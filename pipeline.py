#!/usr/bin/env python
"""
GPI pipeline.
Recommend passing the raw data directory into sanitycheck.py first
(there's a line to do it in this script too).

Command line usage: pipeline.py <rawdata directory> <extra args>

Extra args are used in the recipe creations.
Args starting "0-" are used in the podc creation.
Args starting "1-" are used in the stokes creation.

Recommended extra args:
'0-x0="148"'', '0-y0="140"'      - Change initial guess for PSF centre.
'0-search_window="5"'            - Change search window for PSF centre.
'0-x_off="0.0"', '0-y_off="0.0"' - Change flexure offset.

Comments made by Evan Rich labed with EAR. Comments are mostly to help Evan understand the program
"""
# %%%%%%%%%%% IMPORTS %%%%%%%%%%%
import sys
import os
import numpy as np
import glob, os
from astropy.io import fits
# Anna's functions:
# import centre_dists
import sanity
from gpi_analysis.inputs     import getfilenames, getfitskeywords,\
                                    getfitsdata, savefitsdata,\
                                    getfilenames_l, createrecipe
from gpi_analysis.analysis   import combinestokes, make_radialstokes,\
                                    make_linpolint, getmidvalues, stokesdc_bootstrap
from gpi_analysis.correction import stellar_polarisation_subtraction, polarization_plots, qphi_uphi_subplots
from gpi_analysis.plot       import plotfour, imshow_fancy, get_vlims,\
                                    gather_flexure, plot_fluxchange
from centering_routines import mask_fits, rename_fits, update_cent, plot_center,find_stokesBcent

import datetime
import os.path
from os import path


# %%%%%%%%%%% FUNCTIONS %%%%%%%%%%%
def wait(TIME_FOR_EACH, NUMBER_OF_FILES=1, CHECK=''):
    """
    # If CHECK = '': Wait for the podc creation recipe to run in the GPI pipeline.
    # Should take around 15s to create each podc, 6s per stokesdc.
    """
    import time
    from astropy.io import fits
    if CHECK == '':
        # Time to wait in seconds:
        time_to_wait = NUMBER_OF_FILES*TIME_FOR_EACH
        mins = int(time_to_wait/60.0)
        secs = int(time_to_wait%60.0)
        print('Wait ' + str(mins)+'m'+str(secs)+ 's. ' +\
              'Started at ' + time.strftime('%X %x %Z'))
        time.sleep(time_to_wait)
        print('Wait ended.')
        return
    else: # If check = 'filename', then the loop will look for the file to be created, otheriwise it will wait
        # Until the file has been created or maximum amount of time (time_to_wait) has ellpased.
        print(CHECK)
        stime = time.time()
        time_to_wait = NUMBER_OF_FILES*TIME_FOR_EACH
        dtime = 0.0
        mins = int(time_to_wait/60.0)
        secs = int(time_to_wait%60.0)
        print('Started at ' + time.strftime('%X %x %Z'))
        print('Max waiting time: ' + str(mins)+'m'+str(secs)+ 's. ')
        while dtime < time_to_wait:
            exist = path.exists(CHECK)
            if exist == True:
                print('\n Wait Complete')
                dtime = time_to_wait
            else:
                print('.', end="", flush=True)
                #print('waiting another half minute minute. Started at ' + time.strftime('%X %x %Z') )
                time.sleep(10.0)
                dtime = time.time() - stime

# %%%%%%%%%%% MAIN FUNCTION %%%%%%%%%%%
def main():
    #Import arguments and flags from the origional commandline prompt. EAR
    # Example on how to execute the code from command line:
    #     python3 pipeline.py template.txt /location_of_raw_data_in_raw_directory/
    #     This executes the code, tells the code where it can find the raw data files and
    #     where it can find parameters on how to reduce the data.

    instructions = sys.argv[1]

    rawdata_dir = rawdata_directory + sys.argv[2]

    #This imports keywords that can be appended after the template and raw directory argv.
    #    Current keywords include: first, images, stokes, final
    skip = ''
    if len(sys.argv) > 3:
        skip = sys.argv[3]

    #Sort the keys specified in the template.txt file into lists so they can be inserted into
    #   the appropriate recipe template.
    basic_keys = []
    podc_keys   = []
    stokes_keys = []
    oper_keys = []

    try:
        F = open(rawdata_dir + instructions,'r')
        lines = F.readlines()
        F.close()
        if len(lines) > 0:
            for line in lines:
                if line[0] != '#':
                    keyword = line.strip('\n').split('-')
                    if int(keyword[0])==0:
                        basic_keys.append('-'.join(keyword[1:]))
                    elif int(keyword[0])==1:
                        podc_keys.append('-'.join(keyword[1:]))
                    elif int(keyword[0])==2:
                        oper_keys.append('-'.join(keyword[1:]))
                    else:
                        stokes_keys.append('-'.join(keyword[1:]))
        oper_keys = ':'.join(oper_keys)
        print(basic_keys,podc_keys,stokes_keys,oper_keys)

    except IOError:
        print('Cannot find ' + rawdata_dir + instructions)
        print('exiting pipeline now. Goodbye')
        exit()

    #
    rawdata_files    = getfilenames_l(rawdata_dir,'files.lst',KEEPDIR=0) #Function to return .fits names in text listed file
    target           = getfitskeywords(rawdata_dir+rawdata_files[0], 'OBJECT')
    waveband         = getfitskeywords(rawdata_dir+rawdata_files[0],
                                       'OBSMODE')[0]
    date_obs         = getfitskeywords(rawdata_dir+rawdata_files[-1], 'DATE-OBS')
    date_obs         = ''.join(date_obs.split('-'))

    # Sanity check to check airmass, time, order, and WPANGLE.
    #   Output file of _sanitycheck.txt located in sepcified raw directory
    here_dir = os.path.dirname(os.path.abspath(__file__))
    os.system("python3 "+here_dir+"/sanity.py '"+rawdata_dir+\
              "' 'AIRMASS' 'RMSERR' 'RAWIQ' 'IFSFILT'")
    wait(5)

    #if the first keyword is specified, a recipe is created and run on the first raw data file
    #    named in files.lst
    if skip.find('first') != -1:  # produces the PODC frame of just the first frame
        output_directory = createrecipe(rawdata_files[:1], rawdata_dir,
                                       reduced_directory, queue_directory,
                                       basic_keys, RECIPE=0,
                                       RETURN_OUTPUT_DIRECTORY=1)
        exit() # pipeline ends here if first keyword is used

    #EAR: not sure if this is needed anymore
    if skip.find('images') != -1: #reruns calibration images
        #Interpret the necessary oper_keys
        output_directory = reduced_directory + target + '-' + waveband + '/' + date_obs + '/'
        if oper_keys.find('MINR') != -1:
            shold = oper_keys.split('MINR')[1]
            MINR = float(shold.split('"')[1])
            print(MINR)
        else:
            MINR = 70.0
            print('MINR not specified, using: ' + MINR)
        if oper_keys.find('MAXR') != -1:
            shold = oper_keys.split('MAXR')[1]
            MAXR = float(shold.split('"')[1])
            print(MAXR)
            if MAXR == 140.:
                mask_limit=0.05
            else:
                mask_limit=0.0
        else:
            MAXR = 70.0
            mask_limit = 0.0
        polarization_plots(output_directory, target + '-' + waveband + '_' + str(int(MINR)) + '_' + str(int(MAXR)) + \
                           '_' + str(mask_limit),MINR,MAXR)

        #make subplots of Qphi and Uphi of rstokesdc frames
        qphi_uphi_subplots(output_directory,  target + '-' + waveband + '_' + str(int(MINR)) + '_' + str(int(MAXR)) + \
                           '_' + str(mask_limit))
        exit()

    # continuation of pipeline. Skipped if keywords final or stokes are used.
    elif skip.find('stokes') == -1 and skip.find('final') == -1:
        #Check if the data used a coronograph
        coron = getfitskeywords(rawdata_dir+rawdata_files[0], 'OBSMODE')
        if coron.split('_')[1] == 'coron':
            print('using coronographic reduction')
            # Create podc recipe:
            output_directory = createrecipe(rawdata_files, rawdata_dir,
                                           reduced_directory, queue_directory,
                                           basic_keys, RECIPE=5,
                                           RETURN_OUTPUT_DIRECTORY=1)
            wait(30,len(rawdata_files),CHECK=output_directory+rawdata_files[-1].split('.')[0]+'_podc.fits')

            # Write all of the keys to a log file to keep track of what is done in each reduction
            exist = path.exists(output_directory + 'logfile_pipeline.log')
            if exist == True:
                F = open(output_directory + 'logfile_pipeline.log','r+')
                history = F.readlines()[-1]
            else:
                F = open(output_directory + 'logfile_pipeline.log','w+')
                history = ''
            string = 'execute ' + str(datetime.datetime.now()) + ': '
            F.write(string + ' '.join(sys.argv) + '\n')
            F.close()

            os.system('cp ' + rawdata_dir.replace(' ','\ ') + 'files.lst ' + output_directory.replace(' ','\ ') + 'files.lst')

            output_directory = reduced_directory + target + '-' + waveband + '/' + date_obs + '/'

            bpfix_files       = getfilenames_l(output_directory,'files.lst',
                                             END='_podc.fits',
                                             KEEPDIR=0)
            if oper_keys.find('Bcent') != -1:
                shold = oper_keys.split('Bcent')[1]
                Bcent = (float(shold.split('(')[1].split(',')[0]),float(shold.split(',')[1].split(')')[0]))
                mask_fits(bpfix_files,output_directory,reduced_directory + 'calibrations/',Bcent=Bcent)
            else:
                mask_fits(bpfix_files,output_directory,reduced_directory + 'calibrations/')

            mask_files       = getfilenames_l(output_directory,'files.lst',
                                             END='_podc_mask.fits',
                                             KEEPDIR=0)

            output_directory = createrecipe(mask_files, output_directory,
                                           reduced_directory, queue_directory,
                                           podc_keys, RECIPE=1,
                                           RETURN_OUTPUT_DIRECTORY=1)

            wait(30,len(rawdata_files),CHECK=output_directory+rawdata_files[-1].split('.')[0]+'_podc_mask_.fits')
            PSFCENT_files       = getfilenames_l(output_directory,'files.lst',
                                             END='_podc_mask_.fits',
                                             KEEPDIR=0)
            rename_fits(PSFCENT_files,output_directory,'_PSFCENT.fits')
            update_cent(PSFCENT_files,output_directory)


        else:
            print('using non-coronographic reduction')
            output_directory = createrecipe(rawdata_files, rawdata_dir,
                                           reduced_directory, queue_directory,
                                           basic_keys, RECIPE=5,
                                           RETURN_OUTPUT_DIRECTORY=1)
            wait(30,len(rawdata_files),CHECK=output_directory+rawdata_files[-1].split('.')[0]+'_podc.fits')
            os.system('cp ' + rawdata_dir.replace(' ','\ ') + 'files.lst ' + output_directory.replace(' ','\ ') + 'files.lst')
        # outputs *_podc_centers.png files in target reduced directory
        plot_center(output_directory)

        #Check if Flexure is correct. Outputs *_flexure.png files in target reduced directory
        gather_flexure(rawdata_files,rawdata_dir,output_directory)

        #Check for flux variability and AO performance. First check keys to see if there is
        #   binary. Binary masked if binary. Outputs *_residual.png files in target reduced
        #   directory
        if oper_keys.find('Bcent') != -1:
            shold = oper_keys.split('Bcent')[1]
            Bcent = (float(shold.split('(')[1].split(',')[0]),float(shold.split(',')[1].split(')')[0]))
            print(Bcent)
            plot_fluxchange(rawdata_files,output_directory,Bcent=Bcent)
        else:
            plot_fluxchange(rawdata_files,output_directory)

        print('Plotting flexure and flux change complete')
        # Check PSF centres, update any anomalies:
        #os.system('cp ' + rawdata_dir.replace(' ','\ ') + 'files.lst ' + output_directory.replace(' ','\ ') + 'files.lst')
        # centre_dists.main(output_directory)
        # centers the objects using a radon routine origionally written within the GPI IDL reduction pipeline
        #    now includes masking of the companion
        F = open(output_directory + 'logfile_pipeline.log','a')
        F.write('dir ' +  output_directory + '\n')
        F.close()

    else: #Will read the logfile for the last FILE_OUTPUT_DIRECTORY
        output_directory = reduced_directory + target + '-' + waveband + '/' + date_obs + '/'
        #plot_fluxchange(rawdata_files,output_directory)
        print('Skipping over PODC files. Will use already processed files in : ' + output_directory)

    if skip == '' or skip.find('stokes') != -1:
        # Create stokesdc recipe:
        # To remove PSFcheck step, change to END='podc.fits'
        podc_files       = getfilenames_l(output_directory,'files.lst',
                                        END='_podc.fits',
                                        KEEPDIR=0)

        if skip.find('stokes32') == 0:
            print('Making 1 stokesdc from all podc instead of 1 per HWP cycle.')
            output_directory = createrecipe(podc_files, output_directory,
                                            reduced_directory, queue_directory,
                                            stokes_keys,
                                            RECIPE=2, RETURN_OUTPUT_DIRECTORY=1)
        else:
            for X in range(0,len(podc_files),4):
                output_directory = createrecipe(podc_files[X:X+4], output_directory,
                                                reduced_directory, queue_directory,
                                                stokes_keys,
                                                RECIPE=2, RETURN_OUTPUT_DIRECTORY=1)
        wait(30,len(rawdata_files)*0.25,CHECK= output_directory + podc_files[-1][:-5] + '_stokesdc.fits')

        # Update output_directory here manually if necessary
        # e.g. comment out above and only run stokes part of pipeline.

    else:
        print('Skipping over creation of Stokes files. Will begin correcting already created stokesdc files')

    # Gather stokesdc, stellar polarisation subtraction:
    stokes_files     = getfilenames(output_directory, END='_podc_stokesdc.fits',
                                        KEEPDIR=1)

    multiple_stokes_lists = []
    #Interpret the necessary oper_keys
    if oper_keys.find('MINR') != -1:
        shold = oper_keys.split('MINR')[1]
        MINR = float(shold.split('"')[1])
    else:
        MINR = 70.0
        print('MINR not specified, using: ' + MINR)
    if oper_keys.find('MAXR') != -1:
        shold = oper_keys.split('MAXR')[1]
        MAXR = float(shold.split('"')[1])
    else:
        MAXR = 70.0
        print('MAXR not specified, using: ' + MAXR)
    if oper_keys.find('NBOOT') != -1:
        shold = oper_keys.split('NBOOT')[1]
        NBOOT = int(shold.split('"')[1])
    else:
        NBOOT = 100
    if oper_keys.find('Bcent') != -1:
        shold = oper_keys.split('Bcent')[1]
        Bcent = (float(shold.split('(')[1].split(',')[0]),float(shold.split(',')[1].split(')')[0]))
        Bcent_stokes = find_stokesBcent(rawdata_files,output_directory,stokes_files[0],Bcent)

    for F in range(len(stokes_files)):
        stokes_name  = stokes_files[F]
        stokes       = getfitsdata(stokes_name)
        if oper_keys.find('Bcent') != -1:
            stokes, mask_limit = stellar_polarisation_subtraction(stokes,stokes_name,
                                                            MINR=MINR, MAXR=MAXR, Bcent_stokes=Bcent_stokes)
        else:
            stokes, mask_limit = stellar_polarisation_subtraction(stokes,stokes_name,
                                                            MINR=MINR, MAXR=MAXR)

        i,q,u,v      = stokes
        qphi, uphi   = make_radialstokes(q,u)
        rstokes_list = [i,qphi,uphi,v]
        savefitsdata(output_directory + stokes_name.split('/')[-1], rstokes_list, output_directory +\
                     stokes_name.split('/')[-1].split('_')[0] + '_podc_rstokesdc.fits')
        multiple_stokes_lists.append(stokes)

    #Save the non-bootstrapped files
    combined_stokes_list = combinestokes(multiple_stokes_lists)
    old_stokes       = stokes_files[-1].split('/')[-1]
    new_stokes_stem  = output_directory + target + '-' + waveband +\
                       '_' + old_stokes[:9] + '_' + version + '_m' + str(int(MINR)) + 'M' + str(int(MAXR))
    new_stokes       = new_stokes_stem + '_combined_stokesdc.fits'
    savefitsdata(output_directory+old_stokes, combined_stokes_list, new_stokes)

    i,q,u,v          = combined_stokes_list
    qphi, uphi       = make_radialstokes(q,u)
    rstokes_list     = [i,qphi,uphi,v]
    new_rstokes      = new_stokes_stem + '_combined_rstokesdc.fits'
    savefitsdata(output_directory+old_stokes, rstokes_list, new_rstokes)

    #make QU, PA, and %Pol plots
    polarization_plots(output_directory, target + '-' + waveband + '_' + str(int(MINR)) + '_' + str(int(MAXR)) + \
                       '_' + str(mask_limit),MINR,MAXR)

    #make subplots of Qphi and Uphi of rstokesdc frames
    qphi_uphi_subplots(output_directory,  target + '-' + waveband + '_' + str(int(MINR)) + '_' + str(int(MAXR)) + \
                       '_' + str(mask_limit))

    #Bootstrap the _podc_stokesdc files for future error analysis
    rstokes_files     = getfilenames(output_directory, END='_podc_rstokesdc.fits',
                                    KEEPDIR=1)
    multiple_rstokes_lists = stokesdc_bootstrap(rstokes_files,output_directory,NBOOT,target)

    #Combine the bootstrapped images and save
    combined_rstokes_list = combinestokes(multiple_rstokes_lists)
    old_stokes       = stokes_files[-1].split('/')[-1]
    new_stokes_stem  = output_directory + target + '-' + waveband +\
                       '_' + old_stokes[:14]
    new_rstokes      = new_stokes_stem + '_BScombined_rstokesdc.fits'
    savefitsdata(output_directory+old_stokes, combined_rstokes_list, new_rstokes)

    # Plot output:
    lpi              = make_linpolint(qphi,uphi)
    plotfour([lpi,i,qphi,uphi],  TITLE=target+'-'+waveband,
             SAVENAME=output_directory+target+'-'+waveband+'_maps.png')

    #Check to see if the DRP pipeline inserted any warnings into the headers.
    # Eg. did not preform dark corrections because there was not an appropriate file
    last_file = glob.glob(output_directory + "*_combined_rstokesdc.fits")
    hdul = fits.open(last_file[0])
    hdr = hdul[0].header
    string = str(hdr['History'])
    string = string.split('\n')
    warn = {}
    for i in range(0,len(string)):
        if string[i].find("WARNING") != -1:
            warn.update([(i, string[i] + ' \n' + string[i+1])])
    klaxon = '\n ********************************************** \n '
    if len(warn) == 0:
        print('\n No warnings generated during the reduction. \n')
    else:
        print(klaxon + str(len(warn)) + ' warning(s) generated during the reduction.' + klaxon)
        for item in warn:
            print('Warning in Header History line ' + str(item) + ':')
            print(warn[item])
    #Pipeline complete
    print('Pipeline complete!')
    return


# %%%%%%%%%%% MAIN FUNCTION %%%%%%%%%%%
if __name__ == '__main__':
    try:
        cwd = os.getcwd()
        home = '/' + cwd.split('/')[1] + '/' + cwd.split('/')[2] + '/'
        print(home)
        F = open(home + '.gpienv','r')
        lines = F.readlines()
        F.close()
        line = ';'.join(lines)
        hold = line.split('export GPI_DATA_ROOT=')
        ROOT_dir = hold[1][1:].split(r'"')[0]
        queue_directory = ROOT_dir + 'queue/'
        reduced_directory = ROOT_dir + 'Reduced/'
        rawdata_directory = ROOT_dir + 'Raw/'
        version = 'v' + ROOT_dir.split('/')[-2].split('_')[-1]
        print('Successfuly got directories from .gpienv file')
        print(queue_directory,reduced_directory,rawdata_directory,version)

    except IOError:
        print('Using directories sepcified in pipeline.py')
        # GPI_DRP_QUEUE_DIR and GPI_REDUCED_DATA_DIR from the gpienv.
        queue_directory   = \
        '/Users/earich/work_reduction/queue/'
        #'/Users/al630/gpi_pipeline_original_1.4.0_source/data/queue/'

        reduced_directory = \
        '/Users/earich/work_reduction/Reduced/'
        #'/Users/al630/gpi_pipeline_original_1.4.0_source/data/Reduced/'
        rawdata_directory = '/Users/earich/work_reduction/Raw/'
        print(queue_directory,reduced_directory,rawdata_directory)

    main()
