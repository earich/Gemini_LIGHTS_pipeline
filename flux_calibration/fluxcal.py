# coding=utf-8
"""
Keywords worth changing:
"stellar_flux", "stellar_flux_err", "SecondOrder" (1 for J, 0 for H)
"""
import sys
sys.path.append('../')
import numpy as np
import glob, os
from datetime import datetime, timedelta
from gpi_analysis.inputs import getfitskeywords, createrecipe, getfilenames
from astropy.io import fits
from pipeline import wait

import time
from astropy.io import fits

import sanity

# -------------------------
# FUNCTIONS
# -------------------------
def magtoflux(MAGS,ZERO,ZERO_ERR):
    """
    Convert magnitudes to flux.

    --INPUTS:
    MAGS            --  List of [magnitude, uncertainty in magnitude].
    ZERO            --  Zero flux
    ZERO_ERR        --  Error in zero flux.

    --RETURNS:
    [FLUX,MEANERR]  --  List of corresponding flux and uncertainty.
    """
    # Find upper and lower limit of magnitude, uperr and downerr:
    MAG     = MAGS[0]
    UPERR   = MAGS[0]+MAGS[1]
    DOWNERR = MAGS[0]-MAGS[1]
    # Empty list to contain [mag,uperr,downerr] converted to flux.
    FLUXES=[]
    # List for iterating. Multiply different factor of zeroerr each time.
    MZEROS=[0,-1,1]
    for X,M in enumerate([MAG,UPERR,DOWNERR]):
        FLUX = 10**(-M*0.4)*(ZERO + MZEROS[X]*ZERO_ERR)
        FLUXES.append(FLUX)
    # Calculate mean error, the mean of (max. flux - flux) and
    # (flux - min. flux).
    MEANERR = 0.5*( (FLUXES[0]-FLUXES[1]) + (FLUXES[2]-FLUXES[0]))
    # Final flux and uncertainty to keep:
    FLUX    = round(FLUXES[0],3)
    MEANERR = round(MEANERR,  3)
    return [FLUX, MEANERR]

# def fluxtomag(fluxes,zero,zeroerr):
#     flux    = fluxes[0]
#     uperr   = fluxes[0]+fluxes[1]
#     downerr = fluxes[0]-fluxes[1]
#     mags=[]
#     fzeros=[0,-1,1]
#     for x,f in enumerate([flux,uperr,downerr]):
#         mag = 10**(-f*0.4)*(zero+fzeros[x]*zeroerr) INCOMPLETE
#         mags.append(mag)
#     meanerr = 0.5*( (mags[0]-mags[1]) + (mags[2]-mags[0]))
#     mag = round(mags[0],3)
#     meanerr = round(meanerr,3)
#     return [mag, meanerr]

def fluxcal(REDUCED_DIR, QUEUE_DIR, MAG, MAG_ERR, BAND, MAGTOFLUX=1):
    """
    Run flux cal recipes through the pipeline and tabulate output.

    The GPI data reduction pipeline must be running.
    Recipes are created and saved into the Queue directory to be
    executed immediately.

    --INPUTS:
    PODC_DIR    --  Input directory containing podc files to stack.
    MAG         --  Reference magnitude for calibration.
    MAG_ERR     --  Error in reference magnitude.
    QUEUE_DIR   --  GPI_DRP_QUEUE_DIR
    REDUCED_DIR --  GPI_REDUCED_DATA_DIR
    --OPTIONAL:
    MAGTOFLUX   --  Whether to convert MAG,MAG_ERR to FLUX,FLUX_ERR
                    or just set FLUX=MAG, FLUX_ERR=MAG_ERR.
    OUTPUT_DIR -- Specify the output directory.
                  (new to create_recipe in inputs.py 08/JUL/19)

    --OUTPUTS:
    VALS        -- A list of useful values to keep.
    (inputs):   TARGET,     n_podc,     TIME_INT,   FLUX,   FLUX_ERR,
    (outputs):  CALIBFAC,   CALIBERR,   SCALE,      SCALE_PERC_ERR

    Saved pipeline files in REDUCED_DIR:
    _mean.fits      -- All podc aligned and stacked into one image.
    _mean_phot.fits -- Flux-calibrated file.
    """
    # Get all relevant podc file names
    PODC_NAMES = getfilenames(REDUCED_DIR, END = '_podc.fits',
                                 KEEPDIR = 0)
    #BAND       = getfitskeywords(PODC_DIR+PODC_NAMES[0], 'OBSMODE')[0]
    if MAGTOFLUX>0:
        # Set these values depending on waveband.
        # ZERO      --  Zero magnitude.
        # ZERO_ERR  --  Error in zero magnitude.
        # PERC_ERR  --  Extra percentage error to add in quadrature.
        # ORDER     --  Use first-order or second-order spots.
        #               H uses 1st-order (ORDER=0),
        #               J uses 2nd-order (ORDER=1).
        if BAND=='J':
            ZERO        = 1594
            ZERO_ERR    = 27.8
            PERC_ERR    = 20.0
            ORDER       = 1
        else: # H-band:
            ZERO        = 1024
            ZERO_ERR    = 20.0
            PERC_ERR    = 13.0
            ORDER       = 0
        # Convert given magnitudes to fluxes:
        [FLUX,FLUX_ERR] = magtoflux([MAG,MAG_ERR],ZERO,ZERO_ERR)
    else:
        # Use input "magnitudes" as fluxes.
        FLUX,FLUX_ERR   = MAG,MAG_ERR

    # Keywords for the recipe template:
    KEYS_FLUX  = ['stellar_flux="'    + str(FLUX)     +'"',
                  'stellar_flux_err="'+ str(FLUX_ERR) +'"',
                  'SecondOrder="'     + str(ORDER)    +'"'
                  ]

    # Create recipe to align and stack all podc into one data cube.
    # Output file ends '_mean.fits'
    OUT_DIR    = createrecipe(PODC_NAMES,   REDUCED_DIR,
                              REDUCED_DIR,  QUEUE_DIR,
                              RECIPE=3,     RETURN_OUTPUT_DIRECTORY=1)#,
                              # OUTPUT_DIR=OUTPUT_DIR)
    # Wait for the recipe to run in the GPI pipeline.
    time.sleep(30)
    #wait(10,1,CHECK=REDUCED_DIR
    #
    # Find the name of the output file from the previous recipe.
    MEAN_FILE  = getfilenames(OUT_DIR,      END='_mean.fits',
                              KEEPDIR=0)
    # Run a new recipe to calibrate flux from said file.
    # New output file ends '_mean_phot.fits'.
    createrecipe(MEAN_FILE,         OUT_DIR,
                 REDUCED_DIR,       QUEUE_DIR,
                 KEYS=KEYS_FLUX,    RECIPE=4)#,
                 # OUTPUT_DIR=OUT_DIR)
    # Wait for the recipe to run in the GPI pipeline.
    time.sleep(30)
    #
    # Store values in a list.
    # Now get the relevant keywords from that 'phot' output file:
    PHOT_NAMES     = getfilenames(   OUT_DIR,       END='_mean_phot.fits')
    TARGET         = getfitskeywords(PHOT_NAMES[0], 'OBJECT')+'-'+BAND
    TIME_INT       = getfitskeywords(PHOT_NAMES[0], 'ITIME',    HEADER='SCI')
    CALIBFAC       = getfitskeywords(PHOT_NAMES[0], 'CALIBFAC', HEADER='SCI')
    CALIBERR       = getfitskeywords(PHOT_NAMES[0], 'CALIBERR', HEADER='SCI')
    # SCALE is in units of mJy/arcsec^2/(ADU/coadds/sec)
    SCALE          = round( CALIBFAC * 1000.0 * TIME_INT / (0.0144**2.0) ,3)
    SCALE_PERC_ERR = round( (PERC_ERR**2.0 + CALIBERR**2.0)**0.5 ,3)
    # Store these values in 'VALS' list:
    VALS=[]
    VALS.append([TARGET,    len(PODC_NAMES),    TIME_INT,
                 FLUX,      FLUX_ERR,           CALIBFAC,
                 CALIBERR,  SCALE,              SCALE_PERC_ERR])
    VALS.append(['TARGET',  'n_podc',           'TIME_INT',
                 'FLUX',    'FLUX_ERR',         'CALIBFAC',
                 'CALIBERR','SCALE',            'SCALE_PERC_ERR'])

    #Update the combined rstokes file so that the above values are saved into its header
    #    and its units are now in mJy/arcsec^2
    files = glob.glob(REDUCED_DIR + "*_combined_rstokesdc.fits")
    for file in files:
        hdul = fits.open(file)
        hdr = hdul[1].header
        hdr['CALIBFAC'] = CALIBFAC
        hdr['CALIBERR'] = CALIBERR
        hdr['SCALE'] = SCALE
        hdr['SCALE_PERC_ERR'] = SCALE_PERC_ERR
        hdr['BUNIT'] = 'mJy arcsec^2'
        hdul[1].data = hdul[1].data*float(SCALE)/TIME_INT
        new_file = file.split('.')[0] + '_phot.fits'
        hdul.writeto(new_file, overwrite = True)
        hdul.close()
    
    return VALS

def assumed_scale(REDUCED_DIR,SCALE):
    files = glob.glob(REDUCED_DIR + "*_combined_rstokesdc.fits")
    for file in files:
        hdul = fits.open(file)
        hdr = hdul[1].header
        TIME_INT = hdr['ITIME']
        hdr['CALIBFAC'] = -999
        hdr['CALIBERR'] = -999
        hdr['SCALE'] = SCALE
        hdr['SCALE_PERC_ERR'] = -999
        hdr['BUNIT'] = 'mJy arcsec^2'
        hdul[1].data = hdul[1].data*float(SCALE)/TIME_INT
        hdr['HISTORY'] = 'Scale from mean of observations of same night'
        new_file = file.split('.')[0] + '_phot.fits'
        hdul.writeto(new_file, overwrite = True)
        hdul.close()
    

def find_magnitude(RA,DEC,target,filter,table_name='J_H_flux_table.txt'):

    F = open(table_name)
    lines = F.readlines()[2:]
    F.close()

    i = 0
    resid = 100.
    while resid > 0.003:
        line = lines[i].split(',')
        resid = np.sqrt((float(line[1]) - RA)**2. + (float(line[2]) - DEC)**2.)
        i += 1
    
    if filter == 'J':
        mag, emag = float(line[3]),float(line[4])
    elif filter == 'H':
        mag, emag = float(line[5]),float(line[6])
    else:
        print('A filter other than J or H was specified in the IFSFILT header keyword \n')
        print('exiting pipeline now!')
        exit()
    print('Magnitude: ' + str(mag) + ' +/- ' + str(emag) + ' for object ' + target + '-' + filter)
    return mag, emag

# %%%%%%%%%%% MAIN FUNCTION %%%%%%%%%%%
def main(QUEUE_DIR, REDUCED_DIR):
    """
    Run flux cal recipes through the pipeline and tabulate output.

    The GPI data reduction pipeline must be running.
    Recipes are created and saved into the Queue directory to be
    executed immediately.

    --INPUTS:
    DIRS        --  List of directories containing podc. Each directory
                    is considered separately for the flux calibration.
    MAGS        --  List of reference magnitudes and uncertainties for
                    each star. Must be same length as DIRS.
                    MAGS=[[flux1, err1], [flux2, err2],...,[fluxn,errn]]
    QUEUE_DIR   -- GPI_DRP_QUEUE_DIR
    REDUCED_DIR -- GPI_REDUCED_DATA_DIR

    --OUTPUTS:
    Saved pipeline files in REDUCED_DIR:
    _mean.fits      -- All podc aligned and stacked into one image.
    _mean_phot.fits -- Flux-calibrated file.

    Saved table contents:
    (inputs):  TARGET,     n_podc,     TIME_INT,   FLUX,   FLUX_ERR,
    (outputs): CALIBFAC,   CALIBERR,   SCALE,      SCALE_PERC_ERR
    """
    
    file = glob.glob(REDUCED_DIR + "*_combined_rstokesdc.fits")
    hdul = fits.open(file[0])
    header = hdul[0].header
    RA, DEC, target, filter = header['RA'],header['DEC'], header['OBJECT'], header['IFSFILT'].split('_')[1]
    mag,emag = find_magnitude(RA,DEC,target,filter)
    hdul.close()
    
    #for i in range(len(DIRS)):
    #    print('=========='+DIRS[i].split('/')[-2]+'==========')
    VALS = fluxcal(REDUCED_DIR,QUEUE_DIR, mag, emag,filter, MAGTOFLUX=1)
    print('\n')
    TABLE = sanity.maketable(1, len(VALS[1]), VALS[1])
    TABLE[1,:] = VALS[0]
    SAVE_NAME  = 'fluxcals.txt'
    FMT    = "%15s %6s %9s %6s %9s %13s %13s %8s %15s" #'%s'#"%60s %12s %8.1e %8.1e %8.1e"
    np.savetxt(fname=REDUCED_DIR + SAVE_NAME,X=TABLE,fmt=FMT)
    
    return

# ----------------------------------------------------------------------
# MAIN CONTENT
# ----------------------------------------------------------------------
if __name__ == "__main__":
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
        queue_dir = ROOT_dir + 'queue/'
        reduced_dir = ROOT_dir + 'Reduced/'
        print('Successfuly got directories from .gpienv file')
        print(queue_dir,reduced_dir)
    except IOError:
        print('Using directories sepcified in pipeline.py')
        # GPI_DRP_QUEUE_DIR and GPI_REDUCED_DATA_DIR from the gpienv.
        queue_dir   = \
        '/Users/earich/work_reduction/queue/'
        #'/Users/al630/gpi_pipeline_original_1.4.0_source/data/queue/'
        reduced_dir = \
        '/Users/earich/work_reduction/Reduced/'
        print(queue_dir,reduced_dir)
    
    #Read in arguemnts to find file that needs to be converted
    reduced_dir += sys.argv[1]
    if len(sys.argv) > 2:
        SCALE = sys.argv[2]
        assumed_scale(reduced_dir,SCALE)
    else: 
        main(queue_dir, reduced_dir)

    # dir = '/Users/al630/Documents/GPI/GPI LLP data paper/'+\
    #       'Reduction details/data paper reductions/'
    # podc_dirs = [
    # 'FU Ori-J',
    # 'Hen 3-365-J',
    # 'hd 139614-J',
    # 'hd45677-J',
    # 'hd45677-H',
    # 'hd 139614-H',
    # 'mwc 789-H'
    # ]
    # for i in range(len(podc_dirs)):
    #     podc_dirs[i] = dir+podc_dirs[i]+'/'
    #
    # mags = [
    # [6.519, 0.023], # FU Ori J
    # [6.217, 0.037], # Hen 3-365 J
    # [7.669, 0.026], # HD 139614 J
    # [7.242, 0.026], # HD 45677 J
    # [6.347, 0.023], # HD 45677 H
    # [7.333, 0.040], # HD 139614 H
    # [7.528, 0.026], # MWC 789 H
    # ]


