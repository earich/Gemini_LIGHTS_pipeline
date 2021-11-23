#!/usr/bin/env python
"""
Routines that help with centering.

The centering technique is to use the RADON transformation that is in the GPI IDL pipeline.
However, we will first mask the image where we take a mask that only leaves regions that are
part of the AO waffle of the primary. Additionally, any binary star is masked out with a circular mask.

These masked .fits files are saved and then passed through the the GPI IDL RADON transform routine.
"""

import numpy as np
from astropy.io import fits
import scipy.ndimage
import matplotlib.pyplot as plt
import os

from gpi_analysis.inputs import getfitskeywords, getfilenames, getfilenames_l, getfitsdata

def circular_mask(image, center, radius=15., filler=0.0):
    Y, X = np.ogrid[:len(image[:,0]),:len(image[0,:])]
    dist_from_center = np.sqrt((X - center[0])**2 + (Y-center[1])**2)

    dist_from_center[dist_from_center < radius] = filler
    dist_from_center[dist_from_center >= radius] = 1.0
    image = dist_from_center*image
    return image

def rename_fits(files,directory,END):
    for name in files:
        os.system('mv ' + directory+name + ' ' + directory+name.split('_')[0] + END)

def update_cent(files,directory):
    for name in files:
        hdu1 = fits.open(directory + name.split('_')[0] + '_podc.fits')
        HDUL3 = hdu1
        #hdu2 = fits.open(directory + name.split('_')[0] + '_PSFCENT.fits')
        XYarray = get_psfcoords(directory + name.split('_')[0] + '_PSFCENT.fits')
        HDUL3['sci'].header['PSFCENTX'] = 148.0 #XYarray[0]
        HDUL3['sci'].header['PSFCENTY'] = 145.0 #XYarray[1]
        HDUL3.writeto(directory + name.split('_')[0] + '_podc.fits' ,overwrite=True)
        hdu1.close()
        HDUL3.close()

def get_psfcoords(FITS_NAME):
    """
    Get coords from podc files:
    """
    XCENT = getfitskeywords(FITS_NAME,'PSFCENTX',HEADER='SCI')
    YCENT = getfitskeywords(FITS_NAME,'PSFCENTY',HEADER='SCI')
    return np.array([XCENT,YCENT])

def plot_center(dir):
    files       = getfilenames_l(dir, 'files.lst', END='_podc.fits')
    coords      = []
    for p in files:
        psfcent = get_psfcoords(p)
        coords.append(psfcent)
    coords      = np.array(coords)
    for i in range(0,len(files)):
        image = getfitsdata(files[i])
        print(files[i],coords[i][0],coords[i][1])
        x,y = int(coords[i][0]), int(coords[i][1])
        ds = 25
        dx, dy = x - ds, y - ds
        fig = plt.figure()
        img1 = plt.imshow(image[0,y-ds:y+ds,x-ds:x+ds],origin='lower')
        for A in coords[int(i//4*4):int(i//4*4)+4]:
            plt.plot(A[0]-dx,A[1]-dy,'x',color='w')
        plt.plot(coords[i][0]-dx,coords[i][1]-dy,'o',color='w')
        plt.savefig(files[i][:-5] + '_centers.png')

def mask_fits(files,directory,reduced_directroy,Bcent=(0,0),Pcent = (148,146)): #Using the AO spot mask is no longer used
    PSFCENTX, PSFCENTY = Pcent[0], Pcent[1]
    if Bcent != (0,0):
        print('Binary star specified at: ',Bcent)
        hdu = fits.open(directory + files[0])
        header = hdu[1].header
        PA = header['AVPARANG']
        PA_first = PA
        theta = np.arctan2(Bcent[1] - PSFCENTY,Bcent[0] - PSFCENTX)
        r = np.sqrt((PSFCENTX-Bcent[0])**2. + (PSFCENTY-Bcent[1])**2.)
        hdu.close()
    else:
        print('No Binary star specified')
    
    """ #Masking only the AO spots is no longer used.
    try:
        hduT = fits.open(reduced_directroy + 'mask_multcirc.fits')
        Tmask = hduT[0].data
        #Tmask = np.ones(Tmask.shape)  #REMOVE THIS LINE
    except IOError:
        print('Could not find the specified mask at: ' + reduced_directroy + mask)
        print('Closing program now')
        exit()
    """
    
    for i in range(0,len(files)):
        hdu1 = fits.open(directory + files[i])
        header = hdu1[1].header
        data = hdu1[1].data
        Tmask = np.ones(data[0].shape)
        if Bcent != (0,0):
            PA = header['AVPARANG']
            ang = theta - (PA_first-PA)*np.pi/180.
            x2,y2 = r*np.cos(ang)+PSFCENTX, r*np.sin(ang)+PSFCENTY
            mask = circular_mask(np.copy(Tmask),(x2,y2),radius=40.)
        else:
            mask = np.copy(Tmask)

        mask = scipy.ndimage.gaussian_filter(mask*1000, 6.0)/1000.
        data[0], data[1] = data[0]*mask, data[1]*mask
        hdu1[1].data = data
        HDUL2 = hdu1
        HDUL2.writeto(directory + files[i].split('_')[0] + '_podc_mask.fits' ,overwrite=True)
        hdu1.close()
        HDUL2.close()
    return

def find_stokesBcent(rawdata_files,output_directory,stokes_file,Bcent):
    hdu = fits.open(output_directory + rawdata_files[0].split('.fits')[0] + '_podc.fits')
    header = hdu[1].header
    PA = header['AVPARANG']
    PSFCENTX = header['PSFCENTX']
    PSFCENTY = header['PSFCENTY']
    PA_first = PA
    theta = np.arctan2(Bcent[1] - PSFCENTY,Bcent[0] - PSFCENTX)
    r = np.sqrt((PSFCENTX-Bcent[0])**2. + (PSFCENTY-Bcent[1])**2.)
    hdu.close()
    
    hdu = fits.open(stokes_file)
    header = hdu[1].header
    ROTANG = header['ROTANG']
    #fudge factor of 2 degrees is included since I know where the companion is in the first frame
    # but not the 4th frame where the ROTANG is calculated from. A mask of sufficient size will be
    # fine with this estimate.
    ang = theta + (ROTANG-2)*np.pi/180.
    x2,y2 = r*np.cos(ang)+140.0, r*np.sin(ang)+140.0
    x2 = 140.+(140.-x2)
    print('Binary star located at (x,y) = ',x2,y2) 
    return (x2,y2) 

