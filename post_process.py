"""
This python program takes in the bootstrap rstokesdc files and outputs
median maps of I, Qphi, and Uphi
SNR maps of I, Qphi, and Uphi
Radial profilesof I, Qphi, and Uphi
"""

import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import os
import sys
import glob
from gpi_analysis.analysis import radialprofile

def make_SNR_maps(img,error):
    SNR = np.zeros((3,len(img[0,:,0]),len(img[0,0,:])))
    SNR[0,:,:] = np.abs(img[0,:,:])/error[0,:,:]
    SNR[1,:,:] = np.abs(img[1,:,:])/error[1,:,:]
    SNR[2,:,:] = np.abs(img[2,:,:])/error[2,:,:]
    return SNR

def update_flux(name_change, name_phot,dir_change,dir_phot):
    hdu_phot = fits.open(dir_phot + name_phot)
    header_phot = hdu_phot['Sci'].header 
   
    hdu_change = fits.open(dir_change + name_change)
    header_change = hdu_change['Sci'].header
    header_change.set('SATSORDR',header_phot['SATSORDR'])
    header_change.set('CALIBFAC',header_phot['CALIBFAC'])
    header_change.set('CALIBERR',header_phot['CALIBERR'])
    ITIME = float(header_phot['ITIME'])
    header_change['BUNIT'] = 'mJy/arcsec^2'
    image = hdu_change['Sci'].data
    image = image*float(header_phot['CALIBFAC'])*1000./0.01414/0.01414
    hdu_change['Sci'].data = image
    hdu_change.writeto(dir_change + name_change.split('.fits')[0] + '_phot.fits', overwrite = True)
    hdu_change.close()
    hdu_phot.close()

def main():
    directory = sys.argv[1]
    #directory = '/Users/earich/work_reduction/Reduced/190813/FUOri-J/'
    files = glob.glob(directory + '*_rstokesdc_BS_*.fits')
    for i in range(0,len(files)):
        hdul = fits.open(files[i])
        images = hdul['Sci'].data
        if i == 0:
             data = np.zeros((3,len(images[0,:,0]),len(images[0,0,:]),len(files)))
        data[:,:,:,i] = images[:3,:,:]
        hdul.close()
    
    error = np.zeros((3,len(data[0,:,0,0]),len(data[0,0,:,0])))
    img = np.zeros((3,len(data[0,:,0,0]),len(data[0,0,:,0])))
    for i in range(len(data[0,:,0,0])):
        for j in range(len(data[0,0,:,0])):
            error[0,i,j] = np.nanstd(data[0,i,j,:])
            error[1,i,j] = np.nanstd(data[1,i,j,:])
            error[2,i,j] = np.nanstd(data[2,i,j,:])
            img[0,i,j] = np.nanmedian(data[0,i,j,:])
            img[1,i,j] = np.nanmedian(data[1,i,j,:])
            img[2,i,j] = np.nanmedian(data[2,i,j,:])
    
    """
    SNR = make_SNR_maps(img,error)
    HDUO = fits.open(files[0])
    target = HDUO[0].header['OBJECT']
    hdul2 = HDUO
    hdul2['Sci'].data = SNR
    hdul2.writeto(directory + target + '_BS_SNRmaps.fits', overwrite=True)
    HDUO.close()
    fig, ax = plt.subplots(1)
    scale = 0.01414
    Rbins, Qbins = radialprofile(img[0,:,:])
    Rbins, EQbins = radialprofile(error[0,:,:], ERROR = True)
    Rbins, Qbins, EQbins = np.array(Rbins)*scale, np.array(Qbins), np.array(EQbins)
    #plt.plot(np.arrayplt.plot(Rbins,Qbins,'-', color='b', label='I')
    #ax.fill_between(Rbins, Qbins-EQbins,Qbins+EQbins,color='b', alpha=0.1)Rbins,Qbins,'-',\ 
    #                label='I')
    Rbins, Qbins = radialprofile(img[1,:,:])
    Rbins, EQbins = radialprofile(error[1,:,:], ERROR = True)
    Rbins, Qbins, EQbins = np.array(Rbins)*scale, np.array(Qbins), np.array(EQbins)
    plt.plot(np.array(Rbins)*scale,Qbins,plt.plot(Rbins,Qbins,'--', color='r', label='Q$_\phi$')
    ax.fill_between(Rbins, Qbins-EQbins,Qbins+EQbins,color='r', alpha=0.1)'--', label='Q$_\phi$')
    Rbins, Qbins = radialprofile(img[2,:,:])
    plt.plot(np.array(Rbins)*scale,Qbins,'-.', label='U$_\phi$')
    plt.yscale("log")
    plt.xlabel(r'Radial Distance (")')
    plt.ylabel("Flux (counts)")
    plt.legend(loc=1)
    plt.savefig(directory + directory.split('/')[-3] + '_BS_radialprof.png', dpi=100)
    plt.close()
    """

if __name__ == '__main__':
    main()
