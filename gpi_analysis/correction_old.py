# coding=utf-8
"""
GPI analysis functions for pre/post corrections of instrumental polarisation.
"""
import numpy as np
from scipy.optimize import minimize
from . import analysis
from . import inputs
from astropy.io import fits
import glob, os
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import matplotlib
from mpl_toolkits.axes_grid1 import make_axes_locatable

# from analysis import make_radialstokes
# from inputs import getfitsdata, savefitsdata

# -------------------------
# PRECORRECTION
# -------------------------
def precorrection(XY,YY, MAXR=85.0, XTEMP=140.0, YTEMP=140.0):
    """
    Equalise flux in the two polarization states of a podc.

    Input:
    [XY,YY] -- The two polarisation states in the podc.
    MAXR    -- maximum distance from center pixel to compare summed fluxes.
               MAXR is in pixel units, e.g. MAXR=85pixels.
    Returns:
    [XY,YY] -- List of the two corrected ordinary and extraordinary images.
               These can then be saved as an updated podc fits file.
    """
    # Sum up all the flux within some radius of the centre pixel
    XYSUM = radialsum(XY, MAXR, XTEMP, YTEMP)
    YYSUM = radialsum(YY, MAXR, XTEMP, YTEMP)
    # Find the ratio of the summed flux in XY map compared with YY map.
    RATIO = XYSUM / YYSUM
    # If RATIO>1.0 (i.e. XY>YY), this makes YY larger  and XY smaller.
    # If RATIO<1.0 (i.e. XY<YY), this makes YY smaller and XY larger.
    # The post-correction ratio is always 1.0.
    XY_PRECORRECTED = XY * np.sqrt(1.0/RATIO)
    YY_PRECORRECTED = YY * np.sqrt(RATIO)
    XYSUM_PRECORRECTED = radialsum(XY_PRECORRECTED, MAXR)
    YYSUM_PRECORRECTED = radialsum(YY_PRECORRECTED, MAXR)
    # Find the ratio of the summed flux in XY map compared with YY map.
    RATIO_PRECORRECTED = XYSUM_PRECORRECTED / YYSUM_PRECORRECTED
    print(RATIO, RATIO_PRECORRECTED)
    # Return the corrected maps:
    return [XY_PRECORRECTED,YY_PRECORRECTED]

def circular_mask(image, center, radius=15., filler=0.0):
    Y, X = np.ogrid[:len(image[:,0]),:len(image[0,:])]
    dist_from_center = np.sqrt((X - center[0])**2 + (Y-center[1])**2)
    
    dist_from_center[dist_from_center < radius] = filler
    dist_from_center[dist_from_center >= radius] = 1.0
    image = dist_from_center
    return image

# -------------------------
# fQ subtraction correction
# -------------------------
def stellar_polarisation_subtraction(STOKES,stokes_name,MINR=0.0,MAXR=140.0,Bcent_stokes=(0,0)):
    """
    Recommended subtraction area: 70<r<80pix.
    """
    XPIX,YPIX = analysis.getmidvalues(STOKES[1],MINR=MINR,MAXR=MAXR)
    XOUT,YOUT = analysis.getmidvalues(STOKES[1],MINR=80,MAXR=140)
    """
    Test that creates a mask for Q/I values > 0.1. Those pixles will not be included
    in stellar polarisation calculation.
    """
    mask_limit = 0.05
    #mask_limit = 'nan'
    print('mask ratio limit: ', mask_limit)
    if Bcent_stokes[0] != 0:
        maskI = circular_mask(STOKES[0],([int(Bcent_stokes[0]),int(Bcent_stokes[1])]),radius=40,filler=np.nan)
        print('Binary companion masked at: (x,y) = ',Bcent_stokes)
    else:
        maskI = np.ones((STOKES[0].shape))
    ratio = np.abs(STOKES[1]/STOKES[0]*maskI)
    ratioU = np.abs(STOKES[2]/STOKES[0]*maskI)
    total_ratio = np.copy(ratio)
    total_ratio[total_ratio>-1000000]=1.0
    total = np.nansum(total_ratio)
    if mask_limit == 0.05:
        floor = np.nanstd((STOKES[0]*maskI)[XOUT,YOUT])
        print('Minimum counts used when finding ratio: ' + str(floor))
        print('floor used: ' + str(200.))
        floor = 66.67
        ratio[ratio>mask_limit] = np.nan
        ratio[ratioU>mask_limit] = np.nan
        ratio[STOKES[0]<floor*3.] = np.nan

    ratio[ratio>-100000] = 1.0
    print('Number of masked pixels: ', total - np.nansum(ratio))
    
    #Back to your normaly scheduled TV
    MEANI = np.nanmean((STOKES[0]*ratio)[XPIX,YPIX])
    MEANQ = np.nanmean((STOKES[1]*ratio)[XPIX,YPIX])
    MEANU = np.nanmean((STOKES[2]*ratio)[XPIX,YPIX])
    FQ = MEANQ/MEANI
    FU = MEANU/MEANI
    print('Mean FQ: ', FQ, '   Mean FU: ', FU)
    
    #print('Actually use FQ: ', FQ, '   FU: ', FU)
    # make new stokes maps:
    Qstell, Ustell = STOKES[0]*(FQ), STOKES[0]*(FU)   
    QSTAR = STOKES[1] - Qstell
    USTAR = STOKES[2] - Ustell
 
    #Add FQ and FU values to the podc_stokesdc file
    hdul = fits.open(stokes_name)
    hdr = hdul['sci'].header
    if 'FQ' in hdr:
        hdr['FQ'] = (FQ, 'Fraction of Q Stellar Polarization Removed')
        hdr['FU'] = (FU, 'Fraction of Q Stellar Polarization Removed')
        hdr['MeanI'] = (MEANI, 'Mean I intensity used to remove stellar polarization')
        hdr.append(('HISTORY','For Stellar Polarization, used MINR=' + str(MINR) + ' MAXR=' +\
                    str(MAXR) + 'and mask_limit=' + str(mask_limit)))
    else:
        hdr.append(('FQ', FQ, 'Fraction of Q Stellar Polarization Removed'))
        hdr.append(('FU', FU, 'Fraction of U Stellar Polarization Removed'))
        hdr.append(('HISTORY','For Stellar Polarization, used MINR=' + str(MINR) + ' MAXR=' +\
                    str(MAXR) + 'and mask_limit=' + str(mask_limit)))
        hdr.append(('MeanI',MEANI, 'Mean I intensity used to remove stellar polarization'))
    hdul.writeto(stokes_name, overwrite=True)
    hdul.close()
    # new stokes maps:
    
    #SUBT_STOKES = [STOKES[0],STOKES[1],STOKES[2],STOKES[3]]
    SUBT_STOKES = [STOKES[0],QSTAR,USTAR,STOKES[3]]
    return SUBT_STOKES, mask_limit


def polarized_subplots(directory,ax1,ax2,ax3,ax6):
    files = glob.glob(directory + "*_podc_stokesdc.fits")
    files = sorted(files)
    p_avg = 0.0
    pa_avg = 0.0
    I_max = 0.0
    for i in range(0,len(files)):
        hdul = fits.open(files[i])
        hdr = hdul['sci'].header
        fQ = hdr['FQ']
        fU = hdr['FU']
        MeanI = hdr['MeanI']
        if MeanI > I_max:
            I_max = MeanI
        p = np.sqrt(fQ**2. + fU**2.)
        pa = np.arctan2(fU,fQ)*180./np.pi/2.
        p_avg += p*100.
        pa_avg += pa
        ax1.plot(fQ*100.,fU*100.,'o',color='k')
        ax2.plot(i,pa,'o',color='k')
        ax3.plot(i,p*100.,'o',color='k')
        ax6.plot(i,MeanI,'o',color='k')

    ax1.set_ylabel('U%')
    ax1.set_xlabel('Q%')
    ax1.set_xlim((-4.0,4.0))
    ax1.set_ylim((-4.0,4.0))
    p_avg = p_avg/len(files)
    pa_avg = pa_avg/len(files)
    ax1.plot((0,0),(-4.0,4.0),'-',color='k',linewidth=1)
    ax1.plot((-4.0,4.0),(0,0),'-',color='k',linewidth=1)
    ax2.plot((0,3,len(files)),(pa_avg,pa_avg,pa_avg),'--',color='r')
    ax2.text(0,80,'PA Avg = ' + str(round(pa_avg,1)),color='r',Fontsize=12)
    ax3.plot((0,3,len(files)),(p_avg,p_avg,p_avg),'--',color='r')
    ax3.text(0,3.5,'%Pol. Avg = ' + str(round(p_avg,2)),color='r',Fontsize=12)
    ax2.set_xlabel('Cycles') 
    ax2.set_ylabel('Polarization Angle (deg)')
    ax2.set_ylim((-90,90))
    ax3.set_xlabel('Cycles')
    ax3.set_ylabel('% Polarization')
    ax3.set_ylim((0,4))
    ax6.set_ylim((0,I_max + 1000))
    ax6.set_xlabel('Cycles') 
    ax6.set_ylabel('raw Mean I from pol. estimate')
    return

def QUrot_subplots(directory, ax4,ax5,fig):
    file = glob.glob(directory + "*_combined_rstokesdc.fits")
    hdul = fits.open(file[0])
    data = hdul['sci'].data
    hdul.close()

    VMAX = np.nanmax(data[1,:,:])/50
    VMIN = 0.0
    img1 = ax4.imshow(data[1,:,:], origin = 'lower', CMAP='Blues_r', #vmin= -50, vmax = 50,
               norm=matplotlib.colors.SymLogNorm(linthresh=10.0, 
               linscale=1.0, vmin=VMIN, vmax=VMAX))
    divider = make_axes_locatable(ax4)
    cax = divider.append_axes('right', size='5%', pad=0.05)
    fig.colorbar(img1, cax=cax, orientation='vertical')
    ax4.set_xlabel('Combined Qphi image')
    
    VMAX = np.nanmax(data[2,:,:])/100.
    VMIN = -VMAX
    img2 = ax5.imshow(data[2,:,:], origin = 'lower', CMAP='seismic', #vmin= -50, vmax = 50,
               norm=matplotlib.colors.SymLogNorm(linthresh=10.0, 
               linscale=1.0, vmin=VMIN, vmax=VMAX))
    divider = make_axes_locatable(ax5)
    cax = divider.append_axes('right', size='5%', pad=0.05)
    fig.colorbar(img2, cax=cax, orientation='vertical')
    ax5.set_xlabel('Combined Qphi image')

def polarization_plots(directory,name):
    fig = plt.figure(figsize=(19.5,10))
    fig.suptitle(name,fontsize=14)
    gs1 = gridspec.GridSpec(2,3)
    ax1 = plt.subplot(gs1[0])
    ax2 = plt.subplot(gs1[1])
    ax3 = plt.subplot(gs1[2])
    ax4 = plt.subplot(gs1[3])
    ax5 = plt.subplot(gs1[4])
    ax6 = plt.subplot(gs1[5])
    polarized_subplots(directory,ax1,ax2,ax3,ax6)
    QUrot_subplots(directory, ax4,ax5, fig)
    plt.savefig(directory + 'polarization_' + name + '.png',dpi=150)
    print(directory + 'polarization_' + name + '.png')
    plt.close()

"""
def polarized_subplots(directory,gs1):
    files = glob.glob(directory + "*_podc_stokesdc.fits")
    files = sorted(files)
    p_avg = 0.0
    pa_avg = 0.0
    for i in range(0,len(files)):
        hdul = fits.open(files[i])
        hdr = hdul['sci'].header
        fQ = hdr['FQ']
        fU = hdr['FU']
        p = np.sqrt(fQ**2. + fU**2.)
        pa = np.arctan2(fU,fQ)*180./np.pi/2.
        p_avg += p*100.
        pa_avg += pa
        ax1 = plt.subplot(gs1[0])
        ax1.plot(fQ*100.,fU*100.,'o',color='k')
        ax2 = plt.subplot(gs1[1])
        ax2.plot(i,pa,'o',color='k')
        ax3 = plt.subplot(gs1[2])
        ax3.plot(i,p*100.,'o',color='k')

    p_avg = p_avg/len(files)
    pa_avg = pa_avg/len(files)
    ax2.plot((0,3,len(files)),(pa_avg,pa_avg,pa_avg),'--',color='r')
    ax2.text(0,80,'PA Avg = ' + str(round(pa_avg,1)),color='r',Fontsize=12)
    ax3 = plt.subplot(gs1[2])
    ax3.plot((0,3,len(files)),(p_avg,p_avg,p_avg),'--',color='r')
    ax3.text(0,2.5,'%Pol. Avg = ' + str(round(p_avg,2)),color='r',Fontsize=12)
    return

def polarization_plots(directory,description):
    plt.figure(figsize=(15,4))
    gs1 = gridspec.GridSpec(1,3)
    polarized_subplots(directory,gs1)
    ax1 = plt.subplot(gs1[0])
    ax1.set_ylabel('U%')
    ax1.set_xlabel('Q%')
    ax1.set_xlim((-2.0,2.0))
    ax1.set_ylim((-2.0,2.0))
    ax1.plot((0,0),(-2.0,2.0),'-',color='k',linewidth=1)
    ax1.plot((-2.0,2.0),(0,0),'-',color='k',linewidth=1)
    ax2 = plt.subplot(gs1[1])
    ax2.set_xlabel('Cycles') 
    ax2.set_ylabel('Polarization Angle (deg)')
    ax2.set_ylim((-90,90))
    ax3 = plt.subplot(gs1[2])
    ax3.set_xlabel('Cycles')
    ax3.set_ylabel('% Polarization')
    ax3.set_ylim((0,3))
    plt.savefig(directory + 'polarization_' + description + '.png',dpi=150)
    print(directory + 'polarization_' + description + '.png')
    plt.close()
"""

# -------------------------
# FINE CORRECTION - OPTIMISATION FUNCTION
# -------------------------
def weighted_stokesmaps(STOKES_WEIGHTS, STOKES_MAPS):
    """
    Find Q* and U*, the new Q and U weighted by intensity (Avenhaus '18).

    Input:
    STOKES_WEIGHTS -- the values of gamma,c1,c2,c3,c4 to weight Q and U.
    STOKES_MAPS    -- the stokes parameter maps I,Q,U,V.
    """
    # Changeable coefficients:
    if len(STOKES_WEIGHTS)<2:
        GAMMA = STOKES_WEIGHTS[0]
        C1, C2, C3, C4 = 0.0,0.0,0.0,0.0
        print('minimising gamma only', GAMMA)
    else:
        GAMMA, C1, C2, C3, C4 = STOKES_WEIGHTS
    # The input maps:
    I,Q,U,V = STOKES_MAPS
    # Find modified Q and U (Avenhaus 2018's Q* and U*)
    QSTAR = Q + C1*I + C2
    USTAR = U + C3*I + C4
    return QSTAR, USTAR


def find_ur_sum(STOKES_WEIGHTS, EXTRA_ARGS):
    """
    Create weighted Qr and Ur maps, and return sum of counts in Ur.

    Goal: minimise sum(abs(Uphi)) in the region between two circles.

    Input:
    STOKES_WEIGHTS -- the values of gamma,c1,c2,c3,c4 to weight Q and U.
    STOKES_MAPS    -- the stokes parameter maps I,Q,U,V.
    """
    if len(EXTRA_ARGS)==5:
        STOKES_MAPS = EXTRA_ARGS[:4]
        RADIUS = EXTRA_ARGS[4]
        # XCENT=EXTRA_ARGS[-2]
        # YCENT=EXTRA_ARGS[-1]
        YCENT=float(int(0.5*STOKES_MAPS[1].shape[0]))
        XCENT=float(int(0.5*STOKES_MAPS[1].shape[1]))
    else:
        STOKES_MAPS = EXTRA_ARGS
        RADIUS=16.0
        YCENT=float(int(0.5*STOKES_MAPS[1].shape[0]))
        XCENT=float(int(0.5*STOKES_MAPS[1].shape[1]))
    # Find weighted Q* and U* using weights (gamma, c1, c2, c3, c4):
    QSTAR, USTAR = weighted_stokesmaps(STOKES_WEIGHTS, STOKES_MAPS)
    # and use these to find Qr and Ur.
    # n.b. STOKES_WEIGHTS[0] is gamma.
    QR, UR = analysis.make_radialstokes(QSTAR,USTAR,STOKES_WEIGHTS[0],XTEMP=XCENT,YTEMP=YCENT)
    # Find the sum of values in Ur within a radius MAXR of centre pixel.
    USUM = radialsum(UR, MAXR=RADIUS,XTEMP=XCENT,YTEMP=YCENT)
    return USUM


def find_ur_sum_for_centre(CENTS, STOKES_MAPS):
    """
    Create weighted Qr and Ur maps, and return sum of counts in Ur.

    Goal: minimise sum(abs(Uphi)) in the region between two circles.
    edited update this
    Input:
    STOKES_WEIGHTS -- the values of gamma,c1,c2,c3,c4 to weight Q and U.
    STOKES_MAPS    -- the stokes parameter maps I,Q,U,V.
    """
    Q = STOKES_MAPS[1]
    U = STOKES_MAPS[2]
    XCENT=CENTS[0]
    YCENT=CENTS[1]
    QR, UR = analysis.make_radialstokes(Q,U,XTEMP=XCENT,YTEMP=YCENT)
    # Find the sum of values in Ur within a radius MAXR of centre pixel.
    USUM = radialsum(UR, MAXR=85.0)
    return USUM


# def radialsum(MAP, MAXR=10001, XTEMP=10001,YTEMP=10001,ABS=1, MINR=0.0):
#     """
#     Returns sum of fluxes within some distance MAXR from centre of MAP.
#     """
#     X=0
#     if len(MAP.shape)>2:
#         X+=len(MAP.shape)-2#1
#     # Find the distance of each coordinate from the center coordinate.
#     # This makes it easier to define a circle containing values to sum up.
#     # Centre pixel is (XTEMP,YTEMP).
#     if XTEMP>10000 or YTEMP>10000:
#         # i.e. if no sensible input is given:
#         # x/ytemp mark the centre pixel in a stokesdc image.
#         # For standard stokesdc, centre is (140,140).
#         YTEMP=float(int(0.5*MAP.shape[X]))
#         XTEMP=float(int(0.5*MAP.shape[X+1]))
#
#     # Get two grids. YY0 has rows containing -140,-139,...,139,140
#     # and XX0 has columns of the same.
#     XX0, YY0 = np.meshgrid(np.arange(-YTEMP,YTEMP+1.0),
#                            np.arange(-XTEMP,XTEMP+1.0))
#     # Get grid of distances from centre pixel.
#     R0 = (YY0**2.0 + XX0**2.0)**0.5
#     #AS = 1000.0/14.14    # 1 arcsec in pixels (J-band)
#     # Define min and max radii. The equivalent of Avenhaus '18's
#     # postStokesCorr_rInner and postStokesCorr_rOuter from Table 8.
#     # Measured 85 pixels off GPItv for HR 5999. Excludes companion.
#     if MAXR>10000:
#         # i.e. if no sensible input is given:
#         MAXR = (85.0/140.0)*MAP.shape[X]
#
#     # Find coordinates within the allowed range of the centre:
#     XPIX, YPIX = np.where( (R0<MAXR) & (R0>=MINR))
#     # Take an array containing only these values in the allowed range:
#     VALS = MAP[XPIX,YPIX]
#     # Sum up all these values:
#     # (keep 'nansum' to account for 'Not A Number' values outside detector)
#     if ABS>0:
#         SUM = np.nansum(np.abs(VALS))
#     else:
#         SUM = np.nansum(VALS)
#     return SUM


def radialsum(MAP, MAXR=10001, XTEMP=10001,YTEMP=10001,ABS=1, MINR=0.0):
    """
    Returns sum of fluxes within some distance MAXR from centre of MAP.
    """
    # Get values within MINR<=R<=MAXR of MAP centre.
    XPIX,YPIX = getmidvalues(MAP,XTEMP,YTEMP,MAXR,MINR)
    # Take an array containing only these values in the allowed range:
    VALS = MAP[XPIX,YPIX]
    # Sum up all these values:
    # (keep 'nansum' to account for 'Not A Number' values outside detector)
    if ABS>0:
        SUM = np.nansum(np.abs(VALS))
    else:
        SUM = np.nansum(VALS)
    return SUM

def optimise_params(START_PARAMS,EXTRA_ARGS,FUN=find_ur_sum):
    """
    Input:
    start_params -- [gamma,c1,c2,c3,c4]
    stokes_maps  -- [I,Q,U,V]

    Also:
    find_ur_sum  -- finds the sum of the Ur map within some radius of the
                    centre. This is the function to be minimised.

    Output:
    optimised_params  -- the optimised values [gamma,c1,c2,c3,c4]
    optimised_results -- the full results from the optimiser,
                         including the minimised value of sum(Ur).
    """
    # Run the minimization function:
    # minimize is scipy.optimize.minimize
    OPTIMISED_RESULTS = \
    minimize(fun=FUN, x0=START_PARAMS,
                      args=EXTRA_ARGS, method='Nelder-Mead')
    # fun    -- the function to minimise
    # x0     -- array of initial guesses for parameters c1-->c4 and gamma
    # args   -- extra arguments that don't change in minimisation (Q,U,I)
    # method -- use 'Nelder-Mead' to avoid error messages.
    #           Other methods (e.g. the default BFGS) assume that the
    #           function can be differentiated. However we're minimising
    #           a variation of abs(Ur). abs(Ur) has a discontinuity in
    #           gradient at Ur=0 so can't be minimised using those methods.
    # These are the parameters that minimise sum(Ur):
    OPTIMISED_PARAMS = OPTIMISED_RESULTS['x']
    return OPTIMISED_PARAMS, OPTIMISED_RESULTS


# -------------------------
# BIG LOOPING FUNCTION
# -------------------------
# def oldoptimisation(INPUT_DIRECTORY, INPUT_FILENAMES, SHORTFILENAMES, START_PARAMS, OUTPUT_FILE_NAME):
#     """
#     Runs the finecorrection from the Avenhaus '18 paper.
#     Weights Q and U by the intensity I and other coefficients in order to
#     minimise sum(Ur), i.e. minimise the noise in the image.
#
#     Input:
#     filenames        -- Full filenames for each of the input stokesdc.
#                         Includes full path to the file's directory.
#     shortfilenames   -- Shorthand versions of the filenames to use in the
#                         output parameter list and fits files.
#     start_params     -- initial [gamma,c1,c2,c3,c4].
#     output_file_name -- an identifier for the output files.
#
#     Output:
#     stokesdc         -- contains [I,Q*,U*,V] for the newly-found Q* and U*
#                         weighted to minimise sum(Ur).
#     txt file         -- contains the optimal parameters [gamma,c1,c2,c3,c4]
#                         to minimise sum(Ur).
#
#     Returns:
#     --
#     """
#     # Create an empty 2D array to be filled with the new parameters.
#     # I've done it this way to save having to transpose it at the end.
#     # len(shortfilenames) number of rows,
#     # and 6 columns - first is input filename, other five are the parameters
#     # found in the minimisation process.
#     ALLPARAMS= np.zeros((len(SHORTFILENAMES),6),dtype=object)
#     # Fill the first column with the input filenames:
#     ALLPARAMS[:,0] = SHORTFILENAMES
#     # Find the minimised parameters for each file:
#     for F in range(len(INPUT_FILENAMES)):
#         # Import stokes maps using the list of filenames.
#         STOKES_MAPS = inputs.getfitsdata(INPUT_DIRECTORY+INPUT_FILENAMES[F])
#         # EXTRA PART FOR MONNIER 2017:
#         # Mask the values.
#         MASK = np.where( np.abs(STOKES_MAPS[1]) < 10**(-16) )
#         STOKES_MAPS[0][MASK] = np.nan
#         STOKES_MAPS[1][MASK] = np.nan
#         STOKES_MAPS[2][MASK] = np.nan
#         STOKES_MAPS[3][MASK] = np.nan
#         # Run the optimisation, and return the optimised parameters in the list
#         # "params". params = [gamma,c1,c2,c3,c4]
#         PARAMS, OPTIMISED_RESULTS = optimise_params(START_PARAMS,STOKES_MAPS)
#         # Fill the Fth row of the big 2D array with these parameters:
#         # (start at 1st column to avoid overwriting the short file name)
#         ALLPARAMS[F,1:] = PARAMS
#
#         # Save a fits file of the corrected stokes Q and U: Q* and U*.
#         # Make new Q* and U* maps using the optimised parameters in params.
#         QSTAR, USTAR = weighted_stokesmaps(PARAMS, STOKES_MAPS)
#         NEW_STOKES_LIST = [STOKES_MAPS[0],QSTAR,USTAR,STOKES_MAPS[3]]
#         # leeched_file is the name of the fits file to use as the basis for
#         # creating a new fits file. This is a fast way to copy over all the
#         # header information and other extensions.
#         LEECHED_FILE = INPUT_DIRECTORY+INPUT_FILENAMES[F]
#         # new_file_name is name of the newly-created file.
#         NEW_FILE_NAME = INPUT_DIRECTORY + SHORTFILENAMES[F] + '_' + OUTPUT_FILE_NAME +'_stokesdc.fits'
#         inputs.savefitsdata(LEECHED_FILE, NEW_STOKES_LIST, NEW_FILE_NAME)
#         print('Completed file ', F+1 , 'of', len(INPUT_FILENAMES))
#
#     # Save the new-found output parameters to a txt file.
#     # fname is new filename
#     # x is data to be saved
#     # fmt is format (float/integer etc., precision)
#     # header is the line that appears at the top of the file
#     FNAME = INPUT_DIRECTORY + SHORTFILENAMES[0]+'-'+SHORTFILENAMES[-1] +\
#             OUTPUT_FILE_NAME + '_parameters.txt'
#     FMT    = "%25s %14.7e %14.7e %14.7e %14.7e %14.7e"
#     HEADER = 'input file, gamma, c1, c2, c3, c4'
#     np.savetxt(fname=FNAME,X=ALLPARAMS,header=HEADER,fmt=FMT)
#     print('Complete')
#     return


def optimisation(INPUT_DIRECTORY, INPUT_FILENAMES, SHORTFILENAMES, START_PARAMS, OUTPUT_FILE_NAME, XTEMPS=[],YTEMPS=[]):
    """
    Runs the finecorrection from the Avenhaus '18 paper.
    Weights Q and U by the intensity I and other coefficients in order to
    minimise sum(Ur), i.e. minimise the noise in the image.

    Input:
    filenames        -- Full filenames for each of the input stokesdc.
                        Includes full path to the file's directory.
    shortfilenames   -- Shorthand versions of the filenames to use in the
                        output parameter list and fits files.
    start_params     -- initial [gamma,c1,c2,c3,c4].
    output_file_name -- an identifier for the output files.

    Output:
    stokesdc         -- contains [I,Q*,U*,V] for the newly-found Q* and U*
                        weighted to minimise sum(Ur).
    txt file         -- contains the optimal parameters [gamma,c1,c2,c3,c4]
                        to minimise sum(Ur).

    Returns:
    --
    """
    # Create an empty 2D array to be filled with the new parameters.
    # I've done it this way to save having to transpose it at the end.
    # len(shortfilenames) number of rows,
    # and 6 columns - first is input filename, other five are the parameters
    # found in the minimisation process.
    ALLPARAMS= np.zeros((len(SHORTFILENAMES),6),dtype=object)
    # Fill the first column with the input filenames:
    ALLPARAMS[:,0] = SHORTFILENAMES
    # Find the minimised parameters for each file:
    for F in range(len(INPUT_FILENAMES)):
        # Import stokes maps using the list of filenames.
        STOKES_MAPS = inputs.getfitsdata(INPUT_DIRECTORY+INPUT_FILENAMES[F])
        # EXTRA PART FOR MONNIER 2017:
        # Mask the values.
        MASK = np.where( np.abs(STOKES_MAPS[1]) < 10**(-16) )
        STOKES_MAPS[0][MASK] = np.nan
        STOKES_MAPS[1][MASK] = np.nan
        STOKES_MAPS[2][MASK] = np.nan
        STOKES_MAPS[3][MASK] = np.nan
        # Run the optimisation, and return the optimised parameters in the list
        # "params". params = [gamma,c1,c2,c3,c4]
        if len(XTEMPS)<1 and len(YTEMPS)<1:
            EXTRA_ARGS = list(STOKES_MAPS)
        else:
            EXTRA_ARGS = list(STOKES_MAPS)+[XTEMPS[F],YTEMPS[F]]
        PARAMS, OPTIMISED_RESULTS = optimise_params(START_PARAMS,EXTRA_ARGS)
        # Fill the Fth row of the big 2D array with these parameters:
        # (start at 1st column to avoid overwriting the short file name)
        ALLPARAMS[F,1:] = PARAMS

        # Save a fits file of the corrected stokes Q and U: Q* and U*.
        # Make new Q* and U* maps using the optimised parameters in params.
        QSTAR, USTAR = weighted_stokesmaps(PARAMS, STOKES_MAPS)
        NEW_STOKES_LIST = [STOKES_MAPS[0],QSTAR,USTAR,STOKES_MAPS[3]]
        # leeched_file is the name of the fits file to use as the basis for
        # creating a new fits file. This is a fast way to copy over all the
        # header information and other extensions.
        LEECHED_FILE = INPUT_DIRECTORY+INPUT_FILENAMES[F]
        # new_file_name is name of the newly-created file.
        NEW_FILE_NAME = INPUT_DIRECTORY + SHORTFILENAMES[F] + '_' + OUTPUT_FILE_NAME +'_stokesdc.fits'
        inputs.savefitsdata(LEECHED_FILE, NEW_STOKES_LIST, NEW_FILE_NAME)
        print('Completed file ', F+1 , 'of', len(INPUT_FILENAMES))

    # Save the new-found output parameters to a txt file.
    # fname is new filename
    # x is data to be saved
    # fmt is format (float/integer etc., precision)
    # header is the line that appears at the top of the file
    FNAME = INPUT_DIRECTORY + 'newcentres_individual_' + \
            OUTPUT_FILE_NAME + '_parameters.txt'
    # FNAME = INPUT_DIRECTORY + SHORTFILENAMES[0]+'-'+SHORTFILENAMES[-1] +\
    #         OUTPUT_FILE_NAME + '_parameters.txt'
    FMT    = "%25s %14.7e %14.7e %14.7e %14.7e %14.7e"
    HEADER = 'input file, gamma, c1, c2, c3, c4'
    np.savetxt(fname=FNAME,X=ALLPARAMS,header=HEADER,fmt=FMT)
    print('Complete')
    return

if __name__ == '__main__':
    exit()
