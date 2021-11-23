# coding=utf-8
"""
GPI analysis functions for transforming existing Stokes maps.
"""
import numpy as np
# For Gaussian smoothing
from astropy.convolution import convolve, Gaussian2DKernel
from . import plot
from math import sqrt, log
from astropy.io import fits

# -------------------------
# USEFUL FUNCTIONS
# -------------------------
def getmidvalues(MAP,XTEMP=10001,YTEMP=10001,MAXR=10001,MINR=0.0):
    """
    Get coordinates (indices) of the values within some distance R
    (MINR<=R<=MAXR) of the image (MAP) centre (XTEMP,YTEMP).

    Centre value is (XTEMP,YTEMP). Can set this to anything but
    usually want the image centre.
    """
    # Sanity checks:
    if XTEMP>10000 or YTEMP>10000:
        # i.e. if no sensible input is given:
        # Set XTEMP,YTEMP to centre position of image.
        # For standard stokesdc, centre is (140,140).
        YTEMP=0.5*MAP.shape[0]
        XTEMP=0.5*MAP.shape[1]
    if MAXR>10000:
        # i.e. if no sensible input is given:
        # Sets MAXR to this distance from the image centre.
        # Equivalent of radius 80pix of a 280x280 stokesdc.
        # Usually goes near the detector edge but not touching it.
        MAXR = (80.0/280.0)*MAP.shape[0]
        print('getmidvalues: MAXR',MAXR,' MINR',MINR)
    if MINR>MAXR:
        print('getmidvalues: MINR is greater than MAXR')

    # Get grid of distances from centre pixel, R0 (units of pixels).
    # First get two 2D arrays of the right size. YY0 has rows of
    # -YTEMP,-(YTEMP-1),...,+(YTEMP-1),+YTEMP
    # and XX0 has columns of the same for XTEMP.
    XX0, YY0 = np.meshgrid(np.arange(-YTEMP,YTEMP+1.0),
                           np.arange(-XTEMP,XTEMP+1.0))
    R0 = (YY0**2.0 + XX0**2.0)**0.5
    # Find coordinates XPIX,YPIX within allowed range of the centre:
    XPIX, YPIX = np.where( (R0<=MAXR) & (R0>=MINR))
    return XPIX, YPIX

# -------------------------
# MAKING NEW MAPS
# -------------------------
def make_linpolint(Q,U):
    """
    Combine Q and U maps to find linear polarized intensity.

    Returns LP, a map of linear polarized intensity.
    """
    LP = ( Q**2.0 + U**2.0 )**0.5
    return LP


def make_radialstokes(Q, U, GAMMA=0.0, XTEMP=10001, YTEMP=10001):
    """
    Given input Q0 and U0 maps, returns radial Qr and Ur maps.

    Gamma is offset for polar angle - default 0.0.
    XTEMP,YTEMP are the coordinates of the centre of the polar angle
    PHI - usually want to use the image centre.
    """
    if XTEMP>10000 or YTEMP>10000:
        # i.e. if no sensible input is given:
        # Find centre coordinates XTEMP,YTEMP for polar angle.
        # For standard stokesdc, centre is (140,140).
        YTEMP=float(int(0.5*Q.shape[0]))
        XTEMP=float(int(0.5*Q.shape[1]))

    # Make yy0, grid with each row    going 0,1,2,3,...,MAX
    # and  xx0, grid with each column going 0,1,2,3,...,MAX
    MAX = Q.shape[0]
    XX0, YY0 = np.meshgrid(np.arange(0.0,MAX),np.arange(0.0,MAX))
    # Make polar angle grid. Use it to find Qphi, Uphi.
    PHI = np.arctan2((YY0-YTEMP+10.0**-10.0),(XX0-XTEMP)) + GAMMA
    QPHI =  U*np.sin(2.0*PHI) + Q*np.cos(2.0*PHI)
    UPHI = -Q*np.sin(2.0*PHI) + U*np.cos(2.0*PHI)
    return QPHI, UPHI


# -------------------------
# EDIT EXISTING MAPS
# -------------------------
def better_smooth_map(MAP,STD=0.0,MAXR=10001):
    """
    Smooth map, accounting for high values and NaN border.

    Smooths with Gaussian kernal with standard deviation STD pixels.
    Masks large values located outside MAXR.
    """
    if STD==0.0:
        STD = 2.1 / (2.0*sqrt(2.0*log(2.0)))
    # Mask out unusually high values near the map edge.
    # Find values within some distance of the centre.
    XPIX, YPIX = getmidvalues(MAP,MAXR=MAXR)
    # Make a "mask" so values near centre have 1.0, too far have 0.0.
    MASK_MAP = np.full(MAP.shape,0.0)
    MASK_MAP[YPIX,XPIX]=1.0

    # Find colour scale limits using the central quarter of the image:
    VLIMS         = plot.get_vlims(MAP,0.05)
    # Find locations of values far outside this colour range:
    HIGHVALS      = np.where( ((np.abs(MAP)>VLIMS[1]*0.1) & (MASK_MAP<1.0)) )
    # HIGHVALS      = np.where( ((np.abs(MAP)>VLIMS[1]) & (MASK_MAP<1.0)) )
    # Set these values to NaN later. Shrink edge first.
    # (else map would form holes where we set the NaNs).

    # Shrink edge of map.
    # We want to exclude edge values from the later interpolation.
    # Create a new map to populate with the valid old map values.
    NEW_MAP = np.full(MAP.shape,np.NaN)
    # Find non-NaN values in the map.
    NOT_NAN = np.where(np.abs(MAP)>=0.0)
    # Define min and max values to loop over.
    # Avoids sampling outside the array - large NaN border already.
    XMIN = min(NOT_NAN[0])
    XMAX = max(NOT_NAN[0])
    YMIN = min(NOT_NAN[1])
    YMAX = max(NOT_NAN[1])
    for X in range(XMIN,XMAX+1):
        for Y in range(YMIN,YMAX+1):
            # If there is a NaN within C pixels, set value to NaN.
            # Largest YMAX is 286 for Hen3-365, so can't have C above 4.
            C = 4
            # if Y-C<YMIN or Y+C>YMAX or X-C<XMIN or X+C>XMAX:
            #     NEW_MAP[Y,X] = np.NaN#MAP[Y,X]
            if np.isnan(np.max(MAP[Y-C:Y+C,X-C:X+C])) == True:
                # (using 'max' rather than 'nanmax' sets result to nan
                # if there is a nan present in these coordinates).
                NEW_MAP[Y,X] = np.NaN
            else:
                NEW_MAP[Y,X] = MAP[Y,X]

    # Set the high values outside the centre to NaN.
    NEW_MAP[HIGHVALS] = np.nan

    # Now do the proper smoothing of this new map.
    SMOOTHED_MAP = smooth_map(NEW_MAP,STD=STD)
    return SMOOTHED_MAP

def smooth_map(MAP,STD=0.0):
    """
    Makes a smoothed-out version of the input stokesdc MAP.
    The smoothing is done using a 2D Gaussian blur with standard
    deviation STD.

    This uses the astropy function "convolve" and kernel
    "Gaussian2DKernel". (Instead using scipy's ndimage.gaussian_filter
    would cause the sea of NaN at the edges to spread inwards and
    shrink the map.)
    """
    if STD==0.0:
        # Set to 2.1 pixel standard deviation.
        STD = 2.1 / (2.0*sqrt(2.0*log(2.0)))
    SMOOTHED_MAP = convolve(MAP, Gaussian2DKernel(STD),
                            boundary='fill',fill_value=0.0,
                            nan_treatment='fill',preserve_nan=True)
    # Extra kwargs set NaN to 0.0 for the sake of smoothing,
    # then revert those values to NaN for the returned smoothed map.
    # Similarly uses 0.0 as placeholder values when smoothing values
    # that are at the edge of the array.
    return SMOOTHED_MAP


def combinestokes(MULTIPLE_STOKES_LISTS):
    """
    Combines multiple maps into one mean map for each Stokes parameter.

    Input:
    MULTIPLE_STOKES_LISTS -- [ [I1,Q1,U1,V1], [I2,Q2,U2,V2], ...] etc.
                             where each list of Stokes maps has come
                             from a different stokesdc.

    Output:
    COMBINED_STOKES_LIST  -- [I,Q,U,V] for the combined maps.
    """
    # Separate the complete STOKES_LISTS into four lists,
    # one of each type of stokes map I,Q,U,V.
    # Turn the group of lists into an array to select by column.
    MULTIPLE_STOKES_LISTS = np.array(MULTIPLE_STOKES_LISTS)
    MULTIPLE_I = MULTIPLE_STOKES_LISTS[:,0]
    MULTIPLE_Q = MULTIPLE_STOKES_LISTS[:,1]
    MULTIPLE_U = MULTIPLE_STOKES_LISTS[:,2]
    MULTIPLE_V = MULTIPLE_STOKES_LISTS[:,3]
    # Combine multiple maps for each Stokes parameter into one mean map.
    # (Can change "nanmean" to "nanmedian" to instead combine maps
    #  using median, but it makes little difference.)
    COMBINED_I  = np.nanmean(MULTIPLE_I,axis=0)
    COMBINED_Q  = np.nanmean(MULTIPLE_Q,axis=0)
    COMBINED_U  = np.nanmean(MULTIPLE_U,axis=0)
    COMBINED_V  = np.nanmean(MULTIPLE_V,axis=0)
    #COMBINED_I  = np.nanmedian(MULTIPLE_I,axis=0)
    #COMBINED_Q  = np.nanmedian(MULTIPLE_Q,axis=0)
    #COMBINED_U  = np.nanmedian(MULTIPLE_U,axis=0)
    #COMBINED_V  = np.nanmedian(MULTIPLE_V,axis=0)
    COMBINED_STOKES_LIST = [COMBINED_I, COMBINED_Q, COMBINED_U, COMBINED_V]
    return COMBINED_STOKES_LIST


# def combinestokes(STOKES_LISTS):
#     """
#     02/JUL/19 - untested but should work fine - replace above?
#     Combines multiple maps into one mean map for each Stokes parameter.
#
#     (Can change "nanmean" to "nanmedian" to instead combine maps
#     using median, but it makes little difference.)
#
#     Input:
#     MULTIPLE_STOKES_LISTS -- [ [I1,Q1,U1,V1], [I2,Q2,U2,V2], ...] etc.
#                              where each list of Stokes maps has come
#                              from a different stokesdc.
#     Output:
#     COMBINED_STOKES_LIST  -- [I,Q,U,V] for the combined maps.
#     """
#     COMBINED_STOKES = [np.nanmean(STOKES_LISTS[:,X],axis=0) for X in range(0,4)]
#     return COMBINED_STOKES

# -------------------------
# RADIAL PROFILES
# -------------------------
def radialprofile(MAP, SMOOTH=1, BWIDTH=1.0, XTEMP=10001, YTEMP=10001, ERROR=False):
    """
    Find variation in flux with distance from centre of map.

    For each binned distance, average the flux from all pixels
    within a small range of radial distances.

    BWIDTH -- bin width in pixels.
    XTEMP, YTEMP -- centre coords. Take r=0 to be here.

    Returns:
    RSMOOTH -- List of sampled distances from centre.
    QSMOOTH -- List of average values for each distance in RSMOOTH.
    """
    # Centre pixel (XTEMP,YTEMP) is (140,140) for default stokesdc.
    # Mark coordinates as (XTEMP,YTEMP):
    if XTEMP>10000 or YTEMP>10000:
        YTEMP=float(int(0.5*MAP.shape[0]))
        XTEMP=float(int(0.5*MAP.shape[1]))
    # Get two grids of distances. YY0 has rows containing
    # -YTEMP,-(YTEMP-1),...,(YTEMP-1),YTEMP
    # and XX0 has columns of the same.
    XX0, YY0 = np.meshgrid(np.arange(-YTEMP,YTEMP+1.0),
                           np.arange(-XTEMP,XTEMP+1.0))
    # Get grid of distances from centre pixel.
    R0 = (YY0**2.0 + XX0**2.0)**0.5
    # Group together distances from centre in some small range of distances.
    RBINS  = np.arange(0.0,np.max(R0)+BWIDTH,BWIDTH)
    QBINS  = []
    for D in range(len(RBINS)-1):
        # Look at each bin containing distances between r and r+dr
        # where dr=BWIDTH.
        # Find pixels where the distance lies in the bin.
        XPIX, YPIX = np.where( (R0<RBINS[D+1]) & (R0>=RBINS[D]))
        # Take values from the input map at pixels lying in this bin.
        QVALS = MAP[XPIX,YPIX]
        # Find the mean of all of the binned values.
        # Mean essentially normalises by number of pixels in each bin.
        # Take absolute value AFTER finding the mean.
        if ERROR == False:
            QBINS.append(abs(np.mean(QVALS)))
        else:
            QBINS.append(np.sqrt(np.sum(QVALS**2.))/float(len(QVALS)))
        # # NEW 19th June 2019 - abs THEN mean.
        # QBINS.append(np.mean(np.abs(QVALS)))
    # Redefine distance list. Essentially remove the first boundary
    # at distance=0.0. This makes RBINS and QBINS the same length.
    RBINS = RBINS[:-1]+BWIDTH
    # Take running average to smooth out the profile somewhat.
    # (each value becomes mean of itself and the values each side of it).
    if SMOOTH>0:
        QBINS, RBINS = radialprofile_smooth(QBINS,RBINS)
    return RBINS, QBINS

def radialprofile_meanratio(QBINS, UBINS, RBINS, RMIN=6.0, RMAX=18.0):
    """
    Find mean Qr/Ur in subset of Qr and Ur radial profiles.

    RMIN and RMAX (units of pixels) are the min/max radii to consider Qr/Ur.
    Default averages values in the inner emission region (6-->18 pixels).
    """
    # Find ratio Qr/Ur, avoiding division by zero.
    QR_UR_RAT = np.array(QBINS) / ( np.array(UBINS)+10**-10 )
    # Find distances between RMIN and RMAX in the radius bins.
    # IER = Inner Emission Region
    R_IER     = np.where( (np.array(RBINS)>=RMIN) & (np.array(RBINS)<=RMAX) )[0]
    # Find mean value of Qr/Ur within this region.
    MEANRATIO = np.nanmean( QR_UR_RAT[R_IER] )
    return MEANRATIO


def radialprofile_smooth(QBINS,RBINS):
    """
    Returns simple running average of input array QBINS, with RBINS to match.
    """
    Q_SMOOTHED = []
    R_SMOOTHED = []
    for E in range(1,len(QBINS)-1):
        # Take the current value of QBINS and the value before and after it.
        QARR = np.array( [QBINS[E-1], QBINS[E], QBINS[E+1]] )
        # New Q is the mean of these three values.
        # This "if" statement accounts for all elements of QARR being NaN.
        # It doesn't suppress the warning entirely but makes it less alarming.
        if np.isnan(QARR).all()==True:
            RUNNING_Q = np.nan
        else:
            RUNNING_Q = np.nanmean(QARR)
        # State the radius for this Q, and append to list.
        # Makes sure that RUNNING_Q and RUNNING_R are the same length.
        RUNNING_R = RBINS[E]
        Q_SMOOTHED.append(RUNNING_Q)
        R_SMOOTHED.append(RUNNING_R)
    return Q_SMOOTHED, R_SMOOTHED

def radialprofile_errorbars(FILES, FLUXSCALE=1.0, SMOOTH=0):
    """
    Calculate errorbar values for radial profiles I,QPHI,UPHI.

    Inputs:
    FILES     -- List of stokesdc to use in calculating errorbars.
                 Use the interim stokesdc that were averaged to create
                 the combined cube.
    FLUXSCALE -- The flux scale factor in ( mJy/arcsec^2 )/(ADU/coadd)/s.
                 If left default, profiles will retain units of ADU/coadd.
    SMOOTH    -- Whether to smooth each stokesdc prior to profile creation.

    Returns:
    STDS      -- Lists of standard deviations to use for errorbars.
                 Contains a list for each map type I,Qphi,Uphi.
                 Each list has length the same as RBINS - sampled identically.
                 Each value in the list is the standard deviation of values
                 in each radial bin across all interim profiles.
    RBINS     -- The binned radii (helps with later plotting).
    """
    # Empty lists to contain profiles. I,Qphi,Uphi.
    PROFS=[[],[],[]]
    # Open each stokesdc separately.
    # Create radial profiles I,Qphi,Uphi.
    # Store all profiles in the PROFS list.
    for FILE in FILES:
        I,Q,U,V     = getfitsdata(FILE)
        QPHI,UPHI   = make_radialstokes(Q,U)
        for X,MAP in enumerate([I,QPHI,UPHI]):
            MAP = shift_centre(MAP,140.0,140.0)
            if FLUXSCALE!=1.0:
                ITIME = getfitskeywords(FILE, 'ITIME', HEADER='SCI')
                print('Scaling to mJy/arcsec^2')
                print('Fluxscale, int. time:',FLUXSCALE,ITIME)
                MAP = MAP*FLUXSCALE/ITIME
            else:
                print('No flux scaling.')
            if SMOOTH>0:
                # Smooth the J-band images:
                print('Smoothing map before profile.')
                MAP    = better_smooth_map(MAP,   MAXR=80.0)
            RBINS, MAP_BINS = radialprofile(MAP,    SMOOTH=1, BWIDTH=1.0)
            PROFS[X].append(MAP_BINS)
    # For each type of profile (I,Qphi, and Uphi), find standard deviations
    # for errorbars. Find the std of values across all cubes at each
    # binned radius, and output this as list STDS.
    STDS = []
    for TYPE in PROFS:
        STDS.append([np.nanstd(np.array(TYPE)[:,X]) for X in range(len(TYPE[0]))])
    return STDS, RBINS


# -------------------------
# CHANGE SCALES, SIZES
# -------------------------
def map_distancenormalise(MAP):
    """
    Scale each pixel according to its distance from the map centre.

    Returns:
    SCALED_MAP -- MAP * (distance from centre, R)^2.
    """
    # Centre pixel is (140,140) for default stokesdc.
    # Mark coordinates as (XTEMP,YTEMP):
    YTEMP=0.5*MAP.shape[X]
    XTEMP=0.5*MAP.shape[X+1]
    # Need to normalise by distance in arcsec, not pixels.
    # Multiply distances by pixel scale PX.
    PX=0.0141 # GPI pixel scale - 1pixel is 0.0141arcsec.

    # Get two grids. YY0 has rows containing -140,-139,...,139,140
    # and XX0 has columns of the same.
    XX0, YY0 = np.meshgrid(np.arange(-YTEMP,YTEMP+1.0)*PX,
                           np.arange(-XTEMP,XTEMP+1.0)*PX)
    # Get grid of distances from centre pixel.
    R0 = (YY0**2.0 + XX0**2.0)**0.5
    # Edit the value at the centre to prevent multiplication by zero.
    R0[int(YTEMP),int(XTEMP)] = 1.0
    SCALED_MAP = MAP * R0
    return SCALED_MAP



def shift_centre(MAP,XCENT=10001,YCENT=10001,XSHIFT=10001,YSHIFT=10001, NEWCOLS=10):
    """
    Make a new version of MAP with a larger size and different centre.

    By default, adds 5 columns of NaN to each side of MAP.
    Can use this function with no shift to pad edges of MAP.
    Increase NEWCOLS if you need to shift by more than 5 pixels.

    Must pass in either the amount to shift (XSHIFT, YSHIFT)
    OR the new required centre coordinates (XCENT, YCENT).
    Recommend XSHIFT,YSHIFT directly for non-square MAP.
    """
    if XCENT>10000 and YCENT>10000:
        XCENT = int(MAP.shape[1]*0.5)
        YCENT = int(MAP.shape[0]*0.5)
    if XSHIFT>10000 and YSHIFT>10000:
        XCENT = int(round(XCENT,0))
        YCENT = int(round(YCENT,0))
        #MAP[YCENT,XCENT] = np.NaN
        # Current centre pixel:
        XTEMP = int(MAP.shape[1]*0.5)
        YTEMP = int(MAP.shape[0]*0.5)
        # Number of rows/columns to shift by:
        XSHIFT = int(XTEMP - XCENT)
        YSHIFT = int(YTEMP - YCENT)

    # Add NEWCOLS total columns to the size of the array.
    # NEWCOLS must be even.
    # NEWCOLS   = 10
    NEWHEIGHT = MAP.shape[0]+NEWCOLS
    NEWWIDTH  = MAP.shape[1]+NEWCOLS
    NEWARRAY  = np.full((NEWHEIGHT,NEWWIDTH), np.nan)

    LEFTCOLS  = int(0.5*NEWCOLS) + XSHIFT
    RIGHTCOLS = int(0.5*NEWCOLS) - XSHIFT
    UPCOLS    = int(0.5*NEWCOLS) - YSHIFT
    DOWNCOLS  = int(0.5*NEWCOLS) + YSHIFT

    NEWARRAY[DOWNCOLS:-UPCOLS,LEFTCOLS:-RIGHTCOLS] = MAP
    return NEWARRAY

def stokesdc_bootstrap(stokes_files, output_directory,bootstrapN,target):
    for i in range(0,len(stokes_files)):
        hdul = fits.open(stokes_files[i])
        images = hdul['Sci'].data
        if i == 0:
            data = np.zeros((3,len(images[0,:,0]),len(images[0,0,:]),len(stokes_files)))
        data[:,:,:,i] = images[:3,:,:]
        hdul.close()

    HDUO = fits.open(stokes_files[0])
    output = np.zeros((3,len(data[0,:,0,0]),len(data[0,0,:,0])))
    hold = np.zeros((3,len(data[0,:,0,0]),len(data[0,0,:,0]),len(stokes_files)))
    multiple_rstokes_lists = []
    for z in range(0,bootstrapN):
        for i in range(len(stokes_files)):
            rand = int(np.random.rand()*len(stokes_files))
            hold[:,:,:,i] = data[:,:,:,rand]
        output_med = np.nanmedian(np.copy(hold), axis=3)
        hdul2 = HDUO
        hdul2['Sci'].data = output_med
        hdul2.writeto(output_directory + target + '_rstokesdc_BS_median_'
                      + str(z) + '.fits',overwrite=True)
        output = np.nanmean(np.copy(hold), axis=3)
        #hdul2 = fits.PrimaryHDU(output[:,:,:])
        hdul2 = HDUO
        hdul2['Sci'].data = output
        hdul2.writeto(output_directory + target + '_rstokesdc_BS_mean_'
                      + str(z) + '.fits',overwrite=True)
        i,qphi,uphi,v = output[0,:,:],output[1,:,:],output[2,:,:],output[0,:,:]*0.0
        multiple_rstokes_lists.append([i,qphi,uphi,v])
    hdul2.close()
    HDUO.close()
    return multiple_rstokes_lists
