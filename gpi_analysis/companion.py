"""
GPI analysis functions for studying a flux cut across a Stokes map.

Most of the companion-related functions will only be useful in specific cases
where the companion is smeared out horizontally across the stokesdc image.
It's not worth updating this code to apply it to other angles.
"""

import numpy as np
# Use this for fitting a Gaussian model to the binary companion PSF.
from scipy.stats import norm
from scipy.optimize import curve_fit

# -------------------------
# FITTING TO PSF
# -------------------------
def gaussfunction(X, A, MU, STD):
    """
    Gaussian function to fit to companion PSF.

    Input:
    X   -- Range of x-values to create a Gaussian over.
    A   -- Amplitude of Gaussian.
    MU  -- Mean of the Gaussian.
    STD -- Standard deviation of the Gaussian.
    """
    G = A * np.exp( -(X-MU)**2 / (2* STD**2) )
    return G


def gaussmodel(FLUXCUT, X_LIST, MU, STD):
    """
    Fits model Gaussian curve to fluxes across some row of pixels, FLUXCUT.
    Returns the optimal parameters that define this model curve.

    Uses curve_fit from scipy.optimize.

    Inputs:
    FLUXCUT -- list of fluxes across a row of some stokesdc.
    XRANGE  -- the pixel numbers in this row. (e.g. X=range(130,150,1)).
    MU      -- mean of fluxes in FLUXCUT.
    STD     -- standard deviation of fluxes in FLUXCUT.

    Returns:
    HEIGHT -- height of the model curve.
    CENT   -- centre value of x (= centre of model curve).
    SIGMA  -- standard deviation of model curve.
    FWHM   -- full width half maximum of the curve.
    """
    # Fit a curve to the data (x=XRANGE, y=FLUXCUT).
    # p0 list contains initial guess of parameters. These will be varied
    # until the curve is a good fit to the data. Initial guesses:
    # model amplitude = max(FLUXCUT),
    # mean is MU,
    # and standard deviation is STD.
    # The data to be fitted is stored in XRANGE (x pixels)
    # and FLUXCUT (y values).
    POPT, PCOV = curve_fit(gaussfunction, XRANGE, FLUXCUT,
                           p0=[max(FLUXCUT), MU_initial, STD_initial])
    # PCOV is "The estimated covariance of popt." (I'm ignoring it.)
    # POPT is array of the optimal variables for fitting gaussfunction.
    # POPT is a list containing the height of the model curve, the
    # centre value of x, and standard deviation sigma of the curve.
    HEIGHT = POPT[0]
    CENT   = POPT[1]
    SIGMA  = POPT[2]
    # Also find FWHM, which will be useful for comparing the spread of
    # flux across this row in two stokesdc files.
    FWHM   = 2*np.sqrt(2*np.log(2))*abs(POPT[2])
    # (Ok to take abs of sigma. The function is 50/50 about finding a
    # positive value of sigma.)
    # Return model curve parameters:
    return HEIGHT, CENT, SIGMA, FWHM


def fit_gaussian(YCUT, MAP, XRANGE):
    """
    Fit a model Gaussian to the given flux values (MAP[YCUT,XRANGE]).

    Input:
    YCUT   -- y-value of the row of interest.
    MAP    -- map to have its flux cut
    XRANGE -- range of x values over which to examine flux

    Returns:
    A list of the following parameters:
    FLUXCUT -- list of fluxes at the studied pixels
    XRANGE  -- range of x values over which to examine flux
    YGAUSS  -- y-values of model curve
    XGAUSS  -- x-values of model curve
    HEIGHT  -- height of model curve
    CENTRE  -- centre of model curve in x
    SIGMA   -- standard deviation of model curve
    FWHM    -- full width half maximum of model curve
    """
    # Isolate the relevant pixels of the map. Contain the values
    # immediately surrounding the binary companion.
    XMIN = min(XRANGE)
    XMAX = max(XRANGE)
    FLUXCUT = MAP[YCUT,XMIN:XMAX]
    # Find mean and standard deviation of these arrays for model
    # Gaussian fit. These have the same value as instead using
    # np.nanmean(fluxcut) and np.nanstd(fluxcut).
    MEAN, STD = norm.fit(CFLUXCUT)
    # Find the parameters beloning to the optimal Gaussian fit.
    # The model curve's height, x centre, and standard deviation.
    # Also find the FWHM of the curve, for easier comparisons.
    HEIGHT, CENTRE, SIGMA, FWHM = gaussmodel(FLUXCUT, XRANGE, MEAN, STD)
    # Create a model Gaussian curve using these parameters.
    # Set it to be plotted against some x array with lots of values
    # for a smoother curve.
    XGAUSS = np.linspace(XRANGE[0],XRANGE[-1],1000)
    YGAUSS  = gaussfunction(XGAUSS, HEIGHT, CENTRE, SIGMA)
    # Return all the relevant parameters in one list.
    return [FLUXCUT, XRANGE, YGAUSS, XGAUSS, HEIGHT, CENTRE, SIGMA, FWHM]


# -------------------------
# PLOT MODEL CURVES
# -------------------------
def plotfluxfwhm(TITLE, S_FLUXCUT, S_XRANGE, S_YGAUSS, S_XGAUSS,
                 S_HEIGHT, S_CENTRE, S_SIGMA, S_FWHM,
                 C_FLUXCUT, C_XRANGE, C_YGAUSS, C_XGAUSS,
                 C_HEIGHT, C_CENTRE, C_SIGMA, C_FWHM):
    """
    Plots row of fluxes in S/CFLUXCUT, and their model Gaussian fits S/CGAUSS.

    S/CFLUXCUT is plotted against S/CX, and S/CYGAUSS against S/CXGAUSS.
    The model parameters for the curve's HEIGHT, CENTRE, and FWHM
    are used to overlay arrows highlighting the FWHM of the model fit.

    Use this to compare the flux distribution across a binary companion
    for two stokesdc, one created with all 32 podc files, and one made by
    combining eight intermediary stokesdc files.
    """
    # Parameters for overlaying arrows to show FWHM.
    # For single stokesdc curve:
    SARROWLEFT   = S_CENT - 0.5*S_FWHM
    SARROWRIGHT  = S_CENT + 0.5*S_FWHM
    SARROWBOTTOM = 0.5*S_HEIGHT
    SARROWTOP    = 0.5*S_HEIGHT
    # For combined stokesdc curve:
    CARROWLEFT   = C_CENT - 0.5*C_FWHM
    CARROWRIGHT  = C_CENT + 0.5*C_FWHM
    CARROWBOTTOM = 0.5*C_HEIGHT
    CARROWTOP    = 0.5*C_HEIGHT

    # Plot the data points. Scatter points with solid line connecting.
    plt.scatter(S_XRANGE, S_FLUXCUT, marker='o', alpha=0.5,
             color='r',    label='Single stokesdc')
    plt.scatter(C_XRANGE, C_FLUXCUT, marker='s', alpha=0.5,
             color='b',    label='Combined stokesdc')
    # Plot the model Gaussian fit on top.
    plt.plot(S_XGAUSS, S_YGAUSS, color='r')#, label='Gaussian fit to\nsingle stokesdc')
    plt.plot(C_XGAUSS, C_YGAUSS, color='b')#, label='Gaussian fit to\ncombined stokedsc')
    # Plot setup:
    plt.xlabel('Pixel number in x')
    plt.ylabel('Flux (counts)')
    plt.title(TITLE)
    plt.legend(loc='upper left')
    # Overlay arrows showing FWHM.
    # "xy" and "xytext" set the endpoints of the arrows.
    # See here: https://stackoverflow.com/questions/38677467/how-to-annotate-a-range-of-the-x-axis-in-matplotlib
    # The second annotation adds the label off to the side, with FWHM value.
    # FWHM arrow and label for single stokesdc flux:
    plt.annotate('',   xy=(SARROWLEFT,  SARROWBOTTOM),
                 xytext  =(SARROWRIGHT, SARROWTOP),
                 xycoords='data', textcoords='data',
                 arrowprops={'arrowstyle': '<->'})
    plt.annotate('FWHM '+str(round(S_FWHM,2)),
                 xy=(max(S_XRANGE)-5, S_HEIGHT/2.0),
                 ha='center', va='center', color='r')
    # FWHM arrow and label for combined stokesdc flux:
    plt.annotate('', xy=(CARROWLEFT, CARROWBOTTOM),
                 xytext=(CARROWRIGHT, CARROWTOP),
                 xycoords='data', textcoords='data',
                 arrowprops={'arrowstyle': '<->'})
    plt.annotate('FWHM '+str(round(C_FWHM,2)),
                 xy=(max(S_XRANGE)-5, C_HEIGHT/2.0),
                 ha='center', va='center', color='b')
    plt.show()
    return
