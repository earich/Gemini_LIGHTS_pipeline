# coding=utf-8
"""
GPI analysis functions for plotting various things.
"""
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
# For subplots:
import matplotlib.gridspec as gridspec
# For colourbar. Somehow the usual one isn't working today.
from matplotlib.colorbar import Colorbar
from math import pi, sin, cos, atan, sqrt, log
# # Subplots:
# from matplotlib.colorbar import ColorbarBase
# # For ticks:
from matplotlib.ticker import MultipleLocator, FormatStrFormatter, LogFormatter
from gpi_analysis.inputs import getfitskeywords, getfitsdata
from astropy.io import fits
from mpl_toolkits.axes_grid1 import make_axes_locatable

# -------------------------
# DISPLAY MAPS
# -------------------------
def imshow_basic(MAP, VMAX=0.0, VMIN=0.0, LOG=10.0,
                 TITLE='', CMAP='Blues_r', SAVE_NAME=''):
    """
    Basic imshow for displaying maps with colourbar.

    Inputs:
    MAP       -- Map to be shown
    VMIN      -- Minimum value on colourbar
    VMAX      -- Maximum value on colourbar
    LOG       -- 1 for logarithmic colour scale, 0 for linear, -1 for symlog.
    TITLE     -- Title
    SAVE_NAME -- New file name for saving figure
    CMAP      -- Colourmap to use. e.g. 'Blues_r','seismic','gist_heat'
    """
    plt.figure(figsize = (10,10))
    if abs(VMAX)<0.1 and abs(VMIN)<0.1:
        VMAX = np.nanmax(MAP)
        VMIN = 0.001*VMAX
    # Keep origin='lower' to make (0,0) be in the bottom left corner.
    if abs(LOG) < 1.0:
        plt.imshow(MAP, origin='lower', cmap=CMAP, vmin=VMIN,vmax=VMAX),
    elif LOG < 0.0:
        # SymLogNorm is a log colorbar that shows values below zero.
        # (symmetrical log)
        # In SymLogNorm, linthresh is the region that is linearly mapped
        # around zero. Leave linscale as default 1.0 to have the range
        # covered by linthresh shown as one decade.
        plt.imshow(MAP, origin='lower', cmap=CMAP,
                   norm=matplotlib.colors.SymLogNorm(linthresh=10.0,
                   linscale=1.0, vmin=-VMAX, vmax=VMAX) )
    else:
        # Standard log colorbar for all fluxes greater than zero.
        plt.imshow(MAP, origin='lower', cmap=CMAP, norm=LogNorm(VMIN,VMAX))
    plt.colorbar()
    plt.title(TITLE)
    if len(SAVE_NAME)<2:
        plt.show()
    else:
        plt.savefig(SAVE_NAME)
        plt.close()
    return


def imshow_correction_residuals(QR_NOCORRECT,UR_NOCORRECT,
                                QR_CORRECTED,UR_CORRECTED,
                                QR_RESIDUAL, UR_RESIDUAL,
                                VMAX_CORRECTED, VMIN_CORRECTED,
                                VMAX_RESIDUALS, VMIN_RESIDUALS,
                                NEW_FILE_NAME='', CMAP='seismic',
                                TITLE_LEFT='',TITLE_CENT='',TITLE_RIGHT='',
                                PA=10002):
    """
    THIS NEEDS TIDYING

    Plot maps of 1) corrected Qr&Ur with 2) standard Qr&Ur and 3) residuals.

    [LT][CT][] [RT][]
    [LB][CB][] [RB][]
    Left Top   & Left Bottom   - Corrected Qr & Ur
    Centre Top & Centre Bottom - Standard Qr & Ur (no correction)
    Right Top  & Right Bottom  - Residuals, LT-CT & LB-LT.

    Inputs:
    QR_NOCORRECT   -- No correction, Qr map
    UR_NOCORRECT   -- No correction, Ur map
    QR_CORRECTED   -- Corrected Qr map
    UR_CORRECTED   -- Corrected Ur map
    QR_RESIDUAL    -- Residual Qr (corrected - no correction)
    UR_RESIDUAL    -- Residual Ur (corrected - no correction)
    VMAX_CORRECTED -- Colour scale max limit for left four panels.
    VMIN_CORRECTED -- Same but min limit.
    VMAX_RESIDUALS -- Colour scale max limit for right two panels.
    VMIN_RESIDUALS -- Same but min limit.
    NEW_FILE_NAME  -- File name to save the figure as.
    CMAP -- colour map. Use seismic for vmax=-vmin.
    """
    # Set figure size. (8,4) means 800pix wide by 400 tall.
    plt.figure(figsize = (10,6))
    # (3,2) creates a grid of subplots 3 panels tall by 2 wide.
    # Use width ratios to make the colourbars much skinnier than the maps.
    # Extra unused panel (width 8) for a gap down the centre.
    gs = gridspec.GridSpec(2, 6, width_ratios=(16,16,1,8,16,1))
    gs.update(#left=0.0, right=1.0, top=1.0, bottom=0.0,
              hspace=0.0,wspace=0.0)
    # Identify plots by L LEFT, C CENTRE, R RIGHT, T TOP, B BOTTOM.
    # Also Left and Right ColorBAR.
    LT    = plt.subplot(gs[0,0])
    LB    = plt.subplot(gs[1,0],sharex=LT,sharey=LT)
    CT    = plt.subplot(gs[0,1],sharex=LT,sharey=LT)
    CB    = plt.subplot(gs[1,1],sharex=LT,sharey=LT)
    LCBAR = plt.subplot(gs[:,2])
    RT    = plt.subplot(gs[0,4],sharex=LT,sharey=LT)
    RB    = plt.subplot(gs[1,4],sharex=LT,sharey=LT)
    RCBAR = plt.subplot(gs[:,5])
    # Fill in maps to the left four panels and colourbar:
    LT.imshow(QR_NOCORRECT, origin='lower', cmap=CMAP,
              vmin=VMIN_CORRECTED,vmax=VMAX_CORRECTED)
    LB.imshow(UR_NOCORRECT, origin='lower', cmap=CMAP,
              vmin=VMIN_CORRECTED,vmax=VMAX_CORRECTED)
    CAX = CT.imshow(QR_CORRECTED, origin='lower', cmap=CMAP,
                    vmin=VMIN_CORRECTED,vmax=VMAX_CORRECTED)
    CB.imshow(UR_CORRECTED, origin='lower', cmap=CMAP,
              vmin=VMIN_CORRECTED, vmax=VMAX_CORRECTED)
    Colorbar(ax=LCBAR, mappable=CAX, label = 'counts')
    # Fill in maps for the residual maps and colourbar:
    CAX2 = RT.imshow(QR_RESIDUAL, origin='lower', cmap=CMAP,
                     vmin=VMIN_RESIDUALS,vmax=VMAX_RESIDUALS)
    RB.imshow(UR_RESIDUAL, origin='lower', cmap=CMAP,
              vmin=VMIN_RESIDUALS,vmax=VMAX_RESIDUALS)
    Colorbar(ax=RCBAR, mappable=CAX2, label = 'counts')
    LT.set_title(TITLE_LEFT)
    CT.set_title(TITLE_CENT)
    RT.set_title(TITLE_RIGHT)
    # Remove ticks from unseen axes in the center of the panels.
    #LT.set_xticks([])
    #CT.set_xticks([])
    #CT.set_yticks([])
    #CB.set_yticks([])
    #RT.set_xticks([])
    #RT.set_yticks([])
    #RB.set_yticks([])
    plt.setp(LT.get_xticklabels(), visible=False)
    plt.setp(CT.get_xticklabels(), visible=False)
    plt.setp(RT.get_xticklabels(), visible=False)
    plt.setp(CT.get_yticklabels(), visible=False)
    plt.setp(CB.get_yticklabels(), visible=False)
    plt.setp(RB.get_yticklabels(), visible=False)
    # Temporarily use axis titles to label the Qr and Ur maps.
    LT.set_ylabel('Qr')
    LB.set_ylabel('Ur')
    # Create coronagraphic spot.
    # Inner Working Angle IWA is the radius of the occulting spot (pix).
    # Use 6 pixels for J-band, 8 pixels for H-band.
    IWA = 6.0
    SPOTX = []
    SPOTY = []
    for THETA in np.linspace(0,2*pi,100):
        X = 140.0 + IWA*sin(THETA)
        Y = 140.0 + IWA*cos(THETA)
        SPOTX.append(X)
        SPOTY.append(Y)
    LT.plot(SPOTX,SPOTY,color='white')
    LB.plot(SPOTX,SPOTY,color='white')
    CT.plot(SPOTX,SPOTY,color='white')
    CB.plot(SPOTX,SPOTY,color='white')
    RT.plot(SPOTX,SPOTY,color='white')
    RB.plot(SPOTX,SPOTY,color='white')

    # Overplot position angles
    if PA<10001:
        LIMS=[]
        if len(LIMS)<1:
            # LIMS = sharedaxeslimits([QR_NOCORRECT], AXES='pixels')
            # print(LIMS)
            # LIMS = np.array(LIMS)+140
            # print(LIMS)
            LIMS=[0,280,0,280]
        LT = plot_pa_hd45677(PA,LT,LIMS,LENGTH=100,ANNOTATIONCOLOUR='k')
        LB = plot_pa_hd45677(PA,LB,LIMS,LENGTH=100,ANNOTATIONCOLOUR='k')
        CT = plot_pa_hd45677(PA,CT,LIMS,LENGTH=100,ANNOTATIONCOLOUR='k')
        CB = plot_pa_hd45677(PA,CB,LIMS,LENGTH=100,ANNOTATIONCOLOUR='k')
        RT = plot_pa_hd45677(PA,RT,LIMS,LENGTH=100,ANNOTATIONCOLOUR='k')
        RB = plot_pa_hd45677(PA,RB,LIMS,LENGTH=100,ANNOTATIONCOLOUR='k')

    if len(NEW_FILE_NAME)<2:
        plt.show()
    else:
        plt.savefig(NEW_FILE_NAME)
        plt.show()
    return



def sharedaxeslimits(MAP_LIST, LIMS=[], NCOLS=280, AXES='arcsec', SCALE=0.0):
    """
    'lims' will contain the axis limits shared between the subplots.
    Start with the values in input LIMS, then update with the most extreme
    values across all maps in MAP_LIST.
    """
    # MID is the middle pixel of the MAP.
    MID = int(0.5*MAP_LIST[0].shape[0])
    if len(LIMS)<1:
        LIMS = np.array([MID,MID,MID,MID])
    for M in MAP_LIST:
        # Crop off NaN border.
        NOTNAN = np.where(np.abs(M) > 0.0)
        # Axis limits of this subplot only, templims=[xmin,xmax,ymin,ymax].
        MLIMS  = np.array([min(NOTNAN[1]),max(NOTNAN[1]),
                           min(NOTNAN[0]),max(NOTNAN[0])])
        # Keep the lowest values for minima and the highest values for maxima.
        # Minus 140 to use the same condition for all, not having to do half
        # with difference>0 and half with difference<0.
        LIMS   = np.array([max([t,l],key=abs)+MID
                           for t,l in zip(MLIMS-MID,LIMS-MID)])
    if AXES == 'arcsec':
        PX = 0.0141
    elif AXES == 'physical':
        # 1arcsec is 1.0/0.0141 pixels.
        PX = SCALE*0.0141
    else:
        PX = 1.0
    LIMS = (LIMS-MID)*PX
    if NCOLS<280:
        LIMS = LIMS*(NCOLS+1)
    return LIMS



def plotfour(MAPS,VLIMS=[],IWA=0.0,TITLE='',LIMS=[],SAVENAME='',CBARLABEL=''):
    """
    Display four maps together - LPI, I, Qphi, Uphi.

    Inputs:
    MAPS      -- List of maps to plot: [LPI,I,Qphi,Uphi].
                 Each map is a 2D array.
    Optional:
    VLIMS     -- List containing colourbar limits (v limits) for each
                 subplot. VLIMS must be same length as MAPS.
    IWA       -- Inner Working Angle, radius of spot in pixels.
                 (6.0 for J-band, 8.0 for H-band).
    TITLE     -- Title displayed above plots.
    LIMS      -- Axis limits for all plots. [XMIN,XMAX,YMIN,YMAX].
    SAVENAME  -- Save name for the resulting plot.
    CBARLABEL -- Label for the colourbars.
    """
    # Set up labels, colourmaps etc.
    LABELS  = [r'LPI=$\sqrt{Q^2+U^2}$'+'\n(log)',
               'Stokes I \n(log)',
               r'$Q_{\phi}$'+'\n(symlog)',
               r'$U_{\phi}$'+'\n(symlog)']
    CMAPS   = ['inferno', 'inferno', 'PuOr',   'PuOr'  ]
    VSCALES = ['log',     'log',     'symlog', 'symlog']
    COLOURS = ['r',       'r',       'r',      'r'     ]

    # If no colourbar limits are given, calculate them now.
    if len(VLIMS)<1:
        VLIMS=[get_vlims(MAPS[X],SCALE=0.005,TYPE=VSCALES[X])\
               for X in range(len(MAPS))]


    # --- FIGURE SETUP ---
    # Set figure size. (8,7) means 800pix wide by 700 tall.
    plt.figure(figsize = (8,7))
    # plt.title(TITLE)
    plt.suptitle(TITLE)
    # (3,2) creates a grid of subplots 3 panels tall by 2 wide.
    gs = gridspec.GridSpec(2,4, width_ratios=(1,32,32,1))
    gs.update(left  =0.09, right =0.91,
              top   =0.9,  bottom=0.01,
              wspace=0.0,  hspace=0.0)
    # Identify axes by L LEFT, R RIGHT, T TOP, B BOTTOM, C COLOURBAR.
    LT  = plt.subplot(gs[0,1])
    LB  = plt.subplot(gs[1,1], sharex=LT, sharey=LT)
    RT  = plt.subplot(gs[0,2], sharex=LT, sharey=LT)
    RB  = plt.subplot(gs[1,2], sharex=LT, sharey=LT)
    LTC = plt.subplot(gs[0,0])
    LBC = plt.subplot(gs[1,0])
    RTC = plt.subplot(gs[0,3])
    RBC = plt.subplot(gs[1,3])
    #  ---   --- ---   ---
    # |LTC| |LT |RT | |RTC|
    #  ---   --- ---   ---
    # |LBC| |LB |RB | |RBC|
    #  ---   --- ---   ---
    # Maps and corresponding colourbars in order [LPI, I, Qphi, Uphi]:
    GRID    = [LT,RT,LB,RB,LTC,RTC,LBC,RBC]

    for X in range(4):
        GRID[X], CAX = imshow_fancy(GRID[X],MAPS[X],
                                    LIMS=LIMS,IWA=IWA,
                                    LABEL=LABELS[X],
                                    ANNOTATIONCOLOUR=COLOURS[X],BG='w')
        # Remove tick labels and marks and axis labels:
        GRID[X].set_xticks([])
        GRID[X].set_yticks([])
        GRID[X].tick_params(top=False,right=False,
                            bottom=False,left=False,which='both')
        GRID[X].set_xlabel('')
        GRID[X].set_ylabel('')
        # Make colourbar.
        CAX.set_cmap(CMAPS[X])
        CAX,CBAR = scale_colourbar(CAX,VLIMS[X],GRID[X+4],
                    VSCALE=VSCALES[X],
                    ORIENTATION='vertical',LABEL=CBARLABEL)
    # Change colourbar label position:
    LTC.yaxis.set_ticks_position('left')
    LTC.yaxis.set_label_position('left')
    LBC.yaxis.set_ticks_position('left')
    LBC.yaxis.set_label_position('left')

    # --- SAVE OR DISPLAY THE FIGURE ---
    if len(SAVENAME) > 0:
        plt.savefig(SAVENAME)
        plt.close()
    else:
        plt.show()
    return

def get_vlims(MAP,SCALE=0.001,TYPE='log'):
    """
    Get colourbar limits (vlims) for plotting MAP.

    Inputs:
    MAP       -- 2D array for which we're finding colourbar limits.
    SCALE     -- VMIN will be set to this fraction of VMAX.
    TYPE      -- Type of colourbar limits required.
                 If it's for a symlog plot, colourbar limits will be
                 calculated using absolute values and VMIN=-VMAX.

    Returns:
    VMIN,VMAX -- The calculated min/max colourbar limits.
    """
    print('Calculating vlims.')
    # Define slices for the central quarter of the image.
    # (don't use the full image to avoid high edge values)
    Q1 = int(MAP.shape[0]*0.25)
    Q2 = int(MAP.shape[1]*0.25)
    # Find vlims from this smaller map MCOPY:
    MCOPY = np.copy(MAP[Q1:-Q1,Q2:-Q2])
    MCOPY = np.abs(MCOPY) if TYPE=='symlog' else MCOPY
    # Maximum value for colour scale:
    VMAX  = np.nanmean(MCOPY) + 4.0*np.nanstd(MCOPY)
    # If VMAX is off the scale, set VMAX to max instead.
    if VMAX > np.nanmax(np.abs(MCOPY)):
        print('vmax adjusted')
        VMAX = np.nanmax(np.abs(MCOPY))
    VMIN = -VMAX if TYPE=='symlog' else SCALE*VMAX
    print('vmin, vmax:',VMIN,VMAX)
    return VMIN, VMAX


def imshow_fancy(AX, MAP, LIMS=[], IWA=0.0, LABEL='', SCALE=0.0, ANNOTATIONCOLOUR='k', AXES='pixels', XSPOTCENT=0.0, YSPOTCENT=0.0, BG='',CMAP = 'inferno'):
    """
    A fancy way to display a MAP in an axis object AX.

    This was designed with stokesdc in mind - check that it works with podc.

    Inputs:
    AX    -- the axis object where the image will be placed.
    MAP   -- the image to be displayed, e.g. a Stokes Q map.

    Optional:
    LIMS  -- axes limits [xmin,xmax,ymin,ymax] in pixel scale.
    IWA   -- inner working angle - the radius of the coronagraphic
             spot. IWA=6.0 for J-band and 8.0 for H-band.
    LABEL -- annotations for the corner of the subplot.
    SCALE -- the physical scale of 1arcsec for this map, to be
             plotted with a scale arrow in the subplot.
    AXES  -- Axis scale. 'arcsec', 'physical', or else gives pixels.
    BG    -- Background colour for region outside detector.
    CMAP  -- Colour map.

    Outputs:
    AX    -- the axis object, now containing the displayed image.
    CAX   -- cax object for the map subplot for use with
             colourbars.
    """
    # --- AXIS SCALES ---
    # Set the pixel scale. This later affects actual axis limits, axis
    # labels, and number and location of axis ticks.
    # Also set unit name for axis labels.
    if AXES == 'arcsec':
        PX      = 0.0141 #1 arcsec is 1000.0/14.14 pixels
        UNIT    = r'$\prime\prime$'
    elif AXES == 'physical':
        PX      = SCALE*0.0141
        UNIT    = 'AU'
    else:
        PX      = 1.0
        UNIT    = 'pixels'

    # --- AXIS LIMITS ---
    if len(LIMS)<1:
        # If no limits are given, calculate them here.
        L = MAP.shape[0]*0.5
        LIMS = PX*np.array([-L,L,-L,L])
    AX.set_xlim(LIMS[0],LIMS[1])
    AX.set_ylim(LIMS[2],LIMS[3])
    AX.set_xlabel(r'$\Delta \alpha $'+'('+UNIT+')')
    AX.set_ylabel(r'$\Delta \delta $'+'('+UNIT+')')

    # --- DISPLAY MAP ---
    # Find extent of image - the values at the edges of the array.
    LIMIT  = MAP.shape[0]*0.5*PX
    EXTENT = [-LIMIT,LIMIT,-LIMIT,LIMIT]
    # imshow to display the map within axes AX.
    # Keep origin='lower' to make (0,0) be in the bottom left corner.
    # Note that v limits must be updated outside this function.
    CAX = AX.imshow(MAP, extent=EXTENT,     origin='lower',
                    interpolation='none',   cmap=CMAP)

    # --- SET BAD COLOUR ---
    # Set colour of pixels outside the colourscale.
    # Specified with keyword else use minimium colour of colour map.
    CAX.cmap.set_bad(BG) if len(BG)>0 else\
    CAX.cmap.set_bad( matplotlib.cm.get_cmap(CMAP)(0.0) )

    # ----- TICKS -----
    # Add ticks to the top and right-hand sides.
    AX.tick_params(top=True,bottom=True,left=True,right=True,which='both')
    if AXES=='arcsec':
        # Edit ticks.
        # Set major ticks to every 0.5 arcsec, with labels to 1dp.
        # Set minor ticks to every 0.25 arcsec.
        # For smaller plots, change below to 0.25, 3.2f, 0.25*0.25
        # (also change SMALLARROW below).
        majorLocator   = MultipleLocator(0.5)
        majorFormatter = FormatStrFormatter('%3.1f')
        minorLocator   = MultipleLocator(0.25)
        # Apply these tick parameters to the subplot.
        AX.xaxis.set_major_locator(  majorLocator)
        AX.xaxis.set_major_formatter(majorFormatter)
        AX.xaxis.set_minor_locator(  minorLocator)
        AX.yaxis.set_major_locator(  majorLocator)
        AX.yaxis.set_major_formatter(majorFormatter)
        AX.yaxis.set_minor_locator(  minorLocator)

    if AXES=='arcsec' or AXES=='physical':
        # Reverse the x-axis tick labels for RA.
        labels = AX.get_xticks().tolist()
        # If uneven number of labels, when reversing list the 0.0 tick
        # will have shifted from 0.0. Make sure list is balanced, then
        # update labels:
        if abs(labels[0]) > abs(labels[-1]):
            labels.append(-labels[0])
            labels=labels[1:-1]
        elif abs(labels[-1]) > abs(labels[0]):
            labels.insert(0, -labels[-1])
            labels=labels[1:-1]
        ylabels = labels
        labels  = labels[::-1]
        if AXES=='physical':
            labels  = [int(i) for i in labels]
            ylabels = [int(i) for i in ylabels]
            AX.set_yticklabels(ylabels)
            # Set minor ticks - three between each pair of major ticks.
            min_majortick = int(0.5*len(labels))-1
            print(min_majortick)
            # x = 0
            # while labels[min_majortick] == 0:
            #     x+=1
            #     min_majortick = int(0.5*len(labels))-1-x
            majorLocator = MultipleLocator(     labels[min_majortick])
            minorLocator = MultipleLocator(0.25*labels[min_majortick])
            AX.xaxis.set_major_locator(majorLocator)
            AX.yaxis.set_major_locator(majorLocator)
            AX.xaxis.set_minor_locator(minorLocator)
            AX.yaxis.set_minor_locator(minorLocator)
        AX.set_xticklabels(labels)

    # ----- ANNOTATIONS -----
    # Add annotations (object name, waveband, map type).
    AX.annotate(LABEL,
                xy=(0.95*LIMS[0], 0.95*LIMS[3]),
                ha='left', va='top',
                color=ANNOTATIONCOLOUR, size=10)
    if SCALE>0.0 and UNIT!='AU':
        # Overlay arrow for scale in arcsec.
        LABEL = '%s' % int(float('%.2g' % SCALE))+'AU'
        # See here: https://stackoverflow.com/questions/38677467/
        #           how-to-annotate-a-range-of-the-x-axis-in-matplotlib
        # Arrow spatial parameters:
        SARROWRIGHT  = LIMS[1]*0.95
        SARROWLEFT   = SARROWRIGHT - 1.0
        SARROWBOTTOM = LIMS[2]*0.90
        SARROWTOP    = SARROWBOTTOM

        SMALLARROW=0
        if SMALLARROW>0:
            SARROWLEFT += 0.5
            SCALE=0.5*SCALE
            LABEL = '%s' % int(float('%.2g' % SCALE))+'AU'
        # Draw the arrow.
        # "xy" and "xytext" set the endpoints of the arrow.
        AX.annotate('',
                    xy      =(SARROWLEFT,  SARROWBOTTOM),
                    xytext  =(SARROWRIGHT, SARROWTOP),
                    xycoords='data', textcoords='data',
                    arrowprops=dict(arrowstyle= '<|-|>',
                                    color=ANNOTATIONCOLOUR,
                                    lw=1.0, ls='-'))
        # This second annotation adds the label by the arrow.
        AX.annotate(LABEL,
                    xy=(((SARROWLEFT+SARROWRIGHT)*0.5), LIMS[2]*0.8),
                    ha='center', va='center',
                    color=ANNOTATIONCOLOUR, size=10)

    # ----- SPOT -----
    # If IWA is given, draw on a coronagraphic spot.
    if IWA>0.0:
        # Create coordinates of spot, then plot as line.
        SPOTX, SPOTY = coronaspot(IWA=IWA, PX=PX)
        AX.plot(SPOTX,SPOTY,color='k')
    return AX, CAX


def scale_colourbar(CAX,VLIMS,CBAR_AXES,VSCALE='log',
                    ORIENTATION='vertical',LABEL='Surface brightness'):
    """

    # Update CAX so that the displayed map colour scale also updates.

    Inputs:
    CAX         -- Cax object from imshow.
    VLIMS       -- V limits for the colourbar.
    CBAR_AXES   -- Axis object. Place the cbar in here.
    VSCALE      -- Scaling - 'log', 'symlog', else gives linear.
    ORIENTATION -- Cbar orientation: 'vertical' or 'horizontal'.
    LABEL       -- Cbar label.

    Outputs:
    CAX         -- Updated CAX.
    CBAR        -- Colourbar in the CBAR_AXES.
    """
    # Update the colour scale.
    # --- LOG ---
    if VSCALE == 'log':
        # Log colourbar.
        CAX.set_norm(LogNorm(VLIMS[0],VLIMS[1]))
        CBAR = plt.colorbar(cax=CBAR_AXES, mappable=CAX,
                            label=LABEL,
                            orientation=ORIENTATION,
                            format=LogFormatter(10, labelOnlyBase=False))
    # --- SYMLOG ---
    elif VSCALE == 'symlog':
        if VLIMS[1] < 10.0:
            VLIMS = [-15.0,15.0]
            print('Updated symlog colour limits: ',VLIMS)
        # SymLogNorm is a log colorbar that shows values below zero.
        # (symmetrical log)
        # In SymLogNorm, linthresh is the region that is linearly mapped
        # around zero. Leave linscale as default 1.0 to have the range
        # covered by linthresh shown as one decade.
        CAX.set_norm(matplotlib.colors.SymLogNorm(linthresh=1.0,
                     linscale=1.0, vmin=-VLIMS[1], vmax=VLIMS[1]))
        # Change to a diverging colourmap for symmetry:
        CAX.set_cmap('PuOr')
        # Set bad colour for NaN to be middle colour of divering cmap.
        CAX.cmap.set_bad( matplotlib.cm.get_cmap('PuOr')(0.5) )
        # CAX.cmap.set_bad( 'w' )
        CBAR = plt.colorbar(cax=CBAR_AXES, mappable=CAX,
                            label=LABEL, orientation=ORIENTATION,
                            format = LogFormatter(10,labelOnlyBase=False))
        # Using the same formatter as for normal log results in
        # negative Uphi values being labelled as positive on the colorbar.
        # Manually set half of them to be negative with the following:
        if ORIENTATION=='horizontal':
            labels = CBAR_AXES.get_xticklabels()
        else:
            labels = CBAR_AXES.get_yticklabels()
        newlabels=[]
        for a in range(len(labels)):
            x = labels[a].get_text()
            if len(x)>0:
                # i.e. if there is a label other than empty '':
                # Convert math mode to a float:
                x = eval(x)
                # If the number is in the first half of the labels,
                # convert it to a negative value.
                if a < 0.5*len(labels):
                    x = -x
            newlabels.append(x)
        if ORIENTATION=='horizontal':
            CBAR_AXES.set_xticklabels(newlabels)
        else:
            CBAR_AXES.set_yticklabels(newlabels)
    # --- LINEAR ---
    else:
        # Linear colourbar.
        CAX.set_clim(vmin=VLIMS[0], vmax=VLIMS[1])
        CBAR = plt.colorbar(cax=CBAR_AXES, mappable=CAX,
                            label=LABEL,
                            orientation=ORIENTATION)

    if ORIENTATION == 'horizontal':
        CBAR.ax.xaxis.set_ticks_position('top')
        CBAR.ax.xaxis.set_label_position('top')
    return CAX, CBAR


def coronaspot(IWA=6.0, PX=1.0, XCENT=0.0, YCENT=0.0):
    """
    Create coordinates of the coronagraphic spot for plotting.

    Inputs:
    IWA          -- Inner Working Angle, the radius of the occulting
                    spot in pixels. Use 6 pixels for J-band, 8 pixels
                    for H-band.
    PX           -- Pixel scale. PX=1.0 for coordinates in pixels
                    or PX=0.0141 for coords in arcsec.
    XCENT,YCENT  -- Centre coordinates of the ring.

    Returns:
    SPOTX, SPOTY -- lists of X and Y coordinates to make a ring.
    """
    SPOTX = []
    SPOTY = []
    for THETA in np.linspace(0,2*pi,100):
        X = PX*(IWA*sin(THETA) + XCENT)
        Y = PX*(IWA*cos(THETA) + YCENT)
        SPOTX.append(X)
        SPOTY.append(Y)
    return SPOTX, SPOTY




# -------------------------
# DISPLAY RADIAL PROFILE
# -------------------------
def radialprofile_plot(LINES, RBINS,ERRS=[],TITLE='', NEW_FILE_NAME='n',LABELS=''):
    """
    Plot radial profiles.

    Inputs:
    LINES         -- List of profiles, typically [I,Qphi,Uphi].
    RBINS         -- Binned distances from centre pixel to match the
                     profiles in LINES.

    Optional:
    IWA           -- Value of IWA (region behind spot).
                     Draws on a shaded region up to this value.
    ERRS          -- List of errorbar-related profiles.
                     ERRS[0] - errorbar sizes.
                     ERRS[1] - the profile that the errorbars belong to.
                     ERRS[2] - radial bins for plotting the errorbars.
    TITLE         -- Optional title to print at top of figure.
    NEW_FILE_NAME -- File name to save the figure.
    LABELS        -- Labels for each profile for the legend.
    """
    plt.figure()
    plt.yscale('log')
    # Plot line graphs of variations.
    STYLES = ['-','--','-.']
    for X in range(len(LINES)):
        # LT.semilogy(RBINS,LINES[X],color='k',label=LABELS[X],linestyle=STYLES[X])
        plt.plot(RBINS,LINES[X],color='k',label=LABELS[X],linestyle=STYLES[X])
        if len(ERRS)>0:
            STDS = np.array(ERRS[0][X][::2])
            ERR_QS = np.array(ERRS[1][X][::2])
            ERR_RBINS = np.array(ERRS[2][::2])
            # Errorbar +/- are relative to the point, so subtract the
            # main data points.
            MAXS = ERR_QS+STDS
            MINS = ERR_QS-STDS
            plt.fill_between(ERR_RBINS, MINS, MAXS,
                        alpha=0.5, edgecolor='gray', facecolor='Gainsboro')
    # Setup:
    plt.legend(loc='upper right')
    plt.xlabel('Distance')
    plt.ylabel('Flux')
    plt.title(TITLE)
    if len(NEW_FILE_NAME)>2:
        fig = plt.gcf()
        fig.savefig(NEW_FILE_NAME)
    else:
        plt.show()
    return


# -------------------------
# FLUXCUT
# -------------------------
def plotfluxcut(SFLUXCUT, CFLUXCUT, SMEAN, SSTD, CMEAN, CSTD, TITLE=''):
    """
    Compare how two rows of flux (SFLUXCUT,CFLUXCUT) vary across the map.

    Use this to compare the same row of fluxes as it appears in a stokesdc
    made from 32 podc files, and in a stokesdc made by combining eight
    intermediary stokesdc files.

    Plots the changing flux FLUXCUT across some row of pixels, one set each
    for combined stokesdc and single stokesdc. Also shows the mean flux and
    the standard deviation of each data set.
    """
    # Find the X values to plot against:
    # (assumes that the full row is plotted, i.e. x-values 0 to 280)
    X = range(len(SFLUXCUT))
    # Add dashed lines to mark the mean values of the fluxes:
    plt.axhline(SMEAN, linestyle='-.', color='LightCoral')
                #, label='single mean: '+str(round(s_mean,2)))
    plt.axhline(CMEAN, linestyle=':',  color='CornflowerBlue')
                #, label='combined mean: '+str(round(c_mean,2)))
    # Add a coloured background in the region covered by mean +/- std:
    plt.axhspan(SMEAN-SSTD, SMEAN+SSTD, alpha=0.5, color='LightCoral')
    plt.axhspan(CMEAN-CSTD, CMEAN+CSTD, alpha=1.0, color='LightSteelBlue')
    # Now plot the actual data on top:
    plt.plot(X,SFLUXCUT, color='r', label='Single stokesdc')
    plt.plot(X,CFLUXCUT, color='b', label='Combined stokesdc')
    # Plot setup:
    plt.xlabel('Pixel number in x')
    plt.ylabel('Flux (counts)')
    plt.title(TITLE)
    plt.legend(loc='upper left')
    plt.show()
    plt.close()
    return


#-----------------------------
# Plot flexure
#-----------------------------
def plot_flexure(image,loc,name,sdx,sdy,reduced_dir):
    """
    Plot the measured flexure measured and found in the headers of the podc files onto the raw files
    """
    c = (1920,1335) #(1915,1339)
    du = 25
    fig = plt.figure(figsize=(6,6))
    img1 = image[int(c[1])-du:int(c[1])+du,int(c[0])-du:int(c[0])+du]
    vmin, vmax = np.median(img1)*0.01, np.median(img1)*10.
    plt.imshow(img1, origin='lower', vmin=vmin, vmax=vmax)
    for i in range(0,len(loc[0,:,0,0])):
        for j in range(0,len(loc[0,0,:,0])):
            if str(loc[0,i,j,0]) != 'nan':
                x, y = loc[0,i,j,0]+sdx, loc[0,i,j,1]+sdy
                if np.abs(x - c[0]) < du-1 and np.abs(y - c[1]) < du-1:
                    plt.plot(x-int(c[0])+du,y-int(c[1])+du,'x',color='w',ms=8,mew=3)
                    x2, y2 = loc[1,i,j,0]+sdx, loc[1,i,j,1]+sdy
                    if np.abs(x2 - c[0]) < du-1 and np.abs(y2 - c[1]) < du-1:
                        plt.plot(x2-int(c[0])+du,y2-int(c[1])+du,'x',color='r',ms=3,mew=2)
                        dx, dy = -int(c[0])+du, -int(c[1])+du
                        plt.plot((x +dx,x2 +dx),(y+dy,y2+dy),'--',c='r')

    plt.title(name + '_flexure.png')
    plt.savefig(reduced_dir + name + '_flexure.png')
    plt.close()
    return

#-----------------------------------------
# Gather flexure corrections for plotting
#----------------------------------------
def gather_flexure(raw_files,rawdata_dir,reduced_dir):
    podc_file = raw_files[0][:-5] + '_podc.fits'
    hold = getfitskeywords(reduced_dir+podc_file, 'HISTORY')
    print(hold)
    #hold = str(header['HISTORY']).split('\n')
    seg = 'GPI_LOAD_POLARIMETRY_SPOT_CALIBRATION: /Users/'
    for i in range(0,len(hold)):
        line = hold[i]
        loc = line.find(seg)
        if loc != -1:
            if line.find('.fits') == -1:
                calib = line+hold[i+1]
                if calib[-5:]!='.fits':
                    extra = hold[i+2].split(' ')[0]
                    calib+=extra
                if calib[-5:]!='.fits': #if the above didn't fix it...
                    print('there is a problem in gather_flexure')
            else:
                calib = line
            calib = calib.split(':')[-1].strip(' ')
    print(calib)
    hdu = fits.open(calib) #('/Users/earich/work_reduction2/Reduced/calibrations/S20171231S0432_J-polcal.fits')
    loc = hdu[1].data
    for i in range(0,len(raw_files)):
        podc_file = raw_files[i][:-5] + '_podc.fits'
        images = getfitsdata(rawdata_dir+raw_files[i])
        SPOT_DX = getfitskeywords(reduced_dir+podc_file, 'SPOT_DX')
        SPOT_DY = getfitskeywords(reduced_dir+podc_file, 'SPOT_DY')
        plot_flexure(images,loc,raw_files[i].split('/')[-1][:-5],SPOT_DX,SPOT_DY,reduced_dir)
    return

#--------------------------------------
# Plot change in Flux distribution from cycle to cycle
#-----------------------------------------
def plot_fluxchange(files,directory,Bcent=(-1,-1), du = 20):
    if Bcent != (-1,-1):
        hdu = fits.open(directory + files[0].strip('\n')[:-5] + '_podc.fits')
        header = hdu[1].header
        PA, PSFCENTX, PSFCENTY = header['AVPARANG'],header['PSFCENTX'],header['PSFCENTY']
        PA_first = PA
        theta = np.arctan2(Bcent[1] - PSFCENTY,Bcent[0] - PSFCENTX)
        r = np.sqrt((PSFCENTX-Bcent[0])**2. + (PSFCENTY-Bcent[1])**2.)
    for i in range(0,len(files)-4):
        file1, file2 = files[i].strip('\n')[:-5] + '_podc.fits', files[i+4].strip('\n')[:-5] + '_podc.fits'
        hdu1 = fits.open(directory + file1)
        data1 = hdu1[1].data
        hdu2 = fits.open(directory + file2)
        data2 = hdu2[1].data
        header11 = hdu1[1].header
        header10 = hdu1[0].header
        RMSERR = header10['RMSERR']
        if Bcent != (-1,-1):
            PA, PSFCENTX, PSFCENTY = header11['AVPARANG'],header11['PSFCENTX'],header11['PSFCENTY']
            ang = theta - (PA_first-PA)*np.pi/180.
            x2,y2 = r*np.cos(ang)+PSFCENTX, r*np.sin(ang)+PSFCENTY
            data1[0][int(y2)-du:int(y2)+du,int(x2)-du:int(x2)+du] = np.nan
            data2[0][int(y2)-du:int(y2)+du,int(x2)-du:int(x2)+du] = np.nan
        hdu1.close()
        hdu2.close()
        fig = plt.figure(figsize=(8,8))
        gs1 = gridspec.GridSpec(2,2)
        gs1.update(wspace=0.25, hspace=0.25)
        ax1 = plt.subplot(gs1[0])
        ax3 = plt.subplot(gs1[2])
        ax4 = plt.subplot(gs1[3])
        ax1.title.set_text('frame ' + str(i) + ' - frame ' + str(i+4))
        img = ax1.imshow((data1[0]-data2[0])/data2[0], origin='lower', vmin = -0.3, vmax=0.3)
        divider = make_axes_locatable(ax1)
        cax = divider.append_axes('right', size='5%', pad=0.05)
        fig.colorbar(img, cax=cax, orientation='vertical')
        ax3.title.set_text('cycle ' + str(i//4) + ' frame ' + str(i))
        ax4.title.set_text('cycle ' + str((i+4)//4) + ' frame ' + str(i+4))
        residual = np.nansum(np.abs(data1[0]-data2[0])/data2[0])
        ax1.text(140,25,'Resid. Sum: ' + str(int(round(residual,0))), Fontsize = 10)
        ax1.text(140,10,'RMS W. Err: ' + str(int(RMSERR)), Fontsize = 10)
        ax1.text(10,250,'cycle: ' + str(i//4) + ' - cycle: ' + str((i+4)//4))
                
        img3 = ax3.imshow(data1[0], origin='lower',#CMAP='Blues_r',
               norm=matplotlib.colors.SymLogNorm(linthresh=10.0, 
               linscale=1.0, vmin=50.0, vmax=10000.0)) # vmin = 0.0, vmax=10000.)
        divider = make_axes_locatable(ax3)
        cax3 = divider.append_axes('right', size='5%', pad=0.05)
        fig.colorbar(img3, cax=cax3, orientation='vertical')
        img4 = ax4.imshow(data2[0], origin='lower', #CMAP='Blues_r',
               norm=matplotlib.colors.SymLogNorm(linthresh=10.0, 
               linscale=1.0, vmin=50.0, vmax=10000.0)) #vmin = 0.0, vmax=10000.)
        divider = make_axes_locatable(ax4)
        cax4 = divider.append_axes('right', size='5%', pad=0.05)
        fig.colorbar(img4, cax=cax4, orientation='vertical')
        
        plt.savefig(directory + file1.split('_')[0] + '_' + file2.split('_')[0] + '_residual.png')
        print(directory + file1.split('_')[0] + '_' + file2.split('_')[0] + '_residual.png')
        plt.close()
    return

""" #Depreciated. If above function works, needs to be deleted.
def plot_fluxchange(files,directory,Bcent=(-1,-1), du = 20):
    if Bcent != (-1,-1):
        hdu = fits.open(directory + files[0].strip('\n')[:-5] + '_podc.fits')
        header = hdu[1].header
        PA, PSFCENTX, PSFCENTY = header['AVPARANG'],header['PSFCENTX'],header['PSFCENTY']
        PA_first = PA
        theta = np.arctan2(Bcent[1] - PSFCENTY,Bcent[0] - PSFCENTX)
        r = np.sqrt((PSFCENTX-Bcent[0])**2. + (PSFCENTY-Bcent[1])**2.)
    for i in range(0,len(files)-4):
        file1, file2 = files[i].strip('\n')[:-5] + '_podc.fits', files[i+4].strip('\n')[:-5] + '_podc.fits'
        hdu1 = fits.open(directory + file1)
        data1 = hdu1[1].data
        hdu2 = fits.open(directory + file2)
        data2 = hdu2[1].data
        header11 = hdu1[1].header
        header10 = hdu1[0].header
        RMSERR = header10['RMSERR']
        if Bcent != (-1,-1):
            PA, PSFCENTX, PSFCENTY = header11['AVPARANG'],header11['PSFCENTX'],header11['PSFCENTY']
            ang = theta - (PA_first-PA)*np.pi/180.
            x2,y2 = r*np.cos(ang)+PSFCENTX, r*np.sin(ang)+PSFCENTY
            data1[0][int(y2)-du:int(y2)+du,int(x2)-du:int(x2)+du] = np.nan
            data2[0][int(y2)-du:int(y2)+du,int(x2)-du:int(x2)+du] = np.nan
        hdu1.close()
        hdu2.close()
        fig = plt.figure(figsize=(8,8))
        img = plt.imshow((data1[0]-data2[0])/data2[0], origin='lower', vmin = -0.3, vmax=0.3)
        plt.colorbar()
        plt.title('(' + file1.split('_')[0] + ' - ' + file2.split('_')[0] + ')/ ' + file2.split('_')[0])
        residual = np.nansum(np.abs(data1[0]-data2[0])/data2[0])
        plt.text(165,20,'Residual Sum: ' + str(round(residual,1)))
        plt.text(165,10,'RMS Wavefront Error' + str(RMSERR))
        plt.text(10,250,'cycle: ' + str(i//4) + ' - cycle: ' + str((i+4)//4))
        plt.savefig(directory + file1.split('_')[0] + '_' + file2.split('_')[0] + '_residual.png')
        print(directory + file1.split('_')[0] + '_' + file2.split('_')[0] + '_residual.png')
        plt.close()
    return
"""
