"""
Demo of some functions that are useful outside the reduction pipeline.

Written for python 3.
"""
import numpy as np
import matplotlib.pyplot as plt

# For fancier plot:
import matplotlib.gridspec as gridspec
from matplotlib.colors import LogNorm
import matplotlib.cm

# Custom functions from our pypeline:
from gpi_analysis.inputs import getfitsdata, getfitskeywords
from gpi_analysis.analysis import make_radialstokes, make_linpolint
from gpi_analysis.plot import imshow_fancy, get_vlims, scale_colourbar


# Import an rstokesdc:
filename = '/Users/al630/Documents/GPI/GPI LLP data paper/Reduction details/stokesdc/Hen3-365-J_S20170406S0176_combined_rstokesdc.fits'
i,qphi,uphi,v = getfitsdata(filename)

# # Alternative:
# # Import a stokesdc:
# f_stokes = '/Users/al630/Documents/GPI/GPI LLP data paper/Reduction details/stokesdc/Hen3-365-J_S20170406S0176_combined_stokesdc.fits'
# i,q,u,v = getfitsdata(f_stokes)
# # n.b. for Evan's reduced files, this q,u, is equivalent to Q*,U* in Laws et al. 2020.
# # for files straight out of the unmodified DRP, this q,u is Q,U in Laws+2020.
# qphi, uphi = make_radialstokes(q,u)
# lpi = make_linpolint(q,u)

# Get useful keywords from fits header:
target  = getfitskeywords(filename, 'OBJECT')
itime   = getfitskeywords(filename, 'ITIME', HEADER='SCI')
print('target, itime', target,itime)

# Flux conversion:
# rstokesdc has qphi,uphi in ADU/coadd. Convert to mJy/arcsec^2.
# fluxscale is different for every observation,
# even for the same object in the same band observed on different days.
fluxscale = 2.6 # units of {mJy arcsec−2}/ {(ADU/coadd)/s}
i,qphi,uphi,v = np.array([i,qphi,uphi,v])*fluxscale/itime

exit()

# Some plot examples:
# # Plot 1 - Basic plot:
plt.imshow(np.log(qphi), origin='lower')
plt.show()
# # n.b. causes gaps and ugly images where values are 0 or less.
# exit()


# Plot 2 - basic plot, but less ugly:
cax = plt.imshow( qphi, origin='lower')
# This is a quick way to get sensible colour scale limits:
vlims = get_vlims(qphi,TYPE='log')
# Update the colour scale.
cax.set_norm(LogNorm(vlims[0],vlims[1]))
# Set values outside colour scale to something other than white:
cax.cmap.set_bad( matplotlib.cm.get_cmap('viridis')(0.0) )
plt.colorbar(label='log label')
plt.show()
# exit()


# Plot 3 - symmetrical log scale (e.g. for Uphi):
cax = plt.imshow(uphi, origin='lower')
# Update colour scale limits:
vlims = get_vlims(uphi,TYPE='symlog')
# Set to symlog scale:
cax.set_norm(matplotlib.colors.SymLogNorm(linthresh=1.0,
             linscale=1.0, vmin=-vlims[1], vmax=vlims[1]))
# Change to a diverging colourmap for symmetry:
cax.set_cmap('PuOr')
# Set bad colour for NaN to be middle colour of divering cmap.
cax.cmap.set_bad( matplotlib.cm.get_cmap('PuOr')(0.5) )
plt.colorbar(label='symlog label')
plt.show()
# exit()


# Plot 4 - more fancy plot:
# imshow_fancy() is a bit of a Frankenstein. It has lots of
# random small details that have been added over the years.
# Some of it might be useful to copy elsewhere,
# even if not calling the full function.
plt.figure()#figsize=(3.0,3.0))
gs = gridspec.GridSpec(1,2,width_ratios=[32,1],wspace=0.0,
    left=0.18,right=0.8,top=0.95, bottom=0.15)
# the left,right etc. are for whitespace around the edge.
# These specific numbers are arbitrary.
gs0 = plt.subplot(gs[0,0])
gs1 = plt.subplot(gs[0,1])

gs0,cax = imshow_fancy(gs0,qphi,CMAP='Reds',AXES='arcsec')
gs0.annotate('an annotation',xy=(-0.6,0.58),va='top')
vlims = get_vlims(qphi,TYPE='log')
cax, cbar = scale_colourbar(cax,vlims,CBAR_AXES=gs1,
                            VSCALE='log',
                            LABEL='cbar label')
plt.show()
# exit()


# Plot 5 - detailed imshow_fancy:
plt.figure()#figsize=(3.0,3.0))
gs = gridspec.GridSpec(2,1,height_ratios=[1,32],wspace=0.0,
    left=0.18,right=0.99,top=0.9, bottom=0.15)
gs0 = plt.subplot(gs[1,0])
gs1 = plt.subplot(gs[0,0])

gs0,cax = imshow_fancy(
    AX=gs0, # axes to plot on
    MAP=qphi, # the map (numpy array) to plot
    # LIMS=[-1.2,1.2,-0.8,0.8], # axis limits xmin,xmax,ymin,ymax
    IWA=8.0, # draws on a circle in the image centre (mark coronagraph). Units of pixels.
    LABEL='some \nwords', # adds an annotation
    SCALE=2000.0, # physical scale - how many AU are in 1 arcsec. Adds a scale arrow to the image.
    ANNOTATIONCOLOUR='w', # colour of labels
    AXES='arcsec', # axis scale. Can be 'arcsec' or else for pixels. ('physical' is broken today.)
    # XSPOTCENT=0.0, # centre x coordinate of IWA circle
    # YSPOTCENT=0.0, # centre y coordinate of IWA circle
    BG='LimeGreen', # set background colour (default is based on colour map)
    CMAP = 'inferno', # colour map.
    # FONTSIZE=10, # font size for axes, labels etc.
    PX=0.0141 # pixel scale - for GPI, 1 pixel is 0.0141 arcsec.
    )

vlims = get_vlims(qphi,TYPE='log')
cax, cbar = scale_colourbar(cax,vlims,CBAR_AXES=gs1,
                            VSCALE='log',
                            LABEL='cbar label',FONTSIZE=10,
                            ORIENTATION='horizontal')
plt.show()
