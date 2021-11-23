"""
Calculate distances between PSF centres of podc to find outliers.
"""
import numpy as np
import os
from gpi_analysis.inputs import getfitskeywords, getfilenames, getfilenames_l, getfitsdata
from astropy.io import fits
import matplotlib.pyplot as plt

def get_psfcoords(FITS_NAME):
    """
    Get coords from podc files:
    """
    XCENT = getfitskeywords(FITS_NAME,'PSFCENTX',HEADER='SCI')
    YCENT = getfitskeywords(FITS_NAME,'PSFCENTY',HEADER='SCI')
    return np.array([XCENT,YCENT])

def updatefitsheader(FILENAME, KEYS, VALUES, HEADER='PRIMARY', NEW_FILE_NAME=''):
    """
    Save a new fits of FILENAME, overwriting KEYS in HEADER with VALUES.

    Example usage: to overwrite PSFCENTX, PSFCENTY with calculated values.
    Each element of VALUES is a tuple (new value, comment about old value).
    """
    if len(NEW_FILE_NAME)<1:
        NEW_FILE_NAME=FILENAME
    # Open the fits file.
    HDUL = fits.open(FILENAME)
    HDR = HDUL[HEADER].header
    for K in range(len(KEYS)):
        HDR[KEYS[K]] = VALUES[K]
    HDUL2 = HDUL
    # Save a new file with the new data but old headers etc.
    HDUL2.writeto(NEW_FILE_NAME)
    # Close the fits file again.
    # Anything that happens before this close statement will be saved.
    # The old file won't be overwritten unless explicitly told to
    # (.e.g with writeto())
    HDUL.close()
    print('Saved file: '+NEW_FILE_NAME)
    return

def dists_array(coords):
    """
    Get distances from each array element to each other.
    """
    # Cycle through y and x to get n arrays each with n distances.
    # e.g. four values, four arrays of four distances each.
    # Each distance appears in two arrays (n1->n2 and n2->n1).
    y=1
    x=0
    dists = [[]]
    while x < len(coords):
        dist = np.linalg.norm(coords[x]-coords[y])
        dists[x].append(dist)
        y+=1
        if y > len(coords)-1:
            y=0
        if y==x:
            x+=1
            y+=2
            if x < len(coords):
                dists.append([])
        if y > len(coords)-1:
            y=0
    lencheck = [len(d) for d in dists]
    if np.mean(lencheck) != lencheck[0]:
        print('maybe weird list lengths? Check.')
        print(lencheck)
    return dists

# def anom_coords(COORDS, FILES):
#     """
#     Outdated - use mean and std to decide if value is anomalous.
#     """
#     print('Initial coords:')
#     for a in range(len(FILES)):
#         print(FILES[a].split('/')[-1], COORDS[a])
#
#     DISTS = dists_array(COORDS)
#     MEANS = [np.mean(X) for X in DISTS]
#     MEAN  = np.mean(MEANS)
#     STD   = np.std(MEANS)
#     print('Mean dists: ', np.round(MEANS,2))
#     print('Mean: ',round(MEAN,2), 'std: ',round(STD,2))
#     # FIND POINTS WHERE MEAN DISTANCE IS GREATER THAN MEAN+/-STD.
#     # ANOM = np.where( abs(MEANS-MEAN) > abs(STD))[0]
#     # FINE = np.where( abs(MEANS-MEAN) <= abs(STD))[0]
#     ANOM = np.where( abs(MEANS-MEAN) > 1.0)[0]
#     FINE = np.where( abs(MEANS-MEAN) <= 1.0)[0]
#     print(len(ANOM), 'anomalous coords:')
#     if len(ANOM)>0:
#         for A in ANOM:
#             print(FILES[A].split('/')[-1], COORDS[A])
#         X = np.mean(COORDS[FINE,0])
#         Y = np.mean(COORDS[FINE,1])
#         print('Updated to coords:      ', [round(X,3),round(Y,3)])
#         for I in ANOM:
#             COORDS[ANOM] = [round(X,3),round(Y,3)]
#     print('Final coords: ')
#     for a in range(len(FILES)):
#         print(FILES[a].split('/')[-1], COORDS[a])
#     print('')
#     return COORDS

def far_check(DISTS):
    """
    Check if coord is more than one pixel (dist=1.41) away from all others.
    """
    # Indices of good/bad coords:
    GOODS    = []
    BADS     = []
    for D in range(len(DISTS)):
        FARS = np.where( np.array(DISTS[D]) > 2**0.5 )[0]
        # Mark BAD if it's over the distance from all other points.
        if len(FARS) > len(DISTS[D])-1:
            BADS.append(D)
        else:
            GOODS.append(D)
    return GOODS, BADS

def far_coords(COORDS, FILES, TEXT=''):
    """
    Return updated coordinates with anomalies replaced with averages.
    """
    TEXT+='Initial coords:'
    TEXT+='\n'
    for A in range(len(FILES)):
        TEXT+=FILES[A].split('/')[-1]+ ' ' + str(COORDS[A])
        TEXT+='\n'
    # Find distances from each point to each other:
    DISTS            = dists_array(COORDS)
    # Find locations of "bad" anomalous values to be replaced.
    GOODS, BADS      = far_check(DISTS)
    TEXT+=str(len(BADS))+' anomalous coords.'
    TEXT+='\n'
    for A in BADS:
        # Overwrite "bad" coords with average of "good" coords.
        TEXT+=FILES[A].split('/')[-1]+ ' ' + str(COORDS[A])
        TEXT+='\n'
        XMEAN        = np.mean(COORDS[GOODS][:,0])
        YMEAN        = np.mean(COORDS[GOODS][:,1])
        COORDS[BADS] = [round(XMEAN,3),round(YMEAN,3)]
        TEXT+='Updated to coords:       '+ str([round(XMEAN,3),round(YMEAN,3)])
        TEXT+='\n'
    if len(BADS)>0:
        TEXT+='Final coords: '
        TEXT+='\n'
        for A in range(len(FILES)):
            TEXT+=FILES[A].split('/')[-1]+ ' ' + str(COORDS[A])
            TEXT+='\n'
    TEXT+='\n'
    return COORDS, TEXT

def plot_centre(coords, files,dir):
    for i in range(0,len(files)):
        image = getfitsdata(files[i])
        print(files[i],coords[i][0],coords[i][1])
        x,y = int(coords[i][0]), int(coords[i][1])
        ds = 25
        dx, dy = x - ds, y - ds
        fig = plt.figure()
        img1 = plt.imshow(image[0,y-ds:y+ds,x-ds:x+ds],origin='lower')
        for A in coords:
            plt.plot(A[0]-dx,A[1]-dy,'x',color='w')
        plt.plot(coords[i][0]-dx,coords[i][1]-dy,'o',color='w')
        plt.savefig(files[i][:-5] + '_centers.png')

def main(dir):
    text=''
    text+=dir+'\n'
    #files       = getfilenames(dir, END='podc.fits')
    files       = getfilenames_l(dir, 'files.lst', END='_podc.fits')
    coords      = []
    for p in files:
        psfcent = get_psfcoords(p)
        coords.append(psfcent)
    coords      = np.array(coords)
    #
    # Run coords check in cycles of four podc:
    i=1
    new_coords  = np.empty(coords.shape)
    for x in range(0,len(files),4):
        text+='================ Cycle '+str(i)+' ================'+'\n'
        new_coords[x:x+4],text = far_coords(np.copy(coords[x:x+4]), files[x:x+4],text)
        i+=1
        #print(coords[x:x+4],new_coords[x:x+4])
        plot_centre(new_coords[x:x+4],files[x:x+4],dir)
    #
    # Save text file of useful output:
    text_file_name = '/'.join(files[0].split('/')[:-1])+'/'+\
                     'PSFcheck.txt'
    file = open(text_file_name, 'w')
    file.write(text)
    file.close()
    print(text)

    # Save new podc with updated coords in header:
    keys        = ['PSFCENTX', 'PSFCENTY']
    for a in range(len(files)):
        values  = []
        for k in range(len(keys)):
            value      = (new_coords[a][k],
                          'Pipeline-calculated value: '+str(coords[a][k]) )
            values.append(value)
        old_file_name  = files[a].split('/')[-1]
        new_file_name  = '/'.join(files[a].split('/')[:-1])+'/'+\
                         old_file_name.split('_')[0]+\
                         '_PSFchecked_'+\
                         old_file_name.split('_')[-1]
        updatefitsheader(files[a], keys, values, HEADER='SCI',
                         NEW_FILE_NAME=new_file_name)
    return



# %%%%%%%%%% MAIN %%%%%%%%%%%%%
if __name__ == '__main__':
    dir = '/Users/earich/work_reduction/Reduced/190715/FUOri-J/' #'/Users/al630/Documents/GPI/Results/Minimisation tests/stokesdc/pipeline correction/FU Ori/'#'/Users/al630/Documents/GPI/Results/Minimisation tests/stokesdc/pipeline correction/FU Ori/'
    main(dir)
