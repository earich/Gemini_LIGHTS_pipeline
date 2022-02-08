# coding=utf-8
"""
GPI analysis functions related to im/exporting data from fits files.
"""
from astropy.io import fits
import numpy as np
import os


# -------------------------
# CHECK FILES
# -------------------------
def getfilenames(DIRECTORY, END='.fits', KEEPDIR=1):
    """
    Originally in plotfourtest - find permanent home
    """
    NAMES     = os.listdir(DIRECTORY)
    LEN       = len(END)
    if KEEPDIR==1:
        FITSNAMES = [ DIRECTORY+A for A in NAMES if A[-LEN:]==END ]
    else:
        FITSNAMES = [ A for A in NAMES if A[-LEN:]==END ]
    return FITSNAMES

# -------------------------
# return Filenames listed in a parameter file
# -------------------------
def getfilenames_l(DIRECTORY,filename, KEEPDIR=1, END='.fits'):
    F = open(DIRECTORY + filename) #open command automatically adds '\' after spaces!!!!!
    #if DIRECTORY.find('\ ') == -1:
    #    DIRECTORY = DIRECTORY.replace(' ','\ ')
    lines = F.readlines()
    F.close()
    list1 = []
    if KEEPDIR==1:
        for line in lines:
            line = DIRECTORY + line.strip('\n')[:-5] + END
            list1.append(line)
        return list1
    else:
        for line in lines:
            line = line.strip('\n')[:-5] + END
            list1.append(line)
        return list1


# -------------------------
# IMPORTING DATA
# -------------------------
def getfitsdata(FILENAME):
    """
    Get the science data out of the fits file FILENAME.

    For input stokesdc, the output is a list with maps [I,Q,U,V].
    For input podc, the output is a list with maps [XY,YY].
    """
    # Open the fits file.
    HDUL = fits.open(FILENAME)
    # First extension hdul[1] has the science data. 'SCI' extension.
    # For stokesdc, SCI_DATA[0] is I, [1] is Q, [2] U, [3] V.
    SCI_DATA = HDUL['SCI'].data
    # Close the fits file again.
    HDUL.close()
    return SCI_DATA


def getfitskeywords(FITS_NAME, KEYWORD, HEADER='PRIMARY'):
    """
    Returns the value of some KEYWORD in the fits header of FITS_NAME.

    Example keywords: 'WPANGLE' - HWP angle, 'OBJECT' - target name.
    Can use this for quickly viewing the HWP angles of multiple podc.
    """
    if KEYWORD == 'OBJECT':
        HDUL = fits.open(FITS_NAME)
        # Primary HDU hdul[0]
        # First extension hdul[1] has the science data. 'SCI' extension.
        HDUL[HEADER].header.set('OBJECT',HDUL[HEADER].header[KEYWORD].replace(' ',''))
        KEYWORD_VALUE = HDUL[HEADER].header[KEYWORD]
        HDUL.writeto(FITS_NAME,overwrite=True)
        HDUL.close()
    else:
        HDUL = fits.open(FITS_NAME)
        # Primary HDU hdul[0]
        # First extension hdul[1] has the science data. 'SCI' extension.
        KEYWORD_VALUE = HDUL[HEADER].header[KEYWORD]
        HDUL.close()
    return KEYWORD_VALUE


# -------------------------
# EXPORTING DATA
# -------------------------
def savefitsdata(FILENAME, NEW_DATA, NEW_FILE_NAME):
    """
    Overwrite a fits file FILENAME with NEW_DATA, keeping old headers.

    Example usage: to save a fits version of a combined stokesdc map,
    with NEWDATA as a list of maps [I, Q, U, V].
    """
    # Open the fits file.
    HDUL = fits.open(FILENAME)
    # First extension hdul[1] has the science data. 'SCI' extension.
    # For stokesdc, SCI_DATA[0] is I, [1] is Q, [2] U, [3] V.
    HDUL['SCI'].data = NEW_DATA
    HDUL2 = HDUL#['SCI']   <-- check out why I had 'SCI' here...
    # Save a new file with the new data but old headers etc.
    HDUL2.writeto(NEW_FILE_NAME, overwrite = True)
    # Close the fits file again.
    # Anything that happens before this close statement will be saved.
    # The old file won't be overwritten unless explicitly told to
    # (.e.g with writeto())
    HDUL.close()
    print('Saved file: '+NEW_FILE_NAME)
    return

# -------------------------
# CREATE RECIPES
# -------------------------
def createrecipe(INPUT_FILES, INPUT_DIRECTORY, REDUCED_DIRECTORY,
                 QUEUE_DIRECTORY, KEYS=[], RECIPE=0, PRINT_USER_INFO='', RETURN_OUTPUT_DIRECTORY=0):
    """
    Create a recipe for the GPI pipeline and pass it into the DRP queue.

    RECIPE          - Set to 0 for basicpolcube (to create podc files),
                      or to 1 for polsequencefromcubes (create stokesdc
                      from podc).
    INPUT_FILES     - The input files. Can be podc or stokesdc.
    INPUT_DIRECTORY - Directory containing input files. Include path.
    QUEUE_DIRECTORY - Directory for queue. Should be in data/queue,
                      near data/reduced directory with results.
    PRINT_USER_INFO - To bypass printing user info (username, pc name)
                      in the recipe file.

    Returns:
    FILE_OUTPUT_DIRECTORY - This output directory will be the input
                            directory for later stokesdc creation recipes.
                            Export it here to use again later.
    """
    # Locate the relevant recipe template and define the output recipe
    # file name. Keep 'waiting' in output file name so the recipe will
    # run immediately when placed in queue.
    DIR = os.path.dirname(os.path.abspath(__file__))
    if RECIPE == 0:
        # podc creation
        TYPE = 'podc'
        RECIPE_TEMPLATE = 'basicpolcube_nocent_template.txt'
        RECIPE_TYPE = '_basicpolcube_recipe.waiting.xml'
    elif RECIPE == 1:
        # podc creation
        TYPE = 'podc'
        RECIPE_TEMPLATE = 'basicpolcube_onlycent_template.txt'
        RECIPE_TYPE = '_basicpolcube_recipe.waiting.xml'
    elif RECIPE == 2:
        # stokesdc creation from podc
        TYPE = 'stokesdc'
        RECIPE_TEMPLATE = 'polsequencefromcubes_template.txt'
        RECIPE_TYPE = '_polsequencefromcubes_recipe.waiting.xml'
        FILE_OUTPUT_DIRECTORY = INPUT_DIRECTORY
    elif RECIPE == 3:
        # combine/align for flux cal
        TYPE = 'combine/align'
        RECIPE_TEMPLATE = 'combine_align_3dcube_template.txt'
        RECIPE_TYPE = '_combinealign_recipe.waiting.xml'
        FILE_OUTPUT_DIRECTORY = INPUT_DIRECTORY
    elif RECIPE == 4:
        # flux calibration
        TYPE = 'fluxcal'
        RECIPE_TEMPLATE = 'calibrate_flux_pol_template.txt'
        RECIPE_TYPE = '_polcalibrateflux_recipe.waiting.xml'
        FILE_OUTPUT_DIRECTORY = INPUT_DIRECTORY
    elif RECIPE == 5:
        # podc creation
        TYPE = 'podc'
        RECIPE_TEMPLATE = 'basicpolcube_template.txt'
        RECIPE_TYPE = '_basicpolcube_recipe.waiting.xml'
    else:
        print('Wrong recipe attempted')
        exit()
    RECIPE_TEMPLATE = DIR + '/'+RECIPE_TEMPLATE

    # Copy all the text from the existing recipe file.
    f = open(RECIPE_TEMPLATE, 'r')
    TEMPLATE_FULL_STRING = f.read()
    f.close()
    # TEMPLATE_SPLIT_STRINGS[0] is all the data before the inserted lines.
    # TEMPLATE_SPLIT_STRINGS[1] is all the data following it.
    TEMPLATE_SPLIT_STRINGS = TEMPLATE_FULL_STRING.split('<!--')

    # GET TARGET NAME
    TARGET_NAME = getfitskeywords(INPUT_DIRECTORY+INPUT_FILES[-1], 'OBJECT')
    DATE_OBS = getfitskeywords(INPUT_DIRECTORY+INPUT_FILES[-1], 'DATE-OBS')
    DATE_OBS = ''.join(DATE_OBS.split('-'))
    WAVEBAND    = getfitskeywords(INPUT_DIRECTORY+INPUT_FILES[-1], 'OBSMODE')[0]

    # ----------
    # GET USER INFO STRING
    # Get current date and time in UTC. Used to define output directory.
    from datetime import datetime
    DATE      = str(datetime.utcnow()) # ='2011-05-03 17:45:35.177000'
    SHORTDATE = DATE[2:10][0:2]+DATE[2:10][3:5]+DATE[2:10][6:8] #='110503'
    # I anticipate this causing compatibility issues, so made it optional:
    if len(PRINT_USER_INFO)<1:
        import getpass
        import socket
        # Get username
        USERNAME = getpass.getuser()
        # Get computer name
        PCNAME = socket.gethostname()
        # Example full author string:
        #<!-- recipe written by al630 on peter.lan at 2018-05-23T22:03:23 UTC -->
        AUTHOR_STRING = '<!-- recipe written by ' + USERNAME + ' on ' +\
                       PCNAME + ' at ' + DATE + ' UTC -->'
    else:
        AUTHOR_STRING = '<!-- recipe written by ? on ? at ' + DATE + ' UTC -->'
    
    # Return the file output directory.
    # The podc output directory is also the stokesdc input directory.
    #FILE_OUTPUT_DIRECTORY = REDUCED_DIRECTORY + SHORTDATE + '/' +\
    #                        TARGET_NAME + '-' + WAVEBAND + '/'
    
    if TYPE == 'podc':
        FILE_OUTPUT_DIRECTORY = REDUCED_DIRECTORY + TARGET_NAME + '-' + WAVEBAND + '/' +\
                                DATE_OBS + '/'

    # ----------
    # BUILD DATASET STRING
    # This will add all the input filenames to the recipe.
    #DATASET_STRING =  '\n<dataset InputDir="' + INPUT_DIRECTORY + \
    #                  '" OutputDir="' + REDUCED_DIRECTORY + '/' + \
    #                  SHORTDATE + '/' + TARGET_NAME + '-' + WAVEBAND + '/">'
    DATASET_STRING =  '\n<dataset InputDir="' + INPUT_DIRECTORY + \
                      '" OutputDir="' + FILE_OUTPUT_DIRECTORY + '">'
    #                  '" OutputDir="' + REDUCED_DIRECTORY + '/' + \
    #                  TARGET_NAME + '-' + WAVEBAND + '/' + DATE_OBS + '/">'
    for F in range(len(INPUT_FILES)):
        DATASET_STRING = DATASET_STRING+ \
                         '\n   <fits FileName="'+INPUT_FILES[F]+'" />'

    # ----------
    # INSERT OPTIONAL KEYWORDS
    PRIMITIVE_STRINGS = TEMPLATE_SPLIT_STRINGS[1]
    if len(KEYS)>0:
        print('Amended Keys from Template File: ', KEYS)
        for KEY in KEYS:
            # Isolate e.g. 'x0' from 'x0="140"'
            SPLIT_KEY = KEY.split('=')
            # Keep full string up to 'x0':
            # (keep the 'join' statement in case multiple keywords share
            # a string, e.g. 'stellar_flux' and 'stellar_flux_err')
            SPLIT_PRIMITIVE_STRINGS = PRIMITIVE_STRINGS.split(SPLIT_KEY[0])
            SPLIT_PRIMITIVE_STRINGS = [SPLIT_PRIMITIVE_STRINGS[0],
                                       SPLIT_KEY[0].join(SPLIT_PRIMITIVE_STRINGS[1:])]
            # Split by spaces to remove the '="140"' before:
            SPACE_SPLIT = SPLIT_PRIMITIVE_STRINGS[1].split(' ')
            # Make one long string out of everything after '="140"' exclusive:
            AFTER_KEY = ' '.join(SPACE_SPLIT[1:])
            PRIMITIVE_STRINGS = SPLIT_PRIMITIVE_STRINGS[0]+\
                                SPLIT_KEY[0]+'='+SPLIT_KEY[1]+\
                                ' '+AFTER_KEY

    # ----------
    # BUILD FINAL OUTPUT STRING
    # The first few lines of the template, plus our additions,
    # then the final lines of the template.
    INSERTED_STRING = AUTHOR_STRING + DATASET_STRING
    OUTPUT_STRING   = TEMPLATE_SPLIT_STRINGS[0] + INSERTED_STRING +\
                      PRIMITIVE_STRINGS

    # ----------
    # Save the recipe as a .xml in the queue directory.
    # Turn the output string into an array for the savetxt() function.
    OUTPUT = np.array([OUTPUT_STRING])
    # Create a new output filename based on the input files.
    FILES_IN_RECIPE = INPUT_FILES[0][:14]
    if len(INPUT_FILES)>1:
        FILES_IN_RECIPE = FILES_IN_RECIPE + '-' + INPUT_FILES[-1][:14]
    OUTPUT_FILENAME  = QUEUE_DIRECTORY + FILES_IN_RECIPE + RECIPE_TYPE
    # fname is new filename, X is data to be saved, fmt is data format.
    np.savetxt(fname=OUTPUT_FILENAME, X=OUTPUT, fmt='%s')

    print('Created recipe for '+TYPE+'.')
    print('Output directory: ', FILE_OUTPUT_DIRECTORY)
    if RETURN_OUTPUT_DIRECTORY > 0:
        return FILE_OUTPUT_DIRECTORY
    else:
        return
