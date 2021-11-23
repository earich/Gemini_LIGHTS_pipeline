#!/usr/bin/env python
"""
Sanity check that the directory contains suitable files for the
pipeline.

Command line usage: sanity.py <directory> <optional args>

Returns a table of the raw data files in the directory,
their observation times and half-wave plate angles,
and any optional args too.
Sanity check is passed if time difference is regular
and angles have regular intervals (e.g. 45deg difference between files).

Recommended optional args: 'AIRMASS', 'RAWIQ', 'RMSERR', 'DIMMSEE'.
"""
# %%%%%%%%%%% IMPORTS %%%%%%%%%%%
import sys
import os
from gpi_analysis.inputs import getfitskeywords, getfilenames, getfilenames_l
from datetime import datetime, timedelta
import numpy as np

# %%%%%%%%%%% FUNCTIONS %%%%%%%%%%%
# def getfilenames(DIRECTORY, END='.fits'):
#     """
#     Originally in plotfourtest - find permanent home
#     """
#     NAMES     = os.listdir(DIRECTORY)
#     LEN       = len(END)
#     FITSNAMES = [ DIRECTORY+A for A in NAMES if A[-LEN:]==END ]
#     return FITSNAMES

def maketable(NROWS, NCOLS, HEADINGS):
    """
    Make a blank table for filling with info.
    """
    TABLE      = np.full((NROWS+1,NCOLS),-100,dtype=object)
    TABLE[0,:] = HEADINGS
    return TABLE

def find_obs_midpoint(TABLE, WHICH):
    """
    Find midpoint time of observations.
    """
    # Find UT (or other time) column in table:
    COL      = np.where(TABLE[0]==WHICH)
    COL      = COL[0][0]
    # Make a copy of this column. (Else values may be overwritten)
    VALS     = np.copy(TABLE[:,COL])
    FMT      = '%H:%M:%S'
    START    = datetime.strptime(VALS[1][:8],  FMT)
    END      = datetime.strptime(VALS[-1][:8], FMT)
    # Midpoint in datetime.datetime format, M:
    M        = START + (END - START)/2
    # Convert M to '%H:%M:%S' string, MIDPOINT:
    MIDPOINT = M.strftime('%X %x %Z')[:8]
    print('Midpoint of observations: ', MIDPOINT)
    return MIDPOINT

def savetable(RAWDATA_DIRECTORY, TABLE, KEYWORDS,
              STATUS_UT, STATUS_WPANGLE,
              TARGET,    BAND, MIDPOINT):
    """
    Save a table of file names, obs times, hwp angles, and any
    optional keywords.

    Header contains target, band, midpoint time, and sanity check
    statuses.
    """
    # Save txt - include status as comment lines.
    FNAME  = RAWDATA_DIRECTORY + TARGET+'-'+BAND[0] + '_sanitycheck.txt'
    HEADER = TARGET+'-'+BAND[0]                  + '\n' +\
             'Observation midpoint: ' + MIDPOINT + '\n' +\
             STATUS_UT                           + '\n' +\
             STATUS_WPANGLE
    FMT    = "%30s %12s %8s %8s %13s" #'%s'#"%60s %12s %8.1e %8.1e %8.1e"
    if len(KEYWORDS)>3:
        for K in range(3,len(KEYWORDS)):
            FMT+=' %13s'
    np.savetxt(fname=FNAME,X=TABLE,header=HEADER,fmt=FMT)
    return

def sanitycheck(TABLE, WHICH, LIMIT):
    """
    Sanity check - do values in column have regular intervals?
    Returns the table updated with the intervals and a pass/fail string.

    Recommend limit = 80 seconds different.
    """
    # Isolate the relevant column in the table.
    COL        = np.where(TABLE[0]==WHICH)
    if len(COL)<1:
        STATUS = 'No sanity check performed for '+WHICH+' - not in table.'
        print(STATUS)
        return TABLE, STATUS
    # Make a copy of this column. (Else values may be overwritten)
    COL        = COL[0][0]
    VALS       = np.copy(TABLE[:,COL])
    #
    # Find difference between successive values in column.
    # Make a new column to store these diffs:
    DIFFS      = np.empty(len(VALS),dtype=object)
    DIFFS[0]   = WHICH+'_DIFFS'
    #
    X          = 1
    INSANE     = 0
    while X<len(VALS):
        if X+4>=len(VALS):
            CYCLE = VALS[X:]
        # Only need to sanity check in blocks of four
        # (four files is one cycle --> one stokesdc).
        else:
            CYCLE = VALS[X:X+4]
        if isinstance(VALS[1],str):
            # If in datetime format, change to seconds for calculation.
            CYCLE_DIFFS=[]
            for Y in range(1,len(CYCLE)):
                FMT = '%H:%M:%S'
                TIME_DIFF = datetime.strptime(CYCLE[Y][:8],  FMT) -\
                            datetime.strptime(CYCLE[Y-1][:8],FMT)
                TIME_DIFF = timedelta.total_seconds(TIME_DIFF)
                CYCLE_DIFFS.append(TIME_DIFF)
            CYCLE_DIFFS   = np.array(CYCLE_DIFFS)
        else:
            # Set any HWP angles of 360.0 to 0.0.
            ZEROS = np.where(CYCLE > 359.0)
            CYCLE[ZEROS] = CYCLE[ZEROS]-360.0
            CYCLE_DIFFS  = np.diff(CYCLE)
        # Perform sanity check.
        # Subtract mean value from all differences.
        # If any remaining value is above the limit, flag as insane.
        # For sane HWP angles, difference will be constant -> all values 0.
        if len(CYCLE)>1:
            CHECK_DIFFS  = np.abs(  CYCLE_DIFFS - np.median(CYCLE_DIFFS))
            BEYOND_LIMIT = np.where(CHECK_DIFFS > LIMIT)
            if len(BEYOND_LIMIT[0])>0:
                INSANE  += 1
        else:
            print('Cycle length less than 2.')
            INSANE += 1
        # Store results.
        # Add a zero for the first difference (VALS[X]-VALS[X]=0)
        # for the sake of the output table.
        CYCLE_DIFFS  = np.append([0.0],CYCLE_DIFFS)
        DIFFS[X:X+4] = CYCLE_DIFFS
        X=X+4

    # Add the DIFFS column to the existing table
    # next to its relevant column.
    DIFFS      = np.transpose(DIFFS)
    TABLE      = np.concatenate( ( TABLE[:,:(COL+1)],
                                  DIFFS[:,None],
                                  TABLE[:,(COL+1):]
                                  ), axis=1 )
    # Make a copy of the pass/fail status to print and save to file.
    if INSANE>0:
        STATUS = 'Failed sanity check: '+WHICH+\
                 '. Limit exceeds '+str(LIMIT)+\
                 ' in '+str(INSANE)+' cycle(s).'
    else:
        STATUS = 'Passed sanity check: '+WHICH+' with limit '+str(LIMIT)+'.'
    print(STATUS)
    return TABLE, STATUS


# %%%%%%%%%%% MAIN FUNCTION %%%%%%%%%%%
def main():
    rawdata_directory    = sys.argv[1]
    keywords             = ['FILE', 'UT', 'WPANGLE']
    if len(sys.argv)>2:
        keywords        += sys.argv[2:]

    #rawdata_names        = getfilenames(rawdata_directory)
    rawdata_names        = getfilenames_l(rawdata_directory,'files.lst')
    target               = getfitskeywords(rawdata_names[0], 'OBJECT')
    band                 = getfitskeywords(rawdata_names[0], 'OBSMODE')
    print(target+'-'+band[0])

    table                = maketable( len(rawdata_names),
                                      len(keywords),
                                      HEADINGS=keywords )
    # Populate table:
    for x, rawdata_file in enumerate(rawdata_names,start=1):
        table[x,0]        = rawdata_file.split('/')[-1]
        for y, keyword in enumerate(keywords[1:], start=1):
            table[x,y]    = getfitskeywords(rawdata_file, keyword)

    table, status_ut      = sanitycheck(table, 'UT',      LIMIT=80)
    table, status_wpangle = sanitycheck(table, 'WPANGLE', LIMIT=5.0)
    midpoint              = find_obs_midpoint(table, 'UT')

    savetable(rawdata_directory, table, keywords,
              status_ut, status_wpangle,
              target,    band,   midpoint)


# %%%%%%%%%%% MAIN FUNCTION %%%%%%%%%%%
if __name__ == '__main__':
    main()
