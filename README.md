**Dependencies**
python3
astropy
matplotlib
numpy
scipy
The IDL Data Reduction Pipeline (DRP). To download the most recent version here: https://github.com/geminiplanetimager/gpi_pipeline

**Installation**
Download git repository. This code acts as a wrapper for the IDL Data Reduction Pipeline (DRP) for Gemini Planet Imager (GPI). http://docs.planetimager.org/pipeline/ . This wrapper was written to work with version 1.5.0 of the DRP to reduce GPI1.0 polarized light data. Older versions of DRP will be incompatable.
Follow the instructions to install the DRP.

**Gemini-LIGHTS pipeline wrapper**
This wrapper functions to modify template recipies. These recipies are normally modified (when only using the IDL DRP code) with the Recipe Editor. This pipeline allows to bypass the GUI and automate the process. The Gemini-LIGHTS pipeline has also implemented sanity checks such as centering, flexure, stellar/insturmental polarization to ensure accurate and reproducible results.

The python wrapper works by taking key parameters indicated in the template file in the Raw directory. 

**Usage**
In this example, I will use the example for the target MWC 275 which was observed in 2014. I will walk you step by step and point out files that are important to look at during each step.

1) Start the DRP with IDL using the gpi-pipeline command in the terminal.
2) Place necessary calibration files into the calibration directory
3) Update the list of calibration files by clicking Rescan Calib. DB in the bottom left corner of the GPI DRP Status Console
4) list names of all files into a files.lst file located in the targets Raw directory. This can easily be generated with the terminal command: ls *.fits > files.lst
5) execute python3 python3 pipeline.py mwc275-J.txt mwc275/20140424-J/ first
6) The first keyword will:
  A) Perform sanity check on the files within files.lst to ensure the files are in the right time order and the waveplate angle WPANGLE is in the proper order.
  B) Create the first PODC file based on the first file in files.lst
  C) end the pipeline
7) Inspect the sanitycheck text file located in the targets raw directory (mwc275-J_sanitycheck.txt). Search for files in the wrong order, and WPANGLE observations that are in the right order and that there are 4 rotations of the half wave plate (0,22.5,45,67.5) to make an individular stokes.
8) Inspect the podc file (suffix _bpfix.fits). Search for bright binaries or failure of the flexure. Adjust these parameters in the template file (mwc275-J.txt) as needed.
9) execute python3 python3 pipeline.py mwc275-J.txt mwc275/20140424-J/
10) No keyword will execute the entire pipeline producing _combined_rstokesdc.fits files. 
11) 
python3 pipeline.py template.txt /location_in_Raw_directory/ optional_key_word
where template.txt is a text file listing the requested reduction parameters 


**Outline of file outputs** 
Diagnostic Images:

Output FITS files:

