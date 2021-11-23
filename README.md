# Gemini_LIGHTS_pipeline
Pipeline to reduce GPI polarized light data. Program uses the IDL DRP created for GPI

  This program will help reduce high-contrast polarized light imagery from the Gemini Planet Imager. 
This code acts as a wrapper around the IDL Data Reduction Pipline (DRP) written to reduce GPI data,
thus the IDL DRP pipeline must be installed for this Gemini-LIGHTS pipeline to opportate. You can 
find a description of the DRP pipeline here: http://docs.planetimager.org/pipeline/
This program automates the IDL DRP pipeline, makes adjustments to it's reduction, and allows for
greater reproducablility.

Dependencies:
numpy, astropy, matplotlib
Requires the instillation of the IDL DRP pipeline. Written for verison 1.5 which the latest version
can be found here: (https://github.com/geminiplanetimager/gpi_pipeline)

This program was written by Anna Laws and Evan Rich and utilized in Laws et al. 2020 and 
Rich et al. In Prep. Please cite both works when utilizing this program along with the relevent
IDL DRP papers described here: http://docs.planetimager.org/pipeline/
ADS: https://ui.adsabs.harvard.edu/abs/2020ApJ...888....7L/abstract


To execute the program, first start IDL DRP with the gpi-pipeline command. To run, execute:
python3 pipeline.py mwc275-J.txt mwc275/20140424-J/
where pipeline.py is the main program, 'mwc275-J.txt' is the parameter file for the reduction, and 
'mwc275/20140424-J/' is the sub-directory within the Raw directory where the raw data files for 
mwc275 and the parameter file are located. All output files are saved within the Reduced directory
with the sub-directory labed based on the target name, flux band observed, and observational date.

This readme file is incomplete. Further discripton is needed of the flux calibration, the parameter file structure, and output files.

