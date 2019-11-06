'''
Python Routine to wrap together Temperature Map Fitting given a binning scheme
What should you already have?
    - binning scheme (WVT perhaps)
    - blank sky background
    - Fits and Image file used for spectral extraction
'''
#-----------------------------------------------IMPORTS--------------------------------------------------#
import os
import sys
from read_input import read_input_file
from binned_spectra import create_spectra
from Fits import PrimeFitting
from Plots import plot_data
from diffuse_specextract_blank import main_extract
from deproject_mod import deproj_final
#-----------------------------------------------READ IN--------------------------------------------------#
inputs = read_input_file(sys.argv[1])
base = inputs['base_dir']+inputs['Name']
#----------------------------------------------------SPECTRA---------------------------------------------------#
if inputs['extract_spectrum'].lower() == 'true':
    for reg_file in inputs['reg_files']:
        main_extract(base,base+'/regions',inputs['ObsIDs'],reg_file)
num_bins = len(inputs['reg_files'])
#----------------------------------------------DEPROJECTION--------------------------------------------#
if num_bins > 1:
    for obsid_ in inputs['dir_list']:
        prefix = inputs['base_dir']+'/'+obsid_+'/repro/el'
        deproj_final(prefix,'.pi',1,num_bins,0.0,prefix,'.deproj')
#-----------------------------------------------FIT SPECTRA-------------------------------------------#
if inputs['fit_spec'].lower() == 'true':
    PrimeFitting( inputs['base_dir']+'/'+inputs['Name'],inputs['ObsIDs'],inputs['source_file'],int(num_bins),inputs['redshift'],inputs['n_H'],inputs['Temp_Guess'],inputs['Temp_data'],inputs['multi'])
#-----------------------------------------------PLOT FITS-----------------------------------------------#
if inputs['plot'].lower() == 'true':
    plot_data(base+'/'+inputs['Temp_data'],base+'/regions/',inputs['reg_files'],base)
    if num_bins > 1:
        plot_data(base+'/'+inputs['Temp_data']+'_deproj',base+'/regions/',inputs['reg_files'],base)
