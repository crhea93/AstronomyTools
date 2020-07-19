'''
Python Routine to wrap together Temperature Map Fitting given a binning scheme
'''
#-----------------------------------------------IMPORTS--------------------------------------------------#
import os
import sys
from read_input import read_input_file
from binned_spectra import create_spectra
from Fits import Fitting
from Plots import plot_Bins, plot_Ab
from WVT import wvt_main
from Hardness_Ratio import hardness_ratio
#-----------------------------------------------READ IN--------------------------------------------------#
inputs = read_input_file(sys.argv[1])
base = inputs['base_dir']+'/'+inputs['Name']
#----------------------------------------------------WVT---------------------------------------------------#
if inputs['wvt'].lower() == 'true':
    wvt_main(inputs,base,inputs['WVT_data'])
#-----------------------------------------------BIN SPECTRA-------------------------------------------#
num_bins = 0
if inputs['bin_spec'].lower() == 'true':
    num_bins = create_spectra(inputs['base_dir']+'/'+inputs['Name'],inputs['WVT_data'],inputs['ObsIDs'],inputs['source_file'],inputs['output_dir'],inputs['WVT_data'])
if inputs['bin_spec'].lower() == 'false':
    num_bins = inputs['num_bins']
#-----------------------------------------------FIT SPECTRA-------------------------------------------#
if inputs['fit_spec'].lower() == 'true':
    Fitting( inputs['base_dir']+'/'+inputs['Name'],inputs['ObsIDs'],inputs['source_file'],int(num_bins),
        inputs['redshift'],inputs['n_H'],inputs['Temp_Guess'],inputs['Temp_data'],inputs['output_dir'],inputs['multi'])
#-----------------------------------------------PLOT FITS-----------------------------------------------#
if inputs['plot'].lower() == 'true':
   plot_Bins(base+'/'+inputs['WVT_data']+'.txt',base+'/'+inputs['Temp_data'],base,inputs['Name'],inputs['Colormap'],
        inputs['stn_target'],base+'/'+inputs['source_file']+'.img')
   plot_Ab(base+'/'+inputs['WVT_data']+'.txt',base+'/'+inputs['Temp_data'],base,inputs['Name'],inputs['Colormap'],
       inputs['stn_target'],base+'/'+inputs['source_file']+'.img')
   hardness_ratio(base,inputs['Name'],inputs['ObsIDs'][0],inputs['stn_target'],num_bins,inputs['output_dir'],inputs['Colormap'],
        base+'/'+inputs['source_file']+'.fits',base+'/'+inputs['WVT_data']+'.txt')
