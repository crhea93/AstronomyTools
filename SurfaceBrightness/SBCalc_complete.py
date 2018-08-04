'''
Calculate Surface Brightness from Scratch

INPUTS:
    chandra_dir -- full path to data directory (e.g. '/home/user/Documents/Data')
    evt_file -- name of event file without extension (e.g. 'acisf#####_repro_evt2')
    energy_range -- energy range in electron volts (e.g. 500:2000)
    region -- name of region file of interest without .reg extension (e.g. 'simple') Set None if for entire event
    background -- name of background region file without .reg extension (e.g. 'simple_background')
    confidence -- confidence level (e.g. 0.9)
    exposure -- Boolean determining method to calculate Net Energy Flux. See
        Documentation for more information. (e.g. True)

OUTPUTS:
    .par file containing aprates solutions meaning all counts/rates/flux info (e.g. aprates+region.par)
'''
import os
from astropy.io import fits
from ciao_contrib.runtool import *
from ToolBox import calc_effenergy, calc_flux
#------------------INPUTS------------------------------------------------------#
chandra_dir = '%%%'
evt_file = '%%%'
energy_range = '%%%'
region = '%%%'
background = '%%%'
exposure = True
#------------------------------------------------------------------------------#
os.chdir(chandra_dir)


run = input("Do you wish to run SBCalc_complete.py as a standalone app: ")
if run.lower() == 'yes':
    calc_flux(evt_file,energy_range,region,background,exposure)
else:
    print("Not running SBCalc_complete as a standalone app :)")
