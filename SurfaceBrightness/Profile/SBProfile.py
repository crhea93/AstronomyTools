'''
python program to tie together three main components in the creation of a
surface brightess profile
'''
import os
import numpy as np
from Tools.LSCalc import calc_scale
from Tools.fit_rprofile import profiles
from Tools.radial_prof import calc_profs
from Tools.annuli_create import create_ann
#--------------------------INPUTS-----------------------------#
merged_dir = '/home/%%%/Documents/Data/Merged'
ra = '##:##:##.###'
dec = '+/-##:##:##.###'
z = 1
Flux = True
model_type = 'single'
#Defaults
evt_file = 'merged_evt.fits'
exposure_map = 'broad_thresh.expmap'
bkg_region = 'simple_merged_bkg.reg'
#-------------------------------------------------------------#
os.chdir(merged_dir)
#---------------------Create Annuli---------------------------#
create_ann(ra,dec)
#---------------------Create Profile--------------------------#
calc_profs(evt_file,exposure_map,bkg_region)
#---------------------PostProcess-----------------------------#
scaling = calc_scale(z)
profiles(scaling,Flux)
