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
merged_dir = '/home/carterrhea/Desktop/test'
ra = '15:32:53.758'
dec = '+30:20:58.52'
z = 0.3621
Flux = True
model_type = 'single'
#Defaults
evt_file = 'merged_evt.fits'
exposure_map = 'broad_thresh.expmap'
bkg_region = ''
#-------------------------------------------------------------#
os.chdir(merged_dir)
scaling = calc_scale(z)
#---------------------Create Annuli---------------------------#
create_ann(ra,dec)
#---------------------Create Profile--------------------------#
calc_profs(evt_file,exposure_map,bkg_region,scaling)
#---------------------PostProcess-----------------------------#
profiles(scaling,Flux)
