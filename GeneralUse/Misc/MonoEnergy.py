'''
Quick script to calculate monochromatic energy
'''
import os
from ToolBox import calc_effenergy
#---------------INPUTS------------------------#
chandra_dir = '%%%'
region = '40kpc_merged'
energy_range = '0.5:2.0'
#---------------------------------------------#
os.chdir(chandra_dir)
mono = calc_effenergy(region,energy_range)
print("The effective monochromatic energy is: "+str(mono))
