'''
Small script to effective energy
'''
import os
import sys
sys.path.append('../GeneralUse')
from ToolBox import calc_effenergy
from ciao_contrib.runtool import *

#----------------------INPUTS--------------------------------------------------#
chandra_dir = '%%%'
os.chdir(chandra_dir)
ra_central = '%%%'
dec_central = '%%%'
region_names = ['%%%']
background_name = '%%%'
energy_range = '%%%' #(e.g. 0.5:2.0)
#------------------------------------------------------------------------------#

def flux_calc(event_file,region_name,ra_central,dec_central,background_name,energy_range,eff_energy):
    srcflux.infile = event_file
    srcflux.pos = ra_central+" "+dec_central
    srcflux.outroot = region_name+"/"
    srcflux.bands = energy_range+":"+str(eff_energy)
    srcflux.srcreg = region_name+".reg"
    srcflux.bkgreg = background_name+".reg"
    srcflux.clobber = True
    print(srcflux())

    return None

def calc_sb(region_names,background_name,energy_range,ra_central,dec_central):
    event_file = None
    for file in os.listdir(os.getcwd()):
        if file.endswith("_evt2.fits"):
            event_file = file
    for region_name in region_names:
        effen = calc_effenergy(region_name,energy_range)
        flux_calc(event_file,region_name,ra_central,dec_central,background_name,energy_range,effen)

    return None

def main():
    calc_sb(region_names,background_name,energy_range,ra_central,dec_central)
    return None
main()
