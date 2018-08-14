'''
Calculate PSF Fraction in region
Inputs:
    chandra_dir - Full path to chandra directory (e.g. '/home/user/Documents/Chandra')
    evt_file - name of event file (e.g. 'acsif#####_repro_evt2')
    region_file  - name of .reg file (e.g. 'simple')
    psf_file - name of .psf file (e.g. 'simple')
    eff_en - Effective monochromatic  energy (e.g. 1.53)
'''
import os
from ciao_contrib.runtool import *
from SBCalc_complete import calc_effenergy,calc_flux

#------------INPUTS------------------------------------------------------------#
evt_file = '%%%'
region_file = '%%%'
psf_file = '%%%'
eff_en = '%%%'



def calculate_PSFfraction(evt_file,region_file,psf_file,eff_en):
    src_psffrac.punlearn()
    src_psffrac.infile = evt_file+'.fits'
    src_psffrac.region = region_file+'.reg'
    src_psffrac.outfile = 'out_temp.fits'
    src_psffrac.energy = eff_en
    src_psffrac.psffile = psf_file+'.psf'
    src_psffrac()

    temp_b = os.popen('dmlist out_temp.fits header | grep PSFFRAC').readlines()
    psfFraction = temp_b[0]
    print("The PSF fraction is "+str(psfFraction))
    return psfFraction
