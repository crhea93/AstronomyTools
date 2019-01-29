'''
This python script contains a full reprocessing suite to clean especially faint and extended X-ray observations

INPUTS:
chandra_dir -- path to chandra directory (i.e. '/home/usr/Documents/Data')
OBSID -- OBSID of interest (i.e. '#####')
source -- region file containing source -- only used for astrometric corrections (i.e. 'source.reg')
source_ra -- right ascension of source (i.e. '##:##:##.#')
source_dec -- declination of source (i.e. '##:##:##.#')
output_dir -- name of output directory (i.e. 'repro')
flare_gti -- gti file created from background flare cleaning (i.e. 'ccd2_bkg_clean.gti')

NOTES:
If running the Flares module assumes that we have already created a background light Curve
    and create a gti file from it in the primary directory

'''

import os
from shutil import copyfile
from ciao_contrib.runtool import *

from Astrometric import Astrometric
from Destreak import Destreak
from BadPixel import BadPixel
from Flares import Flares
from Process import Process
#----------------INPUTS-------------------#
chandra_dir = '/home/carterrhea/Documents/Data'
OBSIDS = ['20940']#,'20941','21129']
source = 'extended.reg'
source_ra = '10:49:22.767'
source_dec = '56:40:26.49'
output_dir = 'repro'
flare_gti = 'ccd2_bkg_clean.gti'
#-----------------------------------------#
def get_filenames():
    filenames = dict()
    num_of_biases = 0
    for file in os.listdir(os.getcwd()+'/primary'):
        if file.endswith("_evt2.fits"):
            filenames['evt2'] = os.getcwd()+'/primary/'+ file
        if file.endswith("_asol1.fits"):
            filenames['asol1'] = os.getcwd()+'/primary/'+ file
    for file in os.listdir(os.getcwd()+'/secondary'):
        if file.endswith("_msk1.fits"):
            filenames['msk1'] = os.getcwd()+'/secondary/'+ file
        if file.endswith("_stat1.fits"):
            filenames['stat1'] = os.getcwd()+'/secondary/'+ file
        if file.endswith("_pbk0.fits"):
            filenames['pbk0'] = os.getcwd()+'/secondary/'+ file
        if file.endswith("_evt1.fits"):
            filenames['evt1'] = os.getcwd()+'/secondary/'+ file
        if file.endswith("_mtl1.fits"):
            filenames['mtl1'] = os.getcwd()+'/secondary/'+ file
        if file.endswith("_flt1.fits"):
            filenames['flt1'] = os.getcwd()+'/secondary/'+ file
        if file.endswith("_bias0.fits"):
            CCD_number = file.split("_")[1]
            filenames[CCD_number+'_bias0'] = os.getcwd()+'/secondary/'+ file
            num_of_biases += 1
    return filenames,num_of_biases



def main():
	for OBSID in OBSIDS:
		base_dir = chandra_dir+'/'+OBSID
		os.chdir(base_dir)
		if not os.path.exists(os.getcwd()+'/'+output_dir):
		    os.makedirs(os.getcwd()+'/'+output_dir)
		if os.path.isfile(os.getcwd()+'/'+source):
		    copyfile(os.getcwd()+'/'+source,os.getcwd()+'/'+output_dir+'/'+source)
		if not os.path.isfile(os.getcwd()+'/'+source):
		    print("Please create source.reg file")
		    print("Exiting program")
		    exit()
		filenames,num_of_biases = get_filenames()
		os.chdir(base_dir+'/'+output_dir)

		print("Appling Astrometric Corrections...")
		Astrometric(OBSID,filenames,source,source_ra,source_dec)
		print("Appling Background Flare Information...")
		os.chdir(base_dir)
		Flares(flare_gti,base_dir,output_dir,filenames)
		os.chdir(base_dir+'/'+output_dir)
		print("Destreaking Event File...")
		Destreak(base_dir,output_dir,filenames)
		print("Creating New Badpixel File...")
		BadPixel(base_dir,output_dir,OBSID,filenames,num_of_biases)
		print("Apply GTI and Completing Reprocessing...")
		Process(filenames,OBSID)
		print()
	return None


main()
