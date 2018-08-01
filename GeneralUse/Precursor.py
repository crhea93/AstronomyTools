'''
File to create files necessary for WVT algorithm
We only look in the following energy range:
energy=500:10000

'''
import os
from ciao_contrib.runtool import *

def get_filenames():
    filenames = dict()
    for file in os.listdir(os.getcwd()+'/repro'):
        if file.endswith("_evt2.fits"):
            filenames['evt2'] = file
        if file.endswith("_bpix1.fits"):
            filenames['bpix1'] = file
    for file in os.listdir(os.getcwd()+'/primary'):
        if file.endswith("_asol1.fits"):
            filenames['asol1'] = file
    for file in os.listdir(os.getcwd()+'/secondary'):
        if file.endswith("_msk1.fits"):
            filenames['msk1'] = file
    os.chdir('repro')
    return filenames

def dmcopy_func(filenames,region_name):
    dmcopy.infile = filenames['evt2']+"[sky=region("+region_name+".reg)]"
    dmcopy.outfile = region_name+".fits"
    dmcopy.clobber = True
    dmcopy()
    dmcopy.infile = region_name+".fits[filter energy=500:10000]" #Only look in this energy range
    dmcopy.outfile = region_name+".fits"
    dmcopy.clobber = True
    dmcopy()
    dmcopy.infile = region_name+'.fits'
    dmcopy.outfile = region_name+"_image.fits"
    dmcopy.option = 'image'
    dmcopy.clobber = True
    dmcopy()
    return None

def specextract_func(filenames,region_name):
    specextract.infile = filenames['evt2']+"[sky=region("+region_name+".reg)]"
    specextract.outroot = region_name
    specextract.bkgfile = filenames['evt2']+"[sky=region("+region_name+"_bkg.reg)]"
    specextract.asp = '../primary/'+filenames['asol1']
    specextract.mskfile = '../secondary/'+filenames['msk1']
    specextract.badpixfile = filenames['bpix1']
    specextract.grouptype = 'NUM_CTS'
    specextract.binspec = '1'
    specextract.clobber = True
    specextract()
    return None

def exposure_map_func(region_name):
    fluximage.infile = region_name+".fits"
    fluximage.outroot = "flux"
    fluximage.bands = "broad"
    fluximage.binsize = "1"
    fluximage.units = "area"
    fluximage.clobber = True
    fluximage.cleanup = False
    fluximage()
    return None


def process_data(region_name):
    filenames = get_filenames()
    dmcopy_func(filenames,region_name)
    specextract_func(filenames,region_name)
    exposure_map_func(region_name)
    return None



def main():
    chandra_dir = '%%%'
    region_name = '%%%'
    os.chdir(chandra_dir)
    process_data(region_name)
    return None
main()
