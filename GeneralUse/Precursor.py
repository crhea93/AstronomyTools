'''
Run basic specextract and create simple fits images
We only look at the broad energy range (0.5-10keV)

INPUTS:
    chandra_dir -- full path to chandra data including obsid (e.g. '/home/user/Documents/ChandraData/20940')
    repro_dir -- name of reprocessed data directory (e.g. 'repro')
    region_name -- name of region for specextract (e.g. '400kpc')
    background_name -- name of background region (e.g. 'simple')


'''
import os
from ciao_contrib.runtool import *

def get_filenames(repro_dir):
    filenames = dict()
    for file in os.listdir(os.getcwd()+'/'+repro_dir):
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
    os.chdir(repro_dir)
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
    dmcopy.infile = filenames['evt2']
    dmcopy.outfile = "image.fits"
    dmcopy.option = 'image'
    dmcopy.clobber = True
    dmcopy()
    return None

def specextract_func(filenames,region_name,background_name):
    specextract.infile = filenames['evt2']+"[sky=region("+region_name+".reg)]"
    specextract.outroot = region_name
    specextract.bkgfile = filenames['evt2']+"[sky=region("+background_name+"_bkg.reg)]"
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
    fluximage.cleanup = True
    fluximage()
    return None


def process_data(repro_dir,region_name,background_name):
    filenames = get_filenames(repro_dir)
    dmcopy_func(filenames,region_name)
    specextract_func(filenames,region_name,background_name)
    exposure_map_func(region_name)
    return None



def main():
    chandra_dir = '%%%'
    repro_dir = '%%%'
    region_name = '%%%'
    background_name = '%%%'
    os.chdir(chandra_dir)
    process_data(repro_dir,region_name,background_name)
    return None
main()
