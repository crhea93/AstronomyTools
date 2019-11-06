'''
Create a spectrum with a corrected ARF for a diffuse low-count object
'''
import os
import shutil
from ciao_contrib.runtool import *
#-------------------------------INPUTS-----------------------------------------#
energy_range = "500:7000"
#------------------------------------------------------------------------------#
def spec_basic(evt_file,src_reg,obsid):
    #Create basic spectrum with normal arf file (un-corrected)
    specextract.punlearn()
    specextract.infile = evt_file+'[sky=region('+src_reg+'.reg)]'#[energy='+energy_range+']'
    specextract.outroot = src_reg
    specextract.bkgfile = obsid+'_blank.evt[sky=region('+src_reg+'.reg)]'#[energy='+energy_range+']'
    specextract.weight = True
    specextract.weight_rmf = True
    specextract.energy_wmap = '500:10000'#energy_range
    specextract.bkgresp = False
    specextract.clobber = True
    specextract.verbose = '5'
    specextract.grouptype = 'NONE'
    specextract.binspec = 'NONE'
    specextract()
    return None

def fits_and_img(evt_file,src_reg):
    dmcopy.punlearn()
    dmcopy.infile = evt_file+'[sky=region('+src_reg+'.reg)]'
    dmcopy.outfile = src_reg+'.fits'
    dmcopy.clobber = True
    dmcopy()
    #image
    dmcopy.punlearn()
    dmcopy.infile = src_reg+'.fits'
    dmcopy.outfile = src_reg+'.img'
    dmcopy.option = 'image'
    dmcopy.clobber = True
    dmcopy()


def main_extract(chandra_dir,source_dir,OBSIDS,source_reg):
    '''


    '''
    reproccesed_dir = 'repro'
    for obsid in OBSIDS:
        print("Calculating Spectrum for %s"%obsid)
        os.chdir(chandra_dir+'/'+obsid+'/'+reproccesed_dir)
        #Make sure region file is in reprocessed directory
        shutil.copy(source_dir+'/'+source_reg+'.reg',os.getcwd()+'/'+source_reg+'.reg')
        #Get event file
        if len(obsid) == 4:
            evt_file = "acisf0"+obsid+"_repro_evt2.fits"
        else:
            evt_file = "acisf"+obsid+"_repro_evt2.fits"
        #Check that blank sky file is there. If not, make on
        if os.path.exists(os.getcwd()+'/'+obsid+'_blank.evt'):
            print(' Blanksky File Exists')
        else:
                print(' Blanksky file does not exist... making one now...')
                blanksky.punlearn()
                blanksky.evtfile = evt_file
                blanksky.outfile = obsid+'_blank.evt'
                blanksky()
        print(' Extracting Spectra')
        spec_basic(evt_file,source_reg,obsid)
        #Make fits and image files for region
        fits_and_img(evt_file,source_reg)
    return None
