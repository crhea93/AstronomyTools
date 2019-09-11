'''
Create Radial Surface Brightness Profile

Be sure to have already merged using Merge.py, created annuli using annuli_create.py,
and chosen a background region using ds9

OUTPUTS:
    rprofile_nonincl_r_data.fits - Fits file containing SB information of independent annuli
    rprofile_r_data.fits- Fits file containing SB information of compounded annuli
'''
import os
from astropy.io import fits
from ciao_contrib.runtool import *

#---------------INPUTS------------------#
chandra_dir = '/home/carterrhea/Documents/Data/Merged/Bin1/Soft'
evt_file = 'merged_evt.fits'
exposure_map = 'broad_thresh.expmap'
bkg_region = 'simple_merged_bkg.reg' #Must have this created!
#---------------------------------------#

def calc_profs(evt_file,exposure_map,bkg_region,scaling):
    dmextract.punlearn()
    dmextract.infile = evt_file+"[bin sky=@annuli.reg]"
    dmextract.outfile = "rprofile.fits"
    #dmextract.bkg = evt_file+"[bin sky=region("+bkg_region+")]"
    dmextract.exp = exposure_map
    #dmextract.bkgexp = exposure_map
    dmextract.opt = 'generic'
    dmextract.clobber = True
    dmextract()

    dmtcalc.punlearn()
    dmtcalc.infile = "rprofile.fits"
    dmtcalc.outfile = "rprofile_rmid.fits"
    dmtcalc.expression = "rmid=0.5*(R[0]+R[1])"
    dmtcalc.clobber = True
    dmtcalc()


    dmcopy.punlearn()
    dmcopy.infile = "rprofile_rmid.fits[cols rmid,sur_bri,sur_bri_err]"
    dmcopy.outfile = "rprofile_rmid_data.fits"
    dmcopy.clobber = True
    dmcopy()

    d_ = fits.open("rprofile_rmid_data.fits")
    for i in range(len(d_[1].data[:])):
        d_[1].data[:][i][0] = d_[1].data[:][i][0]*scaling
    d_.writeto('rprofile_rmid_data.fits',overwrite=True)    

    return None
