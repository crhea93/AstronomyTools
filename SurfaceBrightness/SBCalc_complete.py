'''
Calculate Surface Brightness from Scratch

INPUTS:
    chandra_dir -- full path to data directory (e.g. '/home/user/Documents/Data')
    evt_file -- name of event file without extension (e.g. 'acisf#####_repro_evt2')
    energy_range -- energy range in electron volts (e.g. 500:2000)
    region -- name of region file of interest without .reg extension (e.g. 'simple')
    background -- name of background region file without .reg extension (e.g. 'simple_background')
    confidence -- confidence level (e.g. 0.9)
    exposure -- Boolean determining method to calculate Net Energy Flux. See
        Documentation for more information. (e.g. True)

OUTPUTS:
    .par file containing aprates solutions meaning all counts/rates/flux info (e.g. aprates+region.par)
'''
import os
from astropy.io import fits
from ciao_contrib.runtool import *


def calc_effenergy(region,energy_range2):
    dmtcalc.infile = region+'.arf'
    dmtcalc.outfile = "arf_weights"+str(region)
    dmtcalc.expression = "mid_energy=(energ_lo+energ_hi)/2.0;weights=(mid_energy*specresp)"
    dmtcalc.clobber =True
    dmtcalc()
    dmstat.infile = "arf_weights"+str(region)+"[mid_energy="+str(energy_range2)+"][cols weights]"
    dmstat.verbose = True
    dmstat()
    weight_sum = float(dmstat.out_sum)
    dmstat.infile = "arf_weights"+str(region)+"[mid_energy="+str(energy_range2)+"][cols specresp]"
    dmstat.verbose = True
    dmstat()
    specresp_sum = float(dmstat.out_sum)
    eff_energy = weight_sum/specresp_sum
    print("Our effective energy is: "+str(eff_energy))
    return eff_energy

def calc_flux(evt_file,energy_range,region,background,confidence = 0.9,exposure = False):
    #Rearrange energy ranges
    energies = [float(x) for x in energy_range.split(':')]
    energy_range2 = str(energies[0]/1000)+':'+str(energies[1]/1000) #for effective energy (eV)
    energy_range3 = str(energies[0]/1000)+'-'+str(energies[1]/1000)  #For average effective exposures (eV)
    #Get counts for region and background
    print("Calculating all data needed to calculate flux")
    dmextract.infile = evt_file+".fits[energy="+energy_range+"][bin sky=region("+region+".reg)]"
    dmextract.outfile = region+'_counts.fits'
    dmextract.opt = 'generic'
    dmextract.bkg = evt_file+".fits[energy="+energy_range+"][bin sky=region("+background+".reg)]"
    dmextract.clobber = True
    dmextract()
    dmstat.infile = region+'_counts.fits[cols counts]'
    dmstat()
    counts = float(dmstat.out_sum)
    dmstat.infile = region+'_counts.fits[cols area]'
    dmstat()
    area = float(dmstat.out_sum)
    dmstat.infile = region+'_counts.fits[cols bg_counts]'
    dmstat()
    bg_counts = float(dmstat.out_sum)
    dmstat.infile = region+'_counts.fits[cols bg_area]'
    dmstat()
    bg_area = float(dmstat.out_sum)

    #Set PSF elements
    alpha = 1 #PSF fraction in source aperature; 1-perfect
    beta = 0 #PSF fraction in background aperature; 0-perfect

    #Exposure Time
    hdu = fits.open(evt_file+'.fits')
    hdr = hdu[0].header
    T_s = hdr['TSTOP']-hdr['TSTART']
    T_b = T_s
    hdu.close()

    #Calculate exposure maps
    effen = calc_effenergy(region,energy_range2)
    #Create Exposure Map for the band of interest
    fluximage.infile = evt_file+".fits"
    fluximage.outroot = "exp"
    fluximage.bands = energy_range2+":"+str(effen)
    fluximage.binsize = "1"
    fluximage.units = "default"
    fluximage.clobber = True
    fluximage.cleanup = True
    fluximage()

    #Calculate average effective exposures
    dmstat.punlearn()
    dmstat.infile = 'exp_'+energy_range3+'_thresh.expmap[sky=region('+region+'.reg)]'
    dmstat.centroid = False
    dmstat()
    E_s = dmstat.out_mean
    dmstat.punlearn()
    dmstat.infile = 'exp_'+energy_range3+'_thresh.expmap[sky=region('+background+'.reg)]'
    dmstat.centroid = False
    dmstat()
    E_b = dmstat.out_mean

    #Calculate average photon energies in source and background aperature
    if exposure == False:
        dmtcalc.punlearn()
        dmtcalc.infile = evt_file+".fits[energy="+energy_range+",sky=region("+region+".reg)]"
        dmtcalc.outfile = region+"_source_energy.fits"
        dmtcalc.expression = 'energy=1.6e-12*energy' #Convert to ergs
        dmtcalc.clobber = True
        dmtcalc()
        dmstat.punlearn()
        dmstat.infile = region+'_source_energy.fits[cols energy]'
        dmstat()
        eng_s = dmstat.out_mean
        dmtcalc.punlearn()
        dmtcalc.infile = evt_file+".fits[energy="+energy_range+",sky=region("+background+".reg)]"
        dmtcalc.outfile = region+"_background_energy.fits"
        dmtcalc.expression = 'energy=1.6e-12*energy' #Convert to ergs
        dmtcalc.clobber = True
        dmtcalc()
        dmstat.punlearn()
        dmstat.infile = region+'_background_energy.fits[cols energy]'
        dmstat()
        eng_b = dmstat.out_mean

        #set flux_s,flux_b to zero to ignore exposure
        flux_s = 1; flux_b = 1


    if exposure == True:
        eff2evt.punlearn()
        eff2evt.infile = evt_file+".fits[energy="+energy_range+"][sky=region("+region+".reg)]"
        eff2evt.outfile = region+"_source_effexp.fits"
        eff2evt.clobber = True
        eff2evt()
        dmstat.punlearn()
        dmstat.infile = region+'_source_effexp.fits[cols flux]'
        dmstat()
        flux_s = dmstat.out_mean
        eff2evt.punlearn()
        eff2evt.infile = evt_file+".fits[energy="+energy_range+"][sky=region("+background+".reg)]"
        eff2evt.outfile = region+"_background_effexp.fits"
        eff2evt.clobber = True
        eff2evt()
        dmstat.punlearn()
        dmstat.infile = region+'_background_effexp.fits[cols flux]'
        dmstat()
        flux_b = dmstat.out_mean

        #Conversely set eng_s,eng_b to one to signify we are using effective exposure
        eng_s = 1; eng_b = 1



    #Calculate energy flux and bounds
    print("Setting aprates values")
    aprates.punlearn()
    aprates.conf = confidence
    aprates.n = counts
    aprates.m = bg_counts
    aprates.A_s = area
    aprates.A_b = bg_area
    aprates.alpha = alpha
    aprates.beta = beta
    aprates.T_s = T_s
    aprates.T_b = T_b
    aprates.E_s = E_s
    aprates.E_b = E_b
    aprates.eng_s = eng_s
    aprates.eng_b = eng_b
    aprates.flux_s = flux_s
    aprates.flux_b = flux_b
    aprates.outfile = 'aprates_'+region+'.par'
    aprates.clobber = True
    aprates.pdf = 'alternate'
    print("Running aprates for flux value")
    aprates()

    return None


def main():
    chandra_dir = '%%%'
    os.chdir(chandra_dir)
    evt_file = '%%%'
    energy_range = '%%%' #in electron volts
    region = '%%%'
    background = '%%%'
    confidence = 0.9
    exposure = True

    calc_flux(evt_file,energy_range,region,background,confidence,exposure)

main()