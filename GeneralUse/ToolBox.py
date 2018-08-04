'''
Collection of helpful subroutines used in an assortment of calculations
SUBROUTINES: (in no particular order)
    calc_effenergy - calculate effective monochromatic energy
    calc_flux - calculate flux parameters and bounds using aprates
    calc_bounds - calculate energy flux bounds from .par file
    simulatePSF - simulate psf on event file to create .psf file
'''
import os
from ciao_contrib.runtool import *
#-------------------------------------------------#
#-------------------------------------------------#
'''
calculate effective monochromatic energy
    parameter:
        region - region of interest (e.g. 'simple')
        energy_range2 - energy range in kiloelectron volts (e.g. '0.5:2.0')
'''
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
#-------------------------------------------------#
#-------------------------------------------------#
'''
Calculate various quantities considered surface brightness such as:
    - net counts
    - net count Rate
    - net photon flux
    - net energy flux (two options)
        see further documentation
    parameters:
        evt_file - classic event fits file (e.g. 'acsif_#####_repro_evt2')
        energy_range - energy range in electron volts (e.g. 500:2000)
        region - region of interest (e.g. 'simple')
        background - background .reg file without extension (e.g. 'simple_bkg')
        exposure - boolean to use exposure fluxes (e.g. True) (See documentation)
    outputs:
        .par file containing all calculated quantities (.e.g. 'aprates_'+region+'.par')
'''
def calc_flux(evt_file,energy_range,region,background,exposure = False):
    #Rearrange energy ranges
    energies = [float(x) for x in energy_range.split(':')]
    energy_range2 = str(energies[0]/1000)+':'+str(energies[1]/1000) #for effective energy (eV)
    energy_range3 = str(energies[0]/1000)+'-'+str(energies[1]/1000)  #For average effective exposures (eV)
    #Get counts for region and background
    print("Calculating all data needed to calculate flux")
    if region != None:
        dmextract.infile = evt_file+".fits[energy="+energy_range+"][bin sky=region("+region+".reg)]"
    if region == None:
        dmextract.infile = evt_file+".fits[energy="+energy_range+"[bin sky=region("+evt_file+".fits)]"
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
    if region != None:
        dmstat.infile = 'exp_'+energy_range3+'_thresh.expmap[sky=region('+region+'.reg)]'
    if region == None:
        dmstat.infile = 'exp_'+energy_range3+'_thresh.expmap[sky=region('+evt_file+'.fits)]'
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
        if region != None:
            dmtcalc.infile = evt_file+".fits[energy="+energy_range+",sky=region("+region+".reg)]"
        if region == None:
            dmtcalc.infile = evt_file+".fits[energy="+energy_range+",sky=region("+evt_file+".fits)]"
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
        if region != None:
            eff2evt.infile = evt_file+".fits[energy="+energy_range+"][sky=region("+region+".reg)]"
        if region == None:
            eff2evt.infile = evt_file+".fits[energy="+energy_range+"][sky=region("+evt_file+".fits)]"
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
    aprates.conf = 0.90
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
#-------------------------------------------------#
#-------------------------------------------------#
'''
Simulate PSF using marx
    parameters:
        evt_file - associated event file (string)
        psf_file - associated psf fits file (string)
        outname - name of .psf file to be created (string)
    outputs:
        Creates a .psf for event
'''
def simulatePSF(evt_file,psf_file,outname):
    #Get source RA and dec from PSF
    temp_a = os.popen('dmlist '+psf_file+' header | grep SRC_RA').readlines()
    ra = temp_a[0].split("[deg]")[0].split("SRC_RA")[1].strip()
    temp_b = os.popen('dmlist '+psf_file+' header | grep SRC_DEC').readlines()
    dec = temp_b[0].split("[deg]")[0].split("SRC_DEC")[1].strip()
    #Simulate PSF
    simulate_psf.punlearn()
    simulate_psf.infile = evt_file+'.fits'
    simulate_psf.outroot = outname
    simulate_psf.ra = ra
    simulate_psf.dec = dec
    simulate_psf.simulator = 'file'
    simulate_psf.rayfile = psf_file
    simulate_psf.projector = 'marx'
    simulate_psf()
    return None
#-------------------------------------------------#
#-------------------------------------------------#
'''
 Calculate Energy Flux Bounds from .par file created from aprates
   parameters:
       region - region of interest for which the .par file was created (string)
       quantity_to_calc - Net Counts/Rates/Flux option (string)
            Options for quantity to calculate:
                NC - Net Counts
                NCR - Net Count Rates
                NPF - Net Photon Flux
                NEFA - Net Energy Flux option A
                NEFB - Net Energy Flux option B
   outputs:
       val - calculated value of parameter of interest (float)
       lower - lower confidence bound (float)
       upper - upper confidence bound (float)
'''
def calc_bounds(region,quantity_to_calc):
    with open('aprates_'+region+'.par') as f:
        data = []
        count = 0
        for line in f:
            if count < 35:
                data.append(line.split(',')[3])
            count += 1
    if quantity_to_calc == 'NC':
        val = float(data[0])
        lower = float(data[1])
        upper = float(data[2])
        print("Net Counts is calculated at %.2E with an lower bound of %.2E and an upper bound of %.2E"%(val,lower,upper))
    if quantity_to_calc == 'NCR':
        val = float(data[7])
        lower = float(data[8])
        upper = float(data[9])
        print("Net Count Rate is calculated at %.2E with an lower bound of %.2E and an upper bound of %.2E"%(val,lower,upper))
    if quantity_to_calc == 'NPF':
        val = float(data[15])
        lower = float(data[16])
        upper = float(data[17])
        print("Net Photon Flux is calculated at %.2E with an lower bound of %.2E and an upper bound of %.2E"%(val,lower,upper))
    if quantity_to_calc == 'NEFA' or quantity_to_calc == 'NEFB':
        val = float(data[28])
        src_val = float(data[7])
        src_rates_lower = float(data[8])
        src_rates_upper = float(data[9])
        lower = val*(src_rates_lower/src_val)
        upper = val*(src_rates_upper/src_val)
        print("Net Energy Flux is calculated at %.2E with an lower bound of %.2E and an upper bound of %.2E"%(val,lower,upper))

    return val,lower,upper
#-------------------------------------------------#
#-------------------------------------------------#
