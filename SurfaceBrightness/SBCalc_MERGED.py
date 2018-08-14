'''
Calculate Surface Brightness from Scratch for MERGED Images

This involves created merged folders for each region and energy range

PLEASE RUN SPECEXTRACT ON EACH OBSERVATION FOR EACH REGION BEFORE RUNNING

INPUTS:
    chandra_dir -- full path to data directory (e.g. '/home/user/Documents/Data')
    output_dir -- name of output directory (e.g. 'Merged_Soft')
    obs_to_merge -- list of OBSIDs for which we want to calculate the SB (e.g. ['#####'])
    repro_dir -- name of reprocessed data directory (e.g. 'repro')
    evt_file -- name of event file without extension (e.g. 'acisf#####_repro_evt2')
            Also used to calculate on the fly exposure map
    energy_range -- energy range in electron volts (e.g. '500:2000')
    region -- name of region file of interest without .reg extension (e.g. 'simple')
    background -- name of background region file without .reg extension (e.g. 'simple_background')
    exposure -- Boolean determining method to calculate Net Energy Flux. See
        Documentation for more information. (e.g. True)


OUTPUTS:
    .par file containing aprates solutions meaning all counts/rates/flux info (e.g. aprates+region.par)
'''
import os
from shutil import copyfile
from astropy.io import fits
from ciao_contrib.runtool import *
#------------------INPUTS------------------------------------------------------#
chandra_dir = '%%%'
output_dir = '%%%'
obs_to_merge = ['#####','#####']
repro_dir = '%%%'
evt_file = '%%%'
energy_range = '####:####' #in electron volts
regions = ['%%%','%%%'] #set None if for entire image
background = '%%%'
exposure = False
#------------------------------------------------------------------------------#
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
            if merged ('merged_evt')
        energy_range - energy range in electron volts (e.g. 500:2000)
        region - region of interest (e.g. 'simple')
        background - background .reg file without extension (e.g. 'simple_bkg')
        exposure - boolean to use exposure fluxes (e.g. True) (See documentation)
        merged - boolean for merged data set or not (e.g. True)
    outputs:
        .par file containing all calculated quantities (.e.g. 'aprates_'+region+'.par')
    Notes:
    Usually we use the region name along with the arf files to calculate the monochromatic
        energy, but if the data set is merged then we must use the evt_file name (see documentation).
        This is handled in the code but be sure to name things appropriately!
'''
def calc_flux(evt_file,energy_range,region,background,exposure = False,merged_obs = ['']):
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
    T_s = 0
    T_b = 0
    for obsid in  merged_obs:
        hdu = fits.open(obsid+'.fits')
        hdr = hdu[0].header
        T_s += hdr['TSTOP']-hdr['TSTART']
        T_b += T_s
        hdu.close()

    #Calculate average effective exposures
        dmstat.punlearn()
        dmstat.infile = energy_range3+'_thresh.expmap[sky=region('+region+'.reg)]'
        dmstat.centroid = False
        dmstat()
        E_s = dmstat.out_mean
        dmstat.punlearn()
        dmstat.infile = energy_range3+'_thresh.expmap[sky=region('+background+'.reg)]'
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

def create_arf(obs_to_merge,region,repro_dir):
    #Create arf files
    arf_files = ''
    pi_files = ''
    for obsid in obs_to_merge:
        arf_files += obsid+'/'+repro_dir+'/'+region+'.arf,'
        pi_files += obsid+'/'+repro_dir+'/'+region+'.pi,'
    arf_files = arf_files[:-1]#get rid of final comma
    pi_files = pi_files[:-1]
    addresp.punlearn()
    addresp.infile = ''
    addresp.arffile = arf_files
    addresp.phafile = pi_files
    addresp.outfile = ''
    addresp.outarf = region+'_merged.arf'
    addresp.clobber = True
    addresp()

def merge_observations(obs_to_merge,output_dir,repro_dir,energy_range2,mono_energy):
    #Merge individual region files
    merging_files = ''
    for obsid in obs_to_merge:
        merging_files += obsid+'/'+repro_dir+'/acisf'+obsid+'_repro_evt2.fits,'
    merging_files = merging_files[:-1]
    merge_obs.punlearn()
    merge_obs.infile = merging_files
    merge_obs.outroot = output_dir+'/'
    merge_obs.bands = energy_range2+":"+str(mono_energy)
    merge_obs.clobber = True
    merge_obs()

def main():
    os.chdir(chandra_dir)
    arfs = input('Do we need to create merged ARF files: ')
    if arfs.lower() == 'yes' or arfs == '':
        print("Combining ARF files")
        for region in regions:
            create_arf(obs_to_merge,region,repro_dir)
    if arfs.lower() != 'yes' and arfs != '':
        print("Combined ARFs not being created")
    energies = [float(x) for x in energy_range.split(':')]
    energy_range2 = str(energies[0]/1000)+':'+str(energies[1]/1000)
    mono_energy = calc_effenergy(region+'_merged',energy_range2)
    print("")
    print("We must now created a merged observation file for this energy band...")
    merge_observations(obs_to_merge,output_dir,repro_dir,energy_range2,mono_energy)
    #We need to copy the region files over AND each individual event file
    for region in regions:
        copyfile(chandra_dir+'/'+obs_to_merge[0]+'/repro/'+region+'.reg',chandra_dir+'/'+output_dir+'/'+region+'.reg')
    copyfile(chandra_dir+'/'+obs_to_merge[0]+'/repro/'+background+'.reg',chandra_dir+'/'+output_dir+'/'+background+'.reg')
    for obser in obs_to_merge:
        copyfile(chandra_dir+'/'+obser+'/repro/acisf'+obser+'_repro_evt2.fits',chandra_dir+'/'+output_dir+'/'+obser+'.fits')
    os.chdir(chandra_dir+'/'+output_dir)
    for region in regions:
        print("Calculating flux for "+region)
        calc_flux(evt_file,energy_range,region,background,exposure,obs_to_merge)
main()
