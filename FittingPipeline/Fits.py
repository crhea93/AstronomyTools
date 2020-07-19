'''
------------------------------------------------------
GOAL:
    Step through bins (spectra) and calculate the temperature value
    of each bin using XSPEC
------------------------------------------------------
INPUTS:
    dir - Full Path to Main Directory (e.g. '/home/user/Documents/Chandra/12833/repro/binned/')
    file_name - FIlename of PI/PHA spectrum (e.g. 'imageA')
    output_file - Filename for output containing temperature information (e.g. 'Temp_bin')
    num_files - number of bins (e.g. 100)
    redshift - redshift of object (e.g. 0.659)
    n_H - Hydrogen Equivalent Column Density in units of 10^{22} atoms cm^{-2} (e.g. 3.6e-2)
    Temp_guess - Guess for Temperature value in units of KeV (e.g. 5)
------------------------------------------------------
OUTPUTS:
    A file which contains the bin number and associated
    temperature and reduced chi-squared
------------------------------------------------------
ADDITIONAL NOTES:
    Be sure to run heainit first
------------------------------------------------------
'''
#from astropy.io import fits
import os
import subprocess
from sherpa.optmethods import LevMar
from sherpa.stats import LeastSq
from sherpa.plot import DataPlot
from sherpa.astro.xspec import *
from sherpa.astro.all import *
from sherpa.astro.ui import *
from pychips.all import *
from sherpa.all import *
from sherpa_contrib.xspec.xsconvolve import load_xscflux
load_xscflux("cflux")
#TURN OFF ON-SCREEN OUTPUT FROM SHERPA
import logging
logger = logging.getLogger("sherpa")
logger.setLevel(logging.WARN)
logger.setLevel(logging.ERROR)
def set_log_sherpa():
    p = get_data_plot_prefs()
    p["xlog"] = True
    p["ylog"] = True
    return None

def isFloat(string):
    if string == None:
        return False
    try:
        float(string)
        return True
    except ValueError:
        return False
    except TypeError:
        return False

#------------------------------INPUTS------------------------------------------#

#------------------------------------------------------------------------------#
#Dynamically set source for OBSID
def obsid_set(src_model_dict,bkg_model_dict,obsid,bkg_src, obs_count,redshift,nH_val,Temp_guess):
    '''
    Function to set the source and background model for the observation
    :Params
    :src_model_dict : Dictionary of all the source models set for particular region
    :bkg_model_dict : Dido except for the background models
    :obsid : The source pha/pi file
    :bkg_src : The background pha/pi file
    :obs_count : The count of the current observation

    '''
    load_pha(obs_count,obsid,use_errors=True) #Read in
    if obs_count == 1:
        src_model_dict[obsid] = xsphabs('abs'+str(obs_count)) * xsapec('apec'+str(obs_count)) #set model and name
        # Change src model component values
        get_model_component('apec' + str(obs_count)).kT = Temp_guess
        get_model_component('apec' + str(obs_count)).redshift = redshift  # need to tie all together
        get_model_component('apec' + str(obs_count)).Abundanc = 0.3
        thaw(get_model_component('apec' + str(obs_count)).Abundanc)
        get_model_component('abs1').nH = nH_val  # change preset value
        freeze(get_model_component('abs1'))
    else:
        src_model_dict[obsid] = get_model_component('abs1') * xsapec('apec' + str(obs_count))
        get_model_component('apec'+str(obs_count)).kT = get_model_component('apec1').kT #link to first kT
        get_model_component('apec' + str(obs_count)).redshift = redshift
        get_model_component('apec' + str(obs_count)).Abundanc = get_model_component('apec1').Abundanc  # link to first kT
    set_source(obs_count, src_model_dict[obsid]) #set model to source
    set_bkg(obs_count, unpack_pha(bkg_src))
    bkg_model_dict[obsid] = xsapec('bkgApec'+str(obs_count))+get_model_component('abs1')*xsbremss('brem'+str(obs_count))
    set_bkg_model(obs_count,bkg_model_dict[obsid])
    #Change bkg model component values
    et_model_component('bkgApec' + str(obs_count)).kT = 0.18
    freeze(get_model_component('bkgApec'+str(obs_count)).kT)
    get_model_component('brem' + str(obs_count)).kT = 40.0
    freeze(get_model_component('brem' + str(obs_count)).kT)

    return None
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
def FitXSPEC(spectrum_files,background_files,redshift,n_H,Temp_guess,grouping,spec_count,plot_dir):
    """
    Function to fit spectra using sherpa and XSPEC
    Args:
        spectrum_files: List of spectrum files for each ObsID
        background_files: List of background files for each ObsID
        redshift: Redshift of cluster
        n_H: Column density
        temp_guess: Initial temperature guess
        grouping: Number of counts to bin in sherpa fit
        spec_count: Bin number
        plot_dir: Path to plot directory
    Return:
        Spectral fit parameters and their errors
    """
    set_stat('chi2gehrels')
    set_method('levmar')
    hdu_number = 1  #Want evnts so hdu_number = 1
    src_model_dict = {}; bkg_model_dict = {}
    obs_count = 1
    for spec_pha in spectrum_files:
        obsid_set(src_model_dict, bkg_model_dict, spec_pha,background_files[int(obs_count-1)], obs_count, redshift, n_H, Temp_guess)
        obs_count += 1
    for ob_num in range(obs_count-1):
        if deproj == False:
            group_counts(ob_num+1,grouping)
        notice_id(ob_num+1,0.5,8.0)
    fit()
    set_log_sherpa()
    set_covar_opt("sigma",1)
    covar(get_model_component('apec1').kT,get_model_component('apec1').Abundanc)
    mins = list(get_covar_results().parmins)
    maxes = list(get_covar_results().parmaxes)
    for val in range(len(mins)):
        if isFloat(mins[val]) == False:
            mins[val] = 0.0
        if isFloat(maxes[val]) == False:
            maxes[val] = 0.0
        else:
            pass
    #Get important values
    Temperature = apec1.kT.val
    Temp_min = Temperature+mins[0]
    Temp_max = Temperature+maxes[0]
    Abundance = apec1.Abundanc.val
    Ab_min = Abundance+mins[1]
    Ab_max = Abundance+maxes[1]
    #Calculate norm as average value
    Norm = 0; Norm_min = 0; Norm_max = 0
    #for id_ in src_ids:
    Norm += get_model_component('apec1').norm.val #add up values
    #get errors
    covar(get_model_component('apec1').norm)
    mins = list(get_covar_results().parmins)
    maxes = list(get_covar_results().parmaxes)
    for val in range(len(mins)):
        if isFloat(mins[val]) == False:
            mins[val] = 0.0
        if isFloat(maxes[val]) == False:
            maxes[val] = 0.0
        else:
            pass
        Norm_min += mins[0]
        Norm_max += maxes[0]
    Norm = Norm/len(spectrum_files)
    Norm_min = Norm+Norm_min/len(spectrum_files)
    Norm_max = Norm+Norm_max/len(spectrum_files)
    f = get_fit_results()
    reduced_chi_sq = f.rstat
    #---------Set up Flux Calculation----------#
    flux_calculation = sample_flux(get_model_component('apec1'), 0.01, 50.0, num=1000, confidence=90)[0]
    Flux = flux_calculation[0]
    Flux_min = flux_calculation[1]
    Flux_max = flux_calculation[2]
    reset(get_model())
    reset(get_source())
    clean()
    return Temperature,Temp_min,Temp_max,Abundance,Ab_min,Ab_max,Norm,Norm_min,Norm_max,reduced_chi_sq,Flux,Flux_min,Flux_max



def Fitting(base_directory,dir,file_name,num_files,redshift,n_H,Temp_guess,output_file:
    """
    Fit each region's spectra and create a text file containing the spectral
    fit information for each bin
    Args:
        base_directory: Path to main Directory
        dir: ObsID
        file_name: Root name of PI/PHA file
        num_files: Number of spatial bins
        redshift: Redshift of cluster
        n_H: Column density
        temp_guess: Initial temperature guess
        output_file: Text file containing each bin's spectral fit information
        bin_spec_dir: Path to extracted spectra for each bin within an ObsID

    Return:
        None
    """
    energy_min = 0.5
    energy_max = 8.0
    grouping = 5
    plot_dir = base_directory+'/FitPlots/'
    output_file = output_file.split('.')[0]
    os.chdir(base_directory)
    if plot_dir != '':
        if not os.path.exists(plot_dir):
            os.makedirs(plot_dir)
    if os.path.isfile(file_name) == True:
        os.remove(file_name) #remove it
    # Non Deprojected Fits
    file_to_write = open(output_file+".txt",'w+')
    file_to_write.write("BinNumber Temperature Temp_min Temp_max Abundance Ab_min Ab_max Norm Norm_min Norm_max ReducedChiSquare Flux Flux_min Flux_max\n")
    for i in range(num_files):
        print("Fitting model to spectrum number "+str(i+1))
        spectrum_files = []
        background_files = []
        arf_files = []
        resp_file = []
        for directory in dir:
            try:
                if num_files == 1:
                    spectrum_files.append(directory+'/repro/'+file_name+".pi")
                    background_files.append(directory+'/repro/'+file_name+"_bkg.pi")
                else:
                    spectrum_files.append(directory+'/repro/'+file_name+'_'+str(i)+".pi")
                    background_files.append(directory+'/repro/'+file_name+'_'+str(i)+"_bkg.pi")
            except:
                pass
        try:
            Temperature,Temp_min,Temp_max,Abundance,Ab_min,Ab_max,Norm,Norm_min,Norm_max,reduced_chi_sq,Flux = FitXSPEC(spectrum_files,background_files,redshift,n_H,Temp_guess,grouping,i,plot_dir)
            file_to_write.write("%i %.2E %.2E %.2E %.2E %.2E %.2E %.2E %.2E %.2E %.2E %.2E\n"%(i,Temperature,Temp_min,Temp_max,Abundance,Ab_min,Ab_max,Norm,Norm_min,Norm_max,reduced_chi_sq,Flux,Flux_min,Flux_max))
        except:
            print("No spectra was fit for bin number %i!"%i)
    file_to_write.close()
