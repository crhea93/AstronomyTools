'''
------------------------------------------------------
GOAL:
    Step through bins (spectra) and calculate the temperature value
    of each bin using XSPEC
------------------------------------------------------
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
from multiprocessing import Process, JoinableQueue
from joblib import Parallel, delayed
from tqdm import tqdm
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
def obsid_set(src_model_dict,bkg_model_dict,obsid, bkg_spec,obs_count,redshift,nH_val,Temp_guess):
    load_pha(obs_count,obsid) #Read in
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

    bkg_model_dict[obsid] = xsapec('bkgApec'+str(obs_count))+get_model_component('abs1')*xsbremss('brem'+str(obs_count))
    set_bkg(obs_count, unpack_pha(bkg_spec))
    set_source(obs_count, src_model_dict[obsid]) #set model to source
    set_bkg_model(obs_count,bkg_model_dict[obsid])
    #Change bkg model component values
    get_model_component('bkgApec' + str(obs_count)).kT = 0.18
    freeze(get_model_component('bkgApec'+str(obs_count)).kT)
    get_model_component('brem' + str(obs_count)).kT = 40.0
    freeze(get_model_component('brem' + str(obs_count)).kT)
    return None
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
def FitXSPEC(spectrum_files,background_files,redshift,n_H,temp_guess,grouping,spec_count,plot_dir):

    set_stat('chi2gehrels')
    set_method('levmar')
    hdu_number = 1  #Want evnts so hdu_number = 1
    src_model_dict = {}; bkg_model_dict = {}
    obs_count = 1
    for spec_pha in spectrum_files:  # Set up model for spectra and background
        obsid_set(src_model_dict, bkg_model_dict, spec_pha, background_files[int(obs_count-1)], obs_count, redshift, n_H, temp_guess)
        obs_count += 1
    for ob_num in range(obs_count-1):  # Group and set range
        group_counts(ob_num+1,grouping)
        notice_id(ob_num+1,0.5,8.0)
    fit()  # FIT!
    set_log_sherpa()
    set_covar_opt("sigma",2)
    covar(get_model_component('apec1').kT,get_model_component('apec1').Abundanc)
    #----------Calculate min/max values---------#
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
    Abundance = apec1.Abundanc.val;
    Ab_min = Abundance+mins[1];
    Ab_max = Abundance+maxes[1]
    #Calculate norm as average value
    Norm = 0; Norm_min = 0; Norm_max = 0
    for id_ in range(len(spectrum_files)):
        Norm += get_model_component('apec'+str(id_+1)).norm.val #add up values
        #get errors
        covar(get_model_component('apec'+str(id_+1)).norm)
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

    return Temperature,Temp_min,Temp_max,Abundance,Ab_min,Ab_max,Norm,Norm_min,Norm_max,reduced_chi_sq, Flux, Flux_min, Flux_max

#---------------------------------------------------------#
def fit_loop(dir,bin_spec_dir,file_name,redshift,n_H,temp_guess,grouping,plot_dir,base_directory,i):
    """
    Loop for fitting spectra.
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
        grouping: Number of counts to bin in sherpa fit
        i: Bin number
        plot_dir: Path to plots

    Return:
        None
    """
    os.chdir(base_directory)
    spectrum_files = []
    background_files = []
    for directory in dir:  # Step through each ObsID
        try:  # Collect spectrum if it exists
            spectrum_files.append(directory+'/repro/'+bin_spec_dir+file_name+"_"+str(i)+".pi")
            background_files.append(directory+'/repro/'+bin_spec_dir+file_name+"_"+str(i)+"_bkg.pi")
        except:
            pass
    Temperature,Temp_min,Temp_max,Abundance,Ab_min,Ab_max,Norm,Norm_min,Norm_max,reduced_chi_sq,Flux,Flux_min,Flux_max =\
         FitXSPEC(spectrum_files,background_files,redshift,n_H,temp_guess,grouping,i,plot_dir)
    with open('temp_'+str(i)+'.txt','w+') as out_temp:
        out_temp.write("%i %.2E %.2E %.2E %.2E %.2E %.2E %.2E %.2E %.2E %.2E %.2E %.2E %.2E\n"% \
                (i,Temperature,Temp_min,Temp_max,Abundance,Ab_min,Ab_max,Norm,Norm_min,Norm_max,reduced_chi_sq,Flux,Flux_min,Flux_max))

    return None
#--------------------------------------------------------------------#
def concat_temp_data(num_spec, output_file):
    """
    Concatenate temperature information. Used if running on several processors
    Args:
        num_spec: Number of spectra (or bins)
        output_file: File name of final concatenated temperature data
    """
    file_to_write = open(output_file+".txt",'w+')
    file_to_write.write("BinNumber Temperature Temp_min Temp_max Abundance Ab_min Ab_max Norm Norm_min Norm_max ReducedChiSquare Flux Flux_min Flux_max\n")
    # Get the data for each temporary file and then delete
    for spec_i in range(num_spec):
        with open('temp_'+str(spec_i)+'.txt', 'r') as temp_file:
            for line in temp_file.readlines():
                file_to_write.write(line)
    file_to_write.close()
    return None
#--------------------------------------------------------------------#
#--------------------------------------------------------------------#
def Fitting(base_directory,dir,file_name,num_files,redshift,n_H,temp_guess,output_file,bin_spec_dir):
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
    # Make sure plotting directory exists
    if plot_dir != '':
        if not os.path.exists(plot_dir):
            os.makedirs(plot_dir)
    if os.path.isfile(file_name) == True:
        os.remove(file_name) #remove it
    Parallel(n_jobs=1,prefer="processes")(delayed(fit_loop)\
            (dir,bin_spec_dir,file_name,redshift,n_H,temp_guess,grouping,plot_dir,
                base_directory,i) for i in tqdm(range(num_files)))
    concat_temp_data(num_files, output_file)
