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
from keras.backend import clear_session
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
import tensorflow as tf
from keras.models import Sequential
from keras.layers import Dense, InputLayer, Flatten, Dropout
from keras.layers.convolutional import Conv1D
from keras.layers.convolutional import MaxPooling1D
from keras.optimizers import Adam
from keras.callbacks import EarlyStopping, ReduceLROnPlateau
from keras.utils import to_categorical
from ciao_contrib.runtool import *
from astropy.io import fits
from sklearn import preprocessing
from pickle import dump, load
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
# Disable Python # WARNING:
import warnings
warnings.filterwarnings("ignore")
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'

# Check if value is a float
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
# MACHINE LEARNING VARIABLES
ML_model_name = 'ML_model_Low_NGC4636.h5'
chan_min = 0  # Minimum X-ray Channel
chan_max = 600  # Maximum X-ray Channel
activation = 'relu'  # activation function
initializer = 'he_normal'  # model initializer
input_shape = (None, chan_max-chan_min, 1)  # shape of input spectra for the input layer
num_filters = [4,16]  # number of filters in the convolutional layers
filter_length = [8,4]  # length of the filters
pool_length = 4  # length of the maxpooling
num_hidden = [256,128]  # number of nodes in the hidden layers
batch_size = 8  # number of data fed into model at once
max_epochs = 10  # maximum number of interations
lr = 0.0007  # initial learning rate
beta_1 = 0.9  # exponential decay rate  - 1st
beta_2 = 0.999  # exponential decay rate  - 2nd
optimizer_epsilon = 1e-08  # For the numerical stability
early_stopping_min_delta = 0.0001
early_stopping_patience = 4
reduce_lr_factor = 0.5
reuce_lr_epsilon = 0.009
reduce_lr_patience = 2
reduce_lr_min = 0.00008
loss_function = 'mean_squared_error'
metrics = ['accuracy', 'mae']
temp_min = 0  # Set in the training
temp_max = 8
# num_labels = (temp_max-temp_min)*10+1  # How many temperatures possible
# temp_labels = np.round(np.linspace(temp_min, temp_max, num_labels), 1)  # Actual Possible temp labels
#------------------------------FUNCTIONS---------------------------------------#
#Dynamically set source for OBSID
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
#Dynamically set source for OBSID
def obsid_set2(src_model_dict,bkg_model_dict,obsid, obs_count,redshift,nH_val,Temp_guess):
    '''
    Add two thermal emission models
    '''
    load_pha(obs_count,obsid) #Read in
    if obs_count == 1:
        src_model_dict[obsid] = xsphabs('abs'+str(obs_count)) * (xsapec('apec1_'+str(obs_count)) + xsapec('apec2_'+str(obs_count))) #set model and name
        # Change src model component values
        get_model_component('apec1_' + str(obs_count)).kT = 1
        get_model_component('apec1_' + str(obs_count)).redshift = redshift  # need to tie all together
        get_model_component('apec1_' + str(obs_count)).Abundanc = 0.3
        thaw(get_model_component('apec1_' + str(obs_count)).Abundanc)
        get_model_component('apec2_' + str(obs_count)).kT = 2.0
        get_model_component('apec2_' + str(obs_count)).redshift = redshift  # need to tie all together
        get_model_component('apec2_' + str(obs_count)).Abundanc = 0.3#get_model_component('apec1_1').Abundanc
        thaw(get_model_component('apec2_' + str(obs_count)).Abundanc)
        get_model_component('abs1').nH = nH_val  # change preset value
        freeze(get_model_component('abs1'))
    else:
        src_model_dict[obsid] = get_model_component('abs1') * (xsapec('apec1_'+str(obs_count)) + xsapec('apec2_'+str(obs_count)))
        get_model_component('apec1_'+str(obs_count)).kT = get_model_component('apec1_1').kT #link to first kT
        get_model_component('apec1_' + str(obs_count)).redshift = redshift
        get_model_component('apec1_' + str(obs_count)).Abundanc = get_model_component('apec1_1').Abundanc  # link to first kT
        get_model_component('apec2_'+str(obs_count)).kT = get_model_component('apec2_1').kT #link to first kT (second thermal)
        get_model_component('apec2_' + str(obs_count)).redshift = redshift
        get_model_component('apec2_' + str(obs_count)).Abundanc = get_model_component('apec2_1').Abundanc  # link to first kT (second thermal)

    bkg_model_dict[obsid] = xsapec('bkgApec'+str(obs_count))+get_model_component('abs1')*xsbremss('brem'+str(obs_count))
    set_source(obs_count, src_model_dict[obsid]) #set model to source
    set_bkg_model(obs_count,bkg_model_dict[obsid])
    #Change bkg model component values
    get_model_component('bkgApec' + str(obs_count)).kT = 0.18
    freeze(get_model_component('bkgApec'+str(obs_count)).kT)
    get_model_component('brem' + str(obs_count)).kT = 40.0
    freeze(get_model_component('brem' + str(obs_count)).kT)
    return None
#------------------------------------------------------------------------------#
#FitXSPEC
# Fit spectra
#   parameters:
#       spectrum_file = Name of combined spectra File
#       background_file = Name of associated background file
#       arf_file = Name of associated arf file
#       resp_file = Name of associated rmf file
#       redshift = redshift of object
#       n_H = Hydrogen Equivalent Column Density
#       Temp_guess = Guess for Temperature value
def FitXSPEC(spectrum_files,background_files,redshift,n_H,Temp_guess,grouping,spec_count,plot_dir,multi=False):
    #FIX HEADER
    set_stat('chi2gehrels')
    set_method('levmar')
    hdu_number = 1  #Want evnts so hdu_number = 1
    src_model_dict = {}; bkg_model_dict = {}
    obs_count = 1
    for spec_pha in spectrum_files:
        obsid_set(src_model_dict, bkg_model_dict, spec_pha, background_files[int(obs_count-1)], obs_count, redshift, n_H, Temp_guess)
        obs_count += 1
    for ob_num in range(obs_count-1):
        group_counts(ob_num+1,grouping)
        notice_id(ob_num+1,0.5,8.0)
    fit()
    set_log_sherpa()
    #plot("fit", 1, "fit", 2)
    #print_window(plot_dir+"%s.ps"%spec_count,['clobber','yes'])
    #set_covar_opt("sigma",3)
    #covar(get_model_component('apec1').kT,get_model_component('apec1').Abundanc)
    #with open(os.getcwd()+'/Fits/Params/%s_err_agn_%s.out'%(spec_count,agn_ct),'w+') as res_out:
    #    res_out.write(str(get_covar_results()))
    #----------Calculate min/max values---------#
    '''mins = list(get_covar_results().parmins)
    maxes = list(get_covar_results().parmaxes)
    for val in range(len(mins)):
        if isFloat(mins[val]) == False:
            mins[val] = 0.0
        if isFloat(maxes[val]) == False:
            maxes[val] = 0.0
        else:
            pass'''
    #Get important values
    Temperature = apec1.kT.val
    print(Temp_guess, Temperature)
    Temp_min = Temperature#+mins[0]
    Temp_max = Temperature#+maxes[0]
    Abundance = apec1.Abundanc.val;
    Ab_min = Abundance#+mins[1];
    Ab_max = Abundance#+maxes[1]
    #Calculate norm as average value
    Norm = 0; Norm_min = 0; Norm_max = 0
    '''for id_ in range(len(spectrum_files)):
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
    Norm_max = Norm+Norm_max/len(spectrum_files)'''
    f = get_fit_results()
    reduced_chi_sq = f.rstat
    reset(get_model())
    reset(get_source())
    clean()

    return Temperature,Temp_min,Temp_max,Abundance,Ab_min,Ab_max,Norm,Norm_min,Norm_max,reduced_chi_sq

def FitXSPEC_multi(spectrum_files,background_files,redshift,n_H,Temp_guess,grouping,spec_count,plot_dir):
    #FIX HEADER
    set_stat('cstat')
    set_method('levmar')
    hdu_number = 1  #Want evnts so hdu_number = 1
    src_model_dict = {}; bkg_model_dict = {}
    obs_count = 1
    for spec_pha in spectrum_files:
        obsid_set2(src_model_dict, bkg_model_dict, spec_pha, obs_count, redshift, n_H, Temp_guess)
        obs_count += 1
    for ob_num in range(obs_count-1):
        group_counts(ob_num+1,grouping)
        notice_id(ob_num+1,0.1,8.0)
    fit()
    Temperature1 = apec1_1.kT.val
    Temperature2 = apec2_1.kT.val
    Abundance1 = apec1_1.Abundanc.val
    Abundance2 = apec1_1.Abundanc.val
    f = get_fit_results()
    reduced_chi_sq = f.rstat
    reset(get_model())
    reset(get_source())
    clean()
    return Temperature1,Temperature2,Abundance1,Abundance2,reduced_chi_sq
#------------------------------------------------------------------------------#
def ML_model():
    # Define the ML model
    clear_session()
    model = Sequential([
        InputLayer(batch_input_shape=input_shape),
        Conv1D(kernel_initializer=initializer, activation=activation, padding="same", filters=num_filters[0], kernel_size=filter_length[0]),
        Conv1D(kernel_initializer=initializer, activation=activation, padding="same", filters=num_filters[1], kernel_size=filter_length[1]),
        MaxPooling1D(pool_size=pool_length),
        Flatten(),
        Dropout(0.2),
        Dense(units=num_hidden[0], kernel_initializer=initializer, activation=activation),
        Dense(units=num_hidden[1], kernel_initializer=initializer, activation=activation),
        Dense(1, activation="linear"),
    ])
    optimizer = Adam(lr=lr, beta_1=beta_1, beta_2=beta_2, epsilon=optimizer_epsilon, decay=0.0)
    early_stopping = EarlyStopping(monitor='val_loss', min_delta=early_stopping_min_delta,
                                           patience=early_stopping_patience, verbose=2, mode='min')
    reduce_lr = ReduceLROnPlateau(monitor='loss', factor=0.5, epsilon=reuce_lr_epsilon,
                                      patience=reduce_lr_patience, min_lr=reduce_lr_min, mode='min', verbose=2)
    model.compile(optimizer=optimizer, loss=loss_function, metrics=metrics)
    # Load the ML model
    model.load_weights(ML_model_name)
    scaler = load(open('Scaler_model_Low.pkl', 'rb'))
    #scaler = None
    return model, scaler
#------------------------------------------------------------------------------#
def ML_fit(spectrum_files,model,scaler,ML_temps,i):
    '''
    Subroutine to apply our machine learning algorithm to the combined and
    normalized spectra in order to estimate a temperature
    '''
    temp_guess = 0.0
    try:
    # Combine spectra
        spectrum_files_str = ''
        for spectrum_file in spectrum_files:
            spectrum_files_str += spectrum_file+','
        spectrum_files_str = spectrum_files_str[:-1]
        combine_spectra.punlearn()
        combine_spectra.src_spectra = spectrum_files_str
        combine_spectra.outroot = 'combined_'+str(i)
        combine_spectra.clobber = True
        combine_spectra()
        # Get combined in correct python form using astropy
        spec_ = fits.open('combined_'+str(i)+'_src.pi')[1].data
        #spec_ = fits.open(spectrum_files[1])[1].data
        # counts = [chan[3] for chan in spec_[chan_min:chan_max]]  # We only take the first 600 channels (i.e. <8 keV)
        vals = list(spec_[chan_min:chan_max][:])
        counts = np.array([list(list(zip(*vals))[3])])
        # Now to normalize the counts
        #scaled_counts = np.array(scaler.transform(counts))
        counts_max = np.max(counts)
        scaled_counts = np.array([count/counts_max for count in counts])
        scaled_counts = scaled_counts.reshape(scaled_counts.shape[0], scaled_counts.shape[1], 1)
        predictions = model.predict(scaled_counts)
        temp_guess = predictions[0][0]
        # met_guess = predictions[0,1]
        # most_prob_ind = np.argmax(predictions)
        # temp_guess = temp_labels[most_prob_ind]  # returns string
        ML_temps.write('%i %s\n'%(i, temp_guess))  # , met_guess))
    except:
        pass
    return temp_guess
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
def fit_loop(dir,bin_spec_dir,file_name,redshift,n_H,Temp_guess,grouping,plot_dir,base_directory,model,scaler,ML_temps,i):
    '''
    Parallelized loop to fit spectra using sherpa and xspec
    '''
    os.chdir(base_directory)
    spectrum_files = []
    background_files = []
    arf_files = []
    resp_file = []
    # Collect data to be fit
    for directory in dir:
        try:
            spectrum_files.append(directory+'/repro/'+bin_spec_dir+file_name+"_"+str(i)+".pi")
            background_files.append(directory+'/repro/'+bin_spec_dir+file_name+"_"+str(i)+"_bkg.pi")
        except:
            pass
    # Obtain initial guess by passing data to a machine learning algorithm
    Temp_guess = ML_fit(spectrum_files,model,scaler,ML_temps,i)
    # Fit data in parallel using standard sherpa fitting processes
    try:
    #if multi.lower() == 'false':
        Temperature,Temp_min,Temp_max,Abundance,Ab_min,Ab_max,Norm,Norm_min,Norm_max,reduced_chi_sq = FitXSPEC(spectrum_files,background_files,redshift,n_H,Temp_guess,grouping,i,plot_dir)
        with open('temp_'+str(i)+'.txt','w+') as out_temp:
            out_temp.write("%i %f %f %f %f %f %f %f %f %f %f\n"%(i,Temperature,Temp_min,Temp_max,Abundance,Ab_min,Ab_max,Norm,Norm_min,Norm_max,reduced_chi_sq))
    #else:
    #    Temperature1,Temperature2,Abundance1,Abundance2,reduced_chi_sq = FitXSPEC_multi(spectrum_files,background_files,redshift,n_H,Temp_guess,grouping,i,plot_dir)
    #    out_temp.write("%i %f %f %f %f %f \n"%(i, Temperature1, Temperature2, Abundance1, Abundance2))
    except:
        with open('temp_'+str(i)+'.txt','w+') as out_temp:
            out_temp.write("%i %f %f %f %f %f %f %f %f %f %f\n"%(i,0,0,0,0,0,0,0,0,0,0))
        print("No spectra was fit")
    return None
#--------------------------------------------------------------------#
def concat_temp_data(num_spec, output_file):
    file_to_write = open(output_file+".txt",'w+')
    file_to_write.write("BinNumber Temperature Temp_min Temp_max Abundance Ab_min Ab_max Norm Norm_min Norm_max ReducedChiSquare \n")
    # Get the data for each temporary file and then delete
    for spec_i in range(num_spec):
        with open('temp_'+str(spec_i)+'.txt', 'r') as temp_file:
            for line in temp_file.readlines():
                file_to_write.write(line)
        os.remove(os.path.join(os.getcwd(),'temp_'+str(spec_i)+'.txt'))
    file_to_write.close()
    # Delete the combined spectra...
    '''for fname in os.listdir(os.getcwd()):
        if fname.startswith("combined"):
            os.remove(os.path.join(os.getcwd(), fname))'''
    return None
#--------------------------------------------------------------------#

#--------------------------------------------------------------------#
#PrimeFitting
# Step through spectra to fit
#   parameters:
#       dir = main Directory
#       file_name = FIlename of PI/PHA spectrum
#       output_file = Filename for output containing temperature information
#       num_files = number of bins
#       redshift = redshift of object
#       n_H = Hydrogen Equivalent Column Density
#       Temp_guess = Guess for Temperature value
def PrimeFitting(base_directory,dir,file_name,num_files,redshift,n_H,Temp_guess,output_file,bin_spec_dir,multi=False):
    energy_min = 0.5
    energy_max = 8.0
    grouping = 1
    plot_dir = base_directory+'/FitPlots/'
    output_file = output_file.split('.')[0]
    os.chdir(base_directory)
    if plot_dir != '':
        if not os.path.exists(plot_dir):
            os.makedirs(plot_dir)
    if os.path.isfile(file_name) == True:
        os.remove(file_name) #remove it
    # Define the ML Model
    model, scaler = ML_model()
    ML_temps = open(base_directory+'/ml_temps.txt','w+')
    ML_temps.write('Spectrum Temp\n')
    # Fit data in parallel loop
    Parallel(n_jobs=1,prefer="processes")(delayed(fit_loop)(dir,bin_spec_dir,file_name,redshift,n_H,Temp_guess,grouping,plot_dir,base_directory,model,scaler,ML_temps,i) for i in tqdm(range(num_files)))
    concat_temp_data(num_files, output_file)
    ML_temps.close()
