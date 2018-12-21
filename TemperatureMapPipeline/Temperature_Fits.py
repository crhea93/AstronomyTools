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

#TURN OFF ON-SCREEN OUTPUT FROM SHERPA
import logging
logger = logging.getLogger("sherpa")
logger.setLevel(logging.WARN)
logger.setLevel(logging.ERROR)
#------------------------------INPUTS------------------------------------------#
base_directory = '/home/carterrhea/Documents/'
fold_ext = 'repro/binned'
dir = [base_directory+'323/'+fold_ext,base_directory+'324/'+fold_ext]
file_name = 'center'
output_file = 'Temp_bin'
num_files = 177
redshift = 0.003129
n_H = 1.91-2
Temp_guess = 1
energy_min = 0.5
energy_max = 8.0
grouping = 10
statistic = 'chi2gehrels'
plot_dir = base_directory+'FitPlots/'
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
def FitXSPEC(spectrum_files,background_files,redshift,n_H,Temp_guess,spec_count,plot_dir):
    #FIX HEADER
    set_stat(statistic)
    hdu_number = 1  #Want evnts so hdu_number = 1
    count = 0
    pha1 = load_pha(1, spectrum_files[0], use_errors=True)
    pha2 = load_pha(2, spectrum_files[1], use_errors=True)
    bkg1 = load_bkg(1, background_files[0])
    bkg2 = load_bkg(2, background_files[1])
    ignore(0,energy_min)
    ignore(energy_max,)
    group_counts(1,grouping);group_counts(2,grouping)
    #Set source with background
    set_source(1 , xsphabs.abs1*(xsapec.apec1))
    set_source(2 , abs1*(xsapec.apec2))
    abs1.nH = 1.91e-2
    freeze(abs1.nH)
    apec1.kT = Temp_guess
    apec2.kT = apec1.kT
    apec1.redshift = redshift
    apec2.redshift = redshift
    fit()
    plot("fit", 1, "fit", 2)
    print_window(plot_dir+"%s.ps"%spec_count,['clobber','yes'])
    Temperature = apec1.kT.val
    f = get_fit_results()
    reduced_chi_sq = f.rstat
    reset(get_model())
    reset(get_source())
    clean()
    return Temperature,reduced_chi_sq

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
def PrimeFitting(dir,file_name,output_file,num_files,redshift,n_H,Temp_guess,plot_dir):
    os.chdir(base_directory)
    if plot_dir != '':
        if not os.path.exists(plot_dir):
            os.makedirs(plot_dir)
    if os.path.isfile(file_name) == True:
        os.remove(file_name) #remove it
    file_to_write = open(output_file+".txt",'w+')
    file_to_write.write("BinNumber Temperature ReducedChiSquare \n")

    for i in range(num_files):
        print("Fitting model to spectrum number "+str(i))
        spectrum_files = []
        background_files = []
        arf_files = []
        resp_file = []
        for directory in dir:
            spectrum_files.append(directory+'/'+file_name+"_"+str(i)+".pi")
            background_files.append(directory+'/'+file_name+"_"+str(i)+"_bkg.pi")
            #resp_file = file_name+"_"+str(i)+".rmf"
        Temperature,reduced_chi_sq = FitXSPEC(spectrum_files,background_files,redshift,n_H,Temp_guess,i,plot_dir)
        file_to_write.write(str(i) + " " + str(Temperature) + " " + str(reduced_chi_sq) + " \n")
    file_to_write.close()

def main():
    print("--------------------------------FITTING DATA------------------------#")
    PrimeFitting(dir,file_name,output_file,num_files,redshift,n_H,Temp_guess,plot_dir)
main()
