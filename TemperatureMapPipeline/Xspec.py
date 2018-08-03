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
from astropy.io import fits
import os
import subprocess
import matplotlib.pyplot as plt
from xspec import *
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
def FitXSPEC(spectrum_file,background_file, arf_file,resp_file,redshift,n_H,Temp_guess,spec_count,plot_dir):
    #FIX HEADER
    hdu_number = 1  #Want evnts so hdu_number = 1
    header = fits.getheader(spectrum_file, hdu_number)
    header['BACKFILE'] = background_file
    header['RESPFILE'] = resp_file
    header['ANCRFILE'] = arf_file
    #Kill the output
    Xset.chatter = 0
    #XSPEC GOODNESS
    s1 = Spectrum(spectrum_file)
    m1 = Model("zphabs*(apec+pow)")
    c1 = m1.zphabs
    c2 = m1.apec

    c1.nH = n_H
    c1.nH.frozen = True
    c1.Redshift = redshift
    c2.kT = Temp_guess
    c2.Redshift = redshift


    s1.ignore('0-0.3 10.0-**')
    #AllData.ignore('bad')
    Plot('data')
    Fit.nIterations = 1000
    Fit.statMethod = "cstat"
    Fit.perform()
    Plot('data')
    reduced_chi_sq = Fit.testStatistic/Fit.dof
    Temperature = c2.kT.values[0]

    #PLOT

    Plot.device = '/Null'
    Plot.xAxis = "keV"
    Plot('data')
    Fit.perform()
    # Get coordinates from plot:
    chans = Plot.x()
    rates = Plot.y()
    folded = Plot.model()
    # Plot using Matplotlib:
    if plot_dir != None:
        plt.plot(chans, rates, '+', chans, folded)
        plt.xlabel('KeV')
        plt.yscale('log')
        plt.ylabel(r'counts/cm$^2$/sec/chan')
        plt.title(r'Spectra '+str(spec_count))
        plt.savefig(plot_dir+'spectrum_'+str(spec_count)+'.pdf')

    AllData.clear()
    AllModels.clear()

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
    os.chdir(dir)
    if plot_dir != '':
        if not os.path.exists(plot_dir):
            os.makedirs(plot_dir)
    if os.path.isfile(file_name) == True:
        os.remove(file_name) #remove it
    file_to_write = open(output_file+".txt",'w+')
    file_to_write.write("BinNumber Temperature ReducedChiSquare \n")
    for i in range(num_files):
        print("Fitting model to spectrum number "+str(i))
        spectrum_file = file_name+"_"+str(i)+".pi"
        background_file = file_name+"_"+str(i)+"_bkg.pi"
        arf_file = file_name+"_"+str(i)+".arf"
        resp_file = file_name+"_"+str(i)+".rmf"
        Temperature,reduced_chi_sq = FitXSPEC(spectrum_file,background_file, arf_file,resp_file,redshift,n_H,Temp_guess,i,plot_dir)
        file_to_write.write(str(i) + " " + str(Temperature) + " " + str(reduced_chi_sq) + " \n")
    file_to_write.close()

def main():
    dir = '/home/crhea/Documents/TestData/12036/repro/binned/'
    file_name = 'simple'
    output_file = 'Temperature_bins'
    num_files = 17
    redshift = 0.01790
    n_H = 1.38e-1
    Temp_guess = 5
    plot_dir = '/home/crhea/Desktop/spec/'
    print("--------------------------------FITTING DATA------------------------#")
    PrimeFitting(dir,file_name,output_file,num_files,redshift,n_H,Temp_guess,plot_dir)
main()
