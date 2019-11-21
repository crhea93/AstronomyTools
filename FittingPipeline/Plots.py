'''
This program will plot the final temperature map using the temperatures calculated from sherpa and the binnings from WVT
'''
import os
import numpy as np
import matplotlib as mpl
from astropy.io import fits
from astropy.wcs import WCS
from LSCalc import ls_calc
import matplotlib.pyplot as plt
import matplotlib.colorbar as cbar
from matplotlib.collections import PatchCollection
from Classes import Annulus
import matplotlib
matplotlib.rcParams.update({'errorbar.capsize': 2})
#-----------------------------------INPUTS------------------------------------------#
num_ticks = 5
#-------------------------------------------------#
def plot_data(temp_data,reg_dir,regions,base_dir,outfile_ext,redshift):
    temp_d = open(temp_data); next(temp_d)
    bins = []
    pixels = []
    pixel_num = 0
    #Create bins and add pixels
    #Add temperatures and Reduced Chi Squares to bins
    r_in = []
    r_out = []
    r_mid = []
    temps = []
    temps_min = []
    temps_max = []
    abs = []
    abs_min = []
    abs_max = []
    norm = []
    norm_min = []
    norm_max = []
    tcool = []
    tcool_min = []
    tcool_max = []
    press = []
    press_min = []
    press_max = []
    dens = []
    dens_min = []
    dens_max = []
    stat = []
    flux = []
    # Read in data
    for line in temp_d:
        temps.append(float(line.split(" ")[1]))
        temps_min.append(float(line.split(" ")[1])-float(line.split(" ")[2]))
        temps_max.append(float(line.split(" ")[3])-float(line.split(" ")[1]))
        abs.append(float(line.split(" ")[4]))
        abs_min.append(float(line.split(" ")[4])-float(line.split(" ")[5]))
        abs_max.append(float(line.split(" ")[6])-float(line.split(" ")[4]))
        stat.append(float(line.split(" ")[10]))
        norm.append(float(line.split(" ")[7]))
        norm_min.append(float(line.split(" ")[7])-float(line.split(" ")[8]))
        norm_max.append(float(line.split(" ")[9])-float(line.split(" ")[7]))
        flux.append(float(line.split(" ")[11]))
    temp_d.close()
    # Read in radius data
    for region in regions:
        with open(reg_dir+region+'.reg') as reg_:
            reg_data = reg_.readlines()[3].split(')')[0].split('(')[1]
            r_in_ = ls_calc(redshift,float(reg_data.split(',')[3].strip('"')))
            r_in.append(r_in_)
            r_out_ = ls_calc(redshift,float(reg_data.split(',')[2].strip('"')))
            r_out.append(r_out_)
            r_mid.append((r_in_+r_out_)/2)
    # Pass to Annuli class for the rest of the calculations
    Annuli = []
    for i in range(len(temps)):
        Annulus_ = Annulus(r_in[i],r_out[i])
        Annulus_.add_fit_data(temps[i],temps_min[i],temps_max[i],abs[i],abs_min[i],abs_max[i],norm[i],norm_min[i],norm_max[i],flux[i],stat[i],False,False,redshift)
        tcool.append(Annulus_.t_cool[1]); tcool_min.append(Annulus_.t_cool[0]); tcool_max.append(Annulus_.t_cool[2])
        press.append(Annulus_.press[1]); press_min.append(Annulus_.press[0]); press_max.append(Annulus_.press[2])
        dens.append(Annulus_.dens[1]); dens_min.append(Annulus_.dens[0]); dens_max.append(Annulus_.dens[2])
    #------------------------Plotting---------------------------#
    # TEMPERATURE
    plt.errorbar(r_mid,temps,yerr=[temps_min,temps_max], fmt='o')
    plt.xlabel('Radius (kpc)')
    plt.ylabel('Temperature (keV)')
    plt.savefig(base_dir+'/'+outfile_ext+'_temp.png')
    plt.clf()
    # Abundance
    plt.errorbar(r_mid,abs,yerr=[abs_min,abs_max], fmt='o')
    plt.xlabel('Radius (kpc)')
    plt.ylabel(r'Abundance (Z$_{\odot}$)')
    plt.savefig(base_dir+'/'+outfile_ext+'_abund.png')
    plt.clf()
    # Cooling Time
    plt.errorbar(r_mid,tcool,yerr=[tcool_min,tcool_max], fmt='o')
    plt.xlabel('Radius (kpc)')
    plt.ylabel(r'Cooling Time (Gyr)')
    plt.yscale('log')
    plt.savefig(base_dir+'/'+outfile_ext+'_Tcool.png')
    plt.clf()
    # Pressure
    plt.errorbar(r_mid,press,yerr=[press_min,press_max], fmt='o')
    plt.xlabel('Radius (kpc)')
    plt.ylabel(r'Pressure (keV cm$^{-3}$)')
    plt.yscale('log')
    plt.savefig(base_dir+'/'+outfile_ext+'_press.png')
    plt.clf()
    # Density
    plt.errorbar(r_mid,dens,yerr=[dens_min,dens_max], fmt='o')
    plt.xlabel('Radius (kpc)')
    plt.ylabel(r'Density (cm$^{-3}$)')
    plt.yscale('log')
    plt.savefig(base_dir+'/'+outfile_ext+'_dens.png')
    plt.clf()
    # Save T_cool and Pressure info
    with open(base_dir+'/AdditionalParams_'+outfile_ext+'.txt', 'w+') as file_out:
        file_out.write('BinNumber R_in R_out Tcool Tcool_min Tcool_max Press Press_min Press_max Dens Dens_min Dens_max\n')
        for i in range(len(temps)):
            file_out.write('%i %f %f %.2E %.2E %.2E %.2E %.2E %.2E %.2E %.2E %.2E\n'%(i,r_in[i],r_out[i],tcool[i],tcool_min[i],tcool_max[i],press[i],press_min[i],press_max[i],dens[i],dens_min[i],dens_max[i]))
    return None
#-------------------------------------------------#
#-------------------------------------------------#
