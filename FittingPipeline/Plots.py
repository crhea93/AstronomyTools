'''
This program will plot the final temperature map using the temperatures calculated from sherpa and the binnings from WVT
'''
import os
import numpy as np
import matplotlib as mpl
from astropy.io import fits
from astropy.wcs import WCS
import matplotlib.pyplot as plt
import matplotlib.colorbar as cbar
from matplotlib.collections import PatchCollection
#-----------------------------------INPUTS------------------------------------------#
num_ticks = 5
#-------------------------------------------------#
def plot_data(temp_data,reg_dir,regions,base_dir):
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
    stat = []
    # Read in temperature data
    for line in temp_d:
        print(line)
        temps.append(float(line.split(" ")[1]))
        temps_min.append(float(line.split(" ")[2])-float(line.split(" ")[1]))
        temps_max.append(float(line.split(" ")[3])-float(line.split(" ")[1]))
        abs.append(float(line.split(" ")[4]))
        abs_min.append(float(line.split(" ")[5]))
        abs_max.append(float(line.split(" ")[6]))
        stat.append(float(line.split(" ")[3]))
    temp_d.close()
    # Read in radius data
    for region in regions:
        with open(reg_dir+region+'.reg') as reg_:
            reg_data = reg_.readlines()[3].split(')')[0].split('(')[1]
            r_in.append(reg_data.split(',')[3])
            r_out.append(reg_data.split(',')[2])
            r_mid.append(float(reg_data.split(',')[2])-float(reg_data.split(',')[3])/2)

    # Plotting
    plt.errorbar(r_mid,temps,yerr=(temps_min,temps_max),fmt='o')
    plt.show()
    return None
#-------------------------------------------------#
#-------------------------------------------------#
