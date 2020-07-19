'''
Calculate the Hardness Ratio for each bin and plot
HR =(H-S)/(H+S)

Where H is the number of counts in the hard band (2.5-8.0 keV) and
S is the number of counts in the soft band (0.5-2.0 keV)


We must have already created a set of region files for each bin
'''
#----------------------IMPORTS-----------------------------#
import os
import shutil
import numpy as np
import matplotlib as mpl
from astropy.io import fits
from astropy.wcs import WCS
import matplotlib.pyplot as plt
import matplotlib.colorbar as cbar
from ciao_contrib.runtool import *
from crates_contrib.utils import *
from matplotlib.collections import PatchCollection

num_ticks = 5
#----------------------CLASSES-----------------------------#
class Bin:
    def __init__(self, number):
        self.bin_number = number
        self.pixels = []
        self.hr = 0
    def add_pixel(self,pixel):
        self.pixels.append(pixel)
    def add_hr(self,hr):
        self.hr = hr

class Pixel:
    def __init__(self, number, pix_x, pix_y):
        self.pix_number = number
        self.pix_x = pix_x
        self.pix_y = pix_y
#---------------------FUNCTIONS----------------------------#

def calc_vals(wcs_fits,reg_file,region):
    #Hard then soft
    dmextract.infile = wcs_fits+"[energy=2500:8000][bin sky=region("+region+".reg)]"
    dmextract.outfile = region+'_Hardcounts.fits'
    dmextract.opt = 'generic'
    dmextract.clobber = True
    dmextract()
    dmstat.infile = region+'_Hardcounts.fits[cols counts]'
    dmstat()
    hard = float(dmstat.out_sum)

    dmextract.infile = wcs_fits+"[energy=500:2000][bin sky=region("+region+".reg)]"
    dmextract.outfile = region+'_Softcounts.fits'
    dmextract.opt = 'generic'
    dmextract.clobber = True
    dmextract()
    dmstat.infile = region+'_Softcounts.fits[cols counts]'
    dmstat()
    soft = float(dmstat.out_sum)

    return hard,soft

def hardness_plot(bin_file,hr_file,file_dir,filename,color_map,stn,wcs_fits):
    '''
    Plot hardness ratio
    '''
    #Get WCS information from header of original image
    hdu = fits.open(wcs_fits)[0]
    wcs = WCS(hdu.header)
    Bins, x_min, x_max, y_min, y_max = read_in(bin_file,hr_file)
    fig = plt.figure()
    fig.set_size_inches(7, 7)
    ax = plt.axes(xlim=(x_min,x_max), ylim=(y_min,y_max),projection=wcs)
    #plt.axes(projection=wcs)
    N = len(Bins)
    cmap = mpl.cm.get_cmap(color_map)
    max_hr = max([bin.hr for bin in Bins])
    median_hr = np.median([bin.hr for bin in Bins])
    std_hr = np.std([bin.hr for bin in Bins])
    Bins_flush = [] #Bins that passed the sigma screening
    Bins_fail = [] #Bins that failed the sigma screening
    step_val = 1
    max_hr = max([bin.hr for bin in Bins])
    min_hr = min([bin.hr for bin in Bins])
    hr_list = []
    hr_norm_list = []
    bin_nums = []
    #Set up color map
    for bin in Bins:
        hr_norm = (bin.hr - min_hr) / (max_hr - min_hr)
        hr_list.append(bin.hr)
        hr_norm_list.append(hr_norm)

    colors = cmap(hr_norm_list)


    #Create rectanges
    rect_step = 0
    for bin in Bins:
        patches = []
        bin_nums.append(bin.bin_number)
        #Color based on signal-to-noise value
        c = colors[rect_step]
        for pixel in bin.pixels:
            x_coord = pixel.pix_x
            y_coord = pixel.pix_y
            #Shift because x_coord,y_coord are the center points
            rectangle = plt.Rectangle((x_coord,y_coord),1,1, color=c)
            #ax.add_patch(rectangle)
            patches.append(rectangle)
            ax.add_patch(rectangle)
        rect_step += 1
    colors = np.linspace(min(hr_norm_list),max(hr_norm_list),N)


    plt.xlabel("RA")
    plt.ylabel("DEC")
    norm = mpl.colors.Normalize(min(hr_norm_list),max(hr_norm_list))
    cax, _ = cbar.make_axes(ax)
    cb2 = cbar.ColorbarBase(cax, cmap=cmap, norm=norm)
    cb2.set_label('Hardness Ratio')
    tick_list = np.linspace(min(hr_norm_list),max(hr_norm_list),num_ticks)
    ticklabel_list = np.linspace(min_hr,max_hr,num_ticks)
    print(max_hr)
    ticklabel_list = [np.round(val,2) for val in ticklabel_list]
    cb2.set_ticks(tick_list)
    cb2.set_ticklabels(ticklabel_list)
    cb2.update_ticks()
    plt.savefig(file_dir+'/'+filename+"_"+str(stn)+"_HR.png")
    return ax

#-------------------------------------------------#
#-------------------------------------------------#
def read_in(bin_data,hardness_data):
    bin_d = open(bin_data); next(bin_d); next(bin_d) #first two lines are header info
    hardness_d = open(hardness_data); next(hardness_d)
    bins = []
    pixels = []
    pixel_num = 0
    #Create bins and add pixels
    for line in bin_d:
        if int(line.split(" ")[2]) not in [bin.bin_number for bin in bins]:
            bins.append(Bin(int(line.split(" ")[2])))
        pixel_ = Pixel(pixel_num,int(line.split(" ")[0]),int(line.split(" ")[1]))
        pixels.append(pixel_)
        bins[int(int(line.split(" ")[2]))].add_pixel(pixel_)
        pixel_num += 1
    #Add temperatures and Reduced Chi Squares to bins
    for line in hardness_d:
        bins[int(int(line.split(",")[0]))].add_hr(float(line.split(",")[3]))
    hardness_d.close()
    bin_d.close()
    min_x = np.min([pixel.pix_x for pixel in pixels])
    max_x = np.max([pixel.pix_x for pixel in pixels])
    min_y = np.min([pixel.pix_y for pixel in pixels])
    max_y = np.max([pixel.pix_y for pixel in pixels])
    return bins, min_x, max_x, min_y, max_y



#-----------------------------MAIN FUNCTION----------------------------#
def hardness_ratio(base_dir,name,obsid_0,stn_target,num_bins,output_dir,color_map,wcs_fits,bin_file):
    '''
    Calculate hardness ratio and plot
    params:
    base_dir -- path to base directory
    name -- name of cluster
    obsid_0 -- first obsid (doesn't actually matter which)
    stn_target -- target signal-to-noise ratio
    num_bins -- number of bins from WVT
    output_dir -- name of folder which contains region files for bins
    '''
    if os.path.exists(base_dir+'/HR'):
        os.chdir(base_dir+'/HR')
    else:
        os.mkdir(base_dir+'/HR')
        os.chdir(base_dir+'/HR')
    hr_file = open('HR_'+str(stn_target)+'.txt','w+')
    hr_file.write('Bin,Hard_Counts,Soft_Counts,HR\n')

    #Calculate for each bin
    for bin_i in range(int(num_bins)):
        reg_file = base_dir+'/'+obsid_0+'/repro/'+output_dir+str(bin_i)+'.reg'
        shutil.copy(reg_file,str(bin_i)+'.reg')
        reg_file = str(bin_i)+'.reg'
        hard,soft = calc_vals(wcs_fits,reg_file,str(bin_i))
        HR = (hard-soft)/(hard+soft)
        hr_file.write(str(bin_i)+','+str(hard)+','+str(soft)+','+str(HR)+'\n')
    hr_file.close()
    hardness_plot(bin_file,'HR_'+str(stn_target)+'.txt',base_dir,name,color_map,stn_target,wcs_fits)
    return None
