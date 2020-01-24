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
#-----------------------------------------------------------------------------------#
#-------------------------------------------------#
#-------------------------------------------------#
class Bin:
    def __init__(self, number):
        self.bin_number = number
        self.pixels = []
        self.temp = 0
        self.stat = 0
        self.abund = 0
    def add_pixel(self,pixel):
        self.pixels.append(pixel)
    def add_temp(self,temp):
        self.temp = temp
    def add_stat(self,stat):
        self.stat = stat
    def add_abund(self, abund):
        self.abund = abund


class Pixel:
    def __init__(self, number, pix_x, pix_y):
        self.pix_number = number
        self.pix_x = pix_x
        self.pix_y = pix_y
#-------------------------------------------------#
#-------------------------------------------------#
# Plot plot_Bins
#   parameters:
#       bin_file - WVT bin data file
#       temp_file - temperature file created during fits
#       file_dir - Full Path to location of new image
#       fileame - name of file to be printed
def plot_Ab(bin_file,temp_file,file_dir,filename,color_map,stn,wcs_image):
    #Get WCS information from header of original image
    hdu = fits.open(wcs_image)[0]
    wcs = WCS(hdu.header)
    Bins, x_min, x_max, y_min, y_max = read_in(bin_file,temp_file)
    fig = plt.figure()
    fig.set_size_inches(7, 7)
    ax = plt.axes(xlim=(x_min,x_max), ylim=(y_min,y_max),projection=wcs)
    #plt.axes(projection=wcs)
    N = len(Bins)
    cmap = mpl.cm.get_cmap(color_map)
    max_temp = max([bin.abund for bin in Bins])
    median_temp = np.median([bin.abund for bin in Bins])
    std_temp = np.std([bin.abund for bin in Bins])
    Bins_flush = [] #Bins that passed the sigma screening
    Bins_fail = [] #Bins that failed the sigma screening
    step_val = 1
    for bin in Bins:
        if bin.abund < 2.5:#median_temp +step_val*std_temp and bin.temp > median_temp - step_val*std_temp:
            Bins_flush.append(bin)
        else:
            Bins_fail.append(bin)
    max_temp = max([bin.abund for bin in Bins_flush])
    temp_list = []
    temp_norm_list = []
    bin_nums = []
    #Set up color map
    for bin in Bins_flush:
        temp_norm = bin.abund / max_temp
        temp_list.append(bin.abund)
        temp_norm_list.append(temp_norm)

    colors = cmap(temp_norm_list)


    #Create rectanges
    rect_step = 0
    for bin in Bins_flush:
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
    colors = np.linspace(min(temp_norm_list),max(temp_norm_list),N)

    #Set failed bins to black
    for bin in Bins_fail:
        patches = []
        bin_nums.append(bin.bin_number)
        #Color based on signal-to-noise value
        for pixel in bin.pixels:
            x_coord = pixel.pix_x
            y_coord = pixel.pix_y
            #Shift because x_coord,y_coord are the center points
            rectangle = plt.Rectangle((x_coord,y_coord),1,1, color='black')
            #ax.add_patch(rectangle)
            patches.append(rectangle)
            ax.add_patch(rectangle)


    plt.xlabel("RA")
    plt.ylabel("DEC")
    plt.title("Metallicity Map for "+filename)
    norm = mpl.colors.Normalize(min(temp_norm_list),max(temp_norm_list))
    cax, _ = cbar.make_axes(ax)
    cb2 = cbar.ColorbarBase(cax, cmap=cmap, norm=norm)
    cb2.set_label(r'Abundance [$Z/Z_{\odot}$]')
    tick_list = np.linspace(min(temp_norm_list),max(temp_norm_list),num_ticks)
    ticklabel_list = np.linspace(min(temp_list),max(temp_list),num_ticks)
    ticklabel_list = [np.round(val,1) for val in ticklabel_list]
    cb2.set_ticks(tick_list)
    cb2.set_ticklabels(ticklabel_list)
    cb2.update_ticks()
    plt.savefig(file_dir+'/'+filename+"_"+str(stn)+"_Met.png")
    return ax
#-------------------------------------------------#
#-------------------------------------------------#
# Plot plot_Bins
#   parameters:
#       bin_file - WVT bin data file
#       temp_file - temperature file created during fits
#       file_dir - Full Path to location of new image
#       fileame - name of file to be printed
def plot_Bins(bin_file,temp_file,file_dir,filename,color_map,stn,wcs_image):
    #Get WCS information from header of original image
    hdu = fits.open(wcs_image)[0]
    wcs = WCS(hdu.header)
    Bins, x_min, x_max, y_min, y_max = read_in(bin_file,temp_file)
    create_image_fits(wcs_image,file_dir,x_min,x_max,y_min,y_max,Bins,filename,stn)
    fig = plt.figure()
    fig.set_size_inches(7, 7)
    ax = plt.axes(xlim=(x_min,x_max), ylim=(y_min,y_max),projection=wcs)
    #plt.axes(projection=wcs)
    N = len(Bins)
    cmap = mpl.cm.get_cmap(color_map)
    max_temp = max([bin.temp for bin in Bins])
    median_temp = np.median([bin.temp for bin in Bins])
    std_temp = np.std([bin.temp for bin in Bins])
    Bins_flush = [] #Bins that passed the sigma screening
    Bins_fail = [] #Bins that failed the sigma screening
    step_val = 1
    for bin in Bins:
        if bin.temp < 15:#median_temp +step_val*std_temp and bin.temp > median_temp - step_val*std_temp:
            Bins_flush.append(bin)
        else:
            Bins_fail.append(bin)
    max_temp = max([bin.temp for bin in Bins_flush])
    temp_list = []
    temp_norm_list = []
    bin_nums = []
    #Set up color map
    for bin in Bins_flush:
        temp_norm = bin.temp / max_temp
        temp_list.append(bin.temp)
        temp_norm_list.append(temp_norm)
    colors = cmap(temp_norm_list)
    #Create rectanges
    rect_step = 0
    for bin in Bins_flush:
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
    colors = np.linspace(min(temp_norm_list),max(temp_norm_list),N)

    #Set failed bins to black
    for bin in Bins_fail:
        patches = []
        bin_nums.append(bin.bin_number)
        #Color based on signal-to-noise value
        for pixel in bin.pixels:
            x_coord = pixel.pix_x
            y_coord = pixel.pix_y
            #Shift because x_coord,y_coord are the center points
            rectangle = plt.Rectangle((x_coord,y_coord),1,1, color='black')
            #ax.add_patch(rectangle)
            patches.append(rectangle)
            ax.add_patch(rectangle)


    plt.xlabel("RA")
    plt.ylabel("DEC")
    plt.title("Temperature Map for "+filename)
    norm = mpl.colors.Normalize(min(temp_norm_list),max(temp_norm_list))
    cax, _ = cbar.make_axes(ax)
    cb2 = cbar.ColorbarBase(cax, cmap=cmap, norm=norm)
    cb2.set_label('Temperature [KeV]')
    tick_list = np.linspace(min(temp_norm_list),max(temp_norm_list),num_ticks)
    ticklabel_list = np.linspace(min(temp_list),max(temp_list),num_ticks)
    ticklabel_list = [np.round(val,1) for val in ticklabel_list]
    cb2.set_ticks(tick_list)
    cb2.set_ticklabels(ticklabel_list)
    cb2.update_ticks()
    plt.savefig(file_dir+'/'+filename+"_"+str(stn)+".png")
    return ax

#-------------------------------------------------#
#-------------------------------------------------#
# Plot plot_Bins_ML
#   parameters:
#       bin_file - WVT bin data file
#       temp_file - temperature file created during fits
#       file_dir - Full Path to location of new image
#       fileame - name of file to be printed
def plot_Bins_ML(bin_file,temp_file,file_dir,filename,color_map,stn,wcs_image):
    #Get WCS information from header of original image
    hdu = fits.open(wcs_image)[0]
    wcs = WCS(hdu.header)
    Bins, x_min, x_max, y_min, y_max = read_in(bin_file,temp_file)
    create_image_fits(wcs_image,file_dir,x_min,x_max,y_min,y_max,Bins,filename,stn)
    fig = plt.figure()
    fig.set_size_inches(7, 7)
    ax = plt.axes(xlim=(x_min,x_max), ylim=(y_min,y_max),projection=wcs)
    #plt.axes(projection=wcs)
    N = len(Bins)
    cmap = mpl.cm.get_cmap(color_map)
    max_temp = max([bin.temp for bin in Bins])
    median_temp = np.median([bin.temp for bin in Bins])
    std_temp = np.std([bin.temp for bin in Bins])
    Bins_flush = [] #Bins that passed the sigma screening
    Bins_fail = [] #Bins that failed the sigma screening
    step_val = 1
    for bin in Bins:
        if bin.temp < 15:#median_temp +step_val*std_temp and bin.temp > median_temp - step_val*std_temp:
            Bins_flush.append(bin)
        else:
            Bins_fail.append(bin)
    max_temp = max([bin.temp for bin in Bins_flush])
    temp_list = []
    temp_norm_list = []
    bin_nums = []
    #Set up color map
    for bin in Bins_flush:
        temp_norm = bin.temp / max_temp
        temp_list.append(bin.temp)
        temp_norm_list.append(temp_norm)
    colors = cmap(temp_norm_list)
    #Create rectanges
    rect_step = 0
    for bin in Bins_flush:
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
    colors = np.linspace(min(temp_norm_list),max(temp_norm_list),N)

    #Set failed bins to black
    for bin in Bins_fail:
        patches = []
        bin_nums.append(bin.bin_number)
        #Color based on signal-to-noise value
        for pixel in bin.pixels:
            x_coord = pixel.pix_x
            y_coord = pixel.pix_y
            #Shift because x_coord,y_coord are the center points
            rectangle = plt.Rectangle((x_coord,y_coord),1,1, color='black')
            #ax.add_patch(rectangle)
            patches.append(rectangle)
            ax.add_patch(rectangle)


    plt.xlabel("RA")
    plt.ylabel("DEC")
    plt.title("Temperature Map for "+filename)
    norm = mpl.colors.Normalize(min(temp_norm_list),max(temp_norm_list))
    cax, _ = cbar.make_axes(ax)
    cb2 = cbar.ColorbarBase(cax, cmap=cmap, norm=norm)
    cb2.set_label('Temperature [KeV]')
    tick_list = np.linspace(min(temp_norm_list),max(temp_norm_list),num_ticks)
    ticklabel_list = np.linspace(min(temp_list),max(temp_list),num_ticks)
    ticklabel_list = [np.round(val,1) for val in ticklabel_list]
    cb2.set_ticks(tick_list)
    cb2.set_ticklabels(ticklabel_list)
    cb2.update_ticks()
    plt.savefig(file_dir+'/'+filename+"_"+str(stn)+"_ML.png")
    return ax

#-------------------------------------------------#

#-------------------------------------------------#
def read_in(bin_data,temp_data):
    bin_d = open(bin_data); next(bin_d); next(bin_d) #first two lines are header info
    temp_d = open(temp_data); next(temp_d)
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
    for line in temp_d:
        bins[int(int(line.split(" ")[0]))].add_temp(float(line.split(" ")[1]))
        try:
            bins[int(int(line.split(" ")[0]))].add_abund(float(line.split(" ")[4]))
            bins[int(int(line.split(" ")[0]))].add_stat(float(line.split(" ")[-1]))
        except Exception:
            pass
    temp_d.close()
    bin_d.close()
    min_x = np.min([pixel.pix_x for pixel in pixels])
    max_x = np.max([pixel.pix_x for pixel in pixels])
    min_y = np.min([pixel.pix_y for pixel in pixels])
    max_y = np.max([pixel.pix_y for pixel in pixels])
    return bins, min_x, max_x, min_y, max_y
#-------------------------------------------------#
#-------------------------------------------------#
def create_image_fits(fits_img,outroot,min_x,max_x,min_y,max_y,bins,filename,stn):
    # Create image array
    x_len = int(max_x-min_x)
    y_len = int(max_y-min_y)
    temp_array = np.zeros((x_len,y_len))
    abund_array = np.zeros((x_len,y_len))
    for bin in bins:
        for pixel in bin.pixels:
            temp_array[int(pixel.pix_x-1),int(pixel.pix_y-1)] = bin.temp
            abund_array[int(pixel.pix_x-1),int(pixel.pix_y-1)] = bin.abund
    # Copy header
    fits_ = fits.open(fits_img)
    hdr = header=fits_[0].header
    # Change image
    hdu = fits.PrimaryHDU(temp_array)
    hdul = fits.HDUList([hdu])
    fits.writeto(outroot+'/temp_'+filename+'_'+str(stn)+'.fits', temp_array.T, hdr, overwrite=True)
    hdu = fits.PrimaryHDU(abund_array)
    fits.writeto(outroot+'/abund_'+filename+'_'+str(stn)+'.fits', abund_array.T, hdr, overwrite=True)
