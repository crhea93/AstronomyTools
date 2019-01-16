'''
This program will plot the final temperature map using the temperatures calculated from sherpa and the binnings from WVT
'''
import os
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colorbar as cbar
from matplotlib.collections import PatchCollection
#-----------------------------------INPUTS------------------------------------------#
base_dir = '/home/user/Documents/NGC4636/'
bin_file = 'Merged_unbinned/WVT_data.txt'
temp_file = 'Temp_bin.txt'
output_dir = base_dir
output_name = 'NGC4636_temperature'
color_map = 'cool'
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
    def add_pixel(self,pixel):
        self.pixels.append(pixel)
    def add_temp(self,temp):
        self.temp = temp
    def add_stat(self,stat):
        self.stat = stat

class Pixel:
    def __init__(self, number, pix_x, pix_y):
        self.pix_number = number
        self.pix_x = pix_x
        self.pix_y = pix_y
#-------------------------------------------------#
#-------------------------------------------------#
# Plot plot_Bins
#   parameters:
#       Bin - list of bin objects
#       x_min - minimum x value of pixel
#       x_max - maximum x value of pixel
#       y_min - minimum y value of pixel
#       x_max - maximum y value of pixel
#       file_dir - Full Path to location of new image
#       fileame - name of file to be printed
def plot_Bins(Bins,x_min,x_max,y_min,y_max,file_dir,filename,color_map):
    fig = plt.figure()
    fig.set_size_inches(7, 7)
    ax = plt.axes(xlim=(x_min,x_max), ylim=(y_min,y_max))
    N = len(Bins)
    cmap = mpl.cm.get_cmap(color_map)
    max_temp = max([bin.temp for bin in Bins])
    median_temp = np.median([bin.temp for bin in Bins])
    temp_list = []
    temp_norm_list = []
    bin_nums = []
    #Set up color map
    for bin in Bins:
        temp_norm = bin.temp / max_temp
        temp_list.append(bin.temp)
        temp_norm_list.append(temp_norm)

    colors = cmap(temp_norm_list)
    #colors.set_clim(min(temp_norm_list),max(temp_norm_list))
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
    colors = np.linspace(min(temp_norm_list),max(temp_norm_list),N)

    #p = PatchCollection(patches, cmap=cmap)
    #p.set_array(np.array(colors))
    #ax.add_collection(p)
    #cbar = fig.colorbar(p, ax=ax)

    plt.xlabel("X")
    plt.ylabel("Y")
    plt.title("Temperature Map for NGC 4636")
    norm = mpl.colors.Normalize(min(temp_norm_list),max(temp_norm_list))
    cax, _ = cbar.make_axes(ax)
    cb2 = cbar.ColorbarBase(cax, cmap=cmap, norm=norm)
    cb2.set_label('Temperature (KeV)')
    tick_list = np.linspace(min(temp_norm_list),max(temp_norm_list),num_ticks)
    ticklabel_list = np.linspace(min(temp_list),max(temp_list),num_ticks)
    ticklabel_list = [np.round(val,1) for val in ticklabel_list]
    cb2.set_ticks(tick_list)
    cb2.set_ticklabels(ticklabel_list)
    cb2.update_ticks()

    plt.savefig(file_dir+'/'+filename+".png")
    return None

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
        bins[int(int(line.split(" ")[0]))].add_stat(float(line.split(" ")[2]))
    temp_d.close()
    bin_d.close()
    min_x = np.min([pixel.pix_x for pixel in pixels])
    max_x = np.max([pixel.pix_x for pixel in pixels])
    min_y = np.min([pixel.pix_y for pixel in pixels])
    max_y = np.max([pixel.pix_y for pixel in pixels])
    return bins, min_x, max_x, min_y, max_y

#-------------------------------------------------#
#-------------------------------------------------#
#-------------------------------------------------#
#-------------------------------------------------#
def main():
    os.chdir(base_dir)
    bins, min_x, max_x, min_y, max_y = read_in(bin_file,temp_file)
    plot_Bins(bins, min_x, max_x, min_y, max_y,output_dir,output_name,color_map)
    return None

main()
