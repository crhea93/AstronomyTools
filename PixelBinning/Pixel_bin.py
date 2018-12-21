'''
DEFINE INPUTS IN MAIN FUNCTION AT END OF FILE
This file will create a bin map for each pixel. In reality we are simply creating a region about each pixel
with a predetermined signal-to-noise value.
---------------------------------------------------
---------------------------------------------------
Inputs:
    image_fits - Name of Image Fits File (e.g. "simple_image.fits")
    exposure_map - Name of exposure map -- optional -- (e.g. 'flux_broad.expmap')
    StN_Target - Target Signal-to-Noise (e.g. 50)
    pixel_size - Pixel radius in degrees (e.g. 0.492 for CHANDRA AXIS I)
    ToL - Tolerance Level (e.g. 1e-16)
    roundness_crit - Critical Value Describing Roundness of Bin (e.g. 0.3)
    home_dir - Full Path to Home Directory (e.g. '/home/user/Documents/ChandraData/12833/repro')
    image_dir - Full Path to Image Directory (e.g. '/home/user/Desktop')
    output_dir - Full Path to Output Directory (e.g. '/home/user/Documents/ChandraData/12833')
---------------------------------------------------
List of Functions (In order of appearance):
    plot_Bins --> Plot bins with five alternating color schemes
    plot_Bins_SN --> Plot bins with color based of Signal-to-Noise value
    read_in --> Read in flux file and gather flux information
    Nearest_Neighbors --> Calculate nearest neighbors using builtin sklearn algorithm
    dist --> Calculates Euclidean distances
    closest_node --> Calculates closest node to bin centroid
    adjacency --> Determines if pixels are adjacent (diagonals don't count)
    Roundness --> Calculates roundness parameter of bin
    Potential_SN --> Calculates Potential New Signal-to-Noise if pixel is added
    Bin_data --> Creates data for Chandra with only unique pixels from read_in
    bin_num_pixel --> Calculates bixel's associated bin (i.e. in which bin is my pixel??)
    assigned_missing_pixels --> Assign pixels which weren't read in due to 0 Signal-to-Noise
    Bin_Acc --> Perform Bin Accretion Algorithm
    converged_met --> Check that all bins have converged in WVT
    Rebin_Pixels --> Rebin pixels in WVT algorithm
    WVT --> Perform Weighted Voronoi Tessellation Algorithm
---------------------------------------------------
Outputs:
    -- Creates bin plots for bin accretion and WVT algorithms
    -- Creates text files for each bin containing chip information and bin number
                    pixel (x,y) are CENTER coordinates!
---------------------------------------------------
TO DO:
    -- Find a better way to get the x and y pixel coordinates
        in bounding box. Current method SOMETIME (RARELY) incure
        errors. The errors aren't a big deal and don't actually change
        the results, but an error is an error!!
---------------------------------------------------
---------------------------------------------------
Carter Rhea
https://carterrhea.com
carterrhea93@gmail.com
'''

#-----------------INPUTS--------------------------#
import os
import sys
import numpy as np
import statistics as stats
from astropy.io import fits
from sklearn.neighbors import NearestNeighbors
from matplotlib.collections import PatchCollection
import matplotlib.pyplot as plt
import matplotlib as mpl
#-------------------------------------------------#
#-------------------------------------------------#
#-------------------------------------------------#
#-------------------------------------------------#
# Plot Bins
#   parameters:
#       Bin - list of bin objects
#       x_min - minimum x value of pixel
#       x_max - maximum x value of pixel
#       y_min - minimum y value of pixel
#       x_max - maximum y value of pixel
#       file_dir - Full Path to location of new image
#       fileame - name of file to be printed
def plot_Bins(Bins,x_min,x_max,y_min,y_max,StN_Target,file_dir,filename):
    fig = plt.figure()
    #fig.set_size_inches(7, 6.5)
    ax = plt.axes(xlim=(x_min,x_max), ylim=(y_min,y_max))
    N = len(Bins)
    StN_list = []
    SNR_list = []
    bin_nums = []
    max_StN = max([bin.StN[0] for bin in Bins])
    StN_list = [bin.StN[0] for bin in Bins]
    median_StN = np.median(StN_list)
    stand_dev = stats.stdev(StN_list)
    mini_pallete = ['mediumspringgreen','salmon','cyan','orchid','yellow','blue','red','magenta','black','white']
    binNumber  = 0
    for bin in Bins:
        bin_nums.append(bin.bin_number)
        SNR = bin.StN[0]/median_StN
        SNR_list.append(SNR)
        for pixel in bin.pixels:
            x_coord = pixel.pix_x
            y_coord = pixel.pix_y
            #patches.append(Rectangle((x_coord,y_coord),1,1))
            if binNumber%10 == 0:
                color = mini_pallete[0]
            if binNumber%10 == 1:
                color = mini_pallete[1]
            if binNumber%10 == 2:
                color = mini_pallete[2]
            if binNumber%10 == 3:
                color = mini_pallete[3]
            if binNumber%10 == 4:
                color = mini_pallete[4]
            if binNumber%10 == 5:
                color = mini_pallete[5]
            if binNumber%10 == 6:
                color = mini_pallete[6]
            if binNumber%10 == 7:
                color = mini_pallete[7]
            if binNumber%10 == 8:
                color = mini_pallete[8]
            if binNumber%10 == 9:
                color = mini_pallete[9]
            #Shift because x_coord,y_coord are the center points
            rectangle = plt.Rectangle((x_coord,y_coord),1,1, fc=color)
            ax.add_patch(rectangle)
        binNumber += 1
    centroids_x = [Bins[i].centroidx[0] for i in range(len(Bins))]
    centroids_y = [Bins[i].centroidy[0] for i in range(len(Bins))]
    ax.scatter(centroids_x,centroids_y,marker='+',c="black")
    plt.xlabel("X")
    plt.ylabel("Y")
    plt.title("Bin Mosaic")
    plt.savefig(file_dir+'/'+filename+".png")
    plt.clf()
    plt.hist(StN_list)
    plt.xlim((median_StN-3*stand_dev,median_StN+3*stand_dev))
    plt.ylabel("Number of Bins")
    plt.xlabel("Signal-to-Noise")
    plt.title("Signal-to-Noise per Bin")
    plt.savefig(file_dir+'/histograms/'+filename+".png")
    plt.clf()
    SNR_std = stats.stdev(SNR_list)
    SNR_med = np.median(SNR_list)
    SNR_nel = len(SNR_list)
    n_el = SNR_nel #Just change this
    plt.scatter(np.arange(n_el),SNR_list, marker = '+', color='salmon',label="Data Points")
    plt.plot(np.arange(n_el),[SNR_med for i in range(n_el)], linestyle='solid', color= 'forestgreen', label="Median")
    plt.plot(np.arange(n_el),[SNR_med+SNR_std for i in range(n_el)], linestyle='--', color= 'black')
    plt.plot(np.arange(n_el),[SNR_med-SNR_std for i in range(n_el)], linestyle='--', color= 'black', label="sigma")
    plt.title("Signal to Noise Ratio")
    plt.ylim((min(SNR_list),max(SNR_list)))
    plt.xlabel("Bin Number")
    plt.ylabel("Signal-to-Noise Normalized by Median Value")
    plt.legend(loc='center left',bbox_to_anchor=(1.0, 0.5),
          ncol=1, fancybox=True, shadow=True)
    plt.savefig(file_dir+'/'+filename+"_scatter.png", bbox_inches="tight")
    plt.clf()
    return None
#-------------------------------------------------#
#-------------------------------------------------#
# Plot plot_Bins_SN
#   parameters:
#       Bin - list of bin objects
#       x_min - minimum x value of pixel
#       x_max - maximum x value of pixel
#       y_min - minimum y value of pixel
#       x_max - maximum y value of pixel
#       file_dir - Full Path to location of new image
#       fileame - name of file to be printed
def plot_Bins_SN(Bins,x_min,x_max,y_min,y_max,file_dir,filename):
    fig = plt.figure()
    fig.set_size_inches(7, 7)
    ax = plt.axes(xlim=(x_min,x_max), ylim=(y_min,y_max))
    N = len(Bins)
    cmap = mpl.cm.get_cmap('viridis')
    max_StN = max([bin.StN[0] for bin in Bins])
    median_StN = np.median([bin.StN[0] for bin in Bins])
    patches = []
    SN_list = []
    SNR_list = []
    bin_nums = []
    for bin in Bins:
        bin_nums.append(bin.bin_number)
        #Color based on signal-to-noise value
        SNR = bin.StN[0]/median_StN
        SN_list.append(bin.StN[0])
        SNR_list.append(SNR)
        for pixel in bin.pixels:
            x_coord = pixel.pix_x
            y_coord = pixel.pix_y
            #Shift because x_coord,y_coord are the center points
            rectangle = plt.Rectangle((x_coord,y_coord),1,1, fc=cmap(SNR))
            #ax.add_patch(rectangle)
            patches.append(rectangle)
    colors = np.linspace(min(SNR_list),max(SNR_list),N)
    p = PatchCollection(patches, cmap=cmap)
    p.set_array(np.array(colors))
    ax.add_collection(p)
    cbar = fig.colorbar(p, ax=ax)
    plt.xlabel("X")
    plt.ylabel("Y")
    plt.title("Signal to Noise per Bin")
    cbar.set_label('Normalized Signal to Noise', rotation=270, labelpad=20)
    plt.savefig(file_dir+'/'+filename+".png")
    plt.clf()
    return None
#-------------------------------------------------#
#-------------------------------------------------#
# Bin Class Information
# Everything in here should be self-explanatory... if not let me know and I
# will most happily comment it! :)
class Bin:
    def __init__(self, number):
        self.bin_number = number
        self.pixels = []
        self.pixel_neighbors = []
        self.centroidx = [0] #So I can pass by value
        self.centroidy = [0]
        self.centroidx_prev = [0]
        self.centroidy_prev = [0]
        self.StN = [0]
        self.StN_prev = [0]
        self.Signal = [0]
        self.Noise = [0]
        self.Area = [0]
        self.Area_prev = [0]
        self.scale_length = [0]
        self.scale_length_prev = [0]
        self.successful = False
        self.WVT_successful = False
        self.avail_reassign = True #Can a pixel be reassigned to you?

    def add_pixel(self, Pixel):
        self.StN[0] = 0
        self.pixels.append(Pixel)
        self.Signal[0] += Pixel.Signal
        self.Noise[0] += Pixel.Noise
        if self.Noise[0] != 0:
            self.StN[0] = self.Signal[0]/(np.sqrt(self.Noise[0]))
        for neigh in Pixel.neighbors:
            if (neigh not in self.pixel_neighbors):
                self.pixel_neighbors.append(neigh)

    def clear_pixels(self):
        self.pixels = []
        self.StN[0] = 0
        self.Signal[0] = 0
        self.Noise[0] = 0
        self.successful = False
        self.WVT_successful = False
        self.avail_reassign = True

    def remove_pixel(self, Pixel):
        self.pixels.remove(Pixel)
        self.StN[0] = 0
        self.Signal[0] -= Pixel.Signal
        self.Noise[0] -= Pixel.Noise
        if self.Noise[0] != 0:
            self.StN[0] = self.Signal[0]/(np.sqrt(self.Noise[0]))

    def success(self):
        self.successful = True

    def WVT_success(self):
        self.WVT_successful = True

    def availabilty(self):
        self.avail_reassign = False

    def update_StN_prev(self):
        self.StN_prev[0] = self.StN[0]

    def CalcCentroid(self):
        self.centroidx_prev[0] = self.centroidx[0]
        self.centroidy_prev[0] = self.centroidy[0]
        self.centroidx[0] = 0
        self.centroidy[0] = 0
        n_cent = 0
        for pixel in self.pixels:
            self.centroidx[0] += pixel.pix_x
            self.centroidy[0] += pixel.pix_y
            n_cent += 1
        self.centroidx[0] *= (1/n_cent)
        self.centroidy[0] *= (1/n_cent)

    def CalcArea(self,pixel_length):
        self.Area_prev[0] = self.Area[0]
        self.Area[0] = 0
        self.Area[0] += pixel_length**2*len(self.pixels)

    def CalcScaleLength(self,StN_Target):
        self.scale_length_prev[0] = self.scale_length[0]
        self.scale_length[0] = 0
        self.scale_length[0] += np.sqrt((self.Area[0]/np.pi)*(StN_Target/self.StN[0]))
#-------------------------------------------------#
#-------------------------------------------------#
# Pixel Class Information
# Ditto à propos the documentation for this class
class Pixel:
    def __init__(self, number, pix_x, pix_y, signal, noise):
        self.pix_number = number
        self.pix_x = pix_x
        self.pix_y = pix_y
        self.Signal = signal
        self.Noise = noise
        self.StN = 0
        if self.Noise != 0:
            self.StN = self.Signal/np.sqrt(self.Noise)
        self.neighbors = []
        self.neighbors_x = []
        self.neighbors_y = []
        self.assigned_to_bin = False
        self.assigned_bin = None
    def add_neighbor(self, pixel,x_pixel,y_pixel):
        self.neighbors.append(pixel)
        self.neighbors_x.append(x_pixel)
        self.neighbors_y.append(y_pixel)
    def add_to_bin(self,bin):
        self.assigned_to_bin = True
        self.assigned_bin = bin
    def clear_bin(self):
        self.assigned_to_bin = False
        self.assigned_bin = None
#-------------------------------------------------#
#-------------------------------------------------#
# READ IN
# we first must read in our data from the Chandra file
# pixel (x,y) are CENTER coordinates!
#   parameters:
#       fits_file = fits file in string format
#       image_fits = fits image file in string format
#       image_coord = reg file containing center of box of interest and sizes
#       exposure_map = exposure map file in string format
def read_in(image_fits,exposure_map = None):
    #Collect Pixel Data
    hdu_list = fits.open(image_fits, memmap=True)
    exposure_time = float(hdu_list[0].header["TSTOP"]) - float(hdu_list[0].header["TSTART"])
    counts = hdu_list[0].data
    y_len = counts.shape[0]
    x_len = counts.shape[1]
    hdu_list.close()
    x_min = 0; y_min = 0;
    x_max = x_len; y_max = y_len;
    if exposure_map != None:
        hdu_list = fits.open(exposure_map, memmap=True)
        exposure = hdu_list[0].data
        hdu_list.close()
    #Currently Do not bother reading background information
    '''bkg_hdu = fits.open(bkg_image_fits, memmap=True)
    bkg_counts = bkg_hdu[0].data
    bkg_hdu.close()
    avg_bkg_counts = np.mean(bkg_counts)
    bkg_sigma = np.std(bkg_counts)'''
    Pixels = []
    pixel_count = 0
    for col in range(int(x_len)):
        for row in range(int(y_len)):
            if exposure_map == None:
                flux = counts[row][col]
                vari = counts[row][col]
            else:
                flux = counts[row][col]/(exposure[row][col]*exposure_time) #- avg_bkg_counts/(exposure[row][col]*exposure_time)
                vari = counts[row][col]/(exposure[row][col]**2*exposure_time**2) #+ bkg_sigma
            Pixels.append(Pixel(pixel_count,x_min+col,y_min+row,flux,vari)) #Bottom Left Corner!
            pixel_count += 1
    #print("We have "+str(pixel_count)+" Pixels! :)")
    return Pixels, x_min, x_max, y_min, y_max
#-------------------------------------------------#
#-------------------------------------------------#
# CALCULATE NearestNeighbors -- REALLY AN ADJACENCY LIST
#   http://scikit-learn.org/stable/modules/neighbors.html
#   parameters:
#       pixel_list - list of pixel objects
def Nearest_Neighbors(pixel_list):
    print("Running Nearest Neighbor Algorithm")
    xvals = []
    yvals = []
    num_neigh = 100
    for pixel in pixel_list:
        xvals.append(pixel.pix_x)
        yvals.append(pixel.pix_y)
    X = np.column_stack((xvals,yvals))
    nbrs = NearestNeighbors(n_neighbors=num_neigh, algorithm='ball_tree').fit(X)
    distances, indices = nbrs.kneighbors(X)
    pix_num = 0
    for pixel in pixel_list:
        for j in range(num_neigh-1):
            if distances[pix_num][j+1] == 1:
                index = indices[pix_num][j+1]
                pixel.add_neighbor(pixel_list[index],xvals[index],yvals[index])
            else:
                pass #not adjacent
        pix_num += 1
    print("Finished Nearest Neighbor Algorithm")
    return None
#-------------------------------------------------#
# CALCULATE NearestNeighbors -- REALLY AN ADJACENCY LIST
#   http://scikit-learn.org/stable/modules/neighbors.html
#   parameters:
#       pixel_list - list of pixel objects
def Nearest_Neighbors_ind(pixel_list,num_neigh,pix_id):
    print("Running Nearest Neighbor Algorithm")
    xvals = []
    yvals = []
    for pixel in pixel_list:
        xvals.append(pixel.pix_x)
        yvals.append(pixel.pix_y)
    X = np.column_stack((xvals,yvals))
    nbrs = NearestNeighbors(n_neighbors=num_neigh, algorithm='ball_tree').fit(X)
    distances, indices = nbrs.kneighbors(X)
    new_neighbors = []
    for j in range(num_neigh-1):
        if distances[pix_id][j+1] == 1:
            index = indices[pix_id][j+1]
            new_neighbors.append(pixel_list[index],xvals[index],yvals[index])
        else:
            pass #not adjacent
    print("Finished Nearest Neighbor Algorithm")
    return new_neighbors


#-------------------------------------------------#
#-------------------------------------------------#
# CALCULATE StN_Potential
#   parameters:
#       Current_bin - current bin object
#       closest_pix - closest pixel to bin
def Potential_SN(Current_bin,closest_pix):
    Current_bin.add_pixel(closest_pix)
    new_StN = Current_bin.StN[0]
    Current_bin.remove_pixel(closest_pix)
    return new_StN
#-------------------------------------------------#
#-------------------------------------------------#
# Create bin information for chandra
#   prameters:
#       Bins - list of final bins
#       min_x - minimum pixel x value
#       min_y - minimum pixel y value
#       output_directory - directory where txt file will be located
#       filename - name of text file
def Bin_data(Bins,missing_pixels,min_x,min_y, output_directory, filename):
    Bins.sort(key=lambda bin: bin.bin_number)
    file = open(output_directory+'/'+filename+'.txt',"w+")
    file.write("This text file contains information necessary for chandra to bin the pixels appropriately for image.fits \n")
    file.write("pixel_x pixel_y bin \n")
    file2 = open(output_directory+'/'+filename+'_bins.txt','w+')
    file2.write("Bin data for paraview script to plot Weighted Voronoi Diagram \n")
    file2.write("centroidx centroidy weight \n")
    binCount = 0
    for bin in Bins:
        for pixel in bin.pixels:
            file.write(str(pixel.pix_x-min_x)+" "+str(pixel.pix_y-min_y)+" "+str(binCount)+' \n')
        file2.write(str(bin.centroidx[0])+" "+str(bin.centroidy[0])+" "+str(bin.scale_length[0])+" \n")
        binCount += 1
    file.close()
    file2.close()
    return None
#-------------------------------------------------#
#-------------------------------------------------#
# Find bin number given pixel
#   parameters:
#       binList - list of bin objects
#       currentPixel  - current pixel object
def bin_num_pixel(binList,currentPixel):
    bin_number_interest = None
    for binNumber in range(len(binList)):
        if currentPixel in binList[binNumber].pixels:
            bin_number_interest = binNumber
            break
    return bin_number_interest
#-------------------------------------------------#
#-------------------------------------------------#
# Bin_Creation Algorithm
#   parameters:
#       Pixels - list of unique pixels
#       pixel_length - physical size of pixel
#       StN_Target - Target value of Signal-to_noise
#       roundness_crit - Critical roundness value (upper bound)
def Bin_Creation(Pixels,pixel_length,StN_Target,roundness_crit):
    #step 1:setup list of bin objects
    print("Starting Bin Accretion Algorithm")
    binCount = 0
    Bin_list = []
    StN_found = False
    for pixel_ in Pixels[:]:
        Current_bin = Bin(binCount)
        Bin_list.append(Current_bin)
        Current_bin.add_pixel(pixel_)
        #Lets add pixels to the bin until the Signal to Noise is reached
        neighbors_unsearched = pixel_.neighbors
        neigh_num = 0
        i = 1  #how many times we have calculated nearest neighbors
        while StN_found == False:
            if len(neighbors_unsearched) == 0:
                neighbors_unsearched = Nearest_Neighbors_ind(Pixels[:],10**(i+1),pixel_.pix_number)
                i += 1 #So we can search further!
            Current_bin.add_pixel(neighbors_unsearched[neigh_num])
            if Current_bin.StN > StN_Target:
                StN_found = False

            neigh_num += 1
        binCount += 1


    print("Completed Bin Accretion Algorithm")
    print("There are a total of "+str(len(Bin_list)+1)+" bins!")
    return Bin_list


#-------------------------------------------------#
#-------------------------------------------------#
#Read input file
#   parameters:
#       input file - .i input file
def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False
def read_input_file(input_file,number_of_inputs):
    inputs = {}
    with open(input_file) as f:
        for line in f:
            if '=' in line:
                inputs[line.split("=")[0].strip().lower()] = line.split("=")[1].strip()
            else: pass
        if len(inputs) != number_of_inputs:
            print("Please recheck the input file since some parameter is missing...")
            print("Exiting program...")
            exit()
        else:
            print("Successfully read in input file")
            for key,val in inputs.items():
                if is_number(val) == True:
                    inputs[key] = float(val)
        return inputs
#-------------------------------------------------#
#-------------------------------------------------#
# #-------------------------------------------------#
#-------------------------------------------------#
def main():
    inputs = read_input_file(sys.argv[1],float(sys.argv[2]))
    os.chdir(inputs['home_dir'])
    if os.path.isdir(inputs['output_dir']+'/histograms') == False:
        os.mkdir(inputs['output_dir']+'/histograms')
    else:
        for graph in os.listdir(inputs['output_dir']+'/histograms'):
            if graph.endswith(".png"):
                os.remove(os.path.join(inputs['output_dir']+'/histograms', graph))
    for key,val in inputs.items():
        print("     "+key+" = "+str(val))
    pixel_size = inputs['pixel_radius']*2
    print("#----------------Algorithm Part 1----------------#")
    Pixels,min_x,max_x,min_y,max_y = read_in(inputs['image_fits'],inputs['exposure_map'])
    Nearest_Neighbors(Pixels)
    Bins = Bin_Creation(Pixels,pixel_size,inputs['stn_target'],inputs['roundness_crit'])
    plot_Bins(Bins,min_x,max_x,min_y,max_y,inputs['stn_target'],inputs['image_dir'],"bin_acc")
    print("#----------------Algorithm Complete--------------#")
    Bin_data(Bins,Pixels,min_x,min_y,inputs['output_dir'],"WVT_data")
    print("#----------------Information Stored--------------#")
main()
