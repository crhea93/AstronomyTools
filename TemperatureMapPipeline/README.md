For the mechanics of what is going on, please see the the documentation :)

This program will results in:
- A Weighted Voronoi Tessellation 
- Spectra for each bin
- Temperature and Abundance Maps

To run the program please supply the following items in an input file (see example)

#----------------------------WVT-----------------------------#
image_fits = source.img #Name of source fits image for WVT algorithm
exposure_map = none #If there is an exposure map, please include it here
stn_target = 75 #Signal-to-Noise target
pixel_radius = 0.5 #Pixel radius for WVT (don't change for Chandra)
tol = 1e-4 #Tolerance of WVT fitting (Shouldn't need to change)
roundness_crit = 0.3 #Criteria of WVT bin (Standard value)
#--------------------------FITTING---------------------------#
#----------INPUT DATA------------#
base_dir = /home/user/Documents/M87_Xray/M87 #Path to obsid directories
Name = M87_2 #Name of output directory
ObsIDs = 5826,5827 #List of observation ids
WVT_data = WVT_data.txt #WVT output file (Default is WVT_data.txt)
source_file = source #Name of source fits file
output_dir = binned/ #Output directory name 
#----------FIT INFO--------------#
redshift = 0.00428 #Redshift of cluster
n_H = 0.0194 #Column density
Temp_Guess = 3.0 #Guess temperature (default is 3.0)
#----------CHOICES---------------#
wvt = False #Would you like to construct a WVT map? (Default True)
bin_spec = False #Do you need to create the spectra for each bin? (Default True)
num_bins = 0 #Only supply if bin_spec is False
fit_spec = False #Do you need to fit the spectra? (Default True)
plot = True #Would you like to plot the maps? (Default True)
Colormap = inferno #Matplotlib colormap
