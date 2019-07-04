#----------------------------WVT-----------------------------#
image_fits = source.img
exposure_map = none
stn_target = 80
pixel_radius = 0.5
tol = 1e-4
roundness_crit = 0.3
#--------------------------FITTING---------------------------#
#----------INPUT DATA------------#
base_dir = /home/carterrhea/Desktop/PipelineClusters/Data/M87
Name = M87_2
#ObsIDs = 2707,3717,5826,5827,5828,6186,7210,7211,7212
ObsIDs = 5826,5827
WVT_data = WVT_data.txt
source_file = source
output_dir = binned/
#----------FIT INFO--------------#
redshift = 0.00428
n_H = 0.0194
Temp_Guess = 2.0
#----------CHOICES---------------#
wvt = False
bin_spec = False
num_bins = 0
fit_spec = False
plot = True
Colormap = inferno

