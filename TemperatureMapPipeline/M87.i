#----------------------------WVT-----------------------------#
image_fits = reg1_ext.img
exposure_map = none
stn_target = 80
pixel_radius = 0.5
tol = 1e-4
roundness_crit = 0.3
WVT_data = WVT_data_reg1_ext_stn80
#--------------------------FITTING---------------------------#
#----------INPUT DATA------------#
base_dir = /mnt/1895813a-f52b-4ccc-9bab-1ee15fee024b/carterrhea/M87/X-ray
Name = Data
ObsIDs = 2707,3717,5826,5827,5828,6186,7210,7211,7212
#ObsIDs = 5826,5827
source_file = reg1_ext
output_dir = binned_reg1_ext_80/
Temp_data = Temp_bin_reg1_ext_stn80.txt
multi = False
#----------FIT INFO--------------#
redshift = 0.00428
n_H = 0.0194
Temp_Guess = 2.0
#----------CHOICES---------------#
wvt = False
bin_spec = False
num_bins = 108
fit_spec = True
plot = True
Colormap = inferno

