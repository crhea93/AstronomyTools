#----------------------------WVT-----------------------------#
image_fits = center.img
exposure_map = none
stn_target = 50
pixel_radius = 0.5
tol = 1e-4
roundness_crit = 0.3
WVT_data = WVT_data_center_stn50
#--------------------------FITTING---------------------------#
#----------INPUT DATA------------#
base_dir = /home/carterrhea/Desktop/Tests
Name = NGC4636
ObsIDs = 323,324
source_file = center 
output_dir = binned/
Temp_data = Temp_bin_center_stn50.txt
multi = False
#----------FIT INFO--------------#
redshift = 0.003129
n_H = 1.91e-2
Temp_guess = 1.0
#----------CHOICES---------------#
wvt = True
bin_spec = True
num_bins = 0
fit_spec = True
plot = True
Colormap = RdYlGn


