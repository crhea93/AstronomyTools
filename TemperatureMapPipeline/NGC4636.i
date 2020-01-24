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
base_dir = /mnt/1895813a-f52b-4ccc-9bab-1ee15fee024b/carterrhea/Pumpkin/Test
Name = NGC4636
ObsIDs = 323,324,3926,4415
source_file = center 
output_dir = binned_50/
Temp_data = Temp_bin_center_stn50.txt
multi = False
#----------FIT INFO--------------#
redshift = 0.003129
n_H = 1.84e-2
Temp_Guess = 2.0
#----------CHOICES---------------#
wvt = False
bin_spec = False
num_bins = 150
fit_spec = True
plot = True
Colormap = viridis


