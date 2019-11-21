#----------------------------WVT-----------------------------#
image_fits = North_filament.img
exposure_map = none
stn_target = 80
pixel_radius = 0.5
tol = 1e-4
roundness_crit = 0.3
WVT_data = WVT_data_NorthFilament_stn80
#--------------------------FITTING---------------------------#
#----------INPUT DATA------------#
base_dir = /mnt/1895813a-f52b-4ccc-9bab-1ee15fee024b/carterrhea/Perseus/Xray
Name = 
ObsIDs = 11713,12025,12033,12036,11715,11716,12037,11714
source_file = North_filament
output_dir = binned_source_80/
Temp_data = Temp_bin_NorthFilament_stn80.txt
multi = False
#----------FIT INFO--------------#
redshift = 0.0179
n_H = 0.137
Temp_Guess = 2.0
#----------CHOICES---------------#
wvt = False
bin_spec = True
num_bins = 2
fit_spec = True
plot = True
Colormap = inferno

