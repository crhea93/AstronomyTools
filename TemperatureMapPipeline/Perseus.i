#----------------------------WVT-----------------------------#
image_fits = PC.img
exposure_map = none
stn_target = 500
pixel_radius = 0.5
tol = 1e-4
roundness_crit = 0.3
WVT_data = WVT_data_PC_stn200
#--------------------------FITTING---------------------------#
#----------INPUT DATA------------#
base_dir = /mnt/1895813a-f52b-4ccc-9bab-1ee15fee024b/carterrhea/Perseus/Xray
Name = 
ObsIDs = 11713,12025,12033,12036,11715,11716,12037,11714
source_file = PC
output_dir = binned_source_PC_200/
Temp_data = Temp_bin_PC_stn200.txt
multi = False
#----------FIT INFO--------------#
redshift = 0.0179
n_H = 0.137
Temp_Guess = 2.0
#----------CHOICES---------------#
wvt = True
bin_spec = False
num_bins = 132
fit_spec = False
plot = False
Colormap = inferno

