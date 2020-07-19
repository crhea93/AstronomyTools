#----------------------------WVT-----------------------------#
image_fits = source.img
exposure_map = none
stn_target = 50
pixel_radius = 0.5
tol = 1e-4
roundness_crit = 0.3
WVT_data = WVT_data_source_stn50
#--------------------------FITTING---------------------------#
#----------INPUT DATA------------#
base_dir = /export/carterrhea/Documents/Myriam/MACSJ1447
Name = 
ObsIDs = 10481,17233,18825
source_file = source
output_dir = binned_source_50/
Temp_data = Temp_bin_source_stn50.txt
#----------FIT INFO--------------#
redshift = 0.3755
n_H = 2.27e-2
Temp_Guess = 3.0
multi = False
#----------CHOICES---------------#
wvt = False
bin_spec = False
num_bins = 58
fit_spec = True
plot = True
Colormap = RdYlGn

