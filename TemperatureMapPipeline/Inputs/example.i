#----------------------------WVT-----------------------------#
image_fits = source.img
exposure_map = none
stn_target = 150
pixel_radius = 0.5
tol = 1e-4
roundness_crit = 0.3
WVT_data = WVT_data_source_stn150
#--------------------------FITTING---------------------------#
#----------INPUT DATA------------#
base_dir = /home/user/Documents/Cluster_Data/Perseus
Name =
ObsIDs = 3209,4289
source_file = source
output_dir = binned_source_150/
Temp_data = Temp_bin_source_stn150.txt
#----------FIT INFO--------------#
redshift = 0.0179
n_H = 0.137
Temp_Guess = 2.0
#----------CHOICES---------------#
wvt = True
bin_spec = True
# Only set num_bins if bin_spec = False
num_bins = 0
fit_spec = True
plot = True
Colormap = magma
