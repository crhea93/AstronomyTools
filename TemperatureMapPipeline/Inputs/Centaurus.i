#----------------------------WVT-----------------------------#
image_fits = source.img
exposure_map = none
stn_target = 60
pixel_radius = 0.5
tol = 1e-4
roundness_crit = 0.3
#--------------------------FITTING---------------------------#
#----------INPUT DATA------------#
base_dir = /home/carterrhea/Desktop/PipelineClusters/Data/Pheonix
Name = Pheonix
ObsIDs = 16135,16545
WVT_data = WVT_data.txt
source_file = source
output_dir = binned/
#----------FIT INFO--------------#
redshift = 0.0104
n_H = 0.12
Temp_Guess = 3.0
#----------CHOICES---------------#
wvt = True
bin_spec = True
num_bins = 0
fit_spec = True
plot = True
Colormap = inferno

