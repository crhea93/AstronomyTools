#----------------------------WVT-----------------------------#
image_fits = source.img
exposure_map = none
stn_target = 50
pixel_radius = 0.5
tol = 1e-4
roundness_crit = 0.3
#--------------------------FITTING---------------------------#
#----------INPUT DATA------------#
base_dir = /home/carterrhea/Desktop/PipelineClusters/Data/MACSJ1447
Name = Merge_unbinned
ObsIDs = 10481,17233,18825
WVT_data = WVT_data.txt
source_file = source
output_dir = binned/
#----------FIT INFO--------------#
redshift = 0.3755
n_H = 2.27e-2
Temp_Guess = 2.0
#----------CHOICES---------------#
wvt = False
bin_spec = False
num_bins = 61
fit_spec = False
plot = True
Colormap = RdYlGn

