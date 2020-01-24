#--------------------------SPECEXTRACT-----------------------------#
reg_files = North_filament_0,North_filament_1,North_filament_2
#--------------------------FITTING---------------------------#
#----------INPUT DATA------------#
base_dir = /home/carterrhea/Documents/Perseus/Xray
Name = 
ObsIDs = 11713,11714,11715,11716,12025,12033,12036,12037
WVT_data = WVT_data.txt
source_file = North_filament
output_dir = binned_NFs/
Temp_data = Temp_bin_NFs.txt
multi = false
#----------FIT INFO--------------#
redshift = 0.0179
n_H = 0.137
Temp_Guess = 2.0
#----------CHOICES---------------#
extract_spectrum = False
fit_spec = True
plot = True


