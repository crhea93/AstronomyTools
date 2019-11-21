--------------------------SPECEXTRACT-----------------------------#
reg_files = 600kpc
#--------------------------FITTING---------------------------#
#----------INPUT DATA------------#
base_dir = /mnt/1895813a-f52b-4ccc-9bab-1ee15fee024b/carterrhea/Pipeline-Clusters/Data/A907
Name = A907
ObsIDs = 3185,3205
source_file = 600kpc
output_dir = binned_600/
Temp_data = Temp_bin_600.txt
multi = false
#----------FIT INFO--------------#
redshift = 0.153
n_H = 0.0528
Temp_Guess = 3.0
#----------CHOICES---------------#
extract_spectrum = True
fit_spec = True
plot = True
