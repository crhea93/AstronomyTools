#--------------------------SPECEXTRACT-----------------------------#
reg_files = el_0,el_1,el_2,el_3,el_4,el_5,el_6,el_7
#--------------------------FITTING---------------------------#
#----------INPUT DATA------------#
base_dir = /mnt/1895813a-f52b-4ccc-9bab-1ee15fee024b/carterrhea/Pipeline-Clusters/Data/MACSJ1447
Name = 
ObsIDs = 10481,17233,18825
WVT_data = WVT_data.txt
source_file = el_
output_dir = binned_el/
Temp_data = Temp_bin.txt
multi = false
#----------FIT INFO--------------#
redshift = 0.3755
n_H = 2.27e-2
Temp_Guess = 6.0
#----------CHOICES---------------#
extract_spectrum = False
fit_spec = True
plot = True


