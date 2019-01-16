'''
This python file will unzip a tar file and then continue to unzip the residual
gunzipped files

This is to be used with data recovered from the Chandra telescope

parameters:
    chandra_dir = chandra_dir directory (e.g. "/home/usrname/Documents/ChandraData"); String
    tar_file = filename of .tar file to be untarred; String

Carter Rhea
06/22/18
'''
import os
import gzip
import shutil
#------------------------INPUTS-------------------------#
chandra_dir = '/home/carterrhea/Desktop/WVT_test'
dir_list = ['323','324']
#-------------------------------------------------------#

#We need to unzip the files with .gz
def gunzip(file):
    with gzip.open(file+'.fits.gz', 'rb') as f_in:
        with open(file+'.fits', 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)

def unzip(chandra_dir,dir_list):
    for dir in dir_list:
        path = os.getcwd()+'/'+dir
        for dirname in ['primary','secondary']:
            os.chdir(chandra_dir+'/'+dir+'/'+dirname)
            for file in os.listdir(os.getcwd()):
                if file.endswith(".fits.gz"):
                    gunzip(file.split(".")[0])
    return None



def main():
    os.chdir(chandra_dir)
    unzip(chandra_dir,dir_list)
    return None

main()
