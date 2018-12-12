'''
Small script to merge objects

INPUTS:
    chandra_dir -- chandra data directory (e.g. '/home/user/Documents/ChandraData/')
    OBSIDS -- list of OBSIDS to merge (e.g. '[11111,22222]')
    reproccesed_dir -- name of reprocessed directory containing event files(e.g. 'repro')
    Merged_Folder -- Name of Folder to be created with merged data (e.g. 'Merged')
    clean -- Boolean to clean files created from 'merge_obs' command (e.g. 'yes')
'''
import os
from astropy.io import fits
from ciao_contrib.runtool import *

#-------------------INPUTS---------------------#
chandra_dir = '%%%'
OBSIDS = ['###','###','###','###']
reproccesed_dir =  'repro_vfaint'
Merged_Folder = 'Merged_all'
clean = 'yes' # yes or no
#----------------------------------------------#

def merge_objects(Obsids,reproccesed_dir,Merged_Folder,clean):
    #Lets create the merged folder if necessary
    if os.path.exists(os.getcwd()+"/"+Merged_Folder):
        pass
    else:
        id_string = ''
        id_hyphen = ''
        for obsid in Obsids:
            id_string += obsid+'/'+reproccesed_dir+"/acisf"+obsid+"_repro_evt2.fits,"
            id_hyphen += obsid+"-"
        os.system("merge_obs '"+id_string+"' "+Merged_Folder+"/ cleanup="+clean )
        #Get rid of last hyphen...
        id_string = id_string[:-1]
        id_hyphen = id_hyphen[:-1]
        #We need to update the HEADER name
        hdu = fits.open(Merged_Folder+"/merged_evt.fits",mode='update')
        hdr = hdu[0].header
        hdr['OBS_ID'] = id_hyphen
        hdr['SEQ_NUM'] = id_hyphen
        hdr['DS_IDENT'] = id_hyphen
        hdu.close()





    return None

def main():
    os.chdir(chandra_dir)
    merge_objects(OBSIDS,reproccesed_dir,Merged_Folder,clean)
    return None

main()
