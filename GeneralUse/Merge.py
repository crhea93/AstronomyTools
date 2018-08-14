'''
Small script to merge objects

INPUTS:
    chandra_dir -- full path to main chandra data (e.g. '/home/user/Documents/ChandraData')
    OBSIDS -- list of obsids to merge (e.g. ['#####','#####'])
    reproccesed_dir -- name of directory containing reprocessed data (e.g. 'repro')
    Merged_Folder -- name for new merged folder (e.g. 'Merged')
'''
import os
from astropy.io import fits
from ciao_contrib.runtool import *


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
    chandra_dir = '%%%'
    os.chdir(chandra_dir)
    OBSIDS = ['%%%%']
    reproccesed_dir =  '%%%'
    Merged_Folder = '%%%'
    clean = 'no' # yes or no
    merge_objects(OBSIDS,reproccesed_dir,Merged_Folder,clean)
    return None

main()
