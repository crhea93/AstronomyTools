'''
Small script to merge objects

INPUTS:
    chandra_dir -- full path to data directory (e.g. '/home/user/Documents/Data')
    OBSIDS -- list of obs_ids to merge (e.g. ['#####','#####'])
    Merged_Folder -- name of new merged folder (e.g. 'Merged')
'''
import os
from astropy.io import fits
from ciao_contrib.runtool import *


def merge_objects(Obsids,Merged_Folder,filenames):
    #Lets create the merged folder if necessary
    if os.path.exists(os.getcwd()+"/"+Merged_Folder):
        pass
    else:
        id_string = ''
        id_hyphen = ''
        for obsid in Obsids:
            id_string += obsid+" "
            id_hyphen += obsid+"-"
        os.system("merge_obs '"+id_string+"' "+Merged_Folder+"/")
        #Get rid of last hyphen...
        id_hyphen = id_hyphen[:-1]
        #We need to update the HEADER name
        hdu = fits.open(Merged_Folder+"/merged_evt.fits",mode='update')
        hdr = hdu[0].header
        print(hdr)
        hdr['OBS_ID'] = id_hyphen
        hdr['SEQ_NUM'] = id_hyphen
        hdr['DS_IDENT'] = id_hyphen
        print(hdr)
        hdu.close()



    return None

def main():
    chandra_dir = '%%%'
    os.chdir(chandra_dir)
    OBSIDS = ['%%%','%%%']
    Merged_Folder = '%%%'
    merge_objects(OBSIDS,Merged_Folder)
    return None

main()
