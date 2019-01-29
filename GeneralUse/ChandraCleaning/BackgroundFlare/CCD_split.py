'''
Create ccd specific fits files
'''

import os
from ciao_contrib.runtool import *

#-----------------INPUTS------------------#
chandra_dir = '/home/carterrhea/Documents/Data'
dir_to_split = ['20528','20941','20940','21129']
ccds = ['0','1','2','3']
background_dir = 'Background'
#-----------------------------------------#

def main():
    os.chdir(chandra_dir)
    if not os.path.exists(chandra_dir+'/'+background_dir):
        os.mkdir(chandra_dir+'/'+background_dir)
    for dir in dir_to_split:
        if not os.path.exists(dir+'/'+background_dir):
            os.makedirs(dir+'/'+background_dir)
        os.chdir(os.getcwd()+'/'+dir+'/primary')
        evt2_file = None
        for file in os.listdir(os.getcwd()):
            if file.endswith("_evt2.fits"):
                evt2_file = os.getcwd() + '/'+ file
        os.chdir(chandra_dir+'/'+dir+'/'+background_dir)
        for ccd in ccds:
            dmcopy.punlearn()
            dmcopy.infile = evt2_file+'[ccd_id='+ccd+']'
            dmcopy.outfile = 'ccd'+ccd+'.fits'
            dmcopy.clobber = True
            dmcopy()
        os.chdir(chandra_dir)
    return None

main()
