import os
from ciao_contrib.runtool import *

def Flares(flare_gti,base_dir,output_dir,filenames):
    os.chdir(os.getcwd()+'/Background')
    #Determine Goodtime intervals
    evt1_name = filenames['evt1'].split('.')[0]
    dmcopy.punlearn()
    dmcopy.clobber = True
    dmcopy.infile = evt1_name+'.fits[@'+flare_gti+']'
    dmcopy.outfile = base_dir+'/'+output_dir+'/'+evt1_name.split('/')[-1]+'_deflared.fits'
    dmcopy()
    #filenames['evt1'].split('.')[1] == .fits
    filenames['evt1_deflared'] = base_dir+'/'+output_dir+'/'+evt1_name.split('/')[-1]+'_deflared.'+filenames['evt1'].split('.')[1]
    #Clear status bits for future cleaning
    os.system("acis_clear_status_bits " + filenames['evt1_deflared'])
    return None
