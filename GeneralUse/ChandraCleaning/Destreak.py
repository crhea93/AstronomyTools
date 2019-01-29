from ciao_contrib.runtool import *

def Destreak(base_dir,output_dir,filenames):
    evt1_name = filenames['evt1_deflared'].split('.')[0].split("/")[-1]
    destreak.punlearn()
    destreak.infile = filenames['evt1_deflared']
    destreak.outfile = base_dir+'/'+output_dir+'/'+evt1_name+'_dstrk.'+filenames['evt1'].split('.')[1]
    destreak.mask = None #All events are evaluated as potential
    destreak.filter = False #Flags but does not remove streak events
    destreak.clobber = True
    destreak()
    #filenames['evt1'].split('.')[1] == .fits
    filenames['evt1_dstrk'] = base_dir+'/'+output_dir+'/'+evt1_name+'_dstrk.'+filenames['evt1'].split('.')[1]
    return None
