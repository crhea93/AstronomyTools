'''
Apply processing and GTI and finish cleaning
'''
from ciao_contrib.runtool import *
def Process(filenames,OBSID):
    #Process event
    acis_process_events.punlearn()
    acis_process_events.infile = filenames['evt1_dstrk']
    acis_process_events.outfile = 'acis_new_evt1.fits'
    acis_process_events.badpixfile = filenames['bpix_repro']
    acis_process_events.acaofffile = filenames['asol1']
    acis_process_events.mtlfile = filenames['mtl1']
    acis_process_events.eventdef = ")stdlev1" #Since VFAINT
    acis_process_events.check_vf_pha = True
    acis_process_events.clobber = True
    acis_process_events()

    #filter bad grades
    dmcopy.punlearn()
    dmcopy.infile = 'acis_new_evt1.fits[EVENTS][grade=0,2,3,4,6,status=0]'
    dmcopy.outfile = 'acis_flt_evt1.fits'
    dmcopy.clobber = True
    dmcopy()

    #Apply GTI and Complete Cleaning!
    dmcopy.punlearn()
    dmcopy.infile = 'acis_flt_evt1.fits[EVENTS][@'+filenames['flt1']+'][cols -phas]'
    dmcopy.outfile = 'acisf'+OBSID+'_repro_evt2.fits'
    dmcopy.clobber = True
    dmcopy()

    return None
