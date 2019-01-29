from ciao_contrib.runtool import *
def Astrometric(OBSID,filenames,source,source_ra,source_dec):
    #Create broad band image
    fluximage.punlearn() #Commented out for now!
    fluximage.infile = filenames['evt2']
    fluximage.outroot = OBSID
    fluximage.binsize = 1
    fluximage.bands = 'broad'
    fluximage.clobber = True
    fluximage()
    #Calculate physical coordinates of source
    dmstat.punlearn()
    dmstat.infile = OBSID+'_broad_thresh.img[sky=region('+source+')]'
    dmstat()
    cntrd_phys = dmstat.out_cntrd_phys.split(',')
    phys_x_obs = float(cntrd_phys[0])
    phys_y_obs = float(cntrd_phys[1])
    #Calculate physical coordinates for actual location
    dmcoords.punlearn()
    dmcoords.infile = OBSID+'_broad_thresh.img'
    dmcoords.option = 'cel'
    dmcoords.celfmt = 'hms'
    dmcoords.ra = source_ra
    dmcoords.dec = source_dec
    dmcoords()
    phys_x_ref = dmcoords.x
    phys_y_ref = dmcoords.y
    diff_x = phys_x_ref-phys_x_obs
    diff_y = phys_y_ref-phys_y_obs
    #Make copies of event file and asol file
    dmcopy.punlearn()
    dmcopy.infile = filenames['evt2']
    dmcopy.outfile = filenames['evt2']+'.corrected'
    dmcopy.clobber = True
    dmcopy()
    dmcopy.punlearn()
    dmcopy.infile = filenames['asol1']
    dmcopy.outfile = filenames['asol1']+'.corrected'
    dmcopy.clobber = True
    dmcopy()
    #Correct event file and asol file
    wcs_update.punlearn()
    wcs_update.infile = filenames['evt2']+'.corrected'
    wcs_update.outfile = ''
    wcs_update.wcsfile = OBSID+'_broad_thresh.img'
    wcs_update.deltax = diff_x
    wcs_update.deltay = diff_y
    wcs_update.clobber = True
    wcs_update()
    wcs_update.punlearn()
    wcs_update.infile = filenames['asol1']+'.corrected'
    wcs_update.outfile = 'acisf'+OBSID+'_new_asol.fits'
    wcs_update.wcsfile = OBSID+'_broad_thresh.img'
    wcs_update.deltax = diff_x
    wcs_update.deltay = diff_y
    wcs_update.clobber = True
    wcs_update()
    return None
