from ciao_contrib.runtool import *

def BadPixel(base_dir,output_dir,OBSID,filenames,num_of_biases):
    #Make parameter list
    par_file = OBSID+"_obs.par"
    dmmakepar.punlearn()
    dmmakepar.input = filenames['evt1_dstrk']
    dmmakepar.output = par_file
    dmmakepar.clobber = True
    dmmakepar()
    #Create bias list
    lis_file = base_dir+'/'+output_dir+'/bias.lis'
    with open(lis_file,'w+') as lis:
        for i in range(num_of_biases):
            lis.write(filenames[str(i+1)+'_bias0']+'\n')

    #Now create initial badpix file
    bpix_file1 = 'acisf'+OBSID+'_abb1_bpix.fits'
    acis_build_badpix.punlearn()
    acis_build_badpix.obsfile = par_file
    acis_build_badpix.pbkfile = filenames['pbk0']
    acis_build_badpix.biasfile = '@'+lis_file
    acis_build_badpix.outfile = bpix_file1
    acis_build_badpix.bitflag = '00000000000000022221100020022212'
    acis_build_badpix.calibfile = 'CALDB'
    acis_build_badpix.clobber = True
    acis_build_badpix()

    #Find afterglow events
    aglow_file = 'acisf'+OBSID+'_aglow_bpix1.fits'
    acis_find_afterglow.punlearn()
    acis_find_afterglow.infile = filenames['evt1_dstrk']
    acis_find_afterglow.outfile = aglow_file
    acis_find_afterglow.badpixfile = bpix_file1
    acis_find_afterglow.maskfile = filenames['msk1']
    acis_find_afterglow.statfile = filenames['stat1']
    acis_find_afterglow.clobber = True
    acis_find_afterglow()

    #Complete bpix file
    bpix_repro = 'acisf'+OBSID+'_repro_bpix.fits'
    acis_build_badpix.obsfile = par_file
    acis_build_badpix.pbkfile = filenames['pbk0']
    acis_build_badpix.biasfile = None
    acis_build_badpix.outfile = bpix_repro
    acis_build_badpix.calibfile = aglow_file
    acis_build_badpix.procbias = False
    acis_build_badpix.clobber = True
    acis_build_badpix()

    filenames['bpix_repro'] = bpix_repro
    return None
