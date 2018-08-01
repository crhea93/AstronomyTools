#!/bin/bash

#Runs specextract on data

punlearn specextract
pset specextract infile='repro/acisf21129_repro_evt2.fits[sky=region(repro/simple.reg)]'
pset specextract outroot=simple
pset specextract bkgfile='repro/acisf21129_repro_evt2.fits[sky=region(repro/simple_bkg.reg)]'
pset specextract asp=primary/pcadf648564790N001_asol1.fits
pset specextract mskfile=secondary/acisf21129_000N001_msk1.fits
pset specextract badpixfile=repro/acisf21129_repro_bpix1.fits
pset specextract grouptype=NUM_CTS
pset specextract binspec=100
pset specextract clobber=yes
specextract mode=h
