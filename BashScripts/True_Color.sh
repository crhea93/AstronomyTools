#!/bin/bash

punlearn dmcopy 
pset dmcopy infile="repro/acisf21129_repro_evt2.fits[energy=200:1500][bin x=::1,y=::1]"
pset dmcopy outfile=soft_img.fits
pset dmcopy clobber=yes
dmcopy

punlearn dmcopy 
pset dmcopy infile="repro/acisf21129_repro_evt2.fits[energy=1500:2500][bin x=::1,y=::1]"
pset dmcopy outfile=med_img.fits
pset dmcopy clobber=yes
dmcopy

punlearn dmcopy 
pset dmcopy infile="repro/acisf21129_repro_evt2.fits[energy=2500:8000][bin x=::1,y=::1]"
pset dmcopy outfile=hard_img.fits
pset dmcopy clobber=yes
dmcopy

punlearn dmimg2jpg
pset dmimg2jpg infile=soft_img.fits
pset dmimg2jpg greenfile=med_img.fits
pset dmimg2jpg bluefile=hard_img.fits
pset dmimg2jpg outfile=truecolor.jpg
pset dmimg2jpg clobber=yes
dmimg2jpg
