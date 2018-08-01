#!bin/bash
punlearn dmcopy 
pset dmcopy infile="acisf21129_repro_evt2.fits[sky=region(40kpc.reg)][energy=500:2000]"
pset dmcopy outfile=40kpc_soft.fits
pset dmcopy clobber=yes
dmcopy

punlearn dmcopy 
pset dmcopy infile="acisf21129_repro_evt2.fits[sky=region(400kpc.reg)][energy=500:2000]"
pset dmcopy outfile=400kpc_soft.fits
pset dmcopy clobber=yes
dmcopy

punlearn srcflux
pset srcflux infile="40kpc_soft.fits"
pset srcflux pos="10:49:22.6 +56:40:32.5"
pset srcflux bands="soft,medium"
pset srcflux outroot=40kpc/
pset srcflux clobber=yes
srcflux

punlearn srcflux
pset srcflux infile="400kpc_soft.fits"
pset srcflux pos="10:49:22.6 +56:40:32.5"
pset srcflux bands="soft,medium"
pset srcflux outroot=400kpc/
pset srcflux clobber=yes
srcflux
