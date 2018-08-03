#!bin/bash

punlearn specextract
pset specextract infile='20528/repro/acisf20528_repro_evt2.fits[sky=region(40kpc.reg)],20940/repro/acisf20940_repro_evt2.fits[sky=region(40kpc.reg)]'
pset specextract outroot='combined/com'
pset specextract weight='yes'
pset specextract combine='yes'
specextract
