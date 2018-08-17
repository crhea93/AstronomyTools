'''
Simulate PSF using marx
    parameters:
        evt_file - associated event file (string)
        psf_file - associated psf fits file (string)
        outname - name of .psf file to be created (string)
    outputs:
        Creates a .psf for event
'''

import os
from ciao_contrib.runtool import *


def simulatePSF(evt_file,psf_file,outname):
    #Get source RA and dec from PSF
    temp_a = os.popen('dmlist '+psf_file+' header | grep SRC_RA').readlines()
    ra = temp_a[0].split("[deg]")[0].split("SRC_RA")[1].strip()
    temp_b = os.popen('dmlist '+psf_file+' header | grep SRC_DEC').readlines()
    dec = temp_b[0].split("[deg]")[0].split("SRC_DEC")[1].strip()
    #Simulate PSF
    simulate_psf.punlearn()
    simulate_psf.infile = evt_file+'.fits'
    simulate_psf.outroot = outname
    simulate_psf.ra = ra
    simulate_psf.dec = dec
    simulate_psf.simulator = 'file'
    simulate_psf.rayfile = psf_file
    simulate_psf.projector = 'marx'
    simulate_psf()
    return None
#-------------------------------------------------#
#-------------------------------------------------#