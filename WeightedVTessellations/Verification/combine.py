'''
Short python script to combine bins in order to validate program
'''
from ciao_contrib.runtool import *
import os

dir = '/home/crhea/Documents/TestData/12036/repro/binned'
os.chdir(dir)
spec_to_comb = ''
for file in os.listdir(dir):
    if file.endswith(".pi"):
        spec_to_comb += file+','
combine_spectra.clobber = True
combine_spectra.verbose = 5
combine_spectra(src_spectra=spec_to_comb, outroot = 'combined')
