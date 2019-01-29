'''
Python script to create lightcurves
'''
import os
from ciao_contrib.runtool import *
from pychips import *
from pychips.hlui import *
from pycrates import *
from lightcurves import *
import numpy as np
#------------------------------INPUTS-------------------------------#
chandra_dir = '/home/usr/Documents/Data/21129/Background'
file_lc = ['ccd2_bkg']
time_bin = '200'
#-------------------------------------------------------------------#

def main():
    if not os.path.exists(chandra_dir):
        print("Please Create a Directory Called Background in the OBSID of Interest...")
        print("This Folder Needs to Contain the fits Files with Sources Excluded...")
        print("The Program is Now Exiting...")
        exit()
    for file_ in file_lc:
        if not os.path.isfile(chandra_dir+'/'+file_+'.fits'):
            print("Missing File: "+file_)
            print("Please Add File to Background Directory")
            print("The Program is Now Exiting...")
            exit()

    os.chdir(chandra_dir)
    for file in file_lc:
        #Create Lightcurve
        dmextract.punlearn()
        dmextract.infile = file+'.fits[bin time=::'+time_bin+']'
        dmextract.outfile = file+'.lc'
        dmextract.opt = 'ltc1'
        dmextract.clobber = True
        dmextract()
        #Plot Lightcurve using CHIPS
        make_figure(file+".lc[cols dt, count_rate]")
        set_curve(["symbol.style", "none"])
        set_plot_title("Light Curve")
        set_plot_xlabel(r"\Delta T (s)")
        set_plot_ylabel("Rate (count s^{-1})")
        set_preference("export.clobber", "yes")
        print_window(file+"_lc.pdf")
        clear()
        add_window()
        '''bkg = read_file(file+".lc")
        dt = get_colvals(bkg,'dt')
        count_rate = get_colvals(bkg,'count_rate')'''
        lc_sigma_clip(file+'.lc',file+"_clean.gti",sigma=10,pattern="none")
        print_window(file + '_cleanedLC.pdf')
        clear()
        add_window()
    return None

main()
