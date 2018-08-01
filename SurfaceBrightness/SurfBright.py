'''
Small script to effective energy
'''
import os
from ciao_contrib.runtool import *

def calc_effenergy(region_name,background_name,energy_range,ra_central,dec_central):
    dmtcalc.infile = region_name+'.arf'
    dmtcalc.outfile = "arf_weights"+str(region_name)
    dmtcalc.expression = "mid_energy=(energ_lo+energ_hi)/2.0;weights=(mid_energy*specresp)"
    dmtcalc.clobber =True
    dmtcalc()
    dmstat.infile = "arf_weights"+str(region_name)+"[mid_energy="+str(energy_range)+"][cols weights]"
    dmstat.verbose = True
    dmstat()
    weight_sum = float(dmstat.out_sum)
    dmstat.infile = "arf_weights"+str(region_name)+"[mid_energy="+str(energy_range)+"][cols specresp]"
    dmstat.verbose = True
    dmstat()
    specresp_sum = float(dmstat.out_sum)
    eff_energy = weight_sum/specresp_sum
    print("Our effective energy is: "+str(eff_energy))
    for file in os.listdir(os.getcwd()):
        if file.endswith("_evt2.fits"):
            event_file = file
    srcflux.infile = event_file
    srcflux.pos = ra_central+" "+dec_central
    srcflux.outroot = region_name+"/"
    srcflux.bands = energy_range+":"+str(eff_energy)
    srcflux.srcreg = region_name+".reg"
    srcflux.bkgreg = background_name+".reg"
    srcflux.clobber = True
    print(srcflux())

    return None


def main():
    chandra_dir = '%%%'
    os.chdir(chandra_dir)
    ra_central = '%%%'
    dec_central = '%%%'
    region_name = '%%%'
    background_name = '%%%'
    energy_range = '%%%'
    calc_effenergy(region_name,background_name,energy_range,ra_central,dec_central)
    return None
main()
