'''
Analysis of aprates data to calculate surface brightness concentration bounds

INPUTS:
    chandra_dir -- full path to directory contianing data (i.e. '/home/user/Documents/Data')
    region -- name of region file of interest without .reg extension (e.g. 'simple')
    quantity_to_calc -- string acronym for quantity to calculate (e.g. 'NEFA')

Options for quantity to calculate:
    NC - Net Counts
    NCR - Net Count Rates
    NPF - Net Photon Flux
    NEFA - Net Energy Flux option A
    NEFB - Net Energy Flux option B

OUTPUTS:
    Prints value and confidence bounds in terminal
'''
import os

def calc_bounds(region,quantity_to_calc):
    with open('aprates_'+region+'.par') as f:
        data = []
        count = 0
        for line in f:
            if count < 35:
                data.append(line.split(',')[3])
            count += 1
    if quantity_to_calc == 'NC':
        val = float(data[0])
        lower = float(data[1])
        upper = float(data[2])
        print("Net Counts is calculated at %.2E with an lower bound of %.2E and an upper bound of %.2E"%(val,lower,upper))
    if quantity_to_calc == 'NCR':
        val = float(data[7])
        lower = float(data[8])
        upper = float(data[9])
        print("Net Count Rate is calculated at %.2E with an lower bound of %.2E and an upper bound of %.2E"%(val,lower,upper))
    if quantity_to_calc == 'NPF':
        val = float(data[15])
        lower = float(data[16])
        upper = float(data[17])
        print("Net Photon Flux is calculated at %.2E with an lower bound of %.2E and an upper bound of %.2E"%(val,lower,upper))
    if quantity_to_calc == 'NEFA' or quantity_to_calc == 'NEFB':
        val = float(data[28])
        src_val = float(data[7])
        src_rates_lower = float(data[8])
        src_rates_upper = float(data[9])
        lower = val*(src_rates_lower/src_val)
        upper = val*(src_rates_upper/src_val)
        print("Net Energy Flux is calculated at %.2E with an lower bound of %.2E and an upper bound of %.2E"%(val,lower,upper))

    return None

def main():
    chandra_dir = '%%%'
    os.chdir(chandra_dir)
    region = '%%%'
    quantity_to_calc = 'NEFB'
    calc_bounds(region,quantity_to_calc)
main()
