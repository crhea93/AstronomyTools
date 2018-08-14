'''
Analysis of aprates data to calculate surface brightness concentration bounds

INPUTS:
    chandra_dir -- full path to directory contianing data (i.e. '/home/user/Documents/Data')
    OBSIDs -- list of obsids to loop through to calculate CSB bounds (e.g. ['#####'])
    repro_dirs -- list of reprocessed directories to calc CSB bounds (e.g. ['repro','repro_vfaint'])
    region1 -- name of first region file of interest without .reg extension (e.g. 'simple')
    region2 -- name of second region file of interest without .reg extension (e.g. 'simple2')
    quantities_to_calc -- list of string acronyms for quantity to calculate (e.g. ['NEFA'])

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

#------------------INPUTS------------------------------------------------------#
chandra_dir = '%%%'
OBSIDs = ['%%%']
repro_dirs = ['%%%']
region1 = '%%'
region2 = '%%'
quantities_to_calc = ['%%%','%%%']
#------------------------------------------------------------------------------#

def isfloat(value):
  try:
    float(value)
    return True
  except ValueError:
    return False
def div_with_err(val1,val2):
    try:
        value = val1/val2
        return value
    except ZeroDivisionError:
        return 0.0
def calc_bounds(region1,region2,quantity_to_calc,fout,fout_all):
    with open('aprates_'+region1+'.par') as f:
        data = []
        count = 0
        for line in f:
            if count < 35:
                val = line.split(',')[3]
                if isfloat(val) == False:
                    data.append(0)
                else:
                    data.append(val)
            count += 1
    with open('aprates_'+region2+'.par') as f:
        data2 = []
        count = 0
        for line in f:
            if count < 35:
                val = line.split(',')[3]
                if isfloat(val) == False:
                    data2.append(0)
                else:
                    data2.append(val)
            count += 1
    if quantity_to_calc == 'NC':
        val_1 = float(data[0])
        lower_1 = float(data[1])
        upper_1 = float(data[2])
        val_2 = float(data2[0])
        lower_2 = float(data2[1])
        upper_2 = float(data2[2])
        csb_val = div_with_err(val_1,val_2)
        csb_lower = div_with_err(lower_1,upper_2)
        csb_upper = div_with_err(upper_1,lower_2)
        #print("Net Counts Concentration is calculated at %.2E with an lower bound of %.2E and an upper bound of %.2E"%(csb_val,csb_lower,csb_upper))
        fout.write("Net Counts Concentration is calculated at %.2E with an lower bound of %.2E and an upper bound of %.2E "%(csb_val,csb_lower,csb_upper))
        fout.write("\n")
        fout_all.write("    Net Counts Concentration is calculated at %.2E with an lower bound of %.2E and an upper bound of %.2E "%(csb_val,csb_lower,csb_upper))
        fout_all.write("\n")
    if quantity_to_calc == 'NCR':
        val_1 = float(data[7])
        lower_1 = float(data[8])
        upper_1 = float(data[9])
        val_2 = float(data2[7])
        lower_2 = float(data2[8])
        upper_2 = float(data2[9])
        csb_val = div_with_err(val_1,val_2)
        csb_lower = div_with_err(lower_1,upper_2)
        csb_upper = div_with_err(upper_1,lower_2)
        #print("Net Count Rate Concentration is calculated at %.2E with an lower bound of %.2E and an upper bound of %.2E"%(csb_val,csb_lower,csb_upper))
        fout.write("Net Count Rate Concentration is calculated at %.2E with an lower bound of %.2E and an upper bound of %.2E"%(csb_val,csb_lower,csb_upper))
        fout.write('\n')
        fout_all.write("    Net Count Rate Concentration is calculated at %.2E with an lower bound of %.2E and an upper bound of %.2E"%(csb_val,csb_lower,csb_upper))
        fout_all.write("\n")
    if quantity_to_calc == 'NPF':
        val_1 = float(data[14])
        lower_1 = float(data[15])
        upper_1 = float(data[16])
        val_2 = float(data2[14])
        lower_2 = float(data2[15])
        upper_2 = float(data2[16])
        csb_val = div_with_err(val_1,val_2)
        csb_lower = div_with_err(lower_1,upper_2)
        csb_upper = div_with_err(upper_1,lower_2)
        #print("Net Photon Flux Concentration is calculated at %.2E with an lower bound of %.2E and an upper bound of %.2E"%(csb_val,csb_lower,csb_upper))
        fout.write("Net Photon Flux Concentration is calculated at %.2E with an lower bound of %.2E and an upper bound of %.2E"%(csb_val,csb_lower,csb_upper))
        fout.write('\n')
        fout_all.write("    Net Photon Flux Concentration is calculated at %.2E with an lower bound of %.2E and an upper bound of %.2E"%(csb_val,csb_lower,csb_upper))
        fout_all.write("\n")
    if quantity_to_calc == 'NEFA':
        val_1 = float(data[21])
        lower_1 = float(data[22])
        upper_1 = float(data[23])
        val_2 = float(data2[21])
        lower_2 = float(data2[22])
        upper_2 = float(data2[23])
        csb_val = div_with_err(val_1,val_2)
        csb_lower = div_with_err(lower_1,upper_2)
        csb_upper = div_with_err(upper_1,lower_2)
        #print("Net Energy Flux Concentration is calculated at %.2E with an lower bound of %.2E and an upper bound of %.2E"%(csb_val,csb_lower,csb_upper))
        fout.write("Net Energy Flux Concentration is calculated at %.2E with an lower bound of %.2E and an upper bound of %.2E"%(csb_val,csb_lower,csb_upper))
        fout.write('\n')
        fout_all.write("    Net Energy Flux Concentration is calculated at %.2E with an lower bound of %.2E and an upper bound of %.2E"%(csb_val,csb_lower,csb_upper))
        fout_all.write("\n")
    if  quantity_to_calc == 'NEFB':
        val_1 = float(data[28])
        src_val_1 = float(data[7])
        src_rates_lower_1 = float(data[8])
        src_rates_upper_1 = float(data[9])
        lower_1 = val_1*(src_rates_lower_1/src_val_1)
        upper_1 = val_1*(src_rates_upper_1/src_val_1)
        val_2 = float(data2[28])
        src_val_2 = float(data2[7])
        src_rates_lower_2 = float(data2[8])
        src_rates_upper_2 = float(data2[9])
        lower_2 = val_2*(src_rates_lower_2/src_val_2)
        upper_2 = val_2*(src_rates_upper_2/src_val_2)
        csb_val = div_with_err(val_1,val_2)
        csb_lower = div_with_err(lower_1,upper_2)
        csb_upper = div_with_err(upper_1,lower_2)
        #print("Net Energy Flux Concentration is calculated at %.2E with an lower bound of %.2E and an upper bound of %.2E"%(csb_val,csb_lower,csb_upper))
        fout.write("Net Energy Flux Concentration is calculated at %.2E with an lower bound of %.2E and an upper bound of %.2E"%(csb_val,csb_lower,csb_upper))
        fout.write('\n')
        fout_all.write("    Net Energy Flux Concentration is calculated at %.2E with an lower bound of %.2E and an upper bound of %.2E"%(csb_val,csb_lower,csb_upper))
        fout_all.write("\n")

    return None

fout_all = open('CSB.txt',"w")
for OBSID in OBSIDs:
    print('#-------------------------------------------------------------#')
    fout_all.write('#----------------------------OBSID '+OBSID+'---------------------------------# \n')
    print("We are on OBSID: "+OBSID)
    for repro_dir in repro_dirs:
        print("We are on reprocessed directory: "+repro_dir)
        fout_all.write("We are on reprocessed directory: "+repro_dir+" \n")
        os.chdir(chandra_dir+'/'+OBSID+'/'+repro_dir)
        fout = open('CSB.txt','w')
        for quantity_to_calc in quantities_to_calc:
            calc_bounds(region1,region2,quantity_to_calc,fout,fout_all)
        fout.close()
