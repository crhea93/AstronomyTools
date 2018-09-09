'''
This python file has the simple task of combining ds9 regions into a single file!
INPUTS:
    chandra_dir -- chandra data directory (e.g. '/home/user/Documents/ChandraData/')
    regs_add -- list of region files to add (e.g. ['100kpc','200kpc'])
    regs_sub -- list of region files to subtract (e.g. ['50kpc'])
    combined_name -- name for combined region file (e.g. 'Combined')
'''
import os
#------------------INPUTS---------------------#
chandra_dir = '/home/user/Documents/Data'
regs_add = ['100kpc', '200kpc']
regs_sub = ['5kpc']
combined_name = 'combined'
#---------------------------------------------#
#Read in first file and keep exact format
os.chdir(chandra_dir)
read_list = []
with open(regs_add[0]+'.reg','r') as file:
    for line in file:
        read_list.append(file)

for region in regs_add[1:]:
    with open(region+'.reg','r') as file:
        for line in file:
            if 'circle' in line:
                read_list.append(line)

for region in regs_sub:
    with open(region+'.reg','r') as file:
        for line in file:
            if 'circle' in line:
                read_list.append('-'+line)

with open(combined_name+'.reg','w+') as fileout:
    for line in read_list:
        fileout.write(line+'\n')
