'''
This small script is used to create a conglomeration of annuli centered about a point
for use in constructing a radial profile of surface brightness
'''

import os
import numpy as np


def create_ann(ra,dec):
    with open('annuli_nonincl.reg','w+') as file:
        file.write("# Region file format: DS9 version 4.1 \n")
        for i in np.arange(0,50,1):
            if i==1:
                file.write("annulus("+ra+","+dec+",%s,%s) \n"%(str(0)+'"',str(1.0)+'"'))
            else:
                file.write("annulus("+ra+","+dec+",%s,%s) \n"%(str(1.0+i)+'"',str(1.0+(i+1))+'"'))

    with open('annuli.reg','w+') as file:
        file.write("# Region file format: DS9 version 4.1 \n")
        for i in np.arange(0,40,1):
            file.write("annulus("+ra+","+dec+",%s,%s) \n"%(str(0)+'"',str(4.16+(i*4.16))+'"'))
