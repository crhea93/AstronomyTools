'''
Small script to calculate and do analysis for Concentration surface brightness
'''
import numpy as np

flux_40 = 9.77E-16
flux_40_min = 2.14E-16
flux_40_max = 2.3E-15

flux_400 = 1.14E-14
flux_400_min = 5.09E-15
flux_400_max = 1.76E-14

CSB = flux_40/flux_400
CSB_min = flux_40_min/flux_400_min
CSB_max = flux_40_max/flux_400_max

print("The Coefficient of Surface Brightness is: "+str(CSB))
print("It has a lower bound of: "+str(CSB_min))
print("It has an upper bound of: "+str(CSB_max))
