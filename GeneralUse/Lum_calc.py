'''
Program to verify that we are recovering an appropriate luminosity value
'''
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import scipy.integrate as spi
import scipy.constants as spc
Omega_mass = 0.3
Omega_lam = 0.7
Hubble_const = 70
Omega_K = 1 - Omega_mass - Omega_lam

#------------------INPUTS----------------------#
z = 1.7
flux_log = -15.7
#----------------------------------------------#

def Energy_func_inv(z,Omega_mass,Omega_lam,Omega_k):
    return 1/(np.sqrt(Omega_mass*(1+z)**3+Omega_k*(1+z)**2+Omega_lam))
def calc_size(z,Omega_mass,Omega_lam,Hubble_const):
    d_H = spc.c/Hubble_const
    d_C = d_H*spi.quad(Energy_func_inv,0,z,args=(Omega_mass,Omega_lam,Omega_K))[0]
    if Omega_K > 0:
        d_M = (d_H/np.sqrt(Omega_K))*np.sinh(np.sqrt(Omega_K)*d_C/d_H)
    if Omega_K == 0:
        d_M = d_C
    if Omega_K < 0:
        d_M = (d_H/np.sqrt(np.abs(Omega_K)))*np.sin(np.sqrt(np.abs(Omega_K))*d_C/d_H)
    d_A = d_M/(1+z)
    return d_A


def ds_calc(z):
    d_A = calc_size(z,Omega_mass,Omega_lam,Hubble_const)
    d_l = d_A*(1+z)**2
    d_l = d_l* 3.086e21  # conversion to cm from kpc
    return d_l


def Flux_to_Lum(flux,z):
    #Luminosity/E(z)
    return 10 ** (flux) * 4 * np.pi * ds_calc(z) ** 2



Lum = Flux_to_Lum(flux_log,z)
print("For redshift %.1f we have a luminosity of %.2E"%(z,Lum))
