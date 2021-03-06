'''
Luminosity Calculator
'''
'''
Python Script to calculate the linear size (the size of the object at a given distance)
given simply the angular size and redshift
INPUTS:
    z -- redshift value
    theta -- angular seperation of image in arcseconds
'''
import sys
import numpy as np
import scipy.integrate as spi
import scipy.constants as spc
from astropy.cosmology import Planck13 as cosmo


def Energy_func_inv(z,Omega_mass,Omega_lam,Omega_k):
    '''
    Classic cosmological energy function. Note we simply invert it here.
    PARAMETERS:
        z - redshift
        Omega_mass - relative mass density in universe
        Omega_lam - relative dark energy density in universe
        Omega_k - 1-(Omega_mass+Omega_lam)
    '''
    return 1/(np.sqrt(Omega_mass*(1+z)**3+Omega_k*(1+z)**2+Omega_lam))
def calc_size(z,Omega_mass,Omega_lam,Hubble_const):
    '''
    Calculate the size of an object at a certain redshift
    PARAMETERS:
        z - redshift
        Omega_mass - relative mass density in universe
        Omega_lam - relative dark energy density in universe
        Omega_k - 1-(Omega_mass+Omega_lam)
    '''
    Omega_K = 1-Omega_mass-Omega_lam
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
    '''
    Calculate comoving distance
    PARAMETERS:
        z - redshift
    '''
    Omega_mass = cosmo.Om0
    Omega_lam = cosmo.Ode0
    Hubble_const = cosmo.H0
    d_A = calc_size(z,Omega_mass,Omega_lam,Hubble_const.value)
    d_l = d_A*(1+z)**2
    return  d_l


def calc_lum(flux,z):
    '''
    Calculate Luminosity given flux and redshift
    PARAMETERS:
        flux - Log of Flux in ergs/s
        z - redshift
    '''
    lum_dist = ds_calc(z)*3.086e21 #kpc to cm
    lum = 10**(flux)*4*np.pi*lum_dist**2
    print("The Calculated Luminosty is %.2E"%lum)
    return lum

def main():
    calc_lum(float(sys.argv[1]),float(sys.argv[2]))

main()
