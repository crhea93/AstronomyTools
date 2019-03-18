'''
Python Script to calculate the linear size (the size of the object at a given distance)
given simply the angular size and redshift
INPUTS:
    z -- redshift value
    theta -- angular seperation of image in arcseconds
'''
import numpy as np
import scipy.integrate as spi
import scipy.constants as spc
from astropy.cosmology import Planck13 as cosmo


def Energy_func_inv(z,Omega_rel,Omega_mass,Omega_lam,Omega_k):
    return 1/(np.sqrt(Omega_rel*(1+z)**4+Omega_mass*(1+z)**3+Omega_k*(1+z)**2+Omega_lam))
def calc_size(z,theta_rad,Omega_rel,Omega_mass,Omega_lam,Hubble_const):
    Omega_K = 1-Omega_mass-Omega_lam
    d_H = spc.c/Hubble_const
    d_C = d_H*spi.quad(Energy_func_inv,0,z,args=(Omega_rel,Omega_mass,Omega_lam,Omega_K))[0]
    if Omega_K > 0:
        d_M = (d_H/np.sqrt(Omega_K))*np.sinh(np.sqrt(Omega_K)*d_C/d_H)
    if Omega_K == 0:
        d_M = d_C
    if Omega_K < 0:
        d_M = (d_H/np.sqrt(np.abs(Omega_K)))*np.sin(np.sqrt(np.abs(Omega_K))*d_C/d_H)
    d_A = d_M/(1+z)
    l = d_A*theta_rad
    print(l)
    return l


def calc_scale(z):
    Omega_rel = cosmo.Onu0
    Omega_mass = cosmo.Om0
    Omega_lam = cosmo.Ode0
    Hubble_const = cosmo.H0

    theta_rad = (np.pi/648000) #radians per theta

    l = calc_size(z,theta_rad,Omega_rel,Omega_mass,Omega_lam,Hubble_const.value)
    return l
