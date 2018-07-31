'''
Python Script to calculate the angular given simply the linear size and redshift
INPUTS:
    z -- redshift value
    l -- angular seperation in kpc
'''
import sys
import numpy as np
import scipy.integrate as spi
import scipy.constants as spc
from astropy.cosmology import Planck13 as cosmo


def Energy_func_inv(z,Omega_rel,Omega_mass,Omega_lam,Omega_k):
    return 1/(np.sqrt(Omega_rel*(1+z)**4+Omega_mass*(1+z)**3+Omega_k*(1+z)**2+Omega_lam))
def calc_ang(z,l,Omega_rel,Omega_mass,Omega_lam,Hubble_const):
    Omega_K = 1-Omega_mass-Omega_lam
    d_H = spc.c/Hubble_const
    #print(spi.quad(Energy_func_inv,0,z,args=(Omega_rel,Omega_mass,Omega_lam,Omega_K))[0])
    d_C = d_H*spi.quad(Energy_func_inv,0,z,args=(Omega_rel,Omega_mass,Omega_lam,Omega_K))[0]
    if Omega_K > 0:
        d_M = (d_H/np.sqrt(Omega_K))*np.sinh(np.sqrt(Omega_K)*d_C/d_H)
    if Omega_K == 0:
        d_M = d_C
    if Omega_K < 0:
        d_M = (d_H/np.sqrt(np.abs(Omega_K)))*np.sin(np.sqrt(np.abs(Omega_K))*d_C/d_H)
    d_A = d_M/(1+z)
    theta_rad = l/d_A
    theta = theta_rad*(648000/np.pi) #degree
    return theta


def main():
    z = input("Please input a redshift value: ")#1.7089
    z = float(z)
    l = input("Please input a linear distance in kpc: ") #kpc
    l = float(l)

    Omega_rel = cosmo.Onu0
    Omega_mass = cosmo.Om0
    Omega_lam = cosmo.Ode0
    Hubble_const = cosmo.H0.value


    theta = calc_ang(z,l,Omega_rel,Omega_mass,Omega_lam,Hubble_const)

    print("Our assumed cosmology is as follows:")
    print("   Omega_rel = "+str(Omega_rel))
    print("   Omega_mass = "+str(Omega_mass))
    print("   Omega_lam = "+str(Omega_lam))
    print("   Hubble_const = "+str(Hubble_const))
    print()
    print("We have an angular distance of theta = "+str(theta))
main()
