'''
Cooling time calculation
Given: Temperature, Luminosity, and Volume
'''

import numpy as np
import scipy.integrate as spi
import scipy.constants as spc
import astropy.constants as cst
from astropy.cosmology import Planck13 as cosmo

#--------------------------INPUT PARAMETERS-----------------------------#
r_in = 0 #in arcsec
r_out = 42 #in arcsec
Temp = 6.07
Norm = 7.7e-3 #average of apec norms
Flux = -12.3#log10 of flux
redshift = 0.702
#-----------------------------------------------------------------------#

#------------------------------CLASS------------------------------------#
class fit:
    def __init__(self,r_in,r_out,temp,temp_min,temp_max,norm,norm_min,norm_max,flux):
        self.r_in = r_in
        self.r_out = r_out
        self.temp = [temp+temp_min,temp,temp+temp_max]
        self.temp_ergs = [val*1.60218e-9 for val in self.temp]
        self.norm = [norm+norm_min,norm,norm+norm_max]
        self.flux = flux
        self.lum = 0
        self.dens = []
        self.press = []
        self.entropy = []
        self.t_cool = []
        self.vol = 0
        self.Da = 0
        self.dl = 0
        self.m_acc_dot = 0
    def calc_D(self,z):
        '''
        Calculate comoving distance in kpc
        '''
        dist = ds_calc(z)
        self.Da = dist[0]*3.086e21 #conversion to cm from kpc
        self.dl = dist[1]*3.086e21 #conversion to cm from kpc
    def calc_lum(self):
        '''
        Calculate luminosity from flux
        '''
        self.lum = 10**(self.flux)*4*np.pi*self.dl**2
    def calc_vol(self,z):
        '''
        Calculate volume given inner and outer radii
        '''
        dist_out = ls_calc(z,self.r_out)*3.086e21 #conversion to cm from kpc
        dist_in = ls_calc(z,self.r_in)*3.086e21 #conversion to cm from kpc
        self.vol = (4/3)*np.pi*(dist_out**3-dist_in**3)
    def calc_dens(self,z):
        '''
        Calculate electron density and errors
        '''
        const = 1e7*np.sqrt(4*np.pi)
        red_dep = self.Da*(1+z)
        for i in range(3):
            norm_vol = np.sqrt((1.2*self.norm[i])/self.vol)
            self.dens.append(const*red_dep*norm_vol)
    def calc_press(self):
        '''
        Calculate pressure and errors
        '''
        self.press = [2*self.temp_ergs[i]*self.dens[i] for i in range(3)]
    def calc_entropy(self):
        '''
        Calculate entropy and errors
        '''
        self.entropy = [self.temp[i]*self.dens[i]**(-2/3) for i in range(3)]
    def calc_tcool(self):
        '''
        Calculate cooling time and errors
        '''
        for i in range(3):
            t_sec = (5/2)*((1.91*self.dens[i]*self.temp_ergs[i])/(self.lum/self.vol))
            self.t_cool.append(t_sec*(1/3.15e16))#3.17098e-8*1e-9) #year -> Gigayears
    def calc_M_acc_dot(self,z):
        '''
        Calculate Accreiton Rate
        '''
        pho = 1.91*self.dens[1]*0.6*(1.00784*cst.u.value)
        alpha = 1.4 #from 1904.08942
        dist_out = ls_calc(z,self.r_out)*3.086e21 #conversion to cm from kpc
        dist_in = ls_calc(z,self.r_in)*3.086e21 #conversion to cm from kpc
        num = (4*np.pi*(dist_out**3-dist_in**3)*pho)/cst.M_sun.value #needs to be in solar masses
        den = alpha*(self.t_cool[1]*3.17098e8) #needs to be in years
        macc = (num/den)
        self.m_acc_dot = num/den


#----------------------------------FUNCTIONS---------------------------------#
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


def ls_calc(z,theta):
    '''
    Calculate the object's physical size which is the angular diameter distance * angular size of object from earth
    PARAMETERS:
        z - redshift
        theta - angle
    '''
    Omega_mass = 0.3#cosmo.Om0
    Omega_lam = 0.7#cosmo.Ode0
    Hubble_const = 70#cosmo.H0
    theta_rad = theta*(np.pi/648000)
    d_A = calc_size(z,Omega_mass,Omega_lam,Hubble_const)
    return d_A*theta_rad

def ds_calc(z):
    '''
    Calculate angular diameter distance and luminosity distance
    PARAMETERS:
        z - redshift
    '''
    Omega_mass = 0.3#cosmo.Om0
    Omega_lam = 0.7#cosmo.Ode0
    Hubble_const = 70#cosmo.H0
    d_A = calc_size(z,Omega_mass,Omega_lam,Hubble_const)
    d_l = d_A*(1+z)**2
    return d_A, d_l

def calc_params(r_in,r_out,Temp,Norm,Flux):
    '''
    Calculate all related parameters for ICM
    '''
    best_fit = fit(r_in,r_out,Temp,0,0,Norm,0,0,Flux)
    best_fit.calc_D(redshift)
    best_fit.calc_lum()
    best_fit.calc_vol(redshift)
    best_fit.calc_dens(redshift)
    best_fit.calc_press()
    best_fit.calc_entropy()
    best_fit.calc_tcool()
    best_fit.calc_M_acc_dot(redshift)
    return best_fit

def main():
    print("Using the bolometric Luminosity we find:")
    best_fit = calc_params(r_in,r_out,Temp,Norm,Flux)
    print(" The Electron Density is %.2E"%best_fit.dens[1])
    print(" The Entropy is %.2E"%best_fit.entropy[1])
    print(" The Luminosity is %.2E"%best_fit.lum)
    print(" The cooling time is %.2E Gyr"%best_fit.t_cool[1])
    print(" The Mass Accretion Rate is %.2E M_sol/yr"%best_fit.m_acc_dot)


    return None
main()
