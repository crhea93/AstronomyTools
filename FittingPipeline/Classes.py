'''
Annulus class definition
'''
import numpy as np
from scipy import interpolate
from LSCalc import ls_calc, ds_calc

class Annulus:
    '''
    Class for a region including the inner/outer radius and all relevant parameters
    :param r_in - inner radius
    :param r_out - outer radius
    :param temp - temperature
    :param temp_min - min temp value
    :param temp_max - max temp value
    :param Ab - abundace
    :param Ab_min - min ab value
    :param Ab_max - max ab value
    :param norm - temperature model normalization value
    :param norm_min - min norm value
    :param norm_max - max norm value
    :param flux - flux value
    :param agn_act - notate the use of AGN in fit or no
    '''
    def __init__(self,r_in,r_out):
        self.r_in = r_in
        self.r_out = r_out
        self.temp = []
        self.temp_ergs = []
        self.abund = []
        self.norm = []
        self.flux = 0
        self.agn_act = False
        self.deproj = True
        self.lum = 0
        self.dens = []
        self.press = []
        self.entropy = []
        self.t_cool = []
        self.vol = 0
        self.Da = 0
        self.dl = 0
    def add_fit_data(self,temp,temp_min,temp_max,abund,abund_min,abund_max,norm,norm_min,norm_max,flux,reduced_chi_sq,agn,deproj,redshift):
        self.temp = [temp_min,temp,temp_max]
        self.temp_ergs = [val*1.60218e-9 for val in self.temp]
        self.abund = [abund_min,abund,abund_max]
        self.norm = [abs(norm_min),norm,abs(norm_max)]
        self.flux = flux
        self.agn_act = agn
        self.deproj = deproj
        self.calc_all(redshift)
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
        Must multiply radii by 0.492 to convert from physical to WCS coordinates (Obsolete -- done earlier)
        '''
        dist_out = self.r_out*3.086e21#ls_calc(z,self.r_out*0.492)*3.086e21 #conversion to cm from kpc
        dist_in = self.r_in*3.086e21#ls_calc(z,self.r_in*0.492)*3.086e21 #conversion to cm from kpc
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
        print(self.dens)
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
            t_sec = (5/2)*((1.91*self.dens[i]*self.temp_ergs[i]*self.vol)/self.lum)
            self.t_cool.append(t_sec*(1/3.15e16)) #sec -> Gigayears

    def calc_all(self,redshift):
        '''
        Calculate all additional PARAMETERS
        '''
        self.calc_D(redshift)
        self.calc_lum()
        self.calc_vol(redshift)
        self.calc_dens(redshift)
        self.calc_press()
        self.calc_entropy()
        self.calc_tcool()

    def save_data(self,file_to_write,file_min,file_max):
        '''
        save data to csv files
        '''
        #now write to files
        agn_corr = 'no'
        if self.agn_act == True:
            agn_corr = 'yes'
        file_to_write.write(str(self.r_in)+","+ str(self.r_out) + "," + str(self.temp[1]) + "," + str(self.abund[1])+ "," + str(self.dens[1])+ "," + str(self.press[1])
            + "," + str(self.entropy[1]) + "," + str(self.t_cool[1])+ ','+ str(agn_corr)+"\n")
        file_min.write(str(self.r_in)+","+ str(self.r_out) + "," + str(self.temp[0]) + "," + str(self.abund[0])+ "," + str(self.dens[0])+ "," + str(self.press[0])
            + "," + str(self.entropy[0]) + "," + str(self.t_cool[0])+  ','+ str(agn_corr)+"\n")
        file_max.write(str(self.r_in)+","+ str(self.r_out) + "," + str(self.temp[2]) + "," + str(self.abund[2])+ "," + str(self.dens[2])+ "," + str(self.press[2])
            + "," + str(self.entropy[2]) + "," + str(self.t_cool[2])+  ','+ str(agn_corr)+"\n")
