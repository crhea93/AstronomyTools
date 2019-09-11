'''
Python Script that will calculate and apply the necessary astrometic corrections given the reference and target sources
'''

import numpy as np
import astropy.units as u
from astropy.coordinates import SkyCoord,FK5
from astropy.io import fits
import matplotlib.pyplot as plt
#Read in coordinates from text file




def read_coord(filename):
    '''
    Read in coordinates (RA/DEC) from text file
    :param filename: Path to fie containing coordinate information
    :return: numpy array of coordinates in the following form -- ['RA DEC',...]
    '''
    coords = []
    with open(filename) as file:
        for line in file.readlines():
            ra = line.split(' ')[0]
            dec = line.split(' ')[1].strip('\n')
            coords.append(ra+' '+dec)
    c = SkyCoord(coords, unit=(u.hourangle, u.deg), frame=FK5, obstime="J2000")
    return c


def correction(ref_file,target_file,target_fits,outfile):
    ref_coords = read_coord(ref_file)
    tar_coords = read_coord(target_file)
    idx, d2d, d3d = tar_coords.match_to_catalog_sky(ref_coords)
    ra_offset = np.mean(tar_coords.ra.deg)-np.mean(ref_coords.ra.deg)
    dec_offset = np.mean(tar_coords.dec.deg) - np.mean(ref_coords.dec.deg)
    #Print out some important information about the offset
    print('The mean RA offset in arcseconds is %.2f'%(np.mean(tar_coords.ra.arcsec)-np.mean(ref_coords.ra.arcsec)))
    print('The mean DEC offset in arcseconds is %.2f' % (np.mean(tar_coords.dec.arcsec) - np.mean(ref_coords.dec.arcsec)))
    #Create new fits file with the offset applied
    new_fits = fits.open(target_fits)
    #update position of center
    new_fits[0].header['CRVAL1'] -= ra_offset
    new_fits[0].header['CRVAL2'] -= dec_offset
    new_fits[0].header['HISTORY'] = 'Shifted image by %s arcsec in ra and %s arcsec in dec'%(ra_offset,dec_offset)
    new_fits.writeto(outfile+'.fits',overwrite=True)
    return None


def select_coords(fits_file):
    '''
    Allow the user to interactively select the point sources using matplotlib
    :param fits_file: input fits file
    :return: sky coordinates of selected sources
    '''
    coords = []
    fits_ = fits.open(fits_file)
    plt.imshow(fits_[0].data)
    plt.show()
    return coords



correction('ref_file.txt','target_file.txt','broad_flux.img','shifted')
