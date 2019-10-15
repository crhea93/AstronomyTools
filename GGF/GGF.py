'''
Construct gaussian gradient filtered images of Chandra data
'''

from astropy.io import fits
import os
import matplotlib.pyplot as plt
import numpy as np
import scipy.ndimage as scim
import cv2
from pylab import rcParams
def ggf1(infile,outfile,sigma):
    '''
    Creates a log-scaled, smoothed, gaussian gradient filtered image (in that order) from a fits file
    :param infile: fits image file to read in
    :param outfile: fits file to create
    :param sigma: sigma value for gaussian used in filtering
    :return: both fits file and png image
    '''
    #First we will just read in the data and get the information we need
    file_temp = fits.open(infile)
    data = file_temp[0].data
    data[data<=0] = np.mean(data)
    data_log =  np.arcsinh(data)#np.log(data) #log data
    #data_log = np.ma.array(data_log,mask=np.isnan(data_log))
    data_filtered = scim.gaussian_gradient_magnitude(data_log, sigma) #filter data
    #Update fits file to save new filtered data
    file_temp[0].data = data_filtered
    if os.path.isfile(outfile+'.img'):
        os.remove(outfile+'.img')
        os.remove(outfile + '.png')
    file_temp.writeto(outfile+'.img')
    plt.imshow(data_filtered)
    plt.colorbar()
    #plt.savefig(outfile+'.png')
    plt.clf()
    return data_filtered

def combine_ggf(img_dir,infiles,radius_bins,weight_bins,center_pixel,fits_file):
    '''
    Combine GGF plots. For each radial region, we choose a weight for the image.
    1. Create mask for each image based of weight_bins and radius_bins
    2. Add weighted images together for each bin
    3. Reconstruct complete weighted image by recombining weighted, binned image
    The image files all need to be the same size and of the same region!
    :param img_dir: Full path to image files
    :param infiles: List of input files -- GGF
    :param radius_bins: List of radii used for binning
    :param weight_bins: List of bin weights corresponding to each GGF image
    :param center_pixel: (X,Y) tuple of central pixel from which the radial bins expand
    :param fits_file: Input Fits Image File for header info
    :return: Reconstructed weighted-GGF image and fits file
    '''
    masked_images = [] #List of 2d numpy arrays. Each item in the list will be one of the masked images
    ct = 0
    (X, Y) = infiles[0].shape
    final_img = np.zeros((X, Y))
    im_blend = np.zeros((X,Y))
    for infile,weight_bin in zip(infiles,weight_bins): #Step through each GGF image, mask it, and add to list
        #Create mask where positive weights map to False because we want to keep them :)
        weight_mask = []
        for weight in weight_bin:
            if weight != 0:
                weight_mask.append(False)
            else:
                weight_mask.append(True)
        #apply weight mask to radius values to only get relevant radii
        radius_masked = np.ma.array(radius_bins,mask=weight_mask)
        #Get rid of nonvalid entries
        radius_masked = radius_masked[~radius_masked.mask]
        #print(radius_masked)
        radii = (radius_masked[0],radius_masked[-1])
        #print(radii)
        rad_mask = make_radial_mask(infile.T,radii,center_pixel)
        masked_images.append(infile*rad_mask)
        plt.imshow(masked_images[ct])
        plt.colorbar()
        '''plt.savefig(img_dir+'img'+str(ct)+'.png')
        if ct >= 1: #Poisson image blending
            im_blend = cv2.seamlessClone(infile, im_blend, rad_mask, center_pixel, cv2.MIXED_CLONE)
        else: #initialize
            im_blend = masked_images[ct]'''
        ct += 1
        #cv2.imshow('stuff',np.matrix(infile))
        plt.clf()
    #Get total sum of weights for each radius
    #plt.imshow(infiles[1].T*rad_mask)
    #plt.savefig(img_dir + 'coadded.png')
    weight_sum = [sum(x) for x in zip(*weight_bins)] #need * to SPLAT list of weights :)
    #Apply weighted sum for each radial bin

    for rad_ct in range(len(radius_bins)):
        radius_weights = [weight_bin[rad_ct]/weight_sum[rad_ct] for weight_bin in weight_bins]
        for im_ct in range(len(masked_images)):
            weighted_img_at_rad = masked_images[im_ct]*radius_weights[im_ct] #Weighted image for a given radius
            final_img += weighted_img_at_rad#weighted_img_at_rad.T
            #final_img = np.add(final_img,weighted_img_at_rad.T,out=final_img,casting='unsafe') #now add to final image
    plt.imshow(final_img.T)
    plt.colorbar()
    plt.savefig(img_dir+'coadded.png')
    #Save as fits file
    fits_temp = fits.open(fits_file)
    fits_temp[0].data = final_img.T
    fits_temp.writeto(img_dir+'GGF.fits',overwrite=True)
    return None

def make_radial_mask(img_array,radii,center_pixel):
    '''
    Create mask for image based off weight bins, center pixel, and radial values
    :param img_array: numpy array from reading in image
    :param radii: (R_in,R_out) tuple of inner and outer radius WITHIN mask
    :param center_pixel: (X,Y) tuple of central pixel from which the radial bins expand
    :return: radial mask for image as np array of booleans
    '''
    #Step through entire image array
    (len_x,len_y) = img_array.shape
    mask = np.empty((len_y,len_x))
    for i in range(len_x):
        i_pos = i - center_pixel[0]
        for j in range(len_y):
            j_pos = j - center_pixel[1]
            if ((i_pos**2+j_pos**2>=radii[0]) and (i_pos**2+j_pos**2<=radii[1])):
                mask[j,i] = True #Pixel within radii bounds
            else:
                mask[j,i] = False #Pixel outside of radii bounds
    return mask

def main():
    img_dir = '/home/carterrhea/Documents/M87/X-ray/'
    infile = '/home/carterrhea/Documents/M87/X-ray/broad_f_unbinned.img'
    rcParams['figure.figsize'] = 32,16
    filtered = []
    radii_ = np.array([0,250,1000,2500,4000])/0.492
    center_ = (1394,1451) #broad_f_unbinned
    #center_ = (836,832) #broad_un_cen.img
    for sigma in [1,2,4,6]:
        outfile = '/home/carterrhea/Documents/M87/X-ray/M87_sig'+str(sigma)
        filtered.append(ggf1(infile,outfile,sigma))
    combine_ggf(img_dir,filtered, radii_, [[1,1,1,0,0],[0,4,4,4,0],[0,0,8,8,8],[0,0,0,32,32]], center_, infile)
main()
