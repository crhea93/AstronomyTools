'''
File to create the correctly binned pi files
This is only important if you are running X-ray data analysis from Chandra
---------------------------------------------------
Goal: Create binned spectra from Chandra data given
    the WVT of the pixels
---------------------------------------------------
INPUTS:
    filename - WVT output (e.g.'/home/user/Desktop/WVT_data.txt')
    dir - Directory with Chandra data (e.g.'/home/usr/CHANDRA/12833/repro/'
    file_to_split - File to read in and used to create bins in WVT (e.g.'imageA')
    output_dir - Output directory concatenated with dir  (e.g.'binned/')
---------------------------------------------------
List of Functions (In order of appearance):
    specextract_run --> Creates PI file for individually binned images
    split_fits --> Break input fits by WVT bins into single pixel pha files
    combine_pha --> Combine the single pixel phas into combined spectra based off WVT bins
    create_spectra --> Main function to create binned spectra
---------------------------------------------------
OUTPUTS:
    -- A combined spectra for each bin as designated by the WVT.
    -- This is to be used for spectral fitting (we'll that's why I made this program)
    -- File put in /PathToChandraData/OBSID/repro/binned
---------------------------------------------------
Additional Notes:
    As mentioned, the program was designed to generate combinned-binned-spectra
    so that I could generate temperature maps...
    The program can easily be canabilized for other uses or specifications
---------------------------------------------------
---------------------------------------------------
Carter Rhea
carterrhea93@gmail.com
https://carterhea93.com
'''
import numpy as np
import os
from ciao_contrib.runtool import *
from crates_contrib.utils import *
# Get necessary filenames such as evt2,asol1,bpix1, and msk1
#   parameters:
#
def get_filenames():
    filenames = dict()
    for file in os.listdir(os.getcwd()+'/repro'):
        if file.endswith("repro_evt2.fits"):
            filenames['evt2'] = file
        if file.endswith("repro_bpix1.fits"):
            filenames['bpix1'] = file
    for file in os.listdir(os.getcwd()+'/primary'):
        if file.endswith("_asol1.fits"):
            filenames['asol1'] = file
    for file in os.listdir(os.getcwd()+'/secondary'):
        if file.endswith("_msk1.fits"):
            filenames['msk1'] = file
    os.chdir('repro')
    return filenames
#---------------------------------------------------#
#---------------------------------------------------#
# specextract_run
# apply specextract functio from ciao
#   parameters:
#       filenames = list of files necessary for specextract
#       file_to_convert = fits file in string format
#       outfile_from_convert = pha outroot in string format
#       output_dir = directory for output file
def specextract_run(filenames,file_to_convert,outfile_from_convert,output_dir):
    specextract.infile = file_to_convert+"[sky=region("+str(output_dir)+"temp.reg)]" #1 to make sure we have enough space. serious overkill!
    specextract.outroot = outfile_from_convert
    specextract.bkgfile = filenames['evt2'] +'[sky=region(simple_bkg.reg)]'
    specextract.asp = '../primary/'+filenames['asol1']
    specextract.mskfile = '../secondary/'+filenames['msk1']
    specextract.badpixfile = filenames['bpix1']
    specextract.grouptype = 'NUM_CTS'
    specextract.binspec = 1
    specextract.clobber = True
    specextract.energy_wmap = '100:14000'
    specextract()
    return True
#---------------------------------------------------#
#---------------------------------------------------#
#create temporary region file with regions in bin
#   parameters:
#       output_dir = output directory
#       reigons = list of regions to add
def create_reg(output_dir,regions):
    #Create temporary reg file for split
    with open(str(output_dir)+"temp.reg","w+") as file:
        file.write("# Region file format: CIAO version 1.0 \n")
        for region in regions:
            file.write(region+" \n")
#---------------------------------------------------#
#---------------------------------------------------#
# Split up fits files into pi/pha files
#   parameters:
#       filenames = list of files necessary for specextract
#       file_to_split = fits file in string format
#       output_dir = directory for output
#       output_split = pha outroot in string format
#       x_center = pixel x-position center in physical (sky) coordinates
#       y_center = pixel y-position center in physical (sky) coordinates
#       width = width of pixel in sky coordinates
#       height = height of pixel in sky coordinates
#       pix_in_bin_num = pixel number relative to bin (0-max(bin.pixels))
#       bin_number = As suspected its the bin number :):):)
def split_fits(filenames,file_to_split,output_dir,x_center,y_center,width,height,pix_in_bin_num,bin_number):
    output=  output_dir+file_to_split+"_"+str(bin_number)+"_"
    regions = []
    for pix in range(pix_in_bin_num):
        regions.append("box("+str(x_center[pix])+","+str(y_center[pix])+","+str(width[pix])+","+str(height[pix])+")")
    file_to_convert = file_to_split+'.fits'
    create_reg(output_dir,regions)
    specRun = specextract_run(filenames,file_to_convert,output,output_dir)
    return None
#---------------------------------------------------#
#---------------------------------------------------#
#Change from image or logical coordinate system to physical or sky coordinate system
#   parameters:
#       pixel_x  = pixels x coordinate in image coordinates
#       pixel_y  = pixels y coordinate in image coordinates
#       file_to_split = file which we will split later (needed to get new coordinates)
def coord_trans(pixel_x,pixel_y,file_to_split):
    pixel_x = pixel_x+0.5
    pixel_y = pixel_y+0.5
    tr = SimpleCoordTransform(file_to_split+"_image.fits")
    x_center,y_center = tr.convert("image", "physical", pixel_x, pixel_y)
    x_min,y_min = tr.convert("image", "physical", pixel_x-0.5, pixel_y-0.5)
    x_max,y_max = tr.convert("image", "physical", pixel_x+0.5, pixel_y+0.5)
    width = x_max-x_min
    height = y_max-y_min
    return x_center,y_center,width,height
#---------------------------------------------------#
#---------------------------------------------------#
# Main program to create bin pi/pha files
#   parameters:
#       filename = name of file to read in WVT bin data
#       dir = Directory for Chandra OBSID
#       file_to_split = Name of file to split in repro directory
#       output_dir = Output path for binned files
def create_spectra(filename,dir,file_to_split,output_dir):
    os.chdir(dir)
    if not os.path.exists(dir+"/repro/"+output_dir):
        os.makedirs(dir+"/repro/"+output_dir)
    if len(os.listdir(dir+"/repro/"+output_dir)) != 0:
        print("Cleaning output directory of  files")
        for item in os.listdir(dir+"/repro/"+output_dir):
            os.remove(os.path.join(dir+"/repro/"+output_dir, item))
    filenames = get_filenames()
    pixel_x = []
    pixel_y = []
    bin = []
    with open(filename) as f:
        next(f)
        next(f) #skip first two lines
        for line in f:
            pixel_x.append(float(line.split(" ")[0]))
            pixel_y.append(float(line.split(" ")[1]))
            bin.append(int(line.split(" ")[2]))
    #Get unique bin numbers and order bins ascending
    bin_sorted = sorted(bin)
    bin_unique = set(bin_sorted)
    total_num = 0
    x_pixels = []
    y_pixels = []
    widths = []
    heights = []
    for bin_i in bin_unique: #for each unique_bin
        pix_in_bin_num = 0
        spec_to_combine = []
        print("We are combining bin number "+str(bin_i+1)+" of "+str(len(bin_unique)))
        for chip_set_num in range(bin.count(bin_i)): #for number of said unique bin
            x_center,y_center,width,height = coord_trans(pixel_x[total_num],pixel_y[total_num],file_to_split)
            x_pixels.append(x_center)
            y_pixels.append(y_center)
            widths.append(width)
            heights.append(height)
            pix_in_bin_num += 1
        specRun = split_fits(filenames,file_to_split,output_dir,x_pixels,y_pixels,widths,heights,pix_in_bin_num,bin_i)
        for item in os.listdir(dir+'/repro/'+output_dir):
            if item.endswith(("_temp.reg")):
                os.remove(os.path.join(dir+'/repro/'+output_dir, item))
def main():
    filename = '/home/crhea/Documents/TestData/WVT_data.txt'
    dir = '/home/crhea/Documents/TestData/12036'
    file_to_split = 'simple'
    output_dir = 'binned/'
    create_spectra(filename,dir,file_to_split,output_dir)

main()
