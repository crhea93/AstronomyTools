'''
File to create the correctly binned pi files
This is only important if you are running X-ray data analysis from Chandra
---------------------------------------------------
Goal: Create binned spectra from Chandra data given
    the WVT of the pixels
---------------------------------------------------
INPUTS:
    filename - WVT output (e.g.'/home/user/Desktop/WVT_data.txt')
    base_directory - Directory with Chandra data (e.g.'/home/usr/CHANDRA')
	source_file - File to read in and used to create bins in WVT (e.g.'source')
    background - name of background file (e.g. 'background_simple')
        Set to 'blank' if using a blank-sky background file.
        Have the file named "obsid_blank.evt"
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
import time
from ciao_contrib.runtool import *
from crates_contrib.utils import *
import shutil
#-------------------------------------------------#
#-------------------------------------------------#
# Bin Class Information
# Everything in here should be self-explanatory... if not let me know and I
# will most happily comment it! :)
class Bin:
    def __init__(self, number):
        self.bin_number = number
        self.pixels = []
        self.total_pixels = 0
    def add_pixel(self, Pixel):
        self.pixels.append(Pixel)
        self.total_pixels += 1
#-------------------------------------------------#
#-------------------------------------------------#
# Pixel Class Information
# Ditto à propos the documentation for this class
class Pixel:
    def __init__(self, number, pix_x, pix_y, width, height):
        self.pix_number = number
        self.pix_x = pix_x
        self.pix_y = pix_y
        self.width = width
        self.height = height
    def update(self, pix_x, pix_y, width, height):
        self.pix_x = pix_x
        self.pix_y = pix_y
        self.width = width
        self.height = height
#-------------------------------------------------#
#-------------------------------------------------#
# Get necessary filenames such as evt2,asol1,bpix1, and msk1
#   parameters:
#
def get_filenames(dir):
    filenames = dict()
    for file in os.listdir(dir):
        if file.endswith("evt2.fits"):
            filenames['evt2'] = file
        if file.endswith("bpix1.fits"):
            filenames['bpix1'] = file
    return filenames

def update_header(ann,ann_values):
    '''
    Must update header for DSDEPROJ
    '''
    dmhedit.punlearn()
    dmhedit.infile = 'Annuli/Annulus_'+str(ann)+'.pi'
    dmhedit.filelist = None
    dmhedit.operation = "add"
    dmhedit.key = "XFLT0001"
    dmhedit.value = ann_values[ann-1]
    dmhedit()
    return None

#---------------------------------------------------#
#---------------------------------------------------#
# specextract_run
# apply specextract functio from ciao
#   parameters:
#       filenames = list of files necessary for specextract
#       file_to_convert = fits file in string format
#       outfile_from_convert = pha outroot in string format
#       output_dir = directory for output file
def specextract_run(obsid,filenames,file_to_convert,outfile_from_convert,output_dir):
    specextract.punlearn()
    specextract.infile = file_to_convert+"[sky=region("+str(output_dir)+"temp.reg)]" #1 to make sure we have enough space. serious overkill!
    specextract.outroot = outfile_from_convert
    specextract.bkgfile = obsid+'_blank.evt[sky=region('+str(output_dir)+'temp.reg)]'
    specextract.bkgresp = False #Necessary if using blank sky
    specextract.badpixfile = filenames['bpix1']
    specextract.grouptype = 'NUM_CTS'
    specextract.binspec = 1
    specextract.clobber = True
    specextract.energy_wmap = '500:14000'
    print("     Running Specextract...")
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
        file.write("# Region file format: DS9 version 4.1 \n")
        file.write("global color=green dashlist=8 3 width=1 font='helvetica 10 normal' select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1 \n")
        file.write("physical \n")
        count = 0
        for region in regions:
            if count == 0:
                file.write(region+" \n")
            else:
                file.write('+'+region+" \n")
            count += 1
#---------------------------------------------------#
#---------------------------------------------------#
#create event file with regions in bin
#   parameters:
#
def create_evt(file_to_split,bin_number,output_dir):
    dmcopy.punlearn()
    dmcopy.infile = file_to_split+'.fits[sky=region('+output_dir+'temp.reg)]'
    dmcopy.outfile = output_dir+'bin_'+str(bin_number)+'.fits'
    dmcopy.clobber = True
    dmcopy()
#---------------------------------------------------#
#---------------------------------------------------#
# Combine individuals pixels into larger boxes to
# reduce spexectract time.
def create_reg_comb(pix_in_bin,file_to_split,bin_number,output_dir):
    points = {} #(x,[y]) pairs
    for pixel in pix_in_bin:
        if pixel.pix_x not in points.keys():
            points[pixel.pix_x] = [pixel.pix_y]
        else:
            points[pixel.pix_x].append(pixel.pix_y)
    pixel_x_list = [pixel.pix_x for pixel in pix_in_bin]
    pixel_y_list = [pixel.pix_y for pixel in pix_in_bin]
    x_min = np.min(pixel_x_list)
    y_min = np.min(pixel_y_list)
    x_max = np.max(pixel_x_list)
    y_max = np.max(pixel_y_list)
    pixels_used = {}

    for x in range(int(x_min),int(x_max)+1):
        if x not in pixels_used.keys():
            pixels_used[x] = []
    with open(output_dir+'temp.reg','w+') as file:#will create a temporary file used to create spectra
        file.write("# Region file format: CIAO version 1.0 \n")
        file_phys = open(output_dir+'temp_phys.reg','w+')
        file_phys.write("# Region file format: CIAO version 1.0 \n")
        #Step through each value to build up rectangle!
        in_count = 0
        for y in range(int(y_min),int(y_max)+1):
            for x in range(int(x_min),int(x_max)+1):
                xwi = 1 #Let us try and add a pixel to our box
                ywi = 1
                if y in points[x]: #only check point if it is in our bin
                    new_box = False #will we create a new region box?
                    continue_extension = True #we want to initialize the while loop for the box
                    if y in pixels_used[x]: #if we already assigned the pixel then dont contine
                        continue_extension = False
                        new_box = False
                    if y not in pixels_used[x]: #if not assigned yet we should assign
                        pixels_used[x].append(y)
                        new_box = True #that means we need a new box!
                    in_count += 1
                    while continue_extension == True: #Lets see how large we can make this box
                        if x+xwi > x_max: #Are we within bounds of bin?
                            break
                        if y+ywi > y_max: #ditto wrt y...
                            break
                        yextend = True
                        for xi in range(x,x+xwi+1): #can we add a y for all x-values currently used
                            if y+ywi in pixels_used[xi] or y+ywi not in points[xi]:
                                yextend = False
                                break
                        xextend = True
                        for yi in range(y,y+ywi+1): #ditto for x
                            if yi not in points[x+xwi] or yi in pixels_used[x+xwi]:
                                xextend = False
                                break
                        xandyextend = False #ok can we add the corner xandy extension point?
                        if x+xwi <= x_max: #lets make sure there is no conflicts
                            xandyextend = (y+ywi not in pixels_used[x+xwi] and y+ywi<=y_max)
                        if x+xwi > x_max:
                            pass
                        if (yextend and xextend and xandyextend): #ok were all good to extend our box!
                            for xi in range(x,x+xwi+1): #gotta get all those new x-values
                                if xi<=x_max:
                                    for yi in range(y,y+ywi+1): #and dont forget those y-values :)
                                        if yi not in pixels_used[xi]:
                                            pixels_used[xi].append(yi) #if it hasn't already been noted lets add it
                            xwi += 1 #ok now to make our box potentially bigger....
                            ywi += 1
                        if not (yextend and xextend and xandyextend):
                            continue_extension = False #fail!!!!
                    if new_box == True:#now we make a new box!
                        x_coord1,y_coord1,ra_center,dec_center = coord_trans((x+.5)+xwi*.5,(y+.5)+ywi*.5,file_to_split)
                        file.write("box(%.16f,%.16f,%i,%i) \n"%(x_coord1,y_coord1,xwi,ywi)) #Have to subtract 0.25 from x to make it align -- non permanent
                        file_phys.write("box(%s,%s,%i,%i) \n"%(ra_center,dec_center,xwi,ywi))
                else:
                    #value not in list of pixels
                    pass
        file_phys.close()
    shutil.copy(output_dir+'temp_phys.reg',output_dir+str(bin_number)+'.reg')
    return None
#---------------------------------------------------#
#---------------------------------------------------#
# Split up fits files into pi/pha files
#   parameters:
#       filenames = list of files necessary for specextract
#       file_to_split = fits file in string format
#       output_file = directory for output
#       output_split = pha outroot in string format
#       x_center = pixel x-position center in physical (sky) coordinates
#       y_center = pixel y-position center in physical (sky) coordinates
#       width = width of pixel in sky coordinates
#       height = height of pixel in sky coordinates
#       pix_in_bin_num = pixel number relative to bin (0-max(bin.pixels))
#       bin_number = As suspected its the bin number :):):)
def split_fits(obsid,filenames,file_to_split,output_file,output_dir,pix_in_bin_num,bin_number):
    output =  output_file+"_"+str(bin_number)+"_"
    regions = []
    file_to_convert = file_to_split+'.fits'
    #create_reg_comb(pix_in_bin_num,file_to_split,bin_number,output_dir)
    create_evt(file_to_split,bin_number,output_dir)
    specextract_run(obsid,filenames,file_to_convert,output,output_dir)
    update_header()
    return None
#---------------------------------------------------#
#---------------------------------------------------#
#Change from image or logical coordinate system to physical or sky coordinate system
#   parameters:
#       pixel_x  = pixels x coordinate in image coordinates
#       pixel_y  = pixels y coordinate in image coordinates
#       file_to_split = file which we will split later (needed to get new coordinates)
def coord_trans(pixel_x,pixel_y,file_to_split):
    tr = SimpleCoordTransform(file_to_split+'.img')
    x_center,y_center = tr.convert("image", "physical", pixel_x, pixel_y)
    dmcoords.punlearn()
    dmcoords.infile = file_to_split+'.img'  # OBSID+'_broad_thresh.img'
    dmcoords.option = 'logical'
    dmcoords.logicalx = pixel_x
    dmcoords.logicaly = pixel_y
    dmcoords()
    ra = dmcoords.ra
    dec = dmcoords.dec
    return x_center,y_center,ra,dec
#---------------------------------------------------#
#---------------------------------------------------#
# Main program to create bin pi/pha files
#   parameters:
#       base_directory = Directory containing Chandra data
#       filename = name of file to read in WVT bin data
#       dir = Directory for Chandra OBSID
#       file_to_split = Name of file to split in repro directory
#       output_dir = Output path for binned files
def create_spectra(base_directory,filename,OBSIDS,source_file,output_dir,wvt_output):
    print("Starting to bin spectra...")
    for obsid in OBSIDS:
        print("We are on obsid %s"%obsid)
        print("#-------------------------------------------------------------------------------------#")
        directory = base_directory+'/'+obsid+'/repro'
        #Just making sure the directory exists and is empty :)
        if not os.path.exists(directory+"/"+output_dir):
            os.makedirs(directory+"/"+output_dir)
        if len(os.listdir(directory+"/"+output_dir)) != 0:
            print("Cleaning output directory of  files")
            for item in os.listdir(directory+"/"+output_dir):
                os.remove(os.path.join(directory+"/"+output_dir, item))
        os.chdir(directory)
        output_file = directory+'/'+output_dir+'/'+source_file
        file_to_split = directory+'/'+source_file
        #snagging some file names for later...
        filenames = get_filenames(directory)
        bins = []
        number_bins = -1
        pix_num = 0
        #Lets get the data from our WVT file
        with open(base_directory+'/'+filename+'.txt') as f:
            next(f)
            next(f) #skip first two lines
            for line in f:
                if not int(line.split(" ")[2]) in [bins[i].bin_number for i in range(len(bins))]:
                    # create new Bin instance if doesnt exist already
                    new_bin = Bin(int(line.split(" ")[2]))
                    bins.append(new_bin)
                    number_bins += 1
                # Create new pixel
                new_pix = Pixel(pix_num,float(line.split(" ")[0]),float(line.split(" ")[1]),0,0) #set width and height to zero
                pix_num += 1
                # Add pixel to current bin
                bins[number_bins].add_pixel(new_pix)

        #Get unique bin numbers and order bins ascending
        for bin_i in bins:#bin_unique: #for each unique_bin
            pix_in_bin = bin_i.pixels
            spec_to_combine = []
            print("  We are combining bin number "+str(bin_i.bin_number+1)+" of "+str(len(bins)))
            print("     We have %i pixels"%bin_i.total_pixels)
            start = time.time() #lets do the actual splitting algorithm
            try:
                split_fits(obsid,filenames,file_to_split,output_file,output_dir,pix_in_bin,bin_i.bin_number)
            except:
                print("     Not enough counts in the region to create a spectrum")
            print("     The creation of the spectrum took %.5f seconds"%(time.time()-start))
            print()
            for item in os.listdir(directory+'/'+output_dir):#we cab get rid of temporary files
                if item.endswith(("_temp.reg")):
                    os.remove(os.path.join(directory+'/'+output_dir, item))
    return len(bins)
