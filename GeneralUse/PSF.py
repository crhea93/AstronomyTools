'''
Create PSF

Gotta run marx  bash script first

First must approximate Photon Flux using ideal PSF -- needed for ChaRT input...


INPUTS:
    chandra_dir -- full path to data directory (e.g. '/home/user/Documents/Data')
    evt_file -- name of event file without extension (e.g. 'acisf#####_repro_evt2')
    region -- name of region file of interest without .reg extension (e.g. 'simple')
    background -- name of background region file without .reg extension (e.g. 'simple_background')
    exposure -- Boolean determining method to calculate Net Energy Flux. See
        Documentation for more information. (e.g. True)
    source_ra -- ra of source (e.g. "00:00:00.0")
    source_dec -- dec of source (e.g. '+00:00:00.0')
    energy_range -- energy range in electron volts (e.g. '500:2000')
    OBSID -- Chandra Obsid (e.g. '#####')
    gmail_info_file -- full path to gmail info file (e.g. '/home/user/Documents/secretFolder/gmail_info.txt')
        This file needs to be 2 lines only of the following format:
        gmail_account = 'account_name'
        gmail_password = 'password'
'''

import os
import wget
import time
import email
import imaplib
import tarfile
import selenium
from selenium import webdriver
from ciao_contrib.runtool import *
from ToolBox import calc_effenergy,calc_flux,calc_bounds,simulatePSF,read_gmail_info

#---------------------INPUTS---------------------------------------------#
chandra_dir = '%%%'
evt_file = '%%%'
region = '%%%'
background = '%%%'
exposure = False
source_ra = '%%%'
source_dec = '%%%'
energy_range = '%%%' #in electron volts
gmail_info_file = '%%%'
#------------------------------------------------------------------------#



def prepare(evt_file,region,background,exposure,energy_range,source_ra,source_dec):
    energies = [float(x) for x in energy_range.split(':')]
    energy_range2 = str(energies[0]/1000)+':'+str(energies[1]/1000) #for effective energy (eV)
    #Calculate Effective Energy
    effen = calc_effenergy(region,energy_range2)
    #Calculation of Source Flux
    print("Calculating Net Energy Flux and Obtaining Bounds")
    calc_flux(evt_file,energy_range,region,background,exposure)
    #Now read flux
    energy_flux,energy_flux_lower,energy_flux_upper = calc_bounds(region,'NEFB')
    return effen,energy_flux

def submit_form(source_ra,source_dec,eff_en,energy_flux,gmail_account,OBS_ID):
    url = 'http://cxc.harvard.edu/ciao/PSFs/chart2/runchart.html'
    data_dict = {}
    data_dict['email'] = gmail_account
    data_dict['ra'] = source_ra
    data_dict['dec'] = source_dec
    data_dict['mono'] = eff_en
    data_dict['flux'] = energy_flux
    data_dict['obsid'] = OBS_ID
    #data = urlencode(data_dict).encode("utf-8")
    #r = requests.post(url, data=data_dict)
    driver = webdriver.Chrome('/usr/lib/chromium-browser/chromedriver')
    # Open the website
    driver.get(url)
    #Fill in form
    for key,val in data_dict.items():
        id_box = driver.find_element_by_name(key)
        id_box.send_keys(str(val))
    # Find login button
    login_button = driver.find_element_by_name('submit')
    # Click login
    login_button.click()
    driver.close()
    return None

def download_data(gmail_account,gmail_password):
    recieved = False
    mail = imaplib.IMAP4_SSL('imap.gmail.com')
    mail.login(gmail_account,gmail_password)
    mail.select('Inbox')
    #Search for specific email
    typ,data = mail.search(None, '(SUBJECT "ChaRT rays ready for download")')
    ids = data[0] # data is a list.
    print()
    id_list = ids.split() # ids is a space separated string
    if len(id_list) == 0:
        print("Ray File Not yet recieved")
    else:
        print("Recieved!")
        recieved = True
        latest_email_id = id_list[-1] # get the latest

        result, data = mail.fetch(latest_email_id, "(RFC822)") # fetch the email body (RFC822)             for the given ID

        raw_email = data[0][1]
        #continue inside the same for loop as above
        raw_email_string = raw_email.decode('utf-8')
        # converts byte literal to string removing b''
        email_message = email.message_from_string(raw_email_string)
        # this will loop through all the available multiparts in mail
        for part in email_message.walk():
            if part.get_content_type() == "text/plain": # ignore attachments/html
                body = part.get_payload(decode=True)
                body = body.decode('utf-8')
                body = body.split("\n")
                weburl = body[4].split("\t")
                weburl = weburl[0].strip()
        #Download tar file
        wget.download(weburl)
        filename_wget = weburl.split("/")[-1]
        #untar file
        tar = tarfile.open(filename_wget)
        tar.extractall()
        tar.close()
        #remove tar file
        os.remove(filename_wget)
    return recieved



def main():
    gmail_account,gmail_password = read_gmail_info(gmail_info_file)
    prepare_bool = input("Do you need to prepare the Chart run: ") #Yes or No
    os.chdir(chandra_dir)
    if prepare_bool.lower() == 'yes':
        eff_en,energy_flux = prepare(evt_file,region,background,exposure,energy_range,source_ra,source_dec)
    if prepare_bool.lower() != 'yes':
        eff_en = float(input("Please input Effective Energy: "))
        energy_flux = float(input("Please input Net Energy Flux: "))
    if energy_flux < 1E-9:
        print("The Net Energy Flux is too low so we cannot proceed with using ChaRT")
        print("The programming is exiting")
        exit()
    print("Submitting Form")
    submit_form(source_ra,source_dec,eff_en,energy_flux,gmail_account,OBS_ID)
    recieved = False
    print("Waiting to recieve email")
    while recieved == False:
        time.sleep(10)
        recieved = download_data(gmail_account,gmail_password)
    for file in os.listdir(os.getcwd()):
        if file.endswith("rays.fits"):
            psf_file = file
    print("Simulating PSF")
    simulatePSF(evt_file,psf_file)

main()
