'''
Obtain gmail account amd password from secret folder containing such information
parameters:
    gmail_info_file -- full path to gmail info file (e.g. '/home/user/Documents/secretFolder/gmail_info.txt')
        This file needs to be 2 lines only of the following format:
        gmail_account = account_name
        gmail_password = password
outputs:
    gmail account name (string)
    gmail password (string)
'''
def read_gmail_info(gmail_info_file):
    with open(gmail_info_file) as f:
        contents = []
        for line in f:
            contents.append(line)
    gmail_account = contents[0].split("=")[1].strip()
    gmail_password = contents[1].split("=")[1].strip()
    return gmail_account, gmail_password
#-------------------------------------------------#
#-------------------------------------------------#