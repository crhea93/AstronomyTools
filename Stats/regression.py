'''
A simple python script for regression fitting data from a mysql database
'''

##---------------------------IMPORTS----------------------------------------##
import pandas as pd
import numpy as np
from scipy import odr
import mysql.connector
import pandas.io.sql as psql
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
##--------------------------Plotting Database-------------------------------##
plot_names = {'csb_ct':'Coefficient of Surface Brightness (ct/s)',
              'csb_pho':'Coefficient of Surface Brightness (photons/s)',
              'csb_flux':'Coefficient of Surface Brightness (erg/s)',
              'r_cool_3':'Cooling Radius at 3 Gyr (kpc)',
              'r_cool_7':'Cooling Radius at 7 Gyr (kpc)'}
##--------------------------DATABASE----------------------------------------##
def read_database(password,table_name):
    mydb = mysql.connector.connect(
        host="localhost",
        user="carterrhea",
        passwd=password,
        database="Lemur_DB"
    )
    query = "SELECT * FROM " + table_name
    # execute the query and assign it to a pandas dataframe
    df = psql.read_sql(query, con=mydb)
    # close the database connection
    mydb.close()
    return df



def read_password(input_file):
    '''
    read in password for database
    PARAMETERS:
        input_file - location of file with password
    '''
    pword = ''
    with open(input_file) as f:
        #Read file
        for line in f:
            if 'password' in line: #Only read lines with '='
                pword= line.split("=")[1].strip()
            else: pass

    return pword

def line(x,a,b):
    return a*x + b

def f(B, x):
    '''Linear function y = m*x + b'''
    # B is a vector of the parameters.
    # x is an array of the current x values.
    # x is in the same format as the x passed to Data or RealData.
    #
    # Return an array in the same format as y passed to Data or RealData.
    return B[0]*x + B[1]
def plot_lin_reg(X, Y, x_name, y_name, x_err=[], y_err=[]):
    '''
    plot linear regression results
    :param x: data for x axis
    :param y: data for y axis
    :param intercept: returned intercept value from linear regression fit
    :param coeff: returned coefficient value from linear regression fit
    :return: None
    '''
    X = X.to_numpy(); Y = Y.to_numpy()
    #popt, pcov = curve_fit(line, X, Y, xerr=x_err, yerr=y_err)
    linear = odr.Model(f)
    mydata = odr.RealData(X, Y, sx=x_err, sy=y_err)
    myodr = odr.ODR(mydata, linear, beta0=[1., 0.2])
    myoutput = myodr.run()
    linear_fit_x = np.linspace(np.min(X), 1.1*np.max(X),100)
    linear_fit_y = list(map(lambda x: x*myoutput.beta[0] + myoutput.beta[1] , linear_fit_x))
    #plotting
    fig,ax = plt.subplots()
    ax.plot(linear_fit_x,linear_fit_y)
    ax.errorbar(X, Y, xerr=x_err, yerr=y_err, color='coral',  fmt='o', capsize=3, capthick=1, ecolor='coral')
    plt.xlabel(plot_names[x_name])
    plt.ylabel(plot_names[y_name])
    fig.show()

def artificial_errors(err,tol):
    '''
    Add minor value to error is none is included
    :param err:  initial error values
    :param tol: artificial error
    :return: new error values
    '''
    vals = []
    for i in range(len(err)):
        if err[i] == 0.0:
            vals.append(tol)
        else:
            vals.append(err[i])
    return vals

def main():
    x_name = 'csb_ct'
    y_name = 'csb_flux'

    #Set up database info
    cols_to_drop = ['csb_ct','csb_flux','csb_pho','R_cool_3','R_cool_7']
    password = read_password('/home/carterrhea/Desktop/pword_mysql.txt')
    clusters = read_database(password, 'Clusters')
    # remove duplicates from clusters
    clusters = clusters.drop(columns=cols_to_drop)
    csb_data = read_database(password, 'csb')
    rcool_data = read_database(password, 'r_cool')
    clusters = pd.merge(clusters, csb_data, on='Name', how='outer', copy=False)
    clusters = pd.merge(clusters, rcool_data, on='Name', how='outer', copy=False)
    #Set up values and errors
    x = clusters[x_name]
    y = clusters[y_name]
    xerr = artificial_errors(clusters[x_name+'_u'] - x, 1e-8)
    yerr = artificial_errors(clusters[y_name+'_u'] - y, 1e-8)
    #plot
    plot_lin_reg(x, y, x_name, y_name, xerr, yerr)

main()