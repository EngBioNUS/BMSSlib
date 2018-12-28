# -*- coding: utf-8 -*-
"""
Created on Fri Dec 28 17:02:29 2018

@author: jingwui
"""

#!/usr/bin/env python

### This module to read raw data, convert units, and plot the processed data  ###

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def readData(filename, NumDataSet):
    ### Read CSV into raw_data dataframe ###
    input_data = pd.read_csv(filename)
    raw_data_header = input_data.columns[1:NumDataSet+1].tolist()
    

    ### Print Data headers
    print('Input Data File Header:', raw_data_header)

    ### number of columns and rows in raw_data (shape only for dataframe)
    input_data_numrows = input_data.shape[0]
    
    input_data_numcols = input_data.shape[1]- NumDataSet #minus off the std deviation data
    print('Input Data: Rows -', input_data_numrows, ', Columns -', input_data_numcols)

    ### Saving raw_data columns into individual list (rows) in rfp_data numpy
    # Flip row and col from raw_data
    rfp_data = np.zeros((input_data_numcols,input_data_numrows))

    for i in range(0, input_data_numcols): # Iterate through all columns
        rfp_data[i, :] = input_data.iloc[:, i]
    # print(rfp_data)

    # Time + All Inducers
    rfp_data_numrows = rfp_data.shape[0] # or = np.size(rfp_data,0)
    # Number of RFP Data Per Inducer
    rfp_data_numcols = rfp_data.shape[1]

    ### Convert RFP/OD into (mol/OD) ###
    # Starts from 1 because Row 0 is time
    for i in range (1, rfp_data_numrows):
        for j in range (0, rfp_data_numcols):
            rfp_data[i][j] = ODconv(rfp_data[i][j])

    #store standard deviation data
    stddev_data = np.zeros((NumDataSet, input_data_numrows))

    for i in range(input_data_numcols, input_data_numcols + NumDataSet): # Iterate through all columns
        stddev_data[i-input_data_numcols, :] = input_data.iloc[:, i]

    stddev_data_numrows = stddev_data.shape[0]
    stddev_data_numcols = stddev_data.shape[1]

    # convert RFP/OD into Molar/OD
    for i in range(0, stddev_data_numrows):
        for j in range(0, stddev_data_numcols):
            stddev_data[i][j] = ODconv(stddev_data[i][j])

    return raw_data_header, rfp_data, stddev_data

### Convert RFPOD into (Molar/OD) ###
def ODconv(x):
    x_conv = x *(1 / 18.814) * (1 / 30) * (10 ** -6)
    return x_conv

def indconv_Percent(raw_data_header, Molar_Mass):
    ### Convert Inducer from % into (mol/L) ###
    numInd = np.size(raw_data_header)
    inducer_array = np.zeros(numInd)

    Inducer_MM = Molar_Mass

    #Arabinose_MM = 150.13  # g/mol

    # Starts from 1 because Row 0 is time
    inducer_array = raw_data_header[0:numInd]
    print(inducer_array)
    for i in range (0,numInd):
        #read only the digit (exclude the unit)
        inducer_array[i] = ''.join([c for c in inducer_array[i] if c in '1234567890.'])
        inducer_array[i] = float(inducer_array[i])
        inducer_array[i] = ((inducer_array[i] / 0.1) / Inducer_MM)
    return inducer_array

def indconv(raw_data_header, ind_unit):
    ### Convert Inducer from % into (mol/L) ###
    numInd = np.size(raw_data_header)
    inducer_array = np.zeros(numInd)

    # Starts from 1 because Row 0 is time
    inducer_array = raw_data_header[0:numInd]
    for i in range (0,numInd):
        inducer_array[i] = ''.join([c for c in inducer_array[i] if c in '1234567890.'])
        inducer_array[i] = float(inducer_array[i])
        inducer_array[i] = inducer_array[i] * ind_unit
    return inducer_array

def indconv_None(raw_data_header):

    numInd = np.size(raw_data_header)
    inducer_array = np.zeros(numInd)

    # Starts from 1 because Row 0 is time
    inducer_array = raw_data_header[0:numInd]
    for i in range (0,numInd):
        inducer_array[i] = ''.join([c for c in inducer_array[i] if c in '1234567890.'])
        inducer_array[i] = float(inducer_array[i])
    return inducer_array

### Plot CSV Data ###
def plot_inputdata(input_data_header, rfp_data, stddev_data):
    ### Defining time array from Row 0 of rfp_data
    time = rfp_data[0]

    # Time grid == no. of times to report solution
    rfp_data_numrows = np.size(rfp_data, 0)
    # rfp_data_numcols = np.size(rfp_data, 1)

    fig = plt.figure(figsize=(5,3.6))
    ax = fig.add_axes([0.16,0.15,0.8,0.78])
    ### Plot RFP Data vs Time ###
    # plt.figure()
    style = ['^','*','>','D','<','s','p','o','d','+','h','x','v','.','H']
    for i in range (1, rfp_data_numrows):
        plt.errorbar(time, rfp_data[i], yerr = stddev_data[i-1], linestyle='None', marker = style[i-1], capsize=2,  markersize = 4) 

    ### Figure Settings
    #plt.rc('font', size=12)  # controls default text sizes
    rfp_data_legend = input_data_header[0:rfp_data_numrows]
    plt.legend(rfp_data_legend, loc='upper left', prop={'size': 12},frameon=False)
    #plt.title('Experimental Peptide Concentration vs Time')
    plt.xlabel('Time (min)')
    plt.ylabel('RFP/OD (M/OD)')

    # Set Y Axis Ticker to Scientific Style
    
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
    #plt.rcParams["figure.figsize"] = 5,3.5  # Set figure width to 12 and height to 9
    #plt.figure(dpi = 600)
    # Figure border
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    axes = plt.gca()
    ymin, ymax = axes.get_ylim()
    axes.set_ylim([0, ymax+0.25*ymax])

##############################################################################
# Helper function to read raw data, convert units, and plot the processed data  
##############################################################################

def run_readdata(input_data_filename, NumDataSet, inducer_const):
    print('\tFilename:', input_data_filename)
    print('\tInducer Unit:', inducer_const)

    ### Execute readData Function ###
    # data_header: Deader row from the raw data file
    # data_array: Extracted RFP data in a new array
    data_header, data_array, Data_stddev = readData(input_data_filename, NumDataSet)

    ### Find Sample Size ###
    Data_array_numrows = data_array.shape[0]    # Time + All Inducers
    Data_array_numcols = data_array.shape[1]    # Number of RFP Data Per Inducer
    sample_size = (Data_array_numrows - 1) * Data_array_numcols

    ### Convert Inducer Unit ###
    ### Independent variable # Inducer concentration tuple
    if inducer_const == '%':
        inducer_molar_mass = input("Please insert the Inducer Molar Mass (in g/mol): ")
        inducer_molar_mass = float(inducer_molar_mass)
        inducer = indconv_Percent(data_header, inducer_molar_mass)
    elif inducer_const == 'dimensionless':
        inducer = indconv_None(data_header)
    else:   # for cases with 10**-6 or 10**-9
        inducer = indconv(data_header, inducer_const)

    print('Inducers Used:', inducer)  # ind for checking (Values after conversion to M)

    ### Convert inducer to log10
    inducer_log = np.log10(inducer)
    print('Inducers (log10):', inducer_log)
    print(type(inducer_log))
    print(np.shape(inducer_log))

    
    ### Plot Processed Data ###
    
    ### Execute plotData Function ###
    # Figure 1 - Processed Data
    plot_inputdata(data_header, data_array, Data_stddev)

    return(data_header, data_array, Data_stddev, inducer, inducer_log, sample_size)