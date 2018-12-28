# -*- coding: utf-8 -*-
"""
Created on Wed Aug  1 17:02:29 2018

@author: jingwui
"""

#!/usr/bin/env python

### This module serves as a prefitting algorithm to fit the dose-response data ###
### for Inducible System ###


import numpy as np
import matplotlib.pyplot as plt

from scipy.optimize import minimize
from scipy.optimize import differential_evolution

### ODE Model for Inducible Promoter ###
def compute_hilleqn(inducer, param, Inhibition):

    # print('param', param)
    if Inhibition == True:
        n_Km = param[0]
        K_ind_Km = param[1]
        Vm_Km = param[2]
        Kinhmax = param[3]
        ninh = param[4]
        Kinh = param[5]
        Kbasal = 0 #param[6]

        numInd = len(inducer)
        model_peakrfp = np.zeros(numInd)

        for i in range(0,numInd-1):
            ind = inducer[i]
            # print('Ind', ind)
            # model_peakrfp[i][0] = ind
            model_peakrfp[i] = Kbasal + Vm_Km*((ind**n_Km)/(ind**n_Km + K_ind_Km**n_Km))* \
                         (1-Kinhmax*((ind**ninh)/((ind**ninh)+(Kinh**ninh))))
    elif (Inhibition == False) or (Inhibition == None):
        ### Parameters
        n_Km = param[0]
        K_ind_Km = param[1]
        Vm_Km = param[2]
        Kbasal = 0; #param[3] #0

        numInd = len(inducer)
        model_peakrfp = np.zeros(numInd)

        for i in range(0,numInd-1):
            ind = inducer[i]
            # print('Ind', ind)
            # model_peakrfp[i][0] = ind
            model_peakrfp[i] = Kbasal + Vm_Km*((ind**n_Km) / (ind**n_Km + K_ind_Km**n_Km))

    # print('model', model_peakrfp)
    return model_peakrfp

def findSSE_Km(param, inducer, rfp_data, Inhibition):

    # model_peakrfp = np.zeros(numInd)
    model_peakrfp = compute_hilleqn(inducer, param, Inhibition)

    numInd = len(inducer)

    rfp_max = np.zeros(numInd)
    for i in range(0, numInd):
        # +1 because index 0 is time
        rfp_max[i] = max(rfp_data[i + 1])


    ###Calculate SSE

    sse_Km = 0

    for i in range (0,numInd-1):
        sse_Km = sse_Km + (rfp_max[i] - model_peakrfp[i]) ** 2

    # print('sse_Km', sse_Km)

    # update and return to main function

    return sse_Km

### Plot CSV and Model Data ###
def plotData_Km(raw_data_header, rfp_data, data_stddev, inducer, inducer_log, param_optimized, Inhibition):

#    ### Parameters
#    n_Km = param_optimized[0]
#    K_ind_Km = param_optimized[1]
#    Vm_Km = param_optimized[2]

    # n_Km = 1
    # K_ind_Km = 5e-08
    # Vm_Km = 6.44006680841e-06
    numInd = len(inducer)

    rfp_max = np.zeros(numInd)
    for i in range(0, numInd):
        # +1 because index 0 is time
        rfp_max[i] = max(rfp_data[i + 1])

    rfp = np.zeros(numInd)

    rfp = compute_hilleqn(inducer, param_optimized, Inhibition)

#    for i in range(1,numInd):
#        ind = inducer[i]
#        # model_peakrfp[i][0] = ind
#        rfp[i] = (Vm_Km * ((ind**n_Km) / (ind**n_Km + K_ind_Km**n_Km)))

    # print(rfp)
    fig = plt.figure(figsize=(5,3.6))
    ax = fig.add_axes([0.16,0.17,0.8,0.77])
    plt.rc('font', size=16)  # controls default text sizes
    # Hide the right and top spines
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    plt.tight_layout
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    
    print(type(data_stddev))
    print(np.shape(data_stddev))
    
    ax.errorbar(inducer_log, rfp_max, yerr = data_stddev[:,-1], capsize = 2, color = 'blue', linestyle='None', marker = 's', mfc= 'blue', mec= 'blue', markersize = 4, label = 'Experiment')
    #ax.plot(inducer_log, rfp_max, 'bs',markersize=4, label = 'Experiment')
    ax.plot(inducer_log, rfp, 'b', linewidth=2, label = 'Model')
    ax.legend(frameon=False, loc='upper left')
#    axes = plt.gca()
#    ymin, ymax = axes.get_ylim()
#    axes.set_ylim([0, ymax+0.1*ymax])
    #axes.set_ylim([0, ymax+0.1*ymax])
    #plt.title('Modelled Max M/OD vs Conc (Log)')
    plt.xlabel('Inducer Concentration ($log_{10}$ (M))')
    plt.ylabel('Expression Level (M/OD)')
    plt.show()


##############################################################################
# Helper function to run the dose-response prefitting
##############################################################################

def run_DoseRes(data_header, data_array, data_stddev, inducer, inducer_log, Inhibition):

    ### Find rfp Max for each inducer

    numInd = len(inducer)
    rfp_max = np.zeros(numInd)
    
    print(np.size(data_array))
    
    for i in range(0, numInd):
        # +1 because index 0 is time
        rfp_max[i] = max(data_array[i + 1])
    # print('rfp_max', rfp_max)


    ### Run Optimizer ###

    ### Number of Parameters to be optimized

    if Inhibition == True:
        ### In order of (n_Km, K_ind_Km, Vm_Km, Kinhmax, ninh, Kinh)
        numParam_Km = 6

        Vm_max = max(rfp_max)
        #Vm_min = min(rfp_max)
        ind_mean = (inducer_log[0]+inducer_log[-2])/2
        K_ind0 = 10**ind_mean # Convert back from log10
        K_inh0 = 10**inducer_log[0] #[-1]

        param0_global_Km = [(0.1, 4), (K_ind0*0.8, K_ind0*1.2), (Vm_max*0.9, Vm_max*1.1),
                            (0, 1), (0.1, 4), (K_inh0*0.5, K_inh0*1.5)]#, (0, Vm_min)] #Diff Evo
    elif (Inhibition == False) or (Inhibition == None):
        ### In order of (n_Km, K_ind_Km, Vm_Km)
        numParam_Km = 3; #3 # Fixed for aTc Model

        Vm_max = max(rfp_max)
        #Vm_min = min(rfp_max)
        ind_mean = (inducer_log[0]+inducer_log[-2])/2
        K_ind0 = 10**ind_mean # Convert back from log10

        param0_global_Km = [(0.1, 4), (K_ind0*0.8, K_ind0*1.2), (Vm_max*0.9, Vm_max*1.1)] # (0, Vm_min)]
        #param0_global_Km = [(0.1, 4), (K_ind0*0.8, K_ind0*1.2), (Vm_max*0.9, Vm_max*1.1)]#, (0, Vm_min)]

    result_diffevo_Km = differential_evolution(findSSE_Km,
                                            param0_global_Km, args=(inducer, data_array, Inhibition), tol = 1e-12)

    param0_local_Km = np.zeros(numParam_Km)
    for i in range(0, numParam_Km):
        param0_local_Km[i] = result_diffevo_Km.x[i]

    result_NM_Km = minimize(findSSE_Km, param0_local_Km, args=(inducer, data_array, Inhibition),
                      method='nelder-mead', options={'xtol': 1e-15, 'disp': False})

    param_optimized_Km = np.zeros(numParam_Km)
    for i in range(0, numParam_Km):
        param_optimized_Km[i] = result_NM_Km.x[i]

    ### Calculate and plot Model results (param_optimized) ###
    
    # param_optimized_Km = [3.05784966e+00,  4.57058298e-09,   6.28748162e-06]
    plotData_Km(data_header, data_array, data_stddev, inducer, inducer_log, param_optimized_Km, Inhibition)

    
    ### Display plots
    
    # plt.show()
    if Inhibition == True:
        n, K_ind, Kinh_max, ninh, Kinh = param_optimized_Km[0:5]
        return (n, K_ind, Kinh_max, ninh, Kinh)

    elif (Inhibition == False) or (Inhibition == None) :
        n, K_ind = param_optimized_Km[0:2]
        return(n, K_ind)
    else:
        return print('Error in setting Inhibition!')