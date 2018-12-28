# -*- coding: utf-8 -*-
"""
Created on Mon Jun 25 09:50:57 2018

@author: JingWui
"""

import numpy as np
import constrNMPy as cNM
import matplotlib.pyplot as plt

from scipy.integrate import odeint
from scipy.optimize import differential_evolution

class InduciblePromoterLibrary:

    ### ODE Model for Constant Inducible Promoter ###
    def solveODE_ConstantInducer(y, t, inducer, param, Operation = 'Solve'):
        # Dependent variables
        mRNA = y[0] # Col 0 of ODESoln
        Pep = y[1]  # Col 1 of ODESoln

        # Parameters
        n = param[0]
        K_ind = param[1]
        syn_mRNA = param[2]
        deg_mRNA = 0.1386   # Fixed as ln2/5(mins) = 0.1386
        syn_Pep = param[3]
        deg_Pep = param[4]

        # Differential equations
        dmRNA_dt = '(syn_mRNA*((inducer**n)/(inducer**n+K_ind**n)))-(deg_mRNA * mRNA)'
        dPep_dt = '(syn_Pep*mRNA)-(deg_Pep*Pep)'
        
        if Operation == 'Solve':
            return [eval(dmRNA_dt), eval(dPep_dt)]
        elif Operation == 'GetODE':
            return [dmRNA_dt, dPep_dt]
        else:
            print('Error: Please Enter the correct Operation for solveODE function')

        # Return differential equations solution
    
    def solveODE_ConstantInducerKMat(y, t, inducer, param, Operation = 'Solve'):
        # Dependent variables
        mRNA = y[0] # Col 0 of ODESoln
        Pep = y[1]  # Col 1 of ODESoln
        Pepm = y[2]

        # Parameters
        n = param[0]
        K_ind = param[1]
        syn_mRNA = param[2]
        deg_mRNA = 0.1386   # Fixed as ln2/5(mins) = 0.1386
        syn_Pep = param[3]
        deg_Pep = param[4]
        Kmature = param[5]
        #Kmature = 0.0316 # log(2)/21.9 (mRFP1 at 37'C)

        # Differential equations
        dmRNA_dt = '(syn_mRNA*((inducer**n)/(inducer**n+K_ind**n)))-(deg_mRNA * mRNA)'
        dPep_dt = '(syn_Pep*mRNA)-(Kmature*Pep)'
        dPepm_dt = '(Kmature*Pep)-(deg_Pep*Pepm)'
        
        if Operation == 'Solve':
            return [eval(dmRNA_dt), eval(dPep_dt), eval(dPepm_dt)]
        elif Operation == 'GetODE':
            return [dmRNA_dt, dPep_dt, dPepm_dt]
        else:
            print('Error: Please Enter the correct Operation for solveODE function')

        
    ### ODE Model for Degradable Inducible Promoter ###
    def solveODE_DegradationInducer(y, t, inducer, param, Operation = 'Solve'):
        # Dependent variables
        Ind = y[0]  # Col 0 of ODESoln
        mRNA = y[1] # Col 1
        Pep = y[2]  # Col 2

        # Parameters
        n = param[0]
        K_ind = param[1]
        syn_mRNA = param[2]
        deg_mRNA = 0.1386 # Fixed as ln2/5(mins) = 0.1386
        syn_Pep = param[3]
        deg_Pep = param[4]
        deg_Ind = param[5]
        

        # Differential equations
        dInd_dt = '-deg_Ind*Ind'
        dmRNA_dt = '(syn_mRNA*((Ind**n)/(Ind**n+K_ind**n)))-(deg_mRNA*mRNA)'
        dPep_dt = '(syn_Pep*mRNA)-(deg_Pep*Pep)'

        if Operation == 'Solve':
            return [eval(dInd_dt), eval(dmRNA_dt), eval(dPep_dt)]
        elif Operation == 'GetODE':
            return [dInd_dt, dmRNA_dt, dPep_dt]
        else:
            print('Error: Please Enter the correct Operation for solveODE function')
            
    ### ODE Model for Degradable Inducible Promoter ###
    def solveODE_DegradationInducerKMat(y, t, inducer, param, Operation = 'Solve'):
        # Dependent variables
        Ind = y[0]  # Col 0 of ODESoln
        mRNA = y[1] # Col 1
        Pep = y[2]  # Col 2
        Pepm = y[3]

        # Parameters
        n = param[0]
        K_ind = param[1]
        syn_mRNA = param[2]
        deg_mRNA = 0.1386 # Fixed as ln2/5(mins) = 0.1386
        syn_Pep = param[3]
        deg_Pep = param[4]
        deg_Ind = param[5]
        Kmature = param[6]
        

        # Differential equations
        dInd_dt = '-deg_Ind*Ind'
        dmRNA_dt = '(syn_mRNA*((Ind**n)/(Ind**n+K_ind**n)))-(deg_mRNA*mRNA)'
        dPep_dt = '(syn_Pep*mRNA)-(Kmature*Pep)'
        dPepm_dt = '(Kmature*Pep)-(deg_Pep*Pepm)'

        if Operation == 'Solve':
            return [eval(dInd_dt), eval(dmRNA_dt), eval(dPep_dt), eval(dPepm_dt)]
        elif Operation == 'GetODE':
            return [dInd_dt, dmRNA_dt, dPep_dt, dPepm_dt]
        else:
            print('Error: Please Enter the correct Operation for solveODE function')


    ### ODE Model for Delay Inducible feature ###
    def solveODE_DelayInducer(y, t, inducer, param, Operation = 'Solve'):
        # Dependent variables
        Inde = y[0] # Col 0 of ODESoln
        Indi = y[1] # Col 1 of ODESoln
        mRNA = y[2] # Col 2 of ODESoln
        Pep = y[3]  # Col 3 of ODESoln

        # Parameters
        Vm = param[0]
        ntrans = param[1]
        Ktrans = param[2]
        n = param[3]
        K_ind = param[4]
        syn_mRNA = param[5]
        deg_mRNA = 0.1386 # Fixed as ln2/5(mins) = 0.1386
        syn_Pep = param[6]
        deg_Pep = param[7]

        # Differential equations
        dInde_dt = '-Vm*((Inde**ntrans)/(Inde**ntrans+Ktrans**ntrans))'
        dIndi_dt = 'Vm*((Inde**ntrans)/(Inde**ntrans+Ktrans**ntrans))'
        dmRNA_dt = '(syn_mRNA*((Indi**n)/(Indi**n+K_ind**n)))-(deg_mRNA*mRNA)'
        dPep_dt = '(syn_Pep*mRNA)-(deg_Pep*Pep)'

        if Operation == 'Solve':
            return [eval(dInde_dt), eval(dIndi_dt), eval(dmRNA_dt), eval(dPep_dt)]
        elif Operation == 'GetODE':
            return [dInde_dt, dIndi_dt, dmRNA_dt, dPep_dt]
        else:
            print('Error: Please Enter the correct Operation for solveODE function')
    
    ### ODE Model for Delay Inducible feature ###
    def solveODE_DelayInducerKMat(y, t, inducer, param, Operation = 'Solve'):
        # Dependent variables
        Inde = y[0] # Col 0 of ODESoln
        Indi = y[1] # Col 1 of ODESoln
        mRNA = y[2] # Col 2 of ODESoln
        Pep = y[3]  # Col 3 of ODESoln
        Pepm = y[4]

        # Parameters
        Vm = param[0]
        ntrans = param[1]
        Ktrans = param[2]
        n = param[3]
        K_ind = param[4]
        syn_mRNA = param[5]
        deg_mRNA = 0.1386 # Fixed as ln2/5(mins) = 0.1386
        syn_Pep = param[6]
        deg_Pep = param[7]
        Kmature = param[8]

        # Differential equations
        dInde_dt = '-Vm*((Inde**ntrans)/(Inde**ntrans+Ktrans**ntrans))'
        dIndi_dt = 'Vm*((Inde**ntrans)/(Inde**ntrans+Ktrans**ntrans))'
        dmRNA_dt = '(syn_mRNA*((Indi**n)/(Indi**n+K_ind**n)))-(deg_mRNA*mRNA)'
        dPep_dt = '(syn_Pep*mRNA)-(Kmature*Pep)'
        dPepm_dt = '(Kmature*Pep)-(deg_Pep*Pepm)'

        if Operation == 'Solve':
            return [eval(dInde_dt), eval(dIndi_dt), eval(dmRNA_dt), eval(dPep_dt), eval(dPepm_dt)]
        elif Operation == 'GetODE':
            return [dInde_dt, dIndi_dt, dmRNA_dt, dPep_dt, dPepm_dt]
        else:
            print('Error: Please Enter the correct Operation for solveODE function')

    

    ### Single-ODE Model for Constant Inducible Promoter ###
    def solveODE_SingleConstantInducer(y, t, inducer, param, Operation = 'Solve'):
        # Dependent variables
        Pep = y[0]  # Col 1 of ODESoln

        # Parameters
        n = param[0]
        K_ind = param[1]
        syn_Pep = param[2]
        deg_Pep = param[3]

        # Differential equations
        dPep_dt = '(syn_Pep*((inducer**n)/(inducer**n+K_ind**n)))-(deg_Pep*Pep)'

        if Operation == 'Solve':
            return [eval(dPep_dt)]
        elif Operation == 'GetODE':
            return [dPep_dt]
        else:
            print('Error: Please Enter the correct Operation for solveODE function')

    ### Single-ODE Model for Inducer with Degradation ###
    def solveODE_SingleDegradationInducer(y, t, inducer, param, Operation = 'Solve'):
        # Dependent variables
        Ind = y[0]  # Col 0 of ODESoln
        Pep = y[1]  # Col 2

        # Parameters
        n = param[0]
        K_ind = param[1]
        syn_Pep = param[2]
        deg_Pep = param[3]
        deg_Ind = param[4]

        # Differential equations
        dInd_dt = '-deg_Ind*Ind'
        dPep_dt = '(syn_Pep*((Ind**n)/(Ind**n+K_ind**n)))-(deg_Pep*Pep)'

        if Operation == 'Solve':
            return [eval(dInd_dt), eval(dPep_dt)]
        elif Operation == 'GetODE':
            return [dInd_dt, dPep_dt]
        else:
            print('Error: Please Enter the correct Operation for solveODE function')

    ### Single-ODE Model for Delay Inducible feature ###
    def solveODE_SingleDelayInducer(y, t, inducer, param, Operation = 'Solve'):
        # Dependent variables
        Inde = y[0] # Col 0 of ODESoln
        Indi = y[1] # Col 1 of ODESoln
        Pep = y[2]  # Col 3 of ODESoln

        # Parameters
        Vm = param[0]
        ntrans = param[1]
        Ktrans = param[2]
        n = param[3]
        K_ind = param[4]
        syn_Pep = param[5]
        deg_Pep = param[6]
        

        # Differential equations
        dInde_dt = '-Vm*((Inde**ntrans)/(Inde**ntrans+Ktrans**ntrans))'
        dIndi_dt = 'Vm*((Inde**ntrans)/(Inde**ntrans+Ktrans**ntrans))'

        dPep_dt = '(syn_Pep*((Indi**n)/(Indi**n+K_ind**n)))-(deg_Pep*Pep)'

        if Operation == 'Solve':
            return [eval(dInde_dt), eval(dIndi_dt), eval(dPep_dt)]
        elif Operation == 'GetODE':
            return [dInde_dt, dIndi_dt, dPep_dt]
        else:
            print('Error: Please Enter the correct Operation for solveODE function')
            
            
    ### ODE Model for Constant Inducible Promoter ###
    def solveODE_ConstantIndInhibition(y, t, inducer, param, Operation = 'Solve'):
        ### Dependent variables
        mRNA = y[0] # Col 0 of ODESoln
        Pep = y[1]  # Col 1 of ODESoln

        ### Parameters
        n = param[0]
        K_ind = param[1]
        syn_mRNA = param[2]
        deg_mRNA = 0.1386   # Fixed as ln2/5(mins) = 0.1386
        syn_Pep = param[3]
        deg_Pep = param[4]
        Kinhmax = param[5]
        ninh = param[6]
        Kinh = param[7]
        KLeak = param[8]

        ### Differential equations
        dmRNA_dt = 'KLeak + (syn_mRNA * ((inducer**n) / (inducer**n + K_ind**n))* (1-Kinhmax*((inducer**ninh)/((inducer**ninh)+(Kinh**ninh))))) - (deg_mRNA * mRNA)'
        dPep_dt = '(syn_Pep * mRNA) - (deg_Pep * Pep)'

        if Operation == 'Solve':
            return [eval(dmRNA_dt), eval(dPep_dt)]
        elif Operation == 'GetODE':
            return [dmRNA_dt, dPep_dt]
        else:
            print('Error: Please Enter the correct Operation for solveODE function')
            
    ### ODE Model for Constant Inducible Promoter ###
    def solveODE_DelayIndInhibition(y, t, inducer, param, Operation = 'Solve'):
        ### Dependent variables
        Inde = y[0] # Col 0 of ODESoln
        Indi = y[1] # Col 1 of ODESoln
        mRNA = y[2] # Col 2 of ODESoln
        Pep = y[3]  # Col 3 of ODESoln

        ### Parameters
        n = param[0]
        K_ind = param[1]
        syn_mRNA = param[2]
        deg_mRNA = 0.1386   # Fixed as ln2/5(mins) = 0.1386
        syn_Pep = param[3]
        deg_Pep = param[4]
        Kinhmax = param[5]
        ninh = param[6]
        Kinh = param[7]
        KLeak = param[8]
        Vm = param[9]
        ntrans = param[10]
        Ktrans = param[11]
 

        ### Differential equations
        dInde_dt = '-Vm*((Inde**ntrans)/(Inde**ntrans+Ktrans**ntrans))'
        dIndi_dt = 'Vm*((Inde**ntrans)/(Inde**ntrans+Ktrans**ntrans))'
        dmRNA_dt = 'KLeak + (syn_mRNA * ((Indi**n) / (Indi**n + K_ind**n))* (1-Kinhmax*((Inde**ninh)/((Inde**ninh)+(Kinh**ninh))))) - (deg_mRNA * mRNA)'
        dPep_dt = '(syn_Pep * mRNA) - (deg_Pep * Pep)'

        if Operation == 'Solve':
            return [eval(dInde_dt), eval(dIndi_dt), eval(dmRNA_dt), eval(dPep_dt)]
        elif Operation == 'GetODE':
            return [dInde_dt, dIndi_dt, dmRNA_dt, dPep_dt]
        else:
            print('Error: Please Enter the correct Operation for solveODE function')

    ### ODE Model for Delay Degradation Inducible feature ###
    def solveODE_DelayDegradationInducer(y, t, inducer, param, Operation = 'Solve'):
        # Dependent variables
        Inde = y[0] # Col 0 of ODESoln
        Indi = y[1] # Col 1 of ODESoln
        mRNA = y[2] # Col 2 of ODESoln
        Pep = y[3]  # Col 3 of ODESoln

        # Parameters
        Vm = param[0]
        ntrans = param[1]
        Ktrans = param[2]
        n = param[3]
        K_ind = param[4]
        syn_mRNA = param[5]
        deg_mRNA = 0.1386 # Fixed as ln2/5(mins) = 0.1386
        syn_Pep = param[6]
        deg_Pep = param[7]
        deg_Ind = param[8]

        # Differential equations
        dInde_dt = '-Vm*((Inde**ntrans)/(Inde**ntrans+Ktrans**ntrans))'
        dIndi_dt = 'Vm*((Inde**ntrans)/(Inde**ntrans+Ktrans**ntrans)) - deg_Ind*Indi'
        dmRNA_dt = '(syn_mRNA*((Indi**n)/(Indi**n+K_ind**n)))-(deg_mRNA*mRNA)'
        dPep_dt = '(syn_Pep*mRNA)-(deg_Pep*Pep)'

        if Operation == 'Solve':
            return [eval(dInde_dt), eval(dIndi_dt), eval(dmRNA_dt), eval(dPep_dt)]
        elif Operation == 'GetODE':
            return [dInde_dt, dIndi_dt, dmRNA_dt, dPep_dt]
        else:
            print('Error: Please Enter the correct Operation for solveODE function')

    ### ODE Model for Degradation Delay Inducible feature ###
    def solveODE_DegradationDelayInducer(y, t, inducer, param, Operation = 'Solve'):
        # Dependent variables
        Inde = y[0] # Col 0 of ODESoln
        Indi = y[1] # Col 1 of ODESoln
        mRNA = y[2] # Col 2 of ODESoln
        Pep = y[3]  # Col 3 of ODESoln

        # Parameters
        Vm = param[0]
        ntrans = param[1]
        Ktrans = param[2]
        n = param[3]
        K_ind = param[4]
        syn_mRNA = param[5]
        deg_mRNA = 0.1386 # Fixed as ln2/5(mins) = 0.1386
        syn_Pep = param[6]
        deg_Pep = param[7]
        deg_Ind = param[8]
        
#        ParamName = ['Vm','ntrans','Ktrans','n','K_ind','syn_mRNA','syn_Pep','deg_Pep','deg_Ind']
#        ParamUnits = ['molL-1min-1','dimensionless','molL-1','dimensionless', 'molL-1', 'molL-1min-1', 'min-1', 'min-1', 'min-1']

        # Differential equations
        dInde_dt = '-Vm*((Inde**ntrans)/(Inde**ntrans+Ktrans**ntrans)) - deg_Ind*Inde'
        dIndi_dt = 'Vm*((Inde**ntrans)/(Inde**ntrans+Ktrans**ntrans))'
        dmRNA_dt = '(syn_mRNA*((Indi**n)/(Indi**n+K_ind**n)))-(deg_mRNA*mRNA)'
        dPep_dt = '(syn_Pep*mRNA)-(deg_Pep*Pep)'

        if Operation == 'Solve':
            return [eval(dInde_dt), eval(dIndi_dt), eval(dmRNA_dt), eval(dPep_dt)]
        elif Operation == 'GetODE':
            return [dInde_dt, dIndi_dt, dmRNA_dt, dPep_dt]
        else:
            print('Error: Please Enter the correct Operation for solveODE function')


    function_mappings = {
        'solveODE_ConstantInducer': solveODE_ConstantInducer,
        'solveODE_ConstantInducerKMat': solveODE_ConstantInducerKMat,
        'solveODE_DegradationInducer': solveODE_DegradationInducer,
        'solveODE_DegradationInducerKMat': solveODE_DegradationInducerKMat,
        'solveODE_DelayInducer': solveODE_DelayInducer,
        'solveODE_DelayInducerKMat': solveODE_DelayInducerKMat,
        'solveODE_SingleConstantInducer': solveODE_SingleConstantInducer,
        'solveODE_SingleDegradationInducer': solveODE_SingleDegradationInducer,
        'solveODE_SingleDelayInducer': solveODE_SingleDelayInducer,
        'solveODE_ConstantIndInhibition': solveODE_ConstantIndInhibition,
        'solveODE_DelayIndInhibition': solveODE_DelayIndInhibition,
        'solveODE_DelayDegradationInducer': solveODE_DelayDegradationInducer,
        'solveODE_DegradationDelayInducer': solveODE_DegradationDelayInducer,
        
    }

    def select_function(self, SystemTypeString):
        ### Convert SystemType from string to function name
        try:
            return self.function_mappings[SystemTypeString]
        except KeyError:
            print('Invalid function, try again.')

    def ComputeSSE(self, param, y0, inducer, rfp_data, time_int, SystemType, VarIndex, OptimizerType):
        ### Calculate SSE - To be minimized by Optimizer ####

        global weight_global, sse_global

        # Time grid == no. of times to report solution
        rfp_data_numrows = np.size(rfp_data, 0)
        rfp_data_numcols = np.size(rfp_data, 1)
        numInd = np.size(inducer)
        
        print(rfp_data_numrows)
        print(rfp_data_numcols)
        print(numInd)

        t_start = rfp_data[0][0] # First value of Time column
        t_end = rfp_data[0][-1] # Last value of Time column
        dt = 1  # minutes
        timestep = int((t_end / dt) + 1)
        t = np.linspace(t_start, t_end, timestep)

        # Initialize mRNA and Pep nested list
        #mRNA = np.zeros((timestep, numInd), dtype=object)  # timestep (rows) x numInd (cols)
        Pep = np.zeros((timestep, numInd), dtype=object)   # timestep (rows) x numInd (cols)
        # Pep2 is transpose of Pep, used to find Pep max
        Pep2 = np.zeros((numInd, timestep), dtype=object)  # numInd(rows) x timestep (cols)

        solveODE_Name = '_'.join(('solveODE', SystemType))
        solveODEfun = self.select_function(solveODE_Name) #convert string to function name

        if OptimizerType == 'Global':
            # Integrate the ODE equations
            for i in range(0, numInd):  # Iterates through Ind Concs
                if (SystemType == 'DelayInducer') or (SystemType == 'DegradationInducer') or \
                (SystemType == 'SingleDelayInducer') or (SystemType == 'SingleDegradationInducer') or \
                (SystemType == 'DelayDegradationInducer') or (SystemType == 'DegradationDelayInducer') or \
                (SystemType == 'DegradationInducerKMat') or (SystemType == 'DelayInducerKMat'):
                    y0[0] = inducer[i]  # add inducer conc into the initial condition array y0
                else:
                    pass
                ODEsoln = odeint(solveODEfun, y0, t, args=(inducer[i], param))
#                for j in range(0, timestep):  # Iterates through timesteps
#                    #mRNA[j][i] = ODEsoln[j][VarIndex[0]] # mRNA array runs downwards with time
#                    Pep[j][i] = ODEsoln[j][VarIndex[0]]  # Pep array runs downwards with time
#                    Pep2[i][j] = ODEsoln[j][VarIndex[0]]  # Pep2 array runs sideways with time

        elif OptimizerType == 'Local':
            # Integrate the ODE equations
            for i in range(0, numInd):  # Iterates through Ind Concs
                if (SystemType == 'DelayInducer') or (SystemType == 'DegradationInducer') or \
                (SystemType == 'SingleDelayInducer') or (SystemType == 'SingleDegradationInducer') or \
                (SystemType == 'DelayDegradationInducer') or (SystemType == 'DegradationDelayInducer') or \
                (SystemType == 'DegradationInducerKMat') or (SystemType == 'DelayInducerKMat'):
                    y0[0] = inducer[i]  # add inducer conc into the initial condition array y0
                else:
                    pass
                ODEsoln = odeint(solveODEfun, y0, t, args=(inducer[i], param))
                for j in range(0, timestep):  # Iterates through timesteps
                    #mRNA[j][i] = ODEsoln[j][VarIndex[0]] # mRNA array runs downwards with time
                    Pep[j][i] = ODEsoln[j][VarIndex[0]]  # Pep array runs downwards with time
                    Pep2[i][j] = ODEsoln[j][VarIndex[0]]  # Pep2 array runs sideways with time
        else:
            print('No specified Optimizer Type')
        '''
        Calculate SSE_Time
        '''
        # rfp_data runs lengthwise with time
        sse_time = 0
        for i in range(1, rfp_data_numrows):  # Start from 1 because Row 0 is time
            for j in range(0, rfp_data_numcols):
                # time_int * i to find the corresponding model Pep value
                sse_time = sse_time + (Pep[int(time_int/dt)*j, i-1] - rfp_data[i][j])**2

        # Array to store max data and model rfp values for each inducer
        rfp_peak = np.zeros(numInd)
        Pep_peak = np.zeros(numInd)

        for i in range(0, numInd):
            rfp_peak[i] = max(rfp_data[i+1])
            Pep_peak[i] = max(Pep2[i])

        '''
        Calculate SSE_Dose (Dosage Response)
        '''
        # rfp_data runs lengthwise with time
        sse_dose = 0
        for i in range(0, numInd):  # Start from 1 because Row 0 is time
            sse_dose = sse_dose + (Pep_peak[i] - rfp_peak[i]) ** 2

        weight = np.real(sse_time)/np.real(sse_dose)
        sse = sse_time + sse_dose * weight

        # update and return to main function
        sse_global = sse
        weight_global = weight
        if OptimizerType == 'Global':
            print('Model: ', SystemType, '- SSE (Global):', sse)
        elif OptimizerType == 'Local':
            print('Model: ', SystemType, '- SSE (Local):', sse)
        else:
            print('Error in ComputeSSE function')
        return sse

    def Run_InducibleSystem(self, SystemType, n_DoseRes, K_ind_DoseRes, data_header, data_array, inducer, inducer_log):

        Time_interval = data_array[0][1] - data_array[0][0]
        
        ### ODE Input (Inducer Degradation) ###
        # Initial conditions for (mRNA, Pep, Ind) at time = 0
        if (SystemType == 'ConstantInducer'):
            mRNA0 = 0.
            Pep0 = data_array[1][0]
            y0 = [mRNA0, Pep0]  # Initial condition
            VariableName = ['mRNA', 'Pep'] # Variables Name
            VarIndex =[VariableName.index('Pep')]  #get the index for mRNA and RFP

            ### Number of Parameters to be optimized
            numParam = 5 # Fixed for Constant Inducer Model
            
            ParamName = ['n','K_ind','syn_mRNA','syn_Pep','deg_Pep']
            ParamUnits = ['dimensionless', 'molL-1', 'molL-1min-1', 'min-1', 'min-1']

            # Global bounds - In order of (n, Kind, a_mRNA, a_Pep, y_Pep)
            param0_global = [(n_DoseRes*0.8, n_DoseRes*1.2), (K_ind_DoseRes*0.8, K_ind_DoseRes*1.2),
                             (5e-8, 5e-7), (0, 0.02), (0.001, 0.02)] #Diff Evo

            # Local bounds
            LB = [0]*(numParam-1) + [0.001] #[0.001]
            UB = [None]*(numParam)
            
        if (SystemType == 'ConstantInducerKMat'):
            mRNA0 = 0.
            Pep0 = 0.
            Pepm0 = data_array[1][0]
            y0 = [mRNA0, Pep0, Pepm0]  # Initial condition
            VariableName = ['mRNA', 'Pep', 'Pepm'] # Variables Name
            VarIndex =[VariableName.index('Pepm')]  #get the index for mRNA and RFP

            ### Number of Parameters to be optimized
            numParam = 6 # Fixed for Constant Inducer Model
            
            ParamName = ['n','K_ind','syn_mRNA','syn_Pep','deg_Pep','Kmature']
            ParamUnits = ['dimensionless', 'molL-1', 'molL-1min-1', 'min-1', 'min-1', 'min-1']

            # Global bounds - In order of (n, Kind, a_mRNA, a_Pep, y_Pep)
            param0_global = [(n_DoseRes*0.8, n_DoseRes*1.2), (K_ind_DoseRes*0.8, K_ind_DoseRes*1.2),
                             (5e-8, 5e-7), (0, 0.02), (0.001, 0.02), (0.001, 1)] #Diff Evo

            # Local bounds
            LB = [0]*(numParam-2) + [0.001]*2 #[0.001]
            UB = [None]*(numParam-1) + [1]
#            LB = [0]*(numParam-1) + [0.001] #[0.001]
#            UB = [None]*(numParam)

        elif (SystemType == 'DegradationInducer'):
            # Initial conditions for (Ind, mRNA, Pep) at time = 0
            Ind0 = 0. # Will be replaced in findsse
            mRNA0 = 0.
            Pep0 = data_array[1][0]
            y0 = [Ind0, mRNA0, Pep0]  # Initial condition
            VariableName = ['Ind', 'mRNA', 'Pep'] # Variables Name
            #Variable = ['Inducer', 'mRNA', 'Peptide'] # Variables Name
            VarIndex =[VariableName.index('Pep')]  #get the index for mRNA and RFP

            # Number of Parameters to be optimized
            numParam = 6 # Fixed for Inducer Degradation Model (aTc)
            
            ParamName = ['n','K_ind','syn_mRNA','syn_Pep','deg_Pep','deg_Ind']
            ParamUnits = ['dimensionless', 'molL-1', 'molL-1min-1', 'min-1', 'min-1', 'min-1']

            # Global bounds - In order of (n, K_ind, syn_mRNA, syn_Pep, deg_Pep, deg_Ind)
            param0_global = [(n_DoseRes*0.8, n_DoseRes*1.2), (K_ind_DoseRes*0.8, K_ind_DoseRes*1.2),
                             (3e-8, 5e-7), (0, 0.02), (0.001, 0.02), (0.001, 0.02)] #Diff Evo

            # Local bounds
            LB = [0]*(numParam-2)+[0.001] + [0]
            UB = [None]*(numParam)
            
        elif (SystemType == 'DegradationInducerKMat'):
            # Initial conditions for (Ind, mRNA, Pep) at time = 0
            Ind0 = 0. # Will be replaced in findsse
            mRNA0 = 0.
            Pep0 = 0.
            Pepm0 = data_array[1][0]
            y0 = [Ind0, mRNA0, Pep0, Pepm0]  # Initial condition
            VariableName = ['Ind', 'mRNA', 'Pep', 'Pepm'] # Variables Name
            #Variable = ['Inducer', 'mRNA', 'Peptide'] # Variables Name
            VarIndex =[VariableName.index('Pepm')]  #get the index for mRNA and RFP

            # Number of Parameters to be optimized
            numParam = 7 # Fixed for Inducer Degradation Model (aTc)
            
            ParamName = ['n','K_ind','syn_mRNA','syn_Pep','deg_Pep','deg_Ind', 'Kmature']
            ParamUnits = ['dimensionless', 'molL-1', 'molL-1min-1', 'min-1', 'min-1', 'min-1', 'min-1']

            # Global bounds - In order of (n, K_ind, syn_mRNA, syn_Pep, deg_Pep, deg_Ind)
            param0_global = [(n_DoseRes*0.8, n_DoseRes*1.2), (K_ind_DoseRes*0.8, K_ind_DoseRes*1.2),
                             (3e-8, 5e-7), (0, 0.02), (0.001, 0.02), (0.001, 0.02), (0.001, 1)] #Diff Evo

            # Local bounds
            LB = [0]*(numParam-3)+[0.001] + [0]+ [0.001]
            UB = [None]*(numParam-1) + [1]

        elif (SystemType == 'DelayInducer'):
            # Initial conditions for (Inde, Indi, mRNA, Pep) at time = 0
            Inde0 = 0. # Will be replaced in findsse
            Indi0 = 0.
            mRNA0 = 0.
            Pep0 = data_array[1][0]
            y0 = [Inde0, Indi0, mRNA0, Pep0]  # Initial condition
            VariableName = ['Inde', 'Indi', 'mRNA', 'Pep'] # Variables Name
            VarIndex =[VariableName.index('Pep')]  #get the index for mRNA and RFP

            ### Number of Parameters to be optimized
            numParam = 8 # User Input
            
            ParamName = ['Vm','ntrans','Ktrans','n','K_ind','syn_mRNA','syn_Pep','deg_Pep']
            ParamUnits = ['molL-1min-1','dimensionless','molL-1','dimensionless', 'molL-1', 'molL-1min-1', 'min-1', 'min-1']

            # Global bounds - In order of (Vm, ntrans, Ktrans, n, Kind, syn_mRNA, syn_Pep, deg_Pep)
            param0_global = [(0, 0.01), (0, 4), (0, 0.01), (n_DoseRes*0.8, n_DoseRes*1.2),
                             (K_ind_DoseRes*0.8, K_ind_DoseRes*1.2), (5e-8, 5e-7),
                             (0, 0.02),(0.001, 0.02)] #Diff Evo IPTG

            # local bounds
            LB = [0]*(numParam-1) +[0.001] #[0.001]
            #UB = [None]*(numParam-2) + [0.02, None]
            UB = [None]*(numParam)
        
        elif (SystemType == 'DelayInducerKMat'):
            # Initial conditions for (Inde, Indi, mRNA, Pep) at time = 0
            Inde0 = 0. # Will be replaced in findsse
            Indi0 = 0.
            mRNA0 = 0.
            Pep0 = 0.
            Pepm0 = data_array[1][0]
            y0 = [Inde0, Indi0, mRNA0, Pep0, Pepm0]  # Initial condition
            VariableName = ['Inde', 'Indi', 'mRNA', 'Pep', 'Pepm'] # Variables Name
            VarIndex =[VariableName.index('Pepm')]  #get the index for mRNA and RFP

            ### Number of Parameters to be optimized
            numParam = 9 # User Input
            
            ParamName = ['Vm','ntrans','Ktrans','n','K_ind','syn_mRNA','syn_Pep','deg_Pep', 'Kmature']
            ParamUnits = ['molL-1min-1','dimensionless','molL-1','dimensionless', 'molL-1', 'molL-1min-1', 'min-1', 'min-1', 'min-1']

            # Global bounds - In order of (Vm, ntrans, Ktrans, n, Kind, syn_mRNA, syn_Pep, deg_Pep)
            param0_global = [(0, 0.01), (0, 4), (0, 0.01), (n_DoseRes*0.8, n_DoseRes*1.2),
                             (K_ind_DoseRes*0.8, K_ind_DoseRes*1.2), (5e-8, 5e-7),
                             (0, 0.02),(0.001, 0.02), (0.001, 1)] #Diff Evo IPTG

            # local bounds
            LB = [0]*(numParam-2) +[0.001]*2 #[0.001]
            #UB = [None]*(numParam-2) + [0.02, None]
            UB = [None]*(numParam-1) + [1]

        elif (SystemType == 'SingleConstantInducer'):
            Pep0 = data_array[1][0]
            y0 = [Pep0]  # Initial condition
            VariableName = ['Pep'] # Variables Name
            VarIndex =[VariableName.index('Pep')]  #get the index for mRNA and RFP

            ### Number of Parameters to be optimized
            numParam = 4 # Fixed for Constant Inducer Model
            
            ParamName = ['n','K_ind','syn_Pep','deg_Pep']
            ParamUnits = ['dimensionless', 'molL-1', 'min-1', 'min-1']

            # Global bounds - In order of (n, Kind, a_mRNA, a_Pep, y_Pep)
            param0_global = [(n_DoseRes*0.8, n_DoseRes*1.2), (K_ind_DoseRes*0.8, K_ind_DoseRes*1.2),
                             (5e-8, 5e-7), (0.001, 0.02)] #Diff Evo

            # Local bounds
            LB = [0]*(numParam-1) + [0.001] #[0.001]
            UB = [None]*(numParam)

        elif (SystemType == 'SingleDegradationInducer'):
            # Initial conditions for (Ind, mRNA, Pep) at time = 0
            Ind0 = 0. # Will be replaced in findsse
            Pep0 = data_array[1][0]
            y0 = [Ind0, Pep0]  # Initial condition
            VariableName = ['Ind', 'Pep'] # Variables Name
            VarIndex =[VariableName.index('Pep')]  #get the index for mRNA and RFP

            # Number of Parameters to be optimized
            numParam = 5 # Fixed for Inducer Degradation Model (aTc)
            
            ParamName = ['n','K_ind','syn_Pep','deg_Pep','deg_Ind']
            ParamUnits = ['dimensionless', 'molL-1', 'min-1', 'min-1', 'min-1']

            # Global bounds - In order of (n, K_ind, syn_mRNA, syn_Pep, deg_Pep, deg_Ind)
            param0_global = [(n_DoseRes*0.8, n_DoseRes*1.2), (K_ind_DoseRes*0.8, K_ind_DoseRes*1.2),
                             (3e-8, 5e-7), (0.001, 0.02), (0.001, 0.02)] #Diff Evo

            # Local bounds
            LB = [0]*(numParam-2)+[0.001] + [0]
            UB = [None]*(numParam)

        elif (SystemType == 'SingleDelayInducer'):
            # Initial conditions for (Inde, Indi, mRNA, Pep) at time = 0
            Inde0 = 0. # Will be replaced in findsse
            Indi0 = 0.
            Pep0 = data_array[1][0]
            y0 = [Inde0, Indi0, Pep0]  # Initial condition
            VariableName = ['Inde', 'Indi', 'Pep'] # Variables Name
            VarIndex =[VariableName.index('Pep')]  #get the index for mRNA and RFP

            ### Number of Parameters to be optimized
            numParam = 7 # User Input
            
            ParamName = ['Vm','ntrans','Ktrans','n','K_ind','syn_Pep','deg_Pep']
            ParamUnits = ['molL-1min-1','dimensionless','molL-1','dimensionless', 'molL-1', 'molL-1min-1', 'min-1']

            # Global bounds - In order of (Vm, ntrans, Ktrans, n, Kind, syn_mRNA, syn_Pep, deg_Pep)
            param0_global = [(0, 0.01), (0, 4), (0, 0.01), (n_DoseRes*0.8, n_DoseRes*1.2),
                             (K_ind_DoseRes*0.8, K_ind_DoseRes*1.2), (5e-8, 5e-7),
                             (0.001, 0.02)] #Diff Evo IPTG

            # local bounds
            LB = [0]*(numParam-1) +[0.001] #[0.001]
            UB = [None]*(numParam)

        elif (SystemType == 'DelayDegradationInducer') or (SystemType == 'DegradationDelayInducer'):
            # Initial conditions for (Inde, Indi, mRNA, Pep) at time = 0
            Inde0 = 0. # Will be replaced in findsse
            Indi0 = 0.
            mRNA0 = 0.
            Pep0 = data_array[1][0]
            y0 = [Inde0, Indi0, mRNA0, Pep0]  # Initial condition
            VariableName = ['Inde', 'Indi', 'mRNA', 'Pep'] # Variables Name
            VarIndex =[VariableName.index('Pep')]  #get the index for mRNA and RFP

            ### Number of Parameters to be optimized
            numParam = 9 # User Input
            
            ParamName = ['Vm','ntrans','Ktrans','n','K_ind','syn_mRNA','syn_Pep','deg_Pep','deg_Ind']
            ParamUnits = ['molL-1min-1','dimensionless','molL-1','dimensionless', 'molL-1', 'molL-1min-1', 'min-1', 'min-1', 'min-1']

            # Global bounds - In order of (Vm, ntrans, Ktrans, n, Kind, syn_mRNA, syn_Pep, deg_Pep)
            param0_global = [(0, 0.01), (0, 4), (0, 0.01), (n_DoseRes*0.8, n_DoseRes*1.2),
                             (K_ind_DoseRes*0.8, K_ind_DoseRes*1.2), (5e-8, 5e-7),
                             (0, 0.02),(0.001, 0.02),(0.0, 0.02)] #Diff Evo IPTG

            # local bounds
            LB = [0]*(numParam-2) +[0.001, 0] #[0.001]
            #UB = [None]*(numParam-2) + [0.02, None]
            UB = [None]*(numParam)

        elif (SystemType == 'ConstantIndInhibition'):
            mRNA0 = 0.
            Pep0 = data_array[1][0]
            y0 = [mRNA0, Pep0]  # Initial condition
            VariableName = ['mRNA', 'Pep'] # Variables Name
            VarIndex =[VariableName.index('Pep')]

            'Vm', 'ntrans', 'Ktrans'

            ### Number of Parameters to be optimized
            numParam = 9 # Fixed for Inhibition Inducer model
            
            ParamName = ['n','K_ind','syn_mRNA','syn_Pep','deg_Pep','Kinhmax','ninh','Kinh', 'KLeak']
            ParamUnits = ['dimensionless', 'molL-1', 'molL-1min-1', 'min-1', 'min-1', 'dimensionless','dimensionless','molL-1', 'molL-1min-1']

            K_inh0 = 10**inducer_log[0] #[-1]

            ### Global Optimizer (Diff Evo) ###

            ### In order of (n, Kind, a_mRNA, a_Pep, y_Pep)
            param0_global = [(n_DoseRes*0.8, n_DoseRes*1.2), (K_ind_DoseRes*0.8, K_ind_DoseRes*1.2),
                             (5e-8, 5e-7), (0, 0.02), (0.001, 0.02), (0, 1), (0, 4), (K_inh0*0.5, K_inh0*1.5), (1e-10, 5e-6)] #Diff Evo

            # local bounds
            LB = [0]*(numParam-5) + [0.001] +[0]*4 #[0.001]
            UB = [None]*(numParam)
            
        elif (SystemType == 'DelayIndInhibition'):
            Inde0 = 0. # Will be replaced in findsse
            Indi0 = 0.
            mRNA0 = 0.
            Pep0 = data_array[1][0]
            y0 = [Inde0, Indi0, mRNA0, Pep0]  # Initial condition
            VariableName = ['Inde', 'Indi', 'mRNA', 'Pep'] # Variables Name
            VarIndex =[VariableName.index('Pep')]


            ### Number of Parameters to be optimized
            numParam = 12 # Fixed for Inhibition Inducer model
            
            ParamName = ['n','K_ind','syn_mRNA','syn_Pep','deg_Pep','Kinhmax','ninh','Kinh', 'KLeak', 'Vm', 'ntrans', 'Ktrans']
            ParamUnits = ['dimensionless', 'molL-1', 'molL-1min-1', 'min-1', 'min-1', 'dimensionless','dimensionless','molL-1', 'molL-1min-1', 'molL-1min-1','dimensionless','molL-1']

            K_inh0 = 10**inducer_log[0] #[-1]

            ### Global Optimizer (Diff Evo) ###

            ### In order of (n, Kind, a_mRNA, a_Pep, y_Pep)
            param0_global = [(n_DoseRes*0.8, n_DoseRes*1.2), (K_ind_DoseRes*0.8, K_ind_DoseRes*1.2),
                             (5e-8, 5e-7), (0, 0.02), (0.001, 0.02), (0, 1), (0, 4), (K_inh0*0.5, K_inh0*1.5), (1e-10, 5e-6), 
                             (0, 0.01), (0, 4), (0, 0.01)] #Diff Evo

            # local bounds
            LB = [0]*(numParam-8) + [0.001] +[0]*7 #[0.001]
            UB = [None]*(numParam)

        else:
            print('Please choose the correct System Type for Inducible Promoter')

        # run Global optimizer
        OptimizerType1 = 'Global'
        result_diffevo = differential_evolution\
            (self.ComputeSSE, param0_global, args=(y0, inducer, data_array, Time_interval, SystemType, VarIndex, OptimizerType1))


        # run Local Optimizer (Nelder Mead)
        OptimizerType2 = 'Local'
        param0_local = np.zeros(numParam)
        for i in range(0, numParam):
            param0_local[i] = result_diffevo.x[i]

        result_NM = cNM.constrNM(self.ComputeSSE, param0_local,LB,UB,args=(y0, inducer, data_array, Time_interval, SystemType, VarIndex, OptimizerType2),
                                 xtol= 1e-15, full_output=True)

        # Optimized Parameters
        param_optimized = np.zeros(numParam)
        for i in range(0, numParam):
            #param_optimized[i] = result_NM.x[i]
            param_optimized[i] = result_NM['xopt'][i]

        return(param_optimized, sse_global, y0, VariableName, ParamName, ParamUnits)

    # Plot CSV and Model Data #
    def plotData_Combined(self, SystemType, Variable, y0, raw_data_header, rfp_data, Data_stddev, inducer, inducer_log, param):
        ### timespan for Model results (t)
        t_start = rfp_data[0][0]
        t_end = rfp_data[0][-1]
        dt = 1  # minutes
        timestep = int((t_end / dt) + 1)
        t = np.linspace(t_start, t_end, timestep)

        # Time grid == no. of times to report solution
        rfp_data_numrows = np.size(rfp_data, 0)
        rfp_data_numcols = np.size(rfp_data, 1)
        numInd = np.size(inducer)

        Inducer = np.zeros((timestep, numInd), dtype=object)
        Inducer_e = np.zeros((timestep, numInd), dtype=object)
        Inducer_i = np.zeros((timestep, numInd), dtype=object)
        mRNA = np.zeros((timestep, numInd), dtype=object)
        Peptide = np.zeros((timestep, numInd), dtype=object)
        Peptidem = np.zeros((timestep, numInd), dtype=object)
        Pep2 = np.zeros((numInd, timestep), dtype=object)  # numInd (rows) x timestep (cols)

        Variable_mappings = {
                'Inducer': Inducer,
                'Inducer_e': Inducer_e,
                'Inducer_i': Inducer_i,
                'mRNA': mRNA,
                'Peptide': Peptide,
                'Peptidem': Peptidem,
            }

        # initiate empty list to store all the arrays
        VariableMatrix = [None]*len(Variable)

        for i in range(0, len(Variable)):
            VariableMatrix[i] = Variable_mappings[Variable[i]]

        solveODE_Name = '_'.join(('solveODE', SystemType))
        solveODEfun = self.select_function(solveODE_Name) #convert string to function name

        # Integrate the ODE equations
        for i in range(0, numInd):  # Iterates through Ind Concs
            if (SystemType == 'DelayInducer') or (SystemType == 'DegradationInducer') or \
            (SystemType == 'SingleDelayInducer') or (SystemType == 'SingleDegradationInducer') or \
            (SystemType == 'DelayDegradationInducer') or (SystemType == 'DegradationDelayInducer') or \
            (SystemType == 'DegradationInducerKMat') or (SystemType == 'DelayInducerKMat') or \
            (SystemType == 'DelayIndInhibition'):
                y0[0] = inducer[i]  # add inducer conc into the initial condition array y0
            else:
                pass
            ODEsoln = odeint(solveODEfun, y0, t, args=(inducer[i], param))
            for j in range(0, len(VariableMatrix)):  # Iterates through timesteps
                #for k in range(0, timestep):
                VariableMatrix[j][:,i] = ODEsoln[:,j]
                    #VariableMatrix[j][k][i] = ODEsoln[k][j]
            Pep2[i,:] = ODEsoln[:,-1]  # Pep2 array runs sideways with time

        ### Retrive the ODEs in String from the corresponding solveODEfun
        ODEstring = solveODEfun(y0, t[0], inducer[0], param, 'GetODE')
        
        print('ODEs in string', ODEstring)
        
        ### Defining time array from Row 0
        time = rfp_data[0]
        

        ### Plot RFP Data vs Time ###
        fig = plt.figure(figsize=(5,3.6))
        ax = fig.add_axes([0.16,0.16,0.8,0.78])
        ax.set_prop_cycle('color',plt.cm.tab10(np.linspace(0, 1, 10)))
        plt.rc('font', size=16)  # controls default text sizes
        ## CSV Data (time)
        style = ['^','*','>','D','<','d','p','o','h','+','s','x','v','.','H']
        for i in range (1, rfp_data_numrows):
            #plt.plot(time, rfp_data[i], linestyle='None', marker = style[i-1], markersize = 4)
            plt.errorbar(time, rfp_data[i], yerr = Data_stddev[i-1], capsize = 2, linestyle='None', marker = style[i-1], markersize = 3)

        rfp_data_legend = raw_data_header[0:rfp_data_numcols]
        
        plt.legend(rfp_data_legend, ncol=2, loc='upper left', prop={'size': 12},frameon=False, markerscale=1.6)

#        n = 20
#        ax = plt.axes()
#        ax.set_prop_cycle('color',[plt.cm.jet(i) for i in np.linspace(0, 1, n)])

        ### Model Data (t)
        # Resets colour cycle for the second plot over same figure
        ax.set_prop_cycle('color',plt.cm.tab10(np.linspace(0, 1, 10)))
        #plt.gca().set_prop_cycle(None)
        # plot model Pep data   
        if 'Peptidem' not in Variable:
            Peptideid = Variable.index('Peptide') #get the index of mRNA
            plt.plot(t, VariableMatrix[Peptideid][:,:], linewidth=2)  # Pep
        else:
            Peptidemid = Variable.index('Peptidem') #get the index of mRNA
            plt.plot(t, VariableMatrix[Peptidemid][:,:], linewidth=2)  # Pep
        #plt.title('Modelled Peptide Concentration vs Time')
        plt.xlabel('Time (min)')
        plt.ylabel('Expression Level (M/OD)')
        # Set Y Axis Ticker to Scientific Style
        plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
        # Figure border
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        axes = plt.gca()
        ymin, ymax = axes.get_ylim()
        #axes.set_ylim([0, ymax+0.3*ymax]) #fpr other [plots]
        axes.set_ylim([0, ymax+0.1*ymax]) #for pLac
        #axes.set_ylim([0, ymax+0.1*ymax]) 

        ### Protein Data (t)
        fig = plt.figure(figsize=(5,3.6))
        ax = fig.add_axes([0.16,0.16,0.8,0.78])
        if 'Peptidem' not in Variable:
            Peptideid = Variable.index('Peptide') #get the index of mRNA
            plt.plot(t, VariableMatrix[Peptideid][:,:], linewidth=2)  # Pep
        else:
            Peptidemid = Variable.index('Peptidem') #get the index of mRNA
            plt.plot(t, VariableMatrix[Peptidemid][:,:], linewidth=2)  # Pep

        plt.legend(rfp_data_legend, ncol=2, loc='upper left', prop={'size': 12},frameon=False)
        #plt.title('Modelled Peptide Concentration vs Time')
        plt.xlabel('Time (min)')
        plt.ylabel('Expression Level (M/OD)')
        # Set Y Axis Ticker to Scientific Style
        plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
        # Figure border
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        axes = plt.gca()
        ymin, ymax = axes.get_ylim()
        axes.set_ylim([0, ymax+0.35*ymax])

        '''
        ### Dosage Response ###
        '''
        rfp_peak = np.zeros(numInd)
        Pep_peak = np.zeros(numInd)
        for i in range(0, numInd):
            rfp_peak[i] = max(rfp_data[i + 1])
            Pep_peak[i] = max(Pep2[i])

        fig = plt.figure(figsize=(5,3.6))
        ax = fig.add_axes([0.16,0.17,0.8,0.77])
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        plt.tight_layout
        plt.rc('font', size=16)  # controls default text sizes
        ax.errorbar(inducer_log, rfp_peak, yerr = Data_stddev[:,-1], capsize = 2, color = 'black', linestyle='None', marker = 's', mfc= 'black', mec= 'black', markersize = 4, label = 'Experiment')
        #ax.plot(inducer_log, rfp_peak, 'ks',markersize=4, label = 'Experiment')
        ax.plot(inducer_log, Pep_peak, 'k', linewidth=2, label = 'Model')
        ax.legend(frameon=False, loc='upper left')
        plt.xlabel('Inducer Concentration ($log_{10}$ (M))')
        plt.ylabel('Expression Level (M/OD)')
        # Set Y Axis Ticker to Scientific Style
        plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
        axes = plt.gca()
        ymin, ymax = axes.get_ylim()
        #axes.set_ylim([0, ymax+0.1*ymax])
        #axes.set_ylim([0, ymax+0.2*ymax])

        if 'mRNA' in Variable:
            ### mRNA Data (t)
            fig = plt.figure(figsize=(5,3.6))
            ax = fig.add_axes([0.16,0.16,0.8,0.78])
            plt.rc('font', size=16)  # controls default text sizes
            mRNAid = Variable.index('mRNA') #get the index of mRNA
            plt.plot(t, VariableMatrix[mRNAid][:,:], linewidth=2)  # mRNA
            #plt.title('Modelled mRNA Concentration vs Time')
            plt.xlabel('Time (min)')
            plt.ylabel('mRNA Concentration (M)')
            plt.legend(rfp_data_legend, ncol=3, loc='upper left', prop={'size': 12},frameon=False)
            # Set Y Axis Ticker to Scientific Style
            plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
            # Figure border
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            axes = plt.gca()
            ymin, ymax = axes.get_ylim()
            #axes.set_ylim([0, ymax+0.45*ymax])
            axes.set_ylim([0, ymax+0.45*ymax])

        if 'Inducer' in Variable:
            ### Inducer Data (t)
            fig = plt.figure(figsize=(5,3.6))
            ax = fig.add_axes([0.16,0.16,0.8,0.78])
            plt.rc('font', size=16)  # controls default text sizes
            # Hide the right and top spines
            ax.spines['right'].set_visible(False)
            ax.spines['top'].set_visible(False)
            plt.tight_layout
            Inducerid = Variable.index('Inducer') #get the index of Inducer
            plt.plot(t, VariableMatrix[Inducerid][:,:], linewidth=2)  # Inducer
            plt.xlabel('Time (min)')
            plt.ylabel('Inducer Concentration (M)')
            plt.legend(rfp_data_legend, ncol=3, loc='upper left', prop={'size': 12},frameon=False)
            # Set Y Axis Ticker to Scientific Style
            plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
            # Figure border
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            axes = plt.gca()
            ymin, ymax = axes.get_ylim()
            axes.set_ylim([0, ymax+0.35*ymax])

        if 'Inducer_e' in Variable:
            ### Extracellular Inducer Data (t)
            fig = plt.figure(figsize=(5,3.6))
            ax = fig.add_axes([0.16,0.16,0.8,0.78])
            plt.rc('font', size=16)  # controls default text sizes
            # Hide the right and top spines
            ax.spines['right'].set_visible(False)
            ax.spines['top'].set_visible(False)
            plt.tight_layout
            Inducer_eid = Variable.index('Inducer_e') #get the index of Extracellular Inducer
            plt.plot(t, VariableMatrix[Inducer_eid][:,:], linewidth=2)  # Inducer
            plt.xlabel('Time (min)')
            plt.ylabel('Inducer$_{ex}$ Concentration (M)')
            plt.legend(rfp_data_legend, ncol=3, loc='upper left', prop={'size': 12},frameon=False)
            # Set Y Axis Ticker to Scientific Style
            plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
            # Figure border
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            axes = plt.gca()
            ymin, ymax = axes.get_ylim()
            axes.set_ylim([0, ymax+0.35*ymax])

        if 'Inducer_i' in Variable:
            ### Intracellular Inducer Data (t)
            fig = plt.figure(figsize=(5,3.6))
            ax = fig.add_axes([0.16,0.16,0.8,0.78])
            plt.rc('font', size=16)  # controls default text sizes
            # Hide the right and top spines
            ax.spines['right'].set_visible(False)
            ax.spines['top'].set_visible(False)
            plt.tight_layout
            Inducer_iid = Variable.index('Inducer_i') #get the index of Intracellular Inducer
            plt.plot(t, VariableMatrix[Inducer_iid][:,:], linewidth=2)  # Inducer_i
            plt.xlabel('Time (min)')
            plt.ylabel('Inducer$_{in}$ Concentration (M)')
            plt.legend(rfp_data_legend, ncol=3, loc='upper left', prop={'size': 12},frameon=False)
            # Set Y Axis Ticker to Scientific Style
            plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
            # Figure border
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            axes = plt.gca()
            ymin, ymax = axes.get_ylim()
            #axes.set_ylim([0, ymax+0.35*ymax])
            axes.set_ylim([0, ymax+0.55*ymax])
            
        return t, VariableMatrix, rfp_data_legend, ODEstring


    def Run_InduciblePlot(self, SystemType, y0, data_header, data_array, Data_stddev, inducer, inducer_log, param_optimized):
        if (SystemType == 'ConstantInducer') or (SystemType == 'ConstantIndInhibition'):
            VariablePlot = ['mRNA', 'Peptide'] # Variables Name
            
        elif (SystemType == 'ConstantInducerKMat'):
            VariablePlot = ['mRNA', 'Peptide', 'Peptidem'] # Variables Name

        elif (SystemType == 'DegradationInducer'):
            VariablePlot = ['Inducer', 'mRNA', 'Peptide'] # Variables Name
            
        elif (SystemType == 'DegradationInducerKMat'):
            VariablePlot = ['Inducer', 'mRNA', 'Peptide', 'Peptidem'] # Variables Name

        elif (SystemType == 'DelayInducer') or (SystemType == 'DelayIndInhibition'):
            VariablePlot = ['Inducer_e', 'Inducer_i', 'mRNA', 'Peptide'] # Variables Name
            
        elif (SystemType == 'DelayInducerKMat'):
            VariablePlot = ['Inducer_e', 'Inducer_i', 'mRNA', 'Peptide', 'Peptidem'] # Variables Name

        elif (SystemType == 'SingleConstantInducer'):
            VariablePlot = ['Peptide'] # Variables Name

        elif (SystemType == 'SingleDegradationInducer'):
            VariablePlot = ['Inducer', 'Peptide'] # Variables Name

        elif (SystemType == 'SingleDelayInducer'):
            VariablePlot = ['Inducer_e', 'Inducer_i', 'Peptide'] # Variables Name

        elif (SystemType == 'DelayDegradationInducer'):
            VariablePlot = ['Inducer_e', 'Inducer_i', 'mRNA', 'Peptide'] # Variables Name

        elif (SystemType == 'DegradationDelayInducer'):
            VariablePlot = ['Inducer_e', 'Inducer_i', 'mRNA', 'Peptide'] # Variables Name
        else:
            print('Error in Plotting module')


        ### Calculate and plot Model results (param_optimized) ###
        t, VariableMatrix, rfp_data_legend, ODEstring = self.plotData_Combined(SystemType, VariablePlot, y0, data_header, data_array, Data_stddev, inducer, inducer_log, param_optimized)
        
        return t, VariableMatrix, rfp_data_legend, ODEstring

    def __del__(self):
      class_name = self.__class__.__name__
      print(class_name, "destroyed")
