# -*- coding: utf-8 -*-
"""
Created on Wed Dec 26 15:21:42 2018

@author: JingWui
"""

#############################################################################
# The Inducible System class
#############################################################################

### Scripts imported ###
import Read_Data as ReadData
import DoseResponse as DoseRes
import ComputeAIC_ as ComputeAIC
import Txtfilename as Txtfilename
### Import Class from module
from InduciblePromoterLibrary_ import InduciblePromoterLibrary

#import numpy as np

### Python Packages Imported ###
import matplotlib.pyplot as plt
import time
import csv
import scipy.stats as ss
from tabulate import tabulate

class InducibleSystem:
    def __init__(self):
        
        plt.close('all')
        
        plt.rcParams['font.family'] = 'Arial' #'sans-serif'
        plt.rcParams['font.weight'] = 'bold'
        plt.rcParams['font.size'] = 16
        plt.rcParams['axes.labelsize'] = 16
        plt.rcParams['axes.labelweight'] = 'bold'
        plt.rcParams['axes.linewidth'] = 2
        plt.rcParams['xtick.labelsize'] = 16
        #plt.rcParams['xtick.labelweight'] = 'medium'
        plt.rcParams['xtick.major.width'] = 2
        plt.rcParams['ytick.labelsize'] = 16
        #plt.rcParams['xtick.labelweight'] = 'medium'
        plt.rcParams['ytick.major.width'] = 2
        plt.rcParams['legend.handletextpad'] = 0.5
        
        ### MSS Main Code ###
        self.starttime = time.time()
        
        #create the first object of InduciblePromoterLibrary class
        self.IndPLib1 = InduciblePromoterLibrary()
        
        self.SystemOptList = ['all', 'constind', 'degradeind', 'delayind', 'inhibition']
        
    def DataReader(self, Input_filename, NumDataSet, Inducer_unit):
        
        ### Read Data ###
        print('Read Data')
        self.Data_header1, self.Data_array1, self.Data_stddev1, self.Inducer1, self.Inducer_log1, self.Sample_size \
            = ReadData.run_readdata(Input_filename, NumDataSet, Inducer_unit)
            
        self.Input_filename_ = Input_filename
        self.NumDataSet_ = NumDataSet
            
    def RunDoseResponsePrefitting(self, OptInhibition = None):
        
        # To calculate initial guess for Param n and K_ind
        print('Dose Response')
        
        self.OptInhibition_ = OptInhibition
        
        if OptInhibition == True:
            # n_DoseRes, K_ind_DoseRes, Kinh_max_DoseRes, ninh_DoseRes, Kinh_DoseRes, Kbasal_DoseRes
            # = DoseRes.run_DoseRes(Data_header1, Data_array1, Inducer1, Inducer_log1, OptInhibition)
            self.n_DoseRes, self.K_ind_DoseRes, Kinh_max_DoseRes, ninh_DoseRes, Kinh_DoseRes \
            = DoseRes.run_DoseRes(self.Data_header1, self.Data_array1, self.Data_stddev1, self.Inducer1, self.Inducer_log1, OptInhibition)
            print('Param estimate for n:', self.n_DoseRes)
            print('Param estimate for Kind :', self.K_ind_DoseRes)
            print('Param estimate for Kinhmax:', Kinh_max_DoseRes)
            print('Param estimate for ninh :', ninh_DoseRes)
            print('Param estimate for Kinh :', Kinh_DoseRes)
            #print('Param estimate for Kbasal :', Kbasal_DoseRes)
        else:
            OptInhibition == None     #default setting
            #n_DoseRes, K_ind_DoseRes, Kbasal_DoseRes = DoseRes.run_DoseRes(Data_header1, Data_array1, Inducer1, Inducer_log1, OptInhibition)
            self.n_DoseRes, self.K_ind_DoseRes = DoseRes.run_DoseRes(self.Data_header1, self.Data_array1, self.Data_stddev1, \
                                                           self.Inducer1, self.Inducer_log1, OptInhibition)
            print('Param estimate for n:', self.n_DoseRes)
            print('Param estimate for Kind :', self.K_ind_DoseRes)
            #print('Param estimate for Kbasal :', Kbasal_DoseRes)
            
    def RunModels(self, SystemOpt = 'ALL', iteration = 1):
        
        # -------------------------------------------------------------------------------- #

        if (SystemOpt.casefold() not in self.SystemOptList):
            raise Exception('Error in the selected System Option')
            
        else:

            ### Define Lists Used ###
            self.Model_List = []     # To store names of Models Tested
            self.SSE_Combined = []   # To store values of SSE
            self.Num_Param_List = [] # To store param length of each model tested
            self.ParamName_List = [] # To store param Name of each model tested
            self.Param_List = [] # To store fitted param values of each model tested
            self.ParamUnits_List = [] # To store the units of the params of each model tested
            self.VarName_List = [] # To store the Variables name of each model tested
            self.y0_List = [] # To store the initial values for the variables of each model tested
            
            # -------------------------------------------------------------------------------- #
            
            for loop in range(0, iteration):
                
                if (SystemOpt.casefold() == 'constind') or (SystemOpt.casefold() == 'all'): 
                
                    ### Model 1 - Constant Inducer ###
                    SystemType = 'ConstantInducer'
                    Param_ConstantInducer, SSE_ConstantInducer, y0_ConstantInducer, VarName_ConstantInducer, ParamName_ConstantInducer, ParamUnits_ConstantInducer \
                        = self.IndPLib1.Run_InducibleSystem(SystemType, self.n_DoseRes, self.K_ind_DoseRes, self.Data_header1, self.Data_array1, self.Inducer1, self.Inducer_log1)
                
                    self.SSE_Combined.append(SSE_ConstantInducer)
                    self.ParamName_List.append(ParamName_ConstantInducer)
                    self.Param_List.append(Param_ConstantInducer)
                    self.ParamUnits_List.append(ParamUnits_ConstantInducer)
                    self.VarName_List.append(VarName_ConstantInducer)
                    self.y0_List.append(y0_ConstantInducer)
                    self.Num_Param_List.append(len(Param_ConstantInducer))
                    self.Model_List.append('Model 1.0 - ConstantInducer')
                
                    ### Model 1.1 - Constant Inducer with Protein Maturation Kinetics ###
                    SystemType = 'ConstantInducerKMat'
                    Param_ConstantInducerKMat, SSE_ConstantInducerKMat, y0_ConstantInducerKMat, VarName_ConstantInducerKMat, ParamName_ConstantInducerKMat, ParamUnits_ConstantInducerKMat \
                        = self.IndPLib1.Run_InducibleSystem(SystemType, self.n_DoseRes, self.K_ind_DoseRes, self.Data_header1, self.Data_array1, self.Inducer1, self.Inducer_log1)
                
                    self.SSE_Combined.append(SSE_ConstantInducerKMat)
                    self.ParamName_List.append(ParamName_ConstantInducerKMat)
                    self.Param_List.append(Param_ConstantInducerKMat)
                    self.ParamUnits_List.append(ParamUnits_ConstantInducerKMat)
                    self.VarName_List.append(VarName_ConstantInducerKMat)
                    self.y0_List.append(y0_ConstantInducerKMat)
                    self.Num_Param_List.append(len(Param_ConstantInducerKMat))
                    self.Model_List.append('Model 1.1 - ConstantInducerKMat')
                    
                    ## Model 1.2 - Single-ODE Constant Inducer ###
                    SystemType = 'SingleConstantInducer'
                    Param_SingleConstantInducer, SSE_SingleConstantInducer, y0_SingleConstantInducer, VarName_SingleConstantInducer, ParamName_SingleConstantInducer, ParamUnits_SingleConstantInducer \
                        = self.IndPLib1.Run_InducibleSystem(SystemType, self.n_DoseRes, self.K_ind_DoseRes, self.Data_header1, self.Data_array1, self.Inducer1, self.Inducer_log1)
                
                    self.SSE_Combined.append(SSE_SingleConstantInducer)
                    self.ParamName_List.append(ParamName_SingleConstantInducer)
                    self.Param_List.append(Param_SingleConstantInducer)
                    self.ParamUnits_List.append(ParamUnits_SingleConstantInducer)
                    self.VarName_List.append(VarName_SingleConstantInducer)
                    self.y0_List.append(y0_SingleConstantInducer)
                    self.Num_Param_List.append(len(Param_SingleConstantInducer))
                    self.Model_List.append('Model 1.2 - SingleConstantInducer')
                
                # -------------------------------------------------------------------------------- #
                
                if (SystemOpt.casefold() == 'degradeind') or (SystemOpt.casefold() == 'all'):
            
                    ### Model 2 - Inducer Degradation ###
                    SystemType = 'DegradationInducer'
                    Param_DegradationInducer, SSE_DegradationInducer, y0_DegradationInducer, VarName_DegradationInducer, ParamName_DegradationInducer, ParamUnits_DegradationInducer \
                        = self.IndPLib1.Run_InducibleSystem(SystemType, self.n_DoseRes, self.K_ind_DoseRes, self.Data_header1, self.Data_array1, self.Inducer1, self.Inducer_log1)
                
                    self.SSE_Combined.append(SSE_DegradationInducer)
                    self.ParamName_List.append(ParamName_DegradationInducer)
                    self.Param_List.append(Param_DegradationInducer)
                    self.ParamUnits_List.append(ParamUnits_DegradationInducer)
                    self.VarName_List.append(VarName_DegradationInducer)
                    self.y0_List.append(y0_DegradationInducer)
                    self.Num_Param_List.append(len(Param_DegradationInducer))
                    self.Model_List.append('Model 2.0 - DegradationInducer')
                    
                    ### Model 2.1 - Inducer Degradation with Protein Maturation Kinetics ###
                    SystemType = 'DegradationInducerKMat'
                    Param_DegradationInducerKMat, SSE_DegradationInducerKMat, y0_DegradationInducerKMat, VarName_DegradationInducerKMat, ParamName_DegradationInducerKMat, ParamUnits_DegradationInducerKMat \
                        = self.IndPLib1.Run_InducibleSystem(SystemType, self.n_DoseRes, self.K_ind_DoseRes, self.Data_header1, self.Data_array1, self.Inducer1, self.Inducer_log1)
                
                    self.SSE_Combined.append(SSE_DegradationInducerKMat)
                    self.ParamName_List.append(ParamName_DegradationInducerKMat)
                    self.Param_List.append(Param_DegradationInducerKMat)
                    self.ParamUnits_List.append(ParamUnits_DegradationInducerKMat)
                    self.VarName_List.append(VarName_DegradationInducerKMat)
                    self.y0_List.append(y0_DegradationInducerKMat)
                    self.Num_Param_List.append(len(Param_DegradationInducerKMat))
                    self.Model_List.append('Model 2.1 - DegradationInducerKMat')
                    
                    ### Model 2.2 - Single-ODE Inducer Degradation ###
                    SystemType = 'SingleDegradationInducer'
                    Param_SingleDegradationInducer, SSE_SingleDegradationInducer, y0_SingleDegradationInducer, VarName_SingleDegradationInducer, ParamName_SingleDegradationInducer, ParamUnits_SingleDegradationInducer \
                        = self.IndPLib1.Run_InducibleSystem(SystemType, self.n_DoseRes, self.K_ind_DoseRes, self.Data_header1, self.Data_array1, self.Inducer1, self.Inducer_log1)
                
                    self.SSE_Combined.append(SSE_SingleDegradationInducer)
                    self.ParamName_List.append(ParamName_SingleDegradationInducer)
                    self.Param_List.append(Param_SingleDegradationInducer)
                    self.ParamUnits_List.append(ParamUnits_SingleDegradationInducer)
                    self.VarName_List.append(VarName_SingleDegradationInducer)
                    self.y0_List.append(y0_SingleDegradationInducer)
                    self.Num_Param_List.append(len(Param_SingleDegradationInducer))
                    self.Model_List.append('Model 2.2 - SingleDegradationInducer')
            
                # -------------------------------------------------------------------------------- #
                
                if (SystemOpt.casefold() == 'delayind') or (SystemOpt.casefold() == 'all'):
            
                    ### Model 3 - Inducer Transport Delay ###
                    SystemType = 'DelayInducer'
                    Param_DelayInducer, SSE_DelayInducer, y0_DelayInducer, VarName_DelayInducer, ParamName_DelayInducer, ParamUnits_DelayInducer \
                        = self.IndPLib1.Run_InducibleSystem(SystemType, self.n_DoseRes, self.K_ind_DoseRes, self.Data_header1, self.Data_array1, self.Inducer1, self.Inducer_log1)
                
                    self.SSE_Combined.append(SSE_DelayInducer)
                    self.ParamName_List.append(ParamName_DelayInducer)
                    self.Param_List.append(Param_DelayInducer)
                    self.ParamUnits_List.append(ParamUnits_DelayInducer)
                    self.VarName_List.append(VarName_DelayInducer)
                    self.y0_List.append(y0_DelayInducer)
                    self.Num_Param_List.append(len(Param_DelayInducer))
                    self.Model_List.append('Model 3.0 - DelayInducer')
                    
                    ### Model 3.1 - Inducer Transport Delay with Protein Maturation Kinetics ###
                    SystemType = 'DelayInducerKMat'
                    Param_DelayInducerKMat, SSE_DelayInducerKMat, y0_DelayInducerKMat, VarName_DelayInducerKMat, ParamName_DelayInducerKMat, ParamUnits_DelayInducerKMat \
                        = self.IndPLib1.Run_InducibleSystem(SystemType, self.n_DoseRes, self.K_ind_DoseRes, self.Data_header1, self.Data_array1, self.Inducer1, self.Inducer_log1)
                
                    self.SSE_Combined.append(SSE_DelayInducerKMat)
                    self.ParamName_List.append(ParamName_DelayInducerKMat)
                    self.Param_List.append(Param_DelayInducerKMat)
                    self.ParamUnits_List.append(ParamUnits_DelayInducerKMat)
                    self.VarName_List.append(VarName_DelayInducerKMat)
                    self.y0_List.append(y0_DelayInducerKMat)
                    self.Num_Param_List.append(len(Param_DelayInducerKMat))
                    self.Model_List.append('Model 3.1 - DelayInducerKMat')
                    
                    ## Model 3.2 - Single-ODE Inducer Transport Delay ###
                    SystemType = 'SingleDelayInducer'
                    Param_SingleDelayInducer, SSE_SingleDelayInducer, y0_SingleDelayInducer, VarName_SingleDelayInducer, ParamName_SingleDelayInducer, ParamUnits_SingleDelayInducer \
                        = self.IndPLib1.Run_InducibleSystem(SystemType, self.n_DoseRes, self.K_ind_DoseRes, self.Data_header1, self.Data_array1, self.Inducer1, self.Inducer_log1)
                
                    self.SSE_Combined.append(SSE_SingleDelayInducer)
                    self.ParamName_List.append(ParamName_SingleDelayInducer)
                    self.Param_List.append(Param_SingleDelayInducer)
                    self.ParamUnits_List.append(ParamUnits_SingleDelayInducer)
                    self.VarName_List.append(VarName_SingleDelayInducer)
                    self.y0_List.append(y0_SingleDelayInducer)
                    self.Num_Param_List.append(len(Param_SingleDelayInducer))
                    self.Model_List.append('Model 3.2 - SingleDelayInducer')
              
                # -------------------------------------------------------------------------------- #
            
                if (SystemOpt.casefold() == 'inhibition') or (SystemOpt.casefold() == 'all') or (self.OptInhibition_ == 'true'):
                
                    ## Model 9 - Inducer Inhibition at High Concentration ###
                    SystemType = 'ConstantIndInhibition'
                    Param_ConstantIndInhibition, SSE_ConstantIndInhibition, y0_ConstantIndInhibition, VarName_ConstantIndInhibition, ParamName_ConstantIndInhibition, ParamUnits_ConstantIndInhibition \
                        = self.IndPLib1.Run_InducibleSystem(SystemType, self.n_DoseRes, self.K_ind_DoseRes, self.Data_header1, self.Data_array1, self.Inducer1, self.Inducer_log1)
                
                    self.SSE_Combined.append(SSE_ConstantIndInhibition)
                    self.ParamName_List.append(ParamName_ConstantIndInhibition)
                    self.Param_List.append(Param_ConstantIndInhibition)
                    self.ParamUnits_List.append(ParamUnits_ConstantIndInhibition)
                    self.VarName_List.append(VarName_ConstantIndInhibition)
                    self.y0_List.append(y0_ConstantIndInhibition)
                    self.Num_Param_List.append(len(Param_ConstantIndInhibition))
                    self.Model_List.append('Model 4.0 - ConstantIndInhibition')
                    
                    ## Model 9.1 - Inducer with Delay response and Inhibition at High Concentration ###
                    SystemType = 'DelayIndInhibition'
                    Param_DelayIndInhibition, SSE_DelayIndInhibition, y0_DelayIndInhibition, VarName_DelayIndInhibition, ParamName_DelayIndInhibition, ParamUnits_DelayIndInhibition \
                        = self.IndPLib1.Run_InducibleSystem(SystemType, self.n_DoseRes, self.K_ind_DoseRes, self.Data_header1, self.Data_array1, self.Inducer1, self.Inducer_log1)
                
                    self.SSE_Combined.append(SSE_DelayIndInhibition)
                    self.ParamName_List.append(ParamName_DelayIndInhibition)
                    self.Param_List.append(Param_DelayIndInhibition)
                    self.ParamUnits_List.append(ParamUnits_DelayIndInhibition)
                    self.VarName_List.append(VarName_DelayIndInhibition)
                    self.y0_List.append(y0_DelayIndInhibition)
                    self.Num_Param_List.append(len(Param_DelayIndInhibition))
                    self.Model_List.append('Model 4.1 - DelayIndInhibition')
                    
                # -------------------------------------------------------------------------------- #
                    
                if (SystemOpt.casefold() == 'all'):
            
                    ## Model 7 - Inducer Delay with Degradation ###
                    SystemType = 'DelayDegradationInducer'
                    Param_DelayDegradationInducer, SSE_DelayDegradationInducer, y0_DelayDegradationInducer, VarName_DelayDegradationInducer, ParamName_DelayDegradationInducer, ParamUnits_DelayDegradationInducer \
                        = self.IndPLib1.Run_InducibleSystem(SystemType, self.n_DoseRes, self.K_ind_DoseRes, self.Data_header1, self.Data_array1, self.Inducer1, self.Inducer_log1)
                
                    self.SSE_Combined.append(SSE_DelayDegradationInducer)
                    self.ParamName_List.append(ParamName_DelayDegradationInducer)
                    self.Param_List.append(Param_DelayDegradationInducer)
                    self.ParamUnits_List.append(ParamUnits_DelayDegradationInducer)
                    self.VarName_List.append(VarName_DelayDegradationInducer)
                    self.y0_List.append(y0_DelayDegradationInducer)
                    self.Num_Param_List.append(len(Param_DelayDegradationInducer))
                    self.Model_List.append('Model 5.0 - DelayDegradationInducer')
                    
                    
                
                    ## Model 8 - Inducer Degradation with Delay ###
                    SystemType = 'DegradationDelayInducer'
                    Param_DegradationDelayInducer, SSE_DegradationDelayInducer, y0_DegradationDelayInducer, VarName_DegradationDelayInducer, ParamName_DegradationDelayInducer, ParamUnits_DegradationDelayInducer \
                        = self.IndPLib1.Run_InducibleSystem(SystemType, self.n_DoseRes, self.K_ind_DoseRes, self.Data_header1, self.Data_array1, self.Inducer1, self.Inducer_log1)
                
                    self.SSE_Combined.append(SSE_DegradationDelayInducer)
                    self.ParamName_List.append(ParamName_DegradationDelayInducer)
                    self.Param_List.append(Param_DegradationDelayInducer)
                    self.ParamUnits_List.append(ParamUnits_DegradationDelayInducer)
                    self.VarName_List.append(VarName_DegradationDelayInducer)
                    self.y0_List.append(y0_DegradationDelayInducer)
                    self.Num_Param_List.append(len(Param_DegradationDelayInducer))
                    self.Model_List.append('Model 5.1 - DegradationDelayInducer')
        
            
            # -------------------------------------------------------------------------------- #
        
    def RunModelSelection(self):
        # -------------------------------------------------------------------------------- #

        ### AIC ###
        self.AIC_Results = ComputeAIC.run_AIC(self.SSE_Combined, self.Sample_size, self.Num_Param_List)
        
        # -------------------------------------------------------------------------------- #
        
        ### Find Min SSE and AIC ###
        min_SSE = min(self.SSE_Combined)
        min_SSE_index = self.SSE_Combined.index(min_SSE)
        
        self.min_AIC = min(self.AIC_Results)
        self.min_AIC_index = self.AIC_Results.index(self.min_AIC)
        
        # -------------------------------------------------------------------------------- #
        
        ### Model Recommendation and Print out###
        
        self.Best_Model = self.Model_List[self.min_AIC_index]
        
        print('\n')
        print('Models Tested:', self.Model_List)
        print('Number of Data points:', self.Sample_size)
        print('NumParam:', self.Num_Param_List)
        print('SSE Results:', self.SSE_Combined)
        print('AIC Results:', self.AIC_Results)
        
        print('Lowest SSE Value and Index:', min_SSE, ',', min_SSE_index)
        print('Lowest AIC Value and Index:', self.min_AIC, ',', self.min_AIC_index)
        
        print('\nBased on the lowest AIC value, the recommended characterization model is:')
        print('>', self.Best_Model)
        
        # -------------------------------------------------------------------------------- #
        
        ### Plot Results of Best Model ###

        # Initiate lists to store strings
        Param_String = []
        ParamName_String = []
        ParamUnits_String = []
        VarName_String = []
        y0_String = []
        
        # store StringName and VariableName into dict for mapping purpose
        for m in range(0, len(self.Model_List)):
            Param_String.append('_'.join(('Param', self.Model_List[m].split()[3])))
            ParamName_String.append('_'.join(('ParamName', self.Model_List[m].split()[3])))
            ParamUnits_String.append('_'.join(('ParamUnits', self.Model_List[m].split()[3])))
            VarName_String.append('_'.join(('VarName', self.Model_List[m].split()[3])))
            y0_String.append('_'.join(('y0', self.Model_List[m].split()[3])))
        
        Param_dict = dict(zip(Param_String, self.Param_List))
        ParamName_dict = dict(zip(ParamName_String, self.ParamName_List))
        ParamUnits_dict = dict(zip(ParamUnits_String, self.ParamUnits_List))
        VarName_dict = dict(zip(VarName_String, self.VarName_List))
        y0_dict = dict(zip(y0_String, self.y0_List))
        
        # Get the best model and update the parameters for txt file and SBML outputs
        BestSystemType = self.Best_Model.split()[3]
        Param_ = '_'.join(('Param', BestSystemType))
        self.FittedParams = Param_dict[Param_] #convert string to variable name
        ParamName_ = '_'.join(('ParamName', BestSystemType))
        self.ParamName = ParamName_dict[ParamName_] #convert string to variable name
        ParamUnits_ = '_'.join(('ParamUnits', BestSystemType))
        self.ParamUnits = ParamUnits_dict[ParamUnits_] #convert string to variable name
        VarName_ = '_'.join(('VarName', BestSystemType))
        self.VarName = VarName_dict[VarName_] #convert string to variable name
        y0_ = '_'.join(('y0', BestSystemType))
        self.Init = y0_dict[y0_] #convert string to variable name
        
        
        self.Time, self.VariableMatrix, self.DataLegend, self.ODEstring = self.IndPLib1.Run_InduciblePlot\
                (BestSystemType, self.Init, self.Data_header1, self.Data_array1, self.Data_stddev1, self.Inducer1, self.Inducer_log1, self.FittedParams)
                
        plt.show()     
        
        # -------------------------------------------------------------------------------- #
        
        ### Calculate Time Elapsed ###
        endtime = time.time()
        self.elapsedtime = endtime - self.starttime
        print('Overall Time taken =', self.elapsedtime, 's')
        
        # -------------------------------------------------------------------------------- #

    def CreateOutputTextFile(self):
        
        ### Rank the AIC results
        Rank = ss.rankdata(self.AIC_Results)
        
        ### create a table with (Model, SSE, AIC, Rank)
        TableData = []
        Header = ["Model", "SSE", "AIC", "\u0394AIC", "Evidence", "Rank"]
        
        i = 0
        dAIC = []
        Evidence = []
        for x in self.AIC_Results:
            dAIC.append(x - self.min_AIC)
            if (dAIC[i] == 0):
                Evidence.append('-')
            elif (dAIC[i] > 0) and (dAIC[i] <= 2):
                Evidence.append('Substantial Support')
            elif (dAIC[i] > 10):
                Evidence.append('No Support')
            else:
                Evidence.append('Weak Support')
            i += 1

            
#        #decide if the model with rank 1 is a better model with confidence
#        Rank_ = [int(x) for x in Rank]
#        BestEvidence = Evidence[Rank_.index(2)]
        
        if 'Substantial Support' in Evidence:
            Count = Evidence.count('Substantial Support')
            BestEvidence = 'low confidence. There are ' + str(Count) + ' other comparably good models'
        else:
            BestEvidence = 'confidence'
        
        
        for i in range(0, len(self.Model_List)):
            TableData.append([str(self.Model_List[i]), str(self.SSE_Combined[i]), str(self.AIC_Results[i]), str(dAIC[i]), Evidence[i], str(int(Rank[i]))])
        Table = tabulate(TableData, Header, tablefmt='orgtbl')
        print(Table)
        
        ### Create and write to Txt File ###
        Txtfilename1, DateTimenow = Txtfilename.gettxtfilename()
        print('\nText File Generated:', Txtfilename)
        
        Txtpath = "Results\\" + Txtfilename1
        
        f = open(Txtpath,encoding = 'utf-8', mode ="a+")
        
        f.write('Input File name: '+self.Input_filename_+'\n')
        f.write('\n')
        f.write('Models Tested: '+str(self.Model_List)+'\n')
        f.write('\n')
        f.write(str(Table))
        f.write('\n')
        f.write('Recommended Model: '+self.Best_Model + ' with ' + BestEvidence +'\n')
        f.write('\n')
        f.write('Optimized Parameters:\n')
        
        
        for i in self.ParamName:
            f.write('\t'+ i +' = '+ str(self.FittedParams[self.ParamName.index(i)])+'\n')
        f.write('\tdeg_mRNA = 0.1386\n')
        f.write('ODE:\n')
        for j in self.VarName:
            f.write('d'+ j + 'dt'+' = '+ self.ODEstring[self.VarName.index(j)]+'\n')
        f.write('\n')
        
        f.write('Number of Data points: '+str(self.Sample_size)+'\n')
        f.write('SSE of Ideal Model: '+str(self.SSE_Combined[self.min_AIC_index])+'\n')
        f.write('AIC of Ideal Model: '+str(self.min_AIC)+'\n')
        f.write('\n')
        f.write('Time taken: '+str(self.elapsedtime)+ 's\n')
        f.write('\n')
        f.write('Date and Time: '+DateTimenow+'\n')
        f.close()
        
    def ExportModelDataFile(self):
        ### To export Model Data in CSV file
        ExportDataFile = input("Please insert 'yes/no' to export Model data file): \n")
        
        while not ((ExportDataFile.casefold() == 'yes') or (ExportDataFile.casefold() == 'no')):
             ExportDataFile = input("Error: Incorrect Choice! Please insert either yes or no only:\n")
                
                
        if ExportDataFile.casefold() == 'yes':
            CSVfileName = Txtfilename.getcsvfilename()
            VariableMatrixData = self.VariableMatrix[-1][:,:].tolist()
            
            CSVfilePath = "Results/" + CSVfileName
            with open(CSVfilePath, 'w', newline='') as csvfile:
                CF = csv.writer(csvfile, delimiter=',')
                CF.writerow(['Time(min)'] + self.DataLegend)
                for i in range(0, len(self.Time)):
                    CF.writerow([self.Time[i]] + VariableMatrixData[i])
        else:
            print('No result data file is exported\n')
            
    def ExportSBMLFile(self):
        
        ### To export SBML file in .xml
        from exportsbml import exportsbml

        ODE = self.ODEstring;
        
        print(ODE)
        
        Variable = []
        for v in self.VarName:
            Variable.append(v)  #to convert to Molar instead of mole
        
        VarInit = self.Init;
        
        #'molL-1s-1', 's-1', 'molL-1', 'dimensionless'
        
        for u in self.ParamUnits:
            if u == 'molL-1min-1':
                self.ParamUnits[self.ParamUnits.index(u)] = 'molperLmin'
            elif u == 'molL-1':
                self.ParamUnits[self.ParamUnits.index(u)] = 'molperL'
            elif u == 's-1':
                self.ParamUnits[self.ParamUnits.index(u)] = 'per_second'
            elif u == 'min-1':
                self.ParamUnits[self.ParamUnits.index(u)] = 'per_min'
            elif u == 'dimensionless':
                self.ParamUnits[self.ParamUnits.index(u)] = 'Dimension_less'
            else:
                print('Error in the defined units for parameters')
        
        ParName = self.ParamName + ['deg_mRNA', 'inducer'];
        Params = self.FittedParams.tolist() + [0.1386, 1];
        ParamsUnit = self.ParamUnits + ['per_min', 'molperL'];
        
        ODE1 = []
        
        for o in range(0, len(ODE)):
            ODE1.append(ODE[o].replace("**", "^"))
            
        print(ODE1)
        
        exportsbml(ODE1, Variable, VarInit, ParName, Params, ParamsUnit)
        
    def PlotSBOLGraphics(self):
        ### To plot and visualize the gene circuit in SBOL visual compliant diagram
        import PlotCircuit as pc
        
        PlotCircuitFile = input("Please insert 'yes/no' to set and visualize the gene circuit diagram:\n")
        
        while not ((PlotCircuitFile == 'yes') or (PlotCircuitFile == 'no')):
             PlotCircuitFile = input("Error: Incorrect Choice! Please insert either yes or no only:\n")
        
        if PlotCircuitFile == 'yes':
            Reporter = input('Please insert Reporter type (RFP/GFP/YFP/BFP/...): ')
            print('Reporter: ', Reporter)
            if Reporter == 'RFP':
                ReporterColor = "red."
            elif Reporter == 'GFP':
                ReporterColor = "green."
            elif Reporter == 'YFP':
                ReporterColor = "yellow."
            elif Reporter == 'BFP':
                ReporterColor = "blue."
            else:
                ReporterColor = "orange."
                
            Origin = []
            PlasmidNum = input('Please insert the number of plasmids: ')
            
            while not PlasmidNum.isdigit():
                PlasmidNum = input('Error: Incorrect input! Please insert the right number: \n')
            
            for pl in range(0, int(PlasmidNum)):
                Origin.append (input('Please insert the Name of Origin '+str(pl+1)+': '))
            print('Origin: ', Origin)
        
            Gene = []
            PartNum = input('Please insert the number of parts: ')
            
            while not PartNum.isdigit():
                PartNum = input('Error: Incorrect input! Please insert the right number: \n')
            
            for pt in range(0, int(PartNum)):
                Gene.append(input('Please insert the Name of gene '+str(pt+1)+': '))
            print('Gene: ', Gene)
            
            PromoterName = input('Please insert the Name of the Inducible Promoter: \n')
            
            if PlasmidNum == '1':
                if PartNum == '1':
                    
                    
                    Input = "p.black."+ PromoterName + " r.black c." + ReporterColor + Reporter+ " " + "t.black o.black." + Origin[0]
        
                    Regulations = [{'type': 'Activation', 'from_part': {'start': 7, 'end': 7},
                         'to_part': {'start': 7,'end': 7, 'fwd': True},
                         'opts': {'color': (1.00,0.50,0.00), 'linewidth': 1.5}}]
            
                elif PartNum == '2':
                    Input = "p.black r.black c.orange."+ Gene[0] + " t.black " + \
                    "p.black." + PromoterName + " r.black c." + ReporterColor + Reporter+ " " + "t.black o.black." + Origin[0]
        
                    Regulations = [{'type': 'Repression', 'from_part': {'start': 43, 'end': 43},
                         'to_part': {'start': 78,'end': 78, 'fwd': True},
                         'opts': {'color': (1.00,0.50,0.00), 'linewidth': 1.5, 'first_arc_y_offset': -3.5}}, 
                        {'type': 'Repression', 'from_part': {'start': 66, 'end': 66},
                         'to_part': {'start': 66,'end': 66, 'fwd': True},
                         'opts': {'color': (0.12,0.47,0.71), 'linewidth': 1.5, 'second_arc_y_offset': 7}}]
                else:
                    print('Error: The numbers of part and plasmid are not compatible')
                    
            elif PlasmidNum == '2':
                if PartNum == '2':
                    Input = "p.black r.black c.orange."+ Gene[0] + " t.black o.black." + Origin[0] +\
                        " =.white " + "p.black." + PromoterName + " r.black c."+ ReporterColor + Reporter + " " + "t.black o.black." + Origin[1]
        
                    Regulations = [{'type': 'Repression', 'from_part': {'start': 43, 'end': 43},
                         'to_part': {'start': 105,'end': 105, 'fwd': True},
                         'opts': {'color': (1.00,0.50,0.00), 'linewidth': 1.5, 'first_arc_y_offset': -3.5}}, 
                        {'type': 'Repression', 'from_part': {'start': 80, 'end': 80},
                         'to_part': {'start': 80,'end': 80, 'fwd': True},
                         'opts': {'color': (0.12,0.47,0.71), 'linewidth': 1.5, 'second_arc_y_offset': 7}}]
                else:
                    print('Error: The numbers of part and plasmid are not compatible')
            else:
                print('The inserted plasmid number is too high')
            
            pc.Run_PlotCircuit(Input, Regulations)
        else:
            pass
    
    ####################################################################            
    # Helper function to simplify the whole process   
    ####################################################################    
        
    def AutoRunInducibleSystem(self, Input_filename, NumDataSet, Inducer_unit, OptInhibition, SystemOpt):
        
        self.DataReader(Input_filename, NumDataSet, Inducer_unit)
        self.RunDoseResponsePrefitting(OptInhibition)
        self.RunModels(SystemOpt)
        self.RunModelSelection()
        self.CreateOutputTextFile()
        self.ExportModelDataFile()
        self.ExportSBMLFile()
        self.PlotSBOLGraphics()
        
    def __del__(self):
        class_name = self.__class__.__name__
        print(class_name, "destroyed")