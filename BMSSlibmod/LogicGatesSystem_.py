# -*- coding: utf-8 -*-
"""
Created on Wed Dec 26 17:19:05 2018

@author: JingWui
"""

#############################################################################
# The LogicGates System class
#############################################################################

### Scripts imported ###
import Read_Data as ReadData
import ComputeAIC_ as ComputeAIC
import Txtfilename as Txtfilename
### Import Class from module
from LogicGatesLibrary_ import LogicGatesLibrary

### Python Packages Imported ###
#from matplotlib.font_manager import FontProperties
import matplotlib.pyplot as plt
import time
import csv
import scipy.stats as ss
from tabulate import tabulate



class LogicGatesSystem:
    def __init__(self):
        
        #close all figures from previous run
        plt.close('all')
        
        plt.rcdefaults()
        #font0 = FontProperties()
        #font = font0.copy()
        #font.set_family('sans-serif')
        #font.set_style('normal')
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
        plt.rcParams['legend.handletextpad'] = 0.8
        
        ### MSS Main Code ###
        self.starttime = time.time()
        
        #create the first object of InduciblePromoterLibrary class
        self.LogicGLib1 = LogicGatesLibrary()
        
        
    def DataReader(self, Input_filename, NumState):
        ### Read Data ###
        print('Read Data')
        self.Data_header1, self.Data_array1, self.Data_stddev1 = ReadData.readData(Input_filename, NumState)
        self.Input_filename_ = Input_filename
        self.NumState_ = NumState
        
        ReadData.plot_inputdata(self.Data_header1, self.Data_array1, self.Data_stddev1)
        
        ### Find Sample Size ###
        Data_array_numrows = self.Data_array1.shape[0]    # Time + All Inducers
        Data_array_numcols = self.Data_array1.shape[1]    # Number of RFP Data Per Inducer
        self.Sample_size = (Data_array_numrows - 1) * Data_array_numcols
        
    def RunModels(self, SystemOpt = 'NOT', iteration = 1):
        
        # -------------------------------------------------------------------------------- #

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
            
            if (SystemOpt.casefold() == 'not'):
        
                ### Model 1.0 - NOTgate ###
                SystemType = 'NOTgate'
                Param_NOTgate, SSE_NOTgate, y0_NOTgate, VarName_NOTgate, ParamName_NOTgate, ParamUnits_NOTgate\
                    = self.LogicGLib1.Run_LogicGatesSystem(SystemType, self.Data_header1, self.Data_array1, self.NumState_)
            
                self.SSE_Combined.append(SSE_NOTgate)
                self.ParamName_List.append(ParamName_NOTgate)
                self.Param_List.append(Param_NOTgate)
                self.ParamUnits_List.append(ParamUnits_NOTgate)
                self.VarName_List.append(VarName_NOTgate)
                self.y0_List.append(y0_NOTgate)
                self.Num_Param_List.append(len(Param_NOTgate))
                self.Model_List.append('Model 1.0 - NOTgate')
            
                ### Model 1.1 - NOTgateKMat ###
                SystemType = 'NOTgateKMat'
                Param_NOTgateKMat, SSE_NOTgateKMat, y0_NOTgateKMat, VarName_NOTgateKMat, ParamName_NOTgateKMat, ParamUnits_NOTgateKMat\
                    = self.LogicGLib1.Run_LogicGatesSystem(SystemType, self.Data_header1, self.Data_array1, self.NumState_)
            
                self.SSE_Combined.append(SSE_NOTgateKMat)
                self.ParamName_List.append(ParamName_NOTgateKMat)
                self.Param_List.append(Param_NOTgateKMat)
                self.ParamUnits_List.append(ParamUnits_NOTgateKMat)
                self.VarName_List.append(VarName_NOTgateKMat)
                self.y0_List.append(y0_NOTgateKMat)
                self.Num_Param_List.append(len(Param_NOTgateKMat))
                self.Model_List.append('Model 1.1 - NOTgateKMat')
            
            
                # -------------------------------------------------------------------------------- #
            
                ### Model 1.2 - NOTgateSingle ###
                SystemType = 'NOTgateSingle'
                Param_NOTgateSingle, SSE_NOTgateSingle, y0_NOTgateSingle, VarName_NOTgateSingle, ParamName_NOTgateSingle, ParamUnits_NOTgateSingle\
                    = self.LogicGLib1.Run_LogicGatesSystem(SystemType, self.Data_header1, self.Data_array1, self.NumState_)
            
                self.SSE_Combined.append(SSE_NOTgateSingle)
                self.ParamName_List.append(ParamName_NOTgateSingle)
                self.Param_List.append(Param_NOTgateSingle)
                self.ParamUnits_List.append(ParamUnits_NOTgateSingle)
                self.VarName_List.append(VarName_NOTgateSingle)
                self.y0_List.append(y0_NOTgateSingle)
                self.Num_Param_List.append(len(Param_NOTgateSingle))
                self.Model_List.append('Model 1.2 - NOTgateSingle')
            
                ### Model 1.3 - NOTgateSingleKMat ###
                SystemType = 'NOTgateSingleKMat'
                Param_NOTgateSingleKMat, SSE_NOTgateSingleKMat, y0_NOTgateSingleKMat, VarName_NOTgateSingleKMat, ParamName_NOTgateSingleKMat, ParamUnits_NOTgateSingleKMat\
                    = self.LogicGLib1.Run_LogicGatesSystem(SystemType, self.Data_header1, self.Data_array1, self.NumState_)
            
                self.SSE_Combined.append(SSE_NOTgateSingleKMat)
                self.ParamName_List.append(ParamName_NOTgateSingleKMat)
                self.Param_List.append(Param_NOTgateSingleKMat)
                self.ParamUnits_List.append(ParamUnits_NOTgateSingleKMat)
                self.VarName_List.append(VarName_NOTgateSingleKMat)
                self.y0_List.append(y0_NOTgateSingleKMat)
                self.Num_Param_List.append(len(Param_NOTgateSingleKMat))
                self.Model_List.append('Model 1.3 - NOTgateSingleKMat')
            
        
            # -------------------------------------------------------------------------------- #
            
            elif SystemOpt.casefold() == 'and':
                ### Model 2.0 - ANDgate ###
                SystemType = 'ANDgate'
                Param_ANDgate, SSE_ANDgate, y0_ANDgate, VarName_ANDgate, ParamName_ANDgate, ParamUnits_ANDgate \
                    = self.LogicGLib1.Run_LogicGatesSystem(SystemType, self.Data_header1, self.Data_array1, self.NumState_)
            
                self.SSE_Combined.append(SSE_ANDgate)
                self.ParamName_List.append(ParamName_ANDgate)
                self.Param_List.append(Param_ANDgate)
                self.ParamUnits_List.append(ParamUnits_ANDgate)
                self.VarName_List.append(VarName_ANDgate)
                self.y0_List.append(y0_ANDgate)
                self.Num_Param_List.append(len(Param_ANDgate))
                self.Model_List.append('Model 2.0 - ANDgate')
            
                ### Model 2.1 - ANDgateBLeak1 ###
                SystemType = 'ANDgateBLeak1'
                Param_ANDgateBLeak1, SSE_ANDgateBLeak1, y0_ANDgateBLeak1, VarName_ANDgateBLeak1, ParamName_ANDgateBLeak1, ParamUnits_ANDgateBLeak1 \
                    = self.LogicGLib1.Run_LogicGatesSystem(SystemType, self.Data_header1, self.Data_array1, self.NumState_)
            
                self.SSE_Combined.append(SSE_ANDgateBLeak1)
                self.ParamName_List.append(ParamName_ANDgateBLeak1)
                self.Param_List.append(Param_ANDgateBLeak1)
                self.ParamUnits_List.append(ParamUnits_ANDgateBLeak1)
                self.VarName_List.append(VarName_ANDgateBLeak1)
                self.y0_List.append(y0_ANDgateBLeak1)
                self.Num_Param_List.append(len(Param_ANDgateBLeak1))
                self.Model_List.append('Model 2.1 - ANDgateBLeak1')
            
                ### Model 2.2 - ANDgateBLeak2 ###
                SystemType = 'ANDgateBLeak2'
                Param_ANDgateBLeak2, SSE_ANDgateBLeak2, y0_ANDgateBLeak2, VarName_ANDgateBLeak2, ParamName_ANDgateBLeak2, ParamUnits_ANDgateBLeak2 \
                    = self.LogicGLib1.Run_LogicGatesSystem(SystemType, self.Data_header1, self.Data_array1, self.NumState_)
            
                self.SSE_Combined.append(SSE_ANDgateBLeak2)
                self.ParamName_List.append(ParamName_ANDgateBLeak2)
                self.Param_List.append(Param_ANDgateBLeak2)
                self.ParamUnits_List.append(ParamUnits_ANDgateBLeak2)
                self.VarName_List.append(VarName_ANDgateBLeak2)
                self.y0_List.append(y0_ANDgateBLeak2)
                self.Num_Param_List.append(len(Param_ANDgateBLeak2))
                self.Model_List.append('Model 2.2 - ANDgateBLeak2')
            
                ### Model 2.3 - ANDgateBLeak3 ###
                SystemType = 'ANDgateBLeak3'
                Param_ANDgateBLeak3, SSE_ANDgateBLeak3, y0_ANDgateBLeak3, VarName_ANDgateBLeak3, ParamName_ANDgateBLeak3, ParamUnits_ANDgateBLeak3 \
                    = self.LogicGLib1.Run_LogicGatesSystem(SystemType, self.Data_header1, self.Data_array1, self.NumState_)
            
                self.SSE_Combined.append(SSE_ANDgateBLeak3)
                self.ParamName_List.append(ParamName_ANDgateBLeak3)
                self.Param_List.append(Param_ANDgateBLeak3)
                self.ParamUnits_List.append(ParamUnits_ANDgateBLeak3)
                self.VarName_List.append(VarName_ANDgateBLeak3)
                self.y0_List.append(y0_ANDgateBLeak3)
                self.Num_Param_List.append(len(Param_ANDgateBLeak3))
                self.Model_List.append('Model 2.3 - ANDgateBLeak3')
            
                ## Model 2.4 - ANDgateBLeak13 ###
                SystemType = 'ANDgateBLeak13'
                Param_ANDgateBLeak13, SSE_ANDgateBLeak13, y0_ANDgateBLeak13, VarName_ANDgateBLeak13, ParamName_ANDgateBLeak13, ParamUnits_ANDgateBLeak13 \
                    = self.LogicGLib1.Run_LogicGatesSystem(SystemType, self.Data_header1, self.Data_array1, self.NumState_)
            
                self.SSE_Combined.append(SSE_ANDgateBLeak13)
                self.ParamName_List.append(ParamName_ANDgateBLeak13)
                self.Param_List.append(Param_ANDgateBLeak13)
                self.ParamUnits_List.append(ParamUnits_ANDgateBLeak13)
                self.VarName_List.append(VarName_ANDgateBLeak13)
                self.y0_List.append(y0_ANDgateBLeak13)
                self.Num_Param_List.append(len(Param_ANDgateBLeak13))
                self.Model_List.append('Model 2.4 - ANDgateBLeak13')
            
                ### Model 2.5 - ANDgateBLeak13KMat ###
                SystemType = 'ANDgateBLeak13KMat'
                Param_ANDgateBLeak13KMat, SSE_ANDgateBLeak13KMat, y0_ANDgateBLeak13KMat, VarName_ANDgateBLeak13KMat, ParamName_ANDgateBLeak13KMat, ParamUnits_ANDgateBLeak13KMat \
                    = self.LogicGLib1.Run_LogicGatesSystem(SystemType, self.Data_header1, self.Data_array1, self.NumState_)
            
                self.SSE_Combined.append(SSE_ANDgateBLeak13KMat)
                self.ParamName_List.append(ParamName_ANDgateBLeak13KMat)
                self.Param_List.append(Param_ANDgateBLeak13KMat)
                self.ParamUnits_List.append(ParamUnits_ANDgateBLeak13KMat)
                self.VarName_List.append(VarName_ANDgateBLeak13KMat)
                self.y0_List.append(y0_ANDgateBLeak13KMat)
                self.Num_Param_List.append(len(Param_ANDgateBLeak13KMat))
                self.Model_List.append('Model 2.5 - ANDgateBLeak13KMat')
        
            # -------------------------------------------------------------------------------- #
            elif (SystemOpt.casefold() == 'or'):
                
                ### Model 3.0 - ORgate ###
                SystemType = 'ORgate'
                Param_ORgate, SSE_ORgate, y0_ORgate, VarName_ORgate, ParamName_ORgate, ParamUnits_ORgate \
                    = self.LogicGLib1.Run_LogicGatesSystem(SystemType, self.Data_header1, self.Data_array1, self.NumState_)
            
                self.SSE_Combined.append(SSE_ORgate)
                self.ParamName_List.append(ParamName_ORgate)
                self.Param_List.append(Param_ORgate)
                self.ParamUnits_List.append(ParamUnits_ORgate)
                self.VarName_List.append(VarName_ORgate)
                self.y0_List.append(y0_ORgate)
                self.Num_Param_List.append(len(Param_ORgate))
                self.Model_List.append('Model 3.0 - ORgate')
            
                ### Model 3.1 - ORgate ###
                SystemType = 'ORgateDelay'
                Param_ORgateDelay, SSE_ORgateDelay, y0_ORgateDelay, VarName_ORgateDelay, ParamName_ORgateDelay, ParamUnits_ORgateDelay \
                    = self.LogicGLib1.Run_LogicGatesSystem(SystemType, self.Data_header1, self.Data_array1, self.NumState_)
            
                self.SSE_Combined.append(SSE_ORgateDelay)
                self.ParamName_List.append(ParamName_ORgateDelay)
                self.Param_List.append(Param_ORgateDelay)
                self.ParamUnits_List.append(ParamUnits_ORgateDelay)
                self.VarName_List.append(VarName_ORgateDelay)
                self.y0_List.append(y0_ORgateDelay)
                self.Num_Param_List.append(len(Param_ORgateDelay))
                self.Model_List.append('Model 3.1 - ORgateDelay')
            
                ### Model 3.2 - ORgate_Delay ###
                SystemType = 'ORgate_Delay'
                Param_ORgate_Delay, SSE_ORgate_Delay, y0_ORgate_Delay, VarName_ORgate_Delay, ParamName_ORgate_Delay, ParamUnits_ORgate_Delay \
                    = self.LogicGLib1.Run_LogicGatesSystem(SystemType, self.Data_header1, self.Data_array1, self.NumState_)
            
                self.SSE_Combined.append(SSE_ORgate_Delay)
                self.ParamName_List.append(ParamName_ORgate_Delay)
                self.Param_List.append(Param_ORgate_Delay)
                self.ParamUnits_List.append(ParamUnits_ORgate_Delay)
                self.VarName_List.append(VarName_ORgate_Delay)
                self.y0_List.append(y0_ORgate_Delay)
                self.Num_Param_List.append(len(Param_ORgate_Delay))
                self.Model_List.append('Model 3.2 - ORgate_Delay')
            
                ### Model 3.3 - ORgateDegradation ###
                SystemType = 'ORgateDegradation'
                Param_ORgateDegradation, SSE_ORgateDegradation, y0_ORgateDegradation, VarName_ORgateDegradation, ParamName_ORgateDegradation, ParamUnits_ORgateDegradation \
                    = self.LogicGLib1.Run_LogicGatesSystem(SystemType, self.Data_header1, self.Data_array1, self.NumState_)
            
                self.SSE_Combined.append(SSE_ORgateDegradation)
                self.ParamName_List.append(ParamName_ORgateDegradation)
                self.Param_List.append(Param_ORgateDegradation)
                self.ParamUnits_List.append(ParamUnits_ORgateDegradation)
                self.VarName_List.append(VarName_ORgateDegradation)
                self.y0_List.append(y0_ORgateDegradation)
                self.Num_Param_List.append(len(Param_ORgateDegradation))
                self.Model_List.append('Model 3.3 - ORgateDegradation')
            
                ### Model 3.4 - ORgate_Degradation ###
                SystemType = 'ORgate_Degradation'
                Param_ORgate_Degradation, SSE_ORgate_Degradation, y0_ORgate_Degradation, VarName_ORgate_Degradation, ParamName_ORgate_Degradation, ParamUnits_ORgate_Degradation \
                    = self.LogicGLib1.Run_LogicGatesSystem(SystemType, self.Data_header1, self.Data_array1, self.NumState_)
            
                self.SSE_Combined.append(SSE_ORgate_Degradation)
                self.ParamName_List.append(ParamName_ORgate_Degradation)
                self.Param_List.append(Param_ORgate_Degradation)
                self.ParamUnits_List.append(ParamUnits_ORgate_Degradation)
                self.VarName_List.append(VarName_ORgate_Degradation)
                self.y0_List.append(y0_ORgate_Degradation)
                self.Num_Param_List.append(len(Param_ORgate_Degradation))
                self.Model_List.append('Model 3.4 - ORgate_Degradation')
            
                # -------------------------------------------------------------------------------- #
            
                ### Model 3.5 - ORgateDelayDegradation ###
                SystemType = 'ORgateDelayDegradation'
                Param_ORgateDelayDegradation, SSE_ORgateDelayDegradation, y0_ORgateDelayDegradation, \
                VarName_ORgateDelayDegradation, ParamName_ORgateDelayDegradation, ParamUnits_ORgateDelayDegradation \
                    = self.LogicGLib1.Run_LogicGatesSystem(SystemType, self.Data_header1, self.Data_array1, self.NumState_)
            
                self.SSE_Combined.append(SSE_ORgateDelayDegradation)
                self.ParamName_List.append(ParamName_ORgateDelayDegradation)
                self.Param_List.append(Param_ORgateDelayDegradation)
                self.ParamUnits_List.append(ParamUnits_ORgateDelayDegradation)
                self.VarName_List.append(VarName_ORgateDelayDegradation)
                self.y0_List.append(y0_ORgateDelayDegradation)
                self.Num_Param_List.append(len(Param_ORgateDelayDegradation))
                self.Model_List.append('Model 3.5 - ORgateDelayDegradation')
            
                ### Model 3.6 - ORgateDegradationDelay ###
                SystemType = 'ORgateDegradationDelay'
                Param_ORgateDegradationDelay, SSE_ORgateDegradationDelay, y0_ORgateDegradationDelay, \
                VarName_ORgateDegradationDelay, ParamName_ORgateDegradationDelay, ParamUnits_ORgateDegradationDelay \
                    = self.LogicGLib1.Run_LogicGatesSystem(SystemType, self.Data_header1, self.Data_array1, self.NumState_)
            
                self.SSE_Combined.append(SSE_ORgateDegradationDelay)
                self.ParamName_List.append(ParamName_ORgateDegradationDelay)
                self.Param_List.append(Param_ORgateDegradationDelay)
                self.ParamUnits_List.append(ParamUnits_ORgateDegradationDelay)
                self.VarName_List.append(VarName_ORgateDegradationDelay)
                self.y0_List.append(y0_ORgateDegradationDelay)
                self.Num_Param_List.append(len(Param_ORgateDegradationDelay))
                self.Model_List.append('Model 3.6 - ORgateDegradationDelay')
            
                ### Model 3.7 - ORgateDelayDelay ###
                SystemType = 'ORgateDelayDelay'
                Param_ORgateDelayDelay, SSE_ORgateDelayDelay, y0_ORgateDelayDelay, \
                VarName_ORgateDelayDelay, ParamName_ORgateDelayDelay, ParamUnits_ORgateDelayDelay \
                    = self.LogicGLib1.Run_LogicGatesSystem(SystemType, self.Data_header1, self.Data_array1, self.NumState_)
            
                self.SSE_Combined.append(SSE_ORgateDelayDelay)
                self.ParamName_List.append(ParamName_ORgateDelayDelay)
                self.Param_List.append(Param_ORgateDelayDelay)
                self.ParamUnits_List.append(ParamUnits_ORgateDelayDelay)
                self.VarName_List.append(VarName_ORgateDelayDelay)
                self.y0_List.append(y0_ORgateDelayDelay)
                self.Num_Param_List.append(len(Param_ORgateDelayDelay))
                self.Model_List.append('Model 3.7 - ORgateDelayDelay')
            
                #-------------------------------------------------------------------------------- #
            
                ### Model 3.8- OR with Delay Degradation (with resource competition at State11) System ###
                SystemType = 'ORgateDelayDegradeResCompete'
                Param_ORgateDelayDegradeResCompete, SSE_ORgateDelayDegradeResCompete, y0_ORgateDelayDegradeResCompete, \
                VarName_ORgateDelayDegradeResCompete, ParamName_ORgateDelayDegradeResCompete, ParamUnits_ORgateDelayDegradeResCompete \
                    = self.LogicGLib1.Run_LogicGatesSystem(SystemType, self.Data_header1, self.Data_array1, self.NumState_)
            
                self.SSE_Combined.append(SSE_ORgateDelayDegradeResCompete)
                self.ParamName_List.append(ParamName_ORgateDelayDegradeResCompete)
                self.Param_List.append(Param_ORgateDelayDegradeResCompete)
                self.ParamUnits_List.append(ParamUnits_ORgateDelayDegradeResCompete)
                self.VarName_List.append(VarName_ORgateDelayDegradeResCompete)
                self.y0_List.append(y0_ORgateDelayDegradeResCompete)
                self.Num_Param_List.append(len(Param_ORgateDelayDegradeResCompete))
                self.Model_List.append('Model 3.8 - ORgateDelayDegradeResCompete')
            
                ### Model 3.9- OR with Degradation Delay (with resource competition at State11) System ###
                SystemType = 'ORgateDegradeDelayResCompete'
                Param_ORgateDegradeDelayResCompete, SSE_ORgateDegradeDelayResCompete, y0_ORgateDegradeDelayResCompete, \
                VarName_ORgateDegradeDelayResCompete, ParamName_ORgateDegradeDelayResCompete, ParamUnits_ORgateDegradeDelayResCompete \
                    = self.LogicGLib1.Run_LogicGatesSystem(SystemType, self.Data_header1, self.Data_array1, self.NumState_)
            
                self.SSE_Combined.append(SSE_ORgateDegradeDelayResCompete)
                self.ParamName_List.append(ParamName_ORgateDegradeDelayResCompete)
                self.Param_List.append(Param_ORgateDegradeDelayResCompete)
                self.ParamUnits_List.append(ParamUnits_ORgateDegradeDelayResCompete)
                self.VarName_List.append(VarName_ORgateDegradeDelayResCompete)
                self.y0_List.append(y0_ORgateDegradeDelayResCompete)
                self.Num_Param_List.append(len(Param_ORgateDegradeDelayResCompete))
                self.Model_List.append('Model 3.9 - ORgateDegradeDelayResCompete')
            
                ### Model 3.10- OR with Delay Delay (with resource competition at State11) System ###
                SystemType = 'ORgateDelayDelayResCompete'
                Param_ORgateDelayDelayResCompete, SSE_ORgateDelayDelayResCompete, y0_ORgateDelayDelayResCompete, \
                VarName_ORgateDelayDelayResCompete, ParamName_ORgateDelayDelayResCompete, ParamUnits_ORgateDelayDelayResCompete \
                    = self.LogicGLib1.Run_LogicGatesSystem(SystemType, self.Data_header1, self.Data_array1, self.NumState_)
            
                self.SSE_Combined.append(SSE_ORgateDelayDelayResCompete)
                self.ParamName_List.append(ParamName_ORgateDelayDelayResCompete)
                self.Param_List.append(Param_ORgateDelayDelayResCompete)
                self.ParamUnits_List.append(ParamUnits_ORgateDelayDelayResCompete)
                self.VarName_List.append(VarName_ORgateDelayDelayResCompete)
                self.y0_List.append(y0_ORgateDelayDelayResCompete)
                self.Num_Param_List.append(len(Param_ORgateDelayDelayResCompete))
                self.Model_List.append('Model 3.10 - ORgateDelayDelayResCompete')
        
            else: 
                print('Error in the selected logic gate system')
        
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
        
        
        self.Time, self.VariableMatrix, self.DataLegend, self.ODEstring = self.LogicGLib1.Run_LogicGatesPlot\
                (BestSystemType, self.Init, self.Data_header1, self.Data_array1, self.Data_stddev1, self.NumState_, self.FittedParams)
        
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
        
#        #use list comprehension
#        dAIC = [x - self.min_AIC for x in self.AIC_Results]
        
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
            BestEvidence = 'low confidence. There are ' + str(Count) + 'other comparably good models'
        else:
            BestEvidence = 'confidence'
        
        for i in range(0, len(self.Model_List)):
            TableData.append([str(self.Model_List[i]), str(self.SSE_Combined[i]), str(self.AIC_Results[i]), str(dAIC[i]), Evidence[i], str(int(Rank[i]))])
        Table = tabulate(TableData, Header, tablefmt='orgtbl')
        print(Table)
        
        ### Create and write to Txt File ###
        Txtfilename1, DateTimenow = Txtfilename.gettxtfilename()
        print('\nText File Generated:', Txtfilename)
        
        #f = open(Txtfilename1,"a+")
        Txtpath = "Results\\" + Txtfilename1
        
        f = open(Txtpath,encoding = 'utf-8', mode ="a+")
        
        f.write('Input File name: '+self.Input_filename_+'\n')
        f.write('\n')
        f.write('Models Tested: '+str(self.Model_List)+'\n')
        f.write('\n')
        f.write(str(Table))
        f.write('\n\n')
        f.write('Recommended Model: '+self.Best_Model+ ' with ' + BestEvidence +'\n')
        f.write('\n')
        f.write('Optimized Parameters:\n')
        
        
        for i in self.ParamName:
            f.write('\t'+ i +' = '+ str(self.FittedParams[self.ParamName.index(i)])+'\n')
        f.write('\tdeg_mRNA1 = 0.1386\n')
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
        ExportDataFile = input("Please insert 'yes/no' to export Model data file):")
    
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
        
        if self.NumState_ == 2:
            ParName = self.ParamName + ['deg_mRNA', 'state'];
            Params = self.FittedParams.tolist() + [0.1386, 1];
            ParamsUnit = self.ParamUnits + ['per_min', 'Dimension_less'];
        elif self.NumState_ == 4:
            ParName = self.ParamName + ['deg_mRNA', 'state1', 'state2'];
            Params = self.FittedParams.tolist() + [0.1386, 1, 1];
            ParamsUnit = self.ParamUnits + ['per_min', 'Dimension_less', 'Dimension_less'];
        else:
            print('Error in setting NumState for SBML output')
        
        exportsbml(ODE, Variable, VarInit, ParName, Params, ParamsUnit)
        
    def PlotSBOLGraphics(self):
        
        ### To plot and visualize the gene circuit in SBOL visual compliant diagram
        import PlotCircuit as pc
        
        
        PlotCircuitFile = input("Please insert 'yes/no' to set and visualize the gene circuit diagram:\n")
        
        while not ((PlotCircuitFile == 'yes') or (PlotCircuitFile == 'no')):
             PlotCircuitFile = input("Error: Incorrect Choice! Please insert either yes or no only:\n")
            
        
        if PlotCircuitFile == 'yes':
            Reporter = input('Please insert Reporter type (RFP/GFP/YFP/BFP/...): \n')
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
            PlasmidNum = input('Please insert the number of plasmids: \n')
            
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
        
            # NOT gate models
            if 'NOTgate' in self.Best_Model:
                if PlasmidNum == '1':
                    if PartNum == '1':
                        Input = "p.black r.black c." + ReporterColor + Reporter + " " + "t.black o.black." + Origin[0]
        
                        Regulations = [{'type': 'Repression', 'from_part': {'start': 8, 'end': 8},
                             'to_part': {'start': 8,'end': 8, 'fwd': True},
                             'opts': {'color': (0.12,0.47,0.71), 'linewidth': 1.5}}]
                    elif PartNum == '2':
        
        
                        Input = "p.black r.black c.blue."+ Gene[0] + " t.black "+ "p.black r.black c."+ ReporterColor + Reporter + " " + "t.black o.black." + Origin[0]
        
                        Regulations = [{'type': 'Repression', 'from_part': {'start': 40, 'end': 40},
                             'to_part': {'start': 78,'end': 78, 'fwd': True},
                             'opts': {'color': (0.12,0.47,0.71), 'linewidth': 1.5}}]
                    else:
                        print('Error: The numbers of part and plasmid are not compatible')
        
                elif PlasmidNum == '2':
                    if PartNum == '2':
                        Input = "p.black r.black c.blue."+ Gene[0] + " t.black o.black." + Origin[0] + " "+ "=.white "+ "p.black r.black c."+ ReporterColor + Reporter + " " + "t.black o.black." + Origin[1]
        
                        Regulations = [{'type': 'Repression', 'from_part': {'start': 43, 'end': 43},
                             'to_part': {'start': 105,'end': 105, 'fwd': True},
                             'opts': {'color': (0.12,0.47,0.71), 'linewidth': 1.5}}]
                    elif PartNum == '3':
                        Input = "p.black r.black c.orange."+ Gene[0] + " t.black o.black." + Origin[0] + " =.white " + \
                        "p.black r.black c.blue."+ Gene[1] + " t.black " + \
                        "p.black r.black c."+ ReporterColor + Reporter + " " + "t.black o.black." + Origin[1]
        
                        Regulations = [{'type': 'Connection', 'from_part': {'start': 43, 'end': 43},
                             'to_part': {'start': 142,'end': 142, 'fwd': True},
                             'opts': {'color': (1.00,0.50,0.00), 'linewidth': 1.5}},
                            {'type': 'Repression', 'from_part': {'start': 142, 'end': 142},
                             'to_part': {'start': 178,'end': 178, 'fwd': True},
                             'opts': {'color': (0.12,0.47,0.71), 'linewidth': 1.5}}]
                    else:
                        print('Error: The numbers of part and plasmid are not compatible')
                elif PlasmidNum == '3':
                    if PartNum == '3':
                        Input = "p.black r.black c.orange."+ Gene[0] + " t.black o.black." + Origin[0] + " =.white " + \
                        "p.black r.black c.blue."+ Gene[1] + " t.black o.black." + Origin[1] + " =.white " + \
                        "p.black r.black c."+ ReporterColor + Reporter + " " + "t.black o.black." + Origin[2]
        
                        Regulations = [{'type': 'Connection', 'from_part': {'start': 43, 'end': 43},
                             'to_part': {'start': 142,'end': 142, 'fwd': True},
                             'opts': {'color': (1.00,0.50,0.00), 'linewidth': 1.5}},
                            {'type': 'Repression', 'from_part': {'start': 142, 'end': 142},
                             'to_part': {'start': 205,'end': 205, 'fwd': True},
                             'opts': {'color': (0.12,0.47,0.71), 'linewidth': 1.5}}]
        
                    else:
                        print('Error: The numbers of part and plasmid are not compatible')
                else:
                    print('The inserted plasmid number is too high')
        
            # AND gate models
            elif 'ANDgate' in self.Best_Model:
                if PlasmidNum == '1':
                    if PartNum == '1':
                        Input = "p.black r.black c." + ReporterColor + Reporter + " " + "t.black o.black." + Origin[0]
        
                        Regulations = [{'type': 'Activation', 'from_part': {'start': 3, 'end': 3},
                         'to_part': {'start': 3,'end': 3, 'fwd': True},
                         'opts': {'color': (0.12,0.47,0.71), 'linewidth': 1.5}},
                        {'type': 'Activation', 'from_part': {'start': 11, 'end': 11},
                         'to_part': {'start': 11,'end': 11, 'fwd': True},
                         'opts': {'color': (1.00,0.50,0.00), 'linewidth': 1.5}}]
        
                    elif PartNum == '2':
                        Input = "p.black r.black c.orange." + Gene[0] + " t.black " + \
                        "p.black r.black c." + ReporterColor + Reporter + " " + "t.black o.black." + Origin[0]
        
                        Regulations = [{'type': 'Activation', 'from_part': {'start': 75, 'end': 75},
                         'to_part': {'start': 75,'end': 75, 'fwd': True},
                         'opts': {'color': (0.12,0.47,0.71), 'linewidth': 1.5}},
                        {'type': 'Activation', 'from_part': {'start': 43, 'end': 43},
                         'to_part': {'start': 83,'end': 83, 'fwd': True},
                         'opts': {'color': (1.00,0.50,0.00), 'linewidth': 1.5}}]
        
                    elif PartNum == '3':
                        Input = "p.black r.black c.orange." + Gene[0] + " t.black " + \
                        "p.black r.black c.blue." + Gene[1] + " t.black " + \
                        "p.black r.black c." + ReporterColor + Reporter + " " + "t.black o.black." + Origin[0]
        
                        Regulations = [{'type': 'Activation', 'from_part': {'start': 43, 'end': 43},
                         'to_part': {'start': 155,'end': 155, 'fwd': True},
                         'opts': {'color': (1.00,0.50,0.00), 'linewidth': 1.5}},
                         {'type': 'Activation', 'from_part': {'start': 115, 'end': 115},
                         'to_part': {'start': 147,'end': 147, 'fwd': True},
                         'opts': {'color': (0.12,0.47,0.71), 'linewidth': 1.5}}]
        
                    else:
                        print('Error: The numbers of part and plasmid are not compatible')
        
                elif PlasmidNum == '2':
                    if PartNum == '2':
                        Input = "p.black r.black c.orange." + Gene[0] + " t.black o.black." + Origin[0] + " =.white "+ \
                        "p.black r.black c."+ ReporterColor + Reporter + " t.black o.black." + Origin[1]
        
                        Regulations = [{'type': 'Activation', 'from_part': {'start': 43, 'end': 43},
                         'to_part': {'start': 110,'end': 110, 'fwd': True},
                         'opts': {'color': (1.00,0.50,0.00), 'linewidth': 1.5}},
                         {'type': 'Activation', 'from_part': {'start': 102, 'end': 102},
                         'to_part': {'start': 102,'end': 102, 'fwd': True},
                         'opts': {'color': (0.12,0.47,0.71), 'linewidth': 1.5}}]
        
                    elif PartNum == '3':
                        Input = "p.black r.black c.orange." + Gene[0] + " t.black p.black r.black c.blue." + Gene[1] +\
                        " t.black o.black." + Origin[0] + " =.white "+ \
                        "p.black r.black c."+ ReporterColor + Reporter + " t.black o.black." + Origin[1]
        
                        Regulations = [{'type': 'Activation', 'from_part': {'start': 43, 'end': 43},
                         'to_part': {'start': 182,'end': 182, 'fwd': True},
                         'opts': {'color': (1.00,0.50,0.00), 'linewidth': 1.5}},
                         {'type': 'Activation', 'from_part': {'start': 115, 'end': 115},
                         'to_part': {'start': 174,'end': 174, 'fwd': True},
                         'opts': {'color': (0.12,0.47,0.71), 'linewidth': 1.5}}]
        
                    else:
                        print('Error: The numbers of part and plasmid are not compatible')
        
                elif PlasmidNum == '3' :
                    if PartNum == '3':
                        Input = "p.black r.black c.orange." + Gene[0] + " t.black o.black." + Origin[0] + " =.white " \
                        "p.black r.black c.blue." + Gene[1] + " t.black o.black." + Origin[1] + " =.white " \
                        "p.black r.black c." + ReporterColor + Reporter + " t.black o.black." + Origin[2]
        
                        Regulations = [{'type': 'Activation', 'from_part': {'start': 44, 'end': 44},
                         'to_part': {'start': 209,'end': 209, 'fwd': True},
                         'opts': {'color': (1.00,0.50,0.00), 'linewidth': 1.5}},
                         {'type': 'Activation', 'from_part': {'start': 143, 'end': 143},
                         'to_part': {'start': 201,'end': 201, 'fwd': True},
                         'opts': {'color': (0.12,0.47,0.71), 'linewidth': 1.5}}]
                    else:
                        print('Error: The numbers of part and plasmid are not compatible')
        
                else:
                    print('Error: The inserted plasmid number is too high')
        
            elif 'ORgate' in self.Best_Model:
                if PlasmidNum == '1':
                    if PartNum == '1':
                        Input = "p.black p.black r.black c." + ReporterColor + Reporter + " t.black o.black." + Origin[0]
        
                        Regulations = None
        
                    elif PartNum == '2':
                        Input = "p.black r.black c." + ReporterColor + Reporter + " t.black " + \
                        "p.black r.black c." + ReporterColor + Reporter + " t.black o.black." + Origin[0]
        
                        Regulations = None
        
                    else:
                        print('Error: The numbers of part and plasmid are not compatible')
                elif PlasmidNum == '2':
                    if PartNum == '2':
                        Input = "p.black r.black c."+ ReporterColor + Reporter + " t.black o.black." + Origin[0] + \
                        " =.white "+ "p.black r.black c."+ ReporterColor + Reporter + " t.black o.black." + Origin[1]
        
                        Regulations = None
        
                    else:
                        print('Error: The numbers of part and plasmid are not compatible')
                else:
                    print('Error: The inserted plasmid number is too high')
        
            else:
                print('Error in the selected Best_Model')
        
            pc.Run_PlotCircuit(Input, Regulations)
        else:
            pass
        
    ####################################################################            
    # Helper function to simplify the whole process   
    ####################################################################    
        
    def AutoRunLogicGatesSystem(self, Inputfile, SystemOpt, NumState):
        
        self.DataReader(Inputfile, NumState)
        self.RunModels(SystemOpt)
        self.RunModelSelection()
        self.CreateOutputTextFile()
        self.ExportModelDataFile()
        self.ExportSBMLFile()
        self.PlotSBOLGraphics()
        
    def __del__(self):
        class_name = self.__class__.__name__
        print(class_name, "destroyed")