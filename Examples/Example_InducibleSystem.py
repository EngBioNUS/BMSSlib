# -*- coding: utf-8 -*-
"""
Created on Wed Dec 26 16:27:29 2018

@author: JingWui
"""

import sys, os

#sys.path.append("..") # adds higher directory to python modules path to access the key modules
sys.path.append("../BMSSlibmod")

from InducibleSystem_ import InducibleSystem

def main():
    
    # create an object from the class InducibleSystem
    IS = InducibleSystem()    
    
   
    #############################################################
    # User Input 
    
    # Input_filename: The input data file in .csv following the sequences below:
    # 
    # NOT gate: [Time(min) | state=0 | state=1 | SDstate=0 | SDstate=1]
    # AND gate: [Time(min) | state=00 | state=01 | state=10 | state=11 | SDstate=00 | SDstate=01 | SDstate=10 | SDstate=11]
    # OR gate: [Time(min) | state=00 | state=01 | state=10 | state=11 | SDstate=00 | SDstate=01 | SDstate=10 | SDstate=11]
    # 
    # SD refers to the standard deviation data
    #
    # NumDataSet: Number of Data Set (excluding standard deviations)
    # Inducer_unit: The unit for the inducer concentration
    # OptInhibition: True if there is Inhibitory effect at high concentrations,
    #                else None if there is only inducible trend
    # SystemOpt: Option for Inducible system ('ALL', 'ConstInd', 'DegradeInd', 'DelayInd', 'Inhibition')   
    # 
    # Note: the SystemOpt is case-insensitive!
    #############################################################
    
    # Example file
#    Input_filename = 'Inducible_pBADAra_rbspBb_MG.csv' # Arabinose (Data set 1)
#    NumDataSet = 9
#    Inducer_unit = '%' # Enter only (10**-6, 10**-9 or '%')
#    OptInhibition = None
    
    Input_filename = 'pLac_IPTG_E6k.csv' 
    NumDataSet = 9
    Inducer_unit = 10**-6 # Enter only (10**-6, 10**-9 or '%')
    OptInhibition = None
    

    
    #############################################################
    
    # To access the input file stored in the InputData folder
    Inputfile = "InputData\\" + Input_filename
    
    #############################################################
    # The sequences of key functions used 
    #############################################################
    
#    IS.DataReader(Inputfile, NumDataSet, Inducer_unit)
#    IS.RunDoseResponsePrefitting(OptInhibition)
#    IS.RunModels()
#    IS.RunModelSelection()
#    IS.CreateOutputTextFile()
#    IS.ExportModelDataFile()
#    IS.ExportSBMLFile()
#    IS.PlotSBOLGraphics()
        
    #############################################################
    # Or use Helper function to autorun all functions
    #############################################################
    
    IS.AutoRunInducibleSystem(Inputfile, NumDataSet, Inducer_unit, OptInhibition)
        
    sys.path.remove('../BMSSlibmod') #remove the path
    del IS
    
#enable the script to be run from command line
if __name__ == "__main__":
    main()