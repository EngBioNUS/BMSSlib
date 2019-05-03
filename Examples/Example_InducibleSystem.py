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
    # e.g.: [Time(min) | 100uM | 50uM | 10uM | 1uM | 0uM | SD100uM | SD50uM | SD10uM | SD1uM | SD0uM ]
    # 
    # The inducer concentration in the input file must be in units of %, uM, or nM; 
    # Note: No spacing in word between concentration value and unit
    # Note: The concentrations of inducer shall be in descending order. 
    # 
    # SD refers to the standard deviation data
    #
    # NumDataSet: Number of Data Set (excluding standard deviations)
    # Inducer_unit: The unit for the inducer concentration
    # OptInhibition: True if there is Inhibitory effect at high concentrations,
    #           else None if there is only inducible trend
    # SystemOpt: Option for Inducible system ('ALL', 'ConstInd', 'DegradeInd', 'DelayInd', 'Inhibition')   
    # 
    # 'ALL': (default setting) run through all the models available in the library
    # 'ConstInd': run models under the assumption of constant inducer concentration
    # 'DegradeInd': run models under the assumption of inducer with fast degradation behavior
    # 'DelayInd': run models under the assumption of significant initial delayed response
    # 'Inhibition': run models under the assumption of inhibition at high inducer concentrations
    #
    # Note: When the Inducer_unit is in '%', which refers to mass concentration in (g/100 mL), 
    # users will be promted to insert the Inducer Molar Mass (in g/mol) at the console the unit 
    # to convert to the unit of Molar.
    #
    # Note: the SystemOpt is case-insensitive!
    #############################################################
    
    # Example file
    Input_filename = 'Inducible_pBADAra_rbspBb_MG.csv' # Arabinose (Data set 1)
    NumDataSet = 9
    Inducer_unit = '%' # Enter only (10**-6, 10**-9 or '%')
    OptInhibition = None
    SystemOpt = 'ConstInd'
    
#    Input_filename = 'pLac_IPTG_E6k.csv' 
#    NumDataSet = 9
#    Inducer_unit = 10**-6 # Enter only (10**-6, 10**-9 or '%')
#    OptInhibition = None
#    SystemOpt = 'DelayInd'
    

    
    #############################################################
    
    # To access the input file stored in the InputData folder
    Inputfile = "InputData\\" + Input_filename
    
    #############################################################
    # The sequences of key functions used 
    #############################################################
    
#    IS.DataReader(Inputfile, NumDataSet, Inducer_unit)
#    IS.RunDoseResponsePrefitting(OptInhibition)
#    IS.RunModels(SystemOpt)
#    IS.RunModelSelection()
#    IS.CreateOutputTextFile()
#    IS.ExportModelDataFile()
#    IS.ExportSBMLFile()
#    IS.PlotSBOLGraphics()
        
    #############################################################
    # Or use Helper function to autorun all functions
    #############################################################
    
    IS.AutoRunInducibleSystem(Inputfile, NumDataSet, Inducer_unit, OptInhibition, SystemOpt)
        
    sys.path.remove('../BMSSlibmod') #remove the path
    del IS
    
#enable the script to be run from command line
if __name__ == "__main__":
    main()