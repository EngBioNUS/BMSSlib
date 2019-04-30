# -*- coding: utf-8 -*-
"""
Created on Wed Dec 26 18:07:21 2018

@author: JingWui
"""

import sys, os

print(sys.path)

# adds higher level directory to python modules path to access the key modules
sys.path.append("../BMSSlibmod") 

from LogicGatesSystem_ import LogicGatesSystem

def main():
    
    # create an object from the class LogicGatesSystem
    LGS = LogicGatesSystem()  
    
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
    # NumState: Number of Data Set (excluding standard deviations)
    # SystemOpt: Option for Logic gate system ('NOT', 'AND', 'OR')   
    # 
    # Note: the SystemOpt is case-insensitive!
    #############################################################
    
    Input_filename = 'LogicGate_NOTTop10d30LB_1.csv' #NOT gate (Data set 1)
    NumState = 2
    SystemOpt = 'NOT' 
    
    #############################################################
    
    # To access the input file stored in the InputData folder
    Inputfile = "InputData\\" + Input_filename
    
    #############################################################
    # The sequences of key functions used 
    #############################################################
    
#    LGS.DataReader(Inputfile, NumState)
#    LGS.RunModels(SystemOpt)
#    LGS.RunModelSelection()
#    LGS.CreateOutputTextFile()
#    LGS.ExportModelDataFile()
#    LGS.ExportSBMLFile()
#    LGS.PlotSBOLGraphics()
    
    #############################################################
    # Or use Helper function to autorun all functions
    #############################################################
    
    LGS.AutoRunLogicGatesSystem(Inputfile, SystemOpt, NumState)
        
    sys.path.remove('../BMSSlibmod') #remove the path
    del LGS
    
    
#enable the script to be run from command line
if __name__ == "__main__":
    main()