# -*- coding: utf-8 -*-
"""
Created on Wed Dec 26 18:07:21 2018

@author: JingWui
"""

import sys, os
sys.path.append("../BMSSlibmod") # adds same-level directory to python modules path.

from BMSSlibmod.LogicGatesSystem_ import LogicGatesSystem

def main():
    
    # create an object from the class LogicGatesSystem
    LGS = LogicGatesSystem()  
    
    # User Input (Data FileName and Number of Data Set in the file) 
    Input_filename = 'LogicGate_NOTTop10d30LB_1.csv' #NOT gate (Data set 1)
    NumState = 2
    
#    LGS.DataReader(Input_filename, NumState)
#    LGS.RunModels()
#    LGS.RunModelSelection()
#    LGS.CreateOutputTextFile()
#    LGS.ExportModelDataFile()
#    LGS.ExportSBMLFile()
#    LGS.PlotSBOLGraphics()
    
    #############################################################
    # Or use Helper function to autorun all functions
    #############################################################
    
    LGS.AutoRunLogicGatesSystem(Input_filename, NumState)
        
    del LGS
    
#enable the script to be run from command line
if __name__ == "__main__":
    main()