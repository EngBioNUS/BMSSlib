# -*- coding: utf-8 -*-
"""
Created on Wed Dec 26 16:27:29 2018

@author: JingWui
"""

import sys, os
#sys.path.append("..") # adds higher directory to python modules path.
sys.path.append("../BMSSlibmod")

from BMSSlibmod.InducibleSystem_ import InducibleSystem

def main():
    
    # create an object from the class InducibleSystem
    IS = InducibleSystem()    
    
    # User Input (Data FileName and Number of Data Set in the file)
    Input_filename = 'pBAD_Ara_rbspBb_MG_csv.csv' # Arabinose (Data set 1)
    NumDataSet = 9
    Inducer_unit = '%' # Enter only (10**-6, 10**-9 or '%')
    OptInhibition = None
    
#    IS.DataReader(Input_filename, NumDataSet, Inducer_unit)
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
    
    IS.AutoRunInducibleSystem(Input_filename, NumDataSet, Inducer_unit, OptInhibition)
        
    del IS
    
#enable the script to be run from command line
if __name__ == "__main__":
    main()