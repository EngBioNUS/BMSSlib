# -*- coding: utf-8 -*-
"""
Created on Wed Dec 26 14:16:04 2018

@author: JingWui
"""

#import sys,os
#print(sys.path.append(os.getcwd()))

import sys, os
sys.path.append("../BMSSlibmod") # adds same-level directory to python modules path.

## or add the corresponding path to the sys path as follows: 
#sys.path.append(os.getcwd().replace('Examples', 'BMSSlibmod')) 

#print(sys.path)
#print(os.path)

### Import Class from module
from BMSSlibmod.ConstitutiveSystem_ import ConstitutiveSystem
#from ConstitutiveSystem_ import ConstitutiveSystem

def main():
    
    ### create an object from the class ConstitutiveSystem
    CS = ConstitutiveSystem()
    
    ### User Input (Data FileName and Number of Data Set in the file)
#    Input_filename = 'Constitutive_p66.csv'
#    NumDataSet = 1
    
    ### multi data set (varied Promoters with fixed RBS)
    Input_filename = 'Constitutive_4xmultiFixRBS.csv'
    NumDataSet = 4
    
#    CS.DataReader(Input_filename, NumDataSet)
#    CS.RunModels()
#    CS.RunModelSelection()
#    CS.CreateOutputTextFile()
#    CS.ExportModelDataFile()
#    CS.ExportSBMLFile()
#    CS.PlotSBOLGraphics()
    
    #############################################################
    # Or use Helper function to autorun all functions
    #############################################################
    CS.AutoRunConstitutiveSystem(Input_filename, NumDataSet)
    
    del CS
    
#enable the script to be run from command line
if __name__ == "__main__":
    main()