# -*- coding: utf-8 -*-
"""
Released on April 29, 2019

@author: Yeoh Jing Wui <bchyjw@nus.edu.sg>; Poh Chueh Loo <poh.chuehloo@nus.edu.sg>

The code is part of BMSS software.

Copyright (c) 2019, National University of Singapore.

"""


import sys, os

# adds higher level directory to python modules path to access the key modules
sys.path.append("../BMSSlibmod") 


### Import Class from module
from ConstitutiveSystem_ import ConstitutiveSystem
#from ConstitutiveSystem_ import ConstitutiveSystem

def main():
    
    ### create an object from the class ConstitutiveSystem
    CS = ConstitutiveSystem()
    
    #############################################################
    # User Input 
    
    # Input_filename: The input data file in .csv following the sequences below:
    # 
    # Single Constitutive Dataset: [Time(min) | Promoter1 | SDPromoter1 ]
    # Multiple Constitutive Dataset with fixed Promoters ( <= 6 datasets): 
    # [Time(min) | RBS1 | RBS2 | RBS3 | SDRBS1 | SDRBS2 | SDRBS3 ]
    # Multiple Constitutive Dataset with fixed RBSs ( <= 6 datasets): 
    # [Time(min) | Promoter1 | Promoter2 | Promoter3 | SDPromoter1 | SDPromoter2 | SDPromoter3 ]
    #
    # SD refers to the standard deviation data
    #
    # NumState: Number of Data Set (excluding standard deviations)
    # SystemOpt: Options for Constitutive system ('SingleConst', 'MultiFixedRBS', 'MultiFixedPromoter')  
    # 
    # Note: the SystemOpt is case-insensitive!
    #############################################################
    
    ### User Input (Data FileName and Number of Data Set in the file)
#    Input_filename = 'Constitutive_p66.csv'
#    NumDataSet = 1
#    SystemOpt = 'SingleConst'
    
    ### multi data set (varied Promoters with fixed RBS)
    Input_filename = 'Constitutive_4xmultiFixRBS.csv'
    NumDataSet = 4
    SystemOpt = 'MultiFixedRBS' 
    
    ### multi data set (varied RBSs with fixed Promoter)
#    SystemOpt = 'MultiFixedPromoter'
    
    #############################################################
    
    # To access the input file stored in the InputData folder
    Inputfile = "InputData\\" + Input_filename
    
    #############################################################
    # The sequences of key functions used 
    #############################################################
    
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
    
    CS.AutoRunConstitutiveSystem(Inputfile, SystemOpt, NumDataSet)
    
    sys.path.remove('../BMSSlibmod') #remove the path
    del CS
    
#enable the script to be run from command line
if __name__ == "__main__":
    main()
