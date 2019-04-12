# -*- coding: utf-8 -*-
"""
Created on Wed Dec 26 10:55:45 2018

@author: JingWui
"""

#############################################################################
# The Constitutive System class
#############################################################################

### Scripts imported ###
import Read_Data as ReadData
import ComputeAIC_ as ComputeAIC
import Txtfilename as Txtfilename
### Import Class from module
from ConstitutivePromoterLibrary_ import ConstitutivePromoterLibrary

### Python Packages Imported ###
import matplotlib.pyplot as plt
import time
import csv
import scipy.stats as ss
from tabulate import tabulate

class ConstitutiveSystem:
    
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
        self.ConstPLib1 = ConstitutivePromoterLibrary()
        
    def DataReader(self, Input_filename, NumDataSet): 
        
        ### Read Data ###
        print('Read Data')
        self.Data_header1, self.Data_array1, self.Data_stddev1 = ReadData.readData(Input_filename, NumDataSet)
        self.Input_filename_ = Input_filename
        self.NumDataSet_ = NumDataSet
        
        ReadData.plot_inputdata(self.Data_header1, self.Data_array1, self.Data_stddev1)
        
        ### Find Sample Size ###
        Data_array_numrows = self.Data_array1.shape[0]    # Time + All Inducers
        Data_array_numcols = self.Data_array1.shape[1]    # Number of RFP Data Per Inducer
        self.Sample_size = (Data_array_numrows - 1) * Data_array_numcols
        
    def RunModels(self, iteration = 1):
        
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
            
        #    ### Model 1 - Double-ODE Constitutive Promoter System ###
        #    SystemType = 'ConstDouble'
        #    Param_ConstDouble, SSE_ConstDouble, y0_ConstDouble, VarName_ConstDouble, ParamName_ConstDouble, ParamUnits_ConstDouble \
        #        = self.ConstPLib1.Run_ConstitutiveSystem(SystemType, self.Data_header1, self.Data_array1, self.NumDataSet_)
        #
        #    self.SSE_Combined.append(SSE_ConstDouble)
        #    self.ParamName_List.append(ParamName_ConstDouble)
        #    self.Param_List.append(Param_ConstDouble)
        #    self.ParamUnits_List.append(ParamUnits_ConstDouble)
        #    self.VarName_List.append(VarName_ConstDouble)
        #    self.y0_List.append(y0_ConstDouble)
        #    self.Num_Param_List.append(len(Param_ConstDouble))
        #    self.Model_List.append('Model 1 - ConstDouble')
        #    
        #    ### Model 1.1 - Double-ODE Constitutive Promoter System with Protein Maturation kinetics ###
        #    SystemType = 'ConstDoubleKMat'
        #    Param_ConstDoubleKMat, SSE_ConstDoubleKMat, y0_ConstDoubleKMat, VarName_ConstDoubleKMat, ParamName_ConstDoubleKMat, ParamUnits_ConstDoubleKMat \
        #        = self.ConstPLib1.Run_ConstitutiveSystem(SystemType, self.Data_header1, self.Data_array1, self.NumDataSet_)
        #
        #    self.SSE_Combined.append(SSE_ConstDoubleKMat)
        #    self.ParamName_List.append(ParamName_ConstDoubleKMat)
        #    self.Param_List.append(Param_ConstDoubleKMat)
        #    self.ParamUnits_List.append(ParamUnits_ConstDoubleKMat)
        #    self.VarName_List.append(VarName_ConstDoubleKMat)
        #    self.y0_List.append(y0_ConstDoubleKMat)
        #    self.Num_Param_List.append(len(Param_ConstDoubleKMat))
        #    self.Model_List.append('Model 1.1 - ConstDoubleKMat')
        
        #    # -------------------------------------------------------------------------------- #
        
#            ### Model 2 - Single-ODE Constitutive Promoter System ###
#            SystemType = 'ConstSingle'
#            Param_ConstSingle, SSE_ConstSingle, y0_ConstSingle, VarName_ConstSingle, ParamName_ConstSingle, ParamUnits_ConstSingle \
#                = self.ConstPLib1.Run_ConstitutiveSystem(SystemType, self.Data_header1, self.Data_array1, self.NumDataSet_)
#        
#            self.SSE_Combined.append(SSE_ConstSingle)
#            self.ParamName_List.append(ParamName_ConstSingle)
#            self.Param_List.append(Param_ConstSingle)
#            self.ParamUnits_List.append(ParamUnits_ConstSingle)
#            self.VarName_List.append(VarName_ConstSingle)
#            self.y0_List.append(y0_ConstSingle)
#            self.Num_Param_List.append(len(Param_ConstSingle))
#            self.Model_List.append('Model 2 - ConstSingle')
            
        #    ### Model 2.1 - Single-ODE Constitutive Promoter System with Protein Maturation Kinetics ###
        #    SystemType = 'ConstSingleKMat'
        #    Param_ConstSingleKMat, SSE_ConstSingleKMat, y0_ConstSingleKMat, VarName_ConstSingleKMat, ParamName_ConstSingleKMat, ParamUnits_ConstSingleKMat \
        #        = self.ConstPLib1.Run_ConstitutiveSystem(SystemType, self.Data_header1, self.Data_array1, self.NumDataSet_)
        #
        #    self.SSE_Combined.append(SSE_ConstSingleKMat)
        #    self.ParamName_List.append(ParamName_ConstSingleKMat)
        #    self.Param_List.append(Param_ConstSingleKMat)
        #    self.ParamUnits_List.append(ParamUnits_ConstSingleKMat)
        #    self.VarName_List.append(VarName_ConstSingleKMat)
        #    self.y0_List.append(y0_ConstSingleKMat)
        #    self.Num_Param_List.append(len(Param_ConstSingleKMat))
        #    self.Model_List.append('Model 2.1 - ConstSingleKMat')
            
        #    # -------------------------------------------------------------------------------- #
        
            ### Model 3 - Multi Double-ODE Constitutive Promoter System (Fix RBS) ###
            SystemType = 'MultiDoubleFixRBS'
            Param_MultiDoubleFixRBS, SSE_MultiDoubleFixRBS, y0_MultiDoubleFixRBS, VarName_MultiDoubleFixRBS, ParamName_MultiDoubleFixRBS, ParamUnits_MultiDoubleFixRBS \
                = self.ConstPLib1.Run_ConstitutiveSystem(SystemType, self.Data_header1, self.Data_array1, self.NumDataSet_)
        
            self.SSE_Combined.append(SSE_MultiDoubleFixRBS)
            self.ParamName_List.append(ParamName_MultiDoubleFixRBS)
            self.Param_List.append(Param_MultiDoubleFixRBS)
            self.ParamUnits_List.append(ParamUnits_MultiDoubleFixRBS)
            self.VarName_List.append(VarName_MultiDoubleFixRBS)
            self.y0_List.append(y0_MultiDoubleFixRBS)
            self.Num_Param_List.append(len(Param_MultiDoubleFixRBS))
            self.Model_List.append('Model 3 - MultiDoubleFixRBS')
            
            ### Model 3.1 - Multi Double-ODE Constitutive Promoter System (Fix RBS) with Protein Maturation Kinetics ###
            SystemType = 'MultiDoubleFixRBSKMat'
            Param_MultiDoubleFixRBSKMat, SSE_MultiDoubleFixRBSKMat, y0_MultiDoubleFixRBSKMat, VarName_MultiDoubleFixRBSKMat, ParamName_MultiDoubleFixRBSKMat, ParamUnits_MultiDoubleFixRBSKMat \
                = self.ConstPLib1.Run_ConstitutiveSystem(SystemType, self.Data_header1, self.Data_array1, self.NumDataSet_)
        
            self.SSE_Combined.append(SSE_MultiDoubleFixRBSKMat)
            self.ParamName_List.append(ParamName_MultiDoubleFixRBSKMat)
            self.Param_List.append(Param_MultiDoubleFixRBSKMat)
            self.ParamUnits_List.append(ParamUnits_MultiDoubleFixRBSKMat)
            self.VarName_List.append(VarName_MultiDoubleFixRBSKMat)
            self.y0_List.append(y0_MultiDoubleFixRBSKMat)
            self.Num_Param_List.append(len(Param_MultiDoubleFixRBSKMat))
            self.Model_List.append('Model 3.1 - MultiDoubleFixRBSKMat')
            
            ### Model 3.2 - Multi Single-ODE Constitutive Promoter System (Fix RBS) ###
            SystemType = 'MultiSingleFixRBS'
            Param_MultiSingleFixRBS, SSE_MultiSingleFixRBS, y0_MultiSingleFixRBS, VarName_MultiSingleFixRBS, ParamName_MultiSingleFixRBS, ParamUnits_MultiSingleFixRBS \
                = self.ConstPLib1.Run_ConstitutiveSystem(SystemType, self.Data_header1, self.Data_array1, self.NumDataSet_)
        
            self.SSE_Combined.append(SSE_MultiSingleFixRBS)
            self.ParamName_List.append(ParamName_MultiSingleFixRBS)
            self.Param_List.append(Param_MultiSingleFixRBS)
            self.ParamUnits_List.append(ParamUnits_MultiSingleFixRBS)
            self.VarName_List.append(VarName_MultiSingleFixRBS)
            self.y0_List.append(y0_MultiSingleFixRBS)
            self.Num_Param_List.append(len(Param_MultiSingleFixRBS))
            self.Model_List.append('Model 3.2 - MultiSingleFixRBS')
            
            ### Model 3.3 - Multi Single-ODE Constitutive Promoter System (Fix RBS) with Protein Maturation Kinetics ###
            SystemType = 'MultiSingleFixRBSKMat'
            Param_MultiSingleFixRBSKMat, SSE_MultiSingleFixRBSKMat, y0_MultiSingleFixRBSKMat, VarName_MultiSingleFixRBSKMat, ParamName_MultiSingleFixRBSKMat, ParamUnits_MultiSingleFixRBSKMat \
                = self.ConstPLib1.Run_ConstitutiveSystem(SystemType, self.Data_header1, self.Data_array1, self.NumDataSet_)
        
            self.SSE_Combined.append(SSE_MultiSingleFixRBSKMat)
            self.ParamName_List.append(ParamName_MultiSingleFixRBSKMat)
            self.Param_List.append(Param_MultiSingleFixRBSKMat)
            self.ParamUnits_List.append(ParamUnits_MultiSingleFixRBSKMat)
            self.VarName_List.append(VarName_MultiSingleFixRBSKMat)
            self.y0_List.append(y0_MultiSingleFixRBSKMat)
            self.Num_Param_List.append(len(Param_MultiSingleFixRBSKMat))
            self.Model_List.append('Model 3.3 - MultiSingleFixRBSKMat')
            
        #    # -------------------------------------------------------------------------------- #
        
        #    ### Model 4 - Multi Double-ODE Constitutive Promoter System (Fix Promoter) ###
        #    SystemType = 'MultiDoubleFixPromoter'
        #    Param_MultiDoubleFixPromoter, SSE_MultiDoubleFixPromoter, y0_MultiDoubleFixPromoter, VarName_MultiDoubleFixPromoter, ParamName_MultiDoubleFixPromoter, ParamUnits_MultiDoubleFixPromoter \
        #        = self.ConstPLib1.Run_ConstitutiveSystem(SystemType, self.Data_header1, self.Data_array1, self.NumDataSet_)
        #
        #    self.SSE_Combined.append(SSE_MultiDoubleFixPromoter)
        #    self.ParamName_List.append(ParamName_MultiDoubleFixPromoter)
        #    self.Param_List.append(Param_MultiDoubleFixPromoter)
        #    self.ParamUnits_List.append(ParamUnits_MultiDoubleFixPromoter)
        #    self.VarName_List.append(VarName_MultiDoubleFixPromoter)
        #    self.y0_List.append(y0_MultiDoubleFixPromoter)
        #    self.Num_Param_List.append(len(Param_MultiDoubleFixPromoter))
        #    self.Model_List.append('Model 4 - MultiDoubleFixPromoter')
        #    
        #    ### Model 4.1 - Multi Double-ODE Constitutive Promoter System (Fix Promoter) ###
        #    SystemType = 'MultiDoubleFixPromoterKMat'
        #    Param_MultiDoubleFixPromoterKMat, SSE_MultiDoubleFixPromoterKMat, y0_MultiDoubleFixPromoterKMat, VarName_MultiDoubleFixPromoterKMat, ParamName_MultiDoubleFixPromoterKMat, ParamUnits_MultiDoubleFixPromoterKMat \
        #        = self.ConstPLib1.Run_ConstitutiveSystem(SystemType, self.Data_header1, self.Data_array1, self.NumDataSet_)
        #
        #    self.SSE_Combined.append(SSE_MultiDoubleFixPromoterKMat)
        #    self.ParamName_List.append(ParamName_MultiDoubleFixPromoterKMat)
        #    self.Param_List.append(Param_MultiDoubleFixPromoterKMat)
        #    self.ParamUnits_List.append(ParamUnits_MultiDoubleFixPromoterKMat)
        #    self.VarName_List.append(VarName_MultiDoubleFixPromoterKMat)
        #    self.y0_List.append(y0_MultiDoubleFixPromoterKMat)
        #    self.Num_Param_List.append(len(Param_MultiDoubleFixPromoterKMat))
        #    self.Model_List.append('Model 4.1 - MultiDoubleFixPromoterKMat')
        
    def RunModelSelection(self):
        
        ### AIC ###
        self.AIC_Results = ComputeAIC.run_AIC(self.SSE_Combined, self.Sample_size, self.Num_Param_List)
        
        # -------------------------------------------------------------------------------- #
        
        ### Find Min SSE and AIC ###
        min_SSE = min(self.SSE_Combined)
        min_SSE_index = self.SSE_Combined.index(min_SSE)
        
        self.min_AIC = min(self.AIC_Results)
        self.min_AIC_index = self.AIC_Results.index(self.min_AIC)
        
        # -------------------------------------------------------------------------------- #
        
        ### Model Recommendation and Print out ###
        
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
        
        
        self.Time, self.VariableMatrix, self.DataLegend, self.ODEstring = self.ConstPLib1.Run_ConstitutivePlot\
                (BestSystemType, self.Init, self.Data_header1, self.Data_array1, self.Data_stddev1, self.FittedParams, self.NumDataSet_)
        
        
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
        Header = ["Model", "SSE", "AIC", "Rank"]
        for i in range(0, len(self.Model_List)):
            TableData.append([str(self.Model_List[i]), str(self.SSE_Combined[i]), str(self.AIC_Results[i]), str(int(Rank[i]))])
        Table = tabulate(TableData, Header, tablefmt='orgtbl')
        print(Table)
        
        ### Create and write to Txt File ###
        Txtfilename1, DateTimenow = Txtfilename.gettxtfilename()
        print('\nText File Generated:', Txtfilename)
        
        f = open(Txtfilename1,"a+")
        f.write('Input File name: '+self.Input_filename_+'\n')
        f.write('\n')
        f.write('Models Tested: '+str(self.Model_List)+'\n')
        f.write('\n')
        f.write(str(Table))
        f.write('\n')
        f.write('Recommended Model: '+self.Best_Model+'\n')
        f.write('\n')
        f.write('Optimized Parameters:\n')
        
        for i in self.ParamName:
            f.write('\t'+ i +' = '+ str(self.FittedParams[self.ParamName.index(i)])+'\n')
        f.write('\tdeg_mRNA = 0.1386\n')
        f.write('ODE:\n')
        for j in self.VarName:
            f.write('d'+ j + 'dt'+' = '+ self.ODEstring[self.VarName.index(j)]+'\n')
        f.write('\n')
        
        ### create a table to show the Relative Strength for Promoters or RBSs
        if 'Multi' and 'FixRBS' in self.Best_Model:
            TableRStrength = []
            Header = ["Promoter", "Relative Strength"]
            RStrength = []
            if 'Single' in self.Best_Model:        
                for s in self.ParamName:
                    if 'syn_Pep' in s:
                        RStrength.append (self.FittedParams[self.ParamName.index(s)])
            elif 'Double' in self.Best_Model:
                for s in self.ParamName:
                    if 'syn_mRNA' in s:
                        RStrength.append (self.FittedParams[self.ParamName.index(s)])
                
            else:
                print('Error on Tabulating Relative Strengths')
                
            MaxStrength = max(RStrength)
            for i in range(0, len(RStrength)):
                TableRStrength.append([str(self.Data_header1[i]), str(float("{0:.3f}".format(RStrength[i]/MaxStrength)))])
            TableRS = tabulate(TableRStrength, Header, tablefmt='orgtbl')
            print(TableRS)
                
            f.write(str(TableRS))
            f.write('\n\n')
        
        elif 'Multi' and 'FixPromoter' in self.Best_Model:
            TableRStrength = []
            Header = ["RBS", "Relative Strength"]
            RStrength = []   
            for s in self.ParamName:
                if 'syn_Pep' in s:
                    RStrength.append (self.FittedParams[self.ParamName.index(s)])
            MaxStrength = max(RStrength)
            for i in range(0, len(RStrength)):
                TableRStrength.append([str(self.Data_header1[i]), str(float("{0:.3f}".format(RStrength[i]/MaxStrength)))])
            TableRS = tabulate(TableRStrength, Header, tablefmt='orgtbl')
            print(TableRS)  
        
            f.write(str(TableRS))
            f.write('\n\n')     
        else:
            print('No Relative Strength Table available')
            
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
    
        while not ((ExportDataFile == 'yes') or (ExportDataFile == 'no')):
             ExportDataFile = input("Error: Incorrect Choice! Please insert either yes or no only:\n")
    
        if ExportDataFile == 'yes':
            CSVfileName = Txtfilename.getcsvfilename()
            VariableMatrixData = self.VariableMatrix[-1][:,:].tolist()
            with open(CSVfileName, 'w', newline='') as csvfile:
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
        
        ParName = self.ParamName + ['deg_mRNA'];
        Params = self.FittedParams.tolist() + [0.1386];
        ParamsUnit = self.ParamUnits + ['per_min', 'molperL'];
        
        ODE1 = []
        
        for o in range(0, len(ODE)):
            ODE1.append(ODE[o].replace("**", "^"))
            
        print(ODE1)
        
        exportsbml(ODE1, Variable, VarInit, ParName, Params, ParamsUnit)
        
    def PlotSBOLGraphics(self):
        
        ### To plot and visualize the gene circuit in SBOL visual compliant diagram
        import PlotCircuit as pc
        
        PlotCircuitFile = input("Please insert 'yes/no' to set and visualize the gene circuit diagram:")
        
        while not ((PlotCircuitFile == 'yes') or (PlotCircuitFile =='no')):
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
            Origin.append (input('Please insert the Name of Origin: '))
            print('Origin: ', Origin)
        
        #    Gene = []
        #    PartNum = input('Please insert the number of parts: ')
        #    for pt in range(0, int(PartNum)):
        #        Gene.append(input('Please insert the Name of gene '+str(pt+1)+': '))
        #    print('Gene: ', Gene)
            
            # MultiData with Fixed RBS models
            if 'Multi' and 'FixRBS' in self.Best_Model:
                RBS = []
                RBS.append (input('Please insert the Name of RBS: '))
                print('RBS: ', RBS)
                
#                Promoter = []
#          
#                PromoterNum = input('Please insert the number of promoters: ')
#                
#                while not PromoterNum.isdigit():
#                    PromoterNum = input('Error: Incorrect input! Please insert the right number: \n')
#                
#                for pr in range(0, int(PromoterNum)):
#                    Promoter.append (input('Please insert the Name of Promoter '+str(pr+1)+': '))
#                print('Promoters: ', Promoter)
                
                Promoter = self.Data_header1
                PromoterNum = len(Promoter)
                
                for pr in range(0, PromoterNum):
                    Input = "p.black." + Promoter[pr] + ".red" + " r.black." + RBS[0] + " c." + ReporterColor + Reporter+ " " + "t.black o.black." + Origin[0]
            
                    Regulations = None
                    
                    pc.Run_PlotCircuit(Input, Regulations)
                    
            # MultiData with Fixed Promoter models
            elif 'Multi' and 'FixPromoter' in self.Best_Model: 
                Promoter = []
                Promoter.append (input('Please insert the Name of Promoter: '))
                print('Promoter: ', Promoter)
#                RBS = []
#                RBSNum = input('Please insert the number of RBSs: ')
#                
#                while not RBSNum.isdigit():
#                    RBSNum = input('Error: Incorrect input! Please insert the right number: \n')
#                
#                for rb in range(0, int(RBSNum)):
#                    RBS.append (input('Please insert the Name of RBS '+str(rb+1)+': '))
#                print('RBSs: ', RBS)
                
                RBS = self.Data_header1
                RBSNum = len(RBS)
                
                for rb in range(0, RBSNum):
                    Input = "p.black." + Promoter[0] + " r.black." + RBS[rb] + ".red" + " c." + ReporterColor + Reporter+ " " + "t.black o.black." + Origin[0]
            
                    Regulations = None
                    
                    pc.Run_PlotCircuit(Input, Regulations)
                    
            else: 
                Promoter = []
                Promoter.append (input('Please insert the Name of Promoter: '))
                print('Promoter: ', Promoter)
                RBS = []
                RBS.append (input('Please insert the Name of RBS: '))
                print('RBS: ', RBS)
                
                Input = "p.black." + Promoter[0] + " r.black." + RBS[0] + " c." + ReporterColor + Reporter+ " " + "t.black o.black." + Origin[0]
            
                Regulations = None
                    
                pc.Run_PlotCircuit(Input, Regulations)
    
    ####################################################################            
    # Helper function to simplify the whole process   
    ####################################################################
         
    def AutoRunConstitutiveSystem(self, Input_filename, NumDataSet):
        
        self.DataReader(Input_filename, NumDataSet)
        self.RunModels()
        self.RunModelSelection()
        self.CreateOutputTextFile()
        self.ExportModelDataFile()
        self.ExportSBMLFile()
        self.PlotSBOLGraphics()
        
    def __del__(self):
        class_name = self.__class__.__name__
        print(class_name, "destroyed")
        
   