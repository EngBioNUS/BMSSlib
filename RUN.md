# Quick Start 
[Spyder], the Scientific Python Development Environment, is a free integrated development environment (IDE) that is included in Anaconda. The IDE includes editing, interactive testing, debugging and introspection features. Refer to the [Anaconda Documentation] for more details. 

- Launch __Spyder (BMSSenv)__ from your start menu or Anaconda Navigator. 
- In Spyder, check if the right __Python version (3.6.5)__ is used in the IPython console shown at the right bottom.  
- IPython console allows users to run code by line, cell, or file, and render plots right inline.
- Editor is where users can input all their codes for running. 
- Click on the green Run button to execute files. 
- To display the output Graphics in a separate window, select Tools --> Preferences --> IPython console --> Graphics, change the option of Backend from Inline to Qt5. 
-	Click Apply then OK
- Close the Spyder IDE and reopen the IDE from Start menu to apply the changed setting.

## To Run Example File
- To get started, at Spyder, File --> Open to browse and open the example file (__*Example_InducibleSystem, Example_ConstitutiveSystem, or Example_LogicGatesSystem*__) in .py located in the __BMSSlib-master/Examples/InputData__ folder. Select the file based on your system of interest. 
- Click on the green Run button to execute the example file without modification.  
- Proper comments have been included in each of the example file. Refer to the given __*BMSS_User_Manual.pdf*__ for more detailed guidelines.
- The running model and the Sum Squared Error (SSE) followed by the type of optimizer (Global or Local) will be listed and continuously updated at the IPython console while the program is running. 
- After the run has completed, users will be prompted to insert either yes or no to export the model simulation data file (.csv).  
- Next, users will be asked to insert yes or no to set and visualize the corresponding SBOL-compliant gene circuit diagram. 
- The preview of the resulting sample circuits for each of the example files can be seen at the [GitHub Examples] folder.  
- All the output files will be exported to the __BMSSlib-master/Examples/Results__ folder (There will be a total of three files named according to the Year-Month-Day-hour-minute the time the run is completed: *.txt* details the model ranking table and the recommended best model candidate with its model formulation and estimated parameters; *.csv* includes the model simulation results for post-processing; and *.xml* file encodes the details of the best model candidate in Systems Biology Markup Language (SBML) format. The SBML file can be imported to other CAD tools for post-processing or further simulations.   

### The Output Specifications for the Example Files (to be inserted when prompted at console):

<__Note__: Insert only the Italic bold information>

__Inducible System__: 
- Reporter type: __*RFP*__
- Number of Plasmids: __*1*__
- Name of Origin: __*CoIE1*__
- Number of parts: __*2*__
- Name of gene 1: __*AraC*__
- Name of gene 2: __*RFP*__
- Name of Inducible Promoter: __*pBAD*__ 

__Constitutive System__: 
- Reporter type: __*RFP*__
- Name of Origin: __*CoIE1*__
- Name of RBS: __*rbs34*__

__Logic Gate System__:
- Reporter type: __*GFP*__
- Number of Plasmids: __*2*__
- Name of Origin 1: __*CoIE1*__
- Name of Origin 2: __*p15A*__
- Number of parts: __*3*__
- Name of gene 1: __*gRNA*__
- Name of gene 2: __*dCas9*__
- Name of gene 3: __*GFP*__

<__Note__: With proper installation process, you should be able to run the example file without errors. If there is any error shown at the console after running the file, please restart your Spyder IDE and try again. If the error persists, please feel free to email the bugs to us at *EngBioBMSS.help@gmail.com*>

## To Run New User Input File 
- Prepare the processed characterization data file in .csv following the example characterization data files for a particular system (Inducible, Constitutive, Logic Gate System) provided in the __BMSSlib-master/Examples/InputData__ folder. More descriptions are also documented in the example python file. 
- Saved the Input Data file in the __BMSSlib-master/Examples/InputData__ folder. __Note__: Follow the naming rules of starting with letters, then followed by letters, digits, or underscores. No spacing is allowed.   
- At the corresponding Example file in Spyder, users can modify the Input_filename, NumDataSet/NumState, SystemOpt, or/and Inducer_unit and OptInhibition (only applicable for Inducible System) accordingly. 
- Click on the green Run button to execute the Example python file. 
- It is advisable to repeat at least three times (depends on the nature of the characterization data trends, some models are just harder to converge to best fit the experimental data) to verify if any models from the BMSS libraries could well capture the newly imported .csv characterization file. 
- The displayed output graphic figures can be saved by clicking on the Save icon on top of the individual graphic windows.
- Again, check the exported output *.txt* file, *.csv* file, and *.xml* (SBML) file in the __BMSSlib-master/Examples/Results__ folder.  

[Spyder]: <https://docs.spyder-ide.org/overview.html>
[Anaconda Documentation]: <https://docs.anaconda.com/anaconda/user-guide/getting-started/>
[GitHub Examples]: <https://github.com/EngBioNUS/BMSSlib/tree/master/Examples>
