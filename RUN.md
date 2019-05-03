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
- Proper comments have been included in each of the example file. Refer to the given *BMSS_User_Manual.pdf* for more detailed guidelines.
- The running model and the Sum Squared Error (SSE) followed by the type of optimizer (Global or Local) will be listed and continuously updated at the IPython console while the program is running. 
- After the run has completed, users will be prompted to insert either yes or no to export the model simulation data file (.csv).  
- Next, users will be asked to insert yes or no to set and visualize the corresponding SBOL-compliant gene circuit diagram. 
- The preview of the resulting sample circuits for each of the example files can be seen at the [GitHub Examples] folder.  
- All the output files will be exported to the __BMSSlib-master/Examples/Results__ folder (There will be a total of three files named according to the Year-Month-Day-hour-minute the time the run is completed: .txt lists  

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

[Spyder]: <https://docs.spyder-ide.org/overview.html>
[Anaconda Documentation]: <https://docs.anaconda.com/anaconda/user-guide/getting-started/>
[GitHub Examples]: <https://github.com/EngBioNUS/BMSSlib/tree/master/Examples>
