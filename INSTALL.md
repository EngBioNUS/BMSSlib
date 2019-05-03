# Installation: 
- Download the BMSS python [BMSSlib] package zip file from GitHub, and unzip to get the BMSSlib-master folder. 

There are two subfolders inside the BMSSlib-master folder: 
- __*BMSSlibmod*__: contains all the main source files
- __*Examples*__: contains the example main files for the three gene regulatory systems (Inducible System, Constitutive System, and Logic Gate System), and the example user input files (.csv) stored in the __*InputData*__ folder containing the characterization data. The generated output files are exported into the __*Results*__ folder.    

BMSSlib is implemented in Python 3.6.5 environment and can be *soon* installed using pip in Anaconda Prompt/terminal as shown below:
```
pip install BMSSlib
```

The system can be run on different OS systems: Windows, Mac, and Linux

New users can first download and install the [Anaconda Distribution] and choose the Python 3.7 version based on their Operating Systems (Windows, MacOS, Linux). The distribution comes with all the fundamental packages and GUI applications (Spyder, Jupyter Notebook etc.) of Python for scientific computing. The open source packages can be individually installed from the Anaconda repository with the *conda install* command or using the *pip install* command that is installed with Anaconda.   

In view of the potential compatibility issues for different versions of dependencies, it is highly recommended to create a virtual environment of __*Python 3.6.5*__ in Anaconda. 
- Open Anaconda Prompt (functions like command prompt) or Anaconda Navigator (Environments --> base (root) --> Open Terminal) to create a virtual environment as BMSSenv (or any other name) by typing the command below:
```
conda create -n BMSSenv python=3.6.5 anaconda
```
- When prompted with a message to ask whether to proceed with the packages installation or update, just insert y, refers to yes to proceed with the installation.
-	After all the installations have finished, insert the command below to activate the newly created environment
```
conda activate BMSSenv
```

In additional to those fundamental packages available in Anaconda, users are required to install two additional packages using the terminal at Anaconda Navigator or Anaconda Prompt: 
```
pip install tabulate
pip install tesbml
```



[BMSSlib]: <https://github.com/EngBioNUS/BMSSlib>
[Anaconda Distribution]: <https://www.anaconda.com/distribution/>
