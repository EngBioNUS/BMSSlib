```diff
- # <Note: The BMSS system is still under development, and the package will be pip-installable after a while>
```

# BMSSlib

Thank you for downloading BMSSlib! This library supports an automated Bio-Model Selection System which enables users to derive the best mathematical model candidate to capture the transient dynamic gene circuit behaviors in a fast, efficient, and automated way using part/circuit characterization data.

This package supports __*three*__ routinely used __gene regulatory systems__:
 
- *Inducible System*
- *Constitutive System (Single Dataset or Multiple Datasets)*
- *Logic Gate System (NOT, AND, and OR gates)*

## Features

- Read in User Input characterization data in .csv file
- Model Fitting and Model Selection from a library of pre-established models.
- Exhibit the graphical figures of processed experimental data, experimental & Simulation data, and/or steady state results.    
- Export the best model candidate as SBML file (.xml)
- Generate SBOL visual-compliant gene circuit diagrams for the best model candidate 
- Export Model Data in .csv file

## Installation

BMSSlib is implemented in *Python 3.6.5* environment and can be *soon* installed using *pip* as shown below:  

```
pip install BMSSlib
```
New users can first download and install the [Anaconda Distribution] and choose the Python 3.7 version based on their Operating Systems (Windows, MacOS, Linux).  

__*Note*__: In view of the potential compatibility issues for different versions of dependencies, it is highly recommended to create virtual environment of Python 3.6.5 in Anaconda. Open Anaconda Prompt (like command prompt) or Anaconda Navigator (Environments --> base (root) --> Open Terminal) to create a virtual environment as BMSSenv and then activate the environment following the commands below.    
```
conda create -n BMSSenv python=3.6.5 anaconda
conda activate BMSSenv
```
Two required packages, in addition to those fundamental packages available in Anaconda, can then be installed using *pip* as follows:
```
pip install tabulate
pip install tesbml
```

## Getting Started
Please refer to the example files (Inducible, Constitutive, or logic gate system) in __*Examples*__ folder to get you started. Users are required to specify the input data filename (.csv) and provide few information as requested before running the example file. For advance users who wish to append their new models into the model bank, they can refer to the model function (in the form of ordinary differential equations) presented in the particular regulatory system library and modify accordingly. 

## System Overview
![manuscript_figure1](https://user-images.githubusercontent.com/32381993/50499775-58e86a00-0a87-11e9-9993-5ed192d7aec2.png)
The user input characterization data in csv file (including time, fluorescence level relative to cell growth, and their corresponding standard deviations) was read by the system data reader for units conversion and dose-response prefitting (only applicable for inducible system). Numerical integration, optimization, and plotting packages were used to solve the ordinary differential equations, iteratively fit models retrieved from the model bank to experimental data through minimizing the sum squared residuals and plot the graphical results for data visualizations. The model selection algorithm was based on the AIC statistical inference criterion. The output figures displayed the time-response and steady-state response behaviors in conjunction with the imported measured experimental data and their corresponding standard deviations represented in error bars. User-defined customizable SBOL-compliant gene circuit graphics rendering was also exhibited by the BMSS platform. The details of the model ranking and the best model candidate were listed in the output text file with the model simulation data exported in a separate csv file. The SBML representation of the chosen model was also programmed and exported in XML format.

All the pre-established model formulations available in the model bank are described in the __*Models_Documentation.pdf*__ file. 

## BMSS-validated Results
### Inducible System
![github_induciblesystem](https://user-images.githubusercontent.com/32381993/50501327-593a3280-0a92-11e9-9491-4574241672e3.png)

### Constitutive System
![github_constitutivesystem](https://user-images.githubusercontent.com/32381993/50503415-d9b45f80-0aa1-11e9-9692-99b0fffde1a4.png)

### Logic Gate System
![github_logicgatessystem](https://user-images.githubusercontent.com/32381993/50503680-a8d52a00-0aa3-11e9-9bab-680b900a6f4e.png)


## References
* [SimpleSBML]  - supports easy SBML model construction and editing
* [DNAplotlib]  - generates SBOL Visual compliant diagrams
* [constrNMPy]  - package for constrained Nelder-Mead optimization

## License

No license. Permission from authors is required. 

[//]: # (These are reference links used in the body of this note and get stripped out when the markdown processor does its job. There is no need to format nicely because it shouldn't be seen. Thanks SO - http://stackoverflow.com/questions/4823468/store-comments-in-markdown-syntax)


   [SimpleSBML]: <https://github.com/sys-bio/simplesbml>
   [DNAplotlib]: <https://github.com/VoigtLab/dnaplotlib>
   [constrNMPy]: <https://github.com/alexblaessle/constrNMPy>
   [Anaconda Distribution]: <https://www.anaconda.com/distribution/>

