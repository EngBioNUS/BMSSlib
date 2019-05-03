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
Spyder, the Scientific Python Development Environment, is a free integrated development environment (IDE) that is included with Anaconda. It includes editing, interactive testing, debugging and introspection features. Launch Spyder from your start menu or Anaconda Navigator. Refer to the [Anaconda Documentation] for more details.  

To get started, please open the example files (Inducible, Constitutive, or logic gate system) from __*Examples*__ folder in Spyder. Users are required to specify the input data filename (.csv) and provide few information as requested before running the example file. The input data file should follow the format as shown in the example files for the three different regulatory systems. The file should be accessible in the same folder as the example file that you need to run, or else please specify the absolute path to the data file. For advance users who wish to append their new models into the model bank, they can refer to the model function (in the form of ordinary differential equations) presented in the particular regulatory system library and modify accordingly. 

## System Overview
![GitHub_SystemOverview](https://user-images.githubusercontent.com/32381993/57119532-8dee2b80-6d9d-11e9-9126-3208d641b662.png)
The user input characterization data in csv file (including time, fluorescence level relative to cell growth, and their corresponding standard deviations) was read by the system data reader for units conversion and dose-response prefitting (only applicable for inducible system). Numerical integration, optimization, and plotting packages were used to solve the ordinary differential equations, iteratively fit models retrieved from the model bank to experimental data through minimizing the sum squared residuals and plot the graphical results for data visualizations. The model selection algorithm was based on the AIC statistical inference criterion. The output figures displayed the time-response and steady-state response behaviors in conjunction with the imported measured experimental data and their corresponding standard deviations represented in error bars. User-defined customizable SBOL-compliant gene circuit graphics rendering was also exhibited by the BMSS platform. The details of the model ranking and the best model candidate were listed in the output text file with the model simulation data exported in a separate csv file. The SBML representation of the chosen model was also programmed and exported in XML format.

All the pre-established model formulations available in the model bank are described in the __*Models_Documentation.pdf*__ file. 

## BMSS-validated Results
### Inducible System
![GitHub_InducibleSystem1](https://user-images.githubusercontent.com/32381993/57119877-081faf80-6da0-11e9-9398-e7fcd25db339.png)

### Constitutive System
![github_constitutivesystem](https://user-images.githubusercontent.com/32381993/50503415-d9b45f80-0aa1-11e9-9692-99b0fffde1a4.png)

### Logic Gate System
![github_logicgatessystem](https://user-images.githubusercontent.com/32381993/50503680-a8d52a00-0aa3-11e9-9bab-680b900a6f4e.png)


## References
* [SimpleSBML]  - supports easy SBML model construction and editing
* [DNAplotlib]  - generates SBOL Visual compliant diagrams
* [constrNMPy]  - package for constrained Nelder-Mead optimization

## License

Copyright __*2019 EngBioNUS*__

Licensed under the __Apache License, Version 2.0__ (the "License"); you may not use this file except in compliance with the License.
You may obtain a copy of the License at
```
    http://www.apache.org/licenses/LICENSE-2.0
```
Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the License.

[//]: # (These are reference links used in the body of this note and get stripped out when the markdown processor does its job. There is no need to format nicely because it shouldn't be seen. Thanks SO - http://stackoverflow.com/questions/4823468/store-comments-in-markdown-syntax)


   [SimpleSBML]: <https://github.com/sys-bio/simplesbml>
   [DNAplotlib]: <https://github.com/VoigtLab/dnaplotlib>
   [constrNMPy]: <https://github.com/alexblaessle/constrNMPy>
   [Anaconda Distribution]: <https://www.anaconda.com/distribution/>
   [Anaconda Documentation]: <https://docs.anaconda.com/anaconda/user-guide/getting-started/>

