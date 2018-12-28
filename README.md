# BMSSlib

Thank you for downloading BMSSlib! This library supports an automated Bio-Model Selection System which enables users to derive the best mathematical model candidate to capture the transient dynamic gene circuit behaviors in a fast, efficient, and automated way using part/circuit characterization data.

This package supports __*three*__ routinely used __gene regulatory systems__:
 
- *Inducible Systems*
- *Constitutive Systems (Single Dataset or Multiple Datasets)*
- *Logic Gates Systems (NOT, AND, and OR gates)*

## Features

- Read in User Input characterization data in .csv file
- Model Fitting and Model Selection from a library of pre-established models.
- Exhibit the graphical figures of processed experimental data, experimental & Simulation data, and/or steady state results.    
- Export the best model candidate as SBML file (.xml)
- Generate SBOL visual-compliant gene circuit diagrams for the best model candidate 
- Export Model Data in .csv file

## Installation

BMSSlib is implemented in *Python 3.6.5* environment and can be installed using *pip* as shown below:  

```
pip install BMSSlib
```


__*Note*__: In view of the potential compatibility issues for the different versions of dependencies, it is highly recommended to create virtual environment of Python 3.6.5 in Anaconda. 
```
conda create -n BMSSenv python=3.6.5 anaconda
conda activate BMSSenv
```
Two required packages can then be installed using *pip* as follows:
```
pip install tabulate
pip install tesbml
```

## Getting Started
Please refer to the example files for each of the three gene regulatory systems to get you started. For advance users who wish to append their new models into the model bank, they can refer to the model function presented in the particular system library and modify accordingly. 

## System Overview
![manuscript_figure1](https://user-images.githubusercontent.com/32381993/50499775-58e86a00-0a87-11e9-9993-5ed192d7aec2.png)

## BMSS-validated Results
### Inducible System
![github_induciblesystem](https://user-images.githubusercontent.com/32381993/50501327-593a3280-0a92-11e9-9491-4574241672e3.png)

### Constitutive System
![github_constitutivesystem](https://user-images.githubusercontent.com/32381993/50503415-d9b45f80-0aa1-11e9-9692-99b0fffde1a4.png)

### Logic Gates System
![github_logicgatessystem](https://user-images.githubusercontent.com/32381993/50503680-a8d52a00-0aa3-11e9-9bab-680b900a6f4e.png)


## References
* [SimpleSBML]  - supports easy SBML model construction and editing
* [DNAplotlib]  - generates SBOL Visual compliant diagrams
* [constrNMPy]  - package for constrained Nelder-Mead optimization

## License

MIT

[//]: # (These are reference links used in the body of this note and get stripped out when the markdown processor does its job. There is no need to format nicely because it shouldn't be seen. Thanks SO - http://stackoverflow.com/questions/4823468/store-comments-in-markdown-syntax)


   [SimpleSBML]: <https://github.com/sys-bio/simplesbml>
   [DNAplotlib]: <https://github.com/VoigtLab/dnaplotlib>
   [constrNMPy]: <https://github.com/alexblaessle/constrNMPy>

