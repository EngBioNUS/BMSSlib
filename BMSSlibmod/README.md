## This folder lists all the python source files necessary to run the BMSS system

The following source files describe all the pre-established kinetic models for the three different regulatory systems: 
- ***InduciblePromoterLibrary_.py***
- ***ConstitutivePromoterLibrary_.py***
- ***LogicGatesLibrary_.py*** 

__*Note__: New models can be added in these files.



The three files below include the main operation functions of the BMSS:
- ***InducibleSystem_.py***
- ***ConstitutiveSystem_.py***
- ***LogicGatesSystem_.py***

__*Note__: Users can select the models to be tested under the *function* **RunModels** by commenting out those unwanted models.

| System Type | Number of Models Available |
| --- | --- |
| Inducible | 13 |
| Constitutive (Single Data Set) | 4 |
| Constitutive (Multiple Data Sets: Fixed RBS) | 4 |
| Constitutive (Multiple Data Sets: Fixed Promoter) | 2 |
| Logic Gate (NOT Gate) | 4 | 
| Logic Gate (AND Gate) | 6 |
| Logic Gate (OR Gate) | 11 |
