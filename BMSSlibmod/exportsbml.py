# -*- coding: utf-8 -*-
"""
Created on Mon Aug  6 14:31:01 2018

@author: jingwui
"""

### This module is to export the model in SBML format (.xml) ###
#input argument:
#    ODE: list of ODE equations in string format
#    variable: list of variables in string format
#    Init: list of Initial values for variables
#    paramName: list of parameters names in string format
#    param: list of parameters values
#    paramUnit: list of parameters Units
#output:
#    SBML file in .xml


import simplesbml #use simplesbml module
import Txtfilename

def exportsbml(ODE, variable, Init, paramName, param, paramUnit):

    model = simplesbml.sbmlModel();

    Variable = [0]*len(variable);

    for i in range(0, len(variable)):
        Variable[i] = '['+variable[i]+']';

    for s in range(0, len(variable)):
        model.addSpecies(Variable[s], Init[s]);

    for p in range(0, len(param)):
        model.addParameter(paramName[p], param[p], paramUnit[p]);

    for r in range(0, len(ODE)):
        model.addRateRule(variable[r], ODE[r]);


    Model = model.toSBML()

    XMLfileName = Txtfilename.getxmlfilename()

    print(Model, file=open(XMLfileName, "w"))


###---------------------------------------------------------------------------#
# the main file to test running the exportsbml function (Test function)

def Main_exportsbml():

    dx = 'sigma *(y - x)'
    dy = 'x*(rho -z) - y'
    dz = 'x*y - beta *z'

    ODE = [dx, dy, dz];
    variable = ['x', 'y','z'];
    Init = [0.8, 1.8, 19];
    paramName = ['sigma', 'rho', 'beta'];
    param = [20.0, 28, 3];
    paramUnit = ['Dimension_less', 'Dimension_less', 'Dimension_less'];

    exportsbml(ODE, variable, Init, paramName, param, paramUnit)
    
# Enable the script to be run from the command line
if __name__ == "__main__":
	Main_exportsbml()