# -*- coding: utf-8 -*-
"""
Created on Wed July 18 09:50:57 2018

@author: JingWui
"""

import numpy as np
import constrNMPy as cNM
import matplotlib.pyplot as plt

from scipy.integrate import odeint
from scipy.optimize import differential_evolution

#import sympy as sym   #use symbolic mathematics to replace eval()

# Note: Please ignore the warning given by code analysis in all the solveODE functions
# Those local variables will be evaluated in eval();
# Using symsy instead of eval() would critically increase the computational time.

class LogicGatesLibrary:

#    mRNA1, Pep1, mRNA2, Pep2, mRNA3, Pep3 = sym.symbols(('mRNA1', 'Pep1', \
#                                                            'mRNA2', 'Pep2',\
#                                                            'mRNA3', 'Pep3'))
#
#    state1, state2 = sym.symbols(('state1','state2'))
#
#    syn_mRNA1, syn_mRNA2, syn_mRNA3, deg_mRNA, syn_Pep, deg_Pep, Pep1max, Pep2max \
#        = sym.symbols(('syn_mRNA1', 'syn_mRNA2','syn_mRNA3','deg_mRNA',\
#                       'syn_Pep','deg_Pep', 'Pep1max', 'Pep2max'))

    ### ODE Model for NOT logic gate (one INPUT, STATES: 0, 1) ###
    def solveODE_NOTgate(y, t, state, param, Operation = 'Solve'):
        # Dependent variables
        mRNA1 = y[0] # Col 0 of ODESoln
        Pep1 = y[1]  # Col 1 of ODESoln
        mRNA2 = y[2]
        Pep2 = y[3]

        # Parameters
        syn_mRNA1 = param[0]
        syn_mRNA2 = param[1]
        deg_mRNA = 0.1386   # Fixed as ln2/5(mins) = 0.1386
        syn_Pep = param[2]
        deg_Pep = param[3]
        Kmaxrep = param[4]
        Pepmax = param[5]

        # Differential equations in string
        dmRNA1_dt = 'syn_mRNA1*(state) - (deg_mRNA *mRNA1)'
        dPep1_dt = '(syn_Pep*mRNA1) - (deg_Pep*Pep1)'
        dmRNA2_dt = '(syn_mRNA2*(1-Kmaxrep*(Pep1/Pepmax)))-(deg_mRNA * mRNA2)'
        dPep2_dt = '(syn_Pep*mRNA2)-(deg_Pep*Pep2)'

        if Operation == 'Solve':
            # Return differential equations solution
            #return [dmRNA1_dt, dPep1_dt, dmRNA2_dt, dPep2_dt]
            return [eval(dmRNA1_dt), eval(dPep1_dt), eval(dmRNA2_dt), eval(dPep2_dt)]
            #return [dmRNA1_dt_, dPep1_dt_, dmRNA2_dt_, dPep2_dt_, dmRNA3_dt_, dPep3_dt_]
        elif Operation == 'GetODE':
            return [dmRNA1_dt, dPep1_dt, dmRNA2_dt, dPep2_dt]
        else:
            print('Error: Please Enter the correct Operation for solveODE function')

    def solveODE_NOTgateKMat(y, t, state, param, Operation = 'Solve'):
        # Dependent variables
        mRNA1 = y[0] # Col 0 of ODESoln
        Pep1 = y[1]  # Col 1 of ODESoln
        mRNA2 = y[2]
        Pep2 = y[3]
        Pep2m = y[4]

        # Parameters
        syn_mRNA1 = param[0]
        syn_mRNA2 = param[1]
        deg_mRNA = 0.1386   # Fixed as ln2/5(mins) = 0.1386
        syn_Pep = param[2]
        deg_Pep = param[3]
        Kmaxrep = param[4]
        Pepmax = param[5]
        Kmature = param[6]
        #Kmature = 0.154 # log(2)/4.5 min

        # Differential equations in string
        dmRNA1_dt = 'syn_mRNA1*(state) - (deg_mRNA *mRNA1)'
        dPep1_dt = '(syn_Pep*mRNA1) - (deg_Pep*Pep1)'
        dmRNA2_dt = '(syn_mRNA2*(1-Kmaxrep*(Pep1/Pepmax)))-(deg_mRNA * mRNA2)'
        dPep2_dt = '(syn_Pep*mRNA2) - (Kmature*Pep2)'
        dPep2m_dt = '(Kmature*Pep2)-(deg_Pep*Pep2m)'

        if Operation == 'Solve':
            # Return differential equations solution
            #return [dmRNA1_dt, dPep1_dt, dmRNA2_dt, dPep2_dt, dPep2m_dt]
            return [eval(dmRNA1_dt), eval(dPep1_dt), eval(dmRNA2_dt), eval(dPep2_dt), eval(dPep2m_dt)]
        elif Operation == 'GetODE':
            return [dmRNA1_dt, dPep1_dt, dmRNA2_dt, dPep2_dt, dPep2m_dt]
        else:
            print('Error: Please Enter the correct Operation for solveODE function')

    ### ODE Model for NOT logic gate repressible promoter only one part (one INPUT, STATES: 0, 1) ###
    def solveODE_NOTgateSingle(y, t, state, param, Operation = 'Solve'):
        # Dependent variables
        mRNA1 = y[0] # Col 0 of ODESoln
        Pep1 = y[1]  # Col 1 of ODESoln

        # Parameters
        syn_mRNA = param[0]
        deg_mRNA = 0.1386   # Fixed as ln2/5(mins) = 0.1386
        syn_Pep = param[1]
        deg_Pep = param[2]
        Kmaxrep = param[3]

        # Differential equations in string
        dmRNA1_dt = '(syn_mRNA*(1-Kmaxrep*(state)))-(deg_mRNA * mRNA1)'
        dPep1_dt = '(syn_Pep*mRNA1)-(deg_Pep*Pep1)'

        if Operation == 'Solve':
            # Return differential equations solution
            #return [dmRNA1_dt, dPep1_dt]
            return [eval(dmRNA1_dt), eval(dPep1_dt)]
        elif Operation == 'GetODE':
            return [dmRNA1_dt, dPep1_dt]
        else:
            print('Error: Please Enter the correct Operation for solveODE function')

    ### ODE Model for NOT logic gate repressible promoter only one part (one INPUT, STATES: 0, 1) ###
    def solveODE_NOTgateSingleKMat(y, t, state, param, Operation = 'Solve'):
        # Dependent variables
        mRNA1 = y[0] # Col 0 of ODESoln
        Pep1 = y[1]  # Col 1 of ODESoln
        Pep1m = y[2]

        # Parameters
        syn_mRNA = param[0]
        deg_mRNA = 0.1386   # Fixed as ln2/5(mins) = 0.1386
        syn_Pep = param[1]
        deg_Pep = param[2]
        Kmaxrep = param[3]
        Kmature = param[4]

        # Differential equations in string
        dmRNA1_dt = '(syn_mRNA*(1-Kmaxrep*(state)))-(deg_mRNA * mRNA1)'
        dPep1_dt = '(syn_Pep*mRNA1)-(Kmature*Pep1)'
        dPep1m_dt = '(Kmature*Pep1)-(deg_Pep*Pep1m)'

        if Operation == 'Solve':
            # Return differential equations solution
            #return [dmRNA1_dt, dPep1_dt]
            return [eval(dmRNA1_dt), eval(dPep1_dt), eval(dPep1m_dt)]
        elif Operation == 'GetODE':
            return [dmRNA1_dt, dPep1_dt, dPep1m_dt]
        else:
            print('Error: Please Enter the correct Operation for solveODE function')

    ### ODE Model for AND gate (two INPUT, STATES: 00,01,10,11)  ###
    def solveODE_ANDgate(y, t, state, param, Operation = 'Solve'):

        # Dependent variables
        mRNA1 = y[0] # Col 1
        Pep1 = y[1]  # Col 2
        mRNA2 = y[2]
        Pep2 = y[3]
        mRNA3 = y[4]
        Pep3 = y[5]

        state1 = state[0]
        state2 = state[1]

        # Parameters
        syn_mRNA1 = param[0]
        syn_mRNA2 = param[1]
        syn_mRNA3 = param[2]
        deg_mRNA = 0.1386 # Fixed as ln2/5(mins) = 0.1386
        syn_Pep = param[3]
        deg_Pep = param[4]
        Pep1max = param[5]
        Pep2max = param[6]

#        mRNA1, Pep1, mRNA2, Pep2, mRNA3, Pep3 = sym.symbols(('mRNA1', 'Pep1', \
#                                                            'mRNA2', 'Pep2',\
#                                                            'mRNA3', 'Pep3'))
#
#        state1, state2 = sym.symbols(('state1','state2'))
#
#        syn_mRNA1, syn_mRNA2, syn_mRNA3, deg_mRNA, syn_Pep, deg_Pep, Pep1max, Pep2max \
#        = sym.symbols(('syn_mRNA1', 'syn_mRNA2','syn_mRNA3','deg_mRNA',\
#                       'syn_Pep','deg_Pep', 'Pep1max', 'Pep2max'))

        # ODEs in String
        dmRNA1_dt = 'syn_mRNA1*(state1) - (deg_mRNA *mRNA1)'
        dPep1_dt = '(syn_Pep*mRNA1) - (deg_Pep*Pep1)'
        dmRNA2_dt = 'syn_mRNA2*(state2) - (deg_mRNA *mRNA2)'
        dPep2_dt = '(syn_Pep*mRNA2) - (deg_Pep*Pep2)'
        dmRNA3_dt = '(syn_mRNA3*(Pep1/Pep1max)*(Pep2/Pep2max))-(deg_mRNA * mRNA3)'
        dPep3_dt = '(syn_Pep*mRNA3)-(deg_Pep*Pep3)'

#        # Differential equations
#        dmRNA1_dt = syn_mRNA1*(state1) - (deg_mRNA *mRNA1)
#        dPep1_dt = (syn_Pep*mRNA1) - (deg_Pep*Pep1)
#        dmRNA2_dt = syn_mRNA2*(state2) - (deg_mRNA *mRNA2)
#        dPep2_dt = (syn_Pep*mRNA2) - (deg_Pep*Pep2)
#        dmRNA3_dt = (syn_mRNA3*(Pep1/Pep1max)*(Pep2/Pep2max))-(deg_mRNA * mRNA3)
#        dPep3_dt = (syn_Pep*mRNA3)-(deg_Pep*Pep3)
#
#        dmRNA1_dt_ = dmRNA1_dt.subs([(syn_mRNA1, syn_mRNA1_),(state1, state1_), \
#                                     (deg_mRNA, deg_mRNA_), (mRNA1, mRNA1_)])
#        dPep1_dt_ = dPep1_dt.subs([(syn_Pep, syn_Pep_),(mRNA1, mRNA1_), \
#                                     (deg_Pep, deg_Pep_), (Pep1, Pep1_)])
#        dmRNA2_dt_ = dmRNA2_dt.subs([(syn_mRNA2, syn_mRNA2_),(state2, state2_), \
#                                     (deg_mRNA, deg_mRNA_), (mRNA2, mRNA2_)])
#        dPep2_dt_ = dPep2_dt.subs([(syn_Pep, syn_Pep_),(mRNA2, mRNA2_), \
#                                     (deg_Pep, deg_Pep_), (Pep2, Pep2_)])
#        dmRNA3_dt_ = dmRNA3_dt.subs([(syn_mRNA3, syn_mRNA3_),(Pep1, Pep1_), \
#                                     (Pep1max, Pep1max_), (Pep2, Pep2_), \
#                                     (Pep2max, Pep2max_), (deg_mRNA, deg_mRNA_), \
#                                     (mRNA3, mRNA3_)])
#        dPep3_dt_ = dPep3_dt.subs([(syn_Pep, syn_Pep_),(mRNA3, mRNA3_), \
#                                     (deg_Pep, deg_Pep_), (Pep3, Pep3_)])

        if Operation == 'Solve':
            # Return differential equations solution
            return [eval(dmRNA1_dt), eval(dPep1_dt), eval(dmRNA2_dt), eval(dPep2_dt), eval(dmRNA3_dt), eval(dPep3_dt)]
            #return [dmRNA1_dt_, dPep1_dt_, dmRNA2_dt_, dPep2_dt_, dmRNA3_dt_, dPep3_dt_]
        elif Operation == 'GetODE':
            return [dmRNA1_dt, dPep1_dt, dmRNA2_dt, dPep2_dt, dmRNA3_dt, dPep3_dt]
        else:
            print('Error: Please Enter the correct Operation for solveODE function')

    ### ODE Model for AND gate (two INPUT, STATES: 00,01,10,11)  ###
    def solveODE_ANDgateBLeak1(y, t, state, param, Operation = 'Solve'):
        # Dependent variables
        mRNA1 = y[0] # Col 1
        Pep1 = y[1]  # Col 2
        mRNA2 = y[2]
        Pep2 = y[3]
        mRNA3 = y[4]
        Pep3 = y[5]

        state1 = state[0]
        state2 = state[1]

        # Parameters
        syn_mRNA1 = param[0]
        syn_mRNA2 = param[1]
        syn_mRNA3 = param[2]
        deg_mRNA = 0.1386 # Fixed as ln2/5(mins) = 0.1386
        syn_Pep = param[3]
        deg_Pep = param[4]
        Pep1max = param[5]
        Pep2max = param[6]
        Kleak = param[7]

        # Differential equations
        dmRNA1_dt = 'Kleak + syn_mRNA1*(state1) - (deg_mRNA *mRNA1)'
        dPep1_dt = '(syn_Pep*mRNA1) - (deg_Pep*Pep1)'
        dmRNA2_dt = 'syn_mRNA2*(state2) - (deg_mRNA *mRNA2)'
        dPep2_dt = '(syn_Pep*mRNA2) - (deg_Pep*Pep2)'
        dmRNA3_dt = '(syn_mRNA3*(Pep1/Pep1max)*(Pep2/Pep2max))-(deg_mRNA * mRNA3)'
        dPep3_dt = '(syn_Pep*mRNA3)-(deg_Pep*Pep3)'

        if Operation == 'Solve':
            # Return differential equations solution
            #return [dmRNA1_dt, dPep1_dt, dmRNA2_dt, dPep2_dt, dmRNA3_dt, dPep3_dt]
            return [eval(dmRNA1_dt), eval(dPep1_dt), eval(dmRNA2_dt), eval(dPep2_dt), eval(dmRNA3_dt), eval(dPep3_dt)]
            #return [dmRNA1_dt_, dPep1_dt_, dmRNA2_dt_, dPep2_dt_, dmRNA3_dt_, dPep3_dt_]
        elif Operation == 'GetODE':
            return [dmRNA1_dt, dPep1_dt, dmRNA2_dt, dPep2_dt, dmRNA3_dt, dPep3_dt]
        else:
            print('Error: Please Enter the correct Operation for solveODE function')

    ### ODE Model for AND gate (two INPUT, STATES: 00,01,10,11)  ###
    def solveODE_ANDgateBLeak2(y, t, state, param, Operation = 'Solve'):
        # Dependent variables
        mRNA1 = y[0] # Col 1
        Pep1 = y[1]  # Col 2
        mRNA2 = y[2]
        Pep2 = y[3]
        mRNA3 = y[4]
        Pep3 = y[5]

        state1 = state[0]
        state2 = state[1]

        # Parameters
        syn_mRNA1 = param[0]
        syn_mRNA2 = param[1]
        syn_mRNA3 = param[2]
        deg_mRNA = 0.1386 # Fixed as ln2/5(mins) = 0.1386
        syn_Pep = param[3]
        deg_Pep = param[4]
        Pep1max = param[5]
        Pep2max = param[6]
        Kleak = param[7]

        # Differential equations
        dmRNA1_dt = 'syn_mRNA1*(state1) - (deg_mRNA *mRNA1)'
        dPep1_dt = '(syn_Pep*mRNA1) - (deg_Pep*Pep1)'
        dmRNA2_dt = 'Kleak + syn_mRNA2*(state2) - (deg_mRNA *mRNA2)'
        dPep2_dt = '(syn_Pep*mRNA2) - (deg_Pep*Pep2)'
        dmRNA3_dt = '(syn_mRNA3*(Pep1/Pep1max)*(Pep2/Pep2max))-(deg_mRNA * mRNA3)'
        dPep3_dt = '(syn_Pep*mRNA3)-(deg_Pep*Pep3)'

        if Operation == 'Solve':
            # Return differential equations solution
            #return [dmRNA1_dt, dPep1_dt, dmRNA2_dt, dPep2_dt, dmRNA3_dt, dPep3_dt]
            return [eval(dmRNA1_dt), eval(dPep1_dt), eval(dmRNA2_dt), eval(dPep2_dt), eval(dmRNA3_dt), eval(dPep3_dt)]
            #return [dmRNA1_dt_, dPep1_dt_, dmRNA2_dt_, dPep2_dt_, dmRNA3_dt_, dPep3_dt_]
        elif Operation == 'GetODE':
            return [dmRNA1_dt, dPep1_dt, dmRNA2_dt, dPep2_dt, dmRNA3_dt, dPep3_dt]
        else:
            print('Error: Please Enter the correct Operation for solveODE function')

    ### ODE Model for AND gate (two INPUT, STATES: 00,01,10,11)  ###
    def solveODE_ANDgateBLeak3(y, t, state, param, Operation = 'Solve'):
        # Dependent variables
        mRNA1 = y[0] # Col 1
        Pep1 = y[1]  # Col 2
        mRNA2 = y[2]
        Pep2 = y[3]
        mRNA3 = y[4]
        Pep3 = y[5]

        state1 = state[0]
        state2 = state[1]

        # Parameters
        syn_mRNA1 = param[0]
        syn_mRNA2 = param[1]
        syn_mRNA3 = param[2]
        deg_mRNA = 0.1386 # Fixed as ln2/5(mins) = 0.1386
        syn_Pep = param[3]
        deg_Pep = param[4]
        Pep1max = param[5]
        Pep2max = param[6]
        Kleak = param[7]

        # Differential equations
        dmRNA1_dt = 'syn_mRNA1*(state1) - (deg_mRNA *mRNA1)'
        dPep1_dt = '(syn_Pep*mRNA1) - (deg_Pep*Pep1)'
        dmRNA2_dt = 'syn_mRNA2*(state2) - (deg_mRNA *mRNA2)'
        dPep2_dt = '(syn_Pep*mRNA2) - (deg_Pep*Pep2)'
        dmRNA3_dt = 'Kleak+ (syn_mRNA3*(Pep1/Pep1max)*(Pep2/Pep2max))-(deg_mRNA * mRNA3)'
        dPep3_dt = '(syn_Pep*mRNA3)-(deg_Pep*Pep3)'

        if Operation == 'Solve':
            # Return differential equations solution
            #return [dmRNA1_dt, dPep1_dt, dmRNA2_dt, dPep2_dt, dmRNA3_dt, dPep3_dt]
            return [eval(dmRNA1_dt), eval(dPep1_dt), eval(dmRNA2_dt), eval(dPep2_dt), eval(dmRNA3_dt), eval(dPep3_dt)]
            #return [dmRNA1_dt_, dPep1_dt_, dmRNA2_dt_, dPep2_dt_, dmRNA3_dt_, dPep3_dt_]
        elif Operation == 'GetODE':
            return [dmRNA1_dt, dPep1_dt, dmRNA2_dt, dPep2_dt, dmRNA3_dt, dPep3_dt]
        else:
            print('Error: Please Enter the correct Operation for solveODE function')

    ### ODE Model for AND gate (two INPUT, STATES: 00,01,10,11)  ###
    def solveODE_ANDgateBLeak13(y, t, state, param, Operation = 'Solve'):
        # Dependent variables
        mRNA1 = y[0] # Col 1
        Pep1 = y[1]  # Col 2
        mRNA2 = y[2]
        Pep2 = y[3]
        mRNA3 = y[4]
        Pep3 = y[5]

        state1 = state[0]
        state2 = state[1]

        # Parameters
        syn_mRNA1 = param[0]
        syn_mRNA2 = param[1]
        syn_mRNA3 = param[2]
        deg_mRNA = 0.1386 # Fixed as ln2/5(mins) = 0.1386
        syn_Pep = param[3]
        deg_Pep = param[4]
        Pep1max = param[5]
        Pep2max = param[6]
        Kleak = param[7]
        Kleak1 = param[8]

        # Differential equations
        dmRNA1_dt = 'Kleak1 + syn_mRNA1*(state1) - (deg_mRNA *mRNA1)'
        dPep1_dt = '(syn_Pep*mRNA1) - (deg_Pep*Pep1)'
        dmRNA2_dt = 'syn_mRNA2*(state2) - (deg_mRNA *mRNA2)'
        dPep2_dt = '(syn_Pep*mRNA2) - (deg_Pep*Pep2)'
        dmRNA3_dt = 'Kleak+ (syn_mRNA3*(Pep1/Pep1max)*(Pep2/Pep2max))-(deg_mRNA * mRNA3)'
        dPep3_dt = '(syn_Pep*mRNA3)-(deg_Pep*Pep3)'

        if Operation == 'Solve':
            # Return differential equations solution
            #return [dmRNA1_dt, dPep1_dt, dmRNA2_dt, dPep2_dt, dmRNA3_dt, dPep3_dt]
            return [eval(dmRNA1_dt), eval(dPep1_dt), eval(dmRNA2_dt), eval(dPep2_dt), eval(dmRNA3_dt), eval(dPep3_dt)]
            #return [dmRNA1_dt_, dPep1_dt_, dmRNA2_dt_, dPep2_dt_, dmRNA3_dt_, dPep3_dt_]
        elif Operation == 'GetODE':
            return [dmRNA1_dt, dPep1_dt, dmRNA2_dt, dPep2_dt, dmRNA3_dt, dPep3_dt]
        else:
            print('Error: Please Enter the correct Operation for solveODE function')

    ### ODE Model for AND gate (two INPUT, STATES: 00,01,10,11)  ###
    def solveODE_ANDgateBLeak13KMat(y, t, state, param, Operation = 'Solve'):
        # Dependent variables
        mRNA1 = y[0] # Col 1
        Pep1 = y[1]  # Col 2
        mRNA2 = y[2]
        Pep2 = y[3]
        mRNA3 = y[4]
        Pep3 = y[5]
        Pep3m = y[6]

        state1 = state[0]
        state2 = state[1]

        # Parameters
        syn_mRNA1 = param[0]
        syn_mRNA2 = param[1]
        syn_mRNA3 = param[2]
        deg_mRNA = 0.1386 # Fixed as ln2/5(mins) = 0.1386
        syn_Pep = param[3]
        deg_Pep = param[4]
        Pep1max = param[5]
        Pep2max = param[6]
        Kleak = param[7]
        Kleak1 = param[8]
        Kmature = param[9]

        # Differential equations
        dmRNA1_dt = 'Kleak1 + syn_mRNA1*(state1) - (deg_mRNA *mRNA1)'
        dPep1_dt = '(syn_Pep*mRNA1) - (deg_Pep*Pep1)'
        dmRNA2_dt = 'syn_mRNA2*(state2) - (deg_mRNA *mRNA2)'
        dPep2_dt = '(syn_Pep*mRNA2) - (deg_Pep*Pep2)'
        dmRNA3_dt = 'Kleak+ (syn_mRNA3*(Pep1/Pep1max)*(Pep2/Pep2max))-(deg_mRNA * mRNA3)'
        dPep3_dt = '(syn_Pep*mRNA3)-(Kmature*Pep3)'
        dPep3m_dt = '(Kmature*Pep3)-(deg_Pep*Pep3m)'

        if Operation == 'Solve':
            # Return differential equations solution
            #return [dmRNA1_dt, dPep1_dt, dmRNA2_dt, dPep2_dt, dmRNA3_dt, dPep3_dt]
            return [eval(dmRNA1_dt), eval(dPep1_dt), eval(dmRNA2_dt), eval(dPep2_dt), eval(dmRNA3_dt), eval(dPep3_dt), eval(dPep3m_dt)]
            #return [dmRNA1_dt_, dPep1_dt_, dmRNA2_dt_, dPep2_dt_, dmRNA3_dt_, dPep3_dt_]
        elif Operation == 'GetODE':
            return [dmRNA1_dt, dPep1_dt, dmRNA2_dt, dPep2_dt, dmRNA3_dt, dPep3_dt, dPep3m_dt]
        else:
            print('Error: Please Enter the correct Operation for solveODE function')

    ### ODE Model for OR gate (two INPUT, STATES: 00,01,10,11) ###
    def solveODE_ORgate(y, t, state, param, Operation = 'Solve'):
        # Dependent variables
        mRNA1 = y[0] # Col 1
        Pep1 = y[1]  # Col 2
        mRNA2 = y[2]
        Pep2 = y[3]
        mRNA3 = y[4]
        Pep3 = y[5]

        state1 = state[0]
        state2 = state[1]

        # Parameters
        syn_mRNA1 = param[0]
        syn_mRNA2 = param[1]
        syn_mRNA3 = param[2]
        deg_mRNA = 0.1386 # Fixed as ln2/5(mins) = 0.1386
        syn_Pep = param[3]
        deg_Pep = param[4]
        Pepmax = param[5]

        # Differential equations
        dmRNA1_dt = 'syn_mRNA1*(state1) - (deg_mRNA *mRNA1)'
        dPep1_dt = '(syn_Pep*mRNA1) - (deg_Pep*Pep1)'
        dmRNA2_dt = 'syn_mRNA2*(state2) - (deg_mRNA *mRNA2)'
        dPep2_dt = '(syn_Pep*mRNA2) - (deg_Pep*Pep2)'
        dmRNA3_dt = '(syn_mRNA3*((Pep1+Pep2)/Pepmax))-(deg_mRNA *mRNA3)'
        dPep3_dt = '(syn_Pep*mRNA3)-(deg_Pep*Pep3)'

        if Operation == 'Solve':
            # Return differential equations solution
            #return [dmRNA1_dt, dPep1_dt, dmRNA2_dt, dPep2_dt, dmRNA3_dt, dPep3_dt]
            return [eval(dmRNA1_dt), eval(dPep1_dt), eval(dmRNA2_dt), eval(dPep2_dt), eval(dmRNA3_dt), eval(dPep3_dt)]
            #return [dmRNA1_dt_, dPep1_dt_, dmRNA2_dt_, dPep2_dt_, dmRNA3_dt_, dPep3_dt_]
        elif Operation == 'GetODE':
            return [dmRNA1_dt, dPep1_dt, dmRNA2_dt, dPep2_dt, dmRNA3_dt, dPep3_dt]
        else:
            print('Error: Please Enter the correct Operation for solveODE function')


    def solveODE_ORgateDelay(y, t, state, param, Operation = 'Solve'):
        # Dependent variables
        Inde = y[0]
        Indi = y[1]
        mRNA1 = y[2] # Col 1
        Pep1 = y[3]  # Col 2
        mRNA2 = y[4]
        Pep2 = y[5]
        mRNA3 = y[6]
        Pep3 = y[7]

        state1 = state[0]
        state2 = state[1]

        # Parameters
        syn_mRNA1 = param[0]
        syn_mRNA2 = param[1]
        syn_mRNA3 = param[2]
        deg_mRNA = 0.1386 # Fixed as ln2/5(mins) = 0.1386
        syn_Pep = param[3]
        deg_Pep = param[4]
        Pepmax = param[5]
        Km = param[6]

        # Differential equations
        dInde_dt = '-(Inde/(Inde+Km))*Inde'
        dIndi_dt = '(Inde/(Inde+Km))*Inde'
        dmRNA1_dt = 'syn_mRNA1*(Indi)*(state1) - (deg_mRNA *mRNA1)'
        dPep1_dt = '(syn_Pep*mRNA1) - (deg_Pep*Pep1)'
        dmRNA2_dt = 'syn_mRNA2*(state2) - (deg_mRNA *mRNA2)'
        dPep2_dt = '(syn_Pep*mRNA2) - (deg_Pep*Pep2)'
        dmRNA3_dt = '(syn_mRNA3*((Pep1+Pep2)/Pepmax))-(deg_mRNA *mRNA3)'
        dPep3_dt = '(syn_Pep*mRNA3)-(deg_Pep*Pep3)'

        if Operation == 'Solve':
            # Return differential equations solution
            #return [dInde1_dt, dIndi1_dt, dInd2_dt, dmRNA1_dt, dPep1_dt, dmRNA2_dt, dPep2_dt, dmRNA3_dt, dPep3_dt]
            return [eval(dInde_dt), eval(dIndi_dt), eval(dmRNA1_dt), \
                    eval(dPep1_dt), eval(dmRNA2_dt), eval(dPep2_dt), eval(dmRNA3_dt), eval(dPep3_dt)]
            #return [dmRNA1_dt_, dPep1_dt_, dmRNA2_dt_, dPep2_dt_, dmRNA3_dt_, dPep3_dt_]
        elif Operation == 'GetODE':
            return [dInde_dt, dIndi_dt, dmRNA1_dt, dPep1_dt, dmRNA2_dt, dPep2_dt, dmRNA3_dt, dPep3_dt]
        else:
            print('Error: Please Enter the correct Operation for solveODE function')

    def solveODE_ORgate_Delay(y, t, state, param, Operation = 'Solve'):
        # Dependent variables
        Inde = y[0]
        Indi = y[1]
        mRNA1 = y[2] # Col 1
        Pep1 = y[3]  # Col 2
        mRNA2 = y[4]
        Pep2 = y[5]
        mRNA3 = y[6]
        Pep3 = y[7]

        state1 = state[0]
        state2 = state[1]

        # Parameters
        syn_mRNA1 = param[0]
        syn_mRNA2 = param[1]
        syn_mRNA3 = param[2]
        deg_mRNA = 0.1386 # Fixed as ln2/5(mins) = 0.1386
        syn_Pep = param[3]
        deg_Pep = param[4]
        Pepmax = param[5]
        Km = param[6]


        # Differential equations
        dInde_dt = '-(Inde/(Inde+Km))*Inde'
        dIndi_dt = '(Inde/(Inde+Km))*Inde'
        dmRNA1_dt = 'syn_mRNA1*(state1) - (deg_mRNA *mRNA1)'
        dPep1_dt = '(syn_Pep*mRNA1) - (deg_Pep*Pep1)'
        dmRNA2_dt = 'syn_mRNA2*(Indi)*(state2) - (deg_mRNA *mRNA2)'
        dPep2_dt = '(syn_Pep*mRNA2) - (deg_Pep*Pep2)'
        dmRNA3_dt = '(syn_mRNA3*((Pep1+Pep2)/Pepmax))-(deg_mRNA *mRNA3)'
        dPep3_dt = '(syn_Pep*mRNA3)-(deg_Pep*Pep3)'

        if Operation == 'Solve':
            # Return differential equations solution
            return [eval(dInde_dt), eval(dIndi_dt), eval(dmRNA1_dt), \
                    eval(dPep1_dt), eval(dmRNA2_dt), eval(dPep2_dt), eval(dmRNA3_dt), eval(dPep3_dt)]
            #return [dmRNA1_dt_, dPep1_dt_, dmRNA2_dt_, dPep2_dt_, dmRNA3_dt_, dPep3_dt_]
        elif Operation == 'GetODE':
            return [dInde_dt, dIndi_dt, dmRNA1_dt, dPep1_dt, dmRNA2_dt, dPep2_dt, dmRNA3_dt, dPep3_dt]
        else:
            print('Error: Please Enter the correct Operation for solveODE function')

    def solveODE_ORgateDegradation(y, t, state, param, Operation = 'Solve'):
        # Dependent variables
        Ind = y[0]
        mRNA1 = y[1] # Col 1
        Pep1 = y[2]  # Col 2
        mRNA2 = y[3]
        Pep2 = y[4]
        mRNA3 = y[5]
        Pep3 = y[6]

        state1 = state[0]
        state2 = state[1]

        # Parameters
        syn_mRNA1 = param[0]
        syn_mRNA2 = param[1]
        syn_mRNA3 = param[2]
        deg_mRNA = 0.1386 # Fixed as ln2/5(mins) = 0.1386
        syn_Pep = param[3]
        deg_Pep = param[4]
        Pepmax = param[5]
        deg_Ind = param[6]

        # Differential equations
        dInd_dt = '-deg_Ind*Ind'
        dmRNA1_dt = 'syn_mRNA1*(Ind)*(state1) - (deg_mRNA *mRNA1)'
        dPep1_dt = '(syn_Pep*mRNA1) - (deg_Pep*Pep1)'
        dmRNA2_dt = 'syn_mRNA2*(state2) - (deg_mRNA *mRNA2)'
        dPep2_dt = '(syn_Pep*mRNA2) - (deg_Pep*Pep2)'
        dmRNA3_dt = '(syn_mRNA3*((Pep1+Pep2)/Pepmax))-(deg_mRNA *mRNA3)'
        dPep3_dt = '(syn_Pep*mRNA3)-(deg_Pep*Pep3)'

        if Operation == 'Solve':
            # Return differential equations solution
            return [eval(dInd_dt), eval(dmRNA1_dt), \
                    eval(dPep1_dt), eval(dmRNA2_dt), eval(dPep2_dt), eval(dmRNA3_dt), eval(dPep3_dt)]
            #return [dmRNA1_dt_, dPep1_dt_, dmRNA2_dt_, dPep2_dt_, dmRNA3_dt_, dPep3_dt_]
        elif Operation == 'GetODE':
            return [dInd_dt, dmRNA1_dt, dPep1_dt, dmRNA2_dt, dPep2_dt, dmRNA3_dt, dPep3_dt]
        else:
            print('Error: Please Enter the correct Operation for solveODE function')

    def solveODE_ORgate_Degradation(y, t, state, param, Operation = 'Solve'):
        # Dependent variables
        Ind = y[0]
        mRNA1 = y[1] # Col 1
        Pep1 = y[2]  # Col 2
        mRNA2 = y[3]
        Pep2 = y[4]
        mRNA3 = y[5]
        Pep3 = y[6]

        state1 = state[0]
        state2 = state[1]

        # Parameters
        syn_mRNA1 = param[0]
        syn_mRNA2 = param[1]
        syn_mRNA3 = param[2]
        deg_mRNA = 0.1386 # Fixed as ln2/5(mins) = 0.1386
        syn_Pep = param[3]
        deg_Pep = param[4]
        Pepmax = param[5]
        deg_Ind = param[6]

        # Differential equations
        dInd_dt = '-deg_Ind*Ind'
        dmRNA1_dt = 'syn_mRNA1*(state1) - (deg_mRNA *mRNA1)'
        dPep1_dt = '(syn_Pep*mRNA1) - (deg_Pep*Pep1)'
        dmRNA2_dt = 'syn_mRNA2*Ind*(state2) - (deg_mRNA *mRNA2)'
        dPep2_dt = '(syn_Pep*mRNA2) - (deg_Pep*Pep2)'
        dmRNA3_dt = '(syn_mRNA3*((Pep1+Pep2)/Pepmax))-(deg_mRNA *mRNA3)'
        dPep3_dt = '(syn_Pep*mRNA3)-(deg_Pep*Pep3)'

        if Operation == 'Solve':
            # Return differential equations solution
            return [eval(dInd_dt), eval(dmRNA1_dt), \
                    eval(dPep1_dt), eval(dmRNA2_dt), eval(dPep2_dt), eval(dmRNA3_dt), eval(dPep3_dt)]
            #return [dmRNA1_dt_, dPep1_dt_, dmRNA2_dt_, dPep2_dt_, dmRNA3_dt_, dPep3_dt_]
        elif Operation == 'GetODE':
            return [dInd_dt, dmRNA1_dt, dPep1_dt, dmRNA2_dt, dPep2_dt, dmRNA3_dt, dPep3_dt]
        else:
            print('Error: Please Enter the correct Operation for solveODE function')

    def solveODE_ORgateDelayDegradation(y, t, state, param, Operation = 'Solve'):
        # Dependent variables
        Inde = y[0]
        Indi = y[1]
        Ind = y[2]
        mRNA1 = y[3] # Col 1
        Pep1 = y[4]  # Col 2
        mRNA2 = y[5]
        Pep2 = y[6]
        mRNA3 = y[7]
        Pep3 = y[8]

        state1 = state[0]
        state2 = state[1]

        # Parameters
        syn_mRNA1 = param[0]
        syn_mRNA2 = param[1]
        syn_mRNA3 = param[2]
        deg_mRNA = 0.1386 # Fixed as ln2/5(mins) = 0.1386
        syn_Pep = param[3]
        deg_Pep = param[4]
        Pepmax = param[5]
        Km = param[6]
        deg_Ind = param[7]

        # Differential equations
        dInde_dt = '-(Inde/(Inde+Km))*Inde'
        dIndi_dt = '(Inde/(Inde+Km))*Inde'
        dInd_dt = '-deg_Ind*Ind'
        dmRNA1_dt = 'syn_mRNA1*(Indi)*(state1) - (deg_mRNA *mRNA1)'
        dPep1_dt = '(syn_Pep*mRNA1) - (deg_Pep*Pep1)'
        dmRNA2_dt = 'syn_mRNA2*(Ind)*(state2) - (deg_mRNA *mRNA2)'
        dPep2_dt = '(syn_Pep*mRNA2) - (deg_Pep*Pep2)'
        dmRNA3_dt = '(syn_mRNA3*((Pep1+Pep2)/Pepmax))-(deg_mRNA *mRNA3)'
        dPep3_dt = '(syn_Pep*mRNA3)-(deg_Pep*Pep3)'

        if Operation == 'Solve':
            # Return differential equations solution
            #return [dInde1_dt, dIndi1_dt, dInd2_dt, dmRNA1_dt, dPep1_dt, dmRNA2_dt, dPep2_dt, dmRNA3_dt, dPep3_dt]
            return [eval(dInde_dt), eval(dIndi_dt), eval(dInd_dt), eval(dmRNA1_dt), \
                    eval(dPep1_dt), eval(dmRNA2_dt), eval(dPep2_dt), eval(dmRNA3_dt), eval(dPep3_dt)]
            #return [dmRNA1_dt_, dPep1_dt_, dmRNA2_dt_, dPep2_dt_, dmRNA3_dt_, dPep3_dt_]
        elif Operation == 'GetODE':
            return [dInde_dt, dIndi_dt, dInd_dt, dmRNA1_dt, dPep1_dt, dmRNA2_dt, dPep2_dt, dmRNA3_dt, dPep3_dt]
        else:
            print('Error: Please Enter the correct Operation for solveODE function')

    def solveODE_ORgateDegradationDelay(y, t, state, param, Operation = 'Solve'):
        # Dependent variables
        Inde = y[0]
        Indi = y[1]
        Ind = y[2]
        mRNA1 = y[3] # Col 1
        Pep1 = y[4]  # Col 2
        mRNA2 = y[5]
        Pep2 = y[6]
        mRNA3 = y[7]
        Pep3 = y[8]

        state1 = state[0]
        state2 = state[1]

        # Parameters
        syn_mRNA1 = param[0]
        syn_mRNA2 = param[1]
        syn_mRNA3 = param[2]
        deg_mRNA = 0.1386 # Fixed as ln2/5(mins) = 0.1386
        syn_Pep = param[3]
        deg_Pep = param[4]
        Pepmax = param[5]
        Km = param[6]
        deg_Ind = param[7]

        # Differential equations
        dInd_dt = '-deg_Ind*Ind'
        dInde_dt = '-(Inde/(Inde+Km))*Inde'
        dIndi_dt = '(Inde/(Inde+Km))*Inde'
        dmRNA1_dt = 'syn_mRNA1*(Ind)*(state1) - (deg_mRNA *mRNA1)'
        dPep1_dt = '(syn_Pep*mRNA1) - (deg_Pep*Pep1)'
        dmRNA2_dt = 'syn_mRNA2*(Indi)*(state2) - (deg_mRNA *mRNA2)'
        dPep2_dt = '(syn_Pep*mRNA2) - (deg_Pep*Pep2)'
        dmRNA3_dt = '(syn_mRNA3*((Pep1+Pep2)/Pepmax))-(deg_mRNA *mRNA3)'
        dPep3_dt = '(syn_Pep*mRNA3)-(deg_Pep*Pep3)'

        if Operation == 'Solve':
            # Return differential equations solution
            return [eval(dInde_dt), eval(dIndi_dt), eval(dInd_dt), eval(dmRNA1_dt), \
                    eval(dPep1_dt), eval(dmRNA2_dt), eval(dPep2_dt), eval(dmRNA3_dt), eval(dPep3_dt)]
            #return [dmRNA1_dt_, dPep1_dt_, dmRNA2_dt_, dPep2_dt_, dmRNA3_dt_, dPep3_dt_]
        elif Operation == 'GetODE':
            return [dInde_dt, dIndi_dt, dInd_dt, dmRNA1_dt, dPep1_dt, dmRNA2_dt, dPep2_dt, dmRNA3_dt, dPep3_dt]
        else:
            print('Error: Please Enter the correct Operation for solveODE function')


    def solveODE_ORgateDelayDelay(y, t, state, param, Operation = 'Solve'):
        # Dependent variables
        Inde1 = y[0]
        Indi1 = y[1]
        Inde2 = y[2]
        Indi2 = y[3]
        mRNA1 = y[4] # Col 1
        Pep1 = y[5]  # Col 2
        mRNA2 = y[6]
        Pep2 = y[7]
        mRNA3 = y[8]
        Pep3 = y[9]

        state1 = state[0]
        state2 = state[1]

        # Parameters
        syn_mRNA1 = param[0]
        syn_mRNA2 = param[1]
        syn_mRNA3 = param[2]
        deg_mRNA = 0.1386 # Fixed as ln2/5(mins) = 0.1386
        syn_Pep = param[3]
        deg_Pep = param[4]
        Pepmax = param[5]
        Km1 = param[6]
        Km2 = param[7]

        # Differential equations
        dInde1_dt = '-(Inde1/(Inde1+Km1))*Inde1'
        dIndi1_dt = '(Inde1/(Inde1+Km1))*Inde1'
        dInde2_dt = '-(Inde2/(Inde2+Km2))*Inde2'
        dIndi2_dt = '(Inde2/(Inde2+Km2))*Inde2'
        dmRNA1_dt = 'syn_mRNA1*(Indi1)*(state1) - (deg_mRNA *mRNA1)'
        dPep1_dt = '(syn_Pep*mRNA1) - (deg_Pep*Pep1)'
        dmRNA2_dt = 'syn_mRNA2*(Indi2)*(state2) - (deg_mRNA *mRNA2)'
        dPep2_dt = '(syn_Pep*mRNA2) - (deg_Pep*Pep2)'
        dmRNA3_dt = '(syn_mRNA3*((Pep1+Pep2)/Pepmax))-(deg_mRNA *mRNA3)'
        dPep3_dt = '(syn_Pep*mRNA3)-(deg_Pep*Pep3)'

        if Operation == 'Solve':
            # Return differential equations solution
            return [eval(dInde1_dt), eval(dIndi1_dt), eval(dInde2_dt), eval(dIndi2_dt), eval(dmRNA1_dt), \
                    eval(dPep1_dt), eval(dmRNA2_dt), eval(dPep2_dt), eval(dmRNA3_dt), eval(dPep3_dt)]
            #return [dmRNA1_dt_, dPep1_dt_, dmRNA2_dt_, dPep2_dt_, dmRNA3_dt_, dPep3_dt_]
        elif Operation == 'GetODE':
            return [dInde1_dt, dIndi1_dt, dInde2_dt, dIndi2_dt, dmRNA1_dt, dPep1_dt, dmRNA2_dt, dPep2_dt, dmRNA3_dt, dPep3_dt]
        else:
            print('Error: Please Enter the correct Operation for solveODE function')


    def solveODE_ORgateDelayDegradeResCompete(y, t, state, param, Operation = 'Solve'):
        # Dependent variables
        Inde = y[0]
        Indi = y[1]
        Ind = y[2]
        mRNA1 = y[3] # Col 1
        Pep1 = y[4]  # Col 2
        mRNA2 = y[5]
        Pep2 = y[6]
        mRNA3 = y[7]
        Pep3 = y[8]

        state1 = state[0]
        state2 = state[1]

        # Parameters
        syn_mRNA1 = param[0]
        syn_mRNA2 = param[1]
        syn_mRNA3 = param[2]
        deg_mRNA = 0.1386 # Fixed as ln2/5(mins) = 0.1386
        syn_Pep = param[3]
        deg_Pep = param[4]
        Pepmax = param[5]
        Km = param[6]
        deg_Ind = param[7]

        if (state1 == 1) and (state2 == 1):
            syn_Pep3 = syn_Pep*param[8]
        else:
            syn_Pep3 = syn_Pep*1

        # Differential equations
        dInde_dt = '-(Inde/(Inde+Km))*Inde'
        dIndi_dt = '(Inde/(Inde+Km))*Inde'
        dInd_dt = '-deg_Ind*Ind'
        dmRNA1_dt = 'syn_mRNA1*(Indi)*(state1) - (deg_mRNA *mRNA1)'
        dPep1_dt = '(syn_Pep*mRNA1) - (deg_Pep*Pep1)'
        dmRNA2_dt = 'syn_mRNA2*(Ind)*(state2) - (deg_mRNA *mRNA2)'
        dPep2_dt = '(syn_Pep*mRNA2) - (deg_Pep*Pep2)'
        dmRNA3_dt = '(syn_mRNA3*((Pep1+Pep2)/Pepmax))-(deg_mRNA *mRNA3)'
        dPep3_dt = '(syn_Pep3*mRNA3)-(deg_Pep*Pep3)'

        if Operation == 'Solve':
            # Return differential equations solution
            #return [dInde1_dt, dIndi1_dt, dInd2_dt, dmRNA1_dt, dPep1_dt, dmRNA2_dt, dPep2_dt, dmRNA3_dt, dPep3_dt]
            return [eval(dInde_dt), eval(dIndi_dt), eval(dInd_dt), eval(dmRNA1_dt), \
                    eval(dPep1_dt), eval(dmRNA2_dt), eval(dPep2_dt), eval(dmRNA3_dt), eval(dPep3_dt)]
            #return [dmRNA1_dt_, dPep1_dt_, dmRNA2_dt_, dPep2_dt_, dmRNA3_dt_, dPep3_dt_]
        elif Operation == 'GetODE':
            return [dInde_dt, dIndi_dt, dInd_dt, dmRNA1_dt, dPep1_dt, dmRNA2_dt, dPep2_dt, dmRNA3_dt, dPep3_dt]
        else:
            print('Error: Please Enter the correct Operation for solveODE function')

    def solveODE_ORgateDegradeDelayResCompete(y, t, state, param, Operation = 'Solve'):
        # Dependent variables
        Inde = y[0]
        Indi = y[1]
        Ind = y[2]
        mRNA1 = y[3] # Col 1
        Pep1 = y[4]  # Col 2
        mRNA2 = y[5]
        Pep2 = y[6]
        mRNA3 = y[7]
        Pep3 = y[8]

        state1 = state[0]
        state2 = state[1]

        # Parameters
        syn_mRNA1 = param[0]
        syn_mRNA2 = param[1]
        syn_mRNA3 = param[2]
        deg_mRNA = 0.1386 # Fixed as ln2/5(mins) = 0.1386
        syn_Pep = param[3]
        deg_Pep = param[4]
        Pepmax = param[5]
        Km = param[6]
        deg_Ind = param[7]

        if (state1 == 1) and (state2 == 1):
            syn_Pep3 = syn_Pep*param[8]
        else:
            syn_Pep3 = syn_Pep*1

        # Differential equations
        dInd_dt = '-deg_Ind*Ind'
        dInde_dt = '-(Inde/(Inde+Km))*Inde'
        dIndi_dt = '(Inde/(Inde+Km))*Inde'
        dmRNA1_dt = 'syn_mRNA1*(Ind)*(state1) - (deg_mRNA *mRNA1)'
        dPep1_dt = '(syn_Pep*mRNA1) - (deg_Pep*Pep1)'
        dmRNA2_dt = 'syn_mRNA2*(Indi)*(state2) - (deg_mRNA *mRNA2)'
        dPep2_dt = '(syn_Pep*mRNA2) - (deg_Pep*Pep2)'
        dmRNA3_dt = '(syn_mRNA3*((Pep1+Pep2)/Pepmax))-(deg_mRNA *mRNA3)'
        dPep3_dt = '(syn_Pep3*mRNA3)-(deg_Pep*Pep3)'

        if Operation == 'Solve':
            # Return differential equations solution
            #return [dInde1_dt, dIndi1_dt, dInd2_dt, dmRNA1_dt, dPep1_dt, dmRNA2_dt, dPep2_dt, dmRNA3_dt, dPep3_dt]
            return [eval(dInde_dt), eval(dIndi_dt), eval(dInd_dt), eval(dmRNA1_dt), \
                    eval(dPep1_dt), eval(dmRNA2_dt), eval(dPep2_dt), eval(dmRNA3_dt), eval(dPep3_dt)]
            #return [dmRNA1_dt_, dPep1_dt_, dmRNA2_dt_, dPep2_dt_, dmRNA3_dt_, dPep3_dt_]
        elif Operation == 'GetODE':
            return [dInde_dt, dIndi_dt, dInd_dt, dmRNA1_dt, dPep1_dt, dmRNA2_dt, dPep2_dt, dmRNA3_dt, dPep3_dt]
        else:
            print('Error: Please Enter the correct Operation for solveODE function')

    def solveODE_ORgateDelayDelayResCompete(y, t, state, param, Operation = 'Solve'):
        # Dependent variables
        Inde1 = y[0]
        Indi1 = y[1]
        Inde2 = y[2]
        Indi2 = y[3]
        mRNA1 = y[4] # Col 1
        Pep1 = y[5]  # Col 2
        mRNA2 = y[6]
        Pep2 = y[7]
        mRNA3 = y[8]
        Pep3 = y[9]

        state1 = state[0]
        state2 = state[1]

        # Parameters
        syn_mRNA1 = param[0]
        syn_mRNA2 = param[1]
        syn_mRNA3 = param[2]
        deg_mRNA = 0.1386 # Fixed as ln2/5(mins) = 0.1386
        syn_Pep = param[3]
        deg_Pep = param[4]
        Pepmax = param[5]
        Km1 = param[6]
        Km2 = param[7]

        if (state1 == 1) and (state2 == 1):
            syn_Pep3 = syn_Pep*param[8]
        else:
            syn_Pep3 = syn_Pep*1

        # Differential equations
        dInde1_dt = '-(Inde1/(Inde1+Km1))*Inde1'
        dIndi1_dt = '(Inde1/(Inde1+Km1))*Inde1'
        dInde2_dt = '-(Inde2/(Inde2+Km2))*Inde2'
        dIndi2_dt = '(Inde2/(Inde2+Km2))*Inde2'
        dmRNA1_dt = 'syn_mRNA1*(Indi1)*(state1) - (deg_mRNA *mRNA1)'
        dPep1_dt = '(syn_Pep*mRNA1) - (deg_Pep*Pep1)'
        dmRNA2_dt = 'syn_mRNA2*(Indi2)*(state2) - (deg_mRNA *mRNA2)'
        dPep2_dt = '(syn_Pep*mRNA2) - (deg_Pep*Pep2)'
        dmRNA3_dt = '(syn_mRNA3*((Pep1+Pep2)/Pepmax))-(deg_mRNA *mRNA3)'
        dPep3_dt = '(syn_Pep3*mRNA3)-(deg_Pep*Pep3)'

        if Operation == 'Solve':
            # Return differential equations solution
            #return [dInde1_dt, dIndi1_dt, dInd2_dt, dmRNA1_dt, dPep1_dt, dmRNA2_dt, dPep2_dt, dmRNA3_dt, dPep3_dt]
            return [eval(dInde1_dt), eval(dIndi1_dt), eval(dInde2_dt), eval(dIndi2_dt), eval(dmRNA1_dt), \
                    eval(dPep1_dt), eval(dmRNA2_dt), eval(dPep2_dt), eval(dmRNA3_dt), eval(dPep3_dt)]
            #return [dmRNA1_dt_, dPep1_dt_, dmRNA2_dt_, dPep2_dt_, dmRNA3_dt_, dPep3_dt_]
        elif Operation == 'GetODE':
            return [dInde1_dt, dIndi1_dt, dInde2_dt, dIndi2_dt, dmRNA1_dt, dPep1_dt, dmRNA2_dt, dPep2_dt, dmRNA3_dt, dPep3_dt]
        else:
            print('Error: Please Enter the correct Operation for solveODE function')


    ### To map the string name with function name
    function_mappings = {
        'solveODE_NOTgate': solveODE_NOTgate,
        'solveODE_NOTgateKMat': solveODE_NOTgateKMat,
        'solveODE_NOTgateSingle': solveODE_NOTgateSingle,
        'solveODE_NOTgateSingleKMat': solveODE_NOTgateSingleKMat,
        'solveODE_ORgate': solveODE_ORgate,
        'solveODE_ORgateDelay': solveODE_ORgateDelay,
        'solveODE_ORgate_Delay': solveODE_ORgate_Delay,
        'solveODE_ORgateDegradation': solveODE_ORgateDegradation,
        'solveODE_ORgate_Degradation': solveODE_ORgate_Degradation,
        'solveODE_ORgateDelayDegradation': solveODE_ORgateDelayDegradation,
        'solveODE_ORgateDegradationDelay': solveODE_ORgateDegradationDelay,
        'solveODE_ORgateDelayDelay': solveODE_ORgateDelayDelay,
        'solveODE_ORgateDelayDegradeResCompete': solveODE_ORgateDelayDegradeResCompete,
        'solveODE_ORgateDegradeDelayResCompete': solveODE_ORgateDegradeDelayResCompete,
        'solveODE_ORgateDelayDelayResCompete': solveODE_ORgateDelayDelayResCompete,
        'solveODE_ANDgate': solveODE_ANDgate,
        'solveODE_ANDgateBLeak1': solveODE_ANDgateBLeak1,
        'solveODE_ANDgateBLeak2': solveODE_ANDgateBLeak2,
        'solveODE_ANDgateBLeak3': solveODE_ANDgateBLeak3,
        'solveODE_ANDgateBLeak13': solveODE_ANDgateBLeak13,
        'solveODE_ANDgateBLeak13KMat': solveODE_ANDgateBLeak13KMat,


    }

    def select_function(self, SystemTypeString):
        ### Convert SystemType from string to function name
        try:
            return self.function_mappings[SystemTypeString]
        except KeyError:
            print('Invalid function, try again.')

    def ComputeSSE(self, param, y0, numState, rfp_data, time_int, SystemType, VarIndex, OptimizerType):
        ### Calculate SSE - To be minimized by Optimizer ####

        global weight_global, sse_global

        # Time grid == no. of times to report solution
        rfp_data_numrows = np.size(rfp_data, 0)
        rfp_data_numcols = np.size(rfp_data, 1)

        t_start = rfp_data[0][0] # First value of Time column
        t_end = rfp_data[0][-1] # Last value of Time column
        dt = 1  # minutes
        timestep = int((t_end / dt) + 1)
        t = np.linspace(t_start, t_end, timestep)

        # Initialize mRNA and Pep nested list
        #mRNA = np.zeros((timestep, numInd), dtype=object)  # timestep (rows) x numInd (cols)
        Pep = np.zeros((timestep, numState), dtype=object)   # timestep (rows) x numInd (cols)

        solveODE_Name = '_'.join(('solveODE', SystemType))
        solveODEfun = self.select_function(solveODE_Name) #convert string to function name

        if OptimizerType == 'Global':
            # Integrate the ODE equations
            if numState == 2:
                state = [0, 1]
                for i in range(0, numState):  # Iterates through Ind Concs
                    ODEsoln = odeint(solveODEfun, y0, t, args=(state[i], param))
            elif numState == 4:
                state = [[0, 0], [0, 1], [1, 0], [1, 1]]
                for i in range(0, numState):  # Iterates through Ind Concs
                    ODEsoln = odeint(solveODEfun, y0, t, args=(state[i], param))

        elif OptimizerType == 'Local':
            # Integrate the ODE equations
            # Integrate the ODE equations
            if numState == 2:
                state = [0, 1]
            elif numState == 4:
                state = [[0, 0], [0, 1], [1, 0], [1, 1]]

            for i in range(0, numState):  # Iterates through number of States
                ODEsoln = odeint(solveODEfun, y0, t, args=(state[i], param))
                for j in range(0, timestep):  # Iterates through timesteps
                    #mRNA[j][i] = ODEsoln[j][VarIndex[0]] # mRNA array runs downwards with time
                    Pep[j][i] = ODEsoln[j][VarIndex[0]]


        else:
            print('No specified Optimizer Type')
        '''
        Calculate SSE_Time
        '''
        # rfp_data runs lengthwise with time
        sse_time = 0
        for i in range(1, rfp_data_numrows):  # Start from 1 because Row 0 is time
            for j in range(0, rfp_data_numcols):
                # time_int * i to find the corresponding model Pep value
                sse_time = sse_time + (Pep[int(time_int/dt)*j, i-1] - rfp_data[i][j])**2

        sse = sse_time

        # update and return to main function
        sse_global = sse
        if OptimizerType == 'Global':
            print('Model: ', SystemType, '- SSE (Global):', sse)
        elif OptimizerType == 'Local':
            print('Model: ', SystemType, '- SSE (Local):', sse)
        else:
            print('Error in ComputeSSE function')
        return sse

    def Run_LogicGatesSystem(self, SystemType, data_header, data_array, numState):

        Time_interval = data_array[0][1] - data_array[0][0]

        ### ODE Input (Inducer Degradation) ###
        # Initial conditions for (mRNA, Pep, Ind) at time = 0
        if (SystemType == 'NOTgate'):
            mRNA10 = 0.
            Pep10 = 0.
            mRNA20 = 0.
            Pep20 = data_array[1][0]
            y0 = [mRNA10, Pep10, mRNA20, Pep20]  # Initial condition
            VariableName = ['mRNA1', 'Pep1', 'mRNA2', 'Pep2'] # Variables Name
            VarIndex =[VariableName.index('Pep2')]  #get the index for mRNA and RFP

            ### Number of Parameters to be optimized
            numParam = 6

            ParamName = ['syn_mRNA1','syn_mRNA2','syn_Pep','deg_Pep','Kmaxrep','Pepmax']
            ParamUnits = ['molL-1min-1', 'molL-1min-1', 'min-1', 'min-1', 'dimensionless', 'molL-1']

            # Global bounds - In order of [syn_mRNA1, syn_mRNA2, syn_Pep, deg_Pep, Kmrep, Pepmax]
            param0_global = [(1e-9, 5e-6), (1e-9, 5e-6), (0, 0.02), (0.001, 0.02), (0.1, 1), (5e-8, 5e-5)] #Diff Evo

            # Local bounds
            LB = [0]*(numParam-3) + [0.001, 0.1, 0]
            UB = [None]*(numParam-2) + [1, None]

        elif (SystemType == 'NOTgateKMat'):
            mRNA10 = 0.
            Pep10 = 0.
            mRNA20 = 0.
            Pep20 = 0.
            Pep2m0 = data_array[1][0]
            y0 = [mRNA10, Pep10, mRNA20, Pep20, Pep2m0]  # Initial condition
            VariableName = ['mRNA1', 'Pep1', 'mRNA2', 'Pep2', 'Pep2m'] # Variables Name
            VarIndex =[VariableName.index('Pep2m')]  #get the index for mRNA and RFP

            ### Number of Parameters to be optimized
            numParam = 7

            ParamName = ['syn_mRNA1','syn_mRNA2','syn_Pep','deg_Pep','Kmaxrep','Pepmax', 'Kmature']
            ParamUnits = ['molL-1min-1', 'molL-1min-1', 'min-1', 'min-1', 'dimensionless', 'molL-1', 'min-1']

            # Global bounds - In order of [syn_mRNA1, syn_mRNA2, syn_Pep, deg_Pep, Kmrep, Pepmax]
            param0_global = [(5e-8, 5e-6), (5e-8, 5e-6), (0, 0.02), (0.001, 0.02), (0.1, 1), (5e-8, 5e-6), (0.001, 1)] #Diff Evo

            # Local bounds
            LB = [0]*(numParam-4) + [0.001, 0.1, 0, 0.001]
            UB = [None]*(numParam-3) + [1, None, 1]
            
#            #when fixing Kmature
#            LB = [0]*(numParam-3) + [0.001, 0.1, 0]
#            UB = [None]*(numParam-2) + [1, None]

        elif (SystemType == 'NOTgateSingle'):
            mRNA10 = 0.
            Pep10 = data_array[1][0]
            y0 = [mRNA10, Pep10]  # Initial condition
            VariableName = ['mRNA1', 'Pep1'] # Variables Name
            VarIndex =[VariableName.index('Pep1')]  #get the index for mRNA and RFP

            ### Number of Parameters to be optimized
            numParam = 4

            ParamName = ['syn_mRNA', 'syn_Pep','deg_Pep','Kmaxrep']
            ParamUnits = ['molL-1min-1', 'min-1', 'min-1', 'dimensionless']

            # Global bounds - In order of [syn_mRNA, syn_Pep, deg_Pep, Kmrep]
            param0_global = [(5e-8, 5e-6), (0, 0.02), (0.001, 0.02), (0.1, 1)] #Diff Evo

            # Local bounds
            LB = [0]*(numParam-2) + [0.001, 0.1]
            UB = [None]*(numParam-1) + [1]

        elif (SystemType == 'NOTgateSingleKMat'):
            mRNA10 = 0.
            Pep10 = 0.
            Pep1m0 = data_array[1][0]
            y0 = [mRNA10, Pep10, Pep1m0]  # Initial condition
            VariableName = ['mRNA1', 'Pep1', 'Pep1m'] # Variables Name
            VarIndex =[VariableName.index('Pep1m')]  #get the index for mRNA and RFP

            ### Number of Parameters to be optimized
            numParam = 5

            ParamName = ['syn_mRNA', 'syn_Pep','deg_Pep','Kmaxrep', 'Kmature']
            ParamUnits = ['molL-1min-1', 'min-1', 'min-1', 'dimensionless', 'min-1']

            # Global bounds - In order of [syn_mRNA, syn_Pep, deg_Pep, Kmrep]
            param0_global = [(5e-8, 5e-6), (0, 0.02), (0.001, 0.02), (0.1, 1), (0.001, 1)] #Diff Evo

            # Local bounds
            LB = [0]*(numParam-3) + [0.001, 0.1, 0.001]
            UB = [None]*(numParam-2) + [1, 1]

        elif (SystemType == 'ANDgate') or (SystemType == 'ANDgateBLeak3') or (SystemType == 'ANDgateBLeak1') or (SystemType == 'ANDgateBLeak2') or (SystemType == 'ANDgateBLeak13'):
            mRNA10 = 0.
            Pep10 = 0.
            mRNA20 = 0.
            Pep20 = 0.
            mRNA30 = 0.
            Pep30 = data_array[1][0]
            y0 = [mRNA10, Pep10, mRNA20, Pep20, mRNA30, Pep30]  # Initial condition
            VariableName = ['mRNA1', 'Pep1', 'mRNA2', 'Pep2', 'mRNA3', 'Pep3'] # Variables Name
            VarIndex =[VariableName.index('Pep3')]  #get the index for mRNA and RFP

            if (SystemType == 'ANDgateBLeak13'):
                ### Number of Parameters to be optimized
                numParam = 9
                ParamName = ['syn_mRNA1','syn_mRNA2','syn_mRNA3','syn_Pep','deg_Pep','Pep1max','Pep2max','Kleak','Kleak1']
                ParamUnits = ['molL-1min-1', 'molL-1min-1', 'molL-1min-1', 'min-1', 'min-1', 'molL-1', 'molL-1', 'molL-1min-1', 'molL-1min-1']
                # Global bounds - In order of [syn_mRNA1, syn_mRNA2, syn_mRNA3, syn_Pep, deg_Pep, Pep1max, Pep2max, Kleak, Kleak1]
                param0_global = [(1e-9, 5e-6), (1e-9, 5e-6), (1e-9, 5e-6), (0, 0.02), (0.001, 0.02), (5e-8, 5e-5), (5e-8, 5e-5), (1e-10, 5e-6), (1e-10, 5e-6)] #Diff Evo

                # Local bounds
                LB = [0]*(numParam-5)+[0.001, 1e-10, 1e-10] + [0]*2
                UB = [None]*(numParam)

            elif (SystemType == 'ANDgateBLeak3') or (SystemType == 'ANDgateBLeak2') or (SystemType == 'ANDgateBLeak1'):
                ### Number of Parameters to be optimized
                numParam = 8
                ParamName = ['syn_mRNA1','syn_mRNA2','syn_mRNA3','syn_Pep','deg_Pep','Pep1max','Pep2max','Kleak']
                ParamUnits = ['molL-1min-1', 'molL-1min-1', 'molL-1min-1', 'min-1', 'min-1', 'molL-1', 'molL-1', 'molL-1min-1']
                # Global bounds - In order of [syn_mRNA1, syn_mRNA2, syn_mRNA3, syn_Pep, deg_Pep, Pep1max, Pep2max, Kleak]
                param0_global = [(1e-9, 5e-6), (1e-9, 5e-6), (1e-9, 5e-6), (0, 0.02), (0.001, 0.02), (5e-8, 5e-5), (5e-8, 5e-5), (1e-10, 5e-6)] #Diff Evo

                # Local bounds
                LB = [0]*(numParam-4)+[0.001, 1e-10, 1e-10] + [0]*1
                UB = [None]*(numParam)

            elif (SystemType == 'ANDgate'):
                ### Number of Parameters to be optimized
                numParam = 7
                ParamName = ['syn_mRNA1','syn_mRNA2','syn_mRNA3','syn_Pep','deg_Pep','Pep1max','Pep2max']
                ParamUnits = ['molL-1min-1', 'molL-1min-1', 'molL-1min-1', 'min-1', 'min-1', 'molL-1', 'molL-1']
                # Global bounds - In order of [syn_mRNA1, syn_mRNA2, syn_mRNA3, syn_Pep, deg_Pep, Pep1max, Pep2max, Kleak]
                param0_global = [(1e-9, 5e-6), (1e-9, 5e-6), (1e-9, 5e-6), (0, 0.02), (0.001, 0.02), (5e-8, 5e-5), (5e-8, 5e-5)] #Diff Evo

                # Local bounds
                LB = [0]*(numParam-3)+[0.001] + [0]*2
                UB = [None]*(numParam)

            else:
                print('Error in ANDgate System subType')

        elif (SystemType == 'ANDgateBLeak13KMat'):
            mRNA10 = 0.
            Pep10 = 0.
            mRNA20 = 0.
            Pep20 = 0.
            mRNA30 = 0.
            Pep30 = 0.
            Pep3m0 = data_array[1][0]
            y0 = [mRNA10, Pep10, mRNA20, Pep20, mRNA30, Pep30, Pep3m0]  # Initial condition
            VariableName = ['mRNA1', 'Pep1', 'mRNA2', 'Pep2', 'mRNA3', 'Pep3', 'Pep3m'] # Variables Name
            VarIndex =[VariableName.index('Pep3m')]  #get the index for mRNA and RFP

            numParam = 10
            ParamName = ['syn_mRNA1','syn_mRNA2','syn_mRNA3','syn_Pep','deg_Pep','Pep1max','Pep2max','Kleak','Kleak1', 'Kmature']
            ParamUnits = ['molL-1min-1', 'molL-1min-1', 'molL-1min-1', 'min-1', 'min-1', 'molL-1', 'molL-1', 'molL-1min-1', 'molL-1min-1', 'min-1']
            # Global bounds - In order of [syn_mRNA1, syn_mRNA2, syn_mRNA3, syn_Pep, deg_Pep, Pep1max, Pep2max, Kleak, Kleak1]
            param0_global = [(1e-9, 5e-6), (1e-9, 5e-6), (1e-9, 5e-6), (0, 0.02), (0.001, 0.02), (5e-8, 5e-5), (5e-8, 5e-5), (1e-10, 5e-6), (1e-10, 5e-6), (0.001, 1)] #Diff Evo

            # Local bounds
            LB = [0]*(numParam-6)+[0.001, 1e-10, 1e-10] + [0]*2 + [0.001]
            UB = [None]*(numParam-1) +[1]

        elif (SystemType == 'ORgate'):
            mRNA10 = 0.
            Pep10 = 0.
            mRNA20 = 0.
            Pep20 = 0.
            mRNA30 = 0.
            Pep30 = data_array[1][0]
            y0 = [mRNA10, Pep10, mRNA20, Pep20, mRNA30, Pep30]  # Initial condition
            VariableName = ['mRNA1', 'Pep1', 'mRNA2', 'Pep2', 'mRNA3', 'Pep3'] # Variables Name
            VarIndex =[VariableName.index('Pep3')]  #get the index for mRNA and RFP

            ### Number of Parameters to be optimized
            numParam = 6
            ParamName = ['syn_mRNA1','syn_mRNA2','syn_mRNA3','syn_Pep','deg_Pep','Pepmax']
            ParamUnits = ['molL-1min-1', 'molL-1min-1', 'molL-1min-1', 'min-1', 'min-1', 'molL-1']
            # Global bounds - In order of [syn_mRNA1, syn_mRNA2, syn_mRNA3, syn_Pep, deg_Pep, Pepmax]
            param0_global = [(5e-8, 5e-6), (5e-8, 5e-6), (5e-8, 5e-6), (0, 0.02), (0.001, 0.02), (5e-8, 5e-6)] #Diff Evo

            # Local bounds
            LB = [0]*(numParam-2)+[0.001] + [0]
            UB = [None]*(numParam)

        elif (SystemType == 'ORgateDelay') or (SystemType == 'ORgate_Delay'):
            Inde0 = 1.
            Indi0 = 0.
            mRNA10 = 0.
            Pep10 = 0.
            mRNA20 = 0.
            Pep20 = 0.
            mRNA30 = 0.
            Pep30 = data_array[1][0]
            y0 = [Inde0, Indi0, mRNA10, Pep10, mRNA20, Pep20, mRNA30, Pep30]  # Initial condition
            VariableName = ['Inde', 'Indi', 'mRNA1', 'Pep1', 'mRNA2', 'Pep2', 'mRNA3', 'Pep3'] # Variables Name
            VarIndex =[VariableName.index('Pep3')]  #get the index for mRNA and RFP

            ### Number of Parameters to be optimized
            numParam = 7
            ParamName = ['syn_mRNA1','syn_mRNA2','syn_mRNA3','syn_Pep','deg_Pep','Pepmax', 'Km']
            ParamUnits = ['molL-1min-1', 'molL-1min-1', 'molL-1min-1', 'min-1', 'min-1', 'molL-1', 'dimensionless']
            # Global bounds - In order of [syn_mRNA1, syn_mRNA2, syn_mRNA3, syn_Pep, deg_Pep, Pepmax, Km]
            param0_global = [(5e-8, 5e-6), (5e-8, 5e-6), (5e-8, 5e-6), (0, 0.02), (0.001, 0.02), (5e-8, 5e-6), (0, 70)] #Diff Evo

            # Local bounds
            LB = [0]*(numParam-3)+[0.001] + [0]*2
            UB = [None]*(numParam)

        elif (SystemType == 'ORgateDegradation') or (SystemType == 'ORgate_Degradation'):
            Ind0 = 1.
            mRNA10 = 0.
            Pep10 = 0.
            mRNA20 = 0.
            Pep20 = 0.
            mRNA30 = 0.
            Pep30 = data_array[1][0]
            y0 = [Ind0, mRNA10, Pep10, mRNA20, Pep20, mRNA30, Pep30]  # Initial condition
            VariableName = ['Ind', 'mRNA1', 'Pep1', 'mRNA2', 'Pep2', 'mRNA3', 'Pep3'] # Variables Name
            VarIndex =[VariableName.index('Pep3')]  #get the index for mRNA and RFP

            ### Number of Parameters to be optimized
            numParam = 7
            ParamName = ['syn_mRNA1','syn_mRNA2','syn_mRNA3','syn_Pep','deg_Pep','Pepmax', 'deg_Ind']
            ParamUnits = ['molL-1min-1', 'molL-1min-1', 'molL-1min-1', 'min-1', 'min-1', 'molL-1', 'min-1']
            # Global bounds - In order of [syn_mRNA1, syn_mRNA2, syn_mRNA3, syn_Pep, deg_Pep, Pepmax, Km]
            param0_global = [(5e-8, 5e-6), (5e-8, 5e-6), (5e-8, 5e-6), (0, 0.02), (0.001, 0.02), (5e-8, 5e-6), (0.001, 0.02)] #Diff Evo

            # Local bounds
            LB = [0]*(numParam-3)+[0.001] + [0]*2
            UB = [None]*(numParam)

        elif (SystemType == 'ORgateDelayDegradation') or (SystemType == 'ORgateDegradationDelay') or (SystemType == 'ORgateDelayDegradeResCompete') or (SystemType == 'ORgateDegradeDelayResCompete'):
            Inde0 = 1.
            Indi0 = 0.
            Ind0 = 1.
            mRNA10 = 0.
            Pep10 = 0.
            mRNA20 = 0.
            Pep20 = 0.
            mRNA30 = 0.
            Pep30 = data_array[1][0]
            y0 = [Inde0, Indi0, Ind0, mRNA10, Pep10, mRNA20, Pep20, mRNA30, Pep30]  # Initial condition
            VariableName = ['Inde', 'Indi', 'Ind', 'mRNA1', 'Pep1', 'mRNA2', 'Pep2', 'mRNA3', 'Pep3'] # Variables Name
            VarIndex =[VariableName.index('Pep3')]  #get the index for mRNA and RFP

            ### Number of Parameters to be optimized
            if (SystemType == 'ORgateDelayDegradation') or (SystemType == 'ORgateDegradationDelay'):
                numParam = 8
                ParamName = ['syn_mRNA1','syn_mRNA2','syn_mRNA3','syn_Pep','deg_Pep','Pepmax','Km', 'deg_Ind']
                ParamUnits = ['molL-1min-1', 'molL-1min-1', 'molL-1min-1', 'min-1', 'min-1', 'molL-1', 'dimensionless', 'min-1']
                # Global bounds
                param0_global = [(5e-8, 5e-6), (5e-8, 5e-6), (5e-8, 5e-6), (0, 0.02), (0.001, 0.02), (5e-8, 5e-6),
                                 (0, 70), (0.001, 0.02)] #Diff Evo

                # Local bounds
                LB = [0]*(numParam-4)+[0.001] + [0]*3
                UB = [None]*(numParam)
            else:
                numParam = 9
                ParamName = ['syn_mRNA1','syn_mRNA2','syn_mRNA3','syn_Pep','deg_Pep','Pepmax','Km', 'deg_Ind', 'Ratio']
                ParamUnits = ['molL-1min-1', 'molL-1min-1', 'molL-1min-1', 'min-1', 'min-1', 'molL-1', 'dimensionless', 'min-1', 'dimensionless']
                # Global bounds - In order of [syn_mRNA1, syn_mRNA2, syn_mRNA3, syn_Pep, deg_Pep, Pepmax]
                param0_global = [(5e-8, 5e-6), (5e-8, 5e-6), (5e-8, 5e-6), (0, 0.02), (0.001, 0.02), (5e-8, 5e-6),
                                 (0, 70), (0.001, 0.02), (0, 1)] #Diff Evo

                # Local bounds
                LB = [0]*(numParam-5)+[0.001] + [0]*4
                UB = [None]*(numParam-1) +[1]
        elif (SystemType == 'ORgateDelayDelay') or (SystemType == 'ORgateDelayDelayResCompete'):
            Inde10 = 1.
            Indi10 = 0.
            Inde20 = 1.
            Indi20 = 0.
            mRNA10 = 0.
            Pep10 = 0.
            mRNA20 = 0.
            Pep20 = 0.
            mRNA30 = 0.
            Pep30 = data_array[1][0]
            y0 = [Inde10, Indi10, Inde20, Indi20, mRNA10, Pep10, mRNA20, Pep20, mRNA30, Pep30]  # Initial condition
            VariableName = ['Inde1', 'Indi1', 'Inde2', 'Indi2', 'mRNA1', 'Pep1', 'mRNA2', 'Pep2', 'mRNA3', 'Pep3'] # Variables Name
            VarIndex =[VariableName.index('Pep3')]  #get the index for mRNA and RFP

            if (SystemType == 'ORgateDelayDelay'):
                numParam = 8
                ParamName = ['syn_mRNA1','syn_mRNA2','syn_mRNA3','syn_Pep','deg_Pep','Pepmax','Km1', 'Km2']
                ParamUnits = ['molL-1min-1', 'molL-1min-1', 'molL-1min-1', 'min-1', 'min-1', 'molL-1', 'dimensionless', 'dimensionless']
                # Global bounds
                param0_global = [(5e-8, 5e-6), (5e-8, 5e-6), (5e-8, 5e-6), (0, 0.02), (0.001, 0.02), (5e-8, 5e-6),
                                 (0, 70), (0, 70)] #Diff Evo

                # Local bounds
                LB = [0]*(numParam-4)+[0.001] + [0]*3
                UB = [None]*(numParam)

            else:
                numParam = 9
                ParamName = ['syn_mRNA1','syn_mRNA2','syn_mRNA3','syn_Pep','deg_Pep','Pepmax','Km1', 'Km2', 'Ratio']
                ParamUnits = ['molL-1min-1', 'molL-1min-1', 'molL-1min-1', 'min-1', 'min-1', 'molL-1', 'dimensionless', 'dimensionless', 'dimensionless']
                # Global bounds - In order of [syn_mRNA1, syn_mRNA2, syn_mRNA3, syn_Pep, deg_Pep, Pepmax]
                param0_global = [(5e-8, 5e-6), (5e-8, 5e-6), (5e-8, 5e-6), (0, 0.02), (0.001, 0.02), (5e-8, 5e-6),
                                 (0, 70), (0, 70), (0, 1)] #Diff Evo

                # Local bounds
                LB = [0]*(numParam-5)+[0.001] + [0]*4
                UB = [None]*(numParam-1) +[1]

        else:
            print('Please choose the correct System Type for Inducible Promoter')

        # run Global optimizer
        OptimizerType1 = 'Global'
        result_diffevo = differential_evolution\
            (self.ComputeSSE, param0_global, args=(y0, numState, data_array, Time_interval, SystemType, VarIndex, OptimizerType1))


        # run Local Optimizer (Nelder Mead)
        OptimizerType2 = 'Local'
        param0_local = np.zeros(numParam)
        for i in range(0, numParam):
            param0_local[i] = result_diffevo.x[i]

        result_NM = cNM.constrNM(self.ComputeSSE, param0_local,LB,UB,args=(y0, numState, data_array, Time_interval, SystemType, VarIndex, OptimizerType2),
                                 xtol= 1e-15, full_output=True)

        # Optimized Parameters
        param_optimized = np.zeros(numParam)
        for i in range(0, numParam):
            #param_optimized[i] = result_NM.x[i]
            param_optimized[i] = result_NM['xopt'][i]

        return param_optimized, sse_global, y0, VariableName, ParamName, ParamUnits

    # Plot CSV and Model Data #
    def plotData_Combined(self, SystemType, Variable, y0, raw_data_header, rfp_data, Data_stddev, numState, param):
        ### timespan for Model results (t)
        t_start = rfp_data[0][0]
        t_end = rfp_data[0][-1]
        dt = 1  # minutes
        timestep = int((t_end / dt) + 1)
        t = np.linspace(t_start, t_end, timestep)

        # Time grid == no. of times to report solution
        rfp_data_numrows = np.size(rfp_data, 0)
        rfp_data_numcols = np.size(rfp_data, 1)

        mRNA = np.zeros((timestep, numState), dtype=object)
        Peptide = np.zeros((timestep, numState), dtype=object)
        Inducer = np.zeros((timestep, numState), dtype=object)
        Inducer_e = np.zeros((timestep, numState), dtype=object)
        Inducer_i = np.zeros((timestep, numState), dtype=object)
        Pep2 = np.zeros((numState, 1), dtype=object)

        Variable_mappings = {
                'mRNA': mRNA,
                'Peptide': Peptide,
#                'Inducer': Inducer,
#                'Inducer_e': Inducer_e,
#                'Inducer_i': Inducer_i,
            }

        # initiate empty list to store all the arrays
        VariableMatrix = [None]*len(Variable)

        for i in range(0, len(Variable)):
            VariableMatrix[i] = Variable_mappings[Variable[i]]

        solveODE_Name = '_'.join(('solveODE', SystemType))
        solveODEfun = self.select_function(solveODE_Name) #convert string to function name

        # Integrate the ODE equations
        if numState == 2:
            state = [0,1]
        elif numState == 4:
            state = [[0,0],[0,1],[1,0],[1,1]]
        for i in range(0, numState):  # Iterates through Ind Concs
            ODEsoln = odeint(solveODEfun, y0, t, args=(state[i], param))
            for j in range(0, len(VariableMatrix)):  # Iterates through timesteps
                #for k in range(0, timestep):
                VariableMatrix[j][:,i] = ODEsoln[:,-2+j]
            Pep2[i,0] = ODEsoln[-1,-1] #get the last col last row (which is the final Pep)

        ### Retrive the ODEs in String from the corresponding solveODEfun
        ODEstring = solveODEfun(y0, t[0], state[0], param, 'GetODE')

        print('ODEs in string', ODEstring)

        ### Defining time array from Row 0
        time = rfp_data[0]

        ### Plot RFP Data vs Time ###
        fig = plt.figure(figsize=(5,3.6))
        #The dimensions [left, bottom, width, height]
        ax = fig.add_axes([0.16,0.15,0.8,0.78])
#        plt.rc('axes', linewidth=2)
#        plt.rc('font', size=14) # weight='bold')  # controls default text sizes
        ## CSV Data (time)
        style = ['^','*','>','D','<','s','p','o','d','+','h','x','v','.','H']
        for i in range (1, rfp_data_numrows):
            #plt.plot(time, rfp_data[i], linestyle='None', marker = style[i-1], markersize = 4)
            plt.errorbar(time, rfp_data[i], yerr = Data_stddev[i-1], capsize = 2, linestyle='None', marker = style[i-1], markersize = 4)

        rfp_data_legend = raw_data_header[0:rfp_data_numcols]
        plt.legend(rfp_data_legend, loc='upper left', prop={'size': 16},frameon=False)

        ### Model Data (t)
        # Resets colour cycle for the second plot over same figure
        plt.gca().set_prop_cycle(None)
        # plot model Pep data
        Peptideid = Variable.index('Peptide') #get the index of mRNA
#        print(Peptideid)
#        print(VariableMatrix[Peptideid][:,:])
        #for i in range(0, numInd):
        plt.plot(t, VariableMatrix[Peptideid][:,:], linewidth=2)  # Pep
        #plt.title('Modelled Peptide Concentration vs Time')
        plt.xlabel('Time (min)')
        plt.ylabel('Expression Level (M/OD)')
        # Set Y Axis Ticker to Scientific Style
        plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
        # Figure border
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        axes = plt.gca()
        ymin, ymax = axes.get_ylim()
        axes.set_ylim([0, ymax+0.25*ymax])

        ### Protein Data (t)
        fig = plt.figure(figsize=(5,3.6))
        ax = fig.add_axes([0.16,0.15,0.8,0.78])
#        plt.rc('axes', linewidth=1.5)
#        plt.rc('font', size=14)  # controls default text sizes
        Peptideid = Variable.index('Peptide') #get the index of mRNA
        plt.plot(t, VariableMatrix[Peptideid][:,:], linewidth=2)  # Pep
        plt.legend(rfp_data_legend, loc='upper left', prop={'size': 16},frameon=False)
        #plt.title('Modelled Peptide Concentration vs Time')
        plt.xlabel('Time (min)')
        plt.ylabel('Expression Level (M/OD)')
        # Set Y Axis Ticker to Scientific Style
        plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
        # Figure border
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        axes = plt.gca()
        ymin, ymax = axes.get_ylim()
        axes.set_ylim([0, ymax+0.25*ymax])

        '''
        ### For plotting Bar chart (get the final value under different states) ###
        '''
        exp_final = np.zeros(numState)
        experr_final = np.zeros(numState)
        for i in range(0, numState):
            #rfp_peak[i] = max(rfp_data[i + 1])
            exp_final[i] = rfp_data[i+1,-1]
            experr_final[i] = Data_stddev[i, -1]

        fig = plt.figure(figsize=(5,3.6))
        ax = fig.add_axes([0.16,0.15,0.8,0.78])
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        plt.tight_layout
#        plt.rc('axes', linewidth=1.5)
#        plt.rc('font', size=14) # weight='bold')  # controls default text sizes
        if numState == 2:
            Xlabel = ['0', '1']
        elif numState == 4:
            Xlabel = ['00', '01', '10', '11']
        #x_pos = [i for i, _ in enumerate(Xlabel)]
        xnum = np.arange(numState)
        width1 = 0.4
        print('exp_final = ', exp_final)
        print('Pep2 = ', Pep2)
        plt.bar(xnum, exp_final, yerr=experr_final, align = 'center', ecolor='black', capsize=10, width=width1, color='blue', label = 'Experiment')
        plt.bar(xnum+width1, Pep2[:,0], width=width1, color='red', label = 'Model')
        #ax.legend(prop={'size': 12}, frameon=False)
        plt.xlabel('State')
        plt.ylabel('Expression Level (M/OD)')
        plt.xticks(xnum +width1 / 2, Xlabel)
        plt.legend(prop={'size': 16}, frameon=False, loc='best')
        # Set Y Axis Ticker to Scientific Style
        plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0)) #, fontweight = 'bold')
        axes = plt.gca()
        ymin, ymax = axes.get_ylim()
        axes.set_ylim([0, ymax+0.2*ymax])

        if 'mRNA' in Variable:
            ### mRNA Data (t)
            fig = plt.figure(figsize=(5,3.6))
            ax = fig.add_axes([0.16,0.15,0.8,0.78])
#            plt.rc('axes', linewidth=2)
#            plt.rc('font', size=14)  # controls default text sizes
            mRNAid = Variable.index('mRNA') #get the index of mRNA
            plt.plot(t, VariableMatrix[mRNAid][:,:], linewidth=2)  # mRNA
            #plt.title('Modelled mRNA Concentration vs Time')
            plt.xlabel('Time (min)')
            plt.ylabel('mRNA Concentration (M)')
            plt.legend(rfp_data_legend, loc='upper left', prop={'size': 16},frameon=False)
            # Set Y Axis Ticker to Scientific Style
            plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
            # Figure border
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            axes = plt.gca()
            ymin, ymax = axes.get_ylim()
            axes.set_ylim([0, ymax+0.35*ymax])

#        if 'Inducer' in Variable:
#            ### Inducer Data (t)
#            fig = plt.figure(figsize=(5,3.5))
#            ax = fig.add_axes([0.15,0.15,0.8,0.8])
#            plt.rc('font', size=12)  # controls default text sizes
#            Inducerid = Variable.index('Inducer') #get the index of Inducer
#            plt.plot(t, VariableMatrix[Inducerid][:,:], linewidth=1.5)  # Inducer
#            plt.xlabel('Time (min)')
#            plt.ylabel('Inducer Concentration (M)')
#            plt.legend(rfp_data_legend, loc='upper left', prop={'size': 12},frameon=False)
#            # Set Y Axis Ticker to Scientific Style
#            plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
#            # Figure border
#            ax.spines['top'].set_visible(False)
#            ax.spines['right'].set_visible(False)
#            axes = plt.gca()
#            ymin, ymax = axes.get_ylim()
#            axes.set_ylim([0, ymax+0.35*ymax])
#
#        if 'Inducer_e' in Variable:
#            ### Inducer_e Data (t)
#            fig = plt.figure(figsize=(5,3.5))
#            ax = fig.add_axes([0.15,0.15,0.8,0.8])
#            plt.rc('font', size=12)  # controls default text sizes
#            Inducer_eid = Variable.index('Inducer_e') #get the index of Extracellular Inducer
#            plt.plot(t, VariableMatrix[Inducer_eid][:,:], linewidth=1.5)  # Inducer_e
#            plt.xlabel('Time (min)')
#            plt.ylabel('Inducer$_{ex}$ Concentration (M)')
#            plt.legend(rfp_data_legend, loc='upper left', prop={'size': 12},frameon=False)
#            # Set Y Axis Ticker to Scientific Style
#            plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
#            # Figure border
#            ax.spines['top'].set_visible(False)
#            ax.spines['right'].set_visible(False)
#            axes = plt.gca()
#            ymin, ymax = axes.get_ylim()
#            axes.set_ylim([0, ymax+0.35*ymax])
#
#        if 'Inducer_i' in Variable:
#            ### Inducer_i Data (t)
#            fig = plt.figure(figsize=(5,3.5))
#            ax = fig.add_axes([0.15,0.15,0.8,0.8])
#            plt.rc('font', size=12)  # controls default text sizes
#            Inducer_iid = Variable.index('Inducer_i') #get the index of Intracellular Inducer
#            plt.plot(t, VariableMatrix[Inducer_iid][:,:], linewidth=1.5)  # Inducer_i
#            plt.xlabel('Time (min)')
#            plt.ylabel('Inducer$_{in}$ Concentration (M)')
#            plt.legend(rfp_data_legend, loc='upper left', prop={'size': 12},frameon=False)
#            # Set Y Axis Ticker to Scientific Style
#            plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
#            # Figure border
#            ax.spines['top'].set_visible(False)
#            ax.spines['right'].set_visible(False)
#            axes = plt.gca()
#            ymin, ymax = axes.get_ylim()
#            axes.set_ylim([0, ymax+0.35*ymax])

        return t, VariableMatrix, rfp_data_legend, ODEstring


    def Run_LogicGatesPlot(self, SystemType, y0, data_header, data_array, Data_stddev, numState, param_optimized):
        if 'NOTgate' in SystemType:
            VariablePlot = ['mRNA', 'Peptide'] # Variables Name

        elif 'ANDgate' in SystemType:
            VariablePlot = ['mRNA', 'Peptide'] # Variables Name

        elif 'ORgate' in SystemType:
            VariablePlot = ['mRNA', 'Peptide'] # Variables Name

        else:
            print('Error in Plotting module')

        ### Calculate and plot Model results (param_optimized) ###
        VariableMatrix = self.plotData_Combined(SystemType, VariablePlot, y0, data_header, data_array, Data_stddev, numState, param_optimized)

        return VariableMatrix

    def __del__(self):
      class_name = self.__class__.__name__
      print(class_name, "destroyed")
