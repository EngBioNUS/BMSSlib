# -*- coding: utf-8 -*-
"""
Created on Mon Jun 25 09:50:57 2018

@author: JingWui
"""

import numpy as np
import constrNMPy as cNM
import matplotlib.pyplot as plt

from scipy.integrate import odeint
from scipy.optimize import differential_evolution

class ConstitutivePromoterLibrary:

    ### ODE Model for Constitutive Promoter ###
    def solveODE_ConstDouble(y, t, param, TotalDataSet, Operation = 'Solve'): #must include the comma at the end
        # Dependent variables
        mRNA = y[0] # Col 0 of ODESoln
        Pep = y[1]  # Col 1 of ODESoln

        # Parameters
        syn_mRNA = param[0]
        deg_mRNA = 0.1386   # Fixed as ln2/5(mins) = 0.1386
        syn_Pep = param[1]
        deg_Pep = param[2]

        # Differential equations
        dmRNA_dt = '(syn_mRNA)-(deg_mRNA * mRNA)'
        dPep_dt = '(syn_Pep*mRNA)-(deg_Pep*Pep)'

        # Return differential equations solution
        if Operation == 'Solve':
            return [eval(dmRNA_dt), eval(dPep_dt)]
        elif Operation == 'GetODE':
            return [dmRNA_dt, dPep_dt]
        else:
            print('Error: Please Enter the correct Operation for solveODE function')
            
    def solveODE_ConstDoubleKMat(y, t, param, TotalDataSet, Operation = 'Solve'): #must include the comma at the end
        # Dependent variables
        mRNA = y[0] # Col 0 of ODESoln
        Pep = y[1]  # Col 1 of ODESoln
        Pepm = y[2]

        # Parameters
        syn_mRNA = param[0]
        deg_mRNA = 0.1386   # Fixed as ln2/5(mins) = 0.1386
        syn_Pep = param[1]
        deg_Pep = param[2]
        Kmature = param[3]
        #Kmature = 0.0316 #log(2)/21.9 min

        # Differential equations
        dmRNA_dt = '(syn_mRNA)-(deg_mRNA * mRNA)'
        dPep_dt = '(syn_Pep*mRNA)-(Kmature*Pep)'
        dPepm_dt = '(Kmature*Pep)-(deg_Pep*Pepm)'

        # Return differential equations solution
        if Operation == 'Solve':
            return [eval(dmRNA_dt), eval(dPep_dt), eval(dPepm_dt)]
        elif Operation == 'GetODE':
            return [dmRNA_dt, dPep_dt, dPepm_dt]
        else:
            print('Error: Please Enter the correct Operation for solveODE function')

    ### Single-ODE Model for Constitutive Promoter ###
    def solveODE_ConstSingle(y, t, param, TotalDataSet, Operation = 'Solve'):
        # Dependent variables
        Pep = y[0]  # Col 1 of ODESoln

        # Parameters
        syn_Pep = param[0]
        deg_Pep = param[1]

        # Differential equations
        dPep_dt = '(syn_Pep)-(deg_Pep*Pep)'

        # Return differential equations solution
        # Return differential equations solution
        if Operation == 'Solve':
            return [eval(dPep_dt)]
        elif Operation == 'GetODE':
            return [dPep_dt]
        else:
            print('Error: Please Enter the correct Operation for solveODE function')
            
    def solveODE_ConstSingleKMat(y, t, param, TotalDataSet, Operation = 'Solve'):
        # Dependent variables
        Pep = y[0]  # Col 1 of ODESoln
        Pepm = y[1]

        # Parameters
        syn_Pep = param[0]
        deg_Pep = param[1]
        Kmature = param[2]
        #Kmature = 0.0316 # log(2)/21.9 (mRFP1 at 37'C)

        # Differential equations
        dPep_dt = '(syn_Pep)-(Kmature*Pep)'
        dPepm_dt = '(Kmature*Pep)-(deg_Pep*Pepm)'

        # Return differential equations solution
        # Return differential equations solution
        if Operation == 'Solve':
            return [eval(dPep_dt), eval(dPepm_dt)]
        elif Operation == 'GetODE':
            return [dPep_dt, dPepm_dt]
        else:
            print('Error: Please Enter the correct Operation for solveODE function')

    ### Multi-ODE Model (different promoters same RBS) for Constitutive Promoter ###
    def solveODE_MultiDoubleFixRBS(y, t, param, TotalDataSet, Operation = 'Solve'): #must include the comma at the end
        # Dependent variables
        mRNA1 = y[0] # Col 0 of ODESoln
        Pep1 = y[1]  # Col 1 of ODESoln
        mRNA2 = y[2] # Col 0 of ODESoln
        Pep2 = y[3]

        # Parameters
        deg_mRNA = 0.1386   # Fixed as ln2/5(mins) = 0.1386
        syn_Pep = param[0]
        deg_Pep = param[1]
        syn_mRNA1 = param[2]
        syn_mRNA2 = param[3]

        # Differential equations
        dmRNA1_dt = '(syn_mRNA1)-(deg_mRNA*mRNA1)'
        dPep1_dt = '(syn_Pep*mRNA1)-(deg_Pep*Pep1)'
        dmRNA2_dt = '(syn_mRNA2)-(deg_mRNA*mRNA2)'
        dPep2_dt = '(syn_Pep*mRNA2)-(deg_Pep*Pep2)'

        dys = [eval(dmRNA1_dt), eval(dPep1_dt), eval(dmRNA2_dt), eval(dPep2_dt)]
        dy = [dmRNA1_dt, dPep1_dt, dmRNA2_dt, dPep2_dt]

        if TotalDataSet >= 3:
            # Dependent variables
            mRNA3 = y[4] # Col 0 of ODESoln
            Pep3 = y[5]  # Col 1 of ODESoln

            syn_mRNA3 = param[4]

            # Differential equations
            dmRNA3_dt = '(syn_mRNA3)-(deg_mRNA*mRNA3)'
            dPep3_dt = '(syn_Pep*mRNA3)-(deg_Pep*Pep3)'

            dys = dys + [eval(dmRNA3_dt), eval(dPep3_dt)]
            dy = dy + [dmRNA3_dt, dPep3_dt]

            if TotalDataSet >= 4:
                # Dependent variables
                mRNA4 = y[6]
                Pep4 = y[7]

                syn_mRNA4 = param[5]

                # Differential equations
                dmRNA4_dt = '(syn_mRNA4)-(deg_mRNA*mRNA4)'
                dPep4_dt = '(syn_Pep*mRNA4)-(deg_Pep*Pep4)'

                dys = dys + [eval(dmRNA4_dt), eval(dPep4_dt)]
                dy = dy + [dmRNA4_dt, dPep4_dt]

                if TotalDataSet >= 5:
                    # Dependent variables
                    mRNA5 = y[8]
                    Pep5 = y[9]

                    syn_mRNA5 = param[6]

                    # Differential equations
                    dmRNA5_dt = '(syn_mRNA5)-(deg_mRNA*mRNA5)'
                    dPep5_dt = '(syn_Pep*mRNA5)-(deg_Pep*Pep5)'

                    dys = dys + [eval(dmRNA5_dt), eval(dPep5_dt)]
                    dy = dy + [dmRNA5_dt, dPep5_dt]
                    
                    if TotalDataSet >= 6:
                        # Dependent variables
                        mRNA6 = y[10]
                        Pep6 = y[11]
    
                        syn_mRNA6 = param[7]
    
                        # Differential equations
                        dmRNA6_dt = '(syn_mRNA6)-(deg_mRNA*mRNA6)'
                        dPep6_dt = '(syn_Pep*mRNA6)-(deg_Pep*Pep6)'
    
                        dys = dys + [eval(dmRNA6_dt), eval(dPep6_dt)]
                        dy = dy + [dmRNA6_dt, dPep6_dt]
                    
        if Operation == 'Solve': 
            return dys            
        elif Operation == 'GetODE':
            return dy
        else:
            print('Error: Please Enter the correct Operation for solveODE function')
            
    def solveODE_MultiDoubleFixRBSKMat(y, t, param, TotalDataSet, Operation = 'Solve'): #must include the comma at the end
        # Dependent variables
        mRNA1 = y[0] # Col 0 of ODESoln
        Pep1 = y[1]  # Col 1 of ODESoln
        Pepm1 = y[2]
        mRNA2 = y[3] # Col 0 of ODESoln
        Pep2 = y[4]
        Pepm2 = y[5]

        # Parameters
        deg_mRNA = 0.1386   # Fixed as ln2/5(mins) = 0.1386
        syn_Pep = param[0]
        deg_Pep = param[1]
        syn_mRNA1 = param[2]
        syn_mRNA2 = param[3]
        Kmature = param[4]

        # Differential equations
        dmRNA1_dt = '(syn_mRNA1)-(deg_mRNA*mRNA1)'
        dPep1_dt = '(syn_Pep*mRNA1)-(Kmature*Pep1)'
        dPepm1_dt = '(Kmature*Pep1)-(deg_Pep*Pepm1)'
        dmRNA2_dt = '(syn_mRNA2)-(deg_mRNA*mRNA2)'
        dPep2_dt = '(syn_Pep*mRNA2)-(Kmature*Pep2)'
        dPepm2_dt = '(Kmature*Pep2)-(deg_Pep*Pepm2)'

        dys = [eval(dmRNA1_dt), eval(dPep1_dt), eval(dPepm1_dt), eval(dmRNA2_dt), eval(dPep2_dt), eval(dPepm2_dt)]
        dy = [dmRNA1_dt, dPep1_dt, dPepm1_dt, dmRNA2_dt, dPep2_dt, dPepm2_dt]

        if TotalDataSet >= 3:
            # Dependent variables
            mRNA3 = y[6] # Col 0 of ODESoln
            Pep3 = y[7]  # Col 1 of ODESoln
            Pepm3 = y[8]

            syn_mRNA3 = param[5]

            # Differential equations
            dmRNA3_dt = '(syn_mRNA3)-(deg_mRNA*mRNA3)'
            dPep3_dt = '(syn_Pep*mRNA3)-(Kmature*Pep3)'
            dPepm3_dt = '(Kmature*Pep3)-(deg_Pep*Pepm3)'

            dys = dys + [eval(dmRNA3_dt), eval(dPep3_dt), eval(dPepm3_dt)]
            dy = dy + [dmRNA3_dt, dPep3_dt, dPepm3_dt]

            if TotalDataSet >= 4:
                # Dependent variables
                mRNA4 = y[9]
                Pep4 = y[10]
                Pepm4 = y[11]

                syn_mRNA4 = param[6]

                # Differential equations
                dmRNA4_dt = '(syn_mRNA4)-(deg_mRNA*mRNA4)'
                dPep4_dt = '(syn_Pep*mRNA4)-(Kmature*Pep4)'
                dPepm4_dt = '(Kmature*Pep4)-(deg_Pep*Pepm4)'

                dys = dys + [eval(dmRNA4_dt), eval(dPep4_dt), eval(dPepm4_dt)]
                dy = dy + [dmRNA4_dt, dPep4_dt, dPepm4_dt]

                if TotalDataSet >= 5:
                    # Dependent variables
                    mRNA5 = y[12]
                    Pep5 = y[13]
                    Pepm5 = y[14]

                    syn_mRNA5 = param[7]

                    # Differential equations
                    dmRNA5_dt = '(syn_mRNA5)-(deg_mRNA*mRNA5)'
                    dPep5_dt = '(syn_Pep*mRNA5)-(Kmature*Pep5)'
                    dPepm5_dt = '(Kmature*Pep5)-(deg_Pep*Pepm5)'

                    dys = dys + [eval(dmRNA5_dt), eval(dPep5_dt), eval(dPepm5_dt)]
                    dy = dy + [dmRNA5_dt, dPep5_dt, dPepm5_dt]
                    
                    if TotalDataSet >= 6:
                        # Dependent variables
                        mRNA6 = y[15]
                        Pep6 = y[16]
                        Pepm6 = y[17]
    
                        syn_mRNA6 = param[8]
    
                        # Differential equations
                        dmRNA6_dt = '(syn_mRNA6)-(deg_mRNA*mRNA6)'
                        dPep6_dt = '(syn_Pep*mRNA6)-(Kmature*Pep6)'
                        dPepm6_dt = '(Kmature*Pep6)-(deg_Pep*Pepm6)'
    
                        dys = dys + [eval(dmRNA6_dt), eval(dPep6_dt), eval(dPepm6_dt)]
                        dy = dy + [dmRNA6_dt, dPep6_dt, dPepm6_dt]
                    
        if Operation == 'Solve': 
            return dys            
        elif Operation == 'GetODE':
            return dy
        else:
            print('Error: Please Enter the correct Operation for solveODE function')
            
    def solveODE_MultiSingleFixRBS(y, t, param, TotalDataSet, Operation = 'Solve'): #must include the comma at the end
        # Dependent variables
        Pep1 = y[0]  # Col 1 of ODESoln
        Pep2 = y[1]

        # Parameters
        deg_Pep = param[0]
        syn_Pep1 = param[1]
        syn_Pep2 = param[2]
        
        # Differential equations
        dPep1_dt = '(syn_Pep1)-(deg_Pep*Pep1)'
        dPep2_dt = '(syn_Pep2)-(deg_Pep*Pep2)'

        dys = [eval(dPep1_dt), eval(dPep2_dt)]
        dy = [dPep1_dt, dPep2_dt]

        if TotalDataSet >= 3:
            # Dependent variables
            Pep3 = y[2]  # Col 1 of ODESoln

            syn_Pep3 = param[3]

            # Differential equations
            dPep3_dt = '(syn_Pep3)-(deg_Pep*Pep3)'

            dys = dys + [eval(dPep3_dt)]
            dy = dy + [dPep3_dt]

            if TotalDataSet >= 4:
                # Dependent variables
                Pep4 = y[3]

                syn_Pep4 = param[4]

                # Differential equations
                dPep4_dt = '(syn_Pep4)-(deg_Pep*Pep4)'

                dys = dys + [eval(dPep4_dt)]
                dy = dy + [dPep4_dt]

                if TotalDataSet >= 5:
                    # Dependent variables
                    Pep5 = y[4]

                    syn_Pep5 = param[5]

                    # Differential equations
                    dPep5_dt = '(syn_Pep5)-(deg_Pep*Pep5)'

                    dys = dys + [eval(dPep5_dt)]
                    dy = dy + [dPep5_dt]
                    
                    if TotalDataSet >= 6:
                        # Dependent variables
                        Pep6 = y[5]
    
                        syn_Pep6 = param[6]
    
                        # Differential equations
                        dPep6_dt = '(syn_Pep6)-(deg_Pep*Pep6)'
    
                        dys = dys + [eval(dPep6_dt)]
                        dy = dy + [dPep6_dt]
                    
        if Operation == 'Solve': 
            return dys            
        elif Operation == 'GetODE':
            return dy
        else:
            print('Error: Please Enter the correct Operation for solveODE function')
            
    def solveODE_MultiSingleFixRBSKMat(y, t, param, TotalDataSet, Operation = 'Solve'): #must include the comma at the end
        # Dependent variables
        Pep1 = y[0]  # Col 1 of ODESoln
        Pepm1 = y[1]
        Pep2 = y[2]
        Pepm2 = y[3]

        # Parameters
        deg_Pep = param[0]
        syn_Pep1 = param[1]
        syn_Pep2 = param[2]
        Kmature = param[3]
        
        # Differential equations
        dPep1_dt = '(syn_Pep1)-(Kmature*Pep1)'
        dPepm1_dt = '(Kmature*Pep1)-(deg_Pep*Pepm1)'
        dPep2_dt = '(syn_Pep2)-(Kmature*Pep2)'
        dPepm2_dt = '(Kmature*Pep2)-(deg_Pep*Pepm2)'

        dys = [eval(dPep1_dt), eval(dPepm1_dt), eval(dPep2_dt), eval(dPepm2_dt)]
        dy = [dPep1_dt, dPepm1_dt, dPep2_dt, dPepm2_dt]

        if TotalDataSet >= 3:
            # Dependent variables
            Pep3 = y[4]  # Col 1 of ODESoln
            Pepm3 = y[5]

            syn_Pep3 = param[4]

            # Differential equations
            dPep3_dt = '(syn_Pep3)-(Kmature*Pep3)'
            dPepm3_dt = '(Kmature*Pep3)-(deg_Pep*Pepm3)'

            dys = dys + [eval(dPep3_dt), eval(dPepm3_dt)]
            dy = dy + [dPep3_dt, dPepm3_dt]

            if TotalDataSet >= 4:
                # Dependent variables
                Pep4 = y[6]
                Pepm4 = y[7]

                syn_Pep4 = param[5]

                # Differential equations
                dPep4_dt = '(syn_Pep4)-(Kmature*Pep4)'
                dPepm4_dt = '(Kmature*Pep4)-(deg_Pep*Pepm4)'

                dys = dys + [eval(dPep4_dt), eval(dPepm4_dt)]
                dy = dy + [dPep4_dt, dPepm4_dt]

                if TotalDataSet >= 5:
                    # Dependent variables
                    Pep5 = y[8]
                    Pepm5 = y[9]

                    syn_Pep5 = param[6]

                    # Differential equations
                    dPep5_dt = '(syn_Pep5)-(Kmature*Pep5)'
                    dPepm5_dt = '(Kmature*Pep5)-(deg_Pep*Pepm5)'

                    dys = dys + [eval(dPep5_dt), eval(dPepm5_dt)]
                    dy = dy + [dPep5_dt, dPepm5_dt]
                    
                    if TotalDataSet >= 6:
                        # Dependent variables
                        Pep6 = y[10]
                        Pepm6 = y[11]
    
                        syn_Pep6 = param[7]
    
                        # Differential equations
                        dPep6_dt = '(syn_Pep6)-(Kmature*Pep6)'
                        dPepm6_dt = '(Kmature*Pep6)-(deg_Pep*Pepm6)'
    
                        dys = dys + [eval(dPep6_dt), eval(dPepm6_dt)]
                        dy = dy + [dPep6_dt, dPepm6_dt]
                    
        if Operation == 'Solve': 
            return dys            
        elif Operation == 'GetODE':
            return dy
        else:
            print('Error: Please Enter the correct Operation for solveODE function')

    ### Multi-ODE Model (different RBSs same promoter) for Constitutive Promoter ###
    def solveODE_MultiDoubleFixPromoter(y, t, param, TotalDataSet, Operation = 'Solve'): #must include the comma at the end
        # Dependent variables
        mRNA1 = y[0]
        Pep1 = y[1]
        mRNA2 = y[2]
        Pep2 = y[3]

        # Parameters
        deg_mRNA = 0.1386   # Fixed as ln2/5(mins) = 0.1386
        syn_mRNA = param[0]
        deg_Pep = param[1]
        syn_Pep1 = param[2]
        syn_Pep2 = param[3]

        # Differential equations
        dmRNA1_dt = '(syn_mRNA)-(deg_mRNA*mRNA1)'
        dPep1_dt = '(syn_Pep1*mRNA1)-(deg_Pep*Pep1)'
        dmRNA2_dt = '(syn_mRNA)-(deg_mRNA*mRNA2)'
        dPep2_dt = '(syn_Pep2*mRNA2)-(deg_Pep*Pep2)'

        dys = [eval(dmRNA1_dt), eval(dPep1_dt), eval(dmRNA2_dt), eval(dPep2_dt)]
        dy = [dmRNA1_dt, dPep1_dt, dmRNA2_dt, dPep2_dt]

        if TotalDataSet >= 3:
            # Dependent variables
            mRNA3 = y[4]
            Pep3 = y[5]

            syn_Pep3 = param[4]

            # Differential equations
            dmRNA3_dt = '(syn_mRNA)-(deg_mRNA*mRNA3)'
            dPep3_dt = '(syn_Pep3*mRNA3)-(deg_Pep*Pep3)'
            
            dys = dys + [eval(dmRNA3_dt), eval(dPep3_dt)]
            dy = dy + [dmRNA3_dt, dPep3_dt]

            if TotalDataSet >= 4:
                # Dependent variables
                mRNA4 = y[6]
                Pep4 = y[7]

                syn_Pep4 = param[5]

                # Differential equations
                dmRNA4_dt = '(syn_mRNA)-(deg_mRNA*mRNA4)'
                dPep4_dt = '(syn_Pep4*mRNA4)-(deg_Pep*Pep4)'

                dys = dys + [eval(dmRNA4_dt), eval(dPep4_dt)]
                dy = dy + [dmRNA4_dt, dPep4_dt]

                if TotalDataSet >= 5:
                    # Dependent variables
                    mRNA5 = y[8]
                    Pep5 = y[9]

                    syn_Pep5 = param[6]

                    # Differential equations
                    dmRNA5_dt = '(syn_mRNA)-(deg_mRNA*mRNA5)'
                    dPep5_dt = '(syn_Pep5*mRNA5)-(deg_Pep*Pep5)'

                    dys = dys + [eval(dmRNA5_dt), eval(dPep5_dt)]
                    dy = dy + [dmRNA5_dt, dPep5_dt]
                    
                    if TotalDataSet >= 6:
                        # Dependent variables
                        mRNA6 = y[10]
                        Pep6 = y[11]
    
                        syn_Pep6 = param[7]
    
                        # Differential equations
                        dmRNA6_dt = '(syn_mRNA)-(deg_mRNA*mRNA6)'
                        dPep6_dt = '(syn_Pep6*mRNA6)-(deg_Pep*Pep6)'
    
                        dys = dys + [eval(dmRNA6_dt), eval(dPep6_dt)]
                        dy = dy + [dmRNA6_dt, dPep6_dt]

        if Operation == 'Solve': 
            return dys            
        elif Operation == 'GetODE':
            return dy
        else:
            print('Error: Please Enter the correct Operation for solveODE function')
            
    def solveODE_MultiDoubleFixPromoterKMat(y, t, param, TotalDataSet, Operation = 'Solve'): #must include the comma at the end
        # Dependent variables
        mRNA1 = y[0]
        Pep1 = y[1]
        Pepm1 = y[2]
        mRNA2 = y[3]
        Pep2 = y[4]
        Pepm2 = y[5]

        # Parameters
        deg_mRNA = 0.1386   # Fixed as ln2/5(mins) = 0.1386
        syn_mRNA = param[0]
        deg_Pep = param[1]
        syn_Pep1 = param[2]
        syn_Pep2 = param[3]
        Kmature = param[4]

        # Differential equations
        dmRNA1_dt = '(syn_mRNA)-(deg_mRNA*mRNA1)'
        dPep1_dt = '(syn_Pep1*mRNA1)-(Kmature*Pep1)'
        dPepm1_dt = '(Kmature*Pep1)-(deg_Pep*Pepm1)'
        dmRNA2_dt = '(syn_mRNA)-(deg_mRNA*mRNA2)'
        dPep2_dt = '(syn_Pep2*mRNA2)-(Kmature*Pep2)'
        dPepm2_dt = '(Kmature*Pep2)-(deg_Pep*Pepm2)'

        dys = [eval(dmRNA1_dt), eval(dPep1_dt), eval(dPepm1_dt), eval(dmRNA2_dt), eval(dPep2_dt), eval(dPepm2_dt)]
        dy = [dmRNA1_dt, dPep1_dt, dPepm1_dt, dmRNA2_dt, dPep2_dt, dPepm2_dt]

        if TotalDataSet >= 3:
            # Dependent variables
            mRNA3 = y[6]
            Pep3 = y[7]
            Pepm3 = y[8]

            syn_Pep3 = param[5]

            # Differential equations
            dmRNA3_dt = '(syn_mRNA)-(deg_mRNA*mRNA3)'
            dPep3_dt = '(syn_Pep3*mRNA3)-(Kmature*Pep3)'
            dPepm3_dt = '(Kmature*Pep3)-(deg_Pep*Pepm3)'
            
            dys = dys + [eval(dmRNA3_dt), eval(dPep3_dt), eval(dPepm3_dt)]
            dy = dy + [dmRNA3_dt, dPep3_dt, dPepm3_dt]

            if TotalDataSet >= 4:
                # Dependent variables
                mRNA4 = y[9]
                Pep4 = y[10]
                Pepm4 = y[11]

                syn_Pep4 = param[6]

                # Differential equations
                dmRNA4_dt = '(syn_mRNA)-(deg_mRNA*mRNA4)'
                dPep4_dt = '(syn_Pep4*mRNA4)-(Kmature*Pep4)'
                dPepm4_dt = '(Kmature*Pep4)-(deg_Pep*Pepm4)'

                dys = dys + [eval(dmRNA4_dt), eval(dPep4_dt), eval(dPepm4_dt)]
                dy = dy + [dmRNA4_dt, dPep4_dt, dPepm4_dt]

                if TotalDataSet >= 5:
                    # Dependent variables
                    mRNA5 = y[12]
                    Pep5 = y[13]
                    Pepm5 = y[14]

                    syn_Pep5 = param[7]

                    # Differential equations
                    dmRNA5_dt = '(syn_mRNA)-(deg_mRNA*mRNA5)'
                    dPep5_dt = '(syn_Pep5*mRNA5)-(Kmature*Pep5)'
                    dPepm5_dt = '(Kmature*Pep5)-(deg_Pep*Pepm5)'

                    dys = dys + [eval(dmRNA5_dt), eval(dPep5_dt), eval(dPepm5_dt)]
                    dy = dy + [dmRNA5_dt, dPep5_dt, dPepm5_dt]
                    
                    if TotalDataSet >= 6:
                        # Dependent variables
                        mRNA6 = y[15]
                        Pep6 = y[16]
                        Pepm6 = y[17]
    
                        syn_Pep6 = param[8]
    
                        # Differential equations
                        dmRNA6_dt = '(syn_mRNA)-(deg_mRNA*mRNA6)'
                        dPep6_dt = '(syn_Pep6*mRNA6)-(Kmature*Pep6)'
                        dPepm6_dt = '(Kmature*Pep6)-(deg_Pep*Pepm6)'
    
                        dys = dys + [eval(dmRNA6_dt), eval(dPep6_dt), eval(dPepm6_dt)]
                        dy = dy + [dmRNA6_dt, dPep6_dt, dPepm6_dt]

        if Operation == 'Solve': 
            return dys            
        elif Operation == 'GetODE':
            return dy
        else:
            print('Error: Please Enter the correct Operation for solveODE function')
            
    

    function_mappings = {
        'solveODE_ConstSingle': solveODE_ConstSingle,
        'solveODE_ConstSingleKMat': solveODE_ConstSingleKMat,
        'solveODE_ConstDouble': solveODE_ConstDouble,
        'solveODE_ConstDoubleKMat': solveODE_ConstDoubleKMat,
        'solveODE_MultiDoubleFixRBS': solveODE_MultiDoubleFixRBS,
        'solveODE_MultiDoubleFixRBSKMat': solveODE_MultiDoubleFixRBSKMat,
        'solveODE_MultiSingleFixRBS': solveODE_MultiSingleFixRBS,
        'solveODE_MultiSingleFixRBSKMat': solveODE_MultiSingleFixRBSKMat, 
        'solveODE_MultiDoubleFixPromoter': solveODE_MultiDoubleFixPromoter,
        'solveODE_MultiDoubleFixPromoterKMat': solveODE_MultiDoubleFixPromoterKMat,
    }

    def select_function(self, SystemTypeString):
        ### Convert SystemType from string to function name
        try:
            return self.function_mappings[SystemTypeString]
        except KeyError:
            print('Invalid function, try again.')

    def ComputeSSE(self, param, y0, rfp_data, time_int, SystemType, VarIndex, OptimizerType, TotalDataSet):
        ### Calculate SSE - To be minimized by Optimizer ####

        global sse_global

        # Time grid == no. of times to report solution
        rfp_data_numrows = np.size(rfp_data, 0)
        rfp_data_numcols = np.size(rfp_data, 1)

        t_start = rfp_data[0][0] # First value of Time column
        t_end = rfp_data[0][-1] # Last value of Time column
        dt = 10  # minutes
        timestep = int((t_end / dt) + 1)
        t = np.linspace(t_start, t_end, timestep)

        # Initialize mRNA and Pep nested list
        #mRNA = np.zeros((timestep, numInd), dtype=object)  # timestep (rows) x numInd (cols)
        Pep = np.zeros((timestep, TotalDataSet), dtype=object)   # timestep (rows) x numInd (cols)

        solveODE_Name = '_'.join(('solveODE', SystemType))
        solveODEfun = self.select_function(solveODE_Name) #convert string to function name

        if OptimizerType == 'Global':
            # Integrate the ODE equations
            ODEsoln = odeint(solveODEfun, y0, t, args=(param, TotalDataSet)) #must have more than one additional arg
            #mRNA[j][i] = ODEsoln[j][VarIndex[0]] # mRNA array runs downwards with time
            for i in range(0, len(VarIndex)):
                for j in range(0, timestep):
                    Pep[j][i] = ODEsoln[j][VarIndex[i]]  # Pep array runs downwards with time


        elif OptimizerType == 'Local':
            # Integrate the ODE equations
            ODEsoln = odeint(solveODEfun, y0, t, args=(param, TotalDataSet))
            for i in range(0, len(VarIndex)):
                for j in range(0, timestep):
                #mRNA[j][i] = ODEsoln[j][VarIndex[0]] # mRNA array runs downwards with time
                    Pep[j][i] = ODEsoln[j][VarIndex[i]]  # Pep array runs downwards with time
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

    def Run_ConstitutiveSystem(self, SystemType, data_header, data_array, TotalDataSet):

        Time_interval = data_array[0][1] - data_array[0][0]
        
        ### ODE Input (Inducer Degradation) ###
        # Initial conditions for (mRNA, Pep, Ind) at time = 0
        if (SystemType == 'ConstDouble'):
            mRNA0 = 0.
            Pep0 = data_array[1][0]
            y0 = [mRNA0, Pep0]  # Initial condition
            VariableName = ['mRNA', 'Pep'] # Variables Name
            VarIndex =[VariableName.index('Pep')]  #get the index for mRNA and RFP

            ### Number of Parameters to be optimized
            numParam = 3
            
            ParamName = ['syn_mRNA','syn_Pep','deg_Pep']
            ParamUnits = ['molL-1min-1', 'min-1', 'min-1']

            # Global bounds - In order of (a_mRNA, a_Pep, y_Pep)
            param0_global = [(5e-8, 5e-7), (0, 0.02), (0.001, 0.02)] #Diff Evo

            # Local bounds
            LB = [1e-10, 1e-6, 0.001]
            UB = [5e-7, 0.02, 0.02]
            #LB = [0]*(numParam-1) + [0.001]
            #UB = [None]*(numParam)
            
        elif (SystemType == 'ConstDoubleKMat'):
            mRNA0 = 0.
            Pep0 = 0.
            Pepm0 = data_array[1][0]
            y0 = [mRNA0, Pep0, Pepm0]  # Initial condition
            VariableName = ['mRNA', 'Pep', 'Pepm'] # Variables Name
            VarIndex =[VariableName.index('Pepm')]  #get the index for mRNA and RFP

            ### Number of Parameters to be optimized
            numParam = 4
            
            ParamName = ['syn_mRNA','syn_Pep','deg_Pep', 'Kmature']
            ParamUnits = ['molL-1min-1', 'min-1', 'min-1', 'min-1']

            # Global bounds - In order of (a_mRNA, a_Pep, y_Pep)
            param0_global = [(5e-8, 5e-7), (0, 0.02), (0.001, 0.02), (0.001, 1)] #Diff Evo

            # Local bounds
            LB = [1e-10, 1e-6, 0.001, 0.002]
            UB = [5e-7, 0.02, 0.02, 1]


        elif (SystemType == 'ConstSingle'):
            Pep0 = data_array[1][0]
            y0 = [Pep0]  # Initial condition
            VariableName = ['Pep'] # Variables Name
            VarIndex =[VariableName.index('Pep')]  #get the index for RFP

            ### Number of Parameters to be optimized
            numParam = 2 # Fixed for Constant Inducer Model
 
            ParamName = ['syn_Pep','deg_Pep']
            ParamUnits = ['molL-1min-1', 'min-1']

            # Global bounds - In order of (a_Pep, y_Pep)
            param0_global = [(5e-8, 5e-7), (0.001, 0.02)] #Diff Evo

            # Local bounds
            LB = [1e-10, 0.001]
            UB = [5e-7, 0.02]
            #LB = [0]*(numParam-1) + [0.001]
            #UB = [None]*(numParam)
            
        elif (SystemType == 'ConstSingleKMat'):
            Pep0 = 0.
            Pepm0 = data_array[1][0]
            y0 = [Pep0, Pepm0]  # Initial condition
            VariableName = ['Pep', 'Pepm'] # Variables Name
            VarIndex =[VariableName.index('Pepm')]  #get the index for RFP

            ### Number of Parameters to be optimized
            numParam = 3 # Fixed for Constant Inducer Model
 
            ParamName = ['syn_Pep','deg_Pep', 'Kmature']
            ParamUnits = ['molL-1min-1', 'min-1', 'min-1']

            # Global bounds - In order of (a_Pep, y_Pep)
            param0_global = [(5e-8, 5e-7), (0.001, 0.02), (0.001, 1)] #Diff Evo

            # Local bounds
            LB = [1e-10, 0.001, 0.002]
            UB = [5e-7, 0.02, 1]

        elif (SystemType == 'MultiDoubleFixRBS'):

            mRNA10 = 0.
            Pep10 = data_array[1][0]
            mRNA20 = 0.
            Pep20 = data_array[2][0]
            y0 = [mRNA10, Pep10, mRNA20, Pep20]  # Initial condition
            VariableName = ['mRNA1', 'Pep1', 'mRNA2', 'Pep2'] # Variables Name
            VarIndex =[VariableName.index('Pep1'), VariableName.index('Pep2')]  #get the index for mRNA and RFP

            ### Number of Parameters to be optimized
            numParam = 4
            
            ParamName = ['syn_Pep','deg_Pep','syn_mRNA1', 'syn_mRNA2']
            ParamUnits = ['min-1','min-1', 'molL-1min-1', 'molL-1min-1']
            
            # Global bounds - In order of (syn_Pep, deg_Pep, syn_mRNA1, syn_mRNA2, syn_mRNA3])
            # Global bounds - In order of (a_mRNA, a_Pep, y_Pep)
            param0_global = [(0, 0.02), (0.001, 0.02), (5e-8, 5e-7),
                             (5e-8, 5e-7)] #Diff Evo

            # Local bounds
            LB = [1e-6, 0.001, 1e-10, 1e-10]
            UB = [0.02, 0.02, 5e-7, 5e-7]

            if TotalDataSet >= 3:
                mRNA30 = 0.
                Pep30 = data_array[3][0]
                y0 = y0 + [mRNA30, Pep30]  # Initial condition
                VariableName = VariableName + ['mRNA3', 'Pep3'] # Variables Name
                VarIndex =VarIndex + [VariableName.index('Pep3')]  #get the index for all Peptides

                ### Number of Parameters to be optimized
                numParam = 5
                
                ParamName = ParamName + ['syn_mRNA3']
                ParamUnits = ParamUnits + ['molL-1min-1']
                # Global bounds - In order of (syn_Pep, deg_Pep, syn_mRNA1, syn_mRNA2, syn_mRNA3])
                # Global bounds - In order of (a_mRNA, a_Pep, y_Pep)
                param0_global = param0_global + [(5e-8, 5e-7)] #Diff Evo

                # Local bounds
                LB = LB + [1e-10]
                UB = UB + [5e-7]

                if TotalDataSet >= 4:
                    mRNA40 = 0.
                    Pep40 = data_array[4][0]
                    y0 = y0 + [mRNA40, Pep40]  # Initial condition
                    VariableName = VariableName + ['mRNA4', 'Pep4'] # Variables Name
                    VarIndex = VarIndex + [VariableName.index('Pep4')]  #get the index for all Peptides

                    ### Number of Parameters to be optimized
                    numParam = 6
                    
                    ParamName = ParamName + ['syn_mRNA4']
                    ParamUnits = ParamUnits + ['molL-1min-1']
                    # Global bounds - In order of (syn_Pep, deg_Pep, syn_mRNA1, syn_mRNA2, syn_mRNA3])
                    # Global bounds - In order of (a_mRNA, a_Pep, y_Pep)
                    param0_global = param0_global + [(5e-8, 5e-7)] #Diff Evo

                    # Local bounds
                    LB = LB + [1e-10]
                    UB = UB + [5e-7]

                    if TotalDataSet >= 5:
                        mRNA50 = 0.
                        Pep50 = data_array[5][0]
                        y0 = y0 + [mRNA50, Pep50]  # Initial condition
                        VariableName = VariableName + ['mRNA5', 'Pep5'] # Variables Name
                        VarIndex = VarIndex + [VariableName.index('Pep5')]  #get the index for all Peptides

                        ### Number of Parameters to be optimized
                        numParam = 7
                        
                        ParamName = ParamName + ['syn_mRNA5']
                        ParamUnits = ParamUnits + ['molL-1min-1']
                    
                        # Global bounds - In order of (syn_Pep, deg_Pep, syn_mRNA1, syn_mRNA2, syn_mRNA3])
                        param0_global = param0_global + [(5e-8, 5e-7)] #Diff Evo

                        # Local bounds
                        LB = LB + [1e-10]
                        UB = UB + [5e-7]
                        
                        if TotalDataSet >= 6:
                            mRNA60 = 0.
                            Pep60 = data_array[6][0]
                            y0 = y0 + [mRNA60, Pep60]  # Initial condition
                            VariableName = VariableName + ['mRNA6', 'Pep6'] # Variables Name
                            VarIndex = VarIndex + [VariableName.index('Pep6')]  #get the index for all Peptides
    
                            ### Number of Parameters to be optimized
                            numParam = 8
                            
                            ParamName = ParamName + ['syn_mRNA6']
                            ParamUnits = ParamUnits + ['molL-1min-1']
                        
                            # Global bounds - In order of (syn_Pep, deg_Pep, syn_mRNA1, syn_mRNA2, syn_mRNA3])
                            param0_global = param0_global + [(5e-8, 5e-7)] #Diff Evo
    
                            # Local bounds
                            LB = LB + [1e-10]
                            UB = UB + [5e-7]
                            
        elif (SystemType == 'MultiDoubleFixRBSKMat'):

            mRNA10 = 0.
            Pep10 = 0.
            Pepm10 = data_array[1][0]
            mRNA20 = 0.
            Pep20 = 0.
            Pepm20 = data_array[2][0]
            y0 = [mRNA10, Pep10, Pepm10, mRNA20, Pep20, Pepm20]  # Initial condition
            VariableName = ['mRNA1', 'Pep1', 'Pepm1', 'mRNA2', 'Pep2', 'Pepm2'] # Variables Name
            VarIndex =[VariableName.index('Pepm1'), VariableName.index('Pepm2')]  #get the index for mRNA and RFP

            ### Number of Parameters to be optimized
            numParam = 5
            
            ParamName = ['syn_Pep','deg_Pep','syn_mRNA1', 'syn_mRNA2', 'Kmature']
            ParamUnits = ['min-1','min-1', 'molL-1min-1', 'molL-1min-1', 'min-1']
            
            # Global bounds - In order of (syn_Pep, deg_Pep, syn_mRNA1, syn_mRNA2, syn_mRNA3])
            # Global bounds - In order of (a_mRNA, a_Pep, y_Pep)
            param0_global = [(0, 0.02), (0.001, 0.02), (5e-8, 5e-7),
                             (5e-8, 5e-7), (0.001, 1)] #Diff Evo

            # Local bounds
            LB = [1e-6, 0.001, 1e-10, 1e-10, 0.002]
            UB = [0.02, 0.02, 5e-7, 5e-7, 1]

            if TotalDataSet >= 3:
                mRNA30 = 0.
                Pep30 = 0.
                Pepm30 = data_array[3][0]
                y0 = y0 + [mRNA30, Pep30, Pepm30]  # Initial condition
                VariableName = VariableName + ['mRNA3', 'Pep3', 'Pepm3'] # Variables Name
                VarIndex =VarIndex + [VariableName.index('Pepm3')]  #get the index for all Peptides

                ### Number of Parameters to be optimized
                numParam = 6
                
                ParamName = ParamName + ['syn_mRNA3']
                ParamUnits = ParamUnits + ['molL-1min-1']
                # Global bounds - In order of (syn_Pep, deg_Pep, syn_mRNA1, syn_mRNA2, syn_mRNA3])
                # Global bounds - In order of (a_mRNA, a_Pep, y_Pep)
                param0_global = param0_global + [(5e-8, 5e-7)] #Diff Evo

                # Local bounds
                LB = LB + [1e-10]
                UB = UB + [5e-7]

                if TotalDataSet >= 4:
                    mRNA40 = 0.
                    Pep40 = 0.
                    Pepm40 = data_array[4][0]
                    y0 = y0 + [mRNA40, Pep40, Pepm40]  # Initial condition
                    VariableName = VariableName + ['mRNA4', 'Pep4', 'Pepm4'] # Variables Name
                    VarIndex = VarIndex + [VariableName.index('Pepm4')]  #get the index for all Peptides

                    ### Number of Parameters to be optimized
                    numParam = 7
                    
                    ParamName = ParamName + ['syn_mRNA4']
                    ParamUnits = ParamUnits + ['molL-1min-1']
                    # Global bounds - In order of (syn_Pep, deg_Pep, syn_mRNA1, syn_mRNA2, syn_mRNA3])
                    # Global bounds - In order of (a_mRNA, a_Pep, y_Pep)
                    param0_global = param0_global + [(5e-8, 5e-7)] #Diff Evo

                    # Local bounds
                    LB = LB + [1e-10]
                    UB = UB + [5e-7]

                    if TotalDataSet >= 5:
                        mRNA50 = 0.
                        Pep50 = 0.
                        Pepm50 = data_array[5][0]
                        y0 = y0 + [mRNA50, Pep50, Pepm50]  # Initial condition
                        VariableName = VariableName + ['mRNA5', 'Pep5', 'Pepm5'] # Variables Name
                        VarIndex = VarIndex + [VariableName.index('Pepm5')]  #get the index for all Peptides

                        ### Number of Parameters to be optimized
                        numParam = 8
                        
                        ParamName = ParamName + ['syn_mRNA5']
                        ParamUnits = ParamUnits + ['molL-1min-1']
                    
                        # Global bounds - In order of (syn_Pep, deg_Pep, syn_mRNA1, syn_mRNA2, syn_mRNA3])
                        param0_global = param0_global + [(5e-8, 5e-7)] #Diff Evo

                        # Local bounds
                        LB = LB + [1e-10]
                        UB = UB + [5e-7]
                        
                        if TotalDataSet >= 6:
                            mRNA60 = 0.
                            Pep60 = 0.
                            Pepm60 = data_array[6][0]
                            y0 = y0 + [mRNA60, Pep60, Pepm60]  # Initial condition
                            VariableName = VariableName + ['mRNA6', 'Pep6', 'Pepm6'] # Variables Name
                            VarIndex = VarIndex + [VariableName.index('Pepm6')]  #get the index for all Peptides
    
                            ### Number of Parameters to be optimized
                            numParam = 9
                            
                            ParamName = ParamName + ['syn_mRNA6']
                            ParamUnits = ParamUnits + ['molL-1min-1']
                        
                            # Global bounds - In order of (syn_Pep, deg_Pep, syn_mRNA1, syn_mRNA2, syn_mRNA3])
                            param0_global = param0_global + [(5e-8, 5e-7)] #Diff Evo
    
                            # Local bounds
                            LB = LB + [1e-10]
                            UB = UB + [5e-7]
                        
        elif (SystemType == 'MultiSingleFixRBS'):

            Pep10 = data_array[1][0]
            Pep20 = data_array[2][0]
            y0 = [Pep10, Pep20]  # Initial condition
            VariableName = ['Pep1', 'Pep2'] # Variables Name
            VarIndex =[VariableName.index('Pep1'), VariableName.index('Pep2')]  #get the index for mRNA and RFP

            ### Number of Parameters to be optimized
            numParam = 3
            
            ParamName = ['deg_Pep','syn_Pep1', 'syn_Pep2']
            ParamUnits = ['min-1', 'molL-1min-1', 'molL-1min-1']
            
            # Global bounds - In order of (syn_Pep, deg_Pep, syn_mRNA1, syn_mRNA2, syn_mRNA3])
            # Global bounds - In order of (a_mRNA, a_Pep, y_Pep)
            param0_global = [(0.001, 0.02), (5e-8, 5e-7), (5e-8, 5e-7)] #Diff Evo

            # Local bounds
            LB = [0.001, 1e-10, 1e-10]
            UB = [0.02, 5e-7, 5e-7]

            if TotalDataSet >= 3:
                Pep30 = data_array[3][0]
                y0 = y0 + [Pep30]  # Initial condition
                VariableName = VariableName + ['Pep3'] # Variables Name
                VarIndex =VarIndex + [VariableName.index('Pep3')]  #get the index for all Peptides

                ### Number of Parameters to be optimized
                numParam = 4
                
                ParamName = ParamName + ['syn_Pep3']
                ParamUnits = ParamUnits + ['molL-1min-1']
                # Global bounds - In order of (syn_Pep, deg_Pep, syn_mRNA1, syn_mRNA2, syn_mRNA3])
                # Global bounds - In order of (a_mRNA, a_Pep, y_Pep)
                param0_global = param0_global + [(5e-8, 5e-7)] #Diff Evo

                # Local bounds
                LB = LB + [1e-10]
                UB = UB + [5e-7]

                if TotalDataSet >= 4:
                    Pep40 = data_array[4][0]
                    y0 = y0 + [Pep40]  # Initial condition
                    VariableName = VariableName + ['Pep4'] # Variables Name
                    VarIndex = VarIndex + [VariableName.index('Pep4')]  #get the index for all Peptides

                    ### Number of Parameters to be optimized
                    numParam = 5
                    
                    ParamName = ParamName + ['syn_Pep4']
                    ParamUnits = ParamUnits + ['molL-1min-1']
                    # Global bounds - In order of (syn_Pep, deg_Pep, syn_mRNA1, syn_mRNA2, syn_mRNA3])
                    # Global bounds - In order of (a_mRNA, a_Pep, y_Pep)
                    param0_global = param0_global + [(5e-8, 5e-7)] #Diff Evo

                    # Local bounds
                    LB = LB + [1e-10]
                    UB = UB + [5e-7]

                    if TotalDataSet >= 5:
                        Pep50 = data_array[5][0]
                        y0 = y0 + [Pep50]  # Initial condition
                        VariableName = VariableName + ['Pep5'] # Variables Name
                        VarIndex = VarIndex + [VariableName.index('Pep5')]  #get the index for all Peptides

                        ### Number of Parameters to be optimized
                        numParam = 6
                        
                        ParamName = ParamName + ['syn_Pep5']
                        ParamUnits = ParamUnits + ['molL-1min-1']
                    
                        # Global bounds - In order of (syn_Pep, deg_Pep, syn_mRNA1, syn_mRNA2, syn_mRNA3])
                        param0_global = param0_global + [(5e-8, 5e-7)] #Diff Evo

                        # Local bounds
                        LB = LB + [1e-10]
                        UB = UB + [5e-7]
                        
                        if TotalDataSet >= 6:
                            Pep60 = data_array[6][0]
                            y0 = y0 + [Pep60]  # Initial condition
                            VariableName = VariableName + ['Pep6'] # Variables Name
                            VarIndex = VarIndex + [VariableName.index('Pep6')]  #get the index for all Peptides
    
                            ### Number of Parameters to be optimized
                            numParam = 7
                            
                            ParamName = ParamName + ['syn_Pep6']
                            ParamUnits = ParamUnits + ['molL-1min-1']
                        
                            # Global bounds - In order of (syn_Pep, deg_Pep, syn_mRNA1, syn_mRNA2, syn_mRNA3])
                            param0_global = param0_global + [(5e-8, 5e-7)] #Diff Evo
    
                            # Local bounds
                            LB = LB + [1e-10]
                            UB = UB + [5e-7]
                            
        elif (SystemType == 'MultiSingleFixRBSKMat'):
            Pep10 = 0.
            Pepm10 = data_array[1][0]
            Pep20 = 0.
            Pepm20 = data_array[2][0]
            y0 = [Pep10, Pepm10, Pep20, Pepm20]  # Initial condition
            VariableName = ['Pep1', 'Pepm1', 'Pep2', 'Pepm2'] # Variables Name
            VarIndex =[VariableName.index('Pepm1'), VariableName.index('Pepm2')]  #get the index for mRNA and RFP

            ### Number of Parameters to be optimized
            numParam = 4
            
            ParamName = ['deg_Pep','syn_Pep1', 'syn_Pep2', 'Kmature']
            ParamUnits = ['min-1', 'molL-1min-1', 'molL-1min-1', 'min-1']
            
            # Global bounds - In order of (syn_Pep, deg_Pep, syn_mRNA1, syn_mRNA2, syn_mRNA3])
            # Global bounds - In order of (a_mRNA, a_Pep, y_Pep)
            param0_global = [(0.001, 0.02), (5e-8, 5e-7), (5e-8, 5e-7), (0.001, 1)] #Diff Evo

            # Local bounds
            LB = [0.001, 1e-10, 1e-10, 0.001]
            UB = [0.02, 5e-7, 5e-7, 1]

            if TotalDataSet >= 3:
                Pep30 = 0.
                Pepm30 = data_array[3][0]
                y0 = y0 + [Pep30, Pepm30]  # Initial condition
                VariableName = VariableName + ['Pep3', 'Pepm3'] # Variables Name
                VarIndex =VarIndex + [VariableName.index('Pepm3')]  #get the index for all Peptides

                ### Number of Parameters to be optimized
                numParam = 5
                
                ParamName = ParamName + ['syn_Pep3']
                ParamUnits = ParamUnits + ['molL-1min-1']
                # Global bounds - In order of (syn_Pep, deg_Pep, syn_mRNA1, syn_mRNA2, syn_mRNA3])
                # Global bounds - In order of (a_mRNA, a_Pep, y_Pep)
                param0_global = param0_global + [(5e-8, 5e-7)] #Diff Evo

                # Local bounds
                LB = LB + [1e-10]
                UB = UB + [5e-7]

                if TotalDataSet >= 4:
                    Pep40 = 0.
                    Pepm40 = data_array[4][0]
                    y0 = y0 + [Pep40, Pepm40]  # Initial condition
                    VariableName = VariableName + ['Pep4', 'Pepm4'] # Variables Name
                    VarIndex = VarIndex + [VariableName.index('Pepm4')]  #get the index for all Peptides

                    ### Number of Parameters to be optimized
                    numParam = 6
                    
                    ParamName = ParamName + ['syn_Pep4']
                    ParamUnits = ParamUnits + ['molL-1min-1']
                    # Global bounds - In order of (syn_Pep, deg_Pep, syn_mRNA1, syn_mRNA2, syn_mRNA3])
                    # Global bounds - In order of (a_mRNA, a_Pep, y_Pep)
                    param0_global = param0_global + [(5e-8, 5e-7)] #Diff Evo

                    # Local bounds
                    LB = LB + [1e-10]
                    UB = UB + [5e-7]

                    if TotalDataSet >= 5:
                        Pep50 = 0.
                        Pepm50 = data_array[5][0]
                        y0 = y0 + [Pep50, Pepm50]  # Initial condition
                        VariableName = VariableName + ['Pep5', 'Pepm5'] # Variables Name
                        VarIndex = VarIndex + [VariableName.index('Pepm5')]  #get the index for all Peptides

                        ### Number of Parameters to be optimized
                        numParam = 7
                        
                        ParamName = ParamName + ['syn_Pep5']
                        ParamUnits = ParamUnits + ['molL-1min-1']
                    
                        # Global bounds - In order of (syn_Pep, deg_Pep, syn_mRNA1, syn_mRNA2, syn_mRNA3])
                        param0_global = param0_global + [(5e-8, 5e-7)] #Diff Evo

                        # Local bounds
                        LB = LB + [1e-10]
                        UB = UB + [5e-7]
                        
                        if TotalDataSet >= 6:
                            Pep60 = 0.
                            Pepm60 = data_array[6][0]
                            y0 = y0 + [Pep60, Pepm60]  # Initial condition
                            VariableName = VariableName + ['Pep6', 'Pepm6'] # Variables Name
                            VarIndex = VarIndex + [VariableName.index('Pepm6')]  #get the index for all Peptides
    
                            ### Number of Parameters to be optimized
                            numParam = 8
                            
                            ParamName = ParamName + ['syn_Pep6']
                            ParamUnits = ParamUnits + ['molL-1min-1']
                        
                            # Global bounds - In order of (syn_Pep, deg_Pep, syn_mRNA1, syn_mRNA2, syn_mRNA3])
                            param0_global = param0_global + [(5e-8, 5e-7)] #Diff Evo
    
                            # Local bounds
                            LB = LB + [1e-10]
                            UB = UB + [5e-7]

        elif (SystemType == 'MultiDoubleFixPromoter'):
            mRNA10 = 0.
            Pep10 = data_array[1][0]
            mRNA20 = 0.
            Pep20 = data_array[2][0]
            y0 = [mRNA10, Pep10, mRNA20, Pep20]  # Initial condition
            VariableName = ['mRNA1', 'Pep1', 'mRNA2', 'Pep2'] # Variables Name
            VarIndex =[VariableName.index('Pep1'), VariableName.index('Pep2')]  #get the index for all Peptides

            ### Number of Parameters to be optimized
            numParam = 4
            
            ParamName = ['syn_mRNA', 'deg_Pep', 'syn_Pep1', 'syn_Pep2']
            ParamUnits = ['molL-1min-1', 'min-1', 'min-1', 'min-1']
            
            # Global bounds - In order of (syn_mRNA, deg_Pep, syn_Pep1, syn_Pep2)
            param0_global = [(5e-8, 5e-7), (0.001, 0.02), (0, 0.02), (0, 0.02)] #Diff Evo

            # Local bounds
            LB = [1e-10, 0.001, 1e-6, 1e-6]
            UB = [5e-7, 0.02, 0.02, 0.02]

            if TotalDataSet >= 3:
                mRNA30 = 0.
                Pep30 = data_array[3][0]
                y0 = y0 + [mRNA30, Pep30]  # Initial condition
                VariableName = VariableName + ['mRNA3', 'Pep3'] # Variables Name
                VarIndex = VarIndex + [VariableName.index('Pep3')]  #get the index for all Peptides

                ### Number of Parameters to be optimized
                numParam = 5
                
                ParamName = ParamName + ['syn_Pep3']
                ParamUnits = ParamUnits + ['min-1']
                
                # Global bounds - In order of (syn_mRNA, deg_Pep, syn_Pep1, syn_Pep2, syn_Pep3)
                param0_global = param0_global + [(0, 0.02)] #Diff Evo

                # Local bounds
                LB = LB + [1e-6]
                UB = UB + [0.02]

                if TotalDataSet >= 4:
                    mRNA40 = 0.
                    Pep40 = data_array[4][0]
                    y0 = y0 + [mRNA40, Pep40]  # Initial condition
                    VariableName = VariableName + ['mRNA4', 'Pep4'] # Variables Name
                    VarIndex = VarIndex + [VariableName.index('Pep4')]  #get the index for all Peptides

                    ### Number of Parameters to be optimized
                    numParam = 6
                    
                    ParamName = ParamName + ['syn_Pep4']
                    ParamUnits = ParamUnits + ['min-1']
                    
                    # Global bounds - In order of (syn_mRNA, deg_Pep, syn_Pep1, syn_Pep2, syn_Pep3)
                    param0_global = param0_global + [(0, 0.02)] #Diff Evo

                    # Local bounds
                    LB = LB + [1e-6]
                    UB = UB + [0.02]

                    if TotalDataSet >= 5:
                        mRNA50 = 0.
                        Pep50 = data_array[5][0]
                        y0 = y0 + [mRNA50, Pep50]  # Initial condition
                        VariableName = VariableName + ['mRNA5', 'Pep5'] # Variables Name
                        VarIndex =VarIndex + [VariableName.index('Pep5')]  #get the index for all Peptides

                        ### Number of Parameters to be optimized
                        numParam = 7
                        
                        ParamName = ParamName + ['syn_Pep5']
                        ParamUnits = ParamUnits + ['min-1']
                        
                        # Global bounds - In order of (syn_mRNA, deg_Pep, syn_Pep1, syn_Pep2, syn_Pep3)
                        param0_global = param0_global + [(0, 0.02)] #Diff Evo

                        # Local bounds
                        LB = LB + [1e-6]
                        UB = UB + [0.02]
                        
                        if TotalDataSet >= 6:
                            mRNA60 = 0.
                            Pep60 = data_array[6][0]
                            y0 = y0 + [mRNA60, Pep60]  # Initial condition
                            VariableName = VariableName + ['mRNA6', 'Pep6'] # Variables Name
                            VarIndex =VarIndex + [VariableName.index('Pep6')]  #get the index for all Peptides
    
                            ### Number of Parameters to be optimized
                            numParam = 8
                            
                            ParamName = ParamName + ['syn_Pep6']
                            ParamUnits = ParamUnits + ['min-1']
                            
                            # Global bounds - In order of (syn_mRNA, deg_Pep, syn_Pep1, syn_Pep2, syn_Pep3)
                            param0_global = param0_global + [(0, 0.02)] #Diff Evo
    
                            # Local bounds
                            LB = LB + [1e-6]
                            UB = UB + [0.02]
                            
        elif (SystemType == 'MultiDoubleFixPromoterKMat'):
            mRNA10 = 0.
            Pep10 = 0.
            Pepm10 = data_array[1][0]
            mRNA20 = 0.
            Pep20 = 0.
            Pepm20 = data_array[2][0]
            y0 = [mRNA10, Pep10, Pepm10, mRNA20, Pep20, Pepm20]  # Initial condition
            VariableName = ['mRNA1', 'Pep1', 'Pepm1', 'mRNA2', 'Pep2', 'Pepm2'] # Variables Name
            VarIndex =[VariableName.index('Pepm1'), VariableName.index('Pepm2')]  #get the index for all Peptides

            ### Number of Parameters to be optimized
            numParam = 5
            
            ParamName = ['syn_mRNA', 'deg_Pep', 'syn_Pep1', 'syn_Pep2', 'Kmature']
            ParamUnits = ['molL-1min-1', 'min-1', 'min-1', 'min-1', 'min-1']
            
            # Global bounds - In order of (syn_mRNA, deg_Pep, syn_Pep1, syn_Pep2)
            param0_global = [(5e-8, 5e-7), (0.001, 0.02), (0, 0.02), (0, 0.02), (0.001, 1)] #Diff Evo

            # Local bounds
            LB = [1e-10, 0.001, 1e-6, 1e-6, 0.002]
            UB = [5e-7, 0.02, 0.02, 0.02, 1]

            if TotalDataSet >= 3:
                mRNA30 = 0.
                Pep30 = 0.
                Pepm30 = data_array[3][0]
                y0 = y0 + [mRNA30, Pep30, Pepm30]  # Initial condition
                VariableName = VariableName + ['mRNA3', 'Pep3', 'Pepm3'] # Variables Name
                VarIndex = VarIndex + [VariableName.index('Pepm3')]  #get the index for all Peptides

                ### Number of Parameters to be optimized
                numParam = 6
                
                ParamName = ParamName + ['syn_Pep3']
                ParamUnits = ParamUnits + ['min-1']
                
                # Global bounds - In order of (syn_mRNA, deg_Pep, syn_Pep1, syn_Pep2, syn_Pep3)
                param0_global = param0_global + [(0, 0.02)] #Diff Evo

                # Local bounds
                LB = LB + [1e-6]
                UB = UB + [0.02]

                if TotalDataSet >= 4:
                    mRNA40 = 0.
                    Pep40 = 0.
                    Pepm40 = data_array[4][0]
                    y0 = y0 + [mRNA40, Pep40, Pepm40]  # Initial condition
                    VariableName = VariableName + ['mRNA4', 'Pep4', 'Pepm4'] # Variables Name
                    VarIndex = VarIndex + [VariableName.index('Pepm4')]  #get the index for all Peptides

                    ### Number of Parameters to be optimized
                    numParam = 7
                    
                    ParamName = ParamName + ['syn_Pep4']
                    ParamUnits = ParamUnits + ['min-1']
                    
                    # Global bounds - In order of (syn_mRNA, deg_Pep, syn_Pep1, syn_Pep2, syn_Pep3)
                    param0_global = param0_global + [(0, 0.02)] #Diff Evo

                    # Local bounds
                    LB = LB + [1e-6]
                    UB = UB + [0.02]

                    if TotalDataSet >= 5:
                        mRNA50 = 0.
                        Pep50 = 0.
                        Pepm50 = data_array[5][0]
                        y0 = y0 + [mRNA50, Pep50, Pepm50]  # Initial condition
                        VariableName = VariableName + ['mRNA5', 'Pep5', 'Pepm5'] # Variables Name
                        VarIndex =VarIndex + [VariableName.index('Pepm5')]  #get the index for all Peptides

                        ### Number of Parameters to be optimized
                        numParam = 8
                        
                        ParamName = ParamName + ['syn_Pep5']
                        ParamUnits = ParamUnits + ['min-1']
                        
                        # Global bounds - In order of (syn_mRNA, deg_Pep, syn_Pep1, syn_Pep2, syn_Pep3)
                        param0_global = param0_global + [(0, 0.02)] #Diff Evo

                        # Local bounds
                        LB = LB + [1e-6]
                        UB = UB + [0.02]
                        
                        if TotalDataSet >= 6:
                            mRNA60 = 0.
                            Pep60 = 0.
                            Pepm60 = data_array[6][0]
                            y0 = y0 + [mRNA60, Pep60, Pepm60]  # Initial condition
                            VariableName = VariableName + ['mRNA6', 'Pep6', 'Pepm6'] # Variables Name
                            VarIndex =VarIndex + [VariableName.index('Pepm6')]  #get the index for all Peptides
    
                            ### Number of Parameters to be optimized
                            numParam = 9
                            
                            ParamName = ParamName + ['syn_Pep6']
                            ParamUnits = ParamUnits + ['min-1']
                            
                            # Global bounds - In order of (syn_mRNA, deg_Pep, syn_Pep1, syn_Pep2, syn_Pep3)
                            param0_global = param0_global + [(0, 0.02)] #Diff Evo
    
                            # Local bounds
                            LB = LB + [1e-6]
                            UB = UB + [0.02]

        # run Global optimizer
        OptimizerType1 = 'Global'
        result_diffevo = differential_evolution\
            (self.ComputeSSE, param0_global, args=(y0, data_array, Time_interval, SystemType, VarIndex, OptimizerType1, TotalDataSet))

        # run Local Optimizer (Nelder Mead)
        OptimizerType2 = 'Local'
        param0_local = np.zeros(numParam)
        for i in range(0, numParam):
            param0_local[i] = result_diffevo.x[i]

        result_NM = cNM.constrNM(self.ComputeSSE, param0_local,LB,UB,args=(y0, data_array, Time_interval, SystemType, VarIndex, OptimizerType2, TotalDataSet),
                                 xtol= 1e-15, full_output=True)

        # Optimized Parameters
        param_optimized = np.zeros(numParam)
        for i in range(0, numParam):
            param_optimized[i] = result_NM['xopt'][i]

        return(param_optimized, sse_global, y0, VariableName, ParamName, ParamUnits)

    # Plot CSV and Model Data #
    def plotData_Combined(self, SystemType, Variable, y0, raw_data_header, rfp_data, Data_stddev, param, TotalDataSet):
        ### timespan for Model results (t)
        t_start = rfp_data[0][0]
        t_end = rfp_data[0][-1]
        dt = 10  # minutes
        timestep = int((t_end / dt) + 1)
        t = np.linspace(t_start, t_end, timestep)

        # Time grid == no. of times to report solution
        rfp_data_numrows = np.size(rfp_data, 0)
        rfp_data_numcols = np.size(rfp_data, 1)

        mRNA = np.zeros((timestep, TotalDataSet), dtype=object)
        Peptide = np.zeros((timestep, TotalDataSet), dtype=object)
        Peptidem = np.zeros((timestep, TotalDataSet), dtype=object)

        Variable_mappings = {
                'mRNA': mRNA,
                'Peptide': Peptide,
                'Peptidem': Peptidem,
            }

        # initiate empty list to store all the arrays
        VariableMatrix = [None]*len(Variable)

        for i in range(0, len(Variable)):
            VariableMatrix[i] = Variable_mappings[Variable[i]]

        solveODE_Name = '_'.join(('solveODE', SystemType))
        solveODEfun = self.select_function(solveODE_Name) #convert string to function name

#        NummRNA = sum('mRNA' in s for s in Variable) #count the number of mRNA data
#        NumPep = sum('Peptide' in s for s in Variable) #count the number of Peptide data

        # Integrate the ODE equations
        ODEsoln = odeint(solveODEfun, y0, t, args=(param, TotalDataSet))
        for j in range(0, len(VariableMatrix)):
            for i in range(0, TotalDataSet):  # Iterates through DataSet
                if (len(VariableMatrix) == 1) and (TotalDataSet > 1):
                    VariableMatrix[j][:,i] = ODEsoln[:,j+i]
                elif (len(VariableMatrix) == 3):
                    VariableMatrix[j][:,i] = ODEsoln[:,j+3*i]
                else:
                    VariableMatrix[j][:,i] = ODEsoln[:,j+2*i]
                    #VariableMatrix[j][k][i] = ODEsoln[k][j]
                    
        ### Retrive the ODEs in String from the corresponding solveODEfun
        ODEstring = solveODEfun(y0, t[0], param, TotalDataSet, 'GetODE')

        ### Defining time array from Row 0
        time = rfp_data[0]

        ### Plot RFP Data vs Time ###
        fig = plt.figure(figsize=(5,3.6))
        ax = fig.add_axes([0.16,0.16,0.8,0.78])
        plt.rc('font', size=16)  # controls default text sizes
        ## CSV Data (time)
        #style = ['o']
        style = ['^','*','>','D','<','d','p','o','h','+','s','x','v','.','H']
        for i in range (1, rfp_data_numrows):
            #plt.plot(time, rfp_data[i], linestyle='None', marker = style[i-1], markersize = 4)
            plt.errorbar(time, rfp_data[i], yerr = Data_stddev[i-1], capsize = 2, linestyle='None', marker = style[i-1], markersize = 4)

        rfp_data_legend = raw_data_header[0:rfp_data_numcols]
        plt.legend(rfp_data_legend, ncol=2, loc='upper left', prop={'size': 16},frameon=False)
        axes = plt.gca()
        ymin, ymax = axes.get_ylim()
        #axes.set_ylim([0, ymax+0.2*ymax])
        axes.set_ylim([0, ymax+0.25*ymax])
        #axes.set_ylim([0, ymax+0.1*ymax]) #for multi-4x_fixpromoter_1

        ### Model Data (t)
        # Resets colour cycle for the second plot over same figure
        plt.gca().set_prop_cycle(None)
        # plot model Pep data
        if 'Peptidem' not in Variable:
            Peptideid = Variable.index('Peptide') #get the index of mRNA
            plt.plot(t, VariableMatrix[Peptideid][:,:], linewidth=2)  # Pep
        else:
            Peptidemid = Variable.index('Peptidem') #get the index of mRNA
            plt.plot(t, VariableMatrix[Peptidemid][:,:], linewidth=2)  # Pep
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
        #axes.set_ylim([0, ymax+0.2*ymax])
        axes.set_ylim([0, ymax+0.25*ymax])
        #axes.set_ylim([0, ymax+0.1*ymax]) #for multi-4x_fixpromoter_1

        ### Protein Data (t)
        fig = plt.figure(figsize=(5,3.6))
        ax = fig.add_axes([0.16,0.15,0.8,0.78])
        plt.rc('font', size=16)  # controls default text sizes
        if 'Peptidem' not in Variable:
            Peptideid = Variable.index('Peptide') #get the index of mRNA
            plt.plot(t, VariableMatrix[Peptideid][:,:], linewidth=2)  # Pep
        else:
            Peptidemid = Variable.index('Peptidem') #get the index of mRNA
            plt.plot(t, VariableMatrix[Peptidemid][:,:], linewidth=2)  # Pep
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
        axes.set_ylim([0, ymax+0.2*ymax])

        if 'mRNA' in Variable:
            ### mRNA Data (t)
            fig = plt.figure(figsize=(5,3.6))
            ax = fig.add_axes([0.16,0.15,0.8,0.78])
            plt.rc('font', size=16)  # controls default text sizes
            mRNAid = Variable.index('mRNA') #get the index of mRNA
            plt.plot(t, VariableMatrix[mRNAid][:,:], linewidth=2)  # mRNA
            #plt.title('Modelled mRNA Concentration vs Time')
            plt.xlabel('Time (min)')
            plt.ylabel('mRNA Concentration (M)')
            plt.legend(rfp_data_legend, ncol=2, loc='upper left', prop={'size': 16},frameon=False)
            # Set Y Axis Ticker to Scientific Style
            plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
            # Figure border
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            axes = plt.gca()
            ymin, ymax = axes.get_ylim()
            #axes.set_ylim([0, ymax+0.35*ymax])
            #axes.set_ylim([0, ymax+0.45*ymax])
            axes.set_ylim([0, ymax+0.5*ymax])
            
        return t, VariableMatrix, rfp_data_legend, ODEstring

    def Run_ConstitutivePlot(self, SystemType, y0, data_header, data_array, Data_stddev, param_optimized, TotalDataSet):
        if (SystemType == 'ConstDouble'):
            VariablePlot = ['mRNA', 'Peptide'] # Variables Name
        elif (SystemType == 'ConstDoubleKMat') or (SystemType == 'MultiDoubleFixPromoterKMat') or (SystemType == 'MultiDoubleFixRBSKMat'):
            VariablePlot = ['mRNA', 'Peptide', 'Peptidem'] # Variables Name

        elif (SystemType == 'ConstSingle'):
            VariablePlot = ['Peptide'] # Variables Name
            
        elif (SystemType == 'ConstSingleKMat') or (SystemType == 'MultiSingleFixRBSKMat'):
            VariablePlot = ['Peptide', 'Peptidem'] # Variables Name
            
        elif (SystemType == 'MultiDoubleFixRBS') or (SystemType == 'MultiDoubleFixPromoter'):
            VariablePlot = ['mRNA', 'Peptide'] # Variables Name

        elif (SystemType == 'MultiSingleFixRBS'):
            VariablePlot = ['Peptide'] # Variables Name
            
        else:
            print('Error in Plotting module')

        ### Calculate and plot Model results (param_optimized) ###
        t, VariableMatrix, rfp_data_legend, ODEstring = self.plotData_Combined(SystemType, VariablePlot, y0, data_header, data_array, Data_stddev, param_optimized, TotalDataSet)

        return t, VariableMatrix, rfp_data_legend, ODEstring

    def __del__(self):
      class_name = self.__class__.__name__
      print(class_name, "destroyed")
