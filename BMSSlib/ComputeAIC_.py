# -*- coding: utf-8 -*-
"""
Created on Fri Dec 28 17:02:29 2018

@author: jingwui
"""

#!/usr/bin/env python

### This module is to compute the AIC value based on SSE, sample size, 
### and number of fitted parameters

import numpy as np

def run_AIC(SSE_array, sample_size, num_param):
    SSE_len = len(SSE_array)
    AIC_array = []  #It is initialized as a list, not a numpy array

    for i in range (0,SSE_len):
        AIC_value = sample_size * np.log(SSE_array[i]/sample_size) + 2*num_param[i]
        AIC_array.append(AIC_value)

    return(AIC_array)