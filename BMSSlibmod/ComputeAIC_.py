# -*- coding: utf-8 -*-
"""
Released on April 29, 2019

@author: Yeoh Jing Wui <bchyjw@nus.edu.sg>; Ivan Ng Kai Boon; Poh Chueh Loo <poh.chuehloo@nus.edu.sg>

The code is part of BMSS software.

Copyright (c) 2019, National University of Singapore.

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
