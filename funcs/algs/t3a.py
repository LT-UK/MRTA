#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 18 15:58:04 2020

Task allocation algorithm: 
    T3A: Truncation Threshold Task Allocation

@author: Teng Li
lt.uk@outlook.com
United Kingdom
All Rights Reserved
"""

import numpy as np
import time
import funcs.vfunc.vf as vf   # value function and getMGV

# =============================================================================
# 
# =============================================================================

# TODO: Alg is not ready yet.

def runT3A(Agents, Tasks, k, eps):
    '''
    Truncation Threshold Task Allocation.
    Input:    
        Agents: [list] all agents' ids
        Tasks: [list] all tasks' ids
        k: [float] truncation parameter determaining the truncation position
        eps: [float] epsilon, the parameter of threshold
    Output:
        selected: [list 2D] selected tasks' ids by each agent 
        values: [list] function value for each agent
        total_value: [float] total function value i.e. sum of all individual agent's function value
        dt: [float] consuming time, unit: sec
        consensus_steps: [int] the number of consensus steps
        n_evs: [int] the number of function evaluations
    '''



# =============================================================================
# 
# =============================================================================
print("----- t3a.py is loaded -----")

