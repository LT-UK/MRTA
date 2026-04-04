#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 01 17:27:04 2020

Task allocation algorithm: 
    STBTA: Sample Threshold Bundle Task Allocation

@author: Teng Li
lt.uk@outlook.com
United Kingdom
All Rights Reserved
"""

import numpy as np
import time
import funcs.vfunc.vf as vf   # value function and getMGV

# =============================================================================
#    When applying lazy strategy on bundle building, mgvs are not accurate
#    because tasks in bundles can be removed due to conflict.
#    Therefore, lazy strategy are not suitable for bundle related algorithms.
# =============================================================================

# TODO: Alg is not ready yet.

def runSTBTA(Agents, Tasks, Pr, eps):
    '''
    Sample Threshold Bundle Task Allocation.
    
    If more than one agent bids for a task, the agent with shorter selection will win the auction.
    If agents have the same length of selections, then allocate according to numerical sequence of their ids. 
    
    Input:    
        Agents: [list] all agents' ids
        Tasks: [list] all tasks' ids
        Pr: [float] sampling probabilty
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
print("----- stbta.py is loaded -----")
