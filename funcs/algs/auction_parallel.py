#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 18 15:58:04 2020

Task allocation algorithm: 
    Parallel Auction Algorithm

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
def runAuctionParallel(Agents, Tasks):
    '''
    Sequencial Greedy Algorithm for Task Allocation.
    Input:    
        Agents: [list] all agents' ids
        Tasks: [list] all tasks' ids
    Output:
        selected: [list 2D] selected tasks' ids by each agent 
        values: [list] function value for each agent
        total_value: [float] total function value i.e. sum of all individual agent's function value
        dt: [float] consuming time /sec
        consensus_steps: [int+] the number of consensus steps
        n_evs: [int+] the number of function evaluations
    '''
    # initialization
    start_time = time.time()
    consensus_steps = 0 # count the total consensus times
    n_evs = 0 # count the number of utility function evaluations    

    selected = [[a] for a in Agents] # start from empty set. The first number is agent id.        
    values = [0 for a in Agents] # record the final function values for each agent
    
    MGV_Tasks = [[] for j in Tasks] # Record all mgv for all agents and tasks
    
    for j in Tasks:
        for a in Agents:
            MGV_Tasks[j].append(vf.getMGV(j,selected[a],Tasks))
            n_evs += 1

    # Find the largeset mgv for each task
    for j in Tasks:
        mgv_j_max = np.max(MGV_Tasks[j])
        agent_id = np.argmax(MGV_Tasks[j])

        if values[agent_id] == 0:
            # mgv is still valid
            values[agent_id] = mgv_j_max
        else:
            # Recalculate mgv
            values[agent_id] += vf.getMGV(j,selected[agent_id],Tasks)
        selected[agent_id].append(j)
    
    
    consensus_steps += 1
           
    # Get the total function value of all agents.
    total_value = np.sum(values)   

    end_time = time.time()
    dt = end_time - start_time
        
    return selected, values, total_value, dt, consensus_steps, n_evs   


# =============================================================================
# 
# =============================================================================
print("----- auction_parallel.py is loaded -----")
