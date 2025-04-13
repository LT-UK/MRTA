#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 18 15:58:04 2020

Task allocation algorithm: 
    Sequencial Auction Algorithm

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
def runAuctionSequential(Agents, Tasks):
    '''
    Sequencial Auction Algorithm for Task Allocation.
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
    
    for j in Tasks:
        mgvs = [] #marginal values
        for a in Agents:           
            mgvs.append(vf.getMGV(j,selected[a],Tasks))
            n_evs += 1
            
        global_max_mgv = np.max(mgvs)
        agent_id = np.argmax(mgvs) # find the agent that wins the bid

        if global_max_mgv > 0:
            # Add the winning element to the corresponding agent's selected set.
            selected[agent_id].append(j)
            values[agent_id] += global_max_mgv
            
        else:
            break

        consensus_steps += 1
           
    # Get the total function value of all agents.
    total_value = np.sum(values)   

    end_time = time.time()
    dt = end_time - start_time
        
    return selected, values, total_value, dt, consensus_steps, n_evs   


# =============================================================================
# 
# =============================================================================
print("----- auction_sequential.py is loaded -----")
