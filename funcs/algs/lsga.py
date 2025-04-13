#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 18 15:58:04 2020

Task allocation algorithm: 
    LSGA: Lazy Sequencial Greedy Algorithm

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
def runLSGA(Agents, Tasks):
    '''
    Lazy Sequencial Greedy Algorithm for Task Allocation.
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
    start_time = time.time()
    consensus_steps = 0 # count the total consensus times
    n_evs = 0 # count the number of utility function evaluations
    
    remained = [list(Tasks) for a in Agents] # start from ground set         
    selected = [[i] for i in Agents] # start from empty set. The first number is agent id.
    N_min = len(Tasks)    
    
    sorted_mgvs = []
    sorted_remained = []   
    local_max_mgvs = [] # max marginal values for each agent
    local_max_us = [] # # max elements for all agents
    values = [0 for a in Agents] # record the final function values for each agent
    
    for a in Agents:
        mgvs = [] #marginal values
        for j in remained[a]:
            mgvs.append(vf.getMGV(j,selected[a],Tasks))
            n_evs += 1
        sorted_mgvs.append(list(np.sort(mgvs)[::-1]))    # decreasing
        sorted_indexes = np.argsort(mgvs)[::-1]
        sorted_remained.append([remained[a][j] for j in sorted_indexes]) # rearrange the remained set
        local_max_mgvs.append(sorted_mgvs[a][0])
        local_max_us.append(sorted_remained[a][0])
    
    searching = [0 for a in Agents]
    
    for i in range(N_min):        
        for a in Agents:   
            while searching[a]:
                if len(sorted_remained[a]) > 1:               
                    #recalculate marginal value for the first element
                    new_mgv = vf.getMGV(sorted_remained[a][0],selected[a],Tasks)            
                    n_evs += 1
                    if new_mgv >= sorted_mgvs[a][1]:
                        local_max_mgvs[a] = new_mgv
                        local_max_us[a] = sorted_remained[a][0]
                        searching[a] = 0 # stop the while loop
                    else:
                        new_mgvs = [new_mgv] + sorted_mgvs[a][1:]
                        sorted_mgvs[a] = list(np.sort(new_mgvs)[::-1])
                        sorted_indexes = np.argsort(new_mgvs)[::-1]
                        sorted_remained[a] = [sorted_remained[a][j] for j in sorted_indexes]              
                elif len(sorted_remained[a]) == 1:
                    local_max_mgvs[a] = vf.getMGV(sorted_remained[a][0],selected[a],Tasks)
                    n_evs += 1
                    local_max_us[a] = sorted_remained[a][0]
                    searching[a] = 0 # stop the while loop    
                else:
                    local_max_mgvs[a] = 0
                    local_max_us[a] = "None"
                    searching[a] = 0 # stop the while loop
            # end while
        # end for
        
        global_max_mgv = np.max(local_max_mgvs)
        agent_id = np.argmax(local_max_mgvs) # find the agent that win the bid
        global_max_u = local_max_us[agent_id] # find the element that helps to win the bid
        
        if global_max_mgv > 0:
            # Add the winning element to the corresponding agent's selected set.
            selected[agent_id].append(global_max_u)
            values[agent_id] += global_max_mgv
            # Remove the winning element from all sampled remained sets.
            for a in Agents:
                if global_max_u in sorted_remained[a]:
                    if global_max_u == sorted_remained[a][0]:
                        searching[a] = 1
                    max_index = sorted_remained[a].index(global_max_u)
                    sorted_remained[a].pop(max_index)
                    sorted_mgvs[a].pop(max_index)
        else:
            break
        
        consensus_steps += 1
    # end for
    
    # Get the total function value of all agents.
    total_value = np.sum(values)   

    end_time = time.time()
    dt = end_time - start_time
        
    return selected, values, total_value, dt, consensus_steps, n_evs


# =============================================================================
# 
# =============================================================================
print("----- lsga.py is loaded -----")
