#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 18 15:58:04 2020

Task allocation algorithm: 
    LTGTA: Lazy Truncation Greedy Task Allocation

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
def runLTGTA(Agents, Tasks, k):
    '''
    Lazy Truncation Greedy Task Allocation.
    Input:    
        Agents: [list] all agents' ids
        Tasks: [list] all tasks' ids
        k: [float] truncation parameter determaining the truncation position
    Output:
        selected: [list 2D] selected tasks' ids by each agent 
        values: [list] function value for each agent
        total_value: [float] total function value i.e. sum of all individual agent's function value
        dt: [float] consuming time, unit: sec
        consensus_steps: [int] the number of consensus steps
        n_evs: [int] the number of function evaluations
    '''
    start_time = time.time()
    consensus_steps = 0 # count the total consensus times
    n_evs = 0 # count the number of utility function evaluations
    
    # initialization
    N_min = len(Tasks)
    values = [0 for a in Agents] # record the final function values for each agent
    remained = [list(Tasks) for a in Agents] # start from ground set   
    # Create the empty allocation for each agent start from empty set. The first number is agent id.
    selected = [[a] for a in Agents] 
    sorted_mgvs = []
    sorted_remained = []
    local_max_mgvs = [] # max marginal values for each agent
    local_max_js = [] # # max elements for all agents  # j represents task id    
    mgvs_table = [[] for j in Tasks] # by tasks
    max_mgvs = [] # maximum mgv for each task
    sorted_front_remained = [[] for a in Agents]
    sorted_front_mgvs = [[] for a in Agents]
    
    # initial evaluation of all task agent pairs
    for a in Agents:
        mgvs = []
        for j in remained[a]:
            mgv = vf.getMGV(j,selected[a],Tasks)
            mgvs.append(mgv)
            mgvs_table[j].append(mgv)
            n_evs += 1
        # Sort the mgvs and tasks in a decreasing order.
        sorted_mgvs.append(list(np.sort(mgvs)[::-1]))   # decreasing
        sorted_indexes = np.argsort(mgvs)[::-1]
        sorted_remained.append([remained[a][j] for j in sorted_indexes])
        local_max_mgvs.append(sorted_mgvs[a][0])
        local_max_js.append(sorted_remained[a][0])
    
    # find the truncation position
    for j in Tasks:
        max_mgvs.append(max(mgvs_table[j])) # find the max mgv for each task
    mgv_k = k*min(max_mgvs) # Get the truncation mgv    
#    global_max_mgv = max(max_mgvs)
    
    # Truncate the remained task list for each agent
    for a in Agents:
        for i in range(len(sorted_remained[a])): # do not use len(Tasks) or j in Tasks
            if sorted_mgvs[a][i] >= mgv_k:
                sorted_front_remained[a].append(sorted_remained[a][i])
                sorted_front_mgvs[a].append(sorted_mgvs[a][i])
            else:
                break
    
    searching = [0 for a in Agents]
    
    # main loop
    for i in range(N_min):        
        for a in Agents:   
            while searching[a]:
                if len(sorted_front_remained[a]) > 1:               
                    #recalculate marginal value for the first element
                    new_mgv = vf.getMGV(sorted_front_remained[a][0],selected[a],Tasks)            
                    n_evs += 1
                    if new_mgv >= sorted_front_mgvs[a][1]:
                        local_max_mgvs[a] = new_mgv
                        local_max_js[a] = sorted_front_remained[a][0]
                        searching[a] = 0 # stop the while loop
                    else:
                        new_mgvs = [new_mgv] + sorted_front_mgvs[a][1:]
                        sorted_front_mgvs[a] = list(np.sort(new_mgvs)[::-1])
                        sorted_indexes = np.argsort(new_mgvs)[::-1]
                        sorted_front_remained[a] = [sorted_front_remained[a][j] for j in sorted_indexes]              
                elif len(sorted_front_remained[a]) == 1:
                    local_max_mgvs[a] = vf.getMGV(sorted_front_remained[a][0],selected[a],Tasks)
                    n_evs += 1
                    local_max_js[a] = sorted_front_remained[a][0]
                    searching[a] = 0 # stop the while loop    
                else:
                    local_max_mgvs[a] = 0
                    local_max_js[a] = "None"
                    searching[a] = 0 # stop the while loop
            # end while
        # end for
        
        global_max_mgv = np.max(local_max_mgvs)
        agent_id = np.argmax(local_max_mgvs) # find the agent that win the bid
        global_max_u = local_max_js[agent_id] # find the element that helps to win the bid
        
        if global_max_mgv > 0:
            # Add the winning element to the corresponding agent's selected set.
            selected[agent_id].append(global_max_u)
            values[agent_id] += global_max_mgv
            # Remove the winning element from all sampled remained sets.
            for a in Agents:
                if global_max_u in sorted_front_remained[a]:
                    if global_max_u == sorted_front_remained[a][0]:
                        searching[a] = 1
                    max_index = sorted_front_remained[a].index(global_max_u)
                    sorted_front_remained[a].pop(max_index)
                    sorted_front_mgvs[a].pop(max_index)
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
print("----- ltgta.py is loaded -----")
