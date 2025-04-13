#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 18 15:58:04 2020

Task allocation algorithm: 
    TGTA: Truncation Greedy Task Allocation

@author: Teng Li
tengli@cranfield.ac.uk; nliteng@foxmail.com
Cranfield University, UK
All Rights Reserved
"""

import numpy as np
import time
import funcs.vfunc.vf as vf   # value function and getMGV

# =============================================================================
# 
# =============================================================================
def runTGTA(Agents, Tasks, k):
    '''
    Decreasing Threshold Task Allocation.
    Input:    
        Agents: [list] all agents' ids
        Tasks: [list] all tasks' ids
        k: [float] truncation parameter determaining the truncation position
    Output:
        selected: [list 2D] selected tasks' ids by each agent 
        values: [list] function value for each agent
        total_value: [float] total function value i.e. sum of all individual agent's function value
        dt: [float] consuming time /sec
        consensus_steps: [int] the number of consensus steps
        n_evs: [int] the number of function evaluations
    '''
    start_time = time.time()
    consensus_steps = 0 # count the total consensus times
    n_evs = 0 # count the number of utility function evaluations
    
#    N_min = np.min([Nt, Na*Lt]) # the stop condition
    N_min = len(Tasks)
    
    # start from empty set. The first number is agent id.
    selected = [[a] for a in Agents]
    remained = [list(Tasks) for a in Agents] # start from ground set  
        
    values = [0 for a in Agents] # record the final function values for each agent 
    mgvs_table = [[] for j in Tasks] # by tasks
    max_mgvs = [] # maximum mgv for each task
    sorted_mgvs = []
    sorted_remained = []
    front_remained = [[] for a in Agents]
    
    # Get the truncation position
    for a in Agents:
        mgvs = [] #marginal values
        for j in remained[a]:
            mgv = vf.getMGV(j,selected[a],Tasks)
            mgvs.append(mgv)
            mgvs_table[j].append(mgv)
            n_evs += 1
        # Sort mgvs and tasks in descending order.
        sorted_mgvs.append(list(np.sort(mgvs)[::-1]))    # descending order
        sorted_indexes = np.argsort(mgvs)[::-1]
        sorted_remained.append([remained[a][i] for i in sorted_indexes]) # rearrange the remained set
    for j in Tasks:
        max_mgvs.append(np.max(mgvs_table[j])) # find the max mgv for each task
    mgv_k = k*np.min(max_mgvs) # Get the truncation mgv
    
    # Truncate the remained task list for each agent
    for a in Agents:
        for i in range(len(sorted_remained[a])):
            if sorted_mgvs[a][i] >= mgv_k:
                front_remained[a].append(sorted_remained[a][i])
            else:
                break

    # Allocate tasks 
    for i in range(N_min):
        local_max_mgvs = [] # max marginal values for each agent
        local_max_js = [] # # max elements for all agents
        for a in Agents:
            
            mgvs = [] #marginal values
#            delta_dists = []
            if front_remained[a]:     # len(front_remained[a]) > 0          
                #calculate marginal values for remained elements
                for j in front_remained[a]:
                    mgvs.append(vf.getMGV(j,selected[a],Tasks))
                    n_evs += 1
                local_max_mgvs.append(np.max(mgvs))
                task_index = np.argmax(mgvs)
                local_max_js.append(front_remained[a][task_index])
            else:
                local_max_mgvs.append(0)
                local_max_js.append("None")
        global_max_mgv = np.max(local_max_mgvs)
        agent_id = np.argmax(local_max_mgvs) # find the agent that wins the bid
        global_max_u = local_max_js[agent_id] # find the element that helps to win the bid
        if global_max_mgv > 0:
            # Add the winning element to the corresponding agent's selected set.
            selected[agent_id].append(global_max_u)
            values[agent_id] += global_max_mgv
            # Remove the winning element from all sampled remained sets.
            for a in Agents:
                if global_max_u in front_remained[a]:
                    front_remained[a].remove(global_max_u)
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
print("----- tgta.py is loaded -----")
