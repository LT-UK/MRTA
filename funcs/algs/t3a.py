#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 18 15:58:04 2020

Task allocation algorithm: 
    T3A: Truncation Threshold Task Allocation

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
        dt: [float] consuming time /sec
        consensus_steps: [int] the number of consensus steps
        n_evs: [int] the number of function evaluations
    '''
    start_time = time.time()
    consensus_steps = 0 # count the total consensus times
    n_evs = 0 # count the number of utility function evaluations
    
    r = len(Tasks)
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
        max_mgvs.append(max(mgvs_table[j])) # find the max mgv for each task
    mgv_k = k*min(max_mgvs) # Get the truncation mgv
    
    global_max_mgv = max(max_mgvs)
    
    # Truncate the remained task list for each agent
    for a in Agents:
        for i in range(len(sorted_remained[a])):  # do not use len(Tasks) or j in Tasks
            if sorted_mgvs[a][i] >= mgv_k:
                front_remained[a].append(sorted_remained[a][i])
            else:
                break
    
    # Init threshold variables
    A = []
    J = []
    W = []
    d = global_max_mgv
    theta = d # initial threshold
    final_threshold = eps/r*d
    
    # Allocate tasks 
    # main loop
    while theta >= final_threshold:
        flag_agent = [1 for a in Agents] # stop flag for each agent
        # If no element is larger than the current threshold theta, then stop searching under this threshold.
        while 1 in flag_agent: 
            # Find the eligible task for each agent.
            for a in Agents:
                omega = 0 # store the mgv that is larger than the current threshold
                if flag_agent[a]:  
                    abandon = []
                    for j in front_remained[a]:
                        mgv_temp = vf.getMGV(j,selected[a],Tasks)
                        n_evs += 1
                        if mgv_temp >= theta:
                            omega = mgv_temp
                            j_star = j
                            if j_star not in J:
                                J.append(j_star)
                                A.append(a)
                                W.append(omega)
                            else: # Deal with the conflict.
                                j_star_index = J.index(j_star)
                                if omega > W[j_star_index]:
                                    A[j_star_index] = a
                                    W[j_star_index] = omega
                            break  # If one eligible task is found, then stop the current searching loop.
                        elif mgv_temp < final_threshold:
                            abandon.append(j)  # DO NOT remove j within for loop!!! 
                    for j in abandon:
                        front_remained[a].remove(j)
                if omega == 0: # There is no task whose mgv is larger than the current threshold.
                    flag_agent[a] = 0
            if len(A) > 0:
                consensus_steps += 1 # new consensus round
                # If no agent can provide an eligible task, then consensus is not required.
            # Add the outbit tasks to corresponding agents and remove these tasks from front_remained sets for all agents       
            for a in Agents:
                if a in A:
                    selected[a].append(J[A.index(a)])
                    values[a] += W[A.index(a)]
                for j in J:
                    if j in front_remained[a]:
                        front_remained[a].remove(j)
            # Clear sets for the new iteration.
            A = []
            J = []
            W = []       
        # If the front_remained sets are empty for all agents, stop the main loop.
        stop_flag = 1
        for a in Agents:
            if front_remained[a]:  # If anyone not empty, don't stop main loop.
                stop_flag = 0 
                break
        if stop_flag:
            break # stop the main loop
        # go to next threshold
        theta *= (1 - eps) 
    # end main loop
            
    # Get the total function value of all agents. 
    total_value = np.sum(values)      
    
    end_time = time.time()
    dt = end_time - start_time 
    
    return selected, values, total_value, dt, consensus_steps, n_evs


# =============================================================================
# 
# =============================================================================
print("----- t3a.py is loaded -----")

