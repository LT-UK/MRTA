#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 18 15:58:04 2020

Task allocation algorithm: 
    LT3A: Lazy Truncation Threshold Task Allocation

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
def runLT3A(Agents, Tasks, k, eps):
    '''
    Lazy Truncation Threshold Task Allocation.
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
    
    # initialization
    r = len(Tasks)
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
    global_max_mgv = max(max_mgvs)
    
    # Truncate the remained task list for each agent
    for a in Agents:
        for i in range(len(sorted_remained[a])): # do not use len(Tasks) or j in Tasks
            if sorted_mgvs[a][i] >= mgv_k:
                sorted_front_remained[a].append(sorted_remained[a][i])
                sorted_front_mgvs[a].append(sorted_mgvs[a][i])
            else:
                break
    
    # init threshold variables
    A = []
    J = []
    W = []
    d = global_max_mgv
    theta = d
    final_threshold = eps/r*d
    
    # main loop: LDTTA based on the sorted_front_remained and sorted_front_mgvs
    while theta >= final_threshold:
        flag_agent = [1 for a in Agents] # stop flag for each agent
        # If no element is larger than the current threshold theta, then stop searching under this threshold.
        while 1 in flag_agent: 
            # Find the eligible task for each agent.
            for a in Agents:
                omega = 0 # store the mgv that is larger than the current threshold
                if len(sorted_front_remained[a]) > 0:
                    while sorted_front_mgvs[a][0] >= theta:
                        # Re-calculate marginal values for all agents because all selected tasks in previous  
                        # iterations has been allocated.
                        new_mgv = vf.getMGV(sorted_front_remained[a][0],selected[a],Tasks)
                        n_evs += 1
                        if new_mgv >= theta:
                            omega = new_mgv
                            j_star = sorted_front_remained[a][0]
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
                        elif new_mgv < final_threshold:
                            sorted_front_remained[a].pop(0) 
                            sorted_front_mgvs[a].pop(0)
                            if len(sorted_front_mgvs[a]) < 1:
                                break
                        else: #resort sorted_front_mgvs and sorted_front_remained
                            sorted_front_mgvs[a][0] = new_mgv
                            new_mgvs = sorted_front_mgvs[a]
                            sorted_front_mgvs[a] = list(np.sort(new_mgvs)[::-1]) # list() convert array to list
                            sorted_indexes = np.argsort(new_mgvs)[::-1]
                            sorted_front_remained[a] = [sorted_front_remained[a][i] for i in sorted_indexes]
                if omega == 0: # There is no task whose mgv is larger than the current threshold.
                    flag_agent[a] = 0
            if len(A) > 0: # new consensus round
                consensus_steps += 1 
                # If no agent can provide an eligible task, then consensus is not required.
            # Add the outbit tasks to corresponding agents and remove these tasks from remained sets for all agents       
            for a in Agents:
                if a in A:
                    selected[a].append(J[A.index(a)])
                    values[a] += W[A.index(a)]
                for j in J:
                    if j in sorted_front_remained[a]:
                        j_index = sorted_front_remained[a].index(j)
                        sorted_front_mgvs[a].pop(j_index)
                        sorted_front_remained[a].pop(j_index)
                    
            # Clear sets for the new iteration.
            A = []
            J = []
            W = []       
        # If the remained sets are empty for all agents, stop the main loop.
        stop_flag = 1
        for a in Agents:
            if sorted_front_remained[a]:  # If anyone not empty, don't stop main loop.
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
print("----- lt3a.py is loaded -----")

