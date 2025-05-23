#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 18 15:58:04 2020

Task allocation algorithm: 
    DTTA: Decreasing Threshold Task Allocation

@author: Teng Li
lt.uk@outlook.com;
United Kingdom
All Rights Reserved
"""

import numpy as np
import time
import funcs.vfunc.vf as vf   # value function and getMGV

# =============================================================================
# 
# =============================================================================
def runDTTA(Agents, Tasks, eps):
    '''
    Decreasing Threshold Task Allocation.
    Input:    
        Agents: [list] all agents' ids
        Tasks: [list] all tasks' ids
        eps: [float] epsilon, the parameter of threshold
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
    
    r = len(Tasks)
    values = [0 for a in Agents] # record the final function values for each agent
    
    # initialization
    remained = [list(Tasks) for a in Agents] # start from ground set   
    remained_theta = [list(Tasks) for a in Agents] # Buffer for remained tasks under the current theta   
    Omega = [[] for a in Agents] # Buffer for valid mgvs
    Omega_J = [[] for a in Agents] # Buffer for task ids corresponding to valid mgvs
    # Create the empty allocation for each agent start from empty set. The first number is agent id.
    selected = [[a] for a in Agents] 
    local_max_mgvs = [] # max marginal values for each agent
    local_max_js = []   # max elements for all agents # j represents task id
    for a in Agents:
        mgvs = []
        for j in remained[a]:
            mgvs.append(vf.getMGV(j,selected[a],Tasks))
        local_max_mgvs.append(np.max(mgvs))
        local_max_js.append(remained[a][np.argmax(mgvs)])
    global_max_mgv = np.max(local_max_mgvs)

    d = global_max_mgv
    theta = d # initial threshold
    final_threshold = eps/r*d
    
    # main loop
    while theta >= final_threshold:
        # Init auxiliary sets for the new iteration.
        A = []
        J = []
        W = []    
        # Find the eligible task for each agent.
        for a in Agents:
            omega = 0 # store the mgv that is larger than the current threshold
            abandon = []
            considered = []
            for j in remained_theta[a]:
                mgv_temp = vf.getMGV(j,selected[a],Tasks)
                considered.append(j) # !!! Do NOT remove considered task within for loop.
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
                else:
                    Omega[a].append(mgv_temp)
                    Omega_J[a].append(j)
            for j in abandon:
                remained[a].remove(j)
            for j in considered:
                remained_theta[a].remove(j)
                
        # Update Omega
        for a in Agents:
            for j in J:
                if j in Omega_J[a]:
                    j_idx = Omega_J[a].index(j)
                    Omega[a].pop(j_idx)
                    Omega_J[a].remove(j)
                    
        # Find the largest remaining mgv
        omega_n_max = 0 # final_threshold
        for a in Agents:
            if Omega[a]:
                omega_n_max = max(omega_n_max, max(Omega[a]))

        consensus_steps += 1 # new consensus round
        
        if J: # J not empty         
            # Add the outbit tasks to corresponding agents and remove these tasks from remained sets for all agents       
            for a in Agents:
                if a in A:
                    selected[a].append(J[A.index(a)])
                    values[a] += W[A.index(a)]
                    Omega[a] = []
                    Omega_J[a] = []
                for j in J:
                    if j in remained[a]:
                        remained[a].remove(j)
                    if j in remained_theta[a]:
                        remained_theta[a].remove(j)
                  
        else:
            # If the remained sets are empty for all agents, stop the main loop.
            if omega_n_max == 0:
                break
                 
            while theta > omega_n_max:
                theta *= (1 - eps) # go to next threshold
            for a in Agents:
                if Omega[a] and max(Omega[a]) < theta:
                    remained_theta[a] = []
                else:
                    remained_theta[a] = [j for j in remained[a]]
                    Omega[a] = []
                    Omega_J[a] = []
        
    # end main loop
            
    # Get the total function value of all agents. 
    total_value = np.sum(values)      
    
    end_time = time.time()
    dt = end_time - start_time 
    
    return selected, values, total_value, dt, consensus_steps, n_evs


# =============================================================================
# 
# =============================================================================
print("----- dtta.py is loaded -----")
