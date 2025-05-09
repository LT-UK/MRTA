#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 18 15:58:04 2020

Task allocation algorithm: 
    STTA: Sample Threshold Task Allocation

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
def runSTTA(Agents, Tasks, Pr, eps):
    '''
    Sample Threshold Task Allocation.
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
    start_time = time.time()
    consensus_steps = 0 # count the total consensus times
    n_evs = 0 # count the number of utility function evaluations
    
    r = len(Tasks)
    values = [0 for a in Agents] # record the final function values for each agent
    
    # Sampling
    remained = []    
    for a in Agents:
        samples = []
        for j in Tasks:
            if np.random.choice([0, 1],p=[1-Pr,Pr]):
                samples.append(j)
        remained.append(samples)
    
    # Create the empty allocation for each agent start from empty set. The first number is agent id.
    selected = [[a] for a in Agents] 
    local_max_mgvs = [] # max marginal values for each agent
    local_max_js = []   # max elements for all agents # j represents task id
    for a in Agents:
        mgvs = []
        for j in remained[a]:
            mgvs.append(vf.getMGV(j,selected[a],Tasks))
            n_evs += 1
        if remained[a]: # not empty
            local_max_mgvs.append(np.max(mgvs))
            local_max_js.append(remained[a][np.argmax(mgvs)])
        else: # empty
            local_max_mgvs.append(0)
            local_max_js.append("none")
    global_max_mgv = np.max(local_max_mgvs)
    agent_id = np.argmax(local_max_mgvs) # find the agent that win the bid
    global_max_j = local_max_js[agent_id] # find the element that helps to win the bid
    A = [agent_id]
    J = [global_max_j]
    W = [global_max_mgv]
    d = global_max_mgv
    theta = d # initial threshold
    final_threshold = eps/r*d
    
    # main loop
    # TODO: Update main loop considering buffers
    while theta >= final_threshold:
        flag_agent = [1 for a in Agents] # stop flag for each agent
        # If no element is larger than the current threshold theta, then stop searching under this threshold.
        while 1 in flag_agent: 
            # Find the eligible task for each agent.
            for a in Agents:
                omega = 0 # store the mgv that is larger than the current threshold
                if flag_agent[a]:  
                    abandon = []
                    for j in remained[a]:
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
                    # end for
                    for j in abandon:
                        remained[a].remove(j)
                if omega == 0: # There is no task whose mgv is larger than the current threshold.
                    flag_agent[a] = 0
            # end for
            if len(A) > 0:
                consensus_steps += 1 # new consensus round
                # If no agent can provide an eligible task, then consensus is not required.
            # Add the outbit tasks to corresponding agents and remove these tasks from remained sets for all agents       
            for a in Agents:
                if a in A:
                    selected[a].append(J[A.index(a)])
                    values[a] += W[A.index(a)]
                for j in J:
                    if j in remained[a]:
                        remained[a].remove(j)
            # Clear sets for the new iteration.
            A = []
            J = []
            W = []       
        # end while
        # If the remained sets are empty for all agents, stop the main loop.
        stop_flag = 1
        for a in Agents:
            if remained[a]:  # If anyone not empty, don't stop main loop.
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
print("----- stta.py is loaded -----")
