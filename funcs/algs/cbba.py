#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 18 15:58:04 2020

Task allocation algorithm: 
    CBBA: Consensus Based Bundle Algorithm

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
def runCBBA(Agents, Tasks):
    '''
    Input:    
        Agents: [list] all agents' ids
        Tasks: [list] all tasks' ids        
    Output:
        selected: [list 2D] selected tasks' ids by each agent 
        values: [list] function value for each agent
        total_value: [float] total function value i.e. sum of all individual agent's function value
        dt: [float] consuming time, unit: sec
        consensus_steps: [int+] the number of consensus steps
        n_evs: [int+] the number of function evaluations
    '''
    start_time = time.time()   
    consensus_steps = 0 # count the total consensus times
    n_evs = 0 # count the number of utility function evaluations
    
    remained_total = list(Tasks)
    
    # start from empty set. The first number is agent id.
    selected = [[a] for a in Agents]
        
    flag_invalid = [0 for a in Agents] # indicate the position where the tasks start to be invalid in the boundle
    bundle_mgvs = [[] for a in Agents]
    bundle_us = [[] for a in Agents]
    values = [0 for a in Agents] # record the final function values for each agent
    
    while remained_total:   # len(remained_total) > 0    
        local_max_mgvs = [] # reset max marginal values for all agent
        local_max_us = [] # reset max elements for all agents
        consensus_flag = 0
        # check whether bundle rebuilding is needed, any rebuilding requires a new consensus
        for a in Agents:
            if not flag_invalid[a]: # rebuild the boundle
                consensus_flag = 1 # new consensus is required
                remained_temp = [u for u in remained_total]
                # clear bundles
                bundle_mgvs[a] = []
                bundle_us[a] = []
                while remained_temp: # len(remained_temp) > 0
                    mgvs = [] # marginal values
                    for u in remained_temp:
                        mgvs.append(vf.getMGV(u,selected[a]+bundle_us[a],Tasks))
                        n_evs += 1
                    bundle_mgvs[a].append(np.max(mgvs))
                    task_index = np.argmax(mgvs)
                    bundle_us[a].append(remained_temp[task_index])
                    remained_temp.pop(task_index)
                flag_invalid[a] = len(bundle_us[a]) - 1 # reset the flag
        if consensus_flag: # new consensus is required
            consensus_steps += 1
        # find the global max mgv and corresponding u
        for a in Agents:
            local_max_mgvs.append(bundle_mgvs[a][0])
            local_max_us.append(bundle_us[a][0]) 
        global_max_mgv = np.max(local_max_mgvs)
        agent_id = np.argmax(local_max_mgvs) # find the agent that win the bid
        global_max_u = local_max_us[agent_id] # find the element that helps to win the bid
        # In non-monotone case, check whether the max mgv is positive or not.
        if global_max_mgv > 0:
            # Add the winning element to the corresponding agent's selected set.
            selected[agent_id].append(global_max_u)
            values[agent_id] += global_max_mgv
            bundle_us[agent_id].pop(0)
            bundle_mgvs[agent_id].pop(0)   
            remained_total.remove(global_max_u)
            # Mark the flag for all agents
            for a in Agents:
                if a == agent_id: # winner
                    flag_invalid[a] -= 1 
                else: # losesr
                    max_u_index = bundle_us[a].index(global_max_u)
                    if max_u_index < flag_invalid[a]:
                        flag_invalid[a] = max_u_index
        else:
            break   
    
    # Get the total function value of all agents.    
    total_value = np.sum(values)      
    
    end_time = time.time()
    dt = end_time - start_time 
    
    return selected, values, total_value, dt, consensus_steps, n_evs


# =============================================================================
# 
# =============================================================================
print("----- cbba.py is loaded -----")