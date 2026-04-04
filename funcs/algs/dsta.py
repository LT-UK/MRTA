#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 18 15:58:04 2020

Task allocation algorithm: 
    DSTA: Decentralised Sample based Task Allocation

@author: Teng Li
lt.uk@outlook.com
United Kingdom
All Rights Reserved
"""

import numpy as np
import time
import init
import funcs.vfunc.vf as vf   # value function and getMGV

# =============================================================================
# 
# =============================================================================
def runDSTA(Agents, Tasks, Pr):
    '''
    Decentralised Sample based Task Allocation
    Input:    
        Agents: [list] all agents' ids
        Tasks: [list] all tasks' ids
        Pr: [float] sampling probability
    Output:
        selected: [list 2D] selected tasks' ids by each agent 
        values: [list] function value for each agent
        total_value: [float] total function value i.e. sum of all individual agent's function value
        dt: [float] consuming time, unit: sec
        consensus_steps: [int+] the number of consensus steps
        n_evs: [int+] the number of function evaluations
    '''
    # initialization
    start_time = time.time()
    consensus_steps = 0 # count the total consensus times
    n_evs = 0 # count the number of utility function evaluations
    
    N_min = len(Tasks)
    # N_min = np.min([Nt, Na*Lt]) # the stop condition
    sampled_remained = []
    admissible_tasks = [[] for a in Agents]
    for a in Agents:
        samples = []
        for j in Tasks:
            if init.scenario == 2:
                if init.Dist_AT[a][j] <= init.ro:
                    admissible_tasks[a].append(j)
                    if np.random.choice([0, 1],p=[1-Pr,Pr]):
                        samples.append(j)   
            else:
                if np.random.choice([0, 1],p=[1-Pr,Pr]):
                    samples.append(j)
        sampled_remained.append(samples)
    
    if init.scenario == 2:
        tasks = admissible_tasks
    else:
        tasks = [Tasks for a in Agents]
    
    selected = [[i] for i in Agents]       
    values = [0 for a in Agents] # record the final function values for each agent
    
    for i in range(N_min):
        local_max_mgvs = [] # max marginal values for each agent
        local_max_us = [] # # max elements for all agents
        for a in Agents:
            mgvs = [] #marginal values
            if len(sampled_remained[a]) > 0:               
                #calculate marginal values for remained elements
                for j in sampled_remained[a]:
                    mgvs.append(vf.getMGV(j,selected[a],tasks[a]))
                    n_evs += 1
                local_max_mgvs.append(np.max(mgvs))
                task_index = np.argmax(mgvs)
                local_max_us.append(sampled_remained[a][task_index])
            else:
                local_max_mgvs.append(0)
                local_max_us.append("None")
        global_max_mgv = np.max(local_max_mgvs)
        agent_id = np.argmax(local_max_mgvs) # find the agent that wins the bid
        global_max_u = local_max_us[agent_id] # find the element that helps to win the bid
        if global_max_mgv > 0:
            # Add the winning element to the corresponding agent's selected set.
            selected[agent_id].append(global_max_u)
            values[agent_id] += global_max_mgv
            # Remove the winning element from all sampled remained sets.
            for a in Agents:
                if global_max_u in sampled_remained[a]:
                    sampled_remained[a].remove(global_max_u)
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
print("----- dsta.py is loaded -----")
