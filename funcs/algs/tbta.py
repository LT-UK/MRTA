#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 18 15:58:04 2020

Task allocation algorithm: 
    TBTA: Threshold Bundle Task Allocation

@author: Teng Li
lt.uk@outlook.com
United Kingdom
All Rights Reserved
"""

import numpy as np
import time
import funcs.vfunc.vf as vf   # value function and getMGV

# =============================================================================
#    When applying lazy strategy on bundle building, mgvs are not accurate
#    because tasks in bundles can be removed due to conflict.
#    Therefore, lazy strategy are not suitable for bundle related algorithms.
# =============================================================================

def runTBTA(Agents, Tasks, eps):
    '''
    Threshold Bundle Task Allocation.
    
    If more than one agent bids for a task, the agent with shorter selection will win the auction.
    If agents have the same length of selections, then allocate according to numerical sequence of their ids. 
    
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
    # initialization
    start_time = time.time()
    consensus_steps = 0 # Count the total consensus times
    n_evs = 0 # Count the number of utility function evaluations
    
    r = len(Tasks)
    values = [0 for a in Agents] # Record the final function values for each agent
    
    remained = [list(Tasks) for a in Agents] # Start from ground set   
    # Create the empty allocation for each agent start from empty set. The first number is agent id.
    selected = [[a] for a in Agents] 
    local_max_mgvs = [] # Max marginal values for each agent
    for a in Agents:
        mgvs = []
        for j in remained[a]:
            mgvs.append(vf.getMGV(j,selected[a],Tasks))
            n_evs += 1
        local_max_mgvs.append(np.max(mgvs))
    global_max_mgv = np.max(local_max_mgvs)
    d = global_max_mgv
    theta = d # initial threshold
    final_threshold = eps/r*d
    
    # Initialise bundles.
    label_invalid = [0 for a in Agents] # Indicate the position where the tasks start to be invalid in the boundle
    bundle_mgvs = [[] for a in Agents]
    bundle_js = [[] for a in Agents]
        
    # main loop
    # TODO: Update main loop considering buffers
    while theta >= final_threshold:
        flag_agent = [1 for a in Agents] # Stop flag for each agent with the current theta
        # If no element is larger than theta, then stop searching under this threshold.
        while 1 in flag_agent: 
            flag_bundle = [1 for a in Agents] # Used for checking whether the bundle has elements
            # Build bundles
            for a in Agents:
                if flag_agent[a]: # Agent a still has task whose mgv is larger than theta
                    to_be_removed = []
                    # Initialise and reset bundles
                    bundle_mgvs[a] = []
                    bundle_js[a] = []
                    label_invalid[a] = 0
                    for j in remained[a]:  # Find the eligible task for each agent.
                        mgv = vf.getMGV(j,selected[a]+bundle_js[a],Tasks)
                        n_evs += 1
                        if mgv >= theta:
                            bundle_js[a].append(j)
                            bundle_mgvs[a].append(mgv)
                            to_be_removed.append(j)
                        elif mgv < final_threshold:
                            to_be_removed.append(j) 
                    # DO NOT remove j in the last for loop, because it will impact the remained
                    for j in to_be_removed: 
                        remained[a].remove(j)
                    if bundle_js[a]:  # len(bundle_js[a]) > 0
                        label_invalid[a] = len(bundle_js[a]) # Reset the label
                    else:
                        flag_agent[a] = 0
                        flag_bundle[a] = 0
                else:
                    flag_bundle[a] = 0
            if 1 in flag_agent:
                consensus_steps += 1 # New consensus round
            # Allocate tasks in the bundles to agents and resolve conflicts
            while 1 in flag_bundle: # Some bundles are still not empty
                # Clear J and A
                J = []
                A = []
                # Allocate tasks from the first row of the bundles
                for a in Agents:
                    # Check if the bundle of agent a has elements
                    if flag_bundle[a]: 
                        if bundle_js[a][0] not in J:
                            J.append(bundle_js[a][0])
                            A.append(a)          
                        else: # Auction
                            # Find the agent that has choosen the same task
                            j_index = J.index(bundle_js[a][0])
                            agent_id = A[j_index]
                            # Compare the numbers of selected tasks
                            if len(selected[a]) < len(selected[agent_id]):
                                A[j_index] = a  # replace agent id
                # Allocation
                for a in Agents:
                    # Check whether the bundle of agent a has elements.
                    if flag_bundle[a]:
                        if a in A: # winner
                            # The mgv of the 1st task in the bundle is still valid
                            if label_invalid[a] > 0: 
                                values[a] += bundle_mgvs[a][0]                              
                                label_invalid[a] -= 1
                            else:
                                # Recalculate mgv
                                mgv = vf.getMGV(bundle_js[a][0],selected[a],Tasks)
                                # No need to count evs here, because this is only used for recording values.
                                values[a] += mgv
                            selected[a].append(bundle_js[a][0])
                        else: # loser
                            # All mgvs in the bundle become invalid
                            label_invalid[a] = 0
                        # Remove the first element from the bundles
                        bundle_js[a].pop(0)
                        bundle_mgvs[a].pop(0)
                # Solve conflicts
                for a in Agents:
                    for j in J:
                        # Check remained tasks
                        if j in remained[a]:
                            remained[a].remove(j)
                        # Check bundles
                        if j in bundle_js[a]:
                            j_index = bundle_js[a].index(j)
                            if j_index < label_invalid[a]:
                                label_invalid[a] = j_index
                            bundle_js[a].pop(j_index)
                            bundle_mgvs[a].pop(j_index)
                    # Check bundle length
                    if not bundle_js[a]: # Same functionality with 'if len(bundle_js[a]) == 0:' but faster
                        flag_bundle[a] = 0
            
        # If the remained sets are empty for all agents, stop the main loop.
        stop_flag = 1
        for a in Agents:
            if remained[a]:  # If anyone not empty, don't stop main loop.
                stop_flag = 0 
                break
        if stop_flag:
            break # stop the main loop
        
        # Go to next threshold
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
print("----- tbta.py is loaded -----")
