
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 27 22:34:41 2020

Monte-Carlo Simulation for tradeoff wrt epsilon
Random task positions, task values, task-agent match fitness factors.

@author: Teng Li
lt.uk@outlook.com
United Kingdom
All Rights Reserved
"""

# =====================  Load Built-in Modules  ===============================
import numpy as np
import time
import os
import math
# =============================================================================

# ==========================  Load Supporting Modules  ========================
# Do NOT change the import order.
import init                          # algorithms initialisation
import funcs.vfunc.vf as vf          # load value function and getMGV
import funcs.sys.store as store      # store and restore settings and outputs
import funcs.sys.progressbar as pb   # print progress bar
# =============================================================================

# ==========================  Load TA Algorithms  =============================
# (prefix 'L' represents 'Lazy')
# TTBTA:    Truncation Threshold Bundle Task Allocation
# STTA:     Sample Threshold Task Allocation
# STBTA:    Sample Threshold Bundle Task Allocation

import funcs.algs.sga as sga       # TA alg: SGA
#import funcs.algs.lsga as lsga     # TA alg: LSGA

import funcs.algs.cbba as cbba     # TA alg: CBBA

import funcs.algs.stta as stta     # TA alg: STTA
import funcs.algs.lstta as lstta   # TA alg: LSTTA
import funcs.algs.stbta as stbta   # TA alg: STBTA
# =============================================================================


def runTradeoff_Eps_Vr():
    
    Tasks = [j for j in range(init.Nt)] 
    
    Agents = [a for a in range(init.Na)]
    
    init.Dist_mat = [[0 for j in range(init.Nt)] for i in range(init.Nt)]
    
    # Tasks' and Agents' positions (x,y) /km
    Tasks_pos = [
            [round(x,3) for x in 
             np.random.uniform(low=0.0, high=init.L, size=init.Nt)], 
            [round(y,3) for y in 
             np.random.uniform(low=0.0, high=init.L, size=init.Nt)]]
    
    if init.monotonicity:        
        Agents_pos = [
            [round(x,3) for x in np.random.uniform(low=0.0, high=init.L, size=init.Na)], 
            [round(y,3) for y in np.random.uniform(low=0.0, high=init.L, size=init.Na)]]       
    else:
        Agents_pos = [Tasks_pos[0][0:init.Na+1], Tasks_pos[1][0:init.Na+1]]
    
    # Get distances between agents and tasks.
    init.Dist_AT = [[0 for j in Tasks] for a in Agents]        
    for a in Agents:
        for j in Tasks:
            init.Dist_AT[a][j] = np.round(vf.getDistAT(a, j, Agents_pos, Tasks_pos), 3)      
    
    # get the distance matrix of all tasks
    for i in range(init.Nt-1):
        for j in range(i+1,init.Nt):
            init.Dist_mat[i][j] = np.round(vf.getDistance(i, j, Tasks_pos), 3)
            init.Dist_mat[j][i] = init.Dist_mat[i][j]
    
    # monotone case
    if init.monotonicity:
        # random importance factors
        init.V = [round(v,3) for v in 
                  np.random.uniform(low=init.v_low, 
                                    high=init.v_high, 
                                    size=init.Nt)
                 ]
        # random task-agent match fitness factors
        init.M = [[round(m,3) for m in 
                   np.random.uniform(init.m_low,init.m_high,init.Nt)
                  ] for a in range(init.Na_end)
                 ]  
    # non-monotone case
    else:
        # special tasks
        v1 = [round(v,3) for v in 
              np.random.uniform(low=init.vs_low, 
                                high=init.vs_high, 
                                size=init.Na)
             ]
        # normal tasks
        v2 = [round(v,3) for v in 
              np.random.uniform(low=init.v_low, 
                                high=init.v_high, 
                                size=init.Nt-init.Na)
             ]
        init.V = v1 + v2
        # special tasks
        m1 = init.ms_low*np.ones(init.Na) + \
            (init.ms_high-init.ms_low)*np.eye(init.Na)
        # normal tasks
        m2 = [[round(m,3) for m in 
               np.random.uniform(init.m_low,init.m_high,init.Nt-init.Na)] 
                for a in range(init.Na)
             ]
        init.M = np.concatenate((m1, m2), axis=1) 
        # inter-task dependences Nt*Nt
        init.X = np.zeros((init.Nt, init.Nt))
        for i in range(init.Nt-1):
            for j in range(i+1, init.Nt):
                init.X[i][j] = round(math.e**(init.V[i]*init.V[j]),4)
                init.X[j][i] = init.X[i][j]
    
    
    # Progress bar color options: RED, GREEN, ORANGE, BLUE, VIOLET, CYAN, GREY.
    bar_color = pb.colors.CYAN # (default: CYAN)
    start_time = time.time()
    time_last = start_time
    item_ctr = 0
    total_item = init.MC_ctr*((init.Eps_end-init.Eps_start)//init.Eps_step+1)
    # Print initial status.
    pb.printProgressBar(item_ctr, total_item, start_time, bar_color)
    
    # ************* test ****************
    if init.CBBA.en:    
        selected, values, total_value, dt, steps, evs \
            = cbba.runCBBA(Agents,Tasks)        
        init.CBBA.utilities = [total_value for i in range(9)]
        init.CBBA.dts = [dt for i in range(9)]
        init.CBBA.steps = [steps for i in range(9)]
        init.CBBA.evs = [evs/10**init.evs_scale for i in range(9)]
        
    if init.SGA.en:    
        selected, values, total_value, dt, steps, evs \
            = sga.runSGA(Agents,Tasks)         
        init.SGA.utilities = [total_value for i in range(9)]
        init.SGA.dts = [dt for i in range(9)]
        init.SGA.steps = [steps for i in range(9)]
        init.SGA.evs = [evs/10**init.evs_scale for i in range(9)]
    # ***********************************
    
    # Change the threshold decreasing parameter       
    for eps in np.arange(init.Eps_start, init.Eps_end+init.Eps_step, init.Eps_step):
       
        init.Epses.append(eps)         
# STTA                 
        if init.STTA.en:
            init.STTA.resetTemp()
# LSTTA                             
        if init.LSTTA.en:
            init.LSTTA.resetTemp()
# STBTA                 
        if init.STBTA.en:
            init.STBTA.resetTemp()

        # Monte-Carlo Loop
        for counter in range(init.MC_ctr):                      
            # Run TA algorithms.           
## STTA                   
            if init.STTA.en:    
                selected, values, total_value, dt, steps, evs \
                    = stta.runSTTA(Agents,Tasks,init.Pr,eps)
                init.STTA.utilities_temp.append(total_value)
                init.STTA.dts_temp.append(dt)
                init.STTA.steps_temp.append(steps)
                init.STTA.evs_temp.append(evs)               
## LSTTA                             
            if init.LSTTA.en:    
                selected, values, total_value, dt, steps, evs \
                    = lstta.runLSTTA(Agents,Tasks,init.Pr,eps)
                init.LSTTA.utilities_temp.append(total_value)
                init.LSTTA.dts_temp.append(dt)
                init.LSTTA.steps_temp.append(steps)
                init.LSTTA.evs_temp.append(evs)               
## STBTA                   
            if init.STBTA.en:    
                selected, values, total_value, dt, steps, evs \
                    = stbta.runSTBTA(Agents,Tasks,init.Pr,eps)
                init.STBTA.utilities_temp.append(total_value)
                init.STBTA.dts_temp.append(dt)
                init.STBTA.steps_temp.append(steps)
                init.STBTA.evs_temp.append(evs)

            # Update progress bar
            item_ctr += 1
            time_current = time.time()
            delta_time = time_current - time_last            
            if delta_time > 1:
                time_last = time_current
                pb.printProgressBar(item_ctr, total_item, start_time, bar_color)                     

        # end for

## STTA                
        if init.STTA.en:
            init.STTA.utilities.append(init.STTA.utilities_temp)
            init.STTA.dts.append(init.STTA.dts_temp)
            init.STTA.steps.append(init.STTA.steps_temp)
            init.STTA.evs.append(init.STTA.evs_temp)
## LSTTA                         
        if init.LSTTA.en:
            init.LSTTA.utilities.append(init.LSTTA.utilities_temp)
            init.LSTTA.dts.append(init.LSTTA.dts_temp)
            init.LSTTA.steps.append(init.LSTTA.steps_temp)
            init.LSTTA.evs.append(init.LSTTA.evs_temp)
## STBTA                 
        if init.STBTA.en:
            init.STBTA.utilities.append(init.STBTA.utilities_temp)
            init.STBTA.dts.append(init.STBTA.dts_temp)
            init.STBTA.steps.append(init.STBTA.steps_temp)
            init.STBTA.evs.append(init.STBTA.evs_temp)                 
    # end for 
    
    # Print final progress status.
    pb.printProgressBar(item_ctr, total_item, start_time, bar_color)
    
    # Save settings and algorithms' output.
    if init.monotonicity:
        store_path = os.path.abspath("output") + "/tradeoff/epsilon_vr/monotone"
    else:
        store_path = os.path.abspath("output") + "/tradeoff/epsilon_vr/non_monotone"
    store.storeData(store_path)
    
    
# =============================================================================
# 
# =============================================================================
print("----- tradeoff_Threshold_vr.py is loaded -----")

    