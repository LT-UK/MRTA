
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 27 22:34:41 2020

Monte-Carlo Simulation for tradeoff wrt truncation position parameter \rho
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
# TGTA:     Truncation Greedy Task Allocation
# T3A:      Truncation Threshold Task Allocation

import funcs.algs.tgta as tgta     # TA alg: TGTA
import funcs.algs.ltgta as ltgta   # TA alg: LTGTA

import funcs.algs.t3a as t3a       # TA alg: T3A
import funcs.algs.lt3a as lt3a     # TA alg: LT3A
import funcs.algs.ttbta as ttbta   # TA alg: TTBTA

# =============================================================================


def runTradeoff_Rho():
    
    Tasks = [j for j in range(init.Nt)]   
    Agents = [a for a in range(init.Na)]  # Na: The number of Agents
    
    init.Dist_mat = [[0 for j in range(init.Nt)] for i in range(init.Nt)]
    
    # Progress bar color options: RED, GREEN, ORANGE, BLUE, VIOLET, CYAN, GREY.
    bar_color = pb.colors.CYAN # (default: CYAN)
    start_time = time.time()
    time_last = start_time
    item_ctr = 0
    # There will be an error with calculating 1.0//0.1. It will return 9.0 instead of 10.0. 
    # Solution: use (1.0+0.1/10)//0.1.
    total_item = init.MC_ctr*((init.Rho_end-init.Rho_start+init.Rho_step/10)//init.Rho_step+1)
    # Print initial status.
    pb.printProgressBar(item_ctr, total_item, start_time, bar_color)
    
    # Change the number of agents       
    for rho in np.arange(init.Rho_start, init.Rho_end+init.Rho_step, init.Rho_step):
       
        init.Rhos.append(rho)
         
# TGTA                 
        if init.TGTA.en:
            init.TGTA.resetSum()
# LTGTA                 
        if init.LTGTA.en:
            init.LTGTA.resetSum()         
# T3A                 
        if init.T3A.en:
            init.T3A.resetSum()
# LT3A                 
        if init.LT3A.en:
            init.LT3A.resetSum()
# TTBTA                 
        if init.TTBTA.en:
            init.TTBTA.resetSum()


        # Monte-Carlo Loop
        for counter in range(init.MC_ctr):       
            # Tasks' and Agents' positions (x,y) /km
            Tasks_pos = [
                    [round(x,3) for x in 
                     np.random.uniform(low=0.0, high=init.L, size=init.Nt)], 
                    [round(y,3) for y in 
                     np.random.uniform(low=0.0, high=init.L, size=init.Nt)]]     
            
            # get the distance matrix of all tasks
            for i in range(init.Nt-1):
                for j in range(i+1,init.Nt):
                    init.Dist_mat[i][j] = round(vf.getDistance(i,j,Tasks_pos),3)
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
    #            M = np.random.uniform(m_low,m_high,(Na,Nt))
                init.M = [[round(m,3) for m in 
                           np.random.uniform(init.m_low,init.m_high,init.Nt)] 
                            for a in range(init.Na)
                         ]
            # non-monotone
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
                m1 = init.ms_low*np.ones(init.Na) \
                        + (init.ms_high-init.ms_low)*np.eye(init.Na)
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
            
            # Run TA algorithms.
## TGTA                   
            if init.TGTA.en:    
                selected, values, total_value, dt, steps, evs \
                    = tgta.runTGTA(Agents,Tasks,rho)
                init.TGTA.utility_sum += total_value
                init.TGTA.dt_sum += dt
                init.TGTA.steps_sum += steps
                init.TGTA.evs_sum += evs 
## LTGTA                   
            if init.LTGTA.en:    
                selected, values, total_value, dt, steps, evs \
                    = ltgta.runLTGTA(Agents,Tasks,rho)
                init.LTGTA.utility_sum += total_value
                init.LTGTA.dt_sum += dt
                init.LTGTA.steps_sum += steps
                init.LTGTA.evs_sum += evs 
## T3A                   
            if init.T3A.en:    
                selected, values, total_value, dt, steps, evs \
                    = t3a.runT3A(Agents,Tasks,rho,init.eps)
                init.T3A.utility_sum += total_value
                init.T3A.dt_sum += dt
                init.T3A.steps_sum += steps
                init.T3A.evs_sum += evs 
## LT3A                   
            if init.LT3A.en:    
                selected, values, total_value, dt, steps, evs \
                    = lt3a.runLT3A(Agents,Tasks,rho,init.eps)
                init.LT3A.utility_sum += total_value
                init.LT3A.dt_sum += dt
                init.LT3A.steps_sum += steps
                init.LT3A.evs_sum += evs             
## TTBTA                   
            if init.TTBTA.en:    
                selected, values, total_value, dt, steps, evs \
                    = ttbta.runTTBTA(Agents,Tasks,rho,init.eps)
                init.TTBTA.utility_sum += total_value
                init.TTBTA.dt_sum += dt
                init.TTBTA.steps_sum += steps
                init.TTBTA.evs_sum += evs 


            # Update progress bar
            item_ctr += 1
            time_current = time.time()
            delta_time = time_current - time_last            
            if delta_time > 1:
                time_last = time_current
                pb.printProgressBar(item_ctr, total_item, start_time, bar_color)                     

        # end for

        # Get average #
## TGTA                 
        if init.TGTA.en:
            init.TGTA.utilities.append(init.TGTA.utility_sum / init.MC_ctr)
            init.TGTA.dts.append(init.TGTA.dt_sum / init.MC_ctr)
            init.TGTA.steps.append(init.TGTA.steps_sum / init.MC_ctr)
            init.TGTA.evs.append(init.TGTA.evs_sum / init.MC_ctr)
## LTGTA                 
        if init.LTGTA.en:
            init.LTGTA.utilities.append(init.LTGTA.utility_sum / init.MC_ctr)
            init.LTGTA.dts.append(init.LTGTA.dt_sum / init.MC_ctr)
            init.LTGTA.steps.append(init.LTGTA.steps_sum / init.MC_ctr)
            init.LTGTA.evs.append(init.LTGTA.evs_sum / init.MC_ctr)
## T3A                 
        if init.T3A.en:
            init.T3A.utilities.append(init.T3A.utility_sum / init.MC_ctr)
            init.T3A.dts.append(init.T3A.dt_sum / init.MC_ctr)
            init.T3A.steps.append(init.T3A.steps_sum / init.MC_ctr)
            init.T3A.evs.append(init.T3A.evs_sum / init.MC_ctr)
## LT3A                 
        if init.LT3A.en:
            init.LT3A.utilities.append(init.LT3A.utility_sum / init.MC_ctr)
            init.LT3A.dts.append(init.LT3A.dt_sum / init.MC_ctr)
            init.LT3A.steps.append(init.LT3A.steps_sum / init.MC_ctr)
            init.LT3A.evs.append(init.LT3A.evs_sum / init.MC_ctr)
## TTBTA                 
        if init.TTBTA.en:
            init.TTBTA.utilities.append(init.TTBTA.utility_sum / init.MC_ctr)
            init.TTBTA.dts.append(init.TTBTA.dt_sum / init.MC_ctr)
            init.TTBTA.steps.append(init.TTBTA.steps_sum / init.MC_ctr)
            init.TTBTA.evs.append(init.TTBTA.evs_sum / init.MC_ctr)
            
    # end for 
    
    # Print final progress status.
    pb.printProgressBar(item_ctr, total_item, start_time, bar_color)
    
    # Save settings and algorithms' output.
    if init.monotonicity:
        store_path = os.path.abspath("output") + "/tradeoff/truncation/monotone"
    else:
        # Non-monotone is not available for truncation actually.
        store_path = os.path.abspath("output") + "/tradeoff/truncation/non_monotone"
    store.storeData(store_path)
    
    
# =============================================================================
# 
# =============================================================================
print("----- tradeoff_Truncation.py is loaded -----")

    