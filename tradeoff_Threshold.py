
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
# DTTA:     Decreasing Threshold Task Allocation
# TBTA:     Threshold Bundle Task Allocation
# T3A:      Truncation Threshold Task Allocation
# TTBTA:    Truncation Threshold Bundle Task Allocation
# STTA:     Sample Threshold Task Allocation
# STBTA:    Sample Threshold Bundle Task Allocation


import funcs.algs.dtta as dtta     # TA alg: DTTA
import funcs.algs.ldtta as ldtta   # TA alg: LDTTA
import funcs.algs.tbta as tbta     # TA alg: TBTA

import funcs.algs.t3a as t3a       # TA alg: T3A
import funcs.algs.lt3a as lt3a     # TA alg: LT3A
import funcs.algs.ttbta as ttbta   # TA alg: TTBTA

import funcs.algs.stta as stta     # TA alg: STTA
import funcs.algs.lstta as lstta   # TA alg: LSTTA
import funcs.algs.stbta as stbta   # TA alg: STBTA
# =============================================================================


def runTradeoff_Eps():
    
    Tasks = [j for j in range(init.Nt)]   
    Agents = [a for a in range(init.Na)]  # Na: The number of Agents
    
    init.Dist_mat = [[0 for j in range(init.Nt)] for i in range(init.Nt)]
    
    # Progress bar color options: RED, GREEN, ORANGE, BLUE, VIOLET, CYAN, GREY.
    bar_color = pb.colors.CYAN # (default: CYAN)
    start_time = time.time()
    time_last = start_time
    item_ctr = 0
    total_item = init.MC_ctr*((init.Eps_end-init.Eps_start)//init.Eps_step+1)
    # Print initial status.
    pb.printProgressBar(item_ctr, total_item, start_time, bar_color)
    
    # Change the number of agents       
    for eps in np.arange(init.Eps_start, init.Eps_end+init.Eps_step, init.Eps_step):
       
        init.Epses.append(eps)
         
# DTTA                        
        if init.DTTA.en:
            init.DTTA.resetSum()
# LDTTA                           
        if init.LDTTA.en:
            init.LDTTA.resetSum()
# TBTA                 
        if init.TBTA.en:
            init.TBTA.resetSum()            
# T3A                 
        if init.T3A.en:
            init.T3A.resetSum()
# LT3A                 
        if init.LT3A.en:
            init.LT3A.resetSum()
# TTBTA                 
        if init.TTBTA.en:
            init.TTBTA.resetSum()
# STTA                 
        if init.STTA.en:
            init.STTA.resetSum()
# LSTTA                             
        if init.LSTTA.en:
            init.LSTTA.resetSum()
# STBTA                 
        if init.STBTA.en:
            init.STBTA.resetSum()

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
## DTTA                
            if init.DTTA.en:    
                selected, values, total_value, dt, steps, evs \
                    = dtta.runDTTA(Agents,Tasks,eps)                    
                init.DTTA.utility_sum += total_value
                init.DTTA.dt_sum += dt
                init.DTTA.steps_sum += steps
                init.DTTA.evs_sum += evs 
## LDTTA               
            if init.LDTTA.en:    
                selected, values, total_value, dt, steps, evs \
                    = ldtta.runLDTTA(Agents,Tasks,eps)
                init.LDTTA.utility_sum += total_value
                init.LDTTA.dt_sum += dt
                init.LDTTA.steps_sum += steps
                init.LDTTA.evs_sum += evs 
## T3A                   
            if init.T3A.en:    
                selected, values, total_value, dt, steps, evs \
                    = t3a.runT3A(Agents,Tasks,init.rho,eps)
                init.T3A.utility_sum += total_value
                init.T3A.dt_sum += dt
                init.T3A.steps_sum += steps
                init.T3A.evs_sum += evs 
## LT3A                   
            if init.LT3A.en:    
                selected, values, total_value, dt, steps, evs \
                    = lt3a.runLT3A(Agents,Tasks,init.rho,eps)
                init.LT3A.utility_sum += total_value
                init.LT3A.dt_sum += dt
                init.LT3A.steps_sum += steps
                init.LT3A.evs_sum += evs             
## STTA                   
            if init.STTA.en:    
                selected, values, total_value, dt, steps, evs \
                    = stta.runSTTA(Agents,Tasks,init.Pr,eps)
                init.STTA.utility_sum += total_value
                init.STTA.dt_sum += dt
                init.STTA.steps_sum += steps
                init.STTA.evs_sum += evs                
## LSTTA                             
            if init.LSTTA.en:    
                selected, values, total_value, dt, steps, evs \
                    = lstta.runLSTTA(Agents,Tasks,init.Pr,eps)
                init.LSTTA.utility_sum += total_value
                init.LSTTA.dt_sum += dt
                init.LSTTA.steps_sum += steps
                init.LSTTA.evs_sum += evs                
## TBTA                   
            if init.TBTA.en:    
                selected, values, total_value, dt, steps, evs \
                    = tbta.runTBTA(Agents,Tasks,eps)
                init.TBTA.utility_sum += total_value
                init.TBTA.dt_sum += dt
                init.TBTA.steps_sum += steps
                init.TBTA.evs_sum += evs 
## TTBTA                   
            if init.TTBTA.en:    
                selected, values, total_value, dt, steps, evs \
                    = ttbta.runTTBTA(Agents,Tasks,init.rho,eps)
                init.TTBTA.utility_sum += total_value
                init.TTBTA.dt_sum += dt
                init.TTBTA.steps_sum += steps
                init.TTBTA.evs_sum += evs 
## STBTA                   
            if init.STBTA.en:    
                selected, values, total_value, dt, steps, evs \
                    = stbta.runSTBTA(Agents,Tasks,init.Pr,eps)
                init.STBTA.utility_sum += total_value
                init.STBTA.dt_sum += dt
                init.STBTA.steps_sum += steps
                init.STBTA.evs_sum += evs 

            # Update progress bar
            item_ctr += 1
            time_current = time.time()
            delta_time = time_current - time_last            
            if delta_time > 1:
                time_last = time_current
                pb.printProgressBar(item_ctr, total_item, start_time, bar_color)                     

        # end for

        # Get average #
## DTTA  
        if init.DTTA.en:
            init.DTTA.utilities.append(init.DTTA.utility_sum / init.MC_ctr)
            init.DTTA.dts.append(init.DTTA.dt_sum / init.MC_ctr)
            init.DTTA.steps.append(init.DTTA.steps_sum / init.MC_ctr)
            init.DTTA.evs.append(init.DTTA.evs_sum / init.MC_ctr)
## LDTTA 
        if init.LDTTA.en:
            init.LDTTA.utilities.append(init.LDTTA.utility_sum / init.MC_ctr)
            init.LDTTA.dts.append(init.LDTTA.dt_sum / init.MC_ctr)
            init.LDTTA.steps.append(init.LDTTA.steps_sum / init.MC_ctr)
            init.LDTTA.evs.append(init.LDTTA.evs_sum / init.MC_ctr)
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
## STTA                
        if init.STTA.en:
            init.STTA.utilities.append(init.STTA.utility_sum / init.MC_ctr)
            init.STTA.dts.append(init.STTA.dt_sum / init.MC_ctr)
            init.STTA.steps.append(init.STTA.steps_sum / init.MC_ctr)
            init.STTA.evs.append(init.STTA.evs_sum / init.MC_ctr)
## LSTTA                         
        if init.LSTTA.en:
            init.LSTTA.utilities.append(init.LSTTA.utility_sum / init.MC_ctr)
            init.LSTTA.dts.append(init.LSTTA.dt_sum / init.MC_ctr)
            init.LSTTA.steps.append(init.LSTTA.steps_sum / init.MC_ctr)
            init.LSTTA.evs.append(init.LSTTA.evs_sum / init.MC_ctr)
## TBTA                 
        if init.TBTA.en:
            init.TBTA.utilities.append(init.TBTA.utility_sum / init.MC_ctr)
            init.TBTA.dts.append(init.TBTA.dt_sum / init.MC_ctr)
            init.TBTA.steps.append(init.TBTA.steps_sum / init.MC_ctr)
            init.TBTA.evs.append(init.TBTA.evs_sum / init.MC_ctr)
## TTBTA                 
        if init.TTBTA.en:
            init.TTBTA.utilities.append(init.TTBTA.utility_sum / init.MC_ctr)
            init.TTBTA.dts.append(init.TTBTA.dt_sum / init.MC_ctr)
            init.TTBTA.steps.append(init.TTBTA.steps_sum / init.MC_ctr)
            init.TTBTA.evs.append(init.TTBTA.evs_sum / init.MC_ctr)
## STBTA                 
        if init.STBTA.en:
            init.STBTA.utilities.append(init.STBTA.utility_sum / init.MC_ctr)
            init.STBTA.dts.append(init.STBTA.dt_sum / init.MC_ctr)
            init.STBTA.steps.append(init.STBTA.steps_sum / init.MC_ctr)
            init.STBTA.evs.append(init.STBTA.evs_sum / init.MC_ctr)            
    # end for 
    
    # Print final progress status.
    pb.printProgressBar(item_ctr, total_item, start_time, bar_color)
    
    # Save settings and algorithms' output.
    if init.monotonicity:
        store_path = os.path.abspath("output") + "/tradeoff/epsilon/monotone"
    else:
        store_path = os.path.abspath("output") + "/tradeoff/epsilon/non_monotone"
    store.storeData(store_path)
    
    
# =============================================================================
# 
# =============================================================================
print("----- tradeoff_Threshold.py is loaded -----")

    