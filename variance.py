#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 16 15:29:42 2020

Show variance resulted from random sampling.
Random task positions, task values, task-agent match fitness factors etc. are 
fixed.

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
# SGA:      Sequencial Greedy Algorithm
# CBBA:     Consensus Based Bundle Algorithm
# TGTA:     Truncation Greedy Task Allocation
# DTTA:     Decreasing Threshold Task Allocation
# TBTA:     Threshold Bundle Task Allocation
# T3A:      Truncation Threshold Task Allocation
# TTBTA:    Truncation Threshold Bundle Task Allocation
# DSTA:     Decentralised Sample based Task Allocation
# STTA:     Sample Threshold Task Allocation
# STBTA:    Sample Threshold Bundle Task Allocation

import funcs.algs.sga as sga       # TA alg: SGA
import funcs.algs.lsga as lsga     # TA alg: LSGA

import funcs.algs.cbba as cbba     # TA alg: CBBA

import funcs.algs.tgta as tgta     # TA alg: TGTA
import funcs.algs.ltgta as ltgta   # TA alg: LTGTA

import funcs.algs.dtta as dtta     # TA alg: DTTA
import funcs.algs.ldtta as ldtta   # TA alg: LDTTA
import funcs.algs.tbta as tbta     # TA alg: TBTA

import funcs.algs.t3a as t3a       # TA alg: T3A
import funcs.algs.lt3a as lt3a     # TA alg: LT3A
import funcs.algs.ttbta as ttbta   # TA alg: TTBTA

import funcs.algs.dsta as dsta     # TA alg: DSTA
import funcs.algs.lsta as lsta     # TA alg: LSTA

import funcs.algs.stta as stta     # TA alg: STTA
import funcs.algs.lstta as lstta   # TA alg: LSTTA
import funcs.algs.stbta as stbta   # TA alg: STBTA
# =============================================================================


def Variance():
    
    Tasks = [j for j in range(init.Nt)] 
    
    init.Dist_mat = [[0 for j in range(init.Nt)] for i in range(init.Nt)]
    
    # Tasks' and Agents' positions (x,y) /km
    Tasks_pos = [
            [round(x,3) for x in 
             np.random.uniform(low=0.0, high=init.L, size=init.Nt)], 
            [round(y,3) for y in 
             np.random.uniform(low=0.0, high=init.L, size=init.Nt)]]     
    
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
                                    size=init.Nt)]
        # random task-agent match fitness factors
        init.M = [[round(m,3) for m in 
                   np.random.uniform(init.m_low,init.m_high,init.Nt)
                  ] for a in range(init.Na_end)]
    
    # Progress bar color options: RED, GREEN, ORANGE, BLUE, VIOLET, CYAN, GREY.
    bar_color = pb.colors.CYAN # (default: CYAN)
    start_time = time.time()
    time_last = start_time
    item_ctr = 0
    total_item = init.MC_ctr*((init.Na_end-init.Na_start)//init.Na_step+1)
    # Print initial status.
    pb.printProgressBar(item_ctr, total_item, start_time, bar_color)
    
    # Change the number of agents       
    for Na in range(init.Na_start, init.Na_end+1, init.Na_step):
        # Na: The number of Agents
        init.Nas.append(Na)
        Agents = [a for a in range(Na)] 
        
        # non-monotone case
        if not init.monotonicity:
            # special tasks math.ceil(Na/1.5)
            v1 = [round(v,3) for v in 
                  np.random.uniform(low=init.vs_low, 
                                    high=init.vs_high, 
                                    size=round(Na*1))
                 ]
            # normal tasks
            v2 = [round(v,3) for v in 
                  np.random.uniform(low=init.v_low, 
                                    high=init.v_high, 
                                    size=init.Nt-round(Na*1))
                 ]
            init.V = v1 + v2
            # special tasks
            m1 = init.ms_low*np.ones(Na) + \
                (init.ms_high-init.ms_low)*np.eye(Na)
            # normal tasks
            m2 = [[round(m,3) for m in 
                   np.random.uniform(init.m_low,init.m_high,init.Nt-Na)] 
                    for a in range(Na)
                 ]
            init.M = np.concatenate((m1, m2), axis=1) 
            # inter-task dependences Nt*Nt
            init.X = np.zeros((init.Nt, init.Nt))
            for i in range(init.Nt-1):
                for j in range(i+1, init.Nt):
                    init.X[i][j] = round(math.e**(init.V[i]*init.V[j]),4)
                    init.X[j][i] = init.X[i][j]

# Stochastic Algs
# DSTA            
        if init.DSTA.en:
            init.DSTA.resetTemp()
# LSTA                     
        if init.LSTA.en:
            init.LSTA.resetTemp()
# STTA                 
        if init.STTA.en:
            init.STTA.resetTemp()
# LSTTA                             
        if init.LSTTA.en:
            init.LSTTA.resetTemp()
# STBTA                 
        if init.STBTA.en:
            init.STBTA.resetTemp()

# Deterministic Algs
## SGA  
        if init.SGA.en:    
            selected, values, total_value, dt, steps, evs \
                = sga.runSGA(Agents,Tasks)
            init.SGA.utilities.append(total_value)
            init.SGA.dts.append(dt)
            init.SGA.steps.append(steps)
            init.SGA.evs.append(evs)              
## LSGA  
        if init.LSGA.en:    
            selected, values, total_value, dt, steps, evs \
                = lsga.runLSGA(Agents,Tasks)
            init.LSGA.utilities.append(total_value)
            init.LSGA.dts.append(dt)
            init.LSGA.steps.append(steps)
            init.LSGA.evs.append(evs)   	           
## CBBA 
        if init.CBBA.en:    
            selected, values, total_value, dt, steps, evs \
                = cbba.runCBBA(Agents,Tasks)
            init.CBBA.utilities.append(total_value)
            init.CBBA.dts.append(dt)
            init.CBBA.steps.append(steps)
            init.CBBA.evs.append(evs)   
## TGTA                   
        if init.TGTA.en:    
            selected, values, total_value, dt, steps, evs \
                = tgta.runTGTA(Agents,Tasks,init.rho)
            init.TGTA.utilities.append(total_value)
            init.TGTA.dts.append(dt)
            init.TGTA.steps.append(steps)
            init.TGTA.evs.append(evs)  
## LTGTA                   
        if init.LTGTA.en:    
            selected, values, total_value, dt, steps, evs \
                = ltgta.runLTGTA(Agents,Tasks,init.rho)
            init.LTGTA.utilities.append(total_value)
            init.LTGTA.dts.append(dt)
            init.LTGTA.steps.append(steps)
            init.LTGTA.evs.append(evs)  
## DTTA                
        if init.DTTA.en:    
            selected, values, total_value, dt, steps, evs \
                = dtta.runDTTA(Agents,Tasks,init.eps)                    
            init.DTTA.utilities.append(total_value)
            init.DTTA.dts.append(dt)
            init.DTTA.steps.append(steps)
            init.DTTA.evs.append(evs)  
## LDTTA               
        if init.LDTTA.en:    
            selected, values, total_value, dt, steps, evs \
                = ldtta.runLDTTA(Agents,Tasks,init.eps)
            init.LDTTA.utilities.append(total_value)
            init.LDTTA.dts.append(dt)
            init.LDTTA.steps.append(steps)
            init.LDTTA.evs.append(evs)  
## TBTA                   
        if init.TBTA.en:    
            selected, values, total_value, dt, steps, evs \
                = tbta.runTBTA(Agents,Tasks,init.eps)
            init.TBTA.utilities.append(total_value)
            init.TBTA.dts.append(dt)
            init.TBTA.steps.append(steps)
            init.TBTA.evs.append(evs)
## T3A                   
        if init.T3A.en:    
            selected, values, total_value, dt, steps, evs \
                = t3a.runT3A(Agents,Tasks,init.rho,init.eps)
            init.T3A.utilities.append(total_value)
            init.T3A.dts.append(dt)
            init.T3A.steps.append(steps)
            init.T3A.evs.append(evs)  
## LT3A                   
        if init.LT3A.en:    
            selected, values, total_value, dt, steps, evs \
                = lt3a.runLT3A(Agents,Tasks,init.rho,init.eps)
            init.LT3A.utilities.append(total_value)
            init.LT3A.dts.append(dt)
            init.LT3A.steps.append(steps)
            init.LT3A.evs.append(evs)  
## TTBTA                   
        if init.TTBTA.en:    
            selected, values, total_value, dt, steps, evs \
                = ttbta.runTTBTA(Agents,Tasks,init.rho,init.eps)
            init.TTBTA.utilities.append(total_value)
            init.TTBTA.dts.append(dt)
            init.TTBTA.steps.append(steps)
            init.TTBTA.evs.append(evs) 
        
        # Monte-Carlo Loop
        for counter in range(init.MC_ctr):       
            # Run stochastic TA algorithms.
## DSTA       
            if init.DSTA.en:    
                selected, values, total_value, dt, steps, evs \
                    = dsta.runDSTA(Agents,Tasks,init.Pr)
                init.DSTA.utilities_temp.append(total_value)
                init.DSTA.dts_temp.append(dt)
                init.DSTA.steps_temp.append(steps)
                init.DSTA.evs_temp.append(evs)                 
## LSTA                           
            if init.LSTA.en:    
                selected, values, total_value, dt, steps, evs \
                    = lsta.runLSTA(Agents,Tasks,init.Pr)
                init.LSTA.utilities_temp.append(total_value)
                init.LSTA.dts_temp.append(dt)
                init.LSTA.steps_temp.append(steps)
                init.LSTA.evs_temp.append(evs)                 
## STTA                   
            if init.STTA.en:    
                selected, values, total_value, dt, steps, evs \
                    = stta.runSTTA(Agents,Tasks,init.Pr,init.eps)
                init.STTA.utilities_temp.append(total_value)
                init.STTA.dts_temp.append(dt)
                init.STTA.steps_temp.append(steps)
                init.STTA.evs_temp.append(evs)               
## LSTTA                             
            if init.LSTTA.en:    
                selected, values, total_value, dt, steps, evs \
                    = lstta.runLSTTA(Agents,Tasks,init.Pr,init.eps)
                init.LSTTA.utilities_temp.append(total_value)
                init.LSTTA.dts_temp.append(dt)
                init.LSTTA.steps_temp.append(steps)
                init.LSTTA.evs_temp.append(evs)               
## STBTA                   
            if init.STBTA.en:    
                selected, values, total_value, dt, steps, evs \
                    = stbta.runSTBTA(Agents,Tasks,init.Pr,init.eps)
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

        # end for (Monte-Carlo loop)

        # Stochasitc Algs #
## DSTA       
        if init.DSTA.en:
            init.DSTA.utilities.append(init.DSTA.utilities_temp)
            init.DSTA.dts.append(init.DSTA.dts_temp)
            init.DSTA.steps.append(init.DSTA.steps_temp)
            init.DSTA.evs.append(init.DSTA.evs_temp)
## LSTA                  
        if init.LSTA.en:
            init.LSTA.utilities.append(init.LSTA.utilities_temp)
            init.LSTA.dts.append(init.LSTA.dts_temp)
            init.LSTA.steps.append(init.LSTA.steps_temp)
            init.LSTA.evs.append(init.LSTA.evs_temp)
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
    # end for (Nas)
    
    # Print final progress status.
    pb.printProgressBar(item_ctr, total_item, start_time, bar_color)
    
    # Save settings and algorithms' output.
    if init.monotonicity:
        store_path = os.path.abspath("output") + "/variance/monotone"
    else:
        store_path = os.path.abspath("output") + "/variance/non_monotone"
    store.storeData(store_path)
    
# =============================================================================
# 
# =============================================================================
print("----- variance.py is loaded -----") 
    