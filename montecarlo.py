#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 16 10:57:07 2020

Monte-Carlo Simulation
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
# GA:       Genetic Algorithm
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
# auction_xx: auction based algorithms

import funcs.algs.ga as ga         # TA alg: GA

import funcs.algs.auction_sequential as auction_sequential         # TA alg: Auction_Sequential
import funcs.algs.auction_parallel as auction_parallel         # TA alg: Auction_Parallel
import funcs.algs.auction_g_prim as auction_g_prim         # TA alg: Auction_G_Prim

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


def MonteCarlo():
    
    Tasks = [j for j in range(init.Nt)]    
    
    init.Dist_mat = [[0 for j in range(init.Nt)] for i in range(init.Nt)]
    
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

# GA        
        if init.GA.en:
            init.GA.resetSum()

# Auction_Sequential        
        if init.AuctionSequential.en:
            init.AuctionSequential.resetSum()            

# Auction_Parallel        
        if init.AuctionParallel.en:
            init.AuctionParallel.resetSum()   

# Auction_G_Prim        
        if init.AuctionGPrim.en:
            init.AuctionGPrim.resetSum()              
            
            
# SGA        
        if init.SGA.en:
            init.SGA.resetSum()
# LSGA        
        if init.LSGA.en:
            init.LSGA.resetSum()
# CBBA        
        if init.CBBA.en:
            init.CBBA.resetSum()
# TGTA                 
        if init.TGTA.en:
            init.TGTA.resetSum()
# LTGTA                 
        if init.LTGTA.en:
            init.LTGTA.resetSum()
# DTTA                        
        if init.DTTA.en:
            init.DTTA.resetSum()
# LDTTA                           
        if init.LDTTA.en:
            init.LDTTA.resetSum()
# T3A                 
        if init.T3A.en:
            init.T3A.resetSum()
# LT3A                 
        if init.LT3A.en:
            init.LT3A.resetSum()
# DSTA            
        if init.DSTA.en:
            init.DSTA.resetSum()
# LSTA                     
        if init.LSTA.en:
            init.LSTA.resetSum()
# STTA                 
        if init.STTA.en:
            init.STTA.resetSum()
# LSTTA                             
        if init.LSTTA.en:
            init.LSTTA.resetSum()
# TBTA                 
        if init.TBTA.en:
            init.TBTA.resetSum()
# TTBTA                 
        if init.TTBTA.en:
            init.TTBTA.resetSum()
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
            # Na is not in settings.py, unlike Nt
#            Agents_pos = [
#                    [round(x,3) for x in np.random.uniform(low=0.0, high=init.L, size=Na)], 
#                    [round(y,3) for y in np.random.uniform(low=0.0, high=init.L, size=Na)]]       
            
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
                            for a in range(Na)
                         ]
            # non-monotone
            else: 
                # special tasks
                v1 = [round(v,3) for v in 
                      np.random.uniform(low=init.vs_low, 
                                        high=init.vs_high, 
                                        size=Na)
                      ]
                # normal tasks
                v2 = [round(v,3) for v in 
                      np.random.uniform(low=init.v_low, 
                                        high=init.v_high, 
                                        size=init.Nt-Na)
                      ]
                init.V = v1 + v2
                # special tasks
                m1 = init.ms_low*np.ones(Na) \
                        + (init.ms_high-init.ms_low)*np.eye(Na)
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
            
            # Run TA algorithms.
## GA  
            if init.GA.en:    
                selected, values, total_value, dt, steps, evs \
                    = ga.runGA(Agents,Tasks,init.pop_size, init.max_ite, init.Pr_mu)
                init.GA.utility_sum += total_value
                init.GA.dt_sum += dt
                init.GA.steps_sum += steps
                init.GA.evs_sum += evs 
## Auction_Sequential  
            if init.AuctionSequential.en:    
                selected, values, total_value, dt, steps, evs \
                    = auction_sequential.runAuctionSequential(Agents,Tasks)
                init.AuctionSequential.utility_sum += total_value
                init.AuctionSequential.dt_sum += dt
                init.AuctionSequential.steps_sum += steps
                init.AuctionSequential.evs_sum += evs    
## Auction_Parallel  
            if init.AuctionParallel.en:    
                selected, values, total_value, dt, steps, evs \
                    = auction_parallel.runAuctionParallel(Agents,Tasks)
                init.AuctionParallel.utility_sum += total_value
                init.AuctionParallel.dt_sum += dt
                init.AuctionParallel.steps_sum += steps
                init.AuctionParallel.evs_sum += evs
## Auction_G_Prim  
            if init.AuctionGPrim.en:    
                selected, values, total_value, dt, steps, evs \
                    = auction_g_prim.runAuctionGPrim(Agents,Tasks)
                init.AuctionGPrim.utility_sum += total_value
                init.AuctionGPrim.dt_sum += dt
                init.AuctionGPrim.steps_sum += steps
                init.AuctionGPrim.evs_sum += evs                
## SGA  
            if init.SGA.en:    
                selected, values, total_value, dt, steps, evs \
                    = sga.runSGA(Agents,Tasks)
                init.SGA.utility_sum += total_value
                init.SGA.dt_sum += dt
                init.SGA.steps_sum += steps
                init.SGA.evs_sum += evs                
## LSGA  
            if init.LSGA.en:    
                selected, values, total_value, dt, steps, evs \
                    = lsga.runLSGA(Agents,Tasks)
                init.LSGA.utility_sum += total_value
                init.LSGA.dt_sum += dt
                init.LSGA.steps_sum += steps
                init.LSGA.evs_sum += evs 	           
## CBBA 
            if init.CBBA.en:    
                selected, values, total_value, dt, steps, evs \
                    = cbba.runCBBA(Agents,Tasks)
                init.CBBA.utility_sum += total_value
                init.CBBA.dt_sum += dt
                init.CBBA.steps_sum += steps
                init.CBBA.evs_sum += evs 
## TGTA                   
            if init.TGTA.en:    
                selected, values, total_value, dt, steps, evs \
                    = tgta.runTGTA(Agents,Tasks,init.rho)
                init.TGTA.utility_sum += total_value
                init.TGTA.dt_sum += dt
                init.TGTA.steps_sum += steps
                init.TGTA.evs_sum += evs 
## LTGTA                   
            if init.LTGTA.en:    
                selected, values, total_value, dt, steps, evs \
                    = ltgta.runLTGTA(Agents,Tasks,init.rho)
                init.LTGTA.utility_sum += total_value
                init.LTGTA.dt_sum += dt
                init.LTGTA.steps_sum += steps
                init.LTGTA.evs_sum += evs 
## DTTA                
            if init.DTTA.en:    
                selected, values, total_value, dt, steps, evs \
                    = dtta.runDTTA(Agents,Tasks,init.eps)                    
                init.DTTA.utility_sum += total_value
                init.DTTA.dt_sum += dt
                init.DTTA.steps_sum += steps
                init.DTTA.evs_sum += evs 
## LDTTA               
            if init.LDTTA.en:    
                selected, values, total_value, dt, steps, evs \
                    = ldtta.runLDTTA(Agents,Tasks,init.eps)
                init.LDTTA.utility_sum += total_value
                init.LDTTA.dt_sum += dt
                init.LDTTA.steps_sum += steps
                init.LDTTA.evs_sum += evs 
## T3A                   
            if init.T3A.en:    
                selected, values, total_value, dt, steps, evs \
                    = t3a.runT3A(Agents,Tasks,init.rho,init.eps)
                init.T3A.utility_sum += total_value
                init.T3A.dt_sum += dt
                init.T3A.steps_sum += steps
                init.T3A.evs_sum += evs 
## LT3A                   
            if init.LT3A.en:    
                selected, values, total_value, dt, steps, evs \
                    = lt3a.runLT3A(Agents,Tasks,init.rho,init.eps)
                init.LT3A.utility_sum += total_value
                init.LT3A.dt_sum += dt
                init.LT3A.steps_sum += steps
                init.LT3A.evs_sum += evs 
## DSTA       
            if init.DSTA.en:    
                selected, values, total_value, dt, steps, evs \
                    = dsta.runDSTA(Agents,Tasks,init.Pr)
                init.DSTA.utility_sum += total_value
                init.DSTA.dt_sum += dt
                init.DSTA.steps_sum += steps
                init.DSTA.evs_sum += evs                 
## LSTA                           
            if init.LSTA.en:    
                selected, values, total_value, dt, steps, evs \
                    = lsta.runLSTA(Agents,Tasks,init.Pr)
                init.LSTA.utility_sum += total_value
                init.LSTA.dt_sum += dt
                init.LSTA.steps_sum += steps
                init.LSTA.evs_sum += evs                 
## STTA                   
            if init.STTA.en:    
                selected, values, total_value, dt, steps, evs \
                    = stta.runSTTA(Agents,Tasks,init.Pr,init.eps)
                init.STTA.utility_sum += total_value
                init.STTA.dt_sum += dt
                init.STTA.steps_sum += steps
                init.STTA.evs_sum += evs                
## LSTTA                             
            if init.LSTTA.en:    
                selected, values, total_value, dt, steps, evs \
                    = lstta.runLSTTA(Agents,Tasks,init.Pr,init.eps)
                init.LSTTA.utility_sum += total_value
                init.LSTTA.dt_sum += dt
                init.LSTTA.steps_sum += steps
                init.LSTTA.evs_sum += evs                
## TBTA                   
            if init.TBTA.en:    
                selected, values, total_value, dt, steps, evs \
                    = tbta.runTBTA(Agents,Tasks,init.eps)
                init.TBTA.utility_sum += total_value
                init.TBTA.dt_sum += dt
                init.TBTA.steps_sum += steps
                init.TBTA.evs_sum += evs 
## TTBTA                   
            if init.TTBTA.en:    
                selected, values, total_value, dt, steps, evs \
                    = ttbta.runTTBTA(Agents,Tasks,init.rho,init.eps)
                init.TTBTA.utility_sum += total_value
                init.TTBTA.dt_sum += dt
                init.TTBTA.steps_sum += steps
                init.TTBTA.evs_sum += evs 
## STBTA                   
            if init.STBTA.en:    
                selected, values, total_value, dt, steps, evs \
                    = stbta.runSTBTA(Agents,Tasks,init.Pr,init.eps)
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
## GA  
        if init.GA.en:
            init.GA.utilities.append(init.GA.utility_sum / init.MC_ctr)
            init.GA.dts.append(init.GA.dt_sum / init.MC_ctr)
            init.GA.steps.append(init.GA.steps_sum / init.MC_ctr)
            init.GA.evs.append(init.GA.evs_sum / init.MC_ctr)
## Auction_Sequential  
        if init.AuctionSequential.en:
            init.AuctionSequential.utilities.append(init.AuctionSequential.utility_sum / init.MC_ctr)
            init.AuctionSequential.dts.append(init.AuctionSequential.dt_sum / init.MC_ctr)
            init.AuctionSequential.steps.append(init.AuctionSequential.steps_sum / init.MC_ctr)
            init.AuctionSequential.evs.append(init.AuctionSequential.evs_sum / init.MC_ctr)            
## Auction_Parallel  
        if init.AuctionParallel.en:
            init.AuctionParallel.utilities.append(init.AuctionParallel.utility_sum / init.MC_ctr)
            init.AuctionParallel.dts.append(init.AuctionParallel.dt_sum / init.MC_ctr)
            init.AuctionParallel.steps.append(init.AuctionParallel.steps_sum / init.MC_ctr)
            init.AuctionParallel.evs.append(init.AuctionParallel.evs_sum / init.MC_ctr) 
## Auction_G_Prim  
        if init.AuctionGPrim.en:
            init.AuctionGPrim.utilities.append(init.AuctionGPrim.utility_sum / init.MC_ctr)
            init.AuctionGPrim.dts.append(init.AuctionGPrim.dt_sum / init.MC_ctr)
            init.AuctionGPrim.steps.append(init.AuctionGPrim.steps_sum / init.MC_ctr)
            init.AuctionGPrim.evs.append(init.AuctionGPrim.evs_sum / init.MC_ctr) 
## SGA  
        if init.SGA.en:
            init.SGA.utilities.append(init.SGA.utility_sum / init.MC_ctr)
            init.SGA.dts.append(init.SGA.dt_sum / init.MC_ctr)
            init.SGA.steps.append(init.SGA.steps_sum / init.MC_ctr)
            init.SGA.evs.append(init.SGA.evs_sum / init.MC_ctr)
## LSGA  
        if init.LSGA.en:
            init.LSGA.utilities.append(init.LSGA.utility_sum / init.MC_ctr)
            init.LSGA.dts.append(init.LSGA.dt_sum / init.MC_ctr)
            init.LSGA.steps.append(init.LSGA.steps_sum / init.MC_ctr)
            init.LSGA.evs.append(init.LSGA.evs_sum / init.MC_ctr)
## CBBA  
        if init.CBBA.en:
            init.CBBA.utilities.append(init.CBBA.utility_sum / init.MC_ctr)
            init.CBBA.dts.append(init.CBBA.dt_sum / init.MC_ctr)
            init.CBBA.steps.append(init.CBBA.steps_sum / init.MC_ctr)
            init.CBBA.evs.append(init.CBBA.evs_sum / init.MC_ctr)
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
## DSTA       
        if init.DSTA.en:
            init.DSTA.utilities.append(init.DSTA.utility_sum / init.MC_ctr)
            init.DSTA.dts.append(init.DSTA.dt_sum / init.MC_ctr)
            init.DSTA.steps.append(init.DSTA.steps_sum / init.MC_ctr)
            init.DSTA.evs.append(init.DSTA.evs_sum / init.MC_ctr)
## LSTA                  
        if init.LSTA.en:
            init.LSTA.utilities.append(init.LSTA.utility_sum / init.MC_ctr)
            init.LSTA.dts.append(init.LSTA.dt_sum / init.MC_ctr)
            init.LSTA.steps.append(init.LSTA.steps_sum / init.MC_ctr)
            init.LSTA.evs.append(init.LSTA.evs_sum / init.MC_ctr)
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
        store_path = os.path.abspath("output") + "/montecarlo/monotone"
    else:
        store_path = os.path.abspath("output") + "/montecarlo/non_monotone"
    store.storeData(store_path)
    
    
# =============================================================================
# 
# =============================================================================
print("----- montecarlo.py is loaded -----")

    