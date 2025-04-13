#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 18 15:58:04 2020

This Python project is corresponding to the PhD thesis:
    
    Title: Efficient Decentralised Task Allocation for Multiple Aerial Robots
    
    Academic year: 2017 -- 2020
    
    Candidate: Teng Li
    
    Supervisors: Prof. Hyo-Sang Shin
                 Prof. Antonios Tsourdos
    
    Affiliation: School of Aerospace, Transport and Manufacturing (SATM)
                 Cranfield University

Language:
    Python 3.12



** non-commercial use only **
@author: Teng Li
lt.uk@outlook.com
United Kingdom
All Rights Reserved
"""

# ==========================  Load Supporting Modules  ========================
# Do NOT change the import order.
import init                          # algorithms initialisation
import funcs.sys.progressbar as pb   # print progress bar
import montecarlo as mc
import montecarlo_change_nt as mc_nt
import variance as vr    # Run stochastic algs for multiple rounds, deterministic algs for one round
import tradeoff_Pr
import tradeoff_Threshold
import tradeoff_Threshold_vr  # For stochastic algs: STTA, LSTTA, STBTA 
import tradeoff_Truncation
# =============================================================================

# =====================  Load Plot Modules  ===================================
import funcs.plots.plotall as plotall        # plot figures for all algorithms
import funcs.plots.plot_dsta_mc as plot_dsta_mc      # plot figures for DSTA vs CBBA
import funcs.plots.plot_t3a_mc  as plot_t3a_mc       # plot figures for T3A vs benchmark
import funcs.plots.plot_dtta_mc as plot_dtta_mc      # plot figures for DTTA vs benchmark
import funcs.plots.plot_stta_mc as plot_stta_mc      # plot figures for STTA vs benchmark
import funcs.plots.plot_dsta_vr as plot_dsta_vr      # plot figures for DSTA variance vs CBBA
import funcs.plots.plot_stta_vr as plot_stta_vr      # plot variance figures for STTA vs benchmark 
import funcs.plots.plot_stta_vs_t3a as plot_stta_vs_t3a
import funcs.plots.plot_stta_vs_t3a_vr as plot_stta_vs_t3a_vr
import funcs.plots.plot_stta_vs_dsta_vr as plot_stta_vs_dsta_vr
import funcs.plots.plot_tradeoff_Pr as plot_tradeoff_Pr
import funcs.plots.plot_tradeoff_Threshold as plot_tradeoff_Threshold
import funcs.plots.plot_tradeoff_Threshold_vr as plot_tradeoff_Threshold_vr
import funcs.plots.plot_tradeoff_Truncation as plot_tradeoff_Truncation

import funcs.plots.plotall_nt as plotall_nt        # plot figures for all algorithms with changing Nt
import funcs.plots.plot_dtta_mc_nt as plot_dtta_mc_nt      # plot figures for DTTA vs benchmark with changing Nt
# =============================================================================


# ====================== Set TA parameters between here =======================

init.Nt = 300    # number of tasks 200, 60, 50   200 tasks too slow
init.evs_scale = 3  # scale in plot figures, scale = 10**evs_scale

init.Na = 100 # 12 or 4, fixed number of agents for tradeoffs 15 for DSTA

init.L = 10.0       # /km  The length of the 2-D space.
init.ro = 20.0  # km  The radius of the operation range of satellites in Satellite surveiliance missions

# value factors of the normal tasks
init.v_low = 0.6    # 0.6
init.v_high = 1.0   # 1.0
# value factors of special tasks in the non-monotone case
init.vs_low = 5   # 5
init.vs_high = 6  # 6

# match fitness factors for normal tasks
init.m_low = 0.5    # 0.5
init.m_high = 1.0   # 1.0
# match fitness factors for special tasks in the non-monotone case
init.ms_low = 0.1   # 0.1
init.ms_high = 0.2  # 0.2

init.d_0 = 1        # /km characteristic distance

init.lambda_x = 0.01     # penalty scaling factor in the non-monotone case

init.Pr = 0.5       # sampling probability for all sample related algorithms
init.eps = 0.1      # epsilon, decreasing threshold parameter
init.rho = 1        # truncation parameter

# GA parameters
init.max_ite = 500 # The maximum number of iterations.
init.pop_size = 20  # Population size.
init.Pr_mu = 0.5   # The probability of mutation.

init.MC_ctr = 20     # The number of Monte Carlo loops.

init.Na_start = 100   # 10, 4
init.Na_end = 101    # 51, 21
init.Na_step = 2    # 5,  2

init.Nt_start = 30   # 
init.Nt_end = 71    # 
init.Nt_step = 5    # 

init.Pr_start = 0.1
init.Pr_end = 0.9
init.Pr_step = 0.1

init.Eps_start = 0.1
init.Eps_end = 0.5
init.Eps_step = 0.05

init.Rho_start = 0
init.Rho_end = 1
init.Rho_step = 0.1

# ------------------ Alg Enabler ------------------
# Algorithm enable flags  1:Enable  0:Disable
init.GA.en = 0

init.AuctionSequential.en = 0
init.AuctionParallel.en = 0
init.AuctionGPrim.en = 0


init.SGA.en = 1
init.LSGA.en = 0

init.CBBA.en = 0

init.TGTA.en = 0
init.LTGTA.en = 0

init.DTTA.en = 1
init.LDTTA.en = 0
init.TBTA.en = 0

init.T3A.en = 0
init.LT3A.en = 0
init.TTBTA.en = 0

init.DSTA.en = 1
init.LSTA.en = 0

init.STTA.en = 0
init.LSTTA.en = 0
init.STBTA.en = 0
# ---------------- End Enabler ---------------------

#init.scenario = 1

init.monotonicity = 1    # 1: monotone 0: non-monotone

# Do NOT enable any two of the followings at the same time.
montecarlo_en = 1

montecarlo_nt_en = 0

variance_en = 0

tradeoff_Pr_en = 0

tradeoff_Eps_en = 0

tradeoff_Eps_Vr_en = 0    # Variance

tradeoff_Rho_en = 0
# ================================= and here ==================================


if __name__ == "__main__":
    
    print("================" + 
          pb.colors.VIOLET + ' main ' + pb.colors.END + 
          "================")
    
    # Monte-Carlo Simulation
    if montecarlo_en:
        mc.MonteCarlo()
        # Plot all algorithms. 1: plot; 0: no plot
        plotall.plotAll(benchmark=1, original=1, bundle=1, lazy=1)
        
        # plot_stta_vs_t3a.plotSTTA_vs_T3A()
        # plot_dtta_mc.plotDTTA_MC()
        plot_dtta_mc.plotDTTA_MC_journal_lg()
        # plot_dtta_mc.plotDTTA_MC_journal()
#        plot_t3a_mc.plotT3A_MC()
#        plot_stta_mc.plotSTTA_MC()
    
#        plot_dsta_mc.plotDSTA_MC()
    
    # Monte-Carlo Simulation -- Change Nt
    if montecarlo_nt_en:
        mc_nt.MonteCarlo_ChangeNt()
        
        plotall_nt.plotAllNt()
        plot_dtta_mc_nt.plotDTTA_MC_Nt(col_index=0)
        plot_dtta_mc_nt.plotDTTA_MC_Nt(col_index=-1)

    # Simulation for showing sampling variance
    if variance_en:
        vr.Variance()
        
#        plot_dsta_vr.plotDSTA_Vr()
#        plot_dsta_vr.plotDSTA_Vr_Box()
#        plot_stta_vr.plotSTTA_Vr()
        plot_stta_vr.plotSTTA_Vr_Box()
        
#        plot_stta_vs_t3a_vr.plotSTTA_vs_T3A_Vr()
#        plot_stta_vs_dsta_vr.plotSTTA_DSTA_Vr()
#        plot_stta_vs_dsta_vr.plotSTTA_DSTA_Vr_Box()
    
    if tradeoff_Pr_en:
        tradeoff_Pr.runTradeOff_Pr()
        
#        plot_tradeoff_Pr.plotTradeoff_Pr() # all
#        plot_tradeoff_Pr.plotTradeoff_Pr_Box() # DSTA
        plot_tradeoff_Pr.plotTradeoff_Pr_Box_STTA()
        
    if tradeoff_Eps_en:
        tradeoff_Threshold.runTradeoff_Eps()
        plot_tradeoff_Threshold.plotTradeoff_Eps()
        
    if tradeoff_Eps_Vr_en:
        tradeoff_Threshold_vr.runTradeoff_Eps_Vr()
#        plot_tradeoff_Threshold_vr.plotTradeoff_Eps_Vr()
        plot_tradeoff_Threshold_vr.plotTradeoff_Eps_Box_STTA()
    
    if tradeoff_Rho_en:
        tradeoff_Truncation.runTradeoff_Rho()
        plot_tradeoff_Truncation.plotTradeoff_Rho()
        
    # ----------------- Log data -----------------
    
    # GA
    if init.GA.en:
        GA_steps = init.GA.steps
        GA_evs = init.GA.evs
        GA_utilities = init.GA.utilities
        
    # Auctions
    if init.AuctionSequential.en:
        AuctionSequential_utilities = init.AuctionSequential.utilities
        AuctionSequential_evs = init.AuctionSequential.evs
        AuctionSequential_steps = init.AuctionSequential.steps
        
    if init.AuctionParallel.en:
        AuctionParallel_steps = init.AuctionParallel.steps
        AuctionParallel_evs = init.AuctionParallel.evs
        AuctionParallel_utilities = init.AuctionParallel.utilities
        
    if init.AuctionGPrim.en:    
        AuctionGPrim_utilities = init.AuctionGPrim.utilities
        AuctionGPrim_evs = init.AuctionGPrim.evs
        AuctionGPrim_steps = init.AuctionGPrim.steps
    
    # SGA
    if init.SGA.en:
        SGA_steps = init.SGA.steps
        SGA_evs = init.SGA.evs
        SGA_utilities = init.SGA.utilities
        
    if init.LSGA.en:
        LSGA_steps = init.LSGA.steps
        LSGA_evs = init.LSGA.evs
        LSGA_utilities = init.LSGA.utilities
        
    
    # CBBA
    if init.CBBA.en:
        CBBA_steps = init.CBBA.steps
        CBBA_evs = init.CBBA.evs
        CBBA_utilities = init.CBBA.utilities
        
    
    # TGTA
    if init.TGTA.en:
        TGTA_steps = init.TGTA.steps
        TGTA_evs = init.TGTA.evs
        TGTA_utilities = init.TGTA.utilities
        
    if init.LTGTA.en:
        LTGTA_steps = init.LTGTA.steps
        LTGTA_evs = init.LTGTA.evs
        LTGTA_utilities = init.LTGTA.utilities
        
    
    # DTTA
    if init.DTTA.en:
        DTTA_steps = init.DTTA.steps
        DTTA_evs = init.DTTA.evs
        DTTA_utilities = init.DTTA.utilities
        
    if init.LDTTA.en:
        LDTTA_steps = init.LDTTA.steps
        LDTTA_evs = init.LDTTA.evs
        LDTTA_utilities = init.LDTTA.utilities   
        
    if init.TBTA.en:
        TBTA_steps = init.TBTA.steps
        TBTA_evs = init.TBTA.evs
        TBTA_utilities = init.TBTA.utilities
    
    # T3A
    if init.T3A.en:
        T3A_steps = init.T3A.steps
        T3A_evs = init.T3A.evs
        T3A_utilities = init.T3A.utilities
        
    if init.LT3A.en:
        LT3A_steps = init.LT3A.steps
        LT3A_evs = init.LT3A.evs
        LT3A_utilities = init.LT3A.utilities   
        
    if init.TBTA.en:
        TBTA_steps = init.TBTA.steps
        TBTA_evs = init.TBTA.evs
        TBTA_utilities = init.TBTA.utilities
    
    
    # DSTA
    if init.DSTA.en:
        DSTA_steps = init.DSTA.steps
        DSTA_evs = init.DSTA.evs
        DSTA_utilities = init.DSTA.utilities
    
    if init.LSTA.en:
        LSTA_steps = init.LSTA.steps
        LSTA_evs = init.LSTA.evs
        LSTA_utilities = init.LSTA.utilities


    # STTA
    if init.STTA.en:
        STTA_steps = init.STTA.steps
        STTA_evs = init.STTA.evs
        STTA_utilities = init.STTA.utilities
        
    if init.LSTTA.en:
        LSTTA_steps = init.LSTTA.steps
        LSTTA_evs = init.LSTTA.evs
        LSTTA_utilities = init.LSTTA.utilities   
        
    if init.STBTA.en:
        STBTA_steps = init.STBTA.steps
        STBTA_evs = init.STBTA.evs
        STBTA_utilities = init.STBTA.utilities
        
        
        
        
        