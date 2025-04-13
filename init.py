#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May 24 11:38:56 2020

Initialization for all algorithms.

Caution! If any change occurs in this file, "Run selection or current line" 
will also call the reloading of this file. Therefore, all parameters will be 
reset. To use the output data again, call functions restoreSettings(path) and 
restorePerformance(path) in funcs/sys/store.py.

@author: Teng Li
tengli@cranfield.ac.uk; nliteng@foxmail.com
Cranfield University, UK
All Rights Reserved
"""

import matplotlib.pyplot as plt
import statistics as st

scenario = 1


Nt = 0          # number of tasks
Lt = Nt         # The limitation of task quantity for each agent.
evs_scale = 0   # scale in plot figures, scale = 10**evs_scale

Na = 0          # Fixed number of agents in trade-off.

Nst = 0         # Number of special tasks

# lambda_d = 0.95 # discount factor for total distance of path
# lambda_n = 0.98 # discount factor for total number of tasks

L = 0.0         # km  The length of the 2-D space.

# The matrix containing distances between tasks.
Dist_mat = []

# The matrix containing distances between Agents and Tasks.
Dist_AT = []

# value factors of the tasks
V = []
# normal tasks
v_low = 0.0
v_high = 0.0
# special tasks in the non-monotone case
vs_low = 0.0    
vs_high = 0.0  


# match fitness factors
M = []
# normal tasks
m_low = 0.0
m_high = 0.0
# special tasks in the non-monotone case
ms_low = 0.0
ms_high = 0.0

d_0 = 1     # /km characteristic distance

monotonicity = 1    # 1: monotone 0: non-monotone
lambda_x = 0.01     # penalty scaling factor in the non-monotone case
X = []              # inter-task effect matrix in the non-monotone case, NtxNt

Pr = 0.5    # sampling probability for all sample related algorithms
eps = 0.1   # epsilon, decreasing threshold parameter
rho = 1       # truncation parameter

# GA parameters
max_ite = 100 # The maximum number of iterations.
pop_size = 1  # Population size.
Pr_mu = 0.1   # The probability of mutation.

MC_ctr = 1  # The number of Monte Carlo loops.

# init the agents number container
Nas = []        # contains all Na values, Na is changing number of agents
Na_start = 0
Na_end = 0
Na_step = 0

# init the task number container
Nts = []        # contains all Nt values, Nt is changing number of tasks
Nt_start = 0
Nt_end = 0
Nt_step = 0

# Used for tradeoff with sampling probability
Prs = []
Pr_start = 0
Pr_end = 0
Pr_step = 0

# Used for tradeoff with threshold decreasing parameter epsilon
Epses = []
Eps_start = 0
Eps_end = 0
Eps_step = 0

# Used for tradeoff with truncation parameter k.
Rhos = []
Rho_start = 0
Rho_end = 0
Rho_step = 0

# error bar settings
#eb_fmt = "c.-"
eb_markersize = 3
eb_capsize = 2
eb_capsize_line = 4
eb_capsize_bar = 2
#eb_label = "DSTA"

# =============================================================================
# Template for all algrothms.
# =============================================================================
class Algorithm:
    
    # initialisation
    def __init__(self, linestyle, lbl, eb_style=""):
        self.en = 0
        self.ls = linestyle
        self.utilities = []
        self.dts = []
        self.steps = []
        self.evs = []
        self.label = lbl
        self.utility_sum = 0
        self.dt_sum = 0
        self.steps_sum = 0
        self.evs_sum = 0
        self.utilities_temp = []
        self.dts_temp = []
        self.steps_temp = []
        self.evs_temp = []
        self.upper_error = []
        self.lower_error = []
        self.utilities_avg = []
        self.error_utilities = []
        self.dts_avg = []
        self.error_dts = []
        self.steps_avg = []
        self.error_steps = []
        self.evs_avg = []
        self.error_evs = []
        self.evs_avg_scale = []
        self.error_evs_scale = []
        self.eb_fmt = eb_style
    
    # reset sum variables, used for Monte-Carlo 
    def resetSum(self):
        self.utility_sum = 0
        self.dt_sum = 0
        self.steps_sum = 0
        self.evs_sum = 0
    
    # reset temparary variables, used for Variance
    def resetTemp(self):
        self.utilities_temp = []
        self.dts_temp = []
        self.steps_temp = []
        self.evs_temp = []
    
    # plot Monte-Carlo figures for all algs 
    # -- random tasks and agents.
    def plotUtilities(self):
        if self.en:
            plt.plot(Nas, self.utilities, self.ls, label = self.label) 
    def plotDts(self):
        if self.en:
            plt.plot(Nas, self.dts, self.ls, label = self.label)
    def plotSteps(self):
        if self.en:
            plt.plot(Nas, self.steps, self.ls, label = self.label)
    def plotEvs(self):
        if self.en:
            plt.plot(Nas, [evs_temp/10**evs_scale for evs_temp in self.evs], 
                     self.ls, label = self.label)
            
    def plotUtilities_Nt(self):
        if self.en:
            plt.plot(Nts, self.utilities, self.ls, label = self.label) 
    def plotDts_Nt(self):
        if self.en:
            plt.plot(Nts, self.dts, self.ls, label = self.label)
    def plotSteps_Nt(self):
        if self.en:
            plt.plot(Nts, self.steps, self.ls, label = self.label)
    def plotEvs_Nt(self):
        if self.en:
            plt.plot(Nts, [evs_temp/10**evs_scale for evs_temp in self.evs], 
                     self.ls, label = self.label)
    
    # plot Monte-Carlo figures with errorbar for stochastic algs 
    # -- fixed tasks and agents.
    def plotVariance_Utilities(self):
        if self.en:
            self.utilities_avg = [st.mean(self.utilities[i]) 
                                    for i in range(len(Nas))]
            self.upper_error = [max(self.utilities[i])-self.utilities_avg[i] 
                            for i in range(len(Nas))]
            self.lower_error = [-min(self.utilities[i])+self.utilities_avg[i] 
                            for i in range(len(Nas))]
            self.error_utilities = [self.lower_error, self.upper_error]
            plt.errorbar(Nas, self.utilities_avg,
                         yerr=self.error_utilities,
                         fmt=self.eb_fmt,
                         capsize=eb_capsize,
                         label = self.label)
            max_height = max([max(self.utilities[i]) for i in range(len(Nas))])
            plt.ylim(ymin=0, ymax=max_height*1.1)
    
    def plotVariance_Dts(self):
        if self.en:
            self.dts_avg = [st.mean(self.dts[i]) 
                                    for i in range(len(Nas))]
            self.upper_error = [max(self.dts[i])-self.dts_avg[i] 
                            for i in range(len(Nas))]
            self.lower_error = [-min(self.dts[i])+self.dts_avg[i] 
                            for i in range(len(Nas))]
            self.error_dts = [self.lower_error, self.upper_error]
            plt.errorbar(Nas, self.dts_avg,
                         yerr=self.error_dts,
                         fmt=self.eb_fmt,
                         capsize=eb_capsize,
                         label = self.label)
#            max_height = max([max(self.dts[i]) for i in range(len(Nas))])
#            plt.ylim(ymin=0, ymax=max_height*1.1)
    
    def plotVariance_Steps(self):
        if self.en:
            self.steps_avg = [st.mean(self.steps[i]) 
                                    for i in range(len(Nas))]
            self.upper_error = [max(self.steps[i])-self.steps_avg[i] 
                            for i in range(len(Nas))]
            self.lower_error = [-min(self.steps[i])+self.steps_avg[i] 
                            for i in range(len(Nas))]
            self.error_steps = [self.lower_error, self.upper_error]
            plt.errorbar(Nas, self.steps_avg,
                         yerr=self.error_steps,
                         fmt=self.eb_fmt,
                         capsize=eb_capsize,
                         label = self.label)
#            max_height = max([max(self.steps[i]) for i in range(len(Nas))])
#            plt.ylim(ymin=0, ymax=max_height*1.1)
    
    def plotVariance_Evs(self):
        if self.en:
            self.evs_avg = [st.mean(self.evs[i]) 
                                    for i in range(len(Nas))]
            self.upper_error = [max(self.evs[i])-self.evs_avg[i] 
                            for i in range(len(Nas))]
            self.lower_error = [-min(self.evs[i])+self.evs_avg[i] 
                            for i in range(len(Nas))]
            self.error_evs = [self.lower_error, self.upper_error]
            self.evs_avg_scale = [evs/10**evs_scale for evs in self.evs_avg]
            self.error_evs_scale = [[err/10**evs_scale for err in self.lower_error], 
                                    [err/10**evs_scale for err in self.upper_error]]
            plt.errorbar(Nas, self.evs_avg_scale,
                         yerr=self.error_evs_scale,
                         fmt=self.eb_fmt,
                         capsize=eb_capsize,
                         label = self.label)
#            max_height = max([max(self.evs[i])/10**evs_scale for i in range(len(Nas))])
#            plt.ylim(ymin=0, ymax=max_height*1.1)
    
    # Plot tradeoff of the decreasing threshold figures.
    def plotTradeoff_Eps_Utilities(self):
        if self.en:
            plt.plot(Epses, self.utilities, self.ls, label = self.label)
    def plotTradeoff_Eps_Dts(self):
        if self.en:
            plt.plot(Epses, self.dts, self.ls, label = self.label)
    def plotTradeoff_Eps_Steps(self):
        if self.en:
            plt.plot(Epses, self.steps, self.ls, label = self.label)
    def plotTradeoff_Eps_Evs(self):
        if self.en:
            plt.plot(Epses, [evs_temp/10**evs_scale for evs_temp in self.evs], 
                     self.ls, label = self.label)
            
    # Plot tradeoff of the truncation position figures.
    def plotTradeoff_Rho_Utilities(self):
        if self.en:
            plt.plot(Rhos, self.utilities, self.ls, label = self.label)
    def plotTradeoff_Rho_Dts(self):
        if self.en:
            plt.plot(Rhos, self.dts, self.ls, label = self.label)
    def plotTradeoff_Rho_Steps(self):
        if self.en:
            plt.plot(Rhos, self.steps, self.ls, label = self.label)
    def plotTradeoff_Rho_Evs(self):
        if self.en:
            plt.plot(Rhos, [evs_temp/10**evs_scale for evs_temp in self.evs], 
                     self.ls, label = self.label)
    
    # plot tradeoff of probability figures with errorbar.
    def plotTradeoff_Pr_Utilities(self):
        if self.en:
            self.utilities_avg = [st.mean(self.utilities[i]) 
                                    for i in range(len(Prs))]
            self.upper_error = [max(self.utilities[i])-self.utilities_avg[i] 
                            for i in range(len(Prs))]
            self.lower_error = [-min(self.utilities[i])+self.utilities_avg[i] 
                            for i in range(len(Prs))]
            self.error_utilities = [self.lower_error, self.upper_error]
            plt.errorbar(Prs, self.utilities_avg,
                         yerr=self.error_utilities,
                         fmt=self.ls,
                         capsize=eb_capsize,
                         label = self.label)
            max_height = max([max(self.utilities[i]) for i in range(len(Prs))])
            plt.ylim(ymin=0, ymax=max_height*1.1)
    
    def plotTradeoff_Pr_Dts(self):        
        if self.en:
            self.dts_avg = [st.mean(self.dts[i]) 
                                    for i in range(len(Prs))]
            self.upper_error = [max(self.dts[i])-self.dts_avg[i] 
                            for i in range(len(Prs))]
            self.lower_error = [-min(self.dts[i])+self.dts_avg[i] 
                            for i in range(len(Prs))]
            self.error_dts = [self.lower_error, self.upper_error]
            plt.errorbar(Prs, self.dts_avg,
                         yerr=self.error_dts,
                         fmt=self.ls,
                         capsize=eb_capsize,
                         label = self.label)
            max_height = max([max(self.dts[i]) for i in range(len(Prs))])
            plt.ylim(ymin=0, ymax=max_height*1.1)
    
    def plotTradeoff_Pr_Steps(self): 
        if self.en:
            self.steps_avg = [st.mean(self.steps[i]) 
                                    for i in range(len(Prs))]
            self.upper_error = [max(self.steps[i])-self.steps_avg[i] 
                            for i in range(len(Prs))]
            self.lower_error = [-min(self.steps[i])+self.steps_avg[i] 
                            for i in range(len(Prs))]
            self.error_steps = [self.lower_error, self.upper_error]
            plt.errorbar(Prs, self.steps_avg,
                         yerr=self.error_steps,
                         fmt=self.ls,
                         capsize=eb_capsize,
                         label = self.label)
            max_height = max([max(self.steps[i]) for i in range(len(Prs))])
            plt.ylim(ymin=0, ymax=max_height*1.1)
        
    def plotTradeoff_Pr_Evs(self): 
        if self.en:
            self.evs_avg = [st.mean(self.evs[i]) 
                                    for i in range(len(Prs))]
            self.upper_error = [max(self.evs[i])-self.evs_avg[i] 
                            for i in range(len(Prs))]
            self.lower_error = [-min(self.evs[i])+self.evs_avg[i] 
                            for i in range(len(Prs))]
            self.error_evs = [self.lower_error, self.upper_error]
            self.evs_avg_scale = [evs/10**evs_scale for evs in self.evs_avg]
            self.error_evs_scale = [[err/10**evs_scale for err in self.lower_error], 
                                    [err/10**evs_scale for err in self.upper_error]]
            plt.errorbar(Prs, self.evs_avg_scale,
                         yerr=self.error_evs_scale,
                         fmt=self.ls,
                         capsize=eb_capsize,
                         label = self.label)
#            max_height = max([max(self.evs[i])/10**evs_scale for i in range(len(Prs))])
#            plt.ylim(ymin=0, ymax=max_height*1.1)
    
    # plot tradeoff of epsilon figures with errorbar for stochastic algorithms.
    def plotTradeoff_Eps_Vr_Utilities(self):
        if self.en:
            self.utilities_avg = [st.mean(self.utilities[i]) 
                                    for i in range(len(Epses))]
            self.upper_error = [max(self.utilities[i])-self.utilities_avg[i] 
                            for i in range(len(Epses))]
            self.lower_error = [-min(self.utilities[i])+self.utilities_avg[i] 
                            for i in range(len(Epses))]
            self.error_utilities = [self.lower_error, self.upper_error]
            plt.errorbar(Epses, self.utilities_avg,
                         yerr=self.error_utilities,
                         fmt=self.ls,
                         capsize=eb_capsize,
                         label = self.label)
    
    def plotTradeoff_Eps_Vr_Dts(self):
        if self.en:
            self.dts_avg = [st.mean(self.dts[i]) 
                                    for i in range(len(Epses))]
            self.upper_error = [max(self.dts[i])-self.dts_avg[i] 
                            for i in range(len(Epses))]
            self.lower_error = [-min(self.dts[i])+self.dts_avg[i] 
                            for i in range(len(Epses))]
            self.error_dts = [self.lower_error, self.upper_error]
            plt.errorbar(Epses, self.dts_avg,
                         yerr=self.error_dts,
                         fmt=self.ls,
                         capsize=eb_capsize,
                         label = self.label)
    
    def plotTradeoff_Eps_Vr_Steps(self):
        if self.en:
            self.steps_avg = [st.mean(self.steps[i]) 
                                    for i in range(len(Epses))]
            self.upper_error = [max(self.steps[i])-self.steps_avg[i] 
                            for i in range(len(Epses))]
            self.lower_error = [-min(self.steps[i])+self.steps_avg[i] 
                            for i in range(len(Epses))]
            self.error_steps = [self.lower_error, self.upper_error]
            plt.errorbar(Epses, self.steps_avg,
                         yerr=self.error_steps,
                         fmt=self.ls,
                         capsize=eb_capsize,
                         label = self.label)
    
#    def plotTradeoff_Eps_Vr_Evs(self):
#        if self.en:
#            self.evs_avg = [st.mean(self.evs[i]) 
#                                    for i in range(len(Epses))]
#            self.upper_error = [max(self.evs[i])-self.evs_avg[i] 
#                            for i in range(len(Epses))]
#            self.lower_error = [-min(self.evs[i])+self.evs_avg[i] 
#                            for i in range(len(Epses))]
#            self.error_evs = [self.lower_error, self.upper_error]
#            plt.errorbar(Epses, self.evs_avg,
#                         yerr=self.error_evs,
#                         fmt=self.ls,
#                         capsize=eb_capsize,
#                         label = self.label)
    
    def plotTradeoff_Eps_Vr_Evs(self): 
        if self.en:
            self.evs_avg = [st.mean(self.evs[i]) 
                                    for i in range(len(Epses))]
            self.upper_error = [max(self.evs[i])-self.evs_avg[i] 
                            for i in range(len(Epses))]
            self.lower_error = [-min(self.evs[i])+self.evs_avg[i] 
                            for i in range(len(Epses))]
            self.error_evs = [self.lower_error, self.upper_error]
            self.evs_avg_scale = [evs/10**evs_scale for evs in self.evs_avg]
            self.error_evs_scale = [[err/10**evs_scale for err in self.lower_error], 
                                    [err/10**evs_scale for err in self.upper_error]]
            plt.errorbar(Epses, self.evs_avg_scale,
                         yerr=self.error_evs_scale,
                         fmt=self.ls,
                         capsize=eb_capsize,
                         label = self.label)
    
    # If necessary, the memory of sum and temp can be released.

# =============================================================================
# 
# =============================================================================
# line styles for plot figures
class Linestyle:
    
    ls_GA = 'rx-'
    
    ls_AuctionSequential = 'b-o'
    ls_AuctionParallel = 'c--x'
    ls_AuctionGPrim = 'r-.+'
    
    ls_SGA = 'm-o'
    ls_LSGA = 'm--o'
    
    ls_CBBA = 'r-*'

    ls_TGTA = 'g-x'
    ls_LTGTA = 'g--x'

    ls_DTTA = 'g.-'
    ls_LDTTA = 'g--x'
    ls_TBTA = 'g-.+'
    
    ls_T3A = 'b.-'
    ls_LT3A = 'b--x'
    ls_TTBTA = 'b+-.'
      
    ls_DSTA = 'b.-'    
    ls_LSTA = 'b.--'
    
    ls_STTA = 'c.-'
    ls_LSTTA = 'c--x'
    ls_STBTA = 'c+-.'
    
    # error bar styles
    eb_STTA = 'c.-' 
    eb_LSTTA = 'b.--'
    eb_STBTA = 'g.-.'
    
#    ls_DSTA = 'c-d'    
#    ls_LSTA = 'c--d'
#    
#    ls_STTA = 'g-h'
#    ls_LSTTA = 'g--h'
#    ls_STBTA = 'g-.h'
    

# algorithms' labels
class Label:
    
    lb_GA = 'GA'
    
    lb_AuctionSequential = 'Auction-Sequential'
    lb_AuctionParallel = 'Auction-Parallel'
    lb_AuctionGPrim = 'Auction-G-Prim'
    
    lb_SGA = 'SGA'
    lb_LSGA = 'LSGA'
    
    lb_CBBA = 'CBBA'
    
    lb_TGTA = 'TGTA'
    lb_LTGTA = 'LTGTA'

    lb_DTTA = 'DTTA'
    lb_LDTTA = 'LDTTA'
    lb_TBTA = 'TBTA'
    
    lb_T3A = 'T3A'
    lb_LT3A = 'LT3A'
    lb_TTBTA = 'TTBTA'
      
    lb_DSTA = 'DSTA'    
    lb_LSTA = 'LSTA'
    
    lb_STTA = 'STTA'
    lb_LSTTA = 'LSTTA'
    lb_STBTA = 'STBTA'
    
    
        
# =============================================================================
# Initialise algorithms.
# =============================================================================

GA = Algorithm(Linestyle.ls_GA, Label.lb_GA)

AuctionSequential = Algorithm(Linestyle.ls_AuctionSequential, Label.lb_AuctionSequential)
AuctionParallel = Algorithm(Linestyle.ls_AuctionParallel, Label.lb_AuctionParallel)
AuctionGPrim = Algorithm(Linestyle.ls_AuctionGPrim, Label.lb_AuctionGPrim)

SGA = Algorithm(Linestyle.ls_SGA, Label.lb_SGA)
LSGA = Algorithm(Linestyle.ls_LSGA, Label.lb_LSGA)

CBBA = Algorithm(Linestyle.ls_CBBA, Label.lb_CBBA)

TGTA = Algorithm(Linestyle.ls_TGTA, Label.lb_TGTA)
LTGTA = Algorithm(Linestyle.ls_LTGTA, Label.lb_LTGTA)

DTTA = Algorithm(Linestyle.ls_DTTA, Label.lb_DTTA)
LDTTA = Algorithm(Linestyle.ls_LDTTA, Label.lb_LDTTA)
TBTA = Algorithm(Linestyle.ls_TBTA, Label.lb_TBTA)

T3A = Algorithm(Linestyle.ls_T3A, Label.lb_T3A)
LT3A = Algorithm(Linestyle.ls_LT3A, Label.lb_LT3A)
TTBTA = Algorithm(Linestyle.ls_TTBTA, Label.lb_TTBTA)


DSTA = Algorithm(Linestyle.ls_DSTA, Label.lb_DSTA)
LSTA = Algorithm(Linestyle.ls_LSTA, Label.lb_LSTA)

STTA = Algorithm(Linestyle.ls_STTA, Label.lb_STTA, Linestyle.eb_STTA)
LSTTA = Algorithm(Linestyle.ls_LSTTA, Label.lb_LSTTA, Linestyle.eb_LSTTA)
STBTA = Algorithm(Linestyle.ls_STBTA, Label.lb_STBTA, Linestyle.eb_STBTA)



# =============================================================================
# 
# =============================================================================
print("----- init.py is loaded -----")

