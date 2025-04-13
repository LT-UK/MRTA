#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 16 23:59:58 2020

@author: Teng Li
tengli@cranfield.ac.uk; nliteng@foxmail.com
Cranfield University, UK
All Rights Reserved
"""

import matplotlib.pyplot as plt
import os # find path
import init
import numpy as np
#import statistics as st

# =============================================================================
# 
# =============================================================================
# monotone
output_path_mono = os.path.abspath("output") \
                    + "/tradeoff/probability/monotone/fig/"
# non-monotone
output_path_non = os.path.abspath("output") \
                    + "/tradeoff/probability/non_monotone/fig/"

# %%
def plotTradeoff_Pr():
    # Draw figures
    print("\n======= Plot Tradeoff of Pr =======")
    
# total values             
    plt.figure("utility tradeoff")
    init.DSTA.plotTradeoff_Pr_Utilities()
    init.LSTA.plotTradeoff_Pr_Utilities()
    init.STTA.plotTradeoff_Pr_Utilities()
    init.LSTTA.plotTradeoff_Pr_Utilities()
    init.STBTA.plotTradeoff_Pr_Utilities()
    # test
    plt.plot(init.Prs, init.CBBA.utilities, init.Linestyle.ls_CBBA, label = init.Label.lb_CBBA)
    # Plot figure
    plt.xlabel("$p$")
    plt.ylabel("Function Utility")
    plt.legend(loc = 'lower right')
    plt.ylim(bottom=0)
    if init.monotonicity:
        plt.savefig(output_path_mono+'tradeoff_Pr_mono_utility.pdf')
    else:
        plt.savefig(output_path_non+'tradeoff_Pr_non_utility.pdf')

## consensus steps             
    plt.figure("steps tradeoff")
    init.DSTA.plotTradeoff_Pr_Steps()
    init.LSTA.plotTradeoff_Pr_Steps()
    init.STTA.plotTradeoff_Pr_Steps()
    init.LSTTA.plotTradeoff_Pr_Steps()
    init.STBTA.plotTradeoff_Pr_Steps()
    # test
    plt.plot(init.Prs, init.CBBA.steps, init.Linestyle.ls_CBBA, label = init.Label.lb_CBBA)
    # Plot figure
    plt.xlabel("$p$")
    plt.ylabel("Consensus Steps")
    plt.legend(loc = 'best')
    plt.ylim(bottom=0)
    if init.monotonicity:
        plt.savefig(output_path_mono+'tradeoff_Pr_mono_steps.pdf')
    else:
        plt.savefig(output_path_non+'tradeoff_Pr_non_steps.pdf')
#
#
## number of function evaluations            
    plt.figure("evs tradeoff")
    init.DSTA.plotTradeoff_Pr_Evs()
    init.LSTA.plotTradeoff_Pr_Evs()
    init.STTA.plotTradeoff_Pr_Evs()
    init.LSTTA.plotTradeoff_Pr_Evs()
    init.STBTA.plotTradeoff_Pr_Evs()
    # test
    plt.plot(init.Prs, init.CBBA.evs, init.Linestyle.ls_CBBA, label = init.Label.lb_CBBA)
    # Plot figure
    plt.xlabel("$p$")
    plt.ylabel("Function Evaluations (x$10^"+str(init.evs_scale)+"$)")
    plt.legend(loc = 'best')
    plt.ylim(bottom=0)
    if init.monotonicity:
        plt.savefig(output_path_mono+'tradeoff_Pr_mono_evs.pdf')
    else:
        plt.savefig(output_path_non+'tradeoff_Pr_non_evs.pdf')

#%%
def set_box_color(bp, color):
    plt.setp(bp['boxes'], color=color)
    plt.setp(bp['whiskers'], color=color)
    plt.setp(bp['caps'], color=color)
    plt.setp(bp['medians'], color=color)

#%%
def plotTradeoff_Pr_Box(): # DSTA
    
    # plot utilities
    plt.figure("boxplot_tradeoff_utilities")
    # ticks =  np.arange(0.04, 0.241, 0.04)
    ticks = [round(Pr,2) for Pr in init.Prs]
    bp = plt.boxplot(init.DSTA.utilities, 
                     positions=np.array(range(len(init.DSTA.utilities)))*2.0, 
                     sym='', widths=0.4)
    
    set_box_color(bp, 'c')
    plt.plot([], c='c', label='DSTA')
    plt.legend()
    plt.xticks(range(0, len(ticks) * 2, 2), ticks)
    plt.xlim(-1, len(ticks)*2-1)
    plt.ylim(bottom=0)
#        plt.tight_layout()
    # test
#    plt.plot(init.Prs, init.CBBA.utilities, init.Linestyle.ls_CBBA, label = init.Label.lb_CBBA)
    plt.axhline(init.CBBA.utilities[0], c = "Red", label = init.Label.lb_CBBA)
#    plt.abline(h=init.CBBA.utilities[0], col = "Red")
    plt.xlabel("Sampling Probability $p$")
    plt.ylabel("Function Utility")
    plt.legend(loc = 'lower right')
    plt.ylim(ymin=0)
    if init.monotonicity:
        plt.savefig(output_path_mono+'tradeoff_Pr_mono_utility_box.pdf')
    else:
        plt.savefig(output_path_non+'tradeoff_Pr_non_utility_box.pdf')
    
    
# plot steps
    plt.figure("boxplot_tradeoff_steps")
    # ticks =  np.arange(0.04, 0.241, 0.04)
    ticks = [round(Pr,2) for Pr in init.Prs]
    bp = plt.boxplot(init.DSTA.steps, 
                     positions=np.array(range(len(init.DSTA.steps)))*2.0, 
                     sym='', widths=0.4)
    
    set_box_color(bp, 'c')
    plt.plot([], c='c', label='DSTA')
    plt.legend()
    plt.xticks(range(0, len(ticks) * 2, 2), ticks)
    plt.xlim(-1, len(ticks)*2-1)
    plt.ylim(bottom=0)
#        plt.tight_layout()
    # test
    plt.axhline(init.CBBA.steps[0], c = "Red", label = init.Label.lb_CBBA)
    plt.xlabel("Sampling Probability $p$")
    plt.ylabel("Consensus Steps")
    plt.legend(loc = 'lower right')
    plt.ylim(ymin=0)
    if init.monotonicity:
        plt.savefig(output_path_mono+'tradeoff_Pr_mono_steps_box.pdf')
    else:
        plt.savefig(output_path_non+'tradeoff_Pr_non_steps_box.pdf')


# plot dts
    plt.figure("boxplot_tradeoff_dts")
    # ticks =  np.arange(0.04, 0.241, 0.04)
    ticks = [round(Pr,2) for Pr in init.Prs]
    bp = plt.boxplot(init.DSTA.dts, 
                     positions=np.array(range(len(init.DSTA.dts)))*2.0, 
                     sym='', widths=0.4)
    
    set_box_color(bp, 'c')
    plt.plot([], c='c', label='DSTA')
    plt.legend()
    plt.xticks(range(0, len(ticks) * 2, 2), ticks)
    plt.xlim(-1, len(ticks)*2-1)
    plt.ylim(bottom=0)
#        plt.tight_layout()
    # test
    plt.axhline(init.CBBA.dts[0], c = "Red", label = init.Label.lb_CBBA)
    plt.xlabel("Sampling Probability $p$")
    plt.ylabel("Running Time (sec)")
    plt.legend(loc = 'center right')
    plt.ylim(ymin=0,ymax=1)
    if init.monotonicity:
        plt.savefig(output_path_mono+'tradeoff_Pr_mono_dts_box.pdf')
    else:
        plt.savefig(output_path_non+'tradeoff_Pr_non_dts_box.pdf')        
    
    # plot evs
    plt.figure("boxplot_tradeoff_evs")
    # ticks =  np.arange(0.04, 0.241, 0.04)
    ticks = [round(Pr,2) for Pr in init.Prs]
#    evs_DSTA_scaled = [evs/10**init.evs_scale for evs in init.DSTA.evs]
    evs_DSTA_scaled = [[evs/10**init.evs_scale for evs in init.DSTA.evs[i]] for i in range(len(init.DSTA.evs))]
    bp = plt.boxplot(evs_DSTA_scaled, 
                     positions=np.array(range(len(init.DSTA.evs)))*2.0, 
                     sym='', widths=0.4)
    
    set_box_color(bp, 'c')
    plt.plot([], c='c', label='DSTA')
    plt.legend()
    plt.xticks(range(0, len(ticks) * 2, 2), ticks)
    plt.xlim(-1, len(ticks)*2-1)
    plt.ylim(bottom=0)
#        plt.tight_layout()
    # test
    plt.axhline(init.CBBA.evs[0], c = "Red", label = init.Label.lb_CBBA)
    plt.xlabel("Sampling Probability $p$")
    plt.ylabel("Function Evaluations (x$10^"+str(init.evs_scale)+"$)")
    plt.legend(loc = 'center right')
    plt.ylim(ymin=0,ymax=200)
    if init.monotonicity:
        plt.savefig(output_path_mono+'tradeoff_Pr_mono_evs_box.pdf')
    else:
        plt.savefig(output_path_non+'tradeoff_Pr_non_evs_box.pdf')
    
    
#    # plot utilities
#    plt.figure("boxplot_tradeoff_utilities")
#    # ticks = np.arange(0.04, 0.241, 0.04)
#    ticks = [0.04, 0.08, 0.12, 0.16, 0.20, 0.24]
#    bpl = plt.boxplot(values_threshold_1, positions=np.array(range(len(values_threshold_1)))*2.0-0.3, sym='', widths=0.4)
#    bpr = plt.boxplot(values_threshold_3, positions=np.array(range(len(values_threshold_3)))*2.0+0.3, sym='', widths=0.4)
#    set_box_color(bpl, 'b')
#    set_box_color(bpr, 'r')    
#    plt.plot([], c='b', label='$m_g=2$')
#    plt.plot([], c='r', label='$m_g=5$')
#    plt.legend()
#    plt.xticks(range(0, len(ticks) * 2, 2), ticks)
#    plt.xlim(-1, len(ticks)*2-1)
#    plt.ylim(0, 6800)
##        plt.tight_layout()
#    
#    plt.xlabel("Threshold Decreasing Parameter $\epsilon$")
#    plt.ylabel("Function Value $f(S)$")
#    plt.legend(loc = 'best')
#    plt.ylim(ymin=0)
#    plt.savefig('boxplot_tradeoff_utility.eps')
        
        
        
#%%        
def plotTradeoff_Pr_Box_STTA():
    
    # plot utilities
    plt.figure("boxplot_tradeoff_utilities_STTA")
    # ticks =  np.arange(0.04, 0.241, 0.04)
    ticks = [round(Pr,2) for Pr in init.Prs]
    bp = plt.boxplot(init.STTA.utilities, 
                     positions=np.array(range(len(init.STTA.utilities)))*2.0, 
                     sym='', widths=0.4)
    
    set_box_color(bp, 'c')
    plt.plot([], c='c', label='STTA')
    plt.legend()
    plt.xticks(range(0, len(ticks) * 2, 2), ticks)
    plt.xlim(-1, len(ticks)*2-1)
    plt.ylim(bottom=0)
#        plt.tight_layout()
    # test
#    plt.plot(init.Prs, init.CBBA.utilities, init.Linestyle.ls_CBBA, label = init.Label.lb_CBBA)
    plt.axhline(init.SGA.utilities[0], c = "m", label = init.Label.lb_SGA)
#    plt.axhline(init.CBBA.utilities[0], c = "Red", label = init.Label.lb_CBBA)
#    plt.abline(h=init.CBBA.utilities[0], col = "Red")
    plt.xlabel("Sampling Probability $p$")
    plt.ylabel("Function Utility")
    plt.legend(loc = 'lower right')
    plt.ylim(ymin=0,ymax=40)
    if init.monotonicity:
        plt.savefig(output_path_mono+'STTA_tradeoff_Pr_mono_utility_box.pdf')
    else:
        plt.savefig(output_path_non+'STTA_tradeoff_Pr_non_utility_box.pdf')
        
        
# plot steps
    plt.figure("boxplot_tradeoff_steps_STTA")
    # ticks =  np.arange(0.04, 0.241, 0.04)
    ticks = [round(Pr,2) for Pr in init.Prs]
    bp = plt.boxplot(init.STTA.steps, 
                     positions=np.array(range(len(init.STTA.steps)))*2.0, 
                     sym='', widths=0.4)
    
    set_box_color(bp, 'c')
    plt.plot([], c='c', label='STTA')
    plt.legend()
    plt.xticks(range(0, len(ticks) * 2, 2), ticks)
    plt.xlim(-1, len(ticks)*2-1)
    plt.ylim(bottom=0)
#        plt.tight_layout()
    # test
    plt.axhline(init.SGA.steps[0], c = "m", label = init.Label.lb_SGA)
#    plt.axhline(init.CBBA.steps[0], c = "Red", label = init.Label.lb_CBBA)
    plt.xlabel("Sampling Probability $p$")
    plt.ylabel("Consensus Steps")
    plt.legend(loc = 'lower right')
    plt.ylim(ymin=0,ymax=50)
    if init.monotonicity:
        plt.savefig(output_path_mono+'STTA_tradeoff_Pr_mono_steps_box.pdf')
    else:
        plt.savefig(output_path_non+'STTA_tradeoff_Pr_non_steps_box.pdf')


# plot dts
    plt.figure("boxplot_tradeoff_dts")
    # ticks =  np.arange(0.04, 0.241, 0.04)
    ticks = [round(Pr,2) for Pr in init.Prs]
    bp = plt.boxplot(init.STTA.dts, 
                     positions=np.array(range(len(init.STTA.dts)))*2.0, 
                     sym='', widths=0.4)
    
    set_box_color(bp, 'c')
    plt.plot([], c='c', label='STTA')
    plt.legend()
    plt.xticks(range(0, len(ticks) * 2, 2), ticks)
    plt.xlim(-1, len(ticks)*2-1)
    plt.ylim(bottom=0)
#        plt.tight_layout()
    # test
    plt.axhline(init.SGA.dts[0], c = "m", label = init.Label.lb_SGA)
#    plt.axhline(init.CBBA.dts[0], c = "Red", label = init.Label.lb_CBBA)
    plt.xlabel("Sampling Probability $p$")
    plt.ylabel("Running Time (sec)")
    plt.legend(loc = 'center right')
    plt.ylim(ymin=0)
    if init.monotonicity:
        plt.savefig(output_path_mono+'STTA_tradeoff_Pr_mono_dts_box.pdf')
    else:
        plt.savefig(output_path_non+'STTA_tradeoff_Pr_non_dts_box.pdf')        
    
    # plot evs
    plt.figure("STTA_boxplot_tradeoff_evs")
    # ticks =  np.arange(0.04, 0.241, 0.04)
    ticks = [round(Pr,2) for Pr in init.Prs]
#    evs_DSTA_scaled = [evs/10**init.evs_scale for evs in init.DSTA.evs]
    evs_STTA_scaled = [[evs/10**init.evs_scale for evs in init.STTA.evs[i]] for i in range(len(init.STTA.evs))]
    bp = plt.boxplot(evs_STTA_scaled, 
                     positions=np.array(range(len(init.STTA.evs)))*2.0, 
                     sym='', widths=0.4)
    
    set_box_color(bp, 'c')
    plt.plot([], c='c', label='STTA')
    plt.legend()
    plt.xticks(range(0, len(ticks) * 2, 2), ticks)
    plt.xlim(-1, len(ticks)*2-1)
    plt.ylim(bottom=0)
#        plt.tight_layout()
    # test
    plt.axhline(init.SGA.evs[0], c = "m", label = init.Label.lb_SGA)
#    plt.axhline(init.CBBA.evs[0], c = "Red", label = init.Label.lb_CBBA)
    plt.xlabel("Sampling Probability $p$")
    plt.ylabel("Function Evaluations (x$10^"+str(init.evs_scale)+"$)")
    plt.legend(loc = 'center left')
    plt.ylim(ymin=0,ymax=18)
    if init.monotonicity:
        plt.savefig(output_path_mono+'STTA_tradeoff_Pr_mono_evs_box.pdf')
    else:
        plt.savefig(output_path_non+'STTA_tradeoff_Pr_non_evs_box.pdf')
        
#%%
# =============================================================================
# 
# =============================================================================
print("----- plot_tradeoff_Pr.py is loaded -----")
