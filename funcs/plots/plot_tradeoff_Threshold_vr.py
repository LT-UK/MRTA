#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 16 12:08:04 2020

Plot tradeoffs related to the threshold.

@author: Teng Li
tengli@cranfield.ac.uk; nliteng@foxmail.com
Cranfield University, UK
All Rights Reserved
"""

import matplotlib.pyplot as plt
import os # find path
import init
import numpy as np

# =============================================================================
# 
# =============================================================================

# monotone
output_path_mono = os.path.abspath("output") + "/tradeoff/epsilon_vr/monotone/fig/"
# non-monotone
output_path_non = os.path.abspath("output") + "/tradeoff/epsilon_vr/non_monotone/fig/"

def plotTradeoff_Eps_Vr():
    # Draw figures
    print("\n========= Plot Tradeoff of Threshold =========")
    
# total values             
    plt.figure("tradeoff utility eps vr")
    init.STTA.plotTradeoff_Eps_Vr_Utilities()
    init.LSTTA.plotTradeoff_Eps_Vr_Utilities()
    init.STBTA.plotTradeoff_Eps_Vr_Utilities()
    # Plot figure
    plt.xlabel("$\epsilon$")
    plt.ylabel("Function Utility")
    plt.legend(loc = 'lower right')
    plt.ylim(bottom=0,top=30)
    if init.monotonicity:
        plt.savefig(output_path_mono+'tradeoff_Eps_mono_utility.pdf')
    else:
        plt.savefig(output_path_non+'tradeoff_Eps_non_utility.pdf')
    

# number of consensus steps       
    plt.figure("tradeoff steps eps vr")
    init.STTA.plotTradeoff_Eps_Vr_Steps()
    init.LSTTA.plotTradeoff_Eps_Vr_Steps()
    init.STBTA.plotTradeoff_Eps_Vr_Steps()
    # Plot figure
    plt.xlabel("$\epsilon$")
    plt.ylabel("Consensus Steps")
    plt.legend(loc = 'best')
    plt.ylim(bottom=0,top=45)
#    plt.ylim(bottom=0, top=init.Nt*1.1)
    if init.monotonicity:
        plt.savefig(output_path_mono+'tradeoff_Eps_mono_steps.pdf')
    else:
        plt.savefig(output_path_non+'tradeoff_Eps_non_steps.pdf')
    
    
#number of function evaluations     
    plt.figure("tradeoff evs eps vr")
    init.STTA.plotTradeoff_Eps_Vr_Evs()
    init.LSTTA.plotTradeoff_Eps_Vr_Evs()
    init.STBTA.plotTradeoff_Eps_Vr_Evs()
    # Plot figure    
    plt.xlabel("$\epsilon$")
    plt.ylabel("Function Evaluations (x$10^"+str(init.evs_scale)+"$)")
    plt.legend(loc = 'best')
    plt.ylim(bottom=0,top=3)
    if init.monotonicity:
        plt.savefig(output_path_mono+'tradeoff_Eps_mono_evs.pdf')
    else:
        plt.savefig(output_path_non+'tradeoff_Eps_non_evs.pdf')
    

def set_box_color(bp, color):
    plt.setp(bp['boxes'], color=color)
    plt.setp(bp['whiskers'], color=color)
    plt.setp(bp['caps'], color=color)
    plt.setp(bp['medians'], color=color)
    
    
def plotTradeoff_Eps_Box_STTA():
    
    # plot utilities
    plt.figure("boxplot_tradeoff_Eps_utilities_STTA")
    # ticks =  np.arange(0.04, 0.241, 0.04)
    ticks = [round(Eps,2) for Eps in init.Epses]
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
    plt.xlabel("$\epsilon$")
    plt.ylabel("Function Utility")
    plt.legend(loc = 'lower right')
    plt.ylim(ymin=0,ymax=40)
    if init.monotonicity:
        plt.savefig(output_path_mono+'STTA_tradeoff_Eps_mono_utility_box.pdf')
    else:
        plt.savefig(output_path_non+'STTA_tradeoff_Eps_non_utility_box.pdf')
        
# plot steps
    plt.figure("boxplot_tradeoff_Eps_steps_STTA")
    # ticks =  np.arange(0.04, 0.241, 0.04)
    ticks = [round(Eps,2) for Eps in init.Epses]
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
    plt.xlabel("$\epsilon$")
    plt.ylabel("Consensus Steps")
    plt.legend(loc = 'center right')
    plt.ylim(ymin=0,ymax=50)
    if init.monotonicity:
        plt.savefig(output_path_mono+'STTA_tradeoff_Eps_mono_steps_box.pdf')
    else:
        plt.savefig(output_path_non+'STTA_tradeoff_Eps_non_steps_box.pdf')


# plot dts
    plt.figure("boxplot_tradeoff_Eps_dts")
    # ticks =  np.arange(0.04, 0.241, 0.04)
    ticks = [round(Eps,2) for Eps in init.Epses]
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
    plt.xlabel("$\epsilon$")
    plt.ylabel("Running Time (sec)")
    plt.legend(loc = 'center right')
    plt.ylim(ymin=0)
    if init.monotonicity:
        plt.savefig(output_path_mono+'STTA_tradeoff_Eps_mono_dts_box.pdf')
    else:
        plt.savefig(output_path_non+'STTA_tradeoff_Eps_non_dts_box.pdf')        
    
    # plot evs
    plt.figure("STTA_boxplot_tradeoff_Eps_evs")
    # ticks =  np.arange(0.04, 0.241, 0.04)
    ticks = [round(Eps,2) for Eps in init.Epses]
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
    plt.xlabel("$\epsilon$")
    plt.ylabel("Function Evaluations (x$10^"+str(init.evs_scale)+"$)")
    plt.legend(loc = 'center right')
    plt.ylim(ymin=0,ymax=18)
    if init.monotonicity:
        plt.savefig(output_path_mono+'STTA_tradeoff_Eps_mono_evs_box.pdf')
    else:
        plt.savefig(output_path_non+'STTA_tradeoff_Eps_non_evs_box.pdf')  



# =============================================================================
# 
# =============================================================================
print("----- plot_tradeoff_Threshold.py is loaded -----")
