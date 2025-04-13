#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 16 23:59:58 2020

@author: Teng Li
lt.uk@outlook.com
United Kingdom
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
output_path_mono = os.path.abspath("output") + "/variance/monotone/plot_stta_vs_dsta/"
# non-monotone
output_path_non = os.path.abspath("output") + "/variance/non_monotone/plot_stta_vs_dsta/"



def plotSTTA_DSTA_Vr():
    # Draw figures
    print("\n======= Plot DSTA vs STTA Variance =======")
        
# total values             
    plt.figure("stta utility variance")
    init.SGA.plotUtilities()
    init.CBBA.plotUtilities()  
    init.STTA.plotVariance_Utilities()
    init.DSTA.plotVariance_Utilities()
    # Plot figure
    plt.xlabel("$N_a$")
    plt.ylabel("Function Utility")
    plt.legend(loc = 'lower right')
    plt.ylim(bottom=0)
    if init.monotonicity:
        plt.savefig(output_path_mono+'STTA_DSTA_utility_mono_vr.pdf')
    else:
        plt.savefig(output_path_non+'STTA_DSTA_utility_non_vr.pdf')

# consensus steps             
    plt.figure("stta steps variance")
    init.SGA.plotSteps()
    init.CBBA.plotSteps()  
    init.STTA.plotVariance_Steps()
    init.DSTA.plotVariance_Steps()
    # Plot figure
    plt.xlabel("$N_a$")
    plt.ylabel("Consensus Steps")
#    plt.legend(loc='center left',bbox_to_anchor=(0.76, 0.6))
    plt.ylim(bottom=0, top=init.Nt*1.1)
    if init.monotonicity:
        plt.legend(loc='center left',bbox_to_anchor=(0.76, 0.6))
        plt.savefig(output_path_mono+'STTA_DSTA_steps_mono_vr.pdf')
    else:
        plt.legend(loc='best')
        plt.savefig(output_path_non+'STTA_DSTA_steps_non_vr.pdf')


# number of function evaluations no CBBA           
    plt.figure("stta evs variance no CBBA")
#    init.SGA.plotEvs()
#    init.CBBA.plotEvs()  
    init.STTA.plotVariance_Evs()
    init.DSTA.plotVariance_Evs()

    # Plot figure
    plt.xlabel("$N_a$")
    plt.ylabel("Function Evaluations (x$10^"+str(init.evs_scale)+"$)")
    plt.legend(loc = 'best')
    plt.ylim(bottom=0)
    if init.monotonicity:
        plt.savefig(output_path_mono+'STTA_DSTA_evs_mono_vr_noCBBA.pdf')
    else:
        plt.savefig(output_path_non+'STTA_DSTA_evs_non_vr_noCBBA.pdf')
        
    # number of function evaluations            
    plt.figure("stta evs variance")
    init.SGA.plotEvs()
    init.CBBA.plotEvs()  
    init.STTA.plotVariance_Evs()
    init.DSTA.plotVariance_Evs()
    # Plot figure
    plt.xlabel("$N_a$")
    plt.ylabel("Function Evaluations (x$10^"+str(init.evs_scale)+"$)")
    plt.legend(loc = 'best')
    plt.ylim(bottom=0)
    if init.monotonicity:
        plt.savefig(output_path_mono+'STTA_DSTA_evs_mono_vr.pdf')
    else:
        plt.savefig(output_path_non+'STTA_DSTA_evs_non_vr.pdf')
        
    
def set_box_color(bp, color):
    plt.setp(bp['boxes'], color=color)
    plt.setp(bp['whiskers'], color=color)
    plt.setp(bp['caps'], color=color)
    plt.setp(bp['medians'], color=color)
    
    
def plotSTTA_DSTA_Vr_Box():
    
    # plot utilities
    plt.figure("boxplot_comparison_utilities")
    # ticks =  np.arange(0.04, 0.241, 0.04)
    ticks = [Na for Na in init.Nas]
    
    SGA_utilities = [[init.SGA.utilities[i]] for i in range(len(init.SGA.utilities))]
    bp1 = plt.boxplot(SGA_utilities, 
                     positions=np.array(range(len(init.SGA.utilities)))*2.0, 
                     sym='', widths=0.4)
    set_box_color(bp1, 'm')
    plt.plot([], c='m', label='SGA')
    
    CBBA_utilities = [[init.CBBA.utilities[i]] for i in range(len(init.CBBA.utilities))]
    bp2 = plt.boxplot(CBBA_utilities, 
                     positions=np.array(range(len(init.CBBA.utilities)))*2.0, 
                     sym='', widths=0.4)
    set_box_color(bp2, 'r')
    plt.plot([], c='r', label='CBBA')
    
    bp3 = plt.boxplot(init.DSTA.utilities, 
                     positions=np.array(range(len(init.DSTA.utilities)))*2.0, 
                     sym='', widths=0.4)
    set_box_color(bp3, 'b')
    plt.plot([], c='b', label='DSTA')
    
    bp4 = plt.boxplot(init.STTA.utilities, 
                     positions=np.array(range(len(init.STTA.utilities)))*2.0, 
                     sym='', widths=0.4)
    set_box_color(bp4, 'c')
    plt.plot([], c='c', label='STTA')
    
#    
#    bp5 = plt.boxplot(init.STBTA.utilities, 
#                     positions=np.array(range(len(init.STBTA.utilities)))*2.0, 
#                     sym='', widths=0.4)
#    set_box_color(bp5, 'g')
#    plt.plot([], c='g', label='STBTA')
    
    plt.legend()
    plt.xticks(range(0, len(ticks) * 2, 2), ticks)
    plt.xlim(-1, len(ticks)*2-1)
    plt.ylim(bottom=0)
#        plt.tight_layout()
    # test
#    plt.plot(init.Prs, init.CBBA.utilities, init.Linestyle.ls_CBBA, label = init.Label.lb_CBBA)
#    plt.axhline(init.CBBA.utilities[0], c = "Red", label = init.Label.lb_CBBA)
#    plt.abline(h=init.CBBA.utilities[0], col = "Red")
    plt.xlabel("$N_a$")
    plt.ylabel("Function Utility")
    plt.legend(loc = 'lower right')
    plt.ylim(ymin=0)
    if init.monotonicity:
        plt.savefig(output_path_mono+'STTA_DSTA_utility_mono_vr_box.pdf')
    else:
        plt.savefig(output_path_non+'STTA_DSTA_utility_non_vr_box.pdf')
        
        
    # plot evs
    plt.figure("boxplot_comparison_evs")
    # ticks =  np.arange(0.04, 0.241, 0.04)
    ticks = [Na for Na in init.Nas]
    
    SGA_evs = [[init.SGA.evs[i]/10**init.evs_scale] for i in range(len(init.SGA.evs))]
    bp1 = plt.boxplot(SGA_evs, 
                     positions=np.array(range(len(init.SGA.evs)))*2.0, 
                     sym='', widths=0.4)
    set_box_color(bp1, 'm')
    plt.plot([], c='m', label='SGA')
    
    CBBA_evs = [[init.CBBA.evs[i]/10**init.evs_scale] for i in range(len(init.CBBA.evs))]
    bp2 = plt.boxplot(CBBA_evs, 
                     positions=np.array(range(len(init.CBBA.evs)))*2.0, 
                     sym='', widths=0.4)
    set_box_color(bp2, 'r')
    plt.plot([], c='r', label='CBBA')
    
    DSTA_evs = [[init.DSTA.evs[i][j]/10**init.evs_scale for j in range(len(init.DSTA.evs[i]))] for i in range(len(init.DSTA.evs))]
    bp3 = plt.boxplot(DSTA_evs, 
                     positions=np.array(range(len(init.DSTA.evs)))*2.0, 
                     sym='', widths=0.4)
    set_box_color(bp3, 'b')
    plt.plot([], c='b', label='DSTA')
    
    STTA_evs = [[init.STTA.evs[i][j]/10**init.evs_scale for j in range(len(init.STTA.evs[i]))] for i in range(len(init.STTA.evs))]
    bp4 = plt.boxplot(STTA_evs, 
                     positions=np.array(range(len(init.STTA.evs)))*2.0, 
                     sym='', widths=0.4)
    set_box_color(bp4, 'c')
    plt.plot([], c='c', label='STTA')
    
    
    
    plt.legend()
    plt.xticks(range(0, len(ticks) * 2, 2), ticks)
    plt.xlim(-1, len(ticks)*2-1)
    plt.ylim(bottom=0)
    plt.xlabel("$N_a$")
    plt.ylabel("Function Evaluations (x$10^"+str(init.evs_scale)+"$)")
    plt.legend(loc = 'upper left')
    plt.ylim(ymin=0)
    if init.monotonicity:
        plt.savefig(output_path_mono+'STTA_DSTA_evs_mono_vr_box.pdf')
    else:
        plt.savefig(output_path_non+'STTA_DSTA_evs_non_vr_box.pdf')
        
        
    # plot steps
    plt.figure("boxplot_comparison_steps")
    # ticks =  np.arange(0.04, 0.241, 0.04)
    ticks = [Na for Na in init.Nas]
    
    SGA_steps = [[init.SGA.steps[i]] for i in range(len(init.SGA.steps))]
    bp1 = plt.boxplot(SGA_steps, 
                     positions=np.array(range(len(init.SGA.steps)))*2.0, 
                     sym='', widths=0.4)
    set_box_color(bp1, 'm')
    plt.plot([], c='m', label='SGA')
    
    CBBA_steps = [[init.CBBA.steps[i]] for i in range(len(init.CBBA.steps))]
    bp2 = plt.boxplot(CBBA_steps, 
                     positions=np.array(range(len(init.CBBA.steps)))*2.0, 
                     sym='', widths=0.4)
    set_box_color(bp2, 'r')
    plt.plot([], c='r', label='CBBA')
    
    bp3 = plt.boxplot(init.DSTA.steps, 
                     positions=np.array(range(len(init.DSTA.steps)))*2.0, 
                     sym='', widths=0.4)
    set_box_color(bp3, 'b')
    plt.plot([], c='b', label='DSTA')
    
    bp4 = plt.boxplot(init.STTA.steps, 
                     positions=np.array(range(len(init.STTA.steps)))*2.0, 
                     sym='', widths=0.4)
    set_box_color(bp4, 'c')
    plt.plot([], c='c', label='STTA')
    
    plt.legend()
    plt.xticks(range(0, len(ticks) * 2, 2), ticks)
    plt.xlim(-1, len(ticks)*2-1)
    plt.ylim(bottom=0)
    plt.xlabel("$N_a$")
    plt.ylabel("Consensus Steps")
    plt.legend(loc = 'lower left')
    plt.ylim(ymin=0)
    if init.monotonicity:
        plt.savefig(output_path_mono+'STTA_DSTA_steps_mono_vr_box.pdf')
    else:
        plt.savefig(output_path_non+'STTA_DSTA_steps_non_vr_box.pdf')
        
# =============================================================================
# 
# =============================================================================
print("----- plot_stta_vs_t3a_vr.py is loaded -----")
