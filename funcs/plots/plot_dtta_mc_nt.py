#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
DTTA vs Benchmarks

@author: Teng Li
teng.li@cranfield.ac.uk; lt.uk@outlook.com
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

#init.DTTA.ls = 'g.-'
#init.LDTTA.ls = 'g--x'
#init.TBTA.ls = 'g-.+'

# monotone
output_path_mono = os.path.abspath("output") + "/montecarlo_nt/monotone/plot_dtta_mc_nt/"
# non-monotone
output_path_non = os.path.abspath("output") + "/montecarlo_nt/non_monotone/plot_dtta_mc_nt/"

def plotDTTA_MC_Nt(col_index=0):
    # Draw figures
    print("\n========= Bar plots =========")
    
    # Break y axis        
    # Function value bar chart   Na = 20  
    plt.figure("ratio_DTTA_break")
    bar_groups = ['CBBA', 'DSTA', 'DTTA']
    # set width of bar
    bar_width = 0.25
    # Set position of bar on X axis
    r1 = np.arange(len(bar_groups))
    r2 = [x + bar_width for x in r1]
    r3 = [x + bar_width for x in r2]
#    r4 = [x + bar_width for x in r3]
    # Calculate the ratios compared with CBBA. Na = 10
    # col_index = len(init.SGA.utilities) - 1
    # CBBA
    ratio_CBBA_utility = init.CBBA.utilities[col_index]/init.SGA.utilities[col_index]*100
    ratio_CBBA_ev = init.CBBA.evs[col_index]/init.SGA.evs[col_index]*100
    ratio_CBBA_step = init.CBBA.steps[col_index]/init.SGA.steps[col_index]*100
    # DSTA
    ratio_DSTA_utility = init.DSTA.utilities[col_index]/init.SGA.utilities[col_index]*100
    ratio_DSTA_ev = init.DSTA.evs[col_index]/init.SGA.evs[col_index]*100
    ratio_DSTA_step = init.DSTA.steps[col_index]/init.SGA.steps[col_index]*100
    # DTTA
    ratio_DTTA_utility = init.DTTA.utilities[col_index]/init.SGA.utilities[col_index]*100
    ratio_DTTA_ev = init.DTTA.evs[col_index]/init.SGA.evs[col_index]*100
    ratio_DTTA_step = init.DTTA.steps[col_index]/init.SGA.steps[col_index]*100
    # Set bar height
    bar_utility = [ratio_CBBA_utility, ratio_DSTA_utility, ratio_DTTA_utility]
    bar_ev = [ratio_CBBA_ev, ratio_DSTA_ev, ratio_DTTA_ev]
    bar_step = [ratio_CBBA_step, ratio_DSTA_step, ratio_DTTA_step]     
    
    # set height_ratios or width_ratios
    f, (ax1, ax2) = plt.subplots(2, 1, sharex=True, gridspec_kw={'height_ratios': [1, 2]})
    max_height = max(max(bar_utility),max(bar_ev),max(bar_step))
    ax1.set_ylim([max_height-40, max_height+20])
    ax2.set_ylim([0, 120])
    # Make the plot
    rect_utility = ax1.bar(r1, bar_utility, color='r', width=bar_width, edgecolor='white', label='Function Value')
    rect_ev = ax1.bar(r2, bar_ev, color='g', width=bar_width, edgecolor='white', label='Running Time')
    rect_step = ax1.bar(r3, bar_step, color='b', width=bar_width, edgecolor='white', label='Consensus Steps')
    
    rect_utility = ax2.bar(r1, bar_utility, color='r', width=bar_width, edgecolor='white', label='Function Value')
    rect_ev = ax2.bar(r2, bar_ev, color='g', width=bar_width, edgecolor='white', label='Running Time')
    rect_step = ax2.bar(r3, bar_step, color='b', width=bar_width, edgecolor='white', label='Consensus Steps')
    
    # hide the spines between ax and ax2
    ax1.spines['bottom'].set_visible(False)
    ax2.spines['top'].set_visible(False)
    ax1.xaxis.tick_top()
    ax1.tick_params(labeltop=False)  # don't put tick labels at the top
    ax2.xaxis.tick_bottom()
    
    d = .015  # how big to make the diagonal lines in axes coordinates
    # arguments to pass to plot, just so we don't keep repeating them
    kwargs = dict(transform=ax1.transAxes, color='k', clip_on=False)
    ax1.plot((-d, +d), (-d, +d), **kwargs)        # top-left diagonal
    ax1.plot((1 - d, 1 + d), (-d, +d), **kwargs)  # top-right diagonal
    
    kwargs.update(transform=ax2.transAxes)  # switch to the bottom axes
    ax2.plot((-d, +d), (1 - d, 1 + d), **kwargs)  # bottom-left diagonal
    ax2.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)  # bottom-right diagonal
    
    # Text for all bars
    for i in range(len(bar_groups)):
        plt.text(rect_utility[i].get_x()+rect_utility[i].get_width()/2., 
                  bar_utility[i], 
                  str(round(bar_utility[i],1)), 
                  ha='center', va='bottom')
        if i>0:
            plt.text(rect_ev[i].get_x()+rect_utility[i].get_width()/2., 
                      bar_ev[i], 
                      str(round(bar_ev[i],1)), 
                      ha='center', va='bottom')
        plt.text(rect_step[i].get_x()+rect_utility[i].get_width()/2., 
                  bar_step[i], 
                  str(round(bar_step[i],1)), 
                  ha='center', va='bottom')
    ax1.text(rect_ev[0].get_x()+rect_utility[0].get_width()/2., 
                      bar_ev[0], 
                      str(round(bar_ev[0],1)), 
                      ha='center', va='bottom')
    # Add xticks on the middle of the group bars
    plt.xticks([r + bar_width for r in range(len(bar_groups))], bar_groups)
    plt.ylabel("Percent (%)")

    ax1.legend(loc = 'best')

    if init.monotonicity:
        plt.savefig(output_path_mono+'DTTA_ratio_mono_break_'+str(col_index)+'.pdf')
    else:
        plt.savefig(output_path_non+'DTTA_ratio_non_break_'+str(col_index)+'.pdf')  

# =============================================================================
# 
# =============================================================================
print("----- plot_dtta_mc_nt.py is loaded -----")
