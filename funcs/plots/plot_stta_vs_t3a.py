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
output_path_mono = os.path.abspath("output") + "/montecarlo/monotone/plot_stta_vs_t3a/"
# non-monotone
output_path_non = os.path.abspath("output") + "/montecarlo/non_monotone/plot_stta_vs_t3a/"

# error bar settings
#eb_fmt = "c.-"
#eb_markersize = 3
#eb_capsize = 2
#eb_capsize_line = 4
#eb_capsize_bar = 2
#eb_label = "STTA"
init.STTA.ls = 'c-+'


def plotSTTA_vs_T3A():
    # Draw figures
    print("\n======= Plot T3A vs STTA =======")
        
# total values             
    plt.figure("stta vs t3a utility")
    init.T3A.plotUtilities() 
    init.STTA.plotUtilities()
    # Plot figure
    plt.xlabel("$N_a$")
    plt.ylabel("Function Utility")
    plt.legend(loc = 'lower right')
    plt.ylim(bottom=0,top=120)
    if init.monotonicity:
        plt.savefig(output_path_mono+'STTA_vs_T3A_utility_mono.pdf')
    else:
        plt.savefig(output_path_non+'STTA_vs_T3A_utility_non.pdf')

# consensus steps             
    plt.figure("stta vs t3a steps")
    init.T3A.plotSteps()  
    init.STTA.plotSteps()
    # Plot figure
    plt.xlabel("$N_a$")
    plt.ylabel("Consensus Steps")
#    plt.legend(loc='center left',bbox_to_anchor=(0.76, 0.6))
    plt.ylim(bottom=0,top=35)
    if init.monotonicity:
        plt.legend(loc='lower left')
        plt.savefig(output_path_mono+'STTA_vs_T3A_steps_mono.pdf')
    else:
        plt.legend(loc='best')
        plt.savefig(output_path_non+'STTA_vs_T3A_steps_non.pdf')


# number of function evaluations         
    plt.figure("stta vs t3a evs")
    init.T3A.plotEvs()
    init.STTA.plotEvs()
    # Plot figure
    plt.xlabel("$N_a$")
    plt.ylabel("Function Evaluations (x$10^"+str(init.evs_scale)+"$)")
    plt.legend(loc = 'best')
    plt.ylim(bottom=0,top=7)
    if init.monotonicity:
        plt.savefig(output_path_mono+'STTA_vs_T3A_evs_mono.pdf')
    else:
        plt.savefig(output_path_non+'STTA_vs_T3A_evs_non.pdf')
        
    
    
# Plot ratio and error bar
    plt.figure("ratio_stta_vs_t3a")
    # Calculate the ratios compared with T3A. 
    Na_index_1 = 0
    Na_index_3 = len(init.STTA.utilities)-1
    Na_index_2 = round((Na_index_1+Na_index_3)/2)
    Na_indices = [Na_index_1, Na_index_2, Na_index_3]
    
    bar_groups = ['$N_a$='+str(init.Nas[Na_index_1]), 
                  '$N_a$='+str(init.Nas[Na_index_2]),
                  '$N_a$='+str(init.Nas[Na_index_3])]
    # set width of bar
    bar_width = 0.2
    # Set position of bar on X axis
    r1 = np.arange(len(bar_groups))
    r2 = [x + bar_width for x in r1]
    r3 = [x + bar_width for x in r2]
    # Get bar height
    bar_utilities = []
    bar_evs = []
    bar_steps = []
    for i in range(len(Na_indices)):
        bar_utilities.append(init.STTA.utilities[Na_indices[i]]
                            /init.T3A.utilities[Na_indices[i]]*100)
        bar_evs.append(init.STTA.evs[Na_indices[i]]
                            /init.T3A.evs[Na_indices[i]]*100)
        bar_steps.append(init.STTA.steps[Na_indices[i]]
                            /init.T3A.steps[Na_indices[i]]*100)
    # Plot bars
    rect_utility = plt.bar(r1, bar_utilities, color='r', 
                           width=bar_width, edgecolor='white', 
                           label='Function Utility')
    rect_ev = plt.bar(r2, bar_evs, color='g', 
                      width=bar_width, edgecolor='white', 
                      label='Running Time')
    rect_step = plt.bar(r3, bar_steps, color='b', 
                        width=bar_width, edgecolor='white', 
                        label='Consensus Steps')
    # Text for all bars
    for i in range(len(bar_groups)):
        plt.text(rect_utility[i].get_x()+rect_utility[i].get_width()/2., 
                 bar_utilities[i], str(round(bar_utilities[i],1)), 
                 ha='center', va='bottom')
        plt.text(rect_ev[i].get_x()+rect_utility[i].get_width()/2., 
                 bar_evs[i], str(round(bar_evs[i],1)), 
                 ha='center', va='bottom')
        plt.text(rect_step[i].get_x()+rect_utility[i].get_width()/2., 
                 bar_steps[i], str(round(bar_steps[i],1)), 
                 ha='center', va='bottom')
    
    # Add xticks on the middle of the group bars
    plt.xticks([r + bar_width for r in range(len(bar_groups))], bar_groups)
    plt.ylabel("Percent (%)")
    # Create legend & Show graphic
    max_height = max(max(bar_utilities),
                     max(bar_evs),
                     max(bar_steps))
#    plt.ylim(ymin=0, ymax=max_height*1.1)
    plt.ylim(ymin=0, ymax=120)
#    plt.legend(loc = 'best')
#    plt.legend(bbox_to_anchor=(0., 1.02, 1., .102),loc = 'best',ncol=3, mode="expand",borderaxespad=0.)
    plt.legend(loc = 'best',ncol=3, mode="expand",borderaxespad=0.2)
    if init.monotonicity:
        plt.savefig(output_path_mono+'STTA_vs_T3A_ratio_mono.pdf')
    else:
        plt.savefig(output_path_non+'STTA_vs_T3A_ratio_non.pdf')

# =============================================================================
# 
# =============================================================================
print("----- plot_stta_vs_t3a.py is loaded -----")
