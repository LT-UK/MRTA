#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 18 15:58:04 2020

DSTA vs CBBA

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
output_path_mono = os.path.abspath("output") + "/montecarlo/monotone/plot_dsta_mc/"
# non-monotone
output_path_non = os.path.abspath("output") + "/montecarlo/non_monotone/plot_dsta_mc/"

def plotDSTA_MC():
    # Draw figures
    print("\n========= Plot DSTA vs CBBA =========")
    
    if not (init.DSTA.en and init.CBBA.en):
        print("\n Error: DSTA or CBBA is not enabled!")
        return
    
# total values             
    plt.figure("dsta utility mc")
    init.CBBA.plotUtilities() 
    init.DSTA.plotUtilities() 
    # Plot figure
    plt.xlabel("$N_a$")
    plt.ylabel("Function Utility")
    plt.legend(loc = 'lower right')
    plt.ylim(bottom=0)
    if init.monotonicity:
        plt.savefig(output_path_mono+'utility_mono_mc.pdf')
    else:
        plt.savefig(output_path_non+'utility_non_mc.pdf')
    

# number of consensus steps       
    plt.figure("dsta steps mc")
    init.CBBA.plotSteps() 
    init.DSTA.plotSteps()
    # Plot figure
    plt.xlabel("$N_a$")
    plt.ylabel("Consensus Steps")
    plt.legend(loc = 'best')
    plt.ylim(bottom=0, top=init.Nt*1.1)
    if init.monotonicity:
        plt.savefig(output_path_mono+'steps_mono_mc.pdf')
    else:
        plt.savefig(output_path_non+'steps_non_mc.pdf')
    
    
#number of function evaluations     
    plt.figure("dsta evs mc")
    init.CBBA.plotEvs() 
    init.DSTA.plotEvs()
    # Plot figure    
    plt.xlabel("$N_a$")
    plt.ylabel("Function Evaluations (x$10^"+str(init.evs_scale)+"$)")
    plt.legend(loc = 'best')
    plt.ylim(bottom=0)
    if init.monotonicity:
        plt.savefig(output_path_mono+'evs_mono_mc.pdf')
    else:
        plt.savefig(output_path_non+'evs_non_mc.pdf')
    

# ratio    
    plt.figure("ratio_dsta mc")
    # Calculate the ratios compared with CBBA. 
    Na_index_1 = 0
    Na_index_3 = len(init.DSTA.utilities)-1
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
        bar_utilities.append(init.DSTA.utilities[Na_indices[i]]
                            /init.CBBA.utilities[Na_indices[i]]*100)
        bar_evs.append(init.DSTA.evs[Na_indices[i]]
                            /init.CBBA.evs[Na_indices[i]]*100)
        bar_steps.append(init.DSTA.steps[Na_indices[i]]
                            /init.CBBA.steps[Na_indices[i]]*100)
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
                 bar_utilities[i], 
                 str(round(bar_utilities[i],1)), 
                 ha='center', va='bottom')
        plt.text(rect_ev[i].get_x()+rect_utility[i].get_width()/2., 
                 bar_evs[i], 
                 str(round(bar_evs[i],1)), 
                 ha='center', va='bottom')
        plt.text(rect_step[i].get_x()+rect_utility[i].get_width()/2., 
                 bar_steps[i], 
                 str(round(bar_steps[i],1)), 
                 ha='center', va='bottom')
    # Add xticks on the middle of the group bars
    plt.xticks([r + bar_width for r in range(len(bar_groups))], bar_groups)
    plt.ylabel("Percent (%)")
    # Create legend & Show graphic
    max_height = max(max(bar_utilities),max(bar_evs),max(bar_steps))
    plt.ylim(ymin=0, ymax=max_height*1.1)
    plt.legend(loc = 'best')
#    plt.legend(bbox_to_anchor=(0., 1.02, 1., .102),loc = 'best',ncol=3, mode="expand",borderaxespad=0.)
#    plt.legend(loc = 'best',ncol=3, mode="expand",borderaxespad=0.2)
    if init.monotonicity:
        plt.savefig(output_path_mono+'ratio_mono_mc.pdf')
    else:
        plt.savefig(output_path_non+'ratio_non_mc.pdf')
 

    
# =============================================================================
# 
# =============================================================================
print("----- plot_dsta_mc.py is loaded -----")
