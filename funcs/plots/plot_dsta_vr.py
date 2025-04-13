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
import statistics as st

# =============================================================================
# 
# =============================================================================
# monotone
output_path_mono = os.path.abspath("output") + "/variance/monotone/plot_dsta_vr/"
# non-monotone
output_path_non = os.path.abspath("output") + "/variance/non_monotone/plot_dsta_vr/"

# error bar settings
eb_fmt = "c.-"
eb_markersize = 3
eb_capsize = 2
eb_capsize_line = 4
eb_capsize_bar = 2
eb_label = "DSTA"

def plotDSTA_Vr():
    # Draw figures
    print("\n======= Plot CBBA vs DSTA Variance =======")
        
    if not (init.DSTA.en and init.CBBA.en):
        print("\n Error: DSTA or CBBA is not enabled!")
        return
    
# total values             
    plt.figure("dsta utility variance")
    init.CBBA.plotUtilities()  
    DSTA_utilities_avg = [st.mean(init.DSTA.utilities[i]) 
                            for i in range(len(init.Nas))]
    upper_error = [max(init.DSTA.utilities[i])-DSTA_utilities_avg[i] 
                    for i in range(len(init.Nas))]
    lower_error = [-min(init.DSTA.utilities[i])+DSTA_utilities_avg[i] 
                    for i in range(len(init.Nas))]
    error_utilities = [lower_error, upper_error]
    plt.errorbar(init.Nas,
                 DSTA_utilities_avg,
                 yerr=error_utilities,
                 fmt=eb_fmt,
                 capsize=eb_capsize,
                 label = eb_label)
    # Plot figure
    plt.xlabel("$N_a$")
    plt.ylabel("Function Utility")
    plt.legend(loc = 'lower right')
    plt.ylim(bottom=0)
    if init.monotonicity:
        plt.savefig(output_path_mono+'utility_mono_vr.pdf')
    else:
        plt.savefig(output_path_non+'utility_non_vr.pdf')

# consensus steps             
    plt.figure("dsta steps variance")
    init.CBBA.plotSteps()  
    DSTA_steps_avg = [st.mean(init.DSTA.steps[i]) 
                        for i in range(len(init.Nas))]
    upper_error = [max(init.DSTA.steps[i])-DSTA_steps_avg[i] 
                    for i in range(len(init.Nas))]
    lower_error = [-min(init.DSTA.steps[i])+DSTA_steps_avg[i] 
                    for i in range(len(init.Nas))]
    error_steps = [lower_error, upper_error]
    plt.errorbar(init.Nas,
                 DSTA_steps_avg,
                 yerr=error_steps,
                 fmt=eb_fmt,
                 capsize=eb_capsize,
                 label = eb_label)
    # Plot figure
    plt.xlabel("$N_a$")
    plt.ylabel("Consensus Steps")
    plt.legend(loc = 'lower right')
    plt.ylim(bottom=0, top=init.Nt*1.1)
    if init.monotonicity:
        plt.savefig(output_path_mono+'steps_mono_vr.pdf')
    else:
        plt.savefig(output_path_non+'steps_non_vr.pdf')


# number of function evaluations            
    plt.figure("dsta evs variance")
    init.CBBA.plotEvs()  
    DSTA_evs_avg = [st.mean(init.DSTA.evs[i]) 
                    for i in range(len(init.Nas))]
    upper_error = [max(init.DSTA.evs[i]) - DSTA_evs_avg[i] 
                   for i in range(len(init.Nas))]
    lower_error = [-min(init.DSTA.evs[i]) + DSTA_evs_avg[i] 
                   for i in range(len(init.Nas))]
    error_evs = [lower_error, upper_error]    
    DSTA_evs_avg_scale = [evs/10**init.evs_scale for evs in DSTA_evs_avg]
    error_evs_scale = [[err/10**init.evs_scale for err in lower_error], 
                       [err/10**init.evs_scale for err in upper_error]]
    plt.errorbar(init.Nas,
                 DSTA_evs_avg_scale,
                 yerr=error_evs_scale,
                 fmt=eb_fmt,
                 capsize=eb_capsize,
                 label = eb_label)
    # Plot figure
    plt.xlabel("$N_a$")
    plt.ylabel("Function Evaluations (x$10^"+str(init.evs_scale)+"$)")
    plt.legend(loc = 'best')
    plt.ylim(bottom=0)
    if init.monotonicity:
        plt.savefig(output_path_mono+'evs_mono_vr.pdf')
    else:
        plt.savefig(output_path_non+'evs_non_vr.pdf')
        


# Plot ratio and error bar
    plt.figure("ratio_dsta_vr_errorbar")
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
        bar_utilities.append(DSTA_utilities_avg[Na_indices[i]]
                            /init.CBBA.utilities[Na_indices[i]]*100)
        bar_evs.append(DSTA_evs_avg[Na_indices[i]]
                            /init.CBBA.evs[Na_indices[i]]*100)
        bar_steps.append(DSTA_steps_avg[Na_indices[i]]
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
                 bar_utilities[i], str(round(bar_utilities[i],1)), 
                 ha='center', va='bottom')
        plt.text(rect_ev[i].get_x()+rect_utility[i].get_width()/2., 
                 bar_evs[i], str(round(bar_evs[i],1)), 
                 ha='center', va='bottom')
        plt.text(rect_step[i].get_x()+rect_utility[i].get_width()/2., 
                 bar_steps[i], str(round(bar_steps[i],1)), 
                 ha='center', va='bottom')
    max_err_h = 0
    # Sample ratio error bar
    for i in range(len(Na_indices)):
        # Get errors
        base_utility = 100/init.CBBA.utilities[Na_indices[i]]
        ratio_err_utility = [[error_utilities[0][Na_indices[i]]*base_utility],
                             [error_utilities[1][Na_indices[i]]*base_utility]]
        base_evs = 100/init.CBBA.evs[Na_indices[i]]
        ratio_err_evs = [[error_evs[0][Na_indices[i]]*base_evs],
                         [error_evs[1][Na_indices[i]]*base_evs]]
        base_steps = 100/init.CBBA.steps[Na_indices[i]]
        ratio_err_steps = [[error_steps[0][Na_indices[i]]*base_steps],
                           [error_steps[1][Na_indices[i]]*base_steps]]
        # Plot errorbar                            
        plt.errorbar(r1[i],bar_utilities[i],yerr=ratio_err_utility,
                     fmt=eb_fmt,markersize=eb_markersize,capsize=eb_capsize)
        plt.errorbar(r2[i],bar_evs[i],yerr=ratio_err_evs,
                     fmt=eb_fmt,markersize=eb_markersize,capsize=eb_capsize)
        plt.errorbar(r3[i],bar_steps[i],yerr=ratio_err_steps,
                     fmt=eb_fmt,markersize=eb_markersize,capsize=eb_capsize)
        max_err_temp = max(ratio_err_utility[1][0],
                           ratio_err_evs[1][0],
                           ratio_err_steps[1][0])
        if max_err_h < max_err_temp:
            max_err_h = max_err_temp
    
    # Add xticks on the middle of the group bars
    plt.xticks([r + bar_width for r in range(len(bar_groups))], bar_groups)
    plt.ylabel("Percent (%)")
    # Create legend & Show graphic
    max_height = max(max(bar_utilities),
                     max(bar_evs),
                     max(bar_steps)) + max_err_h
    plt.ylim(ymin=0, ymax=max_height*1.1)
    plt.legend(loc = 'best')
#    plt.legend(bbox_to_anchor=(0., 1.02, 1., .102),loc = 'best',ncol=3, mode="expand",borderaxespad=0.)
#    plt.legend(loc = 'best',ncol=3, mode="expand",borderaxespad=0.2)
    if init.monotonicity:
        plt.savefig(output_path_mono+'ratio_mono_vr_errorbar.pdf')
    else:
        plt.savefig(output_path_non+'ratio_non_vr_errorbar.pdf')
        
        
# Plot ratio
    plt.figure("ratio_dsta_vr")
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
        bar_utilities.append(DSTA_utilities_avg[Na_indices[i]]
                            /init.CBBA.utilities[Na_indices[i]]*100)
        bar_evs.append(DSTA_evs_avg[Na_indices[i]]
                            /init.CBBA.evs[Na_indices[i]]*100)
        bar_steps.append(DSTA_steps_avg[Na_indices[i]]
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
    plt.ylim(ymin=0, ymax=max_height*1.1)
    plt.legend(loc = 'best')
#    plt.legend(bbox_to_anchor=(0., 1.02, 1., .102),loc = 'best',ncol=3, mode="expand",borderaxespad=0.)
#    plt.legend(loc = 'best',ncol=3, mode="expand",borderaxespad=0.2)
    if init.monotonicity:
        plt.savefig(output_path_mono+'ratio_mono_vr.pdf')
    else:
        plt.savefig(output_path_non+'ratio_non_vr.pdf')



def set_box_color(bp, color):
    plt.setp(bp['boxes'], color=color)
    plt.setp(bp['whiskers'], color=color)
    plt.setp(bp['caps'], color=color)
    plt.setp(bp['medians'], color=color)


def plotDSTA_Vr_Box():
    
    # plot utilities
    plt.figure("boxplot_comparison_utilities")
    # ticks =  np.arange(0.04, 0.241, 0.04)
    ticks = [Na for Na in init.Nas]
    bp1 = plt.boxplot(init.DSTA.utilities, 
                     positions=np.array(range(len(init.DSTA.utilities)))*2.0, 
                     sym='', widths=0.4)
    set_box_color(bp1, 'c')
    plt.plot([], c='c', label='DSTA')
    CBBA_utilities = [[init.CBBA.utilities[i]] for i in range(len(init.CBBA.utilities))]
    bp2 = plt.boxplot(CBBA_utilities, 
                     positions=np.array(range(len(init.DSTA.utilities)))*2.0, 
                     sym='', widths=0.4)
    set_box_color(bp2, 'r')
    plt.plot([], c='r', label='CBBA')
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
        plt.savefig(output_path_mono+'DSTA_utility_mono_vr_box.pdf')
    else:
        plt.savefig(output_path_non+'DSTA_utility_non_vr_box.pdf')
    
    
    # plot steps
    plt.figure("boxplot_comparison_steps")
    # ticks =  np.arange(0.04, 0.241, 0.04)
    ticks = [Na for Na in init.Nas]
    bp1 = plt.boxplot(init.DSTA.steps, 
                     positions=np.array(range(len(init.DSTA.steps)))*2.0, 
                     sym='', widths=0.4)
    set_box_color(bp1, 'c')
    plt.plot([], c='c', label='DSTA')
    CBBA_steps = [[init.CBBA.steps[i]] for i in range(len(init.CBBA.steps))]
    bp2 = plt.boxplot(CBBA_steps, 
                     positions=np.array(range(len(init.DSTA.steps)))*2.0, 
                     sym='', widths=0.4)
    set_box_color(bp2, 'r')
    plt.plot([], c='r', label='CBBA')
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
    plt.ylabel("Consensus Steps")
    plt.legend(loc = 'lower right')
    plt.ylim(ymin=0,ymax=65)
    if init.monotonicity:
        plt.savefig(output_path_mono+'DSTA_steps_mono_vr_box.pdf')
    else:
        plt.savefig(output_path_non+'DSTA_steps_non_vr_box.pdf')
        
        
# plot dts
    plt.figure("boxplot_comparison_dts")
    # ticks =  np.arange(0.04, 0.241, 0.04)
    ticks = [Na for Na in init.Nas]
    bp1 = plt.boxplot(init.DSTA.dts, 
                     positions=np.array(range(len(init.DSTA.dts)))*2.0, 
                     sym='', widths=0.4)
    set_box_color(bp1, 'c')
    plt.plot([], c='c', label='DSTA')
    CBBA_dts = [[init.CBBA.dts[i]] for i in range(len(init.CBBA.dts))]
    bp2 = plt.boxplot(CBBA_dts, 
                     positions=np.array(range(len(init.DSTA.dts)))*2.0, 
                     sym='', widths=0.4)
    set_box_color(bp2, 'r')
    plt.plot([], c='r', label='CBBA')
    plt.legend()
    plt.xticks(range(0, len(ticks) * 2, 2), ticks)
    plt.xlim(-1, len(ticks)*2-1)
    plt.ylim(bottom=0)
    plt.xlabel("$N_a$")
    plt.ylabel("Running Time (sec)")
    plt.legend(loc = 'center right')
    plt.ylim(ymin=0)
    if init.monotonicity:
        plt.savefig(output_path_mono+'DSTA_dts_mono_vr_box.pdf')
    else:
        plt.savefig(output_path_non+'DSTA_dts_non_vr_box.pdf')
    
    
    # plot evs
    plt.figure("boxplot_comparison_evs")
    # ticks =  np.arange(0.04, 0.241, 0.04)
    ticks = [Na for Na in init.Nas]
    evs_DSTA_scaled = [[evs/10**init.evs_scale for evs in init.DSTA.evs[i]] for i in range(len(init.DSTA.evs))]
    bp1 = plt.boxplot(evs_DSTA_scaled, 
                     positions=np.array(range(len(init.DSTA.evs)))*2.0, 
                     sym='', widths=0.4)
    set_box_color(bp1, 'c')
    plt.plot([], c='c', label='DSTA')
    evs_CBBA_scaled = [[init.CBBA.evs[i]/10**init.evs_scale] for i in range(len(init.CBBA.evs))]
    bp2 = plt.boxplot(evs_CBBA_scaled, 
                     positions=np.array(range(len(init.DSTA.evs)))*2.0, 
                     sym='', widths=0.4)
    set_box_color(bp2, 'r')
    plt.plot([], c='r', label='CBBA')
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
    plt.ylabel("Function Evaluations (x$10^"+str(init.evs_scale)+"$)")
    plt.legend(loc = 'center right')
    plt.ylim(ymin=0)
    if init.monotonicity:
        plt.savefig(output_path_mono+'DSTA_evs_mono_vr_box.pdf')
    else:
        plt.savefig(output_path_non+'DSTA_evs_non_vr_box.pdf')
    
    



# =============================================================================
# 
# =============================================================================
print("----- plot_dsta_vr.py is loaded -----")
