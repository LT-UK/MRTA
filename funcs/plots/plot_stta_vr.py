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
output_path_mono = os.path.abspath("output") + "/variance/monotone/plot_stta_vr/"
# non-monotone
output_path_non = os.path.abspath("output") + "/variance/non_monotone/plot_stta_vr/"

# error bar settings
eb_fmt = "c.-"
eb_markersize = 3
eb_capsize = 2
eb_capsize_line = 4
eb_capsize_bar = 2
#eb_label = "STTA"

# plt.rcParams.update({
#     'font.size': 14,        # Default text size
#     'axes.titlesize': 16,   # Title font size
#     'axes.labelsize': 14,   # X/Y label font size
#     'xtick.labelsize': 12,  # X-tick labels
#     'ytick.labelsize': 12,  # Y-tick labels
#     'legend.fontsize': 12   # Legend text
# })

plot_size_text   = 16
plot_size_legend = 14
plot_size_axes_title = 14
plot_size_axes_label = 16
plot_size_tick = 16

plt.rcParams.update({
    'font.size': plot_size_text,        # Default text size
    'axes.titlesize': plot_size_axes_title,   # Title font size
    'axes.labelsize': plot_size_axes_label,   # X/Y label font size
    'xtick.labelsize': plot_size_tick,  # X-tick labels
    'ytick.labelsize': plot_size_tick,  # Y-tick labels
    'legend.fontsize': plot_size_legend   # Legend text
})

plt.tight_layout()


def plotSTTA_Vr():
    # Draw figures
    print("\n======= Plot SGA vs CBBA vs STTA vs LSTTA vs STBTA Variance =======")
        
# total values             
    plt.figure("stta utility variance")
    init.SGA.plotUtilities()
    init.CBBA.plotUtilities()  
    init.STTA.plotVariance_Utilities()
    init.LSTTA.plotVariance_Utilities()
    init.STBTA.plotVariance_Utilities()
    # Plot figure
    plt.xlabel("$N_a$")
    # plt.xlabel("$|\mathcal{A}|$")
    plt.ylabel("Function Utility")
    plt.legend(loc = 'lower right')
    plt.ylim(bottom=0)
    if init.monotonicity:
        plt.savefig(output_path_mono+'STTA_utility_mono_vr.pdf')
    else:
        plt.savefig(output_path_non+'STTA_utility_non_vr.pdf')

# consensus steps             
    plt.figure("stta steps variance")
    init.SGA.plotSteps()
    init.CBBA.plotSteps()  
    init.STTA.plotVariance_Steps()
    init.LSTTA.plotVariance_Steps()
    init.STBTA.plotVariance_Steps()
    # Plot figure
    plt.xlabel("$N_a$")
    plt.ylabel("Consensus Steps")
#    plt.legend(loc='center left',bbox_to_anchor=(0.76, 0.6))
    plt.ylim(bottom=0, top=init.Nt*1.1)
    if init.monotonicity:
        plt.legend(loc='center left',bbox_to_anchor=(0.76, 0.6))
        plt.savefig(output_path_mono+'STTA_steps_mono_vr.pdf')
    else:
        plt.legend(loc='best')
        plt.savefig(output_path_non+'STTA_steps_non_vr.pdf')


# number of function evaluations no CBBA           
    plt.figure("stta evs variance no CBBA")
#    init.SGA.plotEvs()
#    init.CBBA.plotEvs()  
    init.STTA.plotVariance_Evs()
    init.LSTTA.plotVariance_Evs()
    init.STBTA.plotVariance_Evs()
    # Plot figure
    plt.xlabel("$N_a$")
    plt.ylabel("Function Evaluations (x$10^"+str(init.evs_scale)+"$)")
    plt.legend(loc = 'best')
    plt.ylim(bottom=0)
    if init.monotonicity:
        plt.savefig(output_path_mono+'STTA_evs_mono_vr_noCBBA.pdf')
    else:
        plt.savefig(output_path_non+'STTA_evs_non_vr_noCBBA.pdf')
        
    # number of function evaluations            
    plt.figure("stta evs variance")
    init.SGA.plotEvs()
    init.CBBA.plotEvs()  
    init.STTA.plotVariance_Evs()
    init.LSTTA.plotVariance_Evs()
    init.STBTA.plotVariance_Evs()
    # Plot figure
    plt.xlabel("$N_a$")
    plt.ylabel("Function Evaluations (x$10^"+str(init.evs_scale)+"$)")
    plt.legend(loc = 'best')
    plt.ylim(bottom=0)
    if init.monotonicity:
        plt.savefig(output_path_mono+'STTA_evs_mono_vr.pdf')
    else:
        plt.savefig(output_path_non+'STTA_evs_non_vr.pdf')

# Plot Histograms
    if not (init.SGA.en and init.CBBA.en and init.STTA.en and init.LSTTA.en and init.STBTA.en):
        print("\n Error: SGA or CBBA or STTA or LSTTA or STBTA is not enabled!")
        return


# Function value bar chart   Na = Na_Start 
    plt.figure("ratio_STTA_4")
    bar_groups = ['CBBA', 'STTA', 'LSTTA', 'STBTA']
    # set width of bar
    bar_width = 0.25
    # Set position of bar on X axis
    r1 = np.arange(len(bar_groups))
    r2 = [x + bar_width for x in r1]
    r3 = [x + bar_width for x in r2]
#    r4 = [x + bar_width for x in r3]
    # Calculate the ratios compared with CBBA. Na = 10
    Na_index = 0
    # CBBA
    ratio_CBBA_utility = init.CBBA.utilities[Na_index]/init.SGA.utilities[Na_index]*100
    ratio_CBBA_ev = init.CBBA.evs[Na_index]/init.SGA.evs[Na_index]*100
    ratio_CBBA_step = init.CBBA.steps[Na_index]/init.SGA.steps[Na_index]*100
    # STTA
    # init.STTA.utilities_avg was calculated through init.STTA.plotVariance_Utilities()
    ratio_STTA_utility = init.STTA.utilities_avg[Na_index]/init.SGA.utilities[Na_index]*100
    ratio_STTA_ev = init.STTA.evs_avg[Na_index]/init.SGA.evs[Na_index]*100
    ratio_STTA_step = init.STTA.steps_avg[Na_index]/init.SGA.steps[Na_index]*100
    # LSTTA
    ratio_LSTTA_utility = init.LSTTA.utilities_avg[Na_index]/init.SGA.utilities[Na_index]*100
    ratio_LSTTA_ev = init.LSTTA.evs_avg[Na_index]/init.SGA.evs[Na_index]*100
    ratio_LSTTA_step = init.LSTTA.steps_avg[Na_index]/init.SGA.steps[Na_index]*100
    # STBTA
    ratio_STBTA_utility = init.STBTA.utilities_avg[Na_index]/init.SGA.utilities[Na_index]*100
    ratio_STBTA_ev = init.STBTA.evs_avg[Na_index]/init.SGA.evs[Na_index]*100
    ratio_STBTA_step = init.STBTA.steps_avg[Na_index]/init.SGA.steps[Na_index]*100
    # Set bar height
    bar_utility = [ratio_CBBA_utility, ratio_STTA_utility, ratio_LSTTA_utility, ratio_STBTA_utility]
    bar_ev = [ratio_CBBA_ev, ratio_STTA_ev, ratio_LSTTA_ev, ratio_STBTA_ev]
    bar_step = [ratio_CBBA_step, ratio_STTA_step, ratio_LSTTA_step, ratio_STBTA_step]     
    # Make the plot
    rect_utility = plt.bar(r1, bar_utility, color='r', width=bar_width, edgecolor='white', label='Function Value')
    rect_ev = plt.bar(r2, bar_ev, color='g', width=bar_width, edgecolor='white', label='Running Time')
    rect_step = plt.bar(r3, bar_step, color='b', width=bar_width, edgecolor='white', label='Consensus Steps')
    
    # Plot error bar        
    base_utility = 100/init.SGA.utilities[Na_index]
    base_evs = 100/init.SGA.evs[Na_index]
    base_steps = 100/init.SGA.steps[Na_index]
    # STTA
    STTA_ratio_err_utility = [[init.STTA.error_utilities[0][Na_index]*base_utility],
                             [init.STTA.error_utilities[1][Na_index]*base_utility]]                            
    plt.errorbar(r1[1],ratio_STTA_utility,yerr=STTA_ratio_err_utility,
                     fmt=eb_fmt,markersize=eb_markersize,capsize=eb_capsize)
    
    STTA_ratio_err_evs = [[init.STTA.error_evs[0][Na_index]*base_evs],
                             [init.STTA.error_evs[1][Na_index]*base_evs]]                           
    plt.errorbar(r2[1],ratio_STTA_ev,yerr=STTA_ratio_err_evs,
                     fmt=eb_fmt,markersize=eb_markersize,capsize=eb_capsize)
    
    STTA_ratio_err_steps = [[init.STTA.error_steps[0][Na_index]*base_steps],
                             [init.STTA.error_steps[1][Na_index]*base_steps]]                            
    plt.errorbar(r3[1],ratio_STTA_step,yerr=STTA_ratio_err_steps,
                     fmt=eb_fmt,markersize=eb_markersize,capsize=eb_capsize)    
    # LSTTA    
    LSTTA_ratio_err_utility = [[init.LSTTA.error_utilities[0][Na_index]*base_utility],
                             [init.LSTTA.error_utilities[1][Na_index]*base_utility]]                           
    plt.errorbar(r1[2],ratio_LSTTA_utility,yerr=LSTTA_ratio_err_utility,
                     fmt=eb_fmt,markersize=eb_markersize,capsize=eb_capsize)
    LSTTA_ratio_err_evs = [[init.LSTTA.error_evs[0][Na_index]*base_evs],
                             [init.LSTTA.error_evs[1][Na_index]*base_evs]]                          
    plt.errorbar(r2[2],ratio_LSTTA_ev,yerr=LSTTA_ratio_err_evs,
                     fmt=eb_fmt,markersize=eb_markersize,capsize=eb_capsize)
    LSTTA_ratio_err_steps = [[init.LSTTA.error_steps[0][Na_index]*base_steps],
                             [init.LSTTA.error_steps[1][Na_index]*base_steps]]                           
    plt.errorbar(r3[2],ratio_LSTTA_step,yerr=LSTTA_ratio_err_steps,
                     fmt=eb_fmt,markersize=eb_markersize,capsize=eb_capsize)    
    # STBTA
    STBTA_ratio_err_utility = [[init.STBTA.error_utilities[0][Na_index]*base_utility],
                             [init.STBTA.error_utilities[1][Na_index]*base_utility]]                           
    plt.errorbar(r1[3],ratio_STBTA_utility,yerr=STBTA_ratio_err_utility,
                     fmt=eb_fmt,markersize=eb_markersize,capsize=eb_capsize)
    STBTA_ratio_err_evs = [[init.STBTA.error_evs[0][Na_index]*base_evs],
                             [init.STBTA.error_evs[1][Na_index]*base_evs]]                          
    plt.errorbar(r2[3],ratio_STBTA_ev,yerr=STBTA_ratio_err_evs,
                     fmt=eb_fmt,markersize=eb_markersize,capsize=eb_capsize)
    STBTA_ratio_err_steps = [[init.STBTA.error_steps[0][Na_index]*base_steps],
                             [init.STBTA.error_steps[1][Na_index]*base_steps]]                           
    plt.errorbar(r3[3],ratio_STBTA_step,yerr=STBTA_ratio_err_steps,
                     fmt=eb_fmt,markersize=eb_markersize,capsize=eb_capsize)
    
    # Text for all bars
    for i in range(len(bar_groups)):
        plt.text(rect_utility[i].get_x()+rect_utility[i].get_width()/2., bar_utility[i], str(round(bar_utility[i],1)), ha='center', va='bottom')
        plt.text(rect_ev[i].get_x()+rect_utility[i].get_width()/2., bar_ev[i], str(round(bar_ev[i],1)), ha='center', va='bottom')
        plt.text(rect_step[i].get_x()+rect_utility[i].get_width()/2., bar_step[i], str(round(bar_step[i],1)), ha='center', va='bottom')
    # Add xticks on the middle of the group bars
    plt.xticks([r + bar_width for r in range(len(bar_groups))], bar_groups)
    plt.ylabel("Percent (%)")
    # Create legend & Show graphic
#    max_height = max(max(bar_utility),max(bar_ev),max(bar_step))
    max_height = STBTA_ratio_err_utility[1][0]+bar_utility[-1]
    plt.ylim(ymin=0, ymax=max_height*1.45)
    plt.legend(loc = 'best')
#    plt.legend(bbox_to_anchor=(0., 1.02, 1., .102),loc = 'best',ncol=3, mode="expand",borderaxespad=0.)
#    plt.legend(loc = 'best',ncol=3, mode="expand",borderaxespad=0.2)
#    plt.savefig('ratio_monotone_comparison_threshold_10.eps')
    if init.monotonicity:
        plt.savefig(output_path_mono+'STTA_ratio_mono_vr_4.pdf')
    else:
        plt.savefig(output_path_non+'STTA_ratio_non_vr_4.pdf')


# Function value bar chart   Na = Na_End 
    plt.figure("ratio_STTA_20")
    bar_groups = ['CBBA', 'STTA', 'LSTTA', 'STBTA']
    # set width of bar
    bar_width = 0.25
    # Set position of bar on X axis
    r1 = np.arange(len(bar_groups))
    r2 = [x + bar_width for x in r1]
    r3 = [x + bar_width for x in r2]
#    r4 = [x + bar_width for x in r3]
    # Calculate the ratios compared with CBBA. Na = 20
    Na_index = len(init.Nas)-1
    # CBBA
    ratio_CBBA_utility = init.CBBA.utilities[Na_index]/init.SGA.utilities[Na_index]*100
    ratio_CBBA_ev = init.CBBA.evs[Na_index]/init.SGA.evs[Na_index]*100
    ratio_CBBA_step = init.CBBA.steps[Na_index]/init.SGA.steps[Na_index]*100
    # STTA
    # init.STTA.utilities_avg was calculated through init.STTA.plotVariance_Utilities()
    ratio_STTA_utility = init.STTA.utilities_avg[Na_index]/init.SGA.utilities[Na_index]*100
    ratio_STTA_ev = init.STTA.evs_avg[Na_index]/init.SGA.evs[Na_index]*100
    ratio_STTA_step = init.STTA.steps_avg[Na_index]/init.SGA.steps[Na_index]*100
    # LSTTA
    ratio_LSTTA_utility = init.LSTTA.utilities_avg[Na_index]/init.SGA.utilities[Na_index]*100
    ratio_LSTTA_ev = init.LSTTA.evs_avg[Na_index]/init.SGA.evs[Na_index]*100
    ratio_LSTTA_step = init.LSTTA.steps_avg[Na_index]/init.SGA.steps[Na_index]*100
    # STBTA
    ratio_STBTA_utility = init.STBTA.utilities_avg[Na_index]/init.SGA.utilities[Na_index]*100
    ratio_STBTA_ev = init.STBTA.evs_avg[Na_index]/init.SGA.evs[Na_index]*100
    ratio_STBTA_step = init.STBTA.steps_avg[Na_index]/init.SGA.steps[Na_index]*100
    # Set bar height
    bar_utility = [ratio_CBBA_utility, ratio_STTA_utility, ratio_LSTTA_utility, ratio_STBTA_utility]
    bar_ev = [ratio_CBBA_ev, ratio_STTA_ev, ratio_LSTTA_ev, ratio_STBTA_ev]
    bar_step = [ratio_CBBA_step, ratio_STTA_step, ratio_LSTTA_step, ratio_STBTA_step]     
    # Make the plot
    rect_utility = plt.bar(r1, bar_utility, color='r', width=bar_width, edgecolor='white', label='Function Value')
    rect_ev = plt.bar(r2, bar_ev, color='g', width=bar_width, edgecolor='white', label='Running Time')
    rect_step = plt.bar(r3, bar_step, color='b', width=bar_width, edgecolor='white', label='Consensus Steps')
    
    # Plot error bar        
    base_utility = 100/init.SGA.utilities[Na_index]
    base_evs = 100/init.SGA.evs[Na_index]
    base_steps = 100/init.SGA.steps[Na_index]
    # STTA
    STTA_ratio_err_utility = [[init.STTA.error_utilities[0][Na_index]*base_utility],
                             [init.STTA.error_utilities[1][Na_index]*base_utility]]                            
    plt.errorbar(r1[1],ratio_STTA_utility,yerr=STTA_ratio_err_utility,
                     fmt=eb_fmt,markersize=eb_markersize,capsize=eb_capsize)
    
    STTA_ratio_err_evs = [[init.STTA.error_evs[0][Na_index]*base_evs],
                             [init.STTA.error_evs[1][Na_index]*base_evs]]                           
    plt.errorbar(r2[1],ratio_STTA_ev,yerr=STTA_ratio_err_evs,
                     fmt=eb_fmt,markersize=eb_markersize,capsize=eb_capsize)
    
    STTA_ratio_err_steps = [[init.STTA.error_steps[0][Na_index]*base_steps],
                             [init.STTA.error_steps[1][Na_index]*base_steps]]                            
    plt.errorbar(r3[1],ratio_STTA_step,yerr=STTA_ratio_err_steps,
                     fmt=eb_fmt,markersize=eb_markersize,capsize=eb_capsize)    
    # LSTTA    
    LSTTA_ratio_err_utility = [[init.LSTTA.error_utilities[0][Na_index]*base_utility],
                             [init.LSTTA.error_utilities[1][Na_index]*base_utility]]                           
    plt.errorbar(r1[2],ratio_LSTTA_utility,yerr=LSTTA_ratio_err_utility,
                     fmt=eb_fmt,markersize=eb_markersize,capsize=eb_capsize)
    LSTTA_ratio_err_evs = [[init.LSTTA.error_evs[0][Na_index]*base_evs],
                             [init.LSTTA.error_evs[1][Na_index]*base_evs]]                          
    plt.errorbar(r2[2],ratio_LSTTA_ev,yerr=LSTTA_ratio_err_evs,
                     fmt=eb_fmt,markersize=eb_markersize,capsize=eb_capsize)
    LSTTA_ratio_err_steps = [[init.LSTTA.error_steps[0][Na_index]*base_steps],
                             [init.LSTTA.error_steps[1][Na_index]*base_steps]]                           
    plt.errorbar(r3[2],ratio_LSTTA_step,yerr=LSTTA_ratio_err_steps,
                     fmt=eb_fmt,markersize=eb_markersize,capsize=eb_capsize)    
    # STBTA
    STBTA_ratio_err_utility = [[init.STBTA.error_utilities[0][Na_index]*base_utility],
                             [init.STBTA.error_utilities[1][Na_index]*base_utility]]                           
    plt.errorbar(r1[3],ratio_STBTA_utility,yerr=STBTA_ratio_err_utility,
                     fmt=eb_fmt,markersize=eb_markersize,capsize=eb_capsize)
    STBTA_ratio_err_evs = [[init.STBTA.error_evs[0][Na_index]*base_evs],
                             [init.STBTA.error_evs[1][Na_index]*base_evs]]                          
    plt.errorbar(r2[3],ratio_STBTA_ev,yerr=STBTA_ratio_err_evs,
                     fmt=eb_fmt,markersize=eb_markersize,capsize=eb_capsize)
    STBTA_ratio_err_steps = [[init.STBTA.error_steps[0][Na_index]*base_steps],
                             [init.STBTA.error_steps[1][Na_index]*base_steps]]                           
    plt.errorbar(r3[3],ratio_STBTA_step,yerr=STBTA_ratio_err_steps,
                     fmt=eb_fmt,markersize=eb_markersize,capsize=eb_capsize)
    
    # Text for all bars
    for i in range(len(bar_groups)):
        plt.text(rect_utility[i].get_x()+rect_utility[i].get_width()/2., bar_utility[i], str(round(bar_utility[i],1)), ha='center', va='bottom')
        plt.text(rect_ev[i].get_x()+rect_utility[i].get_width()/2., bar_ev[i], str(round(bar_ev[i],1)), ha='center', va='bottom')
        plt.text(rect_step[i].get_x()+rect_utility[i].get_width()/2., bar_step[i], str(round(bar_step[i],1)), ha='center', va='bottom')
    # Add xticks on the middle of the group bars
    plt.xticks([r + bar_width for r in range(len(bar_groups))], bar_groups)
    plt.ylabel("Percent (%)")
    # Create legend & Show graphic
    max_height = max(max(bar_utility),max(bar_ev),max(bar_step))
    plt.ylim(ymin=0, ymax=max_height*1.1)
    plt.legend(loc = 'best')
#    plt.legend(bbox_to_anchor=(0., 1.02, 1., .102),loc = 'best',ncol=3, mode="expand",borderaxespad=0.)
#    plt.legend(loc = 'best',ncol=3, mode="expand",borderaxespad=0.2)
#    plt.savefig('ratio_monotone_comparison_threshold_10.eps')
    if init.monotonicity:
        plt.savefig(output_path_mono+'STTA_ratio_mono_vr_20.pdf')
    else:
        plt.savefig(output_path_non+'STTA_ratio_non_vr_20.pdf')


# Break y axis        
# Function value bar chart   Na = 4  
    plt.figure("ratio_STTA_break_4")
    bar_groups = ['CBBA', 'STTA', 'LSTTA', 'STBTA']
    # set width of bar
    bar_width = 0.25
    # Set position of bar on X axis
    r1 = np.arange(len(bar_groups))
    r2 = [x + bar_width for x in r1]
    r3 = [x + bar_width for x in r2]
#    r4 = [x + bar_width for x in r3]
    # Calculate the ratios compared with CBBA. Na = 10
    Na_index = 0
    # CBBA
    ratio_CBBA_utility = init.CBBA.utilities[Na_index]/init.SGA.utilities[Na_index]*100
    ratio_CBBA_ev = init.CBBA.evs[Na_index]/init.SGA.evs[Na_index]*100
    ratio_CBBA_step = init.CBBA.steps[Na_index]/init.SGA.steps[Na_index]*100
    # STTA
    # init.STTA.utilities_avg was calculated through init.STTA.plotVariance_Utilities()
    ratio_STTA_utility = init.STTA.utilities_avg[Na_index]/init.SGA.utilities[Na_index]*100
    ratio_STTA_ev = init.STTA.evs_avg[Na_index]/init.SGA.evs[Na_index]*100
    ratio_STTA_step = init.STTA.steps_avg[Na_index]/init.SGA.steps[Na_index]*100
    # LSTTA
    ratio_LSTTA_utility = init.LSTTA.utilities_avg[Na_index]/init.SGA.utilities[Na_index]*100
    ratio_LSTTA_ev = init.LSTTA.evs_avg[Na_index]/init.SGA.evs[Na_index]*100
    ratio_LSTTA_step = init.LSTTA.steps_avg[Na_index]/init.SGA.steps[Na_index]*100
    # STBTA
    ratio_STBTA_utility = init.STBTA.utilities_avg[Na_index]/init.SGA.utilities[Na_index]*100
    ratio_STBTA_ev = init.STBTA.evs_avg[Na_index]/init.SGA.evs[Na_index]*100
    ratio_STBTA_step = init.STBTA.steps_avg[Na_index]/init.SGA.steps[Na_index]*100
    # Set bar height
    bar_utility = [ratio_CBBA_utility, ratio_STTA_utility, ratio_LSTTA_utility, ratio_STBTA_utility]
    bar_ev = [ratio_CBBA_ev, ratio_STTA_ev, ratio_LSTTA_ev, ratio_STBTA_ev]
    bar_step = [ratio_CBBA_step, ratio_STTA_step, ratio_LSTTA_step, ratio_STBTA_step]     
    
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
# Plot error bar        
    base_utility = 100/init.SGA.utilities[Na_index]
    base_evs = 100/init.SGA.evs[Na_index]
    base_steps = 100/init.SGA.steps[Na_index]
    # STTA
    STTA_ratio_err_utility = [[init.STTA.error_utilities[0][Na_index]*base_utility],
                             [init.STTA.error_utilities[1][Na_index]*base_utility]]                            
    plt.errorbar(r1[1],ratio_STTA_utility,yerr=STTA_ratio_err_utility,
                     fmt=eb_fmt,markersize=eb_markersize,capsize=eb_capsize)
    
    STTA_ratio_err_evs = [[init.STTA.error_evs[0][Na_index]*base_evs],
                             [init.STTA.error_evs[1][Na_index]*base_evs]]                           
    plt.errorbar(r2[1],ratio_STTA_ev,yerr=STTA_ratio_err_evs,
                     fmt=eb_fmt,markersize=eb_markersize,capsize=eb_capsize)
    
    STTA_ratio_err_steps = [[init.STTA.error_steps[0][Na_index]*base_steps],
                             [init.STTA.error_steps[1][Na_index]*base_steps]]                            
    plt.errorbar(r3[1],ratio_STTA_step,yerr=STTA_ratio_err_steps,
                     fmt=eb_fmt,markersize=eb_markersize,capsize=eb_capsize)    
    # LSTTA    
    LSTTA_ratio_err_utility = [[init.LSTTA.error_utilities[0][Na_index]*base_utility],
                             [init.LSTTA.error_utilities[1][Na_index]*base_utility]]                           
    plt.errorbar(r1[2],ratio_LSTTA_utility,yerr=LSTTA_ratio_err_utility,
                     fmt=eb_fmt,markersize=eb_markersize,capsize=eb_capsize)
    LSTTA_ratio_err_evs = [[init.LSTTA.error_evs[0][Na_index]*base_evs],
                             [init.LSTTA.error_evs[1][Na_index]*base_evs]]                          
    plt.errorbar(r2[2],ratio_LSTTA_ev,yerr=LSTTA_ratio_err_evs,
                     fmt=eb_fmt,markersize=eb_markersize,capsize=eb_capsize)
    LSTTA_ratio_err_steps = [[init.LSTTA.error_steps[0][Na_index]*base_steps],
                             [init.LSTTA.error_steps[1][Na_index]*base_steps]]                           
    plt.errorbar(r3[2],ratio_LSTTA_step,yerr=LSTTA_ratio_err_steps,
                     fmt=eb_fmt,markersize=eb_markersize,capsize=eb_capsize)    
    # STBTA
    STBTA_ratio_err_utility = [[init.STBTA.error_utilities[0][Na_index]*base_utility],
                             [init.STBTA.error_utilities[1][Na_index]*base_utility]]                           
    plt.errorbar(r1[3],ratio_STBTA_utility,yerr=STBTA_ratio_err_utility,
                     fmt=eb_fmt,markersize=eb_markersize,capsize=eb_capsize)
    STBTA_ratio_err_evs = [[init.STBTA.error_evs[0][Na_index]*base_evs],
                             [init.STBTA.error_evs[1][Na_index]*base_evs]]                          
    plt.errorbar(r2[3],ratio_STBTA_ev,yerr=STBTA_ratio_err_evs,
                     fmt=eb_fmt,markersize=eb_markersize,capsize=eb_capsize)
    STBTA_ratio_err_steps = [[init.STBTA.error_steps[0][Na_index]*base_steps],
                             [init.STBTA.error_steps[1][Na_index]*base_steps]]                           
    plt.errorbar(r3[3],ratio_STBTA_step,yerr=STBTA_ratio_err_steps,
                     fmt=eb_fmt,markersize=eb_markersize,capsize=eb_capsize)
    
    # Add xticks on the middle of the group bars
    plt.xticks([r + bar_width for r in range(len(bar_groups))], bar_groups)
    plt.ylabel("Percent (%)")
    # Create legend & Show graphic
    ax1.legend(loc = 'best')
    if init.monotonicity:
        plt.savefig(output_path_mono+'STTA_ratio_mono_break_4.pdf')
    else:
        plt.savefig(output_path_non+'STTA_ratio_non_break_4.pdf') 



# Break y axis        
# Function value bar chart   Na = 20 
    plt.figure("ratio_STTA_break_20")
    bar_groups = ['CBBA', 'STTA', 'LSTTA', 'STBTA']
    # set width of bar
    bar_width = 0.25
    # Set position of bar on X axis
    r1 = np.arange(len(bar_groups))
    r2 = [x + bar_width for x in r1]
    r3 = [x + bar_width for x in r2]
#    r4 = [x + bar_width for x in r3]
    # Calculate the ratios compared with CBBA. Na = 20
    Na_index = len(init.Nas)-1
    # CBBA
    ratio_CBBA_utility = init.CBBA.utilities[Na_index]/init.SGA.utilities[Na_index]*100
    ratio_CBBA_ev = init.CBBA.evs[Na_index]/init.SGA.evs[Na_index]*100
    ratio_CBBA_step = init.CBBA.steps[Na_index]/init.SGA.steps[Na_index]*100
    # STTA
    # init.STTA.utilities_avg was calculated through init.STTA.plotVariance_Utilities()
    ratio_STTA_utility = init.STTA.utilities_avg[Na_index]/init.SGA.utilities[Na_index]*100
    ratio_STTA_ev = init.STTA.evs_avg[Na_index]/init.SGA.evs[Na_index]*100
    ratio_STTA_step = init.STTA.steps_avg[Na_index]/init.SGA.steps[Na_index]*100
    # LSTTA
    ratio_LSTTA_utility = init.LSTTA.utilities_avg[Na_index]/init.SGA.utilities[Na_index]*100
    ratio_LSTTA_ev = init.LSTTA.evs_avg[Na_index]/init.SGA.evs[Na_index]*100
    ratio_LSTTA_step = init.LSTTA.steps_avg[Na_index]/init.SGA.steps[Na_index]*100
    # STBTA
    ratio_STBTA_utility = init.STBTA.utilities_avg[Na_index]/init.SGA.utilities[Na_index]*100
    ratio_STBTA_ev = init.STBTA.evs_avg[Na_index]/init.SGA.evs[Na_index]*100
    ratio_STBTA_step = init.STBTA.steps_avg[Na_index]/init.SGA.steps[Na_index]*100
    # Set bar height
    bar_utility = [ratio_CBBA_utility, ratio_STTA_utility, ratio_LSTTA_utility, ratio_STBTA_utility]
    bar_ev = [ratio_CBBA_ev, ratio_STTA_ev, ratio_LSTTA_ev, ratio_STBTA_ev]
    bar_step = [ratio_CBBA_step, ratio_STTA_step, ratio_LSTTA_step, ratio_STBTA_step]     
    
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
# Plot error bar        
    base_utility = 100/init.SGA.utilities[Na_index]
    base_evs = 100/init.SGA.evs[Na_index]
    base_steps = 100/init.SGA.steps[Na_index]
    # STTA
    STTA_ratio_err_utility = [[init.STTA.error_utilities[0][Na_index]*base_utility],
                             [init.STTA.error_utilities[1][Na_index]*base_utility]]                            
    plt.errorbar(r1[1],ratio_STTA_utility,yerr=STTA_ratio_err_utility,
                     fmt=eb_fmt,markersize=eb_markersize,capsize=eb_capsize)
    
    STTA_ratio_err_evs = [[init.STTA.error_evs[0][Na_index]*base_evs],
                             [init.STTA.error_evs[1][Na_index]*base_evs]]                           
    plt.errorbar(r2[1],ratio_STTA_ev,yerr=STTA_ratio_err_evs,
                     fmt=eb_fmt,markersize=eb_markersize,capsize=eb_capsize)
    
    STTA_ratio_err_steps = [[init.STTA.error_steps[0][Na_index]*base_steps],
                             [init.STTA.error_steps[1][Na_index]*base_steps]]                            
    plt.errorbar(r3[1],ratio_STTA_step,yerr=STTA_ratio_err_steps,
                     fmt=eb_fmt,markersize=eb_markersize,capsize=eb_capsize)    
    # LSTTA    
    LSTTA_ratio_err_utility = [[init.LSTTA.error_utilities[0][Na_index]*base_utility],
                             [init.LSTTA.error_utilities[1][Na_index]*base_utility]]                           
    plt.errorbar(r1[2],ratio_LSTTA_utility,yerr=LSTTA_ratio_err_utility,
                     fmt=eb_fmt,markersize=eb_markersize,capsize=eb_capsize)
    LSTTA_ratio_err_evs = [[init.LSTTA.error_evs[0][Na_index]*base_evs],
                             [init.LSTTA.error_evs[1][Na_index]*base_evs]]                          
    plt.errorbar(r2[2],ratio_LSTTA_ev,yerr=LSTTA_ratio_err_evs,
                     fmt=eb_fmt,markersize=eb_markersize,capsize=eb_capsize)
    LSTTA_ratio_err_steps = [[init.LSTTA.error_steps[0][Na_index]*base_steps],
                             [init.LSTTA.error_steps[1][Na_index]*base_steps]]                           
    plt.errorbar(r3[2],ratio_LSTTA_step,yerr=LSTTA_ratio_err_steps,
                     fmt=eb_fmt,markersize=eb_markersize,capsize=eb_capsize)    
    # STBTA
    STBTA_ratio_err_utility = [[init.STBTA.error_utilities[0][Na_index]*base_utility],
                             [init.STBTA.error_utilities[1][Na_index]*base_utility]]                           
    plt.errorbar(r1[3],ratio_STBTA_utility,yerr=STBTA_ratio_err_utility,
                     fmt=eb_fmt,markersize=eb_markersize,capsize=eb_capsize)
    STBTA_ratio_err_evs = [[init.STBTA.error_evs[0][Na_index]*base_evs],
                             [init.STBTA.error_evs[1][Na_index]*base_evs]]                          
    plt.errorbar(r2[3],ratio_STBTA_ev,yerr=STBTA_ratio_err_evs,
                     fmt=eb_fmt,markersize=eb_markersize,capsize=eb_capsize)
    STBTA_ratio_err_steps = [[init.STBTA.error_steps[0][Na_index]*base_steps],
                             [init.STBTA.error_steps[1][Na_index]*base_steps]]                           
    plt.errorbar(r3[3],ratio_STBTA_step,yerr=STBTA_ratio_err_steps,
                     fmt=eb_fmt,markersize=eb_markersize,capsize=eb_capsize)
    
    # Add xticks on the middle of the group bars
    plt.xticks([r + bar_width for r in range(len(bar_groups))], bar_groups)
    plt.ylabel("Percent (%)")
    # Create legend & Show graphic
    ax1.legend(loc = 'best')
    if init.monotonicity:
        plt.savefig(output_path_mono+'STTA_ratio_mono_break_20.pdf')
    else:
        plt.savefig(output_path_non+'STTA_ratio_non_break_20.pdf') 
        
        

def set_box_color(bp, color):
    plt.setp(bp['boxes'], color=color)
    plt.setp(bp['whiskers'], color=color)
    plt.setp(bp['caps'], color=color)
    plt.setp(bp['medians'], color=color)
    
    
def plotSTTA_Vr_Box():
    
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
    
    bp3 = plt.boxplot(init.STTA.utilities, 
                     positions=np.array(range(len(init.STTA.utilities)))*2.0, 
                     sym='', widths=0.3)
    set_box_color(bp3, 'c')
    plt.plot([], c='c', label='STTA')
    
#    bp4 = plt.boxplot(init.LSTTA.utilities, 
#                     positions=np.array(range(len(init.LSTTA.utilities)))*2.0, 
#                     sym='', widths=0.4)
#    set_box_color(bp4, 'b')
#    plt.plot([], c='b', label='LSTTA')
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
    # plt.xlabel("$N_a$")
    plt.xlabel("$|\mathcal{A}|$")
    plt.ylabel("Function Utility")
    plt.legend(loc = 'lower right')
    plt.ylim(ymin=0)
    plt.tight_layout(pad=0.2)
    if init.monotonicity:
        plt.savefig(output_path_mono+'STTA_utility_mono_vr_box.pdf')
    else:
        plt.savefig(output_path_non+'STTA_utility_non_vr_box.pdf')
        
        
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
    
    STTA_evs = [[init.STTA.evs[i][j]/10**init.evs_scale for j in range(len(init.STTA.evs[i]))] for i in range(len(init.STTA.evs))]
    bp3 = plt.boxplot(STTA_evs, 
                     positions=np.array(range(len(init.STTA.evs)))*2.0, 
                     sym='', widths=0.3)
    set_box_color(bp3, 'c')
    plt.plot([], c='c', label='STTA')
    
    plt.legend()
    plt.xticks(range(0, len(ticks) * 2, 2), ticks)
    plt.xlim(-1, len(ticks)*2-1)
    plt.ylim(bottom=0)
    # plt.xlabel("$N_a$")
    plt.xlabel("$|\mathcal{A}|$")
    plt.ylabel("Function Evaluations (x$10^"+str(init.evs_scale)+"$)")
    plt.legend(loc = 'upper left')
    plt.ylim(ymin=0)
    plt.tight_layout(pad=0.2)
    if init.monotonicity:
        plt.savefig(output_path_mono+'STTA_evs_mono_vr_box.pdf')
    else:
        plt.savefig(output_path_non+'STTA_evs_non_vr_box.pdf')
        
        
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
    
    bp3 = plt.boxplot(init.STTA.steps, 
                     positions=np.array(range(len(init.STTA.steps)))*2.0, 
                     sym='', widths=0.3)
    set_box_color(bp3, 'c')
    plt.plot([], c='c', label='STTA')
    
    plt.legend()
    plt.xticks(range(0, len(ticks) * 2, 2), ticks)
    plt.xlim(-1, len(ticks)*2-1)
    plt.ylim(bottom=0)
    # plt.xlabel("$N_a$")
    plt.xlabel("$|\mathcal{A}|$")
    plt.ylabel("Consensus Steps")
    # plt.legend(loc = 'lower left')
    plt.legend(loc = 'lower right')
    plt.ylim(ymin=0)
    plt.tight_layout(pad=0.2)
    if init.monotonicity:
        plt.savefig(output_path_mono+'STTA_steps_mono_vr_box.pdf')
    else:
        plt.savefig(output_path_non+'STTA_steps_non_vr_box.pdf')



# =============================================================================
# 
# =============================================================================
print("----- plot_stta_vr.py is loaded -----")
