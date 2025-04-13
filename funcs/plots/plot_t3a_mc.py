#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 16 12:08:04 2020

SGA vs CBBA vs T3A vs LT3A vs TTBTA

@author: Teng Li
lt.uk@outlook.com
United Kingdom
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
output_path_mono = os.path.abspath("output") + "/montecarlo/monotone/plot_t3a_mc/"
# non-monotone
output_path_non = os.path.abspath("output") + "/montecarlo/non_monotone/plot_t3a_mc/"

def plotT3A_MC():
    # Draw figures
    print("\n========= Plot SGA vs CBBA vs T3A vs LT3A vs TTBTA =========")
    
    if not (init.SGA.en and init.CBBA.en and init.T3A.en and init.LT3A.en and init.TTBTA.en):
        print("\n Error: SGA or CBBA or T3A or LT3A or TTBTA is not enabled!")
        return
    
# total values             
    plt.figure("t3a utility mc")
    init.SGA.plotUtilities()
    init.CBBA.plotUtilities() 
    init.T3A.plotUtilities() 
    init.LT3A.plotUtilities()
    init.TTBTA.plotUtilities()
    # Plot figure
    plt.xlabel("$N_a$")
    plt.ylabel("Function Utility")
    plt.legend(loc = 'lower right')
    plt.ylim(bottom=0)
    if init.monotonicity:
        plt.savefig(output_path_mono+'T3A_utility_mono_mc.pdf')
    else:
        plt.savefig(output_path_non+'T3A_utility_non_mc.pdf')
    

# number of consensus steps       
    plt.figure("t3a steps mc")
    init.SGA.plotSteps()
    init.CBBA.plotSteps() 
    init.T3A.plotSteps() 
    init.LT3A.plotSteps()
    init.TTBTA.plotSteps()
    # Plot figure
    plt.xlabel("$N_a$")
    plt.ylabel("Consensus Steps")
    plt.legend(loc = 'best')
    plt.ylim(bottom=0, top=init.Nt*1.1)
    if init.monotonicity:
        plt.savefig(output_path_mono+'T3A_steps_mono_mc.pdf')
    else:
        plt.savefig(output_path_non+'T3A_steps_non_mc.pdf')
    
    
#number of function evaluations     
    plt.figure("t3a evs mc")
    init.SGA.plotEvs()
    init.CBBA.plotEvs() 
    init.T3A.plotEvs() 
    init.LT3A.plotEvs()
    init.TTBTA.plotEvs()
    # Plot figure    
    plt.xlabel("$N_a$")
    plt.ylabel("Function Evaluations (x$10^"+str(init.evs_scale)+"$)")
    plt.legend(loc = 'best')
    plt.ylim(bottom=0)
    if init.monotonicity:
        plt.savefig(output_path_mono+'T3A_evs_mono_mc.pdf')
    else:
        plt.savefig(output_path_non+'T3A_evs_non_mc.pdf')
    

#number of function evaluations  no CBBA
    plt.figure("t3a evs mc no CBBA")
    init.SGA.plotEvs()
#    init.CBBA.plotEvs() 
    init.T3A.plotEvs() 
    init.LT3A.plotEvs()
    init.TTBTA.plotEvs()
    # Plot figure    
    plt.xlabel("$N_a$")
    plt.ylabel("Function Evaluations (x$10^"+str(init.evs_scale)+"$)")
    plt.legend(loc = 'best')
    plt.ylim(bottom=0)
    if init.monotonicity:
        plt.savefig(output_path_mono+'T3A_evs_mono_mc_noCBBA.pdf')
    else:
        plt.savefig(output_path_non+'T3A_evs_non_mc_noCBBA.pdf')


    # ratio -- comparison
#    print("\n==============================================")
# Function value bar chart   Na = Na_Start 
    plt.figure("ratio_T3A_10")
    bar_groups = ['CBBA', 'T3A', 'LT3A', 'TTBTA']
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
    # T3A
    ratio_T3A_utility = init.T3A.utilities[Na_index]/init.SGA.utilities[Na_index]*100
    ratio_T3A_ev = init.T3A.evs[Na_index]/init.SGA.evs[Na_index]*100
    ratio_T3A_step = init.T3A.steps[Na_index]/init.SGA.steps[Na_index]*100
    # LT3A
    ratio_LT3A_utility = init.LT3A.utilities[Na_index]/init.SGA.utilities[Na_index]*100
    ratio_LT3A_ev = init.LT3A.evs[Na_index]/init.SGA.evs[Na_index]*100
    ratio_LT3A_step = init.LT3A.steps[Na_index]/init.SGA.steps[Na_index]*100
    # TTBTA
    ratio_TTBTA_utility = init.TTBTA.utilities[Na_index]/init.SGA.utilities[Na_index]*100
    ratio_TTBTA_ev = init.TTBTA.evs[Na_index]/init.SGA.evs[Na_index]*100
    ratio_TTBTA_step = init.TTBTA.steps[Na_index]/init.SGA.steps[Na_index]*100
    # Set bar height
    bar_utility = [ratio_CBBA_utility, ratio_T3A_utility, ratio_LT3A_utility, ratio_TTBTA_utility]
    bar_ev = [ratio_CBBA_ev, ratio_T3A_ev, ratio_LT3A_ev, ratio_TTBTA_ev]
    bar_step = [ratio_CBBA_step, ratio_T3A_step, ratio_LT3A_step, ratio_TTBTA_step]     
    # Make the plot
    rect_utility = plt.bar(r1, bar_utility, color='r', width=bar_width, edgecolor='white', label='Function Value')
    rect_ev = plt.bar(r2, bar_ev, color='g', width=bar_width, edgecolor='white', label='Running Time')
    rect_step = plt.bar(r3, bar_step, color='b', width=bar_width, edgecolor='white', label='Consensus Steps')
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
        plt.savefig(output_path_mono+'T3A_ratio_mono_10.pdf')
    else:
        plt.savefig(output_path_non+'T3A_ratio_non_10.pdf')


# Function value bar chart   Na = 20  
    plt.figure("ratio_T3A_20")
    bar_groups = ['CBBA', 'T3A', 'LT3A', 'TTBTA']
    # set width of bar
    bar_width = 0.25
    # Set position of bar on X axis
    r1 = np.arange(len(bar_groups))
    r2 = [x + bar_width for x in r1]
    r3 = [x + bar_width for x in r2]
#    r4 = [x + bar_width for x in r3]
    # Calculate the ratios compared with CBBA. Na = 10
    Na_index = len(init.SGA.utilities) - 1
    # CBBA
    ratio_CBBA_utility = init.CBBA.utilities[Na_index]/init.SGA.utilities[Na_index]*100
    ratio_CBBA_ev = init.CBBA.evs[Na_index]/init.SGA.evs[Na_index]*100
    ratio_CBBA_step = init.CBBA.steps[Na_index]/init.SGA.steps[Na_index]*100
    # T3A
    ratio_T3A_utility = init.T3A.utilities[Na_index]/init.SGA.utilities[Na_index]*100
    ratio_T3A_ev = init.T3A.evs[Na_index]/init.SGA.evs[Na_index]*100
    ratio_T3A_step = init.T3A.steps[Na_index]/init.SGA.steps[Na_index]*100
    # LT3A
    ratio_LT3A_utility = init.LT3A.utilities[Na_index]/init.SGA.utilities[Na_index]*100
    ratio_LT3A_ev = init.LT3A.evs[Na_index]/init.SGA.evs[Na_index]*100
    ratio_LT3A_step = init.LT3A.steps[Na_index]/init.SGA.steps[Na_index]*100
    # TTBTA
    ratio_TTBTA_utility = init.TTBTA.utilities[Na_index]/init.SGA.utilities[Na_index]*100
    ratio_TTBTA_ev = init.TTBTA.evs[Na_index]/init.SGA.evs[Na_index]*100
    ratio_TTBTA_step = init.TTBTA.steps[Na_index]/init.SGA.steps[Na_index]*100
    # Set bar height
    bar_utility = [ratio_CBBA_utility, ratio_T3A_utility, ratio_LT3A_utility, ratio_TTBTA_utility]
    bar_ev = [ratio_CBBA_ev, ratio_T3A_ev, ratio_LT3A_ev, ratio_TTBTA_ev]
    bar_step = [ratio_CBBA_step, ratio_T3A_step, ratio_LT3A_step, ratio_TTBTA_step]     
    # Make the plot
    rect_utility = plt.bar(r1, bar_utility, color='r', width=bar_width, edgecolor='white', label='Function Value')
    rect_ev = plt.bar(r2, bar_ev, color='g', width=bar_width, edgecolor='white', label='Running Time')
    rect_step = plt.bar(r3, bar_step, color='b', width=bar_width, edgecolor='white', label='Consensus Steps')
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
        plt.savefig(output_path_mono+'T3A_ratio_mono_20.pdf')
    else:
        plt.savefig(output_path_non+'T3A_ratio_non_20.pdf')


# Break y axis        
# Function value bar chart   Na = 4  
    plt.figure("ratio_T3A_break_4")
    bar_groups = ['CBBA', 'T3A', 'LT3A', 'TTBTA']
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
    # T3A
    ratio_T3A_utility = init.T3A.utilities[Na_index]/init.SGA.utilities[Na_index]*100
    ratio_T3A_ev = init.T3A.evs[Na_index]/init.SGA.evs[Na_index]*100
    ratio_T3A_step = init.T3A.steps[Na_index]/init.SGA.steps[Na_index]*100
    # LT3A
    ratio_LT3A_utility = init.LT3A.utilities[Na_index]/init.SGA.utilities[Na_index]*100
    ratio_LT3A_ev = init.LT3A.evs[Na_index]/init.SGA.evs[Na_index]*100
    ratio_LT3A_step = init.LT3A.steps[Na_index]/init.SGA.steps[Na_index]*100
    # TTBTA
    ratio_TTBTA_utility = init.TTBTA.utilities[Na_index]/init.SGA.utilities[Na_index]*100
    ratio_TTBTA_ev = init.TTBTA.evs[Na_index]/init.SGA.evs[Na_index]*100
    ratio_TTBTA_step = init.TTBTA.steps[Na_index]/init.SGA.steps[Na_index]*100
    # Set bar height
    bar_utility = [ratio_CBBA_utility, ratio_T3A_utility, ratio_LT3A_utility, ratio_TTBTA_utility]
    bar_ev = [ratio_CBBA_ev, ratio_T3A_ev, ratio_LT3A_ev, ratio_TTBTA_ev]
    bar_step = [ratio_CBBA_step, ratio_T3A_step, ratio_LT3A_step, ratio_TTBTA_step]     
    
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
    # Create legend & Show graphic
#    max_height = max(max(bar_utility),max(bar_ev),max(bar_step))
#    plt.ylim(ymin=0, ymax=max_height*1.1)
    ax1.legend(loc = 'best')
#    plt.legend(bbox_to_anchor=(0., 1.02, 1., .102),loc = 'best',ncol=3, mode="expand",borderaxespad=0.)
#    plt.legend(loc = 'best',ncol=3, mode="expand",borderaxespad=0.2)
#    plt.savefig('ratio_monotone_comparison_threshold_10.eps')
    if init.monotonicity:
        plt.savefig(output_path_mono+'T3A_ratio_mono_break_4.pdf')
    else:
        plt.savefig(output_path_non+'T3A_ratio_non_break_4.pdf')  


        
# Break y axis        
# Function value bar chart   Na = 20  
    plt.figure("ratio_T3A_break_20")
    bar_groups = ['CBBA', 'T3A', 'LT3A', 'TTBTA']
    # set width of bar
    bar_width = 0.25
    # Set position of bar on X axis
    r1 = np.arange(len(bar_groups))
    r2 = [x + bar_width for x in r1]
    r3 = [x + bar_width for x in r2]
#    r4 = [x + bar_width for x in r3]
    # Calculate the ratios compared with CBBA. Na = 10
    Na_index = len(init.SGA.utilities) - 1
    # CBBA
    ratio_CBBA_utility = init.CBBA.utilities[Na_index]/init.SGA.utilities[Na_index]*100
    ratio_CBBA_ev = init.CBBA.evs[Na_index]/init.SGA.evs[Na_index]*100
    ratio_CBBA_step = init.CBBA.steps[Na_index]/init.SGA.steps[Na_index]*100
    # T3A
    ratio_T3A_utility = init.T3A.utilities[Na_index]/init.SGA.utilities[Na_index]*100
    ratio_T3A_ev = init.T3A.evs[Na_index]/init.SGA.evs[Na_index]*100
    ratio_T3A_step = init.T3A.steps[Na_index]/init.SGA.steps[Na_index]*100
    # LT3A
    ratio_LT3A_utility = init.LT3A.utilities[Na_index]/init.SGA.utilities[Na_index]*100
    ratio_LT3A_ev = init.LT3A.evs[Na_index]/init.SGA.evs[Na_index]*100
    ratio_LT3A_step = init.LT3A.steps[Na_index]/init.SGA.steps[Na_index]*100
    # TTBTA
    ratio_TTBTA_utility = init.TTBTA.utilities[Na_index]/init.SGA.utilities[Na_index]*100
    ratio_TTBTA_ev = init.TTBTA.evs[Na_index]/init.SGA.evs[Na_index]*100
    ratio_TTBTA_step = init.TTBTA.steps[Na_index]/init.SGA.steps[Na_index]*100
    # Set bar height
    bar_utility = [ratio_CBBA_utility, ratio_T3A_utility, ratio_LT3A_utility, ratio_TTBTA_utility]
    bar_ev = [ratio_CBBA_ev, ratio_T3A_ev, ratio_LT3A_ev, ratio_TTBTA_ev]
    bar_step = [ratio_CBBA_step, ratio_T3A_step, ratio_LT3A_step, ratio_TTBTA_step]     
    
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
    # Create legend & Show graphic
#    max_height = max(max(bar_utility),max(bar_ev),max(bar_step))
#    plt.ylim(ymin=0, ymax=max_height*1.1)
    ax1.legend(loc = 'best')
#    plt.legend(bbox_to_anchor=(0., 1.02, 1., .102),loc = 'best',ncol=3, mode="expand",borderaxespad=0.)
#    plt.legend(loc = 'best',ncol=3, mode="expand",borderaxespad=0.2)
#    plt.savefig('ratio_monotone_comparison_threshold_10.eps')
    if init.monotonicity:
        plt.savefig(output_path_mono+'T3A_ratio_mono_break_20.pdf')
    else:
        plt.savefig(output_path_non+'T3A_ratio_non_break_20.pdf')    
# =============================================================================
# 
# =============================================================================
print("----- plot_t3a_mc.py is loaded -----")
