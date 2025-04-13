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
output_path_mono = os.path.abspath("output") + "/variance/monotone/plot_stta_vs_t3a/"
# non-monotone
output_path_non = os.path.abspath("output") + "/variance/non_monotone/plot_stta_vs_t3a/"

# error bar settings
eb_fmt = "c.-"
eb_markersize = 3
eb_capsize = 2
eb_capsize_line = 4
eb_capsize_bar = 2
#eb_label = "STTA"


def plotSTTA_vs_T3A_Vr():
    # Draw figures
    print("\n======= Plot T3A vs STTA Variance =======")
        
# total values             
    plt.figure("stta vs t3a utility variance")
    init.T3A.plotUtilities() 
    init.STTA.plotVariance_Utilities()
    # Plot figure
    plt.xlabel("$N_a$")
    plt.ylabel("Function Utility")
    plt.legend(loc = 'lower right')
    plt.ylim(bottom=0)
    if init.monotonicity:
        plt.savefig(output_path_mono+'STTA_vs_T3A_utility_mono_vr.pdf')
    else:
        plt.savefig(output_path_non+'STTA_vs_T3A_utility_non_vr.pdf')

# consensus steps             
    plt.figure("stta vs t3a steps variance")
    init.T3A.plotSteps()  
    init.STTA.plotVariance_Steps()
    # Plot figure
    plt.xlabel("$N_a$")
    plt.ylabel("Consensus Steps")
#    plt.legend(loc='center left',bbox_to_anchor=(0.76, 0.6))
    plt.ylim(bottom=0)
    if init.monotonicity:
        plt.legend(loc='lower left')
        plt.savefig(output_path_mono+'STTA_vs_T3A_steps_mono_vr.pdf')
    else:
        plt.legend(loc='best')
        plt.savefig(output_path_non+'STTA_vs_T3A_steps_non_vr.pdf')


# number of function evaluations no CBBA           
    plt.figure("stta vs t3a evs variance")
    init.T3A.plotEvs()
    init.STTA.plotVariance_Evs()
    # Plot figure
    plt.xlabel("$N_a$")
    plt.ylabel("Function Evaluations (x$10^"+str(init.evs_scale)+"$)")
    plt.legend(loc = 'best')
    plt.ylim(bottom=0)
    if init.monotonicity:
        plt.savefig(output_path_mono+'STTA_vs_T3A_evs_mono_vr.pdf')
    else:
        plt.savefig(output_path_non+'STTA_vs_T3A_evs_non_vr.pdf')
        
    

# =============================================================================
# 
# =============================================================================
print("----- plot_stta_vs_t3a_vr.py is loaded -----")
