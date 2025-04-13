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
output_path_mono = os.path.abspath("output") + "/tradeoff/epsilon/monotone/fig/"
# non-monotone
output_path_non = os.path.abspath("output") + "/tradeoff/epsilon/non_monotone/fig/"

def plotTradeoff_Eps():
    # Draw figures
    print("\n========= Plot Tradeoff of Threshold =========")
    
# total values             
    plt.figure("tradeoff utility eps")
    init.DTTA.plotTradeoff_Eps_Utilities()
    init.LDTTA.plotTradeoff_Eps_Utilities()
    init.TBTA.plotTradeoff_Eps_Utilities()
    init.T3A.plotTradeoff_Eps_Utilities()
    init.LT3A.plotTradeoff_Eps_Utilities()
    init.TTBTA.plotTradeoff_Eps_Utilities()
    init.STTA.plotTradeoff_Eps_Utilities()
    init.LSTTA.plotTradeoff_Eps_Utilities()
    init.STBTA.plotTradeoff_Eps_Utilities()
    # Plot figure
    plt.xlabel("$\epsilon$")
    plt.ylabel("Function Utility")
    plt.legend(loc = 'lower right')
    plt.ylim(bottom=0, top=110)
    if init.monotonicity:
        plt.savefig(output_path_mono+'tradeoff_Eps_mono_utility.pdf')
    else:
        plt.savefig(output_path_non+'tradeoff_Eps_non_utility.pdf')
    

# number of consensus steps       
    plt.figure("tradeoff steps eps")
    init.DTTA.plotTradeoff_Eps_Steps()
    init.LDTTA.plotTradeoff_Eps_Steps()
    init.TBTA.plotTradeoff_Eps_Steps()
    init.T3A.plotTradeoff_Eps_Steps()
    init.LT3A.plotTradeoff_Eps_Steps()
    init.TTBTA.plotTradeoff_Eps_Steps()
    init.STTA.plotTradeoff_Eps_Steps()
    init.LSTTA.plotTradeoff_Eps_Steps()
    init.STBTA.plotTradeoff_Eps_Steps()
    # Plot figure
    plt.xlabel("$\epsilon$")
    plt.ylabel("Consensus Steps")
    plt.legend(loc = 'best')
    # plt.ylim(bottom=0)
    plt.ylim(bottom=0, top=50)
    if init.monotonicity:
        plt.savefig(output_path_mono+'tradeoff_Eps_mono_steps.pdf')
    else:
        plt.savefig(output_path_non+'tradeoff_Eps_non_steps.pdf')
    
    
#number of function evaluations     
    plt.figure("tradeoff evs eps")
    init.DTTA.plotTradeoff_Eps_Evs()
    init.LDTTA.plotTradeoff_Eps_Evs()
    init.TBTA.plotTradeoff_Eps_Evs()
    init.T3A.plotTradeoff_Eps_Evs()
    init.LT3A.plotTradeoff_Eps_Evs()
    init.TTBTA.plotTradeoff_Eps_Evs()
    init.STTA.plotTradeoff_Eps_Evs()
    init.LSTTA.plotTradeoff_Eps_Evs()
    init.STBTA.plotTradeoff_Eps_Evs()
    # Plot figure    
    plt.xlabel("$\epsilon$")
    plt.ylabel("Function Evaluations (x$10^"+str(init.evs_scale)+"$)")
    plt.legend(loc = 'best')
    plt.ylim(bottom=0)
    plt.ylim(bottom=0, top=2.5)
    if init.monotonicity:
        plt.savefig(output_path_mono+'tradeoff_Eps_mono_evs.pdf')
    else:
        plt.savefig(output_path_non+'tradeoff_Eps_non_evs.pdf')
    


# =============================================================================
# 
# =============================================================================
print("----- plot_tradeoff_Threshold.py is loaded -----")
