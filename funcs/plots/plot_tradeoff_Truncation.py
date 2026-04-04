#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 16 12:08:04 2020

Plot tradeoffs related to the truncation position.

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
output_path_mono = os.path.abspath("output") + "/tradeoff/truncation/monotone/fig/"
# non-monotone
output_path_non = os.path.abspath("output") + "/tradeoff/truncation/non_monotone/fig/"

def plotTradeoff_Rho():
    # Draw figures
    print("\n========= Plot Tradeoff of Truncation =========")
    
# total values             
    plt.figure("tradeoff utility rho")
    init.TGTA.plotTradeoff_Rho_Utilities()
    init.LTGTA.plotTradeoff_Rho_Utilities()
    init.T3A.plotTradeoff_Rho_Utilities()
    init.LT3A.plotTradeoff_Rho_Utilities()
    init.TTBTA.plotTradeoff_Rho_Utilities()
    # Plot figure
    plt.xlabel(r"$\rho$")
    plt.ylabel("Function Utility")
    plt.legend(loc = 'lower right')
    plt.ylim(bottom=0, top=110)
    if init.monotonicity:
        plt.savefig(output_path_mono+'tradeoff_Rho_mono_utility.pdf')
    else:
        plt.savefig(output_path_non+'tradeoff_Rho_non_utility.pdf')
    

# number of consensus steps       
    plt.figure("tradeoff steps rho")
    init.TGTA.plotTradeoff_Rho_Steps()
    init.LTGTA.plotTradeoff_Rho_Steps()
    init.T3A.plotTradeoff_Rho_Steps()
    init.LT3A.plotTradeoff_Rho_Steps()
    init.TTBTA.plotTradeoff_Rho_Steps()
    # Plot figure
    plt.xlabel(r"$\rho$")
    plt.ylabel("Consensus Steps")
    plt.legend(loc = 'lower right')
    plt.ylim(bottom=0,top=30)
#    plt.ylim(bottom=0, top=init.Nt*1.1)
    if init.monotonicity:
        plt.savefig(output_path_mono+'tradeoff_Rho_mono_steps.pdf')
    else:
        plt.savefig(output_path_non+'tradeoff_Rho_non_steps.pdf')
    
    
#number of function evaluations     
    plt.figure("tradeoff evs rho")
    init.TGTA.plotTradeoff_Rho_Evs()
    init.LTGTA.plotTradeoff_Rho_Evs()
    init.T3A.plotTradeoff_Rho_Evs()
    init.LT3A.plotTradeoff_Rho_Evs()
    init.TTBTA.plotTradeoff_Rho_Evs()
    # Plot figure    
    plt.xlabel(r"$\rho$")
    plt.ylabel("Function Evaluations (x$10^"+str(init.evs_scale)+"$)")
    plt.legend(loc = 'lower right')
    plt.ylim(bottom=0,top=18)
    if init.monotonicity:
        plt.savefig(output_path_mono+'tradeoff_Rho_mono_evs.pdf')
    else:
        plt.savefig(output_path_non+'tradeoff_Rho_non_evs.pdf')
    


# =============================================================================
# 
# =============================================================================
print("----- plot_tradeoff_Truncation.py is loaded -----")
