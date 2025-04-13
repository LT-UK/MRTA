#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 18 15:58:04 2020

Plot figures for all algorithms.

@author: Teng Li
lt.uk@outlook.com
United Kingdom
All Rights Reserved
"""

import matplotlib.pyplot as plt
import os # find path
import init

# =============================================================================
# 
# =============================================================================
if init.monotonicity:
    output_path = os.path.abspath("output") + "/montecarlo/monotone/plotall/"
else:
    output_path = os.path.abspath("output") + "/montecarlo/non_monotone/plotall/"

def plotAll(benchmark=1, original=1, bundle=1, lazy=1):
    # Draw figures
    print("\n=========== Plot All ===========")
    
# total values             
    plt.figure("utilities")
    if benchmark:
        init.GA.plotUtilities()
        init.AuctionSequential.plotUtilities()
        init.AuctionParallel.plotUtilities()
        init.AuctionGPrim.plotUtilities()
        init.SGA.plotUtilities() 
        init.CBBA.plotUtilities() 
    if original:
        init.TGTA.plotUtilities() 
        init.DTTA.plotUtilities() 
        init.T3A.plotUtilities() 
        init.DSTA.plotUtilities() 
        init.STTA.plotUtilities() 
    if bundle:
        init.TBTA.plotUtilities() 
        init.TTBTA.plotUtilities()
        init.STBTA.plotUtilities()
    if lazy:
        init.LSGA.plotUtilities() 
        init.LTGTA.plotUtilities()
        init.LDTTA.plotUtilities()
        init.LT3A.plotUtilities()
        init.LSTA.plotUtilities()
        init.LSTTA.plotUtilities()
    # Plot figure
    plt.xlabel("$N_a$")
    plt.ylabel("Function Utility")
    plt.legend(loc = 'lower right')
    plt.ylim(bottom=0)
    if original and not lazy:
        plt.savefig(output_path+'/original/utility_monotone_cmp_original.pdf')
    elif lazy and not original:
        plt.savefig(output_path+'/lazy/utility_monotone_cmp_lazy.pdf')
    else:
        plt.savefig(output_path+'/all/utility_monotone_cmp_all.pdf')
    
    
# total running time /sec      
    plt.figure("time")
    if benchmark:
        init.GA.plotDts()
        init.AuctionSequential.plotDts()
        init.AuctionParallel.plotDts()
        init.AuctionGPrim.plotDts()
        init.SGA.plotDts() 
        init.CBBA.plotDts()
    if original:
        init.TGTA.plotDts() 
        init.DTTA.plotDts() 
        init.T3A.plotDts() 
        init.DSTA.plotDts() 
        init.STTA.plotDts() 
    if bundle:
        init.TBTA.plotDts() 
        init.TTBTA.plotDts()
        init.STBTA.plotDts()
    if lazy:
        init.LSGA.plotDts() 
        init.LTGTA.plotDts()
        init.LDTTA.plotDts()
        init.LT3A.plotDts()
        init.LSTA.plotDts()
        init.LSTTA.plotDts()
    # Plot figure
    plt.xlabel("$N_a$")
    plt.ylabel("Total Time (sec)")
    plt.legend(loc = 'best')  
#    sub_axes = plt.axes([.55, .25, .3, .3]) 
#    sub_axes.plot(gp.Nas, gp.dts_LSTA_Opt, 'g-^')
    plt.ylim(bottom=0)
    if original and not lazy:
        plt.savefig(output_path+'/original/dt_monotone_cmp_original.pdf')
    elif lazy and not original:
        plt.savefig(output_path+'/lazy/dt_monotone_cmp_lazy.pdf')
    else:
        plt.savefig(output_path+'/all/dt_monotone_cmp_all.pdf')
    

# total running time /sec   no CBBA    
    plt.figure("time no CBBA")
    if benchmark:
        init.GA.plotDts()
        init.AuctionSequential.plotDts()
        init.AuctionParallel.plotDts()
        init.AuctionGPrim.plotDts()
        init.SGA.plotDts() 
#        init.CBBA.plotDts() 
    if original:
        init.TGTA.plotDts() 
        init.DTTA.plotDts() 
        init.T3A.plotDts() 
        init.DSTA.plotDts() 
        init.STTA.plotDts() 
    if bundle:
        init.TBTA.plotDts()
        init.TTBTA.plotDts()
        init.STBTA.plotDts()
    if lazy:
        init.LSGA.plotDts() 
        init.LTGTA.plotDts()
        init.LDTTA.plotDts()
        init.LT3A.plotDts()
        init.LSTA.plotDts()
        init.LSTTA.plotDts()
    # Plot figure
    plt.xlabel("$N_a$")
    plt.ylabel("Total Time (sec)")
    plt.legend(loc = 'best')  
#    sub_axes = plt.axes([.55, .25, .3, .3]) 
#    sub_axes.plot(gp.Nas, gp.dts_LSTA_Opt, 'g-^')
    plt.ylim(bottom=0)
    if original and not lazy:
        plt.savefig(output_path+'/original/dt_monotone_cmp_noCBBA_original.pdf')
    elif lazy and not original:
        plt.savefig(output_path+'/lazy/dt_monotone_cmp_noCBBA_lazy.pdf')
    else:
        plt.savefig(output_path+'/all/dt_monotone_cmp_noCBBA_all.pdf')


# number of consensus steps       
    plt.figure("steps")
    if benchmark:
        init.GA.plotSteps()
        init.AuctionSequential.plotSteps()
        init.AuctionParallel.plotSteps()
        init.AuctionGPrim.plotSteps()
        init.SGA.plotSteps() 
        init.CBBA.plotSteps() 
    if original:
        init.TGTA.plotSteps() 
        init.DTTA.plotSteps() 
        init.T3A.plotSteps() 
        init.DSTA.plotSteps() 
        init.STTA.plotSteps() 
    if bundle:
        init.TBTA.plotSteps()
        init.TTBTA.plotSteps()
        init.STBTA.plotSteps()
    if lazy:
        init.LSGA.plotSteps() 
        init.LTGTA.plotSteps()
        init.LDTTA.plotSteps()
        init.LT3A.plotSteps()
        init.LSTA.plotSteps()
        init.LSTTA.plotSteps()
    # Plot figure
    plt.xlabel("$N_a$")
    plt.ylabel("Consensus Steps")
    plt.legend(loc = 'best')
    plt.ylim(bottom=0, top=init.Nt*1.1)
    if original and not lazy:
        plt.savefig(output_path+'/original/steps_monotone_cmp_original.pdf')
    elif lazy and not original:
        plt.savefig(output_path+'/lazy/steps_monotone_cmp_lazy.pdf')
    else:
        plt.savefig(output_path+'/all/steps_monotone_cmp_all.pdf')
    
    
    #number of function evaluations     
    plt.figure("evs")
    if benchmark:
        init.GA.plotEvs()
        init.AuctionSequential.plotEvs()
        init.AuctionParallel.plotEvs()
        init.AuctionGPrim.plotEvs()
        init.SGA.plotEvs() 
        init.CBBA.plotEvs() 
    if original:
        init.TGTA.plotEvs() 
        init.DTTA.plotEvs() 
        init.T3A.plotEvs() 
        init.DSTA.plotEvs() 
        init.STTA.plotEvs() 
    if bundle:
        init.TBTA.plotEvs()
        init.TTBTA.plotEvs()
        init.STBTA.plotEvs()
    if lazy:
        init.LSGA.plotEvs() 
        init.LTGTA.plotEvs()
        init.LDTTA.plotEvs()
        init.LT3A.plotEvs()
        init.LSTA.plotEvs()
        init.LSTTA.plotEvs()
    # Plot figure    
    plt.xlabel("$N_a$")
    plt.ylabel("Function Evaluations (x$10^"+str(init.evs_scale)+"$)")
    plt.legend(loc = 'best')
    plt.ylim(bottom=0)
    if original and not lazy:
        plt.savefig(output_path+'/original/evs_monotone_cmp_original.pdf')
    elif lazy and not original:
        plt.savefig(output_path+'/lazy/evs_monotone_cmp_lazy.pdf')
    else:
        plt.savefig(output_path+'/all/evs_monotone_cmp_all.pdf')
 
    
    #number of function evaluations without CBBA
    plt.figure("evs no CBBA")
    if benchmark:
        init.GA.plotEvs() 
        init.AuctionSequential.plotEvs()
        init.AuctionParallel.plotEvs()
        init.AuctionGPrim.plotEvs()
        init.SGA.plotEvs() 
#        init.CBBA.plotEvs()
    if original:
        init.TGTA.plotEvs() 
        init.DTTA.plotEvs() 
        init.T3A.plotEvs() 
        init.DSTA.plotEvs() 
        init.STTA.plotEvs() 
    if bundle:
        init.TBTA.plotEvs()
        init.TTBTA.plotEvs()
        init.STBTA.plotEvs()
    if lazy:
        init.LSGA.plotEvs() 
        init.LTGTA.plotEvs()
        init.LDTTA.plotEvs()
        init.LT3A.plotEvs()
        init.LSTA.plotEvs()
        init.LSTTA.plotEvs()
    # Plot figure 
    plt.xlabel("$N_a$")
    plt.ylabel("Function Evaluations (x$10^"+str(init.evs_scale)+"$)")
    plt.legend(loc = 'best')
    plt.ylim(bottom=0)
    if original and not lazy:
        plt.savefig(output_path+'/original/evs_monotone_cmp_noCBBA_original.pdf')
    elif lazy and not original:
        plt.savefig(output_path+'/lazy/evs_monotone_cmp_noCBBA_lazy.pdf')
    else:
        plt.savefig(output_path+'/all/evs_monotone_cmp_noCBBA_all.pdf')
    
    
# =============================================================================
# 
# =============================================================================
print("----- plotall.py is loaded -----")
