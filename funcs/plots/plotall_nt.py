#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 18 15:58:04 2020

Plot figures for all algorithms.

@author: Teng Li
teng.li@cranfield.ac.uk; lt.uk@outlook.com
Cranfield University, UK
All Rights Reserved
"""

import matplotlib.pyplot as plt
import os # find path
import init

# =============================================================================
# 
# =============================================================================
if init.monotonicity:
    output_path = os.path.abspath("output") + "/montecarlo_nt/monotone/plotall_nt/"
else:
    output_path = os.path.abspath("output") + "/montecarlo_nt/non_monotone/plotall_nt/"

def plotAllNt(benchmark=1, original=1, bundle=1, lazy=1):
    # Draw figures
    print("\n=========== Plot All ===========")
    
# total values             
    plt.figure("utilities")
    if benchmark:
        init.GA.plotUtilities_Nt()
        init.AuctionSequential.plotUtilities_Nt()
        init.AuctionParallel.plotUtilities_Nt()
        init.AuctionGPrim.plotUtilities_Nt()
        init.SGA.plotUtilities_Nt() 
        init.CBBA.plotUtilities_Nt() 
    if original:
        init.TGTA.plotUtilities_Nt() 
        init.DTTA.plotUtilities_Nt() 
        init.T3A.plotUtilities_Nt() 
        init.DSTA.plotUtilities_Nt() 
        init.STTA.plotUtilities_Nt() 
    if bundle:
        init.TBTA.plotUtilities_Nt() 
        init.TTBTA.plotUtilities_Nt()
        init.STBTA.plotUtilities_Nt()
    if lazy:
        init.LSGA.plotUtilities_Nt() 
        init.LTGTA.plotUtilities_Nt()
        init.LDTTA.plotUtilities_Nt()
        init.LT3A.plotUtilities_Nt()
        init.LSTA.plotUtilities_Nt()
        init.LSTTA.plotUtilities_Nt()
    # Plot figure
    plt.xlabel("$N_t$")
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
        init.GA.plotDts_Nt()
        init.AuctionSequential.plotDts_Nt()
        init.AuctionParallel.plotDts_Nt()
        init.AuctionGPrim.plotDts_Nt()
        init.SGA.plotDts_Nt() 
        init.CBBA.plotDts_Nt()
    if original:
        init.TGTA.plotDts_Nt() 
        init.DTTA.plotDts_Nt() 
        init.T3A.plotDts_Nt() 
        init.DSTA.plotDts_Nt() 
        init.STTA.plotDts_Nt() 
    if bundle:
        init.TBTA.plotDts_Nt() 
        init.TTBTA.plotDts_Nt()
        init.STBTA.plotDts_Nt()
    if lazy:
        init.LSGA.plotDts_Nt() 
        init.LTGTA.plotDts_Nt()
        init.LDTTA.plotDts_Nt()
        init.LT3A.plotDts_Nt()
        init.LSTA.plotDts_Nt()
        init.LSTTA.plotDts_Nt()
    # Plot figure
    plt.xlabel("$N_t$")
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
        init.GA.plotDts_Nt()
        init.AuctionSequential.plotDts_Nt()
        init.AuctionParallel.plotDts_Nt()
        init.AuctionGPrim.plotDts_Nt()
        init.SGA.plotDts_Nt() 
#        init.CBBA.plotDts_Nt() 
    if original:
        init.TGTA.plotDts_Nt() 
        init.DTTA.plotDts_Nt() 
        init.T3A.plotDts_Nt() 
        init.DSTA.plotDts_Nt() 
        init.STTA.plotDts_Nt() 
    if bundle:
        init.TBTA.plotDts_Nt()
        init.TTBTA.plotDts_Nt()
        init.STBTA.plotDts_Nt()
    if lazy:
        init.LSGA.plotDts_Nt() 
        init.LTGTA.plotDts_Nt()
        init.LDTTA.plotDts_Nt()
        init.LT3A.plotDts_Nt()
        init.LSTA.plotDts_Nt()
        init.LSTTA.plotDts_Nt()
    # Plot figure
    plt.xlabel("$N_t$")
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
        init.GA.plotSteps_Nt()
        init.AuctionSequential.plotSteps_Nt()
        init.AuctionParallel.plotSteps_Nt()
        init.AuctionGPrim.plotSteps_Nt()
        init.SGA.plotSteps_Nt() 
        init.CBBA.plotSteps_Nt() 
    if original:
        init.TGTA.plotSteps_Nt() 
        init.DTTA.plotSteps_Nt() 
        init.T3A.plotSteps_Nt() 
        init.DSTA.plotSteps_Nt() 
        init.STTA.plotSteps_Nt() 
    if bundle:
        init.TBTA.plotSteps_Nt()
        init.TTBTA.plotSteps_Nt()
        init.STBTA.plotSteps_Nt()
    if lazy:
        init.LSGA.plotSteps_Nt() 
        init.LTGTA.plotSteps_Nt()
        init.LDTTA.plotSteps_Nt()
        init.LT3A.plotSteps_Nt()
        init.LSTA.plotSteps_Nt()
        init.LSTTA.plotSteps_Nt()
    # Plot figure
    plt.xlabel("$N_t$")
    plt.ylabel("Consensus Steps")
    plt.legend(loc = 'best')
    plt.ylim(bottom=0, top=init.Nt_end*1.1)
    if original and not lazy:
        plt.savefig(output_path+'/original/steps_monotone_cmp_original.pdf')
    elif lazy and not original:
        plt.savefig(output_path+'/lazy/steps_monotone_cmp_lazy.pdf')
    else:
        plt.savefig(output_path+'/all/steps_monotone_cmp_all.pdf')
    
    
    #number of function evaluations     
    plt.figure("evs")
    if benchmark:
        init.GA.plotEvs_Nt()
        init.AuctionSequential.plotEvs_Nt()
        init.AuctionParallel.plotEvs_Nt()
        init.AuctionGPrim.plotEvs_Nt()
        init.SGA.plotEvs_Nt() 
        init.CBBA.plotEvs_Nt() 
    if original:
        init.TGTA.plotEvs_Nt() 
        init.DTTA.plotEvs_Nt() 
        init.T3A.plotEvs_Nt() 
        init.DSTA.plotEvs_Nt() 
        init.STTA.plotEvs_Nt() 
    if bundle:
        init.TBTA.plotEvs_Nt()
        init.TTBTA.plotEvs_Nt()
        init.STBTA.plotEvs_Nt()
    if lazy:
        init.LSGA.plotEvs_Nt() 
        init.LTGTA.plotEvs_Nt()
        init.LDTTA.plotEvs_Nt()
        init.LT3A.plotEvs_Nt()
        init.LSTA.plotEvs_Nt()
        init.LSTTA.plotEvs_Nt()
    # Plot figure    
    plt.xlabel("$N_t$")
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
        init.GA.plotEvs_Nt() 
        init.AuctionSequential.plotEvs_Nt()
        init.AuctionParallel.plotEvs_Nt()
        init.AuctionGPrim.plotEvs_Nt()
        init.SGA.plotEvs_Nt() 
#        init.CBBA.plotEvs_Nt()
    if original:
        init.TGTA.plotEvs_Nt() 
        init.DTTA.plotEvs_Nt() 
        init.T3A.plotEvs_Nt() 
        init.DSTA.plotEvs_Nt() 
        init.STTA.plotEvs_Nt() 
    if bundle:
        init.TBTA.plotEvs_Nt()
        init.TTBTA.plotEvs_Nt()
        init.STBTA.plotEvs_Nt()
    if lazy:
        init.LSGA.plotEvs_Nt() 
        init.LTGTA.plotEvs_Nt()
        init.LDTTA.plotEvs_Nt()
        init.LT3A.plotEvs_Nt()
        init.LSTA.plotEvs_Nt()
        init.LSTTA.plotEvs_Nt()
    # Plot figure 
    plt.xlabel("$N_t$")
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
print("----- plotall_nt.py is loaded -----")
