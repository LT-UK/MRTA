#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 21 18:54:21 2021

@author: Teng Li
tengli@cranfield.ac.uk; nliteng@foxmail.com
Cranfield University, UK
All Rights Reserved
"""

import numpy as np
#import time
import os
import math
import statistics as st
import matplotlib.pyplot as plt

import init 
import funcs.vfunc.vf as vf          # load value function and getMGV

import funcs.algs.ga as ga         # TA alg: GA
import funcs.algs.tgta as tgta     # TA alg: TGTA
import funcs.algs.dtta as dtta     # TA alg: DTTA
import funcs.algs.t3a as t3a     # TA alg: T3A
import funcs.algs.dsta as dsta     # TA alg: DSTA
import funcs.algs.stta as stta     # TA alg: DSTA

SMALL_SIZE = 11
MEDIUM_SIZE = 13
BIGGER_SIZE = 14

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=BIGGER_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=MEDIUM_SIZE)  # fontsize of the figure title


output_path = os.path.abspath("output") + "/tune_ga/"

# ====================== Set TA parameters between here =======================

init.Nt = 50       # number of tasks 200, 60, 50   200 tasks too slow
init.evs_scale = 3  # scale in plot figures, scale = 10**evs_scale

init.Na = 12
init.Nst = 15

init.L = 10.0       # /km  The length of the 2-D space.

# value factors of the normal tasks
init.v_low = 0.6    # 0.6
init.v_high = 1.0   # 1.0
# value factors of special tasks in the non-monotone case
init.vs_low = 5   # 5
init.vs_high = 6  # 6

# match fitness factors for normal tasks
init.m_low = 0.5    # 0.5
init.m_high = 1.0   # 1.0
# match fitness factors for special tasks in the non-monotone case
init.ms_low = 0.1   # 0.1
init.ms_high = 0.2  # 0.2

init.lambda_x = 0.01     # penalty scaling factor in the non-monotone case

init.d_0 = 1        # /km characteristic distance

init.eps = 0.1      # threshold decreasing parameter \epsilon
init.rho = 1        # truncation parameter

init.Pr = 0.5        # sampling probability

init.MC_ctr = 10     # The number of Monte Carlo loops.
    
# GA parameters
init.max_ite = 8000 # The maximum number of iterations.
init.pop_size = 20  # 20 Population size.
init.Pr_mu = 0.5   # 0.5 The probability of mutation.

init.TGTA.en = 0

init.DTTA.en = 1

init.T3A.en = 0

init.DSTA.en = 0

init.STTA.en = 0

Random_en = 0

init.monotonicity = 1    # 1: monotone 0: non-monotone


def getSuppass(val_list, value):
    '''
    Return the (index+1) of the 1st element whose value >= "value"
    val_list values are increasing
    '''
    for i in range(len(val_list)):
        if val_list[i] >= value:
            return i+1
    return False

def getConvergence(val_list):
    '''
    Return the (index+1) of the 1st element reaching convergence
    val_list values are increasing
    '''
    for i in range(len(val_list)):
        if val_list[i] >= val_list[-1]:
            return i+1
#    return False
    


def plotTuneGA(Alg, Alg_utility, Alg_evs, max_utilities, GA_evs):

    # Plot utility   
    plt.figure("GA utilities wrt iterations")
    plt.plot(max_utilities, 'r', label='GA')
    plt.plot(Alg_utility, 'g', label=Alg)
    if Random_en:
        plt.plot([total_value_avg for value in range(len(Alg_utility))], 'b', label="Random")
    plt.xlabel("$Iterations$")
    plt.ylabel("Function Utility")
    plt.legend(loc = 'lower right')
#    plt.ylim(bottom=max_utilities[0]-10, top=max_utilities[-1]*1.05)
    plt.ylim(bottom=0, top=max_utilities[-1]*1.05)
#    plt.ylim(bottom=-100, top=50)
#    plt.ylim(bottom=0)
    plt.savefig(output_path+'GA_tune_utility.pdf')
    
    # Plot evs
    plt.figure("GA evs wrt iterations")
    plt.plot(GA_evs, 'r', label='GA')
    plt.plot(Alg_evs, 'g', label=Alg)
    plt.xlabel("$Iterations$")
    plt.ylabel("Function Evaluations (x$10^"+str(init.evs_scale)+"$)")
    plt.legend(loc = 'lower right')
    plt.ylim(bottom=0)
    plt.savefig(output_path+'GA_tune_evs.pdf')
    
    
def plotTuneGA_Integrate_Variance(Alg, Alg_utility, Alg_evs, max_utilities, GA_evs):
    
    # error bar settings
    eb_fmt = "c."
    eb_markersize = 3
    eb_capsize = 2
    eb_capsize_line = 4
    eb_capsize_bar = 2
    eb_label = Alg
    
    # Plot utility   
    plt.figure("GA and Alg utilities vs evs")
    plt.plot(GA_evs, max_utilities, 'r', label='GA')
#    plt.plot(Alg_evs, Alg_utility, 'g-o', label=Alg)
    
    Alg_utility_avg = [st.mean(Alg_utility)]
    Alg_evs_avg = [st.mean(Alg_evs)]
    y_upper_error = [max(Alg_utility)-Alg_utility_avg[0]]
    y_lower_error = [-min(Alg_utility)+Alg_utility_avg[0]]
    yerr_utilities = [y_lower_error, y_upper_error]
    
    x_upper_error = [max(Alg_evs)-Alg_evs_avg[0]]
    x_lower_error = [-min(Alg_evs)+Alg_evs_avg[0]]
    xerr_evs = [x_lower_error, x_upper_error]
    
    plt.errorbar(Alg_evs_avg,
                 Alg_utility_avg,
                 yerr=yerr_utilities,
                 xerr=xerr_evs,
                 fmt=eb_fmt,
                 capsize=eb_capsize,
                 label = eb_label)
    
#    plt.errorbar(Alg_evs_avg,
#                 Alg_utility_avg,
##                 yerr=yerr_utilities,
#                 xerr=xerr_evs,
#                 fmt=eb_fmt,
#                 capsize=eb_capsize,
#                 label = eb_label)

    plt.xlabel("Function Evaluations (x$10^"+str(init.evs_scale)+"$)")
    plt.ylabel("Function Utility")
    plt.legend(loc = 'best')
#    plt.ylim(bottom=max_utilities[0]-10, top=max_utilities[-1]*1.05)
#    plt.ylim(bottom=0, top=max_utilities[-1]*1.05)
    plt.ylim(bottom=0, top=150)
#    plt.ylim(bottom=0)
    plt.savefig(output_path+'GA_tune_utility_vs_evs_vr.pdf')
    
    
def plotTuneGA_Integrate(Alg, Alg_utility, Alg_evs, max_utilities, GA_evs):
    
    # Plot utility   
    plt.figure("GA and Alg utilities vs evs")
    plt.plot(GA_evs, max_utilities, 'r', label='GA')
#    plt.plot(Alg_evs, Alg_utility, 'g-o', label=Alg)

    plt.plot(Alg_evs, Alg_utility, 'gx', label=Alg)

    plt.xlabel("Function Evaluations (x$10^"+str(init.evs_scale)+"$)")
    plt.ylabel("Function Utility")
    plt.legend(loc = 'lower right')
#    plt.ylim(bottom=max_utilities[0]-10, top=max_utilities[-1]*1.05)
#    plt.ylim(bottom=0, top=max_utilities[-1]*1.05)
    plt.ylim(bottom=70, top=100)
#    plt.ylim(bottom=0)
    plt.savefig(output_path+'GA_tune_utility_vs_evs.pdf')

# ================================= and here ==================================
if __name__ == "__main__":

    Tasks = [j for j in range(init.Nt)]    
    init.Dist_mat = [[0 for j in range(init.Nt)] for i in range(init.Nt)]
    
    Agents = [a for a in range(init.Na)]
    
    Tasks_pos = [
            [round(x,3) for x in 
             np.random.uniform(low=0.0, high=init.L, size=init.Nt)], 
            [round(y,3) for y in 
             np.random.uniform(low=0.0, high=init.L, size=init.Nt)]]
    
    # get the distance matrix of all tasks
    for i in range(init.Nt-1):
        for j in range(i+1,init.Nt):
            init.Dist_mat[i][j] = round(vf.getDistance(i,j,Tasks_pos),3)
            init.Dist_mat[j][i] = init.Dist_mat[i][j]
    
    # monotone case
    if init.monotonicity:
        # random importance factors
        init.V = [round(v,3) for v in 
                  np.random.uniform(low=init.v_low, 
                                    high=init.v_high, 
                                    size=init.Nt)
                  ]
        # random task-agent match fitness factors
        #    M = np.random.uniform(m_low,m_high,(Na,Nt))
        init.M = [[round(m,3) for m in 
                   np.random.uniform(init.m_low,init.m_high,init.Nt)] 
                    for a in range(init.Na)
                 ]
#    Non-monotone:
    else:
        # random importance factors
#        init.V = [round(v,3) for v in 
#                  np.random.uniform(low=init.v_low, 
#                                    high=init.v_high, 
#                                    size=init.Nt)
#                  ]
        # random task-agent match fitness factors
        #    M = np.random.uniform(m_low,m_high,(Na,Nt))
#        init.M = [[round(m,3) for m in 
#                   np.random.uniform(init.m_low,init.m_high,init.Nt)] 
#                    for a in range(init.Na)
#                 ]
#        # special tasks
        v1 = [round(v,3) for v in 
              np.random.uniform(low=init.vs_low, 
                                high=init.vs_high, 
                                size=init.Nst)
             ]
#        # normal tasks
        v2 = [round(v,3) for v in 
              np.random.uniform(low=init.v_low, 
                                high=init.v_high, 
                                size=init.Nt-init.Nst)
             ]
        init.V = v1 + v2
#        # special tasks
        m1 = init.ms_low*np.ones(init.Na) + \
            (init.ms_high-init.ms_low)*np.eye(init.Na)
#        # normal tasks
        m2 = [[round(m,3) for m in 
               np.random.uniform(init.m_low,init.m_high,init.Nt-init.Nst)] 
                for a in range(init.Na)
             ]
        init.M = np.concatenate((m1, m2), axis=1) 
        # inter-task dependences Nt*Nt
        init.X = np.zeros((init.Nt, init.Nt))
        for i in range(init.Nt-1):
            for j in range(i+1, init.Nt):
                init.X[i][j] = round(math.e**(init.V[i]*init.V[j]),4)
                init.X[j][i] = init.X[i][j]
            
                    
    selected, values, total_value, dt, steps, evs, max_utilities, evs_ite \
        = ga.runGA_Tune(Agents,Tasks,init.pop_size, init.max_ite, init.Pr_mu)
    GA_evs = [ev/10**init.evs_scale for ev in evs_ite]
    
#    Find GA suppass point, convergence point, (n-th iteration) 
#    Corresponding evs
    
    if Random_en:
        total_values, total_value_avg = ga.runRandom(Agents, Tasks, init.MC_ctr)
    
    
    if init.DTTA.en:
        selected, values, total_value, dt, steps, evs \
                        = dtta.runDTTA(Agents,Tasks,init.eps)
        Alg_utility = [total_value for i in range(len(max_utilities))]
        Alg_evs = [evs/10**init.evs_scale for i in range(len(max_utilities))]
        # plotTuneGA('DTTA', Alg_utility, Alg_evs, max_utilities, GA_evs)
        plotTuneGA_Integrate('DTTA', total_value, evs/10**init.evs_scale, max_utilities, GA_evs)
        utility_avg = total_value
        

    if init.TGTA.en:
        selected, values, total_value, dt, steps, evs \
                        = tgta.runTGTA(Agents,Tasks,init.rho)
        Alg_utility = [total_value for i in range(len(max_utilities))]
        Alg_evs = [evs/10**init.evs_scale for i in range(len(max_utilities))]
        
        # plotTuneGA('TGTA', Alg_utility, Alg_evs, max_utilities, GA_evs)
        plotTuneGA_Integrate('TGTA', total_value, evs/10**init.evs_scale, max_utilities, GA_evs)
        utility_avg = total_value
        
    
    if init.T3A.en:
        selected, values, total_value, dt, steps, evs \
                        = t3a.runT3A(Agents,Tasks,init.rho,init.eps)
        Alg_utility = [total_value for i in range(len(max_utilities))]
        Alg_evs = [evs/10**init.evs_scale for i in range(len(max_utilities))]
#        plotTuneGA('T3A', Alg_utility, Alg_evs, max_utilities, GA_evs)
        plotTuneGA_Integrate('T3A', total_value, evs/10**init.evs_scale, max_utilities, GA_evs)
        
    # Monte-Carlo Loop
    for counter in range(init.MC_ctr):
        if init.DSTA.en:    
            selected, values, total_value, dt, steps, evs \
                = dsta.runDSTA(Agents,Tasks,init.Pr)
            init.DSTA.utilities_temp.append(total_value)
            init.DSTA.dts_temp.append(dt)
            init.DSTA.steps_temp.append(steps)
            init.DSTA.evs_temp.append(evs) 
            
        if init.STTA.en:    
            selected, values, total_value, dt, steps, evs \
                = stta.runSTTA(Agents,Tasks,init.Pr,init.eps)
            init.STTA.utilities_temp.append(total_value)
            init.STTA.dts_temp.append(dt)
            init.STTA.steps_temp.append(steps)
            init.STTA.evs_temp.append(evs) 
    
    if init.DSTA.en:
        utility_max = max(init.DSTA.utilities_temp)
        Alg_utility_max = [utility_max for i in range(len(max_utilities))]
        index_max = init.DSTA.utilities_temp.index(utility_max)
        evs_max = init.DSTA.evs_temp[index_max]
        Alg_evs_max = [evs_max/10**init.evs_scale for i in range(len(max_utilities))]
        steps_max = init.DSTA.steps_temp[index_max]
        Alg_steps_max = [steps_max for i in range(len(max_utilities))]
        
        utility_avg = st.mean(init.DSTA.utilities_temp)
        Alg_utility_avg = [utility_avg for i in range(len(max_utilities))]
        Alg_steps_avg = st.mean(init.DSTA.steps_temp)
        evs_avg = st.mean(init.DSTA.evs_temp)
        Alg_evs_avg = [evs_avg/10**init.evs_scale for i in range(len(max_utilities))]
        
        Alg_utilities = init.DSTA.utilities_temp
        Alg_evs = [evs/10**init.evs_scale for evs in init.DSTA.evs_temp]
        
#        plotTuneGA('DSTA', Alg_utility_max, Alg_evs_max, max_utilities, GA_evs)
        plotTuneGA('DSTA', Alg_utility_avg, Alg_evs_avg, max_utilities, GA_evs)
        plotTuneGA_Integrate_Variance('DSTA', Alg_utilities, Alg_evs, max_utilities, GA_evs)
    
    if init.STTA.en:
        
        Alg_utilities = init.STTA.utilities_temp
        Alg_evs = [evs/10**init.evs_scale for evs in init.STTA.evs_temp]
        
#        plotTuneGA('DSTA', Alg_utility_max, Alg_evs_max, max_utilities, GA_evs)
#        plotTuneGA('STTA', Alg_utility_avg, Alg_evs_avg, max_utilities, GA_evs)
        plotTuneGA_Integrate_Variance('STTA', Alg_utilities, Alg_evs, max_utilities, GA_evs)
    
#    if init.DSTA.en:
#        selected, values, total_value, dt, steps, evs \
#                        = dsta.runDSTA(Agents,Tasks,init.Pr)
#        Alg_utility = [total_value for i in range(len(max_utilities))]
#        Alg_evs = [evs/10**init.evs_scale for i in range(len(max_utilities))]
#        plotTuneGA('DSTA', Alg_utility, Alg_evs, max_utilities, GA_evs)


#    Print utility% and milestones
    if init.DSTA.en or init.STTA.en:
        final_ratio = utility_avg/max_utilities[-1]
        
        suppass = getSuppass(max_utilities, utility_avg)
        convergence = getConvergence(max_utilities)
        GA_evs[suppass]
        GA_evs[convergence]
        Alg_evs_avg[0]
        max_utilities[-1]
    #    Alg_utility_avg[0]
        Alg_utility_avg = [st.mean(Alg_utilities)]
        Alg_evs_avg = [st.mean(Alg_evs)]
        
        suppass = getSuppass(max_utilities, total_value)
    else:
        final_ratio = total_value/max_utilities[-1]
        suppass = getSuppass(max_utilities, total_value)
        convergence = getConvergence(max_utilities)
        suppass
        GA_evs[suppass]
        GA_evs[convergence]
        total_value
        max_utilities[-1]
        final_ratio
        Alg_evs[0]
        Alg_evs[0]/GA_evs[suppass]
