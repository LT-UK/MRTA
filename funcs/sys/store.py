#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 25 14:18:50 2020

@author: Teng Li
tengli@cranfield.ac.uk; nliteng@foxmail.com
Cranfield University, UK
All Rights Reserved
"""

#import os
import init
import pickle # store and restore variables

# =============================================================================
# 
# =============================================================================

def storeSettings(path):
    '''
    Store settings to a pckl file.
    '''    
#    store_path = os.path.abspath("output") + '/data/settings.pckl'
    store_path = path + '/data/settings.pckl'
    
    settingdata = [
        init.Nt, init.evs_scale, init.Na,
        init.L, 
        init.v_low, init.v_high, 
        init.m_low, init.m_high,
        init.vs_low, init.vs_high,
        init.ms_low, init.ms_high,
        init.monotonicity, 
        init.lambda_x, init.d_0,
        init.MC_ctr, 
        init.Na_start, init.Na_end, init.Na_step,
        init.Pr_start, init.Pr_end, init.Pr_step,
        init.Eps_start, init.Eps_end, init.Eps_step,
        init.Rho_start, init.Rho_end, init.Rho_step,
        init.Pr, init.eps, init.rho,
        init.SGA.en, init.LSGA.en, 
        init.CBBA.en,
        init.TGTA.en, init.LTGTA.en,
        init.DTTA.en, init.LDTTA.en, init.TBTA.en,
        init.T3A.en, init.LT3A.en, init.TTBTA.en,
        init.DSTA.en, init.LSTA.en,
        init.STTA.en, init.LSTTA.en, init.STBTA.en]
    
    with open(store_path, 'wb') as f:
        pickle.dump(settingdata, f)
    # file f has been automatically closed


#settings_path = os.path.abspath("output")+'/data/settings.pckl'

def restoreSettings(path):
    '''
    Restore settings from a pckl file.
    '''    
    with open(path, 'rb') as f:
        [init.Nt, init.evs_scale, init.Na,
         init.L, 
         init.v_low, init.v_high, 
         init.m_low, init.m_high,
         init.vs_low, init.vs_high,
         init.ms_low, init.ms_high,
         init.monotonicity, 
         init.lambda_x, init.d_0,
         init.MC_ctr, 
         init.Na_start, init.Na_end, init.Na_step,
         init.Pr_start, init.Pr_end, init.Pr_step,
         init.Eps_start, init.Eps_end, init.Eps_step,
         init.Rho_start, init.Rho_end, init.Rho_step,         
         init.Pr, init.eps, init.rho,
         init.SGA.en, init.LSGA.en, 
         init.CBBA.en,
         init.TGTA.en, init.LTGTA.en,
         init.DTTA.en, init.LDTTA.en, init.TBTA.en,
         init.T3A.en, init.LT3A.en,  init.TTBTA.en,
         init.DSTA.en, init.LSTA.en,
         init.STTA.en, init.LSTTA.en, init.STBTA.en
        ] = pickle.load(f)
    # file f has been automatically closed


# =============================================================================
# 
# =============================================================================

def storePerformance(path):
    '''
    Store algorithm to a pckl file.
    '''    
#    store_path = os.path.abspath("output") + '/data/performance.pckl'
    store_path = path + '/data/performance.pckl'
    
    performance = [
        init.Nas, init.Prs, init.Epses, init.Rhos,
        init.SGA.utilities,   init.SGA.dts,   init.SGA.steps,   init.SGA.evs,
        init.LSGA.utilities,  init.LSGA.dts,  init.LSGA.steps,  init.LSGA.evs,
        init.CBBA.utilities,  init.CBBA.dts,  init.CBBA.steps,  init.CBBA.evs,
        init.TGTA.utilities,  init.TGTA.dts,  init.TGTA.steps,  init.TGTA.evs,
        init.LTGTA.utilities, init.LTGTA.dts, init.LTGTA.steps, init.LTGTA.evs,
        init.DTTA.utilities,  init.DTTA.dts,  init.DTTA.steps,  init.DTTA.evs,
        init.LDTTA.utilities, init.LDTTA.dts, init.LDTTA.steps, init.LDTTA.evs,
        init.TBTA.utilities,  init.TBTA.dts,  init.TBTA.steps,  init.TBTA.evs,
        init.T3A.utilities,   init.T3A.dts,   init.T3A.steps,   init.T3A.evs,
        init.LT3A.utilities,  init.LT3A.dts,  init.LT3A.steps,  init.LT3A.evs,
        init.TTBTA.utilities, init.TTBTA.dts, init.TTBTA.steps, init.TTBTA.evs,
        init.DSTA.utilities,  init.DSTA.dts,  init.DSTA.steps,  init.DSTA.evs,
        init.LSTA.utilities,  init.LSTA.dts,  init.LSTA.steps,  init.LSTA.evs,
        init.STTA.utilities,  init.STTA.dts,  init.STTA.steps,  init.STTA.evs,
        init.LSTTA.utilities, init.LSTTA.dts, init.LSTTA.steps, init.LSTTA.evs,
        init.STBTA.utilities, init.STBTA.dts, init.STBTA.steps, init.STBTA.evs
        ]
    
    with open(store_path, 'wb') as f:
        pickle.dump(performance, f)
    # file f has been automatically closed


#performance_path = os.path.abspath("output")+'/data/performance.pckl'

def restorePerformance(path):
    '''
    Restore algorithms' performance from a pckl file.
    ''' 
    with open(path, 'rb') as f:
        [init.Nas, init.Prs, init.Epses, init.Rhos,
         init.SGA.utilities,   init.SGA.dts,   init.SGA.steps,   init.SGA.evs,
         init.LSGA.utilities,  init.LSGA.dts,  init.LSGA.steps,  init.LSGA.evs,
         init.CBBA.utilities,  init.CBBA.dts,  init.CBBA.steps,  init.CBBA.evs,
         init.TGTA.utilities,  init.TGTA.dts,  init.TGTA.steps,  init.TGTA.evs,
         init.LTGTA.utilities, init.LTGTA.dts, init.LTGTA.steps, init.LTGTA.evs,
         init.DTTA.utilities,  init.DTTA.dts,  init.DTTA.steps,  init.DTTA.evs,
         init.LDTTA.utilities, init.LDTTA.dts, init.LDTTA.steps, init.LDTTA.evs,
         init.TBTA.utilities,  init.TBTA.dts,  init.TBTA.steps,  init.TBTA.evs,
         init.T3A.utilities,   init.T3A.dts,   init.T3A.steps,   init.T3A.evs,
         init.LT3A.utilities,  init.LT3A.dts,  init.LT3A.steps,  init.LT3A.evs,
         init.TTBTA.utilities, init.TTBTA.dts, init.TTBTA.steps, init.TTBTA.evs,
         init.DSTA.utilities,  init.DSTA.dts,  init.DSTA.steps,  init.DSTA.evs,
         init.LSTA.utilities,  init.LSTA.dts,  init.LSTA.steps,  init.LSTA.evs,
         init.STTA.utilities,  init.STTA.dts,  init.STTA.steps,  init.STTA.evs,
         init.LSTTA.utilities, init.LSTTA.dts, init.LSTTA.steps, init.LSTTA.evs,
         init.STBTA.utilities, init.STBTA.dts, init.STBTA.steps, init.STBTA.evs
        ] = pickle.load(f)
    # file f has been automatically closed

# =============================================================================
# 
# =============================================================================
def storeData(path):
    storeSettings(path)
    storePerformance(path)
    

def restoreData(path):
    # settings
    settings_path = path+'/settings.pckl'
    restoreSettings(settings_path)
    # performance
    performance_path = path+'/performance.pckl'
    restorePerformance(performance_path)

