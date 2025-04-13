#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 25 12:12:31 2020

All "import" lines need executing before running this file.
The Spyder execution direct should be at the current address same with "data".

Do NOT run whole file. 

@author: Teng Li
lt.uk@outlook.com
United Kingdom
All Rights Reserved
"""

# =============================================================================
# Do NOT run these lines, they were just for debuging.
import funcs.sys.store as store
#import funcs.plots.plotall as plotall
#import funcs.plots.plotdsta as plotdsta
# =============================================================================


import os

settings_path = os.path.abspath("data")+'/settings.pckl'
performance_path = os.path.abspath("data")+'/performance.pckl'

store.restoreSettings(settings_path)
store.restorePerformance(performance_path)


