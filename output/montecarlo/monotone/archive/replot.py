#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 25 12:12:31 2020

@author: Teng Li
lt.uk@outlook.com
United Kingdom
All Rights Reserved
"""

import os
import init
import funcs.plots.plotall as plotall

settings_path = os.path.abspath("output")+'/data/settings.pckl'
performance_path = os.path.abspath("output")+'/data/performance.pckl'

init.restoreSettings(settings_path)
init.restorePerformance(performance_path)

plotall.plotAll()