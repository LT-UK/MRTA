#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  8 14:35:52 2020

@author: Teng Li
lt.uk@outlook.com
United Kingdom
All Rights Reserved
"""

# This file is used for testing only.

p = np.arange(0.75,1.0,0.05)

X = [(1/p[i]**2-1)*st.mean(init.DSTA.utilities[i+14])**2-Var[i+14] for i in range(len(p))]

(1 - G_a^2)*max[E(f(s))]^2 - Var 

p = np.arange(0.05,0.71,0.05)

Y = [(1- (p[i]*(1-p[i])/(p[i]+max(p[i],1-p[i])))**2) * max(init.DSTA.utilities_avg)**2 -Var[i]
    for i in range(len(p))
]
