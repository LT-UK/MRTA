#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 18 15:58:04 2020

Value Function

@author: Teng Li
lt.uk@outlook.com
United Kingdom
All Rights Reserved
"""

import math
import init

# =============================================================================
# 
# =============================================================================
def getDistance(j_1, j_2, tasks_pos):
    '''
    Get the distance between task j_1 and j_2.
    Input:
        j_1, j_2: [int] task id
        tasks_pos: [list 2D] tasks' positions
    Output:
        dist: [float] distance between task j_1 and j_2
    '''
    dx = tasks_pos[0][j_1] - tasks_pos[0][j_2]
    dy = tasks_pos[1][j_1] - tasks_pos[1][j_2]
    dist = math.sqrt(dx**2+dy**2)    
    return dist

# =============================================================================
# 
# =============================================================================
def valueFunc(S, T):
    '''
    Get the functioin value, monotone
    
    Input:
        S: [list] current selection
        T: [list] all tasks, indices
    Output:
        f_value: [float] function value of the selection S
    '''

    n = len(list(S))
    if n < 2:
        return 0
    
    agent_id = S[0]
    remained = list(T)
    sum_1 = 0
    sum_2 = 0
    
    for j in S[1:]:
        sum_1 += init.M[agent_id][j]*init.V[j]
        remained.remove(j)
    
    for j in remained:
        d_min = min([init.Dist_mat[j][t] for t in S[1:]])
        sum_2 += init.M[agent_id][j]*init.V[j]*math.e**(-d_min/init.d_0)
    
    f_value = sum_1 + sum_2
    
    return f_value


def getDistAT(a, j, agents_pos, tasks_pos):
    '''
    Get the distance between agent a and task j.
    Input:
        a: [int] agent id
        j: [int] task id
        agents_pos: [list 2D] agents' positions
        tasks_pos: [list 2D] tasks' positions
    Output:
        dist: [float] distance between task j_1 and j_2
    '''
    dx = agents_pos[0][a] - tasks_pos[0][j]
    dy = agents_pos[1][a] - tasks_pos[1][j]
    dist = math.sqrt(dx**2+dy**2)    
    return dist


def valueFuncNon(S, T):
    '''
    Get the functioin value, non-monotone
    
    Input:
        S: [list] current selection
        T: [list] all tasks, indices
    Output:
        f_value: [float] function value of the selection S
    '''
    n = len(list(S))
    if n < 2:
        return 0
    
    agent_id = S[0]
    sum_1 = 0
    penalty = 0
    
    for j in S[1:]:
        sum_1 += init.M[agent_id][j]*init.V[j]
    
    if n > 2:
        for i in range(1,n-1):
            for j in range(i+1,n):
                penalty += init.X[S[i]][S[j]]

    f_value = sum_1 - init.lambda_x*penalty
    
    return f_value


# =============================================================================
# original getMGV, slow
# =============================================================================
def getMGV_Original(u, S, T):
    '''
    Get the marginal value of task u given S with valueFunc
    
    Input:
        u: [int] new task id
        S: [list] current selection
        T: [list] all tasks, indices

    Output:
        mgv: [float] marginal value for task u given S
    '''
    selected = list(S)
    new_set = list(S)
    new_set.append(u)
    
    if init.monotonicity:
        mgv = valueFunc(new_set, T) - valueFunc(selected, T)
    else:
        mgv = valueFuncNon(new_set, T) - valueFuncNon(selected, T)
    return mgv


# =============================================================================
# Optimised getMGV, fast
# =============================================================================
def getMGV(u, S, T):
    '''
    Get the marginal value of task u given S in both monotone and non-monotone
    cases. 
    
    The calculation of mgv is simplified.
    
    Input:
        u: [int] new task id
        S: [list] current selection
        T: [list] all tasks, indices

    Output:
        mgv: [float] marginal value for task u given S
    '''
    agent_id = S[0]
    n = len(list(S))
    
    # no task has been selected
    if n < 2:
        if init.monotonicity:
            mgv = valueFunc(S+[u], T)
        else:
            mgv = valueFuncNon(S+[u], T)
        return mgv
    
    if init.monotonicity:    
        remained = list(T)
        for j in S[1:]:
            remained.remove(j)
        remained.remove(u)
        
        d_min_u = min([init.Dist_mat[u][t] for t in S[1:]])
        term_1 = init.M[agent_id][u]*init.V[u]*(1 - math.e**(-d_min_u/init.d_0))
        
        term_2 = 0
        for j in remained:
            d_min_j = min([init.Dist_mat[j][t] for t in S[1:]])
            d_j_u = init.Dist_mat[j][u]
            if d_j_u < d_min_j:
                term_2 += init.M[agent_id][j]*init.V[j] \
                        *(math.e**(-d_j_u/init.d_0) - math.e**(-d_min_j/init.d_0))
        mgv = term_1 + term_2
    else:
        delta_profit = init.M[agent_id][u]*init.V[u]
        # penalty term in the non-monotone case
        delta_penalty = 0
        if not init.monotonicity:
            for i in range(n-1):
                delta_penalty += init.X[S[i+1]][u]
        # final mgv
        mgv = delta_profit - init.lambda_x*delta_penalty
    
    return mgv

# =============================================================================
# 
# =============================================================================
print("----- vfunc.py is loaded -----")
