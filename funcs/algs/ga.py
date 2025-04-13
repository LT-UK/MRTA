#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 20 21:24:52 2021

Genetic Algorithm for Task Allocation.

@author: Teng Li
tengli@cranfield.ac.uk; nliteng@foxmail.com
Cranfield University, UK
All Rights Reserved
"""

import numpy as np
import time
import funcs.vfunc.vf as vf   # value function and getMGV
import init
import statistics as st

def crossover_1point(P1, P2):
    '''
    Input:
        P1: [list] parent1 (one allocation solution) length=|Tasks|
        P2: [list] parent2
    Output:
        O1: [list] offspring1
        O2: [list] offspring2
    '''
    len_chrom = len(P1)
    cut_pos = np.random.randint(0, len_chrom-1, 1)[0]
    O1 = P1[0:cut_pos+1] + P2[cut_pos+1:len_chrom]
    O2 = P2[0:cut_pos+1] + P1[cut_pos+1:len_chrom]
    
    return O1, O2

    
def mutation_swap(individual):
    '''
    Mutation by randomly swap two points.
    '''
#    ind = [i for i in individual]
    ind = list(individual)
    len_chrom = len(ind)
    pos1, pos2 = np.random.randint(0, len_chrom-1, 2)
    ind[pos1], ind[pos2] = ind[pos2], ind[pos1]
    return ind


def getUtility(individual, Agents, Tasks):
    '''
    Input:
        individual: [list] One allocation solution.
        Agents: [list] The list of agent ids.
    Output:
        
    '''
    n_evs = 0
    len_chrom = len(individual) # equal to len(Tasks)
    # Create the empty allocation for each agent start from empty set. The first number is agent id.
    selected_ind = [[a] for a in Agents]
    values_ind = [0 for a in Agents] # Record the final function values for each agent.
    # Get the allocation result.
    for i in range(len_chrom):
        if individual[i] < len(Agents):
            selected_ind[individual[i]].append(i) 
    # Get the function utility for each agent.
    for a in Agents:
        if init.monotonicity: 
            values_ind[a] = vf.valueFunc(selected_ind[a], Tasks)
        else:
            values_ind[a] = vf.valueFuncNon(selected_ind[a], Tasks)
        n_evs += 1
    # Get total value for this individual.
    total_value = np.sum(values_ind)
    
    return values_ind, total_value, n_evs
    

def runRandom(Agents, Tasks, ites):
    '''
    Rully random selection
    '''
    len_chrom = len(Tasks)
    total_values = []
#    values = [0 for a in Agents] # Record the final function values for each agent.    
    # Create the empty allocation for each agent start from empty set. The first number is agent id.
#    selected = [[a] for a in Agents] 
    
    individual = list(np.random.randint(0, len(Agents), len_chrom))
    
    for ite in range(ites):
        values_ind, total_value, n_evs = getUtility(individual, Agents, Tasks)
        total_values.append(total_value)
    total_value_avg = st.mean(total_values)
    return total_values, total_value_avg


def runGA(Agents, Tasks, pop_size, max_ite, Pr_mu): 
    '''
    Genetic algorithm for task allocation.
    Input:
        Agents: [list] All agents' ids
        Tasks: [list] All tasks' ids
        pop_size: [int] The size of the population
        max_ite: [int] The maximum number of iterations
        Pr_mu: [float] (0,1) The probability of mutation
    Output:
        selected: [list 2D] selected tasks' ids by each agent 
        values: [list] function value for each agent
        total_value: [float] total function value i.e. sum of all individual agent's function value
        dt: [float] consuming time /sec
        consensus_steps: [int] the number of consensus steps
        n_evs: [int] the number of function evaluations
    '''
    # initialization
    start_time = time.time()
    consensus_steps = 0 # count the total consensus times
    n_evs = 0 # count the number of utility function evaluations
    
    len_chrom = len(Tasks)
    values = [0 for a in Agents] # Record the final function values for each agent.    
#    remained = [list(Tasks) for a in Agents] # start from ground set   
    # Create the empty allocation for each agent start from empty set. The first number is agent id.
    selected = [[a] for a in Agents] 
    population = [] # 2D list containing all individuals.
    pop_values = [] # 2D List containing the values of each agent for all individuals.
    pop_utilities = [] # 1D List containing the total function utilities for all individuals.
    
    # initial population
    for pop in range(pop_size):
        # Generate a random individual.
        individual = list(np.random.randint(0, len(Agents+1), len_chrom))
        population.append(individual)
        # Get the utility of the allocation result for the current individual.
        values_ind, total_value, evs = getUtility(individual, Agents, Tasks)
        n_evs += evs
        pop_values.append(values_ind) # Record the value of each agent for the current individual.
        pop_utilities.append(total_value) # Record the total value for the current individual.
            
    # Iterations
    for ite in range(max_ite):
        
        # Randomly select two ids of individuals as parents for Crossover
        p1_id, p2_id = np.random.randint(0, pop_size, 2)
        # Get two offspring ids.
        o1, o2 = crossover_1point(population[p1_id], population[p2_id]) 
        # Mutation for o1 and o2 with probability Pr_mu.
        if np.random.choice([0, 1],p=[1-Pr_mu,Pr_mu]):
            o1 = mutation_swap(o1)
        if np.random.choice([0, 1],p=[1-Pr_mu,Pr_mu]):
            o2 = mutation_swap(o2)
        # Get values of o1 and o2 and add them to the population.
        for oi in [o1,o2]:
            # Get function utilities for oi.
            values_ind, total_value, evs = getUtility(oi, Agents, Tasks)
            n_evs += evs
            # Put o1 into population.
            population.append(oi)
            pop_values.append(values_ind) # Record the value of each agent for o1.
            pop_utilities.append(total_value) # Record the total value for o1.
        # Remove two individuals with the last two highest total value 
        # to keep the size of population as initial value.
        for oi in [o1,o2]:
            min_id = np.argmin(pop_utilities)
            population.pop(min_id)
            pop_values.pop(min_id) # Record the value of each agent for o1.
            pop_utilities.pop(min_id)
    
    # Get the best individual with the highest function utility.        
    max_ind_id = np.argmax(pop_utilities)
    # Get the allocation result.
    for i in range(len_chrom):
        selected[population[max_ind_id][i]].append(i) 
    values = pop_values[max_ind_id]
    total_value = pop_utilities[max_ind_id]
    
    end_time = time.time()
    dt = end_time - start_time 
    
    return selected, values, total_value, dt, consensus_steps, n_evs
    

def runGA_Tune(Agents, Tasks, pop_size, max_ite, Pr_mu): 
    '''
    Tune the parameters for the Genetic Algorithm. (Show utility evolution as iteration increases.)
    Input:
        Agents: [list] All agents' ids
        Tasks: [list] All tasks' ids
        pop_size: [int] The size of the population
        max_ite: [int] The maximum number of iterations
        Pr_mu: [float] (0,1) The probability of mutation
    Output:
        selected: [list 2D] selected tasks' ids by each agent 
        values: [list] function value for each agent
        total_value: [float] total function value i.e. sum of all individual agent's function value
        dt: [float] consuming time /sec
        consensus_steps: [int] the number of consensus steps
        n_evs: [int] the number of function evaluations
        max_utilities: [list] The maximum function utilities in each iteration.
    '''
    # initialization
    start_time = time.time()
    consensus_steps = 0 # count the total consensus times
    n_evs = 0 # count the number of utility function evaluations
    evs_ite = [] # Record the evs wrt different iterations.
    max_utilities = [] # Record the maximum function utility in each iteration.
    
    len_chrom = len(Tasks)
    values = [0 for a in Agents] # Record the final function values for each agent.    
#    remained = [list(Tasks) for a in Agents] # start from ground set   
    # Create the empty allocation for each agent start from empty set. The first number is agent id.
    selected = [[a] for a in Agents] 
    population = [] # 2D list containing all individuals.
    pop_values = [] # 2D List containing the values of each agent for all individuals.
    pop_utilities = [] # 1D List containing the total function utilities for all individuals.
    
    # initial population
    for pop in range(pop_size):
        # Generate a random individual.
        individual = list(np.random.randint(0, len(Agents), len_chrom))
        population.append(individual)
        # Get the utility of the allocation result for the current individual.
        values_ind, total_value, evs = getUtility(individual, Agents, Tasks)
        n_evs += evs
        pop_values.append(values_ind) # Record the value of each agent for the current individual.
        pop_utilities.append(total_value) # Record the total value for the current individual.
        # Record the maximum function utility.
#        max_utilities.append(np.max(pop_utilities))
            
    # Iterations
    for ite in range(max_ite):
        
        # Randomly select two ids of individuals as parents for Crossover
        p1_id, p2_id = np.random.randint(0, pop_size, 2)
        # Get two offspring ids.
        o1, o2 = crossover_1point(population[p1_id], population[p2_id]) 
        # Mutation for o1 and o2 with probability Pr_mu.
        if np.random.choice([0, 1],p=[1-Pr_mu,Pr_mu]):
            o1 = mutation_swap(o1)
        if np.random.choice([0, 1],p=[1-Pr_mu,Pr_mu]):
            o2 = mutation_swap(o2)
        # Get values of o1 and o2 and add them to the population.
        for oi in [o1,o2]:
            # Get function utilities for oi.
            values_ind, total_value, evs = getUtility(oi, Agents, Tasks)
            n_evs += evs
            # Put o1 into population.
            population.append(oi)
            pop_values.append(values_ind) # Record the value of each agent for o1.
            pop_utilities.append(total_value) # Record the total value for o1.
        # Remove two individuals with the last two highest total value 
        # to keep the size of population as initial value.
        for oi in [o1,o2]:
            min_id = np.argmin(pop_utilities)
            population.pop(min_id)
            pop_values.pop(min_id) # Record the value of each agent for o1.
            pop_utilities.pop(min_id)
        # Record the maximum function utility.
        max_utilities.append(np.max(pop_utilities))
        # Record the number of function evaluations wrt different iterations.
        evs_ite.append(n_evs)
    
    # Get the best individual with the highest function utility.        
    max_ind_id = np.argmax(pop_utilities)
    # Get the allocation result.
    for i in range(len_chrom):
        selected[population[max_ind_id][i]].append(i) 
    values = pop_values[max_ind_id]
    total_value = pop_utilities[max_ind_id]
    
    end_time = time.time()
    dt = end_time - start_time 
    
    return selected, values, total_value, dt, consensus_steps, n_evs, max_utilities, evs_ite
    
      
# =============================================================================
# 
# =============================================================================
print("----- ga.py is loaded -----")