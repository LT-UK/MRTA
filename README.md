# MRTA
A serial of multi-robot task allocation algorithms for performance comparison through simulations.

Multi-agent systems, task/resource allocation, submodular optimisation.

## Basic infomation

This Python project is originally corresponding to the PhD thesis:
    
    Title: Efficient Decentralised Task Allocation for Multiple Aerial Robots
    
    Academic year: 2017 -- 2020
    
    Author: Teng Li
    
    Supervisors: Prof. Hyo-Sang Shin
                 Prof. Antonios Tsourdos
    
    Affiliation: Centre for AI, Robotics and Space, FEAS,
                 Cranfield University, 
                 United Kingdom

Programming Language:
    Python 3.12


The simulation scenario is based on a multi-target surveillance mission using multiple UAVs where the utility function is submodular. 

The proposed algorithms can provide a theoretical optimality guarantee. They can achieve comparable solution quality but are more efficient than benchmark algorithm.

Please note that the algorithms are upgraded time to time according to reviewers comments. 
Some of the codes/comments are outdated, users of this repo may contact the author if you enconter any bugs.


## Algorithm list:

- GA:       Genetic Algorithm
- SGA:      Sequencial Greedy Algorithm
- CBBA:     Consensus Based Bundle Algorithm
- TGTA:     Truncation Greedy Task Allocation
- DTTA:     Decreasing Threshold Task Allocation
- TBTA:     Threshold Bundle Task Allocation
- T3A:      Truncation Threshold Task Allocation
- TTBTA:    Truncation Threshold Bundle Task Allocation
- DSTA:     Decentralised Sample based Task Allocation
- STTA:     Sample Threshold Task Allocation
- STBTA:    Sample Threshold Bundle Task Allocation
- Auction_xx: Auction based algorithms

The prefix 'L' letter represents 'Lazy' version of these algorithm in this project.

Algorithms proposed by the author: TGTA, DTTA, TBTA, T3A, TTBTA, DSTA, STTA, STBTA and their 'Lazy' versions.

## Copyrights

** non-commercial use only **<br>
@programmer: Teng Li <br>
lt.uk@outlook.com <br>
United Kingdom <br>
All Rights Reserved <br>
