# MRTA
A serial of multi-robot task allocation algorithms for performance comparison through simulations.

Multi-agent systems, task/resource allocation, submodular optimisation.

## Basic infomation

Programming Language:
    Python 3.12+


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

Algorithms to be updated: TGTA, TBTA, T3A, TTBTA, STBTA.

Algorithm files are located in "func/algs/".


## How to Run Simulations
main.py

1. Set simulation scenario parameters such as number of tasks and agents, number of Monte Carlo runs etc.
2. Enable the algorithms for simulation by setting "init.XXX.en" to 1.
3. Select the mission scenario: montecarlo, variance, tradeoff (Only select one).
4. Adjust plot functions in "__main__".
5. Run the entire project.


## Contributors:
Welcom to contribute to this repo. You can create a branch and make updates to codes. 
If you would like to merge your updates into the main branch, please create a pull request. We will review your updates and merge necessary functionalities.

## TODO list
- Use multi-thread and timer for the progress bar updates.
- Update all threshold related algorithms with buffers (refer to DTTA).


## Citations

### DTTA
Li, Teng, Hyo-Sang Shin, and Antonios Tsourdos. "Efficient Decentralised Parallel Task Allocation for Multiple Robots." IEEE Transactions on Robotics (2025).

### DSTA
Shin, H.S., Li, T., Lee, H.I. and Tsourdos, A., 2022. Sample greedy based task allocation for multiple robot systems. Swarm Intelligence, 16(3), pp.233-260.


## Copyrights

** non-commercial use only **<br>
@author: Teng Li <br>
lt.uk@outlook.com <br>
United Kingdom <br>
All Rights Reserved <br>



