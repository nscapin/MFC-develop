#!/bin/bash

# r -> Number of resource sets (RS) per host (node)
# c -> Number of physical cores (CPUs) per RS
# a -> Number of MPI ranks per RS
# g -> Number of GPUs per RS

# -r6 -a7 -c7 -g1 
## Use whole node
## One resource set per GPU (g=1)
## 42/6 = 7 where 42 = number of total cores (c)
## (a == c) usually

# jsrun -r6 -a7 -c7 -g1 ../../src/pre_process_code/pre_process
# jsrun -r6 -a7 -c7 -g1 ../../src/simulation_code/simulation

jsrun -r1 -a1 -c1 -g1 ../../src/pre_process_code/pre_process
jsrun -r1 -a1 -c1 -g1 ../../src/simulation_code/simulation
# jsrun -r1 -a1 -c1 -g1 nvprof ../../src/simulation_code/simulation
