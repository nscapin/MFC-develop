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
export PGI_ACC_NOTIFY=1

jsrun -r1 -a1 -c1 -g1 ../../src/pre_process_code/pre_process

# jsrun -r1 -a1 -c1 -g1 ../../src/simulation_code/simulation
# jsrun -r1 -a1 -c1 -g1 ../../src/simulation_code/simulation
# jsrun -r1 -a1 -c1 -g1 nvprof --analysis-metrics -o output.nvvp  -f ../../src/simulation_code/simulation

# jsrun -r1 -a1 -c1 -g1 ncu --kernel-regex s_weno_alt_307_gpu --launch-skip 3 --launch-count 1 --set full -o output.prof -f  ../../src/simulation_code/simulation


# jsrun -r1 -a1 -c1 -g1 nsys profile -o output-sys.prof --stats=true -t openacc  ../../src/simulation_code/simulation

# jsrun -r1 -a1 -c1 -g1 ../../src/simulation_code/simulation
