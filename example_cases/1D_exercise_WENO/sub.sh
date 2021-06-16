#!/bin/bash
#SBATCH -n 1
#SBATCH --ntasks-per-node=1
#SBATCH -p GPU-shared
#SBATCH -t 1:00:00
#SBATCH --gres=gpu:v100-32:1

set -x

####SBATCH -C PERF

mpirun -n 1 ../../src/pre_process_code/pre_process
# mpirun -n 1 ../../src/simulation_code/simulation
mpirun -n 1 ncu --launch-skip 2 --launch-count 15 --set full -o output-cu.prof -f   ../../src/simulation_code/simulation
