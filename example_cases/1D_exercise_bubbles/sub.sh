#!/bin/bash
#SBATCH --ntasks-per-node=5
#SBATCH -p GPU-shared
#SBATCH -t 0:05:00
#SBATCH --gpus=v100-32:1
## #SBATCH --constraint=PERF

mpirun -n 1 ../../src/pre_process_code/pre_process
mpirun -n 1 nvprof --export-profile output.nvvp -f --analysis-metrics  ../../src/simulation_code/simulation
