#!/bin/bash
#SBATCH -N 2
#SBATCH -p GPU
#SBATCH -t 00:10:00
#SBATCH --gres=gpu:8

# rm -f ex.o*
# mpif90 -acc ex.f90 -o ex.o_f
# mpirun ex.o_f

nvidia-cuda-mps-control  -d
nvidia-cuda-mps-server

# export PGI_ACC_NOTIFY=3
# export NV_ACC_NOTIFY=3
# export NV_ACC_TIME=1
# export NV_ACC_DEBUG=1

mpirun ../../src/pre_process_code/pre_process
mpirun ../../src/simulation_code/simulation
