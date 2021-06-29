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

#SBATCH --gpu-bind=map_gpu:0,1,2,3,4,5,6,7 to bind task 0 to GPU0, task 1 to GPU1, task 2 to GPU2, etc. See the slurm documentation for srun:  https://slurm.schedmd.com/srun.html
# You can also bind several tasks (4 in the following example) to the same GPU with the same gpu-bind parameter:
# #SBATCH --gpu-bind=map_gpu:0,0,0,0,1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4,5,5,5,5,6,6,6,6,7,7,7,7
# I think you can shorten the previous line like this:   #SBATCH --gpu-bind=map_gpu:0*4,1*4,2*4,3*4,4*4,5*4,6*4,7*4

#SBATCH --gpu-mps

nvidia-cuda-mps-control  -d
nvidia-cuda-mps-server

# export PGI_ACC_NOTIFY=3
# export NV_ACC_NOTIFY=3
# export NV_ACC_TIME=1
# export NV_ACC_DEBUG=1

  # set flag for this example

# BINARY=simulation
# EXEC=../../src/simulation_code/${BINARY}
# OUT=hpctoolkit-${BINARY}-acc 
# STRUCT_FILE=$BINARY.hpcstruct

# rm -r hpctoolkit*
# rm -rf ${OUT}.m ${OUT}.d $STRUCT_FILE

# jsrun -r1 -a1 -c1 -g1 ../../src/pre_process_code/pre_process
# ${HPCTOOLKIT_LULESH_ACC_LAUNCH} hpcrun -o $OUT.m -e REALTIME -e gpu=nvidia -t ${EXEC}
# hpcstruct -o $STRUCT_FILE ${EXEC}
# hpcstruct --gpucfg no $OUT.m
# hpcprof -S $STRUCT_FILE -o $OUT.d $OUT.m

# jsrun -r1 -a1 -c1 -g1 hpcprof -S simulation.hpcstruct -I ../../src/simulation_code/+ hpctoolkit-simulation-measurements-*


# jsrun -r1 -a7 -c7 -g1 ../../src/pre_process_code/pre_process
# jsrun -r1 -a7 -c7 -g1 ../../src/simulation_code/simulation

# jsrun -r6 -a7 -c7 -g1 ../../src/pre_process_code/pre_process
# jsrun -r6 -a7 -c7 -g1 ../../src/simulation_code/simulation

# jsrun -r1 -a1 -c1 -g1 ../../src/pre_process_code/pre_process
# jsrun -r1 -a1 -c1 -g1 ../../src/simulation_code/simulation

# jsrun -r1 -a2 -c2 -g1 ../../src/pre_process_code/pre_process
# jsrun -r1 -a2 -c2 -g1 ../../src/simulation_code/simulation

# jsrun -r1 -a1 -c1 -g1 nsys profile -o output-sys.prof --stats=true -t openacc,nvtx  --force-overwrite true ../../src/simulation_code/simulation

mpirun -n 4 ../../src/pre_process_code/pre_process
# mpirun -n 1 ../../src/simulation_code/simulation

# mpirun -n 2 ../../src/pre_process_code/pre_process
# mpirun -n 2 ../../src/simulation_code/simulation

# mpirun -n 4 ../../src/pre_process_code/pre_process
# mpirun -n 4 ../../src/simulation_code/simulation

# mpirun -n 2 nsys profile -o output-sys.prof_%q{OMPI_COMM_WORLD_RANK} --stats=true -t openacc,nvtx  --force-overwrite true ../../src/simulation_code/simulation
mpirun -n 4 nsys profile -o output-sys.prof --stats=true -t openacc,nvtx  --force-overwrite true ../../src/simulation_code/simulation


# mpirun -n 1 ncu --launch-skip 2 --launch-count 15 --set full -o output-cu.prof -f   ../../src/simulation_code/simulation
# mpirun -n 1 ncu --launch-skip 2 --launch-count 15 --set full -o output-cu.prof -f --nvtx  ../../src/simulation_code/simulation


## Profile 
# jsrun -r1 -a7 -c7 -g1 ../../src/pre_process_code/pre_process
# jsrun -r1 -a7 -c7 -g1 nsys profile -o output-sys.prof_%q{OMPI_COMM_WORLD_RANK} --stats=true -t openacc,nvtx  --force-overwrite true ../../src/simulation_code/simulation

# jsrun -r6 -a7 -c7 -g1 ../../src/pre_process_code/pre_process
# jsrun -r6 -a7 -c7 -g1 nsys profile -o output-sys.prof_%q{OMPI_COMM_WORLD_RANK} --stats=true -t openacc,nvtx  --force-overwrite true ../../src/simulation_code/simulation

# jsrun -r1 -a1 -c1 -g1 ../../src/simulation_code/simulation
# jsrun -r1 -a1 -c1 -g1 nvprof --analysis-metrics -o output.nvvp  -f ../../src/simulation_code/simulation

# jsrun -r1 -a1 -c1 -g1 ncu --kernel-regex s_weno_alt_177_gpu --launch-skip 3 --launch-count 9 --set full -o output.prof -f  ../../src/simulation_code/simulation

# jsrun -r1 -a1 -c1 -g1 ncu --launch-skip 2 --launch-count 5 --set full -o output.prof -f  ../../src/simulation_code/simulation


# jsrun -r1 -a7 -c7 -g1 nsys profile -o output-sys.prof_%q{OMPI_COMM_WORLD_RANK} --stats=true -t openacc  --force-overwrite true ../../src/simulation_code/simulation

# jsrun -r1 -a1 -c1 -g1 ../../src/simulation_code/simulation
