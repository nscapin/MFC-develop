#!/bin/sh
#SBATCH -p normal
#SBATCH -J simulation
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH -t 24:00:00
#SBATCH -o simulation.o%j
#SBATCH -e simulation.o%j
#SBATCH --mail-type=all
#SBATCH --mail-user=
echo MFC - Cases - 2D_bubbles: $SLURM_JOB_NAME.o$SLURM_JOB_ID
echo Description: $SLURM_JOB_ID executed on $SLURM_NTASKS processor'(s)'. The
echo '            ' command-line output information may be found below.
echo Start-date: `date +%D`
echo Start-time: `date +%T`
echo
echo
echo '================================ Terminal Output ==============================='
echo
t_start=$(date +%s)
LD_LIBRARY_PATH="$LD_LIBRARY_PATH:/autofs/nccs-svm1_envoy_od/aradhakr34/MFC-develop/src/../build/common/build/lib/" mpirun /autofs/nccs-svm1_envoy_od/aradhakr34/MFC-develop/src/../build/___current___/build/bin/simulation
t_stop=$(date +%s)
echo
echo '================================================================================'
echo
echo
echo End-date: `date +%D`
echo End-time: `date +%T`
echo
echo Total-time: $(expr $t_stop - $t_start)s
rm -f simulation.inp
rm -f simulation.sh