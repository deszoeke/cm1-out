#!/bin/bash
#SBATCH --job-name=cm1_r2
#SBATCH --output=cm1_par_test.o
#SBATCH --error=cm1_par_test.e
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --time=48:00:00
#SBATCH --partition=ceoas

# Set number of OpenMP threads
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

# Run the executable
./cm1.exe >& cm1.print.out

