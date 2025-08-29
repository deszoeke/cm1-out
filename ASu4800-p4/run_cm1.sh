#!/bin/bash
#SBATCH --job-name=cm1_par_test
#SBATCH --output=cm1_par_test.o
#SBATCH --error=cm1_par_test.e
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=96
#SBATCH --time=40:00:00
#SBATCH --partition=ceoas
#SBATCH -w chen

# Set number of OpenMP threads
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

# Run the executable
./cm1.exe >& cm1.print.out

