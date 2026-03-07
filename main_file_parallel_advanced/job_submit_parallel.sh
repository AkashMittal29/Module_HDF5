#!/bin/bash

#SBATCH -J hdf5check
#SBATCH -n 12
#SBATCH -o out.txt
#SBATCH -t 00:05:00
#SBATCH -p mecfd_q
#SBATCH -e error.txt
#SBATCH -C "YEAR2022|YEAR2024"

cd $SLURM_SUBMIT_DIR/

chmod u+x execute_file

module purge
module load gnu/8.5.0
module load openmpi/4.1.0

srun --cpu-bind=cores ./execute_file
