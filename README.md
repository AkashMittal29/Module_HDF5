# Module_HDF5
This is a Fortran wrapper for writing data into .h5 file, which is particularly created for CFD.

## Description of the files
mod_h5_utility.F90\
The file contains a module mod_h5_utility, serving as a wrapper, which can be directly used in any Fortran code.

main.F90\
Various example files are provided under folders main_file*/. One example for serial mode is provided, however, the parallel mode examples can also be executed in serial mode with necessary changes in the main file and job submission file.

makefile\
This file can be used to compile the code.

job_submit_parallel.sh\
In each example folder, SLURM job submission file is provided.


