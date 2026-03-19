# Module_HDF5
This is a Fortran wrapper for writing data into .h5 file, which is particularly created for CFD.

## Description of the files
**mod_h5_utility.F90**\
The file contains a module mod_h5_utility, serving as a wrapper, which can be directly used in any Fortran code.

**main.F90**\
Various example files are provided under example folders main_file*/. One example for serial mode is provided, however, the parallel mode examples can also be executed in serial mode with necessary changes in the main file and job submission file.

**makefile**\
This file can be used to compile the code.

**job_submit_parallel.sh**\
In each example folder for parallel mode, SLURM job submission file is provided.

**reading_h5_data.m**\
This is an example Matlab file to read .h5 file data.


## Steps to compile and run an example
1. Copy mod_h5_utility.F90 to a directory, say test_dir/.
2. Copy any main.F90 and corresponding job_submit_parallel.sh files to test_dir/.
3. Copy makefile to the same directory.
4. In cluster, first load modules as described in the makefile. For example,\
   module load gnu openmpi
5. If particular versions of the modules are used while compiling, update the modules in job_submit_parallel.sh.
6. To compile, type: make
7. To clean directory, type: make clean
8. To submit job through SLURM, type: sbatch \<name of the executable\>.

## Output files
**.h5 files**  : As defined in the main.F90 files.

**.txt files** : output and error files as defined in job_submit_parallel.sh.
