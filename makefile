# NOTE: Make sure that for the commands followed by tabs, the tabs are actually tabs by making special character visible in the text editor.
# Otherwise error is popped as:  *** missing separator.  Stop.

# This file should only be named: makefile (sensitive to case, If doesn't work, change to MakeFile or so)
# Run this file by typing: make  (only first target is executed (also probably those which are used as dependent with the first target)
# To execute a particular target, type: make target_name (ex. make clean)

# FOR SERIAL I/O
#COMPILER=h5fc # for gfortran with links for hdf5 modules. Use h5pfc for parallel I/O.
## one can execute: h5pfc -show :to see the wrapper options
# COMPILER=ifort -I/opt/rcc/intel/include -L/opt/rcc/intel/lib64 -lhdf5hl_fortran -lhdf5_hl -lhdf5_fortran -lhdf5
# COMPILER=gfortran -I/opt/rcc/gnu/include -L/opt/rcc/gnu/lib64 -lhdf5hl_fortran -lhdf5_hl -lhdf5_fortran -lhdf5
## For ifort use (with specific versions)    : module load intel
## For gfortran use (with specific versions) : module load gnu

# FOR PARALLEL I/O
COMPILER=h5pfc
## one can execute: h5pfc -show :to see the wrapper options
# COMPILER=mpif90 -I/opt/rcc/gnu/openmpi/include -L/opt/rcc/gnu/openmpi/lib64 -lhdf5hl_fortran -lhdf5_hl -lhdf5_fortran -lhdf5
# COMPILER=mpif90 -I/opt/rcc/intel/openmpi/include -L/opt/rcc/intel/openmpi/lib64 -lhdf5hl_fortran -lhdf5_hl -lhdf5_fortran -lhdf5
## For ifort use (with specific versions)    : module load intel openmpi
## For gfortran use (with specific versions) : module load gnu openmpi


# prog rule
prog: code_1.o code_2.o
	$(COMPILER) mod_h5_utility.o main.o -o execute_file

# code_1.o rule
code_1.o: mod_h5_utility.F90
	$(COMPILER) -c mod_h5_utility.F90

code_2.o: main.F90	
	$(COMPILER) -c main.F90


# clean rule     
clean: 
	rm -f *.o *.mod # execute_file  