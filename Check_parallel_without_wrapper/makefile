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

# NOTE: since mod_h5_utility.f90 uses module mpi even if it is not utilized in seriall mode, 
# the code has to be compiled after loading openmpi.

# FOR PARALLEL I/O
COMPILER=h5pfc
OPTIONS=-DDEBUGhdf5
# -DDEBUGhdf5 is used for defining preprocessor directive by the name DEBUGhdf5 in mod_h5_utility. (Not used in this check case.)

## one can execute: h5pfc -show :to see the wrapper options
# COMPILER=mpif90 -I/opt/rcc/gnu/openmpi/include -L/opt/rcc/gnu/openmpi/lib64 -lhdf5hl_fortran -lhdf5_hl -lhdf5_fortran -lhdf5
# COMPILER=mpif90 -I/opt/rcc/intel/openmpi/include -L/opt/rcc/intel/openmpi/lib64 -lhdf5hl_fortran -lhdf5_hl -lhdf5_fortran -lhdf5
## For ifort use (with specific versions)    : module load intel openmpi
## For gfortran use (with specific versions) : module load gnu openmpi


# prog rule
prog: code_2.o
	$(COMPILER) $(OPTIONS) main.o -o execute_file

code_2.o: main.F90	
	$(COMPILER) $(OPTIONS) -c main.F90


# clean rule     
clean: 
	rm -f *.o *.mod # execute_file  