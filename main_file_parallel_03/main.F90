!###########################################################
! EXAMPLE CODE PARALLEL MODE 03: Writing only a small region form the global data
! By Akash Kumar Mittal
! Mar 2026
! Mechanical and Aerospace Engineering
! Florida State University, USA

! Run this example with number of ranks/processors = 4.
! In this example, each rank writes different memory hyperslab from a given 
! memory size to a portion in the file space. The mod_h5_utility module ensures
! that the size in memory hyperslab is equal to the size in file hyperslab.
!###########################################################

PROGRAM main
    !###########################################################
    ! Using various modules
    USE mod_h5_utility
    USE mpi
    USE ISO_FORTRAN_ENV, ONLY: REAL64, INT64, error_unit
    IMPLICIT NONE

    !###########################################################
    ! Declaring variables
    INTEGER, PARAMETER :: real = selected_real_kind(15, 307)
    INTEGER :: i, j, k, nx, ny
    INTEGER :: mpierr, mpirank, mpiprocessors, mpi_comm
    REAL(kind=8), ALLOCATABLE, DIMENSION(:,:) ::x
    INTEGER :: file_slab_start(2), file_slab_end(2), file_slab_count(2) 
    INTEGER :: mem_slab_start(2), mem_slab_count(2)
    TYPE(h5_dataset_type) :: griddata

    !###########################################################
    ! Initializing mpi (message passing interface)
    CALL MPI_INIT(mpierr)
    CALL MPI_COMM_RANK(MPI_COMM_WORLD, mpirank, mpierr)
    CALL MPI_COMM_SIZE(MPI_COMM_WORLD, mpiprocessors, mpierr)
    mpi_comm = MPI_COMM_WORLD

    !###########################################################
    ! In the following eg., each rank has unique file-hyperslab and memory-hyperslab.
    nx = 5; ny = 6;
    ALLOCATE(x(nx,ny))

    CALL h5_utility_set_datatype('H5T_STD_I64LE','H5T_NATIVE_DOUBLE')
    ! Setting parallel I/O capability for mod_h5_utility
    h5_utility_mpi = .TRUE.

    ! Manually setting data and hyperslabs
    IF(mpirank==0) THEN
        ! Assigning x
        k = 1
        DO j=1,ny
        DO i=1,nx
            x(i,j) = k*1.0_real
            k=k+1
        END DO
        END DO

        ! File hyperslab
        file_slab_start = (/1,1/)
        file_slab_end   = (/3,2/)
        file_slab_count = file_slab_end-file_slab_start+1 ! (3,2)
        
        ! Memory hyperslab
        mem_slab_start = (/3,5/)
        mem_slab_count = file_slab_count

    ELSE IF(mpirank==1) THEN
        ! Assigning x
        k = 1
        DO j=1,ny
        DO i=1,nx
            x(i,j) = k*10.0_real
            k=k+1
        END DO
        END DO

        ! File hyperslab
        file_slab_start = (/4,1/)
        file_slab_end   = (/7,2/)
        file_slab_count = file_slab_end-file_slab_start+1 ! (4,2)
        
        ! Memory hyperslab
        mem_slab_start = (/1,5/)
        mem_slab_count = file_slab_count

    ELSE IF(mpirank==2) THEN
        ! Assigning x
        k = 1
        DO j=1,ny
        DO i=1,nx
            x(i,j) = k*100.0_real
            k=k+1
        END DO
        END DO

        ! File hyperslab
        file_slab_start = (/1,3/)
        file_slab_end   = (/3,5/)
        file_slab_count = file_slab_end-file_slab_start+1 ! (3,3)
        
        ! Memory hyperslab
        mem_slab_start = (/3,1/)
        mem_slab_count = file_slab_count

    ELSE IF(mpirank==3) THEN
        ! Assigning x
        k = 1
        DO j=1,ny
        DO i=1,nx
            x(i,j) = k*(-1.0_real)
            k=k+1
        END DO
        END DO

        ! File hyperslab
        file_slab_start = (/4,3/)
        file_slab_end   = (/7,5/)
        file_slab_count = file_slab_end-file_slab_start+1 ! (4,3)
        
        ! Memory hyperslab
        mem_slab_start = (/1,1/)
        mem_slab_count = file_slab_count
    END IF

    !###########################################################
    CALL griddata%create( global_array_size = INT((/7,5/), INT64), &
                        & data_address      = 'field.h5/subgrid',  &
                        & dataset_names     = ['x'],               &
                        & restart           = 0,                   & 
                        & restart_ind       = 1_INT64,             &
                        & slab_start_ind    = file_slab_start,     &
                        & slab_end_ind      = file_slab_end,       &
                        & mem_size          = SHAPE(x),            &
                        & mem_start         = mem_slab_start,      &
                        & mpi_comm          = mpi_comm )

    !###########################################################
    CALL griddata%prepare_next_append(n_extend=1)
    CALL griddata%append(x, data_index=1)

    !###########################################################
    CALL griddata%destructor() 

    ! Finalizing mpi
    CALL MPI_FINALIZE(mpierr)

END PROGRAM MAIN


