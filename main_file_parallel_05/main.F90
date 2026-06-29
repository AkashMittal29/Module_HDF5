!###########################################################
! EXAMPLE CODE PARALLEL MODE 05
! By Akash Kumar Mittal
! June 2026
! Mechanical and Aerospace Engineering
! Florida State University, USA

! In this example, two datasets are written at different frequencies and 
! the restart is checked.
! Note: The datasets writing to the same file should have identical restart flag ( 1 or 0).
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
    INTEGER, PARAMETER :: real = selected_real_kind(15, 307) ! Here it is used form module data.
    REAL(kind=real) :: a1(2), a2(2), v1(2), v2(2) ! Sample data
    TYPE(h5_dataset_type) :: iterdata1, valuedata1, iterdata2, valuedata2
    INTEGER :: iter, restart1, restart2, save_freq1, save_freq2, i
    INTEGER :: mpierr, mpirank
    INTEGER :: file_slab_start(1), file_slab_end(1) ! file hyperslab indices for each rank

    !###########################################################
    ! Initializing mpi (message passing interface)
    CALL MPI_INIT(mpierr)
    CALL MPI_COMM_RANK(MPI_COMM_WORLD, mpirank, mpierr)

    !###########################################################
    ! Assigning local variables
    a1 = (/1.0, 2.0/)+2*mpirank ! Sample data
    a2 = (/5.0, 6.0/)+2*mpirank ! Sample data
    v1 = a1
    v2 = a2
    iter       = 7  ! Change this value (eg., 7) to set restart iteration after running for iter=1, restart1=0, restart2=0.    
    restart1   = 1     
    restart2   = restart1 ! Since the datasets write into a common file
    save_freq1 = 3 ! Save after every 3 iteration
    save_freq2 = 5 ! Save after every 5 iterations

    !###########################################################
    ! Setting file hyperslab for each rank. Here, each rank writes 2 elements in the file dataset.
    file_slab_start = (/mpirank*2+1/)
    file_slab_end   = (/mpirank*2+2/)

    !###########################################################
    ! Setting datatype for int and float for hdf5; This should be called once.
    CALL h5_utility_set_datatype('H5T_STD_I64LE','H5T_NATIVE_DOUBLE')
    ! Setting parallel I/O capability for mod_h5_utility
    h5_utility_mpi = .TRUE. 

    !###########################################################
    ! Creation of datasets; This should be called once before writing the data.
    ! Note: the following datasets write into a common file. Hence, they all shall have identical restart value (restart1=restart2).
    ! data 1: /iter (important for a restart) (DATA SHOULD BE INTEGER TYPE, since it will be read as integer)
    CALL iterdata1%create( global_array_size         = (/1_INT64/),     &
                         & data_address              = 'field.h5/',     &
                         & dataset_names             = ['iter1'],       &
                         & restart                   = restart1,        &
                         & restart_step_value        = INT(iter,INT64), &
                         & restart_step_data_address = '/iter1',        &
                         & data_type                 = 'int' )

    ! data 2: /grid/x, /grid/y (both share same dataspace)
    CALL valuedata1%create( global_array_size = (/4_INT64/),        &
                        & data_address      = 'field.h5/',          &
                        & dataset_names     = ['value1'],           &
                        & restart           = restart1,             & 
                        & restart_ind       = iterdata1%rest_index, &
                        & slab_start_ind    = file_slab_start,      &
                        & slab_end_ind      = file_slab_end )

    ! data 1: /iter (important for a restart) (DATA SHOULD BE INTEGER TYPE, since it will be read as integer)
    CALL iterdata2%create( global_array_size         = (/1_INT64/),     &
                         & data_address              = 'field.h5/',     &
                         & dataset_names             = ['iter2'],       &
                         & restart                   = restart2,        &
                         & restart_step_value        = INT(iter,INT64), &
                         & restart_step_data_address = '/iter2',        &
                         & data_type                 = 'int' )

    ! data 2: /grid/x, /grid/y (both share same dataspace)
    CALL valuedata2%create( global_array_size = (/4_INT64/),        &
                        & data_address      = 'field.h5/',          &
                        & dataset_names     = ['value2'],           &
                        & restart           = restart2,             & 
                        & restart_ind       = iterdata2%rest_index, &
                        & slab_start_ind    = file_slab_start,      &
                        & slab_end_ind      = file_slab_end )



    !###########################################################
    ! Writing data every iteration
    DO i=iter,iter+10-1
        WRITE(error_unit,*) 'loop: ', i, ' ', mpirank
        WRITE(*,*) 'loop: ', i, ' ', mpirank

        ! Updating field variables
        v1 = a1+iter;
        v2 = a2+iter;

        ! Preparing append for hdf5 file. This sets the hyperslab, and extends the dataset by n_extend if 
        ! the file dataset size is exhausted along iteration dimension.
        IF(MOD(iter,save_freq1)==0) THEN
            ! Preparing append
            CALL iterdata1%prepare_next_append(n_extend=5)
            CALL valuedata1%prepare_next_append(n_extend=5)
            ! Wrtiting data
            CALL iterdata1%append((/(INT(iter,INT64))/),  data_index=1)
            CALL valuedata1%append(v1, data_index=1)
        END IF

        IF(MOD(iter,save_freq2)==0) THEN
            ! Preparing append
            CALL iterdata2%prepare_next_append(n_extend=5)
            CALL valuedata2%prepare_next_append(n_extend=5)
            ! Wrtiting data
            CALL iterdata2%append((/(INT(iter,INT64))/),  data_index=1)
            CALL valuedata2%append(v2, data_index=1)
        END IF
       
        iter = iter+1
        IF(h5_utility_mpi) CALL MPI_BARRIER(MPI_COMM_WORLD, mpierr)
    END DO

    !###########################################################
    ! Closing all the resources of hdf5. Should be included at the end of the program to avoid memory leak.
    CALL iterdata1%destructor()
    CALL valuedata1%destructor()
    CALL iterdata2%destructor()
    CALL valuedata2%destructor()

    !###########################################################
    ! Finalizing mpi
    CALL MPI_FINALIZE(mpierr)

END PROGRAM MAIN


