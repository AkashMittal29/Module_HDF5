!###########################################################
! EXAMPLE CODE PARALLEL MODE, ADVANCED DATASET 02: Writing only a small region form the global data
! By Akash Kumar Mittal
! Feb 2026
! Mechanical and Aerospace Engineering
! Florida State University, USA
!###########################################################

PROGRAM main
    !###########################################################
    ! Using various modules
    USE hdf5
    USE mpi
    USE ISO_FORTRAN_ENV, ONLY: REAL64, INT64, error_unit
    IMPLICIT NONE

    !###########################################################
    ! Declaring variables
    INTEGER :: npx, npy ! number of processors along each dimension
    INTEGER :: i, j, k, nx, ny
    INTEGER :: mpierr, mpirank, mpiprocessors, ierror, mp_cart, color, mpi_comm, error
    REAL(kind=8), ALLOCATABLE, DIMENSION(:,:) ::x
    INTEGER(kind=HSIZE_T) :: file_slab_start(3), file_slab_end(3), file_slab_count(3) 
    INTEGER(kind=HSIZE_T) :: mem_slab_start(2), mem_slab_count(2), mem_stride(2), mem_block(2)
    INTEGER(kind=hid_t)   :: file_id, space_id, dataset_id, plist_id, memspace_id, transfer_plist_id
    INTEGER(kind=HSIZE_T) :: dims_current(3), dims_chunk(3), dims_mem(2)
    INTEGER(kind=8)       :: n_objects_open

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

    IF(mpirank==0) THEN
        ! Assigning x
        k = 1
        DO j=1,ny
        DO i=1,nx
            x(i,j) = k*1.0
            k=k+1
        END DO
        END DO

        ! File hyperslab
        file_slab_start = (/1,1,1/)
        file_slab_end   = (/3,2,1/)
        file_slab_count = file_slab_end-file_slab_start+1 ! (3,2)
        
        ! Memory hyperslab
        mem_slab_start = (/3,5/)
        mem_slab_count = file_slab_count(1:2)
        mem_stride = 1
        mem_block  = 1

    ELSE IF(mpirank==1) THEN
        ! Assigning x
        k = 1
        DO j=1,ny
        DO i=1,nx
            x(i,j) = k*10.0
            k=k+1
        END DO
        END DO

        ! File hyperslab
        file_slab_start = (/4,1,1/)
        file_slab_end   = (/7,2,1/)
        file_slab_count = file_slab_end-file_slab_start+1 ! (4,2)
        
        ! Memory hyperslab
        mem_slab_start = (/1,5/)
        mem_slab_count = file_slab_count(1:2)
        mem_stride = 1
        mem_block  = 1

    ELSE IF(mpirank==2) THEN
        ! Assigning x
        k = 1
        DO j=1,ny
        DO i=1,nx
            x(i,j) = k*100.0
            k=k+1
        END DO
        END DO

        ! File hyperslab
        file_slab_start = (/1,3,1/)
        file_slab_end   = (/3,5,1/)
        file_slab_count = file_slab_end-file_slab_start+1 ! (3,3)
        
        ! Memory hyperslab
        mem_slab_start = (/3,1/)
        mem_slab_count = file_slab_count(1:2)
        mem_stride = 1
        mem_block  = 1

    ELSE IF(mpirank==3) THEN
        ! Assigning x
        k = 1
        DO j=1,ny
        DO i=1,nx
            x(i,j) = k*(-1.0)
            k=k+1
        END DO
        END DO

        ! File hyperslab
        file_slab_start = (/4,3,1/)
        file_slab_end   = (/7,5,1/)
        file_slab_count = file_slab_end-file_slab_start+1 ! (4,3)
        
        ! Memory hyperslab
        mem_slab_start = (/1,1/)
        mem_slab_count = file_slab_count(1:2)
        mem_stride = 1
        mem_block  = 1

    END IF

    !###########################################################
    ! Initializing hdf5 library. 
    CALL h5open_f(error) 

    ! Setting property list for accessing file with parallel I/O
    CALL h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
    CALL h5pset_fapl_mpio_f(plist_id, mpi_comm, MPI_INFO_NULL, error)

    ! Opening file
    CALL h5fcreate_f('grid_check.h5', H5F_ACC_TRUNC_F, file_id, error, access_prp=plist_id) 
    CALL h5pclose_f(plist_id, error)

    ! Creating file space id
    dims_current = (/7,5,1/)
    CALL h5screate_simple_f(3, dims_current, space_id, error)
       
    ! Creating property list for chunking
    dims_chunk = dims_current
    CALL h5pcreate_f(H5P_DATASET_CREATE_F, plist_id, error)
    CALL h5pset_chunk_f(plist_id, 3, dims_chunk, error) 
    
    ! Creating dataset id
    CALL h5dcreate_f(file_id, 'x', H5T_NATIVE_DOUBLE, space_id, dataset_id, error, plist_id) 

    ! Creating memory space
    dims_mem = (/nx,ny/)
    CALL h5screate_simple_f(2, dims_mem, memspace_id, error) ! Is used for creating memory hyperslab and writing data

    ! Selecting the memory hyperslab. 
    CALL h5sselect_hyperslab_f(memspace_id, H5S_SELECT_SET_F, mem_slab_start-1, mem_slab_count, error, &
                             & stride=mem_stride, block=mem_block)

    ! Creating dataset transfer mode property list
    CALL h5pcreate_f(H5P_DATASET_XFER_F, transfer_plist_id, error) ! if dxpl_mpio (next line) is not set: in serial: no effect, in parallel: default is independent mode.
    CALL h5pset_dxpl_mpio_f(transfer_plist_id, H5FD_MPIO_COLLECTIVE_F, error) 

    ! selecting the file hyperslab (each processor rank will call this with their locaal self%start and self%count)
    CALL h5sselect_hyperslab_f(space_id, H5S_SELECT_SET_F, file_slab_start-1, file_slab_count, error) ! from the given data space, it selects a given region.
         
    ! Writing data
    CALL h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, x, dims_mem, &
                  & error, mem_space_id=memspace_id, file_space_id=space_id, xfer_prp=transfer_plist_id)

    !###########################################################
    ! Closing h5 resources
    CALL h5sclose_f(space_id, error)
    CALL h5pclose_f(plist_id, error)
    CALL h5pclose_f(transfer_plist_id, error)
    CALL h5sclose_f(memspace_id, error)
    CALL h5dclose_f(dataset_id, error)
    CALL h5fget_obj_count_f(file_id, H5F_OBJ_ALL_F, n_objects_open, error)
    WRITE(*,*) 'n_objects_open: ', n_objects_open
    IF(n_objects_open==1) THEN ! will happen only if all the child objects of the file are closed. No. of objects = 1 is the file itself.
        CALL h5fclose_f(file_id, error) ! closing file_id
    END IF
    CALL h5close_f(error)

    !###########################################################
    ! Finalizing mpi
    CALL MPI_FINALIZE(mpierr)

END PROGRAM MAIN


