!###########################################################
! EXAMPLE CODE PARALLEL MODE
! By Akash Kumar Mittal
! Feb 2026
! Mechanical and Aerospace Engineering
! Florida State University, USA
!###########################################################

! Defining a module data to generate an example of a typical dataset encountered in CFD. 
! It is used in the main program below.
MODULE data
    !###########################################################
    ! Declaring variables
    IMPLICIT NONE
    INTEGER :: iter = 0, ndim ! Module variables have saved attribute automatically.
    REAL(kind=8), ALLOCATABLE, DIMENSION(:,:) :: x_grid, y_grid, velocity, pressure
    REAL(kind=8), ALLOCATABLE, DIMENSION(:)   :: marker, x_grid_1d, y_grid_1d
    REAL(kind=8)    :: dx, dy
    INTEGER(kind=8) :: nxg, nyg, npointg ! Global dimensions
    INTEGER :: nx, ny, npoint          ! Dimensions local to each rank

    !###########################################################
    ! Defining subroutines of module data
    CONTAINS

        !###########################################################
        ! Assigns variables, allocates arrays, and defines grid
        SUBROUTINE define_grid(npx, npy, mpi_coord)
            INTEGER, INTENT(IN) :: npx, npy ! number of processors along each dimension
            INTEGER, DIMENSION(:), INTENT(IN) :: mpi_coord
            REAL(KIND=8) :: x_start, y_start
            INTEGER :: i

            ! Assigning variables
            dx      = 0.2
            dy      = 0.1
            nx      = 5
            ny      = 6
            npoint  = 9
            nxg     = npx*nx
            nyg     = npy*ny
            npointg = npx*npoint

            ! Allocating
            IF(.NOT. ALLOCATED(x_grid_1d)) ALLOCATE(x_grid_1d(nx))
            IF(.NOT. ALLOCATED(y_grid_1d)) ALLOCATE(y_grid_1d(ny))
            IF(.NOT. ALLOCATED(x_grid))    ALLOCATE(x_grid(nx,ny))
            IF(.NOT. ALLOCATED(y_grid))    ALLOCATE(y_grid(nx,ny))
            IF(.NOT. ALLOCATED(velocity))  ALLOCATE(velocity(nx,ny))
            IF(.NOT. ALLOCATED(pressure))  ALLOCATE(pressure(nx,ny))
            IF(.NOT. ALLOCATED(marker))    ALLOCATE(marker(npoint))

            ! Assigning grid
            x_start = 0.0 + mpi_coord(2)*nx*dx
            x_grid_1d = (/(x_start+i*dx, i=0,nx-1)/) ! Implied do loop
            y_start = 0.0 + mpi_coord(1)*ny*dy
            y_grid_1d = (/(y_start+i*dy, i=0,ny-1)/)
            DO i=1,nx
                x_grid(i,:) = x_grid_1d(i)
            END DO
            DO i=1,ny
                y_grid(:,i) = y_grid_1d(i)
            END DO
        END SUBROUTINE define_grid

        !###########################################################
        ! Initializes data
        SUBROUTINE initialize_data()
            velocity = 1.0
            pressure = 0.0
            marker   = 0.0
        END SUBROUTINE initialize_data

        !###########################################################
        ! Updates the arrays
        SUBROUTINE update_data()
            INTEGER :: i,j

            ! Assigning data by some operations
            DO j=1,ny,1
            DO i=1,nx,1
                velocity(i,j) = iter*0.01 + x_grid(i,j)**2 + y_grid(i,j)
                pressure(i,j) = iter*0.01 + x_grid(i,j)**3 + y_grid(i,j)**2
            END DO  
            END DO

            DO i=1,npoint,1
                marker(i) = iter*0.01 + SIN(x_grid(1,1)**2) + i
            END DO
        END SUBROUTINE update_data

END MODULE data


PROGRAM main
    !###########################################################
    ! Using various modules
    USE mod_h5_utility
    USE mpi
    USE data
    USE ISO_FORTRAN_ENV, ONLY: REAL64, INT64, error_unit
    IMPLICIT NONE

    !###########################################################
    ! Declaring variables
    INTEGER :: npx, npy ! number of processors along each dimension
    TYPE(h5_dataset_type) :: griddata, iterdata, fielddata, markerdata  
    INTEGER :: restart, i, j
    INTEGER :: mpierr, mpirank, mpiprocessors, ierror, mp_cart, color, mpi_0y_1Dx
    INTEGER, ALLOCATABLE, DIMENSION(:) ::  mpi_coord

    !###########################################################
    ! Initializing mpi (message passing interface)
    CALL MPI_INIT(mpierr)
    CALL MPI_COMM_RANK(MPI_COMM_WORLD, mpirank, mpierr)
    CALL MPI_COMM_SIZE(MPI_COMM_WORLD, mpiprocessors, mpierr)

    !###########################################################
    ! Assigning local variables
    debug   = .TRUE. ! from MODULE mod_h5_utility
    ndim    = 2      ! from MODULE data
    npx     = 3
    npy     = 4
    iter    = 1      ! from MODULE data
    restart = 0      ! from MODULE data

    !###########################################################
    ! Setting mpi Cartesian structure and a sub-communicator
    IF(.NOT. ALLOCATED(mpi_coord)) ALLOCATE(mpi_coord(ndim))

    ! Getting Cartesian coordinates of the processors. Rank indexing and coord indexing is zero based in mpi. Last dimension is filled first.
    CALL MPI_CART_CREATE(MPI_COMM_WORLD, ndim, (/npy, npx/), (/.FALSE., .FALSE./), .FALSE., mp_cart, ierror) !! CALL mpi_cart_CREATE(MPI_COMM_WORLD, dim, nblocks, periodic_BC, reorder, mp_cart(out), ierror(out))
    CALL MPI_CART_COORDS(mp_cart, mpirank, ndim, mpi_coord, ierror) !! CALL MPI_CART_COORDS(comm, rank, no. of dims, coords, ierr)
    ! for (npy, npx)=(2,3): rank:[0,1,2,3,4,5] -> mpi_coord"[(0,0), (0,1), (0,2), (1,0), (1,1), (1,2)] i.e., highest dimension is filled first.
    WRITE(error_unit,*) 'MPI Cartesian Coord for mpirank = ', mpirank, ' is ', mpi_coord
    CALL MPI_BARRIER(MPI_COMM_WORLD, mpierr)

    ! Creating a sub-communicator
    IF(mpi_coord(1)==0) THEN
        color = 1
    ELSE
        color = MPI_UNDEFINED
    END IF
    CALL MPI_COMM_SPLIT(mp_cart, color, mpi_coord(2), mpi_0y_1Dx, ierror)

    !###########################################################
    ! Defining grid, local to the current rank
    CALL define_grid(npx, npy, mpi_coord)

    !###########################################################
    ! Initializing field variables, local to the current rank
    CALL initialize_data()

    !###########################################################
    ! Setting datatype for int and float for hdf5; This should be called once.
    CALL h5_utility_set_datatype('H5T_STD_I64LE','H5T_NATIVE_DOUBLE')
    ! Setting parallel I/O capability for mod_h5_utility
    h5_utility_mpi = .TRUE. 

    !###########################################################
    ! Creation of datasets; This should be called once before writing the data.
    ! data 0: /grid/x, /grid/y (both share same dataspace)
    CALL griddata%create(global_array_size=(/nxg, nyg/), data_address='field.h5/grid', dataset_names=['x','y'], &
                       & restart=restart, restart_ind=1_INT64,                                                  &
                       & slab_start_ind=(/1+nx*mpi_coord(2), 1+ny*mpi_coord(1)/),                               &
                       & slab_end_ind=(/nx*(1+mpi_coord(2)), ny*(1+mpi_coord(1))/) )
    
    ! data 1: /iter (important for a restart) (DATA SHOULD BE INTEGER TYPE, since it will be read as integer)
    CALL iterdata%create(global_array_size=(/1_INT64/), data_address='field.h5/', dataset_names=['iter'],        &
                       & restart=restart, restart_step_value=INT(iter,INT64), restart_step_data_address='/iter', &
                       data_type='int') 

    ! data 2: /field/velocity, /field/pressure (both share same dataspace)
    CALL fielddata%create(global_array_size=(/nxg,nyg/), data_address='field.h5/field', dataset_names=['velocity','pressure'], &
                        & restart=restart, restart_ind=iterdata%rest_index,                                                    &
                        & slab_start_ind=(/1+nx*mpi_coord(2), 1+ny*mpi_coord(1)/),                                             &
                        & slab_end_ind=(/nx*(1+mpi_coord(2)), ny*(1+mpi_coord(1))/) )

    ! data 3: /marker (is written only by the sub-subcommunicator ranks)
    IF(mpi_0y_1Dx /= MPI_COMM_NULL) THEN ! Can not have same h5 file opened by more than one mpi communicator if they have different number of ranks.
        CALL markerdata%create(global_array_size=(/npointg/), data_address='marker.h5/', dataset_names=['points'], &
                            & restart=restart, restart_ind=iterdata%rest_index,                                    &
                            & slab_start_ind=(/1+npoint*mpi_coord(2)/),                                            &
                            & slab_end_ind=(/npoint*(1+mpi_coord(2))/),                                            &
                            & mpi_comm=mpi_0y_1Dx)
    END IF

    !###########################################################
    ! Writing grid data. This data is written only once.
    ! Preparing append for hdf5 file. This sets the hyperslab, and extends the dataset by n_extend if 
    ! the size of the file dataset is exhausted along iteration dimension.
    CALL griddata%prepare_next_append(n_extend=1)
    ! Writing grid data in each dataset of group grid
    CALL griddata%append(x_grid, data_index=1)
    CALL griddata%append(y_grid, data_index=2)

    !###########################################################
    ! Writing data every iteration
    DO i=1,8
        WRITE(error_unit,*) 'loop: ', i, ' ', mpirank
        WRITE(*,*) 'loop: ', i, ' ', mpirank

        ! Updating field variables
        CALL update_data()

        ! Preparing append for hdf5 file. This sets the hyperslab, and extends the dataset by n_extend if 
        ! the file dataset size is exhausted along iteration dimension.
        CALL iterdata%prepare_next_append(n_extend=5)
        CALL fielddata%prepare_next_append(n_extend=5)
        IF(mpi_0y_1Dx /= MPI_COMM_NULL) CALL markerdata%prepare_next_append(n_extend=7) ! Only executed for 0th y-rank and all x-ranks

        ! Appending data to hdf5 file. HDF5 write operation waits until all ranks in the communicator reach it. 
        CALL iterdata%append((/(INT(iter,INT64))/),  data_index=1)
        CALL fielddata%append(velocity, data_index=1)
        CALL fielddata%append(pressure, data_index=2)
        IF(mpi_0y_1Dx /= MPI_COMM_NULL) CALL markerdata%append(marker, data_index=1) ! Only executed for 0th y-rank and all x-ranks
    
        iter = iter+1
        IF(h5_utility_mpi) CALL MPI_BARRIER(MPI_COMM_WORLD, mpierr)
    END DO

    !###########################################################
    ! Closing all the resources of hdf5. Should be included at the end of the program to avoid memory leak.
    IF(mpi_0y_1Dx /= MPI_COMM_NULL) CALL markerdata%destructor() 
    ! The destructor of an object based on sub-communicator should not be called last due to call h5close_f, 
    ! which is to be called by all the ranks together.
    CALL griddata%destructor()
    CALL iterdata%destructor()
    CALL fielddata%destructor()
    
    !###########################################################
    ! Finalizing mpi
    CALL MPI_FINALIZE(mpierr)

END PROGRAM MAIN



