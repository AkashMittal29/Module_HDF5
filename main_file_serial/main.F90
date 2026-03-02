!###########################################################
! EXAMPLE CODE SERIAL MODE
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
    INTEGER(KIND=8) :: nx, ny, npoint

    !###########################################################
    ! Defining subroutines of module data
    CONTAINS

        !###########################################################
        ! Assigns variables, allocates arrays, and defines grid
        SUBROUTINE define_grid(npx, npy)
            INTEGER, INTENT(IN) :: npx, npy ! number of processors along each dimension. Used here only to match the parallel counterpart.
            REAL(KIND=8) :: x_start, y_start
            INTEGER :: i

            ! Assigning variables
            dx      = 0.2
            dy      = 0.1
            nx      = 5
            ny      = 6
            npoint  = 9
            nx     = npx*nx ! This is done to match with the parallel counterpart.
            ny     = npy*ny
            npoint = npx*npoint

            ! Allocating
            IF(.NOT. ALLOCATED(x_grid_1d)) ALLOCATE(x_grid_1d(nx))
            IF(.NOT. ALLOCATED(y_grid_1d)) ALLOCATE(y_grid_1d(ny))
            IF(.NOT. ALLOCATED(x_grid))    ALLOCATE(x_grid(nx,ny))
            IF(.NOT. ALLOCATED(y_grid))    ALLOCATE(y_grid(nx,ny))
            IF(.NOT. ALLOCATED(velocity))  ALLOCATE(velocity(nx,ny))
            IF(.NOT. ALLOCATED(pressure))  ALLOCATE(pressure(nx,ny))
            IF(.NOT. ALLOCATED(marker))    ALLOCATE(marker(npoint))

            ! Assigning grid
            x_start = 0.0
            x_grid_1d = (/(x_start+i*dx, i=0,nx-1)/) ! Implied do loop
            y_start = 0.0
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
    USE data
    USE ISO_FORTRAN_ENV, ONLY: REAL64, INT64, error_unit
    IMPLICIT NONE

    !###########################################################
    ! Declaring variables
    INTEGER :: npx, npy ! number of processors along each dimension (in serial mode, these are only used to calculate nx and ny)
    TYPE(h5_dataset_type) :: griddata, iterdata, fielddata, markerdata  
    INTEGER :: restart, i, j

    !###########################################################
    ! Assigning local variables
    debug   = .TRUE. ! from MODULE mod_h5_utility
    ndim    = 2      ! from MODULE data
    npx     = 3
    npy     = 4
    iter    = 1      ! from MODULE data
    restart = 0      ! from MODULE data

    !###########################################################
    ! Defining grid, local to the current rank
    CALL define_grid(npx, npy)

    !###########################################################
    ! Initializing field variables, local to the current rank
    CALL initialize_data()

    !###########################################################
    ! Setting datatype for int and float for hdf5; This should be called once.
    CALL h5_utility_set_datatype('H5T_STD_I64LE','H5T_NATIVE_DOUBLE')
    ! Setting parallel I/O capability for mod_h5_utility
    h5_utility_mpi = .FALSE. 

    !###########################################################
    ! Creation of datasets; This should be called once before writing the data.
    ! data 0: /grid/x, /grid/y (both share same dataspace)
    CALL griddata%create(global_array_size=(/nx, ny/), data_address='field.h5/grid', dataset_names=['x','y'], &
                       & restart=restart, restart_ind=1_INT64)
    
    ! data 1: /iter (important for a restart) (DATA SHOULD BE INTEGER TYPE, since it will be read as integer)
    CALL iterdata%create(global_array_size=(/1_INT64/), data_address='field.h5/', dataset_names=['iter'],        &
                       & restart=restart, restart_step_value=INT(iter,INT64), restart_step_data_address='/iter', &
                       data_type='int') 

    ! data 2: /field/velocity, /field/pressure (both share same dataspace)
    CALL fielddata%create(global_array_size=(/nx, ny/), data_address='field.h5/field', dataset_names=['velocity','pressure'], &
                        & restart=restart, restart_ind=iterdata%rest_index)

    ! data 3: /marker
    CALL markerdata%create(global_array_size=(/npoint/), data_address='marker.h5/', dataset_names=['points'], &
                        & restart=restart, restart_ind=iterdata%rest_index)
                        
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
        WRITE(error_unit,*) 'loop: ', i
        WRITE(*,*) 'loop: ', i

        ! Updating field variables
        CALL update_data()

        ! Preparing append for hdf5 file. This sets the hyperslab, and extends the dataset by n_extend if 
        ! the file dataset size is exhausted along iteration dimension.
        CALL iterdata%prepare_next_append(n_extend=5)
        CALL fielddata%prepare_next_append(n_extend=5)
        CALL markerdata%prepare_next_append(n_extend=7)

        ! Appending data to hdf5 file. HDF5 write operation waits until all ranks in the communicator reach it. 
        CALL iterdata%append((/(INT(iter,INT64))/),  data_index=1)
        CALL fielddata%append(velocity, data_index=1)
        CALL fielddata%append(pressure, data_index=2)
        CALL markerdata%append(marker, data_index=1)
    
        iter = iter+1
    END DO

    !###########################################################
    ! Closing all the resources of hdf5. Should be included at the end of the program to avoid memory leak.
    CALL markerdata%destructor() 
    CALL griddata%destructor()
    CALL iterdata%destructor()
    CALL fielddata%destructor()

END PROGRAM MAIN


