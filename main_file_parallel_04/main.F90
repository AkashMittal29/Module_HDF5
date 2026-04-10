!###########################################################
! EXAMPLE CODE PARALLEL MODE 04
! By Akash Kumar Mittal
! Mar 2026
! Mechanical and Aerospace Engineering
! Florida State University, USA

! In this example, each rank writes different memory hyperslab from a given 
! memory size to a portion in the file space. The mod_h5_utility module ensures
! that the size in memory hyperslab is equal to the size in file hyperslab.
! Here, the data within a defined box is written. 
!###########################################################


! Defining a module data to generate an example of a typical dataset encountered in CFD. 
! It is used in the main program below.
MODULE data
    !###########################################################
    ! Declaring variables
    IMPLICIT NONE
    INTEGER, PARAMETER :: real = selected_real_kind(15, 307)
    INTEGER :: iter = 0, ndim ! Module variables have saved attribute automatically.
    REAL(kind=real), ALLOCATABLE, DIMENSION(:,:) :: x_grid, y_grid, velocity, pressure
    REAL(kind=real), ALLOCATABLE, DIMENSION(:)   :: marker, x_grid_1d, y_grid_1d
    REAL(kind=real)    :: dx, dy
    INTEGER(kind=8) :: nxg, nyg, npointg ! Global dimensions
    INTEGER :: nx, ny, npoint          ! Dimensions local to each rank
    INTEGER(kind=8), ALLOCATABLE, DIMENSION(:) :: ltog_ind_map_x, ltog_ind_map_y ! Local to global index map

    !###########################################################
    ! Defining subroutines of module data
    CONTAINS

        !###########################################################
        ! Assigns variables, allocates arrays, and defines grid
        SUBROUTINE define_grid(npx, npy, mpi_coord)
            INTEGER, INTENT(IN) :: npx, npy ! number of processors along each dimension
            INTEGER, DIMENSION(:), INTENT(IN) :: mpi_coord
            REAL(KIND=real) :: x_start, y_start
            INTEGER :: i, j

            ! Assigning variables
            dx      = 0.2_real
            dy      = 0.1_real
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
             IF(.NOT. ALLOCATED(ltog_ind_map_x)) ALLOCATE(ltog_ind_map_x(nx))
            IF(.NOT. ALLOCATED(ltog_ind_map_y))  ALLOCATE(ltog_ind_map_y(ny))

            ! Assigning grid
            x_start = mpi_coord(2)*nx*dx ! mpi_coord starts from 0.
            x_grid_1d = (/(x_start+i*dx, i=0,nx-1)/) ! Implied do loop
            y_start = mpi_coord(1)*ny*dy
            y_grid_1d = (/(y_start+i*dy, i=0,ny-1)/)
            DO i=1,nx
            DO j=1,ny
                x_grid(i,j) = x_grid_1d(i)
                ltog_ind_map_x(i) = i+mpi_coord(2)*nx
            END DO
            END DO
            DO i=1,ny
                y_grid(:,i) = y_grid_1d(i)
                ltog_ind_map_y(i) = i+mpi_coord(1)*ny
            END DO
        END SUBROUTINE define_grid

        !###########################################################
        ! Initializes data
        SUBROUTINE initialize_data()
            velocity = -1.0_real
            pressure = -10.0_real
            marker   = 0.0_real
        END SUBROUTINE initialize_data

        !###########################################################
        ! Updates the arrays
        SUBROUTINE update_data()
            INTEGER :: i,j

            ! Assigning data by some operations
            DO j=1,ny,1
            DO i=1,nx,1
                velocity(i,j) = iter*0.01_real + x_grid(i,j)**2 + y_grid(i,j) 
                pressure(i,j) = iter*0.01_real + x_grid(i,j)**3 + y_grid(i,j)**2
            END DO  
            END DO

            DO i=1,npoint,1
                marker(i) = iter*0.01_real + SIN(x_grid(1,1)**2) + i
            END DO
        END SUBROUTINE update_data

        SUBROUTINE get_indices_for_box(x_grid_1d_, y_grid_1d_, x_box, y_box, is_box_present, &
                                     & x_ind, y_ind, nx_local_box, ny_local_box)
            REAL(kind=real), DIMENSION(:),INTENT(IN) :: x_grid_1d_, y_grid_1d_, x_box, y_box
            INTEGER, DIMENSION(2), INTENT(OUT) :: x_ind, y_ind
            LOGICAL, INTENT(OUT) :: is_box_present
            INTEGER, INTENT(OUT)  :: nx_local_box, ny_local_box
            LOGICAL :: x_start_found, x_end_found, y_start_found, y_end_found
            INTEGER :: i

            ! Initializing
            is_box_present = .FALSE. 
            nx_local_box = 0
            ny_local_box = 0

            x_start_found  = .FALSE.  
            DO i=1,nx
                IF(x_grid_1d_(i)>=x_box(1) .AND. x_grid_1d_(i)<=x_box(2)) THEN
                    IF(.NOT. x_start_found) THEN
                        x_start_found = .TRUE.
                        x_ind(1) = i
                    ELSE
                        x_ind(2) = i
                    END IF
                END IF
            END DO

            y_start_found = .FALSE. 
            DO i=1,ny
                IF(y_grid_1d_(i)>=y_box(1) .AND. y_grid_1d_(i)<=y_box(2)) THEN
                    IF(.NOT. y_start_found) THEN
                        y_start_found = .TRUE.
                        y_ind(1) = i
                    ELSE
                        y_ind(2) = i
                    END IF
                END IF
            END DO

            IF(x_start_found .AND. y_start_found) THEN
                is_box_present = .TRUE.
                nx_local_box = x_ind(2)-x_ind(1)+1
                ny_local_box = y_ind(2)-y_ind(1)+1
            ELSE
                is_box_present = .FALSE.
            END IF
        END SUBROUTINE get_indices_for_box


        SUBROUTINE print_array_2d(array)
            REAL(kind=real), DIMENSION(:,:), INTENT(IN) :: array
            INTEGER i1, i2

            DO i1=LBOUND(array,1), UBOUND(array,1)
                DO i2=LBOUND(array,2), UBOUND(array,2)
                    WRITE(*,'(A,F8.2,A)',ADVANCE='NO') ' ', array(i1,i2), ' '
                END DO
                WRITE(*,*)
            END DO
        END SUBROUTINE print_array_2d

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
    ! INTEGER, PARAMETER :: real = selected_real_kind(15, 307) ! Here it is used form module data.
    REAL(kind=real), PARAMETER :: eps = 1e-15
    INTEGER :: npx, npy ! number of processors along each dimension
    TYPE(h5_dataset_type) :: fullgriddata, griddata, iterdata, fielddata, markerdata  
    INTEGER :: restart, i, j
    INTEGER :: mpierr, mpirank, mpiprocessors, ierror, mp_cart, color, mpi_sub
    INTEGER, ALLOCATABLE, DIMENSION(:) ::  mpi_coord
    REAL(kind=real) :: x_box(2) = (/0.4_real-eps, 2.2_real+eps/), y_box(2) = (/0.3_real-eps, 0.9_real+eps/) ! the box coordinates within which data is to be written.
    LOGICAL :: is_box_present
    INTEGER :: x_ind(2), y_ind(2), nxl_box, nyl_box, nxg_box, nyg_box ! local indices for box, and local & global number of nodes in the box.
    INTEGER :: xg_ind_min, yg_ind_min ! min global index values from all ranks.
    INTEGER :: file_slab_start(2), file_slab_end(2), mem_start(2) ! file hyperslab indices for each rank
    INTEGER :: npx_box, npy_box, npx_min, npy_min, npx_max, npy_max

    !###########################################################
    ! Initializing mpi (message passing interface)
    CALL MPI_INIT(mpierr)
    CALL MPI_COMM_RANK(MPI_COMM_WORLD, mpirank, mpierr)
    CALL MPI_COMM_SIZE(MPI_COMM_WORLD, mpiprocessors, mpierr)

    !###########################################################
    ! Assigning local variables
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

    CALL MPI_BARRIER(MPI_COMM_WORLD, mpierr)

    !###########################################################
    ! Defining grid, local to the current rank
    CALL define_grid(npx, npy, mpi_coord)

    ! Creating a sub-communicator by checking if the current rank contains the box
    ! Box is represented by 'o' in the following illustration.
    !  ---------------------------
    ! |. . . . | . . . . | . . . .|
    ! |. . o o | o o o . | . . . .|
    ! |. . o o | o o o . | . . . .|
    !  ---------------------------
    ! |. . o o | o o o . | . . . .|
    ! |. . . . | . . . . | . . . .|    ^ y 
    ! |. . . . | . . . . | . . . .|    |
    !  ---------------------------      --> x
    CALL get_indices_for_box(x_grid_1d, y_grid_1d, x_box, y_box, is_box_present, x_ind, y_ind, nxl_box, nyl_box)
    IF(is_box_present) THEN 
        color = 1
    ELSE
        color = MPI_UNDEFINED
    END IF
    CALL MPI_COMM_SPLIT(mp_cart, color, 0, mpi_sub, ierror) ! key given here is 0, mpi will take care of the repeated key.

    ! Getting total number of nodes in the box in each direction, and global indices of the box (used for file hyperslab)
    IF(mpi_sub /= MPI_COMM_NULL) THEN
        npx_box = mpi_coord(2) ! x-coordinate of the current rank
        npy_box = mpi_coord(1) ! y-coordinate of the current rank
        CALL MPI_ALLREDUCE(npx_box, npx_min, 1, MPI_INTEGER, MPI_MIN, mpi_sub, ierror)
        CALL MPI_ALLREDUCE(npx_box, npx_max, 1, MPI_INTEGER, MPI_MAX, mpi_sub, ierror)
        CALL MPI_ALLREDUCE(npy_box, npy_min, 1, MPI_INTEGER, MPI_MIN, mpi_sub, ierror)
        CALL MPI_ALLREDUCE(npy_box, npy_max, 1, MPI_INTEGER, MPI_MAX, mpi_sub, ierror)
        npx_box = npx_max-npx_min+1 ! Number of processors in x-dir for box
        npy_box = npy_max-npy_min+1 ! Number of processors in y-dir for box

        CALL MPI_ALLREDUCE(nxl_box, nxg_box, 1, MPI_INTEGER, MPI_SUM, mpi_sub, ierror)
        CALL MPI_ALLREDUCE(nyl_box, nyg_box, 1, MPI_INTEGER, MPI_SUM, mpi_sub, ierror)
        nxg_box = nxg_box/npy_box
        nyg_box = nyg_box/npx_box
        xg_ind_min = ltog_ind_map_x(x_ind(1)) ! min x_index of the box local to each rank
        yg_ind_min = ltog_ind_map_y(y_ind(1)) ! min y_index of the box local to each rank

        ! Getting global min indices of the entire box 
        CALL MPI_ALLREDUCE(MPI_IN_PLACE, xg_ind_min, 1, MPI_INTEGER, MPI_MIN, mpi_sub, ierror)
        CALL MPI_ALLREDUCE(MPI_IN_PLACE, yg_ind_min, 1, MPI_INTEGER, MPI_MIN, mpi_sub, ierror)

        ! eg. start and end x-index of the file hyperslab for the current rank:
        ! start: ltog_ind_map_x(x_ind(1))-xg_ind_min+1
        ! end  : ltog_ind_map_x(x_ind(2))-xg_ind_min+1
    END IF

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
    CALL fullgriddata%create( global_array_size = (/nxg, nyg/),                                &
                            & data_address      = 'field_full.h5/grid',                        &
                            & dataset_names     = ['x','y'],                                   &
                            & restart           = restart,                                     &
                            & restart_ind       = 1_INT64,                                     &
                            & slab_start_ind    = (/1+nx*mpi_coord(2), 1+ny*mpi_coord(1)/),    &
                            & slab_end_ind      = (/nx*(1+mpi_coord(2)), ny*(1+mpi_coord(1))/) )

    ! In the following eg., each rank has unique file-hyperslab and memory-hyperslab.
    IF(mpi_sub /= MPI_COMM_NULL) THEN
        ! data 1: /iter (important for a restart) (DATA SHOULD BE INTEGER TYPE, since it will be read as integer)
        CALL iterdata%create( global_array_size         = (/1_INT64/),     &
                            & data_address              = 'field.h5/',     &
                            & dataset_names             = ['iter'],        &
                            & restart                   = restart,         &
                            & restart_step_value        = INT(iter,INT64), &
                            & restart_step_data_address = '/iter',         &
                            & data_type                 = 'int',           &
                            & mpi_comm                  = mpi_sub )

        ! Setting file hyperslab start and end indices for each rank
        file_slab_start = (/ltog_ind_map_x(x_ind(1))-xg_ind_min+1, ltog_ind_map_y(y_ind(1))-yg_ind_min+1/)
        file_slab_end   = (/ltog_ind_map_x(x_ind(2))-xg_ind_min+1, ltog_ind_map_y(y_ind(2))-yg_ind_min+1/)

        ! data 2: /grid/x, /grid/y (both share same dataspace)
        CALL griddata%create( global_array_size = INT((/nxg_box, nyg_box/), INT64), &
                            & data_address      = 'field.h5/subgrid',               &
                            & dataset_names     = ['x','y'],                        &
                            & restart           = restart,                          & 
                            & restart_ind       = 1_INT64,                          &
                            & slab_start_ind    = file_slab_start,                  &
                            & slab_end_ind      = file_slab_end,                    &
                            & mem_size          = SHAPE(x_grid),                    &
                            & mem_start         = (/x_ind(1), y_ind(1)/),           &
                            & mpi_comm          = mpi_sub )

        ! data 3: /field/velocity, /field/pressure (both share same dataspace)
        CALL fielddata%create(global_array_size = INT((/nxg_box, nyg_box/), INT64), &
                            & data_address      = 'field.h5/field',                 &
                            & dataset_names     = ['velocity','pressure'],          &
                            & restart           = restart,                          &
                            & restart_ind       = iterdata%rest_index,              &
                            & slab_start_ind    = file_slab_start,                  &
                            & slab_end_ind      = file_slab_end,                    &
                            & mem_size          = SHAPE(velocity),                  &
                            & mem_start         = (/x_ind(1), y_ind(1)/),           &
                            & mpi_comm          = mpi_sub )
    END IF
    
    !###########################################################
    ! Writing grid data. This data is written only once.
    ! Preparing append for hdf5 file. This sets the hyperslab, and extends the dataset by n_extend if 
    ! the size of the file dataset is exhausted along iteration dimension.
    CALL fullgriddata%prepare_next_append(n_extend=1)
    ! Writing grid data in each dataset of group grid
    CALL fullgriddata%append(x_grid, data_index=1)
    CALL fullgriddata%append(y_grid, data_index=2)

    IF(mpi_sub /= MPI_COMM_NULL) THEN
        CALL griddata%prepare_next_append(n_extend=1)
        ! CALL print_array_2d(y_grid)
        CALL griddata%append(x_grid, data_index=1)
        CALL griddata%append(y_grid, data_index=2)
    END IF

    !###########################################################
    ! Writing data every iteration
    DO i=1,6
        WRITE(error_unit,*) 'loop: ', i, ' ', mpirank
        WRITE(*,*) 'loop: ', i, ' ', mpirank

        ! Updating field variables
        CALL update_data()

        IF(mpi_sub /= MPI_COMM_NULL) THEN
            ! Preparing append for hdf5 file. This sets the hyperslab, and extends the dataset by n_extend if 
            ! the file dataset size is exhausted along iteration dimension.
            CALL iterdata%prepare_next_append(n_extend=5)
            CALL fielddata%prepare_next_append(n_extend=5) ! Both field variables have identical size and indexing. So, one 
            
            ! Appending data to hdf5 file. HDF5 write operation waits until all ranks in the communicator reach it. 
            CALL iterdata%append((/(INT(iter,INT64))/),  data_index=1)
            CALL fielddata%append(velocity, data_index=1) ! As per memory hyperslab, only a region of the field variable is written into fiile.
            CALL fielddata%append(pressure, data_index=2)

            !IF(mpirank==4) WRITE(*,'(A,F20.16)') 'x(1,1) in rank 4: ', x_grid(1,1)
        END IF
       
        iter = iter+1
        IF(h5_utility_mpi) CALL MPI_BARRIER(MPI_COMM_WORLD, mpierr)
    END DO

    !###########################################################
    ! Closing all the resources of hdf5. Should be included at the end of the program to avoid memory leak.
    IF(mpi_sub /= MPI_COMM_NULL) THEN
        CALL griddata%destructor() 
        CALL iterdata%destructor()
        CALL fielddata%destructor()
    END IF
    ! The destructor of an object based on sub-communicator should not be called last due to call h5close_f, 
    ! which is to be called by all the ranks together.
    CALL fullgriddata%destructor()
    
    
    !###########################################################
    ! Finalizing mpi
    CALL MPI_FINALIZE(mpierr)

END PROGRAM MAIN


