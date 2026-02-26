MODULE data
    IMPLICIT NONE
    INTEGER :: iter = 0 ! Module variables have saved attribute automatically.

    CONTAINS
        SUBROUTINE assign_data(velocity)
            IMPLICIT NONE
            REAL(kind=8), INTENT(INOUT), DIMENSION(:,:) :: velocity
            INTEGER :: i,j

            ! Assigning data
            DO j=1,SIZE(velocity,dim=2),1
            DO i=1,SIZE(velocity,dim=1),1
                velocity(i,j) = iter + i+(j-1)*(SIZE(velocity,dim=1))
            END DO  
            END DO
        END SUBROUTINE assign_data

END MODULE data


PROGRAM main
    USE mod_h5_utility
    USE data
    USE iso_fortran_env, ONLY: REAL64, INT64
    IMPLICIT NONE

    REAL(kind=8), ALLOCATABLE, DIMENSION(:,:) :: velocity, pressure
    REAL(kind=8), ALLOCATABLE, DIMENSION(:)   :: marker
    INTEGER(kind=8) :: nx, ny, npoint, i
    TYPE(h5_dataset_type) :: iterdata, fielddata, markerdata  
    INTEGER :: restart

    iter = 15 ! from MODULE data
    nx   = 4
    ny   = 7
    npoint = 5
    debug = .TRUE. ! from MODULE mod_h5_utility
    restart = 1

    IF(.NOT. ALLOCATED(velocity)) ALLOCATE(velocity(nx,ny))
    IF(.NOT. ALLOCATED(pressure)) ALLOCATE(pressure(nx,ny))
    IF(.NOT. ALLOCATED(marker))   ALLOCATE(marker(npoint))

    ! Initializing field variables
    CALL assign_data(velocity)
    CALL assign_data(pressure)
    marker = 0.00

    ! Setting datatype for int and float for hdf5; This should be called once.
    CALL h5_utility_set_datatype('H5T_STD_I64LE','H5T_NATIVE_DOUBLE')

    ! Creation of datasets; This should be called once before writing the data.
    ! data 1: /iter (important for a restart) (should be integer since it will be read as integer)
    CALL iterdata%create(array_size=(/1_INT64/), data_address='field.h5/', dataset_names=['iter'], &
                       & restart=restart, restart_step_value=INT(iter,INT64), restart_step_data_address='/iter', &
                       data_type='int') 
    ! data 2: /field/velocity, /field/pressure (both share same dataspace)
    CALL fielddata%create(array_size=(/nx,ny/), data_address='field.h5/field', dataset_names=['velocity','pressure'], &
                        & restart=restart, restart_ind=iterdata%restart_index)
    ! data 3: /marker
    CALL markerdata%create(array_size=(/npoint/), data_address='field.h5/marker', dataset_names=['points'], &
                         & restart=restart, restart_ind=iterdata%restart_index)

    ! Writing data
    DO i=1,7
        ! Updating field variables
        velocity = velocity+2.0
        pressure = pressure+10.0
        marker = marker+1.0

        ! Preparing append for hdf5 file. This sets the hyperslab, and extends the dataset by n_extend if 
        ! the file dataset size is exhausted along iteration dimension.
        CALL iterdata%prepare_next_append(n_extend=5)
        CALL fielddata%prepare_next_append(n_extend=5)
        CALL markerdata%prepare_next_append(n_extend=5)
        
        ! Appending data to hdf5 file
        CALL iterdata%append((/(INT(iter,INT64))/),  data_index=1)
        CALL fielddata%append(velocity, data_index=1)
        CALL fielddata%append(pressure, data_index=2)
        CALL markerdata%append(marker,  data_index=1)

        WRITE(*,*) 'loop: ', i
        iter = iter+1
    END DO

    ! Closing all the resources of hdf5. Should be included only at the end of the program to avoid memory leak.
    CALL iterdata%destructor()
    CALL fielddata%destructor()
    CALL markerdata%destructor()


END PROGRAM MAIN


