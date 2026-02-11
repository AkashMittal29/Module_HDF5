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
    IMPLICIT NONE

    REAL(kind=8), ALLOCATABLE, DIMENSION(:,:) :: velocity, pressure
    REAL(kind=8), ALLOCATABLE, DIMENSION(:)   :: marker
    INTEGER(kind=8) :: nx, ny, npoint, restart
    TYPE(h5_dataset_type) :: fielddata, markerdata  

    iter = 0 ! from MODULE data
    nx   = 4
    ny   = 7
    npoint = 5
    debug = .TRUE. ! from MODULE mod_h5_utility

    IF(.NOT. ALLOCATED(velocity)) ALLOCATE(velocity(nx,ny))
    IF(.NOT. ALLOCATED(pressure)) ALLOCATE(pressure(nx,ny))
    IF(.NOT. ALLOCATED(marker))   ALLOCATE(marker(nx))

    CALL fielddata%create(array_size=(/nx,ny/), data_address='field.h5/field',dataset_names=['velocity','pressure'], restart=0)
    CALL markerdata%create(array_size=(/npoint/), data_address='field.h5/marker',dataset_names=['points'], restart=0)

    CALL fielddata%append(velocity, index=1)



    WRITE(*,*) "file:",fielddata%file_name, fielddata%file_id, fielddata%dims_max(3)
    WRITE(*,*) "file:",markerdata%file_name, markerdata%file_id, markerdata%dims_max(2)

    WRITE(*,*) "file:",fielddata%file_name, fielddata%group_id
    WRITE(*,*) "file:",markerdata%file_name, markerdata%group_id

    ! IF(restart\=1) THEN
    !     fielddata%create_dataspace(dims_initial=/nx,ny,0/, dims_max=/nx,ny,HUGE(0_HSIZE_T)/)
    !     markerdata%cteate_dataspace(dims_initial=/npoint, 0/, dims_max=/npoint,HUGE(0_HSIZE_T)/)

    ! END IF





    CALL assign_data(velocity)

    CALL fielddata%destructor()
    CALL markerdata%destructor()


END PROGRAM MAIN


