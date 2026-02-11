
MODULE mod_h5_utility
    USE hdf5 ! The inbuilt module of hdf5
    USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY : error_unit
    IMPLICIT NONE

    PRIVATE
    PUBLIC :: h5_dataset_type
    PUBLIC :: debug
    PUBLIC :: hdf5_initialized

    LOGICAL :: debug
    LOGICAL :: hdf5_initialized = .FALSE.
    
    ! Object with this type can have many datasets under a common file, space, group, plist, and memory space.
    ! Another object will be required if a dataset uses another data space. However, file_id can be copied to
    ! keep the dataset under the same file.
    TYPE h5_dataset_type
        INTEGER(kind=hid_t) :: file_id, space_id, group_id, plist_id, memspace_id ! Common for each dataset
        INTEGER(kind=hid_t),     DIMENSION(:), ALLOCATABLE :: dataset_id ! To have multiple datasets
        INTEGER(kind = HSIZE_T), DIMENSION(:), ALLOCATABLE :: dims_initial, dims_max, dims_chunk, start, count, dims_mem
        CHARACTER(LEN=:), ALLOCATABLE :: file_name
        INTEGER :: error, restart

        CONTAINS
            PROCEDURE, PASS(self), PUBLIC  :: create
            PROCEDURE, PASS(self), PRIVATE :: create_file_id
            PROCEDURE, PASS(self), PUBLIC  :: append
            PROCEDURE, PASS(self), PUBLIC  :: destructor
    END TYPE h5_dataset_type

    TYPE h5_dataset_type_ptr
        TYPE(h5_dataset_type), POINTER :: p
    END TYPE h5_dataset_type_ptr

    TYPE(h5_dataset_type_ptr), POINTER :: opened_file_objects(:) ! pointer to array of type h5_dataset_type_ptr; Being just a pointer, allocatable is invalid.
    INTEGER :: count_files=0


    CONTAINS
        ! Procedures for user defined type: h5_dataset_tye
        SUBROUTINE create(self, array_size, data_address, dataset_names, restart)
            CLASS(h5_dataset_type), INTENT(INOUT) :: self
            INTEGER(kind=HSIZE_T), DIMENSION(:), INTENT(IN)  :: array_size
            CHARACTER(LEN=*), INTENT(IN) :: data_address
            CHARACTER(LEN=*), DIMENSION(:), INTENT(IN) :: dataset_names ! the LEN of each string will be that of the longest string. Therefore use TRIM(ADJUSTL()).
            INTEGER, INTENT(IN) :: restart
            CHARACTER(LEN=:), DIMENSION(:), ALLOCATABLE :: address_parts
            INTEGER :: n_dim, i
            INTEGER(kind=hid_t) :: temp_id

            IF(debug) WRITE(error_unit,*) "mod_h5_write/create/{"

            ! Initializing hdf5 library. This will be executed once in the entire run.
            IF(.NOT. hdf5_initialized) THEN
                CALL h5open_f(self%error) ! Initializing HDF5. Only once is preferred. Multiple times is allowed but need to close also.
                hdf5_initialized = .TRUE.
            END IF

            ! Separating the data address
            CALL separate_string(data_address, address_parts) ! adress_parts(1): file name, subsequent parts are considered as group names. 
            self%file_name = TRIM(ADJUSTL(address_parts(1)))

            ! Creating file id
            CALL self%create_file_id(restart)

            ! Allocating dimensions
            ! dims_initial, dims_max, dims_chunk, start, count, dims_mem
            n_dim = SIZE(array_size) ! no. of dimensions in array
            IF(.NOT. ALLOCATED(self%dims_initial)) ALLOCATE( &
                                & self%dims_initial(1+n_dim), &
                                & self%dims_max(1+n_dim),     &
                                & self%dims_chunk(1+n_dim),   &
                                & self%start(1+n_dim),        &
                                & self%count(1+n_dim),        &
                                & self%dims_mem(n_dim)        &
                                & ) ! +1 for iteration

            ! Assigning dimensions
            self%dims_initial(n_dim+1) = 0
            self%dims_max(n_dim+1)     = HUGE(0_HSIZE_T) ! H5S_UNLIMITED_F will not work.
            WRITE(*,*) " large value: ",self%dims_max(n_dim+1), " ", HUGE(0_HSIZE_T)
            self%dims_chunk(n_dim+1)   = 1 ! chunk's last dim. chunk size is only required during creation of dataset.
            self%count(n_dim+1)        = 1 ! Used for hyperslab. One time step data is written in an iteration.
            DO i=1,n_dim
                self%dims_initial(i) = array_size(i)
                self%dims_max(i)     = array_size(i)
                self%dims_chunk(i)   = array_size(i)
                self%start(i)        = 0 ! Used for hyperslab; 0-based indexing in HDF5. start(n_dim+1) will be assigned later.
                self%count(i)        = array_size(i)
                self%dims_mem(i)     = array_size(i)
            END DO

            ! Creating ids for restart not 1
            restart_if: IF(restart/=1) THEN
                ! Creating space id
                CALL h5screate_simple_f(n_dim+1, self%dims_initial, self%space_id, self%error, self%dims_max) ! Creating dataspace, (rank, dims, space_id, error)

                ! Creating group id if data_address contains group name
                IF(SIZE(address_parts)>1) THEN
                     DO i=2,SIZE(address_parts)
                        IF(i==2) THEN
                            CALL h5gcreate_f(self%file_id, TRIM(ADJUSTL(address_parts(i))), self%group_id, self%error)
                            temp_id = self%group_id
                        ELSE
                            CALL h5gcreate_f(temp_id, TRIM(ADJUSTL(address_parts(i))), self%group_id, self%error)
                            CALL h5gclose_f(temp_id, self%error)
                            temp_id = self%group_id
                        END IF
                    END DO
                END IF

                ! Creating property list for chunking
                CALL h5pcreate_f(H5P_DATASET_CREATE_F, self%plist_id, self%error) ! Creating the property list for creating dataset
                CALL h5pset_chunk_f(self%plist_id, n_dim+1, self%dims_chunk, self%error) ! Setting the property for chunking; plist_id, rank of the data, dimension of the chunk, error

                ! Creating dataset id
                IF(.NOT. ALLOCATED(self%dataset_id)) ALLOCATE(self%dataset_id(SIZE(dataset_names)))
                DO i=1,SIZE(dataset_names)
                    CALL h5dcreate_f(self%group_id, dataset_names(i), H5T_NATIVE_DOUBLE, self%space_id, &
                                   & self%dataset_id(i), self%error, self%plist_id) ! Creating dataset, (group_id/file_id, data_name, type (H5 specific), spce_id, dataset_id, error, property list)
                END DO

                ! Creating memory space
                CALL h5screate_simple_f(n_dim, self%dims_mem, self%memspace_id, self%error)
            
            ELSE ! restart_if
                ! getting dataset id's
                DO i=1,SIZE(dataset_names)
                    CALL h5dopen_f(self%file_id, &
                                 & TRIM(ADJUSTL(data_address(LEN(self%file_name)+2:)))//'/'//dataset_names(i), &
                                 & self%dataset_id(i), self%error)
                END DO

                ! getting dataspace id
                CALL h5dget_space_f(self%dataset_id(1), self%space_id, self%error)

                ! getting group id
                ! is not required.

                ! getting and checking current dimensions
                CALL h5sget_simple_extent_ndims_f(self%space_id, n_dim, self%error)
                IF(n_dim/=SIZE(array_size)+1) THEN
                    WRITE(error_unit,*) "mod_h5_utility/create/ dataset dimesion rank mismatch with the read .h5 file."
                    STOP ! or MPI_ABORT
                END IF
                CALL h5sget_simple_extent_dims_f(self%space_id, self%dims_initial, self%dims_max, self%error) ! gives dims_initial and dims_max corresponding to the space_id
                DO i=1,n_dim
                    IF(self%dims_initial(i)/=array_size(i)) THEN
                        WRITE(error_unit,*) "mod_h5_utility/create/ dataset dimesion mismatch with the read .h5 file."
                        STOP ! or MPI_ABORT
                    END IF
                END DO

            END IF restart_if

            IF(debug) WRITE(error_unit,*) "mod_h5_write/create/}"
        END SUBROUTINE create

        
        SUBROUTINE create_file_id(self, restart)
            CLASS(h5_dataset_type), INTENT(INOUT), TARGET  :: self
            INTEGER, INTENT(IN) :: restart
            INTEGER :: index

            IF(debug) WRITE(error_unit,*) "mod_h5_write/create_file_id/{"

            IF(count_files==0) THEN
                count_files = 1
                ALLOCATE(opened_file_objects(1))
                opened_file_objects(1)%p => self
                IF(restart/=1) THEN
                    CALL h5fcreate_f(self%file_name, H5F_ACC_TRUNC_F, self%file_id, self%error) ! (file_name, File_access_flag, file_id, error)
                    IF(self%error /= 0) WRITE(error_unit,*) "Unable to create file."
                ELSE
                    CALL h5fcreate_f(self%file_name, H5F_ACC_RDWR_F, self%file_id, self%error) ! (file_name, File_access_flag, file_id, error)
                    IF(self%error /= 0) WRITE(error_unit,*) "Unable to open file."
                END IF
            ELSE ! Checking if the file is already opened
                index = is_file_already_opened(self%file_name)
                IF(index/=0) THEN
                    self%file_id = opened_file_objects(index)%p%file_id
                ELSE
                    IF(restart/=1) THEN
                        CALL h5fcreate_f(self%file_name, H5F_ACC_TRUNC_F, self%file_id, self%error) ! (file_name, File_access_flag, file_id, error)
                        IF(self%error /= 0) WRITE(error_unit,*) "Unable to create file."
                    ELSE
                        CALL h5fcreate_f(self%file_name, H5F_ACC_RDWR_F, self%file_id, self%error) ! (file_name, File_access_flag, file_id, error)
                        IF(self%error /= 0) WRITE(error_unit,*) "Unable to open file."
                    END IF
                    count_files = count_files+1
                    CALL extend_opened_file_objects()
                    opened_file_objects(count_files)%p => self
                END IF
            END IF

            IF(debug) WRITE(error_unit,*) "mod_h5_write/create_file_id/}"
        END SUBROUTINE create_file_id


        


        SUBROUTINE destructor(self)
            CLASS(h5_dataset_type), INTENT(INOUT) :: self
            LOGICAL :: is_valid

            IF(debug) WRITE(error_unit,*) "mod_h5_write/destructor/{"

            CALL h5iis_valid_f(self%file_id, is_valid, self%error)
            IF(is_valid) CALL h5fclose_f(self%file_id, self%error) ! closing file

            IF(hdf5_initialized) THEN
                CALL h5close_f(self%error) ! closing the h5 instance. Only once is allowed.
                hdf5_initialized = .FALSE.
            END IF
            
            
            IF(debug) WRITE(error_unit,*) "mod_h5_write/destructor/}"
        END SUBROUTINE destructor


        ! Procedures for module: mod_hf_utility
        SUBROUTINE separate_string(text, out_text)
            CHARACTER(LEN=*), INTENT(IN) :: text
            CHARACTER(LEN=:), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: out_text
            INTEGER :: length, count, i, start, end

            IF(debug) WRITE(error_unit,*) "mod_h5_write/separate_string/{"

            length = LEN_TRIM(text)
            IF(text(length:length)=='/') length = length-1 ! Ignoring trailing '/'.

            ! Counting no. of '/'.
            count = 0
            DO i=1,length
               IF(text(i:i)=='/') count = count+1
            END DO

            ! Allocating out_text with assumed llength 30 characters.
            ALLOCATE(CHARACTER(LEN=30) :: out_text(count+1))

            ! Separating each string
            start = 1
            DO i=1,count
                end = INDEX(text(start:),'/')+start-1
                out_text(i) = text(start:end-1)
                start = end+1
            END DO
            out_text(i) = text(start:length)

            IF(debug) WRITE(error_unit,*) "mod_h5_write/separate_string/}"
        END SUBROUTINE separate_string


        FUNCTION is_file_already_opened(file_name) RESULT(index)
            IMPLICIT NONE
            CHARACTER(LEN=*), INTENT(IN) :: file_name
            INTEGER :: index
            INTEGER :: i

            IF(debug) WRITE(error_unit,*) "mod_h5_write/is_file_already_opened/{"

            index = 0
            DO i=1,count_files
                IF(opened_file_objects(i)%p%file_name == file_name) index = i
            END DO  

            IF(debug) WRITE(error_unit,*) "mod_h5_write/is_file_already_opened/}"
        END FUNCTION is_file_already_opened


        SUBROUTINE extend_opened_file_objects()
            TYPE(h5_dataset_type_ptr), POINTER :: temp(:) ! pointer to array of type h5_dataset_type_ptr
            INTEGER :: n, i

            IF(debug) WRITE(error_unit,*) "mod_h5_write/extend_opened_file_objects/{"

            n = SIZE(opened_file_objects)
            ALLOCATE(temp(n))
            DO i=1,n
                temp(i)%p => opened_file_objects(i)%p
            END DO

            DEALLOCATE(opened_file_objects)
            ALLOCATE(opened_file_objects(n+1))

            DO i=1,n
                opened_file_objects(i)%p => temp(i)%p
            END DO
            opened_file_objects(n+1)%p => NULL()
            DEALLOCATE(temp)

            IF(debug) WRITE(error_unit,*) "mod_h5_write/extend_opened_file_objects/}"
        END SUBROUTINE extend_opened_file_objects



END MODULE mod_h5_utility
