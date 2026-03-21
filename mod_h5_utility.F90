!###########################################################
! MODULE mod_h5_utility
! A Fortran wrapper for writing hdf5 files.
! By Akash Kumar Mittal
! Mar 2026
! Mechanical and Aerospace Engineering
! Florida State University, USA
!###########################################################


MODULE mod_h5_utility
    !###########################################################
    ! Using various modules
    USE hdf5 ! The module of hdf5
    USE mpi 
    USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY : error_unit
    USE iso_fortran_env, ONLY: REAL64, INT64
    IMPLICIT NONE

    !###########################################################
    ! Setting default access
    PRIVATE

    !###########################################################
    ! Public interface
    PUBLIC :: h5_dataset_type
    PUBLIC :: h5_utility_set_datatype
    PUBLIC :: hdf5_initialized
    PUBLIC :: flag_datatype_set
    PUBLIC :: h5_utility_mpi
    PUBLIC :: n_h5_dataset_type_objects

    !###########################################################
    ! Declaring module variables
    CHARACTER(LEN=20) :: h5_utility_real, h5_utility_int
    LOGICAL :: flag_datatype_set = .FALSE.
    LOGICAL :: hdf5_initialized = .FALSE.
    logical :: h5_utility_mpi = .FALSE.
    INTEGER :: count_files = 0, n_h5_dataset_type_objects=0
    
    !###########################################################
    ! Defining derived datatype
    ! Object with this type can have many datasets under a common: file, file space & hyperslab, memory space & hyperslab, group, and plist.
    ! Another object will be required if a dataset uses another file data space or memory space (input variable's size and its hyperslab). 
    ! However, file_id is automatically copied if it is created under the same file to avoid opening multiple file handles.
    ! With two different communicators (having different set of ranks), same file can not be opened.
    TYPE h5_dataset_type
        !###########################################################
        ! Declaring variables
        INTEGER(kind=hid_t) :: file_id=-1, space_id=-1, group_id=-1, plist_id=-1, memspace_id=-1, transfer_plist_id=-1 ! Common for each dataset
        INTEGER(kind=hid_t),     DIMENSION(:), ALLOCATABLE :: dataset_id ! To have multiple datasets
        INTEGER(kind = HSIZE_T), DIMENSION(:), ALLOCATABLE :: dims_current, dims_total, dims_max, dims_chunk, &
                                                            & start, count, dims_mem, mem_start, mem_count
        CHARACTER(LEN=:),                      ALLOCATABLE :: file_name
        INTEGER           :: n_dim_full, error, n_extend, extension_track
        INTEGER(KIND=8)   :: rest_index
        LOGICAL           :: extension_exhausted, has_group = .FALSE.
        CHARACTER(LEN=20) :: data_type
        LOGICAL           :: is_created = .FALSE., is_destroyed = .FALSE. ! These ensure that create and destructor subroutines are executed only once.
        INTEGER           :: mpi_comm ! MPI communicator

        !###########################################################
        ! Declaring procedures of the derived datatype
        CONTAINS
            PROCEDURE, PASS(self), PUBLIC  :: create
            PROCEDURE, PASS(self), PRIVATE :: create_file_id
            PROCEDURE, PASS(self), PRIVATE :: get_rest_index
            PROCEDURE, PASS(self), PRIVATE :: create_datasets
            PROCEDURE, PASS(self), PUBLIC  :: prepare_next_append
            PROCEDURE, PASS(self), PUBLIC  :: append_real64_1d
            PROCEDURE, PASS(self), PUBLIC  :: append_real32_1d
            PROCEDURE, PASS(self), PUBLIC  :: append_real64_2d
            PROCEDURE, PASS(self), PUBLIC  :: append_real32_2d
            PROCEDURE, PASS(self), PUBLIC  :: append_int_1d
            PROCEDURE, PASS(self), PUBLIC  :: destructor
            
            GENERIC :: append => append_real64_1d, append_real32_1d, &
                               & append_real64_2d, append_real32_2d, &
                               & append_int_1d ! Procedure overloading
    END TYPE h5_dataset_type

    !###########################################################
    ! Declaring more module variables
    TYPE h5_dataset_type_ptr
        TYPE(h5_dataset_type), POINTER :: p
        INTEGER :: opened_as ! 0 -> with restart = 0, 1 -> with restart = 1: for checking purpose while opening file.
    END TYPE h5_dataset_type_ptr

    TYPE(h5_dataset_type_ptr), POINTER :: opened_file_objects(:) ! pointer to array of type h5_dataset_type_ptr; Being just a pointer, allocatable is invalid.
    

    !###########################################################
    ! Defining derived datatype's procedures, module subroutines, and module functions
    CONTAINS

        !###########################################################
        ! Procedures for derived type: h5_dataset_type

        !###########################################################
        SUBROUTINE create(self, global_array_size, data_address, dataset_names, &
                        & restart, restart_step_value, restart_step_data_address, restart_ind, &
                        & data_type, slab_start_ind, slab_end_ind, mpi_comm, mem_size, mem_start)
            CLASS(h5_dataset_type), INTENT(INOUT), TARGET        :: self
            INTEGER(kind=HSIZE_T),  INTENT(IN), DIMENSION(:) :: global_array_size ! kind=HSIZE_T (INT64) to accommodate large integers for global_array_size
            CHARACTER(LEN=*),       INTENT(IN)               :: data_address
            CHARACTER(LEN=*),       INTENT(IN), DIMENSION(:) :: dataset_names ! the LEN of each string will be that of the longest string. Therefore use TRIM(ADJUSTL()).
            INTEGER,                INTENT(IN)               :: restart
            INTEGER(KIND=8),        INTENT(IN), OPTIONAL     :: restart_step_value ! must be int64
            CHARACTER(LEN=*),       INTENT(IN), OPTIONAL     :: restart_step_data_address
            CHARACTER(LEN=*),       INTENT(IN), OPTIONAL     :: data_type ! By default, real (h5_utility_real_id)
            INTEGER(KIND=8),        INTENT(IN), OPTIONAL     :: restart_ind
            INTEGER, DIMENSION(:),  INTENT(IN), OPTIONAL     :: slab_start_ind, slab_end_ind ! global start index is 1. However, HDF5 uses zero-based indexing, which is taken care in this module.
            INTEGER,                INTENT(IN), OPTIONAL     :: mpi_comm ! MPI communicator   
            INTEGER, DIMENSION(:),  INTENT(IN), OPTIONAL     :: mem_size, mem_start ! To define memory space for the input variable in h5 write command if its size is different from the size of data to be 
                                                                                    ! written into the file hyperslab. Accordingly, memory hyperslab's start index is to be set.
            
            CHARACTER(LEN=:), DIMENSION(:), ALLOCATABLE :: address_parts
            INTEGER,          DIMENSION(:), ALLOCATABLE :: slabstart, slabend
            INTEGER(KIND=8),  DIMENSION(:), ALLOCATABLE :: mem_stride, mem_block
            INTEGER             :: n_dim, i
            INTEGER(kind=hid_t) :: temp_id
            INTEGER(KIND=8)     :: rest_index

#ifdef DEBUGhdf5
WRITE(error_unit,*) "mod_h5_utility/create/{"
#endif
            ! Incrementing object counter
            IF(self%is_created .EQV. .FALSE.) THEN
                n_h5_dataset_type_objects = n_h5_dataset_type_objects+1
                self%is_created = .TRUE.
            ELSE
                WRITE(error_unit, *) 'mod_h5_utility/create/: object is calling subroutine create more than once.'
                CALL abort_program()
            END IF

            ! Initializing hdf5 library. This will be executed once in the entire run.
            IF(.NOT. hdf5_initialized) THEN
                CALL h5open_f(self%error) ! Initializing HDF5. Only once is preferred. Multiple times is allowed but need to close the same times.
                hdf5_initialized = .TRUE.
            END IF

            ! Checking if datatypes are set
            IF (.NOT. flag_datatype_set) THEN 
                WRITE(error_unit,*) 'mod_h5_utility/create/: datatypes are not defined. Call subroutine h5_utility_set_datatype.'
                CALL abort_program()
            END IF

            ! Assigning datatype
            IF(PRESENT(data_type)) THEN
                IF(data_type=='int') THEN
                    self%data_type = h5_utility_int
                ELSE IF(data_type=='real') THEN
                    self%data_type = h5_utility_real
                ELSE 
                    WRITE(error_unit,*) "mod_h5_utility/create/: wrong data_type provided: available are- 'int', 'real'."
                    CALL abort_program()
                END IF
            ELSE ! Setting default type as real
                self%data_type = h5_utility_real
            END IF

            ! Assigning MPI communicator
            IF(h5_utility_mpi) THEN
                IF (PRESENT(mpi_comm)) THEN
                    self%mpi_comm = mpi_comm
                ELSE
                    self%mpi_comm = MPI_COMM_WORLD
                END IF
            END IF

            ! Initializing extension tracking attributes
            self%extension_exhausted = .TRUE.
            self%n_extend = 0
            self%extension_track = 0

            ! Separating the data address
            CALL separate_string(data_address, address_parts) ! adress_parts(1): file name, subsequent parts are considered as group names. 
            self%file_name = TRIM(ADJUSTL(address_parts(1)))

            ! Creating file id
            CALL self%create_file_id(restart)

            ! Allocating dimensions
            ! dims_current, dims_max, dims_chunk, start, count, dims_mem
            n_dim = SIZE(global_array_size)  ! no. of dimensions in array
            self%n_dim_full = n_dim+1 ! including the extendable dimension
            IF(.NOT. ALLOCATED(self%dims_current)) ALLOCATE(  &
                                & self%dims_current(1+n_dim), &
                                & self%dims_total(1+n_dim),   &
                                & self%dims_max(1+n_dim),     &
                                & self%dims_chunk(1+n_dim),   &
                                & self%start(1+n_dim),        &
                                & self%count(1+n_dim),        &
                                & self%dims_mem(n_dim),       &
                                & self%mem_start(n_dim),      &
                                & self%mem_count(n_dim),      &
                                & mem_stride(n_dim),          &
                                & mem_block(n_dim),           &
                                & slabstart(n_dim),           &
                                & slabend(n_dim)              & 
                                & ) ! +1 for iteration

            ! Assigning dimensions
            self%dims_current(n_dim+1) = 0
            self%dims_max(n_dim+1)     = HUGE(0_HSIZE_T) ! H5S_UNLIMITED_F will not work.
            ! WRITE(*,*) " large value: ",self%dims_max(n_dim+1), " ", HUGE(0_HSIZE_T) !***
            self%dims_chunk(n_dim+1)   = 1 ! chunk's last dim. chunk size is only required during creation of dataset.
            self%count(n_dim+1)        = 1 ! Used for hyperslab. One time step data is written in an iteration.
            IF( PRESENT(slab_start_ind) .AND. PRESENT(slab_end_ind) ) THEN
                slabstart = slab_start_ind
                slabend   = slab_end_ind
            ELSE
                slabstart = 1         ! all elements in slabstart are assumed 1
                slabend   = global_array_size 
            END IF
            
            DO i=1,n_dim
                ! (Note: for each processor, every local object in a processor is disctinct from the object created in other processor.)
                self%dims_current(i)  = global_array_size(i)
                self%dims_max(i)      = global_array_size(i)
                self%start(i)         = slabstart(i)-1       ! Used for file hyperslab; HDF5 uses zero-based indexing. start(n_dim+1) will be assigned later.
                self%count(i)         = slabend(i)-slabstart(i)+1
                self%dims_chunk(i)    = global_array_size(i) ! chunk should have the dimension of the global array size, however the extendible dimension may vary. 
                IF(PRESENT(mem_size) .AND. PRESENT(mem_start)) THEN
                    self%dims_mem(i)  = mem_size(i)          ! memory space is used to write dataset with a memory hyperslab.
                    self%mem_start(i) = mem_start(i)-1       ! HDF5 uses zero-based indexing.
                    mem_stride(i)     = 1
                    mem_block(i)      = 1
                    self%mem_count(i) = self%count(i)        ! Should be equal to self%count(i) element-wise.
                ELSE
                    self%dims_mem(i)  = self%count(i)        ! If variable passed to the append subroutine has same size as that of the file hyperslab's size (self%count).
                    self%mem_start(i) = 0                    ! Default value is 0. HDF5 uses zero-based indexing.
                    mem_stride(i)     = 1                    ! Default memory stride (for contgouous memory, such as for the entire array in the memory).
                    mem_block(i)      = 1                    ! Default is 1. Block -> How many jumps for the next element in the memory.
                    self%mem_count(i) = self%count(i)        ! Should be equal to self%count(i) element-wise.
                END IF
                
            END DO
            self%dims_total = self%dims_current
            

            ! Creating ids 
            restart_if: IF(restart/=1) THEN ! for restart not 1
                ! Creating space id for global dataset
                CALL h5screate_simple_f(n_dim+1, self%dims_current, self%space_id, self%error, self%dims_max) ! Creating dataspace, (rank, dims, space_id, error)

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
                    self%has_group = .TRUE.
                ELSE
                    self%group_id = self%file_id
                    self%has_group = .FALSE.
                END IF

                ! Creating property list for chunking
                CALL h5pcreate_f(H5P_DATASET_CREATE_F, self%plist_id, self%error)        ! Creating the property list for creating dataset
                CALL h5pset_chunk_f(self%plist_id, n_dim+1, self%dims_chunk, self%error) ! Setting the property for chunking; plist_id, rank of the data, dimension of the chunk, error

                ! Creating dataset id
                IF(.NOT. ALLOCATED(self%dataset_id)) ALLOCATE(self%dataset_id(SIZE(dataset_names)))
                CALL self%create_datasets(dataset_names)

            ELSE ! restart_if
                ! getting dataset id's
                IF(.NOT. ALLOCATED(self%dataset_id)) ALLOCATE(self%dataset_id(SIZE(dataset_names)))
                DO i=1,SIZE(dataset_names)
                    CALL h5dopen_f(self%file_id, &
                                 & '/'//TRIM(ADJUSTL(data_address(LEN(self%file_name)+2:)))//'/'//dataset_names(i), &
                                 & self%dataset_id(i), self%error)
                END DO

                ! getting dataspace id
                CALL h5dget_space_f(self%dataset_id(1), self%space_id, self%error)

                ! getting group id
                ! is not required.
                
                ! getting and checking total dimensions
                CALL h5sget_simple_extent_ndims_f(self%space_id, n_dim, self%error)
                n_dim = n_dim-1
                IF(n_dim/=SIZE(global_array_size)) THEN
                    WRITE(error_unit,*) "mod_h5_utility/create/ dataset dimension rank mismatch with the read .h5 file."
                    CALL abort_program()
                END IF
                CALL h5sget_simple_extent_dims_f(self%space_id, self%dims_total, self%dims_max, self%error) ! gives dims_total (current) and dims_max corresponding to the space_id
                DO i=1,n_dim
                    IF(self%dims_total(i)/=global_array_size(i)) THEN
                        WRITE(error_unit,*) "mod_h5_utility/create/ dataset dimension mismatch with the read .h5 file."
                        CALL abort_program()
                    END IF
                END DO

                ! Setting where to start writing data in the file space for a restart
                IF(PRESENT(restart_step_value) .AND. PRESENT(restart_step_data_address)) THEN
                    rest_index = self%get_rest_index(restart_step_value, restart_step_data_address)
                ELSE IF(PRESENT(restart_ind)) THEN
                    rest_index = restart_ind ! restart_ind minimum value can be 1.
                ELSE
                    WRITE(error_unit,*) 'mod_h5_utility/create/: with restart=1, restart_step_value and ' // &
                                      & 'restart_step_data_address (automatic search), or restart_ind should' // &
                                      & 'be provided to know from where to start writing data.'
                    CALL abort_program()
                END IF
                IF(rest_index<1) THEN
                    WRITE(error_unit,*) 'mod_h5_utility/create/: restart_ind can not be less than 1.'
                    CALL abort_program()
                END IF

                self%rest_index = rest_index ! Required, so that other bjects of type h5_dataset_type can copy it.
                self%n_extend = self%dims_total(self%n_dim_full)-rest_index+1_INT64 ! Including the rest_index, therefore +1.
                self%extension_exhausted = .FALSE.
                self%dims_current(self%n_dim_full) = rest_index-1 ! From rest_index the data is to be written.
                self%extension_track = 0

            END IF restart_if

            ! Creating memory space
            CALL h5screate_simple_f(n_dim, self%dims_mem, self%memspace_id, self%error) ! Is used for creating memory hyperslab and writing data

            ! Selecting the memory hyperslab. 
            ! This is more efficient than passing a sliced data which is often not contiguous in memory (results in creation of temporary variable in fortran). 
            ! Hence, just pass the entire data array (only address of first element is passed), and h5 can write only the hyperslab from memory to the file hyperslab.
            CALL h5sselect_hyperslab_f(self%memspace_id, H5S_SELECT_SET_F, self%mem_start, self%mem_count, self%error, &
                                     & stride=mem_stride, block=mem_block)

            ! Creating dataset transfer mode property list
            CALL h5pcreate_f(H5P_DATASET_XFER_F, self%transfer_plist_id, self%error) ! if dxpl_mpio (next line) is not set: in serial: no effect, in parallel: default is independent mode.
            IF(h5_utility_mpi) CALL h5pset_dxpl_mpio_f(self%transfer_plist_id, H5FD_MPIO_COLLECTIVE_F, self%error) ! for mpi, setting collective data transfer: i/o is coordinated among ranks -> more efficient than independent mode.
            !IF(h5_utility_mpi) CALL h5pset_dxpl_mpio_f(self%transfer_plist_id, H5FD_MPIO_INDEPENDENT_F, self%error)

            ! Deallocating local arrays
            IF(ALLOCATED(slabstart))  DEALLOCATE(slabstart)
            IF(ALLOCATED(slabend))    DEALLOCATE(slabend)
            IF(ALLOCATED(mem_stride)) DEALLOCATE(mem_stride)
            IF(ALLOCATED(mem_block))  DEALLOCATE(mem_block)

#ifdef DEBUGhdf5
WRITE(error_unit,*) "mod_h5_utility/create/}"
#endif
        END SUBROUTINE create

        
        !###########################################################
        SUBROUTINE create_file_id(self, restart)
            CLASS(h5_dataset_type), INTENT(INOUT), TARGET  :: self
            INTEGER, INTENT(IN) :: restart
            INTEGER :: index

#ifdef DEBUGhdf5
WRITE(error_unit,*) "mod_h5_utility/create_file_id/{"
#endif

            ! Setting property list for accessing file with parallel I/O
            IF(h5_utility_mpi) THEN
                CALL h5pcreate_f(H5P_FILE_ACCESS_F, self%plist_id, self%error)
                CALL h5pset_fapl_mpio_f(self%plist_id, self%mpi_comm, MPI_INFO_NULL, self%error)
            END IF
            
            IF(count_files==0) THEN
                count_files = 1
                ALLOCATE(opened_file_objects(1))
                opened_file_objects(1)%p => self
                IF(restart/=1) THEN
                    IF(h5_utility_mpi) THEN
                        CALL h5fcreate_f(self%file_name, H5F_ACC_TRUNC_F, self%file_id, self%error, access_prp=self%plist_id) ! (file_name, File_access_flag, file_id, error)
                    ELSE
                        CALL h5fcreate_f(self%file_name, H5F_ACC_TRUNC_F, self%file_id, self%error) ! (file_name, File_access_flag, file_id, error)
                    END IF
                    IF(self%error /= 0) THEN
                        WRITE(error_unit,*) "Unable to create file."
                        CALL abort_program()
                    END IF
                    opened_file_objects(count_files)%opened_as = 0 ! 0 -> with restart = 0
                ELSE
                    IF(h5_utility_mpi) THEN
                        CALL h5fopen_f(self%file_name, H5F_ACC_RDWR_F, self%file_id, self%error, access_prp=self%plist_id) ! (file_name, File_access_flag, file_id, error)
                    ELSE
                        CALL h5fopen_f(self%file_name, H5F_ACC_RDWR_F, self%file_id, self%error) ! (file_name, File_access_flag, file_id, error)
                    END IF
                    IF(self%error /= 0) THEN 
                        WRITE(error_unit,*) "Unable to open file."
                        CALL abort_program()
                    END IF
                    opened_file_objects(count_files)%opened_as = 1 ! 1 -> with restart = 1
                END IF
            ELSE ! Checking if the file is already opened
                index = is_file_already_opened(self%file_name)
                IF(index/=0) THEN
                    self%file_id = opened_file_objects(index)%p%file_id
                    IF((restart==1 .AND. opened_file_objects(index)%opened_as/=1) .OR. &
                     & (restart==0 .AND. opened_file_objects(index)%opened_as/=0)) THEN
                        WRITE(error_unit,*) "restart value and .h5 file access specifier mismatch. " // &
                                          & "Check all objects accessing same file have same restart value."
                        CALL abort_program()
                    END IF
        
                    ! Checking MPI communicator
                    IF(h5_utility_mpi .AND. self%mpi_comm/=opened_file_objects(index)%p%mpi_comm) THEN
                        WRITE(error_unit,*) 'mod_h5_utility/create_file_id/: Same file is opened with different '// &
                                          & 'MPI communicators, which is unsafe if the ranks are not identical in both. Aborting.'
                        CALL abort_program()
                    END IF
                ELSE
                    IF(restart/=1) THEN
                        IF(h5_utility_mpi) THEN
                            CALL h5fcreate_f(self%file_name, H5F_ACC_TRUNC_F, self%file_id, self%error, access_prp=self%plist_id) ! (file_name, File_access_flag, file_id, error)
                        ELSE
                            CALL h5fcreate_f(self%file_name, H5F_ACC_TRUNC_F, self%file_id, self%error) ! (file_name, File_access_flag, file_id, error)
                        END IF
                        IF(self%error /= 0) THEN
                            WRITE(error_unit,*) "Unable to create file."
                            CALL abort_program()
                        END IF
                        count_files = count_files+1
                        CALL extend_opened_file_objects()
                        opened_file_objects(count_files)%p => self
                        opened_file_objects(count_files)%opened_as = 0
                    ELSE
                        IF(h5_utility_mpi) THEN
                            CALL h5fopen_f(self%file_name, H5F_ACC_RDWR_F, self%file_id, self%error, access_prp=self%plist_id) ! (file_name, File_access_flag, file_id, error)
                        ELSE
                            CALL h5fopen_f(self%file_name, H5F_ACC_RDWR_F, self%file_id, self%error) ! (file_name, File_access_flag, file_id, error)
                        END IF
                        IF(self%error /= 0) THEN
                            WRITE(error_unit,*) "Unable to open file."
                            CALL abort_program()
                        END IF
                        count_files = count_files+1
                        CALL extend_opened_file_objects()
                        opened_file_objects(count_files)%p => self
                        opened_file_objects(count_files)%opened_as = 1
                    END IF
                    
                END IF
            END IF

            ! Closing plist_id used for accessing file in parallel I/O
            IF(h5_utility_mpi) THEN
                CALL h5pclose_f(self%plist_id, self%error)
                self%plist_id = -1
            END IF

#ifdef DEBUGhdf5
WRITE(error_unit,*) "mod_h5_utility/create_file_id/}"
#endif
        END SUBROUTINE create_file_id


        !###########################################################
        SUBROUTINE create_datasets(self, dataset_names)
            CLASS(h5_dataset_type), INTENT(INOUT) :: self
            CHARACTER(LEN=*), DIMENSION(:), INTENT(IN) :: dataset_names
            INTEGER :: i

#ifdef DEBUGhdf5
WRITE(error_unit,*) "mod_h5_utility/create_datasets/{"
#endif

            DO i=1,SIZE(self%dataset_id)
                ! Creating dataset, (group_id/file_id, data_name, type (H5 specific), spce_id, dataset_id, error, property list)
                ! NOTE: though here one spaceid is used to define size of the dataset, internally 
                ! each dataset will create their own copy of dataspce.
                SELECT CASE (TRIM(ADJUSTL(self%data_type)))
                    CASE('H5T_NATIVE_INTEGER') 
                        CALL h5dcreate_f(self%group_id, TRIM(ADJUSTL(dataset_names(i))), H5T_NATIVE_INTEGER, self%space_id, &
                                    & self%dataset_id(i), self%error, self%plist_id) 
                    CASE('H5T_STD_I16LE')      
                        CALL h5dcreate_f(self%group_id, TRIM(ADJUSTL(dataset_names(i))), H5T_STD_I16LE, self%space_id, &
                                    & self%dataset_id(i), self%error, self%plist_id) 
                    CASE('H5T_STD_I32LE')      
                        CALL h5dcreate_f(self%group_id, TRIM(ADJUSTL(dataset_names(i))), H5T_STD_I32LE, self%space_id, &
                                    & self%dataset_id(i), self%error, self%plist_id) 
                    CASE('H5T_STD_I64LE')      
                        CALL h5dcreate_f(self%group_id, TRIM(ADJUSTL(dataset_names(i))), H5T_STD_I64LE, self%space_id, &
                                    & self%dataset_id(i), self%error, self%plist_id) 
                    CASE('H5T_NATIVE_REAL')   
                        CALL h5dcreate_f(self%group_id, TRIM(ADJUSTL(dataset_names(i))), H5T_NATIVE_REAL, self%space_id, &
                                    & self%dataset_id(i), self%error, self%plist_id) 
                    CASE('H5T_NATIVE_DOUBLE')  
                        CALL h5dcreate_f(self%group_id, TRIM(ADJUSTL(dataset_names(i))), H5T_NATIVE_DOUBLE, self%space_id, &
                                    & self%dataset_id(i), self%error, self%plist_id) 
                    CASE('H5T_IEEE_F32LE')    
                        CALL h5dcreate_f(self%group_id, TRIM(ADJUSTL(dataset_names(i))), H5T_IEEE_F32LE, self%space_id, &
                                    & self%dataset_id(i), self%error, self%plist_id) 
                    CASE('H5T_IEEE_F64LE')     
                        CALL h5dcreate_f(self%group_id, TRIM(ADJUSTL(dataset_names(i))), H5T_IEEE_F64LE, self%space_id, &
                                    & self%dataset_id(i), self%error, self%plist_id) 
                    CASE DEFAULT
                        WRITE(error_unit,*) 'mod_h5_utility/create_datasets/: Unavailable datatype provided.'// &
                                        & 'Check mod_h5_utility/h5_utility_set_datatype call.'
                        CALL abort_program()
                END SELECT
            END DO

#ifdef DEBUGhdf5
WRITE(error_unit,*) "mod_h5_utility/create_datasets/}"
#endif
        END SUBROUTINE create_datasets


        !###########################################################
        FUNCTION get_rest_index(self, restart_step_value, restart_step_data_address) RESULT(index)
            CLASS(h5_dataset_type), INTENT(INOUT) :: self
            INTEGER(KIND=8),        INTENT(IN)    :: restart_step_value
            CHARACTER(LEN=*),       INTENT(IN)    :: restart_step_data_address
            INTEGER(KIND=8)                :: index
            INTEGER(KIND=hid_t)            :: data_id, space_id, dtype_id
            INTEGER(HSIZE_T), DIMENSION(2) :: dims_total, dims_max
            ! INTEGER :: n_dims
            INTEGER(KIND=8),  DIMENSION(:), ALLOCATABLE :: step_values
            INTEGER(HSIZE_T) :: size_bytes
            INTEGER          :: class_id, i
            ! restart_step_values should have (1,nt) dimension in the file with INTEGER of size 8 bytes (INT64).

#ifdef DEBUGhdf5
WRITE(error_unit,*) "mod_h5_utility/get_rest_index/}"
#endif

            ! Initialilizing index to -1
            index = -1

            ! Getting dataset id
            CALL h5dopen_f(self%file_id, TRIM(ADJUSTL(restart_step_data_address)), data_id, self%error)
            
            ! Getting space id
            CALL h5dget_space_f(data_id, space_id, self%error)

            ! getting current and total (max) dimensions
            ! CALL h5sget_simple_extent_ndims_f(space_id, n_dims, self%error)
            CALL h5sget_simple_extent_dims_f(space_id, dims_total, dims_max, self%error) ! gives dims_total and dims_max corresponding to the space_id
            
            ! Checking dataset type
            CALL h5dget_type_f(data_id, dtype_id, self%error)
            CALL h5tget_class_f(dtype_id, class_id, self%error)
            CALL h5tget_size_f(dtype_id, size_bytes, self%error)
            IF(.NOT. ( class_id==H5T_INTEGER_F .AND. size_bytes==8 )) THEN
                WRITE(error_unit,*) 'mod_h5_utility/create/get_rest_index/: The restart_step_values should ' // &
                                  & 'be written with type INTEGER of size 8 bytes (INT64). ' // &
                                  & 'Check the type of the variable that writes restart_step_values to .h5 file.'
                CALL abort_program()
            END IF
            
            ! Reading data
            ALLOCATE(step_values(dims_total(2)))
            CALL h5dread_f(data_id, H5T_STD_I64LE, step_values, dims_total, self%error)
            !                           ^ memory data type (not file storage type)
            ! From file, data is saved into memory, and then assigned to fortran variable.
            ! Type conversion happenss if the file storage, memory types, and Fortran variable types are compatible types.
            ! H5T_NATIVE_INTEGER is not the safest choice since machine can be 32-bit.

            ! Searching required value
            DO i=1,SIZE(step_values)
                IF(step_values(i)==restart_step_value) THEN
                    index = i
                    EXIT
                END IF
            END DO
            IF(index==-1) THEN
                WRITE(error_unit,*) 'mod_h5_utility/create/get_rest_index/: unable to find ' // &
                                  & 'the restart_step_value in the provided dataset address.'
                CALL abort_program()
            END IF
            DEALLOCATE(step_values)

            ! Closing resources
            CALL h5tclose_f(dtype_id, self%error)
            CALL h5sclose_f(space_id, self%error)
            CALL h5dclose_f(data_id, self%error)
            
#ifdef DEBUGhdf5
WRITE(error_unit,*) "mod_h5_utility/get_rest_index/}"
#endif
        END FUNCTION get_rest_index


        !###########################################################
        SUBROUTINE prepare_next_append(self, n_extend)
            ! If no extension is required for a data where overwriting is desired every step, 
            ! call this subroutine once after calling create with restart = 0.
            ! If data is extendable, then call this subroutine before appending every step. 
            CLASS(h5_dataset_type), INTENT(INOUT) :: self
            INTEGER, INTENT(IN) :: n_extend ! No. of extension elements in extendable dimension.
            ! n_extend is assigned to self%n_extend only if extension_exhausted is true.
            INTEGER :: i, mpierror

#ifdef DEBUGhdf5
WRITE(error_unit,*) "mod_h5_utility/prepare_next_append/{"
#endif

            ASSOCIATE(n_dim_full => self%n_dim_full, extension_exhausted => self%extension_exhausted, &
                    & extension_track => self%extension_track)
                
                IF(extension_exhausted) THEN
                    ! Until extension_exhausted is .TRUE., the input n_extend will not be used.
                    ! Increasing the dimension in time or iter for each dataset by given n_extend
                    self%n_extend = n_extend
                    self%dims_total(n_dim_full) = self%dims_total(n_dim_full) + self%n_extend
                    DO i=1,SIZE(self%dataset_id)
                        IF(h5_utility_mpi) CALL MPI_BARRIER(self%mpi_comm, mpierror)
                        CALL h5dset_extent_f(self%dataset_id(i), self%dims_total, self%error)
                    END DO
                    
                    IF(h5_utility_mpi) CALL MPI_BARRIER(self%mpi_comm, mpierror)
                    extension_exhausted = .FALSE.
                    extension_track = 0

                    ! Creating updated space object id; 
                    CALL h5sclose_f(self%space_id, self%error) ! closing space object
                    CALL h5screate_simple_f(n_dim_full, self%dims_total, self%space_id, self%error) 
                    ! or get space id from the current dataset: CALL h5dget_space_f(self%dataset_id(1), self%space_id, self%error)
                    ! This space object has extension for eg. (nx,ny,nz,nt) (i.e., same as the current dataset); dataset extent is (nx,ny,nz,nt)
                    ! NOTE: Since all datasets have same structure, here one space object is created for writing hyperslab into each dataset.
                END IF

                extension_track = extension_track+1

                ! Extending current dimension for writing data
                self%dims_current(n_dim_full) = self%dims_current(n_dim_full)+1 ! extending the last (time/iter) dimension

                ! setting the hyperslab parameters (for extendable dimension, other dimensions for hyperslab are defined in subroutine create)
                self%start(n_dim_full) = self%dims_current(n_dim_full)-1 ! HDF5 uses zero-based indexing
                
                ! selecting the file hyperslab (each processor rank will call this with their locaal self%start and self%count)
                CALL h5sselect_hyperslab_f(self%space_id, H5S_SELECT_SET_F, self%start, self%count, self%error) ! from the given data space, it selects a given region.
                ! Here, extent of hyperslab (based on count) is for eg. (nx,ny,nz,1) (i.e., a window).

                ! Memory hyperslab is already set in subroutine create.

                IF(extension_track==self%n_extend) extension_exhausted = .TRUE.
            END ASSOCIATE

#ifdef DEBUGhdf5
WRITE(error_unit,*) "mod_h5_utility/prepare_next_append/}"
#endif
        END SUBROUTINE prepare_next_append


        !###########################################################
        SUBROUTINE append_real64_1d(self, array, data_index)
            CLASS(h5_dataset_type), INTENT(INOUT) :: self
            REAL(kind=8), DIMENSION(:), INTENT(IN) :: array
            INTEGER :: data_index

#ifdef DEBUGhdf5
WRITE(error_unit,*) "mod_h5_utility/append_real64_1d/{"
#endif

            CALL h5dwrite_f(self%dataset_id(data_index), H5T_NATIVE_DOUBLE, array, self%dims_mem, &
                          & self%error, mem_space_id=self%memspace_id, file_space_id=self%space_id, xfer_prp=self%transfer_plist_id)
            ! For eg.: 
            !   current dataset extent: (nx,ny,nz,nt)
            !   space extent (file space extent): (nx,ny,nz,nt) : should be same as the current dataset extent
            !                        with hyperslab: start:(0,0,0,nt), count:(nx,ny,nz,1)
            !   memspace extent of the array: (nx,ny,nz) 

#ifdef DEBUGhdf5
WRITE(error_unit,*) "mod_h5_utility/append_real64_1d/}"
#endif
        END SUBROUTINE append_real64_1d


        !###########################################################
        SUBROUTINE append_real32_1d(self, array, data_index)
            CLASS(h5_dataset_type), INTENT(INOUT) :: self
            REAL(kind=4), DIMENSION(:), INTENT(IN) :: array
            INTEGER :: data_index

#ifdef DEBUGhdf5
WRITE(error_unit,*) "mod_h5_utility/append_real32_1d/{"
#endif

            CALL h5dwrite_f(self%dataset_id(data_index), H5T_NATIVE_REAL, array, self%dims_mem, &
                          & self%error, mem_space_id=self%memspace_id, file_space_id=self%space_id, xfer_prp=self%transfer_plist_id)

#ifdef DEBUGhdf5
WRITE(error_unit,*) "mod_h5_utility/append_real32_1d/}"
#endif
        END SUBROUTINE append_real32_1d


        !###########################################################
        SUBROUTINE append_real64_2d(self, array, data_index)
            CLASS(h5_dataset_type), INTENT(INOUT) :: self
            REAL(kind=8), DIMENSION(:,:), INTENT(IN) :: array
            INTEGER :: data_index

#ifdef DEBUGhdf5
WRITE(error_unit,*) "mod_h5_utility/append_real64_2d/{"
#endif

            CALL h5dwrite_f(self%dataset_id(data_index), H5T_NATIVE_DOUBLE, array, self%dims_mem, &
                         & self%error, mem_space_id=self%memspace_id, file_space_id=self%space_id, xfer_prp=self%transfer_plist_id)

#ifdef DEBUGhdf5
WRITE(error_unit,*) "mod_h5_utility/append_real64_2d/}"
#endif
        END SUBROUTINE append_real64_2d


        !###########################################################
        SUBROUTINE append_real32_2d(self, array, data_index)
            CLASS(h5_dataset_type), INTENT(INOUT) :: self
            REAL(kind=4), DIMENSION(:,:), INTENT(IN) :: array
            INTEGER :: data_index

#ifdef DEBUGhdf5
WRITE(error_unit,*) "mod_h5_utility/append_real32_2d/{"
#endif

            CALL h5dwrite_f(self%dataset_id(data_index), H5T_NATIVE_REAL, array, self%dims_mem, &
                          & self%error, mem_space_id=self%memspace_id, file_space_id=self%space_id, xfer_prp=self%transfer_plist_id)

#ifdef DEBUGhdf5
WRITE(error_unit,*) "mod_h5_utility/append_real32_2d/}"
#endif
        END SUBROUTINE append_real32_2d


        !###########################################################
        SUBROUTINE append_real64_3d(self, array, data_index)
            CLASS(h5_dataset_type), INTENT(INOUT) :: self
            REAL(kind=8), DIMENSION(:,:,:), INTENT(IN) :: array
            INTEGER :: data_index

#ifdef DEBUGhdf5
WRITE(error_unit,*) "mod_h5_utility/append_real64_3d/{"
#endif

            CALL h5dwrite_f(self%dataset_id(data_index), H5T_NATIVE_DOUBLE, array, self%dims_mem, &
                         & self%error, mem_space_id=self%memspace_id, file_space_id=self%space_id, xfer_prp=self%transfer_plist_id)

#ifdef DEBUGhdf5
WRITE(error_unit,*) "mod_h5_utility/append_real64_3d/}"
#endif
        END SUBROUTINE append_real64_3d


        !###########################################################
        SUBROUTINE append_real32_3d(self, array, data_index)
            CLASS(h5_dataset_type), INTENT(INOUT) :: self
            REAL(kind=4), DIMENSION(:,:,:), INTENT(IN) :: array
            INTEGER :: data_index

#ifdef DEBUGhdf5
WRITE(error_unit,*) "mod_h5_utility/append_real32_3d/{"
#endif

            CALL h5dwrite_f(self%dataset_id(data_index), H5T_NATIVE_REAL, array, self%dims_mem, &
                          & self%error, mem_space_id=self%memspace_id, file_space_id=self%space_id, xfer_prp=self%transfer_plist_id)

#ifdef DEBUGhdf5
WRITE(error_unit,*) "mod_h5_utility/append_real32_3d/}"
#endif
        END SUBROUTINE append_real32_3d


        !###########################################################
        SUBROUTINE append_int_1d(self, array, data_index)
            CLASS(h5_dataset_type), INTENT(INOUT) :: self
            INTEGER(kind=8), DIMENSION(:), INTENT(IN) :: array
            INTEGER :: data_index

#ifdef DEBUGhdf5
WRITE(error_unit,*) "mod_h5_utility/append_int_1d/{"
#endif

            CALL h5dwrite_f(self%dataset_id(data_index), H5T_STD_I64LE, array, self%dims_mem, &
                          & self%error, mem_space_id=self%memspace_id, file_space_id=self%space_id, xfer_prp=self%transfer_plist_id)

#ifdef DEBUGhdf5
WRITE(error_unit,*) "mod_h5_utility/append_int_1d/}"
#endif
        END SUBROUTINE append_int_1d


        !###########################################################
        SUBROUTINE destructor(self)
            CLASS(h5_dataset_type), INTENT(INOUT) :: self
            LOGICAL :: is_valid
            INTEGER :: i, mpierr
            INTEGER(KIND=8) :: n_objects_open

#ifdef DEBUGhdf5
WRITE(error_unit,*) "mod_h5_utility/destructor/{"
#endif

            ! decrementing object counter
            IF(self%is_destroyed .EQV. .FALSE.) THEN
                n_h5_dataset_type_objects = n_h5_dataset_type_objects-1
                self%is_destroyed = .TRUE.
            ELSE
                WRITE(error_unit, *) 'mod_h5_utility/destructor/: object is calling subroutine destructor more than once.'
                CALL abort_program()
            END IF
            
            ! Closing ids. The sequence of closing resources is important. 
            ! First child sources should be closed, then the fie can be closed. Otherwise, program will hang before termination.
            IF(h5_utility_mpi) CALL MPI_BARRIER(self%mpi_comm, mpierr)
            IF(self%space_id/=-1) THEN 
                CALL h5sclose_f(self%space_id, self%error) ! closing space_id
                self%space_id = -1
            END IF
            
            IF(h5_utility_mpi) CALL MPI_BARRIER(self%mpi_comm, mpierr)
            IF(self%plist_id/=-1) THEN 
                CALL h5pclose_f(self%plist_id, self%error) ! closing plist_id
                self%plist_id = -1
            END IF

            IF(h5_utility_mpi) CALL MPI_BARRIER(self%mpi_comm, mpierr)
            IF(self%transfer_plist_id/=-1) THEN 
                CALL h5pclose_f(self%transfer_plist_id, self%error) ! closing transfer_plist_id
                self%transfer_plist_id = -1
            END IF
            
            IF(h5_utility_mpi) CALL MPI_BARRIER(self%mpi_comm, mpierr)
            IF(self%memspace_id/=-1) THEN 
                CALL h5sclose_f(self%memspace_id, self%error) ! closing memspace_id
                self%memspace_id = -1
            END IF

            IF(ALLOCATED(self%dataset_id)) THEN
                DO i=1,SIZE(self%dataset_id,1)
                    IF(h5_utility_mpi) CALL MPI_BARRIER(self%mpi_comm, mpierr)
                    CALL h5iis_valid_f(self%dataset_id(i), is_valid, self%error)
                    IF(is_valid) CALL h5dclose_f(self%dataset_id(i), self%error) ! closing dataset_id  
                END DO
            END IF

            IF(h5_utility_mpi) CALL MPI_BARRIER(self%mpi_comm, mpierr)
            IF(self%group_id/=-1 .AND. self%has_group) THEN 
                CALL h5gclose_f(self%group_id, self%error) ! closing group_id
                self%group_id = -1
            END IF

            IF(h5_utility_mpi) CALL MPI_BARRIER(self%mpi_comm, mpierr)
            IF(self%file_id/=-1) THEN 
                ! Getting number of opened child objects of the file_id. If any child object is still open, 
                ! then skip closing it since other h5_dataset_type's object might be shareing the same file_id,
                ! and hence closing it will give error.
                CALL h5fget_obj_count_f(self%file_id, H5F_OBJ_ALL_F, n_objects_open, self%error)
                IF(n_objects_open==1) THEN ! will happen only if all the child objects of the file are closed. No. of objects = 1 is the file itself.
                    CALL h5fclose_f(self%file_id, self%error) ! closing file_id
                END IF
                self%file_id = -1
            END IF

            ! Deallocating arrays
            IF(ALLOCATED(self%dataset_id))      DEALLOCATE(self%dataset_id)
            IF(ALLOCATED(self%dims_current))    DEALLOCATE(self%dims_current)
            IF(ALLOCATED(self%dims_total))      DEALLOCATE(self%dims_total)
            IF(ALLOCATED(self%dims_max))        DEALLOCATE(self%dims_max)
            IF(ALLOCATED(self%dims_chunk))      DEALLOCATE(self%dims_chunk)
            IF(ALLOCATED(self%start))           DEALLOCATE(self%start)
            IF(ALLOCATED(self%count))           DEALLOCATE(self%count)
            IF(ALLOCATED(self%dims_mem))        DEALLOCATE(self%dims_mem)
            IF(ALLOCATED(self%mem_start))       DEALLOCATE(self%mem_start)
            IF(ALLOCATED(self%mem_count))       DEALLOCATE(self%mem_count)
            IF(ALLOCATED(self%file_name))       DEALLOCATE(self%file_name)
            IF(ASSOCIATED(opened_file_objects)) DEALLOCATE(opened_file_objects)
            
            ! Closing the hdf5 instance
            IF(hdf5_initialized .AND. n_h5_dataset_type_objects==0) THEN
                IF(h5_utility_mpi) CALL MPI_BARRIER(self%mpi_comm, mpierr)
                CALL h5close_f(self%error) ! closing the hdf5 instance. Only once is allowed.
                hdf5_initialized = .FALSE.
                WRITE(error_unit, *) 'Successful final close of HDF5.'
            END IF
            
#ifdef DEBUGhdf5
WRITE(error_unit,*) "mod_h5_utility/destructor/}"
#endif
        END SUBROUTINE destructor


        !###########################################################
        ! Subroutines and functions for module: mod_hf_utility

        !###########################################################
        SUBROUTINE h5_utility_set_datatype(text_int, text_real)
            CHARACTER(LEN=*), INTENT(IN) :: text_int, text_real

#ifdef DEBUGhdf5
WRITE(error_unit,*) "mod_h5_utility/h5_utility_set_datatype{"
#endif

            SELECT CASE (text_int)
                CASE('H5T_NATIVE_INTEGER') 
                    h5_utility_int = 'H5T_NATIVE_INTEGER' ! Integer size is dependent on the platform corresponds to C native types.
                CASE('H5T_STD_I16LE')      
                    h5_utility_int = 'H5T_STD_I16LE'
                CASE('H5T_STD_I32LE')      
                    h5_utility_int = 'H5T_STD_I32LE'
                CASE('H5T_STD_I64LE')      
                    h5_utility_int = 'H5T_STD_I64LE'
                    ! CALL h5tcopy_f(H5T_STD_I64LE, h5_utility_int_id, error)
                    ! WRITE(*,*) 'DEBUG set_datatype:::: ', error, ' ', H5T_NATIVE_INTEGER, ' ', H5T_STD_I64LE, &
                    !          & ' ', KIND(H5T_STD_I64LE)
                CASE DEFAULT
                    WRITE(error_unit,*) 'mod_h5_utility/h5_utility_set_datatype/: Unavailable type provided.'
                    CALL abort_program()
            END SELECT

            SELECT CASE (text_real)
                CASE('H5T_NATIVE_REAL')    
                    h5_utility_real = 'H5T_NATIVE_REAL' ! Float size is dependent on the platform corresponds to C native types.
                CASE('H5T_NATIVE_DOUBLE')  
                    h5_utility_real = 'H5T_NATIVE_DOUBLE'
                CASE('H5T_IEEE_F32LE')    
                     h5_utility_real = 'H5T_IEEE_F32LE'
                CASE('H5T_IEEE_F64LE')     
                    h5_utility_real = 'H5T_IEEE_F64LE'
                CASE DEFAULT
                    WRITE(error_unit,*) 'mod_h5_utility/h5_utility_set_datatype/: Unavailable type provided.'
                    CALL abort_program()
            END SELECT

            flag_datatype_set = .TRUE.
            ! The options can be extended referring the following link for FORTRAN90 API datatypes: 
            ! https://davis.lbl.gov/Manuals/HDF5-1.6.1/PredefDTypes.html

#ifdef DEBUGhdf5
WRITE(error_unit,*) "mod_h5_utility/h5_utility_set_datatype}"
#endif
        END SUBROUTINE h5_utility_set_datatype


        !###########################################################
        SUBROUTINE separate_string(text, out_text)
            CHARACTER(LEN=*), INTENT(IN) :: text
            CHARACTER(LEN=:), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: out_text
            INTEGER :: length, count, i, start, end

#ifdef DEBUGhdf5
WRITE(error_unit,*) "mod_h5_utility/separate_string/{"
#endif

            length = LEN_TRIM(text)
            IF(text(length:length)=='/') length = length-1 ! Ignoring trailing '/'.

            ! Counting no. of '/'.
            count = 0
            DO i=1,length
               IF(text(i:i)=='/') count = count+1
            END DO

            ! Allocating out_text with assumed length 60 characters.
            ALLOCATE(CHARACTER(LEN=60) :: out_text(count+1))

            ! Separating each string
            start = 1
            DO i=1,count
                end = INDEX(text(start:),'/')+start-1
                out_text(i) = text(start:end-1)
                start = end+1
            END DO
            out_text(i) = text(start:length)

#ifdef DEBUGhdf5
WRITE(error_unit,*) "mod_h5_utility/separate_string/}"
#endif
        END SUBROUTINE separate_string


        !###########################################################
        FUNCTION is_file_already_opened(file_name) RESULT(index)
            IMPLICIT NONE
            CHARACTER(LEN=*), INTENT(IN) :: file_name
            INTEGER :: index
            INTEGER :: i

#ifdef DEBUGhdf5
WRITE(error_unit,*) "mod_h5_utility/is_file_already_opened/{"
#endif

            index = 0
            DO i=1,count_files
                IF(opened_file_objects(i)%p%file_name == file_name) index = i
            END DO  

#ifdef DEBUGhdf5
WRITE(error_unit,*) "mod_h5_utility/is_file_already_opened/}"
#endif
        END FUNCTION is_file_already_opened


        !###########################################################
        SUBROUTINE extend_opened_file_objects()
            TYPE(h5_dataset_type_ptr), POINTER :: temp(:) ! pointer to array of type h5_dataset_type_ptr
            INTEGER :: n, i

#ifdef DEBUGhdf5
WRITE(error_unit,*) "mod_h5_utility/extend_opened_file_objects/{"
#endif

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

#ifdef DEBUGhdf5
WRITE(error_unit,*) "mod_h5_utility/extend_opened_file_objects/}"
#endif
        END SUBROUTINE extend_opened_file_objects


        !###########################################################
        SUBROUTINE abort_program()
            INTEGER :: error

            IF(h5_utility_mpi) THEN
                CALL MPI_ABORT(MPI_COMM_WORLD, 1, error)
            ELSE
                STOP
            END IF
        END SUBROUTINE abort_program


END MODULE mod_h5_utility
