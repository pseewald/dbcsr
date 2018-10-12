!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2018  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief DBCSR tensor framework for block-sparse tensor contraction: Types and create/destroy
!>        routines.
!> \author Patrick Seewald
! **************************************************************************************************
MODULE dbcsr_tensor_types

#:include "dbcsr_tensor.fypp"
#:set maxdim = maxrank
#:set ndims = range(2,maxdim+1)

   USE dbcsr_array_list_methods,        ONLY: array_list,&
                                              array_offsets,&
                                              create_array_list,&
                                              destroy_array_list,&
                                              get_array_elements,&
                                              sizes_of_arrays,&
                                              sum_of_arrays,&
                                              array_sublist,&
                                              get_arrays
   USE dbcsr_api,                       ONLY: &
        dbcsr_distribution_get, &
        dbcsr_distribution_type, &
        dbcsr_get_info, dbcsr_type
   USE dbcsr_rectangular_base,          ONLY: &
      dbcsr_r_create, dbcsr_r_distribution_new,&
      dbcsr_r_distribution_destroy, dbcsr_r_finalize, dbcsr_r_get_info, &
      dbcsr_r_destroy, dbcsr_r_get_stored_coordinates
   USE dbcsr_rectangular_base_types,    ONLY: dbcsr_r_type, dbcsr_r_distribution_type
   USE dbcsr_tensor_index,              ONLY: get_2d_indices,&
                                              get_nd_indices,&
                                              create_nd_to_2d_mapping,&
                                              destroy_nd_to_2d_mapping,&
                                              get_mapping_info,&
                                              nd_to_2d_mapping,&
                                              split_index,&
                                              combine_index
   USE dbcsr_kinds,                     ONLY: default_string_length, int_8
   USE dbcsr_mpiwrap,                   ONLY: mp_cart_create,&
                                              mp_cart_rank,&
                                              mp_environ, &
                                              mp_dims_create, &
                                              mp_comm_free
   USE dbcsr_rectangular_global,        ONLY: dbcsr_r_distribution, dbcsr_r_rowcol_data
#include "base/dbcsr_base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbcsr_tensor_types'

   PUBLIC  :: &
      dbcsr_t_blk_offsets, &
      dbcsr_t_blk_sizes, &
      dbcsr_t_create, &
      dbcsr_t_destroy, &
      dbcsr_t_distribution, &
      dbcsr_t_distribution_destroy, &
      dbcsr_t_distribution_new, &
      dbcsr_t_distribution_type, &
      dbcsr_t_get_data_type, &
      dbcsr_t_nd_mp_comm, &
      dbcsr_t_nd_mp_free, &
      dbcsr_t_type, &
      dims_tensor, &
      ndims_tensor, &
      dbcsr_t_get_stored_coordinates

   TYPE dbcsr_t_type
      TYPE(dbcsr_r_type)                     :: matrix_rep
      TYPE(nd_to_2d_mapping)                 :: nd_index_blk
      TYPE(nd_to_2d_mapping)                 :: nd_index
      TYPE(array_list)                       :: blk_sizes
      TYPE(array_list)                       :: blk_offsets
      TYPE(array_list)                       :: nd_dist
      INTEGER                                :: comm_nd
      LOGICAL                                :: comm_prv = .FALSE.
      LOGICAL                                :: valid = .FALSE.
      CHARACTER(LEN=default_string_length)   :: name
   END TYPE dbcsr_t_type

   TYPE dbcsr_t_distribution_type
      TYPE(dbcsr_r_distribution_type) :: dist
      INTEGER                       :: comm_nd
      TYPE(array_list)              :: nd_dist
   END TYPE

   TYPE, EXTENDS(dbcsr_r_distribution) :: dbcsr_r_dist_t
      INTEGER, DIMENSION(:), ALLOCATABLE :: dims ! tensor dimensions only for this matrix dimensions
      INTEGER, DIMENSION(:), ALLOCATABLE :: dims_grid ! grid dimensions only for this matrix dimension
      TYPE(array_list)       :: nd_dist ! dist only for tensor dimensions belonging to this matrix dimension
   CONTAINS
      PROCEDURE :: dist => r_dist_t
      PROCEDURE :: rowcols => r_rowcols_t
   END TYPE

   TYPE, EXTENDS(dbcsr_r_rowcol_data) :: dbcsr_r_blk_size_t
      INTEGER, DIMENSION(:), ALLOCATABLE :: dims ! tensor dimensions only for this matrix dimensions
      TYPE(array_list) :: blk_size
   CONTAINS
      PROCEDURE :: data => r_blk_size_t
   END TYPE

   INTERFACE dbcsr_t_create
      MODULE PROCEDURE dbcsr_t_create_new
      MODULE PROCEDURE dbcsr_t_create_template
   END INTERFACE

   INTERFACE dbcsr_r_dist_t
      MODULE PROCEDURE new_dbcsr_r_dist_t
   END INTERFACE

   INTERFACE dbcsr_r_blk_size_t
      MODULE PROCEDURE new_dbcsr_r_blk_size_t
   END INTERFACE

CONTAINS

   FUNCTION new_dbcsr_r_dist_t(nd_dist, map_blks, map_grid, which_dim)
      TYPE(array_list), INTENT(IN)                   :: nd_dist
      TYPE(nd_to_2d_mapping), INTENT(IN)             :: map_blks, map_grid
      INTEGER, INTENT(IN) :: which_dim
      TYPE(dbcsr_r_dist_t) :: new_dbcsr_r_dist_t
      INTEGER, DIMENSION(2) :: grid_dims
      INTEGER(KIND=int_8), DIMENSION(2) :: matrix_dims
      INTEGER, DIMENSION(:), ALLOCATABLE :: index_map

      IF(which_dim == 1) THEN
         CALL get_mapping_info(map_blks, &
                               dims_2d_i8 = matrix_dims, &
                               map1_2d=index_map, &
                               dims1_2d=new_dbcsr_r_dist_t%dims)
         CALL get_mapping_info(map_grid, &
                               dims_2d = grid_dims, &
                               dims1_2d=new_dbcsr_r_dist_t%dims_grid)
      ELSEIF(which_dim == 2) THEN
         CALL get_mapping_info(map_blks, &
                               dims_2d_i8 = matrix_dims, &
                               map2_2d=index_map, &
                               dims2_2d=new_dbcsr_r_dist_t%dims)
         CALL get_mapping_info(map_grid, &
                               dims_2d = grid_dims, &
                               dims2_2d=new_dbcsr_r_dist_t%dims_grid)
      ELSE
         DBCSR_ABORT("Unknown value for which_dim")
      ENDIF

      new_dbcsr_r_dist_t%nd_dist = array_sublist(nd_dist, index_map)
      new_dbcsr_r_dist_t%nprowcol = grid_dims(which_dim)
      new_dbcsr_r_dist_t%nmrowcol = matrix_dims(which_dim)
   END FUNCTION

   FUNCTION r_dist_t(t, rowcol)
      CLASS(dbcsr_r_dist_t), INTENT(IN) :: t
      INTEGER(KIND=int_8), INTENT(IN) :: rowcol
      INTEGER, DIMENSION(SIZE(t%dims)) :: ind_blk
      INTEGER, DIMENSION(SIZE(t%dims)) :: dist_blk
      INTEGER :: r_dist_t

      ind_blk(:) = split_index(rowcol, t%dims, base=1, col_major=.TRUE.)
      dist_blk(:) = get_array_elements(t%nd_dist, ind_blk)
      r_dist_t = INT(combine_index(dist_blk, t%dims_grid, base=0, col_major=.FALSE.))
   END FUNCTION

   FUNCTION r_rowcols_t(t, dist)
      CLASS(dbcsr_r_dist_t), INTENT(IN) :: t
      INTEGER, INTENT(IN) :: dist
      INTEGER(KIND=int_8), DIMENSION(:), ALLOCATABLE :: r_rowcols_t
      INTEGER, DIMENSION(SIZE(t%dims)) :: dist_blk
      INTEGER, DIMENSION(:), ALLOCATABLE :: ${varlist("dist")}$, ${varlist("blks")}$, blks_tmp, nd_ind
      INTEGER :: ${varlist("i")}$, i, iblk, iblk_count, nblks
      INTEGER(KIND=int_8) :: nrowcols
      TYPE(array_list) :: blks

      dist_blk(:) = split_index(INT(dist, int_8), t%dims_grid, base=0, col_major=.FALSE.)

#:for ndim in ndims
      IF (SIZE(t%dims) == ${ndim}$) THEN
         CALL get_arrays(t%nd_dist, ${varlist("dist", nmax=ndim)}$)
      ENDIF
#:endfor

#:for idim in range(1, maxdim+1)
      IF (SIZE(t%dims) .GE. ${idim}$) THEN
         nblks = SIZE(dist_${idim}$)
         ALLOCATE (blks_tmp(nblks))
         iblk_count = 0
         DO iblk = 1, nblks
            IF(dist_${idim}$(iblk) == dist_blk(${idim}$)) THEN
               iblk_count = iblk_count+1
               blks_tmp(iblk_count) = iblk
            ENDIF
         ENDDO
         ALLOCATE(blks_${idim}$(iblk_count))
         blks_${idim}$(:) = blks_tmp(:iblk_count)
         DEALLOCATE(blks_tmp)
      ENDIF
#:endfor

#:for ndim in ndims
      IF (SIZE(t%dims) == ${ndim}$) THEN
         CALL create_array_list(blks, ${ndim}$, ${varlist("blks", nmax=ndim)}$)
      ENDIF
#:endfor

      nrowcols = PRODUCT(INT(sizes_of_arrays(blks), int_8))
      ALLOCATE (r_rowcols_t(nrowcols))

#:for ndim in ndims
      IF (SIZE(t%dims) == ${ndim}$) THEN
         ALLOCATE(nd_ind(${ndim}$))
         i = 0
#:for idim in range(1,ndim+1)
         DO i_${idim}$ = 1, SIZE(blks_${idim}$)
#:endfor
            i = i + 1

            nd_ind(:) = get_array_elements(blks, [${varlist("i", nmax=ndim)}$])
            r_rowcols_t(i) = combine_index(nd_ind, t%dims, base=1, col_major=.TRUE.)
#:for idim in range(1,ndim+1)
         ENDDO
#:endfor
      ENDIF
#:endfor

   END FUNCTION

   FUNCTION new_dbcsr_r_blk_size_t(blk_size, map_blks, which_dim)
      TYPE(array_list), INTENT(IN)                   :: blk_size
      TYPE(nd_to_2d_mapping), INTENT(IN)             :: map_blks
      INTEGER, INTENT(IN) :: which_dim
      INTEGER(KIND=int_8), DIMENSION(2) :: matrix_dims
      INTEGER, DIMENSION(:), ALLOCATABLE :: index_map
      TYPE(dbcsr_r_blk_size_t) :: new_dbcsr_r_blk_size_t

      IF (which_dim == 1) THEN
         CALL get_mapping_info(map_blks, &
                               dims_2d_i8=matrix_dims, &
                               map1_2d=index_map, &
                               dims1_2d=new_dbcsr_r_blk_size_t%dims)
      ELSEIF (which_dim == 2) THEN
         CALL get_mapping_info(map_blks, &
                               dims_2d_i8=matrix_dims, &
                               map2_2d=index_map, &
                               dims2_2d=new_dbcsr_r_blk_size_t%dims)
      ELSE
         DBCSR_ABORT("Unknown value for which_dim")
      ENDIF

      new_dbcsr_r_blk_size_t%blk_size = array_sublist(blk_size, index_map)
      new_dbcsr_r_blk_size_t%nmrowcol = matrix_dims(which_dim)

      new_dbcsr_r_blk_size_t%nfullrowcol = PRODUCT(INT(sum_of_arrays(new_dbcsr_r_blk_size_t%blk_size), &
                                                   KIND=int_8))
   END FUNCTION

   FUNCTION r_blk_size_t(t, rowcol)
      CLASS(dbcsr_r_blk_size_t), INTENT(IN) :: t
      INTEGER(KIND=int_8), INTENT(IN) :: rowcol
      INTEGER :: r_blk_size_t
      INTEGER, DIMENSION(SIZE(t%dims)) :: ind_blk
      INTEGER, DIMENSION(SIZE(t%dims)) :: blk_size

      ind_blk(:) = split_index(rowcol, t%dims, base=1, col_major=.TRUE.)
      blk_size(:) = get_array_elements(t%blk_size, ind_blk)
      r_blk_size_t = PRODUCT(blk_size)

   END FUNCTION

! **************************************************************************************************
!> \brief Create a default nd process topology that is consistent with a given 2d topology.
!>        Purpose: a nd tensor defined on the returned process grid can be represented as a DBCSR
!>        matrix with the given 2d topology.
!>        This is needed to enable contraction of 2 tensors (must have the same 2d process grid).
!> \param comm_2d communicator with 2-dimensional topology
!> \param map1_2d which nd-indices map to first matrix index and in which order
!> \param map2_2d which nd-indices map to second matrix index and in which order
!> \param dims_nd nd dimensions
!> \return with nd cartesian grid
! **************************************************************************************************
   FUNCTION dbcsr_t_nd_mp_comm(comm_2d, map1_2d, map2_2d, dims_nd)
      INTEGER, INTENT(IN)                               :: comm_2d
      INTEGER, DIMENSION(:), INTENT(IN)                 :: map1_2d, map2_2d
      INTEGER                                           :: ndim1, ndim2
      INTEGER                                           :: numtask
      INTEGER, DIMENSION(2)                             :: dims_2d, task_coor
      INTEGER, DIMENSION(SIZE(map1_2d) + SIZE(map2_2d)), &
         INTENT(IN), OPTIONAL                           :: dims_nd
      INTEGER, DIMENSION(SIZE(map1_2d))                 :: dims1_nd
      INTEGER, DIMENSION(SIZE(map2_2d))                 :: dims2_nd
      INTEGER, DIMENSION(SIZE(map1_2d) + SIZE(map2_2d)) :: dims_nd_prv, pos
      INTEGER                                           :: dbcsr_t_nd_mp_comm

      ndim1 = SIZE(map1_2d); ndim2 = SIZE(map2_2d)

      CALL mp_environ(numtask, dims_2d, task_coor, comm_2d)

      IF(.NOT. PRESENT(dims_nd)) THEN
         dims1_nd = 0; dims2_nd = 0
         CALL mp_dims_create(dims_2d(1), dims1_nd)
         CALL mp_dims_create(dims_2d(2), dims2_nd)
         dims_nd_prv(map1_2d) = dims1_nd
         dims_nd_prv(map2_2d) = dims2_nd
      ELSE
         DBCSR_ASSERT(PRODUCT(dims_nd(map1_2d)) == dims_2d(1))
         DBCSR_ASSERT(PRODUCT(dims_nd(map2_2d)) == dims_2d(2))
         dims_nd_prv = dims_nd
      ENDIF

      CALL mp_cart_create(comm_2d, ndim1+ndim2, dims_nd_prv, pos, dbcsr_t_nd_mp_comm)

   END FUNCTION

! **************************************************************************************************
!> \brief Release the MPI communicator.
!>
! **************************************************************************************************
   SUBROUTINE dbcsr_t_nd_mp_free(mp_comm)
     INTEGER, INTENT(INOUT)                               :: mp_comm

     CALL mp_comm_free(mp_comm)
   END SUBROUTINE dbcsr_t_nd_mp_free

! **************************************************************************************************
!> \brief Create a tensor distribution.
!>
!> comm_nd and comm_2d must be consistent (the nd process grid must map to the 2d process grid
!> according to map1_2d and map2_2d).
!> comm_2d is only needed for an early testing phase of tensor contraction and in the future,
!> comm_2d should not be an input.
!> \param dist ...
!> \param comm_nd communicator with n-dimensional Cartesian topology
!> \param comm_2d communicator with 2-dimensional topology
!> \param map1_2d which nd-indices map to first matrix index and in which order
!> \param map2_2d which nd-indices map to second matrix index and in which order
!> \param nd_dist_1 distribution vector for first dimension
!> \param nd_dist_2 distribution vector for second dimension
!> \param nd_dist_3 ...
! **************************************************************************************************
   SUBROUTINE dbcsr_t_distribution_new(dist, comm_nd, comm_2d, map1_2d, map2_2d, ${varlist("nd_dist")}$)

      TYPE(dbcsr_t_distribution_type), INTENT(OUT)    :: dist
      INTEGER, INTENT(IN)                             :: comm_nd, comm_2d
      INTEGER, DIMENSION(:), INTENT(IN)               :: map1_2d, map2_2d
      INTEGER, DIMENSION(:), INTENT(IN), OPTIONAL     :: ${varlist("nd_dist")}$
      INTEGER                                         :: ndims
      INTEGER, DIMENSION(2)                           :: pdims_2d_check, &
                                                         pdims_2d, task_coor_2d
      INTEGER, DIMENSION(SIZE(map1_2d)+SIZE(map2_2d)) :: dims, nblks_nd, task_coor
      LOGICAL, DIMENSION(SIZE(map1_2d)+SIZE(map2_2d)) :: periods
      LOGICAL, DIMENSION(2)                           :: periods_2d
      TYPE(array_list)                                :: nd_dist
      TYPE(nd_to_2d_mapping)                          :: map_blks, map_grid
      INTEGER                                         :: handle
      TYPE(dbcsr_r_dist_t)                            :: row_dist_obj, col_dist_obj
      CHARACTER(LEN=*), PARAMETER :: routineN = 'dbcsr_t_distribution_new', &
         routineP = moduleN//':'//routineN

      CALL timeset(routineN, handle)
      ndims = SIZE(map1_2d) + SIZE(map2_2d)
      DBCSR_ASSERT(ndims .GE. 2 .AND. ndims .LE. ${maxdim}$)

      CALL create_array_list(nd_dist, ndims, ${varlist("nd_dist")}$)

      nblks_nd(:) = sizes_of_arrays(nd_dist)

      ! check that 2d process topology is consistent with nd topology.
      CALL mp_environ(comm_nd, ndims, dims, task_coor, periods)

      ! process grid index mapping
      CALL create_nd_to_2d_mapping(map_grid, dims, map1_2d, map2_2d, base=0, col_major=.FALSE.)

      ! blk index mapping
      CALL create_nd_to_2d_mapping(map_blks, nblks_nd, map1_2d, map2_2d)

      row_dist_obj = dbcsr_r_dist_t(nd_dist, map_blks, map_grid, 1)
      col_dist_obj = dbcsr_r_dist_t(nd_dist, map_blks, map_grid, 2)

      CALL get_mapping_info(map_grid, dims_2d=pdims_2d)

      CALL mp_environ(comm_2d, 2, pdims_2d_check, task_coor_2d, periods_2d)
      IF (ANY(pdims_2d_check .NE. pdims_2d)) THEN
         DBCSR_ABORT("inconsistent process grid dimensions")
      ENDIF

      CALL dbcsr_r_distribution_new(dist%dist, comm_2d, row_dist_obj, col_dist_obj)

      !! create 2d distribution from nd distribution
      !ALLOCATE (row_dist(bdims_2d(1)), col_dist(bdims_2d(2)))
      !row_dist(:) = -1
      !col_dist(:) = -1

      !DO brow = 1, bdims_2d(1) ! blk rows
      !   DO bcol = 1, bdims_2d(2) ! blk cols

      !      ! 2d blk index to nd blk index
      !      ind(:) = get_nd_indices(map_blks, [brow, bcol])

      !      ! nd distribution to 2d distribution, mapped using proc grid
      !      dist_2d(:) = get_2d_indices(map_grid, get_array_elements(nd_dist, ind))

      !      ! first index maps to columns, second index to rows (row-major)
      !      IF (row_dist(brow) .LT. 0) row_dist(brow) = dist_2d(1)
      !      IF (col_dist(bcol) .LT. 0) col_dist(bcol) = dist_2d(2)
      !   ENDDO
      !ENDDO

      !CALL dbcsr_distribution_new(dist%dist, group=comm_2d, row_dist=row_dist, col_dist=col_dist, reuse_arrays=.TRUE.)
      dist%nd_dist = nd_dist
      dist%comm_nd = comm_nd
      CALL timestop(handle)
   END SUBROUTINE

! **************************************************************************************************
!> \brief Destroy tensor distribution
!> \param dist ...
! **************************************************************************************************
   SUBROUTINE dbcsr_t_distribution_destroy(dist)
      TYPE(dbcsr_t_distribution_type), INTENT(INOUT) :: dist
      INTEGER                                   :: handle
      CHARACTER(LEN=*), PARAMETER :: routineN = 'dbcsr_t_distribution_destroy', &
         routineP = moduleN//':'//routineN

      CALL timeset(routineN, handle)
      !CALL dbcsr_distribution_release(dist%dist)
      CALL dbcsr_r_distribution_destroy(dist%dist)
      CALL destroy_array_list(dist%nd_dist)
      CALL timestop(handle)
   END SUBROUTINE

   FUNCTION dbcsr_t_distribution(tensor)
      TYPE(dbcsr_t_type), INTENT(IN)  :: tensor
      TYPE(dbcsr_t_distribution_type) :: dbcsr_t_distribution

      CALL dbcsr_r_get_info(tensor%matrix_rep, distribution=dbcsr_t_distribution%dist)
      dbcsr_t_distribution%comm_nd = tensor%comm_nd
      dbcsr_t_distribution%nd_dist = tensor%nd_dist
   END FUNCTION

! **************************************************************************************************
!> \brief create a tensor
!> \param tensor ...
!> \param name ...
!> \param dist ...
!> \param map1_2d which nd-indices to map to first 2d index and in which order
!> \param map2_2d which nd-indices to map to first 2d index and in which order
!> \param data_type ...
!> \param blk_size_1 blk sizes in each dimension
!> \param blk_size_2 ...
! **************************************************************************************************
   SUBROUTINE dbcsr_t_create_new(tensor, name, dist, map1_2d, map2_2d, data_type, &
                                 ${varlist("blk_size")}$)
      TYPE(dbcsr_t_type), INTENT(OUT)                   :: tensor
      CHARACTER(len=*), INTENT(IN)                      :: name
      TYPE(dbcsr_t_distribution_type), INTENT(INOUT)    :: dist
      INTEGER, DIMENSION(:), INTENT(IN)                 :: map1_2d, map2_2d
      INTEGER, INTENT(IN), OPTIONAL                     :: data_type
      INTEGER, DIMENSION(:), INTENT(IN), OPTIONAL       :: ${varlist("blk_size")}$
      INTEGER                                           :: ndims
      INTEGER(KIND=int_8), DIMENSION(2)                             :: dims_2d
      INTEGER, DIMENSION(SIZE(map1_2d) + SIZE(map2_2d)) :: dims
      TYPE(dbcsr_r_blk_size_t)                          :: col_blk_size_obj, row_blk_size_obj
      TYPE(array_list)                                  :: blk_size
      TYPE(nd_to_2d_mapping)                            :: map
      INTEGER                                   :: handle
      CHARACTER(LEN=*), PARAMETER :: routineN = 'dbcsr_t_create_new', &
         routineP = moduleN//':'//routineN

      CALL timeset(routineN, handle)
      ndims = SIZE(map1_2d)+SIZE(map2_2d)
      CALL create_array_list(blk_size, ndims, ${varlist("blk_size")}$)
      dims = sizes_of_arrays(blk_size)

      CALL create_nd_to_2d_mapping(map, dims, map1_2d, map2_2d)
      CALL get_mapping_info(map, dims_2d_i8=dims_2d)

      row_blk_size_obj = dbcsr_r_blk_size_t(blk_size, map, 1)
      col_blk_size_obj = dbcsr_r_blk_size_t(blk_size, map, 2)

      CALL dbcsr_r_create(matrix=tensor%matrix_rep, &
                          name=name//"_matrix_rep", &
                          dist=dist%dist, &
                          row_blk_size=row_blk_size_obj, &
                          col_blk_size=col_blk_size_obj, &
                          data_type=data_type)

      tensor%nd_index_blk = map
      tensor%name = name

      CALL dbcsr_r_finalize(tensor%matrix_rep)
      CALL destroy_nd_to_2d_mapping(map)

      ! map element-wise tensor index
      CALL create_nd_to_2d_mapping(map, sum_of_arrays(blk_size), map1_2d, map2_2d)
      tensor%nd_index = map
      tensor%blk_sizes = blk_size
      tensor%nd_dist = dist%nd_dist
      tensor%comm_nd = dist%comm_nd

      CALL array_offsets(tensor%blk_sizes, tensor%blk_offsets)

      tensor%valid = .TRUE.
      CALL timestop(handle)
   END SUBROUTINE

! **************************************************************************************************
!> \brief create a tensor from template
!> \param tensor_in ...
!> \param tensor ...
!> \param name ...
! **************************************************************************************************
   SUBROUTINE dbcsr_t_create_template(tensor_in, tensor, name)
      TYPE(dbcsr_t_type), INTENT(INOUT)      :: tensor_in
      TYPE(dbcsr_t_type), INTENT(OUT)        :: tensor
      CHARACTER(len=*), INTENT(IN), OPTIONAL :: name
      INTEGER                                   :: handle
      CHARACTER(LEN=*), PARAMETER :: routineN = 'dbcsr_t_create_template', &
         routineP = moduleN//':'//routineN

      CALL timeset(routineN, handle)
      CALL dbcsr_r_create(tensor_in%matrix_rep, tensor%matrix_rep)
      CALL dbcsr_r_finalize(tensor%matrix_rep)

      tensor%nd_index_blk = tensor_in%nd_index_blk
      tensor%nd_index = tensor_in%nd_index
      tensor%blk_sizes = tensor_in%blk_sizes
      tensor%blk_offsets = tensor_in%blk_offsets
      tensor%nd_dist = tensor_in%nd_dist
      tensor%comm_nd = tensor_in%comm_nd
      tensor%valid = .TRUE.
      IF (PRESENT(name)) THEN
         tensor%name = name
      ELSE
         tensor%name = tensor_in%name
      ENDIF
      CALL timestop(handle)
   END SUBROUTINE

! **************************************************************************************************
!> \brief Create 2-rank tensor from matrix.
!> \param matrix_in ...
!> \param tensor ...
!> \param order ...
!> \param name ...
! **************************************************************************************************
   SUBROUTINE dbcsr_t_create_matrix(matrix_in, tensor, order, name)
      TYPE(dbcsr_type), INTENT(IN)                :: matrix_in
      TYPE(dbcsr_t_type), INTENT(OUT)             :: tensor
      INTEGER, DIMENSION(2), INTENT(IN), OPTIONAL :: order
      CHARACTER(len=*), INTENT(IN), OPTIONAL      :: name

      CHARACTER(len=default_string_length)        :: name_in
      INTEGER, DIMENSION(2)                       :: order_in
      INTEGER                                     :: comm_nd, comm_2d, data_type
      TYPE(dbcsr_distribution_type)                :: matrix_dist
      TYPE(dbcsr_t_distribution_type)             :: dist
      INTEGER, DIMENSION(:), POINTER              :: row_blk_size, col_blk_size
      INTEGER, DIMENSION(:), POINTER              :: col_dist, row_dist
      INTEGER                                   :: handle
      CHARACTER(LEN=*), PARAMETER :: routineN = 'dbcsr_t_create_matrix', &
         routineP = moduleN//':'//routineN

      CALL timeset(routineN, handle)
      IF (PRESENT(name)) THEN
         name_in = name
      ELSE
         CALL dbcsr_get_info(matrix_in, name=name_in)
      ENDIF

      IF (PRESENT(order)) THEN
         order_in = order
      ELSE
         order_in = [1, 2]
      ENDIF

      CALL dbcsr_get_info(matrix_in, distribution=matrix_dist)
      CALL dbcsr_distribution_get(matrix_dist, group=comm_2d, row_dist=row_dist, col_dist=col_dist)
      comm_nd = dbcsr_t_nd_mp_comm(comm_2d, [order_in(1)], [order_in(2)])

      CALL dbcsr_t_distribution_new( &
         dist, &
         comm_nd, comm_2d, &
         [order_in(1)], [order_in(2)], &
         row_dist, col_dist)

      CALL dbcsr_get_info(matrix_in, &
                          data_type=data_type, &
                          row_blk_size=row_blk_size, &
                          col_blk_size=col_blk_size)

      CALL dbcsr_t_create_new(tensor, name_in, dist, &
                              [order_in(1)], [order_in(2)], &
                              data_type, &
                              row_blk_size, &
                              col_blk_size)
      CALL timestop(handle)
   END SUBROUTINE

! **************************************************************************************************
!> \brief Destroy a tensor
!> \param tensor ...
! **************************************************************************************************
   SUBROUTINE dbcsr_t_destroy(tensor)
      TYPE(dbcsr_t_type), INTENT(INOUT) :: tensor
      INTEGER                                   :: handle
      CHARACTER(LEN=*), PARAMETER :: routineN = 'dbcsr_t_destroy', &
         routineP = moduleN//':'//routineN

      CALL timeset(routineN, handle)
      CALL dbcsr_r_destroy(tensor%matrix_rep)
      CALL destroy_nd_to_2d_mapping(tensor%nd_index_blk)
      CALL destroy_nd_to_2d_mapping(tensor%nd_index)
      CALL destroy_array_list(tensor%blk_sizes)
      CALL destroy_array_list(tensor%blk_offsets)
      CALL destroy_array_list(tensor%nd_dist)
      IF (tensor%comm_prv) THEN
         ! release communicator if it is marked as private (not created by user)
         CALL mp_comm_free(tensor%comm_nd)
      ENDIF
      tensor%comm_prv = .FALSE.
      tensor%valid = .FALSE.
      tensor%name = ""
      CALL timestop(handle)
   END SUBROUTINE

! **************************************************************************************************
!> \brief tensor rank
!> \param tensor ...
!> \return ...
! **************************************************************************************************
   PURE FUNCTION ndims_tensor(tensor)
      TYPE(dbcsr_t_type), INTENT(IN) :: tensor
      INTEGER                        :: ndims_tensor

      ndims_tensor = tensor%nd_index%ndim_nd
   END FUNCTION

! **************************************************************************************************
!> \brief tensor dimensions
!> \param tensor ...
!> \param dims ...
! **************************************************************************************************
   SUBROUTINE dims_tensor(tensor, dims)
      TYPE(dbcsr_t_type), INTENT(IN)              :: tensor
      INTEGER, DIMENSION(ndims_tensor(tensor)), &
         INTENT(OUT)                              :: dims

      DBCSR_ASSERT(tensor%valid)
      dims = tensor%nd_index%dims_nd
   END SUBROUTINE

! **************************************************************************************************
!> \brief tensor data type
!> \param tensor ...
!> \return ...
! **************************************************************************************************
   FUNCTION dbcsr_t_get_data_type(tensor) RESULT(data_type)
      TYPE(dbcsr_t_type), INTENT(IN) :: tensor
      INTEGER                        :: data_type

      CALL dbcsr_r_get_info(tensor%matrix_rep, data_type=data_type)
   END FUNCTION

! **************************************************************************************************
!> \brief Size of tensor block
!> \param tensor ...
!> \param ind block index ...
!> \param blk_size block size ...
! **************************************************************************************************
   SUBROUTINE dbcsr_t_blk_sizes(tensor, ind, blk_size)
      TYPE(dbcsr_t_type), INTENT(IN)              :: tensor
      INTEGER, DIMENSION(ndims_tensor(tensor)), &
         INTENT(IN)                               :: ind
      INTEGER, DIMENSION(ndims_tensor(tensor)), &
         INTENT(OUT)                              :: blk_size

      blk_size(:) = get_array_elements(tensor%blk_sizes, ind)
   END SUBROUTINE

! **************************************************************************************************
!> \brief offset of tensor block
!> \param tensor ...
!> \param ind block index
!> \param blk_offset block offset
! **************************************************************************************************
   SUBROUTINE dbcsr_t_blk_offsets(tensor, ind, blk_offset)
      TYPE(dbcsr_t_type), INTENT(IN)              :: tensor
      INTEGER, DIMENSION(ndims_tensor(tensor)), &
         INTENT(IN)                               :: ind
      INTEGER, DIMENSION(ndims_tensor(tensor)), &
         INTENT(OUT)                              :: blk_offset

      DBCSR_ASSERT(tensor%valid)
      blk_offset(:) = get_array_elements(tensor%blk_offsets, ind)
   END SUBROUTINE

! **************************************************************************************************
!> \brief Generalization of dbcsr_get_stored_coordinates for tensors.
!> \param tensor ...
!> \param ind_nd ...
!> \param processor ...
! **************************************************************************************************
   SUBROUTINE dbcsr_t_get_stored_coordinates(tensor, ind_nd, processor)
      TYPE(dbcsr_t_type), INTENT(IN)               :: tensor
      INTEGER, DIMENSION(ndims_tensor(tensor)), &
         INTENT(IN)                                :: ind_nd
      INTEGER, INTENT(OUT)                         :: processor

      INTEGER(KIND=int_8), DIMENSION(2)                        :: ind_2d

      ind_2d(:) = get_2d_indices(tensor%nd_index_blk, ind_nd)
      CALL dbcsr_r_get_stored_coordinates(tensor%matrix_rep, ind_2d(1), ind_2d(2), processor)
   END SUBROUTINE

END MODULE
