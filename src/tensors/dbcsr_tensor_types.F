!--------------------------------------------------------------------------------------------------!
! Copyright (C) by the DBCSR developers group - All rights reserved                                !
! This file is part of the DBCSR library.                                                          !
!                                                                                                  !
! For information on the license, see the LICENSE file.                                            !
! For further information please visit https://dbcsr.cp2k.org                                      !
! SPDX-License-Identifier: GPL-2.0+                                                                !
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
                                              get_arrays,&
                                              get_ith_array
   USE dbcsr_api,                       ONLY: &
        dbcsr_distribution_get, &
        dbcsr_distribution_type, &
        dbcsr_get_info, dbcsr_type,&
        ${uselist(dtype_float_param)}$
   USE dbcsr_kinds,                     ONLY: &
        ${uselist(dtype_float_prec)}$,&
        default_string_length
   USE dbcsr_rectangular_base,          ONLY: &
      dbcsr_r_create, dbcsr_r_distribution_new,&
      dbcsr_r_distribution_destroy, dbcsr_r_finalize, dbcsr_r_get_info, &
      dbcsr_r_destroy, dbcsr_r_get_stored_coordinates, dbcsr_r_set, dbcsr_r_filter, &
      dbcsr_r_get_num_blocks, dbcsr_r_get_data_size, dbcsr_r_clear
   USE dbcsr_rectangular_base_types,    ONLY: dbcsr_r_type, dbcsr_r_distribution_type, dbcsr_r_split_info
   USE dbcsr_tensor_index,              ONLY: get_2d_indices,&
                                              get_nd_indices,&
                                              create_nd_to_2d_mapping,&
                                              destroy_nd_to_2d_mapping,&
                                              get_mapping_info,&
                                              nd_to_2d_mapping,&
                                              split_index,&
                                              combine_index,&
                                              ndims_mapping
   USE dbcsr_rectangular_split,         ONLY: dbcsr_r_split_rows_or_cols, &
                                              dbcsr_r_release_info
   USE dbcsr_kinds,                     ONLY: default_string_length, int_8
   USE dbcsr_mpiwrap,                   ONLY: mp_cart_create,&
                                              mp_cart_rank,&
                                              mp_environ, &
                                              mp_dims_create, &
                                              mp_comm_free, mp_comm_dup,&
                                              mp_sum, mp_max
   USE dbcsr_rectangular_global,        ONLY: dbcsr_r_distribution, dbcsr_r_rowcol_data
   USE dbcsr_allocate_wrap,             ONLY: allocate_any
#include "base/dbcsr_base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbcsr_tensor_types'

   PUBLIC  :: &
      blk_dims_tensor, &
      dbcsr_t_blk_offsets, &
      dbcsr_t_blk_sizes, &
      dbcsr_t_create, &
      dbcsr_t_destroy, &
      dbcsr_t_distribution, &
      dbcsr_t_distribution_destroy, &
      dbcsr_t_distribution_new, &
      dbcsr_t_distribution_type, &
      dbcsr_t_get_data_type, &
      dbcsr_t_get_stored_coordinates,&
      dbcsr_t_hold,&
      dbcsr_t_nd_mp_comm, &
      dbcsr_t_nd_mp_free, &
      dbcsr_t_type, &
      dims_tensor, &
      ndims_tensor,&
      dbcsr_t_pgrid_type,&
      mp_environ_pgrid,&
      dbcsr_t_pgrid_create,&
      dbcsr_t_pgrid_destroy, &
      dbcsr_t_set,&
      dbcsr_t_filter,&
      dbcsr_t_get_info,&
      dbcsr_t_write_tensor_dist, &
      dbcsr_t_write_tensor_info, &
      dbcsr_t_clear,&
      dbcsr_t_get_num_blocks,&
      dbcsr_t_finalize

   TYPE dbcsr_t_pgrid_type
      TYPE(nd_to_2d_mapping) :: nd_index_grid
      INTEGER                :: mp_comm_2d
   END TYPE

   TYPE dbcsr_t_type
      TYPE(dbcsr_r_type), POINTER            :: matrix_rep => NULL()
      TYPE(nd_to_2d_mapping)                 :: nd_index_blk
      TYPE(nd_to_2d_mapping)                 :: nd_index
      !TYPE(nd_to_2d_mapping)                 :: nd_index_grid
      TYPE(array_list)                       :: blk_sizes
      TYPE(array_list)                       :: blk_offsets
      TYPE(array_list)                       :: nd_dist
      !INTEGER                                :: comm_nd
      !INTEGER                                :: comm_2d
      TYPE(dbcsr_t_pgrid_type)               :: pgrid
      TYPE(array_list)                       :: blks_local
      INTEGER, DIMENSION(:), ALLOCATABLE     :: nblks_local
      INTEGER, DIMENSION(:), ALLOCATABLE     :: nfull_local
      LOGICAL                                :: valid = .FALSE.
      LOGICAL                                :: owns_matrix = .FALSE.
      CHARACTER(LEN=default_string_length)   :: name
      INTEGER, POINTER :: refcount => NULL() ! lightweight reference counting for communicators
   END TYPE dbcsr_t_type

   TYPE dbcsr_t_distribution_type
      TYPE(dbcsr_r_distribution_type) :: dist
      !INTEGER                       :: comm_nd
      !INTEGER                       :: comm_2d
      TYPE(dbcsr_t_pgrid_type)      :: pgrid
      TYPE(array_list)              :: nd_dist
      !TYPE(nd_to_2d_mapping)        :: nd_index_grid
      INTEGER, POINTER :: refcount => NULL() ! lightweight reference counting for communicators
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
      MODULE PROCEDURE dbcsr_t_create_matrix
   END INTERFACE

   INTERFACE dbcsr_r_dist_t
      MODULE PROCEDURE new_dbcsr_r_dist_t
   END INTERFACE

   INTERFACE dbcsr_r_blk_size_t
      MODULE PROCEDURE new_dbcsr_r_blk_size_t
   END INTERFACE

   INTERFACE dbcsr_t_set
#:for dparam, dtype, dsuffix in dtype_float_list
      MODULE PROCEDURE dbcsr_t_set_${dsuffix}$
#:endfor
   END INTERFACE

   INTERFACE dbcsr_t_filter
#:for dparam, dtype, dsuffix in dtype_float_list
      MODULE PROCEDURE dbcsr_t_filter_${dsuffix}$
#:endfor
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

#:for ndim in range(1, maxdim+1)
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

#:for ndim in range(1, maxdim+1)
      IF (SIZE(t%dims) == ${ndim}$) THEN
         CALL create_array_list(blks, ${ndim}$, ${varlist("blks", nmax=ndim)}$)
      ENDIF
#:endfor

      nrowcols = PRODUCT(INT(sizes_of_arrays(blks), int_8))
      ALLOCATE (r_rowcols_t(nrowcols))

#:for ndim in range(1, maxdim+1)
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
!> \param dims1_nd
!> \param dims2_nd
!> \param pdims_2d if comm_2d does not have a cartesian topology associated, can input dimensions with
!>                 pdims_2d (hack)
!> \return with nd cartesian grid
! **************************************************************************************************
   FUNCTION dbcsr_t_nd_mp_comm(comm_2d, map1_2d, map2_2d, dims_nd, dims1_nd, dims2_nd, pdims_2d)
      INTEGER, INTENT(IN)                               :: comm_2d
      INTEGER, DIMENSION(:), INTENT(IN)                 :: map1_2d, map2_2d
      INTEGER, DIMENSION(SIZE(map1_2d) + SIZE(map2_2d)), &
         INTENT(IN), OPTIONAL                           :: dims_nd
      INTEGER, DIMENSION(SIZE(map1_2d)), INTENT(IN), OPTIONAL :: dims1_nd
      INTEGER, DIMENSION(SIZE(map2_2d)), INTENT(IN), OPTIONAL :: dims2_nd
      INTEGER, DIMENSION(2), INTENT(IN), OPTIONAL           :: pdims_2d
      INTEGER                                           :: ndim1, ndim2
      INTEGER                                           :: numtask
      INTEGER, DIMENSION(2)                             :: dims_2d, task_coor

      INTEGER, DIMENSION(SIZE(map1_2d)) :: dims1_nd_prv
      INTEGER, DIMENSION(SIZE(map2_2d)) :: dims2_nd_prv
      INTEGER, DIMENSION(SIZE(map1_2d) + SIZE(map2_2d)) :: dims_nd_prv
      INTEGER                                           :: handle
      CHARACTER(LEN=*), PARAMETER :: routineN = 'dbcsr_t_nd_mp_comm', &
         routineP = moduleN//':'//routineN
      TYPE(dbcsr_t_pgrid_type)                          :: dbcsr_t_nd_mp_comm

      CALL timeset(routineN, handle)

      ndim1 = SIZE(map1_2d); ndim2 = SIZE(map2_2d)

      IF(PRESENT(pdims_2d)) THEN
         dims_2d(:) = pdims_2d
      ELSE
         CALL mp_environ(numtask, dims_2d, task_coor, comm_2d)
      ENDIF

      IF(.NOT. PRESENT(dims_nd)) THEN
         dims1_nd_prv = 0; dims2_nd_prv = 0
         IF(PRESENT(dims1_nd)) THEN
            dims1_nd_prv(:) = dims1_nd
         ELSE
            CALL mp_dims_create(dims_2d(1), dims1_nd_prv)
         ENDIF

         IF(PRESENT(dims2_nd)) THEN
            dims2_nd_prv(:) = dims2_nd
         ELSE
            CALL mp_dims_create(dims_2d(2), dims2_nd_prv)
         ENDIF
         dims_nd_prv(map1_2d) = dims1_nd_prv
         dims_nd_prv(map2_2d) = dims2_nd_prv
      ELSE
         DBCSR_ASSERT(PRODUCT(dims_nd(map1_2d)) == dims_2d(1))
         DBCSR_ASSERT(PRODUCT(dims_nd(map2_2d)) == dims_2d(2))
         dims_nd_prv = dims_nd
      ENDIF

      CALL dbcsr_t_pgrid_create(comm_2d, dims_nd_prv, dbcsr_t_nd_mp_comm, map1_2d, map2_2d)

      CALL timestop(handle)

   END FUNCTION

   !FUNCTION dbcsr_t_2d_mp_comm(comm_nd, map1_2d, map2_2d)
   !   INTEGER, INTENT(IN) :: comm_nd
   !   INTEGER, DIMENSION(:), INTENT(IN) :: map1_2d, map2_2d
   !   INTEGER, DIMENSION(2) :: pos, dims_2d
   !   INTEGER :: dbcsr_t_2d_mp_comm, handle
   !   INTEGER, DIMENSION(SIZE(map1_2d) + SIZE(map2_2d)) :: dims_nd, pcoord
   !   LOGICAL, DIMENSION(SIZE(map1_2d) + SIZE(map2_2d)) :: periods
   !   CHARACTER(LEN=*), PARAMETER :: routineN = 'dbcsr_t_2d_mp_comm', &
   !      routineP = moduleN//':'//routineN

   !   CALL timeset(routineN, handle)

   !   CALL mp_environ(comm_nd, SIZE(dims_nd), dims_nd, pcoord, periods)

   !   dims_2d = [PRODUCT(dims_nd(map1_2d)), PRODUCT(dims_nd(map2_2d))]

   !   CALL mp_cart_create(comm_nd, 2, dims_2d, pos, dbcsr_t_2d_mp_comm)

   !   CALL timestop(handle)

   !END FUNCTION

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
!> \param ngroup How many subgroups (splits) for longest matrix dimension
!> \param dimsplit which dimension to split (1: row, 2: column)
!> \param own_comm whether distribution should own communicator
! **************************************************************************************************
   SUBROUTINE dbcsr_t_distribution_new(dist, pgrid, map1_2d, map2_2d, ${varlist("nd_dist")}$, ngroup, dimsplit, own_comm)

      TYPE(dbcsr_t_distribution_type), INTENT(OUT)    :: dist
      TYPE(dbcsr_t_pgrid_type), INTENT(IN)            :: pgrid
      INTEGER, DIMENSION(:), INTENT(IN)               :: map1_2d, map2_2d
      INTEGER, DIMENSION(:), INTENT(IN), OPTIONAL     :: ${varlist("nd_dist")}$
      INTEGER, INTENT(IN), OPTIONAL                   :: ngroup
      INTEGER, INTENT(IN), OPTIONAL                   :: dimsplit
      LOGICAL, INTENT(IN), OPTIONAL                   :: own_comm
      INTEGER                                         :: ndims, comm_2d
      INTEGER, DIMENSION(2)                           :: pdims_2d_check, &
                                                         pdims_2d, task_coor_2d
      INTEGER, DIMENSION(SIZE(map1_2d)+SIZE(map2_2d)) :: dims, nblks_nd, task_coor
      LOGICAL, DIMENSION(2)                           :: periods_2d
      TYPE(array_list)                                :: nd_dist
      TYPE(nd_to_2d_mapping)                          :: map_blks, map_grid
      INTEGER                                         :: handle
      TYPE(dbcsr_r_dist_t)                            :: row_dist_obj, col_dist_obj
      TYPE(dbcsr_t_pgrid_type)                        :: pgrid_prv
      LOGICAL                                         :: need_pgrid_remap
      INTEGER,DIMENSION(:), ALLOCATABLE               :: map1_2d_check, map2_2d_check
      INTEGER                                         :: igroup, ngroup_prv
      TYPE(dbcsr_r_split_info)                        :: split_info
      CHARACTER(LEN=*), PARAMETER :: routineN = 'dbcsr_t_distribution_new', &
         routineP = moduleN//':'//routineN

      CALL timeset(routineN, handle)
      ndims = SIZE(map1_2d) + SIZE(map2_2d)
      DBCSR_ASSERT(ndims .GE. 2 .AND. ndims .LE. ${maxdim}$)

      CALL create_array_list(nd_dist, ndims, ${varlist("nd_dist")}$)

      nblks_nd(:) = sizes_of_arrays(nd_dist)

      need_pgrid_remap = .TRUE.
      IF(PRESENT(own_comm)) THEN
         CALL get_mapping_info(pgrid%nd_index_grid, map1_2d=map1_2d_check, map2_2d=map2_2d_check)
         IF(own_comm) THEN
            IF(.NOT. array_eq_i(map1_2d_check, map1_2d) .OR. .NOT. array_eq_i(map2_2d_check, map2_2d)) THEN
               DBCSR_ABORT("map1_2d / map2_2d are not consistent with pgrid")
            ENDIF
            pgrid_prv = pgrid
            need_pgrid_remap = .FALSE.
         ENDIF
      ENDIF

      IF(need_pgrid_remap) CALL dbcsr_t_pgrid_remap(pgrid, map1_2d, map2_2d, pgrid_prv)

      ! check that 2d process topology is consistent with nd topology.
      !CALL mp_environ(comm_nd, ndims, dims, task_coor, periods)
      CALL mp_environ_pgrid(pgrid_prv, dims, task_coor)

      ! process grid index mapping
      CALL create_nd_to_2d_mapping(map_grid, dims, map1_2d, map2_2d, base=0, col_major=.FALSE.)

      ! blk index mapping
      CALL create_nd_to_2d_mapping(map_blks, nblks_nd, map1_2d, map2_2d)

      row_dist_obj = dbcsr_r_dist_t(nd_dist, map_blks, map_grid, 1)
      col_dist_obj = dbcsr_r_dist_t(nd_dist, map_blks, map_grid, 2)

      CALL get_mapping_info(map_grid, dims_2d=pdims_2d)

      comm_2d = pgrid_prv%mp_comm_2d !dbcsr_t_2d_mp_comm(comm_nd, map1_2d, map2_2d)

      CALL mp_environ(comm_2d, 2, pdims_2d_check, task_coor_2d, periods_2d)
      IF (ANY(pdims_2d_check .NE. pdims_2d)) THEN
         PRINT *, pdims_2d, pdims_2d_check
         DBCSR_ABORT("inconsistent process grid dimensions")
      ENDIF

      IF(PRESENT(ngroup)) THEN
         DBCSR_ASSERT(PRESENT(dimsplit))

         IF(MOD(pdims_2d(dimsplit), ngroup) .NE. 0) THEN
            DBCSR_ABORT("split factor must divide process grid dimension")
         ENDIF

         igroup = task_coor_2d(dimsplit)/(pdims_2d(dimsplit)/ngroup)

         ngroup_prv = ngroup
         CALL dbcsr_r_split_rows_or_cols(split_info, comm_2d, ngroup_prv, igroup, dimsplit)

         DBCSR_ASSERT(ngroup_prv == ngroup)

         CALL dbcsr_r_distribution_new(dist%dist, comm_2d, row_dist_obj, col_dist_obj, split_info=split_info)
         CALL dbcsr_r_release_info(split_info)
      ELSE
         CALL dbcsr_r_distribution_new(dist%dist, comm_2d, row_dist_obj, col_dist_obj)
      ENDIF

      dist%nd_dist = nd_dist
      !dist%comm_2d = comm_2d
      dist%pgrid = pgrid_prv
      !IF(PRESENT(own_comm)) THEN
      !   IF(.NOT. own_comm) THEN
      !      CALL mp_comm_dup(comm_2d, dist%pgrid_prv%mp_comm_2d)
      !   ENDIF
      !ELSE
      !   CALL mp_comm_dup(comm_2d, dist%pgrid_prv%mp_comm_2d)
      !ENDIF

      !dist%nd_index_grid = map_grid

      ALLOCATE(dist%refcount)
      dist%refcount = 1
      CALL timestop(handle)

      CONTAINS
         PURE FUNCTION array_eq_i(arr1, arr2)
            INTEGER, INTENT(IN), DIMENSION(:) :: arr1
            INTEGER, INTENT(IN), DIMENSION(:) :: arr2
            LOGICAL                           :: array_eq_i

            array_eq_i = .FALSE.
            IF (SIZE(arr1) .EQ. SIZE(arr2)) array_eq_i = ALL(arr1 == arr2)

         END FUNCTION

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
      LOGICAL :: abort


      CALL timeset(routineN, handle)
      CALL dbcsr_r_distribution_destroy(dist%dist)
      CALL destroy_array_list(dist%nd_dist)
      !CALL destroy_nd_to_2d_mapping(dist%nd_index_grid)

      abort = .FALSE.
      IF (.NOT. ASSOCIATED(dist%refcount)) THEN
         abort = .TRUE.
      ELSEIF (dist%refcount < 1) THEN
         abort = .TRUE.
      ENDIF

      IF (abort) THEN
         DBCSR_ABORT("can not destroy non-existing tensor distribution")
      ENDIF

      dist%refcount = dist%refcount - 1

      IF(dist%refcount == 0) THEN
         CALL dbcsr_t_pgrid_destroy(dist%pgrid)
         !CALL mp_comm_free(dist%comm_2d)
         !CALL mp_comm_free(dist%comm_nd)
         DEALLOCATE(dist%refcount)
      ENDIF

      CALL timestop(handle)
   END SUBROUTINE

   SUBROUTINE dbcsr_t_distribution_hold(dist)
      TYPE(dbcsr_t_distribution_type), INTENT(IN) :: dist
      INTEGER, POINTER                            :: ref

      IF(dist%refcount < 1) THEN
         DBCSR_ABORT("can not hold non-existing tensor distribution")
      ENDIF
      ref => dist%refcount
      ref = ref + 1
   END SUBROUTINE

   FUNCTION dbcsr_t_distribution(tensor)
      TYPE(dbcsr_t_type), INTENT(IN)  :: tensor
      TYPE(dbcsr_t_distribution_type) :: dbcsr_t_distribution

      CALL dbcsr_r_get_info(tensor%matrix_rep, distribution=dbcsr_t_distribution%dist)
      !dbcsr_t_distribution%comm_nd = tensor%comm_nd
      !dbcsr_t_distribution%comm_2d = tensor%comm_2d
      dbcsr_t_distribution%pgrid = tensor%pgrid
      dbcsr_t_distribution%nd_dist = tensor%nd_dist
      !dbcsr_t_distribution%nd_index_grid = tensor%nd_index_grid
      dbcsr_t_distribution%refcount => dbcsr_t_distribution%refcount
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
      INTEGER, DIMENSION(SIZE(map1_2d) + SIZE(map2_2d)) :: dims, pdims, task_coor
      TYPE(dbcsr_r_blk_size_t)                          :: col_blk_size_obj, row_blk_size_obj
      TYPE(array_list)                                  :: blk_size, blks_local
      TYPE(nd_to_2d_mapping)                            :: map
      INTEGER                                   :: handle
      CHARACTER(LEN=*), PARAMETER :: routineN = 'dbcsr_t_create_new', &
         routineP = moduleN//':'//routineN
      INTEGER, DIMENSION(:), ALLOCATABLE              :: ${varlist("blks_local")}$
      INTEGER, DIMENSION(:), ALLOCATABLE              :: ${varlist("dist")}$
      INTEGER                                         :: iblk_count, iblk
      INTEGER, DIMENSION(:), ALLOCATABLE              :: nblks_local, nfull_local


      CALL timeset(routineN, handle)
      ndims = SIZE(map1_2d)+SIZE(map2_2d)
      CALL create_array_list(blk_size, ndims, ${varlist("blk_size")}$)
      dims = sizes_of_arrays(blk_size)

      CALL create_nd_to_2d_mapping(map, dims, map1_2d, map2_2d)
      CALL get_mapping_info(map, dims_2d_i8=dims_2d)

      row_blk_size_obj = dbcsr_r_blk_size_t(blk_size, map, 1)
      col_blk_size_obj = dbcsr_r_blk_size_t(blk_size, map, 2)

      ALLOCATE(tensor%matrix_rep)
      CALL dbcsr_r_create(matrix=tensor%matrix_rep, &
                          name=TRIM(name)//" matrix", &
                          dist=dist%dist, &
                          row_blk_size=row_blk_size_obj, &
                          col_blk_size=col_blk_size_obj, &
                          data_type=data_type)
      tensor%owns_matrix = .TRUE.

      tensor%nd_index_blk = map
      tensor%name = name

      CALL dbcsr_r_finalize(tensor%matrix_rep)
      CALL destroy_nd_to_2d_mapping(map)

      ! map element-wise tensor index
      CALL create_nd_to_2d_mapping(map, sum_of_arrays(blk_size), map1_2d, map2_2d)
      tensor%nd_index = map
      tensor%blk_sizes = blk_size

      CALL mp_environ_pgrid(dist%pgrid, pdims, task_coor)

#:for ndim in range(1, maxdim+1)
      IF (ndims == ${ndim}$) THEN
         CALL get_arrays(dist%nd_dist, ${varlist("dist", nmax=ndim)}$)
         !CALL get_arrays(dist%blk_sizes, ${varlist("blk_sizes", nmax=ndim)}$)
      ENDIF
#:endfor

      ALLOCATE(nblks_local(ndims))
      ALLOCATE(nfull_local(ndims))
      nfull_local(:) = 0
#:for idim in range(1, maxdim+1)
      IF (ndims .GE. ${idim}$) THEN
         nblks_local(${idim}$) = COUNT(dist_${idim}$ == task_coor(${idim}$))
         ALLOCATE (blks_local_${idim}$ (nblks_local(${idim}$)))
         iblk_count = 0
         DO iblk = 1, SIZE(dist_${idim}$)
            IF (dist_${idim}$(iblk) == task_coor(${idim}$)) THEN
               iblk_count = iblk_count + 1
               blks_local_${idim}$ (iblk_count) = iblk
               nfull_local(${idim}$) = nfull_local(${idim}$) + blk_size_${idim}$ (iblk)
            ENDIF
         ENDDO
      ENDIF
#:endfor

#:for ndim in range(1, maxdim+1)
      IF (ndims == ${ndim}$) THEN
         CALL create_array_list(blks_local, ${ndim}$, ${varlist("blks_local", nmax=ndim)}$)
      ENDIF
#:endfor

      ALLOCATE(tensor%nblks_local(ndims))
      ALLOCATE(tensor%nfull_local(ndims))
      tensor%nblks_local(:) = nblks_local
      tensor%nfull_local(:) = nfull_local

      tensor%blks_local = blks_local

      !tensor%nd_index_grid = dist%nd_index_grid
      tensor%nd_dist = dist%nd_dist
      tensor%pgrid = dist%pgrid
      !tensor%comm_nd = dist%comm_nd
      !tensor%comm_2d = dist%comm_2d

      CALL dbcsr_t_distribution_hold(dist)
      tensor%refcount => dist%refcount

      CALL array_offsets(tensor%blk_sizes, tensor%blk_offsets)

      tensor%valid = .TRUE.
      CALL timestop(handle)
   END SUBROUTINE

   SUBROUTINE dbcsr_t_hold(tensor)
      TYPE(dbcsr_t_type), INTENT(IN) :: tensor
      INTEGER, POINTER :: ref

      IF(tensor%refcount < 1) THEN
         DBCSR_ABORT("can not hold non-existing tensor")
      ENDIF
      ref => tensor%refcount
      ref = ref + 1

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
      ALLOCATE(tensor%matrix_rep)
      IF(.NOT. PRESENT(name)) THEN
         CALL dbcsr_r_create(tensor_in%matrix_rep, tensor%matrix_rep, name=TRIM(tensor_in%name)//" matrix")
      ELSE
         CALL dbcsr_r_create(tensor_in%matrix_rep, tensor%matrix_rep, name=TRIM(name)//" matrix")
      ENDIF
      tensor%owns_matrix = .TRUE.
      CALL dbcsr_r_finalize(tensor%matrix_rep)

      tensor%nd_index_blk = tensor_in%nd_index_blk
      tensor%nd_index = tensor_in%nd_index
      !tensor%nd_index_grid = tensor_in%nd_index_grid
      tensor%blk_sizes = tensor_in%blk_sizes
      tensor%blk_offsets = tensor_in%blk_offsets
      tensor%nd_dist = tensor_in%nd_dist
      tensor%blks_local = tensor_in%blks_local
      ALLOCATE(tensor%nblks_local(ndims_tensor(tensor_in)))
      tensor%nblks_local(:) = tensor_in%nblks_local
      ALLOCATE(tensor%nfull_local(ndims_tensor(tensor_in)))
      tensor%nfull_local(:) = tensor_in%nfull_local
      !tensor%comm_nd = tensor_in%comm_nd
      !tensor%comm_2d = tensor_in%comm_2d
      tensor%pgrid = tensor_in%pgrid

      tensor%refcount => tensor_in%refcount
      CALL dbcsr_t_hold(tensor)

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
      INTEGER                                     :: comm_2d, data_type
      TYPE(dbcsr_distribution_type)                :: matrix_dist
      TYPE(dbcsr_t_distribution_type)             :: dist
      INTEGER, DIMENSION(:), POINTER              :: row_blk_size, col_blk_size
      INTEGER, DIMENSION(:), POINTER              :: col_dist, row_dist
      INTEGER                                   :: handle
      CHARACTER(LEN=*), PARAMETER :: routineN = 'dbcsr_t_create_matrix', &
         routineP = moduleN//':'//routineN
      TYPE(dbcsr_t_pgrid_type)                  :: comm_nd
      INTEGER, DIMENSION(2)                     :: pdims_2d

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
      CALL dbcsr_distribution_get(matrix_dist, group=comm_2d, row_dist=row_dist, col_dist=col_dist, &
                                  nprows=pdims_2d(1), npcols=pdims_2d(2))
      comm_nd = dbcsr_t_nd_mp_comm(comm_2d, [order_in(1)], [order_in(2)], pdims_2d = pdims_2d)

      CALL dbcsr_t_distribution_new( &
         dist, &
         comm_nd, &
         [order_in(1)], [order_in(2)], &
         row_dist, col_dist, own_comm=.TRUE.)

      CALL dbcsr_get_info(matrix_in, &
                          data_type=data_type, &
                          row_blk_size=row_blk_size, &
                          col_blk_size=col_blk_size)

      CALL dbcsr_t_create_new(tensor, name_in, dist, &
                              [order_in(1)], [order_in(2)], &
                              data_type, &
                              row_blk_size, &
                              col_blk_size)

      CALL dbcsr_t_distribution_destroy(dist)
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
      LOGICAL :: abort

      CALL timeset(routineN, handle)
      IF (tensor%owns_matrix) THEN
         CALL dbcsr_r_destroy(tensor%matrix_rep)
         DEALLOCATE (tensor%matrix_rep)
      ELSE
         NULLIFY(tensor%matrix_rep)
      ENDIF
      tensor%owns_matrix = .FALSE.

      CALL destroy_nd_to_2d_mapping(tensor%nd_index_blk)
      CALL destroy_nd_to_2d_mapping(tensor%nd_index)
      !CALL destroy_nd_to_2d_mapping(tensor%nd_index_grid)
      CALL destroy_array_list(tensor%blk_sizes)
      CALL destroy_array_list(tensor%blk_offsets)
      CALL destroy_array_list(tensor%nd_dist)
      CALL destroy_array_list(tensor%blks_local)

      DEALLOCATE(tensor%nblks_local, tensor%nfull_local)

      abort = .FALSE.
      IF (.NOT. ASSOCIATED(tensor%refcount)) THEN
         abort = .TRUE.
      ELSEIF (tensor%refcount < 1) THEN
         abort = .TRUE.
      ENDIF

      IF (abort) THEN
         DBCSR_ABORT("can not destroy non-existing tensor")
      ENDIF

      tensor%refcount = tensor%refcount - 1

      IF(tensor%refcount == 0) THEN
         CALL dbcsr_t_pgrid_destroy(tensor%pgrid)
         !CALL mp_comm_free(tensor%comm_2d)
         !CALL mp_comm_free(tensor%comm_nd)
         DEALLOCATE(tensor%refcount)
      ELSE
         CALL dbcsr_t_pgrid_destroy(tensor%pgrid, keep_comm=.TRUE.)
      ENDIF

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

   SUBROUTINE blk_dims_tensor(tensor, dims)
      TYPE(dbcsr_t_type), INTENT(IN)              :: tensor
      INTEGER, DIMENSION(ndims_tensor(tensor)), &
         INTENT(OUT)                              :: dims

      DBCSR_ASSERT(tensor%valid)
      dims = tensor%nd_index_blk%dims_nd
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

   ! todo: should map1_2d and map2_2d be optional? It is not said that it is consistent with nd grid
   SUBROUTINE dbcsr_t_pgrid_create(mp_comm, dims, pgrid, map1_2d, map2_2d)
      INTEGER, INTENT(IN) :: mp_comm
      INTEGER, DIMENSION(:), INTENT(INOUT) :: dims
      TYPE(dbcsr_t_pgrid_type), INTENT(OUT) :: pgrid
      INTEGER, DIMENSION(:), INTENT(IN), OPTIONAL :: map1_2d, map2_2d
      INTEGER :: nproc, iproc, ndims, i, handle
      INTEGER, DIMENSION(2) :: pdims_2d, pos
      INTEGER, DIMENSION(:), ALLOCATABLE :: map1_2d_prv, map2_2d_prv

      CHARACTER(LEN=*), PARAMETER :: routineN = 'dbcsr_t_pgrid_create', &
         routineP = moduleN//':'//routineN

      CALL timeset(routineN, handle)

      ndims = SIZE(dims)
      IF(PRESENT(map1_2d) .AND. PRESENT(map2_2d)) THEN
         CALL allocate_any(map1_2d_prv, source=map1_2d)
         CALL allocate_any(map2_2d_prv, source=map2_2d)
      ELSE
         ALLOCATE(map1_2d_prv(ndims/2))
         ALLOCATE(map2_2d_prv(ndims-ndims/2))
         map1_2d_prv(:) = (/(i, i=1,SIZE(map1_2d_prv))/)
         map2_2d_prv(:) = (/(i, i=SIZE(map1_2d_prv) + 1, SIZE(map1_2d_prv) + SIZE(map2_2d_prv))/)
      ENDIF

      CALL mp_environ(nproc, iproc, mp_comm)
      IF(ANY(dims == 0)) CALL mp_dims_create(nproc, dims)
      CALL create_nd_to_2d_mapping(pgrid%nd_index_grid, dims, map1_2d_prv, map2_2d_prv, base=0, col_major=.FALSE.)
      CALL get_mapping_info(pgrid%nd_index_grid, dims_2d=pdims_2d)
      CALL mp_cart_create(mp_comm, 2, pdims_2d, pos, pgrid%mp_comm_2d)

      CALL timestop(handle)
   END SUBROUTINE

   SUBROUTINE dbcsr_t_pgrid_destroy(pgrid, keep_comm)
      TYPE(dbcsr_t_pgrid_type), INTENT(INOUT) :: pgrid
      LOGICAL, INTENT(IN), OPTIONAL           :: keep_comm
      LOGICAL :: keep_comm_prv
      IF(PRESENT(keep_comm)) THEN
         keep_comm_prv = keep_comm
      ELSE
         keep_comm_prv = .FALSE.
      ENDIF
      IF(.NOT. keep_comm_prv) CALL mp_comm_free(pgrid%mp_comm_2d)
      CALL destroy_nd_to_2d_mapping(pgrid%nd_index_grid)
   END SUBROUTINE

   SUBROUTINE dbcsr_t_pgrid_remap(pgrid_in, map1_2d, map2_2d, pgrid_out)
      TYPE(dbcsr_t_pgrid_type), INTENT(IN) :: pgrid_in
      INTEGER, DIMENSION(:), INTENT(IN) :: map1_2d, map2_2d
      TYPE(dbcsr_t_pgrid_type), INTENT(OUT) :: pgrid_out
      INTEGER, DIMENSION(:), ALLOCATABLE :: dims

      ALLOCATE(dims(SIZE(map1_2d) + SIZE(map2_2d)))
      CALL get_mapping_info(pgrid_in%nd_index_grid, dims_nd = dims)
      CALL dbcsr_t_pgrid_create(pgrid_in%mp_comm_2d, dims, pgrid_out, map1_2d, map2_2d)
   END SUBROUTINE

   SUBROUTINE mp_environ_pgrid(pgrid, dims, task_coor)
      TYPE(dbcsr_t_pgrid_type), INTENT(IN) :: pgrid
      INTEGER, DIMENSION(ndims_mapping(pgrid%nd_index_grid)), INTENT(OUT) :: dims
      INTEGER, DIMENSION(ndims_mapping(pgrid%nd_index_grid)), INTENT(OUT) :: task_coor
      INTEGER, DIMENSION(2)                                          :: dims_2d, task_coor_2d
      INTEGER :: nproc

      CALL mp_environ(nproc, dims_2d, task_coor_2d, pgrid%mp_comm_2d)
      CALL get_mapping_info(pgrid%nd_index_grid, dims_nd=dims)
      task_coor = get_nd_indices(pgrid%nd_index_grid, INT(task_coor_2d, KIND=int_8))
   END SUBROUTINE

#:for dparam, dtype, dsuffix in dtype_float_list
   SUBROUTINE dbcsr_t_set_${dsuffix}$(tensor, alpha)
      TYPE(dbcsr_t_type), INTENT(INOUT)                   :: tensor
      ${dtype}$, INTENT(IN)                               :: alpha
      CALL dbcsr_r_set(tensor%matrix_rep, alpha)
   END SUBROUTINE
#:endfor

#:for dparam, dtype, dsuffix in dtype_float_list
   SUBROUTINE dbcsr_t_filter_${dsuffix}$(tensor, eps, method, use_absolute)
      TYPE(dbcsr_t_type), INTENT(INOUT)    :: tensor
      ${dtype}$, INTENT(IN)                :: eps
      INTEGER, INTENT(IN), OPTIONAL        :: method
      LOGICAL, INTENT(IN), OPTIONAL        :: use_absolute

      CALL dbcsr_r_filter(tensor%matrix_rep, eps, method, use_absolute)

   END SUBROUTINE
#:endfor

   SUBROUTINE dbcsr_t_get_info(tensor, nblks_total, &
                               nfull_total, &
                               nblks_local, &
                               nfull_local, &
                               pdims, &
                               my_ploc, &
                               ${varlist("blks_local")}$, &
                               ${varlist("proc_dist")}$, &
                               ${varlist("blk_size")}$, &
                               ${varlist("blk_offset")}$, &
                               distribution, &
                               name, &
                               data_type)
      TYPE(dbcsr_t_type), INTENT(IN) :: tensor
      INTEGER, INTENT(OUT), OPTIONAL, DIMENSION(ndims_tensor(tensor)) :: nblks_total, nfull_total, &
                                                                  nblks_local, nfull_local, &
                                                                  my_ploc, pdims
      INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(OUT), OPTIONAL :: ${varlist("blks_local")}$
      INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(OUT), OPTIONAL :: ${varlist("proc_dist")}$
      INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(OUT), OPTIONAL :: ${varlist("blk_size")}$
      INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(OUT), OPTIONAL :: ${varlist("blk_offset")}$
      TYPE(dbcsr_t_distribution_type), INTENT(OUT), OPTIONAL         :: distribution
      CHARACTER(len=*), INTENT(OUT), OPTIONAL                   :: name
      INTEGER, INTENT(OUT), OPTIONAL                            :: data_type
      INTEGER, DIMENSION(ndims_tensor(tensor))                  :: pdims_tmp, my_ploc_tmp

      IF(PRESENT(nblks_total)) CALL get_mapping_info(tensor%nd_index_blk, dims_nd=nblks_total)
      IF(PRESENT(nfull_total)) CALL get_mapping_info(tensor%nd_index, dims_nd=nfull_total)
      IF(PRESENT(nblks_local)) nblks_local(:) = tensor%nblks_local
      IF(PRESENT(nfull_local)) nfull_local(:) = tensor%nfull_local

      IF(PRESENT(my_ploc) .OR. PRESENT(pdims)) CALL mp_environ_pgrid(tensor%pgrid, pdims_tmp, my_ploc_tmp)
      IF(PRESENT(my_ploc)) my_ploc = my_ploc_tmp
      IF(PRESENT(pdims)) pdims = pdims_tmp

#:for idim in range(1, maxdim+1)
      IF (${idim}$ <= ndims_tensor(tensor)) THEN
         IF (PRESENT(blks_local_${idim}$)) CALL get_ith_array(tensor%blks_local, ${idim}$, blks_local_${idim}$)
         IF (PRESENT(proc_dist_${idim}$)) CALL get_ith_array(tensor%nd_dist, ${idim}$, proc_dist_${idim}$)
         IF (PRESENT(blk_size_${idim}$)) CALL get_ith_array(tensor%blk_sizes, ${idim}$, blk_size_${idim}$)
         IF (PRESENT(blk_offset_${idim}$)) CALL get_ith_array(tensor%blk_offsets, ${idim}$, blk_offset_${idim}$)
      ENDIF
#:endfor

      IF(PRESENT(distribution)) distribution = dbcsr_t_distribution(tensor)
      IF(PRESENT(name)) name = tensor%name
      IF(PRESENT(data_type)) data_type = dbcsr_t_get_data_type(tensor)

   END SUBROUTINE

   PURE FUNCTION dbcsr_t_get_num_blocks(tensor) RESULT(num_blocks)
      TYPE(dbcsr_t_type), INTENT(IN)    :: tensor
      INTEGER                           :: num_blocks
      num_blocks = dbcsr_r_get_num_blocks(tensor%matrix_rep)
   END FUNCTION

   FUNCTION dbcsr_t_get_data_size(tensor) RESULT(data_size)
      TYPE(dbcsr_t_type), INTENT(IN)    :: tensor
      INTEGER                           :: data_size
      data_size = dbcsr_r_get_data_size(tensor%matrix_rep)
   END FUNCTION

   SUBROUTINE dbcsr_t_clear(tensor)
      TYPE(dbcsr_t_type), INTENT(INOUT) :: tensor
      CALL dbcsr_r_clear(tensor%matrix_rep)
   END SUBROUTINE

   SUBROUTINE dbcsr_t_finalize(tensor)
      TYPE(dbcsr_t_type), INTENT(INOUT) :: tensor
      CALL dbcsr_r_finalize(tensor%matrix_rep)
   END SUBROUTINE

   SUBROUTINE dbcsr_t_write_tensor_dist(tensor, output_unit)
      TYPE(dbcsr_t_type), INTENT(IN) :: tensor
      INTEGER, INTENT(IN)            :: output_unit
      INTEGER                        :: nproc, myproc, nblock_max, nelement_max
      INTEGER(KIND=int_8)            :: nblock_sum, nelement_sum, nblock_tot
      INTEGER, ALLOCATABLE, DIMENSION(:) :: nblock, nelement
      INTEGER                        :: mp_comm
      INTEGER, DIMENSION(2)          :: tmp
      INTEGER, DIMENSION(ndims_tensor(tensor)) :: bdims
      REAL(KIND=real_8)              :: occupation

      mp_comm = tensor%pgrid%mp_comm_2d
      CALL mp_environ(nproc, myproc, mp_comm)

      ALLOCATE(nblock(0:nproc-1))
      nblock(:) = 0

      ALLOCATE(nelement(0:nproc-1))
      nelement(:) = 0

      nblock(myproc) = dbcsr_t_get_num_blocks(tensor)
      nelement(myproc) = dbcsr_t_get_data_size(tensor)

      nblock_sum = nblock(myproc)
      nblock_max = nblock(myproc)
      nelement_sum = nelement(myproc)
      nelement_max = nelement(myproc)
      CALL mp_sum(nblock_sum, mp_comm)
      CALL mp_sum(nelement_sum, mp_comm)
      tmp = (/nblock_max, nelement_max/)
      CALL mp_max(tmp, mp_comm)
      nblock_max = tmp(1); nelement_max = tmp(2)

      CALL blk_dims_tensor(tensor, bdims)
      nblock_tot = PRODUCT(INT(bdims, KIND=int_8))

      occupation = -1.0_real_8
      IF(nblock_tot .NE. 0) occupation = 100.0_real_8*REAL(nblock_sum, real_8)/REAL(nblock_tot, real_8)

      IF(output_unit > 0) THEN
         WRITE(output_unit, "(T2,A)") &
            "DISTRIBUTION OF "//TRIM(tensor%name)
         WRITE(output_unit, "(T15,A,T68,I13)") "Number of non-zero blocks:", nblock_sum
         WRITE(output_unit, "(T15,A,T75,F6.2)") "Percentage of non-zero blocks:", occupation
         WRITE(output_unit, "(T15,A,T68,I13)") "Average number of blocks per CPU:", (nblock_sum+nproc-1)/nproc
         WRITE(output_unit, "(T15,A,T68,I13)") "Maximum number of blocks per CPU:", nblock_max
         WRITE(output_unit, "(T15,A,T68,I13)") "Average number of matrix elements per CPU:", (nelement_sum+nproc-1)/nproc
         WRITE(output_unit, "(T15,A,T68,I13)") "Maximum number of matrix elements per CPU:", nelement_max
      ENDIF

   END SUBROUTINE

   SUBROUTINE dbcsr_t_write_tensor_info(tensor, output_unit, full_info)
      TYPE(dbcsr_t_type), INTENT(IN) :: tensor
      INTEGER, INTENT(IN)            :: output_unit
      LOGICAL, OPTIONAL, INTENT(IN)  :: full_info
      INTEGER, DIMENSION(ndims_tensor(tensor)) :: nblks_total, nfull_total, pdims, my_ploc, nblks_local, nfull_local
      INTEGER, DIMENSION(:), ALLOCATABLE :: ${varlist("blks_local")}$
      INTEGER, DIMENSION(:), ALLOCATABLE :: ${varlist("proc_dist")}$
      INTEGER, DIMENSION(:), ALLOCATABLE :: ${varlist("blk_size")}$
      INTEGER, DIMENSION(:), ALLOCATABLE :: ${varlist("blk_offset")}$
      CHARACTER(len=default_string_length)                   :: name
      INTEGER                            :: idim
      INTEGER                            :: iblk

      CALL dbcsr_t_get_info(tensor, nblks_total, nfull_total, nblks_local, nfull_local, pdims, my_ploc, &
                            ${varlist("blks_local")}$, ${varlist("proc_dist")}$, ${varlist("blk_size")}$, &
                            ${varlist("blk_offset")}$, &
                            name=name)

      IF(output_unit > 0) THEN
         WRITE(output_unit, "(T2,A)") &
            "GLOBAL INFO OF "//TRIM(name)
         WRITE(output_unit, "(T4,A,1X)", advance="no") "block dimensions:"
         DO idim = 1, ndims_tensor(tensor)
            WRITE(output_unit, "(I6)", advance="no") nblks_total(idim)
         ENDDO
         WRITE(output_unit, "(/T4,A,1X)", advance="no") "full dimensions:"
         DO idim = 1, ndims_tensor(tensor)
            WRITE(output_unit, "(I8)", advance="no") nfull_total(idim)
         ENDDO
         WRITE(output_unit, "(/T4,A,1X)", advance="no") "process grid dimensions:"
         DO idim = 1, ndims_tensor(tensor)
            WRITE(output_unit, "(I6)", advance="no") pdims(idim)
         ENDDO

         IF(PRESENT(full_info)) THEN
            IF(full_info) THEN
               WRITE(output_unit,'(/T4,A)', advance='no') "Block sizes:"
#:for dim in range(1, maxdim+1)
               IF(ndims_tensor(tensor)>=${dim}$) THEN
                  WRITE(output_unit,'(/T8,A,1X,I1,A,1X)', advance='no') 'Dim', ${dim}$,':'
                  DO iblk=1, SIZE(blk_size_${dim}$)
                     WRITE(output_unit,'(I2,1X)', advance='no') blk_size_${dim}$(iblk)
                  ENDDO
               ENDIF
#:endfor
               WRITE(output_unit,'(/T4,A)', advance='no') "Block distribution:"
#:for dim in range(1, maxdim+1)
               IF(ndims_tensor(tensor)>=${dim}$) THEN
                  WRITE(output_unit,'(/T8,A,1X,I1,A,1X)', advance='no') 'Dim', ${dim}$,':'
                  DO iblk=1, SIZE(proc_dist_${dim}$)
                     WRITE(output_unit,'(I3,1X)', advance='no') proc_dist_${dim}$(iblk)
                  ENDDO
               ENDIF
#:endfor
            ENDIF
            WRITE (output_unit, *)
         ENDIF
      ENDIF

   END SUBROUTINE
END MODULE
