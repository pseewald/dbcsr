!--------------------------------------------------------------------------------------------------!
! Copyright (C) by the DBCSR developers group - All rights reserved                                !
! This file is part of the DBCSR library.                                                          !
!                                                                                                  !
! For information on the license, see the LICENSE file.                                            !
! For further information please visit https://dbcsr.cp2k.org                                      !
! SPDX-License-Identifier: GPL-2.0+                                                                !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief DBCSR tensor framework for block-sparse tensor contraction.
!>
!> Representation of n-rank tensors as DBCSR matrices.
!> Support for arbitrary redistribution between different representations.
!> Support for arbitrary tensor contractions
!> DBCSR routines are generalized to n dimensions by light-weight wrapper routines.
!> \todo implement checks and error messages
!> \author Patrick Seewald
! **************************************************************************************************
MODULE dbcsr_tensor

#:include "dbcsr_tensor.fypp"
#:set maxdim = maxrank
#:set ndims = range(2,maxdim+1)

   USE dbcsr_allocate_wrap,             ONLY: allocate_any
   USE dbcsr_array_list_methods,        ONLY: get_arrays,&
                                              reorder_arrays,&
                                              get_ith_array,&
                                              array_list, &
                                              array_sublist
   USE dbcsr_api,                       ONLY: &
        dbcsr_type, dbcsr_iterator_type, dbcsr_iterator_blocks_left, &
        dbcsr_iterator_next_block, dbcsr_iterator_start, dbcsr_iterator_stop, &
        dbcsr_transpose, dbcsr_no_transpose, dbcsr_scalar, &
        ${uselist(dtype_float_param)}$
   USE dbcsr_rectangular_base, ONLY:&
        dbcsr_r_copy, dbcsr_r_create, dbcsr_r_finalize, dbcsr_r_get_data_type, dbcsr_r_get_info, dbcsr_r_info
   USE dbcsr_rectangular_mm, ONLY: dbcsr_r_multiply
   USE dbcsr_tensor_block,              ONLY: dbcsr_t_iterator_type, &
                                              dbcsr_t_get_block, &
                                              dbcsr_t_put_block, &
                                              dbcsr_t_iterator_start, &
                                              dbcsr_t_iterator_blocks_left, &
                                              dbcsr_t_iterator_stop, &
                                              dbcsr_t_iterator_next_block, &
                                              ndims_iterator,&
                                              dbcsr_t_reserve_blocks, &
                                              block_nd
   USE dbcsr_tensor_index,              ONLY: get_mapping_info,&
                                              nd_to_2d_mapping,&
                                              dbcsr_t_inverse_order,&
                                              permute_index
   USE dbcsr_tensor_types,              ONLY: dbcsr_t_create,&
                                              dbcsr_t_get_data_type,&
                                              dbcsr_t_type,&
                                              ndims_tensor,&
                                              dims_tensor,&
                                              dbcsr_t_distribution_type,&
                                              dbcsr_t_distribution,&
                                              dbcsr_t_nd_mp_comm,&
                                              dbcsr_t_destroy,&
                                              dbcsr_t_distribution_destroy,&
                                              dbcsr_t_distribution_new,&
                                              dbcsr_t_get_stored_coordinates,&
                                              blk_dims_tensor,&
                                              dbcsr_t_hold
   USE dbcsr_kinds,                     ONLY: ${uselist(dtype_float_prec)}$,&
                                              default_string_length
   USE dbcsr_mpiwrap,                   ONLY: mp_environ
   USE dbcsr_toollib,                   ONLY: sort
   USE dbcsr_tensor_reshape,            ONLY: dbcsr_t_reshape

#include "base/dbcsr_base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbcsr_tensor'


   PUBLIC :: &
      dbcsr_t_contract, &
      dbcsr_t_copy, &
      dbcsr_t_dims, &
      dbcsr_t_get_block, &
      dbcsr_t_get_stored_coordinates, &
      dbcsr_t_inverse_order, &
      dbcsr_t_iterator_blocks_left, &
      dbcsr_t_iterator_next_block, &
      dbcsr_t_iterator_start, &
      dbcsr_t_iterator_stop, &
      dbcsr_t_iterator_type, &
      dbcsr_t_ndims, &
      dbcsr_t_put_block, &
      dbcsr_t_reserve_blocks, &
      dbcsr_t_split_blocks

   INTERFACE dbcsr_t_ndims
      MODULE PROCEDURE ndims_tensor
   END INTERFACE

   INTERFACE dbcsr_t_dims
      MODULE PROCEDURE dims_tensor
   END INTERFACE

CONTAINS

! **************************************************************************************************
!> \brief Copy tensor data.
!>        Redistributes tensor data according to distributions of target and source tensor.
!>        Permutes tensor index according to `order` argument (if present).
!>
!> Source and target tensor formats are arbitrary as long as the following requirements are met:
!> * source and target tensors have the same number of blocks in each dimension and the same block sizes.
!>   If `order` argument is present, this must be the case after index permutation.
!> OR
!> * target tensor is not yet created, in this case an exact copy of source tensor is returned.
!> \param tensor_in Source
!> \param tensor_out Target
!> \param order Permutation of target tensor index. Exact same convention as order argument of RESHAPE
!>              intrinsic
! **************************************************************************************************
   SUBROUTINE dbcsr_t_copy(tensor_in, tensor_out, order)
      TYPE(dbcsr_t_type), INTENT(INOUT)              :: tensor_in, tensor_out
      INTEGER, DIMENSION(ndims_tensor(tensor_in)), &
         INTENT(IN), OPTIONAL                        :: order

      TYPE(dbcsr_t_type)                             :: tensor_tmp
      INTEGER                                        :: handle

      CHARACTER(LEN=*), PARAMETER :: routineN = 'dbcsr_t_copy', &
         routineP = moduleN//':'//routineN

      CALL timeset(routineN, handle)

      IF (PRESENT(order)) THEN
         CALL dbcsr_t_permute_index(tensor_in, tensor_tmp, order)
         CALL dbcsr_t_reshape(tensor_tmp, tensor_out)
         CALL dbcsr_t_destroy(tensor_tmp)
      ELSE
         CALL dbcsr_t_reshape(tensor_in, tensor_out)
      ENDIF

      CALL timestop(handle)

   END SUBROUTINE

! **************************************************************************************************
!> \brief copy matrix to tensor.
!> \param matrix_in ...
!> \param tensor_out ...
! **************************************************************************************************
   SUBROUTINE dbcsr_t_copy_matrix(matrix_in, tensor_out)
      TYPE(dbcsr_type), INTENT(IN)                        :: matrix_in
      TYPE(dbcsr_t_type), INTENT(INOUT)                  :: tensor_out

      INTEGER, DIMENSION(2)                              :: ind_2d
      REAL(KIND=real_8), ALLOCATABLE, DIMENSION(:, :)    :: block_arr
      REAL(KIND=real_8), DIMENSION(:, :), POINTER        :: block
      TYPE(dbcsr_iterator_type)                          :: iter

      INTEGER                                            :: handle
      CHARACTER(LEN=*), PARAMETER :: routineN = 'dbcsr_t_copy_matrix', &
         routineP = moduleN//':'//routineN

      CALL timeset(routineN, handle)
      DBCSR_ASSERT(tensor_out%valid)

      CALL dbcsr_t_reserve_blocks(matrix_in, tensor_out)

      CALL dbcsr_iterator_start(iter, matrix_in)
      DO WHILE (dbcsr_iterator_blocks_left(iter))
         CALL dbcsr_iterator_next_block(iter, ind_2d(1), ind_2d(2), block)
         CALL allocate_any(block_arr, source=block)
         CALL dbcsr_t_put_block(tensor_out, ind_2d, SHAPE(block_arr), block_arr)
         DEALLOCATE (block_arr)
      ENDDO
      CALL dbcsr_iterator_stop(iter)

      CALL timestop(handle)

   END SUBROUTINE

! **************************************************************************************************
!> \brief Contract tensors by multiplying matrix representations.
!> \param tensor_1 first tensor (in)
!> \param tensor_2 second tensor (in)
!> \param contract_1 indices of tensor_1 to contract
!> \param notcontract_1 indices of tensor_1 not to contract
!> \param contract_2 indices of tensor_2 to contract (1:1 with contract_1)
!> \param notcontract_2 indices of tensor_2 not to contract
!> \param map_1 which indices of tensor_3 map to non-contracted indices of tensor_1 (1:1 with notcontract_1)
!> \param map_2 which indices of tensor_3 map to non-contracted indices of tensor_2 (1:1 with notcontract_2)
!> \param tensor_3 contracted tensor (out)
!> \todo Generate input arguments (contract_*, notcontract_*, map_*) from Einstein notation
! **************************************************************************************************
   SUBROUTINE dbcsr_t_contract(tensor_1, tensor_2, tensor_3, &
                               contract_1, notcontract_1, &
                               contract_2, notcontract_2, &
                               map_1, map_2)
      TYPE(dbcsr_t_type), INTENT(INOUT)                 :: tensor_1, tensor_2
      INTEGER, DIMENSION(:), INTENT(IN)              :: contract_1, contract_2, map_1, map_2
      INTEGER, DIMENSION(:), INTENT(IN)              :: notcontract_1, notcontract_2
      TYPE(dbcsr_t_type), INTENT(INOUT)              :: tensor_3
      TYPE(dbcsr_t_type)                             :: tensor_contr_1, tensor_contr_2, tensor_contr_3, &
                                                        tensor_tmp_1, tensor_tmp_2, tensor_tmp_3
      LOGICAL                                        :: assert_stmt
      INTEGER                                        :: data_type, handle, max_mm_dim, max_tensor
      !INTEGER, DIMENSION(:), ALLOCATABLE             :: sort_indices
      INTEGER, DIMENSION(SIZE(contract_1))           :: contract_1_mod
      INTEGER, DIMENSION(SIZE(notcontract_1))        :: notcontract_1_mod
      INTEGER, DIMENSION(SIZE(contract_2))           :: contract_2_mod
      INTEGER, DIMENSION(SIZE(notcontract_2))        :: notcontract_2_mod
      INTEGER, DIMENSION(SIZE(map_1))                :: map_1_mod
      INTEGER, DIMENSION(SIZE(map_2))                :: map_2_mod
      CHARACTER(LEN=1)                               :: trans_1, trans_2, trans_3
      LOGICAL                                        :: new_3
      INTEGER, DIMENSION(:), ALLOCATABLE             :: dims1, dims2, dims3

      CHARACTER(LEN=*), PARAMETER :: routineN = 'dbcsr_t_contract', &
         routineP = moduleN//':'//routineN

      CALL timeset(routineN, handle)

      DBCSR_ASSERT(tensor_1%valid)
      DBCSR_ASSERT(tensor_2%valid)
      DBCSR_ASSERT(tensor_3%valid)

      assert_stmt = SIZE(contract_1) .EQ. SIZE(contract_2)
      DBCSR_ASSERT(assert_stmt)

      assert_stmt = SIZE(map_1) .EQ. SIZE(notcontract_1)
      DBCSR_ASSERT(assert_stmt)

      assert_stmt = SIZE(map_2) .EQ. SIZE(notcontract_2)
      DBCSR_ASSERT(assert_stmt)

      assert_stmt = SIZE(notcontract_1) + SIZE(contract_1) .EQ. dbcsr_t_ndims(tensor_1)
      DBCSR_ASSERT(assert_stmt)

      assert_stmt = SIZE(notcontract_2) + SIZE(contract_2) .EQ. dbcsr_t_ndims(tensor_2)
      DBCSR_ASSERT(assert_stmt)

      assert_stmt = SIZE(map_1) + SIZE(map_2) .EQ. dbcsr_t_ndims(tensor_3)
      DBCSR_ASSERT(assert_stmt)

      assert_stmt = dbcsr_t_get_data_type(tensor_1) .EQ. dbcsr_t_get_data_type(tensor_2)
      DBCSR_ASSERT(assert_stmt)

      data_type = dbcsr_t_get_data_type(tensor_1)

      ! align tensor index with data, tensor data is not modified

      CALL align_tensor(tensor_1, contract_1, notcontract_1, &
                        tensor_tmp_1, contract_1_mod, notcontract_1_mod)

      CALL align_tensor(tensor_2, contract_2, notcontract_2, &
                        tensor_tmp_2, contract_2_mod, notcontract_2_mod)

      CALL align_tensor(tensor_3, map_1, map_2, &
                        tensor_tmp_3, map_1_mod, map_2_mod)


      ALLOCATE(dims1(dbcsr_t_ndims(tensor_1)))
      ALLOCATE(dims2(dbcsr_t_ndims(tensor_2)))
      ALLOCATE(dims3(dbcsr_t_ndims(tensor_3)))

      CALL blk_dims_tensor(tensor_1, dims1)
      CALL blk_dims_tensor(tensor_2, dims2)
      CALL blk_dims_tensor(tensor_3, dims3)

      max_mm_dim = MAXLOC([PRODUCT(dims1(notcontract_1)), PRODUCT(dims1(contract_1)), PRODUCT(dims2(notcontract_2))], DIM=1)
      max_tensor = MAXLOC([PRODUCT(dims1), PRODUCT(dims2), PRODUCT(dims3)], DIM=1)


      SELECT CASE(max_mm_dim)
      CASE(1)

         CALL index_linked_sort(contract_1_mod, contract_2_mod)
         CALL index_linked_sort(map_2_mod, notcontract_2_mod)
         SELECT CASE(max_tensor)
         CASE(1)
            CALL index_linked_sort(notcontract_1_mod, map_1_mod)
         CASE(3)
            CALL index_linked_sort(map_1_mod, notcontract_1_mod)
         CASE DEFAULT
            DBCSR_ABORT("should not happen")
         END SELECT

         CALL reshape_mm_compatible(tensor_tmp_1, tensor_tmp_3, tensor_contr_1, tensor_contr_3, &
                                    contract_1_mod, notcontract_1_mod, map_2_mod, map_1_mod, &
                                    trans_1, trans_3, nodata2=.TRUE., new2=new_3)

         CALL reshape_mm_small(tensor_tmp_2, contract_2_mod, notcontract_2_mod, tensor_contr_2, trans_2)

      CASE(2)

         CALL index_linked_sort(notcontract_1_mod, map_1_mod)
         CALL index_linked_sort(notcontract_2_mod, map_2_mod)
         SELECT CASE(max_tensor)
         CASE(1)
            CALL index_linked_sort(contract_1_mod, contract_2_mod)
         CASE(2)
            CALL index_linked_sort(contract_2_mod, contract_1_mod)
         CASE DEFAULT
            DBCSR_ABORT("should not happen")
         END SELECT

         CALL reshape_mm_compatible(tensor_tmp_1, tensor_tmp_2, tensor_contr_1, tensor_contr_2, &
                                    notcontract_1_mod, contract_1_mod, notcontract_2_mod, contract_2_mod, &
                                    trans_1, trans_2)
         CALL invert_transpose_flag(trans_1)

         CALL reshape_mm_small(tensor_tmp_3, map_1_mod, map_2_mod, tensor_contr_3, trans_3, nodata=.TRUE., new=new_3)

      CASE(3)

         CALL index_linked_sort(map_1_mod, notcontract_1_mod)
         CALL index_linked_sort(contract_2_mod, contract_1_mod)
         SELECT CASE(max_tensor)
         CASE(2)
            CALL index_linked_sort(notcontract_2_mod, map_2_mod)
         CASE(3)
            CALL index_linked_sort(map_2_mod, notcontract_2_mod)
         CASE DEFAULT
            DBCSR_ABORT("should not happen")
         END SELECT

         CALL reshape_mm_compatible(tensor_tmp_2, tensor_tmp_3, tensor_contr_2, tensor_contr_3, &
                                    contract_2_mod, notcontract_2_mod, map_1_mod, map_2_mod, &
                                    trans_2, trans_3, nodata2=.TRUE., new2=new_3)

         CALL invert_transpose_flag(trans_2)
         CALL invert_transpose_flag(trans_3)

         CALL reshape_mm_small(tensor_tmp_1, notcontract_1_mod, contract_1_mod, tensor_contr_1, trans_1)

      END SELECT

      CALL dbcsr_t_destroy(tensor_tmp_1)
      CALL dbcsr_t_destroy(tensor_tmp_2)

      CALL dbcsr_r_multiply(trans_1, trans_2, trans_3, dbcsr_scalar(1.0_real_8), &
                            tensor_contr_1%matrix_rep, tensor_contr_2%matrix_rep, &
                            tensor_contr_3%matrix_rep)


      IF(new_3) THEN
         ! need redistribute if we created new tensor for tensor 3
         CALL dbcsr_t_copy(tensor_contr_3, tensor_tmp_3)
         CALL dbcsr_r_copy(tensor_3%matrix_rep, tensor_tmp_3%matrix_rep, shallow_data=.TRUE.)
      ELSE
         ! shallow copy sufficient since no new matrix has been created for tensor 3
         CALL dbcsr_r_copy(tensor_3%matrix_rep, tensor_contr_3%matrix_rep, shallow_data=.TRUE.)
      ENDIF

      CALL dbcsr_t_destroy(tensor_tmp_3)

      CALL dbcsr_t_destroy(tensor_contr_1)
      CALL dbcsr_t_destroy(tensor_contr_2)
      CALL dbcsr_t_destroy(tensor_contr_3)

      CALL timestop(handle)

   CONTAINS

! **************************************************************************************************
!> \brief align tensor index with data
!> \param tensor_in ...
!> \param contract_in ...
!> \param notcontract_in ...
!> \param tensor_out ...
!> \param contract_out ...
!> \param notcontract_out ...
! **************************************************************************************************
      SUBROUTINE align_tensor(tensor_in, contract_in, notcontract_in, &
                              tensor_out, contract_out, notcontract_out)
         TYPE(dbcsr_t_type), INTENT(INOUT)               :: tensor_in
         INTEGER, DIMENSION(:), INTENT(IN)            :: contract_in, notcontract_in
         TYPE(dbcsr_t_type), INTENT(OUT)              :: tensor_out
         INTEGER, DIMENSION(SIZE(contract_in)), &
            INTENT(OUT)                               :: contract_out
         INTEGER, DIMENSION(SIZE(notcontract_in)), &
            INTENT(OUT)                               :: notcontract_out
         INTEGER, DIMENSION(dbcsr_t_ndims(tensor_in)) :: align

         CALL dbcsr_t_align_index(tensor_in, tensor_out, order=align)
         contract_out = align(contract_in)
         notcontract_out = align(notcontract_in)

      END SUBROUTINE

! **************************************************************************************************
!> \brief Prepare tensor for contraction: redistribute to a 2d format which can be contracted by
!>        matrix multiplication.
!> \param tensor_in ...
!> \param contract ...
!> \param notcontract ...
!> \param contract_index ...
!> \param tensor_out ...
!> \param trans1 'N' if tensor 1 out has format (ind linked, ind free), 'T' otherwise
!> \param trans2 'N' if tensor 2 out has format (ind linked, ind free), 'T' otherwise
! **************************************************************************************************
      SUBROUTINE reshape_mm_compatible(tensor1, tensor2, tensor1_out, tensor2_out, ind1_free, ind1_linked, &
                                       ind2_free, ind2_linked, trans1, trans2, nodata1, nodata2, new1, new2)
         TYPE(dbcsr_t_type), INTENT(INOUT)           :: tensor1
         TYPE(dbcsr_t_type), INTENT(INOUT)           :: tensor2
         TYPE(dbcsr_t_type), INTENT(OUT)             :: tensor1_out, tensor2_out
         INTEGER, DIMENSION(:), INTENT(IN)           :: ind1_free, ind2_free
         INTEGER, DIMENSION(:), INTENT(IN)           :: ind1_linked, ind2_linked
         CHARACTER(LEN=1), INTENT(OUT)               :: trans1, trans2
         LOGICAL, INTENT(IN), OPTIONAL               :: nodata1, nodata2
         LOGICAL, INTENT(OUT), OPTIONAL              :: new1, new2
         INTEGER                                     :: ref_tensor, compat1, compat2
         TYPE(array_list)                            :: dist_list
         INTEGER, DIMENSION(:), ALLOCATABLE                        :: mp_dims
         TYPE(dbcsr_t_distribution_type)             :: dist_in

         IF (SIZE(ind1_free) .GE. SIZE(ind2_free)) THEN
            ref_tensor = 1
         ELSE
            ref_tensor = 2
         ENDIF

         compat1 = compat_map(tensor1%nd_index, ind1_linked)
         compat2 = compat_map(tensor2%nd_index, ind2_linked)

         IF(PRESENT(new1)) new1 = .FALSE.
         IF(PRESENT(new2)) new2 = .FALSE.

         IF(compat1 == 0) THEN
            IF(PRESENT(new1)) new1 = .TRUE.
         ENDIF

         IF(compat2 == 0) THEN
            IF(PRESENT(new2)) new2 = .TRUE.
         ENDIF

         IF (ref_tensor == 1) THEN ! tensor 1 is reference and tensor 2 is reshaped compatible with tensor 1
            IF (compat1 == 0) THEN ! tensor 1 is not contraction compatible --> reshape
               CALL dbcsr_t_remap(tensor1, ind1_linked, ind1_free, tensor1_out, nodata=nodata1)
               compat1 = 1
            ELSE
               CALL dbcsr_t_copy(tensor1, tensor1_out)
            ENDIF
            IF (compat2 == 0) THEN ! tensor 2 is not contraction compatible --> reshape
               dist_in = dbcsr_t_distribution(tensor1_out)
               dist_list = array_sublist(dist_in%nd_dist, ind1_linked)
               IF (compat1 == 1) THEN ! linked index is first 2d dimension
                  ! get distribution of linked index, tensor 2 must adopt this distribution
                  ! get grid dimensions of linked index
                  CALL get_mapping_info(dist_in%nd_index_grid, dims1_2d=mp_dims)
                  CALL dbcsr_t_remap(tensor2, ind2_linked, ind2_free, tensor2_out, comm_2d=dist_in%comm_2d, &
                                     dist1=dist_list, mp_dims_1=mp_dims, nodata=nodata2)
               ELSEIF (compat1 == 2) THEN ! linked index is second 2d dimension
                  ! get distribution of linked index, tensor 2 must adopt this distribution
                  ! get grid dimensions of linked index
                  CALL get_mapping_info(dist_in%nd_index_grid, dims2_2d=mp_dims)
                  CALL dbcsr_t_remap(tensor2, ind2_free, ind2_linked, tensor2_out, comm_2d=dist_in%comm_2d, &
                                     dist2=dist_list, mp_dims_2=mp_dims, nodata=nodata2)
               ELSE
                  DBCSR_ABORT("should not happen")
               ENDIF
            ELSE
               CALL dbcsr_t_copy(tensor2, tensor2_out)
            ENDIF
         ELSE ! tensor 2 is reference and tensor 1 is reshaped compatible with tensor 2
            IF (compat2 == 0) THEN ! tensor 2 is not contraction compatible --> reshape
               CALL dbcsr_t_remap(tensor2, ind2_linked, ind2_free, tensor2_out, nodata=nodata2)
               compat2 = 1
            ELSE
               CALL dbcsr_t_copy(tensor2, tensor2_out)
            ENDIF
            IF (compat1 == 0) THEN ! tensor 1 is not contraction compatible --> reshape
               dist_in = dbcsr_t_distribution(tensor2_out)
               dist_list = array_sublist(dist_in%nd_dist, ind2_linked)
               IF (compat2 == 1) THEN
                  CALL get_mapping_info(dist_in%nd_index_grid, dims1_2d=mp_dims)
                  CALL dbcsr_t_remap(tensor1, ind1_linked, ind1_free, tensor1_out, comm_2d=dist_in%comm_2d, &
                                     dist1=dist_list, mp_dims_1=mp_dims, nodata=nodata1)
               ELSEIF (compat1 == 2) THEN
                  CALL get_mapping_info(dist_in%nd_index_grid, dims2_2d=mp_dims)
                  CALL dbcsr_t_remap(tensor1, ind1_free, ind1_linked, tensor1_out, comm_2d=dist_in%comm_2d, &
                                     dist2=dist_list, mp_dims_2=mp_dims, nodata=nodata1)
               ELSE
                  DBCSR_ABORT("should not happen")
               ENDIF
            ELSE
               CALL dbcsr_t_copy(tensor1, tensor1_out)
            ENDIF
         ENDIF

         SELECT CASE(compat1)
         CASE(1)
            trans1 = dbcsr_no_transpose
         CASE(2)
            trans1 = dbcsr_transpose
         CASE DEFAULT
            DBCSR_ABORT("should not happen")
         END SELECT

         SELECT CASE(compat2)
         CASE(1)
            trans2 = dbcsr_no_transpose
         CASE(2)
            trans2 = dbcsr_transpose
         CASE DEFAULT
            DBCSR_ABORT("should not happen")
         END SELECT


      END SUBROUTINE

      SUBROUTINE reshape_mm_small(tensor_in, ind1, ind2, tensor_out, trans, nodata, new)
         ! was called prep_tensor
         TYPE(dbcsr_t_type), INTENT(INOUT)           :: tensor_in
         INTEGER, DIMENSION(:), INTENT(IN)           :: ind1, ind2
         TYPE(dbcsr_t_type), INTENT(OUT)             :: tensor_out
         CHARACTER(LEN=1), INTENT(OUT)               :: trans
         LOGICAL, INTENT(IN), OPTIONAL               :: nodata
         LOGICAL, INTENT(OUT), OPTIONAL              :: new
         INTEGER                                     :: compat1, compat2

         IF(PRESENT(new)) new = .FALSE.
         compat1 = compat_map(tensor_in%nd_index, ind1)
         compat2 = compat_map(tensor_in%nd_index, ind2)
         IF (compat1 == 0 .or. compat2 == 0) THEN ! index mapping not compatible with contract index
            CALL dbcsr_t_remap(tensor_in, ind1, ind2, tensor_out, nodata=nodata)
            compat1 = 1
            compat2 = 2
            IF(PRESENT(new)) new = .TRUE.
         ELSE
            CALL dbcsr_t_copy(tensor_in, tensor_out)
         ENDIF

         IF(compat1 == 1 .AND. compat2 == 2) THEN
            trans = dbcsr_no_transpose
         ELSEIF (compat1 == 2 .AND. compat2 == 1) THEN
            trans = dbcsr_transpose
         ELSE
            DBCSR_ABORT("this should not happen")
         ENDIF

      END SUBROUTINE


! **************************************************************************************************
!> \brief Check if 2d index is compatible with tensor index
!> \param nd_index ...
!> \param compat_ind ...
! **************************************************************************************************
      FUNCTION compat_map(nd_index, compat_ind)
         TYPE(nd_to_2d_mapping), INTENT(IN) :: nd_index
         INTEGER, DIMENSION(:), INTENT(IN)  :: compat_ind
         INTEGER, DIMENSION(:), ALLOCATABLE :: map1, map2
         INTEGER                            :: compat_map

         CALL get_mapping_info(nd_index, map1_2d=map1, map2_2d=map2)

         compat_map = 0
         IF(array_eq_i(map1, compat_ind)) THEN
            compat_map = 1
         ELSEIF(array_eq_i(map2, compat_ind)) THEN
            compat_map = 2
         ENDIF

      END FUNCTION

! **************************************************************************************************
!> \brief Check if 2 arrays are equal
!> \param arr1 ...
!> \param arr2
! **************************************************************************************************
      PURE FUNCTION array_eq_i(arr1, arr2)
         INTEGER, INTENT(IN), DIMENSION(:) :: arr1
         INTEGER, INTENT(IN), DIMENSION(:) :: arr2
         LOGICAL                           :: array_eq_i

         array_eq_i = .FALSE.
         IF (SIZE(arr1) .EQ. SIZE(arr2)) array_eq_i = ALL(arr1 == arr2)

      END FUNCTION

      SUBROUTINE invert_transpose_flag(trans_flag)
         CHARACTER(LEN=1), INTENT(INOUT)                    :: trans_flag

         IF (trans_flag == dbcsr_transpose) THEN
            trans_flag = dbcsr_no_transpose
         ELSEIF (trans_flag == dbcsr_no_transpose) THEN
            trans_flag = dbcsr_transpose
         ENDIF
      END SUBROUTINE

      SUBROUTINE index_linked_sort(ind_ref, ind_dep)
         INTEGER, DIMENSION(:), INTENT(INOUT) :: ind_ref, ind_dep
         INTEGER, DIMENSION(SIZE(ind_ref))    :: sort_indices

         CALL sort(ind_ref, SIZE(ind_ref), sort_indices)
         ind_dep(:) = ind_dep(sort_indices)

      END SUBROUTINE

   END SUBROUTINE

! **************************************************************************************************
!> \brief set up cyclic distribution. This is used for all internally created temporary tensors.
!> \param dist_array ...
!> \param dist_size ...
!> \param nbins ...
! **************************************************************************************************
   SUBROUTINE cyclic_dist(dist_array, dist_size, nbins)
      INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: dist_array
      INTEGER, INTENT(in)                             :: dist_size, nbins
      INTEGER                                         :: i

      CALL allocate_any(dist_array, source=(/(MODULO(nbins-i, nbins), i=1, dist_size)/))

   END SUBROUTINE cyclic_dist

! **************************************************************************************************
!> \brief Copy tensor to tensor with modified index mapping
!> \param tensor_in ...
!> \param map1_2d new index mapping
!> \param map2_2d new index mapping
!> \param tensor_out ...
!> \param name ...
! **************************************************************************************************
   SUBROUTINE dbcsr_t_remap(tensor_in, map1_2d, map2_2d, tensor_out, comm_2d, dist1, dist2, &
                            mp_dims_1, mp_dims_2, name, nodata)
      TYPE(dbcsr_t_type), INTENT(INOUT)      :: tensor_in
      INTEGER, DIMENSION(:), INTENT(IN)      :: map1_2d, map2_2d
      TYPE(dbcsr_t_type), INTENT(OUT)        :: tensor_out
      CHARACTER(len=*), INTENT(IN), OPTIONAL :: name
      LOGICAL, INTENT(IN), OPTIONAL          :: nodata
      INTEGER, INTENT(IN), OPTIONAL          :: comm_2d
      TYPE(array_list), INTENT(IN), OPTIONAL :: dist1, dist2
      INTEGER, DIMENSION(SIZE(map1_2d)), OPTIONAL :: mp_dims_1
      INTEGER, DIMENSION(SIZE(map2_2d)), OPTIONAL :: mp_dims_2
      CHARACTER(len=default_string_length)   :: name_tmp
      INTEGER, DIMENSION(:), ALLOCATABLE     :: ${varlist("blk_sizes")}$,&
                                                ${varlist("nd_dist")}$
      TYPE(dbcsr_t_distribution_type)        :: dist
      INTEGER                                :: comm_2d_prv, comm_nd, handle, i
      INTEGER, DIMENSION(dbcsr_t_ndims(tensor_in)) :: pdims, myploc
      LOGICAL, DIMENSION(dbcsr_t_ndims(tensor_in)) :: periods
      CHARACTER(LEN=*), PARAMETER :: routineN = 'dbcsr_t_remap', &
         routineP = moduleN//':'//routineN
      LOGICAL                               :: nodata_prv

      CALL timeset(routineN, handle)

      IF (PRESENT(name)) THEN
         name_tmp = name
      ELSE
         name_tmp = tensor_in%name
      ENDIF
      IF(PRESENT(dist1)) THEN
         DBCSR_ASSERT(PRESENT(mp_dims_1))
      ENDIF

      IF(PRESENT(dist2)) THEN
         DBCSR_ASSERT(PRESENT(mp_dims_2))
      ENDIF

      IF(PRESENT(comm_2d)) THEN
         comm_2d_prv = comm_2d
      ELSE
         comm_2d_prv = tensor_in%comm_2d
      ENDIF

      comm_nd = dbcsr_t_nd_mp_comm(comm_2d_prv, map1_2d, map2_2d, dims1_nd=mp_dims_1, dims2_nd=mp_dims_2)
      CALL mp_environ(comm_nd, dbcsr_t_ndims(tensor_in), pdims, myploc, periods)

#:for ndim in ndims
      IF (ndims_tensor(tensor_in) == ${ndim}$) THEN
         CALL get_arrays(tensor_in%blk_sizes, ${varlist("blk_sizes", nmax=ndim)}$)
      ENDIF
#:endfor

! todo: first loop may be redundant
#:for ndim in ndims
      IF (ndims_tensor(tensor_in) == ${ndim}$) THEN
#:for idim in range(1, ndim+1)
         IF(PRESENT (dist1)) THEN
            IF(ANY(map1_2d == ${idim}$)) THEN
               i = MINLOC(map1_2d, dim=1, mask=map1_2d==${idim}$) ! i is location of idim in map1_2d
               CALL get_ith_array(dist1, i, nd_dist_${idim}$)
            ENDIF
         ENDIF

         IF(PRESENT (dist2)) THEN
            IF(ANY(map2_2d == ${idim}$)) THEN
               i = MINLOC(map2_2d, dim=1, mask=map2_2d==${idim}$) ! i is location of idim in map2_2d
               CALL get_ith_array(dist2, i, nd_dist_${idim}$)
            ENDIF
         ENDIF

         IF(.NOT. ALLOCATED(nd_dist_${idim}$)) THEN
            CALL cyclic_dist(nd_dist_${idim}$, SIZE(blk_sizes_${idim}$), pdims(${idim}$))
         ENDIF
#:endfor
         CALL dbcsr_t_distribution_new(dist, comm_nd, map1_2d, map2_2d, &
                                       ${varlist("nd_dist", nmax=ndim)}$, own_comm=.TRUE.)
      ENDIF
#:endfor

#:for ndim in ndims
      IF (ndims_tensor(tensor_in) == ${ndim}$) THEN
         CALL dbcsr_t_create(tensor_out, name_tmp, dist, &
                             map1_2d, map2_2d, dbcsr_r_get_data_type(tensor_in%matrix_rep),&
                             ${varlist("blk_sizes", nmax=ndim)}$)
      ENDIF
#:endfor

      IF(PRESENT(nodata)) THEN
         nodata_prv = nodata
      ELSE
         nodata_prv = .FALSE.
      ENDIF

      IF(.NOT. nodata_prv) CALL dbcsr_t_copy(tensor_in, tensor_out)
      CALL dbcsr_t_distribution_destroy(dist)

      CALL timestop(handle)
   END SUBROUTINE

! **************************************************************************************************
!> \brief Align index with data
!> \param tensor_in ...
!> \param tensor_out ...
!> \param order permutation resulting from alignment
! **************************************************************************************************
   SUBROUTINE dbcsr_t_align_index(tensor_in, tensor_out, order)
      TYPE(dbcsr_t_type), INTENT(INOUT)               :: tensor_in
      TYPE(dbcsr_t_type), INTENT(OUT)                 :: tensor_out
      INTEGER, DIMENSION(:), ALLOCATABLE              :: map1_2d, map2_2d
      INTEGER, DIMENSION(ndims_tensor(tensor_in)), &
         INTENT(OUT), OPTIONAL                        :: order
      INTEGER, DIMENSION(ndims_tensor(tensor_in))     :: order_prv
      CHARACTER(LEN=*), PARAMETER :: routineN = 'dbcsr_t_align_index', &
         routineP = moduleN//':'//routineN
      INTEGER                                         :: handle

      CALL timeset(routineN, handle)

      CALL get_mapping_info(tensor_in%nd_index_blk, map1_2d=map1_2d, map2_2d=map2_2d)
      order_prv = dbcsr_t_inverse_order([map1_2d, map2_2d])
      CALL dbcsr_t_permute_index(tensor_in, tensor_out, order=order_prv)

      IF(PRESENT(order)) order = order_prv

      CALL timestop(handle)
   END SUBROUTINE

! **************************************************************************************************
!> \brief Create new tensor by reordering index, data is copied exactly
!> \param tensor_in ...
!> \param tensor_out ...
!> \param order ...
! **************************************************************************************************
   SUBROUTINE dbcsr_t_permute_index(tensor_in, tensor_out, order)
      TYPE(dbcsr_t_type), INTENT(INOUT)                  :: tensor_in
      TYPE(dbcsr_t_type), INTENT(OUT)                 :: tensor_out
      INTEGER, DIMENSION(ndims_tensor(tensor_in)), &
         INTENT(IN)                                   :: order

      TYPE(nd_to_2d_mapping)                          :: nd_index_blk_rs, nd_index_rs, nd_index_grid_rs
      CHARACTER(LEN=*), PARAMETER :: routineN = 'dbcsr_t_permute_index', &
         routineP = moduleN//':'//routineN
      INTEGER                                         :: handle

      CALL timeset(routineN, handle)

      CALL permute_index(tensor_in%nd_index, nd_index_rs, order)
      CALL permute_index(tensor_in%nd_index_blk, nd_index_blk_rs, order)
      CALL permute_index(tensor_in%nd_index_grid, nd_index_grid_rs, order)

      CALL dbcsr_r_create(tensor_in%matrix_rep, tensor_out%matrix_rep)
      CALL dbcsr_r_finalize(tensor_out%matrix_rep)

      CALL dbcsr_r_copy(tensor_out%matrix_rep, tensor_in%matrix_rep, shallow_data=.TRUE.)
      tensor_out%nd_index = nd_index_rs
      tensor_out%nd_index_blk = nd_index_blk_rs
      tensor_out%nd_index_grid = nd_index_grid_rs
      tensor_out%comm_nd = tensor_in%comm_nd
      tensor_out%comm_2d = tensor_in%comm_2d
      tensor_out%refcount => tensor_in%refcount
      CALL dbcsr_t_hold(tensor_out)

      CALL reorder_arrays(tensor_in%blk_sizes, tensor_out%blk_sizes, order)
      CALL reorder_arrays(tensor_in%blk_offsets, tensor_out%blk_offsets, order)
      CALL reorder_arrays(tensor_in%nd_dist, tensor_out%nd_dist, order)
      tensor_out%name = tensor_in%name
      tensor_out%valid = .TRUE.

      CALL timestop(handle)
   END SUBROUTINE

! **************************************************************************************************
!> \brief Split tensor blocks into smaller blocks of maximum size PRODUCT(block_sizes).
!> \param tensor_in Input tensor
!> \param tensor_out Output tensor (splitted blocks)
!> \param block_sizes block sizes for each of the tensor dimensions
! **************************************************************************************************
   SUBROUTINE dbcsr_t_split_blocks(tensor_in, tensor_out, block_sizes)
      TYPE(dbcsr_t_type), INTENT(INOUT)               :: tensor_in
      TYPE(dbcsr_t_type), INTENT(OUT)                 :: tensor_out
      INTEGER, DIMENSION(dbcsr_t_ndims(tensor_in)), &
         INTENT(IN)                                   :: block_sizes

      TYPE(dbcsr_t_distribution_type)                 :: dist_old, dist_split
      TYPE(dbcsr_t_iterator_type)                     :: iter
      INTEGER, DIMENSION(:), ALLOCATABLE              :: ${varlist("nd_dist_split")}$
      INTEGER, DIMENSION(:), ALLOCATABLE              :: ${varlist("nd_blk_size_split")}$
      INTEGER, DIMENSION(:), ALLOCATABLE              :: ${varlist("index_split_offset")}$
      INTEGER                                         :: ${varlist("split_blk")}$
      INTEGER :: idim, i, isplit_sum, blk_remainder, blk, nsplit, isplit
      INTEGER, DIMENSION(:), ALLOCATABLE :: dist_d, blk_size_d, blk_size_split_d, dist_split_d, &
         map1_2d, map2_2d
      INTEGER, DIMENSION(dbcsr_t_ndims(tensor_in)) :: blk_index, blk_size, blk_offset, nsplit_blk, &
         isplit_sum_blk, blk_remainder_nd, blk_shape
      INTEGER, DIMENSION(${maxdim}$) :: split_index

#:for dparam, dtype, dsuffix in dtype_float_list
#:for ndim in ndims
      ${dtype}$, DIMENSION(${shape_colon(n=ndim)}$), ALLOCATABLE :: block_${dsuffix}$_${ndim}$d
#:endfor
#:endfor

      dist_old = dbcsr_t_distribution(tensor_in)

      DO idim = 1, ndims_tensor(tensor_in)
         CALL get_ith_array(dist_old%nd_dist, idim, dist_d)
         CALL get_ith_array(tensor_in%blk_sizes, idim, blk_size_d)

#:for idim in range(1, maxdim+1)
         IF (idim == ${idim}$) THEN
            ALLOCATE(index_split_offset_${idim}$(SIZE(dist_d)))
         ENDIF
#:endfor

         isplit_sum=0
         DO i = 1, SIZE(dist_d)
            nsplit = (blk_size_d(i)+block_sizes(idim)-1)/block_sizes(idim)
            isplit_sum = isplit_sum+nsplit
         ENDDO
         ALLOCATE(dist_split_d(isplit_sum))
         ALLOCATE(blk_size_split_d(isplit_sum))

         isplit_sum = 0
         DO i = 1, SIZE(dist_d)
            nsplit = (blk_size_d(i)+block_sizes(idim)-1)/block_sizes(idim)
            blk_remainder = blk_size_d(i)
#:for idim in range(1, maxdim+1)
            IF(idim == ${idim}$) THEN
               index_split_offset_${idim}$(i) = isplit_sum
            ENDIF
#:endfor
            DO isplit = 1, nsplit
               isplit_sum = isplit_sum + 1
               dist_split_d(isplit_sum) = dist_d(i)
               blk_size_split_d(isplit_sum) = MIN(block_sizes(idim), blk_remainder)
               blk_remainder = blk_remainder - block_sizes(idim)
            ENDDO

         ENDDO

#:for idim in range(1, maxdim+1)
         IF(idim == ${idim}$) THEN
            CALL allocate_any(nd_dist_split_${idim}$, source=dist_split_d)
            CALL allocate_any(nd_blk_size_split_${idim}$, source=blk_size_split_d)
         ENDIF
#:endfor
         DEALLOCATE(dist_split_d)
         DEALLOCATE(blk_size_split_d)

      ENDDO

      CALL get_mapping_info(tensor_in%nd_index_blk, map1_2d=map1_2d, map2_2d=map2_2d)

#:for ndim in ndims
      IF (ndims_tensor(tensor_in) == ${ndim}$) THEN
         CALL dbcsr_t_distribution_new(dist_split, tensor_in%comm_nd, map1_2d, map2_2d, &
                                       ${varlist("nd_dist_split", nmax=ndim)}$)
         CALL dbcsr_t_create(tensor_out, tensor_in%name, dist_split, map1_2d, map2_2d, &
                             dbcsr_t_get_data_type(tensor_in), ${varlist("nd_blk_size_split", nmax=ndim)}$)
      ENDIF
#:endfor

      CALL dbcsr_t_distribution_destroy(dist_split)

      CALL dbcsr_t_iterator_start(iter, tensor_in)

      DO WHILE (dbcsr_t_iterator_blocks_left(iter))
         CALL dbcsr_t_iterator_next_block(iter, blk_index, blk, blk_size=blk_size, blk_offset=blk_offset)
#:for dparam, dtype, dsuffix in dtype_float_list
         IF (dbcsr_t_get_data_type(tensor_in) == ${dparam}$) THEN
#:for ndim in ndims
            IF (dbcsr_t_ndims(tensor_in) == ${ndim}$) THEN
               CALL dbcsr_t_get_block(tensor_in, blk_index, block_${dsuffix}$_${ndim}$d)
            ENDIF
#:endfor
         ENDIF
#:endfor

         nsplit_blk = (blk_size + block_sizes - 1)/block_sizes

#:for ndim in ndims
         IF (dbcsr_t_ndims(tensor_in) == ${ndim}$) THEN
#:for idim in range(1,ndim+1)

            blk_remainder_nd(${idim}$) = blk_size(${idim}$)
            isplit_sum_blk(${idim}$) = 0
            DO split_blk_${idim}$ = 1, nsplit_blk(${idim}$)
               isplit_sum_blk(${idim}$) = isplit_sum_blk(${idim}$) + 1

               blk_shape(${idim}$) = MIN(block_sizes(${idim}$), blk_remainder_nd(${idim}$))
               blk_remainder_nd(${idim}$) = blk_remainder_nd(${idim}$) - block_sizes(${idim}$)
#:endfor

               split_index(1:${ndim}$) = [${", ".join(["index_split_offset_"+str(idim)+"(blk_index("+str(idim)+")) + isplit_sum_blk("+str(idim)+")" for idim in range(1, ndim+1)])}$]

#:for dparam, dtype, dsuffix in dtype_float_list

               IF (dbcsr_t_get_data_type(tensor_in) == ${dparam}$) THEN

                  CALL dbcsr_t_put_block(tensor_out, split_index(1:${ndim}$), &
                                         blk_shape, &
                                         block_${dsuffix}$_${ndim}$d( &
                                         ${", ".join(["(isplit_sum_blk("+str(idim)+") - 1)*block_sizes("+str(idim)+") + 1:MIN(isplit_sum_blk("+str(idim)+")*block_sizes("+str(idim)+"), blk_size("+str(idim)+"))" for idim in range(1, ndim+1)])}$))
               ENDIF
#:endfor

#:for idim in range(1,ndim+1)
            ENDDO
#:endfor

#:for dparam, dtype, dsuffix in dtype_float_list
            IF (dbcsr_t_get_data_type(tensor_in) == ${dparam}$) THEN
               DEALLOCATE(block_${dsuffix}$_${ndim}$d)
            ENDIF
#:endfor
         ENDIF
#:endfor
      ENDDO
      CALL dbcsr_t_iterator_stop(iter)

      CALL dbcsr_r_finalize(tensor_out%matrix_rep)

   END SUBROUTINE

END MODULE
