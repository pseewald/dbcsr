!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2018  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief base routines, low-level code using Fypp preprocessor.
!>
!> Mostly wrappers around existing DBCSR routines.
!> \author Patrick Seewald
! **************************************************************************************************

MODULE dbcsr_rectangular_base_low
#:include "dbcsr_rectangular.fypp"
#include "../base/dbcsr_base_uses.f90"

   USE dbcsr_rectangular_base_types, ONLY: dbcsr_r_type, dbcsr_r_iterator
   USE dbcsr_rectangular_split, ONLY: block_index_local_to_global, block_index_global_to_local
   USE dbcsr_block_access, ONLY: dbcsr_get_block_p, dbcsr_put_block
   USE dbcsr_kinds, ONLY: int_8, real_4, real_8
   USE dbcsr_iterator_operations, ONLY: dbcsr_iterator_next_block

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbcsr_rectangular_base_low'

   PUBLIC :: dbcsr_r_iterator_next_block, dbcsr_r_put_block, dbcsr_r_get_block

   INTERFACE dbcsr_r_iterator_next_block
#:for dparam, dtype, dsuffix in dtype_float_list
      MODULE PROCEDURE dbcsr_r_iterator_next_block_${dsuffix}$
#:endfor
   END INTERFACE

   INTERFACE dbcsr_r_put_block
#:for dparam, dtype, dsuffix in dtype_float_list
      MODULE PROCEDURE dbcsr_r_put_block_${dsuffix}$
#:endfor
   END INTERFACE

   INTERFACE dbcsr_r_get_block
#:for dparam, dtype, dsuffix in dtype_float_list
      MODULE PROCEDURE dbcsr_r_get_block_p_${dsuffix}$
#:endfor
   END INTERFACE

CONTAINS

#:for dparam, dtype, dsuffix in dtype_float_list
! **************************************************************************************************
!> \brief iterate over group local blocks and get index and block
!> \param iterator ...
!> \param row global block row
!> \param column global block column
!> \param block matrix block
!> \param transposed ...
!> \param block_number ...
!> \param row_size ...
!> \param col_size ...
! **************************************************************************************************
   SUBROUTINE dbcsr_r_iterator_next_block_${dsuffix}$(iterator, row, column, block, transposed, block_number, &
                                                         row_size, col_size)
      TYPE(dbcsr_r_iterator), INTENT(INOUT) :: iterator
      INTEGER(KIND=int_8), INTENT(OUT)                     :: row, column
      ${dtype}$, DIMENSION(:, :), &
        POINTER                                :: block
      LOGICAL, INTENT(OUT)                     :: transposed
      INTEGER, INTENT(OUT), OPTIONAL           :: block_number
      INTEGER, INTENT(OUT), OPTIONAL           :: row_size, col_size
      INTEGER                                  :: row_group, column_group
      CHARACTER(LEN=*), PARAMETER :: routineN = 'dbcsr_r_iterator_next_block_${dsuffix}$', &
         routineP = moduleN//':'//routineN
      INTEGER :: handle

      CALL timeset(routineN, handle)
      CALL dbcsr_iterator_next_block(iterator%iter, row_group, column_group, block, transposed, block_number, &
                                     row_size, col_size)

      CALL block_index_local_to_global(iterator%dist%info, iterator%dist, row_group=row_group, column_group=column_group, &
         row=row, column=column)

      CALL timestop(handle)
   END SUBROUTINE
#:endfor

! **************************************************************************************************
!> \brief ...
!> \param matrix ...
!> \param row ...
!> \param col ...
!> \param block ...
!> \param transposed ...
!> \param summation ...
!> \param scale ...
! **************************************************************************************************

#:for dparam, dtype, dsuffix in dtype_float_list
   SUBROUTINE dbcsr_r_put_block_${dsuffix}$(matrix, row, col, block, transposed, summation, &
                                               scale)
      TYPE(dbcsr_r_type), INTENT(INOUT) :: matrix
      ${dtype}$, DIMENSION(:,:), INTENT(IN) :: block
      LOGICAL, INTENT(IN), OPTIONAL :: transposed, summation
      INTEGER(KIND=int_8), INTENT(IN) :: row, col
      ${dtype}$, INTENT(IN), OPTIONAL      :: scale
      INTEGER :: col_group, row_group
      CHARACTER(LEN=*), PARAMETER :: routineN = 'dbcsr_r_put_block_${dsuffix}$', &
         routineP = moduleN//':'//routineN
      INTEGER :: handle

      CALL timeset(routineN, handle)
      CALL block_index_global_to_local(matrix%dist%info, matrix%dist, row=row, column=col, &
                                       row_group=row_group, column_group = col_group)

      CALL dbcsr_put_block(matrix%matrix, row_group, col_group, block, transposed=transposed, summation=summation, scale=scale)

      CALL timestop(handle)
   END SUBROUTINE
#:endfor

! **************************************************************************************************
!> \brief ...
!> \param matrix ...
!> \param row ...
!> \param col ...
!> \param block ...
!> \param transposed ...
!> \param found ...
!> \param row_size ...
!> \param col_size ...
! **************************************************************************************************

#:for dparam, dtype, dsuffix in dtype_float_list
   SUBROUTINE dbcsr_r_get_block_p_${dsuffix}$(matrix, row, col, block, transposed, found, row_size, col_size)
      TYPE(dbcsr_r_type), INTENT(INOUT)           :: matrix
      INTEGER(KIND=int_8), INTENT(IN)                      :: row, col
      ${dtype}$, DIMENSION(:,:), POINTER         :: block
      LOGICAL, INTENT(OUT)                     :: transposed
      LOGICAL, INTENT(OUT)                     :: found
      INTEGER, INTENT(OUT), OPTIONAL           :: row_size, col_size
      INTEGER :: col_group, row_group
      CHARACTER(LEN=*), PARAMETER :: routineN = 'dbcsr_r_get_block_p_${dsuffix}$', &
         routineP = moduleN//':'//routineN
      INTEGER :: handle

      CALL timeset(routineN, handle)

      CALL block_index_global_to_local(matrix%dist%info, matrix%dist, row=row, column=col, &
                                       row_group=row_group, column_group = col_group)

      CALL dbcsr_get_block_p(matrix%matrix, row_group, col_group, block, transposed, found, row_size=row_size, col_size=col_size)

      CALL timestop(handle)
   END SUBROUTINE
#:endfor

END MODULE
