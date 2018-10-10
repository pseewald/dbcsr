!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2018  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief testing infrastructure for rectangular matrices
! **************************************************************************************************

MODULE dbcsr_rectangular_test
   USE dbcsr_types,                       ONLY: dbcsr_type_real_8
   USE dbcsr_data_methods,              ONLY: dbcsr_scalar
   USE dbcsr_methods,                   ONLY: dbcsr_release
   USE dbcsr_multiply_api,              ONLY: dbcsr_multiply
   USE dbcsr_rectangular_base,          ONLY: &
        dbcsr_r_convert_to_dbcsr, dbcsr_r_create, dbcsr_r_destroy, dbcsr_r_distribution_new, &
        dbcsr_r_finalize, dbcsr_r_get_stored_coordinates, dbcsr_r_nblkcols_total, &
        dbcsr_r_nblkrows_total, dbcsr_r_put_block
   USE dbcsr_rectangular_base_types,    ONLY: dbcsr_r_distribution_type,&
                                              dbcsr_r_type
   USE dbcsr_rectangular_global,        ONLY: dbcsr_r_blk_size_arb,&
                                              dbcsr_r_dist_cyclic
   USE dbcsr_rectangular_mm,            ONLY: dbcsr_r_multiply
   USE dbcsr_rectangular_split,         ONLY: dbcsr_r_mp_comm
   USE dbcsr_types,                     ONLY: dbcsr_type
   USE dbcsr_kinds,                           ONLY: int_8,&
                                              real_8
   USE dbcsr_mpiwrap,                 ONLY: mp_environ
#include "../base/dbcsr_base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: dbcsr_r_setup_test_matrix, dbcsr_r_random_bsizes, dbcsr_r_benchmark_mm

CONTAINS
! **************************************************************************************************
!> \brief ...
!> \param matrix ...
!> \param mp_comm_out ...
!> \param mp_comm ...
!> \param nrows ...
!> \param ncols ...
!> \param rbsizes ...
!> \param cbsizes ...
!> \param dist_splitsize ...
!> \param name ...
!> \param sparsity ...
!> \param reuse_comm ...
! **************************************************************************************************
   SUBROUTINE dbcsr_r_setup_test_matrix(matrix, mp_comm_out, mp_comm, nrows, ncols, rbsizes, cbsizes, &
                                        dist_splitsize, name, sparsity, reuse_comm)

      TYPE(dbcsr_r_type), INTENT(OUT)                    :: matrix
      INTEGER, INTENT(OUT)                               :: mp_comm_out
      INTEGER, INTENT(IN)                                :: mp_comm
      INTEGER(KIND=int_8), INTENT(IN)                    :: nrows, ncols
      INTEGER, DIMENSION(nrows), INTENT(IN)              :: rbsizes
      INTEGER, DIMENSION(ncols), INTENT(IN)              :: cbsizes
      INTEGER, DIMENSION(2), INTENT(IN)                  :: dist_splitsize
      CHARACTER(len=*), INTENT(IN)                       :: name
      REAL(KIND=real_8), INTENT(IN)                      :: sparsity
      LOGICAL, INTENT(IN), OPTIONAL                      :: reuse_comm

      INTEGER                                            :: col_size, max_col_size, max_nze, &
                                                            max_row_size, mynode, node_holds_blk, &
                                                            numnodes, nze, row_size
      INTEGER(KIND=int_8)                                :: col, col_s, row, row_s
      INTEGER, DIMENSION(2)                              :: pcoord, pdims
      LOGICAL                                            :: reuse_comm_prv, tr
      REAL(KIND=real_8)                                  :: rn
      REAL(KIND=real_8), ALLOCATABLE, DIMENSION(:, :)    :: values
      TYPE(dbcsr_r_blk_size_arb)                         :: cbsize_obj, rbsize_obj
      TYPE(dbcsr_r_dist_cyclic)                          :: col_dist_obj, row_dist_obj
      TYPE(dbcsr_r_distribution_type)                    :: dist

      IF (PRESENT(reuse_comm)) THEN
         reuse_comm_prv = reuse_comm
      ELSE
         reuse_comm_prv = .FALSE.
      ENDIF

      IF (reuse_comm_prv) THEN
         mp_comm_out = mp_comm
      ELSE
         mp_comm_out = dbcsr_r_mp_comm(mp_comm, nrows, ncols)
      ENDIF

      CALL mp_environ(numnodes, mynode, mp_comm_out)
      CALL mp_environ(numnodes, pdims, pcoord, mp_comm_out)

      row_dist_obj = dbcsr_r_dist_cyclic(dist_splitsize(1), pdims(1), nrows)
      col_dist_obj = dbcsr_r_dist_cyclic(dist_splitsize(2), pdims(2), ncols)

      rbsize_obj = dbcsr_r_blk_size_arb(rbsizes)
      cbsize_obj = dbcsr_r_blk_size_arb(cbsizes)

      CALL dbcsr_r_distribution_new(dist, mp_comm_out, row_dist_obj, col_dist_obj)
      CALL dbcsr_r_create(matrix, name, dist=dist, data_type=dbcsr_type_real_8, &
                          row_blk_size=rbsize_obj, col_blk_size=cbsize_obj, own_dist=.TRUE.)

      max_row_size = MAXVAL(rbsizes)
      max_col_size = MAXVAL(cbsizes)
      max_nze = max_row_size*max_col_size

      ALLOCATE (values(max_row_size, max_col_size))

      DO row = 1, dbcsr_r_nblkrows_total(matrix)
         DO col = 1, dbcsr_r_nblkcols_total(matrix)
            CALL RANDOM_NUMBER(rn)
            IF (rn .GT. sparsity) THEN
               tr = .FALSE.
               row_s = row; col_s = col
               CALL dbcsr_r_get_stored_coordinates(matrix, row_s, col_s, node_holds_blk)

               IF (node_holds_blk .EQ. mynode) THEN
                  row_size = rbsize_obj%data(row_s)
                  col_size = cbsize_obj%data(col_s)
                  nze = row_size*col_size
                  CALL RANDOM_NUMBER(values(1:row_size, 1:col_size))
                  CALL dbcsr_r_put_block(matrix, row_s, col_s, values(1:row_size, 1:col_size))
               ENDIF
            ENDIF
         ENDDO
      ENDDO

      CALL dbcsr_r_finalize(matrix)

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param transa ...
!> \param transb ...
!> \param matrix_a ...
!> \param matrix_b ...
! **************************************************************************************************
   SUBROUTINE dbcsr_r_benchmark_mm(transa, transb, transc, matrix_a, matrix_b, matrix_c, compare_dbcsr, io_unit)
      CHARACTER(LEN=1), INTENT(IN)                       :: transa, transb, transc
      TYPE(dbcsr_r_type), INTENT(INOUT)                  :: matrix_a, matrix_b, matrix_c

      INTEGER, INTENT(IN), OPTIONAL                      :: io_unit
      LOGICAL, INTENT(IN) :: compare_dbcsr

      INTEGER                                            :: handle1, handle2
      TYPE(dbcsr_type)                                   :: dbcsr_a, dbcsr_b, dbcsr_c

      CALL timeset("benchmark_rectangular_mm", handle1)
      CALL dbcsr_r_multiply(transa, transb, transc, dbcsr_scalar(1.0_real_8), matrix_a, matrix_b, matrix_c, io_unit=io_unit)
      CALL timestop(handle1)

      IF (compare_dbcsr) THEN
         CALL dbcsr_r_convert_to_dbcsr(matrix_a, dbcsr_a)
         CALL dbcsr_r_convert_to_dbcsr(matrix_b, dbcsr_b)
         CALL dbcsr_r_convert_to_dbcsr(matrix_c, dbcsr_c)

         CALL timeset("benchmark_dbcsr_mm", handle2)
         CALL dbcsr_multiply(transa, transb, dbcsr_scalar(1.0_real_8), dbcsr_a, dbcsr_b, &
                             dbcsr_scalar(0.0_real_8), dbcsr_c)
         CALL timestop(handle2)

         CALL dbcsr_release(dbcsr_a)
         CALL dbcsr_release(dbcsr_b)
         CALL dbcsr_release(dbcsr_c)
      ENDIF

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param sizes ...
!> \param repeat ...
!> \param block_sizes ...
! **************************************************************************************************
   SUBROUTINE dbcsr_r_random_bsizes(sizes, repeat, block_sizes)
      INTEGER, DIMENSION(:), INTENT(IN)                  :: sizes
      INTEGER, INTENT(IN)                                :: repeat
      INTEGER, DIMENSION(:), INTENT(OUT)                 :: block_sizes

      INTEGER                                            :: d, size_i

      DO d = 1, SIZE(block_sizes)
         size_i = MOD((d-1)/repeat, SIZE(sizes))+1
         block_sizes(d) = sizes(size_i)
      ENDDO
   END SUBROUTINE

END MODULE
