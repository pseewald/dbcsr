!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2018  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief unit testing for rectangular matrices
! **************************************************************************************************

PROGRAM dbcsr_rectangular_unittest
   USE dbcsr_api,                       ONLY: dbcsr_finalize_lib,&
                                              dbcsr_init_lib,&
                                              dbcsr_type_no_symmetry,&
                                              dbcsr_type_real_8
   USE dbcsr_data_methods,              ONLY: dbcsr_scalar
   USE dbcsr_dist_methods,              ONLY: dbcsr_distribution_new,&
                                              dbcsr_distribution_release
   USE dbcsr_methods,                   ONLY: dbcsr_col_block_sizes,&
                                              dbcsr_get_nze,&
                                              dbcsr_nblkcols_total,&
                                              dbcsr_nblkrows_total,&
                                              dbcsr_release,&
                                              dbcsr_row_block_sizes
   USE dbcsr_mp_methods,                ONLY: dbcsr_mp_release
   USE dbcsr_multiply_api,              ONLY: dbcsr_multiply
   USE dbcsr_operations,                ONLY: dbcsr_add,&
                                              dbcsr_frobenius_norm
   USE dbcsr_rectangular_base,          ONLY: dbcsr_r_convert_to_dbcsr,&
                                              dbcsr_r_destroy,&
                                              dbcsr_r_info,&
                                              dbcsr_r_nblkcols_total,&
                                              dbcsr_r_nblkrows_total,&
                                              dbcsr_r_create
   USE dbcsr_rectangular_base_types,    ONLY: dbcsr_r_type
   USE dbcsr_rectangular_mm,            ONLY: dbcsr_r_multiply
   USE dbcsr_rectangular_split,         ONLY: dbcsr_r_get_split_info,&
                                              dbcsr_r_print_info
   USE dbcsr_rectangular_test,          ONLY: dbcsr_r_random_bsizes,&
                                              dbcsr_r_setup_test_matrix
   USE dbcsr_rectangular_util,          ONLY: dbcsr_mp_environ
   USE dbcsr_transformations,           ONLY: dbcsr_complete_redistribute
   USE dbcsr_types,                     ONLY: dbcsr_distribution_obj,&
                                              dbcsr_mp_obj,&
                                              dbcsr_no_transpose,&
                                              dbcsr_transpose,&
                                              dbcsr_type
   USE dbcsr_dist_util,                      ONLY: dbcsr_checksum
   USE dbcsr_work_operations,           ONLY: dbcsr_create,&
                                              dbcsr_finalize
   USE dbcsr_kinds,                           ONLY: int_8,&
                                              real_8
   USE dbcsr_machine,                         ONLY: default_output_unit
   USE dbcsr_mpiwrap,                 ONLY: mp_cart_create,&
                                              mp_comm_free,&
                                              mp_environ,&
                                              mp_world_finalize,&
                                              mp_world_init
#include "base/dbcsr_base_uses.f90"

   IMPLICIT NONE

   INTEGER(KIND=int_8), PARAMETER :: m = 100, k = 20, n = 10
   TYPE(dbcsr_r_type) :: A, B, C, At, Bt, Ct, A_out, B_out, C_out, At_out, Bt_out, Ct_out
   INTEGER, DIMENSION(m) :: bsize_m
   INTEGER, DIMENSION(n) :: bsize_n
   INTEGER, DIMENSION(k) :: bsize_k
   REAL(KIND=real_8), PARAMETER :: sparsity = 0.1
   REAL(KIND=real_8), PARAMETER :: test_tol = 1.0E-10_real_8
   INTEGER :: mp_comm
   INTEGER :: numnodes, mynode
   INTEGER :: io_unit
   INTEGER :: mp_comm_A, mp_comm_At, mp_comm_B, mp_comm_Bt, mp_comm_C, mp_comm_Ct

   CALL mp_world_init(mp_comm)

   CALL dbcsr_init_lib()

   CALL dbcsr_r_random_bsizes([13, 8, 5, 25, 12], 2, bsize_m)
   CALL dbcsr_r_random_bsizes([3, 78, 33, 12, 3, 15], 1, bsize_n)
   CALL dbcsr_r_random_bsizes([9, 64, 23, 2], 3, bsize_k)

   CALL mp_environ(numnodes, mynode, mp_comm)

   io_unit = 0
   IF (mynode .EQ. 0) io_unit = default_output_unit

   CALL dbcsr_r_setup_test_matrix(A, mp_comm_A, mp_comm, m, k, bsize_m, bsize_k, [5, 1], "A", sparsity)
   CALL dbcsr_r_setup_test_matrix(At, mp_comm_At, mp_comm, k, m, bsize_k, bsize_m, [3, 8], "A^t", sparsity)
   CALL dbcsr_r_setup_test_matrix(B, mp_comm_B, mp_comm, n, m, bsize_n, bsize_m, [3, 2], "B", sparsity)
   CALL dbcsr_r_setup_test_matrix(Bt, mp_comm_Bt, mp_comm, m, n, bsize_m, bsize_n, [1, 3], "B^t", sparsity)
   CALL dbcsr_r_setup_test_matrix(C, mp_comm_C, mp_comm, k, n, bsize_k, bsize_n, [5, 7], "C", sparsity)
   CALL dbcsr_r_setup_test_matrix(Ct, mp_comm_Ct, mp_comm, n, k, bsize_n, bsize_k, [1, 1], "C^t", sparsity)

   CALL dbcsr_r_create(A, A_out)
   CALL dbcsr_r_create(At, At_out)
   CALL dbcsr_r_create(B, B_out)
   CALL dbcsr_r_create(Bt, Bt_out)
   CALL dbcsr_r_create(C, C_out)
   CALL dbcsr_r_create(Ct, Ct_out)

   IF (mynode == 0) WRITE (io_unit, '(A, 1X, A, 1X, A, I10, 1X, A, 1X, I10)') "DBCSR RECT |", "split info for matrix", &
      TRIM(A%matrix%name), dbcsr_r_nblkrows_total(A), 'X', dbcsr_r_nblkcols_total(A)
   CALL dbcsr_r_print_info(dbcsr_r_info(A))
   IF (mynode == 0) WRITE (io_unit, '(A, 1X, A, 1X, A, I10, 1X, A, 1X, I10)') "DBCSR RECT |", "split info for matrix", &
      TRIM(At%matrix%name), dbcsr_r_nblkrows_total(At), 'X', dbcsr_r_nblkcols_total(At)
   CALL dbcsr_r_print_info(dbcsr_r_info(At))
   IF (mynode == 0) WRITE (io_unit, '(A, 1X, A, 1X, A, I10, 1X, A, 1X, I10)') "DBCSR RECT |", "split info for matrix", &
      TRIM(B%matrix%name), dbcsr_r_nblkrows_total(B), 'X', dbcsr_r_nblkcols_total(B)
   CALL dbcsr_r_print_info(dbcsr_r_info(B))
   IF (mynode == 0) WRITE (io_unit, '(A, 1X, A, 1X, A, I10, 1X, A, 1X, I10)') "DBCSR RECT |", "split info for matrix", &
      TRIM(Bt%matrix%name), dbcsr_r_nblkrows_total(Bt), 'X', dbcsr_r_nblkcols_total(Bt)
   CALL dbcsr_r_print_info(dbcsr_r_info(Bt))
   IF (mynode == 0) WRITE (io_unit, '(A, 1X, A, 1X, A, I10, 1X, A, 1X, I10)') "DBCSR RECT |", "split info for matrix", &
      TRIM(C%matrix%name), dbcsr_r_nblkrows_total(C), 'X', dbcsr_r_nblkcols_total(C)
   CALL dbcsr_r_print_info(dbcsr_r_info(C))
   IF (mynode == 0) WRITE (io_unit, '(A, 1X, A, 1X, A, I10, 1X, A, 1X, I10)') "DBCSR RECT |", "split info for matrix", &
      TRIM(Ct%matrix%name), dbcsr_r_nblkrows_total(Ct), 'X', dbcsr_r_nblkcols_total(Ct)
   CALL dbcsr_r_print_info(dbcsr_r_info(Ct))

   CALL test_mm('N', 'N', 'N', B, A, Ct_out, default_output_unit)
   CALL test_mm('T', 'N', 'N', Bt, A, Ct_out, default_output_unit)
   CALL test_mm('N', 'T', 'N', B, At, Ct_out, default_output_unit)
   CALL test_mm('T', 'T', 'N', Bt, At, Ct_out, default_output_unit)
   CALL test_mm('N', 'N', 'T', B, A, C_out, default_output_unit)
   CALL test_mm('T', 'N', 'T', Bt, A, C_out, default_output_unit)
   CALL test_mm('N', 'T', 'T', B, At, C_out, default_output_unit)
   CALL test_mm('T', 'T', 'T', Bt, At, C_out, default_output_unit)

   CALL test_mm('N', 'N', 'N', A, C, Bt_out, default_output_unit)
   CALL test_mm('T', 'N', 'N', At, C, Bt_out, default_output_unit)
   CALL test_mm('N', 'T', 'N', A, Ct, Bt_out, default_output_unit)
   CALL test_mm('T', 'T', 'N', At, Ct, Bt_out, default_output_unit)

   CALL test_mm('N', 'N', 'T', A, C, B_out, default_output_unit)
   CALL test_mm('T', 'N', 'T', At, C, B_out, default_output_unit)
   CALL test_mm('N', 'T', 'T', A, Ct, B_out, default_output_unit)
   CALL test_mm('T', 'T', 'T', At, Ct, B_out, default_output_unit)

   CALL test_mm('N', 'N', 'N', C, B, At_out, default_output_unit)
   CALL test_mm('T', 'N', 'N', Ct, B, At_out, default_output_unit)
   CALL test_mm('N', 'T', 'N', C, Bt, At_out, default_output_unit)
   CALL test_mm('T', 'T', 'N', Ct, Bt, At_out, default_output_unit)

   CALL test_mm('N', 'N', 'T', C, B, A_out, default_output_unit)
   CALL test_mm('T', 'N', 'T', Ct, B, A_out, default_output_unit)
   CALL test_mm('N', 'T', 'T', C, Bt, A_out, default_output_unit)
   CALL test_mm('T', 'T', 'T', Ct, Bt, A_out, default_output_unit)

   CALL dbcsr_r_destroy(A)
   CALL dbcsr_r_destroy(At)
   CALL dbcsr_r_destroy(B)
   CALL dbcsr_r_destroy(Bt)
   CALL dbcsr_r_destroy(C)
   CALL dbcsr_r_destroy(Ct)
   CALL dbcsr_r_destroy(A_out)
   CALL dbcsr_r_destroy(At_out)
   CALL dbcsr_r_destroy(B_out)
   CALL dbcsr_r_destroy(Bt_out)
   CALL dbcsr_r_destroy(C_out)
   CALL dbcsr_r_destroy(Ct_out)

   CALL mp_comm_free(mp_comm_A)
   CALL mp_comm_free(mp_comm_At)
   CALL mp_comm_free(mp_comm_B)
   CALL mp_comm_free(mp_comm_Bt)
   CALL mp_comm_free(mp_comm_C)
   CALL mp_comm_free(mp_comm_Ct)

   CALL dbcsr_finalize_lib(mp_comm, io_unit)

   CALL mp_world_finalize()

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param transa ...
!> \param transb ...
!> \param matrix_a ...
!> \param matrix_b ...
!> \param unit_nr ...
! **************************************************************************************************
   SUBROUTINE test_mm(transa, transb, transc, matrix_a, matrix_b, matrix_c, unit_nr)
      CHARACTER(LEN=1), INTENT(IN)                       :: transa, transb, transc
      TYPE(dbcsr_r_type), INTENT(INOUT)                  :: matrix_a, matrix_b, matrix_c
      INTEGER, INTENT(IN)                                :: unit_nr

      CHARACTER(LEN=1)                                   :: transa_prv, transb_prv
      INTEGER                                            :: comm_dbcsr, io_unit, mp_comm, mynode, &
                                                            numnodes
      INTEGER, DIMENSION(2)                              :: myploc, npdims
      INTEGER, DIMENSION(:), POINTER                     :: cd_a, cd_b, cd_c, rd_a, rd_b, rd_c
      REAL(KIND=real_8)                                  :: norm, rc_cs, sq_cs
      TYPE(dbcsr_distribution_obj)                       :: dist_a, dist_b, dist_c
      TYPE(dbcsr_mp_obj)                                 :: mp_environ_tmp
      TYPE(dbcsr_type)                                   :: dbcsr_a, dbcsr_a_mm, dbcsr_b, &
                                                            dbcsr_b_mm, dbcsr_c, dbcsr_c_mm, &
                                                            dbcsr_c_mm_check

      CALL dbcsr_r_get_split_info(dbcsr_r_info(matrix_a), mp_comm=mp_comm)
      CALL mp_environ(numnodes, mynode, mp_comm)
      io_unit = 0
      IF (mynode .EQ. 0) io_unit = unit_nr

      CALL dbcsr_r_multiply(transa, transb, transc, dbcsr_scalar(1.0_real_8), matrix_a, matrix_b, matrix_c)

      CALL dbcsr_r_convert_to_dbcsr(matrix_a, dbcsr_a)
      CALL dbcsr_r_convert_to_dbcsr(matrix_b, dbcsr_b)
      CALL dbcsr_r_convert_to_dbcsr(matrix_c, dbcsr_c)

      npdims(:) = 0
      CALL mp_cart_create(mp_comm, 2, npdims, myploc, comm_dbcsr)

      CALL random_dist(rd_a, dbcsr_nblkrows_total(dbcsr_a), npdims(1))
      CALL random_dist(cd_a, dbcsr_nblkcols_total(dbcsr_a), npdims(2))
      CALL random_dist(rd_b, dbcsr_nblkrows_total(dbcsr_b), npdims(1))
      CALL random_dist(cd_b, dbcsr_nblkcols_total(dbcsr_b), npdims(2))
      CALL random_dist(rd_c, dbcsr_nblkrows_total(dbcsr_c), npdims(1))
      CALL random_dist(cd_c, dbcsr_nblkcols_total(dbcsr_c), npdims(2))

      mp_environ_tmp = dbcsr_mp_environ(comm_dbcsr)
      CALL dbcsr_distribution_new(dist_a, mp_environ_tmp, rd_a, cd_a, reuse_arrays=.TRUE.)
      CALL dbcsr_distribution_new(dist_b, mp_environ_tmp, rd_b, cd_b, reuse_arrays=.TRUE.)
      CALL dbcsr_distribution_new(dist_c, mp_environ_tmp, rd_c, cd_c, reuse_arrays=.TRUE.)
      CALL dbcsr_mp_release(mp_environ_tmp)

      CALL dbcsr_create(matrix=dbcsr_a_mm, name="matrix a", dist=dist_a, matrix_type=dbcsr_type_no_symmetry, &
                        row_blk_size=dbcsr_row_block_sizes(dbcsr_a), col_blk_size=dbcsr_col_block_sizes(dbcsr_a), &
                        data_type=dbcsr_type_real_8)

      CALL dbcsr_create(matrix=dbcsr_b_mm, name="matrix b", dist=dist_b, matrix_type=dbcsr_type_no_symmetry, &
                        row_blk_size=dbcsr_row_block_sizes(dbcsr_b), col_blk_size=dbcsr_col_block_sizes(dbcsr_b), &
                        data_type=dbcsr_type_real_8)

      CALL dbcsr_create(matrix=dbcsr_c_mm, name="matrix c", dist=dist_c, matrix_type=dbcsr_type_no_symmetry, &
                        row_blk_size=dbcsr_row_block_sizes(dbcsr_c), col_blk_size=dbcsr_col_block_sizes(dbcsr_c), &
                        data_type=dbcsr_type_real_8)

      CALL dbcsr_create(matrix=dbcsr_c_mm_check, name="matrix c check", dist=dist_c, matrix_type=dbcsr_type_no_symmetry, &
                        row_blk_size=dbcsr_row_block_sizes(dbcsr_c), col_blk_size=dbcsr_col_block_sizes(dbcsr_c), &
                        data_type=dbcsr_type_real_8)

      CALL dbcsr_finalize(dbcsr_a_mm)
      CALL dbcsr_finalize(dbcsr_b_mm)
      CALL dbcsr_finalize(dbcsr_c_mm)
      CALL dbcsr_finalize(dbcsr_c_mm_check)

      CALL dbcsr_complete_redistribute(dbcsr_a, dbcsr_a_mm)
      CALL dbcsr_complete_redistribute(dbcsr_b, dbcsr_b_mm)
      CALL dbcsr_complete_redistribute(dbcsr_c, dbcsr_c_mm_check)

      transa_prv = transa; transb_prv = transb

      IF (transc == dbcsr_no_transpose) THEN
         CALL dbcsr_multiply(transa_prv, transb_prv, dbcsr_scalar(1.0_real_8), &
                             dbcsr_a_mm, dbcsr_b_mm, dbcsr_scalar(0.0_real_8), dbcsr_c_mm)
      ELSEIF (transc == dbcsr_transpose) THEN
         CALL invert_transpose_flag(transa_prv)
         CALL invert_transpose_flag(transb_prv)
         CALL dbcsr_multiply(transb_prv, transa_prv, dbcsr_scalar(1.0_real_8), &
                             dbcsr_b_mm, dbcsr_a_mm, dbcsr_scalar(0.0_real_8), dbcsr_c_mm)
      ENDIF

      sq_cs = dbcsr_checksum(dbcsr_c_mm)
      rc_cs = dbcsr_checksum(dbcsr_c_mm_check)
      CALL dbcsr_add(dbcsr_c_mm_check, dbcsr_c_mm, -1.0_real_8, 1.0_real_8)
      norm = dbcsr_frobenius_norm(dbcsr_c_mm_check)
      norm = norm/dbcsr_get_nze(dbcsr_c_mm_check)

      IF (io_unit > 0) THEN
      IF (ABS(norm) .GT. test_tol) THEN
         WRITE (io_unit, '(A, A, A, A, A, 1X, A)') TRIM(matrix_a%matrix%name), transa, ' X ', TRIM(matrix_b%matrix%name), &
            transb, 'failed!'
         WRITE (io_unit, *) "checksums", sq_cs, rc_cs
         WRITE (io_unit, *) "difference norm", norm
         DBCSR_ABORT("")
      ELSE
         WRITE (io_unit, '(A, A, A, A, A, 1X, A)') TRIM(matrix_a%matrix%name), transa, ' X ', TRIM(matrix_b%matrix%name), &
            transb, 'passed!'
         WRITE (io_unit, *) "checksums", sq_cs, rc_cs
         WRITE (io_unit, *) "difference norm", norm
      ENDIF
      ENDIF

      CALL dbcsr_release(dbcsr_a)
      CALL dbcsr_release(dbcsr_a_mm)
      CALL dbcsr_release(dbcsr_b)
      CALL dbcsr_release(dbcsr_b_mm)
      CALL dbcsr_release(dbcsr_c)
      CALL dbcsr_release(dbcsr_c_mm)
      CALL dbcsr_release(dbcsr_c_mm_check)

      CALL dbcsr_distribution_release(dist_a)
      CALL dbcsr_distribution_release(dist_b)
      CALL dbcsr_distribution_release(dist_c)

      CALL mp_comm_free(comm_dbcsr)

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param trans_flag ...
! **************************************************************************************************
   SUBROUTINE invert_transpose_flag(trans_flag)
      CHARACTER(LEN=1), INTENT(INOUT)                    :: trans_flag

      IF (trans_flag == dbcsr_transpose) THEN
         trans_flag = dbcsr_no_transpose
      ELSEIF (trans_flag == dbcsr_no_transpose) THEN
         trans_flag = dbcsr_transpose
      ENDIF
   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param dist_array ...
!> \param dist_size ...
!> \param nbins ...
! **************************************************************************************************
   SUBROUTINE random_dist(dist_array, dist_size, nbins)
      INTEGER, DIMENSION(:), INTENT(out), POINTER        :: dist_array
      INTEGER, INTENT(in)                                :: dist_size, nbins

      INTEGER                                            :: i

      ALLOCATE (dist_array(dist_size))
      DO i = 1, dist_size
         dist_array(i) = MODULO(nbins-i, nbins)
      END DO

   END SUBROUTINE random_dist

END PROGRAM
