!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2017  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Matrix multiplication for rectangular matrices (CARMA algorithm).
!> \author Patrick Seewald
! **************************************************************************************************

MODULE dbcsr_rectangular_mm

   USE dbcsr_data_methods, ONLY: dbcsr_scalar_zero
   USE dbcsr_data_types, ONLY: dbcsr_scalar_type
   USE dbcsr_multiply_api, ONLY: dbcsr_multiply
   USE dbcsr_rectangular_base, ONLY: &
      dbcsr_r_create, dbcsr_r_destroy, dbcsr_r_distribution_destroy, dbcsr_r_distribution_new, &
      dbcsr_r_finalize, dbcsr_r_get_data_type, dbcsr_r_info, dbcsr_r_nblkcols_total, &
      dbcsr_r_nblkrows_total, dbcsr_r_write_matrix_info
   USE dbcsr_rectangular_base_types, ONLY: dbcsr_r_distribution_type, &
      dbcsr_r_split_info, &
      dbcsr_r_type
   USE dbcsr_rectangular_global, ONLY: dbcsr_r_dist_cyclic, dbcsr_r_distribution
   USE dbcsr_rectangular_reshape, ONLY: dbcsr_r_merge, &
      dbcsr_r_replicate, &
      dbcsr_r_reshape
   USE dbcsr_rectangular_split, ONLY: colsplit, &
      dbcsr_r_get_split_info, &
      dbcsr_r_split, &
      dbcsr_r_split_rows_or_cols, &
      rowsplit, &
      world_to_group_proc_map, &
      dbcsr_r_mp_comm, &
      dbcsr_r_release_info, &
      dbcsr_r_print_info
   USE dbcsr_rectangular_util, ONLY: swap
   USE dbcsr_types, ONLY: dbcsr_no_transpose, &
      dbcsr_transpose
   USE dbcsr_kinds, ONLY: default_string_length, &
      int_8, &
      real_8
   USE dbcsr_mpiwrap, ONLY: mp_comm_compare, &
      mp_environ, mp_sum, mp_comm_free
   USE dbcsr_io, ONLY: dbcsr_print
   USE dbcsr_operations, ONLY: dbcsr_scale
#include "../base/dbcsr_base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbcsr_rectangular_mm'

   PUBLIC :: dbcsr_r_multiply

CONTAINS

   ! TODO: recursive splitting
! **************************************************************************************************
!> \brief ...
!> \param transa ...
!> \param transb ...
!> \param transc ...
!> \param alpha ...
!> \param matrix_a ...
!> \param matrix_b ...
!> \param matrix_c ...
! **************************************************************************************************
   SUBROUTINE dbcsr_r_multiply(transa, transb, transc, alpha, matrix_a, matrix_b, beta, matrix_c, &
                               optimize_dist, filter_eps, flop, io_unit)
      CHARACTER(LEN=1), INTENT(IN)                       :: transa, transb, transc
      TYPE(dbcsr_scalar_type), INTENT(IN) :: alpha
      TYPE(dbcsr_r_type), TARGET, INTENT(INOUT) :: matrix_a, matrix_b
      TYPE(dbcsr_r_type), POINTER :: matrix_b_rs => NULL(), matrix_a_rs => NULL(), matrix_c_rs => NULL()
      TYPE(dbcsr_r_type) :: matrix_c_rep, matrix_b_rep, matrix_a_rep
      TYPE(dbcsr_scalar_type), INTENT(IN) :: beta
      TYPE(dbcsr_r_type), TARGET, INTENT(INOUT) :: matrix_c
      LOGICAL, INTENT(IN), OPTIONAL             :: optimize_dist
      REAL(KIND=real_8), INTENT(IN), OPTIONAL            :: filter_eps
      INTEGER(KIND=int_8), INTENT(OUT), OPTIONAL        :: flop
      REAL(KIND=real_8)                                  :: filter_eps_prv
      INTEGER(KIND=int_8), DIMENSION(2) :: dims_a, dims_b, dims_c
      INTEGER(KIND=int_8), DIMENSION(3) :: dims
      INTEGER, OPTIONAL, INTENT(IN) :: io_unit
      INTEGER :: dim_split
      CHARACTER(LEN=1) :: tr_case
      INTEGER :: data_type
      INTEGER :: numnodes
      INTEGER :: mp_comm, comm_tmp
      INTEGER, DIMENSION(2) :: pdims, pcoord
      CHARACTER(LEN=1) :: transa_prv, transb_prv, transc_prv
      TYPE(dbcsr_scalar_type) :: zero
      LOGICAL :: new_a, new_b, new_c
      TYPE(dbcsr_r_split_info) :: info
      INTEGER :: numproc, iproc
      CHARACTER(LEN=*), PARAMETER :: routineN = 'dbcsr_r_multiply', &
                                     routineP = moduleN//':'//routineN
      INTEGER :: handle, io_unit_prv, nsplit
      LOGICAL :: debug=.FALSE.

      IF (PRESENT(io_unit)) THEN
         io_unit_prv = io_unit
      ELSE
         io_unit_prv = 0
      ENDIF

      CALL timeset(routineN, handle)

      IF (.NOT. dbcsr_r_get_data_type(matrix_a) .EQ. dbcsr_r_get_data_type(matrix_b)) THEN
         DBCSR_ABORT("matrices must have same datatype")
      ENDIF

      data_type = dbcsr_r_get_data_type(matrix_a)

      transa_prv = transa; transb_prv = transb; transc_prv = transc

      dims_a = [dbcsr_r_nblkrows_total(matrix_a), dbcsr_r_nblkcols_total(matrix_a)]
      dims_b = [dbcsr_r_nblkrows_total(matrix_b), dbcsr_r_nblkcols_total(matrix_b)]
      dims_c = [dbcsr_r_nblkrows_total(matrix_c), dbcsr_r_nblkcols_total(matrix_c)]

      !fdims_a = [dbcsr_r_nfullrows_total(matrix_a), dbcsr_r_nfullcols_total(matrix_a)]
      !fdims_b = [dbcsr_r_nfullrows_total(matrix_b), dbcsr_r_nfullcols_total(matrix_b)]
      !fdims_c = [dbcsr_r_nfullrows_total(matrix_c), dbcsr_r_nfullcols_total(matrix_c)]

      IF (PRESENT(io_unit)) THEN
         IF (io_unit_prv .GT. 0) THEN
            WRITE (io_unit_prv, '(A)') repeat("-", 80)
            WRITE (io_unit_prv, '(A,1X,A,1X,A,1X,A,1X,A,1X,A)') "DBCSR RECTANGULAR MATRIX MULTIPLICATION:", &
               TRIM(matrix_a%matrix%name), 'x', TRIM(matrix_b%matrix%name), '=', TRIM(matrix_c%matrix%name)
            WRITE (io_unit_prv, '(A)') repeat("-", 80)
         ENDIF


         CALL dbcsr_r_write_matrix_info(matrix_a, io_unit_prv, full_info=debug)
         CALL dbcsr_r_write_matrix_info(matrix_b, io_unit_prv, full_info=debug)
         CALL dbcsr_r_write_matrix_info(matrix_c, io_unit_prv, full_info=debug)
         !WRITE(io_unit_prv, '(A,1X,A)') "Matrix A", TRIM(matrix_a%matrix%name)
         !WRITE(io_unit_prv, '(A,1X,A)') "Matrix B", TRIM(matrix_b%matrix%name)
         !WRITE(io_unit_prv, '(A,1X,A)') "Matrix C", TRIM(matrix_c%matrix%name)
         !WRITE (io_unit_prv, "(A, 1X, I8, 1X, I8, 1X, A, 1X, A)") &
         !   "block dims 1", dims_a(1), dims_a(2), "trans 1", transa_prv
         !WRITE (io_unit_prv, "(A, 1X, I8, 1X, I8)") &
         !   "full dims 1", fdims_a(1), fdims_a(2)
         !WRITE (io_unit_prv, "(A, 1X, I8, 1X, I8, 1X, A, 1X, A)") &
         !   "block dims 2", dims_b(1), dims_b(2), "trans 2", transb_prv
         !WRITE (io_unit_prv, "(A, 1X, I8, 1X, I8)") &
         !   "full dims 2", fdims_b(1), fdims_b(2)
         !WRITE (io_unit_prv, "(A, 1X, I8, 1X, I8, 1X, A, 1X, A)") &
         !   "block dims 3", dims_c(1), dims_c(2), "trans 3", transc_prv
         !WRITE (io_unit_prv, "(A, 1X, I8, 1X, I8)") &
         !   "full dims 3", fdims_c(1), fdims_c(2)
      ENDIF

      IF (transa_prv .EQ. dbcsr_transpose) THEN
         CALL swap(dims_a)
      ENDIF

      IF (transb_prv .EQ. dbcsr_transpose) THEN
         CALL swap(dims_b)
      ENDIF

      dims_c = [dims_a(1), dims_b(2)]

      IF (.NOT. (dims_a(2) .EQ. dims_b(1))) THEN
         DBCSR_ABORT("inconsistent matrix dimensions")
      ENDIF

      dims(:) = [dims_a(1), dims_a(2), dims_b(2)]

      dim_split = MAXLOC(dims, 1)

      tr_case = ''

      IF (io_unit_prv > 0) THEN
         WRITE (io_unit_prv, "(T2,A, 1X, I8, 1X, I8, 1X, I8)") "mm dims:", dims(1), dims(2), dims(3)
      ENDIF

      SELECT CASE (dim_split)
      CASE (1)

         CALL reshape_mm_compatible(matrix_a, matrix_c, matrix_a_rs, matrix_c_rs, &
            new_a, new_c, transa_prv, transc_prv, optimize_dist=optimize_dist, nodata2 = .TRUE., comm_new=comm_tmp, unit_nr=io_unit)

         CALL dbcsr_r_get_split_info(dbcsr_r_info(matrix_a_rs), mp_comm=mp_comm)

         ALLOCATE(matrix_b_rs)
         CALL reshape_mm_small(mp_comm, matrix_b, matrix_b_rs, transb_prv==dbcsr_transpose, transb_prv)
         new_b = .TRUE.
         tr_case = transa_prv

         IF (io_unit_prv > 0) THEN
            IF (tr_case == 'N') THEN
               WRITE (io_unit_prv, "(T2,A, 1X, A)") "mm case:", "| x + = |"
            ELSE
               WRITE (io_unit_prv, "(T2,A, 1X, A)") "mm case:", "--T x + = --T"
            ENDIF
         ENDIF

      CASE (2)

         CALL reshape_mm_compatible(matrix_a, matrix_b, matrix_a_rs, matrix_b_rs, new_a, new_b, transa_prv, transb_prv, &
                                    optimize_dist=optimize_dist, comm_new=comm_tmp, unit_nr=io_unit)

         CALL dbcsr_r_get_split_info(dbcsr_r_info(matrix_a_rs), mp_comm=mp_comm)

         ALLOCATE(matrix_c_rs)
         CALL reshape_mm_small(mp_comm, matrix_c, matrix_c_rs, transc_prv==dbcsr_transpose, transc_prv, nodata=.TRUE.)
         new_c = .TRUE.
         tr_case = transa_prv

         IF (io_unit_prv > 0) THEN
            IF (tr_case == 'N') THEN
               WRITE (io_unit_prv, "(T2,A, 1X, A)") "mm case:", "-- x --T = +"
            ELSE
               WRITE (io_unit_prv, "(T2,A, 1X, A)") "mm case:", "|T x | = +"
            ENDIF
         ENDIF

      CASE (3)

         CALL reshape_mm_compatible(matrix_b, matrix_c, matrix_b_rs, matrix_c_rs, new_b, new_c, transb_prv, &
                                    transc_prv, optimize_dist=optimize_dist, nodata2=.TRUE., comm_new=comm_tmp, unit_nr=io_unit)

         CALL dbcsr_r_get_split_info(dbcsr_r_info(matrix_b_rs), mp_comm=mp_comm)

         ALLOCATE(matrix_a_rs)
         CALL reshape_mm_small(mp_comm, matrix_a, matrix_a_rs, transa_prv==dbcsr_transpose, transa_prv)
         new_a = .TRUE.
         tr_case = transb_prv

         IF (io_unit_prv > 0) THEN
            IF (tr_case == 'N') THEN
               WRITE (io_unit_prv, "(T2,A, 1X, A)") "mm case:", "+ x -- = --"
            ELSE
               WRITE (io_unit_prv, "(T2,A, 1X, A)") "mm case:", "+ x |T = |T"
            ENDIF
         ENDIF

      END SELECT

      !IF (io_unit_prv > 0) THEN
      !   WRITE (io_unit_prv, "(A, 1X, A, 1X, A, 1X, A)") "new trans flags:", transa_prv, transb_prv, transc
      !ENDIF

      info = dbcsr_r_info(matrix_a_rs)
      CALL dbcsr_r_get_split_info(info, nsplit=nsplit, mp_comm=mp_comm)
      CALL mp_environ(numnodes, pdims, pcoord, mp_comm)

      IF(PRESENT(filter_eps)) THEN
         filter_eps_prv = filter_eps
         IF(dim_split == 2) THEN
            ! if result matrices are summed up then we need to correctly propagate filter_eps
            filter_eps_prv = REAL(filter_eps_prv, KIND=real_8)/SQRT(REAL(nsplit, KIND=real_8))
         ENDIF
      ELSE
         filter_eps_prv = 0.0_real_8
      ENDIF

      IF(PRESENT(io_unit)) CALL dbcsr_r_print_info(info, " Split info:", io_unit_prv)

      zero = dbcsr_scalar_zero(data_type)
      SELECT CASE (dim_split)
      CASE (1)
         CALL dbcsr_r_replicate(matrix_b_rs%matrix, dbcsr_r_info(matrix_a_rs), matrix_b_rep)
         SELECT CASE (tr_case)
         CASE (dbcsr_no_transpose)
            CALL mp_environ(numproc, iproc, mp_comm)
            CALL dbcsr_multiply(transa=dbcsr_no_transpose, transb=dbcsr_no_transpose, alpha=alpha, &
                                matrix_a=matrix_a_rs%matrix, matrix_b=matrix_b_rep%matrix, beta=beta, matrix_c=matrix_c_rs%matrix, &
                                filter_eps=filter_eps_prv, flop=flop)
         CASE (dbcsr_transpose)
            CALL dbcsr_multiply(transa=dbcsr_transpose, transb=dbcsr_no_transpose, alpha=alpha, &
                                matrix_a=matrix_b_rep%matrix, matrix_b=matrix_a_rs%matrix, beta=beta, matrix_c=matrix_c_rs%matrix, &
                                filter_eps=filter_eps_prv, flop=flop)
         END SELECT
         CALL dbcsr_r_distribution_destroy(matrix_b_rep%dist)
         CALL dbcsr_r_destroy(matrix_b_rep)
      CASE (2)
         CALL dbcsr_r_replicate(matrix_c_rs%matrix, dbcsr_r_info(matrix_a_rs), matrix_c_rep, nodata=.TRUE.)

         CALL dbcsr_multiply(transa=transa_prv, transb=transb_prv, alpha=alpha, matrix_a=matrix_a_rs%matrix, &
                             matrix_b=matrix_b_rs%matrix, beta=beta, matrix_c=matrix_c_rep%matrix, &
                             filter_eps=filter_eps_prv, flop=flop)

         CALL dbcsr_r_merge(matrix_c_rs%matrix, matrix_c_rep)
         CALL dbcsr_r_distribution_destroy(matrix_c_rep%dist)
         CALL dbcsr_r_destroy(matrix_c_rep)
      CASE (3)
         CALL dbcsr_r_replicate(matrix_a_rs%matrix, dbcsr_r_info(matrix_b_rs), matrix_a_rep)
         SELECT CASE (tr_case)
         CASE (dbcsr_no_transpose)
            CALL dbcsr_multiply(transa=dbcsr_no_transpose, transb=dbcsr_no_transpose, alpha=alpha, &
                                matrix_a=matrix_a_rep%matrix, matrix_b=matrix_b_rs%matrix, beta=beta, matrix_c=matrix_c_rs%matrix, &
                                filter_eps=filter_eps_prv, flop=flop)
         CASE (dbcsr_transpose)
            CALL dbcsr_multiply(transa=dbcsr_no_transpose, transb=dbcsr_transpose, alpha=alpha, &
                                matrix_a=matrix_b_rs%matrix, matrix_b=matrix_a_rep%matrix, beta=beta, matrix_c=matrix_c_rs%matrix, &
                                filter_eps=filter_eps_prv, flop=flop)
         END SELECT
         CALL dbcsr_r_distribution_destroy(matrix_a_rep%dist)
         CALL dbcsr_r_destroy(matrix_a_rep)
      END SELECT

      IF (new_a) THEN
         CALL dbcsr_r_destroy(matrix_a_rs)
         DEALLOCATE(matrix_a_rs)
      ENDIF
      IF (new_b) THEN
         CALL dbcsr_r_destroy(matrix_b_rs)
         DEALLOCATE(matrix_b_rs)
      ENDIF
      IF (new_c) THEN
         CALL dbcsr_scale(matrix_c%matrix, beta)
         CALL dbcsr_r_reshape(matrix_c_rs, matrix_c, summation=.TRUE., transposed=transc_prv/=transc)
         CALL dbcsr_r_destroy(matrix_c_rs)
         DEALLOCATE(matrix_c_rs)
      ENDIF

      IF(PRESENT(flop)) THEN
         CALL mp_sum(flop, mp_comm)
         flop = (flop+numnodes-1)/numnodes
      ENDIF

      IF(PRESENT(optimize_dist)) THEN
         IF(optimize_dist) CALL mp_comm_free(comm_tmp)
      ENDIF
      IF(io_unit_prv > 0) THEN
         WRITE (io_unit_prv, '(A)') repeat("-", 80)
         WRITE (io_unit_prv, '(A,1X,A,1X,A,1X,A,1X,A,1X,A)') "RECTANGULAR MATRIX MULTIPLICATION DONE"
         WRITE (io_unit_prv, '(A)') repeat("-", 80)
      ENDIF

      CALL timestop(handle)

   END SUBROUTINE

   SUBROUTINE reshape_mm_small(mp_comm, matrix_in, matrix_out, transposed, trans, nodata)
      INTEGER, INTENT(IN) :: mp_comm
      TYPE(dbcsr_r_type), INTENT(INOUT) :: matrix_in
      TYPE(dbcsr_r_type), INTENT(OUT) :: matrix_out
      LOGICAL, INTENT(IN) :: transposed
      CHARACTER(LEN=1), INTENT(OUT) :: trans
      LOGICAL, INTENT(IN), OPTIONAL :: nodata
      INTEGER :: numnodes
      INTEGER(KIND=int_8), DIMENSION(2) :: dims
      INTEGER, DIMENSION(2) :: pdims, pcoord
      TYPE(dbcsr_r_dist_cyclic) :: new_row_dist, new_col_dist
      TYPE(dbcsr_r_distribution_type) :: dist
      LOGICAL :: nodata_prv
      CHARACTER(LEN=*), PARAMETER :: routineN = 'reshape_mm_small', &
                                     routineP = moduleN//':'//routineN
      INTEGER :: handle

      CALL timeset(routineN, handle)

      IF(PRESENT(nodata)) THEN
         nodata_prv = nodata
      ELSE
         nodata_prv = .FALSE.
      ENDIF

      IF(transposed) THEN
         SELECT CASE(trans)
         CASE(dbcsr_transpose)
            trans = dbcsr_no_transpose
         CASE(dbcsr_no_transpose)
            trans = dbcsr_transpose
         END SELECT
      ENDIF

      CALL mp_environ(numnodes, pdims, pcoord, mp_comm)

      dims = [dbcsr_r_nblkrows_total(matrix_in), dbcsr_r_nblkcols_total(matrix_in)]

      IF (transposed) CALL swap(dims)

      new_row_dist = dbcsr_r_dist_cyclic(1, pdims(1), dims(1))
      new_col_dist = dbcsr_r_dist_cyclic(1, pdims(2), dims(2))

      CALL dbcsr_r_distribution_new(dist, mp_comm, new_row_dist, new_col_dist, nosplit=.TRUE.)
      IF (.NOT. transposed) THEN
         CALL dbcsr_r_create(matrix_out, matrix_in%matrix%name, dist, dbcsr_r_get_data_type(matrix_in), &
                             matrix_in%row_blk_size, matrix_in%col_blk_size, own_dist=.TRUE.)
      ELSE
         CALL dbcsr_r_create(matrix_out, matrix_in%matrix%name, dist, dbcsr_r_get_data_type(matrix_in), &
                             matrix_in%col_blk_size, matrix_in%row_blk_size, own_dist=.TRUE.)

      ENDIF
      IF(.NOT. nodata_prv) CALL dbcsr_r_reshape(matrix_in, matrix_out, transposed=transposed)

      CALL timestop(handle)

   END SUBROUTINE


   SUBROUTINE reshape_mm_compatible(matrix1_in, matrix2_in, matrix1_out, matrix2_out, new1, new2, trans1, trans2, &
                                    optimize_dist, nodata1, nodata2, comm_new, unit_nr)
      TYPE(dbcsr_r_type), TARGET, INTENT(INOUT) :: matrix1_in, matrix2_in
      TYPE(dbcsr_r_type), POINTER, INTENT(OUT) :: matrix1_out, matrix2_out
      LOGICAL, INTENT(OUT) :: new1, new2
      CHARACTER(LEN=1), INTENT(INOUT) :: trans1, trans2
      LOGICAL, INTENT(IN), OPTIONAL :: optimize_dist
      INTEGER, INTENT(OUT), OPTIONAL :: comm_new
      LOGICAL, OPTIONAL, INTENT(IN) :: nodata1, nodata2
      INTEGER, INTENT(IN), OPTIONAL :: unit_nr

      INTEGER(KIND=int_8), DIMENSION(2) :: dims1, dims2
      CHARACTER(LEN=*), PARAMETER :: routineN = 'reshape_mm_compatible', &
                                     routineP = moduleN//':'//routineN
      INTEGER :: handle, mp_comm, numnodes, split_rc_1, split_rc_2, io_unit
      INTEGER, DIMENSION(2) :: pcoord, pdims
      LOGICAL :: optimize_dist_prv, trans1_newdist, trans2_newdist
      TYPE(dbcsr_r_dist_cyclic) :: row_dist_1, col_dist_1, row_dist_2, col_dist_2
      TYPE(dbcsr_r_distribution_type) :: dist_1, dist_2
      TYPE(dbcsr_r_split_info) :: split_info

      IF(PRESENT(unit_nr)) THEN
         io_unit = unit_nr
      ELSE
         io_unit = 0
      ENDIF

      IF(PRESENT(optimize_dist)) THEN
         optimize_dist_prv = optimize_dist
      ELSE
         optimize_dist_prv = .FALSE.
      ENDIF

      IF(optimize_dist_prv) THEN
         DBCSR_ASSERT(PRESENT(comm_new))
      ENDIF

      CALL timeset(routineN, handle)
      new1 = .FALSE.; new2 = .FALSE.

      IF ((.NOT. optimize_dist_prv) .AND. dist_compatible(matrix1_in, matrix2_in)) THEN
         matrix1_out => matrix1_in
         matrix2_out => matrix2_in
         IF(io_unit > 0) THEN
            WRITE(io_unit, "(T2,A,1X,A,1X,A,1X,A)") "No redistribution of", TRIM(matrix1_in%matrix%name), &
               "and", TRIM(matrix2_in%matrix%name)
         ENDIF
      ELSE


         dims1 = [dbcsr_r_nblkrows_total(matrix1_in), dbcsr_r_nblkcols_total(matrix1_in)]
         dims2 = [dbcsr_r_nblkrows_total(matrix2_in), dbcsr_r_nblkcols_total(matrix2_in)]

         IF (optimize_dist_prv) THEN
            IF(io_unit > 0) THEN
               WRITE(io_unit, "(T2,A,1X,A,1X,A)") "Optimizing distribution of", TRIM(matrix1_in%matrix%name), &
                  "and", TRIM(matrix2_in%matrix%name)
            ENDIF
            CALL dbcsr_r_get_split_info(dbcsr_r_info(matrix1_in), mp_comm=mp_comm, split_rowcol=split_rc_1)
            CALL dbcsr_r_get_split_info(dbcsr_r_info(matrix2_in), split_rowcol=split_rc_2)

            trans1_newdist = (split_rc_1 == colsplit)
            trans2_newdist = (split_rc_2 == colsplit)

            IF(trans1_newdist) THEN
               CALL swap(dims1)
               CALL invert_transpose_flag(trans1)
            ENDIF

            IF(trans2_newdist) THEN
               CALL swap(dims2)
               CALL invert_transpose_flag(trans2)
            ENDIF

            IF (PRODUCT(dims1) .GE. PRODUCT(dims2)) THEN
               comm_new = dbcsr_r_mp_comm(mp_comm, dims1(1), dims1(2))
               CALL dbcsr_r_split(split_info, comm_new, dims1(1), dims1(2))
            ELSE
               comm_new = dbcsr_r_mp_comm(mp_comm, dims2(1), dims2(2))
               CALL dbcsr_r_split(split_info, comm_new, dims2(1), dims2(2))
            ENDIF

            CALL mp_environ(numnodes, pdims, pcoord, comm_new)

            ! todo: should be input whether new distribution should be cyclic or contiguous
            row_dist_1=dbcsr_r_dist_cyclic(1, pdims(1), dims1(1))
            col_dist_1=dbcsr_r_dist_cyclic(1, pdims(2), dims1(2))

            row_dist_2=dbcsr_r_dist_cyclic(1, pdims(1), dims2(1))
            col_dist_2=dbcsr_r_dist_cyclic(1, pdims(2), dims2(2))

            CALL dbcsr_r_distribution_new(dist_1, comm_new, row_dist_1, col_dist_1, split_info=split_info)
            CALL dbcsr_r_distribution_new(dist_2, comm_new, row_dist_2, col_dist_2, split_info=split_info)
            CALL dbcsr_r_release_info(split_info)

            ALLOCATE (matrix1_out)
            IF (.NOT. trans1_newdist) THEN
               CALL dbcsr_r_create(matrix1_out, matrix1_in%matrix%name, dist_1, dbcsr_r_get_data_type(matrix1_in), &
                                   matrix1_in%row_blk_size, matrix1_in%col_blk_size, own_dist=.TRUE.)

            ELSE
               CALL dbcsr_r_create(matrix1_out, matrix1_in%matrix%name, dist_1, dbcsr_r_get_data_type(matrix1_in), &
                                   matrix1_in%col_blk_size, matrix1_in%row_blk_size, own_dist=.TRUE.)
            ENDIF

            ALLOCATE (matrix2_out)
            IF (.NOT. trans2_newdist) THEN
               CALL dbcsr_r_create(matrix2_out, matrix2_in%matrix%name, dist_2, dbcsr_r_get_data_type(matrix2_in), &
                                   matrix2_in%row_blk_size, matrix2_in%col_blk_size, own_dist=.TRUE.)
            ELSE
               CALL dbcsr_r_create(matrix2_out, matrix2_in%matrix%name, dist_2, dbcsr_r_get_data_type(matrix2_in), &
                                   matrix2_in%col_blk_size, matrix2_in%row_blk_size, own_dist=.TRUE.)
            ENDIF

            CALL dbcsr_r_reshape(matrix1_in, matrix1_out, transposed=trans1_newdist)
            CALL dbcsr_r_reshape(matrix2_in, matrix2_out, transposed=trans2_newdist)
            new1 = .TRUE.
            new2 = .TRUE.

         ELSE
            IF (PRODUCT(dims1) .GE. PRODUCT(dims2)) THEN
               IF(io_unit > 0) THEN
                  WRITE(io_unit, "(T2,A,1X,A)") "Redistribution of", TRIM(matrix2_in%matrix%name)
               ENDIF

               ALLOCATE (matrix2_out)
               CALL reshape_mm_template(matrix1_in, matrix2_in, matrix2_out, trans2, nodata=nodata2)
               new2 = .TRUE.
               matrix1_out => matrix1_in
            ELSE
               IF(io_unit > 0) THEN
                  WRITE(io_unit, "(T2,A,1X,A)") "Redistribution of", TRIM(matrix1_in%matrix%name)
               ENDIF
               ALLOCATE (matrix1_out)
               CALL reshape_mm_template(matrix2_in, matrix1_in, matrix1_out, trans1, nodata=nodata1)
               new1 = .TRUE.
               matrix2_out => matrix2_in
            ENDIF
         ENDIF
      ENDIF

      CALL timestop(handle)

   CONTAINS

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
!> \param arr1 ...
!> \param arr2 ...
!> \return ...
! **************************************************************************************************
      PURE FUNCTION array_eq_i8(arr1, arr2)
         INTEGER(KIND=int_8), DIMENSION(:), INTENT(IN)      :: arr1, arr2
         LOGICAL                                            :: array_eq_i8

         array_eq_i8 = .FALSE.
         IF (SIZE(arr1) .EQ. SIZE(arr2)) array_eq_i8 = ALL(arr1 == arr2)

      END FUNCTION

! **************************************************************************************************
!> \brief ...
!> \param mat_a ...
!> \param mat_b ...
!> \return ...
! **************************************************************************************************
      FUNCTION dist_compatible(mat_a, mat_b, io_unit)
         TYPE(dbcsr_r_type), INTENT(IN)                     :: mat_a, mat_b
         INTEGER, INTENT(IN), OPTIONAL :: io_unit
         LOGICAL                                            :: dist_compatible

         INTEGER                                            :: res, same_local_rowcols
         TYPE(dbcsr_r_split_info)                           :: info_a, info_b
         INTEGER :: io_unit_prv, numproc, iproc

         IF (PRESENT(io_unit)) THEN
            io_unit_prv = io_unit
         ELSE
            io_unit_prv = 0
         ENDIF

         dist_compatible = .FALSE.

         info_a = dbcsr_r_info(mat_a)
         info_b = dbcsr_r_info(mat_b)

         CALL mp_environ(numproc, iproc, info_a%mp_comm)

         CALL mp_comm_compare(info_a%mp_comm_group, info_b%mp_comm_group, res)
         IF (res .GT. 1) THEN
            RETURN
         ENDIF

         IF (io_unit_prv > 0) THEN
            WRITE (io_unit_prv, *) "mp comm compatible"
         ENDIF

         IF (mat_a%dist%info%split_rowcol == mat_b%dist%info%split_rowcol) THEN

            IF (io_unit_prv > 0) THEN
               WRITE (io_unit_prv, *) "split compatible"
            ENDIF

            same_local_rowcols = MERGE(1, 0, array_eq_i8(mat_a%dist%local_rowcols, mat_b%dist%local_rowcols))
            CALL mp_sum(same_local_rowcols, info_a%mp_comm) ! todo: if inefficient, this can be done on master ranks only

            IF (same_local_rowcols == numproc) THEN
               IF (io_unit_prv > 0) THEN
                  WRITE (io_unit_prv, *) "local rowcols compatible"
               ENDIF
               dist_compatible = .TRUE.
            ELSE
               IF (io_unit_prv > 0) THEN
                  WRITE (io_unit_prv, *) "local rowcols A", mat_a%dist%local_rowcols
                  WRITE (io_unit_prv, *) "local rowcols B", mat_b%dist%local_rowcols
               ENDIF
            ENDIF
         ENDIF

         IF (io_unit_prv > 0) THEN
            WRITE (io_unit_prv, *) "is compatible?", dist_compatible
         ENDIF

      END FUNCTION

      SUBROUTINE reshape_mm_template(template, matrix_in, matrix_out, trans, nodata)
         TYPE(dbcsr_r_type), INTENT(IN) :: template
         TYPE(dbcsr_r_type), INTENT(INOUT) :: matrix_in
         TYPE(dbcsr_r_type), INTENT(OUT) :: matrix_out
         CHARACTER(LEN=1), INTENT(INOUT) :: trans
         LOGICAL, INTENT(IN), OPTIONAL :: nodata
         CLASS(dbcsr_r_distribution), ALLOCATABLE :: row_dist, col_dist

         TYPE(dbcsr_r_distribution_type) :: dist_new
         TYPE(dbcsr_r_split_info) :: info_template, info_matrix
         INTEGER :: mp_comm, dim_split_template, dim_split_matrix
         LOGICAL :: nodata_prv
         INTEGER, DIMENSION(2)                              :: pcoord, pdims
         INTEGER :: numnodes, handle
         LOGICAL :: transposed

         CHARACTER(LEN=*), PARAMETER :: routineN = 'reshape_mm_template', &
                                        routineP = moduleN//':'//routineN

         CALL timeset(routineN, handle)

         IF (PRESENT(nodata)) THEN
            nodata_prv = nodata
         ELSE
            nodata_prv = .FALSE.
         ENDIF

         info_template = dbcsr_r_info(template)
         info_matrix = dbcsr_r_info(matrix_in)

         dim_split_template = info_template%split_rowcol
         dim_split_matrix = info_matrix%split_rowcol

         transposed = dim_split_template .NE. dim_split_matrix
         IF(transposed) THEN
            SELECT CASE(trans)
            CASE(dbcsr_transpose)
               trans = dbcsr_no_transpose
            CASE(dbcsr_no_transpose)
               trans = dbcsr_transpose
            END SELECT
         ENDIF

         CALL mp_environ(numnodes, pdims, pcoord, info_template%mp_comm)

         SELECT CASE (dim_split_template)
         CASE (1)
            IF(.NOT. transposed) THEN
               ALLOCATE (row_dist, source=template%dist%row_dist)
               ALLOCATE (col_dist, source=dbcsr_r_dist_cyclic(1, pdims(2), matrix_in%nblkcols))
            ELSE
               ALLOCATE (row_dist, source=template%dist%row_dist)
               ALLOCATE (col_dist, source=dbcsr_r_dist_cyclic(1, pdims(2), matrix_in%nblkrows))
            ENDIF
         CASE (2)
            IF(.NOT. transposed) THEN
               ALLOCATE (row_dist, source=dbcsr_r_dist_cyclic(1, pdims(1), matrix_in%nblkrows))
               ALLOCATE (col_dist, source=template%dist%col_dist)
            ELSE
               ALLOCATE (row_dist, source=dbcsr_r_dist_cyclic(1, pdims(1), matrix_in%nblkcols))
               ALLOCATE (col_dist, source=template%dist%col_dist)
            ENDIF
         END SELECT

         CALL dbcsr_r_get_split_info(info_template, mp_comm=mp_comm)
         CALL dbcsr_r_distribution_new(dist_new, mp_comm, row_dist, col_dist, split_info=info_template)
         IF(.NOT. transposed) THEN
            CALL dbcsr_r_create(matrix_out, matrix_in%matrix%name, dist_new, dbcsr_r_get_data_type(matrix_in), &
                                matrix_in%row_blk_size, matrix_in%col_blk_size, own_dist=.TRUE.)
         ELSE
            CALL dbcsr_r_create(matrix_out, matrix_in%matrix%name, dist_new, dbcsr_r_get_data_type(matrix_in), &
                                matrix_in%col_blk_size, matrix_in%row_blk_size, own_dist=.TRUE.)
         ENDIF

         IF (.NOT. nodata_prv) CALL dbcsr_r_reshape(matrix_in, matrix_out, transposed=transposed)

         CALL timestop(handle)

      END SUBROUTINE

   END SUBROUTINE
END MODULE
