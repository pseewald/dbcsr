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
      dbcsr_r_nblkrows_total
   USE dbcsr_rectangular_base_types, ONLY: dbcsr_r_distribution_type, &
      dbcsr_r_split_info, &
      dbcsr_r_type
   USE dbcsr_rectangular_global, ONLY: dbcsr_r_dist_cyclic, dbcsr_r_distribution
   USE dbcsr_rectangular_reshape, ONLY: dbcsr_r_merge, &
      dbcsr_r_replicate, &
      dbcsr_r_reshape
   USE dbcsr_rectangular_split, ONLY: colsplit, &
      dbcsr_r_get_info, &
      dbcsr_r_split, &
      dbcsr_r_split_rows_or_cols, &
      rowsplit, &
      world_to_group_proc_map
   USE dbcsr_rectangular_util, ONLY: swap
   USE dbcsr_types, ONLY: dbcsr_no_transpose, &
      dbcsr_transpose
   USE dbcsr_kinds, ONLY: default_string_length, &
      int_8, &
      real_8
   USE dbcsr_mpiwrap, ONLY: mp_cart_create, &
      mp_comm_compare, &
      mp_environ
#include "../base/dbcsr_base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbcsr_rectangular_mm'

   PUBLIC :: dbcsr_r_multiply

CONTAINS

   ! variant for not necessarily square subgroups, probably not needed
   !SUBROUTINE dbcsr_r_split(split_info, mp_comm, nblkrows, nblkcols)
   !   TYPE(dbcsr_r_split_info), INTENT(OUT) :: split_info
   !   INTEGER, INTENT(IN) :: mp_comm
   !   INTEGER(KIND=int_8) :: nblkrows, nblkcols
   !   INTEGER(KIND=int_8), DIMENSION(2) :: mdims
   !   INTEGER :: split_rowcol
   !   INTEGER :: numproc, iproc
   !   REAL(KIND=real_8) :: s_ideal, s_diff
   !   INTEGER, DIMENSION(3) :: sort_indices
   !   INTEGER, DIMENSION(2) :: npdims, myploc
   !   INTEGER :: p1, p2, tmp
   !   INTEGER :: s, mp_comm_cart, pgrid_split_size, igroup
   !   IF(nblkrows .GE. nblkcols) THEN
   !      split_rowcol=rowsplit
   !   ELSE
   !      split_rowcol=colsplit
   !   ENDIF
   !   mdims(:) = [nblkrows, nblkcols]
   !   CALL mp_environ(numproc, iproc, mp_comm)
   !   s_ideal = REAL(MAXVAL(mdims), KIND=real_8)/REAL(MINVAL(mdims), KIND=real_8)
   !   DO p1 = 1, CEILING(SQRT(REAL(numproc, KIND=real_8)))
   !      IF(MOD(numproc, p1) .EQ. 0) THEN
   !         p2 = numproc/p1
   !         IF(ABS(REAL(p2, KIND=real_8)/REAL(p1, KIND=real_8) - s_ideal) .LT. s_diff) THEN
   !            s_diff = ABS(REAL(p2, KIND=real_8)/REAL(p1, KIND=real_8) - s_ideal)
   !            npdims(:) = [p2,p1]
   !         ENDIF
   !      ENDIF
   !   ENDDO
   !   DBCSR_ASSERT(npdims(1) .GT. npdims(2))
   !   s = CEILING(s_ideal)
   !   pgrid_split_size = (p2-1)/s+1
   !   SELECT CASE(split_rowcol)
   !   CASE(rowsplit)
   !   CASE(colsplit)
   !      tmp = npdims(2)
   !      npdims(2) = npdims(1)
   !      npdims(1) = tmp
   !   END SELECT
   !   CALL mp_cart_create(mp_comm, 2, npdims, myploc, mp_comm_cart)
   !   CALL world_to_group_proc_map(iproc, npdims, split_rowcol, pgrid_split_size, igroup)
   !   CALL dbcsr_r_split_rows_or_cols(split_info, mp_comm_cart, s, igroup, split_rowcol)
   !END SUBROUTINE
   ! create splitting based on matrix dimensions of largest matrix
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
   SUBROUTINE dbcsr_r_multiply(transa, transb, transc, alpha, matrix_a, matrix_b, matrix_c, io_unit)
      CHARACTER(LEN=1), INTENT(IN)                       :: transa, transb
      TYPE(dbcsr_scalar_type), INTENT(IN) :: alpha
      TYPE(dbcsr_r_type), INTENT(INOUT) :: matrix_a, matrix_b
      TYPE(dbcsr_r_type) :: matrix_b_rs, matrix_a_rs, matrix_c_rep, matrix_b_rep, matrix_a_rep
      TYPE(dbcsr_r_distribution_type) :: dist_c_rs
      CHARACTER(LEN=1), INTENT(OUT) :: transc
      TYPE(dbcsr_r_type), INTENT(OUT) :: matrix_c
      INTEGER(KIND=int_8), DIMENSION(2) :: dims_a, dims_b, dims_c
      INTEGER(KIND=int_8), DIMENSION(3) :: dims
      INTEGER, OPTIONAL, INTENT(IN) :: io_unit
      INTEGER :: dim_split
      CHARACTER(LEN=1) :: tr_case
      INTEGER :: data_type
      INTEGER :: numnodes
      INTEGER :: mp_comm
      INTEGER, DIMENSION(2) :: pdims, pcoord
      TYPE(dbcsr_r_dist_cyclic) :: new_row_dist_c, new_col_dist_c
      CHARACTER(LEN=1) :: transa_prv, transb_prv
      CHARACTER(LEN=default_string_length) :: name_c
      TYPE(dbcsr_scalar_type) :: zero
      LOGICAL :: new_a, new_b
      TYPE(dbcsr_r_split_info) :: info
      INTEGER :: numproc, iproc
      CHARACTER(LEN=*), PARAMETER :: routineN = 'dbcsr_r_multiply', &
                                     routineP = moduleN//':'//routineN
      INTEGER :: handle, io_unit_prv

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

      transa_prv = transa; transb_prv = transb
      new_a = .FALSE.; new_b = .FALSE.

      dims_a = [dbcsr_r_nblkrows_total(matrix_a), dbcsr_r_nblkcols_total(matrix_a)]
      dims_b = [dbcsr_r_nblkrows_total(matrix_b), dbcsr_r_nblkcols_total(matrix_b)]

      IF (io_unit_prv > 0) THEN
         WRITE (io_unit_prv, *) "DBCSR MM INFO | dims_a", dims_a, "transa ", transa_prv
         WRITE (io_unit_prv, *) "DBCSR MM INFO | dims_b", dims_b, "transb ", transb_prv
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
         WRITE (io_unit_prv, *) "DBCSR MM INFO | mm dims:", dims
      ENDIF

      SELECT CASE (dim_split)
      CASE (1)

         tr_case = transa_prv
         transc = transa_prv

         IF (io_unit_prv > 0) THEN
            IF (tr_case == 'N') THEN
               WRITE (io_unit_prv, *) "DBCSR MM INFO | mm case:", "| x o = |"
            ELSE
               WRITE (io_unit_prv, *) "DBCSR MM INFO | mm case:", "--T x o = --T"
            ENDIF
         ENDIF

         CALL dbcsr_r_get_info(dbcsr_r_info(matrix_a), mp_comm=mp_comm)

         CALL reshape_mm_small(mp_comm, matrix_b, matrix_b_rs, transb_prv==dbcsr_transpose)

         matrix_a_rs = matrix_a

         new_b = .TRUE.

      CASE (2)

         CALL reshape_mm_compatible(matrix_a, matrix_b, matrix_a_rs, matrix_b_rs, new_a, new_b, transa_prv, transb_prv)
         tr_case = transa_prv
         transc = dbcsr_no_transpose

         IF (io_unit_prv > 0) THEN
            IF (tr_case == 'N') THEN
               WRITE (io_unit_prv, *) "DBCSR MM INFO | mm case:", "-- x --T = o"
            ELSE
               WRITE (io_unit_prv, *) "DBCSR MM INFO | mm case:", "|T x | = o"
            ENDIF
         ENDIF

      CASE (3)
         tr_case = transb_prv
         transc = transb_prv

         IF (io_unit_prv > 0) THEN
            IF (tr_case == 'N') THEN
               WRITE (io_unit_prv, *) "DBCSR MM INFO | mm case:", "o x -- = --"
            ELSE
               WRITE (io_unit_prv, *) "DBCSR MM INFO | mm case:", "o x |T = |T"
            ENDIF
         ENDIF

         CALL dbcsr_r_get_info(dbcsr_r_info(matrix_b), mp_comm=mp_comm)

         CALL reshape_mm_small(mp_comm, matrix_a, matrix_a_rs, transa_prv==dbcsr_transpose)

         matrix_b_rs = matrix_b

         new_a = .TRUE.

      END SELECT

      IF (io_unit_prv > 0) THEN
         WRITE (io_unit_prv, *) "DBCSR MM INFO | new trans flags:", transa_prv, transb_prv, transc
      ENDIF

      info = dbcsr_r_info(matrix_a_rs)
      CALL dbcsr_r_get_info(info, mp_comm=mp_comm)
      CALL mp_environ(numnodes, pdims, pcoord, mp_comm)

      name_c = TRIM(matrix_a_rs%matrix%name)//" X "//TRIM(matrix_b_rs%matrix%name)
      SELECT CASE (dim_split)
      CASE (1)
         ! use split from a
         SELECT CASE (tr_case)
         CASE (dbcsr_transpose)
            new_row_dist_c = dbcsr_r_dist_cyclic(1, pdims(1), dims_c(2))
            CALL dbcsr_r_distribution_new(dist_c_rs, mp_comm, &
                                          new_row_dist_c, matrix_a_rs%dist%col_dist, &
                                          split_info=dbcsr_r_info(matrix_a_rs))
            CALL dbcsr_r_create(matrix_c, name_c, dist_c_rs, data_type, &
                                matrix_b_rs%col_blk_size, matrix_a_rs%col_blk_size, own_dist=.TRUE.)

         CASE (dbcsr_no_transpose)
            new_col_dist_c = dbcsr_r_dist_cyclic(1, pdims(2), dims_c(2))
            CALL dbcsr_r_distribution_new(dist_c_rs, mp_comm, &
                                          matrix_a_rs%dist%row_dist, new_col_dist_c, &
                                          split_info=dbcsr_r_info(matrix_a_rs))
            CALL dbcsr_r_create(matrix_c, name_c, dist_c_rs, data_type, &
                                matrix_a_rs%row_blk_size, matrix_b_rs%col_blk_size, own_dist=.TRUE.)

         END SELECT
      CASE (2)
         new_row_dist_c = dbcsr_r_dist_cyclic(1, pdims(1), dims_c(1))
         new_col_dist_c = dbcsr_r_dist_cyclic(1, pdims(2), dims_c(2))
         CALL dbcsr_r_distribution_new(dist_c_rs, mp_comm, new_row_dist_c, new_col_dist_c, nosplit=.TRUE.)
         SELECT CASE (tr_case)
         CASE (dbcsr_transpose)
            CALL dbcsr_r_create(matrix_c, name_c, dist_c_rs, data_type, &
                                matrix_a_rs%col_blk_size, matrix_b_rs%col_blk_size, own_dist=.TRUE.)
         CASE (dbcsr_no_transpose)
            CALL dbcsr_r_create(matrix_c, name_c, dist_c_rs, data_type, &
                                matrix_a_rs%row_blk_size, matrix_b_rs%row_blk_size, own_dist=.TRUE.)

         END SELECT
      CASE (3)
         ! use split from b
         SELECT CASE (tr_case)
         CASE (dbcsr_transpose)
            new_col_dist_c = dbcsr_r_dist_cyclic(1, pdims(2), dims_c(1))
            CALL dbcsr_r_distribution_new(dist_c_rs, mp_comm, &
                                          matrix_b_rs%dist%row_dist, new_col_dist_c, &
                                          split_info=dbcsr_r_info(matrix_b_rs))
            CALL dbcsr_r_create(matrix_c, name_c, dist_c_rs, data_type, &
                                matrix_b_rs%row_blk_size, matrix_a_rs%row_blk_size, own_dist=.TRUE.)
         CASE (dbcsr_no_transpose)
            new_row_dist_c = dbcsr_r_dist_cyclic(1, pdims(1), dims_c(1))
            CALL dbcsr_r_distribution_new(dist_c_rs, mp_comm, &
                                          new_row_dist_c, matrix_b_rs%dist%col_dist, &
                                          split_info=dbcsr_r_info(matrix_b_rs))
            CALL dbcsr_r_create(matrix_c, name_c, dist_c_rs, data_type, &
                                matrix_a_rs%row_blk_size, matrix_b_rs%col_blk_size, own_dist=.TRUE.)
         END SELECT
      END SELECT
      CALL dbcsr_r_finalize(matrix_c)

      zero = dbcsr_scalar_zero(data_type)
      SELECT CASE (dim_split)
      CASE (1)
         CALL dbcsr_r_replicate(matrix_b_rs%matrix, dbcsr_r_info(matrix_a_rs), matrix_b_rep)
         SELECT CASE (tr_case)
         CASE (dbcsr_no_transpose)
            CALL mp_environ(numproc, iproc, mp_comm)
            CALL dbcsr_multiply(transa=dbcsr_no_transpose, transb=dbcsr_no_transpose, alpha=alpha, &
                                matrix_a=matrix_a_rs%matrix, matrix_b=matrix_b_rep%matrix, beta=zero, matrix_c=matrix_c%matrix)
         CASE (dbcsr_transpose)
            CALL dbcsr_multiply(transa=dbcsr_transpose, transb=dbcsr_no_transpose, alpha=alpha, &
                                matrix_a=matrix_b_rep%matrix, matrix_b=matrix_a_rs%matrix, beta=zero, matrix_c=matrix_c%matrix)
         END SELECT
         CALL dbcsr_r_distribution_destroy(matrix_b_rep%dist)
         CALL dbcsr_r_destroy(matrix_b_rep)
      CASE (2)
         CALL dbcsr_r_replicate(matrix_c%matrix, dbcsr_r_info(matrix_a_rs), matrix_c_rep, nodata=.TRUE.)

         CALL dbcsr_multiply(transa=transa_prv, transb=transb_prv, alpha=alpha, matrix_a=matrix_a_rs%matrix, &
                             matrix_b=matrix_b_rs%matrix, beta=zero, matrix_c=matrix_c_rep%matrix)

         CALL dbcsr_r_merge(matrix_c%matrix, matrix_c_rep)
         CALL dbcsr_r_distribution_destroy(matrix_c_rep%dist)
         CALL dbcsr_r_destroy(matrix_c_rep)
      CASE (3)
         CALL dbcsr_r_replicate(matrix_a_rs%matrix, dbcsr_r_info(matrix_b_rs), matrix_a_rep)
         SELECT CASE (tr_case)
         CASE (dbcsr_no_transpose)
            CALL dbcsr_multiply(transa=dbcsr_no_transpose, transb=dbcsr_no_transpose, alpha=alpha, &
                                matrix_a=matrix_a_rep%matrix, matrix_b=matrix_b_rs%matrix, beta=zero, matrix_c=matrix_c%matrix)
         CASE (dbcsr_transpose)
            CALL dbcsr_multiply(transa=dbcsr_no_transpose, transb=dbcsr_transpose, alpha=alpha, &
                                matrix_a=matrix_b_rs%matrix, matrix_b=matrix_a_rep%matrix, beta=zero, matrix_c=matrix_c%matrix)
         END SELECT
         CALL dbcsr_r_distribution_destroy(matrix_a_rep%dist)
         CALL dbcsr_r_destroy(matrix_a_rep)
      END SELECT

      IF (new_a) THEN
         CALL dbcsr_r_destroy(matrix_a_rs)
      ENDIF
      IF (new_b) THEN
         CALL dbcsr_r_destroy(matrix_b_rs)
      ENDIF

      CALL timestop(handle)

   END SUBROUTINE

   SUBROUTINE reshape_mm_small(mp_comm, matrix_in, matrix_out, transposed, nodata)
      INTEGER, INTENT(IN) :: mp_comm
      TYPE(dbcsr_r_type), INTENT(INOUT) :: matrix_in
      TYPE(dbcsr_r_type), INTENT(OUT) :: matrix_out
      LOGICAL, INTENT(IN) :: transposed
      LOGICAL, INTENT(IN), OPTIONAL :: nodata
      INTEGER :: numnodes
      INTEGER(KIND=int_8), DIMENSION(2) :: dims
      INTEGER, DIMENSION(2) :: pdims, pcoord
      TYPE(dbcsr_r_dist_cyclic) :: new_row_dist, new_col_dist
      TYPE(dbcsr_r_distribution_type) :: dist
      LOGICAL :: nodata_prv

      IF(PRESENT(nodata)) THEN
         nodata_prv = nodata
      ELSE
         nodata_prv = .FALSE.
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

   END SUBROUTINE


   SUBROUTINE reshape_mm_compatible(matrix1_in, matrix2_in, matrix1_out, matrix2_out, new1, new2, trans1, trans2, nodata1, nodata2)
      TYPE(dbcsr_r_type), INTENT(INOUT) :: matrix1_in, matrix2_in
      TYPE(dbcsr_r_type), INTENT(OUT) :: matrix1_out, matrix2_out
      LOGICAL, INTENT(OUT) :: new1, new2
      CHARACTER(LEN=1), INTENT(INOUT) :: trans1, trans2
      LOGICAL, OPTIONAL, INTENT(IN) :: nodata1, nodata2

      INTEGER(KIND=int_8), DIMENSION(2) :: dims1, dims2

      new1 = .FALSE.; new2 = .FALSE.

      IF (dist_compatible(matrix1_in, matrix2_in)) THEN
         matrix1_out = matrix1_in
         matrix2_out = matrix2_in
      ELSE

         dims1 = [dbcsr_r_nblkrows_total(matrix1_in), dbcsr_r_nblkcols_total(matrix1_in)]
         dims2 = [dbcsr_r_nblkrows_total(matrix2_in), dbcsr_r_nblkcols_total(matrix2_in)]

         IF (PRODUCT(dims1) .GE. PRODUCT(dims2)) THEN
            CALL reshape_mm_template(matrix1_in, matrix2_in, matrix2_out, trans2, nodata=nodata2)
            new2 = .TRUE.
            matrix1_out = matrix1_in
         ELSE
            CALL reshape_mm_template(matrix2_in, matrix1_in, matrix1_out, trans1, nodata=nodata1)
            new1 = .TRUE.
            matrix2_out = matrix2_in
         ENDIF
      ENDIF

   CONTAINS

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

         INTEGER                                            :: res
         TYPE(dbcsr_r_split_info)                           :: info_a, info_b
         INTEGER :: io_unit_prv

         IF (PRESENT(io_unit)) THEN
            io_unit_prv = io_unit
         ELSE
            io_unit_prv = 0
         ENDIF

         dist_compatible = .FALSE.

         info_a = dbcsr_r_info(mat_a)
         info_b = dbcsr_r_info(mat_b)

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
            IF (array_eq_i8(mat_a%dist%local_rowcols, mat_b%dist%local_rowcols)) THEN
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
         INTEGER :: numnodes
         LOGICAL :: transposed

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

         !SELECT CASE (dim_split_template)
         !CASE (1)
         !   DBCSR_ASSERT(template%nblkrows == matrix_in%nblkrows)
         !CASE (2)
         !   DBCSR_ASSERT(template%nblkcols == matrix_in%nblkcols)
         !END SELECT

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
               !PRINT *, "matrix_in nblkcols", matrix_in%nblkcols, template%nblkrows
               ALLOCATE (row_dist, source=dbcsr_r_dist_cyclic(1, pdims(1), matrix_in%nblkcols))
               ALLOCATE (col_dist, source=template%dist%col_dist)
            ENDIF
         END SELECT

         CALL dbcsr_r_get_info(info_template, mp_comm=mp_comm)
         CALL dbcsr_r_distribution_new(dist_new, mp_comm, row_dist, col_dist, split_info=info_template)
         IF(.NOT. transposed) THEN ! todo: don't need to call create twice
            CALL dbcsr_r_create(matrix_out, matrix_in%matrix%name, dist_new, dbcsr_r_get_data_type(matrix_in), &
                                matrix_in%row_blk_size, matrix_in%col_blk_size, own_dist=.TRUE.)
         ELSE
            CALL dbcsr_r_create(matrix_out, matrix_in%matrix%name, dist_new, dbcsr_r_get_data_type(matrix_in), &
                                matrix_in%col_blk_size, matrix_in%row_blk_size, own_dist=.TRUE.)
         ENDIF

         IF (.NOT. nodata_prv) CALL dbcsr_r_reshape(matrix_in, matrix_out, transposed=transposed)

      END SUBROUTINE

   END SUBROUTINE

   ! recursive subroutine to get optimal grid size, needs more work!
   !SUBROUTINE dbcsr_r_get_pgrid(mdims, nproc, pdims)
   !   INTEGER(KIND=int_8), DIMENSION(3), INTENT(IN) :: mdims
   !   INTEGER(KIND=int_8), DIMENSION(3) :: mdims_r, mdims_s
   !   INTEGER, DIMENSION(3) :: sort_indices
   !   INTEGER, INTENT(IN) :: nproc
   !   INTEGER, DIMENSION(2), INTENT(OUT) :: pdims
   !   INTEGER, DIMENSION(2) :: psplit_12, psplit_23, psplit_13
   !   INTEGER(KIND=int_8), DIMENSION(2) :: mdims_12, mdims_23, mdims_13
   !   INTEGER :: which_split, s
   !   INTEGER :: nsteps, step
   !   !mdims_r(:) = mdims(:)
   !   mdims_s(:) = mdims(:)
   !   CALL sort(mdims_s, 3, sort_indices)
   !   DBCSR_ASSERT(ALL(mdims(sort_indices) == mdims_s))
   !   psplit_12 = [1,1]
   !   psplit_13 = [1,1]
   !   psplit_23 = [1,1]
   !   mdims_12 = [mdims_s(3), mdims_s(2)]
   !   mdims_13 = [mdims_s(3), mdims_s(1)]
   !   mdims_23 = [mdims_s(2), mdims_s(1)]
   !   !nsteps = 10
   !   step = 0
   !   !print*, "mdims_12", mdims_12
   !   ! do recursive
   !   nsteps = 2
   !   DO
   !      step = step+1
   !      mdims_r = [mdims_12(1), mdims_12(2), mdims_13(2)]
   !      mdims_s = mdims_r
   !      CALL sort(mdims_s, 3, sort_indices)
   !      which_split = sort_indices(3)
   !      s = INT((mdims_s(3) - 1_int_8)/mdims_s(2) + 1_int_8)
   !      s = MIN(s, nproc/PRODUCT(psplit_12))
   !      s = MAX(s, 2)
   !      print *, "mdims_s", mdims_s, REAL(mdims_s(3))/REAL(mdims_s(2))
   !      !print *, "s", s
   !      !PRINT *, psplit_12, "|", psplit_13, "|", psplit_23
   !      IF(PRODUCT(psplit_12)*s .GT. nproc) THEN
   !         EXIT
   !      ENDIF
   !      IF(mdims_s(3) ==1) EXIT
   !      IF(step>nsteps) EXIT
   !      SELECT CASE (which_split)
   !      CASE (1)
   !         mdims_12(1) = (mdims_12(1)-1)/s+1
   !         mdims_13(1) = (mdims_13(1)-1)/s+1
   !         psplit_12(1) = psplit_12(1)*s
   !         psplit_13(1) = psplit_13(1)*s
   !         IF (mdims_23(1) .GE. mdims_23(2)) THEN
   !            psplit_23 = [psplit_12(2), psplit_12(1)] ! match 2
   !         ELSE
   !            psplit_23 = [psplit_13(1), psplit_13(2)] ! match 3
   !         ENDIF
   !      CASE (2)
   !         mdims_12(2) = (mdims_12(2)-1)/s+1
   !         mdims_23(1) = (mdims_23(1)-1)/s+1
   !         psplit_12(2) = psplit_12(2)*s
   !         psplit_23(1) = psplit_23(1)*s
   !         IF (mdims_13(1) .GE. mdims_13(2)) THEN
   !            psplit_13 = [psplit_12(1), psplit_12(2)] ! match 1
   !         ELSE
   !            psplit_13 = [psplit_23(1), psplit_23(2)] ! match 3
   !         ENDIF
   !      CASE (3)
   !         mdims_13(2) = (mdims_13(2)-1)/s+1
   !         mdims_23(2) = (mdims_23(2)-1)/s+1
   !         psplit_13(2) = psplit_13(2)*s
   !         psplit_23(2) = psplit_23(2)*s
   !         IF (mdims_12(1) .GE. mdims_12(2)) THEN
   !            psplit_12 = [psplit_13(1), psplit_13(2)] ! match 1
   !         ELSE
   !            psplit_12 = [psplit_23(2), psplit_23(1)] ! match 2
   !         ENDIF
   !      END SELECT
   !      PRINT*, "pdims", [MAXVAL(psplit_12), MINVAL(psplit_12)], REAL(MAXVAL(psplit_12))/MINVAL(psplit_12)
   !   ENDDO
   !   pdims = [MAXVAL(psplit_12), MINVAL(psplit_12)]
   !   !CALL sort(mdims_s, 3 sort_indices)
   !END SUBROUTINE

END MODULE
