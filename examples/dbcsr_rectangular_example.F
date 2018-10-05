!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2018  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief   example for multiplying rectangular matrices
! **************************************************************************************************

PROGRAM dbcsr_rectangular_example

   USE dbcsr_api,                       ONLY: dbcsr_finalize_lib,&
                                              dbcsr_init_lib,&
                                              dbcsr_type_real_8
   USE dbcsr_data_methods,              ONLY: dbcsr_scalar
   USE dbcsr_rectangular_base,          ONLY: &
        dbcsr_r_create, dbcsr_r_destroy, dbcsr_r_distribution_destroy, dbcsr_r_distribution_new, &
        dbcsr_r_finalize, dbcsr_r_get_stored_coordinates, dbcsr_r_nblkcols_total, &
        dbcsr_r_nblkrows_total, dbcsr_r_put_block
   USE dbcsr_rectangular_base_types,    ONLY: dbcsr_r_distribution_type,&
                                              dbcsr_r_type
   USE dbcsr_rectangular_global,        ONLY: dbcsr_r_blk_size_arb,&
                                              dbcsr_r_dist_arb
   USE dbcsr_rectangular_mm,            ONLY: dbcsr_r_multiply
   USE dbcsr_rectangular_split,         ONLY: dbcsr_r_mp_comm
   USE dbcsr_kinds,                           ONLY: int_8,&
                                              real_8
   USE dbcsr_machine,                         ONLY: default_output_unit
   USE dbcsr_mpiwrap,                 ONLY: mp_comm_free,&
                                              mp_environ,&
                                              mp_world_finalize,&
                                              mp_world_init

   IMPLICIT NONE

   TYPE(dbcsr_r_type) :: matrix_A, matrix_B, matrix_C, matrix_MM
   TYPE(dbcsr_r_distribution_type) :: dist_AB, dist_C
   INTEGER :: mp_comm, mp_comm_cart
   INTEGER(KIND=int_8) :: dmax, dmin
   INTEGER(KIND=int_8), DIMENSION(2) :: dim_A, dim_B, dim_C
   INTEGER :: numnodes, mynode
   INTEGER, DIMENSION(2) :: pdims
   INTEGER, DIMENSION(:), ALLOCATABLE :: row_dist_AB, col_dist_AB, row_dist_C, col_dist_C, &
                                         dmin_blk_sizes, dmax_blk_sizes
   INTEGER :: io_unit, max_nze, nze
   INTEGER :: max_row_size, max_col_size
   REAL(KIND=real_8), DIMENSION(:, :), ALLOCATABLE :: values
   INTEGER(KIND=int_8) :: row, col, row_s, col_s
   INTEGER :: node_holds_blk
   LOGICAL :: tr
   CHARACTER(LEN=1) :: transc
   INTEGER, DIMENSION(2) :: pcoord

   TYPE(dbcsr_r_dist_arb), ALLOCATABLE :: row_dist_AB_obj, col_dist_AB_obj, row_dist_C_obj, col_dist_C_obj
   TYPE(dbcsr_r_blk_size_arb), ALLOCATABLE :: dmin_blk_size_obj, dmax_blk_size_obj

   dmin = 30
   dmax = 120

   dim_A = [dmax, dmin]
   dim_B = [dmax, dmin]
   dim_C = [dmin, dmin]

   CALL mp_world_init(mp_comm)
   CALL mp_environ(numnodes, mynode, mp_comm)

   CALL dbcsr_init_lib()

   io_unit = 0
   IF (mynode .EQ. 0) io_unit = default_output_unit

   mp_comm_cart = dbcsr_r_mp_comm(mp_comm, dmax, dmin)

   CALL mp_environ(numnodes, pdims, pcoord, mp_comm_cart)

   CALL random_dist(row_dist_AB, dmax, pdims(1))
   CALL random_dist(col_dist_AB, dmin, pdims(2))

   CALL random_dist(row_dist_C, dmin, pdims(1))
   CALL random_dist(col_dist_C, dmin, pdims(2))

   ALLOCATE (row_dist_AB_obj, source=dbcsr_r_dist_arb(row_dist_AB, pdims(1), dmax))
   ALLOCATE (col_dist_AB_obj, source=dbcsr_r_dist_arb(col_dist_AB, pdims(2), dmin))

   ALLOCATE (row_dist_C_obj, source=dbcsr_r_dist_arb(row_dist_C, pdims(1), dmin))
   ALLOCATE (col_dist_C_obj, source=dbcsr_r_dist_arb(col_dist_C, pdims(2), dmin))

   ALLOCATE (dmin_blk_sizes(dmin))
   ALLOCATE (dmax_blk_sizes(dmax))

   dmin_blk_sizes(:) = 2
   dmax_blk_sizes(:) = 3

   ALLOCATE (dmin_blk_size_obj, source=dbcsr_r_blk_size_arb(dmin_blk_sizes))
   ALLOCATE (dmax_blk_size_obj, source=dbcsr_r_blk_size_arb(dmax_blk_sizes))

   CALL dbcsr_r_distribution_new(dist_AB, mp_comm_cart, row_dist_AB_obj, col_dist_AB_obj)
   CALL dbcsr_r_distribution_new(dist_C, mp_comm_cart, row_dist_C_obj, col_dist_C_obj)

   CALL dbcsr_r_create(matrix_A, "this is my rectangular matrix A", &
                       dist=dist_AB, &
                       data_type=dbcsr_type_real_8, &
                       row_blk_size=dmax_blk_size_obj, &
                       col_blk_size=dmin_blk_size_obj)

   CALL dbcsr_r_create(matrix_B, "this is my rectangular matrix B", &
                       dist=dist_AB, &
                       data_type=dbcsr_type_real_8, &
                       row_blk_size=dmax_blk_size_obj, &
                       col_blk_size=dmin_blk_size_obj)

   CALL dbcsr_r_create(matrix_C, "this is my square matrix C", &
                       dist=dist_C, &
                       data_type=dbcsr_type_real_8, &
                       row_blk_size=dmin_blk_size_obj, &
                       col_blk_size=dmin_blk_size_obj)

   ! fill A
   max_row_size = 3
   max_col_size = 2
   max_nze = max_row_size*max_col_size

   ALLOCATE (values(max_row_size, max_col_size))

   DO row = 1, dbcsr_r_nblkrows_total(matrix_A)
      DO col = 1, MIN(row+1_int_8, dbcsr_r_nblkcols_total(matrix_A))

         tr = .FALSE.
         row_s = row; col_s = col
         CALL dbcsr_r_get_stored_coordinates(matrix_A, row_s, col_s, node_holds_blk)

         IF (node_holds_blk .EQ. mynode) THEN
            nze = dmax_blk_size_obj%data(row_s)*dmin_blk_size_obj%data(col_s)
            CALL RANDOM_NUMBER(values)
            CALL dbcsr_r_put_block(matrix_A, row_s, col_s, &
                                   values(1:dmax_blk_size_obj%data(row_s), &
                                          1:dmin_blk_size_obj%data(col_s)))
         ENDIF

      ENDDO
   ENDDO
   DEALLOCATE (values)
   CALL dbcsr_r_finalize(matrix_A)

   ! fill B
   max_row_size = 3
   max_col_size = 2
   max_nze = max_row_size*max_col_size
   ALLOCATE (values(max_row_size, max_col_size))

   DO row = 1, dbcsr_r_nblkrows_total(matrix_B)
      DO col = 1, MIN(row+1_int_8, dbcsr_r_nblkcols_total(matrix_B))

         tr = .FALSE.
         row_s = row; col_s = col
         CALL dbcsr_r_get_stored_coordinates(matrix_B, row_s, col_s, node_holds_blk)

         IF (node_holds_blk .EQ. mynode) THEN
            nze = dmax_blk_size_obj%data(row_s)*dmin_blk_size_obj%data(col_s)
            CALL RANDOM_NUMBER(values)
            CALL dbcsr_r_put_block(matrix_B, row_s, col_s, &
                                   values(1:dmax_blk_size_obj%data(row_s), &
                                          1:dmin_blk_size_obj%data(col_s)))
         ENDIF

      ENDDO
   ENDDO
   DEALLOCATE (values)
   CALL dbcsr_r_finalize(matrix_B)

   ! fill C
   max_row_size = 2
   max_col_size = 2
   max_nze = max_row_size*max_col_size
   ALLOCATE (values(max_row_size, max_col_size))
   DO row = 1, dbcsr_r_nblkrows_total(matrix_C)
      DO col = 1, MIN(row+1_int_8, dbcsr_r_nblkcols_total(matrix_C))

         tr = .FALSE.
         row_s = row; col_s = col
         CALL dbcsr_r_get_stored_coordinates(matrix_C, row_s, col_s, node_holds_blk)

         IF (node_holds_blk .EQ. mynode) THEN
            nze = dmin_blk_size_obj%data(row_s)*dmin_blk_size_obj%data(col_s)
            CALL RANDOM_NUMBER(values)
            CALL dbcsr_r_put_block(matrix_C, row_s, col_s, &
                                   values(1:dmin_blk_size_obj%data(row_s), &
                                          1:dmin_blk_size_obj%data(col_s)))
         ENDIF

      ENDDO
   ENDDO
   DEALLOCATE (values)
   CALL dbcsr_r_finalize(matrix_C)

   CALL dbcsr_r_multiply('T', 'N', transc, dbcsr_scalar(1.0_real_8), matrix_A, matrix_B, matrix_MM)

   CALL dbcsr_r_destroy(matrix_A)
   CALL dbcsr_r_destroy(matrix_B)
   CALL dbcsr_r_destroy(matrix_C)
   CALL dbcsr_r_destroy(matrix_MM)

   CALL dbcsr_r_distribution_destroy(dist_AB)
   CALL dbcsr_r_distribution_destroy(dist_C)

   CALL mp_comm_free(mp_comm_cart)

   ! finalize libdbcsr
   CALL dbcsr_finalize_lib(mp_comm, io_unit)

   !
   ! finalize mpi
   CALL mp_world_finalize()

   DEALLOCATE (row_dist_AB)
   DEALLOCATE (row_dist_C)
   DEALLOCATE (col_dist_AB)
   DEALLOCATE (col_dist_C)

   DEALLOCATE (dmin_blk_sizes)
   DEALLOCATE (dmax_blk_sizes)

   DEALLOCATE (row_dist_AB_obj, col_dist_AB_obj, row_dist_C_obj, col_dist_C_obj, dmin_blk_size_obj, dmax_blk_size_obj)

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param dist_array ...
!> \param dist_size ...
!> \param nbins ...
! **************************************************************************************************
   SUBROUTINE random_dist(dist_array, dist_size, nbins)
      INTEGER, ALLOCATABLE, DIMENSION(:), INTENT(out)    :: dist_array
      INTEGER(KIND=int_8), INTENT(in)                    :: dist_size
      INTEGER, INTENT(IN)                                :: nbins

      INTEGER(KIND=int_8)                                :: i

      ALLOCATE (dist_array(dist_size))
      DO i = 1, dist_size
         dist_array(i) = INT(MODULO(INT(nbins, KIND=int_8)-i, INT(nbins, KIND=int_8)))
      END DO

   END SUBROUTINE random_dist

END PROGRAM
