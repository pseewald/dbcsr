!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2018  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief routines to reshape/replicate/merge rectangular matrices.
!> \author Patrick Seewald
! **************************************************************************************************
MODULE dbcsr_rectangular_reshape

   USE dbcsr_block_access,              ONLY: dbcsr_put_block,&
                                              dbcsr_reserve_blocks
   USE dbcsr_data_methods,              ONLY: dbcsr_data_clear_pointer,&
                                              dbcsr_data_init,&
                                              dbcsr_data_new,&
                                              dbcsr_data_release,&
                                              dbcsr_type_1d_to_2d
   USE dbcsr_data_methods_low,          ONLY: internal_data_allocate,&
                                              internal_data_deallocate
   USE dbcsr_data_types,                ONLY: dbcsr_data_obj, dbcsr_scalar_type
   USE dbcsr_dist_methods,              ONLY: dbcsr_distribution_col_dist,&
                                              dbcsr_distribution_row_dist
   USE dbcsr_dist_operations,           ONLY: dbcsr_get_stored_coordinates
   USE dbcsr_index_operations,          ONLY: dbcsr_sort_indices
   USE dbcsr_iterator_operations,       ONLY: dbcsr_iterator_blocks_left,&
                                              dbcsr_iterator_next_block,&
                                              dbcsr_iterator_start,&
                                              dbcsr_iterator_stop
   USE dbcsr_methods,                   ONLY: dbcsr_blk_column_size,&
                                              dbcsr_blk_row_size
   USE dbcsr_operations,                ONLY: dbcsr_get_info, dbcsr_clear
   USE dbcsr_rectangular_base,          ONLY: &
        dbcsr_r_blk_sizes, dbcsr_r_create, dbcsr_r_distribution_new, dbcsr_r_finalize, &
        dbcsr_r_get_data_type, dbcsr_r_get_stored_coordinates, dbcsr_r_info, &
        dbcsr_r_iterator_blocks_left, dbcsr_r_iterator_next_block, dbcsr_r_iterator_start, &
        dbcsr_r_iterator_stop, dbcsr_r_put_block, dbcsr_r_reserve_blocks, &
        dbcsr_repl_get_stored_coordinates
   USE dbcsr_rectangular_base_types,    ONLY: dbcsr_r_distribution_type,&
                                              dbcsr_r_iterator,&
                                              dbcsr_r_split_info,&
                                              dbcsr_r_type
   USE dbcsr_rectangular_global,        ONLY: dbcsr_r_blk_size_arb,&
                                              dbcsr_r_blk_size_repl,&
                                              dbcsr_r_dist_arb,&
                                              dbcsr_r_dist_repl,&
                                              dbcsr_r_distribution,&
                                              dbcsr_r_rowcol_data
   USE dbcsr_rectangular_reshape_low,   ONLY: block_buffer_add_block,&
                                              block_buffer_create,&
                                              block_buffer_destroy,&
                                              block_buffer_get_next_block,&
                                              block_buffer_type,&
                                              dbcsr_r_communicate_buffer
   USE dbcsr_rectangular_split,         ONLY: colsplit,&
                                              dbcsr_r_get_split_info,&
                                              rowsplit
   USE dbcsr_types,                     ONLY: dbcsr_distribution_obj,&
                                              dbcsr_iterator,&
                                              dbcsr_type
   USE dbcsr_work_operations,           ONLY: dbcsr_finalize
   USE dbcsr_kinds,                           ONLY: int_8
   USE dbcsr_mpiwrap,                 ONLY: mp_alltoall,&
                                              mp_environ
#include "../base/dbcsr_base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbcsr_rectangular_reshape'

   PUBLIC :: dbcsr_r_reshape, dbcsr_r_merge, dbcsr_r_replicate

CONTAINS

! **************************************************************************************************
!> \brief get all indices from buffer
!> \param buffer ...
!> \param index ...
! **************************************************************************************************
   SUBROUTINE block_buffer_get_index(buffer, index)
      TYPE(block_buffer_type), INTENT(IN)                :: buffer
      INTEGER(KIND=int_8), ALLOCATABLE, &
         DIMENSION(:, :), INTENT(OUT)                    :: index

      INTEGER, DIMENSION(2)                              :: indx_shape
      CHARACTER(LEN=*), PARAMETER :: routineN = 'block_buffer_get_index', &
         routineP = moduleN//':'//routineN
      INTEGER :: handle

      CALL timeset(routineN, handle)

      indx_shape = SHAPE(buffer%indx)-[0, 1]
      ALLOCATE (INDEX(indx_shape(1), indx_shape(2)))
      INDEX(:, :) = buffer%indx(1:indx_shape(1), 1:indx_shape(2))
      CALL timestop(handle)
   END SUBROUTINE

! **************************************************************************************************
!> \brief Reset buffer block iterator
!> \param buffer ...
! **************************************************************************************************
   SUBROUTINE block_buffer_iterator_reset(buffer)
      TYPE(block_buffer_type), INTENT(INOUT)             :: buffer

      buffer%endpos = 0
   END SUBROUTINE

! **************************************************************************************************
!> \brief how many blocks left in iterator
!> \param buffer ...
!> \return ...
! **************************************************************************************************
   PURE FUNCTION block_buffer_blocks_left(buffer)
      TYPE(block_buffer_type), INTENT(IN)                :: buffer
      LOGICAL                                            :: block_buffer_blocks_left

      block_buffer_blocks_left = buffer%endpos .LT. buffer%nblock
   END FUNCTION

! **************************************************************************************************
!> \brief copy data (involves reshape)
!> \param matrix_in ...
!> \param matrix_out ...
!> \param transposed ...
! **************************************************************************************************
   RECURSIVE SUBROUTINE dbcsr_r_reshape(matrix_in, matrix_out, summation, transposed)
      TYPE(dbcsr_r_type), INTENT(INOUT)                  :: matrix_in, matrix_out
      LOGICAL, INTENT(IN), OPTIONAL                      :: summation
      LOGICAL, INTENT(IN), OPTIONAL                      :: transposed

      CHARACTER(LEN=*), PARAMETER :: routineN = 'dbcsr_r_reshape', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, handle2, iproc, mp_comm, mynode, nblock, &
                                                            ndata, numnodes, bcount, nblk
      INTEGER(KIND=int_8), ALLOCATABLE, DIMENSION(:, :)  :: index_recv, blks_to_allocate
      INTEGER(KIND=int_8), DIMENSION(2)                  :: blk_index
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: num_blocks_recv, num_blocks_send, &
                                                            num_entries_recv, num_entries_send, &
                                                            num_rec, num_send
      INTEGER, ALLOCATABLE, DIMENSION(:, :)              :: req_array
      INTEGER, DIMENSION(2)                              :: blk_size
      LOGICAL                                            :: tr, tr_in
      TYPE(block_buffer_type), ALLOCATABLE, DIMENSION(:) :: buffer_recv, buffer_send
      TYPE(dbcsr_data_obj)                               :: block
      TYPE(dbcsr_r_iterator)                             :: iter
      TYPE(dbcsr_r_split_info)                           :: info

      CALL timeset(routineN, handle)

      IF(PRESENT(summation)) THEN
         IF(.NOT. summation) CALL dbcsr_clear(matrix_out%matrix)
      ELSE
         CALL dbcsr_clear(matrix_out%matrix)
      ENDIF

      IF (PRESENT(transposed)) THEN
         tr_in = transposed
      ELSE
         tr_in = .FALSE.
      ENDIF

      IF (.NOT. matrix_out%valid) THEN
         DBCSR_ABORT("can not reshape into invalid matrix")
      ENDIF

      info = dbcsr_r_info(matrix_in)
      mp_comm = info%mp_comm
      CALL mp_environ(numnodes, mynode, mp_comm)
      ALLOCATE (buffer_send(0:numnodes-1))
      ALLOCATE (buffer_recv(0:numnodes-1))
      ALLOCATE (num_blocks_recv(0:numnodes-1))
      ALLOCATE (num_blocks_send(0:numnodes-1))
      ALLOCATE (num_entries_recv(0:numnodes-1))
      ALLOCATE (num_entries_send(0:numnodes-1))
      ALLOCATE (num_rec(0:2*numnodes-1))
      ALLOCATE (num_send(0:2*numnodes-1))
      num_send(:) = 0
      ALLOCATE (req_array(1:numnodes, 4))
      CALL dbcsr_r_iterator_start(iter, matrix_in)

      CALL timeset(routineN//"_get_coord", handle2)
      DO WHILE (dbcsr_r_iterator_blocks_left(iter))
         CALL dbcsr_r_iterator_next_block(iter, blk_index(1), blk_index(2), nblock, transposed=tr, &
                                          row_size=blk_size(1), col_size=blk_size(2))

         IF (tr_in) THEN
            CALL dbcsr_r_get_stored_coordinates(matrix_out, blk_index(2), blk_index(1), iproc)
         ELSE
            CALL dbcsr_r_get_stored_coordinates(matrix_out, blk_index(1), blk_index(2), iproc)
         ENDIF

         num_send(2*iproc) = num_send(2*iproc)+PRODUCT(blk_size)
         num_send(2*iproc+1) = num_send(2*iproc+1)+1
      ENDDO
      CALL dbcsr_r_iterator_stop(iter)
      CALL timestop(handle2)

      CALL timeset(routineN//"_alltoall", handle2)
      CALL mp_alltoall(num_send, num_rec, 2, mp_comm)
      CALL timestop(handle2)


      CALL timeset(routineN//"_buffer_fill", handle2)
      DO iproc = 0, numnodes-1
         num_entries_recv(iproc) = num_rec(2*iproc)
         num_blocks_recv(iproc) = num_rec(2*iproc+1)
         num_entries_send(iproc) = num_send(2*iproc)
         num_blocks_send(iproc) = num_send(2*iproc+1)

         CALL block_buffer_create(buffer_send(iproc), num_blocks_send(iproc), num_entries_send(iproc), &
                                  dbcsr_r_get_data_type(matrix_in))

         CALL block_buffer_create(buffer_recv(iproc), num_blocks_recv(iproc), num_entries_recv(iproc), &
                                  dbcsr_r_get_data_type(matrix_in))

      ENDDO

      CALL dbcsr_data_init(block)
      CALL dbcsr_data_new(block, dbcsr_type_1d_to_2d(dbcsr_r_get_data_type(matrix_in))) ! need to convert to 2d data type
      CALL dbcsr_r_iterator_start(iter, matrix_in)
      DO WHILE (dbcsr_r_iterator_blocks_left(iter))
         CALL dbcsr_r_iterator_next_block(iter, blk_index(1), blk_index(2), block, tr, &
                                          row_size=blk_size(1), col_size=blk_size(2))

         DBCSR_ASSERT(tr .EQV. .FALSE.)

         IF (tr_in) THEN
            CALL dbcsr_r_get_stored_coordinates(matrix_out, blk_index(2), blk_index(1), iproc)
         ELSE
            CALL dbcsr_r_get_stored_coordinates(matrix_out, blk_index(1), blk_index(2), iproc)
         ENDIF
         CALL block_buffer_add_block(buffer_send(iproc), blk_index, block, transposed=tr_in)
      ENDDO
      CALL dbcsr_r_iterator_stop(iter)
      CALL dbcsr_data_clear_pointer(block)

      CALL timestop(handle2)

      CALL timeset(routineN//"_communicate_buffer", handle2)
      CALL dbcsr_r_communicate_buffer(mp_comm, buffer_recv, buffer_send, req_array)
      CALL timestop(handle2)

      CALL timeset(routineN//"_buffer_obtain", handle2)

      nblk = SUM(num_blocks_recv)
      ALLOCATE (blks_to_allocate(nblk, 2))

      bcount = 0
      DO iproc = 0, numnodes - 1
         CALL block_buffer_get_index(buffer_recv(iproc), index_recv)
         blks_to_allocate(bcount + 1:bcount + SIZE(index_recv, 1), :) = index_recv(:, :)
         bcount = bcount + SIZE(index_recv, 1)
         DEALLOCATE(index_recv)
      ENDDO

      CALL dbcsr_r_reserve_blocks(matrix_out, blks_to_allocate(:, 1), blks_to_allocate(:, 2))
      DEALLOCATE(blks_to_allocate)

      DO iproc = 0, numnodes-1
         ! First, we need to get the index to create block
         DO WHILE (block_buffer_blocks_left(buffer_recv(iproc)))
            CALL block_buffer_get_next_block(buffer_recv(iproc), ndata, blk_index)
            CALL dbcsr_r_blk_sizes(matrix_out, blk_index(1), blk_index(2), blk_size(1), blk_size(2))
            CALL internal_data_allocate(block%d, blk_size)
            CALL block_buffer_get_next_block(buffer_recv(iproc), ndata, blk_index, block)
            CALL dbcsr_r_put_block(matrix_out, blk_index(1), blk_index(2), block, summation=summation)
            CALL internal_data_deallocate(block%d)
         ENDDO
         CALL block_buffer_destroy(buffer_recv(iproc))
         CALL block_buffer_destroy(buffer_send(iproc))
      ENDDO
      CALL dbcsr_data_clear_pointer(block)
      CALL dbcsr_data_release(block)

      CALL timestop(handle2)


      CALL dbcsr_r_finalize(matrix_out)

      CALL timestop(handle)
   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param matrix_in ...
!> \param info ...
!> \param matrix_out ...
!> \param nodata ...
! **************************************************************************************************
   SUBROUTINE dbcsr_r_replicate(matrix_in, info, matrix_out, nodata)
      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix_in
      TYPE(dbcsr_r_split_info), INTENT(IN)               :: info
      TYPE(dbcsr_r_type), INTENT(OUT)                    :: matrix_out
      LOGICAL, INTENT(IN), OPTIONAL                      :: nodata

      INTEGER                                            :: data_type, nblkcols, nblkrows
      INTEGER, DIMENSION(2)                              :: pcoord, pdims
      INTEGER, DIMENSION(:), POINTER                     :: col_blk_size, col_dist, row_blk_size, &
                                                            row_dist
      TYPE(dbcsr_distribution_obj)                       :: dbcsr_dist
      TYPE(dbcsr_r_dist_arb), TARGET                     :: dir_dist
      TYPE(dbcsr_r_dist_repl), TARGET                    :: repl_dist

      CLASS(dbcsr_r_distribution), ALLOCATABLE :: col_dist_obj, row_dist_obj
      CLASS(dbcsr_r_rowcol_data), ALLOCATABLE :: row_bsize_obj, col_bsize_obj
      TYPE(dbcsr_r_blk_size_repl), TARGET :: repl_blksize
      TYPE(dbcsr_r_blk_size_arb), TARGET :: dir_blksize
      TYPE(dbcsr_r_distribution_type) :: dist
      INTEGER :: mp_comm, numnodes, mynode
      TYPE(block_buffer_type), ALLOCATABLE, DIMENSION(:) :: buffer_recv, buffer_send
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: num_blocks_recv, num_blocks_send, &
                                                            num_entries_recv, num_entries_send, &
                                                            num_rec, num_send
      INTEGER, ALLOCATABLE, DIMENSION(:, :) :: req_array, blks_to_allocate
      INTEGER, DIMENSION(2) :: blk_size
      INTEGER, DIMENSION(2) :: blk_index
      INTEGER(KIND=int_8), DIMENSION(2) :: blk_index_i8
      TYPE(dbcsr_iterator) :: iter
      INTEGER :: nblock, i, iproc, bcount, nblk
      INTEGER, DIMENSION(:), ALLOCATABLE :: iprocs
      LOGICAL :: tr, nodata_prv
      INTEGER(KIND=int_8), ALLOCATABLE, DIMENSION(:, :) :: index_recv
      INTEGER :: ndata

      TYPE(dbcsr_data_obj) :: block

      CHARACTER(LEN=*), PARAMETER :: routineN = 'dbcsr_r_replicate', &
                                     routineP = moduleN//':'//routineN

      INTEGER :: handle, handle2

      CALL timeset(routineN, handle)

      IF (PRESENT(nodata)) THEN
         nodata_prv = nodata
      ELSE
         nodata_prv = .FALSE.
      ENDIF

      CALL dbcsr_get_info(matrix_in, distribution=dbcsr_dist, data_type=data_type, &
                          nblkrows_total=nblkrows, nblkcols_total=nblkcols, &
                          row_blk_size=row_blk_size, col_blk_size=col_blk_size)
      row_dist => dbcsr_distribution_row_dist(dbcsr_dist)
      col_dist => dbcsr_distribution_col_dist(dbcsr_dist)

      mp_comm = info%mp_comm

      CALL mp_environ(numnodes, mynode, mp_comm)
      CALL mp_environ(numnodes, pdims, pcoord, mp_comm)

      SELECT CASE (info%split_rowcol)
      CASE (rowsplit)
         repl_dist = dbcsr_r_dist_repl(row_dist, pdims(1), nblkrows, info%ngroup, info%pgrid_split_size, info%pgrid_split_size_last)
         dir_dist = dbcsr_r_dist_arb(col_dist, pdims(2), INT(nblkcols, KIND=int_8))
         repl_blksize = dbcsr_r_blk_size_repl(row_blk_size, info%ngroup)
         dir_blksize = dbcsr_r_blk_size_arb(col_blk_size)
         ALLOCATE (row_dist_obj, source=repl_dist)
         ALLOCATE (col_dist_obj, source=dir_dist)
         ALLOCATE (row_bsize_obj, source=repl_blksize)
         ALLOCATE (col_bsize_obj, source=dir_blksize)
      CASE (colsplit)
         dir_dist = dbcsr_r_dist_arb(row_dist, pdims(1), INT(nblkrows, KIND=int_8))
         repl_dist = dbcsr_r_dist_repl(col_dist, pdims(2), nblkcols, info%ngroup, info%pgrid_split_size, info%pgrid_split_size_last)
         dir_blksize = dbcsr_r_blk_size_arb(row_blk_size)
         repl_blksize = dbcsr_r_blk_size_repl(col_blk_size, info%ngroup)
         ALLOCATE (row_dist_obj, source=dir_dist)
         ALLOCATE (col_dist_obj, source=repl_dist)
         ALLOCATE (row_bsize_obj, source=dir_blksize)
         ALLOCATE (col_bsize_obj, source=repl_blksize)
      END SELECT

      CALL dbcsr_r_distribution_new(dist, mp_comm, row_dist_obj, col_dist_obj, split_info=info)
      CALL dbcsr_r_create(matrix_out, TRIM(matrix_in%name)//"_replicated", &
                          dist, data_type, row_bsize_obj, col_bsize_obj)

      IF (nodata_prv) THEN
         CALL dbcsr_r_finalize(matrix_out)
         CALL timestop(handle)
         RETURN
      ENDIF

      ALLOCATE (buffer_send(0:numnodes-1))
      ALLOCATE (buffer_recv(0:numnodes-1))
      ALLOCATE (num_blocks_recv(0:numnodes-1))
      ALLOCATE (num_blocks_send(0:numnodes-1))
      ALLOCATE (num_entries_recv(0:numnodes-1))
      ALLOCATE (num_entries_send(0:numnodes-1))
      ALLOCATE (num_rec(0:2*numnodes-1))
      ALLOCATE (num_send(0:2*numnodes-1))
      num_send(:) = 0
      ALLOCATE (req_array(1:numnodes, 4))

      ALLOCATE (iprocs(info%ngroup))
      CALL dbcsr_iterator_start(iter, matrix_in)
      DO WHILE (dbcsr_iterator_blocks_left(iter))
         CALL dbcsr_iterator_next_block(iter, blk_index(1), blk_index(2), nblock, transposed=tr, &
                                        row_size=blk_size(1), col_size=blk_size(2))
         CALL dbcsr_repl_get_stored_coordinates(matrix_out, blk_index(1), blk_index(2), iprocs)
         DO i = 1, SIZE(iprocs)
            num_send(2*iprocs(i)) = num_send(2*iprocs(i))+PRODUCT(blk_size)
            num_send(2*iprocs(i)+1) = num_send(2*iprocs(i)+1)+1
         ENDDO
      ENDDO
      CALL dbcsr_iterator_stop(iter)

      CALL timeset(routineN//"_alltoall", handle2)
      CALL mp_alltoall(num_send, num_rec, 2, mp_comm)
      CALL timestop(handle2)

      DO iproc = 0, numnodes-1
         num_entries_recv(iproc) = num_rec(2*iproc)
         num_blocks_recv(iproc) = num_rec(2*iproc+1)
         num_entries_send(iproc) = num_send(2*iproc)
         num_blocks_send(iproc) = num_send(2*iproc+1)

         CALL block_buffer_create(buffer_send(iproc), num_blocks_send(iproc), num_entries_send(iproc), &
                                  data_type)

         CALL block_buffer_create(buffer_recv(iproc), num_blocks_recv(iproc), num_entries_recv(iproc), &
                                  data_type)

      ENDDO

      CALL dbcsr_data_init(block)
      CALL dbcsr_data_new(block, dbcsr_type_1d_to_2d(data_type))
      CALL dbcsr_iterator_start(iter, matrix_in)
      DO WHILE (dbcsr_iterator_blocks_left(iter))
         CALL dbcsr_iterator_next_block(iter, blk_index(1), blk_index(2), block, tr, &
                                        row_size=blk_size(1), col_size=blk_size(2))
         CALL dbcsr_repl_get_stored_coordinates(matrix_out, blk_index(1), blk_index(2), iprocs)
         DO i = 1, SIZE(iprocs)
            CALL block_buffer_add_block(buffer_send(iprocs(i)), INT(blk_index, KIND=int_8), block)
         ENDDO
      ENDDO
      CALL dbcsr_iterator_stop(iter)
      CALL dbcsr_data_clear_pointer(block)

      CALL timeset(routineN//"_communicate_buffer", handle2)
      CALL dbcsr_r_communicate_buffer(mp_comm, buffer_recv, buffer_send, req_array)
      CALL timestop(handle2)

      nblk = SUM(num_blocks_recv)
      ALLOCATE (blks_to_allocate(nblk, 2))

      bcount = 0
      DO iproc = 0, numnodes - 1
         CALL block_buffer_get_index(buffer_recv(iproc), index_recv)
         blks_to_allocate(bcount + 1:bcount + SIZE(index_recv, 1), :) = INT(index_recv(:, :))
         bcount = bcount + SIZE(index_recv, 1)
         DEALLOCATE(index_recv)
      ENDDO

      CALL dbcsr_reserve_blocks(matrix_out%matrix, blks_to_allocate(:, 1), blks_to_allocate(:, 2))
      DEALLOCATE(blks_to_allocate)

      DO iproc = 0, numnodes-1
         ! First, we need to get the index to create block
         DO WHILE (block_buffer_blocks_left(buffer_recv(iproc)))
            CALL block_buffer_get_next_block(buffer_recv(iproc), ndata, blk_index_i8)
            CALL dbcsr_r_blk_sizes(matrix_out, blk_index_i8(1), blk_index_i8(2), blk_size(1), blk_size(2))
            CALL internal_data_allocate(block%d, blk_size)
            CALL block_buffer_get_next_block(buffer_recv(iproc), ndata, blk_index_i8, block)
            CALL dbcsr_put_block(matrix_out%matrix, INT(blk_index_i8(1)), INT(blk_index_i8(2)), block)
            CALL internal_data_deallocate(block%d)
         ENDDO

         CALL block_buffer_destroy(buffer_recv(iproc))
         CALL block_buffer_destroy(buffer_send(iproc))
      ENDDO
      CALL dbcsr_data_clear_pointer(block)
      CALL dbcsr_data_release(block)

      CALL dbcsr_r_finalize(matrix_out)

      CALL timestop(handle)

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param matrix_out ...
!> \param matrix_in ...
! **************************************************************************************************
   SUBROUTINE dbcsr_r_merge(matrix_out, matrix_in)
      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix_out
      TYPE(dbcsr_r_type), INTENT(INOUT)                  :: matrix_in

      CHARACTER(LEN=*), PARAMETER :: routineN = 'dbcsr_r_merge', routineP = moduleN//':'//routineN

      INTEGER                                            :: data_type, handle, handle2, iproc, mp_comm, &
                                                            mynode, nblock, ndata, numnodes, nblk, bcount
      INTEGER(KIND=int_8), ALLOCATABLE, DIMENSION(:, :)  :: index_recv
      INTEGER(KIND=int_8), DIMENSION(2)                  :: blk_index_i8
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: iprocs, num_blocks_recv, &
                                                            num_blocks_send, num_entries_recv, &
                                                            num_entries_send, num_rec, num_send
      INTEGER, ALLOCATABLE, DIMENSION(:, :)              :: req_array, blks_to_allocate, blks_to_allocate_u
      INTEGER, DIMENSION(2)                              :: blk_index, blk_size
      LOGICAL                                            :: tr
      TYPE(block_buffer_type), ALLOCATABLE, DIMENSION(:) :: buffer_recv, buffer_send
      TYPE(dbcsr_data_obj)                               :: block
      TYPE(dbcsr_iterator)                               :: iter
      TYPE(dbcsr_r_split_info)                           :: info

      CALL timeset(routineN, handle)

      CALL dbcsr_clear(matrix_out)

      data_type = dbcsr_r_get_data_type(matrix_in)

      info = dbcsr_r_info(matrix_in)
      CALL dbcsr_r_get_split_info(info, mp_comm=mp_comm)
      CALL mp_environ(numnodes, mynode, mp_comm)

      ALLOCATE (buffer_send(0:numnodes-1))
      ALLOCATE (buffer_recv(0:numnodes-1))
      ALLOCATE (num_blocks_recv(0:numnodes-1))
      ALLOCATE (num_blocks_send(0:numnodes-1))
      ALLOCATE (num_entries_recv(0:numnodes-1))
      ALLOCATE (num_entries_send(0:numnodes-1))
      ALLOCATE (num_rec(0:2*numnodes-1))
      ALLOCATE (num_send(0:2*numnodes-1))
      num_send(:) = 0
      ALLOCATE (req_array(1:numnodes, 4))

      ALLOCATE (iprocs(info%ngroup))

      CALL dbcsr_iterator_start(iter, matrix_in%matrix)
      DO WHILE (dbcsr_iterator_blocks_left(iter))
         CALL dbcsr_iterator_next_block(iter, blk_index(1), blk_index(2), nblock, transposed=tr, &
                                        row_size=blk_size(1), col_size=blk_size(2))
         CALL dbcsr_get_stored_coordinates(matrix_out, blk_index(1), blk_index(2), iproc)
         num_send(2*iproc) = num_send(2*iproc)+PRODUCT(blk_size)
         num_send(2*iproc+1) = num_send(2*iproc+1)+1
      ENDDO
      CALL dbcsr_iterator_stop(iter)

      CALL timeset(routineN//"_alltoall", handle2)
      CALL mp_alltoall(num_send, num_rec, 2, mp_comm)
      CALL timestop(handle2)

      DO iproc = 0, numnodes-1
         num_entries_recv(iproc) = num_rec(2*iproc)
         num_blocks_recv(iproc) = num_rec(2*iproc+1)
         num_entries_send(iproc) = num_send(2*iproc)
         num_blocks_send(iproc) = num_send(2*iproc+1)

         CALL block_buffer_create(buffer_send(iproc), num_blocks_send(iproc), num_entries_send(iproc), &
                                  data_type)

         CALL block_buffer_create(buffer_recv(iproc), num_blocks_recv(iproc), num_entries_recv(iproc), &
                                  data_type)

      ENDDO

      CALL dbcsr_data_init(block)
      CALL dbcsr_data_new(block, dbcsr_type_1d_to_2d(data_type))
      CALL dbcsr_iterator_start(iter, matrix_in%matrix)
      DO WHILE (dbcsr_iterator_blocks_left(iter))
         CALL dbcsr_iterator_next_block(iter, blk_index(1), blk_index(2), block, tr, &
                                        row_size=blk_size(1), col_size=blk_size(2))
         CALL dbcsr_get_stored_coordinates(matrix_out, blk_index(1), blk_index(2), iproc)
         CALL block_buffer_add_block(buffer_send(iproc), INT(blk_index, KIND=int_8), block)
      ENDDO

      CALL dbcsr_iterator_stop(iter)

      CALL timeset(routineN//"_communicate_buffer", handle2)
      CALL dbcsr_r_communicate_buffer(mp_comm, buffer_recv, buffer_send, req_array)
      CALL timestop(handle2)

      nblk = SUM(num_blocks_recv)
      ALLOCATE (blks_to_allocate(nblk, 2))

      bcount = 0
      DO iproc = 0, numnodes - 1
         CALL block_buffer_get_index(buffer_recv(iproc), index_recv)
         blks_to_allocate(bcount + 1:bcount + SIZE(index_recv, 1), :) = INT(index_recv(:, :))
         bcount = bcount + SIZE(index_recv, 1)
         DEALLOCATE(index_recv)
      ENDDO

      CALL index_unique(blks_to_allocate, blks_to_allocate_u)

      CALL dbcsr_reserve_blocks(matrix_out, blks_to_allocate_u(:, 1), blks_to_allocate_u(:, 2))
      DEALLOCATE(blks_to_allocate, blks_to_allocate_u)

      DO iproc = 0, numnodes-1
         ! First, we need to get the index to create block
         DO WHILE (block_buffer_blocks_left(buffer_recv(iproc)))
            CALL block_buffer_get_next_block(buffer_recv(iproc), ndata, blk_index_i8)
            blk_size(1) = dbcsr_blk_row_size(matrix_out, INT(blk_index_i8(1)))
            blk_size(2) = dbcsr_blk_column_size(matrix_out, INT(blk_index_i8(2)))
            CALL internal_data_allocate(block%d, blk_size)
            CALL block_buffer_get_next_block(buffer_recv(iproc), ndata, blk_index_i8, block)
            CALL dbcsr_put_block(matrix_out, INT(blk_index_i8(1)), INT(blk_index_i8(2)), block, summation=.TRUE.)
            CALL internal_data_deallocate(block%d)
         ENDDO
         CALL block_buffer_destroy(buffer_recv(iproc))
         CALL block_buffer_destroy(buffer_send(iproc))
      ENDDO
      CALL dbcsr_data_clear_pointer(block)
      CALL dbcsr_data_release(block)

      CALL dbcsr_finalize(matrix_out)

      CALL timestop(handle)
   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param index_in ...
!> \param index_out ...
! **************************************************************************************************
   SUBROUTINE index_unique(index_in, index_out)
      INTEGER, DIMENSION(:, :), INTENT(IN)   :: index_in
      INTEGER, ALLOCATABLE, &
         DIMENSION(:, :), INTENT(OUT)                    :: index_out

      INTEGER                                            :: blk, count, orig_size
      INTEGER, ALLOCATABLE, DIMENSION(:, :)  :: index_tmp
      INTEGER, DIMENSION(2)                  :: prev_index
      INTEGER, DIMENSION(1:SIZE(index_in, 1)&
         , 1:SIZE(index_in, 2))                          :: index_sorted

      orig_size = SIZE(index_in, 1)
      ALLOCATE (index_tmp(orig_size, 2))
      index_sorted(:, :) = index_in(:, :)
      CALL dbcsr_sort_indices(orig_size, index_sorted(:, 1), index_sorted(:, 2))
      count = 0
      prev_index(:) = [0, 0]
      DO blk = 1, orig_size
         IF (ANY(index_sorted(blk, :) .NE. prev_index(:))) THEN
            count = count+1
            index_tmp(count, :) = index_sorted(blk, :)
            prev_index(:) = index_sorted(blk, :)
         ENDIF
      ENDDO

      ALLOCATE (index_out(count, 2))
      index_out(:, :) = index_tmp(1:count, :)

   END SUBROUTINE

END MODULE
