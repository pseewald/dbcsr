!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2018  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief routines to reshape/replicate/merge rectangular matrices,
!>        low-level code using Fypp preprocessor.
!> \author Patrick Seewald
! **************************************************************************************************

MODULE dbcsr_rectangular_reshape_low

#include "../base/dbcsr_base_uses.f90"
#:include "dbcsr_rectangular.fypp"

   USE dbcsr_kinds, ONLY: int_8, real_4, real_8
   USE dbcsr_data_types, ONLY: ${uselist(dtype_float_param)}$
   USE dbcsr_data_methods_low, ONLY: ${uselist(["dbcsr_get_data_p_2d_" + s for s in dtype_float_suffix_dbcsr])}$
   USE dbcsr_data_types, ONLY: dbcsr_data_obj
   USE dbcsr_mpiwrap, ONLY: mp_irecv, mp_isend, mp_environ, mp_waitall
   USE dbcsr_rectangular_util, ONLY: swap
   USE dbcsr_data_methods, ONLY: dbcsr_data_get_sizes, dbcsr_data_set_pointer
   USE dbcsr_data_methods_low, ONLY: internal_data_deallocate

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbcsr_rectangular_reshape_low'

   PUBLIC :: block_buffer_type, block_buffer_add_block, block_buffer_create, block_buffer_destroy, &
             block_buffer_get_next_block, dbcsr_r_communicate_buffer

   TYPE block_buffer_type
      INTEGER :: nblock=-1
      INTEGER(KIND=int_8), DIMENSION(:, :), ALLOCATABLE :: indx
#:for dparam, dtype, dsuffix in dtype_float_list
      ${dtype}$, DIMENSION(:), ALLOCATABLE       :: msg_${dsuffix}$
#:endfor
      INTEGER :: data_type = -1
      INTEGER :: endpos = -1
   END TYPE

   INTERFACE block_buffer_get_next_block
      MODULE PROCEDURE block_buffer_get_next_area_block
#!for dparam, dtype, dsuffix in dtype_float_list
      !MODULE PROCEDURE block_buffer_get_next_block_${dsuffix}$ ! issue:ambigous interface
#!endfor

   END INTERFACE

   INTERFACE block_buffer_add_block
      MODULE PROCEDURE block_buffer_add_area_block
#:for dparam, dtype, dsuffix in dtype_float_list
      MODULE PROCEDURE block_buffer_add_block_${dsuffix}$
#:endfor
   END INTERFACE

CONTAINS

! **************************************************************************************************
!> \brief Create block buffer for MPI communication.
!> \param buffer block buffer
!> \param nblock number of blocks
!> \param ndata total number of block entries
!> \param data_type ...
! **************************************************************************************************
   SUBROUTINE block_buffer_create(buffer, nblock, ndata, data_type)
      TYPE(block_buffer_type), INTENT(OUT) :: buffer
      INTEGER, INTENT(IN) :: nblock, ndata, data_type

      buffer%nblock = nblock
      buffer%data_type = data_type
      buffer%endpos = 0
      SELECT CASE (data_type)
#:for dparam, dtype, dsuffix in dtype_float_list
      CASE (${dparam}$)
         ALLOCATE (buffer%msg_${dsuffix}$(ndata))
#:endfor
      END SELECT
      ALLOCATE (buffer%indx(nblock, 3))
   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param buffer ...
! **************************************************************************************************
   SUBROUTINE block_buffer_destroy(buffer)
      TYPE(block_buffer_type), INTENT(INOUT) :: buffer

      SELECT CASE (buffer%data_type)
#:for dparam, dtype, dsuffix in dtype_float_list
      CASE (${dparam}$)
         DEALLOCATE (buffer%msg_${dsuffix}$)
#:endfor
      END SELECT
      DEALLOCATE (buffer%indx)
      buffer%nblock = -1
      buffer%data_type = -1
      buffer%endpos = -1
   END SUBROUTINE block_buffer_destroy

! **************************************************************************************************
!> \brief ...
!> \param buffer ...
!> \param index ...
!> \param block ...
!> \param transposed ...
! **************************************************************************************************
   SUBROUTINE block_buffer_add_area_block(buffer, index, block, transposed)
      TYPE(block_buffer_type), INTENT(INOUT)      :: buffer
      INTEGER(KIND=int_8), DIMENSION(2), &
         INTENT(IN)                               :: index
      TYPE(dbcsr_data_obj), INTENT(IN)         :: block
      LOGICAL, INTENT(IN), OPTIONAL :: transposed

#:for dparam, dtype, dsuffix in dtype_float_list
      ${dtype}$, DIMENSION(:,:), POINTER :: block_${dsuffix}$
#:endfor

      SELECT CASE(buffer%data_type)
#:for dparam, dtype, dsuffix, dsuffix_dbcsr in dtype_float_list_dbcsr
      CASE(${dparam}$)
         block_${dsuffix}$ => dbcsr_get_data_p_2d_${dsuffix_dbcsr}$(block)
         CALL block_buffer_add_block_${dsuffix}$(buffer, index, block_${dsuffix}$, transposed)
#:endfor
      END SELECT
   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param buffer ...
!> \param ndata ...
!> \param index ...
!> \param block ...
!> \param advance_iter ...
! **************************************************************************************************
   SUBROUTINE block_buffer_get_next_area_block(buffer, ndata, index, block, advance_iter)
      TYPE(block_buffer_type), INTENT(INOUT)      :: buffer
      INTEGER, INTENT(OUT)                        :: ndata
      INTEGER(KIND=int_8), DIMENSION(2), &
         INTENT(OUT)                              :: index
      TYPE(dbcsr_data_obj), INTENT(INOUT), OPTIONAL           :: block
      LOGICAL, INTENT(IN), OPTIONAL               :: advance_iter
      LOGICAL :: valid
      INTEGER, DIMENSION(2) :: sizes
#:for dparam, dtype, dsuffix in dtype_float_list
      ${dtype}$, DIMENSION(:,:), POINTER :: data_${dsuffix}$
#:endfor
      INTEGER :: handle

      CHARACTER(LEN=*), PARAMETER :: routineN = 'block_buffer_get_next_area_block', &
                                     routineP = moduleN//':'//routineN

      CALL timeset(routineN, handle)
      IF (PRESENT(block)) THEN
         CALL dbcsr_data_get_sizes(block, sizes, valid)
         DBCSR_ASSERT(valid)
         CALL internal_data_deallocate(block%d)
      ENDIF

      SELECT CASE(buffer%data_type)
#:for dparam, dtype, dsuffix in dtype_float_list
      CASE(${dparam}$)
         IF(PRESENT(block)) THEN
            ALLOCATE(data_${dsuffix}$(sizes(1), sizes(2)))
            CALL block_buffer_get_next_block_${dsuffix}$(buffer, ndata, index, data_${dsuffix}$, advance_iter=advance_iter)
            CALL dbcsr_data_set_pointer(block, data_${dsuffix}$)
         ELSE
            CALL block_buffer_get_next_block_${dsuffix}$(buffer, ndata, index, advance_iter=advance_iter)
         ENDIF
#:endfor
      END SELECT
      CALL timestop(handle)
   END SUBROUTINE

#:for dparam, dtype, dsuffix in dtype_float_list
! **************************************************************************************************
!> \brief insert a block into block buffer (at current iterator position)
!> \param buffer ...
!> \param index index of block
!> \param block block
!> \param transposed ...
! **************************************************************************************************
   SUBROUTINE block_buffer_add_block_${dsuffix}$(buffer, index, block, transposed)
      TYPE(block_buffer_type), INTENT(INOUT)      :: buffer
      INTEGER(KIND=int_8), DIMENSION(2), &
         INTENT(IN)                               :: index
      ${dtype}$, DIMENSION(:,:), INTENT(IN)                  :: block
      LOGICAL, INTENT(IN), OPTIONAL :: transposed
      INTEGER(KIND=int_8), DIMENSION(2)   :: index_prv

      LOGICAL :: tr
      INTEGER :: p
      INTEGER :: ndata
      INTEGER :: p_data

      IF(PRESENT(transposed)) THEN
         tr = transposed
      ELSE
         tr = .FALSE.
      ENDIF

      index_prv(:) = index(:)
      IF(tr) THEN
         CALL swap(index_prv)
      ENDIF
      ndata = PRODUCT(SHAPE(block))

      DBCSR_ASSERT(buffer%data_type .EQ. ${dparam}$)
      p = buffer%endpos
      IF (p .EQ. 0) THEN
         p_data = 0
      ELSE
         p_data = INT(buffer%indx(p, 3))
      ENDIF

      IF(tr) THEN
         buffer%msg_${dsuffix}$(p_data+1:p_data+ndata) = RESHAPE(TRANSPOSE(block),[ndata])
      ELSE
         buffer%msg_${dsuffix}$(p_data+1:p_data+ndata) = RESHAPE(block,[ndata])
      ENDIF

      buffer%indx(p+1, 1:2) = index_prv(:)
      IF (p > 0) THEN
         buffer%indx(p+1,3) = buffer%indx(p,3)+INT(ndata, KIND=int_8)
      ELSE
         buffer%indx(p+1,3) = INT(ndata, KIND=int_8)
      ENDIF
      buffer%endpos = buffer%endpos+1
   END SUBROUTINE
#:endfor

#:for dparam, dtype, dsuffix in dtype_float_list
! **************************************************************************************************
!> \brief get next block from buffer. Iterator is advanced only if block is retrieved or advance_iter.
!> \param buffer ...
!> \param ndata ...
!> \param index ...
!> \param block ...
!> \param advance_iter ...
! **************************************************************************************************
   SUBROUTINE block_buffer_get_next_block_${dsuffix}$(buffer, ndata, index, block, advance_iter)
      TYPE(block_buffer_type), INTENT(INOUT)      :: buffer
      INTEGER, INTENT(OUT)                        :: ndata
      INTEGER(KIND=int_8), DIMENSION(2), &
         INTENT(OUT)                              :: index
      ${dtype}$, DIMENSION(:,:), INTENT(OUT), OPTIONAL     :: block
      LOGICAL, INTENT(IN), OPTIONAL               :: advance_iter
      INTEGER :: p, p_data
      LOGICAL :: do_advance

      do_advance = .FALSE.
      IF (PRESENT(advance_iter)) THEN
         do_advance = advance_iter
      ELSE IF (PRESENT(block)) THEN
         do_advance = .TRUE.
      ENDIF

      DBCSR_ASSERT(buffer%data_type .EQ. ${dparam}$)

      p = buffer%endpos
      IF (p .EQ. 0) THEN
         p_data = 0
      ELSE
         p_data = INT(buffer%indx(p, 3))
      ENDIF

      IF (p > 0) THEN
         ndata = INT(buffer%indx(p+1, 3)-buffer%indx(p, 3))
      ELSE
         ndata = INT(buffer%indx(p+1, 3))
      ENDIF
      index(:) = buffer%indx(p+1,1:2)

      IF (PRESENT(block)) THEN
         block(:,:) = RESHAPE(buffer%msg_${dsuffix}$(p_data+1:p_data+ndata),SHAPE(block))
      ENDIF

      IF(do_advance) buffer%endpos = buffer%endpos+1
   END SUBROUTINE
#:endfor

! **************************************************************************************************
!> \brief communicate buffer
!> \param mp_comm ...
!> \param buffer_recv ...
!> \param buffer_send ...
!> \param req_array ...
! **************************************************************************************************
   SUBROUTINE dbcsr_r_communicate_buffer(mp_comm, buffer_recv, buffer_send, req_array)
      INTEGER, INTENT(IN)                    :: mp_comm
      TYPE(block_buffer_type), DIMENSION(0:), INTENT(INOUT) :: buffer_recv, buffer_send
      INTEGER, DIMENSION(:, :), INTENT(OUT)               :: req_array

      INTEGER                                :: iproc, mynode, numnodes, rec_counter, &
                                                send_counter
      INTEGER                                   :: handle
      CHARACTER(LEN=*), PARAMETER :: routineN = 'dbcsr_r_communicate_buffer', &
                                     routineP = moduleN//':'//routineN

      CALL timeset(routineN, handle)
      CALL mp_environ(numnodes, mynode, mp_comm)

      IF (numnodes > 1) THEN

         send_counter = 0
         rec_counter = 0

         DO iproc = 0, numnodes-1
            IF (buffer_recv(iproc)%nblock > 0) THEN
               rec_counter = rec_counter+1
               CALL mp_irecv(buffer_recv(iproc)%indx, iproc, mp_comm, req_array(rec_counter, 3), tag=4)
               SELECT CASE (buffer_recv (iproc)%data_type)
#:for dparam, dtype, dsuffix in dtype_float_list
               CASE (${dparam}$)
                  CALL mp_irecv(buffer_recv(iproc)%msg_${dsuffix}$, iproc, mp_comm, req_array(rec_counter, 4), tag=7)
#:endfor
               END SELECT
            END IF
         END DO

         DO iproc = 0, numnodes-1
            IF (buffer_send(iproc)%nblock > 0) THEN
               send_counter = send_counter+1
               CALL mp_isend(buffer_send(iproc)%indx, iproc, mp_comm, req_array(send_counter, 1), tag=4)
               SELECT CASE (buffer_recv (iproc)%data_type)
#:for dparam, dtype, dsuffix in dtype_float_list
               CASE (${dparam}$)
                  CALL mp_isend(buffer_send(iproc)%msg_${dsuffix}$, iproc, mp_comm, req_array(send_counter, 2), tag=7)
#:endfor
               END SELECT
            END IF
         END DO

         IF (send_counter > 0) THEN
            CALL mp_waitall(req_array(1:send_counter, 1:2))
         ENDIF
         IF (rec_counter > 0) THEN
            CALL mp_waitall(req_array(1:rec_counter, 3:4))
         ENDIF

      ELSE
         IF (buffer_recv(0)%nblock > 0) THEN
            buffer_recv(0)%indx(:, :) = buffer_send(0)%indx(:, :)
            SELECT CASE (buffer_recv (0)%data_type)
#:for dparam, dtype, dsuffix in dtype_float_list
            CASE (${dparam}$)
               buffer_recv(0)%msg_${dsuffix}$(:) = buffer_send(0)%msg_${dsuffix}$(:)
#:endfor
            END SELECT
         ENDIF
      ENDIF
      CALL timestop(handle)
   END SUBROUTINE

END MODULE
