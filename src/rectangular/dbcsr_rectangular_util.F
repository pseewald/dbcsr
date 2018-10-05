!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2018  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief often used hacks / utilities
!> \author Patrick Seewald
! **************************************************************************************************

MODULE dbcsr_rectangular_util

   USE dbcsr_mp_methods,                ONLY: dbcsr_mp_new
   USE dbcsr_types,                     ONLY: dbcsr_mp_obj
   USE dbcsr_kinds,                           ONLY: int_8
   USE dbcsr_mpiwrap,                 ONLY: mp_cart_rank,&
                                              mp_environ
#include "../base/dbcsr_base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbcsr_rectangular_util'

   PUBLIC :: dbcsr_mp_environ, swap

   INTERFACE swap
      MODULE PROCEDURE swap_i8
      MODULE PROCEDURE swap_i
   END INTERFACE

CONTAINS
! **************************************************************************************************
!> \brief ...
!> \param mp_comm ...
!> \return ...
! **************************************************************************************************
   FUNCTION dbcsr_mp_environ(mp_comm) ! hack, to be avoided
      INTEGER, INTENT(IN)                                :: mp_comm
      TYPE(dbcsr_mp_obj)                                 :: dbcsr_mp_environ

      INTEGER                                            :: mynode, numnodes, pcol, prow
      INTEGER, ALLOCATABLE, DIMENSION(:, :)              :: pgrid
      INTEGER, DIMENSION(2)                              :: coord, mycoord, npdims

      CALL mp_environ(numnodes, npdims, mycoord, mp_comm)
      CALL mp_environ(numnodes, mynode, mp_comm)
      ALLOCATE (pgrid(0:npdims(1)-1, 0:npdims(2)-1))
      DO prow = 0, npdims(1)-1
         DO pcol = 0, npdims(2)-1
            coord = (/prow, pcol/)
            CALL mp_cart_rank(mp_comm, coord, pgrid(prow, pcol))
         ENDDO
      ENDDO

      DBCSR_ASSERT(mynode == pgrid(mycoord(1), mycoord(2)))

      CALL dbcsr_mp_new(dbcsr_mp_environ, pgrid, mp_comm, mynode, numnodes, mycoord(1), mycoord(2))
   END FUNCTION

! **************************************************************************************************
!> \brief ...
!> \param arr ...
! **************************************************************************************************
   SUBROUTINE swap_i8(arr)
      INTEGER(KIND=int_8), DIMENSION(2), INTENT(INOUT)   :: arr

      INTEGER(KIND=int_8)                                :: tmp

      tmp = arr(1)
      arr(1) = arr(2)
      arr(2) = tmp
   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param arr ...
! **************************************************************************************************
   SUBROUTINE swap_i(arr)
      INTEGER, DIMENSION(2), INTENT(INOUT)               :: arr

      INTEGER                                            :: tmp

      tmp = arr(1)
      arr(1) = arr(2)
      arr(2) = tmp
   END SUBROUTINE

END MODULE
