!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2018  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Global data (distribution and block sizes)
!>
!> Due to the size of the longest matrix dimension (scaling quadratically with system sizes), these
!> are not stored as arrays but as function objects.
!>
!> \author Patrick Seewald
! **************************************************************************************************

MODULE dbcsr_rectangular_global

   USE dbcsr_kinds,                           ONLY: int_8
#include "../base/dbcsr_base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbcsr_rectangular_global'

   PUBLIC :: dbcsr_r_distribution, dbcsr_r_dist_cyclic, dbcsr_r_dist_arb, dbcsr_r_dist_repl, &
             dbcsr_r_rowcol_data, dbcsr_r_blk_size_arb, dbcsr_r_blk_size_repl

   TYPE, ABSTRACT :: dbcsr_r_distribution
      INTEGER :: nprowcol
      INTEGER(KIND=int_8) :: nmrowcol
   CONTAINS
      PROCEDURE(rowcol_dist), deferred :: dist
      PROCEDURE(dist_rowcols), deferred :: rowcols
   END TYPE

   TYPE, EXTENDS(dbcsr_r_distribution) :: dbcsr_r_dist_cyclic
      INTEGER :: split_size
   CONTAINS
      PROCEDURE :: dist => cyclic_dist
      PROCEDURE :: rowcols => cyclic_rowcols
   END TYPE

   TYPE, EXTENDS(dbcsr_r_distribution) :: dbcsr_r_dist_arb
      INTEGER, DIMENSION(:), ALLOCATABLE :: dist_vec
   CONTAINS
      PROCEDURE :: dist => arb_dist
      PROCEDURE :: rowcols => arb_rowcols
      FINAL :: dbcsr_r_release_dist_arb
   END TYPE

   TYPE, EXTENDS(dbcsr_r_distribution) :: dbcsr_r_dist_repl
      INTEGER, DIMENSION(:), ALLOCATABLE :: dist_vec
      INTEGER, DIMENSION(:), ALLOCATABLE :: dist_vec_last
      INTEGER :: nmrowcol_local
      INTEGER :: n_repl
      INTEGER :: dist_size
      INTEGER :: dist_size_last
   CONTAINS
      PROCEDURE :: dist => repl_dist
      PROCEDURE :: rowcols => repl_rowcols
      FINAL :: dbcsr_r_release_dist_repl
   END TYPE

   TYPE, ABSTRACT :: dbcsr_r_rowcol_data
      INTEGER(KIND=int_8) :: nmrowcol
      INTEGER(KIND=int_8) :: nfullrowcol
   CONTAINS
      PROCEDURE(rowcol_data), deferred :: DATA
   END TYPE

   TYPE, EXTENDS(dbcsr_r_rowcol_data) :: dbcsr_r_blk_size_arb
      INTEGER, DIMENSION(:), ALLOCATABLE :: blk_size_vec
   CONTAINS
      PROCEDURE :: DATA => blk_size_arb
      FINAL :: dbcsr_r_release_blk_size_arb
   END TYPE

   TYPE, EXTENDS(dbcsr_r_rowcol_data) :: dbcsr_r_blk_size_repl
      INTEGER, DIMENSION(:), ALLOCATABLE :: blk_size_vec
      INTEGER :: nmrowcol_local
   CONTAINS
      PROCEDURE :: DATA => blk_size_repl
      FINAL :: dbcsr_r_release_blk_size_repl
   END TYPE

   ABSTRACT INTERFACE
! **************************************************************************************************
!> \brief ...
!> \param t ...
!> \param rowcol ...
!> \return ...
! **************************************************************************************************
      FUNCTION rowcol_dist(t, rowcol)
         IMPORT :: dbcsr_r_distribution, int_8
         CLASS(dbcsr_r_distribution), INTENT(IN) :: t
         INTEGER(KIND=int_8), INTENT(IN) :: rowcol
         INTEGER :: rowcol_dist
      END FUNCTION

! **************************************************************************************************
!> \brief ...
!> \param t ...
!> \param dist ...
!> \return ...
! **************************************************************************************************
      FUNCTION dist_rowcols(t, dist)
         IMPORT :: dbcsr_r_distribution, int_8
         CLASS(dbcsr_r_distribution), INTENT(IN) :: t
         INTEGER, INTENT(IN) :: dist
         INTEGER(KIND=int_8), DIMENSION(:), ALLOCATABLE :: dist_rowcols
      END FUNCTION

! **************************************************************************************************
!> \brief ...
!> \param t ...
!> \param rowcol ...
!> \return ...
! **************************************************************************************************
      FUNCTION rowcol_data(t, rowcol)
         IMPORT :: dbcsr_r_rowcol_data, int_8
         CLASS(dbcsr_r_rowcol_data), INTENT(IN) :: t
         INTEGER(KIND=int_8), INTENT(IN) :: rowcol
         INTEGER :: rowcol_data
      END FUNCTION

   END INTERFACE

   INTERFACE dbcsr_r_dist_cyclic
      MODULE PROCEDURE new_dbcsr_r_dist_cyclic
   END INTERFACE

   INTERFACE dbcsr_r_dist_arb
      MODULE PROCEDURE new_dbcsr_r_dist_arb
   END INTERFACE

   INTERFACE dbcsr_r_dist_repl
      MODULE PROCEDURE new_dbcsr_r_dist_repl
   END INTERFACE

   INTERFACE dbcsr_r_blk_size_arb
      MODULE PROCEDURE new_dbcsr_r_blk_size_arb
   END INTERFACE

   INTERFACE dbcsr_r_blk_size_repl
      MODULE PROCEDURE new_dbcsr_r_blk_size_repl
   END INTERFACE

CONTAINS
! **************************************************************************************************
!> \brief ...
!> \param t ...
!> \param rowcol ...
!> \return ...
! **************************************************************************************************
   FUNCTION blk_size_arb(t, rowcol)
      CLASS(dbcsr_r_blk_size_arb), INTENT(IN) :: t
      INTEGER(KIND=int_8), INTENT(IN) :: rowcol
      INTEGER :: blk_size_arb
      blk_size_arb = t%blk_size_vec(rowcol)
   END FUNCTION

! **************************************************************************************************
!> \brief ...
!> \param t ...
!> \param rowcol ...
!> \return ...
! **************************************************************************************************
   FUNCTION blk_size_repl(t, rowcol)
      CLASS(dbcsr_r_blk_size_repl), INTENT(IN) :: t
      INTEGER(KIND=int_8), INTENT(IN) :: rowcol
      INTEGER :: blk_size_repl
      INTEGER :: igroup
      INTEGER :: rowcol_local

      igroup = INT((rowcol-1_int_8)/t%nmrowcol_local)
      rowcol_local = INT(MOD(rowcol-1_int_8, INT(t%nmrowcol_local, KIND=int_8)))+1
      blk_size_repl = t%blk_size_vec(rowcol_local)

   END FUNCTION

! **************************************************************************************************
!> \brief ...
!> \param blk_size_vec ...
!> \return ...
! **************************************************************************************************
   FUNCTION new_dbcsr_r_blk_size_arb(blk_size_vec)
      INTEGER, DIMENSION(:), INTENT(IN)                  :: blk_size_vec
      TYPE(dbcsr_r_blk_size_arb)                         :: new_dbcsr_r_blk_size_arb

      ALLOCATE (new_dbcsr_r_blk_size_arb%blk_size_vec(SIZE(blk_size_vec)))
      new_dbcsr_r_blk_size_arb%blk_size_vec(:) = blk_size_vec(:)
      new_dbcsr_r_blk_size_arb%nmrowcol = SIZE(blk_size_vec)
      new_dbcsr_r_blk_size_arb%nfullrowcol = SUM(blk_size_vec)
   END FUNCTION

! **************************************************************************************************
!> \brief ...
!> \param blk_size_vec ...
!> \param n_repl ...
!> \return ...
! **************************************************************************************************
   FUNCTION new_dbcsr_r_blk_size_repl(blk_size_vec, n_repl)
      INTEGER, DIMENSION(:), INTENT(IN)                  :: blk_size_vec
      INTEGER, INTENT(IN)                                :: n_repl
      TYPE(dbcsr_r_blk_size_repl)                        :: new_dbcsr_r_blk_size_repl

      new_dbcsr_r_blk_size_repl%nmrowcol_local = SIZE(blk_size_vec)
      ALLOCATE (new_dbcsr_r_blk_size_repl%blk_size_vec(new_dbcsr_r_blk_size_repl%nmrowcol_local))
      new_dbcsr_r_blk_size_repl%blk_size_vec(:) = blk_size_vec(:)
      new_dbcsr_r_blk_size_repl%nmrowcol = new_dbcsr_r_blk_size_repl%nmrowcol_local*n_repl
      new_dbcsr_r_blk_size_repl%nfullrowcol = SUM(blk_size_vec)*n_repl
   END FUNCTION

! **************************************************************************************************
!> \brief ...
!> \param t ...
!> \param rowcol ...
!> \return ...
! **************************************************************************************************
   FUNCTION arb_dist(t, rowcol)
      CLASS(dbcsr_r_dist_arb), INTENT(IN) :: t
      INTEGER(KIND=int_8), INTENT(IN) :: rowcol
      INTEGER :: arb_dist

      arb_dist = t%dist_vec(rowcol)
   END FUNCTION

! **************************************************************************************************
!> \brief ...
!> \param t ...
!> \param rowcol ...
!> \return ...
! **************************************************************************************************
   FUNCTION repl_dist(t, rowcol)
      CLASS(dbcsr_r_dist_repl), INTENT(IN) :: t
      INTEGER(KIND=int_8), INTENT(IN) :: rowcol
      INTEGER :: repl_dist
      INTEGER :: rowcol_local
      INTEGER :: igroup

      igroup = INT((rowcol-1_int_8)/t%nmrowcol_local)
      rowcol_local = INT(MOD(rowcol-1_int_8, INT(t%nmrowcol_local, KIND=int_8)))+1

      IF (igroup < t%n_repl-1) THEN
         repl_dist = t%dist_vec(rowcol_local)+igroup*t%dist_size
      ELSEIF (igroup == t%n_repl-1) THEN
         repl_dist = t%dist_vec_last(rowcol_local)+igroup*t%dist_size
      ENDIF

   END FUNCTION

! **************************************************************************************************
!> \brief ...
!> \param t ...
!> \param dist ...
!> \return ...
! **************************************************************************************************
   FUNCTION repl_rowcols(t, dist)
      CLASS(dbcsr_r_dist_repl), INTENT(IN) :: t
      INTEGER, INTENT(IN) :: dist
      INTEGER :: nrowcols
      INTEGER(KIND=int_8), DIMENSION(:), ALLOCATABLE :: repl_rowcols, rowcols_tmp
      INTEGER :: igroup
      INTEGER :: rowcol, count
      LOGICAL :: cond

      igroup = dist/t%dist_size

      nrowcols = t%nmrowcol_local
      count = 0
      ALLOCATE (rowcols_tmp(nrowcols))
      rowcols_tmp(:) = 0
      DO rowcol = 1, nrowcols
         IF (igroup < t%n_repl-1) THEN
            cond = t%dist_vec(rowcol)+igroup*t%dist_size == dist
         ELSE
            cond = t%dist_vec_last(rowcol)+igroup*t%dist_size == dist
         ENDIF

         IF (cond) THEN
            count = count+1
            rowcols_tmp(count) = rowcol
         ENDIF
      ENDDO

      ALLOCATE (repl_rowcols(count))
      repl_rowcols(:) = rowcols_tmp(1:count)+igroup*t%nmrowcol_local

   END FUNCTION

! **************************************************************************************************
!> \brief ...
!> \param t ...
!> \param dist ...
!> \return ...
! **************************************************************************************************
   FUNCTION arb_rowcols(t, dist)
      CLASS(dbcsr_r_dist_arb), INTENT(IN) :: t
      INTEGER, INTENT(IN) :: dist
      INTEGER(KIND=int_8) :: rowcol, nrowcols
      INTEGER(KIND=int_8), DIMENSION(:), ALLOCATABLE :: arb_rowcols, rowcols_tmp
      INTEGER :: count

      nrowcols = t%nmrowcol
      count = 0
      ALLOCATE (rowcols_tmp(nrowcols))
      rowcols_tmp(:) = 0
      DO rowcol = 1, nrowcols
         IF (t%dist_vec(rowcol) == dist) THEN
            count = count+1
            rowcols_tmp(count) = rowcol
         ENDIF
      ENDDO

      ALLOCATE (arb_rowcols(count))
      arb_rowcols(:) = rowcols_tmp(1:count)
   END FUNCTION

! **************************************************************************************************
!> \brief ...
!> \param split_size ...
!> \param nprowcol ...
!> \param nmrowcol ...
!> \return ...
! **************************************************************************************************
   FUNCTION new_dbcsr_r_dist_cyclic(split_size, nprowcol, nmrowcol)
      INTEGER, INTENT(IN)                                :: split_size, nprowcol
      INTEGER(KIND=int_8), INTENT(IN)                    :: nmrowcol
      TYPE(dbcsr_r_dist_cyclic)                          :: new_dbcsr_r_dist_cyclic

      new_dbcsr_r_dist_cyclic%split_size = split_size
      new_dbcsr_r_dist_cyclic%nprowcol = nprowcol
      new_dbcsr_r_dist_cyclic%nmrowcol = nmrowcol
   END FUNCTION

! **************************************************************************************************
!> \brief ...
!> \param dist_vec ...
!> \param nprowcol ...
!> \param nmrowcol ...
!> \return ...
! **************************************************************************************************
   FUNCTION new_dbcsr_r_dist_arb(dist_vec, nprowcol, nmrowcol)
      INTEGER, DIMENSION(:), INTENT(IN)                  :: dist_vec
      INTEGER, INTENT(IN)                                :: nprowcol
      INTEGER(KIND=int_8), INTENT(IN)                    :: nmrowcol
      TYPE(dbcsr_r_dist_arb)                             :: new_dbcsr_r_dist_arb

      ALLOCATE (new_dbcsr_r_dist_arb%dist_vec(nmrowcol))
      new_dbcsr_r_dist_arb%dist_vec(:) = dist_vec(:)
      new_dbcsr_r_dist_arb%nprowcol = nprowcol
      new_dbcsr_r_dist_arb%nmrowcol = nmrowcol
   END FUNCTION

! **************************************************************************************************
!> \brief ...
!> \param dist_vec ...
!> \param nprowcol ...
!> \param nmrowcol ...
!> \param n_repl ...
!> \param dist_size ...
!> \param dist_size_last ...
!> \return ...
! **************************************************************************************************
   FUNCTION new_dbcsr_r_dist_repl(dist_vec, nprowcol, nmrowcol, n_repl, dist_size, dist_size_last)
      INTEGER, DIMENSION(:), INTENT(IN)                  :: dist_vec
      INTEGER, INTENT(IN)                                :: nprowcol, nmrowcol, n_repl, dist_size, &
                                                            dist_size_last
      TYPE(dbcsr_r_dist_repl)                            :: new_dbcsr_r_dist_repl

      new_dbcsr_r_dist_repl%n_repl = n_repl
      new_dbcsr_r_dist_repl%dist_size = dist_size
      ALLOCATE (new_dbcsr_r_dist_repl%dist_vec(nmrowcol))
      ALLOCATE (new_dbcsr_r_dist_repl%dist_vec_last(nmrowcol))
      new_dbcsr_r_dist_repl%dist_vec(:) = MOD(dist_vec(:), dist_size)
      new_dbcsr_r_dist_repl%dist_vec_last(:) = MOD(dist_vec(:), dist_size_last)
      new_dbcsr_r_dist_repl%nprowcol = nprowcol
      new_dbcsr_r_dist_repl%nmrowcol_local = nmrowcol
      new_dbcsr_r_dist_repl%nmrowcol = nmrowcol*n_repl
   END FUNCTION

! **************************************************************************************************
!> \brief ...
!> \param t ...
!> \param rowcol ...
!> \return ...
! **************************************************************************************************
   FUNCTION cyclic_dist(t, rowcol)
      CLASS(dbcsr_r_dist_cyclic), INTENT(IN) :: t
      INTEGER(KIND=int_8), INTENT(IN) :: rowcol
      INTEGER :: cyclic_dist

      cyclic_dist = INT(MOD((rowcol-1)/INT(t%split_size, KIND=int_8), INT(t%nprowcol, KIND=int_8)))

   END FUNCTION

! **************************************************************************************************
!> \brief ...
!> \param t ...
!> \param dist ...
!> \return ...
! **************************************************************************************************
   FUNCTION cyclic_rowcols(t, dist)
      CLASS(dbcsr_r_dist_cyclic), INTENT(IN) :: t
      INTEGER, INTENT(IN) :: dist
      INTEGER(KIND=int_8), DIMENSION(:), ALLOCATABLE :: cyclic_rowcols
      INTEGER :: count, nsplit, isplit, irowcol, max_size
      INTEGER(KIND=int_8) :: rowcol
      INTEGER(KIND=int_8), DIMENSION(:), ALLOCATABLE :: rowcols_tmp

      nsplit = INT((t%nmrowcol-1)/INT(t%split_size, KIND=int_8)+1_int_8)
      max_size = nsplit*t%split_size
      ALLOCATE (rowcols_tmp(max_size))
      rowcols_tmp(:) = 0
      count = 0
      loop: DO isplit = 1, nsplit
         DO irowcol = 1, t%split_size
            rowcol = INT((dist+(isplit-1)*t%nprowcol), KIND=int_8)*INT(t%split_size, KIND=int_8)+ &
                     INT(irowcol, KIND=int_8)
            IF (rowcol > t%nmrowcol) THEN
               EXIT loop
            ELSE
               count = count+1
               rowcols_tmp(count) = rowcol
            ENDIF
         ENDDO
      ENDDO loop

      ALLOCATE (cyclic_rowcols(count))
      cyclic_rowcols(:) = rowcols_tmp(1:count)
   END FUNCTION

! **************************************************************************************************
!> \brief ...
!> \param t ...
! **************************************************************************************************
   SUBROUTINE dbcsr_r_release_blk_size_arb(t)
      TYPE(dbcsr_r_blk_size_arb)                         :: t

      DEALLOCATE (t%blk_size_vec)
   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param t ...
! **************************************************************************************************
   SUBROUTINE dbcsr_r_release_blk_size_repl(t)
      TYPE(dbcsr_r_blk_size_repl)                        :: t

      DEALLOCATE (t%blk_size_vec)
   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param t ...
! **************************************************************************************************
   SUBROUTINE dbcsr_r_release_dist_arb(t)
      TYPE(dbcsr_r_dist_arb)                             :: t

      DEALLOCATE (t%dist_vec)
   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param t ...
! **************************************************************************************************
   SUBROUTINE dbcsr_r_release_dist_repl(t)
      TYPE(dbcsr_r_dist_repl)                            :: t

      DEALLOCATE (t%dist_vec)
      DEALLOCATE (t%dist_vec_last)
   END SUBROUTINE

END MODULE
