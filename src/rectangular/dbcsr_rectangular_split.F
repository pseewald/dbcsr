!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2018  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief methods to split rectangular matrices along longest dimension.
!>
!> Basically, we are splitting MPI grid and each subgrid contains its own DBCSR matrix.
!> todo: may need several split heuristics:
!> - recursive carma: any split (different group sizes)
!> - non-recursive carma: decompose into groups of same size
!> \author Patrick Seewald
! **************************************************************************************************
MODULE dbcsr_rectangular_split

   USE dbcsr_dist_methods, ONLY: dbcsr_distribution_hold
   USE dbcsr_rectangular_base_types, ONLY: dbcsr_r_distribution_type, &
      dbcsr_r_split_info
   USE dbcsr_rectangular_global, ONLY: dbcsr_r_distribution
   USE dbcsr_rectangular_util, ONLY: swap
   USE dbcsr_toollib, ONLY: sort
   USE dbcsr_kinds, ONLY: int_8
   USE dbcsr_machine, ONLY: default_output_unit
   USE dbcsr_mpiwrap, ONLY: mp_bcast, &
      mp_cart_create, &
      mp_comm_dup, &
      mp_comm_free, &
      mp_comm_split_direct, &
      mp_dims_create, &
      mp_environ
#include "../base/dbcsr_base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: dbcsr_r_split, &
             dbcsr_r_split_rows_or_cols, &
             dbcsr_r_release_info, dbcsr_r_print_info, rowsplit, colsplit, dbcsr_r_get_split_info, &
             world_to_group_proc_map, group_to_world_proc_map, group_to_mrowcol, mrowcol_to_group, &
             block_index_global_to_local, block_index_local_to_global, &
             dbcsr_r_mp_comm, dbcsr_r_info_hold, dbcsr_r_copY_distribution

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbcsr_rectangular_split'

   INTEGER, PARAMETER :: rowsplit = 1, colsplit = 2

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param mp_comm ...
!> \param nblkrows_total ...
!> \param nblkcols_total ...
!> \return ...
! **************************************************************************************************
   FUNCTION dbcsr_r_mp_comm(mp_comm, nblkrows_total, nblkcols_total)
      INTEGER, INTENT(IN)                                :: mp_comm
      INTEGER(KIND=int_8), INTENT(IN)                    :: nblkrows_total, nblkcols_total
      INTEGER                                            :: dbcsr_r_mp_comm

      CHARACTER(LEN=*), PARAMETER :: routineN = 'dbcsr_r_mp_comm', &
                                     routineP = moduleN//':'//routineN

      INTEGER                                            :: group_size, handle, iproc, numproc, s, &
                                                            split_rowcol
      INTEGER(KIND=int_8)                                :: d_max, d_min
      INTEGER, DIMENSION(2)                              :: group_dims, mdims_group, myploc, npdims, &
                                                            sort_indices

      CALL timeset(routineN, handle)

      d_max = MAX(nblkrows_total, nblkcols_total)
      d_min = MIN(nblkrows_total, nblkcols_total)

      s = INT((d_max - 1)/d_min + 1)

      CALL mp_environ(numproc, iproc, mp_comm)

      DO WHILE (MOD(numproc, s) .NE. 0)
         s = s - 1
      ENDDO

      group_size = numproc/s
      group_dims(:) = 0

      CALL mp_dims_create(group_size, group_dims)
      CALL sort(group_dims, 2, sort_indices)

      IF (d_max == nblkrows_total) THEN
         split_rowcol = rowsplit
      ELSE
         split_rowcol = colsplit
      ENDIF

      SELECT CASE (split_rowcol)
      CASE (rowsplit)
         mdims_group = [INT(nblkrows_total/s), INT(nblkcols_total)]
      CASE (colsplit)
         mdims_group = [INT(nblkrows_total), INT(nblkcols_total/s)]
      END SELECT

      IF (mdims_group(1) .GT. mdims_group(2)) THEN
         CALL swap(group_dims)
      ENDIF

      SELECT CASE (split_rowcol)
      CASE (rowsplit)
         npdims(:) = [group_dims(1)*s, group_dims(2)]
      CASE (colsplit)
         npdims(:) = [group_dims(1), group_dims(2)*s]
      END SELECT

      CALL mp_cart_create(mp_comm, 2, npdims, myploc, dbcsr_r_mp_comm)

      CALL timestop(handle)
   END FUNCTION

! **************************************************************************************************
!> \brief ...
!> \param split_info ...
!> \param mp_comm ...
!> \param nblkrows_total ...
!> \param nblkcols_total ...
!> \param nosplit ...
! **************************************************************************************************
   SUBROUTINE dbcsr_r_split(split_info, mp_comm, nblkrows_total, nblkcols_total, nosplit)
      TYPE(dbcsr_r_split_info), INTENT(OUT)              :: split_info
      INTEGER, INTENT(IN)                                :: mp_comm
      INTEGER(KIND=int_8), INTENT(IN)                    :: nblkrows_total, nblkcols_total
      LOGICAL, INTENT(IN), OPTIONAL                      :: nosplit

      CHARACTER(LEN=*), PARAMETER :: routineN = 'dbcsr_r_split', routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, igroup, iproc, numproc, s, &
                                                            split_rowcol
      INTEGER(KIND=int_8)                                :: d_max, d_min
      INTEGER, DIMENSION(2)                              :: pcoord, pdims, pdims_group
      LOGICAL                                            :: nosplit_prv

      CALL timeset(routineN, handle)

      IF (PRESENT(nosplit)) THEN
         nosplit_prv = nosplit
      ELSE
         nosplit_prv = .FALSE.
      ENDIF

      IF (nosplit_prv) THEN
         s = 1
         igroup = 0
         split_rowcol = rowsplit

      ELSE

         d_max = MAX(nblkrows_total, nblkcols_total)
         d_min = MIN(nblkrows_total, nblkcols_total)

         s = INT((d_max - 1)/d_min + 1)

         CALL mp_environ(numproc, iproc, mp_comm)
         CALL mp_environ(numproc, pdims, pcoord, mp_comm)

         IF (d_max == nblkrows_total) THEN
            split_rowcol = rowsplit
         ELSE
            split_rowcol = colsplit
         ENDIF

         DO WHILE (MOD(pdims(split_rowcol), s) .NE. 0)
            s = s - 1
         ENDDO

         SELECT CASE (split_rowcol)
         CASE (rowsplit)
            pdims_group = [pdims(1)/s, pdims(2)]
         CASE (colsplit)
            pdims_group = [pdims(1), pdims(2)/s]
         END SELECT

         igroup = pcoord(split_rowcol)/pdims_group(split_rowcol)

      ENDIF
      CALL dbcsr_r_split_rows_or_cols(split_info, mp_comm, s, igroup, split_rowcol)

      CALL timestop(handle)

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param info ...
!> \param mp_comm ...
!> \param nsplit ...
!> \param igroup ...
!> \param mp_comm_group ...
!> \param split_rowcol ...
!> \param pgrid_offset ...
! **************************************************************************************************
   SUBROUTINE dbcsr_r_get_split_info(info, mp_comm, nsplit, igroup, mp_comm_group, split_rowcol, pgrid_offset)
      TYPE(dbcsr_r_split_info), INTENT(IN)               :: info
      INTEGER, INTENT(OUT), OPTIONAL                     :: mp_comm, nsplit, igroup, mp_comm_group, &
                                                            split_rowcol
      INTEGER, DIMENSION(2), INTENT(OUT), OPTIONAL       :: pgrid_offset

      IF (PRESENT(mp_comm)) mp_comm = info%mp_comm
      IF (PRESENT(mp_comm_group)) mp_comm_group = info%mp_comm_group
      IF (PRESENT(split_rowcol)) split_rowcol = info%split_rowcol
      IF (PRESENT(igroup)) igroup = info%igroup
      IF (PRESENT(nsplit)) nsplit = info%ngroup

      IF (PRESENT(pgrid_offset)) THEN
         SELECT CASE (info%split_rowcol)
         CASE (rowsplit)
            pgrid_offset(:) = [info%igroup*info%pgrid_split_size, 0]
         CASE (colsplit)
            pgrid_offset(:) = [0, info%igroup*info%pgrid_split_size]
         END SELECT
      ENDIF

   END SUBROUTINE

! **************************************************************************************************
!> \brief split mpi grid by rows or columns
!> \param split_info ...
!> \param mp_comm global mpi communicator with a 2d cartesian grid
!> \param ngroup number of groups. Note: this number may be diminished on out in order to have a
!>        distribution where all groups except the last group have the same number of processes, and
!>        the last group has less processes than all other groups.
!> \param igroup my group ID
!> \param split_rowcol split rows or columns
! **************************************************************************************************
   SUBROUTINE dbcsr_r_split_rows_or_cols(split_info, mp_comm, ngroup, igroup, split_rowcol)
      TYPE(dbcsr_r_split_info), INTENT(OUT)              :: split_info
      INTEGER, INTENT(IN)                                :: mp_comm
      INTEGER, INTENT(INOUT)                             :: ngroup
      INTEGER, INTENT(IN)                                :: igroup, split_rowcol

      INTEGER :: igroup_check, iproc, iproc_group, iproc_group_check, mp_comm_group, numproc, &
                 numproc_group, numproc_group_check, handle
      INTEGER, DIMENSION(2)                              :: pcoord, pcoord_group, pdims, pdims_group
      LOGICAL                                            :: to_assert
      CHARACTER(LEN=*), PARAMETER :: routineN = 'dbcsr_r_split_rows_or_cols', &
                                     routineP = moduleN//':'//routineN

      CALL timeset(routineN, handle)

      CALL mp_comm_dup(mp_comm, split_info%mp_comm)
      split_info%igroup = igroup
      split_info%split_rowcol = split_rowcol

      CALL mp_comm_split_direct(mp_comm, mp_comm_group, igroup)

      CALL mp_environ(numproc, iproc, mp_comm)
      CALL mp_environ(numproc, pdims, pcoord, mp_comm)

      CALL mp_environ(numproc_group, iproc_group, mp_comm_group)

      IF (iproc == 0) THEN
         to_assert = MOD(numproc_group, pdims(MOD(split_rowcol, 2) + 1)) == 0
         DBCSR_ASSERT(to_assert)
         split_info%pgrid_split_size = numproc_group/pdims(MOD(split_rowcol, 2) + 1)
      ENDIF
      CALL mp_bcast(split_info%pgrid_split_size, 0, split_info%mp_comm)

      ngroup = (pdims(split_rowcol) + split_info%pgrid_split_size - 1)/split_info%pgrid_split_size
      split_info%ngroup = ngroup
      split_info%group_size = split_info%pgrid_split_size*pdims(MOD(split_rowcol, 2) + 1)
      split_info%pgrid_split_size_last = pdims(split_rowcol) - (split_info%ngroup - 1)*split_info%pgrid_split_size
      split_info%group_size_last = split_info%pgrid_split_size_last*pdims(MOD(split_rowcol, 2) + 1)

      IF (split_info%igroup < split_info%ngroup - 1) THEN
         split_info%pgrid_split_size_local = split_info%pgrid_split_size
         split_info%group_size_local = split_info%group_size
      ELSEIF (split_info%igroup == split_info%ngroup - 1) THEN
         split_info%pgrid_split_size_local = split_info%pgrid_split_size_last
         split_info%group_size_local = split_info%group_size_last
      ELSE
         DBCSR_ABORT("inconsistent number of MPI subgroups")
      ENDIF

      CALL world_to_group_proc_map(iproc, pdims, split_rowcol, split_info%pgrid_split_size, igroup_check, &
                                   split_info%pgrid_split_size_local, pdims_group, iproc_group)

      IF (igroup_check .NE. split_info%igroup) THEN
         DBCSR_ABORT('inconsistent subgroups')
      ENDIF

      CALL mp_cart_create(mp_comm_group, 2, pdims_group, pcoord_group, split_info%mp_comm_group)

      CALL mp_environ(numproc_group_check, iproc_group_check, split_info%mp_comm_group)

      DBCSR_ASSERT(iproc_group_check .EQ. iproc_group)

      CALL mp_comm_free(mp_comm_group)

      ALLOCATE (split_info%refcount)
      split_info%refcount = 1

      CALL timestop(handle)

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param dist_in ...
!> \param dist_out ...
!> \param own_dist ...
! **************************************************************************************************
   SUBROUTINE dbcsr_r_copy_distribution(dist_in, dist_out, own_dist)
      TYPE(dbcsr_r_distribution_type), INTENT(INOUT)     :: dist_in
      TYPE(dbcsr_r_distribution_type), INTENT(OUT)       :: dist_out
      LOGICAL, INTENT(IN), OPTIONAL                      :: own_dist

      LOGICAL                                            :: own_dist_prv

      IF (PRESENT(own_dist)) THEN
         own_dist_prv = own_dist
      ELSE
         own_dist_prv = .FALSE.
      ENDIF

      IF (.NOT. own_dist_prv) THEN
         CALL dbcsr_distribution_hold(dist_in%dbcsr_dist)
         CALL dbcsr_r_info_hold(dist_in%info)
      ENDIF

      dist_out = dist_in
   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param split_info ...
! **************************************************************************************************
   SUBROUTINE dbcsr_r_release_info(split_info)
      TYPE(dbcsr_r_split_info), INTENT(INOUT)            :: split_info

      LOGICAL                                            :: abort

      abort = .FALSE.

      IF (.NOT. ASSOCIATED(split_info%refcount)) THEN
         abort = .TRUE.
      ELSEIF (split_info%refcount < 1) THEN
         abort = .TRUE.
      ENDIF

      IF (abort) THEN
         DBCSR_ABORT("can not destroy non-existing split_info")
      ENDIF

      split_info%refcount = split_info%refcount - 1

      IF (split_info%refcount == 0) THEN
         CALL mp_comm_free(split_info%mp_comm_group)
         CALL mp_comm_free(split_info%mp_comm)
         DEALLOCATE (split_info%refcount)
      ENDIF
   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param split_info ...
! **************************************************************************************************
   SUBROUTINE dbcsr_r_info_hold(split_info)
      TYPE(dbcsr_r_split_info), INTENT(IN)               :: split_info

      INTEGER, POINTER                                   :: ref

      IF (split_info%refcount < 1) THEN
         DBCSR_ABORT("can not hold non-existing split_info")
      ENDIF
      ref => split_info%refcount
      ref = ref + 1
   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param info ...
! **************************************************************************************************
   SUBROUTINE dbcsr_r_print_info(info, name)
      TYPE(dbcsr_r_split_info), INTENT(IN)               :: info
      CHARACTER(len=*), INTENT(IN)                       :: name

      INTEGER                                            :: groupsize, igroup, io_unit, mp_comm, &
                                                            mp_comm_group, mynode, nsplit, &
                                                            numnodes, split_rowcol
      INTEGER, DIMENSION(2)                              :: coord, dims, groupcoord, groupdims, &
                                                            pgrid_offset

      CALL dbcsr_r_get_split_info(info, mp_comm, nsplit, igroup, mp_comm_group, split_rowcol, pgrid_offset)
      CALL mp_environ(numnodes, mynode, mp_comm)
      CALL mp_environ(numnodes, dims, coord, mp_comm)
      CALL mp_environ(groupsize, groupdims, groupcoord, mp_comm_group)

      io_unit = 0
      IF (mynode .EQ. 0) io_unit = default_output_unit

      IF (io_unit > 0) THEN
         SELECT CASE (split_rowcol)
         CASE (rowsplit)
            WRITE (io_unit, "(A,I4,1X,A,I4)") "DBCSR RECT | "//TRIM(name)//"   splitting rows by factor", nsplit
         CASE (colsplit)
            WRITE (io_unit, "(A,I4,1X,A,I4)") "DBCSR RECT | "//TRIM(name)//"   splitting columns by factor", nsplit
         END SELECT
         WRITE (io_unit, "(A, I4,A1,I4)") "DBCSR RECT | "//TRIM(name)//"   global grid sizes:", dims(1), "x", dims(2)
      ENDIF

      IF (io_unit > 0) THEN
         WRITE (io_unit, "(A,I4,A1,I4,A1,I4)") &
            "DBCSR RECT | "//TRIM(name)//"   grid sizes on groups 0, ...,", MAX(0, nsplit - 2), ":", &
            groupdims(1), "x", groupdims(2)
         SELECT CASE (split_rowcol)
         CASE (rowsplit)
            WRITE (io_unit, "(A,I4,A1,I4,A1,I4)") "DBCSR RECT | "//TRIM(name)//"   grid sizes on group", nsplit - 1, ":", &
               dims(1) - groupdims(1)*(nsplit - 1), "x", dims(2)
         CASE (colsplit)
            WRITE (io_unit, "(A,I4,A1,I4,A1,I4)") "DBCSR RECT | "//TRIM(name)//"   grid sizes on group", nsplit - 1, ":", &
               dims(1), "x", dims(2) - groupdims(2)*(nsplit - 1)
         END SELECT
      ENDIF

   END SUBROUTINE

! **************************************************************************************************
!> \brief map global process info to group
!> \param iproc global process ID
!> \param pdims global process dimensions
!> \param split_rowcol split rows or column
!> \param pgrid_split_size how many process rows/cols per group (except last group)
!> \param igroup group ID
!> \param pgrid_split_size_local how many process rows/cols in this group
!> \param pdims_group local process grid dimensions
!> \param iproc_group group local process ID
! **************************************************************************************************
   SUBROUTINE world_to_group_proc_map(iproc, pdims, split_rowcol, pgrid_split_size, igroup, pgrid_split_size_local, &
                                      pdims_group, iproc_group)
      INTEGER, INTENT(IN)                                :: iproc
      INTEGER, DIMENSION(2), INTENT(IN)                  :: pdims
      INTEGER, INTENT(IN)                                :: split_rowcol, pgrid_split_size
      INTEGER, INTENT(OUT)                               :: igroup
      INTEGER, INTENT(IN), OPTIONAL                      :: pgrid_split_size_local
      INTEGER, DIMENSION(2), INTENT(OUT), OPTIONAL       :: pdims_group
      INTEGER, INTENT(OUT), OPTIONAL                     :: iproc_group

      INTEGER, DIMENSION(2)                              :: pcoord, pcoord_group

      IF (PRESENT(pdims_group) .OR. PRESENT(iproc_group)) THEN
         DBCSR_ASSERT(PRESENT(pgrid_split_size_local))
      ENDIF

      IF (PRESENT(iproc_group)) THEN
         DBCSR_ASSERT(PRESENT(pdims_group))
      ENDIF

      pcoord = [iproc/pdims(2), MOD(iproc, pdims(2))]

      igroup = pcoord(split_rowcol)/pgrid_split_size

      SELECT CASE (split_rowcol)
      CASE (rowsplit)
         IF (PRESENT(pdims_group)) pdims_group = [pgrid_split_size_local, pdims(2)]
         IF (PRESENT(iproc_group)) pcoord_group = [MOD(pcoord(1), pgrid_split_size), pcoord(2)]
      CASE (colsplit)
         IF (PRESENT(pdims_group)) pdims_group = [pdims(1), pgrid_split_size_local]
         IF (PRESENT(iproc_group)) pcoord_group = [pcoord(1), MOD(pcoord(2), pgrid_split_size)]
      END SELECT
      IF (PRESENT(iproc_group)) iproc_group = pcoord_group(1)*pdims_group(2) + pcoord_group(2)
   END SUBROUTINE

! **************************************************************************************************
!> \brief map local process info to global info
!> \param iproc global process id
!> \param pdims global process grid dimensions
!> \param split_rowcol split rows or column
!> \param pgrid_split_size how many process rows/cols per group (except last group)
!> \param pgrid_split_size_local how many process rows/cols in this group
!> \param igroup group ID
!> \param iproc_group local process ID
! **************************************************************************************************
   SUBROUTINE group_to_world_proc_map(iproc, pdims, split_rowcol, pgrid_split_size, pgrid_split_size_local, &
                                      igroup, iproc_group)
      INTEGER, INTENT(OUT)                               :: iproc
      INTEGER, DIMENSION(2), INTENT(IN)                  :: pdims
      INTEGER, INTENT(IN)                                :: split_rowcol, pgrid_split_size, &
                                                            pgrid_split_size_local, igroup, &
                                                            iproc_group

      INTEGER, DIMENSION(2)                              :: pcoord, pcoord_group, pdims_group

      SELECT CASE (split_rowcol)
      CASE (rowsplit)
         pdims_group = [pgrid_split_size_local, pdims(2)]
      CASE (colsplit)
         pdims_group = [pdims(1), pgrid_split_size_local]
      END SELECT

      pcoord_group = [iproc_group/pdims_group(2), MOD(iproc_group, pdims_group(2))]

      SELECT CASE (split_rowcol)
      CASE (rowsplit)
         pcoord = [igroup*pgrid_split_size + pcoord_group(1), pcoord_group(2)]
      CASE (colsplit)
         pcoord = [pcoord_group(1), igroup*pgrid_split_size + pcoord_group(2)]
      END SELECT
      iproc = pcoord(1)*pdims(2) + pcoord(2)
   END SUBROUTINE

! **************************************************************************************************
!> \brief map group local block index to global matrix index
!> \param info ...
!> \param dist ...
!> \param row_group group local row block index
!> \param column_group group local column block index
!> \param row global block row
!> \param column global block column
! **************************************************************************************************
   SUBROUTINE block_index_local_to_global(info, dist, row_group, column_group, &
                                          row, column)
      TYPE(dbcsr_r_split_info), INTENT(IN)               :: info
      TYPE(dbcsr_r_distribution_type), INTENT(IN)        :: dist
      INTEGER, INTENT(IN), OPTIONAL                      :: row_group, column_group
      INTEGER(KIND=int_8), INTENT(OUT), OPTIONAL         :: row, column

      IF (PRESENT(row)) THEN
         DBCSR_ASSERT(PRESENT(row_group))
      ENDIF
      IF (PRESENT(column)) THEN
         DBCSR_ASSERT(PRESENT(column_group))
      ENDIF

      SELECT CASE (info%split_rowcol)
      CASE (rowsplit)
         ASSOCIATE (rows=>dist%local_rowcols)
            IF (PRESENT(row)) row = rows(row_group)
            IF (PRESENT(column)) column = column_group
         END ASSOCIATE
      CASE (colsplit)
         ASSOCIATE (cols=>dist%local_rowcols)
            IF (PRESENT(row)) row = row_group
            IF (PRESENT(column)) column = cols(column_group)
         END ASSOCIATE
      END SELECT
   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param info ...
!> \param dist ...
!> \param row ...
!> \param column ...
!> \param row_group ...
!> \param column_group ...
! **************************************************************************************************
   SUBROUTINE block_index_global_to_local(info, dist, row, column, row_group, column_group)
      TYPE(dbcsr_r_split_info), INTENT(IN)               :: info
      TYPE(dbcsr_r_distribution_type), INTENT(IN)        :: dist
      INTEGER(KIND=int_8), INTENT(IN), OPTIONAL          :: row, column
      INTEGER, INTENT(OUT), OPTIONAL                     :: row_group, column_group
      CHARACTER(LEN=*), PARAMETER :: routineN = 'block_index_global_to_local', &
                                     routineP = moduleN//':'//routineN

      IF (PRESENT(row_group)) THEN
         DBCSR_ASSERT(PRESENT(row))
      ENDIF
      IF (PRESENT(column_group)) THEN
         DBCSR_ASSERT(PRESENT(column))
      ENDIF

      SELECT CASE (info%split_rowcol)
      CASE (rowsplit)
         IF (PRESENT(row_group)) CALL mrowcol_to_group(dist, row, row_group)
         IF (PRESENT(column_group)) column_group = INT(column)
      CASE (colsplit)
         IF (PRESENT(row_group)) row_group = INT(row)
         IF (PRESENT(column_group)) CALL mrowcol_to_group(dist, column, column_group)
      END SELECT

   END SUBROUTINE

! **************************************************************************************************
!> \brief map matrix rows/columns to group local rows/columns
!> \param dist ...
!> \param rowcol global row/column to map
!> \param rowcol_group group local row/column
! **************************************************************************************************

   SUBROUTINE mrowcol_to_group(dist, rowcol, rowcol_group)

      CLASS(dbcsr_r_distribution_type), INTENT(IN) :: dist
      INTEGER(KIND=int_8), INTENT(IN) :: rowcol
      INTEGER, INTENT(OUT) :: rowcol_group
      CHARACTER(LEN=*), PARAMETER :: routineN = 'mrowcol_to_group', &
                                     routineP = moduleN//':'//routineN

      ASSOCIATE (rowcols=>dist%local_rowcols)
         rowcol_group = i8_bsearch(rowcols, rowcol)
      END ASSOCIATE
   END SUBROUTINE

   FUNCTION i8_isearch(array, el, l_index, u_index) result(res)
      INTEGER(KIND=int_8), intent(in) :: array(:)
      INTEGER(KIND=int_8), intent(in) :: el
      INTEGER, INTENT(in), OPTIONAL :: l_index, u_index
      integer :: res
      INTEGER(KIND=int_8) :: lindex, uindex, aindex

      lindex = 1_int_8
      uindex = INT(size(array), KIND=int_8)
      if (present(l_index)) lindex = INT(l_index, KIND=int_8)
      if (present(u_index)) uindex = INT(u_index, KIND=int_8)

      DO WHILE (lindex <= uindex)
         aindex = lindex + ((el - array(lindex))*(uindex - lindex)/(array(uindex) - array(lindex)))
         IF (array(aindex) < el) THEN
            lindex = aindex + 1
         ELSEIF (el < array(aindex)) THEN
            uindex = aindex - 1
         ELSE
            res = INT(aindex)
            RETURN
         END IF
      ENDDO
      res = INT(lindex)
   END FUNCTION

   FUNCTION i8_bsearch(array, el, l_index, u_index)&
        result(res)
     INTEGER(KIND=int_8), intent(in) :: array(:)
     INTEGER(KIND=int_8), intent(in) :: el
     INTEGER, INTENT(in), OPTIONAL :: l_index, u_index
     integer :: res

     character(len=*), parameter :: routineN='i8_bsearch',&
           routineP=moduleN//':'//routineN
     INTEGER :: lindex,uindex,aindex

     lindex=1
     uindex=size(array)
     if (present(l_index)) lindex=l_index
     if (present(u_index)) uindex=u_index
     DO WHILE (lindex <= uindex)
        aindex=(lindex+uindex)/2
        IF (array(aindex) < el) THEN
           lindex=aindex+1
        ELSE
           uindex=aindex-1
        END IF
     END DO
     res=lindex
   END FUNCTION

! **************************************************************************************************
!> \brief maps a group to matrix row/columns
!> \param info ...
!> \param rowcol_dist ...
!> \param igroup group ID
!> \param rowcols rows/ columns on this group
! **************************************************************************************************
   SUBROUTINE group_to_mrowcol(info, rowcol_dist, igroup, rowcols)
      TYPE(dbcsr_r_split_info), INTENT(IN)               :: info

      CLASS(dbcsr_r_distribution), INTENT(IN) :: rowcol_dist
      INTEGER, INTENT(IN) :: igroup
      INTEGER(KIND=int_8), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: rowcols
      INTEGER, DIMENSION(0:info%pgrid_split_size - 1) :: nrowcols_group
      INTEGER :: pcoord, nrowcols, count, pcoord_group
      INTEGER, DIMENSION(:), ALLOCATABLE :: sort_indices

      nrowcols_group(:) = 0
      DO pcoord = igroup*info%pgrid_split_size, (igroup + 1)*info%pgrid_split_size - 1
         pcoord_group = pcoord - igroup*info%pgrid_split_size
         nrowcols_group(pcoord_group) = SIZE(rowcol_dist%rowcols(pcoord))
      ENDDO
      nrowcols = SUM(nrowcols_group)

      ALLOCATE (rowcols(nrowcols))

      count = 0
      DO pcoord = igroup*info%pgrid_split_size, (igroup + 1)*info%pgrid_split_size - 1
         pcoord_group = pcoord - igroup*info%pgrid_split_size
         rowcols(count + 1:count + nrowcols_group(pcoord_group)) = rowcol_dist%rowcols(pcoord)
         count = count + nrowcols_group(pcoord_group)
      ENDDO

      ALLOCATE (sort_indices(nrowcols))
      CALL sort(rowcols, nrowcols, sort_indices)
   END SUBROUTINE

END MODULE
