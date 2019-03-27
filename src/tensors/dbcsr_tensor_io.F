!--------------------------------------------------------------------------------------------------!
! Copyright (C) by the DBCSR developers group - All rights reserved                                !
! This file is part of the DBCSR library.                                                          !
!                                                                                                  !
! For information on the license, see the LICENSE file.                                            !
! For further information please visit https://dbcsr.cp2k.org                                      !
! SPDX-License-Identifier: GPL-2.0+                                                                !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief DBCSR tensor Input / Output
!> \author Patrick Seewald
! **************************************************************************************************
MODULE dbcsr_tensor_io

#:include "dbcsr_tensor.fypp"
#:set maxdim = maxrank
#:set ndims = range(2,maxdim+1)

   USE dbcsr_tensor_types, ONLY: &
      dbcsr_t_get_info, dbcsr_t_type, ndims_tensor, dbcsr_t_get_num_blocks, dbcsr_t_get_data_size,&
      blk_dims_tensor
   USE dbcsr_kinds, ONLY: default_string_length, int_8, real_8
   USE dbcsr_mpiwrap, ONLY: mp_environ, mp_sum, mp_max

   IMPLICIT NONE
   PRIVATE
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbcsr_tensor_types'

   PUBLIC :: &
      dbcsr_t_write_tensor_info, &
      dbcsr_t_write_tensor_dist

CONTAINS

! **************************************************************************************************
!> \brief Write tensor global info: block dimensions, full dimensions and process grid dimensions
!> \param tensor ...
!> \param output_unit ...
!> \param full_info Whether to print distribution and block size vectors
! **************************************************************************************************
   SUBROUTINE dbcsr_t_write_tensor_info(tensor, output_unit, full_info)
      TYPE(dbcsr_t_type), INTENT(IN) :: tensor
      INTEGER, INTENT(IN)            :: output_unit
      LOGICAL, OPTIONAL, INTENT(IN)  :: full_info
      INTEGER, DIMENSION(ndims_tensor(tensor)) :: nblks_total, nfull_total, pdims, my_ploc, nblks_local, nfull_local
      INTEGER, DIMENSION(:), ALLOCATABLE :: ${varlist("blks_local")}$
      INTEGER, DIMENSION(:), ALLOCATABLE :: ${varlist("proc_dist")}$
      INTEGER, DIMENSION(:), ALLOCATABLE :: ${varlist("blk_size")}$
      INTEGER, DIMENSION(:), ALLOCATABLE :: ${varlist("blk_offset")}$
      CHARACTER(len=default_string_length)                   :: name
      INTEGER                            :: idim
      INTEGER                            :: iblk

      CALL dbcsr_t_get_info(tensor, nblks_total, nfull_total, nblks_local, nfull_local, pdims, my_ploc, &
                            ${varlist("blks_local")}$, ${varlist("proc_dist")}$, ${varlist("blk_size")}$, &
                            ${varlist("blk_offset")}$, &
                            name=name)

      IF (output_unit > 0) THEN
         WRITE (output_unit, "(T2,A)") &
            "GLOBAL INFO OF "//TRIM(name)
         WRITE (output_unit, "(T4,A,1X)", advance="no") "block dimensions:"
         DO idim = 1, ndims_tensor(tensor)
            WRITE (output_unit, "(I6)", advance="no") nblks_total(idim)
         ENDDO
         WRITE (output_unit, "(/T4,A,1X)", advance="no") "full dimensions:"
         DO idim = 1, ndims_tensor(tensor)
            WRITE (output_unit, "(I8)", advance="no") nfull_total(idim)
         ENDDO
         WRITE (output_unit, "(/T4,A,1X)", advance="no") "process grid dimensions:"
         DO idim = 1, ndims_tensor(tensor)
            WRITE (output_unit, "(I6)", advance="no") pdims(idim)
         ENDDO

         IF (PRESENT(full_info)) THEN
            IF (full_info) THEN
               WRITE (output_unit, '(/T4,A)', advance='no') "Block sizes:"
#:for dim in range(1, maxdim+1)
               IF (ndims_tensor(tensor) >= ${dim}$) THEN
                  WRITE (output_unit, '(/T8,A,1X,I1,A,1X)', advance='no') 'Dim', ${dim}$, ':'
                  DO iblk = 1, SIZE(blk_size_${dim}$)
                     WRITE (output_unit, '(I2,1X)', advance='no') blk_size_${dim}$ (iblk)
                  ENDDO
               ENDIF
#:endfor
               WRITE (output_unit, '(/T4,A)', advance='no') "Block distribution:"
#:for dim in range(1, maxdim+1)
               IF (ndims_tensor(tensor) >= ${dim}$) THEN
                  WRITE (output_unit, '(/T8,A,1X,I1,A,1X)', advance='no') 'Dim', ${dim}$, ':'
                  DO iblk = 1, SIZE(proc_dist_${dim}$)
                     WRITE (output_unit, '(I3,1X)', advance='no') proc_dist_${dim}$ (iblk)
                  ENDDO
               ENDIF
#:endfor
            ENDIF
            WRITE (output_unit, *)
         ENDIF
      ENDIF

   END SUBROUTINE

! **************************************************************************************************
!> \brief Write info on tensor distribution & load balance
!> \param tensor ...
!> \param output_unit ...
! **************************************************************************************************
   SUBROUTINE dbcsr_t_write_tensor_dist(tensor, output_unit)
      TYPE(dbcsr_t_type), INTENT(IN) :: tensor
      INTEGER, INTENT(IN)            :: output_unit
      INTEGER                        :: nproc, myproc, nblock_max, nelement_max
      INTEGER(KIND=int_8)            :: nblock_sum, nelement_sum, nblock_tot
      INTEGER, ALLOCATABLE, DIMENSION(:) :: nblock, nelement
      INTEGER                        :: mp_comm
      INTEGER, DIMENSION(2)          :: tmp
      INTEGER, DIMENSION(ndims_tensor(tensor)) :: bdims
      REAL(KIND=real_8)              :: occupation

      mp_comm = tensor%pgrid%mp_comm_2d
      CALL mp_environ(nproc, myproc, mp_comm)

      ALLOCATE (nblock(0:nproc - 1))
      nblock(:) = 0

      ALLOCATE (nelement(0:nproc - 1))
      nelement(:) = 0

      nblock(myproc) = dbcsr_t_get_num_blocks(tensor)
      nelement(myproc) = dbcsr_t_get_data_size(tensor)

      nblock_sum = nblock(myproc)
      nblock_max = nblock(myproc)
      nelement_sum = nelement(myproc)
      nelement_max = nelement(myproc)
      CALL mp_sum(nblock_sum, mp_comm)
      CALL mp_sum(nelement_sum, mp_comm)
      tmp = (/nblock_max, nelement_max/)
      CALL mp_max(tmp, mp_comm)
      nblock_max = tmp(1); nelement_max = tmp(2)

      CALL blk_dims_tensor(tensor, bdims)
      nblock_tot = PRODUCT(INT(bdims, KIND=int_8))

      occupation = -1.0_real_8
      IF (nblock_tot .NE. 0) occupation = 100.0_real_8*REAL(nblock_sum, real_8)/REAL(nblock_tot, real_8)

      IF (output_unit > 0) THEN
         WRITE (output_unit, "(T2,A)") &
            "DISTRIBUTION OF "//TRIM(tensor%name)
         WRITE (output_unit, "(T15,A,T68,I13)") "Number of non-zero blocks:", nblock_sum
         WRITE (output_unit, "(T15,A,T75,F6.2)") "Percentage of non-zero blocks:", occupation
         WRITE (output_unit, "(T15,A,T68,I13)") "Average number of blocks per CPU:", (nblock_sum + nproc - 1)/nproc
         WRITE (output_unit, "(T15,A,T68,I13)") "Maximum number of blocks per CPU:", nblock_max
         WRITE (output_unit, "(T15,A,T68,I13)") "Average number of matrix elements per CPU:", (nelement_sum + nproc - 1)/nproc
         WRITE (output_unit, "(T15,A,T68,I13)") "Maximum number of matrix elements per CPU:", nelement_max
      ENDIF

   END SUBROUTINE

END MODULE
