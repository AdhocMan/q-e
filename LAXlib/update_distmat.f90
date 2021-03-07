!
! Copyright (C) 2003-2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
#define ZERO ( 0.D0, 0.D0 )
#define ONE  ( 1.D0, 0.D0 )
!
SUBROUTINE update_distmat_c( dm, alpha, v, ldv, w, ldw, kdim, idesc, irc_ip, nrc_ip,&
                             rank_ip, nb1)
 USE laxlib_parallel_include
 ! USE util_param,    ONLY : DP
  USE mp_bands_util,    ONLY : intra_bgrp_comm, inter_bgrp_comm, root_bgrp_id,&
                               nbgrp, my_bgrp_id
  USE mp,               ONLY : mp_bcast, mp_root_sum, mp_sum, mp_barrier, &
                               mp_size, mp_type_free, mp_allgather
 !
 IMPLICIT NONE
 INCLUDE 'laxlib_kinds.fh'
 INCLUDE 'laxlib_low.fh'
 INCLUDE 'laxlib_param.fh'
 !
 COMPLEX(DP), INTENT(IN) :: alpha
 COMPLEX(DP), INTENT(INOUT) :: dm( :, : )
 COMPLEX(DP), INTENT(IN) :: v(:,:)
 INTEGER, INTENT(IN) :: ldv
 COMPLEX(DP), INTENT(IN) :: w(:,:)
 INTEGER, INTENT(IN) :: ldw
 INTEGER, INTENT(IN)     :: kdim
 INTEGER, INTENT(IN) :: idesc(:)
 INTEGER, INTENT(IN) :: irc_ip( : )
 INTEGER, INTENT(IN) :: nrc_ip( : )
 INTEGER, INTENT(IN) :: rank_ip( :, : )
 INTEGER, INTENT(IN) :: nb1
 COMPLEX(DP), ALLOCATABLE :: vtmp( :, : )
 INTEGER :: ipc, ipr
 INTEGER :: nx, nr, nc, ir, ic, root, icc, ii
 INTEGER :: ortho_parent_comm

 nx = idesc(LAX_DESC_NRCX)
 CALL laxlib_getval( ortho_parent_comm = ortho_parent_comm )

 ALLOCATE( vtmp( nx, nx ) )
 !
 vtmp = ZERO
 !
 DO ipc = 1, idesc(LAX_DESC_NPC)
    !
    nc = nrc_ip( ipc )
    ic = irc_ip( ipc )
    !
    IF( ic+nc-1 >= nb1 ) THEN
       !
       nc = MIN( nc, ic+nc-1 - nb1 + 1 )
       IF( ic >= nb1 ) THEN
          ii = ic
          icc = 1
       ELSE
          ii = nb1
          icc = nb1-ic+1
       END IF
       !
       ! icc to nc is the local index of the unconverged bands
       ! ii is the global index of the first unconverged bands
       !
       DO ipr = 1, ipc ! idesc(LAX_DESC_NPR) use symmetry
          !
          nr = nrc_ip( ipr )
          ir = irc_ip( ipr )
          !
          root = rank_ip( ipr, ipc )

          CALL ZGEMM( 'C', 'N', nr, nc, kdim, alpha, v(1, ir), &
                      ldv, w(1,ii), ldw, ZERO, vtmp, nx )
          ! IF (ortho_parent_comm.ne.intra_bgrp_comm .and. nbgrp > 1) vtmp = vtmp/nbgrp
          !
          IF(  (idesc(LAX_DESC_ACTIVE_NODE) > 0) .AND. &
               (ipr-1 == idesc(LAX_DESC_MYR)) .AND. (ipc-1 == idesc(LAX_DESC_MYC)) ) THEN
             CALL mp_root_sum( vtmp(:,1:nc), dm(:,icc:icc+nc-1), root, ortho_parent_comm )
          ELSE
             CALL mp_root_sum( vtmp(:,1:nc), dm, root, ortho_parent_comm )
          END IF

       END DO
       !
    END IF
    !
 END DO
 !
 CALL laxlib_zsqmher( idesc(LAX_DESC_N), dm, nx, idesc )
 !
 DEALLOCATE( vtmp )
 RETURN
END SUBROUTINE update_distmat_c
