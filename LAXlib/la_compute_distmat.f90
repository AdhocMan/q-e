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
SUBROUTINE laxlib_compute_distmat_z( dm, kdim, alpha, v, ldv, w, ldw, idesc, &
                                     irc_ip, nrc_ip, rank_ip, nb1)
 USE laxlib_parallel_include
 USE iso_c_binding
#if defined __SPLA
 USE spla
#endif
 USE mp,               ONLY : mp_root_sum
 !
 IMPLICIT NONE
 INCLUDE 'laxlib_kinds.fh'
 INCLUDE 'laxlib_low.fh'
 INCLUDE 'laxlib_param.fh'
 !
 COMPLEX(DP), INTENT(INOUT) :: dm( :, : )
 INTEGER, INTENT(IN)     :: kdim
 COMPLEX(DP), INTENT(IN) :: alpha
 COMPLEX(DP), INTENT(IN) :: v(:,:)
 INTEGER, INTENT(IN) :: ldv
 COMPLEX(DP), INTENT(IN) :: w(:,:)
 INTEGER, INTENT(IN) :: ldw
 INTEGER, INTENT(IN) :: idesc(:)
 INTEGER, INTENT(IN) :: irc_ip( : )
 INTEGER, INTENT(IN) :: nrc_ip( : )
 INTEGER, INTENT(IN) :: rank_ip( :, : )
 INTEGER, INTENT(IN) :: nb1
 COMPLEX(DP), ALLOCATABLE :: vtmp( :, : )
 INTEGER :: ipc, ipr
 INTEGER :: nx, nr, nc, ir, ic, root, icc, ii
 INTEGER :: ierr
 INTEGER :: status_spla
 INTEGER :: ortho_parent_comm
 type(c_ptr) :: mat_dis_spla, ctx_spla
 !
 nx = idesc(LAX_DESC_NRCX)
 CALL laxlib_getval( ortho_parent_comm = ortho_parent_comm, ctx_spla = ctx_spla, &
                     mat_dis_spla = mat_dis_spla )
 !
#if defined __SPLA
 status_spla = spla_mat_dis_set_row_block_size(mat_dis_spla, idesc(LAX_DESC_NRCX))
 IF( status_spla /= SPLA_SUCCESS ) &
   CALL errore( ' laxlib_compute_distmat ',' error when calling SPLA ', ABS(status_spla) )
 status_spla = spla_mat_dis_set_col_block_size(mat_dis_spla, idesc(LAX_DESC_NRCX))
 IF( status_spla /= SPLA_SUCCESS ) &
   CALL errore( ' laxlib_compute_distmat ',' error when calling SPLA ', ABS(status_spla) )
 !
 ! status_spla = spla_pzgemm_ssbtr(idesc(LAX_DESC_N), idesc(LAX_DESC_N) - nb1 + 1, kdim, &
 !          SPLA_OP_CONJ_TRANSPOSE, alpha, v, ldv, w(:, nb1), ldw, ONE, dm, nx, 0, &
 !          nb1 - 1, SPLA_FILL_MODE_UPPER, mat_dis_spla, ctx_spla)
 status_spla = spla_pzgemm_ssbtr(idesc(LAX_DESC_N), idesc(LAX_DESC_N) - nb1 + 1, kdim, &
          SPLA_OP_CONJ_TRANSPOSE, alpha, c_loc(v_d(1,1)), ldv, c_loc(w_d(1, nb1)), ldw,&
          ONE, c_loc(dm(1,1)), nx, 0, nb1 - 1, SPLA_FILL_MODE_UPPER, mat_dis_spla, ctx_spla)
 IF( status_spla /= SPLA_SUCCESS ) &
   CALL errore( ' laxlib_compute_distmat ',' error when calling SPLA ', ABS(status_spla) )
 !
 CALL laxlib_zsqmher( idesc(LAX_DESC_N), dm, nx, idesc )
 !
#else
 !
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
#endif
 !
 RETURN
END SUBROUTINE laxlib_compute_distmat_z


SUBROUTINE laxlib_compute_distmat_z_gpu( dm, kdim, alpha, v_d, ldv, w_d, ldw,&
                                         idesc, irc_ip, nrc_ip, rank_ip, nb1)
 USE laxlib_parallel_include
 USE iso_c_binding
#if defined __SPLA
 USE spla
#endif
 USE mp,               ONLY : mp_root_sum
 !
 IMPLICIT NONE
 INCLUDE 'laxlib_kinds.fh'
 INCLUDE 'laxlib_low.fh'
 INCLUDE 'laxlib_param.fh'
 !
 COMPLEX(DP), INTENT(INOUT) :: dm( :, : )
 INTEGER, INTENT(IN)     :: kdim
 COMPLEX(DP), INTENT(IN) :: alpha
 COMPLEX(DP), INTENT(IN) :: v_d(:,:)
 INTEGER, INTENT(IN) :: ldv
 COMPLEX(DP), INTENT(IN) :: w_d(:,:)
#if defined(__CUDA)
   attributes(DEVICE)   :: v_d, w_d
#endif
 INTEGER, INTENT(IN) :: ldw
 INTEGER, INTENT(IN) :: idesc(:)
 INTEGER, INTENT(IN) :: irc_ip( : )
 INTEGER, INTENT(IN) :: nrc_ip( : )
 INTEGER, INTENT(IN) :: rank_ip( :, : )
 INTEGER, INTENT(IN) :: nb1
 COMPLEX(DP), ALLOCATABLE :: vtmp( :, : )
 INTEGER :: ipc, ipr
 INTEGER :: nx, nr, nc, ir, ic, root, icc, ii
 INTEGER :: ierr
 INTEGER :: status_spla
 INTEGER :: ortho_parent_comm
 type(c_ptr) :: mat_dis_spla, ctx_spla
 !
 nx = idesc(LAX_DESC_NRCX)
 CALL laxlib_getval( ortho_parent_comm = ortho_parent_comm, ctx_spla = ctx_spla, &
                     mat_dis_spla = mat_dis_spla )
 !
#if defined __SPLA
 status_spla = spla_mat_dis_set_row_block_size(mat_dis_spla, idesc(LAX_DESC_NRCX))
 IF( status_spla /= SPLA_SUCCESS ) &
   CALL errore( ' laxlib_compute_distmat ',' error when calling SPLA ', ABS(status_spla) )
 status_spla = spla_mat_dis_set_col_block_size(mat_dis_spla, idesc(LAX_DESC_NRCX))
 IF( status_spla /= SPLA_SUCCESS ) &
   CALL errore( ' laxlib_compute_distmat ',' error when calling SPLA ', ABS(status_spla) )
 !
 status_spla = spla_pzgemm_ssbtr(idesc(LAX_DESC_N), idesc(LAX_DESC_N) - nb1 + 1, kdim, &
          SPLA_OP_CONJ_TRANSPOSE, alpha, c_loc(v_d(1,1)), ldv, c_loc(w_d(1, nb1)), ldw,&
          ONE, c_loc(dm(1,1)), nx, 0, nb1 - 1, SPLA_FILL_MODE_UPPER, mat_dis_spla, ctx_spla)
 IF( status_spla /= SPLA_SUCCESS ) &
   CALL errore( ' laxlib_compute_distmat ',' error when calling SPLA ', ABS(status_spla) )
 !
 CALL laxlib_zsqmher( idesc(LAX_DESC_N), dm, nx, idesc )
 !
#else
 !
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

          CALL ZGEMM( 'C', 'N', nr, nc, kdim, alpha, v_d(1, ir), &
                      ldv, w_d(1,ii), ldw, ZERO, vtmp, nx )
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
#endif
 !
 RETURN
END SUBROUTINE laxlib_compute_distmat_z_gpu

