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
SUBROUTINE laxlib_distmat_refresh_z( mlocal, ndim, kdim, alpha, v, ldv, dm, idesc, beta,&
                                     w, ldw, irc_ip, nrc_ip, rank_ip)
 USE laxlib_parallel_include
 USE iso_c_binding
#if defined __SPLA
 USE spla
#endif
 USE mp,               ONLY : mp_bcast
 !
 IMPLICIT NONE
 INCLUDE 'laxlib_kinds.fh'
 INCLUDE 'laxlib_low.fh'
 INCLUDE 'laxlib_param.fh'
 !
 INTEGER, INTENT(IN)     :: mlocal
 INTEGER, INTENT(IN)     :: ndim
 INTEGER, INTENT(IN)     :: kdim
 COMPLEX(DP), INTENT(IN) :: alpha
 COMPLEX(DP), INTENT(IN), TARGET :: v(:,:)
 INTEGER, INTENT(IN) :: ldv
 COMPLEX(DP), INTENT(INOUT), TARGET :: dm( :, : )
 INTEGER, INTENT(IN) :: idesc(:)
 COMPLEX(DP), INTENT(IN) :: beta
 COMPLEX(DP), INTENT(IN), TARGET :: w(:,:)
 INTEGER, INTENT(IN) :: ldw
 INTEGER, INTENT(IN) :: irc_ip( : )
 INTEGER, INTENT(IN) :: nrc_ip( : )
 INTEGER, INTENT(IN) :: rank_ip( :, : )
 COMPLEX(DP), ALLOCATABLE :: vtmp( :, : )
 INTEGER :: ipc, ipr
 INTEGER :: nx, nr, nc, ir, ic, root, icc, ii
 INTEGER :: ierr
 INTEGER :: status_spla
 INTEGER :: ortho_parent_comm
 COMPLEX(DP) :: beta_loop
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
 status_spla = spla_pzgemm_sbs(mlocal, ndim, kdim, &
                    alpha, c_loc(v(1,1)), ldv, c_loc(dm(1,1)), nx, 0, 0, mat_dis_spla, &
                    beta, c_loc(w(1,1)), ldw, ctx_spla)
 IF( status_spla /= SPLA_SUCCESS ) &
   CALL errore( ' laxlib_compute_distmat ',' error when calling SPLA ', ABS(status_spla) )
 !
#else
 !
 ALLOCATE( vtmp( nx, nx ) )
 !
 DO ipc = 1, idesc(LAX_DESC_NPC)
    !
    nc = nrc_ip( ipc )
    ic = irc_ip( ipc )
    !
    IF( ic <= ndim ) THEN
       !
       nc = min( nc, ndim - ic + 1 )
       !
       beta_loop = beta
       DO ipr = 1, idesc(LAX_DESC_NPR)
          !
          nr = nrc_ip( ipr )
          ir = irc_ip( ipr )
          !
          IF( ir <= kdim ) THEN
             !
             nr = min( nr, kdim - ir + 1 )
             root = rank_ip( ipr, ipc )

             IF( ipr-1 == idesc(LAX_DESC_MYR) .AND. ipc-1 == idesc(LAX_DESC_MYC) &
                 .AND. idesc(LAX_DESC_ACTIVE_NODE) /= 0 ) THEN
                !
                !  this proc sends his block
                ! 
                CALL mp_bcast( dm(:,1:nc), root, ortho_parent_comm )
                CALL ZGEMM( 'N', 'N', mlocal, nc, nr, alpha, &
                         v(1,ir), ldv, dm, nx, beta_loop, w(1,ic), ldw )
             ELSE
                !
                !  all other procs receive
                ! 
                CALL mp_bcast( vtmp(:,1:nc), root, ortho_parent_comm )
                CALL ZGEMM( 'N', 'N', mlocal, nc, nr, alpha, &
                         v(1,ir), ldv, vtmp, nx, beta_loop, w(1,ic), ldw )
             END IF
             ! 
             beta_loop = ONE
          END IF

       END DO
       !
    END IF
    !
 END DO
 !
 DEALLOCATE( vtmp )
#endif
 !
 RETURN
END SUBROUTINE laxlib_distmat_refresh_z
