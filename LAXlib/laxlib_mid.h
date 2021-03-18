!
! Copyright (C) 2003-2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!

INTERFACE laxlib_pdsyevd
SUBROUTINE laxlib_pdsyevd_x( tv, n, idesc, hh, ldh, e )
         IMPLICIT NONE
         include 'laxlib_param.fh'
         include 'laxlib_kinds.fh'
         LOGICAL, INTENT(IN) :: tv
         INTEGER, INTENT(IN) :: n, ldh
         INTEGER, INTENT(IN) :: idesc(LAX_DESC_SIZE)
         REAL(DP) :: hh( ldh, ldh )
         REAL(DP) :: e( n )
END SUBROUTINE
END INTERFACE

INTERFACE laxlib_pzheevd
SUBROUTINE laxlib_pzheevd_x( tv, n, idesc, hh, ldh, e )
   IMPLICIT NONE
   include 'laxlib_param.fh'
   include 'laxlib_kinds.fh'
   LOGICAL, INTENT(IN) :: tv
   INTEGER, INTENT(IN) :: n, ldh
   INTEGER, INTENT(IN) :: idesc(LAX_DESC_SIZE)
   COMPLEX(DP) :: hh( ldh, ldh )
   REAL(DP) :: e( n )
END SUBROUTINE
END INTERFACE

INTERFACE laxlib_pzpotrf
SUBROUTINE laxlib_pzpotrf_x( sll, ldx, n, idesc )
   implicit none
   include 'laxlib_param.fh'
   include 'laxlib_kinds.fh'
   integer :: n, ldx
   integer, INTENT(IN) :: idesc(LAX_DESC_SIZE)
   complex(DP) :: sll( ldx, ldx )
END SUBROUTINE
END INTERFACE

INTERFACE laxlib_pdpotrf
SUBROUTINE laxlib_pdpotrf_x( sll, ldx, n, idesc )
   implicit none
   include 'laxlib_param.fh'
   include 'laxlib_kinds.fh'
   integer  :: n, ldx
   INTEGER, INTENT(IN) :: idesc(LAX_DESC_SIZE)
   REAL(DP) :: sll( ldx, ldx )
END SUBROUTINE
END INTERFACE

INTERFACE laxlib_pztrtri
SUBROUTINE laxlib_pztrtri_x ( sll, ldx, n, idesc )
   implicit none
   include 'laxlib_param.fh'
   include 'laxlib_kinds.fh'
   INTEGER, INTENT( IN ) :: n, ldx
   INTEGER, INTENT(IN) :: idesc(LAX_DESC_SIZE)
   COMPLEX(DP), INTENT( INOUT ) :: sll( ldx, ldx )
END SUBROUTINE
END INTERFACE

INTERFACE laxlib_pdtrtri
SUBROUTINE laxlib_pdtrtri_x ( sll, ldx, n, idesc )
   implicit none
   include 'laxlib_param.fh'
   include 'laxlib_kinds.fh'
   INTEGER, INTENT( IN ) :: n, ldx
   INTEGER, INTENT(IN) :: idesc(LAX_DESC_SIZE)
   REAL(DP), INTENT( INOUT ) :: sll( ldx, ldx )
END SUBROUTINE
END INTERFACE

INTERFACE laxlib_compute_distmat
SUBROUTINE laxlib_compute_distmat_z( dm, kdim, alpha, v, ldv, w, ldw,&
                                     idesc, irc_ip, nrc_ip, rank_ip, nb1)
 IMPLICIT NONE
 INCLUDE 'laxlib_kinds.fh'
 COMPLEX(DP), INTENT(INOUT) :: dm( :, : )
 INTEGER, INTENT(IN) :: kdim
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
END SUBROUTINE
END INTERFACE

INTERFACE laxlib_compute_distmat_gpu
SUBROUTINE laxlib_compute_distmat_gpu_z( dm, kdim, alpha, v_d, ldv, w_d, ldw,&
                                     idesc, irc_ip, nrc_ip, rank_ip, nb1)
 IMPLICIT NONE
 INCLUDE 'laxlib_kinds.fh'
 COMPLEX(DP), INTENT(INOUT) :: dm( :, : )
 INTEGER, INTENT(IN) :: kdim
 COMPLEX(DP), INTENT(IN) :: alpha
 COMPLEX(DP), INTENT(IN) :: v_d(:,:)
 INTEGER, INTENT(IN) :: ldv
 COMPLEX(DP), INTENT(IN) :: w_d(:,:)
 INTEGER, INTENT(IN) :: ldw
 INTEGER, INTENT(IN) :: idesc(:)
 INTEGER, INTENT(IN) :: irc_ip( : )
 INTEGER, INTENT(IN) :: nrc_ip( : )
 INTEGER, INTENT(IN) :: rank_ip( :, : )
 INTEGER, INTENT(IN) :: nb1
#if defined(__CUDA)
 attributes(DEVICE)   :: v_d, w_d
#endif
END SUBROUTINE
END INTERFACE

INTERFACE laxlib_distmat_refresh
SUBROUTINE laxlib_distmat_refresh_z( mdim, ndim, kdim, alpha, v, ldv, dm, idesc, beta,&
                                     w, ldw, irc_ip, nrc_ip, rank_ip)
 IMPLICIT NONE
 INCLUDE 'laxlib_kinds.fh'
 INTEGER, INTENT(IN)     :: mdim
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
END SUBROUTINE
END INTERFACE

INTERFACE laxlib_distmat_refresh_gpu
SUBROUTINE laxlib_distmat_refresh_gpu_z( mdim, ndim, kdim, alpha, v_d, ldv, dm, idesc, beta,&
                                     w_d, ldw, irc_ip, nrc_ip, rank_ip)
 IMPLICIT NONE
 INCLUDE 'laxlib_kinds.fh'
 INTEGER, INTENT(IN)     :: mdim
 INTEGER, INTENT(IN)     :: ndim
 INTEGER, INTENT(IN)     :: kdim
 COMPLEX(DP), INTENT(IN) :: alpha
 COMPLEX(DP), INTENT(IN), TARGET :: v_d(:,:)
 INTEGER, INTENT(IN) :: ldv
 COMPLEX(DP), INTENT(INOUT), TARGET :: dm( :, : )
 INTEGER, INTENT(IN) :: idesc(:)
 COMPLEX(DP), INTENT(IN) :: beta
 COMPLEX(DP), INTENT(IN), TARGET :: w_d(:,:)
 INTEGER, INTENT(IN) :: ldw
 INTEGER, INTENT(IN) :: irc_ip( : )
 INTEGER, INTENT(IN) :: nrc_ip( : )
 INTEGER, INTENT(IN) :: rank_ip( :, : )
#if defined(__CUDA)
 attributes(DEVICE)   :: v_d, w_d
#endif
END SUBROUTINE
END INTERFACE

