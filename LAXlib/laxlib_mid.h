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

INTERFACE update_distmat_new
SUBROUTINE update_distmat_c( dm, alpha, v, ldv, w, ldw, kdim, idesc, irc_ip, nrc_ip,&
                             rank_ip, nb1)
 IMPLICIT NONE
 INCLUDE 'laxlib_kinds.fh'
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
END SUBROUTINE
END INTERFACE

