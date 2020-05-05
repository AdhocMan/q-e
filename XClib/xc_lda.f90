!----------------------------------------------------------------------------
SUBROUTINE xc_lda_l( length, rho_in, ex_out, ec_out, vx_out, vc_out )
  !--------------------------------------------------------------------------
  !! LDA exchange and correlation functionals - Hartree a.u.
  !
  !! * Exchange:
  !!    * Slater;
  !!    * relativistic Slater.
  !! * Correlation:
  !!    * Ceperley-Alder (Perdew-Zunger parameters);
  !!    * Vosko-Wilk-Nusair;
  !!    * Lee-Yang-Parr;
  !!    * Perdew-Wang;
  !!    * Wigner;
  !!    * Hedin-Lundqvist;
  !!    * Ortiz-Ballone (Perdew-Zunger formula);
  !!    * Ortiz-Ballone (Perdew-Wang formula);
  !!    * Gunnarsson-Lundqvist.
  !
  !! NOTE:
  !! $$ E_x = \int E_x(\text{rho}) dr, E_x(\text{rho}) = 
  !!               \text{rho}\epsilon_c(\text{rho})\ . $$
  !! Same for correlation.
  !
  USE kind_l, ONLY: DP
  USE dft_par_mod
  
  USE exch_lda_l
  USE corr_lda_l
  !
  IMPLICIT NONE
  !
  INTEGER,  INTENT(IN) :: length
  !! length of the I/O arrays
  REAL(DP), INTENT(IN),  DIMENSION(length) :: rho_in
  !! Charge density
  REAL(DP), INTENT(OUT), DIMENSION(length) :: ex_out
  !! \(\epsilon_x(rho)\) ( NOT \(E_x(\text{rho})\) )
  REAL(DP), INTENT(OUT), DIMENSION(length) :: vx_out
  !! \(dE_x(\text{rho})/d\text{rho}\)  ( NOT \(d\epsilon_x(\text{rho})/d\text{rho}\) )
  REAL(DP), INTENT(OUT), DIMENSION(length) :: ec_out
  !! \(\epsilon_c(rho)\) ( NOT \(E_c(\text{rho})\) )
  REAL(DP), INTENT(OUT), DIMENSION(length) :: vc_out
  !! \(dE_c(\text{rho})/d\text{rho}\)  ( NOT \(d\epsilon_c(\text{rho})/d\text{rho}\) )
  !
  ! ... local variables
  !
  INTEGER  :: ir
  REAL(DP) :: rho, rs
  REAL(DP) :: ex, ec, ec_
  REAL(DP) :: vx, vc, vc_
  REAL(DP), PARAMETER :: third = 1.0_DP/3.0_DP, &
                         pi34 = 0.6203504908994_DP, e2 = 2.0_DP
  !                      pi34 = (3/4pi)^(1/3)
  !
#if defined(_OPENMP)
  INTEGER :: ntids
  INTEGER, EXTERNAL :: omp_get_num_threads
  !
  ntids = omp_get_num_threads()
#endif
  !
!$omp parallel if(ntids==1)
!$omp do private( rho, rs, ex, ec, ec_, vx, vc, vc_ )
  DO ir = 1, length
     !
     rho = ABS(rho_in(ir))
     !
     ! ... RHO THRESHOLD
     !
     IF ( rho > rho_threshold ) THEN
        rs = pi34 / rho**third
     ELSE
        ex_out(ir) = 0.0_DP  ;  ec_out(ir) = 0.0_DP
        vx_out(ir) = 0.0_DP  ;  vc_out(ir) = 0.0_DP
        CYCLE
     ENDIF
     !
     ! ... EXCHANGE
     !
     SELECT CASE( iexch )
     CASE( 1 )                      ! 'sla'
        !
        CALL slater_l( rs, ex, vx )
        !
     CASE( 2 )                      ! 'sl1'
        !
        CALL slater1_l( rs, ex, vx )
        !
     CASE( 3 )                      ! 'rxc'
        !
        CALL slater_rxc_l( rs, ex, vx )
        !
     CASE( 4, 5 )                   ! 'oep','hf'
        !
        IF ( exx_started ) THEN
           ex = 0.0_DP
           vx = 0.0_DP
        ELSE
           CALL slater_l( rs, ex, vx )
        ENDIF
        !
     CASE( 6, 7 )                   ! 'pb0x' or 'DF-cx-0', or 'DF2-0',
        !                           ! 'B3LYP'
        CALL slater_l( rs, ex, vx )
        IF ( exx_started ) THEN
           ex = (1.0_DP - exx_fraction) * ex
           vx = (1.0_DP - exx_fraction) * vx
        ENDIF
        !
     CASE( 8 )                      ! 'sla+kzk'
        !
        CALL slaterKZK_l( rs, ex, vx, finite_size_cell_volume )
        !
     CASE( 9 )                      ! 'X3LYP'
        !
        CALL slater_l( rs, ex, vx )
        IF ( exx_started ) THEN
           ex = (1.0_DP - exx_fraction) * ex
           vx = (1.0_DP - exx_fraction) * vx
        ENDIF
        !
     CASE DEFAULT
        !
        ex = 0.0_DP
        vx = 0.0_DP
        !
     END SELECT
     !
     !
     ! ... CORRELATION
     !
     SELECT CASE( icorr )
     CASE( 1 )
        !
        CALL pz_l( rs, 1, ec, vc )
        !
     CASE( 2 )
        !
        CALL vwn_l( rs, ec, vc )
        !
     CASE( 3 )
        !
        CALL lyp_l( rs, ec, vc )
        !
     CASE( 4 )
        !
        CALL pw_l( rs, 1, ec, vc )
        !
     CASE( 5 )
        !
        CALL wignerc_l( rs, ec, vc )
        !
     CASE( 6 )
        !
        CALL hl_l( rs, ec, vc )
        !
     CASE( 7 )
        !
        CALL pz_l( rs, 2, ec, vc )
        ! 
     CASE( 8 )
        !
        CALL pw_l( rs, 2, ec, vc )
        !
     CASE( 9 )
        !
        CALL gl_l( rs, ec, vc )
        !
     CASE( 10 )
        !
        CALL pzKZK_l( rs, ec, vc, finite_size_cell_volume )
        !
     CASE( 11 )
        !
        CALL vwn1_rpa_l( rs, ec, vc )
        !
     CASE( 12 )                ! 'B3LYP'
        !
        CALL vwn_l( rs, ec, vc )
        ec = 0.19_DP * ec
        vc = 0.19_DP * vc
        !
        CALL lyp_l( rs, ec_, vc_ )
        ec = ec + 0.81_DP * ec_
        vc = vc + 0.81_DP * vc_
        !
     CASE( 13 )                ! 'B3LYP-V1R'
        !
        CALL vwn1_rpa_l( rs, ec, vc )
        ec = 0.19_DP * ec
        vc = 0.19_DP * vc
        !
        CALL lyp_l( rs, ec_, vc_ )
        ec = ec + 0.81_DP * ec_
        vc = vc + 0.81_DP * vc_
        !
     CASE( 14 )                ! 'X3LYP'
        !
        CALL vwn1_rpa_l( rs, ec, vc )
        ec = 0.129_DP * ec
        vc = 0.129_DP * vc
        !
        CALL lyp_l( rs, ec_, vc_ )
        ec = ec + 0.871_DP * ec_
        vc = vc + 0.871_DP * vc_
        !
     CASE DEFAULT
        !
        ec = 0.0_DP
        vc = 0.0_DP
        !
     END SELECT
     !
     ex_out(ir) = ex  ;  ec_out(ir) = ec
     vx_out(ir) = vx  ;  vc_out(ir) = vc
     !
  ENDDO
!$omp end do
!$omp end parallel
  !
  !
  RETURN
  !
END SUBROUTINE xc_lda_l

