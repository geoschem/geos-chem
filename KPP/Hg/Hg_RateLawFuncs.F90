!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: fullchem_RateLawFuncs
!
! !DESCRIPTION: Provides rate-law functions used by the "fullchem" chemical
!  mechanism.  This will be referenced from within subroutine Update_RCONST.
!\\
!\\
! !INTERFACE:
!
MODULE Hg_RateLawFuncs
!
! !USES:
!
  USE gckpp_Global
  USE gckpp_Parameters
  USE gckpp_Precision
  USE rateLawUtilFuncs

  IMPLICIT NONE
  PUBLIC
!EOP
!------------------------------------------------------------------------------
!BOC
CONTAINS

  !#########################################################################
  !#####          RATE-LAW FUNCTIONS FOR GAS-PHASE REACTIONS           #####
  !#####   Some common functions are defined in rateLawUtilFuncs.F90   #####
  !#########################################################################

  FUNCTION GCJPLPR_abab( a1, b1, a2, b2, fv ) RESULT( k )
    !
    ! Third body effect for pressure dependence of rate coefficients.
    ! a1, b1 are the Arrhenius parameters for the lower-limit rate.
    ! a2, b2 are the Arrhenius parameters for the upper-limit rate.
    ! fv     is the falloff curve paramter, (see ATKINSON ET. AL (1992)
    !        J. Phys. Chem. Ref. Data 21, P. 1145). Usually fv = 0.6.
    !
    ! For these reactions, these Arrhenius law terms evaluate to 1:
    !    EXP(c1/T)
    !    EXP(c2/T)
    ! because c1 = c2 = 0.  Therefore we can skip computing these
    ! terms.  Also, fct1 = fct2 = 0, so we will skip computing these
    ! terms as well.  This is more computationally efficient.
    ! (bmy, 06 Jan 2022)
    !
    REAL(dp), INTENT(IN) :: a1,   b1,    a2,    b2,   fv
    REAL(dp)             :: rlow, rhigh, xyrat, blog, fexp, k
    !
    rlow  = a1 * ( K300_OVER_TEMP**b1 ) * NUMDEN
    rhigh = a2 * ( K300_OVER_TEMP**b2 )
    xyrat = rlow / rhigh
    blog  = LOG10( xyrat )
    fexp  = 1.0_dp / ( 1.0_dp + ( blog * blog ) )
    k     = rlow * ( fv**fexp ) / ( 1.0_dp + xyrat )
  END FUNCTION GCJPLPR_abab

  !#########################################################################
  !#####        RATE-LAW FUNCTIONS FOR HETEROGENEOUS REACTIONS         #####
  !#####   Some common functions are defined in rateLawUtilFuncs.F90   #####
  !#########################################################################

  FUNCTION CloudHet_Hg( H, srMw, gamma ) RESULT( k )
    !
    ! Calculates the loss frequency (1/s) of gas species due to
    ! heterogeneous chemistry on liquid clouds in a partially cloudy
    ! grid cell. The function uses the "entrainment limited uptake"
    ! equations of Holmes et al. (2019).  Modified from C.Holmes's
    ! CloudHet routine by Viral Shah (Oct 2020).
    !
    ! Reference:
    ! Holmes, C.D., Bertram, T. H., Confer, K. L., Ronan, A. C., Wirks,
    !   C. K., Graham, K. A., Shah, V. (2019) The role of clouds in the
    !   tropospheric NOx cycle: a new modeling approach for cloud chemistry
    !   and its global implications, Geophys. Res. Lett. 46, 4980-4990,
    !   https://doi.org/10.1029/2019GL081990
    !
    TYPE(HetState), INTENT(IN) :: H       ! Hetchem State object
    REAL(dp),       INTENT(IN) :: gamma   ! Reaction probability [1]
    REAL(dp),       INTENT(IN) :: srMw    ! SQRT( mol wt in g/mol )
    REAL(dp)                   :: K       ! Reaction rate [1/s]
    !
    REAL(dp)                   :: kic, kIinv, kEinv
    !
    ! If cloud fraction < 0.0001 (0.01%) or there is zero cloud surface
    ! area, then return zero uptake
    IF ( ( H%CldFr < 0.0001_dp ) .or. ( H%aLiq + H%aIce <= 0.0_dp ) ) THEN
       k = 0.0_dp
       RETURN
    ENDIF

    !-----------------------------------------------------------------------
    ! In-cloud loss frequency, 1/s
    !-----------------------------------------------------------------------
    kic = Ars_L1K( H%aLiq, H%rLiq, gamma, srMw )

    !------------------------------------------------------------------------
    ! Grid-average loss frequency;
    ! Add in-cloud and entrainment rates in series
    !
    ! APPROXIMATE expression for entrainment-limited uptake
    ! Approximation error in loss frequency is typically <2% and always <50%.
    !------------------------------------------------------------------------

    ! Entrainment rate, inverse [s]
    ! (Residence time of air in clouds, = 3600 s)
    kEinv = SafeDiv( ( 3600.0_dp * H%ClearFr ) , H%CldFr, 1e+30_dp )

    ! In-cloud loss rate, inverse [s]
    kIinv = SafeDiv( 1.0_dp, ( H%CldFr * kic ), 1e+30_dp )

    ! Overall heterogeneous loss rate, grid average, 1/s
    k     = SafeDiv( 1.0_dp, ( kEinv + kIinv ), 0.0_dp )

  END FUNCTION CloudHet_Hg

  FUNCTION Het_HgIIP_Org( H, gamma, srMw ) RESULT( k )
    !
    ! Computes the heterogeneous chemistry reaction rate [1/s]
    ! for species forming organic HgIIP aerosol in liquid clouds.
    !
    TYPE(HetState), INTENT(IN) :: H       ! Hetchem State object
    REAL(dp),       INTENT(IN) :: gamma   ! Reaction probability [1]
    REAL(dp),       INTENT(IN) :: srMw    ! SQRT( mol wt in g/mol )
    REAL(dp)                   :: K       ! Reaction rate [1/s]
    !
    k = 0.0_dp
    !
    ! Reaction only takes place in cloudy grid boxes
    IF ( H%cloudBox ) THEN
       k = CloudHet_Hg( H, gamma, srMw ) * H%fracOrgAer
    ENDIF
  END FUNCTION Het_HgIIP_Org

  FUNCTION Het_HgIIP_Inorg( H, gamma, srMw ) RESULT( k )
    !
    ! Computes the heterogeneous chemistry reaction rate [1/s]
    ! for species forming inorganic HgIIP aerosol in liquid clouds
    !
    TYPE(HetState), INTENT(IN) :: H       ! Hetchem State object
    REAL(dp),       INTENT(IN) :: gamma   ! Reaction probability [1]
    REAL(dp),       INTENT(IN) :: srMw    ! SQRT( mol wt in g/mol )
    REAL(dp)                   :: K       ! Reaction rate [1/s]
    !
    k = 0.0_dp
    !
    ! Reaction only takes place in cloudy grid boxes
    IF ( H%cloudBox ) THEN
       k = CloudHet_Hg( H, gamma, srMw ) * H%fracInorgAer
    ENDIF
  END FUNCTION Het_HgIIP_Inorg

END MODULE Hg_RateLawFuncs
!EOC
