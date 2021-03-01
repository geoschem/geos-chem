#ifdef BPCH_DIAG
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: initialize.F90
!
! !DESCRIPTION: Subroutine INITIALIZE does the following:
!
! \begin{enumerate}
! \item Zeroes globally defined GEOS-CHEM variables.
! \item Zeroes accumulating diagnostic arrays.
! \item Resets certain year/month/day and counter variables used
!       in GEOS-Chem diagnostic subroutines.
! \end{enumerate}
!
! !INTERFACE:
!
SUBROUTINE INITIALIZE( Input_Opt, State_Grid, IFLAG, RC )
!
! !USES:
!
  USE CMN_DIAG_MOD
  USE DIAG_MOD
  USE DIAG03_MOD
  USE DIAG53_MOD
  USE ErrCode_Mod
  USE ERROR_MOD
  USE Input_Opt_Mod,  ONLY : OptInput
  USE Precision_Mod
  USE State_Grid_Mod, ONLY : GrdState
  USE TIME_MOD

  IMPLICIT NONE
!
! !INPUT PARAMETERS:
!
  TYPE(OptInput), INTENT(IN)  :: Input_Opt   ! Input Options object
  TYPE(GrdState), INTENT(IN)  :: State_Grid  ! Grid State object

  ! If IFLAG=1, zero global CTM arrays
  ! If IFLAG=2, zero accumulating diagnostic arrays
  ! If IFLAG=3, zero accumulating diagnostic counters
  INTEGER,        INTENT(IN)  :: IFLAG
!
! !OUTPUT PARAMETERS:
!
  INTEGER,        INTENT(OUT) :: RC          ! Success or failure?
!
! !REMARKS:
!  Eventually we will fold this into "diag_mod.f" in a cleaner,
!   more consistent fashion.  Think about this later (bmy, 11/14/02)
!
! !REVISION HISTORY:
!  15 Jun 1998 - M. Prather  - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
  !=================================================================
  ! INITIALIZE begins here!
  !=================================================================

  ! Assume success
  RC = GC_SUCCESS

  ! Error condition if IFLAG does not equal 2, or 3!
  IF ( IFLAG < 2 .or. IFLAG > 3 ) THEN
     CALL ERROR_STOP( 'Invalid IFLAG!', 'initialize.f' )
  ENDIF

  !=================================================================
  ! If IFLAG=2 then zero the accumulating arrays
  !=================================================================
  IF ( IFLAG == 2 ) THEN

     ! Allocatable arrays are zeroed only if their
     ! respective diagnostics are turned on (bmy, 2/17/00)

#ifdef TOMAS
     !------------------------------
     ! TOMAS-specific diagnostics
     !------------------------------
     IF ( ND44 > 0 ) AD44 = 0e0 ! netcdf DRYDEP diag group

     IF ( ND59 > 0 ) THEN
        AD59_NUMB = 0e0
        AD59_SULF = 0e0
        AD59_SALT = 0e0
        AD59_ECIL = 0e0
        AD59_ECOB = 0e0
        AD59_OCIL = 0e0
        AD59_OCOB = 0e0
        AD59_DUST = 0e0
     ENDIF

     IF ( ND60 > 0 ) THEN
        AD60_COND  = 0e0
        AD60_COAG  = 0e0
        AD60_NUCL  = 0e0
        AD60_AQOX  = 0e0
        AD60_ERROR = 0e0
        AD60_SOA   = 0e0
     ENDIF

     IF ( ND61 > 0 ) AD61 = 0e0
#endif

#ifdef RRTMG
     !------------------------------
     ! RRTMG-specific diagnostics
     !------------------------------
     IF ( ND72 > 0 ) AD72 = 0e0
#endif

     ! For ND03 - mercury simulations (eck, sas, bmy, 1/20/05)
     IF ( ND03 > 0 ) THEN
        CALL ZERO_DIAG03
     ENDIF

     ! For ND53: POPS simulations (eck, 9/20/10)
     IF ( ND53 > 0 ) THEN
        CALL ZERO_DIAG53
     ENDIF

     ! Echo output
     IF ( Input_Opt%amIRoot ) THEN
        WRITE( 6, '(a)' ) '     - INITIALIZE: Diag arrays zeroed!'
     ENDIF
  ENDIF

  !=================================================================
  ! If IFLAG=3 then zero the counter variables & arrays
  !=================================================================
  IF ( IFLAG == 3 ) THEN

     ! Now reset timesteps here for now
     CALL SET_CT_CHEM( RESET=.TRUE. )
     CALL SET_CT_RAD ( RESET=.TRUE. )
     CALL SET_CT_CONV( RESET=.TRUE. )
     CALL SET_CT_DYN(  RESET=.TRUE. )
     CALL SET_CT_EMIS( RESET=.TRUE. )
     CALL SET_CT_DIAG( RESET=.TRUE. )

     ! Echo output
     IF ( Input_Opt%amIRoot ) THEN
        WRITE( 6, '(a)' ) '     - INITIALIZE: Diag counters zeroed!'
     ENDIF
  ENDIF

END SUBROUTINE INITIALIZE
!EOC
#endif
