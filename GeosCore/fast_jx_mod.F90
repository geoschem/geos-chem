!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: fast_jx_mod.F90
!
! !DESCRIPTION: Module FAST\_JX\_MOD contains routines and variables for
!  calculating photolysis rates using the Fast-JX scheme (Prather et al).
!  Current implementation is version 7.0a.
!\\
!\\
! !INTERFACE:
!
MODULE FAST_JX_MOD
!
! !USES:
!
  USE CMN_FJX_MOD
  USE PRECISION_MOD    ! For GEOS-Chem Precision (fp)

  IMPLICIT NONE

  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC  :: EXITC
  PUBLIC  :: PHOTO_JX
  PUBLIC  :: INIT_FJX
  PUBLIC  :: FAST_JX
  PUBLIC  :: PHOTRATE_ADJ
!
! !PRIVATE MEMBER FUNCTIONS:
!
  PRIVATE :: SOLAR_JX
  PRIVATE :: OPMIE
  PRIVATE :: MIESCT
  PRIVATE :: LEGND0
  PRIVATE :: BLKSLV
  PRIVATE :: GEN_ID
  PRIVATE :: JRATET
  PRIVATE :: X_INTERP
  PRIVATE :: SPHERE2
  PRIVATE :: EXTRAL
  PRIVATE :: RD_PROF_NC
  PRIVATE :: RD_AOD
  PRIVATE :: RD_MIE
  PRIVATE :: RD_XXX
  PRIVATE :: RD_JS_JX
  PRIVATE :: SET_AER
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
  ! Flags for certain photo-reactions that will be adjusted by
  ! subroutine PHOTRATE_ADJ, which is called by FlexChem (bmy 3/29/16)
  INTEGER, PUBLIC :: RXN_O2    = -1   ! O2  + jv --> O   + O
  INTEGER, PUBLIC :: RXN_O3_1  = -1   ! O3  + hv --> O2  + O
  INTEGER, PUBLIC :: RXN_O3_2a = -1   ! O3  + hv --> 2OH         (Tropchem)
  ! O3  + hv --> O2  + O(1D) (UCX #1)
  INTEGER, PUBLIC :: RXN_O3_2b = -1   ! O3  + hv --> O2  + O(1D) (UCX #2)
  INTEGER, PUBLIC :: RXN_H2SO4 = -1   ! SO4 + hv --> SO2 + 2OH
  INTEGER, PUBLIC :: RXN_NO2   = -1   ! NO2 + hv --> NO  + O

  INTEGER, PUBLIC :: RXN_JHNO3  = -1   ! HNO3 + hv --> OH + NO2
  INTEGER, PUBLIC :: RXN_JNITSa = -1   ! NITs  + hv --> HNO2
  INTEGER, PUBLIC :: RXN_JNITSb = -1   ! NITs  + hv --> NO2
  INTEGER, PUBLIC :: RXN_JNITa  = -1   ! NIT + hv --> HNO2
  INTEGER, PUBLIC :: RXN_JNITb  = -1   ! NIT + hv --> NO2

  ! Needed for UCX_MOD
  INTEGER, PUBLIC :: RXN_NO    = -1
  INTEGER, PUBLIC :: RXN_NO3   = -1
  INTEGER, PUBLIC :: RXN_N2O   = -1

  ! Species ID flags
  INTEGER :: id_CH2IBr, id_IBr,  id_CH2ICl, id_ICl,   id_I2
  INTEGER :: id_HOI,    id_IO,   id_OIO,    id_INO,   id_IONO
  INTEGER :: id_IONO2,  id_I2O2, id_CH3I,   id_CH2I2, id_I2O4
  INTEGER :: id_I2O3

  ! Needed for scaling JNIT/JNITs photolysis to JHNO3
  REAL(fp)      :: JscaleNITs, JscaleNIT, JNITChanA, JNITChanB

#ifdef MODEL_GEOS
  ! Diagnostics arrays (ckeller, 5/22/18)
  REAL, ALLOCATABLE, PUBLIC :: EXTRAL_NLEVS(:,:), EXTRAL_NITER(:,:)
#endif

CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: fast_jx
!
! !DESCRIPTION: Subroutine FAST\_JX loops over longitude and latitude, and
!  calls PHOTO\_JX to compute J-Values for each column at every chemistry
!  time-step.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE FAST_JX( WLAOD, Input_Opt, State_Chm, State_Diag, &
                      State_Grid, State_Met, RC )
!
! !USES:
!
    USE CMN_SIZE_MOD,       ONLY : NDUST, NRH
    USE ErrCode_Mod
    USE ERROR_MOD,          ONLY : ERROR_STOP, ALLOC_ERR
    USE ERROR_MOD,          ONLY : DEBUG_MSG
    USE Input_Opt_Mod,      ONLY : OptInput
    USE State_Chm_Mod,      ONLY : ChmState
    USE State_Chm_Mod,      ONLY : Ind_
    USE State_Diag_Mod,     ONLY : DgnState
    USE State_Grid_Mod,     ONLY : GrdState
    USE State_Met_Mod,      ONLY : MetState
    USE TIME_MOD,           ONLY : GET_MONTH, GET_DAY, GET_DAY_OF_YEAR
    USE TIME_MOD,           ONLY : GET_TAU,   GET_YEAR
    USE TOMS_MOD,           ONLY : GET_OVERHEAD_O3

    IMPLICIT NONE

!==============================================================================
! Uncomment the appropriate #define statement to denote which of the
! available cloud overlap options that you wish to use.

!! Linear overlap
!#define USE_LINEAR_OVERLAP 1

! Approximate random overlap (balance between accuracy & speed)
#define USE_APPROX_RANDOM_OVERLAP 1

!! Maximum random cloud overlap (most computationally intensive)
!#define USE_MAXIMUM_RANDOM_OVERLAP 1
!==============================================================================
!
! !INPUT PARAMETERS:
!
    INTEGER,        INTENT(IN)    :: WLAOD       ! AOD calculated how?
                                                 ! (1: 550 nm, 0: 999 nm)
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input options
    TYPE(ChmState), INTENT(IN)    :: State_Chm   ! Chemistry State object
    TYPE(GrdState), INTENT(IN)    :: State_Grid  ! Grid State object
    TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology State object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(DgnState), INTENT(INOUT) :: State_Diag  ! Diagnostics State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!
! !REMARKS:
!
!  Parameter to choose cloud overlap algorithm:
!  ============================================================================
!  (1 ) OVERLAP (INTEGER) : 1 - Linear Approximation (used up to v7-04-12)
!                           2 - Approximate Random Overlap (default)
!                           3 - Maximum Random Overlap (computation intensive)
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER, SAVE :: LASTMONTH = -1
    INTEGER       :: NLON, NLAT, DAY,  MONTH, DAY_OF_YR, L, N
    INTEGER       :: IOPT, LCHEM
    REAL(fp)      :: CSZA, PRES, SFCA, YLAT,  O3_TOMS
    REAL(fp)      :: O3_CTM(State_Grid%NZ+1)
    REAL(fp)      :: T_CTM(State_Grid%NZ+1), OPTD(State_Grid%NZ)
    REAL(fp)      :: OPTDUST(State_Grid%NZ,NDUST)
    REAL(fp)      :: OPTAER(State_Grid%NZ,A_)

    ! Local variables for cloud overlap (hyl, phs)
    INTEGER       :: NUMB, KK, I
    INTEGER       :: INDIC(State_Grid%NZ+1)
    INTEGER       :: INDGEN(State_Grid%NZ+1)! = (/ (i,i=1,State_Grid%NZ+1) /)
    INTEGER       :: KBOT(State_Grid%NZ)
    INTEGER       :: KTOP(State_Grid%NZ)
    INTEGER       :: INDICATOR(State_Grid%NZ+2)
    REAL(fp)      :: FMAX(State_Grid%NZ)    ! maximum cloud fraction
                                              !  in a block, size can be to 
                                              !  FIX(State_Grid%NZ)+1
    REAL(fp)      :: CLDF1D(State_Grid%NZ)
    REAL(fp)      :: ODNEW(State_Grid%NZ)
    REAL(fp)      :: P_CTM(State_Grid%NZ+2)

    LOGICAL       :: AOD999
    LOGICAL, SAVE :: FIRST = .true.
    LOGICAL       :: prtDebug

    ! Species ID flags
    INTEGER, SAVE :: id_O3

    ! Strings
    CHARACTER(LEN=255) :: ErrMsg, ThisLoc

    !=================================================================
    ! FAST_JX begins here!
    !=================================================================

    ! Initialize
    RC        = GC_SUCCESS
    ErrMsg    = ''
    ThisLoc   = ' -> at Fast_JX (in module GeosCore/fast_jx_mod.F)'
    prtDebug  = ( Input_Opt%LPRT .and. Input_Opt%amIRoot)

    ! Get day of year (0-365 or 0-366)
    DAY_OF_YR = GET_DAY_OF_YEAR()

    ! Get current month
    MONTH     = GET_MONTH()

    ! Get day of month
    DAY       = GET_DAY()

    ! Was AOD calculated at 999 nm or reference?
    AOD999    = ( WLAOD == 0 )

    ! Zero diagnostic archival arrays to make sure that we don't have any
    ! leftover values from the last timestep near the top of the chemgrid
    IF ( State_Diag%Archive_UVFluxDiffuse ) THEN
       State_Diag%UVFluxDiffuse = 0.0_f4
    ENDIF
    IF ( State_Diag%Archive_UVFluxDirect ) THEN
       State_Diag%UVFluxDirect = 0.0_f4
    ENDIF
    IF ( State_Diag%Archive_UVFluxNet ) THEN
       State_Diag%UVFluxNet = 0.0_f4
    ENDIF

    !-----------------------------------------------------------------
    ! Special handling for first-time setup
    !-----------------------------------------------------------------
    IF ( FIRST ) THEN

       ! Get the species ID for O3
       id_O3 = Ind_('O3')
       IF ( id_O3 < 0 ) THEN
          ErrMsg = 'O3 is not a defined species!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

#ifdef USE_MAXIMUM_RANDOM_OVERLAP
       ! Special setup only for max random overlap
       DO i = 1,State_Grid%NZ+1
          INDGEN(i) = i       !(/(i,i=1,State_Grid%NZ+1)/)
       ENDDO
#endif

       ! Reset first-time flag
       FIRST = .FALSE.

    ENDIF

    !=================================================================
    ! For each (NLON,NLAT) location, call subroutine PHOTO_JX (in a
    ! parallel loop to compute J-values for the entire column.
    ! J-values will be stored in the common-block variable ZPJ, and
    ! will be later accessed via function FJXFUNC.
    !=================================================================
    !$OMP PARALLEL DO       &
    !$OMP DEFAULT( SHARED ) &
    !$OMP PRIVATE( NLAT,    NLON,   YLAT,      CSZA,    L       ) &
    !$OMP PRIVATE( P_CTM ,  T_CTM,  SFCA,      O3_TOMS, O3_CTM  ) &
    !$OMP PRIVATE( LCHEM,   OPTAER, N,         IOPT             ) &
    !$OMP PRIVATE( OPTDUST, OPTD,   CLDF1D                      ) &
#ifdef USE_MAXIMUM_RANDOM_OVERLAP
    !$OMP PRIVATE( FMAX,    KK,     NUMB,      KBOT             ) &
    !$OMP PRIVATE( KTOP     ODNEW,  INDICATOR, INDIC            ) &
#endif
    !$OMP SCHEDULE( DYNAMIC )

    ! Loop over latitudes and longitudes
    DO NLAT = 1, State_Grid%NY
    DO NLON = 1, State_Grid%NX

       ! Grid box latitude [degrees]
       YLAT = State_Grid%YMid( NLON, NLAT )

       ! Cosine of solar zenith angle [unitless] at (NLON,NLAT)
       CSZA = State_Met%SUNCOSmid(NLON,NLAT)

       ! Define the P array here
       DO L = 1, State_Grid%NZ+1
          P_CTM(L) = State_Met%PEDGE( NLON, NLAT, L )
       ENDDO

       ! Top edge of P_CTM is top of atmosphere (bmy, 2/13/07)
       P_CTM(State_Grid%NZ+2) = 0e+0_fp

       ! Temperature profile [K] at (NLON,NLAT)
       T_CTM(1:State_Grid%NZ) = State_Met%T( NLON, NLAT, 1:State_Grid%NZ)

       ! Top of atmosphere
       T_CTM(State_Grid%NZ+1) = T_CTM(State_Grid%NZ)

       ! Surface albedo [unitless] at (NLON,NLAT)
       SFCA = State_Met%UVALBEDO(NLON,NLAT)

       ! Overhead ozone column [DU] at (NLON, NLAT)
       ! These values are either from the met fields or TOMS/SBUV,
       ! depending on the settings in input.geos
       O3_TOMS = GET_OVERHEAD_O3( NLON, NLAT )

       ! CTM ozone densities (molec/cm3) at (NLON, NLAT)
       O3_CTM = 0e+0_fp
       LCHEM  = State_Met%ChemGridLev(NLON,NLAT)
       DO L = 1, LCHEM
          O3_CTM(L) = State_Chm%Species(NLON,NLAT,L,id_O3)
       ENDDO

       ! Aerosol OD profile [unitless] at (NLON,NLAT)
       ! and at 1000nm, IWV1000 (DAR)
       !OPTAER wants NAER*NRH values but ODAER is now NAER
       !use IRHARR to map to correct OPTAER bin (DAR 08/13)
       OPTAER = 0.0e+0_fp
       DO N = 1, NAER
       DO L = 1, State_Grid%NZ
          IOPT = ( (N-1) * NRH ) + IRHARR(NLON,NLAT,L)
          OPTAER(L,IOPT) = ODAER(NLON,NLAT,L,IWV1000,N)
       ENDDO
       ENDDO
       DO N = 1, NDUST
       DO L = 1, State_Grid%NZ
          OPTDUST(L,N) = ODMDUST(NLON,NLAT,L,IWV1000,N)
       ENDDO
       ENDDO

       ! Mineral dust OD profile [unitless] at (NLON,NLAT)
       ! and at 1000nm, IWV1000 (DAR)
       !OPTDUST = ODMDUST(NLON,NLAT,:,IWV1000,:)

       ! Cloud OD profile [unitless] at (NLON,NLAT)
       OPTD = State_Met%OPTD(NLON,NLAT,1:State_Grid%NZ)

       !-----------------------------------------------------------
       !### If you want to exclude aerosol OD, mineral dust OD,
       !### or cloud OD, then uncomment the following lines:
       !OPTAER  = 0d0
       !OPTDUST = 0d0
       !OPTD(:)    = 0d0
       !-----------------------------------------------------------

#if defined( USE_LINEAR_OVERLAP )
       !===========================================================
       ! %%%% CLOUD OVERLAP: LINEAR ASSUMPTION %%%%
       !
       ! Directly use OPTDEPTH = TAUCLD * CLDTOT
       !===========================================================

       ! Column cloud fraction (not less than zero)
       CLDF1D = State_Met%CLDF(NLON,NLAT,1:State_Grid%NZ)
       WHERE ( CLDF1D < 0e+0_fp ) CLDF1D = 0e+0_fp

       ! NOTE: For GEOS-FP and MERRA-2 met fields, the optical
       ! depth is the in-cloud optical depth.  At this point it has
       ! NOT been multiplied by cloud fraction yet.  Therefore, we can
       ! just apply the ! we can just apply the linear overlap formula
       ! as written above (i.e. multiply by cloud fraction).
       OPTD = OPTD * CLDF1D

       ! Call FAST-JX routines to compute J-values
       CALL PHOTO_JX( CSZA,       SFCA,      P_CTM,  T_CTM,  &
                      O3_CTM,     O3_TOMS,   AOD999, OPTAER, &
                      OPTDUST,    OPTD,      NLON,   NLAT,   &
                      YLAT,       DAY_OF_YR, MONTH,  DAY,    &
                      Input_Opt, State_Diag, State_Grid, State_Met )

#elif defined( USE_APPROX_RANDOM_OVERLAP )
       !===========================================================
       ! %%%% CLOUD OVERLAP: APPROX RANDOM OVERLAP ASSUMPTION %%%%
       !
       ! Use OPTDEPTH = TAUCLD * CLDTOT**1.5
       !===========================================================

       ! Column cloud fraction (not less than zero)
       CLDF1D = State_Met%CLDF(NLON,NLAT,1:State_Grid%NZ)
       WHERE ( CLDF1D < 0e+0_fp ) CLDF1D = 0e+0_fp

       ! NOTE: For GEOS-FP and MERRA-2 met fields, the optical
       ! depth is the in-cloud optical depth.  At this point it has
       ! NOT been multiplied by cloud fraction yet.  Therefore, we can
       ! just apply the approximate random overlap formula as written
       ! above (i.e. multiply by cloud fraction^1.5).
       OPTD = OPTD * ( CLDF1D )**1.5e+0_fp

       ! Call FAST-JX routines to compute J-values
       CALL PHOTO_JX( CSZA,       SFCA,      P_CTM,  T_CTM,  &
                      O3_CTM,     O3_TOMS,   AOD999, OPTAER, &
                      OPTDUST,    OPTD,      NLON,   NLAT,   &
                      YLAT,       DAY_OF_YR, MONTH,  DAY,    &
                      Input_Opt, State_Diag, State_Grid, State_Met )

#elif defined( USE_MAXIMUM_RANDOM_OVERLAP )
       !===========================================================
       ! %%%% CLOUD OVERLAP: MAXIMUM RANDOM OVERLAP %%%%
       !
       ! The Maximum-Random Overlap (MRAN) scheme assumes that
       ! clouds in adjacent layers are maximally overlapped to
       ! form a cloud block and that blocks of clouds separated by
       ! clear layers are randomly overlapped.  A vertical profile
       ! of fractional cloudiness is converted into a series of
       ! column configurations with corresponding fractions
       ! (see Liu et al., JGR 2006; hyl,3/3/04).
       !
       ! For more details about cloud overlap assumptions and
       ! their effect on photolysis frequencies and key oxidants
       ! in the troposphere, refer to the following articles:
       !
       ! (1) Liu, H., et al., Radiative effect of clouds on
       !      tropospheric chemistry in a global three-dimensional
       !      chemical transport model, J. Geophys. Res., vol.111,
       !      D20303, doi:10.1029/2005JD006403, 2006.
       ! (2) Tie, X., et al., Effect of clouds on photolysis and
       !      oxidants in the troposphere, J. Geophys. Res.,
       !      108(D20), 4642, doi:10.1029/2003JD003659, 2003.
       ! (3) Feng, Y., et al., Effects of cloud overlap in
       !      photochemical models, J. Geophys. Res., 109,
       !      D04310, doi:10.1029/2003JD004040, 2004.
       ! (4) Stubenrauch, C.J., et al., Implementation of subgrid
       !      cloud vertical structure inside a GCM and its effect
       !      on the radiation budget, J. Clim., 10, 273-287, 1997.
       !-----------------------------------------------------------
       ! MMRAN needs IN-CLOUD optical depth (ODNEW) as input
       ! Use cloud fraction, instead of OPTD, to form cloud blocks
       ! (hyl,06/19/04)
       !===========================================================

       ! Sort this out later
       CALL ERROR_STOP('MMRAN_16 not yet FJX compatible.', 'fast_jx_mod.F90')

       !! Initialize
       !FMAX(:)   = 0d0  ! max cloud fraction in each cloud block
       !ODNEW(:)  = 0d0  ! in-cloud optical depth
       !CLDF1D    = State_Met%CLDF(1:State_Grid%NZ,NLON,NLAT)
       !INDICATOR = 0
       !
       !! set small negative CLDF or OPTD to zero.
       !! Set indicator vector.
       !WHERE ( CLDF1D <= 0d0 )
       !   CLDF1D               = 0d0
       !   OPTD                 = 0D0
       !ELSEWHERE
       !   INDICATOR(2:State_Grid%NZ+1) = 1
       !ENDWHERE
       !
       !! Prevent negative opt depth
       !WHERE ( OPTD < 0D0 ) OPTD   = 0D0
       !
       !!--------------------------------------------------------
       !! Generate cloud blocks & get their Bottom and Top levels
       !!--------------------------------------------------------
       !INDICATOR = CSHIFT(INDICATOR, 1) - INDICATOR
       !INDIC     = INDICATOR(1:State_Grid%NZ+1)
       !
       !! Number of cloud block
       !NUMB      = COUNT( INDIC == 1 )
       !
       !! Bottom layer of each block
       !KBOT(1:NUMB) = PACK(INDGEN, (INDIC == 1 ) )
       !
       !! Top layer of each block
       !KTOP(1:NUMB) = PACK(INDGEN, (INDIC == -1) ) - 1
       !
       !!--------------------------------------------------------
       !! For each cloud block, get Max Cloud Fractions, and
       !! in-cloud optical depth vertical distribution.
       !!--------------------------------------------------------
       !DO KK = 1, NUMB
       !
       !   ! Max cloud fraction
       !   FMAX(KK) = MAXVAL( CLDF1D(KBOT(KK):KTOP(KK)) )
       !
       !   ! NOTE: for the GEOS-FP and MERRA-2 met fields (i.e. with
       !   ! optical depth & cloud fractions regridded with RegridTau)
       !   ! OPTD is the in-cloud optical depth.  At this point it has
       !   ! NOT been multiplied by cloud fraction yet.  Therefore,
       !   ! we can just set ODNEW = OPTD. (bmy, hyl, 10/24/08)
       !
       !   ! ODNEW is adjusted in-cloud OD vertical distrib.
       !   ODNEW(KBOT(KK):KTOP(KK)) = OPTD(KBOT(KK):KTOP(KK))
       !
       !ENDDO
       !
       !!--------------------------------------------------------
       !! Apply Max RANdom if 1-6 clouds blocks, else use linear
       !!--------------------------------------------------------
       !SELECT CASE( NUMB )
       !
       !CASE( 0,7: )
       !   CALL PHOTOJ( NLON,  NLAT,     YLAT,    DAY_OF_YR,
       !                MONTH, DAY,      CSZA,    TEMP,
       !                SFCA,  OPTD,     OPTDUST, OPTAER,
       !                O3COL )
       !
       !CASE( 1:6 )
       !    CALL MMRAN_16( NUMB,  NLON,  NLAT,      YLAT,
       !                   DAY,   MONTH, DAY_OF_YR, CSZA,
       !                   TEMP,  SFCA,  OPTDUST,   OPTAER,
       !                   State_Grid%NZ, FMAX,  ODNEW,     KBOT,
       !                   KTOP,  O3COL )
       !
       !END SELECT
#endif

    ENDDO
    ENDDO
    !$OMP END PARALLEL DO

    ! Reset first-time flag
    FIRST=.FALSE.

  END SUBROUTINE FAST_JX
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: blkslv
!
! !DESCRIPTION: Subroutine BLKSLV solves the block tri-diagonal system
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE BLKSLV(FJ,POMEGA,FZ,ZTAU,ZFLUX,RFL,PM,PM0,FJTOP,FJBOT,ND)
!
! !USES:
!
!
! !INPUT PARAMETERS:
!
    INTEGER,  INTENT(IN)  :: ND
    REAL(fp), INTENT(IN)  :: POMEGA(M2_,N_,W_)
    REAL(fp), INTENT(IN)  :: FZ(N_,W_)
    REAL(fp), INTENT(IN)  :: ZTAU(N_,W_)
    REAL(fp), INTENT(IN)  :: PM(M_,M2_)
    REAL(fp), INTENT(IN)  :: PM0(M2_)
    REAL(fp), INTENT(IN)  :: RFL(W_)
    REAL(fp), INTENT(IN)  :: ZFLUX(W_)
!
! !OUTPUT PARAMETERS:
!
    REAL(fp), INTENT(OUT) :: FJ(N_,W_)
    REAL(fp), INTENT(OUT) :: FJTOP(W_)
    REAL(fp), INTENT(OUT) :: FJBOT(W_)
!
! !REMARKS:
! The block tri-diagonal system:
!       A(I)*X(I-1) + B(I)*X(I) + C(I)*X(I+1) = H(I)
!
! !REVISION HISTORY:
!  27 Mar 2013 - S. D. Eastham - Copied from GEOS-Chem v9-01-03
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    REAL(fp), DIMENSION(M_,N_,W_)    ::  A,C,H,   RR
    REAL(fp), DIMENSION(M_,M_,N_,W_) ::  B,AA,CC,  DD
    REAL(fp), DIMENSION(M_,M_)       ::  E
    REAL(fp)  SUMB,SUMBX,SUMT
    INTEGER I, J, K, L

    !=================================================================
    ! BLKSLV begins here!
    !=================================================================

    do K = 1,W_
       call GEN_ID (POMEGA(1,1,K),FZ(1,K),ZTAU(1,K),ZFLUX(K),RFL(K), &
                    PM,PM0, B(1,1,1,K),CC(1,1,1,K),AA(1,1,1,K),      &
                    A(1,1,K),H(1,1,K),C(1,1,K), ND)
    enddo

    do K = 1,W_
       ! UPPER BOUNDARY L=1
       L = 1
       do J = 1,M_
       do I = 1,M_
          E(I,J) = B(I,J,1,K)
       enddo
       enddo

       ! setup L & U matrices
       E(2,1) = E(2,1)/E(1,1)
       E(2,2) = E(2,2)-E(2,1)*E(1,2)
       E(2,3) = E(2,3)-E(2,1)*E(1,3)
       E(2,4) = E(2,4)-E(2,1)*E(1,4)
       E(3,1) = E(3,1)/E(1,1)
       E(3,2) = (E(3,2)-E(3,1)*E(1,2))/E(2,2)
       E(3,3) = E(3,3)-E(3,1)*E(1,3)-E(3,2)*E(2,3)
       E(3,4) = E(3,4)-E(3,1)*E(1,4)-E(3,2)*E(2,4)
       E(4,1) = E(4,1)/E(1,1)
       E(4,2) = (E(4,2)-E(4,1)*E(1,2))/E(2,2)
       E(4,3) = (E(4,3)-E(4,1)*E(1,3)-E(4,2)*E(2,3))/E(3,3)
       E(4,4) = E(4,4)-E(4,1)*E(1,4)-E(4,2)*E(2,4)-E(4,3)*E(3,4)
       ! invert L
       E(4,3) = -E(4,3)
       E(4,2) = -E(4,2)-E(4,3)*E(3,2)
       E(4,1) = -E(4,1)-E(4,2)*E(2,1)-E(4,3)*E(3,1)
       E(3,2) = -E(3,2)
       E(3,1) = -E(3,1)-E(3,2)*E(2,1)
       E(2,1) = -E(2,1)
       ! invert U
       E(4,4) = 1.e+0_fp/E(4,4)
       E(3,4) = -E(3,4)*E(4,4)/E(3,3)
       E(3,3) = 1.e+0_fp/E(3,3)
       E(2,4) = -(E(2,3)*E(3,4)+E(2,4)*E(4,4))/E(2,2)
       E(2,3) = -E(2,3)*E(3,3)/E(2,2)
       E(2,2) = 1.e+0_fp/E(2,2)
       E(1,4) = -(E(1,2)*E(2,4)+E(1,3)*E(3,4)+E(1,4)*E(4,4))/E(1,1)
       E(1,3) = -(E(1,2)*E(2,3)+E(1,3)*E(3,3))/E(1,1)
       E(1,2) = -E(1,2)*E(2,2)/E(1,1)
       E(1,1) = 1.e+0_fp/E(1,1)
       ! multiply U-invers * L-inverse
       E(1,1) = E(1,1)+E(1,2)*E(2,1)+E(1,3)*E(3,1)+E(1,4)*E(4,1)
       E(1,2) = E(1,2)+E(1,3)*E(3,2)+E(1,4)*E(4,2)
       E(1,3) = E(1,3)+E(1,4)*E(4,3)
       E(2,1) = E(2,2)*E(2,1)+E(2,3)*E(3,1)+E(2,4)*E(4,1)
       E(2,2) = E(2,2)+E(2,3)*E(3,2)+E(2,4)*E(4,2)
       E(2,3) = E(2,3)+E(2,4)*E(4,3)
       E(3,1) = E(3,3)*E(3,1)+E(3,4)*E(4,1)
       E(3,2) = E(3,3)*E(3,2)+E(3,4)*E(4,2)
       E(3,3) = E(3,3)+E(3,4)*E(4,3)
       E(4,1) = E(4,4)*E(4,1)
       E(4,2) = E(4,4)*E(4,2)
       E(4,3) = E(4,4)*E(4,3)

       do J = 1,M_
          do I = 1,M_
             DD(I,J,1,K) = -E(I,1)*CC(1,J,1,K)-E(I,2)*CC(2,J,1,K) &
                           -E(I,3)*CC(3,J,1,K)-E(I,4)*CC(4,J,1,K)
          enddo
          RR(J,1,K) = E(J,1)*H(1,1,K)+E(J,2)*H(2,1,K) &
                    + E(J,3)*H(3,1,K)+E(J,4)*H(4,1,K)
       enddo

       ! CONTINUE THROUGH ALL DEPTH POINTS ID=2 TO ID=ND-1
       do L = 2,ND-1

          do J = 1,M_
             do I = 1,M_
                B(I,J,L,K) = B(I,J,L,K) + A(I,L,K)*DD(I,J,L-1,K)
             enddo
             H(J,L,K) = H(J,L,K) - A(J,L,K)*RR(J,L-1,K)
          enddo

          do J = 1,M_
             do I = 1,M_
                E(I,J) = B(I,J,L,K)
             enddo
          enddo

          ! setup L & U matrices
          E(2,1) = E(2,1)/E(1,1)
          E(2,2) = E(2,2)-E(2,1)*E(1,2)
          E(2,3) = E(2,3)-E(2,1)*E(1,3)
          E(2,4) = E(2,4)-E(2,1)*E(1,4)
          E(3,1) = E(3,1)/E(1,1)
          E(3,2) = (E(3,2)-E(3,1)*E(1,2))/E(2,2)
          E(3,3) = E(3,3)-E(3,1)*E(1,3)-E(3,2)*E(2,3)
          E(3,4) = E(3,4)-E(3,1)*E(1,4)-E(3,2)*E(2,4)
          E(4,1) = E(4,1)/E(1,1)
          E(4,2) = (E(4,2)-E(4,1)*E(1,2))/E(2,2)
          E(4,3) = (E(4,3)-E(4,1)*E(1,3)-E(4,2)*E(2,3))/E(3,3)
          E(4,4) = E(4,4)-E(4,1)*E(1,4)-E(4,2)*E(2,4)-E(4,3)*E(3,4)
          ! invert L
          E(4,3) = -E(4,3)
          E(4,2) = -E(4,2)-E(4,3)*E(3,2)
          E(4,1) = -E(4,1)-E(4,2)*E(2,1)-E(4,3)*E(3,1)
          E(3,2) = -E(3,2)
          E(3,1) = -E(3,1)-E(3,2)*E(2,1)
          E(2,1) = -E(2,1)
          ! invert U
          E(4,4) = 1.e+0_fp/E(4,4)
          E(3,4) = -E(3,4)*E(4,4)/E(3,3)
          E(3,3) = 1.e+0_fp/E(3,3)
          E(2,4) = -(E(2,3)*E(3,4)+E(2,4)*E(4,4))/E(2,2)
          E(2,3) = -E(2,3)*E(3,3)/E(2,2)
          E(2,2) = 1.e+0_fp/E(2,2)
          E(1,4) = -(E(1,2)*E(2,4)+E(1,3)*E(3,4)+E(1,4)*E(4,4))/E(1,1)
          E(1,3) = -(E(1,2)*E(2,3)+E(1,3)*E(3,3))/E(1,1)
          E(1,2) = -E(1,2)*E(2,2)/E(1,1)
          E(1,1) = 1.e+0_fp/E(1,1)
          ! multiply U-invers * L-inverse
          E(1,1) = E(1,1)+E(1,2)*E(2,1)+E(1,3)*E(3,1)+E(1,4)*E(4,1)
          E(1,2) = E(1,2)+E(1,3)*E(3,2)+E(1,4)*E(4,2)
          E(1,3) = E(1,3)+E(1,4)*E(4,3)
          E(2,1) = E(2,2)*E(2,1)+E(2,3)*E(3,1)+E(2,4)*E(4,1)
          E(2,2) = E(2,2)+E(2,3)*E(3,2)+E(2,4)*E(4,2)
          E(2,3) = E(2,3)+E(2,4)*E(4,3)
          E(3,1) = E(3,3)*E(3,1)+E(3,4)*E(4,1)
          E(3,2) = E(3,3)*E(3,2)+E(3,4)*E(4,2)
          E(3,3) = E(3,3)+E(3,4)*E(4,3)
          E(4,1) = E(4,4)*E(4,1)
          E(4,2) = E(4,4)*E(4,2)
          E(4,3) = E(4,4)*E(4,3)

          do J = 1,M_
             do I = 1,M_
                DD(I,J,L,K) = - E(I,J)*C(J,L,K)
             enddo
             RR(J,L,K) = E(J,1)*H(1,L,K)+E(J,2)*H(2,L,K) &
                       + E(J,3)*H(3,L,K)+E(J,4)*H(4,L,K)
          enddo

       enddo

       ! FINAL DEPTH POINT: L=ND
       L = ND
       do J = 1,M_
          do I = 1,M_
             B(I,J,L,K) = B(I,J,L,K) &
                   + AA(I,1,L,K)*DD(1,J,L-1,K) + AA(I,2,L,K)*DD(2,J,L-1,K) &
                   + AA(I,3,L,K)*DD(3,J,L-1,K) + AA(I,4,L,K)*DD(4,J,L-1,K)
          enddo
          H(J,L,K) = H(J,L,K) &
                   - AA(J,1,L,K)*RR(1,L-1,K) - AA(J,2,L,K)*RR(2,L-1,K) &
                   - AA(J,3,L,K)*RR(3,L-1,K) - AA(J,4,L,K)*RR(4,L-1,K)
       enddo

       do J = 1,M_
          do I = 1,M_
             E(I,J) = B(I,J,L,K)
          enddo
       enddo

       ! setup L & U matrices
       E(2,1) = E(2,1)/E(1,1)
       E(2,2) = E(2,2)-E(2,1)*E(1,2)
       E(2,3) = E(2,3)-E(2,1)*E(1,3)
       E(2,4) = E(2,4)-E(2,1)*E(1,4)
       E(3,1) = E(3,1)/E(1,1)
       E(3,2) = (E(3,2)-E(3,1)*E(1,2))/E(2,2)
       E(3,3) = E(3,3)-E(3,1)*E(1,3)-E(3,2)*E(2,3)
       E(3,4) = E(3,4)-E(3,1)*E(1,4)-E(3,2)*E(2,4)
       E(4,1) = E(4,1)/E(1,1)
       E(4,2) = (E(4,2)-E(4,1)*E(1,2))/E(2,2)
       E(4,3) = (E(4,3)-E(4,1)*E(1,3)-E(4,2)*E(2,3))/E(3,3)
       E(4,4) = E(4,4)-E(4,1)*E(1,4)-E(4,2)*E(2,4)-E(4,3)*E(3,4)
       ! invert L
       E(4,3) = -E(4,3)
       E(4,2) = -E(4,2)-E(4,3)*E(3,2)
       E(4,1) = -E(4,1)-E(4,2)*E(2,1)-E(4,3)*E(3,1)
       E(3,2) = -E(3,2)
       E(3,1) = -E(3,1)-E(3,2)*E(2,1)
       E(2,1) = -E(2,1)
       ! invert U
       E(4,4) = 1.e+0_fp/E(4,4)
       E(3,4) = -E(3,4)*E(4,4)/E(3,3)
       E(3,3) = 1.e+0_fp/E(3,3)
       E(2,4) = -(E(2,3)*E(3,4)+E(2,4)*E(4,4))/E(2,2)
       E(2,3) = -E(2,3)*E(3,3)/E(2,2)
       E(2,2) = 1.e+0_fp/E(2,2)
       E(1,4) = -(E(1,2)*E(2,4)+E(1,3)*E(3,4)+E(1,4)*E(4,4))/E(1,1)
       E(1,3) = -(E(1,2)*E(2,3)+E(1,3)*E(3,3))/E(1,1)
       E(1,2) = -E(1,2)*E(2,2)/E(1,1)
       E(1,1) = 1.e+0_fp/E(1,1)
       ! multiply U-invers * L-inverse
       E(1,1) = E(1,1)+E(1,2)*E(2,1)+E(1,3)*E(3,1)+E(1,4)*E(4,1)
       E(1,2) = E(1,2)+E(1,3)*E(3,2)+E(1,4)*E(4,2)
       E(1,3) = E(1,3)+E(1,4)*E(4,3)
       E(2,1) = E(2,2)*E(2,1)+E(2,3)*E(3,1)+E(2,4)*E(4,1)
       E(2,2) = E(2,2)+E(2,3)*E(3,2)+E(2,4)*E(4,2)
       E(2,3) = E(2,3)+E(2,4)*E(4,3)
       E(3,1) = E(3,3)*E(3,1)+E(3,4)*E(4,1)
       E(3,2) = E(3,3)*E(3,2)+E(3,4)*E(4,2)
       E(3,3) = E(3,3)+E(3,4)*E(4,3)
       E(4,1) = E(4,4)*E(4,1)
       E(4,2) = E(4,4)*E(4,2)
       E(4,3) = E(4,4)*E(4,3)

       do J = 1,M_
          RR(J,L,K) = E(J,1)*H(1,L,K)+E(J,2)*H(2,L,K) &
                    + E(J,3)*H(3,L,K)+E(J,4)*H(4,L,K)
       enddo

       ! BACK SOLUTION
       do L = ND-1,1,-1
          do J = 1,M_
             RR(J,L,K) = RR(J,L,K) &
                       + DD(J,1,L,K)*RR(1,L+1,K) + DD(J,2,L,K)*RR(2,L+1,K) &
                       + DD(J,3,L,K)*RR(3,L+1,K) + DD(J,4,L,K)*RR(4,L+1,K)
          enddo
       enddo

       ! mean J & H
       do L = 1,ND,2
          FJ(L,K) = RR(1,L,K)*WT(1) + RR(2,L,K)*WT(2) &
                  + RR(3,L,K)*WT(3) + RR(4,L,K)*WT(4)
       enddo
       do L = 2,ND,2
          FJ(L,K) = RR(1,L,K)*WT(1)*EMU(1) + RR(2,L,K)*WT(2)*EMU(2) &
                  + RR(3,L,K)*WT(3)*EMU(3) + RR(4,L,K)*WT(4)*EMU(4)
       enddo

       ! FJTOP = scaled diffuse flux out top-of-atmosphere (limit = mu0)
       ! FJBOT = scaled diffuse flux onto surface:
       ! ZFLUX = reflect/(1 + reflect) * mu0 * Fsolar(lower boundary)
       ! SUMBX = flux from Lambert reflected I+
       SUMT = RR(1, 1,K)*WT(1)*EMU(1) + RR(2, 1,K)*WT(2)*EMU(2) &
            + RR(3, 1,K)*WT(3)*EMU(3) + RR(4, 1,K)*WT(4)*EMU(4)
       SUMB = RR(1,ND,K)*WT(1)*EMU(1) + RR(2,ND,K)*WT(2)*EMU(2) &
            + RR(3,ND,K)*WT(3)*EMU(3) + RR(4,ND,K)*WT(4)*EMU(4)
       SUMBX = 4.e+0_fp*SUMB*RFL(K)/(1.0e+0_fp + RFL(K)) + ZFLUX(K)

       FJTOP(K) = 4.e+0_fp*SUMT
       FJBOT(K) = 4.e+0_fp*SUMB - SUMBX

    enddo

  END SUBROUTINE BLKSLV
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: gen_id
!
! !DESCRIPTION: Subroutine GEN generates coefficient matrices for the block
!  tri-diagonal system described in BLKSLV.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE GEN_ID(POMEGA,FZ,ZTAU,ZFLUX,RFL,PM,PM0,B,CC,AA,A,H,C,ND)
!
! !INPUT PARAMETERS:
!
    INTEGER,  INTENT(IN)  :: ND
    REAL(fp), INTENT(IN)  :: POMEGA(M2_,N_)
    REAL(fp), INTENT(IN)  :: PM(M_,M2_)
    REAL(fp), INTENT(IN)  :: PM0(M2_)
    REAL(fp), INTENT(IN)  :: ZFLUX,RFL
    REAL(fp), INTENT(IN), DIMENSION(N_) :: FZ,ZTAU
!
! !OUTPUT PARAMETERS:
!
    REAL(fp), INTENT(OUT),DIMENSION(M_,M_,N_) :: B,AA,CC
    REAL(fp), INTENT(OUT),DIMENSION(M_,N_)    :: A,C,H

!
! !REVISION HISTORY:
!  28 Mar 2013 - S. D. Eastham - Copied from Fast-JX v7.0
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER I, J, K, L1,L2,LL
    REAL(fp)  SUM0, SUM1, SUM2, SUM3
    REAL(fp)  DELTAU, D1, D2, SURFAC

    REAL(fp), DIMENSION(M_,M_) :: S,T,U,V,W

    !=================================================================
    ! GEN_ID begins here!
    !=================================================================

    ! upper boundary:  2nd-order terms
    L1 = 1
    L2 = 2
    do I = 1,M_
       SUM0 = POMEGA(1,L1)*PM(I,1)*PM0(1) + POMEGA(3,L1)*PM(I,3)*PM0(3) &
            + POMEGA(5,L1)*PM(I,5)*PM0(5) + POMEGA(7,L1)*PM(I,7)*PM0(7)
       SUM2 = POMEGA(1,L2)*PM(I,1)*PM0(1) + POMEGA(3,L2)*PM(I,3)*PM0(3) &
            + POMEGA(5,L2)*PM(I,5)*PM0(5) + POMEGA(7,L2)*PM(I,7)*PM0(7)
       SUM1 = POMEGA(2,L1)*PM(I,2)*PM0(2) + POMEGA(4,L1)*PM(I,4)*PM0(4) &
            + POMEGA(6,L1)*PM(I,6)*PM0(6) + POMEGA(8,L1)*PM(I,8)*PM0(8)
       SUM3 = POMEGA(2,L2)*PM(I,2)*PM0(2) + POMEGA(4,L2)*PM(I,4)*PM0(4) &
            + POMEGA(6,L2)*PM(I,6)*PM0(6) + POMEGA(8,L2)*PM(I,8)*PM0(8)
       H(I,L1) = 0.5e+0_fp*(SUM0*FZ(L1) + SUM2*FZ(L2))
       A(I,L1) = 0.5e+0_fp*(SUM1*FZ(L1) + SUM3*FZ(L2))
    enddo

    do I = 1,M_
    do J = 1,I
       SUM0 = POMEGA(1,L1)*PM(I,1)*PM(J,1) + POMEGA(3,L1)*PM(I,3)*PM(J,3) &
            + POMEGA(5,L1)*PM(I,5)*PM(J,5) + POMEGA(7,L1)*PM(I,7)*PM(J,7)
       SUM2 = POMEGA(1,L2)*PM(I,1)*PM(J,1) + POMEGA(3,L2)*PM(I,3)*PM(J,3) &
            + POMEGA(5,L2)*PM(I,5)*PM(J,5) + POMEGA(7,L2)*PM(I,7)*PM(J,7)
       SUM1 = POMEGA(2,L1)*PM(I,2)*PM(J,2) + POMEGA(4,L1)*PM(I,4)*PM(J,4) &
            + POMEGA(6,L1)*PM(I,6)*PM(J,6) + POMEGA(8,L1)*PM(I,8)*PM(J,8)
       SUM3 = POMEGA(2,L2)*PM(I,2)*PM(J,2) + POMEGA(4,L2)*PM(I,4)*PM(J,4) &
            + POMEGA(6,L2)*PM(I,6)*PM(J,6) + POMEGA(8,L2)*PM(I,8)*PM(J,8)
       S(I,J) = - SUM2*WT(J)
       S(J,I) = - SUM2*WT(I)
       T(I,J) = - SUM1*WT(J)
       T(J,I) = - SUM1*WT(I)
       V(I,J) = - SUM3*WT(J)
       V(J,I) = - SUM3*WT(I)
       B(I,J,L1) = - 0.5e+0_fp*(SUM0 + SUM2)*WT(J)
       B(J,I,L1) = - 0.5e+0_fp*(SUM0 + SUM2)*WT(I)
    enddo
    enddo

    do I = 1,M_
       S(I,I)   = S(I,I)    + 1.0e+0_fp
       T(I,I)   = T(I,I)    + 1.0e+0_fp
       V(I,I)   = V(I,I)    + 1.0e+0_fp
       B(I,I,L1)= B(I,I,L1) + 1.0e+0_fp

       C(I,L1)= S(I,1)*A(1,L1)/EMU(1) + S(I,2)*A(2,L1)/EMU(2) &
              + S(I,3)*A(3,L1)/EMU(3) + S(I,4)*A(4,L1)/EMU(4)
    enddo

    do I = 1,M_
    do J = 1,M_
       W(J,I) = S(J,1)*T(1,I)/EMU(1) + S(J,2)*T(2,I)/EMU(2) &
              + S(J,3)*T(3,I)/EMU(3) + S(J,4)*T(4,I)/EMU(4)
       U(J,I) = S(J,1)*V(1,I)/EMU(1) + S(J,2)*V(2,I)/EMU(2) &
              + S(J,3)*V(3,I)/EMU(3) + S(J,4)*V(4,I)/EMU(4)
    enddo
    enddo
    ! upper boundary, 2nd-order, C-matrix is full (CC)
    DELTAU = ZTAU(L2) - ZTAU(L1)
    D2 = 0.25e+0_fp*DELTAU
    do I = 1,M_
       do J = 1,M_
          B(I,J,L1) = B(I,J,L1) + D2*W(I,J)
          CC(I,J,L1) = D2*U(I,J)
       enddo
       H(I,L1) = H(I,L1) + 2.0e+0_fp*D2*C(I,L1)
       A(I,L1) = 0.0e+0_fp
    enddo
    do I = 1,M_
       D1 = EMU(I)/DELTAU
       B(I,I,L1)  = B(I,I,L1) + D1
       CC(I,I,L1) = CC(I,I,L1) - D1
    enddo

    ! intermediate points:  can be even or odd, A & C diagonal
    ! mid-layer h-points, Legendre terms 2,4,6,8
    do LL=2,ND-1,2
       DELTAU = ZTAU(LL+1) - ZTAU(LL-1)
       do I = 1,M_
          A(I,LL) = EMU(I)/DELTAU
          C(I,LL) = -A(I,LL)
          H(I,LL) = FZ(LL)*( &
               POMEGA(2,LL)*PM(I,2)*PM0(2) + POMEGA(4,LL)*PM(I,4)*PM0(4) &
             + POMEGA(6,LL)*PM(I,6)*PM0(6) + POMEGA(8,LL)*PM(I,8)*PM0(8))
       enddo
       do I = 1,M_
       do J=1,I
          SUM0 = POMEGA(2,LL)*PM(I,2)*PM(J,2) + POMEGA(4,LL)*PM(I,4)*PM(J,4) &
               + POMEGA(6,LL)*PM(I,6)*PM(J,6) + POMEGA(8,LL)*PM(I,8)*PM(J,8)
          B(I,J,LL) =  - SUM0*WT(J)
          B(J,I,LL) =  - SUM0*WT(I)
       enddo
       enddo
       do I = 1,M_
          B(I,I,LL) = B(I,I,LL) + 1.0e+0_fp
       enddo
    enddo

    ! odd-layer j-points, Legendre terms 1,3,5,7
    do LL=3,ND-2,2
       DELTAU = ZTAU(LL+1) - ZTAU(LL-1)
       do I = 1,M_
          A(I,LL) = EMU(I)/DELTAU
          C(I,LL) = -A(I,LL)
          H(I,LL) = FZ(LL)*( &
               POMEGA(1,LL)*PM(I,1)*PM0(1) + POMEGA(3,LL)*PM(I,3)*PM0(3) &
             + POMEGA(5,LL)*PM(I,5)*PM0(5) + POMEGA(7,LL)*PM(I,7)*PM0(7))
       enddo
       do I = 1,M_
       do J=1,I
          SUM0 = POMEGA(1,LL)*PM(I,1)*PM(J,1) + POMEGA(3,LL)*PM(I,3)*PM(J,3) &
               + POMEGA(5,LL)*PM(I,5)*PM(J,5) + POMEGA(7,LL)*PM(I,7)*PM(J,7)
          B(I,J,LL) =  - SUM0*WT(J)
          B(J,I,LL) =  - SUM0*WT(I)
       enddo
       enddo
       do I = 1,M_
          B(I,I,LL) = B(I,I,LL) + 1.0e+0_fp
       enddo
    enddo

    ! lower boundary:  2nd-order terms
    L1 = ND
    L2 = ND-1
    do I = 1,M_
       SUM0 = POMEGA(1,L1)*PM(I,1)*PM0(1) + POMEGA(3,L1)*PM(I,3)*PM0(3) &
            + POMEGA(5,L1)*PM(I,5)*PM0(5) + POMEGA(7,L1)*PM(I,7)*PM0(7)
       SUM2 = POMEGA(1,L2)*PM(I,1)*PM0(1) + POMEGA(3,L2)*PM(I,3)*PM0(3) &
            + POMEGA(5,L2)*PM(I,5)*PM0(5) + POMEGA(7,L2)*PM(I,7)*PM0(7)
       SUM1 = POMEGA(2,L1)*PM(I,2)*PM0(2) + POMEGA(4,L1)*PM(I,4)*PM0(4) &
            + POMEGA(6,L1)*PM(I,6)*PM0(6) + POMEGA(8,L1)*PM(I,8)*PM0(8)
       SUM3 = POMEGA(2,L2)*PM(I,2)*PM0(2) + POMEGA(4,L2)*PM(I,4)*PM0(4) &
            + POMEGA(6,L2)*PM(I,6)*PM0(6) + POMEGA(8,L2)*PM(I,8)*PM0(8)
       H(I,L1) = 0.5e+0_fp*(SUM0*FZ(L1) + SUM2*FZ(L2))
       A(I,L1) = 0.5e+0_fp*(SUM1*FZ(L1) + SUM3*FZ(L2))
    enddo

    do I = 1,M_
    do J = 1,I
       SUM0 = POMEGA(1,L1)*PM(I,1)*PM(J,1) + POMEGA(3,L1)*PM(I,3)*PM(J,3) &
            + POMEGA(5,L1)*PM(I,5)*PM(J,5) + POMEGA(7,L1)*PM(I,7)*PM(J,7)
       SUM2 = POMEGA(1,L2)*PM(I,1)*PM(J,1) + POMEGA(3,L2)*PM(I,3)*PM(J,3) &
            + POMEGA(5,L2)*PM(I,5)*PM(J,5) + POMEGA(7,L2)*PM(I,7)*PM(J,7)
       SUM1 = POMEGA(2,L1)*PM(I,2)*PM(J,2) + POMEGA(4,L1)*PM(I,4)*PM(J,4) &
            + POMEGA(6,L1)*PM(I,6)*PM(J,6) + POMEGA(8,L1)*PM(I,8)*PM(J,8)
       SUM3 = POMEGA(2,L2)*PM(I,2)*PM(J,2) + POMEGA(4,L2)*PM(I,4)*PM(J,4) &
            + POMEGA(6,L2)*PM(I,6)*PM(J,6) + POMEGA(8,L2)*PM(I,8)*PM(J,8)
       S(I,J) = - SUM2*WT(J)
       S(J,I) = - SUM2*WT(I)
       T(I,J) = - SUM1*WT(J)
       T(J,I) = - SUM1*WT(I)
       V(I,J) = - SUM3*WT(J)
       V(J,I) = - SUM3*WT(I)
       B(I,J,L1) = - 0.5e+0_fp*(SUM0 + SUM2)*WT(J)
       B(J,I,L1) = - 0.5e+0_fp*(SUM0 + SUM2)*WT(I)
    enddo
    enddo

    do I = 1,M_
       S(I,I)   = S(I,I)   + 1.0e+0_fp
       T(I,I)   = T(I,I)   + 1.0e+0_fp
       V(I,I)   = V(I,I)   + 1.0e+0_fp
       B(I,I,L1)= B(I,I,L1) + 1.0e+0_fp

       C(I,L1)= S(I,1)*A(1,L1)/EMU(1) + S(I,2)*A(2,L1)/EMU(2) &
              + S(I,3)*A(3,L1)/EMU(3) + S(I,4)*A(4,L1)/EMU(4)
    enddo

    do I = 1,M_
    do J = 1,M_
       W(J,I) = S(J,1)*T(1,I)/EMU(1) + S(J,2)*T(2,I)/EMU(2) &
              + S(J,3)*T(3,I)/EMU(3) + S(J,4)*T(4,I)/EMU(4)
       U(J,I) = S(J,1)*V(1,I)/EMU(1) + S(J,2)*V(2,I)/EMU(2) &
              + S(J,3)*V(3,I)/EMU(3) + S(J,4)*V(4,I)/EMU(4)
    enddo
    enddo

    ! lower boundary, 2nd-order, A-matrix is full (AA)
    DELTAU = ZTAU(L1) - ZTAU(L2)
    D2 = 0.25e+0_fp*DELTAU
    SURFAC = 4.0e+0_fp*RFL/(1.0e+0_fp + RFL)
    do I = 1,M_
       D1 = EMU(I)/DELTAU
       SUM0 = D1 + D2*(W(I,1)+W(I,2)+W(I,3)+W(I,4))
       SUM1 = SURFAC*SUM0
       do J = 1,M_
          AA(I,J,L1) = - D2*U(I,J)
          B(I,J,L1) = B(I,J,L1) + D2*W(I,J) - SUM1*EMU(J)*WT(J)
       enddo
       H(I,L1) = H(I,L1) - 2.0e+0_fp*D2*C(I,L1) + SUM0*ZFLUX
    enddo

    do I = 1,M_
       D1 = EMU(I)/DELTAU
       AA(I,I,L1) = AA(I,I,L1) + D1
       B(I,I,L1)  = B(I,I,L1) + D1
       C(I,L1) = 0.0e+0_fp
    enddo

  END SUBROUTINE GEN_ID
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: jratet
!
! !DESCRIPTION: Subroutine JRATET calculates temperature-dependent J-rates.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE JRATET(PPJ,TTJ,FFF,VALJL,LCTM,LCHEM,NJXU)
!
! !USES:
!
!
! !INPUT PARAMETERS:
!
    integer,  intent(in)    :: LCTM,LCHEM,NJXU
    real(fp), intent(in)    ::  PPJ(JXL1_+1),TTJ(JXL1_+1)
!
! !INPUT/OUTPUT PARAMETERS:
!
    real(fp), intent(inout) ::  FFF(W_,LCTM)
!
! !OUTPUT VARIABLES:
!
    real(fp), intent(out), dimension(LCTM,NJXU) ::  VALJL
!
! !REMARKS:
!
! !REVISION HISTORY:
!  28 Mar 2013 - S. D. Eastham - Copied from Fast-JX v7.0
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    real(fp)  VALJ(X_)
    real(fp)  QO2TOT, QO3TOT, QO31DY, QO31D, QQQT, TFACT
    real(fp)  TT,PP,DD,TT200,TFACA,TFAC0,TFAC1,TFAC2
    real(fp)  QQQA,QQ2,QQ1A,QQ1B
    integer   J,K,L, IV

    !=================================================================
    ! JRATET begins here!
    !=================================================================

    if (NJXU .lt. NJX) then
       call EXITC(' JRATET:  CTM has not enough J-values dimensioned')
    endif
    do L = 1,LCTM
       ! need temperature, pressure, and density at mid-layer
       ! (for some quantum yields):
       TT   = TTJ(L)
       if (L .eq. 1) then
          PP = PPJ(1)
       else
          PP  = (PPJ(L)+PPJ(L+1))*0.5e+0_fp
       endif
       DD = 7.24e18*PP/TT

       ! if W_=18/12, must zero bin-11/5 below 100 hPa, since O2 e-fold is
       ! too weak and does not represent the decay of 215.5-221.5 nm sunlight.
       if (PP .gt. 100.e+0_fp) then
          if (W_ .eq. 18) then
             FFF(11,L) = 0.e+0_fp
          elseif (W_ .eq. 12) then
             FFF(5,L) = 0.e+0_fp
          endif
       endif

       do J = 1,NJXU
          VALJ(J) = 0.e+0_fp
       enddo

       do K = 1,W_
          call X_interp (TT,QO2TOT, TQQ(1,1),QO2(K,1), &
                         TQQ(2,1),QO2(K,2), TQQ(3,1),QO2(K,3), LQQ(1))
          call X_interp (TT,QO3TOT, TQQ(1,2),QO3(K,1), &
                         TQQ(2,2),QO3(K,2), TQQ(3,2),QO3(K,3), LQQ(2))
          call X_interp (TT,QO31DY, TQQ(1,3),Q1D(K,1), &
                         TQQ(2,3),Q1D(K,2), TQQ(3,3),Q1D(K,3), LQQ(3))
          QO31D  = QO31DY*QO3TOT
          VALJ(1) = VALJ(1) + QO2TOT*FFF(K,L)
          VALJ(2) = VALJ(2) + QO3TOT*FFF(K,L)
          VALJ(3) = VALJ(3) + QO31D*FFF(K,L)
       enddo

       do J = 4,NJXU
       do K = 1,W_
          ! also need to allow for Pressure interpolation if SQQ(J) = 'p'
          if (SQQ(J) .eq.'p') then
             call X_interp (PP,QQQT, TQQ(1,J),QQQ(K,1,J), &
                            TQQ(2,J),QQQ(K,2,J), TQQ(3,J),QQQ(K,3,J), LQQ(J))
          else
             call X_interp (TT,QQQT, TQQ(1,J),QQQ(K,1,J), &
                            TQQ(2,J),QQQ(K,2,J), TQQ(3,J),QQQ(K,3,J), LQQ(J))
          endif
          VALJ(J) = VALJ(J) + QQQT*FFF(K,L)
       enddo
       enddo

       do J=1,NJXU
          VALJL(L,J) = VALJ(J)
       enddo

    enddo

    ! Zero non-chemistry layers
    if (LCHEM.lt.LCTM) then
       do L=(LCTM+1),LCHEM
       do J=1,NJXU
          VALJL(L,J) = 0.e+0_fp
       enddo
       enddo
    endif

  END SUBROUTINE JRATET
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: x_interp
!
! !DESCRIPTION: Subroutine X\_INTERP is an up-to-three-point linear interp.
!  function for cross-sections.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE X_INTERP (TINT,XINT,T1,X1,T2,X2,T3,X3,L123)
!
! !USES:
!
!
! !INPUT PARAMETERS:
!
    REAL(fp), INTENT(IN)  ::  TINT,T1,T2,T3, X1,X2,X3
    INTEGER,  INTENT(IN)  ::  L123
!
! !OUTPUT VARIABLES:
!
    REAL(fp), INTENT(OUT) ::  XINT
!
! !REMARKS:
!
! !REVISION HISTORY:
!  28 Mar 2013 - S. D. Eastham - Copied from Fast-JX v7.0
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    REAL(fp)  TFACT

    !=================================================================
    ! X_INTERP begins here!
    !=================================================================

    if (L123 .le. 1) then
       XINT = X1
    elseif (L123 .eq. 2) then
       TFACT = max(0.e+0_fp,min(1.e+0_fp,(TINT-T1)/(T2-T1) ))
       XINT = X1 + TFACT*(X2 - X1)
    else
       if (TINT.le. T2) then
          TFACT = max(0.e+0_fp,min(1.e+0_fp,(TINT-T1)/(T2-T1) ))
          XINT = X1 + TFACT*(X2 - X1)
       else
          TFACT = max(0.e+0_fp,min(1.e+0_fp,(TINT-T2)/(T3-T2) ))
          XINT = X2 + TFACT*(X3 - X2)
       endif
    endif

  END SUBROUTINE X_INTERP
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: sphere2
!
! !DESCRIPTION: Subroutine SPHERE2 is an AMF2.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE SPHERE2 (U0,ZHL,AMF2,L1U,LJX1U)
!
! !INPUT PARAMETERS:
!
    INTEGER,  INTENT(IN)  ::   L1U, LJX1U
    REAL(fp), INTENT(IN)  ::   U0,ZHL(L1U+1)
!
! !OUTPUT VARIABLES:
!
    REAL(fp), INTENT(OUT) ::   AMF2(2*LJX1U+1,2*LJX1U+1)
!
! !REMARKS:
! Quoting from the original:
!  New v6.2: does AirMassFactors for mid-layer, needed for SZA ~ 90
!  This new AMF2 does each of the half-layers of the CTM separately,
!     whereas the original, based on the pratmo code did the whole layers
!     and thus calculated the ray-path to the CTM layre edges, NOT the middle.
!  Since fast-JX is meant to calculate the intensity at the mid-layer, the
!     solar beam at low sun (interpolated between layer edges) was incorrect.
!  This new model does make some approximations of the geometry of the layers:
!     the CTM layer is split evenly in mass (good) and in height (approx).
!                                                                             .
!  Calculation of spherical geometry; derive tangent heights, slant path
!  lengths and air mass factor for each layer. Not called when
!  SZA > 98 degrees.  Beyond 90 degrees, include treatment of emergent
!  beam (where tangent height is below altitude J-value desired at).
!                                                                             .
!  ---------------------------------------------------------------------
!  Inputs:
!     U0      cos(solar zenith angle)
!     RAD  radius of Earth mean sea level (cm)
!     ZHL(L)  height (cm) of the bottom edge of CTM level L
!     ZZHT    scale height (cm) used above top of CTM (ZHL(L_+1))
!     L1U     dimension of CTM = levels +1 (L+1 = above-CTM level)
!  Outputs:
!     AMF2(I,J) = air mass factor for CTM level I for sunlight reaching J
!         ( these are calculated for both layer middle and layer edge)
!  ---------------------------------------------------------------------
!
! !REVISION HISTORY:
!  28 Mar 2013 - S. D. Eastham - Copied from Fast-JX v7.0
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
#ifdef MODEL_GEOS
    INTEGER, PARAMETER  ::  LSPH_ = 200
#else
    INTEGER, PARAMETER  ::  LSPH_ = 100
#endif

    ! RZ      Distance from centre of Earth to each point (cm)
    ! RQ      Square of radius ratios
    ! SHADHT  Shadow height for the current SZA
    ! XL      Slant path between points
    INTEGER  I, J, K, II, L2
    REAL(fp)   XMU1,XMU2,XL,DIFF,SHADHT,RZ(LSPH_+1)
    REAL(fp)   RZ2(2*LSPH_+1),RQ2(2*LSPH_+1)

    !=================================================================
    ! SPHERE2 begins here!
    !=================================================================

    ! must have top-of-atmos (NOT top-of-CTM) defined
    !      ZHL(L1U+1) = ZHL(L1U) + ZZHT

    if (L1U .gt. LSPH_) then
       call EXITC(' SPHERE2: temp arrays not large enough')
    endif

    RZ(1) = RAD + ZHL(1)
    do II = 2,L1U+1
       RZ(II)   = RAD + ZHL(II)
    enddo

    ! calculate heights for edges of split CTM-layers
    L2 = 2*L1U
    do II = 2,L2,2
       I = II/2
       RZ2(II-1) = RZ(I)
       RZ2(II) = 0.5e+0_fp*(RZ(I)+RZ(I+1))
    enddo
    RZ2(L2+1) = RZ(L1U+1)
    do II = 1,L2
       RQ2(II) = (RZ2(II)/RZ2(II+1))**2
    enddo

    ! shadow height for SZA > 90
    if (U0 .lt. 0.0e+0_fp)  then
       SHADHT = RZ2(1)/sqrt(1.0e+0_fp - U0**2)
    else
       SHADHT = 0.e+0_fp
    endif

    ! up from the surface calculating the slant paths between each level
    ! and the level above, and deriving the appropriate Air Mass Factor
    AMF2(:,:) = 0.e+0_fp

    do 16 J = 1,2*L1U+1

       ! Air Mass Factors all zero if below the tangent height
       if (RZ2(J) .lt. SHADHT) goto 16

       ! Ascend from layer J calculating AMF2s
       XMU1 = abs(U0)
       do I = J,2*L1U
          XMU2     = sqrt(1.0e+0_fp - RQ2(I)*(1.0e+0_fp-XMU1**2))
          XL       = RZ2(I+1)*XMU2 - RZ2(I)*XMU1
          AMF2(I,J) = XL / (RZ2(I+1)-RZ2(I))
          XMU1     = XMU2
       enddo

       ! fix above top-of-atmos (L=L1U+1), must set DTAU(L1U+1)=0
       AMF2(2*L1U+1,J) = 1.e+0_fp
       
       ! Twilight case - Emergent Beam, calc air mass factors below layer
       if (U0 .ge. 0.0e+0_fp) goto 16

       ! Descend from layer J
       XMU1       = abs(U0)
       do II = J-1,1,-1
          DIFF        = RZ2(II+1)*sqrt(1.0e+0_fp-XMU1**2)-RZ2(II)
          if (II.eq.1)  DIFF = max(DIFF,0.e+0_fp)   ! filter

          ! Tangent height below current level - beam passes through twice
          if (DIFF .lt. 0.0e+0_fp)  then
             XMU2      = sqrt(1.0e+0_fp - (1.0e+0_fp-XMU1**2)/RQ2(II))
             XL        = abs(RZ2(II+1)*XMU1-RZ2(II)*XMU2)
             AMF2(II,J) = 2.e+0_fp*XL/(RZ2(II+1)-RZ2(II))
             XMU1      = XMU2

          ! Lowest level intersected by emergent beam
          else
             XL        = RZ2(II+1)*XMU1*2.0e+0_fp
             AMF2(II,J) = XL/(RZ2(II+1)-RZ2(II))
             goto 16
          endif
       enddo

16  continue

  END SUBROUTINE SPHERE2
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: extral
!
! !DESCRIPTION: Subroutine EXTRAL adds sub-layers to thick cloud/aerosol layers
!  using log-spacing for sub-layers of increasing thickness ATAU.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE EXTRAL (Input_Opt,DTAUX,L1X,L2X,NX,JXTRA,ILON,ILAT)
!
! !USES:
    USE Input_Opt_Mod,      ONLY : OptInput
!
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN) :: Input_Opt   ! Input options
    INTEGER,        INTENT(IN) :: L1X,L2X     !index of cloud/aerosol
    integer,        intent(in) :: NX          !Mie scattering array size
    real(fp),       intent(in) :: DTAUX(L1X)  !cloud+3aerosol OD in each layer
    integer,        intent(in) :: ILON, ILAT  !lon,lat index
!
! !OUTPUT VARIABLES:
!
    integer,        intent(out):: JXTRA(L2X+1)!number of sub-layers to be added
!
! !REMARKS:
!     DTAUX(L=1:L1X) = Optical Depth in layer L (generally 600 nm OD)
!        This can be just cloud or cloud+aerosol, it is used only to set
!        the number in levels to insert in each layer L
!        Set for log-spacing of tau levels, increasing top-down.
!                                                                             .
!     N.B. the TTAU, etc calculated here are NOT used elsewhere
!                                                                             .
!   The log-spacing parameters have been tested for convergence and chosen
!     to be within 0.5% for ranges OD=1-500, rflect=0-100%, mu0=0.1-1.0
!     use of ATAU = 1.18 and min = 0.01, gives at most +135 pts for OD=100
!     ATAU = 1.12 now recommended for more -accurate heating rates (not J's)
!
! !REVISION HISTORY:
!  28 Mar 2013 - S. D. Eastham - Copied from Fast-JX v7.0
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER JTOTL,I,L,L2
    REAL(fp)  TTAU(L2X+1),DTAUJ, ATAU1,ATAULN,ATAUM,ATAUN1

#ifdef MODEL_GEOS
    ! ckeller, 5/21/18
    LOGICAL  :: failed
    INTEGER  :: N, NMAX
    REAL(fp) :: ATAULOC
#endif

    !=================================================================
    ! EXTRAL begins here!
    !=================================================================

#ifdef MODEL_GEOS
    ! This routine now repeats the extra layer computation with an
    ! increased heating rate (ATAU) if there is an array overfloat.
    ! This is repeated maximum 5 times. The diagnostics arrays
    ! EXTRAL_NLEVS and EXTRAL_NITER archive the number of extra layers
    ! and the number of iterations needed to converge to that solution.
    ! Ideally, NITER is 1 and no adjustments to the heating rate are
    ! needed (ckeller, 5/22/18).
    NMAX = MAX(1,Input_Opt%FJX_EXTRAL_ITERMAX)
    DO N=1,NMAX
       ! local heating rate
       ATAULOC = ATAU + (0.06*(N-1))
#endif

       ! Reinitialize arrays
       TTAU(:)  = 0.e+0_fp
       JXTRA(:) = 0

       ! combine these edge- and mid-layer points into grid of size:
       !     L2X+1 = 2*L1X+1 = 2*L_+3
       ! calculate column optical depths above each level, TTAU(1:L2X+1)
       !     note that TTAU(L2X+1)=0 and TTAU(1)=total OD
       !
       ! Divide thick layers to achieve better accuracy in the scattering code
       ! In the original fast-J, equal sub-layers were chosen, this is wasteful
       ! and this new code (ver 5.3) uses log-scale:
       !     Each succesive layer (down) increase thickness by ATAU > 1
       !     e.g., if ATAU = 2, a layer with OD = 15 could be divided into
       !     4 sub-layers with ODs = 1 - 2 - 4 - 8
       ! The key parameters are:
       !     ATAU = factor increase from one layer to the next
       !     ATAUMN = the smallest OD layer desired
       !     JTAUMX = maximum number of divisions (i.e., may not get to ATAUMN)
       ! These are set in CMN_FJX_MOD, and have been tested/optimized

#if defined( MODEL_GEOS )
       ATAU1  = ATAULOC - 1.e+0_fp
       ATAULN = log(ATAULOC)
#else
       ATAU1  = ATAU - 1.e+0_fp
       ATAULN = log(ATAU)
#endif
       TTAU(L2X+1)  = 0.0e+0_fp

       do L2 = L2X,1,-1
          L         = (L2+1)/2
          DTAUJ     = 0.5e+0_fp * DTAUX(L)
          TTAU(L2)  = TTAU(L2+1) + DTAUJ
          ! Now compute the number of log-spaced sub-layers to be added in
          ! the interval TTAU(L2) > TTAU(L2+1)
          ! The objective is to have successive TAU-layers increasing by factor
          ! ATAU >1 the number of sub-layers + 1
          if (TTAU(L2) .lt. ATAU0) then
             JXTRA(L2) = 0
          else
             ATAUM    = max(ATAU0, TTAU(L2+1))
             ATAUN1 = log(TTAU(L2)/ATAUM) / ATAULN
             JXTRA(L2) = min(JTAUMX, max(0, int(ATAUN1 - 0.5e+0_fp)))
          endif
       enddo

       ! check on overflow of arrays, cut off JXTRA at lower L if too many
       ! levels
#ifdef MODEL_GEOS
       failed   = .FALSE.
       JTOTL    = L2X + 2
       do L2 = L2X,1,-1
          JTOTL  = JTOTL + JXTRA(L2)
          if (JTOTL .gt. NX/2)  then
             failed = .TRUE.
             exit
          endif
       enddo

       ! exit loop if not failed
       if ( .not. failed ) exit
    enddo

    ! print error and cut off JXTRAL at lower L if too many levels
    if ( failed ) then
       IF ( Input_Opt%FJX_EXTRAL_ERR ) THEN
          write(6,'(A,2I5,F9.2)') 'N_/L2_/L2-cutoff JXTRA:',NX,L2X,L2
       ENDIF
       do L = L2,1,-1
          JXTRA(L) = 0
       enddo
       !go to 10
    endif
    !enddo
    !10 continue

    ! Fill diagnostics arrays
    EXTRAL_NLEVS(ILON,ILAT) = SUM(JXTRA(:))
    EXTRAL_NITER(ILON,ILAT) = N
#else
    JTOTL    = L2X + 2
    do L2 = L2X,1,-1
       JTOTL  = JTOTL + JXTRA(L2)
       if (JTOTL .gt. NX/2)  then
          write(6,'(A,2I5,F9.2)') 'N_/L2_/L2-cutoff JXTRA:',NX,L2X,L2
          do L = L2,1,-1
             JXTRA(L) = 0
          enddo
          go to 10
       endif
    enddo
10  continue
#endif

  END SUBROUTINE EXTRAL
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: exitc
!
! !DESCRIPTION: Subroutine EXITC forces an error in GEOS-Chem and quits.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE EXITC (T_EXIT)
!
! !USES:
!
    USE ERROR_MOD, ONLY : ERROR_STOP
!
! !INPUT PARAMETERS:
!
    CHARACTER(LEN=*), INTENT(IN) ::  T_EXIT
!
! !REVISION HISTORY:
!  28 Mar 2013 - S. D. Eastham - Copied from Fast-JX v7.0
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    CALL ERROR_STOP( T_EXIT, 'fast_jx_mod.F90' )

  END SUBROUTINE EXITC
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: int_fjx
!
! !DESCRIPTION: Subroutine INIT\_FJX initializes Fast-JX variables.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE INIT_FJX( Input_Opt, State_Chm, State_Diag, State_Grid, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE Input_Opt_Mod,  ONLY : OptInput
    USE inquireMod,     ONLY : findFreeLUN
    USE State_Chm_Mod,  ONLY : ChmState
    USE State_Chm_Mod,  ONLY : Ind_
    USE State_Diag_Mod, ONLY : DgnState
    USE State_Grid_Mod, ONLY : GrdState
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)  :: Input_Opt   ! Input Options object
    TYPE(ChmState), INTENT(IN)  :: State_Chm   ! Chemistry State object
    TYPE(DgnState), INTENT(IN)  :: State_Diag  ! Diagnostics State object
    TYPE(GrdState), INTENT(IN)  :: State_Grid  ! Grid State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT) :: RC          ! Success or failure?
!
! !REVISION HISTORY:
!  28 Mar 2013 - S. D. Eastham - Copied from Fast-JX v7.0
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    LOGICAL            :: notDryRun
    INTEGER            :: JXUNIT, J, NJXX, PhotoId
    REAL(fp)           :: ND64MULT

    ! Strings
    CHARACTER(LEN=6)   :: TITLEJXX(JVN_)
    CHARACTER(LEN=255) :: DATA_DIR
    CHARACTER(LEN=255) :: FILENAME
    CHARACTER(LEN=255) :: ErrMsg, ThisLoc

    !=================================================================
    ! INIT_FJX begins here!
    !=================================================================

    ! Initialize
    RC          = GC_SUCCESS
    notDryRun   = ( .not. Input_Opt%DryRun )
    ErrMsg      = ''
    ThisLoc     = ' -> at Init_FJX (in module GeosCore/fast_jx_mod.F90)'

    ! Skip these opterations when running in dry-run mode
    IF ( notDryRun ) THEN

       ! Define species IDs
       id_CH2IBr   = IND_('CH2IBr'  )
       id_IBr      = IND_('IBr'     )
       id_CH2ICl   = IND_('CH2ICl'  )
       id_ICl      = IND_('ICl'     )
       id_I2       = IND_('I2'      )
       id_HOI      = IND_('HOI'     )
       id_IO       = IND_('IO'      )
       id_OIO      = IND_('OIO'     )
       id_INO      = IND_('INO'     )
       id_IONO     = IND_('IONO'    )
       id_IONO2    = IND_('IONO2'   )
       id_I2O2     = IND_('I2O2'    )
       id_CH3I     = IND_('CH3i'    )
       id_CH2I2    = IND_('CH2I2'   )
       id_I2O4     = IND_('I2O4'    )
       id_I2O3     = IND_('I2O3'    )

       ! Print info
       IF ( Input_Opt%amIRoot ) THEN
          write(6,*) ' Initializing Fast-JX v7.0 standalone CTM code.'

          if (W_.ne.8 .and. W_.ne.12 .and. W_.ne.18) then
             ErrMsg =  ' INIT_FJX: invalid no. wavelengths'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          endif
       ENDIF

       ! Get a free LUN
       JXUNIT = findFreeLUN()

    ENDIF

    ! Define data directory for FAST-JX input
    DATA_DIR = TRIM( Input_Opt%FAST_JX_DIR )

    !=====================================================================
    ! Read in fast-J X-sections (spectral data)
    !=====================================================================
    FILENAME = TRIM( DATA_DIR ) // 'FJX_spec.dat'

    ! Read file, or just print filename if we are in dry-run mode
    CALL RD_XXX( JXUNIT, TRIM( FILENAME ), Input_Opt, RC)

    ! Trap potential errors
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered in FAST-JX routine "RD_XXX"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    ! Compute factors for UV flux diagnostics
    IF ( notDryRun ) THEN
       IF ( State_Diag%Archive_UVFluxNet      .or. &
            State_Diag%Archive_UVFluxDirect   .or. &
            State_Diag%Archive_UVFluxDiffuse ) THEN
          UVXFACTOR = 0e+0_fp
          ND64MULT  = UVXPLANCK*UVXCCONST*1.0e+13_fp
          DO J = 1, W_
             UVXFACTOR(J) = ND64MULT/WL(J)
          ENDDO
       ENDIF
    ENDIF

    !=====================================================================
    ! Read in 5-wavelength scattering data
    ! (or just print file name if in dry-run mode)
    !=====================================================================
    FILENAME = TRIM( DATA_DIR ) // 'jv_spec_mie.dat'

    ! Read data
    CALL RD_MIE( JXUNIT, TRIM( FILENAME ), Input_Opt, RC )

    ! Trap potential errors
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered in FAST-JX routine "RD_MIE"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    !=====================================================================
    ! Read in AOD data
    ! (or just print file name if in dry-run mode)
    !=====================================================================
    CALL RD_AOD( JXUNIT, Input_Opt, RC )

    ! Trap potential errors
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered in FAST-JX routine "RD_AOD"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    ! Set up MIEDX array to interpret between GC and FJX aerosol indexing
    IF ( notDryRun ) THEN
       CALL SET_AER( Input_Opt )
    ENDIF

    !=====================================================================
    ! Read in T & O3 climatology used to fill e.g. upper layers
    ! or if O3 not calc.
    !=====================================================================
    CALL RD_PROF_NC( Input_Opt, RC )

    ! Trap potential errors
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered in "Rd_Prof_Nc"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    ! Skip if not in dry-run mode
    IF ( notDryRun ) THEN
       NJXX = NJX
       do J = 1,NJX
          TITLEJXX(J) = TITLEJX(J)
       enddo
    ENDIF

    !=====================================================================
    ! Read in photolysis rates used in chemistry code and mapping onto
    ! FJX J's CTM call:  read in J-values names and link to fast-JX names
    !=====================================================================
    FILENAME = TRIM( DATA_DIR ) // 'FJX_j2j.dat'

    ! Read mapping information
    CALL RD_JS_JX( JXUNIT, TRIM( FILENAME ), TITLEJXX, NJXX, &
                   Input_Opt, RC )

    ! Trap potential errors
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered in "Rd_Js_Jx"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    ! Skip further processing if we are in dry-run mode
    IF ( notDryRun ) THEN

       ! Get the GEOS-Chem photolysis index for each of the 1...JVN_ entries
       ! in the FJX_j2j.dat file.  We'll use this for the diagnostics.
       DO J = 1, JVN_

          IF ( J == Rxn_O3_2a ) THEN

             !------------------------------------------------------------
             ! O3 + hv = O + O(1D) branch 1
             !
             ! UCX     : Save this as JO3_O1D in the nPhotol+1 slot
             ! non_UCX : Save this as JO3     in the nPhotol+1 slot
             !------------------------------------------------------------
             GC_Photo_Id(J) = State_Chm%nPhotol + 1

          ELSE IF ( J == Rxn_O3_1 ) THEN

             !------------------------------------------------------------
             ! O3 + hv -> O + O
             !
             ! UCX     : Save this as JO3_O3P in the nPhotol+2 slot
             ! non-UCX : undefined
             !-------------------------------------------------------------
             IF ( Input_Opt%LUCX ) THEN
                GC_Photo_Id(J) = State_Chm%nPhotol + 2
             ELSE
                GC_Photo_Id(J) = -999
             ENDIF

          ELSE IF ( J == Rxn_O3_2b ) THEN

             !------------------------------------------------------------
             ! O3 + hv -> O2 + O(1d) branch 2
             !
             ! UCX     : undefined
             ! non-UCX : Save into the nPhotol+2 slot
             !           NOTE: The JPOH rate in the bpch diagnostic will
             !           now be the sum of the nPhotol+1+nPhotol+2 slots!
             !------------------------------------------------------------
             IF ( Input_Opt%LUCX ) THEN
                GC_Photo_Id(J) = -999
             ELSE
                GC_Photo_Id(J) = State_Chm%nPhotol + 2
             ENDIF

          ELSE

             !------------------------------------------------------------
             ! Everything else
             !
             ! Find the matching GEOS-Chem photolysis species number
             !------------------------------------------------------------
             GC_Photo_Id(J) = Ind_( RNAMES(J), 'P' )

          ENDIF

          ! Print the mapping
          IF ( Input_Opt%amIRoot ) THEN
             IF ( GC_Photo_Id(J) > 0 ) THEN
                WRITE(6, 200) RNAMES(J), J, GC_Photo_Id(J), JFACTA(J)
200             FORMAT( a10, ':', i7, 2x, i7, 2x, f7.4 )
             ENDIF
          ENDIF
       ENDDO
    ENDIF

#ifdef MODEL_GEOS
    ! Diagnostics arrays
    ALLOCATE(EXTRAL_NLEVS(State_Grid%NX,State_Grid%NY))
    ALLOCATE(EXTRAL_NITER(State_Grid%NX,State_Grid%NY))
    EXTRAL_NLEVS(:,:) = 0.0
    EXTRAL_NITER(:,:) = 0.0
#endif

  END SUBROUTINE INIT_FJX
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: rd_xxx
!
! !DESCRIPTION: Subroutine RD\_XXX reads in wavelength bins, solar fluxes,
!  Rayleigh and temperature-dependent cross-sections.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE RD_XXX ( NUN, NAMFIL, Input_Opt, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE Input_Opt_Mod, ONLY : OptInput
!
! !INPUT PARAMETERS:
!
    INTEGER,        INTENT(IN)  :: NUN
    CHARACTER(*),   INTENT(IN)  :: NAMFIL
    TYPE(OptInput), INTENT(IN)  :: Input_Opt
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT) :: RC
!
! !REMARKS:
!    NEW v-6.8  now allow 1 to 3 sets of X-sects for T or P
!           LQQ = 1, 2, or 3 to determine interpolation with T or P
!           IF the temperatures TQQQ are <0, then use as pressure interp (hPa)
!           NB - the temperatures and pressures must be increasing
!    NEW v-6.4  changed to collapse wavelengths & x-sections to Trop-only:
!           WX_ = 18 should match the JX_spec.dat wavelengths
!           W_ = 12 (Trop-only) or 18 (std) is set in (CMN_FJX.F).
!       if W_=12 then drop strat wavels, and drop x-sects (e.g. N2O, ...)
!           W_ = 8, reverts to quick fix:  fast-J (12-18) plus bin (5) scaled
!                                                                             .
!   --------------------------------------------------------------------
!     NAMFIL   Name of spectral data file (FJX_spec.dat) >> j2 for fast-J2
!     NUN      Channel number for reading data file
!
!     NJX    Number of species to calculate J-values for
!     NWWW     Number of wavelength bins, from 1:NWWW
!     WBIN     Boundaries of wavelength bins
!     WL       Centres of wavelength bins - 'effective wavelength'
!     FL       Solar flux incident on top of atmosphere (cm-2.s-1)
!     QRAYL    Rayleigh parameters (effective cross-section) (cm2)
!     QO2      O2 cross-sections
!     QO3      O3 cross-sections
!     Q1D      O3 => O(1D) quantum yield
!     TQQ      Temperature for supplied cross sections
!     QQQ      Supplied cross sections in each wavelength bin (cm2)
!   --------------------------------------------------------------------
!
! !REVISION HISTORY:
!  28 Mar 2013 - S. D. Eastham - Copied from Fast-JX v7.0
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    LOGICAL            :: FileExists
    INTEGER            :: I, J, JJ, K, IW, NQRD, NWWW, LQ
    REAL(fp)           :: TQQ2

    ! Arrays
    REAL(fp)           :: QQ2(199)

    ! Strings
    CHARACTER(LEN=255) :: FileMsg, FileStatus
    CHARACTER(LEN=255) :: ErrMsg,  ThisLoc
    CHARACTER(LEN=78)  :: TITLE0
    CHARACTER(LEN=6 )  :: TITLEJ2, TITLEJ3
    CHARACTER(LEN=1 )  :: TSTRAT

    !=================================================================
    ! In dry-run mode, print file path to dryrun log and exit.
    ! Otherwise, print file path to stdout and continue.
    !=================================================================

    ! Initialize
    RC      = GC_SUCCESS
    ErrMsg  = ''
    ThisLoc = '-> at RD_XXX (in module GeosCore/fast_jx_mod.F)'

    ! Test if the file exists
    INQUIRE( FILE=TRIM( NamFil ), EXIST=FileExists )

    ! Test if the file exists and define an output string
    IF ( FileExists ) THEN
       FileMsg    = 'FAST-JX (RD_XXX): Opening'
    ELSE
       FileMsg    = 'FAST-JX (RD_XXX): REQUIRED FILE NOT FOUND'
    ENDIF

    ! Write to stdout for both regular and dry-run simulations
    IF ( Input_Opt%AmIRoot ) THEN
       WRITE( 6, 300 ) TRIM( FileMsg ), TRIM( NamFil )
300    FORMAT( a, ' ', a )
    ENDIF

    ! For dry-run simulations, return to calling program.
    ! For regular simulations, throw an error if we can't find the file.
    IF ( Input_Opt%DryRun ) THEN
       RETURN
    ELSE
       IF ( .not. FileExists ) THEN
          WRITE( ErrMsg, 300 ) TRIM( FileMsg ), TRIM( NamFil )
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
    ENDIF

    !=================================================================
    ! RD_XXX begins here!
    !=================================================================

    ! Initialization
    TQQ(:,:) = 0.e+0_fp

    ! -------spectral data----set for new format data------------------
    !         note that X_ = max # Xsects read in
    !                   NJX = # fast-JX J-values derived from this (.le. X_)

    ! >>>> W_ = 12 <<<< means trop-only, discard WL #1-4 and #9-10, some X-sects

    ! Open file
    open (NUN,FILE=NAMFIL,status='old',form='formatted')
    read (NUN,100) TITLE0

    ! -note that NQRD is not used any more, a read until 'endofJ' is performed
    read (NUN,101) NQRD,NWWW
    NW1 = 1
    NW2 = NWWW
    IF ( Input_Opt%AmIRoot ) THEN
       write(6,'(1x,a)') TITLE0
       write(6,'(i8)') NWWW
    ENDIF
    ! -J-values:  1=O2, 2=O3P,3=O3D 4=readin Xsects
    read (NUN,102) (WL(IW),IW=1,NWWW)
    read (NUN,102) (FL(IW),IW=1,NWWW)
    read (NUN,102) (QRAYL(IW),IW=1,NWWW)

    ! Read O2 X-sects, O3 X-sects, O3=>O(1D) quant yields (each at 3 temps)
    ! NB the O3 and q-O3-O1D are at different temperatures and cannot be combined
    read (NUN,103) TITLEJX(1),TQQ(1,1), (QO2(IW,1),IW=1,NWWW)
    read (NUN,103) TITLEJ2,   TQQ(2,1), (QO2(IW,2),IW=1,NWWW)
    read (NUN,103) TITLEJ3,   TQQ(3,1), (QO2(IW,3),IW=1,NWWW)

    read (NUN,103) TITLEJX(2),TQQ(1,2), (QO3(IW,1),IW=1,NWWW)
    read (NUN,103) TITLEJ2,   TQQ(2,2), (QO3(IW,2),IW=1,NWWW)
    read (NUN,103) TITLEJ3,   TQQ(3,2), (QO3(IW,3),IW=1,NWWW)

    read (NUN,103) TITLEJX(3),TQQ(1,3), (Q1D(IW,1),IW=1,NWWW)
    read (NUN,103) TITLEJ2,   TQQ(2,3), (Q1D(IW,2),IW=1,NWWW)
    read (NUN,103) TITLEJ3,   TQQ(3,3), (Q1D(IW,3),IW=1,NWWW)

    SQQ(1) = ' '
    SQQ(2) = ' '
    SQQ(3) = ' '

    LQQ(1) = 3
    LQQ(2) = 3
    LQQ(3) = 3

    ! Read remaining species:  X-sections at 1-2-3 T_s
    JJ = 3
    do I=1,9999

       ! try to read in 3 X-sects per J-value (JJ)
       read (NUN,104) TITLEJ2,TSTRAT,TQQ2,(QQ2(IW),IW=1,NWWW)
       if (TITLEJ2 .eq. 'endofJ') goto 1

       ! skip stratosphere only J's (denoted by 'x')if W_<18 => trop-only J's
       if (W_.eq.18 .or. TSTRAT.ne.'x') then
          if (TITLEJ2 .ne. TITLEJX(JJ)) then
             JJ = JJ+1

             if (JJ .gt. X_) then
                call EXITC(' RD_XXX:  X_ not large enough for Xsects read in')
             endif

             TITLEJX(JJ) = TITLEJ2
             LQQ(JJ) = 1
             SQQ(JJ) = TSTRAT
             LQ = LQQ(JJ)
             TQQ(LQ,JJ) = TQQ2
             do IW = 1,NWWW
                QQQ(IW,LQ,JJ) = QQ2(IW)
             enddo
          else
             LQQ(JJ) = LQQ(JJ)+1
             if (LQQ(JJ) .le. 3) then
                LQ = LQQ(JJ)
                TQQ(LQ,JJ) = TQQ2
                do IW = 1,NWWW
                   QQQ(IW,LQ,JJ) = QQ2(IW)
                enddo
             endif
          endif
       endif
    enddo
1   continue
    NJX = JJ

    do J = 1,NJX
       if ( Input_Opt%AmIRoot ) then
          write(6,200) J,TITLEJX(J),SQQ(J),LQQ(J),(TQQ(I,J),I=1,LQQ(J))
       endif
       ! need to check that TQQ is monotonically increasing:
       if (LQQ(J) .eq. 3) then
          if (TQQ(2,J) .ge. TQQ(3,J)) then
             call EXITC ('TQQ out of order')
          endif
          if (TQQ(1,J) .ge. TQQ(2,J)) then
             call EXITC ('TQQ out of order')
          endif
       endif
       if (LQQ(J) .eq. 2) then
          if (TQQ(1,J) .ge. TQQ(2,J)) then
             call EXITC ('TQQ out of order')
          endif
       endif
    enddo

    ! check on doingpressure interp
    ! check on consolidating Qo2 and others into
    ! wrte a newFJX_J2J.dat for mapping on fjx Xsects

    ! truncate number of wavelengths to do troposphere-only
    if (W_ .ne. WX_) then
       ! TROP-ONLY
       if (W_ .eq. 12) then
          if ( Input_Opt%AmIRoot ) then
             write(6,'(a)') &
                  ' >>>TROP-ONLY reduce wavelengths to 12, drop strat X-sects'
          endif
          NW2 = 12
          do IW = 1,4
             WL(IW) = WL(IW+4)
             FL(IW) = FL(IW+4)
             QRAYL(IW) = QRAYL(IW+4)
             do K = 1,3
                QO2(IW,K) = QO2(IW+4,K)
                QO3(IW,K) = QO3(IW+4,K)
                Q1D(IW,K) = Q1D(IW+4,K)
             enddo
             do J = 4,NJX
                do LQ=1,LQQ(J)
                   QQQ(IW,LQ,J) = QQQ(IW+4,LQ,J)
                enddo
             enddo
          enddo
          do IW = 5,12
             WL(IW) = WL(IW+6)
             FL(IW) = FL(IW+6)
             QRAYL(IW) = QRAYL(IW+6)
             do K = 1,3
                QO2(IW,K) = QO2(IW+6,K)
                QO3(IW,K) = QO3(IW+6,K)
                Q1D(IW,K) = Q1D(IW+6,K)
             enddo
             do J = 4,NJX
                do LQ=1,LQQ(J)
                   QQQ(IW,LQ,J) = QQQ(IW+6,LQ,J)
                enddo
             enddo
          enddo
          ! TROP-QUICK  (must scale solar flux for W=5)
       elseif (W_ .eq. 8) then
          if ( Input_Opt%amIRoot ) then
             write(6,'(a)') &
                  ' >>>TROP-QUICK reduce wavelengths to 8, drop strat X-sects'
          endif
          NW2 = 8
          do IW = 1,1
             WL(IW) = WL(IW+4)
             FL(IW) = FL(IW+4)  * 2.e+0_fp
             QRAYL(IW) = QRAYL(IW+4)
             do K = 1,3
                QO2(IW,K) = QO2(IW+4,K)
                QO3(IW,K) = QO3(IW+4,K)
                Q1D(IW,K) = Q1D(IW+4,K)
             enddo
             do J = 4,NJX
                do LQ=1,LQQ(J)
                   QQQ(IW,LQ,J) = QQQ(IW+4,LQ,J)
                enddo
             enddo
          enddo
          do IW = 2,8
             WL(IW) = WL(IW+10)
             FL(IW) = FL(IW+10)
             QRAYL(IW) = QRAYL(IW+10)
             do K = 1,3
                QO2(IW,K) = QO2(IW+10,K)
                QO3(IW,K) = QO3(IW+10,K)
                Q1D(IW,K) = Q1D(IW+10,K)
             enddo
             do J = 4,NJX
                do LQ=1,LQQ(J)
                   QQQ(IW,LQ,J) = QQQ(IW+10,LQ,J)
                enddo
             enddo
          enddo

       else
          call EXITC(' no. wavelengths wrong: W_ .ne. 8,12,18')
       endif
    endif

    close(NUN)

100 format(a)
101 format(10x,5i5)
102 format(10x,    6e10.3/(10x,6e10.3)/(10x,6e10.3))
103 format(a6,1x,f3.0,6e10.3/(10x,6e10.3)/(10x,6e10.3))
104 format(a6,a1,f3.0,6e10.3/(10x,6e10.3)/(10x,6e10.3))
200 format(1x,' x-sect:',i3,a10,a4,i5,3(3x,f6.2))
201 format(' Number of x-sections supplied to Fast-J2: ',i3,/, &
           ' Maximum number allowed (X_) only set to: ',i3,    &
           ' - increase in cmn_FJX.f')

  END SUBROUTINE RD_XXX
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: rd_mie
!
! !DESCRIPTION: Subroutine RD\_MIE retrieves aerosol scattering data for FJX.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE RD_MIE( NUN, NAMFIL, Input_Opt, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE Input_Opt_Mod, ONLY : OptInput
!
! !INPUT PARAMETERS:
!
    INTEGER,        INTENT(IN)  :: NUN         ! Logical unit #
    CHARACTER(*),   INTENT(IN)  :: NAMFIL      ! File name
    TYPE(OptInput), INTENT(IN)  :: Input_Opt   ! Input Options object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT) :: RC          ! Success or failure?
!
! !REMARKS:
!   --------------------------------------------------------------------
!     NAMFIL   Name of scattering data file (e.g., FJX_scat.dat)
!     NUN      Channel number for reading data file
!     NAA      Number of categories for scattering phase functions
!     QAA      Aerosol scattering phase functions
!     WAA      5 Wavelengths for the supplied phase functions
!     PAA      Phase function: first 8 terms of expansion
!     RAA      Effective radius associated with aerosol type
!     SAA      Single scattering albedo
!   --------------------------------------------------------------------
!
! !REVISION HISTORY:
!  28 Mar 2013 - S. D. Eastham - Adapted from GEOS-Chem v9-1-3
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER            :: I, J, K, NK
    LOGICAL            :: LBRC, FileExists

    ! Strings
    CHARACTER(LEN=78 ) :: TITLE0
    CHARACTER(LEN=255) :: FileMsg, ErrMsg, ThisLoc

    !=================================================================
    ! In dry-run mode, print file path to dryrun log and exit.
    ! Otherwise, print file path to stdout and continue.
    !=================================================================

    ! Assume success
    RC = GC_SUCCESS
    ErrMsg = ''
    ThisLoc = ' -> at RD_MIE (in module GeosCore/fast_jx_mod.F90)'

    ! Test if the file exists
    INQUIRE( FILE=TRIM( NamFil ), EXIST=FileExists )

    ! Test if the file exists and define an output string
    IF ( FileExists ) THEN
       FileMsg = 'FAST-JX (RD_MIE): Opening'
    ELSE
       FileMsg = 'FAST-JX (RD_MIE): REQUIRED FILE NOT FOUND'
    ENDIF

    ! Write to stdout for both regular and dry-run simulations
    IF ( Input_Opt%amIRoot ) THEN
       WRITE( 6, 300 ) TRIM( FileMsg ), TRIM( NamFil )
300    FORMAT( a, ' ', a )
    ENDIF

    ! For dry-run simulations, return to calling program
    ! For regular simulations, throw an error if we can't find the file.
    IF ( Input_Opt%DryRun ) THEN
       RETURN
    ELSE
       IF ( .not. FileExists ) THEN
          WRITE( ErrMsg, 300 ) TRIM( FileMsg ), TRIM( NamFil )
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
    ENDIF

    !=================================================================
    ! RD_MIE begins here -- read data from file
    !=================================================================

    ! Copy fields from Input_Opt
    LBRC = Input_Opt%LBRC

    ! Open file
    open (NUN,FILE=NAMFIL,status='old',form='formatted')

    ! Read header lines
    READ( NUN,'(A)' ) TITLE0
    IF  ( Input_Opt%AmIRoot ) WRITE( 6, '(1X,A)' ) TITLE0
    READ( NUN,'(A)' ) TITLE0

    !---Read aerosol phase functions:
    read(NUN,'(A10,I5,/)') TITLE0,NAA

    NK=5
    do j=1,NAA
       read(NUN,110) TITLEAA(j)
       do k=1,NK
          read(NUN,*) WAA(k,j),QAA(k,j),RAA(k,j),SAA(k,j), &
                      (PAA(i,k,j),i=1,8)
       enddo
    enddo

    ! Brown carbon option
    IF (LBRC) THEN

       ! Overwrite OC entries (36-42 in jv_spec_mie.dat)
       ! with BR entries at end of file (labeled 57-63)
       do j= 36, 42
          read(NUN,110) TITLEAA(j)
          do k=1,NK
             read(NUN,*) WAA(k,j),QAA(k,j),RAA(k,j),SAA(k,j), &
                         (PAA(i,k,j),i=1,8)
          enddo
       enddo

    ENDIF

    close(NUN)

    IF ( Input_Opt%amIRoot ) THEN
       write(6,'(a,9f8.1)') ' Aerosol optical: r-eff/rho/Q(@wavel):', &
                            (WAA(K,1),K=1,5)
       do J=1,NAA
          write(6,'(1x,A)') TRIM(TITLEAA(J))
          write(6,'(3x,I2,A,9F8.1)') J,'  wavel=',(WAA(K,J),K=1,NK)
          write(6,'(3x,I2,A,9F8.4)') J,'  Qext =',(QAA(K,J),K=1,NK)
       enddo
    ENDIF

110 format(3x,a80)

  END SUBROUTINE RD_MIE
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: rd_aod
!
! !DESCRIPTION: Subroutine RD\_AOD reads aerosol phase functions that are
!  used to scale diagnostic output to an arbitrary wavelengh.  This
!  facilitates comparing with satellite observations.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE RD_AOD( NJ1, Input_Opt, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE Input_Opt_Mod, ONLY : OptInput
!
! !INPUT PARAMETERS:
!
    INTEGER,          INTENT(IN)  :: NJ1         ! Unit # of file to open
    TYPE(OptInput),   INTENT(IN)  :: Input_Opt   ! Input Options object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(OUT) :: RC          ! Success or failure?
!
! !REMARKS:
!  The .dat files for each species contain the optical properties
!  at multiple wavelengths to be used in the online calculation of the aerosol
!  optical depth diagnostics.
!  These properties have been calculated using the same size and optical
!  properties as the FJX_spec.dat file used for the FAST-J photolysis
!  calculations (which is now redundant for aerosols, the values in the .dat
!  files here are now used). The file currently contains 11 wavelengths
!  for Fast-J and other commonly used wavelengths for satellite and
!  AERONET retrievals. 30 wavelengths follow that map onto RRTMG
!  wavebands for radiaitive flux calculations (not used if RRTMG is off).
!  A complete set of optical properties from 250-2000 nm for aerosols is
!  available at:
!  ftp://ftp.as.harvard.edu/geos-chem/data/aerosol_optics/hi_spectral_res
!                                                                             .
!     -- Colette L. Heald, 05/10/10)
!     -- David A. Ridley, 05/10/13 (update for new optics files)
!
! !REVISION HISTORY:
!  10 May 2010 - C. Heald      - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES
!
    ! Scalars
    INTEGER            :: I, J, K, N
    INTEGER            :: IOS
    LOGICAL            :: LBRC, FileExists

    ! Strings
    CHARACTER(LEN=78 ) :: TITLE0
    CHARACTER(LEN=255) :: DATA_DIR
    CHARACTER(LEN=255) :: THISFILE
    CHARACTER(LEN=255) :: FileMsg
    CHARACTER(LEN=255) :: ErrMsg
    CHARACTER(LEN=255) :: ThisLoc

    ! String arrays
    CHARACTER(LEN=30)  :: SPECFIL(6)
    CHARACTER(LEN=30)  :: SPECFIL_UCX(8)

    !================================================================
    ! RD_AOD begins here!
    !================================================================

    ! Initialize
    RC       = GC_SUCCESS
    ErrMsg   = ''
    ThisLoc  = ' -> at RD_AOD (in module GeosCore/fast_jx_mod.F90)'
    LBRC     = Input_Opt%LBRC
    DATA_DIR = TRIM( Input_Opt%FAST_JX_DIR )

    ! IMPORTANT: aerosol_mod.F and dust_mod.F expect aerosols in this order
    DATA SPECFIL /"so4.dat","soot.dat","org.dat", &
                  "ssa.dat","ssc.dat", "dust.dat"/

    ! For UCX simulations:
    !
    ! Extra two LUT dat files for strat H2SO4 and NAT particles
    !
    ! Treating strat sulfate with GADS data but modified to match
    ! the old Fast-J values size (r=0.09um, sg=0.6) - I think there's
    ! evidence that this is too smale and narrow e.g. Deshler et al. 2003
    ! NAT should really be associated with something like cirrus cloud
    ! but for now we are just treating the NAT like the sulfate... limited
    ! info but ref index is similar e.g. Scarchilli et al. (2005)
    !(DAR 05/2015)
    DATA SPECFIL_UCX /"so4.dat","soot.dat","org.dat", &
                      "ssa.dat","ssc.dat",            &
                      "h2so4.dat","h2so4.dat",        &
                      "dust.dat"/

    ! Loop over the array of filenames
    DO k = 1, NSPAA

       ! Choose different set of input files for UCX and tropchem simulations
       IF ( Input_Opt%LUCX) THEN
          THISFILE = TRIM( DATA_DIR ) // TRIM( SPECFIL_UCX(k) )
       ELSE
          THISFILE = TRIM( DATA_DIR ) // TRIM( SPECFIL(k) )
       ENDIF

       !--------------------------------------------------------------
       ! In dry-run mode, print file path to dryrun log and cycle.
       ! Otherwise, print file path to stdout and continue.
       !--------------------------------------------------------------

       ! Test if the file exists
       INQUIRE( FILE=TRIM( ThisFile ), EXIST=FileExists )

       ! Test if the file exists and define an output string
       IF ( FileExists ) THEN
          FileMsg = 'FAST-JX (RD_AOD): Opening'
       ELSE
          FileMsg = 'FAST-JX (RD_AOD): REQUIRED FILE NOT FOUND'
       ENDIF

       ! Write to stdout for both regular and dry-run simulations
       IF ( Input_Opt%amIRoot ) THEN
          WRITE( 6, 300 ) TRIM( FileMsg ), TRIM( ThisFile )
300       FORMAT( a, ' ', a )
       ENDIF

       ! For dry-run simulations, cycle to next file.
       ! For regular simulations, throw an error if we can't find the file.
       IF ( Input_Opt%DryRun ) THEN
          CYCLE
       ELSE
          IF ( .not. FileExists ) THEN
             WRITE( ErrMsg, 300 ) TRIM( FileMsg ), TRIM( ThisFile )
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF
       ENDIF

       !--------------------------------------------------------------
       ! If not a dry-run, read data from each species file
       !--------------------------------------------------------------

       ! Open file
       OPEN( NJ1, FILE=TRIM( THISFILE ), STATUS='OLD', IOSTAT=RC )

       ! Error check
       IF ( RC /= 0 ) THEN
          ErrMsg = 'Error opening file: ' // TRIM( ThisFile )
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       ! Read header lines
       READ(  NJ1, '(A)' ) TITLE0
       IF ( Input_Opt%amIRoot ) WRITE( 6, '(1X,A)' ) TITLE0

       ! Second header line added for more info
       READ(  NJ1, '(A)' ) TITLE0
       IF ( Input_Opt%amIRoot ) WRITE( 6, '(1X,A)' ) TITLE0

       READ(  NJ1, '(A)' ) TITLE0
110    FORMAT( 3x, a20 )

       DO i = 1, NRAA
       DO j = 1, NWVAA

          READ(NJ1,*) WVAA(j,k),RHAA(i,k),NRLAA(j,i,k),NCMAA(j,i,k), &
                      RDAA(i,k),RWAA(i,k),SGAA(i,k),QQAA(j,i,k),   &
                      ALPHAA(j,i,k),REAA(i,k),SSAA(j,i,k),         &
                      ASYMAA(j,i,k),(PHAA(j,i,k,n),n=1,8)

          ! make note of where 1000nm is for FAST-J calcs
          IF (WVAA(j,k).EQ.1000.0) IWV1000=J

       ENDDO
       ENDDO

       ! Close file
       CLOSE( NJ1 )
    ENDDO

    !=================================================================
    ! Only do the following if we are not running in dry-run mode
    !=================================================================
    IF ( .not. Input_Opt%DryRun ) THEN

       IF ( Input_Opt%amIRoot ) THEN
          WRITE( 6, * ) 'Optics read for all wavelengths successfully'
       ENDIF

       ! Now calculate the required wavelengths in the LUT to calculate
       ! the requested AOD
       CALL CALC_AOD( Input_Opt )
    ENDIF

  END SUBROUTINE RD_AOD
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: calc_aod
!
! !DESCRIPTION: Subroutine CALC\_AOD works out the closest tie points
! in the optics LUT wavelengths and the coefficients required to
! calculate the angstrom exponent for interpolating optics to the requested
! wavelength. If the wavelength requested matches a standard wavelength
! in the LUT then we skip the interpolation (DAR 09/2013)
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE CALC_AOD( Input_Opt )
!
! !USES:
!
    USE CMN_FJX_MOD, ONLY : NWVAA, NWVAA0, WVAA
    USE CMN_FJX_MOD, ONLY : IWVSELECT
    USE CMN_FJX_MOD, ONLY : IRTWVSELECT
    USE CMN_FJX_MOD, ONLY : ACOEF_WV, BCOEF_WV, CCOEF_WV
    USE CMN_FJX_MOD, ONLY : ACOEF_RTWV, BCOEF_RTWV, CCOEF_RTWV
    USE CMN_FJX_MOD, ONLY : NWVREQUIRED, IWVREQUIRED
    USE CMN_FJX_MOD, ONLY : NRTWVREQUIRED, IRTWVREQUIRED
    USE Input_Opt_Mod, ONLY : OptInput
#ifdef RRTMG
    USE PARRRTM,     ONLY : NBNDLW
#endif
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput),   INTENT(IN) :: Input_Opt
!
! !REMARKS:
!  Now the user is able to select any 3 wavelengths for optics
!  output in the input.geos file we need to be able to interpolate
!  to those wavelengths based on what is available in the optics
!  look-up table.
!                                                                             .
!  The standard lookup table currently has values for
!  11 common wavelengths followed by 30 that are required by RRTMG.
!  Only those required to interpolate to user requested
!  wavelengths are selected from the standard wavelengths. RRTMG
!  wavelengths are not used in the interpolation for AOD output
!  (DAR 10/2013)
!                                                                             .
!   UPDATE: because the RT optics output doesnt have access to the
!   standard wavelengths we now calculate two sets of values: one
!   for the ND21 and diag3 outputs that use the standard wavelengths
!   and one for ND72 that interpolates the optics from RRTMG
!   wavelengths. Perhaps a switch needs adding to switch off the RT
!   optics output (and interpolation) if this ends up costing too
!   much and is not used, but it is ideal to have an optics output
!   that matches exactly what RRTMG uses to calculate the fluxes
!
! !REVISION HISTORY:
!  18 Jun 2013 - D. Ridley   - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES
!
    INTEGER             :: MINWV, MAXWV, N, N0, N1, W, NSTEP
    REAL(fp)            :: WVDIF

    !================================================================
    ! CALC_AOD begins here!
    !================================================================

    !cycle over standard wavelengths
    N0=1
    N1=NWVAA0
    NSTEP=1
    NWVREQUIRED=0
    DO W=1,Input_Opt%NWVSELECT
       MINWV     = -999
       MAXWV     =  999
       DO N=N0,N1,NSTEP
          WVDIF = WVAA(N,1)-Input_Opt%WVSELECT(W)
          IF ((WVDIF.LE.0).AND.(WVDIF.GT.MINWV)) THEN
             MINWV = WVDIF
             IWVSELECT(1,W)=N
          ENDIF
          IF ((WVDIF.GE.0).AND.(WVDIF.LT.MAXWV)) THEN
             MAXWV = WVDIF
             IWVSELECT(2,W)=N
          ENDIF
       ENDDO
       IF (IWVSELECT(2,W).EQ.IWVSELECT(1,W)) THEN
          !we have a match!
          MINWV=0
          MAXWV=0
          !add this wavelength to those for output
          NWVREQUIRED=NWVREQUIRED+1
          IWVREQUIRED(NWVREQUIRED)=IWVSELECT(1,W)
       ELSE
          !we are going to have to interpolate to the requested wavelength
          NWVREQUIRED=NWVREQUIRED+1
          IWVREQUIRED(NWVREQUIRED)=IWVSELECT(1,W)
          NWVREQUIRED=NWVREQUIRED+1
          IWVREQUIRED(NWVREQUIRED)=IWVSELECT(2,W)
       ENDIF

       !Error check - ensure we have a match or requested wavelength
       !falls within two LUT tie points
       IF (MINWV.EQ.-999) THEN
          ! requested wavelength is shorter than min wv in LUT
          ! set to min
          write(6,*) 'ERROR requested wavelength is too short!!'
          write(6,*) 'Defaulting to LUT min: ',WVAA(1,1)
          IWVSELECT(1,W)=1
          IWVSELECT(2,W)=1 !300nm
          NWVREQUIRED=NWVREQUIRED-1
          IWVREQUIRED(NWVREQUIRED)=IWVSELECT(1,W)
       ENDIF
       IF (MAXWV.EQ.999) THEN
          ! requested wavelength is longer than min wv in LUT
          ! set to max
          write(6,*) 'ERROR requested wavelength is too long!!'
          write(6,*) 'Defaulting to LUT min: ',WVAA(NWVAA0,1)
          IWVSELECT(1,W)=NWVAA0
          IWVSELECT(2,W)=NWVAA0 !1020nm
          NWVREQUIRED=NWVREQUIRED-1
          IWVREQUIRED(NWVREQUIRED)=IWVSELECT(1,W)
       ENDIF

       !now calcualte the angstrom exponent coefs for interpolation -
       !this is done here to save time and repetition in aerosol_mod.F
       IF (IWVSELECT(1,W).NE.IWVSELECT(2,W)) THEN
          ACOEF_WV(W) = WVAA(IWVSELECT(2,W),1)/Input_Opt%WVSELECT(W)
          BCOEF_WV(W) =1.0d0/(LOG(WVAA(IWVSELECT(2,W),1)/ &
                                  WVAA(IWVSELECT(1,W),1)))
          !relative location of selected wavelength between tie points
          !for interpolating SSA and ASYM for output in aerosol_mod.F and
          !dust_mod.F
          CCOEF_WV(W) =(Input_Opt%WVSELECT(W)-WVAA(IWVSELECT(1,W),1))/ &
                      (WVAA(IWVSELECT(2,W),1)-WVAA(IWVSELECT(1,W),1))
       ENDIF
       IF ( Input_Opt%amIRoot ) THEN
          write(6,*) 'N WAVELENGTHS: ',Input_Opt%NWVSELECT
          write(6,*) 'WAVELENGTH REQUESTED:',Input_Opt%WVSELECT(W)
          write(6,*) 'WAVELENGTH REQUIRED:', NWVREQUIRED
          !write(6,*) IWVSELECT(1,W),WVAA(IWVSELECT(1,W),1)
          !write(6,*) IWVSELECT(2,W),WVAA(IWVSELECT(2,W),1)
          !write(6,*) ACOEF_WV(W),BCOEF_WV(W),CCOEF_WV(W)
          write(6,*) '*********************************'
       ENDIF
    ENDDO !Input_Opt%NWVSELECT
#ifdef RRTMG
    !repeat for RRTMG wavelengths to get the closest wavelength
    !indices and the interpolation coefficients
    !Indices are relative to all wavelengths in the LUT i.e. the RRTMG
    !wavelengths start at NWVAA0+1
    N0=NWVAA0+1
    N1=NWVAA
    NSTEP=1
    NRTWVREQUIRED=0
    DO W=1,Input_Opt%NWVSELECT
       MINWV     = -999
       MAXWV     =  999
       DO N=N0,N1,NSTEP
          WVDIF = WVAA(N,1)-Input_Opt%WVSELECT(W)
          IF ((WVDIF.LE.0).AND.(WVDIF.GT.MINWV)) THEN
             MINWV = WVDIF
             IRTWVSELECT(1,W)=N
          ENDIF
          IF ((WVDIF.GE.0).AND.(WVDIF.LT.MAXWV)) THEN
             MAXWV = WVDIF
             IRTWVSELECT(2,W)=N
          ENDIF
       ENDDO
       IF (IRTWVSELECT(2,W).EQ.IRTWVSELECT(1,W)) THEN
          !we have a match!
          MINWV=0
          MAXWV=0
          !add this wavelength to those for output
          NRTWVREQUIRED=NRTWVREQUIRED+1
          IRTWVREQUIRED(NRTWVREQUIRED)=IRTWVSELECT(1,W)
       ELSE
          !we are going to have to interpolate to the requested
          !wavelength
          NRTWVREQUIRED=NRTWVREQUIRED+1
          IRTWVREQUIRED(NRTWVREQUIRED)=IRTWVSELECT(1,W)
          NRTWVREQUIRED=NRTWVREQUIRED+1
          IRTWVREQUIRED(NRTWVREQUIRED)=IRTWVSELECT(2,W)
       ENDIF

       !Error check - ensure we have a match or requested wavelength
       !falls within two LUT tie points
       IF (MINWV.EQ.-999) THEN
          ! requested wavelength is shorter than min wv in LUT
          ! set to min
          write(6,*) 'ERROR requested wavelength is too short!!'
          write(6,*) 'Defaulting to LUT min: ',WVAA(NWVAA-1,1)
          IRTWVSELECT(1,W)=NWVAA-1
          IRTWVSELECT(2,W)=NWVAA-1
          NRTWVREQUIRED=NRTWVREQUIRED-1
          IRTWVREQUIRED(NRTWVREQUIRED)=IRTWVSELECT(1,W)
       ENDIF
       IF (MAXWV.EQ.999) THEN
          ! requested wavelength is longer than min wv in LUT
          ! set to max
          write(6,*) 'ERROR requested wavelength is too long!!'
          write(6,*) 'Defaulting to LUT min: ',WVAA(NWVAA0+1,1)
          IRTWVSELECT(1,W)=NWVAA0+1
          IRTWVSELECT(2,W)=NWVAA0+1
          NRTWVREQUIRED=NRTWVREQUIRED-1
          IRTWVREQUIRED(NRTWVREQUIRED)=IRTWVSELECT(1,W)
       ENDIF

       !now calcualte the angstrom exponent coefs for interpolation -
       !this is done here to save time and repetition in aerosol_mod.F
       IF (IRTWVSELECT(1,W).NE.IRTWVSELECT(2,W)) THEN
          ACOEF_RTWV(W) = WVAA(IRTWVSELECT(2,W),1)/Input_Opt%WVSELECT(W)
          BCOEF_RTWV(W) =1.0d0/(LOG(WVAA(IRTWVSELECT(2,W),1)/ &
                                    WVAA(IRTWVSELECT(1,W),1)))
          !relative location of selected wavelength between tie points
          !for interpolating SSA and ASYM for output in aerosol_mod.F and
          !dust_mod.F
          CCOEF_RTWV(W) =(Input_Opt%WVSELECT(W)-WVAA(IRTWVSELECT(1,W),1))/ &
                      (WVAA(IRTWVSELECT(2,W),1)-WVAA(IRTWVSELECT(1,W),1))
       ENDIF
       !convert wavelength index to that required by rrtmg_rad_transfer
       !i.e. without the standard and LW wavelengths
       IRTWVSELECT(1,W) = IRTWVSELECT(1,W) - NWVAA0 - NBNDLW
       IRTWVSELECT(2,W) = IRTWVSELECT(2,W) - NWVAA0 - NBNDLW
       IF ( Input_Opt%amIRoot ) THEN
          write(6,*) 'N RT WAVELENGTHS: ',Input_Opt%NWVSELECT
          write(6,*) 'RT WAVELENGTH REQUESTED:',Input_Opt%WVSELECT(W)
          write(6,*) 'RT WAVELENGTH REQUIRED:', NRTWVREQUIRED
          write(6,*) IRTWVSELECT(1,W),WVAA(IRTWVSELECT(1,W)+NWVAA0+NBNDLW,1)
          write(6,*) IRTWVSELECT(2,W),WVAA(IRTWVSELECT(2,W)+NWVAA0+NBNDLW,1)
          write(6,*) ACOEF_WV(W),BCOEF_WV(W),CCOEF_WV(W)
          write(6,*) '*********************************'
       ENDIF
    ENDDO !Input_Opt%NWVSELECT
#endif
  END SUBROUTINE CALC_AOD
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: rd_js_jx
!
! !DESCRIPTION: Subroutine RD\_JS\_JX reads in 'FJX\_j2j.dat', which defines
!  the mapping of Fast-JX J's (TITLEJX(1:NJX)) onto the CTM reactions.
!  Reaction number JJ, named T\_REACT, uses Fast-JX's T\_FJX (including scaling
!  factor F\_FJX).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE RD_JS_JX( NUNIT, NAMFIL, TITLEJX, NJXX, Input_Opt, RC )
!
! !USES:
!
    USE Charpak_Mod,   ONLY : CStrip
    USE ErrCode_Mod
    USE Input_Opt_Mod, ONLY : OptInput
!
! !INPUT PARAMETERS:
!
    INTEGER,          INTENT(IN)                  :: NUNIT
    INTEGER,          INTENT(IN)                  :: NJXX
    CHARACTER(LEN=*), INTENT(IN)                  :: NAMFIL
    CHARACTER(LEN=6), INTENT(IN), DIMENSION(NJXX) :: TITLEJX
    TYPE(OptInput),   INTENT(IN)                  :: Input_Opt
!
! !OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(OUT)                 :: RC
!
! !REMARKS:
!  Now flag special reactions that are to be adjusted for FlexChem later.
!
! !REVISION HISTORY:
!  28 Mar 2013 - S. D. Eastham - Copied from Fast-JX v7.0
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    LOGICAL            :: FileExists
    INTEGER            :: J, JJ, K
    REAL(fp)           :: F_FJX

    ! Strings
    CHARACTER(LEN=6  ) :: T_FJX
    CHARACTER(LEN=50 ) :: T_REACT
    CHARACTER(LEN=50 ) :: TEXT
    CHARACTER(LEN=120) :: CLINE
    CHARACTER(LEN=255) :: ErrMsg, ThisLoc, FileMsg

    ! String arrays
    CHARACTER(LEN=6)   :: JMAP(JVN_)

    !=================================================================
    ! RD_JS_JX begins here!
    !=================================================================

    ! Initialize
    RC       = GC_SUCCESS
    ErrMsg   = ''
    ThisLoc  = ' -> at Rd_Js_Jx (in module GeosCore/fast_jx_mod.F90)'

    !=================================================================
    ! In dry-run mode, print file path to dryrun log and exit.
    ! Otherwise, print file path to stdout and continue.
    !=================================================================

    ! Test if the file exists
    INQUIRE( FILE=TRIM( NamFil ), EXIST=FileExists )

    ! Test if the file exists and define an output string
    IF ( FileExists ) THEN
       FileMsg = 'FAST-JX (RD_JS_JX): Opening'
    ELSE
       FileMsg = 'FAST-JX (RD_JS_JX): REQUIRED FILE NOT FOUND'
    ENDIF

    ! Write to stdout for both regular and dry-run simulations
    IF ( Input_Opt%amIRoot ) THEN
       WRITE( 6, 300 ) TRIM( FileMsg ), TRIM( NamFil )
300    FORMAT( a, ' ', a )
    ENDIF

    ! For dry-run simulations, return to calling program.
    ! For regular simulations, throw an error if we can't find the file.
    IF ( Input_Opt%DryRun ) THEN
       RETURN
    ELSE
       IF ( .not. FileExists ) THEN
          WRITE( ErrMsg, 300 ) TRIM( FileMsg ), TRIM( NamFil )
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
    ENDIF

    !=================================================================
    ! Read the FJX_j2j.dat file to map model specific J's onto fast-JX J's
    ! The chemistry code title describes fully the reaction (a50)
    ! Blank (unfilled) chemistry J's are unmapped
    ! The number NRATJ is the last JJ readin that is .le. JVN
    ! include fractional quantum yield for the fast-JX J's
    !=================================================================

    JLABEL(:) = '------'
    JMAP(:)   = '------'
    JFACTA(:) = 0.e+0_fp

    ! Open file
    open (NUNIT,file=NAMFIL,status='old',form='formatted')

    read (NUNIT,'(a)') CLINE
    IF ( Input_Opt%amIRoot ) THEN
       write(6,'(a)') CLINE
    ENDIF
    do J = 1,JVN_
       read (NUNIT,'(i4,1x,a50,4x,f5.3,2x,a6)') JJ,T_REACT,F_FJX,T_FJX
       if (JJ.gt.JVN_) goto 20
       JLABEL(JJ) = T_REACT
       JFACTA(JJ) = F_FJX
       JMAP(JJ) = T_FJX
       NRATJ = JJ
       ! SDE 03/31/13: Check number of branches
       ! Note that the order of the branches in
       ! globchem.dat must match the order in
       ! FJX_j2j.dat
       READ (T_REACT(1:10),"(a10)") RNAMES(JJ)
       RNAMES(JJ) = TRIM(RNAMES(JJ))
       BRANCH(JJ) = 1
       DO K=1,(JJ-1)
          IF (RNAMES(JJ) == RNAMES(K)) THEN
             BRANCH(JJ) = BRANCH(K) + 1
          ENDIF
       ENDDO
    enddo

20  close(NUNIT)

    ! Zero / Set index arrays that map Jvalue(j) onto rates
    do K = 1,NRATJ
       JIND(K) = 0
       do J = 1,NJXX
          T_FJX = TITLEJX(J)
          if (JMAP(K) .eq. TITLEJX(J)) then
             JIND(K)=J
          endif
       enddo
    enddo

    IF ( Input_Opt%amIRoot ) THEN
       write(6,'(a,i4,a)')'Photochemistry Scheme with',NRATJ,' J-values'
    ENDIF
    do K=1,NRATJ
       if (JMAP(K) .ne. '------' ) then
          J = JIND(K)
          IF ( Input_Opt%amIRoot ) THEN
             if (J.eq.0) then
                write(6,'(i5,1x,a50,f6.3,a,1x,a6)') K,JLABEL(K),JFACTA(K), &
                     ' no mapping onto fast-JX',JMAP(K)
             else
                write(6,'(i5,1x,a50,f6.3,a,i4,1x,a6)') K,JLABEL(K),JFACTA(K), &
                     ' mapped to FJX:',J,TITLEJX(J)
             endif
          ENDIF
       endif
    enddo

    !=================================================================
    ! Flag special reactions that will be later adjusted by
    ! routine PHOTRATE_ADJ (called from FlexChem)
    !=================================================================

    ! Loop over all photolysis reactions
    DO K = 1, NRATJ

       ! Strip all blanks from the reactants and products list
       TEXT = JLABEL(K)
       CALL CSTRIP( TEXT )

       !IF ( Input_Opt%amIRoot ) THEN
       !   WRITE(*,*) K, TRIM( TEXT )
       !ENDIF

       ! Look for certain reactions
       SELECT CASE( TRIM( TEXT ) )

       ! O2 + hv -> O + O
       CASE( 'O2PHOTONOO' )
          RXN_O2 = K

       ! O3 + hv -> O2 + O
       CASE( 'O3PHOTONO2O' )
          RXN_O3_1 = K

       ! O3 + hv -> O2 + O(1D)
       CASE( 'O3PHOTONO2O(1D)' )

          ! NOTE: There are 2 reactions of this form.  We shall save
          ! the first one that is encountered in RXN_O3_2a and the
          ! second one in RXN_O3_2b. (bmy, 3/29/16)
          IF ( RXN_O3_2a > 0 ) THEN
             RXN_O3_2b = K
          ELSE
             RXN_O3_2a = K
          ENDIF

       ! SO4 + hv -> SO2 + OH + OH
       CASE( 'SO4PHOTONSO2OHOH' )
          RXN_H2SO4 = K

       ! NO2 + hv -> NO + O
       CASE( 'NO2PHOTONNOO' )
          RXN_NO2 = K

       ! NO + hv -> N + O
       CASE( 'NOPHOTONNO' )
          RXN_NO = K

       ! NO3 + hv -> NO2 + O
       CASE( 'NO3PHOTONNO2O' )
          RXN_NO3 = K

       ! N2O + hv -> N2 + O
       CASE( 'N2OPHOTONN2O' )
          RXN_N2O = K

       ! NITs + hv -> HNO2
       CASE( 'NITsPHOTONHNO2' )
          RXN_JNITSa = K

       ! NITs + hv -> NO2
       CASE( 'NITsPHOTONNO2' )
          RXN_JNITSb = K

       ! NIT + hv -> HNO2
       CASE( 'NITPHOTONHNO2' )
          RXN_JNITa = K

       ! NIT + hv -> NO2
       CASE( 'NITPHOTONNO2' )
          RXN_JNITb = K

       ! HNO3 + hv = OH + NO2
       CASE( 'HNO3PHOTONNO2OH' )
          RXN_JHNO3 = K

       CASE DEFAULT
          ! Nothing
       END SELECT

    ENDDO

    !---------------------------------------------------------------------
    ! Error check the various rxn flags
    !---------------------------------------------------------------------
    IF ( RXN_O2 < 0 ) THEN
       ErrMsg = 'Could not find rxn O2 + hv -> O + O'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    IF ( RXN_O3_1 < 0 ) THEN
       ErrMsg = 'Could not find rxn O3 + hv -> O2 + O'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    IF ( RXN_O3_2a < 0 ) THEN
       ErrMsg = 'Could not find rxn O3 + hv -> O2 + O(1D) #1'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    IF ( RXN_O3_2b  < 0 ) THEN
       ErrMsg = 'Could not find rxn O3 + hv -> O2 + O(1D) #2'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
    ENDIF

    IF ( RXN_NO2 < 0 ) THEN
       ErrMsg = 'Could not find rxn NO2 + hv -> NO + O'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    IF ( RXN_NO2 < 0 ) THEN
       ErrMsg = 'Could not find rxn NO2 + hv -> NO + O'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    IF ( RXN_JNITSa < 0 ) THEN
       ErrMsg = 'Could not find rxn NITS + hv -> HNO2'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    IF ( RXN_JNITSb < 0 ) THEN
       ErrMsg = 'Could not find rxn NITS + hv -> NO2'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    IF ( RXN_JNITa < 0 ) THEN
       ErrMsg = 'Could not find rxn NIT + hv -> HNO2'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    IF ( RXN_JNITb < 0 ) THEN
       ErrMsg = 'Could not find rxn NIT + hv -> NO2'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    !---------------------------------------------------------------------
    ! These reactions are only defined for the UCX mechanism!
    !---------------------------------------------------------------------
    IF ( Input_Opt%LUCX ) THEN

       IF ( RXN_H2SO4  < 0 ) THEN
          ErrMsg = 'Could not find rxn SO4 + hv -> SO2 + OH + OH!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       IF ( RXN_NO3 < 0 ) THEN
          ErrMsg = 'Could not find rxn NO3 + hv -> NO2 + O'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       IF ( RXN_NO < 0 ) THEN
          ErrMsg = 'Could not find rxn NO + hv -> O + N'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       IF ( RXN_N2O < 0 ) THEN
          ErrMsg = 'Could not find rxn N2O + hv -> N2 + O(1D)'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

    ENDIF

    !------------------------------------
    ! Print out saved rxn flags
    !------------------------------------
    IF ( Input_Opt%amIRoot ) THEN
       WRITE( 6, 100 ) REPEAT( '=', 79 )
       WRITE( 6, 110 )
       WRITE( 6, 120 ) RXN_O2
       WRITE( 6, 130 ) RXN_O3_1
       WRITE( 6, 140 ) RXN_O3_2a
       WRITE( 6, 150 ) RXN_O3_2b
       WRITE( 6, 180 ) RXN_JNITSa
       WRITE( 6, 190 ) RXN_JNITSb
       WRITE( 6, 200 ) RXN_JNITa
       WRITE( 6, 210 ) RXN_JNITb
       IF ( Input_Opt%LUCX ) THEN
          WRITE( 6, 160 ) RXN_H2SO4
       ENDIF
       WRITE( 6, 170 ) RXN_NO2
       WRITE( 6, 100 ) REPEAT( '=', 79 )
    ENDIF

    ! FORMAT statements
100 FORMAT( a                                                 )
110 FORMAT( 'Photo rxn flags saved for use in PHOTRATE_ADJ:', / )
120 FORMAT( 'RXN_O2    [ O2  + hv -> O + O         ]  =  ', i5 )
130 FORMAT( 'RXN_O3_1  [ O3  + hv -> O2 + O        ]  =  ', i5 )
140 FORMAT( 'RXN_O3_2a [ O3  + hv -> O2 + O(1D) #1 ]  =  ', i5 )
150 FORMAT( 'RXN_O3_2b [ O3  + hv -> O2 + O(1D) #2 ]  =  ', i5 )
160 FORMAT( 'RXN_H2SO4 [ SO4 + hv -> SO2 + OH + OH ]  =  ', i5 )
170 FORMAT( 'RXN_NO2   [ NO2 + hv -> NO + O        ]  =  ', i5 )
180 FORMAT( 'RXN_JNITSa [ NITS + hv -> HNO2        ]  =  ', i5 )
190 FORMAT( 'RXN_JNITSb [ NITS + hv -> NO2         ]  =  ', i5 )
200 FORMAT( 'RXN_JNITa  [ NIT + hv -> HNO2         ]  =  ', i5 )
210 FORMAT( 'RXN_JNITb  [ NIT + hv -> NO2          ]  =  ', i5 )

  END SUBROUTINE RD_JS_JX
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: solar_jx
!
! !DESCRIPTION: Subroutine SOLAR\_JX handles solar zenith angles.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE SOLAR_JX(NDAY,COSSZA,SZA,SOLFX)
!
! !USES:
!
    USE PhysConstants,ONLY : PI
!
! !INPUT PARAMETERS:
!
    REAL(fp), INTENT(IN)  ::  COSSZA
    INTEGER,  INTENT(IN)  ::  NDAY
!
! !OUTPUT VARIABLES:
!
    REAL(fp), INTENT(OUT) ::  SZA,SOLFX
!
! !REMARKS:
!  ---------------------------------------------------------------------
!     NDAY   = integer day of the year (used for solar lat and declin)
!     SZA = solar zenith angle in degrees
!     COSSZA = U0 = cos(SZA)
!     SOLFX = Solar function
!  ---------------------------------------------------------------------
!
! !REVISION HISTORY:
!  28 Mar 2013 - S. D. Eastham - Adapted from Fast-JX v7.0
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    REAL(fp)  PI180

    !=================================================================
    ! SOLAR_JX begins here!
    !=================================================================

    PI180  = PI/180.e+0_fp
    SZA    = acos(MIN(MAX(COSSZA,-1._fp),1._fp))/PI180

    ! Offset used for GEOS-Chem slightly different
    !SOLFX  = 1.e+0_fp-(0.034e+0_fp*cos(dble(NDAY-186)*2.e+0_fp*PI/365.e+0_fp))
    SOLFX  = 1.e+0_fp-(0.034e+0_fp*cos(dble(NDAY-172) &
            *2.e+0_fp*PI/365.e+0_fp))

  END SUBROUTINE SOLAR_JX
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: photo_jx
!
! !DESCRIPTION: Subroutine PHOTO\_JX is the core subroutine of Fast-JX.
!    calc J's for a single column atmosphere (aka Indep Colm Atmos or ICA)
!    needs P, T, O3, clds, aersls; adds top-of-atmos layer from climatology
!    needs day-fo-year for sun distance, SZA (not lat or long)
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE PHOTO_JX( U0,             REFLB,     P_COL,       &
                       T_COL,          O3_COL,    O3_TOMS,     &
                       AOD999,         ODAER_COL, ODMDUST_COL, &
                       ODCLOUD_COL_IN, ILON,      ILAT,        &
                       YLAT,           DAY_OF_YR, MONTH,       &
                       DAY,            Input_Opt, State_Diag,  &
                       State_Grid,     State_Met )
!
! !USES:
!
    USE CMN_SIZE_Mod,   ONLY : NRH, NRHAER
    USE Input_Opt_Mod,  ONLY : OptInput
    USE State_Diag_Mod, ONLY : DgnState
    USE State_Grid_Mod, ONLY : GrdState
    USE State_Met_Mod,  ONLY : MetState
!
! !INPUT PARAMETERS:
!
    INTEGER,        INTENT(IN)                      :: ILON, ILAT
    INTEGER,        INTENT(IN)                      :: DAY_OF_YR
    INTEGER,        INTENT(IN)                      :: MONTH,DAY
    REAL(fp),       INTENT(IN)                      :: YLAT
    REAL(fp),       INTENT(IN)                      :: U0,REFLB
    REAL(fp),       INTENT(IN), DIMENSION(L1_+1   ) :: P_COL
    REAL(fp),       INTENT(IN), DIMENSION(L1_     ) :: T_COL
    REAL(fp),       INTENT(IN), DIMENSION(L1_     ) :: O3_COL
    REAL(fp),       INTENT(IN)                      :: O3_TOMS
    REAL(fp),       INTENT(IN), DIMENSION(L_,A_   ) :: ODAER_COL
    REAL(fp),       INTENT(IN), DIMENSION(L_,NDUST) :: ODMDUST_COL
    REAL(fp),       INTENT(IN), DIMENSION(L_      ) :: ODCLOUD_COL_IN
    LOGICAL,        INTENT(IN)                      :: AOD999
    TYPE(OptInput), INTENT(IN)                      :: Input_Opt
    TYPE(GrdState), INTENT(IN)                      :: State_Grid
    TYPE(MetState), INTENT(IN)                      :: State_Met
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(DgnState), INTENT(INOUT)                   :: State_Diag
!
! !REMARKS:
!
! !REVISION HISTORY:
!  28 Mar 2013 - S. D. Eastham - Copied from Fast-JX v7.0
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! --------------------------------------------------------------------
    ! key LOCAL atmospheric data needed to solve plane-parallel J----
    ! --these are dimensioned JXL_, and must have JXL_ .ge. State_Grid%NZ
    real(fp), dimension(JXL1_)      :: DDJ,OOJ
    real(fp), dimension(JXL1_+1)    :: PPJ,ZZJ
    integer,dimension(JXL2_+1)      :: JXTRA
    real(fp), dimension(W_)         :: FJTOP,FJBOT,FSBOT,FLXD0,RFL
    real(fp), dimension(JXL_, W_)   :: AVGF, FJFLX
    real(fp), dimension(JXL1_,W_)   :: DTAUX, FLXD
    real(fp), dimension(8,JXL1_,W_) :: POMEGAX

    ! flux/heating arrays (along with FJFLX,FLXD,FLXD0)
    real(fp) :: ODABS,ODRAY
    real(fp) :: RFLECT
    real(fp) :: AMF2(2*JXL1_+1,2*JXL1_+1)

    ! ---------key SCATTERING arrays for clouds+aerosols------------------
    real(fp) :: OD(5,JXL1_),SSA(5,JXL1_),SLEG(8,5,JXL1_)
    real(fp) :: OD600(JXL1_)

    ! ---------key arrays AFTER solving for J's---------------------------
    real(fp) :: FFF(W_,JXL_),VALJ(X_)
    real(fp) :: VALJL(JXL_,X_) !2-D array of J_s returned by JRATET

    integer  :: L2EDGE, I,J,K,L,M,KMIE,IDXAER,IM
    INTEGER  :: KMIE2, IR
    real(fp) :: XQO3,XQO2,WAVE, TTTX
    ! --------------------------------------------------------------------
    ! For compatibility with GEOS-Chem (SDE 03/30/13)
    REAL(fp) :: QSCALING,LOCALOD,LOCALSSA,SZA,SOLF
    REAL(fp) :: AERX_COL(A_,L1_) ! Accumulated aerosol array
    REAL(fp) :: ODCLOUD_COL(L_)
    REAL(fp) :: VALJXX(L_,JVN_)

    ! T_CLIM: Climatology data temperature (K)
    REAL(fp) :: T_CLIM(JXL1_)

    ! T_INPUT: Input temperature (K) with extra layer for compatibility
    REAL(fp) :: T_INPUT(JXL1_+1)

    ! UVFlux* diagnostics
    REAL(fp) :: FDIRECT (JXL1_)
    REAL(fp) :: FDIFFUSE(JXL1_)
    REAL(fp) :: UVX_CONST

    ! Logical flags
    LOGICAL :: IS_HALOGENS

    !Maps the new LUT optics wavelengths on to
    !the 5 jv_spec_mie.dat wavelengths
    ! N.B. currently 200nm and 300nm data is the same in
    ! jv_spec_mie.dat, so we copy from new LUT 300nm for both
    INTEGER :: LUTIDX(5)
    LUTIDX = (/1,1,2,6,10/)

    !=================================================================
    ! PHOTO_JX begins here!
    !=================================================================

    if (State_Grid%NZ+1 .gt. JXL1_) then
       call EXITC(' PHOTO_JX: not enough levels in JX')
    endif

    ! Define logical flags for ND22
    IS_HALOGENS = ( id_CH2IBr > 0 .and. id_IBr   > 0 .and. &
                    id_CH2ICl > 0 .and. id_ICl   > 0 .and. &
                    id_I2     > 0 .and. id_HOI   > 0 .and. &
                    id_IO     > 0 .and. id_OIO   > 0 .and. &
                    id_INO    > 0 .and. id_IONO  > 0 .and. &
                    id_IONO2  > 0 .and. id_I2O2  > 0 .and. &
                    id_CH3I   > 0 .and. id_CH2I2 > 0 .and. &
                    id_I2O4   > 0 .and. id_I2O3  > 0 )

    ! Copy cloud OD data to a variable array
    DO L=1,L_
       ODCLOUD_COL(L) = ODCLOUD_COL_IN(L)
    ENDDO

    ! Input conversion (SDE 03/29/13)
    ! Calculate solar zenith angle (degrees)
    CALL SOLAR_JX(DAY_OF_YR,U0,SZA,SOLF)

    L2EDGE = L_ + L_ + 2
    FFF(:,:) = 0.e+0_fp

    ! check for dark conditions SZA > 98.0 deg => tan ht = 63 km or
    !                                 99.0                 80 km
    if (SZA .gt. 98.e+0_fp) goto 99

    ! Use GEOS-Chem methodology to set vertical profiles of:
    ! Pressure      (PPJ)    [hPa]
    ! Temperature   (T_CLIm) [K]
    ! Path density  (DDJ)    [# molec/cm2]
    ! New methodology for:
    ! Ozone density (OOJ)    [# O3 molec/cm2]
    CALL SET_PROF (YLAT,        MONTH,     DAY,         &
                   T_COL,       P_COL,     ODCLOUD_COL, &
                   ODMDUST_COL, ODAER_COL, O3_COL,      &
                   O3_TOMS,     AERX_COL,  T_CLIM,      &
                   OOJ,         ZZJ,       DDJ,         &
                   Input_Opt,   State_Grid )

    ! Fill out PPJ and TTJ with CTM data to replace fixed climatology
    DO L=1,L1_
       PPJ(L) = P_COL(L)
       T_INPUT(L) = T_COL(L)
    ENDDO

    ! Ensure TOA pressure is zero
    PPJ(L1_+1) = 0.e+0_fp
    T_INPUT(L1_+1) = T_CLIM(L1_)

    ! calculate spherical weighting functions (AMF: Air Mass Factor)
    RFLECT = REFLB

    ! --------------------------------------------------------------------
    call SPHERE2 (U0,ZZJ,AMF2,L1_,JXL1_)
    ! --------------------------------------------------------------------

    !-----------------------------------------------------------------------
    ! Modification for GEOS-Chem: Optical depths are calculated at a single
    ! wavelength for GEOS-Chem, so we perform scaling in this routine.
    ! elsewhere.
    ! (SDE 03/31/13)
    !-----------------------------------------------------------------------
    ! calculate the optical properties (opt-depth, single-scat-alb,
    ! phase-fn(1:8)) at the 5 std wavelengths 200-300-400-600-999 nm
    ! for cloud+aerosols
    do L = 1,L1_
       OD600(L) = 0.e+0_fp
       ! ODs stored with fine-grain (DTAUX) and coarse (OD)
       do K=1,W_
          DTAUX(L,K) = 0.e+0_fp
       enddo
       do K=1,5
          OD(K,L)  = 0.e+0_fp
          SSA(K,L) = 0.e+0_fp
          do I=1,8
             SLEG(I,K,L) = 0.e+0_fp
          enddo
       enddo
    enddo

    ! Clunky fix to accomodate RRTMG and UCX (DAR 01/2015)
    ! Using a combination of old optics LUT format and
    ! new optics format (greater spectral resolution and range
    ! for RRTMG). Clouds, non-species aerosols and strat aerosols
    ! are not incorporated into the new LUT so must still use the
    ! old LUT for these.

    ! Don't bother on extraneous level - leave as zero
    do L = 1,L_
       DO KMIE=1,5
          ! Clouds and non-species aerosols
          DO M=1,3
             IF (AERX_COL(M,L).gt.0e+0_fp) THEN
                IDXAER=MIEDX(M)
                ! Cloud (600 nm scaling)
                QSCALING = QAA(KMIE,IDXAER)/QAA(4,IDXAER)
                LOCALOD = QSCALING*AERX_COL(M,L)
                LOCALSSA = SAA(KMIE,IDXAER)*LOCALOD
                OD(KMIE,L) = OD(KMIE,L) + LOCALOD
                SSA(KMIE,L)= SSA(KMIE,L) + LOCALSSA
                DO I=1,8
                   SLEG(I,KMIE,L) = SLEG(I,KMIE,L) + &
                                    (PAA(I,KMIE,IDXAER)*LOCALSSA)
                ENDDO ! I (Phase function)
             ENDIF
          ENDDO ! M (Aerosol)
          !transpose wavelength indices from the mie LUT
          !to the new speciated LUT
          KMIE2=LUTIDX(KMIE)

          IF ( Input_Opt%LUCX ) THEN

             ! Strat aerosols for UCX simulations
             IM=10+(NRHAER*NRH)+1
             DO M=IM,IM+1
                IDXAER=M-IM+6 !6-STS, 7-NAT

                IF (AERX_COL(M,L).gt.0d0) THEN
                   IF (AOD999) THEN
                      ! Aerosol/dust (999 nm scaling)
                      ! Fixed to dry radius
                      QSCALING = QQAA(KMIE2,1,IDXAER)/QQAA(10,1,IDXAER)
                   ELSE
                      ! Aerosol/dust (550 nm scaling)
                      QSCALING = QQAA(KMIE2,1,IDXAER)/QQAA(5,1,IDXAER)
                   ENDIF
                   LOCALOD    = QSCALING*AERX_COL(M,L)
                   LOCALSSA   = SSAA(KMIE2,1,IDXAER)*LOCALOD
                   OD(KMIE,L) = OD(KMIE,L) + LOCALOD
                   SSA(KMIE,L)= SSA(KMIE,L) + LOCALSSA
                   DO I=1,8
                      SLEG(I,KMIE,L) = SLEG(I,KMIE,L) + &
                                       (PAA(I,KMIE,IDXAER)*LOCALSSA)
                   ENDDO     ! I (Phase function)
                ENDIF

             ENDDO           ! M (Aerosol)

          ENDIF ! LUCX

          ! Mineral dust (from new optics LUT)
          DO M=4,10
             IF (AERX_COL(M,L).gt.0d0) THEN
                IDXAER=NSPAA !dust is last in LUT (6, or 8 if UCX=y)
                IR=M-3
                IF (AOD999) THEN
                   QSCALING = QQAA(KMIE2,IR,IDXAER)/ &
                              QQAA(10,IR,IDXAER) !1000nm in new .dat
                ELSE
                   ! Aerosol/dust (550 nm scaling)
                   QSCALING = QQAA(KMIE2,IR,IDXAER)/ &
                              QQAA(5,IR,IDXAER)  !550nm in new .dat
                ENDIF
                LOCALOD = QSCALING*AERX_COL(M,L)
                LOCALSSA = SSAA(KMIE2,IR,IDXAER)*LOCALOD
                OD(KMIE,L) = OD(KMIE,L) + LOCALOD
                SSA(KMIE,L)= SSA(KMIE,L) + LOCALSSA
                DO I=1,8
                   SLEG(I,KMIE,L) = SLEG(I,KMIE,L) + &
                                    (PHAA(KMIE2,IR,IDXAER,I)*LOCALSSA)
                ENDDO ! I (Phase function)
             ENDIF
          ENDDO ! M (Aerosol)

          ! Other aerosol (from new optics LUT)
          DO M=1,5
             DO IR=1,5
                IDXAER=10+(M-1)*NRH+IR
                IF (AERX_COL(IDXAER,L).gt.0d0) THEN
                   IF (AOD999) THEN
                      QSCALING = QQAA(KMIE2,IR,M)/ &
                                 QQAA(10,IR,M) !1000nm in new .dat
                   ELSE
                      ! Aerosol/dust (550 nm scaling)
                      QSCALING = QQAA(KMIE2,IR,M)/ &
                                 QQAA(5,IR,M)  !550nm in new .dat
                   ENDIF
                   LOCALOD = QSCALING*AERX_COL(IDXAER,L)
                   LOCALSSA = SSAA(KMIE2,IR,M)*LOCALOD
                   OD(KMIE,L) = OD(KMIE,L) + LOCALOD
                   SSA(KMIE,L)= SSA(KMIE,L) + LOCALSSA
                   DO I=1,8
                      SLEG(I,KMIE,L) = SLEG(I,KMIE,L) + &
                                       (PHAA(KMIE2,IR,M,I)*LOCALSSA)
                   ENDDO ! I (Phase function)
                ENDIF
             ENDDO    ! IR (RH bins)
          ENDDO ! M (Aerosol)
       ENDDO ! KMIE (Mie scattering wavelength bin)

       ! Normalize
       DO KMIE=1,5
          IF (OD(KMIE,L).gt.0.e+0_fp) THEN
             SSA(KMIE,L) = SSA(KMIE,L)/OD(KMIE,L)
             DO I=1,8
                SLEG(I,KMIE,L) = SLEG(I,KMIE,L)/OD(KMIE,L)
             ENDDO
          ENDIF
       ENDDO
       ! Retrieve 600 nm OD to determine added layers
       OD600(L) = OD(4,L)
    ENDDO ! L (Layer)

    ! when combining with Rayleigh and O2-O3 abs, remember the SSA and
    !  phase fn SLEG are weighted by OD and OD*SSA, respectively.
    ! Given the aerosol+cloud OD/layer in visible (600 nm) calculate how to add
    !  additonal levels at top of clouds (now uses log spacing)
    ! --------------------------------------------------------------------
    call EXTRAL(Input_Opt,OD600,L1_,L2EDGE,N_,JXTRA,ILON,ILAT)
    ! --------------------------------------------------------------------

    ! set surface reflectance
    RFL(:) = max(0.e+0_fp,min(1.e+0_fp,RFLECT))

    ! --------------------------------------------------------------------
    ! Loop over all wavelength bins to calc mean actinic flux AVGF(L)
    ! --------------------------------------------------------------------
    do K = 1,W_

       WAVE = WL(K)
       ! Pick nearest Mie wavelength to get scattering properites------------
       KMIE=1                             ! use 200 nm prop for <255 nm
       if( WAVE .gt. 255.e+0_fp ) KMIE=2  ! use 300 nm prop for 255-355 nm
       if( WAVE .gt. 355.e+0_fp ) KMIE=3  ! use 400 nm prop for 355-500 nm
       if( WAVE .gt. 500.e+0_fp ) KMIE=4
       if( WAVE .gt. 800.e+0_fp ) KMIE=5

       ! Combine: Rayleigh scatters & O2 & O3 absorbers to get optical
       ! properties values at L1_=L_+1 are a pseudo/climatol layer above
       ! the top CTM layer (L_)
       do L = 1,L1_
          TTTX     = T_CLIM(L) ! Following GEOS-Chem v9-1-3
          call X_interp (TTTX,XQO2, TQQ(1,1),QO2(K,1), TQQ(2,1),QO2(K,2), &
                         TQQ(3,1),QO2(K,3), LQQ(1))
          call X_interp (TTTX,XQO3, TQQ(1,2),QO3(K,1), TQQ(2,2),QO3(K,2), &
                         TQQ(3,2),QO3(K,3), LQQ(2))
          ODABS = XQO3*OOJ(L) + XQO2*DDJ(L)*0.20948e+0_fp
          ODRAY = DDJ(L)*QRAYL(K)

          DTAUX(L,K) = OD(KMIE,L) + ODABS + ODRAY

          ! Aerosols + clouds + O2 + O3
          do I=1,8
             POMEGAX(I,L,K) = SLEG(I,KMIE,L)*OD(KMIE,L)
          enddo
          ! Add Rayleigh scattering effects
          ! Only non-zero for 1st and 3rd phase functions
          POMEGAX(1,L,K) = POMEGAX(1,L,K) + 1.0e+0_fp*ODRAY
          POMEGAX(3,L,K) = POMEGAX(3,L,K) + 0.5e+0_fp*ODRAY
          ! Normalize
          do I=1,8
             POMEGAX(I,L,K) = POMEGAX(I,L,K)/DTAUX(L,K)
          enddo
       enddo
    enddo
    ! --------------------------------------------------------------------
    call OPMIE(DTAUX,POMEGAX,U0,RFL,AMF2,JXTRA, &
               AVGF,FJTOP,FJBOT,FSBOT,FJFLX,FLXD,FLXD0,State_Grid%NZ)

    !! --------------------------------------------------------------------

    do K = 1,W_

       do L = 1,L_
          FFF(K,L) = FFF(K,L) + SOLF*FL(K)*AVGF(L,K)
       enddo

    enddo

    ! Calculate photolysis rates
    call JRATET(PPJ,T_INPUT,FFF, VALJXX,L_,State_Grid%MaxChemLev,NJX)

    ! Fill out common-block array of J-rates
    DO L=1,State_Grid%MaxChemLev
       DO J=1,NRATJ
          IF (JIND(J).gt.0) THEN
             ZPJ(L,J,ILON,ILAT) = VALJXX(L,JIND(J))*JFACTA(J)
          ELSE
             ZPJ(L,J,ILON,ILAT) = 0.e+0_fp
          ENDIF
       ENDDO
    ENDDO

    ! Set J-rates outside the chemgrid to zero
    IF (State_Grid%MaxChemLev.lt.L_) THEN
       DO L=State_Grid%MaxChemLev+1,L_
          DO J=1,NRATJ
             ZPJ(L,J,ILON,ILAT) = 0.e+0_fp
          ENDDO
       ENDDO
    ENDIF

    !=================================================================
    ! UV radiative fluxes (direct, diffuse, net) [W/m2]
    !
    ! Updated for netCDF from nd64 (JMM 2019-09-11)
    ! Use it to calculate fluxes for output if necessary
    !
    ! Get net direct and net diffuse fluxes separately
    ! Order:
    !    1 - Net flux
    !    2 - Direct flux
    !    3 - Diffuse flux
    ! Convention: negative is downwards
    !=================================================================
    IF ( State_Diag%Archive_UVFluxDiffuse .or. &
         State_Diag%Archive_UVFluxDirect .or. &
         State_Diag%Archive_UVFluxNet ) THEN

       ! Loop over wavelength bins
       DO K = 1, W_

          ! Initialize
          FDIRECT  = 0.0_fp
          FDIFFUSE = 0.0_fp

          ! Direct & diffuse fluxes at each level
          FDIRECT(1)  = FSBOT(K)                    ! surface
          FDIFFUSE(1) = FJBOT(K)                    ! surface
          DO L = 2, State_Grid%NZ
             FDIRECT(L) = FDIRECT(L-1) + FLXD(L-1,K)
             FDIFFUSE(L) = FJFLX(L-1,K)
          ENDDO

          ! Constant to multiply UV fluxes at each wavelength bin
          UVX_CONST = SOLF * FL(K) * UVXFACTOR(K)

          ! Archive into diagnostic arrays
          DO L = 1, State_Grid%NZ

             IF ( State_Diag%Archive_UVFluxNet ) THEN
                State_Diag%UVFluxNet(ILON,ILAT,L,K) = &
                State_Diag%UVFluxNet(ILON,ILAT,L,K) + &
                     ( ( FDIRECT(L) + FDIFFUSE(L) ) * UVX_CONST )
             ENDIF

             IF ( State_Diag%Archive_UVFluxDirect ) THEN
                State_Diag%UVFluxDirect(ILON,ILAT,L,K) = &
                State_Diag%UVFluxDirect(ILON,ILAT,L,K) + &
                     ( FDIRECT(L) * UVX_CONST )
             ENDIF

             IF ( State_Diag%Archive_UVFluxDiffuse ) THEN
                State_Diag%UVFluxDiffuse(ILON,ILAT,L,K) = &
                State_Diag%UVFluxDiffuse(ILON,ILAT,L,K) + &
                     ( FDIFFUSE(L) * UVX_CONST )
             ENDIF
          ENDDO
       ENDDO
    ENDIF

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%% COMMENT OUT FAST-J DIAGNOSTICS FOR NOW (bmy, 3/5/14)
!%%% LEAVE CODE HERE SO THAT IT CAN BE RESTORED
!      !---------------------------------------------------------------
!      ! Majority of heating code ignored by GEOS-Chem
!      ! Everything from here until statement number 99 is ignored
!      ! (SDE 03/31/13)
!      !---------------------------------------------------------------
!      IF (LFJXDIAG) THEN
!         DO K=1,W_
!! -direct(DIR) and diffuse(FLX) fluxes at top(UP) (solar = negative by convention)
!! -     also at bottom (DN), does not include diffuse reflected flux.
!            FLXUP(K) =  FJTOP(K)
!            DIRUP(K) = -FLXD0(K)
!            FLXDN(K) = -FJBOT(K)
!            DIRDN(K) = -FSBOT(K)
!
!            FREFI = FREFI + SOLF*FL(K)*FLXD0(K)/WL(K)
!            FREFL = FREFL + SOLF*FL(K)*FJTOP(K)/WL(K)
!            FREFS = FREFS + SOLF*FL(K)/WL(K)
!
!! for each wavelength calculate the flux budget/heating rates:
!c  FLXD(L) = direct flux deposited in layer L  [approx = MU0*(F(L+1) -F(L))]
!c            but for spherical atmosphere!
!c  FJFLX(L) = diffuse flux across top of layer L
!
!! calculate divergence of diffuse flux in each CTM layer (& t-o-a)
!!      need special fix at top and bottom:
!! FABOT = total abs at L.B. &  FXBOT = net diffusive flux at L.B.
!            FABOT = (1.e+0_fp-RFL(K))*(FJBOT(K)+FSBOT(K))
!            FXBOT = -FJBOT(K) + RFL(K)*(FJBOT(K)+FSBOT(K))
!            FLXJ(1) = FJFLX(1,K) - FXBOT
!            do L=2,LU
!               FLXJ(L) = FJFLX(L,K) - FJFLX(L-1,K)
!            enddo
!            FLXJ(LU+1) = FJTOP(K) - FJFLX(LU,K)
!! calculate net flux deposited in each CTM layer (direct & diffuse):
!            FFX0 = 0.e+0_fp
!            do L=1,L1U
!               FFX(K,L) = FLXD(L,K) - FLXJ(L)
!               FFX0 = FFX0 + FFX(K,L)
!            enddo
!
!c  NB: the radiation level ABOVE the top CTM level is included in these budgets
!c      these are the flux budget/heating terms for the column:
!c  FFXNET(K,1) = FLXD0        direct(solar) flux dep into atmos (spherical)
!c  FFXNET(K,2) = FSBOT        direct(solar) flux dep onto LB (surface)
!c  FFXNET(K,3) = FLXD0+FSBOT  TOTAL solar into atmopshere+surface
!c  FFXNET(K,4) = FJTOP        diffuse flux leaving top-of-atmos
!c  FFXNET(K,5) = FFX0         diffuse flux absorbed in atmos
!c  FFXNET(K,6) = FABOT        total (dir+dif) absorbed at LB (surface)
!c       these are surface fluxes to compare direct vs. diffuse:
!c  FFXNET(K,7) = FSBOT        direct flux dep onto LB (surface) - for srf diags
!c  FFXNET(K,8) = FJBOT        diffuse flux dep onto LB (surface)
!
!            FFXNET(K,1) = FLXD0(K)
!            FFXNET(K,2) = FSBOT(K)
!            FFXNET(K,3) = FLXD0(K) + FSBOT(K)
!            FFXNET(K,4) = FJTOP(K)
!            FFXNET(K,5) = FFX0
!            FFXNET(K,6) = FABOT
!            FFXNET(K,7) = FSBOT(K)
!            FFXNET(K,8) = FJBOT(K)
!
!! --------------------------------------------------------------------
!         enddo       ! end loop over wavelength K
!! --------------------------------------------------------------------
!         FREFL = FREFL/FREFS      !calculate reflected flux (energy weighted)
!         FREFI = FREFI/FREFS
!
!! NB UVB = 280-320 = bins 12:15, UVA = 320-400 = bins 16:17, VIS = bin 18 (++)
!
!
!! mapping J-values from fast-JX onto CTM chemistry is done in main
!
!! --------------------------------------------------------------------
!         if ((Input_Opt%amIRoot).and.(LPRTJ)) then
!! diagnostics below are NOT returned to the CTM code
!            write(6,*)'fast-JX-(7.0)---PHOTO_JX internal print:',
!     &               ' Atmosphere---'
!! used last called values of DTAUX and POMEGAX, should be 600 nm
!         do L=1,L1U
!            DTAU600(L) = DTAUX(L,W_)
!            do I=1,8
!               POMG600(I,L) = POMEGAX(I,L,W_)
!            enddo
!         enddo
!
!      !   call JP_ATM(PPJ,TTJ,DDJ,OOJ,ZZJ,DTAU600,POMG600,JXTRA, LU)
!
!! PRINT SUMMARY of mean intensity, flux, heating rates:
!         if (Input_Opt%amIRoot) then
!            write(6,*)
!            write(6,*)'fast-JX(7.0)---PHOTO_JX internal print:',
!     &               ' Mean Intens---'
!            write(6,'(a,5f10.4)')
!     &       ' SUMMARY fast-JX: albedo/SZA/u0/F-incd/F-refl/',
!     &        RFLECT,SZA,U0,FREFI,FREFL
!
!            write(6,'(a5,18i8)')   ' bin:',(K, K=NW2,NW1,-1)
!            write(6,'(a5,18f8.1)') ' wvl:',(WL(K), K=NW2,NW1,-1)
!            write(6,'(a,a)') ' ----  100000=Fsolar   ',
!     &                  'MEAN INTENSITY per wvl bin'
!         endif
!         do L = LU,1,-1
!            do K=NW1,NW2
!               RATIO(K) = (1.d5*FFF(K,L)/FL(K))
!            enddo
!            if (Input_Opt%amIRoot) then
!               write(6,'(i3,2x,18i8)') L,(RATIO(K),K=NW2,NW1,-1)
!            endif
!         enddo
!
!         if (Input_Opt%amIRoot) then
!            write(6,*)
!            write(6,*)'fast-JX(7.0)---PHOTO_JX internal print:',
!                              ' Net Fluxes---'
!            write(6,'(a11,18i8)')   ' bin:',(K, K=NW2,NW1,-1)
!            write(6,'(a11,18f8.1)') ' wvl:',(WL(K), K=NW2,NW1,-1)
!c            write(6,'(a11,18f8.4)') ' sol in atm',(FFXNET(K,1),
!c     &                        K=NW2,NW1,-1)
!c            write(6,'(a11,18f8.4)') ' sol at srf',(FFXNET(K,2),
!c     &                        K=NW2,NW1,-1)
!            write(6,*) ' ---NET FLUXES--- '
!            write(6,'(a11,18f8.4)') ' sol TOTAL ',(FFXNET(K,3),
!     &                        K=NW2,NW1,-1)
!            write(6,'(a11,18f8.4)') ' dif outtop',(FFXNET(K,4),
!     &                        K=NW2,NW1,-1)
!            write(6,'(a11,18f8.4)') ' abs in atm',(FFXNET(K,5),
!     &                        K=NW2,NW1,-1)
!            write(6,'(a11,18f8.4)') ' abs at srf',(FFXNET(K,6),
!     &                        K=NW2,NW1,-1)
!            write(6,*) ' ---SRF FLUXES--- '
!            write(6,'(a11,18f8.4)') ' srf direct',(FFXNET(K,7),
!     &                        K=NW2,NW1,-1)
!            write(6,'(a11,18f8.4)') ' srf diffus',(FFXNET(K,8),
!     &                        K=NW2,NW1,-1)
!            write(6,'(4a)') '  ---NET ABS per layer:',
!     &         '       10000=Fsolar',
!     &         '  [NB: values <0 = numerical error w/clouds',
!     &         ' or SZA>90, colm OK]'
!         endif
!         do L = LU,1,-1
!            do K=NW1,NW2
!               RATIO(K) = 1.d5*FFX(K,L)
!            enddo
!            if (Input_Opt%amIRoot) then
!               write(6,'(i9,2x,18i8)') L,(RATIO(K),K=NW2,NW1,-1)
!            endif
!         enddo
!         if (Input_Opt%amIRoot) then
!            write(6,'(a)')
!            write(6,'(a)') ' fast-JX (7.0)----J-values----'
!            write(6,'(1x,a,72(a6,3x))') 'L=  ',(TITLEJX(K), K=1,NJX)
!            do L = LU,1,-1
!              write(6,'(i3,1p, 72e9.2)') L,(VALJXX(L,K),K=1,NJX)
!            enddo
!         endif
!
!      ENDIF
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

99  continue

  END SUBROUTINE PHOTO_JX
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: opmie
!
! !DESCRIPTION: Subroutine OPMIE is a core Fast-JX scattering subroutine,
!  specifically for Mie scattering.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE OPMIE(DTAUX,POMEGAX,U0,RFL,AMF2,JXTRA, &
                   FJACT,FJTOP,FJBOT,FSBOT,FJFLX,FLXD,FLXD0,LU)
!
! !USES:
!
!
! !INPUT PARAMETERS:
!
    REAL(fp), INTENT(IN)  ::  DTAUX(JXL1_,W_),POMEGAX(8,JXL1_,W_)
    REAL(fp), INTENT(IN)  ::  AMF2(2*JXL1_+1,2*JXL1_+1)
    REAL(fp), INTENT(IN)  ::  U0,RFL(W_)
    INTEGER,  INTENT(IN)  ::  JXTRA(JXL2_+1), LU
!
! !OUTPUT VARIABLES:
!
    REAL(fp), INTENT(OUT) ::  FJACT(JXL_,W_),FJTOP(W_)
    REAL(fp), INTENT(OUT) ::  FJBOT(W_),FSBOT(W_)
    REAL(fp), INTENT(OUT) ::  FJFLX(JXL_,W_),FLXD(JXL1_,W_),FLXD0(W_)
!
! !REVISION HISTORY:
!  28 Mar 2013 - S. D. Eastham - Copied from Fast-JX v7.0
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER JNDLEV(JXL_),JNELEV(JXL1_)
    INTEGER JADDLV(JXL2_+1),JADDTO(JXL2_+1),L2LEV(JXL2_+1)
    INTEGER JTOTL,I,II,J,K,L,LL,IX,JK,   L2,L2L,L22,LZ,LZZ,ND
    INTEGER L1U,L2U,   LZ0,LZ1,LZMID
    REAL(fp)   SUMT,SUMJ

    REAL(fp)  DTAU(JXL1_+1,W_),POMEGAJ(M2_,JXL2_+1,W_)
    REAL(fp)  TTAU(JXL2_+1,W_)
    REAL(fp)  FTAU2(JXL2_+1,W_),POMEGAB(M2_,W_)
    REAL(fp)  ATAUA,ATAUZ,XLTAU,TAUDN,TAUUP,DTAUJ,FJFLX0
    REAL(fp), DIMENSION(W_) :: TAUBTM,TAUTOP,FBTM,FTOP,ZFLUX

    !  variables used in mie code-----------------------------------------
    REAL(fp), DIMENSION(W_)         :: FJT,FJB
    REAL(fp), DIMENSION(N_,W_)      :: FJ,FZ,ZTAU
    REAL(fp), DIMENSION(M2_,N_,W_)  :: POMEGA
    REAL(fp), DIMENSION(2*JXL1_,W_)  :: FLXD2

    !=================================================================
    ! OPMIE begins here!
    !=================================================================

    ! there is a parallel correspondence:
    !  dimension of JX arrays JXL_ .ge. dimension that CTM is using = L_
    !  but calculation is done for L_=LU, L1_=L1U, L2_=L2U lengths of CTM
    !
    !  fast-J Mie code for J_s, only uses 8-term expansion, 4-Gauss pts
    !
    ! in:
    !     DTAUX(1:L1_,1:W_) = optical depth of each layer
    !     POMEGAX(1:8,1:L1_,1:W_) = scattering phase fn (multiplied by s-s albedo)
    !     U0  = cos (SZA)
    !     RFL(1:W_) = Lambertian albedo of surface
    !     AMF2(1:2*L1_+1,1:2*L1_+1) = air mass factor (I,L)=wt of layer-I to layer-L
    !        AMF2 now does both edges and middle of CTM layers
    !     JXTRA(1:L1_) = number 0:J = no. of additional levels to be inserted
    ! out:
    !     FJACT(1:L_,1:W_) = mean actinic flux(diff+direct) at std CTM levels(mid-lyr)
    !  (new ver 5.7 diagnostics for fluxes, deposition)  fluxes 'down' are <0
    !     FJTOP(1:W_) = diffuse flux out top-of-atmosphere (TAU=0 above top model lyr)
    !     FJBOT(1:W_) = diffuse flux onto surface (<0 by definition)
    !     FSBOT(1:W_) = direct/solar flux onto surface  (<0 by definition)
    !     FJFLX(1:L_,1:W_) = diffuse flux across top of model layer L
    !        this connects with FJBOT = FJFLX(0) & FJTOP = FJFLX(L_+1) (not dim!!)
    !     FLXD(1:L_+1,1:W_) = solar flux deposited in layer L (includes lyr above CTM)
    !        this should take into account sphericity, and is not just = mu0
    !     FLXD0(1:W_) = sum of solar flux deposited in atmos
    !        does NOT include flux on lower surface, does NOT mean absorbed!
    ! --------------------------------------------------------------------
    !
    !     DTAU     Local optical depth of each CTM level
    !     TTAU     Optical depth of air vertically above each point (to top of atm)
    !     FTAU2     Attenuation of solar beam
    !     POMEGAJ  Scattering phase function
    !
    ! new ver 5.3 code adds sub-layers (# = JXTRA(L2)) using ATAU as the
    !   factor increase from sub-layer to sub-layer
    !
    ! ------------------SET UP FOR MIE CODE-------------------------------
    !
    ! --------------wavelength independent--------------------------------
    !
    !  Transpose the ascending TTAU grid to a descending ZTAU grid.
    !  Double the resolution - TTAU points become the odd points on the
    !  ZTAU grid, even points needed for asymm phase fn soln, contain 'h'.
    !  Odd point added at top of grid for unattenuated beam   (Z='inf')
    !
    !  The following mapping holds for JADDLV=0
    !        Surface:   TTAU(1)    ==> ZTAU(2*L2_+1)
    !        Top:       TTAU(L2_)  ==> ZTAU(3)
    !        Infinity:     0.0     ==> ZTAU(1)
    !        index: 2*(L2_+1-L2)+1 ==> LZ
    !
    !  Mie scattering code only used from surface to level L2_
    ! ---------------------------------------------------------------------
    !
    ! ---------------------------------------------------------------------
    !  Insert new levels, working downwards from the top of the atmosphere
    !  to the surface (down in 'LZ', up in 'L2'). This allows ztau and pomega
    !  to be incremented linearly, and the flux fz to be attenuated top-down
    !    (avoiding problems where lower level fluxes are zero).
    ! ---------------------------------------------------------------------
    !
    !  Ascend through atmosphere transposing grid and adding extra points
    !  remember L2=1 is surface of CTM, but last layer (LZ) in scattering code.
    !  there are twice the number of layers in the LZ arrays (2*L2_ + 2*JADDTO + 1)
    !    because we need to insert the intermediate layers (even LZ) for the
    !    asymmetric scattering code.
    !
    !  Transfer the L2=1:L2_+1 values (TTAU,FTAU2,POMEGAJ) onto the reverse
    !    order, expanded, doubled-level scatter grid.
    !    Note that we need to deal with the expansion by JADD levels (L2L).
    !      These JADDLV levels are skipped and need to be interpolated later.
    !    Note that only odd LZ levels are filled,
    !
    ! -------------------re-grid data---------------------------------------------
    !  Calculate cumulative total and define levels we want J-values at.
    !  Sum upwards for levels, and then downwards for Mie code readjustments.
    !
    !     JXTRA(L2)  Number of new levels to add between (L2) and (L2+1)
    !           ***JXTRA(1:L2_+1) is calculated based on the aerosol+cloud OD_s
    !     JADDLV(L2)  Number of new levels actually added at each wavelength
    !            where JADDLV = 0 when there is effectively no FTAU2
    !     JADDTO(L2)   Total number of new levels to add to and above level (L2)
    !     JNDLEV(L) = L2 index that maps on CTM mid-layer L
    !
    ! JADDLV(L2=1:L2_) = number of levels to add between TTAU2(L2) and TTAU(L2+1)
    !     JADDLV is taken from JXTRA, which is based on visible OD.
    !     JADDTO(L2=1:L2_+1) is the cumulative number of levels to be added
    ! these should be fixed for all wavelengths to lock-in the array sizes

    if (LU .gt. JXL_) then
       call EXITC (' OPMIE:  JXL_ .lt. L_')
    endif

    L1U = LU + 1
    L2U = 2*LU + 2

    do L2 = 1,L2U,1
       JADDLV(L2) = JXTRA(L2)
    enddo
    JADDTO(L2U+1) = 0
    do L2 = L2U,1,-1
       JADDTO(L2) = JADDTO(L2+1) + JADDLV(L2)
    enddo

    ! expanded grid now included CTM edge and mid layers plus expanded
    !     grid to allow for finer delta-tau at tops of clouds.
    !     DIM of new grid = L2U + JADDTO(1) + 1

    ! L2LEV(L2) = L2-index for old level L2 in expanded J-grid (w/JADDLV)
    !     in absence of JADDLV, L2LEV(L2) = L2
    L2LEV(1)  = 1
    do L2 = 2,L2U+1
       L2LEV(L2) = L2LEV(L2-1) + 1 + JADDLV(L2-1)
    enddo

    ! JNDLEV(L=1:L_) = L2-index in expanded grid for CTM mid-layer L
    ! JNELEV(L=1:L_) = L2-index for top of layer L
    do L = 1,LU
       JNDLEV(L) = L2LEV(2*L)
       JNELEV(L) = L2LEV(2*L+1)
    enddo
    JNELEV(LU+1) = 0  !need to set this to top-of-atmosphere

    ND = 2*L2U + 2*JADDTO(1) + 1

    if(ND .gt. N_) then
       call EXITC (' overflow of scatter arrays: ND > N_')
    endif

    ! -------------begin wavelength dependent set up---------------------------

    ! Reinitialize arrays
    ZTAU(:,:)     = 0.e+0_fp
    FZ(:,:)       = 0.e+0_fp
    POMEGA(:,:,:) = 0.e+0_fp

    do K=1,W_

       ! Set up optical depth DTAU(L)
       do L = 1,L1U
          DTAU(L,K) = DTAUX(L,K)
       enddo
       DTAU(L1U+1,K) = 0.e+0_fp

       ! Define the total scattering phase fn for each CTM layer L=1:L_+1
       !    from a DTAU-wt_d mix of aerosols, cloud & Rayleigh
       ! No. of quadrature pts fixed at 4(M_), expansion of phase fn @ 8
       do L = 1,L1U
          do I = 1,M2_
             POMEGAJ(I,L,K) = POMEGAX(I,L,K)
          enddo
       enddo

       ! Calculate attenuated incident beam exp(-TTAU/U0 = DTAU * AirMassFactor)
       !       at the middle & edges of the CTM layers L=1:2*L1_+1
       !   L1_ is top-edge of CTM (ie, L=38 = 2 hPa) which has TAU > 0
       !   note that DTAU(L1_) is optical depth in the FULL CTM layer just above
       FTAU2(:,:) = 0.e+0_fp
       FTAU2(L2U+1,:) = 1.0e+0_fp
       do LL = 1,2*L1U+1
          L = (LL+1)/2
          if (AMF2(LL,LL) .gt. 0.0e+0_fp) then
             XLTAU = 0.0e+0_fp
             do II = 1,2*L1U+1
                I = (II+1)/2
                XLTAU = XLTAU + 0.5e+0_fp*DTAU(I,K)*AMF2(II,LL)
             enddo
             if (XLTAU .lt. 76.e+0_fp) then   ! zero out flux at 1e-33
                FTAU2(LL,K) = exp(-XLTAU)
             endif
          endif
       enddo

       ! calculate direct solar flux deposited in each CTM half-layer: L=1:L2_
       !      use FSBOT for surface flux, cannot do layer above CTM (L_+1)
       FLXD2(:,:) = 0.e+0_fp
       do LL = 1,2*L1U
          if (AMF2(LL,LL) .gt. 0.e+0_fp) then
             FLXD2(LL,K) = (FTAU2(LL+1,K) - FTAU2(LL,K))/AMF2(LL,LL)
          endif
       enddo
       if (AMF2(1,1) .gt. 0.e+0_fp) then
          FSBOT(K) = FTAU2(1,K)/AMF2(1,1)
       else
          FSBOT(K) = 0.e+0_fp
       endif

       do LL = 2,2*L1U,2
          L=LL/2
          FLXD(L,K) = FLXD2(LL,K)+FLXD2(LL-1,K)
       enddo

       ! integrate solar flux depositied in CTM layers L=1:L_, cannot do top layer
       !   note FLXD0 .ne. (1.e+0_fp - FTAU(L_+1))/AMF(L_+1,L_+1) with spherical atmos
       FLXD0(K) = 0.e+0_fp
       if (AMF2(2*L1U,2*L1U) .gt. 0.e+0_fp) then
          do L=1,L1U
             FLXD0(K) = FLXD0(K) + FLXD(L,K)
          enddo
       endif

       ! ---------------------------------------------------------------------
       !  Take optical properties on CTM layers and convert to a photolysis
       !  level grid corresponding to layer centres and boundaries. This is
       !  required so that J-values can be calculated for the centre of CTM
       !  layers; the index of these layers is kept in the JNDLEV array.
       ! ---------------------------------------------------------------------
       ! Now combine the CTM layer edges (1:L_+2) with the CTM mid-layer
       !     points (1:L_) plus 1 for the mid point of added top layer.
       ! combine these edge- and mid-layer points into grid of size:
       !               L2_+1 = 2*L1_+1 = 2*L_+3
       ! calculate column optical depths above each level, TTAU(1:L2_+1)
       !       note that TTAU(L2_+1)=0 and TTAU(1)=total OD

       TTAU(L2U+1,K) = 0.0e+0_fp
       do L2 = L2U,1,-1
          L          = (L2+1)/2
          DTAUJ      = 0.5e+0_fp * DTAU(L,K)
          TTAU(L2,K)   = TTAU(L2+1,K) + DTAUJ
       enddo

       ! -solar flux incident on lower boundary & Lambertian reflect factor:
       if (FSBOT(K) .gt. 0.e+0_fp) then
          ZFLUX(K) = FSBOT(K)*RFL(K)/(1.e+0_fp+RFL(K))
       else
          ZFLUX(K) = 0.e+0_fp
       endif

       !  Calculate scattering properties, level centres then level boundaries
       !>>>>>be careful of order, we are overwriting/shifting the 'POMEGAJ'
       !     upward in index
       do L2 = L2U,2,-2
          L   = L2/2
          do I = 1,M2_
             POMEGAJ(I,L2,K) = POMEGAJ(I,L,K)
          enddo
       enddo
       ! lower boundary value is set (POMEGAJ(I,1)), but set upper:
       do I = 1,M2_
          POMEGAJ(I,L2U+1,K) = POMEGAJ(I,L2U,K)
       enddo
       ! now have POMEGAJ filled at even points from L2=3:L2_-1
       ! use inverse interpolation for correct tau-weighted values at edges
       do L2 = 3,L2U-1,2
          TAUDN = TTAU(L2-1,K)-TTAU(L2,K)
          TAUUP = TTAU(L2,K)-TTAU(L2+1,K)
          do I = 1,M2_
             POMEGAJ(I,L2,K) = (POMEGAJ(I,L2-1,K)*TAUDN + &
                                POMEGAJ(I,L2+1,K)*TAUUP) / (TAUDN+TAUUP)
          enddo
       enddo

       ! at this point FTAU2(1:L2_+1) and POMEAGJ(1:8, 1:L2_+1)
       !     where FTAU2(L2_+1) = 1.0 = top-of-atmos, FTAU2(1) = surface
       
       do L2 = 1,L2U+1          ! L2 = index of CTM edge- and mid-layers
          L2L = L2LEV(L2)        ! L2L = index for L2 in expanded scale(JADD)
          LZ  = ND + 2 - 2*L2L  ! LZ = index for L2 in scatt arrays
          ZTAU(LZ,K) = TTAU(L2,K)
          FZ(LZ,K)   = FTAU2(L2,K)
          do I=1,M2_
             POMEGA(I,LZ,K) = POMEGAJ(I,L2,K)
          enddo
       enddo

       ! Now go thru the pairs of L2 levels to see if we need JADD levels
       do L2 = 1,L2U             ! L2 = index of CTM edge- and mid-layers
          L2L = L2LEV(L2)        ! L2L = index for L2 in expanded scale(JADD)
          LZ  = ND + 2 - 2*L2L   ! LZ = index for L2 in scatt arrays
          L22 = L2LEV(L2+1) - L2LEV(L2) - 1   ! L22 = 0 if no added levels

          if (L22 .gt. 0) then
             TAUBTM(K) = TTAU(L2,K)
             TAUTOP(K) = TTAU(L2+1,K)
             FBTM(K)   = FTAU2(L2,K)
             FTOP(K)   = FTAU2(L2+1,K)
             do I = 1,M2_
                POMEGAB(I,K) = POMEGAJ(I,L2,K)
             enddo

             ! to fit L22 new layers between TAUBOT > TAUTOP, calculate new
             ! 1/ATAU factor
             ! such that TAU(just above TAU-btm) = ATUAZ * TAUBTM < TAUBTM
             ATAUZ = exp(-log(TAUBTM(K)/max(TAUTOP(K),ATAU0))/float(L22+1))
             do L = 1,L22      ! add odd levels between L2LEV(L2) & L2LEV(L2+1)
                LZZ = LZ - 2*L ! LZZ = index(odd) of added level in scatt arrays
                ZTAU(LZZ,K) = TAUBTM(K) * ATAUZ

                ! fraction from TAUBTM=>TAUTOP
                ATAUA=(TAUBTM(K)-ZTAU(LZZ,K))/(TAUBTM(K)-TAUTOP(K))
                ! solar flux at interp-levels: use exp(TAU/U0) if U0>0.02
                ! (89 deg), else scale by TAU
                if (U0 .gt. 0.02e+0_fp) then
                   FZ(LZZ,K) = FTOP(K) * exp((TAUTOP(K)-ZTAU(LZZ,K))/U0)
                else
                   if (FBTM(K) .lt. 1.d-32) then
                      FZ(LZZ,K) = 0.e+0_fp
                   else
                      FZ(LZZ,K) = FBTM(K) * (FTOP(K)/FBTM(K))**ATAUA
                   endif
                endif
                do I = 1,M2_
                   POMEGA(I,LZZ,K) = POMEGAB(I,K) + &
                        ATAUA*(POMEGAJ(I,L2+1,K)-POMEGAB(I,K))
                enddo
                TAUBTM(K)    = ZTAU(LZZ,K)
                FBTM(K)      = FZ(LZZ,K)
                do I = 1,M2_
                   POMEGAB(I,K) = POMEGA(I,LZZ,K)
                enddo
             enddo
          endif
       enddo

       ! Now fill in the even points with simple interpolation in scatter arrays
       do LZ = 2,ND-1,2
          ZTAU(LZ,K) = 0.5e+0_fp*(ZTAU(LZ-1,K)+ZTAU(LZ+1,K))
          FZ(LZ,K)   = sqrt(FZ(LZ-1,K)*FZ(LZ+1,K))
          do I=1,M2_
             POMEGA(I,LZ,K) = 0.5e+0_fp*(POMEGA(I,LZ-1,K)+POMEGA(I,LZ+1,K))
          enddo
       enddo

    enddo  ! wavelength loop!

    ! --------------------------------------------------------------------
    call MIESCT(FJ,FJT,FJB,POMEGA,FZ,ZTAU,ZFLUX,RFL,U0,ND)
    ! --------------------------------------------------------------------

    ! Move mean intensity from scatter array FJ(LZ=1:ND)
    !               to CTM mid-level array FJACT(L=1:L_)

    do K=1,W_

       ! mean intensity at mid-layer:  4*<I> + solar
       !do L = 1,LU
       ! L2L = JNDLEV(L)
       ! LZ  = ND+2 - 2*L2L
       ! FJACT(L,K) = 4.e+0_fp*FJ(LZ,K) + FZ(LZ,K)
       !enddo

       ! mean intensity averaged throughout layer:
       do L = 1,LU
          LZ0 = ND+2 - 2*JNELEV(L)
          if (L .gt. 1) then
             LZ1 = ND+2 - 2*JNELEV(L-1)
          else
             LZ1 = ND
          endif
          SUMJ = (4.e+0_fp*FJ(LZ0,K)+FZ(LZ0,K))*(ZTAU(LZ0+2,K)-ZTAU(LZ0,K)) &
               + (4.e+0_fp*FJ(LZ1,K)+FZ(LZ1,K))*(ZTAU(LZ1,K)-ZTAU(LZ1-2,K))
          SUMT = ZTAU(LZ0+2,K)-ZTAU(LZ0,K) + ZTAU(LZ1,K)-ZTAU(LZ1-2,K)

          do LZ = LZ0+2,LZ1-2,2
             SUMJ =SUMJ+(4.e+0_fp*FJ(LZ,K)+FZ(LZ,K))*(ZTAU(LZ+2,K)-ZTAU(LZ-2,K))
             SUMT =SUMT + ZTAU(LZ+2,K)-ZTAU(LZ-2,K)
          enddo
          FJACT(L,K) = SUMJ/SUMT

       enddo

       ! mean diffuse flux:  4<I*mu> (not solar) at top of layer L
       !       average (tau-wtd) the h's just above and below the L-edge
       do L = 1,LU
          L2L = JNELEV(L)
          LZ  = ND+2 - 2*L2L
          FJFLX0 = (ZTAU(LZ+1,K)-ZTAU(LZ,K))/(ZTAU(LZ+1,K)-ZTAU(LZ-1,K))
          FJFLX(L,K)=4.e+0_fp*(FJ(LZ-1,K)*FJFLX0 +FJ(LZ+1,K)*(1.e+0_fp-FJFLX0))
       enddo

       ! diffuse fluxes reflected at top, incident at bottom
       FJTOP(K) = FJT(K)
       FJBOT(K) = FJB(K)

    enddo  ! wavelength loop!

  END SUBROUTINE OPMIE
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: miesct
!
! !DESCRIPTION: Subroutine MIESCT is an adaptation of the Prather radiative
!  transfer code (mjp, 10/95).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE MIESCT(FJ,FJT,FJB, POMEGA,FZ,ZTAU,ZFLUX,RFL,U0,ND)
!
! !USES:
!
!
! !INPUT PARAMETERS:
!
    INTEGER,  INTENT(IN)  ::  ND
    REAL(fp), INTENT(IN)  ::  POMEGA(M2_,N_,W_),FZ(N_,W_),ZTAU(N_,W_), &
                              RFL(W_),U0,ZFLUX(W_)
!
! !OUTPUT VARIABLES:
!
    REAL(fp), INTENT(OUT) ::  FJ(N_,W_),FJT(W_),FJB(W_)
!
! !REMARKS:
!     Prather, 1974, Astrophys. J. 192, 787-792.
!         Solution of inhomogeneous Rayleigh scattering atmosphere.
!         (original Rayleigh w/ polarization)
!     Cochran and Trafton, 1978, Ap.J., 219, 756-762.
!         Raman scattering in the atmospheres of the major planets.
!         (first use of anisotropic code)
!     Jacob, Gottlieb and Prather, 1989, J.Geophys.Res., 94, 12975-13002.
!         Chemistry of a polluted cloudy boundary layer,
!         (documentation of extension to anisotropic scattering)
!                                                                             .
!    takes atmospheric structure and source terms from std J-code
!    ALSO limited to 4 Gauss points, only calculates mean field! (M=1)
!
! !REVISION HISTORY:
!  28 Mar 2013 - S. D. Eastham - Copied from Fast-JX v7.0
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    REAL(fp)  PM(M_,M2_),PM0(M2_)
    INTEGER I, IM  ,K

    !=================================================================
    ! MIESCT begins here!
    !=================================================================

    do I = 1,M_
       call LEGND0 (EMU(I),PM0,M2_)
       do IM = 1,M2_
          PM(I,IM) = PM0(IM)
       enddo
    enddo

    call LEGND0 (-U0,PM0,M2_)
    do IM=1,M2_
       PM0(IM) = 0.25e+0_fp*PM0(IM)
    enddo

    ! BLKSLV now called with all the wavelength arrays (K=1:W_)
    call BLKSLV(FJ,POMEGA,FZ,ZTAU,ZFLUX,RFL,PM,PM0,FJT,FJB, ND)

  END SUBROUTINE MIESCT
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: legnd0
!
! !DESCRIPTION: Subroutine LEGND0 calculates ordinary Legendre functions
!  of X (real) from $P[0] = PL(1) = 1, P[1] = X, \dots, P[N-1] = PL(N)$
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE LEGND0 (X,PL,N)
!
! !USES:
!
!
! !INPUT PARAMETERS:
!
    INTEGER,  INTENT(IN)  :: N
    REAL(fp), INTENT(IN)  :: X
!
! !OUTPUT VARIABLES:
!
    REAL(fp), INTENT(OUT) :: PL(N)
!
! !REMARKS:
!
! !REVISION HISTORY:
!  28 Mar 2013 - S. D. Eastham - Copied from Fast-JX v7.0
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER I
    REAL(fp)  DEN

    !=================================================================
    ! LEGND0 begins here!
    !=================================================================

    ! Always does PL(2) = P[1]
    PL(1) = 1.e+0_fp
    PL(2) = X
    do I = 3,N
       DEN = (I-1)
       PL(I) = PL(I-1)*X*(2.e+0_fp-1.0/DEN) - PL(I-2)*(1.e+0_fp-1.e+0_fp/DEN)
    enddo

  END SUBROUTINE LEGND0
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: set_prof
!
! !DESCRIPTION: Subroutine SET\_PROF sets vertical profiles for a given
!  latitude and longitude.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE SET_PROF( YLAT,      MONTH,  DAY,     T_CTM,  P_CTM,    &
                       CLDOD,     DSTOD,  AEROD,   O3_CTM, O3_TOMS,  &
                       AERCOL,    T_CLIM, O3_CLIM, Z_CLIM, AIR_CLIM, &
                       Input_Opt, State_Grid )
!
! !USES:
!
    USE CMN_SIZE_Mod,       ONLY : NAER, NRH
    USE Input_Opt_Mod,      ONLY : OptInput
    USE PhysConstants,      ONLY : AIRMW, AVO, g0, BOLTZ
    USE State_Grid_Mod,     ONLY : GrdState
!
! !INPUT PARAMETERS:
!
    REAL(fp), INTENT(IN)       :: YLAT              ! Latitude (degrees)
    INTEGER,  INTENT(IN)       :: MONTH             ! Month
    INTEGER,  INTENT(IN)       :: DAY               ! Day *of month*
    REAL(fp), INTENT(IN)       :: T_CTM(L1_)        ! CTM temperatures (K)
    REAL(fp), INTENT(IN)       :: O3_TOMS           ! O3 column (DU)
    REAL(fp), INTENT(IN)       :: P_CTM(L1_+1)      ! CTM edge pressures (hPa)
    REAL(fp), INTENT(INOUT)    :: CLDOD(L_)         ! Cloud optical depth
    REAL(fp), INTENT(IN)       :: DSTOD(L_,NDUST)   ! Mineral dust OD
    REAL(fp), INTENT(IN)       :: AEROD(L_,A_)      ! Aerosol OD
    REAL(fp), INTENT(IN)       :: O3_CTM(L1_)       ! CTM ozone (molec/cm3)
    TYPE(OptInput), INTENT(IN) :: Input_Opt         ! Input options
    TYPE(GrdState), INTENT(IN) :: State_Grid        ! Grid State object
!
! !OUTPUT VARIABLES:
!
    REAL(fp), INTENT(OUT)      :: AERCOL(A_,L1_)    ! Aerosol column
    REAL(fp), INTENT(OUT)      :: T_CLIM(L1_)       ! Clim. temperatures (K)
    REAL(fp), INTENT(OUT)      :: Z_CLIM(L1_+1)     ! Edge altitudes (cm)
    REAL(fp), INTENT(OUT)      :: O3_CLIM(L1_)      ! O3 column depth (#/cm2)
    REAL(fp), INTENT(OUT)      :: AIR_CLIM(L1_)     ! O3 column depth (#/cm2)
!
! !REMARKS:
!
! !REVISION HISTORY:
!  30 Mar 2013 - S. D. Eastham - Adapted from J. Mao code
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER                  :: I, K, L, M, N, LCTM
    REAL(fp)                 :: DLOGP,F0,T0,B0,PB,PC,XC,MASFAC,SCALEH
    REAL(fp)                 :: PSTD(52),OREF2(51),TREF2(51)
    REAL(fp)                 :: PROFCOL, ODSUM
    REAL(fp), PARAMETER      :: ODMAX = 200.0e+0_fp

    ! Local variables for quantities from Input_Opt
    LOGICAL :: USE_ONLINE_O3

    !=================================================================
    ! SET_PROF begins here!
    !=================================================================

    ! Copy fields from INPUT_OPT
    USE_ONLINE_O3   = Input_Opt%USE_ONLINE_O3

    ! Zero aerosol column
    DO K=1,A_
       DO I=1,L1_
          AERCOL(K,I) = 0.e+0_fp
       ENDDO
    ENDDO

    ! Scale optical depths to stay within limits
    ODSUM = 0.e+0_fp
    DO I=1,L_
       CLDOD(I) = DBLE(CLDOD(I))
       ODSUM = ODSUM + CLDOD(I)
    ENDDO
    IF (ODSUM.gt.ODMAX) THEN
       ODSUM = ODMAX/ODSUM ! Temporary
       DO I=1,L_
          CLDOD(I) = CLDOD(I)*ODSUM
       ENDDO
       ODSUM = ODMAX
    ENDIF

    !=================================================================
    ! Set up pressure levels for O3/T climatology - assume that value
    ! given for each 2 km z* level applies from 1 km below to 1 km
    ! above, so select pressures at these boundaries. Surface level
    ! values at 1000 mb are assumed to extend down to the actual
    ! surface pressure for this lat/lon.
    !=================================================================
    PSTD(1)  = MAX(P_CTM(1),1000.e+0_fp)
    PSTD(2)  = 1000.e+0_fp * 10.e+0_fp ** (-1.e+0_fp/16.e+0_fp)
    DLOGP    = 10.e+0_fp**(-2.e+0_fp/16.e+0_fp)
    DO I=3,51
       PSTD(I) = PSTD(I-1) * DLOGP
    ENDDO
    PSTD(52) = 0.e+0_fp

    ! Mass factor - delta-Pressure [hPa] to delta-Column [molec/cm2]
    MASFAC = 100.e+0_fp * AVO / ( AIRMW * g0 * 10.e+0_fp )

    ! Select appropriate monthly and latitudinal profiles
    ! Now use YLAT instead of Oliver's YDGRD(NSLAT) (bmy, 9/13/99)
    M = MAX( 1, MIN( 12, MONTH                   ) )
    L = MAX( 1, MIN( 18, ( INT(YLAT) + 99 ) / 10 ) )

    ! Temporary arrays for climatology data
    DO I = 1, 51
       OREF2(I) = OREF(I,L,M)
       TREF2(I) = TREF(I,L,M)
    ENDDO

    ! Apportion O3 and T on supplied climatology z* levels onto CTM levels
    ! with mass (pressure) weighting, assuming constant mixing ratio and
    ! temperature half a layer on either side of the point supplied.
    DO I = 1, L1_
       F0 = 0.e+0_fp
       T0 = 0.e+0_fp
       DO K = 1, 51
          PC = MIN( P_CTM(I),   PSTD(K)   )
          PB = MAX( P_CTM(I+1), PSTD(K+1) )
          IF ( PC .GT. PB ) THEN
             XC = ( PC - PB ) / ( P_CTM(I) - P_CTM(I+1) )
             F0 = F0 + OREF2(K)*XC
             T0 = T0 + TREF2(K)*XC
          ENDIF
       ENDDO
       T_CLIM(I)  = T0
       O3_CLIM(I) = F0 * 1.e-6_fp
    ENDDO

    !=================================================================
    ! Calculate effective altitudes using scale height at each level
    !=================================================================
    Z_CLIM(1) = 0.e+0_fp
    DO I = 1, L_
       SCALEH = BOLTZ * 1.e+4_fp * MASFAC * T_CLIM(I)

       Z_CLIM(I+1) = Z_CLIM(I) - ( LOG( P_CTM(I+1) / P_CTM(I) ) * SCALEH )
    ENDDO
    Z_CLIM(L1_+1)=Z_CLIM(L1_) + ZZHT

    !=================================================================
    ! Add Aerosol Column - include aerosol types here. Currently use
    ! soot water and ice; assume black carbon x-section of 10 m2/g,
    ! independent of wavelength; assume limiting temperature for
    ! ice of -40 deg C.
    !=================================================================
    DO I = 1, L_
       ! Turn off uniform black carbon profile (rvm, bmy, 2/27/02)
       AERCOL(1,I) = 0e+0_fp

       IF ( T_CTM(I) .GT. 233.e+0_fp ) THEN
          AERCOL(2,I) = CLDOD(I)
          AERCOL(3,I) = 0.e+0_fp
       ELSE
          AERCOL(2,I) = 0.e+0_fp
          AERCOL(3,I) = CLDOD(I)
       ENDIF

       ! Also add in aerosol optical depth columns (rvm, bmy, 9/30/00)
       DO N = 1, NDUST
          AERCOL(3+N,I) = DSTOD(I,N)
       ENDDO

       ! Also add in other aerosol optical depth columns (rvm, bmy, 2/27/02)
       DO N = 1, NAER*NRH
          AERCOL(3+N+NDUST,I) = AEROD(I,N)
       ENDDO

    ENDDO

    DO K = 1,(3+NDUST+(NAER))
       AERCOL(K,L1_    ) = 0.e+0_fp
    ENDDO

    !=================================================================
    ! Calculate column quantities for FAST-JX
    !=================================================================
    PROFCOL = 0e+0_fp

    DO I = 1, L1_

       ! Monthly mean air Column [molec/cm2]
       AIR_CLIM(I)  = ( P_CTM(I) - P_CTM(I+1) ) * MASFAC

       ! Monthly mean O3 column [molec/cm2]
       O3_CLIM(I) = O3_CLIM(I) * AIR_CLIM(I)

       ! Monthly mean O3 column [DU]
       PROFCOL = PROFCOL + ( O3_CLIM(I) / 2.69e+16_fp )
    ENDDO

    !! Top values are special (do not exist in CTM data)
    !AIR_CLIM(L1_)     = P_CTM(L1_) * MASFAC
    !O3_CLIM(L1_) = O3_CLIM(L1_) * AIR_CLIM(L1_)

    !=================================================================
    ! Now weight the O3 column by the observed monthly mean TOMS.
    ! Missing data is denoted by the flag -999. (mje, bmy, 7/15/03)
    !
    ! TOMS/SBUV MERGED TOTAL OZONE DATA, Version 8, Revision 3.
    ! Resolution:  5 x 10 deg.
    !
    ! Methodology (bmy, 2/12/07)
    ! ----------------------------------------------------------------
    ! FAST-J comes with its own default O3 column climatology (from
    ! McPeters 1992 & Nagatani 1991), which is stored in the input
    ! file "jv_atms.dat".  These "FAST-J default" O3 columns are used
    ! in the computation of the actinic flux and other optical
    ! quantities for the FAST-J photolysis.
    !
    ! The TOMS/SBUV O3 columns and 1/2-monthly O3 trends (contained
    ! in the TOMS_200701 directory) are read into GEOS-Chem by routine
    ! READ_TOMS in "toms_mod.f".  Missing values (i.e. locations where
    ! there are no data) in the TOMS/SBUV O3 columns are defined by
    ! the flag -999.
    !
    ! After being read from disk in routine READ_TOMS, the TOMS/SBUV
    ! O3 data are then passed to the FAST-J routine "set_prof.f".  In
    ! "set_prof.f", a test is done to make sure that the TOMS/SBUV O3
    ! columns and 1/2-monthly trends do not have any missing values
    ! for (lat,lon) location for the given month.  If so, then the
    ! TOMS/SBUV O3 column data is interpolated to the current day and
    ! is used to weight the "FAST-J default" O3 column.  This
    ! essentially "forces" the "FAST-J default" O3 column values to
    ! better match the observations, as defined by TOMS/SBUV.
    !
    ! If there are no TOMS/SBUV O3 columns (and 1/2-monthly trends)
    ! at a (lat,lon) location for given month, then FAST-J will revert
    ! to its own "default" climatology for that location and month.
    ! Therefore, the TOMS O3 can be thought of as an  "overlay" data
    ! -- it is only used if it exists.
    !
    ! Note that there are no TOMS/SBUV O3 columns at the higher
    ! latitudes.  At these latitudes, the code will revert to using
    ! the "FAST-J default" O3 columns.
    !
    ! As of February 2007, we have TOMS/SBUV data for 1979 thru 2005.
    ! 2006 TOMS/SBUV data is incomplete as of this writing.  For years
    ! 2006 and onward, we use 2005 TOMS O3 columns.
    !
    ! This methodology was originally adopted by Mat Evans.  Symeon
    ! Koumoutsaris was responsible for creating the downloading and
    ! processing the TOMS O3 data files from 1979 thru 2005 in the
    ! TOMS_200701 directory.
    !=================================================================

    ! Updated with UCX
    ! Since we now have stratospheric ozone calculated online, use
    ! this instead of archived profiles for all chemistry-grid cells
    ! The variable O3_CTM is obtained from State_Met%Species, and will be 0
    ! outside the chemgrid (in which case we use climatology)

    ! Scale monthly O3 profile to the daily O3 profile (if available)
    DO I = 1, L1_

       ! Use online O3 values in the chemistry grid if selected
       IF ( (USE_ONLINE_O3) .and. &
            (I <= State_Grid%MaxChemLev) .and. &
            (O3_CTM(I) > 0e+0_fp) ) THEN

          ! Convert from molec/cm3 to molec/cm2
          O3_CLIM(I) = O3_CTM(I) * (Z_CLIM(I+1)-Z_CLIM(I))

       ! Otherwise, use O3 values from the met fields or TOMS/SBUV
       ELSEIF (O3_TOMS > 0e+0_fp) THEN

          O3_CLIM(I) = O3_CLIM(I) * ( O3_TOMS / PROFCOL )

       ENDIF

    ENDDO

  END SUBROUTINE SET_PROF
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: set_aer
!
! !DESCRIPTION: Subroutine SET\_AER fills out the array MIEDX.
!  Each entry connects a GEOS-Chem aerosol to its Fast-JX counterpart:
!  MIEDX(Fast-JX index) = (GC index)
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE SET_AER( Input_Opt )
!
! !USES:
!
    USE CMN_SIZE_Mod,  ONLY : NRHAER, NRH
    USE Input_Opt_Mod, ONLY : OptInput
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN) :: Input_Opt ! Input options
!
! !REVISION HISTORY:
!  31 Mar 2013 - S. D. Eastham - Adapted from J. Mao FJX v6.2 implementation
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER               :: I, J, K
    INTEGER               :: IND(NRHAER)

    !=================================================================
    ! SER_AER begins here!
    !=================================================================

    ! Taken from aerosol_mod.F
    IND = (/22,29,36,43,50/)

    DO I=1,AN_
       MIEDX(I) = 0
    ENDDO

    ! Select Aerosol/Cloud types to be used - define types here
    ! Each of these types must be listed in the order used by OPMIE.F

    ! Clouds
    MIEDX(1)  =  3   !  Black carbon absorber
    MIEDX(2)  = 10   !  Water Cloud (Deirmenjian 8 micron)
    MIEDX(3)  = 14   !  Irregular Ice Cloud (Mishchenko)

    ! Dust
    MIEDX(4)  = 15   !  Mineral Dust  .15 micron    (rvm, 9/30/00)
    MIEDX(5)  = 16   !  Mineral Dust  .25 micron    (rvm, 9/30/00)
    MIEDX(6)  = 17   !  Mineral Dust  .4  micron    (rvm, 9/30/00)
    MIEDX(7)  = 18   !  Mineral Dust  .8  micron    (rvm, 9/30/00)
    MIEDX(8)  = 19   !  Mineral Dust 1.5  micron    (rvm, 9/30/00)
    MIEDX(9)  = 20   !  Mineral Dust 2.5  micron    (rvm, 9/30/00)
    MIEDX(10) = 21   !  Mineral Dust 4.0  micron    (rvm, 9/30/00)

    ! Aerosols
    DO I=1,NRHAER
       DO J=1,NRH
          MIEDX(10+((I-1)*NRH)+J)=IND(I)+J-1
       ENDDO
    ENDDO

    IF ( Input_Opt%LUCX ) THEN
       ! Stratospheric aerosols - SSA/STS and solid PSCs
       MIEDX(10+(NRHAER*NRH)+1) = 4  ! SSA/LBS/STS
       MIEDX(10+(NRHAER*NRH)+2) = 14 ! NAT/ice PSCs
    endif

    ! Ensure all 'AN_' types are valid selections
    do i=1,AN_
       IF (Input_Opt%amIRoot) write(6,1000) MIEDX(i),TITLEAA(MIEDX(i))
       if (MIEDX(i).gt.NAA.or.MIEDX(i).le.0) then
          if (Input_Opt%amIRoot) then
             write(6,1200) MIEDX(i),NAA
          endif
          CALL EXITC('Bad MIEDX value.')
       endif
    enddo

1000 format('Using Aerosol type: ',i3,1x,a)
1200 format('Aerosol type ',i3,' unsuitable; supplied values must be ', &
            'between 1 and ',i3)

  END SUBROUTINE SET_AER
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: rd_prof_nc
!
! !DESCRIPTION: Subroutine RAD\_PROF\_NC reads in the reference climatology
!  from a NetCDF file rather than an ASCII .dat.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE RD_PROF_NC( Input_Opt, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE Input_Opt_Mod, ONLY : OptInput
    USE m_netcdf_io_open
    USE m_netcdf_io_read
    USE m_netcdf_io_readattr
    USE m_netcdf_io_close
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)  :: Input_Opt   ! Input Options object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT) :: RC          ! Success or failure?
!
! !REMARKS:
!  This file was automatically generated by the Perl scripts in the
!  NcdfUtilities package (which ships w/ GEOS-Chem) and was subsequently
!  hand-edited.
!
! !REVISION HISTORY:
!  19 Apr 2012 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    LOGICAL            :: FileExists          ! Does input file exist?
    INTEGER            :: fId                 ! netCDF file ID

    ! Strings
    CHARACTER(LEN=255) :: nc_dir              ! netCDF directory name
    CHARACTER(LEN=255) :: nc_file             ! netCDF file name
    CHARACTER(LEN=255) :: nc_path             ! netCDF path name
    CHARACTER(LEN=255) :: v_name              ! netCDF variable name
    CHARACTER(LEN=255) :: a_name              ! netCDF attribute name
    CHARACTER(LEN=255) :: a_val               ! netCDF attribute value
    CHARACTER(LEN=255) :: FileMsg
    CHARACTER(LEN=255) :: ErrMsg
    CHARACTER(LEN=255) :: ThisLoc

    ! Arrays
    INTEGER            :: st3d(3), ct3d(3)    ! For 3D arrays

    !=================================================================
    ! RD_PROF_NC begins here!
    !=================================================================

    ! Initialize
    ! Assume success
    RC      = GC_SUCCESS
    ErrMsg  = ''
    ThisLoc = ' -> at RD_PROF_NC (in module GeosCore/fast_jx_mod.F90)'

    ! Directory and file names
    nc_dir  = TRIM( Input_Opt%CHEM_INPUTS_DIR ) // 'FastJ_201204/'
    nc_file = 'fastj.jv_atms_dat.nc'
    nc_path = TRIM( nc_dir ) // TRIM( nc_file )

    !=================================================================
    ! In dry-run mode, print file path to dryrun log and exit.
    ! Otherwise, print file path to stdout and continue.
    !=================================================================

    ! Test if the file exists
    INQUIRE( FILE=TRIM( nc_path ), EXIST=FileExists )

    ! Test if the file exists and define an output string
    IF ( FileExists ) THEN
       FileMsg = 'FAST-JX (RD_PROF_NC): Opening'
    ELSE
       FileMsg = 'FAST-JX (RD_PROF_NC): REQUIRED FILE NOT FOUND'
    ENDIF

    ! Write to stdout for both regular and dry-run simulations
    IF ( Input_Opt%amIRoot ) THEN
       WRITE( 6, 300 ) TRIM( FileMsg ), TRIM( nc_path )
300    FORMAT( a, ' ', a )
    ENDIF

    ! For dry-run simulations, return to calling program.
    ! For regular simulations, throw an error if we can't find the file.
    IF ( Input_Opt%DryRun ) THEN
       RETURN
    ELSE
       IF ( .not. FileExists ) THEN
          WRITE( ErrMsg, 300 ) TRIM( FileMsg ), TRIM( nc_path )
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
    ENDIF

    !=========================================================================
    ! Open and read data from the netCDF file
    !=========================================================================

    ! Open netCDF file
    CALL Ncop_Rd( fId, TRIM(nc_path) )

    ! Echo info to stdout
    IF ( Input_Opt%amIRoot ) THEN
       WRITE( 6, 100 ) REPEAT( '%', 79 )
       WRITE( 6, 110 ) TRIM(nc_file)
       WRITE( 6, 120 ) TRIM(nc_dir)
    ENDIF

    !----------------------------------------
    ! VARIABLE: T
    !----------------------------------------

    ! Variable name
    v_name = "T"

    ! Read T from file
    st3d   = (/  1,  1,  1 /)
    ct3d   = (/ 51, 18, 12 /)
    CALL NcRd( TREF, fId, TRIM(v_name), st3d, ct3d )

    ! Read the T:units attribute
    a_name = "units"
    CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )

    ! Echo info to stdout
    IF ( Input_Opt%amIRoot ) THEN
       WRITE( 6, 130 ) TRIM(v_name), TRIM(a_val)
    ENDIF

    !----------------------------------------
    ! VARIABLE: O3
    !----------------------------------------

    ! Variable name
    v_name = "O3"

    ! Read O3 from file
    st3d   = (/  1,  1,  1 /)
    ct3d   = (/ 51, 18, 12 /)
    CALL NcRd( OREF, fId, TRIM(v_name), st3d, ct3d )

    ! Read the O3:units attribute
    a_name = "units"
    CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )

    ! Echo info to stdout
    IF ( Input_Opt%amIRoot ) THEN
       WRITE( 6, 130 ) TRIM(v_name), TRIM(a_val)
    ENDIF

    !=================================================================
    ! Cleanup and quit
    !=================================================================

    ! Close netCDF file
    CALL NcCl( fId )

    ! Echo info to stdout
    IF ( Input_Opt%amIRoot ) THEN
       WRITE( 6, 140 )
       WRITE( 6, 100 ) REPEAT( '%', 79 )
    ENDIF

    ! FORMAT statements
100 FORMAT( a                                              )
110 FORMAT( '%% Opening file  : ',         a               )
120 FORMAT( '%%  in directory : ',         a, / , '%%'     )
130 FORMAT( '%% Successfully read ',       a, ' [', a, ']' )
140 FORMAT( '%% Successfully closed file!'                 )

  END SUBROUTINE RD_PROF_NC
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: photrate_adj
!
! !DESCRIPTION: Subroutine PHOTRATE\_ADJ adjusts certain photolysis rates
!  for chemistry.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE PHOTRATE_ADJ( Input_Opt, State_Diag,       &
                           I,         J,         L,     &
                           NUMDEN,    TEMP,      C_H2O, &
                           FRAC,      RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE Input_Opt_Mod,  ONLY : OptInput
    USE State_Diag_Mod, ONLY : DgnState
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)    :: Input_Opt  ! Input_Options object
    INTEGER,        INTENT(IN)    :: I, J, L    ! Lon, lat, lev indices
    REAL(fp),       INTENT(IN)    :: NUMDEN     ! Air # density [molec/m3]
    REAL(fp),       INTENT(IN)    :: TEMP       ! Temperature [K]
    REAL(fp),       INTENT(IN)    :: C_H2O      ! H2O conc [molec/cm3]
    REAL(fp),       INTENT(IN)    :: FRAC       ! Result of SO4_PHOTFRAC,
                                                !  called from DO_FLEXCHEM
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(DgnState), INTENT(INOUT) :: State_Diag ! Diagnostics State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC         ! Success or failure
!
! !REMARKS:
!  NOTE: The netCDF diagnostics are attached in DO_FLEXCHEM so that we have
!  access to the adjusted rates.  Only the bpch diagnostics are updated
!  here.
!    -- Bob Yantosca, 19 Dec 2017
!
!  %%%% NOTE: WE SHOULD UPDATE THE COMMENTS TO MAKE SURE THAT WE DO      %%%%
!  %%%% NOT KEEP ANY CONFLICTING OR INCORRECT INFORMATION (bmy, 3/28/16) %%%%
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    LOGICAL  :: DO_ND22
    REAL(fp) :: C_O2,      C_N2,     C_H2, ITEMPK
    REAL(fp) :: RO1DplH2O, RO1DplH2, RO1D

    !=================================================================
    ! PHOTRATE_ADJ begins here!
    !=================================================================

    ! Initialize
    RC      = GC_SUCCESS

    ! For all mechanisms. Set the photolysis rate of NITs and NIT to a
    ! scaled value of JHNO3. NOTE: this is set in input.geos
    IF ( Input_Opt%hvAerNIT ) THEN

       ! Get the photolysis scalars read in from input.geos
       JscaleNITs = Input_Opt%hvAerNIT_JNITs
       JscaleNIT  = Input_Opt%hvAerNIT_JNIT
       ! convert reaction channel % to a fraction
       JNITChanA  = Input_Opt%JNITChanA
       JNITChanB  = Input_Opt%JNITChanB
       JNITChanA  = JNITChanA / 100.0_fp
       JNITChanB  = JNITChanB / 100.0_fp
       ! Set the photolysis rate of NITs
       ZPJ(L,RXN_JNITSa,I,J) = ZPJ(L,RXN_JHNO3,I,J) * JscaleNITs
       ZPJ(L,RXN_JNITSb,I,J) = ZPJ(L,RXN_JHNO3,I,J) * JscaleNITs
       ! Set the photolysis rate of NIT
       ZPJ(L,RXN_JNITa,I,J) = ZPJ(L,RXN_JHNO3,I,J) * JscaleNIT
       ZPJ(L,RXN_JNITb,I,J) = ZPJ(L,RXN_JHNO3,I,J) * JscaleNIT
       ! Adjust to scaling for channels set in input.geos
       ! NOTE: channel scaling is 1 in FJX_j2j.dat, then updated here
       ZPJ(L,RXN_JNITSa,I,J) = ZPJ(L,RXN_JNITSa,I,J) * JNITChanA
       ZPJ(L,RXN_JNITa,I,J) = ZPJ(L,RXN_JNITa,I,J) * JNITChanA
       ZPJ(L,RXN_JNITSb,I,J) = ZPJ(L,RXN_JNITSb,I,J) * JNITChanB
       ZPJ(L,RXN_JNITb,I,J) = ZPJ(L,RXN_JNITb,I,J) * JNITChanB

    ! Gotcha to set JNIT and JNITs to zero if hvAerNIT switch is off
    ELSE

       ! Set the photolysis rate of NITs to zero
       ZPJ(L,RXN_JNITSa,I,J) = 0.0_fp
       ZPJ(L,RXN_JNITSb,I,J) = 0.0_fp
       ! Set the photolysis rate of NIT to zero
       ZPJ(L,RXN_JNITa,I,J) = 0.0_fp
       ZPJ(L,RXN_JNITb,I,J) = 0.0_fp

    ENDIF

    ! Test if the UCX mechanism is being used
    IF ( Input_Opt%LUCX ) THEN

       !==============================================================
       ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       ! %%% FOR MECHANISMS WITH UCX                              %%%
       ! %%% (standard, benchmark, *SOA*, marinePOA, aciduptake)  %%%
       ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       !
       ! SPECIAL TREATMENT FOR H2SO4+hv -> SO2 + 2OH
       !
       ! Only allow photolysis of H2SO4 when gaseous (SDE 04/11/13)
       !==============================================================

       ! Calculate if H2SO4 expected to be gaseous or aqueous
       ! Only allow photolysis above 6 hPa
       ! RXN_H2SO4 specifies SO4 + hv -> SO2 + OH + OH
       ZPJ(L,RXN_H2SO4,I,J) = ZPJ(L,RXN_H2SO4,I,J) * FRAC

       !==============================================================
       ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       ! %%% FOR MECHANISMS WITH UCX                              %%%
       ! %%% (standard, benchmark, *SOA*, marinePOA, aciduptake)  %%%
       ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       !
       ! SPECIAL TREATMENT FOR O3+hv -> O+O2  (UCX simulation)
       !
       ! [O1D]ss=J[O3]/(k[H2O]+k[N2]+k[O2])
       ! SO, THE EFFECTIVE J-VALUE IS J*k[H2O]/(k[H2O]+k[N2]+k[O2])
       !
       ! We don't want to do this if strat-chem is in use, as all
       ! the intermediate reactions are included - this would be
       ! double-counting (SDE 04/01/13)
       !==============================================================

       ! Need to subtract O3->O1D from rate
       ! RXN_O3_1  specifies: O3 + hv -> O2 + O
       ! RXN_O3_2a specifies: O3 + hv -> O2 + O(1D)
       ZPJ(L,RXN_O3_1,I,J) = ZPJ(L,RXN_O3_1,I,J) &
                           - ZPJ(L,RXN_O3_2a,I,J)

    ELSE

       !==============================================================
       ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       ! %%% FOR MECHANISMS WITHOUT UCX %%%
       ! %%% (tropchem)                 %%%
       ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       !
       ! Change rate of O(1D)+ N2 to be 3.1e-11 at 298K rather
       ! than 2.6e-11.  The temperature dependence remains the
       ! same, so the constant changes from 1.8e-11 to 2.14e-11
       ! according to Heard, pers. comm.,2002. (amf, bmy, 1/7/02)
       !==============================================================
       ! Change the rate of O(1D)+H2O from 2.2e-10 to 1.45e-10*
       ! exp(89/temp) on the basis of Dunlea and Ravishankara
       ! 'Measurement of the Rate coefficient for the reaction
       ! of O(1D) with H2O and re-evaluation of the atmospheric
       ! OH Production Rate'.  One of the RSC Journals
       ! (mje 4/5/04)
       !==============================================================
       ! Updated from JPL2006, the difference is pretty small.
       ! (jmao,02/26/2009)
       !==============================================================
       ! Additional update from JPL 10-6 for reaction of
       ! O3 + hv --> HO2 + OH and O1D + H2:
       ! Includes calculation of k[O3], where
       ! k[O3]  = J[O3]*1.2e-10/k[O1D], where
       ! k[O1D] = ([O1D]*[H2]+[O1D]*[H2O])/
       !          ([O1D]*[H2]+[O1D]*[H2O]+[O1D][N2]+[O1D][O2])
       ! (bhh, jmao, eam, 7/18/11)
       !==============================================================

       ! Inverse temperature [K-1]
       ITEMPK    = 1.0_fp / TEMP

       ! Set species concentrations [molec/m3] ???
       C_O2      = 0.2095e+0_fp * NUMDEN
       C_N2      = 0.7808e+0_fp * NUMDEN

       ! Added H2 concentration (bhh, jmao, eam, 7/18/11)
       ! Seasonal variability of H2 may be important,
       ! but not included in this update (bhh, jmao, eam, 7/18/11)
       C_H2      = 0.5000e-6_fp * NUMDEN

       RO1DplH2O = 1.63e-10_fp * EXP(  60.0_fp * ITEMPK ) * C_H2O

       RO1DplH2  = 1.2e-10                                * C_H2

       RO1D      = RO1DplH2O &
                 + RO1DplH2  &
                 + 2.15e-11_fp * EXP( 110.0_fp * ITEMPK ) * C_N2 &
                 + 3.30e-11_fp * EXP(  55.0_fp * ITEMPK ) * C_O2

       ! Prevent div-by-zero
       IF ( RO1D > 0.0_fp ) THEN

          ! RXN_O3_2a specifies: O3 + hv -> O2 + O(1D) #1
          ZPJ(L,RXN_O3_2a,I,J) = ZPJ(L,RXN_O3_2a,I,J) * RO1DplH2O / RO1D

          ! RXN_O3_2b specifies: O3 + hv -> O2 + O(1D) #2
          ZPJ(L,RXN_O3_2b,I,J) = ZPJ(L,RXN_O3_2b,I,J) * RO1DplH2  / RO1D

       ENDIF

    ENDIF

  END SUBROUTINE PHOTRATE_ADJ
!EOC
END MODULE FAST_JX_MOD
