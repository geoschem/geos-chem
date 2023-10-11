!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: fjx_interface_mod.F90
!
! !DESCRIPTION: Module FJX\_INTERFACE\_MOD contains routines and variables
!  for interfacing with the Fast-JX scheme (Prather et al) that calculates
!  photolysis rates. Current implementation is version 7.0a.
!\\
!\\
! !INTERFACE:
!
MODULE FJX_INTERFACE_MOD
!
! !USES:
!
  USE FJX_Mod
  USE PRECISION_MOD    ! For GEOS-Chem Precision (fp)

  IMPLICIT NONE

  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC  :: Init_FastJX
  PUBLIC  :: Run_FASTJX
!
! !PRIVATE MEMBER FUNCTIONS:
!
  PRIVATE :: RD_PROF_NC ! NC-read version of what is in original fast-jx
  PRIVATE :: SET_PROF_FJX ! could consolidate in photolysis_mod perhaps
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC

CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_fastjx
!
! !DESCRIPTION: Subroutine Init\_FastJX initializes Fast-JX variables.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Init_FastJX( Input_Opt, State_Diag, State_Chm, RC )
!
! !USES:
!
    USE CMN_FJX_Mod,    ONLY : JVN_, NJX, NRATJ, W_, WL
    USE CMN_FJX_Mod,    ONLY : TITLEJX, JLABEL, RNAMES, JFACTA
    USE ErrCode_Mod
    USE Input_Opt_Mod,  ONLY : OptInput
    USE inquireMod,     ONLY : findFreeLUN
    USE State_Chm_Mod,  ONLY : ChmState
    USE State_Chm_Mod,  ONLY : Ind_
    USE State_Diag_Mod, ONLY : DgnState
#if defined( MODEL_CESM )
    USE UNITS,          ONLY : freeUnit
#endif

!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
    TYPE(DgnState), INTENT(IN)    :: State_Diag  ! Diagnostics State object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
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
    INTEGER            :: JXUNIT, J, NJXX, K

    ! Strings
    CHARACTER(LEN=6)   :: TITLEJXX(JVN_)
    CHARACTER(LEN=255) :: DATA_DIR
    CHARACTER(LEN=255) :: FILENAME
    CHARACTER(LEN=255) :: ErrMsg, ThisLoc

    !=================================================================
    ! Init_FastJX begins here!
    !=================================================================

    ! Initialize
    RC        = GC_SUCCESS
    notDryRun = ( .not. Input_Opt%DryRun )
    ErrMsg    = ''
    ThisLoc   = ' -> at Init_FastJX (in module GeosCore/fjx_interface_mod.F90)'

    ! Skip these opterations when running in dry-run mode
    IF ( notDryRun ) THEN

       ! Print info
       IF ( Input_Opt%amIRoot ) THEN
          write(6,*) ' Initializing Fast-JX v7.0 standalone CTM code.'

          if (W_.ne.8 .and. W_.ne.12 .and. W_.ne.18) then
             ErrMsg =  'Invalid number of wavelengths (W_)'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          endif
       ENDIF

#if defined( MODEL_CESM )
       IF ( Input_Opt%amIRoot ) THEN
          JXUNIT = findFreeLUN()
       ELSE
          JXUNIT = 0
       ENDIF
#else
       ! Get a free LUN
       JXUNIT = findFreeLUN()
#endif

    ENDIF

    ! Define data directory for FAST-JX input
    DATA_DIR = TRIM( Input_Opt%FAST_JX_DIR )

    !=====================================================================
    ! Read in fast-J X-sections (spectral data)
    !=====================================================================
    FILENAME = TRIM( DATA_DIR ) // 'FJX_spec.dat'

    ! Read file, or just print filename if we are in dry-run mode
    CALL RD_XXX( Input_Opt%amIRoot, Input_Opt%DryRun, JXUNIT, &
                 TRIM( FILENAME ), RC)

    ! Trap potential errors
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered in FAST-JX routine "RD_XXX"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    !=====================================================================
    ! Read in 5-wavelength scattering data
    ! (or just print file name if in dry-run mode)
    !=====================================================================
    FILENAME = TRIM( DATA_DIR ) // 'jv_spec_mie.dat'

    ! Read data
    CALL RD_MIE( Input_Opt%amIRoot, Input_Opt%DryRun, Input_Opt%LBRC, &
                 JXUNIT, TRIM( FILENAME ), RC )

    ! Trap potential errors
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered in FAST-JX routine "RD_MIE"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    !=====================================================================
    ! Read in T & O3 climatology used to fill e.g. upper layers
    ! or if O3 not calc.
    !=====================================================================
    CALL RD_PROF_NC( Input_Opt, State_Chm, RC )

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
    CALL RD_JS_JX( Input_Opt%amIRoot, Input_Opt%DryRun, JXUNIT, &
                   TRIM( FILENAME ), TITLEJXX, NJXX, RC )

    ! Store # of photolysis reactions in state_chm for easy reference
    State_Chm%Phot%nPhotRxns = NRatJ

#if defined( MODEL_CESM )
    IF ( notDryRun .AND. Input_Opt%amIRoot ) THEN
       CALL freeUnit(JXUnit)
    ENDIF
#endif

  END SUBROUTINE Init_FastJX
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: run_fastjx
!
! !DESCRIPTION: Subroutine RUN|_FASTJX loops over longitude and latitude, and
!  calls PHOTO\_JX to compute J-Values for each column at every chemistry
!  time-step.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Run_FastJX( WLAOD, Input_Opt, State_Chm, State_Diag, &
                         State_Grid, State_Met, RC )
!
! !USES:
!
    USE CMN_FJX_Mod,     ONLY : A_, L_, L1_, W_, JVN_, JXL_, JXL1_
    USE CMN_FJX_Mod,     ONLY : NRATJ, JIND, JFACTA, FL
    USE CMN_SIZE_MOD,    ONLY : NDUST, NRH, NAER
    USE ErrCode_Mod
    USE ERROR_MOD,       ONLY : ERROR_STOP, ALLOC_ERR
    USE ERROR_MOD,       ONLY : DEBUG_MSG
    USE Input_Opt_Mod,   ONLY : OptInput
    USE State_Chm_Mod,   ONLY : ChmState
    USE State_Chm_Mod,   ONLY : Ind_
    USE State_Diag_Mod,  ONLY : DgnState
    USE State_Grid_Mod,  ONLY : GrdState
    USE State_Met_Mod,   ONLY : MetState
    USE TIME_MOD,        ONLY : GET_MONTH, GET_DAY, GET_DAY_OF_YEAR
    USE TIME_MOD,        ONLY : GET_TAU,   GET_YEAR
    USE TOMS_MOD,        ONLY : GET_OVERHEAD_O3

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
    INTEGER       :: NLON, NLAT, DAY,  MONTH, DAY_OF_YR, L, N, J
    INTEGER       :: IOPT, LCHEM
    REAL(fp)      :: U0, PRES, RFL, YLAT,  O3_TOMS, SZA, SOLF
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
    REAL(fp)      :: AERX_COL(A_,L1_)
    REAL(fp)      :: T_CLIM(L1_)
    REAL(fp)      :: O3_CLIM(L1_)
    REAL(fp)      :: Z_CLIM(L1_+1)
    REAL(fp)      :: AIR_CLIM(L1_)
    REAL(fp)      :: VALJXX(L_,JVN_)

    LOGICAL       :: AOD999
    LOGICAL, SAVE :: FIRST = .true.
    LOGICAL       :: prtDebug

    ! Species ID flags
    INTEGER, SAVE :: id_O3

    ! Strings
    CHARACTER(LEN=255) :: ErrMsg, ThisLoc

    ! Input to set_prof
    REAL(fp) :: ODCLOUD_COL(L_)

    ! Output from photo_jx for use in diagnostics
    REAL(fp), DIMENSION(W_)         :: FJBOT,FSBOT
    REAL(fp), DIMENSION(JXL1_,W_)   :: FLXD
    REAL(fp), DIMENSION(JXL_, W_)   :: FJFLX

    ! UVFlux* diagnostics
    REAL(fp) :: FDIRECT (JXL1_)
    REAL(fp) :: FDIFFUSE(JXL1_)
    REAL(fp) :: UVX_CONST
    INTEGER  :: S, K

    ! Pointers
    INTEGER,  POINTER :: IRHARR   (:,:,:)
    REAL(fp), POINTER :: UVXFACTOR(:)
    REAL(fp), POINTER :: ZPJ      (:,:,:,:)
    REAL(fp), POINTER :: ODAER    (:,:,:,:,:)
    REAL(fp), POINTER :: ODMDUST  (:,:,:,:,:)

    !=================================================================
    ! Run_FastJX begins here!
    !=================================================================

    ! Initialize
    RC        = GC_SUCCESS
    ErrMsg    = ''
    ThisLoc   = ' -> at Run_FastJX (in module GeosCore/fjx_interface_mod.F90)'
    prtDebug  = Input_Opt%Verbose

    ! Set pointers
    ZPJ       => State_Chm%Phot%ZPJ
    IRHARR    => State_Chm%Phot%IRHARR
    UVXFACTOR => State_Chm%Phot%UVXFACTOR
    ODAER     => State_Chm%Phot%ODAER
    ODMDUST   => State_Chm%Phot%ODMDUST

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
    !$OMP PRIVATE( NLAT,    NLON,   YLAT,      U0,      L       ) &
    !$OMP PRIVATE( P_CTM ,  T_CTM,  RFL,       O3_TOMS, O3_CTM  ) &
    !$OMP PRIVATE( LCHEM,   OPTAER, N,         IOPT,    J       ) &
    !$OMP PRIVATE( OPTDUST, OPTD,   CLDF1D                      ) &
#ifdef USE_MAXIMUM_RANDOM_OVERLAP
    !$OMP PRIVATE( FMAX,    KK,     NUMB,      KBOT             ) &
    !$OMP PRIVATE( KTOP     ODNEW,  INDICATOR, INDIC            ) &
#endif
    !$OMP PRIVATE( SZA, SOLF, ODCLOUD_COL                       ) &
    !$OMP PRIVATE( AERX_COL,  T_CLIM, O3_CLIM, Z_CLIM, AIR_CLIM ) &
    !$OMP PRIVATE( VALJXX,    FSBOT,  FJBOT,   FLXD,   FJFLX    ) &
    !$OMP PRIVATE( FDIRECT,   FDIFFUSE, UVX_CONST, K, S         ) &
    !$OMP SCHEDULE( DYNAMIC )

    ! Loop over latitudes and longitudes
    DO NLAT = 1, State_Grid%NY
    DO NLON = 1, State_Grid%NX

       ! Grid box latitude [degrees]
       YLAT = State_Grid%YMid( NLON, NLAT )

       ! Cosine of solar zenith angle [unitless] at (NLON,NLAT)
       U0 = State_Met%SUNCOSmid(NLON,NLAT)

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
       RFL = State_Met%UVALBEDO(NLON,NLAT)

       ! Overhead ozone column [DU] at (NLON, NLAT)
       ! These values are either from the met fields or TOMS/SBUV,
       ! depending on the settings in geoschem_config.yml
       O3_TOMS = GET_OVERHEAD_O3( State_Chm, NLON, NLAT )

       ! CTM ozone densities (molec/cm3) at (NLON, NLAT)
       O3_CTM = 0e+0_fp
       LCHEM  = State_Met%ChemGridLev(NLON,NLAT)
       DO L = 1, LCHEM
          O3_CTM(L) = State_Chm%Species(id_O3)%Conc(NLON,NLAT,L)
       ENDDO

       ! Aerosol OD profile [unitless] at (NLON,NLAT)
       ! and at 1000nm, IWV1000 (DAR)
       !OPTAER wants NAER*NRH values but ODAER is now NAER
       !use IRHARR to map to correct OPTAER bin (DAR 08/13)
       OPTAER = 0.0e+0_fp
       DO N = 1, NAER
       DO L = 1, State_Grid%NZ
          IOPT = ( (N-1) * NRH ) + IRHARR(NLON,NLAT,L)
          OPTAER(L,IOPT) = ODAER(NLON,NLAT,L,State_Chm%Phot%IWV1000,N)
       ENDDO
       ENDDO
       DO N = 1, NDUST
       DO L = 1, State_Grid%NZ
          OPTDUST(L,N) = ODMDUST(NLON,NLAT,L,State_Chm%Phot%IWV1000,N)
       ENDDO
       ENDDO

       ! Mineral dust OD profile [unitless] at (NLON,NLAT)
       ! and at 1000nm, IWV1000 (DAR)
       !OPTDUST = ODMDUST(NLON,NLAT,:,State_Chm%Phot%IWV1000,:)

       ! Cloud OD profile [unitless] at (NLON,NLAT)
       OPTD = State_Met%OPTD(NLON,NLAT,1:State_Grid%NZ)

       !-----------------------------------------------------------
       !### If you want to exclude aerosol OD, mineral dust OD,
       !### or cloud OD, then uncomment the following lines:
       !OPTAER  = 0d0
       !OPTDUST = 0d0
       !OPTD(:)    = 0d0
       !-----------------------------------------------------------

#if defined( MODEL_GEOS )
       ! Initialize diagnostics arrays
       IF ( State_Diag%Archive_EXTRALNLEVS ) THEN
          State_Diag%EXTRALNLEVS(NLON,NLAT) = 0.0
       ENDIF
       IF ( State_Diag%Archive_EXTRALNITER ) THEN
          State_Diag%EXTRALNITER(NLON,NLAT) = 0.0
       ENDIF
#endif

       if (State_Grid%NZ+1 .gt. JXL1_) then
          ErrMsg = ' PHOTO_JX: not enough levels in JX'
          call Error_Stop( ErrMsg, ThisLoc )
       endif

       ! Input conversion (SDE 03/29/13)
       ! Calculate solar zenith angle (degrees)
       CALL SOLAR_JX(DAY_OF_YR,U0,SZA,SOLF)
       
       ! check for dark conditions SZA > 98.0 deg => tan ht = 63 km or
       !                                 99.0                 80 km
       if (SZA .gt. 98.e+0_fp) cycle

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

#elif defined( USE_MAXIMUM_RANDOM_OVERLAP )
       ! See commented out code in github history
       CALL ERROR_STOP('MMRAN_16 not yet FJX compatible.', &
                       'fjx_interface_mod.F90')


#endif

       ! Copy cloud OD data to a variable array
       DO L=1,L_
          ODCLOUD_COL(L) = OPTD(L)
       ENDDO
       
       ! Use GEOS-Chem methodology to set vertical profiles of:
       ! Pressure      (PPJ)    [hPa]
       ! Temperature   (T_CLIm) [K]
       ! Path density  (DDJ)    [# molec/cm2]
       ! New methodology for:
       ! Ozone density (OOJ)    [# O3 molec/cm2]
       CALL SET_PROF_FJX (YLAT,        MONTH,      DAY,         &
                      T_CTM,       P_CTM,      OPTD,        &
                      OPTDUST,     OPTAER,     O3_CTM,      &
                      O3_TOMS,     AERX_COL,   T_CLIM,      &
                      O3_CLIM,     Z_CLIM,     AIR_CLIM,    &
                      Input_Opt,   State_Grid, State_Chm )

       ! Call FAST-JX routines to compute J-values
       CALL PHOTO_JX( Input_Opt%amIRoot, Input_Opt%DryRun,          &
                      U0,        RFL,        SZA,       SOLF,       &
                      P_CTM,     T_CTM,      AOD999,    NLON,       &
                      NLAT,      AERX_COL,   T_CLIM,    O3_CLIM,    &
                      Z_CLIM,    AIR_CLIM,   State_Grid%maxChemLev, &
                      State_Chm, VALJXX,     FSBOT,     FJBOT,      &
                      FLXD,      FJFLX,      Input_Opt, State_Diag   )

       ! Fill out common-block array of J-rates using PHOTO_JX output
       DO L=1,State_Grid%MaxChemLev
          DO J=1,NRATJ
             IF (JIND(J).gt.0) THEN
                ZPJ(L,J,NLON,NLAT) = VALJXX(L,JIND(J))*JFACTA(J)
             ELSE
                ZPJ(L,J,NLON,NLAT) = 0.e+0_fp
             ENDIF
          ENDDO
       ENDDO

       ! Set J-rates outside the chemgrid to zero
       IF (State_Grid%MaxChemLev.lt.L_) THEN
          DO L=State_Grid%MaxChemLev+1,L_
             DO J=1,NRATJ
                ZPJ(L,J,NLON,NLAT) = 0.e+0_fp
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

       ! TODO TODO TODO - there seems to be a parallelization error here

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
                   S = State_Diag%Map_UvFluxNet%id2slot(K)
                   IF ( S > 0 ) THEN
                      State_Diag%UVFluxNet(NLON,NLAT,L,S) =  &
                      State_Diag%UVFluxNet(NLON,NLAT,L,S) +  &
                           ( ( FDIRECT(L) + FDIFFUSE(L) ) * UVX_CONST )
                   ENDIF
                ENDIF
       
                IF ( State_Diag%Archive_UVFluxDirect ) THEN
                   S = State_Diag%Map_UvFluxDirect%id2slot(K)
                   IF ( S > 0 ) THEN
                      State_Diag%UVFluxDirect(NLON,NLAT,L,S) =  &
                      State_Diag%UVFluxDirect(NLON,NLAT,L,S) +  &
                           ( FDIRECT(L) * UVX_CONST )
                   ENDIF
                ENDIF
       
                IF ( State_Diag%Archive_UVFluxDiffuse ) THEN
                   S = State_Diag%Map_UvFluxDiffuse%id2slot(K)
                   IF ( S > 0 ) THEN
                      State_Diag%UVFluxDiffuse(NLON,NLAT,L,S) =  &
                      State_Diag%UVFluxDiffuse(NLON,NLAT,L,S) +  &
                           ( FDIFFUSE(L) * UVX_CONST )
                   ENDIF
                ENDIF
             ENDDO
          ENDDO
       ENDIF

    ENDDO
    ENDDO
    !$OMP END PARALLEL DO

    ! Reset first-time flag
    FIRST=.FALSE.

    ! Free pointers
    ZPJ       => NULL()
    IRHARR    => NULL()
    UVXFACTOR => NULL()
    ODAER     => NULL()
    ODMDUST   => NULL()

  END SUBROUTINE Run_FastJX
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
  SUBROUTINE RD_PROF_NC( Input_Opt, State_Chm, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE Input_Opt_Mod, ONLY : OptInput
    USE State_Chm_Mod, ONLY : ChmState

#if defined( MODEL_CESM )
    USE CAM_PIO_UTILS,     ONLY : CAM_PIO_OPENFILE
    USE IOFILEMOD,         ONLY : GETFIL
    USE PIO,               ONLY : PIO_CLOSEFILE
    USE PIO,               ONLY : PIO_INQ_DIMID
    USE PIO,               ONLY : PIO_INQ_DIMLEN
    USE PIO,               ONLY : PIO_INQ_VARID
    USE PIO,               ONLY : PIO_GET_VAR
    USE PIO,               ONLY : PIO_NOERR
    USE PIO,               ONLY : PIO_NOWRITE
    USE PIO,               ONLY : FILE_DESC_T
#else
    USE m_netcdf_io_open
    USE m_netcdf_io_read
    USE m_netcdf_io_readattr
    USE m_netcdf_io_close
#endif
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry state object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
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

#if defined( MODEL_CESM )
    type(FILE_DESC_T)  :: ncid
    INTEGER            :: vId, iret
#endif

    ! Pointers
    REAL(fp), POINTER :: OREF(:,:,:)
    REAL(fp), POINTER :: TREF(:,:,:)

    !=================================================================
    ! RD_PROF_NC begins here!
    !=================================================================

    ! Initialize
    ! Assume success
    RC      = GC_SUCCESS
    ErrMsg  = ''
    ThisLoc = ' -> at RD_PROF_NC (in module GeosCore/fjx_interface_mod.F90)'

    ! Set pointers
    OREF => State_Chm%Phot%OREF
    TREF => State_Chm%Phot%TREF

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
#if defined( MODEL_CESM )
    CALL CAM_PIO_OPENFILE( ncid, TRIM(nc_path), PIO_NOWRITE )
#else
    CALL Ncop_Rd( fId, TRIM(nc_path) )
#endif

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
#if defined( MODEL_CESM )
    iret = PIO_INQ_VARID( ncid, trim(v_name), vid  )
    iret = PIO_GET_VAR(   ncid, vid, TREF          )
#else
    CALL NcRd( TREF, fId, TRIM(v_name), st3d, ct3d )

    ! Read the T:units attribute
    a_name = "units"
    CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )

    ! Echo info to stdout
    IF ( Input_Opt%amIRoot ) THEN
       WRITE( 6, 130 ) TRIM(v_name), TRIM(a_val)
    ENDIF
#endif

    !----------------------------------------
    ! VARIABLE: O3
    !----------------------------------------

    ! Variable name
    v_name = "O3"

    ! Read O3 from file
    st3d   = (/  1,  1,  1 /)
    ct3d   = (/ 51, 18, 12 /)
#if defined( MODEL_CESM )
    iret = PIO_INQ_VARID( ncid, trim(v_name), vid  )
    iret = PIO_GET_VAR(   ncid, vid, OREF          )
#else
    CALL NcRd( OREF, fId, TRIM(v_name), st3d, ct3d )

    ! Read the O3:units attribute
    a_name = "units"
    CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )

    ! Echo info to stdout
    IF ( Input_Opt%amIRoot ) THEN
       WRITE( 6, 130 ) TRIM(v_name), TRIM(a_val)
    ENDIF
#endif

    !=================================================================
    ! Cleanup and quit
    !=================================================================

    ! Close netCDF file
#if defined( MODEL_CESM )
    CALL PIO_CLOSEFILE( ncid )
#else
    CALL NcCl( fId )
#endif

    ! Echo info to stdout
    IF ( Input_Opt%amIRoot ) THEN
       WRITE( 6, 140 )
       WRITE( 6, 100 ) REPEAT( '%', 79 )
    ENDIF

    ! Free pointers
    OREF => NULL()
    TREF => NULL()

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
! !IROUTINE: set_prof_fjx
!
! !DESCRIPTION: Subroutine SET\_PROF\_FJX sets vertical profiles for a given
!  latitude and longitude.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE SET_PROF_FJX( YLAT,      MONTH,  DAY,     T_CTM,  P_CTM,    &
                       CLDOD,     DSTOD,  AEROD,   O3_CTM, O3_TOMS,  &
                       AERCOL,    T_CLIM, O3_CLIM, Z_CLIM, AIR_CLIM, &
                       Input_Opt, State_Grid,      State_Chm )
!
! !USES:
!
    USE CMN_FJX_Mod,     ONLY : L_, L1_, A_, ZZHT
    USE CMN_SIZE_Mod,    ONLY : NAER, NRH, NDUST
    USE Input_Opt_Mod,   ONLY : OptInput
    USE PhysConstants,   ONLY : AIRMW, AVO, g0, BOLTZ
    USE State_Chm_Mod,   ONLY : ChmState
    USE State_Grid_Mod,  ONLY : GrdState
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
    TYPE(ChmState), INTENT(IN) :: State_Chm         ! Chemistry State object
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
       OREF2(I) = State_Chm%Phot%OREF(I,L,M)
       TREF2(I) = State_Chm%Phot%TREF(I,L,M)
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

  END SUBROUTINE SET_PROF_FJX
!EOC
END MODULE FJX_INTERFACE_MOD
