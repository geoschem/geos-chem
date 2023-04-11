#ifdef CLOUDJ
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: cldj_interface_mod.F90
!
! !DESCRIPTION: Module CLDJ\_INTERFACE\_MOD contains routines and variables
!  for interfacing with the Cloud-J scheme (Prather et al) that calculates
!  photolysis rates.
!\\
!\\
! !INTERFACE:
!
MODULE CLDJ_INTERFACE_MOD
!
! !USES:
!
  USE PRECISION_MOD
#if defined( MODEL_CESM ) && defined( SPMD )
  USE MPISHORTHAND
  USE SPMD_UTILS
#endif

  IMPLICIT NONE

  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC  :: INIT_CLOUDJ
  PUBLIC  :: RUN_CLOUDJ
!
! !PRIVATE MEMBER FUNCTIONS:
!
  PRIVATE :: SET_PROF_CLOUDJ
  PRIVATE :: SOLAR_JX ! Copy of fjx_mod SOLAR_JX pending looking more at cldj
!
! !REVISION HISTORY:
!  14 Dec 2022 - E. Lundgren - initial version, adapted from fast_jx_mod
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
! !IROUTINE: int_cloudj
!
! !DESCRIPTION: Subroutine INIT\_CLOUDJ initializes Cloud-J variables.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE INIT_CLOUDJ( Input_Opt, State_Diag, State_Chm, RC )
!
! !USES:
!

! ewl: Use, inputs/outputs, and local vars could be slimmed down
    USE Charpak_Mod,    ONLY : CSTRIP
    ! ewl: if these are in cloud-j, why do I need to pass them???
    USE Cldj_Cmn_Mod,   ONLY : JVN_, NJX, NRATJ, W_, WL
    USE Cldj_Cmn_Mod,   ONLY : TITLEJX, JLABEL, JFACTA, RNAMES
    USE Cldj_Cmn_Mod,   ONLY : JIND, BRANCH     ! ewl debugging
    USE Cldj_Init_Mod,  ONLY : Init_CldJ
    USE ErrCode_Mod
    USE Input_Opt_Mod,  ONLY : OptInput
    USE inquireMod,     ONLY : findFreeLUN
    USE State_Chm_Mod,  ONLY : ChmState
    USE State_Diag_Mod, ONLY : DgnState
#if defined( MODEL_CESM )
    USE UNITS,          ONLY : freeUnit
#endif

!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)     :: Input_Opt   ! Input Options object
    TYPE(DgnState), INTENT(IN)     :: State_Diag  ! Diagnostics State object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT)  :: State_Chm   ! Chemistry State object

!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)    :: RC          ! Success or failure?
!
! !REVISION HISTORY:
!  14 Dec 2022 - E. Lundgren - initial version, adapted from fast_jx_mod
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    LOGICAL            :: notDryRun
    INTEGER            :: JXUNIT, J, K, NJXX, PhotoId
    REAL(fp)           :: ND64MULT

    ! Strings
    CHARACTER(LEN=6)   :: TITLEJXX(JVN_)
    CHARACTER(LEN=50 ) :: TEXT
    CHARACTER(LEN=255) :: DATA_DIR
    CHARACTER(LEN=255) :: FILENAME
    CHARACTER(LEN=255) :: ErrMsg, ThisLoc

    !=================================================================
    ! INIT_CLOUDJ begins here!
    !=================================================================

    ! Initialize
    RC          = GC_SUCCESS
    notDryRun   = ( .not. Input_Opt%DryRun )
    ErrMsg      = ''
    ThisLoc     = ' -> at Init_CloudJ (in module GeosCore/cldj_interface_mod.F90)'

    ! Skip these opterations when running in dry-run mode
    IF ( notDryRun ) THEN

       ! Print info
       IF ( Input_Opt%amIRoot ) THEN
          write(6,*) ' Initializing Cloud-J'

          if (W_.ne.8 .and. W_.ne.12 .and. W_.ne.18) then
             ErrMsg =  ' INIT_CLOUDJ: invalid no. wavelengths'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          endif
       ENDIF

    ENDIF

    ! Define data directory for FAST-JX input
    DATA_DIR = TRIM( Input_Opt%CLOUDJ_DIR )

    ! Initialize Cloud-J. Includes reading input data files
    ! FJX_spec.dat (RD_XXX), FJX_scat-aer.dat (RD_MIE), and 
    ! FJX_j2j.dat (RD_JS_JX)
    CALL Init_CldJ(DATA_DIR,TITLEJXX,JVN_,NJXX)

    ! Store # of photolysis reactions in State_Chm object for easy reference
    State_Chm%Phot%nPhotRxns = NRatJ

  END SUBROUTINE INIT_CLOUDJ
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: run_cloudj
!
! !DESCRIPTION: Subroutine RUN\_CLOUDJ loops over longitude and latitude, and
!  calls CLOUD\_JX to compute J-Values for each column at every chemistry
!  time-step.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Run_CloudJ( WLAOD, Input_Opt, State_Chm, State_Diag, &
                         State_Grid, State_Met, RC )
!
! !USES:
!
    ! ewl: if these are in cloud-j, why do we need them here??
    USE Cldj_Cmn_Mod,       ONLY : A_, L_, L1_, W_, S_
    USE Cldj_Cmn_Mod,       ONLY : JVN_, JXL_, JXL1_, AN_, NQD_, W_r
    USE Cldj_Cmn_Mod,       ONLY : JIND, JFACTA, FL 
    USE Cld_Sub_Mod,        ONLY : Cloud_JX
    USE CMN_SIZE_MOD,       ONLY : NDUST, NRH, NAER, NRHAER
    USE ErrCode_Mod
    USE ERROR_MOD,          ONLY : ERROR_STOP, ALLOC_ERR
    USE ERROR_MOD,          ONLY : DEBUG_MSG
    USE Input_Opt_Mod,      ONLY : OptInput
    USE PhysConstants,      ONLY : AVO
    USE State_Chm_Mod,      ONLY : ChmState
    USE State_Chm_Mod,      ONLY : Ind_
    USE State_Diag_Mod,     ONLY : DgnState
    USE State_Grid_Mod,     ONLY : GrdState
    USE State_Met_Mod,      ONLY : MetState
    USE TIME_MOD,           ONLY : GET_MONTH, GET_DAY, GET_DAY_OF_YEAR
    USE TIME_MOD,           ONLY : GET_TAU,   GET_YEAR
    USE TOMS_MOD,           ONLY : GET_OVERHEAD_O3

    IMPLICIT NONE
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
! !REVISION HISTORY:
!  14 Dec 2022 - E. Lundgren - initial version, adapted from fast_jx_mod
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
!ewl: might be able to clean up use statements, inputs/outpus, and local vars
    INTEGER, SAVE :: LASTMONTH = -1
    INTEGER       :: NLON, NLAT, DAY,  MONTH, DAY_OF_YR, L, N, J
    INTEGER       :: IOPT, LCHEM, IWV1000
    REAL(fp)      :: U0, PRES, YLAT,  O3_TOMS, SZA, SOLF
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

    ! More species IDs (added by ewl)
    INTEGER, SAVE :: id_SO4
    INTEGER, SAVE :: id_BCPI
    INTEGER, SAVE :: id_BCPO
    INTEGER, SAVE :: id_OCPI
    INTEGER, SAVE :: id_OCPO
    INTEGER, SAVE :: id_SALA
    INTEGER, SAVE :: id_SALC

    ! Strings
    CHARACTER(LEN=255) :: ErrMsg, ThisLoc

    ! ewl: input to set_prof, but I question whether we need it
    REAL(fp) :: ODCLOUD_COL(L_)

    ! ewl: output from photo_jx for use in diagnostics
    REAL(fp) :: FJBOT(W_)
    REAL(fp) :: FSBOT(W_)
    REAL(fp) :: FLXD(JXL1_,W_)
    REAL(fp) :: FJFLX(JXL_,W_)

    ! UVFlux* diagnostics
    REAL(fp) :: FDIRECT (JXL1_)
    REAL(fp) :: FDIFFUSE(JXL1_)
    REAL(fp) :: UVX_CONST
    INTEGER  :: S, K

    ! ewl: Inputs to Cloud_JX
    integer  :: IRAN
    real(8)  :: CLDCOR !v7.7
    real(8)  :: RFL(5,W_+W_r) ! RFL already exists
    logical  :: LPRTJ
    real(8)  :: PPP(L1_+1)
    real(8)  :: ZZZ(L1_+1)
    real(8)  :: TTT(L1_)
    real(8)  :: HHH(L1_)
    real(8)  :: DDD(L1_)
    real(8)  :: RRR(L1_)
    real(8)  :: OOO(L1_)
    real(8)  :: CCC(L1_)
    real(8)  :: LWP(L1_)
    real(8)  :: IWP(L1_)
    real(8)  :: REFFL(L1_)
    real(8)  :: REFFI(L1_)
    real(8)  :: AERSP(L1_,AN_)
    integer  :: NDXAER(L1_,AN_)
    real(8)  :: CLDF(L1_)
    integer  :: CLDIW(L1_)

    ! ewl: Outputs of Cloud_JX
    !real(8)  :: VALJXX(L1U-1,NJXU) ! VALJXX already exists
    real(8)  :: SKPERD(S_+2,L1_)
    real(8)  :: SWMSQ(6)
    real(8)  :: OD18(L1_)
    real(8)  :: WTQCA(NQD_)
    integer  :: NICA
    integer  :: JCOUNT
    logical  :: LDARK

    ! ewl: added to get concentrations
    integer  :: A
    real(8)  :: MW_kg
    real(8)  :: BoxHt

    !=================================================================
    ! Run_CloudJ begins here!
    !=================================================================

    ! Initialize
    RC        = GC_SUCCESS
    ErrMsg    = ''
    ThisLoc   = ' -> at Run_CloudJ (in module GeosCore/cldj_interface_mod.F90)'
    prtDebug  = Input_Opt%Verbose

    ! Get wavelength index for 1000 hm
    IWV1000 = State_Chm%Phot%IWV1000

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

    ! ewl: set NDXAER to MIEDX duplicated for all levels. If this works
    ! will want to store this elsewhere. Don't want to do this computation
    ! every timestep.
    NDXAER(:,:) = 0.d0
    DO N = 1, AN_
    DO L = 1, L1_
       NDXAER(L,N) = State_Chm%Phot%MIEDX(N)
    ENDDO
    ENDDO

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

       ! Get other species
       id_SO4  = Ind_('SO4')
       id_BCPI = Ind_('BCPI')
       id_BCPO = Ind_('BCPO')
       id_OCPI = Ind_('OCPI')
       id_OCPO = Ind_('OCPO')
       id_SALA = Ind_('SALA')
       id_SALC = Ind_('SALC')

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
    !$OMP PRIVATE( LPRTJ, HHH, RRR, CCC, LWP, IWP, REFFL, REFFI ) &
    !$OMP PRIVATE( CLDF, CLDCOR, CLDIW, AERSP, IRAN, SKPERD     ) &
    !$OMP PRIVATE( SWMSQ, OD18, NICA, JCOUNT, LDARK, WTQCA      ) &
    !$OMP PRIVATE( A, MW_kg, BoxHt                              ) &
    !$OMP PRIVATE( FDIRECT, FDIFFUSE, UVX_CONST, K, S           ) &
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
       RFL(1,:) = State_Met%UVALBEDO(NLON,NLAT)
       RFL(2,:) = State_Met%UVALBEDO(NLON,NLAT)
       RFL(3,:) = State_Met%UVALBEDO(NLON,NLAT)
       RFL(4,:) = State_Met%UVALBEDO(NLON,NLAT)
       RFL(5,:) = State_Met%UVALBEDO(NLON,NLAT)

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
          IOPT = ( (N-1) * NRH ) + State_Chm%Phot%IRHARR(NLON,NLAT,L)
          OPTAER(L,IOPT) = State_Chm%Phot%ODAER(NLON,NLAT,L,IWV1000,N)
       ENDDO
       ENDDO
       DO N = 1, NDUST
       DO L = 1, State_Grid%NZ
          OPTDUST(L,N) = State_Chm%Phot%ODMDUST(NLON,NLAT,L,IWV1000,N)
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
          State_Diag%EXTRALNLEVS(ILON,ILAT) = 0.0
       ENDIF
       IF ( State_Diag%Archive_EXTRALNITER ) THEN
          State_Diag%EXTRALNITER(ILON,ILAT) = 0.0
       ENDIF
#endif

       if (State_Grid%NZ+1 .gt. JXL1_) then
          CALL ERROR_STOP('PHOTO_JX: not enough levels in JX', &
                          'cldj_interface_mod.F90' )
       endif

       ! Input conversion (SDE 03/29/13)
       ! Calculate solar zenith angle (degrees)
       CALL SOLAR_JX(DAY_OF_YR,U0,SZA,SOLF)
       
       ! check for dark conditions SZA > 98.0 deg => tan ht = 63 km or
       !                                 99.0                 80 km
       if (SZA .gt. 98.e+0_fp) cycle

! ewl: better to put these options in geoschem_config.yml?
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

       ! Copy cloud OD data to a variable array (ewl: why???)
       DO L=1,L_
          ODCLOUD_COL(L) = OPTD(L)
       ENDDO
       
       ! Use GEOS-Chem methodology to set vertical profiles of:
       ! Pressure      (PPJ)    [hPa]
       ! Temperature   (T_CLIm) [K]
       ! Path density  (DDJ)    [# molec/cm2]
       ! New methodology for:
       ! Ozone density (OOJ)    [# O3 molec/cm2]
       CALL SET_PROF_CLOUDJ (YLAT,        MONTH,     DAY,         &
                      T_CTM,       P_CTM,     OPTD,        &
                      OPTDUST,     OPTAER,    O3_CTM,      &
                      O3_TOMS,     AERX_COL,  T_CLIM,      &
                      O3_CLIM,     Z_CLIM,    AIR_CLIM,    &
                      Input_Opt,   State_Grid, State_Chm )

!       ! Call FAST-JX routines to compute J-values
!       CALL PHOTO_JX( Input_Opt%amIRoot, Input_Opt%DryRun,          &
!                      U0,        RFL,        SZA,       SOLF,       &
!                      P_CTM,     T_CTM,      AOD999,    NLON,       &
!                      NLAT,      AERX_COL,   T_CLIM,    O3_CLIM,    &
!                      Z_CLIM,    AIR_CLIM,   State_Grid%maxChemLev, &
!                      VALJXX,     FSBOT,     FJBOT,     FLXD,       &
!                      FJFLX                                        )

       ! Call Cloud-JX to compute J-values

       ! Inputs
       LPRTJ = .false. ! Limit Cloud-J debug prints to 1 grid cell
       IF ( NLON .eq. 20 .and. NLAT .eq. 20 ) LPRTJ = .true.
       HHH = 0.5_fp
       CCC = O3_CLIM ! Use O3 for CH4 just to get something going

       ! Need to verify if units okay
       RRR(1:State_Grid%NZ) = State_Met%RH(NLON,NLAT,1:State_Grid%NZ)
       RRR(1:State_Grid%NZ+1) = RRR(State_Grid%NZ)

       ! Need to verify if units okay
       LWP(1:State_Grid%NZ) = State_Met%CLDLIQ(NLON,NLAT,1:State_Grid%NZ)
       LWP(1:State_Grid%NZ+1) = LWP(State_Grid%NZ)

       ! Need to verify if units okay
       IWP(1:State_Grid%NZ) = State_Met%CLDICE(NLON,NLAT,1:State_Grid%NZ)
       IWP(1:State_Grid%NZ+1) = IWP(State_Grid%NZ)

       REFFL = 0.5d0 ! ? dummy
       REFFI = 0.3d0 ! ? dummy
       CLDF = State_Met%CLDFRC(NLON,NLAT)
       CLDCOR = 0.33
       CLDIW = 0 ! Will need to compute this from state_met
                 ! 0 = no cloud, 1 = water cld, 2 = ice cld, 3 = liq+ice cld
       IRAN = 0 ! Dummy since not using CLDFLAG=5

!-------------
       DO L=1,State_Grid%NZ

          ! Layer height [m]
          BoxHt = State_Met%BXHEIGHT(NLON,NLAT,L)

          ! AERSP is concentration in g/m2. Size is L1_, AN_. Will need to
          ! compute this for each lat/lon in this loop. Could be a separate
          ! subroutine, or in set_prof. Best put in set_prof.
          ! AN_ is # of separate aerosols per layer, = 37.
          ! L1_ is # of layer of edges, so 73. Make the top one the same as 72?
          ! Just do that for now, but need to check if okay.
          AERSP(:,:) = 0.d0
          
          ! cloud 1: black carbon absorber. MIEDX(1) = 3 in FJX_scate-aer.dat,
          ! which has all zero phase functions and single scattering albedo.
          AERSP(1:State_Grid%NZ,1) = 0.d0
          
          ! cloud 2: water cloud. MIEDX(2) = 10. Not sure what to put here.
          ! AERCOL for this uses input met cloud optical depth where T > 233
          ! Keeping zero for now.
          AERSP(1:State_Grid%NZ,2) = 0.d0
          
          ! cloud 3: irregular ice cloud. MIEDX(3) = 14. Not sure what to put here.
          ! AERCOL for this uses input met cloud optical depth where T <= 233
          ! Keeping zero for now
          AERSP(1:State_Grid%NZ,3) = 0.d0
          
          ! Mineral dust [kg/m3] -> [g/m2]
          AERSP(L,4)  = State_Chm%SOILDUST(NLON,NLAT,L,1)*BoxHt*1.d3 
          AERSP(L,5)  = State_Chm%SOILDUST(NLON,NLAT,L,2)*BoxHt*1.d3 
          AERSP(L,6)  = State_Chm%SOILDUST(NLON,NLAT,L,3)*BoxHt*1.d3 
          AERSP(L,7)  = State_Chm%SOILDUST(NLON,NLAT,L,4)*BoxHt*1.d3 
          AERSP(L,8)  = State_Chm%SOILDUST(NLON,NLAT,L,5)*BoxHt*1.d3 
          AERSP(L,9)  = State_Chm%SOILDUST(NLON,NLAT,L,6)*BoxHt*1.d3 
          AERSP(L,10) = State_Chm%SOILDUST(NLON,NLAT,L,7)*BoxHt*1.d3 
          
          ! Hygroscopic-growth aerosols, looped over relative humidity?
          ! for now I am using the same concentration for each RH bin
          ! and am adding hydrophilic and hydrophobic.
          ! [molec/cm3]. Will need to convert to g/m2.
          DO A=1,5

             ! Sulfate
             MW_kg = State_Chm%SpcData(id_SO4)%Info%MW_g * 1.e-3_fp
             AERSP(L,10+A) =                                              &
                 ( State_Chm%Species(id_SO4)%Conc(NLON,NLAT,L)            &
                 * 1e+6_fp / ( AVO / MW_kg ) / State_Met%AIRDEN(NLON,NLAT,L) )  &
                 * BoxHt
             ! Black carbon
             MW_kg = State_Chm%SpcData(id_BCPI)%Info%MW_g * 1.e-3_fp
             AERSP(L,10+NRHAER+A) =                                       &
                 ( ( State_Chm%Species(id_BCPI)%Conc(NLON,NLAT,L) +       &
                   State_Chm%Species(id_BCPO)%Conc(NLON,NLAT,L) )         &
                 * 1e+6_fp / ( AVO / MW_kg ) / State_Met%AIRDEN(NLON,NLAT,L) )  &
                 * BoxHt

             ! Organic carbon
             MW_kg = State_Chm%SpcData(id_OCPI)%Info%MW_g * 1.e-3_fp
             AERSP(L,10+NRHAER*2+A) =                                     &
                 ( ( State_Chm%Species(id_OCPI)%Conc(NLON,NLAT,L) +       &
                   State_Chm%Species(id_OCPO)%Conc(NLON,NLAT,L) )         &
                 * 1e+6_fp / ( AVO / MW_kg ) / State_Met%AIRDEN(NLON,NLAT,L) )  &
                 * BoxHt

             ! Seasalt (accum)
             MW_kg = State_Chm%SpcData(id_SALA)%Info%MW_g * 1.e-3_fp
             AERSP(L,10+NRHAER*3+A) =                                     &
                 ( State_Chm%Species(id_SALA)%Conc(NLON,NLAT,L)           &
                 * 1e+6_fp / ( AVO / MW_kg ) / State_Met%AIRDEN(NLON,NLAT,L) )  &
                 * BoxHt

             ! Seasalt (coarse)
             MW_kg = State_Chm%SpcData(id_SALC)%Info%MW_g * 1.e-3_fp
             AERSP(L,10+NRHAER*4+A) =                                     &
                 ( State_Chm%Species(id_SALC)%Conc(NLON,NLAT,L)           &
                 * 1e+6_fp / ( AVO / MW_kg ) / State_Met%AIRDEN(NLON,NLAT,L) )  &
                 * BoxHt

          ENDDO
          
          ! Use sulfate concentration for sulfate stratospheric aerosols
          AERSP(L,36) = AERSP(L,11)
          AERSP(L,37) = AERSP(L,11)

       ENDDO

       ! Set TOA equal to concentration in top level. Not sure if this is
       ! right. We set TOA optical depth to zero when passing to photo_jx
       ! in old fast-jx. 
       AERSP(State_Grid%NZ+1,:) = AERSP(State_Grid%NZ,:)

       ! Will need to convert units!!! What are units above?
       ! Loop over levels
       AERSP(State_Grid%NZ+1,:) = AERSP(State_Grid%NZ,:)

!-------------
       
       ! Outputs
       ! SKPERD
       ! SWMSQ
       ! OD18
       ! NICA
       ! JCOUNT
       ! LDARK
       ! WTQCA

       ! NOTE: there appears to be a parallelization error when computing
       ! the J-values diagnostic, either here or in fullchem_mod.

       ! ewl: deleted arguments CLDFLAG, NRANDO, and LNRG since globally
       ! set in cloud-j init via reading CJ77_inputs file
       CALL Cloud_JX( U0,       SZA,      RFL,      SOLF,     LPRTJ,     &
                      P_CTM,    Z_CLIM,   T_CLIM,   HHH,      AIR_CLIM,  &
                      RRR,      O3_CLIM,  CCC,      LWP,      IWP,       &
                      REFFL,    REFFI,    CLDF,     CLDCOR,   CLDIW,     &
                      AERSP,    NDXAER,   L1_,      AN_,      JVN_,      &
                      VALJXX,   SKPERD,   SWMSQ,    OD18,     IRAN,      &
                      NICA,     JCOUNT,   LDARK,    WTQCA               )

       ! Fill out common-block array of J-rates using PHOTO_JX output
       DO L=1,State_Grid%MaxChemLev
          DO J=1,State_Chm%Phot%nPhotRxns
             IF (JIND(J).gt.0) THEN
                State_Chm%Phot%ZPJ(L,J,NLON,NLAT) = VALJXX(L,JIND(J))*JFACTA(J)
             ELSE
                State_Chm%Phot%ZPJ(L,J,NLON,NLAT) = 0.e+0_fp
             ENDIF
          ENDDO
       ENDDO

       ! Set J-rates outside the chemgrid to zero
       IF (State_Grid%MaxChemLev.lt.L_) THEN
          DO L=State_Grid%MaxChemLev+1,L_
             DO J=1,State_Chm%Phot%nPhotRxns
                State_Chm%Phot%ZPJ(L,J,NLON,NLAT) = 0.e+0_fp
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
       
             ! ewl: this is messed up. FSBOT and FJBOT aren't set.

             ! Direct & diffuse fluxes at each level
             FDIRECT(1)  = FSBOT(K)                    ! surface
             FDIFFUSE(1) = FJBOT(K)                    ! surface
             DO L = 2, State_Grid%NZ
                FDIRECT(L) = FDIRECT(L-1) + FLXD(L-1,K)
                FDIFFUSE(L) = FJFLX(L-1,K)
             ENDDO
       
             ! Constant to multiply UV fluxes at each wavelength bin
             UVX_CONST = SOLF * FL(K) * State_Chm%Phot%UVXFACTOR(K)
       
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

  END SUBROUTINE Run_CloudJ
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: set_prof_cloudj
!
! !DESCRIPTION: Subroutine SET\_PROF\_CLOUDJ sets vertical profiles for a given
!  latitude and longitude.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE SET_PROF_CLOUDJ( YLAT,      MONTH,  DAY,     T_CTM,  P_CTM,    &
                       CLDOD,     DSTOD,  AEROD,   O3_CTM, O3_TOMS,  &
                       AERCOL,    T_CLIM, O3_CLIM, Z_CLIM, AIR_CLIM, &
                       Input_Opt, State_Grid, State_Chm )
!
! !USES:
!
    USE Cldj_Cmn_Mod,       ONLY : L_, L1_, A_, ZZHT
    USE CMN_SIZE_Mod,       ONLY : NAER, NRH, NDUST
    USE Input_Opt_Mod,      ONLY : OptInput
    USE PhysConstants,      ONLY : AIRMW, AVO, g0, BOLTZ
    USE State_Grid_Mod,     ONLY : GrdState
    USE State_Chm_Mod,      ONLY : ChmState
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
!  14 Dec 2022 - E. Lundgren   - Copied from GEOS-Chem fast_jx_mod
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

  END SUBROUTINE SET_PROF_CloudJ
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
    USE PhysConstants, ONLY : PI
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
END MODULE CLDJ_INTERFACE_MOD
#endif
