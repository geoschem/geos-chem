#if defined( ESMF_ )
!------------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 910.1 and      !
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: gigc_chemdr
!
! !DESCRIPTION: Module GC\_CHEMDR is the "bridge" between the ESMF interface
!  to the GEOS-5 GCM and the GEOS-Chem chemistry routines.
!\\
!\\
! !INTERFACE:
!
MODULE GIGC_ChemDr
!
! !USES:
!
  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: GIGC_Do_Chem
!
! !REMARKS:
!  This mo
!
! !REVISION HISTORY: 
!  15 Oct 2012 - M. Long     - Initial version
!  15 Oct 2012 - R. Yantosca - Added ProTeX headers, F90 indentation
!  16 Oct 2012 - R. Yantosca - Renamed GC_STATE to State_Chm
!  16 Oct 2012 - R. Yantosca - Renamed GC_MET to State_Met
!  22 Oct 2012 - R. Yantosca - Renamed to gigc_chemdr.F90
!EOP
!------------------------------------------------------------------------------
!BOC
CONTAINS
!EOC
!------------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 910.1 and      !
!          Harvard University Atmospheric Chemistry Modeling Group
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: gigc_do_chem
!
! !DESCRIPTION: Routine GIGC\_DO\_CHEM is the chemistry driver for the
!  Grid-Independent GEOS-Chem (aka "GIGC") model.  This routine is called by
!  the Run method from the ESMF interface to the GEOS-5 GCM.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE GIGC_Do_Chem( State_Chm, State_Met, am_I_Root, NI, &
                           NJ,        NL,        NCNST,     RC )
!
! !USES:
!
    USE CHEMISTRY_MOD,      ONLY : DO_CHEMISTRY
    USE CMN_DEP_MOD,        ONLY : FRCLND
    USE CMN_SIZE_MOD,       ONLY : DJSIZE, DISIZE, LLPAR, IIPAR, JJPAR
    USE COMODE_MOD,         ONLY : AIRDENS, CSPEC_FULL
    USE COMODE_LOOP_MOD
    USE DAO_MOD
    USE GIGC_State_Chm_Mod, ONLY : ChmState
    USE GIGC_State_Met_Mod, ONLY : MetState
    USE GIGC_Test_Utils
    USE GIGC_Chem_Utils
    USE GRID_MOD,           ONLY : AREA_M2, YEDGE, XEDGE, YMID, XMID
    USE LOGICAL_MOD
    USE PBL_MIX_MOD,        ONLY : PBL_TOP_L, PBL_TOP_M, INIT_PBL_MIX
    USE PRESSURE_MOD,       ONLY : EXTERNAL_PEDGE
    USE TRACER_MOD
    USE TRACERID_MOD     
    USE UVALBEDO_MOD,       ONLY : UVALBEDO
!
! !INPUT PARAMETERS:
!
    INTEGER,        INTENT(IN)    :: NI          ! # of longitudes
    INTEGER,        INTENT(IN)    :: NJ          ! # of latitudes
    INTEGER,        INTENT(IN)    :: NL          ! # of levels
    INTEGER,        INTENT(IN)    :: NCNST       ! # of constituents
    LOGICAL,        INTENT(IN)    :: am_I_Root   ! Are we on the root CPU
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
    TYPE(MetState), INTENT(INOUT) :: State_Met   ! Meteorology State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure
!
! !REMARKS:
!  NOTE: Eventually we will replace all met field arrays with the
!  values from the meteorology state.  For now we need to copy these
!  over. (bmy, 10/15/12)
!
! !REVISION HISTORY:
!  15 Oct 2012 - M. Long     - Initial version
!  15 Oct 2012 - R. Yantosca - Added ProTeX headers, F90 indentation
!  16 Oct 2012 - R. Yantosca - Renamed GC_STATE to State_Chm
!  16 Oct 2012 - R. Yantosca - Renamed GC_MET to State_Met
!  19 Oct 2012 - R. Yantosca - Now reference gigc_state_chm_mod.F90
!  19 Oct 2012 - R. Yantosca - Now reference gigc_state_met_mod.F90
!  22 Oct 2012 - R. Yantosca - Renamed to GIGC_Do_Chem
!  25 Oct 2012 - R. Yantosca - Now convert units of State_Chm%TRACERS from 
!                              [v/v] -> [kg] before calling G-C chemistry 
!                              (and from [kg] -> [v/v] after chemistry)
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!

! TESTING SECTION
!<><><><><><><><><><><><><><><><><><><><><><>
     !CALL GIGC_Dump_Config( am_I_Root )


!<><><><><><><><><><><><><><><><><><><><><><>
! TEMPORARY INITIALIZATION SECTION
      IIPAR      = NI
      JJPAR      = NJ
      LLPAR      = NL
      LLTROP     = NL
      LLTROP_FIX = NL

      RRATE      = 0.E0
      TRATE      = 0.E0

      !======================================================================
      ! Set 2-D variables
      !
      ! NOTE: This is a stopgap measure for testing.  Eventually we will
      ! carry the meteorology state down to all G-C routines.  in order to
      ! test, we need to populate the G-C module arrays from the met state.
      !======================================================================

      ! Met fields
      TO3                = State_Met%TO3       ! Total column O3 [DU]
      ALBD               = State_Met%ALBD      ! visible surface albedo [1]
      AREA_M2            = State_Met%AREA_M2   ! grid box surface area [cm2]
      CLDFRC             = State_Met%CLDFRC    ! column cloud fraction [1]
      FRCLND             = State_Met%FRCLND    ! Olson land fraction [1]
      GWETTOP            = State_Met%GWETTOP   ! top soil moisture [1]
      HFLUX              = State_Met%HFLUX     ! Sensible heat flux [W/m2]
      LWI                = State_Met%LWI       ! Land/water indices [1]
      PARDR              = State_Met%PARDR     ! Direct  PAR [W/m2]
      PARDF              = State_Met%PARDF     ! Diffuse PAR [W/m2]
      PBL_TOP_M          = State_Met%PBLH      ! PBL height [m]
      PRECON             = State_Met%PRECCON   ! Conv  precip @ ground [kg/m2/s]
      PREACC             = State_Met%PRECTOT   ! Total precip @ ground [kg/m2/s]
      RADSWG             = State_Met%RADSWG    ! Solar radiation @ ground [W/m2]
      TSKIN              = State_Met%TS        ! Sea surface temperature [K]
      SUNCOS(1:NI*NJ)    = State_Met%SUNCOS    ! Cosine of solar zenith angle
      SUNCOS_MID(1:NI*NJ)= State_Met%SUNCOS    ! Cosine of solar zenith angle
      TROPP              = State_Met%TROPP     ! Tropopause pressure [hPa]
      TS                 = State_Met%TS        ! Surface temperature [K]
      U10M               = State_Met%U10M      ! E/W wind speed @ 10m [m/s]
      USTAR              = State_Met%USTAR     ! Friction velocity [m/s]
      UVALBEDO           = State_Met%UVALBEDO  ! UV surface albedo [1]
      V10M               = State_Met%V10M      ! N/S wind speed @ 10m [m/s]
      Z0                 = State_Met%Z0        ! Roughness height [m]

      !======================================================================
      ! Set 3-D variables
      !
      ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      ! %%% NOTE: This is a stopgap measure for testing.  Eventually we   %%%
      ! %%% will carry the Meteorology State and Chemistry State objects  %%%
      ! %%% down to all G-C routines.   In order to continue testing the  %%%
      ! %%% GIGC without disrupting existing workflow, we must populate   %%%
      ! %%% G-C module arays from the Meteorology and Chemistry States.   %%%
      ! %%% (bmy, 10/25/12)                                               %%%
      ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      !======================================================================

      ! Met fields
      AD                 = State_Met%AD        ! Air mass [kg]
      AIRDEN             = State_Met%AIRDENS   ! Air density [kg/m3]
      AIRVOL             = State_Met%AIRVOL    ! Grid box volume [m3]
      BXHEIGHT           = State_Met%BXHEIGHT  ! Grid box height [m]
      CLDF               = State_Met%CLDF      ! 3-D cloud fraction
      CMFMC              = State_Met%CMFMC     ! Cloud mass flux [kg/m2/s]
      DELP               = State_Met%DELP      ! Pressure thickness [hPa]
      DQIDTMST           = State_Met%DQIDTMST  ! Ice tendency, mst [kg/kg/s]
      DQLDTMST           = State_Met%DQLDTMST  ! H2O tendency, mst [kg/kg/s]
      DQVDTMST           = State_Met%DQVDTMST  ! Vapor tendency, mst [kg/kg/s]
      DTRAIN             = State_Met%DTRAIN    ! Detrainment flux [kg/m2/s]
      EXTERNAL_PEDGE     = State_Met%PEDGE     ! Pressure @ level edges [hPa]
      MOISTQ             = State_Met%MOISTQ    ! Tendency in sp. hum [kg/kg/s]
      OPTD               = State_Met%OPTD      ! Visible optical depth [1]    
      QI                 = State_Met%QI        ! Ice mixing ratio [kg/kg]
      QL                 = State_Met%QL        ! Water mixing ratio [kg/kg]
      RH                 = State_Met%RH        ! Relative humidity [1]
      SPHU               = State_Met%SPHU      ! Specific humidity [kg/kg]
      T                  = State_Met%T         ! Temperature [K]
      TAUCLI             = State_Met%TAUCLI    ! Opt depth of ice clouds [1]
      TAUCLW             = State_Met%TAUCLW    ! Opt depth of h2o clouds [1]

      ! Constituents
      CSPEC_FULL         = State_Chm%Species   ! Chemical species

      !======================================================================
      ! Call the GEOS-Chem Chemistry routines
      !
      ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      ! %%% NOTE: The values of State_Chm%TRACERS are taken from the ESMF %%%
      ! %%% Internal State, which has units of [v/v].  We need to convert %%%
      ! %%% this to [kg] before passing to the GEOS-Chem chemistry.  The  %%%
      ! %%% GEOS-Chem code expects tracers to be in [kg] at the point     %%%
      ! %%% where the chemistry routines are called. (bmy, 10/25/12)      %%%
      ! %%%                                                               %%%
      ! %%% NOTE: For now we initialize the GEOS-Chem tracer array STT    %%%
      ! %%% with State_Chm%TRACERS within DO_CHEMISTRY (bmy, 10/25/12)    %%%
      ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      !======================================================================

      !### Debug, print values in v/v before chem
      !IF ( am_I_Root ) THEN
      !   WRITE(6,*) '##### GIGC_Do_Chem, TRC_OX before chem [v/v]'
      !   WRITE(6,*) State_Chm%Tracers(1,1,:,2)
      !ENDIF

      ! If we are doing chemistry
      IF ( LCHEM ) THEN

         ! The tracer concentrations in State_Chm%TRACERS state have units 
         ! of [v/v] (since these come from the ESMF internal state).  We need
         ! to convert these to [kg] before calling the GEOS-Chem chemistry.
         CALL Convert_Units( 2, N_TRACERS, TCVV, AD, State_Chm%Tracers )

         ! Call the GEOS-Chem chemistry routines
         CALL Do_Chemistry( am_I_Root, NI, NJ, NL, State_Chm, State_Met, RC )

         ! Convert the tracer concentrations in State_Chm%TRACERS back to
         ! [v/v] so that they can be stored in the ESMF internal state
         ! for the next chemistry timestep.
         CALL Convert_Units( 1, N_TRACERS, TCVV, AD, State_Chm%Tracers )

      ENDIF

      !======================================================================
      ! Reset 3-D variables
      !
      ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      ! %%% NOTE: This is a stopgap measure for testing.  Eventually we   %%%
      ! %%% will carry the Meteorology State and Chemistry State objects  %%%
      ! %%% down to all G-C routines.   In order to continue testing the  %%%
      ! %%% GIGC without disrupting existing workflow, we must populate   %%%
      ! %%% G-C module arays from the Meteorology and Chemistry States.   %%%
      ! %%% (bmy, 10/25/12)                                               %%%
      ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      !======================================================================

      ! Save chemistry output for next timestep
      State_Chm%Species = CSPEC_FULL

      !### Debug
      !IF ( am_I_Root ) THEN
      !   WRITE(6,*) '##### GIGC_Do_Chem, TRC_OX after chem'
      !   WRITE(6,*) State_Chm%Tracers(1,1,:,2)
      !ENDIF

    END SUBROUTINE GIGC_Do_Chem
!EOC
!-----------------------------------------------------------------------------
! NOTE: Preserve version with BCC-specific stuff.  We should probably split
! the BCC code off into another subroutine. (bmy, 10/15/12)
!
!    SUBROUTINE DO_GC_CHEM(State_Chm, State_Met, am_I_Root, NI, NJ, NL, NCNST)
!
!!                    MMR, SPECHUM , T_IN      ,                            &
!!                    PMID_IN     , PDEL_IN    , PINT_IN,                   &
!!                    ALBEDO_IN   , CLTOT_IN, LANDFRAC, PBLHT         ,     &
!!                    PRECL_IN , PRECC_IN, PTROP_IN  , TREFHT         ,     &
!!                    SST_IN   , U_10M   , V_10M     , USTAR_IN       ,     &
!!                    Z0_IN    , CLD     , MC        , DU2            ,     &
!!                    CICEWP   , CLIQWP  , RELHUM    , SHFLX          ,     &
!!                    RHO      , COSLAT  , TAU       ,                      &
!!                    SOLS     , SOLSD   , FSDS      , UV_ALBEDO      ,     &
!!                    DZ       , TEND_CIW, TEND_CLW  , TEND_W         ,     &
!!                    MOISTQ_IN, TS_IN   , GWETTOP_IN,          &
!!                    DLAT, DLON, LAT, LON, ZMID                      ,     &
!!                    NI, NJ, NL, NCNST, State_Chm)!, State_Met)
!
!
!!      USE PPGRID,          ONLY: PCOLS, PVER
!!      USE PHYSICS_TYPES,   ONLY: PHYSICS_STATE, PHYSICS_TEND, PHYSICS_PTEND
!      USE GC_TYPE_MOD,     ONLY: State_Met_LOCAL
!      USE GC_TYPE2_MOD,    ONLY: CHEMSTATE, EXT_STRATOH, EXT_SJVALUE,     &
!                                 EXT_COPROD, EXT_COLOSS
!      USE CHEMISTRY_MOD,   ONLY: DO_CHEMISTRY
!!      USE CHEMISTRY,       ONLY: NCNST
!      USE DAO_MOD,         ONLY: AIRQNT
!!      USE GC_INITIALIZATION_MOD,ONLY: State_Met, State_Chm
!      USE PBL_MIX_MOD,     ONLY: PBL_TOP_L, PBL_TOP_M, INIT_PBL_MIX
!
!! TEMPORARY USE
!      USE GRID_MOD,     ONLY : AREA_M2, YEDGE, XEDGE, YMID, XMID
!      USE CMN_SIZE_MOD, ONLY : DJSIZE, DISIZE, LLPAR, IIPAR, JJPAR
!      ! 3-D
!      USE DAO_MOD,      ONLY : TO3, OPTD, CLDFRC, AIRVOL, BXHEIGHT, &
!                               CLDF, CMFMC, DQIDTMST, DQLDTMST,     &
!                               DQVDTMST, DTRAIN, MOISTQ, OPTDEP,    &
!                               DELP, RH, SPHU, T, TAUCLI, TAUCLW
!      ! 2-D
!      USE DAO_MOD,      ONLY : ALBD, CLDFRC, GWETTOP, HFLUX, LWI,   &
!                               PARDR, PARDF, PRECON, PREACC,        &
!                               RADSWG, TSKIN, SUNCOS, TROPP, TS,    &
!                               U10M, V10M, Z0, USTAR, AD, INIT_DAO, &
!                               SUNCOS_MID, AIRDEN
!      USE CMN_DEP_MOD,  ONLY : FRCLND
!      USE UVALBEDO_MOD, ONLY : UVALBEDO
!      USE COMODE_LOOP_MOD
!      USE COMODE_MOD,   ONLY : AIRDENS, CSPEC_FULL
!      USE TRACER_MOD,   ONLY : TCVV
!      USE TRACERID_MOD                !, ONLY : IDTNOX
!      
!      USE PRESSURE_MOD, ONLY : EXTERNAL_PEDGE !PEDGE, PMID
!
!      USE GC_TEST_UTILS
!      
!      IMPLICIT NONE
!      
!      TYPE(CHEMSTATE),    INTENT(INOUT) :: State_Chm
!      TYPE(State_Met_LOCAL), INTENT(INOUT) :: State_Met
!      
!      INTEGER, INTENT(IN) :: NI, NJ, NL, NCNST
!      LOGICAL, INTENT(IN) :: am_I_Root
!
!      !----------------------------------------
!      ! LINK DATA FIELDS FROM BCC TO GEOS-CHEM
!      ! THESE ARE CALCULATED ONLINE
!      !----------------------------------------
!      ! TRACER MASS MIXING-RATIO ARRAY
!      REAL*8 :: MMR(NI, NJ, NL, NCNST)
!      ! ALBEDO
!      REAL*8 :: ALBEDO_IN(NI, NJ)  ! ALBEDO
!      ! COLUMN-INTEGRATED CLOUD FRACTION
!      REAL*8 :: CLTOT_IN(NI, NJ)
!      ! SPECIFIC HUMIDITY
!      REAL*8 :: SPECHUM(NI, NJ, NL)
!      ! LAND FRACTION
!      REAL*8 :: LANDFRAC(NI, NJ)
!      ! LAND-WATER-ICE FLAGS
!      
!      ! PBL HEIGHT (M)
!      REAL*8 :: PBLHT(NI, NJ)
!      ! CONVECTIVE PRECIP
!      REAL*8 :: PRECC_IN(NI, NJ)
!      ! LS PRECIP
!      REAL*8 :: PRECL_IN(NI, NJ)
!      ! AN PRECIP
!!!!!      REAL*8 :: PTROP(:) ! IS THIS CORRECT?
!      ! TROPOPAUSE PRESSURE
!      REAL*8 :: PTROP_IN(NI, NJ)
!      ! SURFACE TEMPERATURE (2-METER AIR TEMPERATURE??)
!      REAL*8 :: TREFHT(NI, NJ)
!      ! RADIATION @ GROUND (NET, LONG OR SHORTWAVE?)
!      REAL*8 :: FSDS(NI, NJ)
!      ! SEA-SURFACE TEMPERATURE
!      REAL*8 :: SST_IN(NI, NJ)
!      ! 10-METER MERIDIONAL WIND-SPEED
!      REAL*8 :: U_10M(NI, NJ)
!      ! 10-METER ZONAL WIND-SPEED
!      REAL*8 :: V_10M(NI, NJ)
!      ! FRICTION VELOCITY (U*)
!      REAL*8 :: USTAR_IN(NI, NJ)
!      ! ROUGHNESS HEIGHT
!      REAL*8 :: Z0_IN(NI, NJ)
!      ! CLOUD FRACTION
!      REAL*8 :: CLD(NI, NJ, NL)
!      ! CONVECTIVE (CLOUD) MASS FLUX
!      REAL*8 :: MC(NI, NJ, NL)
!      ! DETRAINMENT FLUX (CONVECTIVE: SHALLOW & DEEP?)
!      REAL*8 :: DU2(NI, NJ, NL)
!      ! VISIBLE (ICE) OPTICAL DEPTH
!      REAL*8 :: CICEWP(NI, NJ, NL)
!      ! VISIBLE (LIQUID WATER) OPTICAL DEPTH
!      REAL*8 :: CLIQWP(NI, NJ, NL)
!      ! RELATIVE HUMIDITY
!      REAL*8 :: RELHUM(NI, NJ, NL)
!      ! SENSIBLE HEAT FLUX
!      REAL*8 :: SHFLX(NI, NJ)
!      ! AIR TEMPERATURE
!      REAL*8 :: T_IN(NI, NJ, NL)
!      ! AIR DENSITY
!      REAL*8 :: RHO(NI, NJ, NL)
!      ! COSINE OF SOLAR ZENITH-ANGLE
!      REAL*8 :: COSLAT(NI*NJ)
!      ! MOIST (Q) TENDENCY
!      REAL*8 :: MOISTQ_IN(NI, NJ, NL)
!      ! MODEL LAYER PRESSURE-THICKNESS (PDEL)
!      REAL*8 :: PDEL_IN(NI, NJ, NL)
!      ! AIR PRESSURE (GRIDBOX EDGES)
!      REAL*8 :: PINT_IN(NI, NJ, NL+1)
!      ! MID-LEVEL PERSSURES
!      REAL*8 :: PMID_IN(NI, NJ, NL)
!      ! TOTAL (COLUMN INTEGRATED) OVERHEAD OZONE COLUMN
!!      REAL*8 :: 
!      ! GRID-BOX AREA
!!      REAL*8 :: AREA_M2(NI, NJ)
!      ! TOTAL OPTICAL DEPTH
!      REAL*8 :: TAU(NI, NJ, NL)
!      ! DIRECT PAR (THIS MAY NOT BE THE CORRECT VALUE)
!      REAL*8 :: SOLS(NI, NJ)
!      ! DIFFUSE PAR
!      REAL*8 :: SOLSD(NI, NJ)
!      ! UV ALBEDO
!      REAL*8 :: UV_ALBEDO(NI, NJ)
!      REAL*8 :: DZ(NI, NJ, NL)
!      REAL*8 :: TEND_CIW(NI, NJ, NL)
!      REAL*8 :: TEND_CLW(NI, NJ, NL)
!      REAL*8 :: TEND_W(NI, NJ, NL)
!      REAL*8 :: TS_IN(NI, NJ)
!      REAL*8 :: LAT(NI, NJ), LON(NI, NJ)
!      REAL*8 :: GWETTOP_IN(NI, NJ)
!
!      REAL*8 :: ZMID(NI, NJ, NL)
!      REAL*8 :: DLAT, DLON
!
!
!      ! TEMPORARY ARRAYS - DECLARED AND INITIALIZED HERE
!      ! TO FACILITATE TEST RUNS. THEIR VALUES WILL NEED TO 
!      ! BE PROPERLY DEALT WITH!
!
!
!      !-------------------------------------------------
!      ! THE FOLLOWING FIELDS NEED TO BE READ FROM FILE.
!      ! THE WILL NEED TO BE DISTRIBUTED WITH A
!      ! "SCATTER_FIELD_TO_CHUNK" ROUTINE CALL.
!      !-------------------------------------------------
!      ! IREG - NUMBER OF OLSEN LAND TYPES
!      ! ILAND (1-IREG) OLSEN LAND TYPE INDEX
!      ! IUSE  (1-IREG) FRACTION OF OLSEN LAND TYPE IN BOX
!      ! LAI (1-IREG) LEAF-AREA INDEX FOR LAND TYPE
!      ! STRATOSPHERIC OH
!      ! STRATOSPHERIC PAN J-VALUES
!      ! STRATOSPHERIC H2O2 J-VALUES
!      ! STRATOSEPHRIC ACETONE J-VALUES
!      ! STRATOSPHERIC KETONES J-VALUES
!      ! STRAT. ACETALDEHYDE J-VALUES
!      ! STRAT. ALDEHYDE(>C3) J-VALUES
!      ! STRAT. METHYLVINYLKETONE J-VALUES
!      ! STRAT. METHACROLEIN J-VALUES
!      ! STRAT. ALKYLNITRATE(>C3) J-VALUES
!      ! STRAT. HCHO J-VALUES
!      ! STRAT. N2O5 J-VALUES
!      ! STRAT. HNO4 J-VALUES
!      ! STRAT. CH3OOH J-VALUES
!      ! STRAT. CO PRODUCTION RATES
!      ! STRAT. CO LOSS RATES
!
!! TESTING SECTION
!!<><><><><><><><><><><><><><><><><><><><><><>
!      CALL DUMP_GC_CONFIG(am_I_Root)
!
!
!
!!<><><><><><><><><><><><><><><><><><><><><><>
!! TEMPORARY INITIALIZATION SECTION
!      IIPAR = NI
!      JJPAR = NJ
!      LLPAR = NL
!
!!      CALL INIT_DAO
!!<><><><><><><><><><><><><><><><><><><><><><>
!
!!      CALL BCC_PASS_TO_GC(State_Met, State_Chm)
!
!!      DISIZE = DLON*(360.D0/(2.*3.1416)) ! RADIANS TO DEGREES
!!      DJSIZE = DLAT*(360.D0/(2.*3.1416)) ! RADIANS TO DEGREES
!
!      ! THE FOLLOWING ARE SHORT-TERM FIXES
!!      YEDGE(1:PCOLS) = LAT
!!      YEDGE(PCOLS+1) = LAT(1)
!!      XEDGE(1)       = LON(1)
!!      XEDGE(2)       = LON(2)
!
!!      YMID(:PCOLS)   = LAT
!!      XMID(1)        = LON(1)
!
!      RRATE = 0.E0
!      TRATE = 0.E0
!
!      ! ARBITRARY VALUES PASSED TO SCHEM (STRATPSPHERIC CHEMISTRY)
!      ! THEY ARE GENERALLY READ FROM FILES. THIS SHOULD BE DONE.
!      EXT_STRATOH = 1.E-13
!      EXT_SJVALUE = 1.E-13
!      EXT_COPROD  = 1.E-15
!      EXT_COLOSS  = 0.9E-15
!
!      ! SET 2-D VARS
!      TO3 = State_Met%TO3    ! TOTAL COLUMN O3 BURDEN. SIMPLY A GUESS. SHOULD BE IN DOBSON UNITS
!      ALBD(:,:)    = State_Met%ALBD !ALBEDO_IN   ! VISIBLE SURFACE ALBEDO [UNITLESS]
!      AREA_M2      = State_Met%AREA_M2!AREA_M2     ! GRID BOX SURFACE AREA [CM2]
!      CLDFRC(:,:)  = State_Met%CLDFRC!CLTOT_IN    ! COLUMN CLOUD FRACTION [UNITLESS]
!      FRCLND(:,:)  = State_Met%FRCLND!LANDFRAC    ! OLSON LAND FRACTION [UNITLESS]
!      GWETTOP(:,:) = State_Met%GWETTOP!GWETTOP_IN  ! TOP SOIL MOISTURE [UNITLESS]
!      HFLUX(:,:)   = State_Met%HFLUX!SHFLX       ! SENSIBLE HEAT FLUX [W/M2]
!      LWI          = State_Met%LWI           ! LAND/WATER INDICES [UNITLESS]
!      PARDR(:,:)   = State_Met%PARDR!SOLS        ! DIRECT  PHOTSYN ACTIVE RAD [W/M2]
!      PARDF(:,:)   = State_Met%PARDF!SOLSD       ! DIFFUSE PHOTSYN ACTIVE RAD [W/M2]
!      PBL_TOP_M(:,:)=State_Met%PBLH!PBLHT       ! PBL HEIGHT [M]
!      PRECON(:,:)  = State_Met%PRECCON!PRECC_IN    ! CONV  PRECIP @ GROUND [KG/M2/S]
!      PREACC(:,:)  = State_Met%PRECTOT!PRECL_IN+PRECC_IN    ! TOTAL PRECIP @ GROUND [KG/M2/S]
!      RADSWG(:,:)  = State_Met%RADSWG!SHFLX       ! SOLAR RADIATION @ GROUND [W/M2]
!      TSKIN(:,:)   = State_Met%TS!TS_IN       ! SEA SURFACE TEMPERATURE [K]
!      SUNCOS(1:NI*NJ)    = State_Met%SUNCOS!COSLAT      ! COSINE OF SOLAR ZENITH ANGLE
!      SUNCOS_MID(1:NI*NJ)    = State_Met%SUNCOS!COSLAT      ! COSINE OF SOLAR ZENITH ANGLE
!      TROPP(:,:)   = State_Met%TROPP!PTROP_IN    ! TROPOPAUSE PRESSURE [HPA]
!      TS(:,:)      = State_Met%TS!TREFHT      ! SURFACE TEMPERATURE [K]
!      U10M(:,:)    = State_Met%U10M!U_10M       ! E/W WIND SPEED @ 10M HEIGHT [M/S]
!      USTAR(:,:)   = State_Met%USTAR!USTAR_IN    ! FRICTION VELOCITY [M/S]
!      UVALBEDO(:,:) = State_Met%UVALBEDO!UV_ALBEDO  ! UV SURFACE ALBEDO [UNITLESS]
!      V10M(:,:)    = State_Met%V10M!V_10M       ! N/S WIND SPEED @ 10M HEIGHT [M/S]
!      Z0(:,:)      = State_Met%Z0!Z0_IN       ! ROUGHNESS HEIGHT
!
!! SET 3-D VARS      
!      AD = State_Met%AD        ! AIR MASS [KG]
!!      State_Chm%TRACERS(:,:,:,:) = 1.E-2!MMR(:,:,:,1:NCNST)!(:,:PVER,:)
!      AIRDEN = State_Met%AIRDENS!reshape(State_Met%AIRDENS,(/NI*NJ*NL/)) !RESHAPE(RHO,(/NI*NJ*NL/)) ! AIR DENSITY [KG/M3]
!
!!XXXXXXXXXXXXXXXXXXX
!      AIRVOL   = State_Met%AIRVOL ! GRID BOX VOLUME [M3]
!      BXHEIGHT(:,:,:) = State_Met%BXHEIGHT !DZ  ! GRID BOX HEIGHT [M]
!!XXXXXXXXXXXXXXXXXXX
!
!      CLDF(:,:,:)    = State_Met%CLDF!RESHAPE(CLD, (/NL,NI,NJ/))! 3-D CLOUD FRACTION [UNITLESS]
!      CMFMC(:,:,:)   = State_Met%CMFMC!MC     ! CLOUD MASS FLUX [KG/M2/S]
!
!!XXXXXXXXXXXXXXXXXXX
!      DQIDTMST = State_Met%DQIDTMST ! ICE TENDENCY, MST PROC [KG/KG/S]
!      DQLDTMST = State_Met%DQLDTMST ! H2O TENDENCY, MST PROC [KG/KG/S]
!      DQVDTMST = State_Met%DQVDTMST ! VAPOR TENDENCY, MST PROC [KG/KG/S]
!!XXXXXXXXXXXXXXXXXXX
!
!      DTRAIN(:,:,:)  = State_Met%DTRAIN !DU2   ! DETRAINMENT FLUX [KG/M2/S]
!
!!XXXXXXXXXXXXXXXXXXX
!      MOISTQ(:,:,:)  = State_Met%MOISTQ !RESHAPE(MOISTQ_IN,(/NL,NI,NJ/)) ! TENDENCY IN SP. HUMIDITY [KG/KG/S]
!      OPTD(:,:,:)    = State_Met%OPTD !RESHAPE(TAU      ,(/NL,NI,NJ/)) ! VISIBLE OPTICAL DEPTH [UNITLESS]
!!XXXXXXXXXXXXXXXXXXX
!
!      State_Met%PEDGE   = (State_Met%PEDGE(:,:,NL+1:1:-1))   ! PRESSURE @ LEVEL EDGES [PA]
!      EXTERNAL_PEDGE = State_Met%PEDGE
!!      State_Met%PMID(:,:,:)    =       ! PRESSURE @ LEVEL CENTERS [PA]
!      DELP(:,:,:)    = State_Met%DELP !RESHAPE(PDEL_IN  ,(/NL,NI,NJ/))      !
!      RH(:,:,:)      = State_Met%RH !RELHUM        ! RELATIVE HUMIDITY [UNITLESS]
!      SPHU(:,:,:)    = State_Met%SPHU !SPECHUM      ! SPECIFIC HUMIDITY [KG/KG]
!      T(:,:,:)       = State_Met%T !T_IN         ! TEMPERATURE [K]
!      TAUCLI(:,:,:)  = State_Met%TAUCLI !CICEWP  ! OPT DEPTH OF ICE CLOUDS [UNITLESS]
!      TAUCLW(:,:,:)  = State_Met%TAUCLW !CLIQWP  ! OPT DEPTH OF H2O CLOUDS [UNITLESS]
!
!
!      !CALL CHEMISTRY
!      CSPEC_FULL = State_Chm%CSPEC
!
!      CALL DO_CHEMISTRY(am_I_Root, NI, NJ, NL, State_Chm, State_Met)
!
!      State_Chm%CSPEC = CSPEC_FULL
!      !REINTERFACE WITH GCM: FROM GEOS-CHEM TO BCC
!!      CALL GC_PASS_TO_BCC(State_Chm)
!
!!      MMR(:,:,:,1:NCNST) = State_Chm%TRACERS(:,:,:,1:NCNST)!1.D-12
!
!      RETURN
!    END SUBROUTINE DO_GC_CHEM

END MODULE GIGC_ChemDr
#endif
