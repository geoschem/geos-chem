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
!  08 Nov 2012 - R. Yantosca - Deleted obsolete, commented-out stuff
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
  SUBROUTINE GIGC_Do_Chem( am_I_Root, NI,        NJ,         &
                           NL,        NCNST,     Input_Opt,  &
                           State_Chm, State_Met, RC         )
!
! !USES:
!
    USE CHEMISTRY_MOD,      ONLY : DO_CHEMISTRY
    USE CMN_DEP_MOD,        ONLY : FRCLND
    USE CMN_SIZE_MOD,       ONLY : DJSIZE, DISIZE, LLPAR, IIPAR, JJPAR
    USE COMODE_MOD,         ONLY : AIRDENS, CSPEC_FULL
    USE COMODE_LOOP_MOD
    USE DAO_MOD
    USE GIGC_Input_Opt_Mod, ONLY : OptInput
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
    LOGICAL,        INTENT(IN)    :: am_I_Root   ! Are we on the root CPU
    INTEGER,        INTENT(IN)    :: NI          ! # of longitudes
    INTEGER,        INTENT(IN)    :: NJ          ! # of latitudes
    INTEGER,        INTENT(IN)    :: NL          ! # of levels
    INTEGER,        INTENT(IN)    :: NCNST       ! # of constituents
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
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
!  08 Nov 2012 - R. Yantosca - Now pass the Input_Opt object to DO_CHEMISTRY
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
      IF ( am_I_Root ) THEN
         WRITE(6,*) '##### GIGC_Do_Chem, TRC_OX before chem [v/v]'
         WRITE(6,*) State_Chm%Tracers(1,1,:,2)
         WRITE(6,*) LCHEM
      ENDIF

      ! If we are doing chemistry
      IF ( LCHEM ) THEN

         ! The tracer concentrations in State_Chm%TRACERS state have units 
         ! of [v/v] (since these come from the ESMF internal state).  We need
         ! to convert these to [kg] before calling the GEOS-Chem chemistry.
         CALL Convert_Units( 2, N_TRACERS, TCVV, AD, State_Chm%Tracers )

         ! Call the GEOS-Chem chemistry routines
         CALL Do_Chemistry( am_I_Root, NI,        NJ,        NL,  &
                            Input_Opt, State_Chm, State_Met, RC  )

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
END MODULE GIGC_ChemDr
#endif
