#include "MAPL_Generic.h"
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: gigc_chunk_mod
!
! !DESCRIPTION: Module GC\_CHUNK\_MOD is the module that contains the init,
!  and run methods for the ESMF interface to GEOS-Chem.
!\\
!\\
! !INTERFACE: 
!      
MODULE GIGC_Chunk_Mod
!
! !USES:
!      
  USE MAPL_MOD
  USE ESMF
  USE ErrCode_Mod
  USE Precision_Mod

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: GIGC_Chunk_Init
  PUBLIC :: GIGC_Chunk_Run
!
! !PRIVATE MEMBER FUNCTIONS:
!
  INTEGER  ::  MemDebugLevel
!
! !REVISION HISTORY:
!  22 Jun 2009 - R. Yantosca & P. Le Sager - Chunkized & cleaned up.
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC

CONTAINS
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: gigc_chunk_init
!
! !DESCRIPTION: Subroutine GIGC\_CHUNK\_INIT is the ESMF init method for
!  GEOS-Chem.  This routine calls routines within core GEOS-Chem to allocate 
!  arrays and read input files.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE GIGC_Chunk_Init( nymdB,         nhmsB,      nymdE,           &
                              nhmsE,         tsChem,     tsDyn,           &
                              lonCtr,        latCtr,                      &
#if !defined( MODEL_GEOS )
                              GC,            EXPORT,                      &
#endif
                              Input_Opt,     State_Chm,  State_Diag,      &
                              State_Grid,    State_Met,  HcoConfig,       &
                              HistoryConfig, RC )
!
! !USES:
!
    USE Chemistry_Mod,           ONLY : Init_Chemistry
    USE Emissions_Mod,           ONLY : Emissions_Init
    USE GC_Environment_Mod
    USE GC_Grid_Mod,             ONLY : SetGridFromCtr
    USE GIGC_HistoryExports_Mod, ONLY : HistoryConfigObj
    USE HCO_Types_Mod,           ONLY : ConfigObj
    USE Input_Mod,               ONLY : Read_Input_File
    USE Input_Opt_Mod,           ONLY : OptInput, Set_Input_Opt
    USE Linoz_Mod,               ONLY : Linoz_Read
    USE PBL_Mix_Mod,             ONLY : Init_PBL_Mix
    USE PhysConstants,           ONLY : PI_180
    USE Pressure_Mod,            ONLY : Init_Pressure
    USE Roundoff_Mod,            ONLY : RoundOff
    USE State_Chm_Mod,           ONLY : ChmState
    USE State_Diag_Mod,          ONLY : DgnState
    USE State_Grid_Mod,          ONLY : GrdState, Init_State_Grid
    USE State_Met_Mod,           ONLY : MetState
    USE Strat_Chem_Mod,          ONLY : Init_Strat_Chem
#if defined( MODEL_GEOS )
    USE Tendencies_Mod,          ONLY : TEND_INIT
#endif
    USE Time_Mod,                ONLY : Set_Timesteps
    USE UCX_MOD,                 ONLY : INIT_UCX
    USE UnitConv_Mod,            ONLY : Convert_Spc_Units
!
! !INPUT PARAMETERS:
!
    INTEGER,            INTENT(IN)    :: nymdB       ! YYYYMMDD @ start of run
    INTEGER,            INTENT(IN)    :: nhmsB       ! hhmmss   @ start of run
    INTEGER,            INTENT(IN)    :: nymdE       ! YYYYMMDD @ end of run
    INTEGER,            INTENT(IN)    :: nhmsE       ! hhmmss   @ end of run
    REAL,               INTENT(IN)    :: tsChem      ! Chemistry timestep [s]
    REAL,               INTENT(IN)    :: tsDyn       ! Chemistry timestep [s]
    REAL(ESMF_KIND_R4), INTENT(IN)    :: lonCtr(:,:) ! Lon centers [radians]
    REAL(ESMF_KIND_R4), INTENT(IN)    :: latCtr(:,:) ! Lat centers [radians]
!
! !INPUT/OUTPUT PARAMETERS:
!
#if !defined( MODEL_GEOS )
    TYPE(ESMF_State),    INTENT(INOUT), TARGET :: EXPORT ! Export state object
    TYPE(ESMF_GridComp), INTENT(INOUT)         :: GC     ! Ref to this GridComp
#endif
    TYPE(OptInput),      INTENT(INOUT) :: Input_Opt      ! Input Options object
    TYPE(ChmState),      INTENT(INOUT) :: State_Chm      ! Chem State object 
    TYPE(DgnState),      INTENT(INOUT) :: State_Diag     ! Diag State object
    TYPE(GrdState),      INTENT(INOUT) :: State_Grid     ! Grid State object
    TYPE(MetState),      INTENT(INOUT) :: State_Met      ! Met State object
    TYPE(ConfigObj),     POINTER       :: HcoConfig      ! HEMCO config obj 
    TYPE(HistoryConfigObj), POINTER    :: HistoryConfig  ! History config obj 
!
! !OUTPUT PARAMETERS:
!
    INTEGER,             INTENT(OUT)   :: RC             ! Success or failure?
!
! !REMARKS:
!  Need to add better error checking
!
! !REVISION HISTORY: 
!  18 Jul 2011 - M. Long     - Initial Version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER                        :: I, J, L, STATUS
    CHARACTER(LEN=ESMF_MAXSTR)     :: Iam
    TYPE(ESMF_Config)              :: CF            ! Grid comp config object

    !=======================================================================
    ! GIGC_CHUNK_INIT begins here 
    !=======================================================================

    ! Error trap
    Iam = 'GIGC_CHUNK_INIT (gigc_chunk_mod.F90)'

    ! Assume success
    RC = GC_SUCCESS

#if !defined( MODEL_GEOS )
    ! Get memory debug level
    call ESMF_GridCompGet ( GC, config=CF, RC=STATUS )
    _VERIFY(STATUS)
    call ESMF_ConfigGetAttribute(CF, MemDebugLevel, &
                                 Label="MEMORY_DEBUG_LEVEL:" , RC=STATUS)
    _VERIFY(STATUS)
#endif

    ! Read input.geos at very beginning of simulation on every thread
    CALL Read_Input_File( Input_Opt, State_Grid, RC )
    _ASSERT(RC==GC_SUCCESS, 'Error calling Read_Input_File')

    ! Initialize GEOS-Chem horizontal grid structure
    CALL GC_Init_Grid( Input_Opt, State_Grid, RC )
    _ASSERT(RC==GC_SUCCESS, 'Error calling GC_Init_Grid')

    ! Set maximum number of levels in the chemistry grid
    IF ( Input_Opt%LUCX ) THEN
       State_Grid%MaxChemLev  = State_Grid%MaxStratLev
    ELSE
       State_Grid%MaxChemLev  = State_Grid%MaxTropLev
    ENDIF

    ! In the ESMF/MPI environment, we can get the total overhead ozone
    ! either from the met fields (GIGCsa) or from the Import State (GEOS-5)
    Input_Opt%USE_O3_FROM_MET = .TRUE.
    
    ! Read LINOZ climatology
    IF ( Input_Opt%LLINOZ ) THEN
       CALL Linoz_Read( Input_Opt, RC ) 
       _ASSERT(RC==GC_SUCCESS, 'Error calling Linoz_Read')
    ENDIF

    ! Allocate all lat/lon arrays
    CALL GC_Allocate_All( Input_Opt, State_Grid, RC )
    _ASSERT(RC==GC_SUCCESS, 'Error calling GC_Allocate_All')

    ! Set grid based on passed mid-points
    CALL SetGridFromCtr( Input_Opt, State_Grid, lonCtr, latCtr, RC )
    _ASSERT(RC==GC_SUCCESS, 'Error caling')

    ! Update Input_Opt with timing fields
    Input_Opt%NYMDb   = nymdB
    Input_Opt%NHMSb   = nhmsB
    Input_Opt%NYMDe   = nymdE
    Input_Opt%NHMSe   = nhmsE
    Input_Opt%TS_CHEM = INT( tsChem )   ! Chemistry timestep [sec]
    Input_Opt%TS_EMIS = INT( tsChem )   ! Chemistry timestep [sec]
    Input_Opt%TS_DYN  = INT( tsDyn  )   ! Dynamic   timestep [sec]
    Input_Opt%TS_CONV = INT( tsDyn  )   ! Dynamic   timestep [sec]

    ! Set GEOS-Chem timesteps on all CPUs
    CALL SET_TIMESTEPS( Input_Opt,                                       &
                        Chemistry  = Input_Opt%TS_CHEM,                  &
                        Convection = Input_Opt%TS_CONV,                  &
                        Dynamics   = Input_Opt%TS_DYN,                   &
                        Emission   = Input_Opt%TS_EMIS,                  &
                        Radiation  = Input_Opt%TS_RAD,                   &
                        Unit_Conv  = MAX( Input_Opt%TS_DYN,              &
                                          Input_Opt%TS_CONV ),           &
                        Diagnos    = Input_Opt%TS_DIAG         )

    ! Initialize derived-type objects for met, chem, and diag
    CALL GC_Init_StateObj( HistoryConfig%DiagList, Input_Opt, &
                           State_Chm, State_Diag, State_Grid, State_Met, RC )
    _ASSERT(RC==GC_SUCCESS, 'Error calling GC_Init_StateObj')

    ! Initialize other GEOS-Chem modules
    CALL GC_Init_Extra( HistoryConfig%DiagList, Input_Opt,    &
                        State_Chm, State_Diag, State_Grid, RC ) 
    _ASSERT(RC==GC_SUCCESS, 'Error calling GC_Init_Extra')

    ! Set initial State_Chm%Species units to units expected in transport
# if defined( MODEL_GEOS )
    State_Chm%Spc_Units = 'kg/kg total'
#else
    State_Chm%Spc_Units = 'kg/kg dry'
#endif

    ! Initialize the GEOS-Chem pressure module (set Ap & Bp)
    CALL Init_Pressure( Input_Opt, State_Grid, RC )

    ! Initialize the PBL mixing module
    CALL INIT_PBL_MIX( Input_Opt, State_Grid, RC )
    _ASSERT(RC==GC_SUCCESS, 'Error calling INIT_PBL_MIX')

    ! Initialize chemistry mechanism
    IF ( Input_Opt%ITS_A_FULLCHEM_SIM .OR. Input_Opt%ITS_AN_AEROSOL_SIM ) THEN
       CALL INIT_CHEMISTRY ( Input_Opt,  State_Chm, State_Diag, &
                             State_Grid, RC )
       _ASSERT(RC==GC_SUCCESS, 'Error calling INIT_CHEMISTRY')
    ENDIF

    ! Initialize HEMCO
    CALL EMISSIONS_INIT( Input_Opt, State_Chm, State_Grid, State_Met, RC, &
                         HcoConfig=HcoConfig )
    _ASSERT(RC==GC_SUCCESS, 'Error calling EMISSIONS_INIT')

    ! Stratosphere - can't be initialized without HEMCO because of STATE_PSC
    IF ( Input_Opt%LUCX ) THEN

       ! Initialize stratospheric routines
       CALL INIT_UCX( Input_Opt, State_Chm, State_Diag, State_Grid )

    ENDIF

#if defined( MODEL_GEOS )
    ! Keep commented out line with LLSTRAT as a GEOS-5 option reminder
    !IF ( Input_Opt%LSCHEM .AND. Input_Opt%LLSTRAT < value_LM ) THEN
#endif
     IF ( Input_Opt%LSCHEM ) THEN
       CALL INIT_STRAT_CHEM( Input_Opt, State_Chm, State_Met, State_Grid, RC )
       _ASSERT(RC==GC_SUCCESS, 'Error calling INIT_STRAT_CHEM')
    ENDIF

    !-------------------------------------------------------------------------
    ! Diagnostics and tendencies 
    !-------------------------------------------------------------------------

#if defined( MODEL_GEOS )
    ! The GEOS-Chem diagnostics list, stored in HistoryConfig, is initialized 
    ! during GIGC_INIT_SIMULATION, and corresponding arrays in State_Diag are 
    ! allocated accordingly when initializing State_Diag. Here, we thus 
    ! only need to initialize the tendencies, which have not been initialized
    ! yet (ckeller, 11/29/17). 
    CALL Tend_Init ( Input_Opt, State_Chm, State_Grid, State_Met, RC ) 
    _ASSERT(RC==GC_SUCCESS, 'Error calling Tend_Init')
#endif

#if !defined( MODEL_GEOS )
    ! GCHP only: Convert species units to internal state units (v/v dry)
    CALL Convert_Spc_Units( Input_Opt, State_Chm, State_Grid, State_Met, &
                            'v/v dry', RC )
    _ASSERT(RC==GC_SUCCESS, 'Error calling Convert_Spc_Units')
#endif

    ! Return success
    RC = GC_Success

  END SUBROUTINE GIGC_Chunk_Init
!EOC

!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: gigc_chunk_run
!
! !DESCRIPTION: Subroutine GIGC\_CHUNK\_RUN is the ESMF run method for
!  GEOS-Chem.
!
! !INTERFACE:
!
  SUBROUTINE GIGC_Chunk_Run( GC,                                             &
                             nymd,       nhms,       year,       month,      &
                             day,        dayOfYr,    hour,       minute,     &
                             second,     utc,        hElapsed,   Input_Opt,  &
                             State_Chm,  State_Diag, State_Grid, State_Met,  &
                             Phase,      IsChemTime,                         &
#if defined( MODEL_GEOS )
                             FrstRewind, &
#endif
                             RC )
!
! !USES:
!
    ! GEOS-Chem state objects 
    USE HCO_Interface_Mod,  ONLY : HcoState
    USE Input_Opt_Mod,      ONLY : OptInput
    USE State_Chm_Mod,      ONLY : ChmState
    USE State_Diag_Mod
    USE State_Grid_Mod,     ONLY : GrdState
    USE State_Met_Mod,      ONLY : MetState

    ! GEOS-Chem components
    USE Chemistry_Mod,      ONLY : Do_Chemistry, Recompute_OD
    USE Convection_Mod,     ONLY : Do_Convection
    USE DryDep_Mod,         ONLY : Do_DryDep
    USE Emissions_Mod,      ONLY : Emissions_Run
    USE Mixing_Mod,         ONLY : Do_Tend, Do_Mixing
    USE WetScav_Mod,        ONLY : Setup_WetScav, Do_WetDep

    ! Specialized subroutines
    USE Calc_Met_Mod,       ONLY : AirQnt, Set_Dry_Surface_Pressure
    USE Calc_Met_Mod,       ONLY : GIGC_Cap_Tropopause_Prs
    USE Set_Global_CH4_Mod, ONLY : Set_CH4
    USE MODIS_LAI_Mod,      ONLY : Compute_XLAI
    USE PBL_Mix_Mod,        ONLY : Compute_PBL_Height
    USE Pressure_Mod,       ONLY : Set_Floating_Pressures
    USE TOMS_Mod,           ONLY : Compute_Overhead_O3
    USE UCX_Mod,            ONLY : Set_H2O_Trac

    ! Utilities
    USE ErrCode_Mod
    USE HCO_Error_Mod
    USE HCO_Interface_Mod,  ONLY : SetHcoTime
    USE MAPL_MemUtilsMod
    USE Pressure_Mod,       ONLY : Accept_External_Pedge
    USE State_Chm_Mod,      ONLY : IND_
    USE Time_Mod,           ONLY : Accept_External_Date_Time
    USE UnitConv_Mod,       ONLY : Convert_Spc_Units

    ! Diagnostics
    USE Diagnostics_Mod,    ONLY : Set_Diagnostics_EndofTimestep
    USE Aerosol_Mod,        ONLY : Set_AerMass_Diagnostic

#if defined( MODEL_GEOS )
    USE Derived_Met_Mod,    ONLY : GET_COSINE_SZA
    USE HCOI_GC_MAIN_MOD,   ONLY : HCOI_GC_WriteDiagn
#endif
!
! !INPUT PARAMETERS:
!
    INTEGER,        INTENT(IN)    :: nymd        ! YYYY/MM/DD @ current time
    INTEGER,        INTENT(IN)    :: nhms        ! hh:mm:ss   @ current time
    INTEGER,        INTENT(IN)    :: year        ! UTC year 
    INTEGER,        INTENT(IN)    :: month       ! UTC month
    INTEGER,        INTENT(IN)    :: day         ! UTC day
    INTEGER,        INTENT(IN)    :: dayOfYr     ! UTC day of year
    INTEGER,        INTENT(IN)    :: hour        ! UTC hour
    INTEGER,        INTENT(IN)    :: minute      ! UTC minute
    INTEGER,        INTENT(IN)    :: second      ! UTC second
    REAL*4,         INTENT(IN)    :: utc         ! UTC time [hrs]
    REAL*4,         INTENT(IN)    :: hElapsed    ! Elapsed hours
    INTEGER,        INTENT(IN)    :: Phase       ! Run phase (-1, 1 or 2)
    LOGICAL,        INTENT(IN)    :: IsChemTime  ! Time for chemistry? 
#if defined( MODEL_GEOS )
    LOGICAL,        INTENT(IN)    :: FrstRewind  ! Is it the first rewind? 
#endif
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ESMF_GridComp), INTENT(INOUT) :: GC          ! Ref to this GridComp
    TYPE(OptInput),      INTENT(INOUT) :: Input_Opt   ! Input Options obj
    TYPE(ChmState),      INTENT(INOUT) :: State_Chm   ! Chemistry State obj
    TYPE(DgnState),      INTENT(INOUT) :: State_Diag  ! Diagnostics State obj
    TYPE(GrdState),      INTENT(INOUT) :: State_Grid  ! Grid State obj
    TYPE(MetState),      INTENT(INOUT) :: State_Met   ! Meteorology State obj
!
! !OUTPUT PARAMETERS:
!
    INTEGER,             INTENT(OUT)   :: RC          ! Return code
!
! !REMARKS:
!
! !REVISION HISTORY:
!  18 Jul 2011 - M. Long     - Initial Version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
    TYPE(ESMF_STATE)               :: INTSTATE
    TYPE(MAPL_MetaComp), POINTER   :: STATE
    TYPE(ESMF_VM)                  :: VM            ! ESMF VM object
    TYPE(ESMF_Field)               :: IntField
    REAL*8                         :: DT
    CHARACTER(LEN=ESMF_MAXSTR)     :: Iam, OrigUnit
    INTEGER                        :: STATUS, HCO_PHASE, RST
#if defined( MODEL_GEOS )
    INTEGER                        :: N, I, J, L
#endif

    ! Local logicals to turn on/off individual components
    ! The parts to be executed are based on the input options,
    ! the time step and the phase.
    LOGICAL                        :: DoConv 
    LOGICAL                        :: DoDryDep
    LOGICAL                        :: DoEmis
    LOGICAL                        :: DoTend 
    LOGICAL                        :: DoTurb 
    LOGICAL                        :: DoChem
    LOGICAL                        :: DoWetDep

    ! First call?
    LOGICAL, SAVE                  :: FIRST = .TRUE.

    ! # of times this routine has been called. Only temporary for printing 
    ! processes on the first 10 calls.
    INTEGER, SAVE                  :: NCALLS = 0

    ! Strat. H2O settings 
    LOGICAL                        :: SetStratH2O 
#if defined( MODEL_GEOS )
    LOGICAL, SAVE                  :: LSETH2O_orig
#endif

    ! Whether to scale mixing ratio with meteorology update in AirQnt
    LOGICAL, SAVE                  :: scaleMR = .FALSE.

    !=======================================================================
    ! GIGC_CHUNK_RUN begins here 
    !=======================================================================

    ! Error trap
    Iam = 'GIGC_CHUNK_RUN (gigc_chunk_mod.F90)'

    ! Assume success
    RC = GC_SUCCESS

    ! Get state object (needed for timers)
    CALL MAPL_GetObjectFromGC(GC, STATE, __RC__)

    ! Get the VM for optional memory prints (level >= 2)
    !-----------------------------------
    if ( MemDebugLevel > 0 ) THEN
       call ESMF_VmGetCurrent(VM, RC=STATUS)
       _VERIFY(STATUS)
    endif

    !=======================================================================
    ! Define processes to be covered in this phase
    !
    ! In the standard GEOS-Chem, the following operator sequence is used:
    ! 1. DryDep (kg)
    ! 2. Emissions (kg)
    ! 3. Turbulence (v/v)
    ! 4. Convection (v/v)
    ! 5. Chemistry (kg)
    ! 6. Wetdep (kg)
    !
    ! The GEOS-5 operator sequence is:
    ! 1. Gravity wave drag
    ! 2. Moist (convection)
    ! 3. Chemistry 1 (drydep and emissions)
    ! 4. Surface 1
    ! 5. Turbulence 1
    ! 6. Surface 2
    ! 7. Turbulence 2
    ! 8. Chemistry 2 (chemistry and wet deposition)
    ! 9. Radiation 
    !
    ! Here, we use the following operator sequence:
    ! 
    ! 1.  Convection (v/v) --> Phase 1
    ! 2.  DryDep (kg)      --> Phase 1
    ! 3.  Emissions (kg)   --> Phase 1
    ! 4a. Tendencies (v/v) --> Phase 1
    ! -------------------------------
    ! 4b. Turbulence (v/v) --> Phase 2 
    ! 5.  Chemistry (kg)   --> Phase 2
    ! 6.  WetDep (kg)      --> Phase 2     
    ! 
    ! Any of the listed processes is only executed if the corresponding switch
    ! in the input.geos file is enabled. If the physics component already
    ! covers convection or turbulence, they should not be applied here!
    ! The tendencies are only applied if turbulence is not done within
    ! GEOS-Chem (ckeller, 10/14/14).
    ! 
    ! The standard number of phases in GCHP is 1, set in GCHP.rc, which
    ! results in Phase -1 in gigc_chunk_run. This results in executing
    ! all GEOS-Chem components in a single run rather than splitting up
    ! across two runs as is done in GEOS-5. (ewl, 10/26/18)
    !=======================================================================

    ! By default, do processes as defined in input.geos. DoTend defined below. 
    DoConv   = Input_Opt%LCONV                    ! dynamic time step
    DoDryDep = Input_Opt%LDRYD .AND. IsChemTime   ! chemistry time step
    DoEmis   = IsChemTime                         ! chemistry time step
#if defined( MODEL_GEOS )
    DoTurb   = Input_Opt%LTURB .AND. IsChemTime   ! dynamic time step
#else
    DoTurb   = Input_Opt%LTURB                    ! dynamic time step
#endif
    DoChem   = Input_Opt%LCHEM .AND. IsChemTime   ! chemistry time step
    DoWetDep = Input_Opt%LWETD                    ! dynamic time step 

#if defined( MODEL_GCHP )
    ! Make sure chemistry timestep components are run in first timestep
    ! even if IsChemTime is false. This is a temporary kludge pending a
    ! better solution for the GIGCchem alarm (ewl, 2/21/20)
    IF ( FIRST ) THEN
       DoDryDep = Input_Opt%LDRYD
       DoEmis   = .TRUE.
       DoChem   = Input_Opt%LCHEM
    ENDIF
#endif

    ! If Phase is not -1, only do selected processes for given phases: 
    ! Phase 1: disable turbulence, chemistry and wet deposition.
    IF ( Phase == 1 ) THEN
       DoTurb   = .FALSE.
       DoChem   = .FALSE.
       DoWetDep = .FALSE.

    ! Phase 2: disable convection, drydep and emissions. 
    ELSEIF ( Phase == 2 ) THEN
       DoConv   = .FALSE.
       DoDryDep = .FALSE.
       DoEmis   = .FALSE. 
    ENDIF

    ! Check if tendencies need be applied. The drydep and emission calls
    ! only calculates the emission / drydep rates, but do not apply the
    ! tendencies to the tracer array yet. If turbulence is done as part of
    ! GEOS-5, we need to make sure that these tendencies are applied to the
    ! tracer array. If turbulence is explicitly covered by GEOS-Chem,
    ! however, the tendencies become automatically applied within the PBL
    ! mixing routines (DO_MIXING), so we should never apply the tendencies
    ! in this case.
    DoTend = ( DoEmis .OR. DoDryDep ) .AND. .NOT. Input_Opt%LTURB

    ! testing only
    IF ( Input_Opt%AmIRoot .and. NCALLS < 10 ) THEN 
       write(*,*) 'GEOS-Chem phase ', Phase, ':'
       write(*,*) 'DoConv   : ', DoConv
       write(*,*) 'DoDryDep : ', DoDryDep
       write(*,*) 'DoEmis   : ', DoEmis
       write(*,*) 'DoTend   : ', DoTend
       write(*,*) 'DoTurb   : ', DoTurb
       write(*,*) 'DoChem   : ', DoChem
       write(*,*) 'DoWetDep : ', DoWetDep
       write(*,*) ' '
    ENDIF

    !-------------------------------------------------------------------------
    ! Pre-Run assignments
    !-------------------------------------------------------------------------

    ! Eventually initialize/reset wetdep
    IF ( DoConv .OR. DoChem .OR. DoWetDep ) THEN
       CALL SETUP_WETSCAV( Input_Opt, State_Chm, State_Grid, State_Met, RC )
    ENDIF

    ! Pass time values obtained from the ESMF environment to GEOS-Chem
    CALL Accept_External_Date_Time( value_NYMD     = nymd,       &  
                                    value_NHMS     = nhms,       &  
                                    value_YEAR     = year,       &  
                                    value_MONTH    = month,      &  
                                    value_DAY      = day,        &  
                                    value_DAYOFYR  = dayOfYr,    &  
                                    value_HOUR     = hour,       &  
                                    value_MINUTE   = minute,     &  
                                    value_HELAPSED = hElapsed,   & 
                                    value_UTC      = utc,        &
                                    RC             = RC         )

    ! Set HEMCO time
    CALL SetHcoTime ( DoEmis, RC )

    ! Calculate MODIS leaf area indexes needed for dry deposition
    CALL Compute_XLAI( Input_Opt, State_Grid, State_Met, RC )

    ! Set the pressure at level edges [hPa] from the ESMF environment
    CALL Accept_External_Pedge( State_Met = State_Met,  &
                                RC        = RC         )

    ! Set dry surface pressure (PS1_DRY) from State_Met%PS1_WET
    CALL SET_DRY_SURFACE_PRESSURE( State_Grid, State_Met, 1 )

    ! Set dry surface pressure (PS2_DRY) from State_Met%PS2_WET
    CALL SET_DRY_SURFACE_PRESSURE( State_Grid, State_Met, 2 )

    ! Initialize surface pressures to match the post-advection pressures
    State_Met%PSC2_WET = State_Met%PS1_WET
    State_Met%PSC2_DRY = State_Met%PS1_DRY
    CALL SET_FLOATING_PRESSURES( State_Grid, State_Met, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    ! Define airmass and related quantities
#if defined( MODEL_GEOS )
    CALL AirQnt( Input_Opt, State_Chm, State_Grid, State_Met, RC, .FALSE. )
#else
    ! Scale mixing ratio with changing met only if FV advection is off.
    ! Only do this the first timestep if DELP_DRY found in restart.
    IF ( FIRST .and. .not. Input_Opt%LTRAN ) THEN
       CALL MAPL_Get ( STATE, INTERNAL_ESMF_STATE=INTSTATE, __RC__ )
       CALL ESMF_StateGet( INTSTATE, 'DELP_DRY', IntField, RC=STATUS )
       _VERIFY(STATUS)
       CALL ESMF_AttributeGet( IntField, NAME="RESTART", VALUE=RST, RC=STATUS )
       _VERIFY(STATUS)
       IF ( .not. ( RST == MAPL_RestartBootstrap .OR. &
                    RST == MAPL_RestartSkipInitial ) ) scaleMR = .TRUE.
       CALL AirQnt( Input_Opt, State_Chm, State_Grid, State_Met, RC, scaleMR )
       scaleMR = .TRUE.
    ELSE
       CALL AirQnt( Input_Opt, State_Chm, State_Grid, State_Met, RC, scaleMR )
    ENDIF
#endif

    ! Cap the polar tropopause pressures at 200 hPa, in order to avoid
    ! tropospheric chemistry from happening too high up (cf. J. Logan)
    CALL GIGC_Cap_Tropopause_Prs( Input_Opt      = Input_Opt,  &
                                  State_Grid     = State_Grid, &
                                  State_Met      = State_Met,  &
                                  RC             = RC         )

    ! Call PBL quantities. Those are always needed
    CALL COMPUTE_PBL_HEIGHT( State_Grid, State_Met, RC )
    
    ! Convert to dry mixing ratio
    CALL Convert_Spc_Units ( Input_Opt, State_Chm, State_Grid, State_Met, &
                             'kg/kg dry', RC, OrigUnit=OrigUnit )

    ! SDE 05/28/13: Set H2O to STT if relevant
    IF ( IND_('H2O','A') > 0 ) THEN
       SetStratH2O = .FALSE.
#if defined( MODEL_GEOS )
       !=======================================================================
       ! Tropospheric H2O is always prescribed (using GEOS Q). For strat H2O
       ! there are three options, controlled by toggles 'set initial global MR'
       ! in input.geos and 'Prescribe_strat_H2O' in GEOSCHEMchem_GridComp.rc: 
       ! (A) never prescribe strat H2O -> both toggles off
       ! (B) prescribe strat H2O on init time step -> toggle in input.goes on
       ! (C) always prescribe strat H2O -> toggle in GEOSCHEMchem_GridComp.rc on
       !=======================================================================
       IF ( FIRST ) THEN
          LSETH2O_orig = Input_Opt%LSETH2O
       ENDIF
       IF ( FIRST .OR. FrstRewind ) THEN
          Input_Opt%LSETH2O = LSETH2O_orig
       ENDIF
       IF ( Input_Opt%LSETH2O .OR. .NOT. Input_Opt%LUCX .OR. &
            Input_Opt%AlwaysSetH2O ) THEN
#else
       IF ( Input_Opt%LSETH2O .OR. .NOT. Input_Opt%LUCX ) THEN
#endif
          SetStratH2O = .TRUE.
       ENDIF
       CALL SET_H2O_TRAC( SetStratH2O, Input_Opt, State_Chm, & 
                          State_Grid,  State_Met, RC )
       _ASSERT(RC==GC_SUCCESS, 'Error calling SET_H2O_TRAC')

      ! Only force strat once if using UCX
       IF (Input_Opt%LSETH2O) Input_Opt%LSETH2O = .FALSE.
    ENDIF

#if defined( MODEL_GEOS )
    ! GEOS-5 only: needed in GCHP?
    ! Compute the cosine of the solar zenith angle array:
    !    State_Met%SUNCOS     => COS(SZA) at the current time
    !    State_Met%SUNCOSmid  => COS(SZA) at the midpt of the chem timestep
    !    COS(SZA) at the midpt of the chem timestep 5hrs ago is now
    !    calculated elsewhere, in the HEMCO PARANOx extension
    CALL GET_COSINE_SZA( Input_Opt, State_Grid, State_Met, RC )
    _ASSERT(RC==GC_SUCCESS, 'Error calling GET_COSINE_SZA')
#endif

    !=======================================================================
    ! EMISSIONS. Pass HEMCO Phase 1 which only updates the HEMCO clock
    ! and the HEMCO data list. Should be called every time to make sure 
    ! that the HEMCO clock and the HEMCO data list are up to date.
    !=======================================================================
    HCO_PHASE = 1
    CALL EMISSIONS_RUN( Input_Opt, State_Chm, State_Diag, &
                        State_Grid, State_Met, DoEmis, HCO_PHASE, RC  )
    _ASSERT(RC==GC_SUCCESS, 'Error calling EMISSIONS_RUN')

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!                                PHASE 1 or -1                           !!!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !=======================================================================
    ! 1. Convection
    ! 
    ! Call GEOS-Chem internal convection routines if convection is enabled
    ! in input.geos. This should only be done if convection is not covered
    ! by another gridded component and/or the GC species are not made
    ! friendly to this component!!
    !=======================================================================
    IF ( DoConv ) THEN
       if(Input_Opt%AmIRoot.and.NCALLS<10) write(*,*) ' --- Do convection now'
       CALL MAPL_TimerOn( STATE, 'GC_CONV' )

       CALL DO_CONVECTION ( Input_Opt, State_Chm, State_Diag, &
                            State_Grid, State_Met, RC )
       _ASSERT(RC==GC_SUCCESS, 'Error calling DO_CONVECTION')
 
       CALL MAPL_TimerOff( STATE, 'GC_CONV' )
       if(Input_Opt%AmIRoot.and.NCALLS<10) write(*,*) ' --- Convection done!'
    ENDIF   

    !=======================================================================
    ! 2. Dry deposition
    !
    ! Calculates the deposition rates in [s-1].
    !=======================================================================
    IF ( DoDryDep ) THEN
       if(Input_Opt%AmIRoot.and.NCALLS<10) THEN
          write(*,*) ' --- Do drydep now'
          write(*,*) '     Use FULL PBL: ', Input_Opt%PBL_DRYDEP
       endif
       CALL MAPL_TimerOn( STATE, 'GC_DRYDEP' )
    
       ! Do dry deposition
       CALL Do_DryDep ( Input_Opt, State_Chm, State_Diag, &
                        State_Grid, State_Met, RC ) 
       _ASSERT(RC==GC_SUCCESS, 'Error calling Do_DryDep')

       CALL MAPL_TimerOff( STATE, 'GC_DRYDEP' )
       if(Input_Opt%AmIRoot.and.NCALLS<10) write(*,*) ' --- Drydep done!'
    ENDIF

    !=======================================================================
    ! 3. Emissions (HEMCO)
    !
    ! HEMCO must be called on first time step to make sure that the HEMCO
    ! data lists are all properly set up. 
    !=======================================================================
    IF ( DoEmis ) THEN
#if !defined( MODEL_GEOS )
       ! Optional memory prints (level >= 3)
       if ( MemDebugLevel > 0 ) THEN
          call ESMF_VMBarrier(VM, RC=STATUS)
          _VERIFY(STATUS)
          call MAPL_MemUtilsWrite(VM, &
                  'gigc_chunk_run, before Emissions_Run', RC=STATUS )
          _VERIFY(STATUS)
       endif
#endif

       if(Input_Opt%AmIRoot.and.NCALLS<10) write(*,*) ' --- Do emissions now'
       CALL MAPL_TimerOn( STATE, 'GC_EMIS' )

       ! Do emissions. Pass HEMCO Phase 2 which performs the emissions 
       ! calculations.
       HCO_PHASE = 2
       CALL EMISSIONS_RUN( Input_Opt, State_Chm, State_Diag, &
                           State_Grid, State_Met, DoEmis, HCO_PHASE, RC )
       _ASSERT(RC==GC_SUCCESS, 'Error calling EMISSIONS_RUN')

       CALL MAPL_TimerOff( STATE, 'GC_EMIS' )
       if(Input_Opt%AmIRoot.and.NCALLS<10) write(*,*) ' --- Emissions done!'

       ! Optional memory prints (level >= 3)
       if ( MemDebugLevel > 0 ) THEN
          call ESMF_VMBarrier(VM, RC=STATUS)
          _VERIFY(STATUS)
          call MAPL_MemUtilsWrite(VM,&
                  'gigc_chunk_run, after  Emissions_Run', RC=STATUS )
          _VERIFY(STATUS)
       endif

    ENDIF

    !=======================================================================
    ! If physics covers turbulence, simply add the emission and dry 
    ! deposition fluxes calculated above to the tracer array, without caring
    ! about the vertical distribution. The tracer tendencies are only added
    ! to the tracers array after emissions, drydep. So we need to use the
    ! emissions time step here.
    !=======================================================================
    IF ( DoTend ) THEN 
       if(Input_Opt%AmIRoot.and.NCALLS<10) write(*,*)   &
                           ' --- Add emissions and drydep to tracers'
       CALL MAPL_TimerOn( STATE, 'GC_FLUXES' )

       ! Get emission time step [s]. 
       _ASSERT(ASSOCIATED(HcoState), 'Error: HcoState not associated')
       DT = HcoState%TS_EMIS 

       ! Apply tendencies over entire PBL. Use emission time step.
       CALL DO_TEND( Input_Opt, State_Chm, State_Diag, &
                     State_Grid, State_Met, .FALSE., RC, DT=DT )
       _ASSERT(RC==GC_SUCCESS, 'Error calling DO_TEND')

       ! testing only
       if(Input_Opt%AmIRoot.and.NCALLS<10) write(*,*)   &
                                 '     Tendency time step [s]: ', DT 

       CALL MAPL_TimerOff( STATE, 'GC_FLUXES' )
       if(Input_Opt%AmIRoot.and.NCALLS<10) write(*,*)   &
                                 ' --- Fluxes applied to tracers!' 
    ENDIF ! Tendencies 

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!                              PHASE 2 or -1                             !!!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !=======================================================================
    ! 4. Turbulence
    !
    ! Call GEOS-Chem internal turbulence routines if turbulence is enabled
    ! in input.geos. This should only be done if turbulence is not covered
    ! by another gridded component and/or the GC species are not made
    ! friendly to this component!!
    !=======================================================================
    IF ( DoTurb ) THEN
       if(Input_Opt%AmIRoot.and.NCALLS<10) write(*,*) ' --- Do turbulence now'
       CALL MAPL_TimerOn( STATE, 'GC_TURB' )

       ! Do mixing and apply tendencies. This will use the dynamic time step,
       ! which is fine since this call will be executed on every time step. 
       CALL DO_MIXING ( Input_Opt, State_Chm, State_Diag, &
                        State_Grid, State_Met, RC )
       _ASSERT(RC==GC_SUCCESS, 'Error calling DO_MIXING')

       CALL MAPL_TimerOff( STATE, 'GC_TURB' )
       if(Input_Opt%AmIRoot.and.NCALLS<10) write(*,*) ' --- Turbulence done!'
    ENDIF

    ! Set tropospheric CH4 concentrations and fill species array with
    ! current values. 
#if defined( MODEL_GEOS )
    IF ( .NOT. Input_Opt%LCH4EMIS .AND. ( DoTurb .OR. DoTend ) ) THEN
#else
    IF ( Phase /= 2 .AND. Input_Opt%ITS_A_FULLCHEM_SIM  &
         .AND. IND_('CH4','A') > 0 ) THEN
#endif
       CALL SET_CH4 ( Input_Opt, State_Chm, State_Diag, &
                      State_Grid, State_Met, RC )
       _ASSERT(RC==GC_SUCCESS, 'Error calling SET_CH4')
    ENDIF

    !=======================================================================
    ! 5. Chemistry
    !=======================================================================
    IF ( DoChem ) THEN
       if(Input_Opt%AmIRoot.and.NCALLS<10) write(*,*) ' --- Do chemistry now'
       CALL MAPL_TimerOn( STATE, 'GC_CHEM' )

       ! Calculate TOMS O3 overhead. For now, always use it from the
       ! Met field. State_Met%TO3 is imported from PCHEM (ckeller, 10/21/2014).
       CALL COMPUTE_OVERHEAD_O3( Input_Opt, State_Grid, DAY, .TRUE., &
                                 State_Met%TO3 )

#if !defined( MODEL_GEOS )
       ! Set H2O to species value if H2O is advected
       IF ( IND_('H2O','A') > 0 ) THEN
          CALL SET_H2O_TRAC( (.not. Input_Opt%LUCX), Input_Opt, &
                             State_Chm, State_Grid, State_Met, RC )
       ENDIF
#endif

       ! Optional memory prints (level >= 3)
       if ( MemDebugLevel > 0 ) THEN
          call ESMF_VMBarrier(VM, RC=STATUS)
          _VERIFY(STATUS)
          call MAPL_MemUtilsWrite(VM, &
                  'gigc_chunk_run:, before Do_Chemistry', RC=STATUS )
          _VERIFY(STATUS)
       endif

       ! Do chemistry
       CALL Do_Chemistry( Input_Opt, State_Chm, State_Diag, &
                          State_Grid, State_Met, RC ) 
       _ASSERT(RC==GC_SUCCESS, 'Error calling Do_Chemistr')

       CALL MAPL_TimerOff( STATE, 'GC_CHEM' )
       if(Input_Opt%AmIRoot.and.NCALLS<10) write(*,*) ' --- Chemistry done!'

       ! Optional memory prints (level >= 3)
       if ( MemDebugLevel > 0 ) THEN
          call ESMF_VMBarrier(VM, RC=STATUS)
          _VERIFY(STATUS)
          call MAPL_MemUtilsWrite(VM, &
                  'gigc_chunk_run, after  Do_Chemistry', RC=STATUS )
          _VERIFY(STATUS)
       endif

    ENDIF

    !=======================================================================
    ! 6. Wet deposition
    !=======================================================================
    IF ( DoWetDep ) THEN
       if(Input_Opt%AmIRoot.and.NCALLS<10) write(*,*) ' --- Do wetdep now'
       CALL MAPL_TimerOn( STATE, 'GC_WETDEP' )

       ! Do wet deposition
       CALL DO_WETDEP( Input_Opt, State_Chm, State_Diag, &
                       State_Grid, State_Met, RC )
       _ASSERT(RC==GC_SUCCESS, 'Error calling DO_WETDEP')

       CALL MAPL_TimerOff( STATE, 'GC_WETDEP' )
       if(Input_Opt%AmIRoot.and.NCALLS<10) write(*,*) ' --- Wetdep done!'
    ENDIF

    !=======================================================================
    ! Diagnostics 
    !=======================================================================

    !==============================================================
    !      ***** U P D A T E  O P T I C A L  D E P T H *****          
    !==============================================================
    ! Recalculate the optical depth at the wavelength(s) specified
    ! in the Radiation Menu. This must be done before the call to any
    ! diagnostic and only on a chemistry timestep.
    ! (skim, 02/05/11)
    IF ( DoChem ) THEN
       CALL RECOMPUTE_OD ( Input_Opt, State_Chm, State_Diag, &
                           State_Grid, State_Met, RC )
       _ASSERT(RC==GC_SUCCESS, 'Error calling RECOMPUTE_OD')
    ENDIF

    if(Input_Opt%AmIRoot.and.NCALLS<10) write(*,*) ' --- Do diagnostics now'
    CALL MAPL_TimerOn( STATE, 'GC_DIAGN' )

    ! Set certain diagnostics dependent on state at end of step. This
    ! includes species concentration and dry deposition flux.
    CALL Set_Diagnostics_EndofTimestep( Input_Opt,  State_Chm, State_Diag, &
                                        State_Grid, State_Met, RC )
    _ASSERT(RC==GC_SUCCESS, 'Error calling Set_Diagnostics_EndofTimestep')

    ! Archive aerosol mass and PM2.5 diagnostics
    IF ( State_Diag%Archive_AerMass ) THEN
       CALL Set_AerMass_Diagnostic( Input_Opt,  State_Chm, State_Diag, &
                                    State_Grid, State_Met, RC )
       _ASSERT(RC==GC_SUCCESS, 'Error calling Set_AerMass_Diagnostic')
    ENDIF

    CALL MAPL_TimerOff( STATE, 'GC_DIAGN' )
    if(Input_Opt%AmIRoot.and.NCALLS<10) write(*,*) ' --- Diagnostics done!'

    !=======================================================================
    ! Convert State_Chm%Species units
    !=======================================================================
    CALL Convert_Spc_Units ( Input_Opt, State_Chm, State_Grid, State_Met, &
                             OrigUnit, RC )

#if defined( MODEL_GEOS )
    ! Save specific humidity and dry air mass for total mixing ratio 
    ! adjustment in next timestep, if needed (ewl, 11/8/18)
    State_Met%SPHU_PREV = State_Met%SPHU
#endif

    !=======================================================================
    ! Clean up
    !=======================================================================

    ! testing only
    IF ( PHASE /= 1 .AND. NCALLS < 10 ) NCALLS = NCALLS + 1 

    ! First call is done
    FIRST = .FALSE.

    ! Return success
    RC = GC_SUCCESS

  END SUBROUTINE GIGC_Chunk_Run
!EOC
END MODULE GIGC_Chunk_Mod
