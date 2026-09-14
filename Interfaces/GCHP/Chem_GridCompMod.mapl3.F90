#include "MAPL.h"

!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Model                            !
!------------------------------------------------------------------------------

! !MODULE: GEOSCHEMchem_GridCompMod

! !DESCRIPTION: GEOSCHEMchem_GridComp is an ESMF5 gridded component
! implementing the GEOS-Chem chemistry and related processes, including
! dry deposition, emissions, and wet deposition. In addition, the
! parameterizations for PBL mixing and convection as used in GEOS-Chem
! can be invoked by enabling the corresponding option in the GEOS-Chem
! input file (geoschem_config.yml). In this case, the corresponding GEOS-5
! process must NOT be applied to the GC tracers, i.e. the tracers must
! not be friendly to turbulence (if PBL mixing is used) and/or moist
! (for convection).

! This gridded component contains three run phases:
!
!  -1: Phase -1 is the standard setting in GCHP. It executes all components.
!      Phase is -1 if number of phases is set to 1 in config file GCHP.rc.
!
!   1: Phase 1 is used in GEOS-5. It executes convection, dry deposition,
!      and emissions and should be called before surface processes/turbulence.
!
!   2: Phase 2 is used in GEOS-5. It performs chemistry, and wet deposition,
!      and should be called after turbulence.
!
! GEOS-5 only:
! All GEOS-Chem species are stored in the GEOSCHEMchem internal state object
! in units of kg/kg total.

! !INTERFACE:

MODULE Chem_GridCompMod

  USE CMN_Size_Mod
  USE ESMF                                           ! ESMF library
  USE MAPL ! ewl: replace with only later
  use pFlogger, only: logging, Logger
  USE Charpak_Mod                                    ! String functions
  USE DiagList_Mod                                   ! Internal state prefixes
  USE Hco_Types_Mod, ONLY : ConfigObj
  USE Input_Opt_Mod                                  ! Input Options obj
  USE GCHP_Chunk_Mod                                 ! GCHP IRF methods
  USE GCHP_HistoryExports_Mod
  USE ErrCode_Mod                                    ! Error numbers
  USE State_Chm_Mod                                  ! Chemistry State obj
  USE State_Diag_Mod                                 ! Diagnostics State obj
  USE State_Grid_Mod                                 ! Grid State obj
  USE State_Met_Mod                                  ! Meteorology State obj
  USE Species_Mod,   ONLY : Species

  USE pflogger, ONLY : logger_t => logger

  IMPLICIT NONE
  PRIVATE

  public SetServices
  private Initialize_
  private Run_
  private Finalize_
  private Run1     ! Run wrapper phase 1
  private Run2     ! Run wrapper phase 2
  private Extract_ ! Get values from ESMF

  ! Legacy state
  TYPE GEOSCHEM_State
     PRIVATE
     TYPE(ESMF_Config)             :: myCF           ! Private ESMF Config obj
  END TYPE GEOSCHEM_State

  ! Hook for the ESMF
  TYPE GEOSCHEM_Wrap
     TYPE(GEOSCHEM_State), POINTER :: PTR => null()  ! Ptr to GEOSCHEM_State
  END TYPE GEOSCHEM_Wrap

  ! For passing from internal state to Chm_State and vice versa
  TYPE Int2SpcMap
     CHARACTER(LEN=255)            :: Name
     INTEGER                       :: ID
     REAL(ESMF_KIND_R8), POINTER   :: Internal(:,:,:) => NULL()
  END TYPE Int2SpcMap

  ! Internal run alarms
  TYPE GC_run_alarms
     private
     ! Add alarms here
     type(ESMF_Alarm) :: RRTMG_Alarm
  END TYPE GC_run_alarms

  TYPE GCRA_wrap
     type(GC_run_alarms), pointer  :: ptr
  END TYPE GCRA_wrap

  ! For mapping State_Chm%Tracers/Species arrays onto the internal state.
  TYPE(Int2SpcMap), POINTER        :: Int2Spc(:) => NULL()

  ! Objects for GEOS-Chem
  TYPE(OptInput)                   :: Input_Opt      ! Input Options
  TYPE(ChmState)                   :: State_Chm      ! Chemistry state
  TYPE(DgnState)                   :: State_Diag     ! Diagnostics state
  TYPE(GrdState)                   :: State_Grid     ! Grid state
  TYPE(MetState)                   :: State_Met      ! Meteorology state
  TYPE(Species),          POINTER  :: ThisSpc => NULL()
  TYPE(HistoryConfigObj), POINTER  :: HistoryConfig
  TYPE(ConfigObj),        POINTER  :: HcoConfig
  CLASS(Logger),          POINTER  :: lgr => null()
  LOGICAL                          :: met_wind_is_top_down
  LOGICAL                          :: met_humidity_is_top_down
  LOGICAL                          :: met_nonadv_is_top_down
  LOGICAL                          :: use_extdata2g

  ! Number of run phases, 1 or 2. Set in the rc file; else default is 2.
  INTEGER                          :: NPHASE

  ! Is this being run as a CTM?
  INTEGER                          :: IsCTM

  ! Are we reading in dynamical heating?
  LOGICAL                          :: Read_Dyn_Heating

  ! Memory debug level
  INTEGER                          :: MemDebugLevel


  ! Pointers to import, export and internal state data. Declare them as
  ! module variables so that we have to assign them only on first call.

#include "GCHPchem_DeclarePointer___.h"
  
  ! !REMARKS:
  !  Developed for GEOS-5 release Fortuna 2.0 and later.
  !                                                                             .
  !  NOTES:
  !  - The abbreviation "PET" stands for "Persistent Execution Thread".
  !    It is a synomym for CPU.
  !

CONTAINS

  !------------------------------------------------------------------------------
  !                  GEOS-Chem Global Chemical Model                            !
  !------------------------------------------------------------------------------

  ! !DESCRIPTION: The SetServices routine does the following:
  ! \item Defines the Initialize method for the GEOSCHEMchem gridded component
  ! \item Defines the Run methods for the GEOSCHEMchem gridded component
  !       (phase 1 and phase 2).
  ! \item Defines the Finalize method for the GEOSCHEMchem gridded component
  ! \item Attaches an internal state (which holds a private ESMF Config object)
  !       to the GEOSCHEMchem gridded component.
  ! \end{itemize}

  ! SetServices - Externally visible registration routine
  SUBROUTINE SetServices( GC, RC )

    USE HCOI_ESMF_MOD,        ONLY : HCO_SetServices
    USE GCKPP_Model
    USE CHARPAK_MOD,          ONLY : STRSPLIT, CSTRIP
    USE inquireMod,           ONLY : findFreeLUN
    USE FILE_MOD,             ONLY : IOERROR

    TYPE(ESMF_GridComp)  :: GC ! Composite gridded component
    INTEGER, INTENT(OUT) :: RC ! Error code, 0 all is well

    ! !REMARKS:
    !  ESMF can only attach one Config object per Gridded Component.  The
    !  Config object that is defined from the "MAPL.rc" resource file is
    !  directly attached to the GEOSCHEMchem gridded component.
    !                                                                             .
    !  To attach the Config object defined from the "GEOSCHEMchem_GridComp.rc"
    !  resource file, we must first create a derived type with a pointer to
    !  the Config object, then attach that to the gridded component as an
    !  "internal state" (also called "legacy state").

    TYPE(GEOSCHEM_State), POINTER :: myState       ! Legacy state
    TYPE(GEOSCHEM_Wrap)           :: wrap          ! Wrapper for myState
    CHARACTER(LEN=ESMF_MAXSTR)    :: compName      ! Gridded Component name
    CHARACTER(LEN=ESMF_MAXSTR)    :: COMP_NAME     ! This syntax for mapl_acg.pl
    CHARACTER(LEN=ESMF_MAXSTR)    :: HcoConfigFile ! HEMCO configuration file
    CHARACTER(LEN=ESMF_MAXSTR)    :: SpcName       ! Registered species name
    CHARACTER(LEN=40)             :: AdvSpc(500)
    CHARACTER(LEN=255)            :: LINE, MSG, SUBSTRS(500)
    INTEGER                       :: N, I, J, IU_GEOS, IOS
    INTEGER                       :: Nadv, landTypeInt
    LOGICAL                       :: FOUND
    LOGICAL                       :: EOF
    CHARACTER(LEN=60)             :: landTypeStr, importName, simType
    CHARACTER(LEN=ESMF_MAXPATHLEN):: rstFile
    INTEGER                       :: SpcRestartAttr
    CHARACTER(LEN=ESMF_MAXSTR)    :: HistoryConfigFile ! HISTORY config file
    INTEGER                       :: T

!#ifndef MAPL3
!    TYPE(MAPL_MetaComp),  POINTER :: STATE => NULL()
!#endif ! mapl3 end
    INTEGER                       :: DoIt

    ! Manual internal state entries
    LOGICAL                       :: am_I_Root
    INTEGER                       :: II
    CHARACTER(LEN=2)              :: intStr
    CHARACTER(LEN=ESMF_MAXSTR)    :: myName

    integer :: status
    class(logger_t), pointer :: logger

    call MAPL_GridCompGet(gc, logger=logger, _RC)
    call logger%debug("Chem_GridCompMod.F90::SetServices starting...")

#include "GEOSchem_Import___.h"
#include "GEOSchem_Export___.h"
#include "GEOSchem_Internal___.h"

    ! ewl: this is from advcore, but check the gridcomp that uses fv as a service
    call MAPL_GridCompAddSpec(gc, &
         state_intent=ESMF_STATEINTENT_IMPORT, &
         short_name="TRADV", &
         standard_name="advected_quantities", &
         itemtype=MAPL_STATEITEM_SERVICE, _RC)
    
    call MAPL_GridCompSetEntryPointer(gc, ESMF_Method_Initialize, Initialize, _RC)
    call MAPL_GridCompSetEntryPointer(gc, ESMF_Method_Run, Run, phase_name="Run", _RC)
    call MAPL_GridCompSetEntryPointer(gc, ESMF_Method_Finalize, Finalize, _RC)

    call logger%debug("Chem_GridCompMod.F90::SetServices done")

!    ! Set up traceback info
!    CALL ESMF_GridCompGet( GC, name=compName, __RC__ )
!
!    ! NOTE: We need to use COMP_NAME for mapl_acg.pl script
!    COMP_NAME = TRIM( compName )
!
!    ! Identify this routine to MAPL
!    Iam = TRIM(compName)//'::SetServices'

    ! Root CPU? 
    am_I_Root = MAPL_am_I_Root()    
    
    !=======================================================================
    ! Wrap internal state for storing in this gridded component
    ! Rename this to a "legacy state"
    !=======================================================================
    ALLOCATE( myState, stat=STATUS )
    _VERIFY(STATUS)
    wrap%ptr => myState

    !=======================================================================
    ! Define an ESMF Config object from the Resource file and set it
    ! as an "internal state" of the GEOSCHEMchem gridded component
    !=======================================================================
    !myState%myCF = ESMF_ConfigCreate(__RC__)
    !
    !call ESMF_ConfigLoadFile( myState%myCF, 'GCHP.rc', __RC__)
    !
    !! Get generic state object
    !CALL MAPL_GetObjectFromGC( GC, STATE, __RC__ )
    !call MAPL_GetResource( STATE, IsCTM, label='GEOSChem_CTM:', &
    !                       default=1, rc=status )
    !_VERIFY(STATUS)
    CALL MAPL_GridCompGetResource(gc, "GEOSChem_CTM", IsCTM, default=1, _RC)

!#ifndef MAPL3
!    ! Set the Initialize, Run, Finalize entry points
!    CALL MAPL_GridCompSetEntryPoint( GC, ESMF_METHOD_INITIALIZE,  &
!                                     Initialize_, __RC__ )
!    CALL MAPL_GridCompSetEntryPoint( GC, ESMF_METHOD_RUN, Run2, __RC__ )
!    CALL MAPL_GridCompSetEntryPoint( GC, ESMF_METHOD_FINALIZE,  &
!                                     Finalize_, __RC__ )
!
!    ! Store internal state with Config object in the gridded component
!    CALL ESMF_UserCompSetInternalState( GC, 'GEOSCHEM_State', wrap, STATUS )
!    _VERIFY(STATUS)
!#endif

    !=======================================================================
    ! Get meteorology vertical index orientation and ExtData version.
    ! Whether met will be flipped to be bottom-up is dependent on ExtData
    ! version (1G requires flipping here, while 2G already flipped by MAPL)
    !=======================================================================
!    call ESMF_ConfigGetAttribute(myState%myCF,value=use_extdata2g, &
!         label='USE_EXTDATA2G:', Default=.false., __RC__ )
!
!    call ESMF_ConfigGetAttribute(myState%myCF,value=met_wind_is_top_down, &
!         label='MET_WIND_IS_TOP_DOWN:', Default=.false., __RC__ )
!
!    call ESMF_ConfigGetAttribute(myState%myCF,value=met_humidity_is_top_down, &
!         label='MET_HUMIDITY_IS_TOP_DOWN:', Default=.false., __RC__ )
!
!    call ESMF_ConfigGetAttribute(myState%myCF,value=met_nonadv_is_top_down, &
!         label='MET_NONADVECTION_IS_TOP_DOWN:', Default=.false., __RC__ )
    CALL MAPL_GridCompGetResource(gc, "USE_EXTDATA2G", use_extdata2g, default=.false=., _RC)
    CALL MAPL_GridCompGetResource(gc, "MET_WIND_IS_TOP_DOWN", met_wind_is_top_down, default=.false=., _RC)
    CALL MAPL_GridCompGetResource(gc, "MET_HUMIDITY_IS_TOP_DOWN", met_humidity_is_top_down, default=.false=., _RC)
    CALL MAPL_GridCompGetResource(gc, "MET_NONADVECTION_IS_TOP_DOWN", met_nonadv_is_top_down, default=.false=., _RC)

    ! Print information to log about expectation of vertical direction of met-fields
    if ( use_extdata2g ) then
       call lgr%info('Using MAPL ExtData2G; all ''top-down'' meteorological data is automatically flipped to ''bottom-up'' within MAPL')
    else
       if (met_wind_is_top_down) then
          call lgr%info('Configured to expect ''top-down'' wind met-field imports in Chem_GridCompMod')
       else
          call lgr%info('Configured to expect ''bottom-up'' wind met-field imports in Chem_GridCompMod')
       end if
       if (met_humidity_is_top_down) then
          call lgr%info('Configured to expect ''top-down'' for humidity met-field imports in Chem_GridCompMod')
       else
          call lgr%info('Configured to expect ''bottom-up'' for humidity met-field imports in Chem_GridCompMod')
       end if
       if (met_nonadv_is_top_down) then
          call lgr%info('Configured to expect ''top-down'' for non-advection met-field imports in Chem_GridCompMod')
       else
          call lgr%info('Configured to expect ''bottom-up'' for non-advection met-field imports in Chem_GridCompMod')
       end if
    endif

!------ Species in restart file ------

    ! Determine if all species (SPC_*) are required in initial restart file
#ifdef MAPL3
    CALL ESMF_ConfigGetAttribute( myState%myCF, DoIt, &
                                  Label = "INITIAL_RESTART_SPECIES_REQUIRED:", &
                                  Default = 1, __RC__ )
#else
    CALL MAPL_GridCompGetResource(gc, "INITIAL_RESTART_SPECIES_REQUIRED", DoIt, default=1, _RC)
#endif
    IF ( DoIt == 1 ) THEN
       SpcRestartAttr  = MAPL_RestartRequired
    ELSE
       SpcRestartAttr  = MAPL_RestartOptional
    ENDIF

!-- Read in species from geoschem_config.yml and set FRIENDLYTO

    ! Open geoschem_config.yml to find the sim name and transported species
    IU_GEOS = findFreeLun()
    OPEN( IU_GEOS, FILE='geoschem_config.yml', STATUS='OLD', IOSTAT=IOS )
    IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_GEOS, 'READ_SPECIES_FROM_FILE:1' )
    DO
       READ( IU_GEOS, '(a)', IOSTAT=IOS ) LINE
       IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_GEOS, 'READ_SPECIES_FROM_FILE:2' )
       LINE = ADJUSTL( ADJUSTR( LINE ) )
       IF ( INDEX( LINE, 'name' ) > 0 ) THEN
          CALL STRSPLIT( line, ':', SUBSTRS, N )
          SimType = ADJUSTL( ADJUSTR( SUBSTRS(2) ) )
       ENDIF
       IF ( INDEX( LINE, 'transported_species' ) > 0 ) EXIT
    ENDDO

    ! Read in all advected species names and add them to internal state
    NADV = 0
    DO WHILE ( LEN_TRIM( line ) > 0 )
       READ( IU_GEOS, '(a)', IOSTAT=IOS ) LINE
       EOF = IOS < 0
       IF ( EOF ) EXIT !Simply exit when the file ends (bmy, 12 Jan 2023)
       IF ( IOS > 0 ) CALL IOERROR( IOS, IU_GEOS, 'READ_SPECIES_FROM_FILE:3' )
       LINE = ADJUSTL( ADJUSTR( LINE ) )
       IF ( INDEX( LINE, 'passive_species' ) > 0 ) EXIT
       CALL STRSPLIT( LINE, '-', SUBSTRS, N )
       IF ( INDEX( LINE, '-' ) > 0 ) THEN
          substrs(1) = ADJUSTL( ADJUSTR( substrs(1) ) )

          ! Remove quotes (i.e. 'NO' -> NO)
          J = INDEX( substrs(1), "'" )
          IF ( J > 0 ) THEN
             substrs(1) = substrs(1)(J+1:)
             J = INDEX( substrs(1), "'" )
             IF ( J > 0 ) substrs(1) = substrs(1)(1:J-1)
          ENDIF

          !%%% GEOS-Chem in GCHP ###
          CALL MAPL_AddInternalSpec(GC, &
               SHORT_NAME      = TRIM(SPFX) // TRIM(SUBSTRS(1)),            &
               LONG_NAME       = TRIM(SUBSTRS(1)),                          &
               UNITS           = 'mol mol-1',                               &
               DIMS            = MAPL_DimsHorzVert,                         &
               VLOCATION       = MAPL_VLocationCenter,                      &
               PRECISION       = ESMF_KIND_R8,                              &
               FRIENDLYTO      = 'DYNAMICS:TURBULENCE:MOIST',               &
               RESTART         = SpcRestartAttr,                               &
               RC              = RC                                       )

          ! Add to list of transported speces
          NADV = NADV + 1
          AdvSpc(NADV) = TRIM(SUBSTRS(1))
       ENDIF
    ENDDO
    CLOSE( IU_GEOS )

!-- Add all non-advected species from KPP-based simulations
!-- (but don't add dummy species).  KPP-based simulations now
!-- include fullchem, Hg, and carbon.
    IF ( TRIM( simType ) == 'fullchem'      .or.                            &
         TRIM( simType ) == 'Hg'            .or.                            &
         TRIM( simType ) == 'carbon' ) THEN
       DO I=1,NSPEC
          FOUND = .false.


          ! Skip dummy RR species for prod/loss diagnostic (mps, 8/23/16)
          SpcName = ADJUSTL( Spc_Names(I) )
          IF ( SpcName(1:2) == 'RR' ) CYCLE

          DO J=1,Nadv !Size of AdvSpc
             IF (trim(AdvSpc(J)) .eq. trim(SpcName)) THEN
                FOUND = .true.
                EXIT
             ENDIF
          END DO

          ! Add non-advected species to internal state
          IF ( .NOT. Found ) THEN
             !%%%% GEOS-Chem in GCHP %%%%
             call MAPL_AddInternalSpec(GC, &
                  SHORT_NAME      = TRIM(SPFX) // SpcName,                   &
                  LONG_NAME       = SpcName,                                 &
                  UNITS           = 'mol mol-1',                             &
                  PRECISION       = ESMF_KIND_R8,                            &
                  DIMS            = MAPL_DimsHorzVert,                       &
                  VLOCATION       = MAPL_VLocationCenter,                    &
                  RESTART         = SpcRestartAttr,                          &
                  RC              = STATUS                                  )
          ENDIF
       ENDDO
    ENDIF

    ! Add other internal state variables as real8 for GCHP

    call MAPL_AddInternalSpec(GC, &
       SHORT_NAME         = 'DryDepNitrogen',  &
       LONG_NAME          = 'Dry deposited nitrogen',  &
       UNITS              = 'cm-2s-1', &
       DIMS               = MAPL_DimsHorzOnly,    &
       VLOCATION          = MAPL_VLocationCenter,    &
       PRECISION          = ESMF_KIND_R8, &
       FRIENDLYTO         = trim(COMP_NAME),    &
                                                      RC=STATUS  )
    _VERIFY(STATUS)

    call MAPL_AddInternalSpec(GC, &
       SHORT_NAME         = 'WetDepNitrogen',  &
       LONG_NAME          = 'Wet deposited nitrogen',  &
       UNITS              = 'cm-2s-1', &
       DIMS               = MAPL_DimsHorzOnly,    &
       VLOCATION          = MAPL_VLocationCenter,    &
       PRECISION          = ESMF_KIND_R8, &
       FRIENDLYTO         = trim(COMP_NAME),    &
                                                      RC=STATUS  )
    _VERIFY(STATUS)

    call MAPL_AddInternalSpec(GC, &
       SHORT_NAME         = 'H2O2AfterChem',  &
       LONG_NAME          = 'Soluble fraction H2O2',  &
       UNITS              = 'vv-1', &
       DIMS               = MAPL_DimsHorzVert,    &
       VLOCATION          = MAPL_VLocationCenter,    &
       PRECISION          = ESMF_KIND_R8, &
       FRIENDLYTO         = trim(COMP_NAME),    &
                                                      RC=STATUS  )
    _VERIFY(STATUS)

    call MAPL_AddInternalSpec(GC, &
       SHORT_NAME         = 'SO2AfterChem',  &
       LONG_NAME          = 'Soluble fraction SO2',  &
       UNITS              = 'vv-1', &
       DIMS               = MAPL_DimsHorzVert,    &
       VLOCATION          = MAPL_VLocationCenter,    &
       PRECISION          = ESMF_KIND_R8, &
       FRIENDLYTO         = trim(COMP_NAME),    &
                                                      RC=STATUS  )
    _VERIFY(STATUS)

    call MAPL_AddInternalSpec(GC, &
       SHORT_NAME         = 'KPPHvalue',  &
       LONG_NAME          = 'HSAVE for KPP',  &
       UNITS              = '1', &
       DIMS               = MAPL_DimsHorzVert,    &
       VLOCATION          = MAPL_VLocationCenter,    &
       PRECISION          = ESMF_KIND_R8, &
       FRIENDLYTO         = trim(COMP_NAME),    &
                                                      RC=STATUS  )
    _VERIFY(STATUS)

    ! Sulfur-nitrogen-ammonia water content computed in Isorropia/HETP
    ! after needed in RDAER
    call MAPL_AddInternalSpec(GC, &
       SHORT_NAME         = 'AeroH2O_SNA',  &
       LONG_NAME          = 'Sulfur-nitrogen-ammonia water content',  &
       UNITS              = 'g/m3', &
       DIMS               = MAPL_DimsHorzVert,    &
       VLOCATION          = MAPL_VLocationCenter,    &
       PRECISION          = ESMF_KIND_R8, &
       FRIENDLYTO         = trim(COMP_NAME),    &
                                                      RC=STATUS  )
    _VERIFY(STATUS)

    ! Sesquiterpene mass per grid box
    call MAPL_AddInternalSpec(GC, &
       SHORT_NAME         = 'ORVCSESQ',  &
       LONG_NAME          = 'Sesquiterpenes mass',  &
       UNITS              = 'kg', &
       DIMS               = MAPL_DimsHorzVert,    &
       VLOCATION          = MAPL_VLocationCenter,    &
       PRECISION          = ESMF_KIND_R8, &
       FRIENDLYTO         = trim(COMP_NAME),    &
                                                      RC=STATUS  )
    _VERIFY(STATUS)

    ! Surface J-values for HEMCO
    call MAPL_AddInternalSpec(GC, &
       SHORT_NAME         = 'JOH',  &
       LONG_NAME          = 'Surface J-values for reaction O3 + hv --> O2 + O',  &
       UNITS              = '1', &
       DIMS               = MAPL_DimsHorzOnly,    &
       PRECISION          = ESMF_KIND_R8, &
       FRIENDLYTO         = trim(COMP_NAME),    &
                                                      RC=STATUS  )
    _VERIFY(STATUS)

    call MAPL_AddInternalSpec(GC, &
       SHORT_NAME         = 'JNO2',  &
       LONG_NAME          = 'Surface J-values for reaction NO2 + hv --> NO + O',  &
       UNITS              = '1', &
       DIMS               = MAPL_DimsHorzOnly,    &
       PRECISION          = ESMF_KIND_R8, &
       FRIENDLYTO         = trim(COMP_NAME),    &
                                                      RC=STATUS  )
    _VERIFY(STATUS)

    ! delta dry pressure used to conserve mass across consecutive runs
    call MAPL_AddInternalSpec(GC, &
       SHORT_NAME         = 'DELP_DRY',  &
       LONG_NAME          = 'Delta dry pressure across box',  &
       UNITS              = 'hPa', &
       DIMS               = MAPL_DimsHorzVert,    &
       VLOCATION          = MAPL_VLocationCenter,    &
       PRECISION          = ESMF_KIND_R8, &
       FRIENDLYTO         = trim(COMP_NAME),    &
                                                      RC=STATUS  )
    _VERIFY(STATUS)

    ! Additional outputs useful for unit conversions and post-processing analysis
    call MAPL_AddInternalSpec(GC, &
       SHORT_NAME         = 'AREA',  &
       LONG_NAME          = 'Grid horizontal area',  &
       UNITS              = 'm2', &
       DIMS               = MAPL_DimsHorzOnly,    &
       PRECISION          = ESMF_KIND_R8, &
       FRIENDLYTO         = trim(COMP_NAME),    &
                                                      RC=STATUS  )
    _VERIFY(STATUS)
    call MAPL_AddInternalSpec(GC, &
       SHORT_NAME         = 'BXHEIGHT',  &
       LONG_NAME          = 'Grid box height (w/r/t dry air)',  &
       UNITS              = 'm', &
       DIMS               = MAPL_DimsHorzVert,    &
       VLOCATION          = MAPL_VLocationCenter,    &
       PRECISION          = ESMF_KIND_R8, &
       FRIENDLYTO         = trim(COMP_NAME),    &
                                                      RC=STATUS  )
    _VERIFY(STATUS)

    call MAPL_AddInternalSpec(GC, &
       SHORT_NAME         = 'TropLev',  &
       LONG_NAME          = 'GEOS-Chem level where the tropopause occurs',  &
       UNITS              = '1', &
       DIMS               = MAPL_DimsHorzOnly,    &
       VLOCATION          = MAPL_VLocationCenter,    &
       PRECISION          = ESMF_KIND_R8, &
       FRIENDLYTO         = trim(COMP_NAME),    &
                                                      RC=STATUS  )
    _VERIFY(STATUS)


!
! !EXTERNAL STATE:
!
#   include "GCHPchem_ExportSpec___.h"

    ! Read HISTORY config file and add exports for unique items
    CALL ESMF_ConfigGetAttribute( myState%myCF, HistoryConfigFile, &
                                  Label="HISTORY_CONFIG:",         &
                                  Default="HISTORY.rc", __RC__ )
    CALL HistoryExports_SetServices( MAPL_am_I_Root(), HistoryConfigFile, &
                                     GC, HistoryConfig, __RC__ )

!EOP
!BOC


    ! OLSON
    DO T = 1, NSURFTYPE
       landTypeInt = T-1
       WRITE ( landTypeStr, '(I2.2)' ) landTypeInt
       importName = 'OLSON' // TRIM(landTypeStr)
       CALL MAPL_AddImportSpec(GC,                                  &
          SHORT_NAME         = importName,                          &
          LONG_NAME          = 'OLSON_land_by_type',                &
          UNITS              = 'unitless',                          &
          DIMS               = MAPL_DimsHorzOnly,                   &
          RESTART            = MAPL_RestartSkip,                    &
                                                            __RC__ )
    ENDDO

    ! Set HEMCO services
    ! --------------------
    CALL ESMF_ConfigGetAttribute( myState%myCF, HcoConfigFile, &
                                  Label="HEMCO_CONFIG:", &
                                  Default="HEMCO_Config.rc", __RC__ )
    CALL HCO_SetServices( MAPL_am_I_Root(), GC, HcoConfig,  &
                          TRIM(HcoConfigFile), __RC__ )

    ! Set the Profiling timers
    ! ------------------------
    CALL MAPL_TimerAdd(GC, NAME="INITIALIZE", RC=status)
    _VERIFY(status)
    CALL MAPL_TimerAdd(GC, NAME="RUN", RC=status)
    _VERIFY(status)
    CALL MAPL_TimerAdd(GC, NAME="FINALIZE", RC=status)
    _VERIFY(status)

    CALL MAPL_TimerAdd(GC, NAME="DO_CHEM", RC=status)
    _VERIFY(status)
    CALL MAPL_TimerAdd(GC, NAME="CP_BFRE", RC=status)
    _VERIFY(status)
    CALL MAPL_TimerAdd(GC, NAME="CP_AFTR", RC=status)
    _VERIFY(status)

    ! More timers to be called in gchp_chunk_run
    CALL MAPL_TimerAdd(GC, NAME="GC_CONV"  , __RC__)
    CALL MAPL_TimerAdd(GC, NAME="GC_EMIS"  , __RC__)
    CALL MAPL_TimerAdd(GC, NAME="GC_DRYDEP", __RC__)
    CALL MAPL_TimerAdd(GC, NAME="GC_FLUXES", __RC__)
    CALL MAPL_TimerAdd(GC, NAME="GC_TURB"  , __RC__)
    CALL MAPL_TimerAdd(GC, NAME="GC_CHEM"  , __RC__)
    CALL MAPL_TimerAdd(GC, NAME="GC_WETDEP", __RC__)
    CALL MAPL_TimerAdd(GC, NAME="GC_DIAGN" , __RC__)
    CALL MAPL_TimerAdd(GC, NAME="GC_RAD"   , __RC__)

    ! Generic Set Services
    ! --------------------
    CALL MAPL_GenericSetServices( GC, RC=status )
    _VERIFY(status)

    !=======================================================================
    ! All done
    !=======================================================================
    _RETURN(ESMF_SUCCESS)

  END SUBROUTINE SetServices
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Model                            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialize_
!
! !DESCRIPTION: Initialize_ is the initialize method of the GEOSCHEMchem
!  gridded component.  This is a simple ESMF/MAPL wrapper which calls down
!  to the Initialize method of the GEOS-Chem column chemistry code.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Initialize_( GC, Import, Export, Clock, RC )
!
! !USES:
!
    USE TIME_MOD,  ONLY : GET_TS_CHEM, GET_TS_EMIS
    USE TIME_MOD,  ONLY : GET_TS_DYN,  GET_TS_CONV
    USE TIME_MOD,  ONLY : GET_TS_RAD
    Use pfLogger,  ONLY : Logger
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ESMF_GridComp), INTENT(INOUT)         :: GC       ! Ref to GridComp
    TYPE(ESMF_State),    INTENT(INOUT), TARGET :: Import   ! Import State object
    TYPE(ESMF_State),    INTENT(INOUT), TARGET :: Export   ! Export State object
    TYPE(ESMF_Clock),    INTENT(INOUT)         :: Clock    ! ESMF clock object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,             INTENT(OUT)           :: RC       ! Success or failure?
!
! !REMARKS:
!  We call routine Extract_ to return various values (i.e. grid parameters,
!  start & end dates, PET information, etc.) from the ESMF/MAPL environment.
!  We then pass those to GEOS-Chem via routine GCHP_CHUNK_INIT, which is
!  located in GEOS-Chem module ./GEOS-Chem/ESMF/gchp_chunk_mod.F90.
!
! !REVISION HISTORY:
!  06 Dec 2009 - A. da Silva - Initial version
!  See https://github.com/geoschem/geos-chem for history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! LOCAL VARIABLES:
!
    ! Objects
    TYPE(ESMF_Grid)             :: Grid        ! ESMF Grid object
    TYPE(ESMF_Config)           :: MaplCF      ! ESMF Config obj (MAPL.rc)
    TYPE(ESMF_Config)           :: GeosCF      ! ESMF Config obj (GEOSCHEM*.rc)

    ! Scalars
    LOGICAL                     :: am_I_Root   ! Are we on the root PET?
    INTEGER                     :: myPet       ! # of the PET we are on
    INTEGER                     :: NPES        ! # of total PETs in MPI world
    INTEGER                     :: nymdB       ! GMT date @ start of simulation
    INTEGER                     :: nymdE       ! GMT date @ end of simulation
    INTEGER                     :: nymd        ! GMT date @ current time
    INTEGER                     :: nhmsB       ! GMT time @ start of simulation
    INTEGER                     :: nhmsE       ! GMT time @ end of simulation
    INTEGER                     :: nhms        ! GMT time @ current time
    INTEGER                     :: IM          ! # of longitudes on this PET
    INTEGER                     :: JM          ! # of latitudes  on this PET
    INTEGER                     :: LM          ! # of levels     on this PET
    INTEGER                     :: value_LLSTRAT ! # of strat. levels
    INTEGER                     :: IM_WORLD    ! # of longitudes in global grid
    INTEGER                     :: JM_WORLD    ! # of latitudes  in global grid
    INTEGER                     :: LM_WORLD    ! # of levels     in global grid
    REAL                        :: tsChem      ! Chemistry timestep [s]
    REAL                        :: tsDyn       ! Dynamic timestep [s]
    REAL                        :: tsRad       ! RRTMG timestep [s]
    CHARACTER(LEN=5)            :: petStr      ! String for PET #
    CHARACTER(LEN=ESMF_MAXSTR)  :: compName    ! Name of gridded component

    ! time step error checks
    REAL                         :: ChemTS, EmisTS, RadTS

    ! Pointer arrays
    REAL(ESMF_KIND_R4),  POINTER :: lonCtr(:,:) ! Lon centers on this PET [rad]
    REAL(ESMF_KIND_R4),  POINTER :: latCtr(:,:) ! Lat centers on this PET [rad]

    INTEGER                      :: I, J, nFlds, mpiComm
    TYPE(ESMF_STATE)             :: INTSTATE
    TYPE(ESMF_Field)             :: GcFld

    ! Species information
    TYPE(Species), POINTER       :: SpcInfo

    CHARACTER(LEN=ESMF_MAXSTR)   :: fieldName

    INTEGER                      :: N, trcID
    TYPE(MAPL_MetaComp), POINTER :: STATE => NULL()
    REAL(ESMF_KIND_R8), POINTER  :: Ptr3D(:,:,:) => NULL()
    REAL(ESMF_KIND_R8), POINTER  :: Ptr3D_int(:,:,:) => NULL()
    REAL(ESMF_KIND_R8), POINTER  :: Ptr3D_exp(:,:,:) => NULL()

    ! Internal run alarms
    type(GC_run_alarms), pointer :: GC_alarms
    type(GCRA_wrap)              :: GC_alarm_wrapper
    TYPE(ESMF_Time)              :: currTime       ! Current (start) time
    TYPE(ESMF_Time)              :: ringTime
    type(ESMF_TimeInterval)      :: tsRad_TI
    type(ESMF_TimeInterval)      :: tsChem_TI
    type (ESMF_Calendar)         :: CAL
    INTEGER                      :: yyyy, mm, dd   ! Year, month, day
    INTEGER                      :: h,    m,  s    ! Hour, minute, seconds
    INTEGER                      :: doy

    INTEGER                     :: IL_WORLD, JL_WORLD    ! # lower indices in global grid
    INTEGER                     :: IU_WORLD, JU_WORLD    ! # upper indices in global grid

    !=======================================================================
    ! Initialization
    !=======================================================================

    ! Get my name and set-up traceback handle
    CALL ESMF_GridCompGet( GC, name=compName, __RC__ )

    ! Identify this routine to MAPL
    Iam = TRIM(compName)//'::Initialize_'

    ! Get my MAPL_Generic state
    ! -------------------------
    CALL MAPL_GetObjectFromGC(GC, STATE, RC=STATUS)
    _VERIFY(STATUS)

    !  Start timers
    !  ------------
    CALL MAPL_TimerOn( STATE, "TOTAL")
    CALL MAPL_TimerOn( STATE, "INITIALIZE")

    ! Initialize MAPL Generic
    CALL MAPL_GenericInitialize( GC, Import, Export, Clock, __RC__ )

    ! Get Internal state.
    CALL MAPL_Get ( STATE, INTERNAL_ESMF_STATE=INTSTATE, __RC__ )

    ! Initialize GEOS-Chem Input_Opt fields to zeros or equivalent
    CALL Set_Input_Opt( MAPL_am_I_Root(), Input_Opt, RC )
    _ASSERT(RC==GC_SUCCESS, 'Error calling Set_Input_Opt')

    ! Grab the logger for this component
    call MAPL_GetLogger(GC, Input_Opt%lgr, __RC__)
    Input_Opt%compname = Trim(compname)

    ! Root CPU?
    am_I_Root = MAPL_am_I_Root()

    ! Get various parameters from the ESMF/MAPL framework
    CALL Extract_( GC,                        &  ! Ref to this Gridded Comp
                   Clock,                     &  ! ESMF Clock object
                   Grid        = Grid,        &  ! ESMF Grid object
                   MaplCF      = MaplCF,      &  ! AGCM.rc/GCHP.rc config object
                   GeosCF      = GeosCF,      &  ! GEOSCHEM*.rc Config object
                   IM          = IM,          &  ! # of longitudes on this PET
                   JM          = JM,          &  ! # of latitudes  on this PET
                   LM          = LM,          &  ! # of levels     on this PET
                   IM_WORLD    = IM_WORLD,    &  ! # of lons in global grid
                   JM_WORLD    = JM_WORLD,    &  ! # of lats  in global grid
                   LM_WORLD    = LM_WORLD,    &  ! # of levels in global grid
                   IL_WORLD    = IL_WORLD,    &  ! start index of lons in global grid on this PET
                   IU_WORLD    = IU_WORLD,    &  ! end   index of lons in global grid on this PET
                   JL_WORLD    = JL_WORLD,    &  ! start index of lats in global grid on this PET
                   JU_WORLD    = JU_WORLD,    &  ! end   index of lats in global grid on this PET
                   nymdB       = nymdB,       &  ! YYYYMMDD @ start of sim
                   nhmsB       = nhmsB,       &  ! hhmmss   @ end   of sim
                   nymdE       = nymdE,       &  ! YYYMMDD  @ start of sim
                   nhmsE       = nhmsE,       &  ! hhmmss   @ end   of sim
                   tsChem      = tsChem,      &  ! Chemistry timestep [seconds]
                   tsRad       = tsRad,       &  ! RRTMG timestep [seconds]
                   tsDyn       = tsDyn,       &  ! Dynamics timestep  [seconds]
                   localPet    = myPet,       &  ! PET # that we are on now
                   petCount    = NPES,        &  ! Number of PETs in MPI World
                   mpiComm     = mpiComm,     &  ! MPI Communicator Handle
                   lonCtr      = lonCtr,      &  ! This PET's lon ctrs [radians]
                   latCtr      = latCtr,      &  ! This PET's lat ctrs [radians]
                   __RC__                      )

    ! Set MPI values in Input_Opt
    Input_Opt%thisCPU = myPet
    Input_Opt%MPIComm = mpiComm
    Input_Opt%numCPUs = NPES
    Input_Opt%isMPI   = .true.
    if ( MAPL_am_I_Root() ) Input_Opt%amIRoot = .true.

    ! MSL - shift from 0 - 360 to -180 - 180 degree grid
    where (lonCtr .gt. MAPL_PI ) lonCtr = lonCtr - 2*MAPL_PI

    ! Get the memory debug level
    call ESMF_ConfigGetAttribute(GeosCF, MemDebugLevel, &
                                 Label="MEMORY_DEBUG_LEVEL:" , RC=STATUS)
    _VERIFY(STATUS)

    !=======================================================================
    ! Save values from the resource file (GCHP.rc for GCHP)
    !=======================================================================

    ! # of run phases
    CALL ESMF_ConfigGetAttribute( GeosCF, NPHASE,                   &
                                  Default = 2,                      &
                                  Label   = "RUN_PHASES:",          &
                                  __RC__                           )
    _ASSERT(NPHASE==1.OR.NPHASE==2,'Error calling ESMF_ConfigGetAttribute on RUN_PHASES')

    !=======================================================================
    ! Initialize GEOS-Chem (will also initialize HEMCO)
    !=======================================================================

    ! Initialize fields of the Grid State object
    CALL Init_State_Grid( Input_Opt, State_Grid, RC )
    _ASSERT(RC==GC_SUCCESS,'Error calling Init_State_Grid')

    ! Pass grid information obtained from Extract_ to State_Grid
    State_Grid%NX          = IM            ! # lons   on this PET
    State_Grid%NY          = JM            ! # lats   on this PET
    State_Grid%NZ          = LM            ! # levels on this PET
    State_Grid%GlobalNX    = IM_WORLD      ! # lons   in global grid
    State_Grid%GlobalNY    = JM_WORLD      ! # lats   in global grid
    State_Grid%NativeNZ    = LM_WORLD      ! # levels in global grid
    State_Grid%XMinOffset  = 1             ! X offset from global grid
    State_Grid%XMaxOffset  = State_Grid%NX ! X offset from global grid
    State_Grid%YMinOffset  = 1             ! Y offset from global grid
    State_Grid%YMaxOffset  = State_Grid%NY ! Y offset from global grid

    ! Call the GCHP initialize routine
    CALL GCHP_Chunk_Init( nymdB     = nymdB,      & ! YYYYMMDD @ start of run
                          nhmsB     = nhmsB,      & ! hhmmss   @ start of run
                          nymdE     = nymdE,      & ! YYYYMMDD @ end of run
                          nhmsE     = nhmsE,      & ! hhmmss   @ end of run
                          tsChem    = tsChem,     & ! Chemical timestep [s]
                          tsDyn     = tsDyn,      & ! Dynamic  timestep [s]
                          tsRad     = tsRad,      & ! RRTMG    timestep [s]
                          lonCtr    = lonCtr,     & ! Lon centers [radians]
                          latCtr    = latCtr,     & ! Lat centers [radians]
                          GC        = GC,         & ! Ref to this gridded comp
                          EXPORT    = EXPORT,     & ! Export state object
                          Input_Opt = Input_Opt,  & ! Input Options obj
                          State_Chm = State_Chm,  & ! Chemistry State obj
                          State_Diag= State_Diag, & ! Diagnostics State obj
                          State_Grid= State_Grid, & ! Grid State obj
                          State_Met = State_Met,  & ! Meteorology State obj
                          HcoConfig = HcoConfig,  & ! HEMCO config obj
                          HistoryConfig = HistoryConfig, & ! History Config Obj
                          __RC__                 )

    !=======================================================================
    ! Initialize the Int2Spc object. This is used to copy the tracer arrays
    ! from the internal state to State_Chm%Tracers, and vice versa.
    ! In this step, we also check for the friendlieness of the tracers. If
    ! the GEOS-Chem internal convection/turbulence schemes shall be used
    ! (as specified in geoschem_config.yml), the tracers must not be friendly
    ! to the GEOS-5 moist / turbulence components!
    !=======================================================================
    nFlds = State_Chm%nSpecies
    ALLOCATE( Int2Spc(nFlds), STAT=STATUS )
    _ASSERT(STATUS==0,'Int2Spc could not be allocated')

    ! Do for every tracer in State_Chm
    DO I = 1, nFlds

       SpcInfo => State_Chm%SpcData(I)%Info

       ! Pass tracer name
       Int2Spc(I)%Name = TRIM(SpcInfo%Name)

       ! Get tracer ID
       Int2Spc(I)%ID = IND_( TRIM(Int2Spc(I)%Name) )

       ! If tracer ID is not valid, make sure all vars are at least defined.
       IF ( Int2Spc(I)%ID <= 0 ) THEN
          Int2Spc(I)%Internal => NULL()
          CYCLE
       ENDIF

       ! Get internal state field
       fieldName = TRIM(SPFX)//TRIM(Int2Spc(I)%Name)
       CALL ESMF_StateGet( INTSTATE, TRIM(fieldName), GcFld, RC=STATUS )

       ! This is mostly for testing
       IF ( STATUS /= ESMF_SUCCESS ) THEN
          IF( am_I_Root ) THEN
             WRITE(*,*) 'Cannot find in internal state: ', TRIM(SPFX) &
                        //TRIM(Int2Spc(I)%Name),I
          ENDIF
          Int2Spc(I)%Internal => NULL()
       ENDIF

       ! Get pointer to field
       CALL ESMF_FieldGet( GcFld, 0, Ptr3D, __RC__ )
       Int2Spc(I)%Internal => Ptr3D
       Ptr3D => NULL()
       SpcInfo => NULL()

    ENDDO

    !=======================================================================
    ! Error trap: make sure that chemistry / emission time step are same and
    ! correspond to the chemistry step set in GEOSCHEMchem_GridComp.rc.
    !=======================================================================
    ChemTS = GET_TS_CHEM()
    EmisTS = GET_TS_EMIS()
    IF ( ChemTS /= tsChem .OR. EmisTS /= tsChem ) THEN
       WRITE(*,*) 'GEOS-Chem chemistry and/or emission time step do not'
       WRITE(*,*) 'agree with time step set in GEOSCHEMchem_GridComp.rc'
       WRITE(*,*) 'GEOS-Chem chemistry time step                 : ', ChemTS
       WRITE(*,*) 'GEOS-Chem emission  time step                 : ', EmisTS
       WRITE(*,*) 'CHEMISTRY_TIMESTEP in GCHP.rc                 : ', tsChem
       _ASSERT(.FALSE.,'Error in timesteps')
    ENDIF

    ! Also check for convection and dynamics time step.
    ChemTS = GET_TS_CONV()
    EmisTS = GET_TS_DYN()
    IF ( ChemTS /= tsDyn .OR. EmisTS /= tsDyn ) THEN
       WRITE(*,*) 'GEOS-Chem transport and/or convection time step do not'
       WRITE(*,*) 'agree with time step set in GEOSCHEMchem_GridComp.rc'
       WRITE(*,*) 'GEOS-Chem convection time step                : ', ChemTS
       WRITE(*,*) 'GEOS-Chem dynamics   time step                : ', EmisTS
       WRITE(*,*) 'RUN_DT in CAP.rc                              : ', tsDyn
       _ASSERT(.FALSE.,'Error in timesteps')
    ENDIF

    If (Input_Opt%LRAD) Then
       RadTS  = GET_TS_RAD()
       IF ( RadTS /= tsRad ) THEN
          WRITE(*,*) 'GEOS-Chem radiation time step (for RRTMG) does not'
          WRITE(*,*) 'agree with time step set in GCHP.rc'
          WRITE(*,*) 'GEOS-Chem RRTMG time step                     : ', RadTS
          WRITE(*,*) 'RRTMG_DT in GCHP.rc                           : ', tsRad
          _ASSERT(.FALSE.,'Error in timesteps')
       ENDIF

       ! Redundantly, check that tsRad is a multiple of tsChem
       _ASSERT(MOD(tsRad,tsChem)==0,'Radiation time step must be a multiple of chemistry time step')
    Else
       ! Use chemistry step; this alarm will be ignored, but must be present
       RadTS = ChemTS
    End If

    !=======================================================================
    ! Establish the internal alarms for GEOS-Chem
    !=======================================================================
    allocate(GC_alarms,stat=status)
    _ASSERT(rc==0,'Could not allocate GC alarms')
    GC_alarm_wrapper%ptr => GC_alarms

    call ESMF_UserCompSetInternalState(GC,'gcchem_internal_alarms',GC_alarm_wrapper,status)
    _ASSERT(status==0,'Could not get GEOS-Chem internal alarms')

    ! Get information about/from the clock
    CALL ESMF_ClockGet( Clock,                    &
                        currTime     = currTime,  &
                        calendar     = cal,       &
                        __RC__ )


    ! Set up the radiation alarm
    ! Must ring once per tsRad
    call ESMF_TimeIntervalSet(tsRad_TI, S=nint(tsRad), calendar=cal, RC=STATUS)
    _ASSERT(STATUS==0,'Could not set radiation alarm time interval')

    ! Initialize the ring time to midnight on the starting (current) day
    call ESMF_TimeGet( currTime, YY=yyyy, MM=mm, DD=dd, H=h, M=m, S=s, rc=STATUS )
    _ASSERT(STATUS==0,'Could not extract ESMF clock current time information')
    call ESMF_TimeSet( ringTime, YY=yyyy, MM=mm, DD=dd, H=0, M=0, S=0, rc=STATUS )
    _ASSERT(STATUS==0,'Could not set initial radiation alarm ring time')

    ! Adjust the alarm to go off on the chemistry timestep immediately before the
    ! target output time. This is because RRTMG is run after chemistry.
    call ESMF_TimeIntervalSet(tsChem_TI, S=nint(tsChem), calendar=cal, RC=STATUS)
    _ASSERT(STATUS==0,'Could not set chemistry alarm time interval')
    ringTime = ringTime - tsChem_TI

    ! Advance ring time until it is at or after current time
    do while (ringTime < currTime)
       ringTime = ringTime + tsRad_TI
    end do

    ! Make the alarm 'sticky'. This means it will ring until
    ! the ringer is turned off.
    GC_alarms%RRTMG_alarm = ESMF_AlarmCreate(CLOCK = Clock, &
                            name = "GC_RRTMG_alarm" ,       &
                            RingInterval = tsRad_TI,        &
                            RingTime     = ringTime,        &
!                            Enabled      = .true.   ,       &
                            sticky       = .true.,          &
                            RC           = STATUS      )
    _VERIFY(STATUS)

    ! Start alarm ringing if already reached first alarm time
    if(ringTime == currTime) then
       call ESMF_AlarmRingerOn(GC_alarms%RRTMG_alarm, rc=status)
       _VERIFY(STATUS)
    end if

    ! Stop timers
    ! -----------
    CALL MAPL_TimerOff( STATE, "INITIALIZE")

    CALL MAPL_TimerOff( STATE, "TOTAL")

    ! Successful return
    _RETURN(ESMF_SUCCESS)

    ! Formats
100 FORMAT( '### ',                                           / &
            '### ', a ,                                       / &
            '### ', a, '  |  Initialization on PET # ', i5.5, / &
            '### ' )
200 FORMAT( '### ',                                           / &
            '### ', a, '  |  Execution on PET # ',      i5.5, / &
            '###' )

  END SUBROUTINE Initialize_

!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Model                            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Run1
!
! !DESCRIPTION: Run1 is a wrapper method for the phase 1 run phase of the
!  GEOSCHEMchem gridded component. It calls down to the Run method of the
!  GEOS-Chem column chemistry code.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Run1 ( GC, Import, Export, Clock, RC )
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ESMF_GridComp), INTENT(INOUT) :: GC       ! Ref to this GridComp
    TYPE(ESMF_State),    INTENT(INOUT) :: Import   ! Import State
    TYPE(ESMF_State),    INTENT(INOUT) :: Export   ! Export State
    TYPE(ESMF_Clock),    INTENT(INOUT) :: Clock    ! ESMF Clock object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,             INTENT(OUT)   :: RC       ! Error return code
!
! !REMARKS:
!
! !REVISION HISTORY:
!  22 Sep 2014 - C. Keller   - Initial version.
!  See https://github.com/geoschem/geos-chem for history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    CHARACTER(LEN=ESMF_MAXSTR)  :: compName    ! Name of gridded component
    CHARACTER(LEN=ESMF_MAXSTR)  :: Iam
    INTEGER                     :: STATUS
    INTEGER                     :: PHASE

    !=======================================================================
    ! Run1 starts here
    !=======================================================================

    ! Set up traceback info
    CALL ESMF_GridCompGet( GC, name=compName, __RC__ )

    ! Identify this routine to MAPL
    Iam = TRIM(compName)//'::Run1'

    ! Call run routine stage 1 if more than one phase. If not 2 phases,
    ! such as in GCHP, then we do all chemistry related processes from
    ! Run2 instead.
    IF ( NPHASE == 2 ) THEN
       PHASE = 1
       CALL Run_ ( GC, IMPORT, EXPORT, CLOCK, PHASE, __RC__ )
    ENDIF

    ! Return w/ success
    _RETURN(ESMF_SUCCESS)

  END SUBROUTINE Run1
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Model                            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Run2
!
! !DESCRIPTION: Run2 is a wrapper method for the phase 2 run phase of the
!  GEOSCHEMchem gridded component. It calls down to the Run method of the
!  GEOS-Chem column chemistry code.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Run2 ( GC, Import, Export, Clock, RC )
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ESMF_GridComp), INTENT(INOUT) :: GC       ! Ref to this GridComp
    TYPE(ESMF_State),    INTENT(INOUT) :: Import   ! Import State
    TYPE(ESMF_State),    INTENT(INOUT) :: Export   ! Export State
    TYPE(ESMF_Clock),    INTENT(INOUT) :: Clock    ! ESMF Clock object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,             INTENT(OUT)   :: RC       ! Error return code
!
! !REMARKS:
!
! !REVISION HISTORY:
!  22 Sep 2014 - C. Keller   - Initial version.
!  See https://github.com/geoschem/geos-chem for history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    CHARACTER(LEN=ESMF_MAXSTR)  :: compName    ! Name of gridded component
    CHARACTER(LEN=ESMF_MAXSTR)  :: Iam
    INTEGER                     :: PHASE
    INTEGER                     :: STATUS

    !=======================================================================
    ! Run2 starts here
    !=======================================================================

    ! Set up traceback info
    CALL ESMF_GridCompGet( GC, name=compName, __RC__ )

    ! Identify this routine to MAPL
    Iam = TRIM(compName)//'::Run2'

    ! Set phase number: this is 2 for multi-phase runs (e.g. GEOS-5), and
    ! is -1 for single-phase runs (e.g. GCHP). If set to -1, all processes
    ! are called (drydep, emissions, chemistry, etc.)
    IF ( NPHASE == 1 ) THEN
       PHASE = -1
    ELSE
       PHASE = 2
    ENDIF

    ! Call run routine stage 2
    CALL Run_ ( GC, IMPORT, EXPORT, CLOCK, PHASE, __RC__ )

    ! Optional timer for run 2 (psturm, April 2024)
    ! More realistic timing estimates can detect load imbalances in chemistry
    ! Rather than attributing this to transport operations
    ! This can be set from GEOSCHEMchem_GridComp.rc

    ! Return w/ success
    _RETURN(ESMF_SUCCESS)

  END SUBROUTINE Run2
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Model                            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Run_
!
! !DESCRIPTION: Run_ is the run method of the GEOSCHEMchem gridded component.
!  GC is a simple ESMF/MAPL wrapper which calls down to the Run method of
!  the GEOS-Chem column chemistry code.
!  Note: this routine currently skips the call down to GEOS-Chem on the very
!  first time it is invoked. The reason is that a number of met-variables seem
!  to be undefined still (e.g. BXHEIGHT, T, etc), yielding to seg-faults and/or
!  crazy results.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Run_( GC, Import, Export, Clock, Phase, RC )
!
! !USES:
!
    USE CMN_Size_Mod,            ONLY : NDUST
    USE HCO_State_GC_Mod,        ONLY : HcoState
    USE MAPL_MemUtilsMod
    USE Olson_Landmap_Mod,       ONLY : Compute_Olson_Landmap
    USE Precision_Mod

!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ESMF_GridComp), INTENT(INOUT), TARGET :: GC     ! Ref to this GridComp
    TYPE(ESMF_State),    INTENT(INOUT), TARGET :: Import ! Import State
    TYPE(ESMF_State),    INTENT(INOUT), TARGET :: Export ! Export State
    TYPE(ESMF_Clock),    INTENT(INOUT)         :: Clock  ! ESMF Clock object
    INTEGER,             INTENT(IN   )         :: Phase  ! Run phase (-1/1/2)
!
! !OUTPUT PARAMETERS:
!
    INTEGER,             INTENT(  OUT)         :: RC     ! Error return code
!
! !REMARKS:
!  We call routine Extract_ to return various values (i.e. grid parameters,
!  start & end dates, PET information, etc.) from the ESMF/MAPL environment.
!  We then pass those to GEOS-Chem via routine GCHP_CHUNK_RUN, which is
!  located in GEOS-Chem module ./GEOS-Chem/ESMF/gchp_chunk_mod.F90.

! !REVISION HISTORY:
!  06 Dec 2009 - A. da Silva - Initial version
!  See https://github.com/geoschem/geos-chem for history
!EOP
!------------------------------------------------------------------------------
!BOC

!
! LOCAL VARIABLES:
!
    ! Objects
    TYPE(ESMF_Grid)              :: Grid          ! ESMF Grid object
    TYPE(ESMF_Config)            :: MaplCF        ! Config (MAPL.rc)
    TYPE(ESMF_Config)            :: GeosCF        ! Config (GEOSCHEM*.rc)
    TYPE(ESMF_Alarm)             :: ALARM
    TYPE(ESMF_VM)                :: VM            ! ESMF VM object
    TYPE(ESMF_STATE)             :: INTSTATE

    ! Scalars
    LOGICAL                      :: am_I_Root     ! Are we on the root PET?
    LOGICAL                      :: IsChemTime    ! Chemistry alarm proxy
    LOGICAL                      :: IsRadTime     ! Radiation alarm proxy
    LOGICAL                      :: IsRunTime     ! Time to call GEOS-Chem
    LOGICAL                      :: IsTendTime    ! Time to calculate tendencies
    INTEGER                      :: IND           ! Species or tracer index
    INTEGER                      :: error         ! G-C error return code
    INTEGER(ESMF_KIND_I8)        :: advCount      ! # of clock advances
    INTEGER                      :: nymd          ! YYYY/MM/DD date
    INTEGER                      :: nhms          ! hh:mm:ss time
    INTEGER                      :: myPet         ! PET # we are on now
    INTEGER                      :: nPets         ! Total # of PETs
    INTEGER                      :: I, J, L       ! Loop indices
    INTEGER                      :: IM,JM,LM      ! Grid dimensions
    INTEGER                      :: LR, N         ! Loop indices
    INTEGER                      :: N_TRC         ! Shadow var: # of tracers
    INTEGER                      :: year          ! Current year
    INTEGER                      :: month         ! Current month
    INTEGER                      :: day           ! Current day
    INTEGER                      :: dayOfYr       ! Current day of year
    INTEGER                      :: hour          ! Current hour
    INTEGER                      :: minute        ! Current minute
    INTEGER                      :: second        ! Current second
    REAL                         :: UTC           ! Universal time
    REAL                         :: tsChem        ! Chem timestep [sec]
    REAL                         :: tsRad         ! RRTMG timestep [sec]
    REAL                         :: tsDyn         ! Dynamic timestep [sec]
    REAL                         :: hElapsed      ! Elapsed time [hours]
    REAL*8                       :: lonDeg        ! Longitude [degrees]
    REAL*8                       :: latDeg        ! Latitude [degrees]
    REAL*8                       :: P1, P2        ! Pressure variables
    CHARACTER(LEN=4)             :: petStr        ! String for PET #
    CHARACTER(LEN=ESMF_MAXSTR)   :: compName      ! Gridded Component name

    ! Allocatable local arrays
    REAL,  ALLOCATABLE, TARGET   :: zenith(:,:)   ! Solar zenith angle
    REAL,  ALLOCATABLE, TARGET   :: solar(:,:)    ! Solar insolation

    ! Pointer arrays needed to initialize from imports
    CHARACTER(LEN=2)             :: intStr
    REAL, POINTER                :: Ptr2d   (:,:)   => NULL()
    REAL, POINTER                :: Ptr3d   (:,:,:) => NULL()
    REAL(ESMF_KIND_R8), POINTER  :: Ptr2d_R8(:,:)   => NULL()
    REAL(ESMF_KIND_R8), POINTER  :: Ptr3d_R8(:,:,:) => NULL()

    ! Other pointer arrays
    REAL(ESMF_KIND_R4),  POINTER :: lonCtr  (:,:) ! Lon centers, this PET [rad]
    REAL(ESMF_KIND_R4),  POINTER :: latCtr  (:,:) ! Lat centers, this PET [rad
    TYPE(MAPL_MetaComp), POINTER :: STATE

    ! For CTM Mode
    ! ckeller, 8/22/19: In GEOS, PLE and AIRDENS are from the IMPORT state
    REAL(ESMF_KIND_R8),  POINTER :: PLE(:,:,:)     => NULL() ! INTERNAL: PEDGE

    ! RRTMG FDH needs to be able to read in dynamical heating

    ! Initialize variables used for reading Olson and MODIS LAI imports
    INTEGER            :: TT, VV, landTypeInt
    CHARACTER(len=64)  :: landTypeStr, varName, importName

    ! GCHP only local variables
    INTEGER                      :: trcID, RST
    REAL                         :: COEFF
    CHARACTER(LEN=ESMF_MAXSTR)   :: trcNAME,hcoNAME
    TYPE(ESMF_Field      )       :: trcFIELD
    TYPE(ESMF_FieldBundle)       :: trcBUNDLE
    REAL              , POINTER  :: fPtrArray(:,:,:)
    REAL(ESMF_KIND_R8), POINTER  :: fPtrVal, fPtr1D(:)
    INTEGER                      :: z_lb, z_ub

    ! Alarms
    type(GC_run_alarms), pointer :: GC_alarms
    type(GCRA_wrap)               :: GC_alarm_wrapper

    ! First call?
    LOGICAL, SAVE                :: FIRST = .TRUE.
    INTEGER                      :: NFD, K
    LOGICAL                      :: LAST
    TYPE(ESMF_Time        )      :: currTime, stopTime
    TYPE(ESMF_TimeInterval)      :: tsChemInt
    CHARACTER(len=ESMF_MAXSTR)   :: timestring1, timestring2
#ifdef JACOBIAN
    INTEGER                      :: primarySpcId
    CHARACTER(len=ESMF_MAXSTR)   :: primarySpcName
#endif

    !=======================================================================
    ! Run starts here
    !=======================================================================

    ! Are we on the root PET?
    am_I_Root = MAPL_Am_I_Root()

    ! Set up traceback info
    CALL ESMF_GridCompGet( GC, name=compName, __RC__ )

    ! Identify this routine to MAPL
    Iam = TRIM(compName)//'::Run_'

    ! Get the VM for optional memory prints (level >= 2)
    !-----------------------------------
    if ( MemDebugLevel > 0 ) THEN
       call ESMF_VmGetCurrent(VM, RC=STATUS)
       _VERIFY(STATUS)
    endif

    ! Get my MAPL_Generic state
    ! -------------------------
    CALL MAPL_GetObjectFromGC(GC, STATE, __RC__)

    ! Query the chemistry alarm.
    ! This checks if it's time to do chemistry, based on the time step
    ! set in AGCM.rc (GEOSCHEMCHEM_DT:). If the GEOS-Chem time step is not
    ! specified in AGCM.rc, the heartbeat will be taken (set in MAPL.rc).
    ! ----------------------------------------------------------------------
    CALL MAPL_Get(STATE, RUNALARM=ALARM, __RC__)
    IsChemTime = ESMF_AlarmIsRinging(ALARM, __RC__)

    ! if (am_I_Root) WRITE(*,*) ' Chem clock is reverse? ', ESMF_ClockIsReverse(CLOCK)
    ! Turn off alarm: only if it was on and this is phase 2 (don't turn off
    ! after phase 1 since this would prevent phase 2 from being executed).
    IF ( IsChemTime .AND. PHASE /= 1 .and. .not. ESMF_ClockIsReverse(CLOCK)) THEN
       CALL ESMF_AlarmRingerOff(ALARM, __RC__ )
    ENDIF

    ! Retrieve GEOS-Chem's internal alarms
    call ESMF_UserCompGetInternalState(GC,'gcchem_internal_alarms',GC_alarm_wrapper,status)
    _ASSERT(rc==0,'Could not retrieve radiation alarm')
    GC_alarms => GC_alarm_wrapper%ptr
    ! Query the radiation alarm
    IsRadTime = ESMF_AlarmIsRinging(GC_alarms%RRTMG_alarm,__RC__)

    ! Turn off alarm: only if it was on, chemistry will run, and this is phase 2
    If ( IsRadTime .and. IsChemTime .and. PHASE /= 1 ) Then
       CALL ESMF_AlarmRingerOff(GC_alarms%RRTMG_alarm, __RC__ )
    End If

    ! Get Internal state
    CALL MAPL_Get ( STATE, INTERNAL_ESMF_STATE=INTSTATE, __RC__ )

    ! ----------------------------------------------------------------------
    ! Check if we need to call the GEOS-Chem driver. The GEOS-Chem driver
    ! contains routines for the following processes:
    !
    ! Phase 1:
    ! (1) Convection:     --> Dynamics time step  (optional)
    ! (2) Dry deposition  --> Chemistry time step
    ! (3) Emissions       --> Chemistry time step
    !
    ! Phase 2:
    ! (4) Turbulence      --> Dynamics time step  (optional)
    ! (5) Chemistry       --> Chemistry time step
    ! (6) Wet deposition  --> Dynamics time step
    !
    ! Phase -1:
    ! Includes all of the above
    !
    ! Convection and turbulence are only called if the corresponding
    ! switches are turned on in the GEOS-Chem input file (geoschem_config.yml).
    !
    ! To avoid unnecessary calls to the GEOS-Chem driver routine, we
    ! check here if it's time for any of the processes listed above.
    ! The IsChemTime variable will be passed down to the GEOS-Chem driver
    ! to ensure that chemistry is only executed if it's time to do so.
    !
    ! The O3 and H2O tendencies will only be calculated when doing chemistry
    ! (set to zero otherwise). All other export variables are updated every
    ! time GEOS-Chem is called.
    ! ----------------------------------------------------------------------
    IsRunTime = IsChemTime
    IF ( Input_Opt%LCONV .AND. Phase /= 2 ) IsRunTime = .TRUE.
    IF ( Input_Opt%LTURB .AND. Phase /= 1 ) IsRunTime = .TRUE.
    IF ( Input_Opt%LWETD .AND. Phase /= 1 ) IsRunTime = .TRUE.

    !=======================================================================
    ! Skip GCC during replay, predictor step (posturm and cakelle2)
    !=======================================================================
    IF ( SkipReplayGCC ) THEN
       CALL ESMF_ClockGetAlarm(CLOCK, "PredictorActive", PredictorAlarm, RC=STATUS)
       VERIFY_(STATUS)

       PredictorActive = ESMF_AlarmIsRinging( PredictorAlarm, RC=STATUS )
       VERIFY_(STATUS)

       IF ( PredictorActive ) THEN
          IsRunTime = .FALSE.
          IF ( am_I_root ) write(*,*) '  --- Skipping GCC during Predictor Step '
       END IF
    END IF

    ! Is it time to update tendencies?
    ! Tendencies shall only be updated when chemistry is done, which is
    ! Phase -1 or 2.
    IsTendTime = ( IsChemTime .AND. Phase /= 1 )

    ! Start timers
    ! ------------
    CALL MAPL_TimerOn(STATE, "TOTAL")

    ! Get pointers to fields in import, internal, and export states defined
    ! in the registry file. This has to be done on the first call only.
    IF ( FIRST ) THEN
#      include "GCHPchem_GetPointer___.h"

       !IF ( IsCTM ) THEN
       call MAPL_GetPointer ( IMPORT, PLE,      'PLE',     __RC__ )
       !ENDIF

       ! Pass IMPORT/EXPORT object to HEMCO state object
       !CALL GetHcoState( HcoState )
       _ASSERT(ASSOCIATED(HcoState),'HcoState is not associated')
       HcoState%GRIDCOMP => GC
       HcoState%importState   => IMPORT
       HcoState%exportState   => EXPORT
       !HcoState => NULL()

    ENDIF


    ! Run when it's time to do so
    ! Always run on first call to make sure that all variables become
    ! properly specified and initialized.
    ! ------------------------------------------------------------------
    RunningGEOSChem: IF(IsRunTime .OR. FIRST) THEN

       CALL MAPL_TimerOn(STATE, "RUN"  )

       ! Get various parameters from the ESMF/MAPL framework
       CALL Extract_( GC,                   &  ! Ref to this Gridded Component
                      Clock,                &  ! ESMF Clock object
                      Grid      = Grid,     &  ! ESMF Grid object
                      MaplCf    = MaplCF,   &  ! ESMF Config obj (MAPL*.rc)
                      GeosCf    = GeosCF,   &  ! ESMF Config obj (GEOSCHEM*.rc)
                      tsChem    = tsChem,   &  ! Chemistry timestep [sec]
                      tsRad     = tsRad,    &  ! Radiation timestep [sec]
                      tsDyn     = tsDyn,    &  ! Dynamic timestep [sec]
                      nymd      = nymd,     &  ! Current YYYY/MM/DD date
                      nhms      = nhms,     &  ! Current hh:mm:ss time
                      year      = year,     &  ! Current year
                      month     = month,    &  ! Current month
                      day       = day,      &  ! Current day
                      dayOfYr   = dayOfYr,  &  ! Current day of year
                      hour      = hour,     &  ! Current hour
                      minute    = minute,   &  ! Current minute
                      helapsed  = hElapsed, &  ! Elapsed hours
                      advCount  = advCount, &  ! # of times clock has advanced
                      utc       = utc,      &  ! Universal time [hours]
                      localpet  = myPet,    &  ! # of the PET we are on now
                      petCount  = nPets,    &  ! Total # of PETs
                      __RC__ )

       ! For convenience, set grid dimension variables. These are being
       ! used in Includes_Before_Run.H (ckeller, 8/22/19)
       IM = State_Grid%NX
       JM = State_Grid%NY
       LM = State_Grid%NZ

       ! Allocate GMAO_ZTH (declared at top of module)
       IF ( .not. ALLOCATED( zenith ) ) THEN
          ALLOCATE( zenith(State_Grid%NX,State_Grid%NY), STAT=STATUS)
          _VERIFY(STATUS)
       ENDIF

       ! Allocate GMAO_SLR (declared @ top of module)
       IF ( .not. ALLOCATED( solar ) ) THEN
          ALLOCATE( solar(State_Grid%NX,State_Grid%NY), STAT=STATUS)
          _VERIFY(STATUS)
       ENDIF

       ! Call EXTRACT a second time to get the solar zenith
       ! angle and solar insolation fields
       CALL Extract_( GC,                   &  ! Ref to this Gridded Component
                      Clock,                &  ! ESMF Clock object
                      Grid      = Grid,     &  ! ESMF Grid object
                      MaplCf    = MaplCF,   &  ! ESMF Config obj (MAPL*.rc)
                      GeosCf    = GeosCF,   &  ! ESMF Config obj (GEOSCHEM*.rc)
                      lonCtr    = lonCtr,   &  ! Lon centers on this PET [rad]
                      latCtr    = latCtr,   &  ! Lat centers on this PET [rad]
                      ZTH       = zenith,   &  ! Solar zenith angle
                      SLR       = solar,    &  ! Solar insolation
                      __RC__ )

       ! MSL - shift from 0 - 360 to -180 - 180 degree grid
       where (lonCtr .gt. MAPL_PI ) lonCtr = lonCtr - 2*MAPL_PI

       ! Pass grid area [m2] obtained from dynamics component to State_Grid
       CALL MAPL_GetPointer( IMPORT, Ptr2d_R8, 'AREA', __RC__ )
       State_Grid%Area_M2 = Ptr2d_R8
       Ptr2d_R8 => NULL()

       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       ! KLUDGE (mps, 5/23/19):
       ! Copy to State_Met%AREA_M2 to avoid breaking GCHP benchmarks, which
       ! require the AREA_M2 field saved out to the StateMet diagnostic
       ! collection for things like computing emission totals.
       !
       State_Met%Area_M2 = State_Grid%Area_M2
       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

       !=======================================================================
       ! Prevent use of (occasional) MAPL_UNDEF tropopause pressures
       !=======================================================================

       ! GCCTROPP contains the last valid tropopause pressure
       WHERE ( TROPP /= MAPL_UNDEF ) GCCTROPP = TROPP

       ! If any values in GCCTROPP are undefined, stop the run
       IF ( ANY( GCCTROPP == MAPL_UNDEF ) ) THEN
          PRINT *,TRIM(Iam)//": At least one invalid tropopause pressure."
          STATUS = GC_FAILURE
          _VERIFY(STATUS)
       ENDIF

!       !=======================================================================
!       ! pre-Run method array assignments. This passes the tracer arrays from
!       ! the internal state to State_Chm. On the first call, it also fills the
!       ! internal species arrays in State_Chm with the values read from the
!       ! restart file (and stored in the internal state).
!       !=======================================================================
       CALL MAPL_TimerOn(STATE, "CP_BFRE")
#include "Includes_Before_Run.H"
       CALL MAPL_TimerOff(STATE, "CP_BFRE")

       !=======================================================================
       ! Point GEOS-Chem species concentration arrays to internal state
       !=======================================================================
       DO I = 1, SIZE(Int2Spc,1)
          IF ( Int2Spc(I)%ID <= 0 ) CYCLE
          State_Chm%Species(Int2Spc(I)%ID)%Conc => Int2Spc(I)%Internal(:,:,State_Grid%NZ:1:-1)
       ENDDO


      !=======================================================================
      ! On first call, populate State_Chm%Species(N)%Conc with background
      ! values if the species is missing from the restart file and missing
      ! species are allowed.
      !=======================================================================
       IF ( FIRST ) THEN

          ! Get Generic State
          call MAPL_GetObjectFromGC ( GC, STATE, RC=STATUS)
          _VERIFY(STATUS)
          ! Get Internal state
          CALL MAPL_Get ( STATE, INTERNAL_ESMF_STATE=INTERNAL, __RC__ )

          ! Loop over all species and get info from spc db
          DO N = 1, State_Chm%nSpecies
             ThisSpc => State_Chm%SpcData(N)%Info
             IF ( TRIM(ThisSpc%Name) == '' ) CYCLE
             IND = IND_( TRIM(ThisSpc%Name ) )
             IF ( IND < 0 ) CYCLE

             ! Determine if species in restart file
             CALL ESMF_StateGet( INTERNAL, TRIM(SPFX) // TRIM(ThisSpc%Name),  &
                  trcFIELD, RC=RC )
             CALL ESMF_AttributeGet( trcFIELD, NAME="RESTART",                &
                  VALUE=RST, RC=STATUS )

             ! Set spc conc to background value if rst skipped or var not there
             IF ( ( RC  /= ESMF_SUCCESS           .OR.     &
                    RST == MAPL_RestartBootstrap  .OR.     &
                    RST == MAPL_RestartSkipInitial  )      &
#ifdef JACOBIAN
                    .AND. .NOT. ThisSpc%Is_JacobianTracer  &
#endif
                    ) THEN
                DO L = 1, State_Grid%NZ
                DO J = 1, State_Grid%NY
                DO I = 1, State_Grid%NX
                   IF ( L > State_Met%MaxChemLev .AND. &
                            ( .NOT. ThisSpc%Is_Advected ) ) THEN
                      ! For non-advected spc at L > MaxChemLev, use small number
                      State_Chm%Species(IND)%Conc(I,J,L) = 1.0E-30_FP
                   ELSE
                      ! For all other cases, use the background value in spc db
                      State_Chm%Species(IND)%Conc(I,J,L) = ThisSpc%BackgroundVV
                   ENDIF
                ENDDO
                ENDDO
                ENDDO
                IF ( MAPL_am_I_Root()) THEN
                   WRITE(*,*)  &
                        '   WARNING: using background values from species database'&
                        //' for species '//trim(ThisSpc%Name) 
                ENDIF
             ENDIF

#ifdef JACOBIAN
             ! Do special handling if this is a Jacobian tracer             
             IF ( ThisSpc%Is_JacobianTracer ) THEN
                primarySpcName = ThisSpc%Name(1:LEN(trim(ThisSpc%Name))-8)
                primarySpcId = IND_(trim(primarySpcName))
                State_Chm%Species(IND)%Conc = State_Chm%Species(primarySpcId)%Conc
                IF ( MAPL_am_I_Root()) THEN
                   WRITE(*,*) '   INFO: using the initial concentration of ' &
                        // trim(primarySpcName) //' for the Jacobian tracer ' &
                        // trim(ThisSpc%Name)
                ENDIF
             ENDIF
#endif

             ThisSpc => NULL()
          ENDDO
       ENDIF

       !=======================================================================
       ! On first call, initialize certain State_Chm arrays from
       ! internal state (restart file) if they are found. Do not initialize
       ! from met-fields in restart file that were added for post-processing.
       ! Do not add delta pressure since used for mass conservation scaling
       ! in FV3 prior advection which comes before GEOS-Chem.
       !=======================================================================
       IF ( FIRST ) THEN
          CALL MAPL_GetPointer( INTSTATE, Ptr3d_R8, 'H2O2AfterChem', notFoundOK=.TRUE., __RC__ )
          IF ( ASSOCIATED(Ptr3d_R8) .AND. ASSOCIATED(State_Chm%H2O2AfterChem) ) THEN
             State_Chm%H2O2AfterChem = Ptr3d_R8(:,:,State_Grid%NZ:1:-1)
          ENDIF
          Ptr3d_R8 => NULL()

          CALL MAPL_GetPointer( INTSTATE, Ptr3d_R8, 'SO2AfterChem', notFoundOK=.TRUE., __RC__ )
          IF ( ASSOCIATED(Ptr3d_R8) .AND. ASSOCIATED(State_Chm%SO2AfterChem) ) THEN
             State_Chm%SO2AfterChem = Ptr3d_R8(:,:,State_Grid%NZ:1:-1)
          ENDIF
          Ptr3d_R8 => NULL()

          CALL MAPL_GetPointer( INTSTATE, Ptr2d_R8, 'DryDepNitrogen', notFoundOK=.TRUE., __RC__ )
          IF ( ASSOCIATED(Ptr2d_R8) .AND. ASSOCIATED(State_Chm%DryDepNitrogen) ) THEN
             State_Chm%DryDepNitrogen = Ptr2d_R8
          ENDIF
          Ptr2d_R8 => NULL()

          CALL MAPL_GetPointer( INTSTATE, Ptr2d_R8, 'WetDepNitrogen', notFoundOK=.TRUE., __RC__ )
          IF ( ASSOCIATED(Ptr2d_R8) .AND. ASSOCIATED(State_Chm%WetDepNitrogen) ) THEN
             State_Chm%WetDepNitrogen = Ptr2d_R8
          ENDIF
          Ptr2d_R8 => NULL()

          CALL MAPL_GetPointer( INTSTATE, Ptr3d_R8, 'KPPHvalue', notFoundOK=.TRUE., __RC__ )
          IF ( ASSOCIATED(Ptr3d_R8) .AND. ASSOCIATED(State_Chm%KPPHvalue) ) THEN
             State_Chm%KPPHvalue(:,:,1:State_Met%MaxChemLev) =       &
            Ptr3d_R8(:,:,State_Grid%NZ:State_Grid%NZ-State_Met%MaxChemLev+1:-1)
          ENDIF
          Ptr3d_R8 => NULL()

          CALL MAPL_GetPointer( INTSTATE, Ptr3D, 'STATE_PSC', notFoundOK=.TRUE., __RC__ )
          IF ( ASSOCIATED(Ptr3D) .AND. ASSOCIATED(State_Chm%State_PSC) ) THEN
             State_Chm%State_PSC(:,:,:) = Ptr3D(:,:,LM:1:-1)
          ENDIF
          Ptr3D => NULL()

          CALL MAPL_GetPointer( INTSTATE, Ptr3d_R8, 'AeroH2O_SNA', notFoundOK=.TRUE., __RC__ )
          IF ( ASSOCIATED(Ptr3d_R8) .AND. ASSOCIATED(State_Chm%AeroH2O) ) THEN
             State_Chm%AeroH2O(:,:,1:State_Grid%NZ,NDUST+1) =       &
                                  Ptr3d_R8(:,:,State_Grid%NZ:1:-1)
          ENDIF
          Ptr3d_R8 => NULL()

          CALL MAPL_GetPointer( INTSTATE, Ptr3d_R8, 'ORVCSESQ', notFoundOK=.TRUE., __RC__ )
          IF ( ASSOCIATED(Ptr3d_R8) .AND. ASSOCIATED(State_Chm%ORVCsesq) ) THEN
             State_Chm%ORVCsesq(:,:,1:State_Grid%NZ) =       &
                                  Ptr3d_R8(:,:,State_Grid%NZ:1:-1)
          ENDIF
          Ptr3d_R8 => NULL()

          CALL MAPL_GetPointer( INTSTATE, Ptr2D_R8, 'JOH', notFoundOK=.TRUE., __RC__ )
          IF ( ASSOCIATED(Ptr2D_R8) .AND. ASSOCIATED(State_Chm%JOH) ) THEN
             State_Chm%JOH(:,:) = Ptr2D_R8(:,:)
          ENDIF
          Ptr2D_R8 => NULL()

          CALL MAPL_GetPointer( INTSTATE, Ptr2D_R8, 'JNO2', notFoundOK=.TRUE., __RC__ )
          IF ( ASSOCIATED(Ptr2D_R8) .AND. ASSOCIATED(State_Chm%JNO2) ) THEN
             State_Chm%JNO2(:,:) = Ptr2D_R8(:,:)
          ENDIF
          Ptr2D_R8 => NULL()

          CALL MAPL_GetPointer( INTSTATE, Ptr3d_R8, 'TSTRAT_ADJ', notFoundOK=.TRUE., __RC__ )
          IF ( ASSOCIATED(Ptr3d_R8) .AND. ASSOCIATED(State_Chm%TStrat_Adj) ) THEN
             State_Chm%TStrat_Adj(:,:,1:State_Grid%NZ) =       &
                                  Ptr3d_R8(:,:,State_Grid%NZ:1:-1)
          ENDIF
          Ptr3d_R8 => NULL()
       ENDIF


       !=======================================================================
       ! Set Olson land map types from import of Olson file.
       !=======================================================================
       ! We are currently using land type fractions derived from the 2001
       ! Olson land map instead of GEOS5 vegetation type fractions. Fractions
       ! are calculated by ExtData using conservate fractional regridding of
       ! the native 0.25x0.25 resolution file. (ewl, 11/29/16)
       !
       ! Previous:
       ! Set land types in State_Met from GEOS5 vegetation type fractions or
       ! OLSON land type fractions. For now, the land types are treated as
       ! static and obtained from offline fields. The routine below thus needs
       ! to be called only once.
       ! Once the GEOS-5 land types are dynamic, we should import those from
       ! the surface component (field ITY, or better: vegetation type fractions
       ! per grid box).                                   (ckeller, 01/06/2015)
       !
       !=======================================================================
       IF ( FIRST ) THEN

          ! Set Olson fractional land type from import (ewl)
          If (am_I_Root) Write(*,'(a)') 'Initializing land type ' // &
                           'fractions from Olson imports'
          Ptr2d => NULL()
          DO TT = 1, NSURFTYPE

             ! Create two-char string for land type
             landTypeInt = TT-1
             WRITE ( landTypeStr, '(I2.2)' ) landTypeInt
             importName = 'OLSON' // TRIM(landTypeStr)

             ! Get pointer and populate State_Met variable
             CALL MAPL_GetPointer ( IMPORT, Ptr2D, TRIM(importName),  &
                                    notFoundOK=.TRUE., __RC__ )
             If ( Associated(Ptr2D) ) Then
                State_Met%LandTypeFrac(:,:,TT) = Ptr2D(:,:)
             ELSE
                WRITE(6,*) TRIM(importName) // ' pointer is not associated'
             ENDIF
             Ptr2D => NULL()
          ENDDO

          ! Compute State_Met variables IREG, ILAND, IUSE, and FRCLND
          CALL Compute_Olson_Landmap( Input_Opt, State_Grid, State_Met, RC )
          _ASSERT(RC==GC_SUCCESS,'Error calling Compute_Olson_Landmap')
       ENDIF

       !=======================================================================
       ! Execute GEOS-Chem on multiple PETs
       !=======================================================================

       ! Fix negatives!
       ! These can be brought in as an artifact of convection.
       DO N = 1, State_Chm%nSpecies
          WHERE ( State_Chm%Species(N)%Conc < 0.0e0 )
             State_Chm%Species(N)%Conc = 1.0e-36
          END WHERE 
       ENDDO

       ! Execute GEOS-Chem if it's time to run it
       IF ( IsRunTime ) THEN

          ! This is mostly for testing
          IF ( FIRST ) THEN
             IF ( am_I_Root ) THEN
                WRITE(*,*) ''
                WRITE(*,*) 'Doing warm GEOS-Chem restart'
                WRITE(*,*) ''
             ENDIF
          ENDIF

             ! Optional memory prints (level >= 2)
             if ( MemDebugLevel > 0 ) THEN
                call ESMF_VMBarrier(vm, RC=STATUS)
                _VERIFY(STATUS)
                call MAPL_MemUtilsWrite(VM, &
                  'Chem_GridCompMod, before chunk_run', RC=STATUS )
                _VERIFY(STATUS)
             endif

             CALL MAPL_TimerOn(STATE, "DO_CHEM")

             ! NOTE: Second was not extracted previously; set to 0 for now
             second = 0


             ! Run the GEOS-Chem column chemistry code for the given phase
             CALL GCHP_Chunk_Run( GC         = GC,         & ! Grid comp ref.
                                  nymd       = nymd,       & ! Current YYYYMMDD
                                  nhms       = nhms,       & ! Current hhmmss
                                  year       = year,       & ! Current year
                                  month      = month,      & ! Current month
                                  day        = day,        & ! Current day
                                  dayOfYr    = dayOfYr,    & ! Current doy
                                  hour       = hour,       & ! Current hour
                                  minute     = minute,     & ! Current minute
                                  second     = second,     & ! Current second
                                  utc        = utc,        & ! Current UTC [hrs]
                                  hElapsed   = hElapsed,   & ! Elapsed hours
                                  Input_Opt  = Input_Opt,  & ! Input Options
                                  State_Chm  = State_Chm,  & ! Chemistry State
                                  State_Diag = State_Diag, & ! Diagnostics State
                                  State_Grid = State_Grid, & ! Grid State
                                  State_Met  = State_Met,  & ! Meteorology State
                                  Phase      = Phase,      & ! Run phase
                                  IsChemTime = IsChemTime, & ! Time for chem?
                                  IsRadTime  = IsRadTime,  & ! Time for RRTMG?
                                  __RC__                  )  ! Success or fail?

             CALL MAPL_TimerOff(STATE, "DO_CHEM")

             ! Optional memory prints (level >= 2)
             if ( MemDebugLevel > 0 ) THEN
                call ESMF_VMBarrier(vm, RC=STATUS)
                _VERIFY(STATUS)
                call MAPL_MemUtilsWrite(VM, &
                  'Chem_GridCompMod, after  chunk_run', RC=STATUS )
                _VERIFY(STATUS)
             endif
             where( State_Met%HFLUX .eq. 0.) State_Met%HFLUX = 1e-5

       ENDIF !IsRunTime

       !=======================================================================
       ! post-Run method array assignments. This copies the values back from
       ! the State_Chm tracer arrays to the internal state, so that they can
       ! be seen by other components (moist, turbulence, ...)
       !=======================================================================

       CALL MAPL_TimerOn(STATE, "CP_AFTR")


       CALL MAPL_TimerOff(STATE, "CP_AFTR")

       ! Update non-species dynamic internal state arrays post-run
       ! every timestep, except for area which can be set first run only.
       ! Some of the fields are for post-processing but need to be updated
       ! mid-run for inclusion in mid-run checkpoint files.
       CALL MAPL_GetPointer( INTSTATE, Ptr2d_R8, 'DryDepNitrogen', &
                             notFoundOK=.TRUE., __RC__ )
       IF (ASSOCIATED(Ptr2d_R8) .AND. ASSOCIATED(State_Chm%DryDepNitrogen)) THEN
          Ptr2d_R8 = State_Chm%DryDepNitrogen
       ENDIF
       Ptr2d_R8 => NULL()
       CALL MAPL_GetPointer( INTSTATE, Ptr2d_R8, 'WetDepNitrogen', &
                             notFoundOK=.TRUE., __RC__ )
       IF (ASSOCIATED(Ptr2d_R8) .AND. ASSOCIATED(State_Chm%WetDepNitrogen)) THEN
          Ptr2d_R8 = State_Chm%WetDepNitrogen
       ENDIF
       Ptr2d_R8 => NULL()

       CALL MAPL_GetPointer( INTSTATE, Ptr3d_R8, 'H2O2AfterChem', &
                             notFoundOK=.TRUE., __RC__ )
       IF (ASSOCIATED(Ptr3d_R8) .AND. ASSOCIATED(State_Chm%H2O2AfterChem)) THEN
          Ptr3d_R8(:,:,State_Grid%NZ:1:-1) = State_Chm%H2O2AfterChem
       ENDIF
       Ptr3d_R8 => NULL()

       CALL MAPL_GetPointer( INTSTATE, Ptr3d_R8, 'SO2AfterChem', &
                             notFoundOK=.TRUE., __RC__ )
       IF (ASSOCIATED(Ptr3d_R8) .AND. ASSOCIATED(State_Chm%SO2AfterChem)) THEN
          Ptr3d_R8(:,:,State_Grid%NZ:1:-1) = State_Chm%SO2AfterChem
       ENDIF
       Ptr3d_R8 => NULL()

       CALL MAPL_GetPointer( INTSTATE, Ptr3d_R8, 'KPPHvalue', &
                             notFoundOK=.TRUE., __RC__ )
       IF (ASSOCIATED(Ptr3d_R8) .AND. ASSOCIATED(State_Chm%KPPHvalue)) THEN
          Ptr3d_R8(:,:,1:State_Grid%NZ-State_Met%MaxChemLev) = 0.0
          Ptr3d_R8(:,:,State_Grid%NZ:State_Grid%NZ-State_Met%MaxChemLev+1:-1)=&
             State_Chm%KPPHvalue(:,:,1:State_Met%MaxChemLev)
       ENDIF
       Ptr3d_R8 => NULL()

       CALL MAPL_GetPointer( INTSTATE, Ptr3d_R8, 'AeroH2O_SNA', &
                             notFoundOK=.TRUE., __RC__ )
       IF (ASSOCIATED(Ptr3d_R8) .AND. ASSOCIATED(State_Chm%AeroH2O)) THEN
          Ptr3d_R8(:,:,State_Grid%NZ:1:-1) =  &
                    State_Chm%AeroH2O(:,:,1:State_Grid%NZ,NDUST+1)
       ENDIF
       Ptr3d_R8 => NULL()

       CALL MAPL_GetPointer( INTSTATE, Ptr3d_R8, 'ORVCSESQ', &
                             notFoundOK=.TRUE., __RC__ )
       IF (ASSOCIATED(Ptr3d_R8) .AND. ASSOCIATED(State_Chm%ORVCsesq)) THEN
          Ptr3d_R8(:,:,State_Grid%NZ:1:-1) =  &
                    State_Chm%ORVCsesq(:,:,1:State_Grid%NZ)
       ENDIF
       Ptr3d_R8 => NULL()

       CALL MAPL_GetPointer( INTSTATE, Ptr2D_R8, 'JOH', &
                             notFoundOK=.TRUE., __RC__ )
       IF (ASSOCIATED(Ptr2D_R8) .AND. ASSOCIATED(State_Chm%JOH)) THEN
          Ptr2d_R8(:,:) = State_Chm%JOH(:,:)
       ENDIF
       Ptr2D_R8 => NULL()

       CALL MAPL_GetPointer( INTSTATE, Ptr2D_R8, 'JNO2', &
                             notFoundOK=.TRUE., __RC__ )
       IF (ASSOCIATED(Ptr2D_R8) .AND. ASSOCIATED(State_Chm%JNO2)) THEN
          Ptr2d_R8(:,:) = State_Chm%JNO2(:,:)
       ENDIF
       Ptr2D_R8 => NULL()

       CALL MAPL_GetPointer( INTSTATE, Ptr3D, 'STATE_PSC', &
                             notFoundOK=.TRUE., __RC__ )
       IF (ASSOCIATED(Ptr3D) .AND. ASSOCIATED(State_Chm%State_PSC)) THEN
          Ptr3d(:,:,LM:1:-1) = State_Chm%State_PSC(:,:,:)
       ENDIF
       Ptr3D => NULL()

       CALL MAPL_GetPointer( INTSTATE, Ptr3d_R8, 'DELP_DRY', &
                             notFoundOK=.TRUE., __RC__ )
       IF (ASSOCIATED(Ptr3d_R8) .AND. ASSOCIATED(State_Met%DELP_DRY)) THEN
          Ptr3d_R8(:,:,State_Grid%NZ:1:-1) =  &
                    State_Met%DELP_DRY(:,:,1:State_Grid%NZ)
       ENDIF
       Ptr3d_R8 => NULL()

       CALL MAPL_GetPointer( INTSTATE, Ptr3d_R8, 'TSTRAT_ADJ', &
                             notFoundOK=.TRUE., __RC__ )
       IF (ASSOCIATED(Ptr3d_R8) .AND. ASSOCIATED(State_Met%DELP_DRY) &
           .AND. ASSOCIATED(State_Chm%TStrat_Adj) ) THEN
          Ptr3d_R8(:,:,State_Grid%NZ:1:-1) =  &
                    State_Chm%TStrat_Adj(:,:,1:State_Grid%NZ)
       ENDIF
       Ptr3d_R8 => NULL()

       CALL MAPL_GetPointer( INTSTATE, Ptr3d_R8, 'BXHEIGHT' ,&
                             notFoundOK=.TRUE., __RC__ )
       IF (ASSOCIATED(Ptr3d_R8) .AND. ASSOCIATED(State_Met%BXHEIGHT)) THEN
          Ptr3d_R8(:,:,State_Grid%NZ:1:-1) =  &
                    State_Met%BXHEIGHT(:,:,1:State_Grid%NZ)
       ENDIF
       Ptr3d_R8 => NULL()

       CALL MAPL_GetPointer( INTSTATE, Ptr2d_R8, 'TropLev', &
                             notFoundOK=.TRUE., __RC__ )
       IF (ASSOCIATED(Ptr2d_R8) .AND. ASSOCIATED(State_Met%TropLev)) THEN
          Ptr2d_R8 = State_Met%TropLev
       ENDIF
       Ptr2d_R8 => NULL()

       ! Only update area the first timestep
       IF ( FIRST ) THEN
          CALL MAPL_GetPointer( INTSTATE, Ptr2d_R8, 'AREA', &
                                notFoundOK=.TRUE., __RC__ )
          IF ( ASSOCIATED(Ptr2d_R8) .AND. ASSOCIATED(State_Met%AREA_M2) ) THEN
             Ptr2d_R8 = State_Met%AREA_M2
          ENDIF
          Ptr2d_R8 => NULL()
       ENDIF
       
       ! Stop timer
       ! ----------
       CALL MAPL_TimerOff(STATE, "RUN"  )

    ENDIF RunningGEOSChem

    !=======================================================================
    ! Copy HISTORY.rc diagnostic data to exports. Includes HEMCO emissions
    ! diagnostics but excludes internal state and exports created explicitly
    ! in Chem_GridCompMod. NOTE: Exports created explicitly in Chem_GridCompMod
    ! will eventually be moved elsewhere as diagnostics for use with GEOS-5.
    ! (ewl, 11/2/17)
    !=======================================================================
    IF ( FIRST ) THEN
       CALL HistoryExports_SetDataPointers( am_I_Root,     EXPORT,    &
                                            HistoryConfig, State_Chm, &
                                            State_Diag,    State_Met, &
                                            STATUS )
       _VERIFY(STATUS)
    ENDIF
    CALL CopyGCStates2Exports( am_I_Root, Input_Opt, HistoryConfig, STATUS )
    _VERIFY(STATUS)

    !=======================================================================
    ! All done
    !=======================================================================

    IF ( ALLOCATED( zenith     ) ) DEALLOCATE( zenith     )
    IF ( ALLOCATED( solar      ) ) DEALLOCATE( solar      )

    ! Stop timer
    ! ----------
    CALL MAPL_TimerOff(STATE, "TOTAL")

    ! Update first flags
    FIRST = .FALSE.

    ! Successful return
    _RETURN(ESMF_SUCCESS)

    ! Formats
100 FORMAT( '---> DATE: ', i4.4, '/', i2.2, '/', i2.2,      &
            '  GMT: ', i2.2, ':', i2.2, '  X-HRS: ', f11.3 )
110 FORMAT( 'Box (',i3,',',i3,') on PET ', i3, ' has coords: ', 2f7.2, &
               ' LocT = ', f9.4 )
200 FORMAT( '### ',                                           / &
            '### ', a, '  |  Execution on PET # ',      i5.5, / &
            '###' )

  END SUBROUTINE Run_

  
  !------------------------------------------------------------------------------
  !                  GEOS-Chem Global Chemical Model                            !
  !------------------------------------------------------------------------------
  
  ! !DESCRIPTION: Finalize_ is the finalize method of the GEOSCHEMchem gridded
  !  component.

  SUBROUTINE Finalize_( GC, Import, Export, Clock, RC )

    USE CMN_Size_Mod,          ONLY : NDUST
    USE Input_Opt_Mod,         ONLY : OptInput
    USE Input_Opt_Mod,         ONLY : Cleanup_Input_Opt
    USE State_Chm_Mod,         ONLY : ChmState, Cleanup_State_Chm
    USE State_Diag_Mod,        ONLY : DgnState, Cleanup_State_Diag
    USE State_Grid_Mod,        ONLY : GrdState, Cleanup_State_Grid
    USE State_Met_Mod,         ONLY : MetState, Cleanup_State_Met
    USE HCO_Interface_GC_Mod,  ONLY : HCOI_GC_FINAL

    TYPE(ESMF_GridComp), INTENT(INOUT)         :: GC       ! Ref. to this GC
    TYPE(ESMF_State),    INTENT(INOUT)         :: Import   ! Import State
    TYPE(ESMF_State),    INTENT(INOUT)         :: Export   ! Export State
    TYPE(ESMF_Clock),    INTENT(INOUT)         :: Clock    ! ESMF Clock object

    INTEGER,             INTENT(OUT)           :: RC       ! Success or failure?

    ! !REMARKS:
    !  We call routine Extract_ to return various values (i.e. grid parameters,
    !  start & end dates, PET information, etc.) from the ESMF/MAPL environment.
    !  We then pass those to GEOS-Chem via routine GCHP_CHUNK_FINAL, which is
    !  located in GEOS-Chem module ./GEOS-Chem/ESMF/gchp_chunk_mod.F90.

    ! Objects
    TYPE(ESMF_Grid)            :: Grid        ! ESMF Grid object
    TYPE(ESMF_Config)          :: MaplCF      ! Config (MAPL.rc)
    TYPE(ESMF_Config)          :: GeosCF      ! Config (GEOSCHEM*.rc)

    type(GC_run_alarms), pointer :: GC_alarms
    type(GCRA_wrap)              :: GC_alarm_wrapper

    ! Scalars
    LOGICAL                    :: am_I_Root   ! Are we on the root PET?
    CHARACTER(LEN=ESMF_MAXSTR) :: compName    ! Gridded component name
    INTEGER                    :: error       ! GEOS-Chem error code
    INTEGER                    :: myPet       ! # of PET we are on now
    INTEGER                    :: I,  J,  L   ! Loop indices
    REAL                       :: UTC         ! UTC time [hours]
    
    ! ewl: added for new internal state vars
    CHARACTER(LEN=ESMF_MAXSTR) :: importName, intStr
    INTEGER                    :: LM

    ! Pointers
    TYPE(MAPL_MetaComp), POINTER :: STATE
    TYPE(Species),       POINTER :: ThisSpc

    ! For species copying
    INTEGER                     :: IND
    TYPE(ESMF_STATE)            :: INTSTATE
    REAL, POINTER               :: Ptr2D(:,:)      => NULL()
    REAL, POINTER               :: Ptr3D(:,:,:)    => NULL()
    REAL(ESMF_KIND_R8), POINTER :: Ptr2D_R8(:,:)   => NULL()
    REAL(ESMF_KIND_R8), POINTER :: Ptr3D_R8(:,:,:) => NULL()

    INTEGER                     :: N, K, NFD
    CHARACTER(LEN=ESMF_MAXSTR)  :: TrcName

    ! Are we on the root PET
    am_I_Root = MAPL_Am_I_Root()

    ! Set number of levels
    LM = State_Grid%NZ

    ! Set up traceback info
    CALL ESMF_GridCompGet( GC, name=compName, __RC__ )

    ! Identify this routine to MAPL
    Iam = TRIM(compName)//'::Finalize_'

    ! Get my MAPL_Generic state
    ! -------------------------
    CALL MAPL_GetObjectFromGC(GC, STATE, RC=STATUS)
    _VERIFY(STATUS)

    ! Start timers
    ! ------------
!    CALL MAPL_TimerOn(STATE, "TOTAL")
!    CALL MAPL_TimerOn(STATE, "FINALIZE")

    ! Get various parameters from the ESMF/MAPL framework
    CALL Extract_( GC,                 &    ! Ref to this Gridded Component
                   Clock,              &    ! ESMF Clock object
                   Grid     = Grid,    &    ! ESMF Grid object
                   MaplCF   = MaplCF,  &    ! ESMF Config obj (MAPL.rc)
                   GeosCF   = GeosCF,  &    ! ESMF Config obj (GEOSCHEM*.rc)
                   utc      = utc,     &    ! Universal time [hours]
                   localPET = myPet,   &    ! PET # we are on now
                   __RC__ )

    ! Destroy the internal alarms
    call ESMF_UserCompGetInternalState(GC,'gcchem_internal_alarms',GC_alarm_wrapper,status)
    _ASSERT(status==0,'Could not find GC alarms for destruction')

    GC_alarms => GC_alarm_wrapper%ptr
    call ESMF_AlarmDestroy(GC_alarms%RRTMG_alarm,rc=status)
    _ASSERT(status==0,'Could not destroy radiation alarm')
    deallocate(GC_alarms,stat=status)
    _ASSERT(status==0,'Could not deallocate GC alarms')

    ! Finalize HEMCO
    CALL HCOI_GC_FINAL( .FALSE., RC )
    IF ( Input_Opt%AmIRoot ) THEN
       IF ( RC == GC_SUCCESS ) THEN
          write(*,'(a)') 'HEMCO::Finalize... OK.'
       ELSE
          write(*,'(a)') 'HEMCO::Finalize... FAILURE.'
       ENDIF
    ENDIF

    ! Deallocate fields of the Chemistry State object
    CALL Cleanup_State_Chm( State_Chm, RC )
    IF ( Input_Opt%AmIRoot ) THEN
       IF ( RC == GC_SUCCESS ) THEN
          write(*,'(a)') 'Chem::State_Chm Finalize... OK.'
       ELSE
          write(*,'(a)') 'Chem::State_Chm Finalize... FAILURE.'
       ENDIF
    ENDIF

    ! Deallocate fields of the Diagnostics State object
    CALL Cleanup_State_Diag( State_Diag, RC )
    IF ( Input_Opt%AmIRoot ) THEN
       IF ( RC == GC_SUCCESS ) THEN
          write(*,'(a)') 'Chem::State_Diag Finalize... OK.'
       ELSE
          write(*,'(a)') 'Chem::State_Diag Finalize... FAILURE.'
       ENDIF
    ENDIF

    ! Deallocate fields of the Grid State object
    CALL Cleanup_State_Grid( State_Grid, RC )
    IF ( Input_Opt%AmIRoot ) THEN
       IF ( RC == GC_SUCCESS ) THEN
          write(*,'(a)') 'Chem::State_Grid Finalize... OK.'
       ELSE
          write(*,'(a)') 'Chem::State_Grid Finalize... FAILURE.'
       ENDIF
    ENDIF

    ! Deallocate fields of the Meteorology State object
    CALL Cleanup_State_Met( State_Met, RC )
    IF ( Input_Opt%AmIRoot ) THEN
       IF ( RC == GC_SUCCESS ) THEN
          write(*,'(a)') 'Chem::State_Met Finalize... OK.'
       ELSE
          write(*,'(a)') 'Chem::State_Met Finalize... FAILURE.'
       ENDIF
    ENDIF

    ! Deallocate fields of the Input Options object
    ! The call to Cleanup_Input_Opt causes a memory leak error. Comment
    ! for now (ckeller, 11/29/16).
    ! Does this still cause a memory leak? (ewl, 12/14/18)
     CALL Cleanup_Input_Opt( Input_Opt, RC )
    IF ( Input_Opt%AmIRoot ) THEN
       IF ( RC == GC_SUCCESS ) THEN
          write(*,'(a)') 'Chem::Input_Opt Finalize... OK.'
       ELSE
          write(*,'(a)') 'Chem::Input_Opt Finalize... FAILURE.'
       ENDIF
    ENDIF

    ! Free Int2Spc pointer
    IF ( ASSOCIATED(Int2Spc) ) THEN
       DO I=1,SIZE(Int2Spc,1)
          Int2Spc(I)%Internal => NULL()
       ENDDO
       DEALLOCATE(Int2Spc)
    ENDIF

    ! Deallocate the history interface between GC States and ESMF Exports
    CALL Destroy_HistoryConfig( am_I_Root, HistoryConfig, RC )

    ! Finalize MAPL Generic
    CALL MAPL_GenericFinalize( GC, Import, Export, Clock, __RC__ )

    ! Stop timers
    ! -----------
    !    CALL MAPL_TimerOff(STATE, "FINALIZE")
    !    CALL MAPL_TimerOff(STATE, "TOTAL")

    ! Successful return
    _RETURN(ESMF_SUCCESS)

  END SUBROUTINE Finalize_


  !------------------------------------------------------------------------------
  !                  GEOS-Chem Global Chemical Model                            !
  !------------------------------------------------------------------------------

  ! !DESCRIPTION: GC routine extracts several common quantities from the
  !  ESMF/MAPL environment so that they can be later passed down to the
  !  grid-independent GEOS-Chem code.

  SUBROUTINE Extract_( GC,         Clock,    Grid,    MaplCF, GeosCF,    &
                       localPet,   petCount,                             &
                       IM,         JM,       LM,                         &
                       IM_WORLD,   JM_WORLD, LM_WORLD,                   &
                       IL_WORLD,   IU_WORLD, JL_WORLD, JU_WORLD,         &
                       lonCtr,     latCtr,   advCount,                   &
                       nymdB,      nymdE,    nymd,    nhmsB,  nhmsE,     &
                       nhms,       year,     month,   day,    dayOfYr,   &
                       hour,       minute,   second,  utc,    hElapsed,  &
                       tsChem,     tsDyn,    mpiComm, ZTH,   SLR,        &
                       tsRad,                                            &
                       RC )

    TYPE(ESMF_Clock),    INTENT(IN)            :: Clock       ! ESMF clock obj
    TYPE(ESMF_GridComp), INTENT(INOUT)         :: GC          ! GC grid comp
    TYPE(ESMF_Grid),     INTENT(OUT), OPTIONAL :: Grid        ! ESMF Grid obj
    TYPE(ESMF_Config),   INTENT(OUT), OPTIONAL :: MaplCF      ! AGCM.rc
    TYPE(ESMF_Config),   INTENT(OUT), OPTIONAL :: GeosCF      ! GEOSCHEM*.rc
    INTEGER,             INTENT(OUT), OPTIONAL :: localPet    ! This PET
    INTEGER,             INTENT(OUT), OPTIONAL :: petCount    ! Total # of PETs
    INTEGER,             INTENT(OUT), OPTIONAL :: mpiComm     ! MPI Comm Handle

    ! Local grid coordinates
    ! (defined on the current CPU)
    INTEGER,             INTENT(OUT), OPTIONAL :: IM          ! Total # lons
    INTEGER,             INTENT(OUT), OPTIONAL :: JM          ! Total # lats
    INTEGER,             INTENT(OUT), OPTIONAL :: LM          ! Total # levs
    REAL(ESMF_KIND_R4),  POINTER,     OPTIONAL :: lonCtr(:,:) ! Lon ctrs [rad]
    REAL(ESMF_KIND_R4),  POINTER,     OPTIONAL :: latCtr(:,:) ! Lat ctrs [rad]

    ! Global grid coordinates
    INTEGER,             INTENT(OUT), OPTIONAL :: IM_WORLD    ! Global # lons
    INTEGER,             INTENT(OUT), OPTIONAL :: JM_WORLD    ! Global # lats
    INTEGER,             INTENT(OUT), OPTIONAL :: LM_WORLD    ! Global # levs
    INTEGER,             INTENT(OUT), OPTIONAL :: IL_WORLD    ! Global start lon index on this PET
    INTEGER,             INTENT(OUT), OPTIONAL :: IU_WORLD    ! Global end   lon index on this PET
    INTEGER,             INTENT(OUT), OPTIONAL :: JL_WORLD    ! Global start lat index on this PET
    INTEGER,             INTENT(OUT), OPTIONAL :: JU_WORLD    ! Global end   lat index on this PET

    ! Date and time variables
    INTEGER(ESMF_KIND_I8),INTENT(OUT), OPTIONAL :: advCount    ! # of clock advs
    INTEGER,             INTENT(OUT), OPTIONAL :: nymdB       ! YYYYMMDD @ start
    INTEGER,             INTENT(OUT), OPTIONAL :: nymdE       ! YYYYMMDD @ end
    INTEGER,             INTENT(OUT), OPTIONAL :: nymd        ! YYYYMMDD now
    INTEGER,             INTENT(OUT), OPTIONAL :: nhmsB       ! hhmmss @ start
    INTEGER,             INTENT(OUT), OPTIONAL :: nhmsE       ! hhmmss @ end
    INTEGER,             INTENT(OUT), OPTIONAL :: nhms        ! hhmmss now
    INTEGER,             INTENT(OUT), OPTIONAL :: year        ! UTC year
    INTEGER,             INTENT(OUT), OPTIONAL :: month       ! UTC month
    INTEGER,             INTENT(OUT), OPTIONAL :: day         ! UTC day
    INTEGER,             INTENT(OUT), OPTIONAL :: dayOfYr     ! UTC day of year
    INTEGER,             INTENT(OUT), OPTIONAL :: hour        ! UTC hour
    INTEGER,             INTENT(OUT), OPTIONAL :: minute      ! UTC minute
    INTEGER,             INTENT(OUT), OPTIONAL :: second      ! UTC second
    REAL,                INTENT(OUT), OPTIONAL :: utc         ! UTC time [hrs]
    REAL,                INTENT(OUT), OPTIONAL :: hElapsed    ! Elapsed hours

    ! Timestep variables [seconds]
    REAL,                INTENT(OUT), OPTIONAL :: tsChem      ! Chemistry
    REAL,                INTENT(OUT), OPTIONAL :: tsRad       ! RRTMG
    REAL,                INTENT(OUT), OPTIONAL :: tsDyn       ! Dynamics

    ! Solar parameters
    REAL,                INTENT(OUT), OPTIONAL :: ZTH(:,:)    ! Solar zth angle
    REAL,                INTENT(OUT), OPTIONAL :: SLR(:,:)    ! Insolation

    ! Return code
    INTEGER,             INTENT(OUT), OPTIONAL :: RC          ! 0 = all is well

    ! !REMARKS:
    !  If you need to obtain a quantity not returned by this routine, you can
    !  manually extract it from the MaplCF or GeosCF configuration objects.
!
    TYPE(ESMF_Time)               :: stopTime       ! ESMF stop time obj
    TYPE(ESMF_Time)               :: currTime       ! ESMF current time obj
    TYPE(ESMF_TimeInterval)       :: elapsedTime    ! ESMF elapsed time obj
    TYPE(ESMF_TimeInterval)       :: chemInterval   ! chemistry interval
    TYPE(ESMF_ALARM)              :: ALARM          ! Run alarm
    TYPE(ESMF_VM)                 :: VM             ! ESMF VM object
    TYPE(GEOSCHEM_State), POINTER :: myState        ! Legacy state
    TYPE(GEOSCHEM_Wrap)           :: wrap           ! Wrapper for myState
    TYPE(MAPL_MetaComp),  POINTER :: STATE          ! MAPL MetaComp object
    TYPE(MAPL_SunOrbit)           :: sunOrbit
    CHARACTER(len=ESMF_MAXSTR)    :: compName       ! Gridded component name
    CHARACTER(len=ESMF_MAXSTR)    :: importRstFN    ! Import restart file name
    INTEGER(ESMF_KIND_I8)         :: count          ! # of clock advances
    INTEGER                       :: locDims(3)     ! Array for local dims
    INTEGER                       :: globDims(3)    ! Array for global dims
    INTEGER                       :: doy            ! Day of year (0-365/366)
    INTEGER                       :: yyyy, mm, dd   ! Year, month, day
    INTEGER                       :: h,    m,  s    ! Hour, minute, seconds
    INTEGER                       :: IL,   IU       ! Min/max local lon indices
    INTEGER                       :: JL,   JU       ! Min/max local lat indices
    REAL                          :: elapsedHours   ! Elapsed hours of run
    REAL(ESMF_KIND_R8)            :: dt_r8          ! chemistry timestep
    CHARACTER(len=ESMF_MAXSTR)    :: OUTSTR         ! Parallel write nonsense

    ! Saved variables
    LOGICAL, SAVE                 :: FIRST = .TRUE.
    TYPE(ESMF_Time), SAVE         :: startTime

    !=======================================================================
    ! Initialization
    !=======================================================================

    ! Get my name and set-up traceback handle
    CALL ESMF_GridCompGet( GC, name=compName, vm=VM, __RC__ )

    ! Get the internal state which holds the private Config object
    CALL ESMF_UserCompGetInternalState( GC, 'GEOSCHEM_State', wrap, STATUS )
    _VERIFY(STATUS)
    myState => wrap%ptr

    ! Get generic state object
    CALL MAPL_GetObjectFromGC( GC, STATE, __RC__ )

    ! Assume successful return
    IF ( PRESENT( RC ) ) RC = ESMF_SUCCESS

    ! Zero variables
    locDims  = 0
    globDims = 0
    IL       = 0
    JL       = 0
    IU       = 0
    JU       = 0

    !=======================================================================
    ! Extract information from ESMF VM object
    !=======================================================================


    ! Index of the PET we are on now
    IF ( PRESENT( localPet ) ) THEN
       CALL ESMF_VmGet( VM, localPet=localPet, __RC__ )
    ENDIF

    ! Total # of PETs used by this gridded component
    IF ( PRESENT( petCount ) ) THEN
       CALL ESMF_VmGet( VM, petCount=petCount, __RC__ )
    ENDIF

    ! Global MPI Communicator Handle
    IF ( PRESENT( mpiComm ) ) THEN
       CALL ESMF_VmGet( VM, mpicommunicator=mpiComm, __RC__ )
    ENDIF

    !=======================================================================
    ! Extract information from ESMF Config objects
    !=======================================================================

    ! Get the Config object
    CALL ESMF_GridCompGet( GC, Config=MaplCF, __RC__ )

    ! Get the Config object based on "GEOSCHEMchem_GridComp.rc"
    GeosCF = myState%myCF

    ! Dynamic timestep (in seconds)
    IF ( PRESENT( tsDyn ) ) THEN
       CALL ESMF_ConfigGetAttribute( MaplCF, tsDyn, Default=1800.,        &
                                     Label="RUN_DT:",             __RC__ )
    ENDIF

    ! Radiation timestep (in seconds)
    IF ( PRESENT( tsRad ) ) THEN
       CALL ESMF_ConfigGetAttribute( MaplCF, tsRad, Default=10800.,        &
                                     Label="RRTMG_DT:",             __RC__ )
    ENDIF

    ! Chemistry timestep (in seconds)
    IF ( PRESENT( tsChem ) ) THEN
        CALL MAPL_Get( STATE, RUNALARM=ALARM, __RC__ )
        CALL ESMF_AlarmGet( ALARM, RingInterval=chemInterval, __RC__ )
        CALL ESMF_TimeIntervalGet( chemInterval, s_r8=dt_r8, __RC__ )
        tsChem = real(dt_r8)

        IF(abs(tsChem) < abs(tsDyn)) THEN
           IF( MAPL_AM_I_ROOT() ) THEN
              WRITE(6,*) 'Chem_DT cannot be less than RUN_DT'
           ENDIF
           STATUS = 1
           _VERIFY(STATUS)
        ENDIF
    ENDIF

    If ( PRESENT( tsRad ) .and. PRESENT( tsChem ) ) Then
        _ASSERT(MOD(nint(tsRad),nint(tsChem)) == 0,'RRTMG_DT is not a multiple of GCHPCHEM_DT')
    End If

    !=======================================================================
    ! Extract time/date information
    !=======================================================================

    ! Get the ESMF time object
    CALL ESMF_ClockGet( Clock,                    &
                        stopTime     = stopTime,  &
                        currTime     = currTime,  &
                        advanceCount = count,     &
                         __RC__ )

    !=======================================================================
    ! Current, start, and end times
    !=======================================================================

    ! Get current-time fields from the time object. Set start/end if first.
    CALL ESMF_TimeGet( currTime, yy=yyyy, mm=mm, dd=dd, dayOfYear=doy, &
                                 h=h,     m=m,   s=s,   __RC__ )
    IF ( PRESENT( nymd     ) ) CALL MAPL_PackTime( nymd, yyyy, mm, dd )
    IF ( PRESENT( nhms     ) ) CALL MAPL_PackTime( nhms, h,    m,  s  )
    IF ( PRESENT( year     ) ) year     = yyyy
    IF ( PRESENT( month    ) ) month    = mm
    IF ( PRESENT( day      ) ) day      = dd
    IF ( PRESENT( dayOfYr  ) ) dayOfYr  = doy
    IF ( PRESENT( hour     ) ) hour     = h
    IF ( PRESENT( minute   ) ) minute   = m
    IF ( PRESENT( second   ) ) second   = s
    IF ( PRESENT( utc      ) ) utc      = ( DBLE( h )        ) + &
                                          ( DBLE( m )/60d0   ) + &
                                          ( DBLE( s )/3600d0 )

    ! Simulation start
    IF ( FIRST ) THEN
       startTime = currTime
    ENDIF
    CALL ESMF_TimeGet( startTime, yy=yyyy, mm=mm, dd=dd, dayOfYear=doy, &
                                  h=h,     m=m,   s=s,   __RC__ )
    IF ( PRESENT ( nymdB ) ) CALL MAPL_PackTime( nymdB, yyyy, mm, dd )
    IF ( PRESENT ( nhmsB ) ) CALL MAPL_PackTime( nhmsB, h,    m,  s  )

    ! Simulation end
    CALL ESMF_TimeGet( stopTime, yy=yyyy, mm=mm, dd=dd, dayOfYear=doy, &
                                 h=h,     m=m,   s=s,   __RC__ )
    IF ( PRESENT( nymdE ) ) CALL MAPL_PackTime( nymdE, yyyy, mm, dd )
    IF ( PRESENT( nhmsE ) ) CALL MAPL_PackTime( nhmsE, h,    m,  s  )

    ! # clock steps
    IF ( PRESENT( advCount ) ) advCount = count

    ! Compute elapsed time since start of simulation
    elapsedTime = currTime - startTime

    ! Get time fields from the elapsedTime object
    CALL ESMF_TimeIntervalGet( elapsedTime, h=h, m=m, s=s, __RC__ )

    ! Convert to decimal hours
    elapsedHours = DBLE( h ) + ( DBLE( m )/60d0 ) + ( DBLE( s )/3600d0 )

    ! Save fields for return
    IF ( PRESENT( hElapsed ) ) hElapsed = elapsedHours

    !=======================================================================
    ! Extract grid information
    !=======================================================================
    IF ( PRESENT( Grid ) ) THEN

       ! Get the ESMF grid attached to this gridded component
       CALL ESMF_GridCompGet( GC, grid=Grid, __RC__ )

       ! Get # of dimensions on this pet, and globally
       CALL MAPL_GridGet( Grid,                                        &
                          localCellCountPerDim  = locDims,             &
                          globalCellCountPerDim = globDims,            &
                          __RC__ )

       ! Get the upper and lower bounds of on each PET using MAPL
       CALL MAPL_GridGetInterior( Grid, IL, IU, JL, JU )
       ! if (PRESENT(localPet)) THEN
       !    WRITE (*,1141) localPet, IL, IU, JL, JU
       ! endif

1141   FORMAT(' Process ', i5, ' goes from I = ', i3, ':', i3, '   J = ', i3, ':', i3)

    ENDIF

    ! Save fields for return
    IF ( PRESENT( IM       ) ) IM       = locDims(1)
    IF ( PRESENT( JM       ) ) JM       = locDims(2)
    IF ( PRESENT( LM       ) ) LM       = locDims(3)
    IF ( PRESENT( IM_WORLD ) ) IM_WORLD = globDims(1)
    IF ( PRESENT( JM_WORLD ) ) JM_WORLD = globDims(2)
    IF ( PRESENT( LM_WORLD ) ) LM_WORLD = globDims(3)

    IF ( PRESENT( IL_WORLD ) ) IL_WORLD = IL
    IF ( PRESENT( IU_WORLD ) ) IU_WORLD = IU
    IF ( PRESENT( JL_WORLD ) ) JL_WORLD = JL
    IF ( PRESENT( JU_WORLD ) ) JU_WORLD = JU

    ! Longitude values on this PET
    IF ( PRESENT( lonCtr ) ) THEN
       CALL MAPL_Get( STATE, lons=lonCtr, __RC__ )
    ENDIF

    ! Latitude values on this PET
    IF ( PRESENT( latCtr ) ) THEN
       CALL MAPL_Get( STATE, lats=latCtr, __RC__ )
    ENDIF

    !=======================================================================
    ! Get solar zenith angle enformation
    !=======================================================================
    IF ( PRESENT( ZTH    ) .and. PRESENT( SLR    )  .and. &
         PRESENT( lonCtr ) .and. PRESENT( latCtr ) ) THEN

       ! Get the Orbit object (of type MAPL_SunOrbit),
       ! which is used in the call to MAPL_SunGetInsolation
       CALL MAPL_Get( STATE,                       &
                      LONS      = lonCtr,             &
                      LATS      = latCtr,             &
                      ORBIT     = sunOrbit,           &
                      __RC__                         )

       ! Get the solar zenith angle and solar insolation
       ! NOTE: ZTH, SLR are allocated outside of this routine
       CALL MAPL_SunGetInsolation( LONS  = lonCtr,    &
                                   LATS  = latCtr,    &
                                   ORBIT = sunOrbit,  &
                                   ZTH   = ZTH,       &
                                   SLR   = SLR,       &
                                   CLOCK = Clock,     &
                                   __RC__            )

    ENDIF

    FIRST = .FALSE.
    _RETURN(ESMF_SUCCESS)

  END SUBROUTINE Extract_

END MODULE Chem_GridCompMod
