#include "MAPL_Generic.h"

!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Model                            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: GEOSCHEMchem_GridCompMod
!
! !DESCRIPTION: GEOSCHEMchem_GridComp is an ESMF5 gridded component
! implementing the GEOS-Chem chemistry and related processes, including
! dry deposition, emissions, and wet deposition. In addition, the
! parameterizations for PBL mixing and convection as used in GEOS-Chem
! can be invoked by enabling the corresponding option in the GEOS-Chem
! input file (geoschem_config.yml). In this case, the corresponding GEOS-5
! process must NOT be applied to the GC tracers, i.e. the tracers must
! not be friendly to turbulence (if PBL mixing is used) and/or moist
! (for convection).
!\\
!\\
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
!\\
!\\
! !INTERFACE:
!
#ifdef MODEL_GEOS
MODULE GEOSCHEMchem_GridCompMod
#else
MODULE Chem_GridCompMod
#endif
!
! !USES:
!
  USE CMN_Size_Mod
  USE ESMF                                           ! ESMF library
  USE MAPL_Mod                                       ! MAPL library
  USE MAPL_IOMod
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

#if defined( MODEL_GEOS )
  USE Chem_Mod                            ! Chemistry Base Class (chem_mie?)
  USE Chem_GroupMod                       ! For family transport
  USE PHYSCONSTANTS
#endif

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC   :: SetServices    ! Sets ESMF entry points
!
! !PRIVATE MEMBER FUNCTIONS:
!
  PRIVATE  :: Initialize_    ! Init method
  PRIVATE  :: Run1           ! Run wrapper phase 1
  PRIVATE  :: Run2           ! Run wrapper phase 2
  PRIVATE  :: Run_           ! Run method
  PRIVATE  :: Finalize_      ! Finalize method
  PRIVATE  :: Extract_       ! Get values from ESMF
!
! !PRIVATE TYPES:
!
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
#if defined( MODEL_GEOS )
     REAL, POINTER                 :: Internal(:,:,:) => NULL()
#else
     REAL(ESMF_KIND_R8), POINTER   :: Internal(:,:,:) => NULL()
#endif
  END TYPE Int2SpcMap

  ! Internal run alarms
  TYPE GC_run_alarms
     private
     ! Add alarms here
     type(ESMF_Alarm) :: RRTMG_Alarm
  END TYPE GC_run_alarms

  TYPE GCRA_wrap
     type(GC_run_alarms), pointer  :: ptr
  END TYPE

  ! For mapping State_Chm%Tracers/Species arrays onto the internal state.
  TYPE(Int2SpcMap), POINTER        :: Int2Spc(:) => NULL()
#ifdef ADJOINT
  ! For mapping State_Chm%Tracers/Species arrays onto the internal state.
  TYPE(Int2SpcMap), POINTER        :: Int2Adj(:) => NULL()
#endif

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
  LOGICAL                          :: meteorology_vertical_index_is_top_down

#if defined( MODEL_GEOS )
  ! Is GEOS-Chem the provider for AERO, RATS, and/or Analysis OX?
  LOGICAL                          :: DoAERO

  ! When to do the analysis
  INTEGER                          :: ANAPHASE
  INTEGER, PARAMETER               :: CHEMPHASE = 2
#endif

  ! Number of run phases, 1 or 2. Set in the rc file; else default is 2.
  INTEGER                          :: NPHASE

  ! Is this being run as a CTM?
  INTEGER                          :: IsCTM

  ! Memory debug level
  INTEGER                          :: MemDebugLevel

#if defined( MODEL_GEOS )
  ! GEOS-5 only
  ! Flag to initialize species concentrations from external fields. Read
  ! through GEOSCHEMchem_GridComp.rc. If true, initial species concentrations
  ! are read from an external file instead of taken from the internal state.
  ! The field names in the external file are expected to be 'SPC_<XXX>'.
  ! This option can be used to initialize a simulation using a restart file
  ! from a 'GEOS-Chem classic' CTM simulation.
  LOGICAL                          :: InitFromFile
  LOGICAL                          :: SkipReplayGCC
#endif

  ! Pointers to import, export and internal state data. Declare them as
  ! module variables so that we have to assign them only on first call.

#if defined( MODEL_GEOS )
# include "GEOSCHEMCHEM_DeclarePointer___.h"
#else
# include "GCHPchem_DeclarePointer___.h"
#endif

#if defined( MODEL_GEOS )

  ! Mie table
!  TYPE(Chem_Mie)     :: geoschemMieTable(2)
!  INTEGER, PARAMETER :: instanceComputational = 1
!  INTEGER, PARAMETER :: instanceData          = 2
#endif
!
! !REMARKS:
!  Developed for GEOS-5 release Fortuna 2.0 and later.
!                                                                             .
!  NOTES:
!  - The abbreviation "PET" stands for "Persistent Execution Thread".
!    It is a synomym for CPU.
!
! !REVISION HISTORY:
!  06 Dec 2009 - A. da Silva - Initial version
!  See https://github.com/geoschem/geos-chem for history
!EOP
!------------------------------------------------------------------------------
!BOC
CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Model                            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: SetServices
!
! !DESCRIPTION: The SetServices routine does the following:
!
! \begin{itemize}
! \item Defines the Initialize method for the GEOSCHEMchem gridded component
! \item Defines the Run methods for the GEOSCHEMchem gridded component
!       (phase 1 and phase 2).
! \item Defines the Finalize method for the GEOSCHEMchem gridded component
! \item Attaches an internal state (which holds a private ESMF Config object)
!       to the GEOSCHEMchem gridded component.
! \end{itemize}
!
! !INTERFACE:
!
  SUBROUTINE SetServices( GC, RC )
!
! !USES:
!
    USE HCOI_ESMF_MOD,        ONLY : HCO_SetServices
    USE GCKPP_Model
    USE CHARPAK_MOD,          ONLY : STRSPLIT, CSTRIP
    USE inquireMod,           ONLY : findFreeLUN
    USE FILE_MOD,             ONLY : IOERROR
#if defined( MODEL_GEOS )
    USE CMN_FJX_MOD
    USE GCKPP_Monitor
    USE GCKPP_Parameters
    USE Precision_Mod
    USE GEOS_Analysis,        ONLY : GEOS_AnaInit
    USE GEOS_Interface,       ONLY : MetVars_For_Lightning_Init, &
                                     GEOS_CheckRATSandOx 
    USE GEOS_AeroCoupler,     ONLY : GEOS_AeroSetServices
    USE GEOS_CarbonInterface, ONLY : GEOS_CarbonSetServices
#endif
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ESMF_GridComp), INTENT(INOUT) :: GC       ! Ref to this GridComp
!
! !OUTPUT PARAMETERS:
!
    INTEGER,             INTENT(OUT)   :: RC       ! Success or failure
!
! !REMARKS:
!  ESMF can only attach one Config object per Gridded Component.  The
!  Config object that is defined from the "MAPL.rc" resource file is
!  directly attached to the GEOSCHEMchem gridded component.
!                                                                             .
!  To attach the Config object defined from the "GEOSCHEMchem_GridComp.rc"
!  resource file, we must first create a derived type with a pointer to
!  the Config object, then attach that to the gridded component as an
!  "internal state" (also called "legacy state").
!
! !REVISION HISTORY:
!  06 Dec 2009 - A. da Silva - Initial version
!  See https://github.com/geoschem/geos-chem for history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
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
    INTEGER                       :: restartAttr
    CHARACTER(LEN=ESMF_MAXSTR)    :: HistoryConfigFile ! HISTORY config file
    INTEGER                       :: T

#if !defined( MODEL_GEOS )
    TYPE(MAPL_MetaComp),  POINTER :: STATE => NULL()
#endif

#if defined( MODEL_GEOS )
    CHARACTER(LEN=ESMF_MAXSTR)    :: LongName      ! Long name for diagnostics
    CHARACTER(LEN=ESMF_MAXSTR)    :: ShortName
    CHARACTER(LEN=255)            :: MYFRIENDLIES
    CHARACTER(LEN=127)            :: FullName
    INTEGER                       :: DoIt
    LOGICAL                       :: FriendMoist, SpcInRestart, ReduceSpc
    CHARACTER(LEN=40)             :: SpcsBlacklist(255)
    INTEGER                       :: nBlacklist
    CHARACTER(LEN=ESMF_MAXSTR)    :: Blacklist
#endif
#ifdef ADJOINT
    INTEGER                       :: restartAttrAdjoint
    LOGICAL                       :: useCFMaskFile
#endif

    ! Manual internal state entries
    LOGICAL                       :: am_I_Root
    INTEGER                       :: II
    CHARACTER(LEN=2)              :: intStr
    CHARACTER(LEN=ESMF_MAXSTR)    :: myName

    __Iam__('SetServices')


    lgr => logging%get_logger('GCHPchem')

    !=======================================================================
    ! Set services begins here
    !=======================================================================

    ! Set up traceback info
    CALL ESMF_GridCompGet( GC, name=compName, __RC__ )

    ! NOTE: We need to use COMP_NAME for mapl_acg.pl script
    COMP_NAME = TRIM( compName )

    ! Identify this routine to MAPL
    Iam = TRIM(compName)//'::SetServices'

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
    myState%myCF = ESMF_ConfigCreate(__RC__)

#if defined( MODEL_GEOS )
    call ESMF_ConfigLoadFile( myState%myCF, 'GEOSCHEMchem_GridComp.rc', __RC__)
#else
    call ESMF_ConfigLoadFile( myState%myCF, 'GCHP.rc', __RC__)
#endif

#if !defined( MODEL_GEOS )
    ! Get generic state object
    CALL MAPL_GetObjectFromGC( GC, STATE, __RC__ )
    call MAPL_GetResource( STATE, IsCTM, label='GEOSChem_CTM:', &
                           default=1, rc=status )
    _VERIFY(STATUS)
#endif

    !=======================================================================
    !                 %%% ESMF Functional Services %%%
    !=======================================================================

    ! Set the Initialize, Run, Finalize entry points
    CALL MAPL_GridCompSetEntryPoint( GC, ESMF_METHOD_INITIALIZE,  &
                                     Initialize_, __RC__ )
#if defined( MODEL_GEOS )
    CALL MAPL_GridCompSetEntryPoint( GC, ESMF_METHOD_RUN, Run1, __RC__ )
#endif
    CALL MAPL_GridCompSetEntryPoint( GC, ESMF_METHOD_RUN, Run2, __RC__ )
    CALL MAPL_GridCompSetEntryPoint( GC, ESMF_METHOD_FINALIZE,  &
                                     Finalize_, __RC__ )

    ! Store internal state with Config object in the gridded component
    CALL ESMF_UserCompSetInternalState( GC, 'GEOSCHEM_State', wrap, STATUS )
    _VERIFY(STATUS)

#if defined(MODEL_GEOS)
    CALL GEOS_AeroSetServices  ( GC, DoAERO, __RC__ )
    CALL GEOS_CarbonSetServices( GC, myState%myCF, __RC__ )
#endif

    !=======================================================================
    !                    %%% MAPL Data Services %%%
    !=======================================================================
!EOC
!BOP
!
! !IMPORT STATE:
!
#if defined( MODEL_GEOS )
#   include "GEOSCHEMCHEM_ImportSpec___.h"
#else
#   include "GCHPchem_ImportSpec___.h"
#endif

#if !defined( MODEL_GEOS )
    call MAPL_AddImportSpec(GC, &
       SHORT_NAME         = 'PLE',  &
       LONG_NAME          = 'pressure_level_edges',  &
       UNITS              = 'Pa', &
       PRECISION          = ESMF_KIND_R8, &
       DIMS               = MAPL_DimsHorzVert,    &
       VLOCATION          = MAPL_VLocationEdge,    &
                                                      RC=STATUS  )
    _VERIFY(STATUS)
#endif

    !=======================================================================
    ! Get meteorology vertical index orientation
    !=======================================================================
    call ESMF_ConfigGetAttribute(myState%myCF,value=meteorology_vertical_index_is_top_down, &
    label='METEOROLOGY_VERTICAL_INDEX_IS_TOP_DOWN:', Default=.false., __RC__ )
    if (meteorology_vertical_index_is_top_down) then
       call lgr%info('Configured to expect ''top-down'' meteorological data from ''ExtData''')
    else
       call lgr%info('Configured to expect ''bottom-up'' meteorological data from ''ExtData''')
    end if

#if defined( MODEL_GEOS )
    ! Define imports to fill met fields needed for lightning
    CALL MetVars_For_Lightning_Init( GC, MyState%myCF, __RC__ )
#endif

#ifdef ADJOINT
    CALL ESMF_ConfigGetAttribute( myState%myCF, useCFMaskFile, &
         Label="USE_CF_MASK_FILE:", Default=.false., __RC__ )

    IF (useCFMaskFile) THEN
       call MAPL_AddImportSpec(GC,                    &
            SHORT_NAME         = 'CFN_MASK',            &
            LONG_NAME          = 'cost_function_Mask',  &
            UNITS              = '1',                   &
            DIMS               = MAPL_DimsHorzVert,     &
            VLOCATION          = MAPL_VLocationCenter,  &
            RC=STATUS  )
       _VERIFY(STATUS)
    Endif
#endif

!
! !INTERNAL STATE:
!
#if defined( MODEL_GEOS )
#   include "GEOSCHEMCHEM_InternalSpec___.h"
#else
#   include "GCHPchem_InternalSpec___.h"
#endif

#if !defined( MODEL_GEOS )
    ! Determine if using a restart file for the internal state. Setting
    ! the GCHPchem_INTERNAL_RESTART_FILE to +none in GCHP.rc indicates
    ! skipping the restart file. Species concentrations will be retrieved
    ! from the species database, overwriting MAPL-assigned default values.
    CALL ESMF_ConfigGetAttribute( myState%myCF, rstFile, &
                                  Label = "GCHPchem_INTERNAL_RESTART_FILE:",&
                                  __RC__ )
    IF ( TRIM(rstFile) == '+none' ) THEN
       restartAttr = MAPL_RestartSkipInitial ! file does not exist;
                                             ! use background values
    ELSE
       restartAttr = MAPL_RestartOptional    ! try to read species from file;
                                             ! use background vals if not found
    ENDIF
#ifdef ADJOINT
    restartAttrAdjoint = MAPL_RestartSkip
#endif
#endif

!-- Read in species from geoschem_config.yml and set FRIENDLYTO
    ! ewl TODO: This works but is not ideal. Look into how to remove it.

#if defined( MODEL_GEOS )
    ! Check if species are friendly to moist
    CALL ESMF_ConfigGetAttribute( myState%myCF, DoIt, &
                                  Label = "Species_friendly_to_moist:",&
                                  Default = 0, &
                                  __RC__ )
    FriendMoist = (DoIt==1)
    IF ( MAPL_am_I_Root() ) THEN
       WRITE(*,*) 'GCC species friendly to MOIST: ',FriendMoist
    ENDIF
    ! Determine if non-advected species shall be included in restart file
    CALL ESMF_ConfigGetAttribute( myState%myCF, DoIt, &
                                  Label = "Shortlived_species_in_restart:", &
                                  Default = 1, __RC__ )
    IF ( DoIt==1 ) THEN
       restartAttr  = MAPL_RestartOptional
       SpcInRestart = .TRUE.
    ELSE
       restartAttr  = MAPL_RestartSkip
       SpcInRestart = .FALSE.
    ENDIF
    ! Check if we want to use a reduced set of species for transport
    SpcsBlacklist(:) = ''
    nBlacklist = 0
    CALL ESMF_ConfigGetAttribute( myState%myCF, DoIt, &
                                  Label = "Reduce_transport_species:", &
                                  Default = 0, __RC__ )
    ReduceSpc = ( DoIt==1 )
    ! Get list of blacklisted species
    IF ( ReduceSpc ) THEN
       CALL ESMF_ConfigGetAttribute( myState%myCF, Blacklist, &
                                     Label = "Transport_blacklist:", &
                                     Default = 'CFC11,CFC12', __RC__ )
       IF ( TRIM(ADJUSTL(Blacklist)) /= '' ) THEN
          CALL STRSPLIT( Blacklist, ',', SpcsBlacklist, nBlacklist )
       ENDIF
    ENDIF

    ! Sulfur-nitrogen-ammonia water content computed in Isorropia after needed in RDAER
    call MAPL_AddInternalSpec(GC, &
       SHORT_NAME         = 'AeroH2O_SNA',  &
       LONG_NAME          = 'Sulfur-nitrogen-ammonia water content',  &
       UNITS              = 'g/m3', &
       DIMS               = MAPL_DimsHorzVert,    &
       VLOCATION          = MAPL_VLocationCenter,    &
!       PRECISION          = ESMF_KIND_R8, &
       FRIENDLYTO         = 'GEOSCHEMCHEM',    __RC__ )

#endif

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

#if defined( MODEL_GEOS )
          ! %%% GEOS-Chem in GEOS %%%
          ! Define friendliness to dynamics / turbulence
          MYFRIENDLIES = 'DYNAMICS:TURBULENCE'
          IF ( FriendMoist ) THEN
             MYFRIENDLIES = TRIM(MYFRIENDLIES)//':MOIST'
          ENDIF
          FullName = TRIM(SUBSTRS(1))
          ! Check if this species is blacklisted
          IF ( nBlacklist > 0 ) THEN
             DO I=1,nBlacklist
                IF ( TRIM(SpcsBlacklist(I))==TRIM(FullName) ) THEN
                   MYFRIENDLIES = TRIM(COMP_NAME)
                   EXIT
                ENDIF
             ENDDO
          ENDIF
          CALL MAPL_AddInternalSpec(GC,                                     &
               SHORT_NAME      = TRIM(SPFX)//TRIM(SUBSTRS(1)),              &
               LONG_NAME       = TRIM(FullName)//                           &
                                 ' mass mixing ratio total air',            &
               UNITS           = 'kg kg-1',                                 &
               DIMS            = MAPL_DimsHorzVert,                         &
               VLOCATION       = MAPL_VLocationCenter,                      &
            !!!PRECISION       = ESMF_KIND_R8,                              &
               FRIENDLYTO      = TRIM(MYFRIENDLIES),                        &
               RC              = RC                                        )

          ! Update count of advected species
          Nadv = Nadv + 1
          AdvSpc(Nadv) = TRIM(SUBSTRS(1))

          ! verbose
          if(MAPL_am_I_Root()) write(*,*)                                   &
               'GCC added to internal: '//TRIM(SPFX)//TRIM(SUBSTRS(1)),     &
               '; Friends: ', TRIM(MYFRIENDLIES)
#else
          !%%% GEOS-Chem in GCHP ###
          CALL MAPL_AddInternalSpec(GC, &
               SHORT_NAME      = TRIM(SPFX) // TRIM(SUBSTRS(1)),            &
               LONG_NAME       = TRIM(SUBSTRS(1)),                          &
               UNITS           = 'mol mol-1',                               &
               DIMS            = MAPL_DimsHorzVert,                         &
               VLOCATION       = MAPL_VLocationCenter,                      &
               PRECISION       = ESMF_KIND_R8,                              &
               FRIENDLYTO      = 'DYNAMICS:TURBULENCE:MOIST',               &
               RESTART         = restartAttr,                               &
               RC              = RC                                       )

          ! Add to list of transported speces
          NADV = NADV + 1
          AdvSpc(NADV) = TRIM(SUBSTRS(1))
#endif
#ifdef ADJOINT
          !%%%% GEOS-Chem in GCHP Adjoint %%%%
          if (MAPL_am_I_Root())                                             &
               WRITE(*,*) '  Adding internal spec for '''//TRIM(SPFX) //    &
               TRIM(SUBSTRS(1)) // '_ADJ'''
          call MAPL_AddInternalSpec(GC, &
               SHORT_NAME      = TRIM(SPFX) // TRIM(SUBSTRS(1)) // '_ADJ',  &
               LONG_NAME       = TRIM(SUBSTRS(1)) // ' adjoint variable',   &
               UNITS           = 'mol mol-1',                               &
               DIMS            = MAPL_DimsHorzVert,                         &
               VLOCATION       = MAPL_VLocationCenter,                      &
               PRECISION       = ESMF_KIND_R8,                              &
               FRIENDLYTO      = 'DYNAMICS:TURBULENCE:MOIST',               &
               RESTART         = restartAttrAdjoint,                        &
               RC              = RC                                        )
#endif
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

#if defined( MODEL_GEOS )
          ! Don't need to do anything if short-lived species are not in
          ! restart file
          IF ( .NOT. SpcInRestart ) CYCLE
#endif

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
#if defined( MODEL_GEOS )
             !%%%% GEOS-Chem in GEOS %%%%
             FullName = TRIM(SpcName)

             ! Error trap for POx and LOx. Their species names in the internal
             ! state must be all caps
             ! (ckeller, 3/11/19)
             !IF ( TRIM(SpcName) == 'POx' ) SpcName = 'POX'
             !IF ( TRIM(SpcName) == 'LOx' ) SpcName = 'LOX'

             ! Set some long names manually ...
             SELECT CASE ( TRIM(SpcName) )
                CASE ('OH')
                   FullName = 'Hydroxyl radical (OH, MW = 17.01 g mol-1)'
                CASE ('HO2')
                   FullName = 'Hydroperoxyl radical (HO2, MW = 33.01 g mol-1)'
                CASE ('O')
                   FullName = 'Molecular oxygen (O, MW = 16.01 g mol-1)'
             END SELECT

             CALL MAPL_AddInternalSpec(GC,                                   &
                  SHORT_NAME      = TRIM(SPFX)//TRIM(SpcName),               &
                  LONG_NAME       = TRIM(FullName)//                         &
                                    ' mass mixing ratio total air',          &
                  UNITS           = 'kg kg-1',                               &
               !!!PRECISION       = ESMF_KIND_R8,                            &
                  DIMS            = MAPL_DimsHorzVert,                       &
                  FRIENDLYTO      = COMP_NAME,                               &
                  RESTART         = restartAttr,                             &
                  VLOCATION       = MAPL_VLocationCenter,                    &
                                   __RC__                                   )
             ! verbose
             if(MAPL_am_I_Root()) write(*,*)  &
                  'GCC added to internal: '//TRIM(SPFX)//TRIM(SpcName)
#else
             !%%%% GEOS-Chem in GCHP %%%%
             call MAPL_AddInternalSpec(GC, &
                  SHORT_NAME      = TRIM(SPFX) // SpcName,                   &
                  LONG_NAME       = SpcName,                                 &
                  UNITS           = 'mol mol-1',                             &
                  PRECISION       = ESMF_KIND_R8,                            &
                  DIMS            = MAPL_DimsHorzVert,                       &
                  VLOCATION       = MAPL_VLocationCenter,                    &
                  RESTART         = restartAttr,                             &
                  RC              = STATUS                                  )
#ifdef ADJOINT
             !%%%% GEOS-Chem in GCHP ADJOINT %%%%
             if (MAPL_am_I_Root()) &
                  WRITE(*,*) '  Adding internal spec for '''//TRIM(SPFX) //  &
                             TRIM(SpcName) // '_ADJ'''
             call MAPL_AddInternalSpec(GC,                                   &
                  SHORT_NAME      = TRIM(SPFX) // TRIM(SpcName) // '_ADJ',   &
                  LONG_NAME       = SpcName // ' adjoint variable',          &
                  UNITS           = 'mol mol-1',                             &
                  PRECISION       = ESMF_KIND_R8,                            &
                  DIMS            = MAPL_DimsHorzVert,                       &
                  VLOCATION       = MAPL_VLocationCenter,                    &
                  RESTART         = restartAttrAdjoint,                      &
                  RC              = STATUS                                  )
#endif
#endif
          ENDIF
       ENDDO
    ENDIF

#if !defined( MODEL_GEOS )
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

    ! Sulfur-nitrogen-ammonia water content computed in Isorropia after needed in RDAER
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
#endif

#if defined( MODEL_GEOS )
!-- Add two exra advected species for use in family transport  (Manyin)

          CALL MAPL_AddInternalSpec(GC,                                    &
             SHORT_NAME         = 'SPC_Bry',                               &
             LONG_NAME          = 'Bromine group for use in transport',    &
             UNITS              = 'kg kg-1',                               &
!!!          PRECISION          = ESMF_KIND_R8,                            &
             DIMS               = MAPL_DimsHorzVert,                       &
             FRIENDLYTO         = 'DYNAMICS',                              &
             RESTART            = MAPL_RestartSkip,                        &
             VLOCATION          = MAPL_VLocationCenter,                    &
                                                  __RC__ )
          if(MAPL_am_I_Root()) write(*,*) 'GCC added to internal: SPC_Bry; Friendly to: DYNAMICS'

          CALL MAPL_AddInternalSpec(GC,                                    &
             SHORT_NAME         = 'SPC_Cly',                               &
             LONG_NAME          = 'Chlorine group for use in transport',   &
             UNITS              = 'kg kg-1',                               &
!!!          PRECISION          = ESMF_KIND_R8,                            &
             DIMS               = MAPL_DimsHorzVert,                       &
             FRIENDLYTO         = 'DYNAMICS',                              &
             RESTART            = MAPL_RestartSkip,                        &
             VLOCATION          = MAPL_VLocationCenter,                    &
                                                  __RC__ )
          if(MAPL_am_I_Root()) write(*,*) 'GCC added to internal: SPC_Cly; Friendly to: DYNAMICS'
!
        ! Add additional RATs/ANOX exports
!        call MAPL_AddExportSpec(GC,                                  &
!           SHORT_NAME         = 'OX_TEND',                           &
!           LONG_NAME          = 'tendency_of_odd_oxygen_mixing_ratio_due_to_chemistry', &
!           UNITS              = 'mol mol-1 s-1',                     &
!           DIMS               = MAPL_DimsHorzVert,                   &
!           VLOCATION          = MAPL_VLocationCenter,                &
!                                                     __RC__ )

     call MAPL_AddExportSpec(GC,                                     &
           SHORT_NAME         = 'GCC_O3',                            &
           LONG_NAME          = 'ozone_mass_mixing_ratio_total_air', &
           UNITS              = 'kg kg-1',                           &
           DIMS               = MAPL_DimsHorzVert,                   &
           VLOCATION          = MAPL_VLocationCenter,                &
                                                     __RC__ )

     call MAPL_AddExportSpec(GC,                                       &
           SHORT_NAME         = 'GCC_O3PPMV',                          &
           LONG_NAME          = 'ozone_volume_mixing_ratio_total_air', &
           UNITS              = 'ppmv',                                &
           DIMS               = MAPL_DimsHorzVert,                     &
           VLOCATION          = MAPL_VLocationCenter,                  &
                                                 __RC__  )
#endif

!
! !EXTERNAL STATE:
!
#if defined( MODEL_GEOS )
#   include "GEOSCHEMCHEM_ExportSpec___.h"
#else
#   include "GCHPchem_ExportSpec___.h"
#endif

    ! Read HISTORY config file and add exports for unique items
    CALL ESMF_ConfigGetAttribute( myState%myCF, HistoryConfigFile, &
                                  Label="HISTORY_CONFIG:",         &
                                  Default="HISTORY.rc", __RC__ )
    CALL HistoryExports_SetServices( MAPL_am_I_Root(), HistoryConfigFile, &
                                     GC, HistoryConfig, __RC__ )

!EOP
!BOC

#if defined( MODEL_GEOS )
    ! Add provider services, if any (AERO, RATS, Analysis Ox)
    CALL GEOS_CheckRATSandOx( am_I_Root, GC, __RC__ )

    ! Analysis options
    CALL GEOS_AnaInit( am_I_Root, GC, myState%myCF, ANAPHASE, __RC__ )

#endif

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
#if defined( MODEL_GEOS )
    USE GEOS_INTERFACE,       ONLY : GEOS_AddSpecInfoForMoist
    USE GEOS_AeroCoupler,     ONLY : GEOS_AeroInit
    USE GEOS_CarbonInterface, ONLY : GEOS_CarbonInit
!    USE TENDENCIES_MOD, ONLY : Tend_CreateClass
!    USE TENDENCIES_MOD, ONLY : Tend_Add
#endif
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

#if defined( MODEL_GEOS )
    ! Does GEOS-Chem restart file exist?
    ! Before broadcasting, we check if there is an import restart file for
    ! GEOS-Chem. This variable is then passed to Input_Opt after
    ! initialization (and CPU broadcasting) of all GEOS-Chem variables.
    LOGICAL                      :: haveImpRst

    TYPE(MAPL_MetaComp), POINTER :: STATE

    ! To read various options
    INTEGER                      :: DoIt
    REAL                         :: Val, OzPause

    LOGICAL                      :: DynFriend, IsPresent, FRIENDLY
    REAL, POINTER                :: Ptr3D(:,:,:) => NULL()

#else
    INTEGER                      :: N, trcID
    TYPE(MAPL_MetaComp), POINTER :: STATE => NULL()
    REAL(ESMF_KIND_R8), POINTER  :: Ptr3D(:,:,:) => NULL()
#endif

    ! Internal run alarms
    type(GC_run_alarms), pointer :: GC_alarms
    type(GCRA_wrap)              :: GC_alarm_wrapper
    TYPE(ESMF_Time)              :: startTime      ! Simulation start time
    TYPE(ESMF_Time)              :: currTime
    TYPE(ESMF_Time)              :: ringTime
    type(ESMF_TimeInterval)      :: tsRad_TI
    type(ESMF_TimeInterval)      :: tsChem_TI
    type (ESMF_Calendar)         :: CAL
    INTEGER                      :: yyyy, mm, dd   ! Year, month, day
    INTEGER                      :: h,    m,  s    ! Hour, minute, seconds
    INTEGER                      :: doy

    INTEGER                     :: IL_WORLD, JL_WORLD    ! # lower indices in global grid
    INTEGER                     :: IU_WORLD, JU_WORLD    ! # upper indices in global grid


    __Iam__('Initialize_')

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

#ifdef ADJOINT
    CALL MAPL_GenericStateClockAdd( GC, name='--AdjointCheckpoint', __RC__ )
#endif

    ! Get Internal state.
    CALL MAPL_Get ( STATE, INTERNAL_ESMF_STATE=INTSTATE, __RC__ )

    ! Initialize GEOS-Chem Input_Opt fields to zeros or equivalent
    CALL Set_Input_Opt( MAPL_am_I_Root(), Input_Opt, RC )
    _ASSERT(RC==GC_SUCCESS, 'Error calling Set_Input_Opt')

    ! Root CPU?
    am_I_Root = MAPL_am_I_Root()

    ! Get various parameters from the ESMF/MAPL framework
#if defined( MODEL_GEOS )
    ! Note: variable haveImpRst must not yet be written into Input_Opt yet
    ! since the variables of Input_Opt may be 'erased' during initialization
    ! of GEOS-Chem.
#endif
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
#if defined( MODEL_GEOS )
		   haveImpRst  = haveImpRst,  &  ! Does import restart exist?
#endif
                   __RC__                      )

    ! Set MPI values in Input_Opt
    Input_Opt%thisCPU = myPet
    Input_Opt%MPIComm = mpiComm
    Input_Opt%numCPUs = NPES
    Input_Opt%isMPI   = .true.
    if ( MAPL_am_I_Root() ) Input_Opt%amIRoot = .true.

#if defined( MODEL_GEOS )
    Input_Opt%haveImpRst = haveImpRst
#endif

    ! MSL - shift from 0 - 360 to -180 - 180 degree grid
    where (lonCtr .gt. MAPL_PI ) lonCtr = lonCtr - 2*MAPL_PI

#if !defined( MODEL_GEOS )
    ! Get the memory debug level
    call ESMF_ConfigGetAttribute(GeosCF, MemDebugLevel, &
                                 Label="MEMORY_DEBUG_LEVEL:" , RC=STATUS)
    _VERIFY(STATUS)
#endif

    !=======================================================================
    ! Save values from the resource file (GCHP.rc for GCHP)
    !=======================================================================

    ! # of run phases
    CALL ESMF_ConfigGetAttribute( GeosCF, NPHASE,                   &
                                  Default = 2,                      &
                                  Label   = "RUN_PHASES:",          &
                                  __RC__                           )
    _ASSERT(NPHASE==1.OR.NPHASE==2,'Error calling ESMF_ConfigGetAttribute on RUN_PHASES')

#if defined( MODEL_GEOS )
    ! Top stratospheric level
    CALL ESMF_ConfigGetAttribute( GeosCF, value_LLSTRAT,            &
                                  Default = LM,                     &
                                  Label   = "LLSTRAT:",             &
                                  __RC__                           )
    IF ( Input_Opt%AmIRoot ) THEN
       WRITE(*,*) 'GCC: top strat. level (LLSTRAT) is set to ',     &
                  value_LLSTRAT
    ENDIF

    ! FAST-JX settings: number of levels, number of EXTRAL iterations,
    ! print error if EXTRAL fails?
    ! LLFASTJX: default is 1201 for LM=132, 601 otherwise
    IF ( LM == 132 ) THEN
       I = 1201
    ELSE
       I = 601
    ENDIF
    CALL ESMF_ConfigGetAttribute( GeosCF, Input_Opt%LLFASTJX,     &
                                  Default = I,                    &
                                  Label   = "LLFASTJX:",          &
                                  __RC__                          )

    ! FJX_EXTRAL_ITERMAX: default is 5 for LM=132, 1 otherwise
    IF ( LM == 132 ) THEN
       I = 5
    ELSE
       I = 1
    ENDIF
    CALL ESMF_ConfigGetAttribute( GeosCF, Input_Opt%FJX_EXTRAL_ITERMAX, &
                                  Default = I,                          &
                                  Label   = "FJX_EXTRAL_ITERMAX:",      &
                                  __RC__                                )

    ! FJX_EXTRAL_ERR: default is 1
    I = 1
    CALL ESMF_ConfigGetAttribute( GeosCF, DoIt,   &
                                  Default = I,                        &
                                  Label   = "FJX_EXTRAL_ERR:",        &
                                  __RC__                              )
    Input_Opt%FJX_EXTRAL_ERR = ( DoIt == 1 )
    IF ( Input_Opt%AmIRoot ) THEN
       WRITE(*,*) 'GCC: Fast-JX settings:'
       WRITE(*,*) 'Number of FAST-JX levels      : ',Input_Opt%LLFASTJX
       WRITE(*,*) 'Max. no. of EXTRAL iterations : ',Input_Opt%FJX_EXTRAL_ITERMAX
       WRITE(*,*) 'Show EXTRAL overflow error    : ',Input_Opt%FJX_EXTRAL_ERR
    ENDIF

    ! Stop KPP if integration fails twice
    CALL ESMF_ConfigGetAttribute( GeosCF, DoIt,   &
                                  Default = 1,                       &
                                  Label   = "KPP_STOP_IF_FAIL:",     &
                                  __RC__                              )
    Input_Opt%KppStop = ( DoIt == 1 )
    IF ( Input_Opt%AmIRoot ) THEN
       WRITE(*,*) 'Stop KPP if integration fails: ',Input_Opt%KppStop
    ENDIF

    ! Turn off three heterogenous reactions in stratosphere
    CALL ESMF_ConfigGetAttribute( GeosCF, DoIt, Default = 0, &
                                  Label = "TurnOffHetRates:", __RC__ )
    Input_Opt%TurnOffHetRates = ( DoIt == 1 )
    IF ( Input_Opt%AmIRoot ) THEN
       WRITE(*,*) 'Disable selected het. reactions in stratosphere: ', &
                  Input_Opt%TurnOffHetRates
    ENDIF
#endif

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
    State_Grid%MaxTropLev  = 40            ! # trop. levels
#if defined( MODEL_GEOS )
    State_Grid%MaxStratLev = value_LLSTRAT ! # strat. levels
#else
    State_Grid%MaxStratLev = 59            ! # strat. levels
#endif

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
#if !defined( MODEL_GEOS )
                          GC        = GC,         & ! Ref to this gridded comp
                          EXPORT    = EXPORT,     & ! Export state object
#endif
                          Input_Opt = Input_Opt,  & ! Input Options obj
                          State_Chm = State_Chm,  & ! Chemistry State obj
                          State_Diag= State_Diag, & ! Diagnostics State obj
                          State_Grid= State_Grid, & ! Grid State obj
                          State_Met = State_Met,  & ! Meteorology State obj
                          HcoConfig = HcoConfig,  & ! HEMCO config obj
                          HistoryConfig = HistoryConfig, & ! History Config Obj
                          __RC__                 )

#if defined( MODEL_GEOS )
    !=======================================================================
    ! If GEOS-Chem is the AERO provider, initialize the AERO bundle here.
    ! All GEOS-Chem tracers possibly being added to the AERO bundle are
    ! listed at the beginning of the module. Here, we see which ones of
    ! those are effectively defined and create a field in the bundle for
    ! them. The AERO names are given the names listed at the beginning of
    ! the module.
    ! GEOS-Chem tracers are in mol/mol, whereas the AERO bundle holds
    ! data in kg/kg. We therefore need to copy the data so that we can
    ! change units independently.
    !=======================================================================
    IF ( DoAERO ) THEN
       CALL GEOS_AeroInit( GC, MaplCF, INTSTATE, EXPORT, Grid, __RC__ )
    ENDIF

#endif

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
#if defined( MODEL_GEOS )
       IF ( TRIM(Int2Spc(I)%Name) == 'POX' ) Int2Spc(I)%Name = 'POx'
       IF ( TRIM(Int2Spc(I)%Name) == 'LOX' ) Int2Spc(I)%Name = 'LOx'
#endif

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
#if defined( MODEL_GEOS )
          CYCLE
          _ASSERT(.FALSE.,'Error finding internal state variable')
#endif
       ENDIF

#if defined( MODEL_GEOS )
       ! Check friendliness of field: the field must not be friendly to
       ! moist and/or turbulence if the corresponding GEOS-Chem switches
       ! are turned on!
       DynFriend=.FALSE.
       CALL ESMF_AttributeGet( GcFld, NAME="FriendlyToDYNAMICS", &
                               isPresent=isPresent, RC=STATUS )
       IF ( isPresent ) THEN
          CALL ESMF_AttributeGet( GcFld, NAME="FriendlyToDYNAMICS", &
                               VALUE=DynFriend, RC=STATUS )
       ENDIF
       IF ( DynFriend ) THEN
          ! Check for friendliness to convection: only if GEOS-Chem convection
          ! is enabled
          FRIENDLY=.FALSE.
          CALL ESMF_AttributeGet( GcFld, NAME="FriendlyToMOIST", &
                                  isPresent=isPresent, RC=STATUS )
          IF ( isPresent ) THEN
             CALL ESMF_AttributeGet( GcFld, NAME="FriendlyToMOIST", &
                                     VALUE=FRIENDLY, RC=STATUS )
          ENDIF
          IF ( FRIENDLY .eqv. Input_Opt%LCONV ) THEN
             WRITE(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
             WRITE(*,*) 'GEOS-Chem convection and MOIST friendly on/off ', &
                        'at the same time', FRIENDLY, Input_Opt%LCONV,     &
                        TRIM(Int2Spc(I)%Name)
             WRITE(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
             _ASSERT(.FALSE.,'MOIST friendly error')
          ENDIF

          ! Check for friendliness to turbulence: only if GEOS-Chem turbulence
          ! is enabled
          IF ( Input_Opt%LTURB ) THEN
             FRIENDLY=.FALSE.
             CALL ESMF_AttributeGet( GcFld, NAME="FriendlyToTURBULENCE", &
                                     isPresent=isPresent, RC=STATUS )
             IF ( isPresent ) THEN
                CALL ESMF_AttributeGet( GcFld, NAME="FriendlyToTURBULENCE", &
                                        VALUE=FRIENDLY, RC=STATUS )
             ENDIF
             IF ( FRIENDLY ) THEN
                WRITE(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
                WRITE(*,*) 'GEOS-Chem turbulence is turned on, but tracer is ',&
                           'also friendly to TURBULENCE. Cannot do both: ',    &
                            TRIM(Int2Spc(I)%Name)
                WRITE(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
                _ASSERT(.FALSE.,'Error in Friendly settings')
             ENDIF
          ENDIF
       ENDIF
#endif

       ! Get pointer to field
       CALL ESMF_FieldGet( GcFld, 0, Ptr3D, __RC__ )
       Int2Spc(I)%Internal => Ptr3D
       Ptr3D => NULL()
       SpcInfo => NULL()

    ENDDO

#ifdef ADJOINT
    if (Input_Opt%is_Adjoint) THEN
       ! Now do the same for adjoint variables
       ALLOCATE( Int2Adj(nFlds), STAT=STATUS )
       _ASSERT(STATUS==0,'informative message here')

       ! Do for every tracer in State_Chm
       DO I = 1, nFlds

          ! Get info about this species from the species database
          N = State_Chm%Map_Advect(I)
          ThisSpc => State_Chm%SpcData(N)%Info

          ! Pass tracer name
          Int2Adj(I)%Name = TRIM(ThisSpc%Name)

          ! Get tracer ID
          Int2Adj(I)%ID = IND_( TRIM(Int2Spc(I)%Name) )

          ! If tracer ID is not valid, make sure all vars are at least defined.
          IF ( Int2Spc(I)%ID <= 0 ) THEN
             Int2Spc(I)%Internal => NULL()
             CYCLE
          ENDIF

          ! Get internal state field
          CALL ESMF_StateGet( INTSTATE, TRIM(SPFX) // TRIM(Int2Spc(I)%Name) // '_ADJ', &
               GcFld, RC=STATUS )

          ! This is mostly for testing
          IF ( STATUS /= ESMF_SUCCESS ) THEN
             IF( am_I_Root ) THEN
                WRITE(*,*) 'Cannot find in internal state: ', TRIM(SPFX) &
                     //TRIM(Int2Spc(I)%Name)//'_ADJ',I
             ENDIF
             _ASSERT(.FALSE.,'informative message here')
          ENDIF

          ! Get pointer to field
          CALL ESMF_FieldGet( GcFld, 0, Ptr3D, __RC__ )
          Int2Adj(I)%Internal => Ptr3D

          ! Free pointers
          Ptr3D => NULL()
          ThisSpc  => NULL()

       ENDDO
    ENDIF
#endif

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

    ! Establish the internal alarms for GEOS-Chem
    allocate(GC_alarms,stat=status)
    _ASSERT(rc==0,'Could not allocate GC alarms')
    GC_alarm_wrapper%ptr => GC_alarms

    call ESMF_UserCompSetInternalState(GC,'gcchem_internal_alarms',GC_alarm_wrapper,status)
    _ASSERT(status==0,'Could not get GEOS-Chem internal alarms')

    ! Get information about/from the clock
    CALL ESMF_ClockGet( Clock,                    &
                        startTime    = startTime, &
                        currTime     = currTime,  &
                        calendar     = cal,       &
                        __RC__ )

    ! Set up the radiation alarm
    ! Must ring once per tsRad
    call ESMF_TimeIntervalSet(tsRad_TI, S=nint(tsRad), calendar=cal, RC=STATUS)
    _ASSERT(STATUS==0,'Could not set radiation alarm time interval')

    ! Initially, just set the ring time to be midnight on the starting day
    call ESMF_TimeGet( startTime, YY = yyyy, MM = mm, DD = dd, H=h, M=m, S=s, rc = STATUS )
    _ASSERT(STATUS==0,'Could not extract start time information')
    call ESMF_TimeSet( ringTime,  YY = yyyy, MM = mm, DD = dd, H=0, M=0, S=0, rc = STATUS )
    _ASSERT(STATUS==0,'Could not set initial radiation alarm ring time')

    ! NOTE: RRTMG is run AFTER chemistry has completed. So we actually want the
    ! alarm to go off on the chemistry timestep immediately before the target
    ! output time.
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

#if defined( MODEL_GEOS )
    !=======================================================================
    ! Read GEOSCHEMchem settings
    !=======================================================================
    IF ( am_I_Root ) THEN
       WRITE(*,*) TRIM(Iam), ': options from GEOSCHEMchem_GridComp.rc:'
    ENDIF

    ! Apply correction term to large-scale precipitation that should be
    ! convective?
    CALL ESMF_ConfigGetAttribute( GeosCF, Input_Opt%WETD_CONV_SCAL,          &
                                  Label   = "Convective_precip_correction:", &
                                  Default = 1.0d0,                           &
                                  __RC__                                      )
    IF ( Input_Opt%WETD_CONV_SCAL < 0d0 ) Input_Opt%WETD_CONV_SCAL = 1.0d0
    Input_Opt%WETD_CONV_SCAL = MAX(MIN(1.0d0,Input_Opt%WETD_CONV_SCAL),0.0d0)
    IF ( am_I_Root ) THEN
       WRITE(*,*)   &
       '- Convective precip correction (0=no washout, 1=no correction): ', &
       Input_Opt%WETD_CONV_SCAL
    ENDIF

    ! Use GMI O3 P/L
    CALL ESMF_ConfigGetAttribute( GeosCF, DoIt, Label = "Use_GMI_O3_PL:", &
                                  Default = 0, __RC__ )
    Input_Opt%LSYNOZ = .FALSE.
    IF ( DoIt == 1 ) THEN
       Input_Opt%LGMIOZ = .TRUE.
       Input_Opt%LLINOZ = .FALSE.
    ELSE
       Input_Opt%LGMIOZ = .FALSE.
       !!!Input_Opt%LSYNOZ = .NOT. Input_Opt%LLINOZ
    ENDIF
    IF ( am_I_Root ) THEN
       WRITE(*,*) '- Use GMIOZ: ', Input_Opt%LGMIOZ
       WRITE(*,*) '- Use LINOZ: ', Input_Opt%LLINOZ
       WRITE(*,*) '- Use SYNOZ: ', Input_Opt%LSYNOZ
    ENDIF

    ! Tropopause options
    CALL ESMF_ConfigGetAttribute( GeosCF, DoIt,                    &
                                  Label = "Cap_polar_tropopause:", &
                                  Default = 1, __RC__ )
    Input_Opt%LCAPTROP = ( DoIt == 1 )
    IF ( am_I_Root ) THEN
       WRITE(*,*) '- Cap polar tropopause: ', Input_Opt%LCAPTROP
    ENDIF

    !CALL ESMF_ConfigGetAttribute( GeosCF, Input_Opt%NOx_sensitivity,   &
    !                              Label   = "NOx_sensitivity_factor:", &
    !                              Default = -999.0d0,                  &
    !                              __RC__                                )
    !IF ( am_I_Root ) THEN
    !   WRITE(*,*) '- Use NOx sensitivity factor: ', Input_Opt%NOx_sensitivity
    !ENDIF

    ! Get internal state from external data
    CALL ESMF_ConfigGetAttribute( GeosCF, DoIt, Label = "INIT_SPC_FROM_FILE:", &
                                  Default = 0, __RC__ )
    InitFromFile = ( DoIt == 1 )

    ! Skip GCC during replay predictor step
    CALL ESMF_ConfigGetAttribute( GeosCF, DoIt, Label = "SkipReplayGCC:", &
                                  Default = 0, __RC__ )
    SkipReplayGCC = ( DoIt == 1 )

    ! Always set stratospheric H2O
    CALL ESMF_ConfigGetAttribute( GeosCF, DoIt, &
          Label="Prescribe_strat_H2O:", Default=0, __RC__ )
    Input_Opt%AlwaysSetH2O = ( DoIt == 1 )
    IF ( am_I_Root ) THEN
       WRITE(*,*) '- Prescribe H2O in stratosphere: ', Input_Opt%AlwaysSetH2O
    ENDIF

    ! Compute vertical updraft velocity from online values
    CALL ESMF_ConfigGetAttribute( GeosCF, DoIt, &
          Label="Online_VUD:", Default=0, __RC__ )
    Input_Opt%UseOnlineVUD = ( DoIt == 1 )
    IF ( am_I_Root ) THEN
       WRITE(*,*) '- Compute VUD online: ', Input_Opt%UseOnlineVUD
    ENDIF

    ! Turn on Family Transport
    CALL ESMF_ConfigGetAttribute( GeosCF, DoIt, &
          Label="Bry_Cly_Family_Transport:", Default=1, __RC__ )
    SELECT CASE ( DoIt )
       CASE ( 1 )
          CALL Init_GCC_Chem_Groups()
          IF ( am_I_Root ) WRITE(*,*) 'GCC: Bry and Cly family transport enabled'
       CASE DEFAULT
          IF ( am_I_Root ) WRITE(*,*) 'GCC: Bry and Cly family transport disabled'
    END SELECT

    ! Add Henry law constants and scavenging coefficients to internal state.
    ! These are needed by MOIST for wet scavenging (if this is enabled).
    CALL GEOS_AddSpecInfoForMoist ( am_I_Root, GC, GeosCF, Input_Opt, State_Chm, __RC__ )

    ! Initialize carbon coupling / CO production from CO2 photolysis (if used) 
    CALL GEOS_CarbonInit( GC, GeosCF, State_Chm, State_Grid, __RC__ ) 

    !=======================================================================
    ! All done
    !=======================================================================
#endif

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
#if defined( MODEL_GEOS )
    USE GEOS_Analysis,           ONLY : GEOS_AnaRun
    USE GEOS_Interface,          ONLY : MetVars_For_Lightning_Run, &
                                        GEOS_Diagnostics,          &
                                        GEOS_CalcTotOzone,         &
                                        GEOS_InitFromFile,         &
                                        GEOS_RATSandOxDiags,       &
                                        GEOS_PreRunChecks
    USE GEOS_AeroCoupler,        ONLY : GEOS_FillAeroBundle
    USE GEOS_CarbonInterface,    ONLY : GEOS_CarbonGetConc,        &
                                        GEOS_CarbonRunPhoto
#endif

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
#if !defined( MODEL_GEOS )
    REAL(ESMF_KIND_R8),  POINTER :: PLE(:,:,:)     => NULL() ! INTERNAL: PEDGE
#endif

    ! Initialize variables used for reading Olson and MODIS LAI imports
    INTEGER            :: TT, VV, landTypeInt
    CHARACTER(len=64)  :: landTypeStr, varName, importName

#if defined( MODEL_GEOS )
    ! GEOS-5 only local variables
    INTEGER                      :: LB            ! Loop indices

    ! Some checks for replay runs
    LOGICAL                      :: FIRSTREWIND
    LOGICAL                      :: AFTERREWIND
    LOGICAL                      :: IsFirst
    INTEGER, SAVE                :: pymd = 0         ! previous date
    INTEGER, SAVE                :: phms = 0         ! previous time
    INTEGER, SAVE                :: nnRewind = 0

    ! For skipping GCC during predictor step                                                  
    TYPE(ESMF_Alarm)             :: PredictorAlarm  ! skip GCC during replay
    LOGICAL                      :: PredictorActive ! skip GCC during replay

#else
    ! GCHP only local variables
    INTEGER                      :: trcID, RST
    REAL                         :: COEFF
    CHARACTER(LEN=ESMF_MAXSTR)   :: trcNAME,hcoNAME
    TYPE(ESMF_Field      )       :: trcFIELD
    TYPE(ESMF_FieldBundle)       :: trcBUNDLE
    REAL              , POINTER  :: fPtrArray(:,:,:)
    REAL(ESMF_KIND_R8), POINTER  :: fPtrVal, fPtr1D(:)
    INTEGER                      :: z_lb, z_ub

#endif

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
#ifdef ADJOINT
    LOGICAL                      :: isStartTime
    REAL(ESMF_KIND_r8), POINTER  :: CostFuncMask(:,:,:) => NULL()
#endif

    __Iam__('Run_')

    !=======================================================================
    ! Run starts here
    !=======================================================================

    ! Are we on the root PET?
    am_I_Root = MAPL_Am_I_Root()

    ! Set up traceback info
    CALL ESMF_GridCompGet( GC, name=compName, __RC__ )

    ! Identify this routine to MAPL
    Iam = TRIM(compName)//'::Run_'

#if !defined( MODEL_GEOS )
    ! Get the VM for optional memory prints (level >= 2)
    !-----------------------------------
    if ( MemDebugLevel > 0 ) THEN
       call ESMF_VmGetCurrent(VM, RC=STATUS)
       _VERIFY(STATUS)
    endif
#endif

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
#if defined( MODEL_GEOS )
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
#endif

#ifdef ADJOINT
    if (Input_Opt%is_adjoint .and. first) THEN
       ! the forward model doesn't actually trigger on the final
       ! timestep, so we should skip the first one
       IsRunTime = .false.
    end if
#endif
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
#if defined( MODEL_GEOS )
#      include "GEOSCHEMCHEM_GetPointer___.h"
#else
#      include "GCHPchem_GetPointer___.h"

       !IF ( IsCTM ) THEN
       call MAPL_GetPointer ( IMPORT, PLE,      'PLE',     __RC__ )
       !ENDIF

       ! Pass IMPORT/EXPORT object to HEMCO state object
       !CALL GetHcoState( HcoState )
       _ASSERT(ASSOCIATED(HcoState),'HcoState is not associated')
       HcoState%GRIDCOMP => GC
       HcoState%IMPORT   => IMPORT
       HcoState%EXPORT   => EXPORT
       !HcoState => NULL()

#endif
    ENDIF

#if defined( MODEL_GEOS )
! GCHP ends the (if FIRST) block and then links HEMCO state to GC objects:
!    ENDIF

   ! Link HEMCO state to gridcomp objects
   ASSERT_(ASSOCIATED(HcoState))
   HcoState%GRIDCOMP => GC
   HcoState%IMPORT   => IMPORT
   HcoState%EXPORT   => EXPORT
#endif
#ifdef ADJOINT
       call MAPL_GetPointer( IMPORT, CostFuncMask, &
            'CFN_MASK', notFoundOK=.TRUE.,         &
            __RC__ )
       if (MAPL_Am_I_Root() .and. .not. ASSOCIATED(CostFuncMask)) &
            WRITE(*,*) ' No CFN_MASK import variable found'
#endif

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

#if !defined( MODEL_GEOS )
       ! MSL - shift from 0 - 360 to -180 - 180 degree grid
       where (lonCtr .gt. MAPL_PI ) lonCtr = lonCtr - 2*MAPL_PI
#endif

#if defined( MODEL_GEOS )
       ! Check if this time is before the datetime of the prev timestep, e.g.
       ! if this is after a clock rewind
       AFTERREWIND = .FALSE.
       FIRSTREWIND = .FALSE.
       IF ( nymd < pymd ) THEN
          AFTERREWIND = .TRUE.
       ELSEIF ( (nymd == pymd) .AND. (nhms < phms) ) THEN
          AFTERREWIND = .TRUE.
       ENDIF

       ! If this is after a rewind, check if it's the first rewind. In this
       ! case, we need to re-do some first-time assignments to make sure that
       ! we reset all variables to the initial state!
       IF ( AFTERREWIND ) THEN
          nnRewind = nnRewind + 1
          IF ( nnRewind == 1 ) FIRSTREWIND = .TRUE.
       ENDIF
#endif

       ! Pass grid area [m2] obtained from dynamics component to State_Grid
       State_Grid%Area_M2 = AREA
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

#ifdef ADJOINT
          IF (IsRunTime) THEN
             IF (Input_opt%IS_ADJOINT) THEN
                call WRITE_PARALLEL('  Resetting state from checkpoint file')
                ! call MAPL_GenericRefresh(GC, Import, Export, Clock, RC)
                call Adjoint_StateRefresh( GC, IMPORT, EXPORT, CLOCK, RC )
                ! Loop over all species and get info from spc db
                DO N = 1, State_Chm%nSpecies
                   ThisSpc => State_Chm%SpcData(N)%Info
                   !IF (ThisSpc%Is_Advected) CYCLE
                   IF ( TRIM(ThisSpc%Name) == '' ) CYCLE
                   IND = IND_( TRIM(ThisSpc%Name ) )
                   IF ( IND < 0 ) CYCLE
                   ! Get data from internal state and copy to species array
                   CALL MAPL_GetPointer( INTERNAL, Ptr3D_R8, TRIM(SPFX) //          &
                        TRIM(ThisSpc%Name), notFoundOK=.TRUE.,     &
                        __RC__ )
                   State_Chm%Species(IND)%Conc(:,:,:) = &
                                 Ptr3D_R8(:,:,State_Grid%NZ:1:-1)
                   if ( MAPL_am_I_Root()) WRITE(*,*)                                &
                        'Initialized species from INTERNAL state: ', TRIM(ThisSpc%Name)

                enddo
             ELSE
                call WRITE_PARALLEL('  Recording state to checkpoint file')
                call Adjoint_StateRecord( GC, IMPORT, EXPORT, CLOCK, RC )
                ! call WRITE_PARALLEL('  Done recording state to checkpoint files')
             ENDIF
          ENDIF
#endif

!       !=======================================================================
!       ! pre-Run method array assignments. This passes the tracer arrays from
!       ! the internal state to State_Chm. On the first call, it also fills the
!       ! internal species arrays in State_Chm with the values read from the
!       ! restart file (and stored in the internal state).
!       !=======================================================================
       CALL MAPL_TimerOn(STATE, "CP_BFRE")
#include "Includes_Before_Run.H"
       CALL MAPL_TimerOff(STATE, "CP_BFRE")

#if defined( MODEL_GCHPCTM )
       !=======================================================================
       ! Point GEOS-Chem species concentration arrays to internal state
       !=======================================================================
       DO I = 1, SIZE(Int2Spc,1)
          IF ( Int2Spc(I)%ID <= 0 ) CYCLE
          State_Chm%Species(Int2Spc(I)%ID)%Conc => Int2Spc(I)%Internal(:,:,State_Grid%NZ:1:-1)
       ENDDO

#ifdef ADJOINT
      IF (Input_Opt%Is_Adjoint) THEN
         DO I = 1, SIZE(Int2Adj,1)
            IF ( Int2Adj(I)%ID <= 0 ) CYCLE
            State_Chm%SpeciesAdj(:,:,:,Int2Adj(I)%ID) = Int2Adj(I)%Internal
         ENDDO

         ! Flip in the vertical
         State_Chm%SpeciesAdj = State_Chm%SpeciesAdj( :, :, State_Grid%NZ:1:-1, : )
      ENDIF
#endif
       !=======================================================================
       ! On first call, also need to initialize the species from restart file.
       ! Only need to do this for species that are not advected, i.e. species
       ! that are not tracers (all other species arrays will be filled with
       ! tracer values anyways!).
       ! We only need to do this on the first call because afterwards, species
       ! array already contains values from previous chemistry time step
       ! (advected species will be updated with tracers)
       ! ckeller, 10/27/2014
       !=======================================================================
#ifdef ADJOINT
       IF ( FIRST .or. Input_Opt%IS_ADJOINT) THEN
#else
       IF ( FIRST ) THEN
#endif

          ! Get Generic State
          call MAPL_GetObjectFromGC ( GC, STATE, RC=STATUS)
          _VERIFY(STATUS)
          ! Get Internal state
          CALL MAPL_Get ( STATE, INTERNAL_ESMF_STATE=INTERNAL, __RC__ )

          ! Loop over all species and get info from spc db
          DO N = 1, State_Chm%nSpecies
             ThisSpc => State_Chm%SpcData(N)%Info
             IF (ThisSpc%Is_Advected) CYCLE
             IF ( TRIM(ThisSpc%Name) == '' ) CYCLE
             IND = IND_( TRIM(ThisSpc%Name ) )
             IF ( IND < 0 ) CYCLE

             ! Get data from internal state and copy to species array
             CALL MAPL_GetPointer( INTERNAL, Ptr3D_R8, TRIM(SPFX) //          &
                                   TRIM(ThisSpc%Name), notFoundOK=.TRUE.,     &
                                   __RC__ )
             IF ( .NOT. ASSOCIATED(Ptr3D_R8) ) THEN
                IF ( MAPL_am_I_Root()) WRITE(*,*)                             &
                   'Could not find species in INTERNAL state - will be ' //   &
                   'initialized to zero: ', TRIM(SPFX), TRIM(ThisSpc%Name)
                State_Chm%Species(IND)%Conc(:,:,:) = 1d-26
                CYCLE
             ENDIF

             ! Determine if species in restart file
             CALL ESMF_StateGet( INTERNAL, TRIM(SPFX) // TRIM(ThisSpc%Name),  &
                  trcFIELD, RC=RC )
             CALL ESMF_AttributeGet( trcFIELD, NAME="RESTART",                &
                  VALUE=RST, RC=STATUS )

             ! Set spc conc to background value if rst skipped or var not there
             IF ( RC /= ESMF_SUCCESS .OR. RST == MAPL_RestartBootstrap .OR.   &
                      RST == MAPL_RestartSkipInitial ) THEN
                DO L = 1, State_Grid%NZ
                DO J = 1, State_Grid%NY
                DO I = 1, State_Grid%NX
                   IF ( L > State_Grid%MaxChemLev .AND. &
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
                   '   WARNING: using background values from species database'
                ENDIF
             ENDIF
             ThisSpc => NULL()
          ENDDO
       ENDIF

       !=======================================================================
       ! On first call, initialize certain State_Chm and State_Met arrays from
       ! imports if they are found (ewl, 12/13/18)
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
             State_Chm%KPPHvalue(:,:,1:State_Grid%MaxChemLev) =       &
            Ptr3d_R8(:,:,State_Grid%NZ:State_Grid%NZ-State_Grid%MaxChemLev+1:-1)
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

          CALL MAPL_GetPointer( INTSTATE, Ptr3d_R8, 'DELP_DRY', notFoundOK=.TRUE., __RC__ )
          IF ( ASSOCIATED(Ptr3d_R8) .AND. ASSOCIATED(State_Met%DELP_DRY) ) THEN
             State_Met%DELP_DRY(:,:,1:State_Grid%NZ) =       &
                                  Ptr3d_R8(:,:,State_Grid%NZ:1:-1)
          ENDIF
          Ptr3d_R8 => NULL()

          CALL MAPL_GetPointer( INTSTATE, Ptr3d_R8, 'BXHEIGHT', notFoundOK=.TRUE., __RC__ )
          IF ( ASSOCIATED(Ptr3d_R8) .AND. ASSOCIATED(State_Met%BXHEIGHT) ) THEN
             State_Met%BXHEIGHT(:,:,1:State_Grid%NZ) =       &
                                  Ptr3d_R8(:,:,State_Grid%NZ:1:-1)
          ENDIF
          Ptr3d_R8 => NULL()

          CALL MAPL_GetPointer( INTSTATE, Ptr2d_R8, 'TropLev', notFoundOK=.TRUE., __RC__ )
          IF ( ASSOCIATED(Ptr2d_R8) .AND. ASSOCIATED(State_Met%TropLev) ) THEN
             State_Met%TropLev = Ptr2d_R8
          ENDIF
          Ptr2d_R8 => NULL()
       ENDIF
#endif

#if defined( MODEL_GEOS )
       CALL MetVars_For_Lightning_Run( GC, Import=IMPORT, Export=EXPORT, &
             State_Met=State_Met, State_Grid=State_Grid, __RC__ )

       ! Import CO2 concentrations from external field (e.g., GOCART)
       ! Note: comment out for now and call this below immediately before doing
       ! CO2 photolysis. Reason for this is that CO2 is currently hardcoded to 
       ! 0.0 in fullchem to match the old SMVGEAR code (??)
!       CALL GEOS_CarbonGetConc( Import,    Input_Opt,  State_Chm,        &
!                                State_Met, State_Diag, State_Grid, __RC__ )

       ! Eventually initialize species concentrations from external field.
       IsFirst = ( FIRST .OR. FIRSTREWIND )
       IF ( InitFromFile ) THEN
          CALL GEOS_InitFromFile( GC, Import, INTSTATE, Export, GeosCF, &
                                  Input_Opt, State_Met, State_Chm, Q,  &
                                  PLE, GCCTROPP, IsFirst, __RC__ )
       ENDIF

       ! Initialize RATS and OX tendency diagnostics
       IF ( PHASE == ANAPHASE ) THEN
          CALL GEOS_RATSandOxDiags( GC, INTSTATE, Export, Input_Opt, State_Met, &
                                    State_Chm, State_Grid, Q, 1, tsChem, __RC__ ) 
       ENDIF

       ! Perform GEOS pre-run checks
       CALL GEOS_PreRunChecks( am_I_Root, Input_Opt, State_Met, State_Chm, &
                               GeosCF, IsFirst, __RC__ )
#endif

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
#if defined( MODEL_GEOS )
       IF ( FIRST .OR. FIRSTREWIND ) THEN
#else
       IF ( FIRST ) THEN
#endif

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
       ! Get total ozone column from GEOS-Chem export variable.
       ! Need to calculate from restart variables on first call!
       !=======================================================================
#if defined( MODEL_GEOS )
       IF ( PHASE /= 1 ) THEN
          CALL GEOS_CalcTotOzone( am_I_Root, State_Met, State_Chm, State_Diag, PLE, TROPP, __RC__ )
       ENDIF
#endif

       !=======================================================================
       ! Execute GEOS-Chem on multiple PETs
       !=======================================================================

       ! Fix negatives!
       ! These can be brought in as an artifact of convection.
#ifndef ADJOINT
       DO N = 1, State_Chm%nSpecies
          WHERE ( State_Chm%Species(N)%Conc < 0.0e0 )
             State_Chm%Species(N)%Conc = 1.0e-36
          END WHERE 
       ENDDO
#endif

       ! Execute GEOS-Chem if it's time to run it
       IF ( IsRunTime ) THEN

          ! This is mostly for testing
#if defined( MODEL_GEOS )
          IF ( FIRST .AND. Input_Opt%haveImpRst ) THEN
#else
          IF ( FIRST ) THEN
#endif
             IF ( am_I_Root ) THEN
                WRITE(*,*) ''
                WRITE(*,*) 'Doing warm GEOS-Chem restart'
                WRITE(*,*) ''
             ENDIF
          ENDIF

#if defined( MODEL_GEOS )
          ! Only if restart file exists...
          IF ( Input_Opt%haveImpRst ) THEN
#endif

#if !defined( MODEL_GEOS )
             ! Optional memory prints (level >= 2)
             if ( MemDebugLevel > 0 ) THEN
                call ESMF_VMBarrier(vm, RC=STATUS)
                _VERIFY(STATUS)
                call MAPL_MemUtilsWrite(VM, &
                  'Chem_GridCompMod, before chunk_run', RC=STATUS )
                _VERIFY(STATUS)
             endif
#endif

             CALL MAPL_TimerOn(STATE, "DO_CHEM")

#if !defined( MODEL_GEOS )
             ! NOTE: Second was not extracted previously; set to 0 for now
             second = 0
#endif

#ifdef ADJOINT
             !=======================================================================
             ! If this is an adjoint run, we need to check for the final (first)
             ! timestep and multiply the scaling factor adjoint by the initial concs
             !=======================================================================
             isStartTime = .false.
             IF (Input_Opt%IS_ADJOINT) THEN
                call ESMF_ClockGet(clock, currTime=currTime, startTime=stopTime,  __RC__ )
             else
                call ESMF_ClockGet(clock, currTime=currTime, stopTime=stopTime, __RC__ )
             Endif

                ! call ESMF_TimeIntervalSet(tsChemInt, s_r8=real(-tsChem, 8), __RC__ )
             ! this variable is set to zero but I'm leaving it in case I need this code later
                call ESMF_TimeIntervalSet(tsChemInt, s_r8=real(0, 8), __RC__ )

                call ESMF_TimeGet(currTime + tsChemInt, timeString=timestring1, __RC__ )
                call ESMF_TimeGet(stopTime, timeString=timestring2, __RC__ )

                if (memdebuglevel > 0 .and. am_I_Root) &
                     WRITE(*,*) '   Adjoint checking if ' // trim(timestring1) // ' == ' // trim(timestring2)

                if (currTime + tsChemInt == stopTime) THEN
                   isStartTime = .TRUE.
                ENDIF
#endif

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
#if defined( MODEL_GEOS )
                                  FrstRewind = FirstRewind,& ! First rewind?
#endif
#ifdef ADJOINT
                                  isStartTime = isStartTime, & !back to the first timestep in the reverse run?
#endif
                                  __RC__                  )  ! Success or fail?

             CALL MAPL_TimerOff(STATE, "DO_CHEM")

#if !defined( MODEL_GEOS )
             ! Optional memory prints (level >= 2)
             if ( MemDebugLevel > 0 ) THEN
                call ESMF_VMBarrier(vm, RC=STATUS)
                _VERIFY(STATUS)
                call MAPL_MemUtilsWrite(VM, &
                  'Chem_GridCompMod, after  chunk_run', RC=STATUS )
                _VERIFY(STATUS)
             endif
#endif

#if !defined( MODEL_GEOS )
             where( State_Met%HFLUX .eq. 0.) State_Met%HFLUX = 1e-5
#endif

#if defined( MODEL_GEOS )
          ! Restart file does not exist:
          ELSE
             IF ( am_I_Root ) THEN
                WRITE(*,*) ''
                WRITE(*,*) ' SKIP GEOS-CHEM PHASE ', Phase, &
                           ' BECAUSE IMPORT RESTART FILE IS MISSING'
                WRITE(*,*) ''
             ENDIF
          ENDIF
#endif

       ENDIF !IsRunTime

       !=======================================================================
       ! post-Run method array assignments. This copies the values back from
       ! the State_Chm tracer arrays to the internal state, so that they can
       ! be seen by other components (moist, turbulence, ...)
       !=======================================================================

#if defined( MODEL_GCHPCTM )
       CALL MAPL_TimerOn(STATE, "CP_AFTR")
#endif

#if defined( MODEL_GEOS )
       !=======================================================================
       ! GEOS post-run procedures
       !=======================================================================

       ! Import CO2 concentrations from external field (e.g., GOCART)
       ! Note: this call should be moved to the beginning of the routine once
       ! CO2 is not hardcoded to 0.0 anymore (in fullchem)
       CALL GEOS_CarbonGetConc( Import,    Input_Opt,  State_Chm,        &
                                State_Met, State_Diag, State_Grid, __RC__ )

       IF ( PHASE == CHEMPHASE ) THEN
          ! CO production from CO2 photolysis, using StratChem code
          CALL GEOS_CarbonRunPhoto( Input_Opt,  State_Chm,  State_Met, &
                                    State_Diag, State_Grid, __RC__      )
       ENDIF

       IF ( PHASE == ANAPHASE ) THEN
          ! Call GEOS analysis routine
          CALL GEOS_AnaRun( GC, Import, INTSTATE, Export, Clock, &
                            Input_Opt, State_Met, State_Chm, Q, PLE, TROPP, __RC__ )

          ! GEOS Diagnostics. This includes the 'default' GEOS-Chem diagnostics.
          CALL GEOS_Diagnostics( GC, IMPORT, EXPORT, Clock, Phase, Input_Opt, &
                                 State_Met, State_Chm, State_Diag, State_Grid, __RC__ )

          ! Fill RATS and OX diagnostics
          CALL GEOS_RATSandOxDiags( GC, INTSTATE, Export, Input_Opt, State_Met, &
                                    State_Chm, State_Grid, Q, 2, tsChem, __RC__ ) 
       ENDIF

       ! Update internal state fields 
       CALL MAPL_TimerOn(STATE, "CP_AFTR")
#      include "Includes_After_Run.H"
       CALL MAPL_TimerOff(STATE, "CP_AFTR")

       ! Connect to aerosols - experimental
       IF ( DoAERO ) THEN
          CALL GEOS_FillAeroBundle ( GC, EXPORT, State_Chm, State_Grid, Input_Opt, __RC__ )
       ENDIF

       ! Archive last active time steps
       pymd = nymd
       phms = nhms
#endif

#if defined( MODEL_GCHPCTM )
#ifdef ADJOINT
       IF (Input_Opt%Is_Adjoint) THEN
          State_Chm%SpeciesAdj = State_Chm%SpeciesAdj(:,:,State_Grid%NZ:1:-1,:)

          DO I = 1, SIZE(Int2Adj,1)
             WRITE(*,*) 'Copying adjoint ', Int2Adj(I)%ID, ' to ', I
             IF ( Int2Adj(I)%ID <= 0 ) CYCLE
             Int2Adj(I)%Internal = State_Chm%SpeciesAdj(:,:,:,Int2Adj(I)%ID)
          ENDDO
       ENDIF
#endif
       CALL MAPL_TimerOff(STATE, "CP_AFTR")

       ! Update non-species dynamic internal state arrays post-run
       ! every timestep, except for area which can be set first run only
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
          Ptr3d_R8(:,:,1:State_Grid%NZ-State_Grid%MaxChemLev) = 0.0
          Ptr3d_R8(:,:,State_Grid%NZ:State_Grid%NZ-State_Grid%MaxChemLev+1:-1)=&
             State_Chm%KPPHvalue(:,:,1:State_Grid%MaxChemLev)
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
#endif
       
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

#if defined( MODEL_GCHPCTM )
    !=======================================================================
    ! Nullify GEOS-Chem species concentration pointers
    !=======================================================================
    DO I = 1, SIZE(Int2Spc,1)
       IF ( Int2Spc(I)%ID <= 0 ) CYCLE
       IF ( ASSOCIATED( State_Chm%Species(Int2Spc(I)%ID)%Conc ) ) THEN
          State_Chm%Species(Int2Spc(I)%ID)%Conc => NULL()
       ENDIF
    ENDDO
#endif

    !=======================================================================
    ! All done
    !=======================================================================

    IF ( ALLOCATED( zenith     ) ) DEALLOCATE( zenith     )
    IF ( ALLOCATED( solar      ) ) DEALLOCATE( solar      )

    ! Stop timer
    ! ----------
    CALL MAPL_TimerOff(STATE, "TOTAL")

#if defined( MODEL_GEOS )
    ! The restart file should exist at least after the first full cycle,
    ! e.g. after phase 1 and phase 2 has been called once.
    IF ( FIRST .AND. Phase == 1 ) THEN
       Input_Opt%haveImpRst = Input_Opt%haveImpRst
    ELSE
       Input_Opt%haveImpRst = .TRUE.
    ENDIF
#endif

    ! Update first flags
    FIRST = .FALSE.

#if defined( MODEL_GEOS )
    ! Unlink HEMCO state from gridcomp objects
    HcoState%GRIDCOMP => NULL()
    HcoState%IMPORT   => NULL()
    HcoState%EXPORT   => NULL()
#endif

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
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Model                            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Finalize_
!
! !DESCRIPTION: Finalize_ is the finalize method of the GEOSCHEMchem gridded
!  component.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Finalize_( GC, Import, Export, Clock, RC )
!
! !USES:
!
    USE CMN_Size_Mod,          ONLY : NDUST
    USE Input_Opt_Mod,         ONLY : OptInput
#if !defined( MODEL_GEOS )
    USE Input_Opt_Mod,         ONLY : Cleanup_Input_Opt
#endif
    USE State_Chm_Mod,         ONLY : ChmState, Cleanup_State_Chm
    USE State_Diag_Mod,        ONLY : DgnState, Cleanup_State_Diag
    USE State_Grid_Mod,        ONLY : GrdState, Cleanup_State_Grid
    USE State_Met_Mod,         ONLY : MetState, Cleanup_State_Met
    USE HCO_Interface_GC_Mod,  ONLY : HCOI_GC_FINAL
#if defined( MODEL_GEOS )
    USE HCO_State_GC_Mod,      ONLY : HcoState
    USE GEOS_Analysis,         ONLY : GEOS_AnaFinal
#endif
!
! !INPUT/OUTPUT PARAMETERS:
!
#if defined( MODEL_GEOS )
    TYPE(ESMF_GridComp), INTENT(INOUT), TARGET :: GC       ! Ref. to this GC
    TYPE(ESMF_State),    INTENT(INOUT), TARGET :: Import   ! Import State
    TYPE(ESMF_State),    INTENT(INOUT), TARGET :: Export   ! Export State
#else
    TYPE(ESMF_GridComp), INTENT(INOUT)         :: GC       ! Ref. to this GC
    TYPE(ESMF_State),    INTENT(INOUT)         :: Import   ! Import State
    TYPE(ESMF_State),    INTENT(INOUT)         :: Export   ! Export State
#endif
    TYPE(ESMF_Clock),    INTENT(INOUT)         :: Clock    ! ESMF Clock object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,             INTENT(OUT)           :: RC       ! Success or failure?
!
! !REMARKS:
!  We call routine Extract_ to return various values (i.e. grid parameters,
!  start & end dates, PET information, etc.) from the ESMF/MAPL environment.
!  We then pass those to GEOS-Chem via routine GCHP_CHUNK_FINAL, which is
!  located in GEOS-Chem module ./GEOS-Chem/ESMF/gchp_chunk_mod.F90.
!
! !REVISION HISTORY:
!  01 Dec 2009 - A. Da Silva - Initial version
!  See https://github.com/geoschem/geos-chem for history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! LOCAL VARIABLES:
!
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
#if !defined( MODEL_GEOS )
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
#endif
#ifdef ADJOINT
    ! Finite difference test variables
    INTEGER                        :: IFD, JFD, LFD
    REAL*8                         :: CFN
    CHARACTER(len=ESMF_MAXSTR)     :: FD_SPEC
#endif


    __Iam__('Finalize_')

    !=======================================================================
    ! Initialization
    !=======================================================================

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

    !=======================================================================
    ! Finalize the Gridded Component
    !=======================================================================

#if defined( MODEL_GEOS )
    ! Link HEMCO state to gridcomp objects
    _ASSERT(ASSOCIATED(HcoState),'HcoState is not associated')
    HcoState%GRIDCOMP => GC
    HcoState%IMPORT   => IMPORT
    HcoState%EXPORT   => EXPORT
#endif

#ifdef ADJOINT
    IF (Input_Opt%IS_FD_SPOT_THIS_PET .and. .not. Input_Opt%IS_FD_GLOBAL) THEN
       FD_SPEC = transfer(state_chm%SpcData(Input_Opt%NFD)%Info%Name, FD_SPEC)
       IFD = Input_Opt%IFD
       JFD = Input_Opt%JFD
       LFD = Input_Opt%LFD
       NFD = Input_Opt%NFD
       ! print out the cost function
       WRITE(*,*) ' Computing final cost function'
       CFN = 0d0
       DO L = 1, State_Grid%NZ
       DO J = 1, State_Grid%NY
       DO I = 1, State_Grid%NX
          if (State_Chm%CostFuncMask(I,J,L) > 0d0) THEN
             WRITE (*, 1047) I, J, L, State_Chm%Species(NFD)%conc(I,J,L)
             CFN = CFN + State_Chm%Species(NFD)%Conc(I,J,L)
          endif
       ENDDO
        ENDDO
       ENDDO
       WRITE(*,'(a7, e22.10)') ' CFN = ', CFN
1047   FORMAT('  SPC(', i2, ', ', i2, ', ', i2, ') = ', e22.10)
    ENDIF
#endif

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

#if !defined( MODEL_GEOS )
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
#endif

    ! Free Int2Spc pointer
    IF ( ASSOCIATED(Int2Spc) ) THEN
       DO I=1,SIZE(Int2Spc,1)
          Int2Spc(I)%Internal => NULL()
       ENDDO
       DEALLOCATE(Int2Spc)
    ENDIF

#ifdef ADJOINT
    ! Free Int2Adj pointer
    IF ( ASSOCIATED(Int2Adj) ) THEN
       DO I=1,SIZE(Int2Adj,1)
          Int2Adj(I)%Internal => NULL()
       ENDDO
       DEALLOCATE(Int2Adj)
    ENDIF
#endif

    ! Deallocate the history interface between GC States and ESMF Exports
    CALL Destroy_HistoryConfig( am_I_Root, HistoryConfig, RC )

#if defined( MODEL_GEOS )
    ! Cleanup GEOS analysis module
    CALL GEOS_AnaFinal( __RC__ )
#endif

    ! Finalize MAPL Generic
    CALL MAPL_GenericFinalize( GC, Import, Export, Clock, __RC__ )

    !=======================================================================
    ! All done
    !=======================================================================

    ! Stop timers
    ! -----------
!    CALL MAPL_TimerOff(STATE, "FINALIZE")
!    CALL MAPL_TimerOff(STATE, "TOTAL")

    ! Successful return
    _RETURN(ESMF_SUCCESS)

  END SUBROUTINE Finalize_
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Model                            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Extract_
!
! !DESCRIPTION: GC routine extracts several common quantities from the
!  ESMF/MAPL environment so that they can be later passed down to the
!  grid-independent GEOS-Chem code.
!\\
!\\
! !INTERFACE:
!
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
#if defined( MODEL_GEOS )
                       haveImpRst,                                       &
#endif
                       RC )
!
! !INPUT PARAMETERS:
!
    TYPE(ESMF_Clock),    INTENT(IN)            :: Clock       ! ESMF clock obj
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ESMF_GridComp), INTENT(INOUT)         :: GC          ! GC grid comp
!
! !OUTPUT PARAMETERS:
!
    !----------------------------------
    ! ESMF and/or MAPL quantities
    !----------------------------------
    TYPE(ESMF_Grid),     INTENT(OUT), OPTIONAL :: Grid        ! ESMF Grid obj
    TYPE(ESMF_Config),   INTENT(OUT), OPTIONAL :: MaplCF      ! AGCM.rc
    TYPE(ESMF_Config),   INTENT(OUT), OPTIONAL :: GeosCF      ! GEOSCHEM*.rc
    INTEGER,             INTENT(OUT), OPTIONAL :: localPet    ! This PET
    INTEGER,             INTENT(OUT), OPTIONAL :: petCount    ! Total # of PETs
    INTEGER,             INTENT(OUT), OPTIONAL :: mpiComm     ! MPI Comm Handle

    !----------------------------------
    ! Local grid coordinates
    ! (defined on the current CPU)
    !----------------------------------
    INTEGER,             INTENT(OUT), OPTIONAL :: IM          ! Total # lons
    INTEGER,             INTENT(OUT), OPTIONAL :: JM          ! Total # lats
    INTEGER,             INTENT(OUT), OPTIONAL :: LM          ! Total # levs
    REAL(ESMF_KIND_R4),  POINTER,     OPTIONAL :: lonCtr(:,:) ! Lon ctrs [rad]
    REAL(ESMF_KIND_R4),  POINTER,     OPTIONAL :: latCtr(:,:) ! Lat ctrs [rad]

    !----------------------------------
    ! Global grid coordinates
    !----------------------------------
    INTEGER,             INTENT(OUT), OPTIONAL :: IM_WORLD    ! Global # lons
    INTEGER,             INTENT(OUT), OPTIONAL :: JM_WORLD    ! Global # lats
    INTEGER,             INTENT(OUT), OPTIONAL :: LM_WORLD    ! Global # levs
    INTEGER,             INTENT(OUT), OPTIONAL :: IL_WORLD    ! Global start lon index on this PET
    INTEGER,             INTENT(OUT), OPTIONAL :: IU_WORLD    ! Global end   lon index on this PET
    INTEGER,             INTENT(OUT), OPTIONAL :: JL_WORLD    ! Global start lat index on this PET
    INTEGER,             INTENT(OUT), OPTIONAL :: JU_WORLD    ! Global end   lat index on this PET

    !----------------------------------
    ! Date and time variables
    !----------------------------------
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

    !-----------------------------------
    ! Timestep variables [seconds]
    !-----------------------------------
    REAL,                INTENT(OUT), OPTIONAL :: tsChem      ! Chemistry
    REAL,                INTENT(OUT), OPTIONAL :: tsRad       ! RRTMG
    REAL,                INTENT(OUT), OPTIONAL :: tsDyn       ! Dynamics

    !-----------------------------------
    ! Solar parameters
    !-----------------------------------
    REAL,                INTENT(OUT), OPTIONAL :: ZTH(:,:)    ! Solar zth angle
    REAL,                INTENT(OUT), OPTIONAL :: SLR(:,:)    ! Insolation

#if defined( MODEL_GEOS )
    !-----------------------------------------------
    ! Optional import restart file existence inquiry
    !-----------------------------------------------
    LOGICAL,             INTENT(OUT), OPTIONAL :: haveImpRst ! Import rst exist?
#endif

    !-----------------------------------
    ! Return code
    !-----------------------------------
    INTEGER,             INTENT(OUT), OPTIONAL :: RC          ! 0 = all is well
!
! !REMARKS:
!  If you need to obtain a quantity not returned by this routine, you can
!  manually extract it from the MaplCF or GeosCF configuration objects.
!
! !REVISION HISTORY:
!  01 Dec 2009 - A. Da Silva - Initial version
!  See https://github.com/geoschem/geos-chem for history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! LOCAL VARIABLES:
!
    ! Objects
    TYPE(ESMF_Time)               :: startTime      ! ESMF start time obj
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

    ! Scalars
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

    __Iam__('Extract_')

    !=======================================================================
    ! Initialization
    !=======================================================================

    ! Get my name and set-up traceback handle
    CALL ESMF_GridCompGet( GC, name=compName, vm=VM, __RC__ )
    Iam = TRIM(compName)//'::Extract_'

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
#if defined( MODEL_GEOS )
              WRITE(6,*) 'GEOSCHEMCHEM_DT cannot be less than RUN_DT'
#else
              WRITE(6,*) 'Chem_DT cannot be less than RUN_DT'
#endif
           ENDIF
           STATUS = 1
           _VERIFY(STATUS)
        ENDIF
    ENDIF

    If ( PRESENT( tsRad ) .and. PRESENT( tsChem ) ) Then
        _ASSERT(MOD(nint(tsRad),nint(tsChem)) == 0,'RRTMG_DT is not a multiple of GCHPCHEM_DT')
    End If

#if defined( MODEL_GEOS )
    ! Simulation dates. Legacy stuff, not used.
    ! Set to dummy values (ckeller, 1/18/18)
    !IF ( PRESENT( nymdB ) ) nymdB = 20130701
    !IF ( PRESENT( nhmsB ) ) nhmsB = 000000
    !IF ( PRESENT( nymdE ) ) nymdE = 20130701
    !IF ( PRESENT( nhmsE ) ) nhmsE = 000000

    !=======================================================================
    ! Does the import restart file exist?
    !=======================================================================

    ! Import restart file name
    IF ( PRESENT ( haveImpRst ) ) THEN
       CALL ESMF_ConfigGetAttribute( MaplCF, importRstFN,                      &
                                  DEFAULT = "geoschemchem_import_rst.nc4",     &
                                  LABEL = "GEOSCHEMCHEM_IMPORT_RESTART_FILE:", &
                                  __RC__ )

       ! remove bootstrap parameter
       IF( importRstFN(1:1) == '-' .OR. &
           importRstFN(1:1) == '+'       ) THEN
          importRstFN = importRstFN(2:LEN(TRIM(importRstFN)))
       ENDIF
       INQUIRE( FILE=TRIM( importRstFN ), EXIST=haveImpRst )
       IF( MAPL_AM_I_ROOT() ) THEN
          PRINT *," ",TRIM( importRstFN )," exists: ", haveImpRst
          PRINT *," "
       END IF
    ENDIF
#endif

    !=======================================================================
    ! Extract time/date information
    !=======================================================================

    ! Get the ESMF time object
    CALL ESMF_ClockGet( Clock,                    &
                        startTime    = startTime, &
                        stopTime     = stopTime,  &
                        currTime     = currTime,  &
                        advanceCount = count,     &
                         __RC__ )

    ! Get starting-time fields from the time object
    CALL ESMF_TimeGet( startTime, yy=yyyy, mm=mm, dd=dd, dayOfYear=doy, &
                                  h=h,     m=m,   s=s,   __RC__ )

    ! Save fields for return
    IF ( PRESENT( nymd     ) ) CALL MAPL_PackTime( nymd, yyyy, mm, dd )
    IF ( PRESENT( nhms     ) ) CALL MAPL_PackTime( nhms, h,    m,  s  )

    ! Get ending-time fields from the time object
    CALL ESMF_TimeGet( stopTime, yy=yyyy, mm=mm, dd=dd, dayOfYear=doy, &
                                 h=h,     m=m,   s=s,   __RC__ )

    ! Save packed fields for return
    IF ( PRESENT( nymdE    ) ) CALL MAPL_PackTime( nymdE, yyyy, mm, dd )
    IF ( PRESENT( nhmsE    ) ) CALL MAPL_PackTime( nhmsE, h,    m,  s  )

    IF ( PRESENT( advCount ) ) advCount = count

    !=======================================================================
    ! SDE 2017-01-05: The following calls must be kept as a single block,
    ! or the wrong date/time elements will be returned (the yyyy/mm/dd
    ! etc variables are re-used). Specifically, the output variables must
    ! be set now, before the variables are re-used.
    !=======================================================================
    ! Start of current-time block
    !=======================================================================
    ! Get current-time fields from the time object
    CALL ESMF_TimeGet( currTime, yy=yyyy, mm=mm, dd=dd, dayOfYear=doy, &
                                 h=h,     m=m,   s=s,   __RC__ )

    ! Save packed fields for return
    IF ( PRESENT( nymd     ) ) CALL MAPL_PackTime( nymd, yyyy, mm, dd )
    IF ( PRESENT( nhms     ) ) CALL MAPL_PackTime( nhms, h,    m,  s  )

    ! Save the various extacted current-time fields for return
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

    !=======================================================================
    ! End of current-time block
    !=======================================================================

    CALL ESMF_TimeGet( startTime, yy=yyyy, mm=mm, dd=dd, dayOfYear=doy, &
                                 h=h,     m=m,   s=s,   __RC__ )

    ! Save fields for return
    IF ( PRESENT( nymdB    ) ) CALL MAPL_PackTime( nymdB, yyyy, mm, dd )
    IF ( PRESENT( nhmsB    ) ) CALL MAPL_PackTime( nhmsB, h,    m,  s  )

    CALL ESMF_TimeGet( stopTime, yy=yyyy, mm=mm, dd=dd, dayOfYear=doy, &
                                 h=h,     m=m,   s=s,   __RC__ )

    ! Save fields for return
    IF ( PRESENT( nymdE    ) ) CALL MAPL_PackTime( nymdE, yyyy, mm, dd )
    IF ( PRESENT( nhmsE    ) ) CALL MAPL_PackTime( nhmsE, h,    m,  s  )

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

    !=======================================================================
    ! All done
    !=======================================================================
    _RETURN(ESMF_SUCCESS)

  END SUBROUTINE Extract_
!EOC
#ifdef ADJOINT
   subroutine Adjoint_StateRecord( GC, IMPORT, EXPORT, CLOCK, RC )

     ! !ARGUMENTS:

     type(ESMF_GridComp), intent(inout) :: GC     ! composite gridded component
     type(ESMF_State),    intent(inout) :: IMPORT ! import state
     type(ESMF_State),    intent(inout) :: EXPORT ! export state
     type(ESMF_Clock),    intent(inout) :: CLOCK  ! the clock
     integer, optional,   intent(  out) :: RC     ! Error code:
     ! = 0 all is well
     ! otherwise, error
     !EOPI

     ! LOCAL VARIABLES

     character(len=ESMF_MAXSTR)                  :: IAm
     character(len=ESMF_MAXSTR)                  :: COMP_NAME
     integer                                     :: STATUS

     type (MAPL_MetaComp), pointer               :: STATE
     type (ESMF_State)                           :: INTERNAL
     integer                                     :: hdr
     character(len=ESMF_MAXSTR)                  :: FILETYPE
     character(len=ESMF_MAXSTR)                  :: FNAME, DATESTAMP

     !=============================================================================

     !  Begin...

     _UNUSED_DUMMY(EXPORT)

     Iam = "Adjoint_StateRecord"
     call ESMF_GridCompGet(GC, name=COMP_NAME, RC=STATUS )
     _VERIFY(STATUS)
     Iam = trim(COMP_NAME) // Iam

     ! Get my MAPL_Generic state
     ! -------------------------
     CALL MAPL_GetObjectFromGC(GC, STATE, RC=STATUS)
     _VERIFY(STATUS)

     ! Get Internal State
     call MAPL_Get( STATE, INTERNAL_ESMF_STATE=INTERNAL, __RC__ )

     hdr = 0
     ! call MAPL_GetResource( STATE   , hdr,         &
     !      default=0, &
     !      LABEL="INTERNAL_HEADER:", &
     !      RC=STATUS)
     ! _VERIFY(STATUS)

     call MAPL_DateStampGet(clock, datestamp, __RC__ )

     FILETYPE = 'pnc4'
     FNAME = 'gcadj_import_checkpoint.' // trim(datestamp) // '.nc4'

     call MAPL_CheckpointState(IMPORT, CLOCK, &
          FNAME, &
          FILETYPE, STATE, hdr/=0, &
          RC=STATUS)
     _VERIFY(STATUS)

     FNAME = 'gcadj_internal_checkpoint.' // trim(datestamp) // '.nc4'

     call MAPL_CheckpointState(INTERNAL, CLOCK, &
          FNAME, &
          FILETYPE, STATE, hdr/=0, &
          RC=STATUS)
     _VERIFY(STATUS)


     _RETURN(ESMF_SUCCESS)
   end subroutine Adjoint_StateRecord

   subroutine Adjoint_StateRefresh( GC, IMPORT, EXPORT, CLOCK, RC )

     ! !ARGUMENTS:

     type(ESMF_GridComp), intent(inout) :: GC     ! composite gridded component
     type(ESMF_State),    intent(inout) :: IMPORT ! import state
     type(ESMF_State),    intent(inout) :: EXPORT ! export state
     type(ESMF_Clock),    intent(inout) :: CLOCK  ! the clock
     integer, optional,   intent(  out) :: RC     ! Error code:
     ! = 0 all is well
     ! otherwise, error
     !EOPI

     ! LOCAL VARIABLES

     character(len=ESMF_MAXSTR)                  :: IAm
     character(len=ESMF_MAXSTR)                  :: COMP_NAME
     integer                                     :: STATUS

     type (MAPL_MetaComp), pointer               :: STATE
     type (ESMF_State)                           :: INTERNAL
     integer                                     :: hdr
     integer                                     :: unit

     character(len=ESMF_MAXSTR)                  :: FNAME, datestamp

     !=============================================================================

     _UNUSED_DUMMY(EXPORT)

     !  Begin...

     Iam = "Adjoint_StateRefresh"
     call ESMF_GridCompGet(GC, name=COMP_NAME, RC=STATUS )
     _VERIFY(STATUS)
     Iam = trim(COMP_NAME) // Iam

     ! Get my MAPL_Generic state
     ! -------------------------
     CALL MAPL_GetObjectFromGC(GC, STATE, RC=STATUS)
     _VERIFY(STATUS)

     ! Get Internal state
     CALL MAPL_Get ( STATE, INTERNAL_ESMF_STATE=INTERNAL, __RC__ )

     call MAPL_DateStampGet(clock, datestamp, rc=status)
     _VERIFY(STATUS)

     HDR = 0

     FNAME = 'gcadj_import_checkpoint.' // trim(datestamp) // '.nc4'

     call MAPL_ESMFStateReadFromFile(IMPORT, CLOCK, &
          FNAME, &
          STATE, .FALSE., RC=STATUS)
     _VERIFY(STATUS)
     UNIT = GETFILE(FNAME, RC=STATUS)
     _VERIFY(STATUS)
     call MAPL_DestroyFile(unit = UNIT, rc=STATUS)
     _VERIFY(STATUS)
     CALL FREE_FILE(UNIT, RC=STATUS)
     _VERIFY(STATUS)

     FNAME = 'gcadj_internal_checkpoint.' // trim(datestamp) // '.nc4'

     call MAPL_ESMFStateReadFromFile(INTERNAL, CLOCK, &
          FNAME, &
          STATE, hdr/=0, RC=STATUS)
     _VERIFY(STATUS)
     IF (FNAME(1:1) .eq. '-' .or. &
          FNAME(1:1) .eq. '+') THEN
        UNIT = GETFILE(FNAME(2:), RC=STATUS)
     else
        UNIT = GETFILE(FNAME, RC=STATUS)
     endif
     _VERIFY(STATUS)
     call MAPL_DestroyFile(unit = UNIT, rc=STATUS)
     _VERIFY(STATUS)
     CALL FREE_FILE(UNIT, RC=STATUS)
     _VERIFY(STATUS)

     _RETURN(ESMF_SUCCESS)
   end subroutine Adjoint_StateRefresh

#endif

#ifdef MODEL_GEOS
END MODULE GEOSCHEMchem_GridCompMod
#else
END MODULE Chem_GridCompMod
#endif
