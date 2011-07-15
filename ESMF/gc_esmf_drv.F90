#if defined (ESMF_)
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: gc_chunk_mod
!
! !DESCRIPTION: Module GC\_CHUNK\_MOD is the module that contains 
!  the GEOS-Chem chunk code init, run and finalize methods.
!\\
!\\
! !INTERFACE: 
!      
#include "MAPL_Generic.h"

MODULE gc_esmf_drv
!
! !USES: 
!
  USE ESMF_Mod                                     ! ESMF framework
  USE MAPL_Mod                                     ! MAPL framework
  USE GC_ESMF_COMP, ONLY: SetServices  ! To set IRF methods
#if defined (ESMF_TESTBED_)
  USE TESTBED_VALUE_MOD                            ! GEOS-Chem input values
!  USE GC_ESMF_TESTBED
#endif 


  IMPLICIT NONE
!
! !REMARKS:
!  ut_GEOSCHEM - Simple ESMF/MAPL example demonstrating how to call GEOSCHEM
!                                                                             .
!  It assumes 2 processors, so typically you will run it as
!
!     % mprirun -np 2 ut_GEOSCHEM.x
!                                                                             .
!  Arlindo da Silva <arlindo.dasilva@nasa.gov>, December 2009
!
! !REVISION HISTORY: 
!  01 Dec 2009 - A. Da Silva - Initial version  
!  02 Apr 2010 - R. Yantosca - Modified for GEOS-Chem column code
!  02 Apr 2010 - R. Yantosca - Added ProTex Headers, other cosmetic changes
!  07 Apr 2010 - R. Yantosca - Now populate all import/internal state fields
!  08 Apr 2010 - R. Yantosca - Now reference GEOS-Chem column input values from
!                              module "bmy_GC_Value_Mod.F90".  The header file
!                              "Column_Values_Saved_From_GEOS-Chem.h" is 
!                              now obsolete.
!  01 Jun 2010 - R. Yantosca - Increased output format to preserve precision
!  01 Jul 2010 - R. Yantosca - Now zero D_OH_MASS and D_AIR_MASS
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
  !==========================================================================
  ! Basic ESMF objects being used in this example
  !==========================================================================
  TYPE(ESMF_Grid)         :: Grid            ! Grid
  TYPE(ESMF_VM)           :: VM              ! ESMF Virtual Machine
  TYPE(ESMF_Time)         :: startTime       ! ESMF Time object (start time)
  TYPE(ESMF_Time)         :: stopTime        ! ESMF Time object (stop time)
  TYPE(ESMF_TimeInterval) :: TimeStep        ! ESMF TimeStep object 

  !==========================================================================
  ! Grid component objects
  !==========================================================================
  TYPE(ESMF_GridComp)     :: GrComp          ! ESMF Gridded component
  TYPE(ESMF_State)        :: Import          ! ESMF Import state
  TYPE(ESMF_State)        :: Export
  TYPE(ESMF_Clock)        :: Clock

  !==========================================================================
  ! Basic information about the parallel environment
  ! PET = Persistent Execution Threads
  ! In the current implementation, a PET is equivalent to an MPI process
  !==========================================================================
  INTEGER                 :: myPET           ! The local PET #
  INTEGER                 :: nPET            ! Total # of PETs we are using
  INTEGER                 :: STATUS          ! Status variable
  INTEGER                 :: RC              ! Status variable
  INTEGER                 :: I,  J, N        ! Loop indices 
  INTEGER                 :: IM, JM          ! Loop indices
  INTEGER, PARAMETER      :: NX       = 2    ! Layout: # of PETs in longitude
  INTEGER, PARAMETER      :: NY       = 1    ! Layout: # of PETs in latitude
  INTEGER, PARAMETER      :: IM_WORLD = 72   ! # of longitudes in global grid
  INTEGER, PARAMETER      :: JM_WORLD = 46   ! # of latitudes in global grid
  INTEGER, PARAMETER      :: LM_WORLD = 72   ! # of levels in global grid
  
  !==========================================================================
  ! Character variables
  !==========================================================================
  CHARACTER(LEN=ESMF_MAXSTR)  :: name
  CHARACTER(LEN=*), PARAMETER :: Iam = 'ut_GEOSCHEM'

    CONTAINS
!EOC
!------------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 910.1 and      !
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: GC_Main_Drv
!
! !DESCRIPTION: The Main subroutine does the following:
!
! \begin{enumerate}
! \item Initializes the ESMF and MAPL frameworks
! \item Creates an ESMF Grid object
! \item Creates an ESMF Clock object
! \item Creates and initializes the ESMF Import and Export States
! \item Creates the GEOS-Chem gridded component
! \item Initializes the internal state of the GEOS-Chem gridded component
! \item Runs the GEOS-Chem gridded component
! \item Finalizes the GEOS-Chem gridded component
! \item Finalizes ESMF and MAPL frameworks
! \end{enumerate}
!
! !INTERFACE:
!
  SUBROUTINE gc_main_drv
! !REMARKS:
!`
! !REVISION HISTORY: 
!  12 Jul 2011 - M. Long     - Initial version built from R. Yantosca Sand-box Version 
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    TYPE(ESMF_Config)  :: MaplCF          
    TYPE(ESMF_Config)  :: GeosCF          
    INTEGER            :: runDt,         count
    INTEGER            :: utcStartDate,  utcStartTime
    INTEGER            :: utcEndDate,    utcEndTime
    INTEGER            :: yy0, mm0, dd0, h0, m0, s0
    INTEGER            :: yy1, mm1, dd1, h1, m1, s1
    CHARACTER(LEN=255) :: fileName

    !========================================================================
    ! Initialize
    !========================================================================

    ! Initialize ESMF framework.  For performance reasons, it is
    ! important to turn OFF ESMF's automatic logging feature
    CALL ESMF_Initialize( DefaultLogType=ESMF_LOG_NONE, VM=VM, __RC__ )
    
    ! Get the total # of PETs and local PET index from the ESMF VM 
    CALL ESMF_VMGet( VM, localPET=myPET, PETcount=nPET )  

    ! Quit if we have the wrong # of processors
    IF ( nPET /= 2 ) THEN
       IF ( MAPL_am_I_root() ) THEN
          PRINT*, 'Error: expecting 2 PETs but found ', nPET, 'PETs'
          PRINT*, 'Try:   mpirun -np 2 ut_GEOSCHEM.x'
       ENDIF
       ASSERT_(.FALSE.)
    ENDIF

    ! Echo info if we are on the root PET
    IF ( MAPL_AM_I_ROOT() ) THEN
       WRITE( 6, '(a)' ) REPEAT( '=', 79 )
       WRITE( 6, 100   ) Iam, nPET
100    FORMAT( 'Starting ', a, ' with ', i3, ' PETs ...' )
       WRITE( 6, '(a)' ) REPEAT( '=', 79 )
    ENDIF

    !-----------------------------------
    ! Create the ESMF Grid object
    !-----------------------------------

    ! Create a global 2D Lat-Lon grid on a 2x1 layout
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%% NOTE: This grid does not match up w/ the GEOS-CHEM 4x5 grid.  %%%
    !%%% But this is how the GEOS-5 met fields come in.  Keep it as-is %%%
    !%%% for testing and we'll deal w/ it later. (bmy, 4/12/10)        %%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Grid = MAPL_LatLonGridCreate( name     = 'GEOS-Chem Grid', &
                                  Nx       = NX,               &
                                  Ny       = NY,               &
                                  IM_WORLD = IM_WORLD,         &
                                  JM_WORLD = JM_WORLD,         &
                                  LM_WORLD = LM_WORLD,         &
                                  __RC__ )

    ! Validate the grid
    CALL ESMF_GridValidate( Grid, __RC__ )

    !-----------------------------------
    ! Read start/stop time info
    !-----------------------------------

    ! Create config objects
    MaplCF = ESMF_ConfigCreate( __RC__ ) !MSL NAMELIST OBJECT
    GeosCF = ESMF_ConfigCreate( __RC__ ) !MSL NAMELIST OBJECT

    ! Load config objects   
    CALL ESMF_ConfigLoadFile( MaplCF, 'MAPL.rc',                  __RC__ )
    CALL ESMF_ConfigLoadFile( GeosCF, 'GEOSCHEMchem_GridComp.rc', __RC__ )

    ! Timestep
    CALL ESMF_ConfigGetAttribute( MaplCF, runDt,               &
                                  Label   = "RUN_DT:",         &
                                  Default = 1800,              &
                                  __RC__ )

    ! Start date
    CALL ESMF_ConfigGetAttribute( GeosCF, utcStartDate,        &
                                  Label   = "UTC_START_DATE:", &
                                  Default = 20080101,          &
                                   __RC__ )

    ! Start time
    CALL ESMF_ConfigGetAttribute( GeosCF, utcStartTime,        &
                                  LABEL   = "UTC_START_TIME:", &
                                  Default = 000000,            &
                                   __RC__ )

    ! End date
    CALL ESMF_ConfigGetAttribute( GeosCF, utcEndDate,          &
                                  Label   = "UTC_END_DATE:",   &
                                  Default = 20080101,          &
                                   __RC__ )

    ! End time
    CALL ESMF_ConfigGetAttribute( GeosCF, utcEndTime,          &
                                  LABEL   = "UTC_END_TIME:",   &
                                  Default = 010000,            &
                                   __RC__ )

    ! Split time variables for defining the clock 
    CALL MAPL_UnpackTime( utcStartDate, yy0, mm0, dd0 )
    CALL MAPL_UnpackTime( utcStartTime, h0,  m0,  s0  )
    CALL MAPL_UnpackTime( utcEndDate,   yy1, mm1, dd1 )
    CALL MAPL_UnpackTime( utcEndTime,   h1,  m1,  s1  )

    ! Free the config objects
    CALL ESMF_ConfigDestroy( MaplCF, __RC__ )
    CALL ESMF_ConfigDestroy( GeosCF, __RC__ )

    ! Convert MAPL timestep from seconds to mins
    runDt = runDt / 60

    !-----------------------------------
    ! Create the ESMF Clock object
    !-----------------------------------

    ! Set the calendar type (e.g. "Gregorian")
    CALL ESMF_CalendarSetDefault( ESMF_CAL_GREGORIAN )

    ! Set the starting time (Jul 1, 2008)
    CALL ESMF_TimeSet( startTime, yy=yy0, mm=mm0, dd=dd0, h=h0, m=m0, s=s0 )
    CALL ESMF_TimeSet( stopTime , yy=yy1, mm=mm1, dd=dd1, h=h1, m=m1, s=s1 )

    ! Set the time step (30 mins)
    CALL ESMF_TimeIntervalSet( TimeStep, h=0, m=runDt, s=0, __RC__ )

    ! Create the ESMF clock
    Clock = ESMF_ClockCreate( "Clock",               &
                              timeStep  = TimeStep,  & 
                              startTime = startTime, &
                              stopTime  = stopTime,  &
                               __RC__ )

    !-----------------------------------
    ! States
    !-----------------------------------

    ! Create import state
    Import = ESMF_StateCreate( stateName='impGEOSCHEM', &
                               stateType=ESMF_STATE_IMPORT, __RC__ )

    ! Create export state
    Export = ESMF_StateCreate( stateName='expGEOSCHEM', &
                               stateType=ESMF_STATE_EXPORT, __RC__ )

    !-----------------------------------
    ! Gridded component
    !-----------------------------------

    ! Create gridded component 
    GrComp = ESMF_GridCompCreate( name         = 'GEOSCHEMchem', &
                                  Grid         = Grid,           &
                                  GridCompType = ESMF_ATM,       &
                                  ConfigFile   = 'MAPL.rc',      &
                                  __RC__  )

    ! Set the component services.  This identifies the Initialize,
    ! Run, and Finalize methods to the ESMF and MAPL frameworks.
    CALL ESMF_GridCompSetServices( GrComp, SetServices, __RC__ )

    ! Initialize the gridded component
    CALL ESMF_GridCompInitialize( GrComp, Import, Export, Clock, __RC__ )

#if defined (ESMF_TESTBED_)
    ! Since we are not reading restarts, set the internal state with
    ! reasonable profiles so that we can exercise the code
    ! This is called just once at the start of the simulation
    fileName = 'internal.state.72L.00000'
    CALL Fill_Internal_State_( GrComp, fileName,  __RC__ )
#endif

    !========================================================================
    ! Run
    !========================================================================

    ! timestep counter
    count = 0

    ! Timestepping loop
    DO WHILE ( .not. ESMF_ClockIsStopTime( Clock, RC ) )

#if defined (ESMF_TESTBED_)
       ! Filename for the import state
       WRITE( fileName, 20 ) count
 20    FORMAT( 'import.state.72L.', i5.5 )

       ! Fill the import state w/ values taken from GEOS-Chem
       ! This is called every timestep to evolve the meteorology
       CALL Fill_Import_State_( Import, fileName, __RC__ )
#endif

       ! Call the Run method of the GEOS-Chem gridded component
       CALL ESMF_GridCompRun( GrComp, Import, Export, Clock, __RC__ )

       ! Advance the clock by one timestep
       CALL ESMF_ClockAdvance( Clock, __RC__ )

       ! Counter
       count = count + 1
    ENDDO

    !========================================================================
    ! Finalize
    !========================================================================

    ! Call the Finalize method of the GEOS-Chem gridded component
    CALL ESMF_GridCompFinalize( GrComp, Import, Export, Clock, __RC__ )

    ! Finalize ESMF
    CALL ESMF_Finalize( __RC__ )

    END SUBROUTINE gc_main_drv

!------------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 910.1 and      !
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Fill_Import_State_
!
! !DESCRIPTION: Subroutine Fill_Import_State_ populates the import state of the
!  GEOS-Chem gridded component with reasonable values for the unit test.
!\\
!\\
!  Normally, the GEOS-Chem gridded component would recieve these fields
!  from the GEOS-5 GCM via the import state.  However, to simplify the unit
!  testing, we instead will declare all import state variables here.
!
! !INTERFACE:
!
  SUBROUTINE Fill_Import_State_( Import, inFile, RC )
!
! !USES:
!
#   include "GEOSCHEMchem_DeclarePointer___.h"     ! Pointer decl's to States
!
! !INPUT PARAMETERS:
!
    CHARACTER(LEN=*), INTENT(IN)    :: inFile      ! Input file name
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ESMF_State), INTENT(INOUT) :: Import      ! Import state
!
! !OUTPUT PARAMETERS:
! 
    INTEGER,          INTENT(OUT)   :: RC          ! Return code
!
! !REMARKS:
!  For simplicity, we assign data from a single column to the entire world.
!  Here we are mostly concerned not w/ the performace of the GEOS-Chem code
!  but to see if it will link up with the GEOS-5 GCM.  
!
! !REVISION HISTORY: 
!  01 Dec 2009 - A. Da Silva - Initial version  
!  02 Apr 2010 - R. Yantosca - Modified for GEOS-Chem column code
!  02 Apr 2010 - R. Yantosca - Added ProTex Headers, other cosmetic changes
!  06 Apr 2010 - R. Yantosca - Fill import state fields w/ G-C output values
!  01 Jun 2010 - R. Yantosca - Changed output format to 4(es19.12,1x)
!  04 Jun 2010 - R. Yantosca - Bug fix: now read tag "MET_SPHU"
!  04 Jun 2010 - R. Yantosca - Now init ZPBL with GC_ZPBL field
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER            :: I,  J    ! Loop indices
    INTEGER            :: I1, I2   ! Min & max lon indices on this PET
    INTEGER            :: J1, J2   ! Min & max lat indices on this PET
    INTEGER            :: K1, K2   ! Min & max lev indices on this PET
    INTEGER            :: IOS
    CHARACTER(LEN=255) :: NAME
!
! !DEFINED PARAMETERS:
!
    INTEGER, PARAMETER :: LUN=200

    !========================================================================
    ! Get Pointers to IMPORT state fields
    !========================================================================
    CALL MAPL_GetPointer( Import, ALBEDO,         'ALBEDO',          __RC__ )
    CALL MAPL_GetPointer( Import, CLDFRC,         'CLDFRC',          __RC__ )
    CALL MAPL_GetPointer( Import, SH,             'SH',              __RC__ )
    CALL MAPL_GetPointer( Import, FRLAND,         'FRLAND',          __RC__ )
    CALL MAPL_GetPointer( Import, LWI,            'LWI',             __RC__ )
    CALL MAPL_GetPointer( Import, ZPBL,           'ZPBL',            __RC__ )
    CALL MAPL_GetPointer( Import, CN_PRCP,        'CN_PRCP',         __RC__ )
    CALL MAPL_GetPointer( Import, LS_PRCP,        'LS_PRCP',         __RC__ )
    CALL MAPL_GetPointer( Import, AN_PRCP,        'AN_PRCP',         __RC__ )
    CALL MAPL_GetPointer( Import, TROPP,          'TROPP',           __RC__ )
    CALL MAPL_GetPointer( Import, T2M,            'T2M',             __RC__ )
    CALL MAPL_GetPointer( Import, RADSRF,         'RADSRF',          __RC__ )
    CALL MAPL_GetPointer( Import, TSEA,           'TSEA',            __RC__ )
    CALL MAPL_GetPointer( Import, U10M,           'U10M',            __RC__ )
    CALL MAPL_GetPointer( Import, USTAR,          'USTAR',           __RC__ )
    CALL MAPL_GetPointer( Import, V10M,           'V10M',            __RC__ )
    CALL MAPL_GetPointer( Import, Z0,             'Z0',              __RC__ )
    CALL MAPL_GetPointer( Import, FCLDF,          'FCLDF',           __RC__ )
    CALL MAPL_GetPointer( Import, CNV_FMC,        'CNV_FMC',         __RC__ )
    CALL MAPL_GetPointer( Import, CNV_MFD,        'CNV_MFD',         __RC__ )
    CALL MAPL_GetPointer( Import, TAUCLI,         'TAUCLI',          __RC__ )
    CALL MAPL_GetPointer( Import, TAUCLW,         'TAUCLW',          __RC__ )
    CALL MAPL_GetPointer( Import, RH2,            'RH2',             __RC__ )
    CALL MAPL_GetPointer( Import, Q,              'Q',               __RC__ )
    CALL MAPL_GetPointer( Import, T,              'T',               __RC__ )
    CALL MAPL_GetPointer( Import, AIRDENS,        'AIRDENS',         __RC__ )
    CALL MAPL_GetPointer( Import, COSZ,           'COSZ',            __RC__ )
    CALL MAPL_GetPointer( Import, MOISTQ,         'MOISTQ',          __RC__ )
    CALL MAPL_GetPointer( Import, DELP,           'DELP',            __RC__ )
    CALL MAPL_GetPointer( Import, PLE,            'PLE',             __RC__ )
    CALL MAPL_GetPointer( Import, PL,             'PL',              __RC__ )
    CALL MAPL_GetPointer( Import, TO3,            'TO3',             __RC__ )
    CALL MAPL_GetPointer( Import, UVALBEDO,       'UVALBEDO',        __RC__ )
    CALL MAPL_GetPointer( Import, EMISS_NOX,      'EMISS_NOX',       __RC__ )
    CALL MAPL_GetPointer( Import, EMISS_O3,       'EMISS_O3',        __RC__ )
    CALL MAPL_GetPointer( Import, EMISS_CO,       'EMISS_CO',        __RC__ )
    CALL MAPL_GetPointer( Import, EMISS_ALK4,     'EMISS_ALK4',      __RC__ )
    CALL MAPL_GetPointer( Import, EMISS_ISOP,     'EMISS_ISOP',      __RC__ )
    CALL MAPL_GetPointer( Import, EMISS_HNO3,     'EMISS_HNO3',      __RC__ )
    CALL MAPL_GetPointer( Import, EMISS_ACET,     'EMISS_ACET',      __RC__ )
    CALL MAPL_GetPointer( Import, EMISS_MEK,      'EMISS_MEK',       __RC__ )
    CALL MAPL_GetPointer( Import, EMISS_ALD2,     'EMISS_ALD2',      __RC__ )
    CALL MAPL_GetPointer( Import, EMISS_PRPE,     'EMISS_PRPE',      __RC__ )
    CALL MAPL_GetPointer( Import, EMISS_C3H8,     'EMISS_C3H8',      __RC__ )
    CALL MAPL_GetPointer( Import, EMISS_CH2O,     'EMISS_CH2O',      __RC__ )
    CALL MAPL_GetPointer( Import, EMISS_C2H6,     'EMISS_C2H6',      __RC__ )
    CALL MAPL_GetPointer( Import, EMISS_DMS,      'EMISS_DMS',       __RC__ )
    CALL MAPL_GetPointer( Import, EMISS_SO2,      'EMISS_SO2',       __RC__ )
    CALL MAPL_GetPointer( Import, EMISS_SO4,      'EMISS_SO4',       __RC__ )
    CALL MAPL_GetPointer( Import, EMISS_NH3,      'EMISS_NH3',       __RC__ )
    CALL MAPL_GetPointer( Import, EMISS_BCPI,     'EMISS_BCPI',      __RC__ )
    CALL MAPL_GetPointer( Import, EMISS_OCPI,     'EMISS_OCPI',      __RC__ )
    CALL MAPL_GetPointer( Import, EMISS_BCPO,     'EMISS_BCPO',      __RC__ )
    CALL MAPL_GetPointer( Import, EMISS_OCPO,     'EMISS_OCPO',      __RC__ )
    CALL MAPL_GetPointer( Import, EMISS_ALPH,     'EMISS_ALPH',      __RC__ )
    CALL MAPL_GetPointer( Import, EMISS_LIMO,     'EMISS_LIMO',      __RC__ )
    CALL MAPL_GetPointer( Import, EMISS_ALCO,     'EMISS_ALCO',      __RC__ )
    CALL MAPL_GetPointer( Import, EMISS_DST1,     'EMISS_DST1',      __RC__ )
    CALL MAPL_GetPointer( Import, EMISS_DST2,     'EMISS_DST2',      __RC__ )
    CALL MAPL_GetPointer( Import, EMISS_DST3,     'EMISS_DST3',      __RC__ )
    CALL MAPL_GetPointer( Import, EMISS_DST4,     'EMISS_DST4',      __RC__ )
    CALL MAPL_GetPointer( Import, EMISS_SALA,     'EMISS_SALA',      __RC__ )
    CALL MAPL_GetPointer( Import, EMISS_SALC,     'EMISS_SALC',      __RC__ )
    CALL MAPL_GetPointer( Import, IREG,           'IREG',            __RC__ )
    CALL MAPL_GetPointer( Import, ILAND_01,       'ILAND_01',        __RC__ )
    CALL MAPL_GetPointer( Import, ILAND_02,       'ILAND_02',        __RC__ )
    CALL MAPL_GetPointer( Import, ILAND_03,       'ILAND_03',        __RC__ )
    CALL MAPL_GetPointer( Import, ILAND_04,       'ILAND_04',        __RC__ )
    CALL MAPL_GetPointer( Import, ILAND_05,       'ILAND_05',        __RC__ )
    CALL MAPL_GetPointer( Import, ILAND_06,       'ILAND_06',        __RC__ )
    CALL MAPL_GetPointer( Import, ILAND_07,       'ILAND_07',        __RC__ )
    CALL MAPL_GetPointer( Import, ILAND_08,       'ILAND_08',        __RC__ )
    CALL MAPL_GetPointer( Import, ILAND_09,       'ILAND_09',        __RC__ )
    CALL MAPL_GetPointer( Import, ILAND_10,       'ILAND_10',        __RC__ )
    CALL MAPL_GetPointer( Import, ILAND_11,       'ILAND_11',        __RC__ )
    CALL MAPL_GetPointer( Import, ILAND_12,       'ILAND_12',        __RC__ )
    CALL MAPL_GetPointer( Import, ILAND_13,       'ILAND_13',        __RC__ )
    CALL MAPL_GetPointer( Import, ILAND_14,       'ILAND_14',        __RC__ )
    CALL MAPL_GetPointer( Import, ILAND_15,       'ILAND_15',        __RC__ )
    CALL MAPL_GetPointer( Import, IUSE_01,        'IUSE_01',         __RC__ )
    CALL MAPL_GetPointer( Import, IUSE_02,        'IUSE_02',         __RC__ )
    CALL MAPL_GetPointer( Import, IUSE_03,        'IUSE_03',         __RC__ )
    CALL MAPL_GetPointer( Import, IUSE_04,        'IUSE_04',         __RC__ )
    CALL MAPL_GetPointer( Import, IUSE_05,        'IUSE_05',         __RC__ )
    CALL MAPL_GetPointer( Import, IUSE_06,        'IUSE_06',         __RC__ )
    CALL MAPL_GetPointer( Import, IUSE_07,        'IUSE_07',         __RC__ )
    CALL MAPL_GetPointer( Import, IUSE_08,        'IUSE_08',         __RC__ )
    CALL MAPL_GetPointer( Import, IUSE_09,        'IUSE_09',         __RC__ )
    CALL MAPL_GetPointer( Import, IUSE_10,        'IUSE_10',         __RC__ )
    CALL MAPL_GetPointer( Import, IUSE_11,        'IUSE_11',         __RC__ )
    CALL MAPL_GetPointer( Import, IUSE_12,        'IUSE_12',         __RC__ )
    CALL MAPL_GetPointer( Import, IUSE_13,        'IUSE_13',         __RC__ )
    CALL MAPL_GetPointer( Import, IUSE_14,        'IUSE_14',         __RC__ )
    CALL MAPL_GetPointer( Import, IUSE_15,        'IUSE_15',         __RC__ )
    CALL MAPL_GetPointer( Import, LAI_01,         'LAI_01',          __RC__ )
    CALL MAPL_GetPointer( Import, LAI_02,         'LAI_02',          __RC__ )
    CALL MAPL_GetPointer( Import, LAI_03,         'LAI_03',          __RC__ )
    CALL MAPL_GetPointer( Import, LAI_04,         'LAI_04',          __RC__ )
    CALL MAPL_GetPointer( Import, LAI_05,         'LAI_05',          __RC__ )
    CALL MAPL_GetPointer( Import, LAI_06,         'LAI_06',          __RC__ )
    CALL MAPL_GetPointer( Import, LAI_07,         'LAI_07',          __RC__ )
    CALL MAPL_GetPointer( Import, LAI_08,         'LAI_08',          __RC__ )
    CALL MAPL_GetPointer( Import, LAI_09,         'LAI_09',          __RC__ )
    CALL MAPL_GetPointer( Import, LAI_10,         'LAI_10',          __RC__ )
    CALL MAPL_GetPointer( Import, LAI_11,         'LAI_11',          __RC__ )
    CALL MAPL_GetPointer( Import, LAI_12,         'LAI_12',          __RC__ )
    CALL MAPL_GetPointer( Import, LAI_13,         'LAI_13',          __RC__ )
    CALL MAPL_GetPointer( Import, LAI_14,         'LAI_14',          __RC__ )
    CALL MAPL_GetPointer( Import, LAI_15,         'LAI_15',          __RC__ )
    CALL MAPL_GetPointer( Import, SOX_OH,         'SOX_OH',          __RC__ )
    CALL MAPL_GetPointer( Import, SOX_JV_NOx,     'SOX_JV_NOx',      __RC__ )
    CALL MAPL_GetPointer( Import, SOX_JV_H2O2,    'SOX_JV_H2O2',     __RC__ )
    CALL MAPL_GetPointer( Import, SOX_JV_ACET,    'SOX_JV_ACET',     __RC__ )
    CALL MAPL_GetPointer( Import, SOX_JV_MEK,     'SOX_JV_MEK',      __RC__ )
    CALL MAPL_GetPointer( Import, SOX_JV_ALD2,    'SOX_JV_ALD2',     __RC__ )
    CALL MAPL_GetPointer( Import, SOX_JV_RCHO,    'SOX_JV_RCHO',     __RC__ )
    CALL MAPL_GetPointer( Import, SOX_JV_MVK,     'SOX_JV_MVK',      __RC__ )
    CALL MAPL_GetPointer( Import, SOX_JV_MACR,    'SOX_JV_MACR',     __RC__ )
    CALL MAPL_GetPointer( Import, SOX_JV_R4N2,    'SOX_JV_R4N2',     __RC__ )
    CALL MAPL_GetPointer( Import, SOX_JV_CH2O,    'SOX_JV_CH2O',     __RC__ )
    CALL MAPL_GetPointer( Import, SOX_JV_N2O5,    'SOX_JV_N2O5',     __RC__ )
    CALL MAPL_GetPointer( Import, SOX_JV_HNO4,    'SOX_JV_HNO4',     __RC__ )
    CALL MAPL_GetPointer( Import, SOX_JV_MP,      'SOX_JV_MP',       __RC__ )
    CALL MAPL_GetPointer( Import, SOX_PCO,        'SOX_PCO',         __RC__ )
    CALL MAPL_GetPointer( Import, SOX_LCO,        'SOX_LCO',         __RC__ )

    ! Make sure the pointer is allocated
    ASSERT_( ASSOCIATED( T ) )

    ! Get indices from a 3-D data block
    I1 = LBOUND( T, 1 )
    J1 = LBOUND( T, 2 )
    K1 = LBOUND( T, 3 )
    I2 = UBOUND( T, 1 )
    J2 = UBOUND( T, 2 )
    K2 = UBOUND( T, 3 )

    ! Make sure the horizontal dimensions of the array are OK
    ASSERT_( I1 > 0 )       
    ASSERT_( I2 > 0 )
    ASSERT_( J1 > 0 )
    ASSERT_( J2 > 0 )

    ! Make sure we have 72 levels
    ASSERT_( ( K2 - K1 + 1 ) == 72 )

    !========================================================================
    ! Read data from disk and store in placeholders
    !========================================================================
    PRINT*, '### READING: ', TRIM( inFile ), 'on PET: ', myPet

    ! Open file
    OPEN( LUN, FILE=TRIM( inFile ), STATUS='OLD' )
 10 FORMAT( 4( es19.12, 1x ) )

    DO 
       ! Read the name from the file
       READ( LUN, '(a)', IOSTAT=IOS ) NAME
       IF ( IOS < 0 ) EXIT

       ! Read fields from file
       SELECT CASE( TRIM( NAME ) ) 

          ! Met fields
          CASE( 'MET_ALBD' )   
             READ( LUN, 10 ) GC_ALBD
          CASE( 'MET_AREA_M2' )   
             READ( LUN, 10 ) GC_AREA_M2
          CASE( 'MET_CLDFRC' )    
             READ( LUN, 10 ) GC_CLDFRC
          CASE( 'MET_FRCLND' )  
             READ( LUN, 10 ) GC_FRCLND          
          CASE( 'MET_GWETTOP' )  
             READ( LUN, 10 ) GC_GWETTOP
          CASE( 'MET_HFLUX' )     
             READ( LUN, 10 ) GC_HFLUX
          CASE( 'MET_LWI' )     
             READ( LUN, 10 ) GC_LWI
          CASE( 'MET_PARDR' )     
             READ( LUN, 10 ) GC_PARDR
          CASE( 'MET_PARDF' )     
             READ( LUN, 10 ) GC_PARDF
          CASE( 'MET_PBLH', 'MET_ZPBL' )     
             READ( LUN, 10 ) GC_ZPBL
          CASE( 'MET_PRECCON' )     
             READ( LUN, 10 ) GC_PRECCON
          CASE( 'MET_PRECTOT' )     
             READ( LUN, 10 ) GC_PRECTOT
          CASE( 'MET_RADSWG' )     
             READ( LUN, 10 ) GC_RADSWG
          CASE( 'MET_SST' )     
             READ( LUN, 10 ) GC_SST
          CASE( 'MET_SUNCOS' )     
             READ( LUN, 10 ) GC_COSZ
          CASE( 'MET_TO3' )     
             READ( LUN, 10 ) GC_TO3
          CASE( 'MET_TROPP' )     
             READ( LUN, 10 ) GC_TROPP
          CASE( 'MET_TS' )     
             READ( LUN, 10 ) GC_TS
          CASE( 'MET_U10M' )     
             READ( LUN, 10 ) GC_U10M
          CASE( 'MET_USTAR' )     
             READ( LUN, 10 ) GC_USTAR
          CASE( 'MET_UVALBEDO' )     
             READ( LUN, 10 ) GC_UVALBEDO
          CASE( 'MET_V10M' )     
             READ( LUN, 10 ) GC_V10M
          CASE( 'MET_Z0' )     
             READ( LUN, 10 ) GC_Z0
          CASE( 'MET_AD' )
             READ( LUN, 10 ) GC_AD
          CASE( 'MET_AIRDENS' )
             READ( LUN, 10 ) GC_AIRDENS
          CASE( 'MET_AIRVOL' )
             READ( LUN, 10 ) GC_AIRVOL
          CASE( 'MET_BXHEIGHT' )
             READ( LUN, 10 ) GC_BXHEIGHT
          CASE( 'MET_CLDF' )
             READ( LUN, 10 ) GC_CLDF
          CASE( 'MET_CMFMC' )  
             READ( LUN, 10 ) GC_CMFMC
          CASE( 'MET_DELP' )    
             READ( LUN, 10 ) GC_DELP
          CASE( 'MET_DQIDTMST' )
             READ( LUN, 10 ) GC_DQIDTMST
          CASE( 'MET_DQLDTMST' )
             READ( LUN, 10 ) GC_DQLDTMST
          CASE( 'MET_DQVDTMST' )
             READ( LUN, 10 ) GC_DQVDTMST
          CASE( 'MET_DTRAIN' )  
             READ( LUN, 10 ) GC_DTRAIN
          CASE( 'MET_MOISTQ' )  
             READ( LUN, 10 ) GC_MOISTQ
          CASE( 'MET_OPTD' ) 
             READ( LUN, 10 ) GC_OPTD
          CASE( 'MET_PMID' )    
             READ( LUN, 10 ) GC_PMID
          CASE( 'MET_PEDGE' )   
             READ( LUN, 10 ) GC_PEDGE
          CASE( 'MET_RH' )      
             READ( LUN, 10 ) GC_RH
          CASE( 'MET_Q', 'MET_SPHU' )       
             READ( LUN, 10 ) GC_Q
          CASE( 'MET_T' )       
             READ( LUN, 10 ) GC_T
          CASE( 'MET_TAUCLI' )  
             READ( LUN, 10 ) GC_TAUCLI
          CASE( 'MET_TAUCLW' )  
             READ( LUN, 10 ) GC_TAUCLW

          ! Land & leaf info
          CASE( 'IREG' )     
             READ( LUN, 10 ) GC_IREG
          CASE( 'ILAND01' )     
             READ( LUN, 10 ) GC_ILAND01
          CASE( 'ILAND02'  )     
             READ( LUN, 10 ) GC_ILAND02
          CASE( 'ILAND03'  )     
             READ( LUN, 10 ) GC_ILAND03
          CASE( 'ILAND04'  )     
             READ( LUN, 10 ) GC_ILAND04
          CASE( 'ILAND05'  )     
             READ( LUN, 10 ) GC_ILAND05
          CASE( 'ILAND06'  )     
             READ( LUN, 10 ) GC_ILAND06
          CASE( 'ILAND07'  )     
             READ( LUN, 10 ) GC_ILAND07
          CASE( 'ILAND08'  )     
             READ( LUN, 10 ) GC_ILAND08
          CASE( 'ILAND09'  )     
             READ( LUN, 10 ) GC_ILAND09
          CASE( 'ILAND10'  )     
             READ( LUN, 10 ) GC_ILAND10
          CASE( 'ILAND11'  )     
             READ( LUN, 10 ) GC_ILAND11
          CASE( 'ILAND12'  )     
             READ( LUN, 10 ) GC_ILAND12
          CASE( 'ILAND13'  )     
             READ( LUN, 10 ) GC_ILAND13
          CASE( 'ILAND14'  )     
             READ( LUN, 10 ) GC_ILAND14
          CASE( 'ILAND15'  )     
             READ( LUN, 10 ) GC_ILAND15
          CASE( 'IUSE01'  )     
             READ( LUN, 10 ) GC_IUSE01
          CASE( 'IUSE02'  )     
             READ( LUN, 10 ) GC_IUSE02
          CASE( 'IUSE03'  )     
             READ( LUN, 10 ) GC_IUSE03
          CASE( 'IUSE04' )     
             READ( LUN, 10 ) GC_IUSE04
          CASE( 'IUSE05' )     
             READ( LUN, 10 ) GC_IUSE05
          CASE( 'IUSE06' )     
             READ( LUN, 10 ) GC_IUSE06
          CASE( 'IUSE07' )     
             READ( LUN, 10 ) GC_IUSE07
          CASE( 'IUSE08' )     
             READ( LUN, 10 ) GC_IUSE08
          CASE( 'IUSE09' )     
             READ( LUN, 10 ) GC_IUSE09
          CASE( 'IUSE10')     
             READ( LUN, 10 ) GC_IUSE10
          CASE( 'IUSE11' )     
             READ( LUN, 10 ) GC_IUSE11
          CASE( 'IUSE12' )     
             READ( LUN, 10 ) GC_IUSE12
          CASE( 'IUSE13' )     
             READ( LUN, 10 ) GC_IUSE13
          CASE( 'IUSE14' )     
             READ( LUN, 10 ) GC_IUSE14
          CASE( 'IUSE15' )     
             READ( LUN, 10 ) GC_IUSE15
          CASE( 'LAI01' )      
             READ( LUN, 10 ) GC_LAI01
          CASE( 'LAI02' )     
             READ( LUN, 10 ) GC_LAI02
          CASE( 'LAI03' )     
             READ( LUN, 10 ) GC_LAI03
          CASE( 'LAI04' )     
             READ( LUN, 10 ) GC_LAI04
          CASE( 'LAI05' )     
             READ( LUN, 10 ) GC_LAI05
          CASE( 'LAI06' )     
             READ( LUN, 10 ) GC_LAI06
          CASE( 'LAI07' )     
             READ( LUN, 10 ) GC_LAI07
          CASE( 'LAI08' )     
             READ( LUN, 10 ) GC_LAI08
          CASE( 'LAI09' )     
             READ( LUN, 10 ) GC_LAI09
          CASE( 'LAI10' )     
             READ( LUN, 10 ) GC_LAI10
          CASE( 'LAI11' )     
             READ( LUN, 10 ) GC_LAI11
          CASE( 'LAI12' )     
             READ( LUN, 10 ) GC_LAI12
          CASE( 'LAI13' )     
             READ( LUN, 10 ) GC_LAI13
          CASE( 'LAI14' )     
             READ( LUN, 10 ) GC_LAI14
          CASE( 'LAI15' )     
             READ( LUN, 10 ) GC_LAI15           

          ! Emission data (from SETEMIS)
          CASE( 'EMISSION_NOx', 'EMISSION_NOX' ) 
             READ( LUN, 10 ) GC_EMISS_NOX
         !CASE( 'EMISSION_O3' )    
         !   READ( LUN, 10 ) GC_EMISS_O3
          CASE( 'EMISSION_CO' )    
             READ( LUN, 10 ) GC_EMISS_CO
          CASE( 'EMISSION_ALK4' )  
             READ( LUN, 10 ) GC_EMISS_ALK4
          CASE( 'EMISSION_ISOP  ' )  
             READ( LUN, 10 ) GC_EMISS_ISOP
          CASE( 'EMISSION_HNO3' )  
             READ( LUN, 10 ) GC_EMISS_HNO3
          CASE( 'EMISSION_ACET' )  
             READ( LUN, 10 ) GC_EMISS_ACET
          CASE( 'EMISSION_MEK' )   
             READ( LUN, 10 ) GC_EMISS_MEK
          CASE( 'EMISSION_ALD2' )  
             READ( LUN, 10 ) GC_EMISS_ALD2
          CASE( 'EMISSION_PRPE' )  
             READ( LUN, 10 ) GC_EMISS_PRPE
          CASE( 'EMISSION_C3H8' )  
             READ( LUN, 10 ) GC_EMISS_C3H8
          CASE( 'EMISSION_CH2O' )  
             READ( LUN, 10 ) GC_EMISS_CH2O
          CASE( 'EMISSION_C2H6' )  
             READ( LUN, 10 ) GC_EMISS_C2H6

          ! Strat chemistry
          CASE( 'SCHEM_OH' )       
             READ( LUN, 10 ) GC_SOX_OH
          CASE( 'SCHEM_PCO' )   
             READ( LUN, 10 ) GC_SOX_PCO
          CASE( 'SCHEM_LCO' )      
             READ( LUN, 10 ) GC_SOX_LCO
          CASE( 'SCHEM_NOX' )  
             READ( LUN, 10 ) GC_SOX_JV_NOX
          CASE( 'SCHEM_H2O2' )  
             READ( LUN, 10 ) GC_SOX_JV_H2O2
          CASE( 'SCHEM_ACET' )  
             READ( LUN, 10 ) GC_SOX_JV_ACET
          CASE( 'SCHEM_MEK' )   
             READ( LUN, 10 ) GC_SOX_JV_MEK
          CASE( 'SCHEM_ALD2' )  
             READ( LUN, 10 ) GC_SOX_JV_ALD2
          CASE( 'SCHEM_RCHO' )  
             READ( LUN, 10 ) GC_SOX_JV_RCHO
          CASE( 'SCHEM_MVK' )   
             READ( LUN, 10 ) GC_SOX_JV_MVK
          CASE( 'SCHEM_MACR' )  
             READ( LUN, 10 ) GC_SOX_JV_MACR
          CASE( 'SCHEM_R4N2' )  
             READ( LUN, 10 ) GC_SOX_JV_R4N2
          CASE( 'SCHEM_CH2O' )  
             READ( LUN, 10 ) GC_SOX_JV_CH2O
          CASE( 'SCHEM_N2O5' )  
             READ( LUN, 10 ) GC_SOX_JV_N2O5
          CASE( 'SCHEM_HNO4' )  
             READ( LUN, 10 ) GC_SOX_JV_HNO4
          CASE( 'SCHEM_MP' )     
             READ( LUN, 10 ) GC_SOX_JV_MP
          CASE DEFAULT
             ! Nothing
       END SELECT

    ENDDO

    ! Close file
    CLOSE( LUN )

    !========================================================================
    ! Fill typical values (or read them from disk and assign them here)
    !========================================================================
    DO J = J1, J2
    DO I = I1, I2

       !-----------------------------------------------------
       ! Met fields
       !-----------------------------------------------------

       ! Surface fields
       ALBEDO         (I,J  ) = GC_ALBD
       CLDFRC         (I,J  ) = GC_CLDFRC
       SH             (I,J  ) = GC_HFLUX
       FRLAND         (I,J  ) = GC_FRCLND
       LWI            (I,J  ) = GC_LWI
       ZPBL           (I,J  ) = GC_ZPBL
       CN_PRCP        (I,J  ) = GC_PRECCON
       LS_PRCP        (I,J  ) = GC_PRECTOT-GC_PRECCON
       AN_PRCP        (I,J  ) = 0.0
       TROPP          (I,J  ) = GC_TROPP
       T2M            (I,J  ) = GC_TS
       RADSRF         (I,J  ) = GC_RADSWG
       TSEA           (I,J  ) = GC_SST
       U10M           (I,J  ) = GC_U10M
       USTAR          (I,J  ) = GC_USTAR
       V10M           (I,J  ) = GC_V10M
       Z0             (I,J  ) = GC_Z0
       COSZ           (I,J  ) = GC_COSZ
       TO3            (I,J  ) = GC_TO3
       UVALBEDO       (I,J  ) = GC_UVALBEDO

       ! Column fields
       FCLDF          (I,J,:) = GC_CLDF       (:)
       CNV_FMC        (I,J,:) = GC_CMFMC      (:)
       CNV_MFD        (I,J,:) = GC_DTRAIN     (:)
       TAUCLI         (I,J,:) = GC_TAUCLI     (:)
       TAUCLW         (I,J,:) = GC_TAUCLW     (:)
       RH2            (I,J,:) = GC_RH         (:)
       Q              (I,J,:) = GC_Q          (:)
       T              (I,J,:) = GC_T          (:)
       AIRDENS        (I,J,:) = GC_AIRDENS    (:)
       MOISTQ         (I,J,:) = GC_MOISTQ     (:)
       DELP           (I,J,:) = GC_DELP       (:)
       PLE            (I,J,:) = GC_PEDGE      (:)
       PL             (I,J,:) = GC_PMID       (:)

       !-----------------------------------------------------
       ! Land type & LAI info
       !-----------------------------------------------------
       IREG           (I,J  ) = GC_IREG
       ILAND_01       (I,J  ) = GC_ILAND01
       ILAND_02       (I,J  ) = GC_ILAND02
       ILAND_03       (I,J  ) = GC_ILAND03
       ILAND_04       (I,J  ) = GC_ILAND04
       ILAND_05       (I,J  ) = GC_ILAND05
       ILAND_06       (I,J  ) = GC_ILAND06
       ILAND_07       (I,J  ) = GC_ILAND07
       ILAND_08       (I,J  ) = GC_ILAND08
       ILAND_09       (I,J  ) = GC_ILAND09
       ILAND_10       (I,J  ) = GC_ILAND10
       ILAND_11       (I,J  ) = GC_ILAND11
       ILAND_12       (I,J  ) = GC_ILAND12
       ILAND_13       (I,J  ) = GC_ILAND13
       ILAND_14       (I,J  ) = GC_ILAND14
       ILAND_15       (I,J  ) = GC_ILAND15
       IUSE_01        (I,J  ) = GC_IUSE01
       IUSE_02        (I,J  ) = GC_IUSE02
       IUSE_03        (I,J  ) = GC_IUSE03
       IUSE_04        (I,J  ) = GC_IUSE04
       IUSE_05        (I,J  ) = GC_IUSE05
       IUSE_06        (I,J  ) = GC_IUSE06
       IUSE_07        (I,J  ) = GC_IUSE07
       IUSE_08        (I,J  ) = GC_IUSE08
       IUSE_09        (I,J  ) = GC_IUSE09
       IUSE_10        (I,J  ) = GC_IUSE10
       IUSE_11        (I,J  ) = GC_IUSE11
       IUSE_12        (I,J  ) = GC_IUSE12
       IUSE_13        (I,J  ) = GC_IUSE13
       IUSE_14        (I,J  ) = GC_IUSE14
       IUSE_15        (I,J  ) = GC_IUSE15
       LAI_01         (I,J  ) = GC_LAI01
       LAI_02         (I,J  ) = GC_LAI02
       LAI_03         (I,J  ) = GC_LAI03
       LAI_04         (I,J  ) = GC_LAI04
       LAI_05         (I,J  ) = GC_LAI05
       LAI_06         (I,J  ) = GC_LAI06
       LAI_07         (I,J  ) = GC_LAI07
       LAI_08         (I,J  ) = GC_LAI08
       LAI_09         (I,J  ) = GC_LAI09
       LAI_10         (I,J  ) = GC_LAI10
       LAI_11         (I,J  ) = GC_LAI11
       LAI_12         (I,J  ) = GC_LAI12
       LAI_13         (I,J  ) = GC_LAI13
       LAI_14         (I,J  ) = GC_LAI14
       LAI_15         (I,J  ) = GC_LAI15

       !-----------------------------------------------------
       ! Emissions
       !-----------------------------------------------------
       EMISS_NOX      (I,J,:) = GC_EMISS_NOx  (:)
       EMISS_CO       (I,J,:) = GC_EMISS_CO   (:)
       EMISS_ALK4     (I,J,:) = GC_EMISS_ALK4 (:)
       EMISS_ISOP     (I,J,:) = GC_EMISS_ISOP (:)
       EMISS_HNO3     (I,J,:) = GC_EMISS_HNO3 (:)
       EMISS_ACET     (I,J,:) = GC_EMISS_ACET (:)
       EMISS_MEK      (I,J,:) = GC_EMISS_MEK  (:)
       EMISS_ALD2     (I,J,:) = GC_EMISS_ALD2 (:)
       EMISS_PRPE     (I,J,:) = GC_EMISS_PRPE (:)
       EMISS_C3H8     (I,J,:) = GC_EMISS_C3H8 (:)
       EMISS_CH2O     (I,J,:) = GC_EMISS_CH2O (:)
       EMISS_C2H6     (I,J,:) = GC_EMISS_C2H6 (:)
       !----------------------------------------------------------------------
       ! NOTE: For now the GEOS-Chem column code adds emissions for aerosol
       ! tracers directly into the tracer array.  So leave these zeroed out
       ! for now.  At some future date we will fix it so that the emissions
       ! for aerosol tracers go into these arrays. (bmy, 4/6/10)
       EMISS_O3       (I,J,:) = 0.0
       EMISS_DMS      (I,J,:) = 0.0
       EMISS_SO2      (I,J,:) = 0.0
       EMISS_SO4      (I,J,:) = 0.0
       EMISS_NH3      (I,J,:) = 0.0
       EMISS_BCPI     (I,J,:) = 0.0
       EMISS_OCPI     (I,J,:) = 0.0
       EMISS_BCPO     (I,J,:) = 0.0
       EMISS_OCPO     (I,J,:) = 0.0
       EMISS_ALPH     (I,J,:) = 0.0
       EMISS_LIMO     (I,J,:) = 0.0
       EMISS_ALCO     (I,J,:) = 0.0
       EMISS_DST1     (I,J,:) = 0.0
       EMISS_DST2     (I,J,:) = 0.0
       EMISS_DST3     (I,J,:) = 0.0
       EMISS_DST4     (I,J,:) = 0.0
       EMISS_SALA     (I,J,:) = 0.0
       EMISS_SALC     (I,J,:) = 0.0
       !----------------------------------------------------------------------

       !-----------------------------------------------------
       ! Stratospheric chemistry
       !-----------------------------------------------------
       SOX_OH         (I,J,:) = GC_SOX_OH     (:)
       SOX_JV_NOx     (I,J,:) = GC_SOX_JV_NOx (:)
       SOX_JV_H2O2    (I,J,:) = GC_SOX_JV_H2O2(:)
       SOX_JV_ACET    (I,J,:) = GC_SOX_JV_ACET(:)
       SOX_JV_MEK     (I,J,:) = GC_SOX_JV_MEK (:)
       SOX_JV_ALD2    (I,J,:) = GC_SOX_JV_ALD2(:)
       SOX_JV_RCHO    (I,J,:) = GC_SOX_JV_RCHO(:)
       SOX_JV_MVK     (I,J,:) = GC_SOX_JV_MVK (:)
       SOX_JV_MACR    (I,J,:) = GC_SOX_JV_MACR(:)
       SOX_JV_R4N2    (I,J,:) = GC_SOX_JV_R4N2(:)
       SOX_JV_CH2O    (I,J,:) = GC_SOX_JV_CH2O(:)
       SOX_JV_N2O5    (I,J,:) = GC_SOX_JV_N2O5(:)
       SOX_JV_HNO4    (I,J,:) = GC_SOX_JV_HNO4(:)
       SOX_JV_MP      (I,J,:) = GC_SOX_JV_MP  (:)
       SOX_PCO        (I,J,:) = GC_SOX_PCO    (:)
       SOX_LCO        (I,J,:) = GC_SOX_LCO    (:)
    ENDDO
    ENDDO

    ! Return successfully
    RC = 0

    END SUBROUTINE Fill_Import_State_
!!EOC
!!------------------------------------------------------------------------------
!!     NASA/GSFC, Global Modeling and Assimilation Office, Code 910.1 and      !
!!          Harvard University Atmospheric Chemistry Modeling Group            !
!!------------------------------------------------------------------------------
!!BOP
!!
!! !IROUTINE: Fill_Import_State_
!!
!! !DESCRIPTION: Subroutine Fill_Import_State_ populates the import state of the
!!  GEOS-Chem gridded component with reasonable values for the unit test.
!!\\
!!\\
!!  Normally, the GEOS-Chem gridded component would recieve these fields
!!  from the GEOS-5 GCM via the import state.  However, to simplify the unit
!!  testing, we instead will declare all import state variables here.
!!
!! !INTERFACE:
!!
!  SUBROUTINE Fill_Export_State_( Export, RC )
!!
!! !USES:
!!
!#   include "GEOSCHEMchem_DeclarePointer___.h"     ! Pointer decl's to States
!!
!! !INPUT/OUTPUT PARAMETERS:
!!
!    TYPE(ESMF_State), INTENT(INOUT) :: Export      ! Export state
!!
!! !OUTPUT PARAMETERS:
!! 
!    INTEGER,          INTENT(OUT)   :: RC          ! Return code
!!
!! !REMARKS:
!!  For simplicity, we assign data from a single column to the entire world.
!!  Here we are mostly concerned not w/ the performace of the GEOS-Chem code
!!  but to see if it will link up with the GEOS-5 GCM.  
!!
!! !REVISION HISTORY: 
!!   6 Jul 2010 - R. Yantosca - Initial version
!!EOP
!!------------------------------------------------------------------------------
!!BOC
!!
!! !LOCAL VARIABLES:
!!
!    INTEGER            :: I,  J    ! Loop indices
!    INTEGER            :: I1, I2   ! Min & max lon indices on this PET
!    INTEGER            :: J1, J2   ! Min & max lat indices on this PET
!    INTEGER            :: K1, K2   ! Min & max lev indices on this PET
!    INTEGER            :: IOS
!    CHARACTER(LEN=255) :: NAME
!
!    !========================================================================
!    ! Get Pointers to Export state fields
!    !========================================================================
!    CALL MAPL_GetPointer ( Export,  WD_HNO3,     'WD_HNO3',     __RC__ )
!    CALL MAPL_GetPointer ( Export,  WD_H2O2,     'WD_H2O2',     __RC__ )
!    CALL MAPL_GetPointer ( Export,  WD_CH2O,     'WD_CH2O',     __RC__ )
!    CALL MAPL_GetPointer ( Export,  WD_MP,       'WD_MP',       __RC__ )
!    CALL MAPL_GetPointer ( Export,  WD_SO2,      'WD_SO2',      __RC__ )
!    CALL MAPL_GetPointer ( Export,  WD_SO4,      'WD_SO4',      __RC__ )
!    CALL MAPL_GetPointer ( Export,  WD_SO4s,     'WD_SO4s',     __RC__ )
!    CALL MAPL_GetPointer ( Export,  WD_MSA,      'WD_MSA',      __RC__ )
!    CALL MAPL_GetPointer ( Export,  WD_NH3,      'WD_NH3',      __RC__ )
!    CALL MAPL_GetPointer ( Export,  WD_NH4,      'WD_NH4',      __RC__ )
!    CALL MAPL_GetPointer ( Export,  WD_NIT,      'WD_NIT',      __RC__ )
!    CALL MAPL_GetPointer ( Export,  WD_NITs,     'WD_NITs',     __RC__ )
!    CALL MAPL_GetPointer ( Export,  WD_BCPI,     'WD_BCPI',     __RC__ )
!    CALL MAPL_GetPointer ( Export,  WD_OCPI,     'WD_OCPI',     __RC__ )
!    CALL MAPL_GetPointer ( Export,  WD_BCPO,     'WD_BCPO',     __RC__ )
!    CALL MAPL_GetPointer ( Export,  WD_OCPO,     'WD_OCPO',     __RC__ )
!    CALL MAPL_GetPointer ( Export,  WD_ALPH,     'WD_ALPH',     __RC__ )
!    CALL MAPL_GetPointer ( Export,  WD_LIMO,     'WD_LIMO',     __RC__ )
!    CALL MAPL_GetPointer ( Export,  WD_ALCO,     'WD_ALCO',     __RC__ )
!    CALL MAPL_GetPointer ( Export,  WD_SOG1,     'WD_SOG1',     __RC__ )
!    CALL MAPL_GetPointer ( Export,  WD_SOG2,     'WD_SOG2',     __RC__ )
!    CALL MAPL_GetPointer ( Export,  WD_SOG3,     'WD_SOG3',     __RC__ )
!    CALL MAPL_GetPointer ( Export,  WD_SOG4,     'WD_SOG4',     __RC__ )
!    CALL MAPL_GetPointer ( Export,  WD_SOA1,     'WD_SOA1',     __RC__ )
!    CALL MAPL_GetPointer ( Export,  WD_SOA2,     'WD_SOA2',     __RC__ )
!    CALL MAPL_GetPointer ( Export,  WD_SOA3,     'WD_SOA3',     __RC__ )
!    CALL MAPL_GetPointer ( Export,  WD_SOA4,     'WD_SOA4',     __RC__ )
!    CALL MAPL_GetPointer ( Export,  WD_DST1,     'WD_DST1',     __RC__ )
!    CALL MAPL_GetPointer ( Export,  WD_DST2,     'WD_DST2',     __RC__ )
!    CALL MAPL_GetPointer ( Export,  WD_DST3,     'WD_DST3',     __RC__ )
!    CALL MAPL_GetPointer ( Export,  WD_DST4,     'WD_DST4',     __RC__ )
!    CALL MAPL_GetPointer ( Export,  WD_SALA,     'WD_SALA',     __RC__ )
!    CALL MAPL_GetPointer ( Export,  WD_SALC,     'WD_SALC',     __RC__ )
!    CALL MAPL_GetPointer ( Export,  DD_FX_NOx,   'DD_FX_NOx',   __RC__ )
!    CALL MAPL_GetPointer ( Export,  DD_FX_O3,    'DD_FX_O3',    __RC__ )
!    CALL MAPL_GetPointer ( Export,  DD_FX_PAN,   'DD_FX_PAN',   __RC__ )
!    CALL MAPL_GetPointer ( Export,  DD_FX_HNO3,  'DD_FX_HNO3',  __RC__ )
!    CALL MAPL_GetPointer ( Export,  DD_FX_H2O2,  'DD_FX_H2O2',  __RC__ )
!    CALL MAPL_GetPointer ( Export,  DD_FX_N2O5,  'DD_FX_N2O5',  __RC__ )
!    CALL MAPL_GetPointer ( Export,  DD_FX_PPN,   'DD_FX_PPN',   __RC__ )
!    CALL MAPL_GetPointer ( Export,  DD_FX_R4N2,  'DD_FX_R4N2',  __RC__ )
!    CALL MAPL_GetPointer ( Export,  DD_FX_CH2O,  'DD_FX_CH2O',  __RC__ )
!    CALL MAPL_GetPointer ( Export,  DD_FX_SO2,   'DD_FX_SO2',   __RC__ )
!    CALL MAPL_GetPointer ( Export,  DD_FX_SO4,   'DD_FX_SO4',   __RC__ )
!    CALL MAPL_GetPointer ( Export,  DD_FX_SO4S,  'DD_FX_SO4S',  __RC__ )
!    CALL MAPL_GetPointer ( Export,  DD_FX_MSA,   'DD_FX_MSA',   __RC__ )
!    CALL MAPL_GetPointer ( Export,  DD_FX_NH3,   'DD_FX_NH3',   __RC__ )
!    CALL MAPL_GetPointer ( Export,  DD_FX_NH4,   'DD_FX_NH4',   __RC__ )
!    CALL MAPL_GetPointer ( Export,  DD_FX_NIT,   'DD_FX_NIT',   __RC__ )
!    CALL MAPL_GetPointer ( Export,  DD_FX_NITS,  'DD_FX_NITS',  __RC__ )
!    CALL MAPL_GetPointer ( Export,  DD_FX_BCPI,  'DD_FX_BCPI',  __RC__ )
!    CALL MAPL_GetPointer ( Export,  DD_FX_OCPI,  'DD_FX_OCPI',  __RC__ )
!    CALL MAPL_GetPointer ( Export,  DD_FX_BCPO,  'DD_FX_BCPO',  __RC__ )
!    CALL MAPL_GetPointer ( Export,  DD_FX_OCPO,  'DD_FX_OCPO',  __RC__ )
!    CALL MAPL_GetPointer ( Export,  DD_FX_ALPH,  'DD_FX_ALPH',  __RC__ )
!    CALL MAPL_GetPointer ( Export,  DD_FX_LIMO,  'DD_FX_LIMO',  __RC__ )
!    CALL MAPL_GetPointer ( Export,  DD_FX_ALCO,  'DD_FX_ALCO',  __RC__ )
!    CALL MAPL_GetPointer ( Export,  DD_FX_SOG1,  'DD_FX_SOG1',  __RC__ )
!    CALL MAPL_GetPointer ( Export,  DD_FX_SOG2,  'DD_FX_SOG2',  __RC__ )
!    CALL MAPL_GetPointer ( Export,  DD_FX_SOG3,  'DD_FX_SOG3',  __RC__ )
!    CALL MAPL_GetPointer ( Export,  DD_FX_SOG4,  'DD_FX_SOG4',  __RC__ )
!    CALL MAPL_GetPointer ( Export,  DD_FX_SOA1,  'DD_FX_SOA1',  __RC__ )
!    CALL MAPL_GetPointer ( Export,  DD_FX_SOA2,  'DD_FX_SOA2',  __RC__ )
!    CALL MAPL_GetPointer ( Export,  DD_FX_SOA3,  'DD_FX_SOA3',  __RC__ )
!    CALL MAPL_GetPointer ( Export,  DD_FX_SOA4,  'DD_FX_SOA4',  __RC__ )
!    CALL MAPL_GetPointer ( Export,  DD_FX_DST1,  'DD_FX_DST1',  __RC__ )
!    CALL MAPL_GetPointer ( Export,  DD_FX_DST2,  'DD_FX_DST2',  __RC__ )
!    CALL MAPL_GetPointer ( Export,  DD_FX_DST3,  'DD_FX_DST3',  __RC__ )
!    CALL MAPL_GetPointer ( Export,  DD_FX_DST4,  'DD_FX_DST4',  __RC__ )
!    CALL MAPL_GetPointer ( Export,  DD_FX_SALA,  'DD_FX_SALA',  __RC__ )
!    CALL MAPL_GetPointer ( Export,  DD_FX_SALC,  'DD_FX_SALC',  __RC__ )
!    CALL MAPL_GetPointer ( Export,  DD_FQ_NOx,   'DD_FQ_NOx',   __RC__ )
!    CALL MAPL_GetPointer ( Export,  DD_FQ_O3,    'DD_FQ_O3',    __RC__ )
!    CALL MAPL_GetPointer ( Export,  DD_FQ_PAN,   'DD_FQ_PAN',   __RC__ )
!    CALL MAPL_GetPointer ( Export,  DD_FQ_HNO3,  'DD_FQ_HNO3',  __RC__ )
!    CALL MAPL_GetPointer ( Export,  DD_FQ_H2O2,  'DD_FQ_H2O2',  __RC__ )
!    CALL MAPL_GetPointer ( Export,  DD_FQ_N2O5,  'DD_FQ_N2O5',  __RC__ )
!    CALL MAPL_GetPointer ( Export,  DD_FQ_PPN,   'DD_FQ_PPN',   __RC__ )
!    CALL MAPL_GetPointer ( Export,  DD_FQ_R4N2,  'DD_FQ_R4N2',  __RC__ )
!    CALL MAPL_GetPointer ( Export,  DD_FQ_CH2O,  'DD_FQ_CH2O',  __RC__ )
!    CALL MAPL_GetPointer ( Export,  DD_FQ_SO2,   'DD_FQ_SO2',   __RC__ )
!    CALL MAPL_GetPointer ( Export,  DD_FQ_SO4,   'DD_FQ_SO4',   __RC__ )
!    CALL MAPL_GetPointer ( Export,  DD_FQ_SO4S,  'DD_FQ_SO4S',  __RC__ )
!    CALL MAPL_GetPointer ( Export,  DD_FQ_MSA,   'DD_FQ_MSA',   __RC__ )
!    CALL MAPL_GetPointer ( Export,  DD_FQ_NH3,   'DD_FQ_NH3',   __RC__ )
!    CALL MAPL_GetPointer ( Export,  DD_FQ_NH4,   'DD_FQ_NH4',   __RC__ )
!    CALL MAPL_GetPointer ( Export,  DD_FQ_NIT,   'DD_FQ_NIT',   __RC__ )
!    CALL MAPL_GetPointer ( Export,  DD_FQ_NITS,  'DD_FQ_NITS',  __RC__ )
!    CALL MAPL_GetPointer ( Export,  DD_FQ_BCPI,  'DD_FQ_BCPI',  __RC__ )
!    CALL MAPL_GetPointer ( Export,  DD_FQ_OCPI,  'DD_FQ_OCPI',  __RC__ )
!    CALL MAPL_GetPointer ( Export,  DD_FQ_BCPO,  'DD_FQ_BCPO',  __RC__ )
!    CALL MAPL_GetPointer ( Export,  DD_FQ_OCPO,  'DD_FQ_OCPO',  __RC__ )
!    CALL MAPL_GetPointer ( Export,  DD_FQ_ALPH,  'DD_FQ_ALPH',  __RC__ )
!    CALL MAPL_GetPointer ( Export,  DD_FQ_LIMO,  'DD_FQ_LIMO',  __RC__ )
!    CALL MAPL_GetPointer ( Export,  DD_FQ_ALCO,  'DD_FQ_ALCO',  __RC__ )
!    CALL MAPL_GetPointer ( Export,  DD_FQ_SOG1,  'DD_FQ_SOG1',  __RC__ )
!    CALL MAPL_GetPointer ( Export,  DD_FQ_SOG2,  'DD_FQ_SOG2',  __RC__ )
!    CALL MAPL_GetPointer ( Export,  DD_FQ_SOG3,  'DD_FQ_SOG3',  __RC__ )
!    CALL MAPL_GetPointer ( Export,  DD_FQ_SOG4,  'DD_FQ_SOG4',  __RC__ )
!    CALL MAPL_GetPointer ( Export,  DD_FQ_SOA1,  'DD_FQ_SOA1',  __RC__ )
!    CALL MAPL_GetPointer ( Export,  DD_FQ_SOA2,  'DD_FQ_SOA2',  __RC__ )
!    CALL MAPL_GetPointer ( Export,  DD_FQ_SOA3,  'DD_FQ_SOA3',  __RC__ )
!    CALL MAPL_GetPointer ( Export,  DD_FQ_SOA4,  'DD_FQ_SOA4',  __RC__ )
!    CALL MAPL_GetPointer ( Export,  DD_FQ_DST1,  'DD_FQ_DST1',  __RC__ )
!    CALL MAPL_GetPointer ( Export,  DD_FQ_DST2,  'DD_FQ_DST2',  __RC__ )
!    CALL MAPL_GetPointer ( Export,  DD_FQ_DST3,  'DD_FQ_DST3',  __RC__ )
!    CALL MAPL_GetPointer ( Export,  DD_FQ_DST4,  'DD_FQ_DST4',  __RC__ )
!    CALL MAPL_GetPointer ( Export,  DD_FQ_SALA,  'DD_FQ_SALA',  __RC__ )
!    CALL MAPL_GetPointer ( Export,  DD_FQ_SALC,  'DD_FQ_SALC',  __RC__ )
!    CALL MAPL_GetPointer ( Export,  DD_V_NOx,    'DD_V_NOx',    __RC__ )
!    CALL MAPL_GetPointer ( Export,  DD_V_O3,     'DD_V_O3',     __RC__ )
!    CALL MAPL_GetPointer ( Export,  DD_V_PAN,    'DD_V_PAN',    __RC__ )
!    CALL MAPL_GetPointer ( Export,  DD_V_HNO3,   'DD_V_HNO3',   __RC__ )
!    CALL MAPL_GetPointer ( Export,  DD_V_H2O2,   'DD_V_H2O2',   __RC__ )
!    CALL MAPL_GetPointer ( Export,  DD_V_N2O5,   'DD_V_N2O5',   __RC__ )
!    CALL MAPL_GetPointer ( Export,  DD_V_PPN,    'DD_V_PPN',    __RC__ )
!    CALL MAPL_GetPointer ( Export,  DD_V_R4N2,   'DD_V_R4N2',   __RC__ )
!    CALL MAPL_GetPointer ( Export,  DD_V_CH2O,   'DD_V_CH2O',   __RC__ )
!    CALL MAPL_GetPointer ( Export,  DD_V_SO2,    'DD_V_SO2',    __RC__ )
!    CALL MAPL_GetPointer ( Export,  DD_V_SO4,    'DD_V_SO4',    __RC__ )
!    CALL MAPL_GetPointer ( Export,  DD_V_SO4S,   'DD_V_SO4S',   __RC__ )
!    CALL MAPL_GetPointer ( Export,  DD_V_MSA,    'DD_V_MSA',    __RC__ )
!    CALL MAPL_GetPointer ( Export,  DD_V_NH3,    'DD_V_NH3',    __RC__ )
!    CALL MAPL_GetPointer ( Export,  DD_V_NH4,    'DD_V_NH4',    __RC__ )
!    CALL MAPL_GetPointer ( Export,  DD_V_NIT,    'DD_V_NIT',    __RC__ )
!    CALL MAPL_GetPointer ( Export,  DD_V_NITS,   'DD_V_NITS',   __RC__ )
!    CALL MAPL_GetPointer ( Export,  DD_V_BCPI,   'DD_V_BCPI',   __RC__ )
!    CALL MAPL_GetPointer ( Export,  DD_V_OCPI,   'DD_V_OCPI',   __RC__ )
!    CALL MAPL_GetPointer ( Export,  DD_V_BCPO,   'DD_V_BCPO',   __RC__ )
!    CALL MAPL_GetPointer ( Export,  DD_V_OCPO,   'DD_V_OCPO',   __RC__ )
!    CALL MAPL_GetPointer ( Export,  DD_V_ALPH,   'DD_V_ALPH',   __RC__ )
!    CALL MAPL_GetPointer ( Export,  DD_V_LIMO,   'DD_V_LIMO',   __RC__ )
!    CALL MAPL_GetPointer ( Export,  DD_V_ALCO,   'DD_V_ALCO',   __RC__ )
!    CALL MAPL_GetPointer ( Export,  DD_V_SOG1,   'DD_V_SOG1',   __RC__ )
!    CALL MAPL_GetPointer ( Export,  DD_V_SOG2,   'DD_V_SOG2',   __RC__ )
!    CALL MAPL_GetPointer ( Export,  DD_V_SOG3,   'DD_V_SOG3',   __RC__ )
!    CALL MAPL_GetPointer ( Export,  DD_V_SOG4,   'DD_V_SOG4',   __RC__ )
!    CALL MAPL_GetPointer ( Export,  DD_V_SOA1,   'DD_V_SOA1',   __RC__ )
!    CALL MAPL_GetPointer ( Export,  DD_V_SOA2,   'DD_V_SOA2',   __RC__ )
!    CALL MAPL_GetPointer ( Export,  DD_V_SOA3,   'DD_V_SOA3',   __RC__ )
!    CALL MAPL_GetPointer ( Export,  DD_V_SOA4,   'DD_V_SOA4',   __RC__ )
!    CALL MAPL_GetPointer ( Export,  DD_V_DST1,   'DD_V_DST1',   __RC__ )
!    CALL MAPL_GetPointer ( Export,  DD_V_DST2,   'DD_V_DST2',   __RC__ )
!    CALL MAPL_GetPointer ( Export,  DD_V_DST3,   'DD_V_DST3',   __RC__ )
!    CALL MAPL_GetPointer ( Export,  DD_V_DST4,   'DD_V_DST4',   __RC__ )
!    CALL MAPL_GetPointer ( Export,  DD_V_SALA,   'DD_V_SALA',   __RC__ )
!    CALL MAPL_GetPointer ( Export,  DD_V_SALC,   'DD_V_SALC',   __RC__ )
!
!!    ! Make sure the pointer is allocated
!!    ASSERT_( ASSOCIATED( WD_HNO3 ) )
!!
!!    ! Get indices from a 3-D data block
!!    I1 = LBOUND( WD_HNO3, 1 )
!!    J1 = LBOUND( WD_HNO3, 2 )
!!    K1 = LBOUND( WD_HNO3, 3 )
!!    I2 = UBOUND( WD_HNO3, 1 )
!!    J2 = UBOUND( WD_HNO3, 2 )
!!    K2 = UBOUND( WD_HNO3, 3 )
!!
!!    ! Make sure the horizontal dimensions of the array are OK
!!    ASSERT_( I1 > 0 )       
!!    ASSERT_( I2 > 0 )
!!    ASSERT_( J1 > 0 )
!!    ASSERT_( J2 > 0 )
!!
!!    ! Make sure we have 72 levels
!!    ASSERT_( ( K2 - K1 + 1 ) == 72 )
!
!
!#if defined( DEBUG_PRINT )
!    print*, '### a WD_HNO3: ', associated( WD_HNO3 )
!    print*, '### a WD_H2O2: ', associated( WD_H2O2 )
!#endif
!
!    END SUBROUTINE Fill_Export_State_
!EOC
!------------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 910.1 and      !
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Fill_Internal_State_
!
! !DESCRIPTION: Subroutine Fill_Import_State_ populates the internal state of 
!  the GEOS-Chem gridded component with reasonable values for the unit test.
!\\
!\\
!  Normally, the GEOS-Chem gridded component would read the initial conditions
!  for internal state variables from a restart file.  However, to simplify
!  the unit testing, we instead will declare all internal state variables
!  here.
!
! !INTERFACE:
!
  SUBROUTINE Fill_Internal_State_( GrComp, inFile, RC )
!
! !USES:
!
#   include "GEOSCHEMchem_DeclarePointer___.h"    ! Pointer decl's to States
!
! !INPUT PARAMETERS
!
    CHARACTER(LEN=*),    INTENT(IN)    :: inFile  ! Filename
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ESMF_GridComp), INTENT(INOUT) :: GrComp  ! Gridded component
!
! !OUTPUT PARAMETERS:
! 
    INTEGER,             INTENT(OUT)   :: RC      ! Return code
!
! !REMARKS:
!  For simplicity, we assign data from a single column to the entire world.
!  Here we are mostly concerned not w/ the performace of the GEOS-Chem code
!  but to see if it will link up with the GEOS-5 GCM.  
!
! !REVISION HISTORY: 
!  01 Dec 2009 - A. Da Silva - Initial version  
!  02 Apr 2010 - R. Yantosca - Modified for GEOS-Chem column code
!  02 Apr 2010 - R. Yantosca - Added ProTex Headers, other cosmetic changes
!  07 Apr 2010 - R. Yantosca - Populate internal state fields with values
!                              saved from the GEOS_Chem model
!  07 Apr 2010 - R. Yantosca - Also populate fields for advected tracers
!  01 Jun 2010 - R. Yantosca - Changed output format to 4(es19.12,1x)
!  01 Jul 2010 - R. Yantosca - Now zero D_OH_MASS and D_AIR_MASS
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    TYPE(MAPL_MetaComp), POINTER :: MAPL             ! Pointer to MAPL
    INTEGER                      :: I1, I2           ! Min & max lon indices
    INTEGER                      :: J1, J2           ! Min & max lat indices
    INTEGER                      :: K1, K2           ! Min & max lev indices
    INTEGER                      :: IOS              ! I/O error code
    CHARACTER(LEN=255)           :: NAME             ! Field name in file
    REAL                         :: PEDGE(LMAX+1)    ! Placeholder
    REAL                         :: PMID (LMAX)      ! Placeholder
!
! !DEFINED PARAMETERS:
!
    INTEGER, PARAMETER :: LUN=200

    !========================================================================
    ! Get the internal state from the GEOS-Chem gridded component
    !========================================================================

    ! Get MAPL object attached to ESMF gridded component
    CALL MAPL_GetObjectFromGC( GrComp, MAPL, __RC__ )

    ! Get internal state from MAPL object
    CALL MAPL_Get( MAPL, INTERNAL_ESMF_STATE=Internal, __RC__ )

    !========================================================================
    ! Get Pointers to Internal state
    !========================================================================
    CALL MAPL_GetPointer( Internal, TRC_NOx,       'TRC_NOx',        __RC__ )
    CALL MAPL_GetPointer( Internal, TRC_Ox,        'TRC_Ox',         __RC__ )
    CALL MAPL_GetPointer( Internal, TRC_PAN,       'TRC_PAN',        __RC__ )
    CALL MAPL_GetPointer( Internal, TRC_CO,        'TRC_CO',         __RC__ )
    CALL MAPL_GetPointer( Internal, TRC_ALK4,      'TRC_ALK4',       __RC__ )
    CALL MAPL_GetPointer( Internal, TRC_ISOP,      'TRC_ISOP',       __RC__ )
    CALL MAPL_GetPointer( Internal, TRC_HNO3,      'TRC_HNO3',       __RC__ )
    CALL MAPL_GetPointer( Internal, TRC_H2O2,      'TRC_H2O2',       __RC__ )
    CALL MAPL_GetPointer( Internal, TRC_ACET,      'TRC_ACET',       __RC__ )
    CALL MAPL_GetPointer( Internal, TRC_MEK,       'TRC_MEK',        __RC__ )
    CALL MAPL_GetPointer( Internal, TRC_ALD2,      'TRC_ALD2',       __RC__ )
    CALL MAPL_GetPointer( Internal, TRC_RCHO,      'TRC_RCHO',       __RC__ )
    CALL MAPL_GetPointer( Internal, TRC_MVK,       'TRC_MVK',        __RC__ )
    CALL MAPL_GetPointer( Internal, TRC_MACR,      'TRC_MACR',       __RC__ )
    CALL MAPL_GetPointer( Internal, TRC_PMN,       'TRC_PMN',        __RC__ )
    CALL MAPL_GetPointer( Internal, TRC_PPN,       'TRC_PPN',        __RC__ )
    CALL MAPL_GetPointer( Internal, TRC_R4N2,      'TRC_R4N2',       __RC__ )
    CALL MAPL_GetPointer( Internal, TRC_PRPE,      'TRC_PRPE',       __RC__ )
    CALL MAPL_GetPointer( Internal, TRC_C3H8,      'TRC_C3H8',       __RC__ )
    CALL MAPL_GetPointer( Internal, TRC_CH2O,      'TRC_CH2O',       __RC__ )
    CALL MAPL_GetPointer( Internal, TRC_C2H6,      'TRC_C2H6',       __RC__ )
    CALL MAPL_GetPointer( Internal, TRC_N2O5,      'TRC_N2O5',       __RC__ )
    CALL MAPL_GetPointer( Internal, TRC_HNO4,      'TRC_HNO4',       __RC__ )
    CALL MAPL_GetPointer( Internal, TRC_MP,        'TRC_MP',         __RC__ )
    CALL MAPL_GetPointer( Internal, TRC_DMS,       'TRC_DMS',        __RC__ )
    CALL MAPL_GetPointer( Internal, TRC_SO2,       'TRC_SO2',        __RC__ )
    CALL MAPL_GetPointer( Internal, TRC_SO4,       'TRC_SO4',        __RC__ )
    CALL MAPL_GetPointer( Internal, TRC_SO4s,      'TRC_SO4s',       __RC__ )
    CALL MAPL_GetPointer( Internal, TRC_MSA,       'TRC_MSA',        __RC__ )
    CALL MAPL_GetPointer( Internal, TRC_NH3,       'TRC_NH3',        __RC__ )
    CALL MAPL_GetPointer( Internal, TRC_NH4,       'TRC_NH4',        __RC__ )
    CALL MAPL_GetPointer( Internal, TRC_NIT,       'TRC_NIT',        __RC__ )
    CALL MAPL_GetPointer( Internal, TRC_NITs,      'TRC_NITs',       __RC__ )
    CALL MAPL_GetPointer( Internal, TRC_BCPI,      'TRC_BCPI',       __RC__ )
    CALL MAPL_GetPointer( Internal, TRC_OCPI,      'TRC_OCPI',       __RC__ )
    CALL MAPL_GetPointer( Internal, TRC_BCPO,      'TRC_BCPO',       __RC__ )
    CALL MAPL_GetPointer( Internal, TRC_OCPO,      'TRC_OCPO',       __RC__ )
    CALL MAPL_GetPointer( Internal, TRC_ALPH,      'TRC_ALPH',       __RC__ )
    CALL MAPL_GetPointer( Internal, TRC_LIMO,      'TRC_LIMO',       __RC__ )
    CALL MAPL_GetPointer( Internal, TRC_ALCO,      'TRC_ALCO',       __RC__ )
    CALL MAPL_GetPointer( Internal, TRC_SOG1,      'TRC_SOG1',       __RC__ )
    CALL MAPL_GetPointer( Internal, TRC_SOG2,      'TRC_SOG2',       __RC__ )
    CALL MAPL_GetPointer( Internal, TRC_SOG3,      'TRC_SOG3',       __RC__ )
    CALL MAPL_GetPointer( Internal, TRC_SOG4,      'TRC_SOG4',       __RC__ )
    CALL MAPL_GetPointer( Internal, TRC_SOA1,      'TRC_SOA1',       __RC__ )
    CALL MAPL_GetPointer( Internal, TRC_SOA2,      'TRC_SOA2',       __RC__ )
    CALL MAPL_GetPointer( Internal, TRC_SOA3,      'TRC_SOA3',       __RC__ )
    CALL MAPL_GetPointer( Internal, TRC_SOA4,      'TRC_SOA4',       __RC__ )
    CALL MAPL_GetPointer( Internal, TRC_DST1,      'TRC_DST1',       __RC__ )
    CALL MAPL_GetPointer( Internal, TRC_DST2,      'TRC_DST2',       __RC__ )
    CALL MAPL_GetPointer( Internal, TRC_DST3,      'TRC_DST3',       __RC__ )
    CALL MAPL_GetPointer( Internal, TRC_DST4,      'TRC_DST4',       __RC__ )
    CALL MAPL_GetPointer( Internal, TRC_SALA,      'TRC_SALA',       __RC__ )
    CALL MAPL_GetPointer( Internal, TRC_SALC,      'TRC_SALC',       __RC__ )
    CALL MAPL_GetPointer( Internal, ORVC_TERP,     'ORVC_TERP',      __RC__ )
    CALL MAPL_GetPointer( Internal, ORVC_SESQ,     'ORVC_SESQ',      __RC__ )
    CALL MAPL_GetPointer( Internal, H2O2s,         'H2O2s',          __RC__ )
    CALL MAPL_GetPointer( Internal, SO2s,          'SO2s',           __RC__ )
    CALL MAPL_GetPointer( Internal, A3O2,          'A3O2',           __RC__ )
    CALL MAPL_GetPointer( Internal, ACET,          'ACET',           __RC__ )
    CALL MAPL_GetPointer( Internal, ACTA,          'ACTA',           __RC__ )
    CALL MAPL_GetPointer( Internal, ALD2,          'ALD2',           __RC__ )
    CALL MAPL_GetPointer( Internal, ALK4,          'ALK4',           __RC__ )
    CALL MAPL_GetPointer( Internal, ATO2,          'ATO2',           __RC__ )
    CALL MAPL_GetPointer( Internal, B3O2,          'B3O2',           __RC__ )
    CALL MAPL_GetPointer( Internal, C2H6,          'C2H6',           __RC__ )
    CALL MAPL_GetPointer( Internal, C3H8,          'C3H8',           __RC__ )
    CALL MAPL_GetPointer( Internal, CH2O,          'CH2O',           __RC__ )
    CALL MAPL_GetPointer( Internal, CH4,           'CH4',            __RC__ )
    CALL MAPL_GetPointer( Internal, CO,            'CO',             __RC__ )
    CALL MAPL_GetPointer( Internal, DRYCH2O,       'DRYCH2O',        __RC__ )
    CALL MAPL_GetPointer( Internal, DRYH2O2,       'DRYH2O2',        __RC__ )
    CALL MAPL_GetPointer( Internal, DRYHNO3,       'DRYHNO3',        __RC__ )
    CALL MAPL_GetPointer( Internal, DRYN2O5,       'DRYN2O5',        __RC__ )
    CALL MAPL_GetPointer( Internal, DRYNO2,        'DRYNO2',         __RC__ )
    CALL MAPL_GetPointer( Internal, DRYO3,         'DRYO3',          __RC__ )
    CALL MAPL_GetPointer( Internal, DRYPAN,        'DRYPAN',         __RC__ )
    CALL MAPL_GetPointer( Internal, DRYPMN,        'DRYPMN',         __RC__ )
    CALL MAPL_GetPointer( Internal, DRYPPN,        'DRYPPN',         __RC__ )
    CALL MAPL_GetPointer( Internal, DRYR4N2,       'DRYR4N2',        __RC__ )
    CALL MAPL_GetPointer( Internal, EMISSION,      'EMISSION',       __RC__ )
    CALL MAPL_GetPointer( Internal, EOH,           'EOH',            __RC__ )
    CALL MAPL_GetPointer( Internal, ETO2,          'ETO2',           __RC__ )
    CALL MAPL_GetPointer( Internal, ETP,           'ETP',            __RC__ )
    CALL MAPL_GetPointer( Internal, GCO3,          'GCO3',           __RC__ )
    CALL MAPL_GetPointer( Internal, GLCO3,         'GLCO3',          __RC__ )
    CALL MAPL_GetPointer( Internal, GLP,           'GLP',            __RC__ )
    CALL MAPL_GetPointer( Internal, GLPAN,         'GLPAN',          __RC__ )
    CALL MAPL_GetPointer( Internal, GLYC,          'GLYC',           __RC__ )
    CALL MAPL_GetPointer( Internal, GLYX,          'GLYX',           __RC__ )
    CALL MAPL_GetPointer( Internal, GP,            'GP',             __RC__ )
    CALL MAPL_GetPointer( Internal, GPAN,          'GPAN',           __RC__ )
    CALL MAPL_GetPointer( Internal, H,             'H',              __RC__ )
    CALL MAPL_GetPointer( Internal, H2,            'H2',             __RC__ )
    CALL MAPL_GetPointer( Internal, H2O,           'H2O',            __RC__ )
    CALL MAPL_GetPointer( Internal, H2O2,          'H2O2',           __RC__ )
    CALL MAPL_GetPointer( Internal, HAC,           'HAC',            __RC__ )
    CALL MAPL_GetPointer( Internal, HCOOH,         'HCOOH',          __RC__ )
    CALL MAPL_GetPointer( Internal, HNO2,          'HNO2',           __RC__ )
    CALL MAPL_GetPointer( Internal, HNO3,          'HNO3',           __RC__ )
    CALL MAPL_GetPointer( Internal, HNO4,          'HNO4',           __RC__ )
    CALL MAPL_GetPointer( Internal, HO2,           'HO2',            __RC__ )
    CALL MAPL_GetPointer( Internal, IALD,          'IALD',           __RC__ )
    CALL MAPL_GetPointer( Internal, IAO2,          'IAO2',           __RC__ )
    CALL MAPL_GetPointer( Internal, IAP,           'IAP',            __RC__ )
    CALL MAPL_GetPointer( Internal, INO2,          'INO2',           __RC__ )
    CALL MAPL_GetPointer( Internal, INPN,          'INPN',           __RC__ )
    CALL MAPL_GetPointer( Internal, ISN1,          'ISN1',           __RC__ )
    CALL MAPL_GetPointer( Internal, ISNO3,         'ISNO3',          __RC__ )
    CALL MAPL_GetPointer( Internal, ISNP,          'ISNP',           __RC__ )
    CALL MAPL_GetPointer( Internal, ISOP,          'ISOP',           __RC__ )
    CALL MAPL_GetPointer( Internal, KO2,           'KO2',            __RC__ )
    CALL MAPL_GetPointer( Internal, M,             'M',              __RC__ )
    CALL MAPL_GetPointer( Internal, MACR,          'MACR',           __RC__ )
    CALL MAPL_GetPointer( Internal, MAN2,          'MAN2',           __RC__ )
    CALL MAPL_GetPointer( Internal, MAO3,          'MAO3',           __RC__ )
    CALL MAPL_GetPointer( Internal, MAOP,          'MAOP',           __RC__ )
    CALL MAPL_GetPointer( Internal, MAP,           'MAP',            __RC__ )
    CALL MAPL_GetPointer( Internal, MCO3,          'MCO3',           __RC__ )
    CALL MAPL_GetPointer( Internal, MEK,           'MEK',            __RC__ )
    CALL MAPL_GetPointer( Internal, MGLY,          'MGLY',           __RC__ )
    CALL MAPL_GetPointer( Internal, MNO3,          'MNO3',           __RC__ )
    CALL MAPL_GetPointer( Internal, MO2,           'MO2',            __RC__ )
    CALL MAPL_GetPointer( Internal, MOH,           'MOH',            __RC__ )
    CALL MAPL_GetPointer( Internal, MP,            'MP',             __RC__ )
    CALL MAPL_GetPointer( Internal, MRO2,          'MRO2',           __RC__ )
    CALL MAPL_GetPointer( Internal, MRP,           'MRP',            __RC__ )
    CALL MAPL_GetPointer( Internal, MVK,           'MVK',            __RC__ )
    CALL MAPL_GetPointer( Internal, MVN2,          'MVN2',           __RC__ )
    CALL MAPL_GetPointer( Internal, N2,            'N2',             __RC__ )
    CALL MAPL_GetPointer( Internal, N2O5,          'N2O5',           __RC__ )
    CALL MAPL_GetPointer( Internal, NH2,           'NH2',            __RC__ )
    CALL MAPL_GetPointer( Internal, NH3,           'NH3',            __RC__ )
    CALL MAPL_GetPointer( Internal, NO,            'NO',             __RC__ )
    CALL MAPL_GetPointer( Internal, NO2,           'NO2',            __RC__ )
    CALL MAPL_GetPointer( Internal, NO3,           'NO3',            __RC__ )
    CALL MAPL_GetPointer( Internal, O,             'O',              __RC__ )
    CALL MAPL_GetPointer( Internal, O2,            'O2',             __RC__ )
    CALL MAPL_GetPointer( Internal, O2CH2OH,       'O2CH2OH',        __RC__ )
    CALL MAPL_GetPointer( Internal, O3,            'O3',             __RC__ )
    CALL MAPL_GetPointer( Internal, OH,            'OH',             __RC__ )
    CALL MAPL_GetPointer( Internal, PAN,           'PAN',            __RC__ )
    CALL MAPL_GetPointer( Internal, PMN,           'PMN',            __RC__ )
    CALL MAPL_GetPointer( Internal, PO2,           'PO2',            __RC__ )
    CALL MAPL_GetPointer( Internal, PP,            'PP',             __RC__ )
    CALL MAPL_GetPointer( Internal, PPN,           'PPN',            __RC__ )
    CALL MAPL_GetPointer( Internal, PRN1,          'PRN1',           __RC__ )
    CALL MAPL_GetPointer( Internal, PRPE,          'PRPE',           __RC__ )
    CALL MAPL_GetPointer( Internal, PRPN,          'PRPN',           __RC__ )
    CALL MAPL_GetPointer( Internal, R4N1,          'R4N1',           __RC__ )
    CALL MAPL_GetPointer( Internal, R4N2,          'R4N2',           __RC__ )
    CALL MAPL_GetPointer( Internal, R4O2,          'R4O2',           __RC__ )
    CALL MAPL_GetPointer( Internal, R4P,           'R4P',            __RC__ )
    CALL MAPL_GetPointer( Internal, RA3P,          'RA3P',           __RC__ )
    CALL MAPL_GetPointer( Internal, RB3P,          'RB3P',           __RC__ )
    CALL MAPL_GetPointer( Internal, RCHO,          'RCHO',           __RC__ )
    CALL MAPL_GetPointer( Internal, RCO3,          'RCO3',           __RC__ )
    CALL MAPL_GetPointer( Internal, RCOOH,         'RCOOH',          __RC__ )
    CALL MAPL_GetPointer( Internal, RIO1,          'RIO1',           __RC__ )
    CALL MAPL_GetPointer( Internal, RIO2,          'RIO2',           __RC__ )
    CALL MAPL_GetPointer( Internal, RIP,           'RIP',            __RC__ )
    CALL MAPL_GetPointer( Internal, ROH,           'ROH',            __RC__ )
    CALL MAPL_GetPointer( Internal, RP,            'RP',             __RC__ )
    CALL MAPL_GetPointer( Internal, VRO2,          'VRO2',           __RC__ )
    CALL MAPL_GetPointer( Internal, VRP,           'VRP',            __RC__ )
    CALL MAPL_GetPointer( Internal, DMS,           'DMS',            __RC__ )
    CALL MAPL_GetPointer( Internal, SO2,           'SO2',            __RC__ )
    CALL MAPL_GetPointer( Internal, SO4,           'SO4',            __RC__ )
    CALL MAPL_GetPointer( Internal, MSA,           'MSA',            __RC__ )
    CALL MAPL_GetPointer( Internal, LISOPOH,       'LISOPOH',        __RC__ )
    CALL MAPL_GetPointer( Internal, ALK_EMIS_SALA, 'ALK_EMIS_SALA',  __RC__ )
    CALL MAPL_GetPointer( Internal, ALK_EMIS_SALC, 'ALK_EMIS_SALC',  __RC__ )
    CALL MAPL_GetPointer( Internal, N_DENS_SALA,   'N_DENS_SALA',    __RC__ )
    CALL MAPL_GetPointer( Internal, N_DENS_SALC,   'N_DENS_SALC',    __RC__ )
    CALL MAPL_GetPointer( Internal, GSOAP11,       'GSOAP11',        __RC__ )
    CALL MAPL_GetPointer( Internal, GSOAP21,       'GSOAP21',        __RC__ )
    CALL MAPL_GetPointer( Internal, GSOAP31,       'GSOAP31',        __RC__ )
    CALL MAPL_GetPointer( Internal, GSOAP12,       'GSOAP12',        __RC__ )
    CALL MAPL_GetPointer( Internal, GSOAP22,       'GSOAP22',        __RC__ )
    CALL MAPL_GetPointer( Internal, GSOAP32,       'GSOAP32',        __RC__ )
    CALL MAPL_GetPointer( Internal, GSOAP13,       'GSOAP13',        __RC__ )
    CALL MAPL_GetPointer( Internal, GSOAP23,       'GSOAP23',        __RC__ )
    CALL MAPL_GetPointer( Internal, GSOAP33,       'GSOAP33',        __RC__ )
    CALL MAPL_GetPointer( Internal, GSOAP14,       'GSOAP14',        __RC__ )
    CALL MAPL_GetPointer( Internal, GSOAP24,       'GSOAP24',        __RC__ )
    CALL MAPL_GetPointer( Internal, GSOAP34,       'GSOAP34',        __RC__ )
    CALL MAPL_GetPointer( Internal, GSOAP15,       'GSOAP15',        __RC__ )
    CALL MAPL_GetPointer( Internal, GSOAP25,       'GSOAP25',        __RC__ )
    CALL MAPL_GetPointer( Internal, GSOAP35,       'GSOAP35',        __RC__ )
    CALL MAPL_GetPointer( Internal, GSOAP16,       'GSOAP16',        __RC__ )
    CALL MAPL_GetPointer( Internal, GSOAP26,       'GSOAP26',        __RC__ )
    CALL MAPL_GetPointer( Internal, GSOAP36,       'GSOAP36',        __RC__ )
    CALL MAPL_GetPointer( Internal, ASOAP11,       'ASOAP11',        __RC__ )
    CALL MAPL_GetPointer( Internal, ASOAP21,       'ASOAP21',        __RC__ )
    CALL MAPL_GetPointer( Internal, ASOAP31,       'ASOAP31',        __RC__ )
    CALL MAPL_GetPointer( Internal, ASOAP12,       'ASOAP12',        __RC__ )
    CALL MAPL_GetPointer( Internal, ASOAP22,       'ASOAP22',        __RC__ )
    CALL MAPL_GetPointer( Internal, ASOAP32,       'ASOAP32',        __RC__ )
    CALL MAPL_GetPointer( Internal, ASOAP13,       'ASOAP13',        __RC__ )
    CALL MAPL_GetPointer( Internal, ASOAP23,       'ASOAP23',        __RC__ )
    CALL MAPL_GetPointer( Internal, ASOAP33,       'ASOAP33',        __RC__ )
    CALL MAPL_GetPointer( Internal, ASOAP14,       'ASOAP14',        __RC__ )
    CALL MAPL_GetPointer( Internal, ASOAP24,       'ASOAP24',        __RC__ )
    CALL MAPL_GetPointer( Internal, ASOAP34,       'ASOAP34',        __RC__ )
    CALL MAPL_GetPointer( Internal, ASOAP15,       'ASOAP15',        __RC__ )
    CALL MAPL_GetPointer( Internal, ASOAP25,       'ASOAP25',        __RC__ )
    CALL MAPL_GetPointer( Internal, ASOAP35,       'ASOAP35',        __RC__ )
    CALL MAPL_GetPointer( Internal, ASOAP16,       'ASOAP16',        __RC__ )
    CALL MAPL_GetPointer( Internal, ASOAP26,       'ASOAP26',        __RC__ )
    CALL MAPL_GetPointer( Internal, ASOAP36,       'ASOAP36',        __RC__ )
    CALL MAPL_GetPointer( Internal, D_AIR_MASS,    'D_AIR_MASS',     __RC__ )
    CALL MAPL_GetPointer( Internal, D_OH_MASS,     'D_OH_MASS',      __RC__ )

    !========================================================================
    ! Check boundaries of data block on this PET
    !========================================================================

    ! Check the internal state by seeing if a pointer to it is valid
    ASSERT_( ASSOCIATED( TRC_CO ) )  
    ASSERT_( ASSOCIATED( CO     ) )  

    ! Get min & max lon, lat, lev indices
    I1     = LBOUND( CO, 1 )
    J1     = LBOUND( CO, 2 )
    K1     = LBOUND( CO, 3 )
    I2     = UBOUND( CO, 1 )
    J2     = UBOUND( CO, 2 )
    K2     = UBOUND( CO, 3 )

    ! Make sure the horizontal dimensions of the array are OK
    ASSERT_( I1 > 0 )       
    ASSERT_( I2 > 0 )
    ASSERT_( J1 > 0 )
    ASSERT_( J2 > 0 )

    ! Make sure we have 72 levels
    ASSERT_( ( K2 - K1 + 1 ) == 72 )

    !========================================================================
    ! Read data from disk and store in placeholders
    !========================================================================
    PRINT*, '### READING: ', TRIM( inFile ), 'on PET: ', myPet

    ! Open file
    OPEN( LUN, FILE=TRIM( inFile ), STATUS='OLD' )
 10 FORMAT( 4( es19.12, 1x ) )

    DO 
       ! Read the name from the file
       READ( LUN, '(a)', IOSTAT=IOS ) NAME
       IF ( IOS < 0 ) EXIT

       ! Read fields from file
       SELECT CASE( TRIM( NAME ) ) 

          ! PEDGE, PMID are placeholders
          CASE( 'MET_PEDGE' )   
             READ( LUN, 10 ) PEDGE
          CASE( 'MET_PMID' )   
             READ( LUN, 10 ) PMID

          ! Advected trcers
          CASE( 'TRC_NOx', 'TRC_NOX' )   
             READ( LUN, 10 ) GC_TRC_NOx
          CASE( 'TRC_Ox', 'TRC_OX' )    
             READ( LUN, 10 ) GC_TRC_Ox
          CASE( 'TRC_PAN' )  
             READ( LUN, 10 ) GC_TRC_PAN
          CASE( 'TRC_CO' )  
             READ( LUN, 10 ) GC_TRC_CO
          CASE( 'TRC_ALK4' )     
             READ( LUN, 10 ) GC_TRC_ALK4
          CASE( 'TRC_ISOP' )     
             READ( LUN, 10 ) GC_TRC_ISOP
          CASE( 'TRC_HNO3' )     
             READ( LUN, 10 ) GC_TRC_HNO3
          CASE( 'TRC_H2O2' )     
             READ( LUN, 10 ) GC_TRC_H2O2
          CASE( 'TRC_ACET' )     
             READ( LUN, 10 ) GC_TRC_ACET
          CASE( 'TRC_MEK' )     
             READ( LUN, 10 ) GC_TRC_MEK 
          CASE( 'TRC_ALD2' )     
             READ( LUN, 10 ) GC_TRC_ALD2
          CASE( 'TRC_RCHO' )     
             READ( LUN, 10 ) GC_TRC_RCHO
          CASE( 'TRC_MVK' )      
             READ( LUN, 10 ) GC_TRC_MVK 
          CASE( 'TRC_MACR' )     
             READ( LUN, 10 ) GC_TRC_MACR
          CASE( 'TRC_PMN' )      
             READ( LUN, 10 ) GC_TRC_PMN 
          CASE( 'TRC_PPN' )      
             READ( LUN, 10 ) GC_TRC_PPN
          CASE( 'TRC_R4N2' )     
             READ( LUN, 10 ) GC_TRC_R4N2
          CASE( 'TRC_PRPE' )     
             READ( LUN, 10 ) GC_TRC_PRPE
          CASE( 'TRC_C3H8' )     
             READ( LUN, 10 ) GC_TRC_C3H8
          CASE( 'TRC_CH2O' )     
             READ( LUN, 10 ) GC_TRC_CH2O
          CASE( 'TRC_C2H6' )     
             READ( LUN, 10 ) GC_TRC_C2H6
          CASE( 'TRC_N2O5' )     
             READ( LUN, 10 ) GC_TRC_N2O5
          CASE( 'TRC_HNO4' )     
             READ( LUN, 10 ) GC_TRC_HNO4
          CASE( 'TRC_MP' )       
             READ( LUN, 10 ) GC_TRC_MP  
          CASE( 'TRC_DMS' )      
             READ( LUN, 10 ) GC_TRC_DMS
          CASE( 'TRC_SO2' )      
             READ( LUN, 10 ) GC_TRC_SO2
          CASE( 'TRC_SO4' )      
             READ( LUN, 10 ) GC_TRC_SO4
          CASE( 'TRC_SO4s', 'TRC_SO4S' )     
             READ( LUN, 10 ) GC_TRC_SO4S
          CASE( 'TRC_MSA' )     
             READ( LUN, 10 ) GC_TRC_MSA
          CASE( 'TRC_NH3' )     
             READ( LUN, 10 ) GC_TRC_NH3
          CASE( 'TRC_NH4' )      
             READ( LUN, 10 ) GC_TRC_NH4
          CASE( 'TRC_NIT' )      
             READ( LUN, 10 ) GC_TRC_NIT
          CASE( 'TRC_NITs', 'TRC_NITS'  )     
             READ( LUN, 10 ) GC_TRC_NITS
          CASE( 'TRC_BCPI' )     
             READ( LUN, 10 ) GC_TRC_BCPI
          CASE( 'TRC_OCPI' )     
             READ( LUN, 10 ) GC_TRC_OCPI
          CASE( 'TRC_BCPO' )     
             READ( LUN, 10 ) GC_TRC_BCPO
          CASE( 'TRC_OCPO' )     
             READ( LUN, 10 ) GC_TRC_OCPO
          CASE( 'TRC_ALPH' )     
             READ( LUN, 10 ) GC_TRC_ALPH
          CASE( 'TRC_LIMO' )     
             READ( LUN, 10 ) GC_TRC_LIMO
          CASE( 'TRC_ALCO' )     
             READ( LUN, 10 ) GC_TRC_ALCO
          CASE( 'TRC_SOG1' )     
             READ( LUN, 10 ) GC_TRC_SOG1
          CASE( 'TRC_SOG2' )     
             READ( LUN, 10 ) GC_TRC_SOG2
          CASE( 'TRC_SOG3' )     
             READ( LUN, 10 ) GC_TRC_SOG3
          CASE( 'TRC_SOG4' )     
             READ( LUN, 10 ) GC_TRC_SOG4
          CASE( 'TRC_SOA1' )     
             READ( LUN, 10 ) GC_TRC_SOA1
          CASE( 'TRC_SOA2' )     
             READ( LUN, 10 ) GC_TRC_SOA2
          CASE( 'TRC_SOA3' )     
             READ( LUN, 10 ) GC_TRC_SOA3
          CASE( 'TRC_SOA4' )     
             READ( LUN, 10 ) GC_TRC_SOA4
          CASE( 'TRC_DST1' )     
             READ( LUN, 10 ) GC_TRC_DST1
          CASE( 'TRC_DST2' )     
             READ( LUN, 10 ) GC_TRC_DST2
          CASE( 'TRC_DST3' )     
             READ( LUN, 10 ) GC_TRC_DST3
          CASE( 'TRC_DST4' )     
             READ( LUN, 10 ) GC_TRC_DST4
          CASE( 'TRC_SALA' )     
             READ( LUN, 10 ) GC_TRC_SALA
          CASE( 'TRC_SALC' )     
             READ( LUN, 10 ) GC_TRC_SALC

          ! Chemical species
          CASE( 'A3O2' )
             READ( LUN, 10 ) GC_A3O2
          CASE( 'ACET' ) 
             READ( LUN, 10 ) GC_ACET
          CASE( 'ALD2' ) 
             READ( LUN, 10 ) GC_ALD2
          CASE( 'ALK4' ) 
             READ( LUN, 10 ) GC_ALK4
          CASE( 'ATO2' ) 
             READ( LUN, 10 ) GC_ATO2
          CASE( 'B3O2' ) 
             READ( LUN, 10 ) GC_B3O2
          CASE( 'C2H6' ) 
             READ( LUN, 10 ) GC_C2H6
          CASE( 'C3H8' )
             READ( LUN, 10 ) GC_C3H8
          CASE( 'CH2O' )
             READ( LUN, 10 ) GC_CH2O
          CASE( 'CO' ) 
             READ( LUN, 10 ) GC_CO
          CASE( 'DRYCH2O' )
             READ( LUN, 10 ) GC_DRYCH2O
          CASE( 'DRYH2O2' )
             READ( LUN, 10 ) GC_DRYH2O2
          CASE( 'DRYHNO3' )
             READ( LUN, 10 ) GC_DRYHNO3
          CASE( 'DRYN2O5' )
             READ( LUN, 10 ) GC_DRYN2O5
          CASE( 'DRYNO2' )
             READ( LUN, 10 ) GC_DRYNO2
          CASE( 'DRYO3' )
             READ( LUN, 10 ) GC_DRYO3
          CASE( 'DRYPAN' )
             READ( LUN, 10 ) GC_DRYPAN
          CASE( 'DRYPMN' )
             READ( LUN, 10 ) GC_DRYPMN
          CASE( 'DRYPPN' )
             READ( LUN, 10 ) GC_DRYPPN
          CASE( 'DRYR4N2' )
             READ( LUN, 10 ) GC_DRYR4N2
          CASE( 'ETO2' )
             READ( LUN, 10 ) GC_ETO2
          CASE( 'ETP' )
             READ( LUN, 10 ) GC_ETP
          CASE( 'GCO3' )
             READ( LUN, 10 ) GC_GCO3
          CASE( 'GLYC' )
             READ( LUN, 10 ) GC_GLYC
          CASE( 'GP' )
             READ( LUN, 10 ) GC_GP
          CASE( 'GPAN' )
             READ( LUN, 10 ) GC_GPAN
          CASE( 'H2O2' )
             READ( LUN, 10 ) GC_H2O2
          CASE( 'HAC' )
             READ( LUN, 10 ) GC_HAC
          CASE( 'HNO2' )
             READ( LUN, 10 ) GC_HNO2
          CASE( 'HNO3' )
             READ( LUN, 10 ) GC_HNO3
          CASE( 'HNO4' )
             READ( LUN, 10 ) GC_HNO4
          CASE( 'HO2' )
             READ( LUN, 10 ) GC_HO2
          CASE( 'IALD' )
             READ( LUN, 10 ) GC_IALD
          CASE( 'IAO2' )
             READ( LUN, 10 ) GC_IAO2
          CASE( 'IAP' )
             READ( LUN, 10 ) GC_IAP
          CASE( 'INO2' )
             READ( LUN, 10 ) GC_INO2
          CASE( 'INPN' )
             READ( LUN, 10 ) GC_INPN
          CASE( 'ISN1' )
             READ( LUN, 10 ) GC_ISN1
          CASE( 'ISNP' )
             READ( LUN, 10 ) GC_ISNP
          CASE( 'ISOP' )
             READ( LUN, 10 ) GC_ISOP
          CASE( 'KO2' )
             READ( LUN, 10 ) GC_KO2
          CASE( 'MACR' )
             READ( LUN, 10 ) GC_MACR
          CASE( 'MAN2' )
             READ( LUN, 10 ) GC_MAN2
          CASE( 'MAO3' )
             READ( LUN, 10 ) GC_MAO3
          CASE( 'MAOP' )
             READ( LUN, 10 ) GC_MAOP
          CASE( 'MAP' )
             READ( LUN, 10 ) GC_MAP
          CASE( 'MCO3' )
             READ( LUN, 10 ) GC_MCO3
          CASE( 'MEK' )
             READ( LUN, 10 ) GC_MEK
          CASE( 'MGLY' )
             READ( LUN, 10 ) GC_MGLY
          CASE( 'MO2' )
             READ( LUN, 10 ) GC_MO2
          CASE( 'MP' )
             READ( LUN, 10 ) GC_MP
          CASE( 'MRO2' )
             READ( LUN, 10 ) GC_MRO2
          CASE( 'MRP' )
             READ( LUN, 10 ) GC_MRP
          CASE( 'MVK' )
             READ( LUN, 10 ) GC_MVK
          CASE( 'MVN2' )
             READ( LUN, 10 ) GC_MVN2
          CASE( 'N2O5' )
             READ( LUN, 10 ) GC_N2O5
          CASE( 'NO' )
             READ( LUN, 10 ) GC_NO
          CASE( 'NO2' )
             READ( LUN, 10 ) GC_NO2
          CASE( 'NO3' )
             READ( LUN, 10 ) GC_NO3
          CASE( 'O3' )
             READ( LUN, 10 ) GC_O3
          CASE( 'OH' )
             READ( LUN, 10 ) GC_OH
          CASE( 'PAN' )
             READ( LUN, 10 ) GC_PAN
          CASE( 'PMN' )
             READ( LUN, 10 ) GC_PMN
          CASE( 'PO2' )
             READ( LUN, 10 ) GC_PO2
          CASE( 'PP' )
             READ( LUN, 10 ) GC_PP
          CASE( 'PPN' )
             READ( LUN, 10 ) GC_PPN
          CASE( 'PRN1' )
             READ( LUN, 10 ) GC_PRN1
          CASE( 'PRPE' )
             READ( LUN, 10 ) GC_PRPE
          CASE( 'PRPN' )
             READ( LUN, 10 ) GC_PRPN
          CASE( 'R4N1' )
             READ( LUN, 10 ) GC_R4N1
          CASE( 'R4N2' )
             READ( LUN, 10 ) GC_R4N2
          CASE( 'R4O2' )
             READ( LUN, 10 ) GC_R4O2
          CASE( 'R4P' )
             READ( LUN, 10 ) GC_R4P
          CASE( 'RA3P' )
             READ( LUN, 10 ) GC_RA3P
          CASE( 'RB3P' )
             READ( LUN, 10 ) GC_RB3P
          CASE( 'RCHO' )
             READ( LUN, 10 ) GC_RCHO
          CASE( 'RCO3' )
             READ( LUN, 10 ) GC_RCO3
          CASE( 'RIO1' )
             READ( LUN, 10 ) GC_RIO1
          CASE( 'RIO2' )
             READ( LUN, 10 ) GC_RIO2
          CASE( 'RIP' )
             READ( LUN, 10 ) GC_RIP
          CASE( 'RP' )
             READ( LUN, 10 ) GC_RP
          CASE( 'VRO2' )
             READ( LUN, 10 ) GC_VRO2
          CASE( 'VRP' )
             READ( LUN, 10 ) GC_VRP
          CASE( 'DMS' )
             READ( LUN, 10 ) GC_DMS
          CASE( 'SO2' )
             READ( LUN, 10 ) GC_SO2
          CASE( 'SO4' )
             READ( LUN, 10 ) GC_SO4
          CASE( 'MSA' )
             READ( LUN, 10 ) GC_MSA
          CASE( 'LISOPOH' )
             READ( LUN, 10 ) GC_LISOPOH
          CASE( 'ROH' )
             READ( LUN, 10 ) GC_ROH
          CASE( 'RCOOH' )
             READ( LUN, 10 ) GC_RCOOH
          CASE( 'O2CH2OH' )   
             READ( LUN, 10 ) GC_O2CH2OH
          CASE( 'O2' )   
             READ( LUN, 10 ) GC_O2
          CASE( 'O' )   
             READ( LUN, 10 ) GC_O
          CASE( 'NH3' )   
             READ( LUN, 10 ) GC_NH3
          CASE( 'NH2' )
             READ( LUN, 10 ) GC_NH2
          CASE( 'N2' )   
             READ( LUN, 10 ) GC_N2
          CASE( 'MOH' )   
             READ( LUN, 10 ) GC_MOH
          CASE( 'MNO3' )   
             READ( LUN, 10 ) GC_MNO3
          CASE( 'M' )   
             READ( LUN, 10 ) GC_M
          CASE( 'ISNO3' )   
             READ( LUN, 10 ) GC_ISNO3
          CASE( 'HCOOH' )   
             READ( LUN, 10 ) GC_HCOOH
          CASE( 'H2O' )   
             READ( LUN, 10 ) GC_H2O
          CASE( 'H2' )   
             READ( LUN, 10 ) GC_H2
          CASE( 'H' )   
             READ( LUN, 10 ) GC_H
          CASE( 'GLYX' )   
             READ( LUN, 10 ) GC_GLYX
          CASE( 'GLPAN' )   
             READ( LUN, 10 ) GC_GLPAN
          CASE( 'GLP' )   
             READ( LUN, 10 ) GC_GLP
          CASE( 'GLCO3' )   
             READ( LUN, 10 ) GC_GLCO3
          CASE( 'EOH' )   
             READ( LUN, 10 ) GC_EOH
          CASE( 'EMISSION' )   
             READ( LUN, 10 ) GC_EMISSION
          CASE( 'CH4' )   
             READ( LUN, 10 ) GC_CH4
          CASE( 'ACTA' )   
             READ( LUN, 10 ) GC_ACTA

          ! Other quantities
          CASE( 'H2O2S','H2O2s' ) 
             READ( LUN, 10 ) GC_H2O2s
          CASE( 'SO2S', 'SO2s' ) 
             READ( LUN, 10 ) GC_SO2s
          CASE( 'ORVC_TERP' ) 
             READ( LUN, 10 ) GC_ORVC_TERP
          CASE( 'ORVC_SESQ' )
             READ( LUN, 10 ) GC_ORVC_SESQ
          CASE DEFAULT
             ! Nothing
       END SELECT

    ENDDO

    ! Close file
    CLOSE( LUN )

    !========================================================================
    ! Fill internal state with values saved out from GEOS-Chem
    ! Assign values from the same column (read from disk) everywhere
    !========================================================================
    DO J = J1, J2
    DO I = I1, I2

       !--------------------------------------------------------------------
       ! Advected tracers
       !--------------------------------------------------------------------
       TRC_NOx      (I,J,:) = GC_TRC_NOx  (:)
       TRC_Ox       (I,J,:) = GC_TRC_Ox   (:)
       TRC_PAN      (I,J,:) = GC_TRC_PAN  (:)
       TRC_CO       (I,J,:) = GC_TRC_CO   (:)
       TRC_ALK4     (I,J,:) = GC_TRC_ALK4 (:)
       TRC_ISOP     (I,J,:) = GC_TRC_ISOP (:)
       TRC_HNO3     (I,J,:) = GC_TRC_HNO3 (:)
       TRC_H2O2     (I,J,:) = GC_TRC_H2O2 (:)
       TRC_ACET     (I,J,:) = GC_TRC_ACET (:)
       TRC_MEK      (I,J,:) = GC_TRC_MEK  (:)
       TRC_ALD2     (I,J,:) = GC_TRC_ALD2 (:)
       TRC_RCHO     (I,J,:) = GC_TRC_RCHO (:)
       TRC_MVK      (I,J,:) = GC_TRC_MVK  (:)
       TRC_MACR     (I,J,:) = GC_TRC_MACR (:)
       TRC_PMN      (I,J,:) = GC_TRC_PMN  (:)
       TRC_PPN      (I,J,:) = GC_TRC_PPN  (:)
       TRC_R4N2     (I,J,:) = GC_TRC_R4N2 (:)
       TRC_PRPE     (I,J,:) = GC_TRC_PRPE (:)
       TRC_C3H8     (I,J,:) = GC_TRC_C3H8 (:)
       TRC_CH2O     (I,J,:) = GC_TRC_CH2O (:)
       TRC_C2H6     (I,J,:) = GC_TRC_C2H6 (:)
       TRC_N2O5     (I,J,:) = GC_TRC_N2O5 (:)
       TRC_HNO4     (I,J,:) = GC_TRC_HNO4 (:)
       TRC_MP       (I,J,:) = GC_TRC_MP   (:)
       TRC_DMS      (I,J,:) = GC_TRC_DMS  (:)
       TRC_SO2      (I,J,:) = GC_TRC_SO2  (:)
       TRC_SO4      (I,J,:) = GC_TRC_SO4  (:)
       TRC_SO4s     (I,J,:) = GC_TRC_SO4s (:)
       TRC_MSA      (I,J,:) = GC_TRC_MSA  (:)
       TRC_NH3      (I,J,:) = GC_TRC_NH3  (:)
       TRC_NH4      (I,J,:) = GC_TRC_NH4  (:)
       TRC_NIT      (I,J,:) = GC_TRC_NIT  (:)
       TRC_NITs     (I,J,:) = GC_TRC_NITs (:)
       TRC_BCPI     (I,J,:) = GC_TRC_BCPI (:)
       TRC_OCPI     (I,J,:) = GC_TRC_OCPI (:)
       TRC_BCPO     (I,J,:) = GC_TRC_BCPO (:)
       TRC_OCPO     (I,J,:) = GC_TRC_OCPO (:)
       TRC_ALPH     (I,J,:) = GC_TRC_ALPH (:)
       TRC_LIMO     (I,J,:) = GC_TRC_LIMO (:)
       TRC_ALCO     (I,J,:) = GC_TRC_ALCO (:)
       TRC_SOG1     (I,J,:) = GC_TRC_SOG1 (:)
       TRC_SOG2     (I,J,:) = GC_TRC_SOG2 (:)
       TRC_SOG3     (I,J,:) = GC_TRC_SOG3 (:)
       TRC_SOG4     (I,J,:) = GC_TRC_SOG4 (:)
       TRC_SOA1     (I,J,:) = GC_TRC_SOA1 (:)
       TRC_SOA2     (I,J,:) = GC_TRC_SOA2 (:)
       TRC_SOA3     (I,J,:) = GC_TRC_SOA3 (:)
       TRC_SOA4     (I,J,:) = GC_TRC_SOA4 (:)
       TRC_DST1     (I,J,:) = GC_TRC_DST1 (:)
       TRC_DST2     (I,J,:) = GC_TRC_DST2 (:)
       TRC_DST3     (I,J,:) = GC_TRC_DST3 (:)
       TRC_DST4     (I,J,:) = GC_TRC_DST4 (:)
       TRC_SALA     (I,J,:) = GC_TRC_SALA (:)
       TRC_SALC     (I,J,:) = GC_TRC_SALC (:)

       !--------------------------------------------------------------------
       ! Chemical species (both active & inactive)
       !--------------------------------------------------------------------
       A3O2         (I,J,:) = GC_A3O2     (:)         
       ACET         (I,J,:) = GC_ACET     (:) 
       ALD2         (I,J,:) = GC_ALD2     (:) 
       ALK4         (I,J,:) = GC_ALK4     (:) 
       ATO2         (I,J,:) = GC_ATO2     (:) 
       B3O2         (I,J,:) = GC_B3O2     (:) 
       C2H6         (I,J,:) = GC_C2H6     (:) 
       C3H8         (I,J,:) = GC_C3H8     (:)
       CH2O         (I,J,:) = GC_CH2O     (:)
       CO           (I,J,:) = GC_CO       (:) 
       DRYCH2O      (I,J,:) = GC_DRYCH2O  (:)
       DRYH2O2      (I,J,:) = GC_DRYH2O2  (:)
       DRYHNO3      (I,J,:) = GC_DRYHNO3  (:)
       DRYN2O5      (I,J,:) = GC_DRYN2O5  (:)
       DRYNO2       (I,J,:) = GC_DRYNO2   (:)
       DRYO3        (I,J,:) = GC_DRYO3    (:)
       DRYPAN       (I,J,:) = GC_DRYPAN   (:)
       DRYPMN       (I,J,:) = GC_DRYPMN   (:)
       DRYPPN       (I,J,:) = GC_DRYPPN   (:)
       DRYR4N2      (I,J,:) = GC_DRYR4N2  (:)
       ETO2         (I,J,:) = GC_ETO2     (:)
       ETP          (I,J,:) = GC_ETP      (:)
       GCO3         (I,J,:) = GC_GCO3     (:)
       GLYC         (I,J,:) = GC_GLYC     (:)
       GP           (I,J,:) = GC_GP       (:)
       GPAN         (I,J,:) = GC_GPAN     (:)
       H2O2         (I,J,:) = GC_H2O2     (:)
       HAC          (I,J,:) = GC_HAC      (:)
       HNO2         (I,J,:) = GC_HNO2     (:)
       HNO3         (I,J,:) = GC_HNO3     (:)
       HNO4         (I,J,:) = GC_HNO4     (:)
       HO2          (I,J,:) = GC_HO2      (:)
       IALD         (I,J,:) = GC_IALD     (:)
       IAO2         (I,J,:) = GC_IAO2     (:)
       IAP          (I,J,:) = GC_IAP      (:)
       INO2         (I,J,:) = GC_INO2     (:)
       INPN         (I,J,:) = GC_INPN     (:)
       ISN1         (I,J,:) = GC_ISN1     (:)
       ISNP         (I,J,:) = GC_ISNP     (:)
       ISOP         (I,J,:) = GC_ISOP     (:)
       KO2          (I,J,:) = GC_KO2      (:)
       MACR         (I,J,:) = GC_MACR     (:)
       MAN2         (I,J,:) = GC_MAN2     (:)
       MAO3         (I,J,:) = GC_MAO3     (:)
       MAOP         (I,J,:) = GC_MAOP     (:)
       MAP          (I,J,:) = GC_MAP      (:)
       MCO3         (I,J,:) = GC_MCO3     (:)
       MEK          (I,J,:) = GC_MEK      (:)
       MGLY         (I,J,:) = GC_MGLY     (:)
       MO2          (I,J,:) = GC_MO2      (:)
       MP           (I,J,:) = GC_MP       (:)
       MRO2         (I,J,:) = GC_MRO2     (:)
       MRP          (I,J,:) = GC_MRP      (:)
       MVK          (I,J,:) = GC_MVK      (:)
       MVN2         (I,J,:) = GC_MVN2     (:)
       N2O5         (I,J,:) = GC_N2O5     (:)
       NO           (I,J,:) = GC_NO       (:)
       NO2          (I,J,:) = GC_NO2      (:)
       NO3          (I,J,:) = GC_NO3      (:)
       O3           (I,J,:) = GC_O3       (:)
       OH           (I,J,:) = GC_OH       (:)
       PAN          (I,J,:) = GC_PAN      (:)
       PMN          (I,J,:) = GC_PMN      (:)
       PO2          (I,J,:) = GC_PO2      (:)
       PP           (I,J,:) = GC_PP       (:)
       PPN          (I,J,:) = GC_PPN      (:)
       PRN1         (I,J,:) = GC_PRN1     (:)
       PRPE         (I,J,:) = GC_PRPE     (:)
       PRPN         (I,J,:) = GC_PRPN     (:)
       R4N1         (I,J,:) = GC_R4N1     (:)
       R4N2         (I,J,:) = GC_R4N2     (:)
       R4O2         (I,J,:) = GC_R4O2     (:)
       R4P          (I,J,:) = GC_R4P      (:)
       RA3P         (I,J,:) = GC_RA3P     (:)
       RB3P         (I,J,:) = GC_RB3P     (:)
       RCHO         (I,J,:) = GC_RCHO     (:)
       RCO3         (I,J,:) = GC_RCO3     (:)
       RIO1         (I,J,:) = GC_RIO1     (:)
       RIO2         (I,J,:) = GC_RIO2     (:)
       RIP          (I,J,:) = GC_RIP      (:)
       RP           (I,J,:) = GC_RP       (:)
       VRO2         (I,J,:) = GC_VRO2     (:)
       VRP          (I,J,:) = GC_VRP      (:)
       DMS          (I,J,:) = GC_DMS      (:)
       SO2          (I,J,:) = GC_SO2      (:)
       SO4          (I,J,:) = GC_SO4      (:)
       MSA          (I,J,:) = GC_MSA      (:)
       LISOPOH      (I,J,:) = GC_LISOPOH  (:)
       ROH          (I,J,:) = GC_ROH      (:)
       RCOOH        (I,J,:) = GC_RCOOH    (:)
       O2CH2OH      (I,J,:) = GC_O2CH2OH  (:)   
       O2           (I,J,:) = GC_O2       (:)   
       O            (I,J,:) = GC_O        (:)   
       NH3          (I,J,:) = GC_NH3      (:)   
       NH2          (I,J,:) = GC_NH2      (:)
       N2           (I,J,:) = GC_N2       (:)   
       MOH          (I,J,:) = GC_MOH      (:)   
       MNO3         (I,J,:) = GC_MNO3     (:)   
       M            (I,J,:) = GC_M        (:)   
       ISNO3        (I,J,:) = GC_ISNO3    (:)   
       HCOOH        (I,J,:) = GC_HCOOH    (:)   
       H2O          (I,J,:) = GC_H2O      (:)   
       H2           (I,J,:) = GC_H2       (:)   
       H            (I,J,:) = GC_H        (:)   
       GLYX         (I,J,:) = GC_GLYX     (:)   
       GLPAN        (I,J,:) = GC_GLPAN    (:)   
       GLP          (I,J,:) = GC_GLP      (:)   
       GLCO3        (I,J,:) = GC_GLCO3    (:)   
       EOH          (I,J,:) = GC_EOH      (:)   
       EMISSION     (I,J,:) = GC_EMISSION (:)   
       CH4          (I,J,:) = GC_CH4      (:)   
       ACTA         (I,J,:) = GC_ACTA     (:)   

       !--------------------------------------------------------------------
       ! Other quantities
       !--------------------------------------------------------------------

       ! For wetdep
       H2O2s        (I,J,:) = GC_H2O2s    (:)   
       SO2S         (I,J,:) = GC_SO2s     (:)   

       ! For secondary organic aerosols
       ORVC_TERP    (I,J,:) = GC_ORVC_TERP(:)   
       ORVC_SESQ    (I,J,:) = GC_ORVC_SESQ(:)   

       ! All of the GPROD/APROD parameters are zero at the start of the run
       GSOAP11      (I,J,:) = 0e0
       GSOAP21      (I,J,:) = 0e0
       GSOAP31      (I,J,:) = 0e0
       GSOAP12      (I,J,:) = 0e0
       GSOAP22      (I,J,:) = 0e0
       GSOAP32      (I,J,:) = 0e0
       GSOAP13      (I,J,:) = 0e0
       GSOAP23      (I,J,:) = 0e0
       GSOAP33      (I,J,:) = 0e0
       GSOAP14      (I,J,:) = 0e0
       GSOAP24      (I,J,:) = 0e0
       GSOAP34      (I,J,:) = 0e0
       GSOAP15      (I,J,:) = 0e0
       GSOAP25      (I,J,:) = 0e0
       GSOAP35      (I,J,:) = 0e0
       GSOAP16      (I,J,:) = 0e0
       GSOAP26      (I,J,:) = 0e0
       GSOAP36      (I,J,:) = 0e0
       ASOAP11      (I,J,:) = 0e0
       ASOAP21      (I,J,:) = 0e0
       ASOAP31      (I,J,:) = 0e0
       ASOAP12      (I,J,:) = 0e0
       ASOAP22      (I,J,:) = 0e0
       ASOAP32      (I,J,:) = 0e0
       ASOAP13      (I,J,:) = 0e0
       ASOAP23      (I,J,:) = 0e0
       ASOAP33      (I,J,:) = 0e0
       ASOAP14      (I,J,:) = 0e0
       ASOAP24      (I,J,:) = 0e0
       ASOAP34      (I,J,:) = 0e0
       ASOAP15      (I,J,:) = 0e0
       ASOAP25      (I,J,:) = 0e0
       ASOAP35      (I,J,:) = 0e0
       ASOAP16      (I,J,:) = 0e0
       ASOAP26      (I,J,:) = 0e0
       ASOAP36      (I,J,:) = 0e0

       ! Sea salt parameters are zero at the initial timestep
       ALK_EMIS_SALA(I,J,:) = 0e0
       ALK_EMIS_SALC(I,J,:) = 0e0
       N_DENS_SALA  (I,J,:) = 0e0
       N_DENS_SALC  (I,J,:) = 0e0

       ! Mean OH lifetime arrays are zero at the initial timestep
       D_AIR_MASS   (I,J,:) = 0e0
       D_OH_MASS    (I,J,:) = 0e0

    ENDDO
    ENDDO

  END SUBROUTINE Fill_Internal_State_
!EOC
END MODULE gc_esmf_drv
#endif
