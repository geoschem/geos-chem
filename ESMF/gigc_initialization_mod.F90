#if defined (ESMF_)
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: gigc_initialization_mod
!
! !DESCRIPTION: Module GIGC\_INITIALIZATION\_MOD is the module that
!  the initialize methods for the ESMF interface to the Grid-Independent
!  GEOS-Chem (aka "GIGC").
!\\
!\\
! !INTERFACE: 
!      
MODULE GIGC_Initialization_Mod
!
! !USES:
!      
  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!      
  PUBLIC :: GIGC_Init_Chemistry
  PUBLIC :: GIGC_Emis_Inti
  PUBLIC :: GIGC_Init_Simulation
  PUBLIC :: GIGC_Get_Options
  PUBLIC :: GIGC_Init_Time_Interface
  PUBLIC :: GIGC_Init_Dimensions
  PUBLIC :: GIGC_Allocate_Interface
!
! !REVISION HISTORY: 
!  16 Oct 2012 - M. Long     - Initial version
!  16 Oct 2012 - R. Yantosca - Added ProTeX headers
!  22 Oct 2012 - R. Yantosca - Renamed to gigc_initialization_mod.F90
!  22 Oct 2012 - R. Yantosca - Renamed several routines for better consistency
!  25 Oct 2012 - R. Yantosca - Remove routine GIGC_SetEnv
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
! !IROUTINE: gigc_init_chemistry
!
! !DESCRIPTION: Routine GIGC\_Init\_Chemistry initializes the GEOS-Chem 
!  chemistry module so that it can connect to the ESMF interface.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE GIGC_Init_Chemistry( PLON, PLAT, am_I_Root, RC )
!
! !USES:
!
    USE LOGICAL_MOD,   ONLY : DO_DIAG_WRITE
    USE CMN_SIZE_MOD,  ONLY : JJPAR, IIPAR
!
! !INPUT PARAMETERS: 
!
    INTEGER, INTENT(IN)  :: PLON        ! Number of
    INTEGER, INTENT(IN)  :: PLAT
    LOGICAL, INTENT(IN)  :: am_I_Root   ! Is this the root CPU?
!
! !OUTPUT PARAMETERS:
!
    INTEGER, INTENT(OUT) :: RC          ! Success or failure?  

! 
! !REVISION HISTORY: 
!  15 Oct 2012 - M. Long     - Initial version
!  15 Oct 2012 - R. Yantosca - Added ProTeX Headers, use F90 format/indents
!EOP
!------------------------------------------------------------------------------
!BOC

    ! Don't write diagnostic files
    DO_DIAG_WRITE = .FALSE. 

    !### Debug
    ! SET MODEL DIMENSIONS 
    !IIPAR = 1
    !JJPAR = PCOLS

    !=======================================================================
    ! Read GEOS-Chem emissions data ets
    !=======================================================================
   
    ! For now comment this out
    !CALL GC_EMIS_INTI( PLONL, PLATL, 1 )

  END SUBROUTINE GIGC_Init_Chemistry
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: gigc_emis_inti
!
! !DESCRIPTION: Routine GC\_EMIS\_INTI establishes the emissions fields for
!  each set of geos-chem emissions used.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE GIGC_Emis_Inti( PLONL, PLATL, PPLON )
!
! !USES:
!
!    USE GC_EMIS_EDGAR, ONLY: GC_EDGAR_EMIS_INTI
!
! !INPUT PARAMETERS: 
!
    INTEGER, INTENT(IN) :: PLONL
    INTEGER, INTENT(IN) :: PLATL
    INTEGER, INTENT(IN) :: PPLON
! 
! !REMARKS:
!  This is a stub routine for now.  Will eventually connect with the Emissions
!  Component as being developed by Christoph Keller.
!
! !REVISION HISTORY: 
!  15 Oct 2012 - M. Long     - Initial version
!  15 Oct 2012 - R. Yantosca - Added ProTeX Headers, use F90 format/indents
!  22 Oct 2012 - R. Yantosca - Renamed to GIGC_Emis_Inti
!EOP
!------------------------------------------------------------------------------
!BOC

    ! Comment out for now
    !CALL GC_EDGAR_EMIS_INTI( PLONL, PLATL, PPLON )

  END SUBROUTINE GIGC_Emis_Inti
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: gigc_init_simulation
!
! !DESCRIPTION: Routine GIGC\_INIT\_SIMULATION is the Initialize method for 
!  the ESMF interface that connects the Grid-Independent GEOS-Chem (aka "GIGC")
!  to the GEOS-5 GCM.  Calls to the various GEOS-Chem init routines (which 
!  allocate arrays, etc.) are made from here.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE GIGC_Init_Simulation( State_Met, State_Chm, tsChem,     &
                                   nymd,      nhms,      am_I_Root,  &
                                   RC                               )
!
! !USES:
!
    USE GIGC_Environment_Mod
    USE GIGC_ErrCode_Mod  
    USE GIGC_State_Chm_Mod
    USE GIGC_State_Met_Mod
    USE CMN_SIZE_MOD,         ONLY : IIPAR
    USE CMN_SIZE_MOD,         ONLY : JJPAR
    USE CMN_SIZE_MOD,         ONLY : LLTROP
    USE COMODE_MOD
    USE COMODE_LOOP_MOD       
    USE GCKPP_COMODE_MOD,     ONLY : INIT_GCKPP_COMODE
    USE GRID_MOD,             ONLY : INIT_GRID
    USE DAO_MOD
    USE LOGICAL_MOD
    USE PBL_MIX_MOD,          ONLY : INIT_PBL_MIX
    USE PRESSURE_MOD,         ONLY : INIT_PRESSURE
    USE TRACER_MOD,           ONLY : ITS_A_FULLCHEM_SIM
    USE TRACER_MOD,           ONLY : ITS_AN_AEROSOL_SIM
    USE TRACER_MOD,           ONLY : ID_TRACER
    USE TRACER_MOD,           ONLY : TRACER_NAME
    USE TRACER_MOD,           ONLY : N_TRACERS
    USE TRACERID_MOD,         ONLY : SETTRACE
    USE TOMS_MOD,             ONLY : TO3_DAILY
    USE WETSCAV_MOD,          ONLY : INIT_WETSCAV
    USE ERROR_MOD,            ONLY : DEBUG_MSG

    ! Comment these out for now (bmy, 10/15/12)
    !USE TIME_MANAGER,       ONLY : GET_STEP_SIZE
    !USE GC_GRIDAREA
!
! !INPUT PARAMETERS: 
!
    REAL,           INTENT(IN)    :: tsChem      ! Chemistry timestep [s]
    INTEGER,        INTENT(IN)    :: nymd        ! GMT date (YYYY/MM/DD)
    INTEGER,        INTENT(IN)    :: nhms        ! GMT time (hh:mm:ss)
    LOGICAL,        INTENT(IN)    :: am_I_Root   ! Is this the root CPU?
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
    TYPE(MetState), INTENT(INOUT) :: State_Met   ! Meteorology State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?  
!
! !REMARKS
!  Add other calls to GEOS-Chem init routines as necessary.
!  NOTE: Later on maybe split these init calls among other routines.
!  Also need to add better error trapping
!
! !REVISION HISTORY: 
!  15 Oct 2012 - M. Long     - Initial version
!  15 Oct 2012 - R. Yantosca - Added ProTeX Headers, use F90 format/indents
!  17 Oct 2012 - R. Yantosca - Now initialize the chemistry mechanism
!  19 Oct 2012 - R. Yantosca - Now reference gigc_state_chm_mod.F90
!  19 Oct 2012 - R. Yantosca - Now reference gigc_state_met_mod.F90
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
! 
    INTEGER :: DTIME, K, AS, N, PLON, PLAT

    !=======================================================================
    ! Initialize key GEOS-Chem sections
    !=======================================================================

    ! Initialize
    RC    = GIGC_SUCCESS
    DTIME = tsChem
       
    ! Initialize the timing routines
    CALL GIGC_Init_Time_Interface( DTIME, nymd, nhms, am_I_Root, RC )

    ! Allocate all
    CALL GIGC_Allocate_All( am_I_Root, RC )

    ! Allocate
    CALL GIGC_Allocate_Interface( am_I_Root, RC )

    ! Read options from the GEOS-Chem input file "input.geos"
    CALL GIGC_Get_Options( am_I_Root, RC )

    ! Initialize derived-type objects for meteorology & chemistry states
    CALL GIGC_Init_All( State_Met, State_Chm, am_I_Root, RC )

    ! Save tracer names and ID's into State_Chm
    DO N = 1, N_TRACERS
       State_Chm%TRAC_NAME(N) = 'TRC_' // TRIM( TRACER_NAME(N) )
       State_Chm%TRAC_ID  (N) = ID_TRACER(N)
    ENDDO
       
    ! Allocate and zero GEOS-Chem diagnostic arrays
    CALL NDXX_SETUP

    ! Initialize 
    CALL GIGC_Init_Chemistry( PLON, PLAT, am_I_Root, RC )

    ! Allocate and initialize met field arrays
    CALL INIT_DAO

    ! Initialize the GEOS-Chem pressure module (set Ap & Bp)
    CALL INIT_PRESSURE( am_I_Root )

    ! Initialize the PBL mixing module
    CALL INIT_PBL_MIX
!       
!    ! Initialize the  
!    CALL INIT_WETSCAV

    !=======================================================================
    ! Initialize chemistry mechanism
    !=======================================================================

    ! Set some size variables
    NLAT   = JJPAR
    NLONG  = IIPAR
    NVERT  = IVERT 
    NPVERT = NVERT
    NPVERT = NVERT + IPLUME

!    ! INITIALIZE ALLOCATABLE SMVGEAR/KPP ARRAYS
!    IF ( LEMIS .OR. LCHEM ) THEN

       ! Initialize arrays in comode_mod.F
       CALL INIT_COMODE( am_I_Root )

       ! Initialize KPP (if necessary)
       IF ( LKPP ) THEN
          CALL INIT_GCKPP_COMODE( am_I_Root, IIPAR,   JJPAR, LLTROP,  &
                                  ITLOOP,    NMTRATE, IGAS,  RC      )
       ENDIF
!    ENDIF

    ! Read from data file mglob.dat
    CALL READER( .TRUE., am_I_Root )

    !### Debug
    IF ( LPRT .and. am_I_Root ) THEN
       CALL DEBUG_MSG( '### GIGC_INIT_CHEMISTRY: after READER' )
    ENDIF

    ! Set NCS for urban chemistry only (since that is where we
    ! have defined the GEOS-CHEM mechanism) (bdf, bmy, 4/21/03)
    NCS = NCSURBAN

    !### Debug
    IF ( LPRT .and. am_I_Root ) THEN
       CALL DEBUG_MSG( '### GIGC_INIT_CHEMISTRY: after READER' )        
    ENDIF

    ! Redefine NTLOOP since READER defines it initially (bmy, 9/28/04)
    NLOOP   = NLAT  * NLONG
    NTLOOP  = NLOOP * NVERT
    NTTLOOP = NTLOOP

    ! Read "globchem.dat" chemistry mechanism
    CALL READCHEM( am_I_Root )

    !### Debug
    IF ( LPRT .and. am_I_Root ) THEN
       CALL DEBUG_MSG( '### GIGC_INIT_CHEMISTRY: after READCHEM' )        
    ENDIF

    ! Save Chemical species names ID's into State_Chm
    DO N = 1, IGAS
       IF ( LEN_TRIM( NAMEGAS(N) ) > 0 ) THEN 
          State_Chm%SPEC_NAME(N) = TRIM( NAMEGAS(N) )
          State_Chm%SPEC_ID  (N) = N
       ENDIF
    ENDDO

    ! Set NCS=NCSURBAN here since we have defined our tropospheric
    ! chemistry mechanism in the urban slot of SMVGEAR II (bmy, 4/21/03)
    NCS = NCSURBAN

    !%%%%% FOR NOW HARDWIRE CH4 to 2007 values %%%%%
    ! Get CH4 [ppbv] in 4 latitude bins for each year
    CALL GET_GLOBAL_CH4( 2007,   .TRUE., C3090S, &
                         C0030S, C0030N, C3090N, am_I_Root )

    !### Debug
    IF ( LPRT .and. am_I_Root ) THEN
       CALL DEBUG_MSG( '### GIGC_INIT_CHEMISTRY: after GET_GLOBAL_CH4' )        
    ENDIF

    ! Initialize FAST-J photolysis
    CALL INPHOT( LLPAR, NPHOT, am_I_Root ) 
         
    !### Debug
    IF ( LPRT .and. am_I_Root ) THEN
       CALL DEBUG_MSG( '### GIGC_INIT_CHEMISTRY: after INPHOT' )        
    ENDIF

    ! Flag certain chemical species
    CALL SETTRACE( State_Chm, am_I_Root )

    !### Debug
    IF ( LPRT .and. am_I_Root ) THEN
       CALL DEBUG_MSG( '### GIGC_INIT_CHEMISTRY: after SETTRACE' )
    ENDIF

    ! Flag emission & drydep rxns
    CALL SETEMDEP( N_TRACERS, am_I_Root )

    !### Debug
    IF ( LPRT .and. am_I_Root ) THEN
       CALL DEBUG_MSG( '### GIGC_INIT_CHEMISTRY: after SETEMDEP' )
    ENDIF

    ! Initialize dry deposition (work in progress), add here
    ! Allocate array of overhead O3 columns for TOMS
    ALLOCATE( TO3_DAILY( IIPAR, JJPAR ), STAT=AS )
    TO3_DAILY = 0d0

  END SUBROUTINE GIGC_Init_Simulation
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: gigc_get_options
!
! !DESCRIPTION: Routine GIGC\_GET\_OPTIONS reads options for a GEOS-Chem 
!  simulation from the input.geos\_\_\_.rc input file.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE GIGC_Get_Options( am_I_Root, RC )
!
! !USES:
!
    USE GIGC_ErrCode_Mod
    USE Input_Mod,      ONLY : Read_Input_File
!
! !INPUT PARAMETERS: 
!
    LOGICAL, INTENT(IN)  :: am_I_Root   ! Are we on the root CPU?
!
! !OUTPUT PARAMETERS:
!
    INTEGER, INTENT(OUT) :: RC          ! Success or failure
! 
! !REMARKS:
!  NOTE: We will probably convert the input file into an ESMF resource file
!  in the near future.
!
! !REVISION HISTORY: 
!  15 Oct 2012 - M. Long     - Initial version
!  15 Oct 2012 - R. Yantosca - Added ProTeX Headers, use F90 format/indents
!  22 Oct 2012 - R. Yantosca - Renamed to GIGC_Get_Options
!  22 Oct 2012 - R. Yantosca - Added RC output argument

!EOP
!------------------------------------------------------------------------------
!BOC

    ! Assume succes
    RC = GIGC_Success
    
    ! Read the GEOS-Chem input file here
    CALL Read_Input_File( am_I_Root )

  END SUBROUTINE GIGC_Get_Options
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: gigc_init_time_interface
!
! !DESCRIPTION: Routine GIGC\_INIT\_TIME\_INTERFACE intializes the starting 
!  date of the run as well as the dynamic timestep for the Grid-Independent
!  GEOS-Chem (aka "GIGC") simulation.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE GIGC_Init_Time_Interface( DT, nymd, nhms, am_I_Root, RC )
!
! !USES:
!
    USE GIGC_ErrCode_Mod
    USE JulDay_Mod
    USE Time_Mod
!
! !INPUT PARAMETERS: 
!
    INTEGER, INTENT(IN)  :: DT          ! Timestep [seconds]
    INTEGER, INTENT(IN)  :: nymd        ! GMT date (YYYY/MM/DD)
    INTEGER, INTENT(IN)  :: nhms        ! GMT time (hh:mm:ss)
    LOGICAL, INTENT(IN)  :: am_I_Root   ! Are we on the root CPU
!
! !OUTPUT PARAMETERS:
!
    INTEGER, INTENT(OUT) :: RC          ! Success or failure
! 
! !REVISION HISTORY: 
!  15 Oct 2012 - M. Long     - Initial version
!  15 Oct 2012 - R. Yantosca - Added ProTeX Headers, use F90 format/indents
!  15 Oct 2012 - R. Yantosca - Need to pass am_I_Root to INITIALIZE from here
!  22 Oct 2012 - R. Yantosca - Renamed to GIGC_Init_Time_Interface
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: NYMDB, NHMSB, Y, M, D, TOD, H, MN, S
    INTEGER :: DTM
    REAL*8  :: FRAC_DAY
    
    ! Assume success
    RC = GIGC_SUCCESS

    ! Convert timestep from seconds to minutes
    DTM =  DT / 60

    ! Set GEOS-Chem timesteps (for now, set all to the same value)
    CALL SET_TIMESTEPS( am_I_Root, DTM, DTM, DTM, DTM, DTM, DTM )
       
    ! SET NYMDB & NHMSB
    Y    = nymd / 10000                    ! Year
    M    = ( nymd - Y*10000 ) /100         ! Month
    D    = nymd - Y*10000 - M*100          ! Day
                                           
    H    = nhms/10000                      ! Hour
    MN   = ( nhms - H*10000 ) /100         ! Minute
    S    = nhms - H*10000 - MN*100         ! Seconds
                                           
    TOD = H*3600 + MN*60 + S               ! Seconds since start of day


    ! Fraction of the day that has elapsed
    IF ( TOD .GT. 0 ) THEN
       FRAC_DAY = DBLE( D ) + DBLE( TOD ) / DBLE( 86400 )   
    ELSE
       FRAC_DAY = DBLE( D )
    ENDIF
    
    ! Find the starting date
    CALL CALDATE( JULDAY( Y, M, FRAC_DAY ), NYMDB, NHMSB )

    ! Set the starting date & time of the GEOS-Chem simulation
    CALL SET_BEGIN_TIME( NYMDB, NHMSB )

    ! Set the current time to the starting time
    CALL SET_CURRENT_TIME()

    ! INITIALIZE AND SET TIMESTEP COUNTS USING
    ! GEOS-CHEM'S INITIALIZE (SEE GEOS-CHEM'S MAIN.F)
    CALL INITIALIZE( 2, am_I_Root )
    CALL INITIALIZE( 3, am_I_Root )

  END SUBROUTINE GIGC_Init_Time_Interface
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: gigc_init_dimensions
!
! !DESCRIPTION: Routine GIGC\_INIT\_DIMENSIONS initializes the geospatial
!  dimensions for a Grid-Independent GEOS-Chem (aka "GIGC") simulation 
!  directly from the ESMF interface.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE GIGC_Init_Dimensions( NI, NJ, NL )
!
! !USES:
!
    USE CMN_SIZE_MOD, ONLY: IIPAR, JJPAR, LLPAR, IGLOB, JGLOB, LGLOB
!
! !INPUT PARAMETERS: 
!
    INTEGER, INTENT(IN) :: NI    ! Size of local I dimension
    INTEGER, INTENT(IN) :: NJ    ! Size of local J dimension
    INTEGER, INTENT(IN) :: NL    ! Number of levels
! 
! !REVISION HISTORY: 
!  15 Oct 2012 - M. Long     - Initial version
!  15 Oct 2012 - R. Yantosca - Added ProTeX Headers, use F90 format/indents
!  22 Oct 2012 - R. Yantosca - Renamed to GIGC_Init_Dimensions

!EOP
!------------------------------------------------------------------------------
!BOC
    
    ! Set GEOS-Chem size variables from the locally-defined 
    ! dimensions returned by the ESMF framework
    IGLOB = NI
    JGLOB = NJ
    LGLOB = NL
    IIPAR = NI
    JJPAR = NJ
    LLPAR = NL

  END SUBROUTINE GIGC_Init_Dimensions
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: gigc_allocate_interface
!
! !DESCRIPTION: Routine GIGC\_ALLOCATE\_INTERFACE allocates GEOS-Chem arrays
!  for a Grid-Independent GEOS-Chem (aka "GIGC") simulation.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE GIGC_Allocate_Interface( am_I_Root, RC )
!
! !USES:
!
    USE CMN_SIZE_MOD,     ONLY : IIPAR, JJPAR, LLPAR
    USE GIGC_ErrCode_Mod
    USE ERROR_MOD,        ONLY : ALLOC_ERR
    USE UVALBEDO_MOD,     ONLY : UVALBEDO
   !USE DAO_MOD,          ONLY : TO3
!
! !INPUT PARAMETERS: 
!
    LOGICAL, INTENT(IN)  :: am_I_Root   ! Are we on the root CPU?
!
! !OUTPUT PARAMETERS:
!
    INTEGER, INTENT(OUT) :: RC          ! Success or failure
! 
! !REVISION HISTORY: 
!  15 Oct 2012 - M. Long     - Initial version
!  15 Oct 2012 - R. Yantosca - Added ProTeX Headers, use F90 format/indents
!  22 Oct 2012 - R. Yantosca - Renamed to GIGC_Allocate_Interface

!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: AS

    ! Assume success
    RC = GIGC_SUCCESS

    ! UV albedo for photolysis
    ALLOCATE( UVALBEDO( IIPAR, JJPAR), STAT=RC  )

  END SUBROUTINE GIGC_Allocate_Interface
END MODULE GIGC_Initialization_Mod
#endif
