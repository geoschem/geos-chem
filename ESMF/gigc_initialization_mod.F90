#if defined( ESMF_ )
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
  PUBLIC :: GIGC_Get_Options
  PUBLIC :: GIGC_Init_Simulation
  PUBLIC :: GIGC_Init_Time_Interface
!
! !REVISION HISTORY: 
!  16 Oct 2012 - M. Long     - Initial version
!  16 Oct 2012 - R. Yantosca - Added ProTeX headers
!  22 Oct 2012 - R. Yantosca - Renamed to gigc_initialization_mod.F90
!  22 Oct 2012 - R. Yantosca - Renamed several routines for better consistency
!  25 Oct 2012 - R. Yantosca - Remove routine GIGC_SetEnv
!  03 Dec 2012 - R. Yantosca - Now pass extra arguments to GIGC_Init_Dimensions
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
! !IROUTINE: gigc_get_options
!
! !DESCRIPTION: Routine GIGC\_GET\_OPTIONS reads options for a GEOS-Chem 
!  simulation from the input.geos\_\_\_.rc input file.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE GIGC_Get_Options( am_I_Root, lonCtr, latCtr, Input_Opt, RC )
!
! !USES:
!
    USE CMN_GCTM_Mod       
    USE CMN_SIZE_Mod,       ONLY : dLat
    USE CMN_SIZE_Mod,       ONLY : dLon
    USE GIGC_ErrCode_Mod
    USE GIGC_Input_Opt_Mod, ONLY : OptInput
    USE Input_Mod,          ONLY : Read_Input_File
!
! !INPUT PARAMETERS: 
!
    LOGICAL,        INTENT(IN)    :: am_I_Root   ! Are we on the root CPU?
    REAL*4,         INTENT(IN)    :: lonCtr(:,:)
    REAL*4,         INTENT(IN)    :: latCtr(:,:)
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(INOUT) :: Input_Opt   ! Input Options object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure
! 
! !REMARKS:
!  NOTE: For now assume that GEOS_Chem will always accept a regular 
!  Cartesian grid.  This is more or less dictated by the input data.
!  The GEOS-5 data can be regridded via ESMF from whatever grid it uses.
!  (bmy, 11/30/12)
!
! !REVISION HISTORY: 
!  15 Oct 2012 - M. Long     - Initial version
!  15 Oct 2012 - R. Yantosca - Added ProTeX Headers, use F90 format/indents
!  22 Oct 2012 - R. Yantosca - Renamed to GIGC_Get_Options
!  22 Oct 2012 - R. Yantosca - Added RC output argument
!  01 Nov 2012 - R. Yantosca - Now pass the Input Options object via arg list
!  03 Dec 2012 - R. Yantosca - Reorder subroutines for clarity
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    REAL*8 :: deltaLon, deltaLat
  
    ! Assume success
    RC       = GIGC_SUCCESS

    ! Assume the longitude difference is the same for all grid boxes 
    deltaLon = Roundoff( ( DBLE( lonCtr(2,1) ) / PI_180 ), 4 )  &
             - Roundoff( ( DBLE( lonCtr(1,1) ) / PI_180 ), 4 )

    ! Assume the latitude difference is the same for all grid boxes
    deltaLat = Roundoff( ( DBLE( latCtr(1,2) ) / PI_180 ), 4 )  &
             - Roundoff( ( DBLE( latCtr(1,1) ) / PI_180 ), 4 )

    ! Save into arrays of CMN_SIZE
    dLon     = deltaLon
    dLat     = deltaLat

    ! Read the GEOS-Chem input file here
    CALL Read_Input_File( am_I_Root, Input_Opt, RC )

  END SUBROUTINE GIGC_Get_Options
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
  SUBROUTINE GIGC_Init_Simulation( am_I_Root,       tsChem,          &
                                   nymd,            nhms,            &
                                   lonCtr,          latCtr,          &      
!                                   latEdg,
                                   value_I_LO,      &
                                   value_J_LO,      value_I_HI,      &
                                   value_J_HI,      value_IM,        &
                                   value_JM,        value_LM,        &
                                   value_IM_WORLD,  value_JM_WORLD,  &
                                   value_LM_WORLD,  Input_Opt,       &
                                   State_Chm,       State_Met,       &
                                   mapping,         RC              )      
!
! !USES:
!
    USE GIGC_Environment_Mod
    USE GIGC_ErrCode_Mod  
    USE GIGC_Input_Opt_Mod
    USE GIGC_State_Chm_Mod
    USE GIGC_State_Met_Mod
    USE CMN_GCTM_Mod
    USE CMN_SIZE_MOD
    USE COMODE_MOD
    USE COMODE_LOOP_MOD       
    USE GCKPP_COMODE_MOD,     ONLY : INIT_GCKPP_COMODE
    USE ERROR_MOD,            ONLY : DEBUG_MSG
    USE GRID_MOD,             ONLY : AREA_M2
    USE GRID_MOD,             ONLY : INIT_GRID
    USE GRID_MOD,             ONLY : XEDGE
    USE GRID_MOD,             ONLY : XMID
    USE GRID_MOD,             ONLY : YEDGE
    USE GRID_MOD,             ONLY : YMID
    USE Mapping_Mod,          ONLY : MapWeight
    USE Mapping_Mod,          ONLY : Init_Mapping
    USE Olson_Landmap_Mod,    ONLY : Init_Olson_Landmap
    USE Olson_Landmap_Mod,    ONLY : Compute_Olson_Landmap
    USE Olson_Landmap_Mod,    ONLY : Cleanup_Olson_Landmap
    USE PBL_MIX_MOD,          ONLY : INIT_PBL_MIX
    USE PRESSURE_MOD,         ONLY : INIT_PRESSURE
    USE TRACER_MOD,           ONLY : ITS_A_FULLCHEM_SIM
    USE TRACER_MOD,           ONLY : ITS_AN_AEROSOL_SIM
    USE TRACERID_MOD,         ONLY : SETTRACE
    USE TOMS_MOD,             ONLY : TO3_DAILY
    USE WETSCAV_MOD,          ONLY : INIT_WETSCAV
    USE DRYDEP_MOD,           ONLY : INIT_WEIGHTSS
!
! !INPUT PARAMETERS: 
!
    LOGICAL,         INTENT(IN)    :: am_I_Root        ! Is this the root CPU?
    REAL,            INTENT(IN)    :: tsChem           ! Chemistry timestep [s]
    INTEGER,         INTENT(IN)    :: nymd             ! GMT date (YYYY/MM/DD)
    INTEGER,         INTENT(IN)    :: nhms             ! GMT time (hh:mm:ss)
    REAL*4,  TARGET, INTENT(IN)    :: lonCtr(:,:)      ! Lon centers [radians]
    REAL*4,  TARGET, INTENT(IN)    :: latCtr(:,:)      ! Lat centers [radians]
!    REAL*4,  TARGET, INTENT(IN)    :: latEdg(:,:)      ! Lat centers [radians]
    INTEGER,         INTENT(IN)    :: value_I_LO       ! Min local lon index
    INTEGER,         INTENT(IN)    :: value_J_LO       ! Min local lat index
    INTEGER,         INTENT(IN)    :: value_I_HI       ! Max local lon index
    INTEGER,         INTENT(IN)    :: value_J_HI       ! Max local lat index
    INTEGER,         INTENT(IN)    :: value_IM         ! Local # of lons
    INTEGER,         INTENT(IN)    :: value_JM         ! Local # of lats
    INTEGER,         INTENT(IN)    :: value_LM         ! Local # of levels
    INTEGER,         INTENT(IN)    :: value_IM_WORLD   ! Global # of lons
    INTEGER,         INTENT(IN)    :: value_JM_WORLD   ! Global # of lats
    INTEGER,         INTENT(IN)    :: value_LM_WORLD   ! Global # of levels
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(OptInput),  INTENT(INOUT) :: Input_Opt        ! Input Options
    TYPE(ChmState),  INTENT(INOUT) :: State_Chm        ! Chemistry State
    TYPE(MetState),  INTENT(INOUT) :: State_Met        ! Meteorology State
    TYPE(MapWeight), POINTER       :: mapping(:,:)     ! Olson mapping object
!
!
! !OUTPUT PARAMETERS:
!
    INTEGER,         INTENT(OUT)   :: RC               ! Success or failure?  
!
! !REMARKS
!  Add other calls to G EOS-Chem init routines as necessary.
!  NOTE: Later on maybe split these init calls among other routines.
!  Also need to add better error trapping
!
! !REVISION HISTORY: 
!  15 Oct 2012 - M. Long     - Initial version
!  15 Oct 2012 - R. Yantosca - Added ProTeX Headers, use F90 format/indents
!  17 Oct 2012 - R. Yantosca - Now initialize the chemistry mechanism
!  19 Oct 2012 - R. Yantosca - Now reference gigc_state_chm_mod.F90
!  19 Oct 2012 - R. Yantosca - Now reference gigc_state_met_mod.F90
!  01 Nov 2012 - R. Yantosca - Now reference gigc_input_opt_mod.F90
!  09 Nov 2012 - R. Yantosca - Now use fields from Input Options object
!  13 Nov 2012 - R. Yantosca - Pass Input Options object to routines 
!                              SETEMDEP, INIT_COMODE
!  28 Nov 2012 - R. Yantosca - Remove reference to INIT_DAO, since there are
!                              no more module arrays anymore in dao_mod.F
!  29 Nov 2012 - R. Yantosca - Add lonCtr, latCtr, latEdg as arguments
!  29 Nov 2012 - R. Yantosca - Now pass am_I_Root to Olson landmap routines
!  03 Dec 2012 - R. Yantosca - Now pass value_* arguments to pass dimension
!                              info from ESMF down to lower-level routines
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    LOGICAL :: prtDebug
    INTEGER :: DTIME, K, AS, N, YEAR, I, J, L

    !=======================================================================
    ! Initialize key GEOS-Chem sections
    !=======================================================================

    ! Initialize
    RC    = GIGC_SUCCESS
    DTIME = tsChem

    ! Initialize the timing routines
    CALL GIGC_Init_Time_Interface( DTIME, nymd, nhms, am_I_Root, RC )

    ! Allocate GEOS-Chem module arrays
    CALL GIGC_Allocate_All( am_I_Root,       Input_Opt,       &
                            RC,              value_I_LO,      &
                            value_J_LO,      value_I_HI,      &
                            value_J_HI,      value_IM,        &
                            value_JM,        value_LM,        &
                            value_IM_WORLD,  value_JM_WORLD,  &
                            value_LM_WORLD                   )

    ! Read options from the GEOS-Chem input file "input.geos"
    CALL GIGC_Get_Options( am_I_Root, lonCtr, latCtr, Input_Opt, RC )

    ! Determine if we have to print debug output
    prtDebug = ( Input_Opt%LPRT .and. am_I_Root )

    ! Initialize derived-type objects for meteorology & chemistry states
    CALL GIGC_Init_All( am_I_Root, Input_Opt, State_Chm, State_Met, RC )

    ! Save tracer names and ID's into State_Chm
    DO N = 1, Input_Opt%N_TRACERS
       State_Chm%TRAC_NAME(N) = 'TRC_' // TRIM( Input_Opt%TRACER_NAME(N) )
       State_Chm%TRAC_ID  (N) = Input_Opt%ID_TRACER(N)
    ENDDO
       
    ! Allocate and zero GEOS-Chem diagnostic arrays

    ! Initialize the GEOS-Chem pressure module (set Ap & Bp)
    CALL Init_Pressure( am_I_Root )

    ! Initialize the PBL mixing module
    CALL Init_PBL_Mix()

    ! Initialize arrays SO2s, H2O2s in wetscav_mod.F for use in sulfate chem
    CALL Init_WetScav( State_Met )

    !=======================================================================
    ! Initialize dry deposition 
    !=======================================================================
    IF ( Input_Opt%LDRYD )  THEN

       ! Initialize the derived type object containing
       ! mapping information for the MODIS LAI routines
       IF ( Input_Opt%USE_OLSON_2001 ) THEN
          CALL Init_Mapping( 1440, 720, IIPAR, JJPAR, mapping )
       ELSE
          CALL DEBUG_MSG( '### Init Mapping 05x05')
          CALL Init_Mapping(  720, 360, IIPAR, JJPAR, mapping )
       ENDIF

       ! Compute the Olson land types that occur in each grid box
       ! (i.e. this is a replacement for rdland.F and vegtype.global)
       CALL Init_Olson_Landmap   ( am_I_Root                     )
       CALL Compute_Olson_Landmap( am_I_Root, mapping, State_Met )
       CALL Cleanup_Olson_Landmap( am_I_Root                     )

       !### Debug
       IF ( prtDebug ) THEN
          CALL DEBUG_MSG( '### GIGC_INIT_CHEMISTRY: after OLSON' )
       ENDIF

!-----------------------------------------------------------------------------
! NOTE: Add this soon (bmy, 12/4/12)
!      ! Read drydep inputs from the netCDF file
!      ! Save Olson indices in INDOLSON array, in order to avoid
!      ! confusion w/ previously-assinged variable name IOLSON
!      CALL READ_DRYDEP_INPUTS( am_I_Root,                      &
!     &                            DRYCOEFF, INDOLSON, IDEP,       &
!     &                            IWATER,   NWATER,   IZO,        
!     &                            IDRYDEP,  IRI,      IRLU,     
!     &                            IRAC,     IRGSS,    IRGSO, 
!     &                            IRCLS,    IRCLO,    IVSMAX )
!
!      ! Calls INIT_WEIGHTSS to calculate the volume distribution of 
!      ! sea salt aerosols (jaegle 5/11/11)
       CALL INIT_WEIGHTSS()
!      FIRST = .FALSE.
!-----------------------------------------------------------------------------

    ENDIF

    !=======================================================================
    ! Initialize chemistry mechanism
    !=======================================================================

    ! Set some size variables
    NLAT   = JJPAR
    NLONG  = IIPAR
    NVERT  = IVERT 
    NPVERT = NVERT
    NPVERT = NVERT + IPLUME

    ! If we are doing chemistry ...
    IF ( Input_Opt%LCHEM ) THEN

       ! Initialize arrays in comode_mod.F
       CALL INIT_COMODE( am_I_Root, Input_Opt, RC )

       ! Initialize KPP (if necessary)
       IF ( Input_Opt%LKPP ) THEN
          CALL INIT_GCKPP_COMODE( am_I_Root, IIPAR,   JJPAR, LLTROP,  &
                                  ITLOOP,    NMTRATE, IGAS,  RC      )
       ENDIF
       
       ! Read from data file mglob.dat
       CALL READER( .TRUE., am_I_Root )

       !### Debug
       IF ( prtDebug ) THEN
          CALL DEBUG_MSG( '### GIGC_INIT_CHEMISTRY: after READER' )
       ENDIF

       ! Set NCS for urban chemistry only (since that is where we
       ! have defined the GEOS-CHEM mechanism) (bdf, bmy, 4/21/03)
       NCS = NCSURBAN

       !### Debug
       IF ( prtDebug ) THEN
          CALL DEBUG_MSG( '### GIGC_INIT_CHEMISTRY: after READER' )        
       ENDIF
       
       ! Redefine NTLOOP since READER defines it initially (bmy, 9/28/04)
       NLOOP   = NLAT  * NLONG
       NTLOOP  = NLOOP * NVERT
       NTTLOOP = NTLOOP
       
       ! Read "globchem.dat" chemistry mechanism
       CALL READCHEM( am_I_Root, Input_Opt, RC )
       
       !### Debug
       IF ( prtDebug ) THEN
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

       ! Get CH4 [ppbv] in 4 latitude bins for each year
       YEAR = NYMD / 10000
       CALL GET_GLOBAL_CH4( YEAR,   .TRUE., C3090S, &
                            C0030S, C0030N, C3090N, am_I_Root )

       !### Debug
       IF ( prtDebug ) THEN
          CALL DEBUG_MSG( '### GIGC_INIT_CHEMISTRY: after GET_GLOBAL_CH4' )
       ENDIF

       ! Initialize FAST-J photolysis
       CALL INPHOT( LLPAR, NPHOT, am_I_Root ) 
       
       !### Debug
       IF ( prtDebug ) THEN
          CALL DEBUG_MSG( '### GIGC_INIT_CHEMISTRY: after INPHOT' )        
       ENDIF
       
       ! Flag certain chemical species
       CALL SETTRACE( am_I_Root, Input_Opt, State_Chm, RC )
       
       !### Debug
       IF ( prtDebug ) THEN
          CALL DEBUG_MSG( '### GIGC_INIT_CHEMISTRY: after SETTRACE' )
       ENDIF
       
       ! Flag emission & drydep rxns
       CALL SETEMDEP( am_I_Root, Input_Opt, RC )
       
       !### Debug
       IF ( prtDebug ) THEN
          CALL DEBUG_MSG( '### GIGC_INIT_CHEMISTRY: after SETEMDEP' )
       ENDIF

       ! Initialize dry deposition (work in progress), add here
       ! Allocate array of overhead O3 columns for TOMS
       ALLOCATE( TO3_DAILY( IIPAR, JJPAR ), STAT=AS )
       TO3_DAILY = 0d0

    ENDIF

  END SUBROUTINE GIGC_Init_Simulation

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
! Prior to 12/3/12:
! NOTE: UVALBEDO is now contained in State_Met (bmy, 12/3/12)
!!EOC
!!------------------------------------------------------------------------------
!!          Harvard University Atmospheric Chemistry Modeling Group            !
!!------------------------------------------------------------------------------
!!BOP
!!
!! !IROUTINE: gigc_allocate_interface
!!
!! !DESCRIPTION: Routine GIGC\_ALLOCATE\_INTERFACE allocates GEOS-Chem arrays
!!  for a Grid-Independent GEOS-Chem (aka "GIGC") simulation.
!!\\
!!\\
!! !INTERFACE:
!!
!  SUBROUTINE GIGC_Allocate_Interface( am_I_Root, RC )
!!
!! !USES:
!!
!    USE CMN_SIZE_MOD,     ONLY : IIPAR, JJPAR, LLPAR
!    USE GIGC_ErrCode_Mod
!    USE ERROR_MOD,        ONLY : ALLOC_ERR
!    USE UVALBEDO_MOD,     ONLY : UVALBEDO
!   !USE DAO_MOD,          ONLY : TO3
!!
!! !INPUT PARAMETERS: 
!!
!    LOGICAL, INTENT(IN)  :: am_I_Root   ! Are we on the root CPU?
!!
!! !OUTPUT PARAMETERS:
!!
!    INTEGER, INTENT(OUT) :: RC          ! Success or failure
!! 
!! !REVISION HISTORY: 
!!  15 Oct 2012 - M. Long     - Initial version
!!  15 Oct 2012 - R. Yantosca - Added ProTeX Headers, use F90 format/indents
!!  22 Oct 2012 - R. Yantosca - Renamed to GIGC_Allocate_Interface
!
!!EOP
!!------------------------------------------------------------------------------
!!BOC
!!
!! !LOCAL VARIABLES:
!!
!    INTEGER :: AS
!
!    ! Assume success
!    RC = GIGC_SUCCESS
!
!    ! UV albedo for photolysis
!    ALLOCATE( UVALBEDO( IIPAR, JJPAR), STAT=RC  )
!
!  END SUBROUTINE GIGC_Allocate_Interface
!EOC
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 910.1 and      !
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: RoundOff
!
! !DESCRIPTION: Rounds a number X to N decimal places of precision.
!\\
!\\
! !INTERFACE:
!
  FUNCTION RoundOff( X, N ) RESULT( Y )
!
! !INPUT PARAMETERS:
! 
    REAL*8,  INTENT(IN) :: X   ! Number to be rounded
    INTEGER, INTENT(IN) :: N   ! Number of decimal places to keep
!
! !RETURN VALUE:
!
    REAL*8              :: Y   ! Number rounded to N decimal places
!
! !REMARKS:
!  The algorithm to round X to N decimal places is as follows:
!  (1) Multiply X by 10**(N+1)
!  (2) If X < 0, then add -5 to X; otherwise add 5 to X
!  (3) Round X to nearest integer
!  (4) Divide X by 10**(N+1)
!  (5) Truncate X to N decimal places: INT( X * 10**N ) / 10**N
!                                                                             .
!  Rounding algorithm from: Hultquist, P.F, "Numerical Methods for Engineers 
!   and Computer Scientists", Benjamin/Cummings, Menlo Park CA, 1988, p. 20.
!                                                                             .
!  Truncation algorithm from: http://en.wikipedia.org/wiki/Truncation
!                                                                             .
!  The two algorithms have been merged together for efficiency.
!
! !REVISION HISTORY:
!  14 Jul 2010 - R. Yantosca - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES
!
    ! Round and truncate X to N decimal places
    Y = INT( NINT( X*(10d0**(N+1)) + SIGN( 5d0, X ) ) / 10d0 ) / (10d0**N)

  END FUNCTION RoundOff
!EOC
END MODULE GIGC_Initialization_Mod
#endif
