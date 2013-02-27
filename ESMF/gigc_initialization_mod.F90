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
  PUBLIC  :: GIGC_Get_Options
  PUBLIC  :: GIGC_Init_Simulation
!
! !PRIVATE MEMBER FUNCTIONS:
!
  PRIVATE :: RoundOff
!
! !REVISION HISTORY: 
!  16 Oct 2012 - M. Long     - Initial version
!  16 Oct 2012 - R. Yantosca - Added ProTeX headers
!  22 Oct 2012 - R. Yantosca - Renamed to gigc_initialization_mod.F90
!  22 Oct 2012 - R. Yantosca - Renamed several routines for better consistency
!  25 Oct 2012 - R. Yantosca - Remove routine GIGC_SetEnv
!  03 Dec 2012 - R. Yantosca - Now pass extra arguments to GIGC_Init_Dimensions
!  06 Dec 2012 - R. Yantosca - Now remove routine GIGC_Init_TimeInterface; this
!                              is now superseded by Accept_Date_Time_From_ESMF
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
  SUBROUTINE GIGC_Get_Options( am_I_Root, lonCtr,    latCtr,  &
                               Input_Opt, State_Chm, RC      )
!
! !USES:
!
    USE CMN_GCTM_Mod       
    USE CMN_SIZE_Mod
    USE GIGC_ErrCode_Mod
    USE GIGC_Input_Opt_Mod, ONLY : OptInput
    USE GIGC_State_Chm_Mod, ONLY : ChmState
    USE Input_Mod,          ONLY : Read_Input_File
    USE Input_Mod,          ONLY : Initialize_Geos_Grid
!
! !INPUT PARAMETERS: 
!
    LOGICAL,        INTENT(IN)    :: am_I_Root     ! Are we on the root CPU?
    REAL*4,         INTENT(IN)    :: lonCtr(:,:)   ! Lon ctrs [deg] from ESMF
    REAL*4,         INTENT(IN)    :: latCtr(:,:)   ! Lat ctrs [deg] from ESMF
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(INOUT) :: Input_Opt   ! Input Options object
    TYPE(ChmState), INTENT(INOUT) :: State_Chm
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC            ! Success or failure
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
!  07 Dec 2012 - R. Yantosca - Compute DLON, DLAT more rigorously
!  26 Feb 2013 - M. Long     - Now pass State_Chm as an argument
!  26 Feb 2013 - M. Long     - Read "input.geos" on root CPU only
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER :: I, J, L
  
    ! Assume success
    RC = GIGC_SUCCESS

    !========================================================================
    ! Compute the DLON and DLAT values.  NOTE, this is a kludge since to do
    ! this truly rigorously, we should take the differences between the grid
    ! box edges.  But because I can't seem to find a way to get the grid
    ! box edge information, the next best thing is to take the differences
    ! between the grid box centers.  They should all be the same given that
    ! the GEOS-Chem grid is regular. (bmy, 12/7/12)
    !========================================================================
    DO L = 1, LLPAR
    DO J = 1, JJPAR
    DO I = 1, IIPAR

       ! Compute Delta-Longitudes [degrees]
       IF ( I == IIPAR ) THEN
          dLon(I,J,L) = RoundOff( ( DBLE( lonCtr(IIPAR,  J) ) / PI_180 ), 4 ) &
                      - RoundOff( ( DBLE( lonCtr(IIPAR-1,J) ) / PI_180 ), 4 )
       ELSE
          dLon(I,J,L) = RoundOff( ( DBLE( lonCtr(I+1,    J) ) / PI_180 ), 4 ) &
                      - RoundOff( ( DBLE( lonCtr(I,      J) ) / PI_180 ), 4 )
       ENDIF

       ! Compute Delta-Latitudes [degrees]
       IF ( J == JJPAR ) THEN
          dLat(I,J,L) = RoundOff( ( DBLE( latCtr(I,JJPAR  ) ) / PI_180 ), 4 ) &
                      - RoundOff( ( DBLE( latCtr(I,JJPAR-1) ) / PI_180 ), 4 )
       ELSE
          dLat(I,J,L) = RoundOff( ( DBLE( latCtr(I,J+1    ) ) / PI_180 ), 4 ) &
                      - RoundOff( ( DBLE( latCtr(I,J      ) ) / PI_180 ), 4 )
       ENDIF

    ENDDO
    ENDDO
    ENDDO

    ! Read the GEOS-Chem input file here.  For now only read on the root
    ! CPU so that we can broadcast to other CPUs in GIGC_Init_Simulation
    ! (mlong, bmy, 2/26/13)
    IF ( am_I_Root ) THEN
       CALL Read_Input_File( am_I_Root, Input_Opt, RC )
    ENDIF

    ! We still need to call Initialize_Geos_Grid on all CPUs though.
    ! without having to read the "input.geos" file. (mlong, bmy, 2/26/13)
    IF ( .not. am_I_Root ) THEN
       CALL Initialize_Geos_Grid( am_I_Root, RC )
    ENDIF

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
  SUBROUTINE GIGC_Init_Simulation( am_I_Root,                        &
                                   nymdB,           nhmsB,           &
                                   nymdE,           nhmsE,           &
                                   tsChem,          tsDyn,           &
                                   lonCtr,          latCtr,          &      
                                   value_I_LO,      value_J_LO,      &
                                   value_I_HI,      value_J_HI,      &
                                   value_IM,        value_JM,        &
                                   value_LM,        value_IM_WORLD,  &
                                   value_JM_WORLD,  value_LM_WORLD,  &
                                   Input_Opt,       State_Chm,       &
                                   State_Met,       mapping,         &
                                   RC                               )      
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
    USE GRID_MOD,             ONLY : INIT_GRID
    USE Mapping_Mod,          ONLY : MapWeight
    USE Mapping_Mod,          ONLY : Init_Mapping
    USE Olson_Landmap_Mod,    ONLY : Init_Olson_Landmap
    USE Olson_Landmap_Mod,    ONLY : Compute_Olson_Landmap
    USE Olson_Landmap_Mod,    ONLY : Cleanup_Olson_Landmap
    USE PBL_MIX_MOD,          ONLY : INIT_PBL_MIX
    USE PRESSURE_MOD,         ONLY : INIT_PRESSURE
    USE TRACER_MOD,           ONLY : ITS_A_FULLCHEM_SIM
    USE TRACER_MOD,           ONLY : ITS_AN_AEROSOL_SIM
    USE TRACER_MOD,           ONLY : INIT_TRACER
    USE TRACERID_MOD,         ONLY : SETTRACE
    USE TOMS_MOD,             ONLY : TO3_DAILY
    USE WETSCAV_MOD,          ONLY : INIT_WETSCAV
    USE DRYDEP_MOD,           ONLY : INIT_WEIGHTSS, INIT_DRYDEP
    USE DUST_MOD,             ONLY : INIT_DUST
    USE GIGC_MPI_WRAP
    ! KLUGE - MSL 04 Jan 2013
    USE LOGICAL_MOD,          ONLY : LVARTROP
    USE TIME_MOD,             ONLY : SET_TIMESTEPS
!
! !INPUT PARAMETERS: 
!
    LOGICAL,         INTENT(IN)    :: am_I_Root       ! Is this the root CPU?
    INTEGER,         INTENT(IN)    :: nymdB           ! GMT date (YYYY/MM/DD)
    INTEGER,         INTENT(IN)    :: nhmsB           ! GMT time (hh:mm:ss)
    INTEGER,         INTENT(IN)    :: nymdE           ! GMT date (YYYY/MM/DD)
    INTEGER,         INTENT(IN)    :: nhmsE           ! GMT time (hh:mm:ss)
    REAL,            INTENT(IN)    :: tsChem          ! Chem timestep [seconds]
    REAL,            INTENT(IN)    :: tsDyn           ! Dyn  timestep [seconds]
    REAL*4,  TARGET, INTENT(IN)    :: lonCtr(:,:)     ! Lon centers [radians]
    REAL*4,  TARGET, INTENT(IN)    :: latCtr(:,:)     ! Lat centers [radians]
    INTEGER,         INTENT(IN)    :: value_I_LO      ! Min local lon index
    INTEGER,         INTENT(IN)    :: value_J_LO      ! Min local lat index
    INTEGER,         INTENT(IN)    :: value_I_HI      ! Max local lon index
    INTEGER,         INTENT(IN)    :: value_J_HI      ! Max local lat index
    INTEGER,         INTENT(IN)    :: value_IM        ! # lons on this CPU
    INTEGER,         INTENT(IN)    :: value_JM        ! # lats on this CPU
    INTEGER,         INTENT(IN)    :: value_LM        ! # levs on this CPU
    INTEGER,         INTENT(IN)    :: value_IM_WORLD  ! # lons in whole globe
    INTEGER,         INTENT(IN)    :: value_JM_WORLD  ! # lats in whole globe
    INTEGER,         INTENT(IN)    :: value_LM_WORLD  ! # levs in whole globe
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(OptInput),  INTENT(INOUT) :: Input_Opt       ! Input Options
    TYPE(ChmState),  INTENT(INOUT) :: State_Chm       ! Chemistry State
    TYPE(MetState),  INTENT(INOUT) :: State_Met       ! Meteorology State
    TYPE(MapWeight), POINTER       :: mapping(:,:)    ! Olson mapping object
!
!
! !OUTPUT PARAMETERS:
!
    INTEGER,         INTENT(OUT)   :: RC              ! Success or failure?  
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
!  06 Dec 2012 - R. Yantosca - Now accept start & end dates & times via 
!                              the nymdB, nymdE, nhmsB, nhmsE arguments
!  06 Dec 2012 - R. Yantosca - Remove nymd, nhms arguments, these will be
!                              the same as nymdB, nhmsB (start of run)
!  26 Feb 2013 - M. Long     - Now read ASCII input files on root CPU and 
!                              broadcast to other CPUs.
!  26 Feb 2013 - R. Yantosca - Cosmetic changes
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    LOGICAL :: prtDebug
    INTEGER :: DTIME, K, AS, N, YEAR, I, J, L, TMP

    !=======================================================================
    ! Initialize key GEOS-Chem sections
    !=======================================================================

    ! Initialize
    RC    = GIGC_SUCCESS
    DTIME = tsChem

    if (am_I_Root) write(*,*) 'GRID :: JM_WORLD :: ', value_JM_WORLD
    if (am_I_Root) write(*,*) 'GRID :: IM_WORLD :: ', value_IM_WORLD

    ! Allocate GEOS-Chem module arrays
    CALL GIGC_Allocate_All( am_I_Root,       Input_Opt,       &
                            RC,              value_I_LO,      &
                            value_J_LO,      value_I_HI,      &
                            value_J_HI,      value_IM,        &
                            value_JM,        value_LM,        &
                            value_IM_WORLD,  value_JM_WORLD,  &
                            value_LM_WORLD                   )

    ! Save timing fields in Input_Opt for passing down to module
    ! GeosCore/input_mod.F via routine GIGC_Get_Options (bmy, 12/6/12)
    Input_Opt%NYMDb   = nymdB
    Input_Opt%NHMSb   = nhmsB
    Input_Opt%NYMDe   = nymdE
    Input_Opt%NHMSe   = nhmsE
    Input_Opt%TS_CHEM = INT( tsChem ) / 60   ! Chemistry timestep [min]
    Input_Opt%TS_DYN  = INT( tsDyn  ) / 60   ! Dynamic   timestep [mn]

    ! Read options from the GEOS-Chem input file "input.geos"
    ! And initialize grid
    CALL GIGC_Get_Options( am_I_Root, lonCtr, latCtr, Input_Opt, &
                           State_Chm, RC                         )

    !-----------------------------------------------------------------------
    ! Read "input.geos" on the root CPU and broadcast to other CPUs
    !-----------------------------------------------------------------------

    ! Broadcast "input.geos" options values to all threads with MPI
    CALL GIGC_Input_Bcast( Input_Opt, RC )
    CALL GIGC_IDT_Bcast(   Input_Opt, RC )

    ! Complete initialization ops on all threads
    IF ( .not. am_I_Root ) THEN 

       ! Initialize dry deposition
       CALL INIT_DRYDEP( am_I_Root, Input_Opt, RC )

       ! Initialize tracer quantities
       CALL INIT_TRACER( am_I_Root, Input_Opt, RC )

       ! Working Kluge - MSL; Break this & Fix the result...
       LVARTROP = Input_Opt%LVARTROP 

       ! Set GEOS-Chem timesteps
       CALL SET_TIMESTEPS( am_I_Root,                                  &
                           Input_Opt%TS_CHEM,                          &
                           Input_Opt%TS_CONV,                          &
                           Input_Opt%TS_DYN,                           &
                           Input_Opt%TS_EMIS,                          &
                           MAX( Input_Opt%TS_DYN, Input_Opt%TS_CONV ), &
                           Input_Opt%TS_DIAG )
       
    ENDIF

    !-----------------------------------------------------------------------
    ! Read other ASCII files on the root CPU and broadcast to other CPUs
    !-----------------------------------------------------------------------

    ! Read "mglob.dat"
    IF ( am_I_Root ) THEN

       ! Read from data file mglob.dat
       CALL READER( .TRUE., am_I_Root )

       !### Debug
       IF ( prtDebug ) THEN
          CALL DEBUG_MSG( '### GIGC_INIT_CHEMISTRY: after READER' )
       ENDIF
    ENDIF

    ! Broadcast "mglob.dat"
    CALL GIGC_Reader_Bcast( RC )

    ! Read "globchem.dat" chemistry mechanism
    ! NOTE: for now, read on all CPUs and fix later (bmy, mlong, 2/26/13)
!    IF ( am_I_Root ) THEN
       CALL READCHEM( am_I_Root, Input_Opt, RC )
       
       !### Debug
       IF ( prtDebug ) THEN
          CALL DEBUG_MSG( '### GIGC_INIT_CHEMISTRY: after READCHEM' )        
       ENDIF
!    ENDIF
!
!    ! Broadcast "globchem.dat" to other CPUs
!    CALL GIGC_ReadChem_Bcast( RC )

    ! Initialize FAST-J photolysis
    ! NOTE: for now, read on all CPUs and fix later (bmy, mlong, 2/26/13)   
!    IF ( am_I_Root ) THEN
       CALL INPHOT( LLPAR, NPHOT, Input_Opt, am_I_Root ) 
       
       !### Debug
       IF ( prtDebug ) THEN
          CALL DEBUG_MSG( '### GIGC_INIT_CHEMISTRY: after INPHOT' )        
       ENDIF
!    ENDIF
!
!    Broadcast FAST-J inputs to other CPUs
!    CALL GIGC_Inphot_Bcast(  Input_Opt, RC )

    !-----------------------------------------------------------------------
    ! Continue with GEOS-Chem setup
    !-----------------------------------------------------------------------

    ! Zero diagnostic arrays
    CALL Initialize( 2, am_I_Root )

    ! Zero diagnostic counters
    CALL Initialize( 3, am_I_Root )

    ! Determine if we have to print debug output
    prtDebug = ( Input_Opt%LPRT .and. am_I_Root )

    ! Initialize derived-type objects for meteorology & chemistry states
    CALL GIGC_Init_All( am_I_Root, Input_Opt, State_Chm, State_Met, RC )

    ! Save tracer names and ID's into State_Chm
    DO N = 1, Input_Opt%N_TRACERS
       State_Chm%TRAC_NAME(N) = 'TRC_' // TRIM( Input_Opt%TRACER_NAME(N) )
       State_Chm%TRAC_ID  (N) = Input_Opt%ID_TRACER(N)
    ENDDO
       
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
          CALL Init_Mapping(  720, 360, IIPAR, JJPAR, mapping )
       ENDIF

       ! Compute the Olson land types that occur in each grid box
       ! (i.e. this is a replacement for rdland.F and vegtype.global)
       CALL Init_Olson_Landmap   ( am_I_Root, Input_Opt%DATA_DIR_1x1 )
       CALL Compute_Olson_Landmap( am_I_Root, mapping, State_Met     )
       CALL Cleanup_Olson_Landmap( am_I_Root                         )

       !### Debug
       IF ( prtDebug ) THEN
          CALL DEBUG_MSG( '### GIGC_INIT_CHEMISTRY: after OLSON' )
       ENDIF
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
       YEAR = NYMDb / 10000
       CALL GET_GLOBAL_CH4( YEAR,   .TRUE., C3090S, &
                            C0030S, C0030N, C3090N, am_I_Root )

       !### Debug
       IF ( prtDebug ) THEN
          CALL DEBUG_MSG( '### GIGC_INIT_CHEMISTRY: after GET_GLOBAL_CH4' )
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
