!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: input_mod.F90
!
! !DESCRIPTION: Module INPUT\_MOD contains routines that read the GEOS-Chem
!  input file at the start of the run and pass the information to several
!  other GEOS-Chem F90 modules.
!\\
!\\
! !INTERFACE:
!
MODULE INPUT_MOD
!
! !USES:
!
  USE CharPak_Mod,   ONLY : MAXDIM  => MAXSTRLEN
  USE inquireMod,    ONLY : findFreeLUN
  USE PRECISION_MOD       ! For GEOS-Chem Precision (fp)

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC  :: Read_Input_File
  PUBLIC  :: Do_Error_Checks
  PUBLIC  :: Validate_Directories
!
! !REVISION HISTORY:
!  20 Jul 2004 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !PRIVATE TYPES:
!
  INTEGER            :: IU_GEOS
  INTEGER, PARAMETER :: FIRSTCOL = 26
  INTEGER            :: CT1      = 0
  CHARACTER(LEN=255) :: FILENAME = 'input.geos'

CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_input_file
!
! !DESCRIPTION: Subroutine READ\_INPUT\_FILE is the driver program for
!  reading the GEOS-Chem input file "input.geos" from disk.
!\\
!\\
! In an ESMF environment, all time steps (chemistry, convection, emissions,
! dynamics) must be specified externally before calling this routine. This is
! done in routine GIGC\_Init\_Simulation (gigc\_initialization\_mod.F90).
! The time steps specified in input.geos are ignored.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE READ_INPUT_FILE( Input_Opt, State_Grid, RC )
!
! !USES:
!
    USE CHARPAK_MOD,        ONLY : ReadOneLine, STRREPL
    USE ErrCode_Mod
    USE FILE_MOD,           ONLY : IOERROR
    USE Input_Opt_Mod,      ONLY : OptInput
    USE State_Grid_Mod,     ONLY : GrdState
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(INOUT) :: Input_Opt   ! Input options
    TYPE(GrdState), INTENT(INOUT) :: State_Grid  ! Grid State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure
!
! !REVISION HISTORY:
!  20 Jul 2004 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    LOGICAL               :: EOF
    INTEGER               :: IOS
    CHARACTER(LEN=1)      :: TAB   = ACHAR(9)
    CHARACTER(LEN=1)      :: SPACE = ' '
    CHARACTER(LEN=MAXDIM) :: LINE
    CHARACTER(LEN=255)    :: ErrMsg, ThisLoc

    !=================================================================
    ! READ_INPUT_FILE begins here!
    !=================================================================

    ! Echo output
    IF ( Input_Opt%amIRoot ) THEN
       WRITE( 6, '(a  )' ) REPEAT( '=', 79 )
       WRITE( 6, '(a,/)' ) 'G E O S - C H E M   U S E R   I N P U T'
       WRITE( 6, 100   ) TRIM( FILENAME )
100    FORMAT( 'READ_INPUT_FILE: Opening ./', a )
    ENDIF

    ! Find a free file LUN
    IU_GEOS  = findFreeLUN()

    ! Assume success
    RC       = GC_SUCCESS

    ! For error handling
    ErrMsg   = ''
    ThisLoc  = ' -> at Read_Input_File (in module GeosCore/input_mod.F90)'

    ! Open file
    OPEN( IU_GEOS, FILE=TRIM( FILENAME ), STATUS='OLD', IOSTAT=IOS )
    IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_GEOS, 'read_input_file:1' )

    ! Loop until EOF
    DO

       ! Read a line from the file, exit if EOF
       LINE = ReadOneLine( IU_GEOS, EOF, IOS )
       IF ( EOF     ) EXIT
       IF ( IOS > 0 ) CALL IOERROR( IOS, IU_GEOS, 'read_input_file:2')

       ! Replace tab characters in LINE (if any) w/ spaces
       CALL STRREPL( LINE, TAB, SPACE )

       !=============================================================
       ! Call individual subroutines to read sections of the file
       !
       ! NOTE: You are pretty flexible in setting the order of the
       ! menus in the input file; however, a few guidelines apply:
       !
       ! (1) SIMULATION MENU should be listed first.
       ! (2) TIMESTEP MENU should be listed second.
       ! (3) ADVECTED SPECIES MENU should be listed third.
       ! (4) EMISSIONS, AEROSOL, CHEMISTRY, TRANSPORT, CONVECTION,
       !      and DEPOSITION menus (in any order) should follow.
       ! (5) Diagnostic menus, including OUTPUT, DIAGNOSTIC,
       !      PLANEFLIGHT, ND48, ND49, ND50, ND51, and PROD_LOSS
       !      menus (in any order) should follow next.
       ! (6) The following menus have no other restriction and
       !      can be placed anywhere (but by convention we will
       !      place them after the diagnostic menu): NESTED GRID
       !      UNIX CMDS, ARCHIVED OH, and O3PL menus.
       !=============================================================
       IF      ( INDEX( LINE, 'SIMULATION MENU'  ) > 0 ) THEN
          CALL READ_SIMULATION_MENU( Input_Opt, RC )
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = 'Error in "Read_Simulation_Menu"!'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF

#if !(defined( EXTERNAL_GRID ) || defined( EXTERNAL_FORCING ))
       ELSE IF ( INDEX( LINE, 'GRID MENU' ) > 0 ) THEN
          CALL READ_GRID_MENU( Input_Opt, State_Grid, RC )
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = 'Error in "Read_Grid_Menu"!'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF
#endif

       ELSE IF ( INDEX( LINE, 'TIMESTEP MENU' ) > 0 ) THEN
          CALL READ_TIMESTEP_MENU( Input_Opt, State_Grid, RC )
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = 'Error in "Read_Timestep_Menu"!'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF

       ELSE IF ( INDEX( LINE, 'ADVECTED SPECIES MENU' ) > 0 ) THEN
          CALL READ_ADVECTED_SPECIES_MENU( Input_Opt, RC )
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = 'Error in "Read_Advected_Species_Menu"!'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF

       ELSE IF ( INDEX( LINE, 'AEROSOL MENU'     ) > 0 ) THEN
          CALL READ_AEROSOL_MENU( Input_Opt, RC )
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = 'Error in "Read_Aerosol_Menu"!'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF

       ELSE IF ( INDEX( LINE, 'EMISSIONS MENU'   ) > 0 ) THEN
          CALL READ_EMISSIONS_MENU( Input_Opt, RC )
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = 'Error in "Read_Emissions_Menu"!'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF

       ELSE IF ( INDEX( LINE, 'CHEMISTRY MENU'   ) > 0 ) THEN
          CALL READ_CHEMISTRY_MENU( Input_Opt, RC )
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = 'Error in "Read_Chemistry_Menu"!'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF

       ELSE IF ( INDEX( LINE, 'TRANSPORT MENU'   ) > 0 ) THEN
          CALL READ_TRANSPORT_MENU( Input_Opt, RC )
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = 'Error in "Read_Transport_Menu"!'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF

       ELSE IF ( INDEX( LINE, 'PHOTOLYSIS MENU'   ) > 0 ) THEN
          CALL READ_PHOTOLYSIS_MENU( Input_Opt, RC )
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = 'Error in "Read_Photolysis_Menu"!'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF

       ELSE IF ( INDEX( LINE, 'CONVECTION MENU'  ) > 0 ) THEN
          CALL READ_CONVECTION_MENU( Input_Opt, RC )
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = 'Error in "Read_Convection_Menu"!'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF

       ELSE IF ( INDEX( LINE, 'DEPOSITION MENU'  ) > 0 ) THEN
          CALL READ_DEPOSITION_MENU( Input_Opt, RC )
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = 'Error in "Read_Depositon_Menu"!'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF

       ELSE IF ( INDEX( LINE, 'RADIATION MENU'   ) > 0 ) THEN
          CALL READ_RADIATION_MENU( Input_Opt, RC )
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = 'Error in "Read_Radiation_Menu"!'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF

       ELSE IF ( INDEX( LINE, 'CO SIM MENU'     ) > 0 ) THEN
          CALL READ_CO_SIM_MENU( Input_Opt, RC )
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = 'Error in "Read_CO_Sim_Menu"!'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF

       ELSE IF ( INDEX( LINE, 'CO2 SIM MENU'     ) > 0 ) THEN
          CALL READ_CO2_SIM_MENU( Input_Opt, RC )
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = 'Error in "Read_CO2_Sim_Menu"!'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF

       ELSE IF ( INDEX( LINE, 'POPS MENU'        ) > 0 ) THEN
          CALL READ_POPS_MENU( Input_Opt, RC )
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = 'Error in "Read_POPS_Menu"!'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF

       ELSE IF ( INDEX( LINE, 'MERCURY MENU'     ) > 0 ) THEN
          CALL READ_MERCURY_MENU( Input_Opt, RC  )
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = 'Error in "Read_Mercury_Menu"!'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF

       ELSE IF ( INDEX( LINE, 'CH4 MENU'         ) > 0 ) THEN
          CALL READ_CH4_MENU( Input_Opt, RC )
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = 'Error in "Read_CH4_Menu"!'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF

       ELSE IF ( INDEX( LINE, 'PASSIVE SPECIES' ) > 0 ) THEN
          CALL READ_PASSIVE_SPECIES_MENU( Input_Opt, RC )
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = 'Error in "Read_Passive_Species_Menu"!'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF

#if !(defined( EXTERNAL_GRID ) || defined( EXTERNAL_FORCING ))
       ELSE IF ( INDEX( LINE, 'PLANEFLIGHT MENU' ) > 0 ) THEN
          CALL READ_PLANEFLIGHT_MENU( Input_Opt, RC )
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = 'Error in "Read_Planeflight_Menu"!'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF

#ifdef BPCH_DIAG
       !==============================================================
       ! Skip BPCH-related menus unless compiled with BPCH_DIAG=y
       ! Always skip BPCH-related menus for external ESMs
       !==============================================================
       ELSE IF ( INDEX( LINE, 'GAMAP MENU'       ) > 0 ) THEN
          CALL READ_GAMAP_MENU( Input_Opt, RC )
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = 'Error in "Read_Gamap_Menu"!'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF

       ELSE IF ( INDEX( LINE, 'OUTPUT MENU'      ) > 0 ) THEN
          CALL READ_OUTPUT_MENU( Input_Opt, RC )
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = 'Error in "Read_Output_Menu"!'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF

       ELSE IF ( INDEX( LINE, 'DIAGNOSTIC MENU'  ) > 0 ) THEN
          CALL READ_DIAGNOSTIC_MENU( Input_Opt, RC )
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = 'Error in "Read_Diagnostic_Menu"!'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF

       ELSE IF ( INDEX( LINE, 'ND51 MENU'        ) > 0 ) THEN
          CALL READ_ND51_MENU( Input_Opt, RC )
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = 'Error in "Read_ND51_Menu"!'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF

       ELSE IF ( INDEX( LINE, 'ND51b MENU'       ) > 0 ) THEN
          CALL READ_ND51b_MENU( Input_Opt, RC )
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = 'Error in "Read_ND51b_Menu"!'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF

#if defined TOMAS
       ELSE IF ( INDEX( LINE, 'PROD & LOSS MENU' ) > 0 ) THEN
          CALL READ_PROD_LOSS_MENU( Input_Opt, RC )
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = 'Error in "Read_Prod_Loss_Menu"!'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF
#endif
#endif
#endif

       ELSE IF ( INDEX( LINE, 'END OF FILE'      ) > 0 ) THEN
          EXIT

       ENDIF
    ENDDO

    ! Close input file
    CLOSE( IU_GEOS )

    !=================================================================
    ! Further error-checking and initialization
    !=================================================================

    ! Check GEOS-CHEM timesteps
    CALL CHECK_TIME_STEPS( Input_Opt, State_Grid, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error in "Check_Time_Steps"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

  END SUBROUTINE READ_INPUT_FILE
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: split_one_line
!
! !DESCRIPTION: Subroutine SPLIT\_ONE\_LINE reads a line from the input file
!  (via routine READ\_ONE\_LINE), and separates it into substrings.
!\\
!\\
!  SPLIT\_ONE\_LINE also checks to see if the number of substrings found is
!  equal to the number of substrings that we expected to find.  However, if
!  you don't know a-priori how many substrings to expect a-priori,
!  you can skip the error check.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE SPLIT_ONE_LINE( SUBSTRS, N_SUBSTRS, N_EXP, LOCATION, RC )
!
! !USES:
!
    USE Charpak_Mod, ONLY: ReadOneLine, StrSplit
    Use ErrCode_Mod
    USE File_Mod,    ONLY: IoError
!
! !INPUT PARAMETERS:
!
    INTEGER,            INTENT(IN)  :: N_EXP      ! Expected # of substrs
    CHARACTER(LEN=*),   INTENT(IN)  :: LOCATION   ! Name of calling routine
!
! !OUTPUT PARAMETERS:
!
    ! Array of substrings (separated by " ")
    CHARACTER(LEN=255), INTENT(OUT) :: SUBSTRS(MAXDIM) ! Substrings
    INTEGER,            INTENT(OUT) :: N_SUBSTRS       ! # of substrings
    INTEGER,            INTENT(OUT) :: RC              ! Success/failure?
!
! !REVISION HISTORY:
!  20 Jul 2004 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER               :: IOS
    LOGICAL               :: EOF

    ! Strings
    CHARACTER(LEN=10)     :: IOS_Str
    CHARACTER(LEN=255)    :: ErrMsg, ThisLoc
    CHARACTER(LEN=MAXDIM) :: LINE

    !=================================================================
    ! SPLIT_ONE_LINE begins here!
    !=================================================================

    ! Initialize
    RC      = GC_SUCCESS
    IOS     = 0
    EOF     = .FALSE.
    ErrMsg  = 'Error reading the "input.geos" file!'
    ThisLoc = ' -> at Split_One_Line (in module GeosCore/input_mod.F90)'

    !=================================================================
    ! Read a line from disk
    !=================================================================
    LINE = ReadOneLine( IU_GEOS, EOF, IOS )

    ! Trap premature end-of-file error
    IF ( EOF ) THEN
       ErrMsg = 'Unexpected end-of-file reading input.geos!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    ! Trap other I/O error conditions
    IF ( IOS > 0 ) THEN
       WRITE( IOS_Str, '(i10)' ) IOS
       ErrMsg = 'I/O error number: ' // TRIM( IOS_STR ) // &
                'encountered when readiang "input.geos"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    !=================================================================
    ! Split the lines between spaces -- start at column FIRSTCOL
    !=================================================================
    CALL STRSPLIT( LINE(FIRSTCOL:), ' ', SUBSTRS, N_SUBSTRS )

    ! Sometimes we don't know how many substrings to expect,
    ! if N_EXP is greater than MAXDIM, then skip the error check
    IF ( N_EXP < 0 ) RETURN

    ! Stop if we found the wrong # of substrings
    IF ( N_EXP /= N_SUBSTRS ) THEN
       ErrMsg= 'SPLIT_ONE_LINE: error at ' // TRIM( LOCATION )
       WRITE( 6, '(a)' ) TRIM( ErrMSg )
       WRITE( 6, 100   ) N_EXP, N_SUBSTRS
       WRITE( 6, '(a)' ) 'STOP in SPLIT_ONE_LINE (input_mod.F90)!'
       WRITE( 6, '(a)' ) REPEAT( '=', 79 )
100    FORMAT( 'Expected ',i2, ' substrs but found ',i3 )
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

  END SUBROUTINE SPLIT_ONE_LINE
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_simulation_menu
!
! !DESCRIPTION: Subroutine READ\_SIMULATION\_MENU reads the SIMULATION MENU
!  section of the GEOS-Chem input file.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE READ_SIMULATION_MENU( Input_Opt, RC )
!
! !USES:
!
    USE Charpak_Mod,        ONLY : To_UpperCase
    USE ErrCode_Mod
    USE Input_Opt_Mod,      ONLY : OptInput
    USE TIME_MOD
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(INOUT) :: Input_Opt   ! Input options
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure
!
! !REVISION HISTORY:
!  20 Jul 2004 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER            :: N, C
    REAL(fp)           :: JulianDateStart,  JulianDateEnd
#if defined( ESMF_ ) || defined( MODEL_ )
    INTEGER            :: H,       M,       S
    REAL(f4)           :: init_UTC
#endif

    ! Strings
    CHARACTER(LEN=6)   :: TimeStr
    CHARACTER(LEN=8)   :: DateStr
    CHARACTER(LEN=12)  :: Met
    CHARACTER(LEN=24)  :: Sim
    CHARACTER(LEN=255) :: ErrMsg, ThisLoc

    ! Arrays
    CHARACTER(LEN=255) :: SUBSTRS(MAXDIM)

    !=================================================================
    ! READ_SIMULATION_MENU begins here!
    !=================================================================

    ! Initialize
    RC      = GC_SUCCESS
    ErrMsg  = 'Error reading the "input.geos" file!'
    ThisLoc = ' -> at Read_Simulation_Menu (in GeosCore/input_mod.F90)'

    !-----------------------------------------------------------------
    ! Simulation start and end time
    !-----------------------------------------------------------------

    ! Start YYYYMMDD, HHMMSS
    CALL SPLIT_ONE_LINE( SUBSTRS, N, 2, 'NYMDb, NHMSb', RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    READ( SUBSTRS(1:N), * ) Input_Opt%NYMDb, Input_Opt%NHMSb

    ! Make sure the starting date NYMDb is valid
    IF ( .not. Valid_Date( Input_Opt%NYMDb ) ) THEN
       WRITE( DateStr, '(i8.8)' ) Input_Opt%NYMDb
       ErrMsg = 'Input%Opt%NYMDb = ' // DateStr        // &
                ' is not a valid calendar date!'       // &
                ' Please check your "input.geos" file.'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    ! Make sure the starting time NHMSb is valid
    IF ( .not. Valid_Time( Input_Opt%NHMSb ) ) THEN
       WRITE( TimeStr, '(i6.6)' ) Input_Opt%NHMSb
       ErrMsg = 'Input%Opt%NHMSb = ' // TimeStr        // &
                ' is not a valid clock time!'          // &
                ' Please check your "input.geos" file.'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    ! End YYYYMMDD, HHMMSS
    CALL SPLIT_ONE_LINE( SUBSTRS, N, 2, 'NYMDe, NHMSe', RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    READ( SUBSTRS(1:N), * ) Input_Opt%NYMDe, Input_Opt%NHMSe

    ! Make sure the starting date NYMDb is valid
    IF ( .not. Valid_Date( Input_Opt%NYMDe ) ) THEN
       WRITE( DateStr, '(i8.8)' ) Input_Opt%NYMDe
       ErrMsg = 'Input%Opt%NYMDe = ' // DateStr        // &
                ' is not a valid calendar date!'       // &
                ' Please check your "input.geos" file.'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    ! Make sure the ending time NHMSe is valid
    IF ( .not. Valid_Time( Input_Opt%NHMSe ) ) THEN
       WRITE( TimeStr, '(i6.6)' ) Input_Opt%NHMSe
       ErrMsg = 'Input%Opt%NHMSe = ' // TimeStr        // &
                ' is not a valid clock time!'          // &
                ' Please check your "input.geos" file.'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    ! Compute the length of the simulation, in elapsed seconds
    JulianDateStart = GET_JD( Input_Opt%NymdB, Input_Opt%NhmsB )
    JulianDateEnd   = GET_JD( Input_Opt%NymdE, Input_Opt%NhmsE )
    Input_Opt%SimLengthSec = NINT( ( JulianDateEnd - JulianDateStart ) &
                             * 86400_f8)

    !-----------------------------------------------------------------
    ! Data directories
    !-----------------------------------------------------------------

    ! Run directory
    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'RUN_DIR', RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    READ( SUBSTRS(1:N), '(a)' ) Input_Opt%RUN_DIR

    ! Root data dir
    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'DATA_DIR', RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    READ( SUBSTRS(1:N), '(a)' ) Input_Opt%DATA_DIR

    ! Make sure DATA-DIR ends with a "/" character
    C = LEN_TRIM( Input_Opt%DATA_DIR )
    IF ( Input_Opt%DATA_DIR(C:C) /= '/' ) THEN
       Input_Opt%DATA_DIR = TRIM( Input_Opt%DATA_DIR ) // '/'
    ENDIF

    ! Create CHEM_INPUTS directory
    Input_Opt%CHEM_INPUTS_DIR = TRIM( Input_Opt%DATA_DIR ) // &
                                'CHEM_INPUTS/'

    !-----------------------------------------------------------------
    ! Meteorology field
    !-----------------------------------------------------------------
    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'MetField', RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    READ( SUBSTRS(1:N), '(a)' ) Input_Opt%MetField

    ! Make sure a valid met field is specified
    Met = To_UpperCase( TRIM( Input_Opt%MetField ) )
    SELECT CASE( TRIM( Met ) )
    CASE( 'GEOS-FP', 'GEOSFP' )
       Input_Opt%MetField = 'GEOSFP'
    CASE( 'MERRA-2', 'MERRA2' )
       Input_Opt%MetField = 'MERRA2'
    CASE DEFAULT
       ErrMsg = Trim( Input_Opt%MetField) // ' is not a valid '  // &
                ' met field. Supported met fields are GEOS-FP '   // &
                ' and MERRA-2. Please check your "input.geos" file.'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    END SELECT

    !-----------------------------------------------------------------
    ! Simulation type
    !-----------------------------------------------------------------
    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'SIM_TYPE', RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    READ( SUBSTRS(1:N), * ) Input_Opt%SimulationName

    ! Error check simulation name
    Sim = To_UpperCase( TRIM( Input_Opt%SimulationName ) )
    IF ( TRIM(Sim) /= 'ACIDUPTAKE'       .and. &
         TRIM(Sim) /= 'AEROSOL'          .and. &
         TRIM(Sim) /= 'APM'              .and. &
         TRIM(Sim) /= 'BENCHMARK'        .and. &
         TRIM(Sim) /= 'CH4'              .and. &
         TRIM(Sim) /= 'CO2'              .and. &
         TRIM(Sim) /= 'COMPLEXSOA'       .and. &
         TRIM(Sim) /= 'COMPLEXSOA_SVPOA' .and. &
         TRIM(Sim) /= 'HEMCO'            .and. &
         TRIM(Sim) /= 'HG'               .and. &
         TRIM(Sim) /= 'MARINEPOA'        .and. &
         TRIM(Sim) /= 'POPS'             .and. &
         TRIM(Sim) /= 'RRTMG'            .and. &
         TRIM(Sim) /= 'STANDARD'         .and. &
         TRIM(Sim) /= 'TRANSPORTTRACERS' .and. &
         TRIM(Sim) /= 'TROPCHEM'         .and. &
         TRIM(Sim) /= 'TAGCO'            .and. &
         TRIM(Sim) /= 'TAGCH4'           .and. &
         TRIM(Sim) /= 'TAGHG'            .and. &
         TRIM(Sim) /= 'TAGO3'            .and. &
         TRIM(Sim) /= 'TOMAS12'          .and. &
         TRIM(Sim) /= 'TOMAS15'          .and. &
         TRIM(Sim) /= 'TOMAS30'          .and. &
         TRIM(Sim) /= 'TOMAS40'          ) THEN
       ErrMsg = Trim( Input_Opt%SimulationName) // ' is not a'      // &
                ' valid simulation. Supported simulations are:'     // &
                ' AcidUptake, Aerosol, APM, Benchmark, CH4, CO2,'   // &
                ' ComplexSOA, ComplexSOA_SVPOA, HEMCO, Hg,'         // &
                ' MarinePOA, POPs, RRTMG, Standard,'                // &
                ' TransportTracers, Tropchem, TagCO, TagCH4, TagHg,'// &
                ' TagO3, TOMAS12, TOMAS15, TOMAS30, and TOMAS40.'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    ! Set simulation type flags in Input_Opt
    Input_Opt%ITS_A_CH4_SIM      = ( TRIM(Sim) == 'CH4'              .or. &
                                     TRIM(Sim) == 'TAGCH4'           )
    Input_Opt%ITS_A_CO2_SIM      = ( TRIM(Sim) == 'CO2'              )
    Input_Opt%ITS_A_FULLCHEM_SIM = ( TRIM(Sim) == 'ACIDUPTAKE'       .or. &
                                     TRIM(Sim) == 'APM'              .or. &
                                     TRIM(Sim) == 'BENCHMARK'        .or. &
                                     TRIM(Sim) == 'COMPLEXSOA'       .or. &
                                     TRIM(Sim) == 'COMPLEXSOA_SVPOA' .or. &
                                     TRIM(Sim) == 'HEMCO'            .or. &
                                     TRIM(Sim) == 'MARINEPOA'        .or. &
                                     TRIM(Sim) == 'RRTMG'            .or. &
                                     TRIM(Sim) == 'STANDARD'         .or. &
                                     TRIM(Sim) == 'TROPCHEM'         .or. &
                                     TRIM(Sim) == 'TOMAS12'          .or. &
                                     TRIM(Sim) == 'TOMAS15'          .or. &
                                     TRIM(Sim) == 'TOMAS30'          .or. &
                                     TRIM(Sim) == 'TOMAS40'          )
    Input_Opt%ITS_A_MERCURY_SIM  = ( TRIM(Sim) == 'HG'               .or. &
                                     TRIM(Sim) == 'TAGHG'            )
    Input_Opt%ITS_A_POPS_SIM     = ( TRIM(Sim) == 'POPS'             )
    Input_Opt%ITS_A_RnPbBe_SIM   = ( TRIM(Sim) == 'TRANSPORTTRACERS' )
    Input_Opt%ITS_A_TAGO3_SIM    = ( TRIM(Sim) == 'TAGO3'            )
    Input_Opt%ITS_A_TAGCO_SIM    = ( TRIM(Sim) == 'TAGCO'            )
    Input_Opt%ITS_AN_AEROSOL_SIM = ( TRIM(Sim) == 'AEROSOL'          )

    !-----------------------------------------------------------------
    ! Turn on debug output
    !-----------------------------------------------------------------
    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'Debug output', RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    READ( SUBSTRS(1:N), * ) Input_Opt%LPRT

    !-----------------------------------------------------------------
    ! Turn on GEOS-Chem timers
    !-----------------------------------------------------------------
    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'GEOS-Chem timers', RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    READ( SUBSTRS(1:N), * ) Input_Opt%useTimers

    ! Separator line
    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'separator 1', RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    ! Return success
    RC = GC_SUCCESS

    !=================================================================
    ! Print to screen
    !=================================================================
    IF ( Input_Opt%amIRoot ) THEN
       WRITE( 6, '(/,a)' ) 'SIMULATION MENU'
       WRITE( 6, '(  a)' ) '---------------'
       WRITE( 6, 100 ) 'Start time of run           : ', &
                        Input_Opt%NYMDb, Input_Opt%NHMSb
       WRITE( 6, 100 ) 'End time of run             : ', &
                        Input_Opt%NYMDe, Input_Opt%NHMSe
       WRITE( 6, 110 ) 'Run directory               : ', &
                        TRIM( Input_Opt%RUN_DIR )
       WRITE( 6, 110 ) 'Data Directory              : ', &
                        TRIM( Input_Opt%DATA_DIR )
       WRITE( 6, 110 ) 'CHEM_INPUTS directory       : ', &
                        TRIM( Input_Opt%CHEM_INPUTS_DIR )
       WRITE( 6, 110 ) 'Meteorology field           : ', &
                        TRIM( Input_Opt%MetField )
       WRITE( 6, 110 ) 'Simulation name             : ', &
                        TRIM( Input_Opt%SimulationName )
       WRITE( 6, 120 ) 'Turn on debug output        : ', &
                        Input_Opt%LPRT
       WRITE( 6, 120 ) 'Turn on GEOS-Chem timers    : ', &
                        Input_Opt%useTimers
    ENDIF

    ! Format statements
100 FORMAT( A, I8.8, 1X, I6.6 )
110 FORMAT( A, A              )
120 FORMAT( A, L5             )

    !=================================================================
    ! Call setup routines from other GEOS-CHEM modules
    !=================================================================

#if defined( EXTERNAL_GRID ) || defined( EXTERNAL_FORCING )
    !-----------------------------------------------------------------
    !         %%%%%%% GEOS-Chem HP (with ESMF & MPI) %%%%%%%
    !
    ! If we are connecting to the ESMF interface, we need to take
    ! the start & end dates as defined in the ESMF resource file.
    ! (i.e. GEOSCHEMchem_GridComp_mod.rc) instead of those in
    ! input.geos.  This is because the ESMF Clock object needs to be
    ! defined at the highest level (in the driver routine), before
    ! input.geos is ever read.
    !
    ! Therefore, we will assign the start & end date fields (i.e.
    ! Input_Opt%NYMDb, Input_Opt%NYMDe, Input_Opt%NHMSb, and
    ! Input_Opt%NHMSe) in the Gridded Component module file
    ! GEOSCHEMchem_GridComp_Mod.F90 (i.e. two levels higher
    ! in the code).  We don'need to define those fields here, so
    ! we have bracketed this with an #ifdef.
    !-----------------------------------------------------------------

    ! Get the starting UTC time from Input_Opt%NHMSb for use below
    CALL YMD_Extract( Input_Opt%NHMSb, H, M, S )
    init_UTC = ( H + ( M / 60 ) + ( S / 3600 ) )

    ! Pass the values for the start & end times of the simulation directly
    ! to GeosUtil/time_mod.F90 via subroutine ACCEPT_EXTERNAL_DATE_TIME.
    ! (bmy, 12/6/12)
    CALL Accept_External_Date_Time( value_NYMDb = Input_Opt%NYMDb,   &
                                    value_NHMSb = Input_Opt%NHMSb,   &
                                    value_NYMDe = Input_Opt%NYMDe,   &
                                    value_NHMSe = Input_Opt%NHMSe,   &
                                    value_NYMD  = Input_Opt%NYMDb,   &
                                    value_NHMS  = Input_Opt%NHMSb,   &
                                    value_UTC   = init_UTC,          &
                                    RC          = RC                 )
#else
    !-----------------------------------------------------------------
    !         %%%%%%% GEOS-Chem CLASSIC (with OpenMP) %%%%%%%
    !
    ! If we are not using ESMF, then call the traditional GEOS-Chem
    ! timing routines (from GeosUtil/time_mod.F90) to set the start &
    ! end times of the simulation, as well as the current time.
    ! (bmy, 12/6/12)
    !-----------------------------------------------------------------

    ! Set start time of run in "time_mod.F90"
    CALL SET_BEGIN_TIME( Input_Opt%NYMDb, Input_Opt%NHMSb, RC  )
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered in "Set_Begin_Time"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    ! Set end time of run in "time_mod.F90"
    ErrMsg = 'Error encountered in "Set_Begin_Time"!'
    CALL SET_END_TIME( Input_Opt%NYMDe, Input_Opt%NHMSe, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    ! Set the current time
    CALL SET_CURRENT_TIME()
#endif

#ifdef BPCH_DIAG
    ! Set the start of the 1st diagnostic interval
    CALL SET_DIAGb( GET_TAU() )
#endif

    ! Set menu counter
    CT1 = CT1 + 1

  END SUBROUTINE READ_SIMULATION_MENU
!EOC
#if !(defined( EXTERNAL_GRID ) || defined( EXTERNAL_FORCING ))
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_grid_menu
!
! !DESCRIPTION: Subroutine READ\_GRID\_MENU reads the GRID MENU section
!  of the GEOS-Chem input file.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE READ_GRID_MENU( Input_Opt, State_Grid, RC )
!
! !USES:
!
    USE CharPak_Mod,        ONLY : StrSplit
    USE ErrCode_Mod
    USE Input_Opt_Mod,      ONLY : OptInput
    USE State_Grid_Mod,     ONLY : GrdState
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(INOUT) :: Input_Opt   ! Input options
    TYPE(GrdState), INTENT(INOUT) :: State_Grid  ! Grid State
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure
!
! !REVISION HISTORY:
!  20 Oct 2018 - M. Sulprizio- Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER            :: N, C
    INTEGER            :: nSubStrs

    ! Strings
    CHARACTER(LEN=255) :: ErrMsg, ThisLoc, nLev
    CHARACTER(LEN=10)  :: XMin_Str, XMax_Str
    CHARACTER(LEN=10)  :: YMin_Str, YMax_Str

    ! Arrays
    CHARACTER(LEN=255) :: SUBSTRS(MAXDIM)

    !=================================================================
    ! READ_GRID_MENU begins here!
    !=================================================================

    ! Initialize
    RC      = GC_SUCCESS
    ErrMsg  = 'Error reading the "input.geos" file!'
    ThisLoc = ' -> at Read_Grid_Menu (in GeosCore/input_mod.F90)'

    !-----------------------------------------------------------------
    ! Grid resolution
    !-----------------------------------------------------------------
    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'Grid resolution', RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    READ( SUBSTRS(1:N), * ) State_Grid%GridRes

    ! Split into two values, separated by 'x'
    CALL StrSplit( Trim(State_Grid%GridRes) , 'x', SubStrs, nSubStrs )

    ! Stop with error if there are more than two substrings
    IF ( nSubStrs /= 2 ) THEN
       ErrMsg = 'Error in extracting delta X and Y values from'    // &
                ' State_Grid%GridRes. Values must be separated by' // &
                ' an x. Please check your "input.geos" file.'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    ! Save the delta X and Y values
    READ( SubStrs(1), '(f10.4)' ) State_Grid%DY
    READ( SubStrs(2), '(f10.4)' ) State_Grid%DX

    !-----------------------------------------------------------------
    ! Longitude range
    !-----------------------------------------------------------------

    ! Read in min and max longitude
    CALL SPLIT_ONE_LINE( SUBSTRS, N, 2, 'Lon min/max', RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    READ( SUBSTRS(1:N), * ) State_Grid%XMin, State_Grid%XMax

    ! Make sure values are in valid range
    IF ( State_Grid%XMin >= State_Grid%XMax ) THEN
       WRITE( XMin_Str, '(i10)' ) State_Grid%XMin
       WRITE( XMax_Str, '(i10)' ) State_Grid%XMax
       ErrMsg = 'Lower lon must be smaller than upper lon: ' // &
                TRIM( XMin_Str ) // ' ' // TRIM( XMax_Str )
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    ! Compute X grid dimension
    State_Grid%NX = FLOOR( ( State_Grid%XMax - State_Grid%XMin ) / &
                             State_Grid%DX )

    !-----------------------------------------------------------------
    ! Latitude range
    !-----------------------------------------------------------------

    ! Read in min and max latitude
    CALL SPLIT_ONE_LINE( SUBSTRS, N, 2, 'Lat min/max', RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    READ( SUBSTRS(1:N), * ) State_Grid%YMin, State_Grid%YMax

    ! Make sure values are in valid range
    IF ( State_Grid%YMin >= State_Grid%YMax ) THEN
       WRITE( YMin_Str, '(i10)' ) State_Grid%YMin
       WRITE( YMax_Str, '(i10)' ) State_Grid%YMax
       ErrMsg = 'Lower lat must be smaller than upper lat: ' // &
                TRIM( YMin_Str ) // ' ' // TRIM( YMax_Str )
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    ! Restrict latitude values to -90.0 and 90.0
    IF ( State_Grid%YMin < -90.0_fp ) THEN
       WRITE( YMin_Str, '(i10)' ) State_Grid%YMin
       ErrMsg = 'Lower latitude must be between -90 and 90 degN: ' // &
                TRIM( YMin_Str )
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    IF ( State_Grid%YMax > 90.0_fp ) THEN
       WRITE( YMax_Str, '(i10)' ) State_Grid%YMax
       ErrMsg = 'Upper latitude must be between -90 and 90 degN: ' // &
                TRIM( YMax_Str )
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    ! Compute Y grid dimension
    State_Grid%NY = FLOOR( ( State_Grid%YMax - State_Grid%YMin ) / &
                             State_Grid%DY ) + 1

    ! Use half-sized polar boxes?
    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'Half Polar', RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    READ( SUBSTRS(1:N), * ) State_Grid%HalfPolar

    !-----------------------------------------------------------------
    ! Level range
    !-----------------------------------------------------------------

    ! Read in number of levels
    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'NZ', RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    READ( SUBSTRS(1:N), * ) State_Grid%NZ

    !-----------------------------------------------------------------
    ! Nested grid settings
    !-----------------------------------------------------------------

    ! Is it a nested grid?
    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'Nested grid', RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    READ( SUBSTRS(1:N), * ) State_Grid%NestedGrid

    IF ( State_Grid%NestedGrid ) THEN
       ! Increase NX by 1
       State_Grid%NX        = State_Grid%NX + 1

       ! For now hardcode HalfPolar to false when using a nested grid
       State_Grid%HalfPolar = .FALSE.
    ENDIF

    ! Nested grid transport offsets
    CALL SPLIT_ONE_LINE( SUBSTRS, N, 4, 'Buffers N/S/E/W', RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    READ( SUBSTRS(1:N), * ) State_Grid%NorthBuffer, &
                            State_Grid%SouthBuffer, &
                            State_Grid%EastBuffer,  &
                            State_Grid%WestBuffer

    ! Set buffers to zero for global grids
    IF ( .not. State_Grid%NestedGrid ) THEN
       State_Grid%NorthBuffer = 0
       State_Grid%SouthBuffer = 0
       State_Grid%EastBuffer  = 0
       State_Grid%WestBuffer  = 0
    ENDIF

    ! Separator line
    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'separator 1', RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    ! Return success
    RC = GC_SUCCESS

    !=================================================================
    ! Print to screen
    !=================================================================
    IF ( Input_Opt%amIRoot ) THEN
       WRITE( 6, '(/,a)' ) 'GRID MENU'
       WRITE( 6, '(  a)' ) '---------------'
       WRITE( 6, 100 ) 'Grid resolution             : ', &
                        TRIM( State_Grid%GridRes )
       WRITE( 6, 110 ) 'Min/max longitude           : ', &
                        State_Grid%XMin, State_Grid%XMax
       WRITE( 6, 110 ) 'Min/max longitude           : ', &
                        State_Grid%YMin, State_Grid%YMax
       WRITE( 6, 120 ) 'X grid dimension            : ', &
                        State_Grid%NX
       WRITE( 6, 120 ) 'Y grid dimension            : ', &
                        State_Grid%NY
       WRITE( 6, 120 ) 'Z grid dimension            : ', &
                        State_Grid%NZ
       WRITE( 6, 130 ) 'Use half-sized polar boxes? : ', &
                        State_Grid%HalfPolar
       WRITE( 6, 130 ) 'Is this a nested-grid sim?  : ', &
                        State_Grid%NestedGrid
       WRITE( 6, 140 ) ' --> Buffer zone (N S E W ) : ', &
                        State_Grid%NorthBuffer,          &
                        State_Grid%SouthBuffer,          &
                        State_Grid%EastBuffer,           &
                        State_Grid%WestBuffer
    ENDIF

    ! Format statements
100 FORMAT( A, A                )
110 FORMAT( A, F10.4, 1X, F10.4 )
120 FORMAT( A, I5               )
130 FORMAT( A, L5               )
140 FORMAT( A, 4I5              )

  END SUBROUTINE READ_GRID_MENU
!EOC
#endif
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_timestep_menu
!
! !DESCRIPTION: Subroutine READ\_TIMESTEP\_MENU reads the TIMESTEP MENU
!  section of the GEOS-Chem input file.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE READ_TIMESTEP_MENU( Input_Opt, State_Grid, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE Input_Opt_Mod,      ONLY : OptInput
    USE State_Grid_Mod,     ONLY : GrdState
!
! !INPUT PARAMETERS:
!
    TYPE(GrdState), INTENT(IN)    :: State_Grid  ! Grid State object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(INOUT) :: Input_Opt   ! Input options
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure
!
! !REVISION HISTORY:
!  09 Aug 2017 - M. Sulprizio- Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER            :: N

    ! Strings
    CHARACTER(LEN=255) :: ErrMsg, ThisLoc

    ! Arrays
    CHARACTER(LEN=255) :: SUBSTRS(MAXDIM)

    !=================================================================
    ! READ_TIMESTEP_MENU begins here!
    !=================================================================

    ! Initialize
    RC      = GC_SUCCESS
    ErrMsg  = 'Error reading the "input.geos" file!'
    ThisLoc = ' -> at Read_Timestep_Menu (in module GeosCore/input_mod.F90)'

    ! Transport/convection timestep
    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'TS_DYN', RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    READ( SUBSTRS(1:N), * ) Input_Opt%TS_DYN
    Input_Opt%TS_CONV = Input_Opt%TS_DYN

    ! Chemistry/emissions timestep
    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'TS_CHEM', RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    READ( SUBSTRS(1:N), * ) Input_Opt%TS_CHEM
    Input_Opt%TS_EMIS = Input_Opt%TS_CHEM

    ! Separator line
    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'separator 1', RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    !=================================================================
    ! Error checks
    !=================================================================

    IF ( Input_Opt%SimLengthSec < Input_Opt%TS_DYN .or. &
         Input_Opt%SimLengthSec < Input_Opt%TS_CHEM ) THEN
       IF ( Input_Opt%amIRoot ) THEN

          WRITE( 6, *   ) ''
          WRITE( 6, *   ) 'The length of the simulation is shorter '
          WRITE( 6, *   ) 'than the transport and/or chemistry '
          WRITE( 6, *   ) 'timesteps. Check the settings in '
          WRITE( 6, *   ) 'the "input.geos" file.'
          WRITE( 6, *   ) ''
          WRITE( 6, 100 ) 'Transport/Convection [sec]  : ', &
                           Input_Opt%TS_DYN
          WRITE( 6, 100 ) 'Chemistry/Emissions  [sec]  : ', &
                           Input_Opt%TS_CHEM
          WRITE( 6, 100 ) 'Simulation duration  [sec]  : ', &
                           Input_Opt%SimLengthSec
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN

       ENDIF
    ENDIF

    IF ( TRIM(Input_Opt%MetField) == 'MERRA2' .and. &
         TRIM(State_Grid%GridRes) == '0.5x0.625' ) THEN
       IF ( Input_Opt%ITS_A_CH4_SIM .or. Input_Opt%ITS_A_CO2_SIM ) THEN
          IF ( Input_Opt%TS_DYN > 300 .or. Input_Opt%TS_CHEM > 600 ) THEN
             IF ( Input_Opt%amIRoot ) THEN

                WRITE( 6,* ) ''
                WRITE( 6,* ) 'It has been noted that MERRA-2 nested grid'
                WRITE( 6,* ) ' simulations can have very high species'
                WRITE( 6,* ) ' concentrations in the stratosphere caused'
                WRITE( 6,* ) ' by a violation of the CFL condition due to'
                WRITE( 6,* ) ' strong stratospheric winds. This is'
                WRITE( 6,* ) ' especially problematic when using total'
                WRITE( 6,* ) ' column concentrations. To avoid the issue,'
                WRITE( 6,* ) ' a timestep of 5/10 instead of 10/20 is'
                WRITE( 6,* ) ' recommended for CH4 and CO2 simulations.'
                WRITE( 6,* ) ''
                WRITE( 6,* ) 'You may remove this trap at your own peril,'
                WRITE( 6,* ) ' by commenting out the call to GC_ERROR in'
                WRITE( 6,* ) ' GeosCore/input_mod.F90. '
                WRITE( 6,* ) ''
                WRITE( 6,* ) 'See the MERRA-2 implementation details page'
                WRITE( 6,* ) ' on the GEOS-Chem wiki for details'
                CALL GC_Error( ErrMsg, RC, ThisLoc )
                RETURN

             ENDIF
          ENDIF
       ENDIF
    ENDIF

    ! Return success
    RC = GC_SUCCESS

    !=================================================================
    ! Print to screen
    !=================================================================
    IF ( Input_Opt%amIRoot ) THEN
       WRITE( 6, '(/,a)' ) 'TIMESTEP MENU'
       WRITE( 6, '(  a)' ) '---------------'
       WRITE( 6, 100     ) 'Transport/Convection [sec]  : ', &
                            Input_Opt%TS_DYN
       WRITE( 6, 100     ) 'Chemistry/Emissions  [sec]  : ', &
                            Input_Opt%TS_CHEM
    ENDIF

    ! Format statements
100 FORMAT( A, I5  )

  END SUBROUTINE READ_TIMESTEP_MENU
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_advected_species_menu
!
! !DESCRIPTION: Subroutine READ\_ADVECTED\_SPECIES\_MENU reads the ADVECTED
!  SPECIES MENU section of the GEOS-Chem input file.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE READ_ADVECTED_SPECIES_MENU( Input_Opt, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE Input_Opt_Mod, ONLY : OptInput, Set_Input_Opt_Advect
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(INOUT) :: Input_Opt   ! Input options
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure
!
! !REVISION HISTORY:
!  20 Jul 2004 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER            :: N, T

    ! Strings
    CHARACTER(LEN=255) :: ErrMsg, ThisLoc

    ! Arrays
    CHARACTER(LEN=255) :: SUBSTRS(MAXDIM)

    !=================================================================
    ! READ_ADVECTED_SPECIES_MENU begins here!
    !
    ! Get the simulation type
    !=================================================================

    ! Initialize
    RC      = GC_SUCCESS
    ErrMsg  = 'Error reading the "input.geos" file!'
    ThisLoc = ' -> at Read_Advected_Species_Menu (in GeosCore/input_mod.F90)'

    !=================================================================
    ! Read advected species name
    !=================================================================

    ! Initialize counter
    T = 0

    ! Do until exit
    DO

       ! Read species name
       CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'species', RC )
       IF ( RC /= GC_SUCCESS ) THEN
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       ! Exit here if separator line is encountered
       IF ( INDEX ( TRIM(SUBSTRS(1)), '-----' ) > 0 ) EXIT

       ! Increase number of advected species by one
       T = T + 1

       ! Stop simulation and print warning if we exceed maximum number
       ! of advected species hardcoded in input_opt_mod.F90
       IF ( T > Input_Opt%Max_AdvectSpc ) THEN
          ErrMsg = 'Number of advected species exceeds maximum. ' // &
                   'This value can be modified in input_opt_mod.F90.'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       ! Save advected species name
       Input_Opt%AdvectSpc_Name(T)  = TRIM( SUBSTRS(1) )

    ENDDO

    ! Total number of advected species
    Input_Opt%N_ADVECT = T

    !=================================================================
    ! Call setup routines from other F90 modules
    !=================================================================

    ! Split into tagged species (turn off for full-chemistry)
    IF ( Input_Opt%ITS_A_FULLCHEM_SIM ) THEN

       ! There are no tagged species for fullchem
       Input_Opt%LSPLIT = .FALSE.

    ELSE IF ( Input_Opt%ITS_A_MERCURY_SIM ) THEN

       ! Need Hg0, Hg2, HgP for tagged Mercury
       Input_Opt%LSPLIT = ( Input_Opt%N_ADVECT > 3 )

    ELSE

       ! All other simulations: tagged if more than 1 species listed
       Input_Opt%LSPLIT = ( Input_Opt%N_ADVECT > 1 )

    ENDIF

    ! Initialize arrays in Input_Opt that depend on N_ADVECT
    CALL Set_Input_Opt_Advect( Input_Opt, RC )

    ! Set menu counter
    CT1 = CT1 + 1

  END SUBROUTINE READ_ADVECTED_SPECIES_MENU
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_aerosol_menu
!
! !DESCRIPTION: Subroutine READ\_AEROSOL\_MENU reads the AEROSOL MENU
!  section of the GEOS-Chem input file.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE READ_AEROSOL_MENU( Input_Opt, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE Input_Opt_Mod, ONLY : OptInput
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(INOUT) :: Input_Opt   ! Input options
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure
!
! !REMARKS:
!  Move error checks that depend on species indices to the subroutine
!  DO_ERROR_CHECKS.  This is now called from GC_INIT_EXTRA, after the
!  initialization of the species database.
!
! !REVISION HISTORY:
!  20 Jul 2004 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER            :: N, T

    ! Strings
    CHARACTER(LEN=255) :: ErrMsg, ThisLoc

    ! Arrays
    CHARACTER(LEN=255) :: SUBSTRS(MAXDIM)

    !=================================================================
    ! READ_AEROSOL_MENU begins here!
    !=================================================================

    ! Initialize
    RC      = GC_SUCCESS
    ErrMsg  = 'Error reading the "input.geos" file!'
    ThisLoc = ' -> at Read_Aerosol_Menu (in module GeosCore/input_mod.F90)'

    ! Error check
    IF ( CT1 /= 2 ) THEN
       ErrMsg = 'SIMULATION MENU & ADVECTED SPECIES MENU ' // &
                'must be read in first!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    ! Use online sulfate aerosols?
    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'LSULF', RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    READ( SUBSTRS(1:N), * ) Input_Opt%LSULF

    ! Use metal catalyzed oxidation of SO2?
    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'LMETALCATSO2', RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    READ( SUBSTRS(1:N), * ) Input_Opt%LMETALCATSO2

    ! Use online carbon aerosols?
    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'LCARB', RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    READ( SUBSTRS(1:N), * ) Input_Opt%LCARB

    ! Use brown carbon aerosols?
    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'LBRC', RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    READ( SUBSTRS(1:N), * ) Input_Opt%LBRC

    ! Use secondary organic aerosols?
    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'LSOA', RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    READ( SUBSTRS(1:N), * ) Input_Opt%LSOA

    ! SOAupdate: Add Semi-volatile POA switch (hotp 8/9/09)
    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'LSVPOA', RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    READ( SUBSTRS(1:N), * ) Input_Opt%LSVPOA

    ! Use online dust aerosols ?
    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'LDUST', RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    READ( SUBSTRS(1:N), * ) Input_Opt%LDUST

    ! Use SO2 and HNO3 uptake on dust aerosols
    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'LDSTUP', RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    READ( SUBSTRS(1:N), * ) Input_Opt%LDSTUP

    ! Use online sea-salt aerosols?
    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'LSSALT', RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    READ( SUBSTRS(1:N), * ) Input_Opt%LSSALT

    ! Accum mode seasalt radii bin edges [um]
    CALL SPLIT_ONE_LINE( SUBSTRS, N, 2, 'SALA_REDGE_um', RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    DO T = 1, N
       READ( SUBSTRS(T), * ) Input_Opt%SALA_REDGE_um(T)
    ENDDO

    ! Coarse mode seasalt radii bin edges [um]
    CALL SPLIT_ONE_LINE( SUBSTRS, N, 2, 'SALC_REDGE_um', RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    DO T = 1, N
       READ( SUBSTRS(T), * ) Input_Opt%SALC_REDGE_um(T)
    ENDDO

    ! Use marine organic aerosols?
    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'LMPOA', RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    READ( SUBSTRS(1:N), * ) Input_Opt%LMPOA

    ! Apply gravitational settling in stratosphere?
    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'LGRAVSTRAT', RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    READ( SUBSTRS(1:N), * ) Input_Opt%LGRAVSTRAT

    ! Use solid polar stratospheric clouds (PSCs)?
    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'LSOLIDPSC', RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    READ( SUBSTRS(1:N), * ) Input_Opt%LSOLIDPSC

    ! Allow homogeneous nucleation of NAT?
    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'LHOMNUCNAT', RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    READ( SUBSTRS(1:N), * ) Input_Opt%LHOMNUCNAT

    ! NAT supercooling requirement (K)
    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'T_NAT_SUPERCOOL', RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    READ( SUBSTRS(1:N), * ) Input_Opt%T_NAT_SUPERCOOL

    ! Ice supersaturation ratio requirement
    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'P_ICE_SUPERSAT', RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    READ( SUBSTRS(1:N), * ) Input_Opt%P_ICE_SUPERSAT

    ! Perform PSC-related heterogeneous chemistry?
    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'LPSCCHEM', RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    READ( SUBSTRS(1:N), * ) Input_Opt%LPSCCHEM

    ! Include stratospheric aerosols optical depths?
    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'LSTRATOD', RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    READ( SUBSTRS(1:N), * ) Input_Opt%LSTRATOD

    ! Include BC absorption enhancement due to coating? (xnw, 8/24/15)
    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'LBCAE', RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    READ( SUBSTRS(1:N), * ) Input_Opt%LBCAE

    ! Define BC absorption enhancement (xnw, 8/24/15)
    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'BCAE_1', RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    READ( SUBSTRS(1:N), * ) Input_Opt%BCAE_1

    ! Define BC absorption enhancement (xnw, 8/24/15)
    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'BCAE_2', RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    READ( SUBSTRS(1:N), * ) Input_Opt%BCAE_2

    ! Photoylse nitrate aerosol? (TMS, 23/08/18)
    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'hvAerNIT', RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    READ( SUBSTRS(1:N), * ) Input_Opt%hvAerNIT

    ! scalar for JHNO3 for photoylsing NITs aerosol (TMS, 23/08/18)
    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'hvAerNIT_JNITs', RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    READ( SUBSTRS(1:N), * ) Input_Opt%hvAerNIT_JNITs

    ! scalar for JHNO3 for photoylsing NIT aerosol (TMS, 23/08/18)
    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'hvAerNIT_JNIT', RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    READ( SUBSTRS(1:N), * ) Input_Opt%hvAerNIT_JNIT

    ! Fraction for JNITS/NIT channel A (HNO2) for NITs photoylsis (TMS, 10/10/18)
    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'JNITChanA', RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    READ( SUBSTRS(1:N), * ) Input_Opt%JNITChanA

    ! Fraction for JNITs/NIT channel B (NO2) for NITs photoylsis (TMS, 10/10/18)
    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'JNITChanB', RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    READ( SUBSTRS(1:N), * ) Input_Opt%JNITChanB

    ! Separator line
    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'separator 2', RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    !=================================================================
    ! Error checks
    !=================================================================

    ! Make sure that SALA, SALC bins are contiguous
    IF ( Input_Opt%SALA_REDGE_um(2) /= &
         Input_Opt%SALC_REDGE_um(1)     ) THEN
       ErrMsg = 'SALA and SALC bin edges are not contiguous!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    ! Turn off switches for simulations that don't use aerosols
    IF ( ( .not. Input_Opt%ITS_A_FULLCHEM_SIM )  .and. &
         ( .not. Input_OPt%ITS_AN_AEROSOL_SIM ) ) THEN
       Input_Opt%LSULF        = .FALSE.
       Input_Opt%LMETALCATSO2 = .FALSE.
       Input_Opt%LCARB        = .FALSE.
       Input_Opt%LBRC         = .FALSE.
       Input_Opt%LSOA         = .FALSE.
       Input_Opt%LDUST        = .FALSE.
       Input_Opt%LSSALT       = .FALSE.
       Input_Opt%LMPOA        = .FALSE.
       Input_Opt%LSVPOA       = .FALSE.
       Input_Opt%LBCAE        = .FALSE.
       Input_Opt%hvAerNIT     = .FALSE.
    ENDIF

    ! Return success
    RC = GC_SUCCESS

    !=================================================================
    ! Print to screen
    !=================================================================
    IF ( Input_Opt%amIRoot ) THEN
       WRITE( 6, '(/,a)' ) 'AEROSOL MENU'
       WRITE( 6, '(  a)' ) '------------'
       WRITE( 6, 100     ) 'Online SULFATE AEROSOLS?    : ', &
                            Input_Opt%LSULF
       WRITE( 6, 100     ) 'Metal catalyzed SO2 ox.?    : ', &
                            Input_Opt%LMETALCATSO2
       WRITE( 6, 100     ) 'Online CARBON AEROSOLS?     : ', &
                            Input_Opt%LCARB
       WRITE( 6, 100     ) 'Brown Carbon Aerosol?       : ', &
                            Input_Opt%LBRC
       WRITE( 6, 100     ) 'Online COMPLEX SOA?         : ', &
                            Input_Opt%LSOA
       WRITE( 6, 100     ) 'Semivolatile POA?           : ', &
                            Input_Opt%LSVPOA
       WRITE( 6, 100     ) 'Online DUST AEROSOLS?       : ', &
                            Input_Opt%LDUST
       WRITE( 6, 100     ) 'Acid uptake on dust?        : ', &
                            Input_Opt%LDSTUP
       WRITE( 6, 100     ) 'Online SEA SALT AEROSOLS?   : ', &
                            Input_Opt%LSSALT
       WRITE( 6, 110     ) 'Accum  SEA SALT radii [um]  : ', &
                            Input_Opt%SALA_REDGE_um(1),      &
                            Input_Opt%SALA_REDGE_um(2)
       WRITE( 6, 110     ) 'Coarse SEA SALT radii [um]  : ', &
                            Input_Opt%SALC_REDGE_um(1),      &
                            Input_Opt%SALC_REDGE_um(2)
       WRITE( 6, 100     ) 'MARINE ORGANIC AEROSOLS?    : ', &
                            Input_Opt%LMPOA
       WRITE( 6, 100     ) 'Settle strat. aerosols?     : ', &
                            Input_Opt%LGRAVSTRAT
       WRITE( 6, 100     ) 'Online SOLID PSC aerosols?  : ', &
                            Input_Opt%LSOLIDPSC
       WRITE( 6, 100     ) 'Allow hom. NAT nucleation?  : ', &
                            Input_Opt%LHOMNUCNAT
       WRITE( 6, 120     ) 'NAT supercooling requirement: ', &
                            Input_Opt%T_NAT_SUPERCOOL
       WRITE( 6, 120     ) 'Ice supersaturation req.    : ', &
                            ((Input_Opt%P_ICE_SUPERSAT-1)*1.e+2_fp)
       WRITE( 6, 100     ) 'Perform PSC het. chemistry? : ', &
                            Input_Opt%LPSCCHEM
       WRITE( 6, 100     ) 'Use strat. aerosol OD?      : ', &
                            Input_Opt%LSTRATOD
       WRITE( 6, 100     ) 'BC Absorption Enhancement?  : ', &
                            Input_Opt%LBCAE
       WRITE( 6, 105     ) 'Hydrophilic BC AE factor    : ', &
                            Input_Opt%BCAE_1
       WRITE( 6, 105     ) 'Hydrophobic BC AE factor    : ', &
                            Input_Opt%BCAE_2
       WRITE( 6, 100     ) 'Photolyse nitrate aerosol?  : ', &
                            Input_Opt%hvAerNIT
       WRITE( 6, 105     ) 'JNITs scaling of JHNO3      : ', &
                            Input_Opt%hvAerNIT_JNITs
       WRITE( 6, 105     ) 'JNIT scaling of JHNO3       : ', &
                            Input_Opt%hvAerNIT_JNIT
       WRITE( 6, 105     ) 'JNIT(s) channel A (HONO)    : ', &
                            Input_Opt%JNITChanA
       WRITE( 6, 105     ) 'JNIT(s) channel B (NO2)     : ', &
                            Input_Opt%JNITChanB
    ENDIF

100 FORMAT( A, L5     )
105 FORMAT( A, f6.2 )
110 FORMAT( A, f6.2, ' - ', f6.2 )
120 FORMAT( A, f6.2, 'K' )

  END SUBROUTINE READ_AEROSOL_MENU
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_emissions_menu
!
! !DESCRIPTION: Subroutine READ\_EMISSIONS\_MENU reads the EMISSIONS MENU
!  section of the GEOS-Chem input file.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE READ_EMISSIONS_MENU( Input_Opt, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE Input_Opt_Mod, ONLY : OptInput
    USE TIME_MOD,      ONLY : GET_YEAR
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(INOUT) :: Input_Opt   ! Input options
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure
!
! !REMARKS:
!  The Ind_() function now defines all species ID's.  It returns -1 if
!  a species cannot be found.  Therefore now test for Ind_() > 0  for a
!  valid species.
!
! !REVISION HISTORY:
!  20 Jul 2004 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER            :: N

    ! Strings
    CHARACTER(LEN=255) :: ErrMsg, ThisLoc

    ! Arrays
    CHARACTER(LEN=255) :: SUBSTRS(MAXDIM)

    !=================================================================
    ! READ_EMISSIONS_MENU begins here!
    !=================================================================

    ! Initialize
    RC      = GC_SUCCESS
    ErrMsg  = 'Error reading the "input.geos" file!'
    ThisLoc = ' -> at Read_Emissions_Menu (in module GeosCore/input_mod.F90)'

    ! Error check
    IF ( CT1 /= 2 ) THEN
       ErrMsg = 'SIMULATION MENU & ADVECTED SPECIES MENU ' // &
                'must be read in first!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    ! Turn on emissions?
    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'LLEMIS', RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    READ( SUBSTRS(1:N), * ) Input_Opt%LEMIS

    ! HEMCO Input file
    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'HcoConfigFile', RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    READ( SUBSTRS(1:N), '(a)' ) Input_Opt%HcoConfigFile

    ! Separator line (start of UCX options)
    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'separator 1', RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    ! Use variable methane emissions?
    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'LCH4EMIS', RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    READ( SUBSTRS(1:N), * ) Input_Opt%LCH4EMIS

    ! Initialize strat H2O to GEOS-Chem baseline?
    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'LSETH2O', RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    READ( SUBSTRS(1:N), * ) Input_Opt%LSETH2O

    ! Separator line
    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'separator 4', RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    !=================================================================
    ! Error check logical flags
    !=================================================================

    ! Turn off full-chem only switches
    IF ( .not. Input_Opt%ITS_A_FULLCHEM_SIM ) THEN
       Input_Opt%LSETH2O = .FALSE.
    ENDIF

    ! Return success
    RC = GC_SUCCESS

    !=================================================================
    ! Print to screen
    !=================================================================
    IF ( Input_Opt%amIRoot ) THEN
       WRITE( 6, '(/,a)' ) 'EMISSIONS MENU'
       WRITE( 6, '(  a)' ) '--------------'
       WRITE( 6, 100 ) 'Turn on emissions?          : ', &
                       Input_Opt%LEMIS
       WRITE( 6, 130 ) 'HEMCO Configuration file    : ', &
                        TRIM( Input_Opt%HcoConfigFile )
       WRITE( 6, 100 ) 'Use CH4 emissions inventory?: ', &
                        Input_Opt%LCH4EMIS
       WRITE( 6, 100 ) 'Set initial strat H2O?      : ', &
                        Input_Opt%LSETH2O
    ENDIF

    ! FORMAT statements
100 FORMAT( A, L5 )
110 FORMAT( A, I5 )
130 FORMAT( A, A  )

  END SUBROUTINE READ_EMISSIONS_MENU
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_co_sim_menu
!
! !DESCRIPTION: Subroutine READ\_CO\_SIM\_MENU reads the CO SIM MENU
!  section of the GEOS-Chem input file.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE READ_CO_SIM_MENU( Input_Opt, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE Input_Opt_Mod, ONLY : OptInput
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(INOUT) :: Input_Opt   ! Input options
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure
!
! !REVISION HISTORY:
!  24 Mar 2017 - J. Fisher - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER            :: N
    CHARACTER(LEN=255) :: MSG, LOC

    ! Arrays
    CHARACTER(LEN=255) :: SUBSTRS(MAXDIM)

    !=================================================================
    ! READ_CO_SIM_MENU begins here!
    !=================================================================

    ! Assume success
    RC = GC_SUCCESS

    ! Location for error messages
    LOC = ' -> at READ_CO_SIM_MENU (in GeosCore/input_mod.F90)'

    ! Error check
    IF ( CT1 /= 2 ) THEN
       MSG = 'SIMULATION MENU & ADVECTED SPECIES MENU ' // &
             'must be read in first!'
       CALL GC_Error( Msg, RC, Loc )
       RETURN
    ENDIF

    ! Use P(CO) from CH4 from full chemistry simulation?
    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'LPCO_CH4', RC )
    READ( SUBSTRS(1:N), * ) Input_Opt%LPCO_CH4

    ! Use P(CO) from CH4 from full chemistry simulation?
    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'LPCO_NMVOC', RC )
    READ( SUBSTRS(1:N), * ) Input_Opt%LPCO_NMVOC

    !=================================================================
    ! Print to screen
    !=================================================================
    IF ( Input_Opt%ITS_A_TAGCO_SIM .and. Input_Opt%amIRoot ) THEN
       WRITE(6,'(/,a)') 'CO SIMULATION MENU ' // &
             '(overwrites any other settings related to CO)'
       WRITE(6,'(  a)') '-------------------------------------'
       WRITE(6,100    ) 'Use full chem. P(CO) from CH4?:', &
                        Input_Opt%LPCO_CH4
       WRITE(6,100    ) '                   from NMVOC?:', &
                        Input_Opt%LPCO_NMVOC
       WRITE(6,'(  a)') '-------------------------------------'
    ENDIF

    ! FORMAT statements
90  FORMAT( A )
100 FORMAT( A, L5 )
110 FORMAT( A, L5, A )

  END SUBROUTINE READ_CO_SIM_MENU
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_co2_sim_menu
!
! !DESCRIPTION: Subroutine READ\_CO2\_SIM\_MENU reads the CO2 SIM MENU
!  section of the GEOS-Chem input file.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE READ_CO2_SIM_MENU( Input_Opt, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE Input_Opt_Mod,      ONLY : OptInput
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(INOUT) :: Input_Opt   ! Input options
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure
!
! !REVISION HISTORY:
!  02 Mar 2009 - R. Nassar   - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER            :: N

    ! Strings
    CHARACTER(LEN=255) :: ErrMsg, ThisLoc

    ! Arrays
    CHARACTER(LEN=255) :: SUBSTRS(MAXDIM)

    !=================================================================
    ! READ_CO2_SIM_MENU begins here!
    !=================================================================

    ! Initialize
    RC      = GC_SUCCESS
    ErrMsg  = 'Error reading the "input.geos" file!'
    ThisLoc = ' -> at Read_CO2_SIM_Menu (in module GeosCore/input_mod.F90)'

    ! Error check
    IF ( CT1 /= 2 ) THEN
       ErrMsg = 'SIMULATION MENU & ADVECTED SPECIES MENU ' // &
                'must be read in first!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    ! Use Fossil Fuel emissions?
    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'LFOSSIL', RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    READ( SUBSTRS(1:N), * ) Input_Opt%LFOSSIL

    ! Use Ocean Exchange?
    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'LOCEAN', RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    READ( SUBSTRS(1:N), * ) Input_Opt%LOCEAN

    ! Turn on (balanced) biosphere with diurnal cycle?
    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'LBIODIURNAL', RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    READ( SUBSTRS(1:N), * ) Input_Opt%LBIODIURNAL

    ! Use Net Terrestrial Exchange Climatology?
    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'LBIONETCLIM', RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    READ( SUBSTRS(1:N), * ) Input_Opt%LBIONETCLIM

    ! Turn on Ship emissions?
    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'LSHIP', RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    READ( SUBSTRS(1:N), * ) Input_Opt%LSHIP

    ! Turn on Aviation emissions?
    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'LPLANE', RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    READ( SUBSTRS(1:N), * ) Input_Opt%LPLANE

    ! Turn on CO2 3D chemical source and surface correction?
    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'LCHEMCO2', RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    READ( SUBSTRS(1:N), * ) Input_Opt%LCHEMCO2

    ! Tagged CO2 Separator line
    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'separator 1', RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    ! Background CO2 (no emissions or exchange) for Tagged-CO2 runs
    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'LFFBKGRD',RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    READ( SUBSTRS(1:N), * ) Input_Opt%LFFBKGRD

    ! Turn on biosphere and ocean exchange region tagged species?
    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'LBIOSPHTAG', RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    READ( SUBSTRS(1:N), * ) Input_Opt%LBIOSPHTAG

    ! Turn on fossil fuel emission region tagged species?
    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'LFOSSILTAG', RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    READ( SUBSTRS(1:N), * ) Input_Opt%LFOSSILTAG

    ! Turn on global ship emissions tagged species?
    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'LSHIPTAG', RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    READ( SUBSTRS(1:N), * ) Input_Opt%LSHIPTAG

    ! Turn on global aircraft emissions tagged species?
    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'LPLANETAG', RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    READ( SUBSTRS(1:N), * ) Input_Opt%LPLANETAG

    !=================================================================
    ! Print to screen
    !=================================================================
    IF ( Input_Opt%ITS_A_CO2_SIM .and. Input_Opt%amIRoot ) THEN
       WRITE(6,'(/,a)') 'CO2 SIMULATION MENU ' // &
               '(overwrites any other settings related to CO2)'
       WRITE(6,'(  a)') '-------------------------------------'
       WRITE(6,100    ) 'National Fossil Fuel Emission :', &
                         Input_Opt%LFOSSIL
       WRITE(6,100    ) 'Ocean CO2 Uptake/Emission     :', &
                         Input_Opt%LOCEAN
       WRITE(6,100    ) 'Biosphere seas/diurnal cycle  :', &
                         Input_Opt%LBIODIURNAL
       WRITE(6,100    ) 'Net Terr Exch - Climatology   :', &
                         Input_Opt%LBIONETCLIM
       WRITE(6,100    ) 'Intl/Domestic Ship emissions  :', &
                         Input_Opt%LSHIP
       WRITE(6,100    ) 'Intl/Domestic Aviation emiss  :', &
                         Input_Opt%LPLANE
       WRITE(6,100    ) 'CO2 from oxidation (CO,CH4,..):', &
                         Input_Opt%LCHEMCO2
       WRITE(6, 90    ) 'Tagged CO2 settings'
       WRITE(6,100    ) '  Save Fossil CO2 in Bckgrnd  :', &
                         Input_Opt%LFFBKGRD
       WRITE(6,100    ) '  Tag Biosphere/Ocean CO2     :', &
                           Input_Opt%LBIOSPHTAG
       WRITE(6,100    ) '  Tag Fossil Fuel CO2         :', &
                           Input_Opt%LFOSSILTAG
       WRITE(6,100    ) '  Tag Global Ship CO2         :', &
                           Input_Opt%LSHIPTAG
       WRITE(6,100    ) '  Tag Global Aviation CO2     :', &
                           Input_Opt%LPLANETAG
       WRITE(6,'(  a)') '-------------------------------------'
    ENDIF

    ! FORMAT statements
90  FORMAT( A )
100 FORMAT( A, L5 )

  END SUBROUTINE READ_CO2_SIM_MENU
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_chemistry_menu
!
! !DESCRIPTION: Subroutine READ\_CHEMISTRY\_MENU reads the CHEMISTRY MENU
!  section of the GEOS-Chem input file.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE READ_CHEMISTRY_MENU( Input_Opt, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE Input_Opt_Mod, ONLY : OptInput
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(INOUT) :: Input_Opt   ! Input options
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure
!
! !REVISION HISTORY:
!  20 Jul 2004 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER            :: N

    ! Strings
    CHARACTER(LEN=255) :: ErrMsg, ThisLoc

    ! Arrays
    CHARACTER(LEN=255) :: SUBSTRS(MAXDIM)

    !=================================================================
    ! READ_CHEMISTRY_MENU begins here!
    !=================================================================

    ! Initialize
    RC      = GC_SUCCESS
    ErrMsg  = 'Error reading the "input.geos" file!'
    ThisLoc = ' -> at Read_Chemistry_Menu (in module GeosCore/input_mod.F90)'

    ! Error check
    IF ( CT1 /= 2 ) THEN
       ErrMsg = 'SIMULATION MENU & ADVECTED SPECIES MENU ' // &
                'must be read in first!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    ! Turn on chemistry?
    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'LCHEM', RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    READ( SUBSTRS(1:N), * ) Input_Opt%LCHEM

    ! Turn on stratospheric chemistry?
    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'LSCHEM', RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    READ( SUBSTRS(1:N), * ) Input_Opt%LSCHEM

    ! Use Linoz for stratospheric ozone? (Otherwise, Synoz is used)
    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'LLINOZ', RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    READ( SUBSTRS(1:N), * ) Input_Opt%LLINOZ

    ! Use Synoz if Linoz is turned off
    IF ( .not. Input_Opt%LLINOZ ) Input_Opt%LSYNOZ = .TRUE.

    ! Turn on unified strat-trop chemistry?
    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'LUCX', RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    READ( SUBSTRS(1:N), * ) Input_Opt%LUCX

    ! Turn on online stratospheric H2O?
    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'LACTIVEH2O', RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    READ( SUBSTRS(1:N), * ) Input_Opt%LACTIVEH2O

    ! Separator line
    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'separator 1', RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    ! Use online ozone in extinction calculations for FAST-JX?
    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'USE_ONLINE_O3', RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    READ( SUBSTRS(1:N), * ) Input_Opt%USE_ONLINE_O3

    ! Use ozone columns from met fields?
    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'USE_O3_FROM_MET', RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    READ( SUBSTRS(1:N), * ) Input_Opt%USE_O3_FROM_MET

    ! Use ozone columns from TOMS?
    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'USE_TOMS_O3', RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    READ( SUBSTRS(1:N), * ) Input_Opt%USE_TOMS_O3

    ! GAMMA HO2 ?
    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'GAMMA_HO2', RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    READ( SUBSTRS(1:N), * ) Input_Opt%GAMMA_HO2

    ! Separator line
    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'separator 2', RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    !=================================================================
    ! Error check settings
    !=================================================================

#ifndef MODEL_GEOS
    ! Cannot use Synoz with linearized mesospheric chemistry
    IF ( Input_Opt%LUCX .and. Input_Opt%LSCHEM ) THEN
       IF (.not.Input_Opt%LLINOZ) THEN
          ErrMsg = 'Cannot use Synoz with linearized meso. chem.!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
    ENDIF
#endif

    ! Cannot have active H2O without stratospheric chemistry
    IF ( (.not.Input_Opt%LUCX) .and. Input_Opt%LACTIVEH2O ) THEN
       ErrMsg = 'Cannot have active H2O without full strat chem!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    ! FAST-JX is only used for fullchem and offline aerosol
    IF ( Input_Opt%ITS_A_FULLCHEM_SIM  .or. &
         Input_Opt%ITS_AN_AEROSOL_SIM  ) THEN

       ! Make sure either O3 from met or TOMS is selected
       IF ( .not. Input_Opt%USE_O3_FROM_MET .and. &
            .not. Input_Opt%USE_TOMS_O3 ) THEN
          ErrMsg = 'Must select either O3 from met or TOMS/SBUV O3' &
                // 'for O3 values above the chemistry grid!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
       IF ( Input_Opt%USE_O3_FROM_MET .and. &
            Input_Opt%USE_TOMS_O3 ) THEN
          ErrMsg = 'Must select either O3 from met or TOMS/SBUV O3' &
                // 'for O3 values above the chemistry grid!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       ! Make sure specialty simulations select O3 from met or TOMS
       IF ( Input_Opt%ITS_AN_AEROSOL_SIM ) THEN
          IF ( Input_Opt%USE_ONLINE_O3 ) THEN
             ErrMsg= 'Cannot use online O3 for specialty simulations! ' &
                  // 'Select O3 from met or TOMS O3 instead.'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF
       ENDIF

    ELSE

       Input_Opt%USE_ONLINE_O3   = .FALSE.
       Input_Opt%USE_O3_FROM_MET = .FALSE.
       Input_Opt%USE_TOMS_O3     = .FALSE.

    ENDIF

    ! Return success
    RC = GC_SUCCESS

    !=================================================================
    ! Print to screen
    !=================================================================
    IF ( Input_Opt%amIRoot ) THEN
       WRITE( 6, '(/,a)' ) 'CHEMISTRY MENU'
       WRITE( 6, '(  a)' ) '--------------'
       WRITE( 6, 100     ) 'Turn on chemistry?          : ', &
                            Input_Opt%LCHEM
       WRITE( 6, 100     ) 'Use linear. strat. chem?    : ', &
                            Input_Opt%LSCHEM
       WRITE( 6, 100     ) ' => Use Linoz for O3?       : ', &
                            Input_Opt%LLINOZ
       WRITE( 6, 100     ) 'Enable UCX?                 : ', &
                            Input_Opt%LUCX
       WRITE( 6, 100     ) 'Online strat. H2O?          : ', &
                            Input_Opt%LACTIVEH2O
       WRITE( 6, 100     ) 'Online ozone for FAST-JX?   : ', &
                            Input_Opt%USE_ONLINE_O3
       WRITE( 6, 100     ) 'Ozone from met for FAST-JX? : ', &
                            Input_Opt%USE_O3_FROM_MET
       WRITE( 6, 100     ) 'TOMS/SBUV ozone for FAST-JX?: ', &
                            Input_Opt%USE_TOMS_O3
       WRITE( 6, 110     ) 'GAMMA HO2                   : ', &
                            Input_Opt%GAMMA_HO2
       IF ( Input_Opt%USE_ONLINE_O3 ) THEN
          WRITE( 6, '(a)' ) ''
          WRITE( 6, '(a)' ) 'NOTE ABOUT OVERHEAD O3 FOR FAST-JX:'
          WRITE( 6, '(a)' ) ' Online O3 from GEOS-Chem will be used'
          WRITE( 6, '(a)' ) ' to weight the O3 column within the'
          WRITE( 6, '(a)' ) ' chemistry grid and O3 from met or TOMS'
          WRITE( 6, '(a)' ) ' will be used outside the chemistry grid.'
       ENDIF
    ENDIF

    ! FORMAT statements
100 FORMAT( A, L5  )
110 FORMAT( A, F4.2 )

  END SUBROUTINE READ_CHEMISTRY_MENU
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_radiation_menu
!
! !DESCRIPTION: Subroutine READ\_RADIATION\_MENU reads the RADIATION
! MENU section of the GEOS-Chem input file.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE READ_RADIATION_MENU( Input_Opt, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE Input_Opt_Mod, ONLY : OptInput
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(INOUT) :: Input_Opt   ! Input options
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure
!
! !REMARKS:
!  Flux outputs are now scheduled in the HISTORY.rc file, and the relevant
!  fields of Input_Opt will be populated in the RRTMG module routine
!  Init_RRTMG_Indices (called at startup).
!
! !REVISION HISTORY:
!  18 Jun 2013 - D. Ridley   - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER            :: N

    ! Strings
    CHARACTER(LEN=255) :: ErrMsg, ThisLoc

    ! Arrays
    CHARACTER(LEN=255) :: SUBSTRS(MAXDIM)

    !=================================================================
    ! READ_RADIATION_MENU begins here!
    !=================================================================

    ! Initialize
    RC      = GC_SUCCESS
    ErrMsg  = 'Error reading the "input.geos" file!'
    ThisLoc = ' -> at Read_Radiation_Menu (in module GeosCore/input_mod.F90)'

    ! Error check
    IF ( CT1 /= 2 ) THEN
       ErrMsg = 'SIMULATION MENU & ADVECTED SPECIES MENU ' // &
                'must be read in first!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    ! AOD wavelength selection?
    CALL SPLIT_ONE_LINE( SUBSTRS, Input_Opt%NWVSELECT, -1, 'Wavelengths', RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    DO N = 1, Input_Opt%NWVSELECT
       READ( SUBSTRS(N), * ) Input_Opt%WVSELECT(N)
       ! save the string version also
       Input_Opt%STRWVSELECT(N) = TRIM(SUBSTRS(N))
    ENDDO

    ! Turn on RRTMG?
    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'LRAD', RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    READ( SUBSTRS(1:N), * ) Input_Opt%LRAD

    ! Turn on LW radiation calculation?
    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'LLWRAD', RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    READ( SUBSTRS(1:N), * ) Input_Opt%LLWRAD

    ! Turn on SW radiation calculation?
    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'LSWRAD', RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    READ( SUBSTRS(1:N), * ) Input_Opt%LSWRAD

    ! Calculate for clear-sky?
    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'Clear sky', RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    READ( SUBSTRS(1:N), * ) Input_Opt%LSKYRAD(1)

    ! Calculate for all-sky?
    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'All Sky', RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    READ( SUBSTRS(1:N), * ) Input_Opt%LSKYRAD(2)

    ! Radiation timestep?
    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'TS_RAD', RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    READ( SUBSTRS(1:N), * ) Input_Opt%TS_RAD

    ! Separator line
    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'separator 1', RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    !=================================================================
    ! Error check settings
    !=================================================================

    ! Use of RRTMG necessitates recompilation
#ifdef RRTMG
    IF ( .not. Input_Opt%LRAD ) THEN
       ErrMsg = 'LRAD=F but RRTMG is defined at compile time!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
#else
    IF ( Input_Opt%LRAD ) THEN
       ErrMsg = 'LRAD=T but RRTMG is not defined at compile time!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
#endif

    ! Make sure radiation switches are turned off if RRTMG is off
    IF ( ( .not. Input_Opt%LRAD ) .and. Input_Opt%LLWRAD ) THEN
       ErrMsg = 'Cannot have LW fluxes turned on without RRTMG'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    IF ( ( .not. Input_Opt%LRAD ) .and. Input_Opt%LSWRAD ) THEN
       ErrMsg = 'Cannot have SW fluxes turned on without RRTMG'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    IF ( ( .not. Input_Opt%LRAD ) .and. Input_Opt%LSKYRAD(1) ) THEN
       ErrMsg = 'Cannot have clear-sky flux turned on without RRTMG'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    IF ( ( .not. Input_Opt%LRAD ) .and. Input_Opt%LSKYRAD(2) ) THEN
       ErrMsg = 'Cannot have all-sky flux turned on without RRTMG'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
    ENDIF

    !=================================================================
    ! Print to screen
    !=================================================================
    IF ( Input_Opt%amIRoot ) THEN
       WRITE( 6, '(/,a)' ) 'RADIATION MENU'
       WRITE( 6, '(  a)' ) '--------------'
       DO N=1, Input_Opt%NWVSELECT
          WRITE( 6, 115     ) 'AOD output wavelength (nm)  : ', &
                               Input_Opt%WVSELECT(N)
       ENDDO
       WRITE( 6, 100 ) 'Turn on radiation?          : ', &
                        Input_Opt%LRAD
       WRITE( 6, 100 ) 'Consider LW                 : ', &
                        Input_Opt%LLWRAD
       WRITE( 6, 100 ) 'Consider SW                 : ', &
                        Input_Opt%LSWRAD
       WRITE( 6, 125 ) 'Clear-sky/All-sky           : ', &
                        Input_Opt%LSKYRAD(1), '/',       &
                        Input_Opt%LSKYRAD(2)
       WRITE( 6, 110 ) 'Radiation timestep [sec]    : ', &
                        Input_Opt%TS_RAD
    ENDIF

    ! FORMAT statements
100 FORMAT( A, L5           )
110 FORMAT( A, I5           )
115 FORMAT( A, F7.1         )
120 FORMAT( A, 11I1         )
125 FORMAT( A, L5, A, L5    )
130 FORMAT( A, 12( A2, 1x ) )

  END SUBROUTINE READ_RADIATION_MENU
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_transport_menu
!
! !DESCRIPTION: Subroutine READ\_TRANSPORT\_MENU reads the TRANSPORT MENU
!  section of the GEOS-Chem input file.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE READ_TRANSPORT_MENU( Input_Opt, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE Input_Opt_Mod, ONLY : OptInput
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(INOUT) :: Input_Opt   ! Input options
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure
!
! !REVISION HISTORY:
!  20 Jul 2004 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER            :: N

    ! Strings
    CHARACTER(LEN=255) :: ErrMsg, ThisLoc

    ! Arrays
    CHARACTER(LEN=255) :: SUBSTRS(MAXDIM)

    !=================================================================
    ! READ_TRANSPORT_MENU begins here!
    !=================================================================

    ! Initialize
    RC      = GC_SUCCESS
    ErrMsg  = 'Error reading the "input.geos" file!'
    ThisLoc = ' -> at Read_Transport_Menu (in module GeosCore/input_mod.F90)'

    ! Error check
    IF ( CT1 /= 2 ) THEN
       ErrMSg = 'SIMULATION MENU & ADVECTED SPECIES MENU ' // &
                'must be read in first!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    ! Turn on transport?
    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'LTRAN', RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    READ( SUBSTRS(1:N), * ) Input_Opt%LTRAN

    ! Fill negative values
    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'LFILL', RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    READ( SUBSTRS(1:N), * ) Input_Opt%LFILL

    ! IORD, JORD, KORD
    CALL SPLIT_ONE_LINE( SUBSTRS, N, 3, 'IORD, JORD, KORD', RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    READ( SUBSTRS(1:N), * ) Input_Opt%TPCORE_IORD, &
                            Input_Opt%TPCORE_JORD, &
                            Input_Opt%TPCORE_KORD

    ! Separator line
    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'separator 1', RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    !=================================================================
    ! Print to screen
    !=================================================================
    IF ( Input_Opt%amIRoot ) THEN
       WRITE( 6, '(/,a)' ) 'TRANSPORT MENU'
       WRITE( 6, '(  a)' ) '--------------'
       WRITE( 6, 100     ) 'Turn on transport?          : ', &
                            Input_Opt%LTRAN
       WRITE( 6, 100     ) 'Let TPCORE Fill negatives?  : ', &
                            Input_Opt%LFILL
       WRITE( 6, 110     ) 'IORD, JORD, KORD for TPCORE?: ', &
                            Input_Opt%TPCORE_IORD,           &
                            Input_Opt%TPCORE_JORD,           &
                            Input_Opt%TPCORE_KORD
    ENDIF

    ! FORMAT statements
100 FORMAT( A, L5  )
110 FORMAT( A, 5I5 )

  END SUBROUTINE READ_TRANSPORT_MENU
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_photolysis_menu
!
! !DESCRIPTION: Subroutine READ\_PHOTOLYSIS\_MENU reads the SIMULATION MENU
!  section of the GEOS-Chem input file.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE READ_PHOTOLYSIS_MENU( Input_Opt, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE Input_Opt_Mod, ONLY : OptInput
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(INOUT) :: Input_Opt   ! Input options
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure
!
! !REVISION HISTORY:
!  26 Jul 2019 - T. Sherwen - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER            :: N, C

    ! Strings
    CHARACTER(LEN=255) :: ErrMsg, ThisLoc

    ! Arrays
    CHARACTER(LEN=255) :: SUBSTRS(MAXDIM)

    !=================================================================
    ! READ_PHOTOLYSIS_MENU begins here!
    !=================================================================

    ! Initialize
    RC      = GC_SUCCESS
    ErrMsg  = 'Error reading the "input.geos" file!'
    ThisLoc = ' -> at Read_Photolysis_Menu (in GeosCore/input_mod.F90)'

    !-----------------------------------------------------------------
    ! Photolysis directory
    !-----------------------------------------------------------------

    ! Root data dir
    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'FAST_JX_DIR', RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    READ( SUBSTRS(1:N), '(a)' ) Input_Opt%FAST_JX_DIR

    ! Make sure PHOTOLYSIS-DIR ends with a "/" character
    C = LEN_TRIM( Input_Opt%FAST_JX_DIR )
    IF ( Input_Opt%FAST_JX_DIR(C:C) /= '/' ) THEN
       Input_Opt%FAST_JX_DIR = TRIM( Input_Opt%FAST_JX_DIR ) // '/'
    ENDIF

    ! Separator line
    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'separator 1', RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    ! Return success
    RC = GC_SUCCESS

    !=================================================================
    ! Print to screen
    !=================================================================
    IF ( Input_Opt%amIRoot ) THEN
       WRITE( 6, '(/,a)' ) 'PHOTOLYSIS MENU'
       WRITE( 6, '(  a)' ) '---------------'
       WRITE( 6, 110 ) 'FAST_JX Directory           : ', &
                        TRIM( Input_Opt%FAST_JX_DIR )
    ENDIF

    ! Format statements
100 FORMAT( A, I8.8, 1X, I6.6 )
110 FORMAT( A, A              )

  END SUBROUTINE READ_PHOTOLYSIS_MENU
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_convection_menu
!
! !DESCRIPTION: Subroutine READ\_CONVECTION\_MENU reads the CONVECTION MENU
!  section of the GEOS-Chem input file.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE READ_CONVECTION_MENU( Input_Opt, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE Input_Opt_Mod,      ONLY : OptInput
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(INOUT) :: Input_Opt   ! Input options
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure
!
! !REVISION HISTORY:
!  20 Jul 2004 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER            :: N

    ! Strings
    CHARACTER(LEN=255) :: ErrMsg, ThisLoc

    ! Arrays
    CHARACTER(LEN=255) :: SUBSTRS(MAXDIM)

    !=================================================================
    ! READ_CONVECTION_MENU begins here!
    !=================================================================

    ! Initialize
    RC      = GC_SUCCESS
    ThisLoc = ' -> at Read_Convection_Menu (in module GeosCore/input_mod.F90)'

    ! Error check
    IF ( CT1 /= 2 ) THEN
       ErrMsg = 'SIMULATION MENU & ADVECTED SPECIES MENU ' // &
                'must be read in first!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    ErrMsg  = 'Error reading the "input.geos" file!'

    ! Turn on convection?
    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'LCONV', RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    READ( SUBSTRS(1:N), * ) Input_Opt%LCONV

    ! Turn on BL mixing
    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'LTURB', RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    READ( SUBSTRS(1:N), * ) Input_Opt%LTURB

    ! Turn on non-local PBL scheme (Lin, 03/31/09)
    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'LNLPBL', RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    READ( SUBSTRS(1:N), * ) Input_Opt%LNLPBL

    ! Separator line
    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'separator', RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    ! Return success
    RC = GC_SUCCESS

    !=================================================================
    ! Print to screen
    !=================================================================
    IF ( Input_Opt%amIRoot ) THEN
       WRITE( 6, '(/,a)' ) 'CONVECTION MENU'
       WRITE( 6, '(  a)' ) '----------------'
       WRITE( 6, 100     ) 'Turn on cloud convection?   : ', &
                            Input_Opt%LCONV
       WRITE( 6, 100     ) 'Turn on PBL mixing?         : ', &
                            Input_Opt%LTURB
       WRITE( 6, 100     ) 'Turn on non-local PBL?      : ', &
                            Input_Opt%LNLPBL
    ENDIF

    ! FORMAT statements
100 FORMAT( A, L5 )

  END SUBROUTINE READ_CONVECTION_MENU
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_deposition_menu
!
! !DESCRIPTION: Subroutine READ\_DEPOSITION\_MENU reads the DEPOSITION MENU
!  section of the GEOS-Chem input file.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE READ_DEPOSITION_MENU( Input_Opt, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE Input_Opt_Mod, ONLY : OptInput
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(INOUT) :: Input_Opt   ! Input options
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure
!
! !REVISION HISTORY:
!  20 Jul 2004 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER            :: N

    ! Strings
    CHARACTER(LEN=255) :: ErrMsg, ThisLoc

    ! Arrays
    CHARACTER(LEN=255) :: SUBSTRS(MAXDIM)

    !=================================================================
    ! READ_DEPOSITION_MENU begins here!
    !=================================================================

    ! Initialize
    RC      = GC_SUCCESS
    ErrMsg  = 'Error reading the "input.geos" file!'
    ThisLoc = ' -> at Read_Deposition_Menu (in module GeosCore/input_mod.F90)'

    ! Error check
    IF ( CT1 /= 2 ) THEN
       ErrMsg = 'SIMULATION MENU & ADVECTED SPECIES MENU ' // &
                'must be read in first!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    ! Turn on drydep?
    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'LDRYD', RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    READ( SUBSTRS(1:N), * ) Input_Opt%LDRYD

    ! Turn on wetdep?
    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'LWETD', RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    READ( SUBSTRS(1:N), * ) Input_Opt%LWETD

    ! Turn on CO2 effect on drydep?
    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'CO2_EFFECT', RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    READ( SUBSTRS(1:N), * ) Input_Opt%CO2_EFFECT

    ! CO2 level at simulation
    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'CO2_LEVEL', RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    READ( SUBSTRS(1:N), * ) Input_Opt%CO2_LEVEL

    ! Reference CO2 level
    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'CO2_REF', RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    READ( SUBSTRS(1:N), * ) Input_Opt%CO2_REF

    ! Reference CO2 level
    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'RA_Alt_Above_Sfc', RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    READ( SUBSTRS(1:N), * ) Input_Opt%RA_Alt_Above_Sfc

    ! Separator line
    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'separator', RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    !=================================================================
    ! Error check settings
    !=================================================================

    ! Turn off drydep for simulations that don't need it
    IF ( Input_Opt%ITS_A_TAGCO_SIM   ) Input_Opt%LDRYD = .FALSE.

    ! Turn off wetdep for simulations that don't need it
    IF ( Input_Opt%ITS_A_TAGO3_SIM   ) Input_Opt%LWETD = .FALSE.
    IF ( Input_Opt%ITS_A_TAGCO_SIM   ) Input_Opt%LWETD = .FALSE.
    IF ( Input_Opt%ITS_A_CH4_SIM     ) Input_Opt%LWETD = .FALSE.

    ! Set the PBL drydep flag. This determines if dry deposition is
    ! applied (and drydep frequencies are calculated) over the entire
    ! PBL or the first model layer only. For now, set this value
    ! automatically based upon the selected PBL scheme: 1st model layer
    ! for the non-local PBL scheme, full PBL for the full-mixing scheme.
    IF ( Input_Opt%LNLPBL ) THEN
       Input_Opt%PBL_DRYDEP = .FALSE.
    ELSE
       Input_Opt%PBL_DRYDEP = .TRUE.
    ENDIF

    ! If CO2 effect on RS in turned on, calculate the scaling factor
    ! on Rs based on Franks et al. (2013) (ayhwong, 6/25/2019)
    If (Input_Opt%CO2_EFFECT) THEN
       Input_Opt%RS_SCALE = Input_Opt%CO2_LEVEL / Input_Opt%CO2_REF * &
                           (Input_Opt%CO2_LEVEL + 80.0_fp) *          &
                           (Input_Opt%CO2_REF   - 40.0_fp) /          &
                           (Input_Opt%CO2_LEVEL - 40.0_fp) /          &
                           (Input_Opt%CO2_REF   + 80.0_fp)
    ENDIF

    !=================================================================
    ! Print to screen
    !=================================================================
    IF ( Input_Opt%amIRoot ) THEN
       WRITE( 6, '(/,a)' ) 'DEPOSITION MENU'
       WRITE( 6, '(  a)' ) '---------------'
       WRITE( 6, 100     ) 'Turn on dry deposition?     : ', &
                            Input_Opt%LDRYD
       WRITE( 6, 100     ) 'Dry dep over full PBL?      : ', &
                            Input_Opt%PBL_DRYDEP
       WRITE( 6, 100     ) 'Turn on wet deposition?     : ', &
                            Input_Opt%LWETD
       WRITE( 6, 100     ) 'Turn on CO2 effect?         : ', &
                            Input_Opt%CO2_EFFECT
       WRITE( 6, 110     ) 'CO2 level                   : ', &
                            Input_Opt%CO2_LEVEL
       WRITE( 6, 110     ) 'CO2 reference level         : ', &
                            Input_Opt%CO2_REF
       WRITE( 6, 110     ) 'RIX scaling factor          : ', &
                            Input_Opt%RS_SCALE
    ENDIF

    ! FORMAT statements
100 FORMAT( A, L5 )
110 FORMAT( A, f6.2 )

  END SUBROUTINE READ_DEPOSITION_MENU
!EOC
#if !(defined( EXTERNAL_GRID ) || defined( EXTERNAL_FORCING ))
#ifdef BPCH_DIAG
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_gamap_menu
!
! !DESCRIPTION: Subroutine READ\_GAMAP\_MENU reads the GAMAP MENU section
!  of the GEOS-Chem input file.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE READ_GAMAP_MENU( Input_Opt, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE Input_Opt_Mod, ONLY : OptInput
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(INOUT) :: Input_Opt   ! Input options
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure
!
! !REVISION HISTORY:
!  25 Apr 2005 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER            :: N

    ! Strings
    CHARACTER(LEN=255) :: ErrMsg, ThisLoc

    ! Arrays
    CHARACTER(LEN=255) :: SUBSTRS(MAXDIM)

    !=================================================================
    ! READ_GAMAP_MENU begins here!
    !=================================================================

    ! Initialize
    RC      = GC_SUCCESS
    ErrMsg  = 'Error reading the "input.geos" file!'
    ThisLoc = ' -> at Read_Gamap_Menu (in module GeosCore/input_mod.F90)'

    ! Background
    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'DIAGINFO', RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    READ( SUBSTRS(1:N), '(a)' ) Input_Opt%GAMAP_DIAGINFO

    ! Redirect
    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'TRACERINFO', RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    READ( SUBSTRS(1:N), '(a)' ) Input_Opt%GAMAP_TRACERINFO

    ! Separator line
    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'separator', RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    !=================================================================
    ! Print to screen
    !=================================================================
    IF ( Input_Opt%amIRoot ) THEN
       WRITE( 6, '(/,a)' ) 'GAMAP MENU'
       WRITE( 6, '(  a)' ) '---------------'
       WRITE( 6, '(a,a)' ) 'GAMAP "diaginfo.dat"   file : ', &
                            TRIM( Input_Opt%GAMAP_DIAGINFO   )
       WRITE( 6, '(a,a)' ) 'GAMAP "tracerinfo.dat" file : ', &
                            TRIM( Input_Opt%GAMAP_TRACERINFO )
    ENDIF

  END SUBROUTINE READ_GAMAP_MENU
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_output_menu
!
! !DESCRIPTION: Subroutine READ\_OUTPUT\_MENU reads the OUTPUT MENU section of
!  the GEOS-Chem input file.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE READ_OUTPUT_MENU( Input_Opt, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE Input_Opt_Mod,  ONLY : OptInput
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(INOUT) :: Input_Opt   ! Input options
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure
!
! !REVISION HISTORY:
!  20 Jul 2004 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER            :: IOS

    ! Strings
    CHARACTER(LEN=255) :: ErrMsg, ThisLoc

    !=================================================================
    ! READ_OUTPUT_MENU begins here!
    !=================================================================

    ! Initialize
    RC      = GC_SUCCESS
    ErrMsg  = 'Error reading the "input.geos" file!'
    ThisLoc = ' -> at Read_Output_Menu (in module GeosCore/input_mod.F90)'

    ! Read info
    READ( IU_GEOS, 100, IOSTAT=IOS ) Input_Opt%NJDAY
100 FORMAT( 26x, 31i1, /  26x, 29i1, /, 26x, 31i1, /, 26x, 30i1, /, &
            26x, 31i1, /, 26x, 30i1, /, 26x, 31i1, /, 26x, 31i1, /, &
            26x, 30i1, /  26x, 31i1, /, 26x, 30i1, /, 26x, 31i1 )

    ! Trap potential errors
    IF ( IOS /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    !=================================================================
    ! Print to screen
    !=================================================================
    IF( Input_Opt%amIRoot ) THEN
       WRITE( 6, '(/,a)' ) 'OUTPUT MENU'
       WRITE( 6, '(  a)' ) '-----------'
       WRITE( 6, 110     )
       WRITE( 6, 120     )
       WRITE( 6, 130     )
       WRITE( 6, 140     ) Input_Opt%NJDAY
    ENDIF

    ! FORMAT statements
110 FORMAT( '              1111111111222222222233' )
120 FORMAT( '     1234567890123456789012345678901' )
130 FORMAT( '     -------------------------------' )
140 FORMAT( 'JAN--', 31i1, /, 'FEB--', 29i1, /, 'MAR--', 31i1, /, &
            'APR--', 30i1, /, 'MAY--', 31i1, /, 'JUN--', 30i1, /, &
            'JUL--', 31i1, /, 'AUG--', 31i1, /, 'SEP--', 30i1, /, &
            'OCT--', 31i1, /, 'NOV--', 30i1, /, 'DEC--', 31i1 )

    ! Make sure we have output at end of run
    CALL IS_LAST_DAY_GOOD( Input_Opt, RC )

    ! Trap potential errors
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

  END SUBROUTINE READ_OUTPUT_MENU
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_diagnostic_menu
!
! !DESCRIPTION: Subroutine READ\_DIAGNOSTIC\_MENU reads the DIAGNOSTIC MENU
!  section of the GEOS-Chem input file.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE READ_DIAGNOSTIC_MENU( Input_Opt, RC )
!
! !USES:
!
    USE CMN_DIAG_MOD        ! Needed for timeseries diags (binary only)
    USE CMN_SIZE_MOD,  ONLY : NDSTBIN
    USE BPCH2_MOD,     ONLY : OPEN_BPCH2_FOR_WRITE
    USE DIAG03_MOD,    ONLY : ND03,      PD03,      PD03_PL
    USE DIAG53_MOD,    ONLY : ND53,      PD53
    USE DRYDEP_MOD,    ONLY : NUMDEP
    USE ErrCode_Mod
    USE FILE_MOD,      ONLY : IU_BPCH
    USE Input_Opt_Mod, ONLY : OptInput
    USE TIME_MOD,      ONLY : GET_NYMDb, GET_NHMSb, EXPAND_DATE
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(INOUT) :: Input_Opt   ! Input options
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure
!
! !REVISION HISTORY:
!  20 Jul 2004 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER            :: M, N, N_TMP

    ! Strings
    CHARACTER(LEN=255) :: ErrMsg, ThisLoc, BPCH_FILE

    ! Arrays
    CHARACTER(LEN=255) :: SUBSTRS(MAXDIM)

    !=================================================================
    ! READ_DIAGNOSTIC_MENU begins here!
    !=================================================================

    ! Initialize
    RC      = GC_SUCCESS

    ErrMsg  = 'Error reading the "input.geos" file!'
    ThisLoc = ' -> at Read_Diagnostic_Menu (in module GeosCore/input_mod.F90)'

    ! Error check
    IF ( CT1 /= 2 ) THEN
       ErrMsg= 'SIMULATION MENU & ADVECTED SPECIES MENU ' // &
               'must be read in first!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    IF ( .not. Input_Opt%LEMIS ) THEN
       WRITE( 6, '(a)' ) 'WARNING: Emissions are turned off. The'
       WRITE( 6, '(a)' ) ' following diagnostics will also be turned'
       WRITE( 6, '(a)' ) ' off:  ND06, ND53'
    ENDIF

    ! Binary punch file name
    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1,  'BPCH_FILE', RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    READ( SUBSTRS(1:N), '(a)' ) BPCH_FILE

    ! Separator line
    CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'separator', RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    !--------------------------
    ! ND03: Hg diagnostics
    !--------------------------
    CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'ND03', RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    READ( SUBSTRS(1), * ) ND03
    IF ( .not. Input_Opt%ITS_A_MERCURY_SIM ) ND03 = 0
    CALL SET_TINDEX( Input_Opt, 03, ND03, SUBSTRS(2:N), N-1, PD03 )

    !--------------------------
    ! ND06: Dust emissions
    !--------------------------
    CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'ND06', RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
#ifdef TOMAS
    READ( SUBSTRS(1), * ) ND06
    IF ( .not. Input_Opt%LDUST .or. &
         .not. Input_Opt%LEMIS ) ND06 = 0
    CALL SET_TINDEX( Input_Opt, 06, ND06, SUBSTRS(2:N), N-1, NDSTBIN)
#endif

    !--------------------------
    ! ND44 drydep vel & flux
    !--------------------------

    ! Number of species depends on simulation type
    IF ( Input_Opt%ITS_A_TAGO3_SIM ) THEN
       N_TMP = Input_Opt%N_ADVECT
    ELSE
       N_TMP = NUMDEP
    ENDIF

    CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'ND44', RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
#ifdef TOMAS
    READ( SUBSTRS(1), * ) ND44
    IF ( .not. Input_Opt%LDRYD ) ND44 = 0
    CALL SET_TINDEX( Input_Opt, 44, ND44, SUBSTRS(2:N), N-1, N_TMP )
#endif

    !--------------------------
    ! ND53: POPS
    !--------------------------
    CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'ND53', RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    READ( SUBSTRS(1), * ) ND53
    IF ( .not. Input_Opt%ITS_A_POPS_SIM .or. &
         .not. Input_Opt%LEMIS) ND53 = 0
    CALL SET_TINDEX( Input_Opt, 53, ND53, SUBSTRS(2:N), N-1, PD53 )

    !--------------------------
    ! ND59: TOMAS aerosol emiss
    !--------------------------
    CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'ND59', RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
#ifdef TOMAS
    READ( SUBSTRS(1), * ) ND59
    CALL SET_TINDEX( Input_Opt, 59, ND59, SUBSTRS(2:N), N-1, PD59 )
#endif

    !--------------------------
    ! ND60: Wetland Fraction
    ! ND60: TOMAS rate
    !--------------------------
    CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'ND60', RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
#ifdef TOMAS
    READ( SUBSTRS(1), * ) ND60
    CALL SET_TINDEX( Input_Opt, 60, ND60, SUBSTRS(2:N), N-1, PD60 )
#endif

    !--------------------------
    ! ND61: 3-D TOMAS rate
    !--------------------------
    CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'ND61', RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
#ifdef TOMAS
    READ( SUBSTRS(1), * ) ND61
    CALL SET_TINDEX( Input_Opt, 61, ND61, SUBSTRS(2:N), N-1, PD61 )
#endif

    !--------------------------
    ! ND72: Radiation output
    !--------------------------
    CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'ND72', RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
#ifdef RRTMG
    !output fields are nspecies*nradfields but user can only specify
    !rad fields (e.g. SW TOA ALL-SKY) so we set the max to the total
    !divided by number of allowed species (PD72R)
    READ( SUBSTRS(1), * ) ND72
    CALL SET_TINDEX( Input_Opt, 72, ND72, SUBSTRS(2:N), N-1, PD72R )

    !If LRAD is on then ND72 must be on (so the diagnostic is
    !available to write into). Check for this
    IF ( (Input_Opt%LRAD) .AND. (ND72.EQ.0) ) THEN
       ErrMsg = 'If LRAD is true then ' // &
                'ND72 diagnostic must be switched on'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
#endif

    !=================================================================
    ! %%%%% IF BPCH DIAGNOSTICS ARE ACTIVATED (BPCH_DIAG=y) %%%%%
    !
    ! Copy shadow variables to diagnostic variables and
    ! call various bpch diagnostic setup routines
    !=================================================================
    Input_Opt%ND53  = ND53

#ifdef TOMAS
    Input_Opt%ND06  = ND06
    Input_Opt%ND44  = ND44
    Input_Opt%ND59  = ND59
    Input_Opt%ND60  = ND60
    Input_Opt%ND61  = ND61
#endif

#ifdef RRTMG
    Input_Opt%ND72  = ND72
#endif

    ! Loop over # of diagnostics
    DO M = 1, Input_Opt%Max_BPCH_Diag
       Input_Opt%TCOUNT(M)       = TCOUNT(M)
       Input_Opt%TMAX(M)         = TMAX(M)

       ! Loop over tracers per diagnostic
       DO N = 1, Input_Opt%N_ADVECT
          Input_Opt%TINDEX(M,N)  = TINDEX(M,N)
       ENDDO
    ENDDO

    !=================================================================
    ! Call other setup routines
    !================================================================

    ! Expand YYYYMMDD tokens in the bpch file name
    CALL EXPAND_DATE( BPCH_FILE, GET_NYMDb(), GET_NHMSb() )

    ! Find a free file LUN
    IU_BPCH = findFreeLUN()

    ! Open the binary punch file for output
    CALL OPEN_BPCH2_FOR_WRITE( IU_BPCH, BPCH_FILE )

    ! Return success
    RC = GC_SUCCESS

  END SUBROUTINE READ_DIAGNOSTIC_MENU
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: set_tindex
!
! !DESCRIPTION: Subroutine SET\_TINDEX sets the TINDEX and TMAX arrays,
!  which determine how many tracers to print to the punch file.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE SET_TINDEX(Input_Opt, N_DIAG, L_DIAG, SUBSTRS, N, NMAX)
!
! !USES:
!
#ifdef TOMAS
    USE CHARPAK_MOD,   ONLY : TXTEXT   ! (win, 7/14/09)
#endif
    USE CMN_DIAG_MOD        ! TMAX, TINDEX
    USE Input_Opt_Mod, ONLY : OptInput
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput),     INTENT(IN) :: Input_Opt   ! Input Options object
    INTEGER,            INTENT(IN) :: N_DIAG      ! GEOS-Chem diagnostic #
    INTEGER,            INTENT(IN) :: N           ! # of valid substrs passed
    INTEGER,            INTENT(IN) :: NMAX        ! Max # of tracers allowed
    INTEGER,            INTENT(IN) :: L_DIAG      ! # of levels to save
    CHARACTER(LEN=255), INTENT(IN) :: SUBSTRS(N)  ! Substrs passed from
                                                  !  READ_DIAGNOSTIC_MENU
!
! !REMARKS:
!  NOTE: This routine converts to a stub when BPCH_DIAG=n (bmy, 1/16/18)
!
! !REVISION HISTORY:
!  20 Jul 2004 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    LOGICAL, SAVE       :: FIRST = .TRUE.
    LOGICAL             :: IS_ALL
    INTEGER             :: M
#ifdef TOMAS
    INTEGER             :: NN,     COL,     IFLAG, TC     ! (win, 7/14/09)
    CHARACTER (LEN=255) :: WORD,   SUBWORD, TMP1,  TMP2   ! (win, 7/14/09)
    INTEGER             :: MINTMP, MAXTMP                 ! (win, 7/14/09)
#endif

    !=================================================================
    ! SET_TINDEX begins here!
    !=================================================================

    ! Error check
    IF ( N < 1 ) THEN
       IF ( Input_Opt%amIRoot ) THEN
          WRITE( 6, '(a)' ) 'ERROR: N must be 1 or greater!'
          WRITE( 6, '(a)' ) 'STOP in SET_TINDEX (input_mod.F90)'
          WRITE( 6, '(a)' ) REPEAT( '=', 79 )
       ENDIF
       STOP
    ENDIF

    !=================================================================
    ! If the word "all" is present, then set TMAX, TINDEX to all
    ! available tracers for the given diagnostic.  Otherwise, just
    ! use the tracers that were read in from the line
    !=================================================================
    IF ( TRIM( SUBSTRS(1) ) == 'all'  .or. &
         TRIM( SUBSTRS(1) ) == 'ALL' ) THEN

       ! TMAX is the max # of tracers to print out
       TMAX(N_DIAG) = NMAX

       ! Fill TINDEX with all possible diagnostic tracer numbers
       DO M = 1, TMAX(N_DIAG)
          TINDEX(N_DIAG,M) = M
       ENDDO

       ! Set flag
       IS_ALL = .TRUE.

    ELSE

#ifdef TOMAS
!(win, 7/14/09)  use TXTEXT and split the read in characters by -

       COL      = 1
       NN       = 0
       SUBWORD  = ''
       IFLAG    = 0

       ! Use explicit DO-loop
       DO M = 1, N
          WORD = SUBSTRS(M)

          ! Check if the characters are a range with - in the middle
          CALL TXTEXT ( '-', WORD, COL, SUBWORD, IFLAG )

          ! Found a dash!  Get the numbers on both sides of the dash
          ! since these the min and max of the tracer range
          IF ( IFLAG == 0 ) THEN
             TMP1 = TRIM( WORD(      1:COL-1      ) )
             TMP2 = TRIM( WORD( COL+1:LEN_TRIM( WORD ) ) )

             READ( TMP1, * ) MINTMP
             READ( TMP2, * ) MAXTMP

             DO TC = MINTMP, MAXTMP
                NN = NN + 1
                TINDEX( N_DIAG, NN ) = TC
             ENDDO

          ! If we haven't found a dash, then there is only one number,
          ! so that number is both the min and max of the tracer range
          ELSE IF ( IFLAG == -1 ) THEN
             NN = NN + 1
             TMP1 = TRIM( WORD )
             READ( TMP1, * ) TINDEX( N_DIAG, NN )
          ENDIF

       ENDDO

       ! Set TMAX to the counted # of tracers
       TMAX( N_DIAG ) = NN
#else
       ! Otherwise, set TMAX, TINDEX to the # of tracers
       ! listed in "input.ctm" -- need some error checks too
       TMAX(N_DIAG) = N

       ! Use explicit DO-loop
       DO M = 1, N
          READ( SUBSTRS(M:M), * ) TINDEX(N_DIAG,M)
       ENDDO
#endif
       ! Set flag
       IS_ALL = .FALSE.

    ENDIF

    !=================================================================
    ! Print to screen
    !=================================================================

    ! First-time printing only
    IF ( FIRST ) THEN
       IF( Input_Opt%amIRoot ) THEN
          WRITE( 6,'(/,a)' ) 'DIAGNOSTIC MENU'
          WRITE( 6,'(  a)' ) '---------------'
          WRITE( 6,'(  a)' ) 'Diag    L   Tracers being saved to disk'
       ENDIF
       FIRST = .FALSE.
    ENDIF

    ! Test if all tracers are being printed out
    IF ( IS_ALL ) THEN

       ! Print abbreviated output string
       IF ( L_DIAG > 0 ) THEN
          IF ( Input_Opt%amIRoot ) THEN
             WRITE( 6, 100 ) N_DIAG, L_DIAG, 1, TMAX(N_DIAG)
          ENDIF
100       FORMAT( 'ND', i2.2, 2x, i3, 1x, i3, ' -', i3 )
       ENDIF

    ELSE

       ! Or just list each tracer
       ! Print each diagnostic and # of tracers that will print out
       IF ( L_DIAG > 0 ) THEN
          IF ( Input_Opt%amIRoot ) THEN
             WRITE( 6, 110 ) N_DIAG, L_DIAG, &
                             ( TINDEX(N_DIAG,M), M=1,TMAX(N_DIAG) )
          ENDIF
110       FORMAT( 'ND', i2, 2x, i3, 1x, 100i3 )
       ENDIF

    ENDIF

  END SUBROUTINE SET_TINDEX
!EOC
#endif
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_planeflight_menu
!
! !DESCRIPTION: Subroutine READ\_PLANEFLIGHT\_MENU reads the PLANEFLIGHT MENU
!  section of the GEOS-Chem input file.  This turns on the plane flight track
!  diagnostic.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE READ_PLANEFLIGHT_MENU( Input_Opt, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE Input_Opt_Mod,   ONLY : OptInput
    USE PLANEFLIGHT_MOD, ONLY : SET_PLANEFLIGHT
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(INOUT) :: Input_Opt   ! Input options
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure
!
! !REVISION HISTORY:
!  20 Jul 2004 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER            :: N

    ! Strings
    CHARACTER(LEN=255) :: ErrMsg, ThisLoc

    ! Arrays
    CHARACTER(LEN=255) :: SUBSTRS(MAXDIM)

    !=================================================================
    ! READ_PLANEFLIGHT_MENU begins here!
    !=================================================================

    ! Initialize
    RC      = GC_SUCCESS
    ErrMsg  = 'Error reading the "input.geos" file!'
    ThisLoc = ' -> at Read_Planeflight_Menu (in module GeosCore/input_mod.F90)'

    ! Error check
    IF ( CT1 /= 2 ) THEN
       ErrMsg = 'SIMULATION MENU & ADVECTED SPECIES MENU ' // &
                'must be read in first!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    ! Turn on planeflight diagnostic?
    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'Do_Planeflight', RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    READ( SUBSTRS(1:N), * ) Input_Opt%Do_Planeflight

    ! Input file name (w/ flight track data points)
    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'Planeflight_InFile', RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    READ( SUBSTRS(1:N), '(a)' ) Input_Opt%Planeflight_InFile

    ! Output file name
    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'Planeflight_OutFile', RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    READ( SUBSTRS(1:N), '(a)' ) Input_Opt%Planeflight_OutFile

    !=================================================================
    ! Print to screen
    !=================================================================
    IF( Input_Opt%amIRoot ) THEN
       WRITE( 6, '(/,a)' ) 'PLANEFLIGHT MENU'
       WRITE( 6, '(  a)' ) '----------------'
       WRITE( 6, 100 ) 'Turn on planeflight diag?   : ', &
                        Input_Opt%Do_Planeflight
       WRITE( 6, 110 ) 'Flight track input file     : ', &
                        TRIM( Input_Opt%Planeflight_InFile )
       WRITE( 6, 110 ) 'Output file name            : ', &
                        TRIM( Input_Opt%Planeflight_OutFile )
    ENDIF

    ! FORMAT statements
100 FORMAT( A, L5    )
110 FORMAT( A, A     )

    !=================================================================
    ! Call setup routines from other F90 modules
    !=================================================================

    ! Pass variables to "planeflight_mod.F90"
    CALL SET_PLANEFLIGHT( Input_Opt%Do_Planeflight,    &
                          Input_Opt%Planeflight_InFile, &
                          Input_Opt%Planeflight_OutFile )

  END SUBROUTINE READ_PLANEFLIGHT_MENU
#endif
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_obspack_menu
!
! !DESCRIPTION: Subroutine READ\_OBSPACK\_MENU reads the OBSPACK MENU
!  section of the GEOS-Chem input file.  This turns on the ObsPack diagnostic.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE READ_OBSPACK_MENU( Input_Opt, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE Input_Opt_Mod, ONLY : OptInput
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(INOUT) :: Input_Opt   ! Input options
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure
!
! !REVISION HISTORY:
!  20 Jul 2004 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER            :: N, S

    ! Strings
    CHARACTER(LEN=255) :: ErrMsg, ThisLoc

    ! Arrays
    CHARACTER(LEN=255) :: SUBSTRS(MAXDIM)

    !=================================================================
    ! READ_OBSPACK_MENU begins here!
    !=================================================================

    ! Initialize
    RC      = GC_SUCCESS
    ErrMsg  = 'Error reading the "input.geos" file!'
    ThisLoc = ' -> at Read_ObsPack_Menu (in module GeosCore/input_mod.F90)'

    ! Error check
    IF ( CT1 /= 2 ) THEN
       ErrMsg = 'SIMULATION MENU & ADVECTED SPECIES MENU ' // &
                'must be read in first!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    ! Turn on ObsPack diagnostic?
    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'Do_ObsPack', RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    READ( SUBSTRS(1:N), * ) Input_Opt%Do_ObsPack

    ! ObsPack quiet output?
    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'Do_ObsPack', RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    READ( SUBSTRS(1:N), * ) Input_Opt%ObsPack_Quiet

    ! Input file name (w/ coordinates and sampling strategy)
    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'ObsPack_InputFile', RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    READ( SUBSTRS(1:N), '(a)' ) Input_Opt%ObsPack_InputFile

    ! Output file name
    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'ObsPack_OutputFile', RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    READ( SUBSTRS(1:N), '(a)' ) Input_Opt%ObsPack_OutputFile

    ! Species names
    CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'ObsPack_SpecName', RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    ! Populate the ObsPack species name list
    Input_Opt%ObsPack_nSpc = N
    DO S = 1, Input_Opt%ObsPack_nSpc
       Input_Opt%ObsPack_SpcName(S) = TRIM( SUBSTRS(S) )
    ENDDO

    !=================================================================
    ! Print to screen
    !=================================================================
    IF( Input_Opt%amIRoot ) THEN
       WRITE( 6, '(/,a)' ) 'OBSPACK_MENU'
       WRITE( 6, '(  a)' ) '----------------'
       WRITE( 6, 100 ) 'Turn on ObsPack diagnostic? : ', &
                        Input_Opt%Do_ObsPack
       WRITE( 6, 100 ) 'Suppress logfile output?    : ', &
                        Input_Opt%ObsPack_Quiet
       WRITE( 6, 110 ) 'ObsPack input file          : ', &
                        TRIM( Input_Opt%ObsPack_InputFile  )
       WRITE( 6, 110 ) 'ObsPack output file         : ', &
                        TRIM( Input_Opt%ObsPack_OutputFile )
    ENDIF

    ! FORMAT statements
100 FORMAT( A, L5    )
110 FORMAT( A, A     )

  END SUBROUTINE READ_OBSPACK_MENU
!EOC

#if !(defined( EXTERNAL_GRID ) || defined( EXTERNAL_FORCING ))
#ifdef BPCH_DIAG
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_nd51_menu
!
! !DESCRIPTION: Subroutine READ\_ND51\_MENU reads the ND51 MENU section of
!  the GEOS-Chem input file.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE READ_ND51_MENU( Input_Opt, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE Input_Opt_Mod,      ONLY : OptInput
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(INOUT) :: Input_Opt   ! Input options
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure
!
! !REVISION HISTORY:
!  20 Jul 2004 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER            :: N

    ! Strings
    CHARACTER(LEN=255) :: ErrMsg, ThisLoc

    ! Arrays
    CHARACTER(LEN=255) :: SUBSTRS(MAXDIM)

    !=================================================================
    ! READ_ND51_MENU begins here!
    !=================================================================

    ! Initialize
    RC      = GC_SUCCESS

    ErrMsg  = 'Error reading the "input.geos" file!'
    ThisLoc = ' -> at Read_ND51_Menu (in module GeosCore/input_mod.F90)'

    ! Error check
    IF ( CT1 /= 2 ) THEN
       ErrMsg = 'SIMULATION MENU & ADVECTED SPECIES MENU ' // &
                'must be read in first!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    ! Turn on ND51 diagnostic
    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1,  'read_nd51_menu:1', RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    READ( SUBSTRS(1:N), * ) Input_Opt%DO_ND51

    ! Instantaneous 3-D timeseries file
    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1,  'read_nd51_menu:2', RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    READ( SUBSTRS(1:N), '(a)' ) Input_Opt%ND51_FILE

    ! Output as hdf
    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1,  'read_nd51_menu:3', RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    READ( SUBSTRS(1:N), * ) Input_Opt%LND51_HDF

    ! (bmy, 12/21/09)
    Input_Opt%LND51_HDF = .FALSE.

    ! Tracers to include
    CALL SPLIT_ONE_LINE( SUBSTRS, Input_Opt%N_ND51, -1, &
                         'read_nd51_menu:4', RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    DO N = 1, Input_Opt%N_ND51
       READ( SUBSTRS(N), * ) Input_Opt%ND51_TRACERS(N)
    ENDDO

    ! NHMS_WRITE
    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1,  'read_nd51_menu:6', RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    READ( SUBSTRS(1:N), * ) Input_Opt%ND51_HR_WRITE

    ! Make sure ND51_HR_WRITE is in the range 0-23.999 hrs
    Input_Opt%ND51_HR_WRITE = MOD( Input_Opt%ND51_HR_WRITE, 24e+0_fp )

    ! HR1, HR2
    CALL SPLIT_ONE_LINE( SUBSTRS, N, 2,  'read_nd51_menu:7', RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    READ( SUBSTRS(1:N), * ) Input_Opt%ND51_HR1, Input_Opt%ND51_HR2

    ! IMIN, IMAX
    CALL SPLIT_ONE_LINE( SUBSTRS, N, 2,  'read_nd51_menu:8', RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    READ( SUBSTRS(1:N), * ) Input_Opt%ND51_IMIN, Input_Opt%ND51_IMAX

    ! JMIN, JMAX
    CALL SPLIT_ONE_LINE( SUBSTRS, N, 2,  'read_nd51_menu:9', RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    READ( SUBSTRS(1:N), * ) Input_Opt%ND51_JMIN, Input_Opt%ND51_JMAX

    ! LMIN, LMAX
    CALL SPLIT_ONE_LINE( SUBSTRS, N, 2,  'read_nd51_menu:10', RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    READ( SUBSTRS(1:N), * ) Input_Opt%ND51_LMIN, Input_Opt%ND51_LMAX

    ! Separator line
    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1,  'read_nd51_menu:11', RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    !=================================================================
    ! Print to screen
    !=================================================================
    IF ( Input_Opt%amIRoot ) THEN
       WRITE( 6, '(/,a)' ) 'ND51 MORNING OR AFTERNOON TIMESERIES MENU'
       WRITE( 6, '(  a)' ) '-----------------------------------------'
       WRITE( 6, 100 ) 'Turn on ND51 timeseries?    : ', &
                        Input_Opt%DO_ND51
       WRITE( 6, 110 ) 'ND51 timeseries file name   : ', &
                        TRIM( Input_Opt%ND51_FILE )
       WRITE( 6, 100 ) 'Output as HDF?              : ', &
                        Input_Opt%LND51_HDF
       WRITE( 6, 120 ) 'ND41 timeseries tracers     : ',  &
                        ( Input_Opt%ND51_TRACERS(N), N=1, &
                          Input_Opt%N_ND51 )
       WRITE( 6, 140 ) 'ND51 hour to write to disk  : ', &
                        Input_Opt%ND51_HR_WRITE
       WRITE( 6, 140 ) 'ND51 averaging period [GMT] : ', &
                        Input_Opt%ND51_HR1,  Input_Opt%ND51_HR2
       WRITE( 6, 130 ) 'ND51 longitude limits       : ', &
                        Input_Opt%ND51_IMIN, Input_Opt%ND51_IMAX
       WRITE( 6, 130 ) 'ND51 latitude  limits       : ', &
                        Input_Opt%ND51_JMIN, Input_Opt%ND51_JMAX
       WRITE( 6, 130 ) 'ND51 altitude  limits       : ', &
                        Input_Opt%ND51_LMIN, Input_Opt%ND51_LMAX
    ENDIF

    ! FORMAT statements
100 FORMAT( A, L5    )
110 FORMAT( A, A     )
120 FORMAT( A, 100I4 )
130 FORMAT( A, 2I5   )
140 FORMAT( A, 2F5.1 )

  END SUBROUTINE READ_ND51_MENU
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_nd51b_menu
!
! !DESCRIPTION: Subroutine READ\_ND51b\_MENU reads the ND51 MENU section
!  of the GEOS-Chem input file.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE READ_ND51b_MENU( Input_Opt, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE Input_Opt_Mod, ONLY : OptInput
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(INOUT) :: Input_Opt   ! Input options
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure
!
! !REVISION HISTORY:
!  21 Dec 2009 - Aaron van D - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER              :: N

    ! Strings
    CHARACTER(LEN=255)   :: ErrMsg, ThisLoc

    ! Arrays
    CHARACTER(LEN=255)   :: SUBSTRS(MAXDIM)

    !=================================================================
    ! READ_ND51b_MENU begins here!
    !=================================================================

    ! Initialize
    RC      = GC_SUCCESS

    ErrMsg  = 'Error reading the "input.geos" file!'
    ThisLoc = ' -> at Read_ND51b_MENU  (in module GeosCore/input_mod.F90)'

    ! Error check
    IF ( CT1 /= 2 ) THEN
       ErrMsg = 'SIMULATION MENU & ADVECTED SPECIES MENU ' // &
                'must be read in first!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    ! Turn on ND51b diagnostic
    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1,  'read_nd51b_menu:1', RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    READ( SUBSTRS(1:N), * ) Input_Opt%DO_ND51b

    ! Instantaneous 3-D timeseries file
    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1,  'read_nd51b_menu:2', RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    READ( SUBSTRS(1:N), '(a)' ) Input_Opt%ND51b_FILE

    ! Output as hdf
    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1,  'read_nd51b_menu:3', RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    READ( SUBSTRS(1:N), * ) Input_Opt%LND51b_HDF

    Input_Opt%LND51b_HDF = .FALSE.

    ! Tracers to include
    CALL SPLIT_ONE_LINE( SUBSTRS, Input_Opt%N_ND51b, -1, &
                         'read_nd51b_menu:4',RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    DO N = 1, Input_Opt%N_ND51b
       READ( SUBSTRS(N), * ) Input_Opt%ND51b_TRACERS(N)
    ENDDO

    ! NHMS_WRITE
    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1,  'read_nd51b_menu:5', RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    READ( SUBSTRS(1:N), * ) Input_Opt%ND51b_HR_WRITE

    ! Make sure ND51b_HR_WRITE is in the range 0-23.999 hrs
    Input_Opt%ND51b_HR_WRITE = MOD(Input_Opt%ND51b_HR_WRITE,24e+0_fp)

    ! HR1, HR2
    CALL SPLIT_ONE_LINE( SUBSTRS, N, 2,  'read_nd51b_menu:6', RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    READ( SUBSTRS(1:N), * ) Input_Opt%ND51b_HR1, Input_Opt%ND51b_HR2

    ! IMIN, IMAX
    CALL SPLIT_ONE_LINE( SUBSTRS, N, 2,  'read_nd51b_menu:7', RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    READ( SUBSTRS(1:N), * ) Input_Opt%ND51b_IMIN, Input_Opt%ND51b_IMAX

    ! JMIN, JMAX
    CALL SPLIT_ONE_LINE( SUBSTRS, N, 2,  'read_nd51b_menu:8', RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    READ( SUBSTRS(1:N), * ) Input_Opt%ND51b_JMIN, Input_Opt%ND51b_JMAX

    ! LMIN, LMAX
    CALL SPLIT_ONE_LINE( SUBSTRS, N, 2,  'read_nd51b_menu:9', RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    READ( SUBSTRS(1:N), * ) Input_Opt%ND51b_LMIN, Input_Opt%ND51b_LMAX

    ! Separator line
    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1,  'read_nd51b_menu:10', RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    ! Return success
    RC = GC_SUCCESS

    !=================================================================
    ! Print to screen
    !=================================================================
    IF ( Input_Opt%amIRoot ) THEN
       WRITE( 6, '(/,a)' )'ND51b MORNING OR AFTERNOON TIMESERIES MENU'
       WRITE( 6, '(  a)' )'-----------------------------------------'
       WRITE( 6, 100 ) 'Turn on ND51b timeseries?    : ', &
                        Input_Opt%DO_ND51
       WRITE( 6, 110 ) 'ND51b timeseries file name   : ', &
                        TRIM( Input_Opt%ND51b_FILE )
       WRITE( 6, 100 ) 'Output as HDF?               : ', &
                        Input_Opt%LND51b_HDF
       WRITE( 6, 120 ) 'ND41 timeseries tracers      : ',  &
                        ( Input_Opt%ND51b_TRACERS(N), N=1, &
                          Input_Opt%N_ND51b )
       WRITE( 6, 140 ) 'ND51b hour to write to disk  : ', &
                        Input_Opt%ND51b_HR_WRITE
       WRITE( 6, 140 ) 'ND51b averaging period [GMT] : ', &
                        Input_Opt%ND51b_HR1,  Input_Opt%ND51b_HR2
       WRITE( 6, 130 ) 'ND51b longitude limits       : ', &
                        Input_Opt%ND51b_IMIN, Input_Opt%ND51b_IMAX
       WRITE( 6, 130 ) 'ND51b latitude  limits       : ', &
                        Input_Opt%ND51b_JMIN, Input_Opt%ND51b_JMAX
       WRITE( 6, 130 ) 'ND51b altitude  limits       : ', &
                        Input_Opt%ND51b_LMIN, Input_Opt%ND51b_LMAX
    ENDIF

    ! FORMAT statements
100 FORMAT( A, L5    )
110 FORMAT( A, A     )
120 FORMAT( A, 100I3 )
130 FORMAT( A, 2I5   )
140 FORMAT( A, 2F5.1 )

  END SUBROUTINE READ_ND51b_MENU
#ifdef TOMAS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_prod_loss_menu
!
! !DESCRIPTION: Subroutine READ\_PROD\_LOSS\_MENU reads the PROD AND LOSS MENU
!  section of the GEOS-Chem input file
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE READ_PROD_LOSS_MENU( Input_Opt, RC )
!
! !USES:
!
    USE CMN_DIAG_MOD
    USE ErrCode_Mod
    USE gckpp_Parameters,   ONLY : NFAM
    USE gckpp_Monitor,      ONLY : FAM_NAMES
    USE Input_Opt_Mod,      ONLY : OptInput
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(INOUT) :: Input_Opt   ! Input options
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure
!
! !REVISION HISTORY:
!  20 Jul 2004 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER            :: F, N, N_ADVECT

    ! Strings
    CHARACTER(LEN=255) :: ErrMsg, ThisLoc

    ! Arrays
    CHARACTER(LEN=255) :: SUBSTRS(MAXDIM)

    !=================================================================
    ! READ_PROD_LOSS_MENU begins here!
    !=================================================================

    ! Initialize
    RC      = GC_SUCCESS

    ErrMsg  = 'Error reading the "input.geos" file!'
    ThisLoc = ' -> at Read_Prod_Loss_Menu (in module GeosCore/input_mod.F90)'

    ! Error check
    IF ( CT1 /= 2 ) THEN
       ErrMsg = 'SIMULATION MENU & ADVECTED SPECIES MENU ' // &
                'must be read in first!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    !=================================================================
    ! Read info about prod & loss families
    !=================================================================

    ! Turn on production & loss diagnostic (e.g. ND65 diagnostic)
    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'DO_SAVE_PL', RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    READ( SUBSTRS(1:N), * ) Input_Opt%DO_SAVE_PL

    ! Read number of levels for ND65 diagnostic
    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'ND65', RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    READ( SUBSTRS(1:N), * ) Input_Opt%ND65

    ! Copy field to variable in CMN_DIAG
    ND65 = Input_Opt%ND65

    !=================================================================
    ! Error check families for certain types of simulations
    !=================================================================

    ! Offline aerosol -- turn off DO_SAVE_PL, since we use ND05,
    ! ND06, ND07, ND08, ND13 etc diagnostics instead of ND65
    IF ( Input_Opt%ITS_AN_AEROSOL_SIM ) THEN
       Input_Opt%DO_SAVE_PL    = .FALSE.
       Input_Opt%ND65          = 0
    ENDIF

    !=================================================================
    ! Set fields of Input Options object
    !=================================================================

    ! Number of advected species
    N_ADVECT = Input_Opt%N_ADVECT

    IF ( Input_Opt%DO_SAVE_PL ) THEN
       IF ( Input_Opt%ITS_A_FULLCHEM_SIM ) THEN
          ! Fullchem - Obtain NFAM from KPP
          Input_Opt%NFAM = NFAM
       ELSEIF ( Input_Opt%ITS_A_TAGO3_SIM ) THEN
          ! Tagged O3
          Input_Opt%NFAM = 2*N_ADVECT
       ELSEIF ( Input_Opt%ITS_A_TAGCO_SIM ) THEN
          ! Tagged CO
          IF ( Input_Opt%LPCO_NMVOC ) THEN
             Input_Opt%NFAM = N_ADVECT+2
          ELSE
             Input_Opt%NFAM = N_ADVECT+6
          ENDIF
       ENDIF
    ENDIF

    ! Return if there are no prod/loss families
    ! or if we have turned off this diagnostic
    IF ( .not. ( Input_Opt%DO_SAVE_PL .and. Input_Opt%NFAM > 0 )) THEN
       Input_Opt%DO_SAVE_PL = .FALSE.
       Input_Opt%ND65       = 0
    ENDIF

    ! Loop over families
    DO F = 1, Input_Opt%NFAM

       IF ( Input_Opt%ITS_A_FULLCHEM_SIM ) THEN

          ! Fullchem - Obtain FAM_NAME from KPP
          Input_Opt%FAM_NAME(F) = FAM_NAMES(F)

       ELSEIF ( Input_Opt%ITS_A_TAGO3_SIM ) THEN

          ! Tagged O3
          IF ( F <= N_ADVECT ) THEN
             Input_Opt%FAM_NAME(F) = &
                  'P' // TRIM(Input_Opt%AdvectSpc_Name(F))
          ELSE
             Input_Opt%FAM_NAME(F) = &
                  'L' // TRIM(Input_Opt%AdvectSpc_Name(F-N_ADVECT))
          ENDIF

       ELSEIF ( Input_Opt%ITS_A_TAGCO_SIM ) THEN

          ! Tagged CO
          IF ( F <= N_ADVECT ) THEN
             Input_Opt%FAM_NAME(F) = 'L'//Input_Opt%AdvectSpc_Name(F)
          ELSEIF ( F == N_ADVECT+1 ) THEN
             Input_Opt%FAM_NAME(F) = 'PCO_CH4'
          ELSEIF ( F == N_ADVECT+2 ) THEN
             Input_Opt%FAM_NAME(F) = 'PCO_NMVOC'
          ELSEIF ( F == N_ADVECT+3 ) THEN
             Input_Opt%FAM_NAME(F) = 'PCO_ISOP'
          ELSEIF ( F == N_ADVECT+4 ) THEN
             Input_Opt%FAM_NAME(F) = 'PCO_CH3OH'
          ELSEIF ( F == N_ADVECT+5 ) THEN
             Input_Opt%FAM_NAME(F) = 'PCO_MONO'
          ELSEIF ( F == N_ADVECT+6 ) THEN
             Input_Opt%FAM_NAME(F) = 'PCO_ACET'
          ENDIF

       ENDIF

       ! Get family type as prod or loss
       IF ( Input_Opt%FAM_NAME(F)(1:1) == 'P'   .or. &
            Input_Opt%FAM_NAME(F)(1:1) == 'p' ) THEN
          Input_Opt%FAM_TYPE(F) = 'prod'
       ELSE
          Input_Opt%FAM_TYPE(F) = 'loss'
       ENDIF

    ENDDO

    !=================================================================
    ! Print to screen
    !=================================================================
    IF ( Input_Opt%amIRoot ) THEN
       WRITE( 6, '(/,a)' ) 'PROD & LOSS DIAGNOSTIC MENU'
       WRITE( 6, '(  a)' ) '---------------------------'
       WRITE( 6, 100 ) 'Turn on prod & loss diag?   : ', &
                        Input_Opt%DO_SAVE_PL
       WRITE( 6, 110 ) '# of levels for P/L diag    : ', &
                        Input_Opt%ND65

       ! Loop over families
       DO F = 1, Input_Opt%NFAM

          ! Write family name and type
          WRITE( 6, 120 ) TRIM(Input_Opt%FAM_NAME(F)), &
                          TRIM(Input_Opt%FAM_TYPE(F))

       ENDDO

    ENDIF

    ! FORMAT statements
100 FORMAT( A, L5 )
110 FORMAT( A, I5 )
120 FORMAT( /, 'Family=', A10, '  Type=', A4 )

  END SUBROUTINE READ_PROD_LOSS_MENU
!EOC
#endif
#endif
#endif
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_mercury_menu
!
! !DESCRIPTION: Subroutine READ\_MERCURY\_MENU reads the BENCHMARK MENU
!  section of the GEOS-Chem input file.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE READ_MERCURY_MENU( Input_Opt, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE Input_Opt_Mod, ONLY : OptInput
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(INOUT) :: Input_Opt   ! Input options
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure
!
! !REVISION HISTORY:
!  24 Feb 2006 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    LOGICAL            :: LDYNOCEAN,      LPREINDHG
    LOGICAL            :: LGTMM,          USE_CHECKS
    LOGICAL            :: LARCTICRIV,     LKRedUV
    INTEGER            :: N
    CHARACTER(LEN=255) :: GTMM_RST_FILE

    ! Strings
    CHARACTER(LEN=255) :: ErrMsg,         ThisLoc

    ! Arrays
    CHARACTER(LEN=255) :: SUBSTRS(MAXDIM)

    !=================================================================
    ! READ_MERCURY_MENU begins here!
    !=================================================================

    ! Initialize
    RC      = GC_SUCCESS
    ErrMsg  = 'Error reading the "input.geos" file!'
    ThisLoc = ' -> at Read_Mercury_Menu (in module GeosCore/input_mod.F90)'

    ! Use error check for tag/tot Hg?
    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'USE_CHECKS', RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    READ( SUBSTRS(1:N), * ) Input_Opt%USE_CHECKS

    ! Use dynamic ocean model?
    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'LDYNOCEAN', RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    READ( SUBSTRS(1:N), * ) Input_Opt%LDYNOCEAN

    ! Use preindustrial simulation?
    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'LPREINDHG', RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    READ( SUBSTRS(1:N), * ) Input_Opt%LPREINDHG

    ! Use GTMM ?
    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'LGTMM', RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    READ( SUBSTRS(1:N), * ) Input_Opt%LGTMM

    ! Name of GTMM restart file
    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'GTMM_RST_FILE', RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    READ( SUBSTRS(1:N), * ) Input_Opt%GTMM_RST_FILE

    ! Use Arctic river Hg?
    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'LARCTICRIV', RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    READ( SUBSTRS(1:N), * ) Input_Opt%LARCTICRIV

    ! Tie reducible HgII(aq) to UV-B?
    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'LKREDUV', RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    READ( SUBSTRS(1:N), * ) Input_Opt%LKRedUV

    ! Separator line
    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'separator', RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    ! Check on logical
    IF (.NOT.( Input_Opt%ITS_A_MERCURY_SIM ) ) THEN
       Input_Opt%LGTMM      = .FALSE.
       Input_Opt%LDYNOCEAN  = .FALSE.
       Input_Opt%LARCTICRIV = .FALSE.
       Input_Opt%LKRedUV    = .FALSE.
    ENDIF

    !=================================================================
    ! Print to screen
    !=================================================================
    IF ( Input_Opt%amIRoot ) THEN
       WRITE( 6, '(/,a)' ) 'MERCURY MENU'
       WRITE( 6, '(  a)' ) '------------'
       WRITE( 6, 110 ) 'Error check tag & total Hg? : ', &
                        Input_Opt%USE_CHECKS
       WRITE( 6, 110 ) 'Use dynamic ocean Hg model? : ', &
                        Input_Opt%LDYNOCEAN
       WRITE( 6, 110 ) 'Preindustrial simulation?   : ', &
                        Input_Opt%LPREINDHG
       WRITE( 6, 110 ) 'Use GTMM ?                  : ', &
                        Input_Opt%LGTMM
       WRITE( 6, 120 ) '=> GTMM restart file        : ', &
                        TRIM( Input_Opt%GTMM_RST_FILE )
       WRITE( 6, 110 ) 'Use Arctic river Hg ?       : ', &
                        Input_Opt%LARCTICRIV
       WRITE( 6, 110 ) 'Tie HgII(aq) red. to UV-B?  : ', &
                        Input_Opt%LKRedUV
    ENDIF

    ! FORMAT statements
100 FORMAT( A, I4  )
110 FORMAT( A, L5  )
120 FORMAT( A, A   )

  END SUBROUTINE READ_MERCURY_MENU
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_ch4_menu
!
! !DESCRIPTION: Subroutine READ\_CH4\_MENU reads the CH4 MENU section of
!  the GEOS-Chem input file; this defines emissions options for CH4 tagged
!  simulations.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE READ_CH4_MENU( Input_Opt, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE Input_Opt_Mod, ONLY : OptInput
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(INOUT) :: Input_Opt   ! Input options
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!
! !REVISION HISTORY:
!  03 Aug 2009 - K. Wecht, C. Pickett-Heaps - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER            :: N

    ! Strings
    CHARACTER(LEN=255) :: ErrMsg, ThisLoc

    ! Arrays
    CHARACTER(LEN=255) :: SUBSTRS(MAXDIM)

    !=================================================================
    ! READ_CH4_MENU begins here!
    !=================================================================

    ! Initialize
    RC      = GC_SUCCESS
    ErrMsg  = 'Error reading the "input.geos" file!'
    ThisLoc = ' -> at Read_CH4_Menu (in module GeosCore/input_mod.F90)'

    ! Use GOSAT CH4 observation operator?
    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'GOSAT_CH4_OBS', RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    READ( SUBSTRS(1:N), * ) Input_Opt%GOSAT_CH4_OBS

    ! Use TCCON CH4 observation operator?
    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'TCCON_CH4_OBS', RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    READ( SUBSTRS(1:N), * ) Input_Opt%TCCON_CH4_OBS

    ! Separator line
    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'separator', RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    !=================================================================
    ! Print to screen
    !=================================================================
    IF ( Input_Opt%amIRoot ) THEN
       WRITE( 6, '(/,a)' ) 'CH4 MENU'
       WRITE( 6, '(  a)' ) '-----------'
       WRITE( 6, 100     ) 'Use GOSAT obs operator: ', &
                            Input_Opt%GOSAT_CH4_OBS
       WRITE( 6, 100     ) 'Use TCCON obs operator: ', &
            Input_Opt%TCCON_CH4_OBS
    ENDIF

    ! FORMAT statements
100 FORMAT( A, L5  )

    ! Return success
    RC = GC_SUCCESS

  END SUBROUTINE READ_CH4_MENU
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_pops_menu
!
! !DESCRIPTION: Subroutine READ\_POPS\_MENU reads the POPS MENU section of
!  the GEOS-Chem input file; this defines emissions options for POPs
!  simulations.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE READ_POPS_MENU( Input_Opt, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE Input_Opt_Mod, ONLY : OptInput
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(INOUT) :: Input_Opt   ! Input options
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!
! !REVISION HISTORY:
!  01 Oct 2012 - C. Friedman - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER            :: N

    ! Strings
    CHARACTER(LEN=255) :: ErrMsg, ThisLoc

    ! Arrays
    CHARACTER(LEN=255) :: SUBSTRS(MAXDIM)

    !=================================================================
    ! READ_POPS_MENU begins here!
    !=================================================================

    ! Initialize
    RC      = GC_SUCCESS
    ErrMsg  = 'Error reading the "input.geos" file!'
    ThisLoc = ' -> at Read_POPs_Menu (in module GeosCore/input_mod.F90)'

    ! POP species
    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'POP_TYPE', RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    READ( SUBSTRS(1:N), * ) Input_Opt%POP_TYPE

    ! Dummy for future process logical switches
    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'CHEM_PROCESS', RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    READ( SUBSTRS(1:N), * ) Input_Opt%CHEM_PROCESS

    ! Molecular weight
    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'POP_XMW', RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    READ( SUBSTRS(1:N), * ) Input_Opt%POP_XMW

    ! KOA
    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'POP_KOA', RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    READ( SUBSTRS(1:N), * ) Input_Opt%POP_KOA

    ! KBC
    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'POP_KBC', RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    READ( SUBSTRS(1:N), * ) Input_Opt%POP_KBC

    ! OH oxidation
    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, ' POP_K_POPG_OH', RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    READ( SUBSTRS(1:N), * ) Input_Opt%POP_K_POPG_OH

    ! O3 oxidation 1
    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, ' POP_K_POPP_O3A', RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    READ( SUBSTRS(1:N), * ) Input_Opt%POP_K_POPP_O3A

    ! O3 oxidation 2
    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'POP_K_POPP_O3B', RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    READ( SUBSTRS(1:N), * ) Input_Opt%POP_K_POPP_O3B

    ! H*
    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'POP_HSTAR', RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    READ( SUBSTRS(1:N), * ) Input_Opt%POP_HSTAR

    ! DEL_H
    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'POP_DEL_H', RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    READ( SUBSTRS(1:N), * ) Input_Opt%POP_DEL_H

    ! DEL_Hw
    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'POP_DEL_Hw', RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    READ( SUBSTRS(1:N), * ) Input_Opt%POP_DEL_Hw

    ! Separator line
    CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'separator', RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    !=================================================================
    ! Print to screen
    !=================================================================
    IF ( Input_Opt%amIRoot ) THEN
       WRITE( 6, '(/,a)' ) 'POPS MENU'
       WRITE( 6, '(  a)' ) '------------'
       WRITE( 6, 120     ) 'Species of POP        : ', &
                            Input_Opt%POP_TYPE
       WRITE( 6, 110     ) 'Chemistry on?         : ', &
                            Input_Opt%CHEM_PROCESS
       WRITE( 6, 130     ) 'POP_XMW               : ', &
                            Input_Opt%POP_XMW
       WRITE( 6, 130     ) 'POP_KOA               : ', &
                            Input_Opt%POP_KOA
       WRITE( 6, 130     ) 'POP_KBC               : ', &
                            Input_Opt%POP_KBC
       WRITE( 6, 130     ) 'POP_K_POPG_OH         : ', &
                            Input_Opt%POP_K_POPG_OH
       WRITE( 6, 130     ) 'POP_K_POPP_O3A        : ', &
                            Input_Opt%POP_K_POPP_O3A
       WRITE( 6, 130     ) 'POP_K_POPP_O3B        : ', &
                            Input_Opt%POP_K_POPP_O3B
       WRITE( 6, 130     ) 'POP_HSTAR             : ', &
                            Input_Opt%POP_HSTAR
       WRITE( 6, 130     ) 'POP_DEL_H             : ', &
                            Input_Opt%POP_DEL_H
       WRITE( 6, 130     ) 'POP_DEL_Hw            : ', &
                            Input_Opt%POP_DEL_Hw
    ENDIF

    ! FORMAT statements
110 FORMAT( A, L5  )
120 FORMAT( A, A   )
130 FORMAT( A, ES10.2 )

    ! Return success
    RC = GC_SUCCESS

  END SUBROUTINE READ_POPS_MENU
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_passive_species_menu
!
! !DESCRIPTION: Subroutine READ\_PASSIVE\_SPECIES\_MENU reads the passive
!  species menu section of the GEOS-Chem input file; this defines passive
!  species to be used for this simulation.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE READ_PASSIVE_SPECIES_MENU( Input_Opt, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE Input_Opt_Mod, ONLY : OptInput
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(INOUT) :: Input_Opt   ! Input options
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!
! !REVISION HISTORY:
!  04 Sep 2015 - C. Keller   - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER            :: N, P, D

    ! Strings
    CHARACTER(LEN=255) :: ErrMsg, ThisLoc

    ! Arrays
    CHARACTER(LEN=255) :: SUBSTRS(MAXDIM)

    !=================================================================
    ! READ_PASSIVE_SPECIES_MENU begins here!
    !=================================================================

    ! Initialize
    RC      = GC_SUCCESS
    ErrMsg  = 'Error reading the "input.geos" file!'
    ThisLoc = ' -> at Read_Passive_Species_Menu (in GeosCore/input_mod.F90)'

    ! Initialize
    P = 0
    D = 0

    IF ( Input_Opt%amIRoot ) THEN
       WRITE( 6, '(/,a)' ) 'PASSIVE SPECIES MENU'
       WRITE( 6, '(  a)' ) '---------------------'
    ENDIF

    ! Do until exit
    DO

       ! Read passive species information for each passive species
       ! Every passive species line is expected to have 4 entries:
       ! - Species name
       ! - Species molecular weight
       ! - Atmospheric lifetime (s)
       ! - Initial atmospheric concentration (v/v)
       CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'Passive Species', RC )
       IF ( RC /= GC_SUCCESS ) THEN
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       ! Exit here if separator line is encountered
       IF ( INDEX ( TRIM(SUBSTRS(1)), '-----' ) > 0 ) EXIT

       ! Make sure there are at least 4 entries
       IF ( N < 4 ) THEN
          ErrMsg = 'Each passive species is expected to have '      // &
                   'at least four entries: Name, MW, TAU, initial ' // &
                   'concentration, and long name (optional).'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       ! Increase number of passive species by one
       P = P + 1

       ! Stop simulation and print warning if we exceed maximum number
       ! of passive species hardcoded in input_opt_mod.F90.
       IF ( P > Input_Opt%Max_PassiveSpc ) THEN
          ErrMsg = 'Number of passive species exceeds maximum. ' // &
                   'This value can be modified in input_opt_mod.F90.'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       ! Read and store species information
       Input_Opt%PASSIVE_ID(P)   = P
       Input_Opt%PASSIVE_NAME(P) = TRIM( SUBSTRS(1) )
       READ( SUBSTRS(2), * ) Input_Opt%PASSIVE_MW(P)
       READ( SUBSTRS(3), * ) Input_Opt%PASSIVE_TAU(P)
       READ( SUBSTRS(4), * ) Input_Opt%PASSIVE_INITCONC(P)

       ! Check if optional entry for long name is included
       IF ( N > 4 ) THEN
          READ( SUBSTRS(5), * ) Input_Opt%PASSIVE_LONGNAME(P)
       ENDIF

       ! Determine if the passive species decays (i.e. if it has
       ! an atmospheric lifetime that is not -1).  This will allow
       ! us to skip those passive species that do not decay in
       ! routine CHEM_PASSIVE_SPECIES, to speed up execution.
       IF ( Input_Opt%PASSIVE_TAU(P) > 0.0_fp ) THEN
          D                            = D + 1
          Input_Opt%NPASSIVE_DECAY     = D
          Input_Opt%PASSIVE_DECAYID(D) = P
       ENDIF

       ! Verbose
       IF ( Input_Opt%amIRoot ) THEN
          WRITE( 6, '(a)' ) 'Added passive species: '
          WRITE( 6, 110   ) ' - Species name                : ', &
                TRIM( Input_Opt%PASSIVE_NAME(P) )
          WRITE( 6, 120   ) ' - Molec. weight [g/mol]       : ', &
                Input_Opt%PASSIVE_MW(P)
          WRITE( 6, 130   ) ' - Lifetime [s]                : ', &
                Input_Opt%PASSIVE_TAU(P)
          WRITE( 6, 130   ) ' - Initial concentration [v/v] : ', &
                Input_Opt%PASSIVE_INITCONC(P)
          WRITE( 6, 110   ) ' - Species long name           : ', &
                TRIM( Input_Opt%PASSIVE_LONGNAME(P) )
       ENDIF

    ENDDO

    IF ( P < 0 ) THEN
       ErrMsg = 'Cannot add passive species '     // &
                TRIM(Input_Opt%PASSIVE_NAME(P) ) // &
                ': # of passive species is smaller than 1!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    ! Total number of passive species
    Input_Opt%NPASSIVE = P

110 FORMAT( A, A )
120 FORMAT( A, F10.2  )
130 FORMAT( A, ES10.2 )

    ! Return success
    RC = GC_SUCCESS

  END SUBROUTINE READ_PASSIVE_SPECIES_MENU
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: validate_directories
!
! !DESCRIPTION: Subroutine VALIDATE\_DIRECTORIES makes sure that each of the
!  directories that we have read from the GEOS-Chem input file are valid.
!  Also, trailing separator characters will be added.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE VALIDATE_DIRECTORIES( Input_Opt, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE Input_Opt_Mod, ONLY : OptInput
    USE TIME_MOD,      ONLY : EXPAND_DATE
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(INOUT) :: Input_Opt   ! Input Options object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!
! !REVISION HISTORY:
!  20 Jul 2004 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Strings
    CHARACTER(LEN=255) :: ErrMsg, ThisLoc, Dir

    !=================================================================
    ! VALIDATE_DIRECTORIES begins here!
    !=================================================================

    ! Initialize
    RC      = GC_SUCCESS
    ErrMsg  = 'Invalid directory encountered!'
    ThisLoc = ' -> at Validate_Directories (in module GeosCore/input_mod.F90)'

    ! Skip for dry-runs
    IF ( Input_Opt%DryRun ) RETURN

    ! Check directories
    CALL CHECK_DIRECTORY( Input_Opt, Input_Opt%DATA_DIR, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    CALL CHECK_DIRECTORY( Input_Opt, Input_Opt%CHEM_INPUTS_DIR, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    CALL CHECK_DIRECTORY( Input_Opt, Input_Opt%RUN_DIR, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

  END SUBROUTINE VALIDATE_DIRECTORIES
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: check_directory
!
! !DESCRIPTION: Subroutine CHECK\_DIRECTORY makes sure that the given
!  directory is valid.  Also a trailing slash character will be added if
!  necessary.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE CHECK_DIRECTORY( Input_Opt, DIR, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE FILE_MOD,      ONLY : File_Exists
    USE Input_Opt_Mod, ONLY : OptInput
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(OptInput),   INTENT(INOUT) :: Input_Opt   ! Input Options object
    CHARACTER(LEN=*), INTENT(INOUT) :: DIR         ! Dir to be checked
!
! !OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(OUT)   :: RC          ! Success or failure
!
! !REVISION HISTORY:
!  20 Mar 2003 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER            :: C

    ! Strings
    CHARACTER(LEN=255) :: ErrMsg, ThisLoc

    !=================================================================
    ! CHECK_DIRECTORY begins here!
    !=================================================================

    ! Initialize
    RC      = GC_SUCCESS
    ErrMsg  = ''
    ThisLoc = ' -> at Check_Directory (in module GeosCore/input_mod.F90)'

    ! Locate the last non-white-space character of NEWDIR
    C = LEN_TRIM( DIR )

    ! Add the trailing directory separator if it is not present
    IF ( DIR(C:C) /= '/' ) THEN
       DIR(C+1:C+1) = TRIM( '/' )
    ENDIF

    !=================================================================
    ! Test if the directory actually exists
    !=================================================================

    ! If the directory does not exist then stop w/ an error message
    IF ( .not. FILE_EXISTS( DIR ) ) THEN
       ErrMsg = 'Invalid directory: ' // TRIM( DIR )
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

  END SUBROUTINE CHECK_DIRECTORY
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: check_time_steps
!
! !DESCRIPTION: Subroutine CHECK\_TIME\_STEPS computes the smallest dynamic
!  time step for the model, based on which operation are turned on.  This
!  is called from routine READ\_INPUT\_FILE, after all of the timesteps and
!  logical flags have been read from "input.geos".
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE CHECK_TIME_STEPS( Input_Opt, State_Grid, RC)
!
! !USES:
!
    USE ErrCode_Mod
    USE Input_Opt_Mod,  ONLY : OptInput
    USE State_Grid_Mod, ONLY : GrdState
    USE TIME_MOD,       ONLY : SET_TIMESTEPS
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)  :: Input_Opt   ! Input Options object
    TYPE(GrdState), INTENT(IN)  :: State_Grid  ! Grid State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT) :: RC          ! Success or failure?
!
! !REVISION HISTORY:
!  20 Jul 2004 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    LOGICAL            :: LCONV, LCHEM,       LDRYD
    LOGICAL            :: LEMIS, LTRAN,       LTURB
    INTEGER            :: I,     J,           K
    INTEGER            :: L,     TS_SMALLEST, TS_DIAG
    INTEGER            :: TS_CHEM, TS_EMIS, TS_CONV, TS_DYN
    INTEGER            :: TS_UNIT, TS_RAD,  MAX_DYN

    ! Strings
    CHARACTER(LEN=255) :: ErrMsg, ThisLoc

    !=================================================================
    ! CHECK_TIME_STEPS begins here!
    !=================================================================

    ! Initialize
    RC      = GC_SUCCESS
    ErrMsg  = ''
    ThisLoc = ' -> at Check_Time_Steps (in module GeosCore/input_mod.F90)'

    ! Copy fields from Input_Opt
    LCONV = Input_Opt%LCONV
    LCHEM = Input_Opt%LCHEM
    LDRYD = Input_Opt%LDRYD
    LEMIS = Input_Opt%LEMIS
    LTRAN = Input_Opt%LTRAN
    LTURB = Input_Opt%LTURB

    TS_CHEM = Input_Opt%TS_CHEM
    TS_EMIS = Input_Opt%TS_EMIS
    TS_CONV = Input_Opt%TS_CONV
    TS_DYN  = Input_Opt%TS_DYN
    TS_RAD  = Input_Opt%TS_RAD

    ! NUNIT is time step in minutes for unit conversion
    TS_UNIT = -1

    ! Define maximum timestep for transport
    IF ( TRIM(State_Grid%GridRes) == '4.0x5.0') THEN
       MAX_DYN = 1800
    ELSE IF ( TRIM(State_Grid%GridRes) == '2.0x2.5' ) THEN
       MAX_DYN = 900
    ELSE IF ( TRIM(State_Grid%GridRes) == '0.5x0.625' ) THEN
       MAX_DYN = 600
    ELSE IF ( TRIM(State_Grid%GridRes) == '0.25x0.3125' ) THEN
       MAX_DYN = 300
    ELSE
       MAX_DYN = 3600
    ENDIF

    ! If TS_DYN is greater than MAX_DYN, then stop w/ error
    IF ( .not. Input_Opt%isMPI ) THEN
       IF ( Input_Opt%TS_DYN > MAX_DYN .and. LTRAN ) THEN
          WRITE( ErrMsg, 300 ) 'Transport timestep exceeds max:', &
                                Input_Opt%TS_DYN, MAX_DYN
300       FORMAT( a, i8, ' >', i8 )
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
    ENDIF

    ! Only do unit conversion if necessary
    IF ( LTRAN .or. LCONV .or. LTURB ) THEN
       TS_UNIT = MAX( TS_DYN, TS_CONV )
    ENDIF

    ! Compute NSMALLEST as the minimum of NDYN, NCONV, NSRCE, NCHEM
    I = TS_DYN
    J = TS_CONV
    K = TS_EMIS
    L = TS_CHEM

    ! SDE 2017-02-24: Always use LTRAN on the assumption that it will
    ! be used as a "heartbeat". This ensures that chemistry always
    ! takes place at the same time, regardless of whether or not
    ! transport is enabled.
    !IF ( .not. LTRAN                  ) I = 999999
    IF ( .not. LCONV .and..not. LTURB ) J = 999999
    IF ( .not. LDRYD .and..not. LEMIS ) K = 999999
    IF ( .not. LCHEM                  ) L = 999999

    ! Get the smallest of all of the above
    TS_SMALLEST = MIN( I, J, K, L )

    ! If all of the operators above are turned off,
    ! then set TS_SMALLEST to TS_DYN.
    IF ( TS_SMALLEST == 999999 ) THEN
       TS_SMALLEST = TS_DYN
    ENDIF

    IF ( LTRAN .and. TS_DYN /= TS_SMALLEST ) THEN
       ErrMsg = 'The transport time step should be the smallest one'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    ! If TS_DYN is smaller than TS_SMALLEST, reset TS_DYN
    ! to TS_SMALLEST.
    ! This is useful for runs where transport is turned off,
    ! but where chemistry is turned on.
    IF ( TS_DYN < TS_SMALLEST ) THEN
       TS_DYN = TS_SMALLEST
    ENDIF

    ! Define the largest time step, TS_DIAG, for diagnostics.
    ! Diagnostics should be incremented at the end of multiples of
    ! TS_DIAG, so that the system is at a physical state.
    ! (ccc, 5/13/09)
    IF ( .not. LTRAN                  ) I = -999999
    IF ( .not. LCONV .and..not. LTURB ) J = -999999
    IF ( .not. LDRYD .and..not. LEMIS ) K = -999999
    IF ( .not. LCHEM                  ) L = -999999

    TS_DIAG = MAX( I, J, K, L )

    ! If all the operators are turned off, then set TS_DIAG to TS_CHEM
    ! Usually the chemistry time step is large. (ccc, 5/13/09)
    IF ( TS_DIAG == -999999 ) THEN
       TS_DIAG = TS_CHEM
    ENDIF

    ! Check if all time steps are multiples of the smallest.
    ! (ccc, 5/13/09)
    IF ( L /= -999999 .and. MOD( TS_CHEM, TS_SMALLEST ) /= 0 ) THEN
       WRITE( ErrMsg, 100 ) 'Chemistry', TS_CHEM, TS_SMALLEST
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    IF ( K /= -999999 .and. MOD( TS_EMIS, TS_SMALLEST ) /= 0 ) THEN
       WRITE( ErrMSg, 100 ) 'Emission', TS_EMIS, TS_SMALLEST
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    IF ( J /= -999999 .and. MOD( TS_CONV, TS_SMALLEST ) /= 0 ) THEN
       WRITE( ErrMsg, 100 ) 'Convection', TS_CONV, TS_SMALLEST
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    IF ( I /= -999999 .and. MOD( TS_DYN, TS_SMALLEST ) /= 0 ) THEN
       WRITE( ErrMsg, 100 ) 'Transport', TS_DYN, TS_SMALLEST
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    ! Initialize timesteps in "time_mod.F90"
    CALL SET_TIMESTEPS( Input_Opt,            &
                        CHEMISTRY  = TS_CHEM, &
                        EMISSION   = TS_EMIS, &
                        DYNAMICS   = TS_DYN,  &
                        UNIT_CONV  = TS_UNIT, &
                        CONVECTION = TS_CONV, &
                        DIAGNOS    = TS_DIAG, &
                        RADIATION  = TS_RAD )

100 FORMAT( A, ' time step must be a multiple of the smallest one:', i5, i5 )

  END SUBROUTINE CHECK_TIME_STEPS
!EOC
#if !(defined( EXTERNAL_GRID ) || defined( EXTERNAL_FORCING ))
#ifdef BPCH_DIAG
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: is_last_day_good
!
! !DESCRIPTION: Suborutine IS\_LAST\_DAY\_GOOD tests to see if there is
!  output scheduled on the last day of the run.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE IS_LAST_DAY_GOOD( Input_Opt, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE Input_Opt_Mod, ONLY : OptInput
    USE JULDAY_MOD,    ONLY : JULDAY
    USE TIME_MOD,      ONLY : GET_NYMDe, ITS_A_LEAPYEAR, YMD_EXTRACT
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(INOUT) :: Input_Opt  ! Input options
    INTEGER,        INTENT(OUT)   :: RC         ! Success or failure?
!
! !REVISION HISTORY:
!  20 Jul 2004 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    LOGICAL            :: IS_LEAPYEAR
    INTEGER            :: NYMDe, Y, M, D, LASTDAY
    REAL(fp)           :: JD, JD0

    ! Strings
    CHARACTER(LEN=255) :: ErrMsg, ThisLoc

    !=================================================================
    ! IS_LAST_DAY_GOOD begins here!
    !=================================================================

    ! Initialize
    RC      = GC_SUCCESS
    ErrMsg  = ''
    ThisLoc = ' -> at Is_Last_Day_Good (in module GeosCore/input_mod.F90)'

    ! Astronomical Julian Day corresponding to NYMDe
    NYMDe = GET_NYMDe()
    CALL YMD_EXTRACT( NYMDe, Y, M, D )
    JD = JULDAY( Y, M, DBLE( D ) )

    ! Astronomical Julian Day corresponding to the 1st of the year
    JD0 = JULDAY( Y, 1, 0d0 )

    ! LASTDAY is the day of year corresponding to NYMDe
    LASTDAY = JD - JD0

    ! Skip past the element of NJDAY for Feb 29, if necessary
    IF ( .not. ITS_A_LEAPYEAR( Y, .TRUE. ) .and. LASTDAY > 59 ) THEN
       LASTDAY = LASTDAY + 1
    ENDIF

    ! Exit w/ error if THIS_NJDAY = 0
    IF ( Input_Opt%NJDAY(LASTDAY) == 0 ) THEN
       ErrMsg = 'No output scheduled on last day of run!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

  END SUBROUTINE IS_LAST_DAY_GOOD
!EOC
#endif
#endif
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: do_error_checks
!
! !DESCRIPTION: Makes sure that certain species are defined in order to
!  proceed with a certain option.  Halts the simulation with an error message
!  if incorrect inputs  would have caused  a simulation to crash.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Do_Error_Checks( Input_Opt, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE Input_Opt_Mod, ONLY : OptInput
    USE State_Chm_Mod, ONLY : Ind_
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(INOUT) :: Input_Opt
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC
!
! !REMARKS:
!  These error checks were originally called when the various menus were
!  read in from disk.  However, in order to use the Ind_() function to look
!  up given species indices, we need to call this after the Species Database
!  (which is in State_Chm) is initialized.  Therefore, we have now moved these
!  error checks to this routine, which is now called from GC_Init_Extra.
!                                                                             .
!  The Ind_() function now defines all species ID's.  It returns -1 if
!  a species cannot be found.  The prior behavior was to return 0 if a
!  species wasn't found.  Therefore, in order to preserve the logic of the
!  error checks, we must force any -1's returned by Ind_() to 0's in
!  this subroutine.
!
! !REVISION HISTORY:
!  22 Jun 2016 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER            :: I

    ! Strings
    CHARACTER(LEN=255) :: MSG, LOCATION

    !=================================================================
    ! Initialization
    !=================================================================

    ! Assume success
    RC       = GC_SUCCESS

    ! Define location string
    LOCATION = '-> at Do_Error_Checks (in GeosCore/input_mod.F90)'

    !=================================================================
    ! Error check SEASALT AEROSOLS
    !=================================================================
    I = MAX( Ind_('SALA','A'), 0 ) + MAX( Ind_('SALC','A'), 0 )

    IF ( Input_Opt%LSSALT ) THEN
       IF ( I == 0 ) THEN
          MSG = 'LSSALT=T but ONLINE SEASALT AEROSOLS are undefined!'
          CALL GC_Error( Msg, RC, Location )
          RETURN
       ENDIF
    ELSE
       IF ( I > 0 ) THEN
          MSG = 'Cannot use ONLINE SEASALT AEROSOLS if LSSALT=F!'
          CALL GC_Error( Msg, RC, Location )
          RETURN
       ENDIF
    ENDIF

    !=================================================================
    ! Error check MARINE ORGANIC AEROSOLS
    !=================================================================
    I = MAX( Ind_('MOPO','A'), 0 ) + MAX( Ind_('MOPI','A'), 0 )

    IF ( Input_Opt%LMPOA ) THEN
       IF ( .not. Input_Opt%LSSALT ) THEN
          MSG = 'LMPOA=T but LSSALT=F!'
          CALL GC_Error( Msg, RC, Location )
          RETURN
       ENDIF
       IF ( I == 0 ) THEN
          MSG = 'LMPOA=T but MARINE ORGANIC AEROSOLS are undefined!'
          CALL GC_Error( Msg, RC, Location )
          RETURN
       ENDIF
    ELSE
       IF ( I > 0 ) THEN
          MSG = 'Cannot use MARINE ORGANIC AEROSOLS if LMPOA=F!'
          CALL GC_Error( Msg, RC, Location )
          RETURN
       ENDIF
    ENDIF

    !=================================================================
    ! Error check SULFUR AEROSOLS
    !=================================================================
    I = MAX( Ind_('DMS' ,'A'), 0 ) + &
        MAX( Ind_('SO2' ,'A'), 0 ) + &
        MAX( Ind_('SO4' ,'A'), 0 ) + &
        MAX( Ind_('SO4s','A'), 0 ) + &
        MAX( Ind_('MSA' ,'A'), 0 ) + &
        MAX( Ind_('NH3' ,'A'), 0 ) + &
        MAX( Ind_('NH4' ,'A'), 0 ) + &
        MAX( Ind_('NITs','A'), 0 )

    IF ( Input_Opt%LSULF ) THEN

       ! We now compute the production of SO4s and NITs, so when
       ! LSULF=T, then we must also have LSSALT=T (bec, bmy, 4/13/05)
       IF ( .not. Input_Opt%LSSALT ) THEN
          MSG = 'LSULF=T now also requires LSSALT=T!'
          CALL GC_Error( Msg, RC, Location )
          RETURN
       ENDIF

       ! Stop w/ error if everything is undefined
       IF ( I == 0 ) THEN
          MSG = 'LSULF=T but ONLINE SULFUR AEROSOLS are undefined!'
          CALL GC_Error( Msg, RC, Location )
          RETURN
       ENDIF

    ELSE

       ! If LSULF=F but we have defined species, stop w/ error
       IF ( I > 0 ) THEN
          MSG = 'Cannot use ONLINE SULFUR AEROSOLS if LSULF=F!'
          CALL GC_Error( Msg, RC, Location )
          RETURN
       ENDIF

    ENDIF

    !=================================================================
    ! Error check CARBON AEROSOLS
    !=================================================================

    ! SOAupdate: Add POA (hotp 10/11/09)
    I = MAX( Ind_('BCPO','A'), 0 ) + &
        MAX( Ind_('BCPI','A'), 0 ) + &
        MAX( Ind_('OCPO','A'), 0 ) + &
        MAX( Ind_('OCPI','A'), 0 ) + &
        MAX( Ind_('POA1','A'), 0 )

    IF ( Input_Opt%LCARB ) THEN
       IF ( I == 0 ) THEN
          MSG = 'LCARB=T but ONLINE CARBON AEROSOLS are undefined!'
          CALL GC_Error( Msg, RC, Location )
          RETURN
       ENDIF
    ELSE
       IF ( I > 0 ) THEN
          MSG = 'Cannot use ONLINE CARBON AEROSOLS if LCARB=F!'
          CALL GC_Error( Msg, RC, Location )
          RETURN
       ENDIF
    ENDIF

    IF ( Input_Opt%LSVPOA .and. ( .NOT. Input_Opt%LSOA ) ) THEN
       MSG = 'Semivolatile POA requires COMPLEX SOA (LSOA=T)'
       CALL GC_Error( Msg, RC, Location )
       RETURN
    ENDIF

    ! SOAupdate: Error check (hotp 8/24/09)
    ! OCPI and OCPO are the non-volatile POA species
    ! POA (along w/ POG, OPOA, and OPOG) are the semivol POA species
    ! You can't have both!
    I = MAX( Ind_('OCPI','A'), 0 ) + MAX( Ind_('OCPO','A'), 0 )

    IF ( Ind_('POA1') > 0 ) THEN
       IF ( I > 0 ) THEN
          MSG = 'Semivolatile POA species is defined in addition to ' // &
                 'Nonvolatile POA'
          CALL GC_Error( Msg, RC, Location )
          RETURN
       ENDIF
       IF ( ( .NOT. Input_Opt%LSOA   ) .or. &
            ( .NOT. Input_Opt%LSVPOA ) ) THEN
          MSG = 'Semivolatile POA requires LSOA=T and LSVPOA=T'
          CALL GC_Error( Msg, RC, Location )
          RETURN
       ENDIF
    ENDIF

    ! SOAupdate
    ! Options for organic aerosol species:
    ! IF LSOA = F: only OCPI and OCPO
    ! IF LSOA = T:
    !   OCPI OCPO SOA (non-vol + original traditional)
    !   POA POG OPOA OPOG SOA BTX NAP (semivol + orig trad + IVOC )
    ! NAP emissions are set in input.geos
    ! LSVPOA is just a check (doesn't do anything hotp 7/21/10)
    I = MAX( Ind_('POA1' ,'A'), 0 ) + &
        MAX( Ind_('POA2' ,'A'), 0 ) + &
        MAX( Ind_('POG1' ,'A'), 0 ) + &
        MAX( Ind_('POG2' ,'A'), 0 ) + &
        MAX( Ind_('OPOA1','A'), 0 ) + &
        MAX( Ind_('OPOA2','A'), 0 ) + &
        MAX( Ind_('OPOG1','A'), 0 ) + &
        MAX( Ind_('OPOG2','A'), 0 )

    IF ( Input_Opt%LSVPOA ) THEN
       IF ( I < 8 ) THEN
          MSG = 'Not enough semivolatile POA species!'
          CALL GC_Error( Msg, RC, Location )
          RETURN
       ENDIF
       IF ( Ind_('NAP','A') < 0 ) THEN
          MSG = 'Semivolatile POA requires IVOCs/NAP!'
          CALL GC_Error( Msg, RC, Location )
          RETURN
       ENDIF
    ENDIF

    !=================================================================
    ! Error check SECONDARY ORGANIC AEROSOLS
    !=================================================================

    ! Check for complex SOA species
    I = MAX( Ind_('TSOA1','A'), 0 ) + &
        MAX( Ind_('TSOA2','A'), 0 ) + &
        MAX( Ind_('TSOA3','A'), 0 ) + &
        MAX( Ind_('ASOA1','A'), 0 ) + &
        MAX( Ind_('ASOA2','A'), 0 ) + &
        MAX( Ind_('ASOA3','A'), 0 ) + &
        MAX( Ind_('ASOAN','A'), 0 )

    IF ( Input_Opt%LSOA ) THEN
       IF ( I == 0 ) THEN
          MSG = 'LSOA=T but COMPLEX SOA species are undefined!'
          CALL GC_Error( Msg, RC, Location )
          RETURN
       ENDIF
    ELSE
       IF ( I > 0 ) THEN
          MSG = 'Cannot use COMPLEX SOA species if LSOA=F!'
          CALL GC_Error( Msg, RC, Location )
          RETURN
       ENDIF
    ENDIF

    !=================================================================
    ! Error check DUST AEROSOLS
    !=================================================================

#ifdef TOMAS
    ! For TOMAS only: If DUST1 is present, the other dust species are too
    I = MAX( Ind_('DUST1','A'), 0 )
#else
    ! Non-TOMAS simulations: Need all DST1-DST4 species
    I = MAX( Ind_('DST1','A'), 0 ) + &
        MAX( Ind_('DST2','A'), 0 ) + &
        MAX( Ind_('DST3','A'), 0 ) + &
        MAX( Ind_('DST4','A'), 0 )
#endif

    IF ( Input_Opt%LDUST ) THEN
       IF ( I == 0 ) THEN
          MSG = 'LDUST=T but ONLINE DUST AEROSOLS are undefined!'
          CALL GC_Error( Msg, RC, Location )
          RETURN
       ENDIF
    ELSE
       IF ( I > 0 ) THEN
          MSG = 'Cannot use ONLINE DUST AEROSOLS if LDUST=F!'
          CALL GC_Error( Msg, RC, Location )
          RETURN
       ENDIF
    ENDIF

    !=================================================================
    ! Error check DUST NITRATE    AEROSOLS
    !             DUST SULFATE    AEROSOLS
    !             DUST ALKALINITY AEROSOLS
    !=================================================================
    I = MAX( Ind_('NITd1'  ,'A'), 0 ) + &
        MAX( Ind_('NITd2'  ,'A'), 0 ) + &
        MAX( Ind_('NITd3'  ,'A'), 0 ) + &
        MAX( Ind_('NITd4'  ,'A'), 0 ) + &
        MAX( Ind_('SO4d1'  ,'A'), 0 ) + &
        MAX( Ind_('SO4d2'  ,'A'), 0 ) + &
        MAX( Ind_('SO4d3'  ,'A'), 0 ) + &
        MAX( Ind_('SO4d4'  ,'A'), 0 ) + &
        MAX( Ind_('DSTAL1' ,'A'), 0 ) + &
        MAX( Ind_('DSTAL2' ,'A'), 0 ) + &
        MAX( Ind_('DSTAL3' ,'A'), 0 ) + &
        MAX( Ind_('DSTAL4' ,'A'), 0 )

    IF ( Input_Opt%LDSTUP ) THEN
       IF ( I < 12 ) THEN
          MSG = 'LDSTUP=T but COATED DUST AEROSOLS are undefined!'
          CALL GC_Error( Msg, RC, Location )
          RETURN
       ENDIF
    ELSE
       IF ( I > 0 ) THEN
          MSG = 'Cannot use COATED DUST AEROSOLS if LDSTUP=F!'
          CALL GC_Error( Msg, RC, Location )
          RETURN
       ENDIF
    ENDIF

    !=================================================================
    ! Error check SEASALT AEROSOLS
    !=================================================================
    I = MAX( Ind_('SALA','A'), 0 ) + MAX( Ind_('SALC','A'), 0 )

    IF ( Input_Opt%LSSALT ) THEN
       IF ( I == 0 ) THEN
          MSG = 'LSSALT=T but ONLINE SEASALT AEROSOLS are undefined!'
          CALL GC_Error( Msg, RC, Location )
          RETURN
       ENDIF
    ELSE
       IF ( I > 0 ) THEN
          MSG = 'Cannot use ONLINE SEASALT AEROSOLS if LSSALT=F!'
          CALL GC_Error( Msg, RC, Location )
          RETURN
       ENDIF
    ENDIF

    !=================================================================
    ! Error check stratospheric H2O
    !=================================================================
    IF ( Input_Opt%LUCX .and. Input_Opt%LSETH2O ) THEN
       IF (Ind_('H2O') < 0 ) THEN
          WRITE( 6, '(a)'     ) REPEAT( '=', 79 )
          WRITE( 6, '(/,a,/)' ) 'Warning in input_mod.F90: ' &
               // 'H2O is set but H2O species is undefined.'
          Input_Opt%LSETH2O = .FALSE.
          WRITE( 6, '(a)'     ) REPEAT( '=', 79 )
       ENDIF
    ELSE
       Input_Opt%LSETH2O = .FALSE.
    ENDIF

  END SUBROUTINE DO_ERROR_CHECKS
!EOC
END MODULE INPUT_MOD
