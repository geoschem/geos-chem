!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: planeflight_mod.F90
!
! !DESCRIPTION: Module PLANEFLIGHT\_MOD contains variables and routines which
!  are used to "fly" a plane through the GEOS-Chem model simulation.  This is
!  useful for comparing model results with aircraft observations.
!\\
!\\
! !INTERFACE:
!
MODULE PLANEFLIGHT_MOD
!
! !USES:
!
  USE inquireMod,    ONLY : findFreeLUN
  USE PRECISION_MOD       ! For GEOS-Chem Precision (fp)

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Prior to 5/12/16:
! This routine does not work with the FlexChem implementation. We need to
! rewrite this code to get the necessary information from KPP (mps, 5/12/16)
!  PUBLIC  :: ARCHIVE_RXNS_FOR_PF
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  PUBLIC  :: CLEANUP_PLANEFLIGHT
  PUBLIC  :: PLANEFLIGHT
  PUBLIC  :: SETUP_PLANEFLIGHT
  PUBLIC  :: SET_PLANEFLIGHT
!
! !PRIVATE MEMBER FUNCTIONS:
!
  PRIVATE :: AN_SETUP
  PRIVATE :: INIT_PLANEFLIGHT
  PRIVATE :: NOY_SETUP
  PRIVATE :: READ_VARIABLES
  PRIVATE :: READ_POINTS
  PRIVATE :: RO2_SETUP
  PRIVATE :: TEST_VALID
  PRIVATE :: WRITE_VARS_TO_FILE
!
! !REMARKS:
!  The quantities that are saved to disk by the planeflight diagnostic were
!  requested by GEOS-Chem users.  If you would like to save out a new quantity,
!  then you will have to make your own modifications in this module.
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !PRIVATE TYPES:
!
  !=================================================================
  ! MODULE VARIABLES:
  !
  ! DO_PF       : Turn on the planeflight diagnostic? (T/F)
  ! MAXVARS     : Maximum # of variables allowed
  ! MAXPOINTS   : Maximum # of flight track points allowed
  ! MAXREAC     : Maximum # of reactions allowed
  ! MAXRO2      : Maximum # of RO2 constituents allowed
  ! NPOINTS     : Number of flight track points
  ! PPOINT      : Pointer to last measured output
  ! PDATE       : Array of dates     at each flight point
  ! PTIME       : Array of times     at each flight point
  ! PTAU        : Array of TAU's     at each flight point
  ! PLAT        : Array of latitude  at each flight point
  ! PLON        : Array of longitude at each flight point
  ! PPRESS      : Array of pressure  at each flight point
  ! PTYPE       : Array of ID'#S     at each flight point
  ! NPVAR       : # of var's to be saved at each flight point
  ! PVAR        : Array of variable indices
  ! PNAME       : Array of variable names corresponding to PVAR
  ! NPREAC      : # of variables that are really rxns
  ! PREAC       : Array of rxn index numbers
  ! PRRATE      : Array of rxn rates for each entry in PREAC
  ! NRO2        : # number of RO2 constituents
  ! PRO2        : Array of species that are RO2 const's
  ! INFILENAME  : Name of input file defining the flight track
  ! OUTFILENAME : Name of output file
  !=================================================================

  ! Logicals
  LOGICAL                        :: DO_PF

  ! Parameters
  INTEGER,           PARAMETER   :: MAXVARS   = 800
  INTEGER,           PARAMETER   :: MAXPOINTS = 100000
  INTEGER,           PARAMETER   :: MAXREAC   = 50
  INTEGER,           PARAMETER   :: MAXRO2    = 45
  INTEGER,           PARAMETER   :: MAXAN     = 15
  INTEGER,           PARAMETER   :: MAXNOY    = 15

  ! For specifying flight track points
  INTEGER                        :: NPOINTS
  INTEGER                        :: PPOINT

  ! For specifying date/time
  INTEGER,           ALLOCATABLE :: PDATE(:)
  INTEGER,           ALLOCATABLE :: PTIME(:)
  REAL(fp),          ALLOCATABLE :: PTAU(:)

  ! For specifying lat/lon/alt and ID type
  REAL*4,            ALLOCATABLE :: PLAT(:)
  REAL*4,            ALLOCATABLE :: PLON(:)
  REAL*4,            ALLOCATABLE :: PPRESS(:)
  REAL*4,            ALLOCATABLE :: POBS(:)
  REAL*4,            ALLOCATABLE :: PTAMB(:)
  REAL*4,            ALLOCATABLE :: PH2OMR(:)
  REAL*4,            ALLOCATABLE :: PPOTTEMP(:)
  REAL*4,            ALLOCATABLE :: PRH(:)
  REAL*4,            ALLOCATABLE :: PGPSALT(:)
  CHARACTER(LEN=7),  ALLOCATABLE :: PTYPE(:)

  ! For specifying variables to save at each flight point
  INTEGER                        :: NPVAR
  INTEGER,           ALLOCATABLE :: PVAR(:)
  CHARACTER(LEN=10), ALLOCATABLE :: PNAME(:)

  ! For specifying rxns to save at each flight point
  INTEGER                        :: NPREAC
  INTEGER,           ALLOCATABLE :: PREAC(:)
  REAL(fp),          ALLOCATABLE :: PRRATE(:,:,:,:)

  ! For rate of production
  INTEGER                        :: NPROD
  INTEGER,           ALLOCATABLE :: IPROD(:)

  ! For specifying RO2 constituents at each flight point
  INTEGER                        :: NPRO2
  INTEGER                        :: PRO2(MAXRO2)

  ! For specifying NOY constituents at each flight point
  INTEGER                        :: NPNOY
  INTEGER                        :: PNOY(MAXNOY)

  ! For specifying AN constituents at each flight point
  INTEGER                        :: NPAN
  INTEGER                        :: P_AN(MAXAN)

  ! Input/output file names
  CHARACTER(LEN=255)             :: INFILENAME,  INF
  CHARACTER(LEN=255)             :: OUTFILENAME, OUTF

  ! Logical unit numbers
  INTEGER                        :: IU_FILE
  INTEGER                        :: IU_PLANE

CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: setup_planeflight
!
! !DESCRIPTION: Subroutine SETUP\_PLANEFLIGHT reads information from the
!  input file in order to initialize the planeflight diagnostic.  Also
!  calls INIT\_PLANEFLIGHT to allocate and zero module arrays.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE SETUP_PLANEFLIGHT( Input_Opt, State_Chm, State_Grid, &
                                State_Met, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE FILE_MOD,           ONLY : FILE_EXISTS
    USE FILE_MOD,           ONLY : IOERROR
    USE Input_Opt_Mod,      ONLY : OptInput
    USE State_Chm_Mod,      ONLY : ChmState
    USE State_Grid_Mod,     ONLY : GrdState
    USE State_Met_Mod,      ONLY : MetState
    USE TIME_MOD,           ONLY : EXPAND_DATE
    USE TIME_MOD,           ONLY : GET_NYMD
    USE TIME_MOD,           ONLY : GET_NHMS
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)  :: Input_Opt   ! Input Options object
    TYPE(ChmState), INTENT(IN)  :: State_Chm   ! Chemistry State object
    TYPE(GrdState), INTENT(IN)  :: State_Grid  ! Grid State object
    TYPE(MetState), INTENT(IN)  :: State_Met   ! Meteorology State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT) :: RC          ! Success or failure?
!
! !REMARKS:
!
! !REVISION HISTORY:
!  30 Jul 2002 - M. Evans    - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    LOGICAL, SAVE      :: FIRST = .TRUE.
    LOGICAL            :: IS_OPEN
    INTEGER            :: I,  IP,      N,   TEMP, LENGTH
    INTEGER            :: RN, COUNTER, IOS, NYMD, NHMS
    CHARACTER(LEN=6)   :: TYPE

    !=================================================================
    ! SETUP_PLANEFLIGHT begins here!
    !=================================================================

    ! Assume that there is flight data for today
    DO_PF   = .TRUE.

    ! Find a free file LUN
    IU_FILE = findFreeLun()

    ! Get date & time
    NYMD    = GET_NYMD()
    NHMS    = GET_NHMS()

    ! Copy file names to local variables
    INF     = INFILENAME
    OUTF    = OUTFILENAME

    ! Replace any date & time tokens in the file names
    CALL EXPAND_DATE( INF,  NYMD, NHMS )
    CALL EXPAND_DATE( OUTF, NYMD, NHMS )

    ! If we can't find a flighttrack file for today's date, return
    IF ( .not. FILE_EXISTS( INF ) ) THEN
       DO_PF = .FALSE.
       RETURN
    ENDIF

    ! Echo info
    IF ( Input_Opt%amIRoot ) THEN
       WRITE( 6, '(a)' ) REPEAT( '=', 79 )
       WRITE( 6, 99    )
       WRITE( 6, 100   ) TRIM( INF ), IU_FILE
99     FORMAT( 'P L A N E   F L I G H T   D I A G N O S T I C' )
100    FORMAT( /, 'SETUP_PLANEFLIGHT: Reading ',a, ' on unit ', i4 )
       WRITE( 6, '(a)' )
    ENDIF

    ! Compute # of species and # of points & allocate arrays
    CALL INIT_PLANEFLIGHT( Input_Opt, State_Grid )

    ! Return if there are no flight track points for today
    IF ( NPOINTS == 0 ) THEN
       IF ( Input_Opt%amIRoot ) THEN
          WRITE( 6, '(a)' ) 'No flight track found for today!'
       ENDIF
       DO_PF = .FALSE.
       RETURN
    ENDIF

    !=================================================================
    ! Open file and read info
    !=================================================================

    ! Open input file
    OPEN( IU_FILE, FILE=TRIM( INF ), IOSTAT=IOS )
    IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_FILE, 'setup_planeflight:1')

    ! Read variables to be output -- sort into PVAR array by type
    CALL READ_VARIABLES( Input_Opt, State_Chm, IU_FILE,   RC )

    ! Read information about each point (date/time/lon/lat/alt)
    CALL READ_POINTS( Input_Opt, State_Grid, State_Met, IU_FILE, RC)

    ! Close the file
    CLOSE( IU_FILE )

    ! Set the pointer to the first record
    PPOINT = 1

    !=================================================================
    ! Find the species # for all components of RO2 (fullchem only)
    !=================================================================
    CALL RO2_SETUP( Input_Opt, State_Chm, RC )

    !=================================================================
    ! Find the species # for all components of NOY (fullchem only)
    !=================================================================
    CALL NOY_SETUP( Input_Opt, State_Chm, RC )

    !=================================================================
    ! Find the species # for all components of AN (fullchem only)
    !=================================================================
    CALL AN_SETUP( Input_Opt, State_Chm, RC )

    ! Fancy output
    WRITE( 6, '(a)' ) REPEAT( '=', 79 )

    !=================================================================
    ! Open today's plane.log file and write file header
    !=================================================================

    ! Close previously opened plane.log file
    INQUIRE( IU_PLANE, OPENED=IS_OPEN )
    IF( IS_OPEN ) CLOSE( IU_PLANE )

    ! Close previously-opened file
    CLOSE( IU_PLANE )

    ! Find a free file LUN
    IU_PLANE = findFreeLUN()

    ! Open output file
    OPEN( IU_PLANE, FILE=TRIM( OUTF ), STATUS='UNKNOWN', IOSTAT=IOS )

    ! Error check
    IF ( IOS /= 0 ) THEN
       CALL IOERROR( IOS, IU_PLANE, 'setup_planeflight:1' )
    ENDIF

    ! Write header
    IF ( Input_Opt%amIRoot ) THEN

       WRITE( IU_PLANE, 110 ) 'POINT', 'TYPE',  'YYYYMMDD', 'HHMM', &
                              'LAT',   'LON',   'PRESS',    'OBS', &
                              'T-IND', 'P-IND', 'I-IND',    'J-IND', &
                              ( PNAME(I), I=1,NPVAR )
    ENDIF

    ! FORMAT string
#if   defined( TOMAS )
110 FORMAT( A5,X,A7,X,A8,X,A4,X,A7,X,A7,X,A7,X,A10,X, &
            A9,X,A3,X,A5,X,A5,X,250(a10,x) )
#else
110 FORMAT( A5,X,A7,X,A8,X,A4,X,A7,X,A7,X,A7,X,A10,X, &
            A9,X,A3,X,A5,X,A5,X,200(a11,x) )
#endif

  END SUBROUTINE SETUP_PLANEFLIGHT
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_variables
!
! !DESCRIPTION: Subroutine READ\_VARIABLES reads the list of variables
!  (chemical species, rxn rates, GMAO met fields, or GEOS-Chem species) to be
!  printed out and sorts the information into the appropriate module variables.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE READ_VARIABLES( Input_Opt, State_Chm, IU_FILE, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE ERROR_MOD,  ONLY : GEOS_CHEM_STOP
    USE FILE_MOD,   ONLY : IOERROR
    USE Input_Opt_Mod,      ONLY : OptInput
    USE Species_Mod,        ONLY : Species
    USE State_Chm_Mod,      ONLY : ChmState
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)  :: Input_Opt  ! Input Options object
    TYPE(ChmState), INTENT(IN)  :: State_Chm  ! Chemistry State object
    INTEGER,        INTENT(IN)  :: IU_FILE    ! Logical unit # for ASCII read
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT) :: RC         ! Success or failure?
!
! !REVISION HISTORY:
!  30 Jul 2002 - M. Evans    - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    LOGICAL             :: IS_FULLCHEM
    INTEGER             :: M, N, NUM, R, IK, IOS, nAdvect
    INTEGER             :: PR, J, NF, FM
    CHARACTER(LEN=255)  :: LINE
    CHARACTER(LEN=10)   :: PRODNAME

    ! Objects
    TYPE(Species), POINTER :: SpcInfo

    !=================================================================
    ! READ_VARIABLES begins here!
    !=================================================================

    ! Assume success
    RC = GC_SUCCESS

    ! Get fields from Input_Opt
    IS_FULLCHEM = Input_Opt%ITS_A_FULLCHEM_SIM

    ! Number of advected species
    nAdvect     = State_Chm%nAdvect

    ! Initialize pointer
    SpcInfo => NULL()

    ! Read four lines of header
    DO N = 1, 4
       READ( IU_FILE, '(a)', IOSTAT=IOS )
       IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_FILE, 'read_variables:1')
    ENDDO

    ! Read in the number of species to be output
    ! Read in as I4 now for the number of variables we save (skim, 7/24/13)
    READ( IU_FILE, '(i4)', IOSTAT=IOS ) NPVAR
    IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_FILE, 'read_variables:2' )

    ! Read in a separation line
    READ( IU_FILE, '(a)', IOSTAT=IOS )
    IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_FILE, 'read_variables:3' )

    ! Echo to stdout
    IF ( Input_Opt%amIRoot ) THEN
       WRITE( 6, '(a)' ) '   #    Species       PVAR'
       WRITE( 6, '(a)' ) '-----------------------------'
    ENDIF

    !=================================================================
    ! Sort variables by type; assign indices to PVAR, PREAC arrays
    ! NOTE: Variables for which PVAR(N) = 0 will be skipped!
    !=================================================================

    ! Zero reaction counter
    R = 0

    ! Zero production counter
    PR = 0

    ! Loop over all variables
    DO N = 1, NPVAR

       ! Read each line
       READ( IU_FILE, '(a)', IOSTAT=IOS ) LINE
       IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_FILE, 'read_variables:4')

       ! Save the name of each variable into the global PNAME array
       PNAME(N) = LINE(1:10)

       ! We are searching for a ...
       SELECT CASE ( LINE(1:4) )

       !===========================================================
       ! GEOS-CHEM tracer: listed as "TRA_001", etc.
       ! PVAR offset: 100000
       !===========================================================
       CASE ( 'TRA_' )

          ! Extract tracer # from the string
          READ( LINE(5:14), '(i10)' ) NUM

          ! Make sure the tracer # is valid!
          IF ( NUM < 0 .or. NUM > nAdvect ) THEN
             IF ( Input_Opt%amIRoot ) THEN
                WRITE( 6, 100 ) TRIM( LINE )
100             FORMAT( 'TRACER ', i4, ' is out of range!' )
                WRITE( 6, '(a)' ) 'STOP in SETUP_PLANEFLIGHT!'
                WRITE( 6, '(a)' ) REPEAT( '=', 79 )
             ENDIF
             CALL GEOS_CHEM_STOP
          ENDIF

          ! Save in PVAR -- add offset of 100000
          PVAR(N) = 100000 + NUM

#ifdef TOMAS
       ! Add case matching for TOMAS microphysics rates (win, 7/28/09)
       !===========================================================
       ! GEOS-CHEM tracer: listed as "TMS_001", etc.
       ! PVAR offset: 200000
       !===========================================================
       CASE ( 'TMS_' )

          ! Extract tracer # from the string
          READ( LINE(5:14), '(i10)' ) NUM

          ! Make sure the tracer # is valid!
          IF ( NUM < 0 .or. NUM > nAdvect ) THEN
             IF ( Input_Opt%amIRoot ) THEN
                WRITE( 6, 100 ) TRIM( LINE )
                WRITE( 6, '(a)' ) 'STOP in SETUP_PLANEFLIGHT!'
                WRITE( 6, '(a)' ) REPEAT( '=', 79 )
                CALL GEOS_CHEM_STOP
             ENDIF
          ENDIF

          ! Save in PVAR -- add offset of 200000
          PVAR(N) = 200000 + NUM
#endif

       !===========================================================
       ! GMAO met field: listed as "GMAO_TEMP", etc.
       ! PVAR offset: 1000
       !===========================================================
       CASE ( 'GMAO' )

          IF ( LINE == 'GMAO_TEMP'  ) PVAR(N) = 1001
          IF ( LINE == 'GMAO_ABSH'  ) PVAR(N) = 1002
          IF ( LINE == 'GMAO_SURF'  ) PVAR(N) = 1003
          IF ( LINE == 'GMAO_PSFC'  ) PVAR(N) = 1004
          IF ( LINE == 'GMAO_UWND'  ) PVAR(N) = 1005
          IF ( LINE == 'GMAO_VWND'  ) PVAR(N) = 1006
          IF ( LINE == 'GMAO_IIEV'  ) PVAR(N) = 1007
          IF ( LINE == 'GMAO_JJEV'  ) PVAR(N) = 1008
          IF ( LINE == 'GMAO_LLEV'  ) PVAR(N) = 1009
          IF ( LINE == 'GMAO_RELH'  ) PVAR(N) = 1010
          IF ( LINE == 'GMAO_PVRT'  ) PVAR(N) = 1011
          IF ( LINE == 'GMAO_PSLV'  ) PVAR(N) = 1012
          IF ( LINE == 'GMAO_AVGW'  ) PVAR(N) = 1013
          IF ( LINE == 'GMAO_THTA'  ) PVAR(N) = 1014
          IF ( LINE == 'GMAO_PRES'  ) PVAR(N) = 1015
          IF ( LINE == 'GMAO_ICE00' ) PVAR(N) = 1100
          IF ( LINE == 'GMAO_ICE10' ) PVAR(N) = 1101
          IF ( LINE == 'GMAO_ICE20' ) PVAR(N) = 1102
          IF ( LINE == 'GMAO_ICE30' ) PVAR(N) = 1103
          IF ( LINE == 'GMAO_ICE40' ) PVAR(N) = 1104
          IF ( LINE == 'GMAO_ICE50' ) PVAR(N) = 1105
          IF ( LINE == 'GMAO_ICE60' ) PVAR(N) = 1106
          IF ( LINE == 'GMAO_ICE70' ) PVAR(N) = 1107
          IF ( LINE == 'GMAO_ICE80' ) PVAR(N) = 1108
          IF ( LINE == 'GMAO_ICE90' ) PVAR(N) = 1109

       !===========================================================
       ! Column aerosol optical depths (same order as for FAST-J)
       ! PVAR offset: 2000
       !===========================================================
       CASE ( 'AODC' )

          IF ( LINE == 'AODC_SULF'  ) PVAR(N) = 2001
          IF ( LINE == 'AODC_BLKC'  ) PVAR(N) = 2002
          IF ( LINE == 'AODC_ORGC'  ) PVAR(N) = 2003
          IF ( LINE == 'AODC_SALA'  ) PVAR(N) = 2004
          IF ( LINE == 'AODC_SALC'  ) PVAR(N) = 2005
          IF ( LINE == 'AODC_DUST'  ) PVAR(N) = 2006

       !===========================================================
       ! Aerosol optical depths below the plane
       ! (same order as for FAST-J)  PVAR offset: 3000
       !===========================================================
       CASE ( 'AODB' )

          IF ( LINE == 'AODB_SULF'  ) PVAR(N) = 3001
          IF ( LINE == 'AODB_BLKC'  ) PVAR(N) = 3002
          IF ( LINE == 'AODB_ORGC'  ) PVAR(N) = 3003
          IF ( LINE == 'AODB_SALA'  ) PVAR(N) = 3004
          IF ( LINE == 'AODB_SALC'  ) PVAR(N) = 3005
          IF ( LINE == 'AODB_DUST'  ) PVAR(N) = 3006

       !===========================================================
       ! Hg(II) Partitioning - eds 10/27/11  PVAR offset: 4000
       !===========================================================
       CASE ( 'HG2_' )

          IF ( LINE == 'HG2_FRACG'  ) PVAR(N) = 4001
          IF ( LINE == 'HG2_FRACP'  ) PVAR(N) = 4002

       !===========================================================
       ! ISORROPIA H+, pH, water, and bisulfate (eam, 06/2015)
       !===========================================================
       CASE( 'ISOR' )

          ! ISORROPIA H+ (mol/L):
          IF ( LINE == 'ISOR_HPLUS' ) PVAR(N) = 5001
          ! ISORROPIA pH (non-ideal system, so pH can be negative)
          IF ( LINE == 'ISOR_PH'    ) PVAR(N) = 5002
          ! ISORROPIA aerosol water (ug/m3 air)
          IF ( LINE == 'ISOR_AH2O'  ) PVAR(N) = 5003
          ! ISORROPIA bisulfate (mol/L):
          IF ( LINE == 'ISOR_HSO4'  ) PVAR(N) = 5004

       !===========================================================
       ! Local time (eam, 06/2015)
       !===========================================================
       CASE( 'TIME' )

          ! Local time:
          IF ( LINE == 'TIME_LT' ) PVAR(N) = 6001

       !===========================================================
       ! Uptake coefficient for SOA formation (eam, 06/2015)
       !===========================================================
       CASE( 'GAMM' )

          ! Skip if not full-chemistry
          IF ( IS_FULLCHEM ) THEN

             ! Increment reaction counter:
             R = R + 1

             IF ( LINE == 'GAMM_EPOX' ) THEN
                PVAR(N)  = 22001
                PREAC(R) = 22001

             ELSE IF ( LINE == 'GAMM_IMAE' ) THEN
                PVAR(N)  = 22002
                PREAC(R) = 22002

             ELSE IF ( LINE == 'GAMM_ISOPN' ) THEN
                PVAR(N)  = 22003
                PREAC(R) = 22003

             ELSE IF ( LINE == 'GAMM_DHDN' ) THEN
                PVAR(N)  = 22004
                PREAC(R) = 22004

             ELSE IF ( LINE == 'GAMM_GLYX' ) THEN
                PVAR(N)  = 22005
                PREAC(R) = 22005

             ENDIF

          ENDIF

       !===========================================================
       ! AQUEOUS AEROSOL properties (eam, 08/2015)
       !===========================================================
       CASE( 'AQAE' )

          ! Aqueous aerosol radius (cm):
          IF ( LINE == 'AQAER_RAD'  ) PVAR(N) = 7001
          ! Aqueous aerosol surface area (cm2/cm3):
          IF ( LINE == 'AQAER_SURF' ) PVAR(N) = 7002

       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       ! This code needs to be updated to work with FlexChem (mps, 6/13/17)
       !!===========================================================
       !! Production rate of CSPEC species (molec/cm3/s)
       !! (eam, 08/2015):
       !!===========================================================
       !CASE( 'PROD' )
       !
       !   ! Skip if not full chem:
       !   IF ( IS_FULLCHEM ) THEN
       !
       !      ! Increment prod rate reaction counter:
       !      PR = PR + 1
       !
       !      ! Initialize:
       !      PRODNAME = ''
       !
       !      ! Extract species name:
       !      PRODNAME = LINE(6:15)
       !
       !      ! Get family integer:
       !      DO NF = 1, NFAMILIES
       !         IF ( FAM_NAME(NF) == TRIM(PRODNAME) ) THEN
       !            ! Loop through family members:
       !            DO FM = 1, FAM_NMEM(NF)
       !               ! Loop through JSPEC species:
       !               DO J = 1, NSPEC(NCS)
       !                  IF ( NAMEGAS(J) == FAM_MEMB(FM, NF) ) THEN
       !                     IPROD(PR) = J
       !                  ENDIF
       !               ENDDO   ! CSPEC species
       !            ENDDO   ! family members
       !         ENDIF
       !      ENDDO
       !
       !      PVAR(N) = 8000 + PR
       !
       !   ENDIF
       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       ! Prior to 6/20/16:
       ! This diagnostic does not work with the FlexChem implementation. We need
       ! to rewrite it to get the necessary information from KPP (mps, 6/20/16)
       !!===========================================================
       !! Rxn rate: listed as "REA_001", etc.
       !! PVAR offset: 10000
       !!===========================================================
       !CASE ( 'REA_' )
       !
       !   ! Skip if not full-chemistry
       !   IF ( IS_FULLCHEM ) THEN
       !
       !      ! Increment rxn counter
       !      R = R + 1
       !
       !      IF ( TRIM( LINE ) == 'REA_O1D' ) THEN
       !
       !         ! O1D is a special rxn, give it offset of 20000
       !         PVAR(N)  = 20000
       !         PREAC(R) = 20000
       !
       !      ELSE IF ( TRIM( LINE ) == 'REA_N2O5' ) THEN
       !
       !         ! N2O5 hydrolysis is another special rxn
       !         ! give it an offset of 21000
       !         PVAR(N)  = 21000
       !         PREAC(R) = 21000
       !
       !      ELSE
       !         !==================================================
       !         ! NOTE: the reaction numbers listed in smv2.log
       !         ! aren't really used to index rxns.  The
       !         ! rxns get reordered.  Find the right rxn number,
       !         ! which is stored in NOLDFNEW.  We assume only one
       !         ! chemistry scheme. (mje, bmy, 8/1/03)
       !         !==================================================
       !
       !         ! Extract tracer # from the string
       !         READ( LINE(5:14), '(i10)' ) NUM
       !
       !         ! Initialize
       !         PVAR(N)  = -999
       !         PREAC(R) = -999
       !
       !         ! Search for proper rxn number
       !         DO IK = 1, NMTRATE
       !
       !            ! Offset other reaction rates by 10000
       !            IF ( NOLDFNEW(IK,1) == NUM ) THEN
       !               PVAR(N)  = 10000 + IK
       !               PREAC(R) = 10000 + IK
       !               EXIT
       !            ENDIF
       !         ENDDO
       !
       !         ! Stop w/ error
       !         IF ( PVAR(N) == -999 ) THEN
       !            IF ( Input_Opt%amIRoot ) THEN
       !               WRITE (6,*) 'Cant match up reaction number'
       !               WRITE (6,*) NUM
       !               WRITE (6,*) 'Is it the second line of the'
       !               WRITE (6,*) 'Three body reaction'
       !               WRITE (6,*) 'Stopping'
       !            ENDIF
       !            CALL GEOS_CHEM_STOP
       !         ENDIF
       !      ENDIF
       !   ENDIF
       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

       !===========================================================
       ! Species: listed as "O3", "C2H6", etc.
       ! PVAR offset: 0
       !===========================================================
       CASE DEFAULT

          ! Skip if not full-chemistry
          IF ( IS_FULLCHEM ) THEN

             ! Loop over all species
             ! match w/ species as read from disk
             DO M = 1, State_Chm%nSpecies

                ! Get info about this species from the species database
                SpcInfo => State_Chm%SpcData(M)%Info

                IF ( TRIM( SpcInfo%Name ) == TRIM( LINE ) ) THEN
                   PVAR(N) = M
                   EXIT
                ENDIF

                ! Free pointer
                SpcInfo => NULL()

             ENDDO

             ! Special flag for RO2 species
             IF ( TRIM( LINE ) == 'RO2' ) PVAR(N) = 999

             ! Special flag for AN species FP
             IF ( TRIM( LINE ) == 'AN' ) PVAR(N) = 998

             ! Special flag for NOy species FP
             IF ( TRIM( LINE ) == 'NOy' ) PVAR(N) = 997

             ! Error check
             IF ( PVAR(N) == 0 ) THEN
                IF ( Input_Opt%amIRoot ) THEN
                   WRITE( 6, '(a)' ) 'ERROR: invalid species!'
                   WRITE( 6, 110   ) TRIM( LINE )
110                FORMAT( 'Species ', a, ' not found!' )
                   WRITE( 6, '(a)' ) 'STOP in PLANEFLIGHT!'
                   CALL GEOS_CHEM_STOP
                ENDIF
             ENDIF
          ENDIF

       END SELECT

       ! Echo species names/numbers to screen
       WRITE( 6, 120 ) N, TRIM( LINE ), PVAR(N)
120    FORMAT( i4, 1x, a12, 1x, i10 )

    ENDDO

  END SUBROUTINE READ_VARIABLES
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_points
!
! !DESCRIPTION: Subroutine READ\_POINTS reads the information (ID, date, time,
!  lat, lon, pressure) for each measurement listed in the input file, and
!  sorts these into the appropriate module variables.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE READ_POINTS( Input_Opt, State_Grid, State_Met, IU_FILE, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE ERROR_MOD,          ONLY : GEOS_CHEM_STOP
    USE FILE_MOD,           ONLY : IOERROR
    USE GC_GRID_MOD,        ONLY : GET_IJ
    USE Input_Opt_Mod,      ONLY : OptInput
    USE Ncdf_Mod,           ONLY : GET_TAU0
    USE PhysConstants,      ONLY : g0
    USE State_Grid_Mod,     ONLY : GrdState
    USE State_Met_Mod,      ONLY : MetState
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)  :: Input_Opt   ! Input Options object
    TYPE(GrdState), INTENT(IN)  :: State_Grid  ! Grid State object
    TYPE(MetState), INTENT(IN)  :: State_Met   ! Meteorology State object
    INTEGER,        INTENT(IN)  :: IU_FILE     ! Logical unit # of file
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT) :: RC          ! Success or failure?
!
! !REVISION HISTORY:
!  30 Jul 2002 - M. Evans    - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER           :: N, IOS, QYY, QMM, QDD, QHH, QMN
    REAL*4            :: LAT, LON, PRES, OBS
    REAL*4            :: TAMB, H2OMR, POTTEMP, GPSALT
    CHARACTER(LEN=7)  :: TYPE
    CHARACTER(LEN=7)  :: NAME

    ! ajt for CCGG
    INTEGER           :: IJ(2), L, L_ALT
    REAL(fp)          :: MOD_ELEV

    !=================================================================
    ! READ_POINTS begins here!
    !=================================================================

    ! Assume success
    RC = GC_SUCCESS

    ! Read 4 header lines
    DO N = 1, 4
       READ( IU_FILE, '(a)', IOSTAT=IOS )
       IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_FILE, 'read_points:1' )
    ENDDO

    !=================================================================
    ! Read plane track points -- plane, lat/lon/alt, date/time
    ! We have previously computed NPOINTS in INIT_PLANEFLIGHT
    !=================================================================
    DO N = 1, NPOINTS

       ! Read a line from the file
       READ( IU_FILE, 100, IOSTAT=IOS ) &
            TYPE, QDD, QMM, QYY, QHH, QMN, LAT, LON, PRES, OBS
100    FORMAT( 5x,a7,x,i2,x,i2,x,i4,x,i2,x,i2,x,f7.2,x,f7.2,x,f7.2,x,f10.3 )

       ! Error check
       IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_FILE, 'read_points:2' )

       ! Exit if the word END is found
       IF ( INDEX( TYPE, 'END' ) > 0 ) EXIT

       !==============================================================
       ! Read date and time coordinates -- also do error checks
       !==============================================================

       ! Error check MONTH
       IF ( QMM < 1 .or. QMM > 12 ) THEN
          IF ( Input_Opt%amIRoot ) THEN
             WRITE( 6, 105 ) QMM
             WRITE( 6, 106 )
105          FORMAT( 'ERROR: MONTH out of range: ', f8.3 )
106          FORMAT( 'STOP in READ_POINTS (planeflight_mod.F90)' )
          ENDIF
          CALL GEOS_CHEM_STOP
       ENDIF

       ! Error check DAY
       IF ( QDD < 1 .or. QDD > 31 ) THEN
          IF ( Input_Opt%amIRoot ) THEN
             WRITE( 6, 110 ) QDD
111          WRITE( 6, 106 )
110          FORMAT( 'ERROR: DAY out of range: ', f8.3 )
          ENDIF
          CALL GEOS_CHEM_STOP
       ENDIF

       ! Error check HOUR
       IF ( QHH < 0 .or. QHH > 23 ) THEN
          IF ( Input_Opt%amIRoot ) THEN
             WRITE( 6, 115 ) QHH
             WRITE( 6, 106 )
115          FORMAT( 'ERROR: HOUR out of range: ', f8.3 )
          ENDIF
          CALL GEOS_CHEM_STOP
       ENDIF

       ! Error check MINUTES
       IF ( QMN < 0 .or. QMN > 59 ) THEN
          IF ( Input_Opt%amIRoot ) THEN
             WRITE( 6, 120 ) QMN
             WRITE( 6, 106 )
120          FORMAT( 'ERROR: MINUTES out of range: ', f8.3 )
          ENDIF
          CALL GEOS_CHEM_STOP
       ENDIF

       ! Store type in the global PTYPE array
       PTYPE(N) = TYPE

       ! Store YYYYMMDD in the global PDATE array
       PDATE(N) = ( QYY * 10000 ) + ( QMM * 100 ) + QDD

       ! Store HHMMSS in the global PTIME array
       ! (actaully we read in just HHMM, assume seconds = 00)
       PTIME(N) = ( QHH * 100 ) + QMN

       ! Store TAU (hours since 1 Jan 1985) in the global PTAU array
       PTAU(N)  = GET_TAU0( QMM, QDD, QYY, QHH, QMN, 0 )

       !==============================================================
       ! Read lon/lat/alt coordinates -- also do error checks
       !==============================================================

       ! Put LONGITUDE in the range [-180...180]
       IF ( LON > 180.0 ) LON = LON - 360e0

       ! Error check LONGITUDE
       IF ( LON < -180.0 .OR. LON > 180.0 ) THEN
          IF ( Input_Opt%amIRoot ) THEN
             WRITE( 6, 125 ) LON
             WRITE( 6, 106 ) 'STOP in READ_POINTS (planeflight_mod.F90)'
125          FORMAT( 'ERROR: Longitude out of range: ', f8.3 )
          ENDIF
          CALL GEOS_CHEM_STOP
       ENDIF

       ! Error check LATITUDE
       IF ( LAT < -90.0 .OR. LAT > 90.0 ) THEN
          IF ( Input_Opt%amIRoot ) THEN
             WRITE( 6, 130 ) LAT
             WRITE( 6, 106 ) 'STOP in READ_POINTS (planeflight_mod.F90)'
130          FORMAT( 'ERROR: Latitude out of range: ', f8.3 )
          ENDIF
          CALL GEOS_CHEM_STOP
       ENDIF

       ! Skip observations outside the domain
       IF ( LAT < State_Grid%YMin .OR. LAT > State_Grid%YMax .OR. &
            LON < State_Grid%XMin .OR. LON > State_Grid%XMax ) THEN
          IF ( Input_Opt%amIRoot ) &
               PRINT*, ' Outside nested domain, skipping record ', N
          CYCLE
       ENDIF

       ! Convert from altitude to pressure if we have CCCG data or
       ! tower data
       NAME = ADJUSTL(PTYPE(N))
       IF ( NAME(1:1)  .EQ. 'S'      .OR. & ! NOAA Surface
            TRIM(NAME) .EQ. 'Aacg'   .OR. & ! NOAA Aircraft
            TRIM(NAME) .EQ. 'Abne'   .OR. & ! NOAA Aircraft
            TRIM(NAME) .EQ. 'Acar'   .OR. & ! NOAA Aircraft
            TRIM(NAME) .EQ. 'Acma'   .OR. & ! NOAA Aircraft
            TRIM(NAME) .EQ. 'Acrv'   .OR. & ! NOAA Aircraft
            TRIM(NAME) .EQ. 'Adnd'   .OR. & ! NOAA Aircraft
            TRIM(NAME) .EQ. 'Aesp'   .OR. & ! NOAA Aircraft
            TRIM(NAME) .EQ. 'Aetl'   .OR. & ! NOAA Aircraft
            TRIM(NAME) .EQ. 'Ahil'   .OR. & ! NOAA Aircraft
            TRIM(NAME) .EQ. 'Ahip'   .OR. & ! NOAA Aircraft
            TRIM(NAME) .EQ. 'Alef'   .OR. & ! NOAA Aircraft
            TRIM(NAME) .EQ. 'Anha'   .OR. & ! NOAA Aircraft
            TRIM(NAME) .EQ. 'Apfa'   .OR. & ! NOAA Aircraft
            TRIM(NAME) .EQ. 'Arta'   .OR. & ! NOAA Aircraft
            TRIM(NAME) .EQ. 'Asca'   .OR. & ! NOAA Aircraft
            TRIM(NAME) .EQ. 'Asgp'   .OR. & ! NOAA Aircraft
            TRIM(NAME) .EQ. 'Atgc'   .OR. & ! NOAA Aircraft
            TRIM(NAME) .EQ. 'Athd'   .OR. & ! NOAA Aircraft
            TRIM(NAME) .EQ. 'Awbi'   .OR. & ! NOAA Aircraft
            TRIM(NAME) .EQ. 'Tamt'   .OR. & ! NOAA Tower
            TRIM(NAME) .EQ. 'Tbao'   .OR. & ! NOAA Tower
            TRIM(NAME) .EQ. 'Tcrv'   .OR. & ! NOAA Tower
            TRIM(NAME) .EQ. 'Tlef'   .OR. & ! NOAA Tower
            TRIM(NAME) .EQ. 'Tlew'   .OR. & ! NOAA Tower
            TRIM(NAME) .EQ. 'Tmbo'   .OR. & ! NOAA Tower
            TRIM(NAME) .EQ. 'Tmvy'   .OR. & ! NOAA Tower
            TRIM(NAME) .EQ. 'Tmwo'   .OR. & ! NOAA Tower
            TRIM(NAME) .EQ. 'Tnwr'   .OR. & ! NOAA Tower
            TRIM(NAME) .EQ. 'Tsct'   .OR. & ! NOAA Tower
            TRIM(NAME) .EQ. 'Tsgp'   .OR. & ! NOAA Tower
            TRIM(NAME) .EQ. 'Tstr'   .OR. & ! NOAA Tower
            TRIM(NAME) .EQ. 'Twbi'   .OR. & ! NOAA Tower
            TRIM(NAME) .EQ. 'Twgc'   .OR. & ! NOAA Tower
            TRIM(NAME) .EQ. 'Twkt' ) THEN   ! NOAA Tower
          ! Change units
          L_ALT = 0
          IJ = GET_IJ( LON, LAT, State_Grid )
          DO L = 1, State_Grid%NZ
             MOD_ELEV = State_Met%PHIS(IJ(1),IJ(2))/ g0 &
                        + SUM( State_Met%BXHEIGHT(IJ(1),IJ(2),1:L) )
             IF ( (L_ALT .EQ. 0)  .AND. (MOD_ELEV .GT. PRES) ) THEN
                L_ALT = L
             ENDIF
          ENDDO
          PRES = State_Met%PMID(IJ(1),IJ(2),L_ALT)
       ENDIF

       ! Assign LAT value into global PLAT array
       PLAT(N)   = LAT

       ! Assign LON value into global PLON array
       PLON(N)   = LON

       ! Assign PRES value into global PPRESS array
       PPRESS(N) = PRES

       ! Assign OBS value into global PPRESS array
       POBS(N) = OBS

    ENDDO

    !=================================================================
    ! Echo number of points found and quit
    !=================================================================
    IF ( Input_Opt%amIRoot ) WRITE( 6, 135 ) NPOINTS
135 FORMAT( /, 'Number of flight track points : ', i6 )

  END SUBROUTINE READ_POINTS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ro2_setup
!
! !DESCRIPTION: Subroutine RO2\_SETUP saves the species indices of RO2
!  constituents in the PRO2 array.  Also computes the count NPRO2.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE RO2_SETUP( Input_Opt, State_Chm, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE ERROR_MOD,          ONLY : GEOS_CHEM_STOP
    USE Input_Opt_Mod,      ONLY : OptInput
    USE Species_Mod,        ONLY : Species
    USE State_Chm_Mod,      ONLY : ChmState
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)  :: Input_Opt   ! Input Options object
    TYPE(ChmState), INTENT(IN)  :: State_Chm   ! Chemistry State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT) :: RC          ! Success or failure?
!
! !REVISION HISTORY:
!  01 Aug 2003 - M. Evans    - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: M

    ! Objects
    TYPE(Species), POINTER :: SpcInfo

    !=================================================================
    ! RO2_SETUP begins here!
    !=================================================================

    ! Initialize
    NPRO2   =  0
    SpcInfo => NULL()
    RC = GC_SUCCESS

    ! We only need to proceed for full-chemistry simulations
    IF ( .not. Input_Opt%ITS_A_FULLCHEM_SIM ) RETURN

    !=================================================================
    ! Loop over all species, test for RO2 components
    !=================================================================
    DO M = 1, State_Chm%nSpecies

       ! Get info about this species from the species database
       SpcInfo => State_Chm%SpcData(M)%Info

       ! If we have found an RO2 compoent, add its species # to
       ! the PRO2 global array, and increment counter
       ! NOTE: PO3 was a bug, that should have been PO2 (tmf, 2/10/09)
       SELECT CASE( TRIM( SpcInfo%Name ) )

       CASE ( 'HO2',     'MO2',      'A3O2',     'ATO2',    'B3O2',  &
              'ETO2',    'HPALD1OO', 'HPALD2OO', 'ICHOO',            &
              'GCO3',    'IAO2',     'KO2',      'MCO3',    'PO2',   &
              'ACO3',    'EO2',      'ENCO3',    'ENO2',    'GLCO3', &
              'ICNOO',   'IDHNBOO',  'IDHNDOO1', 'IDNOO',            &
              'IDHNDOO2','IEPOXAOO', 'IEPOXBOO',                     &
              'IHOO1',   'IHOO4',    'IHPNBOO',  'IHPNDOO',          &
              'IHPOO1',  'IHPOO2',   'IHPOO3',   'INO2B',            &
              'INO2D',   'ISOPNOO1', 'ISOPNOO2', 'LIMO2',            &
              'MACR1OO', 'MACRNO2',  'MCROHOO',  'PIO2',             &
              'MVKOHOO', 'R4O2',     'R4N1',     'R4N2',             &
              'C4HVP1',  'C4HVP2',                                   &
              'BRO2',    'TRO2',     'XRO2',     'NRO2',             &
              'NICO3',   'NIO2',     'PYPO2',    'RCO3')            
          NPRO2       = NPRO2 + 1
          PRO2(NPRO2) = M

       CASE DEFAULT
          ! Nothing

       END SELECT

       ! Free pointer
       SpcInfo => NULL()

    ENDDO

    ! Error check
    IF ( NPRO2 > MAXRO2 ) THEN
       IF ( Input_Opt%amIRoot ) THEN
          WRITE( 6, '(a)' ) 'NPRO2 exceeds maximum allowed value!'
          WRITE( 6, '(a)' ) 'STOP in RO2_SETUP (planeflight_mod.F90)'
          WRITE( 6, '(a)' ) REPEAT( '=', 79 )
          CALL GEOS_CHEM_STOP
       ENDIF
    ENDIF

    !=================================================================
    ! Echo number of points found and quit
    !=================================================================
    IF ( Input_Opt%amIRoot ) THEN
       WRITE( 6, 100 ) NPRO2
100    FORMAT( 'Number of RO2 components      : ', i6 )
    ENDIF

  END SUBROUTINE RO2_SETUP
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: noy_setup
!
! !DESCRIPTION: Subroutine NOY\_SETUP saves the species indices of NOy
!  constituents in the PNOY array.  Also computes the count NPNOY.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE NOY_SETUP( Input_Opt, State_Chm, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE ERROR_MOD,          ONLY : GEOS_CHEM_STOP
    USE Input_Opt_Mod,      ONLY : OptInput
    USE Species_Mod,        ONLY : Species
    USE State_Chm_Mod,      ONLY : ChmState
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)  :: Input_Opt   ! Input Options object
    TYPE(ChmState), INTENT(IN)  :: State_Chm   ! Chemistry State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT) :: RC          ! Success or failure?
!
! !REVISION HISTORY:
!  01 Jun 2009 - F. Paulot   - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: M

    ! Objects
    TYPE(Species), POINTER :: SpcInfo

    !=================================================================
    ! NOY_SETUP begins here!
    !=================================================================

    ! Initialize
    NPNOY   =  0
    SpcInfo => NULL()
    RC = GC_SUCCESS

    ! We only need to proceed for full-chemistry simulations
    IF ( .not. Input_Opt%ITS_A_FULLCHEM_SIM ) RETURN

    !=================================================================
    ! Loop over all species, test for NOY components
    !=================================================================
    DO M = 1, State_Chm%nSpecies

       ! Get info about this species from the species database
       SpcInfo => State_Chm%SpcData(M)%Info

       SELECT CASE( TRIM( SpcInfo%Name ) )

       CASE ( 'NO',  'NO2',   'NO3',  'HNO2', 'HNO4', 'HNO3', &
              'PAN', 'PYPAN', 'MPAN', 'PPN')

          NPNOY       = NPNOY + 1
          PNOY(NPNOY) = M

       CASE ( 'N2O5')

          NPNOY       = NPNOY + 1
          PNOY(NPNOY) = M

          NPNOY       = NPNOY + 1
          PNOY(NPNOY) = M

       CASE DEFAULT
          ! Nothing

       END SELECT

       ! Free pointer
       SpcInfo => NULL()

    ENDDO

    ! Error check
    IF ( NPNOY > MAXNOY ) THEN
       IF ( Input_Opt%amIRoot ) THEN
          WRITE( 6, '(a)' ) 'NPNOY exceeds maximum allowed value!'
          WRITE( 6, '(a)' ) 'STOP in NOY_SETUP (planeflight_mod.F90)'
          WRITE( 6, '(a)' ) REPEAT( '=', 79 )
          CALL GEOS_CHEM_STOP
       ENDIF
    ENDIF

    !=================================================================
    ! Echo number of points found and quit
    !=================================================================
    IF ( Input_Opt%amIRoot ) THEN
       WRITE( 6, 100 ) NPNOY
100    FORMAT( 'Number of NOY components      : ', i6 )
    ENDIF

  END SUBROUTINE NOY_SETUP
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: an_setup
!
! !DESCRIPTION: Subroutine AN\_SETUP saves the species indices of AN
!  constituents in the P\_AN array.  Also computes the count NPAN.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE AN_SETUP( Input_Opt, State_Chm, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE ERROR_MOD,          ONLY : GEOS_CHEM_STOP
    USE Input_Opt_Mod,      ONLY : OptInput
    USE Species_Mod,        ONLY : Species
    USE State_Chm_Mod,      ONLY : ChmState
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)  :: Input_Opt   ! Input Options object
    TYPE(ChmState), INTENT(IN)  :: State_Chm   ! Chemistry State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT) :: RC          ! Success or failure?
!
! !REVISION HISTORY:
!  01 Jun 2009 - F. Paulot   - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: M

    ! Objects
    TYPE(Species), POINTER :: SpcInfo

    !=================================================================
    ! AN_SETUP begins here!
    !=================================================================

    ! Initialize
    NPAN    =  0
    SpcInfo => NULL()
    RC = GC_SUCCESS

    ! We only need to proceed for full-chemistry simulations
    IF ( .not. Input_Opt%ITS_A_FULLCHEM_SIM ) RETURN

    !=================================================================
    ! Loop over all species, test for AN components
    !=================================================================
    DO M = 1, State_Chm%nSpecies

       ! Get info about this species from the species database
       SpcInfo => State_Chm%SpcData(M)%Info

       ! If we have found an AN component, add its species # to
       ! the AN global array, and increment counter
       SELECT CASE( TRIM( SpcInfo%Name ) )

       CASE ( 'IHN1',   'IHN2',   'MVKN',  'INPD',   'R4N2',  &
              'INPB',   'PROPNN', 'ETHLN', 'IDN',    'HONIT', &
              'ITCN',   'ITHN',   'MCRHN', 'MCRHNB',          &
              'MONITU', 'MONITS', 'PRPN',  'IHN3',   'IHN4' )

          NPAN       = NPAN + 1
          P_AN(NPAN) = M

       CASE DEFAULT
          ! Nothing

       END SELECT

       ! Free pointer
       SpcInfo => NULL()

    ENDDO

    ! Error check
    IF ( NPAN > MAXAN ) THEN
       IF ( Input_Opt%amIRoot ) THEN
          WRITE( 6, '(a)' ) 'NPAN exceeds maximum allowed value!'
          WRITE( 6, '(a)' ) 'STOP in AN_SETUP (planeflight_mod.F90)'
          WRITE( 6, '(a)' ) REPEAT( '=', 79 )
          CALL GEOS_CHEM_STOP
       ENDIF
    ENDIF

    !=================================================================
    ! Echo number of points found and quit
    !=================================================================
    IF ( Input_Opt%amIRoot ) THEN
       WRITE( 6, 100 ) NPAN
100    FORMAT( 'Number of AN components      : ', i6 )
    ENDIF

  END SUBROUTINE AN_SETUP
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: planeflight
!
! !DESCRIPTION: Subroutine PLANEFLIGHT saves concentrations to disk at
!  locations corresponding to a flight track.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE PLANEFLIGHT( Input_Opt,  State_Chm, State_Diag, &
                          State_Grid, State_Met, RC )
!
! !USES:
!
    USE CMN_FJX_MOD,        ONLY : ODAER, QAA, QAA_AOD, ODMDUST
    USE CMN_FJX_MOD,        ONLY : IWVSELECT, ACOEF_WV, BCOEF_WV
    USE CMN_SIZE_MOD,       ONLY : NDUST, NAER
    USE ErrCode_Mod
    USE ERROR_MOD,          ONLY : GEOS_CHEM_STOP
    USE Input_Opt_Mod,      ONLY : OptInput
    USE Ncdf_Mod,           ONLY : GET_TAU0
    USE OCEAN_MERCURY_MOD,  ONLY : Fg !eds 10/27/11
    USE OCEAN_MERCURY_MOD,  ONLY : OMMFp => Fp
    USE PhysConstants,      ONLY : CONSVAP, AIRMW
    USE State_Chm_Mod,      ONLY : ChmState
    USE State_Diag_Mod,     ONLY : DgnState
    USE State_Grid_Mod,     ONLY : GrdState
    USE State_Met_Mod,      ONLY : MetState
    USE TIME_MOD
    USE UnitConv_Mod,       ONLY : Convert_Spc_Units
#ifdef TOMAS
#ifdef BPCH_DIAG
    USE DIAG_MOD,           ONLY : AD61_INST   ! (win, 7/28/09)
#endif
#endif
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
    TYPE(GrdState), INTENT(IN)    :: State_Grid  ! Grid State object
    TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology State object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
    TYPE(DgnState), INTENT(INOUT) :: State_Diag  ! Diagnostics State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!
! !REVISION HISTORY:
!  08 Jul 2002 - M. Evans - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    LOGICAL, SAVE       :: FIRST = .TRUE.
    LOGICAL             :: IS_FULLCHEM, LINTERP
    INTEGER             :: I, J, L, M, N, R, PR, V
    INTEGER             :: LL     !eds
    INTEGER             :: IWV, ISPC, K
    REAL(fp)            :: TK, PTAUS, PTAUE, CONSEXP, VPRESH2O, SAODnm
    REAL(fp)            :: VARI(NPVAR)
    LOGICAL             :: CHEMSTEP
    REAL*8              :: FLTGMT   ! eam (06/2015)
    REAL*8              :: XRH      ! MET field RH (eam, 08/2015)
    CHARACTER(LEN=63)   :: OrigUnit
    CHARACTER(LEN=7)    :: NAME

    ! Aerosol types: SULF, BLKC, ORGC, SALA, SALC
    INTEGER             :: IND(5) = (/ 22, 29, 36, 43, 50 /)

    ! Allow for more accurate computation of TAU (L. Schiferl, 1/12/15)
    INTEGER             :: YEAR, MONTH, DAY, HOUR, MINUTE

    ! Pointers
    REAL(fp), POINTER   :: Spc(:,:,:,:)
!
! !DEFINED PARAMETERS:
!
    REAL(fp), PARAMETER :: MISSING = -999.99999999e+0_fp  ! Missing data value
    REAL(fp), PARAMETER :: TINY    = 1.e-36_fp            ! arbitary small # to
                                                          !  avoid faulty output
    ! Expand from 4 to 5 for Fast-JX
    INTEGER, PARAMETER  :: IND999  = 4

    REAL*8, PARAMETER   :: CRITRH  = 35.0e+0_fp

    !=================================================================
    ! PLANEFLIGHT begins here!
    !=================================================================

    ! Assume success
    RC = GC_SUCCESS

    ! Get fields from Input_Opt
    IS_FULLCHEM = Input_Opt%ITS_A_FULLCHEM_SIM

    ! Return if there is no flighttrack data for today
    IF ( .not. DO_PF ) RETURN

    ! Update from kyu (03/2015):
    CHEMSTEP = ( MOD(GET_ELAPSED_SEC(), GET_TS_DIAG() ) == ( GET_TS_DIAG() / 2))

    ! Get date & time values
    YEAR   = GET_YEAR()
    MONTH  = GET_MONTH()
    DAY    = GET_DAY()
    HOUR   = GET_HOUR()
    MINUTE = GET_MINUTE()

    !Determine if optical properties need interpolating
    !The LUT wavelengths in IWVSELECT will match if no interpolation
    !is needed. (output is only for the first requested wavelength)
    !(DAR 10/2013)
    IF(IWVSELECT(1,1).EQ.IWVSELECT(2,1)) THEN
       LINTERP=.FALSE.
    ELSE
       LINTERP=.TRUE.
    ENDIF

    ! Loop over all the locations that have not yet been found
    DO M = PPOINT, NPOINTS

       ! Starting & end times of transport timestep
       PTAUE = GET_TAU0( MONTH, DAY, YEAR, HOUR, MINUTE, 0 )
       ! Modification from kyu (eam, 03/2015)
       !PTAUS = PTAUE - ( GET_TS_DIAG() / 3600e+0_fp )
       ! If we just finished a chemistry timestep, write out one full
       ! diagnostic timestep's worth of data
       IF ( CHEMSTEP ) THEN
          PTAUS = PTAUE - ( GET_TS_DIAG() / 3600d0 )
          ! Otherwise write out only half a timestep's worth of data
       ELSE
          ! Flush the last timestep to the output file (bmy, 3/28/19)
          ! NOTE: Luke Schiferl pointed out we should test against
          ! 86400 seconds instead of GET_TS_DIAG*1440.  GET_TS_DIAG
          ! by now is already in seconds and 1440 is minutes/day,
          ! so there's a mismatch.  This caused planeflight data
          ! for the last timestep of a day not to be written out.
          IF (MOD(GET_ELAPSED_SEC(), 86400) == 0) THEN
             PTAUS = PTAUE - ( GET_TS_DIAG() / (2.0 * 3600d0) )
          ELSE
             EXIT
          ENDIF
       ENDIF

       ! Initialize VARI to missing value for this point
       DO V = 1, NPVAR
          VARI(V) = MISSING
       ENDDO

       !==============================================================
       ! We haven't found the first plane point yet...
       !==============================================================
       IF ( PTAU(M) < PTAUS ) THEN

          ! Write all missing values to disk for point #M
          CALL WRITE_VARS_TO_FILE( Input_Opt, State_Grid, State_Met, M, VARI )

          ! Increment pointer
          PPOINT = PPOINT + 1

       !==============================================================
       ! We have already found all of the plane points...
       !==============================================================
       ! Ensure that a model comparison is made for this point
       ! (skim, 7/24/13)
       ELSE IF ( PTAU(M) > PTAUE ) THEN

          ! Exit this loop and the subroutine
          EXIT

       !==============================================================
       ! We have found a plane point at the proper time & location!
       !==============================================================
       ELSE

          ! Print the flight track point number
          WRITE( 6, 100 ) PTYPE(M), PDATE(M), PTIME(M)
100       FORMAT( '     - PLANEFLIGHT: Archived ',a7,1x,i8.8,1x,i4.4 )

          ! Return grid box indices for the chemistry region
          CALL TEST_VALID( M, I, J, L, Input_Opt, State_Grid, State_Met, RC )

          ! If this is a surface observation, set L=1
          !
          ! NOAA Surface observations start with 'S' in Planeflight.dat
          ! Other surface observation strings can be added here
          NAME = ADJUSTL(PTYPE(M))
          IF ( NAME(1:1)  .EQ. 'S' ) THEN

             ! Make sure it is the users intention to set L=1
             IF ( Input_Opt%amIRoot ) THEN
                WRITE( 6, '(a)') 'WARNING: NOAA Surface Observation. '
                WRITE( 6, '(a)') 'Forcing L=1. If this is intended, '
                WRITE( 6, '(a)') 'you may comment out the call to '
                WRITE( 6, '(a)') 'GEOS_CHEM_STOP in routine '
                WRITE( 6, '(a)') 'PLANEFLIGHT (planeflight_mod.F90)'
                WRITE( 6, '(a)') REPEAT( '=', 79 )
             ENDIF
             CALL GEOS_CHEM_STOP

             ! Force L=1
             L = 1

          ENDIF

          ! Initialize reaction counter
          R = 0

          ! Initialize production count:
          PR = 0

          ! Convert species units to [v/v]
          CALL Convert_Spc_Units( Input_Opt, State_Chm, State_Grid, State_Met, &
                                  'v/v dry', RC, OrigUnit=OrigUnit )

          ! Initialize GEOS-Chem species array
          Spc => State_Chm%Species

          ! Loop over all variables to save out
          DO V = 1, NPVAR

             ! Handle each variable
             SELECT CASE ( PVAR(V) )

             !-------------------------
             ! GEOS-Chem Chemical species [molec/cm3]
             !-------------------------
             CASE ( 1:996)

                ! Only archive where chemistry is done
                IF ( State_Met%InChemGrid(I,J,L) ) THEN

                   ! Species concentration [v/v] -> [molec/cm3]
                   VARI(V) = Spc(I,J,L,PVAR(V)) * State_Met%AIRNUMDEN(I,J,L)

                ENDIF

             ! FP 04/01/2010
             !-------------------------
             ! NOy family
             !-------------------------
             CASE ( 997 )

                ! Only archive where chemistry is done
                ! Sum all AN contributions, save as [v/v]
                VARI(V) = 0e+0_fp

                IF ( IS_FULLCHEM .and. State_Met%InChemGrid(I,J,L) ) THEN

                   DO N = 1, NPNOY

                      ! Species concentration [v/v]
                      VARI(V) = VARI(V) + Spc(I,J,L,PNOY(N))

                   ENDDO

                ENDIF

             ! FP 04/01/2010
             !-------------------------
             ! AN family
             !-------------------------
             CASE ( 998 )

                ! Only archive where chemistry is done
                ! Sum all AN contributions, save as [v/v]
                VARI(V) = 0e+0_fp

                IF ( IS_FULLCHEM .and. State_Met%InChemGrid(I,J,L) ) THEN

                   DO N = 1, NPAN

                      ! Species concentration [v/v]
                      VARI(V) = VARI(V) + Spc(I,J,L,P_AN(N))

                   ENDDO

                ENDIF

             !-------------------------
             ! RO2 family
             !-------------------------
             CASE ( 999 )

                ! Only archive where chemistry is done
                ! Sum all RO2 contributions, save as [v/v]
                VARI(V) = 0e+0_fp

                IF ( IS_FULLCHEM .and. State_Met%InChemGrid(I,J,L) ) THEN

                   DO N = 1, NPRO2

                      ! Species concentration [v/v]
                      VARI(V) = VARI(V) + Spc(I,J,L,PRO2(N))

                   ENDDO

                ENDIF

             !--------------------------
             ! GMAO temperature [K]
             !--------------------------
             CASE ( 1001 )
                VARI(V) = State_Met%T(I,J,L)

             !--------------------------
             ! GMAO abs humidity [frac]
             !--------------------------
             CASE ( 1002 )

                ! Only archive where chemistry is done
                IF ( State_Met%InChemGrid(I,J,L) ) THEN
                   VARI(V)  = State_Met%AVGW(I,J,L) * State_Met%AIRNUMDEN(I,J,L)
                   TK       = State_Met%T(I,J,L)
                   CONSEXP  = 17.2693882e+0_fp * &
                              (TK - 273.16e+0_fp) / (TK - 35.86e+0_fp)

                   VPRESH2O = CONSVAP * EXP(CONSEXP) * 1e+0_fp / TK

                   VARI(V)  = VARI(V) * VPRESH2O / State_Met%AIRNUMDEN(I,J,L)
                ENDIF

             !--------------------------
             ! GMAO aerosol sfc area
             !--------------------------
             CASE ( 1003 )

                ! Only archive where chemistry is done
                VARI(V) = 0e+0_fp

                IF ( State_Met%InChemGrid(I,J,L) ) THEN
                   DO N = 1, NDUST + NAER
                      VARI(V) = VARI(V) + State_Chm%AeroArea(I,J,L,N)
                   ENDDO
                ENDIF

             !--------------------------
             ! GMAO sfc pressure [hPa]
             !--------------------------
             CASE ( 1004 )
                VARI(V) = State_Met%PEDGE(I,J,1)

             !-------------------------
             ! GMAO U-wind [m/s]
             !-------------------------
             CASE ( 1005 )
                VARI(V) = State_Met%U(I,J,L)

             !--------------------------
             ! GMAO V-wind [m/s]
             !--------------------------
             CASE ( 1006 )
                VARI(V) = State_Met%V(I,J,L)

             !--------------------------
             ! GEOS-Chem Grid Box I
             !--------------------------
             CASE ( 1007 )
                VARI(V) = I

             !--------------------------
             ! GEOS-Chem Grid Box J
             !--------------------------
             CASE ( 1008 )
                VARI(V) = J

             !--------------------------
             ! GEOS-Chem Grid Box L
             !--------------------------
             CASE ( 1009 )
                VARI(V) = L

             !--------------------------
             ! GEOS-Chem Relative Humidity [%]
             !--------------------------
             CASE ( 1010 )
                VARI(V) = State_Met%RH(I,J,L)

             !--------------------------
             ! GEOS-Chem Ertel's potential vorticity
             !--------------------------
             CASE ( 1011 )
                ! Disable for now. State_Met%PV is not defined.
                !VARI(V) = State_Met%PV(I,J,L)

             !--------------------------
             ! GEOS-Chem Sea Level pressure [hPa]
             !--------------------------
             CASE ( 1012 )
                VARI(V) = State_Met%SLP(I,J)

             !--------------------------
             ! GEOS-Chem Water Vapor
             !  mixing ratio [v/v]
             !--------------------------
             CASE ( 1013 )
                VARI(V) = State_Met%AVGW(I,J,L)

             !--------------------------
             ! GEOS-Chem Potential Temp
             !  (Theta) [k]
             !  (same calc used in diag1.F90)
             !--------------------------
             CASE ( 1014 )
                VARI(V) = State_Met%T(I,J,L) * &
                     ( State_Met%PEDGE(I,J,1) / State_Met%PMID(I,J,L) )**0.286

             !--------------------------
             ! GEOS-Chem Pressure
             ! at center of grid box [hPa]
             !--------------------------
             CASE ( 1015 )
                VARI(V) = State_Met%PMID(I,J,L)

             !--------------------------
             ! GEOS-Chem SEAICE frac's
             !--------------------------
             CASE ( 1100 )
                VARI(V) = State_Met%SEAICE00(I,J)
             CASE ( 1101 )
                VARI(V) = State_Met%SEAICE10(I,J)
             CASE ( 1102 )
                VARI(V) = State_Met%SEAICE20(I,J)
             CASE ( 1103 )
                VARI(V) = State_Met%SEAICE30(I,J)
             CASE ( 1104 )
                VARI(V) = State_Met%SEAICE40(I,J)
             CASE ( 1105 )
                VARI(V) = State_Met%SEAICE50(I,J)
             CASE ( 1106 )
                VARI(V) = State_Met%SEAICE60(I,J)
             CASE ( 1107 )
                VARI(V) = State_Met%SEAICE70(I,J)
             CASE ( 1108 )
                VARI(V) = State_Met%SEAICE80(I,J)
             CASE ( 1109 )
                VARI(V) = State_Met%SEAICE90(I,J)

             !--------------------------
             ! Column aerosol optical
             ! depths [unitless]
             !--------------------------
             CASE ( 2001:2006 )

                ! Remove MISSING flag
                VARI(V) = 0e+0_fp

                ! Aerosol number
                N = PVAR(V) - 2000

                ! Loop over RH bins
                DO  ISPC= 1, NAER
                   IF ( .not. LINTERP ) THEN
                      DO LL = 1, State_Grid%NZ
                         ! Accumulate
                         VARI(V) = VARI(V) + ODAER(I,J,LL,IWVSELECT(1,1),ISPC)
                      ENDDO
                   ELSE
                      DO LL = 1, State_Grid%NZ
                         ! Interpolated using angstrom exponent between
                         ! Closest available wavelengths
                         ! (coefs pre-calculated in CALC_AOD)
                         !catch any zero values before interpolation
                         IF ((ODAER(I,J,LL,IWVSELECT(2,1),ISPC).GT.0).AND. &
                             (ODAER(I,J,LL,IWVSELECT(1,1),ISPC).GT.0)) THEN
                            VARI(V) = VARI(V) + &
                            (ODAER(I,J,LL,IWVSELECT(2,1),ISPC)*ACOEF_WV(1)**   &
                            (BCOEF_WV(1)*LOG(ODAER(I,J,LL,IWVSELECT(1,1),ISPC)/&
                            ODAER(I,J,LL,IWVSELECT(2,1),ISPC))))
                         ENDIF
                      ENDDO
                   ENDIF
                ENDDO

                !now add in the dust
                IF ( .not. LINTERP ) THEN
                   DO LL = 1, State_Grid%NZ
                      DO ISPC = 1, NDUST
                         ! Accumulate
                         VARI(V) = VARI(V) + ODMDUST(I,J,LL,IWVSELECT(1,1),ISPC)
                      ENDDO
                   ENDDO
                ELSE
                   DO LL = 1, State_Grid%NZ
                      ! Interpolated using angstrom exponent between
                      ! Closest available wavelengths
                      ! (coefs pre-calculated in CALC_AOD)
                      !catch any zero values before interpolation
                      DO ISPC = 1, NDUST
                         IF ((ODAER(I,J,LL,IWVSELECT(2,1),ISPC).GT.0).AND. &
                             (ODAER(I,J,LL,IWVSELECT(1,1),ISPC).GT.0)) THEN
                          VARI(V) = VARI(V) + &
                          (ODMDUST(I,J,LL,IWVSELECT(2,1),ISPC)*ACOEF_WV(1)**   &
                          (BCOEF_WV(1)*LOG(ODMDUST(I,J,LL,IWVSELECT(1,1),ISPC)/&
                          ODMDUST(I,J,LL,IWVSELECT(2,1),ISPC))))
                         ENDIF
                      ENDDO
                   ENDDO
                ENDIF

             !--------------------------
             ! Aerosol optical depths
             ! below plane [unitless]
             !--------------------------
             CASE ( 3001:3006 )

                ! Remove MISSING flag
                VARI(V) = 0e+0_fp

                ! Aerosol number
                N = PVAR(V) - 3000

                ! Loop over RH bins
                DO  ISPC= 1, NAER

                   IF ( .not. LINTERP ) THEN
                      DO LL = 1, State_Grid%NZ
                         ! Skip non-tropospheric boxes
                         IF ( .not. State_Met%InTroposphere(I,J,L) ) CYCLE

                         ! Accumulate
                         VARI(V) = VARI(V) + ODAER(I,J,LL,IWVSELECT(1,1),ISPC)
                      ENDDO
                   ELSE
                      DO LL = 1, State_Grid%NZ
                         ! Skip non-tropospheric boxes
                         IF ( .not. State_Met%InTroposphere(I,J,L) ) CYCLE

                         ! Interpolated using angstrom exponent between
                         ! Closest available wavelengths
                         ! (coefs pre-calculated in CALC_AOD)
                         !catch any zero values before interpolation
                         IF ((ODAER(I,J,LL,IWVSELECT(2,1),ISPC).GT.0).AND. &
                             (ODAER(I,J,LL,IWVSELECT(1,1),ISPC).GT.0)) THEN
                           VARI(V) = VARI(V) + &
                            (ODAER(I,J,LL,IWVSELECT(2,1),ISPC)*ACOEF_WV(1)**   &
                            (BCOEF_WV(1)*LOG(ODAER(I,J,LL,IWVSELECT(1,1),ISPC)/&
                            ODAER(I,J,LL,IWVSELECT(2,1),ISPC))))
                         ENDIF
                      ENDDO
                   ENDIF
                ENDDO

                !now add in the dust
                IF ( .not. LINTERP ) THEN
                   DO LL = 1, State_Grid%NZ
                      ! Skip non-tropospheric boxes
                      IF ( .not. State_Met%InTroposphere(I,J,L) ) CYCLE
                      DO ISPC = 1, NDUST
                         ! Accumulate
                         VARI(V) = VARI(V) + ODMDUST(I,J,LL,IWVSELECT(1,1),ISPC)
                      ENDDO
                   ENDDO
                ELSE
                   DO LL = 1, State_grid%NZ
                      ! Skip non-tropospheric boxes
                      IF ( .not. State_Met%InTroposphere(I,J,L) ) CYCLE

                      ! Interpolated using angstrom exponent between
                      ! Closest available wavelengths
                      ! (coefs pre-calculated in CALC_AOD)
                      !catch any zero values before interpolation
                      DO ISPC = 1, NDUST
                         IF ((ODAER(I,J,LL,IWVSELECT(2,1),ISPC).GT.0).AND. &
                             (ODAER(I,J,LL,IWVSELECT(1,1),ISPC).GT.0)) THEN
                          VARI(V) = VARI(V) + &
                          (ODMDUST(I,J,LL,IWVSELECT(2,1),ISPC)*ACOEF_WV(1)**   &
                          (BCOEF_WV(1)*LOG(ODMDUST(I,J,LL,IWVSELECT(1,1),ISPC)/&
                          ODMDUST(I,J,LL,IWVSELECT(2,1),ISPC))))
                         ENDIF
                      ENDDO
                   ENDDO
                ENDIF

             !--------------------------
             ! Hg(II) partitioning eds 10/27/11
             !--------------------------
             CASE ( 4001 )
                VARI(V) = FG(I,J,LL) !L+1 sample 4/24/12

             CASE ( 4002 )
                VARI(V) = OMMFP(I,J,LL) !L+1 sample 4/24/12

             !--------------------------
             ! ISORROPIA H+ and pH (eam, 06/2015)
             !--------------------------
             CASE( 5001 )
                VARI(V) = State_Chm%HplusSav(I,J,L)

             CASE( 5002 )
                VARI(V) = State_Chm%pHSaV(I,J,L)

             CASE( 5003 )
                VARI(V) = State_Chm%WaterSav(I,J,L)

             CASE( 5004 )
                VARI(V) = State_Chm%BisulSav(I,J,L)

             !--------------------------
             ! Local Time (eam, 06/2015)
             !--------------------------
             CASE( 6001 )

                ! Convert GMT from integer to real and
                ! change format from HHMM to HH.MM:
                FLTGMT = REAL(PTIME(M))*1.d-2

                VARI(V) = GET_LOCALTIME(I,J,L,State_Grid,FLTGMT)

             !--------------------------
             ! Aqueous aerosol properties (eam, 08/2015)
             !--------------------------
             ! MET field relative humidity (%):
                XRH = State_Met%RH( I, J, L )
             CASE( 7001 )
                VARI(V) = 0d0
                ! Radius (cm):
                IF ( XRH .gt. CRITRH ) THEN
                   VARI(V) = State_Chm%AeroRadi(I,J,L,8)
                ENDIF

             CASE( 7002 )
                VARI(V) = 0d0
                ! Surface area (cm2/cm3):
                IF ( XRH .gt. CRITRH ) THEN
                   VARI(V) = State_Chm%AeroArea(I,J,L,8)
                ENDIF

             !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
             ! This code needs to be updated to work with FlexChem (mps,6/13/17)
             !!--------------------------
             !! Production rates (eam, 08/2015)
             !!--------------------------
             !CASE( 8001:8999 )
             !
             !   ! Increment reaction count:
             !   PR = PR + 1
             !
             !   ! Production rate in molec/cm3/s:
             !   ! Only archive where SMVGEAR chem is done
             !   IF ( JLOOP /= 0 ) THEN
             !      VARI(V) = CSPEC( JLOOP, IPROD(PR) )/CHEMINTV
             !
             !      ! Make small values as zero:
             !      IF ( VARI(V) < TINY ) VARI(V) = 0d0
             !
             !      IF ( I .ge. 34 .and. I .le. 37 .and. &
             !           J .ge. 60 .and. J .le. 64 .and. &
             !           L .eq. 1 ) THEN
             !         PRINT*, 'VARI = ', VARI(V)
             !         PRINT*, 'CHEMINTV = ', CHEMINTV
             !         PRINT*, 'CSPEC = ', CSPEC(JLOOP,IPROD(PR))
             !         PRINT*, 'JLOOP,IPROD(PR),PR = ',JLOOP,IPROD(PR),PR
             !      ENDIF
             !   ENDIF
             !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
             ! Prior to 5/12/16:
             ! This diagnostic does not work with the FlexChem implementation.
             ! We need to rewrite ARCHIVE_RXNS_FOR_PF to get the necessary
             ! information from KPP (mps, 5/12/16)
             !!--------------------------
             !! Reaction rates
             !!--------------------------
             !CASE ( 10000:99999 )
             !
             !   ! Increment reaction count
             !   R = R + 1
             !             !   ! Only archive where chemistry is done
             !   VARI(V) = PRRATE(I,J,L,R)
             !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

             !--------------------------
             ! GEOS-CHEM advected species [v/v]
             !--------------------------
             CASE( 100000:199999 )

                ! Remove offset from PVAR
                N = PVAR(V) - 100000

                ! Species concentration [v/v]
                VARI(V) = Spc(I,J,L,N)

                IF ( VARI(V) < TINY ) VARI(V) = 0.e+0_fp

#ifdef TOMAS
#ifdef BPCH_DIAG
             !-------------------------------
             ! TOMAS microphysics rate [kg/s] or [no./cm3/s]
             !-------------------------------
             CASE( 200000:299999 )

                ! Remove offset from PVAR
                N = PVAR(V) - 200000

                ! Archive the microphysics rate
                VARI(V) = AD61_INST(I,J,L,N)

                IF ( Input_Opt%amIRoot ) THEN
                   write (6,*) 'ARCHIVE TO PLANEFLIGHT DIAG', &
                               'AD61_INST at',I,J,L,N,'=',AD61_INST(I,J,L,N)
                ENDIF
#endif
#endif

             !--------------------------
             ! Otherwise it's an error!
             !--------------------------
             CASE DEFAULT
                IF ( Input_Opt%amIRoot) THEN
                   WRITE( 6, '(a)' ) REPEAT( '=', 79 )
                   WRITE( 6, '(a)' ) 'PLANEFLIGHT: Bad variable #!'
                   WRITE( 6, '(a)' ) 'STOP in PLANEFLIGHT!'
                   WRITE( 6, '(a)' ) REPEAT( '=', 79 )
                   CALL GEOS_CHEM_STOP
                ENDIF

             END SELECT

          ENDDO

          ! Convert species units back to original unit
          CALL Convert_Spc_Units( Input_Opt, State_Chm, State_Grid, State_Met, &
                                  OrigUnit, RC )

          ! Free pointer
          NULLIFY( Spc )

          ! Write data for the Mth plane point out to disk
          CALL WRITE_VARS_TO_FILE( Input_Opt, State_Grid, State_Met, M, VARI )

          ! Increment the record pointer
          PPOINT = PPOINT + 1

       ENDIF
    ENDDO

  END SUBROUTINE PLANEFLIGHT
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: test_valid
!
! !DESCRIPTION: Subroutine TEST\_VALID tests to see if we are w/in the
!  tropopause, which is where chemistry is done.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE TEST_VALID( IND, I, J, L, Input_Opt, State_Grid, State_Met, RC )
!
! !USES:
!
    USE Input_Opt_Mod,      ONLY : OptInput
    USE State_Grid_Mod,     ONLY : GrdState
    USE State_Met_Mod,      ONLY : MetState
!
! !INPUT PARAMETERS:
!
    INTEGER,        INTENT(IN ) :: IND         ! # of the flight track point
    TYPE(OptInput), INTENT(IN ) :: Input_Opt   ! Input Options object
    TYPE(GrdState), INTENT(IN ) :: State_Grid  ! Grid State object
    TYPE(MetState), INTENT(IN ) :: State_Met   ! Meteorology State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT) :: I           ! GEOS-Chem longitude index
    INTEGER,        INTENT(OUT) :: J           ! GEOS-Chem latitude index
    INTEGER,        INTENT(OUT) :: L           ! GEOS-Chem level index
    INTEGER,        INTENT(OUT) :: RC          ! Success or failure?

!
! !REVISION HISTORY:
!  08 Jul 2002 - M. Evans - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: IL
    LOGICAL :: FOUND

    !=================================================================
    ! TEST_VALID begins here!
    !=================================================================

    ! We have not found a valid point
    FOUND = .FALSE.

    ! Get I corresponding to PLON(IND)
    I = INT( ( PLON(IND) + 180e+0_fp - &
      ( State_Grid%XMinOffset * State_Grid%DX ) ) / State_Grid%DX + 1.5e+0_fp )

    ! Handle date line correctly (bmy, 4/23/04)
    IF ( I > State_Grid%nx ) I = I - State_Grid%NX

    ! Get J corresponding to PLAT(IND)
    J = INT( ( PLAT(IND) +  90e+0_fp - &
      ( State_Grid%YMinOffset * State_Grid%DY ) ) / State_Grid%DY + 1.5e+0_fp )

    ! Get L corresponding to PRESS(IND)
    L = 1
    DO IL = 1, State_Grid%NZ
       IF ( State_Met%PEDGE(I,J,IL) <= PPRESS(IND) .AND..NOT. FOUND ) THEN
          L     = IL-1
          FOUND =.TRUE.
          EXIT
       ENDIF
    ENDDO

    ! Error check: L must be 1 or higher
    IF ( L == 0 ) L = 1

  END SUBROUTINE TEST_VALID
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: write_vars_to_file
!
! !DESCRIPTION: Subroutine WRITE\_VARS\_TO\_FILE writes the values of all
!  the variables for a given flight track point to the output file.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE WRITE_VARS_TO_FILE( Input_Opt, State_Grid, State_Met, IND, VARI )
!
! !USES:
!
    USE FILE_MOD,       ONLY : IOERROR
    USE GC_GRID_MOD,    ONLY : GET_IJ
    USE Input_Opt_Mod,  ONLY : OptInput
    USE State_Grid_Mod, ONLY : GrdState
    USE State_Met_Mod,  ONLY : MetState
    USE TIME_MOD
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN) :: Input_Opt     ! Input options
    TYPE(GrdState), INTENT(IN) :: State_Grid    ! Grid State object
    TYPE(MetState), INTENT(IN) :: State_Met     ! Meteorology State object
    INTEGER,        INTENT(IN) :: IND           ! # of the flight track point
    REAL(fp),       INTENT(IN) :: VARI(NPVAR)   ! Values to print to file
!
! !REVISION HISTORY:
!  08 Jul 2002 - M. Evans - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    LOGICAL, SAVE :: FIRST = .TRUE.
    INTEGER       :: I, IOS
    INTEGER       :: IIJJ(2), II, JJ, L, LL
    REAL*4        :: LON_TMP, LAT_TMP

    !=================================================================
    ! WRITE_VARS_TO_FILE begins here!
    !=================================================================

    LON_TMP = REAL(PLON(IND),4)
    LAT_TMP = REAL(PLAT(IND),4)

    ! Skip observations outside the domain
    IF ( LAT_TMP < State_Grid%YMin .OR. &
         LAT_TMP > State_Grid%YMax .OR. &
         LON_TMP < State_Grid%XMin .OR. &
         LON_TMP > State_Grid%XMax ) THEN
       IF ( Input_Opt%amIRoot ) &
            PRINT*, ' Outside nested domain, skipping record ', IND
       RETURN
    ENDIF

    ! Get Lat, Lon, and Pressure indicies (ajt, 5/26/13)
    IIJJ = GET_IJ( LON_TMP, LAT_TMP, State_Grid )
    II = IIJJ(1)
    JJ = IIJJ(2)
    LL = 0
    IF ( PPRESS(IND) .GT. State_Met%PEDGE(II,JJ,1) ) LL = 1
    DO L = 1, State_Grid%NZ
       IF ( ( PPRESS(IND) .LT. State_Met%PEDGE(II,JJ,L)   ) .AND. &
            ( PPRESS(IND) .GT. State_Met%PEDGE(II,JJ,L+1) ) ) LL = L
    ENDDO
    IF (LL .EQ. 0) LL = State_Grid%NZ

    ! Write data to file
    WRITE( IU_PLANE, 110, IOSTAT=IOS )                            &
           IND, PTYPE(IND), INT( PDATE(IND) ), INT( PTIME(IND) ), &
           PLAT(IND), PLON(IND), PPRESS(IND), POBS(IND),          &
           INT( GET_ELAPSED_SEC() / GET_TS_DYN() ), LL, II, JJ,   &
           ( VARI(I), I=1,NPVAR )

    ! Format string
#ifdef TOMAS
110 FORMAT(I5,X,A7,X,I8.8,X,I4.4,X,F7.2,X,F7.2,X,F7.2,X,F10.3, &
           X,I9.9,X,I3.3,X,I5.5,X,I5.5,X,250(es11.3e3,x))
#else
110 FORMAT(I5,X,A7,X,I8.8,X,I4.4,X,F7.2,X,F7.2,X,F7.2,X,F10.3, &
           X,I9.9,X,I3.3,X,I5.5,X,I5.5,X,200(es11.3e3,x))
#endif

    ! Error check
    IF ( IOS /= 0 ) CALL IOERROR( IOS,IU_PLANE,'write_vars_to_file:1')

    ! Flush the file to disk
    CALL FLUSH( IU_PLANE )

  END SUBROUTINE WRITE_VARS_TO_FILE
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: set_planeflight
!
! !DESCRIPTION: Subroutine SET\_PLANEFLIGHT is used to pass values read in
!  from the GEOS-Chem input file to "planeflight\_mod.F90".
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE SET_PLANEFLIGHT( PF, IN_FILE, OUT_FILE )
!
! !INPUT PARAMETERS:
!
    LOGICAL,            INTENT(IN) :: PF         ! Turn on planeflight diag?
    CHARACTER(LEN=255), INTENT(IN) :: IN_FILE    ! Input file to read
    CHARACTER(LEN=255), INTENT(IN) :: OUT_FILE   ! Output file to write
!
! !REVISION HISTORY:
!  20 Jul 2004 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC

    ! Save arguments to "shadow" module variables
    DO_PF       = PF
    INFILENAME  = TRIM( IN_FILE  )
    OUTFILENAME = TRIM( OUT_FILE )

  END SUBROUTINE SET_PLANEFLIGHT
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_planeflight
!
! !DESCRIPTION: Subroutine INIT\_PLANEFLIGHT reads the input file to compute
!  the number of variables and flight track points to print out.  Also
!  allocates all module arrays.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE INIT_PLANEFLIGHT( Input_Opt, State_Grid )
!
! !USES:
!
    USE ERROR_MOD,      ONLY : ALLOC_ERR
    USE ERROR_MOD,      ONLY : GEOS_CHEM_STOP
    USE FILE_MOD,       ONLY : IOERROR
    USE Input_Opt_Mod,  ONLY : OptInput
    USE State_Grid_Mod, ONLY : GrdState
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN) :: Input_Opt   ! Input Options object
    TYPE(GrdState), INTENT(IN) :: State_Grid  ! Grid State object
!
! !REVISION HISTORY:
!  08 Jul 2002 - M. Evans - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    LOGICAL           :: IS_INIT = .FALSE.
    INTEGER           :: N, AS, IOS
    CHARACTER(LEN=20) :: LINE

    !=================================================================
    ! INIT_PLANEFLIGHT begins here!
    !=================================================================

    ! Find a free file LUN
    IU_FILE = findFreeLUN()

    ! Open file
    OPEN( IU_FILE, FILE=TRIM( INF ), IOSTAT=IOS )
    IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_FILE, 'init_planeflight:1' )

    ! Read four lines of header
    DO N = 1, 4
       READ( IU_FILE, '(a)', IOSTAT=IOS )
       IF ( IOS /= 0 ) CALL IOERROR( IOS,IU_FILE,'init_planeflight:2')
    ENDDO

    !=================================================================
    ! Read in the number of variables to be output -- store in NPVAR
    !=================================================================
    ! Read in as I4 now for the number of variables we save (skim, 7/24/13)
    READ( IU_FILE, '(i4)', IOSTAT=IOS ) NPVAR
    IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_FILE, 'init_planeflight:3' )

    ! Make sure NPVAR is at least 1
    IF ( NPVAR < 1 ) THEN
       IF ( Input_Opt%amIRoot ) THEN
          WRITE( 6, '(a)') 'NPVAR cannot be zero or negative!'
          WRITE( 6, '(a)') 'STOP in INIT_PLANEFLIGHT (planeflight_mod.F90)'
          WRITE( 6, '(a)') REPEAT( '=', 79 )
       ENDIF
       CALL GEOS_CHEM_STOP
    ENDIF

    ! Make sure NPVAR is less than MAXVARS
    IF ( NPVAR > MAXVARS ) THEN
       IF ( Input_Opt%amIRoot ) THEN
          WRITE( 6, '(a)') 'NPVAR exceeds maximum allowed value!'
          WRITE( 6, '(a)') 'STOP in INIT_PLANEFLIGHT (planeflight_mod.F90)'
          WRITE( 6, '(a)') REPEAT( '=', 79 )
       ENDIF
       CALL GEOS_CHEM_STOP
    ENDIF

    ! Read in a separation line
    READ( IU_FILE, '(a)', IOSTAT=IOS )
    IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_FILE, 'init_planeflight:4' )

    ! Initialize chemistry reaction counter
    NPREAC = 0

    ! Initialize prod rate counter:
    NPROD  = 0

    ! Skip past the species declarations
    DO N = 1, NPVAR
       READ( IU_FILE, '(a)', IOSTAT=IOS ) LINE
       IF ( IOS /= 0 ) CALL IOERROR( IOS,IU_FILE,'init_planeflight:5')

       ! Increment number of chemistry reactions found
       IF ( INDEX( LINE, 'REA_' ) > 0 ) NPREAC = NPREAC + 1
       ! Now also add up the number of GAMMA values
       IF ( INDEX( LINE, 'GAMM' ) > 0 ) NPREAC = NPREAC + 1
       ! Count # of production rate outputs:
       IF ( INDEX( LINE, 'PROD' ) > 0 ) NPROD  = NPROD  + 1
    ENDDO

    ! Read 4 header lines
    DO N = 1, 4
       READ( IU_FILE, '(a)', IOSTAT=IOS )
       IF ( IOS /= 0 ) CALL IOERROR( IOS,IU_FILE,'init_planeflight:6')
    ENDDO

    !=================================================================
    ! Read plane track points -- plane, lat/lon/alt, date/time
    !=================================================================
    NPOINTS = 0

    DO

       ! Read a line from the file
       READ( IU_FILE, '(a)', IOSTAT=IOS ) LINE

       ! Exit at end of file
       IF ( IOS < 0 ) EXIT
       IF ( IOS > 0 ) CALL IOERROR( IOS,IU_FILE,'init_planeflight:7' )

       ! Check for END
       IF ( INDEX( LINE, 'END' ) == 0 ) THEN
          NPOINTS = NPOINTS + 1
       ELSE
          EXIT
       ENDIF
    ENDDO

    ! Close file
    CLOSE( IU_FILE )

    ! If there are no flight-track points then just return
    IF ( NPOINTS < 1 ) THEN
       DO_PF = .FALSE.
       RETURN
    ENDIF

    ! Make sure NPOINTS is less than MAXPOINTS
    IF ( NPOINTS > MAXPOINTS ) THEN
       IF ( Input_Opt%amIRoot ) THEN
          WRITE( 6, '(a)') 'NPOINTS exceeds maximum allowed value!'
          WRITE( 6, '(a)') 'STOP in INIT_PLANEFLIGHT (planeflight_mod.F90)'
          WRITE( 6, '(a)') REPEAT( '=', 79 )
       ENDIF
       CALL GEOS_CHEM_STOP
    ENDIF

    !=================================================================
    ! Allocate arrays to maximum sizes
    !
    ! NOTE: To save space, NPREAC is the actual number of reactions
    !       found.  We will worry about this later.  (bmy, 3/25/05)
    !=================================================================
    IF ( .not. IS_INIT ) THEN

       !-------------------------
       ! Arrays of size NPREAC
       !-------------------------
       ALLOCATE( PREAC( MAX( NPREAC, 1 ) ), STAT=AS )
       IF ( AS /= 0 ) CALL ALLOC_ERR( 'PREAC' )

       ALLOCATE( PRRATE( State_Grid%NX, State_Grid%NY, State_Grid%NZ, &
                         MAX( NPREAC, 1 ) ), STAT=AS )
       IF ( AS /= 0 ) CALL ALLOC_ERR( 'PRRATE' )

       ! ---------------------
       ! Arrays of size NPROD
       ! ---------------------
       ALLOCATE( IPROD( MAX( NPROD, 1 ) ), STAT=AS )
       IF ( AS /= 0 ) CALL ALLOC_ERR( 'IPROD' )

       !--------------------------
       ! Arrays of size MAXVARS
       !--------------------------
       ALLOCATE( PVAR( MAXVARS ), STAT=AS )
       IF ( AS /= 0 ) CALL ALLOC_ERR( 'PVAR' )

       ALLOCATE( PNAME( MAXVARS ), STAT=AS )
       IF ( AS /= 0 ) CALL ALLOC_ERR( 'PNAMES' )

       !---------------------------
       ! Arrays of size MAXPOINTS
       !---------------------------
       ALLOCATE( PTYPE( MAXPOINTS ), STAT=AS )
       IF ( AS /= 0 ) CALL ALLOC_ERR( 'PTYPE' )

       ALLOCATE( PDATE( MAXPOINTS ), STAT=AS )
       IF ( AS /= 0 ) CALL ALLOC_ERR( 'PDATE' )

       ALLOCATE( PTIME( MAXPOINTS ), STAT=AS )
       IF ( AS /= 0 ) CALL ALLOC_ERR( 'PTIME' )

       ALLOCATE( PTAU( MAXPOINTS ), STAT=AS )
       IF ( AS /= 0 ) CALL ALLOC_ERR( 'PTAU' )

       ALLOCATE( PLAT( MAXPOINTS ), STAT=AS )
       IF ( AS /= 0 ) CALL ALLOC_ERR( 'PLAT' )

       ALLOCATE( PLON( MAXPOINTS ), STAT=AS )
       IF ( AS /= 0 ) CALL ALLOC_ERR( 'PLON' )

       ALLOCATE( PPRESS( MAXPOINTS ), STAT=AS )
       IF ( AS /= 0 ) CALL ALLOC_ERR( 'PPRESS' )

       ALLOCATE( POBS( MAXPOINTS ), STAT=AS )
       IF ( AS /= 0 ) CALL ALLOC_ERR( 'POBS' )

       ! Reset IS_INIT flag
       IS_INIT = .TRUE.

    ENDIF

    !=================================================================
    ! Initialize arrays
    !=================================================================
    PREAC  = 0
    NPROD  = 0
    IPROD  = 0
    PRRATE = 0e0
    PVAR   = 0
    PNAME  = ''
    PTYPE  = ''
    PDATE  = 0e0
    PTIME  = 0e0
    PTAU   = 0e0
    PLAT   = 0e0
    PLON   = 0e0
    PPRESS = 0e0
    POBS   = 0e0

  END SUBROUTINE INIT_PLANEFLIGHT
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cleanup_planeflight
!
! !DESCRIPTION: Subroutine CLEANUP\_PLANEFLIGHT deallocates all allocatable
!  module arrays.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE CLEANUP_PLANEFLIGHT
!
! !REVISION HISTORY:
!  01 Jul 2001 - M. Evans - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC

    IF ( ALLOCATED( PVAR   ) ) DEALLOCATE( PVAR   )
    IF ( ALLOCATED( PREAC  ) ) DEALLOCATE( PREAC  )
    IF ( ALLOCATED( IPROD  ) ) DEALLOCATE( IPROD  )
    IF ( ALLOCATED( PNAME  ) ) DEALLOCATE( PNAME  )
    IF ( ALLOCATED( PRRATE ) ) DEALLOCATE( PRRATE )
    IF ( ALLOCATED( PTYPE  ) ) DEALLOCATE( PTYPE  )
    IF ( ALLOCATED( PDATE  ) ) DEALLOCATE( PDATE  )
    IF ( ALLOCATED( PTIME  ) ) DEALLOCATE( PTIME  )
    IF ( ALLOCATED( PTAU   ) ) DEALLOCATE( PTAU   )
    IF ( ALLOCATED( PLAT   ) ) DEALLOCATE( PLAT   )
    IF ( ALLOCATED( PLON   ) ) DEALLOCATE( PLON   )
    IF ( ALLOCATED( PPRESS ) ) DEALLOCATE( PPRESS )
    IF ( ALLOCATED( POBS   ) ) DEALLOCATE( POBS   )

  END SUBROUTINE CLEANUP_PLANEFLIGHT
!EOC
END MODULE PLANEFLIGHT_MOD
