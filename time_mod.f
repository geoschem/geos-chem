! $Id: time_mod.f,v 1.3 2003/10/01 20:32:23 bmy Exp $
      MODULE TIME_MOD
!
!******************************************************************************
!  TIME_MOD contains GEOS-CHEM date and time variables and timesteps, and 
!  routines for accessing them. (bmy, 6/21/00, 9/29/03) 
!
!  Module Variables:
!  ============================================================================
!  (1 ) NYMDb       (INTEGER) : YYYYMMDD at beginning of GEOS-CHEM run
!  (2 ) NHMSb       (INTEGER) : HHMMSS   at beginning of GEOS-CHEM run
!  (3 ) NYMDe       (INTEGER) : YYYYMMDD at end       of GEOS-CHEM run
!  (4 ) NHMSe       (INTEGER) : HHMMSS   at end       of GEOS-CHEM run
!  (5 ) NYMD        (INTEGER) : YYYYMMDD at current timestep
!  (6 ) NHMS        (INTEGER) : HHMMSS   at current timestep
!  (7 ) MONTH       (INTEGER) : Current month value (1-12)
!  (8 ) DAY         (INTEGER) : Current day of the month (1-31)
!  (9 ) YEAR        (INTEGER) : Current year (YYYY format)
!  (10) HOUR        (INTEGER) : Current hour of the day (0-23)
!  (11) MINUTE      (INTEGER) : Current minute of the hour (0-59)
!  (12) SECOND      (INTEGER) : Current second of the minute (0-59)
!  (13) NSEASON     (INTEGER) : Season flag (1=DJF, 2=MAM, 3=JJA, 4=SON)
!  (14) DAY_OF_YEAR (INTEGER) : Current day of year (0-365 or 0-366)
!  (15) ELAPSED_MIN (INTEGER) : Elapsed minutes since the end   of the run
!  (16) TAU         (REAL*8 ) : Current TAU value (hours since 0 GMT 1/1/1985)
!  (17) TAUb        (REAL*8 ) : TAU value at beginning of GEOS-CHEM run
!  (18) TAUe        (REAL*8 ) : TAU value at end       of GEOS-CHEM run
!  (19) DIAGb       (REAL*8 ) : TAU value at beginning of diagnostic interval
!  (20) DIAGe       (REAL*8 ) : TAU value at end       of diagnostic interval
!  (21) GMT         (REAL*8 ) : Current Greenwich Mean Time (0.0 - 23.999 hrs)
!  (22) TS_CHEM     (INTEGER) : Chemistry timestep in minutes
!  (23) TS_CONV     (INTEGER) : Convection timestep in minutes
!  (24) TS_DIAG     (INTEGER) : Diagnostic timestep in minutes
!  (25) TS_DYN      (INTEGER) : Dynamic timestep in minutes
!  (26) TS_EMIS     (INTEGER) : Emission timestep in minutes
!  (27) TS_UNIT     (INTEGER) : Unit conversion timestep in minutes
!  (28) CT_CHEM     (INTEGER) : Number of chemistry timesteps executed so far
!  (29) CT_CONV     (INTEGER) : Number of convection timesteps executed so far
!  (30) CT_DYN      (INTEGER) : Number of dynamic timesteps executed so far
!  (31) CT_EMIS     (INTEGER) : Number of emission timesteps executed so far
!  (32) JD85        (REAL*8 ) : Astronomical Julian Day at 0 GMT 1/1/1985 
!
!  Module Routines:
!  ============================================================================
!  (1 ) SET_CURRENT_TIME  : Updates time variables for current timestep
!  (2 ) SET_BEGIN_TIME    : Initializes NYMDb, NHMSb, TAUb variables
!  (3 ) SET_END_TIME      : Initializes NYMDe, NHMSe, TAUe variables
!  (4 ) SET_DIAGb         : Updates DIAGb and DIAGe diagnostic interval times
!  (5 ) SET_DIAGe         : Updates DIAGb and DIAGe diagnostic interval times
!  (6 ) SET_TIMESTEPS     : Updates the elapsed minutes since the start of run
!  (7 ) SET_CT_CHEM       : Increments/resets the chemistry timestep counter
!  (8 ) SET_CT_CONV       : Increments/resets the convection timestep counter
!  (9 ) SET_CT_DYN        : Increments/resets the dynamic timestep counter
!  (10) SET_CT_EMIS       : Increments/resets the emissions timestep counter
!  (11) SET_CT_A3         : Increments/resets the A-3 fields timestep counter
!  (12) SET_CT_A6         : Increments/resets the A-6 fields timestep counter
!  (13) SET_CT_I6         : Increments/resets the I-6 fields timestep counter
!  (14) SET_ELAPSED_MIN   : Updates the elapsed minutes since the start of run
!  (15) GET_JD            : Returns Astronomical Julian Date for NYMD, NHMS
!  (16) GET_ELAPSED_MIN   : Returns the elapsed minutes since the start of run
!  (17) GET_ELAPSED_SEC   : Returns the elapsed seconds since the start of run 
!  (18) GET_NYMDb         : Returns the YYYYMMDD at the beginning of the run
!  (19) GET_NHMSb         : Returns the HHMMSS   at the beginning of the run
!  (20) GET_NYMDe         : Returns the YYYYMMDD at the end of the run
!  (21) GET_NHMSe         : Returns the HHMMSS   at the end of the run
!  (22) GET_NYMD          : Returns the YYYYMMDD at the current time
!  (23) GET_NHMS          : Returns the HHMMSS   at the current time
!  (24) GET_TIME_AHEAD    : Returns the YYYYMMDD, HHMMSS for N_MINS from now
!  (25) GET_MONTH         : Returns the current month (1-12)
!  (26) GET_DAY           : Returns the current day of month (1-31)
!  (27) GET_YEAR          : Returns the current year (YYYY)
!  (28) GET_HOUR          : Returns the current hour (0-23)
!  (29) GET_MINUTE        : Returns the current minute (0-59)
!  (30) GET_SECOND        : Returns the current second (0-59)
!  (31) GET_DAY_OF_YEAR   : Returns the current day of the year (0-366)
!  (32) GET_GMT           : Returns the current GMT (0.0 - 23.999)
!  (33) GET_TAU           : Returns the current TAU value (hrs since 1/1/1985)
!  (34) GET_TAUb          : Returns TAU value at beginning of GEOS-CHEM run
!  (35) GET_TAUe          : Returns TAU value at end of GEOS-CHEM run
!  (36) GET_DIAGb         : Returns TAU value at start of diagnostic interval
!  (37) GET_DIAGe         : Returns TAU value at end of diagnostic interval
!  (38) GET_LOCALTIME     : Returns local time for a grid box (0.0 - 23.999)
!  (39) GET_SEASON        : Returns season flag (1=DJF, 2=MAM, 3=JJA, 4=SON)
!  (40) GET_TS_CHEM       : Returns chemistry timestep in minutes
!  (41) GET_TS_CONV       : Returns convection timestep in minutes
!  (42) GET_TS_DIAG       : Returns diagnostic timestep in minutes
!  (43) GET_TS_DYN        : Returns dynamic timestep in minutes
!  (44) GET_TS_EMIS       : Returns emissions timestep in minutes
!  (45) GET_TS_UNIT       : Returns unit conversion timestep in minutes
!  (46) GET_CT_CHEM       : Returns # of chemistry timesteps already executed
!  (47) GET_CT_CONV       : Returns # of convection timesteps already executed
!  (48) GET_CT_DYN        : Returns # of dynamic timesteps already executed
!  (49) GET_CT_EMIS       : Returns # of emission timesteps already executed
!  (50) GET_CT_A3         : Returns # of times A-3 fields have been read in
!  (51) GET_CT_A6         : Returns # of times A-6 fields have been read in
!  (52) GET_CT_I6         : Returns # of times I-6 fields have been read in
!  (53) GET_A3_TIME       : Returns YYYYMMDD and HHMMSS for the A-3 fields
!  (54) GET_A6_TIME       : Returns YYYYMMDD and HHMMSS for the A-6 fields
!  (55) GET_I6_TIME       : Returns YYYYMMDD and HHMMSS for the I-6 fields
!  (56) GET_FIRST_A3_TIME : Returns YYYYMMDD and HHMMSS for the first A-3 read
!  (57) GET_FIRST_A3_TIME : Returns YYYYMMDD and HHMMSS for the first A-6 read
!  (58) ITS_TIME_FOR_CHEM : Returns TRUE if it is time to do chemistry
!  (59) ITS_TIME_FOR_CONV : Returns TRUE if it is time to do convection
!  (60) ITS_TIME_FOR_DYN  : Returns TRUE if it is time to do dynamics
!  (61) ITS_TIME_FOR_EMIS : Returns TRUE if it is time to do emissions
!  (62) ITS_TIME_FOR_UNIT : Returns TRUE if it is time to do unit conversions
!  (63) ITS_TIME_FOR_DIAG : Returns TRUE if it is time to write diagnostics
!  (64) ITS_TIME_FOR_A3   : Returns TRUE if it is time to read in A-3 fields
!  (65) ITS_TIME_FOR_A6   : Returns TRUE if it is time to read in A-6 fields
!  (66) ITS_TIME_FOR_I6   : Returns TRUE if it is time to read in I-6 fields
!  (67) ITS_TIME_FOR_UNZIP: Returns TRUE if it is the end of the run
!  (68) ITS_TIME_FOR_DEL  : Returns TRUE if it is time to delete temp files
!  (69) ITS_TIME_FOR_EXIT : Returns TRUE if it is the end of the run
!  (70) ITS_A_LEAPYEAR    : Returns TRUE if the current year is a leapyear
!  (71) NYMD_Y2K          : Returns YYMMDD or YYYYMMDD for the proper data set
!  (72) NYMD6_2_NYMD8     : Converts a 6-digit YYMMDD number into YYYYMMDD
!  (73) NYMD_STRING       : ** deprecated, kept for backwards compatibility **
!  (74) DATE_STRING       : Returns a date string in YYMMDD or YYYYMMDD format
!  (75) TIMESTAMP_STRING  : Returns a string "YYYY/MM/DD HH:MM:SS"
!  (76) PRINT_CURRENT_TIME: Prints date time in YYYY/MM/DD, HH:MM:SS format
!  (77) YMD_EXTRACT       : Extracts YYYY, MM, DD from a YYYYMMDD format number
!  (78) EXPAND_DATE       : Replaces date/time tokens w/ actual values
!
!  GEOS-CHEM modules referenced by time_mod.f
!  ============================================================================
!  (1 ) charpak_mod.f   : Module containing string handling routines
!  (2 ) error_mod.f     : Module containing NaN and other error check routines
!  (3 ) grid_mod.f      : Module containing horizontal grid information
!  (4 ) julday_mod.f    : Module containing astronomical Julian date routines
!
!  NOTES:
!  (1 ) Updated comments (bmy, 9/4/01)
!  (2 ) Added routine YMD_EXTRACT.  Also rewrote TIMECHECK using astronomical
!        Julian day routines from "julday_mod.f". (bmy, 11/21/01)
!  (3 ) Eliminated obsolete code (bmy, 2/27/02)
!  (4 ) Updated comments (bmy, 5/28/02)
!  (5 ) Added routine "expand_date".  Also now reference "charpak_mod.f".
!        (bmy, 6/27/02)
!  (6 ) Now references "error_mod.f".  Also added function GET_SEASON, which
!        returns the current season number. (bmy, 10/22/02)
!  (7 ) Now added module variables and various GET_ and SET_ routines to
!        access them.  Now minutes are the smallest timing unit. (bmy, 3/21/03)
!  (8 ) Bug fix in DATE_STRING (bmy, 5/15/03)
!  (9 ) Added GET_FIRST_A3_TIME and GET_FIRST_A6_TIME.  Also added changes for
!        reading fvDAS fields. (bmy, 6/26/03)
!  (10) Now allow ITS_A_LEAPYEAR to take an optional argument.  Bug fix for
!        Linux: must use ENCODE to convert numbers to strings (bmy, 9/29/03)
!******************************************************************************
!
      IMPLICIT NONE

      !=================================================================
      ! MODULE PRIVATE DECLARATIONS -- keep certain internal variables 
      ! and routines from being seen outside "grid_mod.f"
      !=================================================================
      PRIVATE           :: NYMDb,      NHMSb,       NYMDe   
      PRIVATE           :: NHMSe,      NYMD,        NHMS
      PRIVATE           :: MONTH,      DAY,         YEAR
      PRIVATE           :: HOUR,       MINUTE,      SECOND
      PRIVATE           :: NSEASON,    DAY_OF_YEAR, ELAPSED_MIN
      PRIVATE           :: TAU,        TAUb,        TAUe  
      PRIVATE           :: DIAGb,      DIAGe,       GMT
      PRIVATE           :: TS_CHEM,    TS_CONV,     TS_DIAG
      PRIVATE           :: TS_DYN,     TS_EMIS,     TS_UNIT
      PRIVATE           :: CT_CHEM,    CT_CONV,     CT_DYN
      PRIVATE           :: CT_EMIS,    CT_A3,       CT_A6
      PRIVATE           :: CT_I6,      JD85

      !=================================================================
      ! MODULE VARIABLES
      !=================================================================

      ! Date and time variables
      INTEGER           :: NYMDb,      NHMSb,       NYMDe   
      INTEGER           :: NHMSe,      NYMD,        NHMS
      INTEGER           :: MONTH,      DAY,         YEAR
      INTEGER           :: HOUR,       MINUTE,      SECOND
      INTEGER           :: NSEASON,    DAY_OF_YEAR, ELAPSED_MIN
      REAL*8            :: TAU,        TAUb,        TAUe  
      REAL*8            :: GMT,        DIAGb,       DIAGe

      ! Timesteps
      INTEGER           :: TS_CHEM,   TS_CONV,     TS_DIAG
      INTEGER           :: TS_DYN,    TS_EMIS,     TS_UNIT

      ! Timestep counters
      INTEGER           :: CT_CHEM,   CT_CONV,     CT_DYN    
      INTEGER           :: CT_EMIS,   CT_A3,       CT_A6
      INTEGER           :: CT_I6

      ! Astronomical Julian Date at 0 GMT, 1 Jan 1985
      REAL*8, PARAMETER :: JD85 = 2446066.5d0

      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement 
      !=================================================================
      CONTAINS

!------------------------------------------------------------------------------

      SUBROUTINE SET_CURRENT_TIME
!
!******************************************************************************
!  Subroutine SET_TIME takes in the elapsed time in minutes since the start of
!  a GEOS-CHEM simulation and sets the GEOS-CHEM time variables accordingly.
!  (bmy, 2/5/03)
!******************************************************************************
!
      ! References to F90 modules
      USE JULDAY_MOD, ONLY : JULDAY, CALDATE

      ! Local variables
      REAL*4 :: TMP
      REAL*8 :: JD0, JD1
      
      !=================================================================
      ! SET_CURRENT_TIME begins here!
      !=================================================================

      ! JD0: Astronomical Julian Date at start of GEOS-CHEM run
      JD0 = GET_JD( NYMDb, NHMSb )

      ! JD1: Astronomical Julian Date at current time
      JD1  = JD0 + ( DBLE( ELAPSED_MIN ) / 1440d0 )

      ! Call CALDATE to compute the current YYYYMMDD and HHMMSS
      CALL CALDATE( JD1, NYMD, NHMS )

      ! Extract current year, month, day from NYMD
      CALL YMD_EXTRACT( NYMD, YEAR, MONTH, DAY )

      ! Extract current hour, minute, second from NHMS
      CALL YMD_EXTRACT( NHMS, HOUR, MINUTE, SECOND )

      ! Fix minutes & seconds for display purposes (esp. for 1x1)
      IF ( SECOND              == 59 ) SECOND = 0
      IF ( MOD( MINUTE+1, 10 ) == 0  ) MINUTE = MINUTE + 1

      !=================================================================
      ! Compute other GEOS-CHEM timing variables
      !=================================================================

      ! Current Greenwich Mean Time
      GMT         = ( DBLE( HOUR )            ) + 
     &              ( DBLE( MINUTE ) / 60d0   ) + 
     &              ( DBLE( SECOND ) / 3600d0 )

      ! Days elapsed in this year (0-366)
      DAY_OF_YEAR = JD1 - JULDAY( YEAR, 1, 0d0 )

      ! TAU value (# of hours since 1 Jan 1985)
      ! NOTE: TMP is REAL*4 to prevent precision problems
      TMP         = ( JD1 - JD85 ) * 24e0
      TAU         = DBLE( TMP )

      ! Season index (1=DJF, 2=MAM, 3=JJA, 4=SON)
      SELECT CASE ( MONTH )
         CASE ( 12, 1, 2 )
            NSEASON = 1
         CASE ( 3, 4, 5 )
            NSEASON = 2
         CASE ( 6, 7, 8 )
            NSEASON = 3
         CASE ( 9, 10, 11 )
            NSEASON = 4
      END SELECT

      ! Return to calling program
      END SUBROUTINE SET_CURRENT_TIME

!------------------------------------------------------------------------------

      SUBROUTINE SET_BEGIN_TIME( THISNYMDb, THISNHMSb )
!
!******************************************************************************
!  Subroutine SET_BEGIN_TIME initializes NYMDb, NHMSb, and TAUb, which are the
!  YYYYMMDD, HHMMSS, and hours since 1/1/1985 corresponding to the beginning 
!  date and time of a GEOS-CHEM run. (bmy, 2/5/03)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) THISNYMDb (INTEGER) : YYYYMMDD at beginning of GEOS-CHEM run
!  (2 ) THISNHMSb (INTEGER) : HHMMSS   at beginning of GEOS-CHEM run
!
!  NOTES:
!******************************************************************************
!
      ! References to F90 modules
      USE ERROR_MOD,  ONLY : ERROR_STOP

      ! Arguments
      INTEGER, INTENT(IN) :: THISNYMDb, THISNHMSb
      
      ! Local variables
      REAL*4              :: TMP

      !=================================================================
      ! SET_BEGIN_TIME begins here!
      !=================================================================

      ! Make sure THISNYMDb uses 4 digits for the year
      ! and is not less than 1985/01/01
      IF ( THISNYMDb < 19850101 ) THEN
         CALL ERROR_STOP( 'NYMDb must be in the format YYYYMMDD!',
     &                    'SET_BEGIN_TIME (time_mod.f)' )

      ENDIF

      ! Initialize NYMDb, NHMSb
      NYMDb = THISNYMDb
      NHMSb = THISNHMSb

      ! TAUb value (TMP is REAL*4 to prevent precision problems)
      TMP   = ( GET_JD( NYMDb, NHMSb ) - JD85 ) * 24e0
      TAUb  = DBLE( TMP )

      ! Also initialize ELAPSED_MIN
      ELAPSED_MIN = 0

      ! Return to calling program
      END SUBROUTINE SET_BEGIN_TIME

!------------------------------------------------------------------------------

      SUBROUTINE SET_END_TIME( THISNYMDe, THISNHMSe )
!
!******************************************************************************
!  Subroutine SET_END_TIME initializes NYMDe, NHMSe, and TAUe, which are the
!  YYYYMMDD, HHMMSS, and hours since 1/1/1985 corresponding to the ending 
!  date and time of a GEOS-CHEM run. (bmy, 2/5/03)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) THISNYMDe (INTEGER) : YYYYMMDD at end of GEOS-CHEM run
!  (2 ) THISNHMSe (INTEGER) : HHMMSS   at end of GEOS-CHEM run
!
!  NOTES:
!******************************************************************************
!
      ! Arguments
      INTEGER, INTENT(IN) :: THISNYMDe, THISNHMSe
      
      ! Local variables
      REAL*4              :: TMP

      !=================================================================
      ! SET_END_TIME begins here!
      !=================================================================

      ! Initialize NYMDe, NHMSe
      NYMDe = THISNYMDe
      NHMSe = THISNHMSe

      ! TAUe value (TMP is REAL*4 to prevent precision problems)
      TMP   = ( GET_JD( NYMDe, NHMSe ) - JD85 ) * 24e0
      TAUe  = DBLE( TMP )

      ! Return to calling program
      END SUBROUTINE SET_END_TIME

!------------------------------------------------------------------------------

      SUBROUTINE SET_DIAGb( THISDIAGb )
!
!******************************************************************************
!  Subroutine SET_DIAGb initializes DIAGb, the TAU value at the beginning 
!  of the diagnostic averaging interval. (bmy, 3/21/03)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) THISDIAGb (INTEGER) : TAU value at beginning of diagnostic interval
!
!  NOTES:
!******************************************************************************
!
      ! Arguments
      REAL*8, INTENT(IN) :: THISDIAGb

      !=================================================================
      ! SET_DIAGb begins here!
      !=================================================================
      DIAGb = THISDIAGb

      ! Return to calling program
      END SUBROUTINE SET_DIAGb

!------------------------------------------------------------------------------

      SUBROUTINE SET_DIAGe( THISDIAGe )
!
!******************************************************************************
!  Subroutine SET_DIAGe initializes DIAGe, the TAU value at the end
!  of the diagnostic averaging interval. (bmy, 3/21/03)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) THISDIAGe (INTEGER) : TAU value at end of diagnostic interval
!
!  NOTES:
!******************************************************************************
!
      ! Arguments
      REAL*8, INTENT(IN) :: THISDIAGe

      !=================================================================
      ! SET_DIAGe begins here!
      !=================================================================
      DIAGe = THISDIAGe

      ! Return to calling program
      END SUBROUTINE SET_DIAGe

!------------------------------------------------------------------------------
      
      SUBROUTINE SET_TIMESTEPS( CHEMISTRY, CONVECTION, 
     &                          DYNAMICS,  EMISSION,  UNIT_CONV )
!
!******************************************************************************
!  Subroutine SET_TIMESTEPS initializes the timesteps for dynamics, convection,
!  chemistry, and emissions.  Counters are also zeroed. (bmy, 3/21/03)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) CHEMISTRY  (INTEGER) : Chemistry timestep        [minutes]
!  (2 ) CONVECTION (INTEGER) : Convective timestep       [minutes]
!  (3 ) DYNAMICS   (INTEGER) : Dynamic timestep          [minutes]
!  (4 ) EMISSION   (INTEGER) : Emissions timestep        [minutes]
!  (4 ) UNIT_CONV  (INTEGER) : Unit conversion timestep  [minutes]
!
!  NOTES:
!******************************************************************************
!

      ! Arguments
      INTEGER, INTENT(IN) :: CHEMISTRY, CONVECTION, DYNAMICS
      INTEGER, INTENT(IN) :: EMISSION,  UNIT_CONV
      
      !=================================================================
      ! SET_TIMESTEPS begins here!
      !=================================================================

      ! Initialize timesteps
      TS_CHEM = CHEMISTRY
      TS_CONV = CONVECTION
      TS_DYN  = DYNAMICS
      TS_EMIS = EMISSION
      TS_UNIT = UNIT_CONV

      ! Zero timestep counters
      CT_CHEM = 0
      CT_CONV = 0
      CT_DYN  = 0
      CT_EMIS = 0
      CT_A3   = 0
      CT_A6   = 0
      CT_I6   = 0

      ! Echo to stdout
      WRITE( 6, '(a)' ) REPEAT( '=', 79 )
      WRITE( 6, '(a)' ) 'SET_TIMESTEPS: setting GEOS-CHEM timesteps!'
      WRITE( 6, '(a)' )
      WRITE( 6, '(''Chemistry  Timestep [min] : '', i4 )' ) TS_CHEM
      WRITE( 6, '(''Convection Timestep [min] : '', i4 )' ) TS_CONV
      WRITE( 6, '(''Dynamics   Timestep [min] : '', i4 )' ) TS_DYN
      WRITE( 6, '(''Emission   Timestep [min] : '', i4 )' ) TS_EMIS
      WRITE( 6, '(''Unit Conv  Timestep [min] : '', i4 )' ) TS_UNIT
      WRITE( 6, '(a)' ) REPEAT( '=', 79 )

      ! Return to calling program
      END SUBROUTINE SET_TIMESTEPS

!------------------------------------------------------------------------------

      SUBROUTINE SET_CT_CHEM( INCREMENT, RESET )
!
!******************************************************************************
!  Subroutine SET_CT_CHEM increments CT_CHEM, the counter of chemistry 
!  timesteps executed thus far. (bmy, 3/21/03)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) INCREMENT (LOGICAL) : If specified, then will increment counter
!  (2 ) RESET     (LOGICAL) : If specified, then will reset counter to zero
!
!  NOTES:
!******************************************************************************
!
      ! Arguments
      LOGICAL, INTENT(IN), OPTIONAL :: INCREMENT, RESET

      !=================================================================
      ! SET_CT_CHEM begins here!
      !=================================================================
      IF ( PRESENT( INCREMENT ) ) THEN
         CT_CHEM = CT_CHEM + 1
      ELSE IF ( PRESENT( RESET ) ) THEN
         CT_CHEM = 0
      ENDIF

      ! Return to calling program
      END SUBROUTINE SET_CT_CHEM

!------------------------------------------------------------------------------

      SUBROUTINE SET_CT_CONV( INCREMENT, RESET )
!
!******************************************************************************
!  Subroutine SET_CT_CONV increments CT_CONV, the counter of convection 
!  timesteps executed thus far. (bmy, 3/21/03)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) INCREMENT (LOGICAL) : If T, then will increment counter
!  (2 ) RESET     (LOGICAL) : If T, then will reset counter to zero!
!
!  NOTES:
!******************************************************************************
!
      ! Arguments
      LOGICAL, INTENT(IN), OPTIONAL :: INCREMENT, RESET

      !=================================================================
      ! SET_CT_CONV begins here!
      !=================================================================
      IF ( PRESENT( INCREMENT ) ) THEN
         CT_CONV = CT_CONV + 1
      ELSE IF ( PRESENT( RESET ) ) THEN
         CT_CONV = 0 
      ENDIF

      ! Return to calling program
      END SUBROUTINE SET_CT_CONV

!------------------------------------------------------------------------------

      SUBROUTINE SET_CT_DYN( INCREMENT, RESET )
!
!******************************************************************************
!  Subroutine SET_CT_DYN increments CT_DYN, the counter of dynamic
!  timesteps executed thus far. (bmy, 3/21/03)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) INCREMENT (LOGICAL) : If T, then will increment counter
!  (2 ) RESET     (LOGICAL) : If T, then will reset counter to zero!
!
!  NOTES:
!******************************************************************************
!
      ! Arguments
      LOGICAL, INTENT(IN), OPTIONAL :: INCREMENT, RESET

      !=================================================================
      ! SET_CT_DYN begins here!
      !=================================================================
      IF ( PRESENT( INCREMENT ) ) THEN
         CT_DYN = CT_DYN + 1
      ELSE IF ( PRESENT( RESET ) ) THEN
         CT_DYN = 0
      ENDIF

      ! Return to calling program
      END SUBROUTINE SET_CT_DYN

!------------------------------------------------------------------------------

      SUBROUTINE SET_CT_EMIS( INCREMENT, RESET )
!
!******************************************************************************
!  Subroutine SET_CT_EMIS increments CT_EMIS, the counter of emission
!  timesteps executed thus far. (bmy, 3/21/03)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) INCREMENT (LOGICAL) : If T, then will increment counter
!  (2 ) RESET     (LOGICAL) : If T, then will reset counter to zero!
!
!  NOTES:
!******************************************************************************
!      
      ! Arguments
      LOGICAL, INTENT(IN), OPTIONAL :: INCREMENT, RESET

      !=================================================================
      ! SET_CT_EMIS begins here!
      !=================================================================
      IF ( PRESENT( INCREMENT ) ) THEN
         CT_EMIS = CT_EMIS + 1
      ELSE IF ( PRESENT( RESET ) ) THEN
         CT_EMIS = 0
      ENDIF

      ! Return to calling program
      END SUBROUTINE SET_CT_EMIS

!------------------------------------------------------------------------------

      SUBROUTINE SET_CT_A3( INCREMENT, RESET )
!
!******************************************************************************
!  Subroutine SET_CT_A3 increments CT_A3, the counter of the number of times
!  we have read in A-3 fields.  (bmy, 3/21/03)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) INCREMENT (LOGICAL) : If T, then will increment counter
!  (2 ) RESET     (LOGICAL) : If T, then will reset counter to zero!
!
!  NOTES:
!******************************************************************************
!      
      ! Arguments
      LOGICAL, INTENT(IN), OPTIONAL :: INCREMENT, RESET

      !=================================================================
      ! SET_CT_A3 begins here!
      !=================================================================
      IF ( PRESENT( INCREMENT ) ) THEN
         CT_A3 = CT_A3 + 1
      ELSE IF ( PRESENT( RESET ) ) THEN
         CT_A3 = 0
      ENDIF

      ! Return to calling program
      END SUBROUTINE SET_CT_A3

!------------------------------------------------------------------------------

      SUBROUTINE SET_CT_A6( INCREMENT, RESET )
!
!******************************************************************************
!  Subroutine SET_CT_A6 increments CT_A6, the counter of the number of times
!  we have read in A-6 fields.  (bmy, 3/21/03)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) INCREMENT (LOGICAL) : If T, then will increment counter
!  (2 ) RESET     (LOGICAL) : If T, then will reset counter to zero!
!
!  NOTES:
!******************************************************************************
!      
      ! Arguments
      LOGICAL, INTENT(IN), OPTIONAL :: INCREMENT, RESET

      !=================================================================
      ! SET_CT_A3 begins here!
      !=================================================================
      IF ( PRESENT( INCREMENT ) ) THEN
         CT_A6 = CT_A6 + 1
      ELSE IF ( PRESENT( RESET ) ) THEN
         CT_A6 = 0
      ENDIF

      ! Return to calling program
      END SUBROUTINE SET_CT_A6

!------------------------------------------------------------------------------

      SUBROUTINE SET_CT_I6( INCREMENT, RESET )
!
!******************************************************************************
!  Subroutine SET_CT_I6 increments CT_I6, the counter of the number of times
!  we have read in I-6 fields.  (bmy, 3/21/03)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) INCREMENT (LOGICAL) : If T, then will increment counter
!  (2 ) RESET     (LOGICAL) : If T, then will reset counter to zero!
!
!  NOTES:
!******************************************************************************
!      
      ! Arguments
      LOGICAL, INTENT(IN), OPTIONAL :: INCREMENT, RESET

      !=================================================================
      ! SET_CT_I6 begins here!
      !=================================================================
      IF ( PRESENT( INCREMENT ) ) THEN
         CT_I6 = CT_I6 + 1
      ELSE IF ( PRESENT( RESET ) ) THEN
         CT_I6 = 0
      ENDIF

      ! Return to calling program
      END SUBROUTINE SET_CT_I6

!------------------------------------------------------------------------------

      SUBROUTINE SET_ELAPSED_MIN
!
!******************************************************************************
!  Subroutine SET_ELAPSED_MIN increments the number of elapsed minutes by
!  the dynamic timestep TS_DYN. (bmy, 3/21/03)
!******************************************************************************
!
      !=================================================================
      ! SET_ELAPSED_MIN begins here!
      !=================================================================
      ELAPSED_MIN = ELAPSED_MIN + TS_DYN

      ! Return to calling program
      END SUBROUTINE SET_ELAPSED_MIN

!------------------------------------------------------------------------------

      FUNCTION GET_JD( THISNYMD, THISNHMS ) RESULT( THISJD )
!
!******************************************************************************
!  Function GET_JD is a wrapper for the JULDAY routine.  Given the current
!  NYMD and NHMS values, GET_JD will return the current astronomical Julian
!  date. (bmy, 3/21/03)
! 
!  Arguments as Input:
!  ============================================================================
!  (1 ) THISNYMD (INTEGER) : YYYYMMDD value
!  (2 ) THISNHMS (INTEGER) : HHMMSS value
!
!  NOTES:
!******************************************************************************
!
      ! References to F90 m odules
      USE JULDAY_MOD, ONLY : JULDAY

      ! Arguments
      INTEGER, INTENT(IN)  :: THISNYMD, THISNHMS

      ! Local variables
      INTEGER              :: Y, M, D, H, MI, S
      REAL*8               :: DAY

      ! Function variable
      REAL*8               :: THISJD

      !=================================================================
      ! GET_JD begins here!
      !=================================================================

      ! Extract year, month, day from NYMDb
      CALL YMD_EXTRACT( THISNYMD, Y, M, D )
         
      ! Extract hour, minute, second from NHMSb
      CALL YMD_EXTRACT( THISNHMS, H, MI, S )      

      ! Decimal day (including fractional part) 
      DAY  = DBLE( D ) + ( DBLE( H  ) / 24d0    ) + 
     &                   ( DBLE( MI ) / 1440d0  ) +
     &                   ( DBLE( S  ) / 86400d0 ) 
     
      ! Compute astronomical Julian day at start of run
      THISJD = JULDAY( Y, M, DAY )

      ! Return to the calling program
      END FUNCTION GET_JD

!------------------------------------------------------------------------------

      FUNCTION GET_ELAPSED_MIN() RESULT( THIS_ELAPSED_MIN )
!
!******************************************************************************
!  Function GET_ELAPSED_MIN returns the elapsed minutes since the start of
!  a GEOS_CHEM run to the calling program (bmy, 3/21/03) 
!
!  NOTES:
!******************************************************************************
!
      ! Function value
      INTEGER :: THIS_ELAPSED_MIN

      !=================================================================
      ! GET_ELAPSED_MIN begins here!
      !=================================================================     
      THIS_ELAPSED_MIN = ELAPSED_MIN

      ! Return to calling program
      END FUNCTION GET_ELAPSED_MIN

!------------------------------------------------------------------------------

      FUNCTION GET_ELAPSED_SEC() RESULT( THIS_ELAPSED_SEC )
!
!******************************************************************************
!  Function GET_ELAPSED_SEC returns the elapsed minutss since the start of
!  a GEOS_CHEM run to the calling program (bmy, 3/21/03) 
!
!  NOTES:
!******************************************************************************
!
      ! Function value
      INTEGER :: THIS_ELAPSED_SEC

      !=================================================================
      ! GET_ELAPSED_SEC begins here!
      !=================================================================     
      THIS_ELAPSED_SEC = ELAPSED_MIN * 60

      ! Return to calling program
      END FUNCTION GET_ELAPSED_SEC

!------------------------------------------------------------------------------

      FUNCTION GET_NYMDb() RESULT( THISNYMDb )
!
!******************************************************************************
!  Function GET_NYMDb returns the NYMDb value (YYYYMMDD at the beginning of 
!  the run) to the calling program. (bmy, 3/21/03)
! 
!  NOTES:
!******************************************************************************
!
      ! Function value
      INTEGER :: THISNYMDb

      !=================================================================
      ! GET_NYMDb begins here!
      !=================================================================
      THISNYMDb = NYMDb

      ! Return to calling program
      END FUNCTION GET_NYMDb

!------------------------------------------------------------------------------

      FUNCTION GET_NHMSb() RESULT( THISNHMSb )
!
!******************************************************************************
!  Function GET_NHMSb returns the NHMSb value (HHMMSS at the beginning
!  of the run) to the calling program. (bmy, 3/21/03)
! 
!  NOTES:
!******************************************************************************
!
      ! Function value
      INTEGER :: THISNHMSb

      !=================================================================
      ! GET_NHMSb begins here!
      !=================================================================
      THISNHMSb = NHMSb

      ! Return to calling program
      END FUNCTION GET_NHMSb

!------------------------------------------------------------------------------

      FUNCTION GET_NYMDe() RESULT( THISNYMDe )
!
!******************************************************************************
!  Function GET_NYMDe returns the NYMDe value (YYYYMMDD at the end of 
!  the run) to the calling program. (bmy, 3/21/03)
! 
!  NOTES:
!******************************************************************************
!
      ! Function value
      INTEGER :: THISNYMDe

      !=================================================================
      ! GET_NYMDe begins here!
      !=================================================================
      THISNYMDe = NYMDe

      ! Return to calling program
      END FUNCTION GET_NYMDe

!------------------------------------------------------------------------------

      FUNCTION GET_NHMSe() RESULT( THISNHMSe )
!
!******************************************************************************
!  Function GET_NHMSe returns the NHMSe value (HHMMSS at the end
!  of the run) to the calling program. (bmy, 3/21/03)
! 
!  NOTES:
!******************************************************************************
!
      ! Function value
      INTEGER :: THISNHMSe

      !=================================================================
      ! GET_NHMSe begins here!
      !=================================================================
      THISNHMSe = NHMSe

      ! Return to calling program
      END FUNCTION GET_NHMSe

!------------------------------------------------------------------------------

      FUNCTION GET_NYMD() RESULT( THISNYMD )
!
!******************************************************************************
!  Function GET_NYMD returns the current NYMD value (YYYYMMDD) to the 
!  calling program. (bmy, 2/5/03)
! 
!  NOTES:
!******************************************************************************
!
      ! Function value
      INTEGER :: THISNYMD

      !=================================================================
      ! GET_NYMD begins here!
      !=================================================================
      THISNYMD = NYMD

      ! Return to calling program
      END FUNCTION GET_NYMD

!------------------------------------------------------------------------------

      FUNCTION GET_NHMS() RESULT( THISNHMS )
!
!******************************************************************************
!  Function GET_NHMS returns the current NHMS value (HHMMSS) to the 
!  calling program. (bmy, 2/5/03)
! 
!  NOTES:
!******************************************************************************
!
      ! Function value
      INTEGER :: THISNHMS

      !=================================================================
      ! GET_NHMS begins here!
      !=================================================================
      THISNHMS = NHMS

      ! Return to calling program
      END FUNCTION GET_NHMS

!------------------------------------------------------------------------------

      FUNCTION GET_TIME_AHEAD( N_MINS ) RESULT( DATE )
!
!******************************************************************************
!  Function GET_3h_AHEAD returns to the calling program a 2-element vector
!  containing the YYYYMMDD and HHMMSS values at the current time plus N_MINS
!   minutes. (bmy, 3/21/03)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) N_MINS (INTEGER) : Minutes ahead of time to compute YYYYMMDD,HHMMSS
! 
!  NOTES:
!******************************************************************************
!
      ! References to F90 modules
      USE JULDAY_MOD, ONLY : CALDATE

      ! Arguments
      INTEGER, INTENT(IN) :: N_MINS

      ! Local variables
      INTEGER             :: DATE(2)
      REAL*8              :: JD

      !=================================================================
      ! GET_TIME_AHEAD begins here!
      !=================================================================

      ! Astronomical Julian Date at current time + N_MINS
      JD = GET_JD( NYMD, NHMS ) + ( N_MINS / 1440d0 )
 
      ! Call CALDATE to compute the current YYYYMMDD and HHMMSS
      CALL CALDATE( JD, DATE(1), DATE(2) )

      ! Return to calling program
      END FUNCTION GET_TIME_AHEAD

!------------------------------------------------------------------------------

      FUNCTION GET_MONTH() RESULT( THISMONTH )
!
!******************************************************************************
!  Function GET_MONTH returns the current month to the calling program.
!  (bmy, 2/5/03)
! 
!  NOTES:
!******************************************************************************
!
      ! Function value
      INTEGER :: THISMONTH

      !=================================================================
      ! GET_MONTH begins here!
      !=================================================================
      THISMONTH = MONTH

      ! Return to calling program
      END FUNCTION GET_MONTH

!------------------------------------------------------------------------------

      FUNCTION GET_DAY() RESULT( THISDAY )
!
!******************************************************************************
!  Function GET_DAY returns the current day to the calling program.
!  (bmy, 2/5/03)
!
!  NOTES:
!******************************************************************************
!
      ! Function value
      INTEGER :: THISDAY

      !=================================================================
      ! GET_DAY begins here!
      !=================================================================
      THISDAY = DAY

      ! Return to calling program
      END FUNCTION GET_DAY

!------------------------------------------------------------------------------

      FUNCTION GET_YEAR() RESULT( THISYEAR )
!
!******************************************************************************
!  Function GET_YEAR returns the current year to the calling program.
!  (bmy, 2/5/03)
!******************************************************************************
!
      ! Function value
      INTEGER :: THISYEAR

      !=================================================================
      ! GET_YEAR begins here!
      !=================================================================
      THISYEAR = YEAR

      ! Return to calling program
      END FUNCTION GET_YEAR

!------------------------------------------------------------------------------

      FUNCTION GET_HOUR() RESULT( THISHOUR )
!
!******************************************************************************
!  Function GET_HOUR returns the current hour to the calling program.
!  (bmy, 2/5/03)
!
!  NOTES:
!******************************************************************************
!
      ! Function value
      INTEGER :: THISHOUR

      !=================================================================
      ! GET_HOUR begins here!
      !=================================================================
      THISHOUR = HOUR

      ! Return to calling program
      END FUNCTION GET_HOUR

!------------------------------------------------------------------------------

      FUNCTION GET_MINUTE() RESULT( THISMINUTE )
!
!******************************************************************************
!  Function GET_MINUTE returns the current minute to the calling program
!  (bmy, 2/5/03)
! 
!  NOTES:
!******************************************************************************
!
      ! Function value
      INTEGER :: THISMINUTE

      !=================================================================
      ! GET_MINUTE begins here!
      !=================================================================
      THISMINUTE = MINUTE

      ! Return to calling program
      END FUNCTION GET_MINUTE

!------------------------------------------------------------------------------

      FUNCTION GET_SECOND() RESULT( THISSECOND )
!
!******************************************************************************
!  Function GET_SECOND returns the current seconds to the calling program.
!  (bmy, 2/5/03)
!******************************************************************************
!
      ! Function value
      INTEGER :: THISSECOND

      !=================================================================
      ! GET_SECOND begins here!
      !=================================================================
      THISSECOND = SECOND

      ! Return to calling program
      END FUNCTION GET_SECOND

!------------------------------------------------------------------------------

      FUNCTION GET_DAY_OF_YEAR() RESULT( THISDAYOFYEAR )
!
!******************************************************************************
!  Function GET_DAY_OF_YEAR returns the current day of the year (0-365 or
!  0-366 for leap years) to the calling program. (bmy, 2/5/03)
! 
!  NOTES:
!******************************************************************************
!
      ! Function value
      INTEGER :: THISDAYOFYEAR

      !=================================================================
      ! GET_DAY_OF_YEAR begins here!
      !=================================================================
      THISDAYOFYEAR = DAY_OF_YEAR

      ! Return to calling program
      END FUNCTION GET_DAY_OF_YEAR

!------------------------------------------------------------------------------

      FUNCTION GET_GMT() RESULT( THISGMT )
!
!******************************************************************************
!  Function GET_GMT returns the current Greenwich Mean Time to the calling
!  program. (bmy, 2/5/03)
!
!  NOTES:
!******************************************************************************
!
      ! Function value
      REAL*8 :: THISGMT

      !=================================================================
      ! GET_GMT begins here!
      !=================================================================
      THISGMT = GMT

      ! Return to calling program
      END FUNCTION GET_GMT

!------------------------------------------------------------------------------

      FUNCTION GET_TAU() RESULT( THISTAU )
!
!******************************************************************************
!  Function GET_TAU returns the current TAU (# of hours since 1 Jan 1985) 
!  value to the calling program. (bmy, 2/5/03)
!
!  NOTES:
!******************************************************************************
!
      ! Function value
      REAL*8 :: THISTAU

      !=================================================================
      ! GET_TAUb begins here!
      !=================================================================
      THISTAU = TAU

      ! Return to calling program
      END FUNCTION GET_TAU

!------------------------------------------------------------------------------

      FUNCTION GET_TAUb() RESULT( THISTAUb )
!
!******************************************************************************
!  Function GET_TAUb returns TAUb (# of hours since 1 Jan 1985 at the
!  start of a GEOS-CHEM run) to the calling program. (bmy, 2/5/03)
!
!  NOTES:
!******************************************************************************
!
      ! Function value
      REAL*8 :: THISTAUb

      !=================================================================
      ! GET_TAUb begins here!
      !=================================================================
      THISTAUb = TAUb

      ! Return to calling program
      END FUNCTION GET_TAUb

!------------------------------------------------------------------------------

      FUNCTION GET_TAUe() RESULT( THISTAUe )
!
!******************************************************************************
!  Function GET_TAUe returns TAUe (# of hours since 1 Jan 1985 at the 
!  end of a GEOS-CHEM run) to the calling program. (bmy, 2/5/03)
!
!  NOTES:  
!******************************************************************************
!
      ! Function value
      REAL*8 :: THISTAUe

      !=================================================================
      ! GET_TAUe begins here!
      !=================================================================
      THISTAUe = TAUe

      ! Return to calling program
      END FUNCTION GET_TAUe

!------------------------------------------------------------------------------

      FUNCTION GET_DIAGb() RESULT( THISDIAGb )
!
!******************************************************************************
!  Function GET_DIAGb returns DIAGb (# of hours since 1 Jan 1985 at the
!  start of a diagnostic interval) to the calling program. (bmy, 2/5/03)
!
!  NOTES:
!******************************************************************************
!
      ! Function value
      REAL*8 :: THISDIAGb

      !=================================================================
      ! GET_DIAGb begins here!
      !=================================================================
      THISDIAGb = DIAGb

      ! Return to calling program
      END FUNCTION GET_DIAGb

!------------------------------------------------------------------------------

      FUNCTION GET_DIAGe() RESULT( THISDIAGe )
!
!******************************************************************************
!  Function GET_DIAGb returns DIAGe (# of hours since 1 Jan 1985 at the
!  end of a diagnostic interval) to the calling program. (bmy, 2/5/03)
!
!  NOTES:
!******************************************************************************
!
      ! Function value
      INTEGER :: THISDIAGe

      !=================================================================
      ! GET_DIAGe begins here!
      !=================================================================
      THISDIAGe = DIAGe

      ! Return to calling program
      END FUNCTION GET_DIAGe

!------------------------------------------------------------------------------

      FUNCTION GET_LOCALTIME( I ) RESULT( THISLOCALTIME )
!
!******************************************************************************
!  Function GET_LOCALTIME returns the local time of a grid box to the 
!  calling program. (bmy, 2/5/03)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) I (INTEGER) : Grid box longitude index
!
!  NOTES:
!******************************************************************************
!
      ! References to F90 modules
      USE GRID_MOD, ONLY : GET_XMID

      ! Arguments
      INTEGER, INTENT(IN) :: I

      ! Function value
      REAL*8              :: THISLOCALTIME

      !=================================================================
      ! GET_LOCALTIME begins here!
      !=================================================================

      ! Local Time = GMT + ( longitude / 15 ) since each hour of time
      ! corresponds to 15 degrees of longitude on the globe
      THISLOCALTIME = GET_GMT() + ( GET_XMID( I ) / 15d0 )
      
      ! Make sure that THISLOCALTIME is in the range 0-24 hours
      IF ( THISLOCALTIME > 24 ) THISLOCALTIME = THISLOCALTIME - 24d0 
      IF ( THISLOCALTIME < 0  ) THISLOCALTIME = THISLOCALTIME + 24d0

      ! Return to calling program
      END FUNCTION GET_LOCALTIME

!------------------------------------------------------------------------------
      
      FUNCTION GET_SEASON() RESULT( THISSEASON )
!
!******************************************************************************
!  Function GET_SEASON returns the climatological season number 
!  (1=DJF, 2=MAM, 3=JJA, 4=SON) to the calling program. (bmy, 3/21/03)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) THISMONTH (INTEGER) : Current month (1-12)
!
!  NOTES:
!******************************************************************************
!      
      ! Function value
      INTEGER :: THISSEASON

      !=================================================================
      ! GET_SEASON begins here!
      !=================================================================
      THISSEASON = NSEASON

      ! Return to calling program
      END FUNCTION GET_SEASON

!------------------------------------------------------------------------------

      FUNCTION GET_TS_CHEM() RESULT( THIS_TS_CHEM )
!
!******************************************************************************
!  Function GET_TS_CHEM returns the chemistry timestep in minutes to the
!  calling program. (bmy, 3/21/03)
!
!  NOTES:
!******************************************************************************
!
      ! Function value
      INTEGER :: THIS_TS_CHEM

      !=================================================================
      ! GET_TS_CHEM begins here!
      !=================================================================
      THIS_TS_CHEM = TS_CHEM

      ! Return to calling program
      END FUNCTION GET_TS_CHEM

!------------------------------------------------------------------------------

      FUNCTION GET_TS_CONV() RESULT( THIS_TS_CONV )
!
!******************************************************************************
!  Function GET_TS_CONV returns the convection timestep in minutes to the
!  calling program. (bmy, 3/21/03)
!
!  NOTES:
!******************************************************************************
!
      ! Function value
      INTEGER :: THIS_TS_CONV

      !=================================================================
      ! GET_TS_CONV begins here!
      !=================================================================
      THIS_TS_CONV = TS_CONV

      ! Return to calling program
      END FUNCTION GET_TS_CONV

!------------------------------------------------------------------------------

      FUNCTION GET_TS_DIAG() RESULT( THIS_TS_DIAG )
!
!******************************************************************************
!  Function GET_TS_DIAG returns the diagnostic timestep in minutes to the
!  calling program. (bmy, 3/21/03)
!
!  NOTES:
!******************************************************************************
!
      ! Function value
      INTEGER :: THIS_TS_DIAG

      !=================================================================
      ! GET_TS_DIAG begins here!
      !=================================================================
      THIS_TS_DIAG = TS_DIAG

      ! Return to calling program
      END FUNCTION GET_TS_DIAG

!------------------------------------------------------------------------------

      FUNCTION GET_TS_DYN() RESULT( THIS_TS_DYN )
!
!******************************************************************************
!  Function GET_TS_DIAG returns the diagnostic timestep in minutes to the
!  calling program. (bmy, 3/21/03)
!
!  NOTES:
!******************************************************************************
!
      ! Function value
      INTEGER :: THIS_TS_DYN

      !=================================================================
      ! GET_TS_DYN begins here!
      !=================================================================
      THIS_TS_DYN = TS_DYN

      ! Return to calling program
      END FUNCTION GET_TS_DYN

!------------------------------------------------------------------------------

      FUNCTION GET_TS_EMIS() RESULT( THIS_TS_EMIS )
!
!******************************************************************************
!  Function GET_TS_EMIS returns the emission timestep in minutes to the
!  calling program. (bmy, 3/21/03)
!
!  NOTES:
!******************************************************************************
!
      ! Function value
      INTEGER :: THIS_TS_EMIS

      !=================================================================
      ! GET_TS_EMIS begins here!
      !=================================================================
      THIS_TS_EMIS = TS_EMIS

      ! Return to calling program
      END FUNCTION GET_TS_EMIS

!------------------------------------------------------------------------------

      FUNCTION GET_TS_UNIT() RESULT( THIS_TS_UNIT )
!
!******************************************************************************
!  Function GET_TS_EMIS returns the emission timestep in minutes to the
!  calling program. (bmy, 3/21/03)
!
!  NOTES:
!******************************************************************************
!
      ! Function value
      INTEGER :: THIS_TS_UNIT

      !=================================================================
      ! GET_TS_UNIT begins here!
      !=================================================================
      THIS_TS_UNIT = TS_UNIT

      ! Return to calling program
      END FUNCTION GET_TS_UNIT

!------------------------------------------------------------------------------

      FUNCTION GET_CT_CHEM() RESULT( THIS_CT_CHEM )
!
!******************************************************************************
!  Function GET_CT_CHEM returns the chemistry timestep counter to the
!  calling program. (bmy, 3/21/03)
!
!  NOTES:
!******************************************************************************
!
      ! Function value
      INTEGER :: THIS_CT_CHEM

      !=================================================================
      ! GET_CT_CHEM begins here!
      !=================================================================
      THIS_CT_CHEM = CT_CHEM

      ! Return to calling program
      END FUNCTION GET_CT_CHEM

!------------------------------------------------------------------------------

      FUNCTION GET_CT_CONV() RESULT( THIS_CT_CONV )
!
!******************************************************************************
!  Function GET_CT_CONV returns the convection timestep counter to the
!  calling program. (bmy, 3/21/03)
!
!  NOTES:
!******************************************************************************
!
      ! Function value
      INTEGER :: THIS_CT_CONV

      !=================================================================
      ! GET_CT_CONV begins here!
      !=================================================================
      THIS_CT_CONV = CT_CONV

      ! Return to calling program
      END FUNCTION GET_CT_CONV

!------------------------------------------------------------------------------

      FUNCTION GET_CT_DYN() RESULT( THIS_CT_DYN )
!
!******************************************************************************
!  Function GET_CT_CHEM returns the dynamic timestep counter to the
!  calling program. (bmy, 3/21/03)
!
!  NOTES:
!******************************************************************************
!
      ! Function value
      INTEGER :: THIS_CT_DYN

      !=================================================================
      ! GET_CT_DYN begins here!
      !=================================================================
      THIS_CT_DYN = CT_DYN

      ! Return to calling program
      END FUNCTION GET_CT_DYN

!------------------------------------------------------------------------------

      FUNCTION GET_CT_EMIS() RESULT( THIS_CT_EMIS )
!
!******************************************************************************
!  Function GET_CT_CHEM returns the emission timestep counter to the
!  calling program. (bmy, 3/21/03)
!
!  NOTES:
!******************************************************************************
!
      ! Function value
      INTEGER :: THIS_CT_EMIS

      !=================================================================
      ! GET_CT_EMIS begins here!
      !=================================================================
      THIS_CT_EMIS = CT_EMIS

      ! Return to calling program
      END FUNCTION GET_CT_EMIS

!------------------------------------------------------------------------------

      FUNCTION GET_CT_A3() RESULT( THIS_CT_A3 )
!
!******************************************************************************
!  Function GET_CT_CHEM returns the A-3 fields timestep counter to the
!  calling program. (bmy, 3/21/03)
!
!  NOTES:
!******************************************************************************
!
      ! Function value
      INTEGER :: THIS_CT_A3

      !=================================================================
      ! GET_CT_A3 begins here!
      !=================================================================
      THIS_CT_A3 = CT_A3

      ! Return to calling program
      END FUNCTION GET_CT_A3

!------------------------------------------------------------------------------

      FUNCTION GET_CT_A6() RESULT( THIS_CT_A6 )
!
!******************************************************************************
!  Function GET_CT_A6 returns the A-6 fields timestep counter to the
!  calling program. (bmy, 3/21/03)
!
!  NOTES:
!******************************************************************************
!
      ! Function value
      INTEGER :: THIS_CT_A6

      !=================================================================
      ! GET_CT_A6 begins here!
      !=================================================================
      THIS_CT_A6 = CT_A6

      ! Return to calling program
      END FUNCTION GET_CT_A6

!------------------------------------------------------------------------------

      FUNCTION GET_CT_I6() RESULT( THIS_CT_I6 )
!
!******************************************************************************
!  Function GET_CT_I6 returns the I-6 fields timestep counter to the
!  calling program. (bmy, 3/21/03)
!
!  NOTES:
!******************************************************************************
!
      ! Function value
      INTEGER :: THIS_CT_I6

      !=================================================================
      ! GET_CT_I6 begins here!
      !=================================================================
      THIS_CT_I6 = CT_I6

      ! Return to calling program
      END FUNCTION GET_CT_I6

!------------------------------------------------------------------------------

      FUNCTION GET_A3_TIME() RESULT( DATE )
!
!******************************************************************************
!  Function GET_A3_TIME returns the correct YYYYMMDD and HHMMSS values
!  that are needed to read in the next average 3-hour (A-3) fields. 
!  (bmy, 3/21/03, 6/19/03)
!
!  NOTES:
!  (1 ) Now return proper time for GEOS-4/fvDAS fields (bmy, 6/19/03)
!******************************************************************************
!
#     include "define.h"

      ! Local variables
      LOGICAL, SAVE :: FIRST

      ! Function value
      INTEGER :: DATE(2)

      !=================================================================
      ! GET_A3_TIME begins here!
      !=================================================================
#if   defined( GEOS_4 ) 

      ! For GEOS-4/fvDAS, the A-3 fields are timestamped by center time.
      ! Therefore, the difference between the actual time when the fields 
      ! are read and the A-3 timestamp time is 90 minutes.
      DATE = GET_TIME_AHEAD( 90 )   

#else
      
      ! For GEOS-1, GEOS-STRAT, GEOS-3, the A-3 fields are timestamped 
      ! by ending time.  Therefore, the difference between the actual time
      ! when the fields are read and the A-3 timestamp time is 180 minutes.
      DATE = GET_TIME_AHEAD( 180 )      

#endif

      ! Return to calling program
      END FUNCTION GET_A3_TIME
      
!------------------------------------------------------------------------------

      FUNCTION GET_A6_TIME() RESULT( DATE )
!
!******************************************************************************
!  Function GET_A6_TIME returns the correct YYYYMMDD and HHMMSS values
!  that are needed to read in the next average 6-hour (A-6) fields. 
!  (bmy, 3/21/03, 6/26/03)
!
!  NOTES:
!  (1 ) Updated comments (bmy, 6/26/03)
!******************************************************************************
!
      ! Function value
      INTEGER :: DATE(2)

      !=================================================================
      ! GET_A6_TIME begins here!
      !=================================================================

      ! Return the time 3h (180m) from now, since there is a 3-hour 
      ! offset between the actual time when the A-6 fields are read
      ! and the time that the A-6 fields are stamped with.
      DATE = GET_TIME_AHEAD( 180 )      

      ! Return to calling program
      END FUNCTION GET_A6_TIME

!------------------------------------------------------------------------------

      FUNCTION GET_I6_TIME() RESULT( DATE )
!
!******************************************************************************
!  Function GET_I6_TIME returns the correct YYYYMMDD and HHMMSS values
!  that are needed to read in the next instantaneous 6-hour (I-6) fields. 
!  (bmy, 3/21/03)
!
!  NOTES:
!******************************************************************************
!
      ! Arguments
      INTEGER :: DATE(2)

      !=================================================================
      ! GET_I6_TIME begins here!
      !=================================================================

      ! We need to read in the I-6 fields 6h (360 mins) ahead of time
      DATE = GET_TIME_AHEAD( 360 )
      
      ! Return to calling program
      END FUNCTION GET_I6_TIME

!------------------------------------------------------------------------------

      FUNCTION GET_FIRST_A3_TIME() RESULT( DATE )
!
!******************************************************************************
!  Function GET_FIRST_A3_TIME returns the correct YYYYMMDD and HHMMSS values
!  the first time that A-3 fields are read in from disk. (bmy, 6/26/03)
!
!  NOTES:
!******************************************************************************
!
#     include "define.h"

      ! Arguments
      INTEGER :: DATE(2)

      !=================================================================
      ! GET_FIRST_A3_TIME begins here!
      !=================================================================

#if   defined( GEOS_4 )

      ! For GEOS-4, call GET_A3_TIME to return date/time
      ! under which the A-3 fields are timestamped
      DATE = GET_A3_TIME()
      
#else

      ! For GEOS-1, GEOS-STRAT, GEOS-3, return the current date/time
      DATE = (/ NYMD, NHMS /)

#endif

      ! Return to calling program
      END FUNCTION GET_FIRST_A3_TIME

!------------------------------------------------------------------------------

      FUNCTION GET_FIRST_A6_TIME() RESULT( DATE )
!
!******************************************************************************
!  Function GET_FIRST_A6_TIME returns the correct YYYYMMDD and HHMMSS values
!  the first time that A-6 fields are read in from disk. (bmy, 6/26/03)
!
!  NOTES:
!******************************************************************************
!
#     include "define.h"

      ! Arguments
      INTEGER :: DATE(2)

      !=================================================================
      ! GET_FIRST_A6_TIME begins here!
      !=================================================================

#if   defined( GEOS_4 )

      ! For GEOS-4, call GET_A3_TIME to return date/time
      ! under which the A-3 fields are timestamped
      DATE = GET_A6_TIME()
      
#else

      ! For GEOS-1, GEOS-STRAT, GEOS-3, return the current date/time
      DATE = (/ NYMD, NHMS /)

#endif

      ! Return to calling program
      END FUNCTION GET_FIRST_A6_TIME

!------------------------------------------------------------------------------

      FUNCTION ITS_TIME_FOR_CHEM() RESULT( FLAG )
!
!******************************************************************************
!  Function ITS_TIME_FOR_CHEM returns TRUE if it is time to do chemistry
!  and false otherwise. (bmy, 3/21/03)
!
!  NOTES:
!******************************************************************************
!
      ! Function value
      LOGICAL :: FLAG 

      !=================================================================
      ! ITS_TIME_FOR_CHEM begins here!
      !=================================================================
      FLAG = ( MOD( ELAPSED_MIN, TS_CHEM ) == 0 )

      ! Return to calling program
      END FUNCTION ITS_TIME_FOR_CHEM

!------------------------------------------------------------------------------

      FUNCTION ITS_TIME_FOR_CONV() RESULT( FLAG )
!
!******************************************************************************
!  Function ITS_TIME_FOR_CONV returns TRUE if it is time to do chemistry
!  and false otherwise. (bmy, 3/21/03)
!
!  NOTES:
!******************************************************************************
!
      ! Function value
      LOGICAL :: FLAG 

      !=================================================================
      ! ITS_TIME_FOR_CONV begins here!
      !=================================================================
      FLAG = ( MOD( ELAPSED_MIN, TS_CONV ) == 0 )

      ! Return to calling program
      END FUNCTION ITS_TIME_FOR_CONV

!------------------------------------------------------------------------------

      FUNCTION ITS_TIME_FOR_DYN() RESULT( FLAG )
!
!******************************************************************************
!  Function ITS_TIME_FOR_DYN returns TRUE if it is time to do chemistry
!  and false otherwise. (bmy, 3/21/03)
!
!  NOTES:
!******************************************************************************
!
      ! Function value
      LOGICAL :: FLAG 

      !=================================================================
      ! ITS_TIME_FOR_DYN begins here!
      !=================================================================
      FLAG = ( MOD( ELAPSED_MIN, TS_DYN ) == 0 )

      ! Return to calling program
      END FUNCTION ITS_TIME_FOR_DYN

!------------------------------------------------------------------------------

      FUNCTION ITS_TIME_FOR_EMIS() RESULT( FLAG )
!
!******************************************************************************
!  Function ITS_TIME_FOR_EMIS returns TRUE if it is time to do emissions
!  and false otherwise. (bmy, 3/21/03)
!
!  NOTES:
!******************************************************************************
!
      ! Function value
      LOGICAL :: FLAG 

      !=================================================================
      ! ITS_TIME_FOR_EMIS begins here!
      !=================================================================
      FLAG = ( MOD( ELAPSED_MIN, TS_EMIS ) == 0 )

      ! Return to calling program
      END FUNCTION ITS_TIME_FOR_EMIS

!------------------------------------------------------------------------------

      FUNCTION ITS_TIME_FOR_UNIT() RESULT( FLAG )
!
!******************************************************************************
!  Function ITS_TIME_FOR_UNIT returns TRUE if it is time to do unit conversion
!  and false otherwise. (bmy, 3/21/03)
!
!  NOTES:
!******************************************************************************
!
      ! Function value
      LOGICAL :: FLAG 

      !=================================================================
      ! ITS_TIME_FOR_UNIT begins here!
      !=================================================================
      FLAG = ( MOD( ELAPSED_MIN, TS_DYN ) == 0 )

      ! Return to calling program
      END FUNCTION ITS_TIME_FOR_UNIT

!------------------------------------------------------------------------------

      FUNCTION ITS_TIME_FOR_DIAG() RESULT( FLAG )
!
!******************************************************************************
!  Function ITS_TIME_FOR_DIAG returns TRUE if it is time to archive
!  certain diagnostics false otherwise. (bmy, 3/21/03)
!
!  NOTES:
!******************************************************************************
!
      ! Function value
      LOGICAL :: FLAG 

      !=================================================================
      ! ITS_TIME_FOR_DIAG begins here!
      !=================================================================
      FLAG = ( MOD( ELAPSED_MIN, 60 ) == 0 )

      ! Return to calling program
      END FUNCTION ITS_TIME_FOR_DIAG

!------------------------------------------------------------------------------

      FUNCTION ITS_TIME_FOR_A3() RESULT( FLAG )
!
!******************************************************************************
!  Function ITS_TIME_FOR_A3 returns TRUE if it is time to read in A-3
!  (average 3-h fields) and FALSE otherwise. (bmy, 3/21/03)
!
!  NOTES:
!******************************************************************************
!
      ! Function value
      LOGICAL :: FLAG 

      !=================================================================
      ! ITS_TIME_FOR_A3 begins here!
      !=================================================================

      ! We read A-3 fields every 3 hours
      FLAG = ( MOD( NHMS, 030000 ) == 0 )

      ! Return to calling program
      END FUNCTION ITS_TIME_FOR_A3

!------------------------------------------------------------------------------

      FUNCTION ITS_TIME_FOR_A6() RESULT( FLAG )
!
!******************************************************************************
!  Function ITS_TIME_FOR_A6 returns TRUE if it is time to read in A-6
!  (average 6-h fields) and FALSE otherwise. (bmy, 3/21/03, 6/26/03)
!
!  NOTES:
!  (1 ) Now compute when it's time to read in GEOS-4 A-6 fields. (bmy, 6/26/03)
!******************************************************************************
!
      ! Local variables
      INTEGER :: DATE(2)

      ! Function value
      LOGICAL :: FLAG 

      !=================================================================
      ! ITS_TIME_FOR_A6 begins here!
      !=================================================================
      
#if   defined( GEOS_4 )

      ! For GEOS-4, we need to read A-6 fields at 00, 06, 12, 18 GMT
      DATE = GET_TIME_AHEAD( 0 )

#else

      ! For GEOS-1, GEOS-S, GEOS-3, we need to read A-6 fields at
      ! 03, 09, 15, 21 GMT.  There is a 3 hour difference between
      ! the actual time and the A-3 met field timestamp time.  
      DATE = GET_TIME_AHEAD( 180 )
     
#endif
 
      ! If the time now (GEOS-4) or 3 hours from now (GEOS-1, GEOS-S, 
      ! GEOS-3) is 0, 6, 12, 18 GMT, then read A-6 fields 
      FLAG = ( MOD( DATE(2), 060000 ) == 0 ) 

      ! Return to calling program
      END FUNCTION ITS_TIME_FOR_A6

!------------------------------------------------------------------------------

      FUNCTION ITS_TIME_FOR_I6() RESULT( FLAG )
!
!******************************************************************************
!  Function ITS_TIME_FOR_I6 returns TRUE if it is time to read in I-6
!  (instantaneous 6-h fields) and FALSE otherwise. (bmy, 3/21/03)
!
!  NOTES:
!******************************************************************************
!
      ! Function value
      LOGICAL :: FLAG 

      !=================================================================
      ! ITS_TIME_FOR_I6 begins here!
      !=================================================================

      ! We read in I-6 fields every 6 hours
      FLAG = ( MOD( NHMS, 060000 ) == 0 )

      ! Return to calling program
      END FUNCTION ITS_TIME_FOR_I6

!------------------------------------------------------------------------------

      FUNCTION ITS_TIME_FOR_UNZIP() RESULT( FLAG )
!
!******************************************************************************
!  Function ITS_TIME_FOR_UNZIP Treturns TRUE if it is time to unzip
!  the next day's met field files (bmy, 3/21/03)
!
!  NOTES:
!******************************************************************************
!
      ! Local variables
      INTEGER :: DATE(2)

      ! Function value
      LOGICAL :: FLAG 

      !=================================================================
      ! ITS_TIME_FOR_UNZIP begins here!
      !=================================================================

      ! Get YYYYMMDD and HHMMSS 12 hours (720 mins) from now
      DATE = GET_TIME_AHEAD( 720 )

      ! If HHMMSS = 0 then it's time to unzip!
      FLAG = ( DATE(2) == 000000 )

      ! Return to calling program
      END FUNCTION ITS_TIME_FOR_UNZIP     

!------------------------------------------------------------------------------

      FUNCTION ITS_TIME_FOR_DEL() RESULT( FLAG )
!
!******************************************************************************
!  Function ITS_TIME_FOR_DEL returns TRUE if it is time to delete the previous
!  day's met field files in the temporary directory. (bmy, 3/21/03, 6/19/03)
!
!  NOTES:
!  (1 ) Now delete files at 23 GMT each day, since the last fvDAS A-3 field is 
!        22:30 GMT and the last fvDAS A-6 field is 21 GMT. (bmy, 6/19/03)
!******************************************************************************
!
      ! Local variables
      INTEGER :: DATE(2)

      ! Function value
      LOGICAL :: FLAG 

      !=================================================================
      ! ITS_TIME_FOR_DEL begins here!
      !=================================================================

      ! Delete files when it's 23 GMT
      FLAG = ( NHMS == 230000 )

      ! Return to calling program
      END FUNCTION ITS_TIME_FOR_DEL

!------------------------------------------------------------------------------

      FUNCTION ITS_TIME_FOR_EXIT() RESULT( FLAG )
!
!******************************************************************************
!  Function ITS_TIME_FOR_EXIT returns TRUE if it is the end of the run 
!  (i.e. TAU >= TAUe) and false otherwise. (bmy, 3/21/03)
!
!  NOTES:
!******************************************************************************
!
      ! Function value
      LOGICAL :: FLAG 

      !=================================================================
      ! ITS_FOR_EXIT begins here!
      !=================================================================
      FLAG = ( TAU >= TAUe )

      ! Return to calling program
      END FUNCTION ITS_TIME_FOR_EXIT

!------------------------------------------------------------------------------
! Prior to 9/25/03:
! Now allow an argument to be passed to ITS_A_LEAPYEAR (bmy, 9/25/03)
!      FUNCTION ITS_A_LEAPYEAR() RESULT( IS_LEAPYEAR )
!!
!!******************************************************************************
!!  Function ITS_A_LEAPYEAR returns TRUE if YEAR is a leapyear, 
!!  and FALSE otherwise. (bmy, 3/17/99, 3/21/03)
!!
!!  NOTES: 
!!  (1 ) Now remove YEAR from ARG list; use the module variable (bmy, 3/21/03)
!!******************************************************************************
!!
!      ! Function value
!      LOGICAL :: IS_LEAPYEAR
!
!      !=================================================================
!      ! LEAPYEAR begins here!
!      !
!      ! A leap year is:
!      ! (1) evenly divisible by 4 (if not a century year)
!      ! (2) evenly divisible by 4, 100, and 400 (if a century year)
!      !
!      ! EXAMPLES:
!      ! (a) 1992 is a leap year since it is evenly divisible by 4, and 
!      !     is not a century year (i.e. it doesn't end in '00').
!      !
!      ! (b) 1900 is NOT a leap year, since while being evenly divisible 
!      !     by 4 and 100, it is NOT divisible by 400.
!      !
!      ! (c) 2000 is a leap year, since it is divisible by 4, 100, and 
!      !     400.
!      !=================================================================
!      IS_LEAPYEAR = .FALSE.
!
!      IF ( MOD( YEAR, 4 ) == 0 ) THEN
!         IF ( MOD( YEAR, 100 ) == 0 ) THEN
!            IF ( MOD( YEAR, 400 ) == 0 ) THEN
!               IS_LEAPYEAR = .TRUE.
!            ENDIF
!         ELSE
!            IS_LEAPYEAR = .TRUE.
!         ENDIF
!      ENDIF        
!
!      ! Return to calling program
!      END FUNCTION ITS_A_LEAPYEAR

      FUNCTION ITS_A_LEAPYEAR( YEAR_IN ) RESULT( IS_LEAPYEAR )
!
!******************************************************************************
!  Function ITS_A_LEAPYEAR tests to see if a year is really a leapyear. 
!  (bmy, 3/17/99, 9/25/03)
!
!  NOTES: 
!  (1 ) Now remove YEAR from ARG list; use the module variable (bmy, 3/21/03)
!  (2 ) Now add YEAR_IN as an optional argument.  If YEAR_IN is not passed,
!        then test if the current year is a leapyear (bmy, 9/25/03)
!******************************************************************************
!
      ! Arguments
      INTEGER, INTENT(IN), OPTIONAL :: YEAR_IN
      
      ! Local variables
      INTEGER                       :: THISYEAR

      ! Function value
      LOGICAL                       :: IS_LEAPYEAR

      !=================================================================
      ! LEAPYEAR begins here!
      !=================================================================

      ! If YEAR_IN is passed, use that value; otherwise use the value 
      ! of the current year as stored in module variable YEAR.
      IF ( PRESENT( YEAR_IN ) ) THEN
         THISYEAR = YEAR_IN
      ELSE
         THISYEAR = YEAR
      ENDIF

      !=================================================================
      ! A leap year is:
      ! (1) evenly divisible by 4 (if not a century year)
      ! (2) evenly divisible by 4, 100, and 400 (if a century year)
      !
      ! EXAMPLES:
      ! (a) 1992 is a leap year since it is evenly divisible by 4, and 
      !     is not a century year (i.e. it doesn't end in '00').
      !
      ! (b) 1900 is NOT a leap year, since while being evenly divisible 
      !     by 4 and 100, it is NOT divisible by 400.
      !
      ! (c) 2000 is a leap year, since it is divisible by 4, 100, and 
      !     400.
      !=================================================================
      IS_LEAPYEAR = .FALSE.

      IF ( MOD( THISYEAR, 4 ) == 0 ) THEN
         IF ( MOD( THISYEAR, 100 ) == 0 ) THEN
            IF ( MOD( THISYEAR, 400 ) == 0 ) THEN
               IS_LEAPYEAR = .TRUE.
            ENDIF
         ELSE
            IS_LEAPYEAR = .TRUE.
         ENDIF
      ENDIF        

      ! Return to calling program
      END FUNCTION ITS_A_LEAPYEAR

!------------------------------------------------------------------------------

      SUBROUTINE PRINT_CURRENT_TIME
!
!******************************************************************************
!  Subroutine PRINT_CURRENT_TIME prints the date, GMT time, and elapsed
!  hours of a GEOS-CHEM simulation. (bmy, 3/21/03)
!
!  NOTES:
!******************************************************************************
!
      ! Local variables
      REAL*4 :: E_HOURS

      !=================================================================
      ! PRINT_CURRENT_TIME begins here!
      !=================================================================

      ! Hours since start of run
      E_HOURS = REAL( ELAPSED_MIN ) / 60e0 

      ! Write quantities
      WRITE( 6, 100 ) YEAR, MONTH, DAY, HOUR, MINUTE, E_HOURS

      ! Format string
 100  FORMAT( '---> DATE: ', i4.4, '/', i2.2, '/', i2.2, 
     &            '  GMT: ', i2.2, ':', i2.2, '  X-HRS: ', f11.3 )

      ! Return to calling program
      END SUBROUTINE PRINT_CURRENT_TIME

!------------------------------------------------------------------------------

      FUNCTION NYMD_Y2K( NYMD ) RESULT( NYMD_TMP )
!
!******************************************************************************
!  NYMD_Y2K returns the proper form of NYMD for GEOS-1, GEOS-STRAT, 
!  GEOS-2, or GEOS-3 data. (bmy, 6/21/00, 11/21/01)
!                                                                          
!  Arguments as Input:
!  ============================================================================
!  (1 ) NYMD : 8 digit value for the current Year/month/date (YYYYMMDD)
!
!  NOTES:
!  (1 ) GEOS-1     runs from 1985 to 1995, so return a 6-digit YYMMDD value.
!  (2 ) GEOS-STRAT runs from 1995 to 1997, so return a 6-digit YYMMDD value.
!  (3 ) GEOS-2     must be Y2K compliant, so return an 8-digit YYYYMMDD value.
!  (4 ) GEOS-3     must be Y2K compliant, so return an 8-digit YYYYMMDD value.
!  (5 ) Updated comments (bmy, 9/4/01)
!  (6 ) Use #else condition for met fields other than GEOS-1 or GEOS-STRAT.
!        (bmy, 11/21/01)
!******************************************************************************
!
      ! Reference C-preprocessor switches
#     include "define.h"

      ! Arguments
      INTEGER, INTENT(IN) :: NYMD
      
      ! NYMD_TMP is the return value of the function
      INTEGER             :: NYMD_TMP

      !=================================================================
      ! NYMD_Y2K begins here!
      ! 
      ! Return the proper value of NYMD for the given data set.
      !=================================================================
#if   defined( GEOS_1 ) || defined( GEOS_STRAT )
      NYMD_TMP = NYMD - 19000000

#else
      NYMD_TMP = NYMD

#endif

      ! Return to calling program
      END FUNCTION NYMD_Y2K

!------------------------------------------------------------------------------

      FUNCTION NYMD6_2_NYMD8( NYMD ) RESULT( NYMD_TMP )
!
!******************************************************************************
!  Function NYMD6_2_NYMD8 converts a 6-digit YYMMDD string to a 
!  Y2K-compliant 8-digit YYYYMMDD string.  (bmy, 7/17/00, 9/4/01)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) NYMD : 6 digit value for the current Year/month/date (YYYYMMDD)
!
!  NOTES:
!  (1 ) NYMD6_2_NYMD8 operates regardless of which data set we are using.
!  (2 ) Small bug fix -- replace <= by >= in middle IF clause (bmy, 8/29/00)
!  (3 ) Updated comments (bmy, 9/4/01)
!******************************************************************************
!
      ! Arguments
      INTEGER, INTENT(IN) :: NYMD

      ! Return value of the function
      INTEGER             :: NYMD_TMP

      !=================================================================
      ! NYMD_Y2K begins here!
      ! 
      ! If NYMD is lies between 850101 and 991231, then this denotes
      ! the century 1900.  Otherwise, we are in the century 2000.  
      !
      ! If NYMD is already an 8 digit number, return it as is.
      !=================================================================
      IF ( NYMD >= 850101 .and. NYMD <= 991231 ) THEN
         NYMD_TMP = NYMD + 19000000

      ELSE IF ( NYMD >= 000000 .and. NYMD < 850101 ) THEN
         NYMD_TMP = NYMD + 20000000 

      ELSE IF ( NYMD > 19000000 ) THEN
         NYMD_TMP = NYMD

      ENDIF
      
      ! Return to calling program
      END FUNCTION NYMD6_2_NYMD8

!------------------------------------------------------------------------------

      FUNCTION NYMD_STRING( NYMD ) RESULT( NYMD_STR )
!
!******************************************************************************
!  NYMD_STRING returns a string variable of NYMD for GEOS-1, 
!  GEOS-STRAT, GEOS-2, or GEOS-3 data. (bmy, 6/21/00)
!                                                                          
!  Arguments as Input:
!  ============================================================================
!  (1 ) NYMD : 8 digit value for the current Year/month/date (YYYYMMDD)
!
!  NOTES:
!  (1 ) GEOS-1     runs from 1985 to 1995, so return a 6-char YYMMDD string.
!  (2 ) GEOS-STRAT runs from 1995 to 1997, so return a 6-char YYMMDD string.
!  (3 ) GEOS-2     must be Y2K compliant, so return an 8-char YYYYMMDD string.
!  (4 ) GEOS-3     must be Y2K compliant, so return an 8-char YYYYMMDD string.
!  (5 ) Updated comments (bmy, 9/4/01)
!  (6 ) Use #else condition for met fields other than GEOS-1 or GEOS-STRAT.
!        (bmy, 11/21/01)
!******************************************************************************
!     
      ! Reference C-preprocessor switches
#     include "define.h"
      
      ! Arguments
      INTEGER, INTENT(IN) :: NYMD

      ! Local Variables
      INTEGER             :: NYMD_TMP

      ! NYMD_STR holds the return value of the function
#if   defined( GEOS_1 ) || defined( GEOS_STRAT )
      CHARACTER(LEN=6)    :: NYMD_STR

#else
      CHARACTER(LEN=8)    :: NYMD_STR

#endif

      !=================================================================
      ! NYMD_STRING begins here!
      ! 
      ! Call NYMD_TMP to return the proper form of NYMD
      !
      ! Write NYMD_TMP to string variable NYMD_STR, which will be
      ! either 6 characters or 8 characters long, depending on 
      ! the data set being used (see above).
      !=================================================================      
      NYMD_TMP = NYMD_Y2K( NYMD )

#if   defined( GEOS_1 ) || defined( GEOS_STRAT )
      WRITE ( NYMD_STR, '(i6)' ) NYMD_TMP

#else
      WRITE ( NYMD_STR, '(i8)' ) NYMD_TMP

#endif

      ! Return to calling program
      END FUNCTION NYMD_STRING

!------------------------------------------------------------------------------

      FUNCTION DATE_STRING() RESULT( DATE_STR )
!
!******************************************************************************
!  DATE_STRING returns a string variable for the current date: YYMMDD 
!  for GEOS-1 or GEOS-STRAT and YYYYMMDD for GEOS-3 or GEOS-4. 
!  (bmy, 2/5/03, 9/29/03)
!                                                                          
!  Arguments as Input:
!  ============================================================================
!  (1 ) NYMD : 8 digit value for the current Year/month/date (YYYYMMDD)
!
!  NOTES:
!  (1 ) Adapted from routine NYMD_STRING (bmy, 2/5/03)
!  (2 ) Bug fix for GEOS-1/GEOS-S: print YEAR-1900 for 2 digits (bmy, 5/15/03)
!  (3 ) Bug fix for Linux: use ENCODE statement to convert number to string 
!        instead of a F90 internal read. (bmy, 9/29/03)
!******************************************************************************
!     
      ! Reference C-preprocessor switches
#     include "define.h"
      
      !=================================================================
      ! DATE_STRING begins here!
      !=================================================================
#if   defined( GEOS_1 ) || defined( GEOS_STRAT )

      !-----------------------
      ! Write YYMMDD format
      !-----------------------
      CHARACTER(LEN=6) :: DATE_STR
      
      ! Use ENCODE statement for PGI/F90 on Linux, or a regular 
      ! F90 internal read for other platforms (bmy, 9/29/03)
#if   defined( LINUX ) 
      ENCODE( 6, '(3i2.2)', DATE_STR ) YEAR-1900, MONTH, DAY
#else
      WRITE( DATE_STR, '(3i2.2)' ) YEAR-1900, MONTH, DAY
#endif

#else

      !-----------------------
      ! Write YYYYMMDD format
      !-----------------------
      CHARACTER(LEN=8) :: DATE_STR

      ! Use ENCODE statement for PGI/F90 on Linux, or a regular 
      ! F90 internal read for other platforms (bmy, 9/29/03)
#if   defined( LINUX ) 
      ENCODE( 8, '(i4.4,2i2.2)', DATE_STR ) YEAR, MONTH, DAY
#else
      WRITE( DATE_STR, '(i4.4,2i2.2)' ) YEAR, MONTH, DAY
#endif

#endif

      ! Return to calling program
      END FUNCTION DATE_STRING

!------------------------------------------------------------------------------

      FUNCTION TIMESTAMP_STRING() RESULT( TIME_STR )
!
!******************************************************************************
!  TIMESTAMP_STRING returns a formatted string "YYYY/MM/DD HH:MM" for the 
!  current date.  (bmy, 3/21/03, 9/29/03)
!                                                                          
!  NOTES:
!  (1 ) Now use ENCODE statement for PGI/F90 on Linux (bmy, 9/29/03)
!******************************************************************************
!     
#     include "define.h"

      ! Function value
      CHARACTER(LEN=16) :: TIME_STR
      
      !=================================================================
      ! TIMESTAMP_STRING begins here!
      !=================================================================

#if   defined( LINUX ) 
      
      ! For PGI/F90 Linux, we must use the ENCODE command
      ! to convert from numeric to string format (bmy, 9/29/03)
      ENCODE( 16, 100, TIME_STR ) YEAR, MONTH, DAY, HOUR, MINUTE

#else

      ! For other platforms, we can just use an internal read
      WRITE( TIME_STR, 100 ) YEAR, MONTH, DAY, HOUR, MINUTE

#endif

      ! Format statement
 100  FORMAT( i4.4, '/', i2.2, '/', i2.2, ' ', i2.2, ':', i2.2 )

      ! Return to calling program
      END FUNCTION TIMESTAMP_STRING

!------------------------------------------------------------------------------

      SUBROUTINE YMD_EXTRACT( NYMD, Y, M, D )
!
!******************************************************************************
!  Subroutine YMD_EXTRACT extracts the year, month, and date from an integer
!  variable in YYYYMMDD format.  It can also extract the hours, minutes, and
!  seconds from a variable in HHMMSS format. (bmy, 11/21/01)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) NYMD (INTEGER) : Variable in YYYYMMDD (or HHMMSS) format
!  
!  Arguments as Output:
!  ============================================================================
!  (2 ) Y    (INTEGER) : Variable that returns YYYY (or HH - hours  )
!  (3 ) M    (INTEGER) : Variable that returns MM   (or MM - minutes)
!  (4 ) D    (INTEGER) : Variable that returns DD   (or SS - seconds)
!
!  NOTES:
!******************************************************************************
!
      ! Arguments
      INTEGER, INTENT(IN)  :: NYMD
      INTEGER, INTENT(OUT) :: Y, M, D

      ! Local variables
      REAL*8               :: REM

      !=================================================================
      ! YMD_EXTRACT begins here!
      !=================================================================

      ! Extract YYYY from YYYYMMDD 
      Y = INT( DBLE( NYMD ) / 1d4 )

      ! Extract MM from YYYYMMDD
      REM = DBLE( NYMD ) - ( DBLE( Y ) * 1d4 )
      M   = INT( REM / 1d2 )

      ! Extract DD from YYYYMMDD
      REM = REM - ( DBLE( M ) * 1d2 )
      D   = INT( REM )

      ! Return to calling program
      END SUBROUTINE YMD_EXTRACT

!------------------------------------------------------------------------------

      SUBROUTINE EXPAND_DATE( FILENAME, YYYYMMDD, HHMMSS )
!
!******************************************************************************
!  Subroutine EXPAND_DATE replaces "YYYYMMDD" and "hhmmss" tokens within
!  a filename string with the actual values. (bmy, 6/27/02, 9/29/03)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) FILENAME (CHARACTER) : Filename with tokens to replace
!  (2 ) YYYYMMDD (INTEGER  ) : Current Year-Month-Day (must have 8 digits!)
!  (3 ) HHMMSS   (INTEGER  ) : Current Hour-Minute-Seconds
!  
!
!  Arguments as Output:
!  ============================================================================
!  (1 ) FILENAME (CHARACTER) : Modified filename 
!
!  NOTES:
!  (1 ) Bug fix for Linux: use ENCODE statement to convert number to string 
!        instead of F90 internal read. (bmy, 9/29/03)
!******************************************************************************
!      
      ! References to F90 modules
      USE CHARPAK_MOD, ONLY : STRREPL

#     include "define.h"

      ! Arguments
      CHARACTER(LEN=*), INTENT(INOUT) :: FILENAME
      INTEGER,          INTENT(IN)    :: YYYYMMDD, HHMMSS

      ! Local variables
      INTEGER                         :: YYYY, MM, DD, HH, II, SS
      CHARACTER(LEN=2)                :: MM_STR, DD_STR
      CHARACTER(LEN=2)                :: HH_STR, II_STR, SS_STR

#if   defined( GEOS_1 ) || defined( GEOS_STRAT )
      CHARACTER(LEN=2)                :: YYYY_STR
#else
      CHARACTER(LEN=4)                :: YYYY_STR
#endif

      !=================================================================
      ! EXPAND_DATE begins here!
      !=================================================================

      ! Extract today's date into year, month, and day sections
      CALL YMD_EXTRACT( NYMD_Y2K( YYYYMMDD ), YYYY, MM, DD )

      ! Extract today's time into HH, MM, and SS sections
      ! (rename minutes to II so as not to overwrite MM)
      CALL YMD_EXTRACT( HHMMSS, HH, II, SS )

      ! Convert integer values to char strings
#if   defined( GEOS_1 ) || defined( GEOS_STRAT )

      ! Use ENCODE statement for PGI/Linux (bmy, 9/29/03)
#if   defined( LINUX )
      ENCODE( 2, '(i2.2)', YYYY_STR ), YYYY
#else
      WRITE( YYYY_STR, '(i2.2)' ) YYYY
#endif

#else
    
      ! Use ENCODE statement for PGI/Linux (bmy, 9/29/03)
#if   defined( LINUX ) 
      ENCODE( 4, '(i4.4)', YYYY_STR ), YYYY
#else      
      WRITE( YYYY_STR, '(i4.4)' ) YYYY
#endif

#endif

#if   defined( LINUX )
      
      ! Use ENCODE statement for PGI/Linux (bmy, 9/29/03)
      ENCODE( 2, '(i2.2)', MM_STR ) MM
      ENCODE( 2, '(i2.2)', DD_STR ) DD
      ENCODE( 2, '(i2.2)', HH_STR ) HH
      ENCODE( 2, '(i2.2)', II_STR ) II
      ENCODE( 2, '(i2.2)', SS_STR ) SS
#else

      ! For other platforms, use an F90 internal read (bmy, 9/29/03)
      WRITE( MM_STR, '(i2.2)' ) MM
      WRITE( DD_STR, '(i2.2)' ) DD
      WRITE( HH_STR, '(i2.2)' ) HH
      WRITE( II_STR, '(i2.2)' ) II
      WRITE( SS_STR, '(i2.2)' ) SS

#endif

      ! Replace YYYY, MM, DD, HH tokens w/ actual values 
#if   defined( GEOS_1 ) || defined( GEOS_STRAT )
      CALL STRREPL( FILENAME, 'YY',   YYYY_STR )
#else
      CALL STRREPL( FILENAME, 'YYYY', YYYY_STR )
#endif
      CALL STRREPL( FILENAME, 'MM',   MM_STR   )
      CALL STRREPL( FILENAME, 'DD',   DD_STR   )
      CALL STRREPL( FILENAME, 'hh',   HH_STR   )
      CALL STRREPL( FILENAME, 'mm',   II_STR   )
      CALL STRREPL( FILENAME, 'ss',   SS_STR   )

      ! Return to calling program
      END SUBROUTINE EXPAND_DATE

!------------------------------------------------------------------------------

      END MODULE TIME_MOD




