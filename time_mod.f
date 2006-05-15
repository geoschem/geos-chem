! $Id: time_mod.f,v 1.21 2006/05/15 17:52:55 bmy Exp $
      MODULE TIME_MOD
!
!******************************************************************************
!  TIME_MOD contains GEOS-CHEM date and time variables and timesteps, and 
!  routines for accessing them. (bmy, 6/21/00, 4/24/06) 
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
!  (33) NDIAGTIME   (INTEGER) : Time of day (HHMMSS) to write bpch file
!
!  Module Routines:
!  ============================================================================
!  (1 ) SET_CURRENT_TIME  : Updates time variables for current timestep
!  (2 ) SET_BEGIN_TIME    : Initializes NYMDb, NHMSb, TAUb variables
!  (3 ) SET_END_TIME      : Initializes NYMDe, NHMSe, TAUe variables
!  (4 ) SET_NDIAGTIME     : Initializes NDIAGTIME (time to write bpch file)
!  (5 ) SET_DIAGb         : Updates DIAGb and DIAGe diagnostic interval times
!  (6 ) SET_DIAGe         : Updates DIAGb and DIAGe diagnostic interval times
!  (7 ) SET_TIMESTEPS     : Updates the elapsed minutes since the start of run
!  (8 ) SET_CT_CHEM       : Increments/resets the chemistry timestep counter
!  (9 ) SET_CT_CONV       : Increments/resets the convection timestep counter
!  (10) SET_CT_DYN        : Increments/resets the dynamic timestep counter
!  (11) SET_CT_EMIS       : Increments/resets the emissions timestep counter
!  (12) SET_CT_A3         : Increments/resets the A-3 fields timestep counter
!  (13) SET_CT_A6         : Increments/resets the A-6 fields timestep counte
!  (14) SET_CT_I6         : Increments/resets the I-6 fields timestep counter
!  (15) SET_CT_XTRA      : Increments/resets the I-6 fields timestep counter
!  (16) SET_ELAPSED_MIN   : Updates the elapsed minutes since the start of run
!  (17) GET_JD            : Returns Astronomical Julian Date for NYMD, NHMS
!  (18) GET_ELAPSED_MIN   : Returns the elapsed minutes since the start of run
!  (19) GET_ELAPSED_SEC   : Returns the elapsed seconds since the start of run 
!  (20) GET_NYMDb         : Returns the YYYYMMDD at the beginning of the run
!  (21) GET_NHMSb         : Returns the HHMMSS   at the beginning of the run
!  (22) GET_NYMDe         : Returns the YYYYMMDD at the end of the run
!  (23) GET_NHMSe         : Returns the HHMMSS   at the end of the run
!  (24) GET_NYMD          : Returns the YYYYMMDD at the current time
!  (25) GET_NHMS          : Returns the HHMMSS   at the current time
!  (26) GET_NDIAGTIME     : Returns NDIAGTIME (time of day to write bpch file)
!  (27) GET_TIME_AHEAD    : Returns the YYYYMMDD, HHMMSS for N_MINS from now
!  (28) GET_MONTH         : Returns the current month (1-12)
!  (29) GET_DAY           : Returns the current day of month (1-31)
!  (30) GET_YEAR          : Returns the current year (YYYY)
!  (31) GET_HOUR          : Returns the current hour (0-23)
!  (32) GET_MINUTE        : Returns the current minute (0-59)
!  (33) GET_SECOND        : Returns the current second (0-59)
!  (34) GET_DAY_OF_YEAR   : Returns the current day of the year (0-366)
!  (35) GET_DAY_OF_WEEK   : Returns the current day of the week (0-6)
!  (36) GET_GMT           : Returns the current GMT (0.0 - 23.999)
!  (37) GET_TAU           : Returns the current TAU value (hrs since 1/1/1985)
!  (38) GET_TAUb          : Returns TAU value at beginning of GEOS-CHEM run
!  (39) GET_TAUe          : Returns TAU value at end of GEOS-CHEM run
!  (40) GET_DIAGb         : Returns TAU value at start of diagnostic interval
!  (41) GET_DIAGe         : Returns TAU value at end of diagnostic interval
!  (42) GET_LOCALTIME     : Returns local time for a grid box (0.0 - 23.999)
!  (43) GET_SEASON        : Returns season flag (1=DJF, 2=MAM, 3=JJA, 4=SON)
!  (44) GET_TS_CHEM       : Returns chemistry timestep in minutes
!  (45) GET_TS_CONV       : Returns convection timestep in minutes
!  (46) GET_TS_DIAG       : Returns diagnostic timestep in minutes
!  (47) GET_TS_DYN        : Returns dynamic timestep in minutes
!  (48) GET_TS_EMIS       : Returns emissions timestep in minutes
!  (49) GET_TS_UNIT       : Returns unit conversion timestep in minutes
!  (50) GET_CT_CHEM       : Returns # of chemistry timesteps already executed
!  (51) GET_CT_CONV       : Returns # of convection timesteps already executed
!  (52) GET_CT_DYN        : Returns # of dynamic timesteps already executed
!  (53) GET_CT_EMIS       : Returns # of emission timesteps already executed
!  (54) GET_CT_A3         : Returns # of times A-3 fields have been read in
!  (55) GET_CT_A6         : Returns # of times A-6 fields have been read in
!  (56) GET_CT_I6         : Returns # of times I-6 fields have been read in
!  (57) GET_CT_XTRA       : Returns # of times I-6 fields have been read in
!  (58) GET_A3_TIME       : Returns YYYYMMDD and HHMMSS for the A-3 fields
!  (59) GET_A6_TIME       : Returns YYYYMMDD and HHMMSS for the A-6 fields
!  (60) GET_I6_TIME       : Returns YYYYMMDD and HHMMSS for the I-6 fields
!  (61) GET_FIRST_A3_TIME : Returns YYYYMMDD and HHMMSS for the first A-3 read
!  (62) GET_FIRST_A3_TIME : Returns YYYYMMDD and HHMMSS for the first A-6 read
!  (63) ITS_TIME_FOR_CHEM : Returns TRUE if it is time to do chemistry
!  (64) ITS_TIME_FOR_CONV : Returns TRUE if it is time to do convection
!  (65) ITS_TIME_FOR_DYN  : Returns TRUE if it is time to do dynamics
!  (66) ITS_TIME_FOR_EMIS : Returns TRUE if it is time to do emissions
!  (67) ITS_TIME_FOR_UNIT : Returns TRUE if it is time to do unit conversions
!  (68) ITS_TIME_FOR_DIAG : Returns TRUE if it is time to write diagnostics
!  (69) ITS_TIME_FOR_A3   : Returns TRUE if it is time to read in A-3 fields
!  (71) ITS_TIME_FOR_A6   : Returns TRUE if it is time to read in A-6 fields
!  (72) ITS_TIME_FOR_I6   : Returns TRUE if it is time to read in I-6 fields
!  (73) ITS_TIME_FOR_UNZIP: Returns TRUE if it is the end of the run
!  (74) ITS_TIME_FOR_DEL  : Returns TRUE if it is time to delete temp files
!  (75) ITS_TIME_FOR_EXIT : Returns TRUE if it is the end of the run
!  (76) ITS_A_LEAPYEAR    : Returns TRUE if the current year is a leapyear
!  (77) ITS_A_NEW_YEAR    : Returns TRUE if it's a new year
!  (78) ITS_A_NEW_MONTH   : Returns TRUE if it's a new month
!  (79) ITS_MIDMONTH      : Returns TRUE if it's 0 GMT on the 16th of the month
!  (80) ITS_A_NEW_DAY     : Returns TRUE if it's a new day
!  (81) ITS_A_NEW_SEASON  : Returns TRUE if it's a new season
!  (82) TIMESTAMP_STRING  : Returns a string "YYYY/MM/DD HH:MM:SS"
!  (83) PRINT_CURRENT_TIME: Prints date time in YYYY/MM/DD, HH:MM:SS format
!  (84) YMD_EXTRACT       : Extracts YYYY, MM, DD from a YYYYMMDD format number
!  (85) EXPAND_DATE       : Replaces date/time tokens w/ actual values
!  (86) SYSTEM_DATE_TIME  : Returns the system date and time
!  (87) SYSTEM_TIMESTAMP  : Returns a string with the system date and time
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
!  (11) Bug fix in EXPAND_DATE.  Also add optional arguments to function
!        TIMESTAMP_STRNIG. (bmy, 10/28/03)
!  (12) Changed the name of some cpp switches in "define.h" (bmy, 12/2/03)
!  (13) Modified ITS_TIME_FOR_A6 and GET_FIRST_A6_TIME for both GEOS-4 
!        "a_llk_03" and "a_llk_04" data versions. (bmy, 3/22/04)
!  (14) Added routines ITS_A_NEW_MONTH, ITS_A_NEW_YEAR, ITS_A_NEW_DAY.
!        (bmy, 4/1/04)
!  (15) Added routines ITS_A_NEW_SEASON, GET_NDIAGTIME, SET_NDIAGTIME, and
!        variable NDIAGTIME. (bmy, 7/20/04)
!  (17) Added routine GET_DAY_OF_WEEK (bmy, 11/5/04)
!  (18) Removed obsolete FIRST variable in GET_A3_TIME (bmy, 12/10/04)
!  (19) Added routines SYSTEM_DATE_TIME and SYSTEM_TIMESTAMP.  Also modified
!        for GCAP and GEOS-5 met fields. (swu, bmy, 5/3/05)
!  (20) GCAP/GISS met fields don't have leap years (swu, bmy, 8/29/05)
!  (21) Added counter variable & routines for XTRA fields (tmf, bmy, 10/20/05)
!  (22) Bug fix in ITS_A_NEW_YEAR (bmy, 11/1/05)
!  (23) Added function ITS_MIDMONTH.  Also removed obsolete functions
!        NYMD_Y2K, NYMD6_2_NYMD8, NYMD_STRING, DATE_STRING. 
!        (sas, cdh, bmy, 12/15/05)
!  (24) GCAP bug fix: There are no leapyears, so transition from 2/28 to 3/1,
!        skipping 2/29 for all years. (swu, bmy, 4/24/06)
!******************************************************************************
!
      IMPLICIT NONE

      !=================================================================
      ! MODULE PRIVATE DECLARATIONS -- keep certain internal variables 
      ! and routines from being seen outside "time_mod.f"
      !=================================================================

      ! Make everything PUBLIC ...
      PUBLIC

      ! ... except these variables
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
      PRIVATE           :: CT_I6,      CT_XTRA,     JD85
      PRIVATE           :: NDIAGTIME

      !=================================================================
      ! MODULE VARIABLES
      !=================================================================

      ! Date and time variables
      INTEGER           :: NYMDb,      NHMSb,       NYMDe   
      INTEGER           :: NHMSe,      NYMD,        NHMS
      INTEGER           :: MONTH,      DAY,         YEAR
      INTEGER           :: HOUR,       MINUTE,      SECOND
      INTEGER           :: NSEASON,    DAY_OF_YEAR, ELAPSED_MIN
      INTEGER           :: NDIAGTIME
      REAL*8            :: TAU,        TAUb,        TAUe  
      REAL*8            :: GMT,        DIAGb,       DIAGe

      ! Timesteps
      INTEGER           :: TS_CHEM,   TS_CONV,     TS_DIAG
      INTEGER           :: TS_DYN,    TS_EMIS,     TS_UNIT

      ! Timestep counters
      INTEGER           :: CT_CHEM,   CT_CONV,     CT_DYN    
      INTEGER           :: CT_EMIS,   CT_A3,       CT_A6
      INTEGER           :: CT_I6,     CT_XTRA

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
!  Subroutine SET_CURRENT_TIME takes in the elapsed time in minutes since the 
!  start of a GEOS-CHEM simulation and sets the GEOS-CHEM time variables 
!  accordingly. (bmy, 2/5/03, 4/24/06)
!
!  NOTES:
!  (1 ) GCAP/GISS fields don't have leap years, so if JULDAY says it's 
!        Feb 29th, reset MONTH, DAY, JD1 to Mar 1st. (swu, bmy, 8/29/05)
!  (2 ) Now references "define.h".  Now add special handling to skip from
!        Feb 28th to Mar 1st for GCAP model. (swu, bmy, 4/24/06)
!******************************************************************************
!
      ! References to F90 modules
      USE JULDAY_MOD, ONLY : JULDAY, CALDATE

#     include "define.h"

      ! Local variables
      LOGICAL :: IS_LEAPYEAR
      REAL*4  :: TMP
      REAL*8  :: JD0, JD1, JD_JAN_1
      
      !=================================================================
      ! SET_CURRENT_TIME begins here!
      !=================================================================

      ! JD0: Astronomical Julian Date at start of GEOS-CHEM run
      JD0 = GET_JD( NYMDb, NHMSb )

      ! JD1: Astronomical Julian Date at current time
      JD1 = JD0 + ( DBLE( ELAPSED_MIN ) / 1440d0 )

      ! Call CALDATE to compute the current YYYYMMDD and HHMMSS
      CALL CALDATE( JD1, NYMD, NHMS )

      ! Extract current year, month, day from NYMD
      CALL YMD_EXTRACT( NYMD, YEAR, MONTH, DAY )

#if   defined( GCAP ) 

      !-------------------------------
      ! GCAP met fields: no leapyears
      !-------------------------------

      ! Special handling for leap years 
      IF ( ITS_A_LEAPYEAR( YEAR, FORCE=.TRUE. ) ) THEN

         ! Get Astronomical Julian Date on Jan 0th of this year
         JD_JAN_1 = GET_JD( YEAR*10000 + 0101, 000000 )
         
         ! Skip directly from Feb 28 to Mar 1st 
         IF ( JD1 - JD_JAN_1 >= 59d0 ) JD1 = JD1 + 1d0
         
         ! Call CALDATE to recompute YYYYMMDD and HHMMSS
         CALL CALDATE( JD1, NYMD, NHMS )

         ! Extract current year, month, day from NYMD
         CALL YMD_EXTRACT( NYMD, YEAR, MONTH, DAY )
      ENDIF
         
#endif

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
!  date and time of a GEOS-CHEM run. (bmy, 2/5/03, 7/20/04)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) THISNYMDb (INTEGER) : YYYYMMDD at beginning of GEOS-CHEM run
!  (2 ) THISNHMSb (INTEGER) : HHMMSS   at beginning of GEOS-CHEM run
!
!  NOTES:
!  (1 ) Added error check for THISNHMSb (bmy, 7/20/04)
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

      ! Make sure NHMSb is valid
      IF ( THISNHMSb > 235959 ) THEN
         CALL ERROR_STOP( 'NHMSb cannot be greater than 23:59:59!',
     &                    'SET_BEGIN_TIME (time_mod.f)' )
      ENDIF

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
!  date and time of a GEOS-CHEM run. (bmy, 2/5/03, 7/20/04)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) THISNYMDe (INTEGER) : YYYYMMDD at end of GEOS-CHEM run
!  (2 ) THISNHMSe (INTEGER) : HHMMSS   at end of GEOS-CHEM run
!
!  NOTES:
!  (1 ) Added error check for THISNHMSb (bmy, 7/20/04)
!******************************************************************************
!
      ! References to F90 modules
      USE ERROR_MOD,  ONLY : ERROR_STOP

      ! Arguments
      INTEGER, INTENT(IN) :: THISNYMDe, THISNHMSe
      
      ! Local variables
      REAL*4              :: TMP

      !=================================================================
      ! SET_END_TIME begins here!
      !=================================================================

      ! Error check to make sure 
      IF ( THISNHMSe > 235959 ) THEN
         CALL ERROR_STOP( 'NHMSe cannot be greater than 23:59:59!',
     &                    'SET_END_TIME (time_mod.f)' )
      ENDIF

      ! Make sure THISNYMDb uses 4 digits for the year
      ! and is not less than 1985/01/01
      IF ( THISNYMDe < 19850101 ) THEN
         CALL ERROR_STOP( 'NYMDe must be in the format YYYYMMDD!',
     &                    'SET_END_TIME (time_mod.f)' )

      ENDIF

      ! Initialize NYMDe, NHMSe
      NYMDe = THISNYMDe
      NHMSe = THISNHMSe

      ! TAUe value (TMP is REAL*4 to prevent precision problems)
      TMP   = ( GET_JD( NYMDe, NHMSe ) - JD85 ) * 24e0
      TAUe  = DBLE( TMP )

      ! Return to calling program
      END SUBROUTINE SET_END_TIME

!------------------------------------------------------------------------------

      SUBROUTINE SET_NDIAGTIME( THIS_NDIAGTIME )
!
!******************************************************************************
!  Subroutine SET_NDIAGTIME initializes NDIAGTIME, the time of day at which
!  the binary punch file will be written out to disk. (bmy, 7/20/04)
! 
!  NOTES:
!******************************************************************************
!
      ! Arguments
      INTEGER, INTENT(IN) :: THIS_NDIAGTIME

      !=================================================================
      ! SET_NDIAGTIME begins here!
      !=================================================================
      NDIAGTIME = THIS_NDIAGTIME
      
      ! Return to calling program
      END SUBROUTINE SET_NDIAGTIME

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
!  chemistry, and emissions.  Counters are also zeroed. (bmy, 3/21/03,10/20/05)
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
!  (1 ) Suppress some output lines (bmy, 7/20/04)
!  (2 ) Also zero CT_XTRA (tmf, bmy, 10/20/05)
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
      CT_XTRA = 0

      ! Echo to stdout
      WRITE( 6, '(/,a)' ) 'SET_TIMESTEPS: setting GEOS-CHEM timesteps!'
      WRITE( 6, '(  a)' ) '-------------------------------------------'
      WRITE( 6, '(''Chemistry  Timestep [min] : '', i4 )' ) TS_CHEM
      WRITE( 6, '(''Convection Timestep [min] : '', i4 )' ) TS_CONV
      WRITE( 6, '(''Dynamics   Timestep [min] : '', i4 )' ) TS_DYN
      WRITE( 6, '(''Emission   Timestep [min] : '', i4 )' ) TS_EMIS
      WRITE( 6, '(''Unit Conv  Timestep [min] : '', i4 )' ) TS_UNIT

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

      SUBROUTINE SET_CT_XTRA( INCREMENT, RESET )
!
!******************************************************************************
!  Subroutine SET_CT_XTRA increments CT_XTRA, the counter of the number of 
!  times we have read in GEOS-3 XTRA fields.  (tmf, bmy, 10/20/05)
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
         CT_XTRA = CT_XTRA + 1
      ELSE IF ( PRESENT( RESET ) ) THEN
         CT_XTRA = 0
      ENDIF

      ! Return to calling program
      END SUBROUTINE SET_CT_XTRA

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

      FUNCTION GET_NDIAGTIME() RESULT( THIS_NDIAGTIME )
!
!******************************************************************************
!  Subroutine GET_NDIAGTIME returns to the calling program NDIAGTIME, the 
!  time of day at which the binary punch file will be written out to disk. 
!  (bmy, 7/20/04)
! 
!  NOTES:
!******************************************************************************
!
      ! Local variables
      INTEGER :: THIS_NDIAGTIME

      !=================================================================
      ! GET_NDIAGTIME begins here!
      !=================================================================
      THIS_NDIAGTIME = NDIAGTIME
      
      ! Return to calling program
      END FUNCTION GET_NDIAGTIME

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

      FUNCTION GET_DAY_OF_WEEK() RESULT( DAY_NUM )
!
!******************************************************************************
!  Function GET_DAY_OF_WEEK returns the day of the week as a number:
!  Sun=0, Mon=1, Tue=2, Wed=3, Thu=4, Fri=5, Sat=6.  (bmy, 11/5/04)
!
!  Reference:
!  ============================================================================
!  "Practical Astronomy with Your Calculator", 3rd Ed.  Peter Duffett-Smith,
!    Cambridge UP, 1992, p9.
!
!  NOTES:
!******************************************************************************
!
      ! References to F90 modules
      USE JULDAY_MOD, ONLY : JULDAY
      
      ! Return value
      INTEGER :: DAY_NUM

      ! Local variables
      REAL*8  :: A, B, JD, THISDAY

      !=================================================================
      ! GET_DAY_OF_WEEK begins here!
      !=================================================================

      ! Get fractional day
      THISDAY = DAY                 + ( HOUR   / 24d0    ) + 
     &          ( MINUTE / 1440d0 ) + ( SECOND / 86400d0 )

      ! Get current Julian date 
      JD      = JULDAY( YEAR, MONTH, THISDAY )

      ! Add 1.5 to JD and divide by 7
      A       = ( JD + 1.5d0 ) / 7d0 

      ! Take fractional part and multiply by 7
      B       = ( A - INT( A ) ) * 7d0

      ! Round to nearest integer -- this is the day number!
      DAY_NUM = INT( B + 0.5d0 )
      
      ! Return to calling program
      END FUNCTION GET_DAY_OF_WEEK

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

      FUNCTION GET_CT_XTRA() RESULT( THIS_CT_XTRA )
!
!******************************************************************************
!  Function GET_CT_XTRA returns the XTRA fields timestep counter to the
!  calling program. (tmf, bmy, 10/20/05)
!
!  NOTES:
!******************************************************************************
!
      ! Function value
      INTEGER :: THIS_CT_XTRA

      !=================================================================
      ! GET_CT_I6 begins here!
      !=================================================================
      THIS_CT_XTRA = CT_XTRA

      ! Return to calling program
      END FUNCTION GET_CT_XTRA

!------------------------------------------------------------------------------

      FUNCTION GET_A3_TIME() RESULT( DATE )
!
!******************************************************************************
!  Function GET_A3_TIME returns the correct YYYYMMDD and HHMMSS values
!  that are needed to read in the next average 3-hour (A-3) fields. 
!  (bmy, 3/21/03, 5/24/05)
!
!  NOTES:
!  (1 ) Now return proper time for GEOS-4/fvDAS fields (bmy, 6/19/03)
!  (2 ) Remove reference to FIRST variable (bmy, 12/10/04)
!  (3 ) Now modified for GCAP and GEOS-5 met fields (swu, bmy, 5/24/05)
!******************************************************************************
!
#     include "define.h"

      ! Function value
      INTEGER :: DATE(2)

      !=================================================================
      ! GET_A3_TIME begins here!
      !=================================================================

#if   defined( GEOS_1 ) || defined( GEOS_STRAT ) || defined( GEOS_3 )

      ! For GEOS-1, GEOS-STRAT, GEOS-3, the A-3 fields are timestamped 
      ! by ending time.  Therefore, the difference between the actual time
      ! when the fields are read and the A-3 timestamp time is 180 minutes.
      DATE = GET_TIME_AHEAD( 180 )      

#else

      ! For GEOS-4, GEOS-5, or GCAP data: The A-3 fields are timestamped 
      ! by center time.  Therefore, the difference between the actual time 
      ! when the fields are read and the A-3 timestamp time is 90 minutes.
      DATE = GET_TIME_AHEAD( 90 )   

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
!  (bmy, 3/21/03, 4/24/06)
!
!  NOTES:
!  (1 ) Bug fix for GCAP: skip over Feb 29th (no leapyears). (bmy, 4/24/06)
!******************************************************************************
!
      ! Arguments
      INTEGER :: DATE(2)

      !=================================================================
      ! GET_I6_TIME begins here!
      !=================================================================

#if   defined( GCAP ) 

      !-------------------------------
      ! GCAP met fields: no leapyears
      !-------------------------------

      ! If 18 GMT on Feb 28th, the next I-6 time is 0 GMT on Mar 1st
      IF ( MONTH == 2 .and. DAY == 28 .and. HOUR == 18 ) THEN
         DATE = (/ ( YEAR * 10000 ) + 0301, 000000 /)
         RETURN
      ENDIF

#endif

      !-------------------------------
      ! GEOS met fields: w/ leapyears
      !-------------------------------

      ! We need to read in the I-6 fields 6h (360 mins) ahead of time
      DATE = GET_TIME_AHEAD( 360 )

      ! Return to calling program
      END FUNCTION GET_I6_TIME

!------------------------------------------------------------------------------

      FUNCTION GET_FIRST_A3_TIME() RESULT( DATE )
!
!******************************************************************************
!  Function GET_FIRST_A3_TIME returns the correct YYYYMMDD and HHMMSS 
!  values the first time that A-3 fields are read in from disk. 
!  (bmy, 6/26/03, 5/24/05)
!
!  NOTES:
!  (1 ) Now modified for GCAP and GEOS-5 data (swu, bmy, 5/24/05)
!******************************************************************************
!
#     include "define.h"

      ! Arguments
      INTEGER :: DATE(2)

      !=================================================================
      ! GET_FIRST_A3_TIME begins here!
      !=================================================================

#if   defined( GEOS_1 ) || defined( GEOS_STRAT ) || defined( GEOS_3 )

      ! For GEOS-1, GEOS-STRAT, GEOS-3: Return the current date/time
      DATE = (/ NYMD, NHMS /)

#else
  
      ! For GEOS-4, GEOS-5, and GCAP: Call GET_A3_TIME to return 
      ! the date/time under which the A-3 fields are timestamped
      DATE = GET_A3_TIME()
    
#endif

      ! Return to calling program
      END FUNCTION GET_FIRST_A3_TIME

!------------------------------------------------------------------------------

      FUNCTION GET_FIRST_A6_TIME() RESULT( DATE )
!
!******************************************************************************
!  Function GET_FIRST_A6_TIME returns the correct YYYYMMDD & HHMMSS values the
!  first time that A-6 fields are read in from disk. (bmy, 6/26/03, 5/24/05)
!
!  NOTES:
!  (1 ) Now modified for GEOS-4 "a_llk_03" and "a_llk_04" fields (bmy, 3/22/04)
!  (2 ) Modified for GCAP and GEOS-5 met fields (swu, bmy, 5/24/05)
!******************************************************************************
!
#     include "define.h"

      ! Arguments
      INTEGER :: DATE(2)

      !=================================================================
      ! GET_FIRST_A6_TIME begins here!
      !=================================================================

#if   defined( GCAP )

      ! For GCAP data: Call GET_A6_TIME to return date/time
      ! under which the A-6 fields are timestamped
      DATE = GET_A6_TIME()      

#else

      ! For GEOS data: Return the current date/time
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
!  (average 6-h fields) and FALSE otherwise. (bmy, 3/21/03, 5/24/05)
!
!  NOTES:
!  (1 ) Now compute when it's time to read in GEOS-4 A-6 fields. (bmy, 6/26/03)
!  (2 ) Now modified for GEOS-4 "a_llk_03" and "a_llk_04" fields (bmy, 3/22/04)
!  (3 ) Now modified for GCAP and GEOS-5 met fields (swu, bmy, 5/24/05)
!******************************************************************************
!
#     include "define.h"

      ! Local variables
      INTEGER :: DATE(2)

      ! Function value
      LOGICAL :: FLAG 

      !=================================================================
      ! ITS_TIME_FOR_A6 begins here!
      !=================================================================

#if   defined( GCAP )

      ! For GCAP data: We need to read A-6 fields when it 00, 06, 
      ! 12, 18 GMT.  DATE is the current time -- test below.
      DATE = GET_TIME_AHEAD( 0 )

#else

      ! For all GEOS data: We need to read A-6 fields when it is 03, 
      ! 09, 15, 21 GMT.  DATE is the time 3 hours from now -- test below.
      DATE = GET_TIME_AHEAD( 180 )
     
#endif

      ! Test if DATE corresponds to 00, 06, 12, 18 GMT.  
      ! If so, then it is time to read A-6 fields from disk.
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

      ! We read in I-6 fields at 00, 06, 12, 18 GMT
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

      FUNCTION ITS_A_LEAPYEAR( YEAR_IN, FORCE ) RESULT( IS_LEAPYEAR )
!
!******************************************************************************
!  Function ITS_A_LEAPYEAR tests to see if a year is really a leapyear. 
!  (bmy, 3/17/99, 4/24/06)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) YEAR_IN (INTEGER) : (OPTIONAL) Specify a year to test for leapyear
!  (2 ) FORCE   (LOGICAL) : (OPTIONAL) Do not exit if using GCAP met fields  
!
!  NOTES: 
!  (1 ) Now remove YEAR from ARG list; use the module variable (bmy, 3/21/03)
!  (2 ) Now add YEAR_IN as an optional argument.  If YEAR_IN is not passed,
!        then test if the current year is a leapyear (bmy, 9/25/03)
!  (3 ) Now always return FALSE for GCAP (swu, bmy, 8/29/05)
!  (4 ) Now add FORCE argument to force ITS_A_LEAPYEAR to return a value
!        instead of just returning with FALSE for the GCAP met fields.
!        (swu, bmy, 4/24/06)
!******************************************************************************
!
      ! Arguments
      INTEGER, INTENT(IN), OPTIONAL :: YEAR_IN
      LOGICAL, INTENT(IN), OPTIONAL :: FORCE
      
      ! Local variables
      INTEGER                       :: THISYEAR
      LOGICAL                       :: THISFORCE

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

      ! If FORCE is passed, use that value, otherwise default to .FALSE.
      IF ( PRESENT( FORCE ) ) THEN
         THISFORCE = FORCE
      ELSE
         THISFORCE = .FALSE.
      ENDIF

      !=================================================================
      ! A leap year is:
      ! (1) evenly divisible by 4 (if not a century year)
      ! (2) evenly divisible by 4, 100, and 400 (if a century year)
      !
      ! EXAMPLES:
      ! (a) 1992 is a leap year since it is evenly divisible by 4, 
      !     and is not a century year (i.e. it doesn't end in '00').
      !
      ! (b) 1900 is NOT a leap year, since while being evenly divisible 
      !     by 4 and 100, it is NOT divisible by 400.
      !
      ! (c) 2000 is a leap year, since it is divisible by 
      !     4, 100, and 400.
      !=================================================================
      IS_LEAPYEAR = .FALSE.

#if   defined( GCAP )
      !-----------------------------------------------------------------------
      ! Prior to 4/24/06:
      ! For GCAP/GISS met fields, there are no leap years (swu, bmy, 8/29/05)
      !RETURN
      !-----------------------------------------------------------------------

      ! For GCAP met fields, there are no leap years.  However, sometimes
      ! we need to test to see if it would be a leap year so that we can
      ! tell the GEOS-Chem timing functions to skip past Feb 29th.  If 
      ! argument FORCE=T, then return the value of IS_LEAPYEAR to the
      ! calling program (bmy, 4/24/06)
      IF ( .not. THISFORCE ) RETURN

#endif

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

      FUNCTION ITS_A_NEW_YEAR() RESULT( IS_NEW_YEAR )
!
!******************************************************************************
!  Function ITS_A_NEW_YEAR returns TRUE if it's the first of a new month
!  (it also returns TRUE on the first timestep of the run).  This is useful
!  for setting flags for reading in data. (bmy, 4/1/04)
!
!  NOTES:
!  (1 ) Bug fix: Need month & day to be 1 (bmy, 11/1/05)
!******************************************************************************
!
      ! Function value
      LOGICAL :: IS_NEW_YEAR
      
      !=================================================================
      ! ITS_A_NEW_YEAR begins here!
      !=================================================================
      IF ( MONTH == 1 .and. DAY == 1 .and. NHMS == 000000 ) THEN

         ! A new year is Jan 1 at 0 GMT
         IS_NEW_YEAR = .TRUE.

      ELSE IF ( NYMD == NYMDb .and. NHMS == NHMSb ) THEN

         ! Also return TRUE if it's the start of the run
         ! (since files will need to be read in from disk)
         IS_NEW_YEAR = .TRUE.

      ELSE

         ! Otherwise, it's not a new year
         IS_NEW_YEAR = .FALSE.
         
      ENDIF

      ! Return to calling program
      END FUNCTION ITS_A_NEW_YEAR

!------------------------------------------------------------------------------

      FUNCTION ITS_A_NEW_MONTH() RESULT( IS_NEW_MONTH )
!
!******************************************************************************
!  Function ITS_A_NEW_MONTH returns TRUE if it's the first of a new month
!  (it also returns TRUE on the first timestep of the run).  This is useful
!  for setting flags for reading in data. (bmy, 4/1/04)
!
!  NOTES:
!******************************************************************************
!
      ! Function value
      LOGICAL :: IS_NEW_MONTH
      
      !=================================================================
      ! ITS_A_NEW_MONTH begins here!
      !=================================================================
      IF ( DAY == 1 .and. NHMS == 000000 ) THEN

         ! Test for the 1st of the month at 0 GMT
         IS_NEW_MONTH = .TRUE.

      ELSE IF ( NYMD == NYMDb .and. NHMS == NHMSb ) THEN

         ! Also return TRUE if it's the start of the run
         ! (since files will need to be read in from disk)
         IS_NEW_MONTH = .TRUE.

      ELSE

         ! Otherwise, it's not a new year
         IS_NEW_MONTH = .FALSE.
         
      ENDIF

      ! Return to calling program
      END FUNCTION ITS_A_NEW_MONTH

!------------------------------------------------------------------------------

      FUNCTION ITS_MIDMONTH() RESULT( IS_MIDMONTH )
!
!******************************************************************************
!  Function ITS_MIDMONTH returns TRUE if it's the middle of a month
!  -sas 10/10/05
!
!  NOTES:
!******************************************************************************
!
      ! Function value
      LOGICAL :: IS_MIDMONTH
      
      !=================================================================
      ! ITS_MIDMONTH begins here!
      !=================================================================

      ! Test for the 16th of the month at 0 GMT
      IS_MIDMONTH = ( DAY == 16 .and. NHMS == 000000 )

      ! Return to calling program
      END FUNCTION ITS_MIDMONTH

!------------------------------------------------------------------------------

      FUNCTION ITS_A_NEW_DAY( ) RESULT( IS_NEW_DAY )
!
!******************************************************************************
!  Function ITS_A_NEW_DAY returns TRUE if it's the first timestep of a new
!  day (it also returns TRUE on the first timestep of the run).  This is 
!  useful for setting flags for reading in data. (bmy, 4/1/04)
!
!  NOTES:
!******************************************************************************
!
      ! Arguments
      LOGICAL :: IS_NEW_DAY
      
      !=================================================================
      ! ITS_A_NEW_DAY begins here!
      !=================================================================
      IF ( NHMS == 000000 ) THEN

         ! Test if it's 0 GMT
         IS_NEW_DAY = .TRUE.

      ELSE IF ( NYMD == NYMDb .and. NHMS == NHMSb ) THEN

         ! Also return TRUE if it's the start of the run
         ! (since files will need to be read in from disk)
         IS_NEW_DAY = .TRUE.

      ELSE

         ! Otherwise, it's not a new year
         IS_NEW_DAY = .FALSE.
         
      ENDIF

      ! Return to calling program
      END FUNCTION ITS_A_NEW_DAY

!------------------------------------------------------------------------------

      FUNCTION ITS_A_NEW_SEASON( ) RESULT( IS_NEW_SEASON )
!
!******************************************************************************
!  Function ITS_A_NEW_SEASON returns TRUE if it's a new season or FALSE
!  if it's not a new season.  Seasons are (1=DJF, 2=MAM, 3=JJA, 4=SON).
!  (bmy, 7/20/04)
!
!  NOTES:
!******************************************************************************
!
      ! Function value
      LOGICAL       :: IS_NEW_SEASON
      
      ! Local variables
      INTEGER, SAVE :: LAST_SEASON = -1
      
      !=================================================================
      ! ITS_A_NEW_SEASON begins here!
      !=================================================================
      IF ( NSEASON /= LAST_SEASON ) THEN
         IS_NEW_SEASON = .TRUE.
         LAST_SEASON   = NSEASON
      ELSE
         IS_NEW_SEASON = .FALSE.
      ENDIF

      ! Return to calling program
      END FUNCTION ITS_A_NEW_SEASON

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

      FUNCTION TIMESTAMP_STRING( YYYYMMDD, HHMMSS ) RESULT( TIME_STR )
!
!******************************************************************************
!  TIMESTAMP_STRING returns a formatted string "YYYY/MM/DD HH:MM" for the a
!  date and time specified by YYYYMMDD and HHMMSS.  If YYYYMMDD and HHMMSS are
!  omitted, then TIMESTAMP_STRING will create a formatted string for the 
!  current date and time. (bmy, 3/21/03, 12/2/03)
!                                                                          
!  NOTES:
!  (1 ) Now use ENCODE statement for PGI/F90 on Linux (bmy, 9/29/03)
!  (2 ) Now add optional arguments YYYYMMDD and HHMMSS (bmy, 10/27/03)
!  (3 ) Renamed LINUX to LINUX_PGI (bmy, 12/2/03)
!******************************************************************************
!     
#     include "define.h"

      ! Arguments
      INTEGER, INTENT(IN), OPTIONAL :: YYYYMMDD, HHMMSS 

      ! Local variables
      INTEGER                       :: THISYEAR, THISMONTH,  THISDAY
      INTEGER                       :: THISHOUR, THISMINUTE, THISSECOND

      ! Function value
      CHARACTER(LEN=16)             :: TIME_STR
      
      !=================================================================
      ! TIMESTAMP_STRING begins here!
      !=================================================================

      ! If YYYYMMDD is passed, then use that date.  Otherwise use the 
      ! current date stored in global variables YEAR, MONTH, DAY.
      IF ( PRESENT( YYYYMMDD ) ) THEN
         CALL YMD_EXTRACT( YYYYMMDD, THISYEAR, THISMONTH, THISDAY )
      ELSE
         THISYEAR  = YEAR
         THISMONTH = MONTH
         THISDAY   = DAY
      ENDIF
         
      ! If HHMMSS is passed, then use that time.  Otherwise use the 
      ! current time stored in global variables HOUR and MINUTE.
      IF ( PRESENT( HHMMSS ) ) THEN
         CALL YMD_EXTRACT( HHMMSS, THISHOUR, THISMINUTE, THISSECOND )
      ELSE
         THISHOUR   = HOUR
         THISMINUTE = MINUTE
      ENDIF

#if   defined( LINUX_PGI ) 
      
      ! For PGI/F90 Linux, we must use the ENCODE command
      ! to convert from numeric to string format (bmy, 9/29/03)
      ENCODE( 16, 100, TIME_STR ) THISYEAR, THISMONTH, 
     &                            THISDAY,  THISHOUR, THISMINUTE

#else

      ! For other platforms, we can just use a FORTRAN internal write
      WRITE( TIME_STR, 100 ) THISYEAR, THISMONTH, 
     &                       THISDAY,  THISHOUR, THISMINUTE

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
!  a filename string with the actual values. (bmy, 6/27/02, 7/20/04)
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
!  (2 ) Now replace 2 and 4 digit year strings for all models (bmy, 10/23/03)
!  (3 ) Renamed LINUX to LINUX_PGI (bmy, 12/2/03)
!  (4 ) Now do not replace "ss" with seconds, as the smallest GEOS-CHEM
!        timestep is in minutes. (bmy, 7/20/04)
!******************************************************************************
!      
      ! References to F90 modules
      USE CHARPAK_MOD, ONLY : STRREPL

#     include "define.h"

      ! Arguments
      CHARACTER(LEN=*), INTENT(INOUT) :: FILENAME
      INTEGER,          INTENT(IN)    :: YYYYMMDD, HHMMSS

      ! Local variables
      INTEGER                         :: YYYY, YY, MM, DD, HH, II, SS
      CHARACTER(LEN=2)                :: MM_STR, DD_STR
      CHARACTER(LEN=2)                :: HH_STR, II_STR, SS_STR
      CHARACTER(LEN=2)                :: YY_STR
      CHARACTER(LEN=4)                :: YYYY_STR

      !=================================================================
      ! EXPAND_DATE begins here!
      !=================================================================

      ! Extract today's date into year, month, and day sections
      CALL YMD_EXTRACT( YYYYMMDD, YYYY, MM, DD )

      ! Extract today's time into HH, MM, and SS sections
      ! (rename minutes to II so as not to overwrite MM)
      CALL YMD_EXTRACT( HHMMSS, HH, II, SS )

      ! 2-digit year number (e.g. "97" instead of "1997")
      YY = YYYY - 1900
      IF ( YY >= 100 ) YY = YY - 100

#if   defined( LINUX_PGI )
      
      ! Use ENCODE statement for PGI/Linux (bmy, 9/29/03)
      ENCODE( 4, '(i4.4)', YYYY_STR ) YYYY
      ENCODE( 2, '(i2.2)', YY_STR   ) YY
      ENCODE( 2, '(i2.2)', MM_STR   ) MM
      ENCODE( 2, '(i2.2)', DD_STR   ) DD
      ENCODE( 2, '(i2.2)', HH_STR   ) HH
      ENCODE( 2, '(i2.2)', II_STR   ) II

#else

      ! For other platforms, use an F90 internal write (bmy, 9/29/03)
      WRITE( YYYY_STR, '(i4.4)' ) YYYY
      WRITE( YY_STR,   '(i2.2)' ) YY
      WRITE( MM_STR,   '(i2.2)' ) MM
      WRITE( DD_STR,   '(i2.2)' ) DD
      WRITE( HH_STR,   '(i2.2)' ) HH
      WRITE( II_STR,   '(i2.2)' ) II

#endif

      ! Replace YYYY, MM, DD, HH tokens w/ actual values 
      CALL STRREPL( FILENAME, 'YYYY', YYYY_STR )
      CALL STRREPL( FILENAME, 'YY',   YY_STR   )
      CALL STRREPL( FILENAME, 'MM',   MM_STR   )
      CALL STRREPL( FILENAME, 'DD',   DD_STR   )
      CALL STRREPL( FILENAME, 'hh',   HH_STR   )
      CALL STRREPL( FILENAME, 'mm',   II_STR   )

      ! Return to calling program
      END SUBROUTINE EXPAND_DATE

!------------------------------------------------------------------------------

      SUBROUTINE SYSTEM_DATE_TIME( SYS_NYMD, SYS_NHMS )
!
!******************************************************************************
!  Subroutine SYSTEM_DATE_TIME returns the actual local date and time 
!  (as opposed to the model date and time).  (bmy, 5/2/05)
!
!  Arguments as Output:
!  ============================================================================
!  (1 ) SYS_NYMD (INTEGER) : System date (local time) in YYYYMMDD format
!  (2 ) SYS_NHMS (INTEGER) : System time (local time) in HHMMSS   format
!
!  NOTES: 
!******************************************************************************
!
      ! Arguments
      INTEGER, INTENT(OUT) :: SYS_NYMD
      INTEGER, INTENT(OUT) :: SYS_NHMS

      ! Local variables
      INTEGER              :: V(8)
      CHARACTER(LEN=8)     :: D
      CHARACTER(LEN=10)    :: T

      !=================================================================
      ! SYSTEM_DATE_TIME begins here!
      !=================================================================

      ! Initialize
      D = 'ccyymmdd'
      T = 'hhmmss.sss'

      ! Call the F90 intrinsic routine DATE_AND_TIME
      ! Return values are (/YYYY, MM, DD, GMT_MIN, HH, MM, SS, MSEC/)
      CALL DATE_AND_TIME( DATE=D, TIME=T, VALUES=V )

      ! Save to YYYYMMDD and HHMMSS format
      SYS_NYMD = ( V(1) * 10000 ) + ( V(2) * 100 ) + V(3) 
      SYS_NHMS = ( V(5) * 10000 ) + ( V(6) * 100 ) + V(7)

      ! Return to calling program
      END SUBROUTINE SYSTEM_DATE_TIME

!------------------------------------------------------------------------------

      FUNCTION SYSTEM_TIMESTAMP() RESULT( STAMP )
!
!******************************************************************************
!  Function SYSTEM_TIMESTAMP returns a 16 character string with the system
!  date and time in YYYY/MM/DD HH:MM format. (bmy, 5/3/05)
!
!  NOTES:
!******************************************************************************
!
      ! Local variables
      INTEGER           :: SYS_NYMD, SYS_NHMS
      CHARACTER(LEN=16) :: STAMP

      !=================================================================
      ! SYSTEM_TIMESTAMP begins here!
      !=================================================================

      ! Get system date and time
      CALL SYSTEM_DATE_TIME( SYS_NYMD, SYS_NHMS )

      ! Create a string w/ system date & time
      STAMP = TIMESTAMP_STRING( SYS_NYMD, SYS_NHMS )

      ! Return to calling program 
      END FUNCTION SYSTEM_TIMESTAMP

!------------------------------------------------------------------------------

      END MODULE TIME_MOD




