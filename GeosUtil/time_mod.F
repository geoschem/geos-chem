!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: time_mod.F
!
! !DESCRIPTION: Module TIME\_MOD contains GEOS-Chem date and time variables
!  and timesteps, and routines for accessing them.
!\\
!\\
! !INTERFACE: 
!
      MODULE TIME_MOD
! 
! !USES:
!
      USE PRECISION_MOD    ! For GEOS-Chem Precision (fp)

      IMPLICIT NONE
!
! !PRIVATE TYPES:
!
      PRIVATE
      PRIVATE :: HG2_DIAG !hma, 25 Oct 2011
!
! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC  :: SET_CURRENT_TIME
      PUBLIC  :: SET_BEGIN_TIME
      PUBLIC  :: SET_END_TIME
      PUBLIC  :: SET_DIAGb
      PUBLIC  :: SET_DIAGe
      PUBLIC  :: SET_TIMESTEPS
      PUBLIC  :: SET_CT_CHEM
      PUBLIC  :: SET_CT_CONV
      PUBLIC  :: SET_CT_DYN
      PUBLIC  :: SET_CT_EMIS
      PUBLIC  :: SET_CT_DIAG
      PUBLIC  :: SET_CT_RAD
      PUBLIC  :: SET_CT_A1 
      PUBLIC  :: SET_CT_A3
      PUBLIC  :: SET_CT_I3
      PUBLIC  :: SET_CT_XTRA
      PUBLIC  :: SET_ELAPSED_SEC
      PUBLIC  :: GET_JD
      PUBLIC  :: GET_ELAPSED_SEC
      PUBLIC  :: GET_NYMDb
      PUBLIC  :: GET_NHMSb
      PUBLIC  :: GET_NYMDe
      PUBLIC  :: GET_NHMSe
      PUBLIC  :: GET_NYMD
      PUBLIC  :: GET_NHMS
      PUBLIC  :: GET_TIME_AHEAD
      PUBLIC  :: GET_MONTH
      PUBLIC  :: GET_DAY
      PUBLIC  :: GET_YEAR
      PUBLIC  :: GET_HOUR
      PUBLIC  :: GET_MINUTE
      PUBLIC  :: GET_SECOND
      PUBLIC  :: GET_DAY_OF_YEAR
      PUBLIC  :: GET_DAY_OF_WEEK
      PUBLIC  :: GET_DAY_OF_WEEK_LT
      PUBLIC  :: GET_GMT
      PUBLIC  :: GET_TAU
      PUBLIC  :: GET_TAUb
      PUBLIC  :: GET_TAUe
      PUBLIC  :: GET_DIAGb
      PUBLIC  :: GET_DIAGe
      PUBLIC  :: GET_LOCALTIME
      PUBLIC  :: GET_SEASON
      PUBLIC  :: GET_TS_CHEM
      PUBLIC  :: GET_TS_CONV
      PUBLIC  :: GET_TS_DIAG
      PUBLIC  :: GET_TS_DYN
      PUBLIC  :: GET_TS_EMIS
      PUBLIC  :: GET_TS_UNIT
      PUBLIC  :: GET_TS_RAD
      PUBLIC  :: GET_CT_CHEM
      PUBLIC  :: GET_CT_CONV
      PUBLIC  :: GET_CT_DYN
      PUBLIC  :: GET_CT_EMIS
      PUBLIC  :: GET_CT_RAD
      PUBLIC  :: GET_CT_A1
      PUBLIC  :: GET_CT_A3
      PUBLIC  :: GET_CT_I3
      PUBLIC  :: GET_CT_XTRA
      PUBLIC  :: GET_CT_DIAG
      PUBLIC  :: GET_A1_TIME
      PUBLIC  :: GET_A3_TIME
      PUBLIC  :: GET_I3_TIME
      PUBLIC  :: GET_FIRST_A1_TIME
      PUBLIC  :: GET_FIRST_A3_TIME
      PUBLIC  :: GET_FIRST_I3_TIME
      PUBLIC  :: ITS_TIME_FOR_CHEM
      PUBLIC  :: ITS_TIME_FOR_CONV
      PUBLIC  :: ITS_TIME_FOR_DYN
      PUBLIC  :: ITS_TIME_FOR_EMIS
      PUBLIC  :: ITS_TIME_FOR_EXCH
      PUBLIC  :: ITS_TIME_FOR_UNIT
      PUBLIC  :: ITS_TIME_FOR_DIAG
      PUBLIC  :: ITS_TIME_FOR_RT
      PUBLIC  :: ITS_TIME_FOR_SURFACE_RAD
      PUBLIC  :: ITS_TIME_FOR_A1
      PUBLIC  :: ITS_TIME_FOR_A3
      PUBLIC  :: ITS_TIME_FOR_I3
      PUBLIC  :: ITS_TIME_FOR_EXIT
      PUBLIC  :: ITS_TIME_FOR_BPCH
      PUBLIC  :: ITS_A_LEAPYEAR
      PUBLIC  :: ITS_A_NEW_YEAR
      PUBLIC  :: ITS_A_NEW_MONTH
      PUBLIC  :: ITS_MIDMONTH
      PUBLIC  :: ITS_A_NEW_DAY
      ! Add for FLAMBE compatibility (skim, 6/21/13)
      PUBLIC  :: ITS_A_NEW_HOUR
      PUBLIC  :: ITS_A_NEW_SEASON
      PUBLIC  :: PRINT_CURRENT_TIME
      PUBLIC  :: TIMESTAMP_STRING
      PUBLIC  :: YMD_EXTRACT
      PUBLIC  :: EXPAND_DATE
      PUBLIC  :: Valid_Date
      PUBLIC  :: Valid_Time
      PUBLIC  :: SYSTEM_DATE_TIME
      PUBLIC  :: SYSTEM_TIMESTAMP
      PUBLIC  :: TIMESTAMP_DIAG
      PUBLIC  :: GET_NYMD_DIAG
      PUBLIC  :: GET_Hg2_DIAG !hma, 25 Oct 2011
      PUBLIC  :: SET_Hg2_DIAG !hma, 25 Oct 2011
![eml
      PUBLIC  :: SET_HISTYR
      PUBLIC  :: GET_HISTYR
!eml]
#if defined( ESMF_ ) || defined( MODEL_ )
      PUBLIC  :: Accept_External_Date_Time
#endif
!
! !REMARKS:
!  References:
!  ---------------------------------------------------------------------------
!  (1) "Practical Astronomy with Your Calculator", 3rd Ed.
!        Peter Duffett-Smith, Cambridge UP, 1992, p9.
!  (2) Rounding algorithm from: Hultquist, P.F, "Numerical Methods for
!        Engineers and Computer Scientists", Benjamin/Cummings, Menlo Park CA,
!        1988, p. 20.
!  (3) Truncation algorithm from: http://en.wikipedia.org/wiki/Truncation
!
! !REVISION HISTORY:
!  21 Jun 2000 - R. Yantosca - Initial version
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
!  (25) Remove support for GEOS-1 and GEOS-STRAT met fields (bmy, 8/4/06)
!  (26) Further bug fix to skip over Feb 29th in GCAP (phs, bmy, 10/3/06)
!  (27) Moved ITS_TIME_FOR_BPCH here from "main.f" (bmy, 2/2/07)
!  (28) Add TS_DIAG and CT_DIAG variables to correctly output diagnostics 
!        (good time step).
!       Add SET_CT_DIAG and GET_CT_DIAG to implement TS_DIAG correctly.
!       (ccc, 5/21/09)
!  (29) Add NYMD_DIAG, GET_NYMD_DIAG, TIMESTAMP_DIAG to get the good timestamp
!       for diagnostic filenames (ccc, 8/12/09)
!  15 Jan 2010 - R. Yantosca - Added ProTeX headers
!  27 Apr 2010 - R. Yantosca - Added OFFSET argument to GET_LOCALTIME
!  27 Apr 2010 - R. Yantosca - Added TS_SUN_2 to hold 1/2 of the interval
!                              for computing SUNCOS.
!  27 Apr 2010 - R. Yantosca - Added public routine GET_TS_SUN_2
!  19 Aug 2010 - R. Yantosca - Added variable CT_A1 and routine SET_CT_A1
!  20 Aug 2010 - R. Yantosca - Added function ITS_TIME_FOR_A1
!  27 Sep 2010 - R. Yantosca - Added function GET_FIRST_I6_TIME
!  17 Dec 2010 - R. Yantosca - Bug fix for HHMMSS=240000 in GET_TIME_AHEAD
!  27 Mar 2011 - R. Yantosca - Bug fix for GCAP leap year problem
!  29 Jul 2011 - R. Yantosca - Add LEAP_YEAR_DAYS as a SAVEd module variable
!  17 Feb 2011 - R. Yantosca - Added ITS_TIME_FOR_A6UPDATE for APM (G. Luo)
!  07 Oct 2011 - M. Payer    - Modifications for central chemistry timestep
!  07 Oct 2011 - R. Yantosca - Remove obsolete TS_SUN_2, GET_TS_SUN_2
!  07 Oct 2011 - R. Yantosca - Remove obsolete OFFSET argument to GET_LOCALTIME
!  12 Oct 2011 - R. Yantosca - Modified ITS_A_NEW_MONTH for central chem step
!  25 Oct 2011 - H. Amos     - bring Hg2 gas-particle partitioning code into
!                              v9-01-02
!  02 Feb 2012 - R. Yantosca - Added modifications for GEOS-5.7.x met fields
!  03 Feb 2012 - R. Yantosca - Added I3 fields timestep variable & routines
!  28 Feb 2012 - R. Yantosca - Removed support for GEOS-3
!  01 Mar 2012 - R. Yantosca - GET_LOCALTIME now takes (I,J,L,GMT) as args
!  14 Jun 2013 - R. Yantosca - Now compute day of week in SET_CURRENT_TIME
!  14 Jun 2013 - R. Yantosca - Added comments to module variable declarations
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!  18 Sep 2013 - M. Long     - Now use #if defined( ESMF_ ) for HPC code
!  02 Dec 2014 - M. Yannetti - Added PRECISION_MOD
!  17 Dec 2014 - R. Yantosca - Need to keep time variables as 8-byte precision
!  06 Jan 2015 - R. Yantosca - Add ITS_TIME_FOR_EXCH function for 2-way nests
!  23 Jun 2016 - R. Yantosca - Remove references to APM code; it is no longer
!                              compatible with the FlexChem implementation
!  23 Jun 2016 - R. Yantosca - Remove references to APM code; it is no longer
!                              compatible with the FlexChem implementation
!  29 Nov 2016 - R. Yantosca - grid_mod.F90 is now gc_grid_mod.F90
!  24 Aug 2017 - M. Sulprizio- Remove support for GCAP, GEOS-4, GEOS-5 and
!                              MERRA; Remove obsolete A6 and I6 routines
!  06 Feb 2018 - E. Lundgren - Change timestep units from min to secf
!  06 Jul 2018 - R. Yantosca - Add Valid_Date and Valid_Time functions
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !PRIVATE TYPES:
!
      !---------------------------------
      ! GMT dates @ start & end of run
      !---------------------------------
      INTEGER             :: NYMDb        ! YYYY/MM/DD on which to start run
      INTEGER             :: NHMSb        ! hh:mm:ss   on which to start run
      REAL(f8)            :: TAUb         ! Hours since 1/1/1985 @ start of run

      INTEGER             :: NYMDe        ! YYYY/MM/DD on which to end run
      INTEGER             :: NHMSe        ! hh:mm:ss   on which to end run
      REAL(f8)            :: TAUe         ! Hours since 1/1/1985 @ end of run

      !---------------------------------
      ! Current GMT date/time values
      !---------------------------------
      INTEGER             :: NYMD         ! YYYY/MM/DD date
      INTEGER             :: NHMS         ! hh:mm:ss time
      INTEGER             :: MONTH        ! Month (1-12)
      INTEGER             :: DAY          ! Day of month (1-31)
      INTEGER             :: YEAR         ! Year (YYYY format)
      INTEGER             :: HOUR         ! hour (0-23)
      INTEGER             :: MINUTE       ! Minutes (0-59)
      INTEGER             :: SECOND       ! Seconds (0-59)
      INTEGER             :: NSEASON      ! Season (1=DJF,2=MAM,3=JJA,4=SON)
      INTEGER             :: DAY_OF_YEAR  ! Day of year (0-365 or 0-366 if LY)
      INTEGER             :: DAY_OF_WEEK  ! Day of week (0=Sun,1=Mon,.., 6=Sat)
      REAL(f8)            :: TAU          ! Hours since 1/1/1985
      REAL(f8)            :: GMT          ! Current GMT time [decimal hours]

      !---------------------------------
      ! Elapsed time values
      !---------------------------------
      INTEGER             :: ELAPSED_SEC    ! # of seconds into the run
      REAL(f8)            :: LEAP_YEAR_DAYS ! # of leap yrs since start of run

      !---------------------------------
      ! Diagnostic date/time values
      !---------------------------------
      INTEGER             :: NYMD_DIAG    ! Diagnostic date for output
      REAL(f8)            :: DIAGb        ! Start of diag averaging interval
      REAL(f8)            :: DIAGe        ! End   of diag averaging interval

      !---------------------------------
      ! Timestep values
      !---------------------------------
      INTEGER             :: TS_CHEM      ! Chemistry timestep       [sec]
      INTEGER             :: TS_CONV      ! Convection timestep      [sec]
      INTEGER             :: TS_DIAG      ! Diagnostic timestep      [sec]
      INTEGER             :: TS_DYN       ! Dynamic timestep         [sec]
      INTEGER             :: TS_EMIS      ! Emissions timestep       [sec]
      INTEGER             :: TS_UNIT      ! Unit conversion timestep [sec]
      INTEGER             :: TS_RAD       ! Radiation timestep       [sec]

      !---------------------------------
      ! Counter values
      !---------------------------------
      INTEGER             :: CT_CHEM      ! # of elapsed chemistry timesteps
      INTEGER             :: CT_CONV      ! # of elapsed convection timesteps
      INTEGER             :: CT_DYN       ! # of elapsed dynamic timesteps
      INTEGER             :: CT_EMIS      ! # of elapsed emission timesteps
      INTEGER             :: CT_RAD       ! # of elapsed radiation timesteps
      INTEGER             :: CT_A1        ! # of elapsed A1 met field reads
      INTEGER             :: CT_A3        ! # of elapsed A3 met field reads
      INTEGER             :: CT_I3        ! # of elapsed I3 met field reads
      INTEGER             :: CT_XTRA      ! # of elapsed XTRA met field reads
      INTEGER             :: CT_DIAG      ! # of elapsed diagnostic timesteps
      INTEGER             :: Hg2_DIAG     ! # of elapsed Hg2 diag outputs

      !---------------------------------
      ! For historical emissions
      !---------------------------------
      INTEGER             :: HISTYRA      ! Year for historical emissions
!
! !DEFINED PARAMETERS:
!
      ! Astronomical Julian Date at 0 GMT, 1 Jan 1985
      REAL(f8), PARAMETER :: JD85 = 2446066.5e+0_f8

      CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Set_current_time
!
! !DESCRIPTION: Subroutine SET\_CURRENT\_TIME takes in the elapsed time in 
!  minutes since the start of a GEOS-Chem simulation and sets the GEOS-Chem 
!  time variables accordingly.  NOTE: All time variables are returned
!  w/r/t Greenwich Mean Time (aka Universal Time).
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE SET_CURRENT_TIME
!
! !USES:
!
      USE JULDAY_MOD, ONLY : JULDAY, CALDATE
!
! !REMARKS:
!  The GEOS met fields are assimilated data, and therefore contain data on
!  the leap-year days.  SET_CURRENT_TIME computes the days according to the
!  Astronomical Julian Date algorithms (in "julday_mod.f"), which contain
!  leap-year days.
! 
! !REVISION HISTORY: 
!  05 Feb 2006 - R. Yantosca - Initial Version
!  (1 ) GCAP/GISS fields don't have leap years, so if JULDAY says it's 
!        Feb 29th, reset MONTH, DAY, JD1 to Mar 1st. (swu, bmy, 8/29/05)
!  (2 ) Now references "define.h".  Now add special handling to skip from
!        Feb 28th to Mar 1st for GCAP model. (swu, bmy, 4/24/06)
!  (3 ) Fix bug in case of GCAP fields for runs that start during leap year
!       and after February 29 (phs, 9/27/06)  
!  15 Jan 2010 - R. Yantosca - Added ProTeX headers
!  29 Jul 2011 - R. Yantosca - Bug fix: For GCAP, we need to skip over the
!                              # of leap-year-days that have already occurred
!                              when going from Julian date to Y/M/D date
!  14 Jun 2013 - R. Yantosca - Now move the day of week computation here
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!  06 Feb 2018 - E. Lundgren - Update Julian Data calc to use new minimum time
!                              unit of seconds (previously minutes)
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL, SAVE :: FIRST = .TRUE.
      REAL(f4)      :: TMP
      REAL(f8)      :: A, B, JD, JD0, JD1, THISDAY

      !=================================================================
      ! SET_CURRENT_TIME begins here!
      !=================================================================

      ! Initialize LEAP_YEAR_DAYS
      IF ( FIRST ) THEN 
         LEAP_YEAR_DAYS = 0e+0_f8
         FIRST          = .FALSE.
      ENDIF

      ! JD0: Astronomical Julian Date at start of GEOS-Chem run
      JD0 = GET_JD( NYMDb, NHMSb )

      ! JD1: Astronomical Julian Date at current time
      JD1 = JD0 + ( DBLE( ELAPSED_SEC ) / 86400e+0_f8 )

      ! Call CALDATE to compute the current YYYYMMDD and HHMMSS
      CALL CALDATE( DBLE(JD1), NYMD, NHMS )

      ! Extract current year, month, day from NYMD
      CALL YMD_EXTRACT( NYMD, YEAR, MONTH, DAY )

      ! Extract current hour, minute, second from NHMS
      CALL YMD_EXTRACT( NHMS, HOUR, MINUTE, SECOND )

      ! Fix minutes & seconds for display purposes (esp. for 1x1)
      IF ( SECOND              == 59 ) SECOND = 0
      IF ( MOD( MINUTE+1, 10 ) == 0  ) MINUTE = MINUTE + 1

      !=================================================================
      ! Compute other GEOS-Chem timing variables
      !=================================================================

      ! Current Greenwich Mean Time
      GMT         = ( DBLE( HOUR )                ) + 
     &              ( DBLE( MINUTE ) /   60e+0_f8 ) + 
     &              ( DBLE( SECOND ) / 3600e+0_f8 )

      ! Days elapsed in this year (0-366)
      DAY_OF_YEAR = JD1 - JULDAY( YEAR, 1, 0d0 ) ! has to be REAL*8

      ! TAU value (# of hours since 1 Jan 1985)
      ! NOTE: TMP is REAL*4 to prevent precision problems
      TMP         = ( JD1 - JD85 ) * 24e+0_f4
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

      !=================================================================
      ! Compute day of week w/r/t the GMT date
      ! (moved here from routine GET_DAY_OF_WEEK)
      !=================================================================

      ! Get fractional GMT day
      THISDAY     = DAY + ( GMT / 24e+0_f8 ) ! real*8

      ! Get current Julian date 
      JD          = JULDAY( YEAR, MONTH, DBLE(THISDAY) ) ! real*8

      ! Add 1.5 to JD and divide by 7
      A           = ( JD + 1.5e+0_f8 ) / 7e+0_f8

      ! Take fractional part and multiply by 7
      B           = ( A - INT( A ) ) * 7e+0_f8

      ! NOTE: This is not in the Duffett-Smith book, but when using
      ! REAL(f8), we have to round B to 4 decimal places in order
      ! to avoid floating-point errors. (bmy, 6/13/13)
      B           = INT( NINT( B*1e+5_f8 + SIGN( 5e+0_f8, B ) ) 
     &            / 10e+0_f8 ) / 1e+4_f8

      ! Discard the fractional part of B
      DAY_OF_WEEK = INT( B )

      END SUBROUTINE SET_CURRENT_TIME
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Set_begin_time
!
! !DESCRIPTION: Subroutine SET\_BEGIN\_TIME initializes NYMDb, NHMSb, and TAUb,
!  which are the YYYYMMDD, HHMMSS, and hours since 1/1/1985 corresponding to 
!  the beginning date and time of a GEOS-Chem run.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE SET_BEGIN_TIME( THISNYMDb, THISNHMSb, RC )
!
! !USES:
!
      USE ErrCode_Mod
!
! !INPUT PARAMETERS: 
!
      INTEGER, INTENT(IN)  :: THISNYMDb   ! YYYYMMDD @ start of G-C simulation
      INTEGER, INTENT(IN)  :: THISNHMSb   ! HHMMSS   @ start of G-C simulation
!
! !OUTPUT PARAMETERS:
!
      INTEGER, INTENT(OUT) :: RC          ! Success or failure
! 
! !REVISION HISTORY: 
!  20 Jul 2004 - R. Yantosca - Initial Version
!  15 Jan 2010 - R. Yantosca - Added ProTeX headers
!  16 Dec 2010 - R. Yantosca - Updated error check for THISNYMDe, since 
!                              MERRA met data goes back prior to 1985
!  08 Nov 2017 - R. Yantosca - Return error condition to calling program
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL(f4) :: TMP

      ! Strings
      CHARACTER(LEN=255) :: ErrMsg, ThisLoc

      !=================================================================
      ! SET_BEGIN_TIME begins here!
      !=================================================================

      ! Initialize
      RC      = GC_SUCCESS
      ErrMsg  = ''
      ThisLoc = 
     &   ' -> at Set_Begin_Time (in module GeosUtil/time_mod.F)'

      ! Error check THISNHMSb
      IF ( THISNHMSb > 235959 ) THEN
         ErrMsg = 'NHMSb cannot be greater than 23:59:59!'
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF

      ! Error check THISNYMDb
      IF ( THISNYMDb < 19000101 ) THEN
         ErrMsg = 'NYMDb must be in the format YYYYMMDD!'
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF

      ! Initialize NYMDb, NHMSb
      NYMDb = THISNYMDb
      NHMSb = THISNHMSb

      ! TAUb value (TMP is REAL*4 to prevent precision problems)
      TMP   = ( GET_JD( NYMDb, NHMSb ) - JD85 ) * 24e+0_f4
      TAUb  = DBLE( TMP )

      ! Also initialize ELAPSED_SEC
      ELAPSED_SEC = 0

      ! Return to calling program
      END SUBROUTINE SET_BEGIN_TIME
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Set_end_time
!
! !DESCRIPTION: Subroutine SET\_END\_TIME initializes NYMDe, NHMSe, and TAUe, 
!  which are the YYYYMMDD, HHMMSS, and hours since 1/1/1985 corresponding to 
!  the ending date and time of a GEOS-Chem run.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE SET_END_TIME( THISNYMDe, THISNHMSe, RC )
!
! !USES:
!
      USE ErrCode_Mod
!
! !INPUT PARAMETERS: 
!
      INTEGER, INTENT(IN)  :: THISNYMDe   ! YYYYMMDD @ end of G-C simulation
      INTEGER, INTENT(IN)  :: THISNHMSe   ! HHMMSS   @ end of G-C simulation
!
! !OUTPUT PARAMETERS:
!
      INTEGER, INTENT(OUT) :: RC          ! Success or failure
! 
! !REVISION HISTORY: 
!  20 Jul 2004 - R. Yantosca - Initial Version
!  15 Jan 2010 - R. Yantosca - Added ProTeX headers
!  16 Dec 2010 - R. Yantosca - Updated error check for THISNYMDe, since 
!                              MERRA met data goes back prior to 1985
!  08 Nov 2017 - R. Yantosca - Return error condition to calling program
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! Scalars
      REAL(f4) :: TMP

      ! Strings
      CHARACTER(LEN=255) :: ErrMsg, ThisLoc

      !=================================================================
      ! SET_END_TIME begins here!
      !=================================================================

      ! Initialize
      RC      = GC_SUCCESS
      ErrMsg  = ''
      ThisLoc = 
     &   ' -> at Set_End_Time (in module GeosUtil/time_mod.F)'

      ! Error check THISNHMS
      IF ( THISNHMSe > 235959 ) THEN
         ErrMsg = 'NHMSe cannot be greater than 23:59:59!'
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF

      ! Error check THISNYMDe
      IF ( THISNYMDe < 19000101 ) THEN
         ErrMsg = 'NYMDe must be in the format YYYYMMDD!'
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF

      ! Initialize NYMDe, NHMSe
      NYMDe = THISNYMDe
      NHMSe = THISNHMSe

      ! TAUe value (TMP is REAL*4 to prevent precision problems)
      TMP   = ( GET_JD( NYMDe, NHMSe ) - JD85 ) * 24e+0_f4
      TAUe  = DBLE( TMP )

      ! Return to calling program
      END SUBROUTINE SET_END_TIME
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Set_diagb
!
! !DESCRIPTION: Subroutine SET\_DIAGb initializes DIAGb, the TAU value at the 
!  start of the diagnostic averaging interval.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE SET_DIAGb( THISDIAGb )
!
! !INPUT PARAMETERS: 
!
      REAL(f8), INTENT(IN) :: THISDIAGb  ! Initial DIAGb value [hrs from 1/1/85]
! 
! !REVISION HISTORY: 
!  21 Mar 2003 - R. Yantosca - Initial Version
!  15 Jan 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
      DIAGb = THISDIAGb

      END SUBROUTINE SET_DIAGb
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Set_diage
!
! !DESCRIPTION: Subroutine SET\_DIAGe initializes DIAGe, the TAU value at the 
!  end of the diagnostic averaging interval.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE SET_DIAGe( THISDIAGe )
!
! !INPUT PARAMETERS: 
!
      REAL(f8), INTENT(IN) :: THISDIAGe  ! Initial DIAGe value [hrs from 1/1/85]
! 
! !REVISION HISTORY: 
!  21 Mar 2003 - R. Yantosca - Initial Version
!  15 Jan 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
      DIAGe = THISDIAGe

      END SUBROUTINE SET_DIAGe
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Set_timesteps
!
! !DESCRIPTION: Subroutine SET\_TIMESTEPS initializes the timesteps for 
!  dynamics, convection, chemistry, emissions, and diagnostics.  
!  Counters are also zeroed. 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE SET_TIMESTEPS( am_I_Root,
     &                          CHEMISTRY, CONVECTION, DYNAMICS,  
     &                          EMISSION,  UNIT_CONV,  DIAGNOS,
     &                          RADIATION )
!
! !INPUT PARAMETERS:
!
      LOGICAL, INTENT(IN) :: am_I_Root    ! Is this the root CPU?
      INTEGER, INTENT(IN) :: CHEMISTRY    ! Chemistry  timestep [sec]
      INTEGER, INTENT(IN) :: CONVECTION   ! Convection timestep [sec]
      INTEGER, INTENT(IN) :: DYNAMICS     ! Dynamic    timestep [sec]
      INTEGER, INTENT(IN) :: EMISSION     ! Emission   timestep [sec]
      INTEGER, INTENT(IN) :: UNIT_CONV    ! Unit conve timestep [sec]
      INTEGER, INTENT(IN) :: DIAGNOS      ! Diagnostic timestep [sec]
      INTEGER, INTENT(IN) :: RADIATION    ! Radiation  timestep [sec]
! 
! !REVISION HISTORY: 
!  05 Feb 2003 - R. Yantosca - Initial Version
!  (1 ) Suppress some output lines (bmy, 7/20/04)
!  (2 ) Also zero CT_XTRA (tmf, bmy, 10/20/05)
!  (3 ) Add TS_DIAG as the diagnostic timestep. (ccc, 5/13/09)
!  15 Jan 2010 - R. Yantosca - Added ProTeX headers
!  27 Apr 2010 - R. Yantosca - Now add SUNCOS argument to set 1/2 of the
!                              interval for computing the cosine of the
!                              solar zenith angle.
!  07 Oct 2011 - R. Yantosca - Remove obsolete SUNCOS argument
!  30 Jul 2012 - R. Yantosca - Now accept am_I_Root as an argument when
!                              running with the traditional driver main.F
!  06 Feb 2018 - E. Lundgren - Update timestep unit from min to sec
!EOP
!------------------------------------------------------------------------------
!BOC
      ! Initialize timesteps
      TS_CHEM  = CHEMISTRY
      TS_CONV  = CONVECTION
      TS_DYN   = DYNAMICS
      TS_EMIS  = EMISSION
      TS_UNIT  = UNIT_CONV
      TS_DIAG  = DIAGNOS
      TS_RAD   = RADIATION

      ! Zero timestep counters
      CT_CHEM  = 0
      CT_CONV  = 0
      CT_DYN   = 0
      CT_EMIS  = 0
      CT_RAD   = 0
      CT_A3    = 0
      CT_XTRA  = 0
      CT_DIAG  = 0
      Hg2_DIAG = 0

      ! Echo to stdout
      IF ( am_I_Root ) THEN
       WRITE( 6, '(/,a)' ) 'SET_TIMESTEPS: setting GEOS-Chem timesteps!'
       WRITE( 6, '(  a)' ) '-------------------------------------------'
       WRITE( 6, '(''Chemistry  Timestep [sec] : '', i6)' ) TS_CHEM
       WRITE( 6, '(''Convection Timestep [sec] : '', i6)' ) TS_CONV
       WRITE( 6, '(''Dynamics   Timestep [sec] : '', i6)' ) TS_DYN
       WRITE( 6, '(''Emission   Timestep [sec] : '', i6)' ) TS_EMIS
       WRITE( 6, '(''Unit Conv  Timestep [sec] : '', i6)' ) TS_UNIT
       WRITE( 6, '(''Diagnostic Timestep [sec] : '', i6)' ) TS_DIAG
       WRITE( 6, '(''Radiation  Timestep [sec] : '', i6)' ) TS_RAD
      ENDIF

      END SUBROUTINE SET_TIMESTEPS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Set_ct_chem
!
! !DESCRIPTION: Subroutine SET\_CT\_CHEM increments CT\_CHEM, the counter
!  of chemistry timesteps executed thus far. 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE SET_CT_CHEM( INCREMENT, RESET )
!
! !INPUT PARAMETERS:
!
      LOGICAL, INTENT(IN), OPTIONAL :: INCREMENT  ! Increment counter?
      LOGICAL, INTENT(IN), OPTIONAL :: RESET      ! Reset counter?
! 
! !REVISION HISTORY: 
!  21 Mar 2009 - R. Yantosca - Initial Version
!  15 Jan 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
      IF ( PRESENT( INCREMENT ) ) THEN
         CT_CHEM = CT_CHEM + 1
      ELSE IF ( PRESENT( RESET ) ) THEN
         CT_CHEM = 0
      ENDIF

      END SUBROUTINE SET_CT_CHEM
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Set_ct_rad
!
! !DESCRIPTION: Subroutine SET\_CT\_RAD increments CT\_RAD, the
! counter of radiation timesteps executed thus far.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE SET_CT_RAD( INCREMENT, RESET )
!
! !INPUT PARAMETERS:
!
      LOGICAL, INTENT(IN), OPTIONAL :: INCREMENT  ! Increment counter?
      LOGICAL, INTENT(IN), OPTIONAL :: RESET      ! Reset counter?
!
! !REVISION HISTORY:
!  06 Oct 2012 - D. Ridley   - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
      IF ( PRESENT( INCREMENT ) ) THEN
         CT_RAD = CT_RAD + 1
      ELSE IF ( PRESENT( RESET ) ) THEN
         CT_RAD = 0
      ENDIF

      END SUBROUTINE SET_CT_RAD
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Set_hg2_diag
!
! !DESCRIPTION: Subroutine SET\_Hg2\_DIAG increments Hg2\_DIAG, the counter for
!  the number of times AAD03\_Fg and AD03\_Fp are recorded. (hma 20100218)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE SET_Hg2_DIAG( INCREMENT, RESET )
!
! !INPUT PARAMETERS:
!
      LOGICAL, INTENT(IN), OPTIONAL :: INCREMENT  ! Increment counter?
      LOGICAL, INTENT(IN), OPTIONAL :: RESET      ! Reset counter?
! 
! !REVISION HISTORY:
!  18 Feb 2012 - H. Amos     - Initial version 
!  07 Mar 2012 - M. Payer    - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
      IF ( PRESENT( INCREMENT ) ) THEN
         Hg2_DIAG = Hg2_DIAG + 1
      ELSE IF ( PRESENT( RESET ) ) THEN
         Hg2_DIAG = 0
      ENDIF

      END SUBROUTINE SET_Hg2_DIAG

!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Set_ct_conv
!
! !DESCRIPTION: Subroutine SET\_CT\_CONV increments CT\_CONV, the counter
!  of convection timesteps executed thus far. 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE SET_CT_CONV( INCREMENT, RESET )
!
! !INPUT PARAMETERS:
!
      LOGICAL, INTENT(IN), OPTIONAL :: INCREMENT  ! Increment counter?
      LOGICAL, INTENT(IN), OPTIONAL :: RESET      ! Reset counter?
! 
! !REVISION HISTORY: 
!  21 Mar 2009 - R. Yantosca - Initial Version
!  15 Jan 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
      IF ( PRESENT( INCREMENT ) ) THEN
         CT_CONV = CT_CONV + 1
      ELSE IF ( PRESENT( RESET ) ) THEN
         CT_CONV = 0 
      ENDIF

      END SUBROUTINE SET_CT_CONV
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Set_ct_dyn
!
! !DESCRIPTION: Subroutine SET\_CT\_DYN increments CT\_DYN, the counter
!  of dynamical timesteps executed thus far. 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE SET_CT_DYN( INCREMENT, RESET )
!
! !INPUT PARAMETERS:
!
      LOGICAL, INTENT(IN), OPTIONAL :: INCREMENT  ! Increment counter?
      LOGICAL, INTENT(IN), OPTIONAL :: RESET      ! Reset counter?
! 
! !REVISION HISTORY: 
!  21 Mar 2009 - R. Yantosca - Initial Version
!  15 Jan 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
      IF ( PRESENT( INCREMENT ) ) THEN
         CT_DYN = CT_DYN + 1
      ELSE IF ( PRESENT( RESET ) ) THEN
         CT_DYN = 0
      ENDIF

      END SUBROUTINE SET_CT_DYN
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Set_ct_emis
!
! !DESCRIPTION: Subroutine SET\_CT\_EMIS increments CT\_EMIS, the counter
!  of emission timesteps executed thus far. 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE SET_CT_EMIS( INCREMENT, RESET )
!
! !INPUT PARAMETERS:
!
      LOGICAL, INTENT(IN), OPTIONAL :: INCREMENT  ! Increment counter?
      LOGICAL, INTENT(IN), OPTIONAL :: RESET      ! Reset counter?
! 
! !REVISION HISTORY: 
!  21 Mar 2009 - R. Yantosca - Initial Version
!  15 Jan 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
      IF ( PRESENT( INCREMENT ) ) THEN
         CT_EMIS = CT_EMIS + 1
      ELSE IF ( PRESENT( RESET ) ) THEN
         CT_EMIS = 0
      ENDIF

      END SUBROUTINE SET_CT_EMIS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Set_ct_diag
!
! !DESCRIPTION: Subroutine SET\_CT\_DIAG increments CT\_DIAG, the counter
!  of largest timesteps executed thus far. 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE SET_CT_DIAG( INCREMENT, RESET )
!
! !INPUT PARAMETERS:
!
      LOGICAL, INTENT(IN), OPTIONAL :: INCREMENT  ! Increment counter?
      LOGICAL, INTENT(IN), OPTIONAL :: RESET      ! Reset counter?
! 
! !REVISION HISTORY: 
!  13 May 2009 - C. Carouge  - Initial Version
!  15 Jan 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
      IF ( PRESENT( INCREMENT ) ) THEN
         CT_DIAG = CT_DIAG + 1
      ELSE IF ( PRESENT( RESET ) ) THEN
         CT_DIAG = 0
      ENDIF

      END SUBROUTINE SET_CT_DIAG
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Set_ct_a1
!
! !DESCRIPTION: Subroutine SET\_CT\_A1 increments CT\_A1, the counter of the 
!  number of times we have read in A1 fields.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE SET_CT_A1( INCREMENT, RESET )
!
! !INPUT PARAMETERS:
!
      LOGICAL, INTENT(IN), OPTIONAL :: INCREMENT  ! Increment counter?
      LOGICAL, INTENT(IN), OPTIONAL :: RESET      ! Reset counter?
! 
! !REVISION HISTORY: 
!  19 Aug 2010 - R. Yantosca - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
      IF ( PRESENT( INCREMENT ) ) THEN
         CT_A1 = CT_A1 + 1
      ELSE IF ( PRESENT( RESET ) ) THEN
         CT_A1 = 0
      ENDIF

      END SUBROUTINE SET_CT_A1
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Set_ct_a3
!
! !DESCRIPTION: Subroutine SET\_CT\_A3 increments CT\_A3, the counter of the 
!  number of times we have read in A-3 fields.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE SET_CT_A3( INCREMENT, RESET )
!
! !INPUT PARAMETERS:
!
      LOGICAL, INTENT(IN), OPTIONAL :: INCREMENT  ! Increment counter?
      LOGICAL, INTENT(IN), OPTIONAL :: RESET      ! Reset counter?
! 
! !REVISION HISTORY: 
!  21 Mar 2003 - R. Yantosca - Initial Version
!  15 Jan 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
      IF ( PRESENT( INCREMENT ) ) THEN
         CT_A3 = CT_A3 + 1
      ELSE IF ( PRESENT( RESET ) ) THEN
         CT_A3 = 0
      ENDIF

      END SUBROUTINE SET_CT_A3
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Set_ct_i3
!
! !DESCRIPTION: Subroutine SET\_CT\_I3 increments CT\_I3, the counter of the 
!  number of times we have read in I-3 fields.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE SET_CT_I3( INCREMENT, RESET )
!
! !INPUT PARAMETERS:
!
      LOGICAL, INTENT(IN), OPTIONAL :: INCREMENT  ! Increment counter?
      LOGICAL, INTENT(IN), OPTIONAL :: RESET      ! Reset counter?
! 
! !REVISION HISTORY:
!  03 Feb 2012 - R. Yantosca - Initial version, for GEOS-5.7.2
!EOP
!------------------------------------------------------------------------------
!BOC
      IF ( PRESENT( INCREMENT ) ) THEN
         CT_I3 = CT_I3 + 1
      ELSE IF ( PRESENT( RESET ) ) THEN
         CT_I3 = 0
      ENDIF

      END SUBROUTINE SET_CT_I3
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Set_ct_xtra
!
! !DESCRIPTION: Subroutine SET\_CT\_XTRA increments CT\_XTRA, the counter of 
!  the number of times we have read in GEOS-3 XTRA fields. 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE SET_CT_XTRA( INCREMENT, RESET )
!
! !INPUT PARAMETERS: 
!
      LOGICAL, INTENT(IN), OPTIONAL :: INCREMENT  ! Increment counter?
      LOGICAL, INTENT(IN), OPTIONAL :: RESET      ! Reset counter?
! 
! !REVISION HISTORY: 
!  20 Oct 2009 - T-M Fu, R. Yantosca - Initial Version
!  15 Jan 2010 -         R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
      IF ( PRESENT( INCREMENT ) ) THEN
         CT_XTRA = CT_XTRA + 1
      ELSE IF ( PRESENT( RESET ) ) THEN
         CT_XTRA = 0
      ENDIF

      END SUBROUTINE SET_CT_XTRA
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Set_elapsed_sec
!
! !DESCRIPTION: Subroutine SET\_ELAPSED\_SEC increments the number of elapsed 
!  seconds by the dynamic timestep TS\_DYN.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE SET_ELAPSED_SEC
! 
! !REVISION HISTORY: 
!  05 Feb 2003 - R. Yantosca - Initial Version
!  15 Jan 2010 - R. Yantosca - Added ProTeX headers
!  06 Feb 2018 - E. lundgren - Rename from SET_ELAPSED_MIN to SET_ELAPSED_SEC;
!                              update minimum time unit from min to sec
!EOP
!------------------------------------------------------------------------------
!BOC
      ELAPSED_SEC = ELAPSED_SEC + TS_DYN

      END SUBROUTINE SET_ELAPSED_SEC
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get_jd
!
! !DESCRIPTION: Function GET\_JD is a wrapper for the JULDAY routine.  Given 
!  the current NYMD and NHMS values, GET\_JD will return the current 
!  astronomical Julian date.
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_JD( THISNYMD, THISNHMS ) RESULT( THISJD )
!
! !USES:
!
      USE JULDAY_MOD, ONLY : JULDAY
!
! !INPUT PARAMETERS:
!
      INTEGER, INTENT(IN)  :: THISNYMD   ! YYYY/MM/DD value
      INTEGER, INTENT(IN)  :: THISNHMS   ! hh:mm:ss   value
!
! !RETURN VALUE:
!
      REAL(f8)             :: THISJD     ! Output value
! 
! !REVISION HISTORY: 
!  05 Feb 2003 - R. Yantosca - Initial Version
!  15 Jan 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER  :: Y, M, D, H, MI, S
      REAL(f8) :: DAY

      !=================================================================
      ! GET_JD begins here!
      !=================================================================

      ! Extract year, month, day from NYMDb
      CALL YMD_EXTRACT( THISNYMD, Y, M, D )
         
      ! Extract hour, minute, second from NHMSb
      CALL YMD_EXTRACT( THISNHMS, H, MI, S )      

      ! Decimal day (including fractional part) 
      DAY  = DBLE( D ) + ( DBLE( H  ) /    24e+0_f8 ) + 
     &                   ( DBLE( MI ) /  1440e+0_f8 ) +
     &                   ( DBLE( S  ) / 86400e+0_f8 ) 
     
      ! Compute astronomical Julian day at start of run
      THISJD = JULDAY( Y, M, DAY )

      ! Return to the calling program
      END FUNCTION GET_JD
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get_elapsed_sec
!
! !DESCRIPTION: Function GET\_ELAPSED\_SEC returns the elapsed seconds since 
!  the start of a GEOS-Chem run to the calling program.
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_ELAPSED_SEC() RESULT( THIS_ELAPSED_SEC )
!
! !RETURN VALUE:
!
      INTEGER :: THIS_ELAPSED_SEC
! 
! !REVISION HISTORY: 
!  05 Feb 2003 - R. Yantosca - Initial Version
!  15 Jan 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
      THIS_ELAPSED_SEC = ELAPSED_SEC

      END FUNCTION GET_ELAPSED_SEC
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get_nymdb
!
! !DESCRIPTION: Function GET\_NYMDb returns the NYMDb value (YYYYMMDD at the 
!  beginning of the run).
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_NYMDb() RESULT( THISNYMDb )
!
! !RETURN VALUE:
!
      INTEGER :: THISNYMDb
! 
! !REVISION HISTORY: 
!  05 Feb 2003 - R. Yantosca - Initial Version
!  15 Jan 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
      THISNYMDb = NYMDb

      END FUNCTION GET_NYMDb
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get_nhmsb
!
! !DESCRIPTION: Function GET\_NHMSb returns the NHMSb value (HHMMSS at the 
!  beginning of the run) to the calling program. (bmy, 3/21/03)
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_NHMSb() RESULT( THISNHMSb )
!
! !RETURN VALUE:
!
      INTEGER :: THISNHMSb
! 
! !REVISION HISTORY: 
!  05 Feb 2003 - R. Yantosca - Initial Version
!  15 Jan 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
      THISNHMSb = NHMSb

      END FUNCTION GET_NHMSb
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get_nymde
!
! !DESCRIPTION: Function GET\_NYMDe returns the NYMDe value (YYYYMMDD at the 
!  end of the run) to the calling program. (bmy, 3/21/03)
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_NYMDe() RESULT( THISNYMDe )
!
! !RETURN VALUE:
!
      INTEGER :: THISNYMDe
! 
! !REVISION HISTORY: 
!  05 Feb 2003 - R. Yantosca - Initial Version
!  15 Jan 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
      THISNYMDe = NYMDe

      END FUNCTION GET_NYMDe
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get_nhmse
!
! !DESCRIPTION: Function GET\_NHMSe returns the NHMSe value (HHMMSS at the end
!  of the run).
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_NHMSe() RESULT( THISNHMSe )
!
! !RETURN VALUE:
!
      INTEGER :: THISNHMSe
! 
! !REVISION HISTORY: 
!  05 Feb 2003 - R. Yantosca - Initial Version
!  15 Jan 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
      THISNHMSe = NHMSe

      END FUNCTION GET_NHMSe
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get_nymd
!
! !DESCRIPTION: Function GET\_NYMD returns the current NYMD value (YYYYMMDD).
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_NYMD() RESULT( THISNYMD )
!
! !RETURN VALUE:
!
      INTEGER :: THISNYMD
! 
! !REVISION HISTORY: 
!  05 Feb 2003 - R. Yantosca - Initial Version
!  15 Jan 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
      THISNYMD = NYMD

      END FUNCTION GET_NYMD
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get_nhms
!
! !DESCRIPTION: Function GET\_NHMS returns the current NHMS value (HHMMSS).
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_NHMS() RESULT( THISNHMS )
!
! !RETURN VALUE:
!
      INTEGER :: THISNHMS
! 
! !REVISION HISTORY: 
!  05 Feb 2003 - R. Yantosca - Initial Version
!  15 Jan 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
      THISNHMS = NHMS

      END FUNCTION GET_NHMS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get_time_ahead
!
! !DESCRIPTION: Function GET\_3h\_AHEAD returns to the calling program a 
!  2-element vector containing the YYYYMMDD and HHMMSS values at the current 
!  time plus N\_SECS seconds.
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_TIME_AHEAD( N_SECS ) RESULT( DATE )
!
! !USES:
!
      USE JULDAY_MOD, ONLY : CALDATE
!
! !INPUT PARAMETERS:
!
      INTEGER, INTENT(IN) :: N_SECS   ! Seconds ahead to compute date & time
!
! !RETURN VALUE:
!
      INTEGER             :: DATE(2)  ! Date & time output
! 
! !REVISION HISTORY: 
!  21 Mar 2003 - R. Yantosca - Initial Version
!  (1 ) Bug fix for GCAP leap year case (phs, bmy, 12/8/06)
!  15 Jan 2010 - R. Yantosca - Added ProTeX headers
!  17 Dec 2010 - R. Yantosca - Added fix in case HHMMSS is returned as 240000
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!  06 Feb 2018 - E. Lundgren - Change passed units from minutes to seconds
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
! 
      INTEGER  :: THISYEAR, THISMONTH, THISDAY, TMP
      REAL(f8) :: JD

      !=================================================================
      ! GET_TIME_AHEAD begins here!
      !=================================================================

      ! Astronomical Julian Date at current time + N_SECS
      JD = GET_JD( NYMD, NHMS ) + ( N_SECS / 86400e+0_f8 )

      ! Call CALDATE to compute the current YYYYMMDD and HHMMSS
      CALL CALDATE( DBLE(JD), DATE(1), DATE(2) )

      ! Check to see if HHMMSS is 240000.  This may occur due to a
      ! roundoff error in CALDATE.  If this is the case, then add 1 
      ! to the date and then set HHMMSS = 0.  Use the GET_JD and 
      ! CALDATE functions to do this computation rigorously.
      IF ( DATE(2) == 240000 ) THEN

         ! Split the date into Y/M/D variables
         CALL YMD_EXTRACT( DATE(1), THISYEAR, THISMONTH, THISDAY )

         ! Increment the Astronomical Julian Date by 1 day
         TMP = THISYEAR*10000 + THISMONTH*100 + THISDAY
         JD  = GET_JD( TMP, 000000 ) + 1

         ! Convert to YYYY/MM/DD and hh:mm:ss
         CALL CALDATE( DBLE(JD), DATE(1), DATE(2) )
      ENDIF

      END FUNCTION GET_TIME_AHEAD
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get_month
!
! !DESCRIPTION: Function GET\_MONTH returns the current GMT month.
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_MONTH() RESULT( THISMONTH )
!
! !RETURN VALUE:
!
      INTEGER :: THISMONTH
! 
! !REVISION HISTORY: 
!  05 Feb 2003 - R. Yantosca - Initial Version
!  15 Jan 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
      THISMONTH = MONTH

      END FUNCTION GET_MONTH
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get_day
!
! !DESCRIPTION: Function GET\_DAY returns the current GMT day.
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_DAY() RESULT( THISDAY )
!
! !RETURN VALUE:
!
      INTEGER :: THISDAY
! 
! !REVISION HISTORY: 
!  05 Feb 2003 - R. Yantosca - Initial Version
!  15 Jan 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
      THISDAY = DAY

      END FUNCTION GET_DAY
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get_year
!
! !DESCRIPTION: Function GET\_YEAR returns the current GMT year.
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_YEAR() RESULT( THISYEAR )
!
! !RETURN VALUE:
!
      INTEGER :: THISYEAR
! 
! !REVISION HISTORY: 
!  05 Feb 2003 - R. Yantosca - Initial Version
!  15 Jan 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
      THISYEAR = YEAR

      END FUNCTION GET_YEAR
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Set_histyr
!
! !DESCRIPTION: Function SET\_HISTYR returns the year stored in HISTYR w/o
!  needing to include the CMN\_O3 common block (eml, 8/20/08)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE SET_HISTYR( YEARIN )
!
! !INPUT PARAMETERS:
!
      INTEGER, INTENT(IN) :: YEARIN
! 
! !REVISION HISTORY:
!  20 Aug 2008 - E. Leibensperger - Initial version
!  07 Mar 2012 - M. Payer         - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
      HISTYRA = YEARIN
      
      END SUBROUTINE SET_HISTYR
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get_histyr
!
! !DESCRIPTION: Function GET\_HISTYR returns the year stored in HISTYR w/o
!  needing to include the CMN\_O3 common block (eml, 8/20/08)
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_HISTYR() RESULT( HISTYEAR )
!
! !RETURN VALUE:
!
      INTEGER :: HISTYEAR
! 
! !REVISION HISTORY:
!  20 Aug 2008 - E. Leibensperger - Initial version
!  07 Mar 2012 - M. Payer         - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
      HISTYEAR = HISTYRA

      END FUNCTION GET_HISTYR
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get_hour
!
! !DESCRIPTION: Function GET\_HOUR returns the current GMT hour.
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_HOUR() RESULT( THISHOUR )
!
! !RETURN VALUE:
!
      INTEGER :: THISHOUR
! 
! !REVISION HISTORY: 
!  05 Feb 2003 - R. Yantosca - Initial Version
!  15 Jan 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
      THISHOUR = HOUR

      END FUNCTION GET_HOUR
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get_minute
!
! !DESCRIPTION: Function GET\_MINUTE returns the current GMT minutes.
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_MINUTE() RESULT( THISMINUTE )
!
! !RETURN VALUE:
!
      INTEGER :: THISMINUTE
! 
! !REVISION HISTORY: 
!  05 Feb 2003 - R. Yantosca - Initial Version
!  15 Jan 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
      THISMINUTE = MINUTE

      END FUNCTION GET_MINUTE
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get_second
!
! !DESCRIPTION: Function GET\_SECOND returns the current GMT seconds.
!  calling program.
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_SECOND() RESULT( THISSECOND )
!
! !RETURN VALUE:
!
      INTEGER :: THISSECOND
! 
! !REVISION HISTORY: 
!  05 Feb 2003 - R. Yantosca - Initial Version
!  15 Jan 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
      THISSECOND = SECOND

      END FUNCTION GET_SECOND
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get_day_of_year
!
! !DESCRIPTION: Function GET\_DAY\_OF\_YEAR returns the current day of the
!  year (0-365 or 0-366 for leap years) to the calling program.
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_DAY_OF_YEAR() RESULT( THISDAYOFYEAR )
!
! !RETURN VALUE:
!
      INTEGER :: THISDAYOFYEAR  ! Day of year
! 
! !REVISION HISTORY: 
!  05 Feb 2003 - R. Yantosca - Initial Version
!  15 Jan 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
      THISDAYOFYEAR = DAY_OF_YEAR

      END FUNCTION GET_DAY_OF_YEAR
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get_day_of_week
!
! !DESCRIPTION: Function GET\_DAY\_OF\_WEEK returns the day of the week (with
!  respect to GREENWICH MEAN TIME) as a number: Sun=0, Mon=1, Tue=2, Wed=3,
!  Thu=4, Fri=5, Sat=6.
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_DAY_OF_WEEK() RESULT( DAY_NUM )
!
! !USES:
!
      USE JULDAY_MOD, ONLY : JULDAY
!
! !RETURN VALUE:
!
      INTEGER :: DAY_NUM   ! Day number of week
! 
! !REVISION HISTORY: 
!  05 Feb 2003 - R. Yantosca - Initial Version
!  15 Jan 2010 - R. Yantosca - Added ProTeX headers
!  14 Jun 2013 - R. Yantosca - Now move computation to SET_CURRENT_TIME
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
! 
      ! Return the day of week w/r/t GMT date
      ! (stored in module variable DAY_OF_WEEK)
      DAY_NUM = DAY_OF_WEEK

      END FUNCTION GET_DAY_OF_WEEK
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get_day_of_week_lt
!
! !DESCRIPTION: Function GET\_DAY\_OF\_WEEK\_LT returns the day of the week
!  (with repect to the SOLAR LOCAL TIME AT GRID BOX [I,J,L]) as a number:
!  Sun=0, Mon=1, Tue=2, Wed=3, Thu=4, Fri=5, Sat=6.
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_DAY_OF_WEEK_LT( I, J, L ) RESULT( DAY_NUM )
!
! !USES:
!
      USE GC_GRID_MOD, ONLY : GET_XMID
!
! !INPUT PARAMETERS:
!
      INTEGER, INTENT(IN) :: I         ! Grid box lon index
      INTEGER, INTENT(IN) :: J         ! Grid box lat index
      INTEGER, INTENT(IN) :: L         ! Grid box level index
!
! !RETURN VALUE:
!
      INTEGER             :: DAY_NUM   ! Day of week, w/r/t local time
! 
! !REMARKS:
!  This routine is used by various emissions routines, in order to determine
!  whether weekday or weekend emissions need to be applied.
!
! !REVISION HISTORY: 
!  13 Jun 2013 - R. Yantosca - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL(f8) :: DOW_LOCAL, LOCAL_TIME

      !=================================================================
      ! GET_DAY_OF_WEEK begins here!
      !=================================================================

      ! Get the local time [hrs] at grid box (I,J,L),
      ! but do not force it into the range of [0..24] hours.
      LOCAL_TIME = GMT + ( GET_XMID( I, J, L ) / 15e+0_f8 )

      ! Add the local fraction of day elapsed to the day of week 
      ! obtained w/r/t the GREENWICH MEAN TIME.  This accounts for 
      ! the distance of grid box (I,J,L) from the Greenwich meridian.
      DOW_LOCAL  = DBLE( DAY_OF_WEEK ) + ( LOCAL_TIME / 24e+0_f8 )

      ! If DOW_LOCAL is negative, then we need to add 7.
      ! This means that the day of week w/r/t GMT is Sunday
      ! but that it is still Saturday at grid box (I,J,L).
      IF ( DOW_LOCAL < 0e+0_f8 ) THEN
         DOW_LOCAL = DOW_LOCAL + 7e+0_f8
      ENDIF

      ! Force the day of week w/r/t LOCAL TIME to fall 
      ! within the range 0 (Sunday) thru 6 (Saturday).
      DOW_LOCAL  = MOD( DOW_LOCAL, 7e+0_f8 ) 

      ! Cast to integer and return result
      DAY_NUM    = INT( DOW_LOCAL )
    
      END FUNCTION GET_DAY_OF_WEEK_LT
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get_gmt
!
! !DESCRIPTION:  Function GET\_GMT returns the current Greenwich Mean Time 
!  to the calling program.
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_GMT() RESULT( THISGMT )
!
! !RETURN VALUE:
!
      REAL(f8) :: THISGMT   ! Greenwich mean time [hrs]
! 
! !REVISION HISTORY: 
!  05 Feb 2003 - R. Yantosca - Initial Version
!  15 Jan 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
      THISGMT = GMT

      END FUNCTION GET_GMT
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get_tau
!
! !DESCRIPTION: Function GET\_TAU returns TAU (hours since 1 Jan 
!  1985 at the start of a GEOS-Chem run) to the calling program. 
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_TAU() RESULT( THISTAU )
!
! !RETURN VALUE:
!
      REAL(f8) :: THISTAU  ! TAUb [hrs since 1/1/1985]
! 
! !REVISION HISTORY: 
!  05 Feb 2003 - R. Yantosca - Initial Version
!  15 Jan 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
      THISTAU = TAU

      END FUNCTION GET_TAU
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get_taub
!
! !DESCRIPTION: Function GET\_TAUb returns TAUb (hours since 1 Jan 1985 
!  at the start of a GEOS-Chem run) to the calling program. 
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_TAUb() RESULT( THISTAUb )
!
! !RETURN VALUE:
!
      REAL(f8) :: THISTAUb  ! TAUb [hrs since 1/1/1985]
! 
! !REVISION HISTORY: 
!  05 Feb 2003 - R. Yantosca - Initial Version
!  15 Jan 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
      THISTAUb = TAUb

      END FUNCTION GET_TAUb
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get_taue
!
! !DESCRIPTION: Function GET\_TAUe returns TAUe (hours since 1 Jan 1985 
!  at the end of a GEOS-Chem run) to the calling program. 
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_TAUe() RESULT( THISTAUe )
!
! !RETURN VALUE:
!
      REAL(f8) :: THISTAUe  ! TAUe [hrs since 1/1/1985]
! 
! !REVISION HISTORY: 
!  05 Feb 2003 - R. Yantosca - Initial Version
!  15 Jan 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
      THISTAUe = TAUe

      END FUNCTION GET_TAUe
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get_diagb
!
! !DESCRIPTION: Function GET\_DIAGb returns DIAGb (hours since 1 Jan 1985
!  at the start of a diagnostic interval) to the calling program.
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_DIAGb() RESULT( THISDIAGb )
!
! !RETURN VALUE:
!
      REAL(f8) :: THISDIAGb   ! DIAGb [hrs sincd 1/1/1985]
! 
! !REVISION HISTORY: 
!  05 Feb 2003 - R. Yantosca - Initial Version
!  15 Jan 2010 - R. Yantosca - Added ProTeX headers
!  29 Jul 2014 - R. Yantosca - Bug fix: THISDIAGb needs to be REAL(fp)
!EOP
!------------------------------------------------------------------------------
!BOC
      THISDIAGb = DIAGb

      END FUNCTION GET_DIAGb
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get_diage
!
! !DESCRIPTION: Function GET\_DIAGe returns DIAGe (hours since 1 Jan 1985
!  at the end of a diagnostic interval) to the calling program.
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_DIAGe() RESULT( THISDIAGe )
!
! !RETURN VALUE:
!
      REAL(f8) :: THISDIAGe   ! DIAGe [hrs sincd 1/1/1985]
! 
! !REVISION HISTORY: 
!  05 Feb 2003 - R. Yantosca - Initial Version
!  15 Jan 2010 - R. Yantosca - Added ProTeX headers
!  29 Jul 2014 - R. Yantosca - Bug fix: THISDIAGe needs to be REAL(fp)
!EOP
!------------------------------------------------------------------------------
!BOC
      THISDIAGe = DIAGe

      END FUNCTION GET_DIAGe
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get_localtime
!
! !DESCRIPTION: Function GET\_LOCALTIME returns the local time of a grid
!  box to the calling program. (bmy, 2/5/03)
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_LOCALTIME( I, J, L, GMT ) RESULT( THISLOCALTIME )
!
! !USES:
!
      USE GC_GRID_MOD, ONLY : GET_XMID
!
! !INPUT PARAMETERS:
!
      INTEGER,  INTENT(IN)           :: I              ! Longitude index
      INTEGER,  INTENT(IN)           :: J              ! Latitude index
      INTEGER,  INTENT(IN)           :: L              ! Level index
      REAL(f8), INTENT(IN), OPTIONAL :: GMT            ! GMT time of day [hrs]
!
! !RETURN VALUE:
!
      REAL(f8)                       :: THISLOCALTIME  ! Local time [hrs]
! 
! !REMARKS:
!  Local Time = GMT + ( longitude / 15 ) since each hour of time
!  corresponds to 15 degrees of longitude on the globe
!
! !REVISION HISTORY: 
!  05 Feb 2003 - R. Yantosca - Initial Version
!  15 Jan 2010 - R. Yantosca - Added ProTeX headers
!  27 Apr 2010 - R. Yantosca - Add OFFSET to argument list, to allow the
!                              local time to be computed at an arbitrary time
!                              (e.g. at the halfway point of an interval)
!  05 Oct 2011 - R. Yantosca - Now add GMT as an optional argument
!  07 Oct 2011 - R. Yantosca - Removed obsolete OFFSET argument
!  01 Mar 2012 - R. Yantosca - Now use GET_XMID(I,J,L) from grid_mod.F90, and
!                              add J, L indices to the argument list
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL(f8) :: GMT_HRS

      ! Save the value of the argument GMT in a local variable
      ! If not passed, then use the current GMT time from time_mod.f
      IF ( PRESENT( GMT ) ) THEN
         GMT_HRS = GMT
      ELSE
         GMT_HRS = GET_GMT()
      ENDIF

      ! Local time  = GMT time [hrs] + longitude / 15
      THISLOCALTIME = GMT_HRS + ( GET_XMID( I, J, L ) / 15e+0_f8 ) 

      ! Make sure that THISLOCALTIME is in the range 0-24 hours
      IF ( THISLOCALTIME > 24 ) THISLOCALTIME = THISLOCALTIME - 24e+0_f8
      IF ( THISLOCALTIME < 0  ) THISLOCALTIME = THISLOCALTIME + 24e+0_f8

      ! Return to calling program
      END FUNCTION GET_LOCALTIME
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get_season
!
! !DESCRIPTION: Function GET\_SEASON returns the climatological season number 
!  (1=DJF, 2=MAM, 3=JJA, 4=SON) to the calling program.
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_SEASON() RESULT( THISSEASON )
!
! !RETURN VALUE:
!
      INTEGER :: THISSEASON   ! Current season
! 
! !REVISION HISTORY: 
!  21 Mar 2003 - R. Yantosca - Initial Version
!  15 Jan 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
      THISSEASON = NSEASON

      END FUNCTION GET_SEASON
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get_ts_chem
!
! !DESCRIPTION: Function GET\_TS\_CHEM returns the chemistry timestep in 
!  seconds.
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_TS_CHEM() RESULT( THIS_TS_CHEM )
!
! !RETURN VALUE:
! 
      INTEGER :: THIS_TS_CHEM   ! ! Chemistry timestep [sec]
! 
! !REVISION HISTORY: 
!  21 Mar 2003 - R. Yantosca - Initial Version
!  15 Jan 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
      THIS_TS_CHEM = TS_CHEM

      END FUNCTION GET_TS_CHEM
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get_ts_rad
!
! !DESCRIPTION: Function GET\_TS\_RAD returns the radiation timestep in
!  seconds.
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_TS_RAD() RESULT( THIS_TS_RAD )
!
! !RETURN VALUE:
!
      INTEGER :: THIS_TS_RAD   ! ! Radiation timestep [sec]
!
! !REVISION HISTORY:
!  06 Oct 2012 - D. Ridley   - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
      THIS_TS_RAD = TS_RAD

      END FUNCTION GET_TS_RAD
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get_ts_conv
!
! !DESCRIPTION: Function GET\_TS\_CONV returns the convection timestep in 
!  seconds.
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_TS_CONV() RESULT( THIS_TS_CONV )
!
! !RETURN VALUE:
!
      INTEGER :: THIS_TS_CONV   ! Convective timestep [sec]
! 
! !REVISION HISTORY: 
!  21 Mar 2003 - R. Yantosca - Initial Version
!  15 Jan 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
      THIS_TS_CONV = TS_CONV

      END FUNCTION GET_TS_CONV
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get_ts_diag
!
! !DESCRIPTION: Function GET\_TS\_DIAG returns the diagnostic timestep in 
!  seconds.
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_TS_DIAG() RESULT( THIS_TS_DIAG )
!
! !RETURN VALUE:
!
      INTEGER :: THIS_TS_DIAG   ! Diagnostic timestep [sec]
! 
! !REVISION HISTORY: 
!  21 Mar 2003 - R. Yantosca - Initial Version
!  15 Jan 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
      THIS_TS_DIAG = TS_DIAG

      END FUNCTION GET_TS_DIAG
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get_ts_dyn
!
! !DESCRIPTION: Function GET\_TS\_DIAG returns the diagnostic timestep in 
!  seconds.
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_TS_DYN() RESULT( THIS_TS_DYN )
!
! !RETURN VALUE:
!
      INTEGER :: THIS_TS_DYN    ! Dynamic timestep [sec]
! 
! !REVISION HISTORY: 
!  21 Mar 2003 - R. Yantosca - Initial Version
!  15 Jan 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
      THIS_TS_DYN = TS_DYN

      END FUNCTION GET_TS_DYN
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get_ts_emis
!
! !DESCRIPTION: Function GET\_TS\_EMIS returns the emission timestep in 
!  seconds.
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_TS_EMIS() RESULT( THIS_TS_EMIS )
!
! !RETURN VALUE:
!
      INTEGER :: THIS_TS_EMIS   ! Emissions timestep [sec]
! 
! !REVISION HISTORY: 
!  21 Mar 2003 - R. Yantosca - Initial Version
!  15 Jan 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
      THIS_TS_EMIS = TS_EMIS

      END FUNCTION GET_TS_EMIS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get_ts_unit
!
! !DESCRIPTION: Function GET\_TS\_UNIT returns the unit-conversion timestep
!  in seconds.
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_TS_UNIT() RESULT( THIS_TS_UNIT )
!
! !RETURN VALUE:
!
      INTEGER :: THIS_TS_UNIT   ! Unit conversion timestep [sec]
! 
! !REVISION HISTORY: 
!  21 Mar 2003 - R. Yantosca - Initial Version
!  15 Jan 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
      THIS_TS_UNIT = TS_UNIT

      END FUNCTION GET_TS_UNIT
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get_ct_chem
!
! !DESCRIPTION: Function GET\_CT\_CHEM returns the chemistry timestep counter 
!  to the calling program.
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_CT_CHEM() RESULT( THIS_CT_CHEM )
!
! !RETURN VALUE:
!
      INTEGER :: THIS_CT_CHEM
! 
! !REVISION HISTORY: 
!  21 Mar 2003 - R. Yantosca - Initial Version
!  15 Jan 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
      THIS_CT_CHEM = CT_CHEM

      END FUNCTION GET_CT_CHEM
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get_ct_rad
!
! !DESCRIPTION: Function GET\_CT\_RAD returns the radiation timestep counter
!  to the calling program.
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_CT_RAD() RESULT( THIS_CT_RAD )
!
! !RETURN VALUE:
!
      INTEGER :: THIS_CT_RAD
!
! !REVISION HISTORY:
!  06 Oct 2012 - D. Ridley   - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
      THIS_CT_RAD = CT_RAD

      END FUNCTION GET_CT_RAD
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get_ct_conv
!
! !DESCRIPTION: Function GET\_CT\_CONV returns the convection timestep 
!  counter to the calling program.
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_CT_CONV() RESULT( THIS_CT_CONV )
!
! !RETURN VALUE:
!
      INTEGER :: THIS_CT_CONV   ! # of convection timesteps
! 
! !REVISION HISTORY: 
!  21 Mar 2003 - R. Yantosca - Initial Version
!  15 Jan 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
      THIS_CT_CONV = CT_CONV

      END FUNCTION GET_CT_CONV
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get_ct_dyn
!
! !DESCRIPTION: Function GET\_CT\_CHEM returns the dynamic timestep counter 
!  to the calling program.
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_CT_DYN() RESULT( THIS_CT_DYN )
!
! !RETURN VALUE:
!
      INTEGER :: THIS_CT_DYN   ! # of dynamics timesteps
! 
! !REVISION HISTORY: 
!  21 Mar 2003 - R. Yantosca - Initial Version
!  15 Jan 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
      THIS_CT_DYN = CT_DYN

      END FUNCTION GET_CT_DYN
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get_ct_emis
!
! !DESCRIPTION: Function GET\_CT\_CHEM returns the emissions timestep counter
!  to the calling program.
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_CT_EMIS() RESULT( THIS_CT_EMIS )
!
! !RETURN VALUE:
!
      INTEGER :: THIS_CT_EMIS  ! # of emissions timesteps
!
! 
! !REVISION HISTORY: 
!  21 Mar 2003 - R. Yantosca - Initial Version
!  15 Jan 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
      THIS_CT_EMIS = CT_EMIS

      END FUNCTION GET_CT_EMIS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get_ct_a1
!
! !DESCRIPTION: Function GET\_CT\_A1 returns the A1 fields timestep 
!  counter to the calling program.
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_CT_A1() RESULT( THIS_CT_A1 )
!
! !RETURN VALUE:
!
      INTEGER :: THIS_CT_A1   ! # of A-3 timesteps
! 
! !REVISION HISTORY: 
!  19 Aug 2010 - R. Yantosca - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
      THIS_CT_A1 = CT_A1

      END FUNCTION GET_CT_A1
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get_ct_a3
!
! !DESCRIPTION: Function GET\_CT\_A3 returns the A-3 fields timestep 
!  counter to the calling program.
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_CT_A3() RESULT( THIS_CT_A3 )
!
! !RETURN VALUE:
!
      INTEGER :: THIS_CT_A3   ! # of A-3 timesteps
! 
! !REVISION HISTORY: 
!  21 Mar 2003 - R. Yantosca - Initial Version
!  15 Jan 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
      THIS_CT_A3 = CT_A3

      END FUNCTION GET_CT_A3
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get_ct_i3
!
! !DESCRIPTION: Function GET\_CT\_I3 returns the I-3 fields timestep counter
!  to the calling program
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_CT_I3() RESULT( THIS_CT_I3 )
!
! !RETURN VALUE:
!
      INTEGER :: THIS_CT_I3   ! # of I-6 timesteps
! 
! !REVISION HISTORY: 
!  03 Feb 2012 - R. Yantosca - Initial version, for GEOS-5.7.2
!EOP
!------------------------------------------------------------------------------
!BOC
      THIS_CT_I3 = CT_I3

      END FUNCTION GET_CT_I3
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get_ct_xtra
!
! !DESCRIPTION: Function GET\_CT\_XTRA returns the XTRA fields timestep 
!  counter to the calling program.
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_CT_XTRA() RESULT( THIS_CT_XTRA )
!
! !RETURN VALUE:
!
      INTEGER :: THIS_CT_XTRA    ! # of XTRA timesteps
! 
! !REVISION HISTORY: 
!  20 Oct 2005 - T-M Fu, R. Yantosca - Initial Version
!  15 Jan 2010 -         R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
      THIS_CT_XTRA = CT_XTRA

      END FUNCTION GET_CT_XTRA
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get_ct_diag
!
! !DESCRIPTION: Function GET\_CT\_DIAG returns the DIAG timestep counter to the
!  calling program.
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_CT_DIAG() RESULT( THIS_CT_DIAG )
!
! !RETURN VALUE:
!
      INTEGER :: THIS_CT_DIAG   ! # of diagnostic timesteps
! 
! !REVISION HISTORY: 
!  21 May 2009 - C. Carouge  - Initial Version
!  15 Jan 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
      THIS_CT_DIAG = CT_DIAG

      END FUNCTION GET_CT_DIAG
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get_hg2_diag
!
! !DESCRIPTION: Function GET\_Hg2\_DIAG returns the DIAG timestep counter to
!  the calling program. (hma 20100218)
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_Hg2_DIAG() RESULT( THIS_Hg2_DIAG )
!
! !RETURN VALUE:
!
      INTEGER :: THIS_Hg2_DIAG  ! # of diagnostic timesteps
! 
! !REVISION HISTORY:
!  18 Feb 2012 - H. Amos     - Initial version 
!  07 Mar 2012 - M. Payer    - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
      THIS_Hg2_DIAG = Hg2_DIAG

      END FUNCTION GET_Hg2_DIAG
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get_a1_time
!
! !DESCRIPTION: Function GET\_A1\_TIME returns the correct YYYYMMDD and HHMMSS 
!  values that are needed to read in the next average 1-hour (A-1) fields. 
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_A1_TIME() RESULT( DATE )
!
! !RETURN VALUE:
!
      INTEGER :: DATE(2)   ! YYYYMMDD and HHMMSS values
! 
! !REVISION HISTORY: 
!  19 Aug 2010 - R. Yantosca - Initial version
!  02 Feb 2012 - R. Yantosca - Added modifications for GEOS-5.7.x met fields
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!  26 Sep 2013 - R. Yantosca - Renamed GEOS_57 Cpp switch to GEOS_FP
!  11 Aug 2015 - R. Yantosca - Return DATE for MERRA2 the same as for GEOS-FP
!  06 Feb 2017 - E. Lundgren - Update GET_TIME_AHEAD arg to pass sec not min
!EOP
!------------------------------------------------------------------------------
!BOC

      ! GEOS-FP and MERRA-2 met fields are 1-hour time-averages, timestamped 
      ! at the center of the averaging periods (00:30, 01:30, 02:30 ... 23:30)
      DATE = GET_TIME_AHEAD( 1800 )

      END FUNCTION GET_A1_TIME
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get_a3_time
!
! !DESCRIPTION: Function GET\_A3\_TIME returns the correct YYYYMMDD and HHMMSS 
!  values that are needed to read in the next average 3-hour (A-3) fields. 
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_A3_TIME() RESULT( DATE )
!
! !RETURN VALUE:
!
      INTEGER :: DATE(2)   ! YYYYMMDD and HHMMSS values
! 
! !REVISION HISTORY: 
!  21 Mar 2003 - R. Yantosca - Initial Version
!  (1 ) Now return proper time for GEOS-4/fvDAS fields (bmy, 6/19/03)
!  (2 ) Remove reference to FIRST variable (bmy, 12/10/04)
!  (3 ) Now modified for GCAP and GEOS-5 met fields (swu, bmy, 5/24/05)
!  (4 ) Remove support for GEOS-1 and GEOS-STRAT met fields (bmy, 8/4/06)
!  15 Jan 2010 - R. Yantosca - Added ProTeX headers
!  28 Feb 2012 - R. Yantosca - Removed support for GEOS-3
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!  06 Feb 2017 - E. Lundgren - Update GET_TIME_AHEAD arg to pass sec not min
!EOP
!------------------------------------------------------------------------------
!BOC

      ! The A-3 fields are timestamped by center time.
      ! Therefore, the difference between the actual time when the fields
      ! are read and the A-3 timestamp time is 90 minutes.
      DATE = GET_TIME_AHEAD( 5400 )   

      END FUNCTION GET_A3_TIME
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get_i3_time
!
! !DESCRIPTION: Function GET\_I3\_TIME returns the correct YYYYMMDD and 
!  HHMMSS values that are needed to read in the next instantaneous 3-hour 
!  (I-3) fields. 
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_I3_TIME() RESULT( DATE )
!
! !RETURN VALUE:
!
      INTEGER :: DATE(2)   ! YYYYMMDD and HHMMSS values
!  
! !REMARKS:
!  Modified for start times other than 0 GMT.
!
! !REVISION HISTORY: 
!  06 Feb 2012 - R. Yantosca - Initial version
!  06 Feb 2018 - E. Lundgren - Update GET_TIME_AHEAD arg to pass sec not min
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
! 
      LOGICAL, SAVE :: FIRST = .TRUE.
      INTEGER       :: HH, MM, SS, SECS, OFFSET

      !=================================================================
      ! ALL MET FIELDS:
      !=================================================================
 
      IF ( FIRST ) THEN

         !--------------------------------------------------------------
         ! FIRST-TIME ONLY!  Get the proper # of hours until the next 
         ! I3 time.  Also works for start times other than 0 GMT.
         !--------------------------------------------------------------

         ! Split NHMS into hours, mins, seconds
         CALL YMD_EXTRACT( NHMS, HH, MM, SS )

         ! Compute seconds elapsed in the 3-hour interval
         SECS   = MOD( HH, 3 )*3600 + MM*60 + SS

         ! Compute offset to next I-3 time
         OFFSET = 10800 - SECS

         ! Get YYYY/MM/DD and hh:mm:ss to next I3 time
         DATE   = GET_TIME_AHEAD( OFFSET ) 

         ! Reset first-time flag
         FIRST = .FALSE.

      ELSE

         !--------------------------------------------------------------
         ! Other than the 1st time: Search 180 mins ahead
         !--------------------------------------------------------------

         ! We need to read in the I-3 fields 3h (180 mins, or 10800 secs) 
         ! ahead of time
         DATE = GET_TIME_AHEAD( 10800 )

      ENDIF

      END FUNCTION GET_I3_TIME
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get_first_a1_time
!
! !DESCRIPTION: Function GET\_FIRST\_A1\_TIME returns the correct YYYYMMDD
!  and HHMMSS values the first time that A-3 fields are read in from disk. 
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_FIRST_A1_TIME() RESULT( DATE )
!
! !RETURN VALUE:
!
      INTEGER :: DATE(2)   ! YYYYMMDD and HHMMSS values
! 
! !REVISION HISTORY: 
!  26 Jun 2003 - R. Yantosca - Initial Version
!  (1 ) Now modified for GCAP and GEOS-5 data (swu, bmy, 5/24/05) 
!  (2 ) Remove support for GEOS-1 and GEOS-STRAT met fields (bmy, 8/4/06)
!  15 Jan 2010 - R. Yantosca - Added ProTeX headers
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!EOP
!------------------------------------------------------------------------------
!BOC
!
      ! Return the first A-1 time
      DATE = GET_A1_TIME()
  
      END FUNCTION GET_FIRST_A1_TIME
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get_first_a3_time
!
! !DESCRIPTION: Function GET\_FIRST\_A3\_TIME returns the correct YYYYMMDD
!  and HHMMSS values the first time that A-3 fields are read in from disk. 
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_FIRST_A3_TIME() RESULT( DATE )
!
! !RETURN VALUE:
!
      INTEGER :: DATE(2)   ! YYYYMMDD and HHMMSS values
! 
! !REVISION HISTORY: 
!  26 Jun 2003 - R. Yantosca - Initial Version
!  (1 ) Now modified for GCAP and GEOS-5 data (swu, bmy, 5/24/05) 
!  (2 ) Remove support for GEOS-1 and GEOS-STRAT met fields (bmy, 8/4/06)
!  15 Jan 2010 - R. Yantosca - Added ProTeX headers
!  27 Sep 2010 - R. Yantosca - Modified for start times other than 0 GMT
!  28 Feb 2012 - R. Yantosca - Removed support for GEOS-3
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!  06 Feb 2018 - E. Lundgren - Update GET_TIME_AHEAD arg to pass sec not min
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER :: HH, MM, SS, SECS, OFFSET

      !==================================================================
      ! A3 fields are indexed at the midpoint of the 3-hr interval
      !==================================================================

      ! Split NYMS into hours, mins, seconds
      CALL YMD_EXTRACT( NHMS, HH, MM, SS )

      ! Compute seconds elapsed in the 3-hour interval
      SECS   = MOD( HH, 3 )*3600 + MM*60 + SS
      
      ! Compute offset to midpoint of 3hr interval
      OFFSET = 10800 - ( SECS + 5400 )

      ! Get YYYY/MM/DD and hh:mm:ss at midpoint of 3hr interval
      DATE   = GET_TIME_AHEAD( OFFSET )

      END FUNCTION GET_FIRST_A3_TIME
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get_first_i3_time
!
! !DESCRIPTION: Function GET\_FIRST\_I3\_TIME returns the correct YYYYMMDD and
!  HHMMSS values the first time that I-3 fields are read in from disk.
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_FIRST_I3_TIME() RESULT( DATE )
!
! !RETURN VALUE:
!
      INTEGER :: DATE(2)    ! YYYYMMDD, HHMMSS values
! 
! !REVISION HISTORY: 
!  03 Feb 2012 - R. Yantosca - Initial version, for GEOS-5.7.2
!  06 Feb 2018 - E. Lundgren - Update GET_TIME_AHEAD arg to pass sec not min
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER :: HH, MM, SS, SECS, OFFSET

      !==================================================================
      ! Compute first I-6 time for all met field types
      !==================================================================

      ! Split NYMS into hours, mins, seconds
      CALL YMD_EXTRACT( NHMS, HH, MM, SS )        

      ! Compute seconds elapsed in the 3-hour interval
      SECS   = MOD( HH, 3 )*3600 + MM*60 + SS

      ! Compute offset to nearest I-6 time
      OFFSET = -SECS

      ! Get YYYY/MM/DD and hh:mm:ss to nearest I-6 time
      DATE   = GET_TIME_AHEAD( OFFSET ) 

      END FUNCTION GET_FIRST_I3_TIME
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Its_Time_For_chem
!
! !DESCRIPTION: Function ITS\_TIME\_FOR\_CHEM returns TRUE if it is time to do 
!  chemistry, or FALSE otherwise.
!\\
!\\
! !INTERFACE:
!
      FUNCTION ITS_TIME_FOR_CHEM() RESULT( FLAG )
!
! !RETURN VALUE:
!
      LOGICAL :: FLAG 
! 
! !REVISION HISTORY: 
!  21 Mar 2003 - R. Yantosca - Initial Version
!  15 Jan 2010 - R. Yantosca - Added ProTeX headers
!  27 Sep 2011 - M. Payer    - Modifications for centralizing the chemistry
!                              time step (lzh)
!  06 Feb 2018 - E. Lundgren - Change ELAPSED_MIN to ELAPSED_SEC to match new
!                              timestep units of seconds
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! changes for proper chemistry time (lzh, ccc, 3/20/10)
      INTEGER :: M 

      !=================================================================
      ! ITS_TIME_FOR_CHEM begins here!
      !=================================================================

      ! Get multiplier between transport and chemistry:
      M = TS_CHEM/TS_DYN

      ! Divide by 2 (get middle). KEEP INTEGERS!!!!
      M = MAX( M/2, 1 )

      ! Is it time for chemistry?
      ! Chemistry time step is now in the center of transport time steps
      FLAG = ( MOD( ELAPSED_SEC, TS_CHEM ) == (M-1)*TS_DYN )

      END FUNCTION ITS_TIME_FOR_CHEM
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: its_time_for_rt
!
! !DESCRIPTION: Function ITS\_TIME\_FOR\_RT returns TRUE if it is time to do
!  radiative transfer calculations, or FALSE otherwise.
!\\
!\\
! !INTERFACE:
!
      FUNCTION ITS_TIME_FOR_RT() RESULT( FLAG )
!
! !RETURN VALUE:
!
      LOGICAL :: FLAG
!
! !REVISION HISTORY:
!  04 Oct 2012 - D. Ridley   - Initial Version
!  06 Feb 2018 - E. Lundgren - Change ELAPSED_MIN to ELAPSED_SEC to match new
!                              timestep units of seconds
!  18 May 2018 - C. Holmes   - Call radiation the first time after half-interval,
!                              without requiring an exact time match.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER :: M, TMOD, TS

      !=================================================================
      ! ITS_TIME_FOR_RT begins here!
      !=================================================================
      ! Get half a time step value
      M = TS_RAD/2

      ! Elapsed time in this radiation step
      TMOD =  MOD( ELAPSED_SEC, TS_RAD )

      ! Smallest time step
      TS = TS_DYN

      ! It's time for radiation, when the current time is greater than
      ! half step and the previous time was less.
      FLAG = ( TMOD >= M ) .and. ( TMOD-TS < M )

      END FUNCTION ITS_TIME_FOR_RT
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: its_time_for_surface_rad
!
! !DESCRIPTION: Function ITS\_TIME\_FOR\_SURFACE\_RAD returns TRUE if it is
!  time to read surface albedo and emissivity fields, or FALSE otherwise.
!\\
!\\
! !INTERFACE:
!
      FUNCTION ITS_TIME_FOR_SURFACE_RAD() RESULT( FLAG )
!
! !RETURN VALUE:
!
      LOGICAL :: FLAG
!
! !REVISION HISTORY:
!  06 Oct 2012 - D. Ridley   - Initial Version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER :: M

      !=================================================================
      ! ITS_TIME_FOR_SURFACE_RAD begins here!
      !=================================================================
      ! Get half a time step value
      M = GET_DAY_OF_YEAR()

      ! Is it time to read the 8-day file?
      ! files start on 1, then 9, 17... etc
      FLAG = ( MOD( M, 8 ) == 1 )

      END FUNCTION ITS_TIME_FOR_SURFACE_RAD
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Its_Time_For_conv
!
! !DESCRIPTION: Function ITS\_TIME\_FOR\_CONV returns TRUE if it is time to do 
!  convection, or FALSE otherwise.
!\\
!\\
! !INTERFACE:
!
      FUNCTION ITS_TIME_FOR_CONV() RESULT( FLAG )
!
! !RETURN VALUE:
!
      LOGICAL :: FLAG 
! 
! !REVISION HISTORY: 
!  21 Mar 2003 - R. Yantosca - Initial Version
!  15 Jan 2010 - R. Yantosca - Added ProTeX headers
!  06 Feb 2018 - E. Lundgren - Change ELAPSED_MIN to ELAPSED_SEC to match new
!                              timestep units of seconds
!EOP
!------------------------------------------------------------------------------
!BOC
      ! Is it time for convection?
      FLAG = ( MOD( ELAPSED_SEC, TS_CONV ) == 0 )

      END FUNCTION ITS_TIME_FOR_CONV
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Its_Time_For_dyn
!
! !DESCRIPTION: Function ITS\_TIME\_FOR\_DYN returns TRUE if it is time to do 
!  chemistry and false otherwise.
!\\
!\\
! !INTERFACE:
!
      FUNCTION ITS_TIME_FOR_DYN() RESULT( FLAG )
!
! !RETURN VALUE:
!
      LOGICAL :: FLAG 
! 
! !REVISION HISTORY: 
!  21 Mar 2003 - R. Yantosca - Initial Version
!  15 Jan 2010 - R. Yantosca - Added ProTeX headers
!  06 Feb 2018 - E. Lundgren - Change ELAPSED_MIN to ELAPSED_SEC to match new
!                              timestep units of seconds
!EOP
!------------------------------------------------------------------------------
!BOC
      ! Is it time for dynamics?
      FLAG = ( MOD( ELAPSED_SEC, TS_DYN ) == 0 )

      END FUNCTION ITS_TIME_FOR_DYN
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Its_Time_For_emis
!
! !DESCRIPTION: Function ITS\_TIME\_FOR\_EMIS returns TRUE if it is time to do 
!  emissions, or FALSE otherwise. 
!\\
!\\
! !INTERFACE:
!
      FUNCTION ITS_TIME_FOR_EMIS() RESULT( FLAG )
!
! !RETURN VALUE:
!
      LOGICAL :: FLAG 
! 
! !REVISION HISTORY: 
!  21 Mar 2003 - R. Yantosca - Initial Version
!  15 Jan 2010 - R. Yantosca - Added ProTeX headers
!  07 Oct 2011 - R. Yantosca - Modifications for centralizing the chemistry
!                              time step (lzh)
!  06 Feb 2018 - E. Lundgren - Change ELAPSED_MIN to ELAPSED_SEC to match new
!                              timestep units of seconds
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! changes for proper chemistry time (lzh, ccc, 3/20/10)
      INTEGER :: M 

      !=================================================================
      ! ITS_TIME_FOR_EMIS begins here!
      !=================================================================

      ! Get multiplier between transport and chemistry:
      M = TS_EMIS/TS_DYN

      ! Divide by 2 (get middle). KEEP INTEGERS!!!!
      M = MAX( M/2, 1 )

      ! Is it time for emissions?
      ! Emission time step is now in the center of transport time steps
      FLAG = ( MOD( ELAPSED_SEC, TS_EMIS ) == (M-1)*TS_DYN )

      END FUNCTION ITS_TIME_FOR_EMIS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Its_Time_For_exch
!
! !DESCRIPTION: Function ITS\_TIME\_FOR\_EXCH returns TRUE if it is time to do
!  exchange for two-way coupled simulation, or FALSE otherwise.
!\\
!\\
! !INTERFACE:
!
      FUNCTION ITS_TIME_FOR_EXCH() RESULT( FLAG )
!
! !RETURN VALUE:
!
      LOGICAL :: FLAG
!
! !REVISION HISTORY:
!  30 Mar 2014 - Y.Y. Yan    - Initial Version
!  06 Feb 2018 - E. Lundgren - Change ELAPSED_MIN to ELAPSED_SEC to match new
!                              timestep units of seconds
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! changes for proper chemistry time (lzh, ccc, 3/20/10)
      INTEGER :: M

      !=================================================================
      ! ITS_TIME_FOR_CHEM begins here!
      !=================================================================

      ! Get multiplier between transport and chemistry:
      M = TS_CHEM/TS_DYN

      ! Divide by 2 (get middle). KEEP INTEGERS!!!!
      M = MAX( M/2, 1 )

      ! Is it time for exchange?
      ! Chemistry time step is now in the center of transport time steps
      FLAG = ( MOD( ELAPSED_SEC, 10800 ) == 0 )

      END FUNCTION ITS_TIME_FOR_EXCH
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Its_Time_For_unit
!
! !DESCRIPTION: Function ITS\_TIME\_FOR\_UNIT returns TRUE if it is time to do 
!  unit conversion, or FALSE otherwise.
!\\
!\\
! !INTERFACE:
!
      FUNCTION ITS_TIME_FOR_UNIT() RESULT( FLAG )
!
! !RETURN VALUE:
!
      LOGICAL :: FLAG 
! 
! !REVISION HISTORY: 
!  21 Mar 2003 - R. Yantosca - Initial Version
!  15 Jan 2010 - R. Yantosca - Added ProTeX headers
!  06 Feb 2018 - E. Lundgren - Change ELAPSED_MIN to ELAPSED_SEC to match new
!                              timestep units of seconds
!EOP
!------------------------------------------------------------------------------
!BOC

      ! Is it time for unit conversion?
      FLAG = ( MOD( ELAPSED_SEC, TS_DYN ) == 0 )

      END FUNCTION ITS_TIME_FOR_UNIT
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Its_Time_For_diag
!
! !DESCRIPTION: Function ITS\_TIME\_FOR\_DIAG returns TRUE if it is time to 
!  archive certain diagnostics, or FALSE otherwise.
!\\
!\\
! !INTERFACE:
!
      FUNCTION ITS_TIME_FOR_DIAG() RESULT( FLAG )
!
! !RETURN VALUE:
!
      LOGICAL :: FLAG 
! 
! !REVISION HISTORY: 
!  21 Mar 2003 - R. Yantosca - Initial Version
!  20 Jul 2009 - C. Carouge  - Use TS_DIAG now and not 60 minutes
!  15 Jan 2010 - R. Yantosca - Added ProTeX headers
!  06 Feb 2018 - E. Lundgren - Change ELAPSED_MIN to ELAPSED_SEC to match new
!                              timestep units of seconds
!EOP
!------------------------------------------------------------------------------
!BOC
      ! Is it time for diagnostics?
      FLAG = ( MOD( ELAPSED_SEC, TS_DIAG ) == 0 )

      END FUNCTION ITS_TIME_FOR_DIAG
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Its_Time_For_a1
!
! !DESCRIPTION: Function ITS\_TIME\_FOR\_A1 returns TRUE if it is time to read 
!  in A1 (average 1-hr fields) and FALSE otherwise. 
!\\
!\\
! !INTERFACE:
!
      FUNCTION ITS_TIME_FOR_A1() RESULT( FLAG )
!
! !RETURN VALUE:
!
      LOGICAL :: FLAG 
! 
! !REVISION HISTORY: 
!  20 Aug 2010 - R. Yantosca - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
      ! We read A1 fields every 3 hours
      FLAG = ( MOD( NHMS, 010000 ) == 0 )

      END FUNCTION ITS_TIME_FOR_A1
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Its_Time_For_a3
!
! !DESCRIPTION: Function ITS\_TIME\_FOR\_A3 returns TRUE if it is time to read 
!  in A3 (average 3-hr fields) and FALSE otherwise. 
!\\
!\\
! !INTERFACE:
!
      FUNCTION ITS_TIME_FOR_A3() RESULT( FLAG )
!
! !RETURN VALUE:
!
      LOGICAL :: FLAG 
! 
! !REVISION HISTORY: 
!  21 Mar 2003 - R. Yantosca - Initial Version
!  15 Jan 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
      ! We read A3 fields every 3 hours
      FLAG = ( MOD( NHMS, 030000 ) == 0 )

      END FUNCTION ITS_TIME_FOR_A3
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Its_Time_For_i3
!
! !DESCRIPTION: Function ITS\_TIME\_FOR\_I3 returns TRUE if it is time to read 
!  in I2 (instantaneous 3-hr fields) and FALSE otherwise.
!\\
!\\
! !INTERFACE:
!
      FUNCTION ITS_TIME_FOR_I3() RESULT( FLAG )
!
! !RETURN VALUE:
!
      LOGICAL :: FLAG 
! 
! !REVISION HISTORY: 
!  03 Feb 2012 - R. Yantosca - Initial version, for GEOS-5.7.2
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARAIABLES:
!
      LOGICAL, SAVE :: FIRST = .TRUE.
      
      ! We read in I-6 fields at 00, 03, 06, 09, 12, 15, 18, 21 GMT
      FLAG = ( ( MOD( NHMS, 030000 ) == 0 ) .or. FIRST )
      
      FIRST = .FALSE.

      END FUNCTION ITS_TIME_FOR_I3
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Its_Time_For_exit
!
! !DESCRIPTION: Function ITS\_TIME\_FOR\_EXIT returns TRUE if it is the end of
!  the GEOS-Chem simulation (i.e. TAU >= TAUe), or FALSE otherwise.
!\\
!\\
! !INTERFACE:
!
      FUNCTION ITS_TIME_FOR_EXIT() RESULT( FLAG )
!
! !RETURN VALUE:
!
      LOGICAL :: FLAG 
! 
! !REVISION HISTORY: 
!  21 Mar 2003 - R. Yantosca - Initial Version
!  15 Jan 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
      ! Test if it's end of run
      FLAG = ( TAU >= TAUe )

      END FUNCTION ITS_TIME_FOR_EXIT
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Its_Time_For_bpch
!
! !DESCRIPTION: Function ITS\_TIME\_FOR\_BPCH returns TRUE if it's time to 
!  write output to the bpch file, or FALSE otherwise.
!\\
!\\
! !INTERFACE:
!
      FUNCTION ITS_TIME_FOR_BPCH( Input_Opt ) RESULT( DO_BPCH )
!
! !USES:
!
      USE Input_Opt_Mod, ONLY : OptInput
!
! !INPUT PARAMETERS:
!
      TYPE(OptInput), INTENT(IN) :: Input_Opt  ! Input options
!
! !RETURN VALUE:
!
      LOGICAL :: DO_BPCH
! 
! !REMARKS:
!  NJDAY is now located in CMN_SIZE_mod.F, so that we can eventually
!  retire the obsolete CMN_DIAG_mod.F. (bmy, 1/16/18)
!
! !REVISION HISTORY: 
!  02 Feb 2007 - R. Yantosca - Initial Version
!  15 Jan 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC

      ! Local variables
      INTEGER :: DOY, THIS_NJDAY

      !=================================================================
      ! ITS_TIME_FOR_BPCH begins here!
      !================================================================= 

      ! Return FALSE if it's the first timestep
      IF ( TAU == TAUb ) THEN
         DO_BPCH = .FALSE.
         RETURN
      ENDIF

      ! Day of year (0..365 or 0..366 leapyears)
      DOY = DAY_OF_YEAR

      ! Look up appropriate value of NJDAY array.  We may need to add a
      ! day to skip past the Feb 29 element of NJDAY for non-leap-years.
      IF ( .not. ITS_A_LEAPYEAR( FORCE=.TRUE. ) .and. DOY > 59 ) THEN
         THIS_NJDAY = Input_Opt%NJDAY( DOY + 1 ) 
      ELSE
         THIS_NJDAY = Input_Opt%NJDAY( DOY )
      ENDIF

      ! Test if this is the day & time to write to the BPCH file!
      IF ( ( THIS_NJDAY > 0 ) .and. NHMS == NHMSe ) THEN
         DO_BPCH = .TRUE.
      ELSE
         DO_BPCH = .FALSE.
      ENDIF

      END FUNCTION ITS_TIME_FOR_BPCH
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Its_A_LeapYear
!
! !DESCRIPTION: Function ITS\_A\_LEAPYEAR tests to see if a year is really a 
!  leapyear. 
!\\
!\\
! !INTERFACE:
!
      FUNCTION ITS_A_LEAPYEAR( YEAR_IN, FORCE ) RESULT( IS_LEAPYEAR )
!
! !INPUT PARAMETERS: 
!
      INTEGER, INTENT(IN), OPTIONAL :: YEAR_IN   ! Year to test if leapyear
      LOGICAL, INTENT(IN), OPTIONAL :: FORCE     ! Do not exit if using GCAP
!
! !RETURN VALUE:
!
      LOGICAL                       :: IS_LEAPYEAR  ! =T if it's a leapyear
! 
! !REVISION HISTORY: 
!  17 Mar 1999 - R. Yantosca - Initial Version
!  (1 ) Now remove YEAR from ARG list; use the module variable (bmy, 3/21/03)
!  (2 ) Now add YEAR_IN as an optional argument.  If YEAR_IN is not passed,
!        then test if the current year is a leapyear (bmy, 9/25/03)
!  (3 ) Now always return FALSE for GCAP (swu, bmy, 8/29/05)
!  (4 ) Now add FORCE argument to force ITS_A_LEAPYEAR to return a value
!        instead of just returning with FALSE for the GCAP met fields.
!        (swu, bmy, 4/24/06)
!  15 Jan 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER  :: THISYEAR
      LOGICAL  :: THISFORCE

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

      IF ( MOD( THISYEAR, 4 ) == 0 ) THEN
         IF ( MOD( THISYEAR, 100 ) == 0 ) THEN
            IF ( MOD( THISYEAR, 400 ) == 0 ) THEN
               IS_LEAPYEAR = .TRUE.
            ENDIF
         ELSE
            IS_LEAPYEAR = .TRUE.
         ENDIF
      ENDIF        

      END FUNCTION ITS_A_LEAPYEAR
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Its_A_New_Year
!
! !DESCRIPTION: Function ITS\_A\_NEW\_YEAR returns TRUE if it's the first 
!  timestep of the year when we have to read in annual data.
!\\
!\\
! !INTERFACE:
!
      FUNCTION ITS_A_NEW_YEAR( NO_CCTS ) RESULT( IS_NEW_YEAR )
!
! !INPUT PARAMETERS:
!
      LOGICAL, OPTIONAL :: NO_CCTS       ! =T reverts to previous behavior
                                         ! (i.e. w/o using central chem step)
!
! !RETURN VALUE:
!
      LOGICAL           :: IS_NEW_YEAR   ! =T if it's 1st data read of year
! 
! !REMARKS:
!  ITS_A_NEW_YEAR assumes that we are using the central chemistry timestep
!  option (i.e. do chemistry & emissions & related processes at the midpoint
!  of each chemistry timestep).  To revert to the prior behavior, set the
!  optional flag NO_CCTS = .TRUE.
!                                                                             .
!  If we are using the central chemistry timestep option (which is now the
!  default behavior), then we must not read data at 00:00 GMT on the first day 
!  of the year, but at the center of the first chemistry timestep of the 
!  year.  This is because emissions and chemistry are done at the same time.  
!  The proper time of day for reading emissions is determined by function 
!  ITS_TIME_FOR_EMIS, also within time_mod.f.
!                                                                             .
!  Similarly, for simulations that start at an arbitrary midmonth date and 
!  time, we must not read data at the starting date and time of the simulation,
!  but at the midpoint of the first chemistry timestep of the simulation.
!                                                                             .
!  If we are not using the central chemistry timestep option (specified by
!  NO_CCTS=.TRUE.), then the first data read of the month occurs at 00:00 GMT 
!  on the Jan 1st.  Similarly, for those simulations that start at midmonth, 
!  the first data read will occur the starting date and time
!  of the simulation.
!
! !REVISION HISTORY: 
!  01 Apr 2004 - R. Yantosca - Initial Version
!  01 Nov 2005 - R. Yantosca - Bug fix: Need month & day to be 1
!  15 Jan 2010 - R. Yantosca - Added ProTeX headers
!  14 Oct 2011 - R. Yantosca - Modified for central chemistry timestep
!  06 Feb 2018 - E. Lundgren - Updates for new timestep unit ([sec] not [min])
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES
!
      LOGICAL  :: NO_CENTRAL
      INTEGER  :: HH,     MM,   SS
      REAL(f8) :: GMTb,   TS

      !==============================================================
      ! Initialization
      !==============================================================

      ! Save the optional argument NO_CCTS in a local shadow variable
      IF ( PRESENT( NO_CCTS ) ) THEN
         NO_CENTRAL = NO_CCTS
      ELSE
         NO_CENTRAL = .FALSE.
      ENDIF

      ! Emissions timestep [hours]
      TS = DBLE( TS_EMIS ) / 3600e+0_f8
     
      !==============================================================
      ! FOR JANUARY 1st OF THE YEAR
      !==============================================================
      IF ( MONTH == 1 .and. DAY == 1 ) THEN

         IF ( NO_CENTRAL ) THEN

            ! Here, we are not using the central chemistry timestep.  
            ! Therefore, the first data read of the year should occur
            ! at 00:00 GMT on Jan 1st.
            IS_NEW_YEAR = ( NHMS == 000000 )

         ELSE

            ! Here, we are using the central chemistry timestep option.
            ! Therefore, the first data read of the year will occur not at 
            ! 00:00 GMT on the Jan 1st, but offset by a small amount (as 
            ! diagnosed by function ITS_TIME_FOR_EMIS).
            IS_NEW_YEAR = ( GMT < TS .and. ITS_TIME_FOR_EMIS() )

         ENDIF

      ELSE IF ( NYMD == NYMDb ) THEN

         !==============================================================
         ! FOR THE FIRST DAY OF THE SIMULATION
         ! (i.e. for simulations that start at other times of the year)
         !==============================================================
         IF ( NO_CENTRAL ) THEN

            ! Here, we are not using the central chemistry timestep.  
            ! Therefore, the first data read of this year should occur 
            ! at the start time of the simulation.
            IS_NEW_YEAR = ( NHMS == NHMSb )

         ELSE

            ! Split starting time into hour, minute, second
            CALL YMD_EXTRACT( NHMSb, HH, MM, SS )

            ! Compute GMT at the start of the simulation
            GMTb         = DBLE( HH ) + ( DBLE( MM ) / 60e+0_f8 )
            
            ! Here, we are using the central chemistry timestep option.
            ! Therefore, the first data read of the year will occur not 
            ! at  00:00 GMT on Jan 1st, but offset by a small amount (as  
            ! diagnosed by function ITS_TIME_FOR_EMIS).
            IS_NEW_YEAR  = ( GMT < GMTb+TS .and. ITS_TIME_FOR_EMIS() )

         ENDIF

      ELSE

         !==============================================================
         ! FOR ALL OTHER DAYS
         !==============================================================

         ! It isn't time for the first data read of the year; return FALSE
         IS_NEW_YEAR = .FALSE.
         
      ENDIF

      END FUNCTION ITS_A_NEW_YEAR
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Its_A_New_Month
!
! !DESCRIPTION: Function ITS\_A\_NEW\_MONTH returns TRUE if it's the first
!  timestep of the month when we have to read in monthly data.
!\\
!\\
! !INTERFACE:
!
      FUNCTION ITS_A_NEW_MONTH( NO_CCTS ) RESULT( IS_NEW_MONTH )
!
! !INPUT PARAMETERS
!
      LOGICAL, OPTIONAL :: NO_CCTS       ! =T reverts to previous behavior
                                         ! (i.e. w/o using central chem step)
!
! !RETURN VALUE:
!
      LOGICAL           :: IS_NEW_MONTH  ! =T if it's 1st data read of month
!
! !REMARKS:
!  ITS_A_NEW_MONTH assumes that we are using the central chemistry timestep
!  option (i.e. do chemistry & emissions & related processes at the midpoint
!  of each chemistry timestep).  To revert to the prior behavior, set the
!  optional flag NO_CCTS = .TRUE.
!                                                                             .
!  If we are using the central chemistry timestep option (which is now the
!  default behavior), then we must not read data at 00:00 GMT on the first day 
!  of the month, but at the center of the first chemistry timestep of the 
!  month.  This is because emissions and chemistry are done at the same time.  
!  The proper time of day for reading emissions is determined by function 
!  ITS_TIME_FOR_EMIS, also within time_mod.f.
!                                                                             .
!  Similarly, for simulations that start at an arbitrary midmonth date and 
!  time, we must not read data at the starting date and time of the simulation,
!  but at the midpoint of the first chemistry timestep of the simulation.
!                                                                             .
!  If we are not using the central chemistry timestep option (specified by
!  NO_CCTS=.TRUE.), then the first data read of the month occurs at 00:00 GMT 
!  on the first day of the month.  Similarly, for those simulations that start
!  at midmonth, the first data read will occur the starting date and time
!  of the simulation.
! 
! !REVISION HISTORY: 
!  01 Apr 2004 - R. Yantosca - Initial Version
!  15 Jan 2010 - R. Yantosca - Added ProTeX headers
!  12 Oct 2011 - R. Yantosca - Modified for central chemistry timestep option
!  27 Apr 2016 - E. Lundgren - Include minute condition for first day of month
!  06 Feb 2018 - E. Lundgren - Updates for new timestep unit ([sec] not [min])
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES
!
      LOGICAL  :: NO_CENTRAL
      INTEGER  :: HH,   MM,  SS
      REAL(f8) :: GMTb, TS

      !==============================================================
      ! Initialization
      !==============================================================

      ! Save 
      IF ( PRESENT( NO_CCTS ) ) THEN
         NO_CENTRAL = NO_CCTS
      ELSE
         NO_CENTRAL = .FALSE.
      ENDIF

      ! Emissions timestep [hours]
      TS = DBLE( TS_EMIS ) / 3600e+0_f8

      IF ( DAY == 1 .and. HOUR == 0 .and. MINUTE == 0 ) THEN

         !==============================================================
         ! FOR THE FIRST DAY OF THE MONTH
         !==============================================================
         IF ( NO_CENTRAL ) THEN

            ! Here, we are not using the central chemistry timestep.  
            ! Therefore, the first data read of the month should occur 
            ! at 00:00 GMT of the 1st day of the month.
            IS_NEW_MONTH = ( NHMS == 000000 ) 

         ELSE

            ! Here, we are using the central chemistry timestep option.
            ! Therefore, the first data read of the month will occur not at 
            ! 00:00 GMT on the 1st day of the month, but offset by a small
            ! amount (as diagnosed by function ITS_TIME_FOR_EMIS).
            IS_NEW_MONTH = ( GMT < TS .and. ITS_TIME_FOR_EMIS() )

         ENDIF
      
      ELSE IF ( NYMD == NYMDb ) THEN

         !==============================================================
         ! FOR THE FIRST DAY OF THE SIMULATION
         ! (i.e. for simulations that start at other times of the month)
         !==============================================================
         IF ( NO_CENTRAL ) THEN

            ! Here, we are not using the central chemistry timestep.  
            ! Therefore, the first data read of this month should occur 
            ! at the start time of the simulation.
            IS_NEW_MONTH = ( NHMS == NHMSb )

         ELSE

            ! Split starting time into hour, minute, second
            CALL YMD_EXTRACT( NHMSb, HH, MM, SS )

            ! Compute GMT at the start of the simulation
            GMTb         = DBLE( HH ) + ( DBLE( MM ) / 60e+0_f8 )
            
            ! Here, we are using the central chemistry timestep option.
            ! Therefore, the first data read of the month will occur not at 
            ! 00:00 GMT on the 1st day of the month, but offset by a small
            ! amount (as diagnosed by function ITS_TIME_FOR_EMIS).
            IS_NEW_MONTH = ( GMT < GMTb+TS .and. ITS_TIME_FOR_EMIS() )

         ENDIF

      ELSE

         !==============================================================
         ! FOR ALL OTHER DAYS
         !==============================================================

         ! It isn't time for the first data read of the month; return FALSE
         IS_NEW_MONTH = .FALSE.

      ENDIF
      
      END FUNCTION ITS_A_NEW_MONTH
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Its_MidMonth
!
! !DESCRIPTION: Function ITS\_MIDMONTH returns TRUE if it's the middle of a 
!  month.
!\\
!\\
! !INTERFACE:
!
      FUNCTION ITS_MIDMONTH() RESULT( IS_MIDMONTH )
!
! !RETURN VALUE:
!
      LOGICAL :: IS_MIDMONTH
!
! !REVISION HISTORY: 
!  10 Oct 2005 - S. Strode   - Initial Version
!  15 Jan 2010 - R. Yantosca - Added ProTeX headers
!  14 Oct 2011 - R. Yantosca - Modified for central chemistry timestep 
!EOP
!------------------------------------------------------------------------------
!BOC

      ! Test for the 16th of the month at 0 GMT
      IS_MIDMONTH = ( DAY == 16 .and. ITS_A_NEW_DAY() )

      END FUNCTION ITS_MIDMONTH
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Its_A_New_Day
!
! !DESCRIPTION: Function ITS\_A\_NEW\_DAY returns TRUE if it's the first
!  timestep of the day when we have to read in daily data.
!\\
!\\
! !INTERFACE:
!
      FUNCTION ITS_A_NEW_DAY( NO_CCTS ) RESULT( IS_NEW_DAY )
!
! !INPUT PARAMETERS
!
      LOGICAL, OPTIONAL :: NO_CCTS       ! =T reverts to previous behavior
                                         ! (i.e. w/o using central chem step)
!
! !RETURN VALUE:
!
      LOGICAL           :: IS_NEW_DAY    ! =T if it's 1st data read of day
!
! !REMARKS:
!  ITS_A_NEW_DAY assumes that we are using the central chemistry timestep
!  option (i.e. do chemistry & emissions & related processes at the midpoint
!  of each chemistry timestep).  To revert to the prior behavior, set the
!  optional flag NO_CCTS = .TRUE.
!                                                                             .
!  If we are using the central chemistry timestep option (which is now the
!  default behavior), then we must not read data at 00:00 GMT of each day,
!  but at the center of the first chemistry timestep of the day.  This is 
!  because emissions and chemistry are done at the same time.  The proper 
!  time of day for reading emissions is determined by function 
!  ITS_TIME_FOR_EMIS, also within time_mod.f.
!                                                                             .
!  Similarly, for simulations that start at an arbitrary midmonth date and 
!  time, we must not read data at the starting date and time of the simulation,
!  but at the midpoint of the first chemistry timestep of the simulation.
!                                                                             .
!  If we are not using the central chemistry timestep option (specified by
!  NO_CCTS=.TRUE.), then the first data read of the month occurs at 00:00 GMT 
!  each day.  Similarly, for those simulations that start at midmonth, the 
!  first data read will occur the starting date and time of the simulation.
!
! !REVISION HISTORY: 
!  01 Apr 2004 - R. Yantosca - Initial Version
!  15 Jan 2010 - R. Yantosca - Added ProTeX headers
!  14 Oct 2011 - R. Yantosca - Modified for central chemistry timestep
!  06 Feb 2018 - E. Lundgren - Updates for new timestep unit ([sec] not [min])
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES
!
      LOGICAL  :: NO_CENTRAL
      INTEGER  :: HH,     MM,   SS
      REAL(f8) :: GMTb,   TS

      !==============================================================
      ! Initialization
      !==============================================================

      ! Save optional argument NO_CCTS in a local shadow variable
      IF ( PRESENT( NO_CCTS ) ) THEN
         NO_CENTRAL = NO_CCTS
      ELSE
         NO_CENTRAL = .FALSE.
      ENDIF

      ! Emissions timestep [hours]
      TS = DBLE( TS_EMIS ) / 3600e+0_f8
     
      IF ( NYMD == NYMDb ) THEN

         !==============================================================
         ! FOR THE FIRST DAY OF THE SIMULATION
         ! (i.e. for simulations that start at any time of the year)
         !==============================================================
         IF ( NO_CENTRAL ) THEN

            ! Here, we are not using the central chemistry timestep option 
            ! Therefore, the first data read of this day should occur at
            ! the start time of the simulation.
            IS_NEW_DAY = ( NHMS == NHMSb )

         ELSE

            ! Split starting time into hour, minute, second
            CALL YMD_EXTRACT( NHMSb, HH, MM, SS )

            ! Compute GMT at the start of the simulation
            GMTb       = DBLE( HH ) + ( DBLE( MM ) / 60e+0_f8 )
            
            ! Here, we are using the central chemistry timestep option.
            ! Therefore, the first data read of the day will occur not at 
            ! 00:00 GMT, but offset by a small amount (as diagnosed by 
            ! function ITS_TIME_FOR_EMIS).
            IS_NEW_DAY = ( GMT < GMTb+TS .and. ITS_TIME_FOR_EMIS() )

         ENDIF

      ELSE

         !==============================================================
         ! FOR EACH NEW DAY
         !==============================================================
         IF ( NO_CENTRAL ) THEN

            ! Here, we are not using the central chemistry timestep.  
            ! Therefore, the first data read of this day should occur 
            ! at 00:00 GMT.
            IS_NEW_DAY = ( NHMS == 000000 )

         ELSE

            ! Here, we are using the central chemistry timestep option.
            ! Therefore, the first data read of the day will occur not at 
            ! 00:00 GMT, but offset by a small amount (as diagnosed by 
            ! function ITS_TIME_FOR_EMIS).
            IS_NEW_DAY = ( GMT < TS .and. ITS_TIME_FOR_EMIS() )

         ENDIF

      ENDIF

      END FUNCTION ITS_A_NEW_DAY
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Its_A_New_Hour
!
! !DESCRIPTION: Function ITS\_A\_NEW\_HOUR returns TRUE if it's the first
!  timestep of a new hour (it also returns TRUE on the first timestep of the
!  run).  This is useful for setting flags for reading in data. (bmy, 4/1/04)
!\\
!\\
! !INTERFACE:
!
      FUNCTION ITS_A_NEW_HOUR( ) RESULT( IS_NEW_HOUR )
!
! !RETURN VALUE:
!
      LOGICAL :: IS_NEW_HOUR
!
! !REVISION HISTORY: 
!  01 Apr 2004 - R. Yantosca - Initial Version
!  25 Feb 2014 - M. Sulprizio- Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      !=================================================================
      ! ITS_A_NEW_HOUR begins here!
      !=================================================================
      IF ( MOD( NHMS, 10000 ) == 0 ) THEN

         ! Test if it's 0 GMT
         IS_NEW_HOUR = .TRUE.

      ELSE IF ( NYMD == NYMDb .and. NHMS == NHMSb ) THEN

         ! Also return TRUE if it's the start of the run
         ! (since files will need to be read in from disk)
         IS_NEW_HOUR = .TRUE.

      ELSE

         ! Otherwise, it's not a new year
         IS_NEW_HOUR = .FALSE.
         
      ENDIF

      END FUNCTION ITS_A_NEW_HOUR
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Its_A_New_Season
!
! !DESCRIPTION: Function ITS\_A\_NEW\_SEASON returns TRUE if it's a new season 
!  or FALSE if it's not a new season.  Seasons are (1=DJF, 2=MAM, 3=JJA, 
!  4=SON).
!\\
!\\
! !INTERFACE:
!
      FUNCTION ITS_A_NEW_SEASON( ) RESULT( IS_NEW_SEASON )
!
! !RETURN VALUE:
!
      LOGICAL :: IS_NEW_SEASON  
!
! !REVISION HISTORY: 
!  20 Jul 2004 - R. Yantosca - Initial Version
!  15 Jan 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER, SAVE :: LAST_SEASON = -1
      
      IF ( NSEASON /= LAST_SEASON ) THEN
         IS_NEW_SEASON = .TRUE.
         LAST_SEASON   = NSEASON
      ELSE
         IS_NEW_SEASON = .FALSE.
      ENDIF

      END FUNCTION ITS_A_NEW_SEASON
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Print_Current_Time
!
! !DESCRIPTION: Subroutine PRINT\_CURRENT\_TIME prints the date, GMT time, and 
!  elapsed hours of a GEOS-Chem simulation. 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE PRINT_CURRENT_TIME
! 
! !REVISION HISTORY: 
!  21 Mar 2003 - R. Yantosca - Initial Version
!  15 Jan 2010 - R. Yantosca - Added ProTeX headers
!  10 Oct 2017 - R. Yantosca - Now print 6 decimal digits for elapsed time
!  10 Oct 2017 - R. Yantosca - Change "GMT" to "UTC", which is more standard
!  06 Feb 2018 - E. Lundgren - Updates for new timestep unit ([sec] not [min])
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL(f4) :: E_HOURS

      ! Hours since start of run
      E_HOURS = REAL( ELAPSED_SEC ) / 3600e+0_f4

      ! Write quantities
      WRITE( 6, 100 ) YEAR, MONTH, DAY, HOUR, MINUTE, E_HOURS

      ! Format string
 100  FORMAT( '---> DATE: ', i4.4, '/', i2.2, '/', i2.2, 
     &            '  UTC: ', i2.2, ':', i2.2, '  X-HRS: ', f13.6 )

      END SUBROUTINE PRINT_CURRENT_TIME
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Timestamp_String
!
! !DESCRIPTION: Function TIMESTAMP\_STRING returns a formatted string 
!  "YYYY/MM/DD hh:mm" for the a date and time specified by YYYYMMDD and hhmmss.
!  If YYYYMMDD and hhmmss are omitted, then TIMESTAMP\_STRING will create a 
!  formatted string for the current date and time.
!\\
!\\
! !INTERFACE:
!
      FUNCTION TIMESTAMP_STRING( YYYYMMDD, HHMMSS ) RESULT( TIME_STR )
!
! !INPUT PARAMETERS: 
!
      INTEGER, INTENT(IN), OPTIONAL :: YYYYMMDD   ! YYYY/MM/DD date
      INTEGER, INTENT(IN), OPTIONAL :: HHMMSS     ! hh:mm:ss time
!
! !RETURN VALUE:
!
      CHARACTER(LEN=16)             :: TIME_STR
!
! !REVISION HISTORY: 
!  21 Mar 2003 - R. Yantosca - Initial Version
!  (1 ) Now use ENCODE statement for PGI/F90 on Linux (bmy, 9/29/03)
!  (2 ) Now add optional arguments YYYYMMDD and HHMMSS (bmy, 10/27/03)
!  (3 ) Renamed LINUX to LINUX_PGI (bmy, 12/2/03)
!  15 Jan 2010 - R. Yantosca - Added ProTeX headers
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER :: THISYEAR, THISMONTH,  THISDAY
      INTEGER :: THISHOUR, THISMINUTE, THISSECOND

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

      END FUNCTION TIMESTAMP_STRING
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Ymd_Extract
!
! !DESCRIPTION: Subroutine YMD\_EXTRACT extracts the year, month, and date 
!  from an integer variable in YYYYMMDD format.  It can also extract the 
!  hours, minutes, and seconds from a variable in HHMMSS format.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE YMD_EXTRACT( NYMD, Y, M, D )
!
! !INPUT PARAMETERS: 
!
      INTEGER, INTENT(IN)  :: NYMD      ! YYYY/MM/DD format date
!
! !OUTPUT PARAMETERS:
!
      INTEGER, INTENT(OUT) :: Y, M, D   ! Separated YYYY, MM, DD values
! 
! !REVISION HISTORY: 
!  21 Nov 2001 - R. Yantosca - Initial Version
!  15 Jan 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL(f8) :: REM

      ! Extract YYYY from YYYYMMDD 
      Y = INT( DBLE( NYMD ) / 1e+4_f8 )

      ! Extract MM from YYYYMMDD
      REM = DBLE( NYMD ) - ( DBLE( Y ) * 1e+4_f8 )
      M   = INT( REM / 1e+2_f8 )

      ! Extract DD from YYYYMMDD
      REM = REM - ( DBLE( M ) * 1e+2_f8 )
      D   = INT( REM )

      ! Return to calling program
      END SUBROUTINE YMD_EXTRACT
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Expand_Date
!
! !DESCRIPTION: Subroutine EXPAND\_DATE replaces "YYYYMMDD" and "hhmmss" 
!  tokens within a filename string with the actual values.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE EXPAND_DATE( FILENAME, YYYYMMDD, HHMMSS )
!
! !USES:
!
      USE CHARPAK_MOD, ONLY : STRREPL
!
! !INPUT PARAMETERS: 
!
      INTEGER,          INTENT(IN)    :: YYYYMMDD   ! YYYY/MM/DD date
      INTEGER,          INTENT(IN)    :: HHMMSS     ! hh:mm:ss time
!
! !INPUT/OUTPUT PARAMETERS: 
!
      CHARACTER(LEN=*), INTENT(INOUT) :: FILENAME   ! Filename to modify
! 
! !REVISION HISTORY: 
!  27 Jun 2002 - R. Yantosca - Initial Version
!  (1 ) Bug fix for Linux: use ENCODE statement to convert number to string 
!        instead of F90 internal read. (bmy, 9/29/03)
!  (2 ) Now replace 2 and 4 digit year strings for all models (bmy, 10/23/03)
!  (3 ) Renamed LINUX to LINUX_PGI (bmy, 12/2/03)
!  (4 ) Now do not replace "ss" with seconds, as the smallest GEOS-Chem
!        timestep is in minutes. (bmy, 7/20/04)
!  15 Jan 2010 - R. Yantosca - Added ProTeX headers
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER           :: YYYY, YY, MM, DD, HH, II, SS
      CHARACTER(LEN=2)  :: MM_STR, DD_STR
      CHARACTER(LEN=2)  :: HH_STR, II_STR, SS_STR
      CHARACTER(LEN=2)  :: YY_STR
      CHARACTER(LEN=4)  :: YYYY_STR

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

      END SUBROUTINE EXPAND_DATE
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Valid_Date
!
! !DESCRIPTION: Function VALID\_DATE returns TRUE if the input date is
!  a valid calendar date, or FALSE otherwise.
!\\
!\\
! !INTERFACE:
!
      FUNCTION Valid_Date( YYYYMMDD ) RESULT( Is_Valid )
!
! !INPUT PARAMETERS: 
!
      INTEGER, INTENT(IN) :: YYYYMMDD   ! YYYY/MM/DD date
!
! !RETURN VALUE:
!
      LOGICAL             :: Is_Valid   ! =T if YYYYMMDD is a valid date
                                        ! =F otherwise
!
! !REVISION HISTORY: 
!  06 Jul 2018 - R. Yantosca - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL :: Is_Leap
      INTEGER :: YYYY, MM, DD, HH, II, SS, LastDay

      !=================================================================
      ! VALID_DATE begins here!
      !=================================================================

      ! Assume success until otherwise proven
      Is_Valid = .TRUE.
      
      ! Extract the date value into year, month, and day sections
      CALL YMD_EXTRACT( YYYYMMDD, YYYY, MM, DD )

      ! Exit if month is out of range
      IF ( MM < 1 .or. MM > 12 ) THEN
         Is_Valid = .FALSE.
         RETURN
      ENDIF

      ! Get the last day of the month
      SELECT CASE( MM )

         ! Check for leap year day if it's Feb
         CASE( 2 )
            IF ( ITS_A_LEAPYEAR( YYYY ) ) THEN
               LastDay = 29
            ELSE
               LastDay = 28
            ENDIF
    
         ! 30 days hath September, April, June, and November
         CASE( 4, 6, 9, 11 )
            LastDay = 30

         ! All the rest have 31
         CASE DEFAULT
            LastDay = 31
             
      END SELECT

      ! Exit if the day is out of range
      IF ( DD < 1 .or. DD > LastDay ) THEN
         Is_Valid = .FALSE.
         RETURN
      ENDIF

      END FUNCTION Valid_Date
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Valid_Time
!
! !DESCRIPTION: Function VALID\_TIME returns TRUE if the input date is
!  a valid clock time, or FALSE otherwise.
!\\
!\\
! !INTERFACE:
!
      FUNCTION Valid_Time( HHMMSS ) RESULT( Is_Valid )
!
! !INPUT PARAMETERS: 
!
      INTEGER, INTENT(IN) :: HHMMSS     ! HH:MM:SS time
!
! !RETURN VALUE:
!
      LOGICAL             :: Is_Valid   ! =T if HHMMSS is a valid time
                                        ! =F otherwise
!
! !REVISION HISTORY: 
!  06 Jul 2018 - R. Yantosca - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER :: HH, MM, SS

      !=================================================================
      ! VALID_TIME begins here!
      !=================================================================

      ! Assume success until otherwise proven
      Is_Valid = .TRUE.
      
      ! Extract the time value into hour, minute, and second sections
      CALL YMD_EXTRACT( HHMMSS, HH, MM, SS )

      ! Exit if hours are out of range
      IF ( HH < 0 .or. HH > 23 ) THEN
         Is_Valid = .FALSE.
         RETURN
      ENDIF

      ! Exit if minutes are out of range
      IF ( MM < 0 .or. MM > 59 ) THEN
         Is_Valid = .FALSE.
         RETURN
      ENDIF

      ! Exit if seconds are out of range
      IF ( SS < 0 .or. SS > 59 ) THEN
         Is_Valid = .FALSE.
         RETURN
      ENDIF

      END FUNCTION Valid_Time
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: System_Date_Time
!
! !DESCRIPTION: Subroutine SYSTEM\_DATE\_TIME returns the actual local date 
!  and time (as opposed to the model date and time).
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE SYSTEM_DATE_TIME( SYS_NYMD, SYS_NHMS )
!
! !OUTPUT PARAMETERS:
!
      INTEGER, INTENT(OUT) :: SYS_NYMD   ! System date in YYYY/MM/DD format
      INTEGER, INTENT(OUT) :: SYS_NHMS   ! System time in YYYY/MM/DD format
!
! !REMARKS:
!  Uses the F90 intrinsic function DATE_AND_TIME.
!
! !REVISION HISTORY: 
!  02 May 2005 - R. Yantosca - Initial Version
!  15 Jan 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! Arguments

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

      END SUBROUTINE SYSTEM_DATE_TIME
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: System_Timestamp
!
! !DESCRIPTION: Function SYSTEM\_TIMESTAMP returns a 16 character string with 
!  the system date and time in YYYY/MM/DD HH:MM format.
!\\
!\\
! !INTERFACE:
!
      FUNCTION SYSTEM_TIMESTAMP() RESULT( STAMP )
!
! !RETURN VALUE:
!
      CHARACTER(LEN=16) :: STAMP
! 
! !REVISION HISTORY: 
!  03 May 2005 - R. Yantosca - Initial version
!  15 Jan 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER           :: SYS_NYMD, SYS_NHMS

      !=================================================================
      ! SYSTEM_TIMESTAMP begins here!
      !=================================================================

      ! Get system date and time
      CALL SYSTEM_DATE_TIME( SYS_NYMD, SYS_NHMS )

      ! Create a string w/ system date & time
      STAMP = TIMESTAMP_STRING( SYS_NYMD, SYS_NHMS )

      END FUNCTION SYSTEM_TIMESTAMP
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: timestamp_diag
!
! !DESCRIPTION: Subroutine TIMESTAMP\_DIAG save timestamps to be used in 
!  filenames for diagnostics. We do not want the time when the diagnostic 
!  is saved but the time for previous dynamic time step because midnight is
!  considered as the beginning of next day (and not ending of previous day).
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE TIMESTAMP_DIAG
! 
! !REVISION HISTORY: 
!  12 Aug 2009 - C. Carouge  - Initial version
!  15 Jan 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC

      NYMD_DIAG = GET_NYMD()

      END SUBROUTINE TIMESTAMP_DIAG
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get_nymd_diag
!
! !DESCRIPTION: Function GET\_NYMD\_DIAG returns the previous NYMD value 
!  (YYYYMMDD) to the calling program.  Used for diagnostic filenames. 
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_NYMD_DIAG() RESULT( THISNYMD )
!
! !RETURN VALUE:
!
      INTEGER :: THISNYMD
! 
! !REVISION HISTORY:
!  12 Aug 2009 - C. Carouge  - Initial version
!  15 Jan 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC

      THISNYMD = NYMD_DIAG

      END FUNCTION GET_NYMD_DIAG
!EOC
#if defined( ESMF_ ) || defined( MODEL_ )
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Accept_External_Date_Time
!
! !DESCRIPTION: Subroutine ACCEPT\_EXTERNAL\_DATE\_TIME sets the date and
!  time variables in time\_mod.F with the values obtained from an external
!  GCM (such as NASA's GEOS-5 GCM).  The various date \& time values from 
!  the GCM are passed as arguments.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE Accept_External_Date_Time( 
     &   am_I_Root,    value_NYMDb,    value_NYMDe,   value_NYMD,
     &   value_NHMSb,  value_NHMSe,    value_NHMS,    value_YEAR,
     &   value_MONTH,  value_DAY,      value_DAYOFYR, value_HOUR,
     &   value_MINUTE, value_SECOND,   value_UTC,     value_HELAPSED,
     &   value_TS_CHEM, value_TS_CONV, value_TS_DYN,  value_TS_EMIS,
     &   RC                                                           )
!
! !USES:
!
      USE ErrCode_Mod
      USE JULDAY_MOD,       ONLY : JULDAY
      USE JULDAY_MOD,       ONLY : CALDATE
!
! !INPUT PARAMETERS: 
!
      LOGICAL,  INTENT(IN)  :: am_I_Root       ! Are we on the root CPU?
      INTEGER,  OPTIONAL    :: value_NYMDb     ! YYYY/MM/DD @ start of run
      INTEGER,  OPTIONAL    :: value_NYMDe     ! YYYY/MM/DD @ end of run
      INTEGER,  OPTIONAL    :: value_NYMD      ! YYYY/MM/DD @ current time
      INTEGER,  OPTIONAL    :: value_NHMSb     ! hh:mm:ss   @ start of run
      INTEGER,  OPTIONAL    :: value_NHMSe     ! hh:mm:ss   @ end of run
      INTEGER,  OPTIONAL    :: value_NHMS      ! hh:mm:ss   @ current time
      INTEGER,  OPTIONAL    :: value_YEAR      ! UTC year 
      INTEGER,  OPTIONAL    :: value_MONTH     ! UTC month
      INTEGER,  OPTIONAL    :: value_DAY       ! UTC day
      INTEGER,  OPTIONAL    :: value_DAYOFYR   ! UTC day of year
      INTEGER,  OPTIONAL    :: value_HOUR      ! UTC hour
      INTEGER,  OPTIONAL    :: value_MINUTE    ! UTC minute
      INTEGER,  OPTIONAL    :: value_SECOND    ! UTC second
      REAL(f4), OPTIONAL    :: value_UTC       ! UTC time [hrs]
      REAL(f4), OPTIONAL    :: value_HELAPSED  ! Elapsed hours
      INTEGER,  OPTIONAL    :: value_TS_CHEM   ! Chemistry  timestep [sec]
      INTEGER,  OPTIONAL    :: value_TS_CONV   ! Convection timestep [sec]
      INTEGER,  OPTIONAL    :: value_TS_DYN    ! Dynamic    timestep [sec]
      INTEGER,  OPTIONAL    :: value_TS_EMIS   ! Emissions  timestep [sec]
!
! !OUTPUT ARGUMENTS:
!
      INTEGER, INTENT(OUT) :: RC          ! Success or failure?
!
! !REMARKS:
!  The date and time values are obtained via the Extract_ subroutine in 
!  module file GEOSCHEMchem_GridCompMod.F90.      
! 
! !REVISION HISTORY: 
!  06 Dec 2012 - Initial version
!  11 Dec 2012 - R. Yantosca - Renamed to ACCEPT_EXTERNAL_DATE_TIME
!  18 Jun 2013 - R. Yantosca - Now compute day of week w/r/t GMT, which is
!                              the same modification made in SET_CURRENT_TIME
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL(f8) :: A, B, JD, THISDAY, TMP

      !=================================================================
      ! Set time/date values of time_mod.F from the ESMF inputs 
      !=================================================================
      IF ( PRESENT( value_NYMDb    ) ) NYMDb       = value_NYMdb
      IF ( PRESENT( value_NYMDe    ) ) NYMDe       = value_NYMDe
      IF ( PRESENT( value_NYMD     ) ) NYMD        = value_NYMD
      IF ( PRESENT( value_NHMSb    ) ) NHMSb       = value_NHMSb 
      IF ( PRESENT( value_NHMSe    ) ) NHMSe       = value_NHMSe
      IF ( PRESENT( value_NHMS     ) ) NHMS        = value_NHMS 
      IF ( PRESENT( value_YEAR     ) ) YEAR        = value_YEAR 
      IF ( PRESENT( value_MONTH    ) ) MONTH       = value_MONTH
      IF ( PRESENT( value_DAY      ) ) DAY         = value_DAY 
      IF ( PRESENT( value_DAYOFYR  ) ) DAY_OF_YEAR = value_DAYOFYR
      IF ( PRESENT( value_HOUR     ) ) HOUR        = value_HOUR
      IF ( PRESENT( value_MINUTE   ) ) MINUTE      = value_MINUTE
      IF ( PRESENT( value_SECOND   ) ) SECOND      = value_SECOND
      IF ( PRESENT( value_TS_CHEM  ) ) TS_CHEM     = value_TS_CHEM
      IF ( PRESENT( value_TS_CONV  ) ) TS_CONV     = value_TS_CONV
      IF ( PRESENT( value_TS_DYN   ) ) TS_DYN      = value_TS_DYN
      IF ( PRESENT( value_TS_EMIS  ) ) TS_EMIS     = value_TS_EMIS

      ! Special handling for GMT to avoid roundoff error
      IF ( PRESENT( value_UTC  ) ) THEN 
         TMP = value_UTC
         GMT = TMP
      ENDIF

      !=================================================================
      ! Compute various other derived time/date values
      !=================================================================
      
      ! Elapsed seconds since the start of the run
      IF ( PRESENT( value_HELAPSED ) ) THEN
         ELAPSED_SEC = value_HELAPSED * 3600
      ENDIF

      ! TAUb (hours since 0 UTC on 01 Jan 1985) @ start of simulation
      IF ( PRESENT( value_NYMDb ) .and. PRESENT( value_NHMSb ) ) THEN
         TMP  = ( GET_JD( NYMDb, NHMSb ) - JD85 ) * 24e+0_f8
         TAUb = DBLE( TMP )
      ENDIF

      ! TAUe (hours since 0 UTC on 01 Jan 1985) @ end of simulation
      IF ( PRESENT( value_NYMDe ) .and. PRESENT( value_NHMSe ) ) THEN
         TMP  = ( GET_JD( NYMDe, NHMSe ) - JD85 ) * 24e+0_f8
         TAUe = DBLE( TMP )
      ENDIF

      ! TAU value ( hours since 0 UTC on 01 Jan 1985)
      IF ( PRESENT( value_NYMD ) .and. PRESENT( value_NHMS ) ) THEN
         TMP  = ( GET_JD( NYMD, NHMS ) - JD85 ) * 24e+0_f8
         TAU  = DBLE( TMP )
      ENDIF

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

      ! Day of week w/r/t the GMT date  
      ! Use same algorithm as in routine SET_CURRENT_TIME
      THISDAY     = DAY + ( GMT / 24e+0_f8 )
      JD          = JULDAY( YEAR, MONTH, THISDAY )
      A           = ( JD + 1.5e+0_f8 ) / 7e+0_f8
      B           = ( A - INT( A ) ) * 7e+0_f8
      B           = INT( NINT( B*1e+5_f8 + SIGN( 5e+0_f8, B ) ) 
     &            / 10e+0_f8 ) / 1e+4_f8
      DAY_OF_WEEK = INT( B )

      ! Return successfully
      RC = GC_SUCCESS

      END SUBROUTINE Accept_External_Date_Time
!EOC
#endif
      END MODULE TIME_MOD




