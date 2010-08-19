!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: time_mod
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
      IMPLICIT NONE
      PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC  :: SET_CURRENT_TIME
      PUBLIC  :: SET_BEGIN_TIME
      PUBLIC  :: SET_END_TIME
      PUBLIC  :: SET_NDIAGTIME
      PUBLIC  :: SET_DIAGb
      PUBLIC  :: SET_DIAGe
      PUBLIC  :: SET_TIMESTEPS
      PUBLIC  :: SET_CT_CHEM
      PUBLIC  :: SET_CT_CONV
      PUBLIC  :: SET_CT_DYN
      PUBLIC  :: SET_CT_EMIS
      PUBLIC  :: SET_CT_DIAG
      PUBLIC  :: SET_CT_A1
      PUBLIC  :: SET_CT_A3
      PUBLIC  :: SET_CT_A6
      PUBLIC  :: SET_CT_I6
      PUBLIC  :: SET_CT_XTRA
      PUBLIC  :: SET_ELAPSED_MIN
      PUBLIC  :: GET_JD
      PUBLIC  :: GET_ELAPSED_MIN
      PUBLIC  :: GET_ELAPSED_SEC
      PUBLIC  :: GET_NYMDb
      PUBLIC  :: GET_NHMSb
      PUBLIC  :: GET_NYMDe
      PUBLIC  :: GET_NHMSe
      PUBLIC  :: GET_NYMD
      PUBLIC  :: GET_NHMS
      PUBLIC  :: GET_NDIAGTIME
      PUBLIC  :: GET_TIME_AHEAD
      PUBLIC  :: GET_MONTH
      PUBLIC  :: GET_DAY
      PUBLIC  :: GET_YEAR
      PUBLIC  :: GET_HOUR
      PUBLIC  :: GET_MINUTE
      PUBLIC  :: GET_SECOND
      PUBLIC  :: GET_DAY_OF_YEAR
      PUBLIC  :: GET_DAY_OF_WEEK
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
      PUBLIC  :: GET_TS_SUN_2
      PUBLIC  :: GET_CT_CHEM
      PUBLIC  :: GET_CT_CONV
      PUBLIC  :: GET_CT_DYN
      PUBLIC  :: GET_CT_EMIS
      PUBLIC  :: GET_CT_A1
      PUBLIC  :: GET_CT_A3
      PUBLIC  :: GET_CT_A6
      PUBLIC  :: GET_CT_I6
      PUBLIC  :: GET_CT_XTRA
      PUBLIC  :: GET_CT_DIAG
      PUBLIC  :: GET_A1_TIME
      PUBLIC  :: GET_A3_TIME
      PUBLIC  :: GET_A6_TIME
      PUBLIC  :: GET_I6_TIME
      PUBLIC  :: GET_FIRST_A1_TIME
      PUBLIC  :: GET_FIRST_A3_TIME
      PUBLIC  :: GET_FIRST_A6_TIME
      PUBLIC  :: ITS_TIME_FOR_CHEM
      PUBLIC  :: ITS_TIME_FOR_CONV
      PUBLIC  :: ITS_TIME_FOR_DYN
      PUBLIC  :: ITS_TIME_FOR_EMIS
      PUBLIC  :: ITS_TIME_FOR_UNIT
      PUBLIC  :: ITS_TIME_FOR_DIAG
      PUBLIC  :: ITS_TIME_FOR_A3
      PUBLIC  :: ITS_TIME_FOR_A6
      PUBLIC  :: ITS_TIME_FOR_I6
      PUBLIC  :: ITS_TIME_FOR_UNZIP     
      PUBLIC  :: ITS_TIME_FOR_DEL
      PUBLIC  :: ITS_TIME_FOR_EXIT
      PUBLIC  :: ITS_TIME_FOR_BPCH
      PUBLIC  :: ITS_A_LEAPYEAR
      PUBLIC  :: ITS_A_NEW_YEAR
      PUBLIC  :: ITS_A_NEW_MONTH
      PUBLIC  :: ITS_MIDMONTH
      PUBLIC  :: ITS_A_NEW_DAY
      PUBLIC  :: ITS_A_NEW_SEASON
      PUBLIC  :: PRINT_CURRENT_TIME
      PUBLIC  :: TIMESTAMP_STRING
      PUBLIC  :: YMD_EXTRACT
      PUBLIC  :: EXPAND_DATE
      PUBLIC  :: SYSTEM_DATE_TIME
      PUBLIC  :: SYSTEM_TIMESTAMP
      PUBLIC  :: TIMESTAMP_DIAG
      PUBLIC  :: GET_NYMD_DIAG
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
!EOP
!------------------------------------------------------------------------------
!BOC
      ! Date and time variables
      INTEGER           :: NYMDb,      NHMSb,       NYMDe   
      INTEGER           :: NHMSe,      NYMD,        NHMS
      INTEGER           :: NYMD_DIAG
      INTEGER           :: MONTH,      DAY,         YEAR
      INTEGER           :: HOUR,       MINUTE,      SECOND
      INTEGER           :: NSEASON,    DAY_OF_YEAR, ELAPSED_MIN
      INTEGER           :: NDIAGTIME
      REAL*8            :: TAU,        TAUb,        TAUe  
      REAL*8            :: GMT,        DIAGb,       DIAGe

      ! Timesteps
      INTEGER           :: TS_CHEM,    TS_CONV,     TS_DIAG
      INTEGER           :: TS_DYN,     TS_EMIS,     TS_UNIT
      INTEGER           :: TS_SUN_2

      ! Timestep counters
      INTEGER           :: CT_CHEM,    CT_CONV,     CT_DYN    
      INTEGER           :: CT_EMIS,    CT_A3,       CT_A6
      INTEGER           :: CT_I6,      CT_XTRA,     CT_DIAG
      INTEGER           :: CT_A1

      ! Astronomical Julian Date at 0 GMT, 1 Jan 1985
      REAL*8, PARAMETER :: JD85 = 2446066.5d0

      CONTAINS
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: set_current_time
!
! !DESCRIPTION: Subroutine SET\_CURRENT\_TIME takes in the elapsed time in 
!  minutes since the start of a GEOS-Chem simulation and sets the GEOS-Chem 
!  time variables accordingly.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE SET_CURRENT_TIME
!
! !USES:
!
      USE JULDAY_MOD, ONLY : JULDAY, CALDATE

#     include "define.h"
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
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL :: IS_LEAPYEAR
      REAL*4  :: TMP
      REAL*8  :: JD0, JD1, JD_JAN_1
      
      !=================================================================
      ! SET_CURRENT_TIME begins here!
      !=================================================================

      ! JD0: Astronomical Julian Date at start of GEOS-Chem run
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
         IF (  ( JD1 - JD_JAN_1 >= 59d0 )  .and.
     &         ( JD0 - JD_JAN_1 <= 59d0 ) ) THEN
            JD1 = JD1 + 1d0
         ENDIF

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
      ! Compute other GEOS-Chem timing variables
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
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: set_begin_time
!
! !DESCRIPTION: Subroutine SET\_BEGIN\_TIME initializes NYMDb, NHMSb, and TAUb,
!  which are the YYYYMMDD, HHMMSS, and hours since 1/1/1985 corresponding to 
!  the beginning date and time of a GEOS-Chem run.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE SET_BEGIN_TIME( THISNYMDb, THISNHMSb )
!
! !USES:
!
      USE ERROR_MOD,  ONLY : ERROR_STOP
!
! !INPUT PARAMETERS: 
!
      INTEGER, INTENT(IN) :: THISNYMDb   ! YYYYMMDD @ start of G-C simulation
      INTEGER, INTENT(IN) :: THISNHMSb   ! HHMMSS   @ start of G-C simulation
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
      REAL*4 :: TMP

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
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: set_end_time
!
! !DESCRIPTION: Subroutine SET\_END\_TIME initializes NYMDe, NHMSe, and TAUe, 
!  which are the YYYYMMDD, HHMMSS, and hours since 1/1/1985 corresponding to 
!  the ending date and time of a GEOS-Chem run.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE SET_END_TIME( THISNYMDe, THISNHMSe )
!
! !USES:
!
      USE ERROR_MOD,  ONLY : ERROR_STOP
!
! !INPUT PARAMETERS: 
!
      INTEGER, INTENT(IN) :: THISNYMDe   ! YYYYMMDD @ end of G-C simulation
      INTEGER, INTENT(IN) :: THISNHMSe   ! HHMMSS   @ end of G-C simulation
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
      REAL*4 :: TMP

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
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: set_ndiagtime
!
! !DESCRIPTION:Subroutine SET\_NDIAGTIME initializes NDIAGTIME, the time of 
!  day at which the binary punch file will be written out to disk.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE SET_NDIAGTIME( THIS_NDIAGTIME )
!
! !INPUT PARAMETERS: 
!
      INTEGER, INTENT(IN) :: THIS_NDIAGTIME  ! Initial NDIAGTIMEe [hrs]
! 
! !REVISION HISTORY: 
!  20 Jul 2004 - R. Yantosca - Initial Version
!  15 Jan 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
      NDIAGTIME = THIS_NDIAGTIME

      END SUBROUTINE SET_NDIAGTIME
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: set_diagb
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
      REAL*8, INTENT(IN) :: THISDIAGb  ! Initial DIAGb value [hrs from 1/1/85]
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
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: set_diage
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
      REAL*8, INTENT(IN) :: THISDIAGe  ! Initial DIAGe value [hrs from 1/1/85]
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
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: set_timesteps
!
! !DESCRIPTION: Subroutine SET\_TIMESTEPS initializes the timesteps for 
!  dynamics, convection, chemistry, emissions, and diagnostics.  
!  Counters are also zeroed. 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE SET_TIMESTEPS( CHEMISTRY, CONVECTION, DYNAMICS,  
     &                          EMISSION,  UNIT_CONV,  DIAGNOS,
     &                          SUNCOS )
!
! !INPUT PARAMETERS:
!
      INTEGER, INTENT(IN) :: CHEMISTRY    ! Chemistry  timestep [min]
      INTEGER, INTENT(IN) :: CONVECTION   ! Convection timestep [min]
      INTEGER, INTENT(IN) :: DYNAMICS     ! Dynamic    timestep [min]
      INTEGER, INTENT(IN) :: EMISSION     ! Emission   timestep [min]
      INTEGER, INTENT(IN) :: UNIT_CONV    ! Unit conve timestep [min]
      INTEGER, INTENT(IN) :: DIAGNOS      ! Diagnostic timestep [min]
      INTEGER, INTENT(IN) :: SUNCOS       ! 1/2 of timestep for SUNCOS [min]
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
      TS_SUN_2 = SUNCOS

      ! Zero timestep counters
      CT_CHEM = 0
      CT_CONV = 0
      CT_DYN  = 0
      CT_EMIS = 0
      CT_A3   = 0
      CT_A6   = 0
      CT_I6   = 0
      CT_XTRA = 0
      CT_DIAG = 0

      ! Echo to stdout
      WRITE( 6, '(/,a)' ) 'SET_TIMESTEPS: setting GEOS-Chem timesteps!'
      WRITE( 6, '(  a)' ) '-------------------------------------------'
      WRITE( 6, '(''Chemistry  Timestep [min] : '', i4 )' ) TS_CHEM
      WRITE( 6, '(''Convection Timestep [min] : '', i4 )' ) TS_CONV
      WRITE( 6, '(''Dynamics   Timestep [min] : '', i4 )' ) TS_DYN
      WRITE( 6, '(''Emission   Timestep [min] : '', i4 )' ) TS_EMIS
      WRITE( 6, '(''Unit Conv  Timestep [min] : '', i4 )' ) TS_UNIT
      WRITE( 6, '(''Diagnostic Timestep [min] : '', i4 )' ) TS_DIAG
      WRITE( 6, '(''Offset for SUNCOS   [min] : '', i4 )' ) TS_SUN_2

      ! Return to calling program
      END SUBROUTINE SET_TIMESTEPS
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: set_ct_chem
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
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: set_ct_conv
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
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: set_ct_dyn
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
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: set_ct_emis
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
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: set_ct_diag
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
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: set_ct_a1
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
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: set_ct_a3
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
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: set_ct_a6
!
! !DESCRIPTION: Subroutine SET\_CT\_A6 increments CT\_A6, the counter of the 
!  number of times we have read in A-6 fields.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE SET_CT_A6( INCREMENT, RESET )
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
         CT_A6 = CT_A6 + 1
      ELSE IF ( PRESENT( RESET ) ) THEN
         CT_A6 = 0
      ENDIF

      END SUBROUTINE SET_CT_A6
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: set_ct_i6
!
! !DESCRIPTION: Subroutine SET\_CT\_I6 increments CT\_I6, the counter of the 
!  number of times we have read in I-6 fields.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE SET_CT_I6( INCREMENT, RESET )
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
         CT_I6 = CT_I6 + 1
      ELSE IF ( PRESENT( RESET ) ) THEN
         CT_I6 = 0
      ENDIF

      END SUBROUTINE SET_CT_I6
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: set_ct_xtra
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
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: set_elapsed_min
!
! !DESCRIPTION: Subroutine SET\_ELAPSED\_MIN increments the number of elapsed 
!  minutes by the dynamic timestep TS\_DYN.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE SET_ELAPSED_MIN
! 
! !REVISION HISTORY: 
!  05 Feb 2003 - R. Yantosca - Initial Version
!  15 Jan 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
      ELAPSED_MIN = ELAPSED_MIN + TS_DYN

      END SUBROUTINE SET_ELAPSED_MIN
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_jd
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
      REAL*8               :: THISJD     ! Output value
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
      INTEGER :: Y, M, D, H, MI, S
      REAL*8  :: DAY

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
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_elapsed_min
!
! !DESCRIPTION: Function GET\_ELAPSED\_MIN returns the elapsed minutes since 
!  the start of a GEOS-chem run.
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_ELAPSED_MIN() RESULT( THIS_ELAPSED_MIN )
!
! !RETURN VALUE:
!
      INTEGER :: THIS_ELAPSED_MIN
! 
! !REVISION HISTORY: 
!  05 Feb 2003 - R. Yantosca - Initial Version
!  15 Jan 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
      THIS_ELAPSED_MIN = ELAPSED_MIN

      END FUNCTION GET_ELAPSED_MIN
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_elapsed_sec
!
! !DESCRIPTION: Function GET\_ELAPSED\_SEC returns the elapsed minutes since 
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
      THIS_ELAPSED_SEC = ELAPSED_MIN * 60

      END FUNCTION GET_ELAPSED_SEC
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_nymdb
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
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_nhmsb
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
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_nymde
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
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_nhmse
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
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_nymd
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
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_nhms
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
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_ndiagtime
!
! !DESCRIPTION: Subroutine GET\_NDIAGTIME returns to the calling program 
!  NDIAGTIME, the time of day at which the binary punch file will be written 
!  out to disk. 
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_NDIAGTIME() RESULT( THIS_NDIAGTIME )
!
! !RETURN VALUE:
!
      INTEGER :: THIS_NDIAGTIME
! 
! !REVISION HISTORY: 
!  20 Jul 2004 - R. Yantosca - Initial Version
!  15 Jan 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
      THIS_NDIAGTIME = NDIAGTIME

      END FUNCTION GET_NDIAGTIME
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_time_ahead
!
! !DESCRIPTION: Function GET\_3h\_AHEAD returns to the calling program a 
!  2-element vector containing the YYYYMMDD and HHMMSS values at the current 
!  time plus N\_MINS minutes.
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_TIME_AHEAD( N_MINS ) RESULT( DATE )
!
! !USES:
!
      USE JULDAY_MOD, ONLY : CALDATE

#     include "define.h"   ! C-preprocessor flags
!
! !INPUT PARAMETERS:
!
      INTEGER, INTENT(IN) :: N_MINS   ! Minutes ahead to compute date & time
!
! !RETURN VALUE:
!
      INTEGER             :: DATE(2)  ! Date & time output
! 
! !REVISION HISTORY: 
!  21 Mar 2003 - R. Yantosca - Initial Version
!  (1 ) Bug fix for GCAP leap year case (phs, bmy, 12/8/06)
!  15 Jan 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
! 
      INTEGER :: THISYEAR, THISMONTH, THISDAY
      REAL*8  :: JD

      !=================================================================
      ! GET_TIME_AHEAD begins here!
      !=================================================================

      ! Astronomical Julian Date at current time + N_MINS
      JD = GET_JD( NYMD, NHMS ) + ( N_MINS / 1440d0 )

      ! Call CALDATE to compute the current YYYYMMDD and HHMMSS
      CALL CALDATE( JD, DATE(1), DATE(2) )

#if   defined( GCAP )

      !-------------------------------
      ! GCAP met fields: no leapyears
      !-------------------------------

      ! Extract current year, month, day from DATE(1)
      CALL YMD_EXTRACT( DATE(1), THISYEAR, THISMONTH, THISDAY )

      ! Special handling for leap years
      IF ( ITS_A_LEAPYEAR( THISYEAR, FORCE=.TRUE. )  .AND.
     &     THISMONTH == 2                            .AND.
     &     THISDAY   == 29 ) THEN 
           DATE(1) = ( THISYEAR * 10000 ) + 0301
        ENDIF

#endif

      ! Return to calling program
      END FUNCTION GET_TIME_AHEAD
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_month
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
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_day
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
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_year
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
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_hour
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
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_minute
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
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_second
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
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_day_of_year
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
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_day_of_week
!
! !DESCRIPTION: Function GET\_DAY\_OF\_WEEK returns the day of the week as a 
!  number: Sun=0, Mon=1, Tue=2, Wed=3, Thu=4, Fri=5, Sat=6.
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
! !REMARKS:
!  Reference:
!  ----------
!  "Practical Astronomy with Your Calculator", 3rd Ed.  Peter Duffett-Smith,
!    Cambridge UP, 1992, p9.
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
      
      END FUNCTION GET_DAY_OF_WEEK
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_gmt
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
      REAL*8 :: THISGMT   ! Greenwich mean time [hrs]
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
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_tau
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
      REAL*8 :: THISTAU  ! TAUb [hrs since 1/1/1985]
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
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_taub
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
      REAL*8 :: THISTAUb  ! TAUb [hrs since 1/1/1985]
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
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_taue
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
      REAL*8 :: THISTAUe  ! TAUe [hrs since 1/1/1985]
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
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_diagb
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
      INTEGER :: THISDIAGb   ! DIAGb [hrs sincd 1/1/1985]
! 
! !REVISION HISTORY: 
!  05 Feb 2003 - R. Yantosca - Initial Version
!  15 Jan 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
      THISDIAGb = DIAGb

      END FUNCTION GET_DIAGb
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_diage
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
      INTEGER :: THISDIAGe   ! DIAGe [hrs sincd 1/1/1985]
! 
! !REVISION HISTORY: 
!  05 Feb 2003 - R. Yantosca - Initial Version
!  15 Jan 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
      THISDIAGe = DIAGe

      END FUNCTION GET_DIAGe
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_localtime
!
! !DESCRIPTION: Function GET\_LOCALTIME returns the local time of a grid
!  box to the calling program. (bmy, 2/5/03)
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_LOCALTIME( I, OFFSET ) RESULT( THISLOCALTIME )
!
! !USES:
!
      USE GRID_MOD, ONLY : GET_XMID
!
! !INPUT PARAMETERS:
!
      INTEGER, INTENT(IN)           :: I        ! Longitude index
      REAL*8,  INTENT(IN), OPTIONAL :: OFFSET   ! Offset to apply to GMT [hrs]
!
! !RETURN VALUE:
!
      REAL*8                        :: THISLOCALTIME  ! Local time [hrs]
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
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL*8 :: OFF 

      ! Save the value of OFFSET in a local variable
      IF ( PRESENT( OFFSET ) ) THEN 
         OFF = OFFSET
      ELSE
         OFF = 0d0
      ENDIF

      ! Local time  =   GMT time [hrs]    +   longitude     / 15
      THISLOCALTIME = ( GET_GMT() + OFF ) + ( GET_XMID( I ) / 15d0 ) 

      ! Make sure that THISLOCALTIME is in the range 0-24 hours
      IF ( THISLOCALTIME > 24 ) THISLOCALTIME = THISLOCALTIME - 24d0 
      IF ( THISLOCALTIME < 0  ) THISLOCALTIME = THISLOCALTIME + 24d0

      ! Return to calling program
      END FUNCTION GET_LOCALTIME
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_season
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
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_ts_chem
!
! !DESCRIPTION: Function GET\_TS\_CHEM returns the chemistry timestep in 
!  minutes.
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_TS_CHEM() RESULT( THIS_TS_CHEM )
!
! !RETURN VALUE:
! 
      INTEGER :: THIS_TS_CHEM   ! ! Chemistry timestep [min]
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
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_ts_conv
!
! !DESCRIPTION: Function GET\_TS\_CONV returns the convection timestep in 
!  minutes.
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_TS_CONV() RESULT( THIS_TS_CONV )
!
! !RETURN VALUE:
!
      INTEGER :: THIS_TS_CONV   ! Convective timestep [min]
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
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_ts_diag
!
! !DESCRIPTION: Function GET\_TS\_DIAG returns the diagnostic timestep in 
!  minutes.
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_TS_DIAG() RESULT( THIS_TS_DIAG )
!
! !RETURN VALUE:
!
      INTEGER :: THIS_TS_DIAG   ! Diagnostic timestep [min]
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
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_ts_dyn
!
! !DESCRIPTION: Function GET\_TS\_DIAG returns the diagnostic timestep in 
!  minutes.
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_TS_DYN() RESULT( THIS_TS_DYN )
!
! !RETURN VALUE:
!
      INTEGER :: THIS_TS_DYN    ! Dynamic timestep [min]
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
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_ts_emis
!
! !DESCRIPTION: Function GET\_TS\_EMIS returns the emission timestep in 
!  minutes.
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_TS_EMIS() RESULT( THIS_TS_EMIS )
!
! !RETURN VALUE:
!
      INTEGER :: THIS_TS_EMIS   ! Emissions timestep [min]
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
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_ts_unit
!
! !DESCRIPTION: Function GET\_TS\_UNIT returns the unit-conversion timestep
!  in minutes.
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_TS_UNIT() RESULT( THIS_TS_UNIT )
!
! !RETURN VALUE:
!
      INTEGER :: THIS_TS_UNIT   ! Unit conversion timestep [min]
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
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_ts_sun_2
!
! !DESCRIPTION: Function GET\_TS\_SUN\_2 returns TS_SUN_2, which is 1/2 of
!  the interval at which we are computing the cosine of the solar zenith
!  angle, aka SUNCOS.  This is required to move the time at which we compute
!  SUNCOS to the middle of the chemistry timestep interval.
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_TS_SUN_2() RESULT( THIS_TS_SUN_2 )
!
! !RETURN VALUE:
!
      INTEGER :: THIS_TS_SUN_2  ! 1/2 of SUNCOS interval [min]
! 
! !REVISION HISTORY: 
!  27 Apr 2010 - R. Yantosca - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
      THIS_TS_SUN_2 = TS_SUN_2

      END FUNCTION GET_TS_SUN_2
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_ct_chem
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
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_ct_conv
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
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_ct_dyn
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
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_ct_emis
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
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_ct_a1
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
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_ct_a3
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
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_ct_a6
!
! !DESCRIPTION: Function GET\_CT\_A6 returns the A-6 fields timestep counter 
!  to the calling program.
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_CT_A6() RESULT( THIS_CT_A6 )
!
! !RETURN VALUE:
!
      INTEGER :: THIS_CT_A6   ! # of A-6 timesteps
! 
! !REVISION HISTORY: 
!  21 Mar 2003 - R. Yantosca - Initial Version
!  15 Jan 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
      THIS_CT_A6 = CT_A6

      END FUNCTION GET_CT_A6
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_ct_i6
!
! !DESCRIPTION: Function GET\_CT\_I6 returns the I-6 fields timestep counter
!  to the calling program
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_CT_I6() RESULT( THIS_CT_I6 )
!
! !RETURN VALUE:
!
      INTEGER :: THIS_CT_I6   ! # of I-6 timesteps
! 
! !REVISION HISTORY: 
!  21 Mar 2003 - R. Yantosca - Initial Version
!  15 Jan 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
      THIS_CT_I6 = CT_I6

      END FUNCTION GET_CT_I6
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_ct_xtra
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
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_ct_diag
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
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_a1_time
!
! !DESCRIPTION: Function GET\_A1\_TIME returns the correct YYYYMMDD and HHMMSS 
!  values that are needed to read in the next average 1-hour (A-1) fields. 
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_A1_TIME() RESULT( DATE )
!
! !USES:
!
#     include "define.h"
!
! !RETURN VALUE:
!
      INTEGER :: DATE(2)   ! YYYYMMDD and HHMMSS values
! 
! !REVISION HISTORY: 
!  19 Aug 2010 - R. Yantosca - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC

#if   defined( MERRA )

      ! MERRA met fields are 1-hour time-averages, timestamped at the
      ! center of the averaging periods (00:30, 01:30, 02:30 ... 23:30)
      DATE = GET_TIME_AHEAD( 30 )

#else
      
      ! Otherwise return the current time
      DATE = GET_TIME_AHEAD( 0 )

#endif

      END FUNCTION GET_A1_TIME
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_a3_time
!
! !DESCRIPTION: Function GET\_A3\_TIME returns the correct YYYYMMDD and HHMMSS 
!  values that are needed to read in the next average 3-hour (A-3) fields. 
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_A3_TIME() RESULT( DATE )
!
! !USES:
!
#     include "define.h"
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
!EOP
!------------------------------------------------------------------------------
!BOC

#if   defined( GEOS_3 )

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

      END FUNCTION GET_A3_TIME
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_a6_time
!
! !DESCRIPTION: Function GET\_A6\_TIME returns the correct YYYYMMDD and HHMMSS 
!  values that are needed to read in the next average 6-hour (A-6) fields. 
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_A6_TIME() RESULT( DATE )
!
! !RETURN VALUE:
!
      INTEGER :: DATE(2)   ! YYYYMMDD and HHMMSS time
! 
! !REVISION HISTORY: 
!  21 Mar 2003 - R. Yantosca - Initial Version
!  15 Jan 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC

      ! Return the time 3h (180m) from now, since there is a 3-hour 
      ! offset between the actual time when the A-6 fields are read
      ! and the time that the A-6 fields are stamped with.
      DATE = GET_TIME_AHEAD( 180 )      

      END FUNCTION GET_A6_TIME
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_i6_time
!
! !DESCRIPTION: Function GET\_I6\_TIME returns the correct YYYYMMDD and 
!  HHMMSS values that are needed to read in the next instantaneous 6-hour 
!  (I-6) fields. 
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_I6_TIME() RESULT( DATE )
!
! !RETURN VALUE:
!
      INTEGER :: DATE(2)   ! YYYYMMDD and HHMMSS values
!  
! !REVISION HISTORY: 
!  21 Mar 2003 - R. Yantosca - Initial Version
!  (1 ) Bug fix for GCAP: skip over Feb 29th (no leapyears). (bmy, 4/24/06)
!  15 Jan 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC

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

      END FUNCTION GET_I6_TIME
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_first_a1_time
!
! !DESCRIPTION: Function GET\_FIRST\_A1\_TIME returns the correct YYYYMMDD
!  and HHMMSS values the first time that A-3 fields are read in from disk. 
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_FIRST_A1_TIME() RESULT( DATE )
!
! !USES:
!
#     include "define.h"
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
!EOP
!------------------------------------------------------------------------------
!BOC
!
      ! Return the first A-1 time
      DATE = GET_A1_TIME()
  
      END FUNCTION GET_FIRST_A1_TIME
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_first_a3_time
!
! !DESCRIPTION: Function GET\_FIRST\_A3\_TIME returns the correct YYYYMMDD
!  and HHMMSS values the first time that A-3 fields are read in from disk. 
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_FIRST_A3_TIME() RESULT( DATE )
!
! !USES:
!
#     include "define.h"
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
!EOP
!------------------------------------------------------------------------------
!BOC
!
#if   defined( GEOS_3 )

      ! For GEOS-1, GEOS-STRAT, GEOS-3: Return the current date/time
      DATE = (/ NYMD, NHMS /)

#else
  
      ! For GEOS-4, GEOS-5, and GCAP: Call GET_A3_TIME to return 
      ! the date/time under which the A-3 fields are timestamped
      DATE = GET_A3_TIME()
    
#endif

      END FUNCTION GET_FIRST_A3_TIME
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_first_a6_time
!
! !DESCRIPTION: Function GET\_FIRST\_A6\_TIME returns the correct YYYYMMDD and
!  HHMMSS values the first time that A-6 fields are read in from disk.
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_FIRST_A6_TIME() RESULT( DATE )
!
! !USES:
! 
#     include "define.h"
!
! !RETURN VALUE:
!
      INTEGER :: DATE(2)    ! YYYYMMDD, HHMMSS values
! 
! !REVISION HISTORY: 
!  26 Jun 2003 - R. Yantosca - Initial Version
!  (1 ) Now modified for GEOS-4 "a_llk_03" and "a_llk_04" fields (bmy, 3/22/04)
!  (2 ) Modified for GCAP and GEOS-5 met fields (swu, bmy, 5/24/05)
!  15 Jan 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC

#if   defined( GCAP )

      ! For GCAP data: Call GET_A6_TIME to return date/time
      ! under which the A-6 fields are timestamped
      DATE = GET_A6_TIME()      

#else

      ! For GEOS data: Return the current date/time
      DATE = (/ NYMD, NHMS /)

#endif

      END FUNCTION GET_FIRST_A6_TIME
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: its_time_for_chem
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
!EOP
!------------------------------------------------------------------------------
!BOC
      ! Is it time for chemistry?
      FLAG = ( MOD( ELAPSED_MIN, TS_CHEM ) == 0 )

      END FUNCTION ITS_TIME_FOR_CHEM
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: its_time_for_conv
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
!EOP
!------------------------------------------------------------------------------
!BOC
      ! Is it time for convection?
      FLAG = ( MOD( ELAPSED_MIN, TS_CONV ) == 0 )

      END FUNCTION ITS_TIME_FOR_CONV
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: its_time_for_dyn
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
!EOP
!------------------------------------------------------------------------------
!BOC
      ! Is it time for dynamics?
      FLAG = ( MOD( ELAPSED_MIN, TS_DYN ) == 0 )

      END FUNCTION ITS_TIME_FOR_DYN
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: its_time_for_emis
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
!EOP
!------------------------------------------------------------------------------
!BOC
      ! Is it time for emissions?
      FLAG = ( MOD( ELAPSED_MIN, TS_EMIS ) == 0 )

      END FUNCTION ITS_TIME_FOR_EMIS
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: its_time_for_unit
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
!EOP
!------------------------------------------------------------------------------
!BOC

      ! Is it time for unit conversion?
      FLAG = ( MOD( ELAPSED_MIN, TS_DYN ) == 0 )

      END FUNCTION ITS_TIME_FOR_UNIT
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: its_time_for_diag
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
!EOP
!------------------------------------------------------------------------------
!BOC
      ! Is it time for diagnostics?
      FLAG = ( MOD( ELAPSED_MIN, TS_DIAG ) == 0 )

      END FUNCTION ITS_TIME_FOR_DIAG
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: its_time_for_a3
!
! !DESCRIPTION: Function ITS\_TIME\_FOR\_A3 returns TRUE if it is time to read 
!  in A-3 (average 3-h fields) and FALSE otherwise. 
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
      ! We read A-3 fields every 3 hours
      FLAG = ( MOD( NHMS, 030000 ) == 0 )

      END FUNCTION ITS_TIME_FOR_A3
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: its_time_for_a6
!
! !DESCRIPTION: Function ITS\_TIME\_FOR\_A6 returns TRUE if it is time to read 
!  in A-6 (average 6-h fields) and FALSE otherwise.
!\\
!\\
! !INTERFACE:
!
      FUNCTION ITS_TIME_FOR_A6() RESULT( FLAG )
!
! !USES:
!
#     include "define.h"
!
! !RETURN VALUE:
!
      LOGICAL :: FLAG 
! 
! !REVISION HISTORY: 
!  21 Mar 2003 - R. Yantosca - Initial Version
!  (1 ) Now compute when it's time to read in GEOS-4 A-6 fields. (bmy, 6/26/03)
!  (2 ) Now modified for GEOS-4 "a_llk_03" and "a_llk_04" fields (bmy, 3/22/04)
!  (3 ) Now modified for GCAP and GEOS-5 met fields (swu, bmy, 5/24/05)
!  15 Jan 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER :: DATE(2)

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

      END FUNCTION ITS_TIME_FOR_A6
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: its_time_for_i6
!
! !DESCRIPTION: Function ITS\_TIME\_FOR\_I6 returns TRUE if it is time to read 
!  in I-6 (instantaneous 6-h fields) and FALSE otherwise.
!\\
!\\
! !INTERFACE:
!
      FUNCTION ITS_TIME_FOR_I6() RESULT( FLAG )
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
      ! We read in I-6 fields at 00, 06, 12, 18 GMT
      FLAG = ( MOD( NHMS, 060000 ) == 0 )

      END FUNCTION ITS_TIME_FOR_I6
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: its_time_for_unzip
!
! !DESCRIPTION: Function ITS\_TIME\_FOR\_UNZIP Treturns TRUE if it is time to 
!  unzip the next day's met field files, or FALSE otherwise.
!\\
!\\
! !INTERFACE:
!
      FUNCTION ITS_TIME_FOR_UNZIP() RESULT( FLAG )
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
! !LOCAL VARIABLES:
!
      INTEGER :: DATE(2)

      ! Get YYYYMMDD and HHMMSS 12 hours (720 mins) from now
      DATE = GET_TIME_AHEAD( 720 )

      ! If HHMMSS = 0 then it's time to unzip!
      FLAG = ( DATE(2) == 000000 )

      END FUNCTION ITS_TIME_FOR_UNZIP     
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: its_time_for_del
!
! !DESCRIPTION: Function ITS\_TIME\_FOR\_DEL returns TRUE if it is time to 
!  delete the previous day's met field files in the temporary directory. 
!\\
!\\
! !INTERFACE:
!
      FUNCTION ITS_TIME_FOR_DEL() RESULT( FLAG )
!
! !RETURN VALUE:
!
      LOGICAL :: FLAG 

! 
! !REVISION HISTORY: 
!  21 Mar 2003 - R. Yantosca - Initial Version
!  19 Jun 2003 - R. Yantosca - Now delete files at 23 GMT each day, since the 
!                              last fvDAS A-3 field is 22:30 GMT and the last 
!                              fvDAS A-6 field is 21 GMT
!  15 Jan 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER :: DATE(2)

      ! Delete files when it's 23 GMT
      FLAG = ( NHMS == 230000 )

      END FUNCTION ITS_TIME_FOR_DEL
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: its_time_for_exit
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
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: its_time_for_bpch
!
! !DESCRIPTION: Function ITS\_TIME\_FOR\_BPCH returns TRUE if it's time to 
!  write output to the bpch file, or FALSE otherwise.
!\\
!\\
! !INTERFACE:
!
      FUNCTION ITS_TIME_FOR_BPCH() RESULT( DO_BPCH )
!
! !USES:
!
#     include "CMN_SIZE"  ! Size parameters
#     include "CMN_DIAG"  ! NJDAY
!
! !RETURN VALUE:
!
      LOGICAL :: DO_BPCH
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
         THIS_NJDAY = NJDAY( DOY + 1 ) 
      ELSE
         THIS_NJDAY = NJDAY( DOY )
      ENDIF

      ! Test if this is the day & time to write to the BPCH file!
      IF ( ( THIS_NJDAY > 0 ) .and. NHMS == NDIAGTIME ) THEN
         DO_BPCH = .TRUE.
      ELSE
         DO_BPCH = .FALSE.
      ENDIF

      END FUNCTION ITS_TIME_FOR_BPCH
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: its_a_leapyear
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

#if   defined( GCAP )
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

      END FUNCTION ITS_A_LEAPYEAR
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: its_a_new_year
!
! !DESCRIPTION: Function ITS\_A\_NEW\_YEAR returns TRUE if it's the first 
!  of a new month (it also returns TRUE on the first timestep of the run).  
!  This is useful for setting flags for reading in data.
!\\
!\\
! !INTERFACE:
!
      FUNCTION ITS_A_NEW_YEAR() RESULT( IS_NEW_YEAR )
!
! !RETURN VALUE:
!
      LOGICAL :: IS_NEW_YEAR
! 
! !REVISION HISTORY: 
!  01 Apr 2004 - R. Yantosca - Initial Version
!  01 Nov 2005 - R. Yantosca - Bug fix: Need month & day to be 1
!  15 Jan 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
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

      END FUNCTION ITS_A_NEW_YEAR
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: its_a_new_month
!
! !DESCRIPTION: Function ITS\_A\_NEW\_MONTH returns TRUE if it's the first 
!  of a new month (it also returns TRUE on the first timestep of the run).  
!  This is useful for setting flags for reading in data. 
!\\
!\\
! !INTERFACE:
!
      FUNCTION ITS_A_NEW_MONTH() RESULT( IS_NEW_MONTH )
!
! !RETURN VALUE:
!
      LOGICAL :: IS_NEW_MONTH
! 
! !REVISION HISTORY: 
!  01 Apr 2004 - R. Yantosca - Initial Version
!  15 Jan 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
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

      END FUNCTION ITS_A_NEW_MONTH
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: its_midmonth
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
!EOP
!------------------------------------------------------------------------------
!BOC

      ! Test for the 16th of the month at 0 GMT
      IS_MIDMONTH = ( DAY == 16 .and. NHMS == 000000 )

      END FUNCTION ITS_MIDMONTH
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: its_a_new_season
!
! !DESCRIPTION: Function ITS\_A\_NEW\_DAY returns TRUE if it's the first 
!  timestep of a new day (it also returns TRUE on the first timestep of the
!  run).  This is useful for setting flags for reading in data.
!\\
!\\
! !INTERFACE:
!
      FUNCTION ITS_A_NEW_DAY( ) RESULT( IS_NEW_DAY )
!
! !RETURN VALUE:
!
      LOGICAL :: IS_NEW_DAY
!
! !REVISION HISTORY: 
!  01 Apr 2004 - R. Yantosca - Initial Version
!  15 Jan 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
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

      END FUNCTION ITS_A_NEW_DAY
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: its_a_new_season
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
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: print_current_time
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
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL*4 :: E_HOURS

      ! Hours since start of run
      E_HOURS = REAL( ELAPSED_MIN ) / 60e0 

      ! Write quantities
      WRITE( 6, 100 ) YEAR, MONTH, DAY, HOUR, MINUTE, E_HOURS

      ! Format string
 100  FORMAT( '---> DATE: ', i4.4, '/', i2.2, '/', i2.2, 
     &            '  GMT: ', i2.2, ':', i2.2, '  X-HRS: ', f11.3 )

      END SUBROUTINE PRINT_CURRENT_TIME
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: timestamp_string
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
! !USES:
!
#     include "define.h"
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
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ymd_extract
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
      REAL*8 :: REM

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
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: expand_date
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

#     include "define.h"
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
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: system_date_time
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
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: system_timestamp
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
!          Harvard University Atmospheric Chemistry Modeling Group            !
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
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_nymd_diag
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


      END MODULE TIME_MOD




