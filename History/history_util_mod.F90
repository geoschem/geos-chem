!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: history_util_mod.F90
!
! !DESCRIPTION: Contains defined parameters and utility routines for
!  the GEOS-Chem History Component.
!\\
!\\
! !INTERFACE:
!
MODULE History_Util_Mod
!
! !USES:
!
  USE Precision_Mod

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: Compute_Julian_Date
  PUBLIC :: Compute_Elapsed_Time
  PUBLIC :: Compute_DeltaYmdHms_For_End
!
! !DEFINED PARAMETERS:
!
  !-------------------------------------------------------------------------
  ! MISSING DATA VALUES:
  !
  ! Specify missing data values for various numeric types.
  !-------------------------------------------------------------------------
  INTEGER,          PARAMETER, PUBLIC :: UNDEFINED_INT      = -999
  REAL(f4),         PARAMETER, PUBLIC :: UNDEFINED          = -1.0e+31_f4
  REAL(f8),         PARAMETER, PUBLIC :: UNDEFINED_DBL      = -1.0e+31_f8
  CHARACTER(LEN=9), PARAMETER, PUBLIC :: UNDEFINED_STR      = 'not found'

  !-------------------------------------------------------------------------
  ! OPERATION CODES:
  !
  ! 0 = Copy       data from source pointer to the HISTORY ITEM data array
  ! 1 = Accumulate data from source pointer to the HISTORY ITEM data array
  !-------------------------------------------------------------------------
  INTEGER,          PARAMETER, PUBLIC :: COPY_FROM_SOURCE   = 0
  INTEGER,          PARAMETER, PUBLIC :: ACCUM_FROM_SOURCE  = 1

  !-------------------------------------------------------------------------
  ! ROUNDING AND NUMERICAL TESTING PARAMETRS
  ! Specifies the number of decimal digits for rounding, as well as an
  ! epsilon value that can be used for floating point equality testing.
  !-------------------------------------------------------------------------
  INTEGER,          PARAMETER, PUBLIC :: ROUNDOFF_DECIMALS  = 4
  REAL(f8),         PARAMETER, PUBLIC :: EPS                = 1e-5_f8

  !-------------------------------------------------------------------------
  ! TIME CONVERSION PARAMETERS
  ! Specifies the number of minutes and seconds per day, etc.
  !-------------------------------------------------------------------------
  REAL(f8),         PARAMETER, PUBLIC :: HOURS_PER_DAY      = 24.0_f8
  REAL(f8),         PARAMETER, PUBLIC :: MINUTES_PER_DAY    = 1440.0_f8
  REAL(f8),         PARAMETER, PUBLIC :: MINUTES_PER_HOUR   = 60.0_f8
  REAL(f8),         PARAMETER, PUBLIC :: SECONDS_PER_DAY    = 86400.0_f8
  REAL(f8),         PARAMETER, PUBLIC :: SECONDS_PER_HOUR   = 3600.0_f8
  REAL(f8),         PARAMETER, PUBLIC :: SECONDS_PER_MINUTE = 60.0_f8
!
! !REVISION HISTORY:
!  16 Jun 2017 - R. Yantosca - Initial version
!  03 Aug 2017 - R. Yantosca - Add operation code parameters
!  04 Aug 2017 - R. Yantosca - Rename operation code parameters
!  09 Aug 2017 - R. Yantosca - Add UNDEFINED_INT
!  15 Aug 2017 - R. Yantosca - Add UNDEFINED_STR
!  16 Aug 2017 - R. Yantosca - Added ACTION_* parameters
!  17 Aug 2017 - R. Yantosca - Renamed to history_util_mod.F90; added routine
!                              Compute_Julian_Date
!  21 Aug 2017 - R. Yantosca - Removed some parameters that weren't needed
!EOP
!------------------------------------------------------------------------------
!BOC
CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Compute_Julian_Date
!
! !DESCRIPTION: Computes the Astronomical Julian Date corresponding to a
!  given date and time.  This is useful for computing elapsed times.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Compute_Julian_Date( yyyymmdd, hhmmss, Jd )
!
! !USES:
!
    USE Julday_Mod, ONLY : Julday
    USE Time_Mod,   ONLY : Ymd_Extract
!
! !INPUT PARAMETERS:
!
    INTEGER,  INTENT(IN)  :: yyyymmdd   ! Current Year/month/day
    INTEGER,  INTENT(IN)  :: hhmmss     ! Current hour/minute/second
!
! !OUTPUT PARAMETERS:
!
    REAL(f8), INTENT(OUT) :: Jd         ! Astronomical Julian date
!
! !REMARKS:
!  This is a convenience wrapper for the JULDAY routine, which is located
!  in GeosUtil/julday_mod.F.
!
! !REVISION HISTORY:
!  06 Jan 2015 - R. Yantosca - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER  :: Year, Month, Day, Hour, Minute, Second
    REAL(f8) :: FracDay

    ! Extract year/month/day and hour/minute/seconds from the time
    CALL Ymd_Extract( yyyymmdd, Year, Month,  Day    )
    CALL Ymd_Extract( hhmmss,   Hour, Minute, Second )

    ! Compute the fractional day
    FracDay = DBLE( Day ) + ( DBLE( Hour   ) / HOURS_PER_DAY   )  +          &
                            ( DBLE( Minute ) / MINUTES_PER_DAY )  +          &
                            ( DBLE( Second ) / SECONDS_PER_DAY )

    ! Return the Astronomical Julian Date
    Jd = JulDay( Year, Month, FracDay )

  END SUBROUTINE Compute_Julian_Date
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Compute_Elapsed_Time
!
! !DESCRIPTION: Computes elapsed time in minutes, given the current
!  Astronomical Julian Date value, plus a reference Julian Date.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Compute_Elapsed_Time( CurrentJsec, TimeBaseJsec, ElapsedSec )
!
! !INPUT PARAMETERS:
!
    REAL(f8), INTENT(IN)  :: CurrentJsec  ! Current astronomical Julian date
    REAL(f8), INTENT(IN)  :: TimeBaseJsec ! Reference astronomical Julian date
!
! !OUTPUT PARAMETERS:
!
    REAL(f8), INTENT(OUT) :: ElapsedSec   ! Elapsed time [seconds]
!
! !REMARKS:
!  The netCDF file reference date and time are given by the ReferenceYmd,
!  ReferenceHms, and ReferenceJd fields of the Container object.  This
!  denotes the simulation date & time when the netCDF file was created.
!
! !REVISION HISTORY:
!  18 Aug 2017 - R. Yantosca - Initial version
!  13 Sep 2017 - R. Yantosca - Avoid roundoff error; return integral minutes
!  18 Sep 2017 - R. Yantosca - Now return elapsed seconds to avoid roundoff
!  29 Sep 2017 - R. Yantosca - Use NINT instead of INT to avoid roundoff
!  11 Jul 2018 - R. Yantosca - Now accept input arguments in seconds
!EOP
!------------------------------------------------------------------------------
!BOC

    !=======================================================================
    ! Compute elapsed time in minutes
    !=======================================================================

    ! Compute elapsed minutes since start of simulation
    ElapsedSec = ( CurrentJsec - TimeBaseJsec )

    ! Just keep the integer part, since we are dealing in integral seconds
    ! NINT ensures that we round up in case there is underflow
    ElapsedSec = NINT( ElapsedSec )

  END SUBROUTINE Compute_Elapsed_Time
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Compute_DeltaYmdHms_For_End
!
! !DESCRIPTION: Returns the DeltaYMD and DeltaHMS parameters for collections
!  that are specified with "End".
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Compute_DeltaYmdHms_For_End( yyyymmdd,     hhmmss,              &
                                          yyyymmdd_end, hhmmss_end,          &
                                          deltaYmd,     deltaHms            )
!
! !USES:
!
    USE Time_Mod, ONLY : Ymd_Extract
!
! !INPUT PARAMETERS:
!
    INTEGER, INTENT(IN)  :: yyyymmdd       ! Simulation start date
    INTEGER, INTENT(IN)  :: hhmmss         ! Simulation start time
    INTEGER, INTENT(IN)  :: yyyymmdd_end   ! Simulation end date
    INTEGER, INTENT(IN)  :: hhmmss_end     ! Simulation end time
!
! !OUTPUT PARAMETERS:
!
    INTEGER, INTENT(OUT) :: deltaYmd       ! YMD interval for the collection
    INTEGER, INTENT(OUT) :: deltaHms       ! HMS interval for the collection
!
! !REMARKS:
!  NOTE: This algorithm should work for most typical model start and end dates
!  (which are usually integral intervals of months, days, hours, or minutes.
!  There may be some edge cases that will cause this to fail.  But it is an
!  improvement over the prior situation. (bmy, 2/26/19)
!
! !REVISION HISTORY:
!  26 Feb 2019 - R. Yantosca - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: Year0, Month0, Day0, Hour0, Minute0, Second0
    INTEGER :: Year1, Month1, Day1, Hour1, Minute1, Second1
    INTEGER :: dYear, dMonth, dDay, dHour, dMinute, dSecond

    !=======================================================================
    ! GetRestartDeltaYmdHms begins here!
    !=======================================================================

    ! Split starting date
    CALL Ymd_Extract( yyyymmdd,     Year0, Month0,  Day0    )
    CALL Ymd_Extract( hhmmss,       Hour0, Minute0, Second0 )

    ! Split ending date
    CALL Ymd_Extract( yyyymmdd_end, Year1, Month1,  Day1    )
    CALL Ymd_Extract( hhmmss_end,   Hour1, Minute1, Second1 )

    ! Compute intervals
    dYear    = Year1   - Year0
    dMonth   = Month1  - Month0
    dDay     = Day1    - Day0
    dHour    = Hour1   - Hour0
    dMinute  = Minute1 - Minute0
    dSecond  = Second1 - Second0

    ! Adjust intervals (NOTE: Might not be as robust, more testing needed)
    IF ( dSecond < 0 ) THEN
       dSecond = dSecond + 60
       dMinute = MAX( dMinute - 1, 0 )
    ENDIF

    IF ( dMinute < 0 ) THEN
       dMinute = dMinute + 60
       dHour   = MAX( dHour - 1, 0 )
    ENDIF

    IF ( dHour < 0 ) THEN
       dHour   = dHour + 24
       dDay    = MAX( dDay - 1, 0 )
    ENDIF

    IF ( dDay < 0 ) THEN
       dDay    = dDay + Day0
       dMonth  = MAX( dMonth - 1, 0 )
    ENDIF

    IF ( dMonth < 0 ) THEN
       dMonth  = dMonth + 12
       dYear   = MAX( dYear - 1, 0 )
    ENDIF

    ! Construct the YMD and HMS intervals
    deltaYmd = ( dYear * 10000 ) + ( dMonth  * 100 ) + dDay
    deltaHms = ( dHour * 10000 ) + ( dMinute * 100 ) + dSecond

  END SUBROUTINE Compute_DeltaYmdHms_For_End
!EOC
END MODULE History_Util_Mod
