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
END MODULE History_Util_Mod
