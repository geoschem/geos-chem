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
!
! !DEFINED PARAMETERS:
!
  !-------------------------------------------------------------------------
  ! MISSING DATA VALUES:
  !
  ! Specify missing data values for various numeric types.
  !-------------------------------------------------------------------------
  INTEGER,          PARAMETER, PUBLIC :: UNDEFINED_INT     = -999
  REAL(f4),         PARAMETER, PUBLIC :: UNDEFINED         = -1.0e+31_f4
  REAL(f8),         PARAMETER, PUBLIC :: UNDEFINED_DBL     = -1.0e+31_f8
  CHARACTER(LEN=9), PARAMETER, PUBLIC :: UNDEFINED_STR     = 'not found'

  !-------------------------------------------------------------------------
  ! OPERATION CODES:
  !
  ! 0 = Copy       data from source pointer to the HISTORY ITEM data array
  ! 1 = Accumulate data from source pointer to the HISTORY ITEM data array
  !-------------------------------------------------------------------------
  INTEGER,          PARAMETER, PUBLIC :: COPY_FROM_SOURCE  = 0
  INTEGER,          PARAMETER, PUBLIC :: ACCUM_FROM_SOURCE = 1

  !-------------------------------------------------------------------------
  ! ALARM CODES:
  ! Used to define the alarms that denote when to perform a given action
  ! 
  ! 0 = Update (aka archive) data from the source into the HISTORY ITEM
  ! 1 = Write data to the netCDF file
  ! 2 = Close the current netCDF file and reopen it for the next interval
  !-------------------------------------------------------------------------
  INTEGER,          PARAMETER, PUBLIC :: ALARM_UPDATE      = 0
  INTEGER,          PARAMETER, PUBLIC :: ALARM_FILE_WRITE  = 1
  INTEGER,          PARAMETER, PUBLIC :: ALARM_FILE_CLOSE  = 2

  !-------------------------------------------------------------------------
  ! ELAPSED TIME CODES: 
  ! Define the base time from which to compute the time elapsed in minutes.  
  !
  ! 0 = Since start of run
  ! 1 = Since netCDF file creation (i.e. ReferenceYmd, ReferenceHms)
  !-------------------------------------------------------------------------
  INTEGER,          PARAMETER, PUBLIC :: FROM_START_OF_SIM = 0
  INTEGER,          PARAMETER, PUBLIC :: FROM_FILE_CREATE  = 1

  !-------------------------------------------------------------------------
  ! ROUNDOFF_DECIMALS:
  ! Specifies the number of decimal digits for rounding, 
  ! when converting Julian Dates to elapsed time values.
  !-------------------------------------------------------------------------
  INTEGER,          PARAMETER, PUBLIC :: ROUNDOFF_DECIMALS = 5

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
!  18 Aug 2017 - R. Yantosca - Renamed ACTION_* parameter to ALARM_*
!  18 Aug 2017 - R. Yantosca - Add FROM_* parameters for elapsed time
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
    FracDay = DBLE( Day ) + ( DBLE( Hour   ) / 24.0_f8    )  +               & 
                            ( DBLE( Minute ) / 1440.0_f8  )  +               &
                            ( DBLE( Second ) / 86400.0_f8 ) 

    ! Return the Astronomical Julian Date
    Jd = JulDay( Year, Month, FracDay )

  END SUBROUTINE Compute_Julian_Date
!EOC
END MODULE History_Util_Mod
