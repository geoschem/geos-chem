!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: hco_timeshift_mod.F90
!
! !DESCRIPTION: Module hco\_timeshift\_mod.F90 contains routines to shift the
! file reference time by a given value. Time stamps shifts can be provided as
! optional fifth element to the time stamp attribute in the HEMCO configuration
! file.
!\\
!\\
! For instance, consider the case where 3-hourly averages are provided in
! individual files with centered time stamps, e.g.: file.yyyymmdd\_0130z.nc,
! file.yyyymmdd\_0430z.nc, ..., file.yyymmdd\_2230z.nc
! To read these files *at the beginning* of their time intervals, the time
! stamp can be shifted by 90 minutes, e.g. the file name, variable, and time
! attribute section reads:
! ... file.\$yyyy$mm$dd\_\$hh\$mnz.nc VARNAME 2000-2016/1-12/1-31/0-23/+90minutes ...
!\\
!\\
! At time 00z, HEMCO will then read file 0130z and keep using this file until
! 03z, when it switches to file 0430z. Similarly, it is possible to shift the
! file reference time by any number of years, months, days, or hours. Time
! shifts can be forward or backward in time (use - sign to shift backwards).
!\\
!\\
! This module contains subroutines to determine the time to be shifted (stored
! in filedata variable tshift) and shifts the desired reference time as needed.
!\\
!\\
! !INTERFACE:
!
MODULE HCO_TIMESHIFT_MOD
!
! !USES:
!
  USE HCO_Error_Mod
  USE HCO_State_Mod,  ONLY : HCO_State

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: TimeShift_Set
  PUBLIC :: TimeShift_Apply
!
! !REMARKS:
!
! !REVISION HISTORY:
!  29 Feb 2016 - C. Keller   - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !PRIVATE TYPES:
!
CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: TimeShift_Set
!
! !DESCRIPTION: Subroutine TimeShift\_Set sets the time shift values. The
! time shift attribute tshift contains two entries: the first entry denotes
! the number of months to be shifted (integer value), while the second entry
! denotes then number of seconds to be shifted.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE TimeShift_Set( HcoConfig, Dta, shift, RC )
!
! !USES:
!
    USE HCO_TYPES_MOD,    ONLY : ListCont
    USE HCO_TYPES_MOD,    ONLY : FileData
    USE HCO_TYPES_MOD,    ONLY : ConfigObj
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ConfigObj),  POINTER        :: HcoConfig    ! HEMCO config
    TYPE(FileData),   POINTER        :: Dta       ! file container
    CHARACTER(LEN=*), INTENT(IN   )  :: shift     ! time shift
    INTEGER,          INTENT(INOUT)  :: RC        ! Return code
!
! !REVISION HISTORY:
!  29 Feb 2016 - C. Keller - Initial version
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER                        :: SHFT, IDX
    CHARACTER(LEN=255)             :: MSG
    CHARACTER(LEN=255)             :: iShift
    CHARACTER(LEN=255)             :: tShift
    CHARACTER(LEN=255), PARAMETER  :: LOC = 'TimeShift_Set (hco_tShift_mod.F90)'

    !======================================================================
    ! TimeShift_Set begins here!
    !======================================================================

    ! Assume success until otherwise
    RC = HCO_SUCCESS

    ! Init
    Dta%tShift(:) = 0.0_hp

    ! Mirror variable
    iShift = TRIM(ADJUSTL(Shift))

    ! Parse time iShift and write to desired slot:
    IDX = INDEX( TRIM(iShift), 'year' )
    IF( IDX < 0 ) IDX = INDEX( TRIM(iShift), 'yr' )
    IF ( IDX > 1 ) THEN
       READ( iShift(1:(IDX-1)), * ) SHFT
       Dta%tShift(1) = Dta%tShift(1) + SHFT
       WRITE(tShift,*) Dta%tShift(1), ' years'
    ENDIF

    IDX = INDEX( TRIM(iShift), 'month' )
    IF( IDX < 0 ) IDX = INDEX( TRIM(iShift), 'mt' )
    IF ( IDX > 1 ) THEN
       READ( iShift(1:(IDX-1)), * ) SHFT
       Dta%tShift(2) = Dta%tShift(2) + SHFT
       WRITE(tShift,*) Dta%tShift(2), ' months'
    ENDIF

    IDX = INDEX( TRIM(iShift), 'day' )
    IF( IDX < 0 ) IDX = INDEX( TRIM(iShift), 'dy' )
    IF ( IDX > 1 ) THEN
       READ( iShift(1:(IDX-1)), * ) SHFT
       Dta%tShift(3) = Dta%tShift(3) + SHFT
       WRITE(tShift,*) Dta%tShift(3), ' days'
    ENDIF

    IDX = INDEX( TRIM(iShift), 'hour' )
    IF( IDX < 0 ) IDX = INDEX( TRIM(iShift), 'hr' )
    IF ( IDX > 1 ) THEN
       READ( iShift(1:(IDX-1)), * ) SHFT
       Dta%tShift(4) = Dta%tShift(4) + SHFT
       WRITE(tShift,*) Dta%tShift(4), ' hours'
    ENDIF

    IDX = INDEX( TRIM(iShift), 'min' )
    IF( IDX < 0 ) IDX = INDEX( TRIM(iShift), 'mn' )
    IF ( IDX > 1 ) THEN
       READ( iShift(1:(IDX-1)), * ) SHFT
       Dta%tShift(5) = Dta%tShift(5) + SHFT
       WRITE(tShift,*) Dta%tShift(5), ' minutes'
    ENDIF

    IDX = INDEX( TRIM(iShift), 'sec' )
    IF ( IDX > 1 ) THEN
       READ( iShift(1:(IDX-1)), * ) SHFT
       Dta%tShift(6) = Dta%tShift(6) + SHFT
       WRITE(tShift,*) Dta%tShift(6), ' seconds'
    ENDIF

    ! Error: cannot shift data if we use weekday data
    IF ( Dta%tShift(1) /= 0 .OR. Dta%tShift(2) /= 0 .or. &
         Dta%tShift(3) /= 0 .OR. Dta%tShift(4) /= 0 .or. &
         Dta%tShift(5) /= 0 .OR. Dta%tShift(6) /= 0 ) THEN
       IF (   Dta%ncDys(1) == -10 .OR. &
            ( Dta%ncDys(1) == 1 .AND. Dta%ncDys(2) == 7 ) ) THEN
          WRITE(MSG,*) 'Time shift not supported for weekday data: ', &
             TRIM(Dta%ncFile)
          CALL HCO_ERROR( HcoConfig%Err, MSG, RC, THISLOC=LOC )
          RETURN
       ENDIF
    ENDIF

    ! verbose mode
    IF ( HCO_IsVerb(HcoConfig%Err,2) ) THEN
       WRITE(MSG,*) 'Will shift time stamp of field ', TRIM(Dta%ncPara), &
                    ': ', TRIM(tShift)
       CALL HCO_MSG(HcoConfig%Err,MSG)
    ENDIF

    RC = HCO_SUCCESS

  END SUBROUTINE TimeShift_Set
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: TimeShift_Apply
!
! !DESCRIPTION: Subroutine TimeShift\_Apply shifts the reference time
! (provided through arguments yr, mt, dy, hr, and mn, by the time shift
! specified in the HEMCO configuration file (if specified at all).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE TimeShift_Apply( HcoState, Lct, &
                              yr, mt, dy, hr, mn,  RC )
!
! !USES:
!
    USE Julday_Mod
    USE HCO_TYPES_MOD,    ONLY : ListCont
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(HCO_State),  POINTER        :: HcoState  ! Hemco state
    TYPE(ListCont),   POINTER        :: Lct       ! List container
    INTEGER,          INTENT(INOUT)  :: yr        ! year
    INTEGER,          INTENT(INOUT)  :: mt        ! month
    INTEGER,          INTENT(INOUT)  :: dy        ! day
    INTEGER,          INTENT(INOUT)  :: hr        ! hour
    INTEGER,          INTENT(INOUT)  :: mn        ! minute
    INTEGER,          INTENT(INOUT)  :: RC        ! Return code
!
! !REVISION HISTORY:
!  29 Feb 2016 - C. Keller - Initial version
!  19 Nov 2018 - C. Keller - Add option TimeShiftCap
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER                        :: nYr, nMt, nDy, nHr, nMn
    INTEGER                        :: oYr, oMt, oDy, oHr, oMn
    INTEGER                        :: SHFT, IDX
    INTEGER                        :: YYYYMMDD, HHMMSS
    REAL(dp)                       :: TimeShift, DAY, UTC, JD

    CHARACTER(LEN=255)             :: MSG
    CHARACTER(LEN=255), PARAMETER  :: LOC = 'TimeShift_Apply (hco_tShift_mod.F90)'

    !======================================================================
    ! TimeShift_Apply begins here!
    !======================================================================

    ! Assume success until otherwise
    RC = HCO_SUCCESS

    ! Nothing to do if time shift is zero
    IF ( Lct%Dct%Dta%tShift(1) == 0 .AND. &
         Lct%Dct%Dta%tShift(2) == 0 .AND. &
         Lct%Dct%Dta%tShift(3) == 0 .AND. &
         Lct%Dct%Dta%tShift(4) == 0 .AND. &
         Lct%Dct%Dta%tShift(5) == 0 .AND. &
         Lct%Dct%Dta%tShift(6) == 0 ) RETURN

    ! Mirror values
    nYr = MAX(Yr,1)
    nMt = MAX(Mt,1)
    nDy = MAX(Dy,1)
    nHr = MAX(Hr,0)
    nMn = MAX(Mn,0)

    ! Archive original values
    oYr = nYr
    oMt = nMt
    oDy = nDy
    oHr = nHr
    oMn = nMn

    ! Apply minute time shift
    IF ( ABS(Lct%Dct%Dta%tShift(5)) > 0 ) THEN
       nMn = nMn + Lct%Dct%Dta%tShift(5)
    ENDIF

    ! Make sure new minute value is within 0-59 and adjust hour accordingly
    IF ( nMn > 59 ) THEN
       DO WHILE ( nMn > 59 )
          nMn = nMn - 60
          nHr = nHr + 1
       ENDDO
    ELSEIF( nMn < 0 ) THEN
       DO WHILE ( nMn < 1 )
          nMn = nMn + 60
          nHr = nHr - 1
       ENDDO
    ENDIF

    ! Apply hour time shift
    IF ( ABS(Lct%Dct%Dta%tShift(4)) > 0 ) THEN
       nHr = nHr + Lct%Dct%Dta%tShift(4)
    ENDIF

    ! Make sure new hour value is within 0-23 and adjust day accordingly
    IF ( nHr > 23 ) THEN
       DO WHILE ( nHr > 23 )
          nHr = nHr - 24
          nDy = nDy + 1
       ENDDO
    ELSEIF( nHr < 0 ) THEN
       DO WHILE ( nHr < 1 )
          nHr = nHr + 24
          nDy = nDy - 1
       ENDDO
    ENDIF

    ! Apply day time shift
    IF ( ABS(Lct%Dct%Dta%tShift(3)) > 0 ) THEN
       nDy = nDy + Lct%Dct%Dta%tShift(3)
    ENDIF

    ! Make sure new day is within ndays in month and adjust month accordingly
    IF ( nDy > 28 .AND. oMt == 2 .AND. MOD( oYr, 4 ) /= 0 ) THEN
       nDy = nDy - 28
       nMt = nMt + 1
    ELSEIF ( nDy > 29 .AND. oMt == 2 .AND. MOD( oYr, 4 ) == 0) THEN
       nDy = nDy - 29
       nMt = nMt + 1
    ELSEIF ( nDy > 30 .AND. &
           ( oMt == 4 .OR. oMt == 6 .OR. oMt == 9 .OR. oMt == 11 ) ) THEN
       nDy = nDy - 30
       nMt = nMt + 1
    ELSEIF ( nDy > 31 .AND. &
           ( oMt == 1 .OR. oMt == 3  .OR. oMt == 5 .OR. oMt == 7 .OR. &
             oMt == 8 .OR. oMt == 10 .OR. oMt == 12 ) ) THEN
       nDy = nDy - 31
       nMt = nMt + 1
    ELSEIF ( nDy <= 0 ) THEN
       IF ( oMt == 3 ) THEN
          IF ( MOD( oYr, 4 ) == 0 ) THEN
             nDy = nDy + 29
          ELSE
             nDy = nDy + 28
          ENDIF
          nMt = nMt - 1
       ELSEIF ( oMt == 5 .OR. oMt == 7 .OR. oMt == 10 .OR. oMt == 12 ) THEN
          nDy = nDy + 30
          nMt = nMt - 1
       ELSEIF ( oMt == 1. .OR. oMt == 2 .OR. oMt == 4 .OR. oMt == 6 .OR. &
                oMt == 8  .OR. oMt == 9 .OR. oMt == 11 ) THEN
          nDy = nDy + 31
          nMt = nMt - 1
       ENDIF
    ENDIF

    ! Apply month time shift
    IF ( ABS(Lct%Dct%Dta%tShift(2)) > 0 ) THEN
       nMt = nMt + Lct%Dct%Dta%tShift(2)
    ENDIF

    ! Make sure new month is within 1-12 and adjust year accordingly
    IF ( nMt > 12 ) THEN
       DO WHILE ( nMt > 12 )
          nMt = nMt - 12
          nYr = nYr + 1
       ENDDO
    ELSEIF( nMt <= 0 ) THEN
       DO WHILE ( nMt < 1 )
          nMt = nMt + 12
          nYr = nYr - 1
       ENDDO
    ENDIF

    ! Apply year time shift
    IF ( ABS(Lct%Dct%Dta%tShift(1)) > 0 ) THEN
       nYr = nYr + Lct%Dct%Dta%tShift(1)
    ENDIF

    ! Eventually cap time shift for the same day.
    IF ( HcoState%Options%TimeShiftCap ) THEN
       IF ( ( nYr < oYr                  ) .OR. &
            ( nMt < oMt .AND. nYr == oYr ) .OR. &
            ( nDy < oDy .AND. nMt == oMt )       ) THEN
          nYr = oYr
          nMt = oMt
          nDy = oDy
          nHr = 0
          nMn = 0
          IF ( HCO_IsVerb(HcoState%Config%Err,3) ) THEN
             MSG = 'Options set to cap time shift - set to low bound'
             CALL HCO_MSG(HcoState%Config%Err,MSG)
          ENDIF
       ELSEIF ( nYr > oYr .OR. nMt > oMt .OR. nDy > oDy ) THEN
          nYr = oYr
          nMt = oMt
          nDy = oDy
          nHr = 23
          nMn = 59
          IF ( HCO_IsVerb(HcoState%Config%Err,3) ) THEN
             MSG = 'Options set to cap time shift - set to high bound'
             CALL HCO_MSG(HcoState%Config%Err,MSG)
          ENDIF
       ENDIF
    ENDIF

    ! verbose mode
    IF ( HCO_IsVerb(HcoState%Config%Err,3) ) THEN
       WRITE(MSG,*) 'Adjusted time stamp of field ', TRIM(Lct%Dct%cName)
       CALL HCO_MSG(HcoState%Config%Err,MSG)
       WRITE(MSG,*) 'Time shift (YMDhms): ', Lct%Dct%Dta%tShift
       CALL HCO_MSG(HcoState%Config%Err,MSG)
       WRITE(MSG,'(a27,i4.4,a1,i2.2,a1,i2.2,a1,i2.2,a1,i2.2)') 'Original Yr/Mt/Dy-Hr:Mn = ',oYr,'/',oMt,'/',oDy,'-',oHr,':',oMn
       CALL HCO_MSG(HcoState%Config%Err,MSG)
    ENDIF

    ! Add back to output values
    Yr = nYr
    Mt = nMt
    Dy = nDy
    Hr = nHr
    Mn = nMn

    ! verbose mode
    IF ( HCO_IsVerb(HcoState%Config%Err,3) ) THEN
       WRITE(MSG,'(a27,i4.4,a1,i2.2,a1,i2.2,a1,i2.2,a1,i2.2)') 'Adjusted Yr/Mt/Dy-Hr:Mn = ',Yr,'/',Mt,'/',Dy,'-',Hr,':',Mn
       CALL HCO_MSG(HcoState%Config%Err,MSG)
    ENDIF

    ! Return w/ success
    RC = HCO_SUCCESS

  END SUBROUTINE TimeShift_Apply
!EOC
END MODULE HCO_TIMESHIFT_MOD
