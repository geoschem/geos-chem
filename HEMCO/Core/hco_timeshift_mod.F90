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
    ! Year and month are placed in slot 1. Years are 
    ! converted to months.
    IDX = INDEX( TRIM(iShift), 'year' )
    IF( IDX < 0 ) IDX = INDEX( TRIM(iShift), 'yr' )
    IF ( IDX > 1 ) THEN
       READ( iShift(1:(IDX-1)), * ) SHFT
       Dta%tShift(1) = Dta%tShift(1) + &
                             ( SHFT * 12 )
    ENDIF

    IDX = INDEX( TRIM(iShift), 'month' )
    IF( IDX < 0 ) IDX = INDEX( TRIM(iShift), 'mt' )
    IF ( IDX > 1 ) THEN
       READ( iShift(1:(IDX-1)), * ) SHFT
       Dta%tShift(1) = Dta%tShift(1) + SHFT 
    ENDIF

    ! Days, Hours, and minutes are placed in slot 2.
    ! All values are converted to seconds. 
    IDX = INDEX( TRIM(iShift), 'day' )
    IF( IDX < 0 ) IDX = INDEX( TRIM(iShift), 'dy' )
    IF ( IDX > 1 ) THEN
       READ( iShift(1:(IDX-1)), * ) SHFT
       Dta%tShift(2) = Dta%tShift(2) + ( SHFT * 86400 )
                             
    ENDIF

    IDX = INDEX( TRIM(iShift), 'hour' )
    IF( IDX < 0 ) IDX = INDEX( TRIM(iShift), 'hr' )
    IF ( IDX > 1 ) THEN
       READ( iShift(1:(IDX-1)), * ) SHFT
       Dta%tShift(2) = Dta%tShift(2) + ( SHFT * 3600 )
                             
    ENDIF

    IDX = INDEX( TRIM(iShift), 'min' )
    IF( IDX < 0 ) IDX = INDEX( TRIM(iShift), 'mn' )
    IF ( IDX > 1 ) THEN
       READ( iShift(1:(IDX-1)), * ) SHFT
       Dta%tShift(2) = Dta%tShift(2) + ( SHFT * 60 )
                             
    ENDIF

    ! Error: cannot shift data if we use weekday data
    IF ( Dta%tShift(1) /= 0 .OR. Dta%tShift(2) /= 0 ) THEN
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
       WRITE(MSG,*) 'Will shift time stamp of field ', TRIM(Dta%ncFile), ': '
       CALL HCO_MSG(HcoConfig%Err,MSG)
       WRITE(MSG,*) 'Time shift (months, seconds): ', Dta%tShift(1), &
                                                      Dta%tShift(2)
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
  SUBROUTINE TimeShift_Apply( am_I_Root, HcoState, Lct, &
                              yr, mt, dy, hr, mn,  RC )
!
! !USES:
!
    USE Julday_Mod
    USE HCO_TYPES_MOD,    ONLY : ListCont
!
! !INPUT/OUTPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN   )  :: am_I_Root ! Root CPU
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
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER                        :: nYr, nMt, nDy, nHr, nMn
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
         Lct%Dct%Dta%tShift(2) == 0        ) RETURN

    ! Mirror values
    nYr = MAX(Yr,1)
    nMt = MAX(Mt,1)
    nDy = MAX(Dy,1)
    nHr = MAX(Hr,0)
    nMn = MAX(Mn,0)

    ! Get time shift in days
    TimeShift = REAL(Lct%Dct%Dta%tShift(2),dp) / 86400.0_dp

    ! Add monthly time shift
    IF ( Lct%Dct%Dta%tShift(1) > 0 ) THEN
       nMt = nMt + Lct%Dct%Dta%tShift(1)
       ! Make sure new value is within 1-12. Also 
       ! adjust year accordingly
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
    ENDIF

    ! Get current date as Julian day.
    UTC = ( REAL(nHr,dp) / 24.0_dp    ) + &
          ( REAL(nMn,dp) / 1440.0_dp  ) + &
          ( REAL(0  ,dp) / 86400.0_dp )
    DAY = REAL(nDy,dp) + UTC
    JD  = JULDAY( nYr, nMt, DAY )

    ! Add time shift in seconds
    JD = JD + TimeShift

    ! Translate back into dates.
    CALL CALDATE( JD, YYYYMMDD, HHMMSS )
    nYr = FLOOR ( MOD( YYYYMMDD, 100000000) / 1.0e4_dp )
    nMt = FLOOR ( MOD( YYYYMMDD, 10000    ) / 1.0e2_dp )
    nDy = FLOOR ( MOD( YYYYMMDD, 100      ) / 1.0e0_dp )

    nHr = FLOOR ( MOD(   HHMMSS, 1000000  ) / 1.0e4_dp )
    nMn = FLOOR ( MOD(   HHMMSS, 10000    ) / 1.0e2_dp )

    ! verbose mode
    IF ( HCO_IsVerb(HcoState%Config%Err,3) ) THEN
       WRITE(MSG,*) 'Adjusted time stamp of field ', TRIM(Lct%Dct%cName), ': '
       CALL HCO_MSG(HcoState%Config%Err,MSG)
       WRITE(MSG,*) 'Time shift (months, seconds): ', Lct%Dct%Dta%tShift(1), &
                                                      Lct%Dct%Dta%tShift(2)
       CALL HCO_MSG(HcoState%Config%Err,MSG)
       WRITE(MSG,'(a27,i4.4,a1,i2.2,a1,i2.2,a1,i2.2,a1,i2.2)') 'Original Yr/Mt/Dy-Hr:Mn = ',Yr,'/',Mt,'/',Dy,'-',Hr,':',Mn
       CALL HCO_MSG(HcoState%Config%Err,MSG)
    ENDIF 

    ! Add back to output values
    IF ( Yr >  0 ) Yr = nYr
    IF ( Mt >  0 ) Mt = nMt
    IF ( Dy >  0 ) Dy = nDy
    IF ( Hr > -1 ) Hr = nHr
    IF ( Mn > -1 ) Mn = nMn

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
