! $Id: soiltype.f,v 1.1 2009/09/16 14:06:05 bmy Exp $
      SUBROUTINE SOILTYPE

!**********************************************************************
!                                                                     *
!  HARVARD ATMOSPHERIC CHEMISTRY MODELING GROUP                       *
!  MODULE FOR SOIL NOx EMISSIONS                                      *
!  by Yuhang Wang, Gerry Gardner and Prof. Daniel Jacob               *
!  (Release V2.1)                                                     *
!                                                                     *
!  Contact person: Bob Yantosca (bmy@io.harvard.edu)                  *
!                                                                     *
!**********************************************************************
! Be sure to force double precision with the DBLE function            *
! and the "D" exponent, wherever necessary (bmy, 10/6/99)             *
!**********************************************************************

      ! References to F90 modules (bmy, 2/11/03)
      USE TIME_MOD, ONLY : GET_MONTH, GET_DAY_OF_YEAR

      USE CMN_SIZE_MOD
      USE COMMSOIL_MOD

      IMPLICIT NONE
#     include "define.h"

!**********************************************************************
!  SOILTYPE DETERMINES WHETHER SOIL IS DRY OR WET                     *
!  UPDATED DAILY.                                                     *
!**********************************************************************
! SOILPREP  = Two month observed precip (mm/day/box                   *
!             (divide by # of days in month))                         *
! JENDDAY   = Julian ending day of previous month                     *
! WETSOIL   = Criteria for wet soil mm                                *
! LENGTHDAY = Number of days for pulse                                *
! MONTHDAY  = Day of the month                                        *
! NCURRENT  = Number of days in current  month                        *
! NPREV     = Number of days in previous month                        *
! JDAY      = Julian day                                              *
! MONTH     = Month number                                            *
! RAIN      = Total rain                                              *
! NPULSE    = Number of types of pulsing                              *
! NLAND     = Total number of land boxes                              *
! SOILPULS  = Tracking of wet/dry & three types of pulsing (Y&L, 94)  *
!**********************************************************************
!

      ! Now make JDAY, MONTH local variables
      INTEGER :: JDAY, MONTH

      INTEGER LENGTHDAY,JDAYSAVE,M,K,MONTHDAY,NCURRENT,NPREV
      
      REAL*8  WETSOIL,RAIN

      REAL*8 JENDDAY(12)
      DATA JENDDAY /0,31,59,90,120,151,181,212,243,273,304,334/
      DATA WETSOIL /10.D0/        !ABOVE 10 MM FOR TWO WEEKS
      DATA LENGTHDAY /14/
      DATA JDAYSAVE /0/ 
      
      !=================================================================
      ! SOILTYPE begins here
      !=================================================================

      ! Get month and day of year
      MONTH = GET_MONTH()
      JDAY  = GET_DAY_OF_YEAR()

      ! If it's a new day...
      IF (JDAYSAVE.NE.JDAY) THEN
         JDAYSAVE=JDAY
         MONTHDAY=JDAY-JENDDAY(MONTH)
         NCURRENT=MIN0(LENGTHDAY,MONTHDAY)
         NPREV=MAX0(0,LENGTHDAY-NCURRENT)
         
         DO M=1,NLAND
! For each land grid-box
	    RAIN=SOILPREP(1,M)*DBLE(NPREV)+SOILPREP(2,M)*
     *           DBLE(NCURRENT)
	    IF (RAIN.GT.WETSOIL) THEN
! WET
               SOILPULS(1,M)=-1.D0
               DO K=1,NPULSE
                  SOILPULS(1+K,M)=0.D0
               END DO
	    ELSE
! DRY
               SOILPULS(1,M)=1.D0
	    END IF
         END DO
      END IF

      RETURN
      END
