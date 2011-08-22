C $Id: soiltype.f,v 1.1 2009/09/16 14:06:05 bmy Exp $
      SUBROUTINE SOILTYPE

C**********************************************************************
C                                                                     *
C  HARVARD ATMOSPHERIC CHEMISTRY MODELING GROUP                       *
C  MODULE FOR SOIL NOx EMISSIONS                                      *
C  by Yuhang Wang, Gerry Gardner and Prof. Daniel Jacob               *
C  (Release V2.1)                                                     *
C                                                                     *
C  Contact person: Bob Yantosca (bmy@io.harvard.edu)                  *
C                                                                     *
C**********************************************************************
C Be sure to force double precision with the DBLE function            *
C and the "D" exponent, wherever necessary (bmy, 10/6/99)             *
C**********************************************************************

      ! References to F90 modules (bmy, 2/11/03)
      USE TIME_MOD, ONLY : GET_MONTH, GET_DAY_OF_YEAR

      USE CMN_SIZE_MOD
      USE COMMSOIL_MOD

      IMPLICIT NONE
#     include "define.h"

C**********************************************************************
C  SOILTYPE DETERMINES WHETHER SOIL IS DRY OR WET                     *
C  UPDATED DAILY.                                                     *
C**********************************************************************
C SOILPREP  = Two month observed precip (mm/day/box                   *
C             (divide by # of days in month))                         *
C JENDDAY   = Julian ending day of previous month                     *
C WETSOIL   = Criteria for wet soil mm                                *
C LENGTHDAY = Number of days for pulse                                *
C MONTHDAY  = Day of the month                                        *
C NCURRENT  = Number of days in current  month                        *
C NPREV     = Number of days in previous month                        *
C JDAY      = Julian day                                              *
C MONTH     = Month number                                            *
C RAIN      = Total rain                                              *
C NPULSE    = Number of types of pulsing                              *
C NLAND     = Total number of land boxes                              *
C SOILPULS  = Tracking of wet/dry & three types of pulsing (Y&L, 94)  *
C**********************************************************************
C

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
C For each land grid-box
	    RAIN=SOILPREP(1,M)*DBLE(NPREV)+SOILPREP(2,M)*
     *           DBLE(NCURRENT)
	    IF (RAIN.GT.WETSOIL) THEN
C WET
               SOILPULS(1,M)=-1.D0
               DO K=1,NPULSE
                  SOILPULS(1+K,M)=0.D0
               END DO
	    ELSE
C DRY
               SOILPULS(1,M)=1.D0
	    END IF
         END DO
      END IF

      RETURN
      END
