! $Id: interpOH.f,v 1.1 2003/06/30 20:26:08 bmy Exp $
C**********************************************************************
      SUBROUTINE INTERPOH
C**********************************************************************
      ! References to F90 modules
      USE ERROR_MOD, ONLY : ERROR_STOP
      IMPLICIT NONE
#     include "CMN_OH"
      INTEGER I,J,K
C**********************************************************************
C Created by Bryan Duncan.     
C**********************************************************************
C In subdomains where OH is very low or negligible due to low
C  sunlight (e.g., high latitudes in winter), concentrations of OH are
C  set to climatological mean values as a function of latitude,
C  altitude and season. This SR picks the appropriate average OH
C  field described in Spivakovsky et al., "Three-dimensional climatological
C  distribution of tropospheric OH: update and evaluation", accepted
C  to JGR, 1999.  The fields are stored in array avgOH and read in
C  in SR readavgOH.
C
C avgOH = array containing the climatological OH values.
C
C NCMSALTS = number of altitude levels of climatology
C
C NCMSLATS = number of latitude bands of climatology
C
C**********************************************************************
C 
      I=OH_SEASON
C
      DO J=1,NCMSLATS
C
       DO K=NCMSALTS,1,-1
C
         IF(OH_LAT.GE.CMSLATS(J)) THEN         
C
           IF(CMSALTS(K).GE.OH_PRESS) THEN
              PARAMOH=avgOH(I,J,K)*1.D5
              GOTO 2
           ENDIF
C
           IF(K.EQ.1) THEN
              PARAMOH=avgOH(I,J,K)*1.D5
              GOTO 2
           ENDIF
C
         ENDIF
C
       ENDDO
C
      ENDDO
C
C Error Check.
      !-----------------------------------------
      ! Prior to 10/15/02:
      !PRINT*,'STOPPED IN SR interpOH!'
      !PRINT*,'Point lies nowhere!'
      !STOP
      !-----------------------------------------
      CALL ERROR_STOP( 'Point lies nowhere', 'interpOH.f' )

C
 2    CONTINUE
C
      RETURN
      END
