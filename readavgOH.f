C $Id: readavgOH.f,v 1.1 2003/06/30 20:26:05 bmy Exp $
C**********************************************************************
      SUBROUTINE readavgOH
C**********************************************************************
      IMPLICIT NONE
#     include "CMN_OH"
      INTEGER I,J,K
C**********************************************************************
C Created by Bryan Duncan.
C**********************************************************************
C In subdomains where OH is very low or negligible due to low 
C  sunlight (e.g., high latitudes in winter), concentrations of OH are 
C  set to climatological mean values as a function of latitude,
C  altitude and season. This SR reads in the average OH
C  fields from Spivakovsky et al., "Three-dimensional climatological 
C  distribution of tropospheric OH: update and evaluation", accepted
C  to JGR, 1999.
C**********************************************************************
C List of Variables & Arrays
C
C avgOH = array containing the climatological OH values. 
C
C NCMSALTS = number of altitude levels of climatology.
C
C NCMSLATS = number of latitude bands of climatology.
C
C NSEAS    = number of seasons of climatlogy.
C
C NGENERICFILE   = filenumber specified by user in CMN_OH.
C**********************************************************************
C 
      print*,'Reading in OH climatlogy of cms in readavgOH.'
C
      OPEN(UNIT=NGENERICFILE,FILE='avgOH',STATUS='OLD')
C
      READ(NGENERICFILE,2)
 2    FORMAT(//)
C
      DO I=1,NSEAS
C
      READ(NGENERICFILE,3)
 3    FORMAT(////)
C
      DO J=1,NCMSLATS
        READ(NGENERICFILE,4)(avgOH(I,J,K),K=1,NCMSALTS)
      ENDDO
C
      ENDDO
C
 4    FORMAT(5x,7(F6.2))
C
      CLOSE(NGENERICFILE)
C
C**********************************************************************
      RETURN
      END
