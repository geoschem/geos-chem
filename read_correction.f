C**********************************************************************
      SUBROUTINE read_correction
C**********************************************************************
      IMPLICIT NONE
#     include "CMN_OH"
      INTEGER I,J,K
C**********************************************************************
C Created by Bryan Duncan.
C**********************************************************************
c This SR reads in correction factors to remove the problem with high
c background aerosol (tau=0.1) in the OH parameterization.
c This file contains ratios of OH calculated assuming a uniform
c tau=0.01 over the globe and assuming tau=0.1. OH calculated by CMS, June, 2001.
c Ratio created by bnd, June, 2001.
c*************************************************************************
c      INTEGER NCMSALTS2,NCMSLATS2,NSEAS
c      PARAMETER(NSEAS=4,NCMSALTS2=9,NCMSLATS2=24)
c      REAL*4 CMSALTS2(NCMSALTS2),CMSLATS2(NCMSLATS2)
c      REAL*4 correction(NSEAS,NCMSLATS2,NCMSALTS2)
c      DATA CMSALTS2/1000.,900.,800.,700.,500.,300.,200.,150.,100./
c      DATA CMSLATS2/89.,84.,76.,68.,60.,52.,44.,36.,28.,20.,12.,4.,
c     * -4.,-12.,-20.,-28.,-36.,-44.,-52.,-60.,-68.,-76.,-84.,-89./
c*************************************************************************
C List of Variables & Arrays
C
C correction = correction factor 
C
C NCMSALTS2 = number of altitude levels.
C
C NCMSLATS2 = number of latitude bands.
C
C NSEAS    = number of seasons.
C
C NGENERICFILE   = filenumber specified by user in CMN_OH.
C**********************************************************************
C      print*,'Reading in OH correction factors in SR read_correction.'
C              CALL FLUSH( 6 )
C
      OPEN(UNIT=NGENERICFILE,FILE='CorrectRatio',STATUS='OLD')
C
      correction(:,:,:)=1.
C
      READ(NGENERICFILE,2)
 2    FORMAT(//////////////)
C
      DO I=1,NSEAS
C
      READ(NGENERICFILE,3)
 3    FORMAT()
C
      DO J=1,NCMSLATS2
        READ(NGENERICFILE,4)(correction(I,J,K),K=1,NCMSALTS2)
      ENDDO
C
      ENDDO
C
 4    FORMAT(4x,9(f5.3,2x)) 
C
      CLOSE(NGENERICFILE)
C
      DO I=1,NSEAS
      DO J=1,NCMSLATS2
      DO K=1,NCMSALTS2
      IF(correction(I,J,K).LT.0..OR.correction(I,J,K).GT.2.)THEN
           PRINT*,'Correction factor < 0!'
           PRINT*,I,J,K,correction(I,J,K)
           PRINT*,'Stopped in SR READ_CORRECTION.'
           STOP
        ENDIF
      ENDDO
      ENDDO
      ENDDO
C
C    print*,'Leaving SR read_correction.'
C              CALL FLUSH( 6 )
C**********************************************************************
      RETURN
      END
