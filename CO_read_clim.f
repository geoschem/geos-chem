      SUBROUTINE CO_read_clim()
C-----------------------------------------------------------------------
C  Routine to input O3 reference profiles from Jennifer's
C  2d stratospheric climatologies (monthly) and 3d tropopsheric
C  climatologies (seasonal).
C-----------------------------------------------------------------------
C
      IMPLICIT NONE
#     include "CMN_SIZE"
#     include "CMN"
#     include "CMN_SETUP"
#     include "CMN_CO"
C
      REAL*8 aread(NNMONTHS+2)
      REAL*8 tropO3(IIPAR,JJPAR,NLAYO3)
      INTEGER II,JJ,KK,I,J,L
C      DATA ptropclim /1000.,900.,800.,700.,500.,300.,200.,150.,100./
C
      print*,'Read in jals O3 climatologies for total O3 column'
      print*,'from SR CO_read_clim.'
C*********************
C A) Stratosphere
C*********************
C
C Calculate Z* and corresponding pressure.
C zstarstratclim (cm)
C pstratclim(JJ) (mb)
C
       DO JJ=1,NNALTITUDES
         zstarstratclim(JJ)=2.*real(JJ)-2.
         pstratclim(JJ)=1013./(10**(zstarstratclim(JJ)/16.))
         zstarstratclim(JJ)=zstarstratclim(JJ)*1.D5
       ENDDO
C
      OPEN(UNIT=NO3READ,FILE='O3.profiles.year',STATUS='OLD')
C
      READ(NO3READ,88)
 88   FORMAT(//)
C
C stratclimO3 (ppmv)
C
      DO KK=1,NNLATITUDES
       DO JJ=1,NNALTITUDES
        READ(NO3READ,*)aread
         DO II=3,NNMONTHS+2
           stratclimO3(KK,JJ,II-2)=aread(II)*1.D-6
         ENDDO
       ENDDO
      ENDDO
C
      CLOSE(NO3READ)
C
C*********************
C Troposphere
C*********************
C
C tropO3 (ppbv)
C ptropclim (mb)
C
      OPEN(UNIT=NO3READ,FILE='o3.win',STATUS='OLD')
C
      DO I=1,IIPAR
      DO J=1,JJPAR
      DO L=1,NLAYO3
           READ(NO3READ,*) tropO3(I,J,L)
      ENDDO
      ENDDO
      ENDDO
C
      tropclimO3(:,:,:,1)=tropO3(:,:,:)
C
      CLOSE(NO3READ)
C
      OPEN(UNIT=NO3READ,FILE='o3.spr',STATUS='OLD')
C
      DO I=1,IIPAR
      DO J=1,JJPAR
      DO L=1,NLAYO3
           READ(NO3READ,*) tropO3(I,J,L)
      ENDDO
      ENDDO
      ENDDO
C
      tropclimO3(:,:,:,2)=tropO3(:,:,:)
C
      CLOSE(NO3READ)
C
      OPEN(UNIT=NO3READ,FILE='o3.sum',STATUS='OLD')
C
      DO I=1,IIPAR
      DO J=1,JJPAR
      DO L=1,NLAYO3
           READ(NO3READ,*) tropO3(I,J,L)
      ENDDO
      ENDDO
      ENDDO
C
      tropclimO3(:,:,:,3)=tropO3(:,:,:)
C
      CLOSE(NO3READ)
C
      OPEN(UNIT=NO3READ,FILE='o3.aut',STATUS='OLD')
C
      DO I=1,IIPAR
      DO J=1,JJPAR
      DO L=1,NLAYO3
           READ(NO3READ,*) tropO3(I,J,L)
      ENDDO
      ENDDO
      ENDDO
C
      tropclimO3(:,:,:,4)=tropO3(:,:,:)
C
      CLOSE(NO3READ)
C
      RETURN
      END

