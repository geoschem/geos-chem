C $Id: avgrefl.f,v 1.1 2003/06/30 20:26:03 bmy Exp $
      SUBROUTINE avgrefl(FIRSTCHEM,OPTD,NTDT,NHMSb,NSEC,XLON,NMIN)
C*****************************************************************************
      IMPLICIT NONE
#     include "CMN_SIZE"  ! Size parameters
#     include "CMN"       ! STT, LPAUSE
#     include "CMN_CO"    ! CO arrays
      LOGICAL FIRSTCHEM
      INTEGER I,J,K,L,M,NTDT,NTIMES,MNDT,MOUTPUTON,NMIN
      INTEGER IJLOOP,IREF,JREF,LHR0,NHMSb,NSEC,NFNUM
      REAL*8  OPTD(LLPAR,IIPAR,JJPAR),Rout(IIPAR,JJPAR,LLPAR),
     &     SZA(IIPAR,JJPAR),templat,templon,temp,
     &     SUNCOS(MAXIJ),XLON(IGLOB)
     
cbnd
      NFNUM = 65
cbnd

C*****************************************************************************

C  NEED TO SEE IF DAYLIGHT AVERAGE VS 6-HR AVG AROUND
C  NOON PRODUCES DIFFERENT RESULTS!!!!!!! 

C*****************************************************************************
C
C Designed to convert optical depths to reflectivities
C   and calculate a daylight average which will be used
C   for the CO/OH parameterization option (see SR chemco).
C   Cannot average optical depths, but can average reflectivities.
C Created by Bryan Duncan (1/99).
C
C The "daylight average" can be misleading since, for OH,
C   the optical depths around the solar apex (i.e., noon)
C   are more important in determining its concentration than
C   around dawn or sunset.  Therefore, areas with morning
C   fog (such as marine environments) or areas affected by
C   tropical convection may have a "skewed" daylight average
C   reflectivities. (see Rossow and Schiffer, ISCCP Cloud Data 
C   Products, p.15)
C
C NCHEM is chemical time step (i.e., 24 hours).
C NDT   is meteorological data time step (i.e., 6 hours).
C NTDT  is dynamical time step (i.e., 1/2 hours).
C
C*****************************************************************************
C The annual mean tropopause is stored in the LPAUSE array 
C (from header file "CMN").  LPAUSE is defined such that: 
C
C Levels            1 <= L <= LPAUSE(I,J) - 1 are tropospheric
C         LPAUSE(I,J) <= L <= LLPAR           are stratospheric
C
C (bmy, 4/18/00)  
C*****************************************************************************
C
      MNDT=NTDT/60
      NTIMES=NCHEM/MNDT

      IF(NMIN.LE.NCHEM) NTIMES=NTIMES+1
C     
      IF (NMIN.EQ.0) THEN
         MNEW=0
         MNCOUNT=0
      ENDIF
C
C If a new 24-hr period, set Ravg = 0.
C
      IF (MNEW.EQ.0) THEN
         Rcount(:,:,:)=0.
         Ravg(:,:,:)=0.
         Ravga(:,:,:)=0.
         Ravgb(:,:,:)=0.
         MNEW=1
         MNCOUNT=0
      ENDIF
C*****************************************************************************
C Read in reflectivity table (Srefl.table)
C Table data used to convert user's optical depths to
C   flux fractions - see get_fluxfraction.f.
C*****************************************************************************
      IF(FIRSTCHEM) THEN
C
         OPEN (UNIT=NFNUM,FILE='Srefl.table',STATUS='old')
         REWIND(NFNUM)
C     
         READ(NFNUM,9)
 9       FORMAT(////)
         
         DO 90 M=1,NREFL
 90      READ(NFNUM,100)OPTDEPTH(M),(RFLC(M,J),J=1,5)
C     
 100     FORMAT(3X,F6.3,5(F8.3))
C     
         CLOSE(NFNUM)
            
      ENDIF
C*****************************************************************************
C
C Calculate solar zenith angle (SZA).
C
      CALL COSSZA(JDAY,NHMSb,NSEC,RLAT,XLON,SUNCOS)
C
C Convert SUNCOS (radians) to SZA (degrees).
C
      DO J=1,JJPAR
      DO I=1,IIPAR
         SZA(I,J)=ACOS(SUNCOS((J-1)*IIPAR+I))*180./3.14
      ENDDO
      ENDDO
C       
C SR get_refl converts optical depths to reflectivities
C   and then to flux fractions for a given solar zenith angle.
C
      Rout(:,:,:)=0.
      CALL get_refl(OPTD,Rout,SZA)
C
C Sum reflectivities in Ravga() over daylight hours.
C If Rout = -999 then dark, so don't include in sum.
C
      WHERE(Rout(:,:,:).GE.0.)Ravg(:,:,:)=Ravg(:,:,:)+
     &     Rout(:,:,:)
C
C*****************************************************************************
C Keep track to see if at end of NCHEM time step.
C*****************************************************************************
C
      MNCOUNT=MNCOUNT+1
C
      IF (MNCOUNT.EQ.NTIMES) THEN
C
C Average Ravg over NCHEM time step.
C
         WHERE(Rcount(:,:,:).GT.0)
     &        Ravg(:,:,:)=Ravg(:,:,:)/REAL(Rcount(:,:,:))
C
C If Rcount = 0 and Ravg < 0 then polar night!
C
         WHERE(Rcount(:,:,:).LE.0) Ravg(:,:,:) = -999.
C 
C Convert average reflectivities to flux fractions above
C   and below.
C     
         Ravga(:,:,:)=0.
         Ravgb(:,:,:)=0.
C
C above:
C
         WHERE(Ravg(:,:,:).GT.0.) Ravga(:,:,:)=Ravg(:,:,:)
         WHERE(Ravg(:,:,:).GT.-1..AND.Ravg(:,:,:).LT.0.)
     &           Ravga(:,:,:)=0.
         WHERE(Ravg(:,:,:).LT.-1.) Ravga(:,:,:)=-999.
C
C below:
C
         WHERE(Ravga(:,:,:).LT.0.) Ravgb(:,:,:)=-999.

         ! Now L goes up to the annual mean tropopause (bmy, 4/17/00)       
         DO J = 1, JJPAR
         DO I = 1, IIPAR
         DO L = 1, LPAUSE(I,J) - 1
            IF(Ravga(I,J,L).GT.0.) Ravgb(I,J,L)=Ravg(I,J,1)-
     &           Ravga(I,J,L)
         ENDDO
         ENDDO
         ENDDO

         WHERE(Ravgb(:,:,:).GT.-1.AND.Ravgb(:,:,:).LT.0.)
     &        Ravgb(:,:,:)=0.
C
         MNEW=0
C*****************************************************************************
         MOUTPUTON=0
         IF(MOUTPUTON.EQ.1) THEN
            open(NFNUM,file='bryanoutput',status='unknown')
            rewind(NFNUM)
            
            DO I=1,IIPAR
            DO J=1,JJPAR
            ! Now L goes up to the annual mean tropopause (bmy, 4/17/00)
            DO L = 1, LPAUSE(I,J) - 1
               write(NFNUM,*)Ravg(I,J,L)
            ENDDO
            ENDDO
            ENDDO

            DO I=1,IIPAR
            DO J=1,JJPAR
            ! Now L goes up to the annual mean tropopause (bmy, 4/17/00)
            DO L = 1, LPAUSE(I,J) - 1
               write(NFNUM,*)Ravga(I,J,L)
            ENDDO
            ENDDO
            ENDDO

            DO I=1,IIPAR
            DO J=1,JJPAR
            ! Now L goes up to the annual mean tropopause (bmy, 4/17/00)
            DO L = 1, LPAUSE(I,J) - 1
               write(NFNUM,*)Ravgb(I,J,L)
            ENDDO
            ENDDO
            ENDDO

            DO I=1,IIPAR
            DO J=1,JJPAR
            ! Now L goes up to the annual mean tropopause (bmy, 4/17/00)
            DO L = 1, LPAUSE(I,J) - 1 
               if(Ravg(I,J,L).gt.0.) THen
                  temp=(Ravgb(I,J,L)+Ravga(I,J,L))/Ravg(I,J,1)
               else
                  temp=-999
               endif
               write(NFNUM,*)temp
            ENDDO
            ENDDO
            ENDDO

            ! Now L goes up to the annual mean tropopause (bmy, 4/17/00)
            DO I=1,IIPAR
            DO J=1,JJPAR
            DO L = 1, LPAUSE(I,J) - 1
               write(NFNUM,*)OPTD(L,I,J)
            ENDDO
            ENDDO
            ENDDO

            DO I=1,IIPAR
            DO J=1,JJPAR
            DO L = 1, LPAUSE(I,J)
               write(NFNUM,*)Rcount(I,J,L)
            ENDDO
            ENDDO
            ENDDO
               
            CLOSE(NFNUM)
            STOP
         ENDIF
C*****************************************************************************
C
C*****************************************************************************
      ENDIF
C*****************************************************************************
C
      RETURN
      END
