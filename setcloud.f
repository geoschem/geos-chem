! $Id: setcloud.f,v 1.1 2003/06/30 20:26:01 bmy Exp $
      SUBROUTINE SETCLOUD(CLOUDS,
     C                    GMU,CLOUDREF,IJLOOP,ALT,LDEBUG,
     C                    PSURF,SIG, IJWINDOW)
!
! NOTE: SETCLOUD is only used for SLOW-J photolysis.
! 

      ! References to F90 modules
      USE ERROR_MOD, ONLY : GEOS_CHEM_STOP

      IMPLICIT NONE

C The values in CLOUDREF are those from AVGF. 1=surface albedo,
C 2-5=fractions of flux through top of atmosphere that are reflected
C by clouds at 2=800mb, 3=500mb, 4=200mb, 5=100mb.

C The values CLOUD(L) are cloud reflectivities at the bottom of standard
C level L (for L=1, it is the surface albedo). So, for L>1, CLOUD(L) is
C at Z=(2*L-3)km.
!
! Remove PTOP from the arg list.  PTOP is now a parameter in 
! "CMN_SIZE". (bmy, 2/10/00)

#     include "CMN_SIZE"
#     include "comode.h"
#     include "CMN_CLD"
      INTEGER IJLOOP,IK
      INTEGER L,M,LL,IALTL, IJWINDOW
      REAL*8 CLOUDS(MAXIJ,11)
      REAL*8 SIG(LGLOB)
      REAL*8 CLOUDP(3),CLOUDREF(5),PSURF
      REAL*8 GMU,ALT(MAXIJ,NPVERT)
      REAL*8 PLEV1,ZLEV1,PLEV2,ZLEV2,CLOUDSIG,SUMALB,WGT,ZCLOUD
      LOGICAL LDEBUG

      DATA CLOUDP /800.0,500.0,200.0/


      DO L=1,11
            CLOUDS(IJLOOP,L)=0.0
      ENDDO

      IF (GMU.LE.0.0) RETURN

C Find the standard level corresponding to 800mb,500mb,200mb.
C (Note that PSURF is really P(surface)-P(top))
      DO 100 M=1,3
         CLOUDSIG = (CLOUDP(M)-PTOP)/PSURF
         LL = 0
         DO L=2,NVERT
            IF (SIG(L-1).GE.CLOUDSIG .AND. SIG(L).LT.CLOUDSIG) LL=L
         ENDDO
         IF (CLOUDSIG.GE.SIG(1)) LL=1
C The cloud layer is between ctm layers LL-1 and LL.
         PLEV1 = SIG(LL)*PSURF+PTOP
         ZLEV1 = ALT(IJWINDOW,LL)
         IF (LL.GT.1) THEN
            PLEV2 = SIG(LL-1)*PSURF+PTOP
            ZLEV2 = ALT(IJWINDOW,LL-1)
         ENDIF
C Interpolate cloud height using the ctm layer heights already
C calculated. Do the interpolation based on log-P.
         IF (LL.GT.1) THEN
            WGT = LOG(CLOUDP(M)/PLEV2) / LOG(PLEV1/PLEV2)
            ZCLOUD = ZLEV2 + WGT * (ZLEV1-ZLEV2)
         ELSE
            ZCLOUD = ZLEV1
         ENDIF
         IALTL = NINT(ZCLOUD/2.E5+1.5)
         IALTL = MAX(IALTL,2)

         IF(IALTL .GT.11) THEN
          WRITE(*,*) 'IALTL greater than 11, IALTL = ', IALTL
          print*, 'ZLEV1, ZLEV2', ZLEV1, ZLEV2
          print*, 'PLEV1, PLEV2', PLEV1, PLEV2
          CALL GEOS_CHEM_STOP
         ENDIF

         CLOUDS(IJLOOP,IALTL) = CLOUDREF(M+1)
      IF (LDEBUG) THEN
        WRITE(*,*) 'CLOUDS:',M,CLOUDREF(M+1),IALTL, 
     C         ZCLOUD/1.e5,PLEV1,ZLEV1/1.e5,PLEV2,ZLEV2/1.e5
      ENDIF
         IF (CLOUDS(IJLOOP,IALTL).LT.0.0) CLOUDS(IJLOOP,IALTL)=0.0
 100  CONTINUE

C surface albedo
      CLOUDS(IJLOOP,1) = CLOUDREF(1)

C make sure we don't have too much reflection
      SUMALB = 0.0
      DO L=2,11
            SUMALB = SUMALB + CLOUDS(IJLOOP,L)
      ENDDO

      IF (LDEBUG) WRITE (*,*) 'SUMALB=',SUMALB !***LWH***

      IF (SUMALB.GT.1.0) THEN
            WRITE(*,*) 'in SETCLOUD: too much reflection: CLOUDS='
            WRITE(*,*) (CLOUDS(IJLOOP,IK),IK=1,11)
            DO L=2,11
                  CLOUDS(IJLOOP,L)=CLOUDS(IJLOOP,L)/SUMALB
            ENDDO
      ENDIF

      RETURN
      END
