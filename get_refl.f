C $Id: get_refl.f,v 1.1 2003/06/30 20:26:07 bmy Exp $
      SUBROUTINE get_refl(OPTD,Rout,SZA)
C*****************************************************************************
C Using optical depth and the cos(SZA), use reflectivity
C   table to find the reflectivity above each point.  
C   To calculate the flux fraction at a point:
C    flux fraction=(flux above)/(flux total from TOA to surface)
C Adapated & modified from chem1d by Bryan Duncan (1/99).
C*****************************************************************************
C The annual mean tropopause is stored in the LPAUSE array 
C (from header file "CMN").  LPAUSE is defined such that: 
C
C Levels            1 <= L <= LPAUSE(I,J) - 1 are tropospheric
C         LPAUSE(I,J) <= L <= LLPAR           are stratospheric
C
C We now use LPAUSE instead of NSKIPL to denote the strat/trop boundary. 
C (bmy, 4/18/00)  
C*****************************************************************************
      IMPLICIT NONE
#     include "CMN_SIZE"  ! Size parameters
#     include "CMN"       ! LPAUSE
#     include "CMN_CO"    ! CO arrays
C     
      INTEGER NMAX,MMAX,JA,KA
      PARAMETER (NMAX=5,MMAX=5)
      REAL*8 YMTMP(MMAX),YNTMP(NMAX)
      INTEGER NREFL,NREFLCOL,ICOL,IINT,JINT,IROW,II,JJ,Q
      INTEGER I,J,K,L
      
      ! Remove duplicate definition of SZA (bmy, 11/15/01)
      REAL*8 Y,DY,RFLCVAL,PI,TORAD,BCOSSZA(IIPAR,JJPAR),
     &     REFINT(3,3),OPTDINT(3),COSTHETAINT(3),
     &     PTRFLCABOVE,PTRFLCBELOW,PTRFLC,
     &     PTOPTDEPTH(IIPAR,JJPAR,LLPAR),
     &     OPTD(LLPAR,IIPAR,JJPAR),SZA(IIPAR,JJPAR),
     &     Rout(IIPAR,JJPAR,LLPAR),
     &     COSTHETA(NREFLCOL),WGHTS(NREFLCOL)
      DATA COSTHETA /0.0469,0.2308,0.5,0.7692,0.9531/
      DATA WGHTS/0.1185,0.2393,0.2844,0.2393,0.1185/
C     
C  To calculate the flux fraction above (Ravga) each point,
C    first sum up the optical depths above. (Optical depths
C    are additive.)  Make the assumption that the point never
C    lies inside a cloud.
C  If the optical depth is greater than 90, then set it =90,
C  which is the maximum optical depth on the chart.
C  Compute cos(SZA) for the point.
C  If the cos(SZA) is beyond the limits of the reflectivity
C  table, set it equal to the limit.
C
      PI=3.14
      TORAD=PI/180.
      BCOSSZA(:,:)=0.	
      PTOPTDEPTH(:,:,:)=0.
C
      DO I=1,IIPAR
      DO J=1,JJPAR
      ! L now goes up to the annual mean tropopause (bmy, 4/18/00)
      !DO L=1,NSKIPL-1
      DO L = 1, LPAUSE(I,J) - 1
C
         PTOPTDEPTH(I,J,L)=SUM(OPTD(L+1:LLPAR,I,J))
         IF(PTOPTDEPTH(I,J,L).GT.90.) PTOPTDEPTH(I,J,L)=90.
C     
      ENDDO
      ENDDO
      ENDDO
C
      BCOSSZA(:,:)=COS(SZA(:,:)*TORAD)
      WHERE(BCOSSZA(:,:).LT.0.0469)BCOSSZA(:,:)=0.0469
      WHERE(BCOSSZA(:,:).GT.0.9531)BCOSSZA(:,:)=0.9531
C
C*****************************************************************************
C Begin looping over IIPAR and JJPAR
C*****************************************************************************
C
C  Now we have both optical depth and cos(SZA), and can use cms'
C  reflectivity table to find the reflectivity of point.  note:
C  need to interpolate in two directions in the reflectivity table.
C
      DO I=1,IIPAR
      DO J=1,JJPAR
      ! L now goes up to the annual mean tropopause (bmy, 4/18/00)
      !DO L=1,NSKIPL-1
      DO L = 1, LPAUSE(I,J) - 1
C
C  If the solar zenith angle is > 85, then set PTRFLC=0.      
C
         IF(SZA(I,J).GT.85.) THEN
            PTRFLC=-999.
            GOTO 200 
         ENDIF
C
C  Calculated reflectivity is stored in PTRFLC.
C
         PTRFLC=0.
C
C  If the solar zenith angle is < 85, then update Rcount.
C
         Rcount(I,J,L)=Rcount(I,J,L)+1
C
C  Loop across columns, looking for the correct cos(SZA) values
C  to interpolate between.  note:  the columns correspond to the
C  quadrature points.
C     
         DO ICOL=1,NREFLCOL-1
            IF ((BCOSSZA(I,J).GE.COSTHETA(ICOL)).AND.
     &           (BCOSSZA(I,J).LE.COSTHETA(ICOL
     &           +1))) THEN 
               JINT=ICOL
               GOTO 230
            ENDIF
         ENDDO
C     
 230     IF((BCOSSZA(I,J)-COSTHETA(NREFLCOL-1).LE.(COSTHETA(NREFLCOL)
     &        -BCOSSZA(I,J))).AND.(JINT.GT.1)) THEN
            JINT=JINT-1
         ENDIF
         IF(JINT.GT.3) JINT=3
C     
C Look down columns for correct optical depth to interpolate between.
C
         DO IROW=1,NREFL-1
            IF((PTOPTDEPTH(I,J,L).GE.OPTDEPTH(IROW)).AND.
     &           (PTOPTDEPTH(I,J,L).LE.OPTDEPTH
     &           (IROW+1))) THEN
               IINT=IROW
               GOTO 330
            ENDIF 
         ENDDO
C
 330     IF(((PTOPTDEPTH(I,J,L)-OPTDEPTH(NREFL-1)).LE.
     &        (OPTDEPTH(NREFL)-
     &        PTOPTDEPTH(I,J,L))).AND.(IINT.GT.1)) THEN 
            IINT=IINT-1
         ENDIF
         IF(IINT.GT.(NREFL-2)) THEN
            IINT=NREFL-2
         ENDIF 
C     
         DO II=1,3
            OPTDINT(II)=OPTDEPTH(IINT+II-1)
            COSTHETAINT(II)=COSTHETA(JINT+II-1)
            DO JJ=1,3
               REFINT(II,JJ)=RFLC(IINT+II-1,JINT+JJ-1)
            ENDDO
         ENDDO
C     
C ****************************************************************************
C Interpolation
C  Calculates an interpolated function value y
C  and an accuracy indication dy.
C
         DO 12 JA=1,3
            DO 11 KA=1,3
               YNTMP(KA)=REFINT(JA,KA)
 11         CONTINUE
            
            CALL POLINT(COSTHETAINT,YNTMP,3,BCOSSZA(I,J),YMTMP(JA),DY)
            
 12      CONTINUE

         CALL POLINT(OPTDINT,YMTMP,3,PTOPTDEPTH(I,J,L),Y,DY)
C
C ****************************************************************************
C
         RFLCVAL=Y
C
C Now use the quadrature weights to get the final reflectivity.
C
         DO 120 Q=1,5
 120     PTRFLC=PTRFLC+WGHTS(Q)*RFLCVAL
C
 200     CONTINUE
C*****************************************************************************
C Reflectivity above is stored in RAVGA.
C  The total reflectivity from the surface to TOA is stored in
C  ROUT(I,J,1).
         Rout(I,J,L)=PTRFLC
C*****************************************************************************
      ENDDO
      ENDDO
      ENDDO
C*****************************************************************************
C
      RETURN
      END
