C $Id: scattr.f,v 1.1 2003/06/30 20:26:03 bmy Exp $
      SUBROUTINE SCATTR(IJLOOP, GMU, CLOUDS, ZZHT)

      IMPLICIT NONE

C*
C**********************************************************************
C                                                                     *
C  HARVARD TROPOSPHERIC CHEMISTRY MODULE FOR 3-D APPLICATIONS         *
C  by Larry Horowitz, Jinyou Liang, Gerry Gardner, Prof. Daniel Jacob *
C  of HARVARD UNIVERSITY    (Release V2.0)                            *
C**********************************************************************
C*
C**********************************************************************
C*
C*    Compute the effective transmission of the atmosphere, given TAU(),
C*    the optical depth, PITAU(), the single scattering albedo, and 
C*    RFLECT, the Lambert reflectivity of ground; assume TAU(1)=0.
C*    Atmosphere stops at I=NTT; usually NC+1.
C*              
C**********************************************************************
C*                     
#     include "CMN_SIZE"
#     include "comode.h"
      INTEGER IJLOOP
      REAL*8 CLOUDS(MAXIJ,11)
      REAL*8 GMU, ZZHT
      REAL*8 B(3,3)
      INTEGER III,ICLD,NTTDO,K,NTT1,J,II,NWT,NMU,IMU

      REAL*8 WTCLD,DTEMP,CTEMP,ATEMP,BTEMP,DTC,FACTOR,XMEAN,DTA,DTAU,
     1  QTEMP,WTING,SCALEH,TRANS0,TTFBOT,AIRMAS,U0,GMU0
C*                  
C*    TTA(J) = EXP(-TAU(J)/GMU), where GMU is corrected for the effects
C*             of curvature of the atmosphere with the function AIRMAS
C*    TTU(II)  -- the cosines of the quadrature angles
C*           (6-stream calculation)
C*    TTU2(II) -- TTU(II)**2
C*    TTP0(II) -- the light scattered from the incident beam into the 
C*                quadrature angle II.  For isotropic scattering by
C*                aerosols, TTP0(II) = 1.  For Rayleigh-phase scattering
C*                TTPO(II) is calculated from the P0 matrix as given by
C*                CHANDRASEKHAR and PRATHER(1974), assuming that the
C*                incident and scattered lights are natural, i.e.,
C*                I(TOTAL) = 1/2 (I= + IT), which reduces the PO matrix
C*                to a single number.  Since we look at the radiation
C*                averaged over all azimuthal angles, the P1 and P2
C*                matrices are 0 .
C*    TTP(II,III) -- the light scattered from the quadrature angle III
C*                   into the quadrature angle II.  Calculated in the
C*                   same way as TTP0(II).
C*    TTW(II)     -- the weights assigned to each quadrature angle
C*                   (cf. CHANDRASEKHAR, CHAP.2)
C*
      REAL*8 TTA(66),TTC(3,66),TTD(3,3,66),TTBOT(10),TTEX(66)
      REAL*8 TTQ(3),TTU(3),TTU2(3),TTW(3),TTP(3,3),TTP0(3),TTF(3,66)
C*                 
      LOGICAL LREDO,LDIFUS
C*                    
      REAL*8 TAU,PITAU,FLTAU,RFLECT,AIRM
      INTEGER NTT
      COMMON /CCSCAT/ TAU(66),PITAU(66),FLTAU(66),RFLECT,NTT,LDIFUS

      !=================================================================
      ! /CCSCAT/ needs to be declared THREADPRIVATE for the OpenMP
      ! parallelization for all platforms (bmy, 3/23/03)
      !=================================================================
!$OMP THREADPRIVATE( /CCSCAT/ )
C*                     
      DATA  TTP/0.3099042361D0,0.2861275609D0,0.2322916667D0,
     1          0.4578040972D0,0.4479166666D0,0.4255292360D0,
     2          0.2322916667D0,0.2659557725D0,0.3421790972D0/
      DATA  TTU/0.1127016654D0,0.5000000000D0,0.8872983346D0/
      DATA  TTW/0.2777777778D0,0.4444444444D0,0.2777777778D0/
      DATA TTU2/0.0127016654D0,0.2500000000D0,0.7872983346D0/
      DATA NWT/3/
C*                    
C*    If NWT = 4; then  
C*      DATA TTU/0.0694318442,0.3300094782,0.6699905218,0.9305681558/
C*      DATA TTW/0.1739274226,0.3260725774,0.3260725774,0.1739274226/
C*          
C*    ZZHT, scale height of the atmosphere; RAD, radius of the earth. 
C*    TAU(NTT) is the optical depth at the bottom of the atmosphere.
C*   
      SCALEH = ZZHT/6345.65D5
      IF (TAU(NTT) .GT. 1.0D-4) GO TO 110
C*  
C*    Optically thin limit----no allowance for clouds

      TRANS0 = 1.0D0
      IF (LDIFUS) TRANS0 = 1.0D0 + GMU*2.0D0*RFLECT

      DO II=1,NTT
         FLTAU(II) = TRANS0
      ENDDO
      GO TO 460
C* 
C*    Standard transmission problem
  110 NTT1 = NTT - 1
      DO II=1,10
         TTBOT(II) = 0.0D0
      ENDDO
C*
C*    NWT = 3;  [data statement at the beginning of subroutine]
      DO II=1,66
         DO J=1,NWT
            TTF(J,II) = 0.0D0
         ENDDO
         FLTAU(II) = 0.0D0
      ENDDO
C*
C*    Calculate the attenuated incident beam EXP(-TAU/U0), averaged 
C*    over the appropriate angles.  The following section corrects the 
C*    optical path through each layer for the earth's curvature.  The
C*    correction (function AIRMAS) is very small, except for large 
C*    zenith angles.  Without the correction, AIRMAS = 1/GMU
C*                      
      WTING = 1.0D0
      TTFBOT = 0.0D0
C*
      NMU = 1
      GMU0= 0.D0
C*                 
      DO 200 IMU=1,NMU
         IF (GMU .LT. GMU0) GO TO 210
         U0 = GMU

         AIRM = AIRMAS(U0,SCALEH)
C*                       
C*    Rayleigh phase source function, calculated as indicated in a
C*    comment above.
         DO J=1,NWT
            TTP0(J) = 0.375D0*(3.0D0 - TTU2(J) - U0*U0 +
     x                              3.0D0*TTU2(J)*U0*U0)
         ENDDO
C*            
C*    TTA(1) is EXP(optical depth) at the top of the atmosphere = 1
C*    TTA(II) is then calculated going downward. 
C*    TTA(61) is the value at the bottom of the atmosphere.
C     Instead of using our complicated ATEMP/BTEMP/...
C     scheme to calculate ATEMP, just use the straight-
C     forward EXP(-TAU*AIRMASS), should now agree w/ rjs
C*
      DO 23 II=1,NTT
C*    NOTE: TTA(I) should be equal to EXP(-TAU(I)*AIRMASS)
C*    TTF(J,II) is the normalized intensity of the beam scattered from
C*    the incident beam into quadrature angle J.
         TTA(II) = EXP(-TAU(II)*AIRM)
         DO J=1,NWT
            TTF(J,II) = TTF(J,II) + WTING*TTA(II)*TTP0(J)
         ENDDO
         FLTAU(II) = FLTAU(II) + WTING*TTA(II)
 23   CONTINUE
      DO III=1,10
         TTBOT(III) = TTBOT(III) + WTING*U0*TTA(NTT+1-III)
      ENDDO

  200 CONTINUE
  210 CONTINUE
C*                    
C*    The optical path has been correctly calculated, giving us the
C*    value of the attenuated direct beam.  
C*    FLTAU -- averaged direct beam 
C*    TTF   -- averaged Rayleigh-phase angle-weighted source
      IF (.NOT. LDIFUS) GO TO 460
      LREDO = .FALSE.
C*    
C*    Determine if multiple level clouds are to be included:
      WTCLD = 1.D0
      DO ICLD=2,11
         WTCLD = WTCLD - CLOUDS(IJLOOP,ICLD)
      ENDDO
      DO II=1,66
         TTEX(II) = 0.0D0
      ENDDO
      DO 440 ICLD=1,11
         IF (ICLD .GT. 1) WTCLD = CLOUDS(IJLOOP,ICLD)
* for first layer, if cumulative cloud reflectivity >99.9%
         IF (WTCLD .LT. 0.001D0) GO TO 440
         IF(ICLD.LT.11) THEN
            TTFBOT = TTBOT(ICLD)
         ELSE
            TTFBOT = 0.D0
         ENDIF
         IF (ICLD .GT. 1) RFLECT = 1.00D0
         NTTDO = NTT + 1 - ICLD
C*    
C*    Include the diffusely scattered field:
C*    TTP(J,K) = P(U(J),U(K)) : the scattering matrix 
C*             = WT(K)        : isotropic
C*             = (3/8)(3 - UJ2 - UK2 + 3*UJ2*UK2)*WT(K) : Rayleigh phase
C*    TTP0(J) = P(U(J),U0) : incident source scattering
C*            = 1          : isotropic
         NTT1 = NTTDO - 1
C*    
C     Second order upper boundary condition
         DTAU = TAU(2)-TAU(1)
         IF(LREDO) GO TO 270
         DO 250 J=1,NWT
C*    
C*    TTU(J) are the cosines of angles used in the quadrature
            ATEMP = TTU(J)/DTAU
            BTEMP = 0.50D0/ATEMP
C*    
C*    PITAU(J) is the contribution to TAU(J) from layer J
C*    TTQ is DJ(O,TTU) = (.5*DTAU/TTU)*J(0,TTU),
C*       where J(0,TTU) is the source function from the incident flux,
C*           (1/4)*TTP0(J)*TTA(1)=.25*TTF(J,1)
C*    Normally, PITAU(1)=1.
C*    
C*    Set up 3 2nd-order O.D.E.'s with coefficients B(J,K) for each 
C*    altitude layer (layer 1 is the highest); solve with some matrix
C*    inversion method.
C*        
            TTQ(J) = BTEMP*0.25D0*PITAU(1)*TTF(J,1)
            DO K=1,NWT
               B(J,K) = -PITAU(1)*TTP(J,K)*BTEMP
            ENDDO
            B(J,J) = B(J,J) + 1.0D0 + ATEMP + BTEMP
  250    CONTINUE
C**********************
         CALL MATINV(B)
C**********************
         DO K=1,NWT
            ATEMP = TTU(K)/DTAU
            TTC(K,1) = 0.0D0
            DO J=1,NWT
               TTD(J,K,1) = B(J,K)*ATEMP
               TTC(K,1) = TTC(K,1) + B(K,J)*TTQ(J)
            ENDDO
         ENDDO
         GO TO 290
C*    Redo option
  270    QTEMP = 0.25D0*PITAU(1)*0.5D0*DTAU*DTAU
         DO J=1,NWT
            TTC(J,1) = 0.0D0
            DO K=1,NWT
               TTC(J,1) = TTC(J,1) + TTD(J,K,1)*QTEMP*TTF(K,1)/TTU2(K)
            ENDDO
         ENDDO
C*                          
C*    Continue through all depth points
  290    CONTINUE
         DO 360 II=2,NTT1
            DTEMP = 2.0D0/(TAU(II+1)-TAU(II-1))
            DTA = DTEMP/(TAU(II)-TAU(II-1))
            DTC = DTEMP/(TAU(II+1)-TAU(II))
            IF (LREDO) GO TO 330
            DO J=1,NWT
               ATEMP = TTU2(J)*DTA
               CTEMP = TTU2(J)*DTC
               TTQ(J) = 0.25D0*PITAU(II)*TTF(J,II) + ATEMP*TTC(J,II-1)
               DO K=1,NWT
                  B(J,K) = -ATEMP*TTD(J,K,II-1) - PITAU(II)*TTP(J,K)
               ENDDO
               B(J,J) = B(J,J) + 1.0D0 + ATEMP + CTEMP
            ENDDO
C**********************
            CALL MATINV(B)
C**********************
            DO J=1,NWT
               CTEMP = TTU2(J)*DTC
               TTC(J,II) = 0.0D0
               DO K=1,NWT
                  TTC(J,II) = TTC(J,II) + B(J,K)*TTQ(K)
                  TTD(K,J,II) = B(K,J)*CTEMP
               ENDDO
            ENDDO
            GO TO 360
C*    Redo option
  330       CONTINUE
            DO J=1,NWT
               TTQ(J) = (0.25D0*PITAU(II)*TTF(J,II)/TTU2(J) +
     1                                 DTA*TTC(J,II-1))/DTC
            ENDDO
            DO J=1,NWT
               TTC(J,II) = 0.0D0
               DO K=1,NWT
                  TTC(J,II) = TTC(J,II) + TTD(J,K,II)*TTQ(K) 
               ENDDO
            ENDDO
  360    CONTINUE
C*                      
C*    Lower boundary condition is second order
         II = NTTDO
         DTAU = TAU(II) - TAU(II-1)
         FACTOR = 4.0D0*RFLECT/(1.0D0+RFLECT)
         CTEMP = FACTOR*0.25D0*TTFBOT
         DO 380 J=1,NWT
            ATEMP = TTU(J)/DTAU
            BTEMP = 0.50D0/ATEMP
            TTQ(J) = ATEMP*TTC(J,II-1) +
     x               BTEMP*0.25D0*PITAU(II)*TTF(J,II) + CTEMP
            DO K=1,NWT
               B(J,K) = - ATEMP*TTD(J,K,II-1) - BTEMP*PITAU(II)*TTP(J,K)
     1                  - FACTOR*TTU(K)*TTW(K)
            ENDDO
            B(J,J) = B(J,J) + ATEMP + BTEMP + 1.0D0
  380    CONTINUE
C**********************
         CALL MATINV(B)
C**********************
         DO J=1,NWT
            TTC(J,II) = 0.0D0
            DO K=1,NWT
               TTC(J,II) = TTC(J,II) + B(J,K)*TTQ(K)
            ENDDO
         ENDDO
C*    Do back solution to compute the specific intensities (SI)
         DO III=1,NTT1
            II = NTTDO-III
            DO J=1,NWT
               DO K=1,NWT
                  TTC(J,II) = TTC(J,II) + TTD(J,K,II)*TTC(K,II+1)
               ENDDO
            ENDDO
         ENDDO
         LREDO = .TRUE.
C*    Add on the diffuse field to the direct beam
         DO II=1,NTTDO
            XMEAN = FLTAU(II)
            DO J=1,NWT
               XMEAN = XMEAN + 4.0D0*TTW(J)*TTC(J,II)
            ENDDO
            TTEX(II) = TTEX(II) + WTCLD*XMEAN
         ENDDO
  440 CONTINUE
      DO II=1,NTT
         FLTAU(II) = TTEX(II)
      ENDDO
  460 CONTINUE
C*         
C*    Return to SOL
      RETURN
      END
