#ifdef APM
! !MODULE: apm_radf_mod
!
! !DESCRIPTION: Module APM\_RADF\_MOD solves 1-column radiative forcing.
!\\
!\\
! !INTERFACE:
!
      MODULE APM_RADF_MOD
!
! !USES:
      IMPLICIT NONE
      PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC  :: APM_RADF
!
! !PRIVATE MEMBER FUNCTIONS:
!
! HISTORY:
!
!  MAR 17, 2011 - X. Ma      - preparation for call to radiation
!                              transfer model
!
       CONTAINS
!
      SUBROUTINE APM_RADF(JJ,II,LAY,LEV,NBS,NBL,
     1                 IAER,EXT,OMGA,G,ABSO,
     2                 GT,ALB,SCOS,TEMP,WATER,OZONE,PRES,AIRD,
     3                 CLD,CLDLIQ,CLDICE,RCLDL,YCOD,
     +                 RWCL,RWCI,
     4                 RCST,RCSB,RFST,RFSG,RFSA,RHRS,
     5                 RCLT,RCLB,RFUL,RFDL,RFLA,RHRL)
!
      INTEGER, PARAMETER :: ILG = 1
      INTEGER, PARAMETER :: IL1 = 1
      INTEGER, PARAMETER :: IL2 = 1
!
      INTEGER :: I,J,K,IB,LAY,LEV,NBS,NBL
      INTEGER :: II,JJ
      INTEGER :: IAER
!
      REAL  :: GT,SCOS
!      REAL  :: ALB
      REAL  :: ALB(NBS)
!
      REAL  :: TEMP(LAY),WATER(LAY),OZONE(LAY)
      REAL  :: PRES(LEV)
      REAL  :: CLD(LAY),CLDLIQ(LAY),CLDICE(LAY)
      REAL  :: AIRD(LAY)
      REAL  :: RCLDL(LAY)  !Yu+
      REAL  :: YCOD(LAY)  !Yu+
!
      REAL  :: SHTJ(ILG,LEV)
      REAL  :: T   (ILG,LAY)
      REAL  :: Q   (ILG,LAY)
      REAL  :: O3  (ILG,LAY)
!
      REAL  :: EXT (LAY,NBS)
      REAL  :: OMGA(LAY,NBS)
      REAL  :: G   (LAY,NBS)
      REAL  :: ABSO(LAY,NBL)
!
      REAL  :: GTEMP (ILG)
      REAL  :: SALB  (ILG,NBS)
!
      REAL  :: CCLD  (ILG,LAY),CLDLWC(ILG,LAY),CLDIWC(ILG,LAY)
      REAL  :: RHOP  (ILG,LAY)
!
      REAL  :: EXTA  (ILG,LAY,NBS)
      REAL  :: EXOMA (ILG,LAY,NBS)
      REAL  :: EXOMGA(ILG,LAY,NBS)
      REAL  :: FA    (ILG,LAY,NBS)
      REAL  :: ABSA  (ILG,LAY,NBL)
!
! output
!
      REAL  :: CST(ILG),CSB(ILG)
      REAL  :: FST(ILG),FSG(ILG),FSA(ILG)
      REAL  :: HRS(ILG,LAY),HRL(ILG,LAY)
      REAL  :: CLT(ILG),CLB(ILG)
      REAL  :: FUL(ILG),FDL(ILG),FLA(ILG)
      REAL  :: WCL(ILG),WCI(ILG)
!
      REAL  :: RCST,RCSB,RFST,RFSG,RFSA
      REAL  :: RHRS(LAY),RHRL(LAY)
      REAL  :: RCLT,RCLB,RFUL,RFDL,RFLA
      REAL  :: RWCL,RWCI
!
      DO IB = 1, NBS
      DO K  = 1, LAY
      DO I  = IL1,IL2
       EXTA  (I,K,IB) = EXT(K,IB)
       EXOMA (I,K,IB) = EXT(K,IB)*OMGA(K,IB)
       EXOMGA(I,K,IB) = EXT(K,IB)*OMGA(K,IB)*G(K,IB)
       FA    (I,K,IB) = EXT(K,IB)*OMGA(K,IB)*G(K,IB)*G(K,IB)
      ENDDO
      ENDDO
      ENDDO

      DO IB = 1, NBL
      DO K  = 1, LAY
      DO I  = IL1,IL2
       ABSA  (I,K,IB) = ABSO(K,IB)
      ENDDO
      ENDDO
      ENDDO
!
      DO IB  = 1, NBS
      DO I   = IL1, IL2
!       SALB(I,IB) = ALB
       SALB(I,IB) = ALB(IB)
      ENDDO
      ENDDO
!
      DO K  = 1, LAY
      DO I  = IL1, IL2
       SHTJ(I,K) = PRES (K)   !Yu+ in mb
       T   (I,K) = TEMP (K)
       Q   (I,K) = WATER(K)
       O3  (I,K) = OZONE(K)
      ENDDO
      ENDDO
!
      DO I  = IL1, IL2
       SHTJ(I,LEV) = PRES (LAY+1)
      ENDDO
!
!      write(*,*) 'Profile',IAER,II,JJ
      DO K  = 1, LAY
      DO I  = IL1, IL2
       CCLD (I,K)   = CLD (K)
       CLDLWC (I,K) = CLDLIQ (K)
       CLDIWC (I,K) = CLDICE (K)
       RHOP (I,K) = AIRD (K)
!      write(*,900) K,T(I,K),SHTJ(I,K),Q(I,K),O3(I,K),CCLD(I,K)
900   format(i5,2F10.2,2E10.2,F10.2)
      ENDDO
      ENDDO
!
      DO I = IL1,IL2
       GTEMP(I) = GT
      ENDDO
!
      CALL RAD1D(JJ,II,LAY,LEV,ILG,IL1,IL2,NBS,NBL,
     1           IAER,EXTA,EXOMA,EXOMGA,FA,ABSA,
     2           SHTJ,T,Q,O3,GTEMP,SALB,SCOS,
     3           CCLD,RHOP,CLDLWC,CLDIWC,
     +           WCL,WCI,RCLDL,YCOD,
     4           CST,CSB,FST,FSG,FSA,HRS,
     5           CLT,CLB,FUL,FDL,FLA,HRL)
!
      DO I = IL1, IL2
       RCST = CST(I)
       RCSB = CSB(I)
       RFST = FST(I)
       RFSG = FSG(I)
       RFSA = FSA(I)
       RCLT = CLT(I)
       RCLB = CLB(I)
       RFUL = FUL(I)
       RFDL = FDL(I)
       RFLA = FLA(I)
       RWCL = WCL(I)
       RWCI = WCI(I)
      ENDDO
!
      DO K = 1, LAY
      DO I = IL1, IL2
       RHRS(K) = HRS(I,K)
       RHRL(K) = HRL(I,K)
      ENDDO
      ENDDO 
!
      END SUBROUTINE APM_RADF

      END MODULE APM_RADF_MOD
C

        SUBROUTINE RAD1D(JJ,II,LAY,LEV,ILG,IL1,IL2,NBS,NBL,
     1                   IAER,EXTA,EXOMA,EXOMGA,FA,ABSA,
     2                   SHTJ,T,Q,O3,GT,SALB,SCOS,
     3                   CCLD,RHO,CLDLWC,CLDIWC,
     +                   WCL,WCI,RCLDL,YCOD,
     4                   CST,CSB,FST,FSG,FSA,HRS,
     5                   CLT,CLB,FUL,FDL,FLA,HRL)

C
!      USE APM_INIT_MOD, ONLY : YRLIQ,YRICE

      IMPLICIT REAL (A-H,O-Z),
     +INTEGER (I-N)
C Yu++
C
C     * INPUT ARRAYS FOR AEROSOL OPTICAL PROPERTIES
C
      REAL EXTA(ILG,LAY,NBS), EXOMA(ILG,LAY,NBS), EXOMGA(ILG,LAY,NBS),
     1     FA(ILG,LAY,NBS)
C
      REAL ABSA(ILG,LAY,NBL)
C
C     * INPUT/OUTPUT ARRAYS FOR PHYSICS.
C
C mxy+
      REAL FSA(ILG),DSHJ(ILG,LAY),HRSM(ILG,LAY)
      REAL FLA(ILG),HRLM(ILG,LAY)
      REAL WCL(ILG),WCI(ILG)
      REAL RCLDL(LAY)   !Yu++
      REAL YCOD(LAY)   !Yu++
C mxy+
      REAL FST(ILG),RHO(ILG,LAY)
      REAL FSG(ILG), FSD(ILG), FSF(ILG), FSV(ILG), FSI(ILG),
     1     ALBPLA(ILG), FDL(ILG), FUL(ILG), HRS(ILG,LAY), HRL(ILG,LAY),
     2     CST(ILG), CSB(ILG), CLT(ILG), CLB(ILG), PAR(ILG)
C
      REAL PRESSG(ILG), SHTJ(ILG,LEV), SHJ(ILG,LAY), TFULL(ILG,LEV), 
     1     T(ILG,LAY), GT(ILG), O3(ILG,LAY), O3TOP(ILG), Q(ILG,LAY),
     2     RMU(ILG)
C
      REAL SALB(ILG,NBS), EM0(ILG)
C
      REAL TAUCS(ILG,LAY,NBS), OMCS(ILG,LAY,NBS), GCS(ILG,LAY,NBS),
     1     TAUCL(ILG,LAY,NBL), OMCL(ILG,LAY,NBL), GCL(ILG,LAY,NBL),
     2     CLD(ILG,LAY),lwf(ilg,lev), FLXU(ilg,lev), FLXD(ilg,lev)
C
       REAL FSLO(ILG), FSAMOON(ILG)
C
      REAL P     (ILG,LAY), PFULL  (ILG,LEV)
      LOGICAL LCSW, LCLW
C
      REAL CCLD(ILG,LAY), CLDT(ILG)
      REAL CLDLWC(ILG,LAY), CLDIWC(ILG,LAY)
      REAL WCDW(ILG,LAY), WCDI(ILG,LAY), WCLW(ILG,LAY), WCLI(ILG,LAY),
     1     RADEQV(ILG,LAY), RADEQVW(ILG,LAY), RADEQVI(ILG,LAY)
C
      REAL N2O_PPM
      COMMON /RADCON/ SOLAR_C
      COMMON /TRACE / RMCO2, RMCH4, RMN2O, RMO2, RMF11, RMF12, RMF113,
     1                RMF114

      DATA LCSW /.true./, LCLW /.true./
C
C----------------------------------------------------------------------C
C     DIFFERENT ATMOSPHERE PROFILE                                     C
C----------------------------------------------------------------------C
C
      DO I = IL1, IL2
       TFULL(I,LEV) = GT(I)
      ENDDO
C
C----------------------------------------------------------------------C
C     CALCULATION OF THE MEAN O3 MIXING RATIO ABOVE MODEL TOP LEVEL    C
C----------------------------------------------------------------------C
C    
      DO K = 1, LAY 
      DO I = IL1, IL2
c       O3(I,K) =OOO (k)
c       Q (I,K) =WWW (k)
c      if(ii.eq.60.and.jj.eq.33) then
c       write(*,*) 'O3(I,K)',K,OOO(k),O3(I,K)
c      endif
      ENDDO 
      ENDDO
C
      DO I = IL1, IL2
        O3TOP(I) = 0.0
      ENDDO
C
C----------------------------------------------------------------------C
C     INITIATION OF THE INPUT                                          C
C----------------------------------------------------------------------C
C  
      DO I = IL1, IL2

        PRESSG(I) = 100.0
        RMU(I)    = MAX(SCOS,0.0) 
C
c... surface albedo
        DO J = 1, 4
C          SALB(I,J) = 0.1
        ENDDO
c... surface emissivity
        EM0(I) = 1.0
      ENDDO
C
      DO K = 1, LAY
        DO I = IL1, IL2
          SHJ(I,K)   = SQRT(SHTJ(I,K) * SHTJ(I,K+1))
          DSHJ(I,K)  = SHTJ(I,K+1) - SHTJ(I,K)
          IF (K .EQ. 1)                                             THEN
            TFULL(I,K) = 0.5 * (3.0 * T(I,1) - T(I,2))
          ELSE
            TFULL(I,K) = 0.5 * (T(I,K-1) + T(I,K))
          ENDIF
        ENDDO
      ENDDO
      DO K = 1, LAY
        DO I = IL1, IL2
          X      =  0.01 * PRESSG(I)
          P(I,K) =  SHJ (I,K) * X
        ENDDO
      ENDDO
C
      SOLAR_C  = 1365.6417 
      CO2_PPM  = 350. * 1.E-06 
      CH4_PPM  = 2.7 * 1.E-06
      N2O_PPM  = 1.2700 * 1.E-06
c      F11_PPM  = 0.E-3 * 1.E-06
c      F12_PPM  = 0.3618E-3 * 1.E-06
      F11_PPM  = 0.268E-3 * 1.E-06
      F12_PPM  =  0.535E-3 * 1.E-06
      F113_PPM = 0.084E-3 * 1.E-06
      F114_PPM = 0.015E-3 * 1.E-06
C
C----------------------------------------------------------------------C
C     INPUT TRACE GAS CONCENTRATIONS IN UNIT PPMV,                     C
C     PARTS PER MILLION BY VOLUME, TRANSFORM TO MASS MIXING RATIO.     C
C     THE SAME AS WATER VAPOR AND OZONE.                               C
C     1.5188126 = 44.    / 28.97                                       C
C     0.5522955 = 16.    / 28.97                                       C
C     1.5188126 = 44.    / 28.97                                       C
C     O2 INPUT AS A CONSTANT, UNIT MIXING RATIO BY MASS                C
C     4.7418019 = 137.37 / 28.97                                       C
C     4.1736279 = 120.91 / 28.97                                       C
C     5.2440456 = 151.92 / 28.97                                       C
C     5.8301691 = 168.90 / 28.97                                       C
C     28.97 MOLECULAR WEIGHT OF AIR, E-06 PER MILLION                  C
C----------------------------------------------------------------------C
C
      RMCO2  = 1.5188126 * CO2_PPM
      RMCH4  = 0.5522955 * CH4_PPM
      RMN2O  = 1.5188126 * N2O_PPM
      RMO2   = 0.2315
      RMF11  = 4.7418019 * F11_PPM
      RMF12  = 4.1736279 * F12_PPM
      RMF113 = 5.2440456 * F113_PPM
      RMF114 = 5.8301691 * F114_PPM
C
      DO I = IL1, IL2
        DO K = 1, LAY
          RADEQVW(I,K) = 0.0
          RADEQVI(I,K) = 0.0
          WCLW(I,K)    = 0.0
          WCLI(I,K)    = 0.0
        ENDDO
      ENDDO
C
C----------------------------------------------------------------------C
C     CLOUD INPUT INFORMATION                                          C
C     RADEQVW: EFFECTIVE RADIUS FOR WATER CLOUD                        C
C     WCLW: LIQUID WATER PATH LENGTH (IN GRAM / M^2)                   C
C     RADEQVI: EFFECTIVE RADIUS FOR ICE CLOUD                          C
C     WCLI: ICE WATER PATH LENGTH (IN GRAM / M^2)                      C
C     CCLD: CLOUD FRACTION                                             C
C     FRACW: FRACTION OF LIQUID WATER TO LIQUID + ICE WATER            C
C                                                                      C
C     ANU: FACTOR OF GAMMA DISTRIBUTION FOR INHOMOGENOUS CLOUD         C
C     IMPORTANT NOTE: ANU IS GIVEN IN CLDIFM                           C  
C     IN THIS MODEL THE CLOUD TOP LEVEL HAS TO BE LOWER THAN 1 MB      C
C----------------------------------------------------------------------C
C
       RGAS = 287.04
       GRAV = 9.8
       RGOCP = 2./7.
C
c.............................................
C
      DO I = IL1, IL2
       WCL(I) =0.0
       WCI(I) =0.0
      ENDDO
C
      DO K = 1, LAY
      DO I = IL1, IL2

        DISTANCE     = (RGAS/GRAV) * T(I,K) *
     1                     LOG(SHTJ(I,K+1)/SHTJ(I,K))
!Yu+
!        IF((MOD(II,10).EQ.1).and.(MOD(JJ,10).EQ.1).and.(MOD(K,5).EQ.1))
!     &    THEN
!          WRITE(112,108)II,JJ,K,SHTJ(I,K),CCLD(I,K),CLDLWC(I,K),
!     &     RCLDL(K)
! 108      FORMAT(I3,I3,I3,15(1PE9.2))
!        ENDIF

c        RADEQVW(I,K) = 10.0 
!Yu        RADEQVW(I,K) = 5.89
!!        RADEQVW(I,K) = YRLIQ
!         RADEQVW(I,K) = 3.00  !ok
!        RADEQVW(I,K) = 2.5  !ok
!        RADEQVW(I,K) = 2.2   !ok
!        RADEQVW(I,K) = 2.1  !ok
!        RADEQVW(I,K) = 2.0  !ok
!        RADEQVW(I,K) = MAX(1.,RCLDL(K))   !not ok
!        RADEQVW(I,K) = MAX(2.,RCLDL(K))   !ok
        RADEQVW(I,K) = MAX(1.5,RCLDL(K))   !ok
!        RADEQVW(I,K) = MAX(1.2,RCLDL(K))   ! not ok
!        RADEQVW(I,K) = MAX(1.3,RCLDL(K))   ! not ok
!        RADEQVW(I,K) = MAX(1.4,RCLDL(K))   ! not ok
!        RADEQVW(I,K) = MIN(20.,RADEQVW(I,K))
        RADEQVW(I,K) = MIN(30.,RADEQVW(I,K))

        WCLW(I,K)    = CLDLWC(I,K) * DISTANCE
C
        RADEQVI(I,K) = 41.5
!        RADEQVI(I,K) = YRICE
        WCLI(I,K)    = CLDIWC(I,K) * DISTANCE
C
        WCL(I) = WCL(I) + WCLW(I,K)
        WCI(I) = WCI(I) + WCLI(I,K)

C set cloud cover as zero
C
c        CCLD(I,K) =0.0
C
C      if(ii.eq.60.and.jj.eq.33) then
C       write(*,800) K,T(I,K),SHTJ(I,K),
C     &    DISTANCE,CLDIWC(I,K),WCLW(I,K),WCLI(I,K),CCLD(I,K)
800    format(I5,2F10.2,4F12.2,3X,F12.4)
C      endif

       ENDDO
       ENDDO
C
      CALL CLOUDS (TAUCS, OMCS, GCS, TAUCL, OMCL, GCL, CCLD, CLD, 
     1             CLDT, WCLW, WCLI, RADEQVW, RADEQVI, WCDW, WCDI,
     2             IL1, IL2, ILG, LAY, LEV)
C
!Yu+
!      if(ii.eq.10.and.jj.eq.6) then
      I = 1
      DO K = 1, LAY
       YCOD(LAY-K+1)=TAUCS(1,K,1)*CCLD(1,K)
!       if(MOD(ii,20).eq.1.and.MOD(jj,10).eq.1.and.mod(K,2).eq.1) then
!        DISTANCE     = (RGAS/GRAV) * T(I,K) *
!     1                     LOG(SHTJ(I,K+1)/SHTJ(I,K))
!        write(112,801)ii,jj, K,T(I,K),SHTJ(I,K),DISTANCE,CCLD(I,K),
!     &    CLDLWC(I,K),WCLW(I,K),
!     &    CLDIWC(I,K),WCLI(I,K),TAUCS(1,K,1),TAUCS(1,K,2),
!     &    TAUCLW(ii,jj,LAY-K+1),TAUCLI(ii,jj,LAY-K+1)
!       endif
      ENDDO
     
801   format(I3,I3,I3,20(1PE9.2))


      CALL RADDRIV (FSG, FSD, FSF, FSV, FSI, ALBPLA, FDL, FUL,
     1              HRS, HRL, CST, CSB, CLT, CLB, PAR, PRESSG,
     2              SHTJ, SHJ, TFULL, T, GT, O3, O3TOP, Q, RMU,
     3              SALB, EM0, TAUCS, OMCS, GCS, TAUCL, OMCL, GCL,
     4              CLD, EXTA, EXOMA, EXOMGA, FA, ABSA, 
     8              FSLO, FSAMOON, LCSW, LCLW,
     +              FST, RHO, II, JJ, IAER,
     6              IL1, IL2, ILG, LAY, LEV, FLXU, FLXD)
C
C Aerosol absorption
C ------------------------------------ mxy +
      DO L = 1, LAY
      DO I = IL1,IL2
       HRSM(I,L) = MAX(HRS(I,L),0.)
       HRLM(I,L) = MAX(HRL(I,L),0.)
      ENDDO
      ENDDO
C
      DO I = IL1, IL2
       FSA(I) = DSHJ(I,1)*HRSM(I,1) 
       FLA(I) = DSHJ(I,1)*HRLM(I,1) 
      ENDDO
C
      DO L = 2, LAY
      DO I = IL1,IL2
       FSA(I) = FSA(I) + DSHJ(I,L)*HRSM(I,L)
       FLA(I) = FLA(I) + DSHJ(I,L)*HRLM(I,L)
      ENDDO
      ENDDO
C
C   * CONVERT ABSORBED ATMOSPHERIC RADIATIVE FLUX TO W/M2 AND ACCUMULATE.
C
      CPRES=RGAS/RGOCP
      CPOGRAV=CPRES/GRAV
      DO I = IL1, IL2
       CONST  = CPOGRAV*PRESSG(I)
       FSAM   = MAX(FSAMOON(I),0.)
       FSA(I) = CONST*FSA(I) + FSAM
       FLA(I) = CONST*FLA(I)
C       print*,'FSA',II,JJ,FSA(I)
      ENDDO
C-------------------------------------  mxy +
c      print*,'99999',IAER,II,JJ
      DO I = IL1, IL2
C        WRITE(*,201)
C        WRITE(*,200) CST(I), CSB(I)
C
c        WRITE(*,*) 'FST',FST(I)
c        WRITE(*,200)
c        WRITE(*,202)
c        WRITE(*,200) FSG(I),FSD(I),FSF(I),FSV(I),FSI(I),ALBPLA(I)

      ENDDO
c      WRITE(*,200)

      DO I = IL1, IL2
c        WRITE(*,203)
c        WRITE(*,200) - CLT(I), - CLB(I)
C
c        WRITE(*,200)
c        WRITE(*,204)
c        WRITE(*,200) FUL(I), FDL(I)
c        WRITE(*,200)
      ENDDO


  200 FORMAT(3X, F9.3, 2X, F9.3, 2X, F9.3, 2X, F9.3, 2X, F9.3, 2X, F9.3)
c      WRITE(*,205)
c      WRITE(*,206)

      DO K = 1, LAY
        DO I = IL1, IL2
c         if(K.ge.30) then
c          WRITE(*,207)  P(I,K), HRS(I,K) * 86400., HRL(I,K) * 86400.
c         endif
        ENDDO
      ENDDO
C
    
  201 FORMAT(7X,'CST', 8X, 'CSB')
  202 FORMAT(7X,'FSG', 8X, 'FSD', 8X, 'FSF',8X, 'FSV',8X,
     1      'FSI',7X, 'ALBPLA')
  203 FORMAT(7X,'CLT', 8X, 'CLB')
  204 FORMAT(7X,'FUL', 8X, 'FDL')
  205 FORMAT(5X, 'PRESSUE', 2X, 'SOLAR HEATING RATE', 2X,
     1           'INFRARED COOLING RATE' )
  206 FORMAT(6X, '(MB)', 9X, '(K/DAY)', 14X, '(K/DAY)' )
  207 FORMAT(1X, F10.5, 5X, F9.3, 12X, F9.3)
  208 FORMAT(1X, F10.5, 5X, F9.3, 12X, F9.3,12X,F9.3)

      RETURN
      END


      SUBROUTINE RADDRIV (FSG, FSD, FSF, FSV, FSI, ALBPLA, FDL, FUL,
     1                    HRS, HRL, CST, CSB, CLT, CLB, PAR, PRESSG,
     2                    SHTJ, SHJ, TFULL, T, GT, O3, O3TOP, Q, RMU,
     3                    SALB, EM0, TAUCS, OMCS, GCS, TAUCL, OMCL, GCL,
     4                    CLD, EXTA, EXOMA, EXOMGA, FA, ABSA, 
     8                    FSLO, FSAMOON, LCSW, LCLW,
     +                    FST,RHO,II,JJ,IAER,
     6                    IL1, IL2, ILG, LAY, LEV, FLXU, FLXD)
C
C----------------------------------------------------------------------C
C     FSG:         DOWNWARD FLUX ABSORBED BY GROUND.                   C
C     FSD:         DIRECT DOWNWARD FLUX AT THE SURFACE.                C
C     FSF:         DIFFUSE DOWNWARD FLUX AT THE SURFACE.               C
C     FSV:         VISIBLE DOWNWARD FLUX AT THE SURFACE.               C
C     FSI:         NEAR INFRARED DOWNWARD FLUX AT THE SURFACE.         C
C     PAR:         PHOTOSYNTHETIC ACTIVE RADIATION.                    C
C     ALBPLA:      PLANETARY ALBEDO.                                   C
C     CST/CSB:     NET CLEAR SKY SOLAR FLUX AT TOP / SURFACE           C
C     FUL/FDL:     UPWARD LW FLUX AT THE TOP / SURFACE                 C
C     CLT/CLB:     NET CLEAR SKY LONGWAVE FLUX AT TOP/SURFACE          C
C     HRS/HRL:     SOLAR HEATING RATE / LONGWAVE COOLING RATE          C
C     PRESSG:      PRESSURE AT GROUND IN UNIT PA                       C
C     SHTJ:        PRESSURE AT MODEL LEVELS                            C
C     SHJ:         PRESSURE AT MODEL LAYER CENTER                      C
C     TFULL/T:     TEMPERATURE AT MODEL LEVEL / LAYER CENTER           C
C     GT:          GROUND TEMPERATURE                                  C
C     O3:          OZONE MASS MIXING RATIO IN (G/G)                    C
C     O3TOP:       ACCUMULATED OZONE MASS ABOVE THE MODEL TOP          C
C     Q:           WATER VAPOR SPECIFIC HUMIDITY (MASS MIXING RATIO IN C
C                  SOME VERSIONS)                                      C
C     RMU:         COSINE OF SOLAR ZENITH ANGLE                        C
C     SALB:        SURFACE ALBEDO                                      C
C     TAUCS/TAUCL: CLOUD OPTICAL DEPTH FOR SOLAR/LONGWAVE              C
C     OMCS/OMCL:   CLOUD SINGLE SCATTERING ALBEDO FOR SOLAR / LONGWAVE C
C     GCS/GCL:     CLOUD ASYMMETRY FACTOR FOR SOLAR/LONGWAVE           C
C     CLD:         CLOUD FRACTION                                      C
C     RADJ:        LATITUDE                                            C
C     GC:          INDEX FOR OCAN / LAND / SEA ICE                     C
C     PBLT:        BOUNDARY LAYER HEIGHT                               C
C     EXTA:        EXTINCTION COEFFICIENT                              C
C     EXOMA:       EXTINCTION COEFFICIENT TIMES SINGLE SCATTERING      C 
C                  ALBEDO                                              C
C     EXOMGA: EXOMA TIMES ASYMMETRY FACTOR                             C
C     FA:     SQUARE OF ASYMMETRY FACTOR                               C
C     ABSA:   ABSORPTION COEFFICIENT                                   C
C     FSLO:        SOLAR INCOMING FLUX AT INFRARED RANGE (0-2500CM-1)  C
C     FSAMOON:     THE ENERGY ABSORBED BETWEEN TOA AND MODEL TOP LEVEL C
C----------------------------------------------------------------------C
C
      IMPLICIT REAL (A-H,O-Z),
     +INTEGER (I-N)
C
      PARAMETER (NBS = 4, NBL = 9)
C
C     * INPUT/OUTPUT ARRAYS FOR PHYSICS.
C
      REAL FST(ILG),RHO(ILG,LAY)
      REAL FSG(ILG), FSD(ILG), FSF(ILG), FSV(ILG), FSI(ILG),
     1     ALBPLA(ILG), FDL(ILG), FUL(ILG), HRS(ILG,LAY), HRL(ILG,LAY),
     2     CST(ILG), CSB(ILG), CLT(ILG), CLB(ILG), PAR(ILG)
      REAL CSALB(ILG), BUOY(ILG,LAY), CLDIWC(ILG,LAY)
C
      REAL PRESSG(ILG), SHTJ(ILG,LEV), SHJ(ILG,LAY), 
     1     TFULL(ILG,LEV), T(ILG,LAY), GT(ILG), O3(ILG,LAY),
     2     O3TOP(ILG), Q(ILG,LAY), RMU(ILG)
C
      REAL SALB(ILG,NBS)
C
      REAL TAUCS(ILG,LAY,NBS), OMCS(ILG,LAY,NBS), GCS(ILG,LAY,NBS),
     1     TAUCL(ILG,LAY,NBL), OMCL(ILG,LAY,NBL), GCL(ILG,LAY,NBL),
     2     CLD  (ILG,LAY), GC(ILG)
C
      REAL  FSLO(ILG), FSAMOON(ILG)
C
C   * INPUT ARRAYS FOR AEROSOL OPTICAL PROPERTIES FOR SHORTWAVE
C
      REAL EXTA(ILG,LAY,NBS), EXOMA(ILG,LAY,NBS), EXOMGA(ILG,LAY,NBS),
     1     FA(ILG,LAY,NBS)
C
C   * INPUT ARRAYS FOR AEROSOL OPTICAL PROPERTIES FOR LONGWAVE
C
      REAL ABSA(ILG,LAY,NBL)
C
      LOGICAL LCSW, LCLW
C
C     * GENERAL WORK ARRAYS.
C
      REAL PG     (ILG,LAY), QG    (ILG,LAY), FLXU  (ILG,LEV),
     1     FLXD   (ILG,LEV), P     (ILG,LAY), DP    (ILG,LAY),
     2     TAUR   (ILG,LAY), TAUG  (ILG,LAY), TAUA  (ILG,LAY),
     3     PFULL  (ILG,LEV), F1    (ILG,LAY), F2    (ILG,LAY),
     4     ANU    (ILG,LAY), URBF  (ILG,LAY), TAUOMA(ILG,LAY),
     5     TAUOMGA(ILG,LAY), DIP   (ILG,LAY), DIPG  (ILG,LAY),
     6     DT     (ILG,LAY)
      REAL REFL   (ILG,2,LEV), TRAN  (ILG,2,LEV) 
C
      REAL C1(ILG), C2(ILG), BS(ILG)
C
C     * GATHERED AND OTHER WORK ARRAYS USED GENERALLY BY SOLAR.
C
      REAL A1(ILG,12), A1G(ILG,12), CUMDTR(ILG,4,LEV)
C
      REAL TAUCSG(ILG,LAY), TAUOMC(ILG,LAY), TAUOMGC(ILG,LAY),
     1     PFULLG(ILG,LEV), O3G(ILG,LAY), CLDG(ILG,LAY), CLDMG(ILG,LAY),
     2     TG(ILG,LAY)
C
      REAL O3TOPG(ILG), ALBSUR(ILG), EM0(ILG), RMUG(ILG)
C
      INTEGER INPTG(ILG,LAY), INPTMG(ILG,LAY), NBLK(ILG,LAY)
      INTEGER ISUN(ILG)
C
C     * WORK ARRAYS USED GENERALLY BY LONGWAVE.
C
      REAL TAUCI(ILG,LAY), OMCI(ILG,LAY), GCI(ILG,LAY), BF(ILG,LEV), 
     1     CLDM(ILG,LAY)
C
      INTEGER INPT(ILG,LAY), INPTM(ILG,LAY), INPR(ILG,LAY), 
     1        NCD(ILG,LAY), NCU(ILG,LAY)
      INTEGER NCT(ILG), NCTG(ILG), NCUM(LAY), NCDM(LAY)
C
C     * BAND INFORMATION.
C
      REAL SFINPTL(NBL)
      INTEGER KGS(NBS), KGSGH(NBS), KGL(NBL), KGLGH(NBL)
C
      LOGICAL GH
C
      COMMON /ECCENT/ RRSQ
      COMMON /RADCON/ SOLAR_C
      COMMON /TRACE / RMCO2, RMCH4, RMN2O, RMO2, RMF11, RMF12, RMF113,
     1                RMF114
C
C----------------------------------------------------------------------C
C     FOR HRCOEF, 9.80665 / 1004.64 / 100 = 9.761357E-05, IN (K / SEC),C
C     SINCE WE USE DP (DIFF IN PRESSURE) INSTEAD OF DIFF IN METER,     C
C     THERE IS A FACTOR 1.02. THUS 9.761357E-05 * 1.02 = 9.9565841E-05 C
C     UU3 = 3 * U * U, U = 1 / e^0.5                                   C
C----------------------------------------------------------------------C
C
      DATA HRCOEF, UU3, CUT / 9.9565841E-05, 1.1036383, 0.001 /
C
C     VARIABLES TO PASS TO INTRINSICS SUCH AS MIN/MAX:
C
      DATA QMIN     /1.E-7 /
      DATA OPEM10   /1.E-10/
C----------------------------------------------------------------------C
      PI                        =  3.1415926
C
C----------------------------------------------------------------------C
C     SCALE MEAN (ANNUAL) VALUE OF SOLAR CONSTANT BY RRSQ ACCOUNTING   C
C     FOR ECCENTRICITY (PASSED THROUGH COMMON BLOCK "ECCENT" - SEE     C
C     ROUTINE SDET2). THE SPECTRAL IRRADIANCE FOR MODEL IS 1367.9396   C
C     W / M^2  WHICH IS THE SOLAR ENERGY CONTAINED IN THE SPECTRAL     C
C     REGION 0.2 - 10 UM (50000 - 1000 CM) (BASED ON KURUCZ SOLAR DATA)C
C     FOR LONGWAVE, FROM BAND1 TO BAND9, THE SOLAR AND INFRARED        C
C     INTERACTION IS CONSIDERED. THE TOTAL SOLAR ENERGY CONSIDERED IN  C
C     THE INFRARED REGION IS 11.9096 W / M^2. SFINPTL IS THE INPUT     C
C     SOLAR FLUX IN EACH LONGWAVE BAND                                 C
C     THE SOLAR INPUT IN SHORTWAVE REGION IS 1367.9396 - 11.9096 =     C
C     1356.0300, THE SOLAR FRACTIONS FOR EACH BAND ARE SET IN GASOPTS  C
C----------------------------------------------------------------------C
C
C----------------------------------------------------------------------C
C     RRSQ IS THE SOLAR ORBIT PARAMETER, ADJUST IT IN ANY GLOBAL MODEL C
C----------------------------------------------------------------------C
C 
      DO IB = 1,NBS
      DO K =1,LAY
      DO I =IL1,IL2
c      if(ii.eq.60.and.jj.eq.33) then
c      WRITE(1003,*) EXTA(I,K,IB),EXOMA(I,K,IB),EXOMGA(I,K,IB),FA(I,K,IB)
c      endif
      ENDDO
      ENDDO
      ENDDO
C
      RRSQ = 1
      SOLARC                    =  SOLAR_C
      FRACS                     =  RRSQ * SOLARC / 1367.9396
      X                         =  FRACS / PI
C
      SFINPTL(1)                =  3.67839 * X
      SFINPTL(2)                =  2.79694 * X
      SFINPTL(3)                =  3.20284 * X
      SFINPTL(4)                =  1.13984 * X
      SFINPTL(5)                =  0.31893 * X
      SFINPTL(6)                =  0.35404 * X
      SFINPTL(7)                =  0.29578 * X
      SFINPTL(8)                =  0.99624E-01 * X
      SFINPTL(9)                =  0.23220E-01 * X
C
C----------------------------------------------------------------------C
C     INITIALIZATION                                                   C
C----------------------------------------------------------------------C
C
      DO 20 I = IL1, IL2
        FST(I)                  =  0.0
        FSG(I)                  =  0.0
        FSD(I)                  =  0.0
        FSF(I)                  =  0.0
        FSI(I)                  =  0.0
        FSV(I)                  =  0.0
        CST(I)                  =  0.0
        CSB(I)                  =  0.0
        PAR(I)                  =  0.0
        FSAMOON(I)              =  0.0
        FSLO(I)                 =  11.9096 * RMU(I) * FRACS 
        ALBPLA(I)               =  0.0
        PFULL(I,LEV)            =  0.01 * PRESSG(I) * SHTJ(I,LEV)
        C1(I)                   =  0.0
c        print*,'FSLO',I,FSLO(I),PRESSG(I),SHTJ(I,LEV),RMU(I)
   20 CONTINUE
C
      DO 30 K = 1, LAY
      DO 30 I = IL1, IL2
        HRS(I,K)                =  0.0
        X                       =  0.01 * PRESSG(I)
        P(I,K)                  =  SHJ (I,K) * X
        PFULL(I,K)              =  SHTJ(I,K) * X
C
C----------------------------------------------------------------------C
C     DETERMINE LEV1 FOR SOLAR RADIATION, LEV1 IS THE LEVEL DETERMININGC
C     HIGHEST LEVEL FOR FULL RADIATION CALCULATION, ABOVE LEV1 THE     C
C     RADIATION IS CALCULATED IN A SIMPLE WAY, LEV1 SHOULD BE BETWEEN  C
C     1 TO 0.5 MB                                                      C
C----------------------------------------------------------------------C
C
        PRESLEV1 =  0.5
        IF (PFULL(I,K) .GE. PRESLEV1)                               THEN
          C1(I)                 =  C1(I) + 1.0
          IF (C1(I) .EQ. 1.0)  INPT(I,1) =  K
        ENDIF
C
   30 CONTINUE
C
      LEV1                      =  LEV
C
      DO 31 I = IL1, IL2
        LEV1                    =  MIN (LEV1, INPT(I,1))
   31 CONTINUE
C
C----------------------------------------------------------------------C
C     THIS CODE CAN BE EXTENDED TO ABOUT 100 KM, IF THE MODEL TOP LEVELC
C     IS LOWER THAN THE MAXIMUM HEIGHT, THE CALCULATION CAN BE         C
C     SIMPLIFIED WITH LESS NUMBERS OF KGSGH AND KGLGH ACCOUNTED        C
C----------------------------------------------------------------------C
C
      DATA KGS   / 6, 4, 6, 4 /
      DATA KGL   / 1, 1, 2, 5, 2, 3, 3, 6, 4 /
C
      PTOP =  0.5
      DO 32 I = IL1, IL2
        PTOP = MIN (PFULL(I,1), PTOP)
   32 CONTINUE
C
      IF (PTOP .GT. 1.0)                                            THEN
        KGSGH(1) =  3
        KGSGH(2) =  3
        KGSGH(3) =  3
        KGSGH(4) =  6
C
        KGLGH(1) =  2
        KGLGH(2) =  1  
        KGLGH(3) =  2
        KGLGH(4) =  5
        KGLGH(5) =  3
        KGLGH(6) =  0 
        KGLGH(7) =  6
        KGLGH(8) =  2
        KGLGH(9) =  3
      ELSE
        KGSGH(1) =  3
        KGSGH(2) =  4
        KGSGH(3) =  4
        KGSGH(4) =  9
C
        KGLGH(1) =  5
        KGLGH(2) =  1
        KGLGH(3) =  3
        KGLGH(4) =  5
        KGLGH(5) =  4
        KGLGH(6) =  0 
        KGLGH(7) =  7
        KGLGH(8) =  3
        KGLGH(9) =  6
      ENDIF
C
C----------------------------------------------------------------------C
C     CALCULATE THE CLOUD PARAMETERS FOR SWTRAN AND LWTRAN             C
C----------------------------------------------------------------------C
C
      CALL CLDIFM (CLDM, TAUOMGC, ANU, A1, NCD, NCU, NBLK, NCT, NCUM,
     1             NCDM, CLD, CUT, MAXC, BUOY, CLDIWC, 
     2             IL1, IL2, ILG, LAY, LEV)
C
C----------------------------------------------------------------------C
C     DETERMINATION OF THE INTERPRETATION POINTS IN PRESSURE. INPT FOR C
C     28 REFERENCE LEVELS AND INPTM FOR 18 LEVELS                      C
C----------------------------------------------------------------------C
C
      CALL PREINTP (INPT, INPTM, DIP, A1(1,12), P, IL1, IL2, ILG, LAY)
C
      IF (LCSW)                                                     THEN
C
C----------------------------------------------------------------------C
C     DETERMINE WHETHER GRID POINTS ARE IN DAYLIGHT. GATHER THE        C
C     REQUIRED FIELD FOR DAYLIGHT REGION                               C
C----------------------------------------------------------------------C
C
      JYES = 0
      DO 200 I = IL1, IL2
        IF (RMU(I) .GT. 0.001)                                      THEN
          JYES                  =  JYES + 1
          ISUN(JYES)            =  I
        ENDIF
  200 CONTINUE
      LENGATH = JYES
C
C----------------------------------------------------------------------C
C     SKIP UNNECESSARY SOLAR                                           C
C----------------------------------------------------------------------C
C
      IF (LENGATH .EQ. 0) GO TO 499      
C
      DO 230 I = 1, LENGATH
        J = ISUN(I)
        O3TOPG(I)               =  O3TOP(J)
C
C----------------------------------------------------------------------C
C     C1 AND C2 ARE COEFFICIENTS FOR SWTRAN                            C
C                                                                      C
C     RMUG IS THE EFFECTIVE SOLAR ZENITH ANGLE, IT BECOMES IMPORTANT   C
C     WHEN SOLAR ZENITH ANGLE IS CLOSE TO 90 DEGREE (LI & SHIBATA JAS  C
C     2006), IT IS FOR GLOBAL MODEL, SET RMUG = RMU FOR PURE SINGLE    C
C     COLUMN MODEL AND RELATED INTERCOMPARISION                        C
C----------------------------------------------------------------------C
C
        RMUG(I)                 = (2.0 *  RMU(J) + 
     1                             SQRT(498.5225 * RMU(J) * RMU(J) + 
     2                             1.0)) / 24.35
C
        C1(I)                   =  0.75 * RMUG(I)
        C2(I)                   =  2.0 * C1(I) * RMUG(I)
C
        A1G(I,1)                =  A1(J,1)
        A1G(I,2)                =  A1(J,2)
        A1G(I,3)                =  A1(J,3)
        A1G(I,4)                =  A1(J,4)
        A1G(I,5)                =  A1(J,5)
        A1G(I,6)                =  A1(J,6)
        A1G(I,7)                =  1.0 - A1G(I,1) - A1G(I,2) - A1G(I,3)
        IF (A1G(I,2) .GE . CUT)                                     THEN
          A1G(I,8)              =  A1G(I,4) / A1G(I,2)
        ELSE
          A1G(I,8)              =  0.0
        ENDIF
C
        A1G(I,9)                =  0.0
        A1G(I,10)               =  0.0
        A1G(I,11)               =  0.0
        X                       =  A1G(I,3) + A1G(I,5) + A1G(I,6)
        IF (X .GE . CUT)                                            THEN
          IF (A1G(I,1) .GE . CUT)                                   THEN
            A1G(I,9)            =  A1G(I,6) / (X * A1G(I,1))
          ENDIF
          IF (A1G(I,2) .GE . CUT)                                   THEN
            A1G(I,10)           =  A1G(I,5) / (X * A1G(I,2))
          ENDIF
          A1G(I,11)             =  A1G(I,3) / X
        ENDIF
C
        A1G(I,12)               =  A1(J,12)
        NCTG(I)                 =  NCT(J)
        FLXU(I,LEV)             =  0.0
        FLXD(I,LEV)             =  0.0
        PFULLG(I,LEV)           =  PFULL(J,LEV)
C
C----------------------------------------------------------------------C
C     USING URBF SPACE FOR SCALING FACTOR USED IN SOLAR                C
C----------------------------------------------------------------------C
C
        URBF(I,1)               = (2.0 - RMUG(I)) ** 0.40
  230 CONTINUE
C
C----------------------------------------------------------------------C
C     USING A1(I,3) FOR RMU3                                           C
C----------------------------------------------------------------------C
C
      DO 240 I = 1, LENGATH
        X                       =  1.0 - RMUG(I)
        A1(I,3)                 =  X * X * X
C
        A1(I,4)                 =  0.0
  240 CONTINUE
C
      DO 255 K = 1, LAY
        KP1 = K + 1
        DO 250 I = 1, LENGATH
          J = ISUN(I)
          FLXU(I,K)             =  0.0
          FLXD(I,K)             =  0.0
          PFULLG(I,K)           =  PFULL(J,K)
C
C----------------------------------------------------------------------C
C     CONVERT FROM SPECIFIC HUMIDITY TO MIXING RATIO.  PLEASE NOTE IF  C
C     Q IS IN UNIT OF MIXING RATIO, TAKE OUT IT                        C
C     REUSING A1(I,5) FOR DT0                                          C
C----------------------------------------------------------------------C
C
          QMR                   =  Q(J,K) / (1.0 - Q(J,K))
          QG(I,K)               =  MAX (QMR, QMIN)
C
          CLDMG(I,K)            =  TAUOMGC(J,K)
          CLDG(I,K)             =  CLD(J,K)
C
          O3G(I,K)              =  O3(J,K)
          TG(I,K)               =  T(J,K)
          DT(I,K)               =  TG(I,K) - 250.0
          A1(I,5)               =  2.0 * TG(I,1) - TG(I,2) - 250.0
          PG(I,K)               =  P(J,K)
          DIPG(I,K)             =  DIP(J,K) 

          INPTG(I,K)            =  INPT(J,K)
          INPTMG(I,K)           =  INPTM(J,K)
c          print*,'TAUOMGC',J,K,Q(J,K),TAUOMGC(J,K),
c     &     CLD(J,K),O3(J,K),T(J,K),P(J,K)
C
C----------------------------------------------------------------------C
C     HERE DP = DIFP / G = RHO * DZ, WHERE DIFP IS THE LAYER PRESSURE  C
C     DIFFERENCE (IN MB), G IS THE GRAVITY CONSTANT, RHO IS AIR        C
C     DENSITY, AND DZ IS LAYER THICKNESS (IN CM). THEREFORE GAS MIXING C
C     RATIO * DP = GAS MASS * DZ. OR WE CAN CALL DP AS THE AIR MASS    C
C     PATH FOR A MODEL LAYER.                                          C
C     0.0102 = 1.02 * 0.01                                             C
C     1MB = 100 PASCAL = 1000 DYNES / CM^2,                            C
C     1.02 = (1000 DYNES / CM^2) / (980 CM / (SECOND^2)).              C
C     PRESSG, SURFACE PRESSURE IN UNIT PASCAL, SO WITH 0.01 FACTOR     C
C                                                                      C
C     REUSING BF AS A FACTOR FOR CLOUD SUBGRID VARIABILITY IN SOLAR    C
C----------------------------------------------------------------------C
C
          DP(I,K)               =  0.0102 * PRESSG(J) * 
     1                            (SHTJ(J,KP1) - SHTJ(J,K))
C          if(ii.eq.60.and.jj.eq.33) then
C           print*,'DP',K,DP(I,K)
C          endif
          IF (CLDG(I,K) .LT. CUT)                                   THEN
            BF(I,K)             =  0.0
          ELSE
            BF(I,K)             =  1.0 / (1.0 + 5.68 * ANU(J,K) ** 1.4)
          ENDIF
  250   CONTINUE
  255 CONTINUE
C
C----------------------------------------------------------------------C
C     SOLAR: 4 BAND FOR CLOUD, AEROSOL, AND RAYLEIGH,                  C
C     20 + 15 (20) MONOCHROMATIC CALCULATIONS FOR GAS AND RADIATIVE    C
C     TRANSFER                                                         C
C                                                                      C
C     FLXU:   ALL SKY SW UPWARD FLUX.                                  C
C     FLXD:   ALL SKY SW DOWNWARD FLUX.                                C
C     FSG:    DOWNWARD FLUX ABSORBED BY GROUND.                        C
C     FSD:    DIRECT DOWNWARD FLUX AT THE SURFACE.                     C
C     FSF:    DIFFUSE DOWNWARD FLUX AT THE SURFACE.                    C
C     FSV:    VISIBLE DOWNWARD FLUX AT THE SURFACE.                    C
C     FSI:    NEAR INFRARED DOWNWARD FLUX AT THE SURFACE.              C
C     PAR:    PHOTOSYNTHETIC ACTIVE RADIATION.                         C
C     ALBPLA: PLANETARY ALBEDO.                                        C
C     CST:    NET CLEAR SKY FLUX AT TOP.                               C
C     CSB:    NET CLEAR SKY FLUX AT SURFACE.                           C
C----------------------------------------------------------------------C
C
      DO 480 IB = 1, NBS  
C
C----------------------------------------------------------------------C
C     TWO SURFACE ALBEDO FOR ALL SKY AND CLEAR SKY PORTIONS IN A GRID  C
C     CELL, USER CAN USE TWO SEPARATE INPUT OF SALB                    C
C----------------------------------------------------------------------C
C 
        DO 300 I = 1, LENGATH
          J = ISUN(I)
          ALBSUR(I)             =  SALB(J,IB)
          CSALB(I)              =  SALB(J,IB)
  300   CONTINUE
C
C----------------------------------------------------------------------C
C     AEROSOL OPTICAL PROPERTIES. TAUA IS AEROSOL OPTICAL DEPTH        C
C----------------------------------------------------------------------C
C
        AOD = 0.
        DO 310 K = 1, LAY
        DO 310 I = 1, LENGATH
          J = ISUN(I)
          DZ                   =  DP(I,K) / RHO(I,K)   !mxy
          TAUA(I,K)             =  EXTA(J,K,IB) * DZ
          TAUOMA(I,K)           =  EXOMA(J,K,IB) * DZ
          TAUOMGA(I,K)          =  EXOMGA(J,K,IB) * DZ
          F1(I,K)               =  FA(J,K,IB)
          AOD = AOD+TAUA(1,K)
C----------------------------------------------------------------------C
C     SCALING THE CLOUD OPTICAL PROPERTIES DUE TO SUBGRID VARIABILITY  C
C     AND STANDARD SCALING FOR RADIATIVE TRANSFER                      C
C----------------------------------------------------------------------C
C
          IF (CLDG(I,K) .GE. CUT)                                   THEN
            IF (K .EQ. 1)                                           THEN
              TAUCI(I,K)        =  TAUCS(J,K,IB)
              X                 =  TAUCS(J,K,IB) +
     1                             9.2 * SQRT(TAUCS(J,K,IB)) 
            ELSE
              TAUCI(I,K)        =  TAUCI(I,K-1) + TAUCS(J,K,IB)
              X                 =  TAUCS(J,K,IB) +
     1                             9.2 * SQRT(TAUCI(I,K-1)) 
            ENDIF
C
            TAUCSG(I,K)         =  TAUCS(J,K,IB) / (1.0 + 0.185 *
     1                             X * URBF(I,1) * BF(I,K))
C
            C20                 =  TAUCSG(I,K) * OMCS(J,K,IB)
            TAUOMC(I,K)         =  TAUOMA(I,K) + C20
C
            C30                 =  C20 * GCS(J,K,IB)
            TAUOMGC(I,K)        =  TAUOMGA(I,K) + C30
            F2(I,K)             =  F1(I,K) + C30 * GCS(J,K,IB)
          ELSE
            TAUCI(I,K)          =  0.0
            TAUCSG(I,K)         =  0.0
            TAUOMC(I,K)         =  0.0
            TAUOMGC(I,K)        =  0.0
            F2(I,K)             =  0.0
          ENDIF
C 
  310   CONTINUE
c        if(ii.eq.60.and.jj.eq.33) then
c        print*,'AOD',AOD
c        endif
C
C----------------------------------------------------------------------C
C     RAYLEI, NEAR-IR RAYLEIGH SCATTERING, IT IS INDEPENDENT OF IG.    C
C----------------------------------------------------------------------C
C
        IF (IB .NE. 1)                                              THEN
          CALL RAYLEI (TAUR, IB, DP, 1, LENGATH, ILG, LAY)
        ENDIF
C
        GH = .FALSE.
C
        DO 400 IG = 1, KGS(IB) 
C
          IF (IB .EQ. 1)                                            THEN
C
C----------------------------------------------------------------------C
C     RAYLEV, VISIBLE RAYLEIGH SCATTERING, IT IS DEPENDANT ON IG.      C
C----------------------------------------------------------------------C
C
            CALL RAYLEV (TAUR, IG, DP, A1(1,3), 1, LENGATH, ILG, LAY)
C
C----------------------------------------------------------------------C
C     SOLAR ATTENUATION ABOVE THE MODEL TOP LAY. ONLY APPLY TO BAND    C
C     ONE FOR O3 AND O2. THIS IS TRUE ONLY FOR MODEL TOP LEVEL ABOVE   C
C     ABOUT 1 MB, WATER VAPOR CONTRIBUTION IS SMALL.                   C
C     REUSING A1(I,1) FOR MOON LAYER ATTENUATION                       C
C----------------------------------------------------------------------C
C
            CALL SATTENU (A1, IB, IG, RMUG, O3TOPG, QG, PFULLG, 
     1                    A1G(1,12), DT, A1(1,5), INPTG, GH, 1, LENGATH,
     2                    ILG)
          ELSE
            DO 320 I = 1, LENGATH
              A1(I,1)           =  1.0
  320       CONTINUE
          ENDIF
C
C----------------------------------------------------------------------C
C     DOWNWARD FLUX ABOVE 1 MB, FURTHER FLUX ATTENUATION FACTOR FOR    C
C     THE LOWER REGION                                                 C
C----------------------------------------------------------------------C
C
          IF (LEV1 .GT. 1)                                          THEN
            CALL STRANDN (TRAN, BS, A1, RMUG, DP, DT, O3G, A1(1,3), IB,
     1                    IG, LEV1, 1, LENGATH, ILG, LAY, LEV)
          ELSE
            DO 330 I = 1, LENGATH
              BS(I)             =  A1(I,1)
  330       CONTINUE
          ENDIF
C
          CALL GASOPTS (TAUG, GW, DP, IB, IG, O3G, QG, INPTMG, DIPG,
     1                  DT, A1(1,3), LEV1, GH, 1, LENGATH, ILG, LAY)
C mxy+
          DO K = 1,LAY
          DO I = IL1,IL2
           TAUG(I,K) = MAX(TAUG(I,K),1.E-20)
          ENDDO
          ENDDO
C mxy+
C
          CALL SWTRAN (REFL, TRAN, CUMDTR, BS, TAUA, TAUR, TAUG, TAUOMA,
     1                 TAUOMGA, F1, F2, TAUCSG, TAUOMC, TAUOMGC, CLDG,
     2                 CLDMG, A1G, RMUG, C1, C2, ALBSUR, CSALB,
     1                 NBLK, NCTG,II,JJ,IAER,
     3                 CUT, LEV1, 1, LENGATH, ILG, LAY, LEV)
C
          IF (LEV1 .GT. 1)                                          THEN
            CALL STRANUP (REFL, DP, DT, O3G, IB, IG, LEV1, 1, LENGATH,
     1                    ILG, LAY, LEV)
          ENDIF
C
C----------------------------------------------------------------------C
C     GATHER BACK THE REQUIRED FIELDS                                  C
C----------------------------------------------------------------------C
C
          RGW = GW * FRACS
          DO 350 I = 1, LENGATH
            J = ISUN(I)
            X                   =  A1G(I,7) * CUMDTR(I,1,LEV) +
     1                             A1G(I,1) * CUMDTR(I,2,LEV) +
     2                             A1G(I,2) * CUMDTR(I,3,LEV) +
     3                             A1G(I,3) * CUMDTR(I,4,LEV)
            A1(I,2)             =  RGW * RMUG(I)
            FSD(J)              =  FSD(J) + X * BS(I) * A1(I,2)
            CST(J)              =  CST(J) + (1.0 - REFL(I,1,1) *
     1                             A1(I,1)) * A1(I,2)
            CSB(J)              =  CSB(J) + (TRAN(I,1,LEV) - 
     1                             REFL(I,1,LEV)) * A1(I,2)
C
c          IF(II.eq.16.and.JJ.eq.12) THEN
c          print*,'before',II,JJ,IAER,FLXU(I,1),FLXD(I,1)
c          endif
            FLXU(I,1)           =  FLXU(I,1) + REFL(I,2,1) * A1(I,2)
            FLXD(I,1)           =  FLXD(I,1) + TRAN(I,2,1) * A1(I,2)
c          if(II.eq.16.and.JJ.eq.12) THEN
c          print*,'after',II,JJ,IAER,REFL(I,2,1),TRAN(I,2,1),A1(I,2)
c          endif
  350     CONTINUE
C
C----------------------------------------------------------------------C
C     HEATING RATE CALCULATION, FOR STABILITY IN CALCULATION, EACH IG  C
C     IS DONE SEPARATELY. HEATING RATE IN (K / SEC),                   C
C----------------------------------------------------------------------C
C
          DO 375 K = 1, LAY
            KP1 = K + 1
            DO 370 I = 1, LENGATH
              J = ISUN(I)
              DFNET             = (TRAN(I,2,K) - TRAN(I,2,KP1) -
     1                             REFL(I,2,K) + REFL(I,2,KP1)) * 
     2                             A1(I,2)
              HRS(J,K)          =  HRS(J,K) + HRCOEF * DFNET / DP(I,K)
C
              FLXU(I,KP1)       =  FLXU(I,KP1) + REFL(I,2,KP1) * A1(I,2)
              FLXD(I,KP1)       =  FLXD(I,KP1) + TRAN(I,2,KP1) * A1(I,2)
c          if(II.eq.16.and.JJ.eq.12) THEN
c          print*,'FLXUD',II,JJ,IAER,K,FLXU(I,KP1),FLXD(I,KP1)
c          endif
  370       CONTINUE
  375     CONTINUE
C
C----------------------------------------------------------------------C
C     FSAMOON IS THE ENERGY ABSORBED BETWEEN TOA AND MODEL TOP LEVEL.  C
C     A1(I,4) IS THE ADJUSTMENT FOR UPWARD FLUX FROM MODEL TOP LEVEL   C
C     TO TOA USED FOR PLANETARY ALBEDO                                 C
C----------------------------------------------------------------------C
C
          IF (IB .EQ. 1)                                            THEN
            DO 380 I = 1, LENGATH
              J = ISUN(I)
              X                 = (1.0 - A1(I,1)) * A1(I,2)
              FSAMOON(J)        =  FSAMOON(J) + X * (1.0 + REFL(I,2,1)) 
              A1(I,4)           =  A1(I,4) - X * REFL(I,2,1)
  380       CONTINUE
          ENDIF
C
          IF (IB .EQ. 1 .AND. IG .EQ. 2)                            THEN
            DO 390 I = 1, LENGATH
              PAR(ISUN(I))      =  FLXD(I,LEV)
  390       CONTINUE
          ENDIF
C
  400   CONTINUE
C
C----------------------------------------------------------------------C
C     IN ACCUMULATED SPACE WITH INTERVAL CLOSE TO 1, THE EXTINCTION    C
C     COEFFICIENTS IS EXTREMELY LARGE, THE CALCULATION PROCESS CAN BE  C
C     SIMPLIFIED BY IGNORING SCATTERING, REFLECTION, CLOUD AND AEROSOL.C
C----------------------------------------------------------------------C
C
        GH = .TRUE.
C
        DO 450 IG = 1, KGSGH(IB)
C
          CALL SATTENU (A1, IB, IG, RMUG, O3TOPG, QG, PFULLG, A1G(1,12),
     1                  DT, A1(1,5), INPTG, GH, 1, LENGATH, ILG)
C
          CALL STRANDNGH (TRAN, GWGH, A1, TAUA, TAUOMA, TAUCSG, TAUOMC,
     1                    CLDG, RMUG, DP, O3G, QG, IB, IG, INPTG, DIPG,
     2                    DT, LEV1, GH, CUT, 1, LENGATH, ILG, LAY, LEV)
C
          RGW = GWGH * FRACS
C
          DO 430 I = 1, LENGATH
            J = ISUN(I)
            A1(I,2)             =  RGW * RMUG(I)
            CST(J)              =  CST(J) + A1(I,2)
            CSB(J)              =  CSB(J) + TRAN(I,1,LEV) * A1(I,2)
C
            FSAMOON(J)          =  FSAMOON(J) + 
     1                             A1(I,2) * (1.0 - TRAN(I,2,1))
            FLXD(I,1)           =  FLXD(I,1) + A1(I,2) * TRAN(I,2,1)
  430     CONTINUE

          DO 445 K = 1, LAY
            KP1 = K + 1
            DO 440 I = 1, LENGATH
              J = ISUN(I)
              FLXD(I,KP1)       =  FLXD(I,KP1) + A1(I,2) * TRAN(I,2,KP1)
              HRS(J,K)          =  HRS(J,K) + HRCOEF * A1(I,2) * 
     1                            (TRAN(I,2,K) - TRAN(I,2,KP1)) / 
     2                             DP(I,K)
  440       CONTINUE
  445     CONTINUE
C
  450   CONTINUE
C
        IF (IB .EQ. 1)                                              THEN
          DO 460 I = 1, LENGATH
            FSV(ISUN(I))        =  FLXD(I,LEV)
  460     CONTINUE
        ENDIF
C
  480 CONTINUE
C
C----------------------------------------------------------------------C
C     GATHER BACK REQUIRED FIELD. FOR PLANETARY ALBEDO THE INCOMING    C
C     ENERGY OF 11.9096 * FRACS IS TOTALLY ABSORBED IN LONGWAVE PART   C
C----------------------------------------------------------------------C
C
      RSOLARC = RRSQ * SOLARC
      DO 490 I = 1, LENGATH
        J = ISUN(I)
        FST(J)                  =  FLXD(I,1) - FLXU(I,1)
        FSG(J)                  =  FLXD(I,LEV) - FLXU(I,LEV)
        FSI(J)                  =  FLXD(I,LEV) - FSV(J)
        FSF(J)                  =  FLXD(I,LEV) - FSD(J)
C
!20120713        CST(J)                  =  CST(J) + FSLO(J)
        CST(J)                  =  CST(J) 
        ALBPLA(J)               = (FLXU(I,1) + A1(I,4)) / 
     1                            (RSOLARC * RMUG(I))
  490 CONTINUE
C
  499 CONTINUE
C        
      ENDIF
C     (LCSW)
C
C----------------------------------------------------------------------C
C     LONGWAVE: 9 BAND FOR CLOUD, AEROSOL, CONTINUUM, AND PLANCK.      C
C     24+22 MONOCHROMATIC CALCULATIONS FOR GAS AND RADIATIVE TRANSFER  C
C                                                                      C
C     FLXU: ALL SKY LW UPWARD FLUX.                                    C
C     FLXD: ALL SKY LW DOWNWARD FLUX.                                  C
C     FUL:  UPWARD LW FLUX AT THE TOP.                                 C
C     FDL:  DOWN LW FLUX RECEIVED AT THE GROUND.                       C
C     CLT:  NET CLEAR SKY UPWARD FLUX AT THE TOP.                      C
C     CLB:  NET CLEAR SKY DOWNWARD FLUX AT THE SURFACE.                C
C----------------------------------------------------------------------C
C
      IF (LCLW)                                                     THEN
C
C----------------------------------------------------------------------C
C     CONVERT FROM SPECIFIC HUMIDITY TO MIXING RATIO (BOUNDED) AND     C
C     BOUND TEMPERATURE FOR PLANCK CALCULATION.                        C
C----------------------------------------------------------------------C
C
      DO 505 K = 1, LAY
        KP1 = K + 1
        DO 500 I = IL1, IL2
C
C----------------------------------------------------------------------C
C     CONVERT SPECIFIC HUMIDITY TO MIXING RATIO. TAKE OUT IT IF NOT IN C
C     UNIT OF SPECIFIC HUMIDITY                                        C
C----------------------------------------------------------------------C
C
          QMR                   =  Q(I,K) / (1.0 - Q(I,K))
          QG(I,K)               =  MAX (QMR, QMIN)
C
          HRL(I,K)              =  0.0
          DP(I,K)               =  0.0102 * PRESSG(I) * 
     1                            (SHTJ(I,KP1) - SHTJ(I,K))
          DT(I,K)               =  T(I,K) - 250.0 
  500   CONTINUE
  505 CONTINUE
C
      DO 510 I = IL1, IL2
        CLT(I)                  =  0.0
        CLB(I)                  =  0.0
C mxy+
        FUL(I)                  = 0.0
        FDL(I)                  = 0.0
C mxy+
        INPTG(I,1)              =  0
        ISUN(I)                 =  1 
        A1(I,5)                 =  2.0 * T(I,1) - T(I,2) - 250.0
  510 CONTINUE
C
C----------------------------------------------------------------------C
C     MCONT IS THE HIGHEST PRESSURE LEVEL FOR PERFORM THE CONTINUUM    C
C     CALCULATIONS (> 138.9440 MB). REUSING SPACES OF INPTG AND ISUN.  C
C----------------------------------------------------------------------C
C
      DO 520 K = 1, LEV
      DO 520 I = IL1, IL2
        FLXU(I,K)               =  0.0
        FLXD(I,K)               =  0.0
C
        IF (PFULL(I,K) .GE. 138.9440)                               THEN
          INPTG(I,1)            =  INPTG(I,1) + 1
          IF (INPTG(I,1) .EQ. 1) ISUN(I) =  K
        ENDIF
  520 CONTINUE
C
      MCONT = LEV
C
      DO 530 I = IL1, IL2
        MCONT                   =  MIN (ISUN(I), MCONT)
  530 CONTINUE
      MCONT = MCONT - 1
C
C----------------------------------------------------------------------C
C     DETERMINATION OF THE INTERPRETATION POINTS IN THE RATIO OF CO2   C
C     TO WATER VAPOR FOR TLINEHC. REUSE THE SPACE OF PG FOR DIR        C
C----------------------------------------------------------------------C
C
      CALL PREINTR (INPR, PG, QG, IL1, IL2, ILG, LAY)
C
      DO 900 IB =  1, NBL  !    
C
C----------------------------------------------------------------------C
C     USING C1 SPACE FOR SLWF WHICH IS THE INPUT SOLAR ENERGY IN THE   C
C     INFRARED REGION. TOTAL 11.9096 W / M^2 FROM STANDARD             C
C     CALCULATION                                                      C
C     SCALING CLOUD OPTICAL PROPERTIES FOR IR SCATTERING CALCULATION   C
C----------------------------------------------------------------------C
C
        DO 605 I = IL1, IL2
          IF (RMU(I) .GT. 0.0)                                      THEN
            C1(I)               =  RMU(I) * SFINPTL(IB)
          ELSE
            C1(I)               =  0.0
          ENDIF
 605    CONTINUE
C
        DO 610 K = 1, LAY
        DO 610 I = IL1, IL2
          DZ                   =  DP(I,K) / RHO(I,K)  !mxy
          TAUA(I,K)            =  ABSA(I,K,IB) * DZ   !mxy 
          TAUCI(I,K)            =  0.0 
          OMCI(I,K)             =  0.0 
          GCI(I,K)              =  0.0 
          F2(I,K)               =  0.0 
C
          IF (CLD(I,K) .GE. CUT)                                    THEN
            TAUCI(I,K)          =  TAUCL(I,K,IB)
            OMCI(I,K)           =  OMCL(I,K,IB) * TAUCI(I,K)
            F2(I,K)             =  GCL(I,K,IB) * GCL(I,K,IB)
            GCI(I,K)            = (GCL(I,K,IB) - F2(I,K)) / 
     1                            (1.0 - F2(I,K))
            GCI(I,K)            =  - 0.5 * (1.0 - UU3 * GCI(I,K))
          ENDIF
  610   CONTINUE
C      
C----------------------------------------------------------------------C
C     REUSING SPACE O3G FOR DBF                                        C
C----------------------------------------------------------------------C
C
        CALL PLANCK (BF, BS, URBF, A1(1,2), A1(1,3), O3G, TFULL, GT, IB,
     1               IL1, IL2, ILG, LAY, LEV)
C
        GH = .FALSE.
C
        DO 700 IG = 1, KGL(IB)
C
          CALL GASOPTL (TAUG, GW, DP, IB, IG, O3, QG, INPR, INPTM, 
     1                  MCONT, PG, DIP, DT, LEV1, GH, 
     2                  IL1, IL2, ILG, LAY)
C mxy+
          DO K = 1,LAY
          DO I = IL1,IL2
           TAUG(I,K) = MAX(TAUG(I,K),1.E-20)
          ENDDO
          ENDDO
C mxy+
C
          CALL LWTRAN (REFL, TRAN, C1, TAUCI, OMCI, GCI, F2, TAUA,
     1                 TAUG, BF, BS, URBF, O3G, EM0, CLD, CLDM, ANU, 
     2                 NCT, NCD, NCU, NCUM, NCDM, LEV1, MAXC,
     3                 CUT, IL1, IL2, ILG, LAY, LEV,
     +                 II,JJ)
C
          PGW = PI * GW
          DO 650 K = LEV1, LAY
            KP1 = K + 1
            DO 600 I = IL1, IL2
              FLXU(I,K)         =  FLXU(I,K) + REFL(I,2,K) * PGW
              FLXD(I,K)         =  FLXD(I,K) + TRAN(I,2,K) * PGW
C
              DFNET             =  TRAN(I,2,K) - TRAN(I,2,KP1) -
     1                             REFL(I,2,K) + REFL(I,2,KP1)
              HRL(I,K)          =  HRL(I,K) + 
     1                             HRCOEF * DFNET / DP(I,K) * PGW
             if(FLXU(I,K).gt.2000. .or. FLXU(I,K).lt.-1000.) then
               print*,'FLXU(I,K)',I,K,FLXU(I,K),REFL(I,2,K),PGW
             endif
             if(FLXD(I,K).gt.2000. .or. FLXD(I,K).lt.-1000.) then
               print*,'FLXD(I,K)',I,K,FLXD(I,K),TRAN(I,2,K),PGW
             endif
  600       CONTINUE
  650     CONTINUE
C
          DO 660 I = IL1, IL2
            FLXU(I,LEV)         =  FLXU(I,LEV) + REFL(I,2,LEV) * PGW
            FLXD(I,LEV)         =  FLXD(I,LEV) + TRAN(I,2,LEV) * PGW
            if(FLXU(I,LEV).gt.2000. .or. FLXU(I,LEV).lt.-1000.) then
             print*,'FLXU(I,LEV)',I,FLXU(I,LEV),REFL(I,2,LEV),PGW
            endif
            if(FLXD(I,LEV).gt.2000. .or. FLXD(I,LEV).lt.-1000.) then
             print*,'FLXD(I,LEV)',I,FLXD(I,LEV),TRAN(I,2,LEV),PGW
            endif
C
            CLT(I)              =  CLT(I) - REFL(I,1,LEV1) * PGW
            CLB(I)              =  CLB(I) - 
     1                            (REFL(I,1,LEV) - TRAN(I,1,LEV)) * PGW
           if(CLT(I).gt.2000. .or. CLT(I).lt.-1000.) then
            print*,'CLT(I)',I,CLT(I),REFL(I,1,LEV1),PGW
           endif
           if(CLB(I).gt.2000. .or. CLB(I).lt.-1000.) then
            print*,'CLB(I)',I,CLB(I),REFL(I,1,LEV),TRAN(I,1,LEV)
           endif
  660     CONTINUE
C
          IF (LEV1 .GT. 1)                                          THEN
            DO 680 K = LEV1 - 1, 1, - 1
              KP1 =  K + 1
              DO 670 I = IL1, IL2
                FLXU(I,K)       =  FLXU(I,K) + REFL(I,2,LEV1) * PGW
                FLXD(I,K)       =  FLXD(I,K) + C1(I) * PGW
  670         CONTINUE
  680       CONTINUE
          ENDIF
C
  700   CONTINUE
C
        IF (IB .NE. 6)                                              THEN
C
          GH = .TRUE.
C
          DO 800 IG = 1, KGLGH(IB)
C
            CALL GASOPTLGH (TAUG, GWGH, DP, IB, IG, O3, QG, INPT, 
     1                      MCONT, PG, DIP, DT, LEV1, GH, 
     2                      IL1, IL2, ILG, LAY)
C mxy+
          DO K = 1,LAY
          DO I = IL1,IL2
           TAUG(I,K) = MAX(TAUG(I,K),1.E-20)
          ENDDO
          ENDDO
C mxy+
C
C----------------------------------------------------------------------C
C     CONSIDER THE ATTENUATION FOR THE DOWNWARD FLUX ABOVE THE MODEL   C
C     TOP LEVEL. THIS IS IMPORTANT TO GET THE CORRECT COOLING RATE. IF C
C     THE MODEL TOP LEVEL PRESSURE IS LOWER THAN 0.01. THIS IS NOT     C
C     NECESSARY                                                        C
C----------------------------------------------------------------------C
C
            CALL LATTENU (A1, IB, IG, O3TOP, QG, PFULL, A1(1,12), DT, 
     1                    A1(1,5), INPT, IL1, IL2, ILG)
C
            DO 710 I = IL1, IL2
              TRAN0             =  EXP(- A1(I,1))
              if(TRAN0.gt.2000. .or. TRAN0.lt.-1000.) then
               print*,'TRAN0',I,TRAN0,A1(I,1)
              endif
C
              IF (PFULL(I,1) .GT. 0.001)                            THEN
                X               =  MAX(A1(I,1), OPEM10)
                ETA             =  1.0 / (1.0 + A1(I,3) / X)
                IF (ABS(ETA) .LT. 100.)                             THEN
                  C2(I)         =  C1(I) * TRAN0 + ETA *
     1                            (BF(I,1) - A1(I,2) * TRAN0)
                ELSE
                  C2(I)         =  C1(I) * TRAN0 + X * A1(I,2) * TRAN0
                ENDIF
              ELSE
                C2(I)           =  C1(I) * TRAN0
              ENDIF
  710       CONTINUE           
C
            CALL LWTRAGH (REFL, TRAN, C2, TAUCI, OMCI, TAUA, TAUG, BF, 
     1                    URBF, CLD, EM0, BS, CUT, 
     2                    IL1, IL2, ILG, LAY, LEV)
C
            PGW = PI * GWGH
            DO 740 K = 1, LAY
              KP1 = K + 1
              DO 730 I = IL1, IL2
                FLXU(I,K)       =  FLXU(I,K) + REFL(I,2,K) * PGW
                FLXD(I,K)       =  FLXD(I,K) + TRAN(I,2,K) * PGW
C
                DFNET           =  TRAN(I,2,K) - TRAN(I,2,KP1) -
     1                             REFL(I,2,K) + REFL(I,2,KP1)
                HRL(I,K)        =  HRL(I,K) +
     1                             HRCOEF * DFNET / DP(I,K) * PGW
  730         CONTINUE
  740       CONTINUE
C
C----------------------------------------------------------------------C
C     THE ATTENUATION FOR THE UPWARD FLUX ABOVE THE MODEL TOP IS NOT   C
C     CONSIDERED, SINCE THE IMPACT ON UPWARD FLUX IS VERY SMALL IF THE C
C     MODEL TOP IS ABOUT 1 MB OR HIGHER                                C
C----------------------------------------------------------------------C
C
            DO 750 I = IL1, IL2
              FLXU(I,LEV)       =  FLXU(I,LEV) + REFL(I,2,LEV) * PGW
              FLXD(I,LEV)       =  FLXD(I,LEV) + TRAN(I,2,LEV) * PGW
              CLT(I)            =  CLT(I) -  REFL(I,1,1) * PGW
              CLB(I)            =  CLB(I) - 
     1                            (REFL(I,1,LEV) - TRAN(I,1,LEV)) * PGW
             if(CLT(I).gt.2000. .or. CLT(I).lt.-1000. .or.
     &          CLB(I).gt.2000. .or. CLB(I).lt.-1000.) then
             print*,'CLT(I),CLB(I)',I,CLT(I),CLB(I),REFL(I,1,1),PGW,
     &        REFL(I,1,LEV),TRAN(I,1,LEV),FLXU(I,LEV),REFL(I,2,LEV),
     &        FLXD(I,LEV),TRAN(I,2,LEV)
             endif

  750       CONTINUE
C
  800     CONTINUE
C
        ENDIF
  900 CONTINUE  
C
      DO 950 I = IL1, IL2
        FDL(I)                  =  FLXD(I,LEV)
        FUL(I)                  =  FLXU(I,1)
  950 CONTINUE
C
      ENDIF
C     (LCLW)
C
      RETURN 
      END

      SUBROUTINE ATTENUE (ATTEN, COEF1, O3, Q, DP, DIP, DT, DT0, RMU,
     1                    INPT, N, NG, ISL, IL1, IL2, ILG)
C
C----------------------------------------------------------------------C
C     THIS SUBROUTINE CALCULATES THE DOWNWARD FLUX ATTENUATION ABOVE   C
C     THE MODEL TOP LEVEL                                              C
C     ISL = 1 FOR SOLAR, ISL = 2 FOR INFRARED.                         C
C     NG = 1, H2O; NG = 2, O3; NG = 3, CO2; NG = 6, O2                 C
C     ASSUMING THE TEMPERATURE AT 0.0005 MB IS 210 K                   C
C                                                                      C
C     ATTEN: FOR SOLAR: THE ATTENUATION FACTOR FOR DOWNWARD FLUX FROM  C
C            TOA TO THE MODEL TOP LEVEL; FOR LONGWAVE: THE OPTICAL     C
C            / DIFFUSE FACTOR                                          C
C     DP:    AIR MASS PATH FOR A MODEL LAYER (EXLAINED IN RADDRIV).    C
C     O3:    O3 MASS MIXING RATIO                                      C
C     Q:     WATER VAPOR MASS MIXING RATIO                             C
C     DIP:   INTERPRETATION FACTOR FOR PRESSURE BETWEEN TWO            C
C            NEIGHBORING STANDARD INPUT DATA PRESSURE LEVELS           C
C     DT:    LAYER TEMPERATURE - 250 K                                 C
C     DT0:   TEMPERATURE IN MOON LAYER - 250 K                         C
C     RMU:   COS OF SOLAR ZENITH ANGLE                                 C
C----------------------------------------------------------------------C
C
      IMPLICIT REAL (A-H,O-Z),
     +INTEGER (I-N)
C
      REAL ATTEN(ILG), COEF1(5,N), DP(ILG), O3(ILG), Q(ILG), DIP(ILG),
     1     DT(ILG), DT0(ILG), RMU(ILG), S1(ILG)
      INTEGER INPT(ILG)
C
      COMMON /TRACE / RMCO2, RMCH4, RMN2O, RMO2, RMF11, RMF12, RMF113,
     1                RMF114
C
      DATA RU / 1.6487213 /
      DATA N1 / 1 /
C=======================================================================
      IF (NG .EQ. 1)                                                THEN
        DO 100 I = IL1, IL2
          S1(I)    =  1.02 * Q(I)
 100    CONTINUE        
      ELSE IF (NG .EQ. 2)                                           THEN
        DO 200 I = IL1, IL2
          S1(I)    =  1.02 * O3(I)
 200    CONTINUE
      ELSE IF (NG .EQ. 3)                                           THEN
        DO 300 I = IL1, IL2
          S1(I)    =  1.02 * RMCO2
 300    CONTINUE
      ELSE IF (NG .EQ. 4)                                           THEN
        DO 400 I = IL1, IL2
          S1(I)    =  1.02 * RMCH4
 400    CONTINUE
      ELSE IF (NG .EQ. 6)                                           THEN
        DO 600 I = IL1, IL2
          S1(I)    =  1.02 * RMO2
 600    CONTINUE
      ENDIF 
C
      IF (ISL .EQ. N1)                                              THEN
       IF (INPT(1) .LT. 950)                                        THEN
        DO 1000 I = IL1, IL2
          N =  INPT(I)
          NM1 =  MAX (N - 1, N1)
          X1       =   COEF1(1,NM1) + DT0(I) * (COEF1(2,NM1) + DT0(I) *
     1                (COEF1(3,NM1) + DT0(I) * (COEF1(4,NM1) + DT0(I) *
     2                 COEF1(5,1))))
          X2       =   COEF1(1,N) + DT(I) * (COEF1(2,N) + DT(I) *
     1                (COEF1(3,N) + DT(I) * (COEF1(4,N) + DT(I) *
     2                 COEF1(5,N))))
C
          TAU      =  (X1 + (X2 - X1) * DIP(I)) * S1(I) * DP(I)
C
          ATTEN(I) =   EXP( - TAU / RMU(I))
 1000   CONTINUE
       ELSE
        N =  INPT(1) - 1000
        NM1 =  MAX (N - 1, N1)
        DO 1002 I = IL1, IL2
          X1       =   COEF1(1,NM1) + DT0(I) * (COEF1(2,NM1) + DT0(I) *
     1                (COEF1(3,NM1) + DT0(I) * (COEF1(4,NM1) + DT0(I) *
     2                 COEF1(5,1))))
          X2       =   COEF1(1,N) + DT(I) * (COEF1(2,N) + DT(I) *
     1                (COEF1(3,N) + DT(I) * (COEF1(4,N) + DT(I) *
     2                 COEF1(5,N))))
C
          TAU      =  (X1 + (X2 - X1) * DIP(I)) * S1(I) * DP(I)
C
          ATTEN(I) =   EXP( - TAU / RMU(I))
 1002   CONTINUE
       ENDIF
C
      ELSE
       IF (INPT(1) .LT. 950)                                        THEN
        DO 2000 I = IL1, IL2
          N =  INPT(I) - 1000
          NM1 =  MAX (N - 1, N1)
          X1       =   COEF1(1,NM1) + DT0(I) * (COEF1(2,NM1) + DT0(I) *
     1                (COEF1(3,NM1) + DT0(I) * (COEF1(4,NM1) + DT0(I) *
     2                 COEF1(5,1))))
          X2       =   COEF1(1,N) + DT(I) * (COEF1(2,N) + DT(I) *
     1                (COEF1(3,N) + DT(I) * (COEF1(4,N) + DT(I) *
     2                 COEF1(5,N))))
C
          TAU      =  (X1 + (X2 - X1) * DIP(I)) * S1(I) * DP(I)
C
          ATTEN(I) =   RU * TAU
 2000   CONTINUE
       ELSE
        N =  INPT(1) - 1000
        NM1 =  MAX (N - 1, N1)
        DO 2002 I = IL1, IL2
          X1       =   COEF1(1,NM1) + DT0(I) * (COEF1(2,NM1) + DT0(I) *
     1                (COEF1(3,NM1) + DT0(I) * (COEF1(4,NM1) + DT0(I) *
     2                 COEF1(5,1))))
          X2       =   COEF1(1,N) + DT(I) * (COEF1(2,N) + DT(I) *
     1                (COEF1(3,N) + DT(I) * (COEF1(4,N) + DT(I) *
     2                 COEF1(5,N))))
C
          TAU      =  (X1 + (X2 - X1) * DIP(I)) * S1(I) * DP(I)
C
          ATTEN(I) =   RU * TAU
 2002   CONTINUE
       ENDIF 
      ENDIF
C
      RETURN
      END

      SUBROUTINE CLDIFM (CLDMIN, CLDMAX, ANU, A1, NCD, NCU, NBLK, NCT,
     1                   NCUM, NCDM, CLD, CUT, MAXC,
     2                   BUOY, CLDIWC, IL1, IL2, ILG, LAY, LEV) 
C
      IMPLICIT REAL (A-H,O-Z),
     +INTEGER (I-N)
C
      REAL CLDMIN(ILG,LAY), CLDMAX(ILG,LAY), ANU(ILG,LAY), A1(ILG,10),
     1     CLD(ILG,LAY), BUOY(ILG,LAY), CLDIWC(ILG,LAY)
C
      INTEGER NCD(ILG,LAY), NCU(ILG,LAY), NBLK(ILG,LAY), NCT(ILG), 
     1        NCUM(LAY), NCDM(LAY), LEVC(ILG,LAY), INTG1(ILG), 
     2        INTG2(ILG)
      PARAMETER ( ALF=-3.E-05 )
C
      DATA ANUMIN, ANUMAX /0.5, 4./
C
C----------------------------------------------------------------------C
C     THIS SUBROUTINE DETERMINES THE INFO FOR CLOUD AND LEVEL INFO     C
C     FOR GASEOUS CALCULATION                                          C
C----------------------------------------------------------------------C
C
      DO 10 I = IL1, IL2
        INTG1(I)                =  0
        INTG2(I)                =  0
        NCT(I)                  =  LEV
   10 CONTINUE
C
C----------------------------------------------------------------------C
C     DETERMINE THE HIGHEST CLOUD LOCATION. NCT IS THE UPPER LEVEL OF  C
C     THE HIGHEST CLOUD,                                               C
C     DETERMINE THE NU (ANU) FACTOR FOR CLOUD SUB-GRID VARIABILITY     C
C     BASED ON CLOUD FRACTION.                                         C
C----------------------------------------------------------------------C
C
      MAXC = LEV
C
      DO 25 K = 1, LAY
        KM1 = K - 1
        DO 20 I = IL1, IL2
c         IF (CLDIWC(I,K) .GT. 0.)                                 THEN
            IF (CLD(I,K) .LE. 0.9)                                 THEN
              ANU(I,K)          =  1.0
            ELSEIF (CLD(I,K) .GT. 0.9 .AND. CLD(I,K) .LT. 1.0)     THEN
              ANU(I,K)          =  2.0
            ELSE
              ANU(I,K)          =  4.0
            ENDIF
c          ELSE
c            IF (BUOY(I,K) .LT. 0.)                                 THEN
c              ANU(I,K)          =  ALF / BUOY(I,K)
c              ANU(I,K)          =  MIN (MAX (ANUMIN, ANU(I,K)), ANUMAX)
c            ELSE
c              ANU(I,K)          =  4.
c            END IF
c          ENDIF
C
C----------------------------------------------------------------------C
C     MINIMUM ANU, IT IS EXTREMELY IMPORTANT TO ENSURE CONSISTENCY     C
C     BETWEEN THE DEFINITIONS OF ANU HERE AND THEIR SUBSEQUENT USE IN  C
C     LWTRAN,                                                          C
C     CLDMAX THE MAXIMUM CLOUD FRACTION FOR EACH CLOUD BLOCK.          C
C     CLDMIN THE MINIMUM CLOUD FRACTION FOR EACH CLOUD BLOCK.          C
C----------------------------------------------------------------------C
C
          IF (CLD(I,K) .LT. CUT)                                    THEN
            ANU(I,K)            =  1000.0
            CLDMAX(I,K)         =  0.0
          ELSE
            IF (K .EQ. 1)                                           THEN
              CLDMAX(I,K)       =  CLD(I,K)
            ELSE
              ANU(I,K)          =  MIN (ANU(I,KM1), ANU(I,K))
              CLDMAX(I,K)       =  MAX (CLDMAX(I,KM1), CLD(I,K))
            ENDIF
C
            INTG2(I)            =  INTG2(I) + 1
            IF (INTG2(I) .EQ. 1) NCT(I) = K
          ENDIF
C
   20   CONTINUE
   25 CONTINUE
C
      DO 40 I = IL1, IL2
        MAXC                    =  MIN (NCT(I), MAXC)
        LEVC(I,1)               =  0
        LEVC(I,2)               =  0
        A1(I,1)                 =  0.
        A1(I,2)                 =  0.
        A1(I,3)                 =  0.
   40 CONTINUE
C
C----------------------------------------------------------------------C
C     DETERMINE THE LAYER ORDER FOR EACH CLOUD BLOCK THROUGH DOWN AND  C
C     UP PATHS, NCD AND NCU.                                           C
C     DETERMINE THE TOTAL CLOUD FRACTIONS LOOKING FROM TOP AND SURFACE C
C     FOR ONE CLOUD BLOCK (A CLOUD OCCUPY SEVERAL LAYERS, CHOOSE THE   C
C     MINIMUM VALUE OF NU FOR THE BLOCK.                               C
C     NCT IS THE TOP LEVEL NUMBER FOR THE HIGHEST CLOUD                C
C     DETERMINE THE MINIMUM ANU                                        C
C----------------------------------------------------------------------C
C
      DO 65 K = 2, LEV
        KM1 = K - 1
        L = LEV - K + 1
        LP1 = L + 1
        DO 60 I = IL1, IL2
          IF (CLD(I,KM1) .LT. CUT)                                  THEN
            LEVC(I,1)           =  0
            NCD(I,KM1)          =  0
          ELSE
            LEVC(I,1)           =  LEVC(I,1) + 1
            NCD(I,KM1)          =  LEVC(I,1)
          ENDIF
C
          IF (CLD(I,L) .GE. CUT .AND. L .LT. LAY)                   THEN
            ANU(I,L)            =  MIN (ANU(I,LP1), ANU (I,L))
            CLDMAX(I,L)         =  MAX (CLDMAX(I,LP1), CLDMAX(I,L))
          ENDIF
   60   CONTINUE
   65 CONTINUE
C
      DO 75 L = LAY, 1, -1
        LP1 = L + 1
        DO 70 I = IL1, IL2
          IF (CLD(I,L) .LT. CUT)                                    THEN
            LEVC(I,2)           =  0
            NCU(I,L)            =  0
            NBLK(I,L)           =  0
            CLDMIN(I,L)         =  1.
          ELSE
            LEVC(I,2)           =  LEVC(I,2) + 1
            NCU(I,L)            =  LEVC(I,2)
            IF (NCU(I,L) .EQ. 1)                                    THEN
              INTG1(I)          =  INTG1(I) + 1
              NBLK(I,L)         =  INTG1(I)
              IF (NBLK(I,L) .GT. 3)  NBLK(I,L) =  3
              IF (NBLK(I,L) .EQ. 1)  A1(I,1)   =  CLDMAX(I,L)
              IF (NBLK(I,L) .EQ. 2)  A1(I,2)   =  CLDMAX(I,L)
              IF (NBLK(I,L) .EQ. 3)  A1(I,3)   =
     1                               MAX (A1(I,3), CLDMAX(I,L))
            ELSE
              NBLK(I,L)         =  NBLK(I,LP1)
            ENDIF
C
            IF (NCU(I,L) .EQ. 1)                                    THEN
              CLDMIN(I,L)       =  CLD(I,L)
            ELSE
              CLDMIN(I,L)       =  MIN (CLDMIN(I,LP1), CLD(I,L))
            ENDIF
          ENDIF
   70   CONTINUE
   75 CONTINUE
C
      DO 80 I = IL1, IL2
        X                       =  A1(I,3) * (1.0 - A1(I,1)) *
     1                            (1.0 - A1(I,2))
        A1(I,4)                 =  A1(I,1) * A1(I,2)
        A1(I,1)                 =  A1(I,1) * (1.0 - A1(I,2))
        IF (A1(I,3) .GE. X + A1(I,2))                               THEN
          Y                     =  A1(I,2)
          Z                     =  A1(I,3) - X - A1(I,2)
        ELSE
          Y                     =  A1(I,3) - X
          Z                     =  0.
        ENDIF
C
        IF (A1(I,3) .GE. X + A1(I,1))                               THEN
          A1(I,6)               =  A1(I,1)
          A1(I,5)               =  A1(I,3) - X - A1(I,6)
        ELSE
          A1(I,6)               =  A1(I,3) - X
          A1(I,5)               =  0.
        ENDIF
        A1(I,3)                 =  X
        A1(I,5)                 =  0.5 * (A1(I,5) + Y)
        A1(I,6)                 =  0.5 * (A1(I,6) + Z)
   80 CONTINUE
C
C----------------------------------------------------------------------C
C     DETERMINE THE MAXIMUM PORTION IN A CLOUD BLOCK                   C
C     DETERMINE THE MAXIMUM NUMBER FOR NCD AND NCU, FOR ITERATION IN   C
C     LONGWAVE                                                         C
C----------------------------------------------------------------------C
C
      DO 105 K = 1, LAY
        KM1 = K - 1
        NCUM(K)                 =  0
        NCDM(K)                 =  0
        DO 100 I = IL1, IL2
          IF (NCD(I,K) .GT. 1)                                      THEN
            CLDMIN(I,K)         =  MIN (CLDMIN(I,KM1), CLDMIN(I,K))
          ENDIF
C
          NCUM(K)               =  MAX (NCU(I,K), NCUM(K))
          NCDM(K)               =  MAX (NCD(I,K), NCDM(K))
  100   CONTINUE
  105 CONTINUE
C
      RETURN
      END

      SUBROUTINE GASOPTL (TAUG, GW, DP, IB, IG, O3, Q, INPTR, INPT, 
     1                    MCONT, DIR, DIP, DT, LEV1, GH, 
     2                    IL1, IL2, ILG, LAY)
C
C----------------------------------------------------------------------C
C     CALCULATION OF THE OPTICAL DEPTHS DUE TO NONGRAY GASEOUS         C
C     ABSORPTION FOR THE INFRARED, IN EACH LAYER FOR A GIVEN BAND IB   C
C     AND CUMULATIVE PROBABILITY GW.                                   C
C     FROM BAND1 TO BAND4, THE SOLAR AND INFRARED INTERACTION IS       C
C     CONSIDERED. THE TOTAL SOLAR ENERGY CONSIDERED IN THE INFRARED    C
C     REGION IS 11.9096 W / M^2                                        C
C     FOR GASES WITH CONSTANT MIXING RATIO:                            C
C     1 = H2O                                                          C
C     2 = O3                                                           C
C     3 = CO2                                                          C
C     4 = CH4                                                          C
C     5 = N2O                                                          C
C     6 = O2                                                           C
C     7 = CFC11                                                        C
C     8 = CFC12                                                        C
C     TLINE, ETC., DEAL WITH LINE ABSORPTION AND TCONTL AND TCONTHL    C
C     DEAL WITH WATER VAPOR CONTINUUM                                  C
C                                                                      C
C     TAUG:  GASEOUS OPTICAL DEPTH                                     C
C     DP:    AIR MASS PATH FOR A MODEL LAYER (EXLAINED IN RADDRIV).    C
C     O3:    O3 MASS MIXING RATIO                                      C
C     Q:     WATER VAPOR MASS MIXING RATIO                             C
C     DIR:   INTERPRETATION FACTOR FOR MASS RATIO OF H2O / CO2         C
C            BETWEEN TWO NEIGHBORING STANDARD INPUT RATIOS             C
C     DIP:   INTERPRETATION FACTOR FOR PRESSURE BETWEEN TWO            C
C            NEIGHBORING STANDARD INPUT DATA PRESSURE LEVELS           C
C     DT:    LAYER TEMPERATURE - 250 K                                 C
C     INPR:  NUMBER OF THE RATIO LEVEL FOR THE STANDARD 5 RATIOS       C
C     INPT:  NUMBER OF THE LEVEL FOR THE STANDARD INPUT DATA PRESSURES C
C     MCONT: THE HIGHEST LEVEL FOR WATER VAPOR CONTINUUM CALCULATION   C
C----------------------------------------------------------------------C
      IMPLICIT REAL (A-H,O-Z),
     +INTEGER (I-N)
C
      REAL TAUG(ILG,LAY)
C
      REAL DP(ILG,LAY), O3(ILG,LAY), Q(ILG,LAY), DIR(ILG,LAY), 
     1     DIP(ILG,LAY), DT(ILG,LAY)
      INTEGER INPTR(ILG,LAY), INPT(ILG,LAY)
      LOGICAL GH
C
      COMMON /TRACE / RMCO2, RMCH4, RMN2O, RMO2, RMF11, RMF12, RMF113,
     1                RMF114
C
      COMMON /BANDL1/  GW1(1),  CL1H2O(5,18), CL1CO2(5,18)
      COMMON /BANDL2/  GW2(1),  CL2H2O(5,18), 
     1                          CL2CS(5,3), CL2CF(5,3)
      COMMON /BANDL3/  GW3(2),  CL3H2O(5,18,2), CL3CH4(2),
     1                          CL3CS(5,4,2), CL3CF(5,4,2)
      COMMON /BANDL4/  GW4(5),  CL4H2O(5,18,5), CL4N2O(5,18,5),
     1                          CL4CH4(5,18,5),
     2                          CL4CS(5,4,5), CL4CF(5,4,5),
     3                          CL4F12(5)
      COMMON /BANDL5/  GW5(2),  CL5H2O(5,18,2), CL5O3(5,18,2),
     1                          CL5CS(5,4,2), CL5CF(5,4,2),
     2                          CL5F11(2), CL5F12(2)
      COMMON /BANDL6/  GW6(3),  CL6H2O(5,18,3), CL6F11(5,18,3),
     1                          CL6F12(5,18,3), CL6CS(5,4), CL6CF(5,4)
      COMMON /BANDL7/  GW7(3),  CL7H2OU(5,11,3), CL7H2OD(5,5,7,3),
     1                          CL7CO2U(5,11,3), CL7CO2D(5,5,7,3),
     2                          CL7N2O(3), CL7O3(3),
     3                          CL7CS(5,5,4,3), CL7CF(5,5,4,3)
      COMMON /BANDL8/  GW8(6),  CL8H2O(5,18,6),
     1                          CL8CS(5,6,4), CL8CF(5,6,4)
      COMMON /BANDL9/  GW9(4),  CL9H2O(5,18,4),
     1                          CL9CS(5,6,4), CL9CF(5,6,4)
C
C     * NUMBER OF VERTICAL LEVELS IN ABSORBER PRESSURE-BASED COEFFICIENT
C     * ARRAY ("M" REFERENCES NON-SATURATED BANDS ACTIVE BELOW 1 MB ONLY).
C
      DATA MTL /18/
C=======================================================================
      IF (IB .EQ. 1)                                                THEN
C
C----------------------------------------------------------------------C
C     BAND (2500 - 2200 CM^-1), NONGRAY GASEOUS ABSORPTION OF H2O AND  C
C     CO2.                                                             C
C----------------------------------------------------------------------C
C
      NG2 =  3
      CALL TLINE2 (TAUG, CL1H2O, CL1CO2, Q, O3, NG2, DP, DIP, DT, INPT,
     1             LEV1, GH, MTL, IL1, IL2, ILG, LAY)
C
C----------------------------------------------------------------------C
C     SIMPLY ADD THE N2O EFFECT                                        C
C----------------------------------------------------------------------C
C
      DO 100 K = LEV1, LAY
      DO 100 I = IL1, IL2
        CRMN2O    =  2.0437 * 1.E-07 * SQRT(RMN2O * 1.E+07) 
        FACT      =  Q(I,K) / (Q(I,K) + 8.E+04 * CRMN2O)
        TAUG(I,K) =  TAUG(I,K) + (754.9786 + 10141.5049 * FACT * FACT) *
     1               CRMN2O * DP(I,K)
  100 CONTINUE
C
      GW =  GW1(IG)
C
      ELSE IF (IB .EQ. 2)                                           THEN
C
C----------------------------------------------------------------------C
C     BAND (2200 - 1900 CM^-1), NONGRAY GASEOUS ABSORPTION OF H2O + N2OC
C----------------------------------------------------------------------C
C
      NG =  1
      INIT = 2
      CALL TLINE1 (TAUG, CL2H2O, Q, NG, DP, DIP, DT, INPT, LEV1, GH, 
     1             MTL, INIT, IL1, IL2, ILG, LAY)
C
      LC =  3
      CALL TCONTL1 (TAUG, CL2CS, CL2CF, Q, DP, DIP, DT, LC, INPT, MCONT,
     1              GH, IL1, IL2, ILG, LAY)
C
C----------------------------------------------------------------------C
C     SIMPLY ADD THE N2O EFFECT                                        C
C----------------------------------------------------------------------C
C
      DO 200 K = LEV1, LAY
      DO 200 I = IL1, IL2
        CRMN2O    =  2.0437 * 1.E-07 * SQRT(RMN2O * 1.E+07)
        FACT      =  Q(I,K) / (Q(I,K) + 72000. * CRMN2O)
        TAUG(I,K) =  TAUG(I,K) + (93. + 3500. * FACT * FACT) *
     1               CRMN2O * DP(I,K)
  200 CONTINUE
C
      GW =  GW2(IG)
C
      ELSE IF (IB .EQ. 3)                                           THEN
C
C----------------------------------------------------------------------C
C     BAND (1900 - 1400 CM^-1), NONGRAY GASEOUS ABSORPTION OF H2O.     C
C----------------------------------------------------------------------C
C
      NG =  1
      INIT = 2
      CALL TLINE1 (TAUG, CL3H2O(1,1,IG), Q, NG, DP, DIP, DT, INPT,
     1             LEV1, GH, MTL, INIT, IL1, IL2, ILG, LAY)
C
      LC =  4
      CALL TCONTL1 (TAUG, CL3CS(1,1,IG), CL3CF(1,1,IG), Q, DP, DIP, DT, 
     1              LC, INPT, MCONT, GH, IL1, IL2, ILG, LAY)
C
      DO 300 K = LEV1, LAY
      DO 300 I = IL1, IL2
        TAUG(I,K) =  TAUG(I,K) + CL3CH4(IG) * RMCH4 * DP(I,K)
  300 CONTINUE
C
      GW =  GW3(IG)
C
      ELSE IF (IB .EQ. 4)                                           THEN
C
C----------------------------------------------------------------------C
C     BAND3 (1100 - 1400 CM^-1), OVERLAPPING ABSORPTION OF H2O, N2O,   C
C     CH4 AND CFC12. DIRECT MAPPING METHOD FOR H2O AND CH4 AND N2O     C
C     CFC ARE CONSIDERED AS MINOR GASES                                C
C----------------------------------------------------------------------C
C
      NG2 =  4
      NG3 =  5
      CALL TLINE3 (TAUG, CL4H2O(1,1,IG), CL4CH4(1,1,IG), 
     1             CL4N2O(1,1,IG), Q, NG2, NG3, DP, DIP, DT, INPT, 
     2             LEV1, GH, MTL, IL1, IL2, ILG, LAY)
C
      LC =  4
      CALL TCONTL1 (TAUG, CL4CS(1,1,IG), CL4CF(1,1,IG), Q, DP, DIP, DT, 
     1              LC, INPT, MCONT, GH, IL1, IL2, ILG, LAY)
C
C----------------------------------------------------------------------C
C     SIMPLY ADD THE CFC EFFECT                                        C
C----------------------------------------------------------------------C
C
        DO 400 K = LEV1, LAY
        DO 400 I = IL1, IL2
          TAUG(I,K) =  TAUG(I,K) + (CL4F12(IG) * RMF12 + 
     1                 1037.3 * RMF113 + 1426.9 * RMF114) * DP(I,K)
  400   CONTINUE
C
      GW =  GW4(IG)
C
      ELSE IF (IB .EQ. 5)                                           THEN
C
C----------------------------------------------------------------------C
C     BAND5 (980 - 1100 CM^-1), OVERLAPPING ABSORPTION OF H2O AND O3   C
C     DIRECT MAPPING METHOD. CO2 AND CFC ARE SIMPLY ADDED              C
C----------------------------------------------------------------------C
C
      NG2 =  2
      CALL TLINE2 (TAUG, CL5H2O(1,1,IG), CL5O3(1,1,IG), Q, O3, NG2, DP,
     1             DIP, DT, INPT, LEV1, GH, MTL, IL1, IL2, ILG, LAY)
C
      LC =  4
      CALL TCONTL1 (TAUG, CL5CS(1,1,IG), CL5CF(1,1,IG), Q, DP, DIP, DT, 
     1              LC, INPT, MCONT, GH, IL1, IL2, ILG, LAY)
C
C----------------------------------------------------------------------C
C     SIMPLY ADD THE CO2 + CFC EFFECT                                  C
C     SINCE THE INTERACTION OF CO2 AND H2O, Q(I,K) APPEARS IN CO2      C
C     EFFECT                                                           C
C----------------------------------------------------------------------C
C
      DO 500 K = LEV1, LAY
      DO 500 I = IL1, IL2
        CRMCO2    =  2.3056E-04 * SQRT(RMCO2 * 1.E+04)
        TAUG(I,K) =  TAUG(I,K) + ( (0.009 +  0.093 * Q(I,K) / (Q(I,K) +
     1               2.1 * CRMCO2)) * RMCO2 + CL5F11(IG) * RMF11 +
     2               CL5F12(IG) * RMF12 + 
     3               1687.4 * RMF113 + 2924.1 * RMF114 ) * DP(I,K)
  500 CONTINUE
C
      GW =  GW5(IG)
C
      ELSE IF (IB .EQ. 6)                                           THEN
C
C----------------------------------------------------------------------C
C     BAND (800 - 980 CM^-1), NONGRAY GASEOUS ABSORPTION OF H2O.       C
C     + CFC11 AND CFC12                                                C
C----------------------------------------------------------------------C
C
      NG2 =  7
      NG3 =  8
      CALL TLINE3 (TAUG, CL6H2O(1,1,IG), CL6F11(1,1,IG), CL6F12(1,1,IG),
     1             Q, NG2, NG3, DP, DIP, DT, INPT, LEV1, GH, MTL, 
     2             IL1, IL2, ILG, LAY)
C
C----------------------------------------------------------------------C
C     SIMPLY ADD THE CO2 + CFC EFFECT                                  C
C----------------------------------------------------------------------C
C
        DO 600 K = LEV1, LAY
        DO 600 I = IL1, IL2
          TAUG(I,K) =  TAUG(I,K) + ( (0.0074 + 0.0396 * Q(I,K) /
     1                (Q(I,K) + 2.8 * RMCO2)) * RMCO2 +
     2                 1191.3 * RMF113 + 1098.2 * RMF114 ) * DP(I,K)
  600   CONTINUE
C
        IF (IG .EQ. 1)                                              THEN
          LC =  4
          CALL TCONTL1 (TAUG, CL6CS, CL6CF, Q, DP, DIP, DT, LC, INPT,
     1                  MCONT, GH, IL1, IL2, ILG, LAY)
C
        ENDIF
      GW =  GW6(IG)
C
      ELSE IF (IB .EQ. 7)                                           THEN
C
C----------------------------------------------------------------------C
C     BAND6 (540 - 800 CM^-1), OVERLAPPING ABSORPTION OF H2O AND CO2   C
C     EXACT MAPPING METHOD FOR H2O AND CO2, DIRECT MAPPING FOR N2O     C
C     O3 EFFECT IS SIMPLY ADDED                                        C
C----------------------------------------------------------------------C
C
      CALL TLINEHC (TAUG, CL7H2OU(1,1,IG), CL7H2OD(1,1,1,IG),
     1              CL7CO2U(1,1,IG), CL7CO2D(1,1,1,IG), Q, DP, DIP,
     2              DIR, DT, INPTR, INPT, LEV1, IL1, IL2, ILG, LAY)
C
      CALL TCONTHL (TAUG, CL7CS(1,1,1,IG), CL7CF(1,1,1,IG), Q, DP, DIP,
     1              DIR, DT, INPTR, INPT, MCONT, IL1, IL2, ILG, LAY)
C
C----------------------------------------------------------------------C
C     SIMPLY ADD THE O3 EFFECT                                         C
C----------------------------------------------------------------------C
C
      DO 700 K = LEV1, LAY
      DO 700 I = IL1, IL2
        CRMN2O    =  2.0437 * 1.E-07 * SQRT(RMN2O * 1.E+07)
        TAUG(I,K) =  TAUG(I,K) + (CL7O3(IG) * O3(I,K) + 
     1               CL7N2O(IG) * CRMN2O) * DP(I,K)
  700 CONTINUE
C
      GW =  GW7(IG)
C
      ELSE IF (IB .EQ. 8)                                           THEN
C
C----------------------------------------------------------------------C
C     BAND (340 - 540 CM^-1), NONGRAY GASEOUS ABSORPTION OF H2O.       C
C----------------------------------------------------------------------C
C
      NG =  1 
      INIT = 2
      CALL TLINE1 (TAUG, CL8H2O(1,1,IG), Q, NG, DP, DIP, DT, INPT, 
     1             LEV1, GH, MTL, INIT, IL1, IL2, ILG, LAY)
C
      IF (IG .LE. 4)                                                THEN
        LC =  6
        CALL TCONTL1 (TAUG, CL8CS(1,1,IG), CL8CF(1,1,IG), Q, DP, DIP,
     1                DT, LC, INPT, MCONT, GH, IL1, IL2, ILG, LAY)
      ENDIF
C
      GW =  GW8(IG)
C
      ELSE IF (IB .EQ. 9)                                           THEN
C
C----------------------------------------------------------------------C
C     BAND (0 - 340 CM^-1), NONGRAY GASEOUS ABSORPTION OF H2O.         C
C----------------------------------------------------------------------C
C
      NG =  1 
      INIT = 2
      CALL TLINE1 (TAUG, CL9H2O(1,1,IG), Q, NG, DP, DIP, DT, INPT,
     1             LEV1, GH, MTL, INIT, IL1, IL2, ILG, LAY)
C
      LC =  6
      CALL TCONTL1 (TAUG, CL9CS(1,1,IG), CL9CF(1,1,IG), Q, DP, DIP, DT,
     1              LC, INPT, MCONT, GH, IL1, IL2, ILG, LAY)
C
      GW =  GW9(IG)
C
      ENDIF
C
      RETURN
      END

      SUBROUTINE GASOPTLGH (TAUG, GWGH, DP, IB, IG, O3, Q, INPT, MCONT, 
     1                      DIR, DIP, DT, LEV1, GH, IL1, IL2, ILG, LAY)
C
C----------------------------------------------------------------------C
C     THE SAME AS GASOPTL BUT FOR INTERVALS CLOSE TO 1 IN THE          C
C     ACCUMULATED PROBABILITY SPACE                                    C
C     1 = H2O                                                          C
C     2 = O3                                                           C
C     3 = CO2                                                          C
C     4 = CH4                                                          C
C     5 = N2O                                                          C
C     6 = O2                                                           C
C     7 = CFC11                                                        C
C     8 = CFC12                                                        C
C     TLINE, ETC., DEAL WITH LINE ABSORPTION AND TCONTL AND TCONTHL    C
C     DEAL WITH WATER VAPOR CONTINUUM                                  C
C                                                                      C
C     TAUG: GASEOUS OPTICAL DEPTH                                      C
C     DP:   AIR MASS PATH FOR A MODEL LAYER (EXLAINED IN RADDRIV).     C
C     O3:   O3 MASS MIXING RATIO                                       C
C     Q:    WATER VAPOR MASS MIXING RATIO                              C
C     DIP:  INTERPRETATION FACTOR FOR PRESSURE BETWEEN TWO             C
C           NEIGHBORING STANDARD INPUT DATA PRESSURE LEVELS            C
C     DT:   LAYER TEMPERATURE - 250 K                                  C
C     INPT: NUMBER OF THE LEVEL FOR THE STANDARD INPUT DATA PRESSURES  C
C----------------------------------------------------------------------C
      IMPLICIT REAL (A-H,O-Z),
     +INTEGER (I-N)
C
      REAL TAUG(ILG,LAY)
C
      REAL DP(ILG,LAY), O3(ILG,LAY), Q(ILG,LAY), DIR(ILG,LAY), 
     1     DIP(ILG,LAY), DT(ILG,LAY)
      INTEGER INPT(ILG,LAY)
      LOGICAL GH
C
      COMMON /BANDL1GH/  GW1GH(5),  CL1CO2GH(5,28,5)
      COMMON /BANDL2GH/  GW2GH(1),  CL2H2OGH(5,28), 
     1                              CL2CSGH(5,3), CL2CFGH(5,3)
      COMMON /BANDL3GH/  GW3GH(3),  CL3H2OGH(5,28,3),
     1                              CL3CSGH(5,4), CL3CFGH(5,4) 
      COMMON /BANDL4GH/  GW4GH(5),  CL4H2OGH(5,28,5), CL4N2OGH(5,28,5),
     1                              CL4CH4GH(5,28,5)
      COMMON /BANDL5GH/  GW5GH(4),  CL5H2OGH(5,28,4), CL5O3GH(5,28,4),
     1                              CL5CSGH(5,4,2), CL5CFGH(5,4,2)
C     (NO BAND 6 FOR GH)
      COMMON /BANDL7GH/  GW7GH(7),  CL7H2OGH(5,28,4), CL7CO2GH(5,28,7),
     1                              CL7O3GH(2)
      COMMON /BANDL8GH/  GW8GH(3),  CL8H2OGH(5,28,3)
      COMMON /BANDL9GH/  GW9GH(6),  CL9H2OGH(5,28,6)
      DATA NTL /28/
C=======================================================================
      IF (IB .EQ. 1)                                                THEN
C
C----------------------------------------------------------------------C
C     BAND (2500 - 2200 CM^-1), NONGRAY GASEOUS ABSORPTION OF CO2.     C
C----------------------------------------------------------------------C
C
      NG =  3
      INIT = 2
      CALL TLINE1 (TAUG, CL1CO2GH(1,1,IG), Q, NG, DP, DIP, DT, INPT, 
     1             LEV1, GH, NTL, INIT, IL1, IL2, ILG, LAY)
C
      GWGH =  GW1GH(IG)
C
      ELSE IF (IB .EQ. 2)                                           THEN
C
C----------------------------------------------------------------------C
C     BAND (2200 - 1900 CM^-1), NONGRAY GASEOUS ABSORPTION OF H2O      C
C----------------------------------------------------------------------C
C
      NG =  1
      INIT = 2
      CALL TLINE1 (TAUG, CL2H2OGH, Q, NG, DP, DIP, DT, INPT, LEV1, GH,
     1             NTL, INIT, IL1, IL2, ILG, LAY)
C
      LC =  3
      CALL TCONTL1 (TAUG, CL2CSGH, CL2CFGH, Q, DP, DIP, DT, LC, INPT, 
     1              MCONT, GH, IL1, IL2, ILG, LAY)
C
      GWGH =  GW2GH(IG)
C
      ELSE IF (IB .EQ. 3)                                           THEN
C
C----------------------------------------------------------------------C
C     BAND (1900 - 1400 CM^-1), NONGRAY GASEOUS ABSORPTION OF H2O.     C
C----------------------------------------------------------------------C
C
      NG =  1
      INIT = 2
      CALL TLINE1 (TAUG, CL3H2OGH(1,1,IG), Q, NG, DP, DIP, DT, INPT,
     1             LEV1, GH, NTL, INIT, IL1, IL2, ILG, LAY)
C
      IF (IG .EQ. 1)                                                THEN
        LC =  4
        CALL TCONTL1 (TAUG, CL3CSGH, CL3CFGH, Q, DP, DIP, DT, LC, INPT,
     1                MCONT, GH, IL1, IL2, ILG, LAY)
C
      ENDIF
C
      GWGH =  GW3GH(IG)
C
      ELSE IF (IB .EQ. 4)                                           THEN
C
C----------------------------------------------------------------------C
C     BAND3 (1100 - 1400 CM^-1), OVERLAPPING ABSORPTION OF H2O, N2O,   C
C     AND CH4. DIRECT MAPPING METHOD FOR H2O AND CH4 AND N2O           C
C----------------------------------------------------------------------C
C
      NG2 =  4
      NG3 =  5 
      CALL TLINE3 (TAUG, CL4H2OGH(1,1,IG), CL4CH4GH(1,1,IG), 
     1             CL4N2OGH(1,1,IG), Q, NG2, NG3, DP, DIP, DT, INPT,
     2             LEV1, GH, NTL, IL1, IL2, ILG, LAY)
C
      GWGH =  GW4GH(IG)
C
      ELSE IF (IB .EQ. 5)                                           THEN
C
C----------------------------------------------------------------------C
C     BAND5 (980 - 1100 CM^-1), OVERLAPPING ABSORPTION OF H2O AND O3   C
C     DIRECT MAPPING METHOD                                            C
C----------------------------------------------------------------------C
C
      NG2 =  2
      CALL TLINE2 (TAUG, CL5H2OGH(1,1,IG), CL5O3GH(1,1,IG), Q, O3, NG2,
     1             DP, DIP, DT, INPT, LEV1, GH, NTL, IL1, IL2, ILG, LAY)
C
      IF (IG .LE. 2)                                                THEN
        LC =  4
        CALL TCONTL1 (TAUG, CL5CSGH(1,1,IG), CL5CFGH(1,1,IG), Q, DP,
     1                DIP, DT, LC, INPT, MCONT, GH, IL1, IL2, ILG, LAY)
      ENDIF
C
      GWGH =  GW5GH(IG)
C
C----------------------------------------------------------------------C
C     BAND (800 - 980 CM^-1), NO GH                                    C
C----------------------------------------------------------------------C
C
      ELSE IF (IB .EQ. 7)                                           THEN
C
C----------------------------------------------------------------------C
C     BAND6 (540 - 800 CM^-1), OVERLAPPING ABSORPTION OF H2O AND CO2   C
C     DIRECT MAPPING METHOD. FOR IG > 4, THE CONTRIBUTION BY H2O IS    C
C     VERY SMALL.                                                      C
C----------------------------------------------------------------------C
C
      IF (IG .LE. 4)                                                THEN
        NG2 =  3
        CALL TLINE2 (TAUG, CL7H2OGH(1,1,IG), CL7CO2GH(1,1,IG), Q, O3, 
     1               NG2, DP, DIP, DT, INPT, LEV1, GH, NTL,  
     2               IL1, IL2, ILG, LAY)

C----------------------------------------------------------------------C
C     SIMPLY ADD THE O3 EFFECT                                         C
C----------------------------------------------------------------------C
C
        IF (IG .LE. 2)                                              THEN
          DO 700 K = 1, LAY
          DO 700 I = IL1, IL2
            TAUG(I,K) =  TAUG(I,K) + CL7O3GH(IG) * O3(I,K) * DP(I,K)
  700     CONTINUE
        ENDIF
      ELSE
C
        NG =  3
        INIT = 2
        CALL TLINE1 (TAUG, CL7CO2GH(1,1,IG), Q, NG, DP, DIP, DT, INPT,
     1               LEV1, GH, NTL, INIT, IL1, IL2, ILG, LAY)
      ENDIF
C
      GWGH =  GW7GH(IG)
C
      ELSE IF (IB .EQ. 8)                                           THEN
C
C----------------------------------------------------------------------C
C     BAND (340 - 540 CM^-1), NONGRAY GASEOUS ABSORPTION OF H2O.       C
C----------------------------------------------------------------------C
C
      NG =  1
      INIT = 2
      CALL TLINE1 (TAUG, CL8H2OGH(1,1,IG), Q, NG, DP, DIP, DT, INPT,
     1             LEV1, GH, NTL, INIT, IL1, IL2, ILG, LAY)
C
      GWGH =  GW8GH(IG)
C
      ELSE IF (IB .EQ. 9)                                           THEN
C
C----------------------------------------------------------------------C
C     BAND (0 - 340 CM^-1), NONGRAY GASEOUS ABSORPTION OF H2O.         C
C----------------------------------------------------------------------C
C
      NG =  1
      INIT = 2
      CALL TLINE1 (TAUG, CL9H2OGH(1,1,IG), Q, NG, DP, DIP, DT, INPT,
     1             LEV1, GH, NTL, INIT, IL1, IL2, ILG, LAY)
C
      GWGH =  GW9GH(IG)
C
      ENDIF
C
      RETURN
      END

      SUBROUTINE GASOPTS (TAUG, GW, DP, IB, IG, O3, Q, INPT, DIP, DT,
     1                    RMU3, LEV1, GH, IL1, IL2, ILG, LAY)
C
C----------------------------------------------------------------------C
C     CALCULATION OF THE OPTICAL DEPTHS DUE TO NONGRAY GASEOUS         C
C     ABSORPTION FOR THE SOLAR, IN EACH LAYER FOR A GIVEN BAND IB AND  C
C     CUMULATIVE PROBABILITY GW.                                       C
C     RELATIVE SOLAR ENERGY IN EACH SOLAR BAND ARE                     C
C     BAND 1:   622.8483                                               C
C     BAND 1GH:   7.5917                                               C
C     BAND 2:   430.0919                                               C
C     BAND 2GH:   8.9036                                               C
C     BAND 3:   238.6979                                               C
C     BAND 3GH:   7.4453                                               C
C     BAND 4:    33.4129                                               C
C     BAND 4GH:   7.0384                                               C
C                                                                      C
C     TOTAL RELATIVE SOLAR ENERGY IN FROM 0.2 - 4 UM IS                C
C     1356.0300 W / M^2, PLUS 11.9096 W / M^2 IN 4 - 10 UM.            C
C     TOTAL  1367.9396 W / M^2                                         C
C                                                                      C
C     MINOR GAS:                                                       C
C     3 = CO2                                                          C
C     4 = CH4                                                          C
C     6 = O2                                                           C
C     THIS SUBROUTINE ONLY CALCULATES TAUG BELOW 1 MB                  C
C                                                                      C
C     TAUG:  GASEOUS OPTICAL DEPTH                                     C
C     DP:    AIR MASS PATH FOR A MODEL LAYER (EXLAINED IN RADDRIV).    C
C     O3:    O3 MASS MIXING RATIO                                      C
C     Q:     WATER VAPOR MASS MIXING RATIO                             C
C     DIP:   INTERPRETATION FACTOR FOR PRESSURE BETWEEN TWO            C
C            NEIGHBORING STANDARD INPUT DATA PRESSURE LEVELS           C
C     DT:    LAYER TEMPERATURE - 250 K                                 C
C     INPT:  NUMBER OF THE LEVEL FOR THE STANDARD INPUT DATA PRESSURES C
C     RMU3:  A FACTOR OF SOLAR ZENITH ANGLE, GIVEN IN RADDRIV          C
C----------------------------------------------------------------------C
      IMPLICIT REAL (A-H,O-Z),
     +INTEGER (I-N)
C
      REAL TAUG(ILG,LAY)
C
      REAL DP(ILG,LAY), O3(ILG,LAY), Q(ILG,LAY), DIP(ILG,LAY),
     1     DT(ILG,LAY), RMU3(ILG), S(ILG,LAY)
      INTEGER INPT(ILG,LAY)
      LOGICAL GH
C
      COMMON /BANDS1/ GW1(6), CS1O3(3,6), CS1O21
      COMMON /BANDS2/ GW2(4), CS2H2O(5,18,4), CS2O2(5,18,2), CS2O3(4)
      COMMON /BANDS3/ GW3(6), CS3H2O(5,18,6), CS3CO2(5,18,6),
     1                        CS3CH4(5,18,2)
      COMMON /BANDS4/ GW4(4), CS4H2O(5,18,4), CS4CO2(5,18,4),
     1                        CS4CH4(5,18,2)
C
      COMMON /TRACE / RMCO2, RMCH4, RMN2O, RMO2, RMF11, RMF12, RMF113,
     1                RMF114
C
C     * NUMBER OF VERTICAL LEVELS IN ABSORBER PRESSURE-BASED COEFFICIENT
C     * ARRAY ("M" REFERENCES NON-SATURATED BANDS ACTIVE BELOW 1 MB ONLY).
C
      DATA MTL /18/
C=======================================================================
      IF (IB .EQ. 1)                                                THEN
C
C----------------------------------------------------------------------C
C     BAND (14500 - 50000 CM^-1), NONGRAY GASEOUS ABSORPTION OF O3,    C
C     H2O AND O2.                                                      C
C     RELATIVE SOLAR ENERGY 630.4401 WM^-2.                            C
C     IG9 (50000-43000)  UVC                           1.21100 (W M^-2)C
C     IG8 (43000-37500)  UVC                           3.17570         C
C     IG7 (37500-35700)  UVC                           3.20501         C
C     UVC ALL INCLUDED IN GH PART                                      C
C                                                                      C
C     IG6 (35700-34200)  UVB                           4.73084         C
C     IG5 (34200-32185)  UVB                          10.14919         C
C     IG4 (32185-31250)  UVB  J VALUE: 32185 CM^-1     6.70594         C
C                                                                      C
C     IG3 (31250-25000)  UVA                          83.43346         C
C                                                                      C
C     IG2 (25000-19000)  PAR                         236.97212         C
C     IG1 (19000-14500)  PAR                         280.85678         C
C     PAR: PHOTOSYNTHETIC ACTIVE RADIATION                             C
C     NOTE THE SPECTRAL STRUCTURE IS SLIGHTLY DIFF FROM LI & BARKER    C
C     (2005 JAS)                                                       C
C                                                                      C
C     THE EFFECT OF H2O AND O2 IS ADDED WITH SIMPLE METHOD             C
C----------------------------------------------------------------------C
C
        IF (IG .EQ. 1)                                              THEN
          DO 110 K = LEV1, LAY
          DO 110 I = IL1, IL2
            IF (INPT(1,K) .LT. 950)                                 THEN
              M =  INPT(I,K)
            ELSE 
              M =  INPT(I,K) - 1000
            ENDIF
C
            IF (M .LT. 7)                                           THEN
              X       = (CS1O21 - 0.881E-05 * RMU3(I)) * RMO2
            ELSE
              X       = (0.108E-04 - 0.881E-05 * RMU3(I)) * RMO2
            ENDIF
C
            IF (M .LT. 15)                                          THEN
              X       =  X + (0.199E-02 - 0.952E-03 * RMU3(I)) * Q(I,K) 
            ELSE
              X       =  X + (0.208E-02 - 0.952E-03 * RMU3(I)) * Q(I,K)
            ENDIF
C
            DTO3      =  DT(I,K) + 23.13
            TAUG(I,K) = ((CS1O3(1,IG) + DTO3 * (CS1O3(2,IG) +
     1                   DTO3 * CS1O3(3,IG))) * O3(I,K) + X) * DP(I,K)
  110     CONTINUE
        ELSE
          DO 120 K = LEV1, LAY
          DO 120 I = IL1, IL2
            DTO3      =  DT(I,K) + 23.13
            TAUG(I,K) = (CS1O3(1,IG) + DTO3 * (CS1O3(2,IG) +
     1                   DTO3 * CS1O3(3,IG))) * O3(I,K) * DP(I,K)
  120     CONTINUE
        ENDIF
C
        GW =  GW1(IG)
C
      ELSE IF (IB .EQ. 2)                                           THEN
C
C----------------------------------------------------------------------C
C     BAND (8400 - 14500 CM^-1), NONGRAY GASEOUS ABSORPTION OF H2O,    C
C     O2 AND O3                                                        C
C     RELATIVE SOLAR ENERGY 430.0919 W M^-2                            C
C----------------------------------------------------------------------C
C
        IF (IG .LE. 2)                                              THEN
          NG2 =  6
          CALL TLINE2 (TAUG, CS2H2O(1,1,IG), CS2O2(1,1,IG), Q, O3, NG2,
     1                 DP, DIP, DT, INPT, LEV1, GH, MTL, 
     2                 IL1, IL2, ILG, LAY)
        ELSE 
          NG =  1
          INIT = 2
          CALL TLINE1 (TAUG, CS2H2O(1,1,IG), Q, NG, DP, DIP, DT, INPT,
     1                 LEV1, GH, MTL, INIT, IL1, IL2, ILG, LAY)
        ENDIF
C
C----------------------------------------------------------------------C
C     SIMPLY ADD O3 EFFECT                                             C
C----------------------------------------------------------------------C
C
        DO 200 K = LEV1, LAY
        DO 200 I = IL1, IL2
          TAUG(I,K)   =  TAUG(I,K) + CS2O3(IG) * O3(I,K) * DP(I,K)
  200   CONTINUE
C
        GW =  GW2(IG)
C
      ELSE IF (IB .EQ. 3)                                           THEN
C
C----------------------------------------------------------------------C
C     BAND (4200 - 8400 CM^-1), NONGRAY GASEOUS ABSORPTION OF H2O, CO2 C
C     AND CH4                                                          C
C     RELATIVE SOLAR ENERGY 238.6979 W M^-2                            C
C----------------------------------------------------------------------C
C
        IF (IG .LE. 2)                                              THEN
          NG2 =  3
          NG3 =  4
          CALL TLINE3 (TAUG, CS3H2O(1,1,IG), CS3CO2(1,1,IG),
     1                 CS3CH4(1,1,IG), Q, NG2, NG3, DP, DIP, DT, INPT,
     2                 LEV1, GH, MTL, IL1, IL2, ILG, LAY)
        ELSE
          NG2 =  3
          CALL TLINE2 (TAUG, CS3H2O(1,1,IG), CS3CO2(1,1,IG), Q, O3, 
     1                 NG2, DP, DIP, DT, INPT, LEV1, GH, MTL, 
     2                 IL1, IL2, ILG, LAY)
        ENDIF
C
        GW =  GW3(IG)
C
      ELSE IF (IB .EQ. 4)                                           THEN
C
C----------------------------------------------------------------------C
C     BAND (2500 - 4200 CM^-1), NONGRAY GASEOUS ABSORPTION OF H2O      C
C     AND CO2                                                          C
C     RELATIVE SOLAR ENERGY 33.4129 W M^-2                             C
C----------------------------------------------------------------------C
C
        IF (IG .LE. 2)                                              THEN
          NG2 =  3
          NG3 =  4
          CALL TLINE3 (TAUG, CS4H2O(1,1,IG), CS4CO2(1,1,IG),
     1                 CS4CH4(1,1,IG), Q, NG2, NG3, DP, DIP, DT, INPT,
     2                 LEV1, GH, MTL, IL1, IL2, ILG, LAY)
        ELSE
          NG2 =  3
          CALL TLINE2 (TAUG, CS4H2O(1,1,IG), CS4CO2(1,1,IG), Q, O3, 
     1                 NG2, DP, DIP, DT, INPT, LEV1, GH, MTL, 
     2                 IL1, IL2, ILG, LAY)
        ENDIF
C
        GW =  GW4(IG)
C
      ENDIF
C
      RETURN
      END

      SUBROUTINE LATTENU (ATTEN, IB, IG, O3, Q, DP, DIP, DT, DT0, INPT,
     1                    IL1, IL2, ILG)
C
C----------------------------------------------------------------------C
C     CALCULATION OF THE ATTENUATION FOR THE DOWNWARD FLUX ABOVE THE   C
C     MODEL TOP LEVEL. SINCE THE TEMPERATURE AT 0.005 MB IS UNKNOWN WE C
C     ASSUME IT IS THE SAME AS THAT OF MODEL TOP LEVEL                 C
C                                                                      C
C     ATTEN: FOR SOLAR: THE ATTENUATION FACTOR FOR DOWNWARD FLUX FROM  C
C            TOA TO THE MODEL TOP LEVEL; FOR LONGWAVE: THE OPTICAL     C
C            / DIFFUSE FACTOR                                          C
C     DP:    HERE DP IS ONLY THE PRESSURE DIFFERENCE, DIFFERENT FROM   C
C            THAT DEFINED IN RADDRIV. SO THERE IS A FACTOR 1.02        C
C     O3:    O3 MASS MIXING RATIO                                      C
C     Q:     WATER VAPOR MASS MIXING RATIO                             C
C     DIP:   INTERPRETATION FACTOR FOR PRESSURE BETWEEN TWO            C
C            NEIGHBORING STANDARD INPUT DATA PRESSURE LEVELS           C
C     DT:    LAYER TEMPERATURE - 250 K                                 C
C     DT0:   TEMPERATURE IN MOON LAYER - 250 K                         C
C----------------------------------------------------------------------C
      IMPLICIT REAL (A-H,O-Z),
     +INTEGER (I-N)
C
      REAL ATTEN(ILG), O3(ILG), Q(ILG), DP(ILG), DIP(ILG), DT(ILG), 
     1     DT0(ILG), RMU(ILG)
      INTEGER INPT(ILG)
C
      COMMON /BANDL1GH/  GW1GH(5),  CL1CO2GH(5,28,5)
      COMMON /BANDL3GH/  GW3GH(3),  CL3H2OGH(5,28,3), CL3CSGH(5,4), 
     1                              CL3CFGH(5,4)
      COMMON /BANDL5GH/  GW5GH(4),  CL5H2OGH(5,28,4), CL5O3GH(5,28,4),
     1                              CL5CSGH(5,4,2), CL5CFGH(5,4,2)
      COMMON /BANDL7GH/  GW7GH(7),  CL7H2OGH(5,28,4), CL7CO2GH(5,28,7),
     1                              CL7O3GH(2)
      COMMON /BANDL8GH/  GW8GH(3),  CL8H2OGH(5,28,3)
      COMMON /BANDL9GH/  GW9GH(6),  CL9H2OGH(5,28,6)
C
C     * NUMBER OF VERTICAL LEVELS IN ABSORBER PRESSURE-BASED COEFFICIENT
C     * ARRAY.
C
      DATA NTL /28/
C=======================================================================
      IF (IB .EQ. 1)                                                THEN
        NG = 3
        ISL = 2
        CALL ATTENUE (ATTEN, CL1CO2GH(1,1,IG), O3, Q, DP, DIP, DT, DT0, 
     1                RMU, INPT, NTL, NG, ISL, IL1, IL2, ILG)
C
      ELSE IF (IB .EQ. 2)                                           THEN
        DO 200 I = IL1, IL2
          ATTEN(I) =  0.0  
  200   CONTINUE
C
      ELSE IF (IB .EQ. 3)                                           THEN
        NG = 1
        ISL = 2
        CALL ATTENUE (ATTEN, CL3H2OGH(1,1,IG), O3, Q, DP, DIP, DT, DT0, 
     1                RMU, INPT, NTL, NG, ISL, IL1, IL2, ILG)
C
      ELSE IF (IB .EQ. 4)                                           THEN
        DO 400 I = IL1, IL2
          ATTEN(I) =  0.0
  400   CONTINUE
C
      ELSE IF (IB .EQ. 5)                                           THEN
        NG = 2
        ISL = 2
        CALL ATTENUE (ATTEN, CL5O3GH(1,1,IG), O3, Q, DP, DIP, DT, DT0, 
     1                RMU, INPT, NTL, NG, ISL, IL1, IL2, ILG)
C
      ELSE IF (IB .EQ. 7)                                           THEN
        NG = 3
        ISL = 2
        CALL ATTENUE (ATTEN, CL7CO2GH(1,1,IG), O3, Q, DP, DIP, DT, DT0, 
     1                RMU, INPT, NTL, NG, ISL, IL1, IL2, ILG)
C
      ELSE IF (IB .EQ. 8)                                           THEN
        NG = 1
        ISL = 2
        CALL ATTENUE (ATTEN, CL8H2OGH(1,1,IG), O3, Q, DP, DIP, DT, DT0, 
     1                RMU, INPT, NTL, NG, ISL, IL1, IL2, ILG)
C
      ELSE IF (IB .EQ. 9)                                           THEN
        NG = 1
        ISL = 2
        CALL ATTENUE (ATTEN, CL9H2OGH(1,1,IG), O3, Q, DP, DIP, DT, DT0,
     1                RMU, INPT, NTL, NG, ISL, IL1, IL2, ILG)
C
      ENDIF
C
      RETURN
      END

      SUBROUTINE PLANCK (BF, BS, URBF, BF0, URBF0, DBF, TFULL, GT, IB,
     1                   IL1, IL2, ILG, LAY, LEV)
C
C----------------------------------------------------------------------C
C     CALCULATION OF PLANCK FUNCTION IN VALID RANGE 120 - 360 K        C
C                                                                      C
C     BF:    BLACKBODY INTENSITY INTEGRATED OVER EACH BAND AT EACH     C
C            LEVEL IN UNITS W / M^2 / SR.                              C
C     BS:    THE BLACKBODY INTENSITY AT THE SURFACE.                   C
C     BF0:   THE BLACKBODY INTENSITY AT THE TOA (ASSUME 210 K).        C
C     TFULL: TEMPERATURE AT EACH LEVEL                                 C
C     GT:    TEMPERATURE AT GROUND                                     C
C     U:     1 / DIFFUSIVITY FACTOR                                    C
C     URBF:  U TIMES THE DIFFERENCE OF LOG(BF) FOR TWO NEIGHBOR LEVELS C
C            USED FOR EXPONENTIAL SOURCE FUNCTION (LI, 2002 JAS P3302) C
C     URBF0: U TIMES THE DIFFERENCE OF LOG(BF) FOR TOA AND LEVEL 1     C
C     DBF:   DIFFERENCE OF BF FOR TWO NEIGHBOR LEVELS USED FOR LINEAR  C
C            SOURCE FUNCTION (LI, 2002 JAS P3302)                      C
C                                                                      C
C     0.0040816327 = 1 / 245 (245 THE STANDARD TEMPERATURE FOR POLY.   C
C     FIT)                                                             C
C----------------------------------------------------------------------C
      IMPLICIT REAL (A-H,O-Z),
     +INTEGER (I-N)
C
      REAL  BF(ILG,LEV), BS(ILG), BF0(ILG), URBF(ILG,LAY), URBF0(ILG),
     1      DBF(ILG,LAY), TFULL(ILG,LEV), GT(ILG), XX(ILG,LAY)
      REAL  XP(6,9)
      DATA  U / 0.60653066 /, RTSTAND / 0.0040816327 /
      DATA ((XP(I,J), I = 1, 6), J = 1, 9) /
     1  -2.9876423E+00,    1.3660089E+01,   -1.2944461E+01,
     1   1.1775748E+01,   -1.9236798E+01,    2.3584435E+01,
     2  -1.6414103E+00,    1.1898535E+01,   -1.1262182E+01,
     2   1.0236863E+01,   -1.6677772E+01,    2.0423136E+01,
     3   6.5215205E-01,    9.2657366E+00,   -8.5872301E+00,
     3   7.6765044E+00,   -1.2287254E+01,    1.4990547E+01,
     4   1.5442143E+00,    7.2253228E+00,   -6.7811515E+00,
     4   6.1572299E+00,   -9.8725011E+00,    1.1997278E+01,
     5   1.2777580E+00,    6.1257638E+00,   -5.7906013E+00,
     5   5.3296782E+00,   -8.7529282E+00,    1.0741367E+01,
     6   2.1005257E+00,    5.2376301E+00,   -4.8915631E+00,
     6   4.5030997E+00,   -7.3199981E+00,    8.9204038E+00,
     7   2.9091223E+00,    3.9860795E+00,   -3.5829565E+00,
     7   3.2692193E+00,   -5.1799711E+00,    6.2157752E+00,
     8   2.7856424E+00,    2.8179582E+00,   -2.3780464E+00,
     8   2.1432949E+00,   -3.4540206E+00,    4.1814100E+00,
     9   2.4623332E+00,    1.8731841E+00,   -1.3659538E+00,
     9   1.1484948E+00,   -1.5975564E+00,    1.7791135E+00 /
C
      DO 100 I = IL1, IL2
        DT          =  GT(I) * RTSTAND - 1.0
        BS(I)       =  EXP( XP(1,IB) +
     1                            DT * (XP(2,IB) + DT * (XP(3,IB) +
     2                            DT * (XP(4,IB) + DT * (XP(5,IB) +
     3                            DT *  XP(6,IB) )))) )
C
        DT          = (2. * TFULL(I,1) - TFULL(I,2)) * RTSTAND - 1.0
        XXT         =  XP(1,IB) + DT * (XP(2,IB) + DT * (XP(3,IB) +
     2                            DT * (XP(4,IB) + DT * (XP(5,IB) +
     3                            DT *  XP(6,IB) ))))
C
        DT          =  TFULL(I,1) * RTSTAND - 1.0
        XX0         =  XP(1,IB) + DT * (XP(2,IB) + DT * (XP(3,IB) +
     1                            DT * (XP(4,IB) + DT * (XP(5,IB) +
     2                            DT *  XP(6,IB) ))))
        DT          =  TFULL(I,2) * RTSTAND - 1.
        XX(I,1)     =  XP(1,IB) + DT * (XP(2,IB) + DT * (XP(3,IB) +
     1                            DT * (XP(4,IB) + DT * (XP(5,IB) +
     2                            DT *  XP(6,IB) ))))
C
        BF0(I)      =  EXP(XXT)
        URBF0(I)    =  U * (XX0 - XXT)
        BF(I,1)     =  EXP(XX0)
        BF(I,2)     =  EXP(XX(I,1))
        DBF(I,1)    =  BF(I,2) - BF(I,1)
        URBF(I,1)   =  U * (XX(I,1) - XX0)
  100 CONTINUE
C
      DO 205 K = 2, LAY
        KM1 = K - 1
        KP1 = K + 1
        DO 200 I = IL1, IL2
          DT        =  TFULL(I,KP1) * RTSTAND - 1.0
          XX(I,K)   =  XP(1,IB) + DT * (XP(2,IB) + DT * (XP(3,IB) +
     1                            DT * (XP(4,IB) + DT * (XP(5,IB) +
     2                            DT *  XP(6,IB) ))))
C
          BF(I,KP1) =  EXP(XX(I,K))
          DBF(I,K)  =  BF(I,KP1) - BF(I,K)
          URBF(I,K) =  U * (XX(I,K) - XX(I,KM1))
  200     CONTINUE
  205   CONTINUE
C
      RETURN
      END

      SUBROUTINE PREINTP (INPT, INPTM, DIP, DIP0, P, IL1, IL2, ILG,
     1                    LAY)
C
C----------------------------------------------------------------------C
C     THIS SUBROUTINE DETERMINES THE PRESSURE INTERPRETATION POINTS    C
C                                                                      C
C     INPT:  NUMBER OF THE LEVEL FOR THE STANDARD INPUT DATA PRESSURES C
C            (FOR 28 INTERPRETATION LEVELS)                            C
C     INPTM: NUMBER OF THE LEVEL FOR THE STANDARD INPUT DATA PRESSURES C
C            (FOR 18 INTERPRETATION LEVELS BELOW 1 MB)                 C
C     P:     PRESSURE AT MIDDLE OF EACH LAYER                          C
C     DIP:   INTERPRETATION FACTOR FOR PRESSURE BETWEEN TWO            C
C            NEIGHBORING STANDARD INPUT DATA PRESSURE LEVELS           C
C     DIP0:  INTERPRETATION FACTOR FOR PRESSURE ABOVE MODEL TOP LEVEL  C
C----------------------------------------------------------------------C
      IMPLICIT REAL (A-H,O-Z),
     +INTEGER (I-N)
C
      REAL DIP(ILG,LAY), DIP0(ILG)
      REAL P(ILG,LAY), STANDP(28)
      INTEGER INPT(ILG,LAY), INPTM(ILG,LAY)
C
      DATA STANDP / 5.0000E-04, 1.4604E-03, 2.9621E-03, 6.0080E-03,
     1              1.2186E-02, 2.4717E-02, 5.0134E-02, 1.0169E-01,
     2              2.0625E-01, 4.1834E-01,
     3              1.2180, 1.8075, 2.6824, 3.9806, 5.9072, 8.7662,
     4              13.0091, 19.3054, 28.6491, 42.5151, 63.0922,
     5              93.6284, 138.9440, 206.1920, 305.9876, 454.0837,
     6              673.8573, 1000.0000 /
C
      JENDS = 27
      DO 600 K = 1, LAY
        DO 100 I = IL1, IL2
          INPT(I,K)   =  0
  100   CONTINUE
C
        DO 200 J = 1, JENDS
        DO 200 I = IL1, IL2
          IF (P(I,K) .GT. STANDP(J))                                THEN
            INPT(I,K) =  INPT(I,K) + 1
          ENDIF
  200   CONTINUE
C
C----------------------------------------------------------------------C
C     CALCULATE ARRAYS DIP AND DIT REQUIRED LATER FOR GASOPT ROUTINES. C
C     ALSO, SET VALUES OF INPT FOR A GIVEN LEVEL TO BE NEGATIVE IF ALL C
C     LONGITUDE VALUES ARE THE SAME. THIS IS ALSO USED IN THE GASOPT   C
C     ROUTINES TO IMPROVE PERFORMANCE BY ELIMINATING THE UNNECESSARY   C
C     INDIRECT-ADDRESSING IF INPL IS NEGATIVE FOR A GIVEN LEVEL.       C
C     NOTE THAT FOR INPT=0, IT IS ASSUMED THAT LEVELS ARE MORE OR      C
C     LESS HORIZONTAL IN PRESSURE, SO SCALING BY -1 STILL PRESERVES    C
C     THE VALUE OF ZERO AND NO INDIRECT-ADDRESSING IS DONE IN THE      C
C     GASOPT ROUTINES.                                                 C
C----------------------------------------------------------------------C
C
        INPDIF =  0
        INP1   =  INPT(1,K)
        DO 300 I = IL1, IL2
          IF (INPT(I,K) .NE. INP1)  INPDIF = 1
          M  =  INPT(I,K)
          N  =  M + 1
          IF (M .GT. 0)                                             THEN
            DIP(I,K)  = (P(I,K) - STANDP(M)) / (STANDP(N) - STANDP(M))
          ELSE
            DIP(I,K)  =  P(I,K) / STANDP(1)
          ENDIF
  300   CONTINUE
C
        IF (INPDIF .EQ. 0)                                          THEN
          DO 400 I = IL1, IL2
            INPT(I,K) =  INPT(I,K) + 1000
  400     CONTINUE
        ENDIF
C
        DO 500 I = IL1, IL2
          INPTM(I,K)  =  INPT(I,K) - 10
  500   CONTINUE
  600 CONTINUE
C
C----------------------------------------------------------------------C
C     INTERPRETATION FACTOR FOR LATTENU AND SATTENU (ATTENUATION ABOVE C
C     MODEL TOP                                                        C
C----------------------------------------------------------------------C
C
      PM =  P(1,1)  
      DO 700 I = IL1, IL2
        PM          =  MIN (PM, P(I,1))
  700 CONTINUE
C
      IF (PM .LE. 0.0005)                                           THEN
        DO 800 I = IL1, IL2
          DIP0(I)   =  0.0
  800   CONTINUE
      ELSE
        DO 900 I = IL1, IL2
          P0        =  P(I,1) * P(I,1) / P(I,2)
          X         =  SQRT (P0 * P(I,1))
          DIP0(I)   = (X - P(I,1)) / (P0 - P(I,1))
  900   CONTINUE
      ENDIF
C
      RETURN
      END

      SUBROUTINE PREINTR (INPR, DIR, Q, IL1, IL2, ILG, LAY)
C
C----------------------------------------------------------------------C
C     THIS SUBROUTINE DETERMINES THE INTERPRETION POINTS FOR THE RATIO C
C     OF H2O AND CO2.                                                  C
C                                                                      C
C     INPR:  NUMBER OF THE RATIO LEVEL FOR THE STANDARD 5 RATIOS       C
C     DIR:   INTERPRETATION FACTOR FOR MASS RATIO OF H2O / CO2         C
C            BETWEEN TWO NEIGHBORING STANDARD INPUT RATIOS             C
C     Q:     WATER VAPOR MASS MIXING RATIO                             C
C     RHC:  THE RATIO OF THE H2O MASS MIXING TO CO2 MASS MIXING        C
C----------------------------------------------------------------------C
      IMPLICIT REAL (A-H,O-Z),
     +INTEGER (I-N)
C
      REAL DIR(ILG,LAY), Q(ILG,LAY), RHC(ILG,LAY), STANDR(5)
      INTEGER INPR(ILG,LAY)
C
      COMMON /TRACE / RMCO2, RMCH4, RMN2O, RMO2, RMF11, RMF12, RMF113,
     1                RMF114
C
      DATA STANDR / .06,  .24,  1., 4., 16. /
C=======================================================================
      RRMCO2 =  1.0 / (RMCO2 + 1.E-10)
      JENDR  =  5
      DO 400 K = 1, LAY
C
        DO 100 I = IL1, IL2
          INPR(I,K)   =  0
          RHC(I,K)    =  Q(I,K) * RRMCO2
  100   CONTINUE
C
        DO 200 J = 1, JENDR
        DO 200 I = IL1, IL2
          IF (RHC(I,K) .GT. STANDR(J))                              THEN
            INPR(I,K) =  INPR(I,K) + 1
          ENDIF
  200   CONTINUE
C
        DO 300 I = IL1, IL2
          L   =  INPR(I,K)
          LP1 =  L + 1
          IF (L .GE. 1 .AND. L .LT. 5)                              THEN
            DIR(I,K)  = (RHC(I,K) - STANDR(L)) / 
     1                  (STANDR(LP1) - STANDR(L))
          ELSE
C
C----------------------------------------------------------------------C
C     DIR IS NOT USED WITH VALUES OF {0,5} IN TLINEHC, BUT WE          C
C     INITIALIZE HERE TO AVOID PROBLEMS WITH NaN WHEN USED             C
C     IN MULTITASKING MODE.                                            C
C----------------------------------------------------------------------C
C
            DIR(I,K)  =  0.0
          ENDIF
  300   CONTINUE
C
  400 CONTINUE
C
      RETURN
      END

      SUBROUTINE RAYLEI (TAUR, IB, DP, IL1, IL2, ILG, LAY)
C
C----------------------------------------------------------------------C
C     RAYLEIGH SCATTERING FOR BANDS2-BANDS4, NEAR INFRARED REGION      C
C                                                                      C
C     TAUR: RAYLEIGH OPTICAL DEPTH                                     C
C     DP:   AIR MASS PATH FOR A MODEL LAYER (EXLAINED IN RADDRIV).     C
C----------------------------------------------------------------------C
      IMPLICIT REAL (A-H,O-Z),
     +INTEGER (I-N)
C
      REAL RI(3), TAUR(ILG,LAY), DP(ILG,LAY)
      DATA RI / .16305E-04, .17997E-05, .13586E-06 /
C=======================================================================
      IBM1 = IB - 1
      DO 100 K = 1, LAY
      DO 100 I = IL1, IL2
        TAUR(I,K) =  RI(IBM1) * DP(I,K)
  100 CONTINUE
C
      RETURN
      END

      SUBROUTINE RAYLEV (TAUR, IG, DP, RMU3, IL1, IL2, ILG, LAY)
C
C----------------------------------------------------------------------C
C     RAYLEIGH SCATTERING FOR EACH SUB-BAND IN BANDS1, VISIBLE REGION  C
C     TAUR IS THE OPTICAL DEPTH RAYLEIGH SCATTERING FOR A GIVEN LAYER  C
C     FOR UVC (35700 - 50000 CM^-1), SINCE THE OPTICAL DEPTH OF O3 AND C
C     O2 ARE VERY LARGE, RAYLEIGH SCATTERING EFFECT IS NEGLECTED, IT   C
C     IS SHOWN EVEN FOR 10% O3 AMOUNT OF THE STANDARD ATMO, THE        C
C     RAYLEIGH SCATTERING FOR UVC STILL CAN BE NEGLECTED.              C
C     FOR PAR AND UVA, SINCE THEIR SPECTRAL RANGES ARE VERY WIDE, SMALLC
C     ERRORS COULD OCCUR FOR LARGE ZENITH ANGLE, SLIGHTLY ADJUSTMENT   C
C     IS NEEDED, THIS DOES MEAN THE RAYLEIGH OPTICAL DEPTH IS RELATED  C
C     SOLAR ZENITH ANGLE FOR MULTIPLE SCATTERING PROCESS IN SWTRAN.    C
C                                                                      C
C     TAUR: RAYLEIGH OPTICAL DEPTH                                     C
C     DP:   AIR MASS PATH FOR A MODEL LAYER (EXLAINED IN RADDRIV).     C
C     RMU3:  A FACTOR OF SOLAR ZENITH ANGLE, GIVEN IN RADDRIV          C
C----------------------------------------------------------------------C
      IMPLICIT REAL (A-H,O-Z),
     +INTEGER (I-N)
C
      REAL RI0(6), RI2(3), TAUR(ILG,LAY), DP(ILG,LAY), RMU3(ILG)
C
      DATA RI0 / .67758E-04, .19425E-03, .52145E-03, .88082E-03,
     1           .10541E-02, .12835E-02 /
      DATA RI2 / .28000E-05, .21300E-04, .36000E-04 /
C=======================================================================
      IF (IG .LE. 3)                                                THEN
        DO 100 K = 1, LAY
        DO 100 I = IL1, IL2
          TAUR(I,K) = (RI0(IG) - RI2(IG) * RMU3(I)) * DP(I,K)
  100   CONTINUE
      ELSE
        DO 200 K = 1, LAY
        DO 200 I = IL1, IL2
          TAUR(I,K) =  RI0(IG) * DP(I,K)
  200   CONTINUE
      ENDIF
C
      RETURN
      END

      SUBROUTINE SATTENU (ATTEN, IB, IG, RMU, O3, Q, DP, DIP, DT, DT0,
     1                    INPT, GH, IL1, IL2, ILG)
C
C----------------------------------------------------------------------C
C     CALCULATION OF SOLAR ATTENUATION ABOVE THE MODEL TOP LEVEL. FOR  C
C     BAND1 ONLY O3 AND O2 ARE CONSIDERED, THE CONTRIBUTION OF OTHER   C
C     GASES IS SMALL. FOR BAND 3 AND 4, CO2 IS CONSIDERED FOR GH       C
C                                                                      C
C     ATTEN: ATTENUATION FACTOR FOR DOWNWARD FLUX FROM TOA TO THE      C
C            MODEL TOP LEVEL                                           C
C     O3:    O3 MASS MIXING RATIO                                      C
C     Q:     WATER VAPOR MASS MIXING RATIO                             C
C     DP:    HERE DP IS ONLY THE PRESSURE DIFFERENCE, DIFFERENT FROM   C
C            THAT DEFINED IN RADDRIV. SO THERE IS A FACTOR 1.02        C
C     DIP:   INTERPRETATION FACTOR FOR PRESSURE                        C
C     DT:    LAYER TEMPERATURE - 250 K                                 C
C     DT0:   TEMPERATURE IN MOON LAYER - 250 K                         C
C----------------------------------------------------------------------C
      IMPLICIT REAL (A-H,O-Z),
     +INTEGER (I-N)
C
      REAL ATTEN(ILG), RMU(ILG), O3(ILG), Q(ILG), DP(ILG), DIP(ILG),
     1     DT(ILG), DT0(ILG)
      INTEGER INPT(ILG)
      LOGICAL GH
C
      COMMON /TRACE / RMCO2, RMCH4, RMN2O, RMO2, RMF11, RMF12, RMF113,
     1                RMF114
      COMMON /BANDS1/ GWS1(6), CS1O3(3,6), CS1O21
C
      COMMON /BANDS1GH/ GWS1GH(3), CS1O3GH(3,3),     CS1O2GH3
      COMMON /BANDS2GH/ GWS2GH(4), CS2H2OGH(5,28),   CS2O2GH(5,28,3)
      COMMON /BANDS3GH/ GWS3GH(4), CS3H2OGH(5,28,2), CS3CO2GH(5,28,4)
      COMMON /BANDS4GH/ GWS4GH(9), CS4H2OGH(5,28,6), CS4CO2GH(5,28,6),
     1                             CS4CH4GH(5,28)
C
C     * NUMBER OF VERTICAL LEVELS IN ABSORBER PRESSURE-BASED COEFFICIENT
C     * ARRAY.
C
      DATA NTL /28/
C=======================================================================
      IF (IB .EQ. 1)                                                THEN
        IF (GH)                                                     THEN
          IF (IG .EQ. 3)                                            THEN
            DO 100 I = IL1, IL2
              DTO3     =  DT(I) + 23.13
              TAU      =  1.02 * ((CS1O3GH(1,IG) + DTO3 *
     1                   (CS1O3GH(2,IG) + DTO3 * CS1O3GH(3,IG))) *
     2                    O3(I) +  CS1O2GH3 * RMO2) * DP(I)
              ATTEN(I) =  EXP( - TAU / RMU(I))
  100       CONTINUE
          ELSE
            DO 110 I = IL1, IL2
              DTO3     =  DT(I) + 23.13
              TAU      =  1.02 * (CS1O3GH(1,IG) + DTO3 *
     1                   (CS1O3GH(2,IG) + DTO3 * CS1O3GH(3,IG))) *
     2                    O3(I) * DP(I)
              ATTEN(I) =  EXP( - TAU / RMU(I))
  110       CONTINUE
          ENDIF
        ELSE
C
          IF (IG .EQ. 1)                                            THEN
            DO 120 I = IL1, IL2
              DTO3     =  DT(I) + 23.13
              TAU      =  1.02 * ((CS1O3(1,IG) + DTO3 * (CS1O3(2,IG) +
     1                    DTO3 * CS1O3(3,IG))) * O3(I) +
     2                    CS1O21 * RMO2) * DP(I)
              ATTEN(I) =  EXP( - TAU / RMU(I))
  120       CONTINUE
          ELSE
            DO 130 I = IL1, IL2
              DTO3     =  DT(I) + 23.13
              TAU      =  1.02 * (CS1O3(1,IG) + DTO3 * (CS1O3(2,IG) +
     1                    DTO3 * CS1O3(3,IG))) * O3(I) * DP(I)
              ATTEN(I) =  EXP( - TAU / RMU(I))
  130       CONTINUE
          ENDIF
        ENDIF
C
      ELSE IF (IB .EQ. 2)                                           THEN
        IF (IG .EQ. 1)                                              THEN
          DO 200 I = IL1, IL2
            ATTEN(I)   =  1.0
  200     CONTINUE
        ELSE
          NG = 6
          IM = IG - 1
          ISL = 1
          CALL ATTENUE (ATTEN, CS2O2GH(1,1,IM), O3, Q, DP, DIP, DT,
     1                   DT0, RMU, INPT, NTL, NG, ISL, IL1, IL2, ILG)
        ENDIF
C
      ELSE IF (IB .EQ. 3)                                           THEN
        NG = 3
        ISL = 1
        CALL ATTENUE (ATTEN, CS3CO2GH(1,1,IG), O3, Q, DP, DIP, DT, DT0, 
     1                 RMU, INPT, NTL, NG, ISL, IL1, IL2, ILG)
C
      ELSE IF (IB .EQ. 4)                                           THEN
        NG = 3
        IF (IG .NE. 4 .AND. IG .NE. 6 .AND. IG .NE. 8)              THEN
          IF (IG .LE. 3)  IM =  IG 
          IF (IG .EQ. 5)  IM =  IG - 1
          IF (IG .EQ. 7)  IM =  IG - 2
          IF (IG .EQ. 9)  IM =  IG - 3
          ISL = 1
          CALL ATTENUE (ATTEN, CS4CO2GH(1,1,IM), O3, Q, DP, DIP, DT, 
     1                  DT0, RMU, INPT, NTL, NG, ISL, IL1, IL2, ILG)
        ELSE IF (IG .EQ. 4)                                         THEN
          NG = 4
          CALL ATTENUE (ATTEN, CS4CH4GH, O3, Q, DP, DIP, DT, DT0, RMU,
     1                  INPT, NTL, NG, 1, IL1, IL2, ILG)
        ELSE
          DO 400 I = IL1, IL2
            ATTEN(I)   =  1.0
  400     CONTINUE
        ENDIF
      ENDIF
C
      RETURN
      END

      SUBROUTINE STRANDNGH (TRAN, GWGH, ATTEN, TAUA, TAUOMA, TAUCS, 
     1                      TAUOMC, CLD, RMU, DP, O3, Q, IB, IG, INPT,
     2                      DIP, DT, LEV1, GH, CUT, 
     3                      IL1, IL2, ILG, LAY, LEV)
C
C----------------------------------------------------------------------C
C     CALCULATION OF THE DOWNWARD SOLAR FLUX UNDER THE CONDITION THAT  C
C     THE EXTINCTION COEFFICIENT OF GAS IS VERY LARGE, THE SCATTERING  C
C     EFFECTS CAN BE NEGLECTED. THE CLOUD OPTICAL DEPTH IS MUCH SMALLERC
C     THAN THE GASEOUS OPTICAL DEPTH, THE CLOUD EFFECT IS VERY SMALL   C
C     AND BE TREATED SIMPLY                                            C
C                                                                      C
C     HRS:   SOLAR HEATING RATE (K / SEC)                              C
C     TRAN:  DOWNWARD FLUX                                             C
C     ATTEN: ATTENUATION FACTOR FOR DOWNWARD FLUX FROM TOA TO THE      C
C            MODEL TOP LEVEL                                           C
C     TAUCS: CLOUD OPTICAL DEPTH                                       C
C     CLD:   CLOUD FRACTION                                            C
C     RMU:   COS OF SOLAR ZENITH ANGLE                                 C
C     DP:    AIR MASS PATH FOR A MODEL LAYER (EXLAINED IN RADDRIV).    C
C     O3:    O3 MASS MIXING RATIO                                      C
C     Q:     WATER VAPOR MASS MIXING RATIO                             C
C     INPT:  NUMBER OF THE LEVEL FOR THE STANDARD INPUT DATA PRESSURES C
C     DIP:   INTERPRETATION FACTOR FOR PRESSURE BETWEEN TWO            C
C            NEIGHBORING STANDARD INPUT DATA PRESSURE LEVELS           C
C     DT:    LAYER TEMPERATURE - 250 K                                 C
C----------------------------------------------------------------------C
      IMPLICIT REAL (A-H,O-Z),
     +INTEGER (I-N)
C
      REAL TRAN(ILG,2,LEV), ATTEN(ILG), TAUA(ILG,LAY), TAUOMA(ILG,LAY),
     1     TAUCS(ILG,LAY), TAUOMC(ILG,LAY), CLD(ILG,LAY),RMU(ILG), 
     2     DP(ILG,LAY), O3(ILG,LAY), Q(ILG,LAY), DIP(ILG,LAY), 
     3     DT(ILG,LAY), TAUG(ILG,LAY)
      INTEGER INPT(ILG,LAY)
      LOGICAL GH
C
      COMMON /BANDS1GH/ GWS1GH(3), CS1O3GH(3,3),     CS1O2GH3
      COMMON /BANDS2GH/ GWS2GH(4), CS2H2OGH(5,28),   CS2O2GH(5,28,3)
      COMMON /BANDS3GH/ GWS3GH(4), CS3H2OGH(5,28,2), CS3CO2GH(5,28,4)
      COMMON /BANDS4GH/ GWS4GH(9), CS4H2OGH(5,28,6), CS4CO2GH(5,28,6),
     1                             CS4CH4GH(5,28)
C
      COMMON /TRACE / RMCO2, RMCH4, RMN2O, RMO2, RMF11, RMF12, RMF113,
     1                RMF114
C
C     * NUMBER OF VERTICAL LEVELS IN ABSORBER PRESSURE-BASED COEFFICIENT
C     * ARRAY.
C
      DATA NTL /28/
C=======================================================================
      DO 10 I = IL1, IL2
        TRAN(I,1,1)           =  ATTEN(I)
        TRAN(I,2,1)           =  ATTEN(I)
   10 CONTINUE
C
      IF (IB .EQ. 1)                                                THEN
C
C----------------------------------------------------------------------C
C     BAND1 FOR UVC (35700 - 50000 CM^-1), NONGRAY GASEOUS ABSORPTION  C
C     OF O2  AND O3. SOLAR ENERGY  7.5917 W M^-2                       C
C----------------------------------------------------------------------C
C
        IF (IG .EQ. 3)                                              THEN
          DO 105 K = 1, LAY
            KP1 = K + 1
            DO 100 I = IL1, IL2
              DTO3            =  DT(I,K) + 23.13
              TAU             = ((CS1O3GH(1,IG) +
     1                           DTO3 * (CS1O3GH(2,IG) +
     2                           DTO3 * CS1O3GH(3,IG))) * O3(I,K) +
     3                           CS1O2GH3 * RMO2) * DP(I,K) + TAUA(I,K)
              DTR1            =  EXP( - (TAU - TAUOMA(I,K)) / RMU(I))
              TRAN(I,1,KP1)   =  TRAN(I,1,K) * DTR1
C
              IF (CLD(I,K) .LT. CUT)                                THEN
                TRAN(I,2,KP1) =  TRAN(I,2,K) * DTR1
              ELSE
                ABSC          = (1.0 - CLD(I,K)) * DTR1 + CLD(I,K) *
     1                           EXP( - (TAU + TAUCS(I,K) - TAUOMC(I,K))
     2                           / RMU(I))
                TRAN(I,2,KP1) =  TRAN(I,2,K) * ABSC
              ENDIF
  100       CONTINUE
  105     CONTINUE
        ELSE
          DO 115 K = 1, LAY
            KP1 = K + 1
            DO 110 I = IL1, IL2
              DTO3            =  DT(I,K) + 23.13
              TAU             = (CS1O3GH(1,IG) + DTO3 * (CS1O3GH(2,IG) +
     1                           DTO3 * CS1O3GH(3,IG))) * O3(I,K) *
     2                           DP(I,K) + TAUA(I,K)
              DTR1            =  EXP( - (TAU - TAUOMA(I,K)) / RMU(I))
              TRAN(I,1,KP1)   =  TRAN(I,1,K) * DTR1
C
              IF (CLD(I,K) .LT. CUT)                                THEN
                TRAN(I,2,KP1) =  TRAN(I,2,K) * DTR1
              ELSE
                ABSC          = (1.0 - CLD(I,K)) * DTR1 + CLD(I,K) *
     1                           EXP( - (TAU + TAUCS(I,K) - TAUOMC(I,K))
     2                           / RMU(I))
                TRAN(I,2,KP1) =  TRAN(I,2,K) * ABSC
              ENDIF
  110       CONTINUE
  115     CONTINUE
        ENDIF
      GWGH =  GWS1GH(IG)
C
      ELSE IF (IB .EQ. 2)                                           THEN
C
C----------------------------------------------------------------------C
C     BAND (8400 - 14500 CM^-1), NONGRAY GASEOUS ABSORPTION OF O2      C
C     AND O3. SOLAR ENERGY  8.9036 W M^-2                              C
C----------------------------------------------------------------------C
C
        IF (IG .EQ. 1)                                              THEN
          NG =  1
          INIT = 2
          CALL TLINE1 (TAUG, CS2H2OGH(1,1), Q, NG, DP, DIP, DT, INPT,
     1                 LEV1, GH, NTL, INIT, IL1, IL2, ILG, LAY)
        ELSE
          IM =  IG - 1
          NG =  6
          INIT = 2
          CALL TLINE1 (TAUG, CS2O2GH(1,1,IM), Q, NG, DP, DIP, DT, INPT,
     1                 LEV1, GH, NTL, INIT, IL1, IL2, ILG, LAY)
        ENDIF
C
        DO 205 K = 1, LAY
          KP1 = K + 1
          DO 200 I = IL1, IL2
            TAU               =  TAUG(I,K) + TAUA(I,K)
            DTR1              =  EXP( - (TAU - TAUOMA(I,K)) / RMU(I))
            TRAN(I,1,KP1)     =  TRAN(I,1,K) * DTR1
C
            IF (CLD(I,K) .LT. CUT)                                  THEN
              TRAN(I,2,KP1)   =  TRAN(I,2,K) * DTR1
            ELSE
              ABSC            = (1.0 - CLD(I,K)) * DTR1 + CLD(I,K) *
     1                           EXP( - (TAU + TAUCS(I,K) - TAUOMC(I,K))
     2                           / RMU(I))
              TRAN(I,2,KP1)   =  TRAN(I,2,K) * ABSC
            ENDIF
  200     CONTINUE
  205   CONTINUE
C
        GWGH =  GWS2GH(IG)
C
      ELSE IF (IB .EQ. 3)                                           THEN
C
C----------------------------------------------------------------------C
C     BAND (4200 - 8400 CM^-1), NONGRAY GASEOUS ABSORPTION OF H2O AND  C
C     CO2. SOLAR ENERGY 7.4453 W M^-2                                  C
C----------------------------------------------------------------------C
C
        IF (IG .LE. 2)                                              THEN
          NG2 =  3
          CALL TLINE2 (TAUG, CS3H2OGH(1,1,IG), CS3CO2GH(1,1,IG), Q, O3,
     1                 NG2, DP, DIP, DT, INPT, LEV1, GH, NTL, 
     2                 IL1, IL2, ILG, LAY)

        ELSE
          NG =  3
          INIT = 2
          CALL TLINE1 (TAUG, CS3CO2GH(1,1,IG), Q, NG, DP, DIP, DT, 
     1                 INPT, LEV1, GH, NTL, INIT, IL1, IL2, ILG, LAY)
        ENDIF
C
        DO 305 K = 1, LAY
          KP1 = K + 1
          DO 300 I = IL1, IL2
            TAU               =  TAUG(I,K) + TAUA(I,K)
            DTR1              =  EXP( - (TAU - TAUOMA(I,K)) / RMU(I))
            TRAN(I,1,KP1)     =  TRAN(I,1,K) * DTR1
C
            IF (CLD(I,K) .LT. CUT)                                  THEN
              TRAN(I,2,KP1)   =  TRAN(I,2,K) * DTR1
            ELSE
              ABSC            = (1.0 - CLD(I,K)) * DTR1 + CLD(I,K) *
     1                           EXP( - (TAU + TAUCS(I,K) - TAUOMC(I,K))
     2                           / RMU(I))
              TRAN(I,2,KP1)   =  TRAN(I,2,K) * ABSC
            ENDIF
  300     CONTINUE
  305   CONTINUE
C
        GWGH =  GWS3GH(IG)
C
      ELSE IF (IB .EQ. 4)                                           THEN
C
C----------------------------------------------------------------------C
C     BAND (2500 - 4200 CM^-1), NONGRAY GASEOUS ABSORPTION OF H2O      C
C     AND CO2. SOLAR ENERGY 7.0384 W M^-2                              C
C----------------------------------------------------------------------C
C
        IF (IG .LE. 3)                                              THEN
          NG2 =  3
          CALL TLINE2 (TAUG, CS4H2OGH(1,1,IG), CS4CO2GH(1,1,IG), Q, O3,
     1                 NG2, DP, DIP, DT, INPT, LEV1, GH, NTL,
     2                 IL1, IL2, ILG, LAY)
        ELSE IF (IG .EQ. 6 .OR. IG .EQ. 8)                          THEN
          NG =  1
          IF (IG .EQ. 6)  IM = 5
          IF (IG .EQ. 8)  IM = 6
          INIT = 2
          CALL TLINE1 (TAUG, CS4H2OGH(1,1,IM), Q, NG, DP, DIP, DT, 
     1                 INPT, LEV1, GH, NTL, INIT, IL1, IL2, ILG, LAY)
        ELSE IF (IG .EQ. 5 .OR. IG .EQ. 7 .OR. IG .EQ. 9)           THEN
          NG =  3
          IF (IG .EQ. 5)  IM = 4
          IF (IG .EQ. 7)  IM = 5
          IF (IG .EQ. 9)  IM = 6
          INIT = 2
          CALL TLINE1 (TAUG, CS4CO2GH(1,1,IM), Q, NG, DP, DIP, DT, 
     1                 INPT, LEV1, GH, NTL, INIT, IL1, IL2, ILG, LAY)
        ELSE
          NG2 =  4
          IM = 4
          CALL TLINE2 (TAUG, CS4H2OGH(1,1,IM), CS4CH4GH, Q, O3,
     1                 NG2, DP, DIP, DT, INPT, LEV1, GH, NTL,
     2                 IL1, IL2, ILG, LAY)
        ENDIF
C
        DO 405 K = 1, LAY
          KP1 = K + 1
          DO 400 I = IL1, IL2
            TAU               =  TAUG(I,K) + TAUA(I,K)
            DTR1              =  EXP( - (TAU - TAUOMA(I,K)) / RMU(I))
            TRAN(I,1,KP1)     =  TRAN(I,1,K) * DTR1
C
            IF (CLD(I,K) .LT. CUT)                                  THEN
              TRAN(I,2,KP1)   =  TRAN(I,2,K) * DTR1
            ELSE
              ABSC            = (1.0 - CLD(I,K)) * DTR1 + CLD(I,K) *
     1                           EXP( - (TAU + TAUCS(I,K) - TAUOMC(I,K))
     2                           / RMU(I))
              TRAN(I,2,KP1)   =  TRAN(I,2,K) * ABSC
            ENDIF
  400     CONTINUE
  405   CONTINUE
C
        GWGH =  GWS4GH(IG)
C
      ENDIF
C
      RETURN
      END

      SUBROUTINE STRANDN (TRAN, ATTN, ATTNTOP, RMU, DP, DT, O3, RMU3,
     1                    IB, IG, LEV1, IL1, IL2, ILG, LAY, LEV)
C
C----------------------------------------------------------------------C
C     CALCULATION OF THE DOWNWARD FLUX FROM TOP LEVEL TO 1 MB, NO      C
C     SCATTERING EFFECT IS CONSIDERED                                  C
C                                                                      C
C     TRAN:    TRANSMISIVITY                                           C
C     ATTN:    ATTENUATION FACTOR FOR REDUCING SOLAR FLUX FROM MODEL   C
C              TOP LEVEL TO 1 MB                                       C
C     ATTNTOP: ATTENUATION FACTOR FOR REDUCING SOLAR FLUX FROM TOA TO  C
C              MODEL TOP LEVEL                                         C
C     RMU:     COSINE OF SOLAR ZENITH ANGLE                            C
C     DP:      AIR MASS PATH FOR A MODEL LAYER (EXLAINED IN RADDRIV).  C
C     O3:      O3 MASS MIXING RATIO                                    C
C----------------------------------------------------------------------C
      IMPLICIT REAL (A-H,O-Z),
     +INTEGER (I-N)
C
      REAL TRAN(ILG,2,LEV), ATTN(ILG), ATTNTOP(ILG), RMU(ILG), 
     1     DP(ILG,LAY), DT(ILG,LAY), O3(ILG,LAY), RMU3(ILG)
C
      COMMON /TRACE / RMCO2, RMCH4, RMN2O, RMO2, RMF11, RMF12, RMF113,
     1                RMF114
      COMMON /BANDS1/ GW1(6), CS1O3(3,6), CS1O21
C=======================================================================
      LEV1M1 =  LEV1 - 1
C
      IF (IB .EQ. 1)                                                THEN
C
C----------------------------------------------------------------------C
C     BAND (14500 - 50000 CM^-1), NONGRAY GASEOUS ABSORPTION OF O3     C
C     ANS O2.                                                          C
C----------------------------------------------------------------------C
C
        DO 10 I = IL1, IL2
          TRAN(I,1,1)       =  ATTNTOP(I)
          TRAN(I,2,1)       =  ATTNTOP(I)
   10   CONTINUE
C
        IF (IG .EQ. 1)                                              THEN
          DO 150 K = 1, LEV1M1
            KP1 = K + 1
            DO 100 I = IL1, IL2
              DTO3          =  DT(I,K) + 23.13
              X             = (CS1O21 - 0.881E-05 * RMU3(I)) * RMO2
              TAU           = ((CS1O3(1,IG) + DTO3 * (CS1O3(2,IG) +
     1                         DTO3 * CS1O3(3,IG))) * O3(I,K) + X) *
     2                         DP(I,K)
              TRAN(I,1,KP1) =  TRAN(I,1,K) * EXP( - TAU / RMU(I))
              TRAN(I,2,KP1) =  TRAN(I,1,KP1)
  100       CONTINUE
  150     CONTINUE
C
        ELSE 
          DO 250 K = 1, LEV1M1
            KP1 = K + 1
            DO 200 I = IL1, IL2
              DTO3          =  DT(I,K) + 23.13
              TAU           = (CS1O3(1,IG) + DTO3 * (CS1O3(2,IG) +
     1                         DTO3 * CS1O3(3,IG))) * O3(I,K) * DP(I,K)
              TRAN(I,1,KP1) =  TRAN(I,1,K) * EXP( - TAU / RMU(I))
              TRAN(I,2,KP1) =  TRAN(I,1,KP1)
  200       CONTINUE
  250     CONTINUE
        ENDIF
C
C----------------------------------------------------------------------C
C     FLUX ADJUSTMENT FOR REGION BELOW 1 MB                            C
C----------------------------------------------------------------------C
C
        DO 400 I = IL1, IL2
          ATTN(I)           =  TRAN(I,1,LEV1)
  400   CONTINUE
C
      ELSE
C
        DO 500 K = 1, LEV1
        DO 500 I = IL1, IL2
          TRAN(I,1,K)       =  1.0
          TRAN(I,2,K)       =  1.0
  500   CONTINUE
C
        DO 600 I = IL1, IL2
          ATTN(I)           =  1.0
  600   CONTINUE
      ENDIF
C
      RETURN
      END

      SUBROUTINE STRANUP (REFL, DP, DT, O3, IB, IG, LEV1, 
     1                    IL1, IL2, ILG, LAY, LEV)
C
C----------------------------------------------------------------------C
C     CALCULATION OF THE UPWARD SOLAR FLUX ABOVE 1 MB, NO SCATTERING   C
C     EFFECT IS CONSIDERED                                             C
C                                                                      C
C     REFL:  REFLECTIVITY                                              C
C     TAU:   OPTICAL DEPTH                                             C
C     DTR:   DIRECT TRANSMISSION FUNCTION                              C
C     1.6487213 DIFFUSIVITY FACTOR (LI, 2000 JAS P753)                 C
C----------------------------------------------------------------------C
C
      IMPLICIT REAL (A-H,O-Z),
     +INTEGER (I-N)
C
      REAL REFL(ILG,2,LEV), DP(ILG,LAY), DT(ILG,LAY), O3(ILG,LAY)
C
      COMMON /TRACE / RMCO2, RMCH4, RMN2O, RMO2, RMF11, RMF12, RMF113,
     1                RMF114
      COMMON /BANDS1/ GWS1(6), CS1O3(3,6), CS1O21
C=======================================================================
C
      LEV1M1 =  LEV1 - 1
C
      IF (IB .EQ. 1)                                                THEN
C
        IF (IG .EQ. 1)                                              THEN
          DO 105 K = LEV1M1, 1, - 1
            KP1 = K + 1
            DO 100 I = IL1, IL2
              DTO3        =  DT(I,K) + 23.13
              TAU         = ((CS1O3(1,IG) + DTO3 * (CS1O3(2,IG) +
     1                       DTO3 * CS1O3(3,IG))) * O3(I,K) + 
     2                       CS1O21 * RMO2) * DP(I,K)
              DTR         =  EXP( - 1.6487213 * TAU)
              REFL(I,1,K) =  REFL(I,1,KP1) * DTR
              REFL(I,2,K) =  REFL(I,2,KP1) * DTR
  100       CONTINUE
  105     CONTINUE
C
        ELSE
          DO 115 K = LEV1M1, 1, - 1
            KP1 = K + 1
            DO 110 I = IL1, IL2
              DTO3        =  DT(I,K) + 23.13
              TAU         = (CS1O3(1,IG) + DTO3 * (CS1O3(2,IG) +
     1                       DTO3 * CS1O3(3,IG))) * O3(I,K) * DP(I,K)
              DTR         =  EXP( - 1.6487213 * TAU)
              REFL(I,1,K) =  REFL(I,1,KP1) * DTR
              REFL(I,2,K) =  REFL(I,2,KP1) * DTR
  110       CONTINUE
  115     CONTINUE
C
        ENDIF
C
      ELSE
C
        DO 250 K = LEV1M1, 1, - 1
          KP1 = K + 1
          DO 200 I = IL1, IL2
            REFL(I,1,K)   =  REFL(I,1,KP1)
            REFL(I,2,K)   =  REFL(I,2,KP1)
  200     CONTINUE
  250   CONTINUE
C
      ENDIF
C
      RETURN
      END

      SUBROUTINE SWTRAN (REFL, TRAN, CUMDTR, TRAN0, TAUA, TAUR, TAUG,
     1                   TAUOMA, TAUOMGA, F1, F2, TAUCS, TAUOMC, 
     2                   TAUOMGC, CLD, CLDM, A1, RMU, C1, C2, ALBSUR,
     3                   CSALB, NBLK, NCT,II,JJ,IAER, 
     4                   CUT, LEV1, IL1, IL2, ILG, LAY, LEV)
C
C----------------------------------------------------------------------C
C     DELTA-EDDINGTON APPROXIMATION AND ADDING PROCESS FOR CLEAR AND   C
C     ALL SKY, THE ADDING METHOD BY COAKLEY ET AL (1983). THIS CODE    C
C     CAN DEAL WITH SOLAR RADIATIVE TRANSFER THROUGH ATMOSPHERE WITH   C
C     PROPER TREATMENT OF CLOUD OVERLAP (RANDOM + MAXIMUM OR RANDOM +  C
C     SLANTWISE) AND CLOUD SUB-GRID VARIABILITY. THE THEORY FOR ADDING,C
C     CLOUD OVERLAP LI AND DOBBIE (2003). CLOUD SUB-GRID VARIABILITY   C
C     SIMILAR TO WITH ADJUSTMENT OF CLOUD OPTICAL DEPTH                C
C                                                                      C
C     REFL:    REFLECTIVITY (1) CLEAR SKY; (2) ALL SKY                 C
C     TRAN:    TRANSMITIVITY                                           C
C     CUMDTR:  DIRECT TRANSMISSION FOR MULT-LAYERS                     C
C     TAUA:    AEROSOL OPTICAL DEPTH                                   C
C     TAUR:    RAYLEIGH OPTICAL DEPTH                                  C
C     TAUG:    GASEOUS OPTICAL DEPTH                                   C
C     TAUOMA:  AEROSOL OPTICAL DEPTH TIMES AEROSOL SINGLE SCATTERING   C
C              ALBEDO                                                  C
C     TAUOMGA: TAUOMA TIMES AEROSOL ASYMMETRY FACTOR                   C
C     F1:      SQUARE OF AEROSOL ASYMMETRY FACTOR                      C
C     F2:      SQUARE OF CLOUD ASYMMETRY FACTOR                        C
C     TAUCS:   CLOUD OPTICAL DEPTH                                     C
C     TAUOMC:  CLOUD OPTICAL DEPTH TIMES CLOUD SINGLE SCATTERING ALBEDOC
C     TAUOMGC: TAUOMC TIMES CLOUD ASYMMETRY FACTOR                     C
C     CLD:     CLOUD FRACTION                                          C
C     CLDM:    MAXIMUM PORTION IN EACH CLOUD BLOCK, IN WHICH THE EXACT C
C              SOLUTION FOR SUBGRID VARIABILITY IS APPLIED             C
C     A1:      VARIOUS RELATION FOR CLOUD OVERLAP                      C
C     RMU:     COS OF SOLAR ZENITH ANGLE                               C
C     C1 @ C2: TWO FACTORS NOT DEPENDENT ON IB AND IG CALCULATED       C
C              OUTSIDE FOR EFFICIENCY                                  C
C     ALBSUR:  SURFACE ALBEDO                                          C
C     CSALB:   CLEAR-SKY SURFACE ALBEDO                                C
C     NBLK:    NUMBER OF CLOUD BLOCKS ACCOUNTED FROM SURFACE           C
C     NCT:     THE HIGHEST CLOUD TOP LEVEL FOR THE LONGITUDE AND       C
C              LATITUDE LOOP (ILG)                                     C
C     RDF:     LAYER DIFFUSE REFLECTION                                C
C     TDF:     LAYER DIFFUSE TRANSMISSION                              C
C     RDR:     LAYER DIRECT REFLECTION                                 C
C     TDR:     LAYER DIRECT TRANSMISSION                               C
C     DTR:     DIRECT TRANSMISSION                                     C
C     RMDF:    BLOCK DIFFUSE REFLECTION FROM MODEL TOP LEVEL           C
C     TMDR:    BLOCK DIRECT TRANSMISSION FROM MODEL TOP LEVEL          C
C     RMUR:    BLOCK DIRECT REFLECTION FROM MODEL BOTTOM LEVEL         C
C     RMUF:    BLOCK DIFFUSE REFLECTION FROM MODEL BOTTOM LEVEL        C
C----------------------------------------------------------------------C
      IMPLICIT REAL (A-H,O-Z),
     +INTEGER (I-N)
C
      REAL REFL(ILG,2,LEV), TRAN(ILG,2,LEV), CUMDTR(ILG,4,LEV), 
     1     TRAN0(ILG)
      REAL TAUA(ILG,LAY), TAUR(ILG,LAY), TAUG(ILG,LAY), TAUOMA(ILG,LAY),
     1     TAUOMGA(ILG,LAY), F1(ILG,LAY), F2(ILG,LAY), TAUCS(ILG,LAY),
     2     TAUOMC(ILG,LAY), TAUOMGC(ILG,LAY), CLD(ILG,LAY), 
     3     CLDM(ILG,LAY), A1(ILG,11), RMU(ILG), C1(ILG), C2(ILG),
     4     ALBSUR(ILG), CSALB(ILG)
C
C mxy+
C CHANGE ARRAYS TO REAL *8 
C
      REAL*8 RDF (ILG,4,LAY), TDF (ILG,4,LAY), RDR (ILG,4,LAY),
     1     TDR (ILG,4,LAY), DTR (ILG,4,LAY), RMDF(ILG,4,LEV),
     2     TMDR(ILG,4,LEV), RMUR(ILG,4,LEV), RMUF(ILG,4,LEV)
C mxy+
      INTEGER NBLK(ILG, LAY), NCT(ILG)
C mxy+
c
c     * scalars promoted to 64-bit in loop 200 to avoid pole singularity
c     * from two-stream calculation in clear-sky.
c
      real*8 extopt,omars,ssalb,sf1,sf,tau1,om1,cow
      real*8 ssgas1,cowg,alamd,u,u2,uu,efun,efun2,x1,rn,x2
      real*8 y,yx,dm,gscw,appgm,apmgm,omarcs,sf2
      real*8 tau2,om2,ssgas2,sdtr,x3,x4,srdf,stdf,srdr,stdr
      real*8 x5,x6,x7,x,lowbnd
c
      data lowbnd /1.e-15/
C mxy+
C
C----------------------------------------------------------------------C
C     COMBINE THE OPTICAL PROPERTIES FOR SOLAR,                        C
C     1, AEROSOL + RAYLEIGH + GAS; 2, CLOUD + AEROSOL + RAYLEIGH + GAS C
C     CALCULATE THE DIRECT AND DIFFUSE REFLECTION AND TRANSMISSION IN  C
C     THE SCATTERING LAYERS USING THE DELTA-EDDINGTON METHOD.          C
C----------------------------------------------------------------------C
C
      DO 200 K = LEV1, LAY
      DO 200 I = IL1, IL2
C mxy!!
        if(TAUA(I,K).lt.0. .or. 
     &     TAUR(I,K).lt.0. .or.
     &     TAUG(I,K).lt.0) then
c         print*,'negative optical depth',II,JJ,K,
c     &    TAUA(I,K),TAUR(I,K),TAUG(I,K)
        endif
        TAUA(I,K) = MAX(TAUA(I,K),1.e-20)
        TAUR(I,K) = MAX(TAUR(I,K),1.e-20)
        TAUG(I,K) = MAX(TAUG(I,K),1.e-20)
C mxy!!
        EXTOPT                    =  TAUA(I,K) + TAUR(I,K) + TAUG(I,K)
        OMARS                     =  TAUOMA(I,K) + TAUR(I,K)
        SSALB                     =  OMARS / (EXTOPT + 1.E-20)
        SF1                       =  F1(I,K) / OMARS
        SF                        =  SSALB * SF1
        TAU1                      =  EXTOPT * (1.0 - SF)
        OM1                       = (SSALB - SF) / (1.0 - SF)
        COW                       =  1.0 - OM1 + 1.E-10
        SSGAS1                    = (TAUOMGA(I,K) / OMARS - SF1) / 
     1                              (1.0 - SF1)
        COWG                      =  1.0 - OM1 * SSGAS1
C
        DTR(I,1,K)                =  EXP( - TAU1 / RMU(I))
        ALAMD                     =  SQRT(3.0 * COW * COWG)
        U                         =  1.50 * COWG / ALAMD
        U2                        =  U + U
        UU                        =  U * U
        EFUN                      =  EXP(- ALAMD * TAU1)
        EFUN2                     =  EFUN * EFUN
        X1                        = (UU - U2 + 1.0) * EFUN2
        RN                        =  1.0 / (UU + U2 + 1.0 - X1)
        X2                        =  ALAMD * RMU(I) 
        Y                         =  1.0 - X2 * X2
        YX                        =  SIGN( MAX( ABS(Y), LOWBND), Y )
        DM                        =  OM1 / YX
        GSCW                      =  SSGAS1 * COW
        APPGM                     = (C1(I) + 0.50 + 
     1                               GSCW * (C1(I) + C2(I))) * DM
        APMGM                     = (C1(I) - 0.50 + 
     1                               GSCW * (C1(I) - C2(I))) * DM
        RDF(I,1,K)                = (UU - 1.0) * (1.0  - EFUN2) * RN 
        TDF(I,1,K)                = (U2 + U2) * EFUN * RN
        RDR(I,1,K)                =  APPGM * RDF(I,1,K) + APMGM *
     1                              (TDF(I,1,K) * DTR(I,1,K) - 1.0)
        TDR(I,1,K)                =  APPGM * TDF(I,1,K) +
     1                              (APMGM * RDF(I,1,K) - APPGM + 1.0) *
     2                               DTR(I,1,K)
C
        IF (CLD(I,K) .LT. CUT)                                      THEN
          RDF(I,2,K)              =  RDF(I,1,K)
          TDF(I,2,K)              =  TDF(I,1,K)
          RDR(I,2,K)              =  RDR(I,1,K)
          TDR(I,2,K)              =  TDR(I,1,K)
          DTR(I,2,K)              =  DTR(I,1,K)
          RDF(I,3,K)              =  RDF(I,1,K)
          TDF(I,3,K)              =  TDF(I,1,K)
          RDR(I,3,K)              =  RDR(I,1,K)
          TDR(I,3,K)              =  TDR(I,1,K)
          DTR(I,3,K)              =  DTR(I,1,K)
          RDF(I,4,K)              =  RDF(I,1,K)
          TDF(I,4,K)              =  TDF(I,1,K)
          RDR(I,4,K)              =  RDR(I,1,K)
          TDR(I,4,K)              =  TDR(I,1,K)
          DTR(I,4,K)              =  DTR(I,1,K)
        ELSE
          EXTOPT                  =  TAUCS(I,K) + EXTOPT
          OMARCS                  =  TAUOMC(I,K) + TAUR(I,K)
          SSALB                   =  OMARCS / EXTOPT
          SF2                     =  F2(I,K) / OMARCS
          SF                      =  SSALB * SF2
          TAU2                    =  EXTOPT * (1.0 - SF)  
          OM2                     = (SSALB - SF) / (1.0 - SF)
          COW                     =  1.0 - OM2
          SSGAS2                  = (TAUOMGC(I,K) / OMARCS - SF2) / 
     1                              (1.0 - SF2)
          COWG                    =  1.0 - OM2 * SSGAS2
          ALAMD                   =  SQRT(3.0 * COW * COWG)
          U                       =  1.50 * COWG / ALAMD
          U2                      =  U + U
          UU                      =  U * U
          SDTR                    =  EXP(- TAU2 / RMU(I))
          EFUN                    =  EXP(- ALAMD * TAU2)
          EFUN2                   =  EFUN * EFUN
          X3                      = (UU - U2 + 1.0) * EFUN2
          RN                      =  1.0 / (UU + U2 + 1.0 - X3)
          X4                      =  ALAMD * RMU(I)
          Y                       =  1.0 - X4 * X4
          YX                      =  SIGN( MAX( ABS(Y), LOWBND), Y )
          DM                      =  OM2 / (1.0 - X4 * X4)
          GSCW                    =  SSGAS2 * COW
          APPGM                   = (C1(I) + 0.50 +
     1                               GSCW * (C1(I) + C2(I))) * DM
          APMGM                   = (C1(I) - 0.50 +
     1                               GSCW * (C1(I) - C2(I))) * DM
          SRDF                    = (UU - 1.0) * (1.0 - EFUN2) * RN
          STDF                    = (U2 + U2) * EFUN * RN 
          SRDR                    =  APPGM * SRDF + APMGM *
     1                              (STDF * SDTR - 1.0) 
          STDR                    =  APPGM * STDF + (APMGM * SRDF - 
     1                               APPGM + 1.0) * SDTR
          IF (NBLK(I,K) .EQ. 3)                                     THEN
            X5                    =  A1(I,9) * CLD(I,K)
            RDF(I,2,K)            =  RDF(I,1,K) + X5 * 
     1                              (SRDF - RDF(I,1,K))
            TDF(I,2,K)            =  TDF(I,1,K) + X5 * 
     1                              (STDF - TDF(I,1,K))
            RDR(I,2,K)            =  RDR(I,1,K) + X5 * 
     1                              (SRDR - RDR(I,1,K))
            TDR(I,2,K)            =  TDR(I,1,K) + X5 * 
     1                              (STDR - TDR(I,1,K))
            DTR(I,2,K)            =  DTR(I,1,K) + X5 * 
     1                              (SDTR - DTR(I,1,K))
            X6                    =  A1(I,10) * CLD(I,K)
            RDF(I,3,K)            =  RDF(I,1,K) + X6 *
     1                              (SRDF - RDF(I,1,K))
            TDF(I,3,K)            =  TDF(I,1,K) + X6 *
     1                              (STDF - TDF(I,1,K))
            RDR(I,3,K)            =  RDR(I,1,K) + X6 *
     1                              (SRDR - RDR(I,1,K))
            TDR(I,3,K)            =  TDR(I,1,K) + X6 *
     1                              (STDR - TDR(I,1,K))
            DTR(I,3,K)            =  DTR(I,1,K) + X6 *
     1                              (SDTR - DTR(I,1,K))
            X7                    =  A1(I,11) * CLD(I,K)
            RDF(I,4,K)            =  RDF(I,1,K) + X7 *
     1                              (SRDF - RDF(I,1,K))
            TDF(I,4,K)            =  TDF(I,1,K) + X7 *
     1                              (STDF - TDF(I,1,K))
            RDR(I,4,K)            =  RDR(I,1,K) + X7 *
     1                              (SRDR - RDR(I,1,K))
            TDR(I,4,K)            =  TDR(I,1,K) + X7 *
     1                              (STDR - TDR(I,1,K))
            DTR(I,4,K)            =  SDTR
          ELSE IF (NBLK(I,K) .EQ. 1)                                THEN
            RDF(I,4,K)            =  RDF(I,1,K)
            TDF(I,4,K)            =  TDF(I,1,K)
            RDR(I,4,K)            =  RDR(I,1,K)
            TDR(I,4,K)            =  TDR(I,1,K)
            DTR(I,4,K)            =  DTR(I,1,K)
            X                     =  CLD(I,K) / CLDM(I,K)
            RDF(I,2,K)            =  RDF(I,1,K) + X * 
     1                              (SRDF - RDF(I,1,K))
            TDF(I,2,K)            =  TDF(I,1,K) + X *
     1                              (STDF - TDF(I,1,K))
            RDR(I,2,K)            =  RDR(I,1,K) + X *
     1                              (SRDR - RDR(I,1,K))
            TDR(I,2,K)            =  TDR(I,1,K) + X *
     1                              (STDR - TDR(I,1,K))
            DTR(I,2,K)            =  SDTR 
            IF (A1(I,2) .GE. CUT)                                   THEN
              Y                   =  X * A1(I,8)
              RDF(I,3,K)          =  RDF(I,1,K) + Y *
     1                              (SRDF - RDF(I,1,K))
              TDF(I,3,K)          =  TDF(I,1,K) + Y *
     1                              (STDF - TDF(I,1,K))
              RDR(I,3,K)          =  RDR(I,1,K) + Y *
     1                              (SRDR - RDR(I,1,K))
              TDR(I,3,K)          =  TDR(I,1,K) + Y *
     1                              (STDR - TDR(I,1,K))
              DTR(I,3,K)          =  DTR(I,1,K) + Y *
     1                              (SDTR - DTR(I,1,K))
            ELSE
              RDF(I,3,K)          =  RDF(I,1,K)
              TDF(I,3,K)          =  TDF(I,1,K)
              RDR(I,3,K)          =  RDR(I,1,K)
              TDR(I,3,K)          =  TDR(I,1,K)
              DTR(I,3,K)          =  DTR(I,1,K)
            ENDIF
          ELSE IF (NBLK(I,K) .EQ. 2)                                THEN
            RDF(I,2,K)            =  RDF(I,1,K)
            TDF(I,2,K)            =  TDF(I,1,K)
            RDR(I,2,K)            =  RDR(I,1,K)
            TDR(I,2,K)            =  TDR(I,1,K)
            DTR(I,2,K)            =  DTR(I,1,K)
            RDF(I,4,K)            =  RDF(I,1,K)
            TDF(I,4,K)            =  TDF(I,1,K)
            RDR(I,4,K)            =  RDR(I,1,K)
            TDR(I,4,K)            =  TDR(I,1,K)
            DTR(I,4,K)            =  DTR(I,1,K)
            X8                    =  CLD(I,K) / CLDM(I,K)
            RDF(I,3,K)            =  RDF(I,1,K) + X8 *
     1                              (SRDF - RDF(I,1,K))
            TDF(I,3,K)            =  TDF(I,1,K) + X8 *
     1                              (STDF - TDF(I,1,K))
            RDR(I,3,K)            =  RDR(I,1,K) + X8 *
     1                              (SRDR - RDR(I,1,K))
            TDR(I,3,K)            =  TDR(I,1,K) + X8 *
     1                              (STDR - TDR(I,1,K))
            DTR(I,3,K)            =  SDTR
          ENDIF
        ENDIF
  200 CONTINUE
C
      DO 300 I = IL1, IL2
C
C----------------------------------------------------------------------C
C     INITIALIZATION FOR THE FIRST LEVEL (LEV1).                       C
C----------------------------------------------------------------------C
C
        ATRAN0                    =  1.0 - TRAN0(I)
        TMDR(I,1,LEV1)            =  TRAN0(I)
        RMDF(I,1,LEV1)            =  ATRAN0
        CUMDTR(I,1,LEV1)          =  TRAN0(I)
        TMDR(I,2,LEV1)            =  TRAN0(I)
        RMDF(I,2,LEV1)            =  ATRAN0
        CUMDTR(I,2,LEV1)          =  TRAN0(I)
        TMDR(I,3,LEV1)            =  TRAN0(I)
        RMDF(I,3,LEV1)            =  ATRAN0
        CUMDTR(I,3,LEV1)          =  TRAN0(I)
        TMDR(I,4,LEV1)            =  TRAN0(I)
        RMDF(I,4,LEV1)            =  ATRAN0
        CUMDTR(I,4,LEV1)          =  TRAN0(I)
C
C----------------------------------------------------------------------C
C     INITIALIZATION FOR THE GROUND LAYER.                             C
C----------------------------------------------------------------------C
C
        RMUR(I,1,LEV)             =  CSALB(I)
        RMUF(I,1,LEV)             =  CSALB(I)
        RMUR(I,2,LEV)             =  ALBSUR(I)
        RMUF(I,2,LEV)             =  ALBSUR(I)
        RMUR(I,3,LEV)             =  ALBSUR(I)
        RMUF(I,3,LEV)             =  ALBSUR(I)
        RMUR(I,4,LEV)             =  ALBSUR(I)
        RMUF(I,4,LEV)             =  ALBSUR(I)
C
        CUMDTR(I,2,LEV)           =  0.0
        CUMDTR(I,3,LEV)           =  0.0
        CUMDTR(I,4,LEV)           =  0.0
  300 CONTINUE
C
C----------------------------------------------------------------------C
C     ADD THE LAYERS DOWNWARD FROM THE SECOND LAYER TO THE SURFACE.    C
C----------------------------------------------------------------------C
C
      DO 450 K = LEV1 + 1, LEV  
        KM1 = K - 1
        L = LEV - K + LEV1
        LP1 = L + 1
        DO 400 I = IL1, IL2
          DMM                     =  TDF(I,1,KM1) / 
     1                              (1.0 - RDF(I,1,KM1) * RMDF(I,1,KM1))
          FMM                     =  RMDF(I,1,KM1) * DMM
          TMDR(I,1,K)             =  CUMDTR(I,1,KM1) * (TDR(I,1,KM1) +
     1                               RDR(I,1,KM1) * FMM) + 
     2                              (TMDR(I,1,KM1) - CUMDTR(I,1,KM1)) * 
     3                               DMM
          RMDF(I,1,K)             =  RDF(I,1,KM1) + TDF(I,1,KM1) * FMM
          CUMDTR(I,1,K)           =  CUMDTR(I,1,KM1) * DTR(I,1,KM1)
C
          IF (A1(I,1) .GE. CUT)                                     THEN
            IF (K .LE. NCT(I))                                      THEN
              TMDR(I,2,K)         =  TMDR(I,1,K)
              RMDF(I,2,K)         =  RMDF(I,1,K)
              CUMDTR(I,2,K)       =  CUMDTR(I,1,K)
            ELSE
              DPP                 =  TDF(I,2,KM1) /
     1                              (1.0 - RMDF(I,2,KM1) * RDF(I,2,KM1))
              FPP                 =  RMDF(I,2,KM1) * DPP
              TMDR(I,2,K)         =  CUMDTR(I,2,KM1) * (TDR(I,2,KM1) +
     1                               RDR(I,2,KM1) * FPP) +
     2                              (TMDR(I,2,KM1) - CUMDTR(I,2,KM1)) *
     3                               DPP
              RMDF(I,2,K)         =  RDF(I,2,KM1) + TDF(I,2,KM1) * FPP
              CUMDTR(I,2,K)       =  CUMDTR(I,2,KM1) * DTR(I,2,KM1)
            ENDIF
          ELSE
            TMDR(I,2,K)           =  1.0
            RMDF(I,2,K)           =  0.0
            CUMDTR(I,2,K)         =  0.0
          ENDIF
C
          IF (A1(I,2) .GE. CUT)                                     THEN
            IF (K .LE. NCT(I))                                      THEN
              TMDR(I,3,K)         =  TMDR(I,1,K)
              RMDF(I,3,K)         =  RMDF(I,1,K)
              CUMDTR(I,3,K)       =  CUMDTR(I,1,K)
            ELSE
              DPP                 =  TDF(I,3,KM1) /
     1                              (1.0 - RMDF(I,3,KM1) * RDF(I,3,KM1))
              FPP                 =  RMDF(I,3,KM1) * DPP
              TMDR(I,3,K)         =  CUMDTR(I,3,KM1) * (TDR(I,3,KM1) +
     1                               RDR(I,3,KM1) * FPP) +
     2                              (TMDR(I,3,KM1) - CUMDTR(I,3,KM1)) *
     3                               DPP
              RMDF(I,3,K)         =  RDF(I,3,KM1) + TDF(I,3,KM1) * FPP
              CUMDTR(I,3,K)       =  CUMDTR(I,3,KM1) * DTR(I,3,KM1)
            ENDIF
C
            IF (A1(I,3) .GE. CUT)                                   THEN
              IF (K .LE. NCT(I))                                    THEN
                TMDR(I,4,K)       =  TMDR(I,1,K)
                RMDF(I,4,K)       =  RMDF(I,1,K)
                CUMDTR(I,4,K)     =  CUMDTR(I,1,K)
              ELSE
                DPP               =  TDF(I,4,KM1) /
     1                              (1.0 - RMDF(I,4,KM1) * RDF(I,4,KM1))
                FPP               =  RMDF(I,4,KM1) * DPP
                TMDR(I,4,K)       =  CUMDTR(I,4,KM1) * (TDR(I,4,KM1) +
     1                               RDR(I,4,KM1) * FPP) +
     2                              (TMDR(I,4,KM1) - CUMDTR(I,4,KM1)) *
     3                               DPP
                RMDF(I,4,K)       =  RDF(I,4,KM1) + TDF(I,4,KM1) * FPP
                CUMDTR(I,4,K)     =  CUMDTR(I,4,KM1) * DTR(I,4,KM1)
              ENDIF
            ELSE
              TMDR(I,4,K)         =  1.0
              RMDF(I,4,K)         =  0.0
              CUMDTR(I,4,K)       =  0.0
            ENDIF
          ELSE
            TMDR(I,3,K)           =  1.0
            RMDF(I,3,K)           =  0.0
            CUMDTR(I,3,K)         =  0.0
            TMDR(I,4,K)           =  1.0
            RMDF(I,4,K)           =  0.0
            CUMDTR(I,4,K)         =  0.0
          ENDIF
C
C----------------------------------------------------------------------C
C     ADD THE LAYERS UPWARD FROM ONE LAYER ABOVE SURFACE TO THE LEV1.  C
C----------------------------------------------------------------------C
C
          UMM                     =  TDF(I,1,L) /
     1                              (1.0 - RDF(I,1,L) * RMUF(I,1,LP1))
          FMM                     =  RMUF(I,1,LP1) * UMM
          RMUR(I,1,L)             =  RDR(I,1,L) + DTR(I,1,L) *
     1                               RMUR(I,1,LP1) * UMM + (TDR(I,1,L) -
     2                               DTR(I,1,L)) * FMM
          RMUF(I,1,L)             =  RDF(I,1,L) + TDF(I,1,L) * FMM
C
          IF (A1(I,1) .GE. CUT)                                     THEN
            UPP                   =  TDF(I,2,L) /
     1                              (1.0 - RMUF(I,2,LP1) * RDF(I,2,L))
            FPP                   =  RMUF(I,2,LP1) * UPP
            RMUR(I,2,L)           =  RDR(I,2,L) + DTR(I,2,L) *
     1                               RMUR(I,2,LP1) * UPP + (TDR(I,2,L) -
     2                               DTR(I,2,L)) * FPP
            RMUF(I,2,L)           =  RDF(I,2,L) + TDF(I,2,L) * FPP
          ELSE
            RMUR(I,2,L)           =  0.0
            RMUF(I,2,L)           =  0.0
          ENDIF
C
          IF (A1(I,2) .GE. CUT)                                     THEN
            UPP                   =  TDF(I,3,L) /
     1                              (1.0 - RMUF(I,3,LP1) * RDF(I,3,L))
            FPP                   =  RMUF(I,3,LP1) * UPP
            RMUR(I,3,L)           =  RDR(I,3,L) + DTR(I,3,L) *
     1                               RMUR(I,3,LP1) * UPP + (TDR(I,3,L) -
     2                               DTR(I,3,L)) * FPP
            RMUF(I,3,L)           =  RDF(I,3,L) + TDF(I,3,L) * FPP
C
            IF (A1(I,3) .GE. CUT)                                   THEN
              UPP                 =  TDF(I,4,L) /
     1                              (1.0 - RMUF(I,4,LP1) * RDF(I,4,L))
              FPP                 =  RMUF(I,4,LP1) * UPP
              RMUR(I,4,L)         =  RDR(I,4,L) + DTR(I,4,L) *
     1                               RMUR(I,4,LP1) * UPP + (TDR(I,4,L) -
     2                               DTR(I,4,L)) * FPP
              RMUF(I,4,L)         =  RDF(I,4,L) + TDF(I,4,L) * FPP
            ELSE
              RMUR(I,4,L)         =  0.0
              RMUF(I,4,L)         =  0.0
            ENDIF
          ELSE
            RMUR(I,3,L)           =  0.0
            RMUF(I,3,L)           =  0.0
            RMUR(I,4,L)           =  0.0
            RMUF(I,4,L)           =  0.0
          ENDIF
  400   CONTINUE
  450 CONTINUE
C
C----------------------------------------------------------------------C
C     ADD DOWNWARD TO CALCULATE THE RESULTANT REFLECTANCE AND          C
C     TRANSMITTANCE AT FLUX LEVELS.                                    C
C----------------------------------------------------------------------C
C
      DO 550 K = LEV1, LEV
        KM1 = K - 1
        DO 500 I = IL1, IL2
          DMM                     =  1.0 /
     1                              (1.0 - RMUF(I,1,K) * RMDF(I,1,K))
          X                       =  CUMDTR(I,1,K) * RMUR(I,1,K)
          Y                       =  TMDR(I,1,K) - CUMDTR(I,1,K)
          TRAN(I,1,K)             =  CUMDTR(I,1,K) +
     1                              (X * RMDF(I,1,K) + Y) * DMM
          REFL(I,1,K)             = (X + Y * RMUF(I,1,K)) * DMM
C
          IF (A1(I,1) .GE. CUT)                                     THEN
            DPP                   =  1.0 /
     1                              (1.0 - RMUF(I,2,K) * RMDF(I,2,K))
            X                     =  CUMDTR(I,2,K) * RMUR(I,2,K)
            Y                     =  TMDR(I,2,K) - CUMDTR(I,2,K)
            TRAN(I,2,K)           =  A1(I,1) * (CUMDTR(I,2,K) +
     1                              (X * RMDF(I,2,K) + Y) * DPP) +
     2                               A1(I,7) * TRAN(I,1,K)
            REFL(I,2,K)           =  A1(I,1) * (X + Y * RMUF(I,2,K)) *
     1                               DPP + A1(I,7) * REFL(I,1,K)
          ELSE
            TRAN(I,2,K)           =  A1(I,7) * TRAN(I,1,K)
            REFL(I,2,K)           =  A1(I,7) * REFL(I,1,K)
          ENDIF
C
          IF (A1(I,2) .GE. CUT)                                     THEN
            DPP                   =  1.0 /
     1                              (1.0 - RMUF(I,3,K) * RMDF(I,3,K))
            X                     =  CUMDTR(I,3,K) * RMUR(I,3,K)
            Y                     =  TMDR(I,3,K) - CUMDTR(I,3,K)
            TRANPP                =  CUMDTR(I,3,K) +
     1                              (X * RMDF(I,3,K) + Y) * DPP
            REFLPP                = (X + Y * RMUF(I,3,K)) * DPP
            TRAN(I,2,K)           =  A1(I,2) * TRANPP + TRAN(I,2,K)
            REFL(I,2,K)           =  A1(I,2) * REFLPP + REFL(I,2,K)
C
            IF (A1(I,3) .GE. CUT)                                   THEN
              DPP                 =  1.0 /
     1                              (1.0 - RMUF(I,4,K) * RMDF(I,4,K))
              X                   =  CUMDTR(I,4,K) * RMUR(I,4,K)
              Y                   =  TMDR(I,4,K) - CUMDTR(I,4,K)
              TRANPP              =  CUMDTR(I,4,K) +
     1                              (X * RMDF(I,4,K) + Y) * DPP
              REFLPP              = (X + Y * RMUF(I,4,K)) * DPP
              TRAN(I,2,K)         =  A1(I,3) * TRANPP + TRAN(I,2,K)
              REFL(I,2,K)         =  A1(I,3) * REFLPP + REFL(I,2,K)
            ENDIF
          ENDIF
c          if(II.eq.16.and.JJ.eq.12) then
c           print*,'SWTRAN',II,JJ,K,TRAN(I,2,K),REFL(I,2,K),A1(I,3)
c          endif
  500   CONTINUE
  550 CONTINUE
C
      RETURN
      END

      SUBROUTINE LWTRAN (FU, FD, SLWF, TAUCI, OMCI, GCI, FL, TAUAL,
     1                   TAUG, BF, BS, URBF, DBF, EM0, CLD, CLDM, ANU,
     2                   NCT, NCD, NCU, NCUM, NCDM, LEV1, MAXC,
     3                   CUT, IL1, IL2, ILG, LAY, LEV,
     +                   II,JJ)
C
C----------------------------------------------------------------------C
C     CALCULATION OF LONGWAVE RADIATIVE TRANSFER USING ABSORPTION      C
C     APPROXIMATION. THE FINITE CLOUD EFFECT IS PROPERLY CONSIDERED    C
C     WITH RANDOM AND FULL OVERLAP ASSUMPTION. CLOUD SUBGRID           C
C     VARIABILITY IS INCLUDED (BASED ON LI, 2002 JAS P3302; LI AND     C
C     BARKER JAS P3321).                                               C
C                                                                      C
C     FU:     UPWARD INFRARED FLUX                                     C
C     FD:     DOWNWARD INFRARED FLUX                                   C
C     SLWF:   INPUT SOLAR FLUX AT MODEL TOP LEVEL FOR EACH BAND        C
C     TAUCI:  CLOUD OPTICAL DEPTH FOR THE INFRARED                     C
C     OMCI:   CLOUD SINGLE SCATTERING ALBEDO TIMES OPTICAL DEPTH       C
C     GCI:    CLOUD ASYMMETRY FACTOR TIMES OMCI                        C
C     FL:     SQUARE OF CLOUD ASYMMETRY FACTOR                         C
C     TAUAL:  AEROSOL OPTICAL DEPTH FOR THE INFRARED                   C
C     TAUG:   GASEOUS OPTICAL DEPTH FOR THE INFRARED                   C
C     BF:     BLACKBODY INTENSITY INTEGRATED OVER EACH BAND AT EACH    C
C             LEVEL IN UNITS W / M^2 / SR. THEREFOR A PI FACTOR NEEDED C
C             FOR FLUX                                                 C
C     BS:     THE BLACKBODY INTENSITY AT THE SURFACE.                  C
C     URBF:   U TIMES THE DIFFERENCE OF LOG(BF) FOR TWO NEIGHBOR       C
C             LEVELS USED FOR EXPONENTIAL SOURCE FUNCTION              C
C     DBF:    DIFFERENCE OF BF FOR TWO NEIGHBOR LEVELS USED FOR        C
C             LINEAR SOURCE FUNCTION                                   C
C     EM0:    SURFACE EMISSION                                         C
C     CLD:    CLOUD FRACTION                                           C
C     CLDM:   MAXIMUM PORTION IN EACH CLOUD BLOCK, IN WHICH THE EXACT  C
C             SOLUTION FOR SUBGRID VARIABILITY IS APPLIED LI AND       C
C             BARKER JAS P3321).                                       C
C     ANU:    NU FACTOR FOR CLOUD SUBGRID VARIABILITY                  C
C     S:      TOTAL SCATTERING                                         C
C     SCATBK: BACKWARD SCATTERING                                      C
C     SCATFW: FORWARD SCATTERING                                       C
C     SCATSM: INTERNAL SCATTERING                                      C
C     TAUM:   TAUM(1) A FACTOR RELATED TO TAU IN ZETA FACTOR FOR       C
C             LINEAR SOURCE TERM; TAUM(2) THE CUMULATED TAUM(1) FOR    C
C             SUBGRID VARIABILITY CALCULATION                          C
C     XU:     THE EMISSION PART IN THE UPWARD FLUX TRANSMISSION        C
C             (LI, 2002 JAS P3302)                                     C
C     XD:     THE EMISSION PART IN THE DOWNWARD FLUX TRANSMISSION      C
C     DTR:    DIRECT TRANSMISSION                                      C
C     FY:     UPWARD FLUX FOR PURE CLEAR PORTION (1) AND PURE CLOUD    C
C             PORTION (2)                                              C
C     FX:     THE SAME AS FY BUT FOR THE DOWNWARD FLUX                 C
C     FW:     A TERM FOR TRANSFER WITHIN CLDM                          C
C     NCT:    THE HIGHEST CLOUD TOP LEVEL FOR THE LONGITUDE AND        C
C             LATITUDE LOOP (ILG)                                      C
C     NCD:    LAYER INSIDE A CLOUD BLOCK ACCOUNTED FROM CLOUD TOP      C
C     NCU:    LAYER INSIDE A CLOUD BLOCK ACCOUNTED FROM CLOUD BOTTOM   C
C     NCUM:   MAXIMUM LOOP NUMBER CLOUD VERTICAL CORRELATION ACCOUNTED C
C             FROM LOWER LEVEL TO HIGHER LEVEL                         C
C     NCDM:   MAXIMUM LOOP NUMBER CLOUD VERTICAL CORRELATION ACCOUNTED C
C             FROM HIGHER LEVEL TO LOWER LEVEL                         C
C----------------------------------------------------------------------C
      IMPLICIT REAL (A-H,O-Z),
     +INTEGER (I-N)
C
      REAL FU(ILG,2,LEV), FD(ILG,2,LEV)
      REAL SLWF(ILG), TAUCI(ILG,LAY), OMCI(ILG,LAY), GCI(ILG,LAY),
     1     FL(ILG,LAY), TAUAL(ILG,LAY), TAUG(ILG,LAY), BF(ILG,LEV), 
     2     BS(ILG), URBF(ILG,LAY), DBF(ILG,LAY), EM0(ILG), CLD(ILG,LAY),
     3     CLDM(ILG,LAY), ANU(ILG,LAY)
      REAL S(ILG,LAY), TERM1(ILG), EMISW(ILG,LAY), SCATBK(ILG,LAY),
     1     SCATFW(ILG,LAY), SCATSM(ILG,4,LAY), TAUM(ILG,4,LAY),
     2     XD(ILG,4,LAY), XU(ILG,4,LAY), DTR(ILG,4,LAY), FX(ILG,4,LEV),
     3     FY(ILG,4,LEV), FW(ILG,4,LEV), S1(ILG), T(ILG)
      INTEGER NCT(ILG), NCD(ILG,LAY), NCU(ILG,LAY), NCUM(LAY), NCDM(LAY)
C
      DATA  RU / 1.6487213 /
C
C----------------------------------------------------------------------C
C     INITIALIZATION FOR FIRST LAYER. CALCULATE THE DOWNWARD FLUX IN   C
C     THE SECOND LAYER                                                 C
C     COMBINE THE OPTICAL PROPERTIES FOR THE INFRARED,                 C
C     1, AEROSOL + GAS; 2, CLOUD + AEROSOL + GAS.                      C
C     FD (FU) IS DOWN (UPWARD) FLUX, FX (FY) IS THE INCIDENT FLUX      C
C     ABOVE (BELOW) THE CONSIDERED LAYER.                              C
C     GAUSSIAN INTEGRATION AND DIFFUSIVITY FACTOR, RU (LI JAS 2000)    C
C     ABOVE MAXC, EXPONENTIAL SOURCE FUNCTION IS USED                  C
C     BELOW MAXC, LINEAR SOURCE FUNCTION IS USED                       C
C----------------------------------------------------------------------C
C
C mxy+
      DO K = 1, LAY 
      DO I = IL1, IL2
       if(TAUAL(I,K).lt.0) then
c        print*,'negative TAUAL',II,JJ,K,TAUAL(I,K)
        TAUAL(I,K)=MAX(TAUAL(I,K),1.e-20)
       endif
       if(TAUG(I,K).lt.0) then
        print*,'negative TAUG',II,JJ,K,TAUG(I,K)
        TAUG(I,K)=MAX(TAUG(I,K),1.e-20)
       endif
      ENDDO
      ENDDO
C mxy+ 
      L1 =  LEV1
      L2 =  LEV1 + 1
      DO 100 I = IL1, IL2
        FD(I,1,LEV1)              =  SLWF(I)
        FD(I,2,LEV1)              =  SLWF(I)
        FX(I,1,LEV1)              =  FD(I,1,LEV1)
        FX(I,2,LEV1)              =  FD(I,2,LEV1)
  100 CONTINUE
C
      DO 150 K = L2, MAXC
        KM1 = K - 1
        DO 125 I = IL1, IL2
          TAUL1                   =  TAUAL(I,KM1) + TAUG(I,KM1)
          RTAUL1                  =  TAUL1 * RU
          DTR(I,1,KM1)            =  EXP (- RTAUL1)
          UBETA                   =  URBF(I,KM1) / (TAUL1 + 1.E-20)
          EPSD                    =  UBETA + 1.0
          EPSU                    =  UBETA - 1.0
C
          IF (ABS(EPSD) .GT. 0.001)                                 THEN
            XD(I,1,KM1)           = (BF(I,K) - BF(I,KM1) * 
     1                               DTR(I,1,KM1)) / EPSD
          ELSE
            XD(I,1,KM1)           =  RTAUL1 * BF(I,KM1) * DTR(I,1,KM1)
          ENDIF
C
          IF (ABS(EPSU) .GT. 0.001)                                 THEN
            XU(I,1,KM1)           = (BF(I,K) * DTR(I,1,KM1) - 
     1                               BF(I,KM1)) / EPSU
          ELSE
            XU(I,1,KM1)           =  RTAUL1 * BF(I,K) * DTR(I,1,KM1)
          ENDIF
C
          FD(I,1,K)               =  FD(I,1,KM1) * DTR(I,1,KM1) + 
     1                               XD(I,1,KM1)
          FD(I,2,K)               =  FD(I,2,KM1) * DTR(I,1,KM1) + 
     1                               XD(I,1,KM1)
          FX(I,1,K)               =  FD(I,2,K)
          FX(I,2,K)               =  FD(I,2,K)
          if(XD(I,1,KM1).gt.2000. .or. XD(I,1,KM1).lt.-1000. .or.
     &       XU(I,1,KM1).gt.2000. .or. XU(I,1,KM1).lt.-1000.)then
          print*,'XU(I,1,KM1)',I,XU(I,1,KM1),XD(I,1,KM1),RTAUL1,
     &     BF(I,KM1),DTR(I,1,KM1),EPSU,TAUL1,TAUAL(I,KM1),TAUG(I,KM1)
          endif
          if(FD(I,2,K).gt.2000. .or. FD(I,2,K).lt.-1000.) then
           print*,'FD(I,2,K)',I,K,FD(I,2,K),FD(I,2,KM1),DTR(I,1,KM1),
     &       XD(I,1,KM1),RTAUL1,BF(I,KM1),EPSU,UBETA,UBETA,
     &       URBF(I,KM1),TAUL1,DTR(I,1,KM1),RTAUL1,RU,TAUAL(I,KM1),
     &       TAUG(I,KM1)
          endif
          if(FD(I,1,K).gt.2000. .or. FD(I,1,K).lt.-1000.) then
           print*,'FD(I,1,K)',I,K,FD(I,1,K),FD(I,1,KM1),DTR(I,1,KM1),
     &       XD(I,1,KM1),BF(I,K),BF(I,KM1),EPSD,UBETA,RTAUL1,TAUL1,
     &       RU,TAUAL(I,KM1),TAUG(I,KM1)
          endif
  125   CONTINUE
  150 CONTINUE 
C
C----------------------------------------------------------------------C
C     ADD THE LAYERS DOWNWARD FROM THE SECOND LAYER TO THE SURFACE.    C
C     DETERMINE THE XU FOR THE UPWARD PATH.                            C
C     USING EXPONENTIAL SOURCE FUNCTION FOR CLR FLUX CALCULATION AND   C
C     ALSO FOR ALL SKY FLUX IN CLOUD FREE LAYERS.                      C
C----------------------------------------------------------------------C
C
      IF (MAXC .LT. LEV)                                            THEN
        DO 250 K = MAXC + 1, LEV
          KM1 = K - 1
          KM2 = KM1 - 1
          DO 225 I = IL1, IL2
            TAUL1                 =  TAUAL(I,KM1) + TAUG(I,KM1)
            RTAUL1                =  TAUL1 * RU
            DTR(I,1,KM1)          =  EXP (- RTAUL1)
            if(RTAUL1.lt.0) then
             print*,'RTAUL1',I,K,RTAUL1,DTR(I,1,KM1),TAUL1,
     &        TAUAL(I,KM1),TAUG(I,KM1)
            endif
C
            UBETA                 =  URBF(I,KM1) / (TAUL1 + 1.E-20)
            EPSD                  =  UBETA + 1.0
            EPSU                  =  UBETA - 1.0
C
            IF (ABS(EPSD) .GT. 0.001)                               THEN
              XD(I,1,KM1)         = (BF(I,K) - BF(I,KM1) * 
     1                               DTR(I,1,KM1)) / EPSD
            ELSE
              XD(I,1,KM1)         =  RTAUL1 * BF(I,KM1) * DTR(I,1,KM1)
            ENDIF
C
            IF (ABS(EPSU) .GT. 0.001)                               THEN
              XU(I,1,KM1)         = (BF(I,K) * DTR(I,1,KM1) - 
     1                               BF(I,KM1)) / EPSU
            ELSE
              XU(I,1,KM1)         =  RTAUL1 * BF(I,K) * DTR(I,1,KM1)
            ENDIF
C
            FD(I,1,K)             =  FD(I,1,KM1) * DTR(I,1,KM1) +
     1                               XD(I,1,KM1)
            IF (CLD(I,KM1) .LT. CUT)                                THEN
              FD(I,2,K)           =  FD(I,2,KM1) * DTR(I,1,KM1) + 
     1                               XD(I,1,KM1)
              FX(I,1,K)           =  FD(I,2,K)
              FX(I,2,K)           =  FD(I,2,K)
            ELSE
              TAUL2               =  TAUCI(I,KM1) + TAUG(I,KM1)
              SSALB               =  OMCI(I,KM1) / (TAUL2 + 1.E-20)
              SF                  =  SSALB * FL(I,KM1)
              W                   = (SSALB - SF) / (1.0 - SF)
              COW                 =  1.0 - W
              TAUM(I,1,KM1)       = (COW * TAUL2 * (1.0 - SF) + 
     1                               TAUAL(I,KM1)) * RU
              ZETA                =  DBF(I,KM1) / TAUM(I,1,KM1)
              TAU2                =  TAUM(I,1,KM1) + TAUM(I,1,KM1)
C
              SANU                =  ANU(I,KM1)
              X                   =  SANU / (SANU + TAUM(I,1,KM1))
              Y                   =  SANU / (SANU + TAU2)
C
              IF (SANU .LE. 0.50)                                   THEN
                DTR(I,2,KM1)      =  SQRT(X)
                DTR2              =  SQRT(Y)
                EMISW(I,KM1)      =  ZETA * (SQRT(1.0 + TAU2) - 1.0)
                EMBK              =  1.0 - SQRT(1.0 + TAU2 + TAU2)
               if(DTR(I,2,KM1).gt.2000. .or. DTR(I,2,KM1).lt.-1000..or.
     &            DTR2.gt.2000. .or. DTR2.lt.-1000. .or.
     &            EMISW(I,KM1).gt.2000. .or. EMISW(I,KM1).lt.-1000..or.
     &            EMBK.gt.2000. .or. EMBK.lt.-1000.) then
                print*,'DTR(I,2,KM1),DTR2',I,K,DTR(I,2,KM1),DTR2,X,Y,
     &           EMISW(I,KM1),EMBK,ZETA,TAU2,SANU,TAU2,TAUM(I,1,KM1),
     &           TAUCI(I,KM1),TAUG(I,KM1),OMCI(I,KM1),FL(I,KM1),
     &           TAUAL(I,KM1),DBF(I,KM1),FD(I,1,KM1),DTR(I,1,KM1),
     &           XD(I,1,KM1),FD(I,2,K)
               endif 
              ELSE IF (SANU .GT. 0.50 .AND. SANU .LE. 1.0)          THEN
                WT                =  2.0 * SANU - 1.0
                SX                =  SQRT(X)
                SY                =  SQRT(Y)
                DTR(I,2,KM1)      =  SX + (X - SX) * WT
                DTR2              =  SY + (Y - SY) * WT
                P1                =  SQRT(1.0 + TAU2) - 1.0
                EMISW(I,KM1)      =  ZETA * (P1 + (LOG(1.0 + 
     1                               TAUM(I,1,KM1)) - P1) * WT)
                P2                =  1.0 - SQRT(1.0 + TAU2 + TAU2)
                EMBK              =  P2 - (LOG(1.0 + TAU2) + P2) * WT
                if(SX.gt.2000. .or. SX.lt.-1000. .or.
     &             SY.gt.2000. .or. SY.lt.-1000. .or.
     &             P1.GT.2000. .or. P1.lt.-1000. .or.
     &             EMISW(I,KM1).gt.2000. .or. EMISW(I,KM1).lt.-1000..or.
     &             P2.gt.2000. .or. P2.lt.-1000. .or.
     &             EMBK.gt.2000. .or. EMBK.lt.-1000.) then
                print*,'SX,SY',I,K,SX,SY,P1,P2,EMISW(I,KM1),EMBK,
     &           SANU,X,Y,DTR(I,2,KM1),DTR2,TAU2,ZETA,TAUM(I,1,KM1)
                endif 
              ELSE IF (SANU .GT. 1.0 .AND. SANU .LE. 2.0)           THEN
                WT                =  SANU - 1.0
                DTR(I,2,KM1)      =  X + (X * X - X) * WT
                DTR2              =  Y + (Y * Y - Y) * WT
                Z                 =  SANU / (SANU - 1.0)
                P1                =  LOG(1.0 + TAUM(I,1,KM1))
                EMISW(I,KM1)      =  ZETA * (P1 + (Z * (1.0 - X) - P1) *
     1                               WT)
                P2                =  - LOG(1.0 + TAU2)
                EMBK              =  P2 + (Z * (Y - 1.0) - P2) * WT
               if(P1.gt.2000. or. P1.lt.-1000. .or.
     &            P2.gt.2000. or. P2.lt.-1000.) then
               print*,'P1,P2',I,K,P1,P2,SANU,DTR(I,2,KM1),TAUM(I,1,KM1),
     &          EMISW(I,KM1),TAU2
               endif
              ELSE IF (SANU .GT. 2.0 .AND. SANU .LE. 3.0)           THEN
                X2                =  X * X
                Y2                =  Y * Y
                WT                =  SANU - 2.0
                DTR(I,2,KM1)      =  X2 + (X * X2 - X2) * WT
                DTR2              =  Y2 + (Y * Y2 - Y2) * WT
                Z                 =  SANU / (SANU - 1.0)
                EMISW(I,KM1)      =  Z * ZETA *
     1                              (1.0 - X + (X - X2) * WT)
                EMBK              =  Z * (Y - 1.0 + (Y2 - Y) * WT)
              ELSE IF (SANU .GT. 3.0 .AND. SANU .LE. 4.0)           THEN
                X2                =  X * X
                Y2                =  Y * Y
                X3                =  X2 * X
                Y3                =  Y2 * Y
                WT                =  SANU - 3.0
                DTR(I,2,KM1)      =  X3 + (X2 * X2 - X3) * WT
                DTR2              =  Y3 + (Y2 * Y2 - Y3) * WT
                Z                 =  SANU / (SANU - 1.0)
                EMISW(I,KM1)      =  Z * ZETA *
     1                              (1.0 - X2 + (X2 - X3) * WT)
                EMBK              =  Z * (Y2 - 1.0 + (Y3 - Y2) * WT)
C
C----------------------------------------------------------------------C
C     FOR ANU > 4, THE INHOMOENEITY EFFECT IS VERY WEAK, FOR SAVING    C
C     THE INTEGER ANU IS ASSUMED. FOR ANU > 20, HOMOGENOUS IS ASSUMED  C 
C----------------------------------------------------------------------C
C
              ELSE IF (SANU .GT. 4.0 .AND. SANU .LE. 20.0)          THEN
                NANU              =  INT(SANU)
                DTR(I,2,KM1)      =  X ** NANU
                DTR2              =  Y ** NANU
                Z                 =  SANU / (SANU - 1.0)
                EMISW(I,KM1)      =  Z * ZETA * (1.0 - DTR(I,2,KM1) / X)
                EMBK              =  Z * (DTR2 / Y - 1.0)
              ELSE
                EMISW(I,KM1)      =  ZETA * (1.0 - EXP(- TAUM(I,1,KM1)))
                EMBK              = (EXP(- TAU2) - 1.0)
               if(EMISW(I,KM1).gt.2000. .or. EMISW(I,KM1).lt.-1000.
     &          .or. EMBK.gt.2000. .or. EMBK.lt.-1000.) then
                print*,'EMISW(I,KM1),EMBK',I,KM1,EMISW(I,KM1),EMBK,
     &           TAUM(I,1,KM1),TAU2
               endif
              ENDIF
C
              XD(I,2,KM1)         =  BF(I,K) - BF(I,KM1) * 
     1                               DTR(I,2,KM1) - EMISW(I,KM1) 
              XU(I,2,KM1)         =  BF(I,KM1) - BF(I,K) * 
     1                               DTR(I,2,KM1) + EMISW(I,KM1)
C
              WGRCOW              =  W * GCI(I,KM1) / COW
              TAUDTR              =  TAUM(I,1,KM1) * DTR(I,2,KM1)
C
              SCATFW(I,KM1)       =  WGRCOW * X * TAUDTR
              SCATBK(I,KM1)       =  0.5 * WGRCOW * (DTR2 - 1.0)
C
              X                   =  WGRCOW * (2.0 * EMISW(I,KM1) +
     1                              (0.5 * EMBK - TAUDTR) * ZETA)
              SCATSM(I,1,KM1)     =  - SCATBK(I,KM1) * BF(I,K) -
     1                               SCATFW(I,KM1) * BF(I,KM1) - X
              SCATSM(I,2,KM1)     =  - SCATBK(I,KM1) * BF(I,KM1) -
     1                               SCATFW(I,KM1) * BF(I,K) + X
C
              IF (K .EQ. L2)                                        THEN
                FX(I,1,K)         =  FX(I,1,KM1) * DTR(I,1,KM1) + 
     1                               XD(I,1,KM1)
                FX(I,2,K)         =  FX(I,2,KM1) * DTR(I,2,KM1) + 
     1                               XD(I,2,KM1)
              ELSE IF (CLD(I,KM1) .LE. CLD(I,KM2))                  THEN
                FX(I,1,K)         = ( FX(I,2,KM1) + (1.0 - CLD(I,KM2)) /
     1                              (1.0 - CLD(I,KM1) - 1.E-10) * 
     2                              (FX(I,1,KM1) - FX(I,2,KM1)) ) * 
     3                               DTR(I,1,KM1) + XD(I,1,KM1)
                FX(I,2,K)         =  FX(I,2,KM1) * DTR(I,2,KM1) + 
     1                               XD(I,2,KM1)
              ELSE IF (CLD(I,KM1) .GT. CLD(I,KM2))                  THEN
                FX(I,1,K)         =  FX(I,1,KM1) * DTR(I,1,KM1) + 
     1                               XD(I,1,KM1)
                FX(I,2,K)         = (FX(I,1,KM1) + 
     1                               CLD(I,KM2) / CLD(I,KM1) * 
     2                              (FX(I,2,KM1) - FX(I,1,KM1))) * 
     3                               DTR(I,2,KM1) + XD(I,2,KM1)
              ENDIF
C
              FD(I,2,K)           =  FX(I,1,K) + CLD(I,KM1) * 
     1                              (FX(I,2,K) - FX(I,1,K))
            ENDIF   
  225     CONTINUE
  250   CONTINUE
      ENDIF
C
C----------------------------------------------------------------------C
C     INITIALIZATION FOR SURFACE                                       C
C----------------------------------------------------------------------C
      K = LEV - 1
      DO 300 I = IL1, IL2
C mxy+
        EMBS                      =  EM0(I) * BS(I)
        ABSE0                     =  1.0 - EM0(I)
        FU(I,1,LEV)               =  EMBS + ABSE0 * FD(I,1,LEV)
        FY(I,1,LEV)               =  EMBS + ABSE0 * FX(I,1,LEV)
        FY(I,2,LEV)               =  EMBS + ABSE0 * FX(I,2,LEV)
        FU(I,2,LEV)               =  FY(I,1,LEV) + CLD(I,K) *
     1                              (FY(I,2,LEV) - FY(I,1,LEV))
        FW(I,2,LEV)               =  FY(I,2,LEV)
C mxy+
C        FU(I,1,LEV)               =  FD(I,1,LEV) + EM0(I) *
c     1                              (BS(I) - FD(I,1,LEV))
c        FY(I,1,LEV)               =  FX(I,1,LEV) + EM0(I) *
c     1                              (BS(I) - FX(I,1,LEV))
c        FY(I,2,LEV)               =  FX(I,2,LEV) + EM0(I) *
c     1                              (BS(I) - FX(I,2,LEV))
c        FU(I,2,LEV)               =  FY(I,1,LEV) + CLD(I,K) *
c     1                              (FY(I,2,LEV) - FY(I,1,LEV))
c        FW(I,2,LEV)               =  FY(I,2,LEV)
C
C----------------------------------------------------------------------C
C     DETERMINING THE UPWARD FLUX FOR THE FIRST LAY ABOVE SURFACE      C
C----------------------------------------------------------------------C
C
        FU(I,1,K)                 =  FU(I,1,LEV) * DTR(I,1,K) +
     1                               XU(I,1,K)
C
        FY(I,1,K)                 =  FY(I,1,LEV) * DTR(I,1,K) +
     1                               XU(I,1,K)
        IF (CLD(I,K) .LT. CUT)                                      THEN
          S(I,K)                  =  0.0
          FY(I,2,K)               =  FY(I,2,LEV) * DTR(I,1,K) +
     1                               XU(I,1,K)
          FU(I,2,K)               =  FY(I,1,K)
          FW(I,2,K)               =  FY(I,2,K)
          TAUM(I,2,K)             =  0.0
        ELSE
          S(I,K)                  =  SCATBK(I,K) * FX(I,2,K) +
     1                               SCATFW(I,K) * FY(I,2,LEV) +
     2                               SCATSM(I,2,K)
          FY(I,2,K)               =  FY(I,2,LEV) * DTR(I,2,K) +
     1                               XU(I,2,K) + S(I,K)
          FU(I,2,K)               =  FY(I,1,K) + CLD(I,K) *
     1                              (FY(I,2,K) - FY(I,1,K))
          FW(I,2,K)               =  FY(I,2,K)
          FW(I,2,K)               =  FY(I,2,K)
          TAUM(I,2,K)             =  TAUM(I,1,K)
        ENDIF
  300 CONTINUE
C
C----------------------------------------------------------------------C
C     ADD THE LAYERS UPWARD FROM THE SECOND LAYER TO MAXC              C
C     SCATTERING EFFECT FOR UPWARD PATH IS INCLUDED                    C
C----------------------------------------------------------------------C
C
      DO 450 K = LEV - 2, MAXC, - 1
        KP1 = K + 1
        KP2 = K + 2
        DO 400 I = IL1, IL2
          IF (K .GE. NCT(I))                                        THEN
            FU(I,1,K)             =  FU(I,1,KP1) * DTR(I,1,K) + 
     1                               XU(I,1,K)
C
            IF (CLD(I,K) .LT. CUT)                                  THEN
              FU(I,2,K)           =  FU(I,2,KP1) * DTR(I,1,K) + 
     1                               XU(I,1,K)
              FY(I,1,K)           =  FU(I,2,K)
              FY(I,2,K)           =  FU(I,2,K)
              FW(I,2,K)           =  FU(I,2,K)
              TAUM(I,2,K)         =  0.0
            ELSE
C
C----------------------------------------------------------------------C
C     FY(I,2,K) CONTAINS UNPERTURBED + BACKWARD SCATTERING EFFECT +    C
C     FORWARD SCATTERING EFFECT + INTERNAL SCATTERING EFFECT           C
C    (LI AND FU, JAS 2000)                                             C
C----------------------------------------------------------------------C
C
              IF (CLD(I,K) .LE. CLD(I,KP1) .OR.
     1                          CLD(I,K) - CLDM(I,K) .LT. CUT)      THEN
C
                FY(I,1,K)         = ( FY(I,2,KP1) + (1.0 - CLD(I,KP1)) /
     1                              (1.0 - CLD(I,K) - 1.E-10) *
     2                              (FY(I,1,KP1) - FY(I,2,KP1)) ) * 
     3                               DTR(I,1,K) + XU(I,1,K)
                T(I)              =  FY(I,2,KP1)
              ELSE
                FY(I,1,K)         =  FY(I,1,KP1) * DTR(I,1,K) + 
     1                               XU(I,1,K)
                T(I)              =  FY(I,1,KP1) + 
     1                              (CLD(I,KP1) - CLDM(I,KP1)) / 
     2                              (CLD(I,K) - CLDM(I,K)) *
     3                              (FY(I,2,KP1) - FY(I,1,KP1))
              ENDIF
C
              BKINS               =  SCATBK(I,K) * FX(I,2,K) + 
     1                               SCATSM(I,2,K)
              FY(I,2,K)           =  T(I) * (DTR(I,2,K) + SCATFW(I,K)) +
     1                               XU(I,2,K) + BKINS 
              TAUM(I,2,K)         =  TAUM(I,2,KP1) + TAUM(I,1,K)
              S1(I)               =  0.0
              S(I,K)              =  BKINS + SCATFW(I,K) * FY(I,2,KP1) 
              TERM1(I)            =  0.0 
C
              IF (NCU(I,K) .GT. 1)                                  THEN
                KX = K + NCU(I,K)
                KXM = KX - 1
C
                SANU              =  ANU(I,KXM)
                ANUTAU            =  SANU / (SANU + TAUM(I,2,K))
                IF (SANU .LE. 0.50)                                 THEN
                  DTRGW           =  SQRT(ANUTAU)
                ELSE IF (SANU .GT. 0.50 .AND. SANU .LE. 1.0)        THEN
                  X               =  SQRT(ANUTAU)
                  DTRGW           =  X + 2.0 * (SANU - 0.50) *
     1                              (ANUTAU - X)
                ELSE IF (SANU .GT. 1.0 .AND. SANU .LE. 2.0)         THEN
                  DTRGW           =  ANUTAU + (SANU - 1.0) * ANUTAU *
     1                              (ANUTAU - 1.0)
                ELSE IF (SANU .GT. 2.0 .AND. SANU .LE. 3.0)         THEN
                  X               =  ANUTAU * ANUTAU
                  DTRGW           =  X + (SANU - 2.0) * X *
     1                              (ANUTAU - 1.0)
                ELSE IF (SANU .GT. 3.0 .AND. SANU .LE. 4.0)         THEN
                  X               =  ANUTAU * ANUTAU * ANUTAU
                  DTRGW           =  X + (SANU - 3.0) * X *
     1                              (ANUTAU - 1.0)
                ELSE IF (SANU .GT. 4.0 .AND. SANU .LE. 20.0)        THEN
                  DTRGW           =  ANUTAU ** (INT(SANU))
                ELSE
                  DTRGW           =  EXP(- TAUM(I,2,K))
                 if(DTRGW.gt.2000. .or. DTRGW.lt.-1000.) then
                 print*,'DTRGW',DTRGW,TAUM(I,2,K)
                 endif
                ENDIF
C
                TERM1(I)          = (FW(I,2,KX) - BF(I,KX)) * DTRGW
                S1(I)             = (EMISW(I,KP1) + S(I,KP1)) * 
     1                               DTR(I,2,K)
              ENDIF
            ENDIF
          ENDIF 
  400   CONTINUE
C
C----------------------------------------------------------------------C
C     DETERMINING THE TERMS GOING INTO THE CORRELATION CALCULATIONS    C
C     FOR SUBGRID VARIABILITY FOR CLDM PORTION.                        C
C----------------------------------------------------------------------C
C
        IF (NCUM(K) .GT. 2)                                         THEN
          DO 420 KK = KP2, K + NCUM(K) - 1
          DO 420 I = IL1, IL2
            IF (K .GE. NCT(I) .AND. CLD(I,K) .GE. CUT .AND.
     1          NCU(I,K) .GT. 2 .AND. KK .LE. K + NCU(I,K) - 1)     THEN
C
              SANU                =  ANU(I,KK)
              ANUTAU              =  SANU / (SANU + 
     1                               TAUM(I,2,K) - TAUM(I,2,KK))
              IF (SANU .LE. 0.50)                                   THEN
                DTRGW             =  SQRT(ANUTAU)
              ELSE IF (SANU .GT. 0.50 .AND. SANU .LE. 1.0)          THEN
                X                 =  SQRT(ANUTAU)
                DTRGW             =  X + 2.0 * (SANU - 0.50) *
     1                              (ANUTAU - X)
              ELSE IF (SANU .GT. 1.0 .AND. SANU .LE. 2.0)           THEN
                DTRGW             =  ANUTAU + (SANU - 1.0) * ANUTAU *
     1                              (ANUTAU - 1.0)
              ELSE IF (SANU .GT. 2.0 .AND. SANU .LE. 3.0)           THEN
                X                 =  ANUTAU * ANUTAU
                DTRGW             =  X + (SANU - 2.0) * X *
     1                              (ANUTAU - 1.0)
              ELSE IF (SANU .GT. 3.0 .AND. SANU .LE. 4.0)           THEN
                X                 =  ANUTAU * ANUTAU * ANUTAU
                DTRGW             =  X + (SANU - 3.0) * X *
     1                              (ANUTAU - 1.0)
              ELSE IF (SANU .GT. 4.0 .AND. SANU .LE. 20.0)          THEN
                DTRGW             =  ANUTAU ** (INT(SANU))
              ELSE
                DTRGW             =  EXP(- TAUM(I,2,KK) + TAUM(I,2,KK))
               if(DTRGW.gt.2000. .or. DTRGW.lt.-1000.) then
               print*,'DTRGW',DTRGW,TAUM(I,2,KK),TAUM(I,2,KK)
               endif
              ENDIF
C
              S1(I)               =  S1(I) + 
     1                              (EMISW(I,KK) + S(I,KK)) * DTRGW
            ENDIF
  420     CONTINUE
        ENDIF
C
C----------------------------------------------------------------------C
C     IN CLDM REGION CONSIDER THE CORRELATION BETWEEN DIFFERENT LAYERS C
C----------------------------------------------------------------------C
C
        DO 430 I = IL1, IL2
          IF (K .GE. NCT(I))                                        THEN
            IF (CLD(I,K) .GE. CUT)                                  THEN
              IF (NCU(I,K) .EQ. 1)                                  THEN
                FW(I,2,K)         =  FY(I,2,K) 
                FU(I,2,K)         =  FY(I,1,K) + CLD(I,K) * (FY(I,2,K) -
     1                               FY(I,1,K)) 
              ELSE
                FW(I,2,K)         =  TERM1(I) + S1(I) + BF(I,K) + 
     1                               EMISW(I,K) + S(I,K)
                FU(I,2,K)         =  CLDM(I,K) * (FW(I,2,K) - 
     1                               FY(I,2,K)) + FY(I,1,K) + 
     2                               CLD(I,K) * (FY(I,2,K) - FY(I,1,K))
              ENDIF
            ENDIF
          ENDIF
  430   CONTINUE        
  450 CONTINUE
C
C----------------------------------------------------------------------C
C     ADD THE LAYERS UPWARD ABOVE THE HIGHEST CLOUD  TO THE TOA, NO    C
C     SCATTERING                                                       C
C----------------------------------------------------------------------C
C
      DO 550 K = LEV - 1, L1, - 1
        KP1 = K + 1
C
        DO 500 I = IL1, IL2
          IF (KP1 .LE. NCT(I))                                      THEN
            FU(I,1,K)             =  FU(I,1,KP1) * DTR(I,1,K) + 
     1                               XU(I,1,K)
            FU(I,2,K)             =  FU(I,2,KP1) * DTR(I,1,K) + 
     1                               XU(I,1,K)
          if(FU(I,2,K).gt.2000. .or. FU(I,2,K).lt.-1000.) then
           print*,'FU(I,2,K)',FU(I,2,KP1),DTR(I,1,K),XU(I,1,K),
     &        FU(I,1,K),FU(I,1,KP1)
          endif
          ENDIF
C
C----------------------------------------------------------------------C
C     SCATTERING EFFECT FOR DOWNWARD PATH AT THE TOP LAYER OF THE      C
C     HIGHEST CLOUD                                                    C
C----------------------------------------------------------------------C
C
          IF (K .EQ. NCT(I))                                        THEN
            FW(I,1,K)             =  FX(I,1,K)
            FWINS                 =  SCATSM(I,1,K) + 
     1                               SCATFW(I,K) * FX(I,2,K)
            FMBK                  =  FX(I,2,K) * DTR(I,2,K) + 
     1                               XD(I,2,K) + FWINS
            FX(I,2,KP1)           =  FMBK + SCATBK(I,K) * FY(I,2,KP1)
            TAUM(I,2,K)           =  TAUM(I,1,K)
            S(I,K)                =  SCATBK(I,K) * FW(I,2,KP1) + FWINS
C
            FW(I,1,KP1)           =  FMBK + SCATBK(I,K) * FW(I,2,KP1)
            FD(I,2,KP1)           =  FX(I,1,KP1) + CLD(I,K) *
     1                              (FX(I,2,KP1) - FX(I,1,KP1))
          ENDIF
  500   CONTINUE
  550 CONTINUE
C
C----------------------------------------------------------------------C
C     SCATTERING EFFECT FOR DOWNWARD PATH IN FROM MAXC TO THE SURFACE  C
C----------------------------------------------------------------------C
C
      DO 750 K = MAXC + 2, LEV
        KM1 = K - 1
        KM2 = K - 2
        KM3 = K - 3  
        DO 700 I = IL1, IL2
          IF (KM2 .GE. NCT(I))                                      THEN
            IF (CLD(I,KM1) .LT. CUT)                                THEN
              FD(I,2,K)           =  FD(I,2,KM1) * DTR(I,1,KM1) + 
     1                               XD(I,1,KM1)
              FX(I,1,K)           =  FD(I,2,K)
              FX(I,2,K)           =  FD(I,2,K)
              FW(I,1,K)           =  FD(I,2,K)
              TAUM(I,2,KM1)       =  0.0
            ELSE 
              IF (CLD(I,KM1) .LE. CLD(I,KM2) .OR. 
     1                        CLD(I,KM1) - CLDM(I,KM1) .LT. CUT)    THEN
C
                FX(I,1,K)         = (FX(I,2,KM1) + (1.0 - CLD(I,KM2)) /
     1                              (1.0 - CLD(I,KM1) - 1.E-10) * 
     2                              (FX(I,1,KM1) - FX(I,2,KM1))) * 
     3                               DTR(I,1,KM1) + XD(I,1,KM1)
                T(I)              =  FX(I,2,KM1)
              ELSE
                FX(I,1,K)         =  FX(I,1,KM1) * DTR(I,1,KM1) + 
     1                               XD(I,1,KM1)
                T(I)              =  FX(I,1,KM1) + 
     1                              (CLD(I,KM2) - CLDM(I,KM2)) / 
     2                              (CLD(I,KM1) - CLDM(I,KM1)) * 
     3                              (FX(I,2,KM1) -  FX(I,1,KM1))
              ENDIF
C
              FX(I,2,K)           =  T(I) * DTR(I,2,KM1) + XD(I,2,KM1) +
     1                               SCATBK(I,KM1) * FY(I,2,K) + 
     2                               SCATFW(I,KM1) * T(I) + 
     3                               SCATSM(I,1,KM1)
C
              TAUM(I,2,KM1)       =  TAUM(I,2,KM2) + TAUM(I,1,KM1)
              S1(I)               =  0.0
              S(I,KM1)            =  SCATBK(I,KM1) * FW(I,2,K) + 
     1                               SCATFW(I,KM1) * FW(I,1,KM1) + 
     2                               SCATSM(I,1,KM1)
              TERM1(I)            =  0.0 
C
              IF (NCD(I,KM1) .GT. 1)                                THEN
                KX = K - NCD(I,KM1)
                SANU              =  ANU(I,KX)
                ANUTAU            =  SANU / (SANU + TAUM(I,2,KM1))
                IF (SANU .LE. 0.50)                                 THEN
                  DTRGW           =  SQRT(ANUTAU)
                ELSE IF (SANU .GT. 0.50 .AND. SANU .LE. 1.0)        THEN
                  X               =  SQRT(ANUTAU)
                  DTRGW           =  X + 2.0 * (SANU - 0.50) *
     1                              (ANUTAU - X)
                ELSE IF (SANU .GT. 1.0 .AND. SANU .LE. 2.0)         THEN
                  DTRGW           =  ANUTAU + (SANU - 1.0) * ANUTAU *
     1                              (ANUTAU - 1.0)
                ELSE IF (SANU .GT. 2.0 .AND. SANU .LE. 3.0)         THEN
                  X               =  ANUTAU * ANUTAU
                  DTRGW           =  X + (SANU - 2.0) * X *
     1                              (ANUTAU - 1.0)
                ELSE IF (SANU .GT. 3.0 .AND. SANU .LE. 4.0)         THEN
                  X               =  ANUTAU * ANUTAU * ANUTAU
                  DTRGW           =  X + (SANU - 3.0) * X *
     1                              (ANUTAU - 1.0)
                ELSE IF (SANU .GT. 4.0 .AND. SANU .LE. 20.0)        THEN
                  DTRGW           =  ANUTAU ** (INT(SANU))
                ELSE
                  DTRGW           =  EXP(- TAUM(I,2,KM1))
                 if(DTRGW.gt.2000. .or. DTRGW.lt.-1000.)then
                  print*,'DTRGW,TAUM(I,2,KM1)',DTRGW,TAUM(I,2,KM1)
                 endif
                ENDIF
C
                TERM1(I)          = (FW(I,1,KX) - BF(I,KX)) * DTRGW
                S1(I)             = (S(I,KM2) - EMISW(I,KM2)) * 
     1                               DTR(I,2,KM1)
              ENDIF
            ENDIF
          ENDIF
  700   CONTINUE          
C
C----------------------------------------------------------------------C
C     DETERMINING THE TERMS GOING INTO THE CORRELATION CALCULATIONS    C
C     FOR CLDM PORTION.                                                C
C----------------------------------------------------------------------C
C
        IF (NCDM(KM1) .GT. 2)                                       THEN
C
C----------------------------------------------------------------------C
C     NOTE THAT IN THE FOLLOWING LOOP, "KM1" IS ACTUALLY THE           C
C     REPRESENTATIVE VARIABLE, SO THAT K-NCD(I,KM1) IS ACTUALLY        C
C     KM1-NCD(I,KM1)+1. THE SIMPLER FORM IS USED ONLY FOR              C
C     COMPUTATIONAL EFFICIENCY.                                        C
C----------------------------------------------------------------------C
C
          DO 720 KK = KM3, K - NCDM(KM1), - 1
          DO 720 I = IL1, IL2
            IF (KM2 .GE. NCT(I) .AND. CLD(I,KM1) .GE. CUT .AND.
     1          NCD(I,KM1) .GT. 2 .AND. KK .GE. K - NCD(I,KM1))     THEN
C              
              SANU                =  ANU(I,KK)
              ANUTAU              =  SANU / (SANU + 
     1                               TAUM(I,2,KM1) - TAUM(I,2,KK))
              IF (SANU .LE. 0.50)                                   THEN
                DTRGW             =  SQRT(ANUTAU)
              ELSE IF (SANU .GT. 0.50 .AND. SANU .LE. 1.0)          THEN
                X                 =  SQRT(ANUTAU)
                DTRGW             =  X + 2.0 * (SANU - 0.50) *
     1                              (ANUTAU - X)
              ELSE IF (SANU .GT. 1.0 .AND. SANU .LE. 2.0)           THEN
                DTRGW             =  ANUTAU + (SANU - 1.0) * ANUTAU *
     1                              (ANUTAU - 1.0)
              ELSE IF (SANU .GT. 2.0 .AND. SANU .LE. 3.0)           THEN
                X                 =  ANUTAU * ANUTAU
                DTRGW             =  X + (SANU - 2.0) * X *
     1                              (ANUTAU - 1.0)
              ELSE IF (SANU .GT. 3.0 .AND. SANU .LE. 4.0)           THEN
                X                 =  ANUTAU * ANUTAU * ANUTAU
                DTRGW             =  X + (SANU - 3.0) * X *
     1                              (ANUTAU - 1.0)
              ELSE IF (SANU .GT. 4.0 .AND. SANU .LE. 20.0)          THEN
                DTRGW             =  ANUTAU ** (INT(SANU))
              ELSE
                DTRGW             =  EXP(- TAUM(I,2,KM1) + TAUM(I,2,KK))
               if(DTRGW.gt.2000. .or. DTRGW.lt.-1000.) then
                print*,'DTRGW,TAUM(I,2,KM1),TAUM(I,2,KK)',
     &          DTRGW,TAUM(I,2,KM1),TAUM(I,2,KK)
               endif
              ENDIF
C
              S1(I)               =  S1(I) - 
     1                              (EMISW(I,KK) - S(I,KK)) * DTRGW
            ENDIF
  720     CONTINUE
        ENDIF
C
        DO 730 I = IL1, IL2
          IF (KM2 .GE. NCT(I))                                      THEN
            IF (CLD(I,KM1) .GE. CUT)                                THEN
              IF (NCD(I,KM1) .EQ. 1)                                THEN
                FW(I,1,K)         =  FX(I,2,K)
                FD(I,2,K)         =  FX(I,1,K) + CLD(I,KM1) * 
     1                              (FX(I,2,K) - FX(I,1,K)) 
              ELSE
                FW(I,1,K)         =  TERM1(I) + S1(I) + BF(I,K) - 
     1                               EMISW(I,KM1) + S(I,KM1)
                FD(I,2,K)         =  CLDM(I,KM1) * 
     1                              (FW(I,1,K) - FX(I,2,K)) + 
     2                               FX(I,1,K) + CLD(I,KM1) * 
     3                              (FX(I,2,K) - FX(I,1,K))
              ENDIF
            ENDIF
          ENDIF
  730   CONTINUE
  750 CONTINUE
C
      RETURN
      END

      SUBROUTINE LWTRAGH (FU, FD, SLWF, TAUCI, OMCI, TAUAL, TAUG, BF,
     1                    URBF, CLD, EM0, BS, CUT,
     2                    IL1, IL2, ILG, LAY, LEV)
C
C----------------------------------------------------------------------C
C     IN THE G SPACE WITH INTERVAL CLOSE 1 (VERY LARGE OPTICAL DEPTH)  C
C     OR IN THE CASE WITH CLOUD ABSORPTION IS VERY SMALL OR THE WEIGHT C
C     OF FLUX AND COOLING RATE ARE VERY SMALL. THE CLOUD RADIATIVE     C
C     PROCESS CAN BE HIGHLY SIMPLIFIED. THE ABSORPTION APPROXIMATION   C
C     METHOD IS USED AND CLOUD RANDOM AND MAXIMUM OVERLAP IS           C
C     CONSIDERED, BUT CLOUD SCATTERING AND INHOMOGENEITY ARE IGNORED.  C
C     THE EXPONENTIAL SOURCE PLANCK FUNCTION IS USED WHICH IS MORE     C
C     ACCURATE IN THE REGION ABOVE 200 MB IN COMPARISON WITH LINEAR    C
C     SOURCE FUNCTION                                                  C
C                                                                      C
C     FU:    UPWARD INFRARED FLUX                                      C
C     FD:    DOWNWARD INFRARED FLUX                                    C
C     SLWF:  INPUT SOLAR FLUX AT MODEL TOP LEVEL FOR EACH BAND         C
C     TAUCI: CLOUD OPTICAL DEPTH FOR THE INFRARED                      C
C     OMCI:  CLOUD SINGLE SCATTERING ALBEDO TIMES OPTICAL DEPTH        C
C     TAUAL: AEROSOL OPTICAL DEPTH FOR THE INFRARED                    C
C     TAUG:  GASEOUS OPTICAL DEPTH FOR THE INFRARED                    C
C     BF:    BLACKBODY INTENSITY INTEGRATED OVER EACH BAND AT EACH     C
C            LEVEL IN UNITS W / M^2 / SR.                              C
C     BS:    THE BLACKBODY INTENSITY AT THE SURFACE.                   C
C     URBF:  U TIMES THE DIFFERENCE OF LOG(BF) FOR TWO NEIGHBOR LEVELS C
C            USED FOR EXPONENTIAL SOURCE FUNCTION (LI, 2002 JAS P3302) C
C     CLD:   CLOUD FRACTION                                            C
C     EM0:   SURFACE EMISSION                                          C
C     XU:    THE EMISSION PART IN THE UPWARD FLUX TRANSMISSION         C
C            (LI, 2002 JAS P3302)                                      C
C     XD:    THE EMISSION PART IN THE DOWNWARD FLUX TRANSMISSION       C
C     DTR:   DIRECT TRANSMISSION                                       C
C     FY:    UPWARD FLUX FOR PURE CLEAR PORTION (1) AND PURE CLOUD     C
C            PORTION (2)                                               C
C     FX:    THE SAME AS FY BUT FOR THE DOWNWARD FLUX                  C
C----------------------------------------------------------------------C
C
      IMPLICIT REAL (A-H,O-Z),
     +INTEGER (I-N)

      REAL FU(ILG,2,LEV), FD(ILG,2,LEV)
      REAL SLWF(ILG), TAUCI(ILG,LAY), OMCI(ILG,LAY),
     2     TAUAL(ILG,LAY), TAUG(ILG,LAY), BF(ILG,LEV), URBF(ILG,LAY),
     3     CLD(ILG,LAY), EM0(ILG), BS(ILG)
      REAL XU(ILG,4,LAY), XD(ILG,4,LAY), DTR(ILG,4,LAY), 
     1     FY(ILG,4,LEV), FX(ILG,4,LEV)
C
      DATA  RU / 1.6487213 /
C
C----------------------------------------------------------------------C
C     INITIALIZATION FOR FIRST LAYER. CALCULATE THE DOWNWARD FLUX IN   C
C     THE SECOND LAYER                                                 C
C     COMBINE THE OPTICAL PROPERTIES FOR THE INFRARED,                 C
C     1, AEROSOL + GAS; 2, CLOUD + AEROSOL + GAS.                      C
C     FD (FU) IS DOWN (UPWARD) FLUX                                    C
C     THE OVERLAP BETWEEN SOLAR AND INFRARED IN 4 - 10 UM IS           C
C     CONSIDERED, SLWF IS THE INCOMING SOLAR FLUX                      C
C     SINGULARITY FOR XD AND XU HAS BEEN CONSIDERED AS LI JAS 2002     C
C----------------------------------------------------------------------C
C
      DO 100 I = IL1, IL2
        FD(I,1,1)         =  SLWF(I)
        FD(I,2,1)         =  SLWF(I)
        FX(I,1,1)         =  SLWF(I)
        FX(I,2,1)         =  SLWF(I)
C
        TAUL1             =  TAUAL(I,1) + TAUG(I,1)
        RTAUL1            =  TAUL1 * RU
        DTR(I,1,1)        =  EXP (- RTAUL1)
         if(DTR(I,1,1).gt.2000. .or. DTR(I,1,1).lt.-1000.) then
          print*,'DTR(I,1,1)',I,DTR(I,1,1),RTAUL1,TAUAL(I,1),TAUG(I,1)
         endif
        UBETA             =  URBF(I,1) / (TAUL1 + 1.E-20)
        EPSD              =  UBETA + 1.0
        EPSU              =  UBETA - 1.0
C
        IF (ABS(EPSD) .GT. 0.001)                                   THEN
          XD(I,1,1)       = (BF(I,2) - BF(I,1) * DTR(I,1,1)) / EPSD
        ELSE
          XD(I,1,1)       =  RTAUL1 * BF(I,1) * DTR(I,1,1) 
        ENDIF
        IF (ABS(EPSU) .GT. 0.001)                                   THEN
          XU(I,1,1)       = (BF(I,2) * DTR(I,1,1) - BF(I,1)) / EPSU
        ELSE
          XU(I,1,1)       =  RTAUL1 * BF(I,2) * DTR(I,1,1)
        ENDIF
C
        FD(I,1,2)         =  FD(I,1,1) * DTR(I,1,1) + XD(I,1,1)
C
        IF (CLD(I,1) .LT. CUT)                                      THEN
          FX(I,1,2)       =  FD(I,1,2)
          FX(I,2,2)       =  FD(I,1,2)
          FD(I,2,2)       =  FD(I,1,2)
        ELSE
          TAUL2           =  TAUCI(I,1) + TAUL1
          COW             =  1.0 - OMCI(I,1) / TAUL2
          CTAUL2          =  COW * TAUL2
          CRTAUL2         =  CTAUL2 * RU
          DTR(I,2,1)      =  EXP (- CRTAUL2)
           if(DTR(I,2,1).gt.2000. .or. DTR(I,2,1).lt.-1000.) then
            print*,'DTR(I,2,1)',I,DTR(I,2,1),CRTAUL2,TAUL2,
     &       TAUCI(I,1),TAUL1
           endif
          UBETA           =  URBF(I,1) / (CTAUL2)
          EPSD            =  UBETA + 1.0
          EPSU            =  UBETA - 1.0
C
          IF (ABS(EPSD) .GT. 0.001)                                 THEN
            XD(I,2,1)     = (BF(I,2) - BF(I,1) * DTR(I,2,1)) / EPSD
          ELSE
            XD(I,2,1)     =  CRTAUL2 * BF(I,1) * DTR(I,2,1) 
          ENDIF
          IF (ABS(EPSU) .GT. 0.001)                                 THEN
            XU(I,2,1)     = (BF(I,2) * DTR(I,2,1) - BF(I,1)) / EPSU
          ELSE
            XU(I,2,1)     =  CRTAUL2 * BF(I,2) * DTR(I,2,1)
          ENDIF
C
          FX(I,1,2)       =  FX(I,1,1) * DTR(I,1,1) + XD(I,1,1)
          FX(I,2,2)       =  FX(I,2,1) * DTR(I,2,1) + XD(I,2,1)
          FD(I,2,2)       =  FX(I,1,2) +
     1                       CLD(I,1) * (FX(I,2,2) - FX(I,1,2))
        ENDIF
  100 CONTINUE
C
      DO 250 K = 3, LEV
        KM1 = K - 1
        KM2 = KM1 - 1
        DO 200 I = IL1, IL2
          TAUL1           =  TAUAL(I,KM1) + TAUG(I,KM1)
          RTAUL1          =  TAUL1 * RU
          DTR(I,1,KM1)    =  EXP (- RTAUL1)
          if(DTR(I,1,KM1).gt.2000. .or. DTR(I,1,KM1).lt.-1000.) then
           print*,'DTR(I,1,KM1)',I,DTR(I,1,KM1),RTAUL1,
     &         TAUAL(I,KM1),TAUG(I,KM1)
          endif
          UBETA           =  URBF(I,KM1) / (TAUL1 + 1.E-20)
          EPSD            =  UBETA + 1.0
          EPSU            =  UBETA - 1.0
C
          IF (ABS(EPSD) .GT. 0.001)                                 THEN
            XD(I,1,KM1)   = (BF(I,K) - BF(I,KM1) * DTR(I,1,KM1)) / EPSD
          ELSE
            XD(I,1,KM1)   =  RTAUL1 * BF(I,KM1) * DTR(I,1,KM1)  
          ENDIF
          IF (ABS(EPSU) .GT. 0.001)                                 THEN
            XU(I,1,KM1)   = (BF(I,K) * DTR(I,1,KM1) - BF(I,KM1)) / EPSU
          ELSE
            XU(I,1,KM1)   =  RTAUL1 * BF(I,K) * DTR(I,1,KM1)
          ENDIF
C
          FD(I,1,K)       =  FD(I,1,KM1) * DTR(I,1,KM1) + XD(I,1,KM1)
C
          IF (CLD(I,KM1) .LT. CUT)                                  THEN
            FD(I,2,K)     =  FD(I,2,KM1) * DTR(I,1,KM1) + XD(I,1,KM1)
            FX(I,1,K)     =  FD(I,2,K)
            FX(I,2,K)     =  FD(I,2,K)
          ELSE
            TAUL2         =  TAUCI(I,KM1) + TAUL1
            COW           =  1.0 - OMCI(I,KM1) / TAUL2
            CTAUL2        =  COW * TAUL2
            CRTAUL2       =  CTAUL2 * RU
            DTR(I,2,KM1)  =  EXP (- CRTAUL2)
             if(DTR(I,2,KM1).gt.2000. .or. DTR(I,2,KM1).lt.-1000.) then
             print*,'DTR(I,2,KM1)',I,DTR(I,2,KM1),CRTAUL2,
     &        TAUCI(I,KM1),TAUL1,TAUL2,OMCI(I,KM1),COW
             endif
            UBETA         =  URBF(I,KM1) / (CTAUL2)
            EPSD          =  UBETA + 1.0
            EPSU          =  UBETA - 1.0
C
            IF (ABS(EPSD) .GT. 0.001)                               THEN
              XD(I,2,KM1) = (BF(I,K) - BF(I,KM1) * DTR(I,2,KM1)) / EPSD
            ELSE
              XD(I,2,KM1) =  CRTAUL2 * BF(I,KM1) * DTR(I,2,KM1) 
            ENDIF
            IF (ABS(EPSU) .GT. 0.001)                               THEN
              XU(I,2,KM1) = (BF(I,K) * DTR(I,2,KM1) - BF(I,KM1)) / EPSU
            ELSE
              XU(I,2,KM1) =  CRTAUL2 * BF(I,K) * DTR(I,2,KM1) 
            ENDIF
C
            IF (CLD(I,KM1) .LE. CLD(I,KM2))                         THEN
              FX(I,1,K)   = ( FX(I,2,KM1) + (1.0 - CLD(I,KM2)) /
     1                      (1.0 - CLD(I,KM1) - 1.E-10) *
     2                      (FX(I,1,KM1) - FX(I,2,KM1)) ) *
     3                       DTR(I,1,KM1) + XD(I,1,KM1)
              FX(I,2,K)   =  FX(I,2,KM1) * DTR(I,2,KM1) + XD(I,2,KM1)
            ELSE IF (CLD(I,KM1) .GT. CLD(I,KM2))                    THEN
              FX(I,1,K)   =  FX(I,1,KM1) * DTR(I,1,KM1) + XD(I,1,KM1)
              FX(I,2,K)   = (FX(I,1,KM1) + CLD(I,KM2) / CLD(I,KM1) *
     1                      (FX(I,2,KM1) - FX(I,1,KM1))) *
     2                       DTR(I,2,KM1) + XD(I,2,KM1)
            ENDIF
C
            FD(I,2,K)     =  FX(I,1,K) + CLD(I,KM1) * (FX(I,2,K) -
     1                       FX(I,1,K))
          ENDIF
  200   CONTINUE
  250 CONTINUE
C
      DO 300 I = IL1, IL2
        FU(I,1,LEV)      =  FD(I,1,LEV) + EM0(I) * (BS(I) - FD(I,1,LEV))
        FY(I,1,LEV)      =  FX(I,1,LEV) + EM0(I) * (BS(I) - FX(I,1,LEV))
        FY(I,2,LEV)      =  FX(I,2,LEV) + EM0(I) * (BS(I) - FX(I,2,LEV))
C
        IF (CLD(I,LAY) .GT. CUT)                                    THEN
          FU(I,2,LEV)    =  FY(I,1,LEV) +
     1                      CLD(I,LAY) * (FY(I,2,LEV) - FY(I,1,LEV))
        ELSE
          FU(I,2,LEV)    =  FY(I,2,LEV)
        ENDIF
C
        FU(I,1,LAY)      =  FU(I,1,LEV) * DTR(I,1,LAY) + XU(I,1,LAY)
C
        IF (CLD(I,LAY) .LT. CUT)                                    THEN
          FU(I,2,LAY)    =  FU(I,2,LEV) * DTR(I,1,LAY) + XU(I,1,LAY)
          FY(I,1,LAY)    =  FU(I,2,LAY)
          FY(I,2,LAY)    =  FU(I,2,LAY)
        ELSE
          FY(I,1,LAY)    =  FY(I,1,LEV) * DTR(I,1,LAY) + XU(I,1,LAY)
          FY(I,2,LAY)    =  FY(I,2,LEV) * DTR(I,2,LAY) + XU(I,2,LAY)
          FU(I,2,LAY)    =  FY(I,1,LAY) +
     1                      CLD(I,LAY) * (FY(I,2,LEV) - FY(I,1,LEV))
        ENDIF
  300 CONTINUE
C
      DO 450 K = LEV - 2, 1, - 1
        KP1 = K + 1
        DO 400 I = IL1, IL2
          FU(I,1,K)      =  FU(I,1,KP1) * DTR(I,1,K) + XU(I,1,K)
C
          IF (CLD(I,K) .LT. CUT)                                    THEN
            FU(I,2,K)    =  FU(I,2,KP1) * DTR(I,1,K) + XU(I,1,K)
            FY(I,1,K)    =  FU(I,2,K)
            FY(I,2,K)    =  FU(I,2,K)
          ELSE
            IF (CLD(I,K) .LT. CLD(I,KP1))                           THEN
              FY(I,1,K)  = ( FY(I,2,KP1) + (1.0 - CLD(I,KP1)) /
     1                     (1.0 - CLD(I,K)) * (FY(I,1,KP1) -
     2                      FY(I,2,KP1)) ) * DTR(I,1,K) + XU(I,1,K)
              FY(I,2,K)  =  FY(I,2,KP1) * DTR(I,2,K) + XU(I,2,K)
            ELSE
              FY(I,1,K)  =  FY(I,1,KP1) * DTR(I,1,K) + XU(I,1,K)
              FY(I,2,K)  = ( FY(I,1,KP1) + CLD(I,KP1) / CLD(I,K) *
     1                     (FY(I,2,KP1) - FY(I,1,KP1)) ) * DTR(I,2,K) +
     2                      XU(I,2,K)
            ENDIF
C
            FU(I,2,K)    =  FY(I,1,K) +
     1                      CLD(I,K) * (FY(I,2,K) - FY(I,1,K))
          ENDIF
  400   CONTINUE
  450 CONTINUE
C
      RETURN
      END


      SUBROUTINE TLINE1 (TAUG, COEF1, S, NG, DP, DIP, DT, INPT,  
     1                   LEV1, GH, N, IPLUS, IL1, IL2, ILG, LAY)
C
C----------------------------------------------------------------------C
C     CALCULATION OF OPTICAL DEPTH FOR ONE GAS (LINE CONTRIBUTION ONLY)C
C     THE GASEOUS ABSORPTION COEFFICIENTS IN UNITS OF CM^2 / GRAM.     C
C     S IN MASS MIXING RATIO. ABSORPTION COEFFICIENT ARE CALCULATED    C
C     AT THE TEMPERATURE T FOR THE 18 OR 26 (28) PRESSURE LEVELS.      C
C     ISL = 1 FOR SOLAR, ISL = 2 FOR INFRARED.                         C
C                                                                      C
C     TAUG: GASEOUS OPTICAL DEPTH                                      C
C     S:    INPUT GAS MIXING RATIO FOR EACH LAYER                      C
C     DP:   AIR MASS PATH FOR A MODEL LAYER (EXLAINED IN RADDRIV).     C
C     DIP:  INTERPRETATION FACTOR FOR PRESSURE BETWEEN TWO NEIGHBORING C
C           STANDARD INPUT DATA PRESSURE LEVELS                        C
C     DT:   LAYER TEMPERATURE - 250 K                                  C
C     INPT: NUMBER OF THE LEVEL FOR THE STANDARD INPUT DATA PRESSURES  C
C----------------------------------------------------------------------C
      IMPLICIT REAL (A-H,O-Z),
     +INTEGER (I-N)
C
      REAL TAUG(ILG,LAY), COEF1(5,N), S(ILG,LAY), DP(ILG,LAY),
     1     DIP(ILG,LAY), DT(ILG,LAY), S1(ILG,LAY)
      INTEGER INPT(ILG,LAY)
      LOGICAL GH
C
      COMMON /TRACE / RMCO2, RMCH4, RMN2O, RMO2, RMF11, RMF12, RMF113,
     1                RMF114
C=======================================================================
      IF (GH)                                                       THEN
        LAY1 =  1
      ELSE
        LAY1 =  LEV1
      ENDIF
      LAY2   =  LAY
C
C     * INITIALIZE TAUG IF IPLUS=2.
C
      IF (IPLUS .EQ. 2) THEN
        DO  50 K = LAY1, LAY2
        DO  50 I = IL1, IL2
          TAUG(I,K)     =  0.
   50   CONTINUE
      ENDIF
C
      IF (NG .EQ. 1)                                                THEN
        DO 100 K = LAY1, LAY2
        DO 100 I = IL1, IL2
          S1(I,K)       =  S(I,K)
  100   CONTINUE       
C      
      ELSE IF (NG .EQ. 3)                                           THEN
        DO 300 K = LAY1, LAY2
        DO 300 I = IL1, IL2
          S1(I,K)       =  RMCO2
  300   CONTINUE
C
      ELSE IF (NG .EQ. 5)                                           THEN
        DO 500 K = LAY1, LAY2
        DO 500 I = IL1, IL2
          S1(I,K)       =  RMN2O
  500   CONTINUE
C
      ELSE IF (NG .EQ. 6)                                           THEN
        DO 600 K = LAY1, LAY2
        DO 600 I = IL1, IL2
          S1(I,K)       =  RMO2
  600   CONTINUE
      ENDIF
C
      DO 2000 K = LAY1, LAY2
        IF (INPT(1,K) .LT. 950)                                     THEN
          DO 1000 I = IL1, IL2  
            M  =  INPT(I,K)
            N  =  M + 1
            X2        =  COEF1(1,N) + DT(I,K) * (COEF1(2,N) + 
     1                   DT(I,K) * (COEF1(3,N) + DT(I,K) * 
     2                  (COEF1(4,N) + DT(I,K) * COEF1(5,N))))
            IF (M .GT. 0)                                           THEN
              X1      =  COEF1(1,M) + DT(I,K) * (COEF1(2,M) + 
     1                   DT(I,K) * (COEF1(3,M) + DT(I,K) * 
     2                  (COEF1(4,M) + DT(I,K) * COEF1(5,M))))
            ELSE
              X1      =  0.0
            ENDIF
C
            TAUG(I,K) =  TAUG(I,K) + (X1 + (X2 - X1) * DIP(I,K)) * 
     1                   S1(I,K) * DP(I,K)
 1000     CONTINUE        
        ELSE
          M  =  INPT(1,K) - 1000
          N  =  M + 1
          DO 1500 I = IL1, IL2
            X2        =  COEF1(1,N) + DT(I,K) * (COEF1(2,N) + 
     1                   DT(I,K) * (COEF1(3,N) + DT(I,K) * 
     2                  (COEF1(4,N) + DT(I,K) * COEF1(5,N))))
            IF (M .GT. 0)                                           THEN
              X1      =  COEF1(1,M) + DT(I,K) * (COEF1(2,M) + 
     1                   DT(I,K) * (COEF1(3,M) + DT(I,K) * 
     2                  (COEF1(4,M) + DT(I,K) * COEF1(5,M))))
            ELSE
              X1      =  0.0
            ENDIF
C
            TAUG(I,K) =  TAUG(I,K) + (X1 + (X2 - X1) * DIP(I,K)) * 
     1                     S1(I,K) * DP(I,K)
 1500     CONTINUE
        ENDIF 
 2000 CONTINUE
C
      RETURN
      END

      SUBROUTINE TLINE2 (TAUG, COEF1, COEF2, S1, S, NG2, DP, DIP, DT,
     1                   INPT, LEV1, GH, N, IL1, IL2, ILG, LAY)
C
C----------------------------------------------------------------------C
C     THE SAME AS TLINE1, BUT WITH TWO MIXED GASES, ONE MUST BE H2O    C
C                                                                      C
C     TAUG: GASEOUS OPTICAL DEPTH                                      C
C     S1:   INPUT H2O MIXING RATIO FOR EACH LAYER                      C
C     S:    INPUT ANOTHER GAS MIXING RATIO FOR EACH LAYER              C
C     DP:   AIR MASS PATH FOR A MODEL LAYER (EXLAINED IN RADDRIV).     C
C     DIP:  INTERPRETATION FACTOR FOR PRESSURE BETWEEN TWO NEIGHBORING C
C           STANDARD INPUT DATA PRESSURE LEVELS                        C
C     DT:   LAYER TEMPERATURE - 250 K                                  C
C     INPT: NUMBER OF THE LEVEL FOR THE STANDARD INPUT DATA PRESSURES  C
C----------------------------------------------------------------------C
      IMPLICIT REAL (A-H,O-Z),
     +INTEGER (I-N)
C
      REAL TAUG(ILG,LAY), COEF1(5,N), COEF2(5,N), S1(ILG,LAY), 
     1     S(ILG,LAY), DP(ILG,LAY), DIP(ILG,LAY), DT(ILG,LAY),
     2     S2(ILG,LAY)
      INTEGER INPT(ILG,LAY)
      LOGICAL GH
C
      COMMON /TRACE / RMCO2, RMCH4, RMN2O, RMO2, RMF11, RMF12, RMF113,
     1                RMF114
C=======================================================================
      IF (GH)                                                       THEN
        LAY1 =  1
      ELSE
        LAY1 =  LEV1
      ENDIF
      LAY2   =  LAY
C
      IF (NG2 .EQ. 2)                                               THEN
        DO 200 K = LAY1, LAY2
        DO 200 I = IL1, IL2
          S2(I,K)     =  S(I,K)
  200   CONTINUE
C
      ELSE IF (NG2 .EQ. 3)                                          THEN
        DO 300 K = LAY1, LAY2
        DO 300 I = IL1, IL2
          S2(I,K)     =  RMCO2
  300   CONTINUE
C
      ELSE IF (NG2 .EQ. 4)                                          THEN
          DO 400 K = LAY1, LAY2
          DO 400 I = IL1, IL2
          S2(I,K)     =  RMCH4
  400   CONTINUE
C
      ELSE IF (NG2 .EQ. 5)                                          THEN
        DO 500 K = LAY1, LAY2
        DO 500 I = IL1, IL2
          S2(I,K)     =  RMN2O
  500   CONTINUE
C
      ELSE IF (NG2 .EQ. 6)                                          THEN
        DO 600 K = LAY1, LAY2
        DO 600 I = IL1, IL2
          S2(I,K)     =  RMO2
  600   CONTINUE
      ENDIF
C
      DO 2000 K = LAY1, LAY2
        IF (INPT(1,K) .LT. 950)                                     THEN
          DO 1000 I = IL1, IL2
            M  =  INPT(I,K)
            N  =  M + 1
            X2        =  COEF1(1,N) + DT(I,K) * (COEF1(2,N) + DT(I,K) * 
     1                  (COEF1(3,N) + DT(I,K) * (COEF1(4,N) +
     2                   DT(I,K) * COEF1(5,N))))
C
            Y2        =  COEF2(1,N) + DT(I,K) * (COEF2(2,N) + DT(I,K) * 
     1                  (COEF2(3,N) + DT(I,K) * (COEF2(4,N) +
     2                   DT(I,K) * COEF2(5,N))))
            IF (M .GT. 0)                                           THEN
              X1      =  COEF1(1,M) + DT(I,K) * (COEF1(2,M) + DT(I,K) *
     1                  (COEF1(3,M) + DT(I,K) * (COEF1(4,M) +
     2                   DT(I,K) * COEF1(5,M))))
C
              Y1      =  COEF2(1,M) + DT(I,K) * (COEF2(2,M) + DT(I,K) * 
     1                  (COEF2(3,M) + DT(I,K) * (COEF2(4,M) +
     2                   DT(I,K) * COEF2(5,M))))
            ELSE
              X1      =  0.0
              Y1      =  0.0
            ENDIF
C
            TAUG(I,K) = ( (X1 + (X2 - X1) * DIP(I,K)) * S1(I,K) +
     1                    (Y1 + (Y2 - Y1) * DIP(I,K)) * S2(I,K) ) * 
     1                     DP(I,K)
 1000     CONTINUE        
        ELSE
          M  =  INPT(1,K) - 1000
          N  =  M + 1
          DO 1500 I = IL1, IL2
            X2        =  COEF1(1,N) + DT(I,K) * (COEF1(2,N) + DT(I,K) *
     1                  (COEF1(3,N) + DT(I,K) * (COEF1(4,N) +
     2                   DT(I,K) * COEF1(5,N))))
C
            Y2        =  COEF2(1,N) + DT(I,K) * (COEF2(2,N) + DT(I,K) *
     1                  (COEF2(3,N) + DT(I,K) * (COEF2(4,N) +
     2                   DT(I,K) * COEF2(5,N))))
            IF (M .GT. 0)                                           THEN
              X1      =  COEF1(1,M) + DT(I,K) * (COEF1(2,M) + DT(I,K) *
     1                  (COEF1(3,M) + DT(I,K) * (COEF1(4,M) +
     2                   DT(I,K) * COEF1(5,M))))
C
              Y1      =  COEF2(1,M) + DT(I,K) * (COEF2(2,M) + DT(I,K) *
     1                  (COEF2(3,M) + DT(I,K) * (COEF2(4,M) +
     2                   DT(I,K) * COEF2(5,M))))
            ELSE
              X1      =  0.0
              Y1      =  0.0
            ENDIF
C
            TAUG(I,K) = ( (X1 + (X2 - X1) * DIP(I,K)) * S1(I,K) +
     1                    (Y1 + (Y2 - Y1) * DIP(I,K)) * S2(I,K) ) *
     1                     DP(I,K)
 1500     CONTINUE
        ENDIF 
 2000 CONTINUE
C
      RETURN
      END

      SUBROUTINE TLINE3 (TAUG, COEF1, COEF2, COEF3, S1, NG2, NG3,
     1                   DP, DIP, DT, INPT, LEV1, GH, N,  
     2                   IL1, IL2, ILG, LAY)
C
C----------------------------------------------------------------------C
C     THE SAME AS TLINEL, BUT WITH THREE MIXED GASES. ONE WITH VARYING C
C     MIXING RATIO THE OTHER TWO WITH CONSTANT MIXING RATIO            C
C                                                                      C
C     TAUG: GASEOUS OPTICAL DEPTH                                      C
C     S1:   INPUT H2O MIXING RATIO FOR EACH LAYER                      C
C     DP:   AIR MASS PATH FOR A MODEL LAYER (EXLAINED IN RADDRIV).     C
C     DIP:  INTERPRETATION FACTOR FOR PRESSURE BETWEEN TWO NEIGHBORING C
C           STANDARD INPUT DATA PRESSURE LEVELS                        C
C     DT:   LAYER TEMPERATURE - 250 K                                  C
C     INPT: NUMBER OF THE LEVEL FOR THE STANDARD INPUT DATA PRESSURES  C
C----------------------------------------------------------------------C
      IMPLICIT REAL (A-H,O-Z),
     +INTEGER (I-N)
C
      REAL TAUG(ILG,LAY), COEF1(5,N), COEF2(5,N), COEF3(5,N), 
     1     S1(ILG,LAY), DP(ILG,LAY), DIP(ILG,LAY), DT(ILG,LAY)
      INTEGER INPT(ILG,LAY)
      LOGICAL GH
C
      COMMON /TRACE / RMCO2, RMCH4, RMN2O, RMO2, RMF11, RMF12, RMF113,
     1                RMF114
C=======================================================================
      IF (NG2 .EQ. 3) S2 =  RMCO2
      IF (NG2 .EQ. 4) S2 =  RMCH4
      IF (NG2 .EQ. 5) S2 =  RMN2O
      IF (NG2 .EQ. 6) S2 =  RMO2
      IF (NG2 .EQ. 7) S2 =  RMF11
      IF (NG2 .EQ. 8) S2 =  RMF12
      iF (NG3 .EQ. 7) S3 =  RMF11
      IF (NG3 .EQ. 8) S3 =  RMF12
 
C
      IF (NG3 .EQ. 3) S3 =  RMCO2
      IF (NG3 .EQ. 4) S3 =  RMCH4
      IF (NG3 .EQ. 5) S3 =  RMN2O
      IF (NG3 .EQ. 6) S3 =  RMO2
C
      IF (GH)                                                       THEN
        LAY1 =  1
      ELSE
        LAY1 =  LEV1
      ENDIF
      LAY2   =  LAY
C
      DO 200 K = LAY1, LAY2
        IF (INPT(1,K) .LT. 950)                                     THEN
          DO 100 I = IL1, IL2
            M  =  INPT(I,K)
            N  =  M + 1
            X2        =  COEF1(1,N) + DT(I,K) * (COEF1(2,N) + DT(I,K) *
     1                  (COEF1(3,N) + DT(I,K) * (COEF1(4,N) +
     2                   DT(I,K) * COEF1(5,N))))
C
            Y2        =  COEF2(1,N) + DT(I,K) * (COEF2(2,N) + DT(I,K) *
     1                  (COEF2(3,N) + DT(I,K) * (COEF2(4,N) +
     2                   DT(I,K) * COEF2(5,N))))
C
            Z2        =  COEF3(1,N) + DT(I,K) * (COEF3(2,N) + DT(I,K) * 
     1                  (COEF3(3,N) + DT(I,K) * (COEF3(4,N) +
     2                   DT(I,K) * COEF3(5,N))))
            IF (M .GT. 0)                                           THEN
              X1      =  COEF1(1,M) + DT(I,K) * (COEF1(2,M) + DT(I,K) *
     1                  (COEF1(3,M) + DT(I,K) * (COEF1(4,M) +
     2                   DT(I,K) * COEF1(5,M))))
C
              Y1      =  COEF2(1,M) + DT(I,K) * (COEF2(2,M) + DT(I,K) *
     1                  (COEF2(3,M) + DT(I,K) * (COEF2(4,M) +
     2                   DT(I,K) * COEF2(5,M))))
C
              Z1      =  COEF3(1,M) + DT(I,K) * (COEF3(2,M) + DT(I,K) *
     1                  (COEF3(3,M) + DT(I,K) * (COEF3(4,M) +
     2                   DT(I,K) * COEF3(5,M))))
            ELSE
              X1      =  0.0
              Y1      =  0.0
              Z1      =  0.0
            ENDIF
C
            TAUG(I,K) = ( (X1 + (X2 - X1) * DIP(I,K)) * S1(I,K) +
     1                    (Y1 + (Y2 - Y1) * DIP(I,K)) * S2 +
     2                    (Z1 + (Z2 - Z1) * DIP(I,K)) * S3  ) * DP(I,K)
  100     CONTINUE        
        ELSE
          M  =  INPT(1,K) - 1000
          N  =  M + 1
          DO 150 I = IL1, IL2
            X2        =  COEF1(1,N) + DT(I,K) * (COEF1(2,N) + DT(I,K) * 
     1                  (COEF1(3,N) + DT(I,K) * (COEF1(4,N) +
     2                   DT(I,K) * COEF1(5,N))))
C
            Y2        =  COEF2(1,N) + DT(I,K) * (COEF2(2,N) + DT(I,K) * 
     1                  (COEF2(3,N) + DT(I,K) * (COEF2(4,N) +
     2                   DT(I,K) * COEF2(5,N))))
C
            Z2        =  COEF3(1,N) + DT(I,K) * (COEF3(2,N) + DT(I,K) * 
     1                  (COEF3(3,N) + DT(I,K) * (COEF3(4,N) +
     2                   DT(I,K) * COEF3(5,N))))
            IF (M .GT. 0)                                           THEN
              X1      =  COEF1(1,M) + DT(I,K) * (COEF1(2,M) + DT(I,K) *
     1                  (COEF1(3,M) + DT(I,K) * (COEF1(4,M) +
     2                   DT(I,K) * COEF1(5,M))))
C
              Y1      =  COEF2(1,M) + DT(I,K) * (COEF2(2,M) + DT(I,K) *
     1                  (COEF2(3,M) + DT(I,K) * (COEF2(4,M) +
     2                   DT(I,K) * COEF2(5,M))))
C
              Z1      =  COEF3(1,M) + DT(I,K) * (COEF3(2,M) + DT(I,K) *
     1                  (COEF3(3,M) + DT(I,K) * (COEF3(4,M) +
     2                   DT(I,K) * COEF3(5,M))))
            ELSE
              X1      =  0.0
              Y1      =  0.0
              Z1      =  0.0
            ENDIF
C
            TAUG(I,K) = ( (X1 + (X2 - X1) * DIP(I,K)) * S1(I,K) +
     1                    (Y1 + (Y2 - Y1) * DIP(I,K)) * S2 +
     2                    (Z1 + (Z2 - Z1) * DIP(I,K)) * S3  ) * DP(I,K)

  150     CONTINUE
        ENDIF 
  200 CONTINUE
C
      RETURN
      END

      SUBROUTINE TLINEHC (TAUG, COEF1U, COEF1D, COEF2U, COEF2D, S, DP, 
     1                    DIP, DIR, DT, INPTR, INPT, LEV1,
     2                    IL1, IL2, ILG, LAY)
C
C----------------------------------------------------------------------C
C     THIS SUBROUTINE DETERMINES THE OPTICAL DEPTH FOR H2O AND CO2 IN  C
C     THE REGION OF 540-800 CM^-1                                      C
C                                                                      C
C     TAUG: GASEOUS OPTICAL DEPTH                                      C
C     S:    INPUT H2O MIXING RATIO FOR EACH LAYER                      C
C     DP:   AIR MASS PATH FOR A MODEL LAYER (EXLAINED IN RADDRIV).     C
C     DIP:  INTERPRETATION FACTOR FOR PRESSURE BETWEEN TWO NEIGHBORING C
C           STANDARD INPUT DATA PRESSURE LEVELS                        C
C     DIR:  INTERPRETATION FACTOR FOR MASS RATIO OF H2O / CO2 BETWEEN  C
C           TWO NEIGHBORING STANDARD INPUT RATIOS                      C
C     DT:   LAYER TEMPERATURE - 250 K                                  C
C     INPR: NUMBER OF THE RATIO LEVEL FOR THE STANDARD 5 RATIOS        C
C     INPT: NUMBER OF THE LEVEL FOR THE STANDARD INPUT DATA PRESSURES  C
C----------------------------------------------------------------------C
      IMPLICIT REAL (A-H,O-Z),
     +INTEGER (I-N)
C
      REAL TAUG(ILG,LAY), COEF1U(5,11), COEF1D(5,5,7), 
     1                    COEF2U(5,11), COEF2D(5,5,7)
      REAL S(ILG,LAY), DP(ILG,LAY), DIP(ILG,LAY), DIR(ILG,LAY), 
     1     DT(ILG,LAY)
      INTEGER INPTR(ILG,LAY), INPT(ILG,LAY)
C
      COMMON /TRACE / RMCO2, RMCH4, RMN2O, RMO2, RMF11, RMF12, RMF113,
     1                RMF114
C=======================================================================
      DO 300 K = LEV1, LAY
       IF (INPT(1,K) .LT. 950)                                      THEN
        DO 100 I = IL1, IL2
         M =  INPT(I,K)
         IF (M .LE. 11)                                             THEN
          N     =  M + 1
          IF (M .GT. 0)                                             THEN
            X1  =  COEF1U(1,M) + DT(I,K) * (COEF1U(2,M) + DT(I,K) * 
     1            (COEF1U(3,M) + DT(I,K) * (COEF1U(4,M) +
     2             DT(I,K) * COEF1U(5,M))))
            Y1  =  COEF2U(1,M) + DT(I,K) * (COEF2U(2,M) + DT(I,K) * 
     1            (COEF2U(3,M) + DT(I,K) * (COEF2U(4,M) +
     2             DT(I,K) * COEF2U(5,M))))
          ELSE
            X1  =  0.0
            Y1  =  0.0
          ENDIF
          IF (M .LT. 11)                                            THEN
            X2  =  COEF1U(1,N) + DT(I,K) * (COEF1U(2,N) + DT(I,K) * 
     1            (COEF1U(3,N) + DT(I,K) * (COEF1U(4,N) +
     2             DT(I,K) * COEF1U(5,N))))
            Y2  =  COEF2U(1,N) + DT(I,K) * (COEF2U(2,N) + DT(I,K) * 
     1            (COEF2U(3,N) + DT(I,K) * (COEF2U(4,N) +
     2             DT(I,K) * COEF2U(5,N))))
          ELSE
            X2  =  COEF1D(1,1,1) + DT(I,K) * (COEF1D(2,1,1) + DT(I,K) *
     1            (COEF1D(3,1,1) + DT(I,K) * (COEF1D(4,1,1) +
     2             DT(I,K) * COEF1D(5,1,1))))
            Y2  =  COEF2D(1,1,1) + DT(I,K) * (COEF2D(2,1,1) + DT(I,K) *
     1            (COEF2D(3,1,1) + DT(I,K) * (COEF2D(4,1,1) +
     2             DT(I,K) * COEF2D(5,1,1))))
          ENDIF
         ELSE
          J     =  M - 11
          N     =  J + 1
          L     =  INPTR(I,K)
          IF (L .LT. 1)                                             THEN
            X1  =  COEF1D(1,1,J) + DT(I,K) * (COEF1D(2,1,J) + DT(I,K) * 
     1            (COEF1D(3,1,J) + DT(I,K) * (COEF1D(4,1,J) +
     2             DT(I,K) * COEF1D(5,1,J))))
            X2  =  COEF1D(1,1,N) + DT(I,K) * (COEF1D(2,1,N) + DT(I,K) * 
     1            (COEF1D(3,1,N) + DT(I,K) * (COEF1D(4,1,N) +
     2             DT(I,K) * COEF1D(5,1,N))))
C
            Y1  =  COEF2D(1,1,J) + DT(I,K) * (COEF2D(2,1,J) + DT(I,K) * 
     1            (COEF2D(3,1,J) + DT(I,K) * (COEF2D(4,1,J) +
     2             DT(I,K) * COEF2D(5,1,J))))
            Y2  =  COEF2D(1,1,N) + DT(I,K) * (COEF2D(2,1,N) + DT(I,K) * 
     1            (COEF2D(3,1,N) + DT(I,K) * (COEF2D(4,1,N) +
     2             DT(I,K) * COEF2D(5,1,N))))
C
          ELSE IF (L .LT. 5)                                        THEN
            LP1 =  L + 1
            X11 =  COEF1D(1,L,J) + DT(I,K) * (COEF1D(2,L,J) + DT(I,K) * 
     1            (COEF1D(3,L,J) + DT(I,K) * (COEF1D(4,L,J) +
     2             DT(I,K) * COEF1D(5,L,J))))
            X21 =  COEF1D(1,L,N) + DT(I,K) * (COEF1D(2,L,N) + DT(I,K) * 
     1            (COEF1D(3,L,N) + DT(I,K) * (COEF1D(4,L,N) +
     2             DT(I,K) * COEF1D(5,L,N))))
C
            Y11 =  COEF2D(1,L,J) + DT(I,K) * (COEF2D(2,L,J) + DT(I,K) * 
     1            (COEF2D(3,L,J) + DT(I,K) * (COEF2D(4,L,J) +
     2             DT(I,K) * COEF2D(5,L,J))))
            Y21 =  COEF2D(1,L,N) + DT(I,K) * (COEF2D(2,L,N) + DT(I,K) * 
     1            (COEF2D(3,L,N) + DT(I,K) * (COEF2D(4,L,N) +
     2             DT(I,K) * COEF2D(5,L,N))))
C
            X12 =  COEF1D(1,LP1,J) + DT(I,K) * (COEF1D(2,LP1,J) +
     1                               DT(I,K) * (COEF1D(3,LP1,J) + 
     2                               DT(I,K) * (COEF1D(4,LP1,J) + 
     3                               DT(I,K) * COEF1D(5,LP1,J))))
            X22 =  COEF1D(1,LP1,N) + DT(I,K) * (COEF1D(2,LP1,N) +
     1                               DT(I,K) * (COEF1D(3,LP1,N) + 
     2                               DT(I,K) * (COEF1D(4,LP1,N) + 
     3                               DT(I,K) * COEF1D(5,LP1,N))))
C
            Y12 =  COEF2D(1,LP1,J) + DT(I,K) * (COEF2D(2,LP1,J) +
     1                               DT(I,K) * (COEF2D(3,LP1,J) + 
     2                               DT(I,K) * (COEF2D(4,LP1,J) + 
     3                               DT(I,K) * COEF2D(5,LP1,J))))
            Y22 =  COEF2D(1,LP1,N) + DT(I,K) * (COEF2D(2,LP1,N) +
     1                               DT(I,K) * (COEF2D(3,LP1,N) + 
     2                               DT(I,K) * (COEF2D(4,LP1,N) + 
     3                               DT(I,K) * COEF2D(5,LP1,N))))
C
            X1  =  X11 + (X12 - X11) * DIR(I,K)
            X2  =  X21 + (X22 - X21) * DIR(I,K)
            Y1  =  Y11 + (Y12 - Y11) * DIR(I,K)
            Y2  =  Y21 + (Y22 - Y21) * DIR(I,K)
          ELSE
            X1  =  COEF1D(1,5,J) + DT(I,K) * (COEF1D(2,5,J) + DT(I,K) * 
     1            (COEF1D(3,5,J) + DT(I,K) * (COEF1D(4,5,J) +
     2             DT(I,K) * COEF1D(5,5,J))))
            X2  =  COEF1D(1,5,N) + DT(I,K) * (COEF1D(2,5,N) + DT(I,K) * 
     1            (COEF1D(3,5,N) + DT(I,K) * (COEF1D(4,5,N) +
     2             DT(I,K) * COEF1D(5,5,N))))
            Y1  =  COEF2D(1,5,J) + DT(I,K) * (COEF2D(2,5,J) + DT(I,K) * 
     1            (COEF2D(3,5,J) + DT(I,K) * (COEF2D(4,5,J) +
     2             DT(I,K) * COEF2D(5,5,J))))
            Y2  =  COEF2D(1,5,N) + DT(I,K) * (COEF2D(2,5,N) + DT(I,K) *
     1            (COEF2D(3,5,N) + DT(I,K) * (COEF2D(4,5,N) +
     2             DT(I,K) * COEF2D(5,5,N))))
          ENDIF
         ENDIF
C
         TAUG(I,K) = ( (X1 + (X2 - X1) * DIP(I,K)) * S(I,K) +
     1                 (Y1 + (Y2 - Y1) * DIP(I,K)) * RMCO2 ) * DP(I,K)
  100   CONTINUE
C
       ELSE
        M =  INPT(1,K) - 1000
        DO 200 I = IL1, IL2
         IF (M .LE. 11)                                             THEN
          N     =  M + 1
          IF (M .GT. 0)                                             THEN
            X1  =  COEF1U(1,M) + DT(I,K) * (COEF1U(2,M) + DT(I,K) * 
     1            (COEF1U(3,M) + DT(I,K) * (COEF1U(4,M) +
     2             DT(I,K) * COEF1U(5,M))))
            Y1  =  COEF2U(1,M) + DT(I,K) * (COEF2U(2,M) + DT(I,K) * 
     1            (COEF2U(3,M) + DT(I,K) * (COEF2U(4,M) +
     2             DT(I,K) * COEF2U(5,M))))
          ELSE
            X1  =  0.0
            Y1  =  0.0
          ENDIF
          IF (M .LT. 11)                                            THEN
            X2  =  COEF1U(1,N) + DT(I,K) * (COEF1U(2,N) + DT(I,K) * 
     1            (COEF1U(3,N) + DT(I,K) * (COEF1U(4,N) +
     2             DT(I,K) * COEF1U(5,N))))
            Y2  =  COEF2U(1,N) + DT(I,K) * (COEF2U(2,N) + DT(I,K) * 
     1            (COEF2U(3,N) + DT(I,K) * (COEF2U(4,N) +
     2             DT(I,K) * COEF2U(5,N))))
          ELSE
            X2  =  COEF1D(1,1,1) + DT(I,K) * (COEF1D(2,1,1) + DT(I,K) *
     1            (COEF1D(3,1,1) + DT(I,K) * (COEF1D(4,1,1) +
     2             DT(I,K) * COEF1D(5,1,1))))
            Y2  =  COEF2D(1,1,1) + DT(I,K) * (COEF2D(2,1,1) + DT(I,K) *
     1            (COEF2D(3,1,1) + DT(I,K) * (COEF2D(4,1,1) +
     2             DT(I,K) * COEF2D(5,1,1))))
          ENDIF
         ELSE
          J     =  M - 11
          N     =  J + 1
          L     =  INPTR(I,K)
          IF (L .LT. 1)                                             THEN
            X1  =  COEF1D(1,1,J) + DT(I,K) * (COEF1D(2,1,J) + DT(I,K) * 
     1            (COEF1D(3,1,J) + DT(I,K) * (COEF1D(4,1,J) +
     2             DT(I,K) * COEF1D(5,1,J))))
            X2  =  COEF1D(1,1,N) + DT(I,K) * (COEF1D(2,1,N) + DT(I,K) * 
     1            (COEF1D(3,1,N) + DT(I,K) * (COEF1D(4,1,N) +
     2             DT(I,K) * COEF1D(5,1,N))))
C
            Y1  =  COEF2D(1,1,J) + DT(I,K) * (COEF2D(2,1,J) + DT(I,K) * 
     1            (COEF2D(3,1,J) + DT(I,K) * (COEF2D(4,1,J) +
     2             DT(I,K) * COEF2D(5,1,J))))
            Y2  =  COEF2D(1,1,N) + DT(I,K) * (COEF2D(2,1,N) + DT(I,K) * 
     1            (COEF2D(3,1,N) + DT(I,K) * (COEF2D(4,1,N) +
     2             DT(I,K) * COEF2D(5,1,N))))
C
          ELSE IF (L .LT. 5)                                        THEN
            LP1 =  L + 1
            X11 =  COEF1D(1,L,J) + DT(I,K) * (COEF1D(2,L,J) + DT(I,K) * 
     1            (COEF1D(3,L,J) + DT(I,K) * (COEF1D(4,L,J) +
     2             DT(I,K) * COEF1D(5,L,J))))
            X21 =  COEF1D(1,L,N) + DT(I,K) * (COEF1D(2,L,N) + DT(I,K) * 
     1            (COEF1D(3,L,N) + DT(I,K) * (COEF1D(4,L,N) +
     2             DT(I,K) * COEF1D(5,L,N))))
C
            Y11 =  COEF2D(1,L,J) + DT(I,K) * (COEF2D(2,L,J) + DT(I,K) * 
     1            (COEF2D(3,L,J) + DT(I,K) * (COEF2D(4,L,J) +
     2             DT(I,K) * COEF2D(5,L,J))))
            Y21 =  COEF2D(1,L,N) + DT(I,K) * (COEF2D(2,L,N) + DT(I,K) * 
     1            (COEF2D(3,L,N) + DT(I,K) * (COEF2D(4,L,N) +
     2             DT(I,K) * COEF2D(5,L,N))))
C
            X12 =  COEF1D(1,LP1,J) + DT(I,K) * (COEF1D(2,LP1,J) +
     1                               DT(I,K) * (COEF1D(3,LP1,J) + 
     2                               DT(I,K) * (COEF1D(4,LP1,J) + 
     3                               DT(I,K) * COEF1D(5,LP1,J))))
            X22 =  COEF1D(1,LP1,N) + DT(I,K) * (COEF1D(2,LP1,N) +
     1                               DT(I,K) * (COEF1D(3,LP1,N) + 
     2                               DT(I,K) * (COEF1D(4,LP1,N) + 
     3                               DT(I,K) * COEF1D(5,LP1,N))))
C
            Y12 =  COEF2D(1,LP1,J) + DT(I,K) * (COEF2D(2,LP1,J) +
     1                               DT(I,K) * (COEF2D(3,LP1,J) + 
     2                               DT(I,K) * (COEF2D(4,LP1,J) + 
     3                               DT(I,K) * COEF2D(5,LP1,J))))
            Y22 =  COEF2D(1,LP1,N) + DT(I,K) * (COEF2D(2,LP1,N) +
     1                               DT(I,K) * (COEF2D(3,LP1,N) + 
     2                               DT(I,K) * (COEF2D(4,LP1,N) + 
     3                               DT(I,K) * COEF2D(5,LP1,N))))
C
            X1  =  X11 + (X12 - X11) * DIR(I,K)
            X2  =  X21 + (X22 - X21) * DIR(I,K)
            Y1  =  Y11 + (Y12 - Y11) * DIR(I,K)
            Y2  =  Y21 + (Y22 - Y21) * DIR(I,K)
          ELSE
            X1  =  COEF1D(1,5,J) + DT(I,K) * (COEF1D(2,5,J) + DT(I,K) * 
     1            (COEF1D(3,5,J) + DT(I,K) * (COEF1D(4,5,J) +
     2             DT(I,K) * COEF1D(5,5,J))))
            X2  =  COEF1D(1,5,N) + DT(I,K) * (COEF1D(2,5,N) + DT(I,K) * 
     1            (COEF1D(3,5,N) + DT(I,K) * (COEF1D(4,5,N) +
     2             DT(I,K) * COEF1D(5,5,N))))
            Y1  =  COEF2D(1,5,J) + DT(I,K) * (COEF2D(2,5,J) + DT(I,K) * 
     1            (COEF2D(3,5,J) + DT(I,K) * (COEF2D(4,5,J) +
     2             DT(I,K) * COEF2D(5,5,J))))
            Y2  =  COEF2D(1,5,N) + DT(I,K) * (COEF2D(2,5,N) + DT(I,K) *
     1            (COEF2D(3,5,N) + DT(I,K) * (COEF2D(4,5,N) +
     2             DT(I,K) * COEF2D(5,5,N))))
          ENDIF
         ENDIF
C
         TAUG(I,K) = ( (X1 + (X2 - X1) * DIP(I,K)) * S(I,K) +
     1                 (Y1 + (Y2 - Y1) * DIP(I,K)) * RMCO2 ) * DP(I,K)
  200   CONTINUE
       ENDIF
  300 CONTINUE

      RETURN
      END

      SUBROUTINE TCONTL1 (TAUG, COEF1, COEF2, S, DP, DIP, DT, LC, INPT,
     1                    MCONT, GH, IL1, IL2, ILG, LAY)
C
C----------------------------------------------------------------------C
C     INFRARED WATER VAPOR CONTINUUM, COEF1 IS THE COEFFICIENT FOR     C
C     SELF, COEF2 IS THE COEFFICIENT FOR FOREIGN. THE CONTINUUM ONLY   C
C     APPLIES TO THE LAYERS BELOW 138.9440 MB OR EVEN LOWER REGION     C
C     DEPENDING ON EACH BAND. LC IS NUMBER OF LEVEL FOR STANDARD       C
C     PRESSURE CONSIDERED IN CALCULATING THE CONTINUUM.                C
C     1.608 = 28.97 / 18.016, A FCTR FOR WATER VAPOR PARTIAL PRESSURE  C
C                                                                      C
C     TAUG:  GASEOUS OPTICAL DEPTH                                     C
C     S:     INPUT H2O MIXING RATIO FOR EACH LAYER                     C
C     DP:    AIR MASS PATH FOR A MODEL LAYER (EXLAINED IN RADDRIV).    C
C     DIP:   INTERPRETATION FACTOR FOR PRESSURE BETWEEN TWO            C
C            NEIGHBORING STANDARD INPUT DATA PRESSURE LEVELS           C
C     DT:    LAYER TEMPERATURE - 250 K                                 C
C     INPT:  NUMBER OF THE LEVEL FOR THE STANDARD INPUT DATA PRESSURES C
C     MCONT: THE HIGHEST LEVEL FOR WATER VAPOR CONTINUUM CALCULATION   C
C----------------------------------------------------------------------C
      IMPLICIT REAL (A-H,O-Z),
     +INTEGER (I-N)
C
      REAL TAUG(ILG,LAY), COEF1(5,LC), COEF2(5,LC)
      REAL S(ILG,LAY), DP(ILG,LAY), DIP(ILG,LAY), DT(ILG,LAY)
      INTEGER INPT(ILG,LAY)
      LOGICAL GH
C=======================================================================
      IF (GH)                                                       THEN
        NC =  29 - LC
      ELSE
        NC =  19 - LC
      ENDIF
C 
      DO 200 K = MCONT, LAY
        IF (INPT(1,K) .LT. 950)                                     THEN
          DO 100 I = IL1, IL2
            J =  INPT(I,K)
            IF (J .GE. NC)                                          THEN
              M  =  J - NC + 1
              N  =  M + 1
              X1        =  COEF1(1,M) + DT(I,K) * (COEF1(2,M) +
     1                     DT(I,K) * (COEF1(3,M) + DT(I,K) *
     2                    (COEF1(4,M) + DT(I,K) * COEF1(5,M))))
C
              X2        =  COEF1(1,N) + DT(I,K) * (COEF1(2,N) +
     1                     DT(I,K) * (COEF1(3,N) + DT(I,K) *
     2                    (COEF1(4,N) + DT(I,K) * COEF1(5,N))))
C
              Y1        =  COEF2(1,M) + DT(I,K) * (COEF2(2,M) +
     1                     DT(I,K) * (COEF2(3,M) + DT(I,K) *
     2                    (COEF2(4,M) + DT(I,K) * COEF2(5,M))))
C
              Y2        =  COEF2(1,N) + DT(I,K) * (COEF2(2,N) +
     1                     DT(I,K) * (COEF2(3,N) + DT(I,K) *
     2                    (COEF2(4,N) + DT(I,K) * COEF2(5,N))))
C
              TAUG(I,K) =  TAUG(I,K) +
     1                    ( (X1 - Y1 + (X2 - X1 - Y2 + Y1) *
     2                     DIP(I,K)) * 1.608 * S(I,K) +
     3                     Y1 + (Y2 - Y1) * DIP(I,K) ) *
     4                     S(I,K) * DP(I,K)
            ENDIF
  100     CONTINUE
        ELSE
          J  =  INPT(1,K) - 1000
          M  =  J - NC + 1
          N  =  M + 1
          DO 150 I = IL1, IL2
            IF (J .GE. NC)                                          THEN
              X1        =  COEF1(1,M) + DT(I,K) * (COEF1(2,M) +
     1                     DT(I,K) * (COEF1(3,M) + DT(I,K) *
     2                    (COEF1(4,M) + DT(I,K) * COEF1(5,M))))
C
              X2        =  COEF1(1,N) + DT(I,K) * (COEF1(2,N) +
     1                     DT(I,K) * (COEF1(3,N) + DT(I,K) *
     2                    (COEF1(4,N) + DT(I,K) * COEF1(5,N))))
C
              Y1        =  COEF2(1,M) + DT(I,K) * (COEF2(2,M) +
     1                     DT(I,K) * (COEF2(3,M) + DT(I,K) *
     2                    (COEF2(4,M) + DT(I,K) * COEF2(5,M))))
C
              Y2        =  COEF2(1,N) + DT(I,K) * (COEF2(2,N) +
     1                     DT(I,K) * (COEF2(3,N) + DT(I,K) *
     2                    (COEF2(4,N) + DT(I,K) * COEF2(5,N))))
C
              TAUG(I,K) =  TAUG(I,K) +
     1                    ( (X1 - Y1 + (X2 - X1 - Y2 + Y1) *
     2                     DIP(I,K)) * 1.608 * S(I,K) +
     3                     Y1 + (Y2 - Y1) * DIP(I,K) ) *
     4                     S(I,K) * DP(I,K)
            ENDIF
  150     CONTINUE
        ENDIF
  200 CONTINUE
C
      RETURN
      END

      SUBROUTINE TCONTHL (TAUG, COEF1, COEF2, S, DP, DIP, DIR, DT, 
     1                    INPTR, INPT, MCONT, IL1, IL2, ILG, LAY)
C
C----------------------------------------------------------------------C
C     WATER VAPOR CONTINUUM FOR 540-800 CM-1. DIFFERENT FROM TCONTL,   C
C     VARIATION OF MASS MIXING RATIO FOR H2O AND CO2 IS CONSIDER.      C
C     LC = 4, BUT THE INPTUT DATA ARE WITH 6 GROUP, 2 OF THEM ARE USED C
C     FOR MASS MIXING RATIO CHANGES                                    C
C                                                                      C
C     TAUG:  GASEOUS OPTICAL DEPTH                                     C
C     S:     INPUT H2O MIXING RATIO FOR EACH LAYER                     C
C     DP:    AIR MASS PATH FOR A MODEL LAYER (EXLAINED IN RADDRIV).    C
C     DIP:   INTERPRETATION FACTOR FOR PRESSURE BETWEEN TWO            C
C            NEIGHBORING STANDARD INPUT DATA PRESSURE LEVELS           C
C     DT:    LAYER TEMPERATURE - 250 K                                 C
C     INPT:  NUMBER OF THE LEVEL FOR THE STANDARD INPUT DATA PRESSURES C
C     MCONT: THE HIGHEST LEVEL FOR WATER VAPOR CONTINUUM CALCULATION   C
C----------------------------------------------------------------------C
      IMPLICIT REAL (A-H,O-Z),
     +INTEGER (I-N)
C
      REAL TAUG(ILG,LAY), COEF1(5,5,4), COEF2(5,5,4)
C
      REAL S(ILG,LAY), DP(ILG,LAY), DIP(ILG,LAY), DIR(ILG,LAY),
     1     DT(ILG,LAY)
      INTEGER INPTR(ILG,LAY), INPT(ILG,LAY)
C
      COMMON /TRACE / RMCO2, RMCH4, RMN2O, RMO2, RMF11, RMF12, RMF113,
     1                RMF114 
C
      DATA RMIN,RMAX /0.5, 16.0/
C=======================================================================
      RRMCO2 =  1. / (RMCO2 + 1.E-10)
C
      DO 200 K = MCONT, LAY
        IF (INPT(1,K) .LT. 950)                                     THEN
          DO 100 I = IL1, IL2
            M =  INPT(I,K) - 14
            N =  M + 1
            IF (M .GE. 1)                                           THEN
              L   =  INPTR(I,K)
              IF (L .LT. 1)                                         THEN
                X1  =  COEF1(1,1,M) + DT(I,K) * (COEF1(2,1,M) + 
     1                                DT(I,K) * (COEF1(3,1,M) + 
     2                                DT(I,K) * (COEF1(4,1,M) +
     3                                DT(I,K) * COEF1(5,1,M))))
                X2  =  COEF1(1,1,N) + DT(I,K) * (COEF1(2,1,N) + 
     1                                DT(I,K) * (COEF1(3,1,N) + 
     2                                DT(I,K) * (COEF1(4,1,N) +
     3                                DT(I,K) * COEF1(5,1,N))))
C
                Y1  =  COEF2(1,1,M) + DT(I,K) * (COEF2(2,1,M) + 
     1                                DT(I,K) * (COEF2(3,1,M) + 
     2                                DT(I,K) * (COEF2(4,1,M) +
     3                                DT(I,K) * COEF2(5,1,M))))
                Y2  =  COEF2(1,1,N) + DT(I,K) * (COEF2(2,1,N) +
     1                                DT(I,K) * (COEF2(3,1,N) + 
     2                                DT(I,K) * (COEF2(4,1,N) +
     3                                DT(I,K) * COEF2(5,1,N))))
C
              ELSE IF (L .LT. 5)                                    THEN
                LP1 =  L + 1
                X11 =  COEF1(1,L,M) + DT(I,K) * (COEF1(2,L,M) + 
     1                                DT(I,K) * (COEF1(3,L,M) +
     2                                DT(I,K) * (COEF1(4,L,M) +
     3                                DT(I,K) * COEF1(5,L,M))))
                X21 =  COEF1(1,L,N) + DT(I,K) * (COEF1(2,L,N) + 
     1                                DT(I,K) * (COEF1(3,L,N) +
     2                                DT(I,K) * (COEF1(4,L,N) +
     3                                DT(I,K) * COEF1(5,L,N))))
C
                Y11 =  COEF2(1,L,M) + DT(I,K) * (COEF2(2,L,M) + 
     1                                DT(I,K) * (COEF2(3,L,M) +
     2                                DT(I,K) * (COEF2(4,L,M) +
     3                                DT(I,K) * COEF2(5,L,M))))
                Y21 =  COEF2(1,L,N) + DT(I,K) * (COEF2(2,L,N) +
     1                                DT(I,K) * (COEF2(3,L,N) +
     2                                DT(I,K) * (COEF2(4,L,N) +
     3                                DT(I,K) * COEF2(5,L,N))))
C
                X12 =  COEF1(1,LP1,M) + DT(I,K) * (COEF1(2,LP1,M) +
     1                                  DT(I,K) * (COEF1(3,LP1,M) +
     2                                  DT(I,K) * (COEF1(4,LP1,M) +
     3                                  DT(I,K) * COEF1(5,LP1,M))))
                X22 =  COEF1(1,LP1,N) + DT(I,K) * (COEF1(2,LP1,N) +
     1                                  DT(I,K) * (COEF1(3,LP1,N) +
     2                                  DT(I,K) * (COEF1(4,LP1,N) +
     3                                  DT(I,K) * COEF1(5,LP1,N))))
C
                Y12 =  COEF2(1,LP1,M) + DT(I,K) * (COEF2(2,LP1,M) +
     1                                  DT(I,K) * (COEF2(3,LP1,M) +
     2                                  DT(I,K) * (COEF2(4,LP1,M) +
     3                                  DT(I,K) * COEF2(5,LP1,M))))
                Y22 =  COEF2(1,LP1,N) + DT(I,K) * (COEF2(2,LP1,N) +
     1                                  DT(I,K) * (COEF2(3,LP1,N) +
     2                                  DT(I,K) * (COEF2(4,LP1,N) +
     3                                  DT(I,K) * COEF2(5,LP1,N))))
C
                X1  =  X11 + (X12 - X11) * DIR(I,K)
                X2  =  X21 + (X22 - X21) * DIR(I,K)
                Y1  =  Y11 + (Y12 - Y11) * DIR(I,K)
                Y2  =  Y21 + (Y22 - Y21) * DIR(I,K)
              ELSE
                X1  =  COEF1(1,5,M) + DT(I,K) * (COEF1(2,5,M) + 
     1                                DT(I,K) * (COEF1(3,5,M) + 
     2                                DT(I,K) * (COEF1(4,5,M) +
     3                                DT(I,K) * COEF1(5,5,M))))
                X2  =  COEF1(1,5,N) + DT(I,K) * (COEF1(2,5,N) + 
     1                                DT(I,K) * (COEF1(3,5,N) + 
     2                                DT(I,K) * (COEF1(4,5,N) +
     3                                DT(I,K) * COEF1(5,5,N))))
                Y1  =  COEF2(1,5,M) + DT(I,K) * (COEF2(2,5,M) + 
     1                                DT(I,K) * (COEF2(3,5,M) + 
     2                                DT(I,K) * (COEF2(4,5,M) +
     3                                DT(I,K) * COEF2(5,5,M))))
                Y2  =  COEF2(1,5,N) + DT(I,K) * (COEF2(2,5,N) +
     1                                DT(I,K) * (COEF2(3,5,N) + 
     2                                DT(I,K) * (COEF2(4,5,N) +
     3                                DT(I,K) * COEF2(5,5,N))))
              ENDIF
                TAUG(I,K) =  TAUG(I,K) + ((X1 - Y1 + 
     1                      (X2 - X1 - Y2 + Y1) * DIP(I,K)) * 1.608 * 
     2                       S(I,K) + Y1 + (Y2 - Y1) * DIP(I,K)) * 
     3                       S(I,K) * DP(I,K)
            ENDIF
  100     CONTINUE
        ELSE
          J  =  INPT(1,K) - 1000
          M =  J - 14
          N =  M + 1
          DO 150 I = IL1, IL2
            IF (M .GE. 1)                                           THEN
              L   =  INPTR(I,K)
              IF (L .LT. 1)                                         THEN
                X1  =  COEF1(1,1,M) + DT(I,K) * (COEF1(2,1,M) + 
     1                                DT(I,K) * (COEF1(3,1,M) + 
     2                                DT(I,K) * (COEF1(4,1,M) +
     3                                DT(I,K) * COEF1(5,1,M))))
                X2  =  COEF1(1,1,N) + DT(I,K) * (COEF1(2,1,N) + 
     1                                DT(I,K) * (COEF1(3,1,N) + 
     2                                DT(I,K) * (COEF1(4,1,N) +
     3                                DT(I,K) * COEF1(5,1,N))))
C
                Y1  =  COEF2(1,1,M) + DT(I,K) * (COEF2(2,1,M) + 
     1                                DT(I,K) * (COEF2(3,1,M) + 
     2                                DT(I,K) * (COEF2(4,1,M) +
     3                                DT(I,K) * COEF2(5,1,M))))
                Y2  =  COEF2(1,1,N) + DT(I,K) * (COEF2(2,1,N) +
     1                                DT(I,K) * (COEF2(3,1,N) + 
     2                                DT(I,K) * (COEF2(4,1,N) +
     3                                DT(I,K) * COEF2(5,1,N))))
C
              ELSE IF (L .LT. 5)                                    THEN
                LP1 =  L + 1
                X11 =  COEF1(1,L,M) + DT(I,K) * (COEF1(2,L,M) + 
     1                                DT(I,K) * (COEF1(3,L,M) +
     2                                DT(I,K) * (COEF1(4,L,M) +
     3                                DT(I,K) * COEF1(5,L,M))))
                X21 =  COEF1(1,L,N) + DT(I,K) * (COEF1(2,L,N) + 
     1                                DT(I,K) * (COEF1(3,L,N) +
     2                                DT(I,K) * (COEF1(4,L,N) +
     3                                DT(I,K) * COEF1(5,L,N))))
C
                Y11 =  COEF2(1,L,M) + DT(I,K) * (COEF2(2,L,M) + 
     1                                DT(I,K) * (COEF2(3,L,M) +
     2                                DT(I,K) * (COEF2(4,L,M) +
     3                                DT(I,K) * COEF2(5,L,M))))
                Y21 =  COEF2(1,L,N) + DT(I,K) * (COEF2(2,L,N) +
     1                                DT(I,K) * (COEF2(3,L,N) +
     2                                DT(I,K) * (COEF2(4,L,N) +
     3                                DT(I,K) * COEF2(5,L,N))))
C
                X12 =  COEF1(1,LP1,M) + DT(I,K) * (COEF1(2,LP1,M) +
     1                                  DT(I,K) * (COEF1(3,LP1,M) +
     2                                  DT(I,K) * (COEF1(4,LP1,M) +
     3                                  DT(I,K) * COEF1(5,LP1,M))))
                X22 =  COEF1(1,LP1,N) + DT(I,K) * (COEF1(2,LP1,N) +
     1                                  DT(I,K) * (COEF1(3,LP1,N) +
     2                                  DT(I,K) * (COEF1(4,LP1,N) +
     3                                  DT(I,K) * COEF1(5,LP1,N))))
C
                Y12 =  COEF2(1,LP1,M) + DT(I,K) * (COEF2(2,LP1,M) +
     1                                  DT(I,K) * (COEF2(3,LP1,M) +
     2                                  DT(I,K) * (COEF2(4,LP1,M) +
     3                                  DT(I,K) * COEF2(5,LP1,M))))
                Y22 =  COEF2(1,LP1,N) + DT(I,K) * (COEF2(2,LP1,N) +
     1                                  DT(I,K) * (COEF2(3,LP1,N) +
     2                                  DT(I,K) * (COEF2(4,LP1,N) +
     3                                  DT(I,K) * COEF2(5,LP1,N))))
C
                X1  =  X11 + (X12 - X11) * DIR(I,K)
                X2  =  X21 + (X22 - X21) * DIR(I,K)
                Y1  =  Y11 + (Y12 - Y11) * DIR(I,K)
                Y2  =  Y21 + (Y22 - Y21) * DIR(I,K)
              ELSE
                X1  =  COEF1(1,5,M) + DT(I,K) * (COEF1(2,5,M) + 
     1                                DT(I,K) * (COEF1(3,5,M) + 
     2                                DT(I,K) * (COEF1(4,5,M) +
     3                                DT(I,K) * COEF1(5,5,M))))
                X2  =  COEF1(1,5,N) + DT(I,K) * (COEF1(2,5,N) + 
     1                                DT(I,K) * (COEF1(3,5,N) + 
     2                                DT(I,K) * (COEF1(4,5,N) +
     3                                DT(I,K) * COEF1(5,5,N))))
                Y1  =  COEF2(1,5,M) + DT(I,K) * (COEF2(2,5,M) + 
     1                                DT(I,K) * (COEF2(3,5,M) + 
     2                                DT(I,K) * (COEF2(4,5,M) +
     3                                DT(I,K) * COEF2(5,5,M))))
                Y2  =  COEF2(1,5,N) + DT(I,K) * (COEF2(2,5,N) +
     1                                DT(I,K) * (COEF2(3,5,N) + 
     2                                DT(I,K) * (COEF2(4,5,N) +
     3                                DT(I,K) * COEF2(5,5,N))))
              ENDIF
                TAUG(I,K) =  TAUG(I,K) + ((X1 - Y1 + 
     1                      (X2 - X1 - Y2 + Y1) * DIP(I,K)) * 1.608 * 
     2                       S(I,K) + Y1 + (Y2 - Y1) * DIP(I,K)) * 
     3                       S(I,K) * DP(I,K)
            ENDIF
  150     CONTINUE
        ENDIF
  200 CONTINUE
C
      RETURN
      END

      BLOCK DATA CKDLW
C
      IMPLICIT REAL (A-H,O-Z),
     +INTEGER (I-N)
C
C----------------------------------------------------------------------C
C     Li 2006 REVISED VERSION, COMPARED TO LI / BARKER (2005) THE MINORC
C     GAS, CO2, N2O, CH4 ECT ARE MORE ACCURATELY EVALUED               C
C----------------------------------------------------------------------C
C
C----------------------------------------------------------------------C
C     SPECTRAL REGION: 2200-2500 CM^-1. GAS: H2O, CO2                  C
C     K IN UNIT CM^2 / GRAM FOR ALL DATA IN LONGWAVE                   C
C----------------------------------------------------------------------C
C
      COMMON /BANDL1/ GWL1(1), CL1H2O(5,18), CL1CO2(5,18)
      DATA GWL1 / .4500 / 
      DATA (( CL1H2O(K,J), K = 1, 5), J = 1, 18)             /
     1 .7729E-02, .2411E-03, .2956E-05, .1077E-07,-.4460E-10,
     1 .7769E-02, .2493E-03, .3163E-05, .9180E-08,-.8711E-10,
     1 .8240E-02, .2499E-03, .2716E-05, .9216E-08,-.1298E-10,
     1 .8275E-02, .2642E-03, .3245E-05, .1214E-07,-.3853E-10,
     1 .8342E-02, .2629E-03, .3204E-05, .1246E-07,-.3213E-10,
     1 .8414E-02, .2588E-03, .3128E-05, .1643E-07, .2609E-10,
     1 .8506E-02, .2647E-03, .3275E-05, .1555E-07,-.2425E-11,
     1 .8624E-02, .2658E-03, .3256E-05, .1554E-07, .6619E-12,
     1 .8844E-02, .2684E-03, .3227E-05, .1590E-07, .1128E-10,
     1 .9063E-02, .2741E-03, .3251E-05, .1459E-07,-.3333E-11,
     1 .9189E-02, .2746E-03, .3223E-05, .1586E-07, .1758E-10,
     1 .9306E-02, .2778E-03, .3257E-05, .1565E-07, .1259E-10,
     1 .9324E-02, .2823E-03, .3399E-05, .1497E-07,-.1343E-10,
     1 .9387E-02, .2820E-03, .3363E-05, .1622E-07, .7870E-11,
     1 .9468E-02, .2822E-03, .3352E-05, .1724E-07, .2221E-10,
     1 .9603E-02, .2824E-03, .3283E-05, .1713E-07, .2963E-10,
     1 .9891E-02, .2820E-03, .2973E-05, .1730E-07, .8245E-10,
     1 .9909E-02, .3031E-03, .3554E-05, .1315E-07,-.3465E-10 /
      DATA (( CL1CO2(K,J), K = 1, 5), J = 1, 18)             /
     1 .1225E-02, .2765E-04, .2198E-06, .7150E-09, .1736E-11,
     1 .1702E-02, .3824E-04, .3057E-06, .9857E-09, .1965E-11,
     1 .2359E-02, .5288E-04, .4341E-06, .1416E-08, .1210E-11,
     1 .3290E-02, .7353E-04, .6138E-06, .2027E-08, .2927E-12,
     1 .4627E-02, .1024E-03, .8357E-06, .2952E-08, .4676E-11,
     1 .6529E-02, .1434E-03, .1148E-05, .4100E-08, .9037E-11,
     1 .9233E-02, .2015E-03, .1577E-05, .5542E-08, .1518E-10,
     1 .1303E-01, .2829E-03, .2162E-05, .7338E-08, .2307E-10,
     1 .1829E-01, .3952E-03, .2984E-05, .9599E-08, .2792E-10,
     1 .2544E-01, .5478E-03, .4138E-05, .1292E-07, .3108E-10,
     1 .3508E-01, .7593E-03, .5745E-05, .1694E-07, .3395E-10,
     1 .4732E-01, .1047E-02, .8186E-05, .2136E-07, .2658E-11,
     1 .6245E-01, .1406E-02, .1136E-04, .2811E-07,-.4140E-10,
     1 .8172E-01, .1836E-02, .1490E-04, .3860E-07,-.4997E-10,
     1 .1064E+00, .2378E-02, .1923E-04, .4933E-07,-.8392E-10,
     1 .1375E+00, .3059E-02, .2446E-04, .6294E-07,-.7413E-10,
     1 .1773E+00, .3944E-02, .3123E-04, .7490E-07,-.1274E-09,
     1 .2250E+00, .4978E-02, .3971E-04, .9725E-07,-.2245E-09 /
C
      COMMON /BANDL1GH/ GWL1GH(5), CL1CO2GH(5,28,5)
      DATA GWL1GH / .31000, .13000, .07000, .03000, .01000 /
      DATA ((( CL1CO2GH(K,J,I), K = 1, 5), J = 1, 28), I = 1, 3) /
     2 .6738E-02, .1474E-03, .1914E-05, .2326E-07, .1486E-09,
     2 .1066E-01, .2126E-03, .2442E-05, .2517E-07, .1354E-09,
     2 .1737E-01, .3190E-03, .3098E-05, .2805E-07, .1422E-09,
     2 .2858E-01, .4879E-03, .3877E-05, .3181E-07, .1829E-09,
     2 .4760E-01, .7644E-03, .4907E-05, .3653E-07, .2698E-09,
     2 .7982E-01, .1237E-02, .7031E-05, .3818E-07, .2891E-09,
     2 .1352E+00, .2034E-02, .1055E-04, .3694E-07, .2859E-09,
     2 .2326E+00, .3384E-02, .1584E-04, .3587E-07, .3916E-09,
     2 .4048E+00, .5695E-02, .2485E-04, .3519E-07, .5623E-09,
     2 .7080E+00, .9656E-02, .4384E-04, .3414E-07, .1870E-09,
     2 .1706E+01, .2188E-01, .1009E-03, .3305E-07,-.1328E-08,
     2 .2386E+01, .2974E-01, .1333E-03, .4796E-07,-.1721E-08,
     2 .3350E+01, .4070E-01, .1757E-03, .3596E-07,-.2213E-08,
     2 .4730E+01, .5602E-01, .2264E-03, .2641E-07,-.1543E-08,
     2 .6689E+01, .7768E-01, .3047E-03,-.3627E-07,-.2017E-08,
     2 .9492E+01, .1081E+00, .4080E-03,-.7802E-07,-.1117E-08,
     2 .1348E+02, .1506E+00, .5554E-03,-.1418E-06,-.4081E-09,
     2 .1910E+02, .2106E+00, .7977E-03,-.3192E-06,-.6841E-08,
     2 .2700E+02, .2918E+00, .1121E-02,-.3098E-06,-.1607E-07,
     2 .3803E+02, .4044E+00, .1567E-02,-.4187E-06,-.2876E-07,
     2 .5343E+02, .5610E+00, .2087E-02,-.3239E-06,-.2884E-07,
     2 .7462E+02, .7751E+00, .2679E-02,-.7249E-06,-.1690E-07,
     2 .1027E+03, .1064E+01, .3476E-02,-.2547E-05, .1010E-07,
     2 .1386E+03, .1409E+01, .4877E-02,-.3540E-05,-.2021E-07,
     2 .1857E+03, .1825E+01, .6552E-02,-.2374E-05,-.5292E-07,
     2 .2455E+03, .2303E+01, .8402E-02, .1602E-05,-.9114E-07,
     2 .3170E+03, .2824E+01, .9579E-02, .6663E-05,-.3177E-07,
     2 .3924E+03, .3305E+01, .1073E-01, .1084E-04,-.5578E-08,
     3 .2275E+00, .7262E-02, .1177E-03, .1392E-05, .8708E-08,
     3 .2788E+00, .8069E-02, .1224E-03, .1343E-05, .8082E-08,
     3 .3628E+00, .9253E-02, .1305E-03, .1318E-05, .7473E-08,
     3 .4993E+00, .1100E-01, .1428E-03, .1297E-05, .6672E-08,
     3 .7244E+00, .1350E-01, .1597E-03, .1314E-05, .6102E-08,
     3 .1101E+01, .1740E-01, .1916E-03, .1325E-05, .3991E-08,
     3 .1763E+01, .2312E-01, .2301E-03, .1408E-05, .2465E-08,
     3 .2935E+01, .3173E-01, .2773E-03, .1579E-05, .2010E-08,
     3 .5070E+01, .4504E-01, .3339E-03, .1827E-05, .2292E-08,
     3 .9005E+01, .6610E-01, .3775E-03, .2286E-05, .1025E-07,
     3 .2241E+02, .1307E+00, .4512E-03, .2587E-05, .3188E-07,
     3 .3170E+02, .1724E+00, .4976E-03, .2542E-05, .4311E-07,
     3 .4493E+02, .2297E+00, .7073E-03, .2388E-05, .3032E-07,
     3 .6390E+02, .3110E+00, .1110E-02, .1629E-05,-.1046E-07,
     3 .9126E+02, .4127E+00, .1572E-02, .2561E-05,-.5125E-07,
     3 .1308E+03, .5465E+00, .2060E-02, .4614E-05,-.9379E-07,
     3 .1880E+03, .7338E+00, .2436E-02, .8014E-05,-.1026E-06,
     3 .2711E+03, .9960E+00, .2260E-02, .1285E-04,-.6623E-08,
     3 .3909E+03, .1362E+01, .1590E-02, .1661E-04, .1523E-06,
     3 .5621E+03, .1903E+01, .6027E-04, .8748E-05, .4488E-06,
     3 .8007E+03, .2645E+01, .7700E-03,-.1335E-04, .4367E-06,
     3 .1131E+04, .3452E+01, .5129E-02,-.1868E-04,-.1365E-06,
     3 .1592E+04, .4300E+01, .7911E-02, .1566E-04,-.4882E-06,
     3 .2230E+04, .5343E+01, .3924E-02, .7204E-04,-.2967E-06,
     3 .3083E+04, .7011E+01,-.7174E-02, .9605E-04, .2091E-06,
     3 .4174E+04, .9686E+01,-.1919E-01, .1037E-03,-.6807E-07,
     3 .5522E+04, .1421E+02,-.4063E-01, .1461E-04, .2124E-05,
     3 .7210E+04, .2017E+02,-.3837E-01,-.2109E-03, .1876E-05,
     4 .1108E+02, .3429E+00, .4345E-02, .2875E-04, .9755E-07,
     4 .1144E+02, .3418E+00, .4266E-02, .2884E-04, .1077E-06,
     4 .1214E+02, .3422E+00, .4194E-02, .2895E-04, .1098E-06,
     4 .1333E+02, .3516E+00, .4189E-02, .2782E-04, .9596E-07,
     4 .1534E+02, .3690E+00, .4142E-02, .2645E-04, .9214E-07,
     4 .1851E+02, .3987E+00, .4196E-02, .2525E-04, .8338E-07,
     4 .2363E+02, .4459E+00, .4358E-02, .2358E-04, .6343E-07,
     4 .3217E+02, .5172E+00, .4574E-02, .2220E-04, .5251E-07,
     4 .4634E+02, .6332E+00, .5193E-02, .1942E-04,-.4200E-08,
     4 .7052E+02, .7966E+00, .6123E-02, .1998E-04,-.5578E-07,
     4 .1470E+03, .1198E+01, .7539E-02, .2687E-04,-.1083E-07,
     4 .1974E+03, .1421E+01, .8434E-02, .3125E-04,-.1393E-07,
     4 .2690E+03, .1711E+01, .8646E-02, .3705E-04, .1207E-06,
     4 .3707E+03, .2104E+01, .8274E-02, .3783E-04, .3138E-06,
     4 .5159E+03, .2658E+01, .8280E-02, .2972E-04, .3843E-06,
     4 .7219E+03, .3425E+01, .7638E-02, .1577E-04, .6599E-06,
     4 .1010E+04, .4531E+01, .9475E-02,-.2722E-04, .6672E-06,
     4 .1406E+04, .5896E+01, .1459E-01,-.8524E-04, .4200E-06,
     4 .1942E+04, .7336E+01, .2846E-01,-.1039E-03,-.8586E-06,
     4 .2677E+04, .8625E+01, .4491E-01,-.3324E-04,-.2342E-05,
     4 .3696E+04, .9635E+01, .4063E-01, .1259E-03,-.1060E-05,
     4 .5057E+04, .1136E+02, .2002E-01, .1965E-03, .1828E-05,
     4 .6836E+04, .1477E+02, .4018E-02, .9761E-04, .1983E-05,
     4 .9024E+04, .1972E+02,-.1033E-01,-.1716E-03, .3867E-05,
     4 .1163E+05, .2698E+02, .6287E-03,-.9831E-03, .1049E-04,
     4 .1477E+05, .2452E+02, .1734E+00,-.9441E-03,-.3499E-05,
     4 .1910E+05, .9065E+01, .1955E+00,-.3053E-03,-.2538E-05,
     4 .2461E+05,-.9810E+01, .1550E+00, .1568E-03,-.1033E-05 /
      DATA ((( CL1CO2GH(K,J,I), K = 1, 5), J = 1, 28), I = 4, 5) /
     5 .7986E+03, .1677E+02, .1210E+00, .2314E-03,-.1038E-05,
     5 .7993E+03, .1676E+02, .1211E+00, .2310E-03,-.1046E-05,
     5 .8010E+03, .1675E+02, .1212E+00, .2303E-03,-.1039E-05,
     5 .8054E+03, .1670E+02, .1211E+00, .2327E-03,-.9910E-06,
     5 .8153E+03, .1661E+02, .1210E+00, .2362E-03,-.9439E-06,
     5 .8377E+03, .1643E+02, .1206E+00, .2536E-03,-.1036E-05,
     5 .8863E+03, .1611E+02, .1183E+00, .3027E-03,-.1159E-05,
     5 .9827E+03, .1601E+02, .1064E+00, .3228E-03, .1490E-07,
     5 .1148E+04, .1637E+02, .9050E-01, .2998E-03, .1869E-05,
     5 .1416E+04, .1733E+02, .7984E-01, .2508E-03, .2858E-05,
     5 .2198E+04, .2018E+02, .6760E-01, .1663E-03, .5252E-05,
     5 .2681E+04, .2193E+02, .7743E-01, .7488E-04, .3599E-05,
     5 .3323E+04, .2426E+02, .1092E+00,-.5198E-04,-.1494E-05,
     5 .4212E+04, .2681E+02, .1294E+00,-.1206E-03,-.4404E-05,
     5 .5449E+04, .2865E+02, .1226E+00, .2527E-04,-.3612E-05,
     5 .7140E+04, .3023E+02, .6931E-01, .2631E-03, .5208E-05,
     5 .9356E+04, .3291E+02, .5977E-01, .3857E-03, .6337E-05,
     5 .1233E+05, .3880E+02, .8451E-01, .1655E-03,-.2356E-05,
     5 .1629E+05, .4636E+02, .6859E-02, .6641E-04, .6865E-05,
     5 .2128E+05, .5681E+02,-.4704E-01, .1862E-03, .9793E-05,
     5 .2726E+05, .7843E+02,-.4767E-01,-.1853E-02, .1714E-04,
     5 .3422E+05, .9059E+02, .2994E+00,-.3000E-02,-.9698E-05,
     5 .4300E+05, .6553E+02, .4254E+00,-.8920E-03,-.6760E-06,
     5 .5474E+05, .2816E+02, .4003E+00,-.5322E-03,-.4311E-06,
     5 .6905E+05,-.2155E+02, .3057E+00, .1376E-02,-.1023E-04,
     5 .8411E+05,-.5324E+02,-.1065E+00, .9477E-03, .1867E-04,
     5 .9612E+05,-.6344E+02,-.1826E+00, .6654E-04, .1295E-04,
     5 .1020E+06,-.6770E+02,-.1738E+00, .2740E-04, .7889E-05,
     6 .2751E+06,-.1179E+03,-.1559E+00, .1723E-02, .9732E-05,
     6 .2751E+06,-.1179E+03,-.1596E+00, .1720E-02, .1039E-04,
     6 .2751E+06,-.1179E+03,-.1582E+00, .1726E-02, .1009E-04,
     6 .2751E+06,-.1179E+03,-.1549E+00, .1730E-02, .9570E-05,
     6 .2750E+06,-.1178E+03,-.1565E+00, .1717E-02, .9918E-05,
     6 .2750E+06,-.1178E+03,-.1541E+00, .1702E-02, .9738E-05,
     6 .2750E+06,-.1177E+03,-.1550E+00, .1651E-02, .1055E-04,
     6 .2749E+06,-.1180E+03,-.1415E+00, .1612E-02, .9390E-05,
     6 .2747E+06,-.1189E+03,-.1235E+00, .1598E-02, .7860E-05,
     6 .2744E+06,-.1207E+03,-.1105E+00, .1556E-02, .7872E-05,
     6 .2737E+06,-.1259E+03,-.7812E-01, .1345E-02, .7062E-05,
     6 .2732E+06,-.1290E+03,-.7044E-01, .1259E-02, .8850E-05,
     6 .2725E+06,-.1327E+03,-.8298E-01, .1147E-02, .1445E-04,
     6 .2716E+06,-.1364E+03,-.8634E-01, .9498E-03, .1880E-04,
     6 .2701E+06,-.1385E+03,-.6870E-01, .5579E-03, .1955E-04,
     6 .2680E+06,-.1381E+03,-.2096E-01, .1477E-03, .1196E-04,
     6 .2650E+06,-.1353E+03,-.4792E-01, .4883E-05, .1189E-04,
     6 .2607E+06,-.1311E+03,-.1436E+00, .3124E-03, .2108E-04,
     6 .2547E+06,-.1254E+03,-.1750E+00, .4198E-03, .1541E-04,
     6 .2471E+06,-.1240E+03,-.2163E+00, .6966E-04, .1613E-04,
     6 .2378E+06,-.1394E+03,-.2353E+00, .1544E-02, .7290E-05,
     6 .2266E+06,-.1548E+03,-.5408E+00, .2372E-02, .2895E-04,
     6 .2121E+06,-.1434E+03,-.6396E+00, .4379E-03, .2133E-04,
     6 .1931E+06,-.1265E+03,-.5790E+00, .4702E-03, .1629E-04,
     6 .1698E+06,-.1066E+03,-.4847E+00, .3626E-03, .9090E-05,
     6 .1434E+06,-.8694E+02,-.4484E+00, .6072E-03, .1448E-04,
     6 .1162E+06,-.6672E+02,-.3481E+00, .3352E-03, .8064E-05,
     6 .9060E+05,-.5008E+02,-.2867E+00, .2329E-03, .1025E-04 /
C
C----------------------------------------------------------------------C
C     SPECTRAL REGION: 1900-2200 CM^-1. GAS: H2O                       C
C----------------------------------------------------------------------C
C
      COMMON /BANDL2/ GWL2(1), CL2H2O(5,18), CL2CS(5,3), CL2CF(5,3)
      DATA GWL2 / .65000 /
      DATA (( CL2H2O(K,J), K = 1, 5), J = 1, 18)             /
     1 .1635E-03, .2510E-05, .1165E-07, .6099E-11,-.6637E-13,
     1 .2428E-03, .3731E-05, .1727E-07, .8565E-11,-.8956E-13,
     1 .3603E-03, .5548E-05, .2564E-07, .1195E-10,-.1220E-12,
     1 .5345E-03, .8248E-05, .3813E-07, .1682E-10,-.1752E-12,
     1 .7928E-03, .1226E-04, .5699E-07, .2270E-10,-.2956E-12,
     1 .1175E-02, .1823E-04, .8542E-07, .3186E-10,-.5231E-12,
     1 .1739E-02, .2706E-04, .1287E-06, .4714E-10,-.9960E-12,
     1 .2570E-02, .4011E-04, .1939E-06, .7463E-10,-.1852E-11,
     1 .3793E-02, .5927E-04, .2896E-06, .1431E-09,-.2750E-11,
     1 .5588E-02, .8725E-04, .4266E-06, .2680E-09,-.3413E-11,
     1 .8207E-02, .1282E-03, .6235E-06, .4425E-09,-.3973E-11,
     1 .1202E-01, .1875E-03, .9017E-06, .7194E-09,-.3595E-11,
     1 .1748E-01, .2737E-03, .1297E-05, .8795E-09,-.2771E-11,
     1 .2519E-01, .3970E-03, .1885E-05, .9664E-09,-.3063E-11,
     1 .3575E-01, .5659E-03, .2843E-05, .1564E-08,-.2169E-10,
     1 .5020E-01, .7858E-03, .4096E-05, .3595E-08,-.4269E-10,
     1 .7024E-01, .1068E-02, .5442E-05, .7079E-08,-.4472E-10,
     1 .9820E-01, .1444E-02, .6762E-05, .9510E-08,-.2463E-10 /
      DATA (( CL2CS(K,J), K = 1, 5), J = 1, 3)               /
     1 .2074E+01,-.3404E-01, .3410E-03,-.2830E-05, .1315E-07,
     1 .3090E+01,-.5095E-01, .5008E-03,-.4236E-05, .2128E-07,
     1 .4599E+01,-.7532E-01, .7173E-03,-.6543E-05, .3852E-07 /
      DATA (( CL2CF(K,J), K = 1, 5), J = 1, 3)               /
     1 .1645E-01,-.1179E-04, .4388E-06,-.9480E-09,-.7306E-10,
     1 .2422E-01,-.1917E-04, .4174E-06,-.1452E-08,-.7633E-10,
     1 .3556E-01,-.1747E-04,-.4582E-07,-.4974E-08, .1378E-10 /
      COMMON /BANDL2GH/ GWL2GH(1), CL2H2OGH(5,28), 
     1                             CL2CSGH(5,3), CL2CFGH(5,3)
      DATA GWL2GH / .35000 /
      DATA (( CL2H2OGH(K,J), K = 1, 5), J = 1, 28)           /
     1 .3402E+01, .2934E-01, .4803E-05,-.2079E-06, .6822E-09,
     1 .3468E+01, .2991E-01, .4902E-05,-.2120E-06, .6954E-09,
     1 .3535E+01, .3049E-01, .5015E-05,-.2159E-06, .7054E-09,
     1 .3602E+01, .3106E-01, .5104E-05,-.2200E-06, .7192E-09,
     1 .3668E+01, .3164E-01, .5168E-05,-.2241E-06, .7385E-09,
     1 .3735E+01, .3221E-01, .5290E-05,-.2280E-06, .7472E-09,
     1 .3802E+01, .3279E-01, .5377E-05,-.2322E-06, .7613E-09,
     1 .3869E+01, .3336E-01, .5486E-05,-.2363E-06, .7722E-09,
     1 .3935E+01, .3394E-01, .5573E-05,-.2405E-06, .7872E-09,
     1 .4002E+01, .3451E-01, .5669E-05,-.2446E-06, .8003E-09,
     1 .4065E+01, .3506E-01, .5779E-05,-.2484E-06, .8106E-09,
     1 .4130E+01, .3564E-01, .5857E-05,-.2525E-06, .8281E-09,
     1 .4197E+01, .3621E-01, .5962E-05,-.2567E-06, .8412E-09,
     1 .4264E+01, .3679E-01, .6058E-05,-.2609E-06, .8562E-09,
     1 .4329E+01, .3736E-01, .6171E-05,-.2652E-06, .8692E-09,
     1 .4394E+01, .3794E-01, .6279E-05,-.2694E-06, .8841E-09,
     1 .4460E+01, .3850E-01, .6409E-05,-.2738E-06, .8969E-09,
     1 .4523E+01, .3907E-01, .6532E-05,-.2780E-06, .9115E-09,
     1 .4587E+01, .3963E-01, .6646E-05,-.2823E-06, .9258E-09,
     1 .4652E+01, .4019E-01, .6765E-05,-.2864E-06, .9355E-09,
     1 .4716E+01, .4075E-01, .6816E-05,-.2905E-06, .9500E-09,
     1 .4781E+01, .4130E-01, .6831E-05,-.2946E-06, .9629E-09,
     1 .4844E+01, .4183E-01, .6790E-05,-.2989E-06, .9745E-09,
     1 .4908E+01, .4236E-01, .6666E-05,-.3028E-06, .9879E-09,
     1 .4970E+01, .4287E-01, .6389E-05,-.3071E-06, .1009E-08,
     1 .5031E+01, .4335E-01, .6006E-05,-.3119E-06, .1030E-08,
     1 .5090E+01, .4381E-01, .5556E-05,-.3174E-06, .1046E-08,
     1 .5145E+01, .4424E-01, .5125E-05,-.3225E-06, .1050E-08 /
      DATA (( CL2CSGH(K,J), K = 1, 5), J = 1, 3)             /
     1 .6173E+01,-.9432E-01, .7384E-03,-.5151E-05, .3060E-07,
     1 .9216E+01,-.1404E+00, .1120E-02,-.7676E-05, .4232E-07,
     1 .1376E+02,-.2109E+00, .1733E-02,-.1102E-04, .4933E-07 /
      DATA (( CL2CFGH(K,J), K = 1, 5), J = 1, 3)             /
     1 .9702E-01,-.5003E-03, .1178E-05,-.7351E-08, .1902E-09,
     1 .1451E+00,-.7394E-03, .2253E-05,-.1085E-07, .2121E-09,
     1 .2170E+00,-.1121E-02, .4710E-05,-.1042E-07, .5378E-10 /
C
C----------------------------------------------------------------------C
C     SPECTRAL REGION: 1400-1900 CM^-1. GAS: H2O                       C
C----------------------------------------------------------------------C
C
      COMMON /BANDL3/ GWL3(2), CL3H2O(5,18,2), CL3CH4(2), CL3CS(5,4,2),
     1                CL3CF(5,4,2)
      DATA GWL3 / .35000, .42800 / 
      DATA CL3CH4 / 3.2742, 5.4269 /
      DATA ((( CL3H2O(K,J,I), K = 1, 5), J = 1, 18), I = 1, 2) /
     1 .1346E-01, .3885E-04,-.1051E-06, .1206E-09,-.3369E-11,
     1 .2004E-01, .5747E-04,-.1486E-06, .2288E-09,-.5499E-11,
     1 .2983E-01, .8477E-04,-.2198E-06, .4393E-09,-.7404E-11,
     1 .4442E-01, .1251E-03,-.3292E-06, .7923E-09,-.9568E-11,
     1 .6618E-01, .1843E-03,-.4857E-06, .1357E-08,-.1418E-10,
     1 .9847E-01, .2718E-03,-.7280E-06, .2189E-08,-.2010E-10,
     1 .1467E+00, .4012E-03,-.1094E-05, .3386E-08,-.2946E-10,
     1 .2183E+00, .5922E-03,-.1624E-05, .5178E-08,-.4821E-10,
     1 .3247E+00, .8738E-03,-.2425E-05, .7887E-08,-.7561E-10,
     1 .4824E+00, .1290E-02,-.3749E-05, .1189E-07,-.9333E-10,
     1 .7162E+00, .1903E-02,-.5877E-05, .1754E-07,-.8788E-10,
     1 .1060E+01, .2814E-02,-.9071E-05, .2461E-07,-.6135E-10,
     1 .1565E+01, .4155E-02,-.1449E-04, .3320E-07, .9926E-10,
     1 .2302E+01, .6147E-02,-.2199E-04, .3954E-07, .2928E-09,
     1 .3371E+01, .9014E-02,-.3310E-04, .5333E-07, .6451E-09,
     1 .4910E+01, .1318E-01,-.5060E-04, .6969E-07, .1479E-08,
     1 .7113E+01, .1895E-01,-.7459E-04, .1205E-06, .2632E-08,
     1 .1025E+02, .2642E-01,-.1077E-03, .2868E-06, .3422E-08,
     2 .1827E+00, .3217E-03, .2040E-06, .3824E-08,-.1900E-10,
     2 .2697E+00, .4479E-03, .1816E-06, .5756E-08,-.2757E-10,
     2 .3983E+00, .6269E-03, .6031E-07, .8292E-08,-.3313E-10,
     2 .5886E+00, .8792E-03,-.2372E-06, .1193E-07,-.3245E-10,
     2 .8702E+00, .1242E-02,-.7474E-06, .1636E-07,-.3505E-10,
     2 .1286E+01, .1761E-02,-.1512E-05, .2243E-07,-.5179E-10,
     2 .1902E+01, .2518E-02,-.2696E-05, .2919E-07,-.7042E-10,
     2 .2808E+01, .3618E-02,-.4325E-05, .3677E-07,-.1088E-09,
     2 .4137E+01, .5203E-02,-.7001E-05, .4220E-07,-.7366E-10,
     2 .6076E+01, .7428E-02,-.9254E-05, .5012E-07,-.1806E-09,
     2 .8884E+01, .1047E-01,-.7983E-05, .7354E-07,-.1034E-08,
     2 .1292E+02, .1458E-01,-.4125E-05, .1232E-06,-.2362E-08,
     2 .1871E+02, .1936E-01, .1851E-05, .3028E-06,-.4561E-08,
     2 .2696E+02, .2490E-01, .3493E-06, .6347E-06,-.7641E-08,
     2 .3871E+02, .3116E-01,-.6083E-04, .1167E-05,-.3467E-08,
     2 .5530E+02, .4078E-01,-.1803E-03, .1381E-05, .2540E-08,
     2 .7833E+02, .5606E-01,-.2409E-03, .1278E-05,-.6412E-08,
     2 .1096E+03, .6881E-01,-.3671E-03, .2262E-05,-.1106E-07 /
      DATA ((( CL3CS(K,J,I), K = 1, 5), J = 1, 4), I = 1, 2) /
     1 .3996E+02,-.5365E+00, .3540E-02,-.2717E-04, .2081E-06,
     1 .5942E+02,-.7948E+00, .5248E-02,-.4151E-04, .3219E-06,
     1 .8827E+02,-.1189E+01, .8359E-02,-.6131E-04, .3842E-06,
     1 .1308E+03,-.1790E+01, .1399E-01,-.8808E-04, .2861E-06,
     2 .7479E+02,-.1212E+01, .1079E-01,-.7085E-04, .2736E-06,
     2 .1124E+03,-.1825E+01, .1624E-01,-.1047E-03, .3862E-06,
     2 .1688E+03,-.2728E+01, .2371E-01,-.1581E-03, .6927E-06,
     2 .2536E+03,-.4072E+01, .3414E-01,-.2392E-03, .1289E-05 /
      DATA ((( CL3CF(K,J,I), K = 1, 5), J = 1, 4), I = 1, 2) /
     1 .1834E+01,-.3554E-02,-.1324E-05,-.1868E-06, .4608E-08,
     1 .2720E+01,-.5179E-02,-.1857E-06,-.3105E-06, .6846E-08,
     1 .4029E+01,-.7783E-02, .3094E-04,-.4834E-06, .4418E-08,
     1 .5946E+01,-.1209E-01, .1296E-03,-.7238E-06,-.9306E-08,
     2 .4059E+01,-.1863E-01, .8584E-04,-.2080E-06,-.2539E-08,
     2 .6099E+01,-.2803E-01, .1283E-03,-.2700E-06,-.4145E-08,
     2 .9158E+01,-.4174E-01, .1571E-03,-.4211E-06, .3909E-09,
     2 .1376E+02,-.6215E-01, .1582E-03,-.6142E-06, .1511E-07 /
      COMMON /BANDL3GH/ GWL3GH(3), CL3H2OGH(5,28,3), CL3CSGH(5,4),
     1                  CL3CFGH(5,4)
      DATA GWL3GH / .21200, .00900, .00100 /
      DATA ((( CL3H2OGH(K,J,I), K = 1, 5), J = 1, 28), I = 1, 3) /
     1 .2818E+01, .4977E-01, .3333E-03, .1033E-05, .5543E-09,
     1 .2852E+01, .5005E-01, .3357E-03, .1037E-05, .6017E-09,
     1 .2904E+01, .5029E-01, .3382E-03, .1042E-05, .6547E-09,
     1 .2990E+01, .5042E-01, .3410E-03, .1044E-05, .7345E-09,
     1 .3143E+01, .5032E-01, .3453E-03, .1065E-05, .4436E-09,
     1 .3436E+01, .4988E-01, .3482E-03, .1138E-05,-.1935E-09,
     1 .3997E+01, .4948E-01, .3348E-03, .1222E-05, .7964E-09,
     1 .5031E+01, .5016E-01, .3080E-03, .1223E-05, .2391E-08,
     1 .6866E+01, .5321E-01, .2905E-03, .1111E-05, .1466E-08,
     1 .1009E+02, .5979E-01, .2750E-03, .9415E-06, .1258E-08,
     1 .2023E+02, .7639E-01, .3233E-03, .1042E-05,-.1104E-07,
     1 .2702E+02, .8421E-01, .3657E-03, .1351E-05,-.1924E-07,
     1 .3660E+02, .9413E-01, .3198E-03, .1687E-05,-.1054E-07,
     1 .5018E+02, .1060E+00, .2284E-03, .2134E-05, .2182E-08,
     1 .6946E+02, .1204E+00, .1469E-03, .2435E-05, .4934E-08,
     1 .9688E+02, .1372E+00, .5108E-04, .2902E-05, .5369E-08,
     1 .1356E+03, .1614E+00,-.4585E-04, .3030E-05,-.1509E-08,
     1 .1902E+03, .1914E+00,-.2006E-03, .3070E-05,-.8193E-08,
     1 .2661E+03, .2191E+00,-.3400E-03, .4267E-05,-.3542E-07,
     1 .3698E+03, .2405E+00,-.7398E-03, .6909E-05,-.3960E-07,
     1 .5068E+03, .2605E+00,-.1477E-02, .1130E-04,-.4522E-07,
     1 .6786E+03, .3004E+00,-.2379E-02, .1174E-04,-.4238E-07,
     1 .8796E+03, .3413E+00,-.3577E-02, .9236E-05, .5275E-07,
     1 .1099E+04, .3193E+00,-.3336E-02, .7836E-05,-.5794E-07,
     1 .1322E+04, .2322E+00,-.3523E-02, .1497E-04,-.1227E-06,
     1 .1525E+04, .1338E+00,-.3548E-02, .1797E-04,-.1815E-06,
     1 .1696E+04,-.1624E-02,-.3572E-02, .1515E-04,-.1152E-06,
     1 .1825E+04,-.1337E+00,-.3511E-02, .5190E-05, .4501E-07,
     2 .1141E+05, .7342E+02,-.1039E+00,-.3141E-03, .3118E-05,
     2 .1148E+05, .7392E+02,-.1046E+00,-.3161E-03, .3141E-05,
     2 .1156E+05, .7441E+02,-.1052E+00,-.3181E-03, .3158E-05,
     2 .1165E+05, .7490E+02,-.1059E+00,-.3203E-03, .3172E-05,
     2 .1173E+05, .7538E+02,-.1065E+00,-.3228E-03, .3195E-05,
     2 .1181E+05, .7583E+02,-.1070E+00,-.3267E-03, .3215E-05,
     2 .1190E+05, .7621E+02,-.1068E+00,-.3309E-03, .3172E-05,
     2 .1202E+05, .7644E+02,-.1061E+00,-.3354E-03, .3129E-05,
     2 .1218E+05, .7635E+02,-.1052E+00,-.3401E-03, .3181E-05,
     2 .1243E+05, .7562E+02,-.1030E+00,-.3404E-03, .3128E-05,
     2 .1320E+05, .7181E+02,-.9546E-01,-.3504E-03, .2937E-05,
     2 .1379E+05, .6938E+02,-.8864E-01,-.3585E-03, .2327E-05,
     2 .1460E+05, .6611E+02,-.8538E-01,-.3954E-03, .2062E-05,
     2 .1570E+05, .6148E+02,-.7656E-01,-.3982E-03, .1060E-05,
     2 .1718E+05, .5490E+02,-.8533E-01,-.3193E-03, .3048E-05,
     2 .1908E+05, .4677E+02,-.8618E-01,-.2462E-03, .4298E-05,
     2 .2138E+05, .3694E+02,-.6726E-01,-.1685E-03, .2112E-05,
     2 .2403E+05, .2616E+02,-.6348E-01, .4946E-05,-.1312E-07,
     2 .2666E+05, .1605E+02,-.3698E-01, .1642E-03,-.5717E-05,
     2 .2888E+05, .6788E+01,-.4987E-01, .4297E-03,-.5067E-05,
     2 .3024E+05, .1303E+01,-.3833E-01, .1349E-03,-.4714E-05,
     2 .3042E+05,-.3760E+01,-.5905E-01,-.8417E-04, .5035E-05,
     2 .2933E+05,-.8881E+01, .1676E-01,-.5683E-04,-.3842E-05,
     2 .2728E+05,-.1272E+02, .1116E-01, .6502E-04,-.8432E-06,
     2 .2440E+05,-.1390E+02, .1179E-01,-.1524E-03, .1268E-05,
     2 .2112E+05,-.1461E+02, .2503E-01,-.3063E-03, .2521E-05,
     2 .1778E+05,-.1432E+02, .3541E-01,-.3318E-03, .1949E-05,
     2 .1459E+05,-.1398E+02, .5595E-01,-.1771E-03,-.2845E-05,
     3 .4620E+06,-.9277E+03, .8396E+00, .3048E-02,-.4573E-04,
     3 .4647E+06,-.9322E+03, .8442E+00, .3071E-02,-.4592E-04,
     3 .4674E+06,-.9377E+03, .8487E+00, .3089E-02,-.4613E-04,
     3 .4701E+06,-.9431E+03, .8526E+00, .3107E-02,-.4631E-04,
     3 .4727E+06,-.9485E+03, .8577E+00, .3119E-02,-.4673E-04,
     3 .4752E+06,-.9530E+03, .8603E+00, .3142E-02,-.4688E-04,
     3 .4776E+06,-.9576E+03, .8591E+00, .3162E-02,-.4652E-04,
     3 .4797E+06,-.9603E+03, .8565E+00, .3201E-02,-.4653E-04,
     3 .4811E+06,-.9603E+03, .8497E+00, .3272E-02,-.4681E-04,
     3 .4811E+06,-.9540E+03, .8312E+00, .3305E-02,-.4627E-04,
     3 .4742E+06,-.9195E+03, .7422E+00, .3348E-02,-.4106E-04,
     3 .4700E+06,-.8967E+03, .6618E+00, .3325E-02,-.3259E-04,
     3 .4628E+06,-.8660E+03, .6419E+00, .3593E-02,-.3209E-04,
     3 .4516E+06,-.8219E+03, .5752E+00, .3479E-02,-.2460E-04,
     3 .4350E+06,-.7585E+03, .6894E+00, .2544E-02,-.4573E-04,
     3 .4122E+06,-.6799E+03, .7230E+00, .1644E-02,-.5840E-04,
     3 .3823E+06,-.5859E+03, .5548E+00, .7847E-03,-.3405E-04,
     3 .3449E+06,-.4838E+03, .5570E+00,-.1030E-02,-.1073E-04,
     3 .3025E+06,-.3884E+03, .3217E+00,-.2968E-02, .5427E-04,
     3 .2574E+06,-.3008E+03, .5530E+00,-.6339E-02, .4849E-04,
     3 .2125E+06,-.2520E+03, .6118E+00,-.4391E-02, .4666E-04,
     3 .1707E+06,-.2130E+03, .1042E+01,-.2282E-02,-.5365E-04,
     3 .1341E+06,-.1744E+03, .5511E+00,-.2045E-02, .1543E-04,
     3 .1017E+06,-.1342E+03, .5526E+00,-.3116E-02, .1288E-04,
     3 .7507E+05,-.1065E+03, .6268E+00,-.2878E-02, .4302E-05,
     3 .5468E+05,-.8312E+02, .5621E+00,-.2137E-02, .2262E-05,
     3 .3924E+05,-.6453E+02, .5006E+00,-.1175E-02,-.4003E-05,
     3 .2792E+05,-.4693E+02, .3507E+00,-.9196E-03, .8430E-05 /
      DATA (( CL3CSGH(K,J), K = 1, 5), J = 1, 4)             /
     1 .8254E+02,-.1308E+01, .1136E-01,-.8548E-04, .4490E-06,
     1 .1227E+03,-.1942E+01, .1684E-01,-.1282E-03, .6873E-06,
     1 .1828E+03,-.2898E+01, .2532E-01,-.1890E-03, .9647E-06,
     1 .2721E+03,-.4314E+01, .3809E-01,-.2830E-03, .1388E-05 /
      DATA (( CL3CFGH(K,J), K = 1, 5), J = 1, 4)             /
     1 .4588E+01,-.1977E-01, .8047E-04,-.3912E-06, .1906E-08,
     1 .6826E+01,-.2936E-01, .1185E-03,-.6154E-06, .3345E-08,
     1 .1018E+02,-.4402E-01, .1899E-03,-.8496E-06, .2045E-08,
     1 .1517E+02,-.6521E-01, .3075E-03,-.1328E-05,-.7347E-09 /
C
C----------------------------------------------------------------------C
C     SPECTRAL REGION: 1100_1400 CM^-1. GAS: H2O, N2O, CH4             C
C----------------------------------------------------------------------C
C
      COMMON /BANDL4/ GWL4(5), CL4H2O(5,18,5), CL4N2O(5,18,5), 
     1                CL4CH4(5,18,5), CL4CS(5,4,5), CL4CF(5,4,5),
     2                CL4F12(5)
      DATA (CL4F12(I), I = 1, 5) / .2137E+04, .5064E+03, 0., 0.,
     1                             .7622E+02                   /
      DATA GWL4 / .41600, .19000, .05600, .06600, .13000 /
      DATA ((( CL4H2O(K,J,I), K = 1, 5), J = 1, 18), I = 1, 5) /
     1 .4522E-04, .8245E-06, .5346E-08, .1280E-10,-.1079E-13,
     1 .6617E-04, .1203E-05, .7622E-08, .1888E-10, .1704E-13,
     1 .9684E-04, .1760E-05, .1107E-07, .2745E-10, .3642E-13,
     1 .1415E-03, .2557E-05, .1568E-07, .4371E-10, .1647E-12,
     1 .2058E-03, .3743E-05, .2343E-07, .6158E-10, .1603E-12,
     1 .2993E-03, .5469E-05, .3457E-07, .8483E-10, .1397E-12,
     1 .4342E-03, .7950E-05, .5162E-07, .1268E-09, .2153E-13,
     1 .6306E-03, .1149E-04, .7380E-07, .1884E-09, .1725E-12,
     1 .9155E-03, .1675E-04, .1097E-06, .2639E-09,-.2506E-12,
     1 .1329E-02, .2417E-04, .1555E-06, .4129E-09, .3781E-12,
     1 .1916E-02, .3479E-04, .2190E-06, .6039E-09, .1572E-11,
     1 .2705E-02, .4936E-04, .3259E-06, .8815E-09, .1638E-12,
     1 .3820E-02, .6863E-04, .4362E-06, .1232E-08, .1997E-11,
     1 .5350E-02, .9586E-04, .6027E-06, .1772E-08, .4391E-11,
     1 .7351E-02, .1320E-03, .8311E-06, .2401E-08, .5945E-11,
     1 .9807E-02, .1753E-03, .1111E-05, .3083E-08, .3626E-11,
     1 .1287E-01, .2297E-03, .1492E-05, .4040E-08,-.2848E-11,
     1 .1727E-01, .3040E-03, .1918E-05, .5621E-08, .5334E-11,
     2 .1656E-02, .2445E-04, .1066E-06, .6257E-10, .5280E-12,
     2 .2420E-02, .3584E-04, .1581E-06, .8450E-10, .6278E-12,
     2 .3532E-02, .5255E-04, .2433E-06, .9942E-10,-.1181E-11,
     2 .5160E-02, .7654E-04, .3459E-06, .1875E-09, .4289E-12,
     2 .7508E-02, .1115E-03, .5156E-06, .3366E-09,-.2368E-13,
     2 .1092E-01, .1620E-03, .7701E-06, .5736E-09,-.2538E-11,
     2 .1585E-01, .2338E-03, .1129E-05, .1361E-08,-.8501E-12,
     2 .2304E-01, .3376E-03, .1616E-05, .2244E-08, .2423E-11,
     2 .3333E-01, .4885E-03, .2410E-05, .3030E-08,-.1164E-10,
     2 .4812E-01, .7065E-03, .3576E-05, .3943E-08,-.3673E-10,
     2 .6940E-01, .1013E-02, .5035E-05, .6085E-08,-.4138E-10,
     2 .9942E-01, .1434E-02, .6652E-05, .9787E-08, .2130E-10,
     2 .1407E+00, .2016E-02, .8991E-05, .1109E-07, .4046E-10,
     2 .1970E+00, .2802E-02, .1200E-04, .1927E-07, .1629E-09,
     2 .2707E+00, .3848E-02, .1668E-04, .2874E-07, .2089E-09,
     2 .3622E+00, .5233E-02, .2438E-04, .3566E-07,-.3670E-10,
     2 .4726E+00, .7117E-02, .3608E-04, .4033E-07,-.2733E-09,
     2 .6076E+00, .9701E-02, .5685E-04, .2606E-07,-.1326E-08,
     3 .2216E-03, .3554E-05, .2075E-07,-.2917E-10,-.9166E-12,
     3 .3235E-03, .5219E-05, .3392E-07,-.2936E-10,-.1760E-11,
     3 .4736E-03, .7562E-05, .4818E-07, .6518E-11,-.1707E-11,
     3 .6990E-03, .1127E-04, .7254E-07, .1496E-10,-.2450E-11,
     3 .1048E-02, .1664E-04, .9840E-07, .4709E-10,-.2066E-11,
     3 .1557E-02, .2410E-04, .1391E-06, .1806E-09,-.2184E-11,
     3 .2336E-02, .3598E-04, .1922E-06, .2000E-09,-.1765E-11,
     3 .3477E-02, .5344E-04, .2640E-06, .4470E-09, .3117E-11,
     3 .5099E-02, .7906E-04, .4082E-06, .3883E-09,-.8146E-12,
     3 .7646E-02, .1168E-03, .4470E-06,-.1348E-09, .1608E-10,
     3 .1126E-01, .1758E-03, .7926E-06,-.1417E-08,-.8058E-11,
     3 .1664E-01, .2664E-03, .1509E-05,-.2640E-08,-.7754E-10,
     3 .2434E-01, .3944E-03, .1962E-05, .7539E-11,-.8640E-11,
     3 .3374E-01, .5451E-03, .2840E-05, .3798E-08, .1227E-10,
     3 .4843E-01, .8062E-03, .4038E-05, .1406E-08, .3070E-10,
     3 .6891E-01, .1094E-02, .6680E-05, .8490E-08,-.1232E-09,
     3 .9395E-01, .1433E-02, .1118E-04, .2712E-07,-.4738E-09,
     3 .1249E+00, .1763E-02, .8218E-05, .5259E-07, .4610E-09,
     4 .2978E-03, .4096E-05, .1821E-07, .1786E-10,-.1332E-12,
     4 .4490E-03, .6184E-05, .2733E-07,-.4189E-11,-.5915E-12,
     4 .6686E-03, .9254E-05, .4312E-07,-.1507E-10,-.1331E-11,
     4 .9926E-03, .1376E-04, .7283E-07,-.3685E-10,-.3920E-11,
     4 .1487E-02, .2058E-04, .1053E-06,-.5939E-10,-.5269E-11,
     4 .2232E-02, .3027E-04, .1280E-06, .3122E-10,-.1858E-11,
     4 .3334E-02, .4576E-04, .2042E-06,-.1316E-09,-.6451E-11,
     4 .4918E-02, .6854E-04, .3483E-06,-.4051E-09,-.1835E-10,
     4 .7333E-02, .9814E-04, .3934E-06, .1784E-09,-.1910E-11,
     4 .1075E-01, .1408E-03, .5224E-06, .9762E-09, .1017E-10,
     4 .1545E-01, .2011E-03, .9206E-06, .2739E-08, .1004E-11,
     4 .2213E-01, .2862E-03, .1320E-05, .4216E-08, .9818E-11,
     4 .3112E-01, .4058E-03, .2414E-05, .4595E-08,-.9702E-10,
     4 .4392E-01, .5939E-03, .3704E-05,-.1661E-08,-.2390E-09,
     4 .6035E-01, .8365E-03, .6483E-05,-.4570E-08,-.5600E-09,
     4 .8357E-01, .1060E-02, .4470E-05, .8835E-08, .1228E-09,
     4 .1158E+00, .1378E-02, .3205E-05, .9850E-08, .4720E-09,
     4 .1582E+00, .1941E-02, .5674E-05, .9262E-08, .4614E-09,
     5 .3080E-01, .4314E-03, .1986E-05, .1956E-08,-.1402E-10,
     5 .4426E-01, .6104E-03, .2697E-05, .2671E-08,-.1020E-10,
     5 .6373E-01, .8681E-03, .3768E-05, .3091E-08,-.2370E-10,
     5 .9181E-01, .1233E-02, .5120E-05, .3968E-08,-.1630E-10,
     5 .1324E+00, .1759E-02, .7044E-05, .4812E-08,-.5648E-11,
     5 .1911E+00, .2512E-02, .9758E-05, .5417E-08, .1050E-11,
     5 .2757E+00, .3591E-02, .1365E-04, .7245E-08, .1285E-10,
     5 .3972E+00, .5141E-02, .1923E-04, .6153E-08,-.1002E-10,
     5 .5701E+00, .7312E-02, .2630E-04, .6808E-08, .8962E-10,
     5 .8135E+00, .1036E-01, .3567E-04,-.4509E-10, .2321E-09,
     5 .1147E+01, .1455E-01, .5332E-04,-.1359E-07,-.4653E-09,
     5 .1591E+01, .2004E-01, .7386E-04,-.1697E-07,-.8798E-09,
     5 .2162E+01, .2714E-01, .1005E-03,-.9313E-08,-.1059E-08,
     5 .2882E+01, .3633E-01, .1426E-03, .8892E-07,-.9844E-09,
     5 .3810E+01, .4930E-01, .2167E-03, .1527E-06,-.3229E-08,
     5 .5089E+01, .6618E-01, .2436E-03, .2950E-06, .4565E-08,
     5 .6791E+01, .9040E-01, .3260E-03, .6389E-08, .2512E-08,
     5 .8947E+01, .1217E+00, .4721E-03,-.2841E-06,-.4937E-08 /
      DATA ((( CL4N2O(K,J,I), K = 1, 5), J = 1, 18), I = 1, 5) /
     1 .1221E+00, .4607E-02,-.1092E-03,-.7950E-06, .2260E-07,
     1 .1736E+00, .4696E-02,-.1028E-03,-.7547E-06, .2142E-07,
     1 .2480E+00, .4554E-02,-.8742E-04,-.6458E-06, .1842E-07,
     1 .3534E+00, .2425E-03, .5366E-04, .3462E-06,-.1017E-07,
     1 .5066E+00, .6379E-03, .6300E-04, .4188E-06,-.1199E-07,
     1 .7217E+00, .7438E-03, .8979E-04, .6128E-06,-.1728E-07,
     1 .1022E+01,-.1727E-02, .1967E-03, .1406E-05,-.3894E-07,
     1 .1423E+01, .2931E-02, .1229E-03, .8570E-06,-.2317E-07,
     1 .1969E+01, .4369E-02, .1433E-03, .1043E-05,-.2714E-07,
     1 .2703E+01, .4148E-02, .2228E-03, .1704E-05,-.4389E-07,
     1 .3660E+01, .1966E-02, .3835E-03, .2971E-05,-.7781E-07,
     1 .4891E+01, .6197E-02, .3848E-03, .3143E-05,-.7857E-07,
     1 .6230E+01, .9593E-02, .5283E-03, .3610E-05,-.1012E-06,
     1 .7923E+01, .1151E-01, .6181E-03, .3950E-05,-.1194E-06,
     1 .9847E+01, .1396E-01, .7092E-03, .4052E-05,-.1370E-06,
     1 .1212E+02, .1507E-01, .7382E-03, .4072E-05,-.1381E-06,
     1 .1381E+02, .1632E-01, .7336E-03, .3528E-05,-.1242E-06,
     1 .1469E+02, .6227E-02, .5287E-03, .6379E-05,-.9610E-07,
     2 .1673E+00, .3458E-03, .1699E-04, .9082E-07,-.3284E-08,
     2 .2408E+00, .4617E-03, .2320E-04, .1292E-06,-.4544E-08,
     2 .3458E+00, .5739E-03, .3052E-04, .1966E-06,-.5963E-08,
     2 .4970E+00, .7537E-03, .4125E-04, .2809E-06,-.8042E-08,
     2 .7095E+00, .1138E-02, .5930E-04, .3617E-06,-.1148E-07,
     2 .1013E+01, .1578E-02, .8125E-04, .4838E-06,-.1564E-07,
     2 .1433E+01, .2548E-02, .1216E-03, .5861E-06,-.2299E-07,
     2 .2019E+01, .3330E-02, .1598E-03, .7950E-06,-.3003E-07,
     2 .2879E+01, .3812E-02, .2007E-03, .1309E-05,-.3749E-07,
     2 .4062E+01, .5861E-02, .2928E-03, .1624E-05,-.5574E-07,
     2 .5640E+01, .9247E-02, .4214E-03, .1720E-05,-.7880E-07,
     2 .7643E+01, .1320E-01, .5633E-03, .1869E-05,-.9837E-07,
     2 .1072E+02, .1005E-01, .5204E-03, .3153E-05,-.9999E-07,
     2 .1472E+02, .1435E-01, .6798E-03, .3103E-05,-.1311E-06,
     2 .1964E+02, .1496E-01, .7514E-03, .4260E-05,-.1423E-06,
     2 .2542E+02, .1663E-03, .4323E-03, .8283E-05,-.1120E-06,
     2 .3044E+02, .8567E-02, .8712E-03, .9115E-05,-.2504E-06,
     2 .3769E+02, .1006E-01, .1067E-02, .1187E-04,-.2966E-06,
     3 .1484E+02, .5334E-01, .2501E-02, .1325E-04,-.4186E-06,
     3 .1941E+02, .6501E-01, .3013E-02, .1519E-04,-.5088E-06,
     3 .2572E+02, .7385E-01, .3493E-02, .1886E-04,-.5901E-06,
     3 .3425E+02, .8483E-01, .4074E-02, .2342E-04,-.6783E-06,
     3 .4614E+02, .9777E-01, .4768E-02, .2844E-04,-.7989E-06,
     3 .6298E+02, .1052E+00, .5435E-02, .3644E-04,-.9412E-06,
     3 .8640E+02, .1317E+00, .6768E-02, .4267E-04,-.1233E-05,
     3 .1173E+03, .1642E+00, .8449E-02, .5211E-04,-.1579E-05,
     3 .1586E+03, .1946E+00, .1018E-01, .6401E-04,-.1946E-05,
     3 .2076E+03, .2239E+00, .1144E-01, .6956E-04,-.2126E-05,
     3 .2711E+03, .2072E+00, .1130E-01, .7784E-04,-.2170E-05,
     3 .3552E+03, .2081E+00, .1209E-01, .8748E-04,-.2514E-05,
     3 .4641E+03, .1634E+00, .1221E-01, .1153E-03,-.2857E-05,
     3 .5812E+03, .1662E+00, .1184E-01, .1027E-03,-.2862E-05,
     3 .6853E+03, .4071E+00, .1649E-01, .3255E-04,-.3037E-05,
     3 .8408E+03, .7898E-01, .2854E-02, .1555E-04,-.4177E-07,
     3 .1031E+04,-.1059E+00,-.1015E-01,-.1071E-03, .2761E-05,
     3 .1245E+04,-.5350E+00,-.2760E-01,-.1375E-03, .6120E-05,
     4 .1791E+01, .4671E-02, .2245E-03, .1154E-05,-.4032E-07,
     4 .2465E+01, .5317E-02, .2662E-03, .1574E-05,-.4716E-07,
     4 .3441E+01, .5277E-02, .2950E-03, .2225E-05,-.5287E-07,
     4 .4766E+01, .7379E-02, .3979E-03, .2820E-05,-.7061E-07,
     4 .6615E+01, .9570E-02, .5076E-03, .3359E-05,-.9344E-07,
     4 .9163E+01, .1539E-01, .7306E-03, .3510E-05,-.1351E-06,
     4 .1279E+02, .1617E-01, .8247E-03, .5140E-05,-.1464E-06,
     4 .1789E+02, .2130E-01, .1031E-02, .5581E-05,-.1820E-06,
     4 .2391E+02, .3796E-01, .1474E-02, .3418E-05,-.2163E-06,
     4 .3294E+02, .3936E-01, .1483E-02, .3084E-05,-.2004E-06,
     4 .4575E+02, .2773E-01, .1137E-02, .5439E-05,-.1215E-06,
     4 .6301E+02, .1832E-01, .7399E-03, .4006E-05,-.5753E-07,
     4 .8609E+02, .8497E-02, .1991E-03,-.1211E-06, .3695E-07,
     4 .1156E+03,-.2822E-01,-.8305E-03, .5836E-05, .1911E-06,
     4 .1653E+03,-.1781E+00,-.4480E-02, .3480E-04, .5037E-06,
     4 .2123E+03, .5571E-01, .1559E-02,-.1667E-04,-.4486E-06,
     4 .2704E+03, .7448E-01, .6277E-03,-.3309E-04, .2755E-06,
     4 .3374E+03,-.7487E-02,-.4234E-02,-.6698E-04, .1356E-05,
     5 .7360E+00, .2716E-02, .1361E-03, .9457E-06,-.2056E-07,
     5 .9957E+00, .3583E-02, .1763E-03, .1155E-05,-.2717E-07,
     5 .1358E+01, .4862E-02, .2315E-03, .1361E-05,-.3651E-07,
     5 .1864E+01, .6831E-02, .3127E-03, .1536E-05,-.5179E-07,
     5 .2562E+01, .8541E-02, .3961E-03, .2102E-05,-.6380E-07,
     5 .3592E+01, .1004E-01, .4876E-03, .2997E-05,-.7786E-07,
     5 .4915E+01, .1421E-01, .6639E-03, .3703E-05,-.1042E-06,
     5 .6845E+01, .1895E-01, .8707E-03, .4489E-05,-.1401E-06,
     5 .9531E+01, .2304E-01, .1086E-02, .5880E-05,-.1806E-06,
     5 .1338E+02, .2529E-01, .1300E-02, .8908E-05,-.2159E-06,
     5 .1859E+02, .2983E-01, .1642E-02, .1271E-04,-.2808E-06,
     5 .2555E+02, .3943E-01, .2143E-02, .1555E-04,-.3852E-06,
     5 .3279E+02, .6610E-01, .3127E-02, .1593E-04,-.5566E-06,
     5 .4081E+02, .8791E-01, .4129E-02, .2120E-04,-.7149E-06,
     5 .5086E+02, .1277E+00, .5650E-02, .2339E-04,-.9697E-06,
     5 .6022E+02, .1590E+00, .7100E-02, .3072E-04,-.1207E-05,
     5 .7874E+02, .8541E-01, .6203E-02, .6140E-04,-.1318E-05,
     5 .9991E+02, .1184E+00, .8159E-02, .7115E-04,-.1868E-05 /
      DATA ((( CL4CH4(K,J,I), K = 1, 5), J = 1, 18), I = 1, 5) /
     1 .9137E-01, .1088E-02, .4972E-05, .1499E-07,-.2950E-10,
     1 .1187E+00, .1335E-02, .5474E-05, .1353E-07,-.4323E-10,
     1 .1551E+00, .1652E-02, .5309E-05, .1271E-07, .1100E-09,
     1 .2047E+00, .2143E-02, .7594E-05,-.6530E-08,-.2620E-09,
     1 .2732E+00, .2688E-02, .7900E-05,-.4138E-08,-.1992E-09,
     1 .3623E+00, .3390E-02, .9013E-05, .6175E-08,-.4666E-10,
     1 .4805E+00, .4327E-02, .1129E-04, .2055E-07, .1086E-09,
     1 .6329E+00, .5551E-02, .1557E-04, .3098E-07,-.8012E-13,
     1 .8298E+00, .7139E-02, .2222E-04, .5144E-07,-.1101E-09,
     1 .1092E+01, .9247E-02, .3166E-04, .4297E-07,-.1151E-08,
     1 .1422E+01, .1168E-01, .3258E-04, .1151E-06, .4173E-09,
     1 .1831E+01, .1444E-01, .3322E-04, .2244E-06, .2416E-08,
     1 .2325E+01, .1759E-01, .4613E-04, .1798E-06,-.1103E-08,
     1 .2957E+01, .2229E-01, .6204E-04, .2746E-07,-.5260E-08,
     1 .3733E+01, .2761E-01, .6807E-04, .7011E-07,-.4948E-08,
     1 .4620E+01, .3617E-01, .1329E-03,-.2251E-06,-.1553E-07,
     1 .5887E+01, .4158E-01,-.3489E-04, .4706E-06, .1494E-07,
     1 .7667E+01, .5108E-01,-.8798E-04, .8921E-06, .1366E-07,
     2 .2207E+00, .2008E-02, .8345E-05,-.3194E-08,-.7102E-10,
     2 .2915E+00, .2405E-02, .8360E-05, .1547E-08, .1615E-09,
     2 .3892E+00, .2974E-02, .8936E-05,-.7955E-08, .3137E-09,
     2 .5225E+00, .3597E-02, .8928E-05, .1703E-07, .8602E-09,
     2 .7027E+00, .4644E-02, .1317E-04,-.1240E-07, .7214E-09,
     2 .9451E+00, .6081E-02, .2750E-04,-.5610E-07,-.9667E-09,
     2 .1287E+01, .7621E-02, .3806E-04,-.3813E-07,-.1639E-08,
     2 .1759E+01, .9770E-02, .5419E-04,-.5969E-07,-.3062E-08,
     2 .2394E+01, .1202E-01, .7968E-04, .4902E-07,-.4798E-08,
     2 .3271E+01, .1417E-01, .8826E-04, .3344E-06,-.2054E-08,
     2 .4481E+01, .1632E-01, .1096E-03, .7152E-06,-.1640E-08,
     2 .6107E+01, .1992E-01, .1268E-03, .5190E-06,-.4426E-08,
     2 .8397E+01, .2262E-01, .1162E-03, .8987E-06,-.1597E-08,
     2 .1147E+02, .2732E-01, .7073E-04, .1168E-05, .3837E-08,
     2 .1583E+02, .3975E-01,-.1534E-03,-.1913E-06, .3267E-07,
     2 .2112E+02, .5559E-01,-.2737E-03,-.1964E-05, .6435E-07,
     2 .2715E+02, .5229E-01, .3344E-03, .1562E-05,-.5051E-07,
     2 .3467E+02, .3740E-01, .2171E-03, .3982E-05,-.2886E-07,
     3 .2407E+00, .2039E-02, .1470E-04,-.2640E-07,-.8912E-09,
     3 .3188E+00, .2290E-02, .1830E-04, .8208E-08,-.1229E-08,
     3 .4211E+00, .2495E-02, .2941E-04, .8411E-07,-.2567E-08,
     3 .5768E+00, .2748E-02, .3366E-04, .1720E-06,-.2882E-08,
     3 .8098E+00, .3386E-02, .3183E-04, .2231E-06,-.2288E-08,
     3 .1146E+01, .4600E-02, .2915E-04, .2112E-06,-.1932E-08,
     3 .1632E+01, .6820E-02, .5652E-05, .6807E-07, .1466E-08,
     3 .2301E+01, .1010E-01,-.1987E-04,-.1066E-06, .7104E-08,
     3 .3197E+01, .1400E-01,-.1927E-04,-.2818E-06, .9577E-08,
     3 .4407E+01, .1813E-01,-.2739E-07,-.4008E-06, .1006E-07,
     3 .6107E+01, .2184E-01, .3399E-04,-.6081E-06, .7130E-08,
     3 .8543E+01, .2574E-01, .1007E-03,-.4743E-06, .1007E-08,
     3 .1159E+02, .3344E-01, .3544E-03,-.9878E-06,-.2314E-07,
     3 .1563E+02, .3590E-01, .6387E-03,-.5254E-06,-.6602E-07,
     3 .2259E+02, .4160E-01, .6703E-03, .5133E-06,-.8723E-07,
     3 .3362E+02, .5251E-01,-.5808E-03, .1067E-05, .1409E-06,
     3 .4756E+02, .7566E-01,-.4022E-03,-.3155E-05, .8415E-07,
     3 .6182E+02, .8859E-01, .3841E-03,-.6404E-05,-.2222E-07,
     4 .2988E+02, .3101E+00, .2051E-02, .4433E-05,-.9226E-09,
     4 .3625E+02, .3201E+00, .2046E-02, .6589E-05, .1640E-07,
     4 .4474E+02, .3416E+00, .2380E-02, .6462E-05,-.2999E-07,
     4 .5692E+02, .3657E+00, .2588E-02, .6339E-05,-.7359E-07,
     4 .7390E+02, .3828E+00, .2355E-02, .8384E-05,-.3824E-07,
     4 .9619E+02, .4112E+00, .2900E-02, .8277E-05,-.1406E-06,
     4 .1281E+03, .4211E+00, .1764E-02, .1049E-04, .1216E-07,
     4 .1721E+03, .4390E+00,-.1284E-03, .1275E-04, .3038E-06,
     4 .2294E+03, .4536E+00,-.3047E-03, .1296E-04, .2741E-06,
     4 .3065E+03, .4878E+00, .7446E-03, .1574E-04, .6018E-07,
     4 .4118E+03, .5523E+00, .7160E-03, .1813E-04, .4300E-07,
     4 .5493E+03, .6324E+00, .6411E-03, .1278E-04, .1444E-06,
     4 .7262E+03, .5727E+00, .9721E-03, .2628E-04, .4376E-07,
     4 .9323E+03, .4248E+00, .6676E-02, .6926E-04,-.1448E-05,
     4 .1174E+04, .8226E+00,-.1731E-01, .7288E-05, .2157E-05,
     4 .1350E+04, .1466E+01, .1950E-02,-.6319E-04,-.1249E-05,
     4 .1469E+04, .1092E+01, .7563E-03, .4991E-04,-.3726E-06,
     4 .1518E+04, .6355E-01,-.1367E-02, .1037E-03, .1104E-05,
     5 .2036E+01, .1970E-01, .6954E-04, .9777E-06, .1427E-07,
     5 .2489E+01, .2249E-01, .9201E-04, .6975E-06, .9741E-08,
     5 .3087E+01, .2558E-01, .1240E-03, .5592E-06, .4320E-08,
     5 .3903E+01, .2815E-01, .1190E-03, .6690E-06, .7383E-08,
     5 .5022E+01, .3039E-01, .1048E-03, .1061E-05, .1199E-07,
     5 .6720E+01, .3072E-01,-.1784E-04, .1993E-05, .3744E-07,
     5 .8736E+01, .3483E-01, .1329E-03, .2483E-05, .1622E-07,
     5 .1156E+02, .4173E-01, .2808E-03, .2431E-05,-.7176E-08,
     5 .1559E+02, .5133E-01, .5124E-03, .1729E-05,-.5871E-07,
     5 .2104E+02, .6360E-01, .3861E-03, .1025E-05,-.4281E-07,
     5 .2888E+02, .7407E-01,-.2885E-03, .8115E-06, .4623E-07,
     5 .3933E+02, .1041E+00,-.5496E-03, .1333E-05, .2076E-07,
     5 .5151E+02, .1307E+00,-.1016E-02, .3993E-05, .8703E-07,
     5 .6489E+02, .1424E+00,-.1220E-02, .8490E-05, .2009E-06,
     5 .7902E+02, .4233E-01,-.2834E-04, .2733E-04, .5109E-07,
     5 .9987E+02, .5085E-02,-.1784E-02, .2788E-04, .1585E-06,
     5 .1163E+03, .1377E+00,-.5601E-03, .5898E-06, .6894E-07,
     5 .1301E+03, .2690E+00, .7167E-04,-.1505E-04,-.1757E-06 /
      DATA ((( CL4CS(K,J,I), K = 1, 5), J = 1, 4), I = 1, 5) /
     1 .2963E+01,-.6248E-01, .7072E-03,-.6852E-05, .4073E-07,
     1 .4389E+01,-.9282E-01, .1053E-02,-.1018E-04, .6033E-07,
     1 .6504E+01,-.1379E+00, .1563E-02,-.1513E-04, .8989E-07,
     1 .9655E+01,-.2045E+00, .2315E-02,-.2248E-04, .1345E-06,
     2 .5426E+01,-.7237E-01, .5760E-03,-.4182E-05, .2285E-07,
     2 .7880E+01,-.1059E+00, .8600E-03,-.5971E-05, .2769E-07,
     2 .1137E+02,-.1502E+00, .1194E-02,-.9983E-05, .6469E-07,
     2 .1633E+02,-.2158E+00, .1902E-02,-.1555E-04, .7317E-07,
     3 .1842E+01,-.3104E-01, .2853E-03,-.2193E-05, .1164E-07,
     3 .2775E+01,-.4614E-01, .4189E-03,-.3011E-05, .1430E-07,
     3 .4184E+01,-.6835E-01, .6010E-03,-.4193E-05, .1973E-07,
     3 .6221E+01,-.1013E+00, .8702E-03,-.6274E-05, .3281E-07,
     4 .2481E+01,-.3523E-01, .3003E-03,-.1935E-05, .6590E-08,
     4 .3682E+01,-.5227E-01, .4152E-03,-.2940E-05, .1664E-07,
     4 .5486E+01,-.7810E-01, .6293E-03,-.4376E-05, .2069E-07,
     4 .8208E+01,-.1172E+00, .9091E-03,-.6048E-05, .3264E-07,
     5 .8893E+01,-.1096E+00, .7672E-03,-.4571E-05, .1791E-07,
     5 .1335E+02,-.1655E+00, .1103E-02,-.6265E-05, .3180E-07,
     5 .2001E+02,-.2519E+00, .1738E-02,-.8552E-05, .2439E-07,
     5 .3000E+02,-.3746E+00, .2320E-02,-.1363E-04, .9085E-07 /
      DATA ((( CL4CF(K,J,I), K = 1, 5), J = 1, 4), I = 1, 5) /
     1 .7013E-03,-.1603E-05, .1962E-07, .3683E-10,-.3567E-11,
     1 .8985E-03,-.1971E-05, .2133E-07,-.5080E-10,-.2559E-11,
     1 .1188E-02,-.3713E-05, .2277E-07,-.6132E-11,-.1074E-11,
     1 .1773E-02,-.7232E-05, .3104E-08, .2434E-09, .2715E-11,
     2 .9373E-01,-.2093E-03, .3432E-06, .4281E-08, .4751E-10,
     2 .1318E+00,-.2736E-03, .7616E-06, .1319E-07,-.1956E-09,
     2 .1819E+00,-.1819E-03,-.1961E-05,-.1174E-07, .5800E-09,
     2 .2472E+00,-.1279E-04, .3770E-05,-.5185E-07,-.2937E-09,
     3 .6826E-02,-.2275E-04, .8054E-07,-.6566E-09, .1963E-10,
     3 .1086E-01,-.4792E-04, .4297E-06, .6485E-09,-.3405E-10,
     3 .1737E-01,-.8208E-04, .5458E-06, .1082E-08,-.3918E-10,
     3 .2609E-01,-.1179E-03, .1103E-06,-.6427E-09, .5693E-10,
     4 .3075E-01,-.9096E-04, .1324E-05,-.4341E-08,-.1020E-09,
     4 .4570E-01,-.1504E-03, .3913E-06,-.6475E-08, .1799E-09,
     4 .6792E-01,-.2302E-03, .1675E-05,-.1538E-07,-.8774E-12,
     4 .1032E+00,-.3967E-03, .1016E-05,-.6557E-08, .2296E-09,
     5 .2480E+00,-.1007E-02, .4769E-05,-.2681E-07,-.2797E-10,
     5 .3749E+00,-.1603E-02, .4728E-05,-.1799E-07, .3663E-09,
     5 .5657E+00,-.2589E-02, .1026E-04,-.5741E-08,-.3444E-09,
     5 .8529E+00,-.3892E-02, .3271E-05,-.3967E-07, .1784E-08 /
      COMMON /BANDL4GH/ GWL4GH(5), CL4H2OGH(5,28,5), CL4N2OGH(5,28,5),
     1                             CL4CH4GH(5,28,5)
      DATA GWL4GH / .06000, .02000, .04000, .02000, .00200 /
      DATA ((( CL4H2OGH(K,J,I), K = 1, 5), J = 1, 28), I = 1, 3) /
     1 .2389E-01, .8023E-03, .1119E-04, .7886E-07, .2424E-09,
     1 .2524E-01, .8180E-03, .1115E-04, .7834E-07, .2450E-09,
     1 .2783E-01, .8535E-03, .1113E-04, .7728E-07, .2525E-09,
     1 .3203E-01, .9216E-03, .1138E-04, .7569E-07, .2400E-09,
     1 .3887E-01, .1039E-02, .1195E-04, .7361E-07, .2101E-09,
     1 .4984E-01, .1225E-02, .1294E-04, .7330E-07, .1907E-09,
     1 .6746E-01, .1517E-02, .1453E-04, .7410E-07, .1717E-09,
     1 .9657E-01, .1984E-02, .1733E-04, .7454E-07, .7852E-10,
     1 .9675E-01, .1987E-02, .1709E-04, .7457E-07, .1291E-09,
     1 .2319E+00, .3958E-02, .2654E-04, .8626E-07, .1313E-09,
     1 .5036E+00, .7653E-02, .4228E-04, .8702E-07,-.2182E-10,
     1 .6835E+00, .1001E-01, .5237E-04, .8422E-07,-.2534E-09,
     1 .9369E+00, .1319E-01, .6201E-04, .8985E-07, .1192E-09,
     1 .1291E+01, .1752E-01, .7501E-04, .7891E-07, .2683E-09,
     1 .1788E+01, .2328E-01, .9027E-04, .1094E-06, .9234E-09,
     1 .2472E+01, .3122E-01, .1211E-03, .1212E-06,-.3275E-09,
     1 .3430E+01, .4227E-01, .1553E-03, .7843E-07,-.1063E-08,
     1 .4744E+01, .5677E-01, .1867E-03, .5932E-07,-.1628E-09,
     1 .6491E+01, .7539E-01, .2152E-03, .1557E-06, .3990E-08,
     1 .8763E+01, .1009E+00, .3017E-03, .2280E-06, .2903E-08,
     1 .1172E+02, .1371E+00, .4766E-03, .3806E-06,-.2891E-08,
     1 .1582E+02, .1884E+00, .5719E-03, .2031E-06, .9369E-08,
     1 .2147E+02, .2628E+00, .7701E-03,-.8997E-06, .4368E-08,
     1 .2891E+02, .3565E+00, .1105E-02,-.1426E-05,-.1688E-07,
     1 .3862E+02, .4676E+00, .1176E-02,-.1233E-05, .1237E-07,
     1 .5095E+02, .6090E+00, .1076E-02,-.2648E-05, .7162E-07,
     1 .6524E+02, .7809E+00, .1555E-02,-.5509E-05, .2720E-07,
     1 .8082E+02, .9810E+00, .2414E-02,-.8815E-05,-.5367E-07,
     2 .7236E+03, .7962E+01, .1904E-01,-.3996E-04,-.5148E-06,
     2 .7236E+03, .7962E+01, .1896E-01,-.4060E-04,-.4972E-06,
     2 .7236E+03, .7962E+01, .1896E-01,-.4067E-04,-.4960E-06,
     2 .7236E+03, .7962E+01, .1910E-01,-.4064E-04,-.5201E-06,
     2 .7236E+03, .7962E+01, .1936E-01,-.4068E-04,-.5609E-06,
     2 .7236E+03, .7968E+01, .1931E-01,-.4193E-04,-.5658E-06,
     2 .7236E+03, .7956E+01, .1883E-01,-.3833E-04,-.4978E-06,
     2 .7236E+03, .7938E+01, .1889E-01,-.3021E-04,-.4364E-06,
     2 .7242E+03, .7932E+01, .1871E-01,-.2850E-04,-.4095E-06,
     2 .7242E+03, .7896E+01, .2018E-01,-.4137E-04,-.4395E-06,
     2 .7350E+03, .7950E+01, .1039E-01,-.6774E-04, .1403E-05,
     2 .7356E+03, .7944E+01, .1055E-01,-.6654E-04, .1361E-05,
     2 .7434E+03, .7470E+01, .1583E-01, .3154E-04,-.4696E-07,
     2 .7446E+03, .7338E+01, .3541E-01, .6012E-04,-.4080E-05,
     2 .7464E+03, .7386E+01, .3458E-01, .5217E-04,-.3962E-05,
     2 .7488E+03, .7416E+01, .3224E-01, .4928E-04,-.3556E-05,
     2 .8112E+03, .7692E+01,-.2378E-01, .3646E-04, .5880E-05,
     2 .8088E+03, .7818E+01,-.1949E-01, .1360E-04, .5100E-05,
     2 .8052E+03, .9090E+01, .1899E-01,-.2113E-03,-.2146E-05,
     2 .7992E+03, .8844E+01, .1637E-01,-.5873E-04,-.6390E-08,
     2 .7920E+03, .8742E+01, .1591E-01,-.5873E-04, .1426E-07,
     2 .7818E+03, .8598E+01, .1531E-01,-.5706E-04, .4255E-07,
     2 .7674E+03, .8406E+01, .1509E-01,-.5579E-04,-.3805E-07,
     2 .7488E+03, .8136E+01, .1412E-01,-.4969E-04,-.4195E-07,
     2 .7260E+03, .7854E+01, .1159E-01,-.5386E-04, .2095E-06,
     2 .6942E+03, .7494E+01, .1076E-01,-.5536E-04, .1964E-06,
     2 .6540E+03, .7014E+01, .9864E-02,-.5251E-04, .1809E-06,
     2 .6072E+03, .6516E+01, .1027E-01,-.6330E-04,-.1893E-06,
     3 .2108E+01, .3332E-01, .1275E-03, .1178E-05, .1325E-07,
     3 .2108E+01, .3294E-01, .1379E-03, .1257E-05, .1107E-07,
     3 .2108E+01, .3294E-01, .1381E-03, .1268E-05, .1090E-07,
     3 .2126E+01, .3294E-01, .1198E-03, .1265E-05, .1391E-07,
     3 .2156E+01, .3293E-01, .8849E-04, .1269E-05, .1899E-07,
     3 .2117E+01, .3215E-01, .9393E-04, .1425E-05, .1947E-07,
     3 .2106E+01, .3404E-01, .1532E-03, .9540E-06, .1125E-07,
     3 .2108E+01, .3590E-01, .1476E-03,-.5847E-07, .3361E-08,
     3 .2075E+01, .3686E-01, .1700E-03,-.2662E-06, .3063E-10,
     3 .2045E+01, .3572E-01, .1901E-03,-.3326E-07,-.3711E-08,
     3 .1824E+01, .3158E-01, .2442E-03, .4196E-06,-.1189E-07,
     3 .1820E+01, .3198E-01, .2232E-03, .3441E-06,-.7679E-08,
     3 .1878E+01, .3197E-01, .1528E-03, .1701E-06, .2109E-08,
     3 .1815E+01, .3070E-01, .1596E-03, .2971E-06, .1288E-08,
     3 .1808E+01, .2829E-01, .8270E-04, .6339E-06, .1470E-07,
     3 .1766E+01, .2901E-01, .1250E-03, .2424E-06, .3840E-08,
     3 .1745E+01, .3131E-01, .1697E-03,-.9534E-06,-.1533E-07,
     3 .1680E+01, .2975E-01, .1597E-03,-.4634E-06,-.9044E-08,
     3 .1628E+01, .2907E-01, .1676E-03,-.9052E-07,-.6416E-08,
     3 .1541E+01, .2736E-01, .1364E-03, .2546E-06, .2139E-08,
     3 .1424E+01, .2636E-01, .1318E-03, .2340E-06, .1267E-08,
     3 .1267E+01, .2634E-01, .1556E-03, .8799E-07,-.3775E-08,
     3 .1109E+01, .2486E-01, .1188E-03, .4647E-07, .4857E-08,
     3 .8613E+00, .2496E-01, .1294E-03,-.3978E-06, .7470E-08,
     3 .4707E+00, .1772E-01, .3293E-03, .4122E-06,-.2613E-07,
     3 .4857E+00, .1040E-01, .1019E-03, .9338E-06, .6153E-08,
     3 .4899E+00, .1038E-01, .9413E-04, .2505E-06,-.1454E-08,
     3 .4833E+00, .7562E-02, .3546E-04, .3612E-06, .4517E-08 /
      DATA ((( CL4H2OGH(K,J,I), K = 1, 5), J = 1, 28), I = 4, 5) /
     4 .2355E+02, .2819E+00,-.4206E-03,-.6588E-05, .1073E-06,
     4 .2355E+02, .2819E+00,-.4209E-03,-.6590E-05, .1074E-06,
     4 .2355E+02, .2819E+00,-.4208E-03,-.6590E-05, .1074E-06,
     4 .2353E+02, .2819E+00,-.4211E-03,-.6590E-05, .1074E-06,
     4 .2353E+02, .2819E+00,-.4211E-03,-.6590E-05, .1074E-06,
     4 .2353E+02, .2819E+00,-.4217E-03,-.6590E-05, .1075E-06,
     4 .2353E+02, .2819E+00,-.4215E-03,-.6551E-05, .1070E-06,
     4 .2352E+02, .2817E+00,-.4235E-03,-.6551E-05, .1074E-06,
     4 .2352E+02, .2817E+00,-.4236E-03,-.6562E-05, .1073E-06,
     4 .2346E+02, .2916E+00,-.8387E-03,-.3810E-05, .1225E-06,
     4 .2105E+02, .2838E+00, .1490E-02, .1889E-05,-.3219E-06,
     4 .2091E+02, .2828E+00, .1485E-02, .1736E-05,-.3197E-06,
     4 .1859E+02, .3990E+00, .2987E-03,-.2243E-04, .1270E-07,
     4 .1823E+02, .4320E+00,-.4620E-02,-.2981E-04, .1022E-05,
     4 .1755E+02, .4209E+00,-.4274E-02,-.2850E-04, .9660E-06,
     4 .1670E+02, .4073E+00,-.3795E-02,-.2700E-04, .8865E-06,
     4 .5743E+00, .3263E+00, .1010E-01,-.2136E-04,-.1434E-05,
     4 .4236E+00, .2876E+00, .9015E-02,-.1662E-04,-.1252E-05,
     4 .3987E+00,-.4227E-01,-.6518E-03, .3882E-04, .5528E-06,
     4 .4442E+00, .5766E-02, .7944E-05,-.6069E-07, .7719E-09,
     4 .4914E+00, .7053E-02, .1367E-04,-.7704E-07, .1159E-08,
     4 .5693E+00, .9095E-02, .3801E-04,-.1469E-06,-.2042E-08,
     4 .6678E+00, .1032E-01, .2110E-04, .1560E-06, .4224E-08,
     4 .7842E+00, .1593E-01, .6816E-05,-.3554E-06, .1052E-07,
     4 .8238E+00, .2334E-01, .8207E-04,-.1137E-05, .4041E-08,
     4 .3765E+00, .2898E-01, .7379E-03,-.1364E-05,-.9786E-07,
     4 .3693E+00, .2225E-01, .5910E-03, .1003E-05,-.5321E-07,
     4 .4338E+00,-.1136E-04,-.6776E-04, .5445E-05, .7788E-07,
     5 .9687E-02, .1412E-03, .5307E-06,-.5223E-09,-.2049E-11,
     5 .9690E-02, .1413E-03, .5309E-06,-.5186E-09,-.2009E-11,
     5 .9698E-02, .1413E-03, .5304E-06,-.5111E-09,-.1814E-11,
     5 .9711E-02, .1415E-03, .5309E-06,-.4959E-09,-.1613E-11,
     5 .9738E-02, .1418E-03, .5316E-06,-.4677E-09,-.1178E-11,
     5 .9793E-02, .1425E-03, .5325E-06,-.4074E-09,-.2019E-12,
     5 .9906E-02, .1439E-03, .5349E-06,-.2852E-09, .1674E-11,
     5 .1013E-01, .1466E-03, .5391E-06,-.3668E-10, .5639E-11,
     5 .1013E-01, .1466E-03, .5391E-06,-.3675E-10, .5640E-11,
     5 .1153E-01, .1642E-03, .5758E-06, .9417E-09, .2184E-10,
     5 .1511E-01, .2118E-03, .6582E-06, .6767E-09, .3303E-10,
     5 .1770E-01, .2436E-03, .7868E-06, .1764E-08, .3789E-10,
     5 .2153E-01, .2931E-03, .9288E-06, .2837E-08, .5458E-10,
     5 .2721E-01, .3674E-03, .1228E-05, .5044E-08, .4404E-10,
     5 .3557E-01, .4854E-03, .1211E-05, .5428E-08, .1463E-09,
     5 .4920E-01, .7096E-03,-.7649E-06,-.1568E-07, .4641E-09,
     5 .6774E-01, .9642E-03,-.1386E-05,-.1477E-07, .6924E-09,
     5 .8295E-01, .1503E-02, .9383E-05,-.5270E-07,-.9002E-09,
     5 .1053E+00, .2090E-02, .2388E-04,-.6977E-07,-.2912E-08,
     5 .1493E+00, .2139E-02, .4439E-05, .4974E-07, .1430E-08,
     5 .1845E+00, .3423E-02, .5772E-05,-.8297E-07, .1452E-08,
     5 .1401E+00, .2064E-02,-.1348E-04, .1100E-06, .6345E-08,
     5 .1269E+00, .1383E-02, .6630E-05, .3254E-07, .2028E-09,
     5 .1572E+00, .2025E-02, .8321E-05,-.4928E-07,-.4596E-09,
     5 .1524E+00, .2465E-02, .3357E-04, .1213E-07,-.3288E-08,
     5 .1567E+00, .2385E-02, .8036E-05, .1980E-07, .7383E-09,
     5 .1728E+00, .2756E-02, .2018E-04,-.3805E-08,-.8982E-09,
     5 .1614E+00, .1926E-02, .1683E-04, .9033E-07, .2361E-09 /
      DATA ((( CL4N2OGH(K,J,I), K = 1, 5), J = 1, 28), I = 1, 3) /
     1 .4988E-01, .8580E-03, .3692E-04, .2180E-06,-.3876E-08,
     1 .5165E-01, .8462E-03, .3643E-04, .2110E-06,-.3960E-08,
     1 .5233E-01, .9050E-03, .3818E-04, .2030E-06,-.4235E-08,
     1 .5855E-01, .9452E-03, .4001E-04, .2128E-06,-.4518E-08,
     1 .7071E-01, .8644E-03, .3850E-04, .2357E-06,-.4512E-08,
     1 .8943E-01, .9599E-03, .4286E-04, .2450E-06,-.5571E-08,
     1 .1124E+00, .1075E-02, .4820E-04, .2629E-06,-.6742E-08,
     1 .1522E+00, .1398E-02, .6125E-04, .2914E-06,-.9041E-08,
     1 .2362E+00, .1645E-02, .7551E-04, .3925E-06,-.1201E-07,
     1 .3505E+00, .1991E-02, .9643E-04, .6076E-06,-.1485E-07,
     1 .6880E+00, .3824E-02, .1786E-03, .1011E-05,-.2760E-07,
     1 .9256E+00, .4666E-02, .2183E-03, .1189E-05,-.3527E-07,
     1 .1217E+01, .6512E-02, .2931E-03, .1344E-05,-.4883E-07,
     1 .1636E+01, .6174E-02, .3114E-03, .2063E-05,-.5101E-07,
     1 .2280E+01, .6586E-02, .3680E-03, .2917E-05,-.6316E-07,
     1 .3157E+01, .7771E-02, .4553E-03, .3738E-05,-.8330E-07,
     1 .4350E+01, .8261E-02, .5513E-03, .5395E-05,-.1028E-06,
     1 .6032E+01, .1205E-01, .7615E-03, .6747E-05,-.1472E-06,
     1 .8129E+01, .1545E-01, .9619E-03, .8408E-05,-.1837E-06,
     1 .1056E+02, .1961E-01, .1206E-02, .1050E-04,-.2262E-06,
     1 .1323E+02, .3396E-01, .1781E-02, .1175E-04,-.3247E-06,
     1 .1643E+02, .4218E-01, .2290E-02, .1606E-04,-.4258E-06,
     1 .2066E+02, .5072E-01, .2954E-02, .2282E-04,-.5743E-06,
     1 .2703E+02, .6713E-01, .4093E-02, .3230E-04,-.8506E-06,
     1 .3367E+02, .9359E-01, .5513E-02, .4178E-04,-.1121E-05,
     1 .4364E+02, .1033E+00, .6674E-02, .5709E-04,-.1403E-05,
     1 .5596E+02, .1693E+00, .9550E-02, .6365E-04,-.2023E-05,
     1 .6287E+02, .4216E+00, .1751E-01, .3710E-04,-.3394E-05,
     2 .5170E-01, .1044E-02, .4268E-04, .1976E-06,-.4816E-08,
     2 .5375E-01, .1045E-02, .4276E-04, .1985E-06,-.4844E-08,
     2 .5968E-01, .1002E-02, .4182E-04, .2101E-06,-.4757E-08,
     2 .6615E-01, .1034E-02, .4324E-04, .2126E-06,-.5081E-08,
     2 .7531E-01, .1173E-02, .4785E-04, .1966E-06,-.6056E-08,
     2 .1018E+00, .1173E-02, .5008E-04, .2383E-06,-.6684E-08,
     2 .1342E+00, .1309E-02, .5669E-04, .2672E-06,-.8100E-08,
     2 .1917E+00, .1682E-02, .7213E-04, .2842E-06,-.1152E-07,
     2 .2972E+00, .1897E-02, .8747E-04, .4546E-06,-.1414E-07,
     2 .3855E+00, .3031E-02, .1288E-03, .5204E-06,-.1954E-07,
     2 .7423E+00, .4535E-02, .1965E-03, .8849E-06,-.2927E-07,
     2 .9859E+00, .4597E-02, .2158E-03, .1331E-05,-.3080E-07,
     2 .1316E+01, .6987E-02, .3086E-03, .1467E-05,-.4673E-07,
     2 .1809E+01, .6752E-02, .3375E-03, .2355E-05,-.5047E-07,
     2 .2343E+01, .9041E-02, .4397E-03, .3010E-05,-.6184E-07,
     2 .3125E+01, .1392E-01, .6355E-03, .3391E-05,-.9717E-07,
     2 .4146E+01, .1739E-01, .8144E-03, .4677E-05,-.1251E-06,
     2 .5620E+01, .2043E-01, .1038E-02, .7100E-05,-.1679E-06,
     2 .7262E+01, .2520E-01, .1319E-02, .9658E-05,-.2124E-06,
     2 .1011E+02, .4856E-01, .2259E-02, .1040E-04,-.4136E-06,
     2 .1330E+02, .7198E-01, .3203E-02, .1265E-04,-.5738E-06,
     2 .1731E+02, .7242E-01, .3471E-02, .1956E-04,-.5831E-06,
     2 .2076E+02, .1096E+00, .4867E-02, .2097E-04,-.8149E-06,
     2 .2346E+02, .1381E+00, .5856E-02, .2247E-04,-.9158E-06,
     2 .2809E+02, .1022E+00, .5081E-02, .3419E-04,-.7825E-06,
     2 .2414E+02, .1368E+00, .5566E-02, .2155E-04,-.7330E-06,
     2 .1476E+02, .1447E+00, .5366E-02, .1282E-04,-.6336E-06,
     2 .5787E+01, .5165E-02, .1446E-02, .3817E-04, .4241E-07,
     3 .6458E+04,-.4833E-01, .3890E-02, .1125E-03,-.1087E-05,
     3 .6458E+04,-.4827E-01, .3888E-02, .1124E-03,-.1086E-05,
     3 .6458E+04,-.4791E-01, .3887E-02, .1123E-03,-.1084E-05,
     3 .6458E+04,-.4778E-01, .3872E-02, .1119E-03,-.1076E-05,
     3 .6458E+04,-.4840E-01, .3821E-02, .1115E-03,-.1062E-05,
     3 .6458E+04, .2823E-02, .5023E-02, .9978E-04,-.1185E-05,
     3 .6458E+04, .2388E-02, .4853E-02, .9690E-04,-.1138E-05,
     3 .6451E+04,-.7123E-02, .4602E-02, .9830E-04,-.1114E-05,
     3 .6458E+04,-.1205E+00, .5893E-02, .1799E-03,-.2481E-05,
     3 .6458E+04,-.1635E+00, .3339E-02, .1668E-03,-.1773E-05,
     3 .6426E+04, .1113E+01, .3695E-01,-.1197E-03,-.7629E-05,
     3 .6419E+04, .1027E+01, .3404E-01,-.1017E-03,-.6733E-05,
     3 .6383E+04, .1558E+01, .4609E-01,-.2414E-03,-.8186E-05,
     3 .6366E+04, .3765E+00, .2042E-01, .1448E-03,-.3753E-05,
     3 .6340E+04, .7142E-01, .3606E-02, .4442E-04,-.8531E-08,
     3 .6333E+04, .2743E+00, .1009E-01, .1786E-04,-.1330E-05,
     3 .6282E+04, .3212E+00, .4309E-02,-.1311E-03, .2304E-06,
     3 .6232E+04,-.2892E+00,-.1575E-01,-.7142E-04, .4072E-05,
     3 .6172E+04,-.9869E+00,-.3712E-01, .1841E-06, .7296E-05,
     3 .6032E+04,-.5739E+00,-.3114E-01,-.1738E-03, .7066E-05,
     3 .5885E+04,-.5738E+00,-.3492E-01,-.2867E-03, .6931E-05,
     3 .5678E+04,-.4250E+00,-.3835E-01,-.4296E-03, .9120E-05,
     3 .5459E+04,-.8243E+00,-.5504E-01,-.4660E-03, .1238E-04,
     3 .5187E+04,-.9670E+00,-.6470E-01,-.5443E-03, .1471E-04,
     3 .4896E+04,-.1658E+01,-.8486E-01,-.4743E-03, .1706E-04,
     3 .4541E+04,-.1835E+01,-.8762E-01,-.4542E-03, .1573E-04,
     3 .4055E+04,-.1467E+01,-.7584E-01,-.4803E-03, .1400E-04,
     3 .3500E+04,-.1173E+01,-.6168E-01,-.3969E-03, .1165E-04 /
      DATA ((( CL4N2OGH(K,J,I), K = 1, 5), J = 1, 28), I = 4, 5) /
     4 .3893E+03, .1021E+01, .3123E-01,-.1210E-03,-.5144E-05,
     4 .3893E+03, .1021E+01, .3123E-01,-.1210E-03,-.5144E-05,
     4 .3893E+03, .1021E+01, .3123E-01,-.1210E-03,-.5144E-05,
     4 .3894E+03, .1020E+01, .3123E-01,-.1210E-03,-.5144E-05,
     4 .3894E+03, .1020E+01, .3122E-01,-.1210E-03,-.5144E-05,
     4 .3931E+03, .9536E+00, .2956E-01,-.1073E-03,-.4952E-05,
     4 .3934E+03, .9516E+00, .2955E-01,-.1060E-03,-.4956E-05,
     4 .3942E+03, .9580E+00, .2943E-01,-.1122E-03,-.4888E-05,
     4 .3876E+03, .1090E+01, .2708E-01,-.2214E-03,-.3004E-05,
     4 .3873E+03, .1126E+01, .2891E-01,-.2179E-03,-.3568E-05,
     4 .4184E+03,-.5248E+00,-.1714E-01, .1022E-03, .4688E-05,
     4 .4248E+03,-.3435E+00,-.1757E-01,-.3228E-04, .5424E-05,
     4 .5268E+03,-.2199E+01,-.6428E-01, .3501E-03, .1122E-04,
     4 .5412E+03,-.7284E+00,-.3442E-01,-.1727E-03, .6188E-05,
     4 .5608E+03,-.5328E+00,-.1857E-01, .4252E-04, .3948E-05,
     4 .5900E+03,-.1583E+01,-.4900E-01, .1828E-03, .8308E-05,
     4 .6188E+03,-.8484E+00,-.2166E-01, .1429E-03, .1974E-05,
     4 .6592E+03, .6440E-01, .1962E-02,-.1083E-03,-.3201E-05,
     4 .6804E+03, .5872E+00, .1542E-01,-.2383E-03,-.5716E-05,
     4 .7504E+03,-.7692E+00,-.1903E-01, .8080E-04,-.1088E-06,
     4 .7328E+03,-.1586E+00,-.1535E-02, .6512E-04,-.5992E-06,
     4 .7340E+03,-.5016E+00,-.5500E-02, .1926E-03,-.1719E-05,
     4 .6772E+03,-.1717E+00, .5508E-02, .1743E-03,-.3544E-05,
     4 .6228E+03,-.1493E+00, .8392E-02, .2060E-03,-.4472E-05,
     4 .5224E+03, .4880E+00, .2377E-01, .8404E-04,-.5656E-05,
     4 .4116E+03, .5196E+00, .2191E-01, .7456E-04,-.3667E-05,
     4 .3674E+03, .3179E+00, .1746E-01, .1248E-03,-.3306E-05,
     4 .3502E+03, .2198E+00, .8108E-02, .1192E-04,-.1179E-05,
     5 .3956E+04,-.9684E+01,-.2756E+00, .1190E-02, .3309E-04,
     5 .3956E+04,-.9688E+01,-.2756E+00, .1190E-02, .3309E-04,
     5 .3955E+04,-.9684E+01,-.2756E+00, .1190E-02, .3310E-04,
     5 .3955E+04,-.9684E+01,-.2756E+00, .1190E-02, .3309E-04,
     5 .3954E+04,-.9684E+01,-.2756E+00, .1190E-02, .3310E-04,
     5 .3952E+04,-.9684E+01,-.2756E+00, .1190E-02, .3311E-04,
     5 .3949E+04,-.9684E+01,-.2755E+00, .1191E-02, .3312E-04,
     5 .3943E+04,-.9684E+01,-.2756E+00, .1195E-02, .3321E-04,
     5 .3930E+04,-.9680E+01,-.2754E+00, .1195E-02, .3328E-04,
     5 .3903E+04,-.9664E+01,-.2752E+00, .1196E-02, .3339E-04,
     5 .3804E+04,-.9592E+01,-.2737E+00, .1193E-02, .3378E-04,
     5 .3734E+04,-.1061E+02,-.2560E+00, .2086E-02, .2102E-04,
     5 .3032E+04, .9632E+00, .3287E-01,-.2678E-03,-.1206E-04,
     5 .2910E+04, .6236E+00, .2294E-01,-.1914E-03,-.9672E-05,
     5 .2743E+04, .1983E+01, .3851E-01,-.1425E-02,-.2592E-04,
     5 .2118E+04, .9328E+01, .2214E+00,-.2838E-02,-.4424E-04,
     5 .1852E+04, .5196E+00,-.2855E-01,-.9160E-03, .9244E-05,
     5 .1217E+04,-.2219E+01,-.7616E-01, .5096E-03, .2496E-04,
     5 .6204E+03,-.1749E+00,-.1348E-01, .5316E-03, .2213E-04,
     5 .1615E+03, .6936E+01, .1902E+00,-.9124E-03,-.1962E-04,
     5 .1904E+03,-.4192E-01, .2323E-02, .6436E-04,-.1037E-05,
     5 .1386E+03, .2807E+00, .5912E-02,-.4012E-04, .5732E-06,
     5 .1503E+03, .4096E+00, .6616E-02,-.1279E-03, .5024E-06,
     5 .2210E+03, .1447E+00, .1008E-02,-.5940E-04, .8608E-06,
     5 .3346E+03, .6364E-01,-.6224E-02,-.1564E-03, .2112E-05,
     5 .5816E+03,-.6336E+00,-.1241E-01, .1678E-03, .1669E-06,
     5 .7484E+03,-.2662E+00, .1446E-01, .4064E-03,-.6384E-05,
     5 .8964E+03, .5260E+00, .6216E-01, .7264E-03,-.1772E-04 /
      DATA ((( CL4CH4GH(K,J,I), K = 1, 5), J = 1, 28), I = 1, 3) /
     1 .8700E+00, .1887E-01, .1585E-03, .6363E-06, .1081E-08,
     1 .8721E+00, .1898E-01, .1413E-03, .5553E-06, .3840E-08,
     1 .8595E+00, .1880E-01, .1489E-03, .5433E-06, .2285E-08,
     1 .8466E+00, .1902E-01, .1696E-03, .4416E-06,-.1869E-08,
     1 .8448E+00, .1890E-01, .1844E-03, .5346E-06,-.3618E-08,
     1 .8460E+00, .1804E-01, .1724E-03, .7254E-06,-.4830E-09,
     1 .8769E+00, .1775E-01, .1611E-03, .7674E-06, .1459E-08,
     1 .9585E+00, .1811E-01, .1305E-03, .8235E-06, .8208E-08,
     1 .1076E+01, .1884E-01, .1309E-03, .7077E-06, .7203E-08,
     1 .1241E+01, .1978E-01, .1883E-03, .6351E-06,-.3513E-08,
     1 .1919E+01, .2328E-01, .1845E-03, .4407E-06,-.3924E-08,
     1 .2375E+01, .2612E-01, .1809E-03, .1315E-06,-.4521E-08,
     1 .2954E+01, .2696E-01, .1611E-03, .3363E-06,-.1368E-08,
     1 .3774E+01, .2908E-01, .4320E-04, .5202E-06, .2121E-07,
     1 .4824E+01, .3621E-01, .7428E-04,-.4026E-06, .1263E-07,
     1 .6174E+01, .4110E-01, .2357E-03,-.7044E-06,-.1811E-07,
     1 .8148E+01, .4674E-01, .2657E-03,-.1123E-06,-.1929E-07,
     1 .1036E+02, .4995E-01, .6102E-03, .7980E-06,-.7266E-07,
     1 .1414E+02, .6255E-01, .2994E-03, .1214E-06,-.3504E-07,
     1 .1903E+02, .7083E-01,-.1424E-03, .2311E-05, .2283E-07,
     1 .2470E+02, .1071E+00, .4068E-03,-.2388E-05,-.5868E-07,
     1 .3108E+02, .1164E+00, .3711E-03,-.2117E-05, .2126E-07,
     1 .3771E+02, .1077E+00, .1159E-02, .3057E-05,-.9174E-07,
     1 .4809E+02, .1572E+00,-.8061E-03,-.7326E-06, .1206E-06,
     1 .5589E+02, .2635E+00,-.2446E-02,-.9432E-05, .3012E-06,
     1 .6153E+02, .3111E+00,-.5085E-03,-.9444E-05, .6441E-07,
     1 .6975E+02, .3501E+00, .4641E-03,-.2692E-05,-.9003E-07,
     1 .8316E+02, .4320E+00,-.2826E-02,-.1446E-04, .4611E-06,
     2 .8589E+00, .1775E-01, .1834E-03, .1063E-05,-.4194E-09,
     2 .8619E+00, .1764E-01, .1781E-03, .1091E-05, .5043E-09,
     2 .8634E+00, .1822E-01, .1967E-03, .1024E-05,-.2585E-08,
     2 .8694E+00, .1804E-01, .1970E-03, .1136E-05,-.1791E-08,
     2 .8784E+00, .1763E-01, .1835E-03, .1063E-05,-.1280E-08,
     2 .8652E+00, .1762E-01, .2120E-03, .1033E-05,-.6234E-08,
     2 .9228E+00, .1781E-01, .2100E-03, .9363E-06,-.7014E-08,
     2 .1009E+01, .1795E-01, .1923E-03, .8424E-06,-.4398E-08,
     2 .1153E+01, .1739E-01, .1936E-03, .1064E-05,-.4563E-08,
     2 .1369E+01, .1763E-01, .1357E-03, .7248E-06, .3159E-08,
     2 .2197E+01, .2318E-01, .1815E-03,-.4185E-06,-.2029E-07,
     2 .2599E+01, .2272E-01, .1598E-03, .9402E-07,-.1427E-07,
     2 .3441E+01, .2298E-01,-.6693E-04, .9123E-06, .2633E-07,
     2 .4131E+01, .2627E-01, .8970E-04, .1018E-05, .4293E-08,
     2 .5148E+01, .2958E-01, .1367E-03, .1379E-05,-.5712E-08,
     2 .6624E+01, .2888E-01,-.8442E-04, .3177E-05, .3993E-07,
     2 .8817E+01, .5097E-01, .5445E-04,-.1598E-05, .3513E-08,
     2 .1220E+02, .6765E-01,-.1743E-03,-.3681E-05, .3231E-07,
     2 .1527E+02, .6414E-01, .6732E-04, .2051E-05, .3408E-07,
     2 .1929E+02, .7989E-01, .3669E-03,-.2651E-05, .3153E-09,
     2 .2435E+02, .7494E-01, .7125E-05, .2423E-05, .4284E-07,
     2 .3306E+02, .8523E-01,-.1250E-02, .6003E-05, .1263E-06,
     2 .4257E+02, .2278E+00,-.4995E-03,-.2163E-04,-.2819E-07,
     2 .5352E+02, .1195E+00, .1703E-02, .8169E-05,-.4698E-06,
     2 .7062E+02, .4803E-01, .1991E-02, .5082E-04,-.8274E-06,
     2 .9615E+02, .2264E+00,-.1710E-01, .4059E-04, .2430E-05,
     2 .9891E+02, .5631E+00,-.1117E-01,-.3456E-05, .1443E-05,
     2 .9564E+02, .8934E+00, .6780E-03,-.4305E-04,-.6852E-06,
     3 .9288E+00, .2126E-01, .1902E-03, .4353E-06,-.2975E-08,
     3 .9306E+00, .2115E-01, .1935E-03, .4560E-06,-.3603E-08,
     3 .9387E+00, .2134E-01, .1973E-03, .4104E-06,-.4506E-08,
     3 .9636E+00, .2127E-01, .1868E-03, .4221E-06,-.3084E-08,
     3 .9810E+00, .2158E-01, .1865E-03, .3675E-06,-.2958E-08,
     3 .1034E+01, .2176E-01, .1823E-03, .3576E-06,-.2314E-08,
     3 .1118E+01, .2164E-01, .1850E-03, .4491E-06,-.2807E-08,
     3 .1232E+01, .2098E-01, .2039E-03, .7506E-06,-.4689E-08,
     3 .1621E+01, .2367E-01, .1262E-03, .3819E-06, .5298E-08,
     3 .2102E+01, .2708E-01, .1778E-03, .2137E-06,-.5181E-08,
     3 .3567E+01, .3681E-01, .2037E-03,-.3219E-07,-.5022E-08,
     3 .4428E+01, .4692E-01, .3123E-03,-.1232E-05,-.2515E-07,
     3 .5760E+01, .5175E-01, .2055E-03,-.1133E-05,-.3948E-08,
     3 .7287E+01, .5418E-01, .2750E-03,-.4224E-06,-.6207E-08,
     3 .9219E+01, .5823E-01, .5517E-03,-.2557E-06,-.3555E-07,
     3 .1274E+02, .6270E-01, .2239E-03,-.3996E-06, .1998E-07,
     3 .1728E+02, .6105E-01, .3735E-03, .9441E-06, .2895E-08,
     3 .2393E+02, .5451E-01, .6534E-03, .2003E-05,-.7485E-07,
     3 .3444E+02, .4383E-01,-.2811E-03, .5136E-05, .7017E-07,
     3 .4809E+02, .4314E-01,-.7395E-03, .6189E-05, .1073E-06,
     3 .6579E+02, .5343E-01,-.1095E-02, .1158E-04, .8841E-07,
     3 .9036E+02, .9243E-01,-.8157E-03, .1553E-04,-.9471E-07,
     3 .1236E+03, .2937E+00,-.4410E-02,-.4323E-06, .4398E-06,
     3 .1607E+03, .5868E+00,-.7134E-02,-.5088E-04, .1144E-05,
     3 .1815E+03, .5172E+00, .1127E-01,-.3252E-04,-.1679E-05,
     3 .2223E+03,-.1233E-01, .5787E-02, .3966E-04,-.4056E-06,
     3 .2791E+03,-.3414E+00,-.1515E-02, .1201E-04, .7755E-07,
     3 .3342E+03,-.3147E+00,-.2453E-02, .3483E-04,-.1173E-06 /
      DATA ((( CL4CH4GH(K,J,I), K = 1, 5), J = 1, 28), I = 4, 5) /
     4 .4757E+04, .3723E+02,-.1232E-01,-.1715E-03, .2729E-05,
     4 .4757E+04, .3723E+02,-.1232E-01,-.1714E-03, .2731E-05,
     4 .4757E+04, .3723E+02,-.1232E-01,-.1714E-03, .2729E-05,
     4 .4758E+04, .3723E+02,-.1231E-01,-.1713E-03, .2728E-05,
     4 .4759E+04, .3722E+02,-.1230E-01,-.1711E-03, .2727E-05,
     4 .4762E+04, .3720E+02,-.1230E-01,-.1710E-03, .2729E-05,
     4 .4768E+04, .3717E+02,-.1225E-01,-.1714E-03, .2737E-05,
     4 .4779E+04, .3709E+02,-.1190E-01,-.1722E-03, .2726E-05,
     4 .4803E+04, .3690E+02,-.1079E-01,-.1692E-03, .2628E-05,
     4 .4856E+04, .3645E+02,-.9280E-02,-.1562E-03, .2479E-05,
     4 .5073E+04, .3486E+02,-.8470E-02,-.1160E-03, .2350E-05,
     4 .5236E+04, .3382E+02,-.1048E-01,-.9670E-04, .2602E-05,
     4 .5475E+04, .3243E+02,-.1390E-01,-.7469E-04, .3056E-05,
     4 .5813E+04, .3057E+02,-.1208E-01,-.3754E-04, .2472E-05,
     4 .6287E+04, .2813E+02,-.1333E-01, .3365E-04, .2126E-05,
     4 .6931E+04, .2547E+02,-.2359E-01, .2483E-04, .3583E-05,
     4 .7753E+04, .2220E+02,-.2430E-01, .1363E-04, .3552E-05,
     4 .8761E+04, .1796E+02,-.2315E-01, .7604E-04, .3398E-05,
     4 .9912E+04, .1323E+02,-.2065E-01, .2217E-03, .1908E-05,
     4 .1108E+05, .9344E+01,-.3527E-01, .1472E-03, .4259E-05,
     4 .1204E+05, .4914E+01,-.1852E-01, .1686E-03, .2102E-05,
     4 .1264E+05, .4016E+00,-.2033E-01, .2019E-03, .3795E-05,
     4 .1273E+05,-.3879E+01,-.4584E-02, .2450E-03, .1802E-05,
     4 .1231E+05,-.7157E+01,-.1190E-01, .1783E-03, .4332E-05,
     4 .1146E+05,-.1005E+02, .1147E-01, .1396E-03, .1857E-05,
     4 .1039E+05,-.1248E+02, .1441E-01, .1426E-03, .2535E-05,
     4 .9251E+04,-.1326E+02, .1343E-01, .8996E-04, .3177E-05,
     4 .8208E+04,-.1274E+02, .1830E-01, .6832E-04, .2486E-05,
     5 .2579E+06,-.5724E+03,-.1749E+00, .7040E-02, .5157E-04,
     5 .2579E+06,-.5723E+03,-.1754E+00, .7039E-02, .5167E-04,
     5 .2579E+06,-.5723E+03,-.1755E+00, .7036E-02, .5167E-04,
     5 .2579E+06,-.5722E+03,-.1754E+00, .7035E-02, .5164E-04,
     5 .2579E+06,-.5721E+03,-.1755E+00, .7034E-02, .5164E-04,
     5 .2578E+06,-.5718E+03,-.1757E+00, .7031E-02, .5163E-04,
     5 .2577E+06,-.5712E+03,-.1768E+00, .7020E-02, .5172E-04,
     5 .2574E+06,-.5701E+03,-.1776E+00, .7002E-02, .5164E-04,
     5 .2568E+06,-.5679E+03,-.1814E+00, .6972E-02, .5195E-04,
     5 .2557E+06,-.5631E+03,-.1885E+00, .6883E-02, .5277E-04,
     5 .2517E+06,-.5460E+03,-.2030E+00, .6495E-02, .5515E-04,
     5 .2487E+06,-.5342E+03,-.1863E+00, .6232E-02, .5273E-04,
     5 .2444E+06,-.5180E+03,-.1508E+00, .5899E-02, .4738E-04,
     5 .2385E+06,-.4955E+03,-.1888E+00, .5289E-02, .5658E-04,
     5 .2302E+06,-.4653E+03,-.1771E+00, .4150E-02, .6117E-04,
     5 .2188E+06,-.4324E+03,-.3606E-01, .4062E-02, .4050E-04,
     5 .2039E+06,-.3916E+03,-.6078E-02, .3974E-02, .3790E-04,
     5 .1853E+06,-.3383E+03, .1536E-01, .2900E-02, .3585E-04,
     5 .1632E+06,-.2789E+03, .9858E-02, .6367E-03, .5666E-04,
     5 .1386E+06,-.2343E+03, .2216E+00, .1417E-02, .2702E-04,
     5 .1138E+06,-.1868E+03, .1016E+00, .7532E-03, .4567E-04,
     5 .9064E+05,-.1440E+03, .1976E+00, .7790E-04, .2259E-04,
     5 .7047E+05,-.1085E+03, .1522E+00,-.5594E-03, .2874E-04,
     5 .5398E+05,-.8334E+02, .2023E+00,-.6550E-03, .1868E-04,
     5 .4124E+05,-.6157E+02, .1818E+00,-.1532E-02, .2463E-04,
     5 .3195E+05,-.5278E+02, .2321E+00,-.6407E-03, .9920E-05,
     5 .2647E+05,-.4635E+02, .1244E+00,-.3183E-03, .7762E-05,
     5 .2338E+05,-.4074E+02, .5421E-01,-.1589E-03, .9273E-05 /
C
C----------------------------------------------------------------------C
C     SPECTRAL REGION: 980-1100 CM^-1. GAS: H2O, O3                    C
C----------------------------------------------------------------------C
C
      COMMON /BANDL5/ GWL5(2), CL5H2O(5,18,2), CL5O3(5,18,2),
     1                         CL5CS(5,4,2), CL5CF(5,4,2),
     2                         CL5F11(2), CL5F12(2)
      DATA GWL5 / .59500, .27000 /
C
      DATA (CL5F11(I), I = 1, 2) / 1526.4, 152.1 /
      DATA (CL5F12(I), I = 1, 2) / 1404.3, 12.13 /
C
      DATA ((( CL5H2O(K,J,I), K = 1, 5), J = 1, 18), I = 1, 2) /
     1 .1199E-04, .3452E-06, .3698E-08, .1688E-10, .2669E-13,
     1 .1776E-04, .5178E-06, .5692E-08, .2394E-10,-.9638E-15,
     1 .2650E-04, .7678E-06, .8330E-08, .3730E-10, .3919E-13,
     1 .3945E-04, .1140E-05, .1221E-07, .5509E-10, .8179E-13,
     1 .5863E-04, .1697E-05, .1830E-07, .8018E-10, .7755E-13,
     1 .8687E-04, .2499E-05, .2670E-07, .1238E-09, .2242E-12,
     1 .1281E-03, .3732E-05, .4083E-07, .1763E-09, .7845E-13,
     1 .1888E-03, .5512E-05, .6031E-07, .2707E-09, .2675E-12,
     1 .2760E-03, .8048E-05, .8772E-07, .4009E-09, .5227E-12,
     1 .3990E-03, .1166E-04, .1275E-06, .5941E-09, .9192E-12,
     1 .5720E-03, .1687E-04, .1881E-06, .8442E-09, .5834E-12,
     1 .8128E-03, .2387E-04, .2675E-06, .1277E-08, .1619E-11,
     1 .1148E-02, .3365E-04, .3711E-06, .1627E-08, .9590E-12,
     1 .1593E-02, .4704E-04, .5275E-06, .2319E-08, .8217E-12,
     1 .2157E-02, .6330E-04, .7332E-06, .3460E-08, .8971E-12,
     1 .2896E-02, .8347E-04, .9484E-06, .4549E-08, .3177E-11,
     1 .3820E-02, .1071E-03, .1158E-05, .5380E-08, .5630E-11,
     1 .4987E-02, .1364E-03, .1349E-05, .6183E-08, .1832E-10,
     2 .1105E-04, .3303E-06, .3827E-08, .1485E-10,-.3545E-13,
     2 .1649E-04, .4848E-06, .5454E-08, .2438E-10, .7586E-14,
     2 .2440E-04, .7254E-06, .8353E-08, .3456E-10,-.4444E-13,
     2 .3615E-04, .1091E-05, .1294E-07, .5181E-10,-.1218E-12,
     2 .5385E-04, .1606E-05, .1879E-07, .8369E-10,-.4627E-13,
     2 .8044E-04, .2408E-05, .2803E-07, .1162E-09,-.1638E-12,
     2 .1199E-03, .3529E-05, .4009E-07, .1841E-09, .9172E-13,
     2 .1801E-03, .5244E-05, .5750E-07, .2597E-09, .2938E-12,
     2 .2675E-03, .7872E-05, .8875E-07, .4029E-09, .2491E-12,
     2 .3930E-03, .1162E-04, .1349E-06, .6092E-09,-.1131E-12,
     2 .5813E-03, .1676E-04, .1875E-06, .9503E-09, .1692E-11,
     2 .8444E-03, .2467E-04, .2840E-06, .1246E-08,-.8799E-12,
     2 .1225E-02, .3499E-04, .3994E-06, .2118E-08, .3382E-11,
     2 .1822E-02, .5009E-04, .5195E-06, .2921E-08, .1109E-10,
     2 .2649E-02, .7517E-04, .7642E-06, .3027E-08, .3209E-11,
     2 .3729E-02, .1078E-03, .1124E-05, .3997E-08,-.4317E-11,
     2 .4991E-02, .1468E-03, .1554E-05, .6573E-08, .8653E-11,
     2 .6182E-02, .1898E-03, .2226E-05, .9410E-08,-.7934E-11 /
      DATA ((( CL5O3(K,J,I), K = 1, 5), J = 1, 18), I = 1, 2) /
     1 .7227E+01, .6737E-01, .2301E-03, .9636E-07,-.3997E-08,
     1 .9629E+01, .8488E-01, .2596E-03,-.1467E-07,-.5302E-08,
     1 .1295E+02, .1082E+00, .2804E-03,-.1549E-06,-.4808E-08,
     1 .1755E+02, .1392E+00, .2985E-03,-.2718E-06,-.4402E-08,
     1 .2382E+02, .1807E+00, .3017E-03,-.4668E-06,-.1571E-08,
     1 .3224E+02, .2364E+00, .3046E-03,-.9947E-06, .7062E-09,
     1 .4321E+02, .3102E+00, .3024E-03,-.1909E-05, .7458E-08,
     1 .5702E+02, .4064E+00, .4323E-03,-.3237E-05,-.2147E-08,
     1 .7388E+02, .5251E+00, .6050E-03,-.4465E-05,-.9082E-08,
     1 .9364E+02, .6619E+00, .8856E-03,-.4863E-05,-.2401E-07,
     1 .1161E+03, .8160E+00, .1052E-02,-.5650E-05,-.1708E-07,
     1 .1419E+03, .9792E+00, .1153E-02,-.6110E-05,-.8594E-08,
     1 .1716E+03, .1144E+01, .1234E-02,-.5189E-05,-.2551E-07,
     1 .2053E+03, .1309E+01, .1311E-02,-.2540E-05,-.9645E-07,
     1 .2413E+03, .1505E+01, .6268E-03,-.2568E-05,-.5036E-07,
     1 .2762E+03, .1718E+01,-.2570E-03,-.1056E-05,-.1939E-07,
     1 .3056E+03, .1961E+01,-.8495E-03,-.6960E-05, .2950E-07,
     1 .3272E+03, .2169E+01,-.6410E-03,-.1738E-04, .3894E-07,
     2 .1420E+03, .1113E+01, .3418E-02, .5459E-05, .4714E-07,
     2 .1711E+03, .1199E+01, .3319E-02, .3513E-05, .4462E-07,
     2 .2093E+03, .1323E+01, .3506E-02,-.1264E-05, .1046E-07,
     2 .2597E+03, .1463E+01, .3852E-02,-.3441E-05,-.5330E-07,
     2 .3260E+03, .1631E+01, .3808E-02,-.4302E-05,-.4033E-07,
     2 .4122E+03, .1813E+01, .3875E-02,-.4645E-05,-.8157E-07,
     2 .5227E+03, .1998E+01, .3714E-02,-.7030E-05,-.1000E-06,
     2 .6645E+03, .2100E+01, .3737E-02,-.2429E-05,-.1654E-06,
     2 .8466E+03, .2074E+01, .3199E-02, .5591E-05,-.2743E-06,
     2 .1078E+04, .1932E+01,-.5908E-03, .1557E-04,-.2044E-07,
     2 .1353E+04, .1757E+01,-.3764E-02, .1058E-04, .8142E-07,
     2 .1656E+04, .1481E+01,-.5733E-02, .9746E-05, .3976E-07,
     2 .1959E+04, .1076E+01,-.5487E-02, .1454E-04,-.1805E-06,
     2 .2239E+04, .4363E+00,-.6362E-02, .2840E-04,-.1221E-06,
     2 .2475E+04,-.2581E+00,-.5978E-02, .1979E-04,-.1271E-07,
     2 .2664E+04,-.9293E+00,-.9063E-02, .2347E-04, .4743E-06,
     2 .2803E+04,-.1565E+01,-.4884E-02, .3316E-04,-.3786E-06,
     2 .2915E+04,-.2057E+01,-.4014E-02, .4109E-04,-.6772E-06 /
      DATA ((( CL5CS(K,J,I), K = 1, 5), J = 1, 4), I = 1, 2) /
     1 .4632E+01,-.1078E+00, .1342E-02,-.1431E-04, .8921E-07,
     1 .6883E+01,-.1601E+00, .1991E-02,-.2127E-04, .1330E-06,
     1 .1023E+02,-.2378E+00, .2957E-02,-.3161E-04, .1978E-06,
     1 .1520E+02,-.3535E+00, .4395E-02,-.4686E-04, .2923E-06,
     2 .4540E+01,-.1051E+00, .1301E-02,-.1384E-04, .8638E-07,
     2 .6735E+01,-.1564E+00, .1939E-02,-.2051E-04, .1270E-06,
     2 .9997E+01,-.2325E+00, .2885E-02,-.3044E-04, .1882E-06,
     2 .1485E+02,-.3453E+00, .4284E-02,-.4530E-04, .2808E-06 /
      DATA ((( CL5CF(K,J,I), K = 1, 5), J = 1, 4), I = 1, 2) /
     1 .2729E-04,-.1128E-06, .4232E-09,-.1910E-11, .7769E-14,
     1 .4053E-04,-.1672E-06, .6241E-09,-.2850E-11, .1225E-13,
     1 .6017E-04,-.2479E-06, .9307E-09,-.4219E-11, .1739E-13,
     1 .8929E-04,-.3675E-06, .1386E-08,-.6293E-11, .2534E-13,
     2 .2415E-04,-.9788E-07, .3696E-09,-.1732E-11, .7413E-14,
     2 .3575E-04,-.1453E-06, .5643E-09,-.2588E-11, .8712E-14,
     2 .5297E-04,-.2160E-06, .8246E-09,-.3847E-11, .1521E-13,
     2 .7852E-04,-.3215E-06, .1206E-08,-.5637E-11, .2550E-13 /
C
      COMMON /BANDL5GH/ GWL5GH(4), CL5H2OGH(5,28,4), CL5O3GH(5,28,4),
     1                             CL5CSGH(5,4,2), CL5CFGH(5,4,2)
      DATA GWL5GH / .06800, .04000, .02000, .00700 /
      DATA ((( CL5H2OGH(K,J,I), K = 1, 5), J = 1, 28), I = 1, 2) /
     1-.7654E-07,-.4227E-10, .7509E-10, .2732E-13,-.1185E-13,
     1-.2006E-07, .2266E-09, .3211E-10, .4164E-13,-.4422E-14,
     1 .2940E-07, .8358E-09, .8166E-11, .2437E-13,-.2927E-16,
     1 .5849E-07, .1670E-08, .1713E-10, .5471E-13,-.1257E-15,
     1 .1186E-06, .3448E-08, .3590E-10, .9692E-13,-.5131E-15,
     1 .2457E-06, .6849E-08, .6347E-10, .2321E-12, .7383E-15,
     1 .4934E-06, .1358E-07, .1248E-09, .5416E-12, .2615E-14,
     1 .1005E-05, .2888E-07, .2864E-09, .8145E-12,-.1790E-14,
     1 .2063E-05, .5785E-07, .5505E-09, .1925E-11, .3449E-14,
     1 .4126E-05, .1123E-06, .1066E-08, .5289E-11, .2494E-13,
     1 .1213E-04, .3413E-06, .3351E-08, .1183E-10, .8867E-14,
     1 .1827E-04, .5080E-06, .4635E-08, .1788E-10, .7838E-13,
     1 .2686E-04, .7628E-06, .7242E-08, .2259E-10, .2242E-13,
     1 .3961E-04, .1110E-05, .1056E-07, .3696E-10, .6435E-13,
     1 .5720E-04, .1690E-05, .1799E-07, .4245E-10,-.3654E-12,
     1 .8458E-04, .2599E-05, .2917E-07, .6167E-10,-.7969E-12,
     1 .1284E-03, .3760E-05, .3775E-07, .1339E-09, .1346E-12,
     1 .1771E-03, .5522E-05, .6646E-07, .2150E-09,-.1313E-11,
     1 .2619E-03, .8018E-05, .9340E-07, .2641E-09,-.2127E-11,
     1 .3956E-03, .1205E-04, .1389E-06, .3524E-09,-.3647E-11,
     1 .5623E-03, .1707E-04, .1992E-06, .5093E-09,-.5530E-11,
     1 .8228E-03, .2325E-04, .2358E-06, .9902E-09, .1748E-11,
     1 .1191E-02, .3244E-04, .2971E-06, .1749E-08, .1240E-10,
     1 .1610E-02, .4717E-04, .5064E-06, .1976E-08,-.3265E-12,
     1 .2351E-02, .6762E-04, .6379E-06, .2052E-08, .4956E-11,
     1 .3566E-02, .1057E-03, .9799E-06, .1821E-08,-.4604E-11,
     1 .4778E-02, .1530E-03, .1767E-05, .1962E-08,-.7291E-10,
     1 .6378E-02, .2029E-03, .2534E-05, .9835E-08,-.3423E-10,
     2 .4716E+00, .1382E-01, .1508E-03, .6249E-06, .7696E-10,
     2 .4737E+00, .1388E-01, .1515E-03, .6278E-06, .7874E-10,
     2 .4758E+00, .1394E-01, .1521E-03, .6306E-06, .7832E-10,
     2 .4779E+00, .1401E-01, .1528E-03, .6333E-06, .7815E-10,
     2 .4800E+00, .1407E-01, .1535E-03, .6361E-06, .7898E-10,
     2 .4821E+00, .1413E-01, .1541E-03, .6389E-06, .7892E-10,
     2 .4841E+00, .1419E-01, .1548E-03, .6416E-06, .8011E-10,
     2 .4862E+00, .1425E-01, .1555E-03, .6443E-06, .8048E-10,
     2 .4882E+00, .1430E-01, .1560E-03, .6471E-06, .8212E-10,
     2 .4902E+00, .1436E-01, .1567E-03, .6497E-06, .8495E-10,
     2 .4915E+00, .1441E-01, .1572E-03, .6503E-06, .5395E-10,
     2 .4933E+00, .1446E-01, .1577E-03, .6542E-06, .7476E-10,
     2 .4955E+00, .1450E-01, .1578E-03, .6581E-06, .1603E-09,
     2 .4970E+00, .1461E-01, .1602E-03, .6482E-06,-.2112E-09,
     2 .4997E+00, .1461E-01, .1588E-03, .6692E-06, .2870E-09,
     2 .5010E+00, .1462E-01, .1584E-03, .6711E-06, .3862E-09,
     2 .5018E+00, .1474E-01, .1613E-03, .6482E-06,-.2112E-09,
     2 .5019E+00, .1472E-01, .1606E-03, .6691E-06, .1435E-09,
     2 .5009E+00, .1469E-01, .1604E-03, .6667E-06, .1074E-09,
     2 .4986E+00, .1464E-01, .1600E-03, .6653E-06, .9096E-10,
     2 .4943E+00, .1450E-01, .1586E-03, .6593E-06, .6827E-10,
     2 .4860E+00, .1431E-01, .1577E-03, .6476E-06,-.1598E-09,
     2 .4790E+00, .1403E-01, .1519E-03, .6259E-06, .1653E-09,
     2 .4669E+00, .1377E-01, .1509E-03, .6118E-06,-.1115E-09,
     2 .4504E+00, .1330E-01, .1455E-03, .5929E-06,-.2088E-10,
     2 .4307E+00, .1271E-01, .1392E-03, .5669E-06,-.5525E-10,
     2 .4018E+00, .1189E-01, .1308E-03, .5330E-06,-.8155E-10,
     2 .3688E+00, .1094E-01, .1205E-03, .4936E-06,-.6727E-10 /
      DATA ((( CL5H2OGH(K,J,I), K = 1, 5), J = 1, 28), I = 3, 4) /
     3 .1975E-02, .8771E-04, .1749E-05, .1968E-07, .9965E-10,
     3 .1976E-02, .8773E-04, .1749E-05, .1968E-07, .9966E-10,
     3 .1977E-02, .8776E-04, .1749E-05, .1968E-07, .9966E-10,
     3 .1979E-02, .8782E-04, .1749E-05, .1968E-07, .9966E-10,
     3 .1984E-02, .8794E-04, .1750E-05, .1968E-07, .9967E-10,
     3 .1993E-02, .8818E-04, .1752E-05, .1968E-07, .9970E-10,
     3 .2012E-02, .8868E-04, .1754E-05, .1968E-07, .9980E-10,
     3 .2035E-02, .8970E-04, .1776E-05, .1967E-07, .9737E-10,
     3 .2097E-02, .9175E-04, .1805E-05, .1966E-07, .9501E-10,
     3 .2226E-02, .9607E-04, .1862E-05, .1965E-07, .9057E-10,
     3 .2728E-02, .1056E-03, .1768E-05, .2223E-07, .1534E-09,
     3 .3097E-02, .1184E-03, .1895E-05, .1938E-07, .1089E-09,
     3 .2434E-02, .1423E-03, .3039E-05, .1606E-07,-.7396E-10,
     3 .3276E-02, .3815E-04,-.6392E-06, .4120E-07, .6979E-09,
     3 .1292E-02, .1253E-03, .3249E-05, .1574E-08,-.3267E-09,
     3 .1511E-02, .1822E-03, .4974E-05, .7616E-09,-.5322E-09,
     3 .2065E-02, .3523E-05,-.2470E-06, .5055E-07, .7047E-09,
     3 .2866E-02, .8566E-04, .1400E-05, .7920E-08,-.3194E-10,
     3 .4800E-02, .1201E-03, .1055E-05, .1082E-07, .1046E-09,
     3 .6333E-02, .1705E-03, .1582E-05,-.2512E-08,-.1049E-09,
     3 .9130E-02, .2425E-03, .1998E-05, .9971E-08, .1034E-09,
     3 .1613E-01, .3682E-03, .8668E-06, .1063E-07, .4306E-09,
     3 .1386E-01, .4767E-03, .7866E-05, .3283E-07,-.2597E-09,
     3 .1393E-01, .3021E-03, .2470E-05, .2731E-07, .2545E-09,
     3 .1274E-01, .2569E-03, .2541E-05, .2895E-07, .1749E-09,
     3 .4975E-02, .9229E-04, .1107E-05, .3766E-07, .4034E-09,
     3 .6847E-02, .1176E-03, .1078E-05, .3656E-07, .4049E-09,
     3 .8163E-02, .1578E-03, .1467E-05, .8812E-08, .8120E-12,
     4 .1003E-03, .2113E-05, .1250E-07,-.7693E-11,-.1713E-12,
     4 .1003E-03, .2115E-05, .1251E-07,-.7643E-11,-.1707E-12,
     4 .1005E-03, .2119E-05, .1254E-07,-.7440E-11,-.1691E-12,
     4 .1008E-03, .2127E-05, .1260E-07,-.7069E-11,-.1653E-12,
     4 .1015E-03, .2144E-05, .1272E-07,-.6313E-11,-.1584E-12,
     4 .1029E-03, .2177E-05, .1297E-07,-.4854E-11,-.1439E-12,
     4 .1057E-03, .2244E-05, .1347E-07,-.1776E-11,-.1144E-12,
     4 .1114E-03, .2382E-05, .1449E-07, .4365E-11,-.5653E-13,
     4 .1229E-03, .2660E-05, .1656E-07, .1684E-10, .6300E-13,
     4 .1463E-03, .3222E-05, .2065E-07, .4198E-10, .3147E-12,
     4 .2351E-03, .5446E-05, .3916E-07, .1354E-09, .8846E-12,
     4 .3027E-03, .7053E-05, .4966E-07, .2098E-09, .1853E-11,
     4 .4006E-03, .9441E-05, .6829E-07, .3165E-09, .2729E-11,
     4 .5457E-03, .1170E-04, .7246E-07, .1601E-08, .2115E-10,
     4 .6958E-03, .2474E-04, .4081E-06, .6539E-09,-.2852E-10,
     4 .9830E-03, .3472E-04, .5765E-06, .1280E-08,-.3606E-10,
     4 .1416E-02, .5141E-04, .8230E-06, .1365E-08,-.5324E-10,
     4 .2046E-02, .7889E-04, .1119E-05, .1223E-08,-.6010E-10,
     4 .1965E-02, .1172E-03, .2723E-05, .1967E-08,-.2599E-09,
     4 .2618E-02, .2531E-04, .8428E-07, .3269E-07, .4267E-09,
     4 .3062E-02, .8200E-04, .1184E-05, .3648E-08,-.6129E-10,
     4 .1906E-02, .2289E-04, .1368E-06, .1956E-07, .2512E-09,
     4 .2480E-02, .3805E-04, .3607E-06, .9217E-08, .8740E-10,
     4 .3818E-02, .3690E-04, .3426E-06, .7042E-08, .1646E-10,
     4 .7511E-02, .9734E-04, .2720E-07, .4261E-09, .3096E-10,
     4 .1132E-01, .1903E-03, .4145E-06,-.1841E-07,-.1902E-09,
     4 .1293E-01, .2035E-03,-.2079E-06,-.1947E-07,-.9651E-10,
     4 .1252E-01, .1719E-03,-.6563E-06,-.1314E-07, .1552E-10 /
      DATA ((( CL5O3GH(K,J,I), K = 1, 5), J = 1, 28), I = 1, 2) /
     1 .1235E+04, .1437E+02, .3476E-01,-.1760E-03, .2748E-07,
     1 .1244E+04, .1450E+02, .3664E-01,-.1823E-03,-.2837E-06,
     1 .1256E+04, .1459E+02, .3672E-01,-.1763E-03,-.2575E-06,
     1 .1267E+04, .1472E+02, .3858E-01,-.1803E-03,-.5642E-06,
     1 .1277E+04, .1485E+02, .3967E-01,-.1827E-03,-.6947E-06,
     1 .1292E+04, .1494E+02, .3904E-01,-.1772E-03,-.5306E-06,
     1 .1309E+04, .1501E+02, .3750E-01,-.1704E-03,-.1571E-06,
     1 .1330E+04, .1512E+02, .3972E-01,-.1833E-03,-.5111E-06,
     1 .1365E+04, .1513E+02, .3958E-01,-.1923E-03,-.3143E-06,
     1 .1427E+04, .1492E+02, .4083E-01,-.1917E-03,-.2172E-06,
     1 .1668E+04, .1372E+02, .4072E-01,-.1408E-03,-.6785E-06,
     1 .1853E+04, .1292E+02, .3722E-01,-.9009E-04,-.8662E-06,
     1 .2117E+04, .1202E+02, .2847E-01,-.6750E-04,-.3243E-06,
     1 .2484E+04, .1080E+02, .1706E-01,-.2108E-04, .3646E-06,
     1 .2969E+04, .9694E+01, .5055E-02,-.4960E-04, .1103E-05,
     1 .3579E+04, .8069E+01, .6587E-02,-.5096E-04,-.5986E-06,
     1 .4328E+04, .5771E+01,-.7496E-02,-.1220E-04, .2698E-06,
     1 .5148E+04, .3211E+01,-.4069E-02, .4743E-04,-.1966E-05,
     1 .5954E+04, .1135E+01,-.9591E-02, .3693E-05,-.2039E-05,
     1 .6605E+04,-.1389E+01,-.7061E-02, .4273E-04,-.2449E-05,
     1 .7002E+04,-.3843E+01, .5711E-03, .3050E-04,-.2514E-05,
     1 .7122E+04,-.6037E+01,-.1215E-01, .1899E-05, .1128E-05,
     1 .6973E+04,-.7800E+01, .2080E-03,-.3970E-04, .1924E-06,
     1 .6627E+04,-.9133E+01, .6821E-02, .4145E-04,-.5536E-06,
     1 .6188E+04,-.9111E+01, .3886E-02, .2176E-04, .3929E-06,
     1 .5749E+04,-.8867E+01, .7791E-02, .1619E-04, .2636E-06,
     1 .5424E+04,-.9020E+01, .1058E-01, .4350E-04,-.1057E-07,
     1 .5168E+04,-.9070E+01, .9129E-02, .4971E-04, .3295E-06,
     2 .7157E+02, .8460E+00, .6690E-02,-.2680E-05,-.8151E-06,
     2 .7449E+02, .8585E+00, .4937E-02, .1631E-05,-.4676E-06,
     2 .7553E+02, .9243E+00, .5326E-02,-.7469E-05,-.5097E-06,
     2 .7747E+02, .9323E+00, .3551E-02,-.4554E-05,-.1456E-06,
     2 .7949E+02, .9316E+00, .2689E-02,-.3205E-05, .2839E-08,
     2 .7936E+02, .9821E+00, .4394E-02,-.1273E-04,-.2897E-06,
     2 .8052E+02, .1046E+01, .6462E-02,-.2437E-04,-.6985E-06,
     2 .8358E+02, .1014E+01, .4208E-02,-.1568E-04,-.2126E-06,
     2 .8608E+02, .1023E+01, .5603E-02,-.1633E-04,-.4437E-06,
     2 .9581E+02, .1002E+01, .5243E-02,-.1903E-04,-.5671E-06,
     2 .1127E+03, .9163E+00, .4838E-02,-.1098E-04,-.5042E-06,
     2 .1285E+03, .9922E+00, .4078E-02,-.2246E-04,-.2863E-06,
     2 .1508E+03, .9634E+00, .4752E-02,-.5648E-05,-.4935E-06,
     2 .1786E+03, .1134E+01, .5403E-02,-.3748E-04,-.5169E-06,
     2 .2172E+03, .9412E+00, .9844E-02,-.3387E-05,-.1518E-05,
     2 .2727E+03, .8324E+00, .9425E-03, .2089E-04,-.1031E-07,
     2 .3342E+03, .8265E+00, .4453E-02, .4150E-04,-.7344E-06,
     2 .4196E+03, .1106E+01,-.1299E-01,-.7813E-05, .2289E-05,
     2 .5050E+03, .1058E+01,-.1588E-01,-.6986E-05, .2980E-05,
     2 .5807E+03, .1296E+01,-.8803E-02,-.6923E-04, .2207E-05,
     2 .6610E+03, .1030E+01,-.8179E-02,-.4433E-05, .1702E-05,
     2 .7302E+03, .9501E+00, .7137E-02, .4027E-04,-.1734E-05,
     2 .8080E+03, .1145E+01,-.1154E-01, .5962E-04, .2895E-06,
     2 .8706E+03, .2410E+01,-.1576E-01,-.1432E-03, .1547E-05,
     2 .9180E+03, .2359E+01,-.4957E-02,-.8643E-04,-.6621E-06,
     2 .9676E+03, .1970E+01, .1470E-01,-.1135E-03,-.2863E-05,
     2 .1015E+04, .1572E+01,-.5877E-04,-.1230E-03, .9117E-06,
     2 .1059E+04, .9647E+00,-.3888E-02,-.1913E-04, .1558E-05 /
      DATA ((( CL5O3GH(K,J,I), K = 1, 5), J = 1, 28), I = 3, 4) /
     3 .2169E+05, .7335E+02,-.4330E+00, .1545E-03, .1892E-05,
     3 .2185E+05, .7391E+02,-.4363E+00, .1557E-03, .1912E-05,
     3 .2202E+05, .7445E+02,-.4395E+00, .1569E-03, .1923E-05,
     3 .2219E+05, .7499E+02,-.4427E+00, .1580E-03, .1927E-05,
     3 .2236E+05, .7547E+02,-.4456E+00, .1591E-03, .1939E-05,
     3 .2253E+05, .7590E+02,-.4484E+00, .1596E-03, .1956E-05,
     3 .2272E+05, .7617E+02,-.4501E+00, .1599E-03, .1919E-05,
     3 .2290E+05, .7612E+02,-.4499E+00, .1607E-03, .1741E-05,
     3 .2313E+05, .7545E+02,-.4463E+00, .1577E-03, .1503E-05,
     3 .2340E+05, .7348E+02,-.4362E+00, .1346E-03, .1458E-05,
     3 .2398E+05, .6486E+02,-.4009E+00,-.8347E-04, .4101E-05,
     3 .2441E+05, .5864E+02,-.3582E+00,-.1835E-03, .2996E-05,
     3 .2493E+05, .4948E+02,-.2987E+00,-.1936E-03, .5015E-06,
     3 .2553E+05, .3768E+02,-.2377E+00,-.1325E-03,-.2168E-05,
     3 .2617E+05, .2330E+02,-.2049E+00,-.4649E-05, .8632E-06,
     3 .2668E+05, .8568E+01,-.1820E+00, .5630E-04, .4441E-05,
     3 .2677E+05,-.4442E+01,-.1110E+00,-.1121E-03, .5058E-06,
     3 .2623E+05,-.1553E+02,-.8104E-01,-.1024E-03, .1870E-05,
     3 .2487E+05,-.2402E+02,-.4247E-01, .1475E-03,-.6166E-06,
     3 .2274E+05,-.2679E+02,-.3230E-01,-.2455E-04, .7648E-06,
     3 .2003E+05,-.2667E+02,-.2679E-01,-.1149E-03, .2611E-05,
     3 .1711E+05,-.2494E+02, .4987E-02,-.1531E-03, .4351E-06,
     3 .1436E+05,-.2266E+02, .1446E-01,-.1134E-03, .7320E-06,
     3 .1204E+05,-.1987E+02, .1266E-01,-.9089E-04, .1351E-05,
     3 .1034E+05,-.1722E+02, .6935E-02,-.7658E-04, .1111E-05,
     3 .9135E+04,-.1510E+02, .3963E-02,-.1003E-03, .8410E-06,
     3 .8254E+04,-.1232E+02,-.1244E-02,-.2059E-03, .1421E-05,
     3 .7687E+04,-.1028E+02, .2411E-02,-.3823E-03, .2191E-05,
     4 .1576E+06,-.5123E+03, .4692E+00, .1804E-02,-.5527E-05,
     4 .1576E+06,-.5122E+03, .4689E+00, .1804E-02,-.5497E-05,
     4 .1576E+06,-.5122E+03, .4684E+00, .1802E-02,-.5441E-05,
     4 .1576E+06,-.5119E+03, .4682E+00, .1804E-02,-.5473E-05,
     4 .1575E+06,-.5116E+03, .4672E+00, .1804E-02,-.5457E-05,
     4 .1574E+06,-.5110E+03, .4657E+00, .1806E-02,-.5507E-05,
     4 .1572E+06,-.5096E+03, .4611E+00, .1810E-02,-.5291E-05,
     4 .1568E+06,-.5070E+03, .4528E+00, .1807E-02,-.5027E-05,
     4 .1558E+06,-.5016E+03, .4356E+00, .1805E-02,-.4340E-05,
     4 .1539E+06,-.4910E+03, .4092E+00, .1793E-02,-.4169E-05,
     4 .1471E+06,-.4547E+03, .3347E+00, .1817E-02,-.6128E-05,
     4 .1426E+06,-.4299E+03, .2583E+00, .1659E-02,-.2025E-05,
     4 .1362E+06,-.3969E+03, .1668E+00, .1383E-02, .3499E-05,
     4 .1278E+06,-.3562E+03, .9654E-01, .9971E-03, .7272E-05,
     4 .1169E+06,-.3081E+03, .9909E-01, .5652E-03,-.6837E-06,
     4 .1039E+06,-.2565E+03, .1040E+00, .1636E-03,-.5338E-05,
     4 .8923E+05,-.2060E+03, .2627E-01, .1960E-03, .4067E-05,
     4 .7390E+05,-.1598E+03, .5783E-01,-.1043E-03,-.7839E-06,
     4 .5917E+05,-.1206E+03, .4088E-01,-.6056E-03, .6282E-05,
     4 .4623E+05,-.9289E+02, .5606E-01,-.3440E-03, .4447E-05,
     4 .3570E+05,-.7021E+02, .7477E-01,-.2703E-03,-.6300E-06,
     4 .2761E+05,-.5355E+02, .5410E-01,-.1901E-03,-.7134E-06,
     4 .2162E+05,-.4037E+02, .4764E-01,-.3343E-03, .1002E-05,
     4 .1745E+05,-.3323E+02, .4941E-01,-.2557E-03,-.7868E-07,
     4 .1442E+05,-.2881E+02, .4497E-01,-.2317E-03, .1581E-05,
     4 .1220E+05,-.2478E+02, .4946E-01,-.1645E-03,-.3950E-06,
     4 .1060E+05,-.2295E+02, .4247E-01, .6230E-04,-.1118E-05,
     4 .9390E+04,-.2105E+02, .2323E-01, .1997E-03,-.6701E-06 /
      DATA ((( CL5CSGH(K,J,I), K = 1, 5), J = 1, 4), I = 1, 2) /
     1 .4348E+01,-.1005E+00, .1239E-02,-.1314E-04, .8208E-07,
     1 .6415E+01,-.1484E+00, .1839E-02,-.1945E-04, .1202E-06,
     1 .9472E+01,-.2190E+00, .2706E-02,-.2862E-04, .1775E-06,
     1 .1399E+02,-.3238E+00, .3994E-02,-.4211E-04, .2608E-06,
     2 .4281E+01,-.9904E-01, .1234E-02,-.1305E-04, .8016E-07,
     2 .6346E+01,-.1463E+00, .1800E-02,-.1913E-04, .1202E-06,
     2 .9353E+01,-.2154E+00, .2646E-02,-.2801E-04, .1754E-06,
     2 .1375E+02,-.3143E+00, .3878E-02,-.4270E-04, .2755E-06 /
      DATA ((( CL5CFGH(K,J,I), K = 1, 5), J = 1, 4), I = 1, 2) /
     1 .2394E-04,-.9995E-07, .3865E-09,-.1762E-11, .6796E-14,
     1 .3556E-04,-.1490E-06, .5572E-09,-.2571E-11, .1282E-13,
     1 .5287E-04,-.2219E-06, .8228E-09,-.3875E-11, .2042E-13,
     1 .7859E-04,-.3299E-06, .1299E-08,-.5788E-11, .1850E-13,
     2 .2561E-04,-.1053E-06, .3964E-09,-.1842E-11, .6970E-14,
     2 .3804E-04,-.1560E-06, .5589E-09,-.2537E-11, .1297E-13,
     2 .5644E-04,-.2318E-06, .8781E-09,-.3557E-11, .8562E-14,
     2 .8398E-04,-.3427E-06, .1217E-08,-.5530E-11, .2224E-13 /
C
C----------------------------------------------------------------------C
C     SPECTRAL REGION: 800-980 CM^-1. GAS: H2O                         C
C----------------------------------------------------------------------C
C
      COMMON /BANDL6/ GWL6(3), CL6H2O(5,18,3), CL6F11(5,18,3),
     1                         CL6F12(5,18,3), CL6CS(5,4),
     2                         CL6CF(5,4)
      DATA GWL6 / .69000, .26500, .04500 /
      DATA ((( CL6H2O(K,J,I), K = 1, 5), J = 1, 18), I = 1, 3) /
     1 .2034E-05, .5496E-07, .5971E-09, .3049E-11, .5524E-14,
     1 .3018E-05, .8157E-07, .8861E-09, .4526E-11, .8197E-14,
     1 .4478E-05, .1210E-06, .1315E-08, .6716E-11, .1217E-13,
     1 .6646E-05, .1796E-06, .1951E-08, .9964E-11, .1806E-13,
     1 .9864E-05, .2666E-06, .2896E-08, .1479E-10, .2687E-13,
     1 .1463E-04, .3957E-06, .4297E-08, .2194E-10, .3991E-13,
     1 .2171E-04, .5872E-06, .6376E-08, .3255E-10, .5925E-13,
     1 .3222E-04, .8714E-06, .9464E-08, .4829E-10, .8792E-13,
     1 .4780E-04, .1293E-05, .1404E-07, .7164E-10, .1305E-12,
     1 .7092E-04, .1919E-05, .2085E-07, .1063E-09, .1929E-12,
     1 .1051E-03, .2847E-05, .3096E-07, .1577E-09, .2825E-12,
     1 .1558E-03, .4222E-05, .4597E-07, .2340E-09, .4123E-12,
     1 .2304E-03, .6253E-05, .6819E-07, .3478E-09, .6146E-12,
     1 .3397E-03, .9236E-05, .1010E-06, .5169E-09, .9155E-12,
     1 .4988E-03, .1357E-04, .1486E-06, .7680E-09, .1434E-11,
     1 .7275E-03, .1977E-04, .2160E-06, .1137E-08, .2385E-11,
     1 .1051E-02, .2856E-04, .3124E-06, .1635E-08, .3326E-11,
     1 .1499E-02, .4073E-04, .4458E-06, .2336E-08, .4737E-11,
     2 .6224E-04, .1760E-05, .2012E-07, .1084E-09, .2104E-12,
     2 .9204E-04, .2605E-05, .2987E-07, .1609E-09, .3066E-12,
     2 .1361E-03, .3858E-05, .4431E-07, .2388E-09, .4499E-12,
     2 .2013E-03, .5710E-05, .6562E-07, .3547E-09, .6832E-12,
     2 .2977E-03, .8454E-05, .9734E-07, .5264E-09, .1003E-11,
     2 .4397E-03, .1251E-04, .1443E-06, .7814E-09, .1473E-11,
     2 .6485E-03, .1845E-04, .2135E-06, .1164E-08, .2243E-11,
     2 .9544E-03, .2719E-04, .3158E-06, .1734E-08, .3350E-11,
     2 .1400E-02, .3989E-04, .4645E-06, .2574E-08, .5141E-11,
     2 .2041E-02, .5804E-04, .6752E-06, .3809E-08, .8410E-11,
     2 .2949E-02, .8364E-04, .9725E-06, .5532E-08, .1268E-10,
     2 .4221E-02, .1191E-03, .1379E-05, .7881E-08, .1873E-10,
     2 .5988E-02, .1679E-03, .1924E-05, .1100E-07, .2742E-10,
     2 .8412E-02, .2350E-03, .2670E-05, .1506E-07, .3740E-10,
     2 .1165E-01, .3245E-03, .3683E-05, .2055E-07, .4842E-10,
     2 .1586E-01, .4391E-03, .4965E-05, .2770E-07, .6569E-10,
     2 .2116E-01, .5844E-03, .6673E-05, .3649E-07, .6974E-10,
     2 .2767E-01, .7531E-03, .8515E-05, .4688E-07, .9368E-10,
     3 .7985E+00, .2077E-01, .2296E-03, .1338E-05, .3583E-08,
     3 .7979E+00, .2076E-01, .2296E-03, .1338E-05, .3580E-08,
     3 .7972E+00, .2075E-01, .2295E-03, .1337E-05, .3576E-08,
     3 .7962E+00, .2073E-01, .2294E-03, .1337E-05, .3573E-08,
     3 .7949E+00, .2071E-01, .2293E-03, .1337E-05, .3568E-08,
     3 .7932E+00, .2067E-01, .2290E-03, .1335E-05, .3561E-08,
     3 .7911E+00, .2063E-01, .2286E-03, .1333E-05, .3553E-08,
     3 .7885E+00, .2056E-01, .2279E-03, .1329E-05, .3544E-08,
     3 .7852E+00, .2047E-01, .2269E-03, .1324E-05, .3530E-08,
     3 .7805E+00, .2034E-01, .2255E-03, .1315E-05, .3509E-08,
     3 .7743E+00, .2018E-01, .2235E-03, .1304E-05, .3480E-08,
     3 .7658E+00, .1993E-01, .2208E-03, .1288E-05, .3440E-08,
     3 .7539E+00, .1960E-01, .2171E-03, .1266E-05, .3384E-08,
     3 .7376E+00, .1915E-01, .2120E-03, .1238E-05, .3317E-08,
     3 .7158E+00, .1854E-01, .2050E-03, .1201E-05, .3240E-08,
     3 .6872E+00, .1775E-01, .1961E-03, .1151E-05, .3120E-08,
     3 .6506E+00, .1674E-01, .1842E-03, .1090E-05, .3080E-08,
     3 .6052E+00, .1553E-01, .1709E-03, .1015E-05, .2911E-08 /
      DATA ((( CL6F11(K,J,I), K = 1, 5), J = 1, 18), I = 1, 3) /
     1 .3381E+04,-.2732E+01,-.1639E+00, .6997E-03, .4222E-05,
     1 .3231E+04,-.1926E+01,-.1613E+00, .6339E-03, .4676E-05,
     1 .3033E+04,-.1040E+01,-.1545E+00, .5682E-03, .4642E-05,
     1 .2786E+04,-.7138E-01,-.1450E+00, .4882E-03, .4578E-05,
     1 .2506E+04, .8395E+00,-.1358E+00, .4018E-03, .4926E-05,
     1 .2212E+04, .1500E+01,-.1207E+00, .3216E-03, .4524E-05,
     1 .1953E+04, .1769E+01,-.1060E+00, .2625E-03, .4082E-05,
     1 .1766E+04, .1686E+01,-.9499E-01, .2305E-03, .3758E-05,
     1 .1670E+04, .1436E+01,-.8888E-01, .2246E-03, .3547E-05,
     1 .1638E+04, .1264E+01,-.8698E-01, .2271E-03, .3464E-05,
     1 .1632E+04, .1214E+01,-.8635E-01, .2286E-03, .3371E-05,
     1 .1632E+04, .1184E+01,-.8583E-01, .2330E-03, .3256E-05,
     1 .1633E+04, .1164E+01,-.8623E-01, .2375E-03, .3262E-05,
     1 .1634E+04, .1253E+01,-.8736E-01, .2210E-03, .3362E-05,
     1 .1633E+04, .1310E+01,-.8648E-01, .2097E-03, .3159E-05,
     1 .1628E+04, .1227E+01,-.8157E-01, .2246E-03, .2256E-05,
     1 .1631E+04, .1091E+01,-.8327E-01, .2157E-03, .2212E-05,
     1 .1633E+04, .1296E+01,-.9664E-01, .1599E-03, .4437E-05,
     2 .6106E+04,-.1961E+02,-.2294E+00, .9638E-03, .6765E-05,
     2 .5837E+04,-.1801E+02,-.2154E+00, .9796E-03, .4416E-05,
     2 .5476E+04,-.1612E+02,-.1988E+00, .1047E-02, .1772E-05,
     2 .5028E+04,-.1338E+02,-.1840E+00, .9520E-03, .4294E-06,
     2 .4503E+04,-.1058E+02,-.1536E+00, .8490E-03,-.2598E-05,
     2 .3981E+04,-.7747E+01,-.1482E+00, .6647E-03,-.5842E-06,
     2 .3519E+04,-.6007E+01,-.1327E+00, .5541E-03,-.2578E-06,
     2 .3197E+04,-.5197E+01,-.1286E+00, .4610E-03, .1385E-05,
     2 .3028E+04,-.5128E+01,-.1193E+00, .4272E-03, .9513E-06,
     2 .2971E+04,-.5395E+01,-.1172E+00, .4522E-03, .9949E-06,
     2 .2957E+04,-.5178E+01,-.1101E+00, .4019E-03,-.2902E-06,
     2 .2950E+04,-.5062E+01,-.1120E+00, .3947E-03, .2682E-06,
     2 .2941E+04,-.4941E+01,-.1135E+00, .3850E-03, .8468E-06,
     2 .2930E+04,-.4965E+01,-.1080E+00, .3827E-03, .1011E-06,
     2 .2903E+04,-.4839E+01,-.1021E+00, .3474E-03,-.1298E-06,
     2 .2884E+04,-.4852E+01,-.1131E+00, .3768E-03, .2581E-05,
     2 .2849E+04,-.4552E+01,-.9886E-01, .4233E-03, .1158E-05,
     2 .2820E+04,-.4794E+01,-.7237E-01, .5268E-03,-.2885E-05,
     3 .6794E+04,-.2563E+02,-.8233E-01, .2736E-02,-.3357E-04,
     3 .6490E+04,-.2222E+02,-.1202E+00, .2399E-02,-.2502E-04,
     3 .6124E+04,-.1574E+02,-.1942E+00, .1330E-02,-.7309E-05,
     3 .5662E+04,-.1167E+02,-.2352E+00, .1043E-02, .2647E-05,
     3 .5138E+04,-.7515E+01,-.2862E+00, .6414E-03, .1440E-04,
     3 .4524E+04,-.6262E+01,-.2168E+00, .7032E-03, .5112E-05,
     3 .3987E+04,-.4657E+01,-.1870E+00, .5214E-03, .3763E-05,
     3 .3553E+04,-.4860E+01,-.9823E-01, .6298E-03,-.7960E-05,
     3 .3332E+04,-.5270E+01,-.8622E-01, .6832E-03,-.7269E-05,
     3 .3269E+04,-.4862E+01,-.7506E-01, .5417E-03,-.8338E-05,
     3 .3284E+04,-.6640E+01,-.1195E+00, .8530E-03, .3108E-06,
     3 .3312E+04,-.6994E+01,-.1160E+00, .8318E-03,-.1221E-05,
     3 .3348E+04,-.7408E+01,-.1011E+00, .8199E-03,-.4718E-05,
     3 .3396E+04,-.8623E+01,-.1163E+00, .1087E-02,-.1853E-05,
     3 .3569E+04,-.1024E+02,-.1645E+00, .1468E-02, .2608E-05,
     3 .3743E+04,-.8904E+01,-.1749E+00, .1068E-02, .4947E-06,
     3 .3915E+04,-.8567E+01,-.2328E+00, .9282E-03, .9550E-05,
     3 .4051E+04,-.1028E+02,-.1838E+00, .1174E-02,-.7296E-06 /
      DATA ((( CL6F12(K,J,I), K = 1, 5), J = 1, 18), I = 1, 3) /
     1 .3831E+04,-.6616E+01,-.1311E+00, .3745E-03, .5462E-05,
     1 .3690E+04,-.5755E+01,-.1296E+00, .3549E-03, .5348E-05,
     1 .3501E+04,-.4707E+01,-.1258E+00, .3331E-03, .5043E-05,
     1 .3264E+04,-.3453E+01,-.1212E+00, .2984E-03, .4772E-05,
     1 .2987E+04,-.2156E+01,-.1142E+00, .2519E-03, .4460E-05,
     1 .2694E+04,-.1040E+01,-.1052E+00, .2027E-03, .4145E-05,
     1 .2430E+04,-.3488E+00,-.9475E-01, .1635E-03, .3724E-05,
     1 .2238E+04,-.1245E+00,-.8580E-01, .1381E-03, .3398E-05,
     1 .2137E+04,-.2169E+00,-.8068E-01, .1339E-03, .3191E-05,
     1 .2104E+04,-.3494E+00,-.7895E-01, .1348E-03, .3156E-05,
     1 .2099E+04,-.3710E+00,-.7997E-01, .1299E-03, .3388E-05,
     1 .2099E+04,-.3345E+00,-.8031E-01, .1120E-03, .3593E-05,
     1 .2100E+04,-.3338E+00,-.7947E-01, .1153E-03, .3395E-05,
     1 .2103E+04,-.4294E+00,-.7784E-01, .1408E-03, .2912E-05,
     1 .2107E+04,-.6355E+00,-.7374E-01, .1914E-03, .1978E-05,
     1 .2124E+04,-.6312E+00,-.8343E-01, .1962E-03, .3349E-05,
     1 .2142E+04,-.5083E+00,-.9308E-01, .1848E-03, .4644E-05,
     1 .2148E+04,-.5107E+00,-.8249E-01, .1974E-03, .2445E-05,
     2 .3174E+04, .1559E+01,-.1525E+00, .9536E-03,-.5138E-05,
     2 .3059E+04, .2038E+01,-.1502E+00, .8782E-03,-.4550E-05,
     2 .2905E+04, .2563E+01,-.1422E+00, .8081E-03,-.4707E-05,
     2 .2707E+04, .3185E+01,-.1359E+00, .6975E-03,-.3539E-05,
     2 .2474E+04, .3716E+01,-.1232E+00, .6017E-03,-.3297E-05,
     2 .2234E+04, .4198E+01,-.1148E+00, .4547E-03,-.1784E-05,
     2 .2015E+04, .4453E+01,-.1056E+00, .3494E-03,-.9578E-06,
     2 .1855E+04, .4296E+01,-.9284E-01, .3081E-03,-.1551E-05,
     2 .1775E+04, .3852E+01,-.8603E-01, .3381E-03,-.1903E-05,
     2 .1749E+04, .3716E+01,-.8444E-01, .3368E-03,-.2097E-05,
     2 .1749E+04, .3599E+01,-.8786E-01, .3613E-03,-.1568E-05,
     2 .1750E+04, .3597E+01,-.8534E-01, .3824E-03,-.2436E-05,
     2 .1749E+04, .3518E+01,-.8555E-01, .3936E-03,-.2192E-05,
     2 .1742E+04, .3750E+01,-.8911E-01, .3302E-03,-.1102E-05,
     2 .1740E+04, .4170E+01,-.1042E+00, .2203E-03, .1979E-05,
     2 .1700E+04, .4231E+01,-.7888E-01, .1958E-03,-.1699E-05,
     2 .1659E+04, .4097E+01,-.5662E-01, .2025E-03,-.4840E-05,
     2 .1651E+04, .4286E+01,-.8666E-01, .1600E-03, .1194E-05,
     3 .2697E+04, .1630E+01,-.1335E+00, .2269E-03, .3152E-05,
     3 .2583E+04, .2012E+01,-.1201E+00, .2585E-03, .2430E-05,
     3 .2430E+04, .3323E+01,-.1422E+00, .4445E-04, .8851E-05,
     3 .2263E+04, .3939E+01,-.1253E+00, .5915E-05, .6341E-05,
     3 .2072E+04, .4716E+01,-.1336E+00,-.1317E-03, .9209E-05,
     3 .1847E+04, .4478E+01,-.9821E-01, .7760E-04, .2975E-05,
     3 .1656E+04, .3802E+01,-.7643E-01, .1790E-03, .7571E-06,
     3 .1527E+04, .3486E+01,-.9055E-01, .1987E-03, .4812E-05,
     3 .1435E+04, .4339E+01,-.9408E-01,-.5978E-04, .6688E-05,
     3 .1396E+04, .4102E+01,-.9198E-01,-.2017E-04, .7000E-05,
     3 .1346E+04, .4119E+01,-.5129E-01,-.6117E-04, .1038E-06,
     3 .1332E+04, .3472E+01,-.6103E-01, .9062E-04, .2054E-05,
     3 .1331E+04, .3920E+01,-.7302E-01,-.2448E-04, .3725E-05,
     3 .1320E+04, .4012E+01,-.7672E-01,-.4056E-04, .4657E-05,
     3 .1259E+04, .4705E+01,-.5093E-01,-.1710E-03, .8423E-06,
     3 .1247E+04, .4276E+01,-.5120E-01,-.1004E-03, .1440E-05,
     3 .1211E+04, .3173E+01,-.3464E-01, .3911E-04, .1399E-06,
     3 .1187E+04, .2100E+01,-.1986E-01, .9420E-04,-.1724E-05 /
      DATA (( CL6CS(K,J), K = 1, 5), J = 1, 4)               /
     1 .1365E+02,-.3302E+00, .4251E-02,-.4696E-04, .2990E-06,
     1 .2023E+02,-.4893E+00, .6302E-02,-.6968E-04, .4436E-06,
     1 .2998E+02,-.7237E+00, .9309E-02,-.1037E-03, .6659E-06,
     1 .4437E+02,-.1072E+01, .1384E-01,-.1541E-03, .9830E-06 /
      DATA (( CL6CF(K,J), K = 1, 5), J = 1, 4)               /
     1 .4005E-03,-.1364E-05, .4014E-09,-.1065E-10, .3771E-12,
     1 .5898E-03,-.1974E-05, .1180E-08,-.2475E-10, .4965E-12,
     1 .8665E-03,-.2779E-05, .2611E-09,-.6482E-10, .1219E-11,
     1 .1264E-02,-.4046E-05, .7903E-08,-.1068E-09, .8824E-12 /
C
C----------------------------------------------------------------------C
C     SPECTRAL REGION: 540-800 CM^-1. GAS: H2O, CO2, N2O               C
C----------------------------------------------------------------------C
C
      COMMON /BANDL7/ GW(3), CL7H2OU(5,11,3), CL7H2OD(5,5,7,3),
     1                       CL7CO2U(5,11,3), CL7CO2D(5,5,7,3),
     2                       CL7N2O(3), CL7O3(3),
     3                       CL7CS(5,5,4,3), CL7CF(5,5,4,3)
      DATA GW / .14000, .23500, .24000 /
      DATA CL7N2O / 0.8775E+01, 0.6730E+02, 0.6931E+02 /
      DATA CL7O3 / 6.7342, 10.9696, 42.1065 /
      DATA ((( CL7H2OU(K,J,I), K = 1, 5), J = 1, 11), I = 1, 3) /
     1 .3198E-03, .7010E-05, .5087E-07, .1823E-09, .8000E-12,
     1 .4744E-03, .1039E-04, .7481E-07, .2733E-09, .1358E-11,
     1 .7031E-03, .1544E-04, .1106E-06, .4005E-09, .2096E-11,
     1 .1043E-02, .2301E-04, .1664E-06, .5895E-09, .2823E-11,
     1 .1542E-02, .3408E-04, .2504E-06, .8968E-09, .3832E-11,
     1 .2277E-02, .5052E-04, .3783E-06, .1336E-08, .4584E-11,
     1 .3349E-02, .7481E-04, .5774E-06, .1984E-08, .4136E-11,
     1 .4915E-02, .1107E-03, .8788E-06, .2888E-08, .2030E-11,
     1 .7128E-02, .1599E-03, .1302E-05, .4852E-08, .4836E-11,
     1 .1027E-01, .2335E-03, .2002E-05, .6551E-08,-.1648E-10,
     1 .1473E-01, .3365E-03, .2992E-05, .9252E-08,-.4561E-10,
     2 .3860E-02, .6592E-04, .4173E-06, .1467E-08, .2925E-11,
     2 .5671E-02, .9710E-04, .6148E-06, .1939E-08, .1773E-11,
     2 .8281E-02, .1405E-03, .9136E-06, .3551E-08, .6629E-11,
     2 .1224E-01, .2051E-03, .1275E-05, .5589E-08, .2230E-10,
     2 .1786E-01, .2895E-03, .1683E-05, .9969E-08, .6822E-10,
     2 .2603E-01, .4244E-03, .2602E-05, .1281E-07, .5263E-10,
     2 .3822E-01, .6217E-03, .3781E-05, .1670E-07, .4943E-10,
     2 .5569E-01, .8796E-03, .4931E-05, .2712E-07, .1516E-09,
     2 .8122E-01, .1300E-02, .7901E-05, .2939E-07,-.1662E-10,
     2 .1186E+00, .1870E-02, .1121E-04, .3777E-07,-.1080E-09,
     2 .1744E+00, .2704E-02, .1413E-04, .4233E-07,-.2056E-10,
     3 .2535E-01, .4376E-03, .2415E-05, .5321E-08, .2793E-10,
     3 .3648E-01, .6177E-03, .3550E-05, .7499E-08,-.8573E-11,
     3 .5284E-01, .8857E-03, .4607E-05, .7056E-08, .1293E-10,
     3 .7569E-01, .1249E-02, .5929E-05, .1228E-07, .1470E-09,
     3 .1098E+00, .1771E-02, .7478E-05, .2647E-07, .4460E-09,
     3 .1581E+00, .2472E-02, .9066E-05, .4093E-07, .8467E-09,
     3 .2314E+00, .3685E-02, .1369E-04, .1789E-07, .7564E-09,
     3 .3308E+00, .5356E-02, .2346E-04,-.5571E-08, .1889E-09,
     3 .4624E+00, .7352E-02, .3650E-04, .3410E-07,-.1585E-09,
     3 .6519E+00, .1004E-01, .3635E-04, .9474E-07, .2475E-08,
     3 .9176E+00, .1456E-01, .6342E-04,-.2006E-08, .3091E-09 /
      DATA (((( CL7H2OD(K,L,J,I), K = 1, 5), L = 1, 5), J = 1, 7),
     1                                                  I = 1, 2) /
     1 .2122E-01, .4754E-03, .3935E-05, .1405E-07,-.7501E-11,
     1 .1206E-01, .2866E-03, .2068E-05, .3657E-08, .1041E-10,
     1 .6465E-02, .1665E-03, .1387E-05, .1798E-08,-.2136E-10,
     1 .3193E-02, .8096E-04, .7273E-06, .2920E-08, .8332E-11,
     1 .1888E-02, .4747E-04, .3982E-06,-.1080E-09,-.1454E-10,
     1 .3023E-01, .6682E-03, .5530E-05, .2168E-07, .8249E-12,
     1 .1711E-01, .4164E-03, .3124E-05, .3919E-08,-.1459E-10,
     1 .9017E-02, .2384E-03, .2100E-05, .2122E-08,-.5086E-10,
     1 .4516E-02, .1140E-03, .1005E-05, .4126E-08, .1567E-10,
     1 .2661E-02, .6773E-04, .6267E-06, .1011E-09,-.2749E-10,
     1 .4278E-01, .9224E-03, .7342E-05, .3254E-07, .6124E-10,
     1 .2396E-01, .5895E-03, .4432E-05, .5729E-08,-.9743E-11,
     1 .1242E-01, .3331E-03, .3106E-05, .4000E-08,-.8491E-10,
     1 .6415E-02, .1521E-03, .1074E-05, .6962E-08, .8110E-10,
     1 .3739E-02, .9653E-04, .9598E-06, .2371E-09,-.4901E-10,
     1 .6120E-01, .1324E-02, .9951E-05, .3875E-07, .9646E-10,
     1 .3376E-01, .8242E-03, .6056E-05, .9813E-08, .2581E-10,
     1 .1695E-01, .4537E-03, .4180E-05, .6927E-08,-.8408E-10,
     1 .8837E-02, .2112E-03, .1434E-05, .8920E-08, .1169E-09,
     1 .5242E-02, .1346E-03, .1358E-05, .9875E-09,-.6359E-10,
     1 .8754E-01, .1892E-02, .1266E-04, .4015E-07, .2227E-09,
     1 .4642E-01, .1137E-02, .8871E-05, .1865E-07, .1096E-11,
     1 .2276E-01, .6006E-03, .5328E-05, .1010E-07,-.6748E-10,
     1 .1207E-01, .3052E-03, .2477E-05, .7695E-08, .3896E-10,
     1 .7335E-02, .1780E-03, .1611E-05, .3203E-08,-.3384E-10,
     1 .1204E+00, .2525E-02, .1520E-04, .4697E-07, .4019E-09,
     1 .6313E-01, .1528E-02, .1087E-04, .2288E-07, .1257E-09,
     1 .2950E-01, .8041E-03, .7944E-05, .1390E-07,-.2084E-09,
     1 .1603E-01, .3958E-03, .3219E-05, .1115E-07, .5508E-10,
     1 .1018E-01, .2421E-03, .2095E-05, .5623E-08,-.1446E-10,
     1 .1625E+00, .3545E-02, .2246E-04, .1677E-07,-.1586E-09,
     1 .8242E-01, .2086E-02, .1565E-04, .1280E-07,-.1366E-09,
     1 .3755E-01, .9723E-03, .8996E-05, .3195E-07, .2688E-10,
     1 .2144E-01, .5492E-03, .5077E-05, .9584E-08,-.9889E-10,
     1 .1378E-01, .3173E-03, .2807E-05, .1019E-07,-.4310E-11,
     2 .2502E+00, .3828E-02, .2098E-04, .6728E-07,-.2163E-09,
     2 .1351E+00, .2202E-02, .1169E-04, .2027E-07,-.3301E-10,
     2 .7174E-01, .1160E-02, .5519E-05, .1603E-07, .1577E-09,
     2 .3879E-01, .6627E-03, .4357E-05, .1218E-07,-.3509E-10,
     2 .2402E-01, .4150E-03, .2207E-05, .5301E-08, .2091E-10,
     2 .3564E+00, .5352E-02, .2909E-04, .1161E-06,-.1663E-09,
     2 .1949E+00, .3118E-02, .1478E-04, .2446E-07, .1307E-09,
     2 .1038E+00, .1665E-02, .7581E-05, .1688E-07, .1882E-09,
     2 .5596E-01, .9185E-03, .5313E-05, .2023E-07, .7491E-10,
     2 .3445E-01, .5925E-03, .3029E-05, .6692E-08, .3801E-10,
     2 .5083E+00, .7692E-02, .3604E-04, .1088E-06, .1514E-09,
     2 .2792E+00, .4398E-02, .2038E-04, .2957E-07, .4913E-10,
     2 .1498E+00, .2363E-02, .1009E-04, .2250E-07, .2980E-09,
     2 .8048E-01, .1288E-02, .6713E-05, .2666E-07, .1613E-09,
     2 .4856E-01, .8429E-03, .4527E-05, .6941E-08,-.5864E-11,
     2 .7258E+00, .1116E-01, .4425E-04, .7730E-07, .7876E-09,
     2 .3979E+00, .6172E-02, .2732E-04, .4212E-07, .5499E-10,
     2 .2161E+00, .3366E-02, .1360E-04, .2253E-07, .3647E-09,
     2 .1146E+00, .1853E-02, .9905E-05, .2338E-07, .2420E-10,
     2 .6742E-01, .1155E-02, .6134E-05, .1191E-07, .2888E-10,
     2 .1030E+01, .1615E-01, .5639E-04, .1022E-06, .2680E-08,
     2 .5701E+00, .9025E-02, .3124E-04,-.5345E-08, .1152E-08,
     2 .3057E+00, .4857E-02, .2139E-04, .9328E-08, .1029E-09,
     2 .1600E+00, .2521E-02, .1453E-04, .4003E-07,-.1410E-09,
     2 .9516E-01, .1618E-02, .7903E-05, .7070E-08, .2683E-10,
     2 .1405E+01, .2295E-01, .9099E-04,-.2412E-07, .1033E-08,
     2 .7993E+00, .1288E-01, .4976E-04,-.1862E-07, .1021E-08,
     2 .4242E+00, .6613E-02, .2929E-04, .4863E-07, .4837E-09,
     2 .2273E+00, .3429E-02, .1564E-04, .4485E-07, .1848E-09,
     2 .1359E+00, .2303E-02, .1000E-04, .2450E-08, .1606E-09,
     2 .1910E+01, .3128E-01, .1092E-03, .5854E-07, .6144E-08,
     2 .1107E+01, .1746E-01, .6235E-04, .5236E-07, .2686E-08,
     2 .5959E+00, .9343E-02, .3786E-04,-.3029E-07, .1286E-09,
     2 .3203E+00, .4732E-02, .1825E-04, .3272E-07, .2774E-09,
     2 .1935E+00, .3265E-02, .1215E-04,-.1320E-07, .3797E-09 /
      DATA (((( CL7H2OD(K,L,J,I), K = 1, 5), L = 1, 5), J = 1, 7),
     1                                                  I = 3, 3) /
     3 .1288E+01, .1982E-01, .7932E-04, .1506E-06, .2769E-08,
     3 .6732E+00, .1075E-01, .5420E-04,-.1191E-07,-.1135E-08,
     3 .3483E+00, .5282E-02, .2187E-04, .1975E-07, .2887E-09,
     3 .1959E+00, .3065E-02, .1492E-04,-.1861E-08,-.3155E-09,
     3 .1173E+00, .1856E-02, .8988E-05,-.1027E-07,-.2508E-09,
     3 .1793E+01, .2778E-01, .1076E-03, .2464E-07, .2382E-08,
     3 .9374E+00, .1489E-01, .7284E-04,-.2248E-07,-.1285E-08,
     3 .4796E+00, .7298E-02, .3483E-04, .1079E-07,-.6009E-09,
     3 .2739E+00, .4303E-02, .2221E-04,-.6108E-08,-.6761E-09,
     3 .1680E+00, .2540E-02, .1009E-04, .5848E-08, .1792E-09,
     3 .2479E+01, .3818E-01, .1255E-03,-.9731E-07, .6012E-08,
     3 .1302E+01, .1940E-01, .7569E-04, .1666E-06, .2569E-08,
     3 .6647E+00, .9983E-02, .4455E-04, .6638E-08,-.4726E-09,
     3 .3882E+00, .5991E-02, .2640E-04,-.2330E-07,-.4318E-09,
     3 .2409E+00, .3575E-02, .1264E-04, .2199E-07, .6260E-09,
     3 .3448E+01, .5463E-01, .2772E-03,-.6672E-06,-.1532E-07,
     3 .1791E+01, .2708E-01, .1157E-03, .4916E-07, .2512E-09,
     3 .9230E+00, .1424E-01, .7422E-04,-.7568E-07,-.3182E-08,
     3 .5445E+00, .8268E-02, .3229E-04, .3190E-08, .4198E-09,
     3 .3417E+00, .5074E-02, .1897E-04, .2249E-07, .5607E-09,
     3 .4601E+01, .6643E-01, .3080E-03, .1927E-06,-.4891E-08,
     3 .2428E+01, .3692E-01, .1684E-03, .1351E-06, .2882E-09,
     3 .1278E+01, .1882E-01, .8434E-04, .1643E-06, .3421E-09,
     3 .7541E+00, .1167E-01, .5045E-04,-.3429E-07,-.4394E-09,
     3 .4761E+00, .7213E-02, .3041E-04,-.8801E-08,-.1827E-09,
     3 .6015E+01, .8087E-01, .3583E-03, .1619E-05, .7963E-08,
     3 .3308E+01, .4919E-01, .2322E-03,-.2056E-07,-.6857E-08,
     3 .1731E+01, .2596E-01, .1057E-03, .7529E-08, .5775E-09,
     3 .1015E+01, .1601E-01, .8324E-04,-.2715E-07,-.2508E-08,
     3 .6430E+00, .9946E-02, .4431E-04,-.4379E-07,-.8073E-09,
     3 .7170E+01, .9746E-01, .5182E-03,-.5037E-06,-.3281E-07,
     3 .4311E+01, .6002E-01, .2196E-03, .1687E-06, .2666E-08,
     3 .2263E+01, .3324E-01, .1076E-03, .8376E-07, .6195E-08,
     3 .1332E+01, .2104E-01, .1020E-03, .6250E-07,-.5010E-09,
     3 .8303E+00, .1303E-01, .6471E-04,-.3176E-07,-.1833E-08 /
      DATA ((( CL7CO2U(K,J,I), K = 1, 5), J = 1, 11), I = 1, 3) /
     1 .4904E-04, .1609E-05, .2180E-07, .1390E-09, .3177E-12,
     1 .7290E-04, .2394E-05, .3248E-07, .2065E-09, .4603E-12,
     1 .1084E-03, .3562E-05, .4835E-07, .3072E-09, .6805E-12,
     1 .1610E-03, .5291E-05, .7169E-07, .4570E-09, .1049E-11,
     1 .2392E-03, .7867E-05, .1064E-06, .6785E-09, .1588E-11,
     1 .3543E-03, .1167E-04, .1582E-06, .1010E-08, .2369E-11,
     1 .5234E-03, .1727E-04, .2351E-06, .1505E-08, .3491E-11,
     1 .7692E-03, .2538E-04, .3466E-06, .2250E-08, .5443E-11,
     1 .1124E-02, .3697E-04, .5026E-06, .3334E-08, .9072E-11,
     1 .1623E-02, .5275E-04, .7100E-06, .4932E-08, .1625E-10,
     1 .2316E-02, .7453E-04, .9939E-06, .6987E-08, .2446E-10,
     2 .1502E-02, .3778E-04, .3395E-06, .1011E-08,-.1556E-11,
     2 .2199E-02, .5550E-04, .5037E-06, .1546E-08,-.2156E-11,
     2 .3214E-02, .8127E-04, .7393E-06, .2358E-08,-.2031E-11,
     2 .4677E-02, .1182E-03, .1078E-05, .3660E-08,-.2471E-12,
     2 .6791E-02, .1719E-03, .1570E-05, .5437E-08, .1264E-11,
     2 .9799E-02, .2477E-03, .2266E-05, .8090E-08, .4361E-11,
     2 .1406E-01, .3540E-03, .3248E-05, .1193E-07, .7723E-11,
     2 .2007E-01, .5033E-03, .4602E-05, .1701E-07, .1187E-10,
     2 .2856E-01, .7077E-03, .6282E-05, .2422E-07, .4679E-10,
     2 .4034E-01, .9965E-03, .8614E-05, .3165E-07, .7316E-10,
     2 .5625E-01, .1394E-02, .1197E-04, .4058E-07, .7027E-10,
     3 .2543E-01, .5185E-03, .3802E-05, .8871E-08,-.1849E-10,
     3 .3643E-01, .7330E-03, .5302E-05, .1280E-07,-.2020E-10,
     3 .5234E-01, .1039E-02, .7422E-05, .1805E-07,-.2737E-10,
     3 .7547E-01, .1478E-02, .1027E-04, .2461E-07,-.2592E-10,
     3 .1088E+00, .2103E-02, .1412E-04, .3220E-07,-.1723E-10,
     3 .1566E+00, .3000E-02, .1945E-04, .4064E-07,-.3222E-11,
     3 .2244E+00, .4274E-02, .2712E-04, .4915E-07,-.3124E-10,
     3 .3186E+00, .6048E-02, .3849E-04, .5988E-07,-.2055E-09,
     3 .4479E+00, .8441E-02, .5385E-04, .8065E-07,-.4098E-09,
     3 .6230E+00, .1159E-01, .7311E-04, .1143E-06,-.5444E-09,
     3 .8571E+00, .1568E-01, .9611E-04, .1568E-06,-.5048E-09 /
      DATA (((( CL7CO2D(K,L,J,I), K = 1, 5), L = 1, 5), J = 1, 7),
     1                                                  I = 1, 2) /
     1 .3255E-02, .1040E-03, .1387E-05, .9581E-08, .3107E-10,
     1 .4614E-02, .1308E-03, .1606E-05, .1133E-07, .3897E-10,
     1 .8000E-02, .2025E-03, .2015E-05, .1305E-07, .6382E-10,
     1 .1608E-01, .4327E-03, .4148E-05, .8855E-08,-.7839E-10,
     1 .3006E-01, .7894E-03, .7339E-05, .4311E-07, .2559E-09,
     1 .4549E-02, .1443E-03, .1895E-05, .1271E-07, .3962E-10,
     1 .6526E-02, .1809E-03, .2151E-05, .1538E-07, .5875E-10,
     1 .1145E-01, .2827E-03, .2617E-05, .1811E-07, .1158E-09,
     1 .2264E-01, .6134E-03, .5976E-05, .1229E-07,-.1286E-09,
     1 .4252E-01, .1115E-02, .9811E-05, .5431E-07, .3797E-09,
     1 .6295E-02, .1988E-03, .2587E-05, .1663E-07, .4580E-10,
     1 .9195E-02, .2505E-03, .2872E-05, .1991E-07, .7763E-10,
     1 .1628E-01, .3991E-03, .3504E-05, .2247E-07, .1547E-09,
     1 .3086E-01, .8906E-03, .1019E-04, .1209E-07,-.4987E-09,
     1 .5926E-01, .1522E-02, .1236E-04, .7739E-07, .7260E-09,
     1 .8550E-02, .2658E-03, .3428E-05, .2218E-07, .6261E-10,
     1 .1271E-01, .3464E-03, .3894E-05, .2481E-07, .8621E-10,
     1 .2267E-01, .5632E-03, .5213E-05, .2786E-07, .1182E-09,
     1 .4263E-01, .1203E-02, .1377E-04, .2406E-07,-.6019E-09,
     1 .8037E-01, .2070E-02, .1694E-04, .9654E-07, .8452E-09,
     1 .1140E-01, .3506E-03, .4551E-05, .2893E-07, .6920E-10,
     1 .1743E-01, .4719E-03, .5160E-05, .3007E-07, .9235E-10,
     1 .3140E-01, .7874E-03, .7397E-05, .3604E-07, .1089E-09,
     1 .5772E-01, .1532E-02, .1545E-04, .5578E-07,-.1105E-09,
     1 .1060E+00, .2889E-02, .2743E-04, .1007E-06, .1655E-09,
     1 .1514E-01, .4601E-03, .5967E-05, .3708E-07, .7576E-10,
     1 .2354E-01, .6282E-03, .6910E-05, .3900E-07, .9360E-10,
     1 .4317E-01, .1064E-02, .8955E-05, .4370E-07, .2462E-09,
     1 .7690E-01, .2072E-02, .2034E-04, .6934E-07,-.8692E-10,
     1 .1381E+00, .3797E-02, .3719E-04, .1078E-06,-.3219E-09,
     1 .1979E-01, .5740E-03, .7208E-05, .4868E-07, .1467E-09,
     1 .3168E-01, .8015E-03, .8039E-05, .5036E-07, .2217E-09,
     1 .5674E-01, .1496E-02, .1408E-04, .3696E-07,-.1705E-09,
     1 .9986E-01, .2547E-02, .2090E-04, .1060E-06, .8137E-09,
     1 .1824E+00, .5133E-02, .4785E-04, .9002E-07,-.5910E-09,
     2 .7750E-01, .1926E-02, .1647E-04, .5152E-07, .5789E-10,
     2 .9252E-01, .2176E-02, .1805E-04, .5697E-07, .4090E-10,
     2 .1256E+00, .2757E-02, .2200E-04, .6171E-07,-.7270E-10,
     2 .1997E+00, .3954E-02, .2617E-04, .6572E-07, .1797E-09,
     2 .3415E+00, .6298E-02, .4389E-04, .1357E-06, .8247E-10,
     2 .1051E+00, .2611E-02, .2185E-04, .6654E-07, .1234E-09,
     2 .1264E+00, .2955E-02, .2427E-04, .7624E-07, .6777E-10,
     2 .1737E+00, .3761E-02, .2920E-04, .8313E-07,-.2525E-10,
     2 .2793E+00, .5549E-02, .3784E-04, .6902E-07,-.2349E-09,
     2 .4851E+00, .8731E-02, .5836E-04, .1805E-06, .1464E-09,
     2 .1393E+00, .3468E-02, .2892E-04, .8188E-07, .1027E-09,
     2 .1688E+00, .3947E-02, .3179E-04, .9314E-07, .8401E-10,
     2 .2355E+00, .5063E-02, .3872E-04, .1004E-06,-.1287E-09,
     2 .3871E+00, .7623E-02, .5141E-04, .8296E-07,-.4458E-09,
     2 .6934E+00, .1205E-01, .7024E-04, .2347E-06, .1293E-08,
     2 .1789E+00, .4491E-02, .3824E-04, .9838E-07,-.7182E-10,
     2 .2207E+00, .5185E-02, .4150E-04, .1067E-06,-.1754E-10,
     2 .3149E+00, .6732E-02, .5003E-04, .1184E-06,-.1814E-09,
     2 .5358E+00, .1028E-01, .6192E-04, .1167E-06, .4598E-09,
     2 .9829E+00, .1710E-01, .9858E-04, .2107E-06, .2801E-09,
     2 .2253E+00, .5681E-02, .4891E-04, .1097E-06,-.3620E-09,
     2 .2837E+00, .6655E-02, .5329E-04, .1255E-06,-.2267E-09,
     2 .4190E+00, .8860E-02, .6180E-04, .1329E-06, .3621E-10,
     2 .7333E+00, .1419E-01, .8173E-04, .7045E-07, .3811E-09,
     2 .1347E+01, .2269E-01, .1314E-03, .3965E-06, .8550E-09,
     2 .2805E+00, .7112E-02, .6151E-04, .1254E-06,-.6195E-09,
     2 .3568E+00, .8458E-02, .6798E-04, .1356E-06,-.5463E-09,
     2 .5474E+00, .1172E-01, .8464E-04, .1108E-06,-.1226E-08,
     2 .9699E+00, .1876E-01, .1179E-03, .1239E-06,-.7090E-09,
     2 .1815E+01, .2923E-01, .1676E-03, .5291E-06,-.3456E-09,
     2 .3383E+00, .8644E-02, .7789E-04, .1622E-06,-.1140E-08,
     2 .4362E+00, .1044E-01, .8731E-04, .1638E-06,-.1184E-08,
     2 .6993E+00, .1458E-01, .1057E-03, .2338E-06,-.7636E-09,
     2 .1276E+01, .2441E-01, .1581E-03, .1820E-06,-.1661E-08,
     2 .2438E+01, .3783E-01, .2135E-03, .7464E-06,-.1470E-08 /
      DATA (((( CL7CO2D(K,L,J,I), K = 1, 5), L = 1, 5), J = 1, 7),
     1                                                  I = 3, 3) /
     3 .1167E+01, .2108E-01, .1237E-03, .1746E-06,-.5044E-09,
     3 .1275E+01, .2281E-01, .1314E-03, .1947E-06,-.2232E-09,
     3 .1486E+01, .2643E-01, .1517E-03, .1999E-06,-.5959E-09,
     3 .1876E+01, .3241E-01, .1728E-03, .2616E-06, .1049E-08,
     3 .2672E+01, .4497E-01, .2296E-03, .3947E-06, .2269E-08,
     3 .1565E+01, .2793E-01, .1598E-03, .1901E-06,-.6590E-09,
     3 .1712E+01, .3031E-01, .1698E-03, .2092E-06,-.3258E-09,
     3 .2002E+01, .3520E-01, .1958E-03, .2216E-06,-.5193E-09,
     3 .2531E+01, .4308E-01, .2237E-03, .3119E-06, .1414E-08,
     3 .3605E+01, .6133E-01, .3406E-03, .2296E-06,-.4799E-08,
     3 .2059E+01, .3596E-01, .2083E-03, .2687E-06,-.1859E-08,
     3 .2258E+01, .3921E-01, .2210E-03, .2619E-06,-.1513E-08,
     3 .2658E+01, .4551E-01, .2468E-03, .3205E-06,-.4218E-09,
     3 .3371E+01, .5604E-01, .2877E-03, .4227E-06, .1310E-08,
     3 .4826E+01, .8129E-01, .4636E-03,-.4164E-07,-.1377E-07,
     3 .2713E+01, .4529E-01, .2301E-03, .3232E-06, .9465E-09,
     3 .2984E+01, .4978E-01, .2522E-03, .2712E-06, .4486E-10,
     3 .3523E+01, .5791E-01, .2817E-03, .3919E-06, .2269E-08,
     3 .4482E+01, .7323E-01, .3860E-03, .2529E-06,-.5461E-08,
     3 .6493E+01, .1068E+00, .5665E-03,-.1502E-06,-.1384E-07,
     3 .3571E+01, .5797E-01, .2707E-03, .1700E-06, .1403E-09,
     3 .3929E+01, .6345E-01, .3014E-03, .1354E-06,-.1783E-08,
     3 .4642E+01, .7465E-01, .3501E-03, .1434E-06,-.1217E-08,
     3 .5959E+01, .9345E-01, .4426E-03, .4426E-06,-.1279E-08,
     3 .8688E+01, .1389E+00, .6834E-03, .9559E-07,-.8553E-08,
     3 .4648E+01, .7368E-01, .3137E-03,-.4722E-07,-.3984E-10,
     3 .5107E+01, .7983E-01, .3357E-03, .7199E-07, .1513E-08,
     3 .6056E+01, .9394E-01, .3957E-03, .1422E-06, .1723E-08,
     3 .7848E+01, .1185E+00, .4563E-03, .5001E-06, .1056E-07,
     3 .1153E+02, .1781E+00, .8249E-03, .7612E-06,-.1393E-08,
     3 .5932E+01, .9222E-01, .3615E-03,-.1798E-06, .1687E-08,
     3 .6475E+01, .9998E-01, .4053E-03,-.2179E-06,-.1145E-08,
     3 .7702E+01, .1173E+00, .4919E-03,-.2205E-06,-.5892E-08,
     3 .9986E+01, .1478E+00, .5509E-03, .9207E-07, .5696E-08,
     3 .1490E+02, .2225E+00, .9115E-03, .1485E-05, .1755E-07 /
      DATA (((( CL7CS(K,L,J,I), K = 1, 5),  L = 1, 5), J = 1, 4),
     1                          I = 1, 3)                    /
     1 .1396E+02,-.2899E+00, .3004E-02,-.2032E-04, .7639E-07,
     1 .1353E+02,-.2584E+00, .2485E-02,-.2471E-04, .1733E-06,
     1 .1257E+02,-.2425E+00, .2678E-02,-.2986E-04, .1957E-06,
     1 .1220E+02,-.2489E+00, .2698E-02,-.4082E-04, .3670E-06,
     1 .1250E+02,-.2966E+00, .4165E-02,-.4513E-04, .2460E-06,
     1 .2076E+02,-.4239E+00, .4211E-02,-.2978E-04, .1392E-06,
     1 .1980E+02,-.3811E+00, .3767E-02,-.3469E-04, .2172E-06,
     1 .1832E+02,-.3536E+00, .3858E-02,-.4568E-04, .3246E-06,
     1 .1781E+02,-.3669E+00, .4217E-02,-.6134E-04, .5105E-06,
     1 .1853E+02,-.4439E+00, .6060E-02,-.6637E-04, .3895E-06,
     1 .3051E+02,-.6205E+00, .6107E-02,-.4250E-04, .1959E-06,
     1 .2883E+02,-.5561E+00, .5442E-02,-.4706E-04, .2801E-06,
     1 .2641E+02,-.5175E+00, .5890E-02,-.6928E-04, .4738E-06,
     1 .2575E+02,-.5494E+00, .6783E-02,-.9038E-04, .6733E-06,
     1 .2749E+02,-.6635E+00, .8844E-02,-.9781E-04, .6089E-06,
     1 .4463E+02,-.8951E+00, .8654E-02,-.6029E-04, .2823E-06,
     1 .4156E+02,-.7814E+00, .7533E-02,-.7295E-04, .4848E-06,
     1 .3793E+02,-.7457E+00, .8307E-02,-.1079E-03, .8356E-06,
     1 .3743E+02,-.8170E+00, .1082E-01,-.1344E-03, .8735E-06,
     1 .4082E+02,-.9914E+00, .1288E-01,-.1441E-03, .9500E-06,
     2 .3576E+02,-.7950E+00, .9379E-02,-.1004E-03, .6593E-06,
     2 .3695E+02,-.8510E+00, .1037E-01,-.9624E-04, .5167E-06,
     2 .3768E+02,-.8574E+00, .1002E-01,-.9300E-04, .5332E-06,
     2 .3653E+02,-.8198E+00, .9766E-02,-.8351E-04, .3882E-06,
     2 .3441E+02,-.7423E+00, .8193E-02,-.7873E-04, .4911E-06,
     2 .5297E+02,-.1176E+01, .1408E-01,-.1506E-03, .9642E-06,
     2 .5513E+02,-.1255E+01, .1508E-01,-.1464E-03, .8476E-06,
     2 .5620E+02,-.1271E+01, .1481E-01,-.1387E-03, .8010E-06,
     2 .5408E+02,-.1219E+01, .1447E-01,-.1224E-03, .5719E-06,
     2 .5082E+02,-.1091E+01, .1205E-01,-.1176E-03, .7357E-06,
     2 .7874E+02,-.1749E+01, .2121E-01,-.2258E-03, .1406E-05,
     2 .8230E+02,-.1860E+01, .2225E-01,-.2246E-03, .1361E-05,
     2 .8369E+02,-.1886E+01, .2197E-01,-.2060E-03, .1185E-05,
     2 .8068E+02,-.1803E+01, .2095E-01,-.1809E-03, .9093E-06,
     2 .7519E+02,-.1597E+01, .1758E-01,-.1773E-03, .1140E-05,
     2 .1172E+03,-.2615E+01, .3219E-01,-.3423E-03, .2103E-05,
     2 .1231E+03,-.2805E+01, .3371E-01,-.3281E-03, .1912E-05,
     2 .1248E+03,-.2812E+01, .3300E-01,-.3006E-03, .1622E-05,
     2 .1202E+03,-.2669E+01, .3020E-01,-.2658E-03, .1463E-05,
     2 .1111E+03,-.2337E+01, .2557E-01,-.2674E-03, .1780E-05,
     3 .3100E+02,-.6945E+00, .8494E-02,-.8672E-04, .5088E-06,
     3 .3181E+02,-.7074E+00, .8455E-02,-.8803E-04, .5462E-06,
     3 .3326E+02,-.7418E+00, .8828E-02,-.8823E-04, .5238E-06,
     3 .3536E+02,-.7839E+00, .9161E-02,-.8986E-04, .5339E-06,
     3 .3711E+02,-.8247E+00, .9616E-02,-.9420E-04, .5646E-06,
     3 .4620E+02,-.1045E+01, .1277E-01,-.1280E-03, .7438E-06,
     3 .4726E+02,-.1059E+01, .1273E-01,-.1299E-03, .7868E-06,
     3 .4936E+02,-.1106E+01, .1326E-01,-.1290E-03, .7332E-06,
     3 .5266E+02,-.1160E+01, .1345E-01,-.1343E-03, .8198E-06,
     3 .5515E+02,-.1228E+01, .1445E-01,-.1392E-03, .8026E-06,
     3 .6904E+02,-.1562E+01, .1886E-01,-.1902E-03, .1144E-05,
     3 .7050E+02,-.1592E+01, .1925E-01,-.1900E-03, .1102E-05,
     3 .7362E+02,-.1645E+01, .1947E-01,-.1898E-03, .1105E-05,
     3 .7804E+02,-.1715E+01, .1994E-01,-.2001E-03, .1228E-05,
     3 .8155E+02,-.1825E+01, .2167E-01,-.2044E-03, .1135E-05,
     3 .1034E+03,-.2340E+01, .2796E-01,-.2808E-03, .1696E-05,
     3 .1056E+03,-.2374E+01, .2853E-01,-.2854E-03, .1686E-05,
     3 .1099E+03,-.2444E+01, .2875E-01,-.2822E-03, .1675E-05,
     3 .1153E+03,-.2530E+01, .2960E-01,-.2994E-03, .1835E-05,
     3 .1201E+03,-.2705E+01, .3262E-01,-.2964E-03, .1525E-05 /
      DATA (((( CL7CF(K,L,J,I), K = 1, 5), L = 1, 5), J = 1, 4),
     1                          I = 1, 3)                    /
     1 .2721E-01,-.5029E-04,-.7229E-06, .1343E-07,-.1062E-09,
     1 .2194E-01, .7539E-04,-.1816E-05, .7955E-09, .1054E-09,
     1 .1344E-01, .1345E-03,-.6031E-06,-.4300E-08,-.3154E-10,
     1 .6618E-02, .9637E-04,-.1803E-05,-.6878E-08, .3614E-09,
     1 .1711E-02, .3592E-04, .1396E-05,-.3379E-08,-.1991E-09,
     1 .4055E-01,-.5527E-04,-.1803E-05, .1873E-07,-.4836E-10,
     1 .3142E-01, .1157E-03,-.2265E-05, .3978E-08, .5934E-10,
     1 .1852E-01, .1935E-03,-.1201E-05,-.7214E-08, .4518E-10,
     1 .8721E-02, .1498E-03,-.1686E-05,-.1406E-07, .3448E-09,
     1 .2531E-02, .2676E-04, .1288E-05,-.7368E-09,-.1513E-09,
     1 .5893E-01,-.6185E-04,-.2884E-05, .2840E-07,-.6746E-10,
     1 .4433E-01, .1868E-03,-.3518E-05, .1134E-07, .5021E-10,
     1 .2432E-01, .2755E-03,-.9316E-06,-.1096E-07,-.1363E-10,
     1 .1031E-01, .1906E-03,-.8226E-06,-.1849E-07, .2199E-09,
     1 .3749E-02, .7716E-05, .9610E-06, .4006E-08,-.5038E-10,
     1 .8478E-01,-.2040E-04,-.4858E-05, .3770E-07,-.7883E-10,
     1 .6069E-01, .3674E-03,-.5458E-05, .3683E-08, .1553E-09,
     1 .3107E-01, .3821E-03,-.2662E-05,-.1749E-07, .3548E-09,
     1 .1220E-01, .2782E-03, .1685E-05,-.3191E-07,-.2148E-09,
     1 .5567E-02,-.3440E-04, .3315E-07, .1223E-07, .1653E-09,
     2 .7821E-01,-.5032E-03, .1815E-05,-.1798E-07, .2725E-09,
     2 .8295E-01,-.6342E-03, .3803E-05,-.6655E-09,-.1459E-09,
     2 .8570E-01,-.6197E-03, .1921E-05, .2405E-08, .8082E-10,
     2 .8002E-01,-.4733E-03, .2838E-05, .3014E-08,-.3074E-09,
     2 .7042E-01,-.3041E-03,-.4697E-06, .1385E-08, .1470E-09,
     2 .1157E+00,-.7264E-03, .3596E-05,-.3150E-07, .2669E-09,
     2 .1245E+00,-.9118E-03, .4689E-05,-.9272E-08,-.2104E-10,
     2 .1285E+00,-.8969E-03, .2715E-05,-.2958E-09, .1281E-09,
     2 .1180E+00,-.7211E-03, .3951E-05, .1018E-07,-.3941E-09,
     2 .1035E+00,-.4199E-03,-.6634E-06,-.3520E-08, .1953E-09,
     2 .1723E+00,-.1062E-02, .6825E-05,-.5327E-07, .1619E-09,
     2 .1871E+00,-.1329E-02, .6384E-05,-.2777E-07, .1529E-09,
     2 .1920E+00,-.1307E-02, .4023E-05,-.4661E-08, .1563E-09,
     2 .1769E+00,-.1041E-02, .4307E-05, .1321E-07,-.3727E-09,
     2 .1526E+00,-.5553E-03,-.1259E-05,-.1704E-07, .3511E-09,
     2 .2566E+00,-.1644E-02, .1261E-04,-.7588E-07, .5754E-10,
     2 .2816E+00,-.2071E-02, .1010E-04,-.2053E-07, .1099E-09,
     2 .2876E+00,-.1907E-02, .7436E-05,-.7658E-08,-.1866E-09,
     2 .2645E+00,-.1549E-02, .2845E-05, .2773E-07, .1846E-10,
     2 .2245E+00,-.7345E-03,-.2852E-05,-.4225E-07, .7253E-09,
     3 .5417E-01,-.2528E-03, .2038E-05,-.1069E-07,-.9237E-10,
     3 .5733E-01,-.2676E-03, .1178E-05,-.1113E-07, .6626E-10,
     3 .6326E-01,-.3189E-03, .1382E-05,-.8199E-08, .2497E-11,
     3 .7192E-01,-.3792E-03, .1203E-05,-.6301E-08, .3709E-10,
     3 .7876E-01,-.4561E-03, .1581E-05,-.7160E-08, .7944E-10,
     3 .8131E-01,-.4257E-03, .2970E-05,-.9726E-08,-.1144E-09,
     3 .8526E-01,-.4267E-03, .2139E-05,-.1287E-07, .3529E-10,
     3 .9381E-01,-.4804E-03, .2577E-05,-.8269E-08,-.1095E-09,
     3 .1076E+00,-.5418E-03, .1434E-05,-.1212E-07, .1226E-09,
     3 .1173E+00,-.6765E-03, .3044E-05,-.9734E-08,-.1606E-10,
     3 .1230E+00,-.6759E-03, .3263E-05,-.9553E-08, .4656E-10,
     3 .1282E+00,-.6832E-03, .3999E-05,-.9516E-08,-.1161E-09,
     3 .1408E+00,-.7172E-03, .3002E-05,-.8244E-08,-.4856E-10,
     3 .1592E+00,-.7890E-03, .2233E-05,-.1915E-07, .2194E-09,
     3 .1729E+00,-.1016E-02, .5181E-05,-.9115E-08,-.1553E-09,
     3 .1872E+00,-.1039E-02, .3510E-05,-.1482E-07, .1817E-09,
     3 .1948E+00,-.1027E-02, .5437E-05,-.2273E-07,-.9952E-10,
     3 .2118E+00,-.1057E-02, .3727E-05,-.1216E-07, .8117E-10,
     3 .2342E+00,-.1121E-02, .4051E-05,-.3507E-07, .2993E-09,
     3 .2527E+00,-.1504E-02, .9940E-05,-.1262E-08,-.7234E-09 /
C
      COMMON /BANDL7GH/ GWL7GH(7), CL7H2OGH(5,28,4),
     1                             CL7CO2GH(5,28,7), CL7O3GH(2)
      DATA GWL7GH / .18200, .13600, .05000, .01270, .00300,
     1              .00100, .00030 /
      DATA CL7O3GH / 40.73849, 43.96434 /
      DATA ((( CL7H2OGH(K,J,I), K = 1, 5), J = 1, 28), I = 1, 3) /
     4 .3979E-03, .1309E-04, .2319E-06, .2275E-08, .9044E-11,
     4 .5975E-03, .1801E-04, .3143E-06, .2752E-08, .6934E-11,
     4 .9489E-03, .2967E-04, .5084E-06, .3076E-08,-.5444E-11,
     4 .1616E-02, .4305E-04, .6158E-06, .3792E-08, .1404E-11,
     4 .2657E-02, .6118E-04, .7583E-06, .5784E-08, .2172E-10,
     4 .4296E-02, .8016E-04, .8953E-06, .8844E-08, .4875E-10,
     4 .7236E-02, .1351E-03, .1344E-05, .8520E-08, .3701E-10,
     4 .1304E-01, .1880E-03, .1277E-05, .1610E-07, .1506E-09,
     4 .2337E-01, .3183E-03, .1052E-05, .1674E-07, .3188E-09,
     4 .4154E-01, .5183E-03, .2804E-05, .2479E-07, .2173E-09,
     4 .1062E+00, .1367E-02, .7187E-05,-.9725E-08,-.2056E-09,
     4 .1507E+00, .1994E-02, .1516E-04,-.3444E-07,-.1388E-08,
     4 .2212E+00, .2635E-02, .1612E-04, .2766E-08,-.9770E-09,
     4 .3162E+00, .3718E-02, .1361E-04, .2238E-07, .6228E-09,
     4 .4491E+00, .5330E-02, .1967E-04, .1130E-07, .9380E-09,
     4 .6229E+00, .6924E-02, .2318E-04, .7405E-07, .1877E-08,
     4 .8901E+00, .1036E-01, .5791E-04,-.1261E-06,-.4170E-08,
     4 .1264E+01, .1329E-01, .5246E-04, .8264E-07, .8006E-09,
     4 .1837E+01, .1838E-01, .1531E-04, .1192E-06, .9628E-08,
     4 .2553E+01, .2565E-01, .1345E-03, .2501E-06,-.7128E-08,
     4 .3616E+01, .3494E-01, .1811E-03, .1876E-06,-.9945E-08,
     4 .5193E+01, .4251E-01, .2234E-03, .1113E-05,-.1956E-07,
     4 .7291E+01, .5548E-01, .2691E-03, .2257E-05,-.2823E-07,
     4 .9882E+01, .8229E-01,-.3997E-03, .2539E-05, .7505E-07,
     4 .1040E+02, .1364E+00, .2396E-03,-.2653E-05,-.8985E-08,
     4 .9481E+01, .1431E+00, .1698E-03,-.1884E-05, .4727E-07,
     4 .7102E+01, .1382E+00, .8000E-03,-.4169E-05,-.5986E-07,
     4 .4817E+01, .1006E+00, .4971E-03, .1741E-05, .3680E-07,
     5 .3073E-01, .1060E-02, .1685E-04, .2057E-06, .1376E-08,
     5 .3250E-01, .1081E-02, .1674E-04, .2063E-06, .1410E-08,
     5 .3729E-01, .1104E-02, .1512E-04, .2029E-06, .1591E-08,
     5 .4358E-01, .1219E-02, .1567E-04, .1915E-06, .1457E-08,
     5 .5140E-01, .1576E-02, .2414E-04, .1799E-06, .3830E-09,
     5 .6843E-01, .1849E-02, .2645E-04, .1623E-06,-.7235E-10,
     5 .9013E-01, .2220E-02, .2932E-04, .1516E-06,-.2603E-09,
     5 .1477E+00, .2575E-02, .1313E-04, .2283E-06, .3220E-08,
     5 .2275E+00, .3596E-02, .1003E-04, .2908E-06, .5134E-08,
     5 .3341E+00, .5577E-02, .2436E-04, .2137E-06, .4207E-08,
     5 .6167E+00, .9654E-02, .8056E-04, .2617E-06,-.1574E-08,
     5 .8408E+00, .1139E-01, .7003E-04, .4346E-06, .2950E-08,
     5 .1122E+01, .1629E-01, .1129E-03, .1767E-06,-.1711E-08,
     5 .1604E+01, .2096E-01, .1567E-03, .4008E-06,-.8837E-08,
     5 .2256E+01, .3159E-01, .2720E-03, .2558E-06,-.2096E-07,
     5 .3228E+01, .3908E-01, .2333E-03, .1666E-05, .1700E-08,
     5 .4433E+01, .5401E-01, .8836E-04, .1209E-05, .3277E-07,
     5 .6214E+01, .7139E-01, .2844E-03, .8659E-06,-.8975E-08,
     5 .8883E+01, .1214E+00, .4692E-03,-.3952E-05,-.4179E-07,
     5 .1070E+02, .1283E+00, .1178E-02, .1746E-05,-.1318E-06,
     5 .1423E+02, .1890E+00, .2108E-04,-.1841E-05, .4531E-07,
     5 .1308E+02, .2516E+00, .1488E-02,-.8041E-05,-.1593E-06,
     5 .1082E+02, .2402E+00, .4970E-03,-.3117E-05, .1181E-06,
     5 .5873E+01, .2170E+00, .2438E-02,-.6079E-05,-.1877E-06,
     5 .2984E+01, .1027E+00, .1260E-02, .7025E-05, .2700E-07,
     5 .1567E+01, .6182E-01, .1249E-02, .4508E-05,-.6663E-07,
     5 .1627E+01, .2983E-01, .1840E-03, .5014E-05, .6238E-07,
     5 .1814E+01, .4301E-01, .3540E-03,-.6079E-06,-.2215E-07,
     6 .2679E+01, .8688E-01, .1266E-02, .5835E-05,-.2292E-07,
     6 .2674E+01, .9032E-01, .1363E-02, .5089E-05,-.4266E-07,
     6 .2907E+01, .9037E-01, .1116E-02, .5038E-05,-.1392E-08,
     6 .3053E+01, .8953E-01, .9663E-03, .5114E-05, .2324E-07,
     6 .3393E+01, .9585E-01, .9985E-03, .3540E-05, .1842E-08,
     6 .3641E+01, .8739E-01, .8533E-03, .6089E-05, .3778E-07,
     6 .3979E+01, .8723E-01, .8792E-03, .6312E-05, .2811E-07,
     6 .4731E+01, .1070E+00, .1184E-02, .4556E-05,-.2335E-07,
     6 .5910E+01, .1202E+00, .1488E-02, .3533E-05,-.9056E-07,
     6 .8394E+01, .1517E+00, .1013E-02, .3045E-05, .9213E-08,
     6 .1445E+02, .2547E+00, .9205E-03,-.6289E-05, .2503E-07,
     6 .1641E+02, .2691E+00, .2174E-02,-.5779E-05,-.2074E-06,
     6 .1824E+02, .2805E+00, .2827E-02, .1470E-05,-.1362E-06,
     6 .2484E+02, .2895E+00, .2182E-02, .2012E-04, .9200E-07,
     6 .3228E+02, .2715E+00,-.7519E-03, .2783E-04, .4698E-06,
     6 .4015E+02, .4094E+00, .2015E-03, .1272E-04, .1433E-06,
     6 .4570E+02, .6543E+00, .2055E-02,-.6647E-05,-.8492E-07,
     6 .4541E+02, .7479E+00, .4208E-03,-.8566E-05, .3976E-06,
     6 .3539E+02, .6533E+00, .2515E-02, .8626E-05, .2095E-06,
     6 .2673E+02, .6560E+00, .2446E-02, .1054E-05, .3383E-06,
     6 .1179E+02, .4299E+00, .5160E-02, .1176E-04,-.1008E-06,
     6 .6849E+01, .1998E+00, .1027E-02, .2448E-04, .4590E-06,
     6 .2467E+01, .1354E+00, .3275E-02, .7239E-05,-.2512E-06,
     6 .2487E+01, .4085E-01, .3215E-03, .1424E-04, .1790E-06,
     6 .2652E+01, .5485E-01, .4623E-03, .9730E-06,-.4515E-10,
     6 .2757E+01, .5373E-01, .5512E-03, .1454E-05,-.1937E-07,
     6 .2587E+01, .5036E-01, .5962E-03, .2420E-05,-.1868E-07,
     6 .1840E+01, .2508E-01, .5036E-03, .6535E-05, .1483E-07 /
      DATA ((( CL7H2OGH(K,J,I), K = 1, 5), J = 1, 28), I = 4, 4) /
     7 .2162E+03, .3719E+01, .1760E-01, .6433E-04, .9080E-06,
     7 .2162E+03, .3706E+01, .1721E-01, .6730E-04, .9864E-06,
     7 .2152E+03, .3705E+01, .1822E-01, .6752E-04, .8199E-06,
     7 .2146E+03, .3707E+01, .1880E-01, .6733E-04, .7229E-06,
     7 .2131E+03, .3677E+01, .1858E-01, .7369E-04, .8188E-06,
     7 .2121E+03, .3715E+01, .1902E-01, .5693E-04, .7750E-06,
     7 .2106E+03, .4204E+01, .3293E-01,-.4423E-04,-.2055E-05,
     7 .2228E+03, .3967E+01, .2013E-01,-.7162E-05,-.1466E-06,
     7 .2174E+03, .3899E+01, .1884E-01,-.3389E-05, .1150E-06,
     7 .2271E+03, .3934E+01, .2140E-01, .1457E-04,-.1974E-06,
     7 .2004E+03, .3456E+01, .2094E-01, .5202E-04,-.1874E-06,
     7 .1905E+03, .3366E+01, .1588E-01, .4846E-04, .7028E-06,
     7 .1803E+03, .3253E+01, .1270E-01, .2200E-04, .4645E-06,
     7 .1488E+03, .3143E+01, .1470E-01,-.5471E-04,-.3892E-06,
     7 .1114E+03, .3069E+01, .2493E-01,-.8397E-04,-.1775E-05,
     7 .6775E+02, .2405E+01, .2132E-01,-.4006E-04,-.7314E-06,
     7 .2845E+02, .1198E+01, .1481E-01, .4502E-04,-.6451E-07,
     7 .3587E+01, .5628E+00, .1903E-01, .5347E-04,-.1599E-05,
     7 .3853E+01, .2854E+00, .8946E-02, .3495E-04,-.6220E-06,
     7 .4115E+01, .3230E-01,-.1339E-03, .9949E-06, .1713E-07,
     7 .4514E+01, .4306E-01, .1045E-03,-.3784E-06,-.3653E-07,
     7 .5037E+01, .5328E-01,-.1313E-03,-.1217E-05, .4063E-08,
     7 .5058E+01, .6589E-01,-.4212E-04,-.1356E-05, .1604E-07,
     7 .4487E+01, .6287E-01, .1931E-03, .7326E-06,-.1275E-08,
     7 .4122E+01, .6398E-01,-.1124E-04,-.3712E-06, .1871E-07,
     7 .3504E+01, .7436E-01,-.1268E-03,-.3277E-05, .5389E-07,
     7 .2483E+01, .7784E-01, .3708E-03,-.4513E-05,-.1066E-07,
     7 .3159E+01, .7310E-01,-.2678E-03,-.6098E-05, .5793E-07 /
      DATA ((( CL7CO2GH(K,J,I), K = 1, 5), J = 1, 28), I = 1, 3) /
     4 .7533E-03, .1759E-04, .2086E-06, .1701E-08, .7958E-11,
     4 .1240E-02, .2724E-04, .3004E-06, .1972E-08, .6215E-11,
     4 .2071E-02, .4233E-04, .4310E-06, .2389E-08, .4750E-11,
     4 .3501E-02, .6658E-04, .6299E-06, .3096E-08, .3672E-11,
     4 .6005E-02, .1055E-03, .9128E-06, .4268E-08, .6298E-11,
     4 .1045E-01, .1707E-03, .1361E-05, .5896E-08, .7869E-11,
     4 .1848E-01, .2775E-03, .2017E-05, .8748E-08, .1520E-10,
     4 .3312E-01, .4599E-03, .3097E-05, .1267E-07, .9837E-11,
     4 .6036E-01, .7744E-03, .4852E-05, .1960E-07,-.1083E-10,
     4 .1122E+00, .1343E-02, .7612E-05, .2884E-07,-.5713E-10,
     4 .2939E+00, .3218E-02, .1513E-04, .5520E-07,-.1293E-09,
     4 .4231E+00, .4509E-02, .1931E-04, .7009E-07,-.8888E-10,
     4 .6096E+00, .6369E-02, .2573E-04, .8445E-07,-.1876E-09,
     4 .8783E+00, .8991E-02, .3577E-04, .1071E-06,-.5840E-09,
     4 .1263E+01, .1267E-01, .4923E-04, .1465E-06,-.1119E-08,
     4 .1814E+01, .1777E-01, .6644E-04, .2139E-06,-.1677E-08,
     4 .2595E+01, .2469E-01, .8332E-04, .3488E-06,-.1476E-08,
     4 .3695E+01, .3432E-01, .1039E-03, .4822E-06,-.1712E-08,
     4 .5231E+01, .4768E-01, .1220E-03, .5797E-06,-.1247E-08,
     4 .7349E+01, .6585E-01, .1372E-03, .6522E-06, .3123E-10,
     4 .1022E+02, .8967E-01, .1715E-03, .7535E-06,-.1828E-08,
     4 .1416E+02, .1201E+00, .1571E-03, .8937E-06, .1947E-08,
     4 .1954E+02, .1614E+00, .9713E-04, .7591E-06, .7562E-08,
     4 .2670E+02, .2189E+00, .1658E-04, .3104E-06, .1403E-07,
     4 .3583E+02, .2967E+00,-.1940E-03,-.3044E-06, .3986E-07,
     4 .4628E+02, .4074E+00,-.2678E-03,-.3593E-05, .4834E-07,
     4 .5626E+02, .5270E+00, .2664E-03,-.7034E-05, .6535E-08,
     4 .6676E+02, .6160E+00, .3844E-03,-.6498E-05, .1904E-07,
     5 .5221E-01, .1805E-02, .2934E-04, .2960E-06, .1506E-08,
     5 .6099E-01, .1958E-02, .3034E-04, .2919E-06, .1428E-08,
     5 .7430E-01, .2183E-02, .3210E-04, .2900E-06, .1305E-08,
     5 .9564E-01, .2507E-02, .3425E-04, .2939E-06, .1267E-08,
     5 .1308E+00, .2982E-02, .3636E-04, .3019E-06, .1392E-08,
     5 .1887E+00, .3740E-02, .4048E-04, .3122E-06, .1437E-08,
     5 .2868E+00, .4941E-02, .4705E-04, .3218E-06, .1375E-08,
     5 .4563E+00, .6876E-02, .5772E-04, .3267E-06, .1108E-08,
     5 .7576E+00, .9984E-02, .7128E-04, .3349E-06, .1147E-08,
     5 .1301E+01, .1510E-01, .9173E-04, .3467E-06, .1288E-08,
     5 .3084E+01, .3015E-01, .1499E-03, .3736E-06, .1305E-08,
     5 .4287E+01, .3959E-01, .1842E-03, .3846E-06, .1300E-08,
     5 .5989E+01, .5211E-01, .2223E-03, .4357E-06, .2183E-08,
     5 .8378E+01, .6927E-01, .2770E-03, .4164E-06, .2439E-08,
     5 .1174E+02, .9184E-01, .3309E-03, .4309E-06, .5361E-08,
     5 .1644E+02, .1224E+00, .3954E-03, .2477E-06, .8943E-08,
     5 .2294E+02, .1624E+00, .5228E-03,-.4672E-07, .4298E-08,
     5 .3187E+02, .2121E+00, .6170E-03,-.3082E-06, .1099E-07,
     5 .4400E+02, .2719E+00, .6996E-03,-.8778E-06, .2000E-07,
     5 .6025E+02, .3486E+00, .8769E-03,-.2275E-05, .1487E-07,
     5 .8148E+02, .4338E+00, .1275E-02,-.2286E-05,-.2183E-07,
     5 .1091E+03, .5291E+00, .1301E-02,-.1092E-05,-.1145E-07,
     5 .1433E+03, .6458E+00, .1046E-02,-.2897E-05, .3588E-07,
     5 .1827E+03, .7671E+00, .1326E-02,-.8701E-05, .4600E-07,
     5 .2277E+03, .8028E+00, .2318E-02,-.7950E-05,-.1857E-07,
     5 .2850E+03, .7545E+00, .2500E-02,-.8734E-05, .3775E-07,
     5 .3618E+03, .6011E+00, .2124E-02,-.7617E-05, .1066E-06,
     5 .4619E+03, .3803E+00, .2990E-02,-.7075E-05,-.9234E-07,
     6 .4780E+01, .1369E+00, .1601E-02, .8635E-05, .1564E-07,
     6 .4916E+01, .1363E+00, .1576E-02, .8748E-05, .1836E-07,
     6 .5166E+01, .1368E+00, .1545E-02, .8726E-05, .2062E-07,
     6 .5579E+01, .1392E+00, .1525E-02, .8545E-05, .1990E-07,
     6 .6258E+01, .1443E+00, .1500E-02, .8234E-05, .1968E-07,
     6 .7359E+01, .1546E+00, .1495E-02, .7704E-05, .1741E-07,
     6 .9122E+01, .1720E+00, .1526E-02, .7079E-05, .1276E-07,
     6 .1195E+02, .1972E+00, .1548E-02, .6922E-05, .1710E-07,
     6 .1663E+02, .2384E+00, .1626E-02, .6359E-05, .1751E-07,
     6 .2454E+02, .3044E+00, .1758E-02, .5424E-05, .1714E-07,
     6 .4848E+02, .4761E+00, .2260E-02, .4878E-05, .3592E-09,
     6 .6383E+02, .5736E+00, .2499E-02, .4589E-05,-.7953E-08,
     6 .8489E+02, .6976E+00, .2668E-02, .3466E-05,-.5880E-08,
     6 .1134E+03, .8629E+00, .2916E-02, .3468E-06, .2279E-08,
     6 .1520E+03, .1073E+01, .3294E-02,-.2331E-05, .2682E-07,
     6 .2033E+03, .1319E+01, .4290E-02,-.2585E-05,-.5028E-07,
     6 .2706E+03, .1592E+01, .5192E-02, .7664E-06,-.1258E-06,
     6 .3575E+03, .1906E+01, .4917E-02, .9025E-05,-.6084E-07,
     6 .4645E+03, .2358E+01, .4227E-02, .5232E-05,-.2967E-07,
     6 .5874E+03, .2909E+01, .2957E-02,-.1132E-04, .1926E-06,
     6 .7166E+03, .3404E+01, .3177E-02,-.3915E-04, .4958E-06,
     6 .8538E+03, .3568E+01, .8872E-02,-.6325E-04, .7559E-07,
     6 .1026E+04, .3023E+01, .1190E-01,-.2772E-04,-.1500E-06,
     6 .1264E+04, .2224E+01, .8202E-02,-.7830E-05, .6620E-07,
     6 .1558E+04, .1497E+01, .5592E-02,-.1036E-04,-.8107E-07,
     6 .1876E+04, .6267E+00, .2353E-02, .7971E-05,-.2378E-07,
     6 .2159E+04,-.1868E+00,-.2069E-03, .1924E-04,-.8951E-07,
     6 .2348E+04,-.8938E+00,-.4795E-02, .1456E-04, .5410E-06 /
      DATA ((( CL7CO2GH(K,J,I), K = 1, 5), J = 1, 28), I = 4, 6) /
     7 .2615E+03, .5304E+01, .3627E-01, .5074E-04,-.3389E-06,
     7 .2617E+03, .5303E+01, .3639E-01, .5039E-04,-.3520E-06,
     7 .2623E+03, .5296E+01, .3652E-01, .5041E-04,-.3596E-06,
     7 .2639E+03, .5276E+01, .3662E-01, .5097E-04,-.3578E-06,
     7 .2679E+03, .5237E+01, .3681E-01, .5271E-04,-.3778E-06,
     7 .2771E+03, .5164E+01, .3685E-01, .5950E-04,-.4241E-06,
     7 .2978E+03, .5056E+01, .3443E-01, .7421E-04,-.1410E-06,
     7 .3368E+03, .5065E+01, .3200E-01, .7296E-04, .2043E-07,
     7 .4073E+03, .5236E+01, .2858E-01, .5791E-04, .2245E-06,
     7 .5211E+03, .5557E+01, .3079E-01, .5554E-04,-.4996E-06,
     7 .8666E+03, .6459E+01, .2403E-01, .7457E-04, .3639E-06,
     7 .1085E+04, .7166E+01, .1887E-01, .4120E-04, .8569E-06,
     7 .1368E+04, .8129E+01, .1834E-01, .8641E-05, .6606E-06,
     7 .1734E+04, .9322E+01, .2297E-01, .3645E-05,-.7712E-06,
     7 .2203E+04, .1051E+02, .1850E-01, .7040E-04,-.9868E-06,
     7 .2771E+04, .1259E+02, .4016E-02,-.1628E-04, .7452E-06,
     7 .3414E+04, .1509E+02,-.3010E-02,-.1970E-03, .2855E-05,
     7 .4135E+04, .1642E+02, .2275E-01,-.3437E-03, .4329E-06,
     7 .5048E+04, .1487E+02, .3174E-01,-.2763E-03, .5903E-06,
     7 .6264E+04, .1120E+02, .3591E-01,-.1358E-03,-.9681E-06,
     7 .7830E+04, .6923E+01, .1644E-01, .2735E-04,-.1368E-05,
     7 .9576E+04, .3529E+01,-.2416E-01, .1569E-03, .4027E-06,
     7 .1108E+05, .2838E+01,-.4676E-01, .2050E-04, .5772E-06,
     7 .1185E+05, .2034E+01,-.4604E-01, .4265E-04, .1925E-06,
     7 .1169E+05, .6925E+00,-.5102E-01, .2157E-04, .1562E-05,
     7 .1066E+05,-.1466E+00,-.3748E-01,-.8924E-05, .9231E-06,
     7 .9231E+04,-.8276E+00,-.3130E-01, .5709E-04, .1299E-05,
     7 .7788E+04, .9468E+00,-.1446E-01,-.1953E-04, .2703E-06,
     8 .4476E+04, .5788E+02, .1918E+00,-.2188E-03,-.1428E-05,
     8 .4477E+04, .5788E+02, .1916E+00,-.2194E-03,-.1407E-05,
     8 .4479E+04, .5788E+02, .1915E+00,-.2201E-03,-.1383E-05,
     8 .4482E+04, .5787E+02, .1913E+00,-.2214E-03,-.1367E-05,
     8 .4488E+04, .5785E+02, .1912E+00,-.2248E-03,-.1350E-05,
     8 .4500E+04, .5781E+02, .1915E+00,-.2427E-03,-.1268E-05,
     8 .4525E+04, .5766E+02, .1989E+00,-.2801E-03,-.2176E-05,
     8 .4590E+04, .5682E+02, .2049E+00,-.2793E-03,-.2457E-05,
     8 .4758E+04, .5483E+02, .2150E+00,-.1972E-03,-.3812E-05,
     8 .5198E+04, .5185E+02, .1836E+00,-.6698E-04,-.7395E-06,
     8 .6666E+04, .4913E+02, .1324E+00,-.1569E-03,-.7295E-06,
     8 .7490E+04, .4931E+02, .1489E+00,-.2815E-03,-.5456E-05,
     8 .8585E+04, .4861E+02, .9724E-01,-.2913E-03, .3049E-05,
     8 .9957E+04, .4754E+02, .7270E-01,-.7409E-03, .8337E-05,
     8 .1172E+05, .4459E+02, .1294E+00,-.1470E-02, .5003E-05,
     8 .1419E+05, .3088E+02, .2149E+00,-.9643E-03,-.3364E-05,
     8 .1783E+05, .8358E+01, .1741E+00, .9848E-04,-.2619E-05,
     8 .2259E+05,-.1332E+02, .6490E-01, .5333E-03, .1723E-05,
     8 .2785E+05,-.2933E+02,-.2199E-01, .5109E-03, .3031E-05,
     8 .3240E+05,-.3835E+02,-.5586E-01, .4911E-03,-.6266E-06,
     8 .3522E+05,-.4264E+02,-.1567E+00, .6481E-03, .6558E-05,
     8 .3537E+05,-.3794E+02,-.1826E+00, .4845E-03, .5126E-05,
     8 .3304E+05,-.2947E+02,-.1975E+00, .6451E-03, .4877E-05,
     8 .2958E+05,-.1716E+02,-.1930E+00, .3573E-03, .2223E-05,
     8 .2679E+05,-.5498E+01,-.1766E+00, .2773E-03, .1864E-05,
     8 .2521E+05, .1577E+01,-.1613E+00, .1568E-03, .3465E-05,
     8 .2448E+05, .4568E+01,-.1126E+00,-.1242E-03, .2870E-06,
     8 .2433E+05,-.1987E+01,-.1213E+00, .1448E-03, .3506E-05,
     9 .4412E+05, .3472E+03,-.2846E+00,-.6166E-02, .8124E-05,
     9 .4413E+05, .3471E+03,-.2844E+00,-.6166E-02, .8085E-05,
     9 .4414E+05, .3471E+03,-.2841E+00,-.6166E-02, .8041E-05,
     9 .4415E+05, .3469E+03,-.2836E+00,-.6166E-02, .7975E-05,
     9 .4419E+05, .3466E+03,-.2833E+00,-.6166E-02, .7981E-05,
     9 .4426E+05, .3460E+03,-.2828E+00,-.6177E-02, .8074E-05,
     9 .4441E+05, .3448E+03,-.2810E+00,-.6204E-02, .8058E-05,
     9 .4472E+05, .3423E+03,-.2720E+00,-.6254E-02, .6958E-05,
     9 .4535E+05, .3371E+03,-.2607E+00,-.6429E-02, .7590E-05,
     9 .4661E+05, .3246E+03,-.2078E+00,-.6677E-02, .7502E-05,
     9 .5167E+05, .2660E+03,-.5957E-01,-.6248E-02, .8927E-05,
     9 .5572E+05, .2209E+03,-.9796E-01,-.5173E-02, .2121E-04,
     9 .6149E+05, .1612E+03,-.9377E-01,-.3597E-02, .1533E-04,
     9 .6919E+05, .9042E+02,-.1655E+00,-.2015E-02, .1561E-04,
     9 .7826E+05, .1096E+02,-.2393E+00, .1775E-02,-.5753E-05,
     9 .8635E+05,-.3982E+02,-.4361E+00, .1924E-02, .1590E-04,
     9 .9097E+05,-.6474E+02,-.3228E+00, .6485E-04,-.1016E-04,
     9 .9124E+05,-.7755E+02,-.5427E+00,-.6166E-03, .3452E-04,
     9 .8580E+05,-.7678E+02,-.3057E+00,-.1002E-02, .1055E-04,
     9 .7662E+05,-.8212E+02,-.3641E+00,-.2757E-03, .3060E-04,
     9 .6534E+05,-.7810E+02,-.1041E+00,-.4618E-03, .8734E-05,
     9 .5522E+05,-.8662E+02, .5203E-01,-.1841E-03, .1274E-05,
     9 .4858E+05,-.1007E+03, .9983E-01, .5084E-04,-.8371E-06,
     9 .4488E+05,-.1080E+03, .6578E-02, .3111E-03, .1139E-04,
     9 .4150E+05,-.1065E+03, .1077E+00,-.9938E-04, .4886E-06,
     9 .3871E+05,-.1083E+03, .1047E+00, .5064E-03,-.4781E-06,
     9 .3633E+05,-.1004E+03, .1903E-01, .5945E-03, .5293E-05,
     9 .3363E+05,-.8354E+02,-.1288E-01, .4747E-03, .1656E-05 /
      DATA ((( CL7CO2GH(K,J,I), K = 1, 5), J = 1, 28), I = 7, 7) /
     a .5140E+06,-.1411E+04,-.1048E+01, .1654E-01, .1947E-04,
     a .5140E+06,-.1411E+04,-.1047E+01, .1652E-01, .1922E-04,
     a .5140E+06,-.1411E+04,-.1047E+01, .1654E-01, .1930E-04,
     a .5140E+06,-.1410E+04,-.1050E+01, .1654E-01, .1964E-04,
     a .5136E+06,-.1409E+04,-.1054E+01, .1654E-01, .2006E-04,
     a .5132E+06,-.1407E+04,-.1056E+01, .1656E-01, .2031E-04,
     a .5120E+06,-.1402E+04,-.1066E+01, .1660E-01, .2123E-04,
     a .5100E+06,-.1392E+04,-.1083E+01, .1665E-01, .2266E-04,
     a .5056E+06,-.1372E+04,-.1124E+01, .1678E-01, .2670E-04,
     a .4968E+06,-.1331E+04,-.1149E+01, .1651E-01, .2540E-04,
     a .4664E+06,-.1185E+04,-.1219E+01, .1504E-01, .1425E-04,
     a .4456E+06,-.1085E+04,-.1202E+01, .1347E-01, .5864E-05,
     a .4172E+06,-.9484E+03,-.1014E+01, .1018E-01,-.1949E-04,
     a .3803E+06,-.7892E+03,-.8996E+00, .8948E-02,-.2216E-04,
     a .3355E+06,-.6040E+03,-.1064E+01, .2726E-02, .5016E-04,
     a .2868E+06,-.4540E+03,-.9408E+00, .9592E-03, .2099E-04,
     a .2383E+06,-.3363E+03,-.9800E+00, .2306E-02, .5124E-04,
     a .1916E+06,-.2437E+03,-.3055E+00, .3395E-02,-.4324E-04,
     a .1512E+06,-.1799E+03,-.5340E+00, .3604E-02, .6332E-06,
     a .1164E+06,-.1153E+03,-.2437E+00, .9536E-03,-.1570E-04,
     a .8864E+05,-.8796E+02,-.1052E+00,-.3615E-03,-.8576E-05,
     a .6732E+05,-.7572E+02,-.5404E-01,-.9768E-03, .1252E-04,
     a .5116E+05,-.7624E+02, .1354E+00,-.2120E-02, .2058E-04,
     a .4036E+05,-.1004E+03, .5168E+00,-.1731E-02,-.5504E-05,
     a .3519E+05,-.1283E+03, .3094E+00, .2044E-03, .1278E-04,
     a .3240E+05,-.1182E+03, .1915E+00,-.1668E-03, .1033E-05,
     a .2966E+05,-.9956E+02, .9680E-01, .1245E-03, .1818E-05,
     a .2662E+05,-.7916E+02, .2546E-01, .3327E-03, .1057E-05 /
C
C----------------------------------------------------------------------C
C     SPECTRAL REGION: 340-540 CM^-1. GAS: H2O                         C
C----------------------------------------------------------------------C
C
      COMMON /BANDL8/ GWL8(6), CL8H2O(5,18,6), CL8CS(5,6,4),
     1                         CL8CF(5,6,4)
      DATA GWL8 / .09800, .27720, .22500, .20600, .12500, .05300 /
      DATA ((( CL8H2O(K,J,I), K = 1, 5), J = 1, 18), I = 1, 4) /
     1 .3229E-03, .4497E-05, .2238E-07, .1178E-09, .6098E-12,
     1 .4900E-03, .6826E-05, .3397E-07, .1787E-09, .9256E-12,
     1 .7434E-03, .1036E-04, .5152E-07, .2711E-09, .1406E-11,
     1 .1127E-02, .1571E-04, .7807E-07, .4108E-09, .2145E-11,
     1 .1708E-02, .2381E-04, .1181E-06, .6211E-09, .3287E-11,
     1 .2588E-02, .3608E-04, .1787E-06, .9398E-09, .5027E-11,
     1 .3918E-02, .5465E-04, .2702E-06, .1421E-08, .7680E-11,
     1 .5931E-02, .8277E-04, .4081E-06, .2146E-08, .1183E-10,
     1 .8965E-02, .1252E-03, .6172E-06, .3235E-08, .1795E-10,
     1 .1354E-01, .1893E-03, .9381E-06, .4878E-08, .2651E-10,
     1 .2044E-01, .2859E-03, .1427E-05, .7366E-08, .3840E-10,
     1 .3081E-01, .4309E-03, .2152E-05, .1118E-07, .5852E-10,
     1 .4640E-01, .6479E-03, .3221E-05, .1699E-07, .9206E-10,
     1 .6971E-01, .9730E-03, .4778E-05, .2551E-07, .1489E-09,
     1 .1044E+00, .1458E-02, .7169E-05, .3796E-07, .2206E-09,
     1 .1558E+00, .2172E-02, .1060E-04, .5747E-07, .3551E-09,
     1 .2314E+00, .3216E-02, .1582E-04, .8670E-07, .5066E-09,
     1 .3426E+00, .4692E-02, .2161E-04, .1347E-06, .1031E-08,
     2 .2104E-02, .3167E-04, .1722E-06, .2633E-09,-.1578E-11,
     2 .3167E-02, .4767E-04, .2592E-06, .3958E-09,-.2383E-11,
     2 .4764E-02, .7173E-04, .3905E-06, .5944E-09,-.3642E-11,
     2 .7163E-02, .1079E-03, .5888E-06, .8916E-09,-.5679E-11,
     2 .1077E-01, .1623E-03, .8888E-06, .1334E-08,-.9037E-11,
     2 .1618E-01, .2442E-03, .1342E-05, .2002E-08,-.1440E-10,
     2 .2430E-01, .3669E-03, .2023E-05, .3027E-08,-.2209E-10,
     2 .3649E-01, .5507E-03, .3042E-05, .4617E-08,-.3314E-10,
     2 .5474E-01, .8258E-03, .4574E-05, .7060E-08,-.4968E-10,
     2 .8204E-01, .1238E-02, .6861E-05, .1084E-07,-.7308E-10,
     2 .1228E+00, .1851E-02, .1022E-04, .1654E-07,-.9934E-10,
     2 .1835E+00, .2766E-02, .1527E-04, .2466E-07,-.1446E-09,
     2 .2736E+00, .4125E-02, .2279E-04, .3703E-07,-.2095E-09,
     2 .4068E+00, .6129E-02, .3408E-04, .5640E-07,-.3322E-09,
     2 .6032E+00, .9050E-02, .5009E-04, .8847E-07,-.4121E-09,
     2 .8918E+00, .1329E-01, .7219E-04, .1369E-06,-.3589E-09,
     2 .1312E+01, .1948E-01, .1064E-03, .1939E-06,-.7568E-09,
     2 .1923E+01, .2839E-01, .1537E-03, .2804E-06,-.1043E-08,
     3 .7446E-02, .1111E-03, .6046E-06, .8289E-09,-.4900E-11,
     3 .1119E-01, .1671E-03, .9134E-06, .1258E-08,-.7887E-11,
     3 .1682E-01, .2510E-03, .1373E-05, .1945E-08,-.1124E-10,
     3 .2526E-01, .3770E-03, .2061E-05, .2983E-08,-.1581E-10,
     3 .3794E-01, .5659E-03, .3090E-05, .4602E-08,-.2153E-10,
     3 .5695E-01, .8488E-03, .4629E-05, .7068E-08,-.2907E-10,
     3 .8540E-01, .1272E-02, .6943E-05, .1085E-07,-.4262E-10,
     3 .1279E+00, .1903E-02, .1036E-04, .1674E-07,-.5470E-10,
     3 .1915E+00, .2846E-02, .1543E-04, .2488E-07,-.7320E-10,
     3 .2860E+00, .4248E-02, .2294E-04, .3714E-07,-.9169E-10,
     3 .4261E+00, .6327E-02, .3436E-04, .5581E-07,-.1570E-09,
     3 .6330E+00, .9376E-02, .5132E-04, .8718E-07,-.2525E-09,
     3 .9376E+00, .1384E-01, .7570E-04, .1350E-06,-.3252E-09,
     3 .1385E+01, .2027E-01, .1078E-03, .2157E-06, .1516E-09,
     3 .2038E+01, .2969E-01, .1553E-03, .2960E-06, .2982E-09,
     3 .2984E+01, .4326E-01, .2229E-03, .3967E-06, .3153E-09,
     3 .4333E+01, .6276E-01, .3264E-03, .5118E-06,-.9117E-09,
     3 .6276E+01, .9074E-01, .4698E-03, .6813E-06,-.1851E-08,
     4 .3967E-01, .6059E-03, .3125E-05, .2384E-08,-.3299E-10,
     4 .5947E-01, .9066E-03, .4649E-05, .3729E-08,-.4392E-10,
     4 .8910E-01, .1356E-02, .6930E-05, .5633E-08,-.6395E-10,
     4 .1335E+00, .2026E-02, .1038E-04, .8408E-08,-.1025E-09,
     4 .1998E+00, .3029E-02, .1547E-04, .1266E-07,-.1495E-09,
     4 .2991E+00, .4525E-02, .2298E-04, .1888E-07,-.2122E-09,
     4 .4470E+00, .6756E-02, .3422E-04, .2749E-07,-.3131E-09,
     4 .6667E+00, .1007E-01, .5103E-04, .4005E-07,-.4745E-09,
     4 .9922E+00, .1498E-01, .7611E-04, .6001E-07,-.7425E-09,
     4 .1472E+01, .2220E-01, .1135E-03, .8729E-07,-.1228E-08,
     4 .2175E+01, .3269E-01, .1657E-03, .1423E-06,-.1472E-08,
     4 .3201E+01, .4791E-01, .2404E-03, .2029E-06,-.1948E-08,
     4 .4684E+01, .6978E-01, .3502E-03, .2804E-06,-.3323E-08,
     4 .6810E+01, .1010E+00, .5049E-03, .3511E-06,-.5520E-08,
     4 .9836E+01, .1453E+00, .7112E-03, .4555E-06,-.6331E-08,
     4 .1413E+02, .2078E+00, .9918E-03, .6837E-06,-.4149E-08,
     4 .2014E+02, .2960E+00, .1430E-02, .7609E-06,-.1040E-07,
     4 .2842E+02, .4146E+00, .2034E-02, .1088E-05,-.2282E-07 /
      DATA ((( CL8H2O(K,J,I), K = 1, 5), J = 1, 18), I = 5, 6) /
     5 .3133E+00, .4523E-02, .1887E-04, .6784E-09,-.1336E-09,
     5 .4679E+00, .6729E-02, .2771E-04, .4290E-09,-.1794E-09,
     5 .6987E+00, .1002E-01, .4070E-04,-.3248E-08,-.2409E-09,
     5 .1043E+01, .1493E-01, .5981E-04,-.8776E-08,-.3005E-09,
     5 .1555E+01, .2224E-01, .8850E-04,-.2206E-07,-.4725E-09,
     5 .2316E+01, .3304E-01, .1310E-03,-.2990E-07,-.6487E-09,
     5 .3444E+01, .4901E-01, .1929E-03,-.3852E-07,-.7544E-09,
     5 .5109E+01, .7250E-01, .2833E-03,-.6448E-07,-.1062E-08,
     5 .7551E+01, .1068E+00, .4147E-03,-.9567E-07,-.1432E-08,
     5 .1111E+02, .1566E+00, .6031E-03,-.1274E-06,-.2022E-08,
     5 .1629E+02, .2283E+00, .8610E-03,-.2567E-06,-.2455E-08,
     5 .2374E+02, .3310E+00, .1219E-02,-.5411E-06,-.3350E-08,
     5 .3439E+02, .4767E+00, .1697E-02,-.8358E-06, .3353E-09,
     5 .4947E+02, .6840E+00, .2392E-02,-.1672E-05,-.1353E-08,
     5 .7058E+02, .9732E+00, .3392E-02,-.3087E-05,-.1572E-07,
     5 .9985E+02, .1362E+01, .4450E-02,-.5045E-05,-.8712E-08,
     5 .1385E+03, .1876E+01, .5845E-02,-.7922E-05,-.8375E-08,
     5 .1862E+03, .2526E+01, .7113E-02,-.1443E-04, .1020E-06,
     6 .4125E+01, .6340E-01, .2952E-03, .1065E-06,-.1185E-08,
     6 .6031E+01, .9171E-01, .4119E-03, .7714E-07,-.8214E-09,
     6 .8805E+01, .1326E+00, .5964E-03, .6948E-07,-.3100E-08,
     6 .1286E+02, .1918E+00, .8694E-03, .6308E-07,-.8493E-08,
     6 .1881E+02, .2766E+00, .1203E-02, .1540E-06,-.7264E-08,
     6 .2745E+02, .3997E+00, .1704E-02, .9639E-07,-.9601E-08,
     6 .3991E+02, .5767E+00, .2445E-02,-.1334E-06,-.2008E-07,
     6 .5792E+02, .8290E+00, .3401E-02,-.6853E-06,-.2564E-07,
     6 .8363E+02, .1189E+01, .4774E-02,-.1636E-05,-.3909E-07,
     6 .1201E+03, .1697E+01, .6539E-02,-.3083E-05,-.2372E-07,
     6 .1711E+03, .2415E+01, .9242E-02,-.6990E-05,-.6165E-07,
     6 .2423E+03, .3393E+01, .1259E-01,-.1152E-04,-.7322E-07,
     6 .3363E+03, .4678E+01, .1776E-01,-.1541E-04,-.2280E-06,
     6 .4560E+03, .6299E+01, .2211E-01,-.1754E-04, .2881E-07,
     6 .5923E+03, .8286E+01, .3171E-01,-.3971E-04,-.4238E-06,
     6 .7342E+03, .1019E+02, .3986E-01,-.2432E-04,-.4088E-06,
     6 .8676E+03, .1207E+02, .4605E-01,-.3246E-04,-.1283E-06,
     6 .9753E+03, .1367E+02, .5930E-01,-.3330E-04,-.1258E-05 /
      DATA ((( CL8CS(K,J,I), K = 1, 5), J = 1, 6), I = 1, 4) /
     1 .1958E+02,-.3684E+00, .3532E-02,-.3134E-04, .2012E-06,
     1 .2926E+02,-.5510E+00, .5292E-02,-.4673E-04, .2975E-06,
     1 .4370E+02,-.8226E+00, .7933E-02,-.6982E-04, .4380E-06,
     1 .6527E+02,-.1228E+01, .1187E-01,-.1047E-03, .6543E-06,
     1 .9733E+02,-.1834E+01, .1790E-01,-.1567E-03, .9477E-06,
     1 .1451E+03,-.2736E+01, .2682E-01,-.2359E-03, .1418E-05,
     2 .3586E+02,-.6042E+00, .5778E-02,-.4229E-04, .1921E-06,
     2 .5342E+02,-.9044E+00, .8585E-02,-.6263E-04, .2978E-06,
     2 .7956E+02,-.1359E+01, .1262E-01,-.9174E-04, .4884E-06,
     2 .1182E+03,-.2029E+01, .1914E-01,-.1372E-03, .6773E-06,
     2 .1757E+03,-.3030E+01, .2900E-01,-.2041E-03, .9247E-06,
     2 .2602E+03,-.4514E+01, .4452E-01,-.3055E-03, .1167E-05,
     3 .5001E+02,-.8113E+00, .6838E-02,-.4663E-04, .2328E-06,
     3 .7453E+02,-.1204E+01, .1025E-01,-.7046E-04, .3336E-06,
     3 .1113E+03,-.1778E+01, .1539E-01,-.1082E-03, .4612E-06,
     3 .1660E+03,-.2640E+01, .2245E-01,-.1629E-03, .7655E-06,
     3 .2476E+03,-.3897E+01, .3187E-01,-.2468E-03, .1418E-05,
     3 .3716E+03,-.5824E+01, .4584E-01,-.3619E-03, .2344E-05,
     4 .6183E+02,-.9639E+00, .7972E-02,-.5661E-04, .2868E-06,
     4 .9171E+02,-.1422E+01, .1168E-01,-.8534E-04, .4554E-06,
     4 .1360E+03,-.2116E+01, .1763E-01,-.1259E-03, .6335E-06,
     4 .2019E+03,-.3151E+01, .2653E-01,-.1855E-03, .9002E-06,
     4 .2994E+03,-.4715E+01, .4125E-01,-.2715E-03, .1029E-05,
     4 .4415E+03,-.7029E+01, .6260E-01,-.4001E-03, .1380E-05 /
      DATA ((( CL8CF(K,J,I), K = 1, 5), J = 1, 6), I = 1, 4) /
     1 .7421E-01,-.3174E-03,-.5260E-06,-.6845E-08, .2739E-09,
     1 .1109E+00,-.4750E-03,-.7115E-06,-.9775E-08, .3891E-09,
     1 .1654E+00,-.7022E-03,-.8033E-06,-.1512E-07, .5245E-09,
     1 .2468E+00,-.1038E-02,-.1078E-05,-.2435E-07, .7631E-09,
     1 .3669E+00,-.1544E-02,-.2402E-06,-.3858E-07, .8937E-09,
     1 .5456E+00,-.2285E-02, .4113E-06,-.6651E-07, .1243E-08,
     2 .1866E+00,-.7713E-03, .5633E-05,-.1263E-07,-.1169E-09,
     2 .2774E+00,-.1212E-02, .7367E-05,-.7826E-08, .5910E-10,
     2 .4119E+00,-.1990E-02, .7048E-05, .2037E-07, .8969E-09,
     2 .6094E+00,-.2981E-02, .1407E-04, .2701E-07, .6895E-09,
     2 .9005E+00,-.4492E-02, .2643E-04, .3653E-07,-.9039E-10,
     2 .1322E+01,-.6669E-02, .5356E-04, .3722E-07,-.2687E-08,
     3 .3532E+00,-.2106E-02, .5403E-05,-.5615E-08, .2486E-09,
     3 .5256E+00,-.3043E-02, .9467E-05,-.2657E-07, .4013E-10,
     3 .7843E+00,-.4225E-02, .1746E-04,-.9593E-07,-.8149E-09,
     3 .1169E+01,-.6121E-02, .1969E-04,-.1723E-06,-.3055E-09,
     3 .1741E+01,-.8663E-02, .1045E-04,-.2977E-06, .2924E-08,
     3 .2618E+01,-.1294E-01,-.2893E-05,-.3646E-06, .6870E-08,
     4 .5280E+00,-.3273E-02, .1221E-04,-.8300E-07, .5134E-09,
     4 .7823E+00,-.4740E-02, .1606E-04,-.1407E-06, .1188E-08,
     4 .1158E+01,-.7080E-02, .2810E-04,-.1972E-06, .1164E-08,
     4 .1717E+01,-.1065E-01, .4674E-04,-.2611E-06, .1267E-08,
     4 .2543E+01,-.1624E-01, .9264E-04,-.3372E-06,-.1986E-08,
     4 .3728E+01,-.2473E-01, .1575E-03,-.4391E-06,-.4979E-08 /
C
      COMMON /BANDL8GH/ GWL8GH(3), CL8H2OGH(5,28,3)
      DATA GWL8GH / .01500, .00062, .00018 /
      DATA ((( CL8H2OGH(K,J,I), K = 1, 5), J = 1, 28), I = 1, 3) /
     1 .1238E+03, .3729E+01, .4546E-01, .1800E-03,-.5574E-06,
     1 .1259E+03, .3788E+01, .4616E-01, .1826E-03,-.5685E-06,
     1 .1281E+03, .3850E+01, .4689E-01, .1850E-03,-.5824E-06,
     1 .1308E+03, .3919E+01, .4766E-01, .1871E-03,-.6020E-06,
     1 .1344E+03, .4001E+01, .4850E-01, .1886E-03,-.6330E-06,
     1 .1398E+03, .4109E+01, .4949E-01, .1891E-03,-.6805E-06,
     1 .1489E+03, .4261E+01, .5046E-01, .1896E-03,-.7030E-06,
     1 .1657E+03, .4519E+01, .5168E-01, .1851E-03,-.7526E-06,
     1 .1971E+03, .4993E+01, .5469E-01, .1743E-03,-.1071E-05,
     1 .2595E+03, .5883E+01, .5851E-01, .1511E-03,-.1484E-05,
     1 .4676E+03, .8664E+01, .6072E-01, .1174E-03,-.6198E-06,
     1 .6148E+03, .1073E+02, .6252E-01, .8393E-04, .6668E-08,
     1 .8135E+03, .1345E+02, .6976E-01, .8291E-04, .3736E-06,
     1 .1074E+04, .1702E+02, .8706E-01, .9267E-04,-.4359E-06,
     1 .1424E+04, .2180E+02, .1051E+00, .6227E-04,-.9664E-06,
     1 .1872E+04, .2771E+02, .1266E+00, .1050E-03,-.5847E-06,
     1 .2434E+04, .3583E+02, .1612E+00,-.4468E-04,-.2449E-05,
     1 .3117E+04, .4481E+02, .1813E+00,-.9744E-04,-.1081E-05,
     1 .3859E+04, .5436E+02, .2064E+00,-.7528E-04, .1488E-06,
     1 .4619E+04, .6568E+02, .2608E+00,-.2598E-03,-.2915E-05,
     1 .5327E+04, .7536E+02, .3107E+00,-.1966E-03,-.4716E-05,
     1 .5951E+04, .8367E+02, .3174E+00,-.2324E-03,-.1596E-05,
     1 .6375E+04, .9020E+02, .3439E+00,-.3226E-03,-.2788E-05,
     1 .6576E+04, .9363E+02, .3532E+00,-.2514E-03,-.1234E-06,
     1 .6512E+04, .9342E+02, .3676E+00,-.4642E-04, .1030E-05,
     1 .6245E+04, .9136E+02, .4084E+00, .2778E-05,-.4842E-05,
     1 .5874E+04, .8886E+02, .4384E+00,-.5382E-04,-.9002E-05,
     1 .5436E+04, .8482E+02, .4088E+00,-.7221E-04,-.3601E-05,
     2 .4196E+05, .7703E+03, .4408E+01, .2918E-02,-.5111E-04,
     2 .4255E+05, .7812E+03, .4468E+01, .2964E-02,-.5162E-04,
     2 .4316E+05, .7920E+03, .4527E+01, .3016E-02,-.5174E-04,
     2 .4377E+05, .8029E+03, .4582E+01, .3080E-02,-.5133E-04,
     2 .4441E+05, .8137E+03, .4630E+01, .3168E-02,-.4987E-04,
     2 .4510E+05, .8246E+03, .4665E+01, .3300E-02,-.4627E-04,
     2 .4589E+05, .8362E+03, .4678E+01, .3467E-02,-.3965E-04,
     2 .4689E+05, .8471E+03, .4644E+01, .3832E-02,-.2528E-04,
     2 .4834E+05, .8587E+03, .4478E+01, .4506E-02, .1044E-04,
     2 .5064E+05, .8720E+03, .4122E+01, .5793E-02, .7756E-04,
     2 .5617E+05, .9023E+03, .4322E+01, .8086E-02, .4080E-04,
     2 .5969E+05, .9258E+03, .4981E+01, .1003E-01,-.6854E-04,
     2 .6446E+05, .9628E+03, .5034E+01, .1001E-01,-.7544E-04,
     2 .7071E+05, .1019E+04, .4418E+01, .7884E-02, .1396E-04,
     2 .7731E+05, .1091E+04, .4082E+01, .4058E-02, .5513E-04,
     2 .8343E+05, .1160E+04, .3898E+01,-.4203E-03, .5483E-04,
     2 .8752E+05, .1227E+04, .3976E+01,-.5464E-02, .4088E-04,
     2 .8794E+05, .1239E+04, .4270E+01,-.8613E-02,-.5020E-04,
     2 .8597E+05, .1240E+04, .4687E+01,-.1000E-01,-.1007E-03,
     2 .8158E+05, .1165E+04, .4020E+01,-.6421E-03, .6572E-04,
     2 .7379E+05, .1068E+04, .3222E+01,-.2272E-03, .1897E-03,
     2 .6323E+05, .9605E+03, .4888E+01, .4098E-02,-.5693E-04,
     2 .5446E+05, .8831E+03, .5626E+01, .3937E-02,-.1589E-03,
     2 .4699E+05, .8218E+03, .5864E+01, .1437E-02,-.2067E-03,
     2 .4138E+05, .7694E+03, .5271E+01,-.4315E-02,-.1818E-03,
     2 .3608E+05, .6972E+03, .4007E+01,-.6729E-02,-.3734E-04,
     2 .3122E+05, .5862E+03, .2734E+01,-.5396E-02, .4385E-04,
     2 .2676E+05, .4905E+03, .2289E+01,-.5190E-02,-.1373E-04,
     3 .5876E+06, .8786E+04, .3326E+02,-.5493E-01,-.2312E-03,
     3 .5953E+06, .8901E+04, .3370E+02,-.5574E-01,-.2357E-03,
     3 .6029E+06, .9017E+04, .3415E+02,-.5639E-01,-.2395E-03,
     3 .6106E+06, .9132E+04, .3461E+02,-.5719E-01,-.2462E-03,
     3 .6182E+06, .9248E+04, .3508E+02,-.5798E-01,-.2547E-03,
     3 .6255E+06, .9355E+04, .3559E+02,-.5890E-01,-.2690E-03,
     3 .6325E+06, .9471E+04, .3619E+02,-.6002E-01,-.2983E-03,
     3 .6386E+06, .9570E+04, .3695E+02,-.6154E-01,-.3523E-03,
     3 .6430E+06, .9668E+04, .3804E+02,-.6379E-01,-.4585E-03,
     3 .6441E+06, .9741E+04, .3977E+02,-.6758E-01,-.6683E-03,
     3 .6324E+06, .9675E+04, .3994E+02,-.7624E-01,-.6359E-03,
     3 .6275E+06, .9648E+04, .3841E+02,-.8242E-01,-.3230E-03,
     3 .6164E+06, .9548E+04, .3859E+02,-.8507E-01,-.3512E-03,
     3 .5966E+06, .9336E+04, .4030E+02,-.8175E-01,-.6161E-03,
     3 .5674E+06, .8972E+04, .4104E+02,-.6873E-01,-.7410E-03,
     3 .5276E+06, .8485E+04, .4094E+02,-.5921E-01,-.7990E-03,
     3 .4773E+06, .7741E+04, .3869E+02,-.2918E-01,-.6078E-03,
     3 .4207E+06, .7007E+04, .3671E+02,-.1289E-01,-.4289E-03,
     3 .3608E+06, .6161E+04, .3337E+02,-.8575E-02,-.3672E-03,
     3 .3035E+06, .5357E+04, .3087E+02,-.2419E-01,-.6899E-03,
     3 .2577E+06, .4701E+04, .2882E+02,-.3025E-01,-.9580E-03,
     3 .2206E+06, .4102E+04, .2169E+02,-.4117E-01,-.3724E-03,
     3 .1843E+06, .3408E+04, .1533E+02,-.3191E-01, .1304E-03,
     3 .1503E+06, .2756E+04, .1227E+02,-.2889E-01,-.2608E-04,
     3 .1215E+06, .2193E+04, .9605E+01,-.1932E-01,-.6223E-04,
     3 .9916E+05, .1780E+04, .7285E+01,-.1978E-01,-.7329E-04,
     3 .7702E+05, .1408E+04, .6071E+01,-.1547E-01,-.8183E-04,
     3 .5669E+05, .1049E+04, .4697E+01,-.1063E-01,-.6367E-04 /
C
C----------------------------------------------------------------------C
C     SPECTRAL REGION: 0-340 CM^-1. GAS: H2O                           C
C----------------------------------------------------------------------C
C
      COMMON /BANDL9/ GWL9(4), CL9H2O(5,18,4), CL9CS(5,6,4),
     1                         CL9CF(5,6,4)
      DATA GWL9 / .31030, .21554, .19500, .14000 /
      DATA ((( CL9H2O(K,J,I), K = 1, 5), J = 1, 18), I = 1, 4) /
     1 .4226E-01, .6101E-04,-.5452E-06, .2917E-08,-.8187E-11,
     1 .6312E-01, .9117E-04,-.8123E-06, .4352E-08,-.1261E-10,
     1 .9429E-01, .1361E-03,-.1215E-05, .6531E-08,-.1799E-10,
     1 .1409E+00, .2025E-03,-.1819E-05, .9890E-08,-.2549E-10,
     1 .2103E+00, .3016E-03,-.2717E-05, .1489E-07,-.3764E-10,
     1 .3140E+00, .4486E-03,-.4078E-05, .2248E-07,-.5288E-10,
     1 .4685E+00, .6674E-03,-.6145E-05, .3394E-07,-.7046E-10,
     1 .6987E+00, .9916E-03,-.9316E-05, .5116E-07,-.8603E-10,
     1 .1041E+01, .1477E-02,-.1405E-04, .7613E-07,-.1176E-09,
     1 .1551E+01, .2199E-02,-.2109E-04, .1127E-06,-.1699E-09,
     1 .2305E+01, .3266E-02,-.3159E-04, .1670E-06,-.2355E-09,
     1 .3424E+01, .4807E-02,-.4743E-04, .2530E-06,-.3121E-09,
     1 .5076E+01, .7066E-02,-.7225E-04, .3774E-06,-.1824E-09,
     1 .7508E+01, .1036E-01,-.1111E-03, .5618E-06, .2550E-09,
     1 .1108E+02, .1516E-01,-.1650E-03, .8310E-06, .2776E-09,
     1 .1631E+02, .2183E-01,-.2425E-03, .1277E-05,-.4663E-09,
     1 .2397E+02, .3125E-01,-.3851E-03, .1953E-05, .2786E-08,
     1 .3504E+02, .4459E-01,-.6395E-03, .2839E-05, .1592E-07,
     2 .2546E+00, .1815E-03,-.2208E-05, .2111E-08,-.1561E-10,
     2 .3800E+00, .2695E-03,-.3286E-05, .3314E-08,-.2203E-10,
     2 .5673E+00, .3995E-03,-.4889E-05, .5241E-08,-.3276E-10,
     2 .8468E+00, .5921E-03,-.7336E-05, .8258E-08,-.4111E-10,
     2 .1263E+01, .8762E-03,-.1088E-04, .1268E-07,-.7659E-10,
     2 .1883E+01, .1300E-02,-.1613E-04, .1870E-07,-.1449E-09,
     2 .2807E+01, .1937E-02,-.2426E-04, .2562E-07,-.1964E-09,
     2 .4180E+01, .2891E-02,-.3614E-04, .3357E-07,-.3073E-09,
     2 .6216E+01, .4325E-02,-.5375E-04, .4035E-07,-.4175E-09,
     2 .9238E+01, .6420E-02,-.7875E-04, .4805E-07,-.7080E-09,
     2 .1369E+02, .9538E-02,-.1166E-03, .4148E-07,-.9143E-09,
     2 .2023E+02, .1414E-01,-.1726E-03, .3111E-09,-.8711E-09,
     2 .2978E+02, .2082E-01,-.2275E-03,-.1247E-06,-.4911E-08,
     2 .4370E+02, .2936E-01,-.2826E-03,-.1437E-06,-.1470E-07,
     2 .6396E+02, .3982E-01,-.4087E-03, .1061E-06,-.2130E-07,
     2 .9322E+02, .5274E-01,-.6560E-03, .3482E-06,-.1457E-07,
     2 .1352E+03, .6623E-01,-.1059E-02, .6230E-06,-.2214E-08,
     2 .1954E+03, .9042E-01,-.1357E-02,-.3471E-06,-.3964E-07,
     3 .9075E+00, .2349E-03,-.4012E-05, .1157E-07, .7962E-10,
     3 .1352E+01, .3354E-03,-.6038E-05, .1817E-07, .1270E-09,
     3 .2015E+01, .4793E-03,-.9052E-05, .2765E-07, .1831E-09,
     3 .3002E+01, .6835E-03,-.1364E-04, .4185E-07, .2676E-09,
     3 .4472E+01, .9851E-03,-.2063E-04, .6109E-07, .4010E-09,
     3 .6654E+01, .1432E-02,-.3124E-04, .8690E-07, .6321E-09,
     3 .9898E+01, .2095E-02,-.4652E-04, .1206E-06, .9307E-09,
     3 .1469E+02, .2980E-02,-.6872E-04, .1737E-06, .1423E-08,
     3 .2176E+02, .4122E-02,-.9926E-04, .2522E-06, .1794E-08,
     3 .3212E+02, .5504E-02,-.1445E-03, .3404E-06, .3016E-08,
     3 .4721E+02, .6684E-02,-.1736E-03, .4461E-06,-.7950E-09,
     3 .6916E+02, .4844E-02,-.1773E-03, .1003E-05,-.1243E-07,
     3 .1010E+03,-.3009E-02,-.2229E-03, .2370E-05,-.2420E-07,
     3 .1471E+03,-.2144E-01,-.3871E-03, .4663E-05,-.2618E-07,
     3 .2134E+03,-.5700E-01,-.8083E-03, .7471E-05,-.2514E-08,
     3 .3089E+03,-.1087E+00,-.1370E-02, .9538E-05, .4948E-08,
     3 .4442E+03,-.1908E+00,-.1056E-02, .1380E-04,-.1654E-06,
     3 .6328E+03,-.3826E+00,-.1620E-02, .3299E-04,-.2052E-06,
     4 .3629E+01,-.7787E-03,-.6976E-05, .1193E-06,-.5476E-09,
     4 .5397E+01,-.1244E-02,-.1101E-04, .1756E-06,-.8064E-09,
     4 .8026E+01,-.1968E-02,-.1715E-04, .2604E-06,-.1219E-08,
     4 .1193E+02,-.2999E-02,-.2528E-04, .3702E-06,-.2060E-08,
     4 .1771E+02,-.4487E-02,-.4061E-04, .5268E-06,-.2565E-08,
     4 .2627E+02,-.6705E-02,-.6518E-04, .7423E-06,-.2699E-08,
     4 .3884E+02,-.1013E-01,-.8749E-04, .1037E-05,-.5258E-08,
     4 .5724E+02,-.1537E-01,-.1224E-03, .1402E-05,-.7379E-08,
     4 .8397E+02,-.2451E-01,-.1376E-03, .1964E-05,-.1405E-07,
     4 .1227E+03,-.4169E-01,-.1280E-03, .3161E-05,-.2905E-07,
     4 .1785E+03,-.7371E-01,-.1371E-03, .5716E-05,-.5114E-07,
     4 .2589E+03,-.1257E+00,-.4193E-03, .9398E-05,-.4654E-07,
     4 .3743E+03,-.2045E+00,-.1390E-02, .1334E-04, .5155E-07,
     4 .5382E+03,-.2915E+00,-.2555E-02, .1103E-04, .1227E-06,
     4 .7683E+03,-.4207E+00,-.1761E-02, .1120E-04,-.1680E-06,
     4 .1087E+04,-.7028E+00,-.2430E-02, .3252E-04,-.2404E-06,
     4 .1518E+04,-.1012E+01,-.7249E-02, .3761E-04, .1481E-06,
     4 .2063E+04,-.1267E+01,-.7527E-02, .1644E-04,-.2221E-06 /
      DATA ((( CL9CS(K,J,I), K = 1, 5), J = 1, 6), I = 1, 4) /
     1 .1166E+03,-.1766E+01, .1400E-01,-.9902E-04, .5646E-06,
     1 .1730E+03,-.2609E+01, .2048E-01,-.1482E-03, .8792E-06,
     1 .2571E+03,-.3825E+01, .2895E-01,-.2303E-03, .1611E-05,
     1 .3826E+03,-.5674E+01, .4231E-01,-.3452E-03, .2522E-05,
     1 .5680E+03,-.8350E+01, .6093E-01,-.5280E-03, .4119E-05,
     1 .8396E+03,-.1230E+02, .8937E-01,-.7953E-03, .6324E-05,
     2 .2341E+03,-.4082E+01, .3964E-01,-.3025E-03, .1407E-05,
     2 .3497E+03,-.6110E+01, .5917E-01,-.4487E-03, .2090E-05,
     2 .5217E+03,-.9141E+01, .8852E-01,-.6579E-03, .2952E-05,
     2 .7755E+03,-.1357E+02, .1309E+00,-.9772E-03, .4503E-05,
     2 .1152E+04,-.2029E+02, .1982E+00,-.1404E-02, .5673E-05,
     2 .1716E+04,-.3034E+02, .3007E+00,-.2035E-02, .6973E-05,
     3 .2553E+03,-.4445E+01, .4359E-01,-.4011E-03, .2437E-05,
     3 .3800E+03,-.6617E+01, .6517E-01,-.6000E-03, .3624E-05,
     3 .5659E+03,-.9905E+01, .9936E-01,-.8851E-03, .4962E-05,
     3 .8467E+03,-.1487E+02, .1509E+00,-.1302E-02, .6830E-05,
     3 .1269E+04,-.2221E+02, .2221E+00,-.1981E-02, .1118E-04,
     3 .1898E+04,-.3310E+02, .3270E+00,-.2980E-02, .1758E-04,
     4 .2682E+03,-.4700E+01, .4697E-01,-.4274E-03, .2521E-05,
     4 .4004E+03,-.7039E+01, .7080E-01,-.6350E-03, .3636E-05,
     4 .5980E+03,-.1053E+02, .1062E+00,-.9540E-03, .5472E-05,
     4 .8897E+03,-.1567E+02, .1591E+00,-.1442E-02, .8254E-05,
     4 .1319E+04,-.2341E+02, .2417E+00,-.2111E-02, .1116E-04,
     4 .1955E+04,-.3458E+02, .3462E+00,-.3158E-02, .1903E-04 /
      DATA ((( CL9CF(K,J,I), K = 1, 5), J = 1, 6), I = 1, 4) /
     1 .3948E+01,-.1931E-01, .3552E-04,-.5027E-06, .9660E-08,
     1 .5835E+01,-.2758E-01, .3790E-04,-.8682E-06, .1676E-07,
     1 .8634E+01,-.3825E-01,-.2193E-04,-.1821E-05, .4119E-07,
     1 .1281E+02,-.5589E-01,-.8400E-04,-.2935E-05, .6985E-07,
     1 .1896E+02,-.7819E-01,-.2426E-03,-.5361E-05, .1255E-06,
     1 .2779E+02,-.1116E+00,-.3913E-03,-.8869E-05, .1948E-06,
     2 .1200E+02,-.9434E-01, .6079E-03,-.1841E-05,-.1235E-07,
     2 .1796E+02,-.1419E+00, .8872E-03,-.2595E-05,-.1581E-07,
     2 .2682E+02,-.2126E+00, .1297E-02,-.3450E-05,-.2306E-07,
     2 .3985E+02,-.3150E+00, .1836E-02,-.5104E-05,-.1777E-07,
     2 .5911E+02,-.4745E+00, .2918E-02,-.5347E-05,-.7086E-07,
     2 .8774E+02,-.7101E+00, .4880E-02,-.5460E-05,-.2140E-06,
     3 .1354E+02,-.1043E+00, .5425E-03,-.4045E-05, .2942E-07,
     3 .2015E+02,-.1553E+00, .8368E-03,-.6041E-05, .4049E-07,
     3 .2999E+02,-.2340E+00, .1429E-02,-.8531E-05, .2462E-07,
     3 .4484E+02,-.3526E+00, .2331E-02,-.1189E-04,-.4458E-08,
     3 .6743E+02,-.5260E+00, .3208E-02,-.1884E-04, .5220E-07,
     3 .1012E+03,-.7764E+00, .4308E-02,-.3073E-04, .1633E-06,
     4 .1442E+02,-.1125E+00, .6521E-03,-.3817E-05, .1981E-07,
     4 .2152E+02,-.1690E+00, .1026E-02,-.5605E-05, .2030E-07,
     4 .3214E+02,-.2538E+00, .1542E-02,-.8381E-05, .3521E-07,
     4 .4779E+02,-.3767E+00, .2389E-02,-.1346E-04, .4555E-07,
     4 .7062E+02,-.5673E+00, .4024E-02,-.1808E-04,-.2171E-07,
     4 .1047E+03,-.8481E+00, .5019E-02,-.2529E-04, .1479E-06 /
C
      COMMON /BANDL9GH/ GWL9GH(6), CL9H2OGH(5,28,6)
      DATA GWL9GH / .08000, .05000, .00800, .00094, .00019, .00003 /
      DATA ((( CL9H2OGH(K,J,I), K = 1, 5), J = 1, 28), I = 1, 3) /
     5 .1189E-01, .8033E-05,-.6385E-07, .2642E-09, .3063E-11,
     5 .2393E-01, .1606E-04,-.1259E-06, .3615E-09, .5371E-11,
     5 .4807E-01, .3218E-04,-.1752E-06, .2064E-09,-.2419E-11,
     5 .9646E-01, .6140E-04,-.1156E-06, .5331E-10,-.4231E-10,
     5 .1936E+00, .1096E-03,-.1202E-06, .9577E-09,-.8929E-10,
     5 .3887E+00, .1807E-03,-.2764E-06, .5614E-08,-.1281E-09,
     5 .7812E+00, .2805E-03,-.9748E-06, .1664E-07,-.1473E-09,
     5 .1571E+01, .4071E-03,-.2464E-05, .3803E-07,-.2193E-09,
     5 .3160E+01, .4679E-03,-.5703E-05, .8494E-07,-.4061E-09,
     5 .6356E+01, .3292E-03,-.1362E-04, .1664E-06,-.7436E-09,
     5 .1826E+02,-.9463E-03,-.4659E-04, .4466E-06,-.2670E-08,
     5 .2705E+02,-.1970E-02,-.7248E-04, .6245E-06,-.3709E-08,
     5 .4005E+02,-.3823E-02,-.1166E-03, .8934E-06,-.4207E-08,
     5 .5918E+02,-.6456E-02,-.1708E-03, .1134E-05,-.6114E-08,
     5 .8717E+02,-.1031E-01,-.2243E-03, .1313E-05,-.1105E-07,
     5 .1280E+03,-.1973E-01,-.2872E-03, .1942E-05,-.1862E-07,
     5 .1871E+03,-.3641E-01,-.3462E-03, .2902E-05,-.3248E-07,
     5 .2722E+03,-.7492E-01,-.4067E-03, .5809E-05,-.5395E-07,
     5 .3942E+03,-.1439E+00,-.5441E-03, .1108E-04,-.9498E-07,
     5 .5692E+03,-.2606E+00,-.1439E-02, .1907E-04,-.4731E-07,
     5 .8190E+03,-.4329E+00,-.3509E-02, .2691E-04, .1211E-06,
     5 .1172E+04,-.6517E+00,-.4980E-02, .2611E-04, .1083E-06,
     5 .1665E+04,-.1030E+01,-.3794E-02, .4164E-04,-.4559E-06,
     5 .2345E+04,-.1756E+01,-.7346E-02, .9487E-04,-.4119E-06,
     5 .3246E+04,-.2487E+01,-.1621E-01, .1023E-03, .1653E-06,
     5 .4341E+04,-.3216E+01,-.1522E-01, .7815E-04,-.5541E-06,
     5 .5545E+04,-.4277E+01,-.7660E-02, .1870E-03,-.2722E-05,
     5 .6763E+04,-.4141E+01,-.4154E-01, .8076E-04, .2662E-05,
     6 .4305E+00, .3161E-02, .1746E-04, .9764E-07, .3996E-09,
     6 .7599E+00, .4495E-02, .2078E-04, .1243E-06, .4571E-09,
     6 .1363E+01, .6564E-02, .2335E-04, .1351E-06, .6532E-09,
     6 .2479E+01, .9605E-02, .2449E-04, .1530E-06, .1117E-08,
     6 .4566E+01, .1414E-01, .2753E-04, .1936E-06, .9333E-09,
     6 .8525E+01, .2076E-01, .2243E-04, .2845E-06, .1622E-08,
     6 .1607E+02, .3052E-01, .1895E-04, .3518E-06,-.1335E-08,
     6 .3061E+02, .4421E-01,-.4170E-05, .6186E-06,-.5155E-08,
     6 .5867E+02, .6587E-01,-.1254E-03, .2684E-06, .1537E-08,
     6 .1131E+03, .9347E-01,-.4256E-03,-.3949E-06, .2888E-07,
     6 .3032E+03, .1242E+00,-.2503E-03, .8244E-06,-.6356E-07,
     6 .4377E+03, .1122E+00,-.1440E-03, .4225E-05,-.1285E-06,
     6 .6310E+03, .3900E-01,-.3976E-03, .1431E-04,-.1646E-06,
     6 .9078E+03,-.8215E-01,-.1018E-02, .2577E-04,-.1927E-06,
     6 .1307E+04,-.2744E+00,-.2323E-02, .4191E-04,-.2083E-06,
     6 .1880E+04,-.5997E+00,-.5883E-02, .6281E-04,-.5534E-07,
     6 .2697E+04,-.9574E+00,-.1055E-01, .6398E-04, .7040E-07,
     6 .3849E+04,-.1573E+01,-.1371E-01, .9501E-04,-.3864E-06,
     6 .5464E+04,-.2523E+01,-.2299E-01, .1431E-03,-.2625E-06,
     6 .7645E+04,-.3584E+01,-.3692E-01, .1610E-03,-.2000E-06,
     6 .1044E+05,-.5018E+01,-.3813E-01, .1942E-03,-.2101E-05,
     6 .1379E+05,-.7395E+01,-.5262E-01, .4428E-03,-.3657E-05,
     6 .1752E+05,-.7635E+01,-.1170E+00, .2848E-03, .3984E-05,
     6 .2124E+05,-.8251E+01,-.1104E+00, .1503E-03, .2151E-06,
     6 .2453E+05,-.1158E+02,-.8440E-01, .2243E-03,-.3439E-05,
     6 .2707E+05,-.1349E+02,-.1654E+00,-.1111E-04, .1140E-04,
     6 .2831E+05,-.1490E+02,-.7166E-01,-.1424E-03, .6813E-06,
     6 .2847E+05,-.1807E+02,-.2165E-01,-.2007E-03,-.2902E-05,
     7 .6099E+03, .7817E+01, .4070E-01, .1169E-03,-.3335E-06,
     7 .6402E+03, .7733E+01, .4121E-01, .1208E-03,-.3984E-06,
     7 .6962E+03, .7564E+01, .4096E-01, .1317E-03,-.4466E-06,
     7 .8011E+03, .7361E+01, .3831E-01, .1482E-03,-.4137E-06,
     7 .9790E+03, .7284E+01, .3362E-01, .1448E-03,-.1560E-06,
     7 .1278E+04, .7135E+01, .2383E-01, .1754E-03, .1057E-05,
     7 .1779E+04, .7101E+01, .2331E-01, .1808E-03, .5076E-06,
     7 .2630E+04, .7231E+01, .1379E-01, .2154E-03, .2451E-05,
     7 .4099E+04, .7726E+01,-.9201E-02, .1460E-03, .5437E-05,
     7 .6644E+04, .7222E+01,-.7209E-02, .2517E-03, .2655E-05,
     7 .1469E+05, .4748E+01, .4023E-01, .4281E-03,-.1120E-04,
     7 .2021E+05, .1874E+01, .1824E-02, .6068E-03,-.7271E-05,
     7 .2809E+05,-.1604E+01,-.9092E-01, .7971E-03, .2654E-05,
     7 .3901E+05,-.8209E+01,-.1128E+00, .1302E-02,-.6675E-05,
     7 .5371E+05,-.1471E+02,-.1666E+00, .1328E-02,-.1126E-04,
     7 .7237E+05,-.3177E+02,-.2566E+00, .2354E-02,-.1280E-04,
     7 .9439E+05,-.4523E+02,-.3673E+00, .3028E-02,-.2408E-04,
     7 .1185E+06,-.4731E+02,-.8933E+00, .1770E-02, .5301E-04,
     7 .1421E+06,-.6197E+02,-.6393E+00, .2076E-02,-.7348E-05,
     7 .1627E+06,-.8442E+02,-.6084E+00, .2413E-02,-.8757E-05,
     7 .1778E+06,-.9687E+02,-.8941E+00, .1762E-02, .4770E-04,
     7 .1849E+06,-.9453E+02,-.1786E+00, .4813E-03,-.4317E-04,
     7 .1853E+06,-.9769E+02,-.2466E-01,-.6269E-03,-.4343E-04,
     7 .1787E+06,-.1168E+03,-.4092E+00, .2030E-03, .5138E-04,
     7 .1655E+06,-.1111E+03,-.2212E+00, .1342E-04, .4497E-04,
     7 .1475E+06,-.8839E+02, .7246E+00,-.7165E-03,-.7844E-04,
     7 .1295E+06,-.5339E+02, .2270E+00,-.3102E-02, .1701E-04,
     7 .1112E+06,-.2827E+02, .2234E+00,-.1780E-02,-.5194E-05 /
      DATA ((( CL9H2OGH(K,J,I), K = 1, 5), J = 1, 28), I = 4, 6) /
     8 .2282E+06, .2322E+04, .5714E+01,-.3761E-01,-.3151E-03,
     8 .2296E+06, .2334E+04, .5747E+01,-.3775E-01,-.3161E-03,
     8 .2312E+06, .2345E+04, .5791E+01,-.3783E-01,-.3171E-03,
     8 .2331E+06, .2352E+04, .5865E+01,-.3775E-01,-.3188E-03,
     8 .2359E+06, .2350E+04, .5975E+01,-.3706E-01,-.3216E-03,
     8 .2410E+06, .2333E+04, .6160E+01,-.3621E-01,-.3337E-03,
     8 .2509E+06, .2300E+04, .5864E+01,-.3760E-01,-.2564E-03,
     8 .2712E+06, .2210E+04, .5542E+01,-.3522E-01,-.2389E-03,
     8 .3125E+06, .2047E+04, .3513E+01,-.2898E-01,-.3194E-04,
     8 .3887E+06, .1756E+04, .2574E+01,-.1444E-01,-.1193E-03,
     8 .6162E+06, .1170E+04,-.1032E+01, .1116E-02,-.1318E-03,
     8 .7345E+06, .8614E+03, .1080E+00, .2579E-01,-.4213E-03,
     8 .8733E+06, .6992E+03,-.5564E+01, .4154E-02, .4235E-03,
     8 .1012E+07, .4598E+03,-.6980E+01,-.5912E-02, .6048E-03,
     8 .1131E+07,-.5627E+02,-.4413E+01, .2644E-01, .2096E-03,
     8 .1228E+07,-.2410E+03,-.8442E+01, .3193E-02, .7231E-03,
     8 .1279E+07,-.4451E+03,-.4661E+01, .6349E-02, .2112E-03,
     8 .1276E+07,-.6033E+03, .3384E+01, .1265E-01,-.9195E-03,
     8 .1234E+07,-.4116E+03,-.2972E+01,-.3335E-01, .4665E-03,
     8 .1148E+07,-.4610E+03,-.4323E+01,-.2588E-01, .7519E-03,
     8 .1026E+07,-.4768E+03, .7506E+00,-.1789E-01, .9226E-05,
     8 .8869E+06,-.3862E+03,-.2371E+00,-.2820E-01, .3429E-03,
     8 .7349E+06,-.2595E+03, .2788E+01,-.2175E-01,-.3169E-04,
     8 .5951E+06, .1444E+01, .5039E+01,-.2911E-01,-.4986E-03,
     8 .4756E+06, .1935E+03, .2978E+01,-.3445E-01,-.2903E-03,
     8 .3819E+06, .1285E+03,-.2575E+01,-.1307E-01, .4107E-03,
     8 .3039E+06,-.3683E+02,-.2901E+01, .1221E-01, .3292E-03,
     8 .2372E+06,-.3331E+02,-.6236E+00, .1030E-01, .7012E-04,
     9 .8912E+07, .1888E+05,-.2018E+03, .2427E+00, .2538E-02,
     9 .8961E+07, .1897E+05,-.2030E+03, .2430E+00, .2589E-02,
     9 .9009E+07, .1904E+05,-.2044E+03, .2421E+00, .2702E-02,
     9 .9058E+07, .1907E+05,-.2059E+03, .2410E+00, .2868E-02,
     9 .9106E+07, .1903E+05,-.2078E+03, .2456E+00, .3077E-02,
     9 .9155E+07, .1883E+05,-.2105E+03, .2553E+00, .3488E-02,
     9 .9210E+07, .1827E+05,-.2120E+03, .2897E+00, .3898E-02,
     9 .9252E+07, .1717E+05,-.2138E+03, .3222E+00, .5106E-02,
     9 .9267E+07, .1509E+05,-.2053E+03, .3561E+00, .6120E-02,
     9 .9146E+07, .1137E+05,-.1615E+03, .4325E+00, .4301E-02,
     9 .8346E+07, .6656E+04,-.4017E+02, .3709E-01, .8325E-05,
     9 .7973E+07, .4554E+04,-.4539E+01,-.1428E+00,-.1672E-02,
     9 .7540E+07, .2223E+04, .3524E+02,-.2751E+00,-.5174E-02,
     9 .7117E+07,-.1772E+03, .5044E+02,-.2482E+00,-.7386E-02,
     9 .6698E+07,-.1730E+04, .1359E+02,-.3353E+00,-.2502E-02,
     9 .6106E+07,-.3671E+04, .1312E+02,-.5376E-01,-.2398E-02,
     9 .5318E+07,-.2869E+04, .2268E+02,-.1105E+00,-.1642E-02,
     9 .4486E+07,-.2202E+04, .1947E+02,-.1024E+00,-.4577E-03,
     9 .3682E+07,-.2441E+04, .2543E+02, .7032E-01,-.2329E-02,
     9 .2939E+07,-.7348E+03, .2495E+02,-.6691E-01,-.2004E-02,
     9 .2317E+07, .1401E+03, .1250E+02,-.1234E+00,-.4546E-03,
     9 .1823E+07,-.5861E+02,-.4110E+01,-.4951E-01, .1516E-02,
     9 .1417E+07,-.5177E+03,-.2440E+01, .2654E-01, .6276E-03,
     9 .1067E+07,-.5562E+03, .7921E+01, .4415E-01,-.6922E-03,
     9 .8017E+06,-.4462E+03, .6757E+01, .3992E-01,-.3595E-03,
     9 .5893E+06,-.2117E+01, .1586E+02,-.1702E-02,-.2042E-02,
     9 .4313E+06, .4256E+03, .1230E+02,-.5361E-01,-.1591E-02,
     9 .3157E+06, .4554E+03, .3390E+01,-.4530E-01,-.1535E-03,
     a .6301E+08,-.3184E+06, .1060E+04,-.2921E+00,-.1112E-01,
     a .6337E+08,-.3201E+06, .1066E+04,-.2934E+00,-.1139E-01,
     a .6368E+08,-.3217E+06, .1074E+04,-.2921E+00,-.1189E-01,
     a .6398E+08,-.3228E+06, .1084E+04,-.3109E+00,-.1286E-01,
     a .6422E+08,-.3233E+06, .1095E+04,-.3827E+00,-.1372E-01,
     a .6427E+08,-.3223E+06, .1111E+04,-.5298E+00,-.1554E-01,
     a .6395E+08,-.3184E+06, .1135E+04,-.7931E+00,-.1905E-01,
     a .6311E+08,-.3087E+06, .1165E+04,-.1223E+01,-.2563E-01,
     a .6104E+08,-.2882E+06, .1178E+04,-.1770E+01,-.3596E-01,
     a .5722E+08,-.2463E+06, .9214E+03,-.2642E+01,-.2202E-01,
     a .4622E+08,-.1555E+06, .1512E+03,-.3459E+00, .4365E-02,
     a .4032E+08,-.1143E+06,-.1324E+03, .3485E-01, .2109E-01,
     a .3365E+08,-.8049E+05,-.2025E+03, .1247E+01, .1580E-01,
     a .2677E+08,-.5099E+05,-.2342E+03, .1098E+01, .2610E-01,
     a .2033E+08,-.2543E+05,-.6168E+02, .5573E+00, .1124E-01,
     a .1517E+08,-.7334E+04, .9473E+02,-.6639E+00,-.3270E-02,
     a .1174E+08,-.4546E+04,-.2513E+02,-.5757E+00, .9035E-02,
     a .9169E+07,-.3423E+04,-.9998E+02,-.5366E+00, .1603E-01,
     a .6818E+07,-.2781E+04,-.4965E+01,-.3910E+00, .2521E-02,
     a .4979E+07,-.4028E+04, .5020E+02, .3666E-01,-.6955E-02,
     a .3584E+07,-.3114E+04, .5436E+02, .2033E+00,-.5940E-02,
     a .2553E+07,-.9871E+03, .3126E+02,-.1372E-01,-.2165E-02,
     a .1803E+07, .2082E+03, .3012E+00,-.1819E+00, .2786E-02,
     a .1270E+07, .6351E+03,-.2304E+02,-.2116E+00, .5869E-02,
     a .8807E+06, .5995E+03,-.2112E+02,-.1497E+00, .5022E-02,
     a .6069E+06, .3525E+03,-.1520E+02,-.3614E-01, .4058E-02,
     a .4196E+06, .4041E+03,-.4181E+01,-.1160E-01, .1868E-02,
     a .2939E+06, .3951E+03, .1880E+01,-.3537E-02, .3566E-03 /
C
      END

      BLOCK DATA CKDSW
C
      IMPLICIT REAL (A-H,O-Z),
     +INTEGER (I-N)
C
C----------------------------------------------------------------------C
C     SPECTRAL REGION: 14500-50000 CM^-1. GAS: H2O, O2                 C
C     K IN UNIT CM^2 / GRAM FOR ALL DATA.                              C
C     NOTE: GW = GW * SOLAR ENERGY FOR ALL SOLAR BANDS                 C
C     CS1O2GH3 FOR O2 ABSORPTANCE COEFFICIENT AT 41200 - 50000 CM^-1   C
C     FOR HERZBERG BAND. THE SCHUMAN-RUNGE BAND NOT INCLUDED           C
C     CS1O21 FOR O2 AT 25000 - 14500 CM^-1     LI 2003                 C
C----------------------------------------------------------------------C
C
      COMMON /BANDS1/ GWS1(6), CS1O3(3,6), CS1O21
      DATA GWS1 / 280.85678, 236.97212, 83.43346, 6.70594, 10.14919,
     1              4.73084 /
C
      DATA ((CS1O3(I,J), I = 1, 3), J = 1, 6) /
     1 .4203E+02,-.2043E-02, .5034E-04,
     2 .6390E+01, .9014E-02, .1214E-03,
     3 .2584E+02, .1504E+00, .1025E-02,
     4 .5269E+03, .1181E+01, .2640E-02,
     5 .4034E+04, .5636E+01, .1300E-01,
     6 .2570E+05, .1570E+02, .8769E-01 /
C
      DATA CS1O21 / 0.511E-04 /
C
      COMMON /BANDS1GH/ GWS1GH(3), CS1O3GH(3,3), CS1O2GH3
      DATA GWS1GH / 3.20501, 3.17570, 1.21100 /
C
      DATA ((CS1O3GH(I,J), I = 1, 3), J = 1, 3) /
     7 .8575E+05, .1976E+02, .1222E-01,
     8 .1254E+06,-.1886E+02,-.2423E+00,
     9 .2228E+05,-.1677E+01, .1736E+00 /
C
      DATA CS1O2GH3 / 0.900007E-01 / 
C
C----------------------------------------------------------------------C
C     SPECTRAL REGION: 8200-14500 CM^-1. GAS: H2O, O2, O3              C
C----------------------------------------------------------------------C
C
      COMMON /BANDS2/ GW(4), CS2H2O(5,18,4), CS2O2(5,18,2), CS2O3(4)
      DATA GW / 292.83546, 62.41776, 53.30647, 21.53219 /
      DATA ( CS2O3(I), I = 1,  4) / .2185E+01, .1517E+01, .1684E+01,
     1                              .2268E+01 /
      DATA ((( CS2H2O(K,J,I), K = 1, 5), J = 1, 18), I = 1, 4) /
     1 .5727E-05, .7913E-08,-.2870E-10, .5559E-13,-.3821E-15,
     1 .8450E-05, .1153E-07,-.4306E-10, .7787E-13,-.7543E-15,
     1 .1247E-04, .1684E-07,-.6614E-10, .1037E-12,-.1100E-14,
     1 .1839E-04, .2444E-07,-.9416E-10, .1830E-12,-.2601E-14,
     1 .2714E-04, .3553E-07,-.1543E-09, .3050E-12,-.1829E-14,
     1 .4005E-04, .5192E-07,-.2335E-09, .4800E-12,-.3429E-14,
     1 .5909E-04, .7574E-07,-.3536E-09, .8039E-12,-.5700E-14,
     1 .8719E-04, .1116E-06,-.5451E-09, .1131E-11,-.7054E-14,
     1 .1287E-03, .1631E-06,-.8644E-09, .1835E-11,-.5946E-14,
     1 .1896E-03, .2405E-06,-.1314E-08, .2609E-11,-.9643E-14,
     1 .2790E-03, .3566E-06,-.2043E-08, .3507E-11,-.3798E-14,
     1 .4088E-03, .5305E-06,-.3089E-08, .4182E-11, .3468E-15,
     1 .5956E-03, .7789E-06,-.4303E-08, .5732E-11,-.4682E-13,
     1 .8593E-03, .1129E-05,-.6199E-08, .8142E-11,-.6072E-13,
     1 .1221E-02, .1632E-05,-.8534E-08, .8854E-11,-.5919E-13,
     1 .1703E-02, .2339E-05,-.1107E-07, .1116E-10,-.8467E-13,
     1 .2318E-02, .3352E-05,-.1481E-07, .1337E-10,-.4540E-13,
     1 .3071E-02, .4751E-05,-.2004E-07, .5499E-11, .1391E-12,
     2 .1812E-03, .5011E-06,-.8218E-09, .2032E-11,-.2913E-14,
     2 .2596E-03, .6953E-06,-.1514E-08, .2176E-11, .3824E-13,
     2 .3718E-03, .9660E-06,-.2442E-08, .1820E-11, .7881E-13,
     2 .5323E-03, .1350E-05,-.3699E-08,-.4671E-12, .1219E-12,
     2 .7617E-03, .1879E-05,-.5219E-08,-.3813E-11, .1363E-12,
     2 .1089E-02, .2593E-05,-.7511E-08,-.7619E-11, .1846E-12,
     2 .1556E-02, .3586E-05,-.1149E-07,-.1861E-10, .3683E-12,
     2 .2217E-02, .4926E-05,-.1452E-07,-.3508E-10, .1626E-12,
     2 .3151E-02, .6624E-05,-.1710E-07,-.4111E-10,-.3080E-12,
     2 .4466E-02, .8753E-05,-.2256E-07,-.3110E-10,-.6199E-12,
     2 .6302E-02, .1146E-04,-.3460E-07,-.1332E-10,-.1316E-12,
     2 .8826E-02, .1518E-04,-.4505E-07,-.1953E-10,-.8579E-13,
     2 .1224E-01, .2030E-04,-.5950E-07,-.2591E-10, .3450E-12,
     2 .1675E-01, .2818E-04,-.8032E-07,-.1297E-09, .1356E-11,
     2 .2259E-01, .3962E-04,-.1029E-06,-.3474E-09, .1504E-11,
     2 .2984E-01, .5403E-04,-.1072E-06,-.5094E-09,-.1288E-11,
     2 .3848E-01, .7229E-04,-.1133E-06,-.7213E-09,-.3025E-11,
     2 .4843E-01, .9131E-04,-.1408E-06,-.5242E-09,-.1241E-11,
     3 .5547E-02, .4158E-04, .3280E-07,-.2032E-09, .2258E-11,
     3 .6413E-02, .4296E-04, .2637E-07,-.1630E-09, .1830E-11,
     3 .7564E-02, .4486E-04, .1419E-07,-.1373E-09, .2123E-11,
     3 .9084E-02, .4738E-04,-.1630E-08,-.9672E-10, .3104E-11,
     3 .1109E-01, .5067E-04,-.1787E-07,-.8095E-10, .4110E-11,
     3 .1377E-01, .5523E-04,-.1690E-07,-.1512E-09, .1471E-11,
     3 .1734E-01, .6124E-04,-.2382E-07,-.2526E-09, .5236E-13,
     3 .2213E-01, .6899E-04,-.3343E-07,-.4029E-09,-.1248E-11,
     3 .2858E-01, .7778E-04,-.5925E-07,-.4614E-09,-.8470E-12,
     3 .3720E-01, .8884E-04,-.8387E-07,-.6317E-09,-.1760E-11,
     3 .4869E-01, .1016E-03,-.1577E-06,-.7254E-09, .4759E-11,
     3 .6377E-01, .1176E-03,-.1994E-06,-.9033E-09, .6915E-11,
     3 .8323E-01, .1374E-03,-.1734E-06,-.1207E-08,-.3134E-11,
     3 .1079E+00, .1582E-03,-.1595E-06,-.8371E-09,-.1246E-10,
     3 .1386E+00, .1846E-03,-.2601E-06,-.6825E-10,-.2876E-11,
     3 .1758E+00, .2208E-03,-.3373E-06, .1794E-09,-.8184E-11,
     3 .2194E+00, .2561E-03,-.5008E-06, .1152E-08,-.2540E-13,
     3 .2687E+00, .3031E-03,-.6058E-06, .1136E-08,-.2361E-10,
     4 .2916E+00, .1904E-02, .1167E-05,-.1007E-07, .1006E-10,
     4 .2971E+00, .1880E-02, .1156E-05,-.9691E-08, .1463E-10,
     4 .3052E+00, .1853E-02, .1207E-05,-.9447E-08,-.6028E-12,
     4 .3169E+00, .1820E-02, .1175E-05,-.8540E-08, .5159E-11,
     4 .3333E+00, .1787E-02, .1203E-05,-.9259E-08,-.1682E-11,
     4 .3564E+00, .1746E-02, .1460E-05,-.8860E-08,-.5574E-10,
     4 .3892E+00, .1703E-02, .1277E-05,-.7849E-08,-.5692E-10,
     4 .4342E+00, .1689E-02, .4175E-06,-.1098E-07, .4219E-10,
     4 .4926E+00, .1672E-02, .1925E-06,-.1304E-07, .4876E-10,
     4 .5680E+00, .1621E-02, .3194E-07,-.9251E-08, .2122E-10,
     4 .6644E+00, .1556E-02,-.4294E-07,-.6425E-08,-.2252E-10,
     4 .7851E+00, .1495E-02,-.2301E-06,-.4893E-08,-.9324E-10,
     4 .9304E+00, .1395E-02, .6525E-07, .9691E-08,-.2604E-09,
     4 .1099E+01, .1363E-02,-.1218E-05, .9560E-08,-.1188E-09,
     4 .1281E+01, .1346E-02,-.1086E-05, .3036E-08,-.1302E-09,
     4 .1477E+01, .1293E-02,-.2185E-05, .1361E-08, .4658E-10,
     4 .1675E+01, .1263E-02,-.1153E-05,-.6468E-08,-.9134E-10,
     4 .1877E+01, .1182E-02,-.1615E-05,-.4820E-09, .7266E-10 /
      DATA ((( CS2O2(K,J,I), K = 1, 5), J = 1, 18), I = 1, 2) /
     1 .3214E-07, .1096E-09,-.3452E-12,-.9474E-15, .1212E-16,
     1 .4547E-07, .1549E-09,-.4458E-12,-.1710E-14, .1636E-16,
     1 .6451E-07, .2141E-09,-.6686E-12,-.1753E-14, .4046E-16,
     1 .9148E-07, .2994E-09,-.7666E-12,-.2964E-14, .2290E-16,
     1 .1291E-06, .4214E-09,-.4340E-12,-.5534E-14,-.6547E-16,
     1 .1826E-06, .5775E-09,-.4514E-12,-.8455E-14,-.8933E-16,
     1 .2577E-06, .7778E-09, .1163E-12,-.1062E-13,-.2128E-15,
     1 .3647E-06, .9986E-09,-.5075E-12,-.4184E-14,-.1315E-15,
     1 .5161E-06, .1342E-08,-.1859E-11,-.5709E-14, .1397E-15,
     1 .7339E-06, .1777E-08,-.2617E-11,-.2611E-14, .2048E-15,
     1 .1044E-05, .2394E-08,-.2534E-11,-.2946E-14,-.1257E-16,
     1 .1486E-05, .3198E-08, .9352E-13, .9540E-14,-.7692E-15,
     1 .2129E-05, .4310E-08, .2237E-13, .3356E-13,-.1033E-14,
     1 .3067E-05, .5636E-08,-.1386E-10, .1004E-12, .7494E-15,
     1 .4422E-05, .7600E-08,-.2434E-10, .1677E-12, .6683E-15,
     1 .6379E-05, .1039E-07,-.4778E-10, .2399E-12, .1082E-14,
     1 .9118E-05, .1456E-07,-.6065E-10, .2221E-12,-.2072E-14,
     1 .1297E-04, .2085E-07,-.1135E-09, .2214E-12, .1963E-14,
     2 .1314E-08, .9760E-12, .8168E-13, .6512E-15,-.2238E-16,
     2 .1686E-08,-.3646E-11, .2078E-12, .1897E-14,-.2439E-16,
     2 .2586E-08, .5144E-11,-.3531E-12, .2446E-16, .7444E-16,
     2 .2835E-08, .1092E-10, .6488E-12,-.9216E-15,-.7964E-16,
     2 .7028E-08, .9956E-11,-.1326E-11,-.2252E-14, .1800E-15,
     2 .9576E-08, .5140E-10,-.1396E-11,-.7284E-14, .2245E-15,
     2 .1471E-07, .3408E-11,-.8488E-12, .3609E-14, .1276E-15,
     2 .2307E-07, .1126E-09,-.3998E-11,-.1265E-13, .7476E-15,
     2 .3514E-07, .1715E-09,-.4584E-11,-.1617E-13, .7056E-15,
     2 .4404E-07, .3456E-09,-.1496E-11,-.3771E-13, .1827E-15,
     2 .5772E-07, .6428E-09, .8844E-11,-.8340E-13,-.1589E-14,
     2 .7356E-07, .5488E-09, .7944E-11,-.4340E-13,-.2486E-15,
     2 .1247E-06, .7840E-09,-.6232E-11,-.6144E-13, .2092E-14,
     2 .1674E-06, .1256E-08, .2241E-10,-.1350E-12,-.3142E-14,
     2 .2637E-06, .8356E-09,-.4084E-10,-.1284E-12, .7640E-14,
     2 .2975E-06, .9760E-09,-.6176E-10,-.1859E-12, .1374E-13,
     2 .3092E-06, .3456E-08, .6784E-10,-.3132E-12,-.9412E-14,
     2 .5828E-06, .5524E-08,-.2897E-10,-.1702E-12, .3144E-14 /
C
      COMMON /BANDS2GH/ GWGH(4), CS2H2OGH(5,28), CS2O2GH(5,28,3)
      DATA GWGH / 4.57116, 3.93806, .30706, .08728 /
      DATA (( CS2H2OGH(K,J), K = 1, 5), J = 1, 28)           /
     5 .3219E+02,-.7035E-02, .5768E-04,-.4975E-06, .1310E-08,
     5 .3219E+02,-.7034E-02, .5773E-04,-.4981E-06, .1309E-08,
     5 .3219E+02,-.7032E-02, .5766E-04,-.4985E-06, .1322E-08,
     5 .3219E+02,-.7032E-02, .5774E-04,-.4980E-06, .1303E-08,
     5 .3219E+02,-.7031E-02, .5766E-04,-.4979E-06, .1316E-08,
     5 .3219E+02,-.7027E-02, .5775E-04,-.4980E-06, .1298E-08,
     5 .3219E+02,-.7021E-02, .5762E-04,-.4975E-06, .1306E-08,
     5 .3219E+02,-.7007E-02, .5760E-04,-.4964E-06, .1288E-08,
     5 .3219E+02,-.6979E-02, .5745E-04,-.4947E-06, .1269E-08,
     5 .3217E+02,-.6926E-02, .5686E-04,-.4910E-06, .1287E-08,
     5 .3212E+02,-.6725E-02, .5586E-04,-.4840E-06, .1227E-08,
     5 .3208E+02,-.6588E-02, .5516E-04,-.4807E-06, .1221E-08,
     5 .3203E+02,-.6416E-02, .5448E-04,-.4763E-06, .1202E-08,
     5 .3195E+02,-.6215E-02, .5359E-04,-.4679E-06, .1130E-08,
     5 .3184E+02,-.5987E-02, .5169E-04,-.4484E-06, .1101E-08,
     5 .3169E+02,-.5687E-02, .4792E-04,-.4268E-06, .1274E-08,
     5 .3148E+02,-.5341E-02, .4536E-04,-.3987E-06, .1130E-08,
     5 .3120E+02,-.5072E-02, .4504E-04,-.3405E-06, .4385E-09,
     5 .3082E+02,-.4734E-02, .4023E-04,-.2752E-06, .1470E-09,
     5 .3033E+02,-.4165E-02, .3385E-04,-.2314E-06, .1184E-10,
     5 .2969E+02,-.3520E-02, .2718E-04,-.1814E-06,-.1465E-09,
     5 .2890E+02,-.2987E-02, .2116E-04,-.1337E-06,-.2584E-10,
     5 .2789E+02,-.2448E-02, .1409E-04,-.1685E-06, .8043E-09,
     5 .2669E+02,-.2465E-02, .1763E-04,-.1577E-06, .1860E-09,
     5 .2530E+02,-.2808E-02, .1733E-04,-.1291E-06, .1692E-09,
     5 .2372E+02,-.3219E-02, .2400E-04,-.1230E-06,-.6400E-09,
     5 .2201E+02,-.3839E-02, .2112E-04,-.9221E-07,-.9607E-11,
     5 .2020E+02,-.4393E-02, .2570E-04,-.1237E-06,-.6172E-09 /
      DATA ((( CS2O2GH(K,J,I), K = 1, 5), J = 1, 28), I = 1, 3) /
     6 .2534E-03, .2682E-05, .6224E-08, .9713E-11, .3900E-12,
     6 .2534E-03, .2682E-05, .6225E-08, .9713E-11, .3899E-12,
     6 .2534E-03, .2682E-05, .6227E-08, .9713E-11, .3898E-12,
     6 .2536E-03, .2681E-05, .6232E-08, .9734E-11, .3892E-12,
     6 .2539E-03, .2681E-05, .6245E-08, .9778E-11, .3876E-12,
     6 .2546E-03, .2678E-05, .6251E-08, .9778E-11, .3886E-12,
     6 .2560E-03, .2675E-05, .6262E-08, .9698E-11, .3902E-12,
     6 .2588E-03, .2668E-05, .6284E-08, .9540E-11, .3935E-12,
     6 .2646E-03, .2652E-05, .6311E-08, .9317E-11, .4035E-12,
     6 .2763E-03, .2620E-05, .6450E-08, .9634E-11, .4072E-12,
     6 .3205E-03, .2486E-05, .7034E-08, .1581E-10, .3879E-12,
     6 .3532E-03, .2375E-05, .7308E-08, .2418E-10, .3959E-12,
     6 .4002E-03, .2213E-05, .7531E-08, .4023E-10, .4261E-12,
     6 .4674E-03, .2023E-05, .8093E-08, .6237E-10, .3657E-12,
     6 .5594E-03, .1903E-05, .1140E-07, .7081E-10,-.1552E-12,
     6 .6864E-03, .1912E-05, .1392E-07, .6452E-10,-.3394E-12,
     6 .8669E-03, .2013E-05, .1536E-07, .5262E-10,-.3039E-12,
     6 .1127E-02, .2129E-05, .1106E-07, .4940E-10, .6748E-12,
     6 .1499E-02, .2219E-05, .5355E-08, .4513E-10, .1711E-11,
     6 .2025E-02, .2429E-05, .1135E-07, .1094E-10, .3164E-12,
     6 .2778E-02, .2537E-05, .1981E-07, .1433E-11,-.1488E-11,
     6 .3846E-02, .1850E-05, .2369E-07, .1651E-09,-.2853E-11,
     6 .5333E-02, .1444E-05, .1251E-07, .2508E-09,-.1932E-11,
     6 .7351E-02, .8431E-06,-.7258E-08, .2651E-09, .7819E-12,
     6 .1003E-01, .9900E-06,-.3748E-07, .3186E-10, .6088E-11,
     6 .1343E-01, .1380E-05, .8575E-08,-.3203E-09, .2250E-12,
     6 .1773E-01,-.1315E-05, .1153E-07,-.5970E-10, .6931E-12,
     6 .2259E-01,-.4343E-05, .5704E-07, .2123E-09,-.8021E-11,
     7 .1446E+00, .1154E-02, .4288E-06,-.7609E-08,-.4095E-10,
     7 .1446E+00, .1154E-02, .4288E-06,-.7611E-08,-.4094E-10,
     7 .1446E+00, .1154E-02, .4285E-06,-.7605E-08,-.4087E-10,
     7 .1447E+00, .1154E-02, .4282E-06,-.7603E-08,-.4080E-10,
     7 .1447E+00, .1154E-02, .4284E-06,-.7602E-08,-.4082E-10,
     7 .1447E+00, .1154E-02, .4288E-06,-.7606E-08,-.4093E-10,
     7 .1447E+00, .1154E-02, .4293E-06,-.7600E-08,-.4099E-10,
     7 .1448E+00, .1154E-02, .4310E-06,-.7597E-08,-.4125E-10,
     7 .1449E+00, .1153E-02, .4327E-06,-.7587E-08,-.4145E-10,
     7 .1451E+00, .1152E-02, .4374E-06,-.7583E-08,-.4212E-10,
     7 .1459E+00, .1149E-02, .4512E-06,-.7645E-08,-.4352E-10,
     7 .1464E+00, .1148E-02, .4639E-06,-.7753E-08,-.4495E-10,
     7 .1473E+00, .1144E-02, .4855E-06,-.7968E-08,-.4735E-10,
     7 .1487E+00, .1138E-02, .5143E-06,-.8306E-08,-.4982E-10,
     7 .1507E+00, .1129E-02, .5171E-06,-.8427E-08,-.4663E-10,
     7 .1539E+00, .1109E-02, .5067E-06,-.7951E-08,-.4063E-10,
     7 .1588E+00, .1078E-02, .5143E-06,-.7148E-08,-.3427E-10,
     7 .1658E+00, .1035E-02, .6376E-06,-.6602E-08,-.4582E-10,
     7 .1760E+00, .9738E-03, .8113E-06,-.5868E-08,-.5717E-10,
     7 .1906E+00, .8877E-03, .8450E-06,-.5023E-08,-.3616E-10,
     7 .2109E+00, .7691E-03, .9423E-06,-.1783E-08,-.5308E-10,
     7 .2376E+00, .6315E-03, .1032E-05, .3848E-10,-.8310E-10,
     7 .2707E+00, .5010E-03, .1865E-06,-.5449E-08, .4747E-10,
     7 .3057E+00, .3436E-03, .5605E-06,-.3812E-08,-.3625E-10,
     7 .3370E+00, .1751E-03, .1381E-05, .7458E-08,-.2713E-09,
     7 .3559E+00, .1121E-03,-.9225E-06, .2911E-08, .9792E-10,
     7 .3505E+00, .1067E-03, .4705E-06,-.5564E-08,-.1082E-09,
     7 .3253E+00, .8664E-04,-.8188E-06,-.3145E-08, .1218E-09,
     8 .3026E+01,-.5643E-02, .4530E-05, .1589E-08, .9835E-10,
     8 .3026E+01,-.5643E-02, .4537E-05, .1583E-08, .9716E-10,
     8 .3026E+01,-.5642E-02, .4538E-05, .1545E-08, .9709E-10,
     8 .3026E+01,-.5642E-02, .4539E-05, .1520E-08, .9684E-10,
     8 .3026E+01,-.5642E-02, .4537E-05, .1589E-08, .9710E-10,
     8 .3025E+01,-.5642E-02, .4520E-05, .1585E-08, .9949E-10,
     8 .3025E+01,-.5641E-02, .4533E-05, .1574E-08, .9609E-10,
     8 .3025E+01,-.5639E-02, .4505E-05, .1615E-08, .1000E-09,
     8 .3024E+01,-.5635E-02, .4492E-05, .1654E-08, .9907E-10,
     8 .3022E+01,-.5626E-02, .4456E-05, .1708E-08, .9894E-10,
     8 .3015E+01,-.5595E-02, .4286E-05, .1900E-08, .1066E-09,
     8 .3010E+01,-.5573E-02, .4169E-05, .2084E-08, .1109E-09,
     8 .3003E+01,-.5539E-02, .3995E-05, .2221E-08, .1165E-09,
     8 .2991E+01,-.5490E-02, .3741E-05, .2785E-08, .1268E-09,
     8 .2975E+01,-.5418E-02, .3364E-05, .3507E-08, .1439E-09,
     8 .2951E+01,-.5306E-02, .2991E-05, .3147E-08, .1315E-09,
     8 .2916E+01,-.5141E-02, .2528E-05, .2387E-08, .1005E-09,
     8 .2865E+01,-.4907E-02, .1848E-05, .2590E-08, .7710E-10,
     8 .2792E+01,-.4578E-02, .9081E-06, .2893E-08, .4496E-10,
     8 .2689E+01,-.4127E-02,-.4217E-06, .5231E-08, .4274E-10,
     8 .2546E+01,-.3517E-02,-.2308E-05,-.2749E-08, .2130E-09,
     8 .2355E+01,-.2767E-02,-.3932E-05,-.1626E-07, .4084E-09,
     8 .2111E+01,-.2072E-02,-.8738E-06, .2197E-08,-.1550E-09,
     8 .1827E+01,-.1303E-02,-.1590E-05,-.4515E-08,-.1066E-10,
     8 .1518E+01,-.6041E-03,-.2937E-05,-.3415E-07, .5588E-09,
     8 .1211E+01,-.3898E-03, .3294E-05, .4595E-08,-.5044E-09,
     8 .9397E+00,-.2416E-03,-.2158E-05, .2068E-07, .2944E-09,
     8 .6996E+00,-.4103E-04, .3979E-06,-.6163E-08,-.5776E-10 /
C
C----------------------------------------------------------------------C
C     SPECTRAL REGION: 4200-8400 CM^-1. GAS: H2O, CO2                  C
C----------------------------------------------------------------------C
C
      COMMON /BANDS3/ GWS3(6), CS3H2O(5,18,6), CS3CO2(5,18,6),
     1                         CS3CH4(5,18,2)
      DATA GWS3 / 126.08439, 41.86427, 28.31155, 22.53174, 15.24951,
     1              4.65647 /
      DATA ((( CS3H2O(K,J,I), K = 1, 5), J = 1, 18), I = 1, 3) /
     1 .5243E-05, .3937E-07, .6392E-10,-.7654E-13,-.1331E-14,
     1 .7756E-05, .5803E-07, .7734E-10,-.6975E-13, .8095E-15,
     1 .1147E-04, .8604E-07, .1175E-09,-.1330E-12,-.6863E-16,
     1 .1691E-04, .1277E-06, .1985E-09,-.2722E-12,-.4618E-14,
     1 .2499E-04, .1892E-06, .2877E-09,-.4275E-12,-.5520E-14,
     1 .3690E-04, .2811E-06, .4177E-09,-.9578E-12,-.7466E-14,
     1 .5447E-04, .4168E-06, .5869E-09,-.1550E-11,-.2811E-14,
     1 .8039E-04, .6156E-06, .8308E-09,-.1998E-11, .2946E-14,
     1 .1186E-03, .9149E-06, .1332E-08,-.4103E-11,-.1431E-13,
     1 .1748E-03, .1355E-05, .2219E-08,-.7211E-11,-.6379E-13,
     1 .2566E-03, .2002E-05, .3330E-08,-.1153E-10,-.8348E-13,
     1 .3752E-03, .2913E-05, .4650E-08,-.1277E-10,-.5410E-13,
     1 .5445E-03, .4232E-05, .7270E-08,-.2068E-10,-.1893E-12,
     1 .7817E-03, .6098E-05, .1068E-07,-.2949E-10,-.2664E-12,
     1 .1107E-02, .8662E-05, .1559E-07,-.3382E-10,-.2637E-12,
     1 .1548E-02, .1207E-04, .2050E-07,-.4091E-10,-.8716E-13,
     1 .2119E-02, .1640E-04, .2884E-07,-.3157E-10,-.2281E-12,
     1 .2836E-02, .2195E-04, .3520E-07,-.6713E-10, .2683E-12,
     2 .3170E-03, .2578E-05, .2844E-08,-.2714E-10, .2976E-13,
     2 .4557E-03, .3688E-05, .4038E-08,-.3660E-10, .9022E-13,
     2 .6542E-03, .5284E-05, .6977E-08,-.5251E-10,-.4629E-13,
     2 .9362E-03, .7518E-05, .1126E-07,-.6840E-10,-.1460E-12,
     2 .1339E-02, .1069E-04, .1826E-07,-.9363E-10,-.3904E-12,
     2 .1912E-02, .1510E-04, .3010E-07,-.1283E-09,-.1136E-11,
     2 .2729E-02, .2113E-04, .4543E-07,-.1308E-09,-.1698E-11,
     2 .3902E-02, .2947E-04, .6188E-07,-.1093E-09,-.1660E-11,
     2 .5571E-02, .4130E-04, .9479E-07,-.1062E-09,-.3316E-11,
     2 .7949E-02, .5766E-04, .1278E-06,-.1065E-09,-.3532E-11,
     2 .1131E-01, .8021E-04, .1722E-06,-.3606E-10,-.3278E-11,
     2 .1603E-01, .1108E-03, .2378E-06, .5192E-10,-.5949E-11,
     2 .2257E-01, .1525E-03, .3034E-06, .1086E-09,-.6218E-11,
     2 .3145E-01, .2093E-03, .3541E-06, .1289E-09, .1560E-11,
     2 .4310E-01, .2838E-03, .5302E-06, .2381E-09,-.3963E-11,
     2 .5820E-01, .3758E-03, .6697E-06, .4019E-09,-.4028E-11,
     2 .7692E-01, .4853E-03, .8518E-06, .1027E-08,-.1271E-10,
     2 .9924E-01, .6198E-03, .8679E-06,-.2819E-09, .1600E-10,
     3 .4048E-02, .2608E-04, .1568E-07,-.2063E-09,-.1770E-11,
     3 .5686E-02, .3555E-04, .1363E-07,-.2302E-09,-.1056E-11,
     3 .7993E-02, .4877E-04, .2005E-07,-.3327E-09,-.2128E-11,
     3 .1125E-01, .6660E-04, .1436E-07,-.3712E-09,-.6695E-12,
     3 .1588E-01, .9130E-04,-.7789E-08,-.4876E-09, .3794E-11,
     3 .2240E-01, .1251E-03,-.2941E-08,-.6350E-09, .3433E-11,
     3 .3168E-01, .1712E-03,-.1668E-07,-.6789E-09, .6265E-11,
     3 .4495E-01, .2362E-03,-.3829E-07,-.1035E-08, .7544E-11,
     3 .6387E-01, .3275E-03,-.1625E-07,-.1633E-08,-.1103E-11,
     3 .9082E-01, .4517E-03,-.4195E-07,-.1946E-08,-.5407E-11,
     3 .1289E+00, .6238E-03,-.2491E-06,-.2302E-08, .2520E-10,
     3 .1822E+00, .8568E-03,-.1757E-06,-.1819E-08,-.1196E-12,
     3 .2566E+00, .1170E-02,-.7726E-07,-.2840E-10,-.6146E-10,
     3 .3592E+00, .1603E-02,-.7184E-06, .7504E-09,-.1300E-10,
     3 .4970E+00, .2192E-02,-.7953E-06, .1548E-08,-.8341E-10,
     3 .6774E+00, .2957E-02,-.1090E-05, .3594E-08,-.1495E-09,
     3 .9085E+00, .3914E-02,-.3384E-05, .1239E-07, .2922E-10,
     3 .1195E+01, .5226E-02,-.7062E-05, .2829E-08, .3743E-09 /
      DATA ((( CS3H2O(K,J,I), K = 1, 5), J = 1, 18), I = 4, 6) /
     4 .2654E-01, .1791E-03, .5150E-06, .3716E-09,-.6728E-11,
     4 .3574E-01, .2209E-03, .5579E-06, .1210E-09,-.1067E-10,
     4 .4848E-01, .2743E-03, .6127E-06, .9519E-10,-.1489E-10,
     4 .6618E-01, .3410E-03, .6880E-06, .4258E-09,-.2173E-10,
     4 .9106E-01, .4297E-03, .6362E-06, .3044E-09,-.8938E-11,
     4 .1258E+00, .5429E-03, .5158E-06, .6320E-09, .2041E-10,
     4 .1748E+00, .6866E-03, .4770E-06, .1558E-08, .2773E-10,
     4 .2439E+00, .8774E-03, .5859E-06, .2075E-08, .2138E-11,
     4 .3417E+00, .1140E-02, .1162E-05, .5718E-09,-.1289E-09,
     4 .4813E+00, .1474E-02, .4319E-06, .1747E-08,-.3889E-10,
     4 .6775E+00, .1915E-02, .1935E-06, .3808E-08,-.5681E-10,
     4 .9536E+00, .2473E-02,-.8389E-06, .9849E-08, .8561E-11,
     4 .1337E+01, .3202E-02,-.2386E-05, .1403E-07, .1709E-10,
     4 .1858E+01, .4109E-02,-.4965E-05, .2996E-07, .1398E-09,
     4 .2557E+01, .5331E-02,-.7711E-05, .4100E-07, .9088E-10,
     4 .3477E+01, .6888E-02,-.1284E-04, .3892E-07, .3316E-09,
     4 .4669E+01, .8815E-02,-.1822E-04, .3060E-07, .2788E-09,
     4 .6170E+01, .1123E-01,-.2198E-04, .3073E-07,-.5434E-09,
     5 .4596E+00, .4624E-02, .2257E-04, .2845E-07,-.3622E-09,
     5 .5488E+00, .4860E-02, .2003E-04, .3360E-07,-.3368E-10,
     5 .6664E+00, .5206E-02, .1995E-04, .3132E-07,-.1326E-09,
     5 .8255E+00, .5611E-02, .1856E-04, .2927E-07,-.4059E-10,
     5 .1040E+01, .6175E-02, .1834E-04, .1304E-07,-.1583E-09,
     5 .1330E+01, .6802E-02, .1795E-04, .8607E-08,-.2169E-09,
     5 .1728E+01, .7468E-02, .1507E-04, .2537E-07, .2628E-09,
     5 .2270E+01, .8344E-02, .1880E-04, .2439E-07,-.5943E-09,
     5 .3020E+01, .9317E-02, .2190E-04, .2841E-07,-.1494E-08,
     5 .4049E+01, .1023E-01, .1697E-04, .7381E-07,-.1003E-08,
     5 .5451E+01, .1143E-01, .8803E-05, .9451E-07,-.3631E-09,
     5 .7344E+01, .1282E-01,-.8213E-06, .1144E-06,-.9634E-10,
     5 .9836E+01, .1453E-01,-.1769E-04, .1035E-06, .1590E-08,
     5 .1305E+02, .1635E-01,-.3205E-04, .1156E-06, .2456E-08,
     5 .1707E+02, .1775E-01,-.1976E-04, .1495E-06,-.1671E-08,
     5 .2195E+02, .1731E-01,-.2906E-04, .3484E-06,-.2427E-08,
     5 .2761E+02, .1700E-01,-.6037E-04, .3231E-06, .5769E-09,
     5 .3391E+02, .1467E-01,-.3589E-04, .4394E-06,-.6147E-08,
     6 .2445E+02, .2122E+00, .4472E-03,-.1094E-05,-.1414E-07,
     6 .2537E+02, .2084E+00, .4022E-03,-.9312E-06,-.6623E-08,
     6 .2666E+02, .2046E+00, .4210E-03,-.8415E-06,-.1126E-07,
     6 .2854E+02, .1985E+00, .3881E-03,-.3549E-06,-.1044E-07,
     6 .3121E+02, .1954E+00, .2828E-03,-.8356E-06, .9093E-09,
     6 .3480E+02, .1925E+00, .2948E-03,-.1152E-05,-.8685E-08,
     6 .3976E+02, .1839E+00, .1594E-03,-.4909E-06, .7086E-08,
     6 .4625E+02, .1761E+00, .1330E-03,-.2552E-06, .2167E-08,
     6 .5486E+02, .1651E+00, .2576E-03, .1713E-06,-.3532E-07,
     6 .6619E+02, .1522E+00, .3718E-04, .4508E-06,-.9015E-08,
     6 .8028E+02, .1386E+00,-.3606E-04, .5855E-06,-.8296E-08,
     6 .9724E+02, .1263E+00,-.7633E-04, .2208E-06,-.1605E-07,
     6 .1169E+03, .1075E+00,-.1033E-03, .1081E-05,-.2280E-07,
     6 .1385E+03, .8503E-01,-.3415E-03, .2483E-05, .1085E-07,
     6 .1611E+03, .7763E-01,-.2098E-03, .4268E-06,-.2197E-07,
     6 .1840E+03, .6100E-01,-.4282E-03, .4791E-06, .4600E-08,
     6 .2041E+03, .3765E-01,-.2582E-03, .9662E-06,-.2339E-07,
     6 .2196E+03, .9987E-02,-.5177E-03, .1695E-05, .1555E-07 /
      DATA ((( CS3CO2(K,J,I), K = 1, 5), J = 1, 18), I = 1, 3) /
     1 .1589E-03, .1640E-05, .7753E-08, .1582E-10, .1076E-12,
     1 .1969E-03, .1903E-05, .1209E-07, .7832E-11,-.5050E-12,
     1 .2501E-03, .2122E-05, .1181E-07, .1903E-10,-.1864E-12,
     1 .3257E-03, .2451E-05, .1067E-07, .2329E-10, .2129E-12,
     1 .4251E-03, .2872E-05, .1272E-07, .2731E-10,-.1337E-13,
     1 .5639E-03, .3319E-05, .1416E-07, .5936E-10, .4578E-13,
     1 .7565E-03, .4035E-05, .1631E-07, .5563E-10,-.3285E-12,
     1 .1015E-02, .5052E-05, .2185E-07, .2752E-10,-.1078E-11,
     1 .1378E-02, .6082E-05, .1394E-07, .6454E-10, .3970E-12,
     1 .1878E-02, .7785E-05,-.2238E-08, .7919E-10, .3187E-11,
     1 .2537E-02, .9750E-05,-.3206E-08, .1284E-09, .3600E-11,
     1 .3385E-02, .1406E-04, .1926E-07,-.1047E-09, .5494E-12,
     1 .4507E-02, .1806E-04, .1531E-07,-.1278E-09, .3997E-11,
     1 .5905E-02, .2093E-04, .3526E-07, .3254E-10, .2282E-11,
     1 .7671E-02, .2429E-04, .2548E-07, .1388E-09, .2323E-11,
     1 .9702E-02, .3106E-04, .5056E-07, .7348E-10, .1621E-11,
     1 .1215E-01, .3862E-04, .5259E-07,-.5148E-09, .9232E-11,
     1 .1469E-01, .3781E-04, .2454E-06, .2456E-09,-.1958E-10,
     2 .2105E-03, .2288E-05, .1443E-07, .7527E-11,-.5328E-12,
     2 .2738E-03, .2583E-05, .6724E-08, .1626E-10, .8559E-12,
     2 .3541E-03, .3170E-05, .9352E-08, .3302E-11, .4627E-12,
     2 .4545E-03, .3917E-05, .2305E-07,-.3317E-10,-.2136E-11,
     2 .6052E-03, .4615E-05, .2623E-07, .7403E-11,-.2319E-11,
     2 .8088E-03, .5913E-05, .2950E-07,-.4464E-10,-.2728E-11,
     2 .1097E-02, .7378E-05, .1852E-07,-.3892E-10, .1180E-12,
     2 .1485E-02, .9226E-05, .1140E-07,-.1022E-10, .2273E-11,
     2 .2006E-02, .1190E-04, .9485E-08,-.1098E-09, .4375E-11,
     2 .2678E-02, .1500E-04, .6989E-07,-.1169E-09,-.3755E-11,
     2 .3659E-02, .2072E-04, .1151E-06,-.5751E-09,-.1163E-10,
     2 .4912E-02, .2360E-04, .7839E-07,-.2371E-09, .1976E-12,
     2 .6481E-02, .2805E-04, .9920E-07,-.2629E-09,-.3809E-11,
     2 .8480E-02, .3527E-04, .1910E-06,-.4034E-09,-.2034E-10,
     2 .1086E-01, .4484E-04, .2655E-06,-.6404E-09,-.3327E-10,
     2 .1411E-01, .5367E-04, .4105E-07,-.1058E-08,-.1568E-12,
     2 .1767E-01, .5685E-04, .2495E-07, .1098E-09,-.6352E-11,
     2 .2177E-01, .5640E-04,-.1901E-06, .2717E-08,-.4546E-11,
     3 .1595E-03, .1355E-05,-.1727E-09, .1354E-10, .8839E-12,
     3 .2014E-03, .1621E-05, .2632E-08, .9287E-11, .5670E-12,
     3 .2616E-03, .1946E-05, .4025E-08,-.1691E-11, .3135E-12,
     3 .3473E-03, .2449E-05, .8980E-09,-.1605E-10, .1302E-11,
     3 .4594E-03, .2941E-05, .8607E-08, .2806E-11, .4789E-12,
     3 .6270E-03, .3602E-05, .5125E-08,-.1165E-10, .1216E-11,
     3 .8471E-03, .4769E-05, .1624E-07,-.9563E-10,-.7373E-12,
     3 .1191E-02, .5463E-05, .4851E-08,-.2268E-10, .1258E-11,
     3 .1662E-02, .6062E-05, .5503E-08, .1767E-09, .2756E-12,
     3 .2283E-02, .7847E-05, .2112E-07, .8182E-10,-.2619E-11,
     3 .3137E-02, .1043E-04, .1308E-07, .7045E-10,-.2083E-12,
     3 .4390E-02, .1149E-04, .3458E-07, .2232E-09,-.9584E-11,
     3 .5971E-02, .1299E-04, .3109E-07, .6560E-09,-.1032E-10,
     3 .8066E-02, .2216E-04,-.5196E-07, .1403E-09,-.3315E-12,
     3 .1084E-01, .3230E-04,-.2381E-06, .2725E-09, .3041E-10,
     3 .1452E-01, .4945E-04,-.1681E-06,-.3496E-09, .6259E-11,
     3 .1870E-01, .7177E-04, .1857E-06,-.9727E-10,-.6312E-10,
     3 .2420E-01, .1144E-03,-.3591E-06,-.3906E-08, .4705E-10 /
      DATA ((( CS3CO2(K,J,I), K = 1, 5), J = 1, 18), I = 4, 6) /
     4 .1149E-03, .1470E-05, .1180E-07,-.4522E-11,-.9234E-12,
     4 .1449E-03, .1609E-05, .9413E-08, .1111E-10,-.4038E-12,
     4 .1875E-03, .1964E-05, .1098E-07,-.1928E-10,-.8581E-12,
     4 .2418E-03, .2192E-05, .7245E-08, .5254E-10, .3422E-12,
     4 .3185E-03, .2861E-05,-.2634E-08,-.8221E-11, .2076E-11,
     4 .4013E-03, .3605E-05, .1092E-07,-.3672E-10,-.9761E-13,
     4 .5047E-03, .3985E-05, .2437E-07, .3069E-10,-.1265E-11,
     4 .6761E-03, .4702E-05, .4147E-07, .6096E-10,-.4863E-11,
     4 .9043E-03, .6012E-05, .6513E-07, .4235E-10,-.7505E-11,
     4 .1261E-02, .7848E-05, .6756E-07,-.6823E-10,-.8912E-11,
     4 .1735E-02, .9775E-05, .2034E-07, .1155E-09, .1017E-11,
     4 .2333E-02, .1269E-04, .2980E-07, .1681E-09, .1106E-11,
     4 .3137E-02, .1708E-04, .8946E-07, .1505E-09,-.1089E-10,
     4 .4193E-02, .2075E-04,-.1380E-06, .1071E-09, .3090E-10,
     4 .5201E-02, .2904E-04, .1113E-07, .8138E-10, .1346E-10,
     4 .6600E-02, .3531E-04, .9985E-07, .4916E-09, .9948E-12,
     4 .8448E-02, .4800E-04,-.9133E-07, .4154E-09, .4552E-10,
     4 .1013E-01, .6543E-04,-.8460E-07,-.1921E-08, .6441E-10,
     5 .1043E-03, .1342E-05, .9444E-08, .4421E-10, .1312E-12,
     5 .1341E-03, .1523E-05, .7051E-08, .3751E-10, .3297E-12,
     5 .1711E-03, .1823E-05, .3842E-08, .1601E-10, .7611E-12,
     5 .2039E-03, .2098E-05, .1038E-07, .1185E-10,-.6357E-13,
     5 .2586E-03, .2693E-05, .1116E-07,-.2386E-10, .8329E-13,
     5 .3329E-03, .2991E-05, .5140E-08, .9912E-11, .1518E-11,
     5 .4322E-03, .3438E-05, .1227E-07, .3771E-10,-.1171E-12,
     5 .5531E-03, .4045E-05, .9454E-09, .1085E-09, .3629E-11,
     5 .7282E-03, .7010E-05, .4254E-07,-.2314E-09,-.3931E-11,
     5 .1028E-02, .9036E-05,-.1828E-08,-.2506E-09, .4968E-11,
     5 .1299E-02, .1069E-04, .2735E-07,-.1173E-09, .3509E-11,
     5 .1669E-02, .1086E-04, .8583E-07, .4590E-09,-.7201E-12,
     5 .2288E-02, .1309E-04, .2544E-07, .2997E-09, .6622E-11,
     5 .2989E-02, .1832E-04, .1007E-06, .4619E-09,-.4546E-11,
     5 .3685E-02, .2291E-04, .1443E-06, .7279E-09,-.5163E-11,
     5 .4473E-02, .2992E-04, .1805E-07, .5183E-09, .2606E-10,
     5 .5317E-02, .3785E-04, .1469E-06,-.8436E-09,-.5855E-11,
     5 .5640E-02, .4414E-04, .4487E-06,-.1092E-08,-.3665E-10,
     6 .1271E-03, .2029E-05, .2419E-07,-.9357E-10,-.3749E-11,
     6 .1657E-03, .1906E-05, .3026E-08, .6135E-11, .5054E-12,
     6 .1977E-03, .2272E-05, .1052E-07, .6407E-10, .1615E-12,
     6 .2538E-03, .2609E-05,-.1819E-08, .2386E-10, .1774E-11,
     6 .3033E-03, .3206E-05, .1380E-07,-.3174E-10,-.7921E-12,
     6 .3598E-03, .4512E-05, .3275E-07,-.1842E-09,-.3494E-11,
     6 .4649E-03, .4377E-05,-.1060E-07, .9835E-10, .6697E-11,
     6 .5968E-03, .4370E-05,-.3778E-07, .2550E-09, .1205E-10,
     6 .7410E-03, .5270E-05, .1088E-07, .2263E-09, .2317E-11,
     6 .9389E-03, .5604E-05, .3888E-07, .2878E-09,-.2994E-11,
     6 .1198E-02, .9329E-05, .2683E-07,-.3494E-09,-.9746E-12,
     6 .1567E-02, .1152E-04,-.7298E-07,-.5222E-09, .1277E-10,
     6 .1821E-02, .1332E-04,-.1054E-06, .1829E-09, .2918E-10,
     6 .2061E-02, .1207E-04,-.1083E-06, .6075E-09, .2954E-10,
     6 .1968E-02, .1176E-04, .1794E-07, .7608E-09, .6172E-11,
     6 .1665E-02, .2114E-04, .1158E-06,-.1138E-08,-.1305E-10,
     6 .1355E-02, .1766E-04,-.7097E-07,-.3349E-09, .3554E-10,
     6 .9925E-03, .1501E-04, .1896E-06,-.3547E-09,-.1132E-10 /
      DATA ((( CS3CH4(K,J,I), K = 1, 5), J = 1, 18), I = 1, 2) /
     1 .4760E+01,-.7433E-02,-.8651E-04, .3592E-06, .1382E-07,
     1 .4876E+01,-.8488E-02,-.1163E-03, .5120E-06, .1951E-07,
     1 .5001E+01,-.8892E-02,-.1293E-03, .5431E-06, .2234E-07,
     1 .5122E+01,-.9002E-02,-.1373E-03, .5421E-06, .2384E-07,
     1 .5215E+01,-.9283E-02,-.1133E-03, .5733E-06, .2014E-07,
     1 .5336E+01,-.9117E-02,-.1147E-03, .5152E-06, .2047E-07,
     1 .5443E+01,-.1019E-01,-.1222E-03, .6970E-06, .2244E-07,
     1 .5554E+01,-.9480E-02,-.1159E-03, .5438E-06, .2150E-07,
     1 .5655E+01,-.9119E-02,-.1051E-03, .3740E-06, .1915E-07,
     1 .5715E+01,-.9100E-02,-.4910E-04, .3550E-06, .9870E-08,
     1 .5829E+01,-.9211E-02,-.2367E-04, .3476E-06, .5370E-08,
     1 .5912E+01,-.9355E-02, .1798E-04, .3160E-06,-.1129E-08,
     1 .6047E+01,-.1031E-01, .2895E-04, .4024E-06,-.2268E-08,
     1 .6183E+01,-.1098E-01, .4426E-04, .4399E-06,-.4526E-08,
     1 .6319E+01,-.1119E-01, .4734E-04, .4298E-06,-.5364E-08,
     1 .6440E+01,-.1104E-01, .5316E-04, .2535E-06,-.4902E-08,
     1 .6572E+01,-.1089E-01, .1532E-04, .1447E-06, .1752E-08,
     1 .6705E+01,-.1126E-01, .3651E-05, .1196E-06, .3177E-08,
     2 .2399E+01, .1297E-01, .1701E-03,-.7983E-06,-.3179E-07,
     2 .2449E+01, .1640E-01, .2868E-03,-.1456E-05,-.5384E-07,
     2 .2535E+01, .1706E-01, .3039E-03,-.1278E-05,-.5427E-07,
     2 .2622E+01, .2030E-01, .4541E-03,-.1914E-05,-.8549E-07,
     2 .2879E+01, .2022E-01, .3127E-03,-.1835E-05,-.6273E-07,
     2 .2964E+01, .1708E-01, .4294E-03,-.1479E-05,-.8110E-07,
     2 .3152E+01, .1857E-01, .4116E-03,-.1968E-05,-.7731E-07,
     2 .3412E+01, .1507E-01, .3193E-03,-.1363E-05,-.6283E-07,
     2 .3568E+01, .1307E-01, .3454E-03,-.7835E-06,-.6551E-07,
     2 .3965E+01, .1210E-01, .1274E-03,-.6940E-06,-.3016E-07,
     2 .4150E+01, .1289E-01, .4190E-04,-.1047E-05,-.1544E-07,
     2 .4396E+01, .1023E-01,-.9786E-04,-.6204E-06, .7852E-08,
     2 .4434E+01, .1174E-01,-.1028E-03,-.8347E-06, .7415E-08,
     2 .4478E+01, .1244E-01,-.1047E-03,-.8041E-06, .6674E-08,
     2 .4665E+01, .1052E-01,-.1911E-03,-.6572E-06, .2439E-07,
     2 .4864E+01, .9705E-02,-.2284E-03,-.1546E-06, .2532E-07,
     2 .4958E+01, .7855E-02,-.7171E-04, .1349E-06,-.1875E-08,
     2 .5046E+01, .8271E-02,-.8230E-04, .7807E-07, .3644E-08 /
C
      COMMON /BANDS3GH/ GWS3GH(4),CS3H2OGH(5,28,2),CS3CO2GH(5,28,4) 
      DATA GWS3GH / 5.23113, 1.01931, 1.13661, .05827 /
      DATA ((( CS3H2OGH(K,J,I), K = 1, 5), J = 1, 28), I = 1, 2) /
     7 .4090E-01, .6614E-03, .3852E-05,-.5346E-07,-.5943E-09,
     7 .4155E-01, .6717E-03, .3917E-05,-.5408E-07,-.6069E-09,
     7 .4221E-01, .6816E-03, .3987E-05,-.5455E-07,-.6221E-09,
     7 .4228E-01, .6938E-03, .4479E-05,-.5088E-07,-.7647E-09,
     7 .4289E-01, .7043E-03, .4681E-05,-.5141E-07,-.8060E-09,
     7 .4379E-01, .7093E-03, .4503E-05,-.4840E-07,-.8074E-09,
     7 .4383E-01, .6954E-03, .5037E-05,-.4743E-07,-.8966E-09,
     7 .4450E-01, .7225E-03, .5526E-05,-.4968E-07,-.1017E-08,
     7 .4430E-01, .7241E-03, .7883E-05,-.4343E-07,-.1410E-08,
     7 .4405E-01, .7328E-03, .6082E-05,-.4128E-07,-.9625E-09,
     7 .4844E-01, .4219E-03,-.3186E-05, .1756E-07, .7569E-09,
     7 .4464E-01, .4431E-03, .3590E-05, .9537E-08,-.4138E-09,
     7 .4615E-01, .3056E-03, .4806E-06, .2421E-07, .5553E-10,
     7 .4596E-01, .4087E-03, .3229E-05,-.6091E-10,-.3594E-09,
     7 .5066E-01, .2991E-03, .4483E-05, .4745E-07,-.4432E-09,
     7 .6246E-01, .1451E-03, .8931E-06, .8218E-07, .1532E-09,
     7 .7743E-01, .2978E-03, .3218E-05, .2540E-07,-.8262E-09,
     7 .1141E+00, .1602E-03,-.1284E-04, .3373E-07, .1619E-08,
     7 .1351E+00, .1439E-03,-.2873E-04, .4696E-07, .4679E-08,
     7 .1248E+00, .6543E-03, .7262E-06,-.3586E-07,-.6705E-09,
     7 .1465E+00, .8661E-03,-.2202E-05,-.3825E-07,-.6210E-09,
     7 .1542E+00, .6856E-03,-.3397E-05, .6530E-07, .6882E-09,
     7 .1665E+00, .9407E-03, .1863E-04, .7827E-08,-.3823E-08,
     7 .1934E+00, .1811E-02, .2740E-04,-.1306E-06,-.5262E-08,
     7 .1931E+00, .1152E-02, .9134E-05, .4371E-07,-.1913E-09,
     7 .2062E+00, .1404E-02, .6415E-05, .7587E-08, .6793E-10,
     7 .2071E+00, .1299E-02, .1244E-04,-.6210E-08,-.1224E-08,
     7 .2052E+00, .1138E-02, .6010E-05, .3393E-07,-.1002E-09,
     8 .2692E+04,-.3887E+00,-.6215E-02, .1309E-04,-.2819E-07,
     8 .2692E+04,-.3886E+00,-.6211E-02, .1307E-04,-.2873E-07,
     8 .2692E+04,-.3887E+00,-.6211E-02, .1309E-04,-.2874E-07,
     8 .2692E+04,-.3887E+00,-.6214E-02, .1308E-04,-.2809E-07,
     8 .2692E+04,-.3887E+00,-.6211E-02, .1308E-04,-.2851E-07,
     8 .2692E+04,-.3887E+00,-.6216E-02, .1307E-04,-.2742E-07,
     8 .2691E+04,-.3889E+00,-.6203E-02, .1308E-04,-.2917E-07,
     8 .2691E+04,-.3889E+00,-.6188E-02, .1305E-04,-.3060E-07,
     8 .2690E+04,-.3889E+00,-.6176E-02, .1275E-04,-.2704E-07,
     8 .2688E+04,-.3890E+00,-.6135E-02, .1249E-04,-.2821E-07,
     8 .2682E+04,-.3883E+00,-.5995E-02, .1111E-04,-.2476E-07,
     8 .2678E+04,-.3894E+00,-.5865E-02, .9878E-05,-.2258E-07,
     8 .2672E+04,-.3947E+00,-.5683E-02, .8088E-05,-.1310E-07,
     8 .2663E+04,-.4096E+00,-.5369E-02, .5914E-05,-.3553E-08,
     8 .2650E+04,-.4442E+00,-.4699E-02, .4560E-05,-.2833E-07,
     8 .2632E+04,-.5039E+00,-.3942E-02, .3817E-05,-.4718E-07,
     8 .2608E+04,-.5748E+00,-.2875E-02, .1508E-06,-.1044E-06,
     8 .2575E+04,-.6519E+00,-.1817E-02,-.4788E-05,-.1221E-06,
     8 .2530E+04,-.7391E+00,-.1314E-02,-.9487E-05,-.2249E-07,
     8 .2466E+04,-.8336E+00,-.3735E-03,-.9715E-05,-.1955E-07,
     8 .2381E+04,-.8993E+00, .3273E-03,-.1061E-04,-.2583E-08,
     8 .2268E+04,-.9297E+00, .7265E-03,-.1119E-04, .2304E-07,
     8 .2126E+04,-.9183E+00, .4952E-03,-.1307E-04, .1552E-06,
     8 .1957E+04,-.9063E+00, .1787E-02,-.1557E-04,-.1646E-08,
     8 .1766E+04,-.9301E+00, .1712E-02,-.1116E-04, .5805E-07,
     8 .1553E+04,-.9041E+00, .2892E-02,-.1484E-04,-.3293E-07,
     8 .1334E+04,-.8676E+00, .2671E-02,-.1604E-04, .6876E-07,
     8 .1118E+04,-.7976E+00, .3832E-02,-.2184E-04,-.2388E-07 /
      DATA ((( CS3CO2GH(K,J,I), K = 1, 5), J = 1, 28), I = 1, 2) /
     7 .5998E-01, .8788E-03, .4237E-05, .5526E-08, .2949E-11,
     7 .5999E-01, .8788E-03, .4238E-05, .5521E-08, .2914E-11,
     7 .6001E-01, .8787E-03, .4240E-05, .5506E-08, .2808E-11,
     7 .6004E-01, .8784E-03, .4246E-05, .5475E-08, .2670E-11,
     7 .6012E-01, .8774E-03, .4254E-05, .5466E-08, .2877E-11,
     7 .6033E-01, .8754E-03, .4274E-05, .5455E-08, .1690E-11,
     7 .6083E-01, .8711E-03, .4292E-05, .5532E-08, .1135E-11,
     7 .6195E-01, .8635E-03, .4297E-05, .5841E-08,-.1452E-12,
     7 .6428E-01, .8536E-03, .4284E-05, .6397E-08,-.6745E-11,
     7 .6882E-01, .8496E-03, .4181E-05, .6368E-08,-.8115E-11,
     7 .8333E-01, .8678E-03, .3891E-05, .5286E-08, .2289E-11,
     7 .9233E-01, .8868E-03, .3900E-05, .4514E-08,-.9559E-11,
     7 .1044E+00, .9098E-03, .3980E-05, .4466E-08,-.2903E-10,
     7 .1205E+00, .9404E-03, .3859E-05, .4419E-08,-.4627E-13,
     7 .1424E+00, .9859E-03, .3757E-05, .3067E-08, .1465E-10,
     7 .1718E+00, .1050E-02, .3886E-05, .4607E-09,-.4396E-11,
     7 .2115E+00, .1132E-02, .3780E-05,-.2982E-08, .3386E-10,
     7 .2642E+00, .1188E-02, .3989E-05, .1327E-08, .2682E-10,
     7 .3336E+00, .1260E-02, .4054E-05, .5781E-08, .7471E-10,
     7 .4247E+00, .1381E-02, .4379E-05, .5213E-08, .6535E-10,
     7 .5454E+00, .1566E-02, .4875E-05, .4952E-08,-.4397E-10,
     7 .7041E+00, .1817E-02, .3137E-05,-.7538E-08, .2525E-09,
     7 .8996E+00, .2033E-02, .3527E-05,-.2078E-07, .2541E-09,
     7 .1133E+01, .2093E-02, .5075E-05,-.2984E-07, .2291E-09,
     7 .1410E+01, .2002E-02, .7502E-05,-.5423E-07, .2280E-09,
     7 .1745E+01, .1352E-02, .1268E-04,-.2288E-07,-.3499E-09,
     7 .2142E+01, .4320E-03, .8474E-05, .9677E-08, .1738E-09,
     7 .2599E+01,-.3257E-03, .5604E-05, .2597E-08, .2432E-09,
     8 .4463E-03, .8307E-05,-.3433E-07, .2048E-10, .1430E-10,
     8 .4475E-03, .8326E-05,-.3442E-07, .2052E-10, .1433E-10,
     8 .4489E-03, .8346E-05,-.3451E-07, .2051E-10, .1436E-10,
     8 .4508E-03, .8365E-05,-.3463E-07, .2044E-10, .1441E-10,
     8 .4533E-03, .8385E-05,-.3476E-07, .2021E-10, .1445E-10,
     8 .4575E-03, .8405E-05,-.3499E-07, .1980E-10, .1449E-10,
     8 .4648E-03, .8424E-05,-.3534E-07, .1882E-10, .1456E-10,
     8 .4787E-03, .8444E-05,-.3595E-07, .1675E-10, .1465E-10,
     8 .5057E-03, .8794E-05,-.4191E-07,-.2605E-09, .1871E-10,
     8 .5598E-03, .8807E-05,-.4430E-07,-.2722E-09, .1905E-10,
     8 .7424E-03, .9005E-05,-.4082E-07,-.3493E-09, .1837E-10,
     8 .8900E-03, .9003E-05,-.4636E-07,-.3769E-09, .1903E-10,
     8 .1085E-02, .1098E-04,-.8920E-07,-.8430E-09, .2733E-10,
     8 .1385E-02, .1061E-04,-.8457E-07,-.2774E-09, .1826E-10,
     8 .1800E-02, .2210E-04,-.3587E-06,-.2764E-08, .6647E-10,
     8 .1551E-02, .2254E-04, .2210E-06,-.1893E-08,-.3515E-10,
     8 .1629E-02,-.1009E-04, .3580E-06, .3774E-08,-.7433E-10,
     8 .1398E-02, .2994E-04, .9413E-07,-.4945E-08,-.3497E-10,
     8 .8576E-03, .3710E-04, .1048E-05,-.6150E-08,-.1951E-09,
     8 .1219E-02, .3918E-04, .1070E-05,-.6266E-08,-.1950E-09,
     8 .1731E-02, .7905E-05, .1396E-06,-.1802E-09, .1317E-10,
     8 .2246E-02, .6806E-05, .3558E-06, .2573E-08, .1007E-10,
     8 .2169E-02, .2470E-04,-.4174E-06,-.3887E-08, .1279E-09,
     8 .4884E-03,-.4568E-05,-.6944E-07, .3013E-08, .5515E-10,
     8 .4664E-03,-.1944E-05, .1159E-06, .5154E-11,-.1063E-10,
     8 .6480E-03,-.5752E-05, .1932E-07, .9659E-09, .1025E-10,
     8 .7927E-03,-.4316E-05, .3472E-08, .2729E-09, .4591E-11,
     8 .8797E-03,-.3889E-05, .1663E-07, .2916E-09, .1433E-11 /
      DATA ((( CS3CO2GH(K,J,I), K = 1, 5), J = 1, 28), I = 3, 4) /
     9 .7072E+01, .4635E-01, .6986E-04,-.4740E-06,-.3617E-08,
     9 .7072E+01, .4635E-01, .6986E-04,-.4740E-06,-.3617E-08,
     9 .7072E+01, .4635E-01, .6986E-04,-.4740E-06,-.3617E-08,
     9 .7072E+01, .4635E-01, .6988E-04,-.4740E-06,-.3619E-08,
     9 .7073E+01, .4634E-01, .6989E-04,-.4740E-06,-.3619E-08,
     9 .7075E+01, .4633E-01, .6990E-04,-.4740E-06,-.3615E-08,
     9 .7079E+01, .4631E-01, .6997E-04,-.4744E-06,-.3616E-08,
     9 .7086E+01, .4626E-01, .7017E-04,-.4760E-06,-.3608E-08,
     9 .7102E+01, .4615E-01, .7057E-04,-.4788E-06,-.3561E-08,
     9 .7136E+01, .4586E-01, .7177E-04,-.4800E-06,-.3523E-08,
     9 .7278E+01, .4465E-01, .7563E-04,-.4760E-06,-.3449E-08,
     9 .7391E+01, .4377E-01, .7696E-04,-.4662E-06,-.3351E-08,
     9 .7564E+01, .4264E-01, .7453E-04,-.4606E-06,-.2675E-08,
     9 .7825E+01, .4096E-01, .7118E-04,-.4239E-06,-.1980E-08,
     9 .8214E+01, .3867E-01, .6192E-04,-.3572E-06,-.5044E-09,
     9 .8766E+01, .3552E-01, .5864E-04,-.2787E-06,-.1857E-09,
     9 .9541E+01, .3102E-01, .5137E-04,-.1241E-06, .9859E-09,
     9 .1056E+02, .2573E-01, .5714E-04,-.5411E-07,-.3767E-09,
     9 .1188E+02, .2028E-01, .3384E-04,-.1040E-06, .3031E-08,
     9 .1348E+02, .1391E-01, .1193E-04,-.1586E-07, .3862E-08,
     9 .1523E+02, .6412E-02, .1362E-04, .1643E-06, .8792E-09,
     9 .1690E+02,-.1419E-04, .1943E-04, .1236E-06,-.2510E-09,
     9 .1836E+02,-.5249E-02,-.2073E-04, .8310E-07, .4450E-08,
     9 .1938E+02,-.9898E-02,-.2011E-04, .1999E-06, .2000E-08,
     9 .1991E+02,-.1224E-01,-.3395E-04, .2242E-06, .2341E-08,
     9 .1993E+02,-.1305E-01,-.5027E-04, .1425E-06, .3202E-08,
     9 .1930E+02,-.1236E-01,-.5040E-04, .3909E-07, .2455E-08,
     9 .1797E+02,-.1138E-01,-.4258E-04, .1719E-07, .1772E-08,
     a .5726E+03,-.8229E+00,-.2428E-02, .7345E-05, .1473E-06,
     a .5726E+03,-.8229E+00,-.2428E-02, .7351E-05, .1474E-06,
     a .5726E+03,-.8228E+00,-.2429E-02, .7345E-05, .1475E-06,
     a .5726E+03,-.8229E+00,-.2428E-02, .7365E-05, .1473E-06,
     a .5725E+03,-.8228E+00,-.2431E-02, .7372E-05, .1479E-06,
     a .5725E+03,-.8225E+00,-.2431E-02, .7362E-05, .1479E-06,
     a .5724E+03,-.8221E+00,-.2433E-02, .7368E-05, .1479E-06,
     a .5721E+03,-.8213E+00,-.2435E-02, .7407E-05, .1474E-06,
     a .5715E+03,-.8194E+00,-.2437E-02, .7411E-05, .1464E-06,
     a .5706E+03,-.8156E+00,-.2453E-02, .7463E-05, .1464E-06,
     a .5668E+03,-.8013E+00,-.2493E-02, .7561E-05, .1434E-06,
     a .5641E+03,-.7911E+00,-.2507E-02, .7457E-05, .1413E-06,
     a .5600E+03,-.7787E+00,-.2452E-02, .7399E-05, .1288E-06,
     a .5541E+03,-.7590E+00,-.2351E-02, .6716E-05, .1114E-06,
     a .5453E+03,-.7331E+00,-.2118E-02, .5597E-05, .7853E-07,
     a .5331E+03,-.6971E+00,-.2033E-02, .4291E-05, .7528E-07,
     a .5160E+03,-.6415E+00,-.1815E-02, .1410E-05, .4675E-07,
     a .4933E+03,-.5763E+00,-.1879E-02,-.1765E-06, .7367E-07,
     a .4634E+03,-.5152E+00,-.1351E-02, .6103E-06, .3226E-09,
     a .4254E+03,-.4422E+00,-.8841E-03,-.7826E-06,-.1778E-07,
     a .3806E+03,-.3482E+00,-.9576E-03,-.3944E-05, .4819E-07,
     a .3315E+03,-.2737E+00,-.1019E-02,-.1508E-05, .4752E-07,
     a .2801E+03,-.2108E+00,-.2805E-03, .7385E-06,-.5492E-07,
     a .2314E+03,-.1434E+00,-.4725E-03,-.9929E-06,-.3489E-08,
     a .1861E+03,-.1157E+00,-.4424E-03, .5096E-06,-.5047E-08,
     a .1439E+03,-.8308E-01,-.4975E-03, .2899E-06, .1378E-07,
     a .1073E+03,-.5905E-01,-.2466E-03, .3445E-06,-.8691E-08,
     a .7825E+02,-.4257E-01,-.3012E-03, .3435E-06, .1945E-07 /
C
C----------------------------------------------------------------------C
C     SPECTRAL REGION: 2500-4200 CM^-1. GAS: H2O, CO2, CH4             C
C----------------------------------------------------------------------C
C
      COMMON /BANDS4/ GWS4(4), CS4H2O(5,18,4), CS4CO2(5,18,4),
     1                         CS4CH4(5,18,2)
      DATA GWS4 / 9.28054, 10.50363, 7.63430, 5.99443 /
      DATA ((( CS4H2O(K,J,I), K = 1, 5), J = 1, 18), I = 1, 4) /
     1 .5592E-04, .2577E-06,-.1287E-09, .2232E-11,-.7984E-15,
     1 .8267E-04, .3817E-06,-.2327E-09, .3227E-11, .3594E-14,
     1 .1222E-03, .5646E-06,-.3442E-09, .4849E-11, .3842E-14,
     1 .1806E-03, .8374E-06,-.6381E-09, .6654E-11, .2161E-13,
     1 .2667E-03, .1242E-05,-.1246E-08, .8999E-11, .8306E-13,
     1 .3936E-03, .1851E-05,-.1463E-08, .1037E-10, .5800E-13,
     1 .5805E-03, .2748E-05,-.2492E-08, .1169E-10, .1462E-12,
     1 .8538E-03, .4054E-05,-.3264E-08, .1603E-10, .1681E-12,
     1 .1253E-02, .5945E-05,-.5601E-08, .2541E-10, .4190E-12,
     1 .1832E-02, .8736E-05,-.6827E-08, .3203E-10, .4118E-12,
     1 .2664E-02, .1280E-04,-.5513E-08, .2484E-10,-.9931E-13,
     1 .3845E-02, .1845E-04,-.6037E-08, .5975E-10,-.5198E-13,
     1 .5488E-02, .2612E-04,-.2634E-08, .1229E-09,-.7183E-12,
     1 .7735E-02, .3642E-04,-.7044E-08, .1836E-09,-.1534E-12,
     1 .1074E-01, .5014E-04,-.1123E-07, .2327E-09, .2741E-12,
     1 .1459E-01, .6760E-04,-.2248E-07, .2997E-09, .3149E-11,
     1 .1939E-01, .8921E-04,-.1968E-07, .4573E-09, .5723E-11,
     1 .2519E-01, .1164E-03, .6593E-08, .4594E-09, .6798E-11,
     2 .1737E-02, .1381E-04, .2550E-07,-.9943E-10,-.9423E-13,
     2 .2477E-02, .1942E-04, .3923E-07,-.1386E-09,-.6275E-12,
     2 .3539E-02, .2735E-04, .5841E-07,-.2036E-09,-.1352E-11,
     2 .5059E-02, .3855E-04, .8579E-07,-.2984E-09,-.2417E-11,
     2 .7238E-02, .5414E-04, .1262E-06,-.3918E-09,-.3878E-11,
     2 .1036E-01, .7529E-04, .1817E-06,-.3652E-09,-.5313E-11,
     2 .1485E-01, .1047E-03, .2530E-06,-.2549E-09,-.5940E-11,
     2 .2129E-01, .1464E-03, .3520E-06,-.1399E-09,-.7909E-11,
     2 .3054E-01, .2048E-03, .4639E-06, .6434E-10,-.5954E-11,
     2 .4373E-01, .2880E-03, .6325E-06,-.1731E-10,-.8920E-11,
     2 .6239E-01, .4029E-03, .8870E-06,-.3573E-10,-.1919E-10,
     2 .8854E-01, .5581E-03, .1099E-05, .5495E-09,-.7205E-11,
     2 .1244E+00, .7685E-03, .1531E-05, .1430E-08,-.1672E-10,
     2 .1734E+00, .1060E-02, .2138E-05, .2439E-08,-.2988E-10,
     2 .2409E+00, .1459E-02, .2604E-05, .3538E-08,-.3703E-11,
     2 .3316E+00, .1982E-02, .3492E-05, .6135E-08,-.4184E-10,
     2 .4508E+00, .2678E-02, .4649E-05, .1041E-07,-.5296E-10,
     2 .6073E+00, .3657E-02, .6078E-05, .6189E-08,-.2860E-10,
     3 .3525E-01, .2661E-03, .5611E-06,-.1622E-08,-.6309E-11,
     3 .4817E-01, .3453E-03, .7108E-06,-.1433E-08,-.7615E-11,
     3 .6628E-01, .4523E-03, .9521E-06,-.1454E-08,-.2059E-10,
     3 .9160E-01, .5952E-03, .1258E-05,-.1684E-08,-.4299E-10,
     3 .1272E+00, .7837E-03, .1749E-05,-.2985E-09,-.7960E-10,
     3 .1779E+00, .1034E-02, .1986E-05, .1948E-08,-.8430E-10,
     3 .2496E+00, .1387E-02, .2600E-05, .3641E-08,-.1519E-09,
     3 .3518E+00, .1866E-02, .3210E-05, .8245E-08,-.2363E-09,
     3 .4966E+00, .2572E-02, .3548E-05, .1028E-07,-.2428E-09,
     3 .7026E+00, .3614E-02, .3911E-05, .5702E-08,-.2776E-09,
     3 .9942E+00, .5042E-02, .3764E-05, .7591E-08,-.2550E-09,
     3 .1403E+01, .7044E-02, .1659E-05, .8579E-08, .1295E-09,
     3 .1968E+01, .9904E-02,-.3554E-06,-.5085E-08, .4454E-09,
     3 .2739E+01, .1383E-01, .4070E-06,-.1132E-07, .3275E-09,
     3 .3776E+01, .1933E-01,-.2005E-05,-.4346E-07, .9674E-09,
     3 .5120E+01, .2652E-01,-.2988E-05,-.1151E-06, .9606E-09,
     3 .6760E+01, .3526E-01,-.1111E-04,-.2240E-06, .3458E-08,
     3 .8610E+01, .4511E-01, .3661E-04,-.3240E-06,-.2824E-08,
     4 .8116E+00, .5574E-02, .1987E-04, .3976E-07, .3601E-09,
     4 .1051E+01, .6453E-02, .2090E-04, .2876E-07, .2784E-09,
     4 .1381E+01, .7553E-02, .2324E-04, .1886E-07, .1448E-09,
     4 .1833E+01, .8690E-02, .2617E-04, .5986E-07, .1783E-09,
     4 .2466E+01, .1024E-01, .3010E-04, .1049E-06, .4209E-10,
     4 .3368E+01, .1209E-01, .2591E-04, .1668E-06, .5946E-09,
     4 .4609E+01, .1482E-01, .4553E-04, .1798E-06,-.2642E-08,
     4 .6369E+01, .1790E-01, .4108E-04, .2408E-06,-.1678E-08,
     4 .8861E+01, .2262E-01, .1943E-04, .2076E-06, .8241E-09,
     4 .1233E+02, .2849E-01, .2379E-05, .2958E-06, .2152E-08,
     4 .1709E+02, .3671E-01,-.1159E-04, .2052E-06, .2775E-08,
     4 .2354E+02, .4598E-01,-.4624E-04, .3609E-06, .1024E-07,
     4 .3217E+02, .5730E-01,-.1217E-03, .5018E-07, .1940E-07,
     4 .4332E+02, .6699E-01,-.1205E-03, .3148E-06, .1618E-07,
     4 .5727E+02, .8306E-01,-.6204E-04,-.2697E-06, .5528E-08,
     4 .7448E+02, .9151E-01,-.2017E-03, .1814E-06, .2361E-07,
     4 .9460E+02, .1120E+00,-.5514E-03,-.3244E-05, .9325E-07,
     4 .1175E+03, .1392E+00, .3271E-03,-.5790E-05,-.2586E-07 /
      DATA ((( CS4CO2(K,J,I), K = 1, 5), J = 1, 18), I = 1, 4) /
     1 .6261E-04, .5037E-06, .1374E-09, .1557E-10, .1901E-12,
     1 .7797E-04, .6904E-06, .5510E-08,-.7219E-11,-.8136E-12,
     1 .1001E-03, .7793E-06, .3577E-08,-.7829E-11,-.4701E-12,
     1 .1293E-03, .9082E-06, .3395E-08, .6750E-11,-.7804E-13,
     1 .1740E-03, .1221E-05, .9226E-08,-.1116E-10,-.1233E-11,
     1 .2352E-03, .1396E-05,-.7048E-09,-.1618E-10, .3818E-12,
     1 .2974E-03, .1528E-05, .1026E-07,-.3034E-12,-.1696E-11,
     1 .3921E-03, .1359E-05,-.6074E-08,-.6263E-11,-.1457E-12,
     1 .4343E-03, .1054E-06,-.3559E-08, .5859E-10, .1778E-12,
     1 .4405E-03,-.2235E-06, .2737E-08, .3772E-10,-.1463E-12,
     1 .4435E-03,-.2474E-07, .6131E-08,-.7681E-11,-.4213E-12,
     1 .4494E-03, .1261E-06, .2010E-09,-.6716E-10, .9195E-12,
     1 .4460E-03, .2448E-07, .9120E-08,-.7907E-10,-.6952E-12,
     1 .4511E-03,-.2294E-06, .6211E-08,-.2296E-10,-.3098E-12,
     1 .4536E-03,-.5014E-06, .9430E-08, .3714E-10,-.1378E-11,
     1 .4462E-03,-.4800E-06, .1249E-07, .3812E-10,-.2149E-11,
     1 .4365E-03,-.2881E-06, .8960E-08,-.1185E-10,-.1193E-11,
     1 .4384E-03,-.8316E-07, .3455E-08,-.5502E-10,-.6247E-13,
     2 .2530E-03, .3255E-05, .1424E-07, .4092E-10, .6495E-12,
     2 .3372E-03, .4202E-05, .7776E-08,-.8892E-11, .1922E-11,
     2 .4260E-03, .5392E-05, .3533E-07,-.3847E-10,-.2635E-11,
     2 .5718E-03, .6611E-05, .3602E-07,-.2931E-10,-.2490E-11,
     2 .7723E-03, .8938E-05, .3108E-07,-.1888E-09,-.1775E-11,
     2 .1021E-02, .1150E-04, .3347E-07,-.1745E-09, .1515E-12,
     2 .1377E-02, .1459E-04, .1789E-07,-.2407E-09, .4409E-11,
     2 .1802E-02, .1775E-04, .1230E-07,-.8160E-10, .1156E-10,
     2 .2325E-02, .2092E-04, .1196E-06,-.7975E-11,-.6533E-11,
     2 .3096E-02, .2188E-04, .1242E-07, .4206E-09, .1353E-10,
     2 .3963E-02, .2430E-04,-.5516E-07, .9825E-09, .3105E-10,
     2 .5052E-02, .2564E-04,-.5061E-07, .2064E-08, .3339E-10,
     2 .6190E-02, .2855E-04,-.3880E-07, .2374E-08, .3456E-10,
     2 .7486E-02, .4424E-04,-.1660E-06, .9655E-10, .6113E-10,
     2 .8037E-02, .5116E-04, .6480E-06,-.1248E-08,-.7897E-10,
     2 .9205E-02, .5148E-04,-.8676E-07, .3035E-09, .9142E-10,
     2 .1002E-01, .6223E-04, .1675E-06,-.5249E-09, .5727E-10,
     2 .1273E-01, .3489E-04,-.5331E-06, .2284E-08, .1166E-09,
     3 .1672E-02, .2189E-04, .1305E-06, .3441E-09,-.3714E-11,
     3 .2281E-02, .2848E-04, .1333E-06, .5774E-09, .2344E-11,
     3 .3127E-02, .3732E-04, .1285E-06, .7978E-09, .1024E-10,
     3 .4268E-02, .5143E-04, .2011E-06, .7393E-09, .7058E-11,
     3 .5937E-02, .6833E-04, .1559E-06, .7198E-09, .2469E-10,
     3 .8074E-02, .9408E-04, .3375E-06, .3314E-10, .5666E-11,
     3 .1121E-01, .1259E-03, .3550E-06,-.9904E-10, .2027E-10,
     3 .1564E-01, .1745E-03, .3868E-06,-.8618E-09, .3906E-10,
     3 .2184E-01, .2361E-03, .5712E-06,-.3948E-09, .4233E-10,
     3 .3087E-01, .3176E-03, .3009E-06, .1091E-08, .1374E-09,
     3 .4256E-01, .4410E-03, .1168E-05,-.2772E-09, .3595E-10,
     3 .5871E-01, .6256E-03, .2174E-05,-.5481E-08,-.9776E-10,
     3 .7925E-01, .8404E-03, .3116E-05,-.5324E-08,-.1592E-09,
     3 .1055E+00, .1152E-02, .5549E-05,-.1392E-07,-.4609E-09,
     3 .1402E+00, .1442E-02, .4969E-05,-.9516E-08,-.1955E-09,
     3 .1763E+00, .1727E-02, .5510E-05, .9812E-08, .1655E-09,
     3 .2140E+00, .2214E-02, .1284E-04,-.6666E-08,-.8985E-09,
     3 .2469E+00, .2252E-02, .7551E-05, .3861E-07, .6569E-09,
     4 .2718E-02, .3544E-04, .1086E-06,-.5151E-09, .2852E-13,
     4 .3695E-02, .4840E-04, .1913E-06,-.7000E-09,-.6782E-11,
     4 .5107E-02, .6631E-04, .2449E-06,-.8813E-09,-.4336E-11,
     4 .7083E-02, .9095E-04, .3115E-06,-.1345E-08,-.3109E-11,
     4 .9688E-02, .1250E-03, .6642E-06,-.1337E-08,-.3880E-10,
     4 .1367E-01, .1660E-03, .7248E-06, .1952E-10,-.1212E-10,
     4 .1924E-01, .2299E-03, .9732E-06, .1382E-10,-.1529E-10,
     4 .2692E-01, .3207E-03, .1658E-05,-.2897E-09,-.8347E-10,
     4 .3765E-01, .4610E-03, .2371E-05,-.2885E-08,-.1218E-09,
     4 .5183E-01, .6597E-03, .4395E-05,-.7073E-08,-.3486E-09,
     4 .7253E-01, .8844E-03, .4624E-05,-.3179E-08,-.2014E-09,
     4 .9962E-01, .1151E-02, .5207E-05, .1879E-08,-.7895E-10,
     4 .1364E+00, .1544E-02, .6446E-05, .2586E-08, .1470E-10,
     4 .1855E+00, .1982E-02, .6654E-05, .1394E-07, .3633E-09,
     4 .2462E+00, .2578E-02, .8716E-05, .2148E-07, .5396E-09,
     4 .3209E+00, .3333E-02, .1851E-04,-.2751E-10,-.1142E-08,
     4 .4063E+00, .3815E-02, .1547E-04, .6847E-07, .3526E-09,
     4 .5056E+00, .5021E-02, .2618E-04, .3087E-07,-.1368E-08 /
      DATA ((( CS4CH4(K,J,I), K = 1, 5), J = 1, 18), I = 1, 2) /
     1 .2404E+01, .1922E-01, .1589E-04, .3120E-07, .3130E-08,
     1 .2622E+01, .1985E-01, .8354E-05, .1428E-06, .5497E-08,
     1 .2897E+01, .2082E-01, .1064E-04, .1411E-06, .5144E-08,
     1 .3235E+01, .2195E-01, .2651E-04, .1600E-06, .2624E-08,
     1 .3671E+01, .2315E-01, .4736E-04, .1881E-06,-.1044E-08,
     1 .4261E+01, .2451E-01,-.1116E-06, .2250E-06, .7420E-08,
     1 .4978E+01, .2665E-01, .1400E-04, .1532E-06, .4900E-08,
     1 .5912E+01, .2872E-01, .4555E-05, .1568E-06, .7155E-08,
     1 .7181E+01, .3240E-01, .2954E-04,-.9724E-07, .4783E-09,
     1 .8793E+01, .3380E-01,-.1014E-04, .1399E-06, .8253E-08,
     1 .1079E+02, .3466E-01,-.5630E-04, .6485E-06, .1775E-07,
     1 .1327E+02, .3778E-01, .1983E-05, .4965E-06, .2306E-09,
     1 .1636E+02, .4285E-01,-.6950E-04, .3729E-06, .1141E-07,
     1 .2003E+02, .4600E-01,-.3880E-04, .9760E-06,-.2602E-09,
     1 .2400E+02, .4897E-01,-.1700E-03, .8973E-06, .2581E-07,
     1 .2827E+02, .5257E-01,-.1923E-04, .6579E-06,-.5314E-08,
     1 .3263E+02, .5451E-01,-.1602E-03, .3972E-06, .3653E-07,
     1 .3758E+02, .6678E-01,-.1595E-03,-.8379E-06, .1569E-07,
     2 .1606E+01, .1498E-01, .4753E-04,-.1993E-06,-.1673E-08,
     2 .1784E+01, .1619E-01, .4839E-04,-.2824E-06,-.2573E-08,
     2 .1978E+01, .1730E-01, .4199E-04,-.3251E-06,-.7562E-09,
     2 .2220E+01, .1803E-01, .5258E-04,-.2383E-06,-.1594E-08,
     2 .2543E+01, .1805E-01, .3960E-04, .8035E-07, .2780E-08,
     2 .2928E+01, .1959E-01, .8107E-04, .7686E-07,-.3476E-08,
     2 .3463E+01, .2116E-01, .7426E-04, .2320E-06,-.5591E-10,
     2 .4144E+01, .2401E-01, .1235E-03, .9916E-07,-.8085E-08,
     2 .5009E+01, .2522E-01, .8870E-04, .4139E-06, .4810E-09,
     2 .6113E+01, .2897E-01, .1353E-03, .3888E-06,-.6142E-08,
     2 .7610E+01, .3515E-01, .9484E-04,-.1455E-06,-.5165E-09,
     2 .9302E+01, .3891E-01, .1254E-03, .1830E-06, .5078E-08,
     2 .1161E+02, .4329E-01, .1841E-03, .3182E-06,-.6170E-08,
     2 .1437E+02, .5487E-01, .1053E-03,-.1140E-05, .1311E-07,
     2 .1763E+02, .5832E-01, .3375E-03,-.5222E-06,-.2512E-07,
     2 .2172E+02, .6013E-01, .4312E-04, .4326E-06, .3830E-07,
     2 .2692E+02, .6939E-01, .1721E-03, .3783E-06, .1288E-07,
     2 .3373E+02, .7861E-01, .4114E-03, .1664E-05,-.6396E-09 /
C
      COMMON /BANDS4GH/ GWS4GH(9), CS4H2OGH(5,28,6), CS4CO2GH(5,28,6),
     1                             CS4CH4GH(5,28)
      DATA GWS4GH / 2.28900, 2.80913, .99438, .67932, .13964,
     1               .07300,  .02357, .01957, .01083 /
      DATA ((( CS4H2OGH(K,J,I), K = 1, 5), J = 1, 28), I = 1, 3) /
     5 .2726E+03, .2168E+01,-.3522E-03, .1056E-04, .4138E-06,
     5 .2722E+03, .2206E+01, .5696E-03, .3559E-05, .2220E-06,
     5 .2719E+03, .2320E+01,-.1043E-02,-.2148E-04, .5757E-06,
     5 .2729E+03, .2387E+01,-.4270E-02,-.3509E-04, .1188E-05,
     5 .2707E+03, .2328E+01,-.2208E-02,-.2483E-04, .8379E-06,
     5 .2634E+03, .2324E+01, .4172E-02,-.2349E-04,-.2219E-06,
     5 .2588E+03, .2309E+01, .8432E-02,-.2680E-04,-.9710E-06,
     5 .2561E+03, .2188E+01, .6812E-02,-.5401E-05,-.5616E-06,
     5 .2616E+03, .2152E+01,-.6678E-03,-.1263E-04, .6462E-06,
     5 .2677E+03, .2117E+01,-.2349E-02,-.1610E-04, .7640E-06,
     5 .2754E+03, .2162E+01, .4339E-02,-.2899E-04,-.2183E-06,
     5 .2860E+03, .2122E+01, .2513E-02,-.3093E-04, .1283E-06,
     5 .3020E+03, .2018E+01,-.1983E-04,-.2233E-04, .6163E-06,
     5 .3259E+03, .1856E+01,-.1469E-02,-.4811E-05, .9509E-06,
     5 .3532E+03, .1854E+01, .6609E-02,-.9267E-05,-.5244E-06,
     5 .4054E+03, .1731E+01,-.2606E-02, .1579E-04, .9374E-06,
     5 .4550E+03, .1798E+01, .4725E-02,-.9628E-05,-.1942E-06,
     5 .5207E+03, .1531E+01, .1205E-01, .7999E-05,-.1385E-05,
     5 .6288E+03, .1251E+01, .1979E-02, .5124E-04, .1916E-06,
     5 .7589E+03, .1451E+01, .5982E-02, .1601E-04,-.8078E-06,
     5 .9164E+03, .1546E+01,-.6987E-02,-.1463E-04, .1418E-05,
     5 .1084E+04, .1357E+01, .8108E-02, .3433E-04,-.1229E-05,
     5 .1302E+04, .1082E+01,-.5348E-02,-.5049E-05, .1016E-05,
     5 .1520E+04, .5612E+00,-.2590E-02, .2664E-04, .7822E-06,
     5 .1756E+04, .3074E+00,-.3550E-02,-.2608E-05, .8115E-06,
     5 .1957E+04, .2334E-01, .4211E-02,-.1794E-04,-.4082E-06,
     5 .2114E+04,-.1073E+01,-.2399E-01, .2688E-04, .4753E-05,
     5 .2162E+04,-.1115E+01, .1787E-01,-.8036E-04,-.1357E-05,
     6 .8795E+01, .7441E-01, .2505E-03,-.5124E-06,-.5407E-08,
     6 .8989E+01, .7027E-01, .6147E-04, .4680E-06, .3332E-07,
     6 .9053E+01, .4932E-01, .5292E-03, .5519E-05,-.6462E-07,
     6 .9031E+01, .4000E-01, .1102E-02, .7722E-05,-.1744E-06,
     6 .9681E+01, .5050E-01, .4492E-03, .6094E-05,-.5828E-07,
     6 .1124E+02, .5599E-01,-.7034E-03, .5280E-05, .1298E-06,
     6 .1194E+02, .5853E-01,-.1178E-02, .5596E-05, .2129E-06,
     6 .1263E+02, .8399E-01,-.1040E-02, .9051E-06, .1716E-06,
     6 .1231E+02, .8346E-01, .1726E-03, .3348E-05,-.2068E-07,
     6 .1230E+02, .8893E-01, .3077E-03, .2845E-05,-.1609E-07,
     6 .1444E+02, .6203E-01,-.1181E-02, .7041E-05, .2244E-06,
     6 .1490E+02, .6382E-01,-.5504E-03, .7373E-05, .1169E-06,
     6 .1587E+02, .8529E-01, .4957E-03, .2055E-05,-.7955E-07,
     6 .1762E+02, .8863E-01, .6113E-03, .1139E-05,-.1248E-06,
     6 .2090E+02, .7270E-01,-.2111E-03, .3077E-05, .8283E-08,
     6 .2429E+02, .7391E-01,-.4123E-03, .3008E-06, .5517E-07,
     6 .2976E+02, .7631E-01,-.2990E-03,-.1688E-05,-.2329E-07,
     6 .3674E+02, .1087E+00,-.8494E-03,-.7448E-05, .7146E-07,
     6 .4285E+02, .8397E-01, .9991E-03,-.4659E-05,-.2454E-06,
     6 .5283E+02, .5742E-01,-.1447E-02,-.2204E-05, .1853E-06,
     6 .6551E+02, .4685E-01,-.8482E-03, .1464E-05, .6698E-07,
     6 .8203E+02, .4464E-01,-.2187E-02,-.3503E-05, .2171E-06,
     6 .9722E+02, .7826E-01, .5294E-03, .3086E-05,-.2331E-06,
     6 .1176E+03, .1019E+00,-.6740E-03, .2251E-05,-.7449E-07,
     6 .1338E+03, .1051E+00,-.2564E-02, .6083E-05, .2898E-06,
     6 .1458E+03, .1311E+00,-.4765E-03,-.7380E-06,-.8958E-07,
     6 .1588E+03, .1813E+00, .4481E-02, .9117E-05,-.9424E-06,
     6 .1820E+03, .1178E+00,-.4064E-02, .2665E-04, .3025E-06,
     7 .1350E+02, .1258E+00, .9812E-03,-.7326E-05,-.2142E-06,
     7 .1350E+02, .1220E+00, .1098E-02,-.6730E-05,-.2355E-06,
     7 .1386E+02, .1317E+00, .4720E-03,-.9470E-05,-.1116E-06,
     7 .1389E+02, .1286E+00, .3827E-03,-.9154E-05,-.9002E-07,
     7 .1347E+02, .1296E+00, .1357E-02,-.9382E-05,-.2710E-06,
     7 .1305E+02, .1165E+00, .1601E-02,-.7478E-05,-.2997E-06,
     7 .1400E+02, .1185E+00, .8303E-03,-.6943E-05,-.1579E-06,
     7 .1441E+02, .1034E+00, .1248E-02,-.3646E-05,-.2422E-06,
     7 .1468E+02, .1208E+00, .1120E-02,-.6506E-05,-.2160E-06,
     7 .1486E+02, .1132E+00, .1664E-02,-.3621E-05,-.2897E-06,
     7 .1621E+02, .1545E+00, .2947E-02,-.1319E-04,-.5434E-06,
     7 .1791E+02, .1449E+00, .2009E-02,-.1273E-04,-.3834E-06,
     7 .1990E+02, .9078E-01,-.6946E-04, .2927E-05, .1279E-08,
     7 .2006E+02, .1190E+00, .1243E-02,-.1796E-05,-.2076E-06,
     7 .2288E+02, .1178E+00,-.1404E-02,-.2334E-05, .2517E-06,
     7 .2193E+02, .1351E+00, .2600E-02,-.7836E-05,-.4556E-06,
     7 .2787E+02, .9534E-01,-.4731E-03, .3187E-05, .7110E-07,
     7 .3671E+02, .1126E+00,-.3526E-02, .9112E-05, .5673E-06,
     7 .4533E+02, .2142E+00,-.5035E-02,-.6585E-05, .7927E-06,
     7 .5077E+02, .1978E+00,-.9447E-04,-.7448E-05,-.4328E-08,
     7 .5932E+02, .1885E+00, .3386E-02,-.4461E-05,-.6316E-06,
     7 .7326E+02, .2391E+00, .8387E-03,-.6350E-05, .1854E-07,
     7 .8862E+02, .2698E+00, .9622E-03,-.9800E-05,-.1311E-06,
     7 .1041E+03, .3608E+00, .3039E-03,-.2158E-04,-.1381E-06,
     7 .1156E+03, .4188E+00, .6243E-02,-.2891E-04,-.1261E-05,
     7 .1324E+03, .2857E+00, .1503E-02, .6821E-06,-.2322E-06,
     7 .1474E+03, .3699E+00,-.2621E-03,-.2587E-04,-.2427E-07,
     7 .1640E+03, .3430E+00, .1358E-03,-.1905E-04,-.5765E-07 /
      DATA ((( CS4H2OGH(K,J,I), K = 1, 5), J = 1, 28), I = 4, 6) /
     8 .2220E+00, .4792E-03,-.1820E-04,-.3642E-07, .1455E-09,
     8 .2220E+00, .4792E-03,-.1820E-04,-.3642E-07, .1456E-09,
     8 .2221E+00, .4792E-03,-.1820E-04,-.3643E-07, .1455E-09,
     8 .2221E+00, .4792E-03,-.1820E-04,-.3645E-07, .1452E-09,
     8 .2221E+00, .4792E-03,-.1820E-04,-.3649E-07, .1451E-09,
     8 .2222E+00, .5192E-03,-.1877E-04,-.6892E-07, .6072E-09,
     8 .2224E+00, .4800E-03,-.1821E-04,-.3673E-07, .1438E-09,
     8 .2231E+00, .4804E-03,-.1866E-04,-.3703E-07, .2148E-09,
     8 .1755E+00, .4812E-03, .3045E-04,-.3659E-07,-.7792E-08,
     8 .1769E+00, .3882E-03, .3168E-04, .3942E-07,-.8912E-08,
     8 .1924E+00,-.1481E-02,-.2282E-04, .3498E-06, .2478E-08,
     8 .2002E+00,-.1087E-02,-.1956E-04, .2681E-06, .1574E-08,
     8 .1339E+00,-.1328E-03, .1769E-04, .1060E-06,-.2596E-08,
     8 .1413E+00, .8508E-03,-.7776E-05,-.8956E-07, .2247E-08,
     8 .2017E+00, .7228E-03,-.5996E-04,-.7324E-07, .1063E-07,
     8 .2177E+00, .2714E-03,-.1114E-04, .6768E-07, .1524E-08,
     8 .2314E+00, .3160E-03,-.1386E-04,-.2056E-07, .3066E-08,
     8 .2082E+00,-.1928E-03, .7700E-05, .4832E-07,-.1028E-08,
     8 .1994E+00,-.1195E-02, .5224E-04, .2436E-06,-.9332E-08,
     8 .2374E+00,-.1360E-02, .3485E-04, .2686E-06,-.6628E-08,
     8 .2816E+00,-.1370E-02, .4352E-04, .2655E-06,-.8572E-08,
     8 .2897E+00, .4256E-03,-.2156E-04,-.3002E-06, .3452E-08,
     8 .1988E+00, .2117E-03, .5660E-04,-.2571E-06,-.8588E-08,
     8 .1859E+00, .3160E-03,-.3988E-05, .8252E-08, .3346E-09,
     8 .1848E+00, .3701E-03,-.1159E-05,-.3270E-08,-.1567E-09,
     8 .1856E+00, .3254E-03,-.4432E-06, .4024E-08,-.1848E-09,
     8 .1875E+00, .7164E-04, .1031E-04, .4956E-07,-.2232E-08,
     8 .1900E+00, .4556E-03, .7244E-05,-.3722E-07,-.1551E-08,
     a .2581E+05, .1110E+03,-.3482E+00,-.8481E-03, .7094E-05,
     a .2581E+05, .1110E+03,-.3483E+00,-.8483E-03, .7101E-05,
     a .2582E+05, .1110E+03,-.3482E+00,-.8485E-03, .7088E-05,
     a .2582E+05, .1110E+03,-.3481E+00,-.8485E-03, .7088E-05,
     a .2582E+05, .1110E+03,-.3481E+00,-.8491E-03, .7092E-05,
     a .2582E+05, .1109E+03,-.3480E+00,-.8498E-03, .7088E-05,
     a .2582E+05, .1109E+03,-.3476E+00,-.8508E-03, .7047E-05,
     a .2584E+05, .1107E+03,-.3469E+00,-.8527E-03, .6984E-05,
     a .2586E+05, .1105E+03,-.3451E+00,-.8549E-03, .6805E-05,
     a .2590E+05, .1100E+03,-.3416E+00,-.8621E-03, .6497E-05,
     a .2606E+05, .1078E+03,-.3312E+00,-.8712E-03, .5862E-05,
     a .2618E+05, .1062E+03,-.3231E+00,-.8886E-03, .5577E-05,
     a .2635E+05, .1039E+03,-.3109E+00,-.9216E-03, .5333E-05,
     a .2660E+05, .1004E+03,-.2949E+00,-.1006E-02, .5581E-05,
     a .2697E+05, .9467E+02,-.2786E+00,-.9339E-03, .5645E-05,
     a .2746E+05, .8771E+02,-.2791E+00,-.9254E-03, .1031E-04,
     a .2803E+05, .7802E+02,-.2157E+00,-.9536E-03, .6434E-05,
     a .2884E+05, .6402E+02,-.1774E+00,-.7108E-03, .8067E-05,
     a .2981E+05, .4686E+02,-.9132E-01, .2758E-03,-.5779E-05,
     a .3067E+05, .3637E+02,-.2183E+00,-.6677E-03, .3057E-04,
     a .3129E+05, .1674E+02,-.6318E-01,-.2634E-03, .9962E-05,
     a .3120E+05,-.1628E+01, .1538E+00, .2571E-03,-.2915E-04,
     a .3023E+05,-.8406E+01, .1257E+00,-.2312E-03,-.2068E-04,
     a .2829E+05,-.1159E+02,-.5485E-02,-.5082E-03, .7188E-05,
     a .2542E+05,-.1543E+02, .6479E-01,-.3855E-03, .6156E-06,
     a .2223E+05,-.1859E+02, .7944E-01,-.8370E-04,-.2836E-05,
     a .1900E+05,-.1906E+02, .7370E-01, .1169E-04,-.2807E-05,
     a .1573E+05,-.1866E+02, .8298E-01, .9688E-04,-.6710E-05,
     2 .2764E+06,-.7481E+03, .1670E+01, .1323E-02,-.4270E-04,
     2 .2764E+06,-.7481E+03, .1672E+01, .1325E-02,-.4298E-04,
     2 .2764E+06,-.7481E+03, .1670E+01, .1328E-02,-.4279E-04,
     2 .2764E+06,-.7480E+03, .1671E+01, .1331E-02,-.4298E-04,
     2 .2764E+06,-.7479E+03, .1670E+01, .1329E-02,-.4293E-04,
     2 .2764E+06,-.7477E+03, .1669E+01, .1334E-02,-.4286E-04,
     2 .2763E+06,-.7472E+03, .1665E+01, .1352E-02,-.4263E-04,
     2 .2761E+06,-.7462E+03, .1661E+01, .1371E-02,-.4276E-04,
     2 .2758E+06,-.7442E+03, .1650E+01, .1418E-02,-.4231E-04,
     2 .2751E+06,-.7402E+03, .1630E+01, .1520E-02,-.4198E-04,
     2 .2724E+06,-.7251E+03, .1556E+01, .1831E-02,-.4071E-04,
     2 .2704E+06,-.7140E+03, .1510E+01, .2025E-02,-.4104E-04,
     2 .2676E+06,-.6980E+03, .1443E+01, .2309E-02,-.4147E-04,
     2 .2634E+06,-.6749E+03, .1342E+01, .2605E-02,-.3914E-04,
     2 .2575E+06,-.6412E+03, .1289E+01, .2266E-02,-.4164E-04,
     2 .2494E+06,-.6022E+03, .1275E+01, .2521E-02,-.5625E-04,
     2 .2388E+06,-.5545E+03, .9650E+00, .3021E-02,-.3194E-04,
     2 .2247E+06,-.4850E+03, .8193E+00, .1703E-02,-.3721E-04,
     2 .2066E+06,-.4023E+03, .6581E+00,-.8701E-03,-.2067E-04,
     2 .1855E+06,-.3298E+03, .5099E+00,-.1223E-02,-.1152E-04,
     2 .1615E+06,-.2532E+03, .9934E-01,-.2178E-02, .4408E-04,
     2 .1364E+06,-.1811E+03,-.4081E+00,-.3698E-02, .1433E-03,
     2 .1106E+06,-.1494E+03, .4278E+00,-.1825E-02, .3201E-05,
     2 .8718E+05,-.1274E+03, .4363E+00,-.1988E-03,-.1005E-04,
     2 .6665E+05,-.9968E+02, .4462E+00, .2184E-03,-.2015E-04,
     2 .4978E+05,-.6806E+02, .1529E+00,-.8756E-03, .1803E-04,
     2 .3637E+05,-.4709E+02, .1173E+00,-.1442E-02, .1685E-04,
     2 .2631E+05,-.3286E+02, .1332E+00,-.1131E-02, .6958E-05 /
      DATA ((( CS4CO2GH(K,J,I), K = 1, 5), J = 1, 28), I = 1, 3) /
     5-.3586E-04, .3917E-06, .4169E-07, .9024E-10,-.5660E-11,
     5-.8123E-05, .4542E-06, .3356E-07, .7622E-10,-.4055E-11,
     5 .3542E-04, .8499E-06, .1115E-07, .3252E-10,-.4685E-12,
     5 .5783E-04, .1176E-05, .1073E-07, .2679E-10,-.1933E-12,
     5 .9229E-04, .1592E-05, .9939E-08, .5081E-10, .5633E-12,
     5 .1479E-03, .2427E-05, .1937E-07, .1526E-10,-.1011E-11,
     5 .2523E-03, .3529E-05, .2235E-07, .5476E-10,-.8509E-12,
     5 .4245E-03, .5414E-05, .2106E-07, .1319E-09, .9150E-12,
     5 .7442E-03, .8649E-05, .3172E-08, .2052E-09, .5988E-11,
     5 .1244E-02, .1477E-04, .5220E-07, .5693E-10,-.8819E-12,
     5 .3071E-02, .3552E-04, .4904E-08,-.5837E-09, .1015E-10,
     5 .4289E-02, .5078E-04, .5595E-07,-.1064E-08, .4500E-11,
     5 .6031E-02, .7573E-04, .1599E-06,-.2316E-08,-.5472E-11,
     5 .8409E-02, .1010E-03, .2290E-06,-.1712E-08, .3500E-11,
     5 .1178E-01, .1345E-03, .1636E-06,-.3684E-09, .4749E-10,
     5 .1646E-01, .1936E-03, .6497E-06,-.1793E-08,-.2176E-10,
     5 .2284E-01, .2665E-03, .1479E-05, .2175E-09,-.1045E-09,
     5 .3222E-01, .3552E-03, .1983E-05, .3634E-08,-.1246E-09,
     5 .4591E-01, .4742E-03, .1590E-05, .8026E-08, .1438E-10,
     5 .6593E-01, .6249E-03, .9761E-06, .1699E-07, .1762E-09,
     5 .9540E-01, .9347E-03, .1569E-05, .1082E-07, .1065E-09,
     5 .1386E+00, .1433E-02,-.4810E-06,-.1224E-08, .6702E-09,
     5 .1988E+00, .2207E-02, .2960E-06,-.4062E-07, .6399E-09,
     5 .2802E+00, .3049E-02, .1754E-05,-.3595E-07, .8203E-09,
     5 .3788E+00, .4419E-02, .1878E-04,-.7727E-07,-.1429E-08,
     5 .5249E+00, .5667E-02, .8574E-05,-.3228E-07, .1549E-08,
     5 .6927E+00, .7259E-02, .6496E-05,-.9997E-07, .1757E-08,
     5 .8476E+00, .1036E-01, .3215E-04,-.2214E-06,-.1200E-09,
     6 .2207E-01, .6691E-03, .8401E-05, .5634E-07, .1958E-09,
     6 .2268E-01, .6711E-03, .8375E-05, .5657E-07, .2003E-09,
     6 .2388E-01, .6739E-03, .8340E-05, .5703E-07, .1991E-09,
     6 .2602E-01, .6883E-03, .8201E-05, .5621E-07, .2078E-09,
     6 .2948E-01, .7206E-03, .8133E-05, .5407E-07, .2036E-09,
     6 .3491E-01, .7709E-03, .8185E-05, .5265E-07, .2006E-09,
     6 .4365E-01, .8506E-03, .8392E-05, .5094E-07, .1897E-09,
     6 .5788E-01, .9716E-03, .8844E-05, .4943E-07, .1636E-09,
     6 .8180E-01, .1158E-02, .9439E-05, .4762E-07, .1610E-09,
     6 .1223E+00, .1449E-02, .1069E-04, .4494E-07, .1200E-09,
     6 .2471E+00, .2190E-02, .1432E-04, .4384E-07,-.1022E-10,
     6 .3302E+00, .2595E-02, .1609E-04, .4974E-07,-.3309E-10,
     6 .4469E+00, .3122E-02, .1863E-04, .5974E-07,-.1200E-09,
     6 .6118E+00, .3814E-02, .2179E-04, .7324E-07,-.2581E-09,
     6 .8465E+00, .4770E-02, .2402E-04, .7915E-07,-.1739E-09,
     6 .1180E+01, .6011E-02, .2748E-04, .8669E-07,-.2744E-09,
     6 .1660E+01, .7618E-02, .2701E-04, .1128E-06, .3677E-09,
     6 .2347E+01, .1006E-01, .2664E-04, .9207E-07, .1119E-08,
     6 .3314E+01, .1351E-01, .3423E-04, .8469E-08, .8782E-09,
     6 .4667E+01, .1770E-01, .4059E-04,-.6268E-07, .1788E-08,
     6 .6533E+01, .2290E-01, .7620E-04,-.1053E-06,-.1551E-08,
     6 .9116E+01, .2854E-01, .1024E-03,-.1878E-07,-.2669E-08,
     6 .1261E+02, .3446E-01, .1122E-03, .2696E-06,-.2715E-08,
     6 .1730E+02, .4271E-01, .7788E-04, .5291E-06, .1076E-08,
     6 .2357E+02, .5628E-01, .2377E-04, .2886E-06, .4933E-08,
     6 .3162E+02, .7478E-01,-.4708E-05,-.6628E-06, .1230E-07,
     6 .4170E+02, .8601E-01, .1364E-03,-.6147E-06,-.1891E-08,
     6 .5490E+02, .8546E-01, .9435E-04, .1279E-06, .3128E-08,
     7 .4975E+01, .9627E-01, .6755E-03, .1113E-05,-.9051E-08,
     7 .4993E+01, .9660E-01, .6780E-03, .1116E-05,-.9090E-08,
     7 .5013E+01, .9692E-01, .6806E-03, .1118E-05,-.9121E-08,
     7 .5038E+01, .9717E-01, .6837E-03, .1124E-05,-.9201E-08,
     7 .5072E+01, .9733E-01, .6870E-03, .1132E-05,-.9274E-08,
     7 .5130E+01, .9740E-01, .6902E-03, .1135E-05,-.9347E-08,
     7 .5238E+01, .9723E-01, .6931E-03, .1139E-05,-.9353E-08,
     7 .5449E+01, .9679E-01, .6941E-03, .1159E-05,-.9359E-08,
     7 .5861E+01, .9618E-01, .6908E-03, .1210E-05,-.9779E-08,
     7 .6659E+01, .9633E-01, .6574E-03, .1196E-05,-.7066E-08,
     7 .9172E+01, .1026E+00, .6102E-03, .6211E-06,-.6252E-08,
     7 .1077E+02, .1073E+00, .6224E-03, .3816E-06,-.9373E-08,
     7 .1293E+02, .1124E+00, .6397E-03, .3321E-06,-.1199E-07,
     7 .1584E+02, .1181E+00, .6475E-03, .4139E-06,-.1133E-07,
     7 .1982E+02, .1231E+00, .5717E-03, .9387E-06, .6094E-08,
     7 .2509E+02, .1309E+00, .7233E-03, .1314E-05,-.1761E-07,
     7 .3239E+02, .1440E+00, .6915E-03, .1293E-05,-.9448E-08,
     7 .4240E+02, .1572E+00, .6223E-03, .2577E-05,-.6983E-08,
     7 .5570E+02, .1873E+00, .5946E-03, .1922E-05,-.1659E-07,
     7 .7348E+02, .2291E+00, .4506E-03, .6033E-06,-.1865E-07,
     7 .9553E+02, .2809E+00, .2253E-03,-.2824E-05, .4360E-07,
     7 .1222E+03, .3235E+00, .6616E-03,-.6892E-05, .4119E-07,
     7 .1561E+03, .2815E+00, .1219E-02,-.3312E-05, .7424E-08,
     7 .2002E+03, .1965E+00, .1423E-02,-.2842E-05,-.2594E-07,
     7 .2539E+03, .7295E-01, .1355E-02,-.6921E-06,-.2890E-07,
     7 .3142E+03,-.6172E-01, .6220E-03, .1974E-05, .2867E-07,
     7 .3717E+03,-.1545E+00, .2074E-03,-.9001E-06, .1705E-07,
     7 .4160E+03,-.2249E+00,-.1259E-03,-.9622E-06, .3930E-07 /
      DATA ((( CS4CO2GH(K,J,I), K = 1, 5), J = 1, 28), I = 4, 6) /
     9 .3783E+03, .4530E+01, .1476E-01,-.1375E-04,-.2371E-06,
     9 .3796E+03, .4545E+01, .1481E-01,-.1380E-04,-.2380E-06,
     9 .3810E+03, .4560E+01, .1486E-01,-.1384E-04,-.2387E-06,
     9 .3823E+03, .4575E+01, .1491E-01,-.1390E-04,-.2396E-06,
     9 .3837E+03, .4590E+01, .1496E-01,-.1394E-04,-.2403E-06,
     9 .3853E+03, .4604E+01, .1501E-01,-.1397E-04,-.2407E-06,
     9 .3871E+03, .4618E+01, .1505E-01,-.1403E-04,-.2415E-06,
     9 .3895E+03, .4628E+01, .1513E-01,-.1412E-04,-.2444E-06,
     9 .3933E+03, .4629E+01, .1526E-01,-.1431E-04,-.2472E-06,
     9 .4001E+03, .4609E+01, .1557E-01,-.1388E-04,-.2691E-06,
     9 .4251E+03, .4471E+01, .1562E-01,-.9185E-05,-.2451E-06,
     9 .4451E+03, .4397E+01, .1525E-01,-.7819E-05,-.1994E-06,
     9 .4744E+03, .4338E+01, .1414E-01,-.1445E-04,-.4886E-07,
     9 .5167E+03, .4289E+01, .1234E-01,-.2705E-04, .2173E-06,
     9 .5752E+03, .4246E+01, .1214E-01,-.4528E-04, .1883E-06,
     9 .6570E+03, .4078E+01, .1300E-01,-.4486E-04,-.7793E-07,
     9 .7711E+03, .3761E+01, .1330E-01,-.3409E-04,-.1904E-06,
     9 .9306E+03, .3268E+01, .1076E-01,-.2771E-04, .1199E-06,
     9 .1143E+04, .2625E+01, .6926E-02,-.5638E-04, .8222E-06,
     9 .1407E+04, .1424E+01, .1489E-01,-.6905E-05,-.9028E-06,
     9 .1713E+04, .4176E-01, .2052E-01, .4325E-04,-.2193E-05,
     9 .2065E+04,-.1299E+01,-.2324E-03, .8443E-04, .3259E-06,
     9 .2376E+04,-.1965E+01,-.1994E-01, .7019E-04, .2637E-05,
     9 .2523E+04,-.2029E+01, .1590E-02, .4053E-04,-.9654E-06,
     9 .2516E+04,-.1988E+01,-.4558E-02, .3569E-04, .1418E-06,
     9 .2347E+04,-.1707E+01,-.1043E-01, .3015E-04, .7956E-06,
     9 .2039E+04,-.1174E+01,-.5731E-02,-.9009E-05,-.6796E-08,
     9 .1677E+04,-.9375E+00,-.4102E-02, .9802E-06, .2382E-07,
     1 .1054E+05, .3739E+02,-.2153E+00, .6276E-04, .5170E-05,
     1 .1054E+05, .3739E+02,-.2153E+00, .6278E-04, .5173E-05,
     1 .1054E+05, .3739E+02,-.2153E+00, .6293E-04, .5175E-05,
     1 .1054E+05, .3738E+02,-.2153E+00, .6275E-04, .5172E-05,
     1 .1054E+05, .3737E+02,-.2152E+00, .6297E-04, .5175E-05,
     1 .1054E+05, .3735E+02,-.2152E+00, .6310E-04, .5179E-05,
     1 .1054E+05, .3731E+02,-.2150E+00, .6291E-04, .5170E-05,
     1 .1055E+05, .3723E+02,-.2146E+00, .6173E-04, .5167E-05,
     1 .1055E+05, .3707E+02,-.2139E+00, .5855E-04, .5181E-05,
     1 .1057E+05, .3675E+02,-.2122E+00, .5106E-04, .5185E-05,
     1 .1063E+05, .3540E+02,-.2016E+00, .3207E-04, .4620E-05,
     1 .1067E+05, .3428E+02,-.1933E+00, .2658E-04, .4209E-05,
     1 .1075E+05, .3247E+02,-.1831E+00, .4269E-04, .3981E-05,
     1 .1086E+05, .2960E+02,-.1644E+00, .1152E-03, .2433E-05,
     1 .1104E+05, .2555E+02,-.1477E+00, .2011E-03, .1865E-05,
     1 .1130E+05, .2101E+02,-.1422E+00, .1256E-03, .4224E-05,
     1 .1163E+05, .1545E+02,-.1330E+00, .4341E-05, .6425E-05,
     1 .1195E+05, .9180E+01,-.6567E-01,-.9699E-04,-.1151E-05,
     1 .1225E+05, .2419E+01,-.2109E-01,-.1230E-03,-.3843E-05,
     1 .1243E+05,-.3630E+01,-.3049E-01,-.1076E-04,-.8186E-06,
     1 .1218E+05,-.6972E+01,-.2846E-01, .2003E-04,-.1210E-05,
     1 .1134E+05,-.8201E+01, .1435E-03, .1053E-03,-.4718E-05,
     1 .9937E+04,-.6245E+01,-.2388E-01, .2701E-04, .3522E-06,
     1 .8154E+04,-.4637E+01,-.2236E-01,-.3785E-04, .1863E-06,
     1 .6409E+04,-.3945E+01,-.2824E-01, .3074E-04, .1093E-05,
     1 .4894E+04,-.2857E+01,-.2079E-01, .2944E-04, .9947E-06,
     1 .3685E+04,-.2131E+01,-.9061E-02, .3638E-04,-.6150E-07,
     1 .2747E+04,-.1420E+01,-.1024E-01, .2575E-05, .5252E-06,
     3 .7181E+05,-.2174E+03, .2331E+00, .7570E-04, .7835E-06,
     3 .7181E+05,-.2174E+03, .2327E+00, .7402E-04, .8545E-06,
     3 .7181E+05,-.2173E+03, .2333E+00, .7417E-04, .7244E-06,
     3 .7180E+05,-.2173E+03, .2326E+00, .7553E-04, .8477E-06,
     3 .7179E+05,-.2173E+03, .2328E+00, .7546E-04, .7572E-06,
     3 .7178E+05,-.2172E+03, .2324E+00, .7559E-04, .7768E-06,
     3 .7174E+05,-.2170E+03, .2321E+00, .7744E-04, .6999E-06,
     3 .7167E+05,-.2166E+03, .2312E+00, .8530E-04, .5558E-06,
     3 .7153E+05,-.2157E+03, .2290E+00, .9672E-04, .3421E-06,
     3 .7125E+05,-.2139E+03, .2260E+00, .1186E-03,-.3004E-06,
     3 .7021E+05,-.2076E+03, .2131E+00, .1971E-03,-.2279E-05,
     3 .6945E+05,-.2029E+03, .2037E+00, .2511E-03,-.3651E-05,
     3 .6837E+05,-.1964E+03, .1800E+00, .3249E-03,-.3888E-05,
     3 .6681E+05,-.1870E+03, .1407E+00, .4271E-03,-.3136E-05,
     3 .6459E+05,-.1739E+03, .8587E-01, .4128E-03,-.9457E-07,
     3 .6151E+05,-.1569E+03, .3713E-01, .4919E-03, .1020E-06,
     3 .5737E+05,-.1358E+03,-.1248E-01, .7323E-03,-.7098E-06,
     3 .5210E+05,-.1101E+03,-.1068E+00, .7885E-03, .5417E-05,
     3 .4567E+05,-.8229E+02,-.1881E+00, .8320E-03, .1052E-04,
     3 .3830E+05,-.5555E+02,-.1597E+00, .5039E-03, .5656E-05,
     3 .3073E+05,-.3674E+02,-.1180E+00, .5087E-03, .3420E-05,
     3 .2383E+05,-.2042E+02,-.1343E+00,-.4800E-04, .7065E-05,
     3 .1790E+05,-.1400E+02,-.7039E-01, .8280E-04, .1786E-05,
     3 .1329E+05,-.9804E+01,-.4522E-01, .1105E-03, .1319E-05,
     3 .9701E+04,-.5980E+01,-.3954E-01, .8881E-04, .1671E-05,
     3 .6913E+04,-.3794E+01,-.1293E-01, .3508E-04,-.9205E-06,
     3 .4865E+04,-.2548E+01,-.1945E-01, .2019E-04, .1124E-05,
     3 .3391E+04,-.1770E+01,-.6844E-02, .1841E-04,-.3833E-06 /
      DATA ((CS4CH4GH(K,J), K = 1, 5), J = 1, 28)            /
     8 .1296E+05,-.2822E+01,-.2591E-02, .1616E-03, .3798E-05,
     8 .1296E+05,-.2822E+01,-.2585E-02, .1616E-03, .3797E-05,
     8 .1296E+05,-.2822E+01,-.2556E-02, .1617E-03, .3791E-05,
     8 .1296E+05,-.2822E+01,-.2585E-02, .1617E-03, .3797E-05,
     8 .1296E+05,-.2821E+01,-.2611E-02, .1618E-03, .3800E-05,
     8 .1296E+05,-.2819E+01,-.2595E-02, .1618E-03, .3796E-05,
     8 .1296E+05,-.2815E+01,-.2656E-02, .1618E-03, .3801E-05,
     8 .1296E+05,-.2808E+01,-.2697E-02, .1621E-03, .3804E-05,
     8 .1295E+05,-.2796E+01,-.2721E-02, .1618E-03, .3803E-05,
     8 .1295E+05,-.2783E+01,-.2730E-02, .1619E-03, .3806E-05,
     8 .1293E+05,-.2783E+01,-.2763E-02, .1651E-03, .3779E-05,
     8 .1292E+05,-.2781E+01,-.2802E-02, .1662E-03, .3743E-05,
     8 .1290E+05,-.2776E+01,-.3053E-02, .1685E-03, .3731E-05,
     8 .1288E+05,-.2765E+01,-.3516E-02, .1711E-03, .3718E-05,
     8 .1285E+05,-.2741E+01,-.4156E-02, .1738E-03, .3698E-05,
     8 .1282E+05,-.2701E+01,-.5078E-02, .1773E-03, .3683E-05,
     8 .1278E+05,-.2643E+01,-.6027E-02, .1805E-03, .3636E-05,
     8 .1273E+05,-.2579E+01,-.7229E-02, .1848E-03, .3617E-05,
     8 .1265E+05,-.2519E+01,-.8325E-02, .1876E-03, .3583E-05,
     8 .1258E+05,-.2474E+01,-.9120E-02, .1856E-03, .3549E-05,
     8 .1247E+05,-.2487E+01,-.9268E-02, .1833E-03, .3469E-05,
     8 .1237E+05,-.2559E+01,-.1003E-01, .1797E-03, .3539E-05,
     8 .1224E+05,-.2684E+01,-.9842E-02, .1770E-03, .3504E-05,
     8 .1208E+05,-.2799E+01,-.9249E-02, .1703E-03, .3428E-05,
     8 .1190E+05,-.2897E+01,-.8342E-02, .1695E-03, .3290E-05,
     8 .1171E+05,-.2955E+01,-.8327E-02, .1661E-03, .3156E-05,
     8 .1149E+05,-.3082E+01,-.8293E-02, .1604E-03, .2953E-05,
     8 .1121E+05,-.3367E+01,-.9956E-02, .1624E-03, .3106E-05 /
      END

      SUBROUTINE CLOUDS (TAUCS, OMCS, GCS, TAUCL, OMCL, GCL, CCLD, CLD,
     1                   CLDT, WCLW, WCLI, RADEQVW, RADEQVI, WCDW, WCDI,
     2                   IL1, IL2, ILG, LAY, LEV)
      PARAMETER (NBS = 4, NBL = 9)
C
      IMPLICIT REAL (A-H,O-Z),
     +INTEGER (I-N)
C
      REAL TAUCS(ILG,LAY,NBS), OMCS(ILG,LAY,NBS), GCS(ILG,LAY,NBS),
     1     TAUCL(ILG,LAY,NBL), OMCL(ILG,LAY,NBL), GCL(ILG,LAY,NBL)
C
      REAL CCLD(ILG,LAY), CLD(ILG,LAY), CLDT(ILG)
C
C     * INPUT ARRAYS... 

      REAL WCDW(ILG,LAY), WCDI(ILG,LAY), WCLW(ILG,LAY), WCLI(ILG,LAY)   
      REAL RADEQVW(ILG,LAY), RADEQVI(ILG,LAY)

C     * WORK ARRAY

      REAL REW1(ILG,LAY), REW2(ILG,LAY), REW3(ILG,LAY), DG(ILG,LAY),
     1     DG2(ILG,LAY), DG3(ILG,LAY), WCLWI(ILG,LAY)
C
      DATA CUT, RADW, RADI /0.001, 0.1, 0.1/
C
      COMMON / WATOP / AWS(4,NBS), BWS(4,NBS), CWS(4,NBS),
     1                 AWL(5,NBL), BWL(4,NBL), CWL(4,NBL)
      COMMON / ICEOP / AIS(2,NBS), BIS(4,NBS), CIS(4,NBS),
     1                 AIL(3,NBL), BIL(4,NBL), CIL(4,NBL)
C
C----------------------------------------------------------------------C
C     CLOUD RADIATIVE PROPERTIES FOR RADIATION.                        C
C     TAUCS, OMCS, GCS (TAUCL, OMCL, GCL): OPTICAL DEPTH, SINGLE       C
C     SCATTERING ALBEDO, ASYMMETRY FACTOR FOR SOLAR (INFRARED).        C
C                                                                      C
C     RADEQVW:     EFFECTIVE RADIU (IN MICROMETER) FOR WATER CLOUD     C
C     RADEQVI:     EFFECTIVE RADIU (IN MICROMETER) FOR ICE CLOUD       C
C     DG:          GEOMETRY LENGTH FOR ICE CLOUD                       C
C     WCDW (WCDI): LIQUID WATER (ICE) CONTENT (IN GRAM / M^3)          C
C     WCLW (WCLI): LIQUID WATER (ICE) PATH LENGTH (IN GRAM / M^2)      C
C     WCLWI:       THE SUM OF THE WATER PATH AND ICE PATH              C
C     CCLD:        CLOUD FRACTION                                      C
C     CLD:         EFFECTIVE CLOUD FRACTION (IE THE CLOUD WATER        C
C                  PATH > 0.001)                                       C
C     CLDT:        THE TOTAL CLOUD FRACTION FOR A GRID CELL            C
C                                                                      C
C     PARAMETERIZATION FOR WATER CLOUD:                                C
C     DOBBIE, ETC. 1999, JGR, 104, 2067-2079                           C
C     LINDNER, T. H. AND J. LI., 2000, J. CLIM., 13, 1797-1805.        C
C     PARAMETERIZATION FOR ICE CLOUD:                                  C
C     FU 1996, J. CLIM., 9, 2223-2337.                                 C
C     FU ET AL. 1998 J. CLIM., 11, 2223-2337.                          C
C----------------------------------------------------------------------C
C
C
      DO 100 J = 1, LAY
      DO 100 K = IL1, IL2
        IF (CCLD(K,J) .LE. CUT)                                     THEN
          REW1(K,J)  =  0.
          REW2(K,J)  =  0.
          REW3(K,J)  =  0.
          DG(K,J)    =  0.
          DG2(K,J)   =  0.
          DG3(K,J)   =  0.
          WCLWI(K,J) =  0.
          CLD (K,J)  =  0.
        ELSE
          REW1(K,J)  =  MAX (RADEQVW(K,J), RADW) 
          REW2(K,J)  =  REW1(K,J) * REW1(K,J)
          REW3(K,J)  =  REW2(K,J) * REW1(K,J)
          DG(K,J)    =  1.5396 * MAX (RADEQVI(K,J), RADI)
          DG2(K,J)   =  DG(K,J)  * DG(K,J)
          DG3(K,J)   =  DG2(K,J) * DG(K,J)
          WCLWI(K,J) =  WCLW(K,J) + WCLI(K,J)       
C
C----------------------------------------------------------------------C
C     THE EFFECTIVE CLOUD AMOUNT CLD, TO ELIMINATE THE CLOUDS WITH TOO C
C     WATER PATH                                                       C
C----------------------------------------------------------------------C
C
          IF (WCLWI(K,J) .GT. CUT)                                  THEN
            CLD(K,J) = CCLD(K,J)
          ELSE
            CLD(K,J) = 0.
          ENDIF
C
        ENDIF

  100 CONTINUE
C
      DO 290 I = 1, NBS
        DO 280 J = 1, LAY
          DO 270 K = IL1, IL2
            IF (CCLD(K,J) .LE. CUT)                                 THEN
              TAUCS(K,J,I) =  0.
              OMCS(K,J,I)  =  0.
              GCS(K,J,I)   =  0.
            ELSE
              IF (WCLWI(K,J) .GT. CUT)                              THEN
                TAUSW =  WCLW(K,J) * (AWS(1,I) + AWS(2,I) / REW1(K,J) + 
     2                   AWS(3,I) / REW2(K,J) + AWS(4,I) / REW3(K,J))
                OMSW  =  1.0 - (BWS(1,I) + BWS(2,I) * REW1(K,J) +
     1                   BWS(3,I) * REW2(K,J) + BWS(4,I) * REW3(K,J))
                GSW   =  CWS(1,I) + CWS(2,I) * REW1(K,J) +
     1                   CWS(3,I) * REW2(K,J) + CWS(4,I) * REW3(K,J)
C
                TAUSI =  WCLI(K,J) * (AIS(1,I) + AIS(2,I) / DG(K,J))
                OMSI  =  1.0 - (BIS(1,I) + BIS(2,I) * DG(K,J) +
     1                   BIS(3,I) * DG2(K,J) + BIS(4,I) * DG3(K,J))
                GSI   =  CIS(1,I) + CIS(2,I) * DG(K,J) + 
     1                   CIS(3,I) * DG2(K,J) + CIS(4,I) * DG3(K,J)
              ELSE
                TAUSW =  0.
                OMSW  =  0.
                GSW   =  0.
C
                TAUSI =  0.
                OMSI  =  0.
                GSI   =  0.
              ENDIF
C
!Yu+test
!              TAUSW = 2.*TAUSW
!              TAUSI = 2.*TAUSI
!Yu+test
              TAUCS(K,J,I)  =  TAUSW + TAUSI
              IF (TAUCS(K,J,I) .GT. 0.)                             THEN
                Y1          =  OMSW * TAUSW 
                Y2          =  OMSI * TAUSI 
                OMCS(K,J,I) = (Y1 + Y2) / TAUCS(K,J,I) 
                GCS (K,J,I) = (Y1 * GSW + Y2 * GSI) / (Y1 + Y2)
              ELSE
                OMCS(K,J,I) =  0.
                GCS (K,J,I) =  0.
              ENDIF
C
            OMCS(K,J,I) = MIN(OMCS(K,J,I),0.999)       !!!
            GCS (K,J,I) = MIN(GCS (K,J,I),0.999)       !!!
C
            ENDIF
  270     CONTINUE
  280   CONTINUE
  290 CONTINUE  

      DO 390 I = 1, NBL
        DO 380 J = 1, LAY
          DO 370 K = IL1, IL2
            IF (CCLD(K,J) .LE. CUT)                                 THEN
              TAUCL(K,J,I) =  0.
              OMCL(K,J,I)  =  0.
              GCL(K,J,I)   =  0.
            ELSE
              IF (WCLWI(K,J) .GT. CUT)                              THEN
                TAULW  =  WCLW(K,J) * (AWL(1,I) + AWL(2,I) * REW1(K,J) + 
     1                    AWL(3,I) / REW1(K,J) + AWL(4,I) / REW2(K,J) + 
     2                    AWL(5,I) / REW3(K,J))
                OMLW   =  1.0 - (BWL(1,I) + BWL(2,I) / REW1(K,J) + 
     1                    BWL(3,I) * REW1(K,J) + BWL(4,I) * REW2(K,J))
                GLW    =  CWL(1,I) + CWL(2,I) / REW1(K,J) + 
     1                    CWL(3,I) * REW1(K,J) + CWL(4,I) * REW2(K,J)
C
C----------------------------------------------------------------------C
C     SINCE IN FU ETC. THE PARAM. IS FOR ABSORPTANCE, SO NEED A FACTOR C
C     WCLI(K,J) / TAULI FOR SINGLE SCATTERING ALBEDO                   C
C----------------------------------------------------------------------C
C
                TAULI  =  WCLI(K,J) * (AIL(1,I) + AIL(2,I) / DG(K,J) +
     1                   AIL(3,I) / DG2(K,J))
                IF (TAULI .GT. 0.0)                                 THEN
                  OMLI =  1.0 - (BIL(1,I) / DG(K,J) + BIL(2,I) +
     1                    BIL(3,I) * DG(K,J) + BIL(4,I) * DG2(K,J)) *
     2                    WCLI(K,J) / TAULI
                ELSE
                  OMLI =  0.0
                ENDIF
C
                GLI    =  CIL(1,I) + CIL(2,I) * DG(K,J) + 
     1                    CIL(3,I) * DG2(K,J) + CIL(4,I) * DG3(K,J)
              ELSE
                TAULW  =  0.
                OMLW   =  0.
                GLW    =  0.
C
                TAULI  =  0.
                OMLI   =  0.
                GLI    =  0.
              ENDIF
C
              TAUCL(K,J,I)   =  TAULW + TAULI
              IF (TAUCL(K,J,I) .GT. 0.)                             THEN
                Y1           =  OMLW * TAULW 
                Y2           =  OMLI * TAULI
                OMCL(K,J,I)  = (Y1 + Y2) / TAUCL(K,J,I)
                GCL (K,J,I)  = (GLW * Y1 + GLI * Y2) / (Y1 + Y2)
              ELSE
                OMCL(K,J,I)  =  0.
                GCL (K,J,I)  =  0.
              ENDIF
            ENDIF
  370     CONTINUE
  380   CONTINUE
  390 CONTINUE 
C
C----------------------------------------------------------------------C
C     CALCULATE THE TOTAL DIAGNOSTIC OVERLAPPED CLOUD AMOUNT FOR THE   C
C     GRID CELL BASED ON THE RANDOM / MAXIMUM OVERLAP ASSUMPTION       C
C----------------------------------------------------------------------C
C
      DO 460 K = IL1, IL2
        CLDT(K) = CLD(K,1)
  460 CONTINUE
C
      DO 490 J = 2, LAY
        JM1 = J - 1
        DO 480 K = IL1, IL2
          IF (CLD(K,JM1) .GT. CUT)                                  THEN
            CLDT(K) = MAX(CLDT(K), CLD(K,J))
          ELSE
            CLDT(K) = 1. - (1. - CLDT(K)) * (1. - CLD(K,J))
          ENDIF
  480   CONTINUE
  490 CONTINUE
C
      RETURN
      END

      BLOCK DATA DATCLDOP
      PARAMETER (NBS = 4, NBL = 9)
      IMPLICIT REAL (A-H,O-Z),
     +INTEGER (I-N)

      COMMON / WATOP / AWS(4,NBS), BWS(4,NBS), CWS(4,NBS),
     1                 AWL(5,NBL), BWL(4,NBL), CWL(4,NBL)
      COMMON / ICEOP / AIS(2,NBS), BIS(4,NBS), CIS(4,NBS),
     1                 AIL(3,NBL), BIL(4,NBL), CIL(4,NBL)
C
C     --------------------------------------------------------
C     NEW WATER PROPERTIES FOR SW (DOBBIE, LI, AND CHYLEK 1999, JGR)
C     --------------------------------------------------------
C
        DATA ((AWS(I,J), I = 1, 4), J = 1, NBS)           /
     1  4.554E-04,   1.500E+00,   7.190E-01,  -9.419E-01,
     2  3.859E-04,   1.508E+00,   9.512E-01,  -1.053E+00,
     3 -3.946E-05,   1.538E+00,   1.035E+00,   2.638E-01,
     4  2.936E-04,   1.541E+00,   1.698E-00,   1.521E+00  /
        DATA ((BWS(I,J), I = 1, 4), J = 1, NBS)           /
     1  6.481E-08,   1.553E-07,  -7.755E-10,   7.616E-12,
     2  1.072E-05,   1.345E-05,  -1.799E-08,  -3.146E-11,
     3  4.078E-04,   2.169E-03,  -2.177E-05,   1.506E-07,
     4  2.013E-01,   1.109E-02,  -2.897E-04,   3.055E-06  /
        DATA ((CWS(I,J), I = 1, 4), J = 1, NBS)           /
     1  8.069E-01,   6.188E-03,  -2.065E-04,   2.352E-06,
     2  7.685E-01,   9.337E-03,  -3.101E-04,   3.527E-06,
     3  7.471E-01,   9.440E-03,  -2.616E-04,   2.614E-06,
     4  7.956E-01,   8.138E-03,  -1.861E-04,   1.611E-06  /
C
C     --------------------------------------------------------
C     NEW WATER PROPERTIES FOR LW (LINDNER AND LI, 2000 JCL)
C     --------------------------------------------------------
C
        DATA ((AWL(I,J), I = 1, 5), J = 1, NBL)                      /
     1   -.21671E-01, .79578E-03, .14899E+01, .62606E+01,-.12705E+02,
     2   -.14126E+00, .28208E-02, .35125E+01,-.34541E+01,-.22679E+01,
     3   -.18829E+00, .34065E-02, .46731E+01,-.11664E+02, .87105E+01,
     4   -.16383E+00, .26574E-02, .48670E+01,-.16442E+02, .16128E+02,
     5   -.20294E-01, .85110E-04, .28650E+01,-.11202E+02, .12047E+02,
     6    .28752E-01,-.37315E-03, .14591E+01,-.48788E+01, .49725E+01,
     7   -.40386E-01, .80822E-03, .25318E+01,-.64641E+01, .55609E+01,
     8   -.48716E-01, .81275E-03, .30390E+01,-.97845E+01, .95101E+01,
     9    .64794E-01,-.98530E-03, .12797E+01,-.55272E+01, .62599E+01 /
        DATA ((BWL(I,J), I = 1, 4), J = 1, NBL)          /
     1    .36899E-02,-.54184E-03, .14561E-01,-.18451E-03,
     2    .62141E-02, .61190E-01, .21127E-01,-.29731E-03,
     3    .87326E-01, .29908E+00, .22928E-01,-.35569E-03,
     4   -.37551E-01, .70237E+00, .26945E-01,-.37999E-03,
     5    .51671E-01, .10199E+01, .18296E-01,-.21209E-03,
     6    .52184E+00, .72352E+00,-.48090E-02, .10414E-03,
     7    .57688E+00, .63008E+00,-.56325E-02, .87852E-04,
     8    .50346E+00, .79407E+00,-.13179E-02, .25467E-04,
     9    .67792E+00, .68259E+00,-.12136E-01, .20941E-03 /
        DATA ((CWL(I,J), I = 1, 4), J = 1, NBL)          /
     1    .73147E+00, .11761E+00, .86402E-02,-.10761E-03,
     2    .81284E+00,-.60287E-01, .45367E-02,-.33372E-04,
     3    .92468E+00,-.39653E+00, .30494E-03, .20980E-04,
     4    .10006E+01,-.71422E+00,-.46784E-02, .10114E-03,
     5    .10635E+01,-.10097E+01,-.58726E-02, .97485E-04,
     6    .10762E+01,-.12482E+01,-.40343E-02, .54330E-04,
     7    .97445E+00,-.13875E+01, .79204E-03,-.27995E-04,
     8    .79053E+00,-.13566E+01, .10452E-01,-.18111E-03,
     9    .35512E+00,-.80671E+00, .30384E-01,-.47204E-03 /
C
C----------------------------------------------------------------------C
C    ICE FU 1997 JCL, FU ET AL. 1998 JCL                               C
C----------------------------------------------------------------------C
C
        DATA ((AIS(I,J), I = 1, 2), J = 1, NBS) /
     1 -0.24276E-04, 2.51884E+00,
     2 -0.48500E-04, 2.52275E+00,
     3 -0.98503E-05, 2.52048E+00,
     4  0.24435E-03, 2.49116E+00                /
C
        DATA ((BIS(I,J), I = 1, 4), J = 1, NBS)            /
     1  0.13031E-06, 0.94102E-07,-0.75971E-10, 0.33977E-12,
     2 -0.77603E-06, 0.73420E-05, 0.11514E-09,-0.90818E-12,
     3  0.10007E-02, 0.10992E-02,-0.45043E-05, 0.12637E-07,
     4  0.21201E+00, 0.25713E-02,-0.19228E-04, 0.62183E-07 /
C
        DATA ((CIS(I,J), I = 1, 4), J = 1, NBS)            /
     1  0.74821E+00, 0.92318E-03,-0.72862E-06,-0.95642E-08,
     2  0.75227E+00, 0.10653E-02,-0.24930E-05,-0.29114E-08,
     3  0.75553E+00, 0.17297E-02,-0.87585E-05, 0.19201E-07,
     4  0.84323E+00, 0.20925E-02,-0.18302E-04, 0.60381E-07 /
C
        DATA ((AIL(I,J), I = 1, 3), J = 1, NBL)    /
     1 -.8839455E-03,  .2662598E+01,  .2196338E+01,
     2 -.2066995E-02,  .2787904E+01,  .1397838E+01,
     3 -.3085730E-02,  .2906257E+01, -.1911363E+01,
     4 -.6968920E-02,  .3284275E+01, -.6973825E+01,
     5 -.8372696E-02,  .3455018E+01, -.1516692E+02,
     6 -.1691632E-02,  .2765756E+01, -.8331033E+01,
     7 -.7098616E-02,  .3343404E+01, -.8144649E+01,
     8 -.1041746E-01,  .3824226E+01, -.2255945E+02,
     9  .5689700E-02,  .2285636E+01, -.1430752E+02 /
C
        DATA ((BIL(I,J), I = 1, 4), J = 1, NBL)                   /
     1  .5723611E+00,  .1627863E-01, -.1684272E-03,  .6061332E-06,
     2  .4402328E+00,  .1736939E-01, -.1656608E-03,  .5709622E-06,
     3  .8802908E+00,  .1249744E-01, -.1550609E-03,  .6105065E-06,
     4  .6351165E+00,  .1781519E-01, -.1979682E-03,  .6000892E-06,
     5  .5409536E-00,  .1949649E-01, -.2050908E-03,  .7364680E-06,
     6  .1195515E+01,  .3350616E-02, -.5266996E-04,  .2233377E-06,
     7  .1186334E+01,  .6213290E-02, -.1044277E-03,  .2233377E-06,
     8  .2279562E+00,  .2017007E-01, -.1756872E-03,  .5703918E-06,
     9  .7718967E+00,  .2120626E-01, -.2587649E-03,  .9878070E-06 /
C
        DATA ((CIL(I,J), I = 1, 4), J = 1, NBL)                   /
     1  .7975757E+00,  .3843973E-02, -.3540463E-04,  .1179791E-06,
     2  .7947997E+00,  .3190423E-02, -.2386042E-04,  .6691811E-07,
     3  .8737279E+00,  .2465886E-02, -.2468764E-04,  .8686448E-07,
     4  .8577221E+00,  .2321034E-02, -.1897764E-04,  .8641223E-07,
     5  .8906280E-00,  .1903269E-02, -.1733552E-04,  .5855071E-07,
     6  .8663385E-00,  .2797934E-02, -.3187011E-04,  .1217209E-06,
     7  .7644037E+00,  .4427001E-02, -.4494615E-04,  .1217209E-06,
     8  .7200100E+00,  .3220301E-02, -.2195542E-04,  .6604318E-07,
     9  .5355918E+00,  .1127081E-01, -.1234705E-03,  .4567953E-06 /
C

      END

!      END MODULE APM_RADF_MOD
#endif
