#ifdef APM
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: apm_grow_mod
!
! !DESCRIPTION: Module APM\_GROW\_MOD contains variables and routines for 
!  computing size-resolved particle growth. 
!\\
!\\
! !INTERFACE:
!
      MODULE APM_GROW_MOD

!
! !USES:
!
      IMPLICIT NONE
      PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!  
      PUBLIC :: APM_GROW
      PUBLIC :: APM_MOVEBIN

! !REVISION HISTORY: 
!  23 Aug 2008 - F. Yu       - Initial version  
!  08 Nov 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
      CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: apm_grow
!
! !DESCRIPTION: Subroutine APM\_GROW calculates H2SO4 condensational growth 
!  of aerosol particles
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE APM_GROW( ICOND,    NSO4,   TK,
     &                     PRESS,    CCOND,  PCOND,
     &                     TCSOTHER, DTNG,   XN,
     &                     XVA,      RWETCM, TOTCONDOTHER,
     &                     XMCOND )
!
! !USES:
!
      USE APM_INIT_MOD, ONLY: DENSULF, RDRY
      USE APM_INIT_MOD, ONLY: ONEPI,   BK,  AVG, RGAS
!
! !INPUT/OUTPUT PARAMETERS: 
! 
      INTEGER :: NSO4            ! Dimension variable
      REAL*8  :: TK              ! Temperature (K)       
      REAL*8  :: PRESS           ! Pressure (pa)
      REAL*8  :: PCOND           ! Condensable vapor production rate (#/cm3s) 
      REAL*8  :: DTNG            ! Time detp for growth (s)
      REAL*8  :: XN(NSO4)        ! Number conc of each bin (#/cm3)
      REAL*8  :: XVA(NSO4)       ! Total acid volume of each bin (cm3/cm3)
!
! !INPUT/OUTPUT PARAMETERS: 
! 
      INTEGER :: ICOND
      REAL*8  :: CCOND           ! Condensable vapor concentration (#/cm3)
      REAL*8  :: TCSOTHER
      REAL*8  :: RWETCM(NSO4)
      REAL*8  :: TOTCONDOTHER
      REAL*8  :: XMCOND
!
! !REMARKS:
! Output:
!   updated XMA, CCOND 
! 
! !REVISION HISTORY: 
!  17 Mar 2010 - F. Yu       - Initial version  
!  08 Nov 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER :: N, NMAX
      REAL*8  :: XMA(NSO4)
      REAL*8  :: YF(NSO4),YGR(NSO4),AKELV(NSO4)
      REAL*8  :: CSCOND,YFSUM1
      REAL*8  :: V1COND,VCOND,TEMP0,WTGAS
      REAL*8  :: WTAIR,DIAMAIR,RHOA,CDIFUS,CDIFUS2,DIFUSC,FREEPD
      REAL*8  :: YY,YKN,FCORR,TEMP1,YFSUM,YEVAP
      REAL*8  :: YEVAPV,CCOND1,CCONDA,TOTEVAP,YGRSUM,YDV,TOTCOND
      REAL*8  :: XRCM
      REAL*8  :: DENCOND,AKELV0,SURFT,YTEMP1
      REAL*8  :: DVA(NSO4)

      !=================================================================
      ! APM_GROW begins here!
      !=================================================================
!
      DENCOND = DENSULF
      IF(ICOND.EQ.1) THEN       !H2SO4
         CSCOND = 1.E4          ! no acid evap now
         AKELV0 = CCOND/CSCOND
         IF(AKELV0.LT.1.2) RETURN 
         AKELV =1.
        
!APM2++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      ELSEIF(ICOND.EQ.2) THEN  ! LV-SOG
       SURFT = 30. !erg/cm2
!Yu20111015       CSCOND = 3.3E6  !
       CSCOND = 1.0E6  !
!  RORG0 : organic activation radius (in cm)
!       RORG0 =ALOG(CCOND/CSCOND)*DENCOND*RGAS*TK/(2.*SURFT*XMORG)
!       AKELV0 = exp(2.*SURFT*XMCOND/(DENCOND*RGAS*TK*RORG0))
        AKELV0 = CCOND/CSCOND
       IF(AKELV0.LE.1.2) RETURN
       DO N = 1, NSO4
        XRCM = RWETCM(N)   ! RWETCM in cm
        AKELV(N) = exp(2.*SURFT*XMCOND/(DENCOND*RGAS*TK*XRCM))
       ENDDO
!APM2++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

 
      ELSE
         WRITE(6,*)"NEED to check ICOND value"
         STOP
      ENDIF
      
      V1COND = XMCOND/(AVG*DENCOND) ! volume of one molecule (cm3)
      
      VCOND = SQRT(8.*RGAS*TK/(ONEPI*XMCOND)) ! cm/s
      TEMP0 = VCOND * XMCOND/AVG *ONEPI

! Cal. gas diffusion coef. and mean free path
!
      WTGAS = XMCOND
      WTAIR = 28.966
      DIAMAIR   = 4.5E-08
      RHOA = WTAIR*PRESS*10./(RGAS*TK) ! AIR DENSITY (G CM-3)
      CDIFUS =3.*SQRT(0.5*RGAS*WTAIR/ONEPI)/(AVG*8.*DIAMAIR**2)
      CDIFUS2=CDIFUS * SQRT((WTGAS+WTAIR)/WTGAS)
      DIFUSC=CDIFUS2 * SQRT(TK)/RHOA
      FREEPD = 3.*DIFUSC/VCOND
!
      YF = 0.
      YFSUM = 0.
      YY = PCOND
      DO N = 1, NSO4
         XRCM = RWETCM(N)       ! RWETCM in cm
         YKN = FREEPD/XRCM
         FCORR = YKN/(0.75+YKN)
         TEMP1 = TEMP0*FCORR
         IF(AKELV(N).LT.AKELV0) THEN
             TEMP1 = TEMP1 *(1.-AKELV(N)/AKELV0)
         ELSE
            TEMP1 = 0.          ! no evapor for now
         ENDIF
         YF(N) = TEMP1 * XN(N) *XRCM*XRCM/(DENCOND*V1COND) ! s-1
         YFSUM = YFSUM + YF(N)
!         YEVAP = CSCOND*AKELV*YF(N)    ! in #/cm3s
!         YEVAPV = YEVAP*DTNG*V1COND      ! in cm3/cm3
!         IF(YEVAPV.GT. XVA(N)) THEN
!             YEVAP = XVA(N)/(DTNG*V1COND)
!         ENDIF
!         YY = YY + YEVAP
      ENDDO
      IF(YFSUM.EQ.0.) THEN
         WRITE(6,100)AKELV0,(AKELV(N),N=1,NSO4,5)
         WRITE(6,100)FCORR,(YF(N),N=1,NSO4,5)
         STOP
      ENDIF
! Consider scavenging of H2SO4 vapor by aerosols other than sulfate
      YFSUM1 = YFSUM + TCSOTHER
      CCOND1 = YY/YFSUM1 + (CCOND - YY/YFSUM1)*exp(-YFSUM1*DTNG)
       
      CCONDA = 0.5 * (CCOND1 + CCOND)
      
      TOTEVAP = 0.
!       YGRSUM = 0.
!       DO N = 1, NSO4
!          YGR(N) = YF(N)*(CCONDA - CSCOND*AKELV)
!          IF(YGR(N).LT.0.) THEN
!            YDV = -YGR(N)*DTNG*V1COND    !cm3
!            IF(YDV.GT.XVA(N)) THEN
!               TOTEVAP = TOTEVAP + XVA(N)
!               XVA(N) = 1.E-50
!            ELSE
!               XVA(N) = XVA(N) - YDV
!               TOTEVAP = TOTEVAP + YDV
!            ENDIF
!          ELSE
!            YGRSUM = YGRSUM + YGR(N)
!          ENDIF
!       ENDDO
!
      TOTCOND=PCOND*DTNG+(CCOND-CCOND1)+TOTEVAP/V1COND ! #/cm3
!  H2SO4 vapor condensing on particles other than sulfate
      TOTCONDOTHER = TOTCOND * (1.-YFSUM/YFSUM1)
!  H2SO4 vapor condensing on sulfate
      TOTCOND = TOTCOND * YFSUM/YFSUM1

      YGR = YF
      YGRSUM = YFSUM
       
      DO N = 1, NSO4
         DVA(N) = TOTCOND*V1COND*YGR(N)/YGRSUM
         XVA(N)=XVA(N)+DVA(N)
      ENDDO
       
       
 100  FORMAT(20(1PE9.2))
!
! Update gas concentration
      CCOND = CCOND1

      END SUBROUTINE APM_GROW
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: apm_movebin
!
! !DESCRIPTION: Subroutine APM\_MOVEBIN moves particles across bins after 
!  growth.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE APM_MOVEBIN( NSO4, XN, XVA )
!
! !USES:
!
      USE APM_INIT_MOD, ONLY: VDRY
!
! !INPUT PARAMETERS: 
!
      INTEGER :: NSO4        ! Dimension for arrays
!
! !INPUT/OUTPUT PARAMETERS: 
! 
      REAL*8  :: XN(NSO4)    ! Number conc of each bin (#/cm3)
      REAL*8  :: XVA(NSO4)   ! Total acid mass of each bin (cm3/cm3)   
!
! !REVISION HISTORY: 
!  17 Mar 2010 - F. Yu       - Initial version  
!  08 Nov 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
      ! Number of core component associated with particles, set to 1 now,
      ! may increase later
      INTEGER, PARAMETER :: NC = 1
!
! !LOCAL VARIABLES:
!
      INTEGER :: N, NMAX, IC, J, J0, J1, JS

      REAL*8  :: YN(NSO4),YVA(NSO4,NC),YC(NSO4,NC)
      REAL*8  :: XVT, VDA, YFR, YFV

      NMAX = NSO4
!
!  Move the partilces across bins due to the condensation/evap
!
      DO N=1,NMAX
         XN(N) = MAX(1.D-30,XN(N))
         XVA(N) = MAX(1.D-30,XVA(N))
         YVA(N,1) = XVA(N)      ! to be modified if NC>1
      ENDDO

      DO n=1,NMAX
         YN(n)=1.E-20
         DO IC=1,NC
            YC(n,IC)=1./float(NC)*YN(n)*VDRY(n)*1.E6 ! in cm3/cm3
         ENDDO

!         IF(YN(n).GT.1.E6.or.YN(n).LE.0.) THEN
!           WRITE(6,*)"4 XN=",n,YN(n),XN(n),YC(n,1)
!         ENDIF

      ENDDO

      DO n=1,NMAX
         XVT = 0.
         DO IC=1,NC
            XVT = XVT + YVA(n,IC) ! cm3/cm3
         ENDDO
         VDA = XVT/XN(n)*1.E-6  ! in m3

!         YRATIO(n) = VDA/VDRY(n)

! if particles become smaller than first bin or larger than last bin,
! scale and put in the first or last bin in a fashion that conserve mass
         IF(VDA.LT.VDRY(1)) THEN
            YN(1) = YN(1) + XN(n)*VDA/VDRY(1)  
            DO IC=1,NC          ! move compositions
               YC(1,IC) = YC(1,IC)+YVA(n,IC)
            ENDDO

!            IF(YN(n).GT.1.E6.or.YN(n).LE.0.) THEN
!              WRITE(6,*)"41 XN=",n,YN(n),XN(n),YC(n,1),YVA(n,1),XVT,VDA
!            ENDIF

         ELSEIF(VDA.GE.VDRY(NMAX)) THEN
            YN(NMAX) = YN(NMAX) + XN(n)*VDA/VDRY(NMAX)
            DO IC=1,NC          ! move compositions
               YC(NMAX,IC) = YC(NMAX,IC)+YVA(n,IC)
            ENDDO
!            IF(YN(n).GT.1.E6.or.YN(n).LE.0.) THEN
!              WRITE(6,*)"42 XN=",n,YN(n),XN(n),YC(n,1),
!     &            XVA(n),YVA(n,1),XVT,VDA
!            ENDIF
         ELSE
            IF(VDA.GE.VDRY(n)) THEN
               J0=n
               J1 = NMAX-1
               JS = 1
            ELSE
               J0=n-1
               J1=1
               JS=-1
            ENDIF
            DO J=J0,J1,JS
               IF(VDA.GE.VDRY(J).AND.VDA.LT.VDRY(J+1)) THEN
                  YFR = (VDRY(J+1)-VDA)/(VDRY(J+1)-VDRY(J))
                  YFV = YFR*VDRY(J)/VDA
                  YN(J) = YN(J) + YFR*XN(n) !move #
                  YN(J+1)=YN(J+1)+(1.-YFR)*XN(n)
                  DO IC=1,NC    ! move compositions
                     YC(J,IC) = YC(J,IC)+YFV*YVA(n,IC)
                     YC(J+1,IC)=YC(J+1,IC)
     &                    +(1.-YFV)*YVA(n,IC)
                  ENDDO
                  GOTO 160
               ENDIF
            ENDDO
         ENDIF
 160     CONTINUE
      ENDDO

!  Update the bin values

      DO n=1,NMAX
         
         IF(YN(n).GT.1.E6.or.YN(n).LE.0.) THEN
            WRITE(6,*)"5 XN=",n,YN(n),XN(n),YC(n,1)
         ENDIF
         
         !GanLuo20190501XN(n) = YN(n)
         XN(n) = MAX(1.D-30,YN(n))
         DO IC=1,NC
            !GanLuo20190501YVA(n,IC) =YC(n,IC)
            YVA(n,IC) =MAX(1.D-30,YC(n,IC))
         ENDDO
         XVA(N) = YVA(n,1)      ! to be modified if NC>1

         XN(N) = MAX(0.d0,XN(N))
         XVA(N) = MAX(0.d0,XVA(N))
      ENDDO

      END SUBROUTINE APM_MOVEBIN
!EOC
      END MODULE APM_GROW_MOD
#endif
