#ifdef APM
!************************************************************************
! This is the module to calculate aerosol optical properties based on
! APM simulated particle size distribution, composition, and mixing state.
! Designed and written by
! Fangqun Yu
! SUNY-Albany
! 09/2010-03/2011
!************************************************************************
      MODULE APM_OPTI_MOD
        
      implicit none
       
      ! Make everything PRIVATE ... 
      PRIVATE
      
      ! ... except these variables ...
!      PUBLIC ::
 
      ! ... and these routines
      PUBLIC :: READOPTABLE, READOPTABLE_LW
      PUBLIC :: APM_OPT, APM_OPT_LW
      PUBLIC :: OPTABLE1,OPTABLE2,OPTABLE3   !for comp only

      !=================================================================
      ! MODULE VARIABLES
      !=================================================================
      ! Parameters
!      integer,parameter :: MWL ! Num of  wavelength across solar spectrum
!      integer,parameter :: MDC#  ! Num of Dcore
!      integer,parameter :: MDS  ! Num of Dshell
!      integer,parameter :: MSR  ! Num of real part of shell refindx

      INTEGER, PARAMETER   :: MWL=16,MSR=6,MSI=17  !wavelength,real refindx, imag refindx
!YuNB      INTEGER, PARAMETER   :: MDC1=1, MDC2=21, MDC3=31 !Dcores for type 1, 2, 3
!YuNB      INTEGER, PARAMETER   :: MDS1=31, MDS2=32, MDS3=32 !Dshell for type 1, 2, 3
      INTEGER, PARAMETER   :: MDC1=1, MDC2=21, MDC3=61 !Dcores for type 1, 2, 3
      INTEGER, PARAMETER   :: MDS1=91, MDS2=32, MDS3=17 !Dshell for type 1, 2, 3


      real*8 :: RESR, RESI
      real*8 :: RESC1,RESC2,RESC3,RESS1,RESS2,RESS3

      ! Arrays
      real*8 :: WAVL(MWL),RSR(MSR),RSI(MSI)
      real*8 :: DC1(MDC1),DC2(MDC2),DC3(MDC3)
      real*8 :: DS1(MDS1),DS2(MDS2),DS3(MDS3)

      real*8 :: OPT1EXT(MWL,MDS1,MSR,MSI),OPT1W(MWL,MDS1,MSR,MSI)
      real*8 :: OPT1G(MWL,MDS1,MSR,MSI)

      real*8 :: OPT2EXT(MWL,MDC2,MDS2,MSR,MSI)
      real*8 :: OPT2W(MWL,MDC2,MDS2,MSR,MSI)
      real*8 :: OPT2G(MWL,MDC2,MDS2,MSR,MSI)

      real*8 :: OPT3EXT(MWL,MDC3,MDS3,MSR,MSI)
      real*8 :: OPT3W(MWL,MDC3,MDS3,MSR,MSI)
      real*8 :: OPT3G(MWL,MDC3,MDS3,MSR,MSI)

! Longwave
      INTEGER, PARAMETER   :: MWLL=41
      INTEGER, PARAMETER   :: MDCLW = 31
      real*8 :: WAVLL(MWLL),DC3LW(MDCLW)
      real*8 :: OPT3EXT_LW(MWLL,MDCLW),OPT3W_LW(MWLL,MDCLW)   !2D for dust
      real*8 :: OPT3G_LW(MWLL,MDCLW)

      contains

!************************************************************************
      subroutine APM_OPT(ITEST,ITYPE,MBIN,XDC,XN,ZFV,NWL,WL,XBEXT,XW,XG)

      INTEGER   :: ITEST, ITYPE,MBIN,NWL
      REAL*8    :: XDC(1:MBIN),XN(1:MBIN)
      REAL*8    :: ZFV(6)
      REAL*8    :: WL(NWL),XBEXT(NWL),XW(NWL),XG(NWL)

      REAL*8    :: XDCORE(MBIN),XDWET(MBIN)
      REAL*8    :: XWL,BSCATSUM,BEXTSUM,GSUM,DCORE,DWET,DSHELL
      REAL*8    :: QEXT,YW,YG,BEXT,BSCAT

! REAL PART OF REF INDEX 
! (Table 1 in Aouizerats et al., 2010, 550 nm, consider wavelength dependence later`)
      REAL*8, PARAMETER   :: ASO4 = 1.52, BSO4=5.E-4
      REAL*8, PARAMETER   :: ANH4 = 1.52, BNH4=5.E-4
      REAL*8, PARAMETER   :: ANO3 = 1.53, BNO3=5.E-3
      REAL*8, PARAMETER   :: ASOA = 1.45, BSOA=1.E-3
      REAL*8, PARAMETER   :: AH2O = 1.33, BH2O=1.8E-8
      REAL*8, PARAMETER   :: APOC = 1.45, BPOC=1.E-3
! Table 1.12 in Krekov 1993
      REAL*8, PARAMETER   :: ASALT = 1.45, BSALT=1.5E-4

      REAL*8    :: AS,BS,ZFVSUM

      REAL*8    :: RDRY(40)
      REAL*8, PARAMETER   :: ONEPI = 3.1415926d0

      INTEGER   :: IWL,I,J

      IF(ITYPE.EQ.1) THEN  !SP
       DO I=1,MBIN
         XDCORE(I)=XDC(I)  ! um
!bug20120211         XDWET(I)= XDCORE(I) * (1./ZFV(1))**(1./3.)  !um
         XDWET(I)= XDCORE(I) * (1./(ZFV(1)+ZFV(5)))**(1./3.)  !um
         IF(ITEST.EQ.1)WRITE(1001,201)ITYPE,I,XDCORE(I), XDWET(I)
         XDCORE(I)=0.
       ENDDO
       AS=ASO4*ZFV(1)+ANH4*ZFV(3)+ANO3*ZFV(4)+ASOA*ZFV(5)+AH2O*ZFV(6)
       BS=BSO4*ZFV(1)+BNH4*ZFV(3)+BNO3*ZFV(4)+BSOA*ZFV(5)+BH2O*ZFV(6)
      ELSEIF(ITYPE.EQ.2) THEN   !sea salt
       DO I=1,MBIN
         XDCORE(I)=XDC(I)  ! um
         XDWET(I)= XDCORE(I) * (1./ZFV(1))**(1./3.)  !um
         XDCORE(I)=0.
       ENDDO
       AS=ASALT*ZFV(1)+ASO4*ZFV(2)+ANH4*ZFV(3)+ANO3*ZFV(4)
     &    +ASOA*ZFV(5)+AH2O*ZFV(6)
       BS=BSALT*ZFV(1)+BSO4*ZFV(2)+BNH4*ZFV(3)+BNO3*ZFV(4)
     &    +BSOA*ZFV(5)+BH2O*ZFV(6)
      ELSEIF(ITYPE.EQ.3) THEN    !dust
       DO I=1,MBIN
         XDCORE(I)=XDC(I)  ! um
         XDWET(I)= XDCORE(I) * (1./ZFV(1))**(1./3.)  !um
       ENDDO
       ZFVSUM = SUM(ZFV(2:6))
       IF(ZFVSUM.LT.1.d-4) THEN
         AS = ASO4
         BS = BSO4
       ELSE
         AS=ASO4*ZFV(2)+ANH4*ZFV(3)+ANO3*ZFV(4)+ASOA*ZFV(5)+AH2O*ZFV(6)
         BS=BSO4*ZFV(2)+BNH4*ZFV(3)+BNO3*ZFV(4)+BSOA*ZFV(5)+BH2O*ZFV(6)
         AS=AS/ZFVSUM
         BS=BS/ZFVSUM
       ENDIF
      ELSEIF(ITYPE.EQ.4) THEN  !BC
       DO I=1,MBIN
         XDCORE(I)=XDC(I)  ! um
         XDWET(I)= XDCORE(I) * (1./ZFV(1))**(1./3.)  !um
       ENDDO
       ZFVSUM = SUM(ZFV(2:6))
       IF(ZFVSUM.LT.1.d-4) THEN
         AS = ASO4
         BS = BSO4
       ELSE
         AS=ASO4*ZFV(2)+ANH4*ZFV(3)+ANO3*ZFV(4)+ASOA*ZFV(5)+AH2O*ZFV(6)
         BS=BSO4*ZFV(2)+BNH4*ZFV(3)+BNO3*ZFV(4)+BSOA*ZFV(5)+BH2O*ZFV(6)
         AS=AS/ZFVSUM
         BS=BS/ZFVSUM
       ENDIF
      ELSEIF(ITYPE.EQ.5) THEN   !POC
       DO I=1,MBIN
         XDCORE(I)=XDC(I)  ! um
         XDWET(I)= XDCORE(I) * (1./ZFV(1))**(1./3.)  !um
         XDCORE(I)=0.
       ENDDO
       AS=APOC*ZFV(1)+ASO4*ZFV(2)+ANH4*ZFV(3)+ANO3*ZFV(4)
     &    +ASOA*ZFV(5)+AH2O*ZFV(6)
       BS=BPOC*ZFV(1)+BSO4*ZFV(2)+BNH4*ZFV(3)+BNO3*ZFV(4)
     &    +BSOA*ZFV(5)+BH2O*ZFV(6)
      ELSE
       WRITE(6,*)"STOP: NEED to check ITYPE",ITYPE
       STOP
      ENDIF

      XBEXT = 0.
      XW = 1.
      XG = 0.5
      DO IWL = 1, NWL
        XWL = WL(IWL)
        BSCATSUM = 1.d-20
        BEXTSUM = 1.d-20
        GSUM = 0.5d-20
        DO I = 1, MBIN
         IF(XN(I).GT.1.d-5) THEN
           DCORE = XDCORE(I)    !um
           DWET =  XDWET(I)    !um
           DSHELL = DWET - DCORE    !um

!           IF(ITYPE.EQ.1.or.ITYPE.EQ.2.or.ITYPE.EQ.5) THEN           !SP,salt,POC
!            IF(DSHELL.GT.20.) THEN
!             WRITE(86,201) ITYPE,I, XDC(I), DWET, XN(I),(ZFV(J),J=1,6)
!            ENDIF
!           ELSE
!            IF(DSHELL.GT.1.) THEN
!             WRITE(86,201) ITYPE,I, XDC(I), DWET, XN(I),(ZFV(J),J=1,6)
!            ENDIF
!           ENDIF

           IF(ITYPE.EQ.1.or.ITYPE.EQ.2.or.ITYPE.EQ.5) THEN           !SP,salt,POC
!             CALL OPTABLE1(XWL,DCORE,DSHELL,AS,QEXT,YW,YG)
             CALL OPTABLE1(IWL,DSHELL,AS,BS,QEXT,YW,YG)
           ELSEIF(ITYPE.EQ.4) THEN                                   !BC
             CALL OPTABLE2(IWL,DCORE,DSHELL,AS,BS,QEXT,YW,YG)
           ELSEIF(ITYPE.EQ.3) THEN                                   !dust
             CALL OPTABLE3(IWL,DCORE,DSHELL,AS,BS,QEXT,YW,YG)
           ENDIF

!           IF(YW.GT.1.) THEN
!              WRITE(1001,*)IWL,DCORE,DSHELL,AS,BS,QEXT,YW,YG
!           ENDIF
           BEXT = QEXT*XN(I)*DWET*DWET*1.E-8*ONEPI/4.   !cm-1
           BSCAT = BEXT*YW
           BEXTSUM = BEXTSUM +  BEXT
           BSCATSUM = BSCATSUM + BSCAT
           GSUM = GSUM + BSCAT*YG
         ENDIF
        ENDDO
        XBEXT(IWL) = BEXTSUM             ! cm-1
        XW(IWL) = BSCATSUM/BEXTSUM
        XG(IWL) = GSUM/BSCATSUM

      ENDDO
200   FORMAT(10(1PE10.3))
201   FORMAT(I2,I3,10(1PE10.3))

      return
      end subroutine apm_opt

!----------------------------------------------------------------------------------
!************************************************************************
      subroutine APM_OPT_LW(ITEST,ITYPE,MBIN,XDC,XN,NWL,WL,XBEXT,XW,XG)

      INTEGER   :: ITEST, ITYPE,MBIN,NWL
      REAL*8    :: XDC(1:MBIN),XN(1:MBIN)
      REAL*8    :: WL(NWL),XBEXT(NWL),XW(NWL),XG(NWL)

      REAL*8    :: XDCORE(MBIN),XDWET(MBIN)
      REAL*8    :: XWL,BSCATSUM,BEXTSUM,GSUM,DCORE,DWET,DSHELL
      REAL*8    :: QEXT,YW,YG,BEXT,BSCAT

      REAL*8    :: AS,ZFVSUM

      REAL*8    :: RDRY(40)
      REAL*8, PARAMETER   :: ONEPI = 3.1415926d0

      INTEGER   :: IWL,I,J

       DO I=1,MBIN
         XDCORE(I)=XDC(I)  ! um
!         XDWET(I)= XDCORE(I) * (1./ZFV(1))**(1./3.)  !um
       ENDDO

      XBEXT = 0.
      XW = 1.
      XG = 0.5
      DO IWL = 1, NWL
        XWL = WL(IWL)
        BSCATSUM = 1.d-20
        BEXTSUM = 1.d-20
        GSUM = 0.5d-20
        DO I = 1, MBIN
         IF(XN(I).GT.1.d-5) THEN
           DCORE = XDCORE(I)    !um
           DWET = DCORE
!           DWET =  XDWET(I)    !um
!           DSHELL = DWET - DCORE    !um

           CALL OPTABLE3_LW(XWL,DCORE,QEXT,YW,YG)

           BEXT = QEXT*XN(I)*DWET*DWET*1.E-8*ONEPI/4.   !cm-1 
           BSCAT = BEXT*YW
           BEXTSUM = BEXTSUM +  BEXT
           BSCATSUM = BSCATSUM + BSCAT
           GSUM = GSUM + BSCAT*YG
         ENDIF
        ENDDO
        XBEXT(IWL) = BEXTSUM             ! cm-1
        XW(IWL) = BSCATSUM/BEXTSUM
        XG(IWL) = GSUM/BSCATSUM
      ENDDO
200   FORMAT(10(1PE10.3))
201   FORMAT(I2,I3,10(1PE10.3))

      return
      end subroutine apm_opt_lw

!----------------------------------------------------------------------------------
      subroutine OPTABLE1(IWL,DSHEL,AS,BS,QEXT,YW,YG)
!
!
! INPUT:
        INTEGER :: IWL
        REAL*8  :: DCORE,DSHEL,AS,BS
! OUTPUT:
        REAL*8  :: QEXT, YW, YG

! Local variables
        REAL*8  :: Y,Z,U,V
        REAL*8  :: VOL,FRACT
        REAL*8  :: dy1,dy2,dz1,dz2,du1,du2,dv1,dv2
        REAL*8  :: dy,dz,du,dv

        INTEGER :: IDC1,IDC2,IDS1,IDS2,ISR1,ISR2,ISI1,ISI2
        INTEGER :: IDC, IDS, ISR,ISI
!
! to avoid the input values to be changed due to out of the range reset
!
!        X = XWL    ! um
!        Y = DCORE    ! um
        Z = DSHEL    ! um
        U = AS
        V = abs(BS)
!
!        IF(X.LT.WAVL(1)) THEN
!           WRITE(86,10) X, WAVL(1), WAVL(1)
!           X = WAVL(1)
!        ELSEIF(X.GT.WAVL(MWL)) THEN
!           WRITE(86,11) X, WAVL(MWL), WAVL(MWL)
!           X =WAVL(MWL)
!        ENDIF

!        IF(Y.LT.DC1(1)) THEN
!           WRITE(86,12) Y, DC1(1), DC1(1)
!           Y =DC1(1) 
!        ELSEIF(Y.GT.DC1(MDC1)) THEN
!           WRITE(86,13) Y, DC1(MDC1), DC1(MDC1)
!           Y =DC1(MDC1)
!        ENDIF

        IF(Z.LT.DS1(1)) THEN
!           WRITE(86,14) Z, DS1(1), DS1(1)
           Z =DS1(1)
        ELSEIF(Z.GT.DS1(MDS1)) THEN
!           WRITE(86,15) Z, DS1(MDS1), DS1(MDS1)
           Z =DS1(MDS1)
        ENDIF
!
        IF(U.LT.RSR(1)) THEN
!           WRITE(86,16) U, RSR(1), RSR(1)
           U =RSR(1)
        ELSEIF(U.GT.RSR(MSR)) THEN
!           WRITE(86,17) U, RSR(MSR), RSR(MSR)
           U =RSR(MSR)
        ENDIF
!
        IF(V.LT.RSI(1)) THEN
!           WRITE(86,18) V, RSI(1), RSI(1)
           V =RSI(1)
        ELSEIF(V.GT.RSI(MSI)) THEN
!           WRITE(86,19) V, RSI(MSI), RSI(MSI)
           V =RSI(MSI)
        ENDIF
!
! 10     FORMAT("OPT1 WARNING: INPUTED WAVL=",E9.3,"<",E9.3,
!     &     " set it to ",E9.3)
! 11     FORMAT("OPT1 WARNING: INPUTED WAVL=",E9.3,">",E9.3,
!     &     " set it to ",E9.3)
! 12     FORMAT("OPT1 WARNING: INPUTED DCORE =",E9.3,"<",E9.3,
!     &     " set it to ",E9.3)
! 13     FORMAT("OPT1 WARNING: INPUTED DCORE =",E9.3,">",E9.3,
!     &     " set it to ",E9.3)
! 14     FORMAT("OPT1 WARNING: INPUTED DSHEL =",E9.3,"<",E9.3,
!     &     " set it to ",E9.3)
! 15     FORMAT("OPT1 WARNING: INPUTED DSHEL =",E9.3,">",E9.3,
!     &     " set it to ",E9.3)
! 16     FORMAT("OPT1 WARNING: INPUTED RSR =",E9.3," <",E9.3,
!     &     " set it to ",E9.3)
! 17     FORMAT("OPT1 WARNING: INPUTED RSR =",E9.3," >",E9.3,
!     &     " set it to ",E9.3)
! 18     FORMAT("OPT1 WARNING: INPUTED BS =",E9.3," <",E9.3,
!     &     " set it to ",E9.3)
! 19     FORMAT("OPT1 WARNING: INPUTED BS =",E9.3," >",E9.3,
!     &     " set it to ",E9.3)

!       wavmid(ns) = 0.25 * (0.55/0.25)**(float(ns-1)/4.0)
!        IWL1 =MAX0(INT(1.+4.*LOG10(X/0.25)/LOG10(0.55/WAVL(1))),1)
!        IWL2 = MIN0(IWL1 + 1,MWL)
!        IF(IWL2.EQ.MWL) IWL1=MWL-1
        
!!        DC2(IDC) = 1.0E-2 * 10.**(float(IDC-1)/10.)
!        IDC1 =MAX0(INT(1.+10.*LOG10(Y/DC1(1))),1)
!        IDC2 = MIN0(IDC1 + 1,MDC1)
!        IF(IDC2.EQ.MDC1) IDC1=MDC1-1

!        DS2(IDS) = DSHELLMIN * 10.**(float(IDS-2)/RESS1), IDS>=1
        IDS1 = MAX0(INT(1.+RESS1*LOG10(Z/DS1(1))),1)
        IDS2 = MIN0(IDS1 + 1,MDS1)
        IF(IDS2.EQ.MDS1) IDS1=MDS1-1

!        RSR(ISR) = 1.33 + float(ISR-1)/RESR
        ISR1 = MAX0(INT(1.+(U-RSR(1))*RESR),1)
        ISR2 = MIN0(ISR1 + 1,MSR)
        IF(ISR2.EQ.MSR) ISR1=MSR-1
!
        IF(V.LT.RSI(2)) THEN
          ISI1 = 1
        ELSE
          ISI1 = MAX0(INT(2.+RESI*LOG10(V/RSI(2))),2)
        ENDIF
        ISI2 = MIN0(ISI1 + 1,MSI)
        IF(ISI2.EQ.MSI) ISI1=MSI-1
!
!	dx1 = X-WAVL(IWL1) 
!	dx2 = WAVL(IWL2)-X
!	dy1 = Y-DC1(IDC1)
!	dy2 = DC1(IDC2)-Y
	dz1 = Z-DS1(IDS1)
	dz2 = DS1(IDS2)-Z
        du1 = U - RSR(ISR1)
        du2 = RSR(ISR2) - U
        dv1 = V - RSI(ISI1)
        dv2 = RSI(ISI2) - V
!
        QEXT = 0.
        YW = 0.
        YG = 0.
!
!        VOL = (dx1+dx2)*(dy1+dy2)*(dz1+dz2)*(du1+du2)*(dv1+dv2)
!        VOL = (dx1+dx2)*(dy1+dy2)*(dz1+dz2)*(du1+du2)
!        VOL = (dx1+dx2)*(dz1+dz2)*(du1+du2)
        VOL = (dz1+dz2)*(du1+du2)*(dv1+dv2)
        DO IDS = IDS1,IDS2
          IF(IDS.EQ.IDS1) THEN
            dz = dz2
	  ELSE
            dz = dz1
          ENDIF
!      	  DO IDC = IDC1,IDC2
!            IF(IDC.EQ.IDC1) THEN
!              dy = dy2
!	    ELSE
!              dy = dy1
!            ENDIF
!            DO IWL = IWL1,IWL2
!              IF(IWL.EQ.IWL1) THEN
!                dx = dx2
!	      ELSE
!                dx = dx1
!              ENDIF

 	      DO ISR =ISR1, ISR2
                IF(ISR.EQ.ISR1) THEN
                  du = du2
	        ELSE
                  du = du1
                ENDIF
                DO ISI = ISI1, ISI2
                 IF(ISI.EQ.ISI1) THEN
                  dv = dv2
                 ELSE
                  dv = dv1
                 ENDIF

!                FRACT = dx*dy*dz*du/VOL 
!                FRACT = dx*dz*du/VOL 
                 FRACT = dz*du*dv/VOL 
                 QEXT = QEXT + FRACT*OPT1EXT(IWL,IDS,ISR,ISI)
                 YW = YW + FRACT*OPT1W(IWL,IDS,ISR,ISI)
                 YG = YG + FRACT*OPT1G(IWL,IDS,ISR,ISI)
!                WRITE(6,30)IDS,ISR,ISI,OPT1EXT(IWL,IDS,ISR,ISI),
!     &                    FRACT,
!     &                Z,DS2(IDS1),DS2(IDS2),
!     &                U,RSR(ISR1),RSR(ISR2),
!     &                V,RSI(ISI1),RSI(ISI2)
	        ENDDO
	      ENDDO
!            ENDDO
!	  ENDDO
	ENDDO
!
!
 30    FORMAT(I3, I3, I3,  20(1PE10.3))
 20    FORMAT(10(1PE10.3))

       RETURN

       END SUBROUTINE OPTABLE1

!----------------------------------------------------------------------------------

      subroutine OPTABLE2(IWL,DCORE,DSHEL,AS,BS,QEXT,YW,YG)
!
!
! INPUT:
        INTEGER :: IWL
        REAL*8  :: DCORE,DSHEL,AS,BS
! OUTPUT:
        REAL*8  :: QEXT, YW, YG

! Local variables
        REAL*8  :: Y,Z,U,V
        REAL*8  :: VOL,FRACT
        REAL*8  :: dy1,dy2,dz1,dz2,du1,du2,dv1,dv2
        REAL*8  :: dy,dz,du,dv

        INTEGER :: IDC1,IDC2,IDS1,IDS2,ISR1,ISR2,ISI1,ISI2
        INTEGER :: IDC, IDS, ISR,ISI
!
! to avoid the input values to be changed due to out of the range reset
!
!        X = XWL    ! um
        Y = DCORE    ! um
        Z = DSHEL    ! um
        U = AS
        V = BS
!
!        IF(X.LT.WAVL(1)) THEN
!           WRITE(86,10) X, WAVL(1), WAVL(1)
!           X = WAVL(1)
!        ELSEIF(X.GT.WAVL(MWL)) THEN
!           WRITE(86,11) X, WAVL(MWL), WAVL(MWL)
!           X =WAVL(MWL)
!        ENDIF

        IF(Y.LT.DC2(1)) THEN
!           WRITE(86,12) Y, DC2(1), DC2(1)
           Y =DC2(1)
        ELSEIF(Y.GT.DC2(MDC2)) THEN
!           WRITE(86,13) Y, DC2(MDC2), DC2(MDC2)
           Y =DC2(MDC2)
        ENDIF
!
        IF(Z.LT.DS2(1)) THEN
!           WRITE(86,14) Z, DS2(1), DS2(1)
           Z =DS2(1)
        ELSEIF(Z.GT.DS2(MDS2)) THEN
!           WRITE(86,15) Z, DS2(MDS2), DS2(MDS2)
           Z =DS2(MDS2)
        ENDIF
!
        IF(U.LT.RSR(1)) THEN
!           WRITE(86,16) U, RSR(1), RSR(1)
           U =RSR(1)
        ELSEIF(U.GT.RSR(MSR)) THEN
!           WRITE(86,17) U, RSR(MSR), RSR(MSR)
           U =RSR(MSR)
        ENDIF
!
        IF(V.LT.RSI(1)) THEN
!           WRITE(86,18) V, RSI(1), RSI(1)
           V =RSI(1)
        ELSEIF(V.GT.RSI(MSI)) THEN
!           WRITE(86,19) V, RSI(MSI), RSI(MSI)
           V =RSI(MSI)
        ENDIF
!
!
! 10     FORMAT("OPT2 WARNING: INPUTED WAVL=",E9.3,"<",E9.3,
!     &     " set it to ",E9.3)
! 11     FORMAT("OPT2 WARNING: INPUTED WAVL=",E9.3,">",E9.3,
!     &     " set it to ",E9.3)
! 12     FORMAT("OPT2 WARNING: INPUTED DCORE =",E9.3,"<",E9.3,
!     &     " set it to ",E9.3)
! 13     FORMAT("OPT2 WARNING: INPUTED DCORE =",E9.3,">",E9.3,
!     &     " set it to ",E9.3)
! 14     FORMAT("OPT2 WARNING: INPUTED DSHEL =",E9.3,"<",E9.3,
!     &     " set it to ",E9.3)
! 15     FORMAT("OPT2 WARNING: INPUTED DSHEL =",E9.3,">",E9.3,
!     &     " set it to ",E9.3)
! 16     FORMAT("OPT2 WARNING: INPUTED RSR =",E9.3," <",E9.3,
!     &     " set it to ",E9.3)
! 17     FORMAT("OPT2 WARNING: INPUTED RSR =",E9.3," >",E9.3,
!     &     " set it to ",E9.3)
! 18     FORMAT("OPT1 WARNING: INPUTED BS =",E9.3," <",E9.3,
!     &     " set it to ",E9.3)
! 19     FORMAT("OPT1 WARNING: INPUTED BS =",E9.3," >",E9.3,
!     &     " set it to ",E9.3)
!
!       wavmid(ns) = 0.25 * (0.55/0.25)**(float(ns-1)/4.0)
!        IWL1 =MAX0(INT(1.+4.*LOG10(X/0.25)/LOG10(0.55/WAVL(1))),1)
!        IWL2 = MIN0(IWL1 + 1,MWL)
!        IF(IWL2.EQ.MWL) IWL1=MWL-1
        
!        DC2(IDC) = 1.0E-2 * 10.**(float(IDC-1)/10.)
        IDC1 =MAX0(INT(1.+RESC2*LOG10(Y/DC2(1))),1)
        IDC2 = MIN0(IDC1 + 1,MDC2)
        IF(IDC2.EQ.MDC2) IDC1=MDC2-1

!        DS2(1) =0.
!        DS2(IDS) = DSHELLMIN * 10.**(float(IDS-2)/10.), IDS>1
        IF(Z.LT.DS2(2)) THEN
          IDS1 =1
        ELSE
          IDS1 = MAX0(INT(2.+RESS2*LOG10(Z/DS2(2))),2)
        ENDIF
        IDS2 = MIN0(IDS1 + 1,MDS2)
        IF(IDS2.EQ.MDS2) IDS1=MDS2-1

!        RSR(ISR) = 1.33 + float(ISR-1)*0.01
        ISR1 = MAX0(INT(1.+(U-RSR(1))*RESR),1)
        ISR2 = MIN0(ISR1 + 1,MSR)
        IF(ISR2.EQ.MSR) ISR1=MSR-1

!
        IF(V.LT.RSI(2)) THEN
          ISI1 = 1
        ELSE
          ISI1 = MAX0(INT(2.+RESI*LOG10(V/RSI(2))),2)
        ENDIF
        ISI2 = MIN0(ISI1 + 1,MSI)
        IF(ISI2.EQ.MSI) ISI1=MSI-1
!
!	dx1 = X-WAVL(IWL1) 
!	dx2 = WAVL(IWL2)-X
	dy1 = Y-DC2(IDC1)
	dy2 = DC2(IDC2)-Y
	dz1 = Z-DS2(IDS1)
	dz2 = DS2(IDS2)-Z
        du1 = U - RSR(ISR1)
        du2 = RSR(ISR2) - U
        dv1 = V - RSI(ISI1)
        dv2 = RSI(ISI2) - V

!
        QEXT = 0.
        YW = 0.
        YG = 0.
!
!        VOL = (dx1+dx2)*(dy1+dy2)*(dz1+dz2)*(du1+du2)*(dv1+dv2)
!        VOL = (dx1+dx2)*(dy1+dy2)*(dz1+dz2)*(du1+du2)
        VOL = (dy1+dy2)*(dz1+dz2)*(du1+du2)*(dv1+dv2)
        DO IDS = IDS1,IDS2
          IF(IDS.EQ.IDS1) THEN
            dz = dz2
	  ELSE
            dz = dz1
          ENDIF
      	  DO IDC = IDC1,IDC2
            IF(IDC.EQ.IDC1) THEN
              dy = dy2
	    ELSE
              dy = dy1
            ENDIF
!            DO IWL = IWL1,IWL2
!              IF(IWL.EQ.IWL1) THEN
!                dx = dx2
!	      ELSE
!                dx = dx1
!              ENDIF

 	      DO ISR =ISR1, ISR2
                IF(ISR.EQ.ISR1) THEN
                  du = du2
	        ELSE
                  du = du1
                ENDIF
                DO ISI = ISI1, ISI2
                 IF(ISI.EQ.ISI1) THEN
                  dv = dv2
                 ELSE
                  dv = dv1
                 ENDIF
                 FRACT = dy*dz*du*dv/VOL 
                 QEXT = QEXT + FRACT*OPT2EXT(IWL,IDC,IDS,ISR,ISI)
                 YW = YW + FRACT*OPT2W(IWL,IDC,IDS,ISR,ISI)
                 YG = YG + FRACT*OPT2G(IWL,IDC,IDS,ISR,ISI)
!                 WRITE(6,30)IDC,IDS,ISR,ISI,OPT1EXT(IWL,IDS,ISR,ISI),
!     &                FRACT,Y,DC2(IDC1),DC2(IDC2),
!     &                Z,DS2(IDS1),DS2(IDS2),
!     &                U,RSR(ISR1),RSR(ISR2),
!     &                V,RSI(ISI1),RSI(ISI2)
	        ENDDO
	      ENDDO
!            ENDDO
	  ENDDO
	ENDDO
!
!
 30    FORMAT(I3, I3, I3, I3, 20(1PE10.3))
 20    FORMAT(10(1PE10.3))

       RETURN

       END SUBROUTINE OPTABLE2
!************************************************************************
!----------------------------------------------------------------------------------
      subroutine OPTABLE3(IWL,DCORE,DSHEL,AS,BS,QEXT,YW,YG)
!
!
! INPUT:
        INTEGER :: IWL
        REAL*8  :: DCORE,DSHEL,AS,BS
! OUTPUT:
        REAL*8  :: QEXT, YW, YG

! Local variables
        REAL*8  :: Y,Z,U,V
        REAL*8  :: VOL,FRACT
        REAL*8  :: dy1,dy2,dz1,dz2,du1,du2,dv1,dv2
        REAL*8  :: dy,dz,du,dv

        INTEGER :: IDC1,IDC2,IDS1,IDS2,ISR1,ISR2,ISI1,ISI2
        INTEGER :: IDC, IDS, ISR, ISI
!
! to avoid the input values to be changed due to out of the range reset
!
!        X = XWL    ! um
        Y = DCORE    ! um
        Z = DSHEL    ! um
        U = AS
        V = BS
!
!        IF(X.LT.WAVL(1)) THEN
!           WRITE(86,10) X, WAVL(1), WAVL(1)
!           X = WAVL(1)
!        ELSEIF(X.GT.WAVL(MWL)) THEN
!           WRITE(86,11) X, WAVL(MWL), WAVL(MWL)
!           X =WAVL(MWL)
!        ENDIF

        IF(Y.LT.DC3(1)) THEN
!           WRITE(86,12) Y, DC3(1), DC3(1)
           Y =DC3(1)
        ELSEIF(Y.GT.DC3(MDC3)) THEN
!           WRITE(86,13) Y, DC3(MDC3), DC3(MDC3)
           Y =DC3(MDC3)
        ENDIF
!
        IF(Z.LT.DS3(1)) THEN
!           WRITE(86,14) Z, DS3(1), DS3(1)
           Z =DS3(1)
        ELSEIF(Z.GT.DS3(MDS3)) THEN
!           WRITE(86,15) Z, DS3(MDS3), DS3(MDS3)
           Z =DS3(MDS3)
        ENDIF
!
        IF(U.LT.RSR(1)) THEN
!           WRITE(86,16) U, RSR(1), RSR(1)
           U =RSR(1)
        ELSEIF(U.GT.RSR(MSR)) THEN
!           WRITE(86,17) U, RSR(MSR), RSR(MSR)
           U =RSR(MSR)
        ENDIF
!
        IF(V.LT.RSI(1)) THEN
!           WRITE(86,18) V, RSI(1), RSI(1)
           V =RSI(1)
        ELSEIF(V.GT.RSI(MSI)) THEN
!           WRITE(86,19) V, RSI(MSI), RSI(MSI)
           V =RSI(MSI)
        ENDIF

!       wavmid(ns) = 0.25 * (0.55/0.25)**(float(ns-1)/4.0)
!        IWL1 =MAX0(INT(1.+4.*LOG10(X/0.25)/LOG10(0.55/WAVL(1))),1)
!        IWL2 = MIN0(IWL1 + 1,MWL)
!        IF(IWL2.EQ.MWL) IWL1=MWL-1
        
!        DC2(IDC) = 1.0E-2 * 10.**(float(IDC-1)/10.)
        IDC1 =MAX0(INT(1.+RESC3*LOG10(Y/DC3(1))),1)
        IDC2 = MIN0(IDC1 + 1,MDC3)
        IF(IDC2.EQ.MDC3) IDC1=MDC3-1

!        DS2(1) =0.
!        DS2(IDS) = DSHELLMIN * 10.**(float(IDS-2)/10.), IDS>1
        IF(Z.LT.DS3(2)) THEN
          IDS1 =1
        ELSE
          IDS1 = MAX0(INT(2.+RESS3*LOG10(Z/DS3(2))),2)
        ENDIF
        IDS2 = MIN0(IDS1 + 1,MDS3)
        IF(IDS2.EQ.MDS3) IDS1=MDS3-1

!        RSR(ISR) = 1.33 + float(ISR-1)*0.01
        ISR1 = MAX0(INT(1.+(U-RSR(1))*RESR),1)
        ISR2 = MIN0(ISR1 + 1,MSR)
        IF(ISR2.EQ.MSR) ISR1=MSR-1

!
        IF(V.LT.RSI(2)) THEN
          ISI1 = 1
        ELSE
          ISI1 = MAX0(INT(2.+RESI*LOG10(V/RSI(2))),2)
        ENDIF
        ISI2 = MIN0(ISI1 + 1,MSI)
        IF(ISI2.EQ.MSI) ISI1=MSI-1
!
!	dx1 = X-WAVL(IWL1) 
!	dx2 = WAVL(IWL2)-X
	dy1 = Y-DC3(IDC1)
	dy2 = DC3(IDC2)-Y
	dz1 = Z-DS3(IDS1)
	dz2 = DS3(IDS2)-Z
        du1 = U - RSR(ISR1)
        du2 = RSR(ISR2) - U
        dv1 = V - RSI(ISI1)
        dv2 = RSI(ISI2) - V
!
        QEXT = 0.
        YW = 0.
        YG = 0.
!
!        VOL = (dx1+dx2)*(dy1+dy2)*(dz1+dz2)*(du1+du2)*(dv1+dv2)
!        VOL = (dx1+dx2)*(dy1+dy2)*(dz1+dz2)*(du1+du2)
        VOL = (dy1+dy2)*(dz1+dz2)*(du1+du2)*(dv1+dv2)
        DO IDS = IDS1,IDS2
          IF(IDS.EQ.IDS1) THEN
            dz = dz2
	  ELSE
            dz = dz1
          ENDIF
      	  DO IDC = IDC1,IDC2
            IF(IDC.EQ.IDC1) THEN
              dy = dy2
	    ELSE
              dy = dy1
            ENDIF
!            DO IWL = IWL1,IWL2
!              IF(IWL.EQ.IWL1) THEN
!                dx = dx2
!	      ELSE
!                dx = dx1
!              ENDIF

 	      DO ISR =ISR1, ISR2
                IF(ISR.EQ.ISR1) THEN
                  du = du2
	        ELSE
                  du = du1
                ENDIF
                DO ISI = ISI1, ISI2
                 IF(ISI.EQ.ISI1) THEN
                  dv = dv2
                 ELSE
                  dv = dv1
                 ENDIF

!                FRACT = dx*dy*dz*du/VOL 
                 FRACT = dy*dz*du*dv/VOL 
                 QEXT = QEXT + FRACT*OPT3EXT(IWL,IDC,IDS,ISR,ISI)
                 YW = YW + FRACT*OPT3W(IWL,IDC,IDS,ISR,ISI)
                 YG = YG + FRACT*OPT3G(IWL,IDC,IDS,ISR,ISI)
!                 WRITE(6,30)IDC,IDS,ISR,ISI,OPT1EXT(IWL,IDS,ISR,ISI),
!     &                FRACT,Y,DC3(IDC1),DC3(IDC2),
!     &                Z,DS3(IDS1),DS3(IDS2),
!     &                U,RSR(ISR1),RSR(ISR2),
!     &                V,RSI(ISI1),RSI(ISI2)
	        ENDDO
	      ENDDO
!            ENDDO
	  ENDDO
	ENDDO
!
!
 30    FORMAT(I3, I3, I3, I3, 20(1PE10.3))
 31    FORMAT(I3, I3, 100(1PE10.3))
 20    FORMAT(10(1PE10.3))

       RETURN

       END SUBROUTINE OPTABLE3
!************************************************************************
!----------------------------------------------------------------------------------
!      subroutine OPTABLE3_LW(XWL,DCORE,DSHEL,AS,QEXT,YW,YG)
      subroutine OPTABLE3_LW(XWL,DCORE,QEXT,YW,YG)
!
!
! INPUT:
        REAL*8  :: XWL,DCORE,DSHEL,AS
! OUTPUT:
        REAL*8  :: QEXT, YW, YG

! Local variables
        REAL*8  :: X,Y,Z,U
        REAL*8  :: VOL,FRACT
        REAL*8  :: dx1,dx2,dy1,dy2,dz1,dz2,du1,du2
        REAL*8  :: dx,dy,dz,du

        INTEGER :: IWL1,IWL2,IDC1,IDC2,IDS1,IDS2,ISR1,ISR2
        INTEGER :: IWL, IDC, IDS, ISR
!
! to avoid the input values to be changed due to out of the range reset
!
        X = XWL    ! um
        Y = DCORE    ! um
!        Z = DSHEL    ! um
!        U = AS
!        V = BS
!
        IF(X.LT.WAVLL(1)) THEN
!           WRITE(86,10) X, WAVLL(1), WAVLL(1)
           X = WAVLL(1)
        ELSEIF(X.GT.WAVLL(MWLL)) THEN
!           WRITE(86,11) X, WAVLL(MWLL), WAVLL(MWLL)
           X =WAVLL(MWL)
        ENDIF

        IF(Y.LT.DC3LW(1)) THEN
           WRITE(86,12) Y, DC3LW(1), DC3LW(1)
           Y =DC3LW(1)
        ELSEIF(Y.GT.DC3LW(MDCLW)) THEN
           WRITE(86,13) Y, DC3LW(MDCLW), DC3LW(MDCLW)
           Y =DC3LW(MDCLW)
        ENDIF

!        IF(Z.LT.DS3(1)) THEN
!           WRITE(86,14) Z, DS3(1), DS3(1)
!           Z =DS3(1)
!        ELSEIF(Z.GT.DS3(MDS3)) THEN
!!           WRITE(86,15) Z, DS3(MDS3), DS3(MDS3)
!           Z =DS3(MDS3)
!        ENDIF
!
!        IF(U.LT.RSR(1)) THEN
!           WRITE(86,16) U, RSR(1), RSR(1)
!           U =RSR(1)
!        ELSEIF(U.GT.RSR(MSR)) THEN
!           WRITE(86,17) U, RSR(MSR), RSR(MSR)
!           U =RSR(MSR)
!        ENDIF

!        IF(V.LT.S(1)) THEN
!           WRITE(86,18) V, S(1), S(1)
!           V =S(1)
!        ELSEIF(V.GT.S(MS)) THEN
!           WRITE(86,19) V, S(MS), S(MS)
!           V =S(MS)
!        ENDIF

 10     FORMAT("OPT3 WARNING: INPUTED WAVL=",ES9.3,"<",ES9.3,
     &     " set it to ",ES9.3)
 11     FORMAT("OPT3 WARNING: INPUTED WAVL=",ES9.3,">",ES9.3,
     &     " set it to ",ES9.3)
 12     FORMAT("OPT3 WARNING: INPUTED DCORE =",ES9.3,"<",ES9.3,
     &     " set it to ",ES9.3)
 13     FORMAT("OPT3 WARNING: INPUTED DCORE =",ES9.3,">",ES9.3,
     &     " set it to ",ES9.3)
 14     FORMAT("OPT3 WARNING: INPUTED DSHEL =",ES9.3,"<",ES9.3,
     &     " set it to ",ES9.3)
 15     FORMAT("OPT3 WARNING: INPUTED DSHEL =",ES9.3,">",ES9.3,
     &     " set it to ",ES9.3)
 16     FORMAT("OPT3 WARNING: INPUTED RSR =",ES9.3," <",ES9.3,
     &     " set it to ",ES9.3)
 17     FORMAT("OPT3 WARNING: INPUTED RSR =",ES9.3," >",ES9.3,
     &     " set it to ",ES9.3)
! 18     FORMAT("OPT3 WARNING: INPUTED S =",ES9.3," <",ES9.3,
!     &     " set it to ",ES9.3)
! 19     FORMAT("OPT3 WARNING: INPUTED S =",ES9.3," >",ES9.3,
!     &     " set it to ",ES9.3)

!       wavmid(ns) = 0.25 * (0.55/0.25)**(float(ns-1)/4.0)
!       wavmid(ns) = 4.0 * (10.)**(float(ns-1)/float(nspint-1))  !um

        IWL1 =MAX0(INT(1.+float(MWLL-1)*LOG10(X/WAVLL(1))),1)
        IWL2 = MIN0(IWL1 + 1,MWLL)
        IF(IWL2.EQ.MWLL) IWL1=MWLL-1
        
!        DC2(IDC) = 1.0E-2 * 10.**(float(IDC-1)/10.)
        IDC1 =MAX0(INT(1.+10.*LOG10(Y/DC3LW(1))),1)
        IDC2 = MIN0(IDC1 + 1,MDCLW)
        IF(IDC2.EQ.MDCLW) IDC1=MDCLW-1

!        DS2(1) =0.
!        DS2(IDS) = DSHELLMIN * 10.**(float(IDS-2)/10.), IDS>1
!        IF(Z.LT.DS3(2)) THEN
!          IDS1 =1
!        ELSE
!          IDS1 = MAX0(INT(2.+10.*LOG10(Z/DS3(2))),2)
!        ENDIF
!        IDS2 = MIN0(IDS1 + 1,MDS3)
!        IF(IDS2.EQ.MDS3) IDS1=MDS3-1

!        RSR(ISR) = 1.33 + float(ISR-1)*0.01
!        ISR1 = MAX0(INT(1.+(U-RSR(1))/0.01),1)
!        ISR2 = MIN0(ISR1 + 1,MSR)
!        IF(ISR2.EQ.MSR) ISR1=MSR-1
!
      dx1 = X-WAVLL(IWL1)
      dx2 = WAVLL(IWL2)-X
      dy1 = Y-DC3LW(IDC1)
      dy2 = DC3LW(IDC2)-Y
!      dz1 = Z-DS3(IDS1)
!      dz2 = DS3(IDS2)-Z
!        du1 = U - RSR(ISR1)
!        du2 = RSR(ISR2) - U
!
        QEXT = 0.
        YW = 0.
        YG = 0.
!
!        VOL = (dx1+dx2)*(dy1+dy2)*(dz1+dz2)*(du1+du2)*(dv1+dv2)
!        VOL = (dx1+dx2)*(dy1+dy2)*(dz1+dz2)*(du1+du2)
        VOL = (dx1+dx2)*(dy1+dy2)
!        DO IDS = IDS1,IDS2
!          IF(IDS.EQ.IDS1) THEN
!            dz = dz2
!	  ELSE
!            dz = dz1
!          ENDIF
      	  DO IDC = IDC1,IDC2
            IF(IDC.EQ.IDC1) THEN
              dy = dy2
	    ELSE
              dy = dy1
            ENDIF
            DO IWL = IWL1,IWL2
              IF(IWL.EQ.IWL1) THEN
                dx = dx2
	      ELSE
                dx = dx1
              ENDIF

! 	      DO ISR =ISR1, ISR2
!                IF(ISR.EQ.ISR1) THEN
!                  du = du2
!	        ELSE
!                  du = du1
!                ENDIF
!                FRACT = dx*dy*dz*du/VOL 
!                QEXT = QEXT + FRACT*OPT3EXT(IWL,IDC,IDS,ISR)
!                YW = YW + FRACT*OPT3W(IWL,IDC,IDS,ISR)
!                YG = YG + FRACT*OPT3G(IWL,IDC,IDS,ISR)
                FRACT = dx*dy/VOL 
                QEXT = QEXT + FRACT*OPT3EXT_LW(IWL,IDC)
                YW = YW + FRACT*OPT3W_LW(IWL,IDC)
                YG = YG + FRACT*OPT3G_LW(IWL,IDC)
!                WRITE(6,31)IWL,IDC,OPT3EXT_LW(IWL,IDC),QEXT,
!     &                    OPT3W_LW(IWL,IDC),YW,FRACT
!	      ENDDO
            ENDDO
	  ENDDO
!	ENDDO
!
!
 30    FORMAT(I3, I3, I3, I3, I3, 10(1PE10.3))
 31    FORMAT(I3, I3, 100(1PE10.3))
 20    FORMAT(10(1PE10.3))

       RETURN

       END SUBROUTINE OPTABLE3_LW
!************************************************************************
!------------------------------------------------------------------------------
      subroutine READOPTABLE(ITABLE,DATA_DIR_1x1)

      IMPLICIT NONE

      integer :: ITABLE
      CHARACTER(LEN=255)   :: DATA_DIR_1x1
      CHARACTER*999 YPATH,YPATH1

      INTEGER :: IML,IDC,IDS,ISR,ISI
      INTEGER :: IML1,IDC1,IDS1,ISR1, ISI1,MDC, MDS
      REAL*8  :: WAVL1,RSR1,RSI1
      REAL*8  :: YREAL,YIM
      REAL*8  :: DCa,DCb,DSa,DSb
      REAL*8  :: YTEMP
      REAL*8  :: qextc,w,gscac
      REAL*8  :: RESC, RESS
      REAL*8  :: DCORE1, DSHELL1, DSHELL2   ! in um

      LOGICAL, SAVE :: FIRST = .TRUE.

      CLOSE(88)
      YPATH = TRIM(DATA_DIR_1x1)//'/APM_data_201906/OPTAB20120823/'    !for RRTMG WL

      IF(ITABLE.EQ.1)THEN      !no core
       MDC = MDC1
       MDS = MDS1
       OPEN(88,file=TRIM(YPATH)//'OPTABLE_NoCore.txt',status='old')
       WRITE(6,*)"read OPTABLE_NoCore.txt"
      ELSEIF(ITABLE.EQ.2)THEN   !BC core
       MDC = MDC2
       MDS = MDS2
       OPEN(88,file=TRIM(YPATH)//'OPTABLE_BC.txt',status='old')
       WRITE(6,*)"read OPTABLE_BC.txt"
      ELSEIF(ITABLE.EQ.3)THEN   !DUST core
       MDC = MDC3
       MDS = MDS3
!       OPEN(88,file=TRIM(YPATH1)//'OPTABLE_DUST.txt',status='old')
       OPEN(88,file=TRIM(YPATH)//'OPTABLE_DUST.txt',status='old')
       WRITE(6,*)"read OPTABLE_DUST.txt"
      ENDIF

      READ(88,*)
      READ(88,101)IML,IDC,IDS,ISR,ISI
101   format(3x,I3,3x,I3,2x,I3,2x,I3,2x,I3)

      IF(IML.NE.MWL.or.IDC.NE.MDC.or.IDS.NE.MDS.
     &           or.ISR.NE.MSR.or.ISI.NE.MSI) THEN
         WRITE(6,*)"STOP: NEED to check IML,IDC,IDS,ISR values"
         WRITE(6,*)IML,IDC,IDS,ISR,ISI
         WRITE(6,*)MWL,MDC,MDS,MSR,MSI
         STOP
      ENDIF

!     "*** resolution ***"
      READ(88,*)
      READ(88,102)RESR, RESI, RESC, RESS
102   format(10F6.1)

!     "*** first 2 bin diameter of core and shell ***"
      READ(88,*)
      READ(88,103)DCORE1, DSHELL1, DSHELL2   ! in um
103   format(10(1PE10.2))

      IF(FIRST) THEN   ! RSR and RSI same for all three types, need to modify if not
       DO ISR=1,MSR
        RSR(ISR) = 1.33 + float(ISR-1)/RESR
       ENDDO
       RSI(1)=1.d-6
       DO ISI=2,MSI
        RSI(ISI) = 1.d-3*10.**(float(ISI-2)/RESI)
       ENDDO
       FIRST = .FALSE.
      ENDIF

!     "***wavelength (um)***"
      READ(88,*)
      do IML=1,MWL
       IF(ITABLE.EQ.1) THEN
        read(88,110)IML1,WAVL1  !um
       ELSEIF(ITABLE.EQ.2.or.ITABLE.EQ.3) THEN
        read(88,110)IML1,WAVL1,YREAL,YIM 
       ENDIF

!       YTEMP = abs(WAVL1-WAVL(IML))/WAVL1  !double check
!       IF(YTEMP.GT.0.01) THEN
!         WRITE(6,*)"STOP: NEED to check WAVL",WAVL1,WAVL(IML)
!       ENDIF
      enddo
!
! Define core and shell diameter
! DWET = DSHELL + DCORE

!       "***Dcore (um)***"
      READ(88,*)
      IF(ITABLE.EQ.1) THEN
       RESC1 = RESC
       do IDC = 1, MDC
        DC1(IDC) = DCORE1 * 10.**(float(IDC-1)/RESC1)   !um
       enddo
      ELSEIF(ITABLE.EQ.2) THEN
       RESC2 = RESC
       do IDC = 1, MDC
        DC2(IDC) = DCORE1 * 10.**(float(IDC-1)/RESC2)
       enddo
      ELSEIF(ITABLE.EQ.3) THEN
       RESC3 = RESC
       do IDC = 1, MDC
        DC3(IDC) = DCORE1 * 10.**(float(IDC-1)/RESC3)
       enddo
      ELSE
       WRITE(6,*) "STOP: NEED to check ITABLE.", ITABLE
       STOP
      ENDIF

      do IDC = 1, MDC
       READ(88,110)IDC1,DCa  !um
       IF(MDC.GT.1) THEN
        IF(ITABLE.EQ.1) DCb=DC1(IDC)
        IF(ITABLE.EQ.2) DCb=DC2(IDC)
        IF(ITABLE.EQ.3) DCb=DC3(IDC)
        YTEMP = abs(DCa-DCb)/max(1.d-30,DCa)  !double check
        IF(YTEMP.GT.0.01) THEN
          WRITE(6,*)"STOP: NEED to check DC",DCa,DCb
        ENDIF
       ENDIF
      enddo
!
!      if(MDC.GT.1) THEN
!        READ(88,*)
!      endif

!    "***Dshell (um)***"
      READ(88,*)
      IF(ITABLE.EQ.1) THEN
       RESS1 = RESS
       do IDS = 1, MDS
        DS1(IDS) = DSHELL1 * 10.**(float(IDS-1)/RESS1)   !um
       enddo
      ELSEIF(ITABLE.EQ.2) THEN
       RESS2 = RESS
       DS2(1) = DSHELL1
       do IDS = 2, MDS
        DS2(IDS) = DSHELL2 * 10.**(float(IDS-2)/RESS2)
       enddo
      ELSEIF(ITABLE.EQ.3) THEN
       RESS3 = RESS
       DS3(1) = DSHELL1
       do IDS = 2, MDS
        DS3(IDS) = DSHELL2 * 10.**(float(IDS-2)/RESS3)
       enddo
      ELSE
       WRITE(6,*) "STOP: NEED to check ITABLE.", ITABLE
       STOP
      ENDIF

      do IDS = 1, MDS
        READ(88,110)IDS1,DSa   !um
        IF(IDS1.GT.0.) THEN
         IF(ITABLE.EQ.1) DSb=DS1(IDS)
         IF(ITABLE.EQ.2) DSb=DS2(IDS)
         IF(ITABLE.EQ.3) DSb=DS3(IDS)
         YTEMP = abs(DSa-DSb)/max(1.d-30,DSa)  !double check
         IF(IDS.GT.1) THEN
         IF(YTEMP.GT.0.01) THEN
          WRITE(6,*)"STOP: NEED to check DS",DSa,DSb
         ENDIF
         ENDIF
        ENDIF
      enddo

! "***realshell***"
      READ(88,*)
      do ISR=1,MSR
        READ(88,110)ISR1, RSR1
        YTEMP = abs(RSR1-RSR(ISR))/max(1.d-30,RSR1)  !double check
        IF(YTEMP.GT.0.01) THEN
         WRITE(6,*)"STOP: NEED to check REALSH",RSR1,RSR(ISR)
        ENDIF
      enddo

! "***imagshell***"
      READ(88,*)
      do ISI=1,MSI
        READ(88,110)ISI1, RSI1
        YTEMP = abs(RSI1-RSI(ISI))/max(1.d-30,RSI1)  !double check
        IF(YTEMP.GT.0.01) THEN
         WRITE(6,*)"STOP: NEED to check IMAGSH",RSI1,RSI(ISI)
        ENDIF
      enddo

      ! Begin spectral loop
      do IML=1,MWL
       do IDC=1,MDC
        do IDS = 1, MDS
         do ISR=1,MSR
          do ISI=1,MSI
           READ(88,122)IML1,IDC1,IDS1,ISR1,ISI1,qextc,w,gscac
           IF(IML1.NE.IML.or.IDC1.NE.IDC.or.IDS1.NE.IDS.
     &                 or.ISR1.NE.ISR.or.ISI1.NE.ISI) THEN
            WRITE(6,*)"STOP: Need to check IML1,IDC1,IDS1,ISR1,ISI1"
            STOP
           ENDIF
           IF(ITABLE.EQ.1) THEN
            OPT1EXT(IML,IDS,ISR,ISI) = qextc 
            OPT1W(IML,IDS,ISR,ISI) =  w
            OPT1G(IML,IDS,ISR,ISI) =  gscac
           ELSEIF(ITABLE.EQ.2) THEN
            OPT2EXT(IML,IDC,IDS,ISR,ISI) = qextc 
            OPT2W(IML,IDC,IDS,ISR,ISI) =  w   
            OPT2G(IML,IDC,IDS,ISR,ISI) =  gscac
           ELSEIF(ITABLE.EQ.3) THEN
            OPT3EXT(IML,IDC,IDS,ISR,ISI) = qextc 
            OPT3W(IML,IDC,IDS,ISR,ISI) =  w   
            OPT3G(IML,IDC,IDS,ISR,ISI) =  gscac
           ENDIF
          enddo   !k -- lshell
         enddo   !l -- realshell
        enddo   !MDS
       enddo   !MDC
      enddo    !wavelength

110   FORMAT(I3,10(1PE11.3))
120   FORMAT(10(1PE9.2))
121   FORMAT(2(1PE10.3),10(1PE10.3))
122   FORMAT(I2,4I3,10(1PE10.3))
      RETURN

      END SUBROUTINE READOPTABLE
! *****************************************************************************

!------------------------------------------------------------------------------
      subroutine READOPTABLE_LW(ITABLE,DATA_DIR_1x1)

      IMPLICIT NONE

      integer :: ITABLE
      CHARACTER(LEN=255)   :: DATA_DIR_1x1
      CHARACTER*999 YPATH,YPATH1

      INTEGER :: IML,IDC,IDS,ISR
      INTEGER :: IML1,IDC1,IDS1,ISR1, MDC, MDS
      REAL*8  :: WAVL1,RSR1
      REAL*8  :: DCa,DCb,DSa,DSb
      REAL*8  :: DCOREMIN,DSHELLMIN,YTEMP
      REAL*8  :: qextc,w,gscac

      LOGICAL, SAVE :: FIRST = .TRUE.
      INTEGER, PARAMETER :: MSRL = 1

      IF(FIRST) THEN
       DO IML=1,MWLL
        WAVLL(IML) = 4.0 * (10.)**(float(IML-1)/float(MWLL-1))  ! um
       ENDDO
       DO ISR=1,MSRL
        RSR(ISR) = 1.33 + float(ISR-1)*0.01
       ENDDO
       FIRST = .FALSE.
      ENDIF

      CLOSE(88)
      YPATH1 = TRIM(DATA_DIR_1x1)//'/APM_data_201906/OPTAB20110418/'

      IF(ITABLE.EQ.1)THEN      !no core
!       MDC = MDC1
!       MDS = MDS1
!       DCOREMIN = 0.0E-0   ! um
!       DSHELLMIN = 2.0E-2   ! um
!       OPEN(88,file=TRIM(YPATH)//'OPTABLE_NoCore.txt',status='old')
!       WRITE(6,*)"read OPTABLE_NoCore.txt"
      ELSEIF(ITABLE.EQ.2)THEN   !BC core
!       MDC = MDC2
!       MDS = MDS2
!       DCOREMIN =  1.0E-2   ! um
!       DSHELLMIN = 1.0E-3   ! um
!       OPEN(88,file=TRIM(YPATH)//'OPTABLE_BC.txt',status='old')
!       WRITE(6,*)"read OPTABLE_BC.txt"
      ELSEIF(ITABLE.EQ.3)THEN   !DUST core
       MDC = MDCLW
!       MDS = MDS3
       MDS = 1
       DCOREMIN = 5.0E-2   ! um
       DSHELLMIN = 1.0E-3   ! um
       OPEN(88,file=TRIM(YPATH1)//'OPTABLE_DUST_LW.txt',status='old')
       WRITE(6,*)"read OPTABLE_DUST_LW.txt"
      ENDIF

      READ(88,*)
      READ(88,101)IML,IDC,IDS,ISR
101   format(3x,I3,3x,I3,2x,I3,2x,I3)

      IF(IML.NE.MWLL.or.IDC.NE.MDC.or.IDS.NE.MDS.or.ISR.NE.MSRL) THEN
         WRITE(6,*)"STOP: NEED to check IML,IDC,IDS,ISR values"
         WRITE(6,*)IML,IDC,IDS,ISR
         WRITE(6,*)MWLL,MDC,MDS,MSRL
         STOP
      ENDIF

!     "***wavelength (um)***"
      READ(88,*)
      do IML=1,MWLL
       read(88,110)IML1,WAVL1  !um
       YTEMP = abs(WAVL1-WAVLL(IML))/max(1.d-30,WAVL1)  !double check
       IF(YTEMP.GT.0.01) THEN
         WRITE(6,*)"STOP: NEED to check WAVL",WAVL1,WAVLL(IML)
       ENDIF
      enddo
!
! Define core and shell diameter
! DWET = DSHELL + DCORE

!       "***Dcore (um)***"
      READ(88,*)
      IF(ITABLE.EQ.1) THEN
!       do IDC = 1, MDC
!        DC1(IDC) = DCOREMIN * 10.**(float(IDC-1)/10.)   !um
!       enddo
      ELSEIF(ITABLE.EQ.2) THEN
!       do IDC = 1, MDC
!        DC2(IDC) = DCOREMIN * 10.**(float(IDC-1)/10.)
!       enddo
      ELSEIF(ITABLE.EQ.3) THEN
       do IDC = 1, MDC
        DC3LW(IDC) = DCOREMIN * 10.**(float(IDC-1)/10.)
       enddo
      ELSE
       WRITE(6,*) "STOP: NEED to check ITABLE.", ITABLE
       STOP
      ENDIF

      do IDC = 1, MDC
       READ(88,110)IDC1,DCa  !um
       IF(MDC.GT.1) THEN
        IF(ITABLE.EQ.1) DCb=DC1(IDC)
        IF(ITABLE.EQ.2) DCb=DC2(IDC)
        IF(ITABLE.EQ.3) DCb=DC3LW(IDC)
        YTEMP = abs(DCa-DCb)/DCa  !double check
        IF(YTEMP.GT.0.01) THEN
          WRITE(6,*)"STOP: NEED to check DC",DCa,DCb
          STOP
        ENDIF
       ENDIF
      enddo
!
      if(MDC.GT.1) THEN
        READ(88,*)
      endif

!    "***Dshell (um)***"
      READ(88,*)
      IF(ITABLE.EQ.1) THEN
!       do IDS = 1, MDS
!        DS1(IDS) = DSHELLMIN * 10.**(float(IDS-1)/10.)   !um
!       enddo
      ELSEIF(ITABLE.EQ.2) THEN
!       DS2(1) = 0.
!       do IDS = 2, MDS
!        DS2(IDS) = DSHELLMIN * 10.**(float(IDS-2)/10.)
!       enddo
      ELSEIF(ITABLE.EQ.3) THEN
       DS3(1) = 0.
       do IDS = 2, MDS
        DS3(IDS) = DSHELLMIN * 10.**(float(IDS-2)/10.)
       enddo
      ELSE
       WRITE(6,*) "STOP: NEED to check ITABLE.", ITABLE
       STOP
      ENDIF

      do IDS = 1, MDS
        READ(88,110)IDS1,DSa   !um
        IF(IDS1.GT.0.) THEN
         IF(ITABLE.EQ.1) DSb=DS1(IDS)
         IF(ITABLE.EQ.2) DSb=DS2(IDS)
         IF(ITABLE.EQ.3) DSb=DS3(IDS)
         YTEMP = abs(DSa-DSb)/max(1.d-30,DSa)  !double check
         IF(YTEMP.GT.0.01) THEN
          WRITE(6,*)"STOP: NEED to check DS",DSa,DSb
         ENDIF
        ENDIF
      enddo

! "***realshell***"
      READ(88,*)
      do ISR=1,MSRL
        READ(88,110)ISR1, RSR1
        YTEMP = abs(RSR1-RSR(ISR))/max(1.d-30,RSR1)  !double check
        IF(YTEMP.GT.0.01) THEN
         WRITE(6,*)"STOP: NEED to check REALSH",RSR1,RSR(ISR)
        ENDIF
      enddo

      ! Begin spectral loop
      do IML=1,MWLL
       do IDC=1,MDC
        do IDS = 1, MDS
         do ISR=1,MSRL
          READ(88,122)IML1,IDC1,IDS1,ISR1,qextc,w,gscac
          IF(IML1.NE.IML.or.IDC1.NE.IDC.or.IDS1.NE.IDS.
     &                 or.ISR1.NE.ISR) THEN
            WRITE(6,*)"STOP: Need to check IML1,IDC1,IDS1,ISR1"
            STOP
          ENDIF
          IF(ITABLE.EQ.1) THEN
!            OPT1EXT(IML,IDS,ISR) = qextc 
!            OPT1W(IML,IDS,ISR) =  w
!            OPT1G(IML,IDS,ISR) =  gscac
          ELSEIF(ITABLE.EQ.2) THEN
!            OPT2EXT(IML,IDC,IDS,ISR) = qextc 
!            OPT2W(IML,IDC,IDS,ISR) =  w   
!            OPT2G(IML,IDC,IDS,ISR) =  gscac
          ELSEIF(ITABLE.EQ.3) THEN
            OPT3EXT_LW(IML,IDC) = qextc 
            OPT3W_LW(IML,IDC) =  w   
            OPT3G_LW(IML,IDC) =  gscac
          ENDIF
         enddo   !l -- realshell
        enddo   !MDS
       enddo   !MDC
      enddo    !wavelength

110   FORMAT(I3,10(1PE10.3))
120   FORMAT(100(1PE9.2))
121   FORMAT(2(1PE10.3),10(1PE10.3))
122   FORMAT(I2,3I3,10(1PE10.3))
      RETURN

      END SUBROUTINE READOPTABLE_LW
! *****************************************************************************
      END MODULE APM_OPTI_MOD
!----------------------------------------------------------------------------------
#endif
