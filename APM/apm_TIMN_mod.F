#ifdef APM
! $Id: apm_nucl_mod.f,v 0.0 2008/09/28 11:30:00 fyu $
      MODULE APM_TIMN_MOD
!
!******************************************************************************
!  Module APM_TIMN_MOD contains variables and routines for computing 
!  ternary IMN rates. (fyu,2018)
!
!  Module Variables:
!  ============================================================================
!  Parameters
!  (1 ) MC   : NUMBER OF POINTS IN H2SO4 CONCENTRATION DIMENSION
!  (2 ) MT   : NUMBER OF POINTS IN TEMPERATURE DIMENSION
!  (3 ) MRH  : NUMBER OF POINTS IN RELATIVE HUMIDITY DIMENSION
!  (4 ) MQ   : NUMBER OF POINTS IN IONIZATION RATE DIMENSION
!  (5 ) MS   : NUMBER OF POINTS IN SURFACE AREA DIMENSION
!  (6 ) MB   : NUMBER OF POINTS IN [NH3] DIMENSION
!  Arrays 
!  (1 ) C   : VALUES AT POINTS IN H2SO4 CONCENTRATION DIMENSION
!  (2 ) T   : VALUES AT POINTS IN TEMPERATURE DIMENSION
!  (3 ) RH  : VALUES AT POINTS IN RELATIVE HUMIDITY DIMENSION
!  (4 ) Q   : VALUES AT POINTS IN IONIZATION RATE DIMENSION
!  (5 ) S   : VALUES AT POINTS IN SURFACE AREA DIMENSION
!  (6 ) S   : VALUES AT POINTS IN [NH3] DIMENSION
!
!  (11) XJTIMN : ION-MEDIATED NUCLEATION RATES (cm-3s-1) AT ALL POINTS IN 6-D SPACE
!  (12) XRSTAR : CRITICAL RADIUS (nm) AT ALL POINTS IN 6-DIMENSION SPACE
!
!  Module Routines:
!  ============================================================================
!  (1 ) YUJTIMN1     : INTERPOLAION SCHEME TO FIND JTIMN FROM LOOKUP TABLE
!  (2 ) READJTIMN   : READ IN THE TIMN LOOKUP TABLE 
!
!  NOTES:
!  (1 ) .... 
!******************************************************************************
!
      IMPLICIT NONE

      !=================================================================
      ! MODULE PRIVATE DECLARATIONS -- keep certain internal variables 
      ! and routines from being seen outside "apm_nucl_mod.f"
      !=================================================================

      ! Make everything PRIVATE ...
      PRIVATE

      ! ... except these variables ...
!      PUBLIC :: 

      ! ... and these routines
      PUBLIC :: YUJTIMN1
      PUBLIC :: READJTIMN
      !=================================================================
      ! MODULE VARIABLES
      !=================================================================
      ! Parameters
      INTEGER, PARAMETER   :: MC  = 32
      INTEGER, PARAMETER   :: MT  = 39
      INTEGER, PARAMETER   :: MRH = 26
      INTEGER, PARAMETER   :: MQ  = 9
!      INTEGER, PARAMETER   :: MS  = 12
      INTEGER, PARAMETER   :: MS  = 2
      INTEGER, PARAMETER   :: MB  = 33

      ! Arrays 
      REAL*8        :: C(MC),RH(MRH),T(MT),Q(MQ),S(MS),B(MB)
      REAL*8        :: XJTIMN(MC,MRH,MT,MQ,MS,MB), XRSTAR(MC,MRH,MT,MB)
  
      CHARACTER(LEN=255)   :: DATA_DIR_1x1
      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement
      !=================================================================
      CONTAINS

!------------------------------------------------------------------------------

      SUBROUTINE YUJTIMN1(X0,Y0,Z0,U0,V0,W0,XJBH,XJBIM,XJTH,XJTIM,
     &                    RBH,RBIM,RTH,RTIM)
!
! This subroutine is to calculate rates and critical cluster properties of
! ion-mediated nucleation (IMN) and kinetic binary homogeneous nucleation (KBHN) 
! from lookup tables using multiple-variable interpolation scheme. 
!
! Here the IMN lookup table reported in Yu (JGR,2010) and BHN lookup table 
! reported in Yu (JGR, 2008) have been integrated. YUJTIMN gives KBHN rates
! when the ionization rate is set to 0. The integrated lookup table has a size
! of ~ 170 MB and is designed for 3-D application. For quick application and 
! comparison, one can obtain nucleation rates under specified conditions using 
! an online nucleation rate calculator @ 
! http://www.albany.edu/~yfq/YuOnLineNucleation.html.
!
! The present lookup table should cover almost all the possible conditions in 
! the troposphere relevant to atmospheric nucleation. The range and resolution 
! in each parameter space can be extended in the future if needed.
!
! Written by 
! Fangqun Yu
! Atmospheric Sciences Research Center
! State University of New York at Albany
! E-mail: yfq@asrc.cestm.albany.edu; fangqun.yu@asrc.albany.edu
!
! Original code writted in 2006. Significnat update in 2008 and 2010. Contact
! Yu for future update or if you have questions.
!
! IMN lookup table reference: 
! 1. Yu, F., Ion-mediated nucleation in the atmosphere: Key controlling 
!      parameters, implications, and look-up table, J. Geophy. Res., 115, 
!      D03206, doi:10.1029/2009JD012630, 2010.
!
! IMN model references: 
! 2. Yu, F., From molecular clusters to nanoparticles: Second-generation 
!      ion-mediated nucleation model, Atmos. Chem. Phys., 6, 5193-5211, 2006.
! 3. Yu, F., and R. P. Turco, Ultrafine aerosol formation via ion-mediated 
!      nucleation, Geophys. Res. Lett., 27, 883-886, 2000.
!
! KBHN lookup table reference:
! 4. Yu, F., Updated H2SO4-H2O binary homogeneous nucleation rate look-up 
!      tables, J. Geophy. Res.,113, D24201, doi:10.1029/2008JD010527, 2008.
!
! KBHN model references:
! 5. Yu, F., Improved quasi-unary nucleation model for binary H2SO4-H2O 
!      homogeneous nucleation, J. Chem. Phys., 127, 054301, 2007.
! 6. Yu, F., Quasi-unary homogeneous nucleation of H2SO4-H2O, J. Chem. 
!      Phys., 122, 074501, 2005.
!
! INPUT (valid value range):
! X0 = [H2SO4] in #/cm3  (5E5-5E8,5E9)
! Y0 = RH in % (0.5-99.5)
! Z0 = T (in K) (190-301)
! U0 = Q = ionization rate (ion-pairs/cm3s) (0, 2-22.8,100)
! V0 = S = surface area (um2/cm3) (1, 10-1000)
! W0 = B = NH3 in #/cm3  (1E5,1E8-1E11,1E12)
!
! OUTPUT:
! XJBH: Bionary Homo  Nucleation rate (#/cm3s)
! XJBIM: Bionary Ion-Mediated Nucleation rate (#/cm3s)
! XJTH: Ternary Homo  Nucleation rate (#/cm3s)
! XJ0: Ternary Ion-Mediated Nucleation rate (#/cm3s)
! XR0: Radius of critical cluster (nm)
!
        REAL*8  :: X0,Y0,Z0,U0,V0,W0
        REAL*8  :: XJBH,XJBIM,XJTH,XJTIM
        REAL*8  :: RBH,RBIM,RTH,RTIM
        REAL*8  :: X,Y,Z,U,V,W
        REAL*8  :: VOL,VOL4,FRACT,FRACT4,XR0
        REAL*8  :: dx1,dx2,dy1,dy2,dz1,dz2,du1,du2,dv1,dv2,dw1,dw2
        REAL*8  :: dx,dy,dz,du,dv,dw
        REAL*8  :: XDH,XDT,X1,X2,Y1,Y2,XJ1,XJ2,YJ

        INTEGER :: IC1,IC2,JRH1,JRH2,KT1,KT2,IQ1,IQ2,IS1,IS2,IB1,IB2
        INTEGER :: IC, JRH, KT, IQ, IS,IB
!
! to avoid the input values to be changed due to out of the range reset
!
        X = X0
        Y = Y0
        Z = Z0
        U = U0
        V = V0
        W = W0
!
! The present lookup table should cover almost all the possible conditions in 
! the troposphere relevant to atmospheric nucleation. The range and resolution 
! in each parameter space can be extended in the future if needed.
! If the inputed values are out of the lookup table valid ranges, set them to 
! boundary values for now. Care should be taken if your inputted values are 
! frequently out of the specified ranges.
! 
        IF(U.LE.1.E-20) U=1.E-20    ! i.e., binary homogeneous nucleation 
!
        IF(X.LT.C(1)) THEN
!           WRITE(86,10) X, C(1), C(1)
           X = C(1)
        ELSEIF(X.GT.C(MC)) THEN
!           WRITE(86,11) X, C(MC), C(MC)
           X =C(MC)
        ENDIF

        IF(Y.LT.RH(1)) THEN
!           WRITE(86,12) Y, RH(1), RH(1)
           Y =RH(1) 
        ELSEIF(Y.GT.RH(MRH)) THEN
!           WRITE(86,13) Y, RH(MRH), RH(MRH)
           Y =RH(MRH)
        ENDIF

        IF(Z.LT.T(1)) THEN
!           WRITE(86,14) Z, T(1), T(1)
           Z =T(1)
        ELSEIF(Z.GT.T(MT)) THEN
!           WRITE(86,15) Z, T(MT), T(MT)
           Z =T(MT)
        ENDIF

        IF(U.LT.Q(1)) THEN
!           WRITE(86,16) U, Q(1), Q(1)
           U =Q(1)
        ELSEIF(U.GT.Q(MQ)) THEN
!           WRITE(86,17) U, Q(MQ), Q(MQ)
           U =Q(MQ)
        ENDIF

!        IF(V.LT.S(1)) THEN
!!           WRITE(86,18) V, S(1), S(1)
!           V =S(1)
!        ELSEIF(V.GT.S(MS)) THEN
!!           WRITE(86,19) V, S(MS), S(MS)
!           V =S(MS)
!        ENDIF

        IF(W.LT.B(1)) THEN
!           WRITE(86,20) W, B(1), B(1)
           W =B(1)
        ELSEIF(W.GT.B(MB)) THEN
!           WRITE(86,21) W, B(MB), B(MB)
           W =B(MB)
        ENDIF


 10     FORMAT("IMN WARNING: INPUTED [H2SO4]=",ES9.2,"<",ES9.2,
     &     " set it to ",ES9.2)
 11     FORMAT("IMN WARNING: INPUTED [H2SO4]=",ES9.2,">",ES9.2,
     &     " set it to ",ES9.2)
 12     FORMAT("IMN WARNING: INPUTED RH =",F5.1,"% <",F5.1,
     &     "% set it to ",F5.1,"%")
 13     FORMAT("IMN WARNING: INPUTED RH =",F5.1,"% >",F5.1,
     &     "% set it to ",F5.1,"%")
 14     FORMAT("IMN WARNING: INPUTED T =",F6.1,"K <",F6.1,
     &     "K set it to ",F6.1,"K")
 15     FORMAT("IMN WARNING: INPUTED T =",F6.1,"K >",F6.1,
     &     "K set it to ",F6.1,"K")
 16     FORMAT("IMN WARNING: INPUTED Q =",F6.1," <",F6.1,
     &     " ion-pair/cm3s set it to ",F6.1)
 17     FORMAT("IMN WARNING: INPUTED Q =",F6.1," >",F6.1,
     &     " ion-pair/cm3s set it to ",F6.1)
 18     FORMAT("IMN WARNING: INPUTED S =",F6.1," <",F6.1,
     &     " um2/cm3 set it to ",F6.1)
 19     FORMAT("IMN WARNING: INPUTED S =",F6.1," >",F6.1,
     &     " um2/cm3 set it to ",F6.1)
 20     FORMAT("IMN WARNING: INPUTED [NH3]=",ES9.2,"<",ES9.2,
     &     " set it to ",ES9.2)
 21     FORMAT("IMN WARNING: INPUTED [NH3]=",ES9.2,">",ES9.2,
     &     " set it to ",ES9.2)

        IC1 =MAX0(INT(1.+10.*LOG10(X/5.E5)),1)
        IC2 = MIN0(IC1 + 1,MC)
        IF(IC2.EQ.MC) IC1=MC-1
        
        XDH = 4.
        IF(Y.LT.RH(2)) THEN
           JRH1 = 1.
        ELSE
         JRH1 = MAX0(INT((Y-RH(2))/XDH+2.),2)
        ENDIF
        JRH2 = MIN0(JRH1 + 1,MRH)
        IF(JRH2.EQ.MRH) JRH1=MRH-1

        XDT = 3.0
        KT1 = MAX0(INT((Z-190.0)/XDT)+1,1)
        KT2 = MIN0(KT1 + 1,MT)
        IF(KT2.EQ.MT) KT1=MT-1
!
        
        IF(U.LT.Q(2)) THEN
          IQ1 =1.
        ELSE
          IQ1 = MAX0(INT(2.+LOG10(U/Q(2))/LOG10(1.5)),2)
        ENDIF
        IQ2 = MIN0(IQ1 + 1,MQ)
        IF(IQ2.EQ.MQ) IQ1=MQ-1
!
!        IF(V.LT.10.0) THEN
!          IS1 =1.
!        ELSE
!          IS1 = MAX0(INT(2.+5.*LOG10(V/10.)),2)
!        ENDIF
!        IS2 = MIN0(IS1 + 1,MS)
!        IF(IS2.EQ.MS) IS1=MS-1

        IF(W.LT.B(2)) THEN
          IB1 =1.
        ELSE
          IB1 = MAX0(INT(2.+10.*LOG10(W/B(2))),2)
        ENDIF
        IB2 = MIN0(IB1 + 1,MB)
        IF(IB2.EQ.MB) IB1=MB-1
!
	dx1 = LOG10(X/C(IC1))   ! logJ log[H2SO4] interpolation
	dx2 = LOG10(C(IC2)/X)
	dy1 = LOG10(Y/RH(JRH1))
	dy2 = LOG10(RH(JRH2)/Y)
	dz1 = Z-T(KT1)
	dz2 = T(KT2)-Z

        du1 = U - Q(IQ1)
        du2 = Q(IQ2) - U
!        dv1 = V- S(IS1)
!        dv2 = S(IS2) - V
        dw1 = LOG10(W/B(IB1))   ! logJ log[H2SO4] interpolation
        dw2 = LOG10(B(IB2)/W)
!
!JTIMN
        XJ1 = 0.  
        XJ2 = 0.
        XR0 = 0.
!
        VOL = (dx1+dx2)*(dy1+dy2)*(dz1+dz2)*(du1+du2)*(dw1+dw2)
        VOL4 = (dx1+dx2)*(dy1+dy2)*(dz1+dz2)*(dw1+dw2)

        DO KT = KT1,KT2
          IF(KT.EQ.KT1) THEN
            dz = dz2
	  ELSE
            dz = dz1
          ENDIF
      	  DO JRH = JRH1,JRH2
            IF(JRH.EQ.JRH1) THEN
              dy = dy2
	    ELSE
              dy = dy1
            ENDIF
            DO IC = IC1,IC2
              IF(IC.EQ.IC1) THEN
                dx = dx2
	      ELSE
                dx = dx1
              ENDIF
 	      DO IB =IB1, IB2
                IF(IB.EQ.IB1) THEN
                  dw = dw2
	        ELSE
                  dw = dw1
                ENDIF

                FRACT4 = dx*dy*dz*dw/VOL4
                XR0 = XR0 + FRACT4*XRSTAR(IC,JRH,KT,IB)

	        DO IQ =IQ1, IQ2
                  IF(IQ.EQ.IQ1) THEN
                    du = du2
	          ELSE
                    du = du1
                  ENDIF
                  FRACT = dx*dy*dz*du*dw/VOL 
                  XJ1 = XJ1 + FRACT*XJTIMN(IC,JRH,KT,IQ,1,IB)
                  XJ2 = XJ2 + FRACT*XJTIMN(IC,JRH,KT,IQ,2,IB)
!                WRITE(6,30)IC,JRH,KT,IQ,IB,
!     &           10.**XJTIMN(IC,JRH,KT,IQ,1,IB), 
!     &           10.**XJTIMN(IC,JRH,KT,IQ,2,IB), FRACT
	        ENDDO
	      ENDDO
            ENDDO
	  ENDDO
	ENDDO

!
! Log10J -->J
        XJ1 = 10.**XJ1
        XJ2 = 10.**XJ2
! Inetrpolate to get J at inputed S
        X1 = S(1)
        Y1 = XJ1
        X2 = S(2)
        Y2 = XJ2
        IF(Y1.GT.Y2) THEN
          YJ = Y1*(Y2/Y1)**((X1-V)/(X1-X2))
        ELSE
          YJ=Y2
        ENDIF
        XJTIM = YJ
        RTIM = XR0
!      WRITE(6,40)X1,X2,V,Y1,Y2,YJ


!JBIM -- NH3=0
        XJ1 = 0.  
        XJ2 = 0.
        XR0 = 0.
!
        VOL = (dx1+dx2)*(dy1+dy2)*(dz1+dz2)*(du1+du2)
        VOL4 = (dx1+dx2)*(dy1+dy2)*(dz1+dz2)

        DO KT = KT1,KT2
          IF(KT.EQ.KT1) THEN
            dz = dz2
	  ELSE
            dz = dz1
          ENDIF
      	  DO JRH = JRH1,JRH2
            IF(JRH.EQ.JRH1) THEN
              dy = dy2
	    ELSE
              dy = dy1
            ENDIF
            DO IC = IC1,IC2
              IF(IC.EQ.IC1) THEN
                dx = dx2
	      ELSE
                dx = dx1
              ENDIF

                FRACT4 = dx*dy*dz/VOL4
                XR0 = XR0 + FRACT4*XRSTAR(IC,JRH,KT,1)

	        DO IQ =IQ1, IQ2
                  IF(IQ.EQ.IQ1) THEN
                    du = du2
	          ELSE
                    du = du1
                  ENDIF
                  FRACT = dx*dy*dz*du/VOL 
                  XJ1 = XJ1 + FRACT*XJTIMN(IC,JRH,KT,IQ,1,1)
                  XJ2 = XJ2 + FRACT*XJTIMN(IC,JRH,KT,IQ,2,1)
	        ENDDO
            ENDDO
	  ENDDO
	ENDDO

!
! Log10J -->J
        XJ1 = 10.**XJ1
        XJ2 = 10.**XJ2
! Inetrpolate to get J at inputed S
        X1 = S(1)
        Y1 = XJ1
        X2 = S(2)
        Y2 = XJ2
        IF(Y1.GT.Y2) THEN
          YJ = Y1*(Y2/Y1)**((X1-V)/(X1-X2))
        ELSE   ! When Y1 and Y2 both very small, Y1 can be LT Y2 becasue of cut-off or semi-equli approx, no extrapol
          YJ=Y2
        ENDIF
        XJBIM = YJ
        RBIM = XR0
!      WRITE(6,40)X1,X2,V,Y1,Y2,YJ
!         IF(XJBIM.GT.(1.1*XJTIM)) THEN
!           WRITE(6,40)X,Y,Z,U,V,W,XJBIM,XJTIM,Y1,Y2,(X1-V)/(X1-X2)
!         ENDIF

!JTHN  -- Q=0
        XJ1 = 0.  
        XJ2 = 0.
        XR0 = 0.
!
        VOL = (dx1+dx2)*(dy1+dy2)*(dz1+dz2)*(dw1+dw2)

        DO KT = KT1,KT2
          IF(KT.EQ.KT1) THEN
            dz = dz2
	  ELSE
            dz = dz1
          ENDIF
      	  DO JRH = JRH1,JRH2
            IF(JRH.EQ.JRH1) THEN
              dy = dy2
	    ELSE
              dy = dy1
            ENDIF
            DO IC = IC1,IC2
              IF(IC.EQ.IC1) THEN
                dx = dx2
	      ELSE
                dx = dx1
              ENDIF
 	      DO IB =IB1, IB2
                IF(IB.EQ.IB1) THEN
                  dw = dw2
	        ELSE
                  dw = dw1
                ENDIF
                FRACT = dx*dy*dz*dw/VOL 
                XR0 = XR0 + FRACT*XRSTAR(IC,JRH,KT,IB)
                XJ1 = XJ1 + FRACT*XJTIMN(IC,JRH,KT,1,1,IB)
                XJ2 = XJ2 + FRACT*XJTIMN(IC,JRH,KT,1,2,IB)
	      ENDDO
            ENDDO
	  ENDDO
	ENDDO
!
! Log10J -->J
        XJ1 = 10.**XJ1
        XJ2 = 10.**XJ2
! Inetrpolate to get J at inputed S
        X1 = S(1)
        Y1 = XJ1
        X2 = S(2)
        Y2 = XJ2
        IF(Y1.GT.Y2) THEN
          YJ = Y1*(Y2/Y1)**((X1-V)/(X1-X2))
        ELSE
          YJ=Y2
        ENDIF
        XJTH = YJ
        RTH = XR0

!JBH, Q=0, NH3=0
        XJ1 = 0.  
        XJ2 = 0.
        XR0 = 0.
        VOL = (dx1+dx2)*(dy1+dy2)*(dz1+dz2)
        DO KT = KT1,KT2
          IF(KT.EQ.KT1) THEN
            dz = dz2
	  ELSE
            dz = dz1
          ENDIF
      	  DO JRH = JRH1,JRH2
            IF(JRH.EQ.JRH1) THEN
              dy = dy2
	    ELSE
              dy = dy1
            ENDIF
            DO IC = IC1,IC2
              IF(IC.EQ.IC1) THEN
                dx = dx2
	      ELSE
                dx = dx1
              ENDIF
                FRACT = dx*dy*dz/VOL 
                XR0 = XR0 + FRACT*XRSTAR(IC,JRH,KT,1)
                XJ1 = XJ1 + FRACT*XJTIMN(IC,JRH,KT,1,1,1)
                XJ2 = XJ2 + FRACT*XJTIMN(IC,JRH,KT,1,2,1)
            ENDDO
	  ENDDO
	ENDDO
!
! Log10J -->J
        XJ1 = 10.**XJ1
        XJ2 = 10.**XJ2
! Inetrpolate to get J at inputed S
        X1 = S(1)
        Y1 = XJ1
        X2 = S(2)
        Y2 = XJ2
        IF(Y1.GT.Y2) THEN
          YJ = Y1*(Y2/Y1)**((X1-V)/(X1-X2))
        ELSE
          YJ=Y2
        ENDIF
        XJBH = YJ
        RBH = XR0
!
 30    FORMAT(I3, I3, I3, I3, I3, 10(1PE10.3))
 40    FORMAT(20(1PE10.3))

       END SUBROUTINE YUJTIMN1
!------------------------------------------------------------------------------

! *****************************************************************************
! IMNIMNIMNIMNIMNIMNIMNIMNIMNIMNIMNIMNIMNIMNIMNIMNIMNIMNIMNIMNIMNIMNIMNIMNIMNIM
! *****************************************************************************

        SUBROUTINE READJTIMN(DATA_DIR_1x1a)
!     
!       WRITTEN by Fangqun Yu, SUNY-Albany, 2006 (Revised, 8/2008)
!
! Read in the integrated IMN and KBHN lookup tables.
! IMN lookup table references: 
! 1. Yu, F., Ion-mediated nucleation in the atmosphere: Key controlling 
!      parameters, implications, and look-up table, J. Geophy. Res., 115, 
!      D03206, doi:10.1029/2009JD012630, 2010.
!
! KBHN lookup table reference:
! 2. Yu, F.,Updated H2SO4-H2O binary homogeneous nucleation rate look-up tables, 
!      J. Geophy. Res.,113, D24201, doi:10.1029/2008JD010527, 2008.
!
        CHARACTER(LEN=255)   :: DATA_DIR_1x1a
        INTEGER :: IC, IRH, IT, IQ, IS ,IB
        REAL*8  :: C11,Q11,S11

!       CHARACTER*2 YPATH
        CHARACTER*999 YPATH

        DATA_DIR_1x1= DATA_DIR_1x1a
        YPATH = TRIM(DATA_DIR_1x1)//'APM_data_201906/TIMN_TROP201706/'
!        YPATH = './'
        WRITE(6,*)"Read 6-D IMN look-up tables @ ", TRIM(YPATH)

        open(31,file=TRIM(YPATH)//'TIMN_J6D.txt',form='formatted')
        open(33,file=TRIM(YPATH)//'TIMN_Rstar4D.txt',form='formatted')

        open(41,file=TRIM(YPATH)//'TIMN_1H2SO4.txt',form='formatted')
        open(42,file=TRIM(YPATH)//'TIMN_2RH.txt',form='formatted')
        open(43,file=TRIM(YPATH)//'TIMN_3T.txt',form='formatted')
        open(44,file=TRIM(YPATH)//'TIMN_4Q.txt',form='formatted')
        open(45,file=TRIM(YPATH)//'TIMN_5S.txt',form='formatted')
        open(46,file=TRIM(YPATH)//'TIMN_6B.txt',form='formatted')
!
        READ(41,100)(C(IC),IC=1,MC)
!        WRITE(6,*)"[H2SO4](IC), IC=1, ", MC, ":"
!        WRITE(6,100)(C(IC),IC=1,MC)
!
        READ(42,100)(RH(IRH),IRH=1,MRH)
!        WRITE(6,*)"RH(IRH), IRH=1, ", MRH, ":"
!        WRITE(6,100)(RH(IRH),IRH=1,MRH)
!
        READ(43,100)(T(IT),IT=1,MT)
!        WRITE(6,*)"T(IT), IT=1, ", MT, ":"
!        WRITE(6,100)(T(IT),IT=1,MT)
!
        READ(44,100)(Q(IQ),IQ=1,MQ)
!        WRITE(6,*)"Q(I), I=1, ", MQ, ":"
!        WRITE(6,100)(Q(IQ),IQ=1,MQ)
!
        READ(45,100)(S(IS),IS=1,MS)
!        WRITE(6,*)"S(IS), IS=1, ", MS, ":"
!        WRITE(6,100)(S(IS),IS=1,MS)
!
        READ(46,100)(B(IB),IB=1,MB)

! Use the formula to calculate C and Q to get values with more digits, otherwise
! may cause problem when input C and Q are very clsoe to C(IC),Q(IQ) as
! IC and IQ are decided with formula 
!
        C(1) = 5.0E5
        C(MC) = 5.0E9
        DO IC = 2, MC-2
           C11 = C(IC)                                                          
           C(IC) = C(IC-1)*10.**(0.1)

           IF(abs(1.-C11/C(IC)).GT.0.02) THEN                                  
              write(6,*)"need check JTIMN look-up table inputs"                  
              stop                                                              
           ENDIF                                                                
        ENDDO

        DO IQ = 1, MQ
           Q11 = Q(IQ)                                                          
           IF(IQ.EQ.1) THEN
              Q(1) =1.E-30
           ELSEIF(IQ.EQ.MQ) THEN
              Q(MQ) =100.
           ELSE
              Q(IQ) = 2.0*1.5**float(IQ-2)
           ENDIF
           IF(abs(1.-Q11/Q(IQ)).GT.0.02) THEN
              write(6,*)"need check JTIMN look-up table inputs"
              stop
           ENDIF
        ENDDO

!        DO IS = 1, MS
!           S11 = S(IS)                                                          
!           IF(IS.EQ.1) THEN
!              S(1) =1.0
!           ELSE
!              S(IS) = 10.*100.**(0.1*float(IS-2))
!           ENDIF
!           IF(abs(1.-S11/S(IS)).GT.0.02) THEN
!              write(6,*)"need check JTIMN look-up table inputs"
!              stop
!           ENDIF
!        ENDDO

!
! READ in formatted 6-D Table
! Due to high sensitivity of J to key parameters, use logJ to interpolate
!
        DO IS =1, MS
         DO IT = 1,MT
         DO IRH = 1,MRH
         DO IQ =1, MQ
         DO IB =1, MB 
          READ(31,201)(XJTIMN(IC,IRH,IT,IQ,IS,IB),IC = 1,MC)
          DO IC=1, MC
           XJTIMN(IC,IRH,IT,IQ,IS,IB)=LOG10(XJTIMN(IC,IRH,IT,IQ,IS,IB))
          ENDDO
         ENDDO
         ENDDO
         ENDDO
         ENDDO
        ENDDO
! Critical cluster properties depend on T, RH, [H2SO4], NH3 only
        DO IT = 1,MT
         DO IRH = 1, MRH
         DO IB =1, MB 
           READ(33,203)(XRSTAR(IC,IRH,IT,IB),IC=1,MC)
         ENDDO  ! B
         ENDDO  ! RH
        ENDDO   !T

        CLOSE(31)
        CLOSE(33)
        CLOSE(41)
        CLOSE(42)
        CLOSE(43)
        CLOSE(44)
        CLOSE(45)
        CLOSE(46)
!
 100    FORMAT(100(1PE10.3))
 200    FORMAT(100(1PE9.2))
 201    FORMAT(100(1PE9.2))
 202    FORMAT(100F5.1)
 203    FORMAT(100F5.2)
 204    FORMAT(100F6.3)

      END SUBROUTINE READJTIMN
! *****************************************************************************
! IMNIMNIMNIMNIMNIMNIMNIMNIMNIMNIMNIMNIMNIMNIMNIMNIMNIMNIMNIMNIMNIMNIMNIMNIMNIM
! *****************************************************************************
!
!------------------------------------------------------------------------------

      ! End of module
      END MODULE APM_TIMN_MOD
#endif
