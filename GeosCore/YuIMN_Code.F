#if   defined ( TOMAS )
C IMNIMNIMNIMNIMNIMNIMNIMNIMNIMNIMNIMNIMNIMNIMNIMNIMNIMNIMNIMNIMNIMNIMN
C test program
C         PROGRAM TEST
C         REAL*8 :: X, Y, Z, U, V, XJ0, XI0,XR0,XAMOLF0
C         CALL READJIMN5D   ! The table only needs to be read once at the beginning
C         DO ITRY = 1, 100
C             WRITE(6,*)"Input T(K)=, RH(%)=, CH2SO4(#/cm3), Q=, S="
C             READ(5,*)Z, Y, X, U, V
C             CALL YUJIMN(X,Y,Z,U, V, XJ0,XI0,XR0,XAMOLF0)
C             XIWATER = XI0 * (1.-XAMOLF0)/XAMOLF0
C             WRITE(6,100)Z, Y, X, U, V,XJ0,2.*XR0,XI0,XIWATER
C         ENDDO
C  100    FORMAT(10(1PE10.3))
C
C         STOP
C         END
C
C *********************************************************************
C IMNIMNIMNIMNIMNIMNIMNIMNIMNIMNIMNIMNIMNIMNIMNIMNIMNIMNIMNIMNIMNIMNIMN
C *********************************************************************
C
         SUBROUTINE YUJIMN(X0,Y0,Z0,U0,V0,XJ0,XI0,XR0,XAMOLF0)
C
C       This subroutine is to calculate rates and critical cluster
C       properties of ion-mediated nucleation (IMN) from lookup tables
C       using multiple-variable interpolation scheme
C
C        WRITTEN by Fangqun Yu, SUNY-Albany, 2006; UPDATED 06/2009
C        Email: yfq@asrc.cestm.albany.edu
C
C References:
C 1. Yu, F., Ion-mediated nucleation in the atmosphere: Key controlling
C       parameters, implications, and look-up table, J. Geophy. Res.,
C       ###, 2009.
C 2. Yu, F., From molecular clusters to nanoparticles: Second-generation
C       ion-mediated nucleation model, Atmos. Chem. Phys., 6, 5193-5211, 2006.
C 3. Yu, F., and R. P. Turco, Ultrafine aerosol formation via ion-mediated
C       nucleation, Geophys. Res. Lett., 27, 883-886, 2000.
C
C Input:
C X0 = [H2SO4] in #/cm3 (5E5-5E8)
C Y0 = RH in % (0.5-99.5)
C Z0 = T (in K) (190-302)
C U0 = Q = ionization rate (ion-pairs/cm3s) (0, 1.5-60)
C S0 = S = surface area (um2/cm3) (1-1000)
C
C OUTPUT:
C XJ0: Nucleation rate (#/cm3s)
C XI0: Number of H2SO4 molecules in critical cluster
C XR0: Radius of critical cluster (nm)
C XAMOLF0: Critical cluster H2SO4 mole fraction
C
         PARAMETER (MC=31, MRH=51,MT=57, MQ=9, MS=7)
         COMMON /YYJIMNHT/C(MC),RH(MRH),T(MT),Q(MQ),S(MS),
     &                 XJIMN(MC,MRH,MT,MQ,MS),
     &                 XISTAR(MC,MRH,MT),XRSTAR(MC,MRH,MT),
     &                 XAMOLFSTAR(MC,MRH,MT)
C
        REAL*8 X0,Y0,Z0,U0,V0,XJ0,XI0,XR0,XAMOLF0
C
C to avoid the input values to be changed due to out of the range reset
C
        X = X0
        Y = Y0
        Z = Z0
        U = U0
        V = V0
C
C The lookup table should cover almost all possible conditions in
C ambient troposphere. For the extreme conditions that are out of
C the ranges of the lookup table, we either reset the inputed
C parameters in a way that may underestimate the JIMN values or
C set the nucleation rate to 1.E-20 cm-3s-1.
C
        IF(X.LT.C(1)) THEN
C           WRITE(6,10) X, C(1)
           XJ0 = 1.E-20
           XI0 = 100.
           XR0 = 1.5
           XAMOLF0 = 0.5
           RETURN
        ELSEIF(X.GT.C(MC)) THEN
c           WRITE(6,11) X, C(MC), C(MC)
           X =C(MC)
        ENDIF
        IF(Y.LT.RH(1)) THEN
c           WRITE(6,12) Y, RH(1)
           XJ0 = 1.E-20
           XI0 = 100.
           XR0 = 1.5
           XAMOLF0 = 0.5
           RETURN
        ELSEIF(Y.GT.RH(MRH)) THEN
c           WRITE(6,13) Y, RH(MRH), RH(MRH)
           Y =RH(MRH)
        ENDIF
        IF(Z.LT.T(1)) THEN
c           WRITE(6,14) Z, T(1), T(1)
           Z =T(1)
        ELSEIF(Z.GT.T(MT)) THEN
c           WRITE(6,15) Z, T(MT)
           XJ0 = 1.E-20
           XI0 = 100.
           XR0 = 1.5
           XAMOLF0 = 0.5
           RETURN
        ENDIF
        IF(U.LT.Q(1)) THEN
c           WRITE(6,16) U, Q(1)
           XJ0 = 1.E-20
           XI0 = 100.
           XR0 = 1.5
           XAMOLF0 = 0.5
           RETURN
        ELSEIF(U.GT.Q(MQ)) THEN
c           WRITE(6,17) U, Q(MQ), Q(MQ)
           U =Q(MQ)
       ENDIF
       IF(V.LT.S(1)) THEN
C          WRITE(6,18) V, S(1), S(1)
          V =S(1)
       ELSEIF(V.GT.S(MS)) THEN
c          WRITE(6,19) V, S(MS)
          XJ0 = 1.E-20
          XI0 = 100.
          XR0 = 1.5
          XAMOLF0 = 0.5
          RETURN
       ENDIF
  10   FORMAT("IMN WARNING: INPUTED [H2SO4]=",E9.2,"<",E9.2,
     &    ", set JIMN to 1.E-20 cm-3s-1")
  11   FORMAT("IMN WARNING: INPUTED [H2SO4]=",E9.2,">",E9.2,
     &    " set it to ",E9.2)
  12   FORMAT("IMN WARNING: INPUTED RH =",F5.1,"% <",F5.1,
     &    "%, set JIMN to 1.E-20 cm-3s-1")
  13   FORMAT("IMN WARNING: INPUTED RH =",F5.1,"% >",F5.1,
     &    "% set it to ",F5.1,"%")
  14   FORMAT("IMN WARNING: INPUTED T =",F6.1,"K <",F6.1,
     &    "K set it to ",F6.1,"K")
  15   FORMAT("IMN WARNING: INPUTED T =",F6.1,"K >",F6.1,
     &    "K, set JIMN to 1.E-20 cm-3s-1")
  16   FORMAT("IMN WARNING: INPUTED Q =",F6.1," <",F6.1,
     &    " ion-pair/cm3s , set JIMN to 1.E-20 cm-3s-1")
  17   FORMAT("IMN WARNING: INPUTED Q =",F6.1," >",F6.1,
     &    " ion-pair/cm3s set it to ",F6.1)
  18   FORMAT("IMN WARNING: INPUTED S =",F6.1," <",F6.1,
     &    " um2/cm3 set it to ",F6.1)
  19   FORMAT("IMN WARNING: INPUTED S =",F6.1," >",F6.1,
     &    "um2/cm3, set JIMN to 1.E-20 cm-3s-1")
       IC1 =MAX0(INT(1.+10.*ALOG10(X/C(1))),1)
       IC2 = MIN0(IC1 + 1,MC)
       IF(IC2.EQ.MC) IC1=MC-1
       IF(Y.LT.RH(2)) THEN
          JRH1 = 1.
       ELSE
        JRH1 = MAX0(INT((Y-RH(2))/2.+2.),2)
       ENDIF
       JRH2 = MIN0(JRH1 + 1,MRH)
       IF(JRH2.EQ.MRH) JRH1=MRH-1
       KT1 = MAX0(INT(Z/2.-94.0),1)
       KT2 = MIN0(KT1 + 1,MT)
       IF(KT2.EQ.MT) KT1=MT-1
C
       IQ1 = MAX0(INT(1.+5.*ALOG10(U/Q(1))),1)
       IQ2 = MIN0(IQ1 + 1,MQ)
       IF(IQ2.EQ.MQ) IQ1=MQ-1
C
       IF(V.LT.10.0) THEN
         IS1 =1.
       ELSE
         IS1 = MAX0(INT(2.+2.5*ALOG10(V/10.)),2)
       ENDIF
       IS2 = MIN0(IS1 + 1,MS)
       IF(IS2.EQ.MS) IS1=MS-1
C
      dx1 =  ALOG10(X/C(IC1))    ! logJ log[H2SO4] interpolation
      dx2 =  ALOG10(C(IC2)/X)
      dy1 =  ALOG10(Y/RH(JRH1))
      dy2 =  ALOG10(RH(JRH2)/Y)
      dz1 =  Z-T(KT1)
      dz2 =  T(KT2)-Z
      du1 =  U - Q(IQ1)
      du2 =  Q(IQ2) - U
      dv1 =  V- S(IS1)
      dv2 =  S(IS2) - V
C
      XJ0 = 0.
      XI0 = 0.
      XR0 = 0.
      XAMOLF0 = 0.
C
      VOL = (dx1+dx2)*(dy1+dy2)*(dz1+dz2)*(du1+du2)*(dv1+dv2)
      VOL3 = (dx1+dx2)*(dy1+dy2)*(dz1+dz2)
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
               FRACT3 = dx*dy*dz/VOL3
               XI0 = XI0 + FRACT3*XISTAR(IC,JRH,KT)
               XR0 = XR0 + FRACT3*XRSTAR(IC,JRH,KT)
               XAMOLF0 = XAMOLF0 + FRACT3*XAMOLFSTAR(IC,JRH,KT)
               DO IS =IS1, IS2
                  IF(IS.EQ.IS1) THEN
                     dv = dv2
                  ELSE
                     dv = dv1
                  ENDIF
                  DO IQ =IQ1, IQ2
                     IF(IQ.EQ.IQ1) THEN
                        du = du2
                     ELSE
                        du = du1
                     ENDIF
                     FRACT = dx*dy*dz*du*dv/VOL
                     XJ0 = XJ0 + FRACT*XJIMN(IC,JRH,KT,IQ,IS)
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO
C
C Log10J -->J
         XJ0 = 10.**XJ0
C
  30   FORMAT(I3, I3, I3, I3, I3, 10(1PE10.3))
  20   FORMAT(10(1PE10.3))
        RETURN
        END
C *********************************************************************
        SUBROUTINE READJIMN5D( am_I_Root, Input_Opt, RC )
C
C       WRITTEN by Fangqun Yu, SUNY-Albany, 2006 (Updated, 6/2009)
C
!
! !USES:
!
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)    :: am_I_Root   ! Are we on the root CPU?
      TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!
! !DEFINED PARAMETERS:
!
C Parameters
C (1 ) MC     : NUMBER OF POINTS IN H2SO4 CONCENTRATION DIMENSION
C (2 ) MT     : NUMBER OF POINTS IN TEMPERATURE DIMENSION
C (3 ) MRH : NUMBER OF POINTS IN RELATIVE HUMIDITY DIMENSION
C (4 ) MQ     : NUMBER OF POINTS IN IONIZATION RATE DIMENSION
C (5 ) MS     : NUMBER OF POINTS IN SURFACE AREA DIMENSION
C Arrays
C (6 ) C    : VALUES AT POINTS IN H2SO4 CONCENTRATION DIMENSION
C (7 ) T    : VALUES AT POINTS IN TEMPERATURE DIMENSION
C (8 ) RH : VALUES AT POINTS IN RELATIVE HUMIDITY DIMENSION
C (9 ) Q    : VALUES AT POINTS IN IONIZATION RATE DIMENSION
C (10) S    : VALUES AT POINTS IN SURFACE AREA DIMENSION
        PARAMETER (MC=31, MRH=51,MT=57, MQ=9, MS=7)
        COMMON /YYJIMNHT/C(MC),RH(MRH),T(MT),Q(MQ),S(MS),
     &                 XJIMN(MC,MRH,MT,MQ,MS),
     &                 XISTAR(MC,MRH,MT),XRSTAR(MC,MRH,MT),
     &                 XAMOLFSTAR(MC,MRH,MT)
!
! !REVISION HISTORY
!  30 Jan 2014 - R. Yantosca - Now read big files from shared data space
!
!----------------------------------------------------------
! Prior to 1/30/14:
!        CHARACTER*6 YPATH
!        YPATH = './'
!----------------------------------------------------------
        CHARACTER(LEN=255) :: DATA_DIR 
        CHARACTER(LEN=255) :: FNAME
        
        ! Assume success
        RC       = GIGC_SUCCESS

        ! Data directory path in shared disk space where files live
        DATA_DIR = TRIM( Input_Opt%DATA_DIR_1x1 ) // 'TOMAS_201402/'

        WRITE(6,*) "Read IMN look-up tables"

        FNAME = TRIM( DATA_DIR ) // 'YuIMN_J5D.txt'
        open(31, file=TRIM( FNAME ), form='formatted')

        FNAME = TRIM( DATA_DIR ) // 'YuIMN_Istar3D.txt'
        open(32, file=TRIM( FNAME ), form='formatted')

        FNAME = TRIM( DATA_DIR ) // 'YuIMN_Rstar3D.txt'
        open(33, file=TRIM( FNAME ), form='formatted')

        FNAME = TRIM( DATA_DIR ) // 'YuIMN_AMOLF3D.txt'
        open(34, file=TRIM( FNAME ), form='formatted')

        FNAME = TRIM( DATA_DIR ) // 'YuIMN_1H2SO4.txt'
        open(41, file=TRIM( FNAME ), form='formatted')

        FNAME = TRIM( DATA_DIR ) // 'YuIMN_2RH.txt'
        open(42, file=TRIM( FNAME ), form='formatted')

        FNAME = TRIM( DATA_DIR ) // 'YuIMN_3T.txt'
        open(43, file=TRIM( FNAME ), form='formatted')

        FNAME = TRIM( DATA_DIR ) // 'YuIMN_4Q.txt'
        open(44, file=TRIM( FNAME ), form='formatted')

        FNAME = TRIM( DATA_DIR ) // 'YuIMN_5S.txt'
        open(45, file=TRIM( FNAME ), form='formatted')
C
        READ(41,101)(C(I),I=1,MC)
        WRITE(6,*)"[H2SO4](I), I=1, ", MC, ":"
        WRITE(6,100)(C(I),I=1,MC)
C
        READ(42,102)(RH(J),J=1,MRH)
        WRITE(6,*)"RH(I), I=1, ", MRH, ":"
        WRITE(6,100)(RH(J),J=1,MRH)
C
        READ(43,103)(T(IT),IT=1,MT)
        WRITE(6,*)"T(I), I=1, ", MT, ":"
        WRITE(6,100)(T(IT),IT=1,MT)
C
        READ(44,104)(Q(IQ),IQ=1,MQ)
        WRITE(6,*)"Q(I), I=1, ", MQ, ":"
        WRITE(6,100)(Q(IQ),IQ=1,MQ)
C
        READ(45,105)(S(IS),IS=1,MS)
        WRITE(6,*)"S(I), I=1, ", MS, ":"
        WRITE(6,100)(S(IS),IS=1,MS)
C
C Use the formula to calculate C and Q to get values with more digits, otherwise
C may cause problem when input C and Q are very close to C(IC),Q(IQ) as
C IC and IQ are decided with formula
C 
       C(1) = 5.0E5
        DO IC = 2, MC
           C11 = C(IC)
           RATIO = 10.**(0.1)
           C(IC) = C(IC-1)*RATIO
           IF(abs(1.-C11/C(IC)).GT.0.02) THEN
              write(6,*)"need check JIMN look-up table inputs C"
              stop
           ENDIF
        ENDDO
        DO IQ = 1, MQ
           Q11 = Q(IQ)
           Q(IQ) = 1.5*10.**(0.2*float(IQ-1))
           IF(abs(1.-Q11/Q(IQ)).GT.0.02) THEN
              write(6,*)"need check JIMN look-up table inputs Q"
              stop
           ENDIF
        ENDDO
        DO IS = 1, MS
           S11 = S(IS)
           IF(IS.EQ.1) THEN
              S(1) =1.0
           ELSE
              S(IS) = 10.*100.**(0.2*float(IS-2))
           ENDIF
           IF(abs(1.-S11/S(IS)).GT.0.02) THEN
              write(6,*)"need check JIMN look-up table inputs S"
              stop
           ENDIF
        ENDDO
C
C Formatted 5-D Table
C
        DO IS =1, MS
          DO KT = 1,MT
            DO JRH = 1,MRH
                DO IQ =1, MQ
                  READ(31,101)(XJIMN(IC,JRH,KT,IQ,IS),IC = 1,MC)
                  DO IC=1, MC
C                   IF(XJIMN(IC,JRH,KT,IQ,IS).LT.1.E-20)
C     &                  XJIMN(IC,JRH,KT,IQ,IS)=1.E-20
C Due to high sensitivity of J to key parameters, use logJ to interpolate
                   XJIMN(IC,JRH,KT,IQ,IS)=ALOG10(XJIMN(IC,JRH,KT,IQ,IS))
                  ENDDO
                ENDDO
            ENDDO
          ENDDO
        ENDDO
C Critical cluster properties depend on T, RH, [H2SO4] only
        DO IT = 1,MT
          DO IRH = 1, MRH
            READ(32,202)(XISTAR(IC,IRH,IT),
     &                                 IC=1,MC)
            READ(33,203)(XRSTAR(IC,IRH,IT),
     &                                 IC=1,MC)
            READ(34,204)(XAMOLFSTAR(IC,IRH,IT),
     &                                 IC=1,MC)
          ENDDO ! RH
        ENDDO   !T
        CLOSE(31)
        CLOSE(32)
        CLOSE(33)
        CLOSE(34)
        CLOSE(41)
        CLOSE(42)
        CLOSE(43)
        CLOSE(44)
        CLOSE(45)
C
  100   FORMAT(100E9.2)
  101   FORMAT(31E9.2) ! H2SO4
  102   FORMAT(51E9.2) ! RH
  103   FORMAT(57E9.2) ! T
  104   FORMAT(9E9.2)  ! Q
  105   FORMAT(7E9.2)  ! S
c  100   FORMAT(100(1PE9.2))
c  101   FORMAT(31(1PE9.2)) ! H2SO4
c  102   FORMAT(51(1PE9.2)) ! RH
c  103   FORMAT(57(1PE9.2)) ! T
c  104   FORMAT(9(1PE9.2))  ! Q
c  105   FORMAT(7(1PE9.2))  ! S
  202   FORMAT(31F5.1) ! Istar
  203   FORMAT(31F5.2) ! Rstar
  204   FORMAT(31F6.3) ! AMOLF

      print*,'read Yu inputs'
      print*,C
        RETURN
        END
C *********************************************************************
C IMNIMNIMNIMNIMNIMNIMNIMNIMNIMNIMNIMNIMNIMNIMNIMNIMNIMNIMNIMNIMNIMNIMN
C *********************************************************************
C
#endif
