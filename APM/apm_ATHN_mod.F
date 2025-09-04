#ifdef APM
      MODULE APM_ATHN_MOD
      IMPLICIT NONE

      !=================================================================
      ! MODULE PRIVATE DECLARATIONS -- keep certain internal variables
      ! and routines from being seen outside "apm_nucl_mod.f"
      !=================================================================

!  Module Variables:
!  ============================================================================
!  Parameters
!  (1 ) MC   : NUMBER OF POINTS IN H2SO4 CONCENTRATION DIMENSION
!  (2 ) MRH  : NUMBER OF POINTS IN RELATIVE HUMIDITY DIMENSION
!  (3 ) MT   : NUMBER OF POINTS IN TEMPERATURE DIMENSION
!  (4 ) MD   : NUMBER OF POINTS IN [DMA] DIMENSION
!  (5 ) MS   : NUMBER OF POINTS IN SURFACE AREA DIMENSION

!  Arrays
!  (6 ) C   : VALUES AT POINTS IN H2SO4 CONCENTRATION DIMENSION
!  (7 ) RH  : VALUES AT POINTS IN RELATIVE HUMIDITY DIMENSION
!  (8 ) T   : VALUES AT POINTS IN TEMPERATURE DIMENSION
!  (9 ) D   : VALUES AT POINTS IN [DMA] DIMENSION
!  (10) S   : VALUES AT POINTS IN SURFACE AREA DIMENSION

!  (11) XJATHN : AMINE TERNARY NUCLEATION RATES (cm-3s-1) AT ALL POINTS IN 5-D SPACE

!  Module Routines:
!  ============================================================================
!  (1 ) YUJATHN     : INTERPOLAION SCHEME TO FIND JATHN FROM LOOKUP TABLE
!  (2 ) READJATHN   : READ IN THE ATHN LOOKUP TABLE
!
!  NOTES:
!  (1 ) ....
!******************************************************************************


      ! Make everything PRIVATE ...
      PRIVATE

      ! ... except these variables ...
!      PUBLIC ::

      ! ... and these routines
      PUBLIC :: YUJATHN
      PUBLIC :: READJATHN

      !=================================================================
      ! MODULE VARIABLES
      !=================================================================
      ! Parameters
      INTEGER, PARAMETER   :: MC  = 31
      INTEGER, PARAMETER   :: MRH = 11
      INTEGER, PARAMETER   :: MT  = 31
      INTEGER, PARAMETER   :: MD  = 31
      INTEGER, PARAMETER   :: MS  = 12

      ! Arrays
      REAL*8               :: C(MC),RH(MRH),T(MT),D(MD),S(MS)
      REAL*8               :: XJATHN(MC,MRH,MT,MD,MS)

      CHARACTER(LEN=255)   :: DATA_DIR_1x1
      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement
      !=================================================================
      CONTAINS


! *********************************************************************
! ATHNATHNATHNATHNATHNATHNATHNATHNATHNATHNATHNATHNATHNATHNATHNATHNATHN
! *********************************************************************
!
	SUBROUTINE YUJATHN(X0,Y0,Z0,U0,V0,XJ0)
!
!      This subroutine is to calculate rates of amine ternary homogeneous 
!      nucleation (ATHN) from lookup tables
!      using multiple-variable interpolation scheme
!
!       WRITTEN by Fangqun Yu, SUNY-Albany, 2014
!       Email: fyu@albany.edu
!       Updated 2015
! References: 
!
! Input:
! X0 = [H2SO4] in #/cm3  (5E5-5E8)
! Y0 = RH in % (0.5-99.5)
! Z0 = T (in K) (250-310)
! U0 = [DMA] (ppt) (0.1-100)
! V0 = S (um2/cm3) (1-1000)
!
! OUTPUT:
! XJ0: Nucleation rate (#/cm3s)
!
!
        REAL*8 :: X0,Y0,Z0,U0,V0,XJ0
        REAL*8  :: X,Y,Z,U,V
        REAL*8  :: VOL,FRACT
        REAL*8  :: dx1,dx2,dy1,dy2,dz1,dz2,du1,du2,dv1,dv2
        REAL*8  :: dx,dy,dz,du,dv

        INTEGER :: IC1, IC2, JRH1, JRH2, KT1, KT2, ID1, ID2, IS1,IS2
        INTEGER :: IC, JRH, KT, ID,IS
!
! to avoid the input values to be changed due to out of the range reset
!
        X = X0
        Y = Y0
        Z = Z0
        U = U0
        V = V0
!
! The lookup table should cover almost all possible conditions in
! ambient troposphere. For the extreme conditions that are out of 
! the ranges of the lookup table, we either reset the inputed 
! parameters in a way that may underestmate the JATHN values or
! set the nucleation rate to 1.E-20 cm-3s-1.
!
        IF(X.LT.C(1)) THEN
!           WRITE(6,10) X, C(1)
           XJ0 = 1.E-20
           RETURN
        ELSEIF(X.GT.C(MC)) THEN
!           WRITE(6,11) X, C(MC), C(MC)
           X =C(MC)
        ENDIF

        IF(Y.LT.RH(1)) THEN
!           WRITE(6,12) Y, RH(1)
           XJ0 = 1.E-20
           RETURN
        ELSEIF(Y.GT.RH(MRH)) THEN
!           WRITE(6,13) Y, RH(MRH), RH(MRH)
           Y =RH(MRH)
        ENDIF

        IF(Z.LT.T(1)) THEN
!           WRITE(6,14) Z, T(1), T(1)
           Z =T(1)
        ELSEIF(Z.GT.T(MT)) THEN
!           WRITE(6,15) Z, T(MT)
           XJ0 = 1.E-20
           RETURN
        ENDIF

        IF(U.LT.D(1)) THEN
           WRITE(6,16) U, D(1)
           XJ0 = 1.E-20
           RETURN
        ELSEIF(U.GT.D(MD)) THEN
!           WRITE(6,17) U, D(MD)
           U =D(MD)
        ENDIF

        IF(V.LT.S(1)) THEN 
!           WRITE(86,18) V, S(1), S(1)
           V =S(1)
        ELSEIF(V.GT.S(MS)) THEN
!           WRITE(86,19) V, S(MS), S(MS)
           V =S(MS)
        ENDIF



 10     FORMAT("ATHN WARNING: INPUTED [H2SO4]=",ES9.2,"<",ES9.2,
     &     ", set JATHN to 1.E-20 cm-3s-1")
 11     FORMAT("ATHN WARNING: INPUTED [H2SO4]=",ES9.2,">",ES9.2,
     &     " set it to ",ES9.2)
 12     FORMAT("ATHN WARNING: INPUTED RH =",F5.1,"% <",F5.1,
     &     "%, set JATHN to 1.E-20 cm-3s-1")
 13     FORMAT("ATHN WARNING: INPUTED RH =",F5.1,"% >",F5.1,
     &     "% set it to ",F5.1,"%")
 14     FORMAT("ATHN WARNING: INPUTED T =",F6.1,"K <",F6.1,
     &     "K set it to ",F6.1,"K")
 15     FORMAT("ATHN WARNING: INPUTED T =",F6.1,"K >",F6.1,
     &     "K, set JATHN to 1.E-20 cm-3s-1")
 16     FORMAT("ATHN WARNING: INPUTED D =",F6.1," <",F6.1,
     &     " ppt , set JATHN to 1.E-20 cm-3s-1")
 17     FORMAT("ATHN WARNING: INPUTED D =",F6.1," >",F6.1,
     &     " ppt set it to ",F6.1)
 18     FORMAT("ATHN WARNING: INPUTED S =",F6.1," <",F6.1,
     &     " um2/cm3 set it to ",F6.1)
 19     FORMAT("ATHN WARNING: INPUTED S =",F6.1," >",F6.1,
     &     " um2/cm3 set it to ",F6.1)


        IC1 =MAX0(INT(1.+10.*LOG10(X/C(1))),1)
        IC2 = MIN0(IC1 + 1,MC)
        IF(IC2.EQ.MC) IC1=MC-1
        
        IF(Y.LT.RH(2)) THEN
           JRH1 = 1.
        ELSE
         JRH1 = MAX0(INT((Y-RH(2))/10.+2.),2)
        ENDIF
        JRH2 = MIN0(JRH1 + 1,MRH)
        IF(JRH2.EQ.MRH) JRH1=MRH-1

        KT1 = MAX0(INT(Z/2.-124.0),1)
        KT2 = MIN0(KT1 + 1,MT)
        IF(KT2.EQ.MT) KT1=MT-1
!
        ID1 = MAX0(INT(1.+10.*LOG10(U/D(1))),1)
        ID2 = MIN0(ID1 + 1,MD)
        IF(ID2.EQ.MD) ID1=MD-1

!
        IF(V.LT.10.0) THEN
          IS1 =1.
        ELSE
          IS1 = MAX0(INT(2.+5.*LOG10(V/10.)),2)
        ENDIF
        IS2 = MIN0(IS1 + 1,MS)
        IF(IS2.EQ.MS) IS1=MS-1

!
	dx1 = LOG10(X/C(IC1))   ! logJ log[H2SO4] interpolation
	dx2 = LOG10(C(IC2)/X)
	dy1 = LOG10(Y/RH(JRH1))
	dy2 = LOG10(RH(JRH2)/Y)
	dz1 = Z-T(KT1)
	dz2 = T(KT2)-Z

        du1 = LOG10(U/D(ID1))
        du2 = LOG10(D(ID2)/U)
        dv1 = LOG10(V/S(IS1))
        dv2 = LOG10(S(IS2)/V)
!
        XJ0 = 0.
!
        VOL = (dx1+dx2)*(dy1+dy2)*(dz1+dz2)*(du1+du2)*(dv1+dv2)
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

	      DO ID =ID1, ID2
                IF(ID.EQ.ID1) THEN
                  du = du2
	        ELSE
                  du = du1
                ENDIF
                DO IS =IS1, IS2
                  IF(IS.EQ.IS1) THEN
                    dv = dv2
                  ELSE
                    dv = dv1
                  ENDIF

                  FRACT = dx*dy*dz*du*dv/VOL 
                  XJ0 = XJ0 + FRACT*XJATHN(IC,JRH,KT,ID,IS)
!                WRITE(6,30)IC,JRH,KT,ID,IS,10.**XJATHN(IC,JRH,KT,ID,IS),
!     &                    FRACT

	        ENDDO
	      ENDDO
            ENDDO
	  ENDDO
	ENDDO
!
! Log10J -->J
         XJ0 = 10.**XJ0
!
 30    FORMAT(I3, I3, I3, I3, I3, 10(1PE10.3))
 20    FORMAT(10(1PE10.3))
        RETURN 
        END

! *********************************************************************
        SUBROUTINE READJATHN(DATA_DIR_1x1a)
!     
!       WRITTEN by Fangqun Yu, SUNY-Albany, 2014, 2015
!
!  ============================================================================

        CHARACTER(LEN=255)   :: DATA_DIR_1x1a
        CHARACTER*999 YPATH
        INTEGER :: IC, IRH, IT, ID, IS
        REAL*8  :: C11,D11,S11,ratio
        
        DATA_DIR_1x1= DATA_DIR_1x1a
        YPATH = TRIM(DATA_DIR_1x1)//'/APM_data_201906/ATHN20150328/'

        WRITE(6,*)"Read ATHN look-up tables"

        open(31,file=TRIM(YPATH)//'YuATHN_J5D.txt',form='formatted')
        open(41,file=TRIM(YPATH)//'YuATHN_1H2SO4.txt',form='formatted')
        open(42,file=TRIM(YPATH)//'YuATHN_2RH.txt',form='formatted')
        open(43,file=TRIM(YPATH)//'YuATHN_3T.txt',form='formatted')
        open(44,file=TRIM(YPATH)//'YuATHN_4DMA.txt',form='formatted')
        open(45,file=TRIM(YPATH)//'YuATHN_5S.txt',form='formatted')

!
        READ(41,100)(C(IC),IC=1,MC)
        WRITE(6,*)"[H2SO4](IC), IC=1, ", MC, ":"
        WRITE(6,100)(C(IC),IC=1,MC)
!
        READ(42,100)(RH(IRH),IRH=1,MRH)
        WRITE(6,*)"RH(I), I=1, ", MRH, ":"
        WRITE(6,100)(RH(IRH),IRH=1,MRH)
!
        READ(43,100)(T(IT),IT=1,MT)
        WRITE(6,*)"T(I), I=1, ", MT, ":"
        WRITE(6,100)(T(IT),IT=1,MT)
!
        READ(44,100)(D(ID),ID=1,MD)
        WRITE(6,*)"D(I), I=1, ", MD, ":"
        WRITE(6,100)(D(ID),ID=1,MD)
!
        READ(45,100)(S(IS),IS=1,MS)
        WRITE(6,*)"S(I), I=1, ", MS, ":"
        WRITE(6,100)(S(IS),IS=1,MS)

!
! Use the formula to calculate C and D to get values with more digits, otherwise
! may cause problem when input C and D are very clsoe to C(IC),D(ID) as
! IC and ID are decided with formula
!
        C(1) = 5.0E5
        DO IC = 2, MC
           C11 = C(IC)                                                          
           RATIO = 10.**(0.1)
           C(IC) = C(IC-1)*RATIO

           IF(abs(1.-C11/C(IC)).GT.0.02) THEN                                  
              write(6,*)"need check JATHN look-up table inputs"                  
              stop                                                              
           ENDIF                                                                
        ENDDO

        DO ID = 1, MD
           D11 = D(ID)                                                          
           D(ID) = 0.1*10.**(0.1*float(ID-1))
           IF(abs(1.-D11/D(ID)).GT.0.02) THEN
              write(6,*)"need check JATHN look-up table inputs"
              stop
           ENDIF
        ENDDO

        DO IS = 1, MS
           S11 = S(IS)                                                          
           IF(IS.EQ.1) THEN
              S(1) =1.0 
           ELSE
              S(IS) = 10.*100.**(0.1*float(IS-2))
           ENDIF
           IF(abs(1.-S11/S(IS)).GT.0.02) THEN
              write(6,*)"need check JATHN look-up table inputs"
              stop
           ENDIF
        ENDDO 

!
! Formatted 5-D Table
!
        DO IT = 1,MT
          DO IRH = 1,MRH
            DO IC =1, MC
              DO IS =1, MS
 	        READ(31,201)(XJATHN(IC,IRH,IT,ID,IS),ID =1,MD)
                DO ID=1, MD
! Due to high sensitivity of J to key parameters, use logJ to interpolate

                  XJATHN(IC,IRH,IT,ID,IS)=LOG10(XJATHN(IC,IRH,IT,ID,IS))
                ENDDO
	      ENDDO
	    ENDDO
          ENDDO
	ENDDO

        CLOSE(31)
        CLOSE(41)
        CLOSE(42)
        CLOSE(43)
        CLOSE(44)
        CLOSE(45)
!
 100    FORMAT(100(1PE10.2))
 201    FORMAT(100(1PE9.2))
 202    FORMAT(100F5.1)
 203    FORMAT(100F5.2)
 204    FORMAT(100F6.3)
 
        RETURN
        END      
! *********************************************************************
! ATHNATHNATHNATHNATHNATHNATHNATHNATHNATHNATHNATHNATHNATHNATHNATHNATHN
! *********************************************************************
!
      ! End of module
      END MODULE APM_ATHN_MOD
#endif
