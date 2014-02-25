!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: apm_nucl_mod
!
! !DESCRIPTION: Module APM\_NUCL\_MOD contains variables and routines for 
!  computing nucleation rates and ionization rates.
!\\
!\\
! !INTERFACE:
!
      MODULE APM_NUCL_MOD
!
! !USES:
!
      IMPLICIT NONE
      PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!  
      PUBLIC  :: YUJIMN
      PUBLIC  :: READJIMN
      PUBLIC  :: IONRATE0
!
! !PRIVATE MEMBER FUNCTIONS:
!
      PRIVATE :: READIONRATE
      PRIVATE :: IONSOIL
      PRIVATE :: GEO2MAGLAT
!
! !REVISION HISTORY: 
!  28 Sep 2008 - F. Yu       - Initial version
!  08 Nov 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! ! DEFINED PARAMETERS:
!
      ! MC     : Number of points in H2SO4 concentration dimension
      ! MT     : Number of points in temperature dimension
      ! MRH    : Number of points in relative humidity dimension
      ! MQ     : Number of points in ionization rate dimension
      ! MS     : Number of points in surface area dimension
      INTEGER, PARAMETER   :: MC  = 31
      INTEGER, PARAMETER   :: MT  = 57
      INTEGER, PARAMETER   :: MRH = 51
      INTEGER, PARAMETER   :: MQ  = 18
      INTEGER, PARAMETER   :: MS  = 12
!
! !LOCAL VARIABLES:
!
      ! C      : Values at points in H2SO4 concentration dimension
      ! T      : Values at points in temperature dimension
      ! RH     : Values at points in relative humidity dimension
      ! Q      : Values at points in ionization rate dimension
      ! S      : Values at points in surface area dimension
      ! XJIMN  : ion-mediated nucleation rates (cm-3s-1) 
      !           at all points in 5-d space
      ! XRSTAR : Critical radius (nm) at all points in 5-dimension space
      REAL*8               :: C(MC)
      REAL*8               :: RH(MRH)
      REAL*8               :: T(MT)
      REAL*8               :: Q(MQ)
      REAL*8               :: S(MS)
      REAL*8               :: XJIMN(MC,MRH,MT,MQ,MS)    
      REAL*8               :: XRSTAR(MC,MRH,MT)
 
      ! Data directory
      CHARACTER(LEN=255)   :: DATA_DIR_1x1

      CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: yujimn
!
! !DESCRIPTION: This subroutine is to calculate rates and critical cluster 
!  properties of ion-mediated nucleation (IMN) and kinetic binary homogeneous
!  nucleation (KBHN) from lookup tables using multiple-variable interpolation 
!  scheme. 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE YUJIMN( X0, Y0, Z0, U0, V0, XJ0, XR0 )
!
! !INPUT PARAMETERS: 
!
      REAL*8  :: X0    ! [H2SO4] in #/cm3  (5E5-5E8)
      REAL*8  :: Y0    ! RH in % (0.5-99.5)
      REAL*8  :: Z0    ! T (in K) (190-302)
      REAL*8  :: U0    ! Q = ionization rate (ion-pairs/cm3s) (0, 1.5-60)
      REAL*8  :: V0    ! S = surface area (um2/cm3) (1-1000)
!
! !OUTPUT PARAMETERS:
!
      REAL*8  :: XJ0   ! Nucleation rate (#/cm3s)
      REAL*8  :: XR0   ! Radius of critical cluster (nm)
!
! !REMARKS:
! Here the IMN lookup table reported in Yu (JGR,2010) and BHN lookup table 
! reproted in Yu (JGR, 2008) have been integrated. YUJIMN gives KBHN rates
! when the ionization rate is set to 0. The integrated lookup table has a size
! of ~ 170 MB and is designed for 3-D application. For quick application and 
! comparison, one can obtain nucleation rates under specified conditions using 
! an online nucleation rate calculator @ 
! http://www.albany.edu/~yfq/YuOnLineNucleation.html.
!                                                                             .
! The present lookup table should cover almost all the possible conditions in 
! the troposphere relevant to atmospheric nucleation. The range and resolution 
! in each parameter space can be extended in the future if needed.
!                                                                             .
! Written by 
! Fangqun Yu
! Atmospheric Sciences Research Center
! State University of New York at Albany
! E-mail: yfq@asrc.cestm.albany.edu; fangqun.yu@asrc.albany.edu
!                                                                             .
! Original code writted in 2006. Significnat update in 2008 and 2010. Contact
! Yu for future update or if you have questions.
!                                                                             .
! IMN lookup table reference: 
! 1. Yu, F., Ion-mediated nucleation in the atmosphere: Key controlling 
!      parameters, implications, and look-up table, J. Geophy. Res., 115, 
!      D03206, doi:10.1029/2009JD012630, 2010.
!                                                                             .
! IMN model references: 
! 2. Yu, F., From molecular clusters to nanoparticles: Second-generation 
!      ion-mediated nucleation model, Atmos. Chem. Phys., 6, 5193-5211, 2006.
! 3. Yu, F., and R. P. Turco, Ultrafine aerosol formation via ion-mediated 
!      nucleation, Geophys. Res. Lett., 27, 883-886, 2000.
!                                                                             .
! KBHN lookup table reference:
! 4. Yu, F., Updated H2SO4-H2O binary homogeneous nucleation rate look-up 
!      tables, J. Geophy. Res.,113, D24201, doi:10.1029/2008JD010527, 2008.
!                                                                             .
! KBHN model references:
! 5. Yu, F., Improved quasi-unary nucleation model for binary H2SO4-H2O 
!      homogeneous nucleation, J. Chem. Phys., 127, 054301, 2007.
! 6. Yu, F., Quasi-unary homogeneous nucleation of H2SO4-H2O, J. Chem. 
!      Phys., 122, 074501, 2005.
!
! !REVISION HISTORY: 
!  28 Sep 2008 - F. Yu       - Initial version
!  08 Nov 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL*8  :: X,Y,Z,U,V
      REAL*8  :: VOL,VOL3,FRACT,FRACT3
      REAL*8  :: dx1,dx2,dy1,dy2,dz1,dz2,du1,du2,dv1,dv2
      REAL*8  :: dx,dy,dz,du,dv
      INTEGER :: IC1, IC2, JRH1, JRH2, KT1, KT2, IQ1, IQ2, IS1, IS2
      INTEGER :: IC, JRH, KT, IQ, IS
!
! to avoid the input values to be changed due to out of the range reset
!
      X = X0
      Y = Y0
      Z = Z0
      U = U0
      V = V0
!
! The present lookup table should cover almost all the possible conditions in 
! the troposphere relevant to atmospheric nucleation. The range and resolution 
! in each parameter space can be extended in the future if needed.
! If the inputed values are out of the lookup table valid ranges, set them to 
! boundary values for now. Care should be taken if your inputted values are 
! frequently out of the specified ranges.
! 
      IF(U.LE.1.E-20) U=1.E-20  ! i.e., binary homogeneous nucleation 
!
      IF(X.LT.C(1)) THEN
         WRITE(86,10) X, C(1), C(1)
         X = C(1)
      ELSEIF(X.GT.C(MC)) THEN
         WRITE(86,11) X, C(MC), C(MC)
         X =C(MC)
      ENDIF

      IF(Y.LT.RH(1)) THEN
!        WRITE(86,12) Y, RH(1), RH(1)
         Y =RH(1) 
      ELSEIF(Y.GT.RH(MRH)) THEN
         WRITE(86,13) Y, RH(MRH), RH(MRH)
         Y =RH(MRH)
      ENDIF
      
      IF(Z.LT.T(1)) THEN
!        WRITE(86,14) Z, T(1), T(1)
         Z =T(1)
      ELSEIF(Z.GT.T(MT)) THEN
!        WRITE(86,15) Z, T(MT), T(MT)
         Z =T(MT)
      ENDIF

      IF(U.LT.Q(1)) THEN
         WRITE(86,16) U, Q(1), Q(1)
         U =Q(1)
      ELSEIF(U.GT.Q(MQ)) THEN
         WRITE(86,17) U, Q(MQ), Q(MQ)
         U =Q(MQ)
      ENDIF

      IF(V.LT.S(1)) THEN
!        WRITE(86,18) V, S(1), S(1)
         V =S(1)
      ELSEIF(V.GT.S(MS)) THEN
!        WRITE(86,19) V, S(MS), S(MS)
         V =S(MS)
      ENDIF

 10   FORMAT("IMN WARNING: INPUTED [H2SO4]=",ES9.2,"<",ES9.2,
     &     " set it to ",ES9.2)
 11   FORMAT("IMN WARNING: INPUTED [H2SO4]=",ES9.2,">",ES9.2,
     &     " set it to ",ES9.2)
 12   FORMAT("IMN WARNING: INPUTED RH =",F5.1,"% <",F5.1,
     &     "% set it to ",F5.1,"%")
 13   FORMAT("IMN WARNING: INPUTED RH =",F5.1,"% >",F5.1,
     &     "% set it to ",F5.1,"%")
 14   FORMAT("IMN WARNING: INPUTED T =",F6.1,"K <",F6.1,
     &     "K set it to ",F6.1,"K")
 15   FORMAT("IMN WARNING: INPUTED T =",F6.1,"K >",F6.1,
     &     "K set it to ",F6.1,"K")
 16   FORMAT("IMN WARNING: INPUTED Q =",F6.1," <",F6.1,
     &     " ion-pair/cm3s set it to ",F6.1)
 17   FORMAT("IMN WARNING: INPUTED Q =",F6.1," >",F6.1,
     &     " ion-pair/cm3s set it to ",F6.1)
 18   FORMAT("IMN WARNING: INPUTED S =",F6.1," <",F6.1,
     &     " um2/cm3 set it to ",F6.1)
 19   FORMAT("IMN WARNING: INPUTED S =",F6.1," >",F6.1,
     &     " um2/cm3 set it to ",F6.1)

      IC1 =MAX0(INT(1.+10.*LOG10(X/5.E5)),1)
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
!     
        
      IF(U.LT.Q(2)) THEN
         IQ1 =1.
      ELSE
         IQ1 = MAX0(INT(2.+10.*LOG10(U/Q(2))),2)
      ENDIF
      IQ2 = MIN0(IQ1 + 1,MQ)
      IF(IQ2.EQ.MQ) IQ1=MQ-1
!
      IF(V.LT.10.0) THEN
         IS1 =1.
      ELSE
         IS1 = MAX0(INT(2.+5.*LOG10(V/10.)),2)
      ENDIF
      IS2 = MIN0(IS1 + 1,MS)
      IF(IS2.EQ.MS) IS1=MS-1
!     
      dx1 = LOG10(X/C(IC1))     ! logJ log[H2SO4] interpolation
      dx2 = LOG10(C(IC2)/X)
      dy1 = LOG10(Y/RH(JRH1))
      dy2 = LOG10(RH(JRH2)/Y)
      dz1 = Z-T(KT1)
      dz2 = T(KT2)-Z
      
      du1 = U - Q(IQ1)
      du2 = Q(IQ2) - U
      dv1 = V- S(IS1)
      dv2 = S(IS2) - V
!
      XJ0 = 0.
      XR0 = 0.
!
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
               XR0 = XR0 + FRACT3*XRSTAR(IC,JRH,KT)
               
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
!                    WRITE(6,30)IC,JRH,KT,IQ,IS,10.**XJIMN(IC,JRH,KT,IQ,IS),
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
 30   FORMAT(I3, I3, I3, I3, I3, 10(1PE10.3))
 20   FORMAT(10(1PE10.3))

      END SUBROUTINE YUJIMN
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: readjimn
!
! !DESCRIPTION: Read in the integrated IMN and KBHN lookup tables.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE READJIMN( DATA_DIR_1x1a )
!
! !INPUT PARAMETERS: 
!
      CHARACTER(LEN=255)  :: DATA_DIR_1x1a    ! Data directory
!
! !REMARKS:
!  WRITTEN by Fangqun Yu, SUNY-Albany, 2006 (Revised, 8/2008)
!                                                                             .
!  IMN lookup table references: 
!  1. Yu, F., Ion-mediated nucleation in the atmosphere: Key controlling 
!         parameters, implications, and look-up table, J. Geophy. Res., 115, 
!         D03206, doi:10.1029/2009JD012630, 2010.
!                                                                             .
!  KBHN lookup table reference:
!  2. Yu, F.,Updated H2SO4-H2O binary homogeneous nucleation rate look-up 
!         tables, J. Geophy. Res.,113, D24201, doi:10.1029/2008JD010527, 2008.
! 
! !REVISION HISTORY: 
!  28 Sep 2008 - F. Yu       - Initial version
!  08 Nov 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER :: IC, IRH, IT, IQ, IS 
      REAL*8  :: C11,Q11,S11

!     CHARACTER*2 YPATH
!     YPATH = './'
      CHARACTER*80 YPATH

      DATA_DIR_1x1= DATA_DIR_1x1a
      YPATH = TRIM(DATA_DIR_1x1)//'/APM_201011/IMN_LT/'
      WRITE(6,*)"Read IMN look-up tables"
      
      open(31,file=TRIM(YPATH)//'Yu_IMN_J5D.txt',form='formatted')
      open(33,file=TRIM(YPATH)//'Yu_IMN_Rstar3D.txt',form='formatted')
      
      open(41,file=TRIM(YPATH)//'Yu_IMN_1H2SO4.txt',form='formatted')
      open(42,file=TRIM(YPATH)//'Yu_IMN_2RH.txt',form='formatted')
      open(43,file=TRIM(YPATH)//'Yu_IMN_3T.txt',form='formatted')
      open(44,file=TRIM(YPATH)//'Yu_IMN_4Q.txt',form='formatted')
      open(45,file=TRIM(YPATH)//'Yu_IMN_5S.txt',form='formatted')
!     
      READ(41,100)(C(IC),IC=1,MC)
      WRITE(6,*)"[H2SO4](IC), IC=1, ", MC, ":"
      WRITE(6,100)(C(IC),IC=1,MC)
!     
      READ(42,100)(RH(IRH),IRH=1,MRH)
      WRITE(6,*)"RH(IRH), IRH=1, ", MRH, ":"
      WRITE(6,100)(RH(IRH),IRH=1,MRH)
!     
      READ(43,100)(T(IT),IT=1,MT)
      WRITE(6,*)"T(IT), IT=1, ", MT, ":"
      WRITE(6,100)(T(IT),IT=1,MT)
!     
      READ(44,100)(Q(IQ),IQ=1,MQ)
      WRITE(6,*)"Q(I), I=1, ", MQ, ":"
      WRITE(6,100)(Q(IQ),IQ=1,MQ)
!     
      READ(45,100)(S(IS),IS=1,MS)
      WRITE(6,*)"S(IS), IS=1, ", MS, ":"
      WRITE(6,100)(S(IS),IS=1,MS)
!
! Use the formula to calculate C and Q to get values with more digits, otherwise
! may cause problem when input C and Q are very clsoe to C(IC),Q(IQ) as
! IC and IQ are decided with formula 
!
      C(1) = 5.0E5
      DO IC = 2, MC
         C11 = C(IC)                                                          
         C(IC) = C(IC-1)*10.**(0.1)

         IF(abs(1.-C11/C(IC)).GT.0.02) THEN                                  
            write(6,*)"need check JIMN look-up table inputs"                  
            stop                                                              
         ENDIF                                                                
      ENDDO

      DO IQ = 1, MQ
         Q11 = Q(IQ)                                                          
         IF(IQ.EQ.1) THEN
            Q(1) =1.E-30
         ELSE
            Q(IQ) = 1.5*10.**(0.1*float(IQ-2))
         ENDIF
         IF(abs(1.-Q11/Q(IQ)).GT.0.02) THEN
            write(6,*)"need check JIMN look-up table inputs"
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
            write(6,*)"need check JIMN look-up table inputs"
            stop
         ENDIF
      ENDDO

!
! READ in formatted 5-D Table
!
      DO IS =1, MS
      DO IT = 1,MT
      DO IRH = 1,MRH
      DO IQ =1, MQ
         READ(31,201)(XJIMN(IC,IRH,IT,IQ,IS),IC = 1,MC)

         DO IC=1, MC

! Due to high sensitivity of J to key parameters, use logJ to interpolate

            XJIMN(IC,IRH,IT,IQ,IS)=LOG10(XJIMN(IC,IRH,IT,IQ,IS))
         ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO

! Critical cluster properties depend on T, RH, [H2SO4] only
      DO IT = 1,MT
         DO IRH = 1, MRH
            READ(33,203)(XRSTAR(IC,IRH,IT),
     &           IC=1,MC)
         ENDDO                  ! RH
      ENDDO                     !T

      CLOSE(31)
      CLOSE(33)
      CLOSE(41)
      CLOSE(42)
      CLOSE(43)
      CLOSE(44)
      CLOSE(45)
!     
 100  FORMAT(100(1PE10.3))
 200  FORMAT(100(1PE9.2))
 201  FORMAT(100(1PE9.2))
 202  FORMAT(100F5.1)
 203  FORMAT(100F5.2)
 204  FORMAT(100F6.3)

      END SUBROUTINE READJIMN
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ionrate0 
!
! !DESCRIPTION: Subroutine IONRATE calculate ionization rate 
!  (ZQ: ion-pairs/cm3s) for given surface type, longitude (in degree), 
!  latitude (in degree), and pressure (mb).
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE IONRATE0( ISURF, YPSURF, XLON, XLAT, XP, ZQ )
!
! !INPUT PARAMETERS: 
!
      INTEGER   :: ISURF    ! Surface type (1=land, 0=ocean, ice, and snow)
      REAL*8    :: YPSURF   ! Surface pressure [hPa]
      REAL*8    :: XLON     ! Longitude [degrees]
      REAL*8    :: XLAT     ! Latitude [degrees]
      REAL*8    :: XP       ! Grid box pressure [hPa], from 5-1015 hPa
!
! !OUTPUT PARAMETERS:
!
      REAL*8    :: ZQ       ! Ionization rate [ion pairs/cm3/s]
!
! !REMARKS:
!  Written by Fangqun Yu and Gan Luo, SUNY-Albany, 2010 
!  (yfq@asrc.cestm.albany.edu)
!                                                                             .
!  Ionization lookup table
!  YQ(1,1):    Q for maglat = 0,  p = 5 mb
!  YQ(91,203): Q for maglat = 90, p = 1015 mb
!  YQ(L,K):    Q for maglat = L-1,p = K*5 mb
! 
! !REVISION HISTORY: 
!  28 Sep 2008 - F. Yu       - Initial version
!  08 Nov 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER       :: L, K
      REAL*8        :: MAGLAT, XMAGLAT, YPR, YQSOIL
      REAL*8,  SAVE :: YQ(91,203)  
      LOGICAL, SAVE :: FIRST = .TRUE.
!
! Read in the ionization rate lookup table
!
      IF(FIRST) THEN
         CALL READIONRATE(YQ)
         FIRST = .FALSE.
      ENDIF
!
! Find the magnetic latitude based on (LON, LAT)

      CALL GEO2MAGLAT(XLAT,XLON,XMAGLAT)
      MAGLAT= abs(XMAGLAT) !  magnetic latitude in degree

      L = INT(MAGLAT+0.5) + 1
      K = MIN( MAX ( INT( XP/5. + 0.5 ), 1 ), 203 )
      ZQ = YQ(L,K)   ! GCR ionization rate from the lookup table

      IF(ISURF.EQ.1) THEN  ! Contribution from radioactive material from soil
        CALL IONSOIL(XP,YPSURF,YQSOIL)
      ELSE
        YQSOIL = 0.
      ENDIF
!      WRITE(6,101)XLON,XLAT,XMAGLAT,XP, L, K,YQSOIL,ZQ

      ZQ = ZQ + YQSOIL
 101  FORMAT(4F7.1, I4, I4, 2F7.1)
      RETURN
      END SUBROUTINE IONRATE0
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: readionrate
!
! !DESCRIPTION: Read pre-calculated GCR ionization rate lookup table
!  The lookup table is generated based on the scheme given in
!  Usoskin and Kovaltsov (2006).
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE READIONRATE( YQ )
!
! !OUTPUT PARAMETERS:
!
      ! YQ: ionization rate in ion-pairs/cm3s
      ! maglat is magnitude latitude in degree
      REAL*8 :: YQ(91,203) 
!
! !REMARKS:
!  YQ(1,1):    Q for maglat = 0,  p = 5 mb
!  YQ(91,203): Q for maglat = 90, p = 1015 mb
!  YQ(L,K):    Q for maglat = L-1,p = K*5 mb
! 
! !REVISION HISTORY: 
!  28 Sep 2008 - F. Yu       - Initial version
!  08 Nov 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER   :: L, K, KP
      CHARACTER*80 YPATH

      WRITE(6,*)"Read in the ionization rate lookup table"
      YPATH = TRIM(DATA_DIR_1x1)//'/APM_201011/'
      CLOSE(30)
      OPEN(30,file=TRIM(YPATH)//'YIONRATE.txt',status='old')
      READ(30,*)   ! first line is magnetic latitude in degree
      DO K=1,203      ! KP is pressure in mb, YQ in ion-pairs/cm3s
         READ(30,100)KP,(YQ(L,K),L=1,91)
      ENDDO
 100  FORMAT(I5,91F6.2)
      CLOSE(30)
      RETURN
      END SUBROUTINE READIONRATE
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ionsoil
!
! !DESCRIPTION: Calculate ionization rate (ion-pairs/cm3s) due to Gama rays, 
!  Radon (over land)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE IONSOIL( YPR, YPSURF, YQSOIL )
!
! !INPUT PARAMETERS: 
!
      REAL*8    :: YPR
      REAL*8    :: YPSURF
!
! !OUTPUT PARAMETERS:
!
      REAL*8    :: YQSOIL
!
! !REVISION HISTORY: 
!  28 Sep 2008 - F. Yu       - Initial version
!  08 Nov 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
      INTEGER, PARAMETER   :: MAXH = 21
!
! !LOCAL VARIABLES:
!
      INTEGER   :: IH, K
      REAL*8    :: XH0, XH, YQ, YQGAMA, YQRADON
      REAL*8    :: YH(MAXH), QGAMA(MAXH), QRADON(MAXH)
      REAL*8    :: MAGLAT
      REAL*8    :: XHSURF

      DATA (YH(k),k=1,21)/0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0, ! in km
     &                    2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,15.0/
      DATA (QGAMA(k),k=1,21)/4.5,1.25,0.21,0.0,0.0,0.0,0.0,0.0,0.0,
     &              0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0/
      DATA (QRADON(k),k=1,21)/3.5,3.24,3.0,2.65,2.43,2.19,1.84,1.36,
     &                        0.97, 0.74,0.56,0.13,
     &                        0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0/

! Get altitude (km) from pressure (YPR in mb) based on standard atmosphere
!
      XH0 = 44.3308 - 4.94654*(100.*YPR)**0.190264   ! convert press to alt (km)
! Surface height
      XHSURF = 44.3308 - 4.94654*(100.*YPSURF)**0.190264  

      XH = XH0 - XHSURF


      IF(XH.LT.5) THEN   ! over land, no gama and radon above 5 km
         IF(XH.LT.1) THEN
            IH = INT(XH*10.)+1
         ELSE
            IH = 10.+INT(XH)
         ENDIF

         IH=MIN(IH,20)
         IH=MAX(IH,1)

         YQGAMA=QGAMA(IH)+(XH-YH(IH))*(QGAMA(IH+1)-QGAMA(IH))
         YQRADON=QRADON(IH)+(XH-YH(IH))*(QRADON(IH+1)-QRADON(IH))
         YQSOIL = YQGAMA + YQRADON
!         WRITE(6,*)XH0,XHSURF, IH,YQGAMA,YQRADON
      ELSE 
         YQSOIL = 0.
      ENDIF

      END SUBROUTINE IONSOIL
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: geo2maglat
!
! !DESCRIPTION: Subroutine GEO2MAGLAT finds magnetic latitude from geo 
!  latitude and longitude.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE GEO2MAGLAT( LAT0, LON0, MAGLAT )
!
! !INPUT PARAMETERS: 
!
      REAL*8 :: lat0     ! Geographic latitude [degrees]
      REAL*8 :: lon0     ! Geographic longitude [degrees]
!
! !OUTPUT PARAMETERS:
!
      REAL*8 :: maglat   ! Magnetic latitude [degrees]
! 
! !REVISION HISTORY: 
!  28 Sep 2008 - F. Yu       - Initial version
!  08 Nov 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL*8 :: yz,yp,yx,yy,PI,YD2R,YGLA,YGLO
      REAL*8 :: az,ap,ax,ay,dot,theta,lon,lat

      PI = 3.1415926
      YD2R = PI/180.
       
      YGLA = 80.*YD2R              ! magnetic dipole north pole latitude
      YGLO = -110.*YD2R            ! magnetic dipole north pole longitude
      
      lat = lat0 * YD2R
      lon = lon0 * YD2R
      
      yz = sin(YGLA)
      yp = cos(YGLA)
      yx = yp * cos(YGLO)
      yy = yp * sin(YGLO)
      
      az = sin(lat)
      ap = cos(lat)
      ax = ap * cos(lon)
      ay = ap * sin(lon)
      dot = ax*yx + ay*yy + az*yz
      theta = acos(dot)            ! theta is the magnetic 
                                   ! colatitude of the point a
      
      maglat = (PI/2.-theta)/YD2R
       
      END SUBROUTINE GEO2MAGLAT
!EOC
      END MODULE APM_NUCL_MOD

