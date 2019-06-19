#if defined ( APM )
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

!APM2+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      PUBLIC :: IONRATE
      PUBLIC :: BHN_VEHKAMAKI
      PUBLIC :: BHN_KULMALA98
      PUBLIC :: LovejoyIN
      PUBLIC :: terpara77
      PUBLIC :: PORGANIC
!APM2+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

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
      REAL*8               :: XJIMN8(MC)
      REAL               :: XJIMN(MC,MRH,MT,MQ,MS)    
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
                     XJ0 = XJ0 + FRACT*DBLE(XJIMN(IC,JRH,KT,IQ,IS))
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
      CHARACTER*999 YPATH

      DATA_DIR_1x1= DATA_DIR_1x1a
      YPATH = TRIM(DATA_DIR_1x1)//'/APM_data_201906/IMN_LT/'
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
         READ(31,201)(XJIMN8(IC),IC = 1,MC)


         DO IC=1, MC

! Due to high sensitivity of J to key parameters, use logJ to interpolate

            XJIMN(IC,IRH,IT,IQ,IS)=LOG10(REAL(XJIMN8(IC)))

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
      SUBROUTINE IONRATE0(DATA_DIR_1x1a,IY,ISURF,YPSURF,XLON,XLAT,XP,ZQ)
!
! !INPUT PARAMETERS: 
!
      INTEGER   :: IY
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

      CHARACTER(LEN=255)  :: DATA_DIR_1x1a    ! Data directory
      INTEGER       :: L, K
      REAL*8        :: MAGLAT, XMAGLAT, YPR, YQSOIL
      REAL*8,  SAVE :: YQ(91,203)  
      LOGICAL, SAVE :: FIRST = .TRUE.
!
! Read in the ionization rate lookup table
!
      IF(FIRST) THEN
         CALL READIONRATE(DATA_DIR_1x1a,IY,YQ)
         FIRST = .FALSE.
      ENDIF
!
! Find the magnetic latitude based on (LON, LAT)

      CALL GEO2MAGLAT(XLAT,XLON,XMAGLAT)
      MAGLAT= abs(XMAGLAT) !  magnetic latitude in degree

      L = INT(MAGLAT+0.5) + 1
      K = MIN(MAX(INT(XP/5.+0.5),1),203)
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
      SUBROUTINE READIONRATE(DATA_DIR_1x1a,IY, YQ )

      CHARACTER(LEN=255)  :: DATA_DIR_1x1a    ! Data directory
      INTEGER   :: IY
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
      CHARACTER*999 YPATH

      WRITE(6,*)"Read in the ionization rate lookup table"
      DATA_DIR_1x1= DATA_DIR_1x1a
      YPATH = TRIM(DATA_DIR_1x1)//'/APM_data_201906/APMTABLES/'
      CLOSE(30)
      WRITE(6,*)"YPATH: ",TRIM(YPATH)
      IF(IY.EQ.1) THEN
       OPEN(30,file=TRIM(YPATH)//'YIONRATE1996.txt',status='old')
       WRITE(6,*) "READ YIONRATE1996.txt"
      ELSEIF(IY.EQ.-1)THEN
       OPEN(30,file=TRIM(YPATH)//'YIONRATE1989.txt',status='old')
       WRITE(6,*) "READ YIONRATE1989.txt"
      ELSE
       OPEN(30,file=TRIM(YPATH)//'YIONRATE.txt',status='old')
       WRITE(6,*) "READ YIONRATE.txt"
      ENDIF
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


! APM2+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine IONRATE(DATA_DIR_1x1a,
     &                   ISURF,XLON,XLAT,YPR,IYEAR,IMONTH,YQ)

!******************************************************************************
! Subroutine IONRATE calculate ionization rate (YQ: ion-pairs/cm3s) for given 
! surface type,longitude (in degree), latitude (in degree), pressure (mb), year,
! month  (fyu, 2006; updated 2008)

      CHARACTER(LEN=255)  :: DATA_DIR_1x1a    ! Data directory
      INTEGER, PARAMETER   :: MAXH = 21
      INTEGER   :: ISURF, IYEAR, IMONTH, IH, K
      REAL*8    :: XLAT, XLON, YPR, XH, YQ, YQGAMA, YQRADON, YQGCR
      REAL*8    :: YH(MAXH), QGAMA(MAXH), QRADON(MAXH)

! Ionuzation rate (ion-pairs/cm3s) due to Gama rays, Radon (over land)
!
      DATA (YH(k),k=1,21)/0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0, ! in km
     &                    2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,15.0/
      DATA (QGAMA(k),k=1,21)/4.5,1.25,0.21,0.0,0.0,0.0,0.0,0.0,0.0,
     &              0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0/
      DATA (QRADON(k),k=1,21)/3.5,3.24,3.0,2.65,2.43,2.19,1.84,1.36,
     &                        0.97, 0.74,0.56,0.13,
     &                        0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0/

! Get altitude (km) from pressure (YPR in mb) based on standard atmosphere
!
      XH = 44.3308 - 4.94654*(100.*YPR)**0.190264   ! convert press to alt (km)

      CALL IONGCR(DATA_DIR_1x1a,XLON,XLAT,YPR,IYEAR,IMONTH,YQGCR)

      IF(ISURF.EQ.1.and.XH.LT.5) THEN   ! over land, no gama and radon above 5 km
         IF(XH.LT.1) THEN
            IH = INT(XH*10.)+1
         ELSE
            IH = 10.+INT(XH)
         ENDIF

         !bug gan luo fixed
         IH=MIN(IH,20)
         IH=MAX(IH,1)

         YQGAMA=QGAMA(IH)+(XH-YH(IH))*(QGAMA(IH+1)-QGAMA(IH))
         YQRADON=QRADON(IH)+(XH-YH(IH))*(QRADON(IH+1)-QRADON(IH))
         YQ = YQGCR + YQGAMA + YQRADON
      ELSE  ! over ice or water
         YQ = YQGCR
      ENDIF

      end subroutine IONRATE
!------------------------------------------------------------------------------

       subroutine  IONGCR(DATA_DIR_1x1a,lon,lat,airdin,year,month,qout)

!******************************************************************************
! Subroutine IONGCR is to obtain ionization induced by cosmic ray based on 
! Usoskin and Kovaltsov, JGR,111,D21206, 2006. Coded by Gan Luo, ASRC, Albany, 
! NY, 2006-11-20 (Thanks to professor Fangqun Yu and Dr Usoskin)
!
! Input: longitude (degree), latitude (degree), 
!       atmospheric depth (or pressure, in mb), 
!       month and year
! Output: cosmic ray induced ionization (ion-pairs/cm3s)
!
! Modified by fyu, 11/26/2006, 2008

       CHARACTER(LEN=255)  :: DATA_DIR_1x1a    ! Data directory
       LOGICAL, SAVE    :: FIRST = .TRUE.
       integer :: i,j,k,t,im,jy
       integer :: year, month
       
       real*8 :: lon,lat,maglat,cosmiclat,airdin,airout,ev
       real*8 :: pc,tforp,tfora,pt1,pt2,tnew1,tnew2,yxt_pad(9), 
     &            yxt_aad(9)
       real*8 :: jlisp1,jlisp2,jlisa1,jlisa2,partp1,partp2,parta1,
     &           parta2
       real*8 :: qforp,qfora,qout

       REAL*8, SAVE :: cosmic(14,55)

       real*8 :: aird(21),airro(21),tvalue(9),
     &           yxt_p(21,9),yxt_a(21,9)

!***********************************************************************
       data (aird(k),k=1,21)/25.,75.,125.,175.,225.,275.,325.,375.,425.,
     &    475.,525.,575.,625.,675.,725.,775.,825.,875.,925.,975.,1025./

       data (airro(k),k=1,21)/3.80E-05,1.20E-04,2.00E-04,2.70E-04,
     &  3.50E-04,4.20E-04,4.80E-04,
     &  5.40E-04,5.90E-04,6.50E-04,7.10E-04,7.60E-04,8.20E-04,8.70E-04,
     &  9.20E-04,9.70E-04,1.00E-03,1.10E-03,1.10E-03,1.20E-03,1.20E-03/

       data (tvalue(t),t=1,9)/0.1,0.3,1.,3.,10.,30.,100.,300.,1000./

       data ((yxt_p(k,t),t=1,9),k=1,21)/
     & 3.40E+02,4.10E+05,4.60E+05,6.00E+05,1.30E+06,2.30E+06,4.90E+06,
     & 8.50E+06,1.50E+07,9.80E+01,4.30E+04,3.30E+05,6.30E+05,1.80E+06,
     & 4.20E+06,1.00E+07,2.30E+07,5.70E+07,5.00E+01,4.40E+03,2.10E+05,
     & 5.20E+05,1.60E+06,4.40E+06,1.20E+07,2.90E+07,8.30E+07,2.00E+01,
     & 2.50E+03,1.30E+05,4.00E+05,1.30E+06,3.80E+06,1.10E+07,3.00E+07,
     & 9.30E+07,7.90E+00,1.40E+03,7.40E+04,2.90E+05,9.90E+05,3.10E+06,
     & 9.90E+06,2.80E+07,9.20E+07,4.90E+00,8.50E+02,4.10E+04,2.10E+05,
     & 7.30E+05,2.40E+06,8.10E+06,2.40E+07,8.20E+07,2.00E+00,4.90E+02,
     & 2.10E+04,1.50E+05,5.30E+05,1.80E+06,6.30E+06,2.00E+07,6.90E+07,
     & 8.10E-01,2.90E+02,1.10E+04,1.00E+05,3.70E+05,1.40E+06,5.00E+06,
     & 1.60E+07,5.70E+07,3.10E-01,1.80E+02,6.20E+03,7.20E+04,2.60E+05,
     & 1.00E+06,3.90E+06,1.30E+07,4.50E+07,0.00E+00,1.10E+02,3.90E+03,
     & 5.00E+04,1.90E+05,7.40E+05,3.10E+06,1.10E+07,3.70E+07,0.00E+00,
     & 6.30E+01,2.60E+03,3.50E+04,1.40E+05,5.60E+05,2.40E+06,8.50E+06,
     & 2.90E+07,0.00E+00,3.80E+01,1.50E+03,2.40E+04,9.80E+04,4.40E+05,
     & 1.90E+06,6.90E+06,2.40E+07,0.00E+00,1.80E+01,1.10E+03,1.60E+04,
     & 7.20E+04,3.50E+05,1.60E+06,5.70E+06,2.00E+07,0.00E+00,1.30E+01,
     & 7.30E+02,1.00E+04,5.30E+04,2.70E+05,1.30E+06,4.70E+06,1.70E+07,
     & 0.00E+00,6.30E+00,4.20E+02,6.90E+03,3.80E+04,2.20E+05,1.10E+06,
     & 3.90E+06,1.40E+07,0.00E+00,5.00E+00,3.10E+02,5.10E+03,2.90E+04,
     & 1.80E+05,9.30E+05,3.40E+06,1.20E+07,0.00E+00,4.90E+00,1.80E+02,
     & 3.00E+03,2.10E+04,1.50E+05,8.30E+05,3.00E+06,1.00E+07,0.00E+00,
     & 1.20E+00,1.10E+02,1.70E+03,1.40E+04,1.30E+05,7.30E+05,2.60E+06,
     & 9.10E+06,0.00E+00,2.30E+00,9.90E+01,1.70E+03,1.20E+04,1.10E+05,
     & 6.60E+05,2.40E+06,8.10E+06,0.00E+00,8.90E-01,6.00E+01,8.00E+02,
     & 9.00E+03,9.60E+04,6.00E+05,2.20E+06,7.30E+06,0.00E+00,1.80E-01,
     & 6.80E+01,7.00E+02,7.40E+03,8.50E+04,5.50E+05,2.00E+06,6.70E+06/
       
       data ((yxt_a(k,t),t=1,9),k=1,21)/
     & 1.00E+03,3.70E+05,4.10E+05,5.80E+05,1.30E+06,2.50E+06,5.50E+06,
     & 8.60E+06,2.20E+07,3.40E+02,3.50E+04,3.00E+05,6.20E+05,1.90E+06,
     & 4.60E+06,1.20E+07,2.70E+07,6.80E+07,1.60E+02,9.90E+03,1.90E+05,
     & 5.30E+05,1.80E+06,4.80E+06,1.40E+07,3.50E+07,9.60E+07,6.50E+01,
     & 5.60E+03,1.30E+05,4.10E+05,1.40E+06,4.20E+06,1.30E+07,3.70E+07,
     & 1.10E+08,2.80E+01,3.40E+03,7.80E+04,3.10E+05,1.10E+06,3.40E+06,
     & 1.10E+07,3.50E+07,1.00E+08,2.40E+01,2.00E+03,4.70E+04,2.20E+05,
     & 7.50E+05,2.60E+06,8.90E+06,3.10E+07,9.70E+07,1.40E+01,1.20E+03,
     & 2.80E+04,1.60E+05,5.50E+05,2.00E+06,6.90E+06,2.60E+07,8.30E+07,
     & 2.90E+00,7.50E+02,1.70E+04,1.10E+05,3.90E+05,1.50E+06,5.30E+06,
     & 2.00E+07,6.80E+07,1.90E+00,4.60E+02,1.20E+04,8.00E+04,2.80E+05,
     & 1.10E+06,4.10E+06,1.60E+07,5.40E+07,1.10E+00,2.80E+02,7.80E+03,
     & 5.30E+04,2.00E+05,8.50E+05,3.10E+06,1.20E+07,4.20E+07,5.70E-01,
     & 1.80E+02,5.00E+03,3.70E+04,1.40E+05,6.40E+05,2.50E+06,9.90E+06,
     & 3.30E+07,2.30E-01,1.00E+02,2.80E+03,2.40E+04,1.00E+05,4.80E+05,
     & 2.00E+06,8.00E+06,2.60E+07,9.20E-03,7.70E+01,2.10E+03,1.80E+04,
     & 7.90E+04,3.80E+05,1.60E+06,6.50E+06,2.10E+07,0.00E+00,3.70E+01,
     & 1.50E+03,1.20E+04,5.00E+04,2.80E+05,1.40E+06,5.40E+06,1.70E+07,
     & 0.00E+00,2.00E+01,6.60E+02,7.40E+03,3.60E+04,2.30E+05,1.10E+06,
     & 4.60E+06,1.40E+07,0.00E+00,1.60E+01,6.30E+02,5.50E+03,3.00E+04,
     & 1.80E+05,9.70E+05,4.00E+06,1.20E+07,0.00E+00,8.20E+00,4.40E+02,
     & 4.10E+03,1.90E+04,1.50E+05,8.60E+05,3.40E+06,1.00E+07,0.00E+00,
     & 9.00E+00,2.50E+02,2.40E+03,1.40E+04,1.30E+05,7.50E+05,3.00E+06,
     & 8.80E+06,0.00E+00,6.70E+00,1.80E+02,1.70E+03,1.10E+04,1.10E+05,
     & 6.60E+05,2.70E+06,7.80E+06,0.00E+00,3.60E+00,8.10E+01,1.10E+03,
     & 7.70E+03,9.80E+04,6.10E+05,2.50E+06,7.10E+06,0.00E+00,1.70E+00,
     & 4.70E+01,1.20E+03,6.40E+03,8.80E+04,5.70E+05,2.30E+06,6.50E+06/
       
!***********************************************************************

       IF(FIRST) THEN
         close(100)
         OPEN(100,
     & FILE=TRIM(DATA_DIR_1x1a)//'/APM_data_201906/Phi_mon.txt')
         do j=1,55
            read(100,*)(cosmic(i,j),i=1,14)
         enddo
         close(100)
         FIRST = .FALSE.
       ENDIF
       
!ev=cosmic(14,40)/1000.0 !modulation potential (i,j) i=month+1;j=year-1949

       im = month + 1
       jy = year -1949

       IF(JY.GT.55) JY = 55 ! current scheme only valid for 1950 - 2004
       IF(JY.LT.1) JY = 1 ! current scheme only valid for 1950 - 2004

       ev=cosmic(im,jy)/1000.0 !modulation potential (i,j) i=month+1;j=year-1949

       do k=1,20
        if(airdin>aird(k).and.airdin<=aird(k+1))then
         do t=1,9
          yxt_pad(t)=yxt_p(k,t)+(airdin-aird(k))
     &      *(yxt_p((k+1),t)-yxt_p(k,t))/(aird(k+1)-aird(k))
          yxt_aad(t)=yxt_a(k,t)+(airdin-aird(k))
     &      *(yxt_a((k+1),t)-yxt_a(k,t))/(aird(k+1)-aird(k))
         enddo
         airout=airro(k)+(airdin-aird(k))
     &      *(airro(k+1)-airro(k))/(aird(k+1)-aird(k))
        endif
       enddo
       if(airdin<=aird(1))then
        do t=1,9
         yxt_pad(t)=yxt_p(1,t)
         yxt_aad(t)=yxt_a(1,t)
        enddo
        airout=airro(1)
       endif
       if(airdin>aird(21))then
        do t=1,9
         yxt_pad(t)=yxt_p(21,t)
         yxt_aad(t)=yxt_a(21,t)
        enddo
        airout=airro(21)
       endif
       
       call geo2maglat(lat,lon,maglat)

       cosmiclat=3.14159265*maglat/180.0
       
       pc=1.9*7.8*(cos(cosmiclat)*cos(cosmiclat)
     &      *cos(cosmiclat)*cos(cosmiclat))

       tforp=sqrt(pc*pc+0.938*0.938)-0.938
       tfora=sqrt(0.25*pc*pc+0.938*0.938)-0.938
!***********************************************************************
       qforp=0.0
       qfora=0.0
       qout=0.0
       
       do t=1,8

        if(tvalue(t+1)>=tforp)then

         if(tvalue(t)<tforp)then
          tnew1=tforp+ev
         else
          tnew1=tvalue(t)+ev
         endif
         tnew2=tvalue(t+1)+ev
       
         pt1=sqrt(tnew1*(tnew1+2.0*0.938))
         pt2=sqrt(tnew2*(tnew2+2.0*0.938))
       
         jlisp1=1.9*(pt1**(-2.78))/(1.0+0.4866*(pt1**(-2.51)))
         jlisp2=1.9*(pt2**(-2.78))/(1.0+0.4866*(pt2**(-2.51)))

         if(tvalue(t)<tforp)then
          partp1=tforp*(tforp+2.0*0.938)/(tnew1*(tnew1+2.0*0.938))
         else
          partp1=tvalue(t)*(tvalue(t)+2.0*0.938)
     &                    /(tnew1*(tnew1+2.0*0.938))
         endif
         partp2=tvalue(t+1)*(tvalue(t+1)+2.0*0.938)
     &                   /(tnew2*(tnew2+2.0*0.938))
       
         if(tvalue(t)<tforp)then
          qforp=qforp+0.5*(jlisp2*partp2*yxt_pad(t+1)+jlisp1*partp1
     &        *(yxt_pad(t)+(tforp-tvalue(t))*(yxt_pad(t+1)-yxt_pad(t))
     &        /(tvalue(t+1)-tvalue(t))))*(tvalue(t+1)-tforp)
         else
          qforp=qforp+0.5*(jlisp2*partp2*yxt_pad(t+1)+jlisp1*partp1
     &           *yxt_pad(t))*(tvalue(t+1)-tvalue(t))
         endif
       
        endif

        if(tvalue(t+1)>=tfora)then
       
         if(tvalue(t)<tfora)then
          tnew1=tfora+ev*0.5
         else
          tnew1=tvalue(t)+ev*0.5
         endif
         tnew2=tvalue(t+1)+ev*0.5

         pt1=sqrt(tnew1*(tnew1+2.0*0.938))
         pt2=sqrt(tnew2*(tnew2+2.0*0.938))

         jlisa1=0.1425*(pt1**(-2.78))/(1.0+0.4866*(pt1**(-2.51)))
         jlisa2=0.1425*(pt2**(-2.78))/(1.0+0.4866*(pt2**(-2.51)))

         if(tvalue(t)<tfora)then
          parta1=tfora*(tfora+2.0*0.938)/(tnew1*(tnew1+2.0*0.938))
         else
          parta1=tvalue(t)*(tvalue(t)+2.0*0.938)
     &                   /(tnew1*(tnew1+2.0*0.938))
         endif
         parta2=tvalue(t+1)*(tvalue(t+1)+2.0*0.938)
     &                    /(tnew2*(tnew2+2.0*0.938))

         if(tvalue(t)<tfora)then
          qfora=qfora+0.5*(jlisa2*parta2*yxt_aad(t+1)+jlisa1*parta1*
     &    (yxt_aad(t)+(tfora-tvalue(t))*(yxt_aad(t+1)-yxt_aad(t))
     &      /(tvalue(t+1)-tvalue(t))))*(tvalue(t+1)-tfora)
         else
          qfora=qfora+0.5*(jlisa2*parta2*yxt_aad(t+1)+jlisa1
     &      *parta1*yxt_aad(t))*(tvalue(t+1)-tvalue(t))
         endif
       
        endif

       enddo

!***********************************************************************
       qout=(qforp+qfora)*airout !GCR induced ionization (ion pairs cm-3 sec-1)
!qout=(qforp+qfora) !cosmic ray induced ionization (ion pairs g-1 sec-1)


       end subroutine  IONGCR
!------------------------------------------------------------------------------

!
! APM2+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

! Pure organic nucl based on Kirkby et al., 2016
!
      subroutine PORGANIC(XQ,XCS,CHOM,XJN,XJI)
      REAL*8 :: XQ,XCS,CHOM,XJN,XJI

      REAL*8, PARAMETER   :: a1 = 0.04001
      REAL*8, PARAMETER   :: a2 = 1.848
      REAL*8, PARAMETER   :: a3 = 0.001366
      REAL*8, PARAMETER   :: a4 = 1.566
      REAL*8, PARAMETER   :: a5 = 0.1863
      REAL*8, PARAMETER   :: XREC = 1.6d-6   ! recomb coef

      REAL*8 :: XHOM,XTEMP, XKI, XNI

! XHOM in the unit of 1E7/cm3
      XHOM=CHOM/1.d7
      XJN = a1*XHOM**(a2+a5/XHOM)
      XTEMP = a3*XHOM**(a4+a5/XHOM)
      XKI = XCS + XTEMP
      XNI = (SQRT(XKI*XKI+4.*XREC*XQ) - XKI)/(2.*XREC)
      XJI = 2*XNI*XTEMP

!       WRITE(6,100)XHOM,XKI,XNI,XJN,XJI,XJN+XJI
! 100  FORMAT(10(1PE9.2))

      RETURN
      END SUBROUTINE PORGANIC 

!------------------------------------------------------------------------------

      SUBROUTINE BHN_VEHKAMAKI(t,rh,rhoa,jnuc,x,ntot,rc)

c    Fortran 90 subroutine BHN_VEHKAMAKI
c
c    Calculates parametrized values of nucleation rate,
c    mole fraction of sulphuric acid
c    total number of particles, and the radius of the critical cluster
c    in H2O-H2SO4 system if temperature, saturatio ratio of water and 
c    sulfuric acid concentration  are given. 
c
c    Copyright (C) 2002 Hanna Vehkamäki
c
c    Division of Atmospheric Sciences
c    Department of Physical Sciences
c    P.O. Box 64
c    FIN-00014 University of Helsinki
c    Finland
c    
c    hanna.vehkamaki@helsinki.fi


      IMPLICIT NONE 

      real*8 ::t,rh,rhoa,jnuc,x,ntot,rc,rhotres
                                
c     t: temperature in K (190.15-300.15K)                                  
c     rh: saturatio ratio of water (0.0001-1)                               
c     rhoa: sulfuric acid concentration in 1/cm3 (10^4-10^11 1/cm3)         
c     jnuc: nucleation rate in 1/cm3s (10^-7-10^10 1/cm3s)                  
c     ntot: total number of molecules in the critical cluster (ntot>4)      
c     x: molefraction of H2SO4 in the critical cluster                      
c     rc: radius of the critical cluster in nm                              
c     rhotres: treshold concentration of h2so4 (1/cm^3)                     
c     which produces nucleation rate   1/(cm^3 s) as a function of rh and t 


      if(t .lt. 190.15) then
!         print *,'BHN Warning: temperature < 190.15 K, using 190.15 K'
         t=190.15
      end if

      if(t .gt. 300.15) then
!         print *,'BHN Warning: temperature > 300.15 K, using 300.15 K'
         t=300.15
      end if

      if(rh .lt. 0.0001) then
!         print *,'BHN Warning: saturation ratio of water < 0.0001,',
!     &           ' using 0.0001'
         rh=0.0001
      end if

      if(rh .gt. 1.0) then
!         print *,'BHN Warning: saturation ratioof water > 1 using 1'
         rh=1.0
      end if

      if(rhoa .lt. 1.e4) then
!         print *,'BHN Warning: sulfuric acid concentration <',
!     &      '1e4 1/cm3, using 1e4 1/cm3'
         rhoa=1.e4
      end if

      if(rhoa .gt. 1.e11) then
!         print *,'BHN Warning: sulfuric acid concentration',
!     &    ' > 1e11 1/cm3, using 1e11 1/cm3'
CY         rhoa=1.e11
      end if



      x=  0.7409967177282139 - 0.002663785665140117*t 
     &     + 0.002010478847383187*Log(rh) 
     &     - 0.0001832894131464668*t*Log(rh) 
     &     + 0.001574072538464286*Log(rh)**2       
     &     - 0.00001790589121766952*t*Log(rh)**2 
     &     + 0.0001844027436573778*Log(rh)**3  
     &     -  1.503452308794887e-6*t*Log(rh)**3 
     &     - 0.003499978417957668*Log(rhoa)
     &     + 0.0000504021689382576*t*Log(rhoa)


      jnuc= 0.1430901615568665 + 2.219563673425199*t -
     &     0.02739106114964264*t**2 +  
     &     0.00007228107239317088*t**3 + 5.91822263375044/x +  
     &     0.1174886643003278*Log(rh) + 0.4625315047693772*t*Log(rh) -  
     &     0.01180591129059253*t**2*Log(rh) +  
     &     0.0000404196487152575*t**3*Log(rh) + 
     &     (15.79628615047088*Log(rh))/x -  
     &     0.215553951893509*Log(rh)**2 - 
     &     0.0810269192332194*t*Log(rh)**2 +  
     &     0.001435808434184642*t**2*Log(rh)**2 -  
     &     4.775796947178588e-6*t**3*Log(rh)**2 -  
     &     (2.912974063702185*Log(rh)**2)/x -
     &     3.588557942822751*Log(rh)**3 +  
     &     0.04950795302831703*t*Log(rh)**3 -  
     &     0.0002138195118737068*t**2*Log(rh)**3 +  
     &     3.108005107949533e-7*t**3*Log(rh)**3 -  
     &     (0.02933332747098296*Log(rh)**3)/x +  
     &     1.145983818561277*Log(rhoa) - 
     &     0.6007956227856778*t*Log(rhoa) + 
     &     0.00864244733283759*t**2*Log(rhoa) - 
     &     0.00002289467254710888*t**3*Log(rhoa) - 
     &     (8.44984513869014*Log(rhoa))/x + 
     &     2.158548369286559*Log(rh)*Log(rhoa) + 
     &     0.0808121412840917*t*Log(rh)*Log(rhoa) - 
     &     0.0004073815255395214*t**2*Log(rh)*Log(rhoa) - 
     &     4.019572560156515e-7*t**3*Log(rh)*Log(rhoa) + 
     &     (0.7213255852557236*Log(rh)*Log(rhoa))/x + 
     &     1.62409850488771*Log(rh)**2*Log(rhoa) - 
     &     0.01601062035325362*t*Log(rh)**2*Log(rhoa) + 
     &     0.00003771238979714162*t**2*Log(rh)**2*Log(rhoa) + 
     &     3.217942606371182e-8*t**3*Log(rh)**2*Log(rhoa) - 
     &     (0.01132550810022116*Log(rh)**2*Log(rhoa))/x + 
     &     9.71681713056504*Log(rhoa)**2 - 
     &     0.1150478558347306*t*Log(rhoa)**2 + 
     &     0.0001570982486038294*t**2*Log(rhoa)**2 + 
     &     4.009144680125015e-7*t**3*Log(rhoa)**2 + 
     &     (0.7118597859976135*Log(rhoa)**2)/x - 
     &     1.056105824379897*Log(rh)*Log(rhoa)**2 + 
     &     0.00903377584628419*t*Log(rh)*Log(rhoa)**2 - 
     &     0.00001984167387090606*t**2*Log(rh)*Log(rhoa)**2 + 
     &     2.460478196482179e-8*t**3*Log(rh)*Log(rhoa)**2 - 
     &     (0.05790872906645181*Log(rh)*Log(rhoa)**2)/x - 

     &     0.1487119673397459*Log(rhoa)**3 + 
     &     0.002835082097822667*t*Log(rhoa)**3 - 
     &     9.24618825471694e-6*t**2*Log(rhoa)**3 + 
     &     5.004267665960894e-9*t**3*Log(rhoa)**3 - 
     &     (0.01270805101481648*Log(rhoa)**3)/x
      jnuc=exp(jnuc)            
c     1/(cm3s)

      ntot =-0.002954125078716302 - 0.0976834264241286*t +
     &     0.001024847927067835*t**2 - 2.186459697726116e-6*t**3 - 
     &     0.1017165718716887/x - 0.002050640345231486*Log(rh) -
     &     0.007585041382707174*t*Log(rh) + 
     &     0.0001926539658089536*t**2*Log(rh) -
     &     6.70429719683894e-7*t**3*Log(rh) - 
     &     (0.2557744774673163*Log(rh))/x +
     &     0.003223076552477191*Log(rh)**2 +
     &     0.000852636632240633*t*Log(rh)**2 - 
     &     0.00001547571354871789*t**2*Log(rh)**2 +
     &     5.666608424980593e-8*t**3*Log(rh)**2 + 
     &     (0.03384437400744206*Log(rh)**2)/x +
     &     0.04743226764572505*Log(rh)**3 - 
     &     0.0006251042204583412*t*Log(rh)**3 +
     &     2.650663328519478e-6*t**2*Log(rh)**3 - 
     &     3.674710848763778e-9*t**3*Log(rh)**3 -
     &     (0.0002672510825259393*Log(rh)**3)/x - 
     &     0.01252108546759328*Log(rhoa) +
     &     0.005806550506277202*t*Log(rhoa) - 
     &     0.0001016735312443444*t**2*Log(rhoa) +
     &     2.881946187214505e-7*t**3*Log(rhoa) + 
     &     (0.0942243379396279*Log(rhoa))/x -
     &     0.0385459592773097*Log(rh)*Log(rhoa) - 
     &     0.0006723156277391984*t*Log(rh)*Log(rhoa) +
     &     2.602884877659698e-6*t**2*Log(rh)*Log(rhoa) + 
     &     1.194163699688297e-8*t**3*Log(rh)*Log(rhoa) -
     &     (0.00851515345806281*Log(rh)*Log(rhoa))/x - 
     &     0.01837488495738111*Log(rh)**2*Log(rhoa) +
     &     0.0001720723574407498*t*Log(rh)**2*Log(rhoa) - 
     &     3.717657974086814e-7*t**2*Log(rh)**2*Log(rhoa) - 
     &     5.148746022615196e-10*t**3*Log(rh)**2*Log(rhoa) + 
     &     (0.0002686602132926594*Log(rh)**2*Log(rhoa))/x -
     &     0.06199739728812199*Log(rhoa)**2 + 
     &     0.000906958053583576*t*Log(rhoa)**2 -
     &     9.11727926129757e-7*t**2*Log(rhoa)**2 - 
     &     5.367963396508457e-9*t**3*Log(rhoa)**2 -
     &     (0.007742343393937707*Log(rhoa)**2)/x + 
     &     0.0121827103101659*Log(rh)*Log(rhoa)**2 -
     &     0.0001066499571188091*t*Log(rh)*Log(rhoa)**2 + 
     &     2.534598655067518e-7*t**2*Log(rh)*Log(rhoa)**2 - 
     &     3.635186504599571e-10*t**3*Log(rh)*Log(rhoa)**2 + 
     &     (0.0006100650851863252*Log(rh)*Log(rhoa)**2)/x +
     &     0.0003201836700403512*Log(rhoa)**3 - 
     &     0.0000174761713262546*t*Log(rhoa)**3 +
     &     6.065037668052182e-8*t**2*Log(rhoa)**3 - 
     &     1.421771723004557e-11*t**3*Log(rhoa)**3 +
     &     (0.0001357509859501723*Log(rhoa)**3)/x
      ntot=exp(ntot)

      rc=exp(-1.6524245+0.42316402*x+0.33466487*log(ntot)) 
C     nm


      if(ntot .lt. 4.0 ) then
         if ( t .lt. 195.15) then
C            print *,'Warning: critical cluster < 4 molecules',
C     &              ' and T<195.15 : rough estimate J=1e5/(cm3s)'
Cyu            jnuc=1.e5
            jnuc=1.e0
         else
C            print *,'Warning: critical cluster < 4 molecules:',
C     &              ' accuracy of the parameterisation reduced,',  
C     &              ' J_theor/J_para up to 50'    
         end if
      end if

      if(jnuc .lt. 1.e-7) then
C         print *,'Warning: nucleation rate < 1e-7/cm3s,', 
C     &           ' using 1.0E-7/cm3s'

         jnuc=1.0E-7
      end if

      if(jnuc .gt. 1.e10) then
C         print *,'Warning: nucleation rate > 1e10/cm3s,',
C     &           ' using 1e10/cm3s'
         jnuc=1.e10
      end if



      rhotres=exp( -279.2430007512709 + 11.73439886096903*rh
     &     + 22700.92970508331/t
     &     - (1088.644983466801*rh)/t + 1.144362942094912*t

     &     - 0.03023314602163684*rh*t - 0.001302541390154324*t**2  
     &     - 6.386965238433532*Log(rh) + (854.980361026715*Log(rh))/t 
     &     + 0.00879662256826497*t*Log(rh)) 
c     1/cm3

      RETURN

      END ! SUBROUTINE BHN_VEHKAMAKI
!
!------------------------------------------------------------------------------
      SUBROUTINE BHN_KULMALA98(temp,rh,h2so4,JBHN_Kul)
      
! *** subroutine to implement Kulmala, Laaksonen, Pirjola
!     Parameterization for sulfuric acid/water
!     nucleation rates, J. Geophys. Research (103), pp 8301-8307,
!     April 20, 1998.

! Modified from the code in WRF-Chem (fyu, 12/10/2012)

       IMPLICIT NONE
       
! *** Input:

! ambient temperature [ K ]
      REAL*8 temp
! fractional relative humidity
      REAL*8 rh
! sulfuric acid concentration [ # / cm**3 ]
      REAL*8 h2so4
      
! *** Output:
      REAL*8 JBHN_Kul
      

! *** NOTE, all units are cgs internally.

      REAL*8 ra
! fractional relative acidity
! sulfuric acid vaper concentration [ cm ** -3 ]
      REAL*8 nav
! water vapor concentration   [ cm ** -3 ]
      REAL*8 nwv
! equilibrium sulfuric acid vapor conc. [ cm ** -3 ]
      REAL*8 nav0
                ! to produce a nucleation rate of 1 [ cm ** -3  s ** -1
      REAL*8 nac  
! critical sulfuric acid vapor concentration [ cm ** -3
! mole fractio of the critical nucleus                 
      REAL*8 xal 
      REAL*8 nsulf, delta
      REAL*8 chi, jnuc
! nucleation rate [ cm ** -3  s ** -1 ]
      REAL*8 tt, rr  ! dummy variables for statement functions
      
! avogadro's constant [ 1/mol ]
      REAL*8 avo
      PARAMETER (avo=6.0221367E23)
      
! universal gas constant [ j/mol-k ]
      REAL*8 rgasuniv
      PARAMETER (rgasuniv=8.314510)
      
! 1 atmosphere in pascals
      REAL*8 atm 
      PARAMETER (atm=1013.25E+02)
      
! *** statement functions **************
! arithmetic statement function to compute

      REAL*8 ph2so4, ph2o                  ! for h2so4 and h2o vapor pressures [ Pa ]


! *** both ph2o and ph2so4 are  as in Kulmala et al.  paper
      ph2o(tt)=exp(77.34491296-7235.4246512/tt-8.2*log(tt)
     &                       +tt*5.7113E-03)
      ph2so4(tt)=exp(27.78492066-10156.0/tt)

! *** begin code:

!     The 1.0e-6 factor in the following converts from MKS to cgs units
! *** get water vapor concentration [ molecles / cm **3 ]
      nwv = rh*ph2o(temp)/(rgasuniv*temp)*avo*1.0E-6

! *** calculate the equilibrium h2so4 vapor concentration.
! *** use Kulmala corrections:
! ***
      nav0 = ph2so4(temp)/(rgasuniv*temp)*avo*1.0E-6

      nav = h2so4

! *** calculate critical concentration of sulfuric acid vapor

      nac = exp(-14.5125+0.1335*temp-10.5462*rh+1958.4*rh/temp)

! *** calculate relative acidity

      ra = nav/nav0

! *** calculate temperature correction

      delta = 1.0 + (temp-273.15)/273.14

! *** calculate molar fraction

      xal = 1.2233 - 0.0154*ra/(ra+rh) + 0.0102*log(nav)
     &  - 0.0415*log(nwv) +  0.0016*temp

! *** calculate Nsulf
      nsulf = log(nav/nac)

! *** calculate particle produtcion rate [ # / cm**3 ]

      chi = 25.1289*nsulf - 4890.8*nsulf/temp - 1743.3/temp -
     & 2.2479*delta*nsulf*rh + 7643.4*xal/temp - 1.9712*xal*delta/rh

      jnuc = exp(chi)

! [ # cm-3s-1 ]
      JBHN_Kul = jnuc

      RETURN

      END SUBROUTINE BHN_KULMALA98


!------------------------------------------------------------------------------
C      subroutine terpara77(t,rh,c2,c3,j,ntot,n2,n3,r)     
      subroutine terpara77(TK,YRH,CACID,CNH3,j)     

c     Fortran 90 subroutine terpara
c
c     Calculates parametrized values of nucleation rate,
c     total number of particles, number of sulfuric acid molecules,
c     number of ammonia molecules, and the radius of the critical cluster
c     in H2O-H2SO4-NH3 system if temperature, saturation ratio of water, 
c     sulfuric acid concentration, and ammonia mixing ratio are given. 
c
c     Copyright (C) 2002 Ismo Napari
c
c     Division of Atmospheric Sciences
c     Department of Physical Sciences
c     P.O. Box 64
c     FIN-00014 University of Helsinki
c     Finland
c    
c     ismo.napari@helsinki.fi
 
      implicit none

      real*8 TK,t    
      real*8 YRH,rh    
      real*8 CACID,c2    
      real*8 CNH3,c3    
      real*8 j    

      real*8 ntot 
      real*8 n2   
      real*8 n3   
      real*8 r          
      real*8 lnj  

c     t:    temperature in K (240-300K)
c     rh:   saturation ratio of water(0.05-0.95)
c     c2:   sulfuric acid concentration in 1/cm3 (10^4-10^9 1/cm3)
c     c3:   ammonia mixing ratio in ppt (volume) (0.1-100 ppt)
c     j:    nucleation rate in 1/cm3s (10^-5-10^6 1/cm3s)
c     ntot: total number of particles
c     n2:   number of sulfuric acid molecules
c     n3:   number of ammonia molecules
c     r:    radius of the critical cluster in nm      

      t=TK
      rh = YRH
      c2 = CACID
      c3 = CNH3

      if(t .lt. 240.) then
C         print *,'Warning: temperature < 240 K, using 240 K'
         t=240.

      end if

      if(t .gt. 300.) then
C         print *,'Warning: temperature > 300 K, using 300 K'
         t=300.
      end if

      if(rh .lt. 0.05) then
C         print *,'Warning: saturation ratio of water < 0.05,',
C     &           ' using 0.05'
         rh=0.05
      end if

      if(rh .gt. 0.95) then
C         print *,'Warning: saturation ratio of water > 0.95,',
C     &           ' using 0.95'
         rh=0.95
      end if

      if(c2 .lt. 1.e4) then
C         print *,'Warning: sulfuric acid concentration < 1e4 1/cm3,',
C     &           ' using 1e4 1/cm3'
         c2=1.e4
      end if

      if(c2 .gt. 1.e9) then
C         print *,'Warning: sulfuric acid concentration > 1e9 1/cm3,',
C     &           ' using 1e9 1/cm3'

         c2=1.e9
      end if

      if(c3 .lt. 0.1) then
C         print *,'Warning: ammonia mixing ratio < 0.1 ppt(volume),',
C     &           ' using 0.1 ppt(volume)'
         c3=0.1
      end if

      if(c3 .gt. 100.) then
C         print *,'Warning: ammonia mixing ratio > 100 ppt(volume),',
C     &           ' using 100 ppt(volume)'
         c3=100.
      end if

      lnj = -84.7551114741543 + 0.3117595133628944*rh +
     &       1.640089605712946*rh*t - 
     &       0.003438516933381083*rh*t**2 -
     &       0.00001097530402419113*rh*t**3 - 
     &       0.3552967070274677/Log(c2) -
     &       (0.06651397829765026*rh)/Log(c2) - 
     &       (33.84493989762471*t)/Log(c2) -
     &       (7.823815852128623*rh*t)/Log(c2) + 
     &       (0.3453602302090915*t**2)/Log(c2) + 
     &       (0.01229375748100015*rh*t**2)/Log(c2) - 
     &       (0.000824007160514956*t**3)/Log(c2) + 
     &       (0.00006185539100670249*rh*t**3)/Log(c2) + 
     &       3.137345238574998*Log(c2) +
     &       3.680240980277051*rh*Log(c2) - 
     &       0.7728606202085936*t*Log(c2) -
     &       0.204098217156962*rh*t*Log(c2) + 
     &       0.005612037586790018*t**2*Log(c2) + 
     &       0.001062588391907444*rh*t**2*Log(c2) - 
     &       9.74575691760229e-6*t**3*Log(c2) - 
     &       1.265595265137352e-6*rh*t**3*Log(c2) +
     &       19.03593713032114*Log(c2)**2 - 
     &       0.1709570721236754*t*Log(c2)**2 + 
     &       0.000479808018162089*t**2*Log(c2)**2 - 
     &       4.146989369117246e-7*t**3*Log(c2)**2 +
     &       1.076046750412183*Log(c3) + 
     &       0.6587399318567337*rh*Log(c3) +
     &       1.48932164750748*t*Log(c3) + 
     &       0.1905424394695381*rh*t*Log(c3) -
     &       0.007960522921316015*t**2*Log(c3) - 
     &       0.001657184248661241*rh*t**2*Log(c3) + 
     &       7.612287245047392e-6*t**3*Log(c3) + 
     &       3.417436525881869e-6*rh*t**3*Log(c3) +  
     &       (0.1655358260404061*Log(c3))/Log(c2) +  
     &       (0.05301667612522116*rh*Log(c3))/Log(c2) + 
     &       (3.26622914116752*t*Log(c3))/Log(c2) - 
     &       (1.988145079742164*rh*t*Log(c3))/Log(c2) - 
     &       (0.04897027401984064*t**2*Log(c3))/Log(c2) + 
     &       (0.01578269253599732*rh*t**2*Log(c3))/Log(c2) + 
     &       (0.0001469672236351303*t**3*Log(c3))/Log(c2) - 
     &       (0.00002935642836387197*rh*t**3*Log(c3))/Log(c2) + 
     &       6.526451177887659*Log(c2)*Log(c3) -  
     &       0.2580021816722099*t*Log(c2)*Log(c3) + 
     &       0.001434563104474292*t**2*Log(c2)*Log(c3) - 
     &       2.020361939304473e-6*t**3*Log(c2)*Log(c3) - 
     &       0.160335824596627*Log(c2)**2*Log(c3) + 
     &       0.00889880721460806*t*Log(c2)**2*Log(c3) - 
     &       0.00005395139051155007*t**2*Log(c2)**2*Log(c3) + 
     &       8.39521718689596e-8*t**3*Log(c2)**2*Log(c3) + 
     &       6.091597586754857*Log(c3)**2 +
     &       8.5786763679309*rh*Log(c3)**2 - 
     &       1.253783854872055*t*Log(c3)**2 -
     &       0.1123577232346848*rh*t*Log(c3)**2 + 
     &       0.00939835595219825*t**2*Log(c3)**2 + 
     &       0.0004726256283031513*rh*t**2*Log(c3)**2 - 
     &       0.00001749269360523252*t**3*Log(c3)**2 - 
     &       6.483647863710339e-7*rh*t**3*Log(c3)**2 + 
     &       (0.7284285726576598*Log(c3)**2)/Log(c2) + 
     &       (3.647355600846383*t*Log(c3)**2)/Log(c2) - 
     &       (0.02742195276078021*t**2*Log(c3)**2)/Log(c2) + 
     &       (0.00004934777934047135*t**3*Log(c3)**2)/Log(c2) + 
     &       41.30162491567873*Log(c2)*Log(c3)**2 - 
     &       0.357520416800604*t*Log(c2)*Log(c3)**2 + 
     &       0.000904383005178356*t**2*Log(c2)*Log(c3)**2 - 
     &       5.737876676408978e-7*t**3*Log(c2)*Log(c3)**2 - 
     &       2.327363918851818*Log(c2)**2*Log(c3)**2 + 
     &       0.02346464261919324*t*Log(c2)**2*Log(c3)**2 - 
     &       0.000076518969516405*t**2*Log(c2)**2*Log(c3)**2 + 
     &       8.04589834836395e-8*t**3*Log(c2)**2*Log(c3)**2 - 
     &       0.02007379204248076*Log(rh) -
     &       0.7521152446208771*t*Log(rh) + 
     &       0.005258130151226247*t**2*Log(rh) - 
     &       8.98037634284419e-6*t**3*Log(rh) + 
     &       (0.05993213079516759*Log(rh))/Log(c2) + 
     &       (5.964746463184173*t*Log(rh))/Log(c2) - 
     &       (0.03624322255690942*t**2*Log(rh))/Log(c2) + 
     &       (0.00004933369382462509*t**3*Log(rh))/Log(c2) - 
     &       0.7327310805365114*Log(c3)*Log(rh) - 
     &       0.01841792282958795*t*Log(c3)*Log(rh) + 
     &       0.0001471855981005184*t**2*Log(c3)*Log(rh) - 
     &       2.377113195631848e-7*t**3*Log(c3)*Log(rh)

      j = exp(lnj)

      if(j .lt. 1.e-5) then
C         print *,'Warning: nucleation rate < 1e-5/cm3s,',
C     &           ' using 0.0/cm3s,'
         j=0.0
      end if
 
      if(j .gt. 1.e6) then
C         print *,'Warning: nucleation rate > 1e6/cm3s, using 1e6/cm3s'
         j=1.e6
      end if

      ntot = 79.34842673023872 + 1.738399264246183*lnj + 
     &       0.007114029786143758*lnj**2 - 0.7449934047779641*t - 
     &       0.00820608264651173*lnj*t + 0.001785497681201935*t**2

      n2 = 38.16448247950508 + 0.7741058259731187*lnj + 
     &     0.002988789927230632*lnj**2 - 0.3576046920535017*t - 
     &     0.003663583011953248*lnj*t + 0.000855300153372776*t**2

      n3 = 26.89821829677973 + 0.6829050502793024*lnj + 
     &     0.003575211441884629*lnj**2 - 0.2657482507363107*t - 
     &     0.003418950836780133*lnj*t + 0.000673454337882302*t**2

      r = 0.1410271086638381 - 0.001226253898894878*lnj - 
     &    7.822111731550752e-6*lnj**2 - 0.001567273351921166*t - 
     &    0.00003075996088273962*lnj*t + 0.00001083754117202233*t**2

      return 

      end subroutine terpara77

!------------------------------------------------------------------------------

C      SUBROUTINE LovejoyIN(H2SO4, T, RH, Q, SA, H1, H3, H4, H5)
      SUBROUTINE LovejoyIN(CH2SO4, TK, RH0, Q0, SA0, H1, H3, H4, H5)
C 
C INPUT: Temprature(T, K), Relative Humidity (RH), Sulfuric acid concentration (C, cm-3), ion 
C production rate (q, cm-3 s-1) and  Pre-existing Aerosol Surface Area (SA, micro m2 cm-3)
C
C Valid range: T: 190-300K; RH: 5-95; C:1E5-1E8/cm3; SA: 2-100 um2/cm3; Q: 1-50 pairs/cm3s
C
C Parameterization from Modgil, M. S., S. Kumar, S. N. Tripathi, and E. R. Lovejoy (2005),
C A parameterization of ion-induced nucleation of sulphuric acid and water for atmospheric
C conditions, J. Geophys. Res., 110, D19205, doi:10.1029/2004JD005475.

C Mathematical Expressions for particle nucleation rate (h1,cm-3 s-1), nucleating
C H2SO4 flux (h2, cm-3 s-1), number of H2SO4 in average nucleating cluster (h3),
C  number of H2O in average nucleating cluster (h4), radius of average nucleating
C cluster (h5, nm), and first order loss of H2SO4 to particles (h6) are given below,
C *********************************************************************
C
C Temprature(T, K), Relative Humidity (RH), Sulfuric acid concentration (H2SO4, cm-3), ion
C production rate (q, cm-3 s-1) and  Pre-existing Aerosol Surface Area (SA, micro m2 cm-3)
C

C	double precision t,rh,h2so4,q,sa,h1,h2,h3,h4,h6,h5
C        double precision CH2SO4, TK, RH0, Q0, SA0
	REAL*8 :: t,rh,h2so4,q,sa,h1,h2,h3,h4,h6,h5
        REAL*8 :: CH2SO4, TK, RH0, Q0, SA0
C
        h2so4 = CH2SO4
        t= TK
        rh = RH0
        q = Q0
        sa = SA0
C
        IF(h2so4.LT.1.E5) THEN
C           WRITE(6,10) h2so4
C           h2so4 = 1.E5
            H1=1.E-20
            H3=78.
            H4=200.
            H5=1.5
            RETURN
        ENDIF

        IF(h2so4.GT.1.E8) THEN
C           WRITE(6,11) h2so4
           h2so4 = 1.E8
        ENDIF

        IF(t.LT.190.) THEN
C           WRITE(6,12) t
           t = 190.
        ENDIF

        IF(t.GT.300.) THEN
C           WRITE(6,13) t
           t = 300.
            H1=1.E-20
            H3=78.
            H4=200.
            H5=1.5
            RETURN
        ENDIF

        IF(rh.LT.5.) THEN
C           WRITE(6,14) rh
C           rh = 5.
            H1=1.E-20
            H3=78.
            H4=200.
            H5=1.5
            RETURN
        ENDIF

        IF(rh.GT.95.) THEN
C           WRITE(6,15) rh
           rh = 95.
        ENDIF

        IF(q.LT.1.) THEN
C           WRITE(6,16) 1
C           q = 1.
            H1=1.E-20
            H3=78.
            H4=200.
            H5=1.5
            RETURN
        ENDIF

        IF(q.GT.50.) THEN
C           WRITE(6,17) 1
           q = 50.
        ENDIF

        IF(sa.LT.2.) THEN
C           WRITE(6,18) sa
           sa = 2.
        ENDIF

        IF(sa.GT.100.) THEN
           WRITE(6,19) sa
           sa = 100.
C            H1=1.E-20
C            H3=78.
C            H4=200.
C            H5=1.5
C            RETURN
        ENDIF

 10     FORMAT("IIN WARNING: INPUTED [H2SO4]=",ES9.2,"< 1.E5",
     &     " set it to 1.E5")
 11     FORMAT("IIN WARNING: INPUTED [H2SO4]=",ES9.2,"> 1.E8",
     &     " set it to 1.E8")
 12     FORMAT("IIN WARNING: INPUTED T=",ES9.2,"< 190.",
     &     " set it to 190")
 13     FORMAT("IIN WARNING: INPUTED T=",ES9.2,"> 300.",
     &     " set it to 300.")
 14     FORMAT("IIN WARNING: INPUTED RH=",ES9.2,"< 5",
     &     " set it to 5")
 15     FORMAT("IIN WARNING: INPUTED RH=",ES9.2,"> 95.",
     &     " set it to 95.")
 16     FORMAT("IIN WARNING: INPUTED Q=",ES9.2,"< 1",
     &     " set it to 1")
 17     FORMAT("IIN WARNING: INPUTED Q=",ES9.2,"> 50.",
     &     " set it to 50.")
 18     FORMAT("IIN WARNING: INPUTED SA=",ES9.2,"< 2",
     &     " set it to 1")
 19     FORMAT("IIN WARNING: INPUTED SA=",ES9.2,"> 100.",
     &     " set it to 50.")

C inputed RH in percentage (5-95), rh used in the parameterization in fraction (0.05-0.95)
C
        rh = RH/100.	

	h6=0.000026859579119003205*SA + 1.7477354270484002e-8*q*SA + 
     &   1.5718068902491457e-8*SA**2 + 8.060796806911441e-8*SA*T + 
     &   3.904048293417882e-7*SA*Log(H2SO4) + 
     &   2.727259306977938e-7*SA*Log(RH)
	

  	h3=-198.8039518313554 + 3357.132963009284*h6 - 130510.31325149858*
     &  h6**2 - 0.7093715033997716*q - 10.713505046150196*h6*q + 
     &  1103.4737682776713*h6**2*q + 0.0052565148186649*q**2 - 
     &  0.20195850414426988*h6*q**2 + 10.961027676935213*h6**2*q**2 - 
     &  26.553841269634976*RH + 2913.499196899548*h6*RH - 
     &  7558.996305824136*h6**2*RH + 0.050092880471591994*q*RH + 
     &  0.39840936335061017*h6*q*RH + 16.140386509938388*h6**2*q*RH - 
     &  0.0008159572217147427*q**2*RH - 0.02492462618304389*h6*q**2*RH + 
     &  3.2372842210428825*RH**2 + 1709.7485838150235*h6*RH**2 - 
     &  4016.182638678486*h6**2*RH**2 - 0.022142010235491123*q*RH**2 - 
     &  1.620063009925805*h6*q*RH**2-0.00028477984814528825*q**2*RH**2 + 
     &  22.136724153656015*RH**3 + 170.8982375938333*h6*RH**3 - 
     &  0.01881686723215867*q*RH**3 + 2.6974144456100957*T - 
     &  96.60591604496483*h6*T + 1772.137264721083*h6**2*T - 
     &  0.0009432251807207652*q*T + 0.06072064184950673*h6*q*T - 
     &  2.5196932894429502*h6**2*q*T - 0.000013848768113392552*q**2*T - 
     &  0.0001948394841164792*h6*q**2*T + 0.1828636512279507*RH*T - 
     &  55.135341874839185*h6*RH*T + 164.02631709083576*h6**2*RH*T + 
     &  0.001745921048607296*q*RH*T + 0.035017713828742754*h6*q*RH*T + 
     &  4.057082638293583e-6*q**2*RH*T - 0.3900441693913758*RH**2*T - 
     &  8.955078982582657*h6*RH**2*T+0.00021434974336412074*q*RH**2*T - 
     &  0.14947568974964962*RH**3*T - 0.022748394377623382*T**2 + 
     &  0.7227721843282815*h6*T**2 - 5.386480671871808*h6**2*T**2 - 
     &  0.000035250836279611095*q*T**2-0.0003363405774846326*h6*q*T**2 + 
     &  2.9254973516794257e-8*q**2*T**2 + 0.003994529829421164*RH*T**2 + 
     &  0.2074067980035454*h6*RH*T**2 - 5.136172472264946e-6*q*RH*T**2 + 
     &  0.0020603328018819816*RH**2*T**2+0.000042019279193164354*T**3 - 
     &  0.002100661388749787*h6*T**3 + 7.309966632740304e-8*q*T**3 - 
     &  0.000016323969052607556*RH*T**3+ 12.330627568298462*Log(H2SO4) + 
     &  768.3961789008589*h6*Log(H2SO4)-11568.47324943553*h6**2*
     &  Log(H2SO4)+0.14349043416922366*q*Log(H2SO4)+0.8946851157223353*
     &  h6*q*Log(H2SO4) - 68.46004143191098*h6**2*q*Log(H2SO4) - 
     &  0.0006241121793370407*q**2*Log(H2SO4)+0.011897674833907721*h6*
     &  q**2*Log(H2SO4) + 1.7860574934328677*RH*Log(H2SO4) + 
     &  316.5406316191978*h6*RH*Log(H2SO4) - 2036.825340216443*h6**2*RH*
     &  Log(H2SO4) - 0.026323914507434605*q*RH*Log(H2SO4) - 
     &  0.37505804775954393*h6*q*RH*Log(H2SO4) + 0.00003454867680790666*
     &  q**2*RH*Log(H2SO4) + 2.844877302874606*RH**2*Log(H2SO4) - 
     &  1.5845895178086176*h6*RH**2*Log(H2SO4) + 0.001732608008714275*q*
     &  RH**2*Log(H2SO4) + 0.5611003862827533*RH**3*Log(H2SO4) + 
     &  0.18033151281768975*T*Log(H2SO4) - 7.807090214680351*h6*T*
     &  Log(H2SO4) + 52.76241348342321*h6**2*T*Log(H2SO4) + 
     &  0.0011535134888316242*q*T*Log(H2SO4) + 0.0068466874844708295*h6*
     &  q*T*Log(H2SO4) - 4.231224168766194e-7*q**2*T*Log(H2SO4) - 
     &  0.10505349775719895*RH*T*Log(H2SO4) - 1.9241106452950727*h6*RH*
     &  T*Log(H2SO4) + 2.0451815715440337e-6*q*RH*T*Log(H2SO4) - 
     &  0.015299483302534183*RH**2*T*Log(H2SO4) + 0.000775633115370002*
     &  T**2*Log(H2SO4) + 0.04228608723566267*h6*T**2*Log(H2SO4) - 
     &  9.572221299803945e-7*q*T**2*Log(H2SO4) + 0.0002669785990812474*
     &  RH*T**2*Log(H2SO4) - 1.7595742533055222e-6*T**3*Log(H2SO4) - 
     &  2.7165200489046812*Log(H2SO4)**2 + 1.963036665672805*h6*
     &  Log(H2SO4)**2 + 33.88545004559797*h6**2*Log(H2SO4)**2 - 
     &  0.01647722099703982*q*Log(H2SO4)**2 - 0.08498050322218324*h6*q*
     &  Log(H2SO4)**2 + 0.00003320475358154802*q**2*Log(H2SO4)**2 + 
     &  0.50626436977659*RH*Log(H2SO4)**2 + 3.914586690682404*h6*RH*
     &  Log(H2SO4)**2 + 0.0006654515705980484*q*RH*Log(H2SO4)**2 - 
     &  0.025122425208873058*RH**2*Log(H2SO4)**2 - 0.021380868797539664*
     &  T*Log(H2SO4)**2 - 0.35405514523597137*h6*T*Log(H2SO4)**2 - 
     &  0.000020639290758942666*q*T*Log(H2SO4)**2+0.0003587951811915662*
     &  RH*T*Log(H2SO4)**2 + 7.620708111644729e-6*T**2*Log(H2SO4)**2 + 
     &  0.2196641696573127*Log(H2SO4)**3 + 1.7291708055805226*h6*
     &  Log(H2SO4)**3 + 0.00038146321602426414*q*Log(H2SO4)**3 - 
     &  0.011997306640487447*RH*Log(H2SO4)**3 + 0.0003857955558500776*T*
     &  Log(H2SO4)**3 - 0.004779827937902779*Log(H2SO4)**4

      h3=EXP(h3)


	
	h1=456229.3726785317 - 696754.0061755505/h3 - 
     &  8.954389043957226e7*h6 + (1.4677717736521986e8*h6)/h3 + 
     &  1867.5296995211318*q - (2798.172491398116*q)/h3 + 
     &  1500.05530404756*h6*q - (171625.68387665015*h6*q)/h3 - 
     &  5924.898937400813*T + (7657.054762118453*T)/h3 + 
     &  1.1074613167489376e6*h6*T - (1.556115054313954e6*h6*T)/h3 - 
     &  28.49469750360138*q*T + (50.46283217861981*q*T)/h3 + 
     &  1753.7251886642075*h6*q*T - (975.4050494746682*h6*q*T)/h3 + 
     &  25.45577410331681*T**2 - (26.190051353137182*T**2)/h3 - 
     &  4516.4884641323815*h6*T**2 + (5030.318058537955*h6*T**2)/h3 + 
     &  0.14403307295844858*q*T**2 - (0.2897062226196713*q*T**2)/h3 - 
     &  15.561717210345154*h6*q*T**2+(18.061492698402766*h6*q*T**2)/h3- 
     &  0.03621852394282133*T**3 + (0.02659567895141871*T**3)/h3 + 
     &  6.0696427040893886*h6*T**3 - (4.585112604562606*h6*T**3)/h3 - 
     &  0.00024142548506534756*q*T**3 + 
     &  (0.0005400351895947918*q*T**3)/h3+0.0341853829283136*h6*q*T**3- 
     &  (0.04459905758423761*h6*q*T**3)/h3-81941.3895777097*Log(H2SO4)+ 
     &  (90462.54571516594*Log(H2SO4))/h3 + 
     &  1.6483362741354546e7*h6*Log(H2SO4) - 
     &  (2.3542460879418086e7*h6*Log(H2SO4))/h3 - 
     &  319.8790300996385*q*Log(H2SO4) + 
     &  (226.84336756604796*q*Log(H2SO4))/h3 - 
     &  26455.51148372377*h6*q*Log(H2SO4) + 
     &  (60326.79915791394*h6*q*Log(H2SO4))/h3 + 
     &  1063.8621634534027*T*Log(H2SO4) - 
     &  (891.6122304868614*T*Log(H2SO4))/h3 - 
     &  203718.41653995324*h6*T*Log(H2SO4) + 
     &  (243001.1509313233*h6*T*Log(H2SO4))/h3 + 
     &  4.963508394566835*q*T*Log(H2SO4) - 
     &  (5.401955915120605*q*T*Log(H2SO4))/h3 + 
     &  25.811543727271804*h6*q*T*Log(H2SO4) - 
     &  (355.292262810142*h6*q*T*Log(H2SO4))/h3 - 
     &  4.5691538346459195*T**2*Log(H2SO4) + 
     &  (2.468423948121569*T**2*Log(H2SO4))/h3 + 
     &  830.195816015575*h6*T**2*Log(H2SO4) - 
     &  (745.2009340266254*h6*T**2*Log(H2SO4))/h3 - 
     &  0.025450884793884514*q*T**2*Log(H2SO4) + 
     &  (0.03593327468289806*q*T**2*Log(H2SO4))/h3 + 
     &  1.3302486952619124*h6*q*T**2*Log(H2SO4) - 
     &  (0.320851018492291*h6*q*T**2*Log(H2SO4))/h3 + 
     &  0.006498242231697865*T**3*Log(H2SO4) - 
     &  (0.0013540108127227987*T**3*Log(H2SO4))/h3 - 
     &  1.1150374268629073*h6*T**3*Log(H2SO4) + 
     &  (0.5953715882656012*h6*T**3*Log(H2SO4))/h3 + 
     &  0.000043146778270078744*q*T**3*Log(H2SO4) - 
     &  (0.0000735656618351653*q*T**3*Log(H2SO4))/h3 - 
     &  0.004059155580949069*h6*q*T**3*Log(H2SO4) + 
     &  (0.0028886014749212705*h6*q*T**3*Log(H2SO4))/h3 + 
     &  4912.2066512397205*Log(H2SO4)**2 - 
     &  (3022.0044642322964*Log(H2SO4)**2)/h3 - 
     &  1.0095170363345986e6*h6*Log(H2SO4)**2 + 
     &  (1.1998313113698033e6*h6*Log(H2SO4)**2)/h3 + 
     &  18.018115617541998*q*Log(H2SO4)**2 + 
     &  (7.561352267717425*q*Log(H2SO4)**2)/h3 + 
     &  3133.188395383639*h6*q*Log(H2SO4)**2 - 
     &  (3338.848456178725*h6*q*Log(H2SO4)**2)/h3 - 
     &  63.77612570995598*T*Log(H2SO4)**2 + 
     &  (20.100257708705428*T*Log(H2SO4)**2)/h3 + 
     &  12469.209965937996*h6*T*Log(H2SO4)**2 - 
     &  (11883.666207080118*h6*T*Log(H2SO4)**2)/h3 - 
     &  0.2853187730266578*q*T*Log(H2SO4)**2 + 
     &  (0.03880844916086884*q*T*Log(H2SO4)**2)/h3 - 
     &  21.708027380830988*h6*q*T*Log(H2SO4)**2 + 
     &  (28.111166854183463*h6*q*T*Log(H2SO4)**2)/h3 + 
     &  0.27388018029546984*T**2*Log(H2SO4)**2 + 
     &  (0.005568448436161079*T**2*Log(H2SO4)**2)/h3 - 
     &  50.78398945757657*h6*T**2*Log(H2SO4)**2 + 
     &  (33.25484929904544*h6*T**2*Log(H2SO4)**2)/h3 + 
     &  0.0014875097131454307*q*T**2*Log(H2SO4)**2 - 
     &  (0.0008810798671059988*q*T**2*Log(H2SO4)**2)/h3 + 
     &  0.007045532214236005*h6*q*T**2*Log(H2SO4)**2 - 
     &  (0.05521643331038819*h6*q*T**2*Log(H2SO4)**2)/h3 - 
     &  0.00038943745406964085*T**3*Log(H2SO4)**2 - 
     &  (0.00015308439633013156*T**3*Log(H2SO4)**2)/h3 + 
     &  0.06817761559886615*h6*T**3*Log(H2SO4)**2 - 
     &  (0.019520711865759842*h6*T**3*Log(H2SO4)**2)/h3 - 
     &  2.554240267270117e-6*q*T**3*Log(H2SO4)**2 + 
     &  (2.526957991768314e-6*q*T**3*Log(H2SO4)**2)/h3 + 
     &  0.00011978108924107446*h6*q*T**3*Log(H2SO4)**2 - 
     &  98.09290363960072*Log(H2SO4)**3 + 
     &  (3.8697214915504805*Log(H2SO4)**3)/h3 + 
     &  20559.743351544734*h6*Log(H2SO4)**3 - 
     &  (18687.867717989076*h6*Log(H2SO4)**3)/h3 - 
     &  0.33118998078526996*q*Log(H2SO4)**3 - 
     &  (0.6955616584928909*q*Log(H2SO4)**3)/h3 - 
     &  92.45613590327787*h6*q*Log(H2SO4)**3 + 
     &  (8.948488896745264*h6*q*Log(H2SO4)**3)/h3 + 
     &  1.273863980203035*T*Log(H2SO4)**3 + 
     &  (0.38051947439108685*T*Log(H2SO4)**3)/h3 - 
     &  253.8228380167286*h6*T*Log(H2SO4)**3 + 
     &  (171.19178070061787*h6*T*Log(H2SO4)**3)/h3 + 
     &  0.005376622635136611*q*T*Log(H2SO4)**3 + 
     &  (0.006658042079999949*q*T*Log(H2SO4)**3)/h3 + 
     &  0.8230892248689082*h6*q*T*Log(H2SO4)**3 - 
     &  (0.058021973134188776*h6*q*T*Log(H2SO4)**3)/h3 - 
     &  0.005471105747328322*T**2*Log(H2SO4)**3 - 
     &  (0.003699167378319841*T**2*Log(H2SO4)**3)/h3 + 
     &  1.0332570198673408*h6*T**2*Log(H2SO4)**3 - 
     &  (0.38600322062471826*h6*T**2*Log(H2SO4)**3)/h3 - 
     &  0.000028584962254144285*q*T**2*Log(H2SO4)**3 - 
     &  (0.000016062832568931284*q*T**2*Log(H2SO4)**3)/h3 - 
     &  0.001819022861961261*h6*q*T**2*Log(H2SO4)**3 + 
     &  7.779696142330516e-6*T**3*Log(H2SO4)**3 + 
     &  (8.513832613126146e-6*T**3*Log(H2SO4)**3)/h3 - 
     &  0.0013866861532058*h6*T**3*Log(H2SO4)**3 + 
     &  4.981099027173453e-8*q*T**3*Log(H2SO4)**3 + 
     &  178821.264938151*Log(RH) - (634125.5419575685*Log(RH))/h3 - 
     &  2.310391119522488e7*h6*Log(RH) + 
     &  (3.564836811852359e7*h6*Log(RH))/h3 - 
     &  744.9440395010726*q*Log(RH)+(1019.4951846825716*q*Log(RH))/h3- 
     &  28008.93484108728*h6*q*Log(RH) - 
     &  (7136.405275749431*h6*q*Log(RH))/h3 - 
     &  2275.393426900043*T*Log(RH) + (9231.8244493922*T*Log(RH))/h3 + 
     &  237571.37589421016*h6*T*Log(RH) - 
     &  (295411.0553218543*h6*T*Log(RH))/h3 + 
     &  4.451622613453289*q*T*Log(RH) - 
     &  (4.150562571935793*q*T*Log(RH))/h3 + 
     &  595.8398944235868*h6*q*T*Log(RH) - 
     &  (210.09033695512161*h6*q*T*Log(RH))/h3 + 
     &  9.678343416541734*T**2*Log(RH) - 
     &  (45.69882758602488*T**2*Log(RH))/h3 - 
     &  747.9878102876322*h6*T**2*Log(RH) + 
     &  (649.4424346377774*h6*T**2*Log(RH))/h3 + 
     &  0.005153649541788936*q*T**2*Log(RH) - 
     &  (0.0023655719118541533*q*T**2*Log(RH))/h3 - 
     &  2.848932474900831*h6*q*T**2*Log(RH) + 
     &  (0.2889680076326506*h6*q*T**2*Log(RH))/h3 - 
     &  0.013824587420637588*T**3*Log(RH) + 
     &  (0.07673955581386811*T**3*Log(RH))/h3 + 
     &  0.6813253055163031*h6*T**3*Log(RH) - 
     &  (0.2663343569482817*h6*T**3*Log(RH))/h3 - 
     &  0.00004364745106737275*q*T**3*Log(RH) - 
     &  (9.434880576275503e-6*q*T**3*Log(RH))/h3 + 
     &  0.0030430214338531205*h6*q*T**3*Log(RH) + 
     &  (0.00371157673745342*h6*q*T**3*Log(RH))/h3 - 
     &  31842.277946766168*Log(H2SO4)*Log(RH) + 
     &  (68134.57926920953*Log(H2SO4)*Log(RH))/h3 + 
     &  3.493985612508899e6*h6*Log(H2SO4)*Log(RH) - 
     &  (4.587173661378969e6*h6*Log(H2SO4)*Log(RH))/h3 + 
     &  190.19310993267845*q*Log(H2SO4)*Log(RH) - 
     &  (233.68205749312318*q*Log(H2SO4)*Log(RH))/h3 - 
     &  1547.3248018686509*h6*q*Log(H2SO4)*Log(RH) + 
     &  (4045.133810971437*h6*q*Log(H2SO4)*Log(RH))/h3 + 
     &  405.31725032582267*T*Log(H2SO4)*Log(RH) - 
     &  (1057.655396923827*T*Log(H2SO4)*Log(RH))/h3 - 
     &  33903.90722480216*h6*T*Log(H2SO4)*Log(RH) + 
     &  (32965.76708786528*h6*T*Log(H2SO4)*Log(RH))/h3 - 
     &  1.5034587825111836*q*T*Log(H2SO4)*Log(RH) + 
     &  (1.106270806420077*q*T*Log(H2SO4)*Log(RH))/h3 - 
     &  25.549725610666396*h6*q*T*Log(H2SO4)*Log(RH) + 
     &  (20.09072468984517*h6*q*T*Log(H2SO4)*Log(RH))/h3 - 
     &  1.7253372988643991*T**2*Log(H2SO4)*Log(RH) + 
     &  (5.612615323997249*T**2*Log(H2SO4)*Log(RH))/h3 + 
     &  95.62445916373963*h6*T**2*Log(H2SO4)*Log(RH) - 
     &  (48.96291394056088*h6*T**2*Log(H2SO4)*Log(RH))/h3 + 
     &  0.0019300142376745193*q*T**2*Log(H2SO4)*Log(RH) + 
     &  (0.00015584582354154648*q*T**2*Log(H2SO4)*Log(RH))/h3 + 
     &  0.19177967246162606*h6*q*T**2*Log(H2SO4)*Log(RH) - 
     &  (0.16871358892301552*h6*q*T**2*Log(H2SO4)*Log(RH))/h3 + 
     &  0.00246781658651071*T**3*Log(H2SO4)*Log(RH) - 
     &  (0.01009444649935961*T**3*Log(H2SO4)*Log(RH))/h3 - 
     &  0.06589616928988701*h6*T**3*Log(H2SO4)*Log(RH) - 
     &  (0.01596274185896012*h6*T**3*Log(H2SO4)*Log(RH))/h3 + 
     &  4.109041091272658e-6*q*T**3*Log(H2SO4)*Log(RH) + 
     &  (2.041743340493583e-7*q*T**3*Log(H2SO4)*Log(RH))/h3 - 
     &  0.0001572125066332154*h6*q*T**3*Log(H2SO4)*Log(RH) + 
     &  1932.4011188760162*Log(H2SO4)**2*Log(RH) - 
     &  (781.4228432455719*Log(H2SO4)**2*Log(RH))/h3 - 
     &  181291.39888342103*h6*Log(H2SO4)**2*Log(RH) + 
     &  (210184.24050416853*h6*Log(H2SO4)**2*Log(RH))/h3 - 
     &  13.666906354559005*q*Log(H2SO4)**2*Log(RH) + 
     &  (15.930151624229865*q*Log(H2SO4)**2*Log(RH))/h3 + 
     &  305.7703887043172*h6*q*Log(H2SO4)**2*Log(RH) - 
     &  (391.71164595591085*h6*q*Log(H2SO4)**2*Log(RH))/h3 - 
     &  24.621398404913865*T*Log(H2SO4)**2*Log(RH) + 
     &  (20.079316822912197*T*Log(H2SO4)**2*Log(RH))/h3 + 
     &  1658.5845815289838*h6*T*Log(H2SO4)**2*Log(RH) - 
     &  (1404.5530139788664*h6*T*Log(H2SO4)**2*Log(RH))/h3 + 
     &  0.11823318365930788*q*T*Log(H2SO4)**2*Log(RH) - 
     &  (0.08056635668752316*q*T*Log(H2SO4)**2*Log(RH))/h3 - 
     &  1.0034435154436134*h6*q*T*Log(H2SO4)**2*Log(RH) + 
     &  (1.6021411957310931*h6*q*T*Log(H2SO4)**2*Log(RH))/h3 + 
     &  0.10491613716181539*T**2*Log(H2SO4)**2*Log(RH) - 
     &  (0.1492097929648805*T**2*Log(H2SO4)**2*Log(RH))/h3 - 
     &  4.0934834940398215*h6*T**2*Log(H2SO4)**2*Log(RH) + 
     &  (1.896735315742051*h6*T**2*Log(H2SO4)**2*Log(RH))/h3 - 
     &  0.00022870693934890214*q*T**2*Log(H2SO4)**2*Log(RH) + 
     &  (0.000013175311622036916*q*T**2*Log(H2SO4)**2*Log(RH))/h3- 
     &  0.00231560013471937*h6*q*T**2*Log(H2SO4)**2*Log(RH) - 
     &  0.00015023186988448507*T**3*Log(H2SO4)**2*Log(RH) + 
     &  (0.0003388679191513823*T**3*Log(H2SO4)**2*Log(RH))/h3 + 
     &  0.0015784630518057667*h6*T**3*Log(H2SO4)**2*Log(RH) - 
     &  1.0056424638350849e-7*q*T**3*Log(H2SO4)**2*Log(RH) - 
     &  39.74503165968021*Log(H2SO4)**3*Log(RH) - 
     &  (68.29220045487304*Log(H2SO4)**3*Log(RH))/h3 + 
     &  3236.6794930454566*h6*Log(H2SO4)**3*Log(RH) - 
     &  (2511.4313482742737*h6*Log(H2SO4)**3*Log(RH))/h3 + 
     &  0.3044230466506939*q*Log(H2SO4)**3*Log(RH) - 
     &  (0.321694018076018*q*Log(H2SO4)**3*Log(RH))/h3 - 
     &  7.83632554091337*h6*q*Log(H2SO4)**3*Log(RH) + 
     &  (1.2615699211576334*h6*q*Log(H2SO4)**3*Log(RH))/h3 + 
     &  0.5071743553774433*T*Log(H2SO4)**3*Log(RH) + 
     &  (0.6993184748948513*T*Log(H2SO4)**3*Log(RH))/h3 - 
     &  28.064969097031586*h6*T*Log(H2SO4)**3*Log(RH) + 
     &  (12.010607450658807*h6*T*Log(H2SO4)**3*Log(RH))/h3 - 
     &  0.0027511166431233884*q*T*Log(H2SO4)**3*Log(RH) + 
     &  (0.0015662272224653852*q*T*Log(H2SO4)**3*Log(RH))/h3 + 
     &  0.03815288545975011*h6*q*T*Log(H2SO4)**3*Log(RH) - 
     &  0.0021639947693271114*T**2*Log(H2SO4)**3*Log(RH) - 
     &  (0.0017767382579465334*T**2*Log(H2SO4)**3*Log(RH))/h3 + 
     &  0.059762491440998364*h6*T**2*Log(H2SO4)**3*Log(RH) + 
     &  6.130388542051968e-6*q*T**2*Log(H2SO4)**3*Log(RH) + 
     &  3.1017462628108047e-6*T**3*Log(H2SO4)**3*Log(RH) - 
     &  7689.12786121461*Log(RH)**2-(25856.380361656236*Log(RH)**2)/h3- 
     &  1.3915437178042033e6*h6*Log(RH)**2 + 
     &  (1.9928171556602388e6*h6*Log(RH)**2)/h3 + 
     &  53.79216985973448*q*Log(RH)**2 - 
     &  (134.8787674139501*q*Log(RH)**2)/h3 - 
     &  1440.4837667461468*h6*q*Log(RH)**2 - 
     &  (3141.0130344871345*h6*q*Log(RH)**2)/h3 + 
     &  94.7545268444358*T*Log(RH)**2 + 
     &  (342.5857682014369*T*Log(RH)**2)/h3 + 
     &  11341.294070320842*h6*T*Log(RH)**2 - 
     &  (18070.160441620144*h6*T*Log(RH)**2)/h3 - 
     &  0.44696864139255954*q*T*Log(RH)**2 + 
     &  (1.1199451205710647*q*T*Log(RH)**2)/h3 + 
     &  19.01601691757823*h6*q*T*Log(RH)**2 + 
     &  (31.604766584775273*h6*q*T*Log(RH)**2)/h3 - 
     &  0.3470127737474881*T**2*Log(RH)**2 - 
     &  (1.4877879003131458*T**2*Log(RH)**2)/h3 - 
     &  18.727169614346867*h6*T**2*Log(RH)**2 + 
     &  (42.5448759900235*h6*T**2*Log(RH)**2)/h3 + 
     &  0.0007727701515171791*q*T**2*Log(RH)**2 - 
     &  (0.0032041999804810098*q*T**2*Log(RH)**2)/h3 - 
     &  0.055631814243943*h6*q*T**2*Log(RH)**2 - 
     &  (0.07696588905487615*h6*q*T**2*Log(RH)**2)/h3 + 
     &  0.00033600378567117826*T**3*Log(RH)**2 + 
     &  (0.002069337837641118*T**3*Log(RH)**2)/h3 - 
     &  0.015788853514098967*h6*T**3*Log(RH)**2 - 
     &  (0.009301939734036398*h6*T**3*Log(RH)**2)/h3 + 
     &  5.572064421478949e-7*q*T**3*Log(RH)**2 + 
     &  (5.677455521083783e-6*q*T**3*Log(RH)**2)/h3 - 
     &  0.00003507762649727548*h6*q*T**3*Log(RH)**2 + 
     &  1392.9578972108989*Log(H2SO4)*Log(RH)**2 + 
     &  (1911.734041113941*Log(H2SO4)*Log(RH)**2)/h3 + 
     &  179329.36146873524*h6*Log(H2SO4)*Log(RH)**2 - 
     &  (119328.53674704138*h6*Log(H2SO4)*Log(RH)**2)/h3 - 
     &  6.462104618811676*q*Log(H2SO4)*Log(RH)**2 + 
     &  (11.294074621589449*q*Log(H2SO4)*Log(RH)**2)/h3 + 
     &  13.39033080899747*h6*q*Log(H2SO4)*Log(RH)**2 - 
     &  (8.52831107206456*h6*q*Log(H2SO4)*Log(RH)**2)/h3 - 
     &  16.8958062365253*T*Log(H2SO4)*Log(RH)**2 - 
     &  (24.81735534395592*T*Log(H2SO4)*Log(RH)**2)/h3 - 
     &  1464.394638508739*h6*T*Log(H2SO4)*Log(RH)**2 + 
     &  (981.4137808275404*h6*T*Log(H2SO4)*Log(RH)**2)/h3 + 
     &  0.051952833287610724*q*T*Log(H2SO4)*Log(RH)**2 - 
     &  (0.061982178004231*q*T*Log(H2SO4)*Log(RH)**2)/h3 - 
     &  0.7507183329339466*h6*q*T*Log(H2SO4)*Log(RH)**2 - 
     &  (0.05082083499340442*h6*q*T*Log(H2SO4)*Log(RH)**2)/h3 + 
     &  0.060317519684020285*T**2*Log(H2SO4)*Log(RH)**2 + 
     &  (0.10238074986938017*T**2*Log(H2SO4)*Log(RH)**2)/h3 + 
     &  2.654746963853928*h6*T**2*Log(H2SO4)*Log(RH)**2 - 
     &  (1.9505367684260284*h6*T**2*Log(H2SO4)*Log(RH)**2)/h3 - 
     &  0.00008828022840243558*q*T**2*Log(H2SO4)*Log(RH)**2 - 
     &  (7.843133775464917e-6*q*T**2*Log(H2SO4)*Log(RH)**2)/h3 + 
     &  0.00439346543661198*h6*q*T**2*Log(H2SO4)*Log(RH)**2 - 
     &  0.00005488176049657228*T**3*Log(H2SO4)*Log(RH)**2 - 
     &  (0.00012357549613231235*T**3*Log(H2SO4)*Log(RH)**2)/h3 + 
     &  0.000877040667994904*h6*T**3*Log(H2SO4)*Log(RH)**2 - 
     &  4.394231627146249e-8*q*T**3*Log(H2SO4)*Log(RH)**2 - 
     &  72.11406227014747*Log(H2SO4)**2*Log(RH)**2 - 
     &  (16.86445272085156*Log(H2SO4)**2*Log(RH)**2)/h3 - 
     &  6535.765336704547*h6*Log(H2SO4)**2*Log(RH)**2 + 
     &  (941.2813201741736*h6*Log(H2SO4)**2*Log(RH)**2)/h3 + 
     &  0.2428031989628447*q*Log(H2SO4)**2*Log(RH)**2 - 
     &  (0.3050220578412721*q*Log(H2SO4)**2*Log(RH)**2)/h3 + 
     &  3.3649457326647947*h6*q*Log(H2SO4)**2*Log(RH)**2 + 
     &  (0.9000970248172522*h6*q*Log(H2SO4)**2*Log(RH)**2)/h3 + 
     &  0.8503624487797692*T*Log(H2SO4)**2*Log(RH)**2 + 
     &  (0.2527834506442179*T*Log(H2SO4)**2*Log(RH)**2)/h3 + 
     &  50.73348890297342*h6*T*Log(H2SO4)**2*Log(RH)**2 - 
     &  (4.356161470448808*h6*T*Log(H2SO4)**2*Log(RH)**2)/h3 - 
     &  0.0018342607279861292*q*T*Log(H2SO4)**2*Log(RH)**2 + 
     &  (0.0020365272689109475*q*T*Log(H2SO4)**2*Log(RH)**2)/h3 - 
     &  0.03276629834613374*h6*q*T*Log(H2SO4)**2*Log(RH)**2 - 
     &  0.0028420227040081517*T**2*Log(H2SO4)**2*Log(RH)**2 - 
     &  (0.0009263932662681723*T**2*Log(H2SO4)**2*Log(RH)**2)/h3 - 
     &  0.08964774960777555*h6*T**2*Log(H2SO4)**2*Log(RH)**2 + 
     &  3.036637047548649e-6*q*T**2*Log(H2SO4)**2*Log(RH)**2 + 
     &  2.0594310710626316e-6*T**3*Log(H2SO4)**2*Log(RH)**2 + 
     &  0.9950338851031777*Log(H2SO4)**3*Log(RH)**2 - 
     &  (0.622965174549953*Log(H2SO4)**3*Log(RH)**2)/h3 + 
     &  48.41845482837981*h6*Log(H2SO4)**3*Log(RH)**2 - 
     &  (1.7738934087048919*h6*Log(H2SO4)**3*Log(RH)**2)/h3 - 
     &  0.0022234755298572505*q*Log(H2SO4)**3*Log(RH)**2 - 
     &  (0.002437957119890622*q*Log(H2SO4)**3*Log(RH)**2)/h3 + 
     &  0.08069009533680993*h6*q*Log(H2SO4)**3*Log(RH)**2 - 
     &  0.010935011617216193*T*Log(H2SO4)**3*Log(RH)**2 + 
     &  (0.003456846985810427*T*Log(H2SO4)**3*Log(RH)**2)/h3 - 
     &  0.24482068250824904*h6*T*Log(H2SO4)**3*Log(RH)**2 + 
     &  0.000011334503487127534*q*T*Log(H2SO4)**3*Log(RH)**2 + 
     &  0.000029425270779265584*T**2*Log(H2SO4)**3*Log(RH)**2

	  h1=EXP(h1)


	h2=-32043.03148295406 + 59725.428570008815/h3 +  
     & 7.128537634261564e6*h6 - (1.3833467233343722e7*h6)/h3+
     & 33.63110252227136*q-(48.61215633992165*q)/h3 - 16602.414377611287
     & *h6*q + (40754.788181739124* h6*q)/h3 -2.3397851800516185*q**2 + 
     & (4.426964073992281*q**2)/h3 + 18.971418767591036*h6*q**2 +
     & (60.446718551038344*h6*q**2)/h3+396.33752593131607*T-
     & (650.8601684277011*T)/h3-83753.25337253512*h6*T+
     & (140932.8771905448*h6*T)/h3 + 1.7779514612590905*q*T -  
     & (4.1201474547289845*q*T)/h3 - 14.23100399324848*h6*q*T - 
     & (199.64146093004214*h6*q*T)/h3 + 0.017586270137944494*q**2*T - 
     & (0.040669111407675304*q**2*T)/h3 - 0.5190592639152767*h6*q**2*T + 
     & (0.4152291779457336*h6*q**2*T)/h3 - 1.6239261240452716*T**2 + 
     & (2.290370025362425*T**2)/h3 + 323.71145981054207*h6*T**2 -  
     & (443.70978129586814*h6*T**2)/h3 - 0.016944490878752067*q*T**2 + 
     & (0.03190748579500371*q*T**2)/h3 + 0.9921294086888162*h6*q*T**2 -  
     & (0.017807643647025154*h6*q*T**2)/h3-0.000026182356192850867*q**2*
     & T**2+(0.00018413948144606674*q**2*T**2)/h3+0.001418259144142961*
     & h6*q**2*T**2 - (0.0005314099952660705*h6*q**2*T**2)/h3 + 
     & 0.002203350630670886*T**3-(0.002628197592847762*T**3)/h3 - 
     & 0.41364399538414304*h6*T**3+(0.4236315079925672*h6*T**3)/h3 + 
     & 0.000036034470518775974*q*T**3-(0.00005072877953687334*q*T**3)/h3
     & -0.0025779163615301994*h6*q*T**3-(0.0001801635475946943*h6*q*
     & T**3)/h3 - 1.0397496813871261e-8*q**2*T**3 -  
     & (4.794510915054545e-7*q**2*T**3)/h3 + 1.0037940606886203e-6*h6*
     & q**2*T**3 +4057.4863055406972*Log(H2SO4) - (4857.4895775507075*
     & Log(H2SO4))/h3 - 873923.0975121224*h6*Log(H2SO4) + 
     & (1.3884698205066123e6*h6*Log(H2SO4))/h3-  
     & 16.567334471680468*q*Log(H2SO4)+(39.06695157867744*q*
     & Log(H2SO4))/h3 + 
     & 4033.3161579197213*h6*q*Log(H2SO4) - (4837.3892389479715*h6*q* 
     & Log(H2SO4))/h3 + 0.32643050327864925*q**2*Log(H2SO4) -  
     & (0.29168894526478634*q**2*Log(H2SO4))/h3 +   3.0197469938934343* 
     & h6*q**2*Log(H2SO4) - (13.620934035683936*h6*q**2*Log(H2SO4))/h3 -  
     & 50.64853826213881*T*Log(H2SO4) + (46.190674056264115*T*
     & Log(H2SO4))/h3 +10309.186386370407*h6*T*Log(H2SO4) - 
     & (13393.112078195805*h6*T* 
     & Log(H2SO4))/h3 - 0.05988575301198189*q*T*Log(H2SO4) + 
     & (0.12991016720908877*q*T*Log(H2SO4))/h3-25.211481437032155*h6*q*
     & T*Log(H2SO4) + (25.76780209283505*h6*q*T*Log(H2SO4))/h3 -  
     & 0.0025433151017881413*q**2*T*Log(H2SO4)-(0.00003896417606619728* 
     & q**2*T*Log(H2SO4))/h3+0.018638367565640454*h6*q**2*T*Log(H2SO4)-  
     & (0.009466280650899394*h6*q**2*T*Log(H2SO4))/h3+
     & 0.20947015213812362*T**2*Log(H2SO4)-(0.1282261301435015*T**2*
     & Log(H2SO4))/h3-40.05118528139879*h6*T**2*Log(H2SO4)+
     & (37.82992227601891* h6*T**2*Log(H2SO4))/h3+0.001426085690256897*
     & q*T**2*Log(H2SO4)-(0.0026110978166799066*q*T**2*Log(H2SO4))/h3 -
     &  0.0017889111275040672*h6*q*T**2*Log(H2SO4) - 
     & (0.0010479732179162718* h6*q*T**2*Log(H2SO4))/h3 + 
     & 4.179579131654065e-6 * q**2*T**2*Log(H2SO4) +
     & (7.970218003902365e-6 * q**2 * 
     & T**2*Log(H2SO4))/h3-0.000116249941623789*h6*q**2*T**2*Log(H2SO4)-  
     & 0.0002870253324673844*T**3*Log(H2SO4)+(0.00009186981667688776*
     & T**3*Log(H2SO4))/h3 + 0.05147918448072947*h6*T**3*Log(H2SO4) -  
     & (0.028081229965757286*h6*T**3*Log(H2SO4))/h3 - 
     & 3.560981861287522e-6* q* 
     & T**3*Log(H2SO4) + (4.961243888005889e-6*q*T**3*Log(H2SO4))/h3 +  
     & 0.0001420993695424454*h6*q*T**3*Log(H2SO4)+1.2719812342415112e-9*
     & q**2*T**3*Log(H2SO4) - 124.19545742393447*Log(H2SO4)**2 +   
     & (41.54331926307027*Log(H2SO4)**2)/h3 + 26490.80092234763*h6*
     & Log(H2SO4)**2- (31462.37707525288*h6*Log(H2SO4)**2)/h3 +  
     & 0.8500524958176608*q* Log(H2SO4)**2 - (3.2674516205233317*q*
     & Log(H2SO4)**2)/h3 - 170.8115272895333*h6*q*Log(H2SO4)**2 + 
     & (120.16514929908757*h6*q* Log(H2SO4)**2)/h3 - 
     & 0.011268864148405246*q**2*Log(H2SO4)**2 +(0.019897271615564007*
     & q**2*Log(H2SO4)**2)/h3 -  0.2574159414899909* 
     & h6*q**2*Log(H2SO4)**2 + (0.5058681955079796*h6*q**2*
     & Log(H2SO4)**2)/h3 +1.5616031047481935*T*Log(H2SO4)**2 +  
     & (0.04762580087464638*T* Log(H2SO4)**2)/h3 - 313.57384684195245*
     & h6*T*Log(H2SO4)**2+(270.50302801880633*h6*T*Log(H2SO4)**2)/h3 -
     &  0.0026826956178692776* q* 
     & T*Log(H2SO4)**2 + (0.02263655181544639*q*T*Log(H2SO4)**2)/h3 +   
     & 1.4230688700541085*h6*q*T*Log(H2SO4)**2 - (0.6649649390384627*
     &  h6* q* T*Log(H2SO4)**2)/h3 + 0.00009215121101630449*q**2*T*
     & Log(H2SO4)**2 -  
     & (0.00010215400234653555*q**2*T*Log(H2SO4)**2)/h3 +  
     & 0.0010136511879910042*h6*q**2*T*Log(H2SO4)**2 -  
     & 0.006504837589925847*T**2*Log(H2SO4)**2 -  (0.002702979699750996* 
     & T**2*Log(H2SO4)**2)/h3 + 1.223462360247441*h6*T**2*Log(H2SO4)**2- 
     & (0.5623561810520002*h6*T**2*Log(H2SO4)**2)/h3 - 
     & 0.000024272415265506407*q*T**2*Log(H2SO4)**2 - 
     & (0.0000317372505770328*q*T**2*Log(H2SO4)**2)/h3 -  
     & 0.0028928533476792*h6*q*T**2*Log(H2SO4)**2-1.7589107450560238e-7* 
     & q**2*T**2*Log(H2SO4)**2 + 8.98009092200891e-6*T**3*Log(H2SO4)**2+   
     & (7.173660190761036e-6*T**3*Log(H2SO4)**2)/h3 -  
     & 0.0015801362599946241*h6*T**3*Log(H2SO4)**2+8.220297383016583e-8*
     & q*T**3*Log(H2SO4)**2 - 12589.220049398413*Log(RH)+
     & (71533.62210328173*Log(RH))/h3 + 1.0678104434003264e6*h6*Log(RH)-
     & (1.6102655953002474e6*h6*Log(RH))/h3+51.82639715490156*q*Log(RH)-
     & (38.71320196661908*q*Log(RH))/h3+3767.8963041701336*h6*q*Log(RH)-
     & (2025.2741146328435*h6*q*Log(RH))/h3+0.43253945689885376*q**2*
     & Log(RH) - (0.5285905724764535*q**2*Log(RH))/h3+2.800088035800146*
     & h6*q**2*Log(RH)+(5.830987444741845*h6*q**2*Log(RH))/h3+
     & 156.76303859222733*T*Log(RH) -  (1050.9788795963977*T* 
     & Log(RH))/h3 + 2968.558570914762*h6*T*Log(RH)-  
     & (13144.694726104603*h6*T*Log(RH))/h3-0.3939836553978611*q*T*
     & Log(RH) - (0.15149542273585964*q*T*Log(RH))/h3-52.00514647221519*
     & h6*q*T*Log(RH)+(33.50499777824741*h6*q*T*Log(RH))/h3-
     & 0.0044939367750510195*q**2*T*Log(RH)+  
     & (0.0037387212699065446*q**2*T*Log(RH))/h3 -  
     & 0.019659867844446746*h6*q**2*T*Log(RH) - (0.013728505197733986*
     & h6*q**2*T*Log(RH))/h3-0.6541703069347068*T**2*Log(RH) + 
     & (5.190500347012039*T**2*Log(RH))/h3-  
     & 87.05017506628032*h6*T**2*Log(RH) +  (179.83498689899525*h6*T**2* 
     & Log(RH))/h3 + 0.00042182566738198824*q*T**2*Log(RH) +  
     & (0.004271845253967774*q*T**2*Log(RH))/h3+0.12589872273911631*h6*
     & q* T**2*Log(RH) - (0.08117178992783239*h6*q*T**2*Log(RH))/h3 +   
     & 0.000010942253625927065*q**2*T**2*Log(RH) - 
     & (3.4136292496834907e-6*q**2* 
     & T**2*Log(RH))/h3 + 3.7003812518449353e-6*h6*q**2*T**2*Log(RH) +  
     & 0.0009227426647911644*T**3*Log(RH) - (0.008600685587907114*T**3* 
     & Log(RH))/h3 + 0.23093344962838908*h6*T**3*Log(RH) -  
     & (0.38790097045592176*h6*T**3*Log(RH))/h3 + 1.2941739001871342e-6*
     & q*T**3*Log(RH)-(0.000013459445557723903*q*T**3*Log(RH))/h3 +  
     & 0.00016881487161048804* h6*q*T**3*Log(RH) + 2.609280671396356e-9*
     & q**2*T**3 *Log(RH)+2130.4592878663534*Log(H2SO4)*Log(RH) -  
     & (4313.3562495734495*Log(H2SO4)*Log(RH))/h3 - 
     & 291848.82747478905*h6*Log(H2SO4)*Log(RH) +  
     & (488820.18312982953*h6*Log(H2SO4)*Log(RH))/h3-9.405300971012823*
     & q*Log(H2SO4)*Log(RH) + (12.85263348658218*q*Log(H2SO4)*
     * Log(RH))/h3+102.49140755096978*h6*q*Log(H2SO4)*Log(RH) - 
     & (289.1347848442807*h6*q*Log(H2SO4)*Log(RH))/h3 - 
     & 0.01500268773941549*q**2*Log(H2SO4)*Log(RH) +  
     & (0.018600134825856034*q**2*Log(H2SO4)*Log(RH))/h3 - 
     & 0.07717355560425733* 
     & h6*q**2*Log(H2SO4)*Log(RH)-(0.161938896012148*h6*q**2*Log(H2SO4)* 
     & Log(RH))/h3 - 27.134659704074295*T*Log(H2SO4)*Log(RH) + 
     & (69.2483692103918*T*Log(H2SO4)*Log(RH))/h3 + 1882.179099085481*
     & h6*T*Log(H2SO4)*Log(RH)-(3332.724995554323*h6*T*Log(H2SO4)*
     & Log(RH))/h3 +   0.0786126279415837*q* 
     & T*Log(H2SO4)*Log(RH) - (0.12353477014379059*q*T*Log(H2SO4)*
     & Log(RH))/h3+2.221848005253306*h6*q*T*Log(H2SO4)*Log(RH) +  
     & (0.532770625490143*h6 *q*T*Log(H2SO4)*Log(RH))/h3 + 
     & 0.00020323051093599427*q**2*T* Log(H2SO4)* 
     & Log(RH) - (0.00014374306638210478*q**2*T*Log(H2SO4)*Log(RH))/h3 +   
     & 0.0006995787862129086*h6*q**2*T*Log(H2SO4)*Log(RH) +   
     & 0.11535250161710112*T**2*Log(H2SO4)*Log(RH)-(0.37254120769567123* 
     & T**2*Log(H2SO4)*Log(RH))/h3 +   0.3210533828417492*h6*T**2* 
     & Log(H2SO4)*Log(RH)+(4.899886058456062*h6*T**2*Log(H2SO4)*
     & Log(RH))/h3-0.00013223802083620908*q*T**2*Log(H2SO4)*Log(RH) +  
     & (0.00031364246871883094*q*T**2*Log(H2SO4)*Log(RH))/h3 - 
     & 0.013017818218354784*h6*q*T**2*Log(H2SO4)*Log(RH) - 
     & 6.285230956812909e-7* 
     & q**2*T**2*Log(H2SO4)*Log(RH) -0.0001644893214543928*T**3* 
     & Log(H2SO4)*Log(RH)+(0.0006683187277334718*T**3*Log(H2SO4)*
     & Log(RH))/h3-0.012564516901538326*h6*T**3*Log(H2SO4)*Log(RH)-
     & 1.2575655284907035e-7*q*T**3*Log(H2SO4)*Log(RH) - 
     & 79.49703918541012*Log(H2SO4)**2*Log(RH) -  
     & (101.22507090731304*Log(H2SO4)**2*Log(RH))/h3 + 
     & 12719.306036991808*h6*Log(H2SO4)**2*Log(RH)-(7977.5339222745315*
     & h6*Log(H2SO4)**2*Log(RH))/h3 
     & +0.3744372178099974*q*Log(H2SO4)**2*Log(RH)-(0.025682479965146*q* 
     & Log(H2SO4)**2*Log(RH))/h3-22.19218619809272*h6*q*Log(H2SO4)**2* 
     & Log(RH) + (5.910834971937343*h6*q*Log(H2SO4)**2*Log(RH))/h3 -   
     & 0.00037981877987097496*q**2*Log(H2SO4)**2*Log(RH) +  
     & (0.00032985668220561373*q**2*Log(H2SO4)**2*Log(RH))/h3 - 
     & 0.0003940881289017262*h6*q**2*Log(H2SO4)**2*Log(RH) + 
     & 1.0255615837685275*T*Log(H2SO4)**2*Log(RH) +  
     & (1.0124457331884233*T* 
     & Log(H2SO4)**2*Log(RH))/h3-111.19643387882614*h6*T*Log(H2SO4)**2* 
     & Log(RH) + (37.058545749337576*h6*T*Log(H2SO4)**2*Log(RH))/h3 - 
     & 0.0034094873321621174*q*T*Log(H2SO4)**2*Log(RH) +  
     & (0.000015597253696421123*q*T*Log(H2SO4)**2*Log(RH))/h3 +   
     & 0.10750304600542902*h6*q*T*Log(H2SO4)**2*Log(RH) + 
     &1.8808923253076581e-6* 
     & q**2*T*Log(H2SO4)**2*Log(RH) - 0.0044066707400722644*T**2*
     & Log(H2SO4)**2*Log(RH) -  (0.0025163539282379434*T**2*
     & Log(H2SO4)**2*Log(RH))/h3 + 0.24020207161395796*h6*T**2*
     & Log(H2SO4)**2*Log(RH)+7.634986243805219e-6*q*T**2*Log(H2SO4)**2*
     & Log(RH) + 6.328401246898811e-6*T**3*Log(H2SO4)**2*Log(RH) + 
     & 3630.6862033225625*Log(RH)**2-(1075.4966438125716*Log(RH)**2)/h3+
     & 54546.69751557024*h6*Log(RH)**2-(49710.530231480734*h6*
     & Log(RH)**2)/h3 - 1.5893360668636096*q*Log(RH)**2 +   
     & (7.970183065343727*q*Log(RH)**2)/h3 + 400.88335315886525*h6*q*
     & Log(RH)**2 -(219.0329650802132*h6*q*Log(RH)**2)/h3 +  
     & 0.015009992625400783*q**2*Log(RH)**2 + (0.0036804085191588995*
     & q**2*Log(RH)**2)/h3+0.4203732847045355*h6*q**2*Log(RH)**2 +  
     & (0.20391771437862574*h6*q**2*Log(RH)**2)/h3 - 44.579109086853684*
     & T*Log(RH)**2 + (8.53054910480424*T* 
     & Log(RH)**2)/h3 + 76.04921125760566*h6*T*Log(RH)**2 -  
     & (329.28669019015933*h6*T*Log(RH)**2)/h3 - 0.019940534315381526*q
     & *T*Log(RH)**2 - (0.050640623913744444*q*T*Log(RH)**2)/h3 - 
     & 1.5765203107981385*h6*q*T*Log(RH)**2 +   (0.3694218342124642*h6*
     & q*T* Log(RH)**2)/h3-0.0001222138921423729*q**2*T*Log(RH)**2 - 
     & (0.000049403831741944986*q**2* 
     & T*Log(RH)**2)/h3 - 0.0025870767829077484*h6*q**2*T*Log(RH)**2 +  
     & 0.17463801499094986*T**2*Log(RH)**2 -  (0.02046987546015939*T**2* 
     & Log(RH)**2)/h3 - 1.3780954635968534*h6*T**2*Log(RH)**2 +  
     & (2.0745682824062985*h6*T**2*Log(RH)**2)/h3 +  
     & 0.00015350888582331981* q*T**2*Log(RH)**2 + 
     & (0.00011285744022151752*q*T**2*Log(RH)**2)/h3 +  
     & 0.00004952110557031714*h6*q*T**2*Log(RH)**2+3.968468979853371e-7*
     & q**2*T**2*Log(RH)**2 - 0.00021374165659919817*T**3*Log(RH)**2 +  
     & (0.000023714077438108623*T**3*Log(RH)**2)/h3 -  
     & 0.0010671420677054659* h6* 
     & T**3*Log(RH)**2 - 1.8562790214558605e-7*q*T**3*Log(RH)**2 - 
     & 279.60556688841575*Log(H2SO4)*Log(RH)**2 +  (44.483779255378025* 
     & Log(H2SO4)*Log(RH)**2)/h3 -11085.947063809159*h6*Log(H2SO4)*
     & Log(RH)**2 +(10910.046235786824*h6*Log(H2SO4)*Log(RH)**2)/h3 + 
     & 0.47924678115681674*q*Log(H2SO4)*Log(RH)**2-(0.3759508982199312*
     & q*Log(H2SO4)*Log(RH)**2)/h3-29.071647637053548*h6*q*Log(H2SO4)*
     & Log(RH)**2 +  (8.887445070958586*h6*q*Log(H2SO4)*Log(RH)**2)/h3-
     & 0.00029636660559253115*q**2 * Log(H2SO4)*Log(RH)**2 + 
     & (0.00046301434118339234*q**2*Log(H2SO4)*Log(RH)**2)/h3 +  
     & 0.005479820108166143*h6*q**2*Log(H2SO4)*Log(RH)**2 +  
     & 3.2937887074139884*T*Log(H2SO4)*Log(RH)**2 -  
     & (0.053358185813312024*T*Log(H2SO4)*Log(RH)**2)/h3 + 
     & 29.117876039650476*h6*T*Log(H2SO4)*Log(RH)**2 -  
     & (34.13451610985067*h6*T*Log(H2SO4)*Log(RH)**2)/h3 -  
     & 0.001150661716592889*q*T*Log(H2SO4)*Log(RH)**2 +  
     & (0.0005497364286653244*q*T* Log(H2SO4)* 
     & Log(RH)**2)/h3+0.09171149584822953*h6*q*T*Log(H2SO4)*Log(RH)**2-  
     & 2.9335804996871176e-6*q**2*T*Log(H2SO4)*Log(RH)**2 -  
     & 0.012024558565450824*T**2*Log(H2SO4)*Log(RH)**2 - 
     & (0.0008282043267863183*T**2*Log(H2SO4)*Log(RH)**2)/h3 + 
     & 0.11350434127907331*h6* T**2* 
     &  Log(H2SO4)*Log(RH)**2 - 3.206307709389464e-6*q* T**2*Log(H2SO4)* 
     & Log(RH)**2 + 0.000012855779740625904*T**3*Log(H2SO4)*Log(RH)**2 +  
     & 3.7577104325072703*Log(H2SO4)**2*Log(RH)**2 -(2.6641547962978636* 
     & Log(H2SO4)**2*Log(RH)**2)/h3 +  476.38347866428205*h6* 
     & Log(H2SO4)**2 *Log(RH)**2 - (127.7453688897455*h6*Log(H2SO4)**2*
     & Log(RH)**2)/h3 -0.01970020856731193*q*Log(H2SO4)**2*Log(RH)**2 + 
     & (0.006983199758249699*q*Log(H2SO4)**2*Log(RH)**2)/h3 +
     &  0.3243000572140283*h6*q*Log(H2SO4)**2* 
     & Log(RH)**2+0.00002903250111944733*q**2*Log(H2SO4)**2*Log(RH)**2-  
     & 0.037500815576593044*T*Log(H2SO4)**2*Log(RH)**2 + 
     & (0.013073560271857765*T*Log(H2SO4)**2*Log(RH)**2)/h3 - 
     & 2.2753830913797373*h6*T* Log(H2SO4)**2* 
     & Log(RH)**2 + 0.00007918460467376976*q*T*Log(H2SO4)**2*Log(RH)**2+  
     & 0.00009291148493939081*T**2*Log(H2SO4)**2*Log(RH)**2

	 h2=exp(h2)



	h4=-233.3693139924163 + 3711.127600293859*h6 - 
     &  127375.45943800849*h6**2 - 0.6541599370168311*q - 
     &  8.950348936875036*h6*q + 1420.4060399615116*h6**2*q + 
     &  0.006010885721884837*q**2 - 0.2514391282801529*h6*q**2 + 
     &  11.74107168004114*h6**2*q**2 - 27.242866772851034*RH + 
     &  3230.6550683739456*h6*RH - 7739.349030052802*h6**2*RH + 
     &  0.02657310586465451*q*RH + 0.8072083676135904*h6*q*RH + 
     &  21.19451916114249*h6**2*q*RH - 0.0013789987709190107*q**2*RH - 
     &  0.02985872690339605*h6*q**2*RH + 10.1858919054768*RH**2 + 
     &  1831.5638235525591*h6*RH**2 - 3345.6757256829833*h6**2*RH**2 - 
     &  0.006667965604100408*q*RH**2 - 1.6996949068091805*h6*q*RH**2 - 
     &  0.0002938711589315329*q**2*RH**2 + 20.8476664473087*RH**3 + 
     &  153.46841353011587*h6*RH**3 - 0.023573221099724897*q*RH**3 + 
     &  3.1234153503175706*T - 106.11845318552218*h6*T + 
     &  1616.0925045006472*h6**2*T - 0.000870039433136904*q*T + 
     &  0.07485289069770486*h6*q*T - 2.8492194924162093*h6**2*q*T - 
     &  0.000025041587132070856*q**2*T -
     &  0.00009872033200025474*h6*q**2*T + 0.1139999815738393*RH*T - 
     &  59.3235586131405*h6*RH*T + 188.66428911211182*h6**2*RH*T + 
     &  0.0021232663933975797*q*RH*T + 0.03525906170060478*h6*q*RH*T + 
     &  5.50009952133025e-6*q**2*RH*T - 0.40949236112495135*RH**2*T - 
     &  9.19061496536003*h6*RH**2*T + 
     &  0.00018605332374157369*q*RH**2*T - 
     &  0.15228517482883766*RH**3*T - 0.026779347826859368*T**2 + 
     &  0.7791278560987452*h6*T**2 - 5.635512845808971*h6**2*T**2 - 
     &  0.00004290197829367679*q*T**2 - 
     &  0.00034888935617559023*h6*q*T**2 + 
     &  3.399576950217114e-8*q**2*T**2 + 
     &  0.0047495543899101315*RH*T**2 + 
     &  0.21966980586301704*h6*RH*T**2 - 
     &  5.402900826149558e-6*q*RH*T**2 + 
     &  0.00209481671640675*RH**2*T**2 + 0.00004712913787763615*T**3- 
     &  0.0022395082766246705*h6*T**3 + 7.857921629051029e-8*q*T**3 - 
     &  0.000017637397545586075*RH*T**3 + 
     &  14.487613189612244*Log(H2SO4) + 
     &  834.2954132542228*h6*Log(H2SO4) - 
     &  9818.723117205202*h6**2*Log(H2SO4) + 
     &  0.1293876064755356*q*Log(H2SO4) + 
     &  0.5266286653232314*h6*q*Log(H2SO4) - 
     &  85.4813173795405*h6**2*q*Log(H2SO4) - 
     &  0.0005441269120155648*q**2*Log(H2SO4) + 
     &  0.013203398533271429*h6*q**2*Log(H2SO4) + 
     &  2.9312753776692975*RH*Log(H2SO4) + 
     &  329.1801633051309*h6*RH*Log(H2SO4) - 
     &  2384.7570613779067*h6**2*RH*Log(H2SO4) - 
     &  0.02754902616213419*q*RH*Log(H2SO4) - 
     &  0.39067534860398334*h6*q*RH*Log(H2SO4) + 
     &  0.000059670451934966084*q**2*RH*Log(H2SO4) + 
     &  2.1796062691153963*RH**2*Log(H2SO4) - 
     &  4.33856569910847*h6*RH**2*Log(H2SO4) + 
     &  0.0016947748838741473*q*RH**2*Log(H2SO4) + 
     &  0.7647868504545748*RH**3*Log(H2SO4) + 
     &  0.22920773390043764*T*Log(H2SO4) - 
     &  8.113046295543565*h6*T*Log(H2SO4) + 
     &  68.80414049795816*h6**2*T*Log(H2SO4) + 
     &  0.0014123733231002104*q*T*Log(H2SO4) + 
     &  0.0062552686253774335*h6*q*T*Log(H2SO4) + 
     &  3.138422092020368e-8*q**2*T*Log(H2SO4) - 
     &  0.11854147612732342*RH*T*Log(H2SO4) - 
     &  2.0186043272890317*h6*RH*T*Log(H2SO4) - 
     &  0.000017510448352402354*q*RH*T*Log(H2SO4) - 
     &  0.014662481263341446*RH**2*T*Log(H2SO4) + 
     &  0.001037948515905241*T**2*Log(H2SO4) + 
     &  0.04470222473271311*h6*T**2*Log(H2SO4) - 
     &  7.507637022804968e-7*q*T**2*Log(H2SO4) + 
     &  0.0002763185325508037*RH*T**2*Log(H2SO4) - 
     &  2.046693588795596e-6*T**3*Log(H2SO4) - 
     &  3.3271552306181786*Log(H2SO4)**2 - 
     &  0.6732293679771482*h6*Log(H2SO4)**2 - 
     &  123.7357115914075*h6**2*Log(H2SO4)**2 - 
     &  0.017682151597098173*q*Log(H2SO4)**2 - 
     &  0.06932927383915821*h6*q*Log(H2SO4)**2 + 
     &  0.000028013174553488913*q**2*Log(H2SO4)**2 + 
     &  0.5385357123752617*RH*Log(H2SO4)**2 + 
     &  4.405735442137794*h6*RH*Log(H2SO4)**2 + 
     &  0.0007699683891127761*q*RH*Log(H2SO4)**2 - 
     &  0.01981622019488726*RH**2*Log(H2SO4)**2 - 
     &  0.02863075409562691*T*Log(H2SO4)**2 - 
     &  0.3815368435566718*h6*T*Log(H2SO4)**2 - 
     &  0.000032173458991045656*q*T*Log(H2SO4)**2 + 
     &  0.0006610846557652917*RH*T*Log(H2SO4)**2 + 
     &  5.708155280776268e-6*T**2*Log(H2SO4)**2 + 
     &  0.2841922065650612*Log(H2SO4)**3 + 
     &  1.937916853785423*h6*Log(H2SO4)**3 + 
     &  0.00046645968682534585*q*Log(H2SO4)**3 - 
     &  0.01392240494812053*RH*Log(H2SO4)**3 + 
     &  0.0005604898286672238*T*Log(H2SO4)**3 - 
     &  0.00648946009121241*Log(H2SO4)**4

	  H4=EXP(H4)
	

	h5=68.64045827314231-3277.3575769882523*h6 + 1.0798559249565618*q- 
     &  25.296110707348316*h6*q + 13.398992645698215*RH + 
     &  922.4932305036297*h6*RH - 0.27140107873619296*q*RH + 
     &  20.08312325165439*h6*q*RH + 66.82077511984484*RH**2 + 
     &  1611.1977384351555*h6*RH**2 - 0.02661518788217287*q*RH**2 + 
     &  3.0843227537972138*h6*q*RH**2 - 1.4080258142983926*T - 
     &  1.8568570408634648*h6*T - 0.0037450866352058397*q*T - 
     &  0.2576980690602505*h6*q*T - 1.3810906781490837*RH*T - 
     &  17.730890154257356*h6*RH*T + 0.001076745543266636*q*RH*T - 
     &  0.08909158555166723*h6*q*RH*T - 0.87735596568044*RH**2*T - 
     &  5.603639080904061*h6*RH**2*T + 0.0002224812194282388*q*RH**2*T + 
     &  0.015404311250124585*T**2 + 0.009972776386592903*h6*T**2 - 
     &  7.895537111616416e-7*q*T**2 - 0.0003918985727565079*h6*q*T**2 + 
     &  0.014090613376170114*RH*T**2 + 0.029815313352513736*h6*RH*T**2 - 
     &  3.607688962369381e-6*q*RH*T**2+0.0028248811086430334*RH**2*T**2- 
     &  0.00010105566309790024*T**3 - 0.00008368737713099654*h6*T**3 + 
     &  1.0956121068002908e-8*q*T**3 - 0.000042201374093774165*RH*T**3+ 
     &  2.429718549962505e-7*T**4 - 0.9150072698367441*Log(H2SO4) + 
     &  653.1422380210507*h6*Log(H2SO4) - 0.16169777329907886*q*
     &  Log(H2SO4)+6.860485142132651*h6*q*Log(H2SO4)+13.912279451959678*
     &  RH*Log(H2SO4) - 12.99843208521885*h6*RH*Log(H2SO4) + 
     &  0.021566098556935497*q*RH*Log(H2SO4) - 0.17455209117051437*h6*q*
     &  RH*Log(H2SO4) + 3.4003096116229066*RH**2*Log(H2SO4) - 
     &  36.48723431979357*h6*RH**2*Log(H2SO4) - 0.001787998461490877*q*
     &  RH**2*Log(H2SO4) - 0.11123959974265946*T*Log(H2SO4) +
     &  1.4641341095178422*h6*T*Log(H2SO4) + 0.0005023439404053182*q*T*
     &  Log(H2SO4) + 0.031171816066622334*h6*q*T*Log(H2SO4) - 
     &  0.1666443655044969*RH*T*Log(H2SO4) + 0.9274525071251363*h6*RH*T*
     &  Log(H2SO4) + 0.00003323858217484082*q*RH*T*Log(H2SO4) - 
     &  0.021039758605777333*RH**2*T*Log(H2SO4) + 0.0018827939246916288*
     &  T**2*Log(H2SO4) + 0.0014616406839659635*h6*T**2*Log(H2SO4) - 
     &  2.8887052861955783e-7*q*T**2*Log(H2SO4) + 0.0007036899839828474*
     &  RH*T**2*Log(H2SO4) - 6.077113667352327e-6*T**3*Log(H2SO4) + 
     &  0.8758449734328851*Log(H2SO4)**2 - 61.130559299849466*h6*
     &  Log(H2SO4)**2 + 0.007311434341744119*q*Log(H2SO4)**2 - 
     &  0.4613762652439876*h6*q*Log(H2SO4)**2 + 0.11048996127121676*RH*
     &  Log(H2SO4)**2 - 5.241218416739315*h6*RH*Log(H2SO4)**2 - 
     &  0.0009849834685920686*q*RH*Log(H2SO4)**2 + 0.034417912409269905*
     &  RH**2*Log(H2SO4)**2 - 0.01737690919778605*T*Log(H2SO4)**2 - 
     &  0.10868884345327394*h6*T*Log(H2SO4)**2-0.000013815430650770167*
     &  q*T*Log(H2SO4)**2 - 0.0038903960432723123*RH*T*Log(H2SO4)**2 + 
     &  0.00005879779401379263*T**2*Log(H2SO4)**2 + 0.03393085910585771*
     &  Log(H2SO4)**3 + 2.2221405215918755*h6*Log(H2SO4)**3 - 
     &  0.00009520679794957912*q*Log(H2SO4)**3 + 0.0166181059051036*RH*
     &  Log(H2SO4)**3 - 0.00014027668903237805*T*Log(H2SO4)**3

	
	
C	write(7,100)h1,h2,h3,h4,h5,h6
 100  format(1p,6e12.4,2x)

        RETURN

        END SUBROUTINE LovejoyIN
C
!EOC
      END MODULE APM_NUCL_MOD
#endif
