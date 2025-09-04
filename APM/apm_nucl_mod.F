#ifdef APM
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
      PUBLIC  :: IONRATE0

!APM2+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      PUBLIC :: IONRATE
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
C
!EOC
      END MODULE APM_NUCL_MOD
#endif
