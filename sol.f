! $Id: sol.f,v 1.1 2003/06/30 20:26:03 bmy Exp $
      SUBROUTINE SOL(NPTS,SUNCOS,ALT,SURFALT,TOTO3,CLOUDS,IDXAIR,IDXO3)
      
      ! References to F90 modules (bmy, 10/19/00)
      USE COMODE_MOD, ONLY : IXSAVE, IYSAVE, JLOP
      USE ERROR_MOD,  ONLY : ERROR_STOP

      IMPLICIT NONE
C*
C**********************************************************************
C                                                                     *
C  HARVARD TROPOSPHERIC CHEMISTRY MODULE FOR 3-D APPLICATIONS         *
C  by Larry Horowitz, Jinyou Liang, Gerry Gardner, Prof. Daniel Jacob *
C  of HARVARD UNIVERSITY    (Release V2.0)                            *
C**********************************************************************
C*                                                                    *
C*    Module contains:  SOL, SCATTR, MATINV, AIRMAS                   *
C*                                                                    *
C**********************************************************************
C*
C*    Compute (1) the column densities and (2) the J-values from the
C*    standard atmosphere level IHI down to ILO.  NLBATM is the lower
C*    boundary of the atmosphere (NLBATM .GE. 1).
C*
C***********************************************************************
C*
#     include "CMN_SIZE"
#     include "comode.h"
#     include "comsol.h"
*************************************************************

      INTEGER NPTS,IDXAIR(NLAT),IDXO3(NLAT)

      REAL*8 SUNCOS(MAXIJ),CLOUDS(MAXIJ,11)
      REAL*8 ALT(MAXIJ,NPVERT),SURFALT(MAXIJ)
      REAL*8 DM(NSTDL),DO3INT(NSTDL)
      REAL*8 ZZHT,AIRMAS

      INTEGER NW21,J,IALTL,IB,JJ,INDXO3,INDXAIR,ITEMP2,ITEMP3,JI,JB,
     1        IDOSR,IWL,IB1T6,IBND,NDOWL,ISR,NC,JLOOP,K,NW1,I,IIO2,
     2        IIO3,NW2,NWSRB,NLBATM

      INTEGER KLOOP,IJLOOP,IX,IY, IJWINDOW

      REAL*8 SCALE,GMU,TOTO3(NLAT),TB,WAVE,XQO2,XQO3,XLOW,XQAERT,
     1       XQAERM,TT,XHIG,XQRAY,XQAER,FINT,WSQI,REFRM1,PO2,
     2       ROVMG,TURBD0,TURBDX,CCETA,CCETA2,QNOLYA,QNOD10,QNOD00
C***********************************************************************
C*                
C*    Local absorption given by AERXCT(), column by AERINT(), product
C*    of OZONE column and cross-section given by DXO3IN, temperature 
C*    dependent cross-section by XBQO3(), parameters for fit to Bass
C*    OZONE data by XBASS, temperature dependent N2O5 cross-sections 
C*    by XQN2O5
C*                   
C***********************************************************************
C*
C*    Each S-R O2 band is divided into 6 intervals of constant opacity.
C*    These intervals are 0.00, 0.05, 0.25, 0.50, 0.75, 0.95, 1.00 
C*    A total of 16 S-R bands are used:
C*                              
C*      1.  1775-1783  T = 190 K           9.  1863-1882  T = 270 K
C*      2.  1783-1793  T = 190 K          10.  1882-1902  T = 270 K
C*      3.  1793-1804  T = 190 K          11.  1902-1925* T = 270 K
C*      4.  1804-1816  T = 190 K          12.  1925-1947  T = 230 K
C*      5.  1816-1831  T = 190 K          13.  1947-1972  T = 230 K
C*      6.  1816-1831  T = 270 K          14.  1972-1985  T = 230 K
C*      7.  1831-1846  T = 270 K          15.  1985-2000  T = 230 K
C*      8.  1846-1863  T = 270 K          16.  2000-2025  T = 230 K
C*             ******  Fang, Dalgarno and Wofsy paper ******
C*                          
C***********************************************************************
C*                  
C*    QNOLYA -- no X-section at LYMAN-ALPHA; assume all dissociation
C*    QNOD10 -- effective 15A cross-section (1816-1831A),
C*                 F(BETHKE) = 5.78E-3, DELTA(1,0)
C*    QNOD00 -- effective 23A cross-section (1902-1925A)
C*                 F(BETHKE) = 2.49E-3, DELTA(0,0)
C*    Assume DELTA(0,0) 100% dissociation, rescale DELTA(1,0) for
C*    solar flux asymmetry.
C*                         
C*    CCETA = 1.0E+8*(Speed of light)*(Planck's const)*(Avogadro's #)/
C*               [(Mean mole. wt.)*(Specific heat at constant-P)]
C*
C*
      REAL*8 TRANS(66),ABASS(16),BBASS(16),CBASS(16),AERINT(61),
     1       DXO3IN(61),XBQO3(61),DO3(61),DMINT(61),O2X(6,16),
     2       WK(6),WXSR(5,5),TAU,PITAU,FLTAU,RFLECT
      INTEGER ISRBEG(5),ISRFIN(5),NTT
      LOGICAL LTHIN,LSRBND,LDIFUS
C*
      COMMON /CCSCAT/ TAU(66),PITAU(66),FLTAU(66),RFLECT,NTT,LDIFUS

      !=================================================================
      ! /CCSCAT/ needs to be declared THREADPRIVATE for the 
      ! OpenMP parallelization for all platforms (bmy, 3/23/03)
      !=================================================================
!$OMP THREADPRIVATE( /CCSCAT/ )
C*
      DATA ISRBEG/1, 6, 9, 12, 14/, ISRFIN/5, 9, 11, 14, 16/
      DATA WK/0.05D0, 0.20D0, 0.25D0, 0.25D0, 0.20D0, 0.05D0/
C*                  
      DATA WXSR/0.16D0, 0.20D0, 0.22D0, 0.24D0, 0.18D0,
     1          0.12D0, 0.30D0, 0.34D0, 0.24D0, 0.00D0,
     2          0.14D0, 0.40D0, 0.46D0, 0.00D0, 0.00D0,
     3          0.44D0, 0.50D0, 0.06D0, 0.00D0, 0.00D0,
     4          0.20D0, 0.30D0, 0.50D0, 0.00D0, 0.00D0/
C*                  
      DATA O2X/ 4.33D-21,4.89D-21,6.63D-21,1.60D-20,7.20D-20,1.59D-18,
     1          2.10D-21,2.32D-21,3.02D-21,6.30D-21,3.46D-20,7.52D-19,
     2          5.95D-22,9.75D-22,2.53D-21,7.57D-21,7.38D-20,7.44D-19,
     3          3.33D-22,1.02D-21,4.09D-21,1.63D-20,8.79D-20,3.81D-19,
C*   Correction to original SR-band at (10-0), new F=1.60D-5
     4          5.00D-22,5.00D-22,2.00D-21,5.00D-21,6.00D-20,9.00D-19,
     5          6.50D-22,6.50D-22,3.00D-21,8.00D-21,6.00D-20,9.00D-19,
     6          3.90D-22,4.90D-22,9.49D-22,3.33D-21,2.14D-20,2.39D-19,
     7          1.29D-22,2.18D-22,8.28D-22,3.46D-21,1.94D-20,1.06D-19,
     8          6.26D-23,7.80D-23,2.62D-22,1.83D-21,1.25D-20,6.95D-20,
     9          2.74D-23,3.58D-23,8.64D-23,4.03D-22,2.13D-21,1.95D-20,
     A          1.95D-23,2.44D-23,4.89D-23,2.87D-22,1.95D-21,1.36D-20,
C*    New band model 1950-2025 based on lower continuum cross-sections
     B          1.50D-23,1.62D-23,2.37D-23,8.18D-23,6.45D-22,3.89D-21,
     C          1.00D-23,1.01D-23,1.07D-23,1.89D-23,1.26D-22,1.51D-21,
     D          0.90D-23,0.90D-23,0.92D-23,1.50D-23,4.81D-23,5.62D-22,
     E          0.80D-23,0.80D-23,0.80D-23,0.80D-23,0.80D-23,2.54D-23,
     F          0.80D-23,0.80D-23,0.80D-23,0.80D-23,0.80D-23,1.77D-23/
C*                         
      DATA QNOLYA/2.43D-38/, QNOD10/1.31D-17/, QNOD00/3.50D-18/
      DATA CCETA/4.11D7/, CCETA2/1.987D-8/
C*               
C*********  Start of data input for temperature dependent fit  *********
C*
C*    Data for temperature dependence of OZONE cross-sections from
C*    2675-3425, from RJS analysis of Bass data. Formulation is
C*      Q(L) = ABASS(L) + BBASS(L)*T(DEG.C) + CBASS(L)*T(DEG.C)**2
C*    WAVELENGTHS/2675,2725,2775,2825,2875,2925,2975,3O25,
C*                3075,3125,3175,3225,3275,3325,3375,3425/
C*                   
      DATA ABASS/ 8.6348D-18, 7.9559D-18, 4.7925D-18, 3.0765D-18,
     1            1.8787D-18, 1.0338D-18, 5.3570D-19, 2.7430D-19,
     2            1.3900D-19, 7.0600D-20, 3.5700D-20, 1.7423D-20,
     3            8.3145D-21, 4.0000D-21, 1.9576D-21, 8.8154D-22/
      DATA BBASS/ 2.9036D-22, 1.4341D-21, 1.1244D-21, 1.1463D-21,
     1            1.1953D-21, 9.3861D-22, 7.5860D-22, 4.8831D-22,
     2            3.2141D-22, 1.9707D-22, 1.1714D-22, 7.1821D-23,
     3            4.1105D-23, 2.4335D-23, 1.5765D-23, 1.0052D-23/
      DATA CBASS/-3.7953D-24, 1.1119D-23,-2.5398D-24,-3.9890D-24,
     1            2.7978D-24, 2.2221D-24, 3.6727D-24, 1.9728D-24,
     2            1.4252D-24, 1.0259D-24, 6.4958D-25, 3.7252D-25,
     3            2.0089D-25, 1.3979D-25, 1.0669D-25, 6.1010D-26/
C*                  
C**********  End of data input for temperature dependent fit  **********
C*   
C*    Redefine the column densities.
C*    ZZHT is scale height to give the column above the
C*    highest model layer. DO3INT(J), DMINT(J) and AERINT(J) are the
C*    column densities from layer J above (molecules cm^-2).
C*
C**** Read information on cross sections, actinic fluxes, O3 columns.
!**** Remove timing calls (bmy, 3/16/99)
	NW1        = 1
	NW2        = 70
	NWSRB      = 4
	NLBATM     = 1

C Get IIO2, IIO3 for TAU
	DO 100 I=1, MXSPE
	   IF('O2      '.EQ.CNAME(I)) IIO2=I
	   IF('O3      '.EQ.CNAME(I)) IIO3=I
  100 CONTINUE
C NSTDL = 41 standard layer, 0-2-4-...-80km
	NC= NSTDL
	TURBD0   = AERSOL(5)
	TURBDX   = AERSOL(6)
	PO2      = 0.2095D0
	DENCONS  = 6.02252D0+23.D0/28.966D0
	ROVMG    = 8.3144D0*1000.D0/(28.966D0*9.81D0)
C
C Reset fluxes to zero
C
      DO JLOOP=1, NPTS
	   DO K   = NW1, NW2+1
		ACTFLX(K, JLOOP)  = 0.D0
	   ENDDO
	ENDDO
C
! NOTE: Remove XSECT from SHARED list (bmy, 11/15/01)
!$OMP PARALLEL DO
!$OMP+SHARED( JLOP, SUNCOS, ACTFLX, ALT, NIJLOOP,
!$OMP+        IXSAVE, IYSAVE, TOTO3, IDXAIR, IDXO3, CLOUDS,
!$OMP+        SURFALT, STDO3, STDAIR, AERXCT, ABASS, BBASS, CBASS )
!$OMP+PRIVATE( IX, IY, IJLOOP, KLOOP, IB, JB, INDXO3, INDXAIR, 
!$OMP+         IJWINDOW, I, J, K, JI, ITEMP2, ITEMP3, IB1T6, IDOSR, 
!$OMP+         ISR, NDOWL, IBND, IALTL, NW21, JJ, IWL, LSRBND, LTHIN,
!$OMP+         XHIG, XLOW, TT, XQAERM, XQAERT, XQAER, WSQI, FINT, 
!$OMP+         XQRAY, REFRM1, TB, SCALE, XQO2, XQO3, WAVE, GMU, ZZHT,
!$OMP+         DMINT, DO3INT, AERINT, DO3, DM, XBQO3, DXO3IN, TRANS,
!$OMP+         NLBATM )
!$OMP+SCHEDULE( DYNAMIC )
	DO 250 IJLOOP = 1, NIJLOOP
         IX         = IXSAVE(IJLOOP)
         IY         = IYSAVE(IJLOOP)
c  Get the right index for SUNCOS, ALT, and SURFALT which are
c  calculated outside of chemistry module.
C  (This works for LEMBED= .TRUE. or .FALSE.)
         IJWINDOW   = (IY-1)*IM + IX

C   see if photolysis should be considered.
         IF (SUNCOS(IJWINDOW).LE.0.D0) GOTO 250
C
C ALT     passed in cm
C TOTO3() passed in (molec cm^-3)
C IDXAIR  passed INDEX FOR STANDARD TEMPERATURE PROFILE
C IDXO3   passed INDEX FOR STANDARD OZONE PROFILE
C
C   get indices for standard air density and ozone arrays
         INDXAIR    = IDXAIR(IY)
         INDXO3     = IDXO3(IY)
         ZZHT       = 4.0D5

C Calculate the lowest layer to use in radiation calculation, based
C on surface altitude.
         NLBATM = SURFALT(IJWINDOW)/2000. + 1.5
C        IF(NLBATM .NE. 1) write(6,*) IX,IY,NLBATM,SURFALT(IJWINDOW)

         IALTL      = ALT(IJWINDOW,NVERT)/2.D5 + NLBATM
C Signal the case when layer exceeds 80km
         IF (IALTL.GE. NSTDL) THEN
		   WRITE(*,*) 'LAYER IS TOO HIGH, ',
     +                    'ALT(IJWINDOW,NVERT)= ',
     +                     ALT(IJWINDOW,NVERT)/1.E5, 'KM'
                   CALL ERROR_STOP( 'STOP 81828', 'sol.f' )
         ENDIF

C Get DMINT, DO3INT for TAU, NC = NSTDL now.
         DO I=1,NC
            DO3(I)    = STDO3(I,INDXO3)
            DM (I)    = STDAIR(I,INDXAIR)
         ENDDO

         DMINT(NC)  = DM(NC)*ZZHT
         DO3INT(NC) = DO3(NC)*ZZHT
         AERINT(NC) = 1.0D-32

         DO I   = NC-1, 1, -1
            DMINT(I)  = DMINT(I+1) + 0.5D0*(DM(I)+DM(I+1))*2.D5
            DO3INT(I) = DO3INT(I+1)+ 0.5D0*(DO3(I)+DO3(I+1))*2.D5
            AERINT(I) = AERINT(I+1)+ 0.5D0*(AERXCT(I)+AERXCT(I+1))*2.D5
         ENDDO

	   SCALE      = TOTO3(IY)/DO3INT(NLBATM)
	   DO 130 I     = 1,NC
            DO3(I)    = DO3(I)*SCALE
            DO3INT(I) = DO3INT(I)*SCALE
  130    CONTINUE

C*
C*    K = NW2+1 is the Chappuis band at 6000 A, JNO3 is influenced.
C*	   NW1=1, NW2=70
C For scattr, SUNCOS is passed
	   GMU      = SUNCOS(IJWINDOW)
 
	   NW21 = NW2+1
	   DO 240 K=NW1,NW21
            IF (NW1.EQ.0) WRITE (*,140)
  140       FORMAT(' NW1=0 IN SOLCHEM: CODE WILL BOMB.')
		IF (K .LE. NW2) GO TO 150
C*
C*    Optically thin photolysis in visible is treated as Chappuis-band
C*    XQO2,XQO3 are the absorption cross-sections for O2 and O3.
C*    CLOUDS(IJLOOP,1) is the surface albedo.

		WAVE = 6000.D0 
		XQO2 = 0.0D0
C           XQO3 = 4.50D-21
C Changed XQO3 to represent the average O3 Chappuis-band
C         cross-section from 500-700nm, rather than just using the
C         600nm peak value
            XQO3 = 2.89D-21
		LTHIN = .TRUE.
		GO TO 180

  150       WAVE = WL(K)
C*                    
C*    Temperature dependence of OZONE cross-sections
C*    For wavelengths 2625-3425, based on RJS' analysis of Bass data.
C*    Note: the integral over altitude (OZONE*CROSS-SECTION(T))
C*    replaces DO3INT for these wavelengths.
C*    
		IF ((K.LT.23)  .OR.  (K.GT.38)) GO TO 170
		IB = K-22
C*       
C*    Initialize arrays XBQO3 and DXO3IN at top of atmosphere

		TB = TATM(NC, INDXAIR ) - 273.15D0
		XBQO3(NC) = ABASS(IB) + BBASS(IB)*TB + CBASS(IB)*TB*TB
		DXO3IN(NC) = XBQO3(NC)*DO3INT(NC)
		DO 160 JJ = 2,NC-NLBATM+1
               JB = NC + 1 -JJ
		   TB = TATM(JB, INDXAIR ) - 273.15D0
		   XBQO3(JB) = ABASS(IB) + BBASS(IB)*TB +
     x                              CBASS(IB)*TB*TB
               DXO3IN(JB) = DXO3IN(JB+1) + 0.25D0*(2.D5)*
     x                     (DO3(JB)+DO3(JB+1))*(XBQO3(JB+1)+XBQO3(JB))
  160       CONTINUE
  170       CONTINUE
C************************************************************
            LTHIN = .FALSE.
  180       CONTINUE
C*                    
C*    Rayleigh + Raman cross-section (include for all wavelengths 
C*    except the S-R band)
C*    XQRAY is the Rayleigh absorption optical thickness, calculated 
C*    with a NASA Handbook formula.
C*    
C*    TURBD0 = AERSOL(5) -- total optical depth from aerosols at 310 nm.
C*    TURBDX = AERSOL(6) -- wavelength dependence of aerosol absorption
C*    XQAER = total aerosol absorption at wavelength WAVE,
C*            normalized to the aerosol column.
C*            
		LDIFUS = .TRUE.
		WSQI = 1.D8/(WAVE*WAVE)
		REFRM1 = 1.0D-6*(64.328D0+29498.1D0/(146.D0-WSQI)+
     +		            255.4D0/(41.D0-WSQI))
		XQRAY = 5.40D-21*(REFRM1*WSQI)**2
		XQAER = (TURBD0* (WAVE/3100.D0)**TURBDX )/AERINT(NLBATM)
C*                       
C*    Locate the Oxygen S-R bands
		LSRBND = .FALSE.
		FINT = 1.0D0
		IBND = 0
		NDOWL = 1
C*                
C*    NWSRB is the lowest wavelength of the S-R band, which then extends
C*    over the next 5 wavelength ranges.
C*            
		ISR = K - NWSRB
		IF (ISR .LT. 1  .OR.  ISR .GT. 5) GO TO 190
C*                 
C*    In the S-R band, ignore the effect of clouds (LDIFUS = .FALSE.)
C*    and of Rayleigh scattering.

		LSRBND = .TRUE.
		LDIFUS = .FALSE.
		XQRAY = 0.0D0
		IDOSR = ISRBEG(ISR)
		NDOWL = 6*(ISRFIN(ISR)+1-IDOSR)
C*                  
C*    NDOWL .GT. 1, only for S-R bands.

  190		DO 230 IWL=1,NDOWL
		   IF(.NOT.LSRBND) GO TO 200
		   IB1T6 = MOD(IWL-1,6) + 1
		   IBND = (IWL-1)/6 + IDOSR
		   XQO2 = O2X(IB1T6,IBND)
		   FINT = WXSR(IBND+1-IDOSR,ISR)*WK(IB1T6)
  200		CONTINUE
C*                         
C*    Set up optical depth scale, where TAU(J) is the optical depth at 
C*    altitude J.  TAU(1) = 0.0 at the top of atmosphere, TAU(61) is the
C*    bottom of the atmosphere (ground).
C*  
C*    NLBATM is the lowest layer of the atmosphere, set = 1 above
C*    XQAERT is the optical depth at altitude J from the aerosol column
C*    XQAERM is the contribution of layer J to the aerosol optical depth
C*                    
		   NTT = NC + 2 - NLBATM
		   DO 210 J=NLBATM,NC
			JI = NC + 2 - J
			XQAERT = XQAER*AERINT(J)
			XQAERM = XQAER*AERXCT(J)
			TT     = TATM(J, INDXAIR)
C Get index of O3 and O2 in XSECT, and TT in TARRAY.
C Only 1 T in 8col.dat
			ITEMP3= 1
			ITEMP2= 1
C The IIO3, IIO2, ITEMP3, ITEMP2 are ready. (local variable.)

			IF( K .LE. NW2 .AND. .NOT.LSRBND)
     x			      XQO2  = XSECT(K, IIO2, ITEMP2, 1)
C*              
C*    Temperature dependence on OZONE cross-section  
C*           
C*    TAU(JI) is the sum of the contributions from Rayleigh scattering,
C*    the OZONE column, and the aerosol column above.
C*    PITAU(JI) is the single-scattering albedo.
C*       
		      IF ((K.LT.23)  .OR.  (K.GT.38)) THEN
	               IF ( K .LE. NW2) XQO3  = XSECT(K, IIO3, ITEMP3, 1)
                        TAU(JI) = XQRAY*DMINT(J) + XQO3*DO3INT(J) +
     x                            XQO2*PO2*DMINT(J) + XQAERT
                        PITAU(JI) = (XQRAY*DM(J))/
     x                   (XQRAY*DM(J)+XQAERM+XQO3*DO3(J)+XQO2*PO2*DM(J))
			   ELSE 
                        TAU(JI)=XQRAY*DMINT(J)+DXO3IN(J)+
     x                          XQO2*PO2*DMINT(J)+XQAERT
                        PITAU(JI) = (XQRAY*DM(J))/
     x                   (XQRAY*DM(J)+XQAERM+XBQO3(J)*DO3(J)+
     x                   XQO2*PO2*DM(J))
			   ENDIF 

  210          CONTINUE
		   TAU(1) = 0.0D0
		   PITAU(1) = PITAU(2)
C*              
C*******************************************************************
               RFLECT = CLOUDS(IJLOOP,1)
               CALL SCATTR (IJLOOP, GMU, CLOUDS, ZZHT)
C*******************************************************************
C*                          
C*    Transfer mean radiation field and rescale solar flux   

               DO 220 J=NLBATM,NC
                  JI = NC + 2 - J
                  TRANS(J) = FLTAU(JI)
  220          CONTINUE 
C Save ACTFLX for grid boxes, use linear interpolation vertically.
C The grid box is between standard levels IALTL and IALTL+1.
               DO 225 J=1,NPVERT
                  IALTL      = ALT(IJWINDOW,J)/2.D5 + NLBATM
                  XLOW       = ALT(IJWINDOW,J) - 2.D5*(IALTL-NLBATM)
                  XHIG       = 2.D5 - XLOW
                  KLOOP    = JLOP(IX,IY,J)
                  IF (KLOOP.EQ.0) GOTO 225
                     IF (LTHIN) THEN
                         ACTFLX(K, KLOOP)= (TRANS(IALTL)*XHIG +
     2                                      TRANS(IALTL+1)* XLOW)/2.D5
                     ELSE
                         ACTFLX(K, KLOOP)= ACTFLX(K, KLOOP) +
     2                                     FL(K)*FINT*
     3                   (TRANS(IALTL)*XHIG + TRANS(IALTL+1)*XLOW)/2.D5
                     ENDIF
  225          CONTINUE  !(J)
C
C     TRANS(1) is the ratio of ground 6000A radiation intensity to the 
C     intensity of 6000A solar radiation; save to use when calculating 
C     the light-dependent ISOPRENE emission rate in -SETUPR-
  230       CONTINUE  !(IWL)
  240    CONTINUE  !(K)

  250 CONTINUE  !(IJLOOP)
!$OMP END PARALLEL DO

      RETURN
      END
