! $Id: photoj.f,v 1.1 2003/06/30 20:26:03 bmy Exp $
      SUBROUTINE PHOTOJ( NSLON, NSLAT, YLAT, IDAY, MONTH,   CSZA,  
     &                   P,     T,     OD,   SA,   OPTDUST, OPTAER )
C-----------------------------------------------------------------------
C----jv_trop.f:  new FAST J-Value code, troposphere only (mjprather 6/96)
C----     uses special wavelength quadrature spectral data (jv_spec.dat)
C---      that includes only 289 nm - 800 nm  (later a single 205 nm add-on)
C---      uses special compact Mie code based on Feautrier/Auer/Prather vers.
C-----------------------------------------------------------------------
C  Add the following input variables for CTM interface (bmy, rvm, 9/30/00)
C
C  Variable  Type    Dimension   Units   Description
C  --------  ----    ---------   -----   -----------
C  NSLON     int        -          -     Longitude index
C  NSLAT     int        -          -     Latitude index
C  YLAT      dble       -         [deg]  Latitude value corresponding to NSLAT
C  IDAY      int        -          -     Day of Year (0-365 or 0-366)
C  MONTH     int        -          -     Month of Year (1-12)
C  CSZA      dble       -          -     Cosine of solar zenith angle
C  P         dble       -         [mb]   Surface Pressure - Model top pressure
C  T         dble     [LPAR]      [K]    Vertical temperature profile
C  OD        dble     [LPAR]       -     Vertical optical depth profile
C  SA        dble       -          -     Surface Albedo
C  OPTDUST   dble     [LPAR,NDUST] -     Dust optical depth (rvm, 9/30/00)
C  OPTAER    dble  [LPAR,NAER*NRH] -     Aerosol optical depth (rvm, 2/27/02)
C-----------------------------------------------------------------------
c
c     zpj      External array providing J-values to main CTM code
c     timej    Offset in hours from start of timestep to time J-values
c              required for - take as half timestep for mid-step Js.
c     solf     Solar distance factor, for scaling; normally given by:
c                      1.0-(0.034*cos(real(iday-172)*2.0*pi/365.))
c
C----------basic common blocks:-----------------------------------------
      IMPLICIT NONE

#     include "cmn_fj.h"
#     include "jv_cmn.h"

C=============== INPUT PARAMETERS ======================================
      INTEGER, INTENT(IN)    :: MONTH, NSLAT,   NSLON,  IDAY
      REAL*8,  INTENT(IN)    :: YLAT,  T(LPAR), OD(LPAR) 
      REAL*8,  INTENT(IN)    :: CSZA,  P,       SA
      REAL*8,  INTENT(INOUT) :: OPTDUST(LPAR,NDUST)   !(rvm, bmy, 9/30/00)
      REAL*8,  INTENT(INOUT) :: OPTAER(LPAR,NAER*NRH) !(rvm, bmy, 2/27/02)

C=============== LOCAL VARIABLES =======================================
      INTEGER             :: I, J
      REAL*8, PARAMETER   :: PI = 3.14159265358979324D0

C-----------------------------------------------------------------------
C---Initialize ZJ and ZPJ
      DO I=1,JPNL
         DO J=1,JPPJ
            ZJ(I,J)=0.D0
            ZPJ(I,J,NSLON,NSLAT)=0.D0
         ENDDO
      ENDDO
C
C---Import the cosine of the SZA from the CTM (bmy, 9/10/99)
      U0  = CSZA
      SZA = ACOS(CSZA) * ( 180.0d0 / PI )
C-----------------------------------------------------------------------
C---If you want to set SZA = 0 degrees for testing,
C---then uncomment the following lines (bmy, 9/13/99) 
C      U0  = 1.0d0
C      SZA = 0.0d0
C-----------------------------------------------------------------------
      IF(SZA.GT.SZAMAX) RETURN
C
C---Set up profiles on model levels
      CALL SET_PROF( P, T, OD, SA, MONTH, YLAT, OPTDUST, OPTAER )
C-----------------------------------------------------------------------
      CALL JVALUE( SA )
C-----------------------------------------------------------------------
      CALL JRATET( T, IDAY )
C-----------------------------------------------------------------------
C
C  "zj" updated in JRATET - pass this back to ASAD as "zpj"
      DO I=1,JPNL
         DO J=1,JPPJ
            ZPJ(I,J,NSLON,NSLAT)= ZJ(I,J)
         ENDDO
      ENDDO

      ! Return to calling program
      END SUBROUTINE PHOTOJ
