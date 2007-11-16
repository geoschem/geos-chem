! $Id: photoj.f,v 1.4 2007/11/16 18:47:44 bmy Exp $
!---------------------------------------------------------------------------
! Prior to 2/13/07:
!      SUBROUTINE PHOTOJ( NLON, NLAT, YLAT, DAY_OF_YR, MONTH,   DAY,  
!     &                   CSZA, P,    T,    SA,  OD,   OPTDUST, OPTAER )
!---------------------------------------------------------------------------
      SUBROUTINE PHOTOJ( NLON, NLAT, YLAT, DAY_OF_YR, MONTH,   DAY,  
     &                   CSZA, T,    SA,   OD,        OPTDUST, OPTAER )
!
!******************************************************************************
!  Subroutine PHOTOJ is the driver routine for the FAST-J photolysis package.
!  (Oliver Wild & Michael Prather, 1996, 1999, bmy, 7/18/03, 2/13/07)
!
!  New FAST J-Value code, troposphere only (mjprather 6/96); uses special 
!  wavelength quadrature spectral data (jv_spec.dat) that includes only 
!  289 nm - 800 nm  (later a single 205 nm add-on); uses special compact Mie 
!  code based on Feautrier/Auer/Prather vers.
!
!  Arguments as Input: 
!  ============================================================================
!  (1 ) NLON      (INTEGER) : Grid box longitude index            [unitless]
!  (2 ) NLAT      (INTEGER) : Grid box latitude index             [unitless]
!  (3 ) YLAT      (REAL*8 ) : Grid box latitude                   [degrees]
!  (4 ) DAY_OF_YR (INTEGER) : Current day of year                 [1-366]
!  (5 ) MONTH     (INTEGER) : Current month                       [1-12]
!  (6 ) DAY       (INTEGER) : Current day of month                [
!  (7 ) CSZA      (REAL*8 ) : Cosine of solar zenith angle        [unitless]
!  (8 ) P         (REAL*8 ) : Psurface - PTOP @ (NSLON,NSLAT)     [hPa]
!  (9 ) T         (REAL*8 ) : Vertical temperature profile        [K]
!  (10) SA        (REAL*8 ) : Surface albedo @ (NSLON,NSLAT)      [unitless]
!  (11) OD        (REAL*8 ) : Vertical optical depth profile      [unitless]
!  (12) OPTDUST   (REAL*8 ) : Vertical dust optical depth profile [unitless]
!  (13) OPTAER    (REAL*8 ) : Vertical aerosol opt. depth profile [unitless]
!
!  Important variables from "jv_cmn.h"
!  ============================================================================
!  (1 ) ZJ      (REAL*8 ) : Column array for J-values 
!  (2 ) ZPJ     (REAL*8 ) : Global array for J-values (passed to SMVGEAR)
!  (3 ) JPNL    (INTEGER) : # of GEOS-CHEM layers in which to compute J-values
!  (4 ) JPPJ    (INTEGER) : # of photolysis rxns for FAST-J
!
!  NOTES:
!  (1 ) Renamed NSLON to NLON and NSLAT to NLAT.  Now add DAY_OF_YR 
!        (formerly IDAY) and DAY to the arg list.  Swap places in arg list 
!        of SA and OD.  Now pass NLON, NLAT, DAY_OF_YR and DAY to "set_prof.f".
!        Added standard documentation header; cosmetic changes. (bmy, 7/15/03)
!  (2 ) We don't need to pass "P" via the arg list (bmy, 2/13/07)
!******************************************************************************
!
      IMPLICIT NONE

#     include "cmn_fj.h"     ! IPAR, JPAR, LPAR, JPNL, JPPJ
#     include "jv_cmn.h"     ! ZJ, ZPJ

      ! Arguments
      INTEGER, INTENT(IN)    :: DAY,   DAY_OF_YR, MONTH
      INTEGER, INTENT(IN)    :: NLAT,  NLON
      REAL*8,  INTENT(IN)    :: YLAT,  T(LPAR),   OD(LPAR) 
      !------------------------------------------------------------
      ! Prior to 2/13/07:
      !REAL*8,  INTENT(IN)    :: CSZA,  P,         SA
      !------------------------------------------------------------
      REAL*8,  INTENT(IN)    :: CSZA,  SA
      REAL*8,  INTENT(INOUT) :: OPTDUST(LPAR,NDUST)   !(rvm, bmy, 9/30/00)
      REAL*8,  INTENT(INOUT) :: OPTAER(LPAR,NAER*NRH) !(rvm, bmy, 2/27/02)

      ! Local variables
      INTEGER                :: I, J
      REAL*8, PARAMETER      :: PI = 3.14159265358979324D0

      !=================================================================
      ! PHOTOJ begins here!
      !=================================================================

      ! Zero ZJ (column J-value array) and ZPJ (global J-value array)
      DO I = 1, JPNL
      DO J = 1, JPPJ
         ZJ(I,J)            = 0.D0
         ZPJ(I,J,NLON,NLAT) = 0.D0
      ENDDO
      ENDDO

      ! Import the cosine of the SZA from the CTM (bmy, 9/10/99)
      U0  = CSZA
      SZA = ACOS(CSZA) * ( 180.0d0 / PI )

      !-----------------------------------------------------------------
      !### If you want to set SZA = 0 degrees for testing,
      !### then uncomment the following lines (bmy, 9/13/99) 
      !U0  = 1.0d0
      !SZA = 0.0d0
      !-----------------------------------------------------------------

      ! Return if sun is below the horizon
      IF ( SZA > SZAMAX ) RETURN

      ! Set up Air, O3, BC profiles on GEOS-CHEM vertical levels
!---------------------------------------------------------------------------
! Prior to 2/13/07:
! We no longer need to pass P to SET_PROF, as we have initialized
! the PJ array in "fast_j.f" (bmy, 2/13/07)
!      CALL SET_PROF( NLON, NLAT, YLAT, MONTH, DAY, 
!     &               P,    T,    SA,   OD,    OPTDUST, OPTAER )
!---------------------------------------------------------------------------
      CALL SET_PROF( NLON, NLAT, YLAT, MONTH,   DAY, 
     &               T,    SA,   OD,   OPTDUST, OPTAER )

      ! Compute actinic flux at each GEOS-CHEM vertical level
      CALL JVALUE( SA )
      
      ! Calculate J-values for all species
      CALL JRATET( T, DAY_OF_YR )

      ! ZJ is the J-value array for this column only
      ! Store in ZPJ (global array) for passing to SMVGEAR
      DO I = 1, JPNL
      DO J = 1, JPPJ
         ZPJ(I,J,NLON,NLAT) = ZJ(I,J)
      ENDDO
      ENDDO

      ! Return to calling program
      END SUBROUTINE PHOTOJ
