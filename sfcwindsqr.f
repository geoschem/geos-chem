C $Id: sfcwindsqr.f,v 1.1 2003/06/30 20:26:08 bmy Exp $
      REAL*8 FUNCTION SFCWINDSQR( I, J ) 
!
!******************************************************************************
!  Function SFCWINDSQR computes the surface wind squared from the DAO 
!  U and V winds at 10 m above the surface.  (bmy, 12/21/98, 11/26/01)
!
!  NOTES:
!  (1 ) The old SFCWINDSQR computed the surface wind squared (m/s)^2 from the
!        the Harvard CTM winds (kg/s).  But since the DAO winds are already
!        in units of (m/s) then the previous unit conversion is unnecessary
!        and costly in terms of computer resources.  
!  (2 ) Since GEOS-1 has U and V at 10 m, these are more representative
!        of the surface than UWND(I,J,1) and VWND(I,J,1).  
!  (3 ) Pass GEOS-1 U10M and V10M fields via CMN_UV10M so that the argument 
!        list does not have to be modified in several existing Harvard CTM 
!        subroutines.
!  (4 ) GEOS-STRAT does not store U10M and V10M, so compute 10 m wind speed
!        from UWND(I,J,1) and VWND(I,J,1) in MAKE_WIND10M.
!  (5 ) Now check for NaN's (bmy, 4/27/00)
!  (6 ) Now reference U10M and V10M from "dao_mod.f" instead of from
!        common block header files "CMN_UV10M".  Also extend code
!        to GEOS-2 and GEOS-3 met fields. (bmy, 7/11/00)
!  (7 ) Now use interface IT_IS_NAN (from "error_mod.f") to trap NaN's.
!        This will work on DEC/Compaq and SGI platforms. (bmy, 3/8/01)
!  (8 ) Now call CHECK_VALUE from "error_mod.f".  This will test SFCWINDSQR
!        for NaN or Infinity conditions.  Also updated comments and made
!        cosmetic changes. (bmy, 7/16/01)
!  (9 ) Removed obsolete, commented-out code from 7/01 (bmy, 11/26/01)
!******************************************************************************
!
      ! References to F90 modules
      USE DAO_MOD,   ONLY : U10M, V10M, WIND_10M, USE_WIND_10M
      USE ERROR_MOD, ONLY : CHECK_VALUE

      IMPLICIT NONE
        
#     include "CMN_SIZE"

      ! Arguments
      INTEGER, INTENT(IN) :: I, J

      !=================================================================
      ! SFCWINDSQR begins here!!
      !
      ! If using GEOS-STRAT fields, then square the value of 
      ! WIND_10M as computed from the instantaneous U and V winds 
      ! in MAKE_WIND10M.
      !
      ! Otherwise (i.e. for GEOS-1, GEOS-2, and GEOS-3 met data), 
      ! square the U10M and V10M wind vectors and take that as the 
      ! surface wind speed.
      !=================================================================
#if   defined ( GEOS_STRAT )
      IF ( USE_WIND_10M ) THEN
         SFCWINDSQR = WIND_10M(I,J)**2
      ELSE
         SFCWINDSQR = U10M(I,J)**2 + V10M(I,J)**2
      ENDIF

#else
      SFCWINDSQR = U10M(I,J)**2 + V10M(I,J)**2

#endif

      ! Now check for NaN or Infinity (bmy, 7/16/01)
      CALL CHECK_VALUE( SFCWINDSQR, (/I,J,0,0/), 
     &                 'SFCWINDSQR', 'at sfcwindsqr.f' )

      ! Return to calling program
      END FUNCTION SFCWINDSQR
