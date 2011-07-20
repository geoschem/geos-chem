!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: sfcwindsqr 
!
! !DESCRIPTION: Function SFCWINDSQR computes the surface wind squared 
!  from the U and V winds at 10 m above the surface.
!\\
!\\
! !INTERFACE:
!
      REAL*8 FUNCTION SFCWINDSQR( I, J ) 
!
! !USES:
!
      USE DAO_MOD, ONLY : U10M, V10M

      IMPLICIT NONE
!  
! !INPUT PARAMETERS: 
!
      INTEGER, INTENT(IN) :: I   ! Longitude index
      INTEGER, INTENT(IN) :: J   ! Latitude index
! 
! !REVISION HISTORY: 
!  21 Dec 1998 - R. Yantosca - Initial version
!  (1 ) The old SFCWINDSQR computed the surface wind squared (m/s)^2 from the
!        the Harvard CTM winds (kg/s).  But since the DAO winds are already
!        in units of (m/s) then the previous unit conversion is unnecessary
!        and costly in terms of computer resources.  
!  (2 ) Since GEOS-1 has U and V at 10 m, these are more representative
!        of the surface than UWND(I,J,1) and VWND(I,J,1).  
!  (3 ) Pass GEOS-1 U10M and V10M fields via F77_CMN_UV10M so that the argument 
!        list does not have to be modified in several existing Harvard CTM 
!        subroutines.
!  (4 ) GEOS-STRAT does not store U10M and V10M, so compute 10 m wind speed
!        from UWND(I,J,1) and VWND(I,J,1) in MAKE_WIND10M.
!  (5 ) Now check for NaN's (bmy, 4/27/00)
!  (6 ) Now reference U10M and V10M from "dao_mod.f" instead of from
!        common block header files "F77_CMN_UV10M".  Also extend code
!        to GEOS-2 and GEOS-3 met fields. (bmy, 7/11/00)
!  (7 ) Now use interface IT_IS_NAN (from "error_mod.f") to trap NaN's.
!        This will work on DEC/Compaq and SGI platforms. (bmy, 3/8/01)
!  (8 ) Now call CHECK_VALUE from "error_mod.f".  This will test SFCWINDSQR
!        for NaN or Infinity conditions.  Also updated comments and made
!        cosmetic changes. (bmy, 7/16/01)
!  (9 ) Removed obsolete, commented-out code from 7/01 (bmy, 11/26/01)
!  (10) Remove support for GEOS-1 and GEOS-STRAT met fields.  Also remove
!        call to CHECK_VALUE. (bmy, 8/4/06)
!  08 Dec 2009 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
      !=================================================================
      ! SFCWINDSQR begins here!!
      !=================================================================

      ! Take the 10m wind speed squared as sfc wind speed squared
      SFCWINDSQR = U10M(I,J)**2 + V10M(I,J)**2

      END FUNCTION SFCWINDSQR
!EOC
