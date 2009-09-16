! $Id: boxvl.f,v 1.1 2009/09/16 14:06:39 bmy Exp $
      REAL*8 FUNCTION BOXVL (I,J,L)
!
!*****************************************************************************
!  The new function BOXVL converts the DAO grid box volume values stored
!  in AIRVOL from m^3 to cm^3.  The conversion factor is (100)^3 = 10^6 
!  cm^3 per m^3.  (bmy, 1/30/98, 8/5/02)
!
!  NOTES:
!  (1 ) CMN_VOL is used to pass AIRVOL.
!  (2 ) Use C-preprocessor #include statement to include CMN_SIZE, which 
!        has IIPAR, JJPAR, LLPAR, IGLOB, JGLOB, LGLOB. 
!  (3 ) Now use F90 syntax for declarations (bmy, 10/5/99)
!  (4 ) Now reference AIRVOL from "dao_mod.f" instead of from common
!        block header file "CMN_VOL". (bmy, 6/26/00)
!  (5 ) Removed obsolete code from 6/26/00 (bmy, 8/31/00)
!  (6 ) Updated comments (bmy, 8/5/02)
!*****************************************************************************
!
      ! References to F90 modules
      USE DAO_MOD, ONLY : AIRVOL

      IMPLICIT NONE

      ! Arguments
      INTEGER, INTENT(IN) :: I, J, L

      ! BOXVL begins here!
      BOXVL = AIRVOL(I,J,L) * 1d6

      ! Return to calling program
      END FUNCTION BOXVL

