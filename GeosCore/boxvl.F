!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: boxvl
!
! !DESCRIPTION: The new function BOXVL converts the DAO grid box volume 
!  values stored in AIRVOL from m3 to cm3.  The conversion factor is 
!  (100)**3 = 1e6 cm3 per m3.
!\\
!\\
! !INTERFACE:
!
      REAL*8 FUNCTION BOXVL( I, J, L, State_Met )
!
! !USES:
!
      USE GIGC_State_Met_Mod, ONLY : MetState

      IMPLICIT NONE
!
! !INPUT PARAMETERS: 
!
      INTEGER,        INTENT(IN) :: I           ! Longitude index
      INTEGER,        INTENT(IN) :: J           ! Latitude index
      INTEGER,        INTENT(IN) :: L           ! Level index
      TYPE(MetState), INTENT(IN) :: State_Met   ! Meteorology State object
!
! !REVISION HISTORY: 
!  30 Jan 1998 - R. Yantosca - Initial version
!  (1 ) CMN_VOL is used to pass AIRVOL.
!  (2 ) Use C-preprocessor #include statement to include CMN_SIZE, which 
!        has IIPAR, JJPAR, LLPAR, IIPAR, JJPAR, LGLOB. 
!  (3 ) Now use F90 syntax for declarations (bmy, 10/5/99)
!  (4 ) Now reference AIRVOL from "dao_mod.f" instead of from common
!        block header file "CMN_VOL". (bmy, 6/26/00)
!  (5 ) Removed obsolete code from 6/26/00 (bmy, 8/31/00)
!  (6 ) Updated comments (bmy, 8/5/02)
!  02 Dec 2010 - R. Yantosca - Initial version
!  09 Nov 2012 - M. Payer    - Replaced all met field arrays with State_Met
!                              derived type object
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      !=================================================================
      ! BOXVL begins here!
      !=================================================================
      BOXVL = State_Met%AIRVOL(I,J,L) * 1d6

      END FUNCTION BOXVL
!EOC
