!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: tcorr
!
! !DESCRIPTION: Function TCORR applies the temperature correction for 
!  isoprene emissions, according to Guenther et al.(92)
!\\
!\\
! !INTERFACE:
!
      FUNCTION TCORR( TEMP )
!
! !USES:
!
      IMPLICIT NONE
#     include "define.h"
!
! !INPUT PARAMETERS: 
!
      REAL*8, INTENT(IN) :: TEMP    ! Temperature [K] 
!
! !RETURN VALUE:
!
      REAL*8             :: TCORR   ! Corrected temp for ISOP emissions [K]
!
! !REMARKS:
!  References:
!  ============================================================================
!  Guenther et al, 1992, ... 
!
! !REVISION HISTORY: 
!  15 Nov 1993 - Y. H. Wang  - Initial version
!  (1 ) Removed DATA statements, replaced w/ F90 syntax.  Updated comments
!        and made cosmetic changes (bmy, 4/4/03)
!  19 Nov 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
      REAL*8, PARAMETER  :: R   = 8.314
      REAL*8, PARAMETER  :: CT1 = 95000.
      REAL*8, PARAMETER  :: CT2 = 230000.
      REAL*8, PARAMETER  :: T1  = 303.
      REAL*8, PARAMETER  :: T3  = 314.
      
      !=================================================================
      ! TCORR begins here!
      !=================================================================
      TCORR =
     &     EXP( CT1/(R*T1*TEMP) * (TEMP-T1) ) /
     &     ( 1 + EXP( CT2/(R*T1*TEMP) * (TEMP-T3) ) )

      END FUNCTION TCORR
!EOC
