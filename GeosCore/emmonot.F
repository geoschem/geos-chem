!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: emmonot
!
! !DESCRIPTION: Subroutine EMMONOT computes the BIOGENIC MONOTERPENE 
!  EMISSIONS for each grid box in units of [atoms C/box/step].
!\\
!\\
! !INTERFACE:
!
      FUNCTION EMMONOT( IJLOOP, TMMP, XNUMOL )
!
! !USES:
!
      USE CMN_SIZE_MOD              ! Size parameters
      USE CMN_MONOT_MOD             ! BASEMONOT
      USE CMN_VEL_MOD               ! XYLAI, IJREG, IJLAND, IJUSE

      IMPLICIT NONE
#     include "define.h"

!
! !INPUT PARAMETERS: 
!
      INTEGER, INTENT(IN) :: IJLOOP   ! 1-D grid box index
      REAL*8,  INTENT(IN) :: TMMP     ! Local air temperature (K)
      REAL*8,  INTENT(IN) :: XNUMOL   ! Number of atoms C / kg C 
!
! !RETURN VALUE:
!
      REAL*8              :: EMMONOT 
!
! !REMARKS:
!  Important Common Block Variables:
!  ============================================================================
!  (1 ) XYLAI     (F77_CMN_VEL ) : Leaf Area Index of land type for current MONTH
!  (2 ) IJREG     (F77_CMN_VEL ) : Number of Olson land types per grid box
!  (3 ) IJLAND+1  (F77_CMN_VEL ) : Olson land type index
!  (4 ) IJUSE     (F77_CMN_VEL ) : Olson land type fraction per box (in mils)
!  (5 ) BASEMONOT (F77_CMN_ISOP) : Baseline MONOTERPENE emissions [kg C/box/step]
! 
! !REVISION HISTORY: 
!  04 Sep 2001 - Y. H. Wang, B. Field, R. Yantosca  - Initial version
!  (1 ) Now use F90 syntax.  Use "D" exponents to force double precision.
!        Updated comments, and mad cosmetic changes (bmy, 9/4/01) 
!  (2 ) Removed obsolete, commented-out code from 8/01 (bmy, 11/26/01)
!  02 Dec 2010 - R. Yantosca - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER             :: INVEG
!
! !DEFINED PARAMETERS:
!
      REAL*8,  PARAMETER  :: TS=303d0, BETA=0.09d0

      !=================================================================
      ! EMMONOT begins here!
      !=================================================================

      ! Initialize
      EMMONOT = 0d0

      ! Loop over all Olson land types in this grid box
      DO INVEG = 1, IJREG(IJLOOP)

         ! Compute monoterpene emissions at box IJLOOP in [kg C/box/step].  
         ! Monoterpenes are now scaled to leaf area index XYLAI.  Also 
         ! multiply by the fraction of grid box IJLOOP occupied
         ! by this Olson land type. (bdf, bmy, 8/2/01)
         EMMONOT = EMMONOT + 
     &             ( BASEMONOT(IJLOOP,INVEG) * XYLAI(IJLOOP,INVEG) * 
     &               DBLE( IJUSE(IJLOOP,INVEG) ) / 1000d0 )

      ENDDO

      !=================================================================
      ! Temperature correction from Guenther et al. (1995)
      ! BETA is an empirical coefficient given by Guenther. (.09 K-1)
      ! TS is leaf temperature at standard conditions, (303 K)
      ! foliar density is accounted for in monotemis.table.
      !=================================================================

      ! Temp-corrected MONOTERPENE emissions in [kg C/box/step]
      EMMONOT = EMMONOT * EXP( BETA * ( TMMP - TS ) )

      ! Convert MONOTERPENE emissions to [atoms C/box/step]
      EMMONOT = EMMONOT * XNUMOL

      END FUNCTION EMMONOT
!EOC
