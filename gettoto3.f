! $Id: gettoto3.f,v 1.1 2003/06/30 20:26:06 bmy Exp $
      SUBROUTINE GETTOTO3( TOTO3, IDXAIR, IDXO3 )
!
!******************************************************************************
!  Subroutine GETTOTO3 looks up the appropriate total ozone column in O3DU 
!  (from o3du.dat) for a given latitude and month (interpolating to GEOS-CHEM
!  grid box latitudes).  GETTOTO3 is only called for SLOW-J photolysis.
!  (bmy, 4/1/03)
!
!  Arguments as Output:
!  ============================================================================
!  (1 ) TOTO3  (REAL*8 ) : Ozone column [DU]
!  (2 ) IDXAIR (INTEGER) : Index for standard temperature profile
!  (3 ) IDXO3  (INTEGER) : Index for standard ozone       profile
!
!  NOTES:
!  (1 ) Updated comments, cosmetic changes.  Now references "grid_mod.f".
!        (bmy, 4/1/03)
!******************************************************************************
!
      IMPLICIT NONE
      
      ! References to F90 modules
      USE GRID_MOD, ONLY : GET_YMID, GET_YOFFSET

#     include "CMN_SIZE"   ! Size parameters
#     include "CMN"        ! 
#     include "comode.h"   ! SMVGEAR II arrays
#     include "comsol.h"   ! SLOW-J arrays

      ! Arguments
      REAL*8,  INTENT(OUT) :: TOTO3(NLAT)      
      INTEGER, INTENT(OUT) :: IDXAIR(JJPAR)
      INTEGER, INTENT(OUT) :: IDXO3(JJPAR)

      ! Local variables
      INTEGER              :: J, J0, JREF, LATI, ISEASON, I, I1, JJ

      ! External functions
      INTEGER, EXTERNAL    :: NUMAIR, NUMO3   

      !=================================================================
      ! GETTOTO3 begins here!
      !=================================================================

      ! Get global offset
      J0 = GET_YOFFSET()

      ! Loop over latitudes
      DO J = 1, JJPAR

         ! Get index of season: 1-Sp, 2-Su, 3-Au, 4-Winter. Sothern  H
         ISEASON = MONTH / 3
         IF ( ISEASON == 0 ) ISEASON=4

          ! Get index of season: 1-Sp, 2-Su, 3-Au, 4-Winter. Northern H
         JREF = J + J0
         IF( JREF < JGLOB/2 ) THEN
            ISEASON = ISEASON - 2
            IF ( ISEASON <= 0 ) ISEASON = ISEASON + 4
         ENDIF
         
         ! GEOS-CHEM grid box latitude
         LATI = GET_YMID( J )

         ! Get indices for standard air density and ozone arrays
         IDXAIR(J) = NUMAIR( ISEASON, ABS(LATI) )
         IDXO3(J)  = NUMO3 ( ISEASON, ABS(LATI) )

         I = LATI / 15. + 6
         IF ( I > 11 ) I = 11  
         I1 = I + 1
         IF ( I  < 1  ) I  = 1
         IF ( I1 > 11 ) I1 = 11
         
         ! Find the weights for latitudes I and I1
         JJ       = LATI - ( I - 6 ) * 15
         TOTO3(J) = ( O3DU(I1,MONTH)*JJ + O3DU(I,MONTH)*(15-JJ) ) / 15.

      ENDDO

      ! Return to calling program
      END SUBROUTINE GETTOTO3
