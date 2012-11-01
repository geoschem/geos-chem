! $Id: fertadd.f,v 1.1 2009/09/16 14:06:29 bmy Exp $
      FUNCTION FERTADD( J, M, NN )
!
!******************************************************************************
!  Subroutine FERTADD computes the amount of soil fertilizer released
!  in a particular grid box according to the Yienger & Levy scheme.
!  (yhw, gmg, djj, 1994; bmy, 2/11/03)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) J  (INTEGER) : Grid box latitude index
!  (2 ) M  (INTEGER) : Grid box surface index (M=1,NLAND)
!  (3 ) NN (INTEGER) : Land type index
!
!  References:
!  ============================================================================
!  (1 ) Yienger, J.J, and H. Levy II, "Empirical model of global soil-biogenic
!        NOx emissions", JGR, 100 (D6), pp. 11447-11464, June 20, 1995.
!
!  NOTES:
!  (1 ) Original code by by Yuhang Wang, Gerry Gardner and Prof. Daniel Jacob
!        written in the early 1990's.  Updated and modified for GEOS-CHEM by
!        Bob Yantosca.  Now uses function GET_YMID of "grid_mod.f" to compute
!        grid box latitudes.  Now use function GET_MONTH from "time_mod.f".
!        Removed reference to header file CMN.  Updated comments, 
!        cosmetic changes. (bmy, 2/11/03)
!  (2 ) Add LANTHRO switch to correctly turn off anthropogenic emissions. 
!        (ccc, 4/15/09)
!  01 Mar 2012 - R. Yantosca - Now use GET_AREA_M2(I,J,L) from grid_mod.F90
!******************************************************************************
!
      ! References to F90 modules
      USE GRID_MOD,     ONLY : GET_YMID
      USE TIME_MOD,     ONLY : GET_MONTH
      USE LOGICAL_MOD,  ONLY : LANTHRO

      USE CMN_SIZE_MOD     ! Size parameters
      USE COMMSOIL_MOD   ! SOILFERT

      IMPLICIT NONE
#     include "define.h"

      ! Arguments
      INTEGER, INTENT(IN) :: J, M, NN

      ! Local variables
      REAL*8              :: Y
      REAL*8, PARAMETER   :: UNITCONV = 4.3d9

      ! Function value
      REAL*8              :: FERTADD

      !=================================================================
      ! FERTADD begins here!
      !=================================================================

      ! Initialize
      FERTADD = 0.D0

      ! Return if soil types are not correct
      ! Soil type 8 refers to different kinds of farmland
      ! Soil type 9 refers to rice paddies 
      IF ( NN /= 8 .and. NN /= 9 ) RETURN

      ! Return if anthropogenic emissions are turned off (ccc, 4/15/09)
      IF (.not.LANTHRO) RETURN 

      ! Latitude of grid box [degrees]
      ! NOTE: For now set I and L to 1.  The soil NOx routines will
      ! be replaced shortly. (bmy, 3/2/12)
      Y = GET_YMID( 1, J, 1 )

      !=================================================================
      ! Case 1: Northern Hemisphere midlatitudes ( Y > 28 degrees )
      !=================================================================
      IF ( Y > 28d0 ) THEN

         ! May, June, July, August...
         IF ( GET_MONTH() >= 5 .and. GET_MONTH() <= 8 ) THEN

            ! NH summer: use value from SOILFERT
            FERTADD = SOILFERT(M)                   

         ELSE

            ! NH winter: no soil NOx emissions
            FERTADD = 0.D0
         ENDIF

      !=================================================================
      ! Case 2: Tropics ( -28 <= Y < 28 degrees )
      !=================================================================
      ELSE IF ( Y > -28d0 ) THEN 

         ! Tropics: use value from soilfert
         FERTADD = SOILFERT(M)

      !=================================================================
      ! Case 3: Southern hemisphere midlatitudes  ( Y <= -28 degrees )
      !=================================================================
      ELSE

         ! Jan, Feb, Nov, Dec
         IF ( GET_MONTH() <= 2 .or. GET_MONTH() >= 11 ) THEN

            ! SH summer: use the values from SOILFERT
            FERTADD = SOILFERT(M)
         ELSE

            ! SH winter: no fertilizer emissions
            FERTADD = 0.d0
         ENDIF
      ENDIF

      !=================================================================
      ! Unit conversion
      !=================================================================

      ! Yienger & Levy state that over rice paddies the fertilizer
      ! emissions should be cut by a factor of 30, since the very
      ! wet soil of rice paddies impedes NOx emission.
      IF ( NN == 9 ) FERTADD = FERTADD / 30.D0

      ! Convert [ng N/m2/s] to [molec/cm2/s]
      FERTADD = FERTADD * UNITCONV
      
      ! Return to calling program
      END FUNCTION FERTADD
