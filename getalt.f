! $Id: getalt.f,v 1.1 2003/06/30 20:26:01 bmy Exp $
      SUBROUTINE GETALT( ALT, SURFALT, ROVMG )
!
!******************************************************************************
!  Subroutine GETALT uses the barometric law to determine the altitudes
!  of GEOS-CHEM model layers. (bmy, 4/1/03)
! 
!  Arguments as Output:
!  ============================================================================
!  (1 ) ALT     (INTEGER) : Model layer midpoint altitudes [cm]
!  (2 ) SURFALT (REAL*8 ) : Altitude of surface grid boxes [m]
!  (3 ) ROVMG   (REAL*8 ) : Gas constant ??
!
!  NOTES:
!  (1 ) Updated comments, cosmetic changes (bmy, 4/1/03)
!******************************************************************************
!
      ! References to F90 modules
      USE DAO_MOD,      ONLY : PHIS,      T 
      USE GRID_MOD,     ONLY : GET_YMID,  GET_XOFFSET, GET_YOFFSET
      USE PRESSURE_MOD, ONLY : GET_PEDGE, GET_PCENTER

      IMPLICIT NONE

#     include "CMN_SIZE"  ! Size parameters
#     include "CMN"       ! ???
#     include "comode.h"  ! SMVGEAR II arrays

      ! Arguments
      REAL*8, INTENT(OUT) :: ALT(MAXIJ,NPVERT)
      REAL*8, INTENT(OUT) :: SURFALT(MAXIJ)
      REAL*8, INTENT(OUT) :: ROVMG
      
      ! Local variables
      INTEGER IJLOOP, LATI, I, J, L
      REAL*8  TOTAL,  DSUM ,TTEMP
      
      !=================================================================
      ! GETALT begins here!
      !=================================================================

      ! Gas constant [units ??]
      ROVMG  = 8.3144d0 * 1000.d0 / ( 28.966d0 * 9.81d0 )

      ! 1-D grid box index
      IJLOOP = 0

      ! Loop over latitudes
      DO J = 1, NLAT
         
         ! Grid box latitude [degrees]
         LATI = GET_YMID( J )

         ! Loop over longitudes
         DO I = 1, NLONG

            ! 1-D grid box index
            IJLOOP = IJLOOP + 1

            ! Surface altitude in [m]
            SURFALT(IJLOOP) = PHIS(I,J)

            ! Initialize
            TOTAL = 0.D0

            ! Use barometric law to determine altitudes of CTM model layers.
            ! d(ln P) = -dz/H
            DO L = 1, NVERT
               IF( L == 1 ) THEN
                  DSUM  = GET_PEDGE(I,J,1) / GET_PCENTER(I,J,1)
                  TTEMP = T(I,J,L)
               ELSE
                  DSUM  = GET_PCENTER(I,J,L-1) / GET_PCENTER(I,J,L)
                  TTEMP = 0.5d0 * ( T(I,J,L-1) + T(I,J,L) )
               ENDIF

               ! lwh  Fixed bug in calculating altitude. Should use log, 
               ! not log10.  Multiply by 100 to convert from meters to cm.
               DSUM  = LOG( DSUM ) * TTEMP * ROVMG * 100.
               TOTAL = TOTAL + DSUM

               ! ALT() in cm here, the height of the midpoint(P) of layer L
               ALT(IJLOOP,L) = TOTAL
            ENDDO
         ENDDO
      ENDDO

      ! Return to calling program
      END SUBROUTINE GETALT
