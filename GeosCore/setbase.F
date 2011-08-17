! $Id: setbase.f,v 1.1 2009/09/16 14:06:10 bmy Exp $
      SUBROUTINE SETBASE( CONVERT, GMONOT )
!
!******************************************************************************
!  Subroutine SETBASE computes the baseline emissions for 
!  ISOPRENE, MONOTERPENES, GRASSLAND ISOPRENE, and METHYL BUTENOL.
!  (bdf, bmy, 8/1/01, 2/11/03)
!
!  Baseline emissions are stored in arrays (from F77_CMN_ISOP and F77_CMN_MONOT)
!  BASEISOP, BASEMONOT, BASEGRASS, BASEMB.  Units are [kg C/box/step].
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) CONVERT (REAL*8) : ISOP  emissions by landtype [atoms C/cm2 leaf/s]
!  (2 ) GMONOT  (REAL*8) : MONOT emissions by landtype [atoms C/cm2 leaf/s]
!
!  NOTES:
!  (1 ) Now use F90 syntax.  Updated comments, cosmetic changes.  Moved 
!        everything to within one I-J loop.  Also removed reference to
!        F77_CMN_O3, which is no longer needed. (bdf, bmy, 8/1/01)
!  (2 ) Replaced all instances of IM with IIPAR and JM with JJPAR, in order
!        to prevent namespace confusion for the new TPCORE (bmy, 6/25/02)
!  (3 ) Now replace DXYP(JREF)*1d4 with routine GET_AREA_CM2 from "grid_mod.f".
!        Now use function GET_TS_EMIS from "grid_mod.f". (bmy, 2/11/03)
!******************************************************************************
!
      ! References to F90 modules
      USE GRID_MOD, ONLY : GET_AREA_CM2
      USE TIME_MOD, ONLY : GET_TS_EMIS

      IMPLICIT NONE

#     include "F77_CMN_SIZE"  ! Size parameters 
#     include "CMN"       ! NSRCE
#     include "F77_CMN_VEL"   ! IJREG, IJUSE, IJLAND
#     include "F77_CMN_ISOP"  ! BASEISOP, BASEGRASS, BASEMB
#     include "F77_CMN_MONOT" ! BASEMONOT

      ! Arguments
      REAL*8, INTENT(IN) :: CONVERT(NVEGTYPE), GMONOT(NVEGTYPE)

      ! Local variables
      INTEGER            :: I, J, IJLOOP, K
      REAL*8             :: DTSRCE, FACTOR

      ! Avogadro's Number
      REAL*8, PARAMETER  :: AVO = 6.023D+23

      !=================================================================
      ! SETBASE begins here!
      !=================================================================

      ! Emission timestep [s]
      DTSRCE = GET_TS_EMIS() * 60d0

      !=================================================================
      ! Set up BASEISOP -- baseline ISOPRENE emissions
      ! Now hardwire molecular weight for Carbon = 0.012 kg/mol
      ! ISOPRENE is traced in terms of equivalent C atoms
      !=================================================================
      IJLOOP = 0

      ! Loop over latitudes
      DO J = 1, JJPAR

         ! Conversion factor from [atoms C/cm2/s] to [kg C/box/step]
         FACTOR = 12d-3 * DTSRCE * GET_AREA_CM2( J ) / AVO

         ! Loop over longitudes
         DO I = 1, IIPAR

            ! 1-D grid box index corresponding to (I,J)
            IJLOOP = IJLOOP + 1

            ! Loop over landtypes per (I,J) box
            DO K = 1, IJREG(IJLOOP)

               ! Baseline emissions for ISOPRENE in [kg C/box/step]
               ! IJLAND+1 is the Olson land type index
               BASEISOP(IJLOOP,K) = CONVERT(IJLAND(IJLOOP,K)+1) * FACTOR

               ! Baseline emissions for MONOTERPENES in [kg C/box/step]
               ! IJLAND+1 is the Olson land type index               
               BASEMONOT(IJLOOP,K) = GMONOT(IJLAND(IJLOOP,K)+1) * FACTOR
            ENDDO

            ! Baseline emissions for GRASSLAND ISOPRENE in [kg C/box/step] 
            ! needed for acetone chemistry.  Based on Kirstine et al 1998.
            BASEGRASS(IJLOOP) = 7.25d10 * FACTOR

            ! Baseline emissions for METHYL BUTENOL in [kg C/box/step] 
            ! needed for acetone chemistry.  Based on  3.2 TgC MB 
            ! emissions in N.america from Guenther 2000
            BASEMB(IJLOOP) = 4.37d11 * FACTOR
         ENDDO
      ENDDO

      ! Return to calling program
      END SUBROUTINE SETBASE
