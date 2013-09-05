! $Id: soilbase.f,v 1.1 2009/09/16 14:06:07 bmy Exp $
      REAL*8 FUNCTION SOILBASE(I,J,M,NN,PULSE)

!**********************************************************************
!                                                                     *
!  HARVARD ATMOSPHERIC CHEMISTRY MODELING GROUP                       *
!  MODULE FOR SOIL NOx EMISSIONS                                      *
!  by Yuhang Wang, Gerry Gardner and Prof. Daniel Jacob               *
!  (Release V2.1)                                                     *
!                                                                     *
!  Contact person: Bob Yantosca (bmy@io.harvard.edu)                  *
!                                                                     *
!**********************************************************************
! Be sure to force double precision with the DBLE function            *
! and the "D" exponent, wherever necessary (bmy, 10/6/99)             *
!**********************************************************************

      USE CMN_SIZE_MOD
      USE COMMSOIL_MOD

      IMPLICIT NONE
#     include "define.h"

!**********************************************************************
! SOILBASE  = Emissions                                               *
! SOILAW    = Wet biome coefficient                                   *
! SOILAD    = Dry biome coefficient                                   *
! SOILPREP  = Two month observed precip (mm/day/box                   *
!             (divide by # of days in month))                         *
! NN        = Soil type                                               *
! M         = Index to land box                                       *
! SOILFERT  = Ferterlizers                                            *
! UNITCONV  = Convert from NG_N/(M^2*S) to MOLECULES/CM^2/S           *
!**********************************************************************


      INTEGER I,J,M,NN
      REAL*8  PULSE,UNITCONV
      DATA    UNITCONV /4.3D9/   !NG_N/(M^2*S)->MOLECULES/CM^2/S
      
      IF (NN.EQ.1) THEN
! Desert
         SOILBASE=0.D0

      ELSE IF (NN.EQ.2) THEN
! Tropical rain forest
         IF (SOILPREP(2,M).GT.1.D0) THEN
! WET season
            SOILBASE=SOILAW(2)
         ELSE
! DRY season
            SOILBASE=SOILAD(2)
         END IF

      ELSE IF (NN.EQ.8.OR.NN.EQ.9) THEN

         SOILBASE=SOILAW(NN)
         IF (NN.EQ.9) SOILBASE=SOILBASE/30.D0

      ELSE
! Other
         IF (SOILPULS(1,M).GT.0.D0) THEN
! DRY
            SOILBASE=SOILAD(NN)*PULSE
         ELSE
! WET
            SOILBASE=SOILAW(NN)
         END IF
      END IF
! Convert units
      SOILBASE=SOILBASE*UNITCONV

      RETURN
      END
