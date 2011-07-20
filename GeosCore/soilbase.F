! $Id: soilbase.f,v 1.1 2009/09/16 14:06:07 bmy Exp $
      REAL*8 FUNCTION SOILBASE(I,J,M,NN,PULSE)

C**********************************************************************
C                                                                     *
C  HARVARD ATMOSPHERIC CHEMISTRY MODELING GROUP                       *
C  MODULE FOR SOIL NOx EMISSIONS                                      *
C  by Yuhang Wang, Gerry Gardner and Prof. Daniel Jacob               *
C  (Release V2.1)                                                     *
C                                                                     *
C  Contact person: Bob Yantosca (bmy@io.harvard.edu)                  *
C                                                                     *
C**********************************************************************
C Be sure to force double precision with the DBLE function            *
C and the "D" exponent, wherever necessary (bmy, 10/6/99)             *
C**********************************************************************

      IMPLICIT NONE

C**********************************************************************
C SOILBASE  = Emissions                                               *
C SOILAW    = Wet biome coefficient                                   *
C SOILAD    = Dry biome coefficient                                   *
C SOILPREP  = Two month observed precip (mm/day/box                   *
C             (divide by # of days in month))                         *
C NN        = Soil type                                               *
C M         = Index to land box                                       *
C SOILFERT  = Ferterlizers                                            *
C UNITCONV  = Convert from NG_N/(M^2*S) to MOLECULES/CM^2/S           *
C**********************************************************************

#     include "F77_CMN_SIZE"
#     include "commsoil.h"

      INTEGER I,J,M,NN
      REAL*8  PULSE,UNITCONV
      DATA    UNITCONV /4.3D9/   !NG_N/(M^2*S)->MOLECULES/CM^2/S
      
      IF (NN.EQ.1) THEN
C Desert
         SOILBASE=0.D0

      ELSE IF (NN.EQ.2) THEN
C Tropical rain forest
         IF (SOILPREP(2,M).GT.1.D0) THEN
C WET season
            SOILBASE=SOILAW(2)
         ELSE
C DRY season
            SOILBASE=SOILAD(2)
         END IF

      ELSE IF (NN.EQ.8.OR.NN.EQ.9) THEN

         SOILBASE=SOILAW(NN)
         IF (NN.EQ.9) SOILBASE=SOILBASE/30.D0

      ELSE
C Other
         IF (SOILPULS(1,M).GT.0.D0) THEN
C DRY
            SOILBASE=SOILAD(NN)*PULSE
         ELSE
C WET
            SOILBASE=SOILAW(NN)
         END IF
      END IF
C Convert units
      SOILBASE=SOILBASE*UNITCONV

      RETURN
      END
