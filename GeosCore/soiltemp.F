C $Id: soiltemp.f,v 1.1 2009/09/16 14:06:05 bmy Exp $
      REAL*8 FUNCTION SOILTEMP(I,J,M,NN,TMMP0)

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
C Yienger and Levy [1995] JGR 100, 11447-11464                        *
C**********************************************************************
C NN        = Soil type                                               *
C SOILTEMP  = Temperature factor                                      *
C TMMP0     = Local air temperature (C),                              *
C             include diurnal temp variation                          *
C SOILTA    = Coefficient used to convert from surface temperture to  *
C             soil temperature                                        *
C SOILTB    = Coefficient used to convert from surface temperture to  *
C             soil temperature                                        *
C**********************************************************************

#     include "F77_CMN_SIZE"
#     include "commsoil.h"

      INTEGER I,J,M,NN
      REAL*8  TMMP0,TMMP

      TMMP=TMMP0
C DRY
C SURFACE TEMPERATURE->SOIL TEMPERATURE
C Convert the lowest model level air temperature to soil temperature
C based on observations of Johansson et. al. [1988]
C add 5 degrees C to model temperature
C
      IF (NN.LE.2) THEN
C Desert and rain forest
         SOILTEMP=1.D0
C                                      Agric.      Rice paddies
      ELSE IF (SOILPULS(1,M).GT.0..AND.NN.NE.8.AND.NN.NE.9) THEN
C DRY
         TMMP=TMMP+5.D0
         IF (TMMP.GT.30.D0) THEN
C Optimal
	    SOILTEMP=1.D0
         ELSE IF (TMMP.GT.0.D0) THEN
C Cold-linear
	    SOILTEMP=TMMP/30.D0
         ELSE
	    SOILTEMP=0.D0
         END IF
      ELSE
C WET

C SURFACE TEMPERATURE->SOIL TEMPERATURE
C**********************************************************************
C Convert the lowest model level air temperature to soil temperature  *
C Use the empirical relationships derived by Williams et al. [1992b]  *
C ECO SYSTEM DEPENDENT                                                *
C**********************************************************************

         TMMP=SOILTA(NN)*TMMP+SOILTB(NN)
         IF (TMMP.GE.30.D0) THEN
C Optimal
	    SOILTEMP=21.97D0
         ELSE IF (TMMP.GE.10.D0) THEN
C Exponential
            SOILTEMP=EXP(0.103D0*TMMP)
         ELSE IF (TMMP.GT.0.D0) THEN
C Cold-linear
	    SOILTEMP=0.28D0*TMMP
         ELSE
	    SOILTEMP=0.D0
         END IF
      END IF

      RETURN
      END
