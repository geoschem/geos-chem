! $Id: soiltemp.f,v 1.1 2009/09/16 14:06:05 bmy Exp $
      REAL*8 FUNCTION SOILTEMP(I,J,M,NN,TMMP0)

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
! Yienger and Levy [1995] JGR 100, 11447-11464                        *
!**********************************************************************
! NN        = Soil type                                               *
! SOILTEMP  = Temperature factor                                      *
! TMMP0     = Local air temperature (C),                              *
!             include diurnal temp variation                          *
! SOILTA    = Coefficient used to convert from surface temperture to  *
!             soil temperature                                        *
! SOILTB    = Coefficient used to convert from surface temperture to  *
!             soil temperature                                        *
!**********************************************************************

      INTEGER I,J,M,NN
      REAL*8  TMMP0,TMMP

      TMMP=TMMP0
! DRY
! SURFACE TEMPERATURE->SOIL TEMPERATURE
! Convert the lowest model level air temperature to soil temperature
! based on observations of Johansson et. al. [1988]
! add 5 degrees C to model temperature
!
      IF (NN.LE.2) THEN
! Desert and rain forest
         SOILTEMP=1.D0
!                                      Agric.      Rice paddies
      ELSE IF (SOILPULS(1,M).GT.0..AND.NN.NE.8.AND.NN.NE.9) THEN
! DRY
         TMMP=TMMP+5.D0
         IF (TMMP.GT.30.D0) THEN
! Optimal
	    SOILTEMP=1.D0
         ELSE IF (TMMP.GT.0.D0) THEN
! Cold-linear
	    SOILTEMP=TMMP/30.D0
         ELSE
	    SOILTEMP=0.D0
         END IF
      ELSE
! WET

! SURFACE TEMPERATURE->SOIL TEMPERATURE
!**********************************************************************
! Convert the lowest model level air temperature to soil temperature  *
! Use the empirical relationships derived by Williams et al. [1992b]  *
! ECO SYSTEM DEPENDENT                                                *
!**********************************************************************

         TMMP=SOILTA(NN)*TMMP+SOILTB(NN)
         IF (TMMP.GE.30.D0) THEN
! Optimal
	    SOILTEMP=21.97D0
         ELSE IF (TMMP.GE.10.D0) THEN
! Exponential
            SOILTEMP=EXP(0.103D0*TMMP)
         ELSE IF (TMMP.GT.0.D0) THEN
! Cold-linear
	    SOILTEMP=0.28D0*TMMP
         ELSE
	    SOILTEMP=0.D0
         END IF
      END IF

      RETURN
      END
