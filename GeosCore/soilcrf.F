! $Id: soilcrf.f,v 1.1 2009/09/16 14:06:06 bmy Exp $
      REAL*8 FUNCTION SOILCRF(I,J,IREF,JREF,IJLOOP,M,NN,K,
     &                        WINDSQR,SUNCOS)

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
C Updated comments (bmy, 1/24/03)                                     *
C**********************************************************************

      IMPLICIT NONE

C**********************************************************************
C SOILEXC   = Canopy wind extinction coeff.                           *
C WINDSQR   = Wind speed squared                                      *
C XLAI      = LAI of land type element K                              *
C CANOPYNOX = Deposition rate constant for NOx                        *
C NN        = Soil type                                               *
C K         = Number in vegationtype of the grid                      *
C VFNEW     = Ventilation rate constant for NOx                       *
C SOILCRF   = Canopy reduction factor                                 *
C SUNCOS    = Array of cosine( Solar zenith angle ) for grid boxes    *
C**********************************************************************
C                                                                     *
C Wang et al.: [1998] JGR vol. 103 p10713-10725                       *
C                                                                     *
C**********************************************************************

#     include "F77_CMN_SIZE"    ! Size parameters
#     include "F77_CMN_DEP"     ! CANOPYNOX
#     include "commsoil.h"  ! Soil pulsing & wetness variables

      INTEGER I,J,IREF,JREF,M,NN,K,IJLOOP

      REAL*8  WINDSQR,VFDAY,VFNIGHT,VFNEW,SUNCOS(MAXIJ)
      
C**********************************************************************
C coefficient ALPHA (2.8E-2, 5.6E-3) day, night canopy ventilation    *
C time of 1 hour day, 5 hour night                                    *
C VFDAY,VFNIGHT - alpha scaled                                        *
C**********************************************************************

      DATA VFDAY,VFNIGHT /1.0D-2,0.2D-2/ !VENTILATION VEL. IN DAY&NIGHT M/S
        
C For GEOS-CTM, RADIAT is a 3-hour average field.  Replace the test for 
C RADIAT > 0 with a test for SUNCOS > 0.  SUNCOS is the cosine of the
C solar zenith angle, so SUNCOS > 0 is day and SUNCOS < 0 is night.
C In the GEOS model, SUNCOS is is computed every dynamic timestep 
C (15 or 30 mins), and thus is a better indicator of where the
C day-night terminator falls. (bmy, 10/20/99)
C      IF (RADIAT(IJLOOP).GT.0D0) THEN
      IF ( SUNCOS(IJLOOP) .GT. 0D0 ) THEN
         ! Day
         VFNEW=VFDAY
      ELSE 
         ! Night
         VFNEW=VFNIGHT
      END IF

      IF ((XLAI(IREF,JREF,K).GT.0.D0).AND.
     &    (CANOPYNOX(IJLOOP,K).GT.0.D0))THEN

         VFNEW=VFNEW*SQRT(WINDSQR/9.D0*7.D0/XLAI(IREF,JREF,K))*
     *        (SOILEXC(2)/SOILEXC(NN))
         SOILCRF=CANOPYNOX(IJLOOP,K)/(CANOPYNOX(IJLOOP,K)
     *        +VFNEW)
      ELSE
     
         SOILCRF=0.D0
      END IF
      
      ! Return to calling program
      END FUNCTION SOILCRF
