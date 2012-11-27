! $Id: soilcrf.f,v 1.1 2009/09/16 14:06:06 bmy Exp $
      REAL*8 FUNCTION SOILCRF(I,J,IREF,JREF,IJLOOP,M,NN,K,
     &                        WINDSQR,SUNCOS)

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
! Updated comments (bmy, 1/24/03)                                     *
!**********************************************************************

      USE CMN_SIZE_MOD    ! Size parameters
      USE CMN_DEP_MOD     ! CANOPYNOX
      USE COMMSOIL_MOD  ! Soil pulsing & wetness variables

      IMPLICIT NONE
#     include "define.h"

!**********************************************************************
! SOILEXC   = Canopy wind extinction coeff.                           *
! WINDSQR   = Wind speed squared                                      *
! XLAI      = LAI of land type element K                              *
! CANOPYNOX = Deposition rate constant for NOx                        *
! NN        = Soil type                                               *
! K         = Number in vegationtype of the grid                      *
! VFNEW     = Ventilation rate constant for NOx                       *
! SOILCRF   = Canopy reduction factor                                 *
! SUNCOS    = Array of cosine( Solar zenith angle ) for grid boxes    *
!**********************************************************************
!                                                                     *
! Wang et al.: [1998] JGR vol. 103 p10713-10725                       *
!                                                                     *
!**********************************************************************


      INTEGER I,J,IREF,JREF,M,NN,K,IJLOOP

      REAL*8  WINDSQR,VFDAY,VFNIGHT,VFNEW,SUNCOS(MAXIJ)
      
!**********************************************************************
! coefficient ALPHA (2.8E-2, 5.6E-3) day, night canopy ventilation    *
! time of 1 hour day, 5 hour night                                    *
! VFDAY,VFNIGHT - alpha scaled                                        *
!**********************************************************************

      DATA VFDAY,VFNIGHT /1.0D-2,0.2D-2/ !VENTILATION VEL. IN DAY&NIGHT M/S
        
! For GEOS-CTM, RADIAT is a 3-hour average field.  Replace the test for 
! RADIAT > 0 with a test for SUNCOS > 0.  SUNCOS is the cosine of the
! solar zenith angle, so SUNCOS > 0 is day and SUNCOS < 0 is night.
! In the GEOS model, SUNCOS is is computed every dynamic timestep 
! (15 or 30 mins), and thus is a better indicator of where the
! day-night terminator falls. (bmy, 10/20/99)
!      IF (RADIAT(IJLOOP).GT.0D0) THEN
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
