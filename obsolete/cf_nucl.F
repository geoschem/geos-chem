#if   defined ( TOMAS )
C     **************************************************
C     *  cf_nucl                                     *
C     **************************************************

C     WRITTEN BY Jeff Pierce, April 2007

C     This subroutine calculates the barrierless nucleation rate and radius of the 
C     critical nucleation cluster using the parameterization of...

C     Clement and Ford (1999) Atmos. Environ. 33:489-499

      SUBROUTINE cf_nucl(tempi,rhi,cna,nh3ppt,fn)

      IMPLICIT NONE

C-----INPUTS------------------------------------------------------------

      real tempi                ! temperature of air [K]
      real rhi                  ! relative humidity of air as a fraction
      double precision cna      ! concentration of gas phase sulfuric acid [molec cm-3]
      double precision nh3ppt   ! mixing ratio of ammonia in ppt

C-----OUTPUTS-----------------------------------------------------------

      double precision fn                   ! nucleation rate [cm-3 s-1]
      double precision rnuc                 ! critical cluster radius [nm]

C-----INCLUDE FILES-----------------------------------------------------

C-----ARGUMENT DECLARATIONS---------------------------------------------

C-----VARIABLE DECLARATIONS---------------------------------------------

      double precision temp                 ! temperature of air [K]
      double precision rh                   ! relative humidity of air as a fraction
      double precision alpha1

C-----ADJUSTABLE PARAMETERS---------------------------------------------

C-----CODE--------------------------------------------------------------

      temp=dble(tempi)
      rh=dble(rhi)

      if (nh3ppt .lt. 0.1) then
         alpha1=4.276e-10*sqrt(temp/293.15) ! For sulfuric acid
      else
         alpha1=3.684e-10*sqrt(temp/293.15) ! For ammonium sulfate
      endif
      fn = alpha1*cna**2*3600.
c sensitivity       fn = 1.e-3 * fn ! 10^-3 tuner
      if (fn.gt.1.0e9) fn=1.0e9 ! For numerical conversion

 10   return
      end
#endif
