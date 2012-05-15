! $Id: arsl1k.f,v 1.3 2010/02/02 16:57:55 bmy Exp $
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: arsl1k
!
! !DESCRIPTION: Function ARSL1K calculates the 1st-order loss rate of species 
!  on wet aerosol surface.
!\\
!\\
! !INTERFACE:
!
      REAL*8 FUNCTION ARSL1K( AREA, RADIUS, DENAIR, STKCF, STK, SQM )
!
! !USES:
!
      IMPLICIT NONE
#     include "define.h"
!
! !INPUT PARAMETERS: 
!
      ! Surface  area of wet aerosols/volume of air [cm2/cm3]
      REAL*8, INTENT(IN) :: AREA     

      ! Radius of wet aerosol [cm], order of 0.01-10 um;
      ! Note that radius here is Rd, not Ro
      REAL*8, INTENT(IN) :: RADIUS 
  
      ! Density of air [#/cm3]
      REAL*8, INTENT(IN) :: DENAIR  
 
      ! Sticking coefficient [unitless], order of 0.1
      REAL*8, INTENT(IN) :: STKCF  
  
      ! Square root of temperature [K]
      REAL*8, INTENT(IN) :: STK  
    
      ! Square root of molecular weight [g/mole]
      REAL*8, INTENT(IN) :: SQM      
!
! !REMARKS:
!  The 1st-order loss rate on wet aerosol (Dentener's Thesis, p. 14)
!  is computed as:
!                                                                             .
!      ARSL1K [1/s] = area / [ radius/dfkg + 4./(stkcf * xmms) ]        
!                                                                             .
!  where XMMS = Mean molecular speed [cm/s] = sqrt(8R*TK/pi/M) for Maxwell 
!        DFKG = Gas phase diffusion coeff [cm2/s] (order of 0.1)

! !REVISION HISTORY:
!  01 Jul 1994 - lwh, jyl, gmg, djj - Initial version 
!  04 Apr 2003 - R. Yantosca - Updated comments, cosmetic changes
!  07 Apr 2004 - R. Yantosca - Now return w/ default value if RADIUS is zero 
!                              (i.e. is smaller than a very small number)
!  03 Dec 2009 - R. Yantosca - Prevent div-by-zero errors by returning the
!                              default value if any of the args are zero 
!  03 Dec 2009 - R. Yantosca - Added ProTeX Header
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL*8 :: DFKG

      !=================================================================
      ! ARSL1K begins here!
      !=================================================================
      IF ( AREA < 0d0   .or. DENAIR < 1d-30 .or. RADIUS < 1d-30  .or.
     &     SQM  < 1d-30 .or. STK    < 1d-30 .or. STKCF  < 1d-30 ) THEN

         ! Use default value if any of the above values are zero
         ! This will prevent div-by-zero errors in the eqns below
         ARSL1K = 1.D-30

      ELSE

         ! DFKG = Gas phase diffusion coeff [cm2/s] (order of 0.1)
         DFKG  = 9.45D17/DENAIR * STK * SQRT(3.472D-2 + 1.D0/(SQM*SQM))

         ! Compute ARSL1K according to the formula listed above
         ARSL1K = AREA / ( RADIUS/DFKG + 2.749064E-4*SQM/(STKCF*STK) )

      ENDIF

      ! Return to calling program
      END FUNCTION ARSL1K
!EOC
