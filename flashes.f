C $Id: flashes.f,v 1.1 2003/06/30 20:26:08 bmy Exp $
      SUBROUTINE FLASHES( I, J, HEIGHT, FLASHRATE, XIGCRATIO )
!
!*****************************************************************************
!  Subroutine FLASHES determines the rate of lightning flashes per minute,
!  based on the height of convective cloud tops, and the inter-cloud
!  to cloud-ground strike ratio.  FLASHES has been optimized for the GEOS-CTM.
! 
!  Rewritten by Bob Yantosca for Harvard Atmospheric Sciences 
!  (10/9/97, 6/26/00)
!
!  Arguments as Input:
!  =====================================================================
!  (1-2) I, J      : Longitude/latitude indices
!  (3  ) HEIGHT    : Height of convective cloud tops in meters
!  (4  ) LWI       : DAO Land/water indices (via CMN_LWI)
!
!  Arguments as Output:
!  =====================================================================
!  (4  ) FLASHRATE : Lightning flash rate in flashes/minute
!  (5  ) XICGRATIO : Intercloud (IC) flashes / Cloud-Ground (CG) flashes
!
!  References:
!  =====================================================================
!  (1) Price & Rind (1992), JGR, vol. 97, 9919-9933.
!  (2) Price & Rind (1994), M. Weather Rev, vol. 122, 1930-1939.
!
!  NOTES:
!  (1) FLASHES is written in Fixed-Form Fortran 90.  Also use F90 
!      declaration syntax. (bmy, 6/26/00)
!
!  (2) Eliminate obsolete code from 6/26/00 (bmy, 8/31/00)
!*****************************************************************************
!    
      ! References to F90 modules
      USE DAO_MOD, ONLY : IS_LAND, IS_WATER

      IMPLICIT NONE

#     include "CMN_SIZE"
      
      ! Arguments
      INTEGER, INTENT(IN)  :: I, J
      REAL*8,  INTENT(IN)  :: HEIGHT
      REAL*8,  INTENT(OUT) :: FLASHRATE, XIGCRATIO
     
      ! Local variables
      INTEGER M
      REAL*8  X, FLAND, FOCEAN, FRATIO
!
!*****************************************************************************
!  FLAND  is the lightning flash/minute rate over land
!  FOCEAN is the lightning flash/minute rate over water
!  X      is the convective cloud top height in km.
!
!  These functions are taken from Price & Rind (1992).
!*****************************************************************************
!      
      FLAND (X) = (3.44D-5) * (X**4.9 )   
      FOCEAN(X) = (6.4 D-4) * (X**1.73)  
!
!*****************************************************************************
!  FRATIO is the IC/CG flash ratio as a function of 
!  X      is the lightning flash rate per minute.
!
!  This function was taken from Price & Rind (1992).
!*****************************************************************************
!
      FRATIO(X) = (2.7) * (X**0.5)      
!
!*****************************************************************************
!  Use DAO land-water indices to determine if over land or water, in
!  order to use the appropriate lightning flash/minute function.
!
!  Note...the conversion from m/km is accounted for below.
!
!  LWI(I,J) = 1 is a water point,    LWI(I,J) = 2 is a land point,
!  LWI(I,J) = 3 is a land ice point, LWI(I,J) = 4 is a sea ice point 
!*****************************************************************************
!
      ! Now use functions IS_LAND and IS_WATER to determine if the box (I,J) 
      ! is a land or water point, for all data sets (bmy, 6/26/00)
      IF ( IS_LAND(I,J)  ) THEN
         FLASHRATE = FLAND( HEIGHT * 1.0D-3 )

      ELSE IF ( IS_WATER(I,J) ) THEN 
         FLASHRATE = FOCEAN( HEIGHT * 1.0D-3 )

      ENDIF
!
!*****************************************************************************
!  Compute the IC/CG ratio as a function of flash rate
!*****************************************************************************
!
      XIGCRATIO = FRATIO(FLASHRATE)

      END SUBROUTINE FLASHES







