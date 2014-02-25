!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: biofit
!
! !DESCRIPTION: Function BIOFIT computes the light correction used in the
!  dry deposition and canopy NOx modules.
!\\
!\\
! !INTERFACE:
!
      REAL*8 FUNCTION BIOFIT( COEFF1, XLAI1, SUNCOS1, CFRAC1 )
!
! !USES:
!
      USE CMN_SIZE_MOD, ONLY : NPOLY        ! # of drydep coefficients

      IMPLICIT NONE
!
! !INPUT PARAMETERS: 
!
      REAL*8, INTENT(IN) :: COEFF1(NPOLY)   ! Baldocchi drydep coefficients
      REAL*8, INTENT(IN) :: XLAI1           ! Leaf area index [cm2/cm2]
      REAL*8, INTENT(IN) :: SUNCOS1         ! Cosine( Solar Zenith Angle )
      REAL*8, INTENT(IN) :: CFRAC1          ! Cloud fraction [unitless]
!
! !REMARKS:
!  This routine is ancient code from Yuhang Wang.  It was part of the old
!  Harvard-GISS CTM and was ported into GEOS-Chem.  See this reference for 
!  more information: 
!                                                                             .
!    Wang, Y., D.J. Jacob, and J.A. Logan, "Global simulation of tropospheric 
!     O3-NOx-hydrocarbon chemistry, 1. Model formulation", J. Geophys. Res., 
!     103/D9, 10,713-10,726, 1998. 
! 
! !REVISION HISTORY: 
!  13 Dec 2012 - R. Yantosca - Added ProTeX headers
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
      INTEGER, PARAMETER :: KK = 4 
!
! !LOCAL VARIABLES:
!
      REAL*8             :: TERM(KK)
      REAL*8             :: REALTERM(NPOLY)
      INTEGER            :: K,K1,K2,K3

      !=================================================================
      ! BIOFIT begins here!
      !=================================================================
      TERM(1)=1.
      TERM(2)=XLAI1
      TERM(3)=SUNCOS1
      TERM(4)=CFRAC1
      CALL SUNPARAM(TERM(2))
      K=0
      DO K3=1,KK
        DO K2=K3,KK
          DO K1=K2,KK
            K=K+1
            REALTERM(K)=TERM(K1)*TERM(K2)*TERM(K3)
          END DO
        END DO
      END DO
      BIOFIT=0
      DO K=1,NPOLY
        BIOFIT=BIOFIT+COEFF1(K)*REALTERM(K)
      END DO
      IF (BIOFIT.LT.0.1) BIOFIT=0.1

      END FUNCTION BIOFIT
!EOC
