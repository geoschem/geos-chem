!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: sunparam
!
! !DESCRIPTION: Subroutine SUNPARAM is called by BIOFIT to perform the 
!  light correction used in the dry deposition and canopy NOx modules.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE SUNPARAM( X )
!
! !USES:
!
      IMPLICIT NONE
!
! !DEFINED PARAMETERS:
!
      INTEGER, PARAMETER    :: NN = 3  ! # of variables (LAI, SUNCOS, CLDFRC)
!
! !INPUT/OUTPUT PARAMETERS:
!
      REAL*8, INTENT(INOUT) :: X(NN)   ! LAI, SUNCOS, or cloud fraction
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

      !===============================================
      ! the sequence is lai,suncos,cloud fraction
      !===============================================

      !  ND = scaling factor for each variable
      INTEGER ND(NN),I
      DATA ND /55,20,11/

      !  X0 = maximum for each variable
      REAL*8 X0(NN),XLOW
      DATA X0 /11.,1.,1./

      DO I=1,NN
         X(I)=MIN(X(I),X0(I))
         ! XLOW = minimum for each variable
         IF (I.NE.3) THEN
            XLOW=X0(I)/REAL(ND(I))
         ELSE
            XLOW= 0.
         END IF
         X(I)=MAX(X(I),XLOW)
         X(I)=X(I)/X0(I)
      END DO

      END SUBROUTINE SUNPARAM
!EOC
