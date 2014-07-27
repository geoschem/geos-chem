!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: drydep_toolbox_mod.F90
!
! !DESCRIPTION: Module DryEep\_ToolBox\_Mod contains routines used for dry
! deposition (and soil NOx emissions) calculations, as implemented into!
! the GEOS-Chem model.
!\\
!\\
! !INTERFACE: 
!
MODULE DryDep_ToolBox_Mod
!
! !USES:
!
  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC  :: BioFit
!
! !PRIVATE MEMBER FUNCTIONS:
!
  PRIVATE :: SUNPARAM 
!
! !REVISION HISTORY:
!  14 Nov 2013 - C. Keller   - Created from BIOFIT.F and SUNPARAM.F 
!  09 Jul 2014 - R. Yantosca - Now use F90 free-format indentation
!  09 Jul 2014 - R. Yantosca - Cosmetic changes in ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: BioFit
!
! !DESCRIPTION: Function BioFit computes the light correction used in the
!  dry deposition and canopy NOx modules.
!\\
!\\
! !INTERFACE:
!
  REAL*8 FUNCTION BIOFIT( COEFF1, XLAI1, SUNCOS1, CFRAC1, NPOLY )
!
! !INPUT PARAMETERS: 
!
    REAL*8,  INTENT(IN) :: COEFF1(NPOLY)   ! Baldocchi drydep coefficients
    REAL*8,  INTENT(IN) :: XLAI1           ! Leaf area index [cm2/cm2]
    REAL*8,  INTENT(IN) :: SUNCOS1         ! Cosine( Solar Zenith Angle )
    REAL*8,  INTENT(IN) :: CFRAC1          ! Cloud fraction [unitless]
    INTEGER, INTENT(IN) :: NPOLY           ! # of drydep coefficients
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

  END FUNCTION BioFit
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: SunParam
!
! !DESCRIPTION: Subroutine SUNPARAM is called by BIOFIT to perform the 
!  light correction used in the dry deposition and canopy NOx modules.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE SUNPARAM( X )
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

  END SUBROUTINE SunParam
!EOC
END MODULE DryDep_ToolBox_Mod
