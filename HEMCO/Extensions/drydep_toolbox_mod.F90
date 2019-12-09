!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
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
  USE HCO_ERROR_MOD, ONLY : hp     ! Precision (hp)

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC  :: BIOFIT

  INTERFACE BIOFIT
    MODULE PROCEDURE BIOFIT_R4
    MODULE PROCEDURE BIOFIT_R8
  END INTERFACE

!
! !PRIVATE MEMBER FUNCTIONS:
!
  PRIVATE :: SUNPARAM_R4
  PRIVATE :: SUNPARAM_R8

!
! !REVISION HISTORY:
!  14 Nov 2013 - C. Keller   - Created from BIOFIT.F and SUNPARAM.F
!  09 Jul 2014 - R. Yantosca - Now use F90 free-format indentation
!  09 Jul 2014 - R. Yantosca - Cosmetic changes in ProTeX headers
!  11 Dec 2014 - M. Yannetti - Split BIOFIT into R4 and R8 versions
!                              Split SUNPARAM into R4 and R8 versions
!EOP
!------------------------------------------------------------------------------
!BOC
CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: BIOFIT_R4
!
! !DESCRIPTION: Function BioFit computes the light correction used in the
!  dry deposition and canopy NOx modules.
!\\
!\\
! !INTERFACE:
!
  FUNCTION BIOFIT_R4( COEFF1, XLAI1, SUNCOS1, CFRAC1, NPOLY ) RESULT ( BIO_FIT )
!
! !INPUT PARAMETERS:
!
    REAL*4,  INTENT(IN) :: COEFF1(NPOLY)   ! Baldocchi drydep coefficients
    REAL*4,  INTENT(IN) :: XLAI1           ! Leaf area index [cm2/cm2]
    REAL*4,  INTENT(IN) :: SUNCOS1         ! Cosine( Solar Zenith Angle )
    REAL*4,  INTENT(IN) :: CFRAC1          ! Cloud fraction [unitless]
    INTEGER, INTENT(IN) :: NPOLY           ! # of drydep coefficients
!
! !RETURN VALUE:
!
    REAL*4              :: BIO_FIT         ! Resultant light correction
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
!  09 Dec 2014 - R. Yantosca - Now use BIO_FIT as the return value
!  11 Dec 2014 - M. Yannetti - Split from BIO_FIT
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
    REAL*4             :: TERM(KK)
    REAL*4             :: REALTERM(NPOLY)
    INTEGER            :: K,K1,K2,K3

    !=================================================================
    ! BIOFIT begins here!
    !=================================================================
    TERM(1)=1.
    TERM(2)=XLAI1
    TERM(3)=SUNCOS1
    TERM(4)=CFRAC1
    CALL SUNPARAM_R4(TERM(2))
    K=0
    DO K3=1,KK
       DO K2=K3,KK
          DO K1=K2,KK
             K=K+1
             REALTERM(K)=TERM(K1)*TERM(K2)*TERM(K3)
          END DO
       END DO
    END DO

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%% KLUDGE: UNCOMMENT THIS CODE TO GET IDENTICAL RESULTS w/ v10-01e!
!%%%%% (bmy, myannetti, 12/10/14)
!%%%%%
!%%%%% NOTE: This code did not use Fortran "d" exponents, so it was not
!%%%%% forcing values to REAL*8 precision.  This was not as precise as it
!%%%%% could have been, but this code was like this for many years.
!%%%%%    BIO_FIT=0
!%%%%%    DO K=1,NPOLY
!%%%%%       BIO_FIT=BIO_FIT+COEFF1(K)*REALTERM(K)
!%%%%%    END DO
!%%%%%    IF (BIO_FIT .LT. 0.1 ) BIO_FIT=0.1
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! Now explicitly use REAL*8 precision.  This will cause very small
    ! differences at the level of numerical noise when comparing to
    ! prior states of the code like v10-01e.  But this is something that
    ! we can live with.  Stick with REAL*8 precision for now, but we'll
    ! try to implement flexible precision into this routine at a later
    ! point. (bmy, myannetti, 12/10/14)
    BIO_FIT = 0e0
    DO K = 1, NPOLY
       BIO_FIT = BIO_FIT + COEFF1(K)*REALTERM(K)
    END DO
    IF ( BIO_FIT .LT. 0.1e0 ) BIO_FIT=0.1e0

  END FUNCTION BIOFIT_R4
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: BIOFIT_R8
!
! !DESCRIPTION: Function BioFit computes the light correction used in the
!  dry deposition and canopy NOx modules.
!\\
!\\
! !INTERFACE:
!
  FUNCTION BIOFIT_R8( COEFF1, XLAI1, SUNCOS1, CFRAC1, NPOLY ) RESULT ( BIO_FIT )
!
! !INPUT PARAMETERS:
!
    REAL*8,  INTENT(IN) :: COEFF1(NPOLY)   ! Baldocchi drydep coefficients
    REAL*8,  INTENT(IN) :: XLAI1           ! Leaf area index [cm2/cm2]
    REAL*8,  INTENT(IN) :: SUNCOS1         ! Cosine( Solar Zenith Angle )
    REAL*8,  INTENT(IN) :: CFRAC1          ! Cloud fraction [unitless]
    INTEGER, INTENT(IN) :: NPOLY           ! # of drydep coefficients
!
! !RETURN VALUE:
!
    REAL*8              :: BIO_FIT         ! Resultant light correction
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
!  09 Dec 2014 - R. Yantosca - Now use BIO_FIT as the return value
!  11 Dec 2014 - M. Yannetti - Split from BIO_FIT
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
    CALL SUNPARAM_R8(TERM(2))
    K=0
    DO K3=1,KK
       DO K2=K3,KK
          DO K1=K2,KK
             K=K+1
             REALTERM(K)=TERM(K1)*TERM(K2)*TERM(K3)
          END DO
       END DO
    END DO

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%% KLUDGE: UNCOMMENT THIS CODE TO GET IDENTICAL RESULTS w/ v10-01e!
!%%%%% (bmy, myannetti, 12/10/14)
!%%%%%
!%%%%% NOTE: This code did not use Fortran "d" exponents, so it was not
!%%%%% forcing values to REAL*8 precision.  This was not as precise as it
!%%%%% could have been, but this code was like this for many years.
!%%%%%    BIO_FIT=0
!%%%%%    DO K=1,NPOLY
!%%%%%       BIO_FIT=BIO_FIT+COEFF1(K)*REALTERM(K)
!%%%%%    END DO
!%%%%%    IF (BIO_FIT .LT. 0.1 ) BIO_FIT=0.1
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! Now explicitly use REAL*8 precision.  This will cause very small
    ! differences at the level of numerical noise when comparing to
    ! prior states of the code like v10-01e.  But this is something that
    ! we can live with.  Stick with REAL*8 precision for now, but we'll
    ! try to implement flexible precision into this routine at a later
    ! point. (bmy, myannetti, 12/10/14)
    BIO_FIT = 0d0
    DO K = 1, NPOLY
       BIO_FIT = BIO_FIT + COEFF1(K)*REALTERM(K)
    END DO
    IF ( BIO_FIT .LT. 0.1d0 ) BIO_FIT=0.1d0

  END FUNCTION BIOFIT_R8
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: SunParam_r4
!
! !DESCRIPTION: Subroutine SUNPARAM is called by BIOFIT to perform the
!  light correction used in the dry deposition and canopy NOx modules.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE SUNPARAM_R4( X )
!
! !DEFINED PARAMETERS:
!
    INTEGER, PARAMETER    :: NN = 3  ! # of variables (LAI, SUNCOS, CLDFRC)
!
! !INPUT/OUTPUT PARAMETERS:
!
    REAL*4, INTENT(INOUT) :: X(NN)   ! LAI, SUNCOS, or cloud fraction
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
!  11 Dec 2014 - M. Yannetti - Split into R4 and R8 versions.
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
    REAL*4 X0(NN),XLOW
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

  END SUBROUTINE SUNPARAM_R4
!EOC


!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: SunParam_r8
!
! !DESCRIPTION: Subroutine SUNPARAM is called by BIOFIT to perform the
!  light correction used in the dry deposition and canopy NOx modules.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE SUNPARAM_R8( X )
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
!  11 Dec 2014 - M. Yannetti - Split into R4 and R8 versions.
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

  END SUBROUTINE SUNPARAM_R8
!EOC
END MODULE DryDep_ToolBox_Mod
