!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: fjx_acet_mod
!
! !DESCRIPTION: \subsection*{Overview}
!  This module contains functions used for the new acetone pressure
!  dependency calculation in JRATET.f introduced in FAST-JX version 6.4
!  The temperature interpolation factors and the Xsect are different for
!  both acetone photolysis reactions and interdependant. See use in JRATET.f
!
!\subsection*{Reference}
!  Blitz, M. A., D. E. Heard, M. J. Pilling, S. R. Arnold, M. P. Chipperfield
!   2004: \emph{Pressure and temperature-dependent quantum yields for the 
!   photodissociation of acetone between 279 and 327.5 nm}, 
!   \underline{GRL}, \textbf{31}, 9, L09104.
!\\
!\\
!
! !INTERFACE:
!
      MODULE FJX_ACET_MOD
!
! !USES:
! 
      IMPLICIT NONE
      PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC :: QQ2_F
      PUBLIC :: QQ1_F
      PUBLIC :: TFACA_F
      PUBLIC :: TFAC0_F
      PUBLIC :: TFAC_F
!
! !AUTHOR:
! Original code from Michael Prather.
! Implemented into GEOS-Chem by Claire Carouge (ccarouge@seas.harvard.edu)
!
! !REVISION HISTORY:
!  20 Apr 2009 - C. Carouge - Created the module from fastJX64.f code.
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!EOP
!------------------------------------------------------------------------------
      CONTAINS

!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: tfaca_f
!
! !DESCRIPTION: Calculates temperature interpolation factors for acetone
!\\
!\\
! !INTERFACE:
!
      FUNCTION TFACA_F(TTT, IV)
!
! !USES
!
      USE CMN_FJ_MOD, ONLY : JPMAX, JPPJ
      USE JV_CMN_MOD
!
! !INPUT PARAMETERS:
!
      ! Index of the specie in jv_spec.dat (should be between 4 and NJVAL)
      INTEGER :: IV

      ! Temperature in 1 grid box
      REAL*8  :: TTT
!
! !OUTPUT VALUE:
!
      ! Temperature interpolation factor
      REAL*8  :: TFACA_F 
!                               with the "D" double-precision exponent.
!EOP
!------------------------------------------------------------------------------
!BOC
!
      TFACA_F = (TTT-TQQ(1,IV))/(TQQ(2,IV)-TQQ(1,IV))
      TFACA_F = max(0.d0, min(1.d0, TFACA_F))

      RETURN
      
      END FUNCTION TFACA_F
!EOC 
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: tfac0_f
!
! !DESCRIPTION: Calculates temperature interpolation factors for acetone
!\\
!\\
! !INTERFACE:
!
      FUNCTION TFAC0_F(TTT, IV)
!
! !USES:
!
      USE CMN_FJ_MOD, ONLY : JPMAX, JPPJ
      USE JV_CMN_MOD
!
! !INPUT PARAMETERS: 
!
      ! Index of the specie in jv_spec.dat (should be between 4 and NJVAL)
      INTEGER :: IV

      ! Temperature in 1 grid box
      REAL*8  :: TTT
!
! !OUTPUT VALUE:
!
      ! Temperature interpolation factor
      REAL*8  :: TFAC0_F 
!EOP
!------------------------------------------------------------------------------
!BOC
!
      TFAC0_F = ( (TTT-TQQ(1,IV))/(TQQ(2,IV)-TQQ(1,IV)) )**2
      if (TTT .lt. TQQ(1,IV)) then
         TFAC0_F = (TTT - 210.d0)/(TQQ(1,IV)-210.d0)
      endif
      TFAC0_F = max(0.d0, min(1.d0, TFAC0_F))

      RETURN
      
      END FUNCTION TFAC0_F
!EOC 
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: tfac_f
!
! !DESCRIPTION: Calculates temperature interpolation factors for acetone
!\\
!\\
! !INTERFACE:
!
      FUNCTION TFAC_F(TTT, IV)
!
! !USES:
!
      USE CMN_FJ_MOD, ONLY : JPMAX, JPPJ
      USE JV_CMN_MOD
!
! !INPUT PARAMETERS: 
!
      ! Index of the specie in jv_spec.dat (should be between 4 and NJVAL)
      INTEGER :: IV

      ! Temperature in 1 grid box
      REAL*8  :: TTT

!
! !OUTPUT VALUE:
!
      ! Temperature interpolation factor
      REAL*8  :: TFAC_F 
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL*8  :: TT200

      TT200 = min(300.d0, max(200.d0, TTT))
      TFAC_F = (TT200-TQQ(1,IV))/(TQQ(2,IV)-TQQ(1,IV))

      RETURN
      
      END FUNCTION TFAC_F
!EOC 
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: qq2_f
!
! !DESCRIPTION: This routine computes the cross-section for acetone.
!\\
!\\
! !INTERFACE:
!
      FUNCTION QQ2_F(TFAC0, IV, K, TTT)
!
! !USES:
!
      USE CMN_FJ_MOD, ONLY : JPMAX, JPPJ
      USE JV_CMN_MOD
!
! !INPUT PARAMETERS: 
!
      ! Index of the specie in jv_spec.dat (should be between 4 and NJVAL)
      INTEGER :: IV

      ! Wavelength
      INTEGER :: K

      ! Temperature in 1 grid box
      REAL*8  :: TTT

      ! Temperature interpolation factor from TFAC0_F function
      REAL*8  :: TFAC0
!
! !OUTPUT VALUE:
!
      ! Xsect (total abs) for Acetone
      REAL*8  :: QQ2_F
!
! !NOTES:
!  (1 ) We use IV-3 and not IV because there is no QQQ values for O2, O3 
!        and O1-D. (ccc, 4/20/19)
!EOP
!------------------------------------------------------------------------------
!BOC
!
      QQ2_F  = QQQ(K,1,IV-3) + (QQQ(K,2,IV-3)-QQQ(K,1,IV-3))*TFAC0
      if (TTT .lt. TQQ(1,IV)) then
         QQ2_F = QQQ(K,1,IV-3)*TFAC0
      endif

      RETURN

      END FUNCTION QQ2_F
!EOC 
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: qq1_f
!
! !DESCRIPTION: This routine computes the cross-section for acetone.
!\\
!\\
! !INTERFACE:
!
      FUNCTION QQ1_F(TFAC, IV, K)
!
! !USES:
!
      USE CMN_FJ_MOD, ONLY : JPMAX, JPPJ
      USE JV_CMN_MOD
!
! !INPUT PARAMETERS: 
!
      ! Index of the specie in jv_spec.dat (should be between 4 and NJVAL)
      INTEGER :: IV

      ! Wavelength
      INTEGER :: K

      ! Temperature interpolation factor from TFAC_F function
      REAL*8  :: TFAC
!
! !OUTPUT VALUE:
!
      ! Xsect (total abs) for Acetone
      REAL*8  :: QQ1_F 
!
! !NOTES:
!  (1 ) We use IV-3 and not IV because there is no QQQ values for O2, O3 
!        and O1-D. (ccc, 4/20/19)
!EOP
!------------------------------------------------------------------------------
!BOC
!
      QQ1_F = QQQ(K,1,IV-3) + (QQQ(K,2,IV-3)-QQQ(K,1,IV-3))*TFAC

      RETURN

      END FUNCTION QQ1_F
!EOC
      END MODULE FJX_ACET_MOD 

