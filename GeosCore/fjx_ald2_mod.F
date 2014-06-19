!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: fjx_ald2_mod
!
! !DESCRIPTION: \subsection*{Overview}
!  This module contains functions used for the new acetaldehyde pressure
!  dependency calculation in JRATET.F
!  This is a hack to effectively implement M. Prather's acetaldehyde
!  photolysis into Fast-J.
!
!\\
!\\
!
! !INTERFACE:
!
      MODULE FJX_ALD2_MOD
!
! !USES:
! 
      IMPLICIT NONE
      PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC :: QQA_ALD2
!
! !AUTHOR:
! Original code from Michael Prather.
! Implemented into GEOS-Chem by Eloise Marais (emarais@seas.harvard.edu)
!
! !REVISION HISTORY:
!  13 Jun 2014 - E. Marais   - Created the module from fjx_acet_mod.F code.
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
! !IROUTINE: QQA_ALD2
!
! !DESCRIPTION: 
!\\
!\\
! !INTERFACE:
!
      subroutine QQA_ALD2(PP,QQQT,K)
!
! !USES:
!
      implicit none
!
! !INPUT PARAMETERS: 
!
      real*8,  intent(in)  :: PP
      integer, intent(in)  :: K
!
! !OUTPUT PARAMETERS: 
!
      real*8,  intent(out) :: QQQT
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      logical,save::FIRST=.TRUE.
      real*8,dimension(7,3),save::QQQ
      real*8,dimension(3),save::TQQ

      if (FIRST) then
         FIRST=.false.
         ! Declare arrays
         ! Pressure at which cross-sections calculated
         TQQ = (/177.0d0,566.0d0,999.0d0/)
         ! Taking only the last 7 bins (consistent with ACET P-dependence)
         QQQ(:,1) = (/ 3.732d-20,  2.707d-20,  1.579d-20,  6.566d-21,
     &                 3.883d-22,  5.683d-26,  0.000d0              /)
         QQQ(:,2) = (/ 2.987d-20,  1.923d-20,  9.497d-21,  3.450d-21,
     &                 1.914d-22,  3.762d-26,  0.000d0              /)
         QQQ(:,3) = (/ 2.479d-20,  1.485d-20,  6.739d-21,  2.319d-21,
     &                 1.258d-22,  2.790d-26,  0.000d0              /)
      endif
      call X_interp_FJX (PP,QQQT, TQQ(1),QQQ(K,1),
     &        TQQ(2),QQQ(K,2), TQQ(3),QQQ(K,3), 3)

      end subroutine QQA_ALD2
!EOC 
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: x_interp_fjx
!
! !DESCRIPTION: Up-to-three-point linear interpolation function for X-sections.
!               This should really be a separate module, as it's used here and
!               in fjx_acet_mod.F.
!\\
!\\
! !INTERFACE:
!
      subroutine X_interp_FJX (TINT,XINT, T1,X1, T2,X2, T3,X3, L123)
!
! !USES:
!
      implicit none
!
! !INPUT PARAMETERS: 
!
      real*8, intent(in)::  TINT,T1,T2,T3, X1,X2,X3
      integer,intent(in)::  L123
!
! !OUTPUT PARAMETERS: 
!
      real*8,intent(out)::  XINT
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      real*8  TFACT

      if (L123 .le. 1) then
           XINT = X1
      elseif (L123 .eq. 2) then
             TFACT = max(0.d0,min(1.d0,(TINT-T1)/(T2-T1) ))
           XINT = X1 + TFACT*(X2 - X1)
      else
        if (TINT.le. T2) then
             TFACT = max(0.d0,min(1.d0,(TINT-T1)/(T2-T1) ))
           XINT = X1 + TFACT*(X2 - X1)
        else
             TFACT = max(0.d0,min(1.d0,(TINT-T2)/(T3-T2) ))
           XINT = X2 + TFACT*(X3 - X2)
        endif
      endif

      END SUBROUTINE X_interp_FJX
!EOC
      END MODULE FJX_ALD2_MOD 

