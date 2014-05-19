!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: fjx_acet_mod
!
! !DESCRIPTION: \subsection*{Overview}
!  This module contains functions used for the new acetone pressure
!  dependency calculation in JRATET.f introduced in FAST-JX version 6.7
!  This is a hack to effectively implement Fast-JX v7.0b acetone
!  photolysis into Fast-J. See use in JRATET.f
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
      PUBLIC :: QQA
      PUBLIC :: QQB
!
! !AUTHOR:
! Original code from Michael Prather.
! Implemented into GEOS-Chem by Claire Carouge (ccarouge@seas.harvard.edu)
!
! !REVISION HISTORY:
!  20 Apr 2009 - C. Carouge - Created the module from fastJX64.f code.
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!  19 May 2014 - M. Sulprizio- Update acetone photolysis to Fast-JX v7.0b
!                              (S.D. Eastham)
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
! !IROUTINE: QQA
!
! !DESCRIPTION: 
!\\
!\\
! !INTERFACE:
!
      subroutine QQA(PP,QQQT,K)
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
         ! Taking only the last 7 bins from Fast-JX (!)
         QQQ(:,1) = (/ 1.980d-20,  5.927d-21,  6.000d-22,  5.868d-23,
     &                 5.934d-25,  0.000d0,    0.000d0              /)
         QQQ(:,2) = (/ 1.240d-20,  4.464d-21,  7.146d-22,  1.171d-22,
     &                 2.202d-24,  0.000d0,    0.000d0              /)
         QQQ(:,3) = (/ 9.213d-21,  3.702d-21,  7.100d-22,  1.357d-22,
     &                 3.115d-24,  0.000d0,    0.000d0              /)
      endif
      call X_interp_FJX (PP,QQQT, TQQ(1),QQQ(K,1),
     &        TQQ(2),QQQ(K,2), TQQ(3),QQQ(K,3), 3)

      end subroutine QQA
!EOC 
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: QQB
!
! !DESCRIPTION: 
!\\
!\\
! !INTERFACE:
!
      subroutine QQB(TT,QQQT,K)
!
! !USES:
!
      implicit none
!
! !INPUT PARAMETERS: 
!
      real*8,  intent(in)  :: TT
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
         ! Temperature at which cross-sections calculated
         TQQ = (/235.0d0,260.0d0,298.0d0/)
         ! Taking only the last 7 bins from Fast-JX (!)
         QQQ(:,1) = (/ 1.158d-22,  2.648d-23,  6.014d-24,  1.502d-24,
     &                 4.211d-26,  0.000d0,    0.000d0              /)
         QQQ(:,2) = (/ 5.664d-22,  1.681d-22,  4.919d-23,  1.477d-23,
     &                 5.602d-25,  0.000d0,    0.000d0              /)
         QQQ(:,3) = (/ 2.804d-21,  1.092d-21,  4.079d-22,  1.496d-22,
     &                 7.707d-24,  0.000d0,    0.000d0              /)
      endif
      call X_interp_FJX (TT,QQQT, TQQ(1),QQQ(K,1), 
     &        TQQ(2),QQQ(K,2), TQQ(3),QQQ(K,3), 3)

      end subroutine QQB
!EOC 
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: x_interp_fjx
!
! !DESCRIPTION: Up-to-three-point linear interpolation function for X-sections
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
      END MODULE FJX_ACET_MOD 

