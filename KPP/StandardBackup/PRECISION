
MODULE gckpp_Precision

!
! Definition of different levels of accuracy
! for REAL variables using KIND parameterization
!
! KPP SP - Single precision kind
  INTEGER, PARAMETER :: sp = SELECTED_REAL_KIND(6,30)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%% KLUDGE FOR FLEX PRECISION (bmy, 12/22/16)
!%%% Because the precision of variables is set when you build KPP, in order
!%%% to get KPP to run in single precision, we have to manually toggle the 
!%%% "DP" parameter to a 4-byte real KIND.  This will be alleviated when
!%%% we build KPP on-the-fly in GEOS-Chem. (bmy, 12/22/16)
!----------------------------------------------------------
! Prior to 12/22/16:
! KPP DP - Double precision kind
!  INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(14,300)
!----------------------------------------------------------
#if defined( USE_REAL8 )
  ! Set KPP DP to double precision (the default)
  INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(14,300)
#else
  ! Set KPP DP to single precision (when compiling with PRECISION=4)
  INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(6,30)
#endif
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! KPP QP - Quadruple precision kind
  INTEGER, PARAMETER :: qp = SELECTED_REAL_KIND(18,400)

END MODULE gckpp_Precision


