! $Id: jv_mie.h,v 1.1 2003/06/30 20:26:08 bmy Exp $
!
!----jv_mie.h-----COMMON BLOCKS for FAST-J code: 4x4x85 (prather 4/96)
!
!  Parameters
!  ----------
!
!     NL    Maximum number of levels after insertion of extra Mie levels
!     N__   Number of levels in Mie grid: 2*(2*lpar+2+jaddto(1))+3
!     M__   Number of Gauss points used
!
!
! NOTES:
! (1 ) Changed RCS ID tags to by adding a ! comment character to allow
!       freeform compilation.  Also added & continuation characters in 
!       column 73 to allow header files to be included in F90 freeform files.
!       Also changed comment character from "C" to "!", to allow this
!       file to be inlined into freeform source code. (bmy, 6/25/02)
! (2 ) Now declare common blocks /MIEBLK/ and /MINDEX/ as THREADPRIVATE for
!       all platforms (bmy, 3/23/03)
!-----------------------------------------------------------------------
      INTEGER    NL, N__, M__
!-----------------------------------------------------------------------
!  NL=250 was too small for the GEOS code, so I upped it to 400.
!  Uncomment this line to restore the original definition (bmy, 9/29/99)
!      PARAMETER (NL=250, N__=2*NL, M__=4)
!-----------------------------------------------------------------------
!  NL=400 was too small again, so we upped it to 500. 
!  Uncomment this line to restore the previous definition (bmy, 9/29/99)
!      PARAMETER (NL=400, N__=2*NL, M__=4)
!-----------------------------------------------------------------------
!  NL=500 was too small again, so we upped it to 750. 
!  Uncomment this line to restore the previous definition (mje, 6/14/01)
!      PARAMETER (NL=500, N__=2*NL, M__=4)
!-----------------------------------------------------------------------
      PARAMETER (NL=750, N__=2*NL, M__=4)
      REAL*8 A,B,C1,H,AA,CC,S,W,U1,V1,WT,EMU,PM,PM0,POMEGA
      REAL*8 ZTAU,FZ,FJ,DD,RR,ZREFL,ZFLUX,RADIUS,ZU0
      INTEGER ND,N,M,MFIT
      COMMON/MIEBLK/ A(M__),B(M__,M__),C1(M__),H(M__),AA(M__,M__),      &
     &   CC(M__,M__),S(M__,M__),W(M__,M__),U1(M__,M__),V1(M__),WT(M__), &
     &   EMU(M__),PM(M__,2*M__),PM0(2*M__),POMEGA(2*M__,N__),ZTAU(N__), &
     &   FZ(N__),FJ(N__),DD(M__,M__,N__),RR(M__,N__),                   &
     &   ZREFL,ZFLUX,RADIUS,ZU0
      COMMON/MINDEX/ ND,N,M,MFIT

      !=================================================================
      ! Declare the following common blocks as THREADPRIVATE for the
      ! OpenMP parallelization on all platforms (bmy, 3/23/03)
      !=================================================================
!$OMP THREADPRIVATE( /MIEBLK/ )
!$OMP THREADPRIVATE( /MINDEX/ )
C-----------------------------------------------------------------------
