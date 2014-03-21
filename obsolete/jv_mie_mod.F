!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: jv_mie_mod.F
!
! !DESCRIPTION: This include file contains physical constants for the
!  GEOS-Chem column chemistry code. 
!\\
!\\
! !INTERFACE: 
!
      MODULE JV_MIE_MOD
!
! !USES:
!
      IMPLICIT NONE
      PUBLIC
!
! !DEFINED PARAMETERS:
!
      ! NL    Maximum number of levels after insertion of extra Mie levels
      ! N__   Number of levels in Mie grid: 2*(2*lpar+2+jaddto(1))+3
      ! M__   Number of Gauss points used
!-----------------------------------------------------------------------
!     NL=1500 was too small for dicarbonyls, so we upped it to 2000.
!     Uncomment this line to restore the previous definition (phs, 10/9/09)
!      INTEGER, PARAMETER :: NL  = 1500
!-----------------------------------------------------------------------
      INTEGER, PARAMETER :: NL  = 2000
      INTEGER, PARAMETER :: N__ = 2*NL
      INTEGER, PARAMETER :: M__ = 4 
!
! !PUBLIC DATA MEMBERS:
!
      ! Arrays
      REAL*8  :: A(M__),           B(M__,M__),   C1(M__)
      REAL*8  :: H(M__),           AA(M__,M__),  CC(M__,M__)
      REAL*8  :: S(M__,M__),       W(M__,M__),   U1(M__,M__)
      REAL*8  :: V1(M__),          WT(M__),      EMU(M__)
      REAL*8  :: PM(M__,2*M__),    PM0(2*M__),   POMEGA(2*M__,N__)
      REAL*8  :: ZTAU(N__),        FZ(N__),      FJ(N__)
      REAL*8  :: DD(M__,M__,N__),  RR(M__,N__)
      REAL*8  :: ZREFL,            ZFLUX

      ! Scalars
      REAL*8  :: RADIUS,           ZU0
      INTEGER :: ND,               N
      INTEGER :: M,                MFIT

      !=================================================================
      ! Declare the following global variables as THREADPRIVATE for the
      ! OpenMP parallelization on all platforms (bmy, 3/23/03)
      !=================================================================
!$OMP THREADPRIVATE( A,B,C1,H,AA,CC,S,W,U1,V1,WT,EMU,PM,PM0,POMEGA )
!$OMP THREADPRIVATE( ZTAU,FZ,FJ,DD,RR,ZREFL,ZFLUX,RADIUS,ZU0 )
!$OMP THREADPRIVATE( ND,N,M,MFIT )
!
! !REMARKS:
!  Keep increasing NL if necessary.  This will avoid the "too many levels
!  in photolysis" error.
!
! !REVISION HISTORY:
!  (1 ) Changed RCS ID tags to by adding a ! comment character to allow
!        freeform compilation.  Also added & continuation characters in 
!        column 73 to allow header files to be included in F90 freeform files.
!        Also changed comment character from "C" to "!", to allow this
!        file to be inlined into freeform source code. (bmy, 6/25/02)
!  (2 ) Now declare common blocks /MIEBLK/ and /MINDEX/ as THREADPRIVATE for
!        all platforms (bmy, 3/23/03)
!  (3 ) Set NL to 1000 to avoid SMVGEAR crash with GEOS-5.2.0 on Sept 1st 2008
!  03 Aug 2011 - M. Long     - Converted from Header file to Module
!EOP
!------------------------------------------------------------------------------
!BOC
	END MODULE JV_MIE_MOD
!EOC
