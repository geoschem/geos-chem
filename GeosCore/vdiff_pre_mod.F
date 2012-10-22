!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: vdiff_pre_mod
!
! !DESCRIPTION: Module VDIFF\_PRE\_MOD contains variables used in VDIFF\_MOD.
!\\
!\\
! !INTERFACE: 
!
      MODULE VDIFF_PRE_MOD
! 
! !USES:
!
      USE TRACER_MOD, ONLY : N_TRACERS

#     include "define.h"
      USE CMN_SIZE_MOD
      USE COMODE_LOOP_MOD                    ! IDEMS, NEMIS, NCS
      USE CMN_O3_MOD                      ! EMISRR, EMISRRN
      USE CMN_DIAG_MOD                    ! ND15

      IMPLICIT NONE

      PRIVATE
!
! !PUBLIC DATA MEMBERS:
!
      PUBLIC :: IIPAR, JJPAR, LLPAR         ! from "CMN_SIZE_mod"
      PUBLIC :: IDEMS, NEMIS, NCS, NDRYDEP  ! from "comode_loop_mod"
      PUBLIC :: EMISRR, EMISRRN             ! from "CMN_O3_mod"
      PUBLIC :: ND15, ND44                  ! from "CMN_DIAG_mod"
      PUBLIC :: emis_save
#if defined( DEVEL )
      PUBLIC :: SET_VDIFF_PRE_MOD
#endif

      ! Make sure MAXTRACERS >= N_TRACERS
      INTEGER, PARAMETER :: MAXTRACERS = 100 

#if defined( DEVEL )
      REAL*8, ALLOCATABLE :: emis_save(:,:,:)!(IIPAR, JJPAR, MAXTRACERS)
#else
      REAL*8 :: emis_save(IIPAR, JJPAR, MAXTRACERS) = 0.d0
#endif
!
! !REVISION HISTORY:
!  01 Jun 2009 - C. Carouge & J. Lin - Initial version  
!  07 Oct 2009 - R. Yantosca         - Added CVS Id tag  
!EOP
!------------------------------------------------------------------------------
!BOC
#if defined( DEVEL )
      CONTAINS

      SUBROUTINE SET_VDIFF_PRE_MOD

! This routine simply allocates VDIFF arrays as initializaion

      IMPLICIT NONE

      INTEGER AS

      ALLOCATE(emis_save(IIPAR, JJPAR, MAXTRACERS), STAT=AS)
      emis_save=0.d0
      
      END SUBROUTINE SET_VDIFF_PRE_MOD
#endif

      END MODULE VDIFF_PRE_MOD
!EOC
