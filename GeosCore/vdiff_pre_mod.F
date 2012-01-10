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
      PUBLIC :: IIPAR, JJPAR, LLPAR         ! from "F77_CMN_SIZE"
      PUBLIC :: IDEMS, NEMIS, NCS, NDRYDEP  ! from "comode.h"
      PUBLIC :: EMISRR, EMISRRN             ! from "F77_CMN_O3"
      PUBLIC :: ND15, ND44                  ! from "F77_CMN_DIAG"
      PUBLIC :: emis_save
      PUBLIC :: SET_VDIFF_PRE_MOD

      ! Make sure MAXTRACERS >= N_TRACERS
      INTEGER, PARAMETER :: MAXTRACERS = 100 

      REAL*8, ALLOCATABLE :: emis_save(:,:,:)!(IIPAR, JJPAR, MAXTRACERS)
!
! !REVISION HISTORY:
!  01 Jun 2009 - C. Carouge & J. Lin - Initial version  
!  07 Oct 2009 - R. Yantosca         - Added CVS Id tag  
!EOP
!------------------------------------------------------------------------------
!BOC
      CONTAINS

      SUBROUTINE SET_VDIFF_PRE_MOD

      IMPLICIT NONE

      INTEGER AS

      ALLOCATE(emis_save(IIPAR, JJPAR, MAXTRACERS), STAT=AS)
      emis_save=0.d0
      
      END SUBROUTINE SET_VDIFF_PRE_MOD

      END MODULE VDIFF_PRE_MOD
!EOC
