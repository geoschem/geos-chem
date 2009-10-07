! $Id: vdiff_pre_mod.f,v 1.2 2009/10/07 14:49:47 bmy Exp $
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: VDIFF_PRE_MOD
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

      IMPLICIT NONE

#     include "CMN_SIZE"
#     include "comode.h"                    ! IDEMS, NEMIS, NCS
#     include "CMN_O3"                      ! EMISRR, EMISRRN
#     include "CMN_DIAG"                    ! ND15

      PRIVATE
!
! !PUBLIC DATA MEMBERS
!
      PUBLIC :: IIPAR, JJPAR, LLPAR         ! from "CMN_SIZE"
      PUBLIC :: IDEMS, NEMIS, NCS, NDRYDEP  ! from "comode.h"
      PUBLIC :: EMISRR, EMISRRN             ! from "CMN_O3"
      PUBLIC :: ND15, ND44                  ! from "CMN_DIAG"
      PUBLIC :: emis_save

      ! Make sure MAXTRACERS >= N_TRACERS
      INTEGER, PARAMETER :: MAXTRACERS = 100 

      REAL*8 :: emis_save(IIPAR, JJPAR, MAXTRACERS) = 0.d0

      END MODULE VDIFF_PRE_MOD
!
! !REVISION HISTORY:
!  01 Jun 2009 - C. Carouge & J. Lin - Initial version  
!  07 Oct 2009 - R. Yantosca         - Added CVS Id tag  
!EOP
!------------------------------------------------------------------------------
