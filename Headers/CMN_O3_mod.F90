#ifdef BPCH_DIAG
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: CMN_O3_mod.F90
!
! !DESCRIPTION: Common blocks for anthro emissions (via SMVGEAR!)
!\\
!\\
! !INTERFACE:
!
MODULE CMN_O3_MOD
!
! !USES:
!
  USE PRECISION_MOD

  IMPLICIT NONE
  PUBLIC
!
! !PUBLIC DATA MEMBERS:
!
  ! SAVEEOH = array to save EOH fields (evf, 5/21/13)
  ! SAVEGLYX= array to save GLYX fields (eam, 8/25/14 )
  ! SAVEOA  = array for total organic aerosol (trc 6 of ND42, eam,7/10/14)
  REAL(fp), ALLOCATABLE :: SAVEEOH(:,:,:)
  REAL(fp), ALLOCATABLE :: SAVEGLYX(:,:,:)
  REAL(fp), ALLOCATABLE :: SAVEOA(:,:,:)
!
! !REMARKS:
!  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!  %%%  NOTE: THIS MODULE WILL BE A STUB UNLESS GEOS-Chem IS COMPILED    %%%
!  %%%  WITH THE BPCH_DIAG=y OPTION. (bmy, 10/4/19)                      %%%
!  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! !REVISION HISTORY:
!  23 Aug 2011 - M. Long   - Converted to Module from Header file
!  See https://github.com/geoschem/geos-chem for complete history
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
! !IROUTINE: Init_Cmn_O3
!
! !DESCRIPTION: Subroutine INIT\_CMN\_O3 allocates all module arrays.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Init_CMN_O3( Input_Opt, State_Grid, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE Input_Opt_Mod,  ONLY : OptInput
    USE State_Grid_Mod, ONLY : GrdState
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)  :: Input_Opt   ! Input Options object
    TYPE(GrdState), INTENT(IN)  :: State_Grid  ! Grid State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT) :: RC          ! Success or failure?
!
! !REVISION HISTORY:
!  19 Nov 2012 - R. Yantosca - Added ProTeX headers
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
    ! Assume success
    RC = GC_SUCCESS

    ! Allocate arrays
    ALLOCATE( SAVEEOH (State_Grid%NX,State_Grid%NY,State_Grid%NZ), STAT=RC )
    ALLOCATE( SAVEGLYX(State_Grid%NX,State_Grid%NY,State_Grid%NZ), STAT=RC )
    ALLOCATE( SAVEOA  (State_Grid%NX,State_Grid%NY,State_Grid%NZ), STAT=RC )

    ! Zero arrays
    SAVEEOH    = 0e+0_fp
    SAVEGLYX   = 0e+0_fp
    SAVEOA     = 0e+0_fp

  END SUBROUTINE Init_CMN_O3
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Cleanup_Cmn_O3
!
! !DESCRIPTION: Subroutine CLEANUP\_CMN\_O3 allocates all module arrays.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Cleanup_CMN_O3( RC )
!
! !USES:
!
    USE ErrCode_Mod
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT) :: RC        ! Success or failure?
!
! !REVISION HISTORY:
!  19 Nov 2012 - R. Yantosca - Added ProTeX headers
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
    ! Assume success
    RC = GC_SUCCESS

    ! Deallocate arrays
    IF ( ALLOCATED( SAVEEOH     ) ) DEALLOCATE( SAVEEOH  )
    IF ( ALLOCATED( SAVEGLYX    ) ) DEALLOCATE( SAVEGLYX )
    IF ( ALLOCATED( SAVEOA      ) ) DEALLOCATE( SAVEOA   )

  END SUBROUTINE Cleanup_CMN_O3
!EOC
END MODULE CMN_O3_MOD
#endif

