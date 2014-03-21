!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
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
      USE CMN_SIZE_MOD
      USE COMODE_LOOP_MOD                   ! IDEMS, NEMIS, NCS
      USE CMN_O3_MOD                        ! EMISRR, EMISRRN
      USE CMN_DIAG_MOD                      ! ND15

      IMPLICIT NONE

      PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
! 
      PUBLIC :: Init_VDIFF_PRE
      PUBLIC :: Cleanup_VDIFF_PRE
!
! !PUBLIC DATA MEMBERS:
!
      PUBLIC :: IIPAR, JJPAR, LLPAR         ! from "CMN_SIZE_mod"
      PUBLIC :: IDEMS, NEMIS, NCS, NDRYDEP  ! from "comode_loop_mod"
      PUBLIC :: EMISRR, EMISRRN             ! from "CMN_O3_mod"
      PUBLIC :: ND15, ND44                  ! from "CMN_DIAG_mod"
      PUBLIC :: emis_save

      ! Make sure MAXTRACERS >= N_TRACERS
      INTEGER, PARAMETER :: MAXTRACERS = 100 

      REAL*8, ALLOCATABLE :: emis_save(:,:,:)
!
! !REVISION HISTORY:
!  01 Jun 2009 - C. Carouge & J. Lin - Initial version  
!  07 Oct 2009 - R. Yantosca         - Added CVS Id tag  
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
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
! !IROUTINE: init_vdiff_pre
!
! !DESCRIPTION: Subroutine INIT\_VDIFF\_PRE allocates all module arrays.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE Init_VDIFF_PRE( am_I_Root, RC )
!
! !USES:
!
      USE GIGC_ErrCode_Mod
!
! !INPUT PARAMETERS: 
!
      LOGICAL, INTENT(IN)  :: am_I_Root   ! Are we on the root CPU?
!
! !OUTPUT PARAMETERS:
!
      INTEGER, INTENT(OUT) :: RC          ! Success or failure?
!
! !REMARKS:
!  Need to add error-checking on the allocation statements, so that we
!  exit the code upon error.
! 
! !REVISION HISTORY: 
!  19 Nov 2012 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
      ! Assume success
      RC = GIGC_SUCCESS

      ! Allocate arrays
      ALLOCATE( emis_save( IIPAR, JJPAR, MAXTRACERS ), STAT=RC )
      IF ( RC /= GIGC_SUCCESS ) THEN
         WRITE( 6, '(a)' ) 'Allocation error: EMIS_SAVE'
         RETURN
      ENDIF
      
      ! Zero arrays
      emis_save = 0.d0
      
      END SUBROUTINE Init_VDIFF_PRE
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cleanup_vdiff_pree
!
! !DESCRIPTION: Subroutine CLEANUP\_VDIFF\_PRE deallocates all module arrays.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE Cleanup_VDIFF_PRE( am_I_Root, RC )
!
! !USES:
!
      USE GIGC_ErrCode_Mod

      IMPLICIT NONE
!
! !INPUT PARAMETERS: 
!
      LOGICAL, INTENT(IN)  :: am_I_Root   ! Are we on the root CPU?
!
! !OUTPUT PARAMETERS:
!
      INTEGER, INTENT(OUT) :: RC          ! Success or failure?
!
! !REMARKS:
!
! !REVISION HISTORY: 
!  19 Nov 2012 - R. Yantosca - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
      ! Assume success
      RC = GIGC_SUCCESS

      ! Deallocate arrays
      IF ( ALLOCATED( EMIS_SAVE ) ) DEALLOCATE( EMIS_SAVE )

      END SUBROUTINE Cleanup_VDIFF_PRE
!EOC
      END MODULE VDIFF_PRE_MOD

