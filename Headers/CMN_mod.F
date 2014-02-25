!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: CMN_mod
!
! !DESCRIPTION: Module CMN is the remnant of header file "CMN", which once
!  held many global variables, but now is reduced to only a couple.
!\\
!\\
! !INTERFACE: 
!
      MODULE CMN_MOD
!
! !USES:
!
      USE CMN_SIZE_MOD, ONLY : IIPAR, JJPAR

      IMPLICIT NONE
      PRIVATE
!
! !PUBLIC DATA MEMBERS:
!
      ! LPAUSE  stores the annual mean tropopause (bmy, 12/6/99)
      INTEGER, PUBLIC, ALLOCATABLE :: LPAUSE(:,:)

      ! IFLX stores the flags for ND27 diagnostic (bmy, 12/6/99)
      ! This is mostly obsolete by now (bmy, 8/24/11)
      INTEGER, PUBLIC, ALLOCATABLE :: IFLX(:,:)
!
! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC :: Init_CMN
      PUBLIC :: Cleanup_CMN
!
! !REMARKS:
!  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!  %%% NOTE: THIS MODULE IS DEPRECATED.  IT WILL BE REMOVED IN THE FUTURE  %%% 
!  %%% COMBINE THESE ARRAYS INTO ANOTHER MODULE FILE.  FOR NOW, KEEP IT.   %%%
!  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! !REVISION HISTORY:
!  23 Aug 2011 - M. Long   - Converted to Module from Header file
!EOP
!-----------------------------------------------------------------------------
!BOC
      CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_cmn
!
! !DESCRIPTION: Subroutine INIT\_CMN allocates all module arrays.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE Init_CMN( am_I_Root, RC )
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
! !REVISION HISTORY: 
!  19 Nov 2012 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
      ! Assume success
      RC = GIGC_SUCCESS
      
      ! Allocate arrays
      ALLOCATE( LPAUSE(IIPAR,JJPAR), STAT=RC )
      IF ( RC /= GIGC_SUCCESS ) THEN
         WRITE( 6, '(a)' ) 'Allocation error: LPAUSE (CMN_mod.F)'
         RETURN
      ENDIF

      ALLOCATE( IFLX(IIPAR,JJPAR), STAT=RC ) 
      IF ( RC /= GIGC_SUCCESS ) THEN
         WRITE( 6, '(a)' ) 'Allocation error: IFLX (CMN_mod.F)'
         RETURN
      ENDIF

      ! Zero arrays
      LPAUSE = 0
      IFLX   = 0
      
      END SUBROUTINE Init_CMN
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cleanup_cmn
!
! !DESCRIPTION: Subroutine CLEANUP\_CMN allocates all module arrays.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE Cleanup_CMN( am_I_Root, RC )
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
! !REVISION HISTORY: 
!  19 Nov 2012 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
      ! Assume success
      RC = GIGC_SUCCESS

      ! Deallocate arrays
      IF ( ALLOCATED( LPAUSE ) ) DEALLOCATE( LPAUSE )
      IF ( ALLOCATED( IFLX   ) ) DEALLOCATE( IFLX   )

      END SUBROUTINE Cleanup_CMN
!EOC
      END MODULE CMN_MOD

