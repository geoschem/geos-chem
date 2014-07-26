!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: CMN_NOX_mod
!
! !DESCRIPTION: CMN\_NOX\_mod is the module file for containing NOx from soils.
!\\
!\\
! !INTERFACE: 
!
      MODULE CMN_NOX_MOD
!
! !USES:
!
      USE CMN_SIZE_MOD, ONLY : IIPAR, JJPAR

      IMPLICIT NONE
      PRIVATE
!
! !PUBLIC DATA MEMBERS:
! 
      ! NOTE: Keep this for backwards compatibility for now (bmy, 10/2/07)
      ! GEMISNOX2 = Soil Nox [molec NOx/cm3/s]
      REAL*8, PUBLIC, ALLOCATABLE :: GEMISNOX2(:,:)
!
! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC :: Init_CMN_NOX
      PUBLIC :: Cleanup_CMN_NOX

! !REMARKS:
!  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!  !%%%  NOTE: THIS MODULE IS DEPRECATED.  WE WILL PROBABLY PLACE THE  %%%
!  !%%%  MODULE ARRAYS INTO A NEW DERIVED TYPE FOR THE GIGC EMISSIONS  %%%
!  !%%%  COMPONENT. (bmy, 11/19/12)                                    %%%
!  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! !REVISION HISTORY:
!  05 Mar 1998 - M. Schultz - Initial version
!  (1 ) Changed RCS ID tags from "C" to "!" to allow 
!        freeform compilation.  (bmy, 6/25/02)
!  (2 ) Moved BXHEIGHT to "dao_mod.f".  The fact that BXHEIGHT
!        was in "CMN_NOX" is historical baggage. (bmy, 9/18/02)
!  (3 ) Now everything except GEMISNOX, GEMISNOX2 is in
!        "lightning_mod.f" (bmy, 4/14/04)
!  (4 ) Remove GEMISNOX from common block (ltm, bmy, 10/2/07)
!  23 Aug 2011 - M. Long    - Converted to Module from Header file
!-----------------------------------------------------------------------------
!BOC
      CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_cmn_nox
!
! !DESCRIPTION: Subroutine INIT\_CMN\_NOX allocates all module arrays.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE Init_CMN_NOX( am_I_Root, RC )
!
! !USES:
!
      USE GIGC_ErrCode_Mod
!
! !INPUT PARAMETERS:
!
      LOGICAL, INTENT(IN)  :: am_I_Root ! Are we on the root CPU?
!
! !OUTPUT PARAMETERS:
!
      INTEGER, INTENT(OUT) :: RC ! Success or failure?
! 
! !REVISION HISTORY: 
!  19 Nov 2012 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC

      ! Assume success
      RC = GIGC_SUCCESS
      
      ! Allocate arrays
      ALLOCATE( GEMISNOX2(IIPAR,JJPAR), STAT=RC )
      IF ( RC /= GIGC_SUCCESS ) THEN
         WRITE( 6, '(a)' ) 'Allocation error: GEMISNOX2 (CMN_NOX_mod.F)'
         RETURN
      ENDIF

      ! Zero arrays
      GEMISNOX2 = 0d0

      END SUBROUTINE Init_CMN_NOX
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cleanup_cmn_nox
!
! !DESCRIPTION: Subroutine CLEANUP\_CMN\_NOX allocates all module arrays.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE Cleanup_CMN_NOX( am_I_Root, RC )
!
! !USES:
!
      USE GIGC_ErrCode_Mod
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)  :: am_I_Root ! Are we on the root CPU?
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT) :: RC ! Success or failure?
! 
! !REVISION HISTORY: 
!  19 Nov 2012 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
      ! Assume success
      RC = GIGC_SUCCESS

      ! Deallocate arrays
      IF ( ALLOCATED( GEMISNOX2 ) ) DEALLOCATE( GEMISNOX2 )

      END SUBROUTINE Cleanup_CMN_NOX
!EOC
      END MODULE CMN_NOX_MOD

