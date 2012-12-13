!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: CMN_DEP_mod
!
! !DESCRIPTION: Common blocks for dry deposition.
!\\
!\\
! !INTERFACE: 
!
      MODULE CMN_DEP_MOD
!
! !USES:
!
      USE CMN_SIZE_MOD, ONLY : IIPAR, JJPAR, NTYPE

      IMPLICIT NONE
      PRIVATE
!
! !PUBLIC DATA MEMBERS:
!
      !----------------------------------------------------------------------
      ! Prior to 12/12/12:
      ! These are now included in State_Met (bmy, 12/12/12)
      !! IREG   - Number of landtypes in grid square (I,J)
      !INTEGER, PUBLIC, ALLOCATABLE :: IREG(:,:)                                
      !
      !! ILAND  - Land type ID in grid square (I,J) for IREG landtypes
      !INTEGER, PUBLIC, ALLOCATABLE :: ILAND(:,:,:)
      !
      !! IUSE   - Fraction ((per mil) of gridbox area occupied by land type
      !INTEGER, PUBLIC, ALLOCATABLE :: IUSE(:,:,:)
      !----------------------------------------------------------------------

      !----------------------------------------------------------------------
      ! Prior to 12/4/12:
      ! This is now included in State_Met (bmy, mpayer, 12/4/12)
      !! Fraction of land in grid box
      !REAL*8,  PUBLIC, ALLOCATABLE :: FRCLND(:,:)
      !----------------------------------------------------------------------

      ! XLAI    - Leaf Area Index of land type for current MONTH
      REAL*8,  PUBLIC, ALLOCATABLE :: XLAI(:,:,:)

      ! XLAI2   - Leaf Area Index of land type for following MONTH
      REAL*8,  PUBLIC, ALLOCATABLE :: XLAI2(:,:,:)
!
! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC :: Init_CMN_DEP
      PUBLIC :: Cleanup_CMN_DEP
!
! !REMARKS:
!  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!  %%%  NOTE: THIS MODULE IS DEPRECATED.  IT WILL BE REMOVED EVENTUALLY.  %%%
!  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! !REVISION HISTORY:
!  23 Aug 2011 - M. Long     - Converted to Module from Header file
!  25 Aug 2011 - R. Yantosca - Change IGLOB,JGLOB to IIPAR,JJPAR
!  04 Dec 2012 - R. Yantosca - Remove FRCLND; it is now included in State_Met
!  12 Dec 2012 - R. Yantosca - IREG, ILAND, IUSE are now included in State_Met
!EOP
!-----------------------------------------------------------------------------
!BOC
      CONTAINS
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_cmn_dep
!
! !DESCRIPTION: Routine INIT\_CMN\_DEP allocates the various arrays
!  given the size parameters from CMN_DEP_MOD.  This allows us to set
!  the size parameters locally for each CPU when running in an MPI
!  environment (e.g. as when connecting to the GEOS-5 GCM).
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE Init_CMN_DEP( am_I_Root, RC )
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
      INTEGER, INTENT(OUT) :: RC          ! Return code
!
! !REMARKS:
!  Need to add better error checking later
! 
! !REVISION HISTORY: 
!  16 Nov 2012 - R. Yantosca - Added ProTeX headers
!  19 Nov 2012 - R. Yantosca - Renamed to INIT_CMN_DEP, added error checks
!  04 Dec 2012 - R. Yantosca - Remove FRCLND; it is now included in State_Met
!  12 Dec 2012 - R. Yantosca - IREG, ILAND, IUSE are now included in State_Met
!EOP
!------------------------------------------------------------------------------
!BOC
      !----------------------------------------------------------------------
      ! Prior to 12/12/12:
      ! These are now included in State_Met (bmy, 12/12/12)
      !ALLOCATE( IREG( IIPAR, JJPAR ), STAT=RC )
      !IF ( RC /= GIGC_SUCCESS ) THEN
      !   WRITE( 6, '(a)' ) 'Cannot allocate IREG (CMN_DEP_mod.F)'
      !   RETURN
      !ENDIF
      !
      !ALLOCATE( ILAND ( IIPAR, JJPAR, NTYPE ), STAT=RC )
      !IF ( RC /= GIGC_SUCCESS ) THEN
      !   WRITE( 6, '(a)' ) 'Cannot allocate ILAND (CMN_DEP_mod.F)'
      !   RETURN
      !ENDIF
      !
      !ALLOCATE( IUSE( IIPAR, JJPAR, NTYPE ), STAT=RC )                    
      !IF ( RC /= GIGC_SUCCESS ) THEN
      !   WRITE( 6, '(a)' ) 'Cannot allocate IUSE (CMN_DEP_mod.F)'
      !   RETURN
      !ENDIF
      !----------------------------------------------------------------------
      
      !------------------------------------------------------------------------
      ! Prior to 12/4/12:
      ! This is now included in State_Met (bmy, mpayer, 12/4/12)
      !ALLOCATE( FRCLND( IIPAR, JJPAR ), STAT=RC )
      !IF ( RC /= GIGC_SUCCESS ) THEN
      !   WRITE( 6, '(a)' ) 'Cannot allocate FRCLND (CMN_DEP_mod.F)'
      !   RETURN
      !ENDIF
      !------------------------------------------------------------------------

      ALLOCATE( XLAI  ( IIPAR, JJPAR, NTYPE ), STAT=RC )
      IF ( RC /= GIGC_SUCCESS ) THEN
         WRITE( 6, '(a)' ) 'Cannot allocate XLAI (CMN_DEP_mod.F)'
         RETURN
      ENDIF

      ALLOCATE( XLAI2 ( IIPAR, JJPAR, NTYPE ), STAT=RC )
      IF ( RC /= GIGC_SUCCESS ) THEN
         WRITE( 6, '(a)' ) 'Cannot allocate XLAI2 (CMN_DEP_mod.F)'
         RETURN
      ENDIF
      
      END SUBROUTINE Init_CMN_DEP
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cleanup_cmn_dep_mod
!
! !DESCRIPTION: Subroutine CLEANUP\_CMN\_DEP\_MOD deallocates all 
!  module arrays.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE Cleanup_CMN_DEP( am_I_Root, RC )
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
      INTEGER, INTENT(OUT) :: RC          ! Return code
!
! !REVISION HISTORY: 
!  19 Nov 2012 - R. Yantosca - Initial version
!  04 Dec 2012 - R. Yantosca - Remove FRCLND; it is now included in State_Met
!  12 Dec 2012 - R. Yantosca - IREG, ILAND, IUSE are now included in State_Met
!EOP
!------------------------------------------------------------------------------
!BOC
      ! Assume success
      RC = GIGC_SUCCESS
      
      ! Deallocate arrays
      !----------------------------------------------------------------------
      ! Prior to 12/12/12:
      ! These are now included in State_Met (bmy, 12/12/12)
      !IF ( ALLOCATED( IREG   ) ) DEALLOCATE( IREG   )
      !IF ( ALLOCATED( ILAND  ) ) DEALLOCATE( ILAND  )
      !IF ( ALLOCATED( IUSE   ) ) DEALLOCATE( IUSE   )
      !---------------------------------------------------------------------
      ! Prior to 12/4/12:
      ! This is now included in State_Met (bmy, 12/4/12)
      !IF ( ALLOCATED( FRCLND ) ) DEALLOCATE( FRCLND )
      !---------------------------------------------------------------------
      IF ( ALLOCATED( XLAI   ) ) DEALLOCATE( XLAI   )
      IF ( ALLOCATED( XLAI2  ) ) DEALLOCATE( XLAI2  )

      END SUBROUTINE Cleanup_CMN_DEP
!EOC
      END MODULE CMN_DEP_MOD
