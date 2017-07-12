!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: histitem_mod.F90
!
! !DESCRIPTION: Contains types and methods to create a HISTORY ITEM object.  
!  A HISTORY ITEM represents a single GEOS-Chem diagnostic quantity that 
!  will be archived to netCDF file output.
!\\
!\\
! !INTERFACE:
!
MODULE HistItem_Mod
!
! !USES:
!
  USE History_Params_Mod
  USE Precision_Mod

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: HistItem_Create
  PUBLIC :: HistItem_Print
  PUBLIC :: HistItem_Destroy
!
! !PUBLIC TYPES:
!
  !=========================================================================
  ! This is the derived type for a SINGLE HISTORY ITEM OBJECT, which will
  ! hold a quantity from GEOS-Chem that we want to save to netCDF output.
  !=========================================================================
  TYPE, PUBLIC :: HistItem

     !----------------------------------------------------------------------
     ! Identifying information
     !----------------------------------------------------------------------
     CHARACTER(LEN=255) :: Name                  ! Item name
     INTEGER            :: Id                    ! Item Id
     INTEGER            :: ContainerId           ! Container Id

     !----------------------------------------------------------------------
     ! netCDF variable attributes (for COARDS-compliance)
     !----------------------------------------------------------------------
     CHARACTER(LEN=255) :: LongName              ! Item description
     CHARACTER(LEN=255) :: Units                 ! Units of data
     REAL(f4)           :: AddOffset             
     REAL(f4)           :: MissingValue         
     REAL(f4)           :: ScaleFactor          

     !----------------------------------------------------------------------
     ! Pointers to the data in State_Chm, State_Diag, or State_Met
     !----------------------------------------------------------------------
     REAL(fp), POINTER  :: Source_0d             ! Ptr to 0D flex-prec data
     REAL(f4), POINTER  :: Source_0d_4           ! Ptr to 0D 4-byte    data
     INTEGER,  POINTER  :: Source_0d_I           ! Ptr to 0D integer   data

     REAL(fp), POINTER  :: Source_1d  (:    )    ! Ptr to 1D flex-prec data
     REAL(f4), POINTER  :: Source_1d_4(:    )    ! Ptr to 1D 4-byte    data
     INTEGER,  POINTER  :: Source_1d_I(:    )    ! Ptr to 1D integer   data

     REAL(fp), POINTER  :: Source_2d  (:,:  )    ! Ptr to 2D flex-prec data
     REAL(f4), POINTER  :: Source_2d_4(:,:  )    ! Ptr to 2D 4-byte    data
     INTEGER,  POINTER  :: Source_2d_I(:,:  )    ! Ptr to 2D integer   data

     REAL(fp), POINTER  :: Source_3d  (:,:,:)    ! Ptr to 3D flex-prec data
     REAL(f4), POINTER  :: Source_3d_4(:,:,:)    ! Ptr to 3D 4-byte    data
     INTEGER,  POINTER  :: Source_3d_I(:,:,:)    ! Ptr to 3D integer   data

     !----------------------------------------------------------------------
     ! Data arrays 
     !----------------------------------------------------------------------
     INTEGER            :: SpaceDim              ! # of dims (0-3)
     REAL(f4)           :: Data_0d               ! 0D scalar 
     REAL(f4), POINTER  :: Data_1d(:    )        ! 1D vector
     REAL(f4), POINTER  :: Data_2d(:,:  )        ! 2D array
     REAL(f4), POINTER  :: Data_3d(:,:,:)        ! 3D array

     !----------------------------------------------------------------------
     ! Data archival
     !----------------------------------------------------------------------
     REAL(f4)           :: nUpdates              ! # of times updated
     INTEGER            :: Operation             ! Operation code
                                                 !  0=nothing
                                                 !  1=add
                                                 !  2=multiply
   END TYPE HistItem
!
! !REMARKS:
!  Linked list routines taken from original code (linkedlist.f90) 
!  by Arjen Markus; http://flibs.sourceforge.net/linked_list.html
!
! !REVISION HISTORY:
!  13 Jun 2017 - R. Yantosca - Initial version, based on code by Arjen Markus
!  06 Jul 2017 - R. Yantosca - Add source pointers to 4-byte and integer data
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
! !IROUTINE: HistItem_Create
!
! !DESCRIPTION: Initializes a single history item that will be archived 
!  via History (and eventually sent to netCDF output).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HistItem_Create( am_I_Root,    Item,         Id,          &
                              ContainerId,  Name,         RC,          &
                              LongName,     Units,        SpaceDim,    &
                              NX,           NY,           NZ,          & 
                              AddOffset,    MissingValue, ScaleFactor )
!
! !USES:
!
  USE ErrCode_Mod
!
! !INPUT PARAMETERS: 
!
    ! Required arguments
    LOGICAL,           INTENT(IN)  :: am_I_Root       ! Root CPU?
    INTEGER,           INTENT(IN)  :: Id              ! History item Id #
    INTEGER,           INTENT(IN)  :: ContainerId     ! Container Id #
    CHARACTER(LEN=*),  INTENT(IN)  :: Name            ! Item's short name
    CHARACTER(LEN=*),  INTENT(IN)  :: LongName        ! Item's long name
    CHARACTER(LEN=*),  INTENT(IN)  :: Units           ! Units of the data
    INTEGER,           INTENT(IN)  :: SpaceDim        ! Dimension of data
    
    ! Optional arguments
    INTEGER,           OPTIONAL    :: NX              ! # boxes in X-direction
    INTEGER,           OPTIONAL    :: NY              ! # boxes in Y-direction
    INTEGER,           OPTIONAL    :: NZ              ! # boxes in Z-direction
    REAL(f4),          OPTIONAL    :: AddOffset       ! COARDS-compliant 
    REAL(f4),          OPTIONAL    :: MissingValue    !  attributes for 
    REAL(f4),          OPTIONAL    :: ScaleFactor     !  netCDF output
!
! !INPUT/OUTPUT PARAMETERS: 
!
    TYPE(HistItem),    POINTER     :: Item            ! Item to be archived
!
! !OUTPUT PARAMETERS: 
!
    INTEGER,           INTENT(OUT) :: RC              ! Success or failure
!
! !REMARKS:
!  (1) We need to copy string data to a temporary string of length 255 
!       characters, or else Gfortran will choke.
!
! !REVISION HISTORY:
!  13 Jun 2017 - R. Yantosca - Initial version

!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Strings
    CHARACTER(LEN=255) :: ErrMsg, ThisLoc, TempStr

    !========================================================================
    ! Initialize
    !========================================================================

    ! Assume success
    RC      = GC_SUCCESS

    ! For error output
    ErrMsg  = ''
    ThisLoc = ' -> at HistItem_Create (in History/histitem_mod.F90)'

    ! Allocate the Item object
    ALLOCATE( Item, STAT=RC )
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Cannot allocate the "Item" object!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    
    !========================================================================
    ! Required inputs, handle these first
    !========================================================================

    !--------------------------------------------
    ! ID of this Item
    !--------------------------------------------
    IF ( Id >= 0 ) THEN
       Item%Id = Id
    ELSE
       ErrMsg = '"Id" cannot be negative!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    !--------------------------------------------
    ! Id of the container this item belongs to
    !--------------------------------------------
    IF ( ContainerId >= 0 ) THEN
       Item%ContainerId = ContainerId
    ELSE
       ErrMsg = '"ContainerId" cannot be negative!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    !--------------------------------------------
    ! Name (aka "short name")
    !--------------------------------------------
    IF ( LEN_TRIM( Name ) > 0 ) THEN
       TempStr   = Name
       Item%Name = TempStr
    ELSE
       ErrMsg = 'You must specify a value for "Name"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    !--------------------------------------------
    ! LongName
    !--------------------------------------------
    IF ( LEN_TRIM( LongName ) > 0 ) THEN
       TempStr       = LongName
       Item%LongName = TempStr
    ELSE
       ErrMsg = 'You must specify a value for "LongName"'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    !--------------------------------------------
    ! Units
    !--------------------------------------------
    IF ( LEN_TRIM( Units ) > 0 ) THEN
       TempStr    = Units
       Item%Units = TempStr
    ELSE
       ErrMsg = 'You must specify a value for "Units"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    !=======================================================================
    ! Data fields: Allocate data fields (0-3 dimensions)
    !=======================================================================
    IF ( SpaceDim < 0 .or. SpaceDim > 3 ) THEN
       ErrMsg = 'SpaceDim must be in the range 0-3!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ELSE
       Item%SpaceDim = SpaceDim
    ENDIF

    ! Nullify pointers to source fields
    Item%Source_0d   => NULL()
    Item%Source_0d_4 => NULL()
    Item%Source_0d_I => NULL()
    Item%Source_1d   => NULL()
    Item%Source_1d_4 => NULL()
    Item%Source_1d_I => NULL()
    Item%Source_2d   => NULL()
    Item%Source_2d_4 => NULL()
    Item%Source_2d_I => NULL()
    Item%Source_3d   => NULL()
    Item%Source_3d_4 => NULL()
    Item%Source_3d_I => NULL()
    
    ! Nullify data fields
    Item%Data_0d     =  0.0_f4
    Item%Data_1d     => NULL()
    Item%Data_2d     => NULL()
    Item%Data_3d     => NULL()

    ! Allocate data field, based on SpaceDim
    SELECT CASE( Item%SpaceDim )

      CASE( 3 )
          ALLOCATE( Item%Data_3d( NX, NY, NZ ), STAT=RC )
          IF ( RC == GC_SUCCESS ) THEN
             Item%Data_3d = 0.0_f4
          ELSE
             ErrMsg = 'Could not allocate "Item%Data_3d" array!'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF

       CASE( 2 )
          ALLOCATE( Item%Data_2d( NX, NY ), STAT=RC )
          IF ( RC == GC_SUCCESS ) THEN
             Item%Data_2d = 0.0_f4
          ELSE
             ErrMsg = 'Could not allocate "Item%Data_2d" array!'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN           
          ENDIF

       CASE( 1 )
          ALLOCATE( Item%Data_1d( NX ), STAT=RC )
          IF ( RC == GC_SUCCESS ) THEN
             Item%Data_1d  = 0.0_f4
          ELSE
             ErrMsg = 'Could not allocate "Item%Data_1d" array!'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN             
          ENDIF

       CASE DEFAULT
          ! Do Nothing

    END SELECT

    !========================================================================
    ! Optional inputs, handle these next
    !========================================================================

    !--------------------------------------------
    ! Add_Offset
    !--------------------------------------------
    IF ( PRESENT( AddOffset ) ) THEN
       Item%AddOffset = AddOffset
    ELSE
       Item%AddOffset = 0.0_f4
    ENDIF

    !--------------------------------------------
    ! MissingValue
    !--------------------------------------------
    IF ( PRESENT( MissingValue ) ) THEN
       Item%MissingValue = MissingValue
    ELSE
       Item%MissingValue = UNDEFINED
    ENDIF

    !--------------------------------------------
    ! Scale_Factor
    !--------------------------------------------
    IF ( PRESENT( ScaleFactor ) ) THEN
       Item%ScaleFactor = ScaleFactor
    ELSE
       Item%ScaleFactor = 1.0_f4
    ENDIF

  END SUBROUTINE HistItem_Create
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HistItem_Print
!
! !DESCRIPTION: Prints information contained within a single history item.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HistItem_Print( am_I_Root, Item, RC, ShortFormat )
!
! !USES:
!
    USE ErrCode_Mod
!
! !INPUT PARAMETERS: 
!
    LOGICAL,        INTENT(IN)  :: am_I_Root    ! Are we on the root CPU?
    TYPE(HistItem), POINTER     :: Item         ! History Item
    LOGICAL,        OPTIONAL    :: ShortFormat  ! Print truncated format
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT) :: RC           ! Success or failure?
!
! !REVISION HISTORY:
!  13 Jun 2017 - R. Yantosca - Initial version
!  06 Jul 2017 - R. Yantosca - Add option to print truncated output format
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARAIBLES
!
    ! Scalars
    LOGICAL          :: Use_ShortFormat

    ! String arrays
    CHARACTER(LEN=4) :: DimStr(0:4) = &
                        (/ '0   ', 'x   ', 'xy  ', 'xyz ', 'xyzn' /)

    !=======================================================================
    ! Initialize
    !=======================================================================

    ! Assume success
    RC = GC_SUCCESS

    ! Define
    IF ( PRESENT( ShortFormat ) ) THEN
       Use_ShortFormat = ShortFormat
    ELSE
       Use_ShortFormat = .FALSE.
    ENDIF

    !=======================================================================
    ! Print information about this HISTORY ITEM (only on the root CPU)
    !=======================================================================
    IF ( ASSOCIATED( Item ) .and. am_I_Root ) THEN

       IF ( Use_ShortFormat ) THEN

          !-----------------------------------------------------------------
          ! Use truncated output format
          !-----------------------------------------------------------------
          WRITE( 6, 100 ) TRIM( Item%Name ),     TRIM( Item%LongName ),  &
                          DimStr(Item%SpaceDim), TRIM( Item%Units    )
 100      FORMAT( 2x, a20, ' | ', a40, ' | ', a4, ' | ', a15 )

       ELSE

          !-----------------------------------------------------------------
          ! Use expanded output format
          !-----------------------------------------------------------------
          PRINT*, REPEAT( '-', 70 )
          PRINT*, 'Name           : ', TRIM( Item%Name      )
          PRINT*, 'Long_Name      : ', TRIM( Item%LongName )
          PRINT*, 'Units          : ', TRIM( Item%Units     )
          PRINT*, 'Add_Offset     : ', Item%AddOffset
          PRINT*, 'Missing_Value  : ', Item%MissingValue
          PRINT*, 'Scale_Factor   : ', Item%ScaleFactor
          PRINT*, ''
          PRINT*, 'Id             : ', Item%ID
          PRINT*, 'CollectionId   : ', Item%ContainerId
          PRINT*, ''
          PRINT*, 'nUpdates       : ', Item%nUpdates
          PRINT*, 'Operation      : ', Item%Operation
          PRINT*, ''
          PRINT*, 'SpaceDim       : ', Item%SpaceDim, ' (',         &
                                       DimStr(Item%SpaceDim), ')'
          PRINT*, 'Data_0d        : ', Item%Data_0d
       
          IF ( ASSOCIATED( Item%Data_1d ) ) THEN 
             PRINT*, 'Min   Data_1d  : ', MINVAL( Item%Data_1d    )
             PRINT*, 'Max   Data_1d  : ', MAXVAL( Item%Data_1d    )
             PRINT*, 'Total Data_1d  : ', SUM   ( Item%Data_1d    )
             PRINT*, 'Size  Data_1d  : ', SIZE  ( Item%Data_1d    )
          ENDIF
          
          IF ( ASSOCIATED( Item%Data_2d ) ) THEN
             PRINT*, 'Min   Data_2d  : ', MINVAL( Item%Data_2d    )
             PRINT*, 'Max   Data_2d  : ', MAXVAL( Item%Data_2d    )
             PRINT*, 'Total Data_2d  : ', SUM   ( Item%Data_2d    )
             PRINT*, 'Size  Data_2d  : ', SIZE  ( Item%Data_2d, 1 ), &
                                          SIZE  ( Item%Data_2d, 2 )
          ENDIF

          IF ( ASSOCIATED( Item%Data_3d ) ) THEN
             PRINT*, 'Min   Data_3d  : ', MINVAL( Item%Data_3d    )
             PRINT*, 'Max   Data_3d  : ', MAXVAL( Item%Data_3d    )
             PRINT*, 'Total Data_3d  : ', SUM   ( Item%Data_3d    )
             PRINT*, 'Size  Data_3d  : ', SIZE  ( Item%Data_3d, 1 ), &
                                          SIZE  ( Item%Data_3d, 2 ), &
                                          SIZE  ( Item%Data_3d, 3 )
          ENDIF
       ENDIF
    ENDIF

  END SUBROUTINE HistItem_Print
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP

! !IROUTINE: HistItem_Destroy( Item )
!
! !DESCRIPTION: Deallocates all pointer-based array fields of the history
!  item, then destroys the history item itself.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HistItem_Destroy( am_I_Root, Item, RC )
!
! !USES:
!
  USE ErrCode_Mod
  USE History_Params_Mod
!
! !INPUT PARAMETERS:
!
    LOGICAL,        INTENT(IN)  :: am_I_Root   ! Are we on the root CPU?
!
! !INPUT/OUTPUT PARAMETERS: 
!
    TYPE(HistItem), POINTER     :: Item        ! History item
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT) :: RC          ! Success or failure
!
! !REMARKS:
!
! !REVISION HISTORY:
!  13 Jun 2017 - R. Yantosca - Initial version
!  06 Jul 2017 - R. Yantosca - Nullify source pointers to 4-byte & integer data
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Strings
    CHARACTER(LEN=255) :: ErrMsg, ThisLoc

    !=======================================================================
    ! Initialize
    !=======================================================================
    ErrMsg = ''
    ThisLoc = ' -> at HistItem_Destroy (in History/histitem_mod.F90)'

    !=======================================================================
    ! Nullify fields that are just pointing to other objects   
    !======================================================================
    Item%Source_0d   => NULL()
    Item%Source_0d_4 => NULL()
    Item%Source_0d_I => NULL()
    Item%Source_1d   => NULL()
    Item%Source_1d_4 => NULL()
    Item%Source_1d_I => NULL()
    Item%Source_2d   => NULL()
    Item%Source_2d_4 => NULL()
    Item%Source_2d_I => NULL()
    Item%Source_3d   => NULL()
    Item%Source_3d_4 => NULL()
    Item%Source_3d_I => NULL()

    !=======================================================================
    ! Free allocated fields
    !=======================================================================
    IF ( ASSOCIATED( Item%Data_3d ) ) THEN
       DEALLOCATE( Item%Data_3d, STAT=RC )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Could not deallocate "Item%Data_3d" array'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
    ENDIF 
      
    IF ( ASSOCIATED( Item%Data_2d ) ) THEN
       DEALLOCATE( Item%Data_2d, STAT=RC )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Could not deallocate "Item%Data_2d" array'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
    ENDIF       

    IF ( ASSOCIATED( Item%Data_1d ) ) THEN
       DEALLOCATE( Item%Data_1d, STAT=RC )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Could not deallocate "Item%Data_1d" array'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
    ENDIF

    !=======================================================================
    ! Free the History Item itself
    !=======================================================================
    IF ( ASSOCIATED( Item ) ) THEN
       DEALLOCATE( Item, STAT=RC )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Could not deallocate "Item"'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
    ENDIF

  END SUBROUTINE HistItem_Destroy
!EOC
END MODULE HistItem_Mod

