!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: registry_mod.F90
!
! !DESCRIPTION: Contains derived types and methods to create a registry
!  of each variable contained within a given module.  This will allow the 
!  user to obtain a pointer to each module variable by searching for its name.
!\\
!\\
! !INTERFACE:
!
MODULE Registry_Mod
!
! !USES:
!
  USE Precision_Mod

  IMPLICIT NONE
  PRIVATE
!
! ! PUBLIC TYPES::
!
  !=========================================================================
  ! Derived type for a REGISTRY ITEM (a single registry entry).
  ! This represents a single module variable, plus some metadata.
  !=========================================================================
  TYPE, PUBLIC :: RegItem

     ! Identifying info
     CHARACTER(LEN=67)          :: FullName         ! e.g. "STATE_VARIABLE"
     CHARACTER(LEN=4 )          :: State            ! Name of state
     CHARACTER(LEN=63)          :: Variable         ! Name of variable

     ! Metadata
     CHARACTER(LEN=255)         :: Description      ! Longer description
     REAL(fp)                   :: MemoryInKb       ! Memory use in Kb
     INTEGER                    :: KindVal          ! Numerical KIND value
     INTEGER                    :: Rank             ! Dimensions of data
     CHARACTER(LEN=255)         :: Units            ! Units of data

     ! Pointers to floating-point data
     REAL(fp),          POINTER :: Ptr0d            ! Pointer to 0D data
     REAL(fp),          POINTER :: Ptr1d (:      )  ! Pointer to 1D data
     REAL(fp),          POINTER :: Ptr2d (:,:    )  ! Pointer to 2D dat
     REAL(fp),          POINTER :: Ptr3d (:,:,:  )  ! Pointer to 3D data
     REAL(fp),          POINTER :: Ptr4d (:,:,:,:)  ! Pointer to 4D data

     ! Pointers to integer data
     INTEGER,           POINTER :: Ptr2dI(:,:    )  ! Pointer to 2d int data
     INTEGER,           POINTER :: Ptr3dI(:,:,:  )  ! Pointer to 3d int data

  END TYPE RegItem

  !=========================================================================
  ! Derived type for a METAREGISTRY ITEM (a linked-list of REGISTRY ITEMS)
  !=========================================================================
  TYPE, PUBLIC :: MetaRegItem
     TYPE(MetaRegItem), POINTER :: Next => NULL()   ! Pointer to next node
     TYPE(RegItem    ), POINTER :: Item => NULL()   ! Registry item within
  END TYPE MetaRegItem
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC  :: Registry_AddField
  PUBLIC  :: Registry_Lookup
  PUBLIC  :: Registry_Print
  PUBLIC  :: Registry_Destroy
!
! !PRIVATE MEMBER FUNCTIONS:
!
  PRIVATE :: MetaRegItem_AddNew
  PRIVATE :: MetaRegItem_Create
  PRIVATE :: MetaRegItem_Insert
  PRIVATE :: MetaRegItem_Destroy
  PRIVATE :: Str2Hash
  PRIVATE :: To_UpperCase
!
! !REMARKS:
!  In Fortran 2003, the maximum variable name length is 63 characers, so we
!  have declared the various character name fields of RegItem accordingly.
!
! !REVISION HISTORY:
!  23 Jun 2017 - R. Yantosca - Initial version
!  27 Jun 2017 - R. Yantosca - Added integer data fields, and description
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
! !IROUTINE: Registry_AddField
!
! !DESCRIPTION: Creates a REGISTRY ITEM, which contains information (i.e. 
!  metadata plus a pointer to the data) about a variable within a module.
!  The REGISTRY ITEM will then be added to the METAREGISTRY ITEM, which is 
!  the master list of all variables in the module.  This will allow the user 
!  to obtain a pointer to any variable by searching for its name.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Registry_AddField( am_I_Root, Registry,    State,  Variable,  &
                                RC,        Description, Units,  Data0d,    &
                                Data1d,    Data2d,      Data3d, Data4d,    &
                                Data2dI,   Data3dI                        )
!
! !USES:
!
    USE ErrCode_Mod
!
! !INPUT PARAMETERS:
!
    ! General identifying information
    LOGICAL,           INTENT(IN)       :: am_I_Root        ! Is this root CPU?
    CHARACTER(LEN=*),  INTENT(IN)       :: State            ! Name of state obj
    CHARACTER(LEN=*),  INTENT(IN)       :: Variable         ! Name of variable
    CHARACTER(LEN=*),  OPTIONAL         :: Description      ! Long description
    CHARACTER(LEN=*),  OPTIONAL         :: Units            ! Units of data

    ! Floating-point data targets
    REAL(fp),          OPTIONAL, TARGET :: Data0d           ! Ptr to 0d data
    REAL(fp),          OPTIONAL, TARGET :: Data1d (:      ) ! Ptr to 1d data
    REAL(fp),          OPTIONAL, TARGET :: Data2d (:,:    ) ! Ptr to 2d data
    REAL(fp),          OPTIONAL, TARGET :: Data3d (:,:,:  ) ! Ptr to 3d data 
    REAL(fp),          OPTIONAL, TARGET :: Data4d (:,:,:,:) ! Ptr to 4d data

    ! Integer data targets
    INTEGER,           OPTIONAL, TARGET :: Data2dI(:,:    ) ! Ptr to 2d int data
    INTEGER,           OPTIONAL, TARGET :: Data3dI(:,:,:  ) ! Ptr to 3d int data
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(MetaRegItem), POINTER          :: Registry         ! Registry object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,           INTENT(OUT)      :: RC               ! Success/failure
!
! !REMARKS:
!  Internally, the REGISTRY ITEM will be refered to by its fullname field,
!  which is "STATE_VARIABLE".  Fullname will be defined automatically from
!  the STATE and VARIABLE inputs.
!
! !REVISION HISTORY:
!  23 Jun 2017 - R. Yantosca - Initial version
!  26 Jun 2017 - R. Yantosca - Changed "StateName" to "State", "Name" to
!                              "Variable", and added "MemoryInKb"
!  27 Jun 2017 - R. Yantosca - Now assigns description and KIND value
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!   
    ! Scalars
    REAL(fp)               :: KbPerElement

    ! Strings
    CHARACTER(LEN=255)     :: ErrMsg,   ThisLoc,  TmpFullName
    CHARACTER(LEN=255)     :: TmpState, TmpUnits, TmpVariable, TmpDescription

    ! Objects
    TYPE(RegItem), POINTER :: Item

    !=======================================================================
    ! Initialize
    !=======================================================================
    RC             = GC_SUCCESS
    KbPerElement   = DBLE( fp ) / 1024.0_fp
    ErrMsg         = ''
    ThisLoc        = ' -> at Registry_AddField (in Headers/registry_mod.F90)'
    TmpState       = To_UpperCase( State    )
    TmpVariable    = To_UpperCase( Variable )
    TmpFullname    = TRIM( TmpState ) // '_' // TRIM( TmpVariable )
    TmpDescription = ''
    TmpUnits       = ''
    

    ! Save optional arguments in shadow variables, if passed
    IF ( PRESENT( Description ) ) TmpDescription = Description
    IF ( PRESENT( Units       ) ) TmpUnits       = Units

    !=======================================================================
    ! Allocate the REGISTRY ITEM object, which will hold metadata about
    ! this field, as well as a pointer to the data source, in the registry.
    !=======================================================================
    ALLOCATE( Item, STAT=RC )
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Could not allocate "Item" for variable: ' // TRIM( Variable )
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    
    ! Save identifying info
    Item%FullName      =  TmpFullName
    Item%State         =  TmpState
    Item%Variable      =  TmpVariable
    Item%Description   =  TmpDescription
    Item%Units         =  TmpUnits
    Item%MemoryInKb    =  0.0_fp

    ! Nullify pointers to data etc.
    Item%Ptr0d         => NULL()
    Item%Ptr1d         => NULL()
    Item%Ptr2d         => NULL()
    Item%Ptr3d         => NULL()
    Item%Ptr4d         => NULL()
    Item%Ptr2dI        => NULL()
    Item%Ptr3dI        => NULL()

    ! Depending on the size of the data being passed,
    ! set the rank appropriately and point to the data.
    ! Also compute the size of the array in bytes.
    IF ( PRESENT( Data4d ) ) THEN
       Item%Rank       =  4
       Item%Ptr4d      => Data4d
       Item%MemoryInKb =  KbPerElement * SIZE( Data4d  )
       Item%KindVal    =  KIND( Data4d )
    ELSE IF ( PRESENT( Data3d  ) ) THEN
       Item%Rank       =  3
       Item%Ptr3d      => Data3d
       Item%MemoryInKb =  KbPerElement * SIZE( Data3d  )
       Item%KindVal    =  KIND( Data3d )
    ELSE IF ( PRESENT( Data3dI ) ) THEN
       Item%Rank       =  3
       Item%Ptr3dI     => Data3dI
       Item%MemoryInKb =  KbPerElement * SIZE( Data3dI )
       Item%KindVal    =  KIND( Data3dI )
    ELSE IF ( PRESENT( Data2d  ) ) THEN
       Item%Rank       =  2
       Item%Ptr2d      => Data2d
       Item%MemoryInKb =  KbPerElement * SIZE( Data2d  )
       Item%KindVal    =  KIND( Data2d )
    ELSE IF ( PRESENT( Data2dI ) ) THEN
       Item%Rank       =  2
       Item%Ptr2dI     => Data2dI
       Item%MemoryInKb =  KbPerElement * SIZE( Data2dI )
       Item%KindVal    =  KIND( Data2dI )
    ELSE IF ( PRESENT( Data1d  ) ) THEN
       Item%Rank       =  1
       Item%Ptr1d      => Data1d
       Item%MemoryInKb =  KbPerElement * SIZE( Data1d  )
       Item%KindVal    =  KIND( Data1d )
    ELSE IF ( PRESENT( Data0d  ) ) THEN
       Item%Rank       =  0
       Item%Ptr0d      => Data0d
       Item%MemoryInKb =  KbPerElement
       Item%KindVal    =  KIND( Data1d )
    ELSE
       ErrMsg = 'Need to specify a data source!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    !=======================================================================
    ! Add the REGISTRY ITEM to the METAREGISTRY ITEM, which represents
    ! the list of all data fields contained in a module.
    !=======================================================================
    CALL MetaRegItem_AddNew( am_I_Root, Registry, Item, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Could not add ' // TRIM( TmpVariable ) // ' to the registry!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

  END SUBROUTINE Registry_AddField
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Registry_Lookup
!
! !DESCRIPTION: Get a pointer to any variable in a module (aka "state") by 
!  searching for its name.  Also returns associated metadata.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Registry_Lookup( am_I_Root, Registry,    State,   Variable,    &
                              RC,        Description, KindVal, MemoryInKb,  &
                              Rank,      Units,       Ptr0d,   Ptr1d,       &
                              Ptr2d,     Ptr3d,       Ptr4d,   Ptr2dI,      &
                              Ptr3dI                                       )
!
! !USES:
!
    USE ErrCode_Mod
!
! !INPUT PARAMETERS:
!
    LOGICAL,           INTENT(IN) :: am_I_Root        ! Root CPU?  
    TYPE(MetaRegItem),    POINTER :: Registry         ! Registry obj
    CHARACTER(LEN=*),  INTENT(IN) :: State            ! State name
    CHARACTER(LEN=*),  INTENT(IN) :: Variable         ! Variable name
!                                                     
! !OUTPUT PARAMETERS:                                 
!                                                     
    ! Required outputs
    INTEGER,          INTENT(OUT) :: RC               ! Success or failure

    ! Optional outputs
    CHARACTER(LEN=255),  OPTIONAL :: Description      ! Description of data
    INTEGER,             OPTIONAL :: KindVal          ! Numerical KIND value
    REAL(fp),            OPTIONAL :: MemoryInKb       ! Memory usage
    INTEGER,             OPTIONAL :: Rank             ! Size of data
    CHARACTER(LEN=255),  OPTIONAL :: Units            ! Units of data

    ! Pointers to floating-point data
    REAL(fp),   POINTER, OPTIONAL :: Ptr0d            ! Ptr to 0d data
    REAL(fp),   POINTER, OPTIONAL :: Ptr1d (:      )  ! Ptr to 1d data
    REAL(fp),   POINTER, OPTIONAL :: Ptr2d (:,:    )  ! Ptr to 2d data
    REAL(fp),   POINTER, OPTIONAL :: Ptr3d (:,:,:  )  ! Ptr to 3d data
    REAL(fp),   POINTER, OPTIONAL :: Ptr4d (:,:,:,:)  ! Ptr to 4d data

    ! Pointers to integer data                        
    INTEGER,    POINTER, OPTIONAL :: Ptr2dI(:,:    )  ! Ptr to 2d int data
    INTEGER,    POINTER, OPTIONAL :: Ptr3dI(:,:,:  )  ! Ptr to 3d int data
!
! !REMARKS:
!  Internally, the REGISTRY ITEM will be refered to by its fullname field,
!  which is "STATE_VARIABLE".  Fullname will be defined automatically from
!  the STATE and VARIABLE inputs.
!
! !REVISION HISTORY:
!  23 Jun 2017 - R. Yantosca - Initial version
!  26 Jun 2017 - R. Yantosca - Changed "StateName" to "State", "Name" to
!                              "Variable", and added "MemoryInKb"
!  27 Jun 2017 - R. Yantosca - Also added "Description" and "KindVal" outputs
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    LOGICAL                    :: Is_Description, Is_KindVal
    LOGICAL                    :: Is_MemoryInKb,  Is_Rank
    LOGICAL                    :: Is_Units
    LOGICAL                    :: Is_0d,          Is_0dI
    LOGICAL                    :: Is_1d,          Is_1dI
    LOGICAL                    :: Is_2d,          Is_2dI
    LOGICAL                    :: Is_3d,          Is_3dI
    LOGICAL                    :: Is_4d,          Is_4dI
    INTEGER                    :: FullHash,       ItemHash

    ! Strings
    CHARACTER(LEN=31)          :: FullName31, ItemName31
    CHARACTER(LEN=255)         :: TmpName,    TmpFullName
    CHARACTER(LEN=255)         :: ErrMsg,     ThisLoc

    ! Objects
    TYPE(MetaRegItem), POINTER :: Current

    !=======================================================================
    ! Initialize
    !=======================================================================
    RC             =  GC_SUCCESS
    Current        => NULL()
    ErrMsg         =  ''
    ThisLoc        =  ' -> at Registry_Lookup (in Headers/registry_mod.F90)'

    ! Construct a hash for the full name (i.e. "State_Variable"
    TmpFullName    =  TRIM( State ) // '_' // TRIM( Variable )
    FullName31     =  To_UpperCase( TmpFullName )
    FullHash       =  Str2Hash( FullName31 )

    !=======================================================================   
    ! Test if the optional variables are present outside of the main loop.
    !=======================================================================

    ! Floating-point data
    Is_0d          =  PRESENT( Ptr0d       )
    Is_1d          =  PRESENT( Ptr1d       )
    Is_2d          =  PRESENT( Ptr2d       )
    Is_3d          =  PRESENT( Ptr3d       )
    Is_4d          =  PRESENT( Ptr4d       )

    ! Integer data
    Is_2dI         =  PRESENT( Ptr2dI      )
    Is_3dI         =  PRESENT( Ptr3dI      )

    ! Metadata
    Is_Description =  PRESENT( Description )
    Is_KindVal     =  PRESENT( KindVal     )
    Is_MemoryInKb  =  PRESENT( MemoryInKb  )
    Is_Rank        =  PRESENT( Rank        )
    Is_Units       =  PRESENT( Units       )

    !=======================================================================
    ! Search for the specified field in the Registry
    !=======================================================================

    ! Point to head of linked list
    Current => Registry
    
    ! As long as this entry of the linked list isn't NULL
    DO WHILE( ASSOCIATED( Current ) )

       ! Construct a hash for the full name of this REGISTRY ITEM 
       ItemName31 = Current%Item%FullName
       ItemHash   = Str2Hash( ItemName31 ) 

       ! If the name-hashes match ...
       IF ( FullHash == ItemHash ) THEN

          ! Return rank, units and memory usage if found
          IF ( Is_Description ) Description = Current%Item%Description
          IF ( Is_KindVal     ) KindVal     = Current%Item%KindVal
          IF ( Is_MemoryInKb  ) MemoryInKb  = Current%Item%MemoryInKb
          IF ( Is_Rank        ) Rank        = Current%Item%Rank
          IF ( Is_Units       ) Units       = Current%Item%Units

          ! Then return a pointer to the field
          SELECT CASE( Current%Item%Rank ) 

             ! 4-d data
             CASE( 4 ) 
                IF ( Is_4d ) THEN
                   Ptr4d  => Current%Item%Ptr4d
                   RETURN
                ENDIF

             ! 3-d data
             CASE( 3 ) 
                IF ( Is_3d ) THEN
                   Ptr3d  => Current%Item%Ptr3d
                   RETURN
                ELSE IF  ( Is_3dI ) THEN
                   Ptr3dI => Current%Item%Ptr3dI
                   RETURN
                ENDIF
                   
             ! 2-d data
             CASE( 2 )
                IF ( Is_2d ) THEN
                   Ptr2d  => Current%Item%Ptr2d
                   RETURN
                ELSE IF ( Is_2dI ) THEN
                   Ptr2dI => Current%Item%Ptr2dI
                   RETURN
                ENDIF

             ! 1-d data
             CASE( 1 )
                IF ( Is_1d ) THEN
                   Ptr1d => Current%Item%Ptr1d
                   RETURN
                ENDIF

             ! 0-d data
             CASE( 0 )
                IF ( Is_0d ) THEN
                   Ptr0d => Current%Item%Ptr0d
                   RETURN
                ENDIF

             CASE DEFAULT
                ErrMsg = 'Pointer to data was not passed from calling routine!'
                CALL GC_Error( ErrMsg, RC, ThisLoc )
                RETURN

          END SELECT
       ENDIF
       
       ! Point to next node for next iteration
       Current => Current%Next
    ENDDO

    !=======================================================================
    ! Cleanup and quit
    !=======================================================================
    Current => NULL()

  END SUBROUTINE Registry_Lookup
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Registry_Print
!
! !DESCRIPTION: Prints each REGISTRY ITEM belonging to a METAREGISTRY ITEM.
!  In other words, this prints information about each field contained within
!  a module.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Registry_Print( am_I_Root, Registry, RC, ShortFormat )
!
! !USES:
!
    USE ErrCode_Mod
!
! !INPUT PARAMETERS:
!
    LOGICAL,           INTENT(IN)  :: am_I_Root     ! Are we on the root CPU?
    LOGICAL,           OPTIONAL    :: ShortFormat   ! Print less information
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(MetaRegItem), POINTER     :: Registry      ! Registry of state fields
!
! !OUTPUT PARAMETERS:
!
    INTEGER,           INTENT(OUT) :: RC            ! Success or failure?
!

! !REVISION HISTORY:
!  23 Jun 2017 - R. Yantosca - Initial version
!  26 Jun 2017 - R. Yantosca - Also print memory usage in Kb
!  27 Jun 2017 - R. Yantosca - Now print numeric KIND value and description
!  29 Jun 2017 - R. Yantosca - Add SHORTFORMAT option to print less output
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    LOGICAL                    :: Use_ShortFormat

    ! Strings
    CHARACTER(LEN=255)         :: ErrMsg, ThisLoc

    ! String arrays
    CHARACTER(LEN=4)           :: DimStr(0:4) = &
                                  (/ '0   ', 'x   ', 'xy  ', 'xyz ', 'xyzn' /)

    ! Objects
    TYPE(MetaRegItem), POINTER :: Current
    TYPE(RegItem    ), POINTER :: Item

    !=======================================================================
    ! Initialize
    !=======================================================================

    ! Only print information on the root CPU
    IF ( .not. am_I_Root ) RETURN

    ! Initialize fields
    RC      =  GC_SUCCESS
    ErrMsg  =  ''
    ThisLoc =  ' -> at Registry_Print (in Headers/registry_mod.F90)'
    Current => NULL()
    Item    => NULL()

    ! Save the value of ShortFormat (if passed) in a shadow variable 
    IF ( PRESENT( ShortFormat ) ) THEN
       Use_ShortFormat = ShortFormat
    ELSE
       Use_ShortFormat = .FALSE.
    ENDIF

    ! Print header line
    PRINT*
    PRINT*, 'Registered variables contained within the State_Met object:'
    PRINT*, REPEAT( '=', 79 )

    !=======================================================================
    ! Print information about each state field stored in the Registry
    !=======================================================================

    ! Point to the head node of the Registry
    Current => Registry

    ! As long as the current node isn't NULL
    DO WHILE( ASSOCIATED( Current ) ) 
       
       ! Get the REGISTRY ITEM belonging to this node of the Registry
       Item => Current%Item

       ! Only print on the root CPU
       IF ( ASSOCIATED( Item ) ) THEN

          IF ( Use_ShortFormat ) THEN

             !--------------------------------------------------------------
             ! Just print the name, description, dimension, and units
             !--------------------------------------------------------------
             WRITE( 6, 100 ) Item%FullName,     Item%Description, &
                             DimStr(Item%Rank), Item%Units
  100        FORMAT( 2x, a20, ' | ', a40, ' | ', a4, ' | ', a15 )

          ELSE

             !--------------------------------------------------------------
             ! Print full information about this REGISTRY ITEM
             !--------------------------------------------------------------

             ! Print identifying information
             PRINT*, REPEAT( '-', 70 )
             PRINT*, 'FullName     : ', TRIM( Item%FullName    )
             PRINT*, 'State        : ', TRIM( Item%State       )
             PRINT*, 'Variable     : ', TRIM( Item%Variable    )
             PRINT*, 'Description  : ', TRIM( Item%Description )
             PRINT*, 'Units        : ', TRIM( Item%Units       )
             PRINT*, 'KIND value   : ', Item%KindVal
             PRINT*, 'Memory (Kb)  : ', Item%MemoryInKb
             PRINT*, 'Rank of data : ', Item%Rank, '(', DimStr(Item%Rank), ')'

             ! 0d data
             IF ( ASSOCIATED( Item%Ptr0d ) ) THEN
                PRINT*, 'Value Ptr0d  : ', Item%Ptr0d
             ENDIF

             ! 1d data
             IF ( ASSOCIATED( Item%Ptr1d ) ) THEN
                PRINT*, 'Min value    : ', MINVAL( Item%Ptr1d     )
                PRINT*, 'Max value    : ', MAXVAL( Item%Ptr1d     )
                PRINT*, 'Total        : ', SUM   ( Item%Ptr1d     )
                PRINT*, 'Dimensions   : ', SIZE  ( Item%Ptr1d     )
             ENDIF

             ! 2d data
             IF ( ASSOCIATED( Item%Ptr2d ) ) THEN
                PRINT*, 'Min value    : ', MINVAL( Item%Ptr2d     )
                PRINT*, 'Max value    : ', MAXVAL( Item%Ptr2d     )
                PRINT*, 'Total        : ', SUM   ( Item%Ptr2d     )
                PRINT*, 'Dimensions   : ', SIZE  ( Item%Ptr2d, 1  ), &
                                           SIZE  ( Item%Ptr2d, 2  )
             ENDIF

             ! 2d data -- Integer
             IF ( ASSOCIATED( Item%Ptr2dI ) ) THEN
                PRINT*, 'Min value    : ', MINVAL( Item%Ptr2dI    )
                PRINT*, 'Max value    : ', MAXVAL( Item%Ptr2dI    )
                PRINT*, 'Total        : ', SUM   ( Item%Ptr2dI    )
                PRINT*, 'Dimensions   : ', SIZE  ( Item%Ptr2dI, 1 ), &
                                           SIZE  ( Item%Ptr2dI, 2 )
             ENDIF

             ! 3d data
             IF ( ASSOCIATED( Item%Ptr3d ) ) THEN
                PRINT*, 'Min value    : ', MINVAL( Item%Ptr3d    )
                PRINT*, 'Max value    : ', MAXVAL( Item%Ptr3d    )
                PRINT*, 'Total        : ', SUM   ( Item%Ptr3d    )
                PRINT*, 'Dimensions   : ', SIZE  ( Item%Ptr3d, 1 ), &
                                           SIZE  ( Item%Ptr3d, 2 ), &
                      SIZE  ( Item%Ptr3d, 3 )
             ENDIF

             ! 3d data -- Integer
             IF ( ASSOCIATED( Item%Ptr3dI ) ) THEN
                PRINT*, 'Min value    : ', MINVAL( Item%Ptr3dI    )
                PRINT*, 'Max value    : ', MAXVAL( Item%Ptr3dI    )
                PRINT*, 'Total        : ', SUM   ( Item%Ptr3dI    )
                PRINT*, 'Dimensions   : ', SIZE  ( Item%Ptr3dI, 1 ), &
                                           SIZE  ( Item%Ptr3dI, 2 ), &
                                           SIZE  ( Item%Ptr3dI, 3 )
             ENDIF

             ! 4d data
             IF ( ASSOCIATED( Item%Ptr4d ) ) THEN
                PRINT*, 'Min value    : ', MINVAL( Item%Ptr4d     )
                PRINT*, 'Max value    : ', MAXVAL( Item%Ptr4d     )
                PRINT*, 'Total        : ', SUM   ( Item%Ptr4d     )
                PRINT*, 'Dimensions   : ', SIZE  ( Item%Ptr4d, 1  ), &
                                           SIZE  ( Item%Ptr4d, 2  ), &
                                           SIZE  ( Item%Ptr4d, 3  ), &
                                           SIZE  ( Item%Ptr4d, 4  )
             ENDIF
          ENDIF

       ENDIF

       ! Point to next node of the Registry
       Current => Current%Next
    ENDDO

    !=======================================================================
    ! Cleanup and quit
    !=======================================================================
    Current => NULL()
    Item    => NULL()
       
  END SUBROUTINE Registry_Print
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Registry_Destroy
!
! !DESCRIPTION: Destroys a METAREGISTRY ITEM (i.e. a linked list of REGISTRY
!  ITEMS), each of which contains information (i.e. metadata, plus a pointer
!  to the data source) for a field contained within a module.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Registry_Destroy( am_I_Root, Registry, RC )
!
! !USES:
!
    USE ErrCode_Mod
!
! !INPUT PARAMETERS:
!
    LOGICAL,           INTENT(IN)  :: am_I_Root   ! Are we on the root CPU?
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(MetaRegItem), POINTER     :: Registry    ! Registry of state fields
!
! !OUTPUT PARAMETERS:
!
    INTEGER,           INTENT(OUT) :: RC          ! Success or failure?
!
! !REVISION HISTORY:
!  23 Jun 2017 - R. Yantosca - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
    ! Strings
    CHARACTER(LEN=255) :: ErrMsg, ThisLoc

    !=======================================================================
    ! Initialize
    !=======================================================================
    RC      = GC_SUCCESS
    ErrMsg  = ''
    ThisLoc = ' -> at Registry_Destroy (in Headers/registry_mod.F90)'

    !=======================================================================
    ! Destroy each REGISTRY ITEM contained in the registry,
    ! then destroy the registry itself.
    !=======================================================================
    CALL MetaRegItem_Destroy( am_I_Root, Registry, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Could not destroy the registry object!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

  END SUBROUTINE Registry_Destroy
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: MetaRegItem_AddNew
!
! !DESCRIPTION: Wrapper for methods MetaRegItem_Create and 
!  MetaRegItem_Insert.  Will create a METAREGISTRY ITEM (containing a
!  REGISTRY ITEM) and (1) set it as the head node of a new linked list, or
!  (2) append it to an existing linked list.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE MetaRegItem_AddNew( am_I_Root, Node, Item, RC )
!
! !USES:
!
    USE ErrCode_Mod
!
! !INPUT PARAMETERS: 
!
    LOGICAL,           INTENT(IN)  :: am_I_Root  ! Are we on the root CPU?
    TYPE(RegItem),     POINTER     :: Item       ! REGISTRY ITEM object
!
! !INPUT/OUTPUT PARAMETERS: 
!
    TYPE(MetaRegItem), POINTER     :: Node       ! METAREGISTRY ITEM object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,           INTENT(OUT) :: RC         ! Success or failure
!
! !REVISION HISTORY:
!  23 Jun 2017 - R. Yantosca - Initial version, based on code by Arjen Markus
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
    RC      = GC_SUCCESS
    ErrMsg  = ''
    ThisLoc = ' -> at MetaRegItem_AddNew (in Headers/registry_mod.F90)'
 
    !=======================================================================
    ! Test if the METAREGISTRY ITEM (aka "Node") has been allocated memory
    ! and is therefore part of an existing linked list
    !=======================================================================
    IF ( .not. ASSOCIATED( Node ) ) THEN

       !--------------------------------------------------------------------
       ! If not, then create a new METAREGISTRY ITEM (named "Node"),
       ! and set it at the head of a new linked list
       !--------------------------------------------------------------------
       CALL MetaRegItem_Create( am_I_Root, Node, Item, RC )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Could not create "Node" as the head node of a list!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

    ELSE

       !--------------------------------------------------------------------
       ! Otherwise, create a new METAREGISTRY ITEM (named "Node"),
       ! and append it to the list, immediately following the head node
       !--------------------------------------------------------------------
       CALL MetaRegItem_Insert( am_I_Root, Node, Item, RC )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Could not insert "Node" into an existing linked list!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

    ENDIF

  END SUBROUTINE MetaRegItem_AddNew
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: MetaRegItem_Create
!
! !DESCRIPTION: This method creates a new METAREGISTRY ITEM (to contain the
!  supplied REGISTRY ITEM) and sets it as the head node of a linked list.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE MetaRegItem_Create( am_I_Root, Node, Item, RC )
!
! !USES:
!
    USE ErrCode_Mod
!
! !INPUT PARAMETERS: 
!
    LOGICAL,           INTENT(IN)  :: am_I_Root  ! Are we on the root CPU?
    TYPE(RegItem),     POINTER     :: Item       ! REGISTRY ITEM object
!
! !INPUT/OUTPUT PARAMETERS: 
!
    TYPE(MetaRegItem), POINTER     :: Node       ! METAREGISTRY ITEM object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,           INTENT(OUT) :: RC         ! Success or failure
!
! !REMARKS:
!  This method is not intended to be called directly, but is rather
!  wrapped by the MetaRegItem_AddNew method.
!
! !REVISION HISTORY:
!  23 Jun 2017 - R. Yantosca - Initial version, based on code by Arjen Markus
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
    RC      = GC_SUCCESS
    ErrMsg  = ''
    ThisLoc = ' -> at MetaRegItem_Create (in Headers/registry_mod.F90)'
    
    !=======================================================================
    ! Initialize the METAREGISTRY ITEM itself
    !=======================================================================

    ! Allocate memory
    ALLOCATE( Node, STAT=RC )
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Could not allocate "Node"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
   ENDIF

    ! Nullify the pointer to the next node (it doesn't exist yet)
    Node%Next => NULL()

    !=======================================================================
    ! Initialize the field that will store the REGISTRY ITEM
    !=======================================================================

    ! Because this is the first METAREGISTRY ITEM that is being created,
    ! we can consider this to be the head node of a linked list.
    IF ( .not. ASSOCIATED( Node%Item ) ) THEN 
       ALLOCATE( Node%Item, STAT=RC )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Could not allocate "Node%Item"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
    ENDIF

    ! Attach the given REGISTRY ITEM to the METAREGISTRY ITEM 
    ! (i.e. place it into the head node of a linked list)
    Node%Item = Item

  END SUBROUTINE MetaRegItem_Create
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: MetaRegItem_Insert
!
! !DESCRIPTION: Creates a new METAREGISTRY ITEM (to contain the supplied
!  REGISTRY ITEM), and pops it into an existing linked list, immediately
!  following the head node.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE MetaRegItem_Insert( am_I_Root, Node, Item, RC )
!
! !USES:
!
    USE ErrCode_Mod
!
! !INPUT PARAMETERS: 
!
    LOGICAL,           INTENT(IN)  :: am_I_Root ! Are we on the root CPU?
    TYPE(RegItem),     POINTER     :: Item      ! REGISTRY ITEM object
!
! !INPUT/OUTPUT PARAMETERS: 
!
    TYPE(MetaRegItem), POINTER     :: Node      ! METAREGISTRY ITEM object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,           INTENT(OUT) :: RC        ! Success or failure
!
! !REMARKS:
!  This method is not intended to be called directly, but is rather
!  wrapped by the MetaRegItem_AddNew method.
!
! !REVISION HISTORY:
!  23 Jun 2017 - R. Yantosca - Initial version, based on code by Arjen Markus
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Strings
    CHARACTER(LEN=255)         :: ErrMsg, ThisLoc
    
    ! Objects
    TYPE(MetaRegItem), POINTER :: Next

    !=======================================================================
    ! Initialize
    !=======================================================================
    RC      = GC_SUCCESS
    ErrMsg  = ''
    ThisLoc = ' -> at MetaRegItem_Insert (in Headers/registry_mod.F90)'
    
    !=======================================================================
    ! Initialize a METAREGISTRY ITEM named "Next", which will be inserted 
    ! into the existing list.  "Next" will contain a new REGISTRY ITEM.
    !=======================================================================

    ! Allocate the "Next" object
    ALLOCATE( Next, STAT=RC )
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Could not allocate "Next"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    ! Allocate the "Next%Item" field, which will hold the REGISTRY ITEM
    ALLOCATE( Next%Item, STAT=RC )
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Could not allocate "Next%Item"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    !=======================================================================
    ! Insert "Next" into the existing linked list
    !=======================================================================

    ! Pop the "Next" object in between the current node (i.e. "Node")
    ! and the node that is currently following it (i.e. "Node%Next")
    Next%Next => Node%Next

    ! Now make sure that the current node (i.e. "Node") 
    ! considers that "Next" to be the next node in the list.
    Node%Next => Next

    ! Now that we have inserted the META REGISTRY ITEM "Next" into 
    ! the list, we can save the REGISTRY ITEM into its "Item" field.
    Next%Item =  Item

  END SUBROUTINE MetaRegItem_Insert
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: MetaRegItem_Destroy
!
! !DESCRIPTION:  This method will destroy the REGISTRY ITEM belonging to
!  each METAREGISTRY ITEM (aka node) of a linked list.  It will then destroy
!  each METAREGISTRY ITEM itself.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE MetaRegItem_Destroy( am_I_Root, List, RC )
!
! !USES:
!
    USE ErrCode_Mod
!
! !INPUT PARAMETERS:
!
    LOGICAL,           INTENT(IN)  :: am_I_Root  ! Are we on the root CPU?
!
! !INPUT/OUTPUT PARAMETERS: 
!
    TYPE(MetaRegItem), POINTER     :: List       ! List of METAREGISTRY ITEMS
!
! !OUTPUT PARAMETERS: 
!
    INTEGER,           INTENT(OUT) :: RC         ! Success or failure?
!
! !REMARKS:
!
! !REVISION HISTORY:
!  23 Jun 2017 - R. Yantosca - Initial version, based on code by Arjen Markus
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Strings
    CHARACTER(LEN=255)          :: ErrMsg,  ThisLoc

    ! Objects
    TYPE(MetaRegItem), POINTER  :: Current, Node

    !=======================================================================
    ! Initialize
    !=======================================================================
    RC      =  GC_SUCCESS
    Current => NULL()
    Node    => NULL()
    ErrMsg  =  ''
    ThisLoc =  ' -> at MetaRegItem_Destroy (in Headers/registry_mod.F90)'

    !=======================================================================
    ! Destroy each METAREGISTRY ITEM in the linked list
    !=======================================================================

    ! Point to the first METAREGISTRY ITEM (aka node) in the list
    Node => List

    ! As long as it doesn't point to NULL()
    DO WHILE ( ASSOCIATED( Node ) )

       ! Set the CURRENT pointer to the current METAREGISTRY ITEM
       Current => Node

       ! Free pointers to floating point data in this REGISTRY ITEM
       Current%Item%Ptr0d  => NULL()
       Current%Item%Ptr1d  => NULL()
       Current%Item%Ptr2d  => NULL()
       Current%Item%Ptr3d  => NULL()
       Current%Item%Ptr4d  => NULL()

       ! Free pointers to integer data in this REGISTRY ITEM
       Current%Item%Ptr2dI => NULL()
       Current%Item%Ptr3dI => NULL()
       
       ! Destroy the REGISTRY ITEM itself
       IF ( ASSOCIATED( Current%Item ) ) DEALLOCATE( Current%Item )

       ! Point to the next METAREGISTRY ITEM for the next iteration
       Node => Current%Next

       ! And destroy the current METAREGISTRY ITEM in the list
       DEALLOCATE( Current, STAT=RC )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Cannot deallocate the "Current" METAREGISTRY ITEM!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

    ENDDO

    !=======================================================================
    ! Cleanup and quit
    !=======================================================================
    Current => NULL()
    Node    => NULL()

  END SUBROUTINE MetaRegItem_Destroy
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Str2Hash
!
! !DESCRIPTION: Returns a unique integer hash for a given character string.
!  This allows us to implement a fast name lookup algorithm.
!\\
!\\
! !INTERFACE:
!
  FUNCTION Str2Hash( Str ) RESULT( Hash )
!
! !INPUT PARAMETERS:
!
    CHARACTER(LEN=31), INTENT(IN) :: Str    ! String (31 chars long)
!
! !RETURN VALUE:
!
    INTEGER                       :: Hash   ! Hash value from string
!
! !REMARKS:
!  (1) Algorithm taken from this web page:
!       https://fortrandev.wordpress.com/2013/07/06/fortran-hashing-algorithm/
!
!  (2) For now, we only use the first 31 characers of the character string
!       to compute the hash value.  Most GEOS-Chem variable names only use
!       up to 31 unique characters.  We can change this later if need be.
!
! !REVISION HISTORY:
!  26 Jun 2017 - R. Yantosca - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Initialize
    Hash = 5381

    !-----------------------------------------------------------------------
    ! Construct the hash from the first 31 characters of the string,
    ! which is about the longest variable name for GEOS-Chem.
    !
    ! NOTE: It's MUCH faster to explicitly write these statements
    ! instead of writing them using a DO loop.
    !-----------------------------------------------------------------------
    Hash = ( ISHFT( Hash, 5 ) + Hash ) + ICHAR( Str( 1: 1) )
    Hash = ( ISHFT( Hash, 5 ) + Hash ) + ICHAR( Str( 2: 2) )
    Hash = ( ISHFT( Hash, 5 ) + Hash ) + ICHAR( Str( 3: 3) )
    Hash = ( ISHFT( Hash, 5 ) + Hash ) + ICHAR( Str( 4: 4) )
    Hash = ( ISHFT( Hash, 5 ) + Hash ) + ICHAR( Str( 5: 5) )
    Hash = ( ISHFT( Hash, 5 ) + Hash ) + ICHAR( Str( 6: 6) )
    Hash = ( ISHFT( Hash, 5 ) + Hash ) + ICHAR( Str( 7: 7) )
    Hash = ( ISHFT( Hash, 5 ) + Hash ) + ICHAR( Str( 8: 8) )
    Hash = ( ISHFT( Hash, 5 ) + Hash ) + ICHAR( Str( 9: 9) )
    Hash = ( ISHFT( Hash, 5 ) + Hash ) + ICHAR( Str(10:10) )
    Hash = ( ISHFT( Hash, 5 ) + Hash ) + ICHAR( Str(11:11) )
    Hash = ( ISHFT( Hash, 5 ) + Hash ) + ICHAR( Str(12:12) )
    Hash = ( ISHFT( Hash, 5 ) + Hash ) + ICHAR( Str(13:13) )
    Hash = ( ISHFT( Hash, 5 ) + Hash ) + ICHAR( Str(14:14) )
    Hash = ( ISHFT( Hash, 5 ) + Hash ) + ICHAR( Str(15:15) )
    Hash = ( ISHFT( Hash, 5 ) + Hash ) + ICHAR( Str(16:16) )
    Hash = ( ISHFT( Hash, 5 ) + Hash ) + ICHAR( Str(17:17) )
    Hash = ( ISHFT( Hash, 5 ) + Hash ) + ICHAR( Str(18:18) )
    Hash = ( ISHFT( Hash, 5 ) + Hash ) + ICHAR( Str(19:19) )
    Hash = ( ISHFT( Hash, 5 ) + Hash ) + ICHAR( Str(20:20) )
    Hash = ( ISHFT( Hash, 5 ) + Hash ) + ICHAR( Str(21:21) )
    Hash = ( ISHFT( Hash, 5 ) + Hash ) + ICHAR( Str(22:22) )
    Hash = ( ISHFT( Hash, 5 ) + Hash ) + ICHAR( Str(23:23) )
    Hash = ( ISHFT( Hash, 5 ) + Hash ) + ICHAR( Str(24:24) )
    Hash = ( ISHFT( Hash, 5 ) + Hash ) + ICHAR( Str(25:25) )
    Hash = ( ISHFT( Hash, 5 ) + Hash ) + ICHAR( Str(26:26) )
    Hash = ( ISHFT( Hash, 5 ) + Hash ) + ICHAR( Str(27:27) )
    Hash = ( ISHFT( Hash, 5 ) + Hash ) + ICHAR( Str(28:28) )
    Hash = ( ISHFT( Hash, 5 ) + Hash ) + ICHAR( Str(29:29) )
    Hash = ( ISHFT( Hash, 5 ) + Hash ) + ICHAR( Str(30:30) )
    Hash = ( ISHFT( Hash, 5 ) + Hash ) + ICHAR( Str(31:31) )

  END FUNCTION Str2Hash
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: To_Uppercase
!
! !DESCRIPTION: Converts a string to uppercase, so that we can reliably
!  do string matching.
!\\
!\\
! !INTERFACE:
!
  FUNCTION To_UpperCase( Text ) RESULT( UpCaseText)
!
! !INPUT PARAMETERS: 
!
    CHARACTER(LEN=*), INTENT(IN) :: Text         ! Input test
!
! !RETURN VALUE:
!
    CHARACTER(LEN=255)           :: UpCaseText   ! Output text, uppercase
!
! !REMARKS:
!  Code originally from routine TRANUC (Author: R. D. Stewart, 19 May 1992)
!
! !REVISION HISTORY:
!  26 Jun 2017 - R. Yantosca - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER :: C, Ascii 
    
    !=======================================================================
    ! Convert to uppercase
    !=======================================================================

    ! Initialize
    UpCaseText = Text

    ! Loop over all characters
    DO C = 1, LEN( UpCaseText )

       ! Get the ASCII code for each character
       Ascii = ICHAR( UpCaseText(C:C) )
       
       ! If lowercase, convert to uppercase
       IF ( Ascii > 96 .and. Ascii < 123 ) THEN
          UpCaseText(C:C) = CHAR( Ascii - 32 )
       ENDIF
    ENDDO

  END FUNCTION To_UpperCase
!EOC
END MODULE Registry_Mod
