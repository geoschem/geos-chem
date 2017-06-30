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

     !----------------------------------------------------------------------
     ! Identifying info
     !----------------------------------------------------------------------
     CHARACTER(LEN=67)    :: FullName         ! e.g. "STATE_VARIABLE"
     CHARACTER(LEN=4 )    :: State            ! Name of state
     CHARACTER(LEN=63)    :: Variable         ! Name of variable

     !----------------------------------------------------------------------
     ! Metadata
     !----------------------------------------------------------------------
     CHARACTER(LEN=255)   :: Description      ! Longer description
     REAL(fp)             :: MemoryInKb       ! Memory use in Kb
     INTEGER              :: KindVal          ! Numerical KIND value
     INTEGER              :: Rank             ! Dimensions of data
     CHARACTER(LEN=255)   :: Units            ! Units of data

     !----------------------------------------------------------------------
     ! Pointers to floating point data (flexible precision)
     !----------------------------------------------------------------------
     REAL(fp), POINTER    :: Ptr0d            ! For 0D flex-prec data
     REAL(fp), POINTER    :: Ptr1d  (:    )   ! For 1D flex-prec data
     REAL(fp), POINTER    :: Ptr2d  (:,:  )   ! For 2D flex-prec data
     REAL(fp), POINTER    :: Ptr3d  (:,:,:)   ! For 3D flex-prec data

     !----------------------------------------------------------------------
     ! Pointers to floating point data (4-byte precision)
     !----------------------------------------------------------------------
     REAL(f4), POINTER    :: Ptr0d_4          ! For 0D 4-byte data
     REAL(f4), POINTER    :: Ptr1d_4(:    )   ! For 1D 4-byte data
     REAL(f4), POINTER    :: Ptr2d_4(:,:  )   ! For 2D 4-byte data
     REAL(f4), POINTER    :: Ptr3d_4(:,:,:)   ! For 3D 4-byte data

     !----------------------------------------------------------------------
     ! Pointers to integer data
     !----------------------------------------------------------------------
     INTEGER,  POINTER    :: Ptr0d_I          ! For 0D int data
     INTEGER,  POINTER    :: Ptr1d_I(:    )   ! For 1D int data
     INTEGER,  POINTER    :: Ptr2d_I(:,:  )   ! For 2D int data
     INTEGER,  POINTER    :: Ptr3d_I(:,:,:)   ! For 3D int data

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
!  29 Jun 2017 - R. Yantosca - Added Ptr1DI to type RegItem
!  30 Jun 2017 - R. Yantosca - Add more pointers to 4-byte and integer data
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
  SUBROUTINE Registry_AddField( am_I_Root, Registry,  State,        &
                                Variable,  RC,        Description,  & 
                                Units,     Data0d,    Data1d,       &
                                Data2d,    Data3d,    Data0d_4,     &
                                Data1d_4,  Data2d_4,  Data3d_4,     &
                                Data0d_I,  Data1d_I,  Data2d_I,     &
                                Data3d_I                           )
!
! !USES:
!
    USE ErrCode_Mod
!
! !INPUT PARAMETERS:
!
    ! General identifying information
    LOGICAL,           INTENT(IN)       :: am_I_Root       ! Root CPU?
    CHARACTER(LEN=*),  INTENT(IN)       :: State           ! State name
    CHARACTER(LEN=*),  INTENT(IN)       :: Variable        ! variable
    CHARACTER(LEN=*),  OPTIONAL         :: Description     ! Long description
    CHARACTER(LEN=*),  OPTIONAL         :: Units           ! Units of data

    ! Floating-point data targets (flexible precision)
    REAL(fp),          OPTIONAL, TARGET :: Data0d          ! 0D flex-prec data
    REAL(fp),          OPTIONAL, TARGET :: Data1d  (:    ) ! 1D flex_prec data
    REAL(fp),          OPTIONAL, TARGET :: Data2d  (:,:  ) ! 2D flex-prec data
    REAL(fp),          OPTIONAL, TARGET :: Data3d  (:,:,:) ! 3D flex-prec data 

    ! Floating-point data targets (4-byte precision)
    REAL(f4),          OPTIONAL, TARGET :: Data0d_4        ! 0D 4-byte data
    REAL(f4),          OPTIONAL, TARGET :: Data1d_4(:    ) ! 1D 4-byte data
    REAL(f4),          OPTIONAL, TARGET :: Data2d_4(:,:  ) ! 2D 4-byte data
    REAL(f4),          OPTIONAL, TARGET :: Data3d_4(:,:,:) ! 3D 4-byte data 

    ! Integer data targets
    INTEGER,           OPTIONAL, TARGET :: Data0d_I        ! 1D int data
    INTEGER,           OPTIONAL, TARGET :: Data1d_I(:    ) ! 1D int data
    INTEGER,           OPTIONAL, TARGET :: Data2d_I(:,:  ) ! 2D int data
    INTEGER,           OPTIONAL, TARGET :: Data3d_I(:,:,:) ! 3D int data
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
!  29 Jun 2017 - R. Yantosca - Added Data1dI for 1D integer data
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
    Item%Ptr0d_4       => NULL()
    Item%Ptr1d_4       => NULL()
    Item%Ptr2d_4       => NULL()
    Item%Ptr3d_4       => NULL()
    Item%Ptr0d_I       => NULL()
    Item%Ptr1d_I       => NULL()
    Item%Ptr2d_I       => NULL()
    Item%Ptr3d_I       => NULL()

    !-----------------------------------------------------------------------
    ! Depending on the size of the data being passed,
    ! set the rank appropriately and point to the data.
    ! Also compute the size of the array in bytes.
    !
    ! Assign pointers of the REGISTRY ITEM to flex-precision data targets
    !-----------------------------------------------------------------------
    IF ( PRESENT( Data3d  ) ) THEN
       Item%Rank       =  3
       Item%Ptr3d      => Data3d
       Item%MemoryInKb =  KbPerElement * SIZE( Data3d  )
       Item%KindVal    =  KIND( Data3d )
    ELSE IF ( PRESENT( Data2d  ) ) THEN
       Item%Rank       =  2
       Item%Ptr2d      => Data2d
       Item%MemoryInKb =  KbPerElement * SIZE( Data2d  )
       Item%KindVal    =  KIND( Data2d )
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

    !-----------------------------------------------------------------------
    ! Assign pointers to 4-byte real data targets
    !-----------------------------------------------------------------------
    ELSE IF ( PRESENT( Data3d_4 ) ) THEN
       Item%Rank       =  3
       Item%Ptr3d_4    => Data3d_4
       Item%MemoryInKb =  KbPerElement * SIZE( Data3d_4 )
       Item%KindVal    =  KIND( Data3d_4 )
    ELSE IF ( PRESENT( Data2d_4 ) ) THEN
       Item%Rank       =  2
       Item%Ptr2d_4    => Data2d_4
       Item%MemoryInKb =  KbPerElement * SIZE( Data2d_4 )
       Item%KindVal    =  KIND( Data2d_4 )
    ELSE IF ( PRESENT( Data1d_4 ) ) THEN
       Item%Rank       =  1
       Item%Ptr1d_4    => Data1d_4
       Item%MemoryInKb =  KbPerElement * SIZE( Data1d_4 )
       Item%KindVal    =  KIND( Data1d_4 )
    ELSE IF ( PRESENT( Data0d_4 ) ) THEN
       Item%Rank       =  0
       Item%Ptr0d_4    => Data0d_4
       Item%MemoryInKb =  KbPerElement
       Item%KindVal    =  KIND( Data0d_4 )

    !-----------------------------------------------------------------------
    ! Assign pointers to integer data targets
    !-----------------------------------------------------------------------
    ELSE IF ( PRESENT( Data3d_I ) ) THEN
       Item%Rank       =  3
       Item%Ptr3d_I    => Data3d_I
       Item%MemoryInKb =  KbPerElement * SIZE( Data3d_I )
       Item%KindVal    =  KIND( Data3d_I )
    ELSE IF ( PRESENT( Data2d_I ) ) THEN
       Item%Rank       =  2
       Item%Ptr2d_I    => Data2d_I
       Item%MemoryInKb =  KbPerElement * SIZE( Data2d_I )
       Item%KindVal    =  KIND( Data2d_I )
    ELSE IF ( PRESENT( Data1d_I  ) ) THEN
       Item%Rank       =  1
       Item%Ptr1d_I    => Data1d_I
       Item%MemoryInKb =  KbPerElement * SIZE( Data1d_I )
       Item%KindVal    =  KIND( Data1d_I )
    ELSE IF ( PRESENT( Data0d_I  ) ) THEN
       Item%Rank       =  0
       Item%Ptr0d_I    => Data0d_I
       Item%MemoryInKb =  KbPerElement
       Item%KindVal    =  KIND( Data0d_I )

    !-----------------------------------------------------------------------
    ! Exit with error message if no data target is passed
    !-----------------------------------------------------------------------
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
                              Ptr2d,     Ptr3d,       Ptr0d_4, Ptr1d_4,     &
                              Ptr2d_4,   Ptr3d_4,     Ptr0d_I, Ptr1d_I,     &
                              Ptr2d_I,   Ptr3d_I                           )
!
! !USES:
!
    USE ErrCode_Mod
!
! !INPUT PARAMETERS:
!
    LOGICAL,           INTENT(IN) :: am_I_Root         ! Root CPU?  
    TYPE(MetaRegItem),    POINTER :: Registry          ! Registry obj
    CHARACTER(LEN=*),  INTENT(IN) :: State             ! State name
    CHARACTER(LEN=*),  INTENT(IN) :: Variable          ! Variable name
!                                                     
! !OUTPUT PARAMETERS:                                 
!                                                     
    ! Required outputs
    INTEGER,          INTENT(OUT) :: RC                ! Success or failure

    ! Optional outputs
    CHARACTER(LEN=255),  OPTIONAL :: Description       ! Description of data
    INTEGER,             OPTIONAL :: KindVal           ! Numerical KIND value
    REAL(fp),            OPTIONAL :: MemoryInKb        ! Memory usage
    INTEGER,             OPTIONAL :: Rank              ! Size of data
    CHARACTER(LEN=255),  OPTIONAL :: Units             ! Units of data

    ! Floating-point data pointers (flex-precision)
    REAL(fp),   POINTER, OPTIONAL :: Ptr0d             ! 0D flex-prec data
    REAL(fp),   POINTER, OPTIONAL :: Ptr1d (:    )     ! 1D flex-prec data
    REAL(fp),   POINTER, OPTIONAL :: Ptr2d (:,:  )     ! 2D flex-prec data
    REAL(fp),   POINTER, OPTIONAL :: Ptr3d (:,:,:)     ! 3D flex-prec data

    ! Floating-point data pointers (4-byte precision)
    REAL(f4),   POINTER, OPTIONAL :: Ptr0d_4           ! 0D 4-byte data
    REAL(f4),   POINTER, OPTIONAL :: Ptr1d_4(:    )    ! 1D 4-byte data
    REAL(f4),   POINTER, OPTIONAL :: Ptr2d_4(:,:  )    ! 2D 4-byte data
    REAL(f4),   POINTER, OPTIONAL :: Ptr3d_4(:,:,:)    ! 3D 4-byte data

    ! Integer data pointers
    INTEGER,    POINTER, OPTIONAL :: Ptr0d_I           ! 0D integer data
    INTEGER,    POINTER, OPTIONAL :: Ptr1d_I(:    )    ! 1D integer data
    INTEGER,    POINTER, OPTIONAL :: Ptr2d_I(:,:  )    ! 2D integer data
    INTEGER,    POINTER, OPTIONAL :: Ptr3d_I(:,:,:)    ! 3D integer data
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
!  30 Jun 2017 - R. Yantosca - Added more pointers for 4-byte and integer data
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
    LOGICAL                    :: Is_0d,          Is_0d_4,       Is_0d_I
    LOGICAL                    :: Is_1d,          Is_1d_4,       Is_1d_I
    LOGICAL                    :: Is_2d,          Is_2d_4,       Is_2d_I
    LOGICAL                    :: Is_3d,          Is_3d_4,       Is_3d_I
    INTEGER                    :: FullHash,       ItemHash

    ! Strings
    CHARACTER(LEN=31)          :: FullName31,     ItemName31
    CHARACTER(LEN=255)         :: TmpName,        TmpFullName
    CHARACTER(LEN=255)         :: ErrMsg,         ThisLoc

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

    ! Floating-point (flex-precision) data pointers
    Is_0d          =  PRESENT( Ptr0d       )
    Is_1d          =  PRESENT( Ptr1d       )
    Is_2d          =  PRESENT( Ptr2d       )
    Is_3d          =  PRESENT( Ptr3d       )

    ! Floating-point (4-byte) data pointers
    Is_0d_4        =  PRESENT( Ptr0d_4     )
    Is_1d_4        =  PRESENT( Ptr1d_4     )
    Is_2d_4        =  PRESENT( Ptr2d_4     )
    Is_3d_4        =  PRESENT( Ptr3d_4     )

    ! Integer data pointers
    Is_0d_I        =  PRESENT( Ptr0d_I     )
    Is_1d_I        =  PRESENT( Ptr1d_I     )
    Is_2d_I        =  PRESENT( Ptr2d_I     )
    Is_3d_I        =  PRESENT( Ptr3d_I     )

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

             ! 3-d data
             CASE( 3 ) 
                IF ( Is_3d ) THEN
                   Ptr3d   => Current%Item%Ptr3d
                   RETURN
                ELSE IF ( Is_3d_4 ) THEN 
                   Ptr3d_4 => Current%Item%Ptr3d_4
                   RETURN                  
                ELSE IF  ( Is_3d_I ) THEN
                   Ptr3d_I => Current%Item%Ptr3d_I
                   RETURN
                ENDIF
                   
             ! 2-d data
             CASE( 2 )
                IF ( Is_2d ) THEN
                   Ptr2d   => Current%Item%Ptr2d
                   RETURN
                ELSE IF ( Is_2d_4 ) THEN
                   Ptr2d_4 => Current%Item%Ptr2d_4
                   RETURN
                ELSE IF ( Is_2d_I ) THEN
                   Ptr2d_I => Current%Item%Ptr2d_I
                   RETURN
                ENDIF

             ! 1-d data
             CASE( 1 )
                IF ( Is_1d ) THEN
                   Ptr1d   => Current%Item%Ptr1d
                   RETURN
                ELSE IF ( Is_1d_4 ) THEN
                   Ptr1d_4 => Current%Item%Ptr1d_4
                   RETURN                 
                ELSE IF ( Is_1d_I ) THEN
                   Ptr1d_I => Current%Item%Ptr1d_I
                   RETURN
                ENDIF
                
             ! 0-d data
             CASE( 0 )
                IF ( Is_0d ) THEN
                   Ptr0d   => Current%Item%Ptr0d
                   RETURN
                ELSE IF ( Is_0d_4 ) THEN
                   Ptr0d_4 => Current%Item%Ptr0d_4
                   RETURN
                ELSE IF ( Is_0d_I ) THEN
                   Ptr0d_I => Current%Item%Ptr0d_I
                   RETURN
                ENDIF

             ! Error message
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
!  27 Jun 2017 - R. Yantosca - Also print Ptr1dI, Ptr2dI, Ptr3dI 
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

             !--------------
             ! 3D data
             !--------------

             ! Flexible precision
             IF ( ASSOCIATED( Item%Ptr3d ) ) THEN
                PRINT*, 'Min value    : ', MINVAL( Item%Ptr3d      )
                PRINT*, 'Max value    : ', MAXVAL( Item%Ptr3d      )
                PRINT*, 'Total        : ', SUM   ( Item%Ptr3d      )
                PRINT*, 'Dimensions   : ', SIZE  ( Item%Ptr3d,   1 ), &
                                           SIZE  ( Item%Ptr3d,   2 ), &
                                           SIZE  ( Item%Ptr3d  , 3 )

             ! 4-byte
             ELSE IF ( ASSOCIATED( Item%Ptr3d_4 ) ) THEN
                PRINT*, 'Min value    : ', MINVAL( Item%Ptr3d_4    )
                PRINT*, 'Max value    : ', MAXVAL( Item%Ptr3d_4    )
                PRINT*, 'Total        : ', SUM   ( Item%Ptr3d_4    )
                PRINT*, 'Dimensions   : ', SIZE  ( Item%Ptr3d_4, 1 ), &
                                           SIZE  ( Item%Ptr3d_4, 2 ), &
                                           SIZE  ( Item%Ptr3d_4, 3 )

             ! Integer
             ELSE IF ( ASSOCIATED( Item%Ptr3d_I ) ) THEN
                PRINT*, 'Min value    : ', MINVAL( Item%Ptr3d_I    )
                PRINT*, 'Max value    : ', MAXVAL( Item%Ptr3d_I    )
                PRINT*, 'Total        : ', SUM   ( Item%Ptr3d_I    )
                PRINT*, 'Dimensions   : ', SIZE  ( Item%Ptr3d_I, 1 ), &
                                           SIZE  ( Item%Ptr3d_I, 2 ), &
                                           SIZE  ( Item%Ptr3d_I, 3 )
             !--------------
             ! 2D data
             !--------------

             ! Flexible precision
             ELSE IF ( ASSOCIATED( Item%Ptr2d ) ) THEN
                PRINT*, 'Min value    : ', MINVAL( Item%Ptr2d      )
                PRINT*, 'Max value    : ', MAXVAL( Item%Ptr2d      )
                PRINT*, 'Total        : ', SUM   ( Item%Ptr2d      )
                PRINT*, 'Dimensions   : ', SIZE  ( Item%Ptr2d, 1   ), &
                                           SIZE  ( Item%Ptr2d, 2   )
             ! 4-byte 
             ELSE IF ( ASSOCIATED( Item%Ptr2d_4 ) ) THEN
                PRINT*, 'Min value    : ', MINVAL( Item%Ptr2d_4    )
                PRINT*, 'Max value    : ', MAXVAL( Item%Ptr2d_4    )
                PRINT*, 'Total        : ', SUM   ( Item%Ptr2d_4    )
                PRINT*, 'Dimensions   : ', SIZE  ( Item%Ptr2d_4, 1 ), &
                                           SIZE  ( Item%Ptr2d_4, 2 )

             ! Integer
             ! 2D data -- Integer
             ELSE IF ( ASSOCIATED( Item%Ptr2d_I ) ) THEN
                PRINT*, 'Min value    : ', MINVAL( Item%Ptr2d_I    )
                PRINT*, 'Max value    : ', MAXVAL( Item%Ptr2d_I    )
                PRINT*, 'Total        : ', SUM   ( Item%Ptr2d_I    )
                PRINT*, 'Dimensions   : ', SIZE  ( Item%Ptr2d_I, 1 ), &
                                           SIZE  ( Item%Ptr2d_I, 2 )
             !--------------
             ! 1D data
             !--------------

             ! Flexible precision
             ELSE IF ( ASSOCIATED( Item%Ptr1d ) ) THEN
                PRINT*, 'Min value    : ', MINVAL( Item%Ptr1d      )
                PRINT*, 'Max value    : ', MAXVAL( Item%Ptr1d      )
                PRINT*, 'Total        : ', SUM   ( Item%Ptr1d      )
                PRINT*, 'Dimensions   : ', SIZE  ( Item%Ptr1d      )

             ! 4-byte
             ELSE IF ( ASSOCIATED( Item%Ptr1d_4 ) ) THEN
                PRINT*, 'Min value    : ', MINVAL( Item%Ptr1d_4    )
                PRINT*, 'Max value    : ', MAXVAL( Item%Ptr1d_4    )
                PRINT*, 'Total        : ', SUM   ( Item%Ptr1d_4    )
                PRINT*, 'Dimensions   : ', SIZE  ( Item%Ptr1d_4    )

             ! Integer
             ELSE IF ( ASSOCIATED( Item%Ptr1d_I ) ) THEN
                PRINT*, 'Min value    : ', MINVAL( Item%Ptr1d_I    )
                PRINT*, 'Max value    : ', MAXVAL( Item%Ptr1d_I    )
                PRINT*, 'Total        : ', SUM   ( Item%Ptr1d_I    )
                PRINT*, 'Dimensions   : ', SIZE  ( Item%Ptr1d_I    )

             !--------------
             ! 0D data
             !--------------

             ! Flexible precision
             ELSE IF ( ASSOCIATED( Item%Ptr0d ) ) THEN
                PRINT*, 'Value Ptr0d  : ', Item%Ptr0d

             ! 4-byte
             ELSE IF ( ASSOCIATED( Item%Ptr0d_4 ) ) THEN
                PRINT*, 'Value Ptr0d  : ', Item%Ptr0d_4

             ! Integer
             ELSE IF ( ASSOCIATED( Item%Ptr0d_I ) ) THEN
                PRINT*, 'Value Ptr0d  : ', Item%Ptr0d_I
  
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
!  29 Jun 2017 - R. Yantosca - Now nullify pointers to integer fields
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

       ! Free flexible-precision data pointers in this REGISTRY ITEM
       Current%Item%Ptr0d   => NULL()
       Current%Item%Ptr1d   => NULL()
       Current%Item%Ptr2d   => NULL()
       Current%Item%Ptr3d   => NULL()

       ! Free 4-byte data pointers in this REGISTRY ITEM
       Current%Item%Ptr0d_4 => NULL()
       Current%Item%Ptr1d_4 => NULL()
       Current%Item%Ptr2d_4 => NULL()
       Current%Item%Ptr3d_4 => NULL()

       ! Free integer data pointers in this REGISTRY ITEM
       Current%Item%Ptr0d_I => NULL()
       Current%Item%Ptr1d_I => NULL()
       Current%Item%Ptr2d_I => NULL()
       Current%Item%Ptr3d_I => NULL()
       
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
  FUNCTION To_UpperCase( Text ) RESULT( UpCaseText )
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
