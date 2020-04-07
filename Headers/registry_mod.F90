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
  USE Registry_Params_Mod

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC  :: Registry_AddField
  PUBLIC  :: Registry_Lookup
  PUBLIC  :: Registry_Print
  PUBLIC  :: Registry_Destroy
  PUBLIC  :: Registry_Set_LookupTable
!
! !PRIVATE MEMBER FUNCTIONS:
!
  PRIVATE :: MetaRegItem_AddNew
  PRIVATE :: MetaRegItem_Create
  PRIVATE :: MetaRegItem_Insert
  PRIVATE :: MetaRegItem_Destroy
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
     CHARACTER(LEN=3)     :: DimNames         ! e.g. "xyz", "yz", "y", "t"
     LOGICAL              :: OnLevelEdges     ! Is data on level edges (T/F)?

     !----------------------------------------------------------------------
     ! Pointers to floating point data (flexible precision)
     !----------------------------------------------------------------------
     REAL(fp), POINTER    :: Ptr0d            ! For 0D flex-prec data
     REAL(fp), POINTER    :: Ptr1d  (:    )   ! For 1D flex-prec data
     REAL(fp), POINTER    :: Ptr2d  (:,:  )   ! For 2D flex-prec data
     REAL(fp), POINTER    :: Ptr3d  (:,:,:)   ! For 3D flex-prec data

     !----------------------------------------------------------------------
     ! Pointers to floating point data (8-byte precision)
     !----------------------------------------------------------------------
     REAL(f8), POINTER    :: Ptr0d_8          ! For 0D 8-byte data
     REAL(f8), POINTER    :: Ptr1d_8(:    )   ! For 1D 8-byte data
     REAL(f8), POINTER    :: Ptr2d_8(:,:  )   ! For 2D 8-byte data
     REAL(f8), POINTER    :: Ptr3d_8(:,:,:)   ! For 3D 8-byte data

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
! !DEFINED PARAMETERS:
!
!
! !REMARKS:
!  In Fortran 2003, the maximum variable name length is 63 characers, so we
!  have declared the various character name fields of RegItem accordingly.
!
! !REVISION HISTORY:
!  23 Jun 2017 - R. Yantosca - Initial version
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
  SUBROUTINE Registry_AddField( Input_Opt, Registry,  State,                 &
                                Variable,  RC,        Description,           &
                                Units,     DimNames,  OnLevelEdges,          &
                                Data0d,    Data1d,    Data2d,                &
                                Data3d,    Data0d_8,  Data1d_8,              &
                                Data2d_8,  Data3d_8,  Data0d_4,              &
                                Data1d_4,  Data2d_4,  Data3d_4,              &
                                Data0d_I,  Data1d_I,  Data2d_I,              &
                                Data3d_I                                    )
!
! !USES:
!
    USE CharPak_Mod,   ONLY : To_Uppercase
    USE ErrCode_Mod
    USE Input_Opt_Mod, ONLY : OptInput
!
! !INPUT PARAMETERS:
!
    ! General identifying information
    TYPE(OptInput),    INTENT(IN)       :: Input_Opt       ! Input options
    CHARACTER(LEN=*),  INTENT(IN)       :: State           ! State name
    CHARACTER(LEN=*),  INTENT(IN)       :: Variable        ! variable
    CHARACTER(LEN=*),  OPTIONAL         :: Description     ! Long description
    CHARACTER(LEN=*),  OPTIONAL         :: Units           ! Units of data
    CHARACTER(LEN=*),  OPTIONAL         :: DimNames        ! "xyz", "xy", "t"
    LOGICAL,           OPTIONAL         :: OnLevelEdges    ! Set =T if data
                                                           !  is on level edges

    ! Floating-point data targets (flexible precision)
    REAL(fp),          OPTIONAL, TARGET :: Data0d          ! 0D flex-prec data
    REAL(fp),          OPTIONAL, TARGET :: Data1d  (:    ) ! 1D flex_prec data
    REAL(fp),          OPTIONAL, TARGET :: Data2d  (:,:  ) ! 2D flex-prec data
    REAL(fp),          OPTIONAL, TARGET :: Data3d  (:,:,:) ! 3D flex-prec data

    ! Floating-point data targets (8-byte precision)
    REAL(f8),          OPTIONAL, TARGET :: Data0d_8        ! 0D flex-prec data
    REAL(f8),          OPTIONAL, TARGET :: Data1d_8(:    ) ! 1D flex_prec data
    REAL(f8),          OPTIONAL, TARGET :: Data2d_8(:,:  ) ! 2D flex-prec data
    REAL(f8),          OPTIONAL, TARGET :: Data3d_8(:,:,:) ! 3D flex-prec data

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
!  the STATE and VARIABLE inputs as STATE_VARIABLE, unless variable is in
!  State_Diag, in which case STATE_ is not appended as a prefix.
!
! !REVISION HISTORY:
!  23 Jun 2017 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    LOGICAL                :: IsOnLevelEdges
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
    IsOnLevelEdges = .FALSE.
    KbPerElement   = DBLE( fp ) / 1024.0_fp
    ErrMsg         = ''
    ThisLoc        = ' -> at Registry_AddField (in Headers/registry_mod.F90)'
    TmpState       = To_UpperCase( State    )
    TmpVariable    = To_UpperCase( Variable )
    IF ( TRIM( TmpState ) /= 'DIAG' ) THEN
       TmpFullname    = TRIM( TmpState ) // '_' // TRIM( TmpVariable )
    ELSE
       TmpFullname = TRIM( TmpVariable )
    ENDIF
    TmpDescription = ''
    TmpUnits       = ''

    ! Save optional arguments in shadow variables, if passed
    IF ( PRESENT( Description  ) ) TmpDescription = Description
    IF ( PRESENT( Units        ) ) TmpUnits       = Units
    IF ( PRESENT( OnLevelEdges ) ) IsOnLevelEdges = OnLevelEdges

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
    Item%OnLevelEdges  =  IsOnLevelEdges

    ! Nullify pointers to data etc.
    Item%Ptr0d         => NULL()
    Item%Ptr1d         => NULL()
    Item%Ptr2d         => NULL()
    Item%Ptr3d         => NULL()
    Item%Ptr0d_8       => NULL()
    Item%Ptr1d_8       => NULL()
    Item%Ptr2d_8       => NULL()
    Item%Ptr3d_8       => NULL()
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
       Item%KindVal    =  KINDVAL_FP
    ELSE IF ( PRESENT( Data2d  ) ) THEN
       Item%Rank       =  2
       Item%Ptr2d      => Data2d
       Item%MemoryInKb =  KbPerElement * SIZE( Data2d  )
       Item%KindVal    =  KINDVAL_FP
    ELSE IF ( PRESENT( Data1d  ) ) THEN
       Item%Rank       =  1
       Item%Ptr1d      => Data1d
       Item%MemoryInKb =  KbPerElement * SIZE( Data1d  )
       Item%KindVal    =  KINDVAL_FP
    ELSE IF ( PRESENT( Data0d  ) ) THEN
       Item%Rank       =  0
       Item%Ptr0d      => Data0d
       Item%MemoryInKb =  KbPerElement
       Item%KindVal    =  KINDVAL_FP

    !-----------------------------------------------------------------------
    ! Assign pointers to 8-byte real data targets
    !-----------------------------------------------------------------------
    ELSE IF ( PRESENT( Data3d_8 ) ) THEN
       Item%Rank       =  3
       Item%Ptr3d_8    => Data3d_8
       Item%MemoryInKb =  KbPerElement * SIZE( Data3d_8 )
       Item%KindVal    =  KINDVAL_F8
    ELSE IF ( PRESENT( Data2d_8 ) ) THEN
       Item%Rank       =  2
       Item%Ptr2d_8    => Data2d_8
       Item%MemoryInKb =  KbPerElement * SIZE( Data2d_8  )
       Item%KindVal    =  KINDVAL_F8
    ELSE IF ( PRESENT( Data1d_8 ) ) THEN
       Item%Rank       =  1
       Item%Ptr1d_8    => Data1d_8
       Item%MemoryInKb =  KbPerElement * SIZE( Data1d_8  )
       Item%KindVal    =  KINDVAL_F8
    ELSE IF ( PRESENT( Data0d_8 ) ) THEN
       Item%Rank       =  0
       Item%Ptr0d_8    => Data0d_8
       Item%MemoryInKb =  KbPerElement
       Item%KindVal    =  KINDVAL_F8

    !-----------------------------------------------------------------------
    ! Assign pointers to 4-byte real data targets
    !-----------------------------------------------------------------------
    ELSE IF ( PRESENT( Data3d_4 ) ) THEN
       Item%Rank       =  3
       Item%Ptr3d_4    => Data3d_4
       Item%MemoryInKb =  KbPerElement * SIZE( Data3d_4 )
       Item%KindVal    =  KINDVAL_F4
    ELSE IF ( PRESENT( Data2d_4 ) ) THEN
       Item%Rank       =  2
       Item%Ptr2d_4    => Data2d_4
       Item%MemoryInKb =  KbPerElement * SIZE( Data2d_4 )
       Item%KindVal    =  KINDVAL_F4
    ELSE IF ( PRESENT( Data1d_4 ) ) THEN
       Item%Rank       =  1
       Item%Ptr1d_4    => Data1d_4
       Item%MemoryInKb =  KbPerElement * SIZE( Data1d_4 )
       Item%KindVal    =  KINDVAL_F4
    ELSE IF ( PRESENT( Data0d_4 ) ) THEN
       Item%Rank       =  0
       Item%Ptr0d_4    => Data0d_4
       Item%MemoryInKb =  KbPerElement
       Item%KindVal    =  KINDVAL_F4

    !-----------------------------------------------------------------------
    ! Assign pointers to integer data targets
    !-----------------------------------------------------------------------
    ELSE IF ( PRESENT( Data3d_I ) ) THEN
       Item%Rank       =  3
       Item%Ptr3d_I    => Data3d_I
       Item%MemoryInKb =  KbPerElement * SIZE( Data3d_I )
       Item%KindVal    =  KINDVAL_I4
    ELSE IF ( PRESENT( Data2d_I ) ) THEN
       Item%Rank       =  2
       Item%Ptr2d_I    => Data2d_I
       Item%MemoryInKb =  KbPerElement * SIZE( Data2d_I )
       Item%KindVal    =  KINDVAL_I4
    ELSE IF ( PRESENT( Data1d_I  ) ) THEN
       Item%Rank       =  1
       Item%Ptr1d_I    => Data1d_I
       Item%MemoryInKb =  KbPerElement * SIZE( Data1d_I )
       Item%KindVal    =  KINDVAL_I4
    ELSE IF ( PRESENT( Data0d_I  ) ) THEN
       Item%Rank       =  0
       Item%Ptr0d_I    => Data0d_I
       Item%MemoryInKb =  KbPerElement
       Item%KindVal    =  KINDVAL_I4

    !-----------------------------------------------------------------------
    ! Exit with error message if no data target is passed
    !-----------------------------------------------------------------------
    ELSE
       ErrMsg = 'Need to specify a data source!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    !=======================================================================
    ! Define the "dimnames" field
    !=======================================================================
    IF ( PRESENT( DimNames ) ) THEN

       ! If the DimNames argument is passed, then use it
       Item%DimNames = TRIM( DimNames )

    ELSE

       ! Otherwise, set default DimNames based on the rank
       SELECT CASE( Item%Rank )
          CASE( 3 )
             Item%DimNames = 'xyz'
          CASE( 2 )
             Item%DimNames = 'xy '
          CASE( 1 )
             Item%DimNames = 'x  '
          CASE( 0 )
             Item%DimNames = '-  '
        END SELECT

     ENDIF

    !=======================================================================
    ! Add the REGISTRY ITEM to the METAREGISTRY ITEM, which represents
    ! the list of all data fields contained in a module.
    !=======================================================================
    CALL MetaRegItem_AddNew( Input_Opt, Registry, Item, RC )
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
  SUBROUTINE Registry_Lookup( am_I_Root,    Registry,     RegDict,           &
                              State,        Variable,     RC,                &
                              Description,  Dimensions,   KindVal,           &
                              MemoryInKb,   OnLevelEdges, Rank,              &
                              Units,        DimNames,     Ptr0d,             &
                              Ptr1d,        Ptr2d,        Ptr3d,             &
                              Ptr0d_8,      Ptr1d_8,      Ptr2d_8,           &
                              Ptr3d_8,      Ptr0d_4,      Ptr1d_4,           &
                              Ptr2d_4,      Ptr3d_4,      Ptr0d_I,           &
                              Ptr1d_I,      Ptr2d_I,      Ptr3d_I           )
!
! !USES:
!
    USE Charpak_Mod,   ONLY : To_UpperCase
    USE Dictionary_M,  ONLY : dictionary_t
    USE ErrCode_Mod
!
! !INPUT PARAMETERS:
!
    LOGICAL,           INTENT(IN) :: am_I_Root         ! Root CPU?
    TYPE(MetaRegItem), POINTER    :: Registry          ! Registry obj
    TYPE(dictionary_t)            :: RegDict           ! Registry lookup table
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
    INTEGER,             OPTIONAL :: Dimensions(3)     ! Dimensions of data
    CHARACTER(LEN=255),  OPTIONAL :: Units             ! Units of data
    CHARACTER(LEN=3),    OPTIONAL :: DimNames          ! "xyz", "xz", "t" etc.
    LOGICAL,             OPTIONAL :: OnLevelEdges      ! Is the data defined
                                                       !  on level edges (T/F)

    ! Floating-point data pointers (flex-precision)
    REAL(fp),   POINTER, OPTIONAL :: Ptr0d             ! 0D flex-prec data
    REAL(fp),   POINTER, OPTIONAL :: Ptr1d  (:    )    ! 1D flex-prec data
    REAL(fp),   POINTER, OPTIONAL :: Ptr2d  (:,:  )    ! 2D flex-prec data
    REAL(fp),   POINTER, OPTIONAL :: Ptr3d  (:,:,:)    ! 3D flex-prec data

    ! Floating-point data pointers (4-byte precision)
    REAL(f8),   POINTER, OPTIONAL :: Ptr0d_8           ! 0D 8-byte data
    REAL(f8),   POINTER, OPTIONAL :: Ptr1d_8(:    )    ! 1D 8-byte data
    REAL(f8),   POINTER, OPTIONAL :: Ptr2d_8(:,:  )    ! 2D 8-byte data
    REAL(f8),   POINTER, OPTIONAL :: Ptr3d_8(:,:,:)    ! 3D 8-byte data

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
!  the STATE and VARIABLE inputs as STATE_VARIABLE, unless variable is in
!  State_Diag, in which case STATE_ is not appended as a prefix.
!
! !REVISION HISTORY:
!  23 Jun 2017 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    LOGICAL                    :: Is_OnLevelEdges, Found
    LOGICAL                    :: Is_Description,  Is_Dimensions, Is_KindVal
    LOGICAL                    :: Is_MemoryInKb,   Is_Rank,       Is_Units
    LOGICAL                    :: Is_0d,           Is_0d_8
    LOGICAL                    :: Is_0d_4,         Is_0d_I
    LOGICAL                    :: Is_1d,           Is_1d_8
    LOGICAL                    :: Is_1d_4,         Is_1d_I
    LOGICAL                    :: Is_2d,           Is_2d_8
    LOGICAL                    :: Is_2d_4,         Is_2d_I
    LOGICAL                    :: Is_3d,           Is_3d_8
    LOGICAL                    :: Is_3d_4,         Is_3d_I
    INTEGER                    :: FullHash,        ItemHash,      N

    ! Strings
    CHARACTER(LEN=5)           :: TmpState
    CHARACTER(LEN=67)          :: FullName,        ItemName
    CHARACTER(LEN=255)         :: ErrMsg,          ThisLoc
    CHARACTER(LEN=255)         :: VariableUC

    ! Objects
    TYPE(MetaRegItem), POINTER :: Current

    !=======================================================================
    ! Initialize
    !=======================================================================
    RC              =  GC_SUCCESS
    Current         => NULL()
    ErrMsg          =  ''
    ThisLoc         =  ' -> at Registry_Lookup (in Headers/registry_mod.F90)'
    TmpState        = TRIM( State ) // '_'
    VariableUC      =  To_UpperCase( Variable )

    ! Prefix the the state name (always uppercase) to the variable, unless:
    ! (1) The state name is already part of the variable
    ! (2) If it's a field from State_Diag, which requires no prefix
    IF ( ( TRIM( State ) == 'DIAG' ) .OR.  &
         ( INDEX( VariableUC, TRIM( TmpState ) ) > 0 ) ) THEN
       FullName  = VariableUC
    ELSE
       FullName  = TRIM( TmpState ) // TRIM( VariableUC )
    ENDIF

    ! Construct a hash for the full name (i.e. "State_Variable")
    FullHash =  RegDict%Get( TRIM( FullName ) )

    ! Return with an error if fullname is not found in this registry
    IF ( FullHash == -1 ) THEN
       errMsg = TRIM( fullName ) // ' is not found in the registry for '  // &
                'the ' // TRIM( State ) // ' object!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    ! Set a flag to denote that we've found the field
    Found = .FALSE.

    !=======================================================================
    ! Test if the optional variables are present outside of the main loop.
    !=======================================================================

    ! Floating-point (flex-precision) data pointers
    Is_0d           =  PRESENT( Ptr0d        )
    Is_1d           =  PRESENT( Ptr1d        )
    Is_2d           =  PRESENT( Ptr2d        )
    Is_3d           =  PRESENT( Ptr3d        )

    ! Floating-point (8-byte) data pointers
    Is_0d_8         =  PRESENT( Ptr0d_8      )
    Is_1d_8         =  PRESENT( Ptr1d_8      )
    Is_2d_8         =  PRESENT( Ptr2d_8      )
    Is_3d_8         =  PRESENT( Ptr3d_8      )

    ! Floating-point (4-byte) data pointers
    Is_0d_4         =  PRESENT( Ptr0d_4      )
    Is_1d_4         =  PRESENT( Ptr1d_4      )
    Is_2d_4         =  PRESENT( Ptr2d_4      )
    Is_3d_4         =  PRESENT( Ptr3d_4      )

    ! Integer data pointers
    Is_0d_I         =  PRESENT( Ptr0d_I      )
    Is_1d_I         =  PRESENT( Ptr1d_I      )
    Is_2d_I         =  PRESENT( Ptr2d_I      )
    Is_3d_I         =  PRESENT( Ptr3d_I      )

    ! Metadata
    Is_Description  =  PRESENT( Description  )
    Is_Dimensions   =  PRESENT( Dimensions   )
    Is_KindVal      =  PRESENT( KindVal      )
    Is_MemoryInKb   =  PRESENT( MemoryInKb   )
    Is_Rank         =  PRESENT( Rank         )
    Is_Units        =  PRESENT( Units        )
    Is_OnLevelEdges =  PRESENT( OnLevelEdges )

    !=======================================================================
    ! Nullify all optional pointer arguments that are passed
    !=======================================================================
    If ( Is_0d   ) Ptr0d   => NULL()
    IF ( Is_0d_8 ) Ptr0d_8 => NULL()
    IF ( Is_0d_4 ) Ptr0d_4 => NULL()
    IF ( Is_0d_I ) Ptr0d_I => NULL()
    IF ( Is_1d   ) Ptr1d   => NULL()
    IF ( Is_1d_8 ) Ptr1d_8 => NULL()
    IF ( Is_1d_4 ) Ptr1d_4 => NULL()
    IF ( Is_1d_I ) Ptr1d_I => NULL()
    IF ( Is_2d   ) Ptr2d   => NULL()
    IF ( Is_2d_8 ) Ptr2d_8 => NULL()
    IF ( Is_2d_4 ) Ptr2d_4 => NULL()
    IF ( Is_2d_I ) Ptr2d_I => NULL()
    IF ( Is_3d   ) Ptr3d   => NULL()
    IF ( Is_3d_8 ) Ptr3d_8 => NULL()
    IF ( Is_3d_4 ) Ptr3d_4 => NULL()
    IF ( Is_3d_I ) Ptr3d_I => NULL()

    !=======================================================================
    ! Search for the specified field in the Registry
    !=======================================================================

    ! Point to head of linked list
    Current => Registry

    ! As long as this entry of the linked list isn't NULL
    DO WHILE( ASSOCIATED( Current ) )

       ! Construct a hash for the full name of this REGISTRY ITEM
       ItemHash   = RegDict%Get( TRIM( Current%Item%FullName ) )

       ! If the name-hashes match (and are not missing data "-1")
       IF ( FullHash == ItemHash .and. ItemHash > 0 ) THEN

          ! Return rank, units and memory usage, etc. if found
          IF ( Is_Description  ) Description  = Current%Item%Description
          IF ( Is_KindVal      ) KindVal      = Current%Item%KindVal
          IF ( Is_MemoryInKb   ) MemoryInKb   = Current%Item%MemoryInKb
          IF ( Is_Rank         ) Rank         = Current%Item%Rank
          IF ( Is_Units        ) Units        = Current%Item%Units
          If ( Is_OnLevelEdges ) OnLevelEdges = Current%Item%OnLevelEdges

          ! Then return a pointer to the field
          SELECT CASE( Current%Item%Rank )

             ! Return the appropriate 3D DATA POINTER (and dimensions)
             CASE( 3 )
                IF ( Current%Item%KindVal == KINDVAL_FP ) THEN
                   IF ( Is_3d ) THEN
                      Ptr3d => Current%Item%Ptr3d
                      Found =  .TRUE.
                      IF ( Is_Dimensions ) THEN
                         DO N = 1, Current%Item%Rank
                            Dimensions(N) = SIZE( Ptr3d, N )
                         ENDDO
                      ENDIF
                   ENDIF
                   EXIT
                ELSE IF ( Current%Item%KindVal == KINDVAL_F8 ) THEN
                   IF ( Is_3d_8 ) THEN
                      Ptr3d_8 => Current%Item%Ptr3d_8
                      Found   =  .TRUE.
                      IF ( Is_Dimensions ) THEN
                         DO N = 1, Current%Item%Rank
                            Dimensions(N) = SIZE( Ptr3d_8, N )
                         ENDDO
                      ENDIF
                   ENDIF
                   EXIT
                ELSE IF ( Current%Item%KindVal == KINDVAL_F4 ) THEN
                   IF ( Is_3d_4 ) THEN
                      Ptr3d_4 => Current%Item%Ptr3d_4
                      Found   =  .TRUE.
                      IF ( Is_Dimensions ) THEN
                         DO N = 1, Current%Item%Rank
                            Dimensions(N) = SIZE( Ptr3d_4, N )
                         ENDDO
                      ENDIF
                   ENDIF
                   EXIT
                ELSE IF ( Current%Item%KindVal == KINDVAL_I4 ) THEN
                   IF ( Is_3d_I ) THEN
                      Ptr3d_I => Current%Item%Ptr3d_I
                      Found   =  .TRUE.
                      IF ( Is_Dimensions ) THEN
                         DO N = 1, Current%Item%Rank
                            Dimensions(N) = SIZE( Ptr3d_I, N )
                         ENDDO
                      ENDIF
                   ENDIF
                   EXIT
                ENDIF

             ! Return the appropriate 2D DATA POINTER (and dimensions)
             CASE( 2 )
                IF ( Current%Item%KindVal == KINDVAL_FP ) THEN
                   IF ( Is_2d ) THEN
                      Ptr2d => Current%Item%Ptr2d
                      Found =  .TRUE.
                      IF ( Is_Dimensions ) THEN
                         DO N = 1, Current%Item%Rank
                            Dimensions(N) = SIZE( Ptr2d, N )
                         ENDDO
                      ENDIF
                   ENDIF
                   EXIT
                ELSE IF ( Current%Item%KindVal == KINDVAL_F8 ) THEN
                   IF ( Is_2d_8 ) THEN
                      Ptr2d_8 => Current%Item%Ptr2d_8
                      Found   =  .TRUE.
                      IF ( Is_Dimensions ) THEN
                         DO N = 1, Current%Item%Rank
                            Dimensions(N) = SIZE( Ptr2d_8, N )
                         ENDDO
                      ENDIF
                   ENDIF
                   EXIT
                ELSE IF ( Current%Item%KindVal == KINDVAL_F4 ) THEN
                   IF ( Is_2d_4 ) THEN
                      Ptr2d_4 => Current%Item%Ptr2d_4
                      Found   =  .TRUE.
                      IF ( Is_Dimensions ) THEN
                         DO N = 1, Current%Item%Rank
                            Dimensions(N) = SIZE( Ptr2d_4, N )
                         ENDDO
                      ENDIF
                   ENDIF
                   EXIT
                ELSE IF ( Current%Item%KindVal == KINDVAL_I4 ) THEN
                   IF ( Is_2d_I ) THEN
                      Ptr2d_I => Current%Item%Ptr2d_I
                      Found =    .TRUE.
                      IF ( Is_Dimensions ) THEN
                         DO N = 1, Current%Item%Rank
                            Dimensions(N) = SIZE( Ptr2d_I, N )
                         ENDDO
                      ENDIF
                   ENDIF
                   EXIT
                ENDIF

             ! Return the appropriate 1D DATA POINTER (and dimensions)
             CASE( 1 )
                IF ( Current%Item%KindVal == KINDVAL_FP ) THEN
                   IF ( Is_1d ) THEN
                      Ptr1d => Current%Item%Ptr1d
                      Found =  .TRUE.
                      IF ( Is_Dimensions ) THEN
                         DO N = 1, Current%Item%Rank
                            Dimensions(N) = SIZE( Ptr1d, N )
                         ENDDO
                      ENDIF
                   ENDIF
                   EXIT
                ELSE IF ( Current%Item%KindVal == KINDVAL_F8 ) THEN
                   IF ( Is_1d_8 ) THEN
                      Ptr1d_8 => Current%Item%Ptr1d_8
                      Found   =  .TRUE.
                      IF ( Is_Dimensions ) THEN
                         DO N = 1, Current%Item%Rank
                            Dimensions(N) = SIZE( Ptr1d_8, N )
                         ENDDO
                      ENDIF
                   ENDIF
                   EXIT
                ELSE IF ( Current%Item%KindVal == KINDVAL_F4 ) THEN
                   IF ( Is_1d_4 ) THEN
                      Ptr1d_4 => Current%Item%Ptr1d_4
                      Found   =  .TRUE.
                      IF ( Is_Dimensions ) THEN
                         DO N = 1, Current%Item%Rank
                            Dimensions(N) = SIZE( Ptr1d_4, N )
                         ENDDO
                      ENDIF
                   ENDIF
                   EXIT
                ELSE IF ( Current%Item%KindVal == KINDVAL_I4 ) THEN
                   IF ( Is_1d_I ) THEN
                      Ptr1d_I => Current%Item%Ptr1d_I
                      Found   =  .TRUE.
                      IF ( Is_Dimensions ) THEN
                         DO N = 1, Current%Item%Rank
                            Dimensions(N) = SIZE( Ptr1d_I, N )
                         ENDDO
                      ENDIF
                   ENDIF
                   EXIT
                ENDIF

             ! Return the appropriate 0D DATA POINTER (and dimensions)
             CASE( 0 )
                IF ( Current%Item%KindVal == KINDVAL_FP ) THEN
                   IF ( Is_0d ) THEN
                      Ptr0d => Current%Item%Ptr0d
                      Found =   .TRUE.
                      IF ( Is_Dimensions ) Dimensions = 0
                   ENDIF
                   EXIT
                ELSE IF ( Current%Item%KindVal == KINDVAL_F8 ) THEN
                   IF ( Is_0d_8 ) THEN
                      Ptr0d_8 => Current%Item%Ptr0d_8
                      Found   =  .TRUE.
                      IF ( Is_Dimensions ) Dimensions = 0
                   ENDIF
                   EXIT
                ELSE IF ( Current%Item%KindVal == KINDVAL_F4 ) THEN
                   IF ( Is_0d_4 ) THEN
                      Ptr0d_4 => Current%Item%Ptr0d_4
                      Found   =  .TRUE.
                      IF ( Is_Dimensions ) Dimensions = 0
                   ENDIF
                   EXIT
                ELSE IF ( Current%Item%KindVal == KINDVAL_I4 ) THEN
                   IF ( Is_0d_I ) THEN
                      Ptr0d_I => Current%Item%Ptr0d_I
                      Found   =  .TRUE.
                      IF ( Is_Dimensions ) Dimensions = 0
                   ENDIF
                   EXIT
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
999 CONTINUE

    ! Free pointer
    Current => NULL()

    ! Throw an error if not found
    IF ( .not. Found ) RC = GC_FAILURE

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
  SUBROUTINE Registry_Print( Input_Opt, Registry, RC, ShortFormat )
!
! !USES:
!
    USE Input_Opt_Mod,      ONLY : OptInput
    USE ErrCode_Mod
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput),    INTENT(IN)  :: Input_Opt     ! Input Options object
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
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    LOGICAL                    :: Use_ShortFormat

    ! Strings
    CHARACTER(LEN=1)           :: CellPos
    CHARACTER(LEN=255)         :: ErrMsg,  ThisLoc

    ! Objects
    TYPE(MetaRegItem), POINTER :: Current
    TYPE(RegItem    ), POINTER :: Item

    !=======================================================================
    ! Initialize
    !=======================================================================

    ! Only print information on the root CPU
    IF ( .not. Input_Opt%amIRoot ) RETURN

    ! Initialize fields
    RC      =  GC_SUCCESS
    CellPos =  ''
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

             ! Denote if 3-D data is defined on level edges (E) or centers (C)
             IF ( Item%Rank == 3 ) THEN
                IF ( Item%OnLevelEdges ) THEN
                   CellPos = 'E'
                ELSE
                   CellPos = 'C'
                ENDIF
             ELSE
                CellPos = ''
             ENDIF

             ! Print information
             WRITE( 6, 100 ) Item%FullName,    Item%Description,             &
                             Item%DimNames,    CellPos,                      &
                             TRIM( Item%Units )
  100        FORMAT( 1x, a30, ' | ', a20, ' | ', a3, ' ', a1, ' | ', a )

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
             PRINT*, 'Dim Names    : ', TRIM( Item%DimNames    )
             PRINT*, 'KIND value   : ', Item%KindVal
             PRINT*, 'Memory (Kb)  : ', Item%MemoryInKb
             PRINT*, 'Rank of data : ', Item%Rank, '(', Item%DimNames, ')'
             PRINT*, 'On Edges?    : ', Item%OnLevelEdges

             !--------------
             ! 3D data
             !--------------

             ! Flexible precision
             IF ( ASSOCIATED( Item%Ptr3d ) ) THEN
                PRINT*, 'Min value    : ', MINVAL( Item%Ptr3d      )
                PRINT*, 'Max value    : ', MAXVAL( Item%Ptr3d      )
                PRINT*, 'Total        : ', SUM   ( Item%Ptr3d      )
                PRINT*, 'Dimensions   : ', SIZE  ( Item%Ptr3d,   1 ),        &
                                           SIZE  ( Item%Ptr3d,   2 ),        &
                                           SIZE  ( Item%Ptr3d  , 3 )

             ! 8-byte
             ELSE IF ( ASSOCIATED( Item%Ptr3d_8 ) ) THEN
                PRINT*, 'Min value    : ', MINVAL( Item%Ptr3d_8    )
                PRINT*, 'Max value    : ', MAXVAL( Item%Ptr3d_8    )
                PRINT*, 'Total        : ', SUM   ( Item%Ptr3d_8    )
                PRINT*, 'Dimensions   : ', SIZE  ( Item%Ptr3d_8, 1 ),        &
                                           SIZE  ( Item%Ptr3d_8, 2 ),        &
                                           SIZE  ( Item%Ptr3d_8, 3 )

             ! 4-byte
             ELSE IF ( ASSOCIATED( Item%Ptr3d_4 ) ) THEN
                PRINT*, 'Min value    : ', MINVAL( Item%Ptr3d_4    )
                PRINT*, 'Max value    : ', MAXVAL( Item%Ptr3d_4    )
                PRINT*, 'Total        : ', SUM   ( Item%Ptr3d_4    )
                PRINT*, 'Dimensions   : ', SIZE  ( Item%Ptr3d_4, 1 ),        &
                                           SIZE  ( Item%Ptr3d_4, 2 ),        &
                                           SIZE  ( Item%Ptr3d_4, 3 )

             ! Integer
             ELSE IF ( ASSOCIATED( Item%Ptr3d_I ) ) THEN
                PRINT*, 'Min value    : ', MINVAL( Item%Ptr3d_I    )
                PRINT*, 'Max value    : ', MAXVAL( Item%Ptr3d_I    )
                PRINT*, 'Total        : ', SUM   ( Item%Ptr3d_I    )
                PRINT*, 'Dimensions   : ', SIZE  ( Item%Ptr3d_I, 1 ),        &
                                           SIZE  ( Item%Ptr3d_I, 2 ),        &
                                           SIZE  ( Item%Ptr3d_I, 3 )
             !--------------
             ! 2D data
             !--------------

             ! Flexible precision
             ELSE IF ( ASSOCIATED( Item%Ptr2d ) ) THEN
                PRINT*, 'Min value    : ', MINVAL( Item%Ptr2d      )
                PRINT*, 'Max value    : ', MAXVAL( Item%Ptr2d      )
                PRINT*, 'Total        : ', SUM   ( Item%Ptr2d      )
                PRINT*, 'Dimensions   : ', SIZE  ( Item%Ptr2d, 1   ),        &
                                           SIZE  ( Item%Ptr2d, 2   )

             ! 8-byte
             ELSE IF ( ASSOCIATED( Item%Ptr2d_8 ) ) THEN
                PRINT*, 'Min value    : ', MINVAL( Item%Ptr2d_8    )
                PRINT*, 'Max value    : ', MAXVAL( Item%Ptr2d_8    )
                PRINT*, 'Total        : ', SUM   ( Item%Ptr2d_8    )
                PRINT*, 'Dimensions   : ', SIZE  ( Item%Ptr2d_8, 1 ),        &
                                           SIZE  ( Item%Ptr2d_8, 2 )

             ! 4-byte
             ELSE IF ( ASSOCIATED( Item%Ptr2d_4 ) ) THEN
                PRINT*, 'Min value    : ', MINVAL( Item%Ptr2d_4    )
                PRINT*, 'Max value    : ', MAXVAL( Item%Ptr2d_4    )
                PRINT*, 'Total        : ', SUM   ( Item%Ptr2d_4    )
                PRINT*, 'Dimensions   : ', SIZE  ( Item%Ptr2d_4, 1 ),        &
                                           SIZE  ( Item%Ptr2d_4, 2 )

             ! Integer
             ELSE IF ( ASSOCIATED( Item%Ptr2d_I ) ) THEN
                PRINT*, 'Min value    : ', MINVAL( Item%Ptr2d_I    )
                PRINT*, 'Max value    : ', MAXVAL( Item%Ptr2d_I    )
                PRINT*, 'Total        : ', SUM   ( Item%Ptr2d_I    )
                PRINT*, 'Dimensions   : ', SIZE  ( Item%Ptr2d_I, 1 ),        &
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

             ! 8-byte
             ELSE IF ( ASSOCIATED( Item%Ptr1d_8 ) ) THEN
                PRINT*, 'Min value    : ', MINVAL( Item%Ptr1d_8    )
                PRINT*, 'Max value    : ', MAXVAL( Item%Ptr1d_8    )
                PRINT*, 'Total        : ', SUM   ( Item%Ptr1d_8    )
                PRINT*, 'Dimensions   : ', SIZE  ( Item%Ptr1d_8    )

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
                PRINT*, 'Value        : ', Item%Ptr0d

             ! 8-byte precision
             ELSE IF ( ASSOCIATED( Item%Ptr0d_8 ) ) THEN
                PRINT*, 'Value        : ', Item%Ptr0d_8

             ! 4-byte
             ELSE IF ( ASSOCIATED( Item%Ptr0d_4 ) ) THEN
                PRINT*, 'Value        : ', Item%Ptr0d_4

             ! Integer
             ELSE IF ( ASSOCIATED( Item%Ptr0d_I ) ) THEN
                PRINT*, 'Value        : ', Item%Ptr0d_I

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
! !IROUTINE: Registry_Set_LookupTable
!
! !DESCRIPTION: Defines the lookup table for registry items, using the
!  dictionary_m algorithm.  This will avoid hash collisions.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Registry_Set_LookupTable( Registry, RegDict, RC )
!
! !USES:
!
    USE Dictionary_M
    USE ErrCode_Mod
    USE Input_Opt_Mod, ONLY : OptInput
!
! !INPUT PARAMETERS:
!
    TYPE(MetaRegItem),  POINTER       :: Registry   ! Registry of state fields
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(dictionary_t), INTENT(INOUT) :: RegDict    ! Registry lookup table
!
! !OUTPUT PARAMETERS:
!
    INTEGER,            INTENT(OUT)   :: RC         ! Success or failure?!
!
! !REVISION HISTORY:
!  07 Apr 2020 - R. Yantosca - Initial version
!  See the subsequent Git history with the gitk browser!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER                    :: index
    INTEGER                    :: nDiags


    ! Strings
    CHARACTER(LEN=255)         :: errMsg
    CHARACTER(LEN=255)         :: thisLoc

    ! Objects
    TYPE(MetaRegItem), POINTER :: current

    !=======================================================================
    ! Registry_Set_LookupTable begins here!
    !=======================================================================

    ! Initialize
    RC      = GC_SUCCESS
    index   = 0
    nDiags  = 0
    errMsg  = ''
    thisLoc = &
     ' -> at Registry_Set_LookupTable (in module Headers/registry_mod.F90)'

    !-----------------------------------------------------------------------
    ! First find out how many diagnostics have been registered
    ! so that we can initialize the lookup table accordingly.
    !-----------------------------------------------------------------------
    current => Registry
    DO WHILE( ASSOCIATED( current ) )
       IF ( ASSOCIATED( current%Item ) ) THEN
          nDiags = nDiags + 1
       ENDIF
       current => current%next
    ENDDO
    current => NULL()
    
    ! Initialize the lookup table
    CALL RegDict%Init( nDiags )

    !-----------------------------------------------------------------------
    ! Then populate the lookup table with the index with which each
    ! diagnostic is found in the list.  NOTE: Registry names are
    ! already uppercase, so no need to convert again.
    !-----------------------------------------------------------------------
    current => Registry
    DO WHILE( ASSOCIATED( current ) )
       IF ( ASSOCIATED( current%Item ) ) THEN
          index = index + 1
          CALL RegDict%Set( TRIM( current%Item%fullName ), index )
       ENDIF
       current => current%next
    ENDDO
    current => NULL()

  END SUBROUTINE Registry_Set_LookupTable
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
  SUBROUTINE Registry_Destroy( Registry, RegDict, RC )
!
! !USES:
!
    USE Dictionary_M
    USE ErrCode_Mod
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(MetaRegItem), POINTER     :: Registry   ! Registry of state fields
    TYPE(dictionary_t)             :: RegDict    ! Registry lookup table
!
! !OUTPUT PARAMETERS:
!
    INTEGER,           INTENT(OUT) :: RC         ! Success or failure?
!
! !REVISION HISTORY:
!  23 Jun 2017 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
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
    CALL MetaRegItem_Destroy( Registry, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Could not destroy the registry object!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    !=======================================================================
    ! Destroy the lookup table for this registry
    !=======================================================================
    CALL RegDict%Destroy()

  END SUBROUTINE Registry_Destroy
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: MetaRegItem_AddNew
!
! !DESCRIPTION: Wrapper for methods MetaRegItem\_Create and
!  MetaRegItem\_Insert.  Will create a METAREGISTRY ITEM (containing a
!  REGISTRY ITEM) and (1) set it as the head node of a new linked list, or
!  (2) append it to an existing linked list.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE MetaRegItem_AddNew( Input_Opt, Node, Item, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE Input_Opt_Mod, ONLY : OptInput
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput),    INTENT(IN)  :: Input_Opt  ! Input Options object
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
!  See https://github.com/geoschem/geos-chem for complete history
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
       CALL MetaRegItem_Create( Input_Opt, Node, Item, RC )
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
       CALL MetaRegItem_Insert( Input_Opt, Node, Item, RC )
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
  SUBROUTINE MetaRegItem_Create( Input_Opt, Node, Item, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE Input_Opt_Mod, ONLY : OptInput
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput),    INTENT(IN)  :: Input_Opt  ! Input Options object
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
!  See https://github.com/geoschem/geos-chem for complete history
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
  SUBROUTINE MetaRegItem_Insert( Input_Opt, Node, Item, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE Input_Opt_Mod, ONLY : OptInput
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput),    INTENT(IN)  :: Input_Opt ! Input Options object
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
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Strings
    CHARACTER(LEN=255)         :: ErrMsg, ThisLoc

    ! Objects
    TYPE(MetaRegItem), POINTER :: Head

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

    ! Allocate the "Head" object
    ALLOCATE( Head, STAT=RC )
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Could not allocate "Next"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    ! Allocate the "Head%Item" field, which will hold the REGISTRY ITEM
    ALLOCATE( Head%Item, STAT=RC )
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Could not allocate "Head%Item"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    !=======================================================================
    ! Insert "Head" at the start of the existing linked list
    !=======================================================================

    ! Save the REGISTRY ITEM argument in the "Item" field of "Head"
    Head%Item  =  Item

    ! The "Next" field of "Head" points to the current head of the list
    Head%Next  => Node

    ! Set "Head" as the new head of the linked list
    Node       => Head

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
  SUBROUTINE MetaRegItem_Destroy( List, RC )
!
! !USES:
!
    USE ErrCode_Mod
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
!  See https://github.com/geoschem/geos-chem for complete history
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

       ! Free 8-byte data pointers in this REGISTRY ITEM
       Current%Item%Ptr0d_8 => NULL()
       Current%Item%Ptr1d_8 => NULL()
       Current%Item%Ptr2d_8 => NULL()
       Current%Item%Ptr3d_8 => NULL()

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
#if defined( ESMF_ )
       IF ( ASSOCIATED( Current%Item ) ) NULLIFY( Current%Item )
#else
       IF ( ASSOCIATED( Current%Item ) ) DEALLOCATE( Current%Item )
#endif

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
END MODULE Registry_Mod
