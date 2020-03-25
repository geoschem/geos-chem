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
     INTEGER            :: NcXDimId              ! Id of netCDF X (lon  ) dim
     INTEGER            :: NcYdimId              ! Id of netCDF Y (lat  ) dim
     INTEGER            :: NcZDimId              ! Id of netCDF Z (lev C) dim
     INTEGER            :: NcIDimId              ! ID of netCDF I (lev E) dim
     INTEGER            :: NcTdimId              ! Id of netCDF T (time ) dim
     INTEGER            :: NcVarId               ! netCDF variable ID
     CHARACTER(LEN=255) :: LongName              ! Item description
     CHARACTER(LEN=255) :: Units                 ! Units of data
     REAL(f4)           :: AddOffset             ! Offset and scale factor
     REAL(f4)           :: ScaleFactor           !  for packed data
     REAL(f4)           :: MissingValue          ! Missing value
     CHARACTER(LEN=255) :: AvgMethod             ! Averaging method

     !----------------------------------------------------------------------
     ! Pointers to the data in State_Chm, State_Diag, or State_Met
     !----------------------------------------------------------------------
     INTEGER            :: Source_KindVal        ! Identifies the source type

     REAL(f8), POINTER  :: Source_0d_8           ! Ptr to 0D 8-byte    data

     REAL(fp), POINTER  :: Source_1d  (:    )    ! Ptr to 1D flex-prec data
     REAL(f8), POINTER  :: Source_1d_8(:    )    ! Ptr to 1D 8-byte    data
     REAL(f4), POINTER  :: Source_1d_4(:    )    ! Ptr to 1D 4-byte    data
     INTEGER,  POINTER  :: Source_1d_I(:    )    ! Ptr to 1D integer   data

     REAL(fp), POINTER  :: Source_2d  (:,:  )    ! Ptr to 2D flex-prec data
     REAL(f8), POINTER  :: Source_2d_8(:,:  )    ! Ptr to 2D 8-byte    data
     REAL(f4), POINTER  :: Source_2d_4(:,:  )    ! Ptr to 2D 4-byte    data
     INTEGER,  POINTER  :: Source_2d_I(:,:  )    ! Ptr to 2D integer   data

     REAL(fp), POINTER  :: Source_3d  (:,:,:)    ! Ptr to 3D flex-prec data
     REAL(f8), POINTER  :: Source_3d_8(:,:,:)    ! Ptr to 3D 8-byte    data
     REAL(f4), POINTER  :: Source_3d_4(:,:,:)    ! Ptr to 3D 4-byte    data
     INTEGER,  POINTER  :: Source_3d_I(:,:,:)    ! Ptr to 3D integer   data

     !----------------------------------------------------------------------
     ! Data arrays
     !----------------------------------------------------------------------
     INTEGER            :: SpaceDim              ! # of dims (0-3)
     REAL(f8), POINTER  :: Data_0d               ! 0D scalar
     REAL(f8), POINTER  :: Data_1d(:    )        ! 1D vector
     REAL(f8), POINTER  :: Data_2d(:,:  )        ! 2D array
     REAL(f8), POINTER  :: Data_3d(:,:,:)        ! 3D array
     CHARACTER(LEN=3)   :: DimNames              ! Used to specify if data is
                                                 !  "xyz", "yz", "x", "y" etc.
     INTEGER, POINTER   :: NcChunkSizes(:)       ! Chunk sizes for netCDF
     LOGICAL            :: OnLevelEdges          ! =T if data is defined on
                                                 !    vertical level edges;
                                                 ! =F if on level centers

     !----------------------------------------------------------------------
     ! Data archival
     !----------------------------------------------------------------------
     REAL(f8)           :: nUpdates              ! # of times updated
     INTEGER            :: Operation             ! Operation code
                                                 !  0=copy from source
                                                 !  1=accumulate from source
   END TYPE HistItem
!
! !REMARKS:
!  Linked list routines taken from original code (linkedlist.f90)
!  by Arjen Markus; http://flibs.sourceforge.net/linked_list.html
!
! !REVISION HISTORY:
!  13 Jun 2017 - R. Yantosca - Initial version, based on code by Arjen Markus
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
! !IROUTINE: HistItem_Create
!
! !DESCRIPTION: Initializes a single history item that will be archived
!  via History (and eventually sent to netCDF output).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HistItem_Create( Input_Opt,    Item,           Id,              &
                              ContainerId,  Name,           RC,              &
                              LongName,     Units,          SpaceDim,        &
                              OnLevelEdges, AddOffset,      MissingValue,    &
                              ScaleFactor,  Source_KindVal, Operation,       &
                              DimNames,     Dimensions,     Subset_X,        &
                              Subset_Y,     Subset_Z,       Source_0d_8,     &
                              Source_1d,    Source_1d_8,    Source_1d_4,     &
                              Source_1d_I,  Source_2d,      Source_2d_8,     &
                              Source_2d_4,  Source_2d_I,    Source_3d,       &
                              Source_3d_8,  Source_3d_4,    Source_3d_I     )
!
! !USES:
!
  USE CharPak_Mod,         ONLY : TranLc
  USE ErrCode_Mod
  USE History_Util_Mod
    USE Input_Opt_Mod,     ONLY : OptInput
  USE Registry_Params_Mod
!
! !INPUT PARAMETERS:
!
    ! Required arguments
    TYPE(OptInput),    INTENT(IN)  :: Input_Opt          ! Input Options object
    INTEGER,           INTENT(IN)  :: Id                 ! History item Id #
    INTEGER,           INTENT(IN)  :: ContainerId        ! Container Id #
    CHARACTER(LEN=*),  INTENT(IN)  :: Name               ! Item's short name
    CHARACTER(LEN=*),  INTENT(IN)  :: LongName           ! Item's long name
    CHARACTER(LEN=*),  INTENT(IN)  :: Units              ! Units of the data
    INTEGER,           INTENT(IN)  :: SpaceDim           ! Dimension of data
    INTEGER,           INTENT(IN)  :: Subset_X(2)        ! X0, X1 indices
    INTEGER,           INTENT(IN)  :: Subset_Y(2)        ! Y0, Y1 indices
    INTEGER,           INTENT(IN)  :: Subset_Z(2)        ! Z0, Z1 indices

    ! Optional arguments
    LOGICAL,           OPTIONAL    :: OnLevelEdges       ! =T if data defined
                                                         !  on level edges;
                                                         ! =F if on centers
    REAL(f4),          OPTIONAL    :: AddOffset          ! COARDS-compliant
    REAL(f4),          OPTIONAL    :: MissingValue       !  attributes for
    REAL(f4),          OPTIONAL    :: ScaleFactor        !  netCDF output
    INTEGER,           OPTIONAL    :: Operation          ! Operation code
                                                         !  0=copy  from source
                                                         !  1=accum from source
    CHARACTER(LEN=*),  OPTIONAL    :: DimNames           ! Use this to specify
                                                         !  dimensions of data
                                                         !  ("yz", "z", etc.)

    ! Optional pointers to data targets
    INTEGER,           OPTIONAL    :: Source_KindVal     ! Type of source data
    REAL(fp), POINTER, OPTIONAL    :: Source_0d_8        ! 0D 8-byte    data
    REAL(fp), POINTER, OPTIONAL    :: Source_1d  (:    ) ! 1D flex-prec data
    REAL(fp), POINTER, OPTIONAL    :: Source_1d_8(:    ) ! 1D 8-byte    data
    REAL(f4), POINTER, OPTIONAL    :: Source_1d_4(:    ) ! 1D 4-byte    data
    INTEGER,  POINTER, OPTIONAL    :: Source_1d_I(:    ) ! 1D integer   data
    REAL(fp), POINTER, OPTIONAL    :: Source_2d  (:,:  ) ! 2D flex-prec data
    REAL(f8), POINTER, OPTIONAL    :: Source_2d_8(:,:  ) ! 2D 8-byte    data
    REAL(f4), POINTER, OPTIONAL    :: Source_2d_4(:,:  ) ! 2D 4-byte    data
    INTEGER,  POINTER, OPTIONAL    :: Source_2d_I(:,:  ) ! 2D integer   data
    REAL(fp), POINTER, OPTIONAL    :: Source_3d  (:,:,:) ! 3D flex-prec data
    REAL(f8), POINTER, OPTIONAL    :: Source_3d_8(:,:,:) ! 3D 8-byte    data
    REAL(f4), POINTER, OPTIONAL    :: Source_3d_4(:,:,:) ! 3D 4-byte    data
    INTEGER,  POINTER, OPTIONAL    :: Source_3d_I(:,:,:) ! 3D integer   data
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(HistItem),    POINTER     :: Item               ! HISTORY ITEM object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,           OPTIONAL    :: Dimensions(3)      ! Spatial dims of data
    INTEGER,           INTENT(OUT) :: RC                 ! Success or failure
!
! !REMARKS:
!  (1) We need to copy string data to a temporary string of length 255
!       characters, or else Gfortran will choke.
!
! !REVISION HISTORY:
!  13 Jun 2017 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    LOGICAL            :: Is_DimNames
    LOGICAL            :: Is_0d_8
    LOGICAL            :: Is_1d,      Is_1d_8,  Is_1d_4,  Is_1d_I
    LOGICAL            :: Is_2d,      Is_2d_8,  Is_2d_4,  Is_2d_I
    LOGICAL            :: Is_3d,      Is_3d_8,  Is_3d_4,  Is_3d_I
    INTEGER            :: X0,         X1,       Y0,       Y1
    INTEGER            :: Z0,         Z1,       N

    ! Arrays
    INTEGER            :: Dims(3)

    ! Strings
    CHARACTER(LEN=3  ) :: TmpDimNames
    CHARACTER(LEN=255) :: ErrMsg,     ThisLoc,  TempStr

    !========================================================================
    ! Initialize
    !========================================================================
    RC               =  GC_SUCCESS
    Dims             =  UNDEFINED_INT
    ErrMsg           =  ''
    ThisLoc          =  ' -> at HistItem_Create (in History/histitem_mod.F90)'

    ! Allocate the Item object
    ALLOCATE( Item, STAT=RC )
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Cannot allocate the "Item" object!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    ! Zero the data fields
    Item%Data_0d     => NULL()
    Item%Data_1d     => NULL()
    Item%Data_2d     => NULL()
    Item%Data_3d     => NULL()

    ! Determine if the optional source pointers are passed
    Is_0d_8          =  PRESENT( Source_0d_8 )
    Is_1d            =  PRESENT( Source_1d   )
    Is_1d_8          =  PRESENT( Source_1d_8 )
    Is_1d_4          =  PRESENT( Source_1d_4 )
    Is_1d_I          =  PRESENT( Source_1d_I )
    Is_2d            =  PRESENT( Source_2d   )
    Is_2d_8          =  PRESENT( Source_2d_8 )
    Is_2d_4          =  PRESENT( Source_2d_4 )
    Is_2d_I          =  PRESENT( Source_2d_I )
    Is_3d            =  PRESENT( Source_3d   )
    Is_3d_8          =  PRESENT( Source_3d_8 )
    Is_3d_4          =  PRESENT( Source_3d_4 )
    Is_3d_I          =  PRESENT( Source_3d_I )

    ! Zero optional source pointers
    IF ( Is_0d_8 ) Item%Source_0d_8 => NULL()
    IF ( Is_1d   ) Item%Source_1d   => NULL()
    IF ( Is_1d_8 ) Item%Source_1d_8 => NULL()
    IF ( Is_1d_4 ) Item%Source_1d_4 => NULL()
    IF ( Is_1d_I ) Item%Source_1d_I => NULL()
    IF ( Is_2d   ) Item%Source_2d   => NULL()
    IF ( Is_2d_8 ) Item%Source_2d_8 => NULL()
    IF ( Is_2d_4 ) Item%Source_2d_4 => NULL()
    IF ( Is_2d_I ) Item%Source_2d_I => NULL()
    IF ( Is_3d   ) Item%Source_3d   => NULL()
    IF ( Is_3d_8 ) Item%Source_3d_8 => NULL()
    IF ( Is_3d_4 ) Item%Source_3d_4 => NULL()
    IF ( Is_3d_I ) Item%Source_3d_I => NULL()

    ! Initialize indices
    X0 = UNDEFINED_INT
    X1 = UNDEFINED_INT
    X1 = UNDEFINED_INT
    Y1 = UNDEFINED_INT
    Z1 = UNDEFINED_INT
    Z1 = UNDEFINED_INT

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

    !--------------------------------------------
    ! These won't be defined until we enter
    ! netCDF define mode, so set them undefined
    !--------------------------------------------
    Item%NcXDimId = UNDEFINED_INT
    Item%NcYDimId = UNDEFINED_INT
    Item%NcZDimId = UNDEFINED_INT
    ITem%NcIDimId = UNDEFINED_INT
    Item%NcTdimId = UNDEFINED_INT
    Item%NcVarId  = UNDEFINED_INT

    !========================================================================
    ! Optional inputs, handle these next
    !========================================================================

    !--------------------------------------------
    ! OnLevelEdges
    !--------------------------------------------
    IF ( PRESENT( OnLevelEdges ) ) THEN
       Item%OnLevelEdges = OnLevelEdges
    ELSE
       Item%OnLevelEdges = .FALSE.
    ENDIF

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

    !--------------------------------------------
    ! Source_KindVal
    !--------------------------------------------
    IF ( PRESENT( Source_KindVal ) ) THEN
       Item%Source_KindVal = Source_KindVal
    ELSE
       Item%Source_KindVal = KINDVAL_FP
    ENDIF

    !--------------------------------------------
    ! Source_KindVal
    !--------------------------------------------
    IF ( PRESENT( Source_KindVal ) ) THEN
       Item%Source_KindVal = Source_KindVal
    ELSE
       Item%Source_KindVal = KINDVAL_FP
    ENDIF

    !--------------------------------------------
    ! Operation
    !--------------------------------------------
    IF ( PRESENT( Operation ) ) THEN
       Item%Operation = Operation
    ELSE
       Item%Operation = COPY_FROM_SOURCE
    ENDIF

    !--------------------------------------------
    ! DimNames
    !--------------------------------------------
    IF ( PRESENT( DimNames ) ) THEN
       TmpDimNames = DimNames
       CALL TranLc( TmpDimNames )
    ELSE
       TmpDimNames = '  '
    ENDIF

    !--------------------------------------------
    ! Averaging method (define from Operation)
    !--------------------------------------------
    IF ( Item%Operation == COPY_FROM_SOURCE ) THEN
       TempStr        = 'instantaneous'
       Item%AvgMethod = TempStr
    ELSE
       TempStr        = 'time-averaged'
       Item%AvgMethod = TempStr
    ENDIF

    !========================================================================
    ! Make sure the spatial dimension is in the range 0-3
    !========================================================================
    IF ( SpaceDim < 0 .or. SpaceDim > 3 ) THEN
       ErrMsg = 'SpaceDim must be in the range 0-3!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ELSE
       Item%SpaceDim = SpaceDim
    ENDIF

    !=======================================================================
    ! Set the DimNames field of the HISTORY ITEM
    !=======================================================================
    IF ( PRESENT( DimNames ) ) THEN

       ! If the DimNames argument was passed, use it
       Item%DimNames = DimNames

    ELSE

       ! Set default values if the DimNames argument isn't passed
       ! Most of the time we deal with either xy or xyz spatial data
       SELECT CASE( Item%SpaceDim )
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

    !========================================================================
    ! Attach pointers to the data source.  Also get the values of NX, NY,
    ! and NZ from the relevant source pointer if they were not passed.
    !========================================================================
    SELECT CASE( Item%SpaceDim )

       ! Attach pointer to 3D data source, depending on its type
       CASE( 3 )

          ! Subsets: assume all 3d are xyz
          X0 = Subset_X(1); X1 = Subset_X(2)
          Y0 = Subset_Y(1); Y1 = Subset_Y(2)
          Z0 = Subset_Z(1); Z1 = Subset_Z(2)

          IF ( Item%Source_KindVal == KINDVAL_FP ) THEN
             IF ( Is_3d   ) THEN
                Item%Source_3d => Source_3d(X0:X1, Y0:Y1, Z0:Z1)
                DO N = 1, Item%SpaceDim
                   Dims(N) = SIZE( Source_3d(X0:X1, Y0:Y1, Z0:Z1), N )
                ENDDO
                GOTO 99
             ENDIF
          ELSE IF ( Item%Source_KindVal == KINDVAL_F8 ) THEN
             IF ( Is_3d_8 ) THEN
                Item%Source_3d_8 => Source_3d_8(X0:X1, Y0:Y1, Z0:Z1)
                DO N = 1, Item%SpaceDim
                   Dims(N) = SIZE( Source_3d_8(X0:X1, Y0:Y1, Z0:Z1), N )
                ENDDO
                GOTO 99
             ENDIF
          ELSE IF ( Item%Source_KindVal == KINDVAL_F4 ) THEN
             IF ( Is_3d_4 ) THEN
                Item%Source_3d_4 => Source_3d_4(X0:X1, Y0:Y1, Z0:Z1)
                DO N = 1, Item%SpaceDim
                   Dims(N) = SIZE( Source_3d_4(X0:X1, Y0:Y1, Z0:Z1), N )
                ENDDO
                GOTO 99
             ENDIF
          ELSE IF ( Item%Source_KindVal == KINDVAL_I4 ) THEN
             IF ( Is_3d_I ) THEN
                Item%Source_3d_I => Source_3d_I(X0:X1, Y0:Y1, Z0:Z1)
                DO N = 1, Item%SpaceDim
                   Dims(N) = SIZE( Source_3d_I(X0:X1, Y0:Y1, Z0:Z1), N )
                ENDDO
                GOTO 99
             ENDIF
          ENDIF

       ! Attach pointer to 2D data source, depending on its type
       CASE( 2 )

          ! Subsets: assume all 2D are xy (for now)
          X0 = Subset_X(1); X1 = Subset_X(2)
          Y0 = Subset_Y(1); Y1 = Subset_Y(2)

          IF ( Item%Source_KindVal == KINDVAL_FP ) THEN
             IF ( Is_2d   ) THEN
                Item%Source_2d => Source_2d(X0:X1, Y0:Y1)
                DO N = 1, Item%SpaceDim
                   Dims(N) = SIZE( Source_2d(X0:X1, Y0:Y1), N )
                ENDDO
                GOTO 99
             ENDIF
          ELSE IF ( Item%Source_KindVal == KINDVAL_F8 ) THEN
             IF ( Is_2d_8 ) THEN
                Item%Source_2d_8 => Source_2d_8(X0:X1, Y0:Y1)
                DO N = 1, Item%SpaceDim
                   Dims(N) = SIZE( Source_2d_8(X0:X1, Y0:Y1), N )
                ENDDO
                GOTO 99
             ENDIF
          ELSE IF ( Item%Source_KindVal == KINDVAL_F4 ) THEN
             IF ( Is_2d_4 ) THEN
                Item%Source_2d_4 => Source_2d_4(X0:X1, Y0:Y1)
                DO N = 1, Item%SpaceDim
                   Dims(N) = SIZE( Source_2d_4(X0:X1, Y0:Y1), N )
                ENDDO
                GOTO 99
             ENDIF
          ELSE IF ( Item%Source_KindVal == KINDVAL_I4 ) THEN
             IF ( Is_2d_I ) THEN
                Item%Source_2d_I => Source_2d_I(X0:X1, Y0:Y1)
                DO N = 1, Item%SpaceDim
                   Dims(N) = SIZE( Source_2d_I(X0:X1, Y0:Y1), N )
                ENDDO
                GOTO 99
             ENDIF
          ENDIF

       ! Attach pointer to 1D data source, depending on its type
       CASE( 1 )

          ! Subsets
          SELECT CASE( TRIM( Item%DimNames ) )
             CASE( 'x' )
                X0 = Subset_X(1); X1 = Subset_X(2)
             CASE( 'y' )
                X0 = Subset_Y(1); X1 = Subset_Y(2)
             CASE( 'z' )
                X0 = Subset_Z(1); X1 = Subset_Z(2)
             CASE DEFAULT
                X0 = 1;           X1 = 1
          END SELECT

          IF ( Item%Source_KindVal == KINDVAL_FP ) THEN
             IF ( Is_1d   ) THEN
                Item%Source_1d   => Source_1d(X0:X1)
                DO N = 1, Item%SpaceDim
                   Dims(N) = SIZE( Source_1d(X0:X1), N )
                ENDDO
                GOTO 99
             ENDIF
          ELSE IF ( Item%Source_KindVal == KINDVAL_F8 ) THEN
             IF ( Is_1d_8 ) THEN
                Item%Source_1d_8 => Source_1d_8(X0:X1)
                DO N = 1, Item%SpaceDim
                   Dims(N) = SIZE( Source_1d_8(X0:X1), N )
                ENDDO
                GOTO 99
             ENDIF
          ELSE IF ( Item%Source_KindVal == KINDVAL_F4 ) THEN
             IF ( Is_1d_4 ) THEN
                Item%Source_1d_4 => Source_1d_4(X0:X1)
                DO N = 1, Item%SpaceDim
                   Dims(N) = SIZE( Source_1d_4(X0:X1), N )
                ENDDO
                GOTO 99
             ENDIF
          ELSE IF ( Item%Source_KindVal == KINDVAL_I4 ) THEN
             IF ( Is_1d_I ) THEN
                Item%Source_1d_I => Source_1d_I(X0:X1)
                DO N = 1, Item%SpaceDim
                   Dims(N) = SIZE( Source_1d_I(X0:X1), N )
                ENDDO
                GOTO 99
             ENDIF
          ENDIF

       ! Attach pointer to 0D data source, depending on its type
       CASE( 0 )
          IF ( Item%Source_KindVal == KINDVAL_F8 ) THEN
             IF ( Is_0d_8 ) THEN
                Item%Source_0d_8 => Source_0d_8
                Dims             =  0
             ENDIF
          ENDIF

    END SELECT

    !=======================================================================
    ! Data fields: Allocate data fields (0-3 dimensions)
    !=======================================================================
 99 CONTINUE

    ! Allocate data field, based on SpaceDim
    SELECT CASE( Item%SpaceDim )

       !------------
       ! 3-D data
       !------------
       CASE( 3 )

          ! Allocate the data array
          ALLOCATE( Item%Data_3d( Dims(1), Dims(2), Dims(3) ), STAT=RC )
          IF ( RC == GC_SUCCESS ) THEN
             Item%Data_3d = 0.0_f8
          ELSE
             ErrMsg = 'Could not allocate "Item%Data_3d" array!'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF

          ! Allocate the NcChunkSizes array
          ALLOCATE( Item%NcChunkSizes( 4 ), STAT=rC )
          IF ( RC == GC_SUCCESS ) THEN
             Item%NcChunkSizes = (/ Dims(1), Dims(2), 1, 1 /)
          ELSE
             ErrMsg = 'Could not allocate "Item%NcChunkSizes" array (3d)!'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF

       !------------
       ! 2-D data
       !------------
       CASE( 2 )
          ALLOCATE( Item%Data_2d( Dims(1), Dims(2) ), STAT=RC )
          IF ( RC == GC_SUCCESS ) THEN
             Item%Data_2d = 0.0_f8
          ELSE
             ErrMsg = 'Could not allocate "Item%Data_2d" array!'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF

          ! Allocate the NcChunkSizes array
          ALLOCATE( Item%NcChunkSizes( 3 ), STAT=rC )
          IF ( RC == GC_SUCCESS ) THEN
             IF ( TRIM( Item%DimNames ) == 'xy' ) THEN
                Item%NcChunkSizes = (/ Dims(1), Dims(2), 1 /)   ! xy
             ELSE
                Item%NcChunkSizes = (/ Dims(1), 1,       1 /)   ! xz or yz
             ENDIF

          ELSE
             ErrMsg = 'Could not allocate "Item%NcChunkSizes" array (2d)!'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF

       !------------
       ! 1-D data
       !------------
       CASE( 1 )
          ALLOCATE( Item%Data_1d( Dims(1) ), STAT=RC )
          IF ( RC == GC_SUCCESS ) THEN
             Item%Data_1d  = 0.0_f8
          ELSE
             ErrMsg = 'Could not allocate "Item%Data_1d" array!'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF

          ! Allocate the NcChunkSizes array
          ALLOCATE( Item%NcChunkSizes( 2 ), STAT=rC )
          IF ( RC == GC_SUCCESS ) THEN
             IF ( TRIM( Item%DimNames ) == 'z' ) THEN
                Item%NcChunkSizes = (/ 1,       1 /)   ! z
             ELSE
                Item%NcChunkSizes = (/ Dims(1), 1 /)   ! x or y
             ENDIF
          ELSE
             ErrMsg = 'Could not allocate "Item%NcChunkSizes" array (1d)!'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF

       !------------
       ! 0-D data
       !------------
       CASE( 0 )
          ALLOCATE( Item%Data_0d, STAT=RC )
          IF ( RC == GC_SUCCESS ) THEN
             Item%Data_0d = 0.0_f8
          ELSE
             ErrMsg = 'Could not allocate "Item%Data_0d" variable!'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF

          ! Allocate the NcChunkSizes array
          ALLOCATE( Item%NcChunkSizes( 1 ), STAT=RC )
          IF ( RC == GC_SUCCESS ) THEN
             Item%NcChunkSizes = (/ 1 /)
          ELSE
             ErrMsg = 'Could not allocate "Item%NcChunkSizes" array (0d)!'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF

    END SELECT

    !=======================================================================
    ! Return to the calling program the spatial dimensions of the data
    ! if the optional DIMENSIONS argument has been passed
    !=======================================================================
    IF ( PRESENT( Dimensions ) ) Dimensions = Dims

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
  SUBROUTINE HistItem_Print( Input_Opt, Item, RC, ShortFormat )
!
! !USES:
!
    USE ErrCode_Mod
    USE History_Util_Mod
    USE Input_Opt_Mod,   ONLY : OptInput
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)  :: Input_Opt    ! Input Options object
    TYPE(HistItem), POINTER     :: Item         ! History Item
    LOGICAL,        OPTIONAL    :: ShortFormat  ! Print truncated format
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT) :: RC           ! Success or failure?
!
! !REVISION HISTORY:
!  13 Jun 2017 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARAIBLES
!
    ! Scalars
    LOGICAL          :: Use_ShortFormat

    ! Strings
    CHARACTER(LEN=1) :: CellPos

    !=======================================================================
    ! Initialize
    !=======================================================================

    ! Assume success
    RC      = GC_SUCCESS
    CellPos = ''

    ! Define
    IF ( PRESENT( ShortFormat ) ) THEN
       Use_ShortFormat = ShortFormat
    ELSE
       Use_ShortFormat = .FALSE.
    ENDIF

    !=======================================================================
    ! Print information about this HISTORY ITEM (only on the root CPU)
    !=======================================================================
    IF ( ASSOCIATED( Item ) .and. Input_Opt%amIRoot ) THEN

       IF ( Use_ShortFormat ) THEN

          !-----------------------------------------------------------------
          ! Use truncated output format
          !-----------------------------------------------------------------

          ! Denote if the data is defined on
          ! level edges (E) or centers (C)
          IF ( Item%SpaceDim == 3 ) THEN
             IF ( Item%OnLevelEdges ) THEN
                CellPos = 'E'
             ELSE
                CellPos = 'C'
             ENDIF
          ENDIF

          ! Print information
          WRITE( 6, 100 ) Item%Name,     Item%LongName,                      &
                          Item%DimNames, CellPos,       TRIM( Item%Units )
 100      FORMAT( 2x, a20, ' | ', a38, ' | ', a3, ' ', a1, ' | ', a )

       ELSE

          !-----------------------------------------------------------------
          ! Use expanded output format
          !-----------------------------------------------------------------
          PRINT*, REPEAT( '-', 70 )
          PRINT*, 'Name           : ', TRIM( Item%Name     )
          PRINT*, 'Long_Name      : ', TRIM( Item%LongName )
          PRINT*, 'Units          : ', TRIM( Item%Units    )
          PRINT*, 'OnLevelEdges   : ', Item%OnLevelEdges
          PRINT*, 'Add_Offset     : ', Item%AddOffset
          PRINT*, 'Missing_Value  : ', Item%MissingValue
          PRINT*, 'Scale_Factor   : ', Item%ScaleFactor
          PRINT*, ''
          PRINT*, 'Id             : ', Item%ID
          PRINT*, 'CollectionId   : ', Item%ContainerId
          PRINT*, 'NetCDF var ID  : ', Item%NcVarId
          PRINT*, 'NetCDF xDim Id : ', Item%NcXDimId
          PRINT*, 'NetCDF yDim Id : ', Item%NcYDimId
          PRINT*, 'NetCDF zDim Id : ', Item%NcZDimId
          PRINT*, 'NetCDF iDim Id : ', Item%NcIDimId
          PRINT*, 'NetCDF tDim Id : ', Item%NcTDimId
          PRINT*, ''
          PRINT*, 'nUpdates       : ', Item%nUpdates
          PRINT*, 'Operation      : ', Item%Operation
          PRINT*, ''
          PRINT*, 'SpaceDim       : ', Item%SpaceDim, ' (', Item%DimNames, ')'
          PRINT*, 'NcChunkSizes   : ', Item%NcChunkSizes

          IF ( ASSOCIATED( Item%Data_0d ) ) THEN
             PRINT*, 'Value Data_0d  : ', Item%Data_0d
          ENDIF

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
  SUBROUTINE HistItem_Destroy( Item, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE History_Util_Mod
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
    ErrMsg = ''
    ThisLoc = ' -> at HistItem_Destroy (in History/histitem_mod.F90)'

    !=======================================================================
    ! Nullify fields that are just pointing to other objects
    !======================================================================
    Item%Source_0d_8 => NULL()
    Item%Source_1d   => NULL()
    Item%Source_1d_8 => NULL()
    Item%Source_1d_4 => NULL()
    Item%Source_1d_I => NULL()
    Item%Source_2d   => NULL()
    Item%Source_2d_8 => NULL()
    Item%Source_2d_4 => NULL()
    Item%Source_2d_I => NULL()
    Item%Source_3d   => NULL()
    Item%Source_3d_8 => NULL()
    Item%Source_3d_4 => NULL()
    Item%Source_3d_I => NULL()

    !=======================================================================
    ! Free allocated pointer-based fields
    !=======================================================================
    IF ( ASSOCIATED( Item%Data_3d ) ) THEN
       DEALLOCATE( Item%Data_3d, STAT=RC )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Could not deallocate "Item%Data_3d" array!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
    ENDIF

    IF ( ASSOCIATED( Item%Data_2d ) ) THEN
       DEALLOCATE( Item%Data_2d, STAT=RC )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Could not deallocate "Item%Data_2d" array!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
    ENDIF

    IF ( ASSOCIATED( Item%Data_1d ) ) THEN
       DEALLOCATE( Item%Data_1d, STAT=RC )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Could not deallocate "Item%Data_1d" array!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
    ENDIF

    IF ( ASSOCIATED( Item%Data_0d ) ) THEN
       DEALLOCATE( Item%Data_0d, STAT=RC )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Could not deallocate "Item%Data_0d" variable!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
    ENDIF

    IF ( ASSOCIATED( Item%NcChunkSizes ) ) THEN
       DEALLOCATE( Item%NcChunkSizes, STAT=RC )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Could not deallocate "Item%NcChunkSizes" array!'
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
