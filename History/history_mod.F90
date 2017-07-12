!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: history_mod.F90
!
! !DESCRIPTION: Driver module for the new netCDF diagnostics package.
!\\
!\\
! !INTERFACE:
!
MODULE History_Mod
!
! !USES:
!
  USE History_Params_Mod
  USE Precision_Mod
  USE HistContainer_MOd,     ONLY : HistContainer
  USE MetaHistContainer_Mod, ONLY : MetaHistContainer

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
! 
  PUBLIC  :: History_ReadCollectionNames
  PUBLIC  :: History_ReadCollectionData
  PUBLIC  :: History_Cleanup
!
! PRIVATE MEMBER FUNCTIONS:

  PRIVATE :: CleanText
  PRIVATE :: ReadOneLine
!
! !REMARKS:
!
! !REVISION HISTORY:
!  06 Jan 2015 - R. Yantosca - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
  ! Scalars
  INTEGER                         :: CollectionCount

  ! Strings
  CHARACTER(LEN=255), ALLOCATABLE :: CollectionName      (:)
  CHARACTER(LEN=255), ALLOCATABLE :: CollectionTemplate  (:)
  CHARACTER(LEN=255), ALLOCATABLE :: CollectionSubsetDims(:)
  CHARACTER(LEN=255), ALLOCATABLE :: CollectionFormat    (:)
  CHARACTER(LEN=255), ALLOCATABLE :: CollectionFrequency (:)
  CHARACTER(LEN=255), ALLOCATABLE :: CollectionDuration  (:)
  CHARACTER(LEN=255), ALLOCATABLE :: CollectionMode      (:)

  ! Objects
  TYPE(MetaHistContainer), POINTER     :: CollectionList

CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Read_Collection_Names
!
! !DESCRIPTION: tbd
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE History_ReadCollectionNames( am_I_Root, HistoryInputFile, RC )
!
! !USES:
!
    USE Charpak_Mod
    USE ErrCode_Mod
    USE InquireMod,  ONLY : FindFreeLun
!
! !INPUT PARAMETERS: 
!
    LOGICAL,                 INTENT(IN)  :: am_I_Root
    CHARACTER(LEN=255), INTENT(IN)  :: HistoryInputFile
!
! !OUTPUT PARAMETERS: 
!
    INTEGER,                 INTENT(OUT) :: RC
!
! !REMARKS:
!
! !REVISION HISTORY:
!  16 Jun 2017 - R. Yantosca - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
    ! Ridiculously big number
    INTEGER, PARAMETER      :: MAX_COLLECTIONS = 500
!
! !LOCAL VARIABLES:
!
    ! Scalars
    LOGICAL                 :: EOF
    INTEGER                 :: fId,    IOS      
    INTEGER                 :: N,      nSubs1,  nSubs2

    ! Strings
    CHARACTER(LEN=255)      :: ErrMsg, ThisLoc, Line,  Line2
    
    ! String arrays
    CHARACTER(LEN=255)      :: Subs1(255)
    CHARACTER(LEN=255)      :: Subs2(255)
    CHARACTER(LEN=255)      :: TmpCollectionName(MAX_COLLECTIONS)

    !=======================================================================
    ! Initialize
    !=======================================================================

    ! Assume success
    RC      = GC_SUCCESS

    ! For error output
    ErrMsg  = ''
    ThisLoc = &
     ' -> at History_ReadCollectionNames (in module History/history_mod.F90)'

    ! Zero global variables
    CollectionCount   = 0
    TmpCollectionName = ''

    !=======================================================================
    ! Open the file containing the list of requested diagnostics
    !=======================================================================

    ! Find a free file unit
    fId     = FindFreeLun()

    ! Open the file
    OPEN( fId, FILE=TRIM( HistoryInputFile ), STATUS='OLD', IOSTAT=RC )
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Could not open "HistoryInputFile"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    !=======================================================================
    ! Read data from the file
    !=======================================================================
    DO

       ! Read a single line, and strip leading/trailing spaces
       Line = ReadOneLine( fId, EOF, IOS, Squeeze=.TRUE. )

       ! Exit the loop if it's the end of the file
       IF ( EOF ) EXIT

       ! If it's a real I/O error, quit w/ error message
       IF ( IOS > 0 ) THEN
          ErrMsg = 'Unexpected end-of-file in ' // &
                    TRIM( HistoryInputFile )    //'!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       ! Skip if the line is commented out
       IF ( Line(1:1) == "#" ) CYCLE
      
       !-------------------------------------------------------------------
       ! Get the list of collections
       !-------------------------------------------------------------------
       IF ( Line(1:12) == 'COLLECTIONS:' ) THEN

          ! Start
          CollectionCount                    = CollectionCount + 1
          TmpCollectionName(CollectionCount) = CleanText( Line(14:) )
          
          ! Loop over all collections
          DO

             ! Read next line and strip leading/trailing blanks
             Line2 = ReadOneLine( fId, EOF, IOS, Squeeze=.TRUE. )

             ! Exit the loop if it's the end of the file
             IF ( EOF ) EXIT

             ! If it's a real I/O error, quit w/ error message
             IF ( IOS > 0 ) THEN
                ErrMsg = 'Unexpected erroe in ' // &
                     TRIM( HistoryInputFile )    //'!'
                CALL GC_Error( ErrMsg, RC, ThisLoc )
                RETURN
             ENDIF

             ! Skip if the line is commented out
             IF ( Line2(1:1) == "#" ) CYCLE
      
             ! 2 colons signal the end of the collecto
             IF ( Line2(1:2) == '::' ) EXIT

             CollectionCount                    = CollectionCount + 1
             TmpCollectionName(CollectionCount) = CleanText( Line2 )
          ENDDO
       ENDIF
    ENDDO

    !=======================================================================
    ! Cleanup and quit
    !=======================================================================

    ! Close the file
    CLOSE( fId )

    ! Now that we now the number of diagnostic collections, we can allocate
    ! the array that will hold the name of each diagnostic collection.
    IF ( .not. ALLOCATED( CollectionName ) ) THEN
       ALLOCATE( CollectionName( CollectionCount ), STAT=RC )
       IF ( RC /= GC_SUCCESS ) THEN 
          ErrMsg = 'Could not allocate "CollectionName"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
    ENDIF

    ! Copy the collection names from the temporary array
    DO N = 1, CollectionCount
       CollectionName(N) = TmpCollectionName(N)
       WRITE( 6, '(i3, ":", a )' ) N, TRIM( CollectionName(N) )
    ENDDO

    ! Allocate CollectionTemplate
    IF ( .not. ALLOCATED( CollectionTemplate ) ) THEN
       ALLOCATE( CollectionTemplate( CollectionCount ), STAT=RC )
       IF ( RC /= GC_SUCCESS ) THEN 
          ErrMsg = 'Could not allocate "CollectionTemplate"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
       CollectionTemplate = ''
    ENDIF

    ! Allocate CollectionFormat
    IF ( .not. ALLOCATED( CollectionFormat ) ) THEN
       ALLOCATE( CollectionFormat( CollectionCount ), STAT=RC )
       IF ( RC /= GC_SUCCESS ) THEN 
          ErrMsg = 'Could not allocate "CollectionFormat"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
       CollectionFormat = ''
    ENDIF

    ! Allocate CollectionFrequency
    IF ( .not. ALLOCATED( CollectionFrequency ) ) THEN
       ALLOCATE( CollectionFrequency( CollectionCount ), STAT=RC )
       IF ( RC /= GC_SUCCESS ) THEN 
          ErrMsg = 'Could not allocate "CollectionFrequency"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
       CollectionFrequency = ''
    ENDIF

    ! Allocate CollectionDuration
    IF ( .not. ALLOCATED( CollectionDuration ) ) THEN
       ALLOCATE( CollectionDuration( CollectionCount ), STAT=RC )
       IF ( RC /= GC_SUCCESS ) THEN 
          ErrMsg = 'Could not allocate "CollectionDuration"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
       CollectionDuration = ''
    ENDIF

    ! Allocate CollectionSubsetDims
    IF ( .not. ALLOCATED( CollectionSubsetDims ) ) THEN
       ALLOCATE( CollectionSubSetDims( CollectionCount ), STAT=RC )
       IF ( RC /= GC_SUCCESS ) THEN 
          ErrMsg = 'Could not allocate "CollectionSubsetDims"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
       CollectionSubsetDims = ''
    ENDIF

    ! Allocate CollectionMode
    IF ( .not. ALLOCATED( CollectionMode ) ) THEN
       ALLOCATE( CollectionMode( CollectionCount ), STAT=RC )
       IF ( RC /= GC_SUCCESS ) THEN 
          ErrMsg = 'Could not allocate "CollectionMode"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
       CollectionMode = ''
    ENDIF

  END SUBROUTINE History_ReadCollectionNames
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Read_Collection_Data
!
! !DESCRIPTION: Reads the file 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE History_ReadCollectionData( am_I_Root, HistoryInputFile, RC )
!
! !USES:
!
    USE Charpak_Mod
    USE ErrCode_Mod
    USE InquireMod,            ONLY : FindFreeLun
    USE MetaHistContainer_Mod
    USE HistContainer_Mod
    USE HistContainer_Mod
    USE MetaHistItem_Mod,
    USE HistItem_Mod
!
! !INPUT PARAMETERS: 
!
    LOGICAL,            INTENT(IN)  :: am_I_Root
    CHARACTER(LEN=255), INTENT(IN)  :: HistoryInputFile
!
! !OUTPUT PARAMETERS: 
!
    INTEGER,            INTENT(OUT) :: RC
!
! !REMARKS:
!
! !REVISION HISTORY:
!  16 Jun 2017 - R. Yantosca - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
     ! Scalars
    LOGICAL                      :: EOF
    INTEGER                      :: C,           N,      W    
    INTEGER                      :: fId,         IOS      
    INTEGER                      :: nSubs1,      nSubs2
    INTEGER                      :: Ind1,        Ind2
    INTEGER                      :: ArchivalYmd, ArchivalHms
    INTEGER                      :: ItemCount,   SpaceDim
    

    ! Strings
    CHARACTER(LEN=255)           :: Line,        FileName
    CHARACTER(LEN=255)           :: ErrMsg,      ThisLoc
    CHARACTER(LEN=255)           :: MetaData,    Reference
    CHARACTER(LEN=255)           :: Title,       ItemName
    
    ! Arrays
    INTEGER                      :: SubsetDims(3)
    CHARACTER(LEN=255)           :: Subs1(255)
    CHARACTER(LEN=255)           :: Subs2(255)

    ! Objects
    TYPE(HistContainer), POINTER :: Collection
    TYPE(HistItem),      POINTER :: Item

    !=======================================================================
    ! Initialize
    !=======================================================================

    ! Assume success
    RC          =  GC_SUCCESS

    ! Initialize variables
    ArchivalYmd = 0
    ArchivalHms = 0
    SpaceDim    = 0
    SubsetDims  = 0

    ! Initialize objects
    Collection  => NULL()
    Item        => NULL()

    ! Initialize Strings
    ErrMsg      =  ''
    ThisLoc     =  &
     ' -> at History_ReadCollectionData (in module History/history_mod.F90)'
    Reference   = 'www.geos-chem.org; wiki.geos-chem.org'

    !=======================================================================
    ! Open the file containing the list of requested diagnostics
    !=======================================================================

    ! Find a free file unit
    fId     = FindFreeLun()

    ! Open the file
    OPEN( fId, FILE=TRIM( HistoryInputFile ), STATUS='OLD', IOSTAT=RC )
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Could not open "HistoryInputFile"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    !=======================================================================
    ! Read data from the file
    !=======================================================================
    DO

       ! Read a single line, and strip leading/trailing spaces
       Line = ReadOneLine( fId, EOF, IOS, Squeeze=.TRUE. )

       ! Exit the loop if it's the end of the file
       IF ( EOF ) EXIT

       ! If it's a real I/O error, quit w/ error message
       IF ( IOS > 0 ) THEN
          ErrMsg = 'Unexpected end-of-file in ' // &
                    TRIM( HistoryInputFile )    //'!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       ! Skip if the line is commented out
       IF ( Line(1:1) == "#" ) CYCLE

       !====================================================================
       ! The HISTORY.rc file specifies collection metadata as:
       !
       !   instantaneous.template:  '%y4%m2%d2.nc4',
       !   instantaneous.format:    'CFIO',
       !   instantaneous.frequency:  010000,
       !   instantaneous.duration:   240000
       !   etc.
       !
       ! where in this example, "instantaneous" is the collection name
       ! and "template", "format", "frequency", "duration" are the
       ! metadata fields.  
       !
       ! Get the metadata belonging to each collection and store them
       ! in the proper arrays for later use.  NOTE: this method does not
       ! assume that the collections are in the same order as they
       ! are listed under the COLLECTIONS section.
       !====================================================================
       IF ( INDEX( Line, 'template' ) > 0 ) THEN 
          CALL GetCollectionMetaData( Line, 'template',   MetaData, C )
          CollectionTemplate(C) = Metadata
       ENDIF

       IF ( INDEX( Line, 'format' ) > 0 ) THEN
          CALL GetCollectionMetaData( Line, 'format',     MetaData, C )
          CollectionFormat(C) = Metadata
       ENDIF

       IF ( INDEX( Line, 'frequency' ) > 0 ) THEN
          CALL GetCollectionMetaData( Line, 'frequency',  MetaData, C )
          CollectionFrequency(C) = Metadata
          READ( CollectionFrequency(C), '(i6.6)' ) ArchivalHms
       ENDIF

       IF ( INDEX( Line, 'duration' ) > 0 ) THEN
          CALL GetCollectionMetaData( Line, 'format',     MetaData, C )
          CollectionDuration(C) = Metadata
       ENDIF

       IF ( INDEX( Line, 'mode' ) > 0 ) THEN
          CALL GetCollectionMetaData( Line, 'mode',       MetaData, C )
          CollectionMode(C) = Metadata
       ENDIF

       IF ( INDEX( Line, 'subsetdims' ) > 0 ) THEN
          
          ! First split the line by colon
          CALL StrSplit( Line, ":", Subs1, nSubs1 )
          CollectionSubsetDims(C) = Subs1(2)

          ! Then split by spaces and convert to INTEGER
          CALL StrSplit( CollectionSubsetDims(C), " ", Subs2, nSubs2 )          
          IF ( nSubs2 > 0 ) THEN
             DO N = 1, nSubs2
                READ( Subs2(N), '(i10)' ) SubsetDims(N)
             ENDDO

             ! Define the number of dimensions
             SpaceDim = nSubs2
             IF ( SpaceDim == 2 ) SubsetDims(3) = 1
          ENDIF
       ENDIF

       !====================================================================
       ! NOTE: We assume FIELDS is the last metadata tag for the
       ! collection.  We need to create the collection object and
       ! an object for each history item stored in the collection.
       !====================================================================
       IF ( INDEX( Line, 'fields' ) > 0 ) THEN

          ! Zero the counter of items
          ItemCount = 0

          ! Create title string for collection
          Title = 'GEOS-Chem diagnostic collection: ' //                    &
                   TRIM( CollectionName(C) )
          
          !=================================================================
          ! Create a HISTORY CONTAINER object for this collection
          !=================================================================
          CALL HistContainer_Create( am_I_Root    = am_I_Root,              &
                                     Container    = Collection,             &
                                     Name         = CollectionName(C),      &
                                     Id           = C,                      &
                                     ArchivalYmd  = ArchivalYmd,            &
                                     ArchivalHms  = ArchivalHms,            &
                                     Conventions  = 'COARDS',               &
                                     FileTemplate = CollectionTemplate(C),  & 
                                     NcFormat     = 'netCDF-4',             &
                                     Reference    = Reference,              &
                                     Title        = Title,                  &
                                     RC           = RC                     )

          ! Trap potential error
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = 'Could not create Collection: ' // &
                      TRIM( CollectionName(C) ) 
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF
         
          !=================================================================
          ! Create a list of HISTORY ITEMS that will be contained in this
          ! collection, for each entry under the "fields" tag.
          !=================================================================
          DO 

             IF ( ItemCount == 0 ) THEN

                !----------------------------------------------------------
                ! If we are on the same line as the "fields" tag, then
                ! the name of the HISTORY ITEM will be the first substring
                ! of MetaData (split on spaces).
                !----------------------------------------------------------
                CALL GetCollectionMetaData( Line, 'fields', MetaData, C )
                CALL StrSplit( MetaData, " ", Subs1, nSubs1 )
                ItemName = Subs1(1)
             
             ELSE

                !----------------------------------------------------------
                ! Otherwise, read the next line to get the name for
                ! each subsequent HISTORY ITEM.  The name will be the
                ! first substring of the line.
                !----------------------------------------------------------

                ! Read a single line, and strip leading/trailing spaces
                Line = ReadOneLine( fId, EOF, IOS, Squeeze=.TRUE. )

                ! If it's a real I/O error, quit w/ error message
                IF ( IOS > 0 ) THEN
                   ErrMsg = 'Unexpected end-of-file in ' // &
                        TRIM( HistoryInputFile )    //'!'
                   CALL GC_Error( ErrMsg, RC, ThisLoc )
                   RETURN
                ENDIF
                             
                ! Remove commas, spaces, and tabs
                Line = CleanText( Line )    

                ! Skip if the line is commented out
                IF ( Line(1:1) == "#"  ) CYCLE

                ! 2 colons denotes the end of the "fields" list
                IF ( Line(1:2) == '::' ) EXIT

                ! The HISTORY ITEM name will be the 1st substring of Line
                CALL StrSplit( Line, " ", Subs1, nSubs1 )
                ItemName = Subs1(1)

             ENDIF

! Feature to be added
!            !------------
!            ! wild card
!            !------------
!
!            ! If there is a wild card then lop
!            W = INDEX( ItemName, '*' )
!            IF ( W > 0 ) THEN
!               N = nSpecies
!            ELSE
!               N = 1
!            ENDDO
!
!
!            DO N = 1, nspecies
!
!               IF ( W > 0 ) THEN
!                  ItemName = ItemName(0:N-1) // &
!                             SpeciesName     // &
!                             ItemName(N+1,:)
!               ENDIF

                ! Increment the number of HISTORY items
                ItemCount = ItemCount + 1

             !--------------------------------------------------------------
             ! After all this setup, create the HISTORY ITEM object itself
             !--------------------------------------------------------------
             CALL HistItem_Create( am_I_Root   = am_I_Root,                 &
                                   Item        = Item,                      &
                                   Id          = ItemCount,                 &
                                   ContainerId = C,                         &
                                   Name        = ItemName,                  &
                                   LongName    = ItemName,                  &
                                   Units       = 'TBD',                     &
                                   SpaceDim    = SpaceDim,                  &
                                   NX          = SubsetDims(1),             &
                                   NY          = SubsetDims(2),             &
                                   NZ          = SubsetDims(3),             & 
                                   RC          = RC                        )


             ! Trap potential error
             IF ( RC /= GC_SUCCESS ) THEN
                ErrMsg = 'Could not create Item: ' // TRIM( ItemName )
                CALL GC_Error( ErrMsg, RC, ThisLoc )
                RETURN
             ENDIF

             !--------------------------------------------------------------
             ! Attach this HISTORY ITEM to the METAHISTORY ITEM (aka list
             ! of HISTORY ITEMS) belonging to the HISTORY CONTAINER object
             ! for the given collection.
             !
             ! In other words, this denotes the list of fields that will
             ! be written to the netCDF file, and with the archiving
             ! frequency, specified by this given diagnostic collection.
             !--------------------------------------------------------------
             CALL MetaHistItem_AddNew( am_I_Root = am_I_Root,               &
                                       Node      = Collection%HistItems,    &
                                       Item      = Item,                    &
                                       RC        = RC                      )

             ! Trap potential error
             IF ( RC /= GC_SUCCESS ) THEN
                ErrMsg = 'Could not add Item' // TRIM( ItemName ) //        &
                     ' to ' // TRIM( CollectionName(C) ) // '%HistItems!'
                CALL GC_Error( ErrMsg, RC, ThisLoc )
                RETURN
             ENDIF

          
          ENDDO

          !=================================================================
          ! Add this HISTORY CONTAINER object (i.e. this collection) into
          ! the METAHISTORY OBJECT (i.e. the master list of collections).
          !=================================================================
          CALL MetaHistContainer_AddNew( am_I_Root   = am_I_Root,           &
                                         Node        = CollectionList,      &
                                         Container   = Collection,          &
                                         RC          = RC                  )

          ! Trap potential error
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = 'Could not add Container' //                          &
                      TRIM( CollectionName(C) ) //                          &
                      ' to the list of collections!'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF
       ENDIF

    ENDDO
       
    !=======================================================================
    ! Cleanup and quit
    !=======================================================================

    ! Close the file
    CLOSE( fId )

    DO C = 1, CollectionCount
       print*, 'Collection       ', TRIM( CollectionName      (C) )
       print*, '  -> Template    ', TRIM( CollectionTemplate  (C) )
       print*, '  -> Format      ', TRIM( CollectionFormat    (C) )
       print*, '  -> Frequency   ', TRIM( CollectionFrequency (C) )
       print*, '  -> Duration    ', TRIM( CollectionDuration  (C) )
       print*, '  -> Subset Dims ', TRIM( CollectionSubsetDims(C) )
       print*, '  -> Mode        ', TRIM( CollectionMode      (C) )
    ENDDO

    ! Print information about each diagnostic collection
    CALL MetaHistContainer_Print( am_I_Root, CollectionList, RC )

  END SUBROUTINE History_ReadCollectionData
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: CleanText
!
! !DESCRIPTION: Strips commas, apostrophes, spaces, and tabs from a string.
!\\
!\\
! !INTERFACE:
!
  FUNCTION CleanText( Str ) RESULT( CleanStr )
!
! !USES:
!
    USE Charpak_Mod, ONLY : CStrip, StrRepl, StrSqueeze
!
! !INPUT PARAMETERS: 
!
    CHARACTER(LEN=*), INTENT(IN) :: Str        ! Original string
!
! !RETURN VALUE
!
    CHARACTER(LEN=255)      :: CleanStr   ! Cleaned-up string
!
! !REMARKS:
!
! !REVISION HISTORY:
!  06 Jan 2015 - R. Yantosca - Initial version
!  21 Jun 2017 - R. Yantosca - Now call CSTRIP to remove tabs etc.
!EOP
!------------------------------------------------------------------------------
!BOC

    ! Initialize
    CleanStr = Str

    ! Strip out non-printing characters (e.g. tabs)
    CALL CStrip    ( CleanStr           )

    ! Remove commas and quotes
    CALL StrRepl   ( CleanStr, ",", " " )
    CALL StrRepl   ( CleanStr, "'", " " )

    ! Remove leading and trailing spaces
    CALL StrSqueeze( CleanStr           ) 

  END FUNCTION CleanText
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ReadOneLine
!
! !DESCRIPTION: Subroutine READ\_ONE\_LINE reads a line from the input file.  
!  If the global variable VERBOSE is set, the line will be printed to stdout.  
!  READ\_ONE\_LINE can trap an unexpected EOF if LOCATION is passed.  
!  Otherwise, it will pass a logical flag back to the calling routine, 
!  where the error trapping will be done.
!\\
!\\
! !INTERFACE:
!
  FUNCTION ReadOneLine( fId, EndOfFile, IoStatus, Squeeze ) RESULT( Line )
!
! !USES:
!
    USE Charpak_Mod, ONLY : StrSqueeze
!
! !INPUT PARAMETERS:
!
    INTEGER, INTENT(IN)  :: fId        ! File unit number
    LOGICAL, OPTIONAL    :: Squeeze    ! Call Strsqueeze?
!
! !OUTPUT PARAMETERS:
!
    LOGICAL, INTENT(OUT) :: EndOfFile  ! Denotes EOF condition
    INTEGER, INTENT(OUT) :: IoStatus   ! I/O status code
!
! !RETURN VALUE:
!
    CHARACTER(LEN=255)   :: Line       ! Single line from the input file
    ! 
! !REVISION HISTORY: 
!  16 Jun 2017 - R. Yantosca - Initial version, based on GEOS-Chem
!EOP
!------------------------------------------------------------------------------
!BOC

    !=================================================================
    ! Initialize
    !=================================================================
    EndOfFile = .FALSE.
    IoStatus  = 0
    Line      = ''

    !=================================================================
    ! Read data from the file
    !=================================================================

    ! Read a line from the file
    READ( fId, '(a)', IOSTAT=IoStatus ) Line

    ! IO Status < 0: EOF condition
    IF ( IoStatus < 0 ) THEN 
       EndOfFile = .TRUE.
       RETURN
    ENDIF

    ! If desired, call StrSqueeze to strip leading and trailing blanks
    IF ( PRESENT( Squeeze ) ) THEN
       IF ( Squeeze ) THEN
          CALL StrSqueeze( Line )
       ENDIF
    ENDIF

   END FUNCTION ReadOneLine
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: GetCollectionMetaData
!
! !DESCRIPTION: tbd
!\\
!\\
! !INTERFACE:
!
   SUBROUTINE GetCollectionMetaData( Line, Pattern, MetaData, nCollection )
!
! !USES:
!
     USE Charpak_Mod, ONLY: StrSplit
!
! !INPUT PARAMETERS: 
!
     CHARACTER(LEN=*),   INTENT(IN)  :: Line
     CHARACTER(LEN=*),   INTENT(IN)  :: Pattern
!
! !OUTPUT PARAMETERS:
!
     CHARACTER(LEN=255), INTENT(OUT) :: MetaData
     INTEGER,            INTENT(OUT) :: nCollection
!
! !REMARKS:
!
! !REVISION HISTORY:
!  06 Jan 2015 - R. Yantosca - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
     INTEGER                 :: C, Ind, nSubStr
     CHARACTER(LEN=255) :: SubStr(255)
     
     !======================================================================
     ! Find the metadata for the given collection
     !======================================================================

     ! Loop over all collection names
     DO C = 1, CollectionCount

        ! Then check to see which collection this is in
        Ind = INDEX( Line, TRIM( CollectionName(C) ) )

        ! If we find a match, then return it
        IF ( Ind > 0 ) THEN
           CALL StrSplit( Line, ':', SubStr, nSubStr )
           IF ( nSubStr == 2 ) THEN
              nCollection = C
              MetaData    = CleanText( SubStr(2) )
              EXIT
           ENDIF
        ENDIF
     ENDDO

 END SUBROUTINE GetCollectionMetaData
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: History_Cleanup
!
! !DESCRIPTION: Deallocates all module variables and objects.
!\\
!\\
! !INTERFACE:
!
   SUBROUTINE History_Cleanup( am_I_Root, RC )
!
! !USES:
!
     USE ErrCode_Mod
!
! !INPUT PARAMETERS: 
!
     LOGICAL, INTENT(IN)  :: am_I_Root
!
! !OUTPUT PARAMETERS: 
!
     INTEGER, INTENT(OUT) :: RC
!
! !REVISION HISTORY:
!  16 Jun 2017 - R. Yantosca - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
     ! Strings
     CHARACTER(LEN=255) :: ErrMsg, ThisLoc

     !======================================================================
     ! Initialize
     !======================================================================

     ! Assume success
     RC      = GC_SUCCESS

     ErrMsg  = ''
     ThisLoc = ' -> History_Cleanup (in module History/history_mod.F90'

     !======================================================================
     ! Deallocate module variables
     !======================================================================

     IF ( ASSOCIATED( CollectionList ) ) THEN
        DEALLOCATE( CollectionList, STAT=RC )
        IF ( RC /= GC_SUCCESS ) THEN
           ErrMsg = 'Could not deallocate "CollectionList"!'
           CALL GC_Error( ErrMsg, RC, ThisLoc )
           RETURN
        ENDIF        
     ENDIF

     IF ( ALLOCATED( CollectionName ) ) THEN
        DEALLOCATE( CollectionName, STAT=RC )
        IF ( RC /= GC_SUCCESS ) THEN
           ErrMsg = 'Could not deallocate "CollectionName"!'
           CALL GC_Error( ErrMsg, RC, ThisLoc )
           RETURN
        ENDIF
     ENDIF
     
     IF ( ALLOCATED( CollectionTemplate ) ) THEN
        DEALLOCATE( CollectionTemplate, STAT=RC )
        IF ( RC /= GC_SUCCESS ) THEN
           ErrMsg = 'Could not deallocate "CollectionTemplate"!'
           CALL GC_Error( ErrMsg, RC, ThisLoc )
           RETURN
        ENDIF
     ENDIF

     IF ( ALLOCATED( CollectionSubsetDims ) ) THEN
        DEALLOCATE( CollectionSubsetDims, STAT=RC )
        IF ( RC /= GC_SUCCESS ) THEN
           ErrMsg = 'Could not deallocate "CollectionSubsetDims"!'
           CALL GC_Error( ErrMsg, RC, ThisLoc )
           RETURN
        ENDIF
     ENDIF

     IF ( ALLOCATED( CollectionFormat ) ) THEN
        DEALLOCATE( CollectionFormat, STAT=RC )
        IF ( RC /= GC_SUCCESS ) THEN
           ErrMsg = 'Could not deallocate "CollectionFormat"!'
           CALL GC_Error( ErrMsg, RC, ThisLoc )
           RETURN
        ENDIF
     ENDIF

     IF ( ALLOCATED( CollectionFrequency ) ) THEN
        DEALLOCATE( CollectionFrequency, STAT=RC )
        IF ( RC /= GC_SUCCESS ) THEN
           ErrMsg = 'Could not deallocate "CollectionFrequency"!'
           CALL GC_Error( ErrMsg, RC, ThisLoc )
           RETURN
        ENDIF
     ENDIF

     IF ( ALLOCATED( CollectionDuration ) ) THEN
        DEALLOCATE( CollectionDuration, STAT=RC )
        IF ( RC /= GC_SUCCESS ) THEN
           ErrMsg = 'Could not deallocate "CollectionDuration"!'
           CALL GC_Error( ErrMsg, RC, ThisLoc )
           RETURN
        ENDIF
     ENDIF

     IF ( ALLOCATED( CollectionMode ) ) THEN
        DEALLOCATE( CollectionMode, STAT=RC )
        IF ( RC /= GC_SUCCESS ) THEN
           ErrMsg = 'Could not deallocate "CollectionMode"!'
           CALL GC_Error( ErrMsg, RC, ThisLoc )
           RETURN
        ENDIF
     ENDIF

   END SUBROUTINE History_Cleanup
!EOC
END MODULE History_Mod
