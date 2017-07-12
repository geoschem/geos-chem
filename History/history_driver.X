!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !PROGRAM
!
! !DESCRIPTION: 
!\\
!\\
! !INTERFACE:
!
PROGRAM History_Driver
!
! !USES:
!
  USE CharPak_Mod
  USE ErrCode_Mod
  USE Error_Mod,          ONLY: GC_Error
  USE History_Params_Mod
  USE History_Mod
  USE HistItem_Mod
  USE MetaHistItem_Mod
  USE HistContainer_Mod
  USE MetaHistContainer_Mod
  USE Precision_Mod


  IMPLICIT NONE
!
! !REMARKS:
!
! !REVISION HISTORY:
!  09 Jun 2015 - R. Yantosca - Initial version
!  12 Jun 2017 - R. Yantosca - Now simulates beter 

!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
  INTEGER,               PARAMETER :: N_COLL           = 6  ! Dummy vars for
  INTEGER,               PARAMETER :: N_ITEMS_PER_COLL = 20 ! sake of argument 
!
! !LOCAL VARIABLES:
!
  ! Scalars
  LOGICAL,               PARAMETER :: am_I_Root = .TRUE.    ! Root cpu?
  INTEGER                          :: RC                    ! Success?
  INTEGER                          :: nNodes
  INTEGER                          :: N

  ! Objects
  TYPE(HistItem),          POINTER :: Item                 ! History item
  TYPE(HistContainer),     POINTER :: Container            ! Collection
  TYPE(MetaHistContainer), POINTER :: ContainerList        ! List of coll's

  CHARACTER(LEN=HMAXCHAR)          :: CollName (N_COLL                  )
  INTEGER                          :: FreqYmd  (N_COLL                  )
  INTEGER                          :: FreqHms  (N_COLL                  )
  CHARACTER(LEN=HMAXCHAR)          :: ItemName (N_COLL, N_ITEMS_PER_COLL)
  CHARACTER(LEN=HMAXCHAR)          :: ItemLName(N_COLL, N_ITEMS_PER_COLL)
  CHARACTER(LEN=HMAXCHAR)          :: ItemUnits(N_COLL, N_ITEMS_PER_COLL)
  INTEGER                          :: ItemSpDim(N_COLL, N_ITEMS_PER_COLL)

  CHARACTER(LEN=HMAXCHAR)          :: ErrMsg, ThisLoc, FileName

  CHARACTER(LEN=HMAXCHAR)          :: SubStr(255)

  !=========================================================================
  ! Define a list of collections as if we were reading them from
  ! a file (that's the next step
  !=========================================================================

  ErrMsg           = ''
  ThisLoc          = ' -> in History/history_driver.F90'
 

  CollName         = ''
  FreqYmd          = -1
  FreqHms          = -1
  ItemName         = ''
  ItemLName        = ''
  ItemUnits        = ''
                   
  CollName (1  )   = 'Restart'
  FreqYmd  (1  )   = 000000
  FreqHms  (1  )   = 000000
  ItemName (1,1)   = 'CHEM_NO_CONC'
  ItemLName(1,1)   = 'NO concentration'
  ItemUnits(1,1)   = 'mol/mol'
  ItemSpDim(1,1)   = 3
  ItemName (1,2)   = 'CHEM_CO_CONC'
  ItemLName(1,2)   = 'CO concentration'
  ItemUnits(1,2)   = 'mol/mol'
  ItemSpDim(1,2)   = 3
  ItemName( 1,3)   = 'CHEM_O3_CONC'
  ItemLName(1,3)   = 'O3 concentration'
  ItemUnits(1,3)   = 'mol/mol'
  ItemSpDim(1,3)   = 3
  ItemName (1,4)   = 'CHEM_PAN_CONC'
  ItemLName(1,4)   = 'NO concentration'
  ItemUnits(1,4)   = 'mol/mol'
  ItemSpDim(1,4)   = 3
  ItemName (1,5)   = 'CHEM_ISOP_CONC'
  ItemLName(1,5)   = 'ISOP concentration'
  ItemUnits(1,5)   = 'mol/mol'
  ItemSpDim(1,5)   = 3
  ItemName (1,6)   = 'CHEM_HNO3_CONC'
  ItemLName(1,6)   = 'HNO3 concentration'
  ItemUnits(1,6)   = 'mol/mol'
  ItemSpDim(1,6)   = 3

  CollName (2  )   = 'Hourly'
  FreqYmd  (2  )   = 000000
  FreqHms  (2  )   = 010000
  ItemName (2,1)   = 'EMIS_NO_BIOM'
  ItemLName(2,1)   = 'NO biomass emissions'
  ItemUnits(2,1)   = 'kg/m2'
  ItemSpDim(2,1)   = 2
  ItemName (2,2)   = 'EMIS_NO_BIOF'
  ItemLName(2,2)   = 'NO biofuel emissions'
  ItemUnits(2,2)   = 'kg/m2'
  ItemSpDim(2,2)   = 2
  ItemName (2,3)   = 'EMIS_NO_ANTH'
  ItemLName(2,3)   = 'NO anthropogenic emissions'
  ItemUnits(2,3)   = 'mol/mol'
  ItemSpDim(2,3)   = 2
  ItemName (2,4)   = 'EMIS_NO_SHIP'
  ItemLName(2,4)   = 'NO ship emissions'
  ItemUnits(2,4)   = 'kg/m2'
  ItemSpDim(2,4)   = 2
  ItemName (2,5)   = 'EMIS_NO_AIRC'
  ItemLName(2,5)   = 'NO aircraft emissions'
  ItemUnits(2,5)   = 'kg/m2'
  ItemSpDim(2,5)   = 2

  CollName (3  ) = 'Daily'
  FreqYmd  (3  ) = 000001
  FreqHms  (3  ) = 000000
  ItemName (3,1) = 'MET_PEDGE'
  ItemLName(3,1) = 'Pressure edges'
  ItemUnits(3,1) = 'hPa'
  ItemSpDim(3,1) = 2
  ItemName (3,2) = 'MET_THETA'
  ItemLName(3,2) = 'Potential temperature'
  ItemUnits(3,2) = 'K'
  ItemSpDim(3,2) = 3
  ItemName (3,3) = 'MET_UWND'
  ItemLName(3,3) = 'E-W wind component'
  ItemUnits(3,3) = 'm/s'
  ItemSpDim(3,3) = 3
  ItemName (3,4) = 'CHEM_NO2_JVALUE'
  ItemLName(3,4) = 'Photolyis rate for NO2'
  ItemUnits(3,4) = 's-1'
  ItemSpDim(3,4) = 3
  ItemName (3,5) = 'DRYD_O3_FLUXCHEM'
  ItemLName(3,5) = 'Drydep flux of O3, applied in chemistry'
  ItemUnits(3,5) = 's-1'
  ItemSpDim(3,5) = 2


!  ErrMsg ='CHEM_SPC_*'
!  N = INDEX( ErrMsg, '*' )
!
!  ErrMsg = ErrMsg(1:N-1) // 'NO' // ErrMsg(N+1:)
!  print*, TRIM( ErrMsg )
!  Stop

!  FileName = 'HISTORY_GCC.rc'
  FileName = 'HISTORY.rc'


  ! Get a list of the collection names in the file
  CALL History_ReadCollectionNames( am_I_root, FileName, RC )
  IF ( RC /= GC_SUCCESS ) THEN
     ErrMsg = 'Error parsing input file ' // TRIM( FileName ) //'!'
     CALL GC_Error( ErrMsg, RC, ThisLoc )
     STOP
  ENDIF

  ! Get a list of the collection names in the file
  CALL History_ReadCollectionData( am_I_root, FileName, RC )
  IF ( RC /= GC_SUCCESS ) THEN
     ErrMsg = 'Error parsing input file ' // TRIM( FileName ) //'!'
     CALL GC_Error( ErrMsg, RC, ThisLoc )
     STOP
  ENDIF


! Comment out for now
!
!  ! Create the collection list
!  CALL Create_HistContainerList  ( am_I_Root, ContainerList, RC )
!
  ! Print the list
  CALL MetaHistContainer_Print   ( am_I_Root, ContainerList, RC )
!
!  ! Cleanup the list
!  CALL MetaHistContainer_Destroy ( am_I_Root, ContainerList, RC )


  CALL History_Cleanup( am_I_Root, RC )
  IF ( RC /= GC_SUCCESS ) THEN
     ErrMsg = 'Could not deallocate module variables!'
     CALL GC_Error( ErrMsg, RC, ThisLoc )
     STOP
  ENDIF


CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Create_List_Of_Collections 
!
! !DESCRIPTION: Creates a linked list for the various diagnostic collections.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Create_HistContainerList( am_I_Root, ContainerList, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE Error_Mod,          ONLY : GC_Error
    USE History_Params_Mod
    USE HistItem_Mod
    USE MetaHistItem_Mod
    USE HistContainer_Mod
    USE MetaHistContainer_Mod
!
! !INPUT PARAMETERS:
!
    LOGICAL,                 INTENT(IN) :: am_I_ROot
!
! !INPUT/OUTPUT PARAMETERS: 
!   
    TYPE(MetaHistContainer), POINTER    :: ContainerList  
!
! !OUTPUT PARAMETERS:
!
    INTEGER,                 INTENT(OUT) :: RC
!
! !REMARKS:
!  Create the list of collections.  Fake it as if we are reading
!  the collections and their fields from a file like HISTORY.rc.
!  In the real code, we can just create each collection and item as they
!  are read from an e.g. HISTORY.rc or similar input file.
!
! !REVISION HISTORY:
!  12 Jun 2017 - R. Yantosca - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER                 :: ArchivalYmd
    INTEGER                 :: ArchivalHms
    INTEGER                 :: C
    INTEGER                 :: I
    INTEGER                 :: ISpaceDim
    INTEGER                 :: ILev

    ! Strings
    CHARACTER(LEN=HMAXCHAR) :: Name
    CHARACTER(LEN=HMAXCHAR) :: Prefix
    CHARACTER(LEN=HMAXCHAR) :: Reference
    CHARACTER(LEN=HMAXCHAR) :: Title
    CHARACTER(LEN=HMAXCHAR) :: IName
    CHARACTER(LEN=HMAXCHAR) :: ILongName
    CHARACTER(LEN=HMAXCHAR) :: IUnits

    ! Pointers
    TYPE(HistItem     ), POINTER :: Item
    TYPE(HistContainer), POINTER :: Container

    CHARACTER(LEN=HMAXCHAR) :: ErrMsg
    CHARACTER(LEN=HMAXCHAR) :: ThisLoc

    !=======================================================================
    ! Create_List_Of_Collections begins here!
    !=======================================================================

    ! For error msgs
    ErrMsg  = ''
    ThisLoc = ' -> at History_Driver (in History/history_driver.F90)'

    ! Free pointers
    Item           => NULL()
    Container      => NULL()
    ContainerList  => NULL()
   
    ! The number of collections and fields per collection would come
    ! from the history input file.  For now we are faking this by looping
    ! over preset array values just to test the code.
    DO C = 1, SIZE( CollName )

       ! Exit on the first null string
       IF ( LEN_TRIM( CollName(C) ) == 0 ) EXIT

       ! Define inputs (normally these would come form an input file
       Name        = CollName(C)
       ArchivalYmd = FreqYmd(C)
       ArchivalHms = FreqHms(C)
       Prefix      = 'GEOSChem_History.' // TRIM( Name ) // '.'
       Reference   = 'www.geos-chem.org; wiki.geos-chem.org'
       Title       = 'GEOS-Chem ' // TRIM( Name ) // ' history archival'      

       !===================================================================
       ! Create the collection object
       !===================================================================
       CALL HistContainer_Create( am_I_Root   = am_I_Root,       &
                                  Container   = Container,       &
                                  Name        = Name,            &
                                  Id          = C,               &
                                  ArchivalYmd = ArchivalYmd,     &
                                  ArchivalHms = ArchivalHms,     &
                                  FilePrefix  = Prefix,          &
                                  Conventions = 'COARDS',        &
                                  NcFormat    = 'netCDF-4',      &
                                  Reference   = Reference,       &
                                  Title       = Title,           &
                                  RC          = RC              )
       
       ! Trap error
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Could not create HISTORY CONTAINER ' // TRIM( Name )
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          STOP
       ENDIF

       !===================================================================
       ! Create a history item object for each diagnostic quantity
       ! that belongs to this collection.
       !===================================================================
       
       ! Normally we would read in the items per collection
       ! from a file, but this loop is meant to "fake" that for now
       DO I = 1, SIZE( ItemName, 2 )

          ! Exit on the first null string
          IF ( LEN_TRIM( ItemName(C,I) ) == 0 ) EXIT

          IName     = ItemName (C,I)
          ILongName = ItemLName(C,I)
          IUnits    = ItemUnits(C,I)

          ! Hack number of levels for now
          IF ( ISpaceDim == 3 ) THEN
             ILev = 72 
          ELSE
             ILev = 1
          ENDIF
        
          ! Create the history item (Can add other fields later)
          CALL HistItem_Create( am_I_Root   = am_I_Root,  &
                                Item        = Item,       &
                                Id          = I,          &
                                ContainerId = C,          &
                                Name        = IName,      &
                                LongName    = ILongName,  &
                                Units       = IUnits,     &
                                SpaceDim    = ISpaceDim,  &
                                NX          = 72,         &
                                NY          = 46,         &
                                NZ          = ILev,       & 
                                RC          = RC         )
        
          ! Trap error
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = 'Could not create HISTORY ITEM ' // TRIM( Name )
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             STOP
          ENDIF

          ! Add this HISTORY ITEM to a new METAHISTORY ITEM, and
          ! add that to the METAHISTORY ITEM object belonging to
          ! the HISTORY CONTAINER.  This will represent the list of
          ! individual diagnostic quantities to be archived to 
          ! netCDF output with a specific archiving frequency.
          CALL MetaHistItem_AddNew( am_I_Root, Container%HistItems,  &
                                    Item,      RC                   )
       ENDDO

       !===================================================================
       ! Add the collection into the linked list of collections
       !===================================================================
       CALL MetaHistContainer_AddNew( am_I_Root, ContainerList,  &
                                      Container, RC             )
    ENDDO

  END SUBROUTINE Create_HistContainerList
!EOC
END PROGRAM History_Driver
