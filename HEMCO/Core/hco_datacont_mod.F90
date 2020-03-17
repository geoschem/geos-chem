!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: hco_datacont_mod.F90
!
! !DESCRIPTION: Module HCO\_DATACONT\_MOD contains routines and
! variables to handle the HEMCO data-container (DataCont) and
! correspoding list-container (ListCont) derived type.
!\\
!\\
! DataCont holds all information of an emission field, such as
! emission category, emission hierarchy, scale factors, etc.
! DataCont also contains a pointer to the source data (see
! HCO\_FILEDATA\_MOD) for more information on the file data object.
! A data-container will be created for every emission field
! specified in the HEMCO configuration file.
!\\
!\\
! The ListCont object is a derived type used to create linked lists.
! It contains a pointer to one data container (Dta) and a pointer to
! the next element of the list (NextCont). All HEMCO lists (ConfigList,
! ReadList, ListCont) are built from ListCont elements.
!\\
!\\
! DataCont consists of the following elements:
!
! \begin{itemize}
! \item cName: container name, as set in the configuration file.
! \item cID: container ID, defined by HEMCO.
! \item targetID: target ID of this container. If target ID differs
!       from the container ID, the data will be added to the content
!       of the container with cID = targetID (e.g. data of container
!       1 will be added to container 5 if it has a target ID of 5).
!       Internal use only.
! \item DctType: container type. 1 for base emissions, 2 for scale
!       factors, 3 for masks (set parameter in HCO\_ERROR\_MOD)
! \item SpcName: Species name associated with this data container, as
!       read from the configuration file. Only relevant for base
!       emission arrays.
! \item HcoID: HEMCO species ID corresponding to SpcName.
! \item ExtNr: Extension number. Extension number 0 is reserved for
!       HEMCO core, other extensions can have freely defined extensions
!       number, as specified in the configuration file. Only relevant
!       for base emissions.
! \item Cat: emission category, as set in the configuration file. Only
!       relevant for base emissions.
! \item Hier: emission hierarchy, as set in the configuration file. Only
!       relevant for base emissions.
! \item ScalID: scale factor ID, as set in the configuration file. Only
!       relevant for scale factors and masks.
! \item Oper: mathematical operator applied to scale factor. If 1, the
!       field will be multiplied (E=BxS); if -1, division is applied
!       (E=B/S); if 2, field will be squared (E=BxSxS). For masks,
!       operator 3 can be used to mirror the mask data, i.e. E=Bx(1-S).
!       Only relevant for scale factors and masks.
! \item Scal\_cID: vector of scale factor IDs associated to a base
!       emission field, as specified in the configuration file. Only
!       relevant for base emissions.
! \item Scal\_cID\_set: the Scal\_cID values read from the configuration
!       file are translated to the corresponding container IDs values
!       (the scale IDs are defined in the configuration file, container
!       IDs are automatically set by HEMCO) to optimize container
!       assignment operations. Scal\_cID\_set indicates whether or not
!       the Scal\_cID holds the container IDs or still the original
!       scale factor IDs. For internal use only.
! \item Dta: a file data object, holding information about the source
!       file, update frequency, the data arrays, etc. See
!       HCO\_FILEDATA\_MOD for more information.
! \item DtaHome: a data container only holds a pointer to a file data
!       object, i.e. it is possible that multiple containers share the
!       same file data object. the DtaHome flag is used to determine
!       whether this is the home container of this file data object. For
!       internal use only.
! \end{itemize}
!
! !INTERFACE:
!
MODULE HCO_DataCont_Mod
!
! !USES:
!
  USE HCO_TYPES_MOD
  USE HCO_Error_Mod
  USE HCO_Arr_Mod

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC  :: DataCont_Init
  PUBLIC  :: DataCont_Cleanup
  PUBLIC  :: cIDList_Create
  PUBLIC  :: cIDList_Cleanup
  PUBLIC  :: Pnt2DataCont
  PUBLIC  :: ListCont_NextCont
  PUBLIC  :: ListCont_Find
  PUBLIC  :: ListCont_Length
  PUBLIC  :: ListCont_Cleanup
!
! !REVISION HISTORY:
!  19 Dec 2013 - C. Keller: Initialization
!  01 Jul 2014 - R. Yantosca - Now use F90 free-format indentation
!  01 Jul 2014 - R. Yantosca - Cosmetic changes in ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
  ! Maximum number of scale factor fields per base field
!  INTEGER, PARAMETER,     PUBLIC :: SclMax = 10

  ! Maximum number of emission categories that can be assigned to a
  ! base field. If multiple emission categories are assigned to one
  ! field, a 'shadow' container is created for every additional
  ! emission category. A dummy scale factor of zero is applied to
  ! this shadow container, making sure that no additional emissions
  ! are created by the shadow container.
  INTEGER, PARAMETER,     PUBLIC :: CatMax = 4

  ! Fixed scale factor ID for 'dummy' scale factor of zero.
  ! Internally used to let an emission field cover multiple
  ! emission categories at once. The scale factor here must not
  ! be used in the HEMCO configuration file, otherwise HEMCO will
  ! exit with an error.
  INTEGER, PARAMETER,     PUBLIC :: ZeroScalID = 65123
!
! !PRIVATE TYPES:
!
  !-------------------------------------------------------------------------
  ! Other module variables
  !-------------------------------------------------------------------------

  ! Interface
  INTERFACE ListCont_Find
     MODULE PROCEDURE ListCont_Find_Name
     MODULE PROCEDURE ListCont_Find_ID
  END INTERFACE ListCont_Find

CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: DataCont_Init
!
! !DESCRIPTION: Subroutine DataCont\_Init initializes a new (blank) data
! container Dct.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE DataCont_Init( Dct, cID )
!
! !INPUT PARAMETERS:
!
    TYPE(DataCont),  POINTER       :: Dct
    INTEGER,         INTENT(IN)    :: cID
!
! !REVISION HISTORY:
!  19 Dec 2013 - C. Keller: Initialization
!
!EOP
!------------------------------------------------------------------------------
!BOC

    !======================================================================
    ! DataCont_Init begins here!
    !======================================================================

    ! Allocate the new container
    IF ( .NOT. ASSOCIATED( Dct) ) ALLOCATE( Dct )

    ! Nullify pointers
    Dct%Dta         => NULL()
    Dct%Scal_cID    => NULL()

    ! Set default values
    Dct%DtaHome      = -999
    Dct%DctType      = -999
    Dct%ExtNr        = 0
    Dct%cName        = ''
    Dct%spcName      = ''
    Dct%ScalID       = -999
    Dct%HcoID        = -999
    Dct%Cat          = -999
    Dct%Hier         = -999
    Dct%Oper         = 1
    Dct%levScalID1   = -1
    Dct%levScalID2   = -1
    Dct%nScalID      = 0
    Dct%Scal_cID_set = .FALSE.

    ! Assign container ID.
    ! Set default target ID to cont. ID.
    Dct%cID          = cID
    Dct%targetID     = Dct%cID

  END SUBROUTINE DataCont_Init
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: DataCont_Cleanup
!
! !DESCRIPTION: Subroutine DataCont\_Cleanup cleans up data container Dct.
! If ArrOnly is set to True, this will only cleanup the data array of the
! container but keep all meta-data.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE DataCont_Cleanup( Dct, ArrOnly )
!
! !USES:
!
    USE HCO_FILEDATA_MOD, ONLY : FileData_Cleanup
!
! !ARGUMENTS:
!
    TYPE(DataCont), POINTER               :: Dct
    LOGICAL,        INTENT(IN), OPTIONAL  :: ArrOnly
!
! !REVISION HISTORY:
!  19 Dec 2013 - C. Keller: Initialization
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: I
    LOGICAL :: DeepClean

    !======================================================================
    ! DataCont_Cleanup begins here!
    !======================================================================

    IF ( PRESENT(ArrOnly) ) THEN
       DeepClean = .NOT. ArrOnly
    ELSE
       DeepClean = .TRUE.
    ENDIF

    ! Only if associated...
    IF ( ASSOCIATED( Dct ) ) THEN

       ! Clean up FileData object. If DeepClean is true, this
       ! will entirely erase the file data object. Otherwise, only the
       ! data arrays will be removed.
       ! Note: do only if this is the home container of the file data
       ! object.
       IF ( Dct%DtaHome == 1 ) THEN
          CALL FileData_Cleanup( Dct%Dta, DeepClean )
       ENDIF

       ! Clean up data container if DeepClean option is enabled.
       IF ( DeepClean ) THEN
          Dct%Dta => NULL()
          IF(ASSOCIATED(Dct%Scal_cID)) DEALLOCATE(Dct%Scal_cID)
          DEALLOCATE ( Dct )
       ENDIF

    ENDIF

  END SUBROUTINE DataCont_Cleanup
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ListCont_Cleanup
!
! !DESCRIPTION: Subroutine ListCont\_Cleanup cleans up list List
! The corresponding data container (LstCont%Dct) is also removed if
! RemoveDct is set to true.
!\\
! !INTERFACE:
!
  SUBROUTINE ListCont_Cleanup( List, RemoveDct )
!
! !INPUT PARAMETERS:
!
    TYPE(ListCont), POINTER      :: List
    LOGICAL,        INTENT(IN)   :: RemoveDct
!
! !REVISION HISTORY:
!  19 Dec 2013 - C. Keller: Initialization
!  25 Oct 2016 - R. Yantosca - Do not nullify pointers in declaration stmts
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    TYPE(ListCont), POINTER  :: TmpLct
    TYPE(ListCont), POINTER  :: NxtLct

    !======================================================================
    ! ListCont_Cleanup begins here!
    !======================================================================

    ! Initialize
    TmpLct => NULL()
    NxtLct => NULL()

    ! Walk through entire list and remove all containers
    TmpLct => List
    DO WHILE ( ASSOCIATED( TmpLct ) )

       ! Detach from list
       NxtLct          => TmpLct%NextCont
       TmpLct%NextCont => NULL()

       ! Clean up data container if flag is enabled. Otherwise, just
       ! remove pointer to container!
       IF ( RemoveDct ) THEN
          CALL DataCont_Cleanup ( TmpLct%Dct )
       ELSE
          TmpLct%Dct => NULL()
       ENDIF

       ! Remove
       DEALLOCATE ( TmpLct )

       ! Advance
       TmpLct => NxtLct
    ENDDO

    ! Nullify pointers
    TmpLct => NULL()
    NxtLct => NULL()
    List   => NULL()

  END SUBROUTINE ListCont_Cleanup
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cIDList_Create
!
! !DESCRIPTION: Subroutine cIDList\_Create creates a vector of pointers
! (cIDList) pointing to all available containers of the passed List.
! The vector index of cIDList corresponds to the container cIDs, i.e.
! cIDList(3) will point to data container with cID = 3.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE cIDList_Create( HcoState, List, RC )
!
! !USES:
!
    USE HCO_STATE_MOD, ONLY : HCO_State
!
! !ARGUMENTS:
!
    TYPE(HCO_State), POINTER       :: HcoState
    TYPE(ListCont),  POINTER       :: List
    INTEGER,         INTENT(INOUT) :: RC
!
! !REVISION HISTORY:
!  24 Aug 2012 - C. Keller - Initial Version
!  25 Oct 2016 - R. Yantosca - Do not nullify pointers in declaration stmts
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER                   :: II
    TYPE(ListCont), POINTER   :: TmpLct
    LOGICAL                   :: verbose
    CHARACTER(LEN=255)        :: MSG

    !======================================================================
    ! cIDList_Create begins here
    !======================================================================

    ! Initialize
    TmpLct => NULL()

    ! Enter
    CALL HCO_ENTER( HcoState%Config%Err, 'cIDList_Create (hco_datacont_mod.F)', RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Set verbose flag
    verbose = HCO_IsVerb ( HcoState%Config%Err, 3 )

    ! Set # of data container in list
    HcoState%nnDataCont = ListCont_Length( List )

    ! Eventually cleanup the list
    IF ( ASSOCIATED ( HcoState%cIDList ) ) THEN
       DO II = 1, HcoState%nnDataCont
          HcoState%cIDList(II)%PNT => NULL()
       ENDDO
       DEALLOCATE ( HcoState%cIDList )
    ENDIF

    ! Leave if no emission fields defined
    IF ( HcoState%nnDataCont == 0 ) THEN
       IF ( verbose ) THEN
          WRITE(MSG,*) 'No emission fields defined!'
          CALL HCO_MSG(HcoState%Config%Err,MSG)
       ENDIF
       RC = HCO_SUCCESS
       RETURN
    ENDIF

    ! verbose
    IF ( verbose ) THEN
       WRITE(MSG,*) 'Create cID list: # of fields: ', HcoState%nnDataCont
       CALL HCO_MSG(HcoState%Config%Err,MSG)
    ENDIF

    ! Allocate IDList
    ALLOCATE ( HcoState%cIDList(HcoState%nnDataCont) )

    ! Now set the quicklist pointers
    IILOOP: DO II = 1, HcoState%nnDataCont

       ! Nullify pointer first
       HcoState%cIDList(II)%PNT => NULL()

       ! Set working container to head of emission fields linked list
       TmpLct => List

       DO WHILE ( ASSOCIATED ( TmpLct ) )

          ! Ignore deallocated fields
          IF ( .NOT. ASSOCIATED(TmpLct%Dct)) THEN
             TmpLct => TmpLct%NextCont
             CYCLE
          ENDIF

          ! Check if current field is the one with the correct FID
          IF ( TmpLct%Dct%cID == II ) THEN

             ! Set pointer to emission field
             HcoState%cIDList(II)%PNT => TmpLct%Dct

             ! Advance in loop
             CYCLE IILOOP
          ENDIF

          ! Advance
          TmpLct => TmpLct%NextCont
       ENDDO

    ENDDO IILOOP

    ! Cleanup and leave w/ success
    TmpLct => NULL()
    CALL HCO_LEAVE ( HcoState%Config%Err, RC )

  END SUBROUTINE cIDList_Create
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cIDList_Cleanup
!
! !DESCRIPTION: Subroutine cIDList\_Cleanup cleans up cIDList.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE cIDList_Cleanup( HcoState )
!
! !USES:
!
    USE HCO_STATE_MOD, ONLY : HCO_State
!
! !ARGUMENTS:
!
    TYPE(HCO_State), POINTER       :: HcoState
!
! !REVISION HISTORY:
!  24 Aug 2012 - C. Keller - Initial Version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: I

    !======================================================================
    ! cIDList_Cleanup begins here
    !======================================================================

    ! Remove links to all containers
    IF ( ASSOCIATED ( HcoState%cIDList ) ) THEN
      DO I = 1, HcoState%nnDataCont
        HcoState%cIDList(I)%PNT => NULL()
      ENDDO
      DEALLOCATE( HcoState%cIDList )
    ENDIF
    HcoState%cIDList    => NULL()
    HcoState%nnDataCont =  0

  END SUBROUTINE cIDList_Cleanup
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Pnt2DataCont
!
! !DESCRIPTION: Subroutine Pnt2DataCont returns the data container Dct
! with container ID cID.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Pnt2DataCont( HcoState, cID, Dct, RC )
!
! !USES:
!
    USE HCO_STATE_MOD, ONLY : HCO_State
!
! !INPUT PARAMETERS:
!
    TYPE(HCO_State), POINTER       :: HcoState
    INTEGER,        INTENT(IN)     :: cID
    TYPE(DataCont), POINTER        :: Dct
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(INOUT)  :: RC
!
! !REVISION HISTORY:
!  11 Apr 2012 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    CHARACTER(LEN=255)  :: MSG, LOC

    !======================================================================
    ! Pnt2DataCont begins here!
    !======================================================================

    ! Enter
    LOC = 'Pnt2DataCont (HCO_DATACONT_MOD.F90)'

    ! Check input
    IF ( cID > HcoState%nnDataCont ) THEN
       MSG = 'cID higher than number of containers'
       CALL HCO_ERROR ( HcoState%Config%Err, MSG, RC, THISLOC=LOC)
       RETURN
    ENDIF

    ! Set pointer to container w/ ID cID
    Dct => HcoState%cIDList(cID)%PNT

    ! Check if data container allocated
    IF ( .NOT. ASSOCIATED( Dct ) ) THEN
       MSG = 'Data container is not associated!'
       CALL HCO_ERROR ( HcoState%Config%Err, MSG, RC, THISLOC=LOC)
       RETURN
    ENDIF

    ! Leave
    RC = HCO_SUCCESS

  END SUBROUTINE Pnt2DataCont
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ListCont_NextCont
!
! !DESCRIPTION: Routine ListCont\_NextCont returns container Lct from
! data list List. This is the generic routine for cycling through
! the data container lists.
!\\
!\\
! If Lct is empty (i.e. NULL), the first container of List is returned.
! If Lct already points to a list container, the pointer is advanced
! to the next container in that list (Lct%NextCont). The return flag
! FLAG is set to HCO\_SUCCESS if the return container Lct is defined,
! and to HCO\_FAIL otherwise.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE ListCont_NextCont( List, Lct, FLAG )
!
! !INPUT PARAMETERS:
!
    TYPE(ListCont), POINTER       :: List
    TYPE(ListCont), POINTER       :: Lct
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(INOUT) :: FLAG
!
! !REVISION HISTORY:
!  11 Apr 2012 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC

    !======================================================================
    ! ListCont_NextCont begins here!
    !======================================================================

    ! Point to head of List if passed container pointer is not yet defined.
    IF ( .NOT. ASSOCIATED ( Lct ) ) THEN
       Lct => List

    ! Otherwise, just point to the next container in list
    ELSE
       Lct => Lct%NextCont
    ENDIF

    ! Set return flag
    IF ( .NOT. ASSOCIATED ( Lct ) ) THEN
       FLAG = HCO_FAIL
    ELSE
       FLAG = HCO_SUCCESS
    ENDIF

  END SUBROUTINE ListCont_NextCont
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ListCont_Find_Name
!
! !DESCRIPTION: Subroutine ListCont\_Find\_Name searches for (data)
! container name NME in list List and returns a pointer pointing
! to this container (Lct).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE ListCont_Find_Name( List, NME, FOUND, Lct )
!
! !ARGUMENTS:
!
    TYPE(ListCont),   POINTER               :: List  ! List to be searched
    CHARACTER(LEN=*), INTENT(IN )           :: NME   ! Container name
    LOGICAL,          INTENT(OUT)           :: FOUND ! Container found?
    TYPE(ListCont),   POINTER, OPTIONAL     :: Lct   ! matched list container
!
! !REVISION HISTORY:
!  04 Dec 2012 - C. Keller: Initialization
!  25 Oct 2016 - R. Yantosca - Do not nullify pointers in declaration stmts
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL ARGUMENTS:
!
    TYPE(ListCont), POINTER :: TmpLct

    !======================================================================
    ! ListCont_Find_Name begins here!
    !======================================================================

    ! Initialize
    TmpLct => NULL()
    FOUND  = .FALSE.

    ! Error trap
    IF ( .NOT. ASSOCIATED(List) ) RETURN

    ! Make CurrCnt point to first element of the EMISSIONS linked list
    TmpLct => List

    ! Loop over EMISSIONS linked list
    DO WHILE ( ASSOCIATED ( TmpLct ) )

       ! Eventually skip over empty data containers
       IF ( .NOT. ASSOCIATED(TmpLct%Dct) ) THEN
          TmpLct => TmpLct%NextCont
          CYCLE
       ENDIF

       ! Get the current container or original ID
       ! Check if current field is the wanted one
       IF ( TRIM(TmpLct%Dct%cName) == TRIM(NME) ) THEN
          IF ( PRESENT(Lct) ) Lct => TmpLct
          FOUND = .TRUE.
          RETURN
       ENDIF

       ! Advance to next field otherwise
       TmpLct => TmpLct%NextCont
    ENDDO

    ! Cleanup
    TmpLct => NULL()

  END SUBROUTINE ListCont_Find_Name
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ListCont_Find_ID
!
! !DESCRIPTION: Subroutine ListCont\_Find\_ID searches for (data)
! container cID or ScalID (ID) in list List and returns a pointer
! pointing to this (list) container (Lct).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE ListCont_Find_ID( List, ID, IsScalID, FOUND, Lct )
!
! !INPUT PARAMETERS:
!
    TYPE(ListCont),   POINTER           :: List     ! List to be searched
    INTEGER,          INTENT(IN )       :: ID       ! cID or ScalID
    INTEGER,          INTENT(IN )       :: IsScalID ! 1=ID is ScalID;
                                                      ! else: ID is cID
!
! !OUTPUT PARAMETERS:
!
    LOGICAL,          INTENT(OUT)       :: FOUND    ! Container found?
    TYPE(ListCont),   POINTER, OPTIONAL :: Lct      ! Container w/ ID
!
! !REVISION HISTORY:
!  04 Dec 2012 - C. Keller: Initialization
!  25 Oct 2016 - R. Yantosca - Do not nullify pointers in declaration stmts
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL ARGUMENTS:
!
    TYPE(ListCont), POINTER :: TmpLct
    INTEGER                 :: thisID

    !======================================================================
    ! ListCont_Find_ID begins here!
    !======================================================================

    ! Initialize
    TmpLct => NULL()
    FOUND  = .FALSE.

    ! Error trap
    IF ( .NOT. ASSOCIATED(List) ) RETURN

    ! Make TmpLct point to first element of the EMISSIONS linked list
    TmpLct => List

    ! Loop over EMISSIONS linked list
    DO WHILE ( ASSOCIATED ( TmpLct ) )

       ! Eventually skip over empty data containers
       IF ( .NOT. ASSOCIATED(TmpLct%Dct) ) THEN
          TmpLct => TmpLct%NextCont
          CYCLE
       ENDIF

       ! Get the current container or original ID
       IF ( IsScalID == 1 ) THEN
          thisID = TmpLct%Dct%scalID
       ELSE
          thisID = TmpLct%Dct%cID
       ENDIF

       ! Check if current field is the wanted one
       IF ( thisID == ID ) THEN
          IF ( PRESENT(Lct) ) Lct => TmpLct
          FOUND = .TRUE.
          RETURN
       ENDIF

       ! Advance to next field otherwise
       TmpLct => TmpLct%NextCont
    ENDDO

    ! Cleanup
    TmpLct => NULL()

  END SUBROUTINE ListCont_Find_ID
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ListCont_Length
!
! !DESCRIPTION: Subroutine ListCont\_Length returns the length of the
!  passed list.
!\\
!\\
! !INTERFACE:
!
  FUNCTION ListCont_Length ( List ) RESULT ( nnCont )
!
! !INPUT PARAMETERS:
!
    TYPE(ListCont),  POINTER       :: List
    INTEGER                        :: nnCont
!
! !REVISION HISTORY:
!  15 Feb 2016 - C. Keller: Initial version
!
!EOP
!------------------------------------------------------------------------------
!BOC
    TYPE(ListCont), POINTER  :: TmpLct

    !======================================================================
    ! ListCont_Length begins here!
    !======================================================================

    nnCont = 0
    TmpLct => List
    DO WHILE ( ASSOCIATED( TmpLct ) )
       nnCont = nnCont + 1
       TmpLct => TmpLct%NextCont
    ENDDO
    TmpLct => NULL()

  END FUNCTION ListCont_Length
!EOC
END MODULE HCO_DATACONT_MOD
!EOM
