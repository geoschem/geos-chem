!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: hco_emislist_mod.F90
!
! !DESCRIPTION: Module HCO\_EmisList\_Mod contains routines and variables
! for the HEMCO emissions list EmisList. EmisList is a sorted collection
! of all data containers needed for emission calculation. The containers
! are sorted by data type, species, emission category, and emission
! hierarchy (in this order).
!\\
!\\
! !INTERFACE:
!
MODULE HCO_EMISLIST_MOD
!
! !USES:
!
  USE HCO_ERROR_MOD
  USE HCO_TYPES_MOD
  USE HCO_STATE_MOD,    ONLY : HCO_State

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC  :: HCO_GetPtr
  PUBLIC  :: EmisList_Pass
!  PUBLIC  :: EmisList_Update
!
! !PRIVATE MEMBER FUNCTIONS:
!
  PRIVATE :: EmisList_Add
  PRIVATE :: Add2EmisList
!
! !REVISION HISTORY:
!  04 Dec 2012 - C. Keller   - Initialization
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !PRIVATE TYPES:
!
  INTERFACE HCO_GetPtr
     MODULE PROCEDURE HCO_GetPtr_2D
     MODULE PROCEDURE HCO_GetPtr_3D
  END INTERFACE HCO_GetPtr

CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: EmisList_Add
!
! !DESCRIPTION: Subroutine EmisList\_Add adds the passed data container
! Dct to EmisList. Within EmisList, Dct becomes placed with
! increasing data type, species ID, category and hierarchy (in this order).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE EmisList_Add( Dct, HcoState, RC )
!
! !USES:
!
    USE HCO_DATACONT_MOD,  ONLY : ListCont_Find
    USE HCO_LOGFILE_MOD,   ONLY : HCO_PrintDataCont
!
! !INPUT PARAMETERS:
!
    TYPE(DataCont),    POINTER       :: Dct        ! Data cont.
    TYPE(HCO_State),   POINTER       :: HcoState   ! HEMCO state
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,           INTENT(INOUT) :: RC         ! Return code
!
! !REVISION HISTORY:
!  04 Dec 2012 - C. Keller - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL ARGUMENTS:
!
    TYPE(ListCont), POINTER                 :: Lct
    LOGICAL                                 :: FOUND, VERBOSE, NEW
    CHARACTER(LEN=255)                      :: MSG
    CHARACTER(LEN= 31)                      :: TempRes

    !======================================================================
    ! EmisList_Add begins here!
    !======================================================================

    ! Enter
    CALL HCO_ENTER ( HcoState%Config%Err, &
                     'EmisList_Add (HCO_EMISLL_MOD.F90)', RC )
    IF(RC /= HCO_SUCCESS) RETURN

    ! Set verbose flag
    VERBOSE = HCO_IsVerb( HcoState%Config%Err, 2 )

    ! Init
    Lct => NULL()

    ! ----------------------------------------------------------------
    ! Nothing to do if it's not a new container, i.e. if container
    ! already exists in EmisList.
    ! ----------------------------------------------------------------
    CALL ListCont_Find ( HcoState%EmisList, Dct%cID, 0, FOUND, Lct )
    IF ( FOUND ) THEN
       CALL HCO_LEAVE ( HcoState%Config%Err, RC )
       RETURN
    ENDIF

    ! ----------------------------------------------------------------
    ! Define new list container. This container points to the passed
    ! data container.
    ! ----------------------------------------------------------------
    ALLOCATE ( Lct )
    Lct%NextCont => NULL()
    IF ( .not. ASSOCIATED( Dct ) ) THEN
       PRINT*, '#### DCT is not associated!'
    ENDIF
    Lct%Dct      => Dct

    ! ----------------------------------------------------------------
    ! Add the new container to EmisList. The container will be placed
    ! according to data type, species ID, hierarchy, and category.
    ! ----------------------------------------------------------------
    CALL Add2EmisList ( HcoState, Lct, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! ----------------------------------------------------------------
    ! Verbose mode
    ! ----------------------------------------------------------------
    IF ( VERBOSE ) THEN
       MSG = 'Container added to EmisList:'
       CALL HCO_MSG(HcoState%Config%Err,MSG)
       CALL HCO_PrintDataCont( HcoState, Lct%Dct, 3 )
    ENDIF

    ! Leave w/ success
    CALL HCO_LEAVE ( HcoState%Config%Err, RC )

  END SUBROUTINE EmisList_Add
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Add2EmisList
!
! !DESCRIPTION: Subroutine Add2EmisList adds list container Lct to
! EmisList. Base emission fields (Data type = 1) are sorted based on
! species ID, category and hierarchy (for fields of same category). Scale
! fields and masks are added to the end of EmisList.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Add2EmisList( HcoState, Lct, RC )
!
! !INPUT PARAMETERS:
!
    TYPE(HCO_State), POINTER       :: HcoState   ! HEMCO state
    TYPE(ListCont),  POINTER       :: Lct
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,         INTENT(INOUT) :: RC
!
! !REVISION HISTORY:
!  06 Dec 2012 - C. Keller - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER                   :: NEWCAT, NEWHIR, NEWSPC
    CHARACTER(LEN=255)        :: MSG

    ! Pointers
    TYPE(ListCont), POINTER   :: TmpLct => NULL()

    !======================================================================
    ! Add2EmisList begins here!
    !======================================================================

    ! Enter
    CALL HCO_ENTER ( HcoState%Config%Err, 'Add2EmisList', RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Update number of containers in EmisList
    HcoState%nnEmisCont = HcoState%nnEmisCont + 1

    ! Flag the content of this container as being used in EmisList
    Lct%Dct%Dta%IsInList = .TRUE.

    ! If this is the first container, we can simply place it at the
    ! beginning of the list.
    IF ( HcoState%nnEmisCont == 1 ) THEN
       HcoState%EmisList => Lct
       CALL HCO_LEAVE( HcoState%Config%Err, RC )
       RETURN ! Leave routine
    ENDIF

    ! Special case where the linked list consists of scale factors
    ! only: In this case, we can place the new container at the
    ! beginning no matter of its content!
    IF ( HcoState%EmisList%Dct%DctType /= HCO_DCTTYPE_BASE ) THEN
       Lct%NextCont      => HcoState%EmisList
       HcoState%EmisList => Lct
       CALL HCO_LEAVE( HcoState%Config%Err, RC )
       RETURN
    ENDIF

    ! Get field species ID, category and priority of the new container
    NEWSPC  = Lct%Dct%HcoID
    NEWCAT  = Lct%Dct%Cat
    NEWHIR  = Lct%Dct%Hier

    ! Containers are listed with increasing species ID. If the current
    ! container has lower speciesID than the first container, just add
    ! it at the beginning of the list.
    IF ( (NEWSPC > 0) .AND. (NEWSPC < HcoState%EmisList%Dct%HcoID) ) THEN
       Lct%NextCont      => HcoState%EmisList
       HcoState%EmisList => Lct
       CALL HCO_LEAVE( HcoState%Config%Err, RC )
       RETURN
    ENDIF

    ! In case that the current container has the same species ID
    ! as the first container in the list: If this container has
    ! lower category, or same category and lower hierarchy, place
    ! it before the first container in the list:
    IF ( NEWSPC == HcoState%EmisList%Dct%HcoID ) THEN
       IF ( (HcoState%EmisList%Dct%Cat  >  NEWCAT) .OR.  &
            (HcoState%EmisList%Dct%Cat  == NEWCAT  .AND. &
             HcoState%EmisList%Dct%Hier >  NEWHIR)        ) THEN
          Lct%NextCont      => HcoState%EmisList
          HcoState%EmisList => Lct
          CALL HCO_LEAVE( HcoState%Config%Err, RC )
          RETURN
       ENDIF
    ENDIF

    ! TmpLct is the temporary working pointer, looping through
    ! the entire EmisList until the correct place for the new
    ! container is found.
    TmpLct => HcoState%EmisList

    ! If the new container contains base data (i.e. data type is 1)
    ! we have to move the TmpLct pointer to the position where the
    ! next container is also a base container and one of the
    ! following: (a) the first container with the same species ID
    ! as the new container; (b) a container with higher species ID;
    ! (c) scale factors. From there, we can determine where to place
    ! the container exactly.
    IF ( Lct%Dct%DctType == HCO_DCTTYPE_BASE ) THEN

       ! Loop over list
       DO WHILE ( ASSOCIATED ( TmpLct%NextCont ) )

          ! Check if next container's species ID is higher or if it's a
          ! scale factor, in which case we have to exit.
          IF ( TmpLct%NextCont%Dct%HcoID   >  NEWSPC          .OR. &
               TmpLct%NextCont%Dct%DctType /= HCO_DCTTYPE_BASE      ) THEN
             EXIT
          ENDIF

          ! Check if next container has the same species ID but a
          ! higher category or the same category but higher hierarchy,
          ! in which case we have to exit.
          IF ( TmpLct%NextCont%Dct%HcoID == NEWSPC ) THEN
             IF ( TmpLct%NextCont%Dct%Cat > NEWCAT ) THEN
                EXIT
             ENDIF
             IF ( TmpLct%NextCont%Dct%Cat  == NEWCAT .AND. &
                  TmpLct%NextCont%Dct%Hier >  NEWHIR        ) THEN
                EXIT
             ENDIF
          ENDIF

          ! Advance in list if none of the above checks was true.
          TmpLct => TmpLct%NextCont
       ENDDO

    ! Scale factors and masks are collected at the end of the list.
    ! Hence, make TmpLct pointer point to the last container w/ base
    ! emissions (or the last container in the list).
    ELSE

       ! Loop over list
       DO WHILE ( ASSOCIATED ( TmpLct%NextCont ) )

          ! Check if next container is scale factor
          IF ( TmpLct%NextCont%Dct%DctType /= HCO_DCTTYPE_BASE ) EXIT

          ! Advance in list
          TmpLct => TmpLct%NextCont
       ENDDO

    ENDIF

    ! Add new container AFTER current one
    Lct%NextCont    => TmpLct%NextCont
    TmpLct%NextCont => Lct

    ! Cleanup and leave
    TmpLct => NULL()
    CALL HCO_LEAVE( HcoState%Config%Err, RC )

  END SUBROUTINE Add2EmisList
!EOC
!!------------------------------------------------------------------------------
!!                  Harvard-NASA Emissions Component (HEMCO)                   !
!!------------------------------------------------------------------------------
!!BOP
!!
!! !IROUTINE: EmisList_Update
!!
!! !DESCRIPTION: Subroutine EmisList\_Update makes sure that all containers
!! of the reading list ReadList are correctly referenced in emissions list
!! EmisList. If a container of ReadList does not yet have a corresponding
!! container in EmisList, such a container is created. Also, additive data
!! arrays (i.e. targetID different than container ID) are added to their
!! target array during this call.
!!\\
!!\\
!! !INTERFACE:
!!
!  SUBROUTINE EmisList_Update ( HcoState, ReadList, RC )
!!
!! !USES:
!!
!    USE HCO_FILEDATA_MOD, ONLY : FileData_ArrIsDefined
!!
!! !INPUT PARAMETERS:
!!
!    TYPE(HCO_State), POINTER       :: HcoState   ! Hemco state object
!    TYPE(ListCont),  POINTER       :: ReadList   ! reading list
!!
!! !INPUT/OUTPUT PARAMETERS:
!!
!    INTEGER,         INTENT(INOUT) :: RC         ! Return code
!!
!! !REVISION HISTORY:
!!  20 Apr 2013 - C. Keller - Initial version
!!  See https://github.com/geoschem/geos-chem for complete history
!!EOP
!!------------------------------------------------------------------------------
!!BOC
!!
!! !LOCAL VARIABLES:
!!
!    ! Scalars
!    INTEGER                     :: cID, iScalID
!    INTEGER                     :: I
!    CHARACTER(LEN=31)           :: ScalName
!
!    ! Pointers
!    TYPE(ListCont), POINTER     :: TmpLct => NULL()
!
!    ! ================================================================
!    ! EmisList_Update begins here
!    ! ================================================================
!
!    ! Enter
!    CALL HCO_ENTER ( HcoState%Config%Err, 'EmisList_Update', RC )
!    IF ( RC /= HCO_SUCCESS ) RETURN
!
!    ! Loop over all containers in ReadList
!    TmpLct => ReadList
!    DO WHILE ( ASSOCIATED( TmpLct ) )
!
!       ! only if array is defined...
!       IF ( FileData_ArrIsDefined(TmpLct%Dct%Dta) ) THEN
!
!          ! Pass container to EmisList
!          CALL EmisList_Pass( HcoState, TmpLct, RC )
!          IF ( RC /= HCO_SUCCESS ) RETURN
!       ENDIF
!
!       ! Advance to next container in ReadList
!       TmpLct => TmpLct%NextCont
!    ENDDO
!
!    ! Leave w/ success
!    CALL HCO_LEAVE ( HcoState%Config%Err, RC )
!
!  END SUBROUTINE EmisList_Update
!!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: EmisList_Pass
!
! !DESCRIPTION: Subroutine EmisList\_Pass passes (the ReadList)
! container Lct to EmisList. This routine mostly checks for
! additive arrays, i.e. if arrays from multiple containers have
! to be added together prior to emission calculation (e.g. sectoral
! data).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE EmisList_Pass( HcoState, Lct, RC )
!
! !USES:
!
    USE HCO_DATACONT_MOD, ONLY : ListCont_Find
    USE HCO_FILEDATA_MOD, ONLY : FileData_ArrCheck
!
! !INPUT PARAMETERS:
!
    TYPE(HCO_State),  POINTER       :: HcoState
    TYPE(ListCont),   POINTER       :: Lct        ! list container
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(INOUT) :: RC         ! Success
!
! !REVISION HISTORY:
!  28 Mar 2013 - C. Keller - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Pointers
    TYPE(ListCont),  POINTER  :: TargetLct

    ! Scalars
    INTEGER                   :: I, J, L, T
    LOGICAL                   :: FOUND, verb, Add
    CHARACTER(LEN=255)        :: MSG

    ! ================================================================
    ! EmisList_Pass begins here
    ! ================================================================

    ! Enter
    CALL HCO_ENTER ( HcoState%Config%Err, &
                     'EmisList_Pass (hco_emislist_mod.F90)', RC )
    IF(RC /= HCO_SUCCESS) RETURN

    ! Init
    TargetLct => NULL()

    ! Verbose mode
    verb = HCO_IsVerb( HcoState%Config%Err, 2 )

    ! Initialize Add flag. This fill only be set to FALSE
    ! if the data of the current container is added to the data of
    ! an existing container in EmisList instead of adding the container
    ! alltogether.
    Add = .TRUE.

    ! ----------------------------------------------------------------
    ! Add data arrays if required
    ! ----------------------------------------------------------------

    ! The target ID of a container denotes the ID (cID) of the
    ! container to which the data shall be added. For example, if
    ! container 1 has a target ID of 5, its content will be added to
    ! container 5.
    ! Usually, the targetID is equal to cID and we don't have to do
    ! anything. If tID /= cID, however, the array is added to the
    ! array of the specified target container and removed afterwards
    ! from the original container.
    ! Note: arrays can only be added to each other if they are for the
    ! same species, have same dimensions, update frequencies, scale
    ! factors, categories, hierarchies, data types, etc.
    ! Note2: in an ESMF environment, this option is disabled (targetID
    ! is always equal to cID).
    IF ( Lct%Dct%targetID /= Lct%Dct%cID ) THEN

       ! TargetLct points to the container holding the target array
       CALL ListCont_Find( HcoState%EmisList, Lct%Dct%targetID, &
                           0, FOUND, TargetLct )
       IF ( .NOT. FOUND ) THEN
          MSG = 'Cannot add emissions to target array: error in ' // &
               TRIM(Lct%Dct%cName)
          CALL HCO_ERROR( HcoState%Config%Err, MSG, RC )
          RETURN
       ENDIF

       ! Do not add data if the current data container is not the
       ! 'home' container for the data object Dta. Dta may be used
       ! by multiple containers, and only the home container should
       ! modify its content!
       IF ( Lct%Dct%DtaHome /= 1 ) THEN

          ! Verbose mode
          IF ( verb ) THEN
             WRITE(MSG,*) 'Do not add data of ', TRIM(Lct%Dct%cName), &
                  ' to ', TRIM(TargetLct%Dct%cName), ' because this', &
                  ' is not the file data home container!'
             CALL HCO_MSG(HcoState%Config%Err,MSG)
          ENDIF

       ! Similarly, do not add data to target container if the target
       ! container is being shared by multiple containers.
       ELSEIF ( TargetLct%Dct%Dta%DoShare ) THEN

          ! Verbose mode
          IF ( verb ) THEN
             WRITE(MSG,*) 'Do not add data of ', TRIM(Lct%Dct%cName), &
                  ' to ', TRIM(TargetLct%Dct%cName), ' because the', &
                  ' target is being shared with other fields!'
             CALL HCO_MSG(HcoState%Config%Err,MSG)
          ENDIF

       ELSE

          ! Check extension number
          IF ( Lct%Dct%ExtNr /= TargetLct%Dct%ExtNr ) THEN
             MSG = 'Wrong ext. number: ' // TRIM(Lct%Dct%cName)
             CALL HCO_ERROR( HcoState%Config%Err, MSG, RC )
             RETURN
          ENDIF

          ! Check data type
          IF ( Lct%Dct%DctType /= TargetLct%Dct%DctType ) THEN
             MSG = 'Wrong data type: ' // TRIM(Lct%Dct%cName)
             CALL HCO_ERROR( HcoState%Config%Err, MSG, RC )
             RETURN
          ENDIF

          ! Check species ID
          IF ( Lct%Dct%HcoID /= TargetLct%Dct%HcoID ) THEN
             MSG = 'Wrong species ID: ' // TRIM(Lct%Dct%cName)
             CALL HCO_ERROR( HcoState%Config%Err, MSG, RC )
             RETURN
          ENDIF

          ! Check for array dimensions
          IF ( Lct%Dct%Dta%SpaceDim /= TargetLct%Dct%Dta%SpaceDim ) THEN
             MSG = 'Wrong space dimension: ' // TRIM(Lct%Dct%cName)
             CALL HCO_ERROR( HcoState%Config%Err, MSG, RC )
             RETURN
          ENDIF
          IF ( Lct%Dct%Dta%nt /= TargetLct%Dct%Dta%nt ) THEN
             MSG = 'Wrong time dim: ' // TRIM(Lct%Dct%cName)
             CALL HCO_ERROR( HcoState%Config%Err, MSG, RC )
             RETURN
          ENDIF
          IF ( Lct%Dct%Dta%SpaceDim <= 2) THEN
             I = SIZE(Lct%Dct%Dta%V2(1)%Val,1)
             J = SIZE(Lct%Dct%Dta%V2(1)%Val,2)
             CALL FileData_ArrCheck( HcoState%Config, &
                                     TargetLct%Dct%Dta, I, J, &
                                     Lct%Dct%Dta%nt, RC )
             IF ( RC /= 0 ) THEN
                MSG = 'Wrong 2D array: ' // TRIM(Lct%Dct%cName)
                CALL HCO_MSG(HcoState%Config%Err,MSG)
                RETURN
             ENDIF
          ELSE
             I = SIZE(Lct%Dct%Dta%V3(1)%Val,1)
             J = SIZE(Lct%Dct%Dta%V3(1)%Val,2)
             L = SIZE(Lct%Dct%Dta%V3(1)%Val,3)
             CALL FileData_ArrCheck( HcoState%Config, &
                                     TargetLct%Dct%Dta, I, J, L, &
                                     Lct%Dct%Dta%nt, RC )
             IF ( RC /= 0 ) THEN
                MSG = 'Wrong 3D array: ' // TRIM(Lct%Dct%cName)
                CALL HCO_MSG(HcoState%Config%Err,MSG)
                RETURN
             ENDIF
          ENDIF

          ! Check operator
          IF ( Lct%Dct%Oper /= TargetLct%Dct%Oper ) THEN
             MSG = 'Wrong operator: ' // TRIM(Lct%Dct%cName)
             CALL HCO_ERROR( HcoState%Config%Err, MSG, RC )
             RETURN
          ENDIF

          ! Check category
          IF ( Lct%Dct%Cat /= TargetLct%Dct%Cat ) THEN
             MSG = 'Wrong category: ' // TRIM(Lct%Dct%cName)
             CALL HCO_ERROR( HcoState%Config%Err, MSG, RC )
             RETURN
          ENDIF

          ! Check hierarchy
          IF ( Lct%Dct%Hier /= TargetLct%Dct%Hier ) THEN
             MSG = 'Wrong hierarchy: ' // TRIM(Lct%Dct%cName)
             CALL HCO_ERROR( HcoState%Config%Err, MSG, RC )
             RETURN
          ENDIF

          ! Error check: cannot add masks if operator is 3
          IF ( Lct%Dct%DctType == HCO_DCTTYPE_MASK .AND. &
               Lct%Dct%Oper    == 3                       ) THEN
             MSG = 'Cannot add masks if operator is 3: ' // &
                  TRIM(Lct%Dct%cName)
             CALL HCO_ERROR( HcoState%Config%Err, MSG, RC )
             RETURN
          ENDIF

          ! If all checks were successful, add current array to
          ! target array.
          DO I = 1, TargetLct%Dct%Dta%nt
             IF ( TargetLct%Dct%Dta%SpaceDim <= 2 ) THEN
                TargetLct%Dct%Dta%V2(I)%Val = &
                     TargetLct%Dct%Dta%V2(I)%Val + Lct%Dct%Dta%V2(I)%Val
             ELSE
                TargetLct%Dct%Dta%V3(I)%Val = &
                     TargetLct%Dct%Dta%V3(I)%Val + Lct%Dct%Dta%V3(I)%Val
             ENDIF
          ENDDO
          IF ( verb ) THEN
             WRITE(MSG,*) 'Added data of ',   TRIM(Lct%Dct%cName), &
                  ' to ', TRIM(TargetLct%Dct%cName)
             CALL HCO_MSG(HcoState%Config%Err,MSG)
          ENDIF

          ! This container does not need to be added to the emissions
          ! list
          Add = .FALSE.

       ENDIF
    ENDIF ! cID /= targetID

    ! ----------------------------------------------------------------
    ! Add/update emissions linked list container.
    ! Only add those containers that are effectively used in the
    ! emissions list, i.e. ignore the containers whose content
    ! has been added to another container (targetID /= cID). Those
    ! containers are not needed for emission calculation since its
    ! content is now stored in another container.
    ! The EmisList_Add call will set the IsInList flag of the given
    ! file data object (Lct%Dct%Dta) to TRUE, denoting that this file
    ! data object is used in EmisList. The data arrays of all file
    ! data objects that are not used in EmisList are removed in a
    ! second step of the ReadList_Read call.
    ! ----------------------------------------------------------------
    IF ( Add ) THEN
       CALL EmisList_Add( Lct%Dct, HcoState, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN
    ENDIF

    ! ----------------------------------------------------------------
    ! Return
    ! ----------------------------------------------------------------
    CALL HCO_LEAVE( HcoState%Config%Err, RC )

  END SUBROUTINE EmisList_Pass
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCO_GetPtr_3D
!
! !DESCRIPTION: Subroutine HCO\_GetPtr\_3D returns the 3D data pointer
! Ptr3D of EmisList that is associated with data container DctName. By
! default, the routine returns an error if the given container name is
! not found. This can be avoided by calling the routine with the optional
! argument FOUND, in which case only this argument will be set to FALSE.
! Similarly, the FILLED flag can be used to control the behaviour if the
! data container is found but empty, e.g. no data is associated with it.
!\\
!\\
! This routine returns the unevaluated data field, e.g. no scale factors
! or masking is applied to the data. Use routine HCO\_EvalFld in
! hco\_calc\_mod.F90 to get evaluated fields.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCO_GetPtr_3D( HcoState, DctName, Ptr3D, &
                            RC, TIDX, FOUND, FILLED )
!
! !USES:
!
    USE HCO_DATACONT_MOD, ONLY : ListCont_Find
!
! !INPUT PARAMETERS:
!
    TYPE(HCO_State),  POINTER               :: HcoState       ! HEMCO state obj
    CHARACTER(LEN=*), INTENT(IN   )         :: DctName        ! container name
    INTEGER,          INTENT(IN), OPTIONAL  :: TIDX           ! time index
!                                                             ! (default=1)
! !OUTPUT PARAMETERS:
!
    REAL(sp),         POINTER               :: Ptr3D(:,:,:)   ! output array
    LOGICAL,          INTENT(OUT), OPTIONAL :: FOUND          ! cont. found?
    LOGICAL,          INTENT(OUT), OPTIONAL :: FILLED         ! array filled?
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(INOUT)         :: RC             ! Success/fail
!
! !REVISION HISTORY:
!  04 Sep 2013 - C. Keller - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER                    :: T
    LOGICAL                    :: FND
    CHARACTER(LEN=255)         :: MSG, LOC

    ! Pointers
    TYPE(ListCont), POINTER    :: Lct

    !=================================================================
    ! HCO_GetPtr_3D BEGINS HERE
    !=================================================================

    ! Enter
    LOC = 'HCO_GetPtr_3D (hco_emislist_mod.F90)'

    ! Init
    Lct => NULL()

    ! Define time index to use
    IF ( PRESENT(TIDX) ) THEN
       T = TIDX
    ELSE
       T = 1
    ENDIF

    ! Init
    IF ( PRESENT(FILLED) ) FILLED = .FALSE.

    ! Search for container in emissions linked list
    CALL ListCont_Find ( HcoState%EmisList, TRIM(DctName), FND, Lct )
    IF ( PRESENT(FOUND) ) FOUND = FND

    ! Check if found. If optional argument FOUND is defined, don't
    ! return an error if container not found but only pass the FOUND
    ! argument to the caller routine. Otherwise, exit with error.
    IF ( .NOT. FND ) THEN
       IF ( PRESENT(FOUND) .OR. PRESENT(FILLED) ) THEN
          Ptr3D => NULL()
          RC    = HCO_SUCCESS
          RETURN
       ELSE
          MSG = 'Container not found: ' // TRIM(DctName)
          CALL HCO_ERROR( HcoState%Config%Err, MSG, RC, THISLOC=LOC )
          RETURN
       ENDIF
    ENDIF

    ! Check spatial dimension
    IF ( Lct%Dct%Dta%SpaceDim /= 3 ) THEN
       MSG = 'Container is not 3D: ' // TRIM(DctName)
       CALL HCO_ERROR( HcoState%Config%Err, MSG, RC, THISLOC=LOC )
       RETURN
    ENDIF

    ! Check time dimension
    IF ( Lct%Dct%Dta%nt < T ) THEN
       MSG = 'not enough time slices: ' // TRIM(DctName)
       CALL HCO_ERROR( HcoState%Config%Err, MSG, RC, THISLOC=LOC )
       RETURN
    ENDIF

    IF ( ASSOCIATED( Lct%Dct%Dta%V3 ) ) THEN
       Ptr3D => Lct%Dct%Dta%V3(T)%Val
       IF ( PRESENT( FILLED ) ) FILLED = .TRUE.
    ELSE
       IF ( PRESENT( FILLED ) ) THEN
          Ptr3D  => NULL()
       ELSE
          MSG = 'Container data not filled: ' // TRIM(DctName)
          CALL HCO_ERROR( HcoState%Config%Err, MSG, RC, THISLOC=LOC )
          RETURN
       ENDIF
    ENDIF

    ! Leave w/ success
    RC = HCO_SUCCESS

  END SUBROUTINE HCO_GetPtr_3D
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCO_GetPtr_2D
!
! !DESCRIPTION: Subroutine HCO\_GetPtr\_2D returns the 2D data pointer
! Ptr2D of EmisList that is associated with data container DctName. See
! HCO\_GetPtr\_3D for more details.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCO_GetPtr_2D( HcoState, DctName, Ptr2D, &
                            RC, TIDX, FOUND, FILLED )
!
! !USES:
!
    USE HCO_DATACONT_MOD, ONLY : ListCont_Find
!
! !INPUT PARAMETERS:
!
    TYPE(HCO_State),  POINTER               :: HcoState    ! HEMCO state obj
    CHARACTER(LEN=*), INTENT(IN   )         :: DctName     ! container name
    INTEGER,          INTENT(IN), OPTIONAL  :: TIDX        ! time index
!                                                          ! (default=1)
! !OUTPUT PARAMETERS:
!
    REAL(sp),         POINTER               :: Ptr2D(:,:)  ! output array
    LOGICAL,          INTENT(OUT), OPTIONAL :: FOUND       ! cont. found?
    LOGICAL,          INTENT(OUT), OPTIONAL :: FILLED      ! array filled?
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(INOUT)         :: RC          ! Success/fail
!
! !REVISION HISTORY:
!  04 Sep 2013 - C. Keller - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER                    :: T
    LOGICAL                    :: FND
    CHARACTER(LEN=255)         :: MSG, LOC

    ! Pointers
    TYPE(ListCont), POINTER    :: Lct

    !=================================================================
    ! HCO_GetPtr_2D BEGINS HERE
    !=================================================================

    ! Enter
    LOC = 'HCO_GetPtr_2D (hco_emislist_mod.F90)'
    Lct => NULL()

    ! Define time index to use
    IF ( PRESENT(TIDX) )THEN
       T = TIDX
    ELSE
       T = 1
    ENDIF

    ! Init
    IF ( PRESENT(FILLED) ) FILLED = .FALSE.

    ! Search for container in emissions linked list
    CALL ListCont_Find( HcoState%EmisList, TRIM(DctName), FND, Lct )
    IF ( PRESENT(FOUND) ) FOUND = FND

    ! Check if found. If optional argument FOUND is defined, don't
    ! return an error if container not found but only pass the FOUND
    ! argument to the caller routine. Otherwise, exit with error.
    IF ( .NOT. FND ) THEN
       IF ( PRESENT(FOUND) .OR. PRESENT(FILLED) ) THEN
          Ptr2D => NULL()
          RC    = HCO_SUCCESS
          RETURN
       ELSE
          MSG = 'Container not found: ' // TRIM(DctName)
          CALL HCO_ERROR( HcoState%Config%Err, MSG, RC, THISLOC=LOC )
          RETURN
       ENDIF
    ENDIF

    ! Check spatial dimension
    IF ( (Lct%Dct%Dta%SpaceDim/=2) .AND. &
         (Lct%Dct%Dta%SpaceDim/=1)        ) THEN
       MSG = 'Container is not 2D: ' // TRIM(DctName)
       CALL HCO_ERROR( HcoState%Config%Err, MSG, RC, THISLOC=LOC )
       RETURN
    ENDIF

    ! Check time dimension
    IF ( Lct%Dct%Dta%nt < T ) THEN
       MSG = 'not enough time slices: ' // TRIM(DctName)
       CALL HCO_ERROR( HcoState%Config%Err, MSG, RC, THISLOC=LOC )
       RETURN
    ENDIF

    IF ( ASSOCIATED( Lct%Dct%Dta%V2 ) ) THEN
       Ptr2D => Lct%Dct%Dta%V2(T)%Val
       IF ( PRESENT( FILLED ) ) FILLED = .TRUE.
    ELSE
       IF ( PRESENT( FILLED ) ) THEN
          Ptr2D  => NULL()
       ELSE
          MSG = 'Container data not filled: ' // TRIM(DctName)
          CALL HCO_ERROR( HcoState%Config%Err, MSG, RC, THISLOC=LOC )
          RETURN
       ENDIF
    ENDIF

    ! Return w/ success
    RC = HCO_SUCCESS

  END SUBROUTINE HCO_GetPtr_2D
!EOC
END MODULE HCO_EMISLIST_MOD
