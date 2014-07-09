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
  USE HCO_DATACONT_MOD, ONLY : DataCont, ListCont
  USE HCO_STATE_MOD,    ONLY : HCO_State

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC  :: EmisList_Cleanup
  PUBLIC  :: EmisList_Update
  PUBLIC  :: EmisList_GetDataArr
  PUBLIC  :: EmisList_NextCont
!
! !PRIVATE MEMBER FUNCTIONS:
!
  PRIVATE :: EmisList_Pass
  PRIVATE :: EmisList_Add
  PRIVATE :: Add2EmisList
!
! !REVISION HISTORY:
!  04 Dec 2012 - C. Keller   - Initialization
!  08 Jul 2014 - R. Yantosca - Now use F90 free-format indentation
!  08 Jul 2014 - R. Yantosca - Cosmetic changes in ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !PRIVATE TYPES:
!
  ! EmisList is the emissions linked list
  TYPE(ListCont), POINTER    :: EmisList   => NULL()

  ! # of containers in EmisList
  INTEGER                    :: nnEmisCont =  0       

  INTERFACE EmisList_GetDataArr 
     MODULE PROCEDURE EmisList_GetDataArr_2D
     MODULE PROCEDURE EmisList_GetDataArr_3D
  END INTERFACE EmisList_GetDataArr

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
  SUBROUTINE EmisList_Add( am_I_Root, Dct, HcoState, RC )
!
! !USES:
!
    USE HCO_TIDX_MOD,      ONLY : tIDx_Assign 
    USE HCO_DATACONT_MOD,  ONLY : ListCont_Find, DataCont_Print
!
! !INPUT PARAMETERS:
!
    LOGICAL,           INTENT(IN)    :: am_I_Root  ! Root CPU?
    TYPE(DataCont),    POINTER       :: Dct        ! Data cont.
    TYPE(HCO_State),   POINTER       :: HcoState   ! HEMCO state
!
! !INPUT/OUTPUT PARAMETERS: 
!
    INTEGER,           INTENT(INOUT) :: RC         ! Return code
!
! !REVISION HISTORY:
!  04 Dec 2012 - C. Keller - Initialization
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Pointers
    TYPE(ListCont), POINTER :: Lct => NULL()

    ! Scalars
    LOGICAL                 :: FOUND, VERBOSE, NEW 
    CHARACTER(LEN=255)      :: MSG
    CHARACTER(LEN= 31)      :: TempRes

    !======================================================================
    ! EmisList_Add begins here!
    !======================================================================

    ! Enter
    CALL HCO_ENTER ('EmisList_Add (hco_emisslist_mod.F90)', RC )
    IF(RC /= HCO_SUCCESS) RETURN

    ! Set verbose flag
    VERBOSE = HCO_VERBOSE_CHECK() .AND. am_I_Root

    ! ----------------------------------------------------------------
    ! Nothing to do if it's not a new container, i.e. if container 
    ! already exists in EmisList. 
    ! ----------------------------------------------------------------
    CALL ListCont_Find( EmisList, Dct%cID, 0, FOUND, Lct )
    IF ( FOUND ) THEN
       CALL HCO_LEAVE( RC )
       RETURN
    ENDIF

    ! ----------------------------------------------------------------
    ! Define new list container. This container points to the passed 
    ! data container. 
    ! ----------------------------------------------------------------
    ALLOCATE ( Lct ) 
    Lct%NextCont => NULL()
    Lct%Dct      => Dct

    ! ----------------------------------------------------------------
    ! Set time index pointer tIDx of this data container. tIDx will
    ! be set according to the number of time slices (and the time
    ! interval between them) hold by this data container. For hourly
    ! data (24 time slices), for example, tIDx will point to the 
    ! corresponding 'HOURLY' or 'HOURLY_GRID' time index collection 
    ! type defined in hco\_tidx\_mod. 
    ! ----------------------------------------------------------------
    CALL tIDx_Assign( HcoState, Lct, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! ----------------------------------------------------------------
    ! Add the new container to EmisList. The container will be placed
    ! according to data type, species ID, hierarchy, and category. 
    ! ----------------------------------------------------------------
    CALL Add2EmisList( Lct, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! ----------------------------------------------------------------
    ! Verbose mode 
    ! ----------------------------------------------------------------
    IF ( VERBOSE ) THEN
       MSG = 'Container added to EmisList:'
       CALL HCO_MSG(MSG,SEP1='-')
       CALL DataCont_Print( Lct%Dct, VERBOSE )
    ENDIF

    ! Leave w/ success
    CALL HCO_LEAVE( RC )

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
  SUBROUTINE Add2EmisList( Lct, RC )
!
! !INPUT PARAMETERS:
!
    TYPE(ListCont), POINTER       :: Lct   
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(INOUT) :: RC
!
! !REVISION HISTORY:
!  06 Dec 2012 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER                   :: NEWCAT, NEWHIR, NEWSPC

    ! Pointers
    TYPE(ListCont), POINTER   :: TmpLct => NULL()

    !======================================================================
    ! Add2EmisList begins here!
    !======================================================================

    ! Enter 
    CALL HCO_ENTER ( 'Add2EmisList', RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Update number of containers in EmisList 
    nnEmisCont = nnEmisCont + 1

    ! If this is the first container, we can simply place it at the
    ! beginning of the list.
    IF ( nnEmisCont == 1 ) THEN
       EmisList => Lct
       CALL HCO_LEAVE(RC)
       RETURN ! Leave routine
    ENDIF

    ! Special case where the linked list consists of scale factors
    ! only: In this case, we can place the new container at the 
    ! beginning no matter of its content! 
    IF ( EmisList%Dct%DctType > 1 ) THEN
       Lct%NextCont => EmisList 
       EmisList        => Lct
       CALL HCO_LEAVE( RC )
       RETURN 
    ENDIF

    ! Get field species ID, category and priority of the new container 
    NEWSPC  = Lct%Dct%HcoID 
    NEWCAT  = Lct%Dct%Cat
    NEWHIR  = Lct%Dct%Hier

    ! Containers are listed with increasing species ID. If the current
    ! container has lower speciesID than the first container, just add
    ! it at the beginning of the list.     
    IF ( (NEWSPC > 0) .AND. (NEWSPC < EmisList%Dct%HcoID) ) THEN
       Lct%NextCont => EmisList
       EmisList     => Lct
       CALL HCO_LEAVE(RC)
       RETURN 
    ENDIF

    ! In case that the current container has the same species ID 
    ! as the first container in the list: If this container has
    ! lower category, or same category and lower hierarchy, place 
    ! it before the first container in the list:
    IF ( NEWSPC == EmisList%Dct%HcoID ) THEN
       IF ( (EmisList%Dct%Cat  >  NEWCAT) .OR.  &
            (EmisList%Dct%Cat  == NEWCAT  .AND. &
            EmisList%Dct%Hier >  NEWHIR)        ) THEN
          Lct%NextCont => EmisList
          EmisList     => Lct
          CALL HCO_LEAVE(RC)
          RETURN 
       ENDIF
    ENDIF

    ! TmpLct is the temporary working pointer, looping through
    ! the entire EmisList until the correct place for the new
    ! container is found.
    TmpLct => EmisList

    ! If the new container contains base data (i.e. data type is 1) 
    ! we have to move the TmpLct pointer to the position where the 
    ! next container is also a base container and one of the 
    ! following: (a) the first container with the same species ID 
    ! as the new container; (b) a container with higher species ID; 
    ! (c) scale factors. From there, we can determine where to place 
    ! the container exactly.
    IF ( Lct%Dct%DctType == 1 ) THEN

       ! Loop over list
       DO WHILE ( ASSOCIATED ( TmpLct%NextCont ) )
             
          ! Check if next container's species ID is higher or if it's a
          ! scale factor, in which case we have to exit.
          IF ( TmpLct%NextCont%Dct%HcoID    > NEWSPC .OR. & 
               TmpLct%NextCont%Dct%DctType > 1             ) THEN
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
          IF ( TmpLct%NextCont%Dct%DctType > 1 ) EXIT

          ! Advance in list
          TmpLct => TmpLct%NextCont
       ENDDO

    ENDIF

    ! Add new container AFTER current one
    Lct%NextCont    => TmpLct%NextCont
    TmpLct%NextCont => Lct 

    ! Cleanup and leave
    TmpLct => NULL()
    CALL HCO_LEAVE( RC )

  END SUBROUTINE Add2EmisList
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: EmisList_Cleanup 
!
! !DESCRIPTION: Subroutine EmisList\_Cleanup cleans up EmisList. If 
! RemoveDct is set to true, all data containers associated to the 
! EmisList elements will be removed as well. Otherwise, only the data 
! container pointers will be removed, i.e. nullified.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE EmisList_Cleanup( RemoveDct )
!
! !USES:
!
    USE HCO_DATACONT_MOD,  ONLY : ListCont_Cleanup
!
! !INPUT PARAMETERS:
!
    LOGICAL, INTENT(IN)  :: RemoveDct
!
! !REVISION HISTORY:
!  04 Dec 2012 - C. Keller: Initialization
!EOP
!------------------------------------------------------------------------------
!BOC
     
    !======================================================================
    ! EmisList_Cleanup begins here!
    !======================================================================

    ! Cleanup EmisList
    CALL ListCont_Cleanup( EmisList, RemoveDct )

  END SUBROUTINE EmisList_Cleanup
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: EmisList_Update
!
! !DESCRIPTION: Subroutine EmisList\_Update makes sure that all containers
! of the reading list ReadList are correctly referenced in emissions list
! EmisList. If a container of ReadList does not yet have a corresponding 
! container in EmisList, such a container is created. Also, additive data 
! arrays (i.e. targetID different than container ID) are added to their 
! target array during this call.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE EmisList_Update ( am_I_Root, HcoState, ReadList, RC )
!
! !USES:
!
    USE HCO_FILEDATA_MOD, ONLY : FileData_ArrIsDefined
!
! !INPUT PARAMETERS:
!
    LOGICAL,         INTENT(IN   ) :: am_I_Root  ! Root CPU?  
    TYPE(HCO_State), POINTER       :: HcoState   ! Hemco state object
    TYPE(ListCont),  POINTER       :: ReadList   ! reading list
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,         INTENT(INOUT) :: RC         ! Return code
!
! !REVISION HISTORY:
!  20 Apr 2013 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER                     :: cID, iScalID
    INTEGER                     :: I
    CHARACTER(LEN=31)           :: ScalName

    ! Pointers
    TYPE(ListCont), POINTER     :: TmpLct => NULL()

    ! ================================================================
    ! EmisList_Update begins here
    ! ================================================================

    ! Enter
    CALL HCO_ENTER ( 'EmisList_Update', RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Loop over all containers in ReadList
    TmpLct => ReadList
    DO WHILE ( ASSOCIATED( TmpLct ) )

       ! only if array is defined...
       IF ( FileData_ArrIsDefined(TmpLct%Dct%Dta) ) THEN

          ! Pass container to EmisList 
          CALL EmisList_Pass( am_I_Root, HcoState, TmpLct, RC )
          IF ( RC /= HCO_SUCCESS ) RETURN
       ENDIF

       ! Advance to next container in ReadList 
       TmpLct => TmpLct%NextCont
    ENDDO

    ! Leave w/ success
    CALL HCO_LEAVE ( RC )

  END SUBROUTINE EmisList_Update
!EOC
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
  SUBROUTINE EmisList_Pass( am_I_Root, HcoState, Lct, RC )
!
! !USES:
!
    USE HCO_DATACONT_MOD, ONLY : ListCont_Find
    USE HCO_FILEDATA_MOD, ONLY : FileData_ArrCheck2D
    USE HCO_FILEDATA_MOD, ONLY : FileData_ArrCheck3D
    USE HCO_DATACONT_MOD, ONLY : DataCont_Cleanup
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN)    :: am_I_Root
    TYPE(HCO_State),  POINTER       :: HcoState
    TYPE(ListCont),   POINTER       :: Lct        ! list container 
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(INOUT) :: RC         ! Success 
!
! !REVISION HISTORY:
!  28 Mar 2013 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Pointers
    TYPE(ListCont),  POINTER  :: TargetLct => NULL()

    ! Scalars
    INTEGER                   :: I, J, L, T
    LOGICAL                   :: FOUND, verb
    CHARACTER(LEN=255)        :: MSG
 
    ! ================================================================
    ! EmisList_Pass begins here
    ! ================================================================

    ! Enter
    CALL HCO_ENTER ('EmisList_Pass (hco_emislist_mod.F90)', RC )
    IF(RC /= HCO_SUCCESS) RETURN

    ! Verbose mode
    verb = HCO_VERBOSE_CHECK() .AND. am_I_Root

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
    IF ( Lct%Dct%targetID /= Lct%Dct%cID ) THEN

       ! TargetLct points to the container holding the target array 
       CALL ListCont_Find( EmisList, Lct%Dct%targetID, &
                           0, FOUND, TargetLct )
       IF ( .NOT. FOUND ) THEN
          MSG = 'Cannot add emissions to target array: error in ' // &
               TRIM(Lct%Dct%cName)
          CALL HCO_ERROR( MSG, RC )
          RETURN
       ENDIF

       ! Check extension number 
       IF ( Lct%Dct%ExtNr /= TargetLct%Dct%ExtNr ) THEN
          MSG = 'Wrong ext. number: ' // TRIM(Lct%Dct%cName)
          CALL HCO_ERROR( MSG, RC )
          RETURN
       ENDIF

       ! Check data type 
       IF ( Lct%Dct%DctType /= TargetLct%Dct%DctType ) THEN
          MSG = 'Wrong data type: ' // TRIM(Lct%Dct%cName)
          CALL HCO_ERROR( MSG, RC )
          RETURN
       ENDIF

       ! Check species ID
       IF ( Lct%Dct%HcoID /= TargetLct%Dct%HcoID ) THEN
          MSG = 'Wrong species ID: ' // TRIM(Lct%Dct%cName)
          CALL HCO_ERROR( MSG, RC )
          RETURN
       ENDIF

       ! Check for array dimensions
       IF ( Lct%Dct%Dta%SpaceDim /= TargetLct%Dct%Dta%SpaceDim ) THEN
          MSG = 'Wrong space dimension: ' // TRIM(Lct%Dct%cName)
          CALL HCO_ERROR( MSG, RC )
          RETURN
       ENDIF
       IF ( Lct%Dct%Dta%nt /= TargetLct%Dct%Dta%nt ) THEN 
          MSG = 'Wrong time dim: ' // TRIM(Lct%Dct%cName)
          CALL HCO_ERROR( MSG, RC )
          RETURN
       ENDIF
       IF ( Lct%Dct%Dta%SpaceDim <= 2) THEN
          I = SIZE(Lct%Dct%Dta%V2(1)%Val,1)
          J = SIZE(Lct%Dct%Dta%V2(1)%Val,2)
          CALL FileData_ArrCheck2D( TargetLct%Dct%Dta, I, J, &
                                    Lct%Dct%Dta%nt, RC )
          IF ( RC /= 0 ) THEN
             MSG = 'Wrong 2D array: ' // TRIM(Lct%Dct%cName)
             CALL HCO_MSG(MSG)
             RETURN
          ENDIF
       ELSE
          I = SIZE(Lct%Dct%Dta%V3(1)%Val,1)
          J = SIZE(Lct%Dct%Dta%V3(1)%Val,2)
          L = SIZE(Lct%Dct%Dta%V3(1)%Val,3)
          CALL FileData_ArrCheck3D( TargetLct%Dct%Dta, I, J, L, &
                                    Lct%Dct%Dta%nt, RC )
          IF ( RC /= 0 ) THEN
             MSG = 'Wrong 3D array: ' // TRIM(Lct%Dct%cName)
             CALL HCO_MSG(MSG)
             RETURN
          ENDIF
       ENDIF

       ! Check operator 
       IF ( Lct%Dct%Oper /= TargetLct%Dct%Oper ) THEN
          MSG = 'Wrong operator: ' // TRIM(Lct%Dct%cName)
          CALL HCO_ERROR( MSG, RC )
          RETURN
       ENDIF

       ! Check category 
       IF ( Lct%Dct%Cat /= TargetLct%Dct%Cat ) THEN
          MSG = 'Wrong category: ' // TRIM(Lct%Dct%cName)
          CALL HCO_ERROR ( MSG, RC )
          RETURN
       ENDIF

       ! Check hierarchy 
       IF ( Lct%Dct%Hier /= TargetLct%Dct%Hier ) THEN
          MSG = 'Wrong hierarchy: ' // TRIM(Lct%Dct%cName)
          CALL HCO_ERROR( MSG, RC )
          RETURN
       ENDIF

       ! Error check: cannot add masks if operator is 3
       IF ( Lct%Dct%DctType == 3 .AND. Lct%Dct%Oper == 3 ) THEN
          MSG = 'Cannot add masks if operator is 3: ' // &
               TRIM(Lct%Dct%cName)
          CALL HCO_ERROR( MSG, RC )
          RETURN
       ENDIF

       ! If all checks were successful, add current array to
       ! target array. Don't do this if this is not the home 
       ! data container for the linked file data object Dta!
       ! Dta may be used by multiple containers, and only the
       ! home container should modify it's content!
       IF ( Lct%Dct%DtaHome /= 1 ) THEN
         
          ! Verbose mode
          IF ( verb ) THEN
             WRITE(MSG,*) 'Do not add data of ', TRIM(Lct%Dct%cName), &
                  'to ', TRIM(TargetLct%Dct%cName), ' because this is', &
                  ' not the file data home container!'
             CALL HCO_MSG(MSG)
          ENDIF

       ! Add to array
       ELSE 
    
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
             CALL HCO_MSG(MSG)
          ENDIF
   
          ! Added array is not used anymore and thus can be cleaned up.
          ! Note: never clean up shared arrays!
          IF ( .NOT. Lct%Dct%Dta%DoShare ) THEN
             CALL DataCont_Cleanup( Lct%Dct, ArrOnly=.TRUE. ) 
          ENDIF

       ENDIF

    ! ----------------------------------------------------------------
    ! Add/update emissions linked list container.
    ! Note: only add target containers, i.e. don't add containers 
    ! whose content has been written to another container. Those 
    ! containers are not needed for emission calculation since its
    ! content is now stored in another container.
    ! ----------------------------------------------------------------

    ELSE
       CALL EmisList_Add( am_I_Root, Lct%Dct, HcoState, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN

    ENDIF ! cID /= targetID

    ! ----------------------------------------------------------------
    ! Return 
    ! ----------------------------------------------------------------
    CALL HCO_LEAVE( RC )

  END SUBROUTINE EmisList_Pass
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: EmisList_GetDataArr_3D 
!
! !DESCRIPTION: Subroutine EMISLIST\_GetDataArr\_3D returns the 3D data 
! array Arr3D of EmisList that is associated with data container DctName. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE EmisList_GetDataArr_3D( am_I_Root, DctName, Arr3D, RC, TIDX )
!
! !USES:
!
    USE HCO_DATACONT_MOD, ONLY : ListCont_Find
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN   )        :: am_I_Root      ! root CPU?
    CHARACTER(LEN=*), INTENT(IN   )        :: DctName        ! container name
    INTEGER,          INTENT(IN), OPTIONAL :: TIDX           ! time index
!                                                            ! (default=1)
! !OUTPUT PARAMETERS:
!
    REAL(hp),         POINTER              :: Arr3D(:,:,:)   ! output array
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(INOUT)        :: RC             ! Success/fail
!
! !REVISION HISTORY: 
!  04 Sep 2013 - C. Keller    - Initialization
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER                    :: T
    LOGICAL                    :: FOUND
    CHARACTER(LEN=255)         :: MSG, LOC

    ! Pointers
    TYPE(ListCont), POINTER    :: Lct => NULL()

    !=================================================================
    ! EmisList_GetDataArr_3D BEGINS HERE
    !=================================================================

    ! Enter
    LOC = 'EmisList_GetDataArr_3D (hco_emislist_mod.F90)'

    ! Define time index to use
    IF ( PRESENT(TIDX) ) THEN
       T = TIDX
    ELSE
       T = 1
    ENDIF
 
    ! Search for container in emissions linked list
    CALL ListCont_Find ( EmisList, TRIM(DctName), FOUND, Lct )

    ! Error check
    IF ( .NOT. FOUND ) THEN
       MSG = 'Container not found: ' // TRIM(DctName)
       CALL HCO_ERROR( MSG, RC, THISLOC=LOC )
       RETURN
    ENDIF

    ! Check spatial dimension 
    IF ( Lct%Dct%Dta%SpaceDim /= 3 ) THEN
       MSG = 'Container is not 3D: ' // TRIM(DctName)
       CALL HCO_ERROR( MSG, RC, THISLOC=LOC )
       RETURN
    ENDIF

    ! Check time dimension 
    IF ( Lct%Dct%Dta%nt < T ) THEN
       MSG = 'not enough time slices: ' // TRIM(DctName)
       CALL HCO_ERROR( MSG, RC, THISLOC=LOC )
       RETURN
    ENDIF

    IF ( ASSOCIATED( Lct%Dct%Dta%V3 ) ) THEN
       Arr3D => Lct%Dct%Dta%V3(T)%Val
    ELSE
       MSG = 'Container data not filled: ' // TRIM(DctName)
       CALL HCO_ERROR( MSG, RC, THISLOC=LOC )
       RETURN
    ENDIF

    ! Leave w/ success
    RC = HCO_SUCCESS

  END SUBROUTINE EmisList_GetDataArr_3D
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: EmisList_GetDataArr_2D 
!
! !DESCRIPTION: Subroutine EmisList\_GetDataArr\_2D returns the 2D data 
! array Arr2D of EmisList that is associated with data container DctName. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE EmisList_GetDataArr_2D( am_I_Root, DctName, Arr2D, RC, TIDX )
!
! !USES:
!
    USE HCO_DATACONT_MOD, ONLY : ListCont_Find
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN   )        :: am_I_Root   ! root CPU?
    CHARACTER(LEN=*), INTENT(IN   )        :: DctName     ! container name
    INTEGER,          INTENT(IN), OPTIONAL :: TIDX        ! time index
!                                                         ! (default=1)
! !OUTPUT PARAMETERS:
!
    REAL(hp),         POINTER              :: Arr2D(:,:)  ! output array
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(INOUT)        :: RC          ! Success/fail
!
! !REVISION HISTORY: 
!  04 Sep 2013 - C. Keller    - Initialization
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER                    :: T
    LOGICAL                    :: FOUND
    CHARACTER(LEN=255)         :: MSG, LOC

    ! Pointers
    TYPE(ListCont), POINTER    :: Lct => NULL()

    !=================================================================
    ! EmisList_GetDataArr_2D BEGINS HERE
    !=================================================================

    ! Enter
    LOC = 'EmisList_GetDataArr_2D (hco_emislist_mod.F90)'

    ! Define time index to use
    IF ( PRESENT(TIDX) )THEN
       T = TIDX
    ELSE
       T = 1
    ENDIF
 
    ! Search for container in emissions linked list
    CALL ListCont_Find( EmisList, TRIM(DctName), FOUND, Lct )

    ! Error check
    IF ( .NOT. FOUND ) THEN
       MSG = 'Container not found: ' // TRIM(DctName)
       CALL HCO_ERROR( MSG, RC, THISLOC=LOC )
       RETURN
    ENDIF

    ! Check spatial dimension 
    IF ( (Lct%Dct%Dta%SpaceDim/=2) .AND. &
         (Lct%Dct%Dta%SpaceDim/=1)        ) THEN
       MSG = 'Container is not 2D: ' // TRIM(DctName)
       CALL HCO_ERROR( MSG, RC, THISLOC=LOC )
       RETURN
    ENDIF

    ! Check time dimension 
    IF ( Lct%Dct%Dta%nt < T ) THEN
       MSG = 'not enough time slices: ' // TRIM(DctName)
       CALL HCO_ERROR( MSG, RC, THISLOC=LOC )
       RETURN
    ENDIF

    IF ( ASSOCIATED( Lct%Dct%Dta%V2 ) ) THEN
       Arr2D => Lct%Dct%Dta%V2(T)%Val
    ELSE
       MSG = 'Container data not filled: ' // TRIM(DctName)
       CALL HCO_ERROR( MSG, RC, THISLOC=LOC )
       RETURN
    ENDIF

    ! Leave w/ success
    RC = HCO_SUCCESS

  END SUBROUTINE EmisList_GetDataArr_2D
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: EmisList_NextCont
!
! !DESCRIPTION: Subroutine EmisList\_NextCont is a simple wrapper routine 
! to get the next/first pointer of the ConfigList. See ListCont\_NextCont 
! for more details. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE EmisList_NextCont( Lct, FLAG ) 
!
! !USES:
!
    USE HCO_DATACONT_MOD,  ONLY : ListCont_NextCont
!
! !OUTPUT PARAMETERS:
!
    TYPE(ListCont),   POINTER       :: Lct
    INTEGER,          INTENT(INOUT) :: FLAG
!
! !REVISION HISTORY: 
!  04 Sep 2013 - C. Keller    - Initialization
!EOP
!------------------------------------------------------------------------------
!BOC

    !=================================================================
    ! EmisList_NextCont BEGINS HERE
    !=================================================================

    CALL ListCont_NextCont( EmisList, Lct, FLAG )

  END SUBROUTINE EmisList_NextCont
!EOC
END MODULE HCO_EMISLIST_MOD
