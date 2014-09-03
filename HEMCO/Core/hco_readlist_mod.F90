!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: hco_readlist_mod.F90
!
! !DESCRIPTION: Module HCO\_ReadList\_Mod contains routines and variables
! for the HEMCO ReadList. ReadList is a collection of all data containers
! used by HEMCO. They are categorized according to their reading update
! frequency, i.e all data containers that need to be updated on an annual
! basis are stored in ReadList 'Year', etc. The following reading update
! frequencies are supported:
!
! \begin{itemize}
! \item Year: update every year (annual data) 
! \item Month: update every month (monthly data) 
! \item Day: update every day (daily data) 
! \item Hour: update every hour (hourly data) 
! \item Once: update only once (time-invariant data)
! \item Always: update every time step
! \end{itemize} 
!
! In an ESMF environment - where data reading and time interpolation is 
! done through the ESMF/MAPL software framework - all arrays are added
! to list 'Always' and hence become refreshed on every time step.
!\\
!\\
! !INTERFACE: 
!
MODULE HCO_ReadList_Mod
!
! !USES:
!
  USE HCO_Error_MOD
  USE HCO_DataCont_MOD, ONLY : DataCont, ListCont
  USE HCO_State_MOD,    ONLY : HCO_State
 
  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC  :: ReadList_Init
  PUBLIC  :: ReadList_Read
  PUBLIC  :: ReadList_Set
  PUBLIC  :: ReadList_Print
  PUBLIC  :: ReadList_Cleanup
  PUBLIC  :: ReadList_to_EmisList
!
! !PRIVATE MEMBER FUNCTIONS:
!
  PRIVATE :: DtCont_Add
  PRIVATE :: ReadList_Fill
!
! !REVISION HISTORY:
!  20 Apr 2013 - C. Keller   - Initial version
!  01 Jul 2014 - R. Yantosca - Cosmetic changes in ProTeX headers
!  01 Jul 2014 - R. Yantosca - Now use F90 free-format indentation
!EOP
!-----------------------------------------------------------------------------
!BOC
!
! !PRIVATE TYPES:
!
  ! RdList contains lists to all ReadLists
  TYPE RdList
     TYPE(ListCont), POINTER :: Once
     TYPE(ListCont), POINTER :: Always
     TYPE(ListCont), POINTER :: Year
     TYPE(ListCont), POINTER :: Month
     TYPE(ListCont), POINTER :: Day
     TYPE(ListCont), POINTER :: Hour
  END TYPE RdList

  ! Internally used ReadLists
  TYPE(RdList),      POINTER :: ReadLists => NULL()

CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ReadList_Set
!
! !DESCRIPTION: Subroutine ReadList\_Set places the passed data container
! Dct in one of the reading lists, according to the data update 
! frequency specified in the HEMCO configuration file. Containers are 
! sorted with increasing container ID.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE ReadList_Set( am_I_Root, HcoState, Dct, RC )
!
! !USES:
!
    USE HCO_LOGFILE_MOD, ONLY : HCO_PrintDataCont
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN   )  :: am_I_Root
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(HCO_State),  POINTER        :: HcoState
    TYPE(DataCont),   POINTER        :: Dct
    INTEGER,          INTENT(INOUT)  :: RC
!
! !REVISION HISTORY:
!  20 Apr 2013 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER            :: intv
    LOGICAL            :: verb
    CHARACTER(LEN=255) :: MSG

    ! ================================================================
    ! ReadList_Set begins here
    ! ================================================================

    ! For error handling
    CALL HCO_ENTER ('ReadList_Set (hco_readlist_mod.F90)', RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Verbose mode
    verb = am_I_Root .AND. HCO_VERBOSE_CHECK()

    ! Add container to ReadList according to update freqency.
    ! Fields in list 'Hour' will be updated (i.e. re-read) every hour, 
    ! fields in list 'Day' every day, etc.
    ! If a time range instead of a single time stamp is given,
    ! categorize the field according to the most rapidly changing time
    ! stamp. If no time attribute exist, put the field into the 'Once'
    ! list which reads the file only at the beginning. 
    IF (     Dct%Dta%ncHrs(1) /= Dct%Dta%ncHrs(2) ) THEN
       intv = 1  
    ELSEIF ( Dct%Dta%ncDys(1) /= Dct%Dta%ncDys(2) ) THEN
       intv = 2
    ELSEIF ( Dct%Dta%ncMts(1) /= Dct%Dta%ncMts(2) ) THEN
       intv = 3
    ELSEIF ( Dct%Dta%ncYrs(1) /= Dct%Dta%ncYrs(2) ) THEN
       intv = 4
    ELSE
       intv = 5
    ENDIF

    ! NOTE: In an ESMF environment, data I/O is organized through 
    ! ESMF/MAPL. These routines interpolate between all timesteps,
    ! we hence need to update ReadList (and EmisList) every time!
    ! The only exception to this are the one-time lists, which don't
    ! need to be renewed at all!
    IF ( HcoState%isESMF ) THEN
       IF ( intv /= 5 ) THEN
          CALL DtCont_Add( ReadLists%Always, Dct )
       ELSE
          CALL DtCont_Add( ReadLists%Once,  Dct ) 
       ENDIF
    ELSE 
       IF ( intv == 1 ) THEN 
          CALL DtCont_Add( ReadLists%Hour,  Dct ) 
       ELSEIF ( intv == 2 ) THEN 
          CALL DtCont_Add( ReadLists%Day,   Dct ) 
       ELSEIF ( intv == 3 ) THEN 
          CALL DtCont_Add( ReadLists%Month, Dct ) 
       ELSEIF ( intv == 4 ) THEN 
          CALL DtCont_Add( ReadLists%Year,  Dct ) 
       ELSE
          CALL DtCont_Add( ReadLists%Once,  Dct ) 
       ENDIF
    ENDIF

    ! Verbose
    IF ( Verb ) THEN
       write(MSG,*) 'New container set to ReadList:'
       CALL HCO_MSG(MSG,SEP1='-')
       CALL HCO_PrintDataCont( Dct, Verb )
    ENDIF

    ! Leave w/ success
    CALL HCO_LEAVE ( RC ) 

  END SUBROUTINE ReadList_Set
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ReadList_Read 
!
! !DESCRIPTION: Subroutine ReadList\_Read makes sure that all arrays in the
! reading lists are up to date, i.e. it invokes the data reading calls for 
! those lists that need to be refreshed.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE ReadList_Read( am_I_Root, HcoState, RC ) 
!
! !USES:
!
    USE HCO_CLOCK_MOD, ONLY : HcoClock_First,    HcoClock_NewYear
    USE HCO_CLOCK_MOD, ONLY : HcoClock_NewMonth, HcoClock_NewDay
    USE HCO_CLOCK_MOD, ONLY : HcoClock_NewHour
!
! !INPUT PARAMETERS:
!
    LOGICAL,         INTENT(IN   )  :: am_I_Root  ! root CPU?
!
! !INPUT/OUTPUT PARAMETERS: 
!
    TYPE(HCO_State), POINTER        :: HcoState   ! HEMCO state object
    INTEGER,         INTENT(INOUT)  :: RC         ! Success or failure?
!
! !REVISION HISTORY:
!  20 Apr 2013 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    LOGICAL            :: verb
    CHARACTER(LEN=255) :: MSG

    ! ================================================================
    ! ReadList_Read begins here
    ! ================================================================

    ! For error handling
    CALL HCO_ENTER ('ReadList_Read (hco_readlist_mod.F90)', RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Verbose mode
    verb = am_I_Root .AND. HCO_VERBOSE_CHECK()

    ! Read content from one-time list on the first call 
    IF ( HcoClock_First() ) THEN
       IF ( Verb ) THEN
          WRITE(MSG,*) 'Now reading once list!'
          CALL HCO_MSG(MSG)
       ENDIF
       CALL ReadList_Fill ( am_I_Root, HcoState, ReadLists%Once, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN 
    ENDIF

    ! Read content from year-list if it's a new year
    IF ( HcoClock_NewYear() ) THEN
       IF ( Verb ) THEN
          WRITE(MSG,*) 'Now reading year list!'
          CALL HCO_MSG(MSG)
       ENDIF
       CALL ReadList_Fill ( am_I_Root, HcoState, ReadLists%Year, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN 
    ENDIF

    ! Read content from month-list if it's a new month
    IF ( HcoClock_NewMonth() ) THEN
       IF ( Verb ) THEN
          WRITE(MSG,*) 'Now reading month list!'
          CALL HCO_MSG(MSG)
       ENDIF
       CALL ReadList_Fill ( am_I_Root, HcoState, ReadLists%Month, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN 
    ENDIF

    ! Read content from day-list if it's a new day 
    IF ( HcoClock_NewDay() ) THEN
       IF ( Verb ) THEN
          WRITE(MSG,*) 'Now reading day list!'
          CALL HCO_MSG(MSG)
       ENDIF
       CALL ReadList_Fill ( am_I_Root, HcoState, ReadLists%Day, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN 
    ENDIF

    ! Read content from hour-list if it's a new hour 
    IF ( HcoClock_NewHour() ) THEN
       IF ( Verb ) THEN
          WRITE(MSG,*) 'Now reading hour list!'
          CALL HCO_MSG(MSG)
       ENDIF
       CALL ReadList_Fill ( am_I_Root, HcoState, ReadLists%Hour, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN 
    ENDIF

    ! Always add/update content from always-list
    IF ( Verb ) THEN
       WRITE(MSG,*) 'Now reading always list!'
       CALL HCO_MSG(MSG)
    ENDIF
    CALL ReadList_Fill ( am_I_Root, HcoState, ReadLists%Always, RC ) 
    IF ( RC /= HCO_SUCCESS ) RETURN 

    ! Leave w/ success
    CALL HCO_LEAVE ( RC ) 

  END SUBROUTINE ReadList_Read
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ReadList_Fill
!
! !DESCRIPTION: Subroutine ReadList\_Fill (re-)reads the data from all
! containers of the passed ReadList. In a non-ESMF environment, this
! routine calls the HEMCO generic (netCDF) reading and remapping 
! routines. In an ESMF environment, the arrays are obtained through 
! the ESMF/MAPL software framework.\\
! This routine calls the HEMCO - data reading interface wrapper routine
! for the current environment. See module HCOI\_DATAREAD\_MOD.F90 for
! more details. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE ReadList_Fill( am_I_Root, HcoState, ReadList, RC ) 
!
! !USES:
!
    USE HCOIO_DataRead_Mod, ONLY : HCOIO_DataRead
    USE HCOIO_DataRead_Mod, ONLY : HCOIO_ReadFromConfig
    USE HCO_FileData_Mod,   ONLY : FileData_ArrIsDefined
!
! !INPUT PARAMETERS:
!
    LOGICAL,         INTENT(IN   )  :: am_I_Root  ! root CPU?
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(HCO_State), POINTER        :: HcoState   ! HEMCO state object
    TYPE(ListCont),  POINTER        :: ReadList   ! Current reading list
    INTEGER,         INTENT(INOUT)  :: RC         ! Success or failure?
!
! !REMARKS: 
!  Different HCOI_DATAREAD routines may be invoked depending on the 
!  model environment. 
!
! !REVISION HISTORY:
!  20 Apr 2013 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    TYPE(ListCont), POINTER  :: Lct => NULL()

    ! ================================================================
    ! ReadList_Fill begins here
    ! ================================================================

    ! For error handling
    CALL HCO_ENTER ('ReadList_Fill (hco_readlist_mod.F90)', RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Loop over all containers
    Lct => ReadList
    DO WHILE ( ASSOCIATED ( Lct ) ) 

       ! Check if data has already been read. Multiple data
       ! containers can have the same file data object, and we
       ! only need to read it once. For each file data object,
       ! we assign a 'home container'. Reading will only be
       ! performed on this container. The DtaHome flag is
       ! initialized to -999. Hence, the first time we read data
       ! for a given container, check if this data file object
       ! has not yet been read, in which case we define this 
       ! container as the home container (flag=1). Otherwise, set
       ! flag to 0 (not home), but make sure that the DoShare flag
       ! of the corresponding data file object is enabled.
       IF ( Lct%Dct%DtaHome == -999 ) THEN
          IF ( FileData_ArrIsDefined(Lct%Dct%Dta) ) THEN
             Lct%Dct%DtaHome     = 0
             Lct%Dct%Dta%DoShare = .TRUE.
          ELSE
             Lct%Dct%DtaHome = 1
          ENDIF
       ENDIF

       ! Read if this is the home container
       IF ( Lct%Dct%DtaHome == 1 ) THEN

          ! Read from configuration file if it's not a netCDF file
          IF ( .NOT. Lct%Dct%Dta%NcRead ) THEN
             CALL HCOIO_ReadFromConfig ( am_I_Root, HcoState, Lct, RC )
             IF ( RC /= HCO_SUCCESS ) RETURN

          ! Read from netCDF file otherwise
          ELSE
             CALL HCOIO_DATAREAD ( am_I_Root, HcoState, Lct, RC )
             IF ( RC /= HCO_SUCCESS ) RETURN
          ENDIF

       ENDIF

       ! Point to next container
       Lct => Lct%NextCont
    ENDDO

    ! Cleanup
    Lct => NULL()

    ! Leave with success 
    CALL HCO_LEAVE ( RC )

  END SUBROUTINE ReadList_Fill
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: DtCont_Add
!
! !DESCRIPTION: Subroutine DtCont\_Add adds a new container to the 
! specified reading list.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE DtCont_Add( ReadList, Dct )
!
! !INPUT PARAMETERS:
!
    TYPE(ListCont), POINTER :: ReadList 
    TYPE(DataCont), POINTER :: Dct
!
! !REVISION HISTORY:
!  20 Apr 2013 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    TYPE(ListCont), POINTER :: Lct => NULL()
    TYPE(ListCont), POINTER :: TmpLct => NULL()
    INTEGER                 :: cID

    ! ================================================================
    ! DtCont_Add begins here
    ! ================================================================

    ! Create new container (to be added to ReadList)
    ALLOCATE ( Lct ) 
    Lct%NextCont => NULL()
      
    ! Make container point to passed data container
    Lct%Dct => Dct 

    ! If ReadList is not defined yet, set new container to head of 
    ! list.
    IF ( .NOT. ASSOCIATED ( ReadList ) ) THEN
       ReadList => Lct
     
    ! If list is already defined, place current container according
    ! to its container ID. Containers are sorted with increasing
    ! cID. 
    ELSE 

       ! cID of current container
       cID = Lct%Dct%cID
 
       ! TmpLct is the temporary pointer to the ReadList containers
       TmpLct => ReadList

       ! Check if cID of first container is higher than current ID. 
       ! In this case, we can place current container at beginning 
       ! of list
       IF ( TmpLct%Dct%cID > cID ) THEN
          Lct%NextCont => TmpLct
          ReadList     => Lct
       ELSE         

          ! Loop over containers in list until we encounter first
          ! container where the upcoming container has higher cID 
          ! than currCont. 
          DO WHILE ( ASSOCIATED ( TmpLct%NextCont ) ) 
             IF ( TmpLct%NextCont%Dct%cID > cID ) EXIT
             TmpLct => TmpLct%NextCont
          ENDDO

          ! Now place current container AFTER TmpLct
          Lct%NextCont    => TmpLct%NextCont
          TmpLct%NextCont => Lct
       ENDIF
    ENDIF

    ! Cleanup
    TmpLct => NULL()
    Lct    => NULL()

  END SUBROUTINE DtCont_Add
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ReadList_Init
!
! !DESCRIPTION: Subroutine ReadList\_Init initializes the ReadList. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE ReadList_Init()
!
! !REVISION HISTORY:
!  20 Apr 2013 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC

    ! ================================================================
    ! ReadList_Init begins here
    ! ================================================================
      
    ! Allocate ReadList and all internal lists. Make sure all internal
    ! lists are defined (nullified). 
    ALLOCATE ( ReadLists )

    ALLOCATE ( ReadLists%Once )
    NULLIFY ( ReadLists%Once  )

    ALLOCATE ( ReadLists%Year )
    NULLIFY ( ReadLists%Year  )

    ALLOCATE ( ReadLists%Month )
    NULLIFY ( ReadLists%Month )

    ALLOCATE ( ReadLists%Day )
    NULLIFY ( ReadLists%Day   )

    ALLOCATE ( ReadLists%Hour )
    NULLIFY ( ReadLists%Hour  )

    ALLOCATE ( ReadLists%Always )
    NULLIFY ( ReadLists%Always  )

  END SUBROUTINE ReadList_Init
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ReadList_Print
!
! !DESCRIPTION: Subroutine ReadList\_Print displays the content of
! ReadList.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE ReadList_Print()
!
! !USES:
!
    USE HCO_LOGFILE_MOD,  ONLY : HCO_PrintList
!
! !REVISION HISTORY:
!  20 Apr 2013 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC

    CHARACTER(LEN=255) :: MSG
    LOGICAL            :: verb

    ! ================================================================
    ! ReadList_Print begins here
    ! ================================================================

    verb = .TRUE.

    ! Print content of all lists
    IF ( ASSOCIATED(ReadLists) ) THEN 

       write(MSG,*) 'Content of one-time list:'
       CALL HCO_MSG(MSG,SEP1='=')
       CALL HCO_PrintList ( ReadLists%Once, verb )

       write(MSG,*) 'Content of year-list:'
       CALL HCO_MSG(MSG,SEP1='=')
       CALL HCO_PrintList ( ReadLists%Year, verb )

       write(MSG,*) 'Content of month-list:'
       CALL HCO_MSG(MSG,SEP1='=')
       CALL HCO_PrintList ( ReadLists%Month, verb )

       write(MSG,*) 'Content of day-list:'
       CALL HCO_MSG(MSG,SEP1='=')
       CALL HCO_PrintList ( ReadLists%Day, verb )

       write(MSG,*) 'Content of hour-list:'
       CALL HCO_MSG(MSG,SEP1='=')
       CALL HCO_PrintList ( ReadLists%Hour, verb )

       write(MSG,*) 'Content of always-to-read list:'
       CALL HCO_MSG(MSG,SEP1='=')
       CALL HCO_PrintList ( ReadLists%Always, verb )

    ELSE
       write(MSG,*) 'ReadList not defined yet!!'
       CALL HCO_MSG(MSG,SEP1='=')
    ENDIF

  END SUBROUTINE ReadList_Print
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ReadList_Cleanup
!
! !DESCRIPTION: Subroutine ReadList\_Cleanup removes all content of ReadList.
! If RemoveDct is set to True, the content of the data containers will be
! also removed, otherwise the corresponding pointer is just nullified. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE ReadList_Cleanup( RemoveDct )
!
! !USES:
!
    USE HCO_DataCont_Mod, ONLY : ListCont_Cleanup
!
! !INPUT PARAMETERS:
!
    LOGICAL, INTENT(IN)   :: RemoveDct
!
! !REVISION HISTORY:
!  20 Apr 2013 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC

    ! ================================================================
    ! ReadList_Cleanup begins here
    ! ================================================================
      
    IF ( ASSOCIATED(ReadLists) ) THEN 

       ! Remove all sublists in ReadList
       CALL ListCont_Cleanup ( ReadLists%Once,   RemoveDct )
       CALL ListCont_Cleanup ( ReadLists%Year,   RemoveDct )
       CALL ListCont_Cleanup ( ReadLists%Month,  RemoveDct )
       CALL ListCont_Cleanup ( ReadLists%Day,    RemoveDct )
       CALL ListCont_Cleanup ( ReadLists%Hour,   RemoveDct )
       CALL ListCont_Cleanup ( ReadLists%Always, RemoveDct )

       ! Remove ReadList 
       DEALLOCATE ( ReadLists )
    ENDIF
    ReadLists => NULL()

  END SUBROUTINE ReadList_Cleanup
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ReadList_to_EmisList 
!
! !DESCRIPTION: Subroutine ReadList\_to\_EmisList makes sure that the
! content of EmisList is up-to-date. Since ReadList and EmisList both
! point to the same data containers, this routine primarily makes sure
! that the data of containers with target IDs different that their 
! container IDs are properly set. See HCO\_EmisList\_MOD.F90 for more
! details.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE ReadList_to_EmisList( am_I_Root, HcoState, RC )
!
! !USES:
!
    USE HCO_Clock_Mod,    ONLY : HcoClock_First,    HcoClock_NewYear
    USE HCO_Clock_Mod,    ONLY : HcoClock_NewMonth, HcoClock_NewDay
    USE HCO_Clock_Mod,    ONLY : HcoClock_NewHour
    USE HCO_Emislist_Mod, ONLY : EmisList_Update
!
! !INPUT PARAMETERS:
!
    LOGICAL,         INTENT(IN   )  :: am_I_Root  ! root CPU?
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(HCO_State), POINTER        :: HcoState   ! Hemco state object
    INTEGER,         INTENT(INOUT)  :: RC         ! Return Code 
!
! !REVISION HISTORY:
!  09 Sep 2013 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC

    ! ================================================================
    ! ReadList_to_EmisList begins here
    ! ================================================================

    ! Enter
    CALL HCO_ENTER ('ReadList_to_EmisList (hco_readlist_mod.F90)', RC)
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Always add/update content from always-to-read list
    CALL EmisList_Update( am_I_Root, HcoState, ReadLists%Always, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Add content from one-time list on the first call 
    IF ( HcoClock_First() ) THEN
       CALL EmisList_Update( am_I_Root, HcoState, ReadLists%Once, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN
    ENDIF

    ! Add content from year-list if it's a new year
    IF ( HcoClock_NewYear() ) THEN
       CALL EmisList_Update( am_I_Root, HcoState, ReadLists%Year, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN
    ENDIF

    ! Add content from month-list if it's a new month
    IF ( HcoClock_NewMonth() ) THEN
       CALL EmisList_Update( am_I_Root, HcoState, ReadLists%Month, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN
    ENDIF

    ! Add content from day-list if it's a new day
    IF ( HcoClock_NewDay() ) THEN
       CALL EmisList_Update( am_I_Root, HcoState, ReadLists%Day, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN
    ENDIF

    ! Add content from hour-list if it's a new hour
    IF ( HcoClock_NewHour() ) THEN
       CALL EmisList_Update( am_I_Root, HcoState, ReadLists%Hour, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN
    ENDIF

    ! Leave w/ success
    CALL HCO_LEAVE ( RC )

  END SUBROUTINE ReadList_to_EmisList
!EOC
END MODULE HCO_ReadList_Mod
