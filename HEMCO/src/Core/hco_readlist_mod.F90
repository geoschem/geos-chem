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
!\\
!\\
! \begin{itemize}
! \item Year: update every year (annual data)
! \item Month: update every month (monthly data)
! \item Day: update every day (daily data)
! \item Hour: update every hour (hourly data)
! \item Hour3: update every 3 hours (3-hourly data)
! \item Once: update only once (time-invariant data)
! \item Always: update every time step
! \end{itemize}
!
! !INTERFACE:
!
MODULE HCO_ReadList_Mod
!
! !USES:
!
  USE HCO_Types_Mod
  USE HCO_Error_MOD
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
  PUBLIC  :: ReadList_Remove
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
  SUBROUTINE ReadList_Set( HcoState, Dct, RC )
!
! !USES:
!
    USE HCO_LOGFILE_MOD, ONLY : HCO_PrintDataCont
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(HCO_State),  POINTER        :: HcoState
    TYPE(DataCont),   POINTER        :: Dct
    INTEGER,          INTENT(INOUT)  :: RC
!
! !REVISION HISTORY:
!  20 Apr 2013 - C. Keller - Initial version
!  08 Aug 2018 - C. Keller - Add extra checks for 'exact' and 'range' fields.
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
    CALL HCO_ENTER (HcoState%Config%Err,'ReadList_Set (hco_readlist_mod.F90)', RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Verbose mode
    verb = HCO_IsVerb( HcoState%Config%Err, 2 )

    ! Add container to ReadList according to update freqency.
    ! Fields in list 'Hour' will be updated (i.e. re-read) every hour,
    ! fields in list 'Day' every day, etc.
    ! If a time range instead of a single time stamp is given,
    ! categorize the field according to the most rapidly changing time
    ! stamp.
    ! If no time attribute exist, put the field into the 'Once' list
    ! which reads the file only at the beginning. Similarly, fields
    ! with an update flag of 'always' will be put into the 'Always'
    ! list. The always update flag is set in routine HCO_ExtractTime,
    ! which is called from Config_ReadCont (hco_config_mod.F90).
    IF     ( Dct%Dta%UpdtFlag == HCO_UFLAG_ALWAYS ) THEN
       intv = 1
    ELSEIF ( Dct%Dta%UpdtFlag == HCO_UFLAG_3HR ) THEN
       intv = 7
    ELSEIF ( Dct%Dta%ncHrs(1) /= Dct%Dta%ncHrs(2) ) THEN
       intv = 2
    ELSEIF ( Dct%Dta%ncDys(1) /= Dct%Dta%ncDys(2) ) THEN
       intv = 3
    ELSEIF ( Dct%Dta%ncMts(1) /= Dct%Dta%ncMts(2) ) THEN
       intv = 4
    ELSEIF ( Dct%Dta%ncYrs(1) /= Dct%Dta%ncYrs(2) ) THEN
       intv = 5
    ELSE
       intv = 6
    ENDIF

    ! NOTE: In an ESMF environment, data I/O is organized through
    ! ESMF/MAPL. The hemco reading call (HCOIO_DATAREAD) sets a
    ! pointer of the data container array to the data array provided
    ! by MAPL. These arrays are already interpolated / updated (over
    ! time) by MAPL, and a pointer needs to be established only once.
    ! Hence, make sure that all containers are added to the one-time
    ! reading list!
    IF ( HcoState%Options%isESMF .AND. Dct%Dta%ncRead ) THEN
       intv = 6
    ENDIF

    ! Special handling of data with 'EXACT' or 'RANGE' flag:
    ! These data sets should be evaluated whenever the data
    ! changes (to see if it still falls within the selected
    ! time window).
    IF ( Dct%Dta%CycleFlag == HCO_CFLAG_RANGE .OR. &
         Dct%Dta%CycleFlag == HCO_CFLAG_EXACT ) THEN
       IF ( intv == 6 ) THEN
          IF ( Dct%Dta%ncHrs(1) == -1 .OR. &
               Dct%Dta%ncHrs(1) /= Dct%Dta%ncHrs(2) ) THEN
             intv = 2
          ELSEIF ( Dct%Dta%ncDys(1) == -1 .OR. &
               Dct%Dta%ncDys(1) /= Dct%Dta%ncDys(2) ) THEN
             intv = 3
          ELSEIF ( Dct%Dta%ncMts(1) == -1 .OR. &
               Dct%Dta%ncMts(1) /= Dct%Dta%ncMts(2) ) THEN
             intv = 4
          ELSEIF ( Dct%Dta%ncYrs(1) == -1 .OR. &
               Dct%Dta%ncYrs(1) /= Dct%Dta%ncYrs(2) ) THEN
             intv = 5
          ! Read always
          ELSE
             intv = 1
          ENDIF
       ENDIF
    ENDIF

    ! Add to ReadList according to interval flag
    IF (     intv == 1 ) THEN
       CALL DtCont_Add( HcoState%ReadLists%Always, Dct )
    ELSEIF ( intv == 2 ) THEN
       CALL DtCont_Add( HcoState%ReadLists%Hour,   Dct )
    ELSEIF ( intv == 3 ) THEN
       CALL DtCont_Add( HcoState%ReadLists%Day,    Dct )
    ELSEIF ( intv == 4 ) THEN
       CALL DtCont_Add( HcoState%ReadLists%Month,  Dct )
    ELSEIF ( intv == 5 ) THEN
       CALL DtCont_Add( HcoState%ReadLists%Year,   Dct )
    ELSEIF (intv == 7 ) THEN
       CALL DtCont_Add( HcoState%ReadLists%Hour3,  Dct )
    ELSE
       CALL DtCont_Add( HcoState%ReadLists%Once,   Dct )
    ENDIF

    ! Verbose
    IF ( Verb ) THEN
       WRITE(MSG,*) 'New container set to ReadList:'
       CALL HCO_MSG(HcoState%Config%Err, MSG)
       CALL HCO_PrintDataCont( HcoState, Dct, 3 )
    ENDIF

    ! Leave w/ success
    CALL HCO_LEAVE ( HcoState%Config%Err, RC )

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
  SUBROUTINE ReadList_Read( HcoState, RC, ReadAll )
!
! !USES:
!
    USE HCO_CLOCK_MOD, ONLY : HcoClock_NewYear
    USE HCO_CLOCK_MOD, ONLY : HcoClock_NewMonth
    USE HCO_CLOCK_MOD, ONLY : HcoClock_NewDay
    USE HCO_CLOCK_MOD, ONLY : HcoClock_NewHour
    USE HCO_CLOCK_MOD, ONLY : HcoClock_New3Hour
    USE HCO_CLOCK_MOD, ONLY : HcoClock_First
!
! !INPUT PARAMETERS:
!
    LOGICAL, OPTIONAL, INTENT(IN   )  :: ReadAll    ! read all fields?
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(HCO_State),   POINTER        :: HcoState   ! HEMCO state object
    INTEGER,           INTENT(INOUT)  :: RC         ! Success or failure?
!
! !REVISION HISTORY:
!  20 Apr 2013 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    LOGICAL            :: verb, RdAll
    CHARACTER(LEN=255) :: MSG

    ! ================================================================
    ! ReadList_Read begins here
    ! ================================================================

    ! For error handling
    CALL HCO_ENTER ( HcoState%Config%Err, &
                    'ReadList_Read (hco_readlist_mod.F90)', RC )
    IF ( RC /= HCO_SUCCESS ) THEN
       PRINT *, "Error in HCO_ENTER called from HEMCO ReadList_Read"
       RETURN
    ENDIF

    ! Verbose mode
    verb = HCO_IsVerb( HcoState%Config%Err, 1 )

    ! Read all fields?
    RdAll = .FALSE.
    IF ( PRESENT(ReadAll) ) RdAll = ReadAll
    ! Now use internal counter to determine first-time reading
    ! (ckeller, 02/07/2019).
    !IF ( HcoClock_First( HcoState%Clock, .FALSE. ) ) RdAll = .TRUE.
    IF ( HcoState%ReadLists%Counter == 0 ) RdAll = .TRUE.

    ! Read content from one-time list on the first call
    IF ( RdAll ) THEN
       IF ( Verb ) THEN
          WRITE(MSG,*) 'Now reading once list!'
          CALL HCO_MSG(HcoState%Config%Err,MSG)
       ENDIF
       CALL ReadList_Fill( HcoState, HcoState%ReadLists%Once, RC )
       IF ( RC /= HCO_SUCCESS ) THEN
          PRINT *, "Error in ReadList_Fill (1) called from HEMCO ReadList_Read"
          RETURN
       ENDIF
    ENDIF

    ! Read content from year list if it's a new year
    IF ( HcoClock_NewYear( HcoState%Clock, .FALSE. ) .OR. RdAll ) THEN
       IF ( Verb ) THEN
          WRITE(MSG,*) 'Now reading year list!'
          CALL HCO_MSG(HcoState%Config%Err,MSG)
       ENDIF
       CALL ReadList_Fill( HcoState, HcoState%ReadLists%Year, RC )
       IF ( RC /= HCO_SUCCESS ) THEN
          PRINT *, "Error in ReadList_Fill (2) called from HEMCO ReadList_Read"
          RETURN
       ENDIF
    ENDIF

    ! Read content from month list if it's a new month
    IF ( HcoClock_NewMonth( HcoState%Clock, .FALSE. ) .OR. RdAll ) THEN
       IF ( Verb ) THEN
          WRITE(MSG,*) 'Now reading month list!'
          CALL HCO_MSG(HcoState%Config%Err,MSG)
       ENDIF
       CALL ReadList_Fill( HcoState, HcoState%ReadLists%Month, RC )
       IF ( RC /= HCO_SUCCESS ) THEN
          PRINT *, "Error in ReadList_Fill (3) called from HEMCO ReadList_Read"
          RETURN
       ENDIF
    ENDIF

    ! Read content from day list if it's a new day
    IF ( HcoClock_NewDay( HcoState%Clock, .FALSE. ) .OR. RdAll ) THEN
       IF ( Verb ) THEN
          WRITE(MSG,*) 'Now reading day list!'
          CALL HCO_MSG(HcoState%Config%Err,MSG)
       ENDIF
       CALL ReadList_Fill( HcoState, HcoState%ReadLists%Day, RC )
       IF ( RC /= HCO_SUCCESS ) THEN
          PRINT *, "Error in ReadList_Fill (4) called from HEMCO ReadList_Read"
          RETURN
       ENDIF
    ENDIF

    ! Read content from hour list if it's a new hour
    IF ( HcoClock_NewHour( HcoState%Clock, .FALSE. ) .OR. RdAll ) THEN
       IF ( Verb ) THEN
          WRITE(MSG,*) 'Now reading hour list!'
          CALL HCO_MSG(HcoState%Config%Err,MSG)
       ENDIF
       CALL ReadList_Fill( HcoState, HcoState%ReadLists%Hour, RC )
       IF ( RC /= HCO_SUCCESS ) THEN
          PRINT *, "Error in ReadList_Fill (5) called from HEMCO ReadList_Read"
          RETURN
       ENDIF
    ENDIF

    ! Read content from 3-hour list if it's a new hour
    IF ( HcoClock_New3Hour( HcoState%Clock, .FALSE. ) .OR. RdAll ) THEN
       IF ( Verb ) THEN
          WRITE(MSG,*) 'Now reading 3-hour list!'
          CALL HCO_MSG(HcoState%Config%Err,MSG)
       ENDIF
       CALL ReadList_Fill( HcoState, HcoState%ReadLists%Hour3, RC )
       IF ( RC /= HCO_SUCCESS ) THEN
          PRINT *, "Error in ReadList_Fill (6) called from HEMCO ReadList_Read"
          RETURN
       ENDIF
    ENDIF

    ! Always add/update content from always-list
    IF ( Verb ) THEN
       WRITE(MSG,*) 'Now reading always list!'
       CALL HCO_MSG(HcoState%Config%Err,MSG)
    ENDIF
    CALL ReadList_Fill( HcoState, HcoState%ReadLists%Always, RC )
    IF ( RC /= HCO_SUCCESS ) THEN
       PRINT *, "Error in called ReadList_Fill (7) from HEMCO ReadList_Read"
       RETURN
    ENDIF

    ! Update counter
    HcoState%ReadLists%Counter = HcoState%ReadLists%Counter + 1

    ! Leave w/ success
    CALL HCO_LEAVE ( HcoState%Config%Err, RC )

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
! the ESMF/MAPL software framework. ReadLIst\_Fill provides the
! interface between HEMCO and the data reading interface. See module
! HCOI\_DATAREAD\_MOD.F90 for more details on data reading.
!\\
!\\
! The ReadList containers are added to EmisList immediately after data
! filling. This has the advantage that data arrays are immediately
! available through routine HCO\_GetPtr. This is required for country
! mappings that depend on the country mask input field.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE ReadList_Fill( HcoState, ReadList, RC )
!
! !USES:
!
    USE HCOIO_DataRead_Mod, ONLY : HCOIO_DataRead
    USE HCOIO_Read_Std_Mod, ONLY : HCOIO_ReadOther
    USE HCOIO_Read_Std_Mod, ONLY : HCOIO_CloseAll
    USE HCO_FileData_Mod,   ONLY : FileData_ArrIsDefined
    USE HCO_FileData_Mod,   ONLY : FileData_ArrIsTouched
    USE HCO_EmisList_Mod,   ONLY : EmisList_Pass
    USE HCO_DataCont_Mod,   ONLY : DataCont_Cleanup
    USE HCO_TIDX_MOD,       ONLY : tIDx_Assign
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
!  23 Dec 2014 - C. Keller - Now pass container to EmisList immediately after
!                            reading the data. Added second loop to remove
!                            data arrays that are not used in EmisList.
!  02 Feb 2015 - C. Keller - Now call tIDx_Assign here instead of in
!                            hco_emislist_mod. This way, hco_emislist_mod
!                            can also be used by hco_clock_mod.
!  24 Mar 2015 - C. Keller - Now avoid closing/reopening the same file all
!                            the time.
!  24 Mar 2016 - C. Keller - Remove GetFileLUN and SaveFileLUN. This is now
!                            handled in hcoio_read_std_mod.F90.
!  26 Oct 2016 - R. Yantosca - Don't nullify local ptrs in declaration stmts
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    TYPE(ListCont), POINTER  :: Lct
    LOGICAL                  :: verb
    CHARACTER(LEN=255)       :: MSG

    ! ================================================================
    ! ReadList_Fill begins here
    ! ================================================================

    ! For error handling
    CALL HCO_ENTER (HcoState%Config%Err,&
                   'ReadList_Fill (hco_readlist_mod.F90)', RC )
    IF ( RC /= HCO_SUCCESS ) THEN
       PRINT *, "Error in HCO_ENTER called from HEMCO ReadList_Fill"
       RETURN
    ENDIF

    ! Verbose mode?
    verb = HCO_IsVerb ( HcoState%Config%Err, 2 )

    ! Loop over all containers
    Lct => ReadList
    DO WHILE ( ASSOCIATED ( Lct ) )

       ! Check if data has already been touched. Multiple data
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
       IF ( Lct%Dct%DtaHome < 0 ) THEN
          IF ( FileData_ArrIsTouched(Lct%Dct%Dta) ) THEN
             Lct%Dct%DtaHome     = 0
             Lct%Dct%Dta%DoShare = .TRUE.
          ELSE
             Lct%Dct%DtaHome = 1
          ENDIF
       ENDIF

       ! Read if this is the home container
       IF ( Lct%Dct%DtaHome == 1 ) THEN

          ! Read from other source if it's not a netCDF file
          IF ( .NOT. Lct%Dct%Dta%NcRead ) THEN
             CALL HCOIO_ReadOther( HcoState, Lct, RC )
             IF ( RC /= HCO_SUCCESS ) THEN
                PRINT *, "Error in HCOIO_ReadOther called from HEMCO ReadList_Fill: ", TRIM(Lct%Dct%cname)
                RETURN
             ENDIF

          ! Read from netCDF file otherwise
          ELSE

             ! Read data
             CALL HCOIO_DATAREAD( HcoState, Lct, RC )
             IF ( RC /= HCO_SUCCESS ) THEN
                PRINT *, "Error in HCOIO_DATAREAD called from HEMCO ReadList_Fill: ", TRIM(Lct%Dct%cname)
                RETURN
             ENDIF
          ENDIF

          ! We now have touched this data container
          Lct%Dct%Dta%IsTouched = .TRUE.

       ENDIF

       ! Pass container to EmisList (only if array is defined)
       IF ( FileData_ArrIsDefined(Lct%Dct%Dta) ) THEN

          ! Set time index pointer tIDx of this data container.
          ! tIDx will be set according to the number of time slices
          ! (and the tim einterval between them) hold by this data
          ! container. For hourly data (24 time slices), for example,
          ! tIDx will point to the corresponding 'HOURLY' or
          ! 'HOURLY_GRID' time index collection type defined in
          ! hco_tidx_mod.
          CALL tIDx_Assign ( HcoState, Lct%Dct, RC )
          IF ( RC /= HCO_SUCCESS ) THEN
             PRINT *, "Error in tIDx_Assign called from HEMCO ReadList_Fill: ", TRIM(Lct%Dct%cname)
             RETURN
          ENDIF

          ! Container is now read to be passed to emissions list.
          CALL EmisList_Pass( HcoState, Lct, RC )
          IF ( RC /= HCO_SUCCESS ) THEN
             PRINT *, "Error in EmisList_Pass called from HEMCO ReadList_Fill: ", TRIM(Lct%Dct%cname)
             RETURN
          ENDIF

       ENDIF

       ! Advance to next container
       Lct => Lct%NextCont
    ENDDO

    ! Make sure that all netCDF files are closed
    CALL HCOIO_CloseAll ( HcoState, RC )
    IF ( RC /= HCO_SUCCESS ) THEN
       PRINT *, "Error in HCOIO_CloseAll called from HEMCO ReadList_Fill: ", TRIM(Lct%Dct%cname)
       RETURN
    ENDIF

    ! Second loop to clean up all data that is not used in EmisList.
    ! This cannot be done within EmisList_Pass since it is possible
    ! that some container data is used by multiple containers.
    Lct => ReadList
    DO WHILE ( ASSOCIATED(Lct) )

       ! Remove array if not used in the emissions list
       IF ( .NOT. Lct%Dct%Dta%IsInList ) THEN
          CALL DataCont_Cleanup( Lct%Dct, ArrOnly=.TRUE. )

          ! Verbose mode
          IF ( verb ) THEN
             MSG = 'Remove data array of ' // TRIM(Lct%Dct%cName)
             CALL HCO_MSG( HcoState%Config%Err, MSG )
          ENDIF
       ENDIF

       ! Advance in list
       Lct => Lct%NextCont
    ENDDO

    ! Cleanup
    Lct => NULL()

    ! Leave with success
    CALL HCO_LEAVE ( HcoState%Config%Err, RC )

  END SUBROUTINE ReadList_Fill
!EOC
!------------------------------------------------------------------------------
!         Harvard-NASA Emissions Component (HEMCO)                   !
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
!  26 Oct 2016 - R. Yantosca - Don't nullify local ptrs in declaration stmts
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    TYPE(ListCont), POINTER :: Lct
    TYPE(ListCont), POINTER :: TmpLct
    INTEGER                 :: cID

    ! ================================================================
    ! DtCont_Add begins here
    ! ================================================================

    ! Init
    Lct    => NULL()
    TmpLct => NULL()

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
  SUBROUTINE ReadList_Init( ReadLists, RC )
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(RdList),     POINTER        :: ReadLists
    INTEGER,          INTENT(INOUT)  :: RC
!
! !REVISION HISTORY:
!  20 Apr 2013 - C. Keller - Initial version
!  07 Feb 2019 - C. Keller - Added counter
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

    ALLOCATE ( ReadLists%Hour3 )
    NULLIFY ( ReadLists%Hour3  )

    ALLOCATE ( ReadLists%Always )
    NULLIFY ( ReadLists%Always  )

    ! No file in buffer yet
    ReadLists%FileInArchive = ''
    ReadLists%FileLun       = -1

    ! Initialize counter
    ReadLists%Counter       = 0

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
  SUBROUTINE ReadList_Print( HcoState, ReadLists, verb )
!
! !USES:
!
    USE HCO_LOGFILE_MOD,  ONLY : HCO_PrintList
!
! !INPUT ARGUMENTS
!
    TYPE(HCO_State), POINTER       :: HcoState
    TYPE(RdList),    POINTER       :: ReadLists
    INTEGER,         INTENT(IN)    :: verb   ! verbose number
!
! !REVISION HISTORY:
!  20 Apr 2013 - C. Keller - Initial version
!  15 Mar 2015 - C. Keller - Added verbose number as input argument
!EOP
!------------------------------------------------------------------------------
!BOC

    CHARACTER(LEN=255) :: MSG

    ! ================================================================
    ! ReadList_Print begins here
    ! ================================================================

    ! Nothing to do if HEMCO verbose level is below passed verbose number
    IF ( .NOT. HCO_IsVerb(HcoState%Config%Err,verb) ) RETURN

    ! Print content of all lists
    IF ( ASSOCIATED(ReadLists) .and. HcoState%amIRoot ) THEN

       WRITE(MSG,*) 'Contents of one-time list:'
       CALL HCO_MSG(HcoState%Config%Err,MSG,SEP1='=')
       CALL HCO_PrintList ( HcoState, ReadLists%Once, verb )

       WRITE(MSG,*) 'Contents of year list:'
       CALL HCO_MSG(HcoState%Config%Err,MSG,SEP1='=')
       CALL HCO_PrintList ( HcoState, ReadLists%Year, verb )

       WRITE(MSG,*) 'Contents of month list:'
       CALL HCO_MSG(HcoState%Config%Err,MSG,SEP1='=')
       CALL HCO_PrintList ( HcoState, ReadLists%Month, verb )

       WRITE(MSG,*) 'Contents of day list:'
       CALL HCO_MSG(HcoState%Config%Err,MSG,SEP1='=')
       CALL HCO_PrintList ( HcoState, ReadLists%Day, verb )

       WRITE(MSG,*) 'Contents of 3-hour list:'
       CALL HCO_MSG(HcoState%Config%Err,MSG,SEP1='=')
       CALL HCO_PrintList ( HcoState, ReadLists%Hour3, verb )

       WRITE(MSG,*) 'Contents of hour list:'
       CALL HCO_MSG(HcoState%Config%Err,MSG,SEP1='=')
       CALL HCO_PrintList ( HcoState, ReadLists%Hour, verb )

       WRITE(MSG,*) 'Contents of always-to-read list:'
       CALL HCO_MSG(HcoState%Config%Err,MSG,SEP1='=')
       CALL HCO_PrintList ( HcoState, ReadLists%Always, verb )

    ELSE
       WRITE(MSG,*) 'ReadList not defined yet!!'
       CALL HCO_MSG(HcoState%Config%Err,MSG,SEP1='=')
    ENDIF

  END SUBROUTINE ReadList_Print
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ReadList_Remove
!
! !DESCRIPTION: Subroutine ReadList\_Remove removes the container given by
! name from the ReadList. If no container with the given name exist, nothing
! is done. This routine returns an error if the container already holds data.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE ReadList_Remove( HcoState, cName, RC )
!
! !USES:
!
!
! !INPUT PARAMETERS:
!
    TYPE(HCO_State),  POINTER         :: HcoState
    CHARACTER(LEN=*), INTENT(IN   )   :: cName
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(INOUT)   :: RC
!
! !REVISION HISTORY:
!  13 Jan 2015 - C. Keller - Initial version
!  26 Oct 2016 - R. Yantosca - Don't nullify local ptrs in declaration stmts
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER                   :: I
    TYPE(ListCont), POINTER   :: This
    TYPE(ListCont), POINTER   :: Prev
    TYPE(ListCont), POINTER   :: Next
    LOGICAL                   :: FOUND
    CHARACTER(LEN=255)        :: MSG
    CHARACTER(LEN=255)        :: LOC = 'ReadList_Remove (HCO_ReadList_Mod.F90)'

    ! ================================================================
    ! ReadList_Remove begins here
    ! ================================================================

    ! Assume success until otherwise
    RC = HCO_SUCCESS
    IF ( .NOT. ASSOCIATED(HcoState%ReadLists) ) RETURN

    ! Init
    This => NULL()
    Prev => NULL()
    Next => NULL()

    ! Search for the given container
    DO I = 1,7

       ! Select list to be used
       IF ( I == 1 ) This => HcoState%ReadLists%Once
       IF ( I == 2 ) This => HcoState%ReadLists%Year
       IF ( I == 3 ) This => HcoState%ReadLists%Month
       IF ( I == 4 ) This => HcoState%ReadLists%Day
       IF ( I == 5 ) This => HcoState%ReadLists%Hour
       IF ( I == 6 ) This => HcoState%ReadLists%Always
       IF ( I == 7 ) This => HcoState%ReadLists%Hour3

       ! Initialize working variables
       FOUND = .FALSE.
       Prev => This

       ! Walk through list, looking for this container
       DO WHILE ( ASSOCIATED(This) )

          ! Next container in list
          Next => This%NextCont

          ! Is that the container of interest?
          IF ( TRIM(This%Dct%cName) == TRIM(cName) ) THEN
             FOUND = .TRUE.
             EXIT
          ENDIF

          ! Advance
          Prev => This
          This => Next
       ENDDO !This

       ! Advance if not found
       IF ( .NOT. FOUND ) CYCLE

       ! Check first if data has already been read. In this case, the data home
       ! flag is updated.
       IF ( This%Dct%DtaHome >= 0 ) THEN
          MSG = 'Cannot remove from ReadList. Data has already been read: ' // &
                TRIM(This%Dct%cName)
          CALL HCO_ERROR( HcoState%Config%Err, MSG, RC, THISLOC = LOC )
       ENDIF

       ! Connect previous container to next container in list:
       ! - Special case that this is the first container in the list
       IF ( Prev%Dct%cID == This%Dct%cID ) THEN
          IF ( I == 1 ) HcoState%ReadLists%Once   => Next
          IF ( I == 2 ) HcoState%ReadLists%Year   => Next
          IF ( I == 3 ) HcoState%ReadLists%Month  => Next
          IF ( I == 4 ) HcoState%ReadLists%Day    => Next
          IF ( I == 5 ) HcoState%ReadLists%Hour   => Next
          IF ( I == 6 ) HcoState%ReadLists%Always => Next
          IF ( I == 7 ) HcoState%ReadLists%Hour3  => Next

       ! - Otherwise, just pop out this container from list
       ELSE
          Prev%NextCont => Next
       ENDIF

       ! Remove pointer to data container, detach this container from list
       This%Dct      => NULL()
       This%NextCont => NULL()

       ! Deallocate this container
       DEALLOCATE(This)

       ! If we make it to here, we have successfully removed the container and
       ! don't need to cycle thorugh the loop any more
       EXIT

    ENDDO !

    ! Free pointer
    This => NULL()
    Prev => NULL()
    Next => NULL()

    ! Return w/ success
    RC = HCO_SUCCESS

  END SUBROUTINE ReadList_Remove
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
  SUBROUTINE ReadList_Cleanup( ReadLists, RemoveDct )
!
! !USES:
!
    USE HCO_DataCont_Mod, ONLY : ListCont_Cleanup
!
! !INPUT PARAMETERS:
!

    LOGICAL,          INTENT(IN   )   :: RemoveDct
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(RdList),     POINTER         :: ReadLists
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
       CALL ListCont_Cleanup ( ReadLists%Hour3,  RemoveDct )
       CALL ListCont_Cleanup ( ReadLists%Always, RemoveDct )

       ! Remove ReadList
       DEALLOCATE ( ReadLists )
    ENDIF
    ReadLists => NULL()

  END SUBROUTINE ReadList_Cleanup
!EOC
END MODULE HCO_ReadList_Mod
