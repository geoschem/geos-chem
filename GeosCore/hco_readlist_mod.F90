!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: hco_readlist_mod
!
! !DESCRIPTION: Module HCO\_READLIST\_MOD contains routines, variables
! and types to define, write and modify data for the HEMCO emissions list. 
! The reading information of all used fields (base emissions, scale factors
! and masks) are saved in a reading list. Different lists are available
! for files that have to be updated every year, month, day, or hour, and
! for files that have to be read only once. 
! For each time step and depending on the current simulation time, these
! list information are collected into a temporary reading list. For
! example, if we enter a new month, the temporary list will contain all
! information listed in the monthly, daily and hourly list (but not of
! the annual list!), as well as information in the one-time list.\\
! In a second step, data is read, regridded and added to the emissions
! lists using the information in the temporary lists. This list is
! cleared after each time step, as is the one-time only reading list.
! 
! \\
! !INTERFACE: 
!
      MODULE HCO_READLIST_MOD
!
! !USES:
!
      USE HCO_ERROR_MOD
      USE HCO_TYPE_MOD,          ONLY : HCO_State, RdCont 
 
      IMPLICIT NONE

      PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC :: ReadData
      PUBLIC :: ReadList_Init
      PUBLIC :: ReadList_Set
      PUBLIC :: ReadList_Print
      PUBLIC :: ReadList_Cleanup
      PUBLIC :: Write_EmisLL
!
! !PRIVATE MEMBER FUNCTIONS:
!
      PRIVATE :: RdCont_Init
      PRIVATE :: RdCont_Add
      PRIVATE :: RdCont_Cleanup
      PRIVATE :: ReadList_Fill
      PRIVATE :: List_Print
      PRIVATE :: List_Cleanup
      PRIVATE :: ReadList2EmisLL
      PRIVATE :: AddCont2EmisLL
      PRIVATE :: AssignScalIDs
!
! !REVISION HISTORY:
!  20 Apr 2013 - C. Keller: Initial version
!
!EOP
!-----------------------------------------------------------------------------
!BOC
!
! !MODULE TYPES/VARIABLES:
!
      TYPE RdList
         TYPE(RdCont), POINTER    :: Once
         TYPE(RdCont), POINTER    :: Always
         TYPE(RdCont), POINTER    :: Year
         TYPE(RdCont), POINTER    :: Month
         TYPE(RdCont), POINTER    :: Day
         TYPE(RdCont), POINTER    :: Hour
      END TYPE RdList

      ! Internally used ReadList
      TYPE(RdList), POINTER       :: ReadList => NULL()
!
! !MODULE INTERFACES: 
!
      INTERFACE ReadList_Set
         MODULE PROCEDURE ReadList_Set_Base
         MODULE PROCEDURE ReadList_Set_Scal
      END INTERFACE

      CONTAINS
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: ReadList_Set_Base
!
! !DESCRIPTION: Subroutine ReadList\_Set\_Base sets a new ReadList container 
! entry to ReadList, containing information on a base emissions field. 
! The container attributes are given the passed arguments. Default values are 
! assigned to all attributes that are not explicitly specified.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE ReadList_Set_Base ( am_I_Root,  HcoState,           &
                                     DataType,   cName,  ncFile,      & 
                                     ncPara,     ncYrs,  ncMts,       &
                                     ncDys,      ncHrs,  TrcID,       &
                                     UseScalIDs, Oper,   Cat,         & 
                                     Hier,       ncRead, ESMF_Dim,    & 
                                     ESMF_Unit,  Add,    ExtNr,       &
                                     RC                                )
!
! !USES:
!
      USE HCO_EMISLL_MOD,         ONLY : SclMax
!
! !ARGUMENTS:
!
      LOGICAL,          INTENT(IN   )         :: am_I_Root
      TYPE(HCO_State),  POINTER               :: HcoState
      INTEGER,          INTENT(IN)            :: DataType
      CHARACTER(LEN=*), INTENT(IN)            :: cName
      CHARACTER(LEN=*), INTENT(IN)            :: ncFile
      CHARACTER(LEN=*), INTENT(IN)            :: ncPara
      INTEGER,          INTENT(IN)            :: ncYrs(2)
      INTEGER,          INTENT(IN)            :: ncMts(2)
      INTEGER,          INTENT(IN)            :: ncDys(2)
      INTEGER,          INTENT(IN)            :: ncHrs(2)
      INTEGER,          INTENT(IN)            :: TrcID
      INTEGER,          INTENT(IN)            :: UseScalIDs(SclMax)
      INTEGER,          INTENT(IN)            :: Oper
      INTEGER,          INTENT(IN)            :: Cat
      INTEGER,          INTENT(IN)            :: Hier
      INTEGER,          INTENT(IN)            :: ExtNr
      LOGICAL,          INTENT(IN)            :: ncRead
      CHARACTER(LEN=*), INTENT(IN)            :: ESMF_Dim
      CHARACTER(LEN=*), INTENT(IN)            :: ESMF_Unit
      LOGICAL,          INTENT(IN)            :: Add
      INTEGER,          INTENT(INOUT)         :: RC
!
! !REVISION HISTORY:
!  20 Apr 2013 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      TYPE(RdCont), POINTER     :: newCont    => NULL() 
      INTEGER                   :: append
      LOGICAL                   :: verb

      ! ================================================================
      ! ReadList_Set_Base begins here
      ! ================================================================

      ! Init
      verb = am_I_Root .AND. HcoState%verbose

      ! Set append flag 
      IF ( Add ) THEN
         append = 1 
      ELSE
         append = 0
      ENDIF

      ! Add an empty container to the correct ReadList.
      ! Define here into which ReadList this field belongs.
      ! Fields in list 'Hour' will be updated (i.e. re-read) every hour, 
      ! fields in list 'Day' every day, etc.
      ! If a time range instead of a single time stamp is given,
      ! categorize the field according to the most rapidly changing time
      ! stamp. 
      ! If no dynamic attribute exist, put the field into the 'Once'
      ! list which reads the file and then forgets about it.
#if defined(ESMF_)
      CALL RdCont_Add( ReadList%Always, newCont, append ) 
#else      
      IF (     ncHrs(1) /= ncHrs(2) ) THEN
         CALL RdCont_Add( ReadList%Hour, newCont, append ) 
      ELSEIF ( ncDys(1) /= ncDys(2) ) THEN
         CALL RdCont_Add( ReadList%Day, newCont, append ) 
      ELSEIF ( ncMts(1) /= ncMts(2) ) THEN
         CALL RdCont_Add( ReadList%Month, newCont, append ) 
      ELSEIF ( ncYrs(1) /= ncYrs(2) ) THEN
         CALL RdCont_Add( ReadList%Year, newCont, append ) 
      ELSE
         CALL RdCont_Add( ReadList%Once, newCont, append ) 
      ENDIF
#endif

      ! Pass container arguments
      newCont%DataType  = DataType
      newCont%ExtNr     = ExtNr
      newCont%cName     = cName
      newCont%ncFile    = ncFile
      newCont%ncPara    = ncPara
      newCont%ncYrs     = ncYrs
      newCont%ncMts     = ncMts
      newCont%ncDys     = ncDys
      newCont%ncHrs     = ncHrs
      newCont%TrcID     = TrcID
      newCont%Oper      = Oper
      newCont%Cat       = Cat
      newCont%Hier      = Hier
      newCont%ncRead    = ncRead
      newCont%ESMF_Dim  = ESMF_Dim
      newCont%ESMF_Unit = ESMF_Unit
      newCont%Add       = Add

      ! Init values
      newCont%ScalID    = -1
      newCont%cID       = -1

      ! Pass scale factor IDs to be used
      ALLOCATE ( newCont%UseScalIDs(SclMax) )
      NewCont%UseScalIDs(:) = UseScalIDs(:)

      ! Verbose mode 
      IF ( verb ) THEN
         write(6,*) 'New base data set to ReadList:'
         write(6,*) 'Extension nr  : ', newCont%ExtNr
         write(6,*) 'Container name: ', TRIM(newCont%cName)
         write(6,*) 'ncFile        : ', TRIM(newCont%ncFile)
         write(6,*) 'ncPara        : ', TRIM(newCont%ncPara)
         write(6,*) 'ncYr          : ', newCont%ncYrs
         write(6,*) 'ncMt          : ', newCont%ncMts
         write(6,*) 'ncDy          : ', newCont%ncDys
         write(6,*) 'ncHr          : ', newCont%ncHrs
         write(6,*) 'ncRead?         ', newCont%ncRead
         write(6,*) 'Category      : ', newCont%Cat
         write(6,*) 'Hierarchy     : ', newCont%Hier
         write(6,*) 'Scale factors : ', newCont%UseScalIDs
         write(6,*) 'Add field?    : ', newCont%Add
      ENDIF

      ! Leave w/ success
      RC = HCO_SUCCESS

      END SUBROUTINE ReadList_Set_Base
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: ReadList_Set_Scal
!
! !DESCRIPTION: Subroutine ReadList\_Set\_Scal sets a new ReadList container 
! entry to ReadList, containing information on a scale field (scale
! factor or mask). 
! The container attributes are given the passed arguments. Default values are 
! assigned to all attributes that are not explicitly specified.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE ReadList_Set_Scal ( am_I_Root,  HcoState,  DataType, &
                                     cName,      ncFile,     ncPara,   &
                                     ncYrs,      ncMts,      ncDys,    &
                                     ncHrs,      ScalID,     Oper,     &
                                     ncRead,     ESMF_Dim,   ESMF_Unit,&
                                     Add,        RC,         Array     )
!
! !ARGUMENTS:
!
      LOGICAL,          INTENT(IN   )         :: am_I_Root
      TYPE(HCO_State),  POINTER               :: HcoState
      INTEGER,          INTENT(IN)            :: DataType
      CHARACTER(LEN=*), INTENT(IN)            :: cName
      CHARACTER(LEN=*), INTENT(IN)            :: ncFile
      CHARACTER(LEN=*), INTENT(IN)            :: ncPara
      INTEGER,          INTENT(IN)            :: ncYrs(2)
      INTEGER,          INTENT(IN)            :: ncMts(2)
      INTEGER,          INTENT(IN)            :: ncDys(2)
      INTEGER,          INTENT(IN)            :: ncHrs(2)
      INTEGER,          INTENT(IN)            :: ScalID
      INTEGER,          INTENT(IN)            :: Oper
      LOGICAL,          INTENT(IN)            :: ncRead
      CHARACTER(LEN=*), INTENT(IN)            :: ESMF_Dim
      CHARACTER(LEN=*), INTENT(IN)            :: ESMF_Unit
      LOGICAL,          INTENT(IN)            :: Add
      INTEGER,          INTENT(INOUT)         :: RC
      REAL*8,           POINTER,    OPTIONAL  :: Array(:,:,:,:)
!
! !REVISION HISTORY:
!  20 Apr 2013 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      TYPE(RdCont), POINTER     :: newCont    => NULL() 
      INTEGER                   :: append
      INTEGER                   :: AS
      CHARACTER(LEN=255)        :: MSG, LOC
      LOGICAL                   :: verb

      ! ================================================================
      ! ReadList_Set_Scal begins here
      ! ================================================================

      ! Init
      LOC = 'ReadList_Set_Scal'
      verb = am_I_Root .AND. HcoState%verbose

      ! Set append flag 
      IF ( Add ) THEN
         append = 1 
      ELSE
         append = 0
      ENDIF

      ! Add an empty container to the correct ReadList.
      ! Define here into which ReadList this field belongs.
      ! Fields in list 'Hour' will be updated (i.e. re-read) every hour, 
      ! fields in list 'Day' every day, etc.
      ! If a time range instead of a single time stamp is given,
      ! categorize the field according to the most rapidly changing time
      ! stamp. 
      ! If no dynamic attribute exist, put the field into the 'Once'
      ! list which reads the file and then forgets about it.
#if defined(ESMF_)
      CALL RdCont_Add( ReadList%Always, newCont, append ) 
#else      
      IF (     ncHrs(1) /= ncHrs(2) ) THEN
         CALL RdCont_Add( ReadList%Hour, newCont, append ) 
      ELSEIF ( ncDys(1) /= ncDys(2) ) THEN
         CALL RdCont_Add( ReadList%Day, newCont, append ) 
      ELSEIF ( ncMts(1) /= ncMts(2) ) THEN
         CALL RdCont_Add( ReadList%Month, newCont, append ) 
      ELSEIF ( ncYrs(1) /= ncYrs(2)  ) THEN
         CALL RdCont_Add( ReadList%Year, newCont, append ) 
      ELSE
         CALL RdCont_Add( ReadList%Once, newCont, append ) 
      ENDIF
#endif

      ! Pass container arguments
      newCont%DataType  = DataType 
      newCont%cName     = cName
      newCont%ncFile    = ncFile
      newCont%ncPara    = ncPara
      newCont%ncYrs     = ncYrs
      newCont%ncMts     = ncMts
      newCont%ncDys     = ncDys
      newCont%ncHrs     = ncHrs
      newCont%Oper      = Oper
      newCont%ScalID    = ScalID
      newCont%ncRead    = ncRead
      newCont%Add       = Add
      newCont%ESMF_Dim  = ESMF_Dim
      newCont%ESMF_Unit = ESMF_Unit

      ! Default values
      newCont%ExtNr = -1
      newCont%TrcID = -1
      newCont%Cat   = -999
      newCont%Hier  = -999
      newCont%cID   = -1

      ! For the moment, also allow to pass array directly. In this case, 
      ! the ncRead flag must be disabled!!
      IF ( PRESENT ( Array ) ) THEN
         IF ( newCont%ncRead ) THEN
            MSG = 'Wrong ncRead flag in ' // TRIM(cName)
            CALL HCO_ERROR ( MSG, LOC, RC ); RETURN 
         ENDIF

         ! Allocate array in newCont
         ALLOCATE( newCont%Array( SIZE(Array,1),SIZE(Array,2), & 
                                  SIZE(Array,3),SIZE(Array,4)) )

         ! Pass array
         newCont%Array(:,:,:,:) = Array(:,:,:,:)
      ENDIF

      ! Verbose mode 
      IF ( Verb ) THEN
         write(6,*) 'New scale data set to ReadList:'
         write(6,*) 'Container name: ', TRIM(newCont%cName)
         write(6,*) 'Data type     : ', newCont%DataType
         write(6,*) 'ncFile        : ', TRIM(newCont%ncFile)
         write(6,*) 'ncPara        : ', TRIM(newCont%ncPara)
         write(6,*) 'ncYr          : ', newCont%ncYrs
         write(6,*) 'ncMt          : ', newCont%ncMts
         write(6,*) 'ncDy          : ', newCont%ncDys
         write(6,*) 'ncHr          : ', newCont%ncHrs
         write(6,*) 'ncRead?         ', newCont%ncRead
         write(6,*) 'Scale ID      : ', newCont%ScalID
      ENDIF

      ! Leave w/ success
      RC = HCO_SUCCESS 

      END SUBROUTINE ReadList_Set_Scal
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: ReadData 
!
! !DESCRIPTION: Subroutine ReadData makes sure that all arrays in the
! reading lists are up to date.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE ReadData ( am_I_Root, HcoState, Clock, RC ) 
!
! !USES:
!
      USE HCO_TIME_MOD,     ONLY : HcoClock
!
! !ARGUMENTS:
!
      LOGICAL,         INTENT(IN   )  :: am_I_Root  ! root CPU?
      TYPE(HCO_State), POINTER        :: HcoState
      TYPE(HcoClock),  POINTER        :: Clock      ! Simulation clock 
      INTEGER,         INTENT(INOUT)  :: RC         ! Success or failure?
!
!
! !REVISION HISTORY:
!  20 Apr 2013 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!

      ! ================================================================
      ! ReadData begins here
      ! ================================================================

      ! Add content from one-time list on the first call 
      IF ( Clock%PrevYear < 0 ) THEN
         CALL ReadList_Fill ( am_I_Root,     HcoState, & 
                              ReadList%Once, Clock,     RC )
         IF ( RC /= HCO_SUCCESS ) RETURN 
      ENDIF

      ! Add content from year-list if it's a new year
      IF ( Clock%PrevYear /= Clock%ThisYear ) THEN
         CALL ReadList_Fill ( am_I_Root,     HcoState, &
                              ReadList%Year, Clock,     RC ) 
         IF ( RC /= HCO_SUCCESS ) RETURN 
      ENDIF

      ! Add content from month-list if it's a new month
      IF ( Clock%PrevMonth /= Clock%ThisMonth ) THEN
         CALL ReadList_Fill ( am_I_Root,      HcoState, & 
                              ReadList%Month, Clock,     RC ) 
         IF ( RC /= HCO_SUCCESS ) RETURN 
      ENDIF

      ! Add content from day-list if it's a new year
      IF ( Clock%PrevDay /= Clock%ThisDay ) THEN
         CALL ReadList_Fill ( am_I_Root,     HcoState, &
                              ReadList%Day,  Clock,     RC ) 
         IF ( RC /= HCO_SUCCESS ) RETURN 
      ENDIF

      ! Add content from hour-list if it's a new year
      IF ( Clock%PrevHour /= Clock%ThisHour ) THEN
         CALL ReadList_Fill ( am_I_Root,     HcoState, &
                              ReadList%Hour, Clock,     RC ) 
         IF ( RC /= HCO_SUCCESS ) RETURN 
      ENDIF

      ! Always add/update content from always-to-read list
      CALL ReadList_Fill ( am_I_Root,       HcoState, &
                           ReadList%Always, Clock,     RC ) 
      IF ( RC /= HCO_SUCCESS ) RETURN 

      ! Leave w/ success
      RC = HCO_SUCCESS 

      END SUBROUTINE ReadData
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: ReadList_Fill
!
! !DESCRIPTION: Subroutine ReadList_Fill fills the current reading list by 
! reading all necessary data from netCDF.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE ReadList_Fill ( am_I_Root, HcoState, & 
                                 TempList,  Clock,     RC ) 
!
! !USES:
!
      USE HCOI_DATAREAD_MOD,     ONLY : HCOI_DATAREAD
      USE HCO_TIME_MOD,          ONLY : HcoClock
!
! !ARGUMENTS:
!
      LOGICAL,         INTENT(IN   )  :: am_I_Root  ! root CPU?
      TYPE(HCO_State), POINTER        :: HcoState
      TYPE(RdCont),    POINTER        :: TempList   ! Current reading list
      TYPE(HcoClock),  POINTER        :: Clock      ! Emission time clock 
      INTEGER,         INTENT(INOUT)  :: RC         ! Success or failure?
!
!
! !REVISION HISTORY:
!  20 Apr 2013 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      TYPE(RdCont),     POINTER  :: ThisRdCont => NULL()

      ! ================================================================
      ! ReadList_Fill begins here
      ! ================================================================

      ! Point to first container in TempList 
      ThisRdCont => TempList

      ! Loop over all containers
      DO WHILE ( ASSOCIATED ( ThisRdCont ) ) 

         ! Read (and regrid) data if the ncRead flag is on
         IF ( ThisRdCont%ncRead ) THEN

            ! Read data
            CALL HCOI_DATAREAD ( am_I_Root,  HcoState, &
                                 ThisRdCont, Clock,  RC )
            IF ( RC /= HCO_SUCCESS ) RETURN
         ENDIF

         ! Point to next container
         ThisRdCont => ThisRdCont%NextRdCont
      ENDDO

      ! Cleanup
      ThisRdCont => NULL()

      ! Leave with success 
      RC = HCO_SUCCESS

      END SUBROUTINE ReadList_Fill
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: RdCont_Add
!
! !DESCRIPTION: Subroutine RdCont_Add adds a new container to the specified
! reading list. 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE RdCont_Add ( List, newCont, append )
!
! !USES:
!
! !ARGUMENTS:
!
      TYPE(RdCont),   POINTER      :: List 
      TYPE(RdCont),   POINTER      :: newCont
      INTEGER,        INTENT(IN)   :: append
!
!
! !REVISION HISTORY:
!  20 Apr 2013 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      TYPE(RdCont), POINTER     :: tmpCont

      ! ================================================================
      ! RdCont_Add begins here
      ! ================================================================

      ! Init
      NULLIFY ( tmpCont ) 

      ! Initalize the new pool container
      CALL RdCont_Init ( newCont )

      ! Set new pool container to head if list is not defined yet.
      IF ( .NOT. ASSOCIATED ( List ) ) THEN
         List => newCont
     
      ! If list is already defined 
      ELSE 

         ! Set temporary pointer to head of list to start with
         tmpCont => List

         ! Append container to list if append argument is set to 1
         IF ( append == 1 ) THEN

            ! Advance to end of list
            DO WHILE ( ASSOCIATED ( tmpCont%NextRdCont ) ) 
               tmpCont => tmpCont%NextRdCont
            ENDDO

            ! Make the last next pointer point to the new container.
            ! The new container will now be the last container in the list!
            tmpCont%NextRdCont => newCont

         ! In all other cases, place new container at beginnig of list
         ELSE
            newCont%NextRdCont => List
            List               => newCont
         ENDIF
      ENDIF

      ! Cleanup
      NULLIFY ( tmpCont )

      END SUBROUTINE RdCont_Add
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: RdCont_Init
!
! !DESCRIPTION: Subroutine RdCont_Init initializes the given ReadList
! container 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE RdCont_Init ( ThisRdCont )
!
! !USES:
!
      USE HCO_EMISLL_MOD, ONLY  : SclMax
!
! !ARGUMENTS:
!
      TYPE(RdCont),   POINTER     :: ThisRdCont
!
! !REVISION HISTORY:
!  11 Apr 2012 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC

      !======================================================================
      ! RdCont_Init begins here!
      !======================================================================

      ! Allocate container
      ALLOCATE ( ThisRdCont )

      ! Nullify pointers
      NULLIFY ( ThisRdCont%UseScalIDs )
      NULLIFY ( ThisRdCont%Array      )
      NULLIFY ( ThisRdCont%NextRdCont )

      END SUBROUTINE RdCont_Init
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: ReadList_Init
!
! !DESCRIPTION: Subroutine ReadList_Init initializes the ReadList. 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE ReadList_Init
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
      ! lists are associated (nullified). 
      ALLOCATE ( ReadList )

      ALLOCATE ( ReadList%Once )
      NULLIFY ( ReadList%Once  )

      ALLOCATE ( ReadList%Year )
      NULLIFY ( ReadList%Year  )

      ALLOCATE ( ReadList%Month )
      NULLIFY ( ReadList%Month )

      ALLOCATE ( ReadList%Day )
      NULLIFY ( ReadList%Day   )

      ALLOCATE ( ReadList%Hour )
      NULLIFY ( ReadList%Hour  )

      ALLOCATE ( ReadList%Always )
      NULLIFY ( ReadList%Always  )

      END SUBROUTINE ReadList_Init
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: ReadList_Print
!
! !DESCRIPTION: Subroutine ReadList_Print displays the content of
! ReadList.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE ReadList_Print
!
! !REVISION HISTORY:
!  20 Apr 2013 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC

      ! ================================================================
      ! ReadList_Print begins here
      ! ================================================================
      
      ! Print content of all lists
      IF ( ASSOCIATED(ReadList) ) THEN 

         write(*,*) 'Content of one-time list:'
         CALL List_Print ( ReadList%Once )

         write(*,*) 'Content of year-list:'
         CALL List_Print ( ReadList%Year )

         write(*,*) 'Content of month-list:'
         CALL List_Print ( ReadList%Month )

         write(*,*) 'Content of day-list:'
         CALL List_Print ( ReadList%Day )

         write(*,*) 'Content of hour-list:'
         CALL List_Print ( ReadList%Hour )

         write(*,*) 'Content of always-to-read list:'
         CALL List_Print ( ReadList%Always )

      ELSE
         write(*,*) 'ReadList not defined yet!!'
      ENDIF

      END SUBROUTINE ReadList_Print
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: List_Print
!
! !DESCRIPTION: Subroutine List_Print displays the content of the given
! list.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE List_Print ( List )
!
! !INPUT ARGUMENTS:
!
      TYPE(RdCont), POINTER   :: List
!
! !REVISION HISTORY:
!  20 Apr 2013 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !ARGUMENTS:
!
      TYPE(RdCont), POINTER   :: thisRdCont

      ! ================================================================
      ! List_Print begins here
      ! ================================================================
      
      ! Init
      NULLIFY ( thisRdCont ) 

      ! Point to first element
      thisRdCont => List
   
      ! Loop over all containers
      DO WHILE ( ASSOCIATED(thisRdCont) ) 
         
         ! Print statements
         write(*,*) 'Container name: ', TRIM(thisRdCont%cName)
         write(*,*) 'Data type     : ', thisRdCont%DataType
         write(*,*) 'Source file   : ', TRIM(thisRdCont%ncFile)
         write(*,*) 'Tracer ID     : ', thisRdCont%TrcID
         write(*,*) 'Scal ID       : ', thisRdCont%ScalID
         write(*,*) 'Cat           : ', thisRdCont%Cat
         write(*,*) 'Hier          : ', thisRdCont%Hier
         write(*,*) 'Oper          : ', thisRdCont%Oper
         write(*,*) 'Add           : ', thisRdCont%Add
         write(*,*) 'ncRead        : ', thisRdCont%ncRead
         write(*,*) 'ESMF_Dim      : ', thisRdCont%ESMF_Dim
         write(*,*) 'ESMF_Unit     : ', thisRdCont%ESMF_Unit
 
         ! Advance
         thisRdCont => thisRdCont%nextRdCont

      ENDDO

      END SUBROUTINE List_Print
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: ReadList_Cleanup
!
! !DESCRIPTION: Subroutine ReadList_Cleanup removes all content of ReadList. 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE ReadList_Cleanup
!
! !REVISION HISTORY:
!  20 Apr 2013 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC

      ! ================================================================
      ! ReadList_Cleanup begins here
      ! ================================================================
      
      IF ( ASSOCIATED(ReadList) ) THEN 

         ! Remove all sublists in ReadList
         CALL List_Cleanup ( ReadList%Once   )
         CALL List_Cleanup ( ReadList%Year   )
         CALL List_Cleanup ( ReadList%Month  )
         CALL List_Cleanup ( ReadList%Day    )
         CALL List_Cleanup ( ReadList%Hour   )
         CALL List_Cleanup ( ReadList%Always )

         ! Remove ReadList 
         DEALLOCATE ( ReadList )
      ENDIF

      END SUBROUTINE ReadList_Cleanup
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: List_Cleanup
!
! !DESCRIPTION: Subroutine List_Cleanup deallocates all content in the list.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE List_Cleanup ( List )
!
! !ARGUMENTS
! 
      TYPE(RdCont),    POINTER  :: List
!
! !REVISION HISTORY:
!  20 Apr 2013 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
      TYPE(RdCont), POINTER    :: thisRdCont => NULL()
      TYPE(RdCont), POINTER    :: nextRdCont => NULL()

      ! ================================================================
      ! LIST_CLEANUP begins here
      ! ================================================================

      ! Init
      thisRdCont => List

      ! Loop over all containers and deallocate them
      DO WHILE ( ASSOCIATED ( thisRdCont ) )

         ! Point to next container before removing the current one 
         nextRdCont => thisRdCont%NextRdCont 

         ! Now remove current container
         CALL RdCont_Cleanup ( thisRdCont )

         ! Advance to next container
         thisRdCont => nextRdCont

      ENDDO

      ! Cleanup
      IF ( ASSOCIATED ( nextRdCont ) ) DEALLOCATE ( nextRdCont )
      IF ( ASSOCIATED ( thisRdCont ) ) DEALLOCATE ( thisRdCont )

      END SUBROUTINE List_Cleanup
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: RdCont_Cleanup
!
! !DESCRIPTION: Subroutine RdCont_Cleanup cleans up the given pool
! container 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE RdCont_Cleanup ( ThisRdCont )
!
! !ARGUMENTS:
!
      TYPE(RdCont), POINTER  :: ThisRdCont
!
! !REVISION HISTORY:
!  11 Apr 2012 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC

      !======================================================================
      ! RdCont_Cleanup begins here!
      !======================================================================

      IF ( ASSOCIATED( ThisRdCont ) ) THEN

         ! First remove connection to next field 
         ThisRdCont%NextRdCont => NULL()

         ! Deallocate all defined pointers
         IF ( ASSOCIATED( ThisRdCont%UseScalIDs ) ) THEN
            DEALLOCATE( ThisRdCont%UseScalIDs )
         ENDIF

         IF ( ASSOCIATED( ThisRdCont%Array ) ) THEN
            DEALLOCATE( ThisRdCont%Array )
         ENDIF

         ! Deallocate the container
         DEALLOCATE( ThisRdCont )

      ENDIF

      END SUBROUTINE RdCont_Cleanup
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: Write_EmisLL 
!
! !DESCRIPTION: Subroutine Write\_EmisLL writes updated data from the
! reading list to the emissions linked list. 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE Write_EmisLL ( am_I_Root, HcoState, EmisLL, Clock, RC )
!
! !USES:
!
      USE HCO_TYPE_MOD,     ONLY : HCO_State
      USE HCO_TIME_MOD,     ONLY : HcoClock
      USE HCO_EMISLL_MOD,   ONLY : EmisCont
!
! !ARGUMENTS:
!
      LOGICAL,         INTENT(IN   )  :: am_I_Root  ! root CPU?
      TYPE(HCO_State), POINTER        :: HcoState   ! Hemco state object
      TYPE(EmisCont),  POINTER        :: EmisLL     ! Emissions linked list
      TYPE(HcoClock),  POINTER        :: Clock      ! Simulation clock 
      INTEGER,         INTENT(INOUT)  :: RC         ! Return Code 
!
!
! !REVISION HISTORY:
!  09 Sep 2013 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!

      ! ================================================================
      ! Write_EmisLL begins here
      ! ================================================================

      ! Always add/update content from always-to-read list
      CALL ReadList2EmisLL ( am_I_Root,       HcoState,   &
                             ReadList%Always, EmisLL,    RC )
      IF ( RC /= HCO_SUCCESS ) RETURN

      ! Add content from one-time list on the first call 
      IF ( Clock%PrevYear < 0 ) THEN
         CALL ReadList2EmisLL ( am_I_Root,     HcoState,   &
                                ReadList%Once, EmisLL,    RC )
         IF ( RC /= HCO_SUCCESS ) RETURN
      ENDIF

      ! Add content from year-list if it's a new year
      IF ( Clock%PrevYear /= Clock%ThisYear ) THEN
         CALL ReadList2EmisLL ( am_I_Root,     HcoState,   &
                                ReadList%Year, EmisLL,    RC )
         IF ( RC /= HCO_SUCCESS ) RETURN
      ENDIF

      ! Add content from month-list if it's a new month
      IF ( Clock%PrevMonth /= Clock%ThisMonth ) THEN
         CALL ReadList2EmisLL ( am_I_Root,      HcoState,   &
                                ReadList%Month, EmisLL,    RC )
         IF ( RC /= HCO_SUCCESS ) RETURN
      ENDIF

      ! Add content from day-list if it's a new day
      IF ( Clock%PrevDay /= Clock%ThisDay ) THEN
         CALL ReadList2EmisLL ( am_I_Root,    HcoState,   &
                                ReadList%Day, EmisLL,    RC )
         IF ( RC /= HCO_SUCCESS ) RETURN
      ENDIF

      ! Add content from hour-list if it's a new hour
      IF ( Clock%PrevHour /= Clock%ThisHour ) THEN
         CALL ReadList2EmisLL ( am_I_Root,     HcoState,   &
                                ReadList%Hour, EmisLL,    RC )
         IF ( RC /= HCO_SUCCESS ) RETURN
      ENDIF

      ! Assign scale factor container IDs to base containers (only on
      ! first call)
      IF ( Clock%PrevYear < 0 ) THEN
         CALL AssignScalIDs( am_I_Root, ReadList%Once,   EmisLL, RC )
         IF ( RC /= HCO_SUCCESS ) RETURN
         CALL AssignScalIDs( am_I_Root, ReadList%Always, EmisLL, RC )
         IF ( RC /= HCO_SUCCESS ) RETURN
         CALL AssignScalIDs( am_I_Root, ReadList%Year,   EmisLL, RC )
         IF ( RC /= HCO_SUCCESS ) RETURN
         CALL AssignScalIDs( am_I_Root, ReadList%Month,  EmisLL, RC )
         IF ( RC /= HCO_SUCCESS ) RETURN
         CALL AssignScalIDs( am_I_Root, ReadList%Day,    EmisLL, RC )
         IF ( RC /= HCO_SUCCESS ) RETURN
         CALL AssignScalIDs( am_I_Root, ReadList%Hour,   EmisLL, RC )
         IF ( RC /= HCO_SUCCESS ) RETURN
      ENDIF

      ! Leave w/ success
      RC = HCO_SUCCESS

      END SUBROUTINE Write_EmisLL
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: ReadList2EmisLL
!
! !DESCRIPTION: Subroutine ReadList2EmisLL translates all data arrays in
! the specified ReadList into the emissions list EmisLL.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE ReadList2EmisLL ( am_I_Root, HcoState,   &
                                   TempList,  EmisLL,    RC )
!
! !USES:
!
      USE HCO_EMISLL_MOD,         ONLY : EmisCont
!
! !ARGUMENTS:
!
      LOGICAL,         INTENT(IN   )  :: am_I_Root
      TYPE(HCO_State), POINTER        :: HcoState
      TYPE(RdCont),    POINTER        :: TempList
      TYPE(EmisCont),  POINTER        :: EmisLL
      INTEGER,         INTENT(INOUT)  :: RC
!
! !REVISION HISTORY:
!  20 Apr 2013 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER                     :: cID, iScalID
      INTEGER                     :: I
      TYPE(RdCont), POINTER       :: tmpRdCont => NULL()
      CHARACTER(LEN=31)           :: ScalName

      ! ================================================================
      ! ReadList2EmisLL begins here
      ! ================================================================

      ! Point to header of emissions pool list
      tmpRdCont => TempList

      ! Loop over all containers
      DO WHILE ( ASSOCIATED( tmpRdCont ) )

         ! only if array is defined...
         IF ( ASSOCIATED( tmpRdCont%Array ) ) THEN

            ! Add/update content in emissions linked list (EmisLL).
            ! The array in the ReadList is copied to EmisLL and removed
            ! afterwards from the ReadList.
            CALL AddCont2EmisLL( am_I_Root, HcoState, tmpRdCont, &
                                 EmisLL,    cID,       RC          )
            IF ( RC /= HCO_SUCCESS ) RETURN

            ! Pass container ID
            tmpRdCont%cID = cID
         ENDIF

         ! Point to next scale factor container
         tmpRdCont => tmpRdCont%NextRdCont

      ENDDO

      ! Leave w/ success
      RC = HCO_SUCCESS

      END SUBROUTINE ReadList2EmisLL
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: AssignScalIDs
!
! !DESCRIPTION: Subroutine AssignScalIDs assigns the container IDs of
! the scale factors to be used to base emissions. 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE AssignScalIDs( am_I_Root, TempList, EmisLL, RC )
!
! !USES:
!
      USE HCO_EMISLL_MOD,         ONLY  : EmisCont
      USE HCO_EMISLL_MOD,         ONLY  : Add_SclID
      USE HCO_EMISLL_MOD,         ONLY  : Show_EmisLL 
!
! !ARGUMENTS:
!
      LOGICAL,         INTENT(IN)     :: am_I_Root
      TYPE(RdCont),    POINTER        :: TempList
      TYPE(EmisCont),  POINTER        :: EmisLL
      INTEGER,         INTENT(INOUT)  :: RC
!
! !REVISION HISTORY:
!  20 Apr 2013 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER                     :: cID, iScalID
      INTEGER                     :: I
      TYPE(RdCont), POINTER       :: tmpRdCont => NULL()

      ! ================================================================
      ! ReadList2EmisLL begins here
      ! ================================================================

      ! Point to header of emissions pool list
      tmpRdCont => TempList

      ! Loop over all containers
      DO WHILE ( ASSOCIATED( tmpRdCont ) )

         ! Eventually set scale factors
         IF ( ASSOCIATED( tmpRdCont%UseScalIDs ) ) THEN
               
            ! Get current container ID:
            cID = tmpRdCont%cID

            ! Loop over the maximum number of scale factors 
            DO I = 1, SIZE( tmpRdCont%UseScalIDs ) 

               ! Get current scale field number
               iScalID = tmpRdCont%UseScalIDs(I)

               ! Leave loop if invalid scale field number
               IF ( iScalID <= 0 ) EXIT

               ! Assign scale factor to base emissions
               CALL Add_SclID( cID, iScalID, EmisLL, RC )
               IF ( RC /= HCO_SUCCESS ) RETURN
            ENDDO
         ENDIF

         ! Point to next scale factor container
         tmpRdCont => tmpRdCont%NextRdCont

      ENDDO

      ! Cleanup
      tmpRdCont => NULL()

      ! Leave w/ success
      RC = HCO_SUCCESS

      END SUBROUTINE AssignScalIDs
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: AddCont2EmisLL
!
! !DESCRIPTION: Subroutine Cont2EmisLL adds the specified container
! to the emissions linked list.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE AddCont2EmisLL ( am_I_Root,  HcoState, &
                                  ThisRdCont, EmisLL,    & 
                                  cID,        RC          )
!
! !USES:
!
      USE HCO_EMISLL_MOD,   ONLY : EmisCont
      USE HCO_EMISLL_MOD,   ONLY : Add_EmisCont, Find_EmisCont
!
! !ARGUMENTS:
!
      LOGICAL,          INTENT(IN)     :: am_I_Root
      TYPE(HCO_State),  POINTER        :: HcoState
      TYPE(RdCont),     POINTER        :: ThisRdCont ! ReadList container
      TYPE(EmisCont),   POINTER        :: EmisLL     ! Emissions linked list
      INTEGER,          INTENT(OUT)    :: cID        ! Container ID in LL
      INTEGER,          INTENT(INOUT)  :: RC         ! Success 
!
! !REVISION HISTORY:
!  28 Mar 2013 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      TYPE(EmisCont),  POINTER  :: ThisCont => NULL()
      INTEGER                   :: I, J, L, T
      LOGICAL                   :: FOUND
      CHARACTER(LEN=255)        :: MSG, LOC
      CHARACTER(LEN= 31)        :: TempRes
      REAL*8, POINTER           :: TMPARR(:,:,:,:) => NULL()
      REAL*8, POINTER           :: NEWARR(:,:,:,:) => NULL()
 
      ! ================================================================
      ! AddCont2EmisLL begins here
      ! ================================================================

      ! Enter
      LOC = 'AddCont2EmisLL'

      ! ----------------------------------------------------------------
      ! Create argument list for subroutine CREATE_FIELD
      ! ----------------------------------------------------------------

      ! NEWARR is the new array (to be added to the emissions list)
      NEWARR => ThisRdCont%Array

      ! Get array size
      I = SIZE(NEWARR,1)
      J = SIZE(NEWARR,2)
      L = SIZE(NEWARR,3)
      T = SIZE(NEWARR,4)

      ! ---------------------------------------
      ! Auto detection of temporal resolution
      ! Eventually extend new array so that it
      ! matches required time dimensions.
      ! ---------------------------------------

      ! No internal time variation
      IF ( T == 1 ) THEN
         TempRes = 'CONSTANT'

      ! Hourly varying data
      ELSEIF ( T == 24 ) THEN 

         ! Check if data is gridded or not, set TempRes attribute
         ! accordingly
         IF ( I == 1 .AND. J == 1 ) THEN
            TempRes = 'HOURLY'
         ELSE
            TempRes = 'HOURLY_GRID'
         ENDIF

      ! Monthly varying data
      ! This must be a scalar. Gridded monthly data should be read only
      ! slice by slice!
      ELSEIF ( T == 12 ) THEN 

         IF ( I == 1 .AND. J == 1 ) THEN
            TempRes = 'MONTHLY'
         ELSE
            MSG = 'Monthly data not gridded:' // TRIM(ThisRdCont%cName)
            CALL HCO_ERROR ( MSG, LOC, RC ); RETURN
         ENDIF

      ! Data varying with day of the week 
      ELSEIF ( T == 7 ) THEN 
   
         ! Check if data is gridded or not, set TempRes attribute
         ! accordingly
         IF ( I == 1 .AND. J == 1 ) THEN
            TempRes = 'WEEKDY'
         ELSE
            TempRes = 'WEEKDY_GRID'
         ENDIF

      ELSE
         MSG = 'Invalid time dimension for field ' // &
               TRIM(ThisRdCont%cName)
         CALL HCO_ERROR ( MSG, LOC, RC ); RETURN
      ENDIF

      ! ---------------------------------------
      ! Add fields if required 
      ! ---------------------------------------

      ! Check if field already exists
      CALL Find_EmisCont( TRIM(ThisRdCont%cName), EmisLL, & 
                          FOUND,                  ThisCont )

      ! If this field already exist and the add flag is on, then add 
      ! the emissions to the existing ones. Check all other parameter 
      ! for consistency beforehand
      IF ( FOUND .AND. ThisRdCont%Add ) THEN

         ! Check extension number 
         IF ( ThisRdCont%ExtNr /= ThisCont%ExtNr ) THEN
            MSG = 'Wrong extension number: ' // TRIM(ThisRdCont%cName)
            CALL HCO_ERROR ( MSG, LOC, RC ); RETURN
         ENDIF

         ! Check data type 
         IF ( ThisRdCont%DataType /= ThisCont%DataType ) THEN
            MSG = 'Wrong data type: ' // TRIM(ThisRdCont%cName)
            CALL HCO_ERROR ( MSG, LOC, RC ); RETURN
         ENDIF

         ! Check tracer ID
         IF ( ThisRdCont%TrcID /= ThisCont%TrcID ) THEN
            MSG = 'Wrong tracer ID: ' // TRIM(ThisRdCont%cName)
            CALL HCO_ERROR ( MSG, LOC, RC ); RETURN
         ENDIF

         ! Check for 3rd and 4th dimension
         IF ( L /= SIZE(ThisCont%Array,3) ) THEN
            MSG = 'Wrong 3rd dim: ' // TRIM(ThisRdCont%cName)
            CALL HCO_ERROR ( MSG, LOC, RC ); RETURN
         ENDIF
         IF ( T /= SIZE(ThisCont%Array,4) ) THEN
            MSG = 'Wrong 4th dim: ' // TRIM(ThisRdCont%cName)
            CALL HCO_ERROR ( MSG, LOC, RC ); RETURN
         ENDIF

         ! Check operator 
         IF ( ThisRdCont%Oper /= ThisCont%Oper ) THEN
            MSG = 'Wrong operator: ' // TRIM(ThisRdCont%cName)
            CALL HCO_ERROR ( MSG, LOC, RC ); RETURN
         ENDIF

         ! Check category 
         IF ( ThisRdCont%Cat /= ThisCont%Cat ) THEN
            MSG = 'Wrong category: ' // TRIM(ThisRdCont%cName)
            CALL HCO_ERROR ( MSG, LOC, RC ); RETURN
         ENDIF

         ! Check hierarchy 
         IF ( ThisRdCont%Hier /= ThisCont%Hier ) THEN
            MSG = 'Wrong hierarchy: ' // TRIM(ThisRdCont%cName)
            CALL HCO_ERROR ( MSG, LOC, RC ); RETURN
         ENDIF

         ! If all checks are successfully passed, pass emissions
         ! to TMPARR
         ALLOCATE(TMPARR(I,J,L,T))
         TMPARR(:,:,:,:) = ThisCont%Array(:,:,:,:) + NEWARR(:,:,:,:) 

      ! If we don't add fields, just make TMPARR point to NEWARR.
      ELSE
         TMPARR => NEWARR
      ENDIF

      ! ----------------------------------------------------------------
      ! Add/ update emissions linked list container. 
      ! ----------------------------------------------------------------
      CALL Add_EmisCont ( aIR       = am_I_Root,           & 
                          HcoState = HcoState,           &
                          DataType  = ThisRdCont%DataType, &
                          ExtNr     = ThisRdCont%ExtNr,    &
                          cName     = ThisRdCont%cName,    &
                          Array     = TMPARR,              &
                          scalID    = ThisRdCont%scalID,   &
                          TrcID     = ThisRdCont%TrcID,    &
                          Oper      = ThisRdCont%Oper,     &
                          Cat       = ThisRdCont%Cat,      &
                          Hier      = ThisRdCont%Hier,     &
                          TempRes   = TempRes,             &
                          EmisLL    = EmisLL,              &
                          cID       = cID,                 &
                          RC        = RC )
      IF ( RC /= HCO_SUCCESS ) RETURN

      ! ----------------------------------------------------------------
      ! Cleanup and quit
      ! ----------------------------------------------------------------

      ! Nullify/deallocate TMPARR
      IF ( ASSOCIATED(TMPARR, TARGET=NEWARR) ) THEN
         TMPARR => NULL()
      ELSE 
         IF ( ASSOCIATED(TMPARR) ) DEALLOCATE ( TMPARR ) 
      ENDIF

      ! Nullify/deallocate NEWARR
      IF ( ASSOCIATED(NEWARR, TARGET=ThisRdCont%Array) ) THEN
         NEWARR => NULL()
      ELSE
         IF ( ASSOCIATED(NEWARR) ) DEALLOCATE ( NEWARR )
      ENDIF 

      ! We don't need the array in the ReadList anymore, so remove it!
      IF (ASSOCIATED(ThisRdCont%Array)) DEALLOCATE ( ThisRdCont%Array)

      ! Return w/ success
      RC = HCO_SUCCESS

      END SUBROUTINE AddCont2EmisLL
!EOC
      END MODULE HCO_READLIST_MOD
!EOM
