!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: hco_filedata_mod.F90
!
! !DESCRIPTION: Module HCO\_Filedata\_Mod contains routines and
! variables to handle the HEMCO file data object FileData.
! FileData holds all information of source file data, such as file
! name, update frequency, temporal resolution, the data array itself,
! etc. These values are specified in the HEMCO configuration file. Many
! of these attributes are primarily used for reading/updating the data
! from file using the HEMCO generic file reading routines. Within an
! ESMF environment, these attributes may be obsolete.
!\\
!\\
! FileData consists of the following elements:
!
! \begin{itemize}
! \item ncFile: path and filename to the source file, as specified in
!       the configuration file.
! \item ncPara: file parameter (variable) of interest, as specified in
!       the configuration file.
! \item ncYrs: range of years in the source file, as specified in the
!       configuration file through the timestamp attribute.
! \item ncMts: range of months in the source file, as specified in the
!       configuration file through the timestamp attribute.
! \item ncDys: range of days in the source file, as specified in the
!       configuration file through the timestamp attribute.
! \item ncHrs: range of hours in the source file, as specified in the
!       configuration file through the timestamp attribute.
! \item CycleFlag: determines how to deal with time stamps that do not
!       correspond to one of the source file time slices. If set to 1,
!       the closest available time slice (in the past) is used (or the
!       first available time slice if model time is before first
!       available time slice). If set to 2, the file data is ignored if
!       model time is outside of the source file range. If CycleFlag is
!       set to 3, an error is returned if none of the file time slices
!       matches the model time.
! \item MustFind: if yes, the code returns with an error if no field
!       can be found for the given simulation time (and according to
!       the cycle flag and time attribute settings). Only of relevance
!       for cycle flags range and exact.
! \item UpdtFlag: determines the update frequency of the data. This is
!       currently only used to distinguish containers that are updated
!       on every time step (always) or according to the frequency
!       provided in the HEMCO configuration file via attribute 'srcTime'.
! \item ncRead: logical denoting whether or not we need to read this
!       data container. ncRead is set to false for containers whose
!       data is directly specified in the configuration file. For
!       internal use only.
! \item OrigUnit: original unit of data.
! \item ArbDimName: name of additional (arbitrary) file dimension.
! \item ArbDimVal : desired value of arbitrary dimension.
! \item IsConc: Set to true if data is concentration. Concentration
!       data will be added to the concentration array instead of the
!       emission array.
! \item IsLocTime: Set to true if data is in local time. Defaults to
!       false and becomes only true if data is scalar (e.g. uniform
!       diurnal scale factors), country-specific data (read from
!       ASCII), or weekdaily data.
! \item V3: vector of 3D fields. For 3D-data, this vector will hold
!       the 3D arrays of all time slices kept in memory (e.g. 24
!       elements for hourly data).
! \item V2: vector of 2D fields. For 2D-data, this vector will hold
!       the 2D arrays of all time slices kept in memory (e.g. 24
!       elements for hourly data).
! \item tIDx: derived type used for proper indexing of the time slices
!       in memory. Internal use only.
! \item Cover: data coverage on this CPU: 0=no overlap; 1=full overlap;
!       -1=partial overlap. As determined from the mask regions
!       specified in the configuration file.
! \item SpaceDim: spatial dimension of data array: 1 = spatially uniform
!       (x=y=z=1); 2 = 2D data (x,y); 3 = 3D data (x,y,z).
! \item Levels: handling of vertical levels (3D data only). For internal
!       use only.
! \item nt: time dimension. length of vector V3 or V2. For internal use
!       only.
! \item DeltaT: time interval between time slices. For internal use only.
!       ID i (e.g. cIDList(3) points to data-container w/ cID = 3).
! \item DoShare: will be set to True if this file data object is shared
!       by multiple data containers. For internal use only.
! \item IsInList: will be set to True if this file data object is part
!       of the emissions list EmisList. For internal use only.
! \item IsTouched: will be set to True as soon as the container becomes
!       touched for the first time. For internal use only.
! \end{itemize}
!
! !INTERFACE:
!
MODULE HCO_FileData_Mod
!
! !USES:
!
  USE HCO_TYPES_MOD
  USE HCO_ERROR_MOD
  USE HCO_ARR_MOD

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC  :: FileData_Init
  PUBLIC  :: FileData_Cleanup
  PUBLIC  :: FileData_ArrCheck
  PUBLIC  :: FileData_ArrIsDefined
  PUBLIC  :: FileData_ArrIsTouched
  PUBLIC  :: FileData_ArrInit
!
! !PRIVATE MEMBER FUNCTIONS:
!
  PRIVATE :: FileData_ArrCheck2D
  PRIVATE :: FileData_ArrCheck3D
!
! !REVISION HISTORY:
!  19 Dec 2013 - C. Keller   - Initialization
!  01 Jul 2014 - R. Yantosca - Cosmetic changes in ProTeX headers
!  01 Jul 2014 - R. Yantosca - Now use F90 free-format indentation
!  21 Aug 2014 - C. Keller   - Added concentration
!  23 Dec 2014 - C. Keller   - Added argument IsInList
!  06 Oct 2015 - C. Keller   - Added argument MustFind
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !PRIVATE TYPES:
!
!
! !INTERFACES:
!
  INTERFACE FileData_ArrCheck
     MODULE PROCEDURE FileData_ArrCheck2D
     MODULE PROCEDURE FileData_ArrCheck3D
  END INTERFACE FileData_ArrCheck

  INTERFACE FileData_ArrInit
     MODULE PROCEDURE FileData_ArrInit2D
     MODULE PROCEDURE FileData_ArrInit3D
  END INTERFACE FileData_ArrInit

CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: FileData_Init
!
! !DESCRIPTION: Subroutine FileData\_Init initializes a new (blank) file
! data object.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE FileData_Init( FileDta )
!
! !INPUT PARAMETERS:
!
    TYPE(FileData), POINTER    :: FileDta
!
! !REVISION HISTORY:
!  19 Dec 2013 - C. Keller   - Initialization
!  21 Aug 2014 - C. Keller   - Added concentration
!  26 Oct 2016 - R. Yantosca - Don't nullify local ptrs in declaration stmts
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    TYPE(FileData), POINTER  :: NewFDta

    !======================================================================
    ! FileData_Init begins here!
    !======================================================================

    ! Allocate the new file data object
    ALLOCATE( NewFDta )

    ! Nullify all pointers and initialize variables
    NewFDta%V3          => NULL()
    NewFDta%V2          => NULL()
    NewFDta%tIDx        => NULL()
    NewFDta%ncFile       = ''
    NewFDta%ncPara       = ''
    NewFDta%ncYrs(:)     = -999
    NewFDta%ncMts(:)     = -999
    NewFDta%ncDys(:)     = -999
    NewFDta%ncHrs(:)     = -999
    NewFDta%tShift(:)    = 0
    NewFDta%CycleFlag    = HCO_CFLAG_CYCLE
    NewFDta%UpdtFlag     = HCO_UFLAG_FROMFILE
    NewFDta%MustFind     = .FALSE.
    NewFDta%ncRead       = .TRUE.
    NewFDta%Cover        = -999
    NewFDta%DeltaT       = 0
    NewFDta%nt           = 0
    NewFDta%SpaceDim     = -1
    NewFDta%Levels       = 0
    NewFDta%EmisL1       = 1.0_hp
    NewFDta%EmisL2       = 1.0_hp
    NewFDta%EmisL1Unit   = HCO_EMISL_LEV
    NewFDta%EmisL2Unit   = HCO_EMISL_LEV
    NewFDta%OrigUnit     = ''
    NewFDta%ArbDimName   = 'none'
    NewFDta%ArbDimVal    = ''
    NewFDta%IsLocTime    = .FALSE.
    NewFDta%IsConc       = .FALSE.
    NewFDta%DoShare      = .FALSE.
    NewFDta%IsInList     = .FALSE.
    NewFDta%IsTouched    = .FALSE.

    ! Return
    FileDta => NewFDta

  END SUBROUTINE FileData_Init
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: FileData_Cleanup
!
! !DESCRIPTION: Subroutine FileData\_Cleanup cleans up the file data object
! FileDta. If DeepClean is set to False, only the data arrays will be removed.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE FileData_Cleanup( FileDta, DeepClean )
!
! !INPUT PARAMETERS:
!
    TYPE(FileData), POINTER    :: FileDta
    LOGICAL,        INTENT(IN) :: DeepClean
!
! !REVISION HISTORY:
!  19 Dec 2013 - C. Keller: Initialization
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: I

    !======================================================================
    ! FileData_Cleanup begins here!
    !======================================================================

    ! Only if associated
    IF ( ASSOCIATED( FileDta ) ) THEN

       ! Deallocate data arrays
       CALL HCO_ArrCleanup( FileDta%V3, DeepClean )
       CALL HCO_ArrCleanup( FileDta%V2, DeepClean )
       FileDta%nt = 0

       IF ( DeepClean ) THEN
          FileDta%tIDx => NULL()
          DEALLOCATE ( FileDta )
       ENDIF
    ENDIF

  END SUBROUTINE FileData_Cleanup
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: FileData_ArrCheck2D
!
! !DESCRIPTION: Subroutine FileData\_ArrCheck2D allocates the 2D data array
! vector of the given file data object if it is not yet allocated. If already
! allocated, it compares the array dimensions against the passed dimensions.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE FileData_ArrCheck2D( HcoConfig, FileDta, nx, ny, nt, RC )
!
! !INPUT PARAMETERS:
!
    TYPE(ConfigObj),POINTER       :: HcoConfig ! HEMCO config object
    TYPE(FileData), POINTER       :: FileDta   ! file data object
    INTEGER,        INTENT(IN)    :: nx        ! x-dim
    INTEGER,        INTENT(IN)    :: ny        ! y-dim
    INTEGER,        INTENT(IN)    :: nt        ! time dim => vector length
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(INOUT) :: RC        ! Return code
!
! !REVISION HISTORY:
!  20 Apr 2013 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER            :: I, AS
    CHARACTER(LEN=255) :: MSG

    ! ================================================================
    ! FileData_ArrCheck2D begins here
    ! ================================================================

    ! Assume success until otherwise
    RC = HCO_SUCCESS

    ! Compare dimensions if array already allocated
    IF ( ASSOCIATED(FileDta%V2) ) THEN
       IF ( (               FileDta%nt  /= nt ) .OR. &
            ( SIZE(FileDta%V2(1)%Val,1) /= nx ) .OR. &
            ( SIZE(FileDta%V2(1)%Val,2) /= ny )       ) THEN
          MSG = 'Wrong dimensions: ' // TRIM(FileDta%ncFile)
          CALL HCO_ERROR ( HcoConfig%Err, MSG, RC )
       ENDIF
       RETURN
    ENDIF

    ! If not associated yet:
    ! Initialize vector and corresponding arrays.
    CALL FileData_ArrInit ( FileDta, nt, nx, ny, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

  END SUBROUTINE FileData_ArrCheck2D
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: FileData_ArrCheck3D
!
! !DESCRIPTION: Subroutine FileData\_ArrCheck3D allocates the 3D data array
! vector of the given file data object if it is not yet allocated. If already
! allocated, it compares the array dimensions against the passed dimensions.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE FileData_ArrCheck3D( HcoConfig, FileDta, nx, ny, nz, nt, RC )
!
! !USES:
!
!
! !INPUT PARAMETERS:
!
    TYPE(ConfigObj),POINTER       :: HcoConfig ! HEMCO config object
    TYPE(FileData), POINTER       :: FileDta  ! Container
    INTEGER,        INTENT(IN)    :: nx        ! x-dim
    INTEGER,        INTENT(IN)    :: ny        ! y-dim
    INTEGER,        INTENT(IN)    :: nz        ! z-dim
    INTEGER,        INTENT(IN)    :: nt        ! Time dim => vector length
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(INOUT) :: RC        ! Return code
!
! !REVISION HISTORY:
!  20 Apr 2013 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER            :: I, AS
    CHARACTER(LEN=255) :: MSG

    ! ================================================================
    ! FileData_ArrCheck3D begins here
    ! ================================================================

    ! Assume success until otherwise
    RC = HCO_SUCCESS

    ! Compare dimensions if array already allocated
    IF ( Associated(FileDta%V3) ) THEN
       IF ( (                FileDta%nt /= nt ) .OR. &
            ( SIZE(FileDta%V3(1)%Val,1) /= nx ) .OR. &
            ( SIZE(FileDta%V3(1)%Val,2) /= ny ) .OR. &
            ( SIZE(FileDta%V3(1)%Val,3) /= nz )       ) THEN
          MSG = 'Wrong dimensions: ' // TRIM(FileDta%ncFile)
          CALL HCO_ERROR ( HcoConfig%Err, MSG, RC )
       ENDIF
       RETURN
    ENDIF

    ! If not associated yet:
    ! Initialize vector and corresponding arrays.
    CALL FileData_ArrInit( FileDta, nt, nx, ny, nz, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

  END SUBROUTINE FileData_ArrCheck3D
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: FileData_ArrIsDefined
!
! !DESCRIPTION: Function FileData\_ArrIsDefined returns true if the data
! array of the given file data object is defined.
!\\
!\\
! !INTERFACE:
!
  FUNCTION FileData_ArrIsDefined( FileDta ) RESULT( IsDefined )
!
! !INPUT PARAMETERS:
!
    TYPE(FileData), POINTER :: FileDta  ! Container
!
! !RETURN VALUE:
!
    LOGICAL                 :: IsDefined
!
! !REVISION HISTORY:
!  20 Apr 2013 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC

    ! ================================================================
    ! FileData_ArrIsDefined begins here
    ! ================================================================

    ! Init
    IsDefined = .FALSE.

    ! Return here if passed FileDta object is not defined
    IF ( .NOT. ASSOCIATED( FileDta ) ) RETURN

    ! nt must be larger than zero!
    IF ( FileDta%nt <= 0 ) Return

    ! 2D array
    IF ( (FileDta%SpaceDim<=2) .AND. ASSOCIATED(FileDta%V2) ) THEN
       IsDefined = .TRUE.
    ENDIF

    ! 3D array
    IF ( (FileDta%SpaceDim==3) .AND. ASSOCIATED(FileDta%V3) ) THEN
       IsDefined = .TRUE.
    ENDIF

  END FUNCTION FileData_ArrIsDefined
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: FileData_ArrIsTouched
!
! !DESCRIPTION: Function FileData\_ArrIsTouched returns true if the data
! array of the given file data object has already been touched, e.g. if
! the data has already been read (or at least attempted to being read).
! This information is mostly important for file data objects that are shared
! by multiple data containers. See ReadList\_Fill in hco\_readlist\_mod.F90
! for more details.
!\\
!\\
! !INTERFACE:
!
  FUNCTION FileData_ArrIsTouched( FileDta ) RESULT( IsTouched )
!
! !INPUT PARAMETERS:
!
    TYPE(FileData), POINTER :: FileDta  ! Container
!
! !RETURN VALUE:
!
    LOGICAL                 :: IsTouched
!
! !REVISION HISTORY:
!  17 Mar 2015 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC

    ! ================================================================
    ! FileData_ArrIsTouched begins here
    ! ================================================================

    ! Init
    IsTouched = .FALSE.

    ! Return touched flag
    IF ( ASSOCIATED( FileDta ) ) THEN
       IsTouched = FileDta%IsTouched
    ENDIF


  END FUNCTION FileData_ArrIsTouched
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: FileData_ArrInit2D
!
! !DESCRIPTION: Subroutine FileData\_ArrInit2D is a wrapper routine
! to initialize 2D data arrays of a file data object. To ensure proper
! functioning of the file data object and related routines, this routine
! should always be used to initialize file data arrays (and NOT HCO\_ArrInit
! directly!).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE FileData_ArrInit2D( FileDta, nt, nx, ny, RC )
!
! !USES:
!
!
! !INPUT PARAMETERS:
!
    TYPE(FileData), POINTER       :: FileDta  ! Container
    INTEGER,        INTENT(IN)    :: nt        ! Time dim => vector length
    INTEGER,        INTENT(IN)    :: nx        ! x-dim
    INTEGER,        INTENT(IN)    :: ny        ! y-dim
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(INOUT) :: RC        ! Return code
!
! !REVISION HISTORY:
!  01 Oct 2014 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! ================================================================
    ! FileData_ArrInit2D begins here
    ! ================================================================

    ! Assume success until otherwise
    RC = HCO_SUCCESS

    ! Initialize vector and corresponding arrays.
    CALL HCO_ArrInit( FileDta%V2, nt, nx, ny, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Update nt
    FileDta%nt = nt

    ! Return w/ success
    RC = HCO_SUCCESS

  END SUBROUTINE FileData_ArrInit2D
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: FileData_ArrInit3D
!
! !DESCRIPTION: Subroutine FileData\_ArrInit3D is a wrapper routine
! to initialize 3D data arrays of a file data object. To ensure proper
! functioning of the file data object and related routines, this routine
! should always be used to initialize file data arrays (and NOT HCO\_ArrInit
! directly!).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE FileData_ArrInit3D( FileDta, nt, nx, ny, nz, RC )
!
! !USES:
!
!
! !INPUT PARAMETERS:
!
    TYPE(FileData), POINTER       :: FileDta  ! Container
    INTEGER,        INTENT(IN)    :: nt        ! Time dim => vector length
    INTEGER,        INTENT(IN)    :: nx        ! x-dim
    INTEGER,        INTENT(IN)    :: ny        ! y-dim
    INTEGER,        INTENT(IN)    :: nz        ! z-dim
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(INOUT) :: RC        ! Return code
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
    ! FileData_ArrInit3D begins here
    ! ================================================================

    ! Assume success until otherwise
    RC = HCO_SUCCESS

    ! Initialize vector and corresponding arrays.
    CALL HCO_ArrInit( FileDta%V3, nt, nx, ny, nz, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Update nt
    FileDta%nt = nt

    ! Return w/ success
    RC = HCO_SUCCESS

  END SUBROUTINE FileData_ArrInit3D
!EOC
END MODULE HCO_FileData_Mod
