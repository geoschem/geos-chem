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
! \item ncRead: logical denoting whether or not we need to read this
!       data container. ncRead is set to false for containers whose
!       data is directly specified in the configuration file. For
!       internal use only. 
! \item OrigUnit: original unit of data.
! \item IsConc: Set to true if data is concentration. Concentration 
!       data will be added to the concentration array instead of the
!       emission array.
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
! \item nt: time dimension. length of vector V3 or V2. For internal use 
!       only.
! \item DeltaT: time interval between time slices. For internal use only.
      ! ID i (e.g. cIDList(3) points to data-container w/ cID = 3). 
! \item DoShare: will be set to True if this file data object is shared
!       by multiple data containers. For internal use only. 
! \end{itemize}
!
! !INTERFACE: 
!
MODULE HCO_FileData_Mod 
!
! !USES:
!
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
  PUBLIC  :: FileData_FileRead
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
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !PRIVATE TYPES:
!
  !-------------------------------------------------------------------------
  ! FileData: Derived type definition for HEMCO filetype object
  !-------------------------------------------------------------------------
  TYPE, PUBLIC :: FileData
     CHARACTER(LEN=255)          :: ncFile    ! file path+name
     CHARACTER(LEN= 31)          :: ncPara    ! file parameter
     INTEGER                     :: ncYrs(2)  ! year range
     INTEGER                     :: ncMts(2)  ! month range
     INTEGER                     :: ncDys(2)  ! day range
     INTEGER                     :: ncHrs(2)  ! hour range
     INTEGER                     :: CycleFlag ! cycle flag
     LOGICAL                     :: ncRead    ! read from source?
     TYPE(Arr3D_HP),     POINTER :: V3(:)     ! vector of 3D fields
     TYPE(Arr2D_HP),     POINTER :: V2(:)     ! vector of 2D fields
     TYPE(TimeIdx),      POINTER :: tIDx      ! for time slice indexing 
     CHARACTER(LEN= 31)          :: OrigUnit  ! original data units 
     LOGICAL                     :: IsConc    ! concentration data?
     INTEGER                     :: Cover     ! data coverage
     INTEGER                     :: SpaceDim  ! space dimension: 1, 2 or 3 
     INTEGER                     :: nt        ! time dimension: length of Arr
     INTEGER                     :: DeltaT    ! temp. resolution of array [h]
     LOGICAL                     :: DoShare   ! shared object?
  END TYPE FileData

  !-------------------------------------------------------------------------
  ! TimeIdx: Derived type definition for the object that points to the 
  ! current time slices of data within a file.  Used by hco_tidx_mod.F90.
  !-------------------------------------------------------------------------
  TYPE, PUBLIC :: TimeIdx
     INTEGER                     :: TypeID
     CHARACTER(LEN=31)           :: TempRes
     LOGICAL                     :: LonDependent
     INTEGER,            POINTER :: CurrIDx(:)
  END TYPE TimeIdx
!
! !INTERFACES:
!
  INTERFACE FileData_ArrCheck
     MODULE PROCEDURE FileData_ArrCheck2D
     MODULE PROCEDURE FileData_ArrCheck3D
  END INTERFACE FileData_ArrCheck

CONTAINS
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group
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
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    TYPE(FileData), POINTER  :: NewFDta => NULL()

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
    NewFDta%CycleFlag    = 1
    NewFDta%ncRead       = .TRUE.
    NewFDta%Cover        = -999 
    NewFDta%DeltaT       = 0
    NewFDta%nt           = 0
    NewFDta%SpaceDim     = -1
    NewFDta%OrigUnit     = ''
    NewFDta%IsConc       = .FALSE.
    NewFDta%DoShare      = .FALSE.

    ! Return
    FileDta => NewFDta

  END SUBROUTINE FileData_Init
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group
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
  SUBROUTINE FileData_ArrCheck2D( FileDta, nx, ny, nt, RC )
!
! !INPUT PARAMETERS:
!
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
          CALL HCO_ERROR ( MSG, RC )
       ENDIF
       RETURN
    ENDIF

    ! If not associated yet:
    ! Initialize vector and corresponding arrays.
    CALL HCO_ArrInit ( FileDta%V2, nt, nx, ny, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Update nt
    FileDta%nt = nt

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
  SUBROUTINE FileData_ArrCheck3D( FileDta, nx, ny, nz, nt, RC )
!
! !USES:
!
!
! !INPUT PARAMETERS:
!
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
          CALL HCO_ERROR ( MSG, RC )
       ENDIF
       RETURN
    ENDIF

    ! If not allocated:
    ! Initialize vector and corresponding arrays.
    CALL HCO_ArrInit( FileDta%V3, nt, nx, ny, nz, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Update nt
    FileDta%nt = nt

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
! !IROUTINE: FileData_FileRead
!
! !DESCRIPTION: Subroutine FileData\_FileRead passes scalar scale factors
! read from the configuration file (instead of the file path+name) to a 
! data array.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE FileData_FileRead( am_I_Root, FileDta,    MW_g,    &
                                EmMW_g,    MolecRatio, TS_EMIS, &
                                IsFirst,   RC                    ) 
!
! !USES:
!
    USE HCO_CHARTOOLS_MOD,  ONLY : HCO_CharSplit
    USE HCO_CHARTOOLS_MOD,  ONLY : HCO_WCD, HCO_SEP
    USE HCO_UNIT_MOD,       ONLY : HCO_Unit_Change
!
! !INPUT PARAMETERS:
!
    LOGICAL,         INTENT(IN   ) :: am_I_Root  ! Root CPU?
    TYPE(FileData),  POINTER       :: FileDta    ! File data obj
    REAL(hp),        INTENT(IN   ) :: MW_g       ! input data MW 
    REAL(hp),        INTENT(IN   ) :: EmMW_g     ! emission data MW
    REAL(hp),        INTENT(IN   ) :: MolecRatio ! molec. ratio
    REAL(sp),        INTENT(IN   ) :: TS_EMIS    ! emission time step [s] 
!
! !OUTPUT PARAMETERS:
!
    INTEGER,       INTENT(  OUT) :: IsFirst   ! First read? 
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,       INTENT(INOUT) :: RC        ! Return code
!
! !REVISION HISTORY:
!  10 Jan 2014 - C. Keller - Initialization (update)
!  21 Aug 2014 - C. Keller - Added concentration
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER            :: I, N, AS
    INTEGER            :: AreaFlag, TimeFlag, Check
    CHARACTER(LEN=255) :: MSG, LOC
    REAL(sp)           :: FileVals(24)
    REAL(sp), POINTER  :: FileArr(:,:,:,:) => NULL()
    LOGICAL            :: Verb, IsPerArea

    !======================================================================
    ! FileData_FileRead begins here
    !======================================================================

    ! Enter
    LOC = 'FileData_FileRead (hco_filedata_mod.F90)'
    Verb = HCO_VERBOSE_CHECK() .and. am_I_Root

    ! Don't read twice. Data file objects can belong to multiple 
    ! containers, and this routine may hence be called multiple
    ! times for the same data object. Make sure that we read the
    ! data only on first call! The variable SpaceDim is initialized
    ! to -1 but becomes set to one after reading, hence use this
    ! variable to identify whether or not this is the first call
    ! for this object.
    IF ( FileDta%SpaceDim == 1 ) THEN 
       IsFirst = 0
       RC      = HCO_SUCCESS
       RETURN
    ENDIF

    ! Verbose
    IF ( Verb ) THEN
       WRITE(MSG, *) 'Read from config file: ', TRIM(FileDta%ncFile)
       CALL HCO_MSG(MSG)
    ENDIF

    ! Read data into array
    CALL HCO_CharSplit ( FileDta%ncFile, HCO_SEP(), &
                         HCO_WCD(), FileVals, N, RC ) 
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Return w/ error if no scale factor defined
    IF ( N == 0 ) THEN
       MSG = 'Cannot read data: ' // TRIM(FileDta%ncFile)
       CALL HCO_ERROR ( MSG, RC, THISLOC=LOC)
       RETURN 
    ENDIF

    ! Convert data to HEMCO units 
    ALLOCATE( FileArr(1,1,1,N), STAT=AS )
    IF ( AS /= 0 ) THEN
       MSG = 'Cannot allocate FileArr'
       CALL HCO_ERROR ( MSG, RC, THISLOC=LOC )
       RETURN
    ENDIF
    FileArr(1,1,1,:) = FileVals(:)
    CALL HCO_Unit_Change( Array       = FileArr,                &
                          Units       = TRIM(FileDta%OrigUnit), &
                          MW_IN       = MW_g,                   &
                          MW_OUT      = EmMW_g,                 &
                          MOLEC_RATIO = MolecRatio,             &
                          YYYY        = -999,                   &
                          MM          = -999,                   &
                          AreaFlag    = AreaFlag,               &
                          TimeFlag    = TimeFlag,               &
                          RC          = RC                       )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Data must be ... 
    ! ... concentration ...
    IF ( AreaFlag == 3 .AND. TimeFlag == 0 ) THEN
       FileDta%IsConc = .TRUE.

    ELSEIF ( AreaFlag == 3 .AND. TimeFlag == 1 ) THEN
       FileDta%IsConc = .TRUE.
       FileArr = FileArr * TS_EMIS
       MSG = 'Data converted from kg/m3/s to kg/m3: ' // &
             TRIM(FileDta%ncFile) // ': ' // TRIM(FileDta%OrigUnit)
       CALL HCO_WARNING ( MSG, RC, THISLOC=LOC )

    ! ... emissions or unitless ...
    ELSEIF ( (AreaFlag == -1 .AND. TimeFlag == -1) .OR. &
             (AreaFlag ==  2 .AND. TimeFlag ==  1)       ) THEN
       FileDta%IsConc = .FALSE.

    ! ... invalid otherwise:
    ELSE
       MSG = 'Unit must be unitless, emission or concentration ' // &
             TRIM(FileDta%OrigUnit)
       CALL HCO_ERROR ( MSG, RC, THISLOC=LOC )
       RETURN
    ENDIF

    ! Copy data into array. Assume all data is temporal
    ! dimension!
    CALL FileData_ArrCheck2D( FileDta, 1, 1, N, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN
    DO I = 1, N
       FileDta%V2(I)%Val(1,1) = FileArr(1,1,1,I)
    ENDDO

    ! Make sure that ncRead flag is turned off.
    FileDta%ncRead   = .FALSE.

    ! Data is always 1D.
    FileDta%SpaceDim = 1

    ! Auto-detect delta t [in hours] between time slices.
    ! Scale factors can be:
    ! length 1 : constant
    ! length 7 : weekday factors: Sun, Mon, ..., Sat
    ! length 12: monthly factors: Jan, Feb, ..., Dec
    ! length 24: hourly  factors: 12am, 1am, ... 11pm
    IF ( N == 1 ) THEN
       FileDta%DeltaT = 0
    ELSEIF ( N == 7 ) THEN
       FileDta%DeltaT = 24
    ELSEIF ( N == 12 ) THEN
       FileDta%DeltaT = 720 
    ELSEIF ( N == 24 ) THEN
       FileDta%DeltaT = 1 
    ELSE
       MSG = 'Factor must be of length 1, 7, 12, or 24!' // &
              TRIM(FileDta%ncFile)
       CALL HCO_ERROR ( MSG, RC, THISLOC=LOC)
       RETURN 
    ENDIF

    ! Cleanup
    IF ( ASSOCIATED(FileArr) ) DEALLOCATE(FileArr)
    FileArr => NULL()

    ! Return w/ success
    IsFirst = 1
    RC      = HCO_SUCCESS 

  END SUBROUTINE FileData_FileRead
!EOC
END MODULE HCO_FileData_Mod
