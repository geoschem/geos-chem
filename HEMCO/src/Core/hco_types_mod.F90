!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: hco_types_mod.F90
!
! !DESCRIPTION: Module HCO\_Types\_Mod contains HEMCO derived type definitions
!  and global parameters.
!  The here derived types become all bundled into the HEMCO state object, which
!  is defined in module hco\_state\_mod.F90. Specific routines for the various
!  derived types defined below can be found in the respective HEMCO modules
!  (e.g. hco\_datacont\_mod.F90 for handling of HEMCO data containers and data
!  container lists).
!\\
!\\
!  Prior to HEMCO v2.0, the derived type definitions were placed within the
!  corresponding modules. Collecting them here breaks unnecessary dependencies
!  amongst the various modules and consequently gives greater flexibility in
!  using the modules independently.
!\\
!\\
! !INTERFACE:
!
MODULE HCO_TYPES_MOD
!
! USES:
!
  USE HCO_Error_Mod
  USE HCO_Arr_Mod

  IMPLICIT NONE
  PUBLIC
!
! !PUBLIC PARAMETERS:
!
  ! Maximum length of option name
  INTEGER, PARAMETER       :: OPTLEN = 1023

  ! Cycle flags. Used to determine the temporal behavior of data
  ! fields. For example, data set to 'cycle' will be recycled if
  ! the simulation date is outside of the datetime range of the
  ! source data. See data reading routine (hcoio_dataread_mod.F90)
  ! for more details.
  INTEGER, PARAMETER, PUBLIC  :: HCO_CFLAG_CYCLE    = 1
  INTEGER, PARAMETER, PUBLIC  :: HCO_CFLAG_RANGE    = 2
  INTEGER, PARAMETER, PUBLIC  :: HCO_CFLAG_EXACT    = 3
  INTEGER, PARAMETER, PUBLIC  :: HCO_CFLAG_INTER    = 4
  INTEGER, PARAMETER, PUBLIC  :: HCO_CFLAG_AVERG    = 5
  INTEGER, PARAMETER, PUBLIC  :: HCO_CFLAG_RANGEAVG = 6

  ! Data container update flags. At the moment, those only indicate
  ! if a container gets updated every time step or based upon the
  ! details specifications ('srcTime') in the HEMCO configuration file.
  INTEGER, PARAMETER, PUBLIC  :: HCO_UFLAG_FROMFILE = 1
  INTEGER, PARAMETER, PUBLIC  :: HCO_UFLAG_ALWAYS   = 2
  INTEGER, PARAMETER, PUBLIC  :: HCO_UFLAG_ONCE     = 3
  INTEGER, PARAMETER, PUBLIC  :: HCO_UFLAG_3HR      = 4

  ! Data container types. These are used to distinguish between
  ! base emissions, scale factors and masks.
  INTEGER, PARAMETER, PUBLIC  :: HCO_DCTTYPE_BASE = 1
  INTEGER, PARAMETER, PUBLIC  :: HCO_DCTTYPE_SCAL = 2
  INTEGER, PARAMETER, PUBLIC  :: HCO_DCTTYPE_MASK = 3

  ! PBL level flag
  INTEGER, PARAMETER, PUBLIC  :: HCO_EMISL_PBL = 0
  INTEGER, PARAMETER, PUBLIC  :: HCO_EMISL_LEV = 1
  INTEGER, PARAMETER, PUBLIC  :: HCO_EMISL_M   = 2
!
! !PUBLIC TYPES:
!
  !=========================================================================
  ! HcoSpc: Derived type for HEMCO species
  !
  ! Notes:
  ! **1 The emission molecular weight is the molecular weight of the
  !     emitted compound. This value is only different to MW_g if the
  !     emitted compound does not correspond to the transported species,
  !     e.g. if emissions are in kg C4H10 but the corresponding species
  !     is transported as mass Carbon.
  ! **2 MolecRatio is the ratio between # of species molecules per emitted
  !       molecule, e.g. 4 if emissions are kg C4H10 but model species
  !       are kg C.
  !=========================================================================
  TYPE :: HcoSpc
     INTEGER                 :: HcoID      ! HEMCO species ID
     INTEGER                 :: ModID      ! Model species ID
     CHARACTER(LEN= 31)      :: SpcName    ! species names
     REAL(hp)                :: MW_g       ! species molecular wt.     (g/mol)
     REAL(hp)                :: EmMW_g     ! emission molecular wt.**1 (g/mol)
     REAL(hp)                :: MolecRatio ! molecule emission ratio**2 (-)
     REAL(hp)                :: HenryK0    ! liq. over gas Henry const [M/atm]
     REAL(hp)                :: HenryCR    ! K0 temp. dependency [K]
     REAL(hp)                :: HenryPKA   ! pKa for Henry const. correction
     TYPE(Arr2D_HP), POINTER :: Depv       ! Deposition velocity [1/s]
     TYPE(Arr3D_HP), POINTER :: Emis       ! Emission flux [kg/m2/s]
     TYPE(Arr3D_HP), POINTER :: Conc       ! Concentration [v/v]
  END TYPE HcoSpc

  !=========================================================================
  ! ModSpc: Derived type for model species
  !=========================================================================
  TYPE :: ModSpc
     INTEGER                 :: HcoID      ! HEMCO species ID
     INTEGER                 :: ModID      ! Model species ID
     CHARACTER(LEN= 31)      :: SpcName    ! species names
  END TYPE ModSpc

  !=========================================================================
  ! HcoOpt: Derived type for HEMCO run options
  !=========================================================================
  TYPE :: HcoOpt
     INTEGER  :: ExtNr          ! ExtNr to be used
     INTEGER  :: SpcMin         ! Smallest HEMCO species ID to be considered
     INTEGER  :: SpcMax         ! Highest HEMCO species ID to be considered
     INTEGER  :: CatMin         ! Smallest category to be considered
     INTEGER  :: CatMax         ! Highest category to be considered
     LOGICAL  :: HcoWritesDiagn ! If set to .TRUE., HEMCO will schedule the
                                ! output of the default HEMCO diagnostics
                                ! (in hco_driver_mod.F90).
     LOGICAL  :: AutoFillDiagn  ! Write into AutoFill diagnostics?
     LOGICAL  :: FillBuffer     ! Write calculated emissions into buffer
                                ! instead of emission array?
     INTEGER  :: NegFlag        ! Negative value flag (from configfile):
                                ! 2 = allow negative values
                                ! 1 = set neg. values to zero and prompt warning
                                ! 0 = return w/ error if neg. value
     LOGICAL  :: PBL_DRYDEP     ! If true, dry deposition frequencies will
                                ! be calculated over the full PBL. If false,
                                ! they are calculated over the first layer only.
     REAL(hp) :: MaxDepExp      ! Maximum value of deposition freq. x time step.
     LOGICAL  :: MaskFractions  ! If TRUE, masks are treated as binary, e.g.
                                ! grid boxes are 100% inside or outside of a
                                ! mask.
     LOGICAL  :: Field2Diagn    ! When reading fields from disk, check if there
                                ! is a diagnostics with the same name and write
                                ! field to that diagnostics? Defaults to yes in
                                ! standalone mode and no in other setups.
     LOGICAL  :: VertWeight     ! if spreading 2D fields across multiple vert.
                                ! levels, weight vertical dilution factors based
                                ! upon level depths?
     LOGICAL  :: ScaleEmis      ! Scale emissions by uniform scale factors set
                                ! in HEMCO configuration file? Defaults to yes.
     LOGICAL  :: TimeShiftCap   ! Cap time shift to same day. Defaults to no.
     LOGICAL  :: isESMF         ! Are we using ESMF?
     LOGICAL  :: isDryRun       ! Are we in a dry run?
  END TYPE HcoOpt

  !=========================================================================
  ! VertGrid: Description of vertical grid
  !=========================================================================
  TYPE :: VertGrid
     INTEGER                 :: ZTYPE           ! Grid type
     REAL(hp),       POINTER :: Ap(:) => NULL() ! Hybrid sigma Ap values
     REAL(hp),       POINTER :: Bp(:) => NULL() ! Hybrid sigma Bp values
  END TYPE VertGrid

  !=========================================================================
  ! HcoGrid: Derived type for HEMCO grid. The grid edges are used for data
  ! interpolation. The pressure midpoints are not needed by HEMCO core but
  ! can be specified for the extensions through ExtList.
  !
  ! NOTES:
  ! *  Not used in ESMF environment
  ! ** Only used by some extensions
  !=========================================================================
  TYPE :: HcoGrid
     TYPE(Arr2D_Hp), POINTER :: XMID       ! mid-points in x-direction (lon)
     TYPE(Arr2D_Hp), POINTER :: YMID       ! mid-points in y-direction (lat)
     TYPE(Arr2D_Hp), POINTER :: XEDGE      ! grid edges in x-direction (lon)*
     TYPE(Arr2D_Hp), POINTER :: YEDGE      ! grid edges in y-direction (lat)*
     TYPE(Arr3D_Hp), POINTER :: PEDGE      ! pressure edges (Pa)
     TYPE(Arr2D_Hp), POINTER :: YSIN       ! sin of y-direction grid edges*
     TYPE(Arr2D_Hp), POINTER :: AREA_M2    ! grid box areas (m2)
     TYPE(Arr2D_Hp), POINTER :: ZSFC       ! surface geopotential height (m)**
     TYPE(Arr2D_Hp), POINTER :: PSFC       ! surface pressure (Pa)
     TYPE(Arr2D_Hp), POINTER :: PBLHEIGHT  ! PBL height in m
     TYPE(Arr3D_Hp), POINTER :: BXHEIGHT_M ! grid box heights (m)**
     TYPE(VertGrid), POINTER :: ZGRID      ! vertical grid description
  END TYPE HcoGrid

  !=========================================================================
  ! HcoClock: Derived type definition for the HEMCO clock object
  !=========================================================================
  TYPE :: HcoClock

     ! Current time stamp (UTC)
     INTEGER            :: ThisYear        ! year
     INTEGER            :: ThisMonth       ! month
     INTEGER            :: ThisDay         ! day
     INTEGER            :: ThisHour        ! hour
     INTEGER            :: ThisMin         ! minute
     INTEGER            :: ThisSec         ! second
     INTEGER            :: ThisDOY         ! day of year
     INTEGER            :: ThisWD          ! weekday (0=Sun,...,6=Sat)
     INTEGER            :: MonthLastDay    ! Last day of month: 28,29,30,31

     ! Current local times
     INTEGER            :: ntz             ! number of time zones
     INTEGER,  POINTER  :: ThisLocYear(:)  ! local year
     INTEGER,  POINTER  :: ThisLocMonth(:) ! local month
     INTEGER,  POINTER  :: ThisLocDay(:)   ! local day
     INTEGER,  POINTER  :: ThisLocWD(:)    ! local weekday
     REAL(sp), POINTER  :: ThisLocHour(:)  ! local hour

     ! Previous time stamp (UTC)
     INTEGER            :: PrevYear
     INTEGER            :: PrevMonth
     INTEGER            :: PrevDay
     INTEGER            :: PrevHour
     INTEGER            :: PrevMin
     INTEGER            :: PrevSec
     INTEGER            :: PrevDOY
     INTEGER            :: PrevWD

     ! Current emission time stamp
     INTEGER            :: ThisEYear
     INTEGER            :: ThisEMonth
     INTEGER            :: ThisEDay
     INTEGER            :: ThisEHour
     INTEGER            :: ThisEMin
     INTEGER            :: ThisESec

     ! Previous emission time stamp
     INTEGER            :: PrevEYear
     INTEGER            :: PrevEMonth
     INTEGER            :: PrevEDay
     INTEGER            :: PrevEHour
     INTEGER            :: PrevEMin
     INTEGER            :: PrevESec

     ! Simulation year, month, day, hour, minute, second. Will only
     ! be different from current time stamp in special cases, e.g.
     ! if emission year shall be fixed.
     INTEGER            :: SimYear        ! year
     INTEGER            :: SimMonth       ! month
     INTEGER            :: SimDay         ! day
     INTEGER            :: SimHour        ! hour
     INTEGER            :: SimMin         ! minute
     INTEGER            :: SimSec         ! second

     ! total number of elapsed time steps and emission time steps
     ! LastEStep denotes the last nEmisSteps values for which
     ! emissions have been calculated.
     INTEGER            :: nSteps
     INTEGER            :: nEmisSteps
     INTEGER            :: LastEStep

     ! Pointer to gridded time zones. Will only hold data if a field `TIMEZONES`
     ! is provided in the HEMCO configuration file.
     !
     ! NOTE: This pointer is initialized by a call to HcoClock_InitTzPtr
     ! from the HEMCO run routine (HCO_Run, in hco_driver_mod.F90).
     ! This is necessary to avoid segmentation faults when running with
     ! OpenMP turned on. (bmy, 2/23/15)
     TYPE(Arr2D_Sp), POINTER :: TIMEZONES

     ! Fixed dates to be used for simulation dates
     INTEGER              :: FixYY    = -1
     INTEGER              :: FixMM    = -1
     INTEGER              :: Fixdd    = -1
     INTEGER              :: Fixhh    = -1

     ! Last time step?
     LOGICAL              :: isLast

  END TYPE HcoClock

  !=========================================================================
  ! HcoPhys: Derived type for HEMCO physical constants
  !=========================================================================
  TYPE :: HcoPhys
     REAL(dp) :: Avgdr   ! Avogadro number (mol-1)
     REAL(dp) :: PI      ! Pi
     REAL(dp) :: PI_180  ! Pi / 180
     REAL(dp) :: Re      ! Earth radius [m]
     REAL(dp) :: AIRMW   ! Molecular weight of air (g/mol)
     REAL(dp) :: g0      ! Gravity at surface of earth (m/s2)
     REAL(dp) :: Rd      ! Gas Constant (R) in dry air (J/K/kg)
     REAL(dp) :: Rdg0    ! Rd/g0
     REAL(dp) :: RSTARG  ! Universal gas constant [J/K/mol]
  END TYPE HcoPhys

  !=========================================================================
  ! HcoMicroPhys: Derived type for aerosol microphysics settings
  !=========================================================================
  TYPE :: HcoMicroPhys
     INTEGER           :: nBins              ! # of size-resolved bins
     INTEGER           :: nActiveModebins    ! # of active mode bins
     REAL(dp), POINTER :: BinBound(:)        ! Size bin boundaries
  END TYPE HcoMicroPhys

  ! RdList contains lists to all ReadLists
  TYPE :: RdList
     TYPE(ListCont), POINTER :: Once
     TYPE(ListCont), POINTER :: Always
     TYPE(ListCont), POINTER :: Year
     TYPE(ListCont), POINTER :: Month
     TYPE(ListCont), POINTER :: Day
     TYPE(ListCont), POINTER :: Hour
     TYPE(ListCont), POINTER :: Hour3
     INTEGER                 :: FileLun       = -1  ! LUN of file in archive
     CHARACTER(LEN=2023)     :: FileInArchive = ''  ! name of file in archive
     INTEGER                 :: Counter       =  0  ! ReadList read counter
  END TYPE RdList

  !-------------------------------------------------------------------------
  ! DataCont: Derived type definition for HEMCO data container
  !-------------------------------------------------------------------------
  TYPE :: DataCont

     ! Container information
     CHARACTER(LEN= 63)          :: cName          ! Cont. name
     INTEGER                     :: cID            ! Cont. ID
     INTEGER                     :: targetID       ! target ID
     INTEGER                     :: DctType        ! Data type
     TYPE(FileData),     POINTER :: Dta            ! data information
     INTEGER                     :: DtaHome        ! Home cont for Dta?
     CHARACTER(LEN= 31)          :: SpcName        ! Species Name
     INTEGER                     :: HcoID          ! HEMCO species ID
     INTEGER                     :: ExtNr          ! Extension #
     INTEGER                     :: Cat            ! Category
     INTEGER                     :: Hier           ! Hierarchy
     INTEGER                     :: ScalID         ! Scale factor ID
     INTEGER                     :: Oper           ! Operator
     INTEGER                     :: levScalID1     ! ID of vertical level field
     INTEGER                     :: levScalID2     ! ID of vertical level field
     INTEGER                     :: nScalID        ! # of scale factor IDs
     INTEGER,            POINTER :: Scal_cID(:)    ! assoc. scalefactor IDs
     LOGICAL                     :: Scal_cID_set   ! cIDs or scalIDs
  END TYPE DataCont

  !-------------------------------------------------------------------------
  ! ListCont: Derived type definition for HEMCO list object
  !-------------------------------------------------------------------------
  TYPE :: ListCont
     TYPE(DataCont),     POINTER :: Dct
     TYPE(ListCont),     POINTER :: NextCont
  END TYPE ListCont

  !-------------------------------------------------------------------------
  ! FileData: Derived type definition for HEMCO filetype object
  !-------------------------------------------------------------------------
  TYPE :: FileData
     CHARACTER(LEN=255)          :: ncFile    ! file path+name
     CHARACTER(LEN= 50)          :: ncPara    ! file parameter
     INTEGER                     :: ncYrs(2)  ! year range
     INTEGER                     :: ncMts(2)  ! month range
     INTEGER                     :: ncDys(2)  ! day range
     INTEGER                     :: ncHrs(2)  ! hour range
     INTEGER                     :: tShift(6) ! time stamp shift (YMDhms)
     INTEGER                     :: CycleFlag ! cycle flag
     LOGICAL                     :: MustFind  ! file must be found
     LOGICAL                     :: UseSimYear! use simulation year
     LOGICAL                     :: Discontinuous ! discontinuous dataset?
     INTEGER                     :: UpdtFlag  ! update flag
     LOGICAL                     :: ncRead    ! read from source?
     TYPE(Arr3D_SP),     POINTER :: V3(:)     ! vector of 3D fields
     TYPE(Arr2D_SP),     POINTER :: V2(:)     ! vector of 2D fields
     TYPE(TimeIdx),      POINTER :: tIDx      ! for time slice indexing
     CHARACTER(LEN= 31)          :: OrigUnit  ! original data units
     CHARACTER(LEN= 63)          :: ArbDimName! name of additional dimension
     CHARACTER(LEN= 63)          :: ArbDimVal ! desired value of additional dimension
     INTEGER                     :: Cover     ! data coverage
     INTEGER                     :: SpaceDim  ! space dimension: 1, 2 or 3
     INTEGER                     :: Levels    ! vertical level handling
!     INTEGER                     :: Lev2D     ! level to use for 2D data
     REAL(hp)                    :: EmisL1     ! emission level 1
     REAL(hp)                    :: EmisL2     ! emission level 2
     INTEGER                     :: EmisL1Unit ! emission level 1 unit
     INTEGER                     :: EmisL2Unit ! emission level 2 unit
     INTEGER                     :: nt        ! time dimension: length of Arr
     INTEGER                     :: DeltaT    ! temp. resolution of array [h]
     LOGICAL                     :: IsLocTime ! local time?
     LOGICAL                     :: IsConc    ! concentration data?
     LOGICAL                     :: DoShare   ! shared object?
     LOGICAL                     :: IsInList  ! is in emissions list?
     LOGICAL                     :: IsTouched ! Has container been touched yet?
  END TYPE FileData

  !-------------------------------------------------------------------------
  ! TimeIdx: Derived type definition for the object that points to the
  ! current time slices of data within a file.  Used by hco_tidx_mod.F90.
  !-------------------------------------------------------------------------
  TYPE :: TimeIdx
     INTEGER                     :: TypeID
     CHARACTER(LEN=31)           :: TempRes
  END TYPE TimeIdx

  ! The TimeIdxCollection derived type contains the pointers with the
  ! current valid vector indices for all defined cycling intervals.
  TYPE ::  TimeIdxCollection
     TYPE(TimeIdx), POINTER :: CONSTANT
     TYPE(TimeIdx), POINTER :: HOURLY
     TYPE(TimeIdx), POINTER :: HOURLY_GRID
     TYPE(TimeIdx), POINTER :: WEEKDAY
     TYPE(TimeIdx), POINTER :: MONTHLY
  END TYPE TimeIdxCollection

  !-------------------------------------------------------------------------
  ! cIdListPnt: Derived type definition for pointing to list containers
  !-------------------------------------------------------------------------
  TYPE :: cIDListPnt
     TYPE(DataCont),     POINTER :: PNT ! Pointer to list container
  END TYPE cIDListPnt

  !-------------------------------------------------------------------------
  ! Variables to store (unique) scale factor IDs and species names
  !-------------------------------------------------------------------------
  TYPE :: ScalIDCont
     INTEGER                   :: ScalID
     TYPE(ScalIDCont), POINTER :: NEXT
  END TYPE ScalIDCont

  TYPE :: SpecNameCont
     CHARACTER(LEN=31)           :: SpecName
     TYPE(SpecNameCont), POINTER :: NEXT
  END TYPE SpecNameCont

  !-------------------------------------------------------------------------
  ! Derived type to store options
  !-------------------------------------------------------------------------
  TYPE :: Opt
     CHARACTER(LEN=OPTLEN)   :: OptName
     CHARACTER(LEN=OPTLEN)   :: OptValue
     TYPE(Opt), POINTER      :: NextOpt => NULL()
  END TYPE Opt

  !-------------------------------------------------------------------------
  ! Derived type to manage list of extensions and associated options
  !-------------------------------------------------------------------------
  TYPE :: Ext
     CHARACTER(LEN=255)       :: ExtName  ! Name
     CHARACTER(LEN=OPTLEN)    :: Spcs     ! Species
     INTEGER                  :: ExtNr    ! Ext. number
     TYPE(Opt), POINTER       :: Opts     ! Options linked list
     TYPE(Ext), POINTER       :: NextExt  ! next list item
  END TYPE Ext

  !-------------------------------------------------------------------------
  ! Configuration object
  !-------------------------------------------------------------------------
  TYPE :: ConfigObj
     CHARACTER(LEN=OPTLEN)        :: ROOT           = ''
     CHARACTER(LEN=255)           :: ConfigFileName = ''
     TYPE(ScalIDCont),   POINTER  :: ScalIDList     => NULL()
     TYPE(SpecNameCont), POINTER  :: SpecNameList   => NULL()
     TYPE(ListCont),     POINTER  :: ConfigList     => NULL()
     TYPE(Ext),          POINTER  :: ExtList        => NULL()
     TYPE(HcoErr),       POINTER  :: Err            => NULL()
     TYPE(ModSpc),       POINTER  :: ModelSpc(:)
     INTEGER                      :: nModelSpc
     INTEGER                      :: nModelAdv
     CHARACTER(LEN=255)           :: MetField
     CHARACTER(LEN=255)           :: GridRes
     LOGICAL                      :: ConfigFileRead = .FALSE.
     LOGICAL                      :: amIRoot       ! Is this the root CPU?
  END TYPE ConfigObj

  !------------------------------------------------------------------------
  ! DiagnCont: Diagnostics container derived type declaration
  !------------------------------------------------------------------------
  TYPE :: DiagnCont
     CHARACTER(LEN= 63)          :: cName          ! Cont. name
     CHARACTER(LEN=255)          :: long_name      ! ncdf long_name attribute
     INTEGER                     :: cID            ! Cont. ID
     INTEGER                     :: ExtNr          ! Extension #
     INTEGER                     :: Cat            ! Category
     INTEGER                     :: Hier           ! Hierarchy
     INTEGER                     :: HcoID          ! HEMCO species ID
     INTEGER                     :: AutoFill       ! fill automatically?
     INTEGER                     :: SpaceDim       ! Space dimension (1-3)
     REAL(sp)                    :: Scalar         ! 1D scalar
     TYPE(Arr2D_SP),     POINTER :: Arr2D          ! 2D array
     TYPE(Arr3D_SP),     POINTER :: Arr3D          ! 3D array
     REAL(sp)                    :: Total          ! Diagnostics total
     LOGICAL                     :: DtaIsPtr       ! Is data just a pointer?
     INTEGER                     :: LevIdx         ! Level index to be used
     CHARACTER(LEN= 31)          :: OutUnit        ! Output unit
     INTEGER                     :: AreaFlag       ! 2=per area, 3=per volume, 0 otherwise
     REAL(hp)                    :: AreaScal       ! Scale factor for area
     REAL(hp)                    :: MassScal       ! Scale factor for mass
     REAL(hp)                    :: ScaleFact      ! Uniform scale factor
     INTEGER                     :: TimeAvg        ! Scale flag for time unit
     INTEGER                     :: Counter        ! time steps since
                                                   ! last output
     CHARACTER(LEN= 31)          :: AvgName        ! Output averaging operation
     INTEGER                     :: AvgFlag        ! Averaging flag for
                                                   !  non-standard units
     INTEGER                     :: LastUpdateID   ! Last update time
     INTEGER                     :: nnGetCalls     ! # of Diagn_Get calls w/o update
     LOGICAL                     :: IsOutFormat    ! Data is in output format?
     INTEGER                     :: CollectionID   ! Collection diagnostics belongs to
     TYPE(DiagnCont),    POINTER :: NextCont       ! Ptr to next item in list
  END TYPE DiagnCont

  !------------------------------------------------------------------------
  ! Diagnostcs collection derived type.
  ! DiagnList      : Linked list with all diagnostics container of
  !                  this collection.
  ! nnDiag         : Number of diagnostics in this collection.
  ! AF_LevelDefined: Set to true if there is at least one autofill
  !                  diagnostics at the given level (1-4).
  ! PREFIX         : Prefix to be used for diagnostics output file name.
  ! NX, NY, NZ     : Grid dimensions.
  ! TS             : Time step. This is only of relevance for emission
  !                  diagnostics that are internally converted from
  !                  kg/m2/s to kg/m2.
  ! AREA_M2        : Surface grid box areas. May be required for unit
  !                  conversions.
  !------------------------------------------------------------------------
  TYPE :: DiagnCollection
     TYPE(DiagnCont),       POINTER :: DiagnList          => NULL()
     INTEGER                        :: nnDiagn            =  0
     LOGICAL                        :: AF_LevelDefined(4) =  .FALSE.
     INTEGER                        :: CollectionID       = -1
     CHARACTER(LEN=255)             :: PREFIX             =  ''
     INTEGER                        :: NX                 =  0
     INTEGER                        :: NY                 =  0
     INTEGER                        :: NZ                 =  0
     INTEGER                        :: deltaYMD           =  0
     INTEGER                        :: deltaHMS           =  0
     INTEGER                        :: lastYMD            = -1
     INTEGER                        :: lastHMS            = -1
     INTEGER                        :: OutTimeStamp       = -1
     REAL(sp)                       :: TS                 =  0       ! Time step
     REAL(hp),              POINTER :: AREA_M2(:,:)       => NULL()
     TYPE(DiagnCollection), POINTER :: NextCollection     => NULL()
  END TYPE DiagnCollection

  TYPE :: DiagnBundle
     TYPE(DiagnCollection), POINTER :: Collections       => NULL()
     INTEGER                        :: HcoDiagnIDDefault = -999
     INTEGER                        :: HcoDiagnIDRestart = -999
     INTEGER                        :: HcoDiagnIDManual  = -999
     INTEGER                        :: nnCollections     = 0
  END TYPE DiagnBundle
!
! !REVISION HISTORY:
!  15 Feb 2016 - C. Keller   - Initial version (collected from various modules)
!  12 May 2017 - C. Keller   - Added option ScaleEmis
!  05 Oct 2018 - R. Yantosca - Added HCO_UFLAG_ONCE parameter
!  23 Oct 2018 - M. Sulprizio- Added derived type for external model species
!                              to ConfigObj to facilitate reading GEOS-Chem
!                              restart file via HEMCO.
!  07 Feb 2019 - C. Keller   - Added ReadList read counter.
!EOP
!------------------------------------------------------------------------------
!BOC
!EOC
END MODULE HCO_TYPES_MOD
