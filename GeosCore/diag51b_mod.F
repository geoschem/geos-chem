!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: diag51b_mod 
!
! !DESCRIPTION: Module DIAG51b\_MOD contains variables and routines to 
!  generate save timeseries data where the local time is between two 
!  user-defined limits. This facilitates comparisons with morning or 
!  afternoon-passing satellites such as GOME.
!\\
!\\
! !INTERFACE: 
!
      MODULE DIAG51b_MOD
!
! !USES:
!
      IMPLICIT NONE
      PRIVATE
!
! !PUBLIC DATA MEMBERS:
!
      LOGICAL, PUBLIC :: DO_SAVE_DIAG51b   ! On/off switch for ND51b diagnostic
!
! !PUBLIC MEMBER FUNCTIONS:
! 
      PUBLIC  :: CLEANUP_DIAG51b
      PUBLIC  :: DIAG51b
      PUBLIC  :: INIT_DIAG51b
!
! !PRIVATE MEMBER FUNCTIONS:
! 
      PRIVATE :: ACCUMULATE_DIAG51        
      PRIVATE :: GET_LOCAL_TIME           
      PRIVATE :: ITS_TIME_FOR_WRITE_DIAG51
      PRIVATE :: WRITE_DIAG51      
!
! !REMARKS:
!  ND51b tracer numbers:
!  ============================================================================
!  1 - N_TRACERS : GEOS-CHEM transported tracers            [v/v        ]
!  76            : OH concentration                         [molec/cm3  ]
!  77            : NO2 concentration                        [v/v        ]
!  78            : PBL heights                              [m          ]
!  79            : PBL heights                              [levels     ]
!  80            : Air density                              [molec/cm3  ]
!  81            : 3-D Cloud fractions                      [unitless   ]
!  82            : Column optical depths                    [unitless   ]
!  83            : Cloud top heights                        [hPa        ]
!  84            : Sulfate aerosol optical depth            [unitless   ]
!  85            : Black carbon aerosol optical depth       [unitless   ]
!  86            : Organic carbon aerosol optical depth     [unitless   ]
!  87            : Accumulation mode seasalt optical depth  [unitless   ]
!  88            : Coarse mode seasalt optical depth        [unitless   ]
!  89            : Total dust optical depth                 [unitless   ]
!  90            : Total seasalt tracer concentration       [unitless   ]
!  91            : Pure O3 (not Ox) concentration           [v/v        ]
!  92            : NO concentration                         [v/v        ]
!  93            : NOy concentration                        [v/v        ]
!  94            : Grid box heights                         [m          ]
!  95            : Relative Humidity                        [%          ]
!  96            : Sea level pressure                       [hPa        ]
!  97            : Zonal wind (a.k.a. U-wind)               [m/s        ]
!  98            : Meridional wind (a.k.a. V-wind)          [m/s        ]
!  99            : P(surface) - PTOP                        [hPa        ]
!  100           : Temperature                              [K          ]
!  101           : PAR direct                               [hPa        ]
!  102           : PAR diffuse                              [hPa        ]
!  103           : Daily LAI                                [hPa        ]
!  104           : Temperature at 2m                        [K          ]
!  105           : Isoprene emissions                       [atomC/cm2/s]
!  106           : Total Monoterpene emissions              [atomC/cm2/s]
!  107           : Methyl Butanol emissions                 [atomC/cm2/s]
!  108           : Alpha-Pinene emissions                   [atomC/cm2/s]
!  109           : Beta-Pinene emissions                    [atomC/cm2/s]
!  110           : Limonene emissions                       [atomC/cm2/s]
!  111           : Sabinene emissions                       [atomC/cm2/s]
!  112           : Myrcene emissions                        [atomC/cm2/s]
!  113           : 3-Carene emissions                       [atomC/cm2/s]
!  114           : Ocimene emissions                        [atomC/cm2/s]
!  115-121       : size resolved dust optical depth         [unitless   ]
!
! !REVISION HISTORY:
!  (1 ) Rewritten for clarity (bmy, 7/20/04)
!  (2 ) Added extra counters for NO, NO2, OH, O3.  Also all diagnostic counter
!        arrays are 1-D since they only depend on longitude. (bmy, 10/25/04)
!  (3 ) Bug fix: Now get I0 and J0 properly for nested grids (bmy, 11/9/04)
!  (4 ) Now only archive AOD's once per chemistry timestep (bmy, 1/14/05)
!  (5 ) Now references "pbl_mix_mod.f" (bmy, 2/16/05)
!  (6 ) Now save cld frac and grid box heights (bmy, 4/20/05)
!  (7 ) Remove TRCOFFSET since it's always zero  Also now get HALFPOLAR for
!        both GCAP and GEOS grids.  (bmy, 6/28/05)
!  (8 ) Bug fix: do not save SLP if it's not allocated (bmy, 8/2/05)
!  (9 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (10) Now references XNUMOLAIR from "tracer_mod.f" (bmy, 10/25/05)
!  (11) Modified INIT_DIAG51 to save out transects (cdh, bmy, 11/30/06)
!  (12) Now use 3D timestep counter for full chem in the trop (phs, 1/24/07)
!  (13) Renumber RH in WRITE_DIAG50 (bmy, 2/11/08)
!  (14) Bug fix: replace "PS-PTOP" with "PEDGE-$" (bmy, phs, 10/7/08)
!  (15) Bug fix in GET_LOCAL_TIME (ccc, 12/10/08)
!  (16) Modified to archive O3, NO, NOy as tracers 89, 90, 91  (tmf, 9/26/07)
!  (17) Updates in WRITE_DIAG51b (ccc, tai, bmy, 10/13/09)
!  (18) Updates to AOD output.  Also have the option to write to HDF 
!        (amv, bmy, 12/21/09)
!  (19) Added MEGAN species (mpb, bmy, 12/21/09)
!  (20) Modify AOD output to wavelength specified in jv_spec_aod.dat 
!       (clh, 05/07/10)
!  12 Nov 2010 - R. Yantosca - Now save out PEDGE-$ (pressure at level edges)
!                              rather than Psurface - PTOP
!  03 Feb 2011 - S. Kim      - Now do not scale the AOD output
!                              (recalculated in RDAER AND DUST_MOD)
!  01 Mar 2012 - R. Yantosca - Use updated GET_LOCALTIME from time_mod.F
!  06 Aug 2012 - R. Yantosca - Now make IU_ND51b a local module variable
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !PRIVATE TYPES:
!
      !=================================================================
      ! MODULE VARIABLES
      !
      ! GOOD             : Array denoting grid boxes w/in LT limits
      ! GOOD_CT          : # of "good" times per grid box
      ! GOOD_CT_CHEM     : # of "good" chemistry timesteps
      ! COUNT_CHEM3D     : Counter for 3D chemistry boxes
      ! ND51_HR_WRITE    : Hour at which to save to disk
      ! I0               : Offset between global & nested grid
      ! J0               : Offset between global & nested grid
      ! IOFF             : Longitude offset
      ! JOFF             : Latitude offset
      ! LOFF             : Altitude offset
      ! ND51_HR1         : Starting hour of user-defined LT interval
      ! ND51_HR2         : Ending hour of user-defined LT interval
      ! ND51_IMIN        : Minimum latitude  index for DIAG51 region
      ! ND51_IMAX        : Maximum latitude  index for DIAG51 region
      ! ND51_JMIN        : Minimum longitude index for DIAG51 region
      ! ND51_JMAX        : Maximum longitude index for DIAG51 region
      ! ND51_LMIN        : Minimum altitude  index for DIAG51 region
      ! ND51_LMAX        : Minimum latitude  index for DIAG51 region
      ! ND51_NI          : Number of longitudes in DIAG51 region 
      ! ND51_NJ          : Number of latitudes  in DIAG51 region
      ! ND51_NL          : Number of levels     in DIAG51 region
      ! ND51_N_TRACERS   : Number of tracers for DIAG51
      ! ND51_OUTPUT_FILE : Name of bpch file w  timeseries data
      ! ND51_TRACERS     : Array of DIAG51 tracer numbers
      ! Q                : Accumulator array for various quantities
      ! TAU0             : Starting TAU used to index the bpch file
      ! TAU1             : Ending TAU used to index the bpch file
      ! HALFPOLAR        : Used for bpch file output
      ! CENTER180        : Used for bpch file output
      ! LONRES           : Used for bpch file output
      ! LATRES           : Used for bpch file output
      ! MODELNAME        : Used for bpch file output
      ! RESERVED         : Used for bpch file output
      !=================================================================
      
      ! Scalars
      INTEGER              :: IOFF,           JOFF,    LOFF
      INTEGER              :: I0,             J0
      ! Increased to 120 from 100 (mpb,2009)
      ! Increased to 121 (ccc, 4/20/10)
      INTEGER              :: ND51_N_TRACERS, ND51_TRACERS(121)
      INTEGER              :: ND51_IMIN,      ND51_IMAX
      INTEGER              :: ND51_JMIN,      ND51_JMAX
      INTEGER              :: ND51_LMIN,      ND51_LMAX
      INTEGER              :: ND51_FREQ,      ND51_NI
      INTEGER              :: ND51_NJ,        ND51_NL
      INTEGER              :: HALFPOLAR
      INTEGER, PARAMETER   :: CENTER180=1
      REAL*4               :: LONRES,         LATRES
      REAL*8               :: TAU0,           TAU1
      REAL*8               :: ND51_HR1,       ND51_HR2
      REAL*8               :: ND51_HR_WRITE
      CHARACTER(LEN=20)    :: MODELNAME
      CHARACTER(LEN=40)    :: RESERVED = ''
      CHARACTER(LEN=80)    :: TITLE
      CHARACTER(LEN=255)   :: ND51_OUTPUT_FILE

      ! LUN for ND51b output file
      INTEGER              :: IU_ND51b

      ! Arrays
      INTEGER, ALLOCATABLE :: GOOD(:)
      INTEGER, ALLOCATABLE :: GOOD_CT(:)
      ! For chemistry
      INTEGER, ALLOCATABLE :: GOOD_CHEM(:)
      INTEGER, ALLOCATABLE :: GOOD_CT_CHEM(:)
      INTEGER, ALLOCATABLE :: COUNT_CHEM3D(:,:,:)
      ! For emissions
      INTEGER, ALLOCATABLE :: GOOD_EMIS(:)
      INTEGER, ALLOCATABLE :: GOOD_CT_EMIS(:)
      REAL*8,  ALLOCATABLE :: Q(:,:,:,:)

      !=================================================================
      ! Original code from old DIAG51_MOD.  Leave here as a guide to 
      ! figure out when the averaging periods should be and when to
      ! write to disk (bmy, 9/28/04)
      !
      !! For timeseries between 1300 and 1700 LT, uncomment this code:
      !!
      !! Need to write to the bpch file at 12 GMT, since this covers
      !! an entire day over the US grid (amf, bmy, 12/1/00)
      !!
      !INTEGER, PARAMETER   :: NHMS_WRITE = 120000
      !REAL*8,  PARAMETER   :: HR1        = 13d0
      !REAL*8,  PARAMETER   :: HR2        = 17d0
      !CHARACTER(LEN=255)   :: FILENAME   = 'ts1_4pm.bpch'
      !=================================================================
      ! For timeseries between 1000 and 1200 LT, uncomment this code:
      !
      ! Between 10 and 12 has been chosen because the off-polar orbit 
      ! of GOME traverses (westward) through local times between 12 
      ! and 10 over North America, finally crossing the equator at 
      ! 10.30 (local time).
      !
      ! Need to write to the bpch file at 00 GMT, since we will be 
      ! interested in the whole northern hemisphere (pip, 12/1/00)
      !
      !INTEGER, PARAMETER   :: NHMS_WRITE = 000000
      !REAL*8,  PARAMETER   :: HR1        = 10d0
      !REAL*8,  PARAMETER   :: HR2        = 12d0
      !CHARACTER(LEN=255)   :: FILENAME   ='ts10_12pm.bpch'
      !=================================================================

      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement 
      !=================================================================
      CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: diag51b
!
! !DESCRIPTION:  Subroutine DIAG51 generates time series (averages from !
!  10am - 12pm LT or 1pm - 4pm LT) for the US grid area.  Output is to 
!  binary punch files or HDF5 files.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE DIAG51b( am_I_Root, Input_Opt,
     &                    State_Met, State_Chm, RC )
!
! !USES:
!
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE GIGC_State_Chm_Mod, ONLY : ChmState
      USE GIGC_State_Met_Mod, ONLY : MetState
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)    :: am_I_Root   ! Are we on the root CPU?
      TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
      TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology State object
!
! !INPUT/OUTPUT PARAMETERS:
!
      TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!
! !REVISION HISTORY: 
!  20 Jul 2004 - R. Yantosca - Initial version
!  (1 ) Rewritten for clarity (bmy, 7/20/04)
!  (2 ) Added TAU_W as a local variable (bmy, 9/28/04)
!  02 Dec 2010 - R. Yantosca - Added ProTeX headers
!  25 Mar 2013 - R. Yantosca - Now accept am_I_Root, Input_Opt, State_Chm, RC
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL*8 :: TAU_W

      !=================================================================
      ! DIAG51 begins here!
      !=================================================================
      
      ! Construct array of where local times are between HR1, HR2
      CALL GET_LOCAL_TIME

      ! Accumulate data in the Q array
      CALL ACCUMULATE_DIAG51
     &   ( am_I_Root, Input_Opt, State_Met, State_Chm, RC )

      ! Write data to disk at the proper time
      IF ( ITS_TIME_FOR_WRITE_DIAG51( TAU_W ) .and. 
     &     Input_Opt%DO_DIAG_WRITE ) THEN
         CALL WRITE_DIAG51( am_I_Root, Input_Opt, TAU_W, RC )
      ENDIF

      END SUBROUTINE DIAG51b
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_local_time
!
! !DESCRIPTION: Subroutine GET\_LOCAL\_TIME computes the local time and 
!  returns an array of points where the local time is between two user-defined 
!  limits. 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE GET_LOCAL_TIME
!
! !USES:
!
      USE TIME_MOD, ONLY : GET_LOCALTIME
      USE TIME_MOD, ONLY : GET_TS_DYN

      USE CMN_SIZE_MOD   ! Size parameters
! 
! !REMARKS:
!  For now use GET_LOCALTIME( I, 1, 1 ) which will be independent of J and L
!  for a pure cartesian grid.  This may need to be revisited once G-C is
!  interfaced into a GCM.
!
! !REVISION HISTORY: 
!  20 Jul 2004 - R. Yantosca - Initial version
!  (1 ) The 1d-3 in the computation of XLOCTM is to remove roundoff ambiguity 
!        if a the local time should fall exactly on an hour boundary.
!        (bmy, 11/29/00)
!  (2 ) Bug fix: XMID(I) should be XMID(II).  Also updated comments.
!        (bmy, 7/6/01)
!  (3 ) Updated comments (rvm, bmy, 2/27/02)
!  (4 ) Now uses function GET_LOCALTIME of "time_mod.f" (bmy, 3/27/03) 
!  (5 ) Removed reference to CMN (bmy, 7/20/04)
!  (6 ) Bug fix: LT should be REAL*8 and not INTEGER (ccarouge, 12/10/08)
!  (7 ) We need to substract TS_DYN to the time to get the local time at 
!        the beginning of previous time step. (ccc, 8/11/09)
!  02 Dec 2010 - R. Yantosca - Added ProTeX headers
!  01 Mar 2012 - R. Yantosca - Now use GET_LOCALTIME(I,J,L) from time_mod.F90
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER :: I
      REAL*8  :: LT, TS_DYN

      !=================================================================
      ! GET_LOCAL_TIME begins here!
      !=================================================================
      TS_DYN = GET_TS_DYN()
      TS_DYN = TS_DYN / 60d0

      DO I = 1, IIPAR

         ! Get local time
         LT = GET_LOCALTIME( I, 1, 1 ) - TS_DYN
         IF ( LT < 0  ) LT = LT + 24d0

         ! GOOD indicates which boxes have local times between HR1 and HR2
         IF ( LT >= ND51_HR1 .and. LT <= ND51_HR2 ) THEN
            GOOD(I) = 1
         ENDIF
      ENDDO

      END SUBROUTINE GET_LOCAL_TIME
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: accumulate_diag51
!
! !DESCRIPTION: Subroutine ACCUMULATE\_DIAG51 accumulates tracers into the 
!  Q array. 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE ACCUMULATE_DIAG51( am_I_Root, Input_Opt,
     &                              State_Met, State_Chm, RC )
!
! !USES:
!
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE GIGC_State_Chm_Mod, ONLY : ChmState
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE MODIS_LAI_MOD,      ONLY : ISOLAI => GC_LAI
      USE PBL_MIX_MOD,        ONLY : GET_PBL_TOP_L,    GET_PBL_TOP_m
      USE PRESSURE_MOD,       ONLY : GET_PEDGE
      USE TIME_MOD,           ONLY : GET_ELAPSED_MIN,  GET_TS_CHEM 
      USE TIME_MOD,           ONLY : TIMESTAMP_STRING, GET_TS_DYN
      USE TIME_MOD,           ONLY : GET_TS_DIAG,      GET_TS_EMIS
      USE TRACER_MOD,         ONLY : XNUMOLAIR
      USE TRACERID_MOD,       ONLY : IDTHNO3, IDTHNO4, IDTN2O5, IDTNO
      USE TRACERID_MOD,       ONLY : IDTPAN,  IDTPMN,  IDTPPN,  IDTO3
      USE TRACERID_MOD,       ONLY : IDTR4N2, IDTSALA, IDTSALC, IDTNO2
      USE TROPOPAUSE_MOD,     ONLY : ITS_IN_THE_TROP
                              
      USE CMN_FJ_MOD,         ONLY : JPMAX, JPPJ      ! includes F77_CMN_SIZE
      USE JV_CMN_MOD               ! ODAER,   QAA,    QAA_AOD
      USE CMN_O3_MOD               ! SAVEOH
      USE CMN_GCTM_MOD             ! SCALE_HEIGHT
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)    :: am_I_Root   ! Are we on the root CPU?
      TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
      TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology State object
!
! !INPUT/OUTPUT PARAMETERS:
!
      TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
! 
! !REVISION HISTORY: 
!  20 Jul 2004 - R. Yantosca - Initial version
!  (1 ) Rewrote to remove hardwiring and for better efficiency.  Added extra
!        diagnostics and updated numbering scheme.  Now scale optical depths
!        to 400 nm (which is usually what QAA(2,*) is.  (bmy, 7/20/04) 
!  (2 ) Now reference GET_ELAPSED_MIN and GET_TS_CHEM from "time_mod.f".  
!        Also now all diagnostic counters are 1-D since they only depend on 
!        longitude. Now only archive NO, NO2, OH, O3 on every chemistry 
!        timestep (i.e. only when fullchem is called). (bmy, 10/25/04)
!  (3 ) Only archive AOD's when it is a chem timestep (bmy, 1/14/05)
!  (4 ) Remove reference to "CMN".  Also now get PBL heights in meters and 
!        model layers from GET_PBL_TOP_m and GET_PBL_TOP_L of "pbl_mix_mod.f".
!        (bmy, 2/16/05)
!  (5 ) Now reference CLDF and BXHEIGHT from "dao_mod.f".  Now save 3-D cloud 
!        fraction as tracer #79 and box height as tracer #93.  Now remove 
!        references to CLMOSW, CLROSW, and PBL from "dao_mod.f". (bmy, 4/20/05)
!  (6 ) Remove TRCOFFSET since it's always zero  Also now get HALFPOLAR for
!        both GCAP and GEOS grids.  (bmy, 6/28/05)
!  (7 ) Now do not save SLP data if it is not allocated (bmy, 8/2/05)
!  (8 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (9 ) Now references XNUMOLAIR from "tracer_mod.f" (bmy, 10/25/05)
!  (10) Now account for time spent in the trop for non-tracers (phs, 1/24/07)
!  (11) We determine points corresponding to the time window at each timestep.
!       But accumulate only when it's time for diagnostic (longest t.s.)
!       (ccc, 8/12/09)
!  (12) Add outputs ("DAO-FLDS" and "BIOGSRCE" categories). Add GOOD_EMIS and 
!       GOOD_CT_EMIS to manage emission outputs. (ccc, 11/20/09)
!  (13) Output AOD at 3rd jv_spec.dat row wavelength.  Include all seven dust 
!        bin's individual AOD (amv, bmy, 12/21/09)
!  (12) Added MEGAN species (mpb, bmy, 12/21/09)
!  12 Nov 2010 - R. Yantosca - Now save out PEDGE-$ (pressure at level edges)
!                              rather than Psurface - PTOP
!  11 Apr 2012 - R. Yantosca - Replace lai_mod.F with modis_lai_mod.F
!  09 Nov 2012 - M. Payer    - Replaced all met field arrays with State_Met
!                              derived type object
!  14 Mar 2013 - M. Payer    - Replace NOx and Ox with NO, NO2, and O3 as part
!                              of removal of NOx-Ox partitioning
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL, SAVE     :: FIRST = .TRUE.
      LOGICAL, SAVE     :: IS_FULLCHEM, IS_SEASALT
      LOGICAL, SAVE     :: IS_NO,       IS_O3,     IS_NO2
      LOGICAL, SAVE     :: IS_CLDTOPS,  IS_NOy,    IS_OPTD, IS_SLP
      LOGICAL           :: IS_CHEM,     IS_DIAG,   IS_EMIS
      INTEGER           :: H, I, J, K, L, M, N
      INTEGER           :: PBLINT,  R, X, Y, W, XMIN
      REAL*8            :: C1, C2, PBLDEC, TEMPBL, TMP, SCALEAODnm
      CHARACTER(LEN=16) :: STAMP

      ! Aerosol types (rvm, aad, bmy, 7/20/04)
      INTEGER           :: IND(6) = (/ 22, 29, 36, 43, 50, 15 /)

      ! For fields from Input_Opt
      INTEGER           :: N_TRACERS
      REAL*8            :: TCVV(Input_Opt%N_TRACERS)

      ! Pointers
      REAL*8, POINTER   :: STT(:,:,:,:)
      REAL*8, POINTER   :: AD(:,:,:)

      !=================================================================
      ! ACCUMULATE_DIAG51 begins here!
      !=================================================================

      ! Assume success
      RC        =  GIGC_SUCCESS

      ! Copy values from Input_Opt
      N_TRACERS = Input_Opt%N_TRACERS
      TCVV      = Input_Opt%TCVV(1:N_TRACERS)

      ! Set logical flags on first call
      IF ( FIRST ) THEN
         IS_OPTD     = ASSOCIATED( State_Met%OPTD    )
         IS_CLDTOPS  = ASSOCIATED( State_Met%CLDTOPS )
         IS_SLP      = ASSOCIATED( State_Met%SLP     )
         IS_FULLCHEM = Input_Opt%ITS_A_FULLCHEM_SIM
         IS_SEASALT  = ( IDTSALA > 0 .and. IDTSALC > 0 )
         IS_NO       = ( IS_FULLCHEM .and. IDTNO   > 0 )
         IS_NO2      = ( IS_FULLCHEM .and. IDTNO2  > 0 )
         IS_O3       = ( IS_FULLCHEM .and. IDTO3   > 0 )
         IS_NOy      = ( IS_FULLCHEM .and. IDTNO   > 0 .and.
     &                   IDTNO2  > 0 .and. IDTPAN  > 0 .and.
     &                   IDTHNO3 > 0 .and. IDTPMN  > 0 .and.
     &                   IDTPPN  > 0 .and. IDTR4N2 > 0 .and.
     &                   IDTN2O5 > 0 .and. IDTHNO4 > 0 ) 
         FIRST       = .FALSE.
      ENDIF

      ! Is it a chemistry timestep?
      IS_CHEM = ( MOD( GET_ELAPSED_MIN()-GET_TS_DYN(), 
     &                 GET_TS_CHEM() ) == 0 )

      ! Is it an emissions timestep?
      IS_EMIS = ( MOD( GET_ELAPSED_MIN()-GET_TS_DYN(),
     &                 GET_TS_EMIS() ) == 0 )

      ! Is it time for diagnostic accumulation ?
      IF ( GET_ELAPSED_MIN() == 0 ) THEN
         IS_DIAG = .FALSE.
      ELSE
         IS_DIAG = ( MOD( GET_ELAPSED_MIN(), GET_TS_DIAG() ) == 0 )
      ENDIF

      ! Echo info
      STAMP = TIMESTAMP_STRING()
      WRITE( 6, 100 ) STAMP
 100  FORMAT( '     - DIAG51b: Accumulation at ', a )
      
      !=================================================================
      ! Archive tracers into accumulating array Q 
      !=================================================================

      ! Archive counter array of good points 
      IF ( IS_DIAG ) THEN
         DO X = 1, ND51_NI
            I          = GET_I( X )
            GOOD_CT(X) = GOOD_CT(X) + GOOD(I)
         ENDDO
      ENDIF

      ! Archive counter array of good points for chemistry timesteps only
      IF ( IS_CHEM ) THEN
         DO X = 1, ND51_NI
            I               = GET_I( X )
            GOOD_CT_CHEM(X) = GOOD_CT_CHEM(X) + GOOD(I)
            ! Save GOOD points only for chemistry time steps
            GOOD_CHEM(X) = GOOD(I)
         ENDDO
      ENDIF

      ! Archive counter array of good points for emissions timesteps only
      IF ( IS_EMIS ) THEN
         DO X = 1, ND51_NI
            I               = GET_I( X )
            GOOD_CT_EMIS(X) = GOOD_CT_EMIS(X) + GOOD(I)
            ! Save GOOD points only for emission time steps
            GOOD_EMIS(X) = GOOD(I)
         ENDDO
      ENDIF

      ! Also increment 3-D counter for boxes in the tropopause
      IF ( IS_FULLCHEM .and. IS_CHEM ) THEN
         
         ! Loop over levels
!$OMP PARALLEL DO 
!$OMP+DEFAULT( SHARED ) 
!$OMP+PRIVATE( X, Y, K, I, J, L )
!$OMP+SCHEDULE( DYNAMIC )
         DO K = 1, ND51_NL
            L = LOFF + K

         ! Loop over latitudes 
         DO Y = 1, ND51_NJ
            J = JOFF + Y

         ! Loop over longitudes
         DO X = 1, ND51_NI
            I = GET_I( X )

            ! Only increment if we are in the trop
            IF ( ITS_IN_THE_TROP( I, J, L, State_Met ) ) THEN
               COUNT_CHEM3D(X,Y,K) = COUNT_CHEM3D(X,Y,K) + GOOD(I)
            ENDIF

         ENDDO
         ENDDO
         ENDDO
!$OMP END PARALLEL DO
      ENDIF

      !------------------------
      ! Accumulate quantities
      !------------------------
      IF ( IS_DIAG ) THEN

         ! Initialize pointers
         STT => State_Chm%Tracers
         AD  => State_Met%AD

!$OMP PARALLEL DO 
!$OMP+DEFAULT( SHARED ) 
!$OMP+PRIVATE( W, N, X, Y, K, I, J, L, TMP, H, R, SCALEAODnm ) 
!$OMP+SCHEDULE( DYNAMIC )
         DO W = 1, ND51_N_TRACERS

         ! ND51 Tracer number
         N = ND51_TRACERS(W)

         ! Loop over levels
         DO K = 1, ND51_NL
            L = LOFF + K

         ! Loop over latitudes 
         DO Y = 1, ND51_NJ
            J = JOFF + Y

         ! Loop over longitudes
         DO X = 1, ND51_NI
            I = GET_I( X )

            ! Archive by simulation 
            IF ( N <= N_TRACERS ) THEN

               !--------------------------------------
               ! GEOS-CHEM tracers [v/v]
               !--------------------------------------

               ! Archive afternoon points
               Q(X,Y,K,W) = Q(X,Y,K,W) + 
     &                      ( STT(I,J,L,N) * TCVV(N) / 
     &                        AD(I,J,L)    * GOOD(I) )

            ELSE IF ( N == 91 .and. IS_O3 ) THEN

               !--------------------------------------
               ! Pure O3 [v/v]
               ! NOTE: Only archive at chem timestep
               !--------------------------------------

               ! Accumulate data
               Q(X,Y,K,W) = Q(X,Y,K,W)   + ( STT(I,J,L,IDTO3) *
     &                      TCVV(IDTO3)  /   AD(I,J,L) *
     &                      GOOD_CHEM(X) )

            ELSE IF ( N == 92 .and. IS_NO ) THEN

               !--------------------------------------
               ! NO [v/v]
               ! NOTE: Only archive at chem timestep
               !--------------------------------------
               
               ! Accumulate data
               Q(X,Y,K,W) = Q(X,Y,K,W)   + ( STT(I,J,L,IDTNO) *
     &                      TCVV(IDTNO)  /   AD(I,J,L) *
     &                      GOOD_CHEM(X) )

            ELSE IF ( N == 93 .and. IS_NOy ) THEN

               !--------------------------------------
               ! NOy [v/v]
               !--------------------------------------
  
               ! Temp variable for accumulation
               TMP = 0d0
            
               ! NO
               TMP = TMP + ( TCVV(IDTNO)         * GOOD(I) *
     &                       STT(I,J,L,IDTNO)    / AD(I,J,L) )

               ! NO2
               TMP = TMP + ( TCVV(IDTNO2)        * GOOD(I) *
     &                       STT(I,J,L,IDTNO2)   / AD(I,J,L) )
               ! PAN
               TMP = TMP + ( TCVV(IDTPAN)        * GOOD(I) *
     &                       STT(I,J,L,IDTPAN)   / AD(I,J,L) )

               ! HNO3
               TMP = TMP + ( TCVV(IDTHNO3)       * GOOD(I) *
     &                       STT(I,J,L,IDTHNO3)  / AD(I,J,L) )
            
               ! PMN
               TMP = TMP + ( TCVV(IDTPMN)        * GOOD(I) *
     &                       STT(I,J,L,IDTPMN)   / AD(I,J,L) )

               ! PPN
               TMP = TMP + ( TCVV(IDTPPN)        * GOOD(I) *
     &                       STT(I,J,L,IDTPPN)   / AD(I,J,L) )
 
               ! R4N2
               TMP = TMP + ( TCVV(IDTR4N2)       * GOOD(I) *
     &                       STT(I,J,L,IDTR4N2)  / AD(I,J,L) )
            
               ! N2O5
               TMP = TMP + ( 2d0 * TCVV(IDTN2O5) * GOOD(I) *
     &                       STT(I,J,L,IDTN2O5)  / AD(I,J,L) )
                        
               ! HNO4
               TMP = TMP + ( TCVV(IDTHNO4)       * GOOD(I) *
     &                       STT(I,J,L,IDTHNO4)  / AD(I,J,L) )

               ! Save afternoon points
               Q(X,Y,K,W) = Q(X,Y,K,W) + TMP
    
            ELSE IF ( N == 76 .and. IS_FULLCHEM ) THEN

               !--------------------------------------
               ! OH [molec/cm3]
               ! NOTE: Only archive at chem timestep
               !--------------------------------------

               ! Accumulate data
               Q(X,Y,K,W) = Q(X,Y,K,W) + 
     &              ( SAVEOH(I,J,L) * GOOD_CHEM(X) )
              
            ELSE IF ( N == 77 .and. IS_NO2 ) THEN

               !--------------------------------------
               ! NO2 [v/v]
               ! NOTE: Only archive at chem timestep
               !--------------------------------------     
               Q(X,Y,K,W) = Q(X,Y,K,W)   + ( STT(I,J,L,IDTNO2) *
     &                      TCVV(IDTNO2) /   AD(I,J,L) *
     &                      GOOD_CHEM(X) )
 
            ELSE IF ( N == 78 ) THEN

               !--------------------------------------
               ! PBL HEIGHTS [m] 
               !--------------------------------------
               IF ( K == 1 ) THEN
                  Q(X,Y,K,W) = Q(X,Y,K,W) + 
     &                         ( GET_PBL_TOP_m( I, J ) * GOOD(I) )  
               ENDIF

            ELSE IF ( N == 79 ) THEN

               !--------------------------------------
               ! PBL HEIGHTS [layers] 
               !--------------------------------------
               IF ( K == 1 ) THEN
                  Q(X,Y,K,W) = Q(X,Y,K,W) +
     &                         ( GET_PBL_TOP_L( I, J ) * GOOD(I) )
               ENDIF

            ELSE IF ( N == 80 ) THEN

               !--------------------------------------
               ! AIR DENSITY [molec/cm3] 
               !--------------------------------------
               Q(X,Y,K,W) = Q(X,Y,K,W) + ( State_Met%AIRDEN(L,I,J) *
     &                      XNUMOLAIR  * 1d-6 * GOOD(I) )

            ELSE IF ( N == 81 ) THEN

               !--------------------------------------
               ! 3-D CLOUD FRACTION [unitless]
               !--------------------------------------
               Q(X,Y,K,W) = Q(X,Y,K,W) + ( State_Met%CLDF(L,I,J) *
     &                      GOOD(I) )

            ELSE IF ( N == 82 .and. IS_OPTD ) THEN

               !--------------------------------------
               ! COLUMN OPTICAL DEPTH [unitless]
               !--------------------------------------
               Q(X,Y,1,W) = Q(X,Y,1,W) + ( State_Met%OPTD(L,I,J) *
     &                      GOOD(I) )

            ELSE IF ( N == 83 .and. IS_CLDTOPS ) THEN

               !--------------------------------------
               ! CLOUD TOP HEIGHTS [mb]
               !--------------------------------------
               IF ( K == 1 ) THEN
                  TMP        = GET_PEDGE( I, J, State_Met%CLDTOPS(I,J) )
                  Q(X,Y,K,W) = Q(X,Y,K,W) + ( TMP * GOOD(I) )
               ENDIF

            ELSE IF ( N == 84 ) THEN

               !--------------------------------------
               ! SULFATE AOD @ jv_spec_aod.dat wavelength [unitless]
               ! NOTE: Only archive at chem timestep
               !--------------------------------------
               DO H = 1, NRH
                  
                  ! Scaling factor to 400 nm
                  ! SCALEAODnm = QAA_AOD(IND(1)+H-1) / QAA(4,IND(1)+H-1) 
                  ! We no longer need to scale by wavelength (skim, 02/03/11)
                  SCALEAODnm = 1.0                

                  ! Accumulate
                  Q(X,Y,K,W) = Q(X,Y,K,W) + 
     &                   ( ODAER(I,J,L,H) * SCALEAODnm * GOOD_CHEM(X) )
               ENDDO

            ELSE IF ( N == 85 ) THEN

               !--------------------------------------
               ! BLACK CARBON AOD @ jv_spec_aod.dat wavelength [unitless]
               ! NOTE: Only archive at chem timestep
               !--------------------------------------
               DO R = 1, NRH

                  ! Index for ODAER
                  H          = NRH    + R

                  ! Scaling factor to 400 nm
                  ! SCALEAODnm = QAA_AOD(IND(2)+R-1) / QAA(4,IND(2)+R-1)
                  ! We no longer need to scale by wavelength (skim, 02/03/11)
                  SCALEAODnm = 1.0                

                  ! Accumulate
                  Q(X,Y,K,W) = Q(X,Y,K,W) + 
     &                   ( ODAER(I,J,L,H) * SCALEAODnm * GOOD_CHEM(X) )
               ENDDO

            ELSE IF ( N == 86 ) THEN

               !--------------------------------------
               ! ORG CARBON AOD @ jv_spec_aod.dat wavelength [unitless]
               ! NOTE: Only archive at chem timestep
               !--------------------------------------
               DO R = 1, NRH

                  ! Index for ODAER
                  H          = 2*NRH  + R

                  ! Scaling factor to 400 nm
                  ! SCALEAODnm = QAA_AOD(IND(3)+R-1) / QAA(4,IND(3)+R-1)
                  ! We no longer need to scale by wavelength (skim, 02/03/11)
                  SCALEAODnm = 1.0                

                  ! Accumulate
                  Q(X,Y,K,W) = Q(X,Y,K,W) +
     &                   ( ODAER(I,J,L,H) * SCALEAODnm * GOOD_CHEM(X) )
               ENDDO

            ELSE IF ( N == 87 ) THEN

               !--------------------------------------
               ! ACCUM SEASALT AOD @ jv_spec_aod.dat wavelength [unitless]
               ! NOTE: Only archive at chem timestep
               !--------------------------------------
               DO R = 1, NRH

                  ! Index for ODAER
                  H          = 3*NRH  + R

                  ! Scaling factor to 400 nm
                  ! SCALEAODnm = QAA_AOD(IND(4)+R-1) / QAA(4,IND(4)+R-1)
                  ! We no longer need to scale by wavelength (skim, 02/03/11)
                  SCALEAODnm = 1.0                

                  ! Accumulate
                  Q(X,Y,K,W) = Q(X,Y,K,W) + 
     &                   ( ODAER(I,J,L,H) * SCALEAODnm * GOOD_CHEM(X) ) 
               ENDDO

            ELSE IF ( N == 88 ) THEN

               !--------------------------------------
               ! COARSE SEASALT AOD @ jv_spec_aod.dat wavelength [unitless]
               ! NOTE: Only archive at chem timestep
               !--------------------------------------
               DO R = 1, NRH

                  ! Index for ODAER
                  H          = 4*NRH + R

                  ! Scaling factor to 400 nm
                  ! SCALEAODnm = QAA_AOD(IND(5)+R-1) / QAA(4,IND(5)+R-1)
                  ! We no longer need to scale by wavelength (skim, 02/03/11)
                  SCALEAODnm = 1.0                

                  ! Accumulate
                  Q(X,Y,K,W) = Q(X,Y,K,W) + 
     &                   ( ODAER(I,J,L,H) * SCALEAODnm * GOOD_CHEM(X) )
               ENDDO

            ELSE IF ( N == 89 ) THEN               

               !--------------------------------------
               ! TOTAL DUST OPTD @ jv_spec_aod.dat wavelength  [unitless]
               ! NOTE: Only archive at chem timestep
               !--------------------------------------
               DO R = 1, NDUST 

                  ! Scaling factor to 400 nm
                  ! SCALEAODnm = QAA_AOD(IND(6)+R-1) / QAA(4,IND(6)+R-1)
                  ! We no longer need to scale by wavelength (skim, 02/03/11)
                  SCALEAODnm = 1.0                

                  ! Accumulate
                  Q(X,Y,K,W) = Q(X,Y,K,W) + 
     &                 ( ODMDUST(I,J,L,R) * SCALEAODnm * GOOD_CHEM(X) )
               ENDDO

            ELSE IF ( N > 114 ) THEN

               !--------------------------------------
               ! Dust BINS 1-7 optical depth [unitless]
               ! NOTE: Only archive at chem timestep
               !--------------------------------------
               R = N - 114

               ! Scaling factor to AOD wavelength (clh, 05/09)
               ! SCALEAODnm = QAA_AOD(IND(6)+R-1) / QAA(4,IND(6)+R-1)
               ! We no longer need to scale by wavelength (skim, 02/03/11)
               SCALEAODnm = 1.0

               ! Accumulate
               Q(X,Y,K,W) = Q(X,Y,K,W) +
     &              ( ODMDUST(I,J,L,R) * SCALEAODnm * GOOD_CHEM(X) )

            ELSE IF ( N == 90 .and. IS_SEASALT ) THEN

               !-----------------------------------
               ! TOTAL SEASALT TRACER [v/v]
               !-----------------------------------
               Q(X,Y,K,W) = Q(X,Y,K,W) +
     &                      ( STT(I,J,L,IDTSALA)   + 
     &                        STT(I,J,L,IDTSALC) ) *
     &                        TCVV(IDTSALA)        /
     &                        AD(I,J,L) * GOOD(I)

            ELSE IF ( N == 94 ) THEN

               !-----------------------------------
               ! GRID BOX HEIGHTS [m]
               !-----------------------------------               
               Q(X,Y,K,W) = Q(X,Y,K,W) + ( State_Met%BXHEIGHT(I,J,L) *
     &                      GOOD(I) )

            ELSE IF ( N == 95 ) THEN

               !-----------------------------------
               ! RELATIVE HUMIDITY [%]
               !-----------------------------------               
               Q(X,Y,K,W) = Q(X,Y,K,W) + ( State_Met%RH(I,J,L) *
     &                      GOOD(I) )

            ELSE IF ( N == 96 .and. IS_SLP ) THEN

               !-----------------------------------
               ! SEA LEVEL PRESSURE [hPa]
               !-----------------------------------               
               IF ( K == 1 ) THEN
                  Q(X,Y,K,W) = Q(X,Y,K,W) + ( State_Met%SLP(I,J) *
     &                         GOOD(I) )
               ENDIF

            ELSE IF ( N == 97 ) THEN

               !-----------------------------------
               ! ZONAL (U) WIND [M/S]
               !-----------------------------------               
               Q(X,Y,K,W) = Q(X,Y,K,W) + ( State_Met%U(I,J,L) *
     &                      GOOD(I) )

            ELSE IF ( N == 98 ) THEN

               !-----------------------------------
               ! MERIDIONAL (V) WIND [M/S]
               !-----------------------------------               
               Q(X,Y,K,W) = Q(X,Y,K,W) + ( State_Met%V(I,J,L) *
     &                      GOOD(I) )

            ELSE IF ( N == 99 ) THEN

               !-----------------------------------
               ! PEDGE-$ (prs @ level edges) [hPa]
               !-----------------------------------
               Q(X,Y,K,W) = Q(X,Y,K,W) + ( GET_PEDGE(I,J,K) * GOOD(I) )

            ELSE IF ( N == 100 ) THEN 

               !-----------------------------------
               ! TEMPERATURE [K]
               !-----------------------------------
               Q(X,Y,K,W) = Q(X,Y,K,W) + ( State_Met%T(I,J,L) *
     &                      GOOD(I) )

! =====================================================================
! Added with MEGAN v2.1. (ccc, 11/20/09)

            ELSE IF ( N == 101 ) THEN 

               !-----------------------------------
               ! PAR DR [W/m2] (mpb,2009)
               !-----------------------------------
              
               IF ( K == 1 ) THEN
              
                  Q(X,Y,K,W) = Q(X,Y,K,W) +
     &                       ( State_Met%PARDR(I,J) * GOOD(I) )
               ENDIF

            ELSE IF ( N == 102 ) THEN 

               !-----------------------------------
               ! PAR DF [W/m2] (mpb,2009)
               !-----------------------------------
              
               IF ( K == 1 ) THEN
              
                  Q(X,Y,K,W) =  Q(X,Y,K,W) +
     &                       (  State_Met%PARDF(I,J) * GOOD(I) )
               ENDIF

            ELSE IF ( N == 103 ) THEN 

               !-----------------------------------
               ! DAILY LAI  [cm2/cm2] (mpb,2009)
               !-----------------------------------
              
               IF ( K == 1 ) THEN
              
                  Q(X,Y,K,W) =  Q(X,Y,K,W) +
     &                     ( ISOLAI(I,J) * GOOD(I) ) 
               ENDIF


            ELSE IF ( N == 104 ) THEN

               !-----------------------------------
               ! T at 2m [K] (mpb,2009)
               !----------------------------------- 
              
               IF ( K == 1 ) THEN

                  Q(X,Y,K,W) = Q(X,Y,K,W) +
     &                       ( State_Met%TS(I,J) * GOOD(I) )

               ENDIF

            ELSE IF ( N == 105 ) THEN

               !-----------------------------------
               ! ISOPRENE EMISSIONS [atom C/cm2/s]
               ! (mpb,2008)
               !-----------------------------------

               IF ( K == 1 ) THEN
              
                  Q(X,Y,K,W) =  Q(X,Y,K,W) +
     &                         ( EMISS_BVOC(I,J,1) * GOOD_EMIS(I) ) 

               ENDIF

            ELSE IF ( N == 106 ) THEN

               !------------------------------------
               ! MONOTERPENE EMISSIONS [atomC/cm2/s]
               ! (mpb,2008)
               !------------------------------------

               IF ( K == 1 ) THEN
              
                  Q(X,Y,K,W) =  Q(X,Y,K,W) +
     &                     ( EMISS_BVOC(I,J,2) * GOOD_EMIS(I) ) 
               ENDIF

            ELSE IF ( N == 107 ) THEN

               !-----------------------------------
               ! MBO EMISSIONS [atom C/cm2/s]
               ! (mpb,2008)
               !-----------------------------------

               IF ( K == 1 ) THEN
              
                  Q(X,Y,K,W) =  Q(X,Y,K,W) +
     &                     ( EMISS_BVOC(I,J,3) * GOOD_EMIS(I) ) 
               ENDIF

            ELSE IF ( N == 108 ) THEN

               !-----------------------------------
               ! A-PINENE EMISSIONS [atom C/cm2/s]
               ! (mpb,2009)
               !-----------------------------------

               IF ( K == 1 ) THEN
              
                  Q(X,Y,K,W) =  Q(X,Y,K,W) +
     &                     ( EMISS_BVOC(I,J,4) * GOOD_EMIS(I) ) 
               ENDIF

           ELSE IF ( N == 109 ) THEN

               !-----------------------------------
               ! B-PINENE EMISSIONS [atom C/cm2/s]
               ! (mpb,2009)
               !-----------------------------------

               IF ( K == 1 ) THEN
              
                  Q(X,Y,K,W) =  Q(X,Y,K,W) +
     &                     ( EMISS_BVOC(I,J,5) * GOOD_EMIS(I) ) 
               ENDIF


           ELSE IF ( N == 110 ) THEN

               !-----------------------------------
               ! LIMONENE EMISSIONS [atom C/cm2/s]
               ! (mpb,2009)
               !-----------------------------------

               IF ( K == 1 ) THEN
              
                  Q(X,Y,K,W) =  Q(X,Y,K,W) +
     &                     ( EMISS_BVOC(I,J,6) * GOOD_EMIS(I) ) 
               ENDIF


           ELSE IF ( N == 111 ) THEN

               !-----------------------------------
               ! SABINE EMISSIONS [atom C/cm2/s]
               ! (mpb,2009)
               !-----------------------------------

               IF ( K == 1 ) THEN
              
                  Q(X,Y,K,W) =  Q(X,Y,K,W) +
     &                     ( EMISS_BVOC(I,J,7) * GOOD_EMIS(I) ) 
               ENDIF

           ELSE IF ( N == 112 ) THEN

               !-----------------------------------
               ! MYRCENE EMISSIONS [atom C/cm2/s]
               ! (mpb,2009)
               !-----------------------------------

               IF ( K == 1 ) THEN
              
                  Q(X,Y,K,W) =  Q(X,Y,K,W) +
     &                     ( EMISS_BVOC(I,J,8) * GOOD_EMIS(I) ) 
               ENDIF

           ELSE IF ( N == 113 ) THEN

               !-----------------------------------
               ! 3-CARENE EMISSIONS [atom C/cm2/s]
               ! (mpb,2009)
               !-----------------------------------

               IF ( K == 1 ) THEN
              
                  Q(X,Y,K,W) =  Q(X,Y,K,W) +
     &                     ( EMISS_BVOC(I,J,9) * GOOD_EMIS(I) ) 
               ENDIF

           ELSE IF ( N == 114 ) THEN

               !-----------------------------------
               ! OCIMENE EMISSIONS [atom C/cm2/s]
               ! (mpb,2009)
               !-----------------------------------

               IF ( K == 1 ) THEN
              
                  Q(X,Y,K,W) =  Q(X,Y,K,W) +
     &                     ( EMISS_BVOC(I,J,10) * GOOD_EMIS(I) ) 
               ENDIF

            ENDIF
         ENDDO
         ENDDO
         ENDDO
         ENDDO 
!$OMP END PARALLEL DO
         GOOD(:) = 0

         ! Free pointers
         NULLIFY( STT )
         NULLIFY( AD  )

      ENDIF

      END SUBROUTINE ACCUMULATE_DIAG51
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: its_time_for_write_diag51
!
! !DESCRIPTION: Function ITS\_TIME\_FOR\_WRITE\_DIAG51 returns TRUE if it's 
!  time to write the ND51 bpch file to disk.  We test the time at the next 
!  dynamic timestep so that we can write to disk properly.
!\\
!\\
! !INTERFACE:
!
      FUNCTION ITS_TIME_FOR_WRITE_DIAG51( TAU_W ) RESULT( ITS_TIME )
!
! !USES:
!
      USE TIME_MOD,  ONLY : GET_HOUR
      USE TIME_MOD,  ONLY : GET_MINUTE
      USE TIME_MOD,  ONLY : GET_TAU
      USE TIME_MOD,  ONLY : GET_TAUb
      USE TIME_MOD,  ONLY : GET_TAUe
      USE TIME_MOD,  ONLY : GET_TS_DYN
      USE TIME_MOD,  ONLY : GET_TS_DIAG
      USE ERROR_MOD, ONLY : GEOS_CHEM_STOP
!
! !OUTPUT PARAMETERS:
!
      REAL*8, INTENT(OUT) :: TAU_W   ! TAU at time of disk write
! 
! !REVISION HISTORY: 
!  20 Jul 2004 - R. Yantosca - Initial version
!  (1 ) Added TAU_W so to make sure the timestamp is accurate. (bmy, 9/28/04)
!  (2 ) Add check with TS_DIAG. (ccc, 7/21/09)
!  02 Dec 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL :: ITS_TIME
      REAL*8  :: TAU, HOUR, DYN, TS_DIAG

      !=================================================================
      ! ITS_TIME_FOR_WRITE_DIAG51 begins here!
      !=================================================================

      ! Initialize
      ITS_TIME = .FALSE.

      ! Add a check for the time to save. Must be a multiple of TS_DIAG
      ! (ccc, 7/21/09)
      TS_DIAG = ( GET_TS_DIAG() / 60d0 )
      IF ( MOD(ND51_HR_WRITE, TS_DIAG) /= 0 ) THEN
         WRITE( 6, 100 ) 'ND51', ND51_HR_WRITE, TS_DIAG
 100     FORMAT( 'The ',a,' output frequency must be a multiple '
     &        'of the largest time step:', i5, i5 )
         CALL GEOS_CHEM_STOP
      ENDIF

      ! Current TAU, Hour, and Dynamic Timestep [hrs]
      TAU      = GET_TAU()
      HOUR     = ( GET_MINUTE() / 60d0 ) + GET_HOUR()
      DYN      = ( GET_TS_DYN() / 60d0 )

      ! If first timestep, return FALSE
      IF ( TAU == GET_TAUb() ) RETURN

      ! If the next dyn timestep is the hour of day
      ! when we have to save to disk, return TRUE
      IF ( MOD( HOUR, 24d0 ) == ND51_HR_WRITE ) THEN
         ITS_TIME = .TRUE.
         TAU_W    = TAU + DYN
         RETURN
      ENDIF

      ! If the next dyn timestep is the 
      ! end of the run, return TRUE
      IF ( TAU == GET_TAUe() ) THEN
         ITS_TIME = .TRUE.
         TAU_W    = TAU + DYN
         RETURN
      ENDIF

      END FUNCTION ITS_TIME_FOR_WRITE_DIAG51
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: write_diag51
!
! !DESCRIPTION: Subroutine WRITE\_DIAG51 computes the time-average of 
!  quantities between local time limits ND51\_HR1 and ND51\_HR2 and writes 
!  them to a bpch file or HDF5 file.  Arrays and counters are also zeroed 
!  for the next diagnostic interval.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE WRITE_DIAG51( am_I_Root, Input_Opt, TAU_W, RC )
!
! !USES:
!
      USE BPCH2_MOD,          ONLY : BPCH2
      USE BPCH2_MOD,          ONLY : OPEN_BPCH2_FOR_WRITE
      USE ERROR_MOD,          ONLY : ALLOC_ERR
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE inquireMod,         ONLY : findFreeLUN
      USE TIME_MOD,           ONLY : EXPAND_DATE
      USE TIME_MOD,           ONLY : GET_NYMD_DIAG    
      USE TIME_MOD,           ONLY : GET_NHMS
      USE TIME_MOD,           ONLY : GET_TAU 
      USE TIME_MOD,           ONLY : TIMESTAMP_STRING
      USE TIME_MOD,           ONLY : GET_TS_DYN

#if   defined( USE_HDF5 )
      ! Only include this if we are linking to HDF5 library (bmy, 12/21/09)
      USE HDF_MOD,            ONLY : OPEN_HDF
      USE HDF_MOD,            ONLY : CLOSE_HDF
      USE HDF_MOD,            ONLY : WRITE_HDF
      USE HDF5,               ONLY : HID_T
      INTEGER(HID_T)              :: IU_ND51b_HDF
#endif

      USE CMN_SIZE_MOD   ! Size Parameters
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)  :: am_I_Root   ! Are we on the root CPU?
      TYPE(OptInput), INTENT(IN)  :: Input_Opt   ! Input Options object
      REAL*8,         INTENT(IN)  :: TAU_W       ! TAU value at time of write
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT) :: RC          ! Success or failure?
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) TAU_W (REAL*8) : TAU value at time of writing to disk 
!
!  NOTES:
! !REVISION HISTORY: 
!  20 Jul 2004 - R. Yantosca - Initial version
!  (1 ) Rewrote to` remove hardwiring and for better efficiency.  Added extra
!        diagnostics and updated numbering scheme. (bmy, 7/20/04) 
!  (2 ) Added TAU_W to the arg list.  Now use TAU_W to set TAU0 and TAU0.
!        Also now all diagnostic counters are 1-D since they only depend on 
!        longitude.  Now only archive NO, NO2, OH, O3 on every chemistry
!        timestep (i.e. only when fullchem is called).  Also remove reference
!        to FIRST. (bmy, 10/25/04)
!  (3 ) Now divide tracers 82-87 (i.e. various AOD's) by GOOD_CT_CHEM since
!        these are only updated once per chemistry timestep (bmy, 1/14/05)
!  (4 ) Now save grid box heights as tracer #93.  Now save 3-D cloud fraction 
!        as tracer #79 (bmy, 4/20/05)
!  (5 ) Remove references to TRCOFFSET because it's always zero (bmy, 6/24/05)
!  (6 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (7 ) DIVISOR is now a 3-D array.  Now zero COUNT_CHEM3D.  Now use CASE
!        statement instead of IF statements.  Now zero counter arrays with
!        array broadcast assignments. (phs, 1/24/07)
!  (8 ) RH should be tracer #17 under "TIME-SER" category (bmy, 2/11/08)
!  (9 ) Bug fix: replace "PS-PTOP" with "PEDGE-$" (bmy, phs, 10/7/08)
!  (10) Change timestamp used for filename.  Now save SLP under tracer #18 in 
!       "DAO-FLDS". (ccc, tai, bmy, 10/13/09)
!  (11) Now have the option of saving out to HDF5 format.  NOTE: we have to
!        bracket HDF-specific code with an #ifdef statement to avoid problems
!        if the HDF5 libraries are not installed. (amv, bmy, 12/21/09)
!  (12) Add outputs ("DAO-FLDS" and "BIOGSRCE" categories). Add GOOD_EMIS and 
!       GOOD_CT_EMIS to manage emission outputs. (ccc, 11/20/09)
!  (13) Added MEGAN species (mpb, bmy, 12/21/09)
!  12 Nov 2010 - R. Yantosca - Now save out PEDGE-$ (pressure at level edges)
!                              rather than Psurface - PTOP
!  03 Aug 2012 - R. Yantosca - Move calls to findFreeLUN out of DEVEL block
!  07 Aug 2012 - R. Yantosca - Now print LUN used to open file
!  25 Mar 2013 - R. Yantosca - Now accept am_I_Root, Input_Opt, RC
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER             :: I,   J,  L,  W, N, GMNL, GMTRC
      INTEGER             :: IOS, X, Y, K, NHMS
      CHARACTER(LEN=16)   :: STAMP
      CHARACTER(LEN=40)   :: CATEGORY
      CHARACTER(LEN=40)   :: UNIT 
      CHARACTER(LEN=255)  :: FILENAME

      ! For fields from Input_Opt
      LOGICAL             :: LND51b_HDF
      INTEGER             :: N_TRACERS

      !=================================================================
      ! WRITE_DIAG51 begins here!
      !=================================================================

      ! Assume success
      RC         =  GIGC_SUCCESS

      ! Copy values from Input_Opt
      LND51b_HDF = Input_Opt%LND51b_HDF
      N_TRACERS  = Input_Opt%N_TRACERS

      ! Find a free file LUN
      IU_ND51b = findFreeLUN()

      ! Replace date tokens in FILENAME
      FILENAME = ND51_OUTPUT_FILE

      ! Change to get the good timestamp: day that was run and not next 
      ! day if saved at midnight
      NHMS = GET_NHMS()
      IF ( NHMS == 0 ) NHMS = 240000

      CALL EXPAND_DATE( FILENAME, GET_NYMD_DIAG(), NHMS )
      
      ! Echo info
      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( '     - DIAG51b: Opening file ', a, ' on unit ', i4 ) 

      ! Open output file
      IF ( LND51b_HDF ) THEN
#if   defined( USE_HDF5 )
         ! Only include this if we are linking to HDF5 library (bmy, 12/21/09)
         CALL OPEN_HDF( IU_ND51b_HDF, FILENAME,  ND51_IMAX, ND51_IMIN, 
     &                  ND51_JMAX,    ND51_JMIN, ND51_NI,   ND51_NJ )
#endif
      ELSE
         CALL OPEN_BPCH2_FOR_WRITE( IU_ND51b, FILENAME, TITLE )
      ENDIF

      ! Set ENDING TAU for this bpch write 
      TAU1 = TAU_W
    
      !=================================================================
      ! Compute time-average of tracers between local time limits
      !=================================================================

      ! Echo info
      STAMP = TIMESTAMP_STRING()
      WRITE( 6, 110 ) STAMP
 110  FORMAT( '     - DIAG51b: Saving to disk at ', a ) 

!$OMP PARALLEL DO 
!$OMP+DEFAULT( SHARED ) 
!$OMP+PRIVATE( X, Y, K, W )

      DO W = 1, ND51_N_TRACERS
 
         ! Loop over grid boxes
         DO K = 1, ND51_NL
         DO Y = 1, ND51_NJ
         DO X = 1, ND51_NI

         SELECT CASE( ND51_TRACERS(W) )

            CASE( 91, 92, 76, 77 )
               !--------------------------------------------------------
               ! Avoid div by zero for tracers which are archived each
               ! chem timestep and only available in the troposphere
               !--------------------------------------------------------
               IF ( COUNT_CHEM3D(X,Y,K) > 0 ) THEN
                  Q(X,Y,K,W) = Q(X,Y,K,W) / COUNT_CHEM3D(X,Y,K)
               ELSE
                  Q(X,Y,K,W) = 0d0
               ENDIF

            CASE( 84:89, 115:121 )

               !--------------------------------------------------------
               ! Avoid division by zero for tracers which are archived
               ! on each chem timestep (at trop & strat levels)
               !--------------------------------------------------------
               IF ( GOOD_CT_CHEM(X) > 0 ) THEN
                  Q(X,Y,K,W) = Q(X,Y,K,W) / GOOD_CT_CHEM(X)
               ELSE
                  Q(X,Y,K,W) = 0d0
               ENDIF

            CASE( 105:114 )

               !--------------------------------------------------------
               ! Avoid division by zero for tracers which are archived 
               ! on each EMISSION timestep (at SURFACE)      (mpb,2009)
               !--------------------------------------------------------
               IF ( GOOD_CT_EMIS(X) > 0 ) THEN
                  Q(X,Y,K,W) = Q(X,Y,K,W) / GOOD_CT_EMIS(X) 
               ELSE
                  Q(X,Y,K,W) = 0d0
               ENDIF

            CASE DEFAULT

               !--------------------------------------------------------
               ! Avoid division by zero for all other tracers
               !--------------------------------------------------------
               IF ( GOOD_CT(X) > 0 ) THEN
                  Q(X,Y,K,W) = Q(X,Y,K,W) / GOOD_CT(X) 
               ELSE
                  Q(X,Y,K,W) = 0d0
               ENDIF

            END SELECT

         ENDDO
         ENDDO
         ENDDO
      ENDDO
!$OMP END PARALLEL DO
      
      !=================================================================
      ! Write each tracer from "timeseries.dat" to the timeseries file
      !=================================================================
      DO W = 1, ND51_N_TRACERS

         ! ND51 tracer number
         N = ND51_TRACERS(W)

         ! Save by simulation
         IF ( N <= N_TRACERS ) THEN

            !---------------------
            ! GEOS-CHEM tracers
            !---------------------
            CATEGORY = 'IJ-AVG-$'
            UNIT     = ''              ! Let GAMAP pick unit
            GMNL     = ND51_NL
            GMTRC    = N

         ELSE IF ( N == 91 ) THEN

            !---------------------
            ! Pure O3
            !---------------------
            CATEGORY = 'IJ-AVG-$'
            UNIT     = ''              ! Let GAMAP pick unit
            GMNL     = ND51_NL
            GMTRC    = N_TRACERS + 1

         ELSE IF ( N == 92 ) THEN
            !---------------------
            ! Pure NO [v/v]
            !---------------------
            CATEGORY = 'TIME-SER'
            UNIT     = ''              ! Let GAMAP pick unit
            GMNL     = ND51_NL
            GMTRC    = 9

         ELSE IF ( N == 93 ) THEN
            !---------------------
            ! NOy 
            !---------------------
            CATEGORY = 'TIME-SER'
            UNIT     = ''              ! Let GAMAP pick unit
            GMNL     = ND51_NL
            GMTRC    = 3

         ELSE IF ( N == 76 ) THEN

            !---------------------
            ! OH 
            !---------------------
            CATEGORY  = 'CHEM-L=$'
            UNIT      = 'molec/cm3'
            GMNL      = ND51_NL
            GMTRC     = 1

         ELSE IF ( N == 77 ) THEN

            !---------------------
            ! NO2 
            !---------------------
            CATEGORY = 'TIME-SER'
            UNIT     = ''              ! Let GAMAP pick unit
            GMNL     = ND51_NL
            GMTRC    = 25

         ELSE IF ( N == 78 ) THEN 

            !---------------------
            ! PBL Height [m] 
            !---------------------
            CATEGORY = 'PBLDEPTH'
            UNIT     = 'm'
            GMNL     = 1
            GMTRC    = 1

         ELSE IF ( N == 79 ) THEN

            !---------------------
            ! PBL Height [levels]
            !---------------------
            CATEGORY = 'PBLDEPTH'
            UNIT     = 'levels'
            GMNL     = 1
            GMTRC    = 2

         ELSE IF ( N == 80 ) THEN

            !---------------------
            ! Air Density 
            !---------------------
            CATEGORY = 'TIME-SER'
            UNIT     = 'molec/cm3'
            GMNL     = ND51_NL
            GMTRC    = 22

         ELSE IF ( N == 81 ) THEN

            !---------------------
            ! 3-D Cloud fractions
            !---------------------
            CATEGORY = 'TIME-SER'
            UNIT     = 'unitless'
            GMNL     = ND51_NL
            GMTRC    = 19

         ELSE IF ( N == 82 ) THEN

            !---------------------
            ! Column opt depths 
            !---------------------
            CATEGORY = 'TIME-SER'
            UNIT     = 'unitless'
            GMNL     = 1
            GMTRC    = 20
            
         ELSE IF ( N == 83 ) THEN
        
            !---------------------
            ! Cloud top heights 
            !---------------------
            CATEGORY = 'TIME-SER'
            UNIT     = 'hPa'
            GMNL     = 1
            GMTRC    = 21

         ELSE IF ( N == 84 ) THEN

            !---------------------
            ! Sulfate AOD
            !---------------------            
            CATEGORY = 'OD-MAP-$'
            UNIT     = 'unitless'
            GMNL     = ND51_NL
            GMTRC    = 6

         ELSE IF ( N == 85 ) THEN

            !---------------------
            ! Black Carbon AOD
            !---------------------            
            CATEGORY = 'OD-MAP-$'
            UNIT     = 'unitless'
            GMNL     = ND51_NL
            GMTRC    = 9

         ELSE IF ( N == 86 ) THEN

            !---------------------
            ! Organic Carbon AOD
            !---------------------            
            CATEGORY = 'OD-MAP-$'
            UNIT     = 'unitless'
            GMNL     = ND51_NL
            GMTRC    = 12
            
         ELSE IF ( N == 87 ) THEN

            !---------------------
            ! SS Accum AOD
            !---------------------            
            CATEGORY = 'OD-MAP-$'
            UNIT     = 'unitless'
            GMNL     = ND51_NL
            GMTRC    = 15

         ELSE IF ( N == 88 ) THEN

            !---------------------
            ! SS Coarse AOD
            !---------------------            
            CATEGORY = 'OD-MAP-$'
            UNIT     = 'unitless'
            GMNL     = ND51_NL
            GMTRC    = 18

         ELSE IF ( N == 89 ) THEN

            !---------------------
            ! Total dust OD
            !---------------------   
            CATEGORY = 'OD-MAP-$'
            UNIT     = 'unitless'
            GMNL     = ND51_NL
            GMTRC    = 4

         ELSE IF ( N > 114 ) THEN

            !---------------------
            ! dust OD (bins 1-7)
            !---------------------
            CATEGORY = 'OD-MAP-$'
            UNIT     = 'unitless'
            GMNL     = ND51_NL
            GMTRC    = N - 94

         ELSE IF ( N == 90 ) THEN

            !---------------------
            ! Total seasalt
            !---------------------            
            CATEGORY = 'TIME-SER'
            UNIT     = ''              ! Let GAMAP pick unit
            GMNL     = ND51_NL
            GMTRC    = 24

         ELSE IF ( N == 94 ) THEN

            !---------------------
            ! Grid box heights
            !---------------------            
            CATEGORY = 'BXHGHT-$'
            UNIT     = 'm'
            GMNL     = ND51_NL
            GMTRC    = 1

         ELSE IF ( N == 95 ) THEN

            !---------------------
            ! Relative humidity 
            !---------------------            
            CATEGORY = 'TIME-SER'
            UNIT     = '%'
            GMNL     = ND51_NL
            GMTRC    = 17

         ELSE IF ( N == 96 ) THEN

            !---------------------
            ! Sea level prs
            !---------------------            
            CATEGORY = 'DAO-FLDS'
            UNIT     = 'hPa'
            GMNL     = 1
            GMTRC    = 18

         ELSE IF ( N == 97 ) THEN

            !---------------------
            ! U-wind
            !---------------------            
            CATEGORY = 'DAO-3D-$'
            UNIT     = 'm/s'
            GMNL     = ND51_NL
            GMTRC    = 1

         ELSE IF ( N == 98 ) THEN

            !---------------------
            ! V-wind
            !---------------------
            CATEGORY = 'DAO-3D-$'
            UNIT     = 'm/s'
            GMNL     = ND51_NL
            GMTRC    = 2

         ELSE IF ( N == 99 ) THEN

            !---------------------
            ! Psurface - PTOP 
            !---------------------
            CATEGORY = 'PEDGE-$'
            UNIT     = 'hPa'
            GMNL     = ND51_NL
            GMTRC    = 1

         ELSE IF ( N == 100 ) THEN

            !---------------------
            ! Temperature
            !---------------------
            CATEGORY = 'DAO-3D-$'
            UNIT     = 'K'
            GMNL     = ND51_NL
            GMTRC    = 3
            
! ================================================================
! Added with MEGAN v2.1. (ccc, 11/20/09)
         ELSE IF ( N == 101 ) THEN
            
            !-----------------------------------
            ! PARDR [W/m2] (mpb,2009)
            !-----------------------------------
            CATEGORY = 'DAO-FLDS' 
            UNIT     = 'W/m2'    
            GMNL     = ND51_NL
            GMTRC    = 20

         ELSE IF ( N == 102 ) THEN
            
            !-----------------------------------
            ! PARDF [W/m2] (mpb,2009)
            !-----------------------------------
            CATEGORY = 'DAO-FLDS' 
            UNIT     = 'W/m2'    
            GMNL     = ND51_NL
            GMTRC    = 21

         ELSE IF ( N == 103 ) THEN

            !-----------------------------------
            ! DAILY LAI [W/m2] (mpb,2009)
            !-----------------------------------
            CATEGORY = 'TIME-SER' 
            UNIT     = 'm2/m2'    
            GMNL     = ND51_NL
            GMTRC    = 32

         ELSE IF ( N == 104 ) THEN

            !---------------------
            ! T at 2m
            ! (mpb,2008)
            !---------------------  
            CATEGORY = 'DAO-FLDS'
            UNIT     = 'K'     
            GMNL     = ND51_NL
            GMTRC    = 5

         ELSE IF ( N == 105 ) THEN

            !---------------------
            ! ISOPRENE emissions 
            ! (mpb,2008)
            !---------------------            
            CATEGORY = 'BIOGSRCE'
            UNIT     = 'atomC/cm2/s'
            GMNL     = ND51_NL
            GMTRC    = 1     

         ELSE IF ( N == 106 ) THEN

            !---------------------
            ! MONOTERPENE emissions 
            ! (mpb,2008)
            !---------------------       

            CATEGORY = 'BIOGSRCE'
            UNIT     = 'atomC/cm2/s'
            GMNL     = ND51_NL
            GMTRC    = 4    

         ELSE IF ( N == 107 ) THEN

            !---------------------
            ! MBO emissions 
            ! (mpb,2008)
            !---------------------   

            CATEGORY = 'BIOGSRCE'
            UNIT     = 'atomC/cm2/s'
            GMNL     = ND51_NL
            GMTRC    = 5    


         ELSE IF ( N == 108 ) THEN

            !---------------------
            ! a-pine emissions 
            ! (mpb,2008)
            !---------------------   

            CATEGORY = 'BIOGSRCE'
            UNIT     = 'atomC/cm2/s'
            GMNL     = ND51_NL
            GMTRC    = 7    

         ELSE IF ( N == 109 ) THEN

            !---------------------
            ! b-pine emissions 
            ! (mpb,2008)
            !---------------------   

            CATEGORY = 'BIOGSRCE'
            UNIT     = 'atomC/cm2/s'
            GMNL     = ND51_NL
            GMTRC    = 8   

         ELSE IF ( N == 110 ) THEN

            !---------------------
            ! Limonene emissions 
            ! (mpb,2008)
            !---------------------   

            CATEGORY = 'BIOGSRCE'
            UNIT     = 'atomC/cm2/s'
            GMNL     = ND51_NL
            GMTRC    = 9    


         ELSE IF ( N == 111 ) THEN

            !---------------------
            ! Sabinene emissions 
            ! (mpb,2008)
            !---------------------   

            CATEGORY = 'BIOGSRCE'
            UNIT     = 'atomC/cm2/s'
            GMNL     = ND51_NL
            GMTRC    = 10    

         ELSE IF ( N == 112 ) THEN

            !---------------------
            ! Myrcene emissions 
            ! (mpb,2008)
            !---------------------   

            CATEGORY = 'BIOGSRCE'
            UNIT     = 'atomC/cm2/s'
            GMNL     = ND51_NL
            GMTRC    = 11    


         ELSE IF ( N == 113 ) THEN

            !---------------------
            ! 3-carene emissions 
            ! (mpb,2008)
            !---------------------   

            CATEGORY = 'BIOGSRCE'
            UNIT     = 'atomC/cm2/s'
            GMNL     = ND51_NL
            GMTRC    = 12    

         ELSE IF ( N == 114 ) THEN

            !---------------------
            ! Ocimene emissions 
            ! (mpb,2008)
            !---------------------   

            CATEGORY = 'BIOGSRCE'
            UNIT     = 'atomC/cm2/s'
            GMNL     = ND51_NL
            GMTRC    = 13    

! ================================================================

         ELSE

            ! Otherwise skip
            CYCLE

         ENDIF

         !------------------------
         ! Save to bpch file
         !------------------------
        IF ( LND51b_HDF ) THEN
#if   defined( USE_HDF5 )
            ! Only include this if we are linking to HDF5 library 
            ! (bmy, 12/21/09)
            CALL WRITE_HDF( IU_ND51b_HDF, N,
     &                  CATEGORY,     GMTRC,        UNIT,
     &                  TAU0,         TAU1,         RESERVED,
     &                  ND51_NI,      ND51_NJ,      GMNL,
     &                  ND51_IMIN+I0, ND51_JMIN+J0, ND51_LMIN,
     &                  REAL( Q(1:ND51_NI, 1:ND51_NJ, 1:GMNL, W)))
#else
         ELSE
            CALL BPCH2( IU_ND51b,     MODELNAME,    LONRES,   
     &                  LATRES,       HALFPOLAR,    CENTER180, 
     &                  CATEGORY,     GMTRC,        UNIT,      
     &                  TAU0,         TAU1,         RESERVED,  
     &                  ND51_NI,      ND51_NJ,      GMNL,     
     &                  ND51_IMIN+I0, ND51_JMIN+J0, ND51_LMIN, 
     &                  REAL( Q(1:ND51_NI, 1:ND51_NJ, 1:GMNL, W),4))
#endif
         ENDIF
      ENDDO

      ! Echo info
      WRITE( 6, 120 ) TRIM( FILENAME )
 120  FORMAT( '     - DIAG51b: Closing file ', a )

      ! Close file
      IF ( LND51b_HDF ) THEN
#if   defined( USE_HDF5 )
         ! Only include this if we are linking to HDF5 library (bmy, 12/21/09)
         CALL CLOSE_HDF( IU_ND51b_HDF )
#endif
      ELSE
         CLOSE( IU_ND51b )
      ENDIF

      !=================================================================
      ! Re-initialize quantities for next diagnostic cycle
      !=================================================================

      ! Echo info
      STAMP = TIMESTAMP_STRING()
      WRITE( 6, 130 ) STAMP
 130  FORMAT( '     - DIAG51b: Zeroing arrays at ', a )

      ! Set STARTING TAU for the next bpch write
      TAU0 = TAU_W

      ! Zero accumulating array for tracer
      Q            = 0d0

      ! Zero counter arrays
      COUNT_CHEM3D = 0d0
      GOOD_CT      = 0d0
      GOOD_CHEM    = 0d0
      GOOD_CT_CHEM = 0d0
      GOOD_EMIS    = 0d0
      GOOD_CT_EMIS = 0d0

      END SUBROUTINE WRITE_DIAG51
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_i
!
! !DESCRIPTION: Function GET\_I returns the absolute longitude index (I), 
!  given the relative longitude index (X).
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_I( X ) RESULT( I )
!
! !USES:
!
      USE CMN_SIZE_MOD   ! Size parameters
!
! !INPUT PARAMETERS: 
!
      INTEGER, INTENT(IN) :: X   ! Relative longitude index
!
! !RETURN VALUE:
!
      INTEGER             :: I   ! Absolute longitude index
!
! !REVISION HISTORY: 
!  20 Jul 2004 - R. Yantosca - Initial version
!  02 Dec 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
      !=================================================================
      ! GET_I begins here!
      !=================================================================

      ! Add the offset to X to get I  
      I = IOFF + X

      ! Handle wrapping around the date line, if necessary
      IF ( I > IIPAR ) I = I - IIPAR

      END FUNCTION GET_I
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_diag51
!
! !DESCRIPTION: Subroutine INIT\_DIAG51b allocates and zeroes all module 
!  arrays.  It also gets values for module variables from "input\_mod.f".
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE INIT_DIAG51b( DO_ND51, N_ND51, TRACERS, HR_WRITE, 
     &                         HR1,     HR2,    IMIN,    IMAX,   
     &                         JMIN,    JMAX,   LMIN,    LMAX,  FILE )
!
! !USES:
!
      USE BPCH2_MOD,  ONLY : GET_MODELNAME
      USE BPCH2_MOD,  ONLY : GET_HALFPOLAR
      USE ERROR_MOD,  ONLY : ALLOC_ERR
      USE ERROR_MOD,  ONLY : ERROR_STOP
      USE GRID_MOD,   ONLY : GET_XOFFSET
      USE GRID_MOD,   ONLY : GET_YOFFSET
      USE GRID_MOD,   ONLY : ITS_A_NESTED_GRID
      USE TIME_MOD,   ONLY : GET_TAUb
      USE TRACER_MOD, ONLY : N_TRACERS
  
      USE CMN_SIZE_MOD   ! Size parameters
!
! !INPUT PARAMETERS: 
!
      ! DO_ND51 : Switch to turn on ND51 timeseries diagnostic
      ! N_ND51  : Number of ND51 read by "input_mod.f"
      ! TRACERS : Array w/ ND51 tracer #'s read by "input_mod.f"
      ! HR_WRITE: GMT hour of day at which to write bpch file
      ! HR1     : Lower limit of local time averaging bin
      ! HR2     : Upper limit of local time averaging bin
      ! IMIN    : Min longitude index read by "input_mod.f"
      ! IMAX    : Max longitude index read by "input_mod.f" 
      ! JMIN    : Min latitude index read by "input_mod.f" 
      ! JMAX    : Min latitude index read by "input_mod.f" 
      ! LMIN    : Min level index read by "input_mod.f" 
      ! LMAX    : Min level index read by "input_mod.f" 
      ! FILE    : ND51 output file name read by "input_mod.f"
      LOGICAL,            INTENT(IN) :: DO_ND51
      INTEGER,            INTENT(IN) :: N_ND51, TRACERS(100)
      INTEGER,            INTENT(IN) :: IMIN,   IMAX 
      INTEGER,            INTENT(IN) :: JMIN,   JMAX      
      INTEGER,            INTENT(IN) :: LMIN,   LMAX 
      REAL*8,             INTENT(IN) :: HR1,    HR2
      REAL*8,             INTENT(IN) :: HR_WRITE
      CHARACTER(LEN=255), INTENT(IN) :: FILE
!
! !REVISION HISTORY: 
!  20 Jul 2004 - R. Yantosca - Initial version
!  (1 ) Diagnostic counter arrays are now only 1-D.  Also add GOOD_CT_CHEM
!        which is the counter array of "good" boxes at each chemistry
!        timesteps.  Now allocate GOOD_CT_CHEM. (bmy, 10/25/04)
!  (2 ) Now get I0 and J0 correctly for nested grid simulations (bmy, 11/9/04)
!  (3 ) Now call GET_HALFPOLAR from "bpch2_mod.f" to get the HALFPOLAR flag 
!        value for GEOS or GCAP grids. (bmy, 6/28/05)
!  (4 ) Now allow ND51_IMIN to be equal to ND51_IMAX and ND51_JMIN to be
!        equal to ND51_JMAX.  This will allow us to save out longitude or
!        latitude transects.  Allocate COUNT_CHEM3D. (cdh, bmy, phs, 1/24/07)
!  (5 ) Allocate GOOD_EMIS and GOOD_CT_EMIS (ccc, 12/12/09)
!  02 Dec 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER            :: AS
      CHARACTER(LEN=255) :: LOCATION
      
      !=================================================================
      ! INIT_DIAG51 begins here!
      !=================================================================

      ! Initialize
      LOCATION               = 'INIT_DIAG51b ("diag51_mod.f")'
      ND51_TRACERS(:)        = 0

      ! Get values from "input_mod.f"
      DO_SAVE_DIAG51b        = DO_ND51 
      ND51_N_TRACERS         = N_ND51
      ND51_TRACERS(1:N_ND51) = TRACERS(1:N_ND51)
      ND51_HR_WRITE          = HR_WRITE
      ND51_HR1               = HR1
      ND51_HR2               = HR2
      ND51_IMIN              = IMIN
      ND51_IMAX              = IMAX
      ND51_JMIN              = JMIN
      ND51_JMAX              = JMAX
      ND51_LMIN              = LMIN
      ND51_LMAX              = LMAX
      ND51_OUTPUT_FILE       = TRIM( FILE )

      ! Make sure ND51_HR_WRITE is in the range 0-23.999 hrs
      ND51_HR_WRITE = MOD( ND51_HR_WRITE, 24d0 )

      ! Exit if ND51 is turned off 
      IF ( .not. DO_SAVE_DIAG51b ) RETURN

      !=================================================================
      ! Error check longitude, latitude, altitude limits
      !=================================================================

      ! Get grid offsets
      IF ( ITS_A_NESTED_GRID() ) THEN
         I0 = GET_XOFFSET()
         J0 = GET_YOFFSET()
      ELSE
         I0 = GET_XOFFSET( GLOBAL=.TRUE. )
         J0 = GET_YOFFSET( GLOBAL=.TRUE. )
      ENDIF

      !-----------
      ! Longitude
      !-----------

      ! Error check ND51_IMIN
      IF ( ND51_IMIN+I0 < 1 .or. ND51_IMIN+I0 > IIPAR ) THEN
         CALL ERROR_STOP( 'Bad ND51_IMIN value!', LOCATION )
      ENDIF

      ! Error check ND51_IMAX
      IF ( ND51_IMAX+I0 < 1 .or. ND51_IMAX+I0 > IIPAR ) THEN
         CALL ERROR_STOP( 'Bad ND51_IMAX value!', LOCATION )
      ENDIF

      ! Compute longitude limits to write to disk
      ! Also handle wrapping around the date line
      IF ( ND51_IMAX >= ND51_IMIN ) THEN
         ND51_NI = ( ND51_IMAX - ND51_IMIN ) + 1
      ELSE 
         ND51_NI = ( IIPAR - ND51_IMIN ) + 1 + ND51_IMAX
         WRITE( 6, '(a)' ) 'We are wrapping over the date line!'
      ENDIF

      ! Make sure that ND50_NI <= IIPAR
      IF ( ND51_NI > IIPAR ) THEN
         CALL ERROR_STOP( 'Too many longitudes!', LOCATION )
      ENDIF

      !-----------
      ! Latitude
      !-----------
      
      ! Error check JMIN_AREA
      IF ( ND51_JMIN+J0 < 1 .or. ND51_JMIN+J0 > JJPAR ) THEN
         CALL ERROR_STOP( 'Bad ND51_JMIN value!', LOCATION )
      ENDIF
     
      ! Error check JMAX_AREA
      IF ( ND51_JMAX+J0 < 1 .or.ND51_JMAX+J0 > JJPAR ) THEN
         CALL ERROR_STOP( 'Bad ND51_JMAX value!', LOCATION )
      ENDIF

      ! Compute latitude limits to write to disk (bey, bmy, 3/16/99)
      IF ( ND51_JMAX >= ND51_JMIN ) THEN
         ND51_NJ = ( ND51_JMAX - ND51_JMIN ) + 1
      ELSE
         CALL ERROR_STOP( 'ND51_JMAX < ND51_JMIN!', LOCATION )
      ENDIF     
  
      !-----------
      ! Altitude
      !-----------

      ! Error check ND51_LMIN, ND51_LMAX
      IF ( ND51_LMIN < 1 .or. ND51_LMAX > LLPAR ) THEN 
         CALL ERROR_STOP( 'Bad ND51 altitude values!', LOCATION )
      ENDIF

      ! # of levels to save in ND51 timeseries
      IF ( ND51_LMAX >= ND51_LMIN ) THEN  
         ND51_NL = ( ND51_LMAX - ND51_LMIN ) + 1
      ELSE
         CALL ERROR_STOP( 'ND51_LMAX < ND51_LMIN!', LOCATION )
      ENDIF

      !-----------
      ! Offsets
      !-----------
      IOFF      = ND51_IMIN - 1
      JOFF      = ND51_JMIN - 1
      LOFF      = ND51_LMIN - 1

      !-----------
      ! For bpch
      !-----------
      TAU0      = GET_TAUb()
      TITLE     = 'GEOS-CHEM DIAG51b time series'
      LONRES    = DISIZE
      LATRES    = DJSIZE
      MODELNAME = GET_MODELNAME()
      HALFPOLAR = GET_HALFPOLAR()

      ! Reset offsets to global values for bpch write
      I0        = GET_XOFFSET( GLOBAL=.TRUE. )
      J0        = GET_YOFFSET( GLOBAL=.TRUE. ) 

      !=================================================================
      ! Allocate arrays
      !=================================================================

      ! Array denoting where LT is between HR1 and HR2
      ALLOCATE( GOOD( IIPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'GOOD' )
      GOOD = 0

      ! Counter of "good" times per day at each grid box
      ALLOCATE( GOOD_CT( ND51_NI ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'GOOD_CT' )
      GOOD_CT = 0

      ! Counter of "good" times per day at each grid box for chemistry species
      ALLOCATE( GOOD_CHEM( ND51_NI ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'GOOD_CHEM' )
      GOOD_CHEM = 0

      ! Counter of "good" times per day for each chemistry timestep
      ALLOCATE( GOOD_CT_CHEM( ND51_NI ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'GOOD_CT_CHEM' )
      GOOD_CT_CHEM = 0

      ! Array denoting where LT is between HR1 and HR2 for emissions
      ALLOCATE( GOOD_EMIS( ND51_NI ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'GOOD_EMIS' )
      GOOD_EMIS = 0

      ! Counter of "good" times per day at each grid box for emissions
      ALLOCATE( GOOD_CT_EMIS( ND51_NI ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'GOOD_CT_EMIS' )
      GOOD_CT_EMIS = 0

      ! Accumulating array
      ALLOCATE( Q( ND51_NI, ND51_NJ, ND51_NL, ND51_N_TRACERS), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'Q' )
      Q = 0d0

      ! Accumulating array
      ALLOCATE( COUNT_CHEM3D( ND51_NI, ND51_NJ, ND51_NL ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'COUNT_CHEM3D' )
      COUNT_CHEM3D = 0

      END SUBROUTINE INIT_DIAG51b
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cleanup_diag51
!
! !DESCRIPTION: Subroutine CLEANUP\_DIAG51 deallocates all module arrays. 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CLEANUP_DIAG51b
! 
! !REVISION HISTORY: 
!  20 Jul 2004 - R. Yantosca - Initial version
!  (1 ) Now deallocate GOOD_CT_CHEM (bmy, 10/25/04)
!  (2 ) Also deallocate COUNT_CHEM3D (phs, 1/24/07)
!  (5 ) Also deallocate Allocate GOOD_EMIS and GOOD_CT_EMIS (ccc, 12/12/09)
!  02 Dec 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
      !=================================================================
      ! CLEANUP_DIAG51 begins here!
      !=================================================================
      IF ( ALLOCATED( COUNT_CHEM3D ) ) DEALLOCATE( COUNT_CHEM3D )
      IF ( ALLOCATED( GOOD         ) ) DEALLOCATE( GOOD         )
      IF ( ALLOCATED( GOOD_CT      ) ) DEALLOCATE( GOOD_CT      )
      IF ( ALLOCATED( GOOD_CT_CHEM ) ) DEALLOCATE( GOOD_CT_CHEM )
      IF ( ALLOCATED( GOOD_CHEM    ) ) DEALLOCATE( GOOD_CHEM    )
      IF ( ALLOCATED( GOOD_CT_EMIS ) ) DEALLOCATE( GOOD_CT_EMIS )
      IF ( ALLOCATED( GOOD_EMIS    ) ) DEALLOCATE( GOOD_EMIS    )
      IF ( ALLOCATED( Q            ) ) DEALLOCATE( Q            )

      END SUBROUTINE CLEANUP_DIAG51b
!EOC
      END MODULE DIAG51b_MOD
