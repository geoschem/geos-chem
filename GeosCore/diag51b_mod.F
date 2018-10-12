!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: diag51b_mod.F 
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
      USE PhysConstants, ONLY : AIRMW
      USE PRECISION_MOD       ! For GEOS-Chem Precision (fp)

      IMPLICIT NONE
      PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
! 
      PUBLIC  :: CLEANUP_DIAG51b
      PUBLIC  :: DIAG51b
      PUBLIC  :: INIT_DIAG51b
!
! !PRIVATE MEMBER FUNCTIONS:
! 
      PRIVATE :: ACCUMULATE_DIAG51b
      PRIVATE :: GET_LOCAL_TIME
      PRIVATE :: ITS_TIME_FOR_WRITE_DIAG51b
      PRIVATE :: WRITE_DIAG51b
!
! !REMARKS:
!  NOTE by Melissa Sulprizio, 26 May 2015
!  ----------------------------------------------------------------------------
!  The emission options in the timeseries diagnostics were removed 
!  from v10-01 since HEMCO now handles the emission diagnostics. To 
!  utilize the diagnostics capability from HEMCO and output hourly 
!  isoprene emissions, you can follow these steps:
!
!    1. At the top of your HEMCO_Config.rc file, set DiagnFreq to Hourly 
!       and add a line for DiagnFile:
!
!          DiagnPrefix: HEMCO_Diagnostics
!          DiagnFreq: Hourly
!          DiagnFile: DiagnFile.rc
!
!    2.  Create a new text file in your run directory named DiagnFile.rc 
!        and list the emission fields that you would like to be saved out. 
!        For example:
!
!          # Name        Spec ExtNr Cat Hier Dim OutUnit
!          ISOP_BIOG     ISOP 108    1   1   2   kg/m2/s
!
!        NOTE: The ExtNr, Cat, Hier, and Dim values listed above were 
!        obtained from the MEGAN entries in the HEMCO_Config.rc file.
!
!    3. You can then run GEOS-Chem as usual. HEMCO will write out the 
!       specified diagnostics in a netCDF file named 
!       HEMCO_Diagnostics.YYYYMMDDHHmm.nc. I recommend running a short 
!       1-day simulation to make sure the diagnostic output is what 
!       you expect.
!
!  For more details on the HEMCO diagnostics, please see this post 
!  in the HEMCO User’s Guide:
!
!       http://wiki.geos-chem.org/The_HEMCO_User%27s_Guide#Diagnostics
!
!  ND51b tracer numbers:
!  ============================================================================
!  1 - nAdvect   : GEOS-CHEM advected species               [v/v        ]
!  501           : OH concentration                         [molec/cm3  ]
!  502           : NOy concentration                        [v/v        ]
!  503           : Relative Humidity                        [%          ]
!  504           : 3-D Cloud fractions                      [unitless   ]
!  505           : Column optical depths                    [unitless   ]
!  506           : Cloud top heights                        [hPa        ]
!  507           : Air density                              [molec/cm3  ]
!  508           : Total seasalt tracer concentration       [unitless   ]
!  509           : PBL heights                              [m          ]
!  510           : PBL heights                              [levels     ]
!  511           : Grid box heights                         [m          ]
!  512           : PEDGE-$ (Pressure @ level edges          [hPa        ]
!  513           : Sea level pressure                       [hPa        ]
!  514           : Zonal wind (a.k.a. U-wind)               [m/s        ]
!  515           : Meridional wind (a.k.a. V-wind)          [m/s        ]
!  516           : Temperature                              [K          ]
!  517           : Sulfate aerosol optical depth            [unitless   ]
!  518           : Black carbon aerosol optical depth       [unitless   ]
!  519           : Organic carbon aerosol optical depth     [unitless   ]
!  520           : Accumulation mode seasalt optical depth  [unitless   ]
!  521           : Coarse mode seasalt optical depth        [unitless   ]
!  522           : Total dust optical depth                 [unitless   ]
!  523-529       : Size resolved dust optical depth         [unitless   ]
!  530           : PAR direct                               [hPa        ]
!  531           : PAR diffuse                              [hPa        ]
!  532           : Daily LAI                                [hPa        ]
!  533           : Temperature at 2m                        [K          ]
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
!  19 Mar 2014 - M. Sulprizio- Updated to allow for more than 75 tracers
!  10 Nov 2014 - M. Yannetti - Added PRECISION_MOD
!  17 Dec 2014 - R. Yantosca - Leave time/date variables as 8-byte
!  15 Apr 2015 - M. Sulprizio- Comment out emission diagnostics for now. These
!                              are now handled by the HEMCO.
!  22 Jun 2016 - M. Yannetti - Replace TCVV with spec db and physical const
!  29 Nov 2016 - R. Yantosca - grid_mod.F90 is now gc_grid_mod.F90
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !PRIVATE TYPES:
!
#if defined( BPCH_DIAG ) || defined( BPCH_TIMESER )
      !=================================================================
      ! MODULE VARIABLES
      ! GOOD             : Array denoting grid boxes w/in LT limits
      ! GOOD_CT          : # of "good" times per grid box
      ! GOOD_CT_CHEM     : # of "good" chemistry timesteps
      ! COUNT_CHEM3D     : Counter for 3D chemistry boxes
      ! I0               : Offset between global & nested grid
      ! J0               : Offset between global & nested grid
      ! IOFF             : Longitude offset
      ! JOFF             : Latitude offset
      ! LOFF             : Altitude offset
      ! ND51b_NI         : Number of longitudes in DIAG51b region 
      ! ND51b_NJ         : Number of latitudes  in DIAG51b region
      ! ND51b_NL         : Number of levels     in DIAG51b region
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
      INTEGER              :: IOFF,           JOFF,     LOFF
      INTEGER              :: I0,             J0
      INTEGER              :: ND51b_NI,       ND51b_NJ, ND51b_NL
      INTEGER              :: HALFPOLAR
      INTEGER, PARAMETER   :: CENTER180=1
      REAL*4               :: LONRES,         LATRES
      REAL(f8)             :: TAU0,           TAU1
      CHARACTER(LEN=20)    :: MODELNAME
      CHARACTER(LEN=40)    :: RESERVED = ''
      CHARACTER(LEN=80)    :: TITLE

      ! LUN for ND51b output file
      INTEGER              :: IU_ND51b

      ! Arrays
      INTEGER, ALLOCATABLE :: GOOD(:)
      INTEGER, ALLOCATABLE :: GOOD_CT(:)
      INTEGER, ALLOCATABLE :: COUNT_CHEM3D(:,:,:)
      REAL(fp),ALLOCATABLE :: Q(:,:,:,:)

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
      !REAL(fp),  PARAMETER   :: HR1        = 13d0
      !REAL(fp),  PARAMETER   :: HR2        = 17d0
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
      !REAL(fp),  PARAMETER   :: HR1        = 10d0
      !REAL(fp),  PARAMETER   :: HR2        = 12d0
      !CHARACTER(LEN=255)   :: FILENAME   ='ts10_12pm.bpch'
      !=================================================================
#endif

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
! !DESCRIPTION:  Subroutine DIAG51b generates time series (averages from !
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
      USE ErrCode_Mod
      USE Input_Opt_Mod,      ONLY : OptInput
      USE State_Chm_Mod,      ONLY : ChmState
      USE State_Met_Mod,      ONLY : MetState
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
#if defined( BPCH_DIAG ) || defined( BPCH_TIMESER )
      REAL(fp) :: TAU_W
#endif

      !=================================================================
      ! DIAG51b begins here!
      !=================================================================

      ! Assume success
      RC = GC_SUCCESS

#if defined( BPCH_DIAG ) || defined( BPCH_TIMESER )

      ! Construct array of where local times are between HR1, HR2
      CALL GET_LOCAL_TIME( Input_Opt )

      ! Accumulate data in the Q array
      CALL ACCUMULATE_DIAG51b
     &   ( am_I_Root, Input_Opt, State_Met, State_Chm, RC )

      ! Write data to disk at the proper time
      IF ( ITS_TIME_FOR_WRITE_DIAG51b( Input_Opt, TAU_W ) .and. 
     &     Input_Opt%DO_DIAG_WRITE ) THEN
         CALL WRITE_DIAG51b( am_I_Root, Input_Opt, State_Chm, TAU_W, RC)
      ENDIF
#endif

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
      SUBROUTINE GET_LOCAL_TIME( Input_Opt )
!
! !USES:
!
      USE Input_Opt_Mod, ONLY : OptInput
#if defined( BPCH_DIAG ) || defined( BPCH_TIMESER )
      USE CMN_SIZE_MOD        ! Size parameters
      USE TIME_MOD,      ONLY : GET_LOCALTIME
      USE TIME_MOD,      ONLY : GET_TS_DYN
#endif
!
! !INPUT ARGUMENTS
!
      TYPE(OptInput), INTENT(IN) :: Input_Opt  ! Input opts 
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
!  (6 ) Bug fix: LT should be REAL(fp) and not INTEGER (ccarouge, 12/10/08)
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
#if defined( BPCH_DIAG ) || defined( BPCH_TIMESER )
      INTEGER  :: I
      REAL(fp) :: LT, TS_DYN

      !=================================================================
      ! GET_LOCAL_TIME begins here!
      !=================================================================
      TS_DYN = GET_TS_DYN() / 3600e+0_fp

      DO I = 1, IIPAR

         ! Get local time
         LT = GET_LOCALTIME( I, 1, 1 ) - TS_DYN
         IF ( LT < 0  ) LT = LT + 24e+0_fp

         ! GOOD indicates which boxes have local times between HR1 and HR2
         IF ( LT >= Input_Opt%ND51b_HR1 .and.
     &        LT <= Input_Opt%ND51b_HR2 ) THEN
            GOOD(I) = 1
         ENDIF
      ENDDO
#endif

      END SUBROUTINE GET_LOCAL_TIME
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: accumulate_diag51b
!
! !DESCRIPTION: Subroutine ACCUMULATE\_DIAG51b accumulates tracers into the 
!  Q array. 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE ACCUMULATE_DIAG51b( am_I_Root, Input_Opt,
     &                               State_Met, State_Chm, RC )
!
! !USES:
!
      USE ErrCode_Mod
      USE Input_Opt_Mod,      ONLY : OptInput
      USE State_Chm_Mod,      ONLY : ChmState
      USE State_Chm_Mod,      ONLY : Ind_
      USE State_Met_Mod,      ONLY : MetState
#if defined( BPCH_DIAG ) || defined( BPCH_TIMESER )
      USE CMN_FJX_MOD,        ONLY : ODAER, ODMDUST, NRH, NDUST
      USE CMN_FJX_MOD,        ONLY : IWVSELECT, ACOEF_WV, BCOEF_WV
      USE CMN_O3_MOD               ! SAVEOH
      USE PhysConstants            ! SCALE_HEIGHT, XNUMOLAIR
      USE PBL_MIX_MOD,        ONLY : GET_PBL_TOP_L,    GET_PBL_TOP_m
      USE TIME_MOD,           ONLY : GET_ELAPSED_SEC,  GET_TS_CHEM 
      USE TIME_MOD,           ONLY : TIMESTAMP_STRING, GET_TS_DYN
      USE TIME_MOD,           ONLY : GET_TS_DIAG,      GET_TS_EMIS
#endif
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
!  06 Nov 2014 - R. Yantosca - Now use State_Met%AIRDEN(I,J,L)
!  06 Nov 2014 - R. Yantosca - Now use State_Met%CLDF(I,J,L)
!  06 Nov 2014 - R. Yantosca - Now use State_Met%OPTD(I,J,L)
!  26 Feb 2015 - E. Lundgren - Replace GET_PEDGE with State_Met%PEDGE.
!  24 Mar 2015 - E. Lundgren - Remove dependency on tracer_mod
!  25 Mar 2015 - E. Lundgren - Change tracer units from kg to kg/kg and
!                              remove AD pointer
!  16 Jun 2016 - K. Yu       - Now define species ID's with the Ind_ function
!  17 Jun 2016 - R. Yantosca - Only define species ID's on the first call
!  11 Aug 2016 - R. Yantosca - Remove temporary tracer-removal code
!  06 Feb 2018 - E. Lundgren - Change GET_ELAPSED_MIN to GET_ELAPSED_SEC to
!                              match new timestep unit of seconds
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
#if defined( BPCH_DIAG ) || defined( BPCH_TIMESER )
      ! SAVEd scalars
      LOGICAL, SAVE     :: FIRST = .TRUE.
      LOGICAL, SAVE     :: IS_FULLCHEM, IS_SEASALT
      LOGICAL, SAVE     :: IS_CLDTOPS,  IS_NOy,    IS_OPTD, IS_SLP
      INTEGER, SAVE     :: id_HNO3,     id_HNO4,   id_N2O5, id_NO
      INTEGER, SAVE     :: id_PAN,      id_NPMN,   id_PPN,  id_O3
      INTEGER, SAVE     :: id_R4N2,     id_SALA,   id_SALC, id_NO2
      INTEGER, SAVE     :: id_IPMN,     id_OH

      ! Scalars
      LOGICAL           :: IS_CHEM,     IS_DIAG,   IS_EMIS
      LOGICAL           :: LINTERP
      INTEGER           :: H, I, J, K, L, M, N, ISPC, NA, nAdvect
      INTEGER           :: PBLINT,  R, X, Y, W, XMIN
      INTEGER           :: IWV
      REAL(fp)          :: C1, C2, PBLDEC, TEMPBL, TMP, SCALEAODnm
      CHARACTER(LEN=16) :: STAMP

      ! Aerosol types (rvm, aad, bmy, 7/20/04)
      INTEGER           :: IND(6) = (/ 22, 29, 36, 43, 50, 15 /)

      ! Pointers
      REAL(fp), POINTER :: Spc(:,:,:,:)
#endif

      !=================================================================
      ! ACCUMULATE_DIAG51b begins here!
      !=================================================================

      ! Assume success
      RC        =  GC_SUCCESS

#if defined( BPCH_DIAG ) || defined( BPCH_TIMESER )

      ! Number of advected species
      nAdvect   =  State_Chm%nAdvect

      ! First-time setup
      IF ( FIRST ) THEN

         ! Define species ID flags
         id_HNO3 = Ind_('HNO3')
         id_HNO4 = Ind_('HNO4')
         id_N2O5 = Ind_('N2O5')
         id_NO   = Ind_('NO')
         id_PAN  = Ind_('PAN')
         id_NPMN = Ind_('NPMN')
         id_IPMN = Ind_('IPMN')
         id_PPN  = Ind_('PPN')
         id_O3   = Ind_('O3')
         id_R4N2 = Ind_('R4N2')
         id_SALA = Ind_('SALA')
         id_SALC = Ind_('SALC')
         id_NO2  = Ind_('NO2')
         id_OH   = Ind_('OH')

         ! Set logical flags on first call
         IS_OPTD     = ASSOCIATED( State_Met%OPTD    )
         IS_CLDTOPS  = ASSOCIATED( State_Met%CLDTOPS )
         IS_SLP      = ASSOCIATED( State_Met%SLP     )
         IS_FULLCHEM = Input_Opt%ITS_A_FULLCHEM_SIM
         IS_SEASALT  = ( id_SALA > 0 .and. id_SALC > 0 )
         IS_NOy      = ( IS_FULLCHEM .and. id_NO   > 0 .and.
     &                   id_NO2  > 0 .and. id_PAN  > 0 .and.
     &                   id_HNO3 > 0 .and. id_NPMN > 0 .and.
     &                   id_PPN  > 0 .and. id_R4N2 > 0 .and.
     &                   id_N2O5 > 0 .and. id_HNO4 > 0 .and.
     &                   id_IPMN > 0 ) 

         ! Reset first-time flag
         FIRST       = .FALSE.
      ENDIF

      ! Force accumulation on every "heartbeat" timestep (bmy, 10/4/18)
      IS_DIAG = .TRUE.

      ! Echo info
      STAMP = TIMESTAMP_STRING()
      WRITE( 6, 100 ) STAMP
 100  FORMAT( '     - DIAG51b: Accumulation at ', a )
      
      !=================================================================
      ! Archive tracers into accumulating array Q 
      !=================================================================

      ! Archive counter array of good points 
      IF ( IS_DIAG ) THEN
         DO X = 1, ND51b_NI
            I          = GET_I( X )
            GOOD_CT(X) = GOOD_CT(X) + GOOD(I)
         ENDDO
      ENDIF

      ! Also increment 3-D counter for boxes in the chemistry grid
      IF ( IS_FULLCHEM ) THEN
         
         ! Loop over levels
!$OMP PARALLEL DO 
!$OMP+DEFAULT( SHARED ) 
!$OMP+PRIVATE( X, Y, K, I, J, L )
!$OMP+SCHEDULE( DYNAMIC )
         DO K = 1, ND51b_NL
            L = LOFF + K

         ! Loop over latitudes 
         DO Y = 1, ND51b_NJ
            J = JOFF + Y

         ! Loop over longitudes
         DO X = 1, ND51b_NI
            I = GET_I( X )

            ! Only increment if we are in the 
            IF (  State_Met%InChemGrid(I,J,L) ) THEN
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
         !Determine if optical properties need interpolating
         !the LUT wavelengths in IWVSELECT will match if no interpolation
         !is needed. (output is only for the first requested wavelength)
         IF(IWVSELECT(1,1).EQ.IWVSELECT(2,1)) THEN
            LINTERP=.FALSE.
         ELSE
            LINTERP=.TRUE.
         ENDIF

         ! Initialize pointers
         Spc => State_Chm%Species

!$OMP PARALLEL DO 
!$OMP+DEFAULT( SHARED ) 
!$OMP+PRIVATE( W, N, X, Y, K, I, J, L, TMP, H, R, SCALEAODnm )
!$OMP+SCHEDULE( DYNAMIC )
         DO W = 1, Input_Opt%N_ND51b

         ! ND51b Tracer number
         N = Input_Opt%ND51b_TRACERS(W)

         ! Loop over levels
         DO K = 1, ND51b_NL
            L = LOFF + K

         ! Loop over latitudes 
         DO Y = 1, ND51b_NJ
            J = JOFF + Y

         ! Loop over longitudes
         DO X = 1, ND51b_NI
            I = GET_I( X )

            ! Archive by simulation 
            IF ( N <= nAdvect ) THEN

               !--------------------------------------
               ! GEOS-CHEM advected species [v/v]
               !--------------------------------------

               ! Archive afternoon points
               Q(X,Y,K,W) = Q(X,Y,K,W) + 
     &                      ( Spc(I,J,L,N) * ( AIRMW 
     &                      / State_Chm%SpcData(N)%Info%emMW_g )
     &                        * GOOD(I) )

! NOTE: We can reactivate once we know what the proper unit conversion is
! (bmy, 11/20/17)
!            ELSE IF ( N == 501 .and. IS_FULLCHEM ) THEN
!
!               !--------------------------------------
!               ! OH [molec/cm3]
!               ! NOTE: Only archive at chem timestep
!               !--------------------------------------
!
!               ! Accumulate data
!               Q(X,Y,K,W) = Q(X,Y,K,W) + 
!     &              ( SAVEOH(I,J,L,id_OH) * GOOD(X) )

            ELSE IF ( N == 502 .and. IS_NOy ) THEN

               !--------------------------------------
               ! NOy [v/v]
               !--------------------------------------
  
               ! Temp variable for accumulation
               TMP = 0e+0_fp
            
               ! NO
               TMP = TMP + ( ( AIRMW 
     &                       / State_Chm%SpcData(id_NO)%Info%emMW_g )         
     &                     * GOOD(I) * Spc(I,J,L,id_NO)    )

               ! NO2
               TMP = TMP + ( ( AIRMW 
     &                        / State_Chm%SpcData(id_NO2)%Info%emMW_g )        
     &                       * GOOD(I) * Spc(I,J,L,id_NO2)   )
               ! PAN
               TMP = TMP + ( ( AIRMW 
     &                        / State_Chm%SpcData(id_PAN)%Info%emMW_g )        
     &                       * GOOD(I) * Spc(I,J,L,id_PAN)   )

               ! HNO3
               TMP = TMP + ( ( AIRMW 
     &                        / State_Chm%SpcData(id_HNO3)%Info%emMW_g )       
     &                       * GOOD(I) * Spc(I,J,L,id_HNO3)  )
            
               ! NPMN
               TMP = TMP + ( ( AIRMW 
     &                        / State_Chm%SpcData(id_NPMN)%Info%emMW_g )        
     &                       * GOOD(I) * Spc(I,J,L,id_NPMN)   )

               ! IPMN
               TMP = TMP + ( ( AIRMW 
     &                        / State_Chm%SpcData(id_IPMN)%Info%emMW_g )        
     &                       * GOOD(I) * Spc(I,J,L,id_IPMN)   )

               ! PPN
               TMP = TMP + ( ( AIRMW 
     &                        / State_Chm%SpcData(id_PPN)%Info%emMW_g )       
     &                       * GOOD(I) * Spc(I,J,L,id_PPN)   )
 
               ! R4N2
               TMP = TMP + ( ( AIRMW 
     &                        / State_Chm%SpcData(id_R4N2)%Info%emMW_g )       
     &                       * GOOD(I) * Spc(I,J,L,id_R4N2)  )
            
               ! N2O5
               TMP = TMP + ( 2e+0_fp * ( AIRMW 
     &                        / State_Chm%SpcData(id_N2O5)%Info%emMW_g ) 
     &                       * GOOD(I) * Spc(I,J,L,id_N2O5)  )
                        
               ! HNO4
               TMP = TMP + ( ( AIRMW 
     &                        / State_Chm%SpcData(id_HNO4)%Info%emMW_g )       
     &                       * GOOD(I) * Spc(I,J,L,id_HNO4)  )

               ! Save afternoon points
               Q(X,Y,K,W) = Q(X,Y,K,W) + TMP
    
            ELSE IF ( N == 503 ) THEN

               !---------------------------------------
               ! RELATIVE HUMIDITY [%]
               !---------------------------------------
               Q(X,Y,K,W) = Q(X,Y,K,W) + ( State_Met%RH(I,J,L) *
     &                      GOOD(I) )

            ELSE IF ( N == 504 ) THEN

               !--------------------------------------
               ! 3-D CLOUD FRACTIONS [unitless]
               !--------------------------------------
               Q(X,Y,K,W) = Q(X,Y,K,W) + ( State_Met%CLDF(I,J,L) *
     &                      GOOD(I) )

            ELSE IF ( N == 505 .and. IS_OPTD ) THEN

               !--------------------------------------
               ! COLUMN OPTICAL DEPTHS [unitless]
               !--------------------------------------
               Q(X,Y,1,W) = Q(X,Y,1,W) + ( State_Met%OPTD(I,J,L) *
     &                      GOOD(I) )

            ELSE IF ( N == 506 .and. IS_CLDTOPS ) THEN

               !--------------------------------------
               ! CLOUD TOP HEIGHTS [mb]
               !--------------------------------------
               IF ( K == 1 ) THEN
                  TMP      = State_Met%PEDGE(I,J,State_Met%CLDTOPS(I,J))
                  Q(X,Y,K,W) = Q(X,Y,K,W) + ( TMP * GOOD(I) )
               ENDIF
              
            ELSE IF ( N == 507 ) THEN

               !--------------------------------------
               ! AIR DENSITY [molec/cm3] 
               !--------------------------------------
               Q(X,Y,K,W) = Q(X,Y,K,W) + ( State_Met%AIRDEN(I,J,L) *
     &                      XNUMOLAIR  * 1d-6 * GOOD(I) )
 
            ELSE IF ( N == 508 .and. IS_SEASALT ) THEN

               !--------------------------------------
               ! TOTAL SEASALT TRACER [v/v]
               !--------------------------------------
               Q(X,Y,K,W) = Q(X,Y,K,W) +
     &                      ( Spc(I,J,L,id_SALA)   + 
     &                        Spc(I,J,L,id_SALC) ) *
     &                        ( AIRMW / 
     &                          State_Chm%SpcData(id_SALA)%Info%emMW_g ) 
     &                          * GOOD(I)

            ELSE IF ( N == 509 ) THEN

               !--------------------------------------
               ! PBL HEIGHTS [m] 
               !--------------------------------------
               IF ( K == 1 ) THEN
                  Q(X,Y,K,W) = Q(X,Y,K,W) + 
     &                         ( GET_PBL_TOP_m( I, J ) * GOOD(I) )  
               ENDIF

            ELSE IF ( N == 510 ) THEN

               !--------------------------------------
               ! PBL HEIGHTS [layers] 
               !--------------------------------------
               IF ( K == 1 ) THEN
                  Q(X,Y,K,W) = Q(X,Y,K,W) +
     &                         ( GET_PBL_TOP_L( I, J ) * GOOD(I) )
               ENDIF

            ELSE IF ( N == 511 ) THEN

               !--------------------------------------
               ! GRID BOX HEIGHTS [m]
               !--------------------------------------
               Q(X,Y,K,W) = Q(X,Y,K,W) + ( State_Met%BXHEIGHT(I,J,L) *
     &                      GOOD(I) )

            ELSE IF ( N == 512 ) THEN

               !--------------------------------------
               ! PEDGE-$ (prs @ level edges) [hPa]
               !--------------------------------------
               Q(X,Y,K,W) = Q(X,Y,K,W) + 
     &                      ( State_Met%PEDGE(I,J,K) * GOOD(I) )

            ELSE IF ( N == 513 .and. IS_SLP ) THEN

               !--------------------------------------
               ! SEA LEVEL PRESSURE [hPa]
               !--------------------------------------
               IF ( K == 1 ) THEN
                  Q(X,Y,K,W) = Q(X,Y,K,W) + ( State_Met%SLP(I,J) *
     &                         GOOD(I) )
               ENDIF

            ELSE IF ( N == 514 ) THEN

               !--------------------------------------
               ! ZONAL (U) WIND [M/S]
               !--------------------------------------
               Q(X,Y,K,W) = Q(X,Y,K,W) + ( State_Met%U(I,J,L) *
     &                      GOOD(I) )

            ELSE IF ( N == 515 ) THEN

               !--------------------------------------
               ! MERIDIONAL (V) WIND [M/S]
               !--------------------------------------
               Q(X,Y,K,W) = Q(X,Y,K,W) + ( State_Met%V(I,J,L) *
     &                      GOOD(I) )

            ELSE IF ( N == 516 ) THEN 

               !--------------------------------------
               ! TEMPERATURE [K]
               !--------------------------------------
               Q(X,Y,K,W) = Q(X,Y,K,W) + ( State_Met%T(I,J,L) *
     &                      GOOD(I) )

            ELSE IF ( N == 517 ) THEN

               !--------------------------------------
               ! SULFATE AOD [unitless]
               ! for wavelengths set in Radiation Menu
               !
               ! NOTE: Only archive at chem timestep
               !--------------------------------------
               ISPC = 1 !sulfate
               IF ( .not. LINTERP ) THEN
                  ! Accumulate
                  Q(X,Y,K,W) = Q(X,Y,K,W) +
     &               (ODAER(I,J,L,IWVSELECT(1,1),ISPC) * GOOD(X))
               ELSE
                  ! Interpolated using angstrom exponent between
                  ! Closest available wavelengths
                  ! (coefs pre-calculated in CALC_AOD (RD_AOD.F)
                  !catch any zero values before interpolation
                  IF ((ODAER(I,J,L,IWVSELECT(2,1),ISPC).GT.0).AND.
     &                (ODAER(I,J,L,IWVSELECT(1,1),ISPC).GT.0)) THEN
                    Q(X,Y,K,W) = Q(X,Y,K,W) + GOOD(X)*
     &              (ODAER(I,J,L,IWVSELECT(2,1),ISPC)*ACOEF_WV(1)**
     &              (BCOEF_WV(1)*LOG(ODAER(I,J,L,IWVSELECT(1,1),ISPC)/
     &                               ODAER(I,J,L,IWVSELECT(2,1),ISPC))))
                  ENDIF
               ENDIF

            ELSE IF ( N == 518 ) THEN

               !--------------------------------------
               ! BLACK CARBON AOD [unitless]
               ! for wavelengths set in Radiation Menu
               !
               ! NOTE: Only archive at chem timestep
               !--------------------------------------
               ISPC = 2 !BC
               IF ( .not. LINTERP ) THEN
                  ! Accumulate
                  Q(X,Y,K,W) = Q(X,Y,K,W) +
     &               (ODAER(I,J,L,IWVSELECT(1,1),ISPC) * GOOD(X))
               ELSE
                  ! Interpolated using angstrom exponent between
                  ! Closest available wavelengths
                  ! (coefs pre-calculated in CALC_AOD (RD_AOD.F)
                  !catch any zero values before interpolation
                  IF ((ODAER(I,J,L,IWVSELECT(2,1),ISPC).GT.0).AND.
     &                (ODAER(I,J,L,IWVSELECT(1,1),ISPC).GT.0)) THEN
                    Q(X,Y,K,W) = Q(X,Y,K,W) + GOOD(X)*
     &              (ODAER(I,J,L,IWVSELECT(2,1),ISPC)*ACOEF_WV(1)**
     &              (BCOEF_WV(1)*LOG(ODAER(I,J,L,IWVSELECT(1,1),ISPC)/
     &                               ODAER(I,J,L,IWVSELECT(2,1),ISPC))))
                  ENDIF
               ENDIF

            ELSE IF ( N == 519 ) THEN

               !--------------------------------------
               ! ORG CARBON AOD [unitless]
               ! for wavelengths set in Radiation Menu
               !
               ! NOTE: Only archive at chem timestep
               !--------------------------------------
               ISPC = 3 !OC
               IF ( .not. LINTERP ) THEN
                  ! Accumulate
                  Q(X,Y,K,W) = Q(X,Y,K,W) +
     &               (ODAER(I,J,L,IWVSELECT(1,1),ISPC) * GOOD(X))
               ELSE
                  ! Interpolated using angstrom exponent between
                  ! Closest available wavelengths
                  ! (coefs pre-calculated in CALC_AOD (RD_AOD.F)
                  !catch any zero values before interpolation
                  IF ((ODAER(I,J,L,IWVSELECT(2,1),ISPC).GT.0).AND.
     &                (ODAER(I,J,L,IWVSELECT(1,1),ISPC).GT.0)) THEN
                    Q(X,Y,K,W) = Q(X,Y,K,W) + GOOD(X)*
     &              (ODAER(I,J,L,IWVSELECT(2,1),ISPC)*ACOEF_WV(1)**
     &              (BCOEF_WV(1)*LOG(ODAER(I,J,L,IWVSELECT(1,1),ISPC)/
     &                               ODAER(I,J,L,IWVSELECT(2,1),ISPC))))
                  ENDIF
               ENDIF

            ELSE IF ( N == 520 ) THEN

               !--------------------------------------
               ! ACCUM SEASALT AOD [unitless]
               ! for wavelengths set in Radiation Menu
               !
               ! NOTE: Only archive at chem timestep
               !--------------------------------------
               ISPC = 4 !SSa
               IF ( .not. LINTERP ) THEN
                  ! Accumulate
                  Q(X,Y,K,W) = Q(X,Y,K,W) +
     &               (ODAER(I,J,L,IWVSELECT(1,1),ISPC) * GOOD(X))
               ELSE
                  ! Interpolated using angstrom exponent between
                  ! Closest available wavelengths
                  ! (coefs pre-calculated in CALC_AOD (RD_AOD.F)
                  !catch any zero values before interpolation
                  IF ((ODAER(I,J,L,IWVSELECT(2,1),ISPC).GT.0).AND.
     &                (ODAER(I,J,L,IWVSELECT(1,1),ISPC).GT.0)) THEN
                    Q(X,Y,K,W) = Q(X,Y,K,W) + GOOD(X)*
     &              (ODAER(I,J,L,IWVSELECT(2,1),ISPC)*ACOEF_WV(1)**
     &              (BCOEF_WV(1)*LOG(ODAER(I,J,L,IWVSELECT(1,1),ISPC)/
     &                               ODAER(I,J,L,IWVSELECT(2,1),ISPC))))
                  ENDIF
               ENDIF

            ELSE IF ( N == 521 ) THEN

               !--------------------------------------
               ! COARSE SEASALT AOD [unitless]
               ! for wavelengths set in Radiation Menu
               !
               ! NOTE: Only archive at chem timestep
               !--------------------------------------
               ISPC = 5 !SSc
               IF ( .not. LINTERP ) THEN
                  ! Accumulate
                  Q(X,Y,K,W) = Q(X,Y,K,W) +
     &               (ODAER(I,J,L,IWVSELECT(1,1),ISPC) * GOOD(X))
               ELSE
                  ! Interpolated using angstrom exponent between
                  ! Closest available wavelengths
                  ! (coefs pre-calculated in CALC_AOD (RD_AOD.F)
                  !catch any zero values before interpolation
                  IF ((ODAER(I,J,L,IWVSELECT(2,1),ISPC).GT.0).AND.
     &                (ODAER(I,J,L,IWVSELECT(1,1),ISPC).GT.0)) THEN
                    Q(X,Y,K,W) = Q(X,Y,K,W) + GOOD(X)*
     &              (ODAER(I,J,L,IWVSELECT(2,1),ISPC)*ACOEF_WV(1)**
     &              (BCOEF_WV(1)*LOG(ODAER(I,J,L,IWVSELECT(1,1),ISPC)/
     &                               ODAER(I,J,L,IWVSELECT(2,1),ISPC))))
                  ENDIF
               ENDIF

            ELSE IF ( N == 522 ) THEN               

               !--------------------------------------
               ! TOTAL DUST OPTD [unitless]
               ! for wavelengths set in Radiation Menu
               !
               ! NOTE: Only archive at chem timestep
               !--------------------------------------
               DO R = 1, NDUST

                  IF ( .not. LINTERP ) THEN
                     Q(X,Y,K,W) = Q(X,Y,K,W) + GOOD(X) *
     &                            ODMDUST(I,J,L,IWVSELECT(1,1),R)
                  ELSE
                     ! Interpolated using angstrom exponent between
                     ! Closest available wavelengths
                     ! (coefs pre-calculated in CALC_AOD (RD_AOD.F)
                     ! AOD sometimes zero (if Q zero), must catch this
                     IF ((ODMDUST(I,J,L,IWVSELECT(1,1),R).GT.0).AND.
     &                   (ODMDUST(I,J,L,IWVSELECT(2,1),R).GT.0)) THEN

                        Q(X,Y,K,W) = Q(X,Y,K,W) + GOOD(X) *
     &                               (ODMDUST(I,J,L,IWVSELECT(2,1),R)*
     &                               ACOEF_WV(1)**(BCOEF_WV(1)*LOG(
     &                               ODMDUST(I,J,L,IWVSELECT(1,1),R)/
     &                               ODMDUST(I,J,L,IWVSELECT(2,1),R))))
                     ENDIF
                  ENDIF

               ENDDO

            ELSE IF ( ( N >= 523 ) .and. ( N <= 529) ) THEN

               !--------------------------------------
               ! Dust BINS 1-7 optical depth [unitless]
               ! for wavelengths set in Radiation Menu
               !
               ! NOTE: Only archive at chem timestep
               !--------------------------------------
               R = N - 522

               ! Accumulate
               IF ( .not. LINTERP ) THEN
                  Q(X,Y,K,W) = Q(X,Y,K,W) + GOOD(X) *
     &                         ODMDUST(I,J,L,IWVSELECT(1,1),R)
               ELSE
                  ! Interpolated using angstrom exponent between
                  ! Closest available wavelengths
                  ! (coefs pre-calculated in CALC_AOD (RD_AOD.F)
                  ! AOD sometimes zero (if Q zero), must catch this
                  IF ((ODMDUST(I,J,L,IWVSELECT(1,1),R).GT.0).AND.
     &                (ODMDUST(I,J,L,IWVSELECT(2,1),R).GT.0)) THEN

                     Q(X,Y,K,W) = Q(X,Y,K,W) + GOOD(X) *
     &                            (ODMDUST(I,J,L,IWVSELECT(2,1),R)*
     &                            ACOEF_WV(1)**(BCOEF_WV(1)*LOG(
     &                            ODMDUST(I,J,L,IWVSELECT(1,1),R)/
     &                            ODMDUST(I,J,L,IWVSELECT(2,1),R))))
                  ENDIF
               ENDIF

! =====================================================================
! Added with MEGAN v2.1. (ccc, 11/20/09)

            ELSE IF ( N == 530 ) THEN 

               !--------------------------------------
               ! PAR DR [W/m2] (mpb,2009)
               !--------------------------------------
              
               IF ( K == 1 ) THEN
              
                  Q(X,Y,K,W) = Q(X,Y,K,W) +
     &                       ( State_Met%PARDR(I,J) * GOOD(I) )
               ENDIF

            ELSE IF ( N == 531 ) THEN 

               !--------------------------------------
               ! PAR DF [W/m2] (mpb,2009)
               !--------------------------------------
              
               IF ( K == 1 ) THEN
              
                  Q(X,Y,K,W) =  Q(X,Y,K,W) +
     &                       (  State_Met%PARDF(I,J) * GOOD(I) )
               ENDIF

            ELSE IF ( N == 532 ) THEN 

               !--------------------------------------
               ! DAILY LAI  [cm2/cm2] (mpb,2009)
               !--------------------------------------
              
               IF ( K == 1 ) THEN
              
                  Q(X,Y,K,W) =  Q(X,Y,K,W) +
     &                     ( State_Met%MODISLAI(I,J) * GOOD(I) ) 
               ENDIF


            ELSE IF ( N == 533 ) THEN

               !--------------------------------------
               ! T at 2m [K] (mpb,2009)
               !--------------------------------------
              
               IF ( K == 1 ) THEN

                  Q(X,Y,K,W) = Q(X,Y,K,W) +
     &                       ( State_Met%TS(I,J) * GOOD(I) )

               ENDIF

!##############################################################################
!###  HEMCO NOW ARCHIVES THE BIOGENIC EMISSION TIMESERIES DIAGNOSTICS,      ###
!###  AS WELL AS THE SOIL NOx TIMESERIES DIAGNOSTICS.                       ###
!###  SEE THE NOTE ABOVE FOR INSTRUCTIONS ON HOW TO SAVE THE EQUIVALENT     ###
!###  OUTPUTS TO A netCDF FILE. (bmy, mps, 5/26/15)                         ###
!##############################################################################

            ELSE

               ! Skip other tracers
               CYCLE

            ENDIF
         ENDDO
         ENDDO
         ENDDO
         ENDDO 
!$OMP END PARALLEL DO

         GOOD(:) = 0

         ! Free pointers
         Spc => NULL()

      ENDIF
#endif

      END SUBROUTINE ACCUMULATE_DIAG51b
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: its_time_for_write_diag51b
!
! !DESCRIPTION: Function ITS\_TIME\_FOR\_WRITE\_DIAG51b returns TRUE if it's 
!  time to write the ND51b bpch file to disk.  We test the time at the next 
!  dynamic timestep so that we can write to disk properly.
!\\
!\\
! !INTERFACE:
!
      FUNCTION ITS_TIME_FOR_WRITE_DIAG51b( Input_Opt, TAU_W )
     &   RESULT( ITS_TIME )
!
! !USES:
!
      USE Input_Opt_Mod, ONLY : OptInput
#if defined( BPCH_DIAG ) || defined( BPCH_TIMESER )
      USE ERROR_MOD,     ONLY : GEOS_CHEM_STOP
      USE TIME_MOD,      ONLY : GET_HOUR
      USE TIME_MOD,      ONLY : GET_MINUTE
      USE TIME_MOD,      ONLY : GET_SECOND
      USE TIME_MOD,      ONLY : GET_TAU
      USE TIME_MOD,      ONLY : GET_TAUb
      USE TIME_MOD,      ONLY : GET_TAUe
      USE TIME_MOD,      ONLY : GET_TS_DYN
      USE TIME_MOD,      ONLY : GET_TS_DIAG
#endif
!
! !INPUT PARAMETERS:
!
      TYPE(OptInput), INTENT(IN) :: Input_Opt   ! Input options
!
! !OUTPUT PARAMETERS:
!
      REAL(fp), INTENT(OUT) :: TAU_W   ! TAU at time of disk write
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
      LOGICAL  :: ITS_TIME

#if defined( BPCH_DIAG ) || defined( BPCH_TIMESER )
      REAL(fp) :: TAU, HOUR, DYN, TS_DIAG

      !=================================================================
      ! ITS_TIME_FOR_WRITE_DIAG51b begins here!
      !=================================================================

      ! Initialize
      ITS_TIME = .FALSE.

      ! Add a check for the time to save. Must be a multiple of TS_DIAG
      ! (ccc, 7/21/09)
      TS_DIAG = ( GET_TS_DIAG() / 3600e+0_fp )

      ! NOTE: Change from equality to greater than an a small number, 
      ! because 20min or 40min timesteps are irrational numbers in hours
      IF ( MOD(Input_Opt%ND51b_HR_WRITE, TS_DIAG) > 1e-5_fp ) THEN
         WRITE( 6, 100 ) Input_Opt%ND51b_HR_WRITE, TS_DIAG
 100     FORMAT( 'The ND51b output frequency must be a multiple '
     &        'of the largest time step:', f9.2, f9.2 )
         CALL GEOS_CHEM_STOP
      ENDIF

      ! Current TAU, Hour, and Dynamic Timestep [hrs]
      TAU      = GET_TAU()
      HOUR     = ( GET_SECOND() / 3600e+0_fp ) + 
     &           ( GET_MINUTE() / 60e+0_fp ) + GET_HOUR()
      DYN      = ( GET_TS_DYN() / 3600e+0_fp )

      ! If first timestep, return FALSE
      IF ( TAU == GET_TAUb() ) RETURN

      ! If the next dyn timestep is the hour of day
      ! when we have to save to disk, return TRUE
      IF ( MOD( HOUR, 24e+0_fp ) == Input_Opt%ND51b_HR_WRITE ) THEN
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
#endif

      END FUNCTION ITS_TIME_FOR_WRITE_DIAG51b
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: write_diag51b
!
! !DESCRIPTION: Subroutine WRITE\_DIAG51b computes the time-average of 
!  quantities between local time limits ND51b\_HR1 and ND51b\_HR2 and writes 
!  them to a bpch file or HDF5 file.  Arrays and counters are also zeroed 
!  for the next diagnostic interval.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE WRITE_DIAG51b( am_I_Root, Input_Opt, 
     &                          State_Chm, TAU_W, RC )
!
! !USES:
!

      USE ErrCode_Mod
      USE Input_Opt_Mod,      ONLY : OptInput
      USE State_Chm_Mod,      ONLY : ChmState    
#if defined( BPCH_DIAG ) || defined( BPCH_TIMESER )
      USE BPCH2_MOD,          ONLY : BPCH2
      USE BPCH2_MOD,          ONLY : OPEN_BPCH2_FOR_WRITE
      USE CMN_SIZE_MOD             ! Size Parameters
      USE ERROR_MOD,          ONLY : ALLOC_ERR
      USE inquireMod,         ONLY : findFreeLUN
      USE TIME_MOD,           ONLY : EXPAND_DATE
      USE TIME_MOD,           ONLY : GET_NYMD_DIAG    
      USE TIME_MOD,           ONLY : GET_NHMS
      USE TIME_MOD,           ONLY : GET_TAU 
      USE TIME_MOD,           ONLY : TIMESTAMP_STRING

#if defined( USE_HDF5 )
      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      !%%%%% NOTE: THIS IS PROBABLY BROKEN, WE WILL NOT REPLACE IT
      !%%%%% BECAUSE WE ARE GOING TO MOVE TO netCDF OUTPUT (bmy, 4/11/18)
      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      ! Only include this if we are linking to HDF5 library (bmy, 12/21/09)
      USE HDF_MOD,            ONLY : OPEN_HDF
      USE HDF_MOD,            ONLY : CLOSE_HDF
      USE HDF_MOD,            ONLY : WRITE_HDF
      USE HDF5,               ONLY : HID_T
      INTEGER(HID_T)              :: IU_ND51b_HDF
#endif

#endif
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)    :: am_I_Root   ! Are we on the root CPU?
      TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
      REAL(fp),       INTENT(IN)    :: TAU_W       ! TAU value at time of write
!
! !INPUT/OUTPUT PARAMETERS:
!
      TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) TAU_W (REAL(fp)) : TAU value at time of writing to disk 
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
#if defined( BPCH_DIAG ) || defined( BPCH_TIMESER )
      INTEGER             :: I,   J,  L,  W, N, GMNL, GMTRC
      INTEGER             :: IOS, X, Y, K, NHMS, nAdvect
      CHARACTER(LEN=16)   :: STAMP
      CHARACTER(LEN=40)   :: CATEGORY
      CHARACTER(LEN=40)   :: UNIT 
      CHARACTER(LEN=255)  :: FILENAME
#endif

      !=================================================================
      ! WRITE_DIAG51b begins here!
      !=================================================================

      ! Assume success
      RC         =  GC_SUCCESS

#if defined( BPCH_DIAG ) || defined( BPCH_TIMESER )

      ! Copy values from Input_Opt
      nAdvect    = State_Chm%nAdvect

      ! Find a free file LUN
      IU_ND51b = findFreeLUN()

      ! Replace date tokens in FILENAME
      FILENAME = Input_Opt%ND51b_FILE

      ! Change to get the good timestamp: day that was run and not next 
      ! day if saved at midnight
      NHMS = GET_NHMS()
      IF ( NHMS == 0 ) NHMS = 240000

      CALL EXPAND_DATE( FILENAME, GET_NYMD_DIAG(), NHMS )
      
      ! Echo info
      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( '     - DIAG51b: Opening file ', a, ' on unit ', i4 )

      ! Open output file
      IF ( Input_Opt%LND51b_HDF ) THEN
#if   defined( USE_HDF5 )
         ! Only include this if we are linking to HDF5 library (bmy, 12/21/09)
         CALL OPEN_HDF( IU_ND51b_HDF,         FILENAME,
     &                  Input_Opt%ND51b_IMAX, Input_Opt%ND51b_IMIN, 
     &                  Input_Opt%ND51b_JMAX, Input_Opt%ND51b_JMIN,
     &                  ND51b_NI,             ND51b_NJ )
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

      DO W = 1, Input_Opt%N_ND51b
 
         ! Loop over grid boxes
         DO K = 1, ND51b_NL
         DO Y = 1, ND51b_NJ
         DO X = 1, ND51b_NI

         SELECT CASE( Input_Opt%ND51b_TRACERS(W) )

            CASE( 502 )
               !--------------------------------------------------------
               ! Avoid div by zero for tracers which are archived each
               ! chem timestep and only available in the troposphere
               !--------------------------------------------------------
               IF ( COUNT_CHEM3D(X,Y,K) > 0 ) THEN
                  Q(X,Y,K,W) = Q(X,Y,K,W) / COUNT_CHEM3D(X,Y,K)
               ELSE
                  Q(X,Y,K,W) = 0e+0_fp
               ENDIF

            CASE DEFAULT

               !--------------------------------------------------------
               ! Avoid division by zero for all other tracers
               !--------------------------------------------------------
               IF ( GOOD_CT(X) > 0 ) THEN
                  Q(X,Y,K,W) = Q(X,Y,K,W) / GOOD_CT(X) 
               ELSE
                  Q(X,Y,K,W) = 0e+0_fp
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
      DO W = 1, Input_Opt%N_ND51b

         ! ND51b tracer number
         N = Input_Opt%ND51b_TRACERS(W)

         ! Save by simulation
         IF ( N <= nAdvect ) THEN

            !---------------------
            ! GEOS-CHEM species
            !---------------------
            CATEGORY = 'IJ-AVG-$'
            UNIT     = ''              ! Let GAMAP pick unit
            GMNL     = ND51b_NL
            GMTRC    = N

         ELSE IF ( N == 501 ) THEN

            !---------------------
            ! OH 
            !---------------------
            CATEGORY  = 'CHEM-L=$'
            UNIT      = 'molec/cm3'
            GMNL      = ND51b_NL
            GMTRC     = 1

         ELSE IF ( N == 502 ) THEN

            !---------------------
            ! NOy 
            !---------------------
            CATEGORY = 'TIME-SER'
            UNIT     = ''              ! Let GAMAP pick unit
            GMNL     = ND51b_NL
            GMTRC    = 2

         ELSE IF ( N == 503 ) THEN

            !---------------------
            ! Relative humidity 
            !---------------------            
            CATEGORY = 'TIME-SER'
            UNIT     = '%'
            GMNL     = ND51b_NL
            GMTRC    = 3

         ELSE IF ( N == 504 ) THEN

            !---------------------
            ! 3-D Cloud fractions
            !---------------------
            CATEGORY = 'TIME-SER'
            UNIT     = 'unitless'
            GMNL     = ND51b_NL
            GMTRC    = 4

         ELSE IF ( N == 505 ) THEN

            !---------------------
            ! Column opt depths 
            !---------------------
            CATEGORY = 'TIME-SER'
            UNIT     = 'unitless'
            GMNL     = 1
            GMTRC    = 5
            
         ELSE IF ( N == 506 ) THEN
        
            !---------------------
            ! Cloud top heights 
            !---------------------
            CATEGORY = 'TIME-SER'
            UNIT     = 'hPa'
            GMNL     = 1
            GMTRC    = 6

         ELSE IF ( N == 507 ) THEN

            !---------------------
            ! Air Density 
            !---------------------
            CATEGORY = 'TIME-SER'
            UNIT     = 'molec/cm3'
            GMNL     = ND51b_NL
            GMTRC    = 7

         ELSE IF ( N == 508 ) THEN

            !---------------------
            ! Total seasalt
            !---------------------            
            CATEGORY = 'TIME-SER'
            UNIT     = ''              ! Let GAMAP pick unit
            GMNL     = ND51b_NL
            GMTRC    = 8

         ELSE IF ( N == 509 ) THEN 

            !---------------------
            ! PBL Height [m] 
            !---------------------
            CATEGORY = 'PBLDEPTH'
            UNIT     = 'm'
            GMNL     = 1
            GMTRC    = 1

         ELSE IF ( N == 510 ) THEN

            !---------------------
            ! PBL Height [levels]
            !---------------------
            CATEGORY = 'PBLDEPTH'
            UNIT     = 'levels'
            GMNL     = 1
            GMTRC    = 2

         ELSE IF ( N == 511 ) THEN

            !---------------------
            ! Grid box heights
            !---------------------            
            CATEGORY = 'BXHGHT-$'
            UNIT     = 'm'
            GMNL     = ND51b_NL
            GMTRC    = 1

         ELSE IF ( N == 512 ) THEN

            !---------------------
            ! PEDGE-$
            !---------------------
            CATEGORY = 'PEDGE-$'
            UNIT     = 'hPa'
            GMNL     = ND51b_NL
            GMTRC    = 1

         ELSE IF ( N == 513 ) THEN

            !---------------------
            ! Sea level prs
            !---------------------            
            CATEGORY = 'DAO-FLDS'
            UNIT     = 'hPa'
            GMNL     = 1
            GMTRC    = 18

         ELSE IF ( N == 514 ) THEN

            !---------------------
            ! U-wind
            !---------------------            
            CATEGORY = 'DAO-3D-$'
            UNIT     = 'm/s'
            GMNL     = ND51b_NL
            GMTRC    = 1

         ELSE IF ( N == 515 ) THEN

            !---------------------
            ! V-wind
            !---------------------
            CATEGORY = 'DAO-3D-$'
            UNIT     = 'm/s'
            GMNL     = ND51b_NL
            GMTRC    = 2

         ELSE IF ( N == 516 ) THEN

            !---------------------
            ! Temperature
            !---------------------
            CATEGORY = 'DAO-3D-$'
            UNIT     = 'K'
            GMNL     = ND51b_NL
            GMTRC    = 3


         ELSE IF ( N == 517 ) THEN

            !---------------------
            ! Sulfate AOD
            !---------------------            
            CATEGORY = 'OD-MAP-$'
            UNIT     = 'unitless'
            GMNL     = ND51b_NL
            GMTRC    = 6

         ELSE IF ( N == 518 ) THEN

            !---------------------
            ! Black Carbon AOD
            !---------------------            
            CATEGORY = 'OD-MAP-$'
            UNIT     = 'unitless'
            GMNL     = ND51b_NL
            GMTRC    = 9

         ELSE IF ( N == 519 ) THEN

            !---------------------
            ! Organic Carbon AOD
            !---------------------            
            CATEGORY = 'OD-MAP-$'
            UNIT     = 'unitless'
            GMNL     = ND51b_NL
            GMTRC    = 12
            
         ELSE IF ( N == 520 ) THEN

            !---------------------
            ! SS Accum AOD
            !---------------------            
            CATEGORY = 'OD-MAP-$'
            UNIT     = 'unitless'
            GMNL     = ND51b_NL
            GMTRC    = 15

         ELSE IF ( N == 521 ) THEN

            !---------------------
            ! SS Coarse AOD
            !---------------------            
            CATEGORY = 'OD-MAP-$'
            UNIT     = 'unitless'
            GMNL     = ND51b_NL
            GMTRC    = 18

         ELSE IF ( N == 522 ) THEN

            !---------------------
            ! Total dust OD
            !---------------------   
            CATEGORY = 'OD-MAP-$'
            UNIT     = 'unitless'
            GMNL     = ND51b_NL
            GMTRC    = 4

         ELSE IF ( ( N >= 523 ) .and. ( N <= 529) ) THEN

            !---------------------
            ! dust OD (bins 1-7)
            !---------------------
            CATEGORY = 'OD-MAP-$'
            UNIT     = 'unitless'
            GMNL     = ND51b_NL
            GMTRC    = 21+(N-523)

! ================================================================
! Added with MEGAN v2.1. (ccc, 11/20/09)

         ELSE IF ( N == 530 ) THEN
            
            !---------------------
            ! PARDR [W/m2]
            ! (mpb,2009)
            !---------------------
            CATEGORY = 'DAO-FLDS' 
            UNIT     = 'W/m2'    
            GMNL     = ND51b_NL
            GMTRC    = 20

         ELSE IF ( N == 531 ) THEN
            
            !---------------------
            ! PARDF [W/m2]
            ! (mpb,2009)
            !---------------------
            CATEGORY = 'DAO-FLDS' 
            UNIT     = 'W/m2'    
            GMNL     = ND51b_NL
            GMTRC    = 21

         ELSE IF ( N == 532 ) THEN

            !---------------------
            ! DAILY LAI [W/m2]
            ! (mpb,2009)
            !---------------------
            CATEGORY = 'TIME-SER' 
            UNIT     = 'm2/m2'    
            GMNL     = ND51b_NL
            GMTRC    = 9

         ELSE IF ( N == 533 ) THEN

            !---------------------
            ! T at 2m
            ! (mpb,2008)
            !---------------------
            CATEGORY = 'DAO-FLDS'
            UNIT     = 'K'     
            GMNL     = ND51b_NL
            GMTRC    = 5

         ELSE

            ! Otherwise skip
            CYCLE

         ENDIF

         !------------------------
         ! Save to bpch file
         !------------------------
        IF ( Input_Opt%LND51b_HDF ) THEN
#if   defined( USE_HDF5 )
            ! Only include this if we are linking to HDF5 library 
            ! (bmy, 12/21/09)
            CALL WRITE_HDF( IU_ND51b_HDF, N,
     &                  CATEGORY,     GMTRC,        UNIT,
     &                  TAU0,         TAU1,         RESERVED,
     &                  ND51b_NI,     ND51b_NJ,     GMNL,
     &                  Input_Opt%ND51b_IMIN+I0,
     &                  Input_Opt%ND51b_JMIN+J0,
     &                  Input_Opt%ND51b_LMIN,
     &                  REAL( Q(1:ND51b_NI, 1:ND51b_NJ, 1:GMNL, W)))
#else
         ELSE
            CALL BPCH2( IU_ND51b,     MODELNAME,    LONRES,   
     &                  LATRES,       HALFPOLAR,    CENTER180, 
     &                  CATEGORY,     GMTRC,        UNIT,      
     &                  TAU0,         TAU1,         RESERVED,  
     &                  ND51b_NI,     ND51b_NJ,     GMNL,     
     &                  Input_Opt%ND51b_IMIN+I0,
     &                  Input_Opt%ND51b_JMIN+J0,
     &                  Input_Opt%ND51b_LMIN, 
     &                  REAL( Q(1:ND51b_NI, 1:ND51b_NJ, 1:GMNL, W),4))
#endif
         ENDIF
      ENDDO

      ! Echo info
      WRITE( 6, 120 ) TRIM( FILENAME )
 120  FORMAT( '     - DIAG51b: Closing file ', a )

      ! Close file
      IF ( Input_Opt%LND51b_HDF ) THEN
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
      Q            = 0e+0_fp

      ! Zero counter arrays
      COUNT_CHEM3D = 0e+0_fp
      GOOD_CT      = 0e+0_fp
#endif

      END SUBROUTINE WRITE_DIAG51b
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
#if defined( BPCH_DIAG ) || defined( BPCH_TIMESER )
      ! Add the offset to X to get I  
      I = IOFF + X

      ! Handle wrapping around the date line, if necessary
      IF ( I > IIPAR ) I = I - IIPAR
#endif

      END FUNCTION GET_I
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_diag51b
!
! !DESCRIPTION: Subroutine INIT\_DIAG51b allocates and zeroes all module 
!  arrays.  It also gets values for module variables from "input\_mod.f".
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE INIT_DIAG51b( am_I_Root, Input_Opt, RC )
!
! !USES:
!
      USE ErrCode_MOd
      USE Input_Opt_Mod, ONLY : OptInput
#if defined( BPCH_DIAG ) || defined( BPCH_TIMESER )
      USE BPCH2_MOD,     ONLY : GET_MODELNAME
      USE BPCH2_MOD,     ONLY : GET_HALFPOLAR
      USE CMN_SIZE_MOD        ! Size parameters
      USE ErrCode_Mod
      USE ERROR_MOD,     ONLY : ALLOC_ERR
      USE ERROR_MOD,     ONLY : ERROR_STOP
      USE GC_GRID_MOD,   ONLY : GET_XOFFSET
      USE GC_GRID_MOD,   ONLY : GET_YOFFSET
      USE GC_GRID_MOD,   ONLY : ITS_A_NESTED_GRID
      USE TIME_MOD,      ONLY : GET_TAUb
#endif
!
! !INPUT PARAMETERS: 
!
      LOGICAL,        INTENT(IN)  :: am_I_Root   ! Are we on the root CPU?
      TYPE(OptInput), INTENT(IN)  :: Input_Opt   ! Input Options object
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT) :: RC          ! Success or failure?
!
! !REVISION HISTORY: 
!  20 Jul 2004 - R. Yantosca - Initial version
!  (1 ) Diagnostic counter arrays are now only 1-D.  Also add GOOD_CT_CHEM
!        which is the counter array of "good" boxes at each chemistry
!        timesteps.  Now allocate GOOD_CT_CHEM. (bmy, 10/25/04)
!  (2 ) Now get I0 and J0 correctly for nested grid simulations (bmy, 11/9/04)
!  (3 ) Now call GET_HALFPOLAR from "bpch2_mod.f" to get the HALFPOLAR flag 
!        value for GEOS or GCAP grids. (bmy, 6/28/05)
!  (4 ) Now allow ND51b_IMIN to be equal to ND51b_IMAX and ND51b_JMIN to be
!        equal to ND51b_JMAX.  This will allow us to save out longitude or
!        latitude transects.  Allocate COUNT_CHEM3D. (cdh, bmy, phs, 1/24/07)
!  (5 ) Allocate GOOD_EMIS and GOOD_CT_EMIS (ccc, 12/12/09)
!  02 Dec 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
#if defined( BPCH_DIAG ) || defined( BPCH_TIMESER )
      INTEGER            :: AS
      CHARACTER(LEN=255) :: LOCATION
#endif
      
      !=================================================================
      ! INIT_DIAG51b begins here!
      !=================================================================

      ! Assume success
      RC = GC_SUCCESS

#if defined( BPCH_DIAG ) || defined( BPCH_TIMESER )

      ! Initialize
      LOCATION               = 'INIT_DIAG51b ("diag51b_mod.f")'

      ! Exit if ND51b is turned off 
      IF ( .not. Input_Opt%DO_ND51b ) RETURN

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

      ! Error check ND51b_IMIN
      IF ( Input_Opt%ND51b_IMIN+I0 < 1      .or.
     &     Input_Opt%ND51b_IMIN+I0 > IIPAR ) THEN
         CALL ERROR_STOP( 'Bad ND51b_IMIN value!', LOCATION )
      ENDIF

      ! Error check ND51b_IMAX
      IF ( Input_Opt%ND51b_IMAX+I0 < 1       .or.
     &     Input_Opt%ND51b_IMAX+I0 > IIPAR ) THEN
         CALL ERROR_STOP( 'Bad ND51b_IMAX value!', LOCATION )
      ENDIF

      ! Compute longitude limits to write to disk
      ! Also handle wrapping around the date line
      IF ( Input_Opt%ND51b_IMAX >= Input_Opt%ND51b_IMIN ) THEN
         ND51b_NI = ( Input_Opt%ND51b_IMAX - Input_Opt%ND51b_IMIN ) + 1
      ELSE 
         ND51b_NI = ( IIPAR - Input_Opt%ND51b_IMIN ) + 1 +
     &               Input_Opt%ND51b_IMAX
         WRITE( 6, '(a)' ) 'We are wrapping over the date line!'
      ENDIF

      ! Make sure that ND50_NI <= IIPAR
      IF ( ND51b_NI > IIPAR ) THEN
         CALL ERROR_STOP( 'Too many longitudes!', LOCATION )
      ENDIF

      !-----------
      ! Latitude
      !-----------
      
      ! Error check JMIN_AREA
      IF ( Input_Opt%ND51b_JMIN+J0 < 1       .or.
     &     Input_Opt%ND51b_JMIN+J0 > JJPAR ) THEN
         CALL ERROR_STOP( 'Bad ND51b_JMIN value!', LOCATION )
      ENDIF
     
      ! Error check JMAX_AREA
      IF ( Input_Opt%ND51b_JMAX+J0 < 1       .or.
     &     Input_Opt%ND51b_JMAX+J0 > JJPAR ) THEN
         CALL ERROR_STOP( 'Bad ND51b_JMAX value!', LOCATION )
      ENDIF

      ! Compute latitude limits to write to disk (bey, bmy, 3/16/99)
      IF ( Input_Opt%ND51b_JMAX >= Input_Opt%ND51b_JMIN ) THEN
         ND51b_NJ = ( Input_Opt%ND51b_JMAX - Input_Opt%ND51b_JMIN ) + 1
      ELSE
         CALL ERROR_STOP( 'ND51b_JMAX < ND51b_JMIN!', LOCATION )
      ENDIF     
  
      !-----------
      ! Altitude
      !-----------

      ! Error check ND51b_LMIN, ND51b_LMAX
      IF ( Input_Opt%ND51b_LMIN < 1       .or.
     &     Input_Opt%ND51b_LMAX > LLPAR ) THEN 
         CALL ERROR_STOP( 'Bad ND51b altitude values!', LOCATION )
      ENDIF

      ! # of levels to save in ND51b timeseries
      IF ( Input_Opt%ND51b_LMAX >= Input_Opt%ND51b_LMIN ) THEN  
         ND51b_NL = ( Input_Opt%ND51b_LMAX - Input_Opt%ND51b_LMIN ) + 1
      ELSE
         CALL ERROR_STOP( 'ND51b_LMAX < ND51b_LMIN!', LOCATION )
      ENDIF

      !-----------
      ! Offsets
      !-----------
      IOFF      = Input_Opt%ND51b_IMIN - 1
      JOFF      = Input_Opt%ND51b_JMIN - 1
      LOFF      = Input_Opt%ND51b_LMIN - 1

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
      ALLOCATE( GOOD_CT( ND51b_NI ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'GOOD_CT' )
      GOOD_CT = 0

      ! Accumulating array
      ALLOCATE( Q( ND51b_NI, ND51b_NJ, ND51b_NL, Input_Opt%N_ND51b ),
     &             STAT=AS)
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'Q' )
      Q = 0e+0_fp

      ! Accumulating array
      ALLOCATE( COUNT_CHEM3D( ND51b_NI, ND51b_NJ, ND51b_NL ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'COUNT_CHEM3D' )
      COUNT_CHEM3D = 0
#endif

      END SUBROUTINE INIT_DIAG51b
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cleanup_diag51b
!
! !DESCRIPTION: Subroutine CLEANUP\_DIAG51b deallocates all module arrays. 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CLEANUP_DIAG51b()
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
#if defined( BPCH_DIAG ) || defined( BPCH_TIMESER )
      !=================================================================
      ! CLEANUP_DIAG51b begins here!
      !=================================================================
      IF ( ALLOCATED( COUNT_CHEM3D ) ) DEALLOCATE( COUNT_CHEM3D )
      IF ( ALLOCATED( GOOD         ) ) DEALLOCATE( GOOD         )
      IF ( ALLOCATED( GOOD_CT      ) ) DEALLOCATE( GOOD_CT      )
      IF ( ALLOCATED( Q            ) ) DEALLOCATE( Q            )
#endif

      END SUBROUTINE CLEANUP_DIAG51b
!EOC
      END MODULE DIAG51b_MOD
