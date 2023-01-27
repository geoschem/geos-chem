!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: input_opt_mod.F90
!
! !DESCRIPTION: Module INPUT\_OPT\_MOD contains the derived type for GEOS-Chem
!  options and logical switches.
!\\
!\\
! !INTERFACE:
!
MODULE Input_Opt_Mod
!
! !USES:
!
  USE PRECISION_MOD    ! For GEOS-Chem Precision (fp)

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: Set_Input_Opt
  PUBLIC :: Set_Input_Opt_Advect
  PUBLIC :: Cleanup_Input_Opt
!
! !PUBLIC DATA MEMBERS:
!
  !=========================================================================
  ! Derived type for Input Options
  !=========================================================================
  TYPE, PUBLIC :: OptInput

     !----------------------------------------
     ! General Runtime & Distributed Comp Info
     !----------------------------------------
     INTEGER                     :: numCPUs    ! Number of MPI procs
     INTEGER                     :: thisCPU    ! Local MPI process handle
     INTEGER                     :: MPIComm    ! MPI Communicator Handle
     LOGICAL                     :: isMPI      ! Is this an MPI sim?
     LOGICAL                     :: amIRoot    ! Is this the root cpu?

     !----------------------------------------
     ! Dry run info (print out file names)
     !----------------------------------------
     LOGICAL                     :: DryRun     ! Is this a dry run?

     !----------------------------------------
     ! SIZE PARAMETER fields
     !----------------------------------------
     INTEGER                     :: Max_BPCH_Diag
     INTEGER                     :: Max_Families
     INTEGER                     :: Max_AdvectSpc
     INTEGER                     :: Max_PassiveSpc

     !----------------------------------------
     ! SIMULATION MENU fields
     !----------------------------------------
     INTEGER                     :: NYMDb
     INTEGER                     :: NHMSb
     INTEGER                     :: NYMDe
     INTEGER                     :: NHMSe
     INTEGER                     :: SimLengthSec
     CHARACTER(LEN=255)          :: RUN_DIR
     CHARACTER(LEN=255)          :: DATA_DIR
     CHARACTER(LEN=255)          :: CHEM_INPUTS_DIR
     CHARACTER(LEN=255)          :: MetField
     CHARACTER(LEN=255)          :: SimulationName
     CHARACTER(LEN=255)          :: SpcDatabaseFile
     CHARACTER(LEN=255)          :: SpcMetaDataOutFile
     LOGICAL                     :: ITS_A_CH4_SIM
     LOGICAL                     :: ITS_A_CO2_SIM
     LOGICAL                     :: ITS_A_FULLCHEM_SIM
     LOGICAL                     :: ITS_A_MERCURY_SIM
     LOGICAL                     :: ITS_A_POPS_SIM
     LOGICAL                     :: ITS_A_RnPbBe_SIM
     LOGICAL                     :: ITS_A_TAGO3_SIM
     LOGICAL                     :: ITS_A_TAGCO_SIM
     LOGICAL                     :: ITS_AN_AEROSOL_SIM
     LOGICAL                     :: ITS_A_TRACEMETAL_SIM
     LOGICAL                     :: ITS_A_CARBON_SIM
     LOGICAL                     :: LPRT
     LOGICAL                     :: useTimers

     !----------------------------------------
     ! PASSIVE SPECIES MENU fields
     !----------------------------------------
     INTEGER                     :: NPASSIVE
     INTEGER                     :: NPASSIVE_DECAY
     CHARACTER(LEN=63),  POINTER :: PASSIVE_NAME    (:)
     CHARACTER(LEN=255), POINTER :: PASSIVE_LONGNAME(:)
     INTEGER,            POINTER :: PASSIVE_ID      (:)
     REAL(fp),           POINTER :: PASSIVE_MW      (:)
     REAL(fp),           POINTER :: PASSIVE_TAU     (:)
     REAL(fp),           POINTER :: PASSIVE_INITCONC(:)
     INTEGER,            POINTER :: PASSIVE_DECAYID (:)

     !----------------------------------------
     ! ADVECTED SPECIES MENU fields
     !----------------------------------------
     INTEGER                     :: N_ADVECT
     CHARACTER(LEN=255), POINTER :: AdvectSpc_Name(:)
     LOGICAL                     :: LSPLIT

     !----------------------------------------
     ! AEROSOL MENU fields
     !----------------------------------------
     LOGICAL                     :: LSULF
     LOGICAL                     :: LMETALCATSO2
     LOGICAL                     :: LCARB
     LOGICAL                     :: LBRC
     LOGICAL                     :: LSOA
     LOGICAL                     :: LMPOA
     LOGICAL                     :: LSVPOA
     LOGICAL                     :: LDUST
     LOGICAL                     :: LDEAD
     LOGICAL                     :: LSSALT
     LOGICAL                     :: LDSTUP
     REAL(fp),           POINTER :: SALA_REDGE_um(:)
     REAL(fp),           POINTER :: SALC_REDGE_um(:)
     LOGICAL                     :: LGRAVSTRAT
     LOGICAL                     :: LSOLIDPSC
     LOGICAL                     :: LHOMNUCNAT
     REAL(fp)                    :: T_NAT_SUPERCOOL
     REAL(fp)                    :: P_ICE_SUPERSAT
     LOGICAL                     :: LPSCCHEM
     LOGICAL                     :: LSTRATOD
     !for BC absorption enhancement, (xnw, 8/24/15)
     LOGICAL                     :: LBCAE
     REAL(fp)                    :: BCAE_1
     REAL(fp)                    :: BCAE_2
     ! for nitrate aerosol photolysis (TMS, 23/08/2018)
     LOGICAL                     :: hvAerNIT
     REAL(fp)                    :: hvAerNIT_JNIT
     REAL(fp)                    :: hvAerNIT_JNITs
     REAL(fp)                    :: JNITChanA
     REAL(fp)                    :: JNITChanB

     !----------------------------------------
     ! EMISSIONS fields
     !----------------------------------------
     LOGICAL                     :: DoEmissions
     INTEGER                     :: TS_EMIS
     LOGICAL                     :: LBIOFUEL
     LOGICAL                     :: LOTDLOC
     LOGICAL                     :: LSOILNOX
     LOGICAL                     :: LCH4SBC
     LOGICAL                     :: LSETH2O
     LOGICAL                     :: LStaticH2OBC
     LOGICAL                     :: LHCodedOrgHal
     LOGICAL                     :: LCMIP6OrgHal
     LOGICAL                     :: DoLightNOx ! Shadow for LightNOX extension

     ! For HEMCO "intermediate" grid (hplin, 6/2/20)
     LOGICAL                     :: LIMGRID    ! Use different grid resolution for HEMCO?
     INTEGER                     :: IMGRID_XSCALE
     INTEGER                     :: IMGRID_YSCALE

     !----------------------------------------
     ! CO MENU fields
     !----------------------------------------
     LOGICAL                     :: LPCO_CH4
     LOGICAL                     :: LPCO_NMVOC

     !----------------------------------------
     ! CO2 MENU fields
     !----------------------------------------
     LOGICAL                     :: LFOSSIL
     LOGICAL                     :: LCHEMCO2
     LOGICAL                     :: LBIODIURNAL
     LOGICAL                     :: LBIONETCLIM
     LOGICAL                     :: LOCEAN
     LOGICAL                     :: LSHIP
     LOGICAL                     :: LPLANE
     LOGICAL                     :: LFFBKGRD
     LOGICAL                     :: LBIOSPHTAG
     LOGICAL                     :: LFOSSILTAG
     LOGICAL                     :: LSHIPTAG
     LOGICAL                     :: LPLANETAG

     !----------------------------------------
     ! CHEMISTRY MENU fields
     !----------------------------------------
     LOGICAL                     :: LCHEM
     LOGICAL                     :: LINEAR_CHEM
     LOGICAL                     :: LLINOZ
     LOGICAL                     :: LSYNOZ
     INTEGER                     :: TS_CHEM
     REAL(fp)                    :: GAMMA_HO2
     LOGICAL                     :: LACTIVEH2O
     LOGICAL                     :: LINITSPEC
     LOGICAL                     :: USE_ONLINE_O3
     LOGICAL                     :: USE_O3_FROM_MET
     LOGICAL                     :: USE_TOMS_O3
     LOGICAL                     :: USE_AUTOREDUCE
     LOGICAL                     :: AUTOREDUCE_IS_KEEPACTIVE
     LOGICAL                     :: AUTOREDUCE_IS_KEY_THRESHOLD
     LOGICAL                     :: AUTOREDUCE_IS_PRS_THRESHOLD
     LOGICAL                     :: AUTOREDUCE_IS_APPEND
     REAL(f8)                    :: AUTOREDUCE_THRESHOLD
     REAL(f8)                    :: AUTOREDUCE_TUNING_OH
     REAL(f8)                    :: AUTOREDUCE_TUNING_NO2
#ifdef MODEL_GEOS
     LOGICAL                     :: LGMIOZ
#endif

     !----------------------------------------
     ! PHOTOLYSIS MENU fields
     !----------------------------------------
     CHARACTER(LEN=255)          :: FAST_JX_DIR

     !----------------------------------------
     ! RADIATION MENU fields
     !----------------------------------------
     LOGICAL                     :: LRAD
     LOGICAL                     :: LLWRAD
     LOGICAL                     :: LSWRAD
     LOGICAL,            POINTER :: LSKYRAD(:)
     INTEGER                     :: TS_RAD
     INTEGER                     :: NWVSELECT
     REAL(8),            POINTER :: WVSELECT(:)
     CHARACTER(LEN=5),   POINTER :: STRWVSELECT(:)
     INTEGER                     :: NSPECRADMENU
     INTEGER,            POINTER :: LSPECRADMENU(:)

     !----------------------------------------
     ! TRANSPORT MENU fields
     !----------------------------------------
     LOGICAL                     :: LTRAN
     LOGICAL                     :: LFILL
     INTEGER                     :: TPCORE_IORD
     INTEGER                     :: TPCORE_JORD
     INTEGER                     :: TPCORE_KORD
     INTEGER                     :: TS_DYN

     !----------------------------------------
     ! CONVECTION MENU fields
     !----------------------------------------
     LOGICAL                     :: LCONV
     LOGICAL                     :: LTURB
     LOGICAL                     :: LNLPBL
     INTEGER                     :: TS_CONV

     !----------------------------------------
     ! DEPOSITION MENU fields
     !----------------------------------------
     LOGICAL                     :: LDRYD
     LOGICAL                     :: LWETD
     REAL(fp)                    :: WETD_CONV_SCAL
     LOGICAL                     :: PBL_DRYDEP
     LOGICAL                     :: CO2_EFFECT
     REAL(fp)                    :: CO2_LEVEL
     REAL(fp)                    :: CO2_REF
     REAL(fp)                    :: RS_SCALE
     INTEGER                     :: RA_Alt_Above_Sfc

     !----------------------------------------
     ! GAMAP MENU fields
     !----------------------------------------
     CHARACTER(LEN=255)          :: GAMAP_DIAGINFO
     CHARACTER(LEN=255)          :: GAMAP_TRACERINFO

     !----------------------------------------
     ! OUTPUT MENU fields
     !----------------------------------------
     INTEGER,            POINTER :: NJDAY(:)

     !----------------------------------------
     ! DIAGNOSTIC MENU fields
     !----------------------------------------
     CHARACTER(LEN=255)          :: HistoryInputFile
     INTEGER                     :: ND03   ! Hg
     INTEGER                     :: ND06   ! TOMAS
     INTEGER                     :: ND44   ! TOMAS
     INTEGER                     :: ND53   ! POPs
     INTEGER                     :: ND59   ! TOMAS
     INTEGER                     :: ND60   ! TOMAS
     INTEGER                     :: ND61   ! TOMAS

     INTEGER                     :: TS_DIAG
     INTEGER,            POINTER :: TINDEX(:,:)
     INTEGER,            POINTER :: TCOUNT(:)
     INTEGER,            POINTER :: TMAX(:)
     LOGICAL                     :: DO_DIAG_WRITE

     ! Collection ids
     INTEGER                     :: DIAG_COLLECTION
     INTEGER                     :: GC_RST_COLLECTION ! Used only for NetCDF

     !----------------------------------------
     ! PLANEFLIGHT MENU fields
     !----------------------------------------
     LOGICAL                     :: Do_Planeflight
     CHARACTER(LEN=255)          :: Planeflight_InFile
     CHARACTER(LEN=255)          :: Planeflight_OutFile

     !----------------------------------------
     ! OBSPACK MENU fields
     !----------------------------------------
     LOGICAL                     :: Do_ObsPack
     LOGICAL                     :: ObsPack_Quiet
     CHARACTER(LEN=255)          :: ObsPack_InputFile
     CHARACTER(LEN=255)          :: ObsPack_OutputFile
     INTEGER                     :: ObsPack_nSpc
     CHARACTER(LEN=255), POINTER :: ObsPack_SpcName(:)

     !----------------------------------------
     ! ND51 MENU fields
     !----------------------------------------
     LOGICAL                     :: DO_ND51
     INTEGER                     :: N_ND51
     CHARACTER(LEN=255)          :: ND51_FILE
     INTEGER,            POINTER :: ND51_TRACERS(:)
     REAL(fp)                    :: ND51_HR_WRITE
     REAL(fp)                    :: ND51_HR1
     REAL(fp)                    :: ND51_HR2
     INTEGER                     :: ND51_IMIN
     INTEGER                     :: ND51_IMAX
     INTEGER                     :: ND51_JMIN
     INTEGER                     :: ND51_JMAX
     INTEGER                     :: ND51_LMIN
     INTEGER                     :: ND51_LMAX

     !----------------------------------------
     ! ND51b MENU fields
     !----------------------------------------
     LOGICAL                     :: DO_ND51b
     INTEGER                     :: N_ND51b
     CHARACTER(LEN=255)          :: ND51b_FILE
     INTEGER,            POINTER :: ND51b_TRACERS(:)
     REAL(fp)                    :: ND51b_HR_WRITE
     REAL(fp)                    :: ND51b_HR1
     REAL(fp)                    :: ND51b_HR2
     INTEGER                     :: ND51b_IMIN
     INTEGER                     :: ND51b_IMAX
     INTEGER                     :: ND51b_JMIN
     INTEGER                     :: ND51b_JMAX
     INTEGER                     :: ND51b_LMIN
     INTEGER                     :: ND51b_LMAX

     !----------------------------------------
     ! PROD LOSS MENU fields
     !----------------------------------------
     LOGICAL                     :: DO_SAVE_PL
     INTEGER                     :: ND65, LD65
     INTEGER                     :: NFAM
     CHARACTER(LEN=255), POINTER :: FAM_NAME(:)
     CHARACTER(LEN=255), POINTER :: FAM_TYPE(:)

     !----------------------------------------
     ! BENCHMARK MENU fields
     !----------------------------------------
     LOGICAL                     :: LSTDRUN
     CHARACTER(LEN=255)          :: STDRUN_INIT_FILE
     CHARACTER(LEN=255)          :: STDRUN_FINAL_FILE

     !----------------------------------------
     ! MERCURY MENU fields
     !----------------------------------------
     INTEGER                     :: ANTHRO_Hg_YEAR
     CHARACTER(LEN=255)          :: HG_SCENARIO
     LOGICAL                     :: USE_CHECKS
     LOGICAL                     :: LDYNOCEAN
     LOGICAL                     :: LPREINDHG
     LOGICAL                     :: LGTMM
     CHARACTER(LEN=255)          :: GTMM_RST_FILE
     LOGICAL                     :: LARCTICRIV
     LOGICAL                     :: LKRedUV

     !----------------------------------------
     ! CH4 MENU fields
     !----------------------------------------
     LOGICAL                     :: GOSAT_CH4_OBS
     LOGICAL                     :: AIRS_CH4_OBS
     LOGICAL                     :: TCCON_CH4_OBS
     LOGICAL                     :: AnalyticalInv
     REAL(fp)                    :: PerturbEmis
     INTEGER                     :: StateVectorElement
     LOGICAL                     :: UseEmisSF
     LOGICAL                     :: UseOHSF

     !----------------------------------------
     ! POPS MENU fields
     !----------------------------------------
     CHARACTER(LEN=3)            :: POP_TYPE
     LOGICAL                     :: CHEM_PROCESS
     REAL(fp)                    :: POP_XMW
     REAL(fp)                    :: POP_KOA
     REAL(fp)                    :: POP_KBC
     REAL(fp)                    :: POP_K_POPG_OH
     REAL(fp)                    :: POP_K_POPP_O3A
     REAL(fp)                    :: POP_K_POPP_O3B
     REAL(fp)                    :: POP_HSTAR
     REAL(fp)                    :: POP_DEL_H
     REAL(fp)                    :: POP_DEL_Hw

     !----------------------------------------
     ! Fields for interface to GEOS-5 GCM
     !----------------------------------------
#ifdef MODEL_GEOS
     LOGICAL                     :: LCAPTROP     = .FALSE.
     !REAL(fp)                    :: OZONOPAUSE   = -999.0
     LOGICAL                     :: haveImpRst   = .FALSE.
     LOGICAL                     :: AlwaysSetH2O = .TRUE.
     LOGICAL                     :: UseOnlineVUD = .FALSE.
     INTEGER                     :: LLFASTJX     = 601
     INTEGER                     :: NN_RxnRates             ! # of diagnosed reaction rates
     INTEGER, POINTER            :: RxnRates_IDs(:)         ! Reaction rate numbers to be diagnosed
     INTEGER                     :: NN_RxnRconst            ! # of diagnosed reaction rates
     INTEGER, POINTER            :: RxnRconst_IDs(:)        ! Reaction rate numbers to be diagnosed
     INTEGER                     :: NN_Jvals                ! # of diagnosed Jvalues
     INTEGER, POINTER            :: Jval_IDs(:)             ! J-values to be diagnosed
     INTEGER                     :: FJX_EXTRAL_ITERMAX = 5
     LOGICAL                     :: FJX_EXTRAL_ERR     = .TRUE.
     ! Toggle for het rates. If true, turns off three Cl producing het reactions
     ! in the stratosphere. In MODEL_GEOS, this flag is set in GEOSCHEMchem_GridComp.rc
     LOGICAL                     :: TurnOffHetRates = .FALSE.
#else
     LOGICAL                     :: AlwaysSetH2O
     LOGICAL                     :: TurnOffHetRates
#endif

#if defined( MODEL_GEOS ) || defined( MODEL_WRF )
     LOGICAL                     :: KppStop            = .TRUE. ! Stop KPP if integration fails twice
#endif

#if defined( MODEL_CESM )
     LOGICAL                     :: onlineAlbedo       = .TRUE. ! Use albedo from land model
     LOGICAL                     :: onlineLandTypes    = .TRUE. ! Use land types from land model
     LOGICAL                     :: ddVel_CLM          = .TRUE. ! Use dry deposition velocities as computed by the Community Land Model
     LOGICAL                     :: applyQtend         = .TRUE. ! Apply water vapor tendency to specific humidity
#endif

#ifdef ADJOINT
     !----------------------------------------
     ! GCHP adjoint fields
     !---------------------------------------
     LOGICAL                     :: IS_ADJOINT
     LOGICAL                     :: IS_FD_SPOT, IS_FD_GLOBAL
     INTEGER                     :: FD_STEP
     LOGICAL                     :: IS_FD_SPOT_THIS_PET
     INTEGER                     :: IFD, JFD, NFD, LFD, NFD_ADJ
     INTEGER                     :: CF_IMIN, CF_IMAX
     INTEGER                     :: CF_JMIN, CF_JMAX
     INTEGER                     :: CF_LMIN, CF_LMAX
#endif

     !----------------------------------------
     ! Fields for LINOZ strat chem
     !----------------------------------------
     INTEGER                     :: LINOZ_NLEVELS
     INTEGER                     :: LINOZ_NLAT
     INTEGER                     :: LINOZ_NMONTHS
     INTEGER                     :: LINOZ_NFIELDS
     REAL(fp),           POINTER :: LINOZ_TPARM(:,:,:,:)

  END TYPE OptInput
!
! !REMARKS:
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Set_Input_Opt
!
! !DESCRIPTION: Subroutine SET\_INPUT\_OPT intializes all GEOS-Chem
!  options carried in Input Options derived type object.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Set_Input_Opt( am_I_Root, Input_Opt, RC )
!
! !USES:
!
    USE ErrCode_Mod
!
! !INPUT PARAMETERS:
!
    LOGICAL,        INTENT(IN)    :: am_I_Root   ! Are we on the root CPU?
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(INOUT) :: Input_Opt   ! Input Options object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!
! !REMARKS:
!
! !REVISION HISTORY:
!  01 Nov 2012 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Strings
    CHARACTER(LEN=30) :: arrayId

    !----------------------------------------
    ! Initialize
    ! Set pointers to NULL for safety's sake
    !----------------------------------------
    RC                               =  GC_SUCCESS
    Input_Opt%PASSIVE_NAME           => NULL()
    Input_Opt%PASSIVE_ID             => NULL()
    Input_Opt%PASSIVE_MW             => NULL()
    Input_Opt%PASSIVE_TAU            => NULL()
    Input_Opt%PASSIVE_INITCONC       => NULL()
    Input_Opt%PASSIVE_DECAYID        => NULL()
    Input_Opt%AdvectSpc_Name         => NULL()
    Input_Opt%SALA_REDGE_um          => NULL()
    Input_Opt%SALC_REDGE_um          => NULL()
    Input_Opt%LSKYRAD                => NULL()
    Input_Opt%LSPECRADMENU           => NULL()
    Input_Opt%NJDAY                  => NULL()
    Input_Opt%TINDEX                 => NULL()
    Input_Opt%TCOUNT                 => NULL()
    Input_Opt%TMAX                   => NULL()
    Input_Opt%ND51_TRACERS           => NULL()
    Input_Opt%ND51b_TRACERS          => NULL()
    Input_Opt%FAM_NAME               => NULL()
    Input_Opt%FAM_TYPE               => NULL()
    Input_Opt%LINOZ_TPARM            => NULL()

    !----------------------------------------
    ! General Runtime & Distributed Comp Info
    !----------------------------------------
    Input_Opt%amIRoot                = am_I_Root
    Input_Opt%isMPI                  = .FALSE.
    Input_Opt%numCPUs                = 1
    Input_Opt%thisCPU                = -1
    Input_Opt%MPIComm                = -1

    !----------------------------------------
    ! Dry run info (print out file names)
    !----------------------------------------
    Input_Opt%DryRun                 = .FALSE.

    !----------------------------------------
    ! SIZE PARAMETER fields
    !
    ! Set to large placeholder values
    !----------------------------------------
#ifdef RRTMG
    Input_Opt%Max_BPCH_Diag          = 187 ! Mirror MAX_DIAG in CMN_DIAG_mod.F90
#else
    Input_Opt%Max_BPCH_Diag          = 80  ! Mirror MAX_DIAG in CMN_DIAG_mod.F90
#endif
    Input_Opt%Max_Families           = 250
    Input_Opt%Max_AdvectSpc          = 600
    Input_Opt%Max_PassiveSpc         = 50

    !----------------------------------------
    ! SIMULATION MENU fields
    !----------------------------------------
    Input_Opt%NYMDb                  = 0
    Input_Opt%NHMSb                  = 0
    Input_Opt%NYMDe                  = 0
    Input_Opt%NHMSe                  = 0
    Input_Opt%SimLengthSec           = 0
    Input_Opt%RUN_DIR                = './'
    Input_Opt%DATA_DIR               = './'
    Input_Opt%CHEM_INPUTS_DIR        = './'
    Input_Opt%MetField               = ''
    Input_Opt%SimulationName         = ''
    Input_Opt%SpcDatabaseFile        = ''
    Input_Opt%SpcMetaDataOutFile     = ''
    Input_Opt%ITS_A_CARBON_SIM       = .FALSE.
    Input_Opt%ITS_A_CH4_SIM          = .FALSE.
    Input_Opt%ITS_A_CO2_SIM          = .FALSE.
    Input_Opt%ITS_A_FULLCHEM_SIM     = .FALSE.
    Input_Opt%ITS_A_MERCURY_SIM      = .FALSE.
    Input_Opt%ITS_A_POPS_SIM         = .FALSE.
    Input_Opt%ITS_A_RnPbBe_SIM       = .FALSE.
    Input_Opt%ITS_A_TAGO3_SIM        = .FALSE.
    Input_Opt%ITS_A_TAGCO_SIM        = .FALSE.
    Input_Opt%ITS_AN_AEROSOL_SIM     = .FALSE.
    Input_Opt%ITS_A_TRACEMETAL_SIM   = .FALSE.
    Input_Opt%LPRT                   = .FALSE.
    Input_Opt%useTimers              = .FALSE.

    !----------------------------------------
    ! ADVECTED SPECIES MENU fields
    !----------------------------------------
    arrayId = 'Input_Opt%AdvectSpc_Name'
    ALLOCATE( Input_Opt%AdvectSpc_Name( Input_Opt%Max_AdvectSpc ), STAT=RC )
    CALL GC_CheckVar( arrayId, 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    Input_Opt%N_ADVECT               = 0
    Input_Opt%AdvectSpc_Name         = ''
    Input_Opt%LSPLIT                 = .FALSE.

    !----------------------------------------
    ! PASSIVE SPECIES MENU fields
    !----------------------------------------

    ALLOCATE( Input_Opt%PASSIVE_NAME    ( Input_Opt%Max_PassiveSpc ), STAT=RC )
    ALLOCATE( Input_Opt%PASSIVE_LONGNAME( Input_Opt%Max_PassiveSpc ), STAT=RC )
    ALLOCATE( Input_Opt%PASSIVE_ID      ( Input_Opt%Max_PassiveSpc ), STAT=RC )
    ALLOCATE( Input_Opt%PASSIVE_MW      ( Input_Opt%Max_PassiveSpc ), STAT=RC )
    ALLOCATE( Input_Opt%PASSIVE_TAU     ( Input_Opt%Max_PassiveSpc ), STAT=RC )
    ALLOCATE( Input_Opt%PASSIVE_INITCONC( Input_Opt%Max_PassiveSpc ), STAT=RC )
    ALLOCATE( Input_Opt%PASSIVE_DECAYID ( Input_Opt%Max_PassiveSpc ), STAT=RC )

    Input_Opt%NPASSIVE               = 0
    Input_Opt%NPASSIVE_DECAY         = 0
    Input_Opt%PASSIVE_NAME           = ''
    Input_Opt%PASSIVE_LONGNAME       = ''
    Input_Opt%PASSIVE_ID             = 0
    Input_Opt%PASSIVE_MW             = 0.0_fp
    Input_Opt%PASSIVE_TAU            = 0.0_fp
    Input_Opt%PASSIVE_INITCONC       = 0.0_fp
    Input_Opt%PASSIVE_DECAYID        = 0

    !----------------------------------------
    ! AEROSOL MENU fields
    !----------------------------------------
    arrayId = 'Input_Opt%SALA_REDGE_um'
    ALLOCATE( Input_Opt%SALA_REDGE_um( 2 ), STAT=RC )
    CALL GC_CheckVar( arrayId, 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    arrayId = 'Input_Opt%SALC_REDGE_um'
    ALLOCATE( Input_Opt%SALC_REDGE_um( 2 ), STAT=RC )
    CALL GC_CheckVar( arrayId, 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    Input_Opt%LSULF                  = .FALSE.
    Input_Opt%LMETALCATSO2           = .FALSE.
    Input_Opt%LCARB                  = .FALSE.
    Input_Opt%LBRC                   = .FALSE.
    Input_Opt%LSOA                   = .FALSE.
    Input_Opt%LMPOA                  = .FALSE.
    Input_Opt%LSVPOA                 = .FALSE.
    Input_Opt%LDUST                  = .FALSE.
    Input_Opt%LDEAD                  = .FALSE.
    Input_Opt%LDSTUP                 = .FALSE.
    Input_Opt%LSSALT                 = .FALSE.
    Input_Opt%SALA_REDGE_um          = 0.0_fp
    Input_Opt%SALC_REDGE_um          = 0.0_fp
    Input_Opt%LGRAVSTRAT             = .FALSE.
    Input_Opt%LSOLIDPSC              = .FALSE.
    Input_Opt%LHOMNUCNAT             = .FALSE.
    Input_Opt%T_NAT_SUPERCOOL        = 0.0_fp
    Input_Opt%P_ICE_SUPERSAT         = 0.0_fp
    Input_Opt%LPSCCHEM               = .FALSE.
    Input_Opt%LSTRATOD               = .FALSE.
    Input_Opt%hvAerNIT               = .FALSE.
    Input_Opt%hvAerNIT_JNIT          = 0.0_fp
    Input_Opt%hvAerNIT_JNITs         = 0.0_fp
    Input_Opt%JNITChanA              = 0.0_fp
    Input_Opt%JNITChanB              = 0.0_fp

    !----------------------------------------
    ! EMISSIONS MENU fields
    !----------------------------------------
    Input_Opt%DoEmissions            = .TRUE. ! On by default
    Input_Opt%TS_EMIS                = 0
    Input_Opt%LSOILNOX               = .FALSE.
    Input_Opt%LCH4SBC                = .FALSE.
    Input_Opt%LSETH2O                = .FALSE.
    Input_Opt%LStaticH2OBC           = .FALSE.
    Input_Opt%LHCodedOrgHal          = .FALSE.
    Input_Opt%LCMIP6OrgHal           = .FALSE.
    Input_Opt%DoLightNOx             = .FALSE.
    Input_Opt%LIMGRID                = .FALSE.
    Input_Opt%IMGRID_XSCALE          = 1
    Input_Opt%IMGRID_YSCALE          = 1

    !----------------------------------------
    ! CO MENU fields
    !----------------------------------------
    Input_Opt%LPCO_CH4               = .FALSE.
    Input_Opt%LPCO_NMVOC             = .FALSE.

    !----------------------------------------
    ! CO2 MENU fields
    !----------------------------------------
    Input_Opt%LFOSSIL                = .FALSE.
    Input_Opt%LCHEMCO2               = .FALSE.
    Input_Opt%LBIOFUEL               = .FALSE.
    Input_Opt%LBIODIURNAL            = .FALSE.
    Input_Opt%LBIONETCLIM            = .FALSE.
    Input_Opt%LOCEAN                 = .FALSE.
    Input_Opt%LSHIP                  = .FALSE.
    Input_Opt%LPLANE                 = .FALSE.
    Input_Opt%LFFBKGRD               = .FALSE.
    Input_Opt%LBIOSPHTAG             = .FALSE.
    Input_Opt%LFOSSILTAG             = .FALSE.
    Input_Opt%LSHIPTAG               = .FALSE.
    Input_Opt%LPLANETAG              = .FALSE.

    !----------------------------------------
    ! CHEMISTRY MENU fields
    !----------------------------------------
    Input_Opt%LCHEM                  = .FALSE.
    Input_Opt%LINEAR_CHEM            = .FALSE.
    Input_Opt%LLINOZ                 = .FALSE.
    Input_Opt%LSYNOZ                 = .FALSE.
#ifdef MODEL_GEOS
    Input_Opt%LGMIOZ                 = .FALSE.
#endif
    Input_Opt%TS_CHEM                = 0
    Input_Opt%GAMMA_HO2              = 0.0_fp
    Input_Opt%LACTIVEH2O             = .FALSE.
    Input_Opt%LINITSPEC              = .FALSE.
    Input_Opt%USE_ONLINE_O3          = .FALSE.
    Input_Opt%USE_O3_FROM_MET        = .FALSE.
    Input_Opt%USE_TOMS_O3            = .FALSE.

    Input_Opt%USE_AUTOREDUCE                = .FALSE.
    Input_Opt%AUTOREDUCE_IS_KEY_THRESHOLD   = .TRUE.
    Input_Opt%AUTOREDUCE_TUNING_OH          = 5e-5_fp
    Input_Opt%AUTOREDUCE_TUNING_NO2         = 1e-4_fp
    Input_Opt%AUTOREDUCE_IS_PRS_THRESHOLD   = .TRUE.
    Input_Opt%AUTOREDUCE_IS_KEEPACTIVE      = .FALSE.
    Input_Opt%AUTOREDUCE_IS_APPEND          = .FALSE.

    !----------------------------------------
    ! PHOTOLYSIS MENU fields
    !----------------------------------------
    Input_Opt%FAST_JX_DIR            = './'

    !----------------------------------------
    ! RADIATION MENU fields (for RRTMG only)
    !----------------------------------------
    arrayId = 'Input_Opt%LSKYRAD'
    ALLOCATE( Input_Opt%LSKYRAD( 2 ), STAT=RC )
    CALL GC_CheckVar( arrayId, 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    arrayId = 'Input_Opt%WVSELECT'
    ALLOCATE( Input_Opt%WVSELECT( 3 ), STAT=RC )
    CALL GC_CheckVar( arrayId, 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    arrayId = 'Input_Opt%STRWVSELECT'
    ALLOCATE( Input_Opt%STRWVSELECT( 3 ), STAT=RC )
    CALL GC_CheckVar( arrayId, 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    ! Number of RRTMG outputs (change as necessary)
    Input_Opt%NSpecRadMenu           = 11

    arrayId = 'Input_Opt%LSPECRADMENU'
    ALLOCATE( Input_Opt%LSPECRADMENU( Input_Opt%NSpecRadMenu ), STAT=RC )
    CALL GC_CheckVar( arrayId, 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    Input_Opt%LSpecRadMenu           = 0

    Input_Opt%LRAD                   = .FALSE.
    Input_Opt%LLWRAD                 = .FALSE.
    Input_Opt%LSWRAD                 = .FALSE.
    Input_Opt%LSKYRAD                = .FALSE.
    Input_Opt%TS_RAD                 = 0
    Input_Opt%NWVSELECT              = 0
    Input_Opt%WVSELECT               = 0.0_fp
    Input_Opt%STRWVSELECT            = ''

    !----------------------------------------
    ! TRANSPORT MENU fields
    !----------------------------------------
    Input_Opt%LTRAN                  = .FALSE.
    Input_Opt%LFILL                  = .FALSE.
    Input_Opt%TPCORE_IORD            = 0
    Input_Opt%TPCORE_JORD            = 0
    Input_Opt%TPCORE_KORD            = 0
    Input_Opt%TS_DYN                 = 0

    !----------------------------------------
    ! CONVECTION MENU fields
    !----------------------------------------
    Input_Opt%LCONV                  = .FALSE.
    Input_Opt%LTURB                  = .FALSE.
    Input_Opt%LNLPBL                 = .FALSE.
    Input_Opt%TS_CONV                = 0

    !----------------------------------------
    ! DEPOSITION MENU fields
    !----------------------------------------
    Input_Opt%LDRYD                  = .FALSE.
    Input_Opt%LWETD                  = .FALSE.
    Input_Opt%WETD_CONV_SCAL         = 1.0_fp
    Input_Opt%PBL_DRYDEP             = .FALSE.
    Input_Opt%CO2_LEVEL              = 390.0_fp
    Input_Opt%CO2_REF                = 390.0_fp
    Input_Opt%CO2_EFFECT             = .FALSE.
    Input_Opt%RS_SCALE               = 1.0_fp
    Input_Opt%RA_Alt_Above_Sfc       = 10       ! default height


    !----------------------------------------
    ! GAMAP_MENU fields
    !----------------------------------------
    Input_Opt%GAMAP_DIAGINFO         = ''
    Input_Opt%GAMAP_TRACERINFO       = ''

    !----------------------------------------
    ! OUTPUT MENU fields
    !----------------------------------------
    arrayId = 'Input_Opt%NJDAY'
    ALLOCATE( Input_Opt%NJDAY( 366 ), STAT=RC )
    CALL GC_CheckVar( arrayId, 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    Input_Opt%NJDAY                  = 0

    !----------------------------------------
    ! DIAGNOSTIC MENU fields
    !----------------------------------------
    Input_Opt%HistoryInputFile       = ''
    Input_Opt%DIAG_COLLECTION        = -999
    Input_Opt%TS_DIAG                = 0
    ALLOCATE( Input_Opt%TCOUNT( Input_Opt%Max_BPCH_Diag ), STAT=RC )
    ALLOCATE( Input_Opt%TMAX  ( Input_Opt%Max_BPCH_Diag ), STAT=RC )

    Input_Opt%ND03                   = 0
    Input_Opt%ND06                   = 0
    Input_Opt%ND44                   = 0
    Input_Opt%ND53                   = 0
    Input_Opt%ND59                   = 0
    Input_Opt%ND60                   = 0
    Input_Opt%ND61                   = 0
    Input_Opt%ND65                   = 0
    Input_Opt%TCOUNT(:)              = 0
    Input_Opt%TMAX(:)	             = 0
#if defined( ESMF_ ) || defined( EXTERNAL_GRID ) || defined( EXTERNAL_FORCING )
    ! Need to shut off G-C diagnostics when
    ! connecting to an external GCM (bmy, 3/29/13)
    Input_Opt%DO_DIAG_WRITE          = .FALSE.
#else
    ! For traditional G-C runs, always write diags (bmy, 3/29/13)
    Input_Opt%DO_DIAG_WRITE          = .TRUE.
#endif

    !----------------------------------------
    ! PLANEFLIGHT MENU fields
    !----------------------------------------
    Input_Opt%Do_Planeflight         = .FALSE.
    Input_Opt%Planeflight_InFile     = ''
    Input_Opt%Planeflight_OutFile    = ''

    !----------------------------------------
    ! PLANEFLIGHT MENU fields
    !----------------------------------------
    ALLOCATE( Input_Opt%ObsPack_SpcName( 1000 ), STAT=RC )

    Input_Opt%Do_ObsPack             = .FALSE.
    Input_Opt%ObsPack_Quiet          = .FALSE.
    Input_Opt%ObsPack_InputFile      = ''
    Input_Opt%ObsPack_OutputFile     = ''
    Input_Opt%ObsPack_nSpc           = 0
    Input_Opt%ObsPack_SpcName        = ''

    !----------------------------------------
    ! ND51 MENU fields
    !----------------------------------------
    Input_Opt%DO_ND51                = .FALSE.
    Input_Opt%N_ND51                 = 0
    Input_Opt%ND51_FILE              = ''
    Input_Opt%ND51_HR_WRITE          = 0.0_fp
    Input_Opt%ND51_HR1               = 0.0_fp
    Input_Opt%ND51_HR2               = 0.0_fp
    Input_Opt%ND51_IMIN              = 0
    Input_Opt%ND51_IMAX              = 0
    Input_Opt%ND51_JMIN              = 0
    Input_Opt%ND51_JMAX              = 0
    Input_Opt%ND51_LMIN              = 0

    !----------------------------------------
    ! ND51b MENU fields
    !----------------------------------------
    Input_Opt%DO_ND51b               = .FALSE.
    Input_Opt%N_ND51b                = 0
    Input_Opt%ND51b_FILE             = ''
    Input_Opt%ND51b_HR_WRITE         = 0.0_fp
    Input_Opt%ND51b_HR1              = 0.0_fp
    Input_Opt%ND51b_HR2              = 0.0_fp
    Input_Opt%ND51b_IMIN             = 0
    Input_Opt%ND51b_IMAX             = 0
    Input_Opt%ND51b_JMIN             = 0
    Input_Opt%ND51b_JMAX             = 0
    Input_Opt%ND51b_LMIN             = 0

    !----------------------------------------
    ! PROD LOSS MENU fields
    !---------------------------------------

    arrayId = 'Input_Opt%FAM_NAME'
    ALLOCATE( Input_Opt%FAM_NAME( Input_Opt%Max_Families ), STAT=RC )
    CALL GC_CheckVar( arrayId, 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    arrayId = 'Input_Opt%FAM_TYPE'
    ALLOCATE( Input_Opt%FAM_TYPE( Input_Opt%Max_Families ), STAT=RC )
    CALL GC_CheckVar( arrayId, 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    Input_Opt%DO_SAVE_PL             = .FALSE.
    Input_Opt%ND65                   = 0
    Input_Opt%NFAM                   = 0
    Input_Opt%FAM_NAME               = ''
    Input_Opt%FAM_TYPE               = ''

    !----------------------------------------
    ! MERCURY MENU fields
    !----------------------------------------
    Input_Opt%ANTHRO_Hg_YEAR         = 0
    Input_Opt%HG_SCENARIO            = ''
    Input_Opt%USE_CHECKS             = .FALSE.
    Input_Opt%LDYNOCEAN              = .FALSE.
    Input_Opt%LPREINDHG              = .FALSE.
    Input_Opt%LGTMM                  = .FALSE.
    Input_Opt%GTMM_RST_FILE          = ''

    !----------------------------------------
    ! CH4 MENU fields
    !----------------------------------------
    Input_Opt%GOSAT_CH4_OBS          = .FALSE.
    Input_Opt%AIRS_CH4_OBS           = .FALSE.
    Input_Opt%TCCON_CH4_OBS          = .FALSE.
    Input_Opt%AnalyticalInv          = .FALSE.
    Input_Opt%PerturbEmis            = 1.0
    Input_Opt%StateVectorElement     = 0
    Input_Opt%UseEmisSF              = .FALSE.
    Input_Opt%UseOHSF                = .FALSE.

    !----------------------------------------
    ! POPS MENU fields
    !----------------------------------------
    Input_Opt%POP_TYPE               = ''
    Input_Opt%CHEM_PROCESS           = .FALSE.
    Input_Opt%POP_XMW                = 0.0_fp
    Input_Opt%POP_KOA                = 0.0_fp
    Input_Opt%POP_KBC                = 0.0_fp
    Input_Opt%POP_K_POPG_OH          = 0.0_fp
    Input_Opt%POP_K_POPP_O3A         = 0.0_fp
    Input_Opt%POP_K_POPP_O3B         = 0.0_fp
    Input_Opt%POP_HSTAR              = 0.0_fp
    Input_Opt%POP_DEL_H              = 0.0_fp
    Input_Opt%POP_DEL_Hw             = 0.0_fp

    !----------------------------------------
    ! Fields for interface to GEOS-5 GCM
    !----------------------------------------
#ifdef MODEL_GEOS
!    Input_Opt%OZONOPAUSE             = -999.0
!    Input_Opt%haveImpRst             = .FALSE.
!    Input_Opt%AlwaysSetH2O           = .FALSE.
!    Input_Opt%LLFASTJX               = -999
    Input_Opt%NN_RxnRates            = -999
    Input_Opt%RxnRates_IDs           => NULL()
    Input_Opt%NN_RxnRconst           = -999
    Input_Opt%RxnRconst_IDs          => NULL()
    Input_Opt%NN_Jvals               = -999
    Input_Opt%Jval_IDs               => NULL()
#else
    Input_Opt%AlwaysSetH2O           = .FALSE.
    Input_Opt%TurnOffHetRates        = .FALSE.
#endif

#ifdef ADJOINT
    !----------------------------------------
    ! Fields for adoint
    !---------------------------------------
    Input_Opt%IS_ADJOINT             = .FALSE.
    Input_Opt%IS_FD_SPOT             = .FALSE.
    Input_Opt%IS_FD_GLOBAL           = .FALSE.
    Input_Opt%IS_FD_SPOT_THIS_PET    = .FALSE.
    Input_Opt%FD_STEP                = -999
    Input_Opt%IFD                    = -999
    Input_Opt%JFD                    = -999
    Input_Opt%NFD                    = -999
    Input_Opt%LFD                    = -999
#endif

    !----------------------------------------
    ! Fields for LINOZ strat chem
    !----------------------------------------
    Input_Opt%LINOZ_NLEVELS          = 25
    Input_Opt%LINOZ_NLAT             = 18
    Input_Opt%LINOZ_NMONTHS          = 12
    Input_Opt%LINOZ_NFIELDS          = 7

    arrayId = 'Input_Opt%LINOZ_TPARM'
    ALLOCATE( Input_Opt%LINOZ_TPARM( Input_Opt%LINOZ_NLEVELS,            &
                                     Input_Opt%LINOZ_NLAT,               &
                                     Input_Opt%LINOZ_NMONTHS,            &
                                     Input_Opt%LINOZ_NFIELDS ), STAT=RC )
    CALL GC_CheckVar( arrayId, 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    Input_Opt%LINOZ_TPARM            = 0.0_fp

  END SUBROUTINE Set_Input_Opt
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Set_Input_Opt_Advect
!
! !DESCRIPTION: Subroutine SET\_INPUT\_OPT\_ADVECT intializes all GEOS-Chem
!  options carried in Input Options derived type object that depend on
!  the number of advected species (Input\_Opt%N_ADVECT).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Set_Input_Opt_Advect( Input_Opt, RC )
!
! !USES:
!
    USE ErrCode_Mod
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(INOUT) :: Input_Opt   ! Input Options object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!
! !REMARKS:
!  NOTE: These arrays are all for bpch diagnostics, and will eventually
!  be removed from GEOS-Chem.

! !REVISION HISTORY:
!  26 Jan 2018 - M. Sulprizio- Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Initialize
    RC = GC_SUCCESS

    !=======================================================================
    ! Allocate arrays
    !=======================================================================

    ALLOCATE( Input_Opt%TINDEX(Input_Opt%Max_BPCH_Diag,Input_Opt%N_ADVECT), &
              STAT=RC )
    CALL GC_CheckVar( 'Input_Opt%TINDEX', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    Input_Opt%TINDEX = 0

    ALLOCATE( Input_Opt%ND51_TRACERS (Input_Opt%N_ADVECT+Input_Opt%Max_BPCH_Diag),&
              STAT=RC )
    CALL GC_CheckVar( 'Input_Opt%ND51_TRACERS', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    Input_Opt%ND51_TRACERS = 0

    ALLOCATE( Input_Opt%ND51b_TRACERS(Input_Opt%N_ADVECT+Input_Opt%Max_BPCH_Diag),&
              STAT=RC )
    CALL GC_CheckVar( 'Input_Opt%ND51b_TRACERS', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    Input_Opt%ND51b_TRACERS = 0

  END SUBROUTINE Set_Input_Opt_Advect
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Cleanup_Input_Opt
!
! !DESCRIPTION: Subroutine CLEANUP\_INPUT\_OPT deallocates all
!  allocatable fields of the Input Options object.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Cleanup_Input_Opt( Input_Opt, RC )
!
! !USES:
!
    USE ErrCode_Mod
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(INOUT) :: Input_Opt   ! Input Options object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure
!
! !REVISION HISTORY:
!  02 Nov 2012 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC

    ! Assume success
    RC = GC_SUCCESS

    !======================================================================
    ! Deallocate fields of the Input Options object
    !======================================================================
    IF ( ASSOCIATED( Input_Opt%PASSIVE_NAME ) ) THEN
       DEALLOCATE( Input_Opt%PASSIVE_NAME, STAT=RC )
       CALL GC_CheckVar( 'Input_Opt%PASSIVE_NAME', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       Input_Opt%PASSIVE_NAME => NULL()
    ENDIF

    IF ( ASSOCIATED( Input_Opt%PASSIVE_LONGNAME ) ) THEN
       DEALLOCATE( Input_Opt%PASSIVE_LONGNAME )
    ENDIF

    IF ( ASSOCIATED( Input_Opt%PASSIVE_ID ) ) THEN
       DEALLOCATE( Input_Opt%PASSIVE_ID, STAT=RC )
       CALL GC_CheckVar( 'Input_Opt%PASSIVE_ID', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       Input_Opt%PASSIVE_ID => NULL()
    ENDIF

    IF ( ASSOCIATED( Input_Opt%PASSIVE_MW ) ) THEN
       DEALLOCATE( Input_Opt%PASSIVE_MW, STAT=RC )
       CALL GC_CheckVar( 'Input_Opt%PASSIVE_MW', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       Input_Opt%PASSIVE_MW => NULL()
    ENDIF

    IF ( ASSOCIATED( Input_Opt%PASSIVE_TAU ) ) THEN
       DEALLOCATE( Input_Opt%PASSIVE_TAU, STAT=RC )
       CALL GC_CheckVar( 'Input_Opt%PASSIVE_TAU', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       Input_Opt%PASSIVE_TAU => NULL()
    ENDIF

    IF ( ASSOCIATED( Input_Opt%PASSIVE_INITCONC ) ) THEN
       DEALLOCATE( Input_Opt%PASSIVE_INITCONC, STAT=RC )
       CALL GC_CheckVar( 'Input_Opt%PASSIVE_INITCONC', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       Input_Opt%PASSIVE_INITCONC => NULL()
    ENDIF

    IF ( ASSOCIATED( Input_Opt%PASSIVE_DECAYID ) ) THEN
       DEALLOCATE( Input_Opt%PASSIVE_DECAYID )
    ENDIF

    IF ( ASSOCIATED( Input_Opt%AdvectSpc_Name ) ) THEN
       DEALLOCATE( Input_Opt%AdvectSpc_Name, STAT=RC )
       CALL GC_CheckVar( 'Input_Opt%AdvectSpcName', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       Input_Opt%AdvectSpc_Name => NULL()
    ENDIF

    IF ( ASSOCIATED( Input_Opt%SALA_REDGE_um ) ) THEN
       DEALLOCATE( Input_Opt%SALA_REDGE_um, STAT=RC )
       CALL GC_CheckVar( 'Input_Opt%SALA_REDGE_um', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       Input_Opt%SALA_REDGE_um => NULL()
    ENDIF

    IF ( ASSOCIATED( Input_Opt%SALC_REDGE_um ) ) THEN
       DEALLOCATE( Input_Opt%SALC_REDGE_um, STAT=RC )
       CALL GC_CheckVar( 'Input_Opt%SALC_REDGE_um', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       Input_Opt%SALC_REDGE_um => NULL()
    ENDIF

    IF ( ASSOCIATED( Input_Opt%NJDAY ) ) THEN
       DEALLOCATE( Input_Opt%NJDAY, STAT=RC )
       CALL GC_CheckVar( 'Input_Opt%NJDAY', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       Input_Opt%NJDAY => NULL()
    ENDIF

    IF ( ASSOCIATED( Input_Opt%TINDEX ) ) THEN
       DEALLOCATE( Input_Opt%TINDEX, STAT=RC )
       CALL GC_CheckVar( 'Input_Opt%TINDEX', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       Input_Opt%TINDEX => NULL()
    ENDIF

    IF ( ASSOCIATED( Input_Opt%TCOUNT ) ) THEN
       DEALLOCATE( Input_Opt%TCOUNT, STAT=RC )
       CALL GC_CheckVar( 'Input_Opt%TCOUNT', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       Input_Opt%TCOUNT => NULL()
    ENDIF

    IF ( ASSOCIATED( Input_Opt%TMAX ) ) THEN
       DEALLOCATE( Input_Opt%TMAX, STAT=RC )
       CALL GC_CheckVar( 'Input_Opt%TMAX', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       Input_Opt%TMAX => NULL()
    ENDIF

    IF ( ASSOCIATED( Input_Opt%ND51_TRACERS ) ) THEN
       DEALLOCATE( Input_Opt%ND51_TRACERS, STAT=RC )
       CALL GC_CheckVar( 'Input_Opt%ND51_TRACERS', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       Input_Opt%ND51_TRACERS => NULL()
    ENDIF

    IF ( ASSOCIATED( Input_Opt%ND51b_TRACERS ) ) THEN
       DEALLOCATE( Input_Opt%ND51b_TRACERS, STAT=RC )
       CALL GC_CheckVar( 'Input_Opt%ND51b_TRACERS', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       Input_Opt%ND51b_TRACERS => NULL()
    ENDIF

    IF ( ASSOCIATED( Input_Opt%FAM_NAME ) ) THEN
       DEALLOCATE( Input_Opt%FAM_NAME, STAT=RC )
       CALL GC_CheckVar( 'Input_Opt%FAM_NAME', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       Input_Opt%FAM_NAME => NULL()
    ENDIF

    IF ( ASSOCIATED( Input_Opt%FAM_TYPE ) ) THEN
       DEALLOCATE( Input_Opt%FAM_TYPE, STAT=RC )
       CALL GC_CheckVar( 'Input_Opt%FAM_TYPE', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       Input_Opt%FAM_TYPE => NULL()
    ENDIF

    IF ( ASSOCIATED( Input_Opt%LINOZ_TPARM ) ) THEN
       DEALLOCATE( Input_Opt%LINOZ_TPARM, STAT=RC )
       CALL GC_CheckVar( 'Input_Opt%LINOZ_TPARM', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       Input_Opt%LINOZ_TPARM => NULL()
    ENDIF

    IF ( ASSOCIATED( Input_Opt%LSPECRADMENU ) ) THEN
       DEALLOCATE( Input_Opt%LSPECRADMENU, STAT=RC )
       CALL GC_CheckVar( 'Input_Opt%LSPECRADMENU', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       Input_Opt%LSPECRADMENU => NULL()
    ENDIF

    IF ( ASSOCIATED( Input_Opt%LSKYRAD ) ) THEN
       DEALLOCATE( Input_Opt%LSKYRAD, STAT=RC )
       CALL GC_CheckVar( 'Input_Opt%LSKYRAD', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       Input_Opt%LSKYRAD => NULL()
    ENDIF

    IF ( ASSOCIATED( Input_Opt%WVSELECT ) ) THEN
       DEALLOCATE( Input_Opt%WVSELECT, STAT=RC )
       CALL GC_CheckVar( 'Input_Opt%WVSELECT', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       Input_Opt%WVSELECT => NULL()
    ENDIF

    IF ( ASSOCIATED( Input_Opt%STRWVSELECT ) ) THEN
       DEALLOCATE( Input_Opt%STRWVSELECT, STAT=RC )
       CALL GC_CheckVar( 'Input_Opt%STRWVSELECT', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       Input_Opt%STRWVSELECT => NULL()
    ENDIF

    IF ( ASSOCIATED( Input_Opt%ObsPack_SpcName ) ) THEN
       DEALLOCATE( Input_Opt%ObsPack_SpcName, STAT=RC )
       CALL GC_CheckVar( 'Input_Opt%ObsPack_SpcName', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       Input_Opt%ObsPack_SpcName => NULL()
    ENDIF

#ifdef MODEL_GEOS
    !=======================================================================
    ! These fields of Input_Opt are only finalized when
    ! GEOS-Chem is coupled to the online NASA/GEOS ESM
    !=======================================================================
    IF ( ASSOCIATED( Input_Opt%RxnRconst_IDs ) ) THEN
       DEALLOCATE( Input_Opt%RxnRconst_IDs, STAT=RC )
       CALL GC_CheckVar( 'Input_Opt%RxnRconst_IDs', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       Input_Opt%RxnRconst_IDs => NULL()
    ENDIF

    IF ( ASSOCIATED( Input_Opt%RxnRates_IDs ) ) THEN
       DEALLOCATE( Input_Opt%RxnRates_IDs, STAT=RC )
       CALL GC_CheckVar( 'Input_Opt%RxnRates_IDs', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       Input_Opt%RxnRates_IDs => NULL()
    ENDIF

    IF ( ASSOCIATED( Input_Opt%Jval_IDs ) ) THEN
       DEALLOCATE( Input_Opt%Jval_IDs, STAT=RC )
       CALL GC_CheckVar( 'Input_Opt%Jval_Ids', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       Input_Opt%Jval_Ids => NULL()
    ENDIF
#endif

  END SUBROUTINE Cleanup_Input_Opt
!EOC
END MODULE Input_Opt_Mod
