!BOC
#if defined( ESMF_ ) || defined( EXTERNAL_GRID ) || defined( MODEL_ )
      !-----------------------------------------------------------------
      !         %%%%%%% GEOS-Chem HP (with ESMF & MPI) %%%%%%%
      !        %%%% GEOS-Chem Coupled with External Models %%%%
      !
      ! When GEOS-Chem is connected to an external model or in GCHP,
      ! the GEOS-Chem classic main.F should not be built.
      !-----------------------------------------------------------------
#else
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: main.F 
!
! !DESCRIPTION: Program GEOS\_CHEM is the main level driver program for the 
!  GEOS-Chem model of atmospheric chemistry and composition.
!\\
!\\
! !INTERFACE:
!
      PROGRAM GEOS_Chem
!
! !USES:
!
      !-----------------------------------------------------------------
      ! Parameters to define floating-point variables
      !-----------------------------------------------------------------
      USE PRECISION_MOD,   ONLY : fpp => fp  ! Flexible precision
      USE PRECISION_MOD,   ONLY : f4         ! 4-byte floating point
      USE PRECISION_MOD,   ONLY : f8         ! 8-byte floating point

      !-----------------------------------------------------------------
      ! Basic GEOS-Chem modules
      !-----------------------------------------------------------------
      USE CMN_SIZE_MOD          ! Size parameters
      USE DiagList_Mod          ! Derived type for diagnostics list
      USE Diagnostics_Mod       ! Set select netcdf diagnostics
      USE ErrCode_Mod           ! Error codes for success or failure
      USE ERROR_MOD             ! For error checking
      USE FILE_MOD              ! For file I/O
      USE GEOS_TIMERS_MOD       ! For GEOS-Chem timers (optional)
      USE GC_Environment_Mod    ! For allocating derived type objects
      USE GC_GRID_MOD           ! For defining the lons/lats/areas of the grid
      USE Input_Opt_Mod         ! Derived type for Input Options 
      USE INPUT_MOD             ! For reading settings from "input.geos"
      USE MAPPING_MOD           ! For regridding MODIS LAI
      USE OLSON_LANDMAP_MOD     ! Computes IREG, ILAND, IUSE from Olson map
      USE PhysConstants         ! Physical constants
      USE PRESSURE_MOD          ! For computing pressure at grid boxes
      USE Grid_Registry_Mod     ! Registers horizontal/vertical grid metadata
      USE State_Chm_Mod         ! Derived type for Chemistry State object
      USE State_Diag_Mod        ! Derived type for Diagnostics State object
      USE State_Met_Mod         ! Derived type for Meteorology State object
      USE TIME_MOD              ! For computing date & time 
      USE UnitConv_Mod          ! For species conc unit conversions

      !-----------------------------------------------------------------
      ! GEOS-Chem chemistry modules
      !-----------------------------------------------------------------
      USE AEROSOL_MOD,  ONLY : Set_AerMass_Diagnostic
      USE CARBON_MOD            ! For SOA simulation
      USE CHEMISTRY_MOD         ! Driver routines for chemistry
      USE MERCURY_MOD           ! For offline Hg simulation (driver)
      USE OCEAN_MERCURY_MOD     ! For offline Hg simulation (ocean model)
      USE STRAT_CHEM_MOD        ! For linearized stratospheric chemistry
      USE TOMS_MOD              ! For overhead O3 columns (for FAST-J)
      USE UCX_MOD               ! For unified trop-strat chemistry (SDE)
      USE UVALBEDO_MOD          ! For reading UV albedoes (for FAST-J)
      USE SET_GLOBAL_CH4_MOD    ! For setting global CH4 concentrations

      !-----------------------------------------------------------------
      ! GEOS-Chem deposition modules
      !-----------------------------------------------------------------
      USE DEPO_MERCURY_MOD      ! Deposition for offline Hg simulation
      USE DRYDEP_MOD            ! For dry deposition
      USE WETSCAV_MOD           ! For wet deposition (rainout & washout)

      !-----------------------------------------------------------------
      ! GEOS-Chem diagnostics modules
      !-----------------------------------------------------------------
#if defined( BPCH_DIAG )
      USE BENCHMARK_MOD         ! For the 1-month benchmark simulations
      USE DIAG_MOD              ! G-C diagnostic arrays & counters
      USE DIAG41_MOD            ! For ND41  (afternoon PBL       ) diag
      USE DIAG42_MOD            ! For ND42  (SOA products        ) diag
#endif
#if defined( BPCH_DIAG ) || defined( BPCH_TIMESER )
      USE CMN_DIAG_MOD          ! Logical switches for G-C diagnostics
      USE DIAG48_MOD            ! For ND48  (station timeseries  ) diag
      USE DIAG49_MOD            ! For ND49  (inst. timeseries    ) diag
      USE DIAG50_MOD            ! For ND50  (24h avg timeseries  ) diag
      USE DIAG51_MOD            ! For ND51  (satellite timeseries) diag
      USE DIAG51b_MOD           ! For ND51b (satellite timeseries) diag
      USE DIAG63_MOD            ! For ND63  (PARANOX timeseries  ) diag
#endif
      USE DIAG_OH_MOD           ! For ND43  (OH,HO2,etc. prod    ) diag
      USE HISTORY_MOD           ! Updated netCDF diagnostics
      USE OBSPACK_MOD           ! For ObsPack diagnostics
      USE PLANEFLIGHT_MOD       ! For ND40  (plane flight track  ) diag
      USE GOSAT_CH4_MOD         ! For GOSAT observation operator
      USE TCCON_CH4_MOD         ! For TCCON observation operator
      USE HCOI_GC_MAIN_MOD      ! Writes out HEMCO diagnostics (C. Keller)

      !-----------------------------------------------------------------
      ! GEOS-Chem dynamics modules
      !-----------------------------------------------------------------
      USE CONVECTION_MOD        ! For deep cloud convection
      USE LINOZ_MOD             ! For LINOX linear strat chemistry
      USE PBL_MIX_MOD           ! To compute PBL height 
      USE TPCORE_BC_MOD         ! For nested-grid boundary conditions
      USE TRANSPORT_MOD         ! Driver routines for advection
      USE VDIFF_MOD             ! For non-local PBL mixing (J. Lin)

      !-----------------------------------------------------------------
      ! GEOS-Chem emissions modules
      !-----------------------------------------------------------------
      USE EMISSIONS_MOD         ! For interfacing with HEMCO emissions
      USE MIXING_MOD            ! performs tracer mixing
      USE MODIS_LAI_MOD         ! For MODIS leaf area indices (replacement)

      !-----------------------------------------------------------------
      ! GEOS-Chem met field I/O modules
      !-----------------------------------------------------------------
      USE DAO_MOD               ! Met field definitions
      USE FLEXGRID_READ_MOD     ! For reading FLEXGRID data 
#if defined( EXCHANGE )
      USE EXCHANGE_MOD          ! For two-way coupled simulations
#endif

#if defined( RRTMG )
      !-----------------------------------------------------------------
      ! Radiation modules
      !-----------------------------------------------------------------
      USE RRTMG_RAD_TRANSFER_MOD, ONLY : Do_RRTMG_Rad_Transfer
      USE RRTMG_RAD_TRANSFER_MOD, ONLY : Init_RRTMG_Rad_Transfer
      USE RRTMG_RAD_TRANSFER_MOD, ONLY : Set_SpecMask
      USE RRTMG_LW_Init,          ONLY : RRTMG_LW_Ini
      USE RRTMG_SW_Init,          ONLY : RRTMG_SW_Ini
#endif

      IMPLICIT NONE
!
! !REMARKS:
!                                                                             .
!     GGGGGG  EEEEEEE  OOOOO  SSSSSSS       CCCCCC H     H EEEEEEE M     M    
!    G        E       O     O S            C       H     H E       M M M M    
!    G   GGG  EEEEEE  O     O SSSSSSS      C       HHHHHHH EEEEEE  M  M  M    
!    G     G  E       O     O       S      C       H     H E       M     M    
!     GGGGGG  EEEEEEE  OOOOO  SSSSSSS       CCCCCC H     H EEEEEEE M     M    
!                                                                             .
!                                                                             .
!                 (formerly known as the Harvard-GEOS model)
!           for 4 x 5, 2 x 2.5 global grids and hi-res nested grids
!                                                                             .
!       Contact: GEOS-Chem Support Team (geos-chem-support@as.harvard.edu)
!                                                                     
!                                                                             .
!  See the GEOS-Chem Web Site:
!                                                                             .
!     http://acmg.seas.harvard.edu/geos/
!                                                                             .
!  and the GEOS-Chem User's Guide:
!                                                                             .
!     http://acmg.seas.harvard.edu/geos/doc/man/
!                                                                             .
!  and the GEOS-Chem wiki:
!                                                                             .
!     http://wiki.seas.harvard.edu/geos-chem/
!                                                                             .
!  for the most up-to-date GEOS-Chem documentation on the following topics:
!                                                                             .
!     - installation, compilation, and execution
!     - coding practice and style
!     - input files and met field data files
!     - horizontal and vertical resolution
!     - modification history
!
! !REVISION HISTORY: 
!  13 Aug 2010 - R. Yantosca - Added ProTeX headers
!  13 Aug 2010 - R. Yantosca - Add modifications for MERRA (treat like GEOS-5)
!  19 Aug 2010 - R. Yantosca - Now call MERRA met field reader routines
!  02 Feb 2011 - S. Kim      - Call Compute_OD after wet deposition
!  05 Oct 2011 - R. Yantosca - Now get SUNCOS30 array from routine COSSZA
!  07 Oct 2011 - R. Yantosca - Rename SUNCOS30 to SUNCOS_MID, which is the
!                              cos(SZA) at the midpt of the chemistry timestep
!  02 Feb 2012 - R. Yantosca - Added modifications for GEOS-5.7.x met fields
!  06 Feb 2012 - R. Yantosca - Reorganize USE statements for clarity
!  06 Feb 2012 - R. Yantosca - Renamed NN to NNN to avoid name confusion
!  07 Feb 2012 - R. Yantosca - Split off met field I/O into internal routines
!                              READ_INITIAL_MET_FIELDS and READ_MET_FIELDS
!  07 Feb 2012 - M. Payer    - Replace call to COSSZA with GET_COSINE_SZA
!  28 Feb 2012 - R. Yantosca - Removed support for GEOS-3
!  06 Mar 2012 - R. Yantosca - Now call READ_TOMS every month (this was
!                              formerly done within routine "fast_j.F")
!  06 Mar 2012 - R. Yantosca - Add subroutine GET_OVERHEAD_O3_FOR_FASTJ
!                              which calls COMPUTE_OVERHEAD_O3 (in toms_mod.F)
!                              to pre-compute the overhead O3 columsn for
!                              FAST-J photolysis.  This removes code from
!                              "set_prof.F" to facilitate the GI model.
!  19 Mar 2012 - R. Yantosca - Now call routines from olson_landmap_mod.F90
!                              to read the Olson land map data
!  04 Apr 2012 - R. Yantosca - Now call updated LAI routines from new module
!                              modis_lai_mod.F90.  Retire routine RDLAI.
!  05 Apr 2012 - R. Yantosca - Removed reference to LXTRA, it's obsolete
!  11 Apr 2012 - R. Yantosca - Replace lai_mod.F with modis_lai_mod.F90
!  11 Apr 2012 - R. Yantosca - Now call INIT_MODIS_LAI (in modis_lai_mod.F90)
!                              here so that we don't have to call it from 
!                              megan_mod.F and mercury_mod.F separately.
!  17 Apr 2012 - R. Yantosca - Need to set the mapping variable to NULL()
!  10 Jun 2012 - L. Murray   - Remove references to UPBDFLX_MOD.F
!  31 Jul 2012 - R. Yantosca - Now pass am_I_Root variable to lower-level
!                              routines in order to allow PRINT and WRITE
!                              statements to execute on the root CPU.  This
!                              is needed for compatibility w/ the GEOS-5 GCM.
!  13 Aug 2012 - R. Yantosca - Now call FILL_CHEM_STATE_IDs to populate
!                              the CHEM_STATE object ID and name fields
!  18 Oct 2012 - R. Yantosca - Rename LOCAL_MET object to State_Met
!  18 Oct 2012 - R. Yantosca - Rename CHEM_STATE object to State_Chm
!  18 Oct 2012 - R. Yantosca - Now pass am_I_Root, RC arguments to routines
!                              ALLOCATE_ALL, INIT_ALL when using -DDEVEL
!  19 Oct 2012 - R. Yantosca - Now reference gigc_state_chm_mod.F90
!  19 Oct 2012 - R. Yantosca - Now reference gigc_state_met_mod.F90
!  25 Oct 2012 - R. Yantosca - Define logical doDebugPrt for ND70 output
!  25 Oct 2012 - R. Yantosca - Add descriptive comments for DEVEL #ifdefs
!  25 Oct 2012 - R. Yantosca - Now reference gigc_errcode_mod.F90
!  01 Nov 2012 - R. Yantosca - Now read soil NOx restart file
!  01 Nov 2012 - R. Yantosca - Now reference gigc_input_opt_mod.F90
!  08 Nov 2012 - R. Yantosca - Now pass Input_Opt as an arg to DO_CHEMISTRY
!  01 Nov 2012 - R. Yantosca - Now read soil NOx restart file
!  14 Nov 2012 - R. Yantosca - Add am_I_Root, Input_Opt, RC as arguments
!                              to various subroutines
!  15 Nov 2012 - M. Payer    - Replaced all met field arrays with State_Met
!                              derived type object
!  15 Nov 2012 - R. Yantosca - Bring Input_Opt out of the DEVEL tags
!  26 Feb 2013 - R. Yantosca - Add placeholder tag for Input_Opt%MAX_DEP
!  05 Mar 2013 - R. Yantosca - Now pass am_I_Root, Input_Opt, RC to routine
!                              DO_PBL_MIX_2 (for non-local PBL mixing)
!  15 Mar 2013 - R. Yantosca - Now set Input_Opt%LINOZ_N* fields here
!  26 Mar 2013 - S.D. Eastham- Added initialization of rare tracers
!  29 Mar 2013 - R. Yantosca - Bring code out of DEVEL blocks
!  30 May 2013 - R. Yantosca - Now pass Input_Opt object to STDRUN routine
!  03 Jun 2013 - R. Yanosca - Use routines from updated mercury_mod.F
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!  23 Oct 2013 - R. Yantosca - Now pass am_I_root, Input_Opt, RC to INIT_DAO
!  13 Dec 2013 - M. Sulprizio- Now set USE_O3_FROM_MET logical flag during
!                              initialization stage
!  11 Apr 2014 - R. Yantosca - Now remove call to INIT_GLOBAL_CH4; this is
!                              now done from routine GIGC_Init_Extra
!  19 May 2014 - C. Keller   - Added INIT_CHEMISTRY 
!  19 May 2014 - C. Keller   - Added HEMCO
!  23 Jun 2014 - R. Yantosca - Removed references to unix_cmds_mod.F
!  23 Jun 2014 - R. Yantosca - Removed references to directory_mod.F
!  23 Jun 2014 - R. Yantosca - Removed references to logical_mod.F
!  15 Jul 2014 - R. Yantosca - Now reference grid_mod.F90, regrid_a2a_mod.F90
!  15 Jul 2014 - R. Yantosca - Now call Init_Map_A2A to store shadow variables
!                              within regrid_a2a_mod.F90.  This helps to 
!                              break dependencies for the HEMCO implementation.
!  25 Jul 2014 - R. Yantosca - Remove reference to commsoil_mod.F90
!  22 Aug 2014 - R. Yantosca - Now save areas [m2] in State_Met%AREA_M2
!  15 Dec 2014 - M. Yannetti - Added PRECISION_MOD
!  06 Jan 2015 - R. Yantosca - Added two-way coupled simulation options
!  17 Feb 2015 - E. Lundgren - Remove STT and TCVV pointers
!  25 Feb 2015 - E. Lundgren - Remove MAKE-RH call since now in AIRQNT
!  16 Mar 2015 - E. Lundgren - Change tracer main units from kg to kg/kg
!  24 Mar 2015 - E. Lundgren - Now pass Input_Opt to Check_STT
!  31 Mar 2015 - E. Lundgren - Move post-transport AIRQNT call to transport_mod 
!  16 Apr 2015 - R. Yantosca - Remove call to INIT_DAO; it's obsolete
!  12 Aug 2015 - R. Yantosca - Add support for MERRA2 meteorology
!  03 Sep 2015 - R. Yantosca - Now call SETUP_WETSCAV instead of INIT_WETSCAV
!  25 Jan 2016 - R. Yantosca - Call LINE_BUFFER to force PGI compiler to flush
!                               STDOUT (unit=6) output to disk during a run
!  03 Feb 2016 - E. Lundgren - Use routine MAKE_RESTART_FILES for all GC rsts
!  29 Apr 2016 - R. Yantosca - Don't initialize pointers in declaration stmts
!  02 May 2016 - R. Yantosca - Now pass Input_Opt to cleanup_gigc_state_chm
!  18 May 2016 - M. Sulprizio- Remove call to INIT_COMODE; it's obsolete
!  06 Jun 2016 - M. Sulprizio- Remove call to FILL_CHEM_STATE_IDs; this routine
!                              was made obsolete by the species database
!  22 Jun 2016 - R. Yantosca - Add error checks to prevent calling UCX routines
!                              when we are running specialty simulations
!  12 Jul 2016 - E. Lundgren - Remove binary punch restart file option
!  29 Nov 2016 - R. Yantosca - grid_mod.F90 is now gc_grid_mod.F90
!  16 Jun 2017 - M. Sulprizio- Add call for GOSAT CH4 observation operator
!  17 Aug 2017 - M. Sulprizio- Add call for TCCON CH4 observation operator
!  23 Aug 2017 - R. Yantosca - Now add at timer for the History component
!  24 Aug 2017 - M. Sulprizio- Remove support for GCAP, GEOS-4, GEOS-5 and MERRA
!  30 Aug 2017 - R. Yantosca - Move the HISTORY UPDATE block before call to
!                               SET_CURRENT_TIME, for clarity
!  13 Sep 2017 - M. Sulprizio- Remove Input_Opt%USE_OLSON_2001. Olson 2001 is
!                              now the default.
!  14 Sep 2017 - M. Sulprizio- Remove setting of USE_O3_FROM_MET. This is now
!                              set in the chemistry menu of input.geos.
!  22 Jan 2018 - R. Yantosca - Attach error traps after every routine w/ RC
!  26 Oct 2018 - M. Sulprizio- Remove routines READ_INITIAL_MET_FIELDS and
!                              READ_MET_FIELDS. Met fields are now obtained from
!                              HEMCO in hcoi_gc_main_mod.F90.
!  08 Nov 2018 - R. Yantosca - Update the initializations for RRTMG
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      !----------------------------
      ! Scalars
      !-----------------------------

      ! Logicals
      LOGICAL                  :: ITS_A_FULLCHEM_SIM
      LOGICAL                  :: ITS_A_MERCURY_SIM 
      LOGICAL                  :: ITS_A_TAGCO_SIM
      LOGICAL                  :: ITS_AN_AEROSOL_SIM
      LOGICAL                  :: ITS_A_CH3I_SIM
      LOGICAL                  :: DO_DIAG_WRITE
      LOGICAL                  :: LCHEM
      LOGICAL                  :: LCONV
      LOGICAL                  :: LDRYD
      LOGICAL                  :: LDYNOCEAN
      LOGICAL                  :: LEMIS
      LOGICAL                  :: LGTMM
      LOGICAL                  :: LLINOZ
      LOGICAL                  :: LNLPBL
      LOGICAL                  :: LPRT
      LOGICAL                  :: LSTDRUN
      LOGICAL                  :: LSCHEM
      LOGICAL                  :: LSETH2O
      LOGICAL                  :: LTRAN
      LOGICAL                  :: LTURB
      LOGICAL                  :: LUCX
      LOGICAL                  :: LWETD
      LOGICAL                  :: LWINDO
      LOGICAL                  :: wasModisRead
      LOGICAL                  :: prtDebug
      LOGICAL                  :: TimeForEmis

      ! Integers
      INTEGER                  :: I,             IOS,         J
      INTEGER                  :: K,             L,           N
      INTEGER                  :: JDAY,          N_DYN
      INTEGER                  :: NNN,           N_DYN_STEPS, NSECb 
      INTEGER                  :: N_STEP,        YEAR,        MONTH
      INTEGER                  :: DAY,           DAY_OF_YEAR, SEASON
      INTEGER                  :: NYMD,          NYMDb,       NHMS
      INTEGER                  :: ELAPSED_SEC,   NHMSb,       RC
      INTEGER                  :: ELAPSED_TODAY, HOUR,        MINUTE
      INTEGER                  :: laiYYYY,       laiMM
      INTEGER                  :: id_H2O,        id_CH4,      SECOND

      ! Reals
      REAL(f8)                 :: TAU,           TAUb      
      REAL(fpp)                :: COST_FUNC ! Add here temporarily (mps,6/16/17)

      ! Strings
      CHARACTER(LEN=255)       :: ThisLoc,       ZTYPE
      CHARACTER(LEN=255)       :: historyConfigFile
      CHARACTER(LEN=512)       :: ErrMsg
      CHARACTER(LEN=512)       :: Instr

#if defined( RRTMG )
      !-----------------------------
      ! Scalars specific to RRTMG
      !-----------------------------
      INTEGER                  :: iCld
      INTEGER                  :: iSeed
      INTEGER                  :: iSpecMenu
      INTEGER                  :: iNcDiag
      INTEGER                  :: RADSPEC
      LOGICAL, SAVE            :: FIRST_RT = .TRUE.
#endif

      !-----------------------------
      ! Derived type objects 
      !-----------------------------
      TYPE(MapWeight), POINTER :: mapping(:,:)
      TYPE(OptInput)           :: Input_Opt   ! Input Options object
      TYPE(ChmState)           :: State_Chm   ! Chemistry State object
      TYPE(DgnState)           :: State_Diag  ! Diagnostics State object
      TYPE(MetState)           :: State_Met   ! Meteorology State object
      TYPE(DgnList )           :: Diag_List   ! Diagnostics list object
!
! !DEFINED PARAMETERS:
!
      ! When connecting G-C to an external GCM, we need to only write 
      ! to stdout if we are on the root CPU.  Otherwise this will slow
      ! down the code.  This is why we introduced the am_I_Root logical
      ! variable.
      !
      ! However, if we are using the "traditional" G-C, then we don't
      ! need to restrict I/O to the root CPU.  (This is because each
      ! GEOS-Chem simulation starts out on a single CPU, with other
      ! CPUs joining only within parallel DO loops).  Therefore, we
      ! can just set am_I_Root = .true. here and then have it propagate
      ! down to all of the lower-level routines.  The main.F routine
      ! is not called when connecting G-C to an external GCM.
      ! (mlong, bmy, 7/30/12)
      LOGICAL, PARAMETER       :: am_I_Root = .TRUE. 

      !=================================================================
      ! GEOS-CHEM starts here!                                            
      !=================================================================

      ! Tell PGI compiler to line-buffer STDOUT, so that we can pipe the 
      ! screen output to a log file and see the results while GC is running.
      CALL LINE_BUFFER( 6 )

#if defined( TOMAS )
      !(sfarina, 6/19/2013) It may seem strange, but this welcome message
      !                     fixes an issue where geoschem crashes with a
      !                     sigsegv immediately after starting.
      !                     This happens on ace-net's glooscap cluster with
      !                     ifort (IFORT) 11.1 20101201
      !                     this issue does not appear when running inside
      !                     a debugger, and is probably related to
      !                     some initialization garbage in memory
      !                     when using -O2 optimization
      !(bmy, 1/27/2014)   - Need to "CALL FLUSH(6).  FLUSH needs
      !                     an argument.  Unit 6 is Unix stdout.
      PRINT*, '%%%%% USING TOMAS MICROPHYSICS PACKAGE %%%%%'
      CALL FLUSH(6)
#endif

      ! Assume a successful return until otherwise
      RC                      = GC_SUCCESS

      ! For error trapping
      ErrMsg                  = ''
      ThisLoc                 = 'GEOS-Chem (in GeosCore/main.F)'

      ! Initialize pointers
      mapping => NULL()

      ! Display current grid resolution and data set type
      CALL Display_Grid_and_Model()

      !=================================================================
      !            ***** I N I T I A L I Z A T I O N *****
      !=================================================================

#if defined( USE_TIMERS )
      !-----------------------------------------------------------------
      ! Initialize GEOS-Chem timers (if compiled with TIMERS=1)
      !-----------------------------------------------------------------

      ! Call timer initilization
      CALL GEOS_Timer_Setup( 1 )

      ! Add timers for various operations
      CALL GEOS_Timer_Add( "GEOS-Chem",                            RC )
      CALL GEOS_Timer_Add( "Initialization",                       RC )
      CALL GEOS_Timer_Add( "Timesteps",                            RC )
      CALL GEOS_Timer_Add( "HEMCO",                                RC )
      CALL GEOS_Timer_Add( "All chemistry",                        RC )
      CALL GEOS_Timer_Add( "=> Gas-phase chem",                    RC )
      CALL GEOS_Timer_Add( "=> FAST-JX photolysis",                RC )
      CALL GEOS_Timer_Add( "=> All aerosol chem",                  RC )
      CALL GEOS_Timer_Add( "=> Strat chem",                        RC )
      CALL GEOS_Timer_Add( "=> Unit conversions",                  RC )
      CALL GEOS_Timer_Add( "Transport",                            RC )
      CALL GEOS_Timer_Add( "Convection",                           RC )
      CALL GEOS_Timer_Add( "Boundary layer mixing",                RC )
      CALL GEOS_Timer_Add( "Dry deposition",                       RC )
      CALL GEOS_Timer_Add( "Wet deposition",                       RC )
#if defined( RRTMG )
      CALL GEOS_Timer_Add( "RRTMG",                                RC )
#endif
      CALL GEOS_Timer_Add( "All diagnostics",                      RC )
      CALL GEOS_Timer_Add( "=> HEMCO diagnostics",                 RC )
      CALL GEOS_Timer_Add( "=> Binary punch diagnostics",          RC )
      CALL GEOS_Timer_Add( "=> ObsPack diagnostics",               RC )
      CALL GEOS_Timer_Add( "=> History (netCDF diags)",            RC )
      CALL GEOS_Timer_Add( "Input",                                RC )
      CALL GEOS_Timer_Add( "Output",                               RC )
      CALL GEOS_Timer_Add( "Finalization",                         RC )

      ! Start running the main and initialization timer            
      CALL GEOS_Timer_Start( "GEOS-Chem",                          RC )
      CALL GEOS_Timer_Start( "Initialization",                     RC )
#endif

      !-----------------------------------------------------------------
      ! Read the user-defined options for the simulation, etc.
      !-----------------------------------------------------------------

      ! Initialize fields of the Input Options object
      CALL Set_Input_Opt( am_I_Root, Input_Opt, RC )
      IF ( RC /= GC_SUCCESS ) THEN
         ErrMsg = 'Error encountered within call to "Set_Input_Opt"!'
         CALL Error_Stop( ErrMsg, ThisLoc )
      ENDIF

      ! Read GEOS-Chem input file at very beginning of simulation
      CALL Read_Input_File( am_I_Root, Input_Opt, RC ) 
      IF ( RC /= GC_SUCCESS ) THEN
         ErrMsg = 'Error encountered in "Read_Input_File"!'
         CALL Error_Stop( ErrMsg, ThisLoc )
      ENDIF

      ! Call the routine GC_Allocate_All (located in module file
      ! GeosCore/gc_environment_mod.F90) to allocate all lat/lon
      ! allocatable arrays used by GEOS-Chem.  The IIPAR and JJPAR
      ! dimensions are not declared as PARAMETERs, but are module
      ! variables that are initialized in this call.
      CALL GC_Allocate_All( am_I_root, Input_Opt, RC )
      IF ( RC /= GC_SUCCESS ) THEN
         ErrMsg = 'Error encountered in "GC_Allocate_All"!'
         CALL Error_Stop( ErrMsg, ThisLoc )
      ENDIF

      ! Store shadow copies of am_I_Root, Input_Opt in error_mod.F
      CALL Init_Error( am_I_Root, Input_Opt, RC )
      IF ( RC /= GC_SUCCESS ) THEN
         ErrMsg = 'Error encountered in "Init_Error"!'
         CALL Error_Stop( ErrMsg, ThisLoc )
      ENDIF

      ! Copy values from Input_Opt.  These replace the variables
      ! from logical_mod.F and tracer_mod.F. (bmy, 3/29/13)
      ITS_A_FULLCHEM_SIM  =  Input_Opt%ITS_A_FULLCHEM_SIM
      ITS_A_MERCURY_SIM   =  Input_Opt%ITS_A_MERCURY_SIM
      ITS_A_TAGCO_SIM     =  Input_Opt%ITS_A_TAGCO_SIM 
      ITS_AN_AEROSOL_SIM  =  Input_Opt%ITS_AN_AEROSOL_SIM
      ITS_A_CH3I_SIM      =  Input_Opt%ITS_A_CH3I_SIM
      DO_DIAG_WRITE       =  Input_Opt%DO_DIAG_WRITE
      LCHEM               =  Input_Opt%LCHEM
      LCONV               =  Input_Opt%LCONV
      LDRYD               =  Input_Opt%LDRYD
      LDYNOCEAN           =  Input_Opt%LDYNOCEAN
      LEMIS               =  Input_Opt%LEMIS
      LGTMM               =  Input_Opt%LGTMM
      LLINOZ              =  Input_Opt%LLINOZ
      LNLPBL              =  Input_Opt%LNLPBL
      LPRT                =  Input_Opt%LPRT
      LSCHEM              =  Input_Opt%LSCHEM
      LSETH2O             =  Input_Opt%LSETH2O
      LSTDRUN             =  Input_Opt%LSTDRUN
      LTRAN               =  Input_Opt%LTRAN
      LTURB               =  Input_Opt%LTURB
      LUCX                =  Input_Opt%LUCX
      LWETD               =  Input_Opt%LWETD
      LWINDO              =  Input_Opt%LWINDO

      ! Set a flag to denote if we should print ND70 debug output
      prtDebug            = ( LPRT .and. am_I_Root )
      
      ! Debug output
      IF ( prtDebug ) CALL Debug_Msg( '### MAIN: a READ_INPUT_FILE' )

      !-----------------------------------------------------------------
      ! %%%% REPLICATING GCHP FUNCTIONALITY IN EXISTING GEOS-CHEM %%%%
      !
      ! Initialize the diagnostic list object which contains the
      ! unique entires in the history config file. Note that this is
      ! done in GCHP Set_Services and therefore must be done prior to
      ! initialization of the state objects. Also note that the diag_list
      ! obj may be stored in the HistoryConfig object in GCHP and we may
      ! want to replicate that behavior in GCC in the future. 
      ! (ewl, 9/26/17)
      !-----------------------------------------------------------------
#if defined( USE_TIMERS )
      CALL GEOS_Timer_Start( "All diagnostics",           RC )
      CALL GEOS_Timer_Start( "=> History (netCDF diags)", RC )
#endif      

      historyConfigFile = 'HISTORY.rc' ! InputOpt not yet initialized
      CALL Init_DiagList( am_I_Root, historyConfigFile, Diag_List, RC )
      IF ( RC /= GC_SUCCESS ) THEN
         ErrMsg = 'Error encountered in "Init_DiagList"!'
         CALL Error_Stop( ErrMsg, ThisLoc )
      ENDIF

      !### Print diagnostic list if needed for debugging 
      IF ( prtDebug ) CALL Print_DiagList( am_I_Root, Diag_List, RC )

#if defined( USE_TIMERS )
      CALL GEOS_Timer_End( "All diagnostics",           RC )
      CALL GEOS_Timer_End( "=> History (netCDF diags)", RC )
#endif

      !-----------------------------------------------------------------
      ! %%%% REPLICATING GCHP FUNCTIONALITY IN EXISTING GEOS-CHEM %%%%
      !
      ! To replicate the functionality of the ESMF interface, we must
      ! initialize the Meteorology State (i.e. State_Met) and the
      ! Chemistry State (i.e. State_Chm) objects.  These objects hold
      ! several individual data fields that need to be passed as 
      ! inputs to the chemistry routines.
      !
      ! The Meteorology State has replaced all of the individual
      ! met field arrays contained in module dao_mod.F. Likewise,
      ! the Chemistry State has replaced the STT tracer array
      ! and CSPEC chemical species array.
      !
      ! The Chemistry and Meteorology State objects facilitate using
      ! GEOS-Chem directly from the ESMF interface.  This is the main
      ! reason we are migrating towards used of these objects instead
      ! of the existing ALLOCATABLE module arrays. (bmy, 10/25/12)
      !-----------------------------------------------------------------

      ! Initialize State_Met, State_Chm, and State_Diag objects
      CALL GC_Init_StateObj( am_I_Root, Diag_List,  Input_Opt, 
     &                       State_Chm, State_Diag, State_Met, RC )
      IF ( RC /= GC_SUCCESS ) THEN
         ErrMsg = 'Error encountered in "GC_Init_All!"!'
         CALL Error_Stop( ErrMsg, ThisLoc )
      ENDIF

      ! Initialize GEOS-Chem horizontal grid structure
      CALL GC_Init_Grid( am_I_Root, Input_Opt, RC )
      IF ( RC /= GC_SUCCESS ) THEN
         ErrMsg = 'Error in "GC_Init_Grid"!'
         CALL Error_Stop( ErrMsg, ThisLoc )
      ENDIF

      ! Also allocate arrays in GEOS-Chem modules that have not yet
      ! been initialized (i.e. SEASALT_MOD, CARBON_MOD,  DUST_MOD, 
      ! SULFATE_MOD).  This removes the init calls from the run stage,
      ! which is required when connecting to the GEOS-5 GCM w/ ESMF.
      ! (bmy, 3/4/13)
      CALL GC_Init_Extra( am_I_Root, Diag_List,  Input_Opt, 
     &                    State_Chm, State_Diag, RC         )
      IF ( RC /= GC_SUCCESS ) THEN
         ErrMsg = 'Error encountered in "GC_Init_Extra"!'
         CALL Error_Stop( ErrMsg, ThisLoc )
      ENDIF

      ! Initialize the regridding module by storing shadow copies
      CALL GC_Init_Regridding( am_I_Root, Input_Opt, RC )
      IF ( RC /= GC_SUCCESS ) THEN
         ErrMsg = 'Error encountered in "Initialize_Regridding"!'
         CALL Error_Stop( ErrMsg, ThisLoc )
      ENDIF

      ! Copy surface area info from AREA_M2 array of grid_mod.F90
      ! to the State_Met%AREA_M2 field. (bmy, 8/22/14)
      ! NOTE: Most routines now use either mixing ratio or kg/m2
      ! and are independent of area
      State_Met%AREA_M2 =  AREA_M2

      ! Define advected species ID flags for use below
      id_H2O   = Ind_('H2O'     ,'A')
      id_CH4   = Ind_('CH4'     ,'A')
      !
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      ! OBSPACK Diagnostics: Get information from the species
      ! database for all requested ObsPack output species
      !-----------------------------------------------------------------
      IF ( Input_Opt%Do_ObsPack ) THEN
         CALL ObsPack_SpeciesMap_Init( am_I_Root, Input_Opt, 
     &                                 State_Chm, State_Diag, RC )
         IF ( RC /= GC_SUCCESS ) THEN
            ErrMsg = 'Error encountered in "ObsPack_SpeciesMap_Init"!'
            CALL Error_Stop( ErrMsg, ThisLoc )
         ENDIF
      ENDIF
         
#if defined( RRTMG )
      !-----------------------------------------------------------------
      ! Initializations for the RRTMG radiative transfer model
      !-----------------------------------------------------------------

      ! Initialize module variables
      CALL Init_RRTMG_Rad_Transfer( am_I_Root,  Input_Opt,
     &                              State_Diag, RC         )

      ! Trap potential errors
      IF ( RC /= GC_SUCCESS ) THEN
         ErrMsg = 'Error encountered in "Init_RRTMG_Rad_Transfer"!'
         CALL Error_Stop( ErrMsg, ThisLoc )
      ENDIF

      ! Initialize RRTMG code in the GeosRad folder
      CALL Rrtmg_Lw_Ini()
      CALL Rrtmg_Sw_Ini()

      ! Settings
      iCld  = 0
      iSeed = 10
#endif

#if defined( BPCH_DIAG )
      !-----------------------------------------------------------------
      ! Initialize bpch diagnostic arrays and counters
      !-----------------------------------------------------------------
      CALL Initialize( am_I_Root, Input_Opt, 2, RC )
      CALL Initialize( am_I_Root, Input_Opt, 3, RC )
      IF ( prtDebug ) CALL Debug_Msg( '### MAIN: a INITIALIZE' )
#endif

      !-----------------------------------------------------------------
      ! Initialize the hybrid pressure module.  Define Ap and Bp.
      !-----------------------------------------------------------------
      CALL Init_Pressure( am_I_Root )
      IF ( prtDebug ) CALL Debug_Msg( '### MAIN: a INIT_PRESSURE' )

      !-----------------------------------------------------------------
      ! Register the horizontal and vertical grid information so that 
      ! the History component can use it for netCDF metadata
      !-----------------------------------------------------------------
      CALL Init_Grid_Registry( am_I_Root, RC )
      IF ( RC /= GC_SUCCESS ) THEN
         ErrMsg = 'Error encountered in "Init_Grid_Registry"!'
         CALL Error_Stop( ErrMsg, ThisLoc )
      ENDIF

      !-----------------------------------------------------------------
      ! Added to read input file for linoz strat chem
      !-----------------------------------------------------------------
      IF ( LLINOZ ) THEN
         CALL Linoz_Read( am_I_Root, Input_Opt, RC )

          ! Trap potential errors
         IF ( RC /= GC_SUCCESS ) THEN
            ErrMsg = 'Error encountered in "Linoz_Read"!'
            CALL Error_Stop( ErrMsg, ThisLoc )
         ENDIF
      ENDIF

      ! Define time variables for use below
      NHMS  = GET_NHMS()
      NHMSb = GET_NHMSb()
      NYMD  = GET_NYMD()
      NYMDb = GET_NYMDb()
      TAU   = GET_TAU()
      TAUb  = GET_TAUb()

      !-----------------------------------------------------------------
      !    ***** H I S T O R Y   I N I T I A L I Z A T I O N *****
      !-----------------------------------------------------------------
#if defined( USE_TIMERS )
      CALL GEOS_Timer_Start( "All diagnostics",           RC )
      CALL GEOS_Timer_Start( "=> History (netCDF diags)", RC )
#endif   

      ! For now, just hardwire the input file for the History component
      Input_Opt%HistoryInputFile = './HISTORY.rc'

      ! Initialize the GEOS-Chem history component
      CALL History_Init( am_I_root, Input_Opt,  State_Met,
     &                   State_Chm, State_Diag, RC         )

      ! Trap error
      IF ( RC /= GC_SUCCESS ) THEN
         ErrMsg = 'Error encountered in "History_Init"!'
         CALL Error_Stop( ErrMsg, ThisLoc )
      ENDIF

#if defined( USE_TIMERS )
      CALL GEOS_Timer_End( "All diagnostics",           RC )
      CALL GEOS_Timer_End( "=> History (netCDF diags)", RC )
#endif   

      !-----------------------------------------------------------------
      !        ***** I N I T I A L I Z A T I O N  continued *****
      !-----------------------------------------------------------------

      ! To enable FlexGrid, need to initialize HEMCO and run phase 1
      ! before reading initial metfields.
      ! (Jiawei Zhuang 2017/6)

      ! Initialize HEMCO. This reads the HEMCO configuration file
      ! and creates entries for all data files needed for emission
      ! calculation. Also sets some logical switches in Input_Opt
      ! (e.g. LSOILNOX).
      ! Note: always call HEMCO, even if LEMIS is set to FALSE. This
      ! is to make sure that HEMCO can still be used to read
      ! non-emission data such as stratospheric Bry fields. If LEMIS
      ! is set to FALSE, the emissions driver routines will make sure
      ! that HEMCO does not calculate any emissions (ckeller, 1/12/15).
      ! 
      ! This call will also initialize the three built-in HEMCO 
      ! diagnostics (default, manual, restart).
#if defined( USE_TIMERS )
      CALL GEOS_Timer_Start( "HEMCO", RC )
#endif

      CALL Emissions_Init( am_I_Root, Input_Opt,
     &                     State_Met, State_Chm, RC )

      ! Trap potential errors
      IF ( RC /= GC_SUCCESS ) THEN
         ErrMsg = 'Error encountered in "Emissions_Init"!'
         CALL Error_Stop( ErrMsg, ThisLoc )
      ENDIF
       
      ! Run HEMCO phase 0 as simplfied phase 1 to get initial met fields
      ! and restart file fields
      CALL Emissions_Run( am_I_Root, Input_Opt,   State_Met,
     &                    State_Chm, State_Diag,  TimeForEmis,
     &                    0,         RC )

      ! Trap potential errors
      IF ( RC /= GC_SUCCESS ) THEN
         ErrMsg = 'Error encountered in "Emissions_Run", Phase 0'
         Instr  =
     &     '(1) This error can indicate a missing GEOS-Chem restart '//
     &     'file (GEOSChem.Restart.YYYYMMDD_hhmmz.nc4), or that the '//
     &     'restart file date does not match the start date of the ' //
     &     'simulation.  Please double-check your run directory to ' //
     &     'make sure that there is a valid restart file for the '   //
     &     'simulation you are using.  (2) This error can also '     //
     &     'indicate that a met field file is missing.  Please '     //
     &     'make sure that you have all relevant met field data '    //
     &     'on disk for the start date of your simulation.'
         CALL Error_Stop( ErrMsg, ThisLoc, Instr )
      ENDIF

#if defined( USE_TIMERS )
      CALL GEOS_Timer_End ( "HEMCO", RC )
#endif

      ! Initialize the derived type object containing
      ! mapping information for the MODIS LAI routines
      CALL Init_Mapping( am_I_Root, Input_Opt,
     &                   1440, 720, IIPAR, JJPAR, mapping, RC )

      ! Trap potential errors
      IF ( RC /= GC_SUCCESS ) THEN
         ErrMsg = 'Error encountered in "Init_Mapping"!'
         CALL Error_Stop( ErrMsg, ThisLoc )
      ENDIF

      ! Compute the Olson land types that occur in each grid box
      ! (i.e. this is a replacement for rdland.F and vegtype.global)
      CALL Init_Olson_Landmap( am_I_Root, Input_Opt, RC )

      ! Trap potential errors
      IF ( RC /= GC_SUCCESS ) THEN
         ErrMsg = 'Error encountered in "Init_Olson_Landmap"!'
         CALL Error_Stop( ErrMsg, ThisLoc )
      ENDIF

      ! Continue computing the Olson landmap
      CALL Compute_Olson_Landmap( am_I_Root, mapping,   State_Met )
      CALL Cleanup_Olson_Landmap( am_I_Root                       )

      ! Initialize PBL quantities but do not do mixing
      ! Add option for non-local PBL (Lin, 03/31/09)
      CALL Init_Mixing( am_I_Root, Input_Opt,  State_Met,
     &                  State_Chm, State_Diag, RC         )

      ! Trap potential errors
      IF ( RC /= GC_SUCCESS ) THEN
         ErrMsg = 'Error encountered in Init_Mixing!'
         CALL Error_Stop( ErrMsg, ThisLoc )
      ENDIF

      ! Initialize chemistry
      ! Moved here because some of the variables are used for non-local 
      ! PBL mixing BEFORE the first call of the chemistry routines
      ! (ckeller, 05/19/14).
      IF ( ITS_A_FULLCHEM_SIM .OR. ITS_AN_AEROSOL_SIM ) THEN
         CALL Init_Chemistry( am_I_Root, Input_Opt,
     &                        State_Chm, State_Diag, RC )

         ! Trap potential errors
         IF ( RC /= GC_SUCCESS ) THEN
            ErrMsg = 'Error encountered in "Init_Chemistry"!'
            CALL Error_Stop( ErrMsg, ThisLoc )
         ENDIF
      ENDIF

      !=================================================================
      !       *****  I N I T I A L   C O N D I T I O N S *****
      !=================================================================

      ! Initialize the UCX module
      IF ( LUCX ) THEN
         CALL INIT_UCX( am_I_Root, Input_Opt, State_Chm, State_Diag )
         IF ( prtDebug ) CALL DEBUG_MSG( '### MAIN: a INIT_UCX' )
      ENDIF

      ! Capture initial state of atmosphere for STE flux calc (ltm, 06/10/12)
      ! NOTE: Species concentrations enter the subroutine in [kg/kg dry air]
      ! and are converted locally to [kg] for chemistry (ewl, 9/18/15)
      IF ( LSCHEM ) THEN
         CALL Init_Strat_Chem( am_I_Root, Input_Opt, State_Chm, 
     &                         State_Met, RC                    )

         ! Trap potential errors
         IF ( RC /= GC_SUCCESS ) THEN
            ErrMsg = 'Error encountered in "Init_Strat_Chem"!'
            CALL Error_Stop( ErrMsg, ThisLoc )
         ENDIF
      ENDIF

#if defined( BPCH_DIAG )
      !-----------------------------------------------------------------
      ! Write initial conditions for benchmark simulations
      !-----------------------------------------------------------------
#if defined( USE_TIMERS )
      CALL GEOS_Timer_Start( "All diagnostics",             RC )
      CALL GEOS_Timer_Start( "=> Binary punch diagnostics", RC )
      CALL GEOS_Timer_Start( "Output",                      RC )
#endif 

      ! Save initial tracer mass [kg] to disk for model benchmarking
      ! for Rn-Pb-Be and full-chem benchmark simulations
      ! NOTE: Species concentrations enter the subroutine in [kg/kg dry air]
      ! and are converted locally to [kg] (mps, 11/29/16)
      IF ( LSTDRUN ) THEN 
         CALL StdRun( am_I_Root,    Input_Opt,  State_Met,
     &                State_Chm,    State_Diag, RC,
     &                LBEGIN=.TRUE.                        )

         ! Trap potential errors
         IF ( RC /= GC_SUCCESS ) THEN
            ErrMsg = 'Error encountered in "StdRun"!'
            CALL Error_Stop( ErrMsg, ThisLoc )
         ENDIF
      ENDIF

#if defined( USE_TIMERS )
      CALL GEOS_Timer_End( "All diagnostics",                  RC )
      CALL GEOS_Timer_End( "=> Binary punch diagnostics",      RC )
      CALL GEOS_Timer_End( "Output",                      RC )
#endif 

#endif

!-----------------------------------------------------------------------------
! TWO-WAY NESTING OPTION
! This is only invoked when compiling GEOS-Chem with COUPLE=y
! 
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! %%%%% NOTE: THIS OPTION WAS PROBABLY BROKEN WITH THE ADDITION OF HEMCO %%%%%
! %%%%% AND FLEXCHEM INTO GEOS-CHEM V10-01 AND v11-01.  BUYER BEWARE.    %%%%%
! %%%%%  --  Bob Yantosca (22 Jan 2018)                                  %%%%%
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
#if defined( EXCHANGE )

      ! Initialize the two-way nesting.  This will only get done if you
      ! compile GEOS-Chem with EXCHANGE=yes. (yanyy, 03/28/14)
      CALL INIT_EXCHANGE( am_I_Root )

#if defined( GRID4x5 ) || defined( GRID2x25 )

      ! Initialize exchange of global boundary conditions
      CALL EXCHANGE_GLOBAL_INIT()

#elif defined( NESTED_CH ) || defined( NESTED_NA ) || defined( NESTED_EU )
      ! Initialize exchange of nested boundary conditions
      CALL EXCHANGE_NESTED_INIT()

#endif

#endif
!-----------------------------------------------------------------------------

#if defined( USE_TIMERS )
      CALL GEOS_Timer_End( "Initialization", RC )
#endif

      !=================================================================
      !        ***** O U T E R   T I M E S T E P   L O O P  *****
      !=================================================================      
 
#if defined( USE_TIMERS )
      CALL GEOS_Timer_Start( "Timesteps", RC )
#endif

      ! Echo message before first timestep
      WRITE( 6, '(a)' )
      WRITE( 6, '(a)' ) REPEAT( '*', 44 )
      WRITE( 6, '(a)' ) '* B e g i n   T i m e   S t e p p i n g !! *'
      WRITE( 6, '(a)' ) REPEAT( '*', 44 )
      WRITE( 6, '(a)' ) 

      ! NSTEP is the number of dynamic timesteps w/in the outer loop
      ! Timesteps are now retrieved in seconds (ewl, 2/6/2018)
      N_DYN_STEPS = 10800 / GET_TS_DYN()     ! 3hr interval

      ! Start a new outer loop
      DO 

      ! Compute time parameters at start of 6-h loop
      CALL Set_Current_Time()

      ! NSECb is # of seconds (measured from 00 GMT today) 
      ! at the start of this 6-hr timestepping loop.
      ! NOTE: Assume we start at the head of each minute (i.e. SECONDS=0)
      HOUR   = GET_HOUR()
      HOUR   = ( HOUR / 6 ) * 6
      MINUTE = GET_MINUTE()
      SECOND = GET_SECOND()
      NSECb  = ( HOUR * 3600 ) + ( MINUTE * 60 ) + SECOND

      ! Get dynamic timestep in seconds
      N_DYN  = GET_TS_DYN()

      !=================================================================
      !     ***** D Y N A M I C   T I M E S T E P   L O O P *****
      !     *****    a k a   H E A R T B E A T   L O O P    *****
      !=================================================================
      DO N_STEP = 1, N_DYN_STEPS
    
#if defined( LINUX_PGI )
         ! Flush the buffer each timestep so we get output
         CALL Flush( 6 )
#endif

         ! Compute & print time quantities at start of dyn step
         CALL Set_Current_Time()
         CALL Print_Current_Time()

         ! Set time variables for dynamic loop
         DAY_OF_YEAR   = GET_DAY_OF_YEAR()
         DAY           = GET_DAY()
         ELAPSED_SEC   = GET_ELAPSED_SEC()
         MONTH         = GET_MONTH()
         NHMS          = GET_NHMS()
         NYMD          = GET_NYMD()
         HOUR          = GET_HOUR()
         MINUTE        = GET_MINUTE()
         SECOND        = GET_SECOND()
         TAU           = GET_TAU()
         YEAR          = GET_YEAR()
         SEASON        = GET_SEASON()
         ELAPSED_TODAY = ( HOUR * 3600 ) + ( MINUTE * 60 ) + SECOND

         !### Debug
         IF ( prtDebug ) THEN
            CALL Debug_Msg( '### MAIN: a SET_CURRENT_TIME' )
         ENDIF

         !--------------------------------------------------------------
         ! %%%%% HISTORY (netCDF diagnostics) %%%%%
         !
         ! Certain Hg simulation diagnostics (e.g. deposition of Hg2
         ! and HgP onto snow and ice) need to be zeroed out at the
         ! start each timestep, before operations like drydep, wetdep,
         ! and convection are executed.  Call a routine to do this.
         ! (bmy, 10/25/18)
         !--------------------------------------------------------------
#if defined( USE_TIMERS )
         CALL GEOS_Timer_Start( "All diagnostics",           RC )
         CALL GEOS_Timer_Start( "=> History (netCDF diags)", RC )
#endif

         IF ( ITS_A_MERCURY_SIM ) THEN
            CALL Reset_Hg_Diags( am_I_Root, Input_Opt, State_Diag, RC  )
            IF ( RC /= GC_SUCCESS ) THEN
               ErrMsg = 'Error encountered in "Reset_Hg_Diags!"'
               CALL Error_Stop( ErrMsg, ThisLoc )
            ENDIF
         ENDIF

#if defined( USE_TIMERS )
         CALL GEOS_Timer_End( "All diagnostics",           RC )
         CALL GEOS_Timer_End( "=> History (netCDF diags)", RC )
#endif

         !==============================================================
         !       ***** R U N   H E M C O   P H A S E   1 *****
         ! 
         !    Phase 1 updates the HEMCO clock and the content of the 
         !    HEMCO data list. This should be done before writing the
         !    diagnostics organized in the HEMCO diagnostics structure,
         !    and before using any of the HEMCO data list fields.
         !    (ckeller, 4/1/15)
         !==============================================================
#if defined( USE_TIMERS )
        CALL GEOS_Timer_Start( "HEMCO", RC ) 
        CALL GEOS_Timer_Start( "Input", RC )
#endif
         
         ! Is it time for emissions?
         TimeForEmis = ITS_TIME_FOR_EMIS()

         ! Run HEMCO Phase 1
         CALL Emissions_Run( am_I_Root, Input_Opt,   State_Met, 
     &                       State_Chm, State_Diag,  TimeForEmis, 
     &                       1,         RC )

         ! Trap potential errors
         IF ( RC /= GC_SUCCESS ) THEN
            ErrMsg = 'Error encountered in "Emissions_Run"!'
            CALL Error_Stop( ErrMsg, ThisLoc )
         ENDIF

#if defined( USE_TIMERS )
         CALL GEOS_Timer_End( "Input", RC )
#endif

         !### Debug
         IF ( prtDebug ) THEN
            CALL Debug_Msg( '### MAIN: a HEMCO PHASE 1' )
         ENDIf

         !==============================================================
         !  ***** W R I T E   H E M C O   D I A G N O S T I C S *****
         !==============================================================
#if defined( USE_TIMERS )
         CALL GEOS_Timer_Start( "All diagnostics",       RC )
         CALL GEOS_Timer_Start( "=> HEMCO diagnostics",  RC )
         CALL GEOS_Timer_Start( "Output",                RC )
#endif

         ! Write HEMCO diagnostics (ckeller, 4/1/15)
         CALL HCOI_GC_WriteDiagn( am_I_Root, Input_Opt, .FALSE., RC )

         ! Trap potential errors
         IF ( RC /= GC_SUCCESS ) THEN
            ErrMsg = 'Error encountered in "HCOI_GC_WriteDiagn" ' //
     &               '(writing HEMCO diagnostics)!'
            CALL Error_Stop( ErrMsg, ThisLoc )
         ENDIF

#if defined( USE_TIMERS )
         CALL GEOS_Timer_End( "All diagnostics",      RC ) 
         CALL GEOS_Timer_End( "=> HEMCO diagnostics", RC )
         CALL GEOS_Timer_End( "HEMCO",                RC )
         CALL GEOS_Timer_End( "Output",               RC )
#endif

         !==============================================================
         !      ***** O B S P A C K   D I A G N O S T I C S *****
         !==============================================================
         IF ( Input_Opt%Do_ObsPack .and. ( ELAPSED_TODAY == 0 ) ) THEN

#if defined( USE_TIMERS )
            CALL GEOS_Timer_Start( "All diagnostics",        RC )
            CALL GEOS_Timer_Start( "=> ObsPack diagnostics",  RC )
#endif
            ! Initialize Obspack for the new day
            CALL ObsPack_Init( am_I_Root, NYMD,       NHMS, 
     &                         Input_Opt, State_Diag, RC    )

            ! Trap potential errors
            IF ( RC /= GC_SUCCESS ) THEN
               ErrMsg = 'Error encountered in "ObsPack_Init"!'
               CALL Error_Stop( ErrMsg, ThisLoc )
            ENDIF

#if defined( USE_TIMERS )
            CALL GEOS_Timer_End( "All diagnostics",         RC )
            CALL GEOS_Timer_End( "=> ObsPack diagnostics",  RC )
#endif
         ENDIF

         ! Original diagnostics
         IF ( ITS_TIME_FOR_BPCH( Input_Opt ) ) THEN

#if defined( BPCH_DIAG )
            !===========================================================
            ! *****  W R I T E   B P C H   D I A G N O S T I C S *****
            !===========================================================
#if defined( USE_TIMERS )
            CALL GEOS_Timer_Start( "All diagnostics",              RC ) 
            CALL GEOS_Timer_Start( "=> Binary punch diagnostics",  RC )
            CALL GEOS_Timer_Start( "Output",                       RC )
#endif  

            ! Set time at end of diagnostic timestep
            CALL SET_DIAGe( TAU )

            ! Write bpch file
            IF ( DO_DIAG_WRITE ) THEN

               ! Write data to the "trac_avg." bpch file
               CALL DIAG3
     &              ( am_I_Root, Input_Opt, State_Met, State_Chm, RC )

               ! Flush file units
               CALL CTM_FLUSH()
            ENDIF

#if defined( USE_TIMERS )
            CALL GEOS_Timer_End( "All diagnostics",              RC ) 
            CALL GEOS_Timer_End( "=> Binary punch diagnostics",  RC )
            CALL GEOS_Timer_End( "Output",                       RC )
#endif

#endif

!------------------------------------------------------------------------------
! Prior to 8/9/18:
! Disable CALC_STE routine, which can give misleading STE info (bmy, 8/9/18)
!            !===========================================================
!            !       ***** S T R A T - T R O P   F L U X E S *****
!            !===========================================================
!            IF ( LCHEM .and. LSCHEM ) THEN
!
!#if defined( USE_TIMERS )
!               CALL GEOS_Timer_Start( "=> Strat chem", RC )
!#endif
!
!               ! Determine and print STE fluxes (ltm, 06/10/12)
!               CALL Calc_STE( am_I_Root, Input_Opt, State_Chm, 
!     &                        State_Met, RC )
!
!               ! Trap potential errors
!               IF ( RC /= GC_SUCCESS ) THEN
!                  ErrMsg = 'Error encountered in "Calc_STE"!'
!                  CALL Error_Stop( ErrMsg, ThisLoc )
!               ENDIF
!
!#if defined( USE_TIMERS )
!               CALL GEOS_Timer_End( "=> Strat chem", RC )
!#endif
!            ENDIF
!------------------------------------------------------------------------------

            !===========================================================
            !    ***** W R I T E   H E M C O   R E S T A R T S *****
            !===========================================================
#if defined( USE_TIMERS )
            CALL GEOS_Timer_Start( "Output", RC )
            CALL GEOS_Timer_Start( "HEMCO",  RC )
#endif

            ! Force the output of a HEMCO restart file (ckeller, 4/1/15)
            CALL HCOI_GC_WriteDiagn( am_I_Root, Input_Opt, .TRUE., RC )

            ! Trap potential errors
            IF ( RC /= GC_SUCCESS ) THEN
               ErrMsg = 'Error encountered in "HCOI_GC_WriteDiagn"!'
               CALL Error_Stop( ErrMsg, ThisLoc )
            ENDIF

#if defined( USE_TIMERS )
            CALL GEOS_Timer_End( "HEMCO",  RC )
            CALL GEOS_Timer_End( "Output", RC )
#endif

            ! Set time at beginning of next diagnostic timestep
            CALL SET_DIAGb( TAU )

#if defined( BPCH_DIAG )
            !===========================================================
            !   ***** Z E R O   B P C H   D I A G N O S T I C S *****
            !===========================================================
            CALL INITIALIZE( am_I_Root, Input_Opt, 2, RC )  ! Zero diag arrays
            CALL INITIALIZE( am_I_Root, Input_Opt, 3, RC )  ! Zero counters
#endif

         ENDIF

         !=============================================================
         !   ***** W R I T E   G T M M   R E S T A R T   F I L E *****
         !     ***** MUST be done after call to diag3 ****
         !
         ! %%%%% NOTE: THIS MAY BE BROKEN %%%%
         !=============================================================
         ! Make land restart file: for GTMM runs only, beginning of each 
         ! month but not start of the run.
         IF ( LGTMM .AND. ITS_A_NEW_MONTH() .AND. NYMD /= NYMDb ) THEN
            IF (.NOT.( ITS_TIME_FOR_BPCH( Input_Opt ) )) THEN

               ! Get the species ID (NNN) from the wetdep ID (N)
               N   = 1
               NNN = State_Chm%Map_Wetdep(N)

               DO
                  ! Exit once we encounter Hg2
                  If ( State_Chm%SpcData(NNN)%Info%Is_Hg2 ) THEN
                     EXIT
                  ENDIF

                  ! Get the species ID (NNN) from the wetdep ID (N)
                  N   = N + 1
                  NNN = State_Chm%Map_Wetdep(N)
               ENDDO
               CALL UPDATE_DEP( N )
            ENDIF
            CALL MAKE_GTMM_RESTART( am_I_Root, Input_Opt, NYMD, 
     &                              NHMS,      TAU,       RC    )
         ENDIF

         !==============================================================
         !       ***** T E S T   F O R   E N D   O F   R U N *****
         !==============================================================
         IF ( ITS_TIME_FOR_EXIT() ) GOTO 9999

         !==============================================================
         !              ***** D A I L Y   D A T A *****
         !==============================================================
         laiMM = -1
         wasModisRead = .FALSE.
         IF ( ITS_A_NEW_DAY() ) THEN 

            ! Find the proper month and year for LAI  based on the current 
            ! month & year.  LAI months begin in the middle of the month.
            CALL Find_Lai_Month( day_of_year, month, year, 
     &                                        laiMM, laiYYYY )
            
            ! Read MODIS LAI from disk (if it is a new LAI month)
            ! Logical flag MODISREAD determines if data was just read in
            CALL Read_Modis_Lai( am_I_Root, Input_Opt,    laiYYYY, 
     &                           laiMM,     wasModisRead, RC       )

            ! Trap potential errors
            IF ( RC /= GC_SUCCESS ) THEN
               ErrMsg = 'Error encountered in "Read_Modis_LaiI"!'
               CALL Error_Stop( ErrMsg, ThisLoc )
            ENDIF

            ! Interpolate LAI quantities from native grid to GEOS-Chem grid
            CALL Compute_Modis_Lai( am_I_Root,    Input_Opt, State_Met, 
     &                              day_of_year,  laiMM,     mapping,
     &                              wasModisRead, RC )

            ! Trap potential errors
            IF ( RC /= GC_SUCCESS ) THEN
               ErrMsg = 'Error encountered in "Compute_Modis_Lai"!'
               CALL Error_Stop( ErrMsg, ThisLoc )
            ENDIF

            !### Debug
            IF ( prtDebug ) CALL Debug_Msg( '### MAIN: a DAILY DATA' )
         ENDIF

         !==============================================================
         !   ***** I N T E R P O L A T E   Q U A N T I T I E S *****   
         !==============================================================

         ! Interpolate I-3 fields to current dynamic timestep, 
         ! based on their values at NSEC and NSEC+N_DYN
         CALL Interp( NSECb, ELAPSED_TODAY,
     &                N_DYN, Input_Opt,     State_Met )

         ! If we are not doing transport, then make sure that
         ! the floating pressure is set to PSC2_WET (bdf, bmy, 8/22/02)
         ! Now also includes PSC2_DRY (ewl, 5/4/16)
         IF ( .not. LTRAN ) THEN
            CALL Set_Floating_Pressures( am_I_Root, State_Met, RC )

            ! Trap potential errors
            IF ( RC /= GC_SUCCESS ) THEN
               ErrMsg = 'Error encountered in "Set_Floating_Pressures"!'
               CALL Error_Stop( ErrMsg, ThisLoc )
            ENDIF
         ENDIF

         ! Compute updated airmass quantities at each grid box 
         ! and update tracer concentration to conserve tracer mass
         ! (ewl, 10/28/15)
         CALL AirQnt( am_I_Root, Input_opt, State_Met, State_Chm, RC, 
     &                Update_Mixing_Ratio=.TRUE. )

         ! Trap potential errors
         IF ( RC /= GC_SUCCESS ) THEN
            ErrMsg = 'Error encountered in "AirQnt (call #2)"!'
            CALL Error_Stop( ErrMsg, ThisLoc )
         ENDIF

         ! SDE 05/28/13: Set H2O to State_Chm tracer if relevant and,
         ! if LUCX=T and LSETH2O=F and LACTIVEH2O=T, update specific humidity 
         ! in the stratosphere
         !
         ! NOTE: Specific humidity may change in SET_H2O_TRAC and
         ! therefore this routine may call AIRQNT again to update
         ! air quantities and tracer concentrations (ewl, 10/28/15)
         IF ( ITS_A_FULLCHEM_SIM .and. id_H2O > 0 ) THEN
            CALL SET_H2O_TRAC( am_I_Root, 
     &                         ( ( .not. LUCX ) .or. LSETH2O ),
     &                         Input_Opt, State_Met, State_Chm, RC )

            ! Trap potential errors
            IF ( RC /= GC_SUCCESS ) THEN
               ErrMsg = 'Error encountered in "Set_H2O_Trac" #1!'
               CALL Error_Stop( ErrMsg, ThisLoc )
            ENDIF

            ! Only force strat once if using UCX
            IF (LSETH2O) LSETH2O = .FALSE.
         ENDIF

         ! Compute the cosine of the solar zenith angle array
         ! State_Met%SUNCOS     = at the current time
         ! State_Met%SUNCOSmid  = at the midpt of the chem timestep
         ! State_Met%SUNCOSmid5 = at the midpt of the chem timestep 5hrs ago
         CALL Get_Cosine_SZA( am_I_Root, Input_Opt, State_Met, RC )

         ! Trap potential errors
         IF ( RC /= GC_SUCCESS ) THEN
            ErrMsg = 'Error encountered in "Get_Cosine_SZA"!'
            CALL Error_Stop( ErrMsg, ThisLoc )
         ENDIF

#if defined( BPCH_DIAG )
         !--------------------------------------------------------------
         !   ***** B P C H   D I A G N O S T I C   U P D A T E *****
         !--------------------------------------------------------------
#if defined( USE_TIMERS )
         CALL GEOS_Timer_Start( "All diagnostics",              RC )
         CALL GEOS_Timer_Start( "=> Binary punch diagnostics",  RC )
#endif

         !-------------------------------------------------------------
         ! ND55 (bpch) DIAGNOSTIC
         !
         ! Tropopause level, height [km], and pressure [hPa]
         !
         ! NOTE: Attach here because this is where the call to
         ! obsolete routine DIAG_TROPOPAUSE was.  This should
         ! replicate results w/ prior versions. (bmy, 1/17/18)
         !-------------------------------------------------------------
         IF ( ND55 > 0 ) THEN
            AD55(:,:,1) = AD55(:,:,1) + State_Met%TropLev
            AD55(:,:,2) = AD55(:,:,2) + State_Met%TropHt
            AD55(:,:,3) = AD55(:,:,3) + State_Met%TropP
         ENDIF

         ! Update dynamic timestep (for bpch diagnostics
         CALL SET_CT_DYN( INCREMENT=.TRUE. )

         !-------------------------------------------------------------
         ! Get averaging intervals for local-time diagnostics
         ! (NOTE: maybe improve this later on)
         ! Placed after interpolation to get correct value of TROPP. 
         ! (ccc, 12/9/08)
         !-------------------------------------------------------------
         CALL DIAG_2PM( State_Met )

#if defined( USE_TIMERS )
         CALL GEOS_Timer_End( "All diagnostics",              RC )
         CALL GEOS_Timer_End( "=> Binary punch diagnostics",  RC )
#endif

#endif

         !### Debug
         IF ( prtDebug ) CALL Debug_Msg( '### MAIN: a INTERP, etc' )

         !----------------------------------------------------------
         ! %%% GET SOME NON-EMISSIONS DATA FIELDS VIA HEMCO %%%
         !
         ! HEMCO can track non-emission data fields for chemistry
         ! simulations.  Put these subroutine calls after the 
         ! call to EMISSIONS_RUN, so that the HEMCO data structure
         ! will be initialized. (bmy, 3/20/15)
         !
         ! HEMCO data list is now updated further above, so can
         ! take these calls out of the emissions sequence.
         ! (ckeller, 4/01/15) 
         !----------------------------------------------------------
         IF ( LCHEM .and. ITS_A_NEW_MONTH() ) THEN
 
            ! The following only apply when photolysis is used,
            ! that is for fullchem or aerosol simulations.
            IF ( ITS_A_FULLCHEM_SIM  .or. ITS_AN_AEROSOL_SIM ) THEN

               ! Copy UV Albedo data (for photolysis) into the
               ! State_Met%UVALBEDO field. (bmy, 3/20/15)
               CALL Get_UvAlbedo( am_I_Root, Input_Opt, State_Met, RC )

               ! Trap potential errors
               IF ( RC /= GC_SUCCESS ) THEN
                  ErrMsg = 'Error encountered in "Get_UvAlbedo"!'
                  CALL Error_Stop( ErrMsg, ThisLoc )
               ENDIF

               IF ( Input_Opt%USE_TOMS_O3 ) THEN
                  ! Get TOMS overhead O3 columns for photolysis from
                  ! the HEMCO data structure (bmy, 3/20/15)
                  CALL Read_TOMS( am_I_Root, Input_Opt, RC )

                  ! Trap potential errors
                  IF ( RC /= GC_SUCCESS ) THEN
                     ErrMsg = 'Error encountered in "Read_TOMS"!'
                     CALL Error_Stop( ErrMsg, ThisLoc )
                  ENDIF
               ENDIF

            ENDIF

            ! Read data required for Hg2 gas-particle partitioning 
            ! (H Amos, 25 Oct 2011)
            IF ( ITS_A_MERCURY_SIM ) THEN 
               CALL Read_Hg2_Partitioning
     &            ( am_I_Root, Input_Opt, State_Met, MONTH, RC )

               ! Trap potential errors
               IF ( RC /= GC_SUCCESS ) THEN
                  ErrMsg = 
     &               'Error encountered in "Read_Hg2_Partitioning"!'
                  CALL Error_Stop( ErrMsg, ThisLoc )
               ENDIF

               IF ( prtDebug ) THEN
                 CALL Debug_Msg( '### MAIN: a READ_HG2_PARTITIONING')
               ENDIF
            ENDIF
         ENDIF

         ! Prescribe methane surface concentrations throughout PBL
         IF ( ITS_A_FULLCHEM_SIM .and. id_CH4 > 0 ) THEN

            ! Debug print
            IF ( prtDebug ) THEN 
               CALL DEBUG_MSG( '### MAIN: Setting PBL CH4 conc')
            ENDIF
            
            ! Set CH4 concentrations
            CALL SET_CH4( am_I_Root, Input_Opt,  State_Met,
     &                    State_Chm, State_Diag, RC )

            ! Trap potential errors
            IF ( RC /= GC_SUCCESS ) THEN
               ErrMsg = 'Error encountered in call to "SET_CH4"!'
               CALL Error_Stop( ErrMsg, ThisLoc )
            ENDIF
         ENDIF

         !==============================================================
         !              ***** T R A N S P O R T *****
         !==============================================================
         IF ( ITS_TIME_FOR_DYN() ) THEN

#if defined( USE_TIMERS )
            CALL GEOS_Timer_Start( "Transport", RC )
#endif

#if defined( BPCH_TPBC )
            ! Output BC's
            ! Save boundary conditions (global grid) for future nested run
            IF ( LWINDO ) THEN
               CALL Save_Global_Tpcore_BC
     &            ( am_I_Root, Input_Opt, State_Chm, RC )

               ! Trap potential errors
               IF ( RC /= GC_SUCCESS ) THEN
                  ErrMsg = 
     &               'Error encountered in "Save_Global_Tpcore_BC"!'
                  CALL Error_Stop( ErrMsg, ThisLoc )
               ENDIF

            ENDIF
#endif

!-----------------------------------------------------------------------------
! TWO-WAY NESTING OPTION
! This is only invoked when compiling GEOS-Chem with COUPLE=y
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! %%%%% NOTE: THIS OPTION WAS PROBABLY BROKEN WITH THE ADDITION OF HEMCO %%%%%
! %%%%% AND FLEXCHEM INTO GEOS-CHEM V10-01 AND v11-01.  BUYER BEWARE.    %%%%%
! %%%%%  --  Bob Yantosca (22 Jan 2018)                                  %%%%%
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
#if defined( EXCHANGE )

#if defined( GRID4x5 ) || defined( GRID2x25 )

            ! Exchange the position of POST (global simulations)
            IF ( ITS_TIME_FOR_EXCHANGE() ) THEN 
               CALL EXCHANGE_GLOBAL_POST
     &              ( am_I_Root, Input_Opt, State_Met, State_Chm, RC )
            ENDIF

#elif defined( NESTED_CH ) || defined( NESTED_NA ) || defined( NESTED_EU )

            ! Exchange the position of POST (nested-grid simulations)
            IF ( ITS_TIME_FOR_EXCHANGE() ) THEN
               CALL EXCHANGE_NESTED_POST
     &              ( am_I_Root, Input_Opt, State_Met, State_Chm, RC )
            ENDIF

#endif

#endif
!-----------------------------------------------------------------------------

            ! Call the appropriate version of TPCORE
            IF ( LTRAN ) THEN 
               CALL Do_Transport( am_I_Root, Input_Opt,  State_Met,
     &                            State_Chm, State_Diag, RC         )

               ! Trap potential error
               IF ( RC /= GC_SUCCESS ) THEN
                  ErrMsg = 'Error encountered in "Do_Transport"!'
                  CALL Error_Stop( ErrMsg, ThisLoc )
               ENDIF

               !### Debug
               IF ( prtDebug ) THEN 
                  CALL Debug_Msg( '### MAIN: a DO_TRANSPORT' )
               ENDIF
            ENDIF


            ! Initialize wet scavenging and wetdep fields after
            ! the airmass quantities are reset after transport
#if defined ( TOMAS )
            ! ... TOMAS microphysics: Always call SETUP_WETSCAV ...
            CALL Setup_WetScav( am_I_Root, Input_Opt, 
     &                          State_Met, State_Chm, RC )

               ! Trap potential error
               IF ( RC /= GC_SUCCESS ) THEN
                  ErrMsg = 
     &               'Error encountered in "Setup_WetScav" (TOMAS)!'
                  CALL Error_Stop( ErrMsg, ThisLoc )
               ENDIF

#else 

            ! ... Standard GEOS-Chem: Call INIT_WETSCAV if   ...
            ! ... convection or wet scavenging or chemistry are turned on ...
            IF ( LCONV .or. LWETD .or. LCHEM ) THEN
               CALL Setup_WetScav( am_I_Root, Input_Opt, 
     &                             State_Met, State_Chm, RC )

               ! Trap potential error
               IF ( RC /= GC_SUCCESS ) THEN
                  ErrMsg = 'Error encountered in "Setup_WetScav"!'
                  CALL Error_Stop( ErrMsg, ThisLoc )
               ENDIF

               !### Debug
               IF ( prtDebug ) THEN 
                  CALL Debug_Msg( '### MAIN: a SETUP_WETSCAV' )
               ENDIF
            ENDIF

#endif

#if defined( USE_TIMERS )
            CALL GEOS_Timer_End( "Transport", RC )
#endif

         ENDIF

         ! Update age of air
         CALL Set_Met_AgeOfAir( State_Met )

         !==============================================================
         !     ***** C O M P U T E   P B L   H E I G H T  etc. *****
         !==============================================================
#if defined( USE_TIMERS )
         CALL GEOS_Timer_Start( "Boundary layer mixing", RC )
#endif

         ! Move this call from the PBL mixing routines because the PBL
         ! height is used by drydep and some of the emissions routines.
         ! (ckeller, 3/5/15) 
         CALL Compute_PBL_Height( am_I_Root, State_Met, RC )

         ! Trap potential errors
         IF ( RC /= GC_SUCCESS ) THEN
            ErrMsg = 'Error encountered in "Compute_PBL_Height"!'
            CALL Error_Stop( ErrMsg, ThisLoc )
         ENDIF

#if defined( USE_TIMERS )
         CALL GEOS_Timer_End( "Boundary layer mixing", RC )
#endif

         !### Debug
         IF ( prtDebug ) THEN 
            CALL Debug_Msg( '### MAIN: a COMPUTE_PBL_HEIGHT' )
         ENDIF

         !--------------------------------------------------------------
         ! Test for emission timestep
         ! Now always do emissions here, even for full-mixing
         ! (ckeller, 3/5/15)
         !--------------------------------------------------------------
         IF ( ITS_TIME_FOR_EMIS() ) THEN

#if defined( BPCH_DIAG )
            ! Increment emission counter
            CALL Set_Ct_Emis( INCREMENT=.TRUE. )
#endif

            !===========================================================
            !         ***** D R Y   D E P O S I T I O N *****
            !===========================================================
            IF ( LDRYD ) THEN

#if defined( USE_TIMERS )
               CALL GEOS_Timer_Start( "Dry deposition", RC )
#endif
               ! Compute drydep velocities
               CALL Do_Drydep( am_I_Root, Input_Opt,  State_Met, 
     &                         State_Chm, State_Diag, RC )

               ! Trap potential errors
               IF ( RC /= GC_SUCCESS ) THEN
                  ErrMsg = 'Error encountered in "Do_Drydep!"!'
                  CALL Error_Stop( ErrMsg, ThisLoc )
               ENDIF

#if defined( USE_TIMERS )
               CALL GEOS_Timer_End ( "Dry deposition", RC )
#endif

               !### Debug
               IF ( prtDebug ) THEN 
                  CALL Debug_Msg( '### MAIN: a DO_DRYDEP' )
               ENDIF
            ENDIF

            !===========================================================
            !             ***** E M I S S I O N S *****
            !
            ! NOTE: For a complete description of how emissions from
            ! HEMCO are added into GEOS-Chem (and how they are mixed
            ! into the boundary layer), please see the wiki page:
            !
            ! http://wiki-geos-chem.org/Distributing_emissions_in_the_PBL
            !===========================================================
#if defined( USE_TIMERS )
            CALL GEOS_Timer_Start( "HEMCO", RC )
#endif
            
            ! Is it time for emissions?
            TimeForEmis = ITS_TIME_FOR_EMIS()

            ! EMISSIONS_RUN will call HEMCO run phase 2. HEMCO run phase
            ! only calculates emissions. All data has been read to disk
            ! in phase 1 at the beginning of the time step. 
            ! (ckeller, 4/1/15)
            CALL Emissions_Run( am_I_Root, Input_Opt,   State_Met, 
     &                          State_Chm, State_Diag,  TimeForEmis, 
     &                          2,         RC )

            ! Trap potential errors
            IF ( RC /= GC_SUCCESS ) THEN
               ErrMsg = 
     &           'Error encountered in "Emissions_Run"! after drydep!'
               CALL Error_Stop( ErrMsg, ThisLoc )
            ENDIF

            !### Debug
            IF ( prtDebug ) THEN 
               CALL Debug_Msg( '### MAIN: a HEMCO PHASE 2' )
            ENDIF

#if defined( USE_TIMERS )
            CALL GEOS_Timer_End( "HEMCO", RC )
#endif
         ENDIF

         !-------------------------------
         ! Test for convection timestep
         !-------------------------------
         IF ( ITS_TIME_FOR_CONV() ) THEN
            
#if defined( BPCH_DIAG )
            ! Increment the convection timestep
            CALL Set_Ct_Conv( INCREMENT=.TRUE. )
#endif

            !===========================================================
            !      ***** M I X E D   L A Y E R   M I X I N G *****
            !===========================================================
#if defined( USE_TIMERS )
            CALL GEOS_Timer_Start( "Boundary layer mixing", RC )
#endif

            ! Note: mixing routine expects tracers in v/v
            ! DO_MIXING applies the tracer tendencies (dry deposition,
            ! emission rates) to the tracer arrays and performs PBL
            ! mixing. 
            ! In the non-local PBL scheme, dry deposition and emission
            ! fluxes below the PBL are handled within the PBL mixing
            ! routine. Otherwise, tracer concentrations are first updated
            ! and the full-mixing is then applied.
            ! (ckeller, 3/5/15)
            ! NOTE: Tracer concentration units are converted locally
            ! to [v/v dry air] for mixing. Eventually mixing should
            ! be updated to use [kg/kg total air] (ewl, 9/18/15)
            CALL Do_Mixing( am_I_Root, Input_Opt,  State_Met,
     &                      State_Chm, State_Diag, RC         )

            ! Trap potential errors
            IF ( RC /= GC_SUCCESS ) THEN
               ErrMsg = 'Error encountered in "Do_Mixing"!'
               CALL ERror_Stop( ErrMsg, ThisLoc )
            ENDIF

#if defined( USE_TIMERS )
            CALL GEOS_Timer_End( "Boundary layer mixing", RC )
#endif

            IF ( prtDebug ) CALL Debug_Msg( '### MAIN: a TURBDAY:2' )

            !===========================================================
            !        ***** C L O U D   C O N V E C T I O N *****
            !===========================================================
            IF ( LCONV ) THEN

#if defined( USE_TIMERS )
               CALL GEOS_Timer_Start( "Convection", RC )
#endif

               ! Call the appropriate convection routine
               ! NOTE: Tracer concentration units are converted locally
               ! to [kg/kg total air] for convection (ewl, 9/18/15)
               CALL Do_Convection( am_I_Root, Input_Opt,  State_Met,
     &                             State_Chm, State_Diag, RC         )

               ! Trap potential errors
               IF ( RC /= GC_SUCCESS ) THEN
                  ErrMsg = 'Error encountered in "Do_Convection"!'
                  CALL Error_Stop( ErrMsg, ThisLoc )
               ENDIF

               !### Debug
               IF ( prtDebug ) THEN
                  CALL Debug_Msg( '### MAIN: a CONVECTION' )
               ENDIF

#if defined( USE_TIMERS )
               CALL GEOS_Timer_End( "Convection", RC )
#endif
            ENDIF 

         ENDIF 

         !==============================================================
         !               ***** C H E M I S T R Y *****
         !============================================================== 
#if defined( USE_TIMERS )
         CALL GEOS_Timer_Start( "All chemistry", RC )
#endif

         ! Get the overhead column O3 for use with FAST-J
         ! NOTE: Move to CHEMISTRY section.  This now has to come after
         ! the call to HEMCO emissions driver EMISSIONS_RUN. (bmy, 3/20/15)
         CALL Get_Overhead_O3_For_FastJ( am_I_Root )
         
         ! Every chemistry timestep...
         IF ( ITS_TIME_FOR_CHEM() ) THEN

#if defined( BPCH_DIAG )
            ! Increment chemistry timestep counter
            CALL Set_Ct_Chem( INCREMENT=.TRUE. )
#endif

            ! SDE 05/28/13: Set H2O to State_Chm tracer if relevant
            IF ( ITS_A_FULLCHEM_SIM .and. id_H2O > 0 ) THEN
               CALL Set_H2O_Trac( am_I_Root, (.not. LUCX), Input_Opt ,
     &                            State_Met, State_Chm, RC )

               ! Trap potential errors
               IF ( RC /= GC_SUCCESS ) THEN
                  ErrMsg = 'Error encountered in "Set_H2O_Trac" #2!'
                  CALL Error_Stop( ErrMsg, ThisLoc )
               ENDIF
            ENDIF

            ! Do GEOS-Chem chemistry
            ! NOTE: Tracer concentration units are converted locally
            ! to [kg] for all of chemistry. We will replace use of [kg] 
            ! once FlexChem is implemented (ewl, 9/18/15)
            CALL Do_Chemistry( am_I_Root, Input_Opt,  State_Met, 
     &                         State_Chm, State_Diag, RC         ) 

            ! Trap potential errors
            IF ( RC /= GC_SUCCESS ) THEN
               ErrMsg = 'Error encountered in "Do_Chemistry"!'
               CALL Error_Stop( ErrMsg, ThisLoc )
            ENDIF

         ENDIF 

#if defined( USE_TIMERS )
         CALL GEOS_Timer_End( "All chemistry", RC )
#endif

         !==============================================================
         ! ***** W E T   D E P O S I T I O N  (rainout + washout) *****
         !==============================================================
         IF ( LWETD .and. ITS_TIME_FOR_DYN() ) THEN

#if defined( USE_TIMERS )
            CALL GEOS_Timer_Start( "Wet deposition", RC )
#endif

            ! Do wet deposition
            ! NOTE: Tracer concentration units are converted locally
            ! to [kg/m2] in wet deposition to enable calculations
            ! along the column (ewl, 9/18/15)
            CALL Do_WetDep( am_I_Root, Input_Opt,  State_Met, 
     &                      State_Chm, State_Diag, RC         )

            ! Trap potential errors
            IF ( RC /= GC_SUCCESS ) THEN
               ErrMsg = 'Error encountered in "Do_WetDep"!'
               CALL Error_Stop( ErrMsg, ThisLoc )
            ENDIF

#if defined( USE_TIMERS )
            CALL GEOS_Timer_End( "Wet deposition", RC )
#endif

         ENDIF

         !==============================================================
         !      ***** U P D A T E   O P T I C A L   D E P T H *****          
         !==============================================================
         IF ( ITS_TIME_FOR_CHEM() ) THEN

#if defined( USE_TIMERS )
            CALL GEOS_Timer_Start( "All chemistry",       RC )
            CALL GEOS_Timer_Start( "=> All aerosol chem", RC )
#endif

            ! Recalculate the optical depth at the wavelength(s) specified
            ! in the Radiation Menu. This must be done before the call to any
            ! diagnostic and only on a chemistry timestep.
            ! (skim, 02/05/11)
            CALL Recompute_OD( am_I_Root, Input_Opt,  State_Met,
     &                         State_Chm, State_Diag, RC         )

            ! Trap potential errors
            IF ( RC /= GC_SUCCESS ) THEN
               ErrMsg = 'Error encountered in "Recompute_OD"!'
               CALL Error_Stop( ErrMsg, ThisLoc )
            ENDIF

#if defined( USE_TIMERS )
            CALL GEOS_Timer_End( "All chemistry",       RC )
            CALL GEOS_Timer_End( "=> All aerosol chem", RC )
#endif
         ENDIF

#if defined( RRTMG )
         !==============================================================
         !  ***** R R T M G   R A D I A T I V E   T R A N S F E R *****
         !
         ! NOTE: Tracer concentration units are converted locally to
         ! [kg] in RRTMG. Units should eventually be [kg/kg]  
         ! (ewl, 9/18/15)
         !==============================================================
         IF ( Input_opt%LRAD .and. ITS_TIME_FOR_RT() ) THEN

#if defined( USE_TIMERS )
            CALL GEOS_Timer_Start( "RRTMG", RC )
#endif

            ! Splash page
            IF ( am_I_Root .and. FIRST_RT ) THEN
               WRITE( 6, '(a)' ) REPEAT( '#', 79 )
               WRITE( 6, 500 )     
     &              'R R T M G : Radiative Transfer Model (by AER)'
 500           FORMAT( '#####', 12x, a, 12x, '#####' )
               WRITE( 6, '(a)' ) REPEAT( '#', 79 )
               FIRST_RT = .FALSE.
            ENDIF

            iSeed = iSeed + 15

            !-----------------------------------------------------------
            ! Determine if we are doing clear-sky or all-sky
            !
            ! Clear-sky is output with all-sky, so we just need
            ! to run once regardless of whether both are required
            ! or just one.
            !-----------------------------------------------------------
            IF ( Input_Opt%LSKYRAD(2) ) THEN
               iCld = 1
            ELSE
               iCld = 0         !clouds are on
            ENDIF

            !-----------------------------------------------------------
            ! Calculation for each of the potential output types
            ! See: wiki.geos-chem.org/Coupling_GEOS-Chem_with_RRTMG
            !
            ! Flux outputs (scheduled in HISTORY.rc):
            !  0-BA  1=O3  2=ME  3=SU   4=NI  5=AM
            !  6=BC  7=OA  8=SS  9=DU  10=PM  11=ST (UCX only)
            !-----------------------------------------------------------
            DO N = 1, State_Diag%nRadFlux

               ! Index number for RRTMG (see list above)
               iSpecMenu = State_Diag%RadFluxInd(N)

               ! Slot # of netCDF diagnostic arrays to update
               iNcDiag = N

               ! Echo info
               WRITE( 6, 520 ) State_Diag%RadFluxName(N), iSpecMenu
 520           FORMAT( 5x, '- Calling RRTMG to compute flux: ',
     &                 a2, ' (Index = ', i2.2, ')' )

               ! Generate mask for species in RT
               CALL Set_SpecMask( iSpecMenu )

               ! Compute radiative fluxes for the given output
               CALL Do_RRTMG_Rad_Transfer( am_I_Root  = am_I_Root,
     &                                     ThisDay    = Day,
     &                                     ThisMonth  = Month,
     &                                     iCld       = iCld,
     &                                     iSpecMenu  = iSpecMenu,
     &                                     iNcDiag    = iNcDiag,
     &                                     iSeed      = iSeed,
     &                                     Input_Opt  = Input_Opt,
     &                                     State_Met  = State_Met,
     &                                     State_Chm  = State_Chm,
     &                                     State_Diag = State_Diag,
     &                                     RC         = RC          )

               ! Trap potential errors
               IF ( RC /= GC_SUCCESS ) THEN
                  ErrMsg = 'Error encountered in '      //
     &                     '"Do_RRTMG_Rad_Transfer", '  //
     &                     'for flux output = '         //
     &                 TRIM( State_Diag%RadFluxName(N) )
                  CALL Error_Stop( ErrMsg, ThisLoc )
               ENDIF
            ENDDO

#if defined( BPCH_DIAG )
            ! Increment radiation timestep counter
            CALL Set_Ct_Rad( INCREMENT=.TRUE. )
#endif

            !### Debug
            IF ( prtDebug ) THEN 
               CALL Debug_Msg( '### MAIN: a DO_RRTMG_RAD_TRANSFER' )
            ENDIF

#if defined( USE_TIMERS )
            CALL GEOS_Timer_End( "RRTMG", RC )
#endif
         ENDIF
#endif

#if defined( DIAG_DEVEL )
         ! Eventually write dobson columns into diagnostics array
         CALL CalcDobsonColumn ( am_I_Root, Input_Opt,
     &                           State_Met, State_Chm, RC )
#endif

         !==============================================================
         !      ***** D I A G N O S T I C S   A R C H I V A L *****
         !==============================================================
 
         !--------------------------------------------------------------
         !           ***** H I S T O R Y   U P D A T E  *****
         !--------------------------------------------------------------
#if defined( USE_TIMERS )
         CALL GEOS_Timer_Start( "All diagnostics",           RC ) 
         CALL GEOS_Timer_Start( "Output",                    RC )
         CALL GEOS_Timer_Start( "=> History (netCDF diags)", RC )
#endif

         ! Set State_Diag arrays that rely on state at end of timestep
         CALL Set_Diagnostics_EndofTimestep( am_I_Root,  Input_Opt,  
     &                                       State_Met,  State_Chm, 
     &                                       State_Diag, RC         )

         ! Trap potential errors
         IF ( RC /= GC_SUCCESS ) THEN
            ErrMsg = 'Error encountered in ' // 
     &               '"Set_Diagnostics_EndOfTimestep"!'
            CALL Error_Stop( ErrMsg, ThisLoc )
         ENDIF

         ! Archive aerosol mass and PM2.5 diagnostics
         IF ( State_Diag%Archive_AerMass ) THEN
            CALL Set_AerMass_Diagnostic( am_I_Root,  Input_Opt,  
     &                                   State_Met,  State_Chm, 
     &                                   State_Diag, RC         )
         ENDIF

         ! Trap potential errors
         IF ( RC /= GC_SUCCESS ) THEN
            ErrMsg = 'Error encountered in ' // 
     &               '"Set_AerMass_Diagnostic"!'
            CALL Error_Stop( ErrMsg, ThisLoc )
         ENDIF

         ! Update each HISTORY ITEM from its data source
         CALL History_Update( am_I_Root, RC )

         ! Trap potential errors
         IF ( RC /= GC_SUCCESS ) THEN
            ErrMsg = 'Error encountered in "History_Update"!'
            CALL Error_Stop( ErrMsg, ThisLoc )
         ENDIF

         ! Increment the timestep values by the heartbeat time
         ! This is becasue we need to write out data with the timestamp
         ! at the end of the heartbeat timestep (i.e. at end of run)
         CALL History_SetTime( am_I_Root, RC )

         ! Trap potential errors
         IF ( RC /= GC_SUCCESS ) THEN
            ErrMsg = 'Error encountered in "History_SetTime"!'
            CALL Error_Stop( ErrMsg, ThisLoc )
         ENDIF

#if defined( USE_TIMERS )
         CALL GEOS_Timer_End( "All diagnostics",           RC ) 
         CALL GEOS_Timer_End( "Output",                    RC )
         CALL GEOS_Timer_End( "=> History (netCDF diags)", RC )
#endif

         !--------------------------------------------------------------
         !      ***** O B S P A C K   D I A G N O S T I C S *****
         !--------------------------------------------------------------
         IF ( Input_Opt%Do_ObsPack ) THEN

#if defined( USE_TIMERS )
            CALL GEOS_Timer_Start( "All diagnostics",        RC ) 
            CALL GEOS_Timer_Start( "=> ObsPack diagnostics", RC )
#endif

            ! Sample the observations in today's ObsPack file
            CALL ObsPack_Sample( am_I_Root,  NYMD,      NHMS, 
     &                           Input_Opt,  State_Met, State_Chm, 
     &                           State_Diag, RC                    )

#if defined( USE_TIMERS )
            CALL GEOS_Timer_End( "All diagnostics",        RC ) 
            CALL GEOS_Timer_End( "=> ObsPack diagnostics", RC )
#endif

         ENDIF

#if defined( BPCH_DIAG ) || defined( BPCH_TIMESER )
         !--------------------------------------------------------------
         !   ***** P L A N E F L I G H T   D I A G   S E T U P  *****
         !
         ! %%%%% NOTE: Works whether BPCH_DIAG=y or n        %%%%%
         ! %%%%% but for timing purposes, count this as bpch %%%%%
         ! %%%%% because it is not computed with History.    %%%%%
         !--------------------------------------------------------------
         IF ( ND40 > 0 .and. DO_DIAG_WRITE ) THEN 
            IF ( ITS_A_NEW_DAY() ) THEN

#if defined( USE_TIMERS )
               CALL GEOS_Timer_Start( "All diagnostics",             RC) 
               CALL GEOS_Timer_Start( "Output",                      RC)
               CALL GEOS_Timer_Start( "=> Binary punch diagnostics", RC)
#endif      

               ! Initialize planeflight diagnostic
               CALL Setup_PlaneFlight( am_I_Root, Input_Opt, 
     &                                 State_Chm, State_Met, RC )

               ! Trap potential errors
               IF ( RC /= GC_SUCCESS ) THEN
                  ErrMsg = 'Error encountered in "Setup_Planeflight"!'
                  CALL Error_Stop( ErrMsg, ThisLoc )
               ENDIF

#if defined( USE_TIMERS )
               CALL GEOS_Timer_End( "All diagnostics",             RC ) 
               CALL GEOS_Timer_End( "Output",                      RC )
               CALL GEOS_Timer_End( "=> Binary punch diagnostics", RC )
#endif      

            ENDIF
         ENDIF
#endif

#if defined( BPCH_DIAG )
         !--------------------------------------------------------------
         ! **** C H 4   S I M U L A T I O N   D I A G N O S I C S ****
         !
         ! CH4 columns from the GOSAT instrument (mps, 6/16/17)
         !--------------------------------------------------------------
         IF ( Input_Opt%GOSAT_CH4_OBS ) THEN
            IF ( ITS_A_NEW_HOUR() ) THEN
      
!------------------------------------------------------------------------------
! For adjoint model only:
!            ! Track cost function contributions
!            CF_PRIOR = COST_FUNC
!------------------------------------------------------------------------------

               CALL CALC_GOSAT_CH4_FORCE( Input_Opt, State_Chm,
     &                                    State_Met, COST_FUNC )

!------------------------------------------------------------------------------
! For adjoint model only:
!            ! Track cost function contributions
!            CF_GOSCH4 = CF_GOSCH4 + COST_FUNC - CF_PRIOR
!------------------------------------------------------------------------------
            ENDIF
         ENDIF

         !--------------------------------------------------------------
         ! **** C H 4   S I M U L A T I O N   D I A G N O S I C S ****
         !
         ! CH4 columns from the TCCON instrument (mps, 8/17/17)
         !--------------------------------------------------------------
         IF ( Input_Opt%TCCON_CH4_OBS ) THEN
            IF ( ITS_A_NEW_HOUR() ) THEN
      
!------------------------------------------------------------------------------
! For adjoint model only:
!            ! Track cost function contributions
!            CF_PRIOR = COST_FUNC
!------------------------------------------------------------------------------

               CALL CALC_TCCON_CH4_FORCE( Input_Opt, State_Chm,
     &                                    State_Met, COST_FUNC )

!------------------------------------------------------------------------------
! For adjoint model only:
!            ! Track cost function contributions
!            CF_GOSCH4 = CF_GOSCH4 + COST_FUNC - CF_PRIOR
!------------------------------------------------------------------------------
            ENDIF
         ENDIF

#endif

         !==============================================================
         !   ***** I N C R E M E N T   E L A P S E D   T I M E *****
         !
         ! Moved before diagnostics to count the last timestep as done.
         ! Need to save timestamps for filenames. (ccc, 5/13/09)
         !============================================================== 
         CALL Timestamp_Diag()
         CALL Set_Elapsed_Sec()
         CALL Set_Current_Time()
         IF ( prtDebug ) THEN
            CALL Debug_Msg( '### MAIN: after SET_ELAPSED_SEC' )
         ENDIF

#if defined( BPCH_DIAG ) || defined( BPCH_TIMESER )
         !==============================================================
         !      ***** D I A G N O S T I C S   A R C H I V A L *****
         !
         !                ***** C O N T I N U E D *****
         !==============================================================
#if defined( USE_TIMERS )
         CALL GEOS_Timer_Start( "All diagnostics",              RC )
         CALL GEOS_Timer_Start( "Output",                       RC )
         CALL GEOS_Timer_Start( "=> Binary punch diagnostics",  RC )
#endif

         !--------------------------------------------------------------
         !   ***** P L A N E F L I G H T   D I A G N O S T I C  *****
         !
         ! Plane following diagnostic, moved here to be consistent
         ! with species array and not be biased by emissions, kyu 30
         ! jan 2015
         !
         ! %%%%% NOTE: works whether BPCH_DIAG=y or n %%%%%
         !--------------------------------------------------------------
         IF ( ND40 > 0 ) THEN 
            ! Archive data along the flight track
            CALL PLANEFLIGHT( am_I_Root, Input_Opt,  State_Met,
     &                        State_Chm, State_Diag, RC )

            ! Trap potential errors
            IF ( RC /= GC_SUCCESS ) THEN
               ErrMsg = 'Error encountered in "Planeflight"!'
               CALL Error_Stop( ErrMsg, ThisLoc )
            ENDIF
         ENDIF
         IF ( prtDebug ) CALL Debug_Msg( '### MAIN: after DIAG40' )

         !--------------------------------------------------------------
         !  ***** A R C H I V E   B P C H   D I A G N O S T I C S *****
         !--------------------------------------------------------------
         IF ( ITS_TIME_FOR_DIAG() ) THEN

            !### Debug
            IF ( prtDebug ) CALL Debug_Msg( '### MAIN: b DIAGNOSTICS' )

            ! Accumulate several diagnostic quantities
            CALL Diag1( am_I_Root, Input_Opt, State_Met, State_Chm, RC )

            ! Trap potential errors
            IF ( RC /= GC_SUCCESS ) THEN
               ErrMsg = 'Error encountered in "Diag1"!'
               CALL Error_Stop( ErrMsg, ThisLoc )
            ENDIF

            IF ( prtDebug ) CALL Debug_Msg( '### MAIN: after DIAG1' )

#if defined( BPCH_DIAG )
            ! ND41: save PBL height in 1200-1600 LT (amf)
            ! (for comparison w/ Holzworth, 1967)
            IF ( ND41 > 0 ) CALL DIAG41()
            IF ( prtDebug ) CALL Debug_Msg( '### MAIN: after DIAG41' )

            ! ND42: SOA concentrations [ug/m3]
            IF ( ND42 > 0 ) THEN
               CALL Diag42( am_I_Root, Input_Opt, State_Met,
     &                      State_Chm, RC                    )

               ! Trap potential errors
               IF ( RC /= GC_SUCCESS ) THEN
                  ErrMsg = 'Error encountered in "Diag42"!'
                  CALL Error_Stop( ErrMsg, ThisLoc )
               ENDIF
            ENDIF
            IF ( prtDebug ) CALL Debug_Msg( '### MAIN: after DIAG42' )
#endif

            ! 24-hr timeseries
            IF ( Input_Opt%DO_ND50 ) THEN
               CALL DIAG50
     &            ( am_I_Root, Input_Opt, State_Met, State_Chm, RC )

               ! Trap potential errors
               IF ( RC /= GC_SUCCESS ) THEN
                  ErrMsg = 'Error encountered in "Diag50"!'
                  CALL Error_Stop( ErrMsg, ThisLoc )
               ENDIF
            ENDIF

            ! Increment diagnostic timestep counter. (ccc, 5/13/09)
            CALL Set_Ct_Diag( INCREMENT=.TRUE. )

            ! Planeflight diagnostic moved to be after chemistry, kyu
            !### Debug
            IF ( prtDebug ) CALL Debug_Msg( '### MAIN: a DIAGNOSTICS' )
         ENDIF

         !--------------------------------------------------------------
         !   ***** T I M E S E R I E S   D I A G N O S T I C S  *****
         !
         ! NOTE: Since we are saving soluble tracers, we must move
         !       the ND40, ND49, and ND52 timeseries diagnostics
         !       to after the call to DO_WETDEP (bmy, 4/22/04)
         !--------------------------------------------------------------

         ! Station timeseries
         IF ( ITS_TIME_FOR_DIAG48( Input_Opt ) ) THEN
            CALL DIAG48
     &         ( am_I_Root, Input_Opt, State_Met, State_Chm, RC )

            ! Trap potential errors
            IF ( RC /= GC_SUCCESS ) THEN
               ErrMsg = 'Error encountered in "Diag48"!'
               CALL Error_Stop( ErrMsg, ThisLoc )
            ENDIF
         ENDIF
         IF ( prtDebug ) CALL Debug_Msg( '### MAIN: after DIAG48' )

         ! 3-D timeseries
         IF ( ITS_TIME_FOR_DIAG49( Input_Opt ) ) THEN
            CALL DIAG49( am_I_Root, Input_Opt,  State_Met,
     &                   State_Chm, State_Diag, RC         )

            ! Trap potential errors
            IF ( RC /= GC_SUCCESS ) THEN
               ErrMsg = 'Error encountered in "Diag49"!'
               CALL Error_Stop( ErrMsg, ThisLoc )
            ENDIF
         ENDIF
         IF ( prtDebug ) CALL Debug_Msg( '### MAIN: after DIAG49' )

         ! Ship timeseries
         IF ( ITS_TIME_FOR_DIAG63( Input_Opt ) ) THEN
            CALL DIAG63( am_I_Root, Input_Opt, State_Chm, RC )

            ! Trap potential errors
            IF ( RC /= GC_SUCCESS ) THEN
               ErrMsg = 'Error encountered in "Diag63"!'
               CALL Error_Stop( ErrMsg, ThisLoc )
            ENDIF
         ENDIF
         IF ( prtDebug ) CALL Debug_Msg( '### MAIN: after DIAG63' )

         ! Morning or afternoon timeseries
         IF ( Input_Opt%DO_ND51 ) THEN
            CALL DIAG51
     &         ( am_I_Root, Input_Opt, State_Met, State_Chm, RC )

            ! Trap potential errors
            IF ( RC /= GC_SUCCESS ) THEN
               ErrMsg = 'Error encountered in "Diag51"!'
               CALL Error_Stop( ErrMsg, ThisLoc )
            ENDIF
         ENDIF

         IF ( Input_Opt%DO_ND51b ) THEN
            CALL DIAG51b
     &         ( am_I_Root, Input_Opt, State_Met, State_Chm, RC )

            ! Trap potential errors
            IF ( RC /= GC_SUCCESS ) THEN
               ErrMsg = 'Error encountered in "Diag51b"!'
               CALL Error_Stop( ErrMsg, ThisLoc )
            ENDIF
         ENDIF
         IF ( prtDebug ) CALL Debug_Msg( '### MAIN: after DIAG51' )

         ! Comment out for now 
         !! Column timeseries
         !IF ( ND52 > 0 .and. ITS_TIME_FOR_ND52() ) THEN
         !   CALL DIAG52
         !   IF ( prtDebug ) CALL Debug_Msg( '### MAIN: a ND52' )
         !ENDIF

         !### After diagnostics
         IF ( prtDebug ) CALL Debug_Msg( '### MAIN: after TIMESERIES' )

#if defined( USE_TIMERS )
         CALL GEOS_Timer_End( "All diagnostics",              RC )
         CALL GEOS_Timer_End( "Output",                       RC )
         CALL GEOS_Timer_End( "=> Binary punch diagnostics",  RC )
#endif

#endif

         !--------------------------------------------------------------
         !            ***** H I S T O R Y   W R I T E *****
         !--------------------------------------------------------------
#if defined( USE_TIMERS )
         CALL GEOS_Timer_Start( "All diagnostics",           RC )
         CALL GEOS_Timer_Start( "Output",                    RC )
         CALL GEOS_Timer_Start( "=> History (netCDF diags)", RC )
#endif

         ! Write HISTORY ITEMS in each diagnostic collection to disk
         ! (or skip writing if it is not the proper output time.
         CALL History_Write( am_I_Root, State_Chm%Spc_Units, RC )

         ! Trap potential errors
         IF ( RC /= GC_SUCCESS ) THEN
            ErrMsg = 'Error encountered in "History_Write"!'
            CALL Error_Stop( ErrMsg, ThisLoc )
         ENDIF

#if defined( USE_TIMERS )
         CALL GEOS_Timer_End( "All diagnostics",           RC )
         CALL GEOS_Timer_End( "Output",                    RC )
         CALL GEOS_Timer_End( "=> History (netCDF diags)", RC )
#endif

         !==============================================================
         !  ***** E N D   O F   D Y N A M I C   T I M E S T E P *****
         !==============================================================
      ENDDO

      !=================================================================
      !     ***** C O P Y   I - 3   F I E L D S *****
      !
      !     The I-3 fields at the end of this outer timestep 
      !     become the fields at the beginning of the next timestep
      !=================================================================
      CALL Copy_I3_Fields( State_Met )
      IF ( prtDebug ) THEN 
         CALL Debug_Msg( '### MAIN: after COPY_I3_FIELDS' )
      ENDIF

      ENDDO

      !=================================================================
      !         ***** C L E A N U P   A N D   Q U I T *****
      !=================================================================
 9999 CONTINUE

#if defined( USE_TIMERS )
      CALL GEOS_Timer_End ( "Timesteps", RC )
#endif
      !-----------------------------------------------------------------
      !         ***** O B S P A C K   D I A G N O S T I C S *****
      !
      ! Flush any unwritten ObsPack data to disk and finalize
      !-----------------------------------------------------------------
      IF ( Input_Opt%Do_ObsPack ) THEN

#if defined( USE_TIMERS )
         CALL GEOS_Timer_Start( "All diagnostics",        RC )
         CALL GEOS_Timer_Start( "=> ObsPack diagnostics", RC )
#endif

         ! Write any remaining ObsPack data to disk, and immediately
         ! thereafter free the ObsPack pointer fields of State_Diag
         IF ( ASSOCIATED( State_Diag%ObsPack_id ) ) THEN
            CALL ObsPack_Write_Output( am_I_Root,  Input_Opt,
     &                                 State_Diag, RC )
            IF ( RC /= GC_SUCCESS ) THEN
               ErrMsg = 'Error encountered in "ObsPack_Write_Output"!'
               CALL Error_Stop( ErrMsg, ThisLoc )
            ENDIF
         ENDIF

         ! Free the ObsPack species mapping fields of State_Diag
         CALL ObsPack_SpeciesMap_Cleanup( am_I_Root, Input_Opt,
     &                                    State_Diag, RC         )
         IF ( RC /= GC_SUCCESS ) THEN
            ErrMsg = 
     &         'Error encountered in "ObsPack_SpeciesMap_Cleanup"!'
            CALL Error_Stop( ErrMsg, ThisLoc )
         ENDIF

#if defined( USE_TIMERS )
         CALL GEOS_Timer_End( "All diagnostics",        RC ) 
         CALL GEOS_Timer_End( "=> ObsPack diagnostics", RC )
#endif      
      ENDIF

      !-----------------------------------------------------------------
      ! Print the mass-weighted mean OH concentration (if applicable)
      !-----------------------------------------------------------------
      CALL Print_Diag_OH( am_I_Root, Input_Opt, RC )

      ! Trap potential errors
      IF ( RC /= GC_SUCCESS ) THEN
         ErrMsg = 'Error encountered in "Print_Diag_OH"!'
         CALL Error_Stop( ErrMsg, ThisLoc )
      ENDIF

#if defined( BPCH_DIAG )
      !-----------------------------------------------------------------
      ! Write final conditions for benchmark simulations
      !-----------------------------------------------------------------
#if defined( USE_TIMERS )
      CALL GEOS_Timer_Start( "All diagnostics",             RC )
      CALL GEOS_Timer_Start( "=> Binary punch diagnostics", RC )
      CALL GEOS_Timer_Start( "Output",                      RC )
#endif 

      ! Save final tracer mass [kg] to disk for model benchmarking
      ! for Rn-Pb-Be and full-chem benchmark simulations
      ! NOTE: Species concentrations enter the subroutine in [kg/kg dry air]
      ! and are converted locally to [kg] (mps, 11/29/16)
      IF ( LSTDRUN ) THEN 
         CALL StdRun( am_I_Root, Input_Opt,  State_Met,
     &                State_Chm, State_Diag, RC,        Lbegin=.FALSE. )

         ! Trap potential errors
         IF ( RC /= GC_SUCCESS ) THEN
            ErrMsg = 'Error encountered in "StdRun"!'
            CALL Error_Stop( ErrMsg, ThisLoc )
         ENDIF
      ENDIF

#if defined( USE_TIMERS )
      CALL GEOS_Timer_End( "All diagnostics",             RC )
      CALL GEOS_Timer_End( "=> Binary punch diagnostics", RC )
      CALL GEOS_Timer_End( "Output",                      RC )
#endif 

#endif

      !-----------------------------------------------------------------
      ! Finalize GEOS-Chem
      !-----------------------------------------------------------------
#if defined( USE_TIMERS )
      CALL GEOS_Timer_Start( "Finalization", RC )
#endif 

      ! Close all files
      CALL CLOSE_FILES()
      IF ( prtDebug ) CALL Debug_Msg( '### MAIN: a CLOSE_FILES' )   
 
      !%%% NOTE: Call HISTORY_CLEANUP from cleanup.F.  This will
      !%%% close all netCDF files upon both normal or abnormal exits.

      ! Deallocate fields of the Input Options object
      CALL Cleanup_Input_Opt( am_I_Root, Input_Opt, RC )
      IF ( RC /= GC_SUCCESS ) THEN
         ErrMsg = 'Error encountered in "Cleanup_Input_Opt"!'
         CALL Error_Stop( ErrMsg, ThisLoc )
      ENDIF
      IF ( prtDebug ) CALL Debug_Msg( '### MAIN: a cleanup Input_Opt' )

      ! Deallocate fields of the Chemistry State object
      CALL Cleanup_State_Chm( am_I_Root, State_Chm, RC )
      IF ( RC /= GC_SUCCESS ) THEN
         ErrMsg = 'Error encountered in "Cleanup_State_Chm"!'
         CALL Error_Stop( ErrMsg, ThisLoc )
      ENDIF
      IF ( prtDebug ) CALL Debug_Msg( '### MAIN: a cleanup State_Chm' )

      ! Deallocate fields of the Diagnostics State object
      CALL Cleanup_State_Diag( am_I_Root, State_Diag, RC )
      IF ( RC /= GC_SUCCESS ) THEN
         ErrMsg = 'Error encountered in "Cleanup_State_Diag"!'
         CALL Error_Stop( ErrMsg, ThisLoc )
      ENDIF
      IF ( prtDebug ) CALL Debug_Msg( '### MAIN: a cleanup State_Diag' )

      ! Deallocate fields of the Meteorology State object
      CALL Cleanup_State_Met( am_I_Root, State_Met, RC )
      IF ( RC /= GC_SUCCESS ) THEN
         ErrMsg = 'Error encountered in "Cleanup_State_Met"!'
         CALL Error_Stop( ErrMsg, ThisLoc )
      ENDIF
      IF ( prtDebug ) CALL Debug_Msg( '### MAIN: a cleanup State_Met' )

      ! Deallocate dynamic module arrays
      CALL CleanUp( am_I_Root, .FALSE., RC )
      IF ( RC /= GC_SUCCESS ) THEN
         ErrMsg = 'Error encountered in "Cleanup"!'
         CALL Error_Stop( ErrMsg, ThisLoc )
      ENDIF

      ! Deallocate the derived type object w/ mapping info
      CALL Cleanup_Mapping( mapping )

      ! Deallocate fields of the diagnostics list object
      CALL Cleanup_DiagList( am_I_Root, Diag_List, RC )
      IF ( RC /= GC_SUCCESS ) THEN
         ErrMsg = 'Error encountered in "Cleanup_DiagList"!'
         CALL Error_Stop( ErrMsg, ThisLoc )
      ENDIF
      IF ( prtDebug ) CALL Debug_Msg( '### MAIN: a cleanup Diag_List' )

!-----------------------------------------------------------------------------
! TWO-WAY NESTING OPTION
! This is only invoked when compiling GEOS-Chem with COUPLE=y
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! %%%%% NOTE: THIS OPTION WAS PROBABLY BROKEN WITH THE ADDITION OF HEMCO %%%%%
! %%%%% AND FLEXCHEM INTO GEOS-CHEM V10-01 AND v11-01.  BUYER BEWARE.    %%%%%
! %%%%%  --  Bob Yantosca (22 Jan 2018)                                  %%%%%
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
#if defined( EXCHANGE )
      ! Finalize the two-way nesting
      CALL Cleanup_Exchange( am_I_Root )
#endif
!-----------------------------------------------------------------------------

      ! Free the shadow variables in error_mod.F
      CALL Cleanup_Error()

#if defined( GTMM_Hg )
      ! Deallocate arrays from GTMM model for mercury simulation
      IF ( LGTMM ) CALL CleanupCASAarrays()
#endif

      IF ( prtDebug ) CALL Debug_Msg( '### MAIN: a CLEANUP' )

#if defined( USE_TIMERS )
      ! Stop remaining timers
      CALL GEOS_Timer_End( "Finalization", RC )
      CALL GEOS_Timer_End( "GEOS-Chem",    RC )

      ! Print timer output
      CALL GEOS_Timer_PrintAll( am_I_Root, RC )
#endif

      ! Print ending time of simulation
      CALL Display_End_Time()

      ! Flush the buffer to get output
      CALL Flush( 6 )

      CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: display_grid_and_model
!
! !DESCRIPTION: Internal Subroutine DISPLAY\_GRID\_AND\_MODEL displays the 
!  appropriate messages for the given model grid and machine type.  It also 
!  prints the starting time and date (local time) of the GEOS-Chem simulation.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE Display_Grid_and_Model()
! 
! !REVISION HISTORY: 
!  02 Dec 2003 - R. Yantosca - Initial version
!  13 Aug 2010 - R. Yantosca - Added ProTeX headers
!  13 Aug 2010 - R. Yantosca - Added extra output
!  02 Feb 2012 - R. Yantosca - Added output for GEOS-5.7.x met fields
!  28 Feb 2012 - R. Yantosca - Removed support for GEOS-3
!  19 Mar 2012 - R. Yantosca - Now echo info for 0.25 x 0.3125 runs
!  19 Mar 2012 - R. Yantosca - Now echo info if ISORROPIA is turned off
!  22 Oct 2012 - R. Yantosca - Now echo info if -DDEVEL is used
!  01 Mar 2017 - R. Yantosca - Now echo info if netCDF can do compression
!  16 Mar 2017 - M. Sulprizio- Print model version to log file
!  26 Jul 2018 - R. Yantosca - Now get version # from gc_classic_version.H
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! For system time stamp
      CHARACTER(LEN=16) :: STAMP

#include "gc_classic_version.H"

      !-----------------------------------------------------------------
      ! Print resolution info
      !-----------------------------------------------------------------

#if   defined( GRID4x5   )
      WRITE( 6, '(a)' )                   
     &    REPEAT( '*', 13 )                                          //
     &    '   S T A R T I N G   4 x 5   G E O S - C H E M   '        //
     &    REPEAT( '*', 13 )

#elif defined( GRID2x25  )
      WRITE( 6, '(a)' ) 
     &    REPEAT( '*', 13 )                                          // 
     &    '   S T A R T I N G   2 x 2.5   G E O S - C H E M   '      //
     &    REPEAT( '*', 13 )

#elif defined( GRID05x0625 )
      WRITE( 6, '(a)' ) 
     &    REPEAT( '*', 13 )                                          // 
     &    '   S T A R T I N G   0.5 x 0.625   G E O S - C H E M   '  //
     &    REPEAT( '*', 13 )

#elif defined( GRID025x03125 )
      WRITE( 6, '(a)' ) 
     &    REPEAT( '*', 13 )                                          // 
     &    '   S T A R T I N G  0.25 x 0.3125  G E O S - C H E M   '  //
     &    REPEAT( '*', 13 )

#endif


      !-----------------------------------------------------------------
      ! Mode of simulation
      !-----------------------------------------------------------------
      WRITE( 6, 100 ) 'GEOS-Chem "Classic"'

      !-----------------------------------------------------------------
      ! Print model version
      !-----------------------------------------------------------------
      WRITE( 6, 110 ) TRIM( GC_CLASSIC_VERSION )

      !-----------------------------------------------------------------
      ! Print compiler
      !-----------------------------------------------------------------
#if defined( LINUX_PGI )
      WRITE( 6, 120 ) 'PGI compiler (aka pgfortran)'
#elif defined( LINUX_IFORT )
      WRITE( 6, 120  ) 'Intel Fortran Compiler (aka ifort)'
#elif defined( LINUX_GFORTRAN )
      WRITE( 6, 120 ) 'GNU Fortran compiler (aka gfortran)'
#endif

      !-----------------------------------------------------------------
      ! Print met field info
      !-----------------------------------------------------------------
#if defined( GEOS_FP )

#if defined( GRIDREDUCED ) 
      WRITE( 6, 130 ) 'GMAO GEOS-FP (on reduced 47-layer vertical grid)'
#else
      WRITE( 6, 130 ) 'GMAO GEOS-FP (on native 72-layer vertical grid)'
#endif

#elif defined( MERRA2 )

#if defined( GRIDREDUCED )
      WRITE( 6, 130 ) 'GMAO MERRA-2 (on reduced 47-layer vertical grid)'
#else
      WRITE( 6, 130 ) 'GMAO MERRA-2 (on native 72-layer vertical grid)'
#endif

#endif

      !-----------------------------------------------------------------
      ! Print status of ISORROPIA
      !-----------------------------------------------------------------
#if defined( NO_ISORROPIA ) 
      WRITE( 6, 140 ) 'TURNED OFF FOR TESTING'
#else
      WRITE( 6, 140 ) 'ON'
#endif

      !-----------------------------------------------------------------
      ! Print status of OpenMP
      !-----------------------------------------------------------------
#if defined( NO_OMP ) 
      WRITE( 6, 150 ) 'OFF'
#else
      WRITE( 6, 150 ) 'ON'
#endif

      !-----------------------------------------------------------------
      ! Print status of binary punch (bpch) diagnostics
      !-----------------------------------------------------------------
#if defined( BPCH_DIAG ) 
      WRITE( 6, 160 ) 'ON'
#else
      WRITE( 6, 160 ) 'OFF'
#endif

      !-----------------------------------------------------------------
      ! Print status of netCDF diagnostics (aka History) - always on
      !-----------------------------------------------------------------
      WRITE( 6, 170 ) 'ON'

      !-----------------------------------------------------------------
      ! Print msg if netCDF compression is supported
      !-----------------------------------------------------------------
#if defined( NC_HAS_COMPRESSION )
      WRITE( 6, 180 ) 'SUPPORTED'
#else
      WRITE( 6, 180 ) 
     & 'NOT SUPPORTED (or shut off w/ NC_NODEFLATE=y)'
#endif

      !-----------------------------------------------------------------
      ! System time stamp
      !-----------------------------------------------------------------
      STAMP = SYSTEM_TIMESTAMP()
      WRITE( 6, 190 ) STAMP

      !-----------------------------------------------------------------
      ! Format strings
      !-----------------------------------------------------------------
 100  FORMAT( /, '===> Mode of operation         : ', a             )
 110  FORMAT(    '===> GEOS-Chem version         : ', a             )
 120  FORMAT(    '===> Compiler                  : ', a             )
 130  FORMAT(    '===> Driven by meteorology     : ', a             )
 140  FORMAT(    '===> ISORROPIA ATE package     : ', a             )
 150  FORMAT(    '===> Parallelization w/ OpenMP : ', a             )
 160  FORMAT(    '===> Binary punch diagnostics  : ', a             )
 170  FORMAT(    '===> netCDF diagnostics        : ', a             )
 180  FORMAT(    '===> netCDF file compression   : ', a             )
 190  FORMAT( /, '===> SIMULATION START TIME: ',      a, ' <===', / )

      END SUBROUTINE DISPLAY_GRID_AND_MODEL
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ctm_flush
!
! !DESCRIPTION: Internal subroutine CTM\_FLUSH flushes certain diagnostic
! file buffers to disk. 
!\\
!\\
! CTM\_FLUSH should normally be called after each diagnostic output, so that 
! in case the run dies, the output files from the last diagnostic timestep 
! will not be lost.  
!\\
!\\
! FLUSH is an intrinsic FORTRAN subroutine and takes as input the unit number 
! of the file to be flushed to disk.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CTM_Flush()
! 
! !REVISION HISTORY: 
!  31 Aug 2000 - R. Yantosca - Initial version
!  13 Aug 2010 - R. Yantosca - Added ProTeX headers
!  06 Aug 2012 - R. Yantosca - IU_BPCH is only global file LUN still needed
!EOP
!------------------------------------------------------------------------------
!BOC

#if defined( BPCH_DIAG )
      ! NOTE: Block out with an #ifdef, so that we don't interfere
      ! with file units in GCHP, or when using GEOS-Chem as an
      ! instance in other models like CESM or ECMWF (bmy, 1/12/18)
      CALL FLUSH( IU_BPCH )  
#endif

      END SUBROUTINE CTM_Flush
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: display_end_time
!
! !DESCRIPTION: Internal subroutine DISPLAY\_END\_TIME prints the ending 
!  time of the GEOS-Chem simulation.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE Display_End_Time()
! 
! !REVISION HISTORY: 
!  03 May 2005 - R. Yantosca - Initial version
!  13 Aug 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      CHARACTER(LEN=16) :: STAMP

      ! Print system time stamp
      STAMP = SYSTEM_TIMESTAMP()
      WRITE( 6, 100 ) STAMP
 100  FORMAT( /, '===> SIMULATION END TIME: ', a, ' <===', / )

      ! Echo info
      WRITE ( 6, 3000 ) 
 3000 FORMAT
     &   ( /, '**************   E N D   O F   G E O S -- C H E M   ',
     &        '**************' )

      END SUBROUTINE Display_End_Time
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_overhead_o3_for_fastj
!
! !DESCRIPTION: Internal subroutine GET\_OVERHEAD\_O3\_FOR\_FASTJ 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE Get_Overhead_O3_For_FastJ( am_I_Root )
!
! INPUT ARGUMENTS:
! 
      LOGICAL, INTENT(IN   )   :: am_I_Root   ! Root CPU?
! 
! !REMARKS:
!  This routine makes use of variables declared in above in the main program
!  (which are visible in all sub-programs below the CONTAINS statement).
!                                                                             .
!  The original code was done in FAST-J routine "set_prof.F", but has been
!  split off to facilitate development of the grid-independent model.
!
! !REVISION HISTORY: 
!  07 Mar 2012 - R. Yantosca - Initial version
!  14 Nov 2013 - R. Yantosca - For GEOS-FP, read O3 from met field files
!  13 Dec 2013 - M. Sulprizio- Moved USE_O3_FROM_MET to the Input_Opt object and
!                              set in initialization stage of GEOS_CHEM
!  22 Oct 2014 - C. Keller   - Added am_I_Root
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! FAST-J is only used for fullchem, offline aerosol
      ! and CH3I simulations, so skip otherwise
      IF ( ITS_A_FULLCHEM_SIM  .or.
     &     ITS_AN_AEROSOL_SIM  .or.
     &     ITS_A_CH3I_SIM     ) THEN

         ! Only execute this if we are doing chemistry
         ! and if it we are at a chemistry timestep
         IF ( LCHEM .and. ITS_TIME_FOR_CHEM() ) THEN
             
            ! Get the overhead O3 column for FAST-J.  Take either the
            ! TOMS O3 data or the column O3 directly from the met fields
            CALL Compute_Overhead_O3( am_I_Root,  
     &                                DAY,
     &                                Input_Opt%USE_O3_FROM_MET, 
     &                                State_Met%TO3 )
         ENDIF
      ENDIF

      END SUBROUTINE Get_Overhead_O3_For_FastJ
!EOC
      END PROGRAM GEOS_Chem
#endif
