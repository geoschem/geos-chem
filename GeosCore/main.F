!BOC
#if !defined(ESMF_)
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
      PROGRAM GEOS_CHEM
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
      USE REGRID_A2A_MOD        ! For horizontal regridding
      USE RESTART_MOD           ! For restart file I/O
      USE State_Chm_Mod         ! Derived type for Chemistry State object
      USE State_Met_Mod         ! Derived type for Meteorology State object
      USE TIME_MOD              ! For computing date & time 
      USE UnitConv_Mod          ! For gas concentration unit conversion

      !-----------------------------------------------------------------
      ! GEOS-Chem chemistry modules
      !-----------------------------------------------------------------
      USE CARBON_MOD            ! For SOA simulation
      USE CHEMISTRY_MOD         ! Driver routines for chemistry
      USE GLOBAL_CH4_MOD        ! For offline CH4 simulation
      USE MERCURY_MOD           ! For offline Hg simulation (driver)
      USE OCEAN_MERCURY_MOD     ! For offline Hg simulation (ocean model)
      USE STRAT_CHEM_MOD        ! For linearized stratospheric chemistry
      USE TOMS_MOD              ! For overhead O3 columns (for FAST-J)
      USE UCX_MOD               ! For unified trop-strat chemistry (SDE)
      USE UVALBEDO_MOD          ! For reading UV albedoes (for FAST-J)

      !-----------------------------------------------------------------
      ! GEOS-Chem deposition modules
      !-----------------------------------------------------------------
      USE DEPO_MERCURY_MOD      ! Deposition for offline Hg simulation
      USE DRYDEP_MOD            ! For dry deposition
      USE WETSCAV_MOD           ! For wet deposition (rainout & washout)

      !-----------------------------------------------------------------
      ! GEOS-Chem diagnostics modules
      !-----------------------------------------------------------------
      USE BENCHMARK_MOD         ! For the 1-month benchmark simulations
      USE CMN_DIAG_MOD          ! Logical switches for G-C diagnostics
      USE DIAG_MOD              ! G-C diagnostic arrays & counters
      USE DIAG41_MOD            ! For ND41  (afternoon PBL       ) diag
      USE DIAG42_MOD            ! For ND42  (SOA products        ) diag
      USE DIAG48_MOD            ! For ND48  (station timeseries  ) diag
      USE DIAG49_MOD            ! For ND49  (inst. timeseries    ) diag
      USE DIAG50_MOD            ! For ND50  (24h avg timeseries  ) diag
      USE DIAG51_MOD            ! For ND51  (satellite timeseries) diag
      USE DIAG51b_MOD           ! For ND51b (satellite timeseries) diag
      USE DIAG63_MOD            ! For ND63  (PARANOX timeseries  ) diag
      USE DIAG_OH_MOD           ! For ND43  (OH,HO2,etc. prod    ) diag
      USE DIAGNOSTICS_MOD       ! NetCDF diagnostics module (C. Keller)
      USE PLANEFLIGHT_MOD       ! For ND40  (plane flight track  ) diag

      !-----------------------------------------------------------------
      ! GEOS-Chem dynamics modules
      !-----------------------------------------------------------------
      USE CHEMGRID_MOD          ! For the dynamic tropopause
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
      USE GCAP_READ_MOD         ! For reading GCAP met data
      USE GEOSFP_READ_MOD       ! For reading GEOS-FP data
      USE MERRA2_READ_MOD       ! For reading MERRA2 data
      USE MERRA_A1_MOD          ! For reading MERRA A1 data
      USE MERRA_A3_MOD          ! For reading MERRA A3 data
      USE MERRA_CN_MOD          ! For reading MERRA CN data
      USE MERRA_I6_MOD          ! For reading MERRA I6 data
      USE A3_READ_MOD           ! For reading A3 data (all other met)
      USE A6_READ_MOD           ! For reading A6 data (all other met)
      USE I6_READ_MOD           ! For reading I6 data (all other met)
#if defined( EXCHANGE )
      USE EXCHANGE_MOD          ! For two-way coupled simulations
#endif

#if defined( RRTMG )
      !-----------------------------------------------------------------
      ! Radiation modules
      !-----------------------------------------------------------------
      USE RRTMG_RAD_TRANSFER_MOD, ONLY: DO_RRTMG_RAD_TRANSFER,
     &                                  LW_UFLUX, LW_DFLUX,
     &                                  SW_UFLUX, SW_DFLUX,
     &                                  LW_UFLUXC, LW_DFLUXC,
     &                                  SW_UFLUXC, SW_DFLUXC,
     &                                  INIT_SURFACE_RAD,
     &                                  READ_SURFACE_RAD,
     &                                  INIT_STRAT_CLIM,
     &                                  READ_STRAT_CLIM,
     &                                  INIT_MCICA_CLOUDS,
     &                                  SET_SPECMASK
      USE CMN_FJX_MOD,            ONLY: NSPECRADMENU,
     &                                  LSPECRADMENU
      USE rrtmg_lw_init,  ONLY : rrtmg_lw_ini
      USE rrtmg_sw_init,  ONLY : rrtmg_sw_ini

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
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! Scalars
      LOGICAL                  :: FIRST_RT = .TRUE.
      LOGICAL                  :: wasModisRead
      LOGICAL                  :: prtDebug
      INTEGER                  :: I,             IOS,         J
      INTEGER                  :: K,             L,           N
      INTEGER                  :: JDAY,          NDIAGTIME,   N_DYN
      INTEGER                  :: NNN,           N_DYN_STEPS, NSECb 
      INTEGER                  :: N_STEP,        YEAR,        MONTH
      INTEGER                  :: DAY,           DAY_OF_YEAR, SEASON
      INTEGER                  :: NYMD,          NYMDb,       NHMS
      INTEGER                  :: ELAPSED_SEC,   NHMSb,       RC
      INTEGER                  :: ELAPSED_TODAY, HOUR,        MINUTE
      INTEGER                  :: laiYYYY,       laiMM,       id_CLOCK
      INTEGER                  :: id_H2O
      REAL(f8)                 :: TAU,           TAUb         
      CHARACTER(LEN=255)       :: ZTYPE

      ! Set cloud flag for RRTMG
      INTEGER                  :: ICLD=0,ISEED=10
      INTEGER                  :: RR

      ! Arrays
      INTEGER                  :: DATE(2)

      ! Derived type objects
      TYPE(MapWeight), POINTER :: mapping(:,:)

      !-----------------------------------------------------------------
      ! %%%% REPLICATING GCHP FUNCTIONALITY IN EXISTING GEOS-CHEM %%%%
      !
      ! We must declare the Chemistry State (State_Chm), Meteorology
      ! State (State_Met), and Input Options (Input_Opt) objects
      ! for use with the GEOS-Chem High Performance code. 
      !
      TYPE(OptInput)           :: Input_Opt   ! Input Options object
      TYPE(ChmState)           :: State_Chm   ! Chemistry State object
      TYPE(MetState)           :: State_Met   ! Meteorology State object
      !
      ! Also define NI, NJ, NL for more consistent naming w/r/t the 
      ! ESMF interface.  Initialize these from IIPAR, JJPAR, JJPAR
      ! (bmy, 11/9/12)
      !
      INTEGER                  :: NI          ! # of longitudes 
      INTEGER                  :: NJ          ! # of latitudes
      INTEGER                  :: NL          ! # of levels
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
      !
      LOGICAL, PARAMETER       :: am_I_Root = .TRUE. 
      !
      ! Define shadow variables for fields in Input_Opt (bmy, 6/25/14)
      LOGICAL                  :: ITS_A_FULLCHEM_SIM 
      LOGICAL                  :: ITS_A_CH4_SIM
      LOGICAL                  :: ITS_A_MERCURY_SIM  
      LOGICAL                  :: ITS_A_TAGCO_SIM
      LOGICAL                  :: ITS_AN_AEROSOL_SIM 
      LOGICAL                  :: ITS_A_CH3I_SIM
      LOGICAL                  :: ITS_A_SPECIALTY_SIM
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
      LOGICAL                  :: LUNZIP
      LOGICAL                  :: LVARTROP
      LOGICAL                  :: LWAIT
      LOGICAL                  :: LWETD
      LOGICAL                  :: LWINDO
      LOGICAL                  :: USE_OLSON_2001

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

      !-----------------------------------------------------------------
      ! %%%% REPLICATING GCHP FUNCTIONALITY IN EXISTING GEOS-CHEM %%%%
      !
      ! Assume a successful return until otherwise
      !
      RC                      = GC_SUCCESS
      !
      ! Define NI, NJ, NL dimensions for use below
      !
      NI                      = IIPAR
      NJ                      = JJPAR
      NL                      = LLPAR
      !
      ! Define dimensions for fields of the Input Options object
      ! NOTE: At this point we have not yet read in the "input.geos"
      ! file so these are just placeholder values (bmy, 11/7/12)
      !
      Input_Opt%MAX_DIAG      = MAX_DIAG
      Input_Opt%MAX_TRCS      = NNPAR
      Input_Opt%MAX_MEMB      = 15
      Input_Opt%MAX_FAMS      = MAXFAM
      Input_Opt%MAX_DEP       = NNPAR
      Input_Opt%LINOZ_NLEVELS = 25
      Input_Opt%LINOZ_NLAT    = 18
      Input_Opt%LINOZ_NMONTHS = 12
      Input_Opt%LINOZ_NFIELDS = 7
      Input_Opt%RootCPU       = .true.
      !
      ! Call the routine GC_Allocate_All (located in module file
      ! GeosCore/gc_environment_mod.F90) to allocate all lat/lon
      ! allocatable arrays used by GEOS-Chem, as well as the Input 
      ! Options object.  The IIPAR and JJPAR dimensions are not 
      ! declared as PARAMETERs, but are module variables that are 
      ! initialized in this call.
      !
      CALL GC_Allocate_All( am_I_root, Input_Opt, RC )
      IF ( RC /= GC_SUCCESS ) GOTO 9999
      !-----------------------------------------------------------------
      
      ! Display current grid resolution and data set type
      CALL DISPLAY_GRID_AND_MODEL

      !=================================================================
      !            ***** I N I T I A L I Z A T I O N *****
      !=================================================================

      ! Initialize pointers
      mapping => NULL()

#if defined( USE_TIMERS )
      ! Call timer initilization
      CALL GEOS_Timer_Setup( 1 )
      CALL GEOS_Timer_Add( "GEOS-Chem",           RC )
      CALL GEOS_Timer_Add( "Initialization",      RC )
      CALL GEOS_Timer_Add( "Timesteps",           RC )
      CALL GEOS_Timer_Add( "HEMCO",               RC )
      CALL GEOS_Timer_Add( "All chemistry",       RC )
      CALL GEOS_Timer_Add( "=> Strat chem",       RC )
      CALL GEOS_Timer_Add( "=> Gas-phase chem",   RC )
      CALL GEOS_Timer_Add( "=> All aerosol chem", RC )
      CALL GEOS_Timer_Add( "Transport",           RC )
      CALL GEOS_Timer_Add( "Convection",          RC )
      CALL GEOS_Timer_Add( "Dry deposition",      RC )
      CALL GEOS_Timer_Add( "Wet deposition",      RC )
      CALL GEOS_Timer_Add( "RRTMG",               RC )
      CALL GEOS_Timer_Add( "Diagnostics",         RC )
      CALL GEOS_Timer_Add( "Reading met fields",  RC )
      CALL GEOS_Timer_Add( "Reading restart file",RC )
      CALL GEOS_Timer_Add( "Writing restart file",RC )
      CALL GEOS_Timer_Start( "GEOS-Chem",         RC )
      CALL GEOS_Timer_Start( "Initialization",    RC )
#endif

      ! Read input file and call init routines from other modules
      CALL Read_Input_File( am_I_Root, Input_Opt, RC ) 

      ! Initialize the regridding module by storing shadow copies
      CALL Initialize_Regridding( am_I_Root, Input_Opt, RC )

      ! Store shadow copies of am_I_Root, Input_Opt in error_mod.F
      CALL Init_Error( am_I_Root, Input_Opt, RC )

      ! Copy values from Input_Opt.  These replace the variables
      ! from logical_mod.F and tracer_mod.F. (bmy, 3/29/13)
      USE_OLSON_2001      =  Input_Opt%USE_OLSON_2001
      ITS_A_FULLCHEM_SIM  =  Input_Opt%ITS_A_FULLCHEM_SIM
      ITS_A_CH4_SIM       =  Input_Opt%ITS_A_CH4_SIM
      ITS_A_MERCURY_SIM   =  Input_Opt%ITS_A_MERCURY_SIM
      ITS_A_TAGCO_SIM     =  Input_Opt%ITS_A_TAGCO_SIM 
      ITS_AN_AEROSOL_SIM  =  Input_Opt%ITS_AN_AEROSOL_SIM
      ITS_A_CH3I_SIM      =  Input_Opt%ITS_A_CH3I_SIM
      ITS_A_SPECIALTY_SIM =  Input_Opt%ITS_A_SPECIALTY_SIM
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
      LUNZIP              =  Input_Opt%LUNZIP
      LVARTROP            =  Input_Opt%LVARTROP
      LWAIT               =  Input_Opt%LWAIT
      LWETD               =  Input_Opt%LWETD
      LWINDO              =  Input_Opt%LWINDO
      USE_OLSON_2001      =  Input_Opt%USE_OLSON_2001

      ! Set a flag to denote if we should print ND70 debug output
      prtDebug            = ( LPRT .and. am_I_Root )
      
      ! Debug outputg
      IF ( prtDebug ) CALL DEBUG_MSG( '### MAIN: a READ_INPUT_FILE' )

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
      !
      ! Initialize Input_Opt, State_Chem, State_Met, objects. 
      ! (bmy, 3/4/13)
      !
      CALL GC_Init_All
     &   ( am_I_Root, Input_Opt, State_Chm, State_Met, RC )
      !
      ! Also allocate arrays in GEOS-Chem modules that have not yet
      ! been initialized (i.e. SEASALT_MOD, CARBON_MOD,  DUST_MOD, 
      ! SULFATE_MOD).  This removes the init calls from the run stage,
      ! which is required when connecting to the GEOS-5 GCM w/ ESMF.
      ! (bmy, 3/4/13)
      CALL GC_Init_Extra( am_I_Root, Input_Opt, State_Chm, RC )
      !
      ! Copy surface area info from AREA_M2 array of grid_mod.F90
      ! to the State_Met%AREA_M2 field. (bmy, 8/22/14)
      ! NOTE: Most routines now use either mixing ratio or kg/m2
      ! and are independent of area. Area will eventually be 
      ! completely removed for grid-independence (ewl, 9/10/15)
      State_Met%AREA_M2 =  AREA_M2

      ! Define advected species ID flags for use below
      id_CLOCK = Ind_('CLOCK','A')
      id_H2O   = Ind_('H2O'  ,'A')
      !-----------------------------------------------------------------

#if   defined( RRTMG )
      ! Initialize surface rad fields
      ! from rrtmg_rad_transfer_mod.F
      !WRITE(6,*) 'Calling INIT_SURFACE_RAD'
      CALL INIT_SURFACE_RAD()
      !WRITE(6,*) 'Calling INIT_STRAT_CLIM'
      CALL INIT_STRAT_CLIM()
      CALL READ_STRAT_CLIM( Input_Opt )
      CALL INIT_MCICA_CLOUDS()

      ! Initialization
       CALL RRTMG_LW_INI()
       CALL RRTMG_SW_INI()

#endif

      ! Initialize diagnostic arrays and counters
      CALL INITIALIZE( am_I_Root, Input_Opt, 2, RC )
      CALL INITIALIZE( am_I_Root, Input_Opt, 3, RC )
      IF ( prtDebug ) CALL DEBUG_MSG( '### MAIN: a INITIALIZE' )

      ! Initialize the new hybrid pressure module.  Define Ap and Bp.
      CALL INIT_PRESSURE( am_I_Root )
      IF ( prtDebug ) CALL DEBUG_MSG( '### MAIN: a INIT_PRESSURE' )

      ! Read annual mean tropopause if not a variable tropopause
      ! read_tropopause is obsolete with variable tropopause
      IF ( .not. Input_Opt%LVARTROP ) THEN
         CALL READ_TROPOPAUSE( am_I_Root, Input_Opt, RC )
         IF ( prtDebug ) CALL DEBUG_MSG( '### MAIN: a READ_TROPOPAUSE' )
      ENDIF

      ! Added to read input file for linoz strat (dbj, jliu, bmy, 10/16/09)
      IF ( LLINOZ ) THEN
         CALL LINOZ_READ( am_I_Root, Input_Opt, RC )
      ENDIF

      ! Define time variables for use below
      NHMS  = GET_NHMS()
      NHMSb = GET_NHMSb()
      NYMD  = GET_NYMD()
      NYMDb = GET_NYMDb()
      TAU   = GET_TAU()
      TAUb  = GET_TAUb()

      !=================================================================
      !        ***** I N I T I A L I Z A T I O N  continued *****
      !=================================================================
#if defined( USE_TIMERS )
      CALL GEOS_Timer_Start( "Reading met fields",  RC )
#endif

      ! Read the initial met fields from disk
      CALL READ_INITIAL_MET_FIELDS()

#if defined( USE_TIMERS )
      CALL GEOS_Timer_End( "Reading met fields",  RC )
#endif

      ! Set dry surface pressure (PS1_DRY) from State_Met%PS1_WET
      ! and compute avg surface pressures near polar caps
      CALL SET_DRY_SURFACE_PRESSURE( State_Met, 1 )
      CALL AVGPOLE( State_Met%PS1_DRY )

      ! Compute avg surface pressure near polar caps
      CALL AVGPOLE( State_Met%PS1_WET ) 
      IF ( prtDebug ) CALL DEBUG_MSG( '### MAIN: a AVGPOLE' )

      ! Initialize surface pressures prior to interpolation
      ! to allow initialization of floating pressures
      State_Met%PSC2_WET = State_Met%PS1_WET
      State_Met%PSC2_DRY = State_Met%PS1_DRY
      CALL SET_FLOATING_PRESSURES( am_I_Root, State_Met, RC )
      IF ( prtDebug ) CALL DEBUG_MSG( '### MAIN: a SET_FLT_PRS' )

      ! Call AIRQNT to compute initial air mass quantities
      ! Do not update initial tracer concentrations since not read 
      ! from restart file yet (ewl, 10/28/15)
      CALL AIRQNT( am_I_Root, Input_opt, State_Met, State_Chm, RC, 
     &             update_mixing_ratio=.FALSE. )
      IF ( prtDebug ) CALL DEBUG_MSG( '### MAIN: a AIRQNT' )

      ! Initialize the derived type object containing
      ! mapping information for the MODIS LAI routines
      IF ( USE_OLSON_2001 ) THEN
         CALL Init_Mapping( am_I_Root, Input_Opt,
     &                      1440, 720, IIPAR, JJPAR, mapping, RC )
      ELSE
         CALL Init_Mapping( am_I_Root, Input_Opt,
     &                       720, 360, IIPAR, JJPAR, mapping, RC )
      ENDIF

      ! Compute the Olson land types that occur in each grid box
      ! (i.e. this is a replacement for rdland.F and vegtype.global)
      CALL Init_Olson_Landmap   ( am_I_Root, Input_Opt, RC        )
      CALL Compute_Olson_Landmap( am_I_Root, mapping,   State_Met )
      CALL Cleanup_Olson_Landmap( am_I_Root                       )

       ! Initialize PBL quantities but do not do mixing
       ! Add option for non-local PBL (Lin, 03/31/09) 
      CALL INIT_MIXING ( am_I_Root, Input_Opt, 
     &                   State_Met, State_Chm, RC ) 
 
      ! Initialize chemistry
      ! Moved here (from chemistry_mod.F and chemdr.F) because some
      ! of the variables are used for non-local PBL mixing BEFORE 
      ! the first call of the chemistry routines (ckeller, 05/19/14).
      IF ( ITS_A_FULLCHEM_SIM .OR. ITS_AN_AEROSOL_SIM ) THEN
         CALL INIT_CHEMISTRY ( am_I_Root, Input_Opt, State_Chm, RC )
      ENDIF

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

      CALL EMISSIONS_INIT( am_I_Root, Input_Opt,
     &                     State_Met, State_Chm, RC )

      ! Stop the run if the HEMCO init sequence dies w/ an error
      IF ( RC /= GC_SUCCESS ) THEN
         CALL ERROR_STOP( 'Error in EMISSIONS_INIT!', 'main.F' )
      ENDIF

#if defined( USE_TIMERS )
      CALL GEOS_Timer_End ( "HEMCO", RC )
#endif

#if defined( NC_DIAG )
      ! Initialize diagnostics other than the default HEMCO diagnostics
      ! initialized above (ckeller, 1/13/15).
      CALL Diagnostics_Init( am_I_Root, Input_Opt, 
     &                       State_Met, State_Chm, RC )
#endif

      !=================================================================
      !       *****  I N I T I A L   C O N D I T I O N S *****
      !=================================================================

      ! Read initial species conditions
      ! NOTE: Species concentrations are read in as [v/v dry air] 
      ! and converted to [kg/kg dry air] (ewl, 8/4/16)
#if defined( USE_TIMERS )
      CALL GEOS_Timer_Start( "Reading restart file",  RC )
#endif

      CALL READ_GC_RESTART( am_I_Root, NYMDb,     NHMSb, 
     &                      Input_Opt, State_Met, State_Chm, RC )
      IF ( prtDebug ) CALL DEBUG_MSG( '### MAIN: a READ_GC_RESTART' )

#if defined( USE_TIMERS )
      CALL GEOS_Timer_End( "Reading restart file",  RC )
#endif

      ! Initialize stratospheric routines
      IF ( LUCX ) THEN

         ! Initialize the UCX module
         CALL INIT_UCX( am_I_Root, Input_Opt, State_Chm )
         IF ( prtDebug ) CALL DEBUG_MSG( '### MAIN: a INIT_UCX' )

         ! Set simple initial tracer conditions [kg/kg dry air]
         CALL SET_INITIAL_MIXRATIOS( am_I_Root, Input_Opt, State_Met,
     &                               State_Chm )
         IF ( prtDebug ) THEN
            CALL DEBUG_MSG( '### MAIN: a SET_INITIAL_MIXRATIOS' )
         ENDIF

      ENDIF

      ! Capture initial state of atmosphere for STE flux calc (ltm, 06/10/12)
      ! NOTE: Species concentrations enter the subroutine in [kg/kg dry air]
      ! and are converted locally to [kg] for chemistry (ewl, 9/18/15)
      IF ( LSCHEM ) THEN
         CALL INIT_STRAT_CHEM( am_I_Root, Input_Opt, State_Chm, 
     &                         State_Met, RC )
      ENDIF

      ! Save initial tracer mass [kg] to disk for model benchmarking
      ! for Rn-Pb-Be and full-chem benchmark simulations
      ! NOTE: Species concentrations enter the subroutine in [kg/kg dry air]
      ! and are converted locally to [kg] (mps, 11/29/16)
      IF ( LSTDRUN ) THEN 
         CALL STDRUN( am_I_Root, Input_Opt, State_Met, State_Chm, RC, 
     &                LBEGIN=.TRUE. )
      ENDIF
 
!-----------------------------------------------------------------------------
! TWO-WAY NESTING OPTION
! This is only invoked when compiling GEOS-Chem with COUPLE=y
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

      !=================================================================
      !      ***** INITIALIZE GEOS-CHEM DIAGNOSTICS AND RESTARTS ***** 
      !=================================================================      

#if defined( NC_DIAG )
      ! Initialize diagnostics other than the default HEMCO diagnostics
      ! initialized above (ckeller, 1/13/15).
      CALL Diagnostics_Init( am_I_Root, Input_Opt, 
     &                       State_Met, State_Chm, RC )
#endif 

      CALL INIT_GC_RESTART( am_I_Root, Input_Opt, State_Chm, RC )

#if defined( USE_TIMERS )
      CALL GEOS_Timer_End( "Initialization", RC )
#endif

      !=================================================================
      !      ***** 6 - H O U R   T I M E S T E P   L O O P  *****
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
#if   defined( GEOS_FP ) || defined( MERRA2 )
      N_DYN_STEPS = 180 / GET_TS_DYN()     ! GEOS-FP/MERRA2 have a 3hr interval
#else
      N_DYN_STEPS = 360 / GET_TS_DYN()     ! All other met has a 6hr interval
#endif

      ! Start a new outer loop
      DO 

      ! Compute time parameters at start of 6-h loop
      CALL SET_CURRENT_TIME

      ! NSECb is # of seconds (measured from 00 GMT today) 
      ! at the start of this 6-hr timestepping loop.
      ! NOTE: Assume we start at the head of each minute (i.e. SECONDS=0)
      HOUR   = GET_HOUR()
      HOUR   = ( HOUR / 6 ) * 6
      MINUTE = GET_MINUTE()
      NSECb  = ( HOUR * 3600 ) + ( MINUTE * 60 )

      ! Get dynamic timestep in seconds
      N_DYN  = 60e+0_fpp * GET_TS_DYN()

      !=================================================================
      !     ***** D Y N A M I C   T I M E S T E P   L O O P *****
      !=================================================================
      DO N_STEP = 1, N_DYN_STEPS
    
#if defined( LINUX_PGI )
         ! Flush the buffer each timestep so we get output
         CALL FLUSH( 6 )
#endif

         ! Compute & print time quantities at start of dyn step
         CALL SET_CURRENT_TIME
         CALL PRINT_CURRENT_TIME

         ! Set time variables for dynamic loop
         DAY_OF_YEAR   = GET_DAY_OF_YEAR()
         DAY           = GET_DAY()
         ELAPSED_SEC   = GET_ELAPSED_SEC()
         MONTH         = GET_MONTH()
         NHMS          = GET_NHMS()
         NYMD          = GET_NYMD()
         HOUR          = GET_HOUR()
         MINUTE        = GET_MINUTE()
         TAU           = GET_TAU()
         YEAR          = GET_YEAR()
         SEASON        = GET_SEASON()
         ELAPSED_TODAY = ( HOUR * 3600 ) + ( MINUTE * 60 )

         !### Debug
         IF ( prtDebug ) THEN
            CALL DEBUG_MSG( '### MAIN: a SET_CURRENT_TIME' )
         ENDIf

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
#endif

         CALL EMISSIONS_RUN( am_I_Root, Input_Opt, State_Met, 
     &                       State_Chm, ITS_TIME_FOR_EMIS(), 1, RC )

         !### Debug
         IF ( prtDebug ) THEN
            CALL DEBUG_MSG( '### MAIN: a HEMCO PHASE 1' )
         ENDIf

         !==============================================================
         !   ***** W R I T E   D I A G N O S T I C   F I L E S *****
         !==============================================================
#if defined( USE_TIMERS )
         CALL GEOS_Timer_Start( "Diagnostics", RC )
#endif

         ! Write HEMCO diagnostics as well as new GEOS-Chem diagnostics
         ! (the latter is currently only defined in development mode)
         ! (ckeller, 4/1/15). 
         CALL Diagnostics_Write( am_I_Root, Input_Opt, State_Chm,
     &                           .FALSE.,   RC)

#if defined( USE_TIMERS )
         CALL GEOS_Timer_End( "HEMCO", RC )
#endif
         ! Original diagnostics
         IF ( ITS_TIME_FOR_BPCH() ) THEN

            ! Set time at end of diagnostic timestep
            CALL SET_DIAGe( TAU )

            ! Write bpch file
            IF ( DO_DIAG_WRITE ) THEN

               CALL DIAG3
     &              ( am_I_Root, Input_Opt, State_Met, State_Chm, RC )
 
               ! Flush file units
               CALL CTM_FLUSH

            ENDIF

            !===========================================================
            !    *****  W R I T E   R E S T A R T   F I L E S  *****
            !===========================================================
            ! Determine and print STE fluxes (ltm, 06/10/12)
            IF ( LSCHEM ) THEN
               CALL CALC_STE( am_I_Root, Input_Opt, State_Chm, 
     &                        State_Met, RC )
            ENDIF

#if defined( USE_TIMERS )
            CALL GEOS_Timer_Start( "Writing restart file",  RC )
#endif

            ! Make atmospheric restart file
            CALL WRITE_GC_RESTART( am_I_Root, NYMD, NHMS, Input_Opt,
     &                             State_Chm, State_Met, RC )

            !### Debug
            IF ( prtDebug ) THEN
               CALL DEBUG_MSG( '### MAIN: a WRITE_GC_RESTART' )
            ENDIF

#if defined( USE_TIMERS )
            CALL GEOS_Timer_End( "Writing restart file",  RC )
#endif

#if defined( USE_TIMERS )
            CALL GEOS_Timer_Start( "HEMCO", RC )
#endif

            ! Force the output of a HEMCO restart file. 
            ! (ckeller, 4/1/15)
            CALL Diagnostics_Write( am_I_Root, Input_Opt, State_Chm,
     &                              .TRUE.,    RC)

#if defined( USE_TIMERS )
            CALL GEOS_Timer_End( "HEMCO", RC )
#endif

            ! Set time at beginning of next diagnostic timestep
            CALL SET_DIAGb( TAU )

            !===========================================================
            !        ***** Z E R O   D I A G N O S T I C S *****
            !===========================================================
            CALL INITIALIZE( am_I_Root, Input_Opt, 2, RC )  ! Zero diag arrays

            CALL INITIALIZE( am_I_Root, Input_Opt, 3, RC )  ! Zero counters

         ENDIF

#if defined( USE_TIMERS )
         CALL GEOS_Timer_End( "Diagnostics", RC )
#endif

         !=============================================================
         !   ***** W R I T E   MERCURY RESTART  F I L E *****
         !     ***** MUST be done after call to diag3 *****
         !=============================================================
         ! Make land restart file: for GTMM runs only, beginning of each 
         ! month but not start of the run.
         IF ( LGTMM .AND. ITS_A_NEW_MONTH() .AND. NYMD /= NYMDb ) THEN
            IF (.NOT.( ITS_TIME_FOR_BPCH() )) THEN

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
         !          ****** R E A D   M E T   F I E L D S ******
         !
         !    For clarity, we have split off the code that reads met 
         !    fields into the internal subroutine READ_MET_FIELDS.
         !    Lightning NOx emissions are also computed at the same
         !    time when A3 or A6 met fields are read from disk.
         !==============================================================
#if defined( USE_TIMERS )
         CALL GEOS_Timer_Start( "Reading met fields",  RC )
#endif

         CALL READ_MET_FIELDS()

#if defined( USE_TIMERS )
         CALL GEOS_Timer_End( "Reading met fields",  RC )
#endif

#if   defined( GEOS_5 )

         ! When TOMS data is not available, then we will use
         ! the O3 columns directly from the met fields.  
         ! Currently, this is only done for GEOS-5 met.
         IF ( YEAR > LAST_TOMS_YEAR ) THEN
            Input_Opt%USE_O3_FROM_MET = .TRUE.
         ELSE
            Input_Opt%USE_O3_FROM_MET = .FALSE.
         ENDIF

#elif defined( GEOS_FP ) || defined( MERRA2 )

         ! For GEOS-FP met, always use the O3 columns
         ! from the A1 met field files (bmy, 11/14/13)
         Input_Opt%USE_O3_FROM_MET = .TRUE.

#else

         ! When TOMS data is not available, then we will use
         ! the TOMS data from the last valid year.  This is
         ! done for all met field products other than GEOS-5.
         Input_Opt%USE_O3_FROM_MET = .FALSE.

#endif

         !==============================================================
         !              ***** D A I L Y   D A T A *****
         !==============================================================
         laiMM = -1
         wasModisRead = .FALSE.
         IF ( ITS_A_NEW_DAY() ) THEN 

            !-----------------------------------------------------------
            !    %%%%%%% UPDATED LEAF AREA INDEX ALGORITHM %%%%%%%
            ! 
            ! Note that GEOS-Chem uses LAI data from two separate
            ! sources.  The dry deposition and soil NOx modules rely
            ! on the data from lai*.global ASCII files.  These files
            ! (which are pre-processed offline by IDL codes) are 
            ! generated for each specific GEOS-Chem grid configuration
            ! (e.g. 4x5, 2x25, 0.5x0.666 nested grids).  These files
            ! are read from disk by routine RDLAI, which saves
            ! the LAI data into the XLAI and XYLAI arrays.  XLAI and
            ! XYLAI store the leaf area index as a function of Olson 
            ! land type (cf Olson 1992 land map).
            !
            ! However, the MEGAN biogenic emissions code relies on 
            ! LAI data stored at 1x1 resolution stored in bpch format.
            ! These binary files are read by routine RDISOLAI (and
            ! other underlying routines in lai_mod.F), and are
            ! regridded on-the-fly to the current GEOS-Chem grid
            ! resolution.
            !
            ! Therefore, these two sources of LAI data present an
            ! inconsistency that should be resolved.  Also, for the 
            ! Grid-Indpendent GEOS-Chem project, we must move away
            ! from ASCII files (which prevent interfacing with 
            ! external GCMs).  We also cannot assume any particular
            ! horizontal grid, since that is now to be specified at
            ! the start of the simulation.
            !
            ! Also, to facilitate simulations at ultra-fine horizontal
            ! resolution, we will eventually adopt the Olson 2001 land 
            ! map, which has a native resolution of 0.25 x 0.25 
            ! degrees, and likewise use an updated version of the 
            ! MODIS LAI data at 0.25 x 0.25 resolution.
            !
            ! To resolve these issues, we have created a new module
            ! (modis_lai_mod.F90) which reads from the MODIS LAI data
            ! in netCDF format at the native resolution and then
            ! regrids the LAI data to GEOS-Chem resolution on-the-fly.
            ! The XLAI array is are populated for backwards
            ! compatibility with the existing legacy codes.  The LAI
            ! arrays used for MEGAN (ISOLAI, PMISOLAI, MISOLAI, and
            ! NMISOLAI) are now replaced by arrays GC_LAI, GC_LAI_PM, 
            ! GC_LAI_CM, and GC_LAI_NM) from modis_lai_mod.F.
            !
            ! We have validated that the new scheme generates identical 
            ! XLAI arrays w/r/t the old scheme.  The arrays GC_LAI etc. 
            ! differ from the ISOLAI etc. arrays slightly (but generally 
            ! agree to within 0.001).  This is due to the fact that the 
            ! ISOLAI arrays were regridded from 1 x 1 native resolution, 
            ! but now we are regridding from much finer resolution 
            ! (either 0.5 x 0.5 or 0.25 x 0.25).
            !
            ! NOTES:
            ! (1) At the present time, we have removed all references
            !     to the obsolete XYLAI array and its parent module
            !     CMN_VEL_mod.F.
            ! (2) As of Dec 2012, the XLAI and XLAI2 arrays are now
            !     carried in State_Met.  We can remove the obsolete
            !     module Headers/CMN_DEP_mod.F.
            ! 
            !      -- Bob Yantosca (09 Apr 2012)
            !-----------------------------------------------------------

            ! Find the proper month and year for LAI  based on the current 
            ! month & year.  LAI months begin in the middle of the month.
            CALL Find_Lai_Month( day_of_year, month, year, 
     &                                        laiMM, laiYYYY )
            
            ! Read MODIS LAI from disk (if it is a new LAI month)
            ! Logical flag MODISREAD determines if data was just read in
            CALL Read_Modis_Lai( am_I_Root, Input_Opt,    laiYYYY, 
     &                           laiMM,     wasModisRead, RC       )

            ! Interpolate LAI quantities from native grid to GEOS-Chem grid
            CALL Compute_Modis_Lai( am_I_Root,    Input_Opt, State_Met, 
     &                              day_of_year,  laiMM,     mapping,
     &                              wasModisRead, RC )

            !### Debug
            IF ( prtDebug ) CALL DEBUG_MSG ( '### MAIN: a DAILY DATA' )
         ENDIF

         !==============================================================
         !   ***** I N T E R P O L A T E   Q U A N T I T I E S *****   
         !==============================================================

         ! Interpolate I-6 fields to current dynamic timestep, 
         ! based on their values at NSEC and NSEC+N_DYN
         CALL INTERP( NSECb, ELAPSED_TODAY,
     &                N_DYN, Input_Opt,     State_Met )

         ! Case of variable tropopause:
         ! Check LLTROP and set LMIN, LMAX, and LPAUSE
         ! since this is not done with READ_TROPOPAUSE anymore.
         ! (Need to double-check that LMIN, Lmax are not used before-phs) 
         IF ( LVARTROP ) CALL CHECK_VAR_TROP( State_Met )

         ! If we are not doing transport, then make sure that
         ! the floating pressure is set to PSC2_WET (bdf, bmy, 8/22/02)
         ! Now also includes PSC2_DRY (ewl, 5/4/16)
         IF ( .not. LTRAN ) THEN
            CALL SET_FLOATING_PRESSURES( am_I_Root, State_Met, RC )
         ENDIF

         ! Compute updated airmass quantities at each grid box 
         ! and update tracer concentration to conserve tracer mass
         ! (ewl, 10/28/15)
         CALL AIRQNT( am_I_Root, Input_opt, State_Met, State_Chm, RC, 
     &                update_mixing_ratio=.TRUE. )

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
            ! Only force strat once if using UCX
            IF (LSETH2O) LSETH2O = .FALSE.
         ENDIF

         ! Compute the cosine of the solar zenith angle array
         ! State_Met%SUNCOS     = at the current time
         ! State_Met%SUNCOSmid  = at the midpt of the chem timestep
         ! State_Met%SUNCOSmid5 = at the midpt of the chem timestep 5hrs ago
         CALL GET_COSINE_SZA( am_I_Root, Input_Opt, State_Met, RC )

         ! Compute tropopause height for ND55 diagnostic
         ! Consider updating all Met diagnostics here (ewl, 1/22/16) 
#if defined( BPCH_DIAG )
         IF ( Input_Opt%ND55 > 0 ) THEN
            CALL DIAG_TROPOPAUSE( am_I_root, Input_Opt, State_Met, RC )
         ENDIF
#elif defined( NC_DIAG )
         ! If netcdf diagnostics, output if Met diagnostic group is on
         IF ( Input_Opt%ND68 > 0 ) THEN
            CALL DIAG_TROPOPAUSE( am_I_root, Input_Opt, State_Met, RC )
         ENDIF
#endif

         ! Update dynamic timestep
         CALL SET_CT_DYN( INCREMENT=.TRUE. )

         !### Debug
         IF ( prtDebug ) CALL DEBUG_MSG( '### MAIN: a INTERP, etc' )

         ! Get averaging intervals for local-time diagnostics
         ! (NOTE: maybe improve this later on)
         ! Placed after interpolation to get correct value of TROPP. 
         ! (ccc, 12/9/08)
         CALL DIAG_2PM( State_Met )

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
               CALL GET_UVALBEDO(am_I_Root, Input_Opt, State_Met, RC)

               ! Get TOMS overhead O3 columns for photolysis from
               ! the HEMCO data structure (bmy, 3/20/15)
               CALL READ_TOMS( am_I_Root, Input_Opt, RC )

            ENDIF

            ! Read data required for Hg2 gas-particle partitioning 
            ! (H Amos, 25 Oct 2011)
            IF ( ITS_A_MERCURY_SIM ) THEN 
               CALL READ_Hg2_PARTITIONING
     &            ( am_I_Root, Input_Opt, State_Met, MONTH, RC )
               IF ( prtDebug ) THEN
                 CALL DEBUG_MSG( '### MAIN: a READ_HG2_PARTITIONING')
               ENDIF
            ENDIF
         ENDIF

         !==============================================================
         !              ***** T R A N S P O R T *****
         !==============================================================
         IF ( ITS_TIME_FOR_DYN() ) THEN

#if defined( USE_TIMERS )
            CALL GEOS_Timer_Start( "Transport", RC )
#endif

            ! Output BC's
            ! Save boundary conditions (global grid) for future nested run
            IF ( LWINDO ) THEN
               CALL SAVE_GLOBAL_TPCORE_BC
     &            ( am_I_Root, Input_Opt, State_Chm, RC )
            ENDIF

!-----------------------------------------------------------------------------
! TWO-WAY NESTING OPTION
! This is only invoked when compiling GEOS-Chem with COUPLE=y
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
            ! NOTE: Tracer concentration units are converted locally
            ! within transport to [kg/kg total air] for advection 
            ! (ewl, 9/18/15)
            IF ( LTRAN ) THEN 
               CALL DO_TRANSPORT
     &            ( am_I_Root, Input_Opt, State_Met, State_Chm, RC )
               IF ( RC /= GC_SUCCESS ) THEN
                  CALL ERROR_STOP('Error in DO_TRANSPORT!', 'main.F')
               ENDIF

               !### Debug
               IF ( prtDebug ) THEN 
                  CALL DEBUG_MSG( '### MAIN: a DO_TRANSPORT' )
               ENDIF
            ENDIF

            ! Initialize wet scavenging and wetdep fields after
            ! the airmass quantities are reset after transport
#if defined ( TOMAS )
            ! ... TOMAS microphysics: Always call SETUP_WETSCAV ...
            CALL SETUP_WETSCAV
     &           ( am_I_Root, Input_Opt, State_Met, State_Chm, RC )

#else 

            ! ... Standard GEOS-Chem: Call INIT_WETSCAV if   ...
            ! ... convection or wet scavenging or chemistry are turned on ...
            IF ( LCONV .or. LWETD .or. LCHEM ) THEN
               CALL SETUP_WETSCAV
     &              ( am_I_Root, Input_Opt, State_Met, State_Chm, RC )

               !### Debug
               IF ( prtDebug ) THEN 
                  CALL DEBUG_MSG( '### MAIN: a SETUP_WETSCAV' )
               ENDIF
            ENDIF

#endif

#if defined( USE_TIMERS )
            CALL GEOS_Timer_End( "Transport", RC )
#endif

         ENDIF

         ! SDE 10/16/13: Set clock tracer if relevant
         IF ( id_CLOCK > 0 ) THEN
            CALL SET_CLOCK_TRAC( GET_TS_DYN(), State_Chm )
         ENDIF

         !--------------------------------------------------------------
         ! Compute PBL height
         ! Move this call from the PBL mixing routines because the PBL
         ! height is used by drydep and some of the emissions routines.
         ! (ckeller, 3/5/15) 
         !--------------------------------------------------------------
         CALL COMPUTE_PBL_HEIGHT ( State_Met )

         !### Debug
         IF ( prtDebug ) THEN 
            CALL DEBUG_MSG( '### MAIN: a COMPUTE_PBL_HEIGHT' )
         ENDIF

         !--------------------------------------------------------------
         ! Test for emission timestep
         ! Now always do emissions here, even for full-mixing
         ! (ckeller, 3/5/15)
         !--------------------------------------------------------------
         IF ( ITS_TIME_FOR_EMIS() ) THEN

            ! Increment emission counter
            CALL SET_CT_EMIS( INCREMENT=.TRUE. )
     
            !===========================================================
            !         ***** D R Y   D E P O S I T I O N *****
            !===========================================================
            IF ( LDRYD ) THEN

#if defined( USE_TIMERS )
               CALL GEOS_Timer_Start( "Dry deposition", RC )
#endif

               CALL DO_DRYDEP( am_I_Root, Input_Opt, 
     &                         State_Met, State_Chm, RC )
               IF ( RC /= GC_SUCCESS ) THEN
                  CALL ERROR_STOP('Error in DO_DRYDEP!', 'main.F')
               ENDIF

#if defined( USE_TIMERS )
               CALL GEOS_Timer_End ( "Dry deposition", RC )
#endif

               !### Debug
               IF ( prtDebug ) THEN 
                  CALL DEBUG_MSG( '### MAIN: a DO_DRYDEP' )
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

            ! EMISSIONS_RUN will call HEMCO run phase 2. HEMCO run phase
            ! only calculates emissions. All data has been read to disk
            ! in phase 1 at the beginning of the time step. 
            ! (ckeller, 4/1/15)
            CALL EMISSIONS_RUN( am_I_Root, Input_Opt, State_Met, 
     &                          State_Chm, ITS_TIME_FOR_EMIS(), 2, RC )

            ! Stop the run if the HEMCO run sequence dies w/ an error
            IF ( RC /= GC_SUCCESS ) THEN
               CALL ERROR_STOP( 'Error in EMISSIONS_RUN after drydep!',
     &                          'main.F' )
            ENDIF

            !### Debug
            IF ( prtDebug ) THEN 
               CALL DEBUG_MSG( '### MAIN: a HEMCO PHASE 2' )
            ENDIF

#if defined( USE_TIMERS )
            CALL GEOS_Timer_End( "HEMCO", RC )
#endif
         ENDIF

         !-------------------------------
         ! Test for convection timestep
         !-------------------------------
         IF ( ITS_TIME_FOR_CONV() ) THEN
            
#if defined( USE_TIMERS )
            CALL GEOS_Timer_Start( "Convection", RC )
#endif

            ! Increment the convection timestep
            CALL SET_CT_CONV( INCREMENT=.TRUE. )

            !===========================================================
            !      ***** M I X E D   L A Y E R   M I X I N G *****
            !===========================================================

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
            CALL DO_MIXING ( am_I_Root, Input_Opt, 
     &                       State_Met, State_Chm, RC )
            IF ( RC /= GC_SUCCESS ) THEN
               CALL ERROR_STOP('Error in DO_MIXING!', 'main.F')
            ENDIF

            IF ( RC /= GC_SUCCESS ) THEN
               CALL ERROR_STOP( 'Error in MIX_TRACERS!', 'main.F' )
            ENDIF

            IF ( prtDebug ) CALL DEBUG_MSG( '### MAIN: a TURBDAY:2' )

            !===========================================================
            !        ***** C L O U D   C O N V E C T I O N *****
            !===========================================================

            IF ( LCONV ) THEN

               ! Call the appropriate convection routine
               ! NOTE: Tracer concentration units are converted locally
               ! to [kg/kg total air] for convection (ewl, 9/18/15)
               CALL DO_CONVECTION
     &            ( am_I_Root, Input_Opt, State_Met, State_Chm, RC )
               IF ( RC /= GC_SUCCESS ) THEN
                  CALL ERROR_STOP('Error in DO_CONVECTION!', 'main.F')
               ENDIF

               !### Debug
               IF ( prtDebug ) THEN
                  CALL DEBUG_MSG( '### MAIN: a CONVECTION' )
               ENDIF
            ENDIF 

#if defined( USE_TIMERS )
            CALL GEOS_Timer_End( "Convection", RC )
#endif

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
         CALL GET_OVERHEAD_O3_FOR_FASTJ( am_I_Root )
         
         ! Also need to compute avg P, T for CH4 chemistry (bmy, 1/16/01)
         IF ( ITS_A_CH4_SIM ) CALL CH4_AVGTP( State_Met )

         ! Every chemistry timestep...
         IF ( ITS_TIME_FOR_CHEM() ) THEN

            ! Increment chemistry timestep counter
            CALL SET_CT_CHEM( INCREMENT=.TRUE. )

            ! SDE 05/28/13: Set H2O to State_Chm tracer if relevant
            IF ( ITS_A_FULLCHEM_SIM .and. id_H2O > 0 ) THEN
               CALL SET_H2O_TRAC( am_I_Root, (.not. LUCX), Input_Opt ,
     &                            State_Met, State_Chm, RC )
            ENDIF

            ! Do GEOS-Chem chemistry
            ! NOTE: Tracer concentration units are converted locally
            ! to [kg] for all of chemistry. We will replace use of [kg] 
            ! once FlexChem is implemented (ewl, 9/18/15)
            CALL DO_CHEMISTRY( am_I_Root, Input_Opt, 
     &                         State_Chm, State_Met, RC ) 
            IF ( RC /= GC_SUCCESS ) THEN
               CALL ERROR_STOP('Error in DO_CHEMISTRY!', 'main.F')
            ENDIF

         ENDIF 

         ! Check State_Chm tracers (yxw)   
#if defined( GRID05x0666 ) || defined( GRID05x0625 ) || defined( GRID025x03125 )
         CALL CHECK_SPC_NESTED( State_Chm, Input_Opt,
     &                          'after chemistry' )
#endif 

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
            CALL DO_WETDEP
     &         ( am_I_Root, Input_Opt, State_Met, State_Chm, RC )
            IF ( RC /= GC_SUCCESS ) THEN
               CALL ERROR_STOP('Error in DO_WETDEP!', 'main.F')
            ENDIF

#if defined( USE_TIMERS )
            CALL GEOS_Timer_End( "Wet deposition", RC )
#endif

         ENDIF

         !==============================================================
         !      ***** U P D A T E  O P T I C A L  D E P T H *****          
         !==============================================================
         ! Recalculate the optical depth at the wavelength(s) specified
         ! in the Radiation Menu. This must be done before the call to any
         ! diagnostic and only on a chemistry timestep.
         ! (skim, 02/05/11)
         IF ( ITS_TIME_FOR_CHEM() ) THEN
            CALL RECOMPUTE_OD
     &         ( am_I_Root, Input_Opt, State_Met, State_Chm, RC )
         ENDIF

#if defined( RRTMG )

         !==============================================================
         !  ***** R R T M G   R A D I A T I V E   T R A N S F E R *****
         !==============================================================
         ! NOTE: Tracer concentration units are converted locally to
         ! [kg] in RRTMG. Units should eventually be [kg/kg]  (ewl, 9/18/15)
         IF ( Input_opt%LRAD ) THEN 
            IF ( ITS_TIME_FOR_RT() ) THEN

#if defined( USE_TIMERS )
               CALL GEOS_Timer_Start( "RRTMG", RC )
#endif

               ! Splash page
               IF ( am_I_Root .and. FIRST_RT ) THEN
                  WRITE( 6, '(a)' ) REPEAT( '#', 79 )
                  WRITE( 6, 500 )     
     &              'R R T M G : Radiative Transfer Model (by AER)'
 500              FORMAT( '#####', 12x, a, 12x, '#####' )
                  WRITE( 6, '(a)' ) REPEAT( '#', 79 )
                  FIRST_RT = .FALSE.
               ENDIF

               iseed=iseed+15

               !clear-sky is output with all-sky, so we just need to run once 
               !regardless of whether both are required or just one
               IF ( Input_Opt%LSKYRAD(2) ) THEN
                  !Determine if we are doing CLEAR-SKY (1) or ALL-SKY (2)
                  ICLD=1
               ELSE
                  ICLD=0 !clouds are on
               ENDIF

               !set all values to one for baseline RT calc
               CALL SET_SPECMASK(0)

               !First call for baseline (all aerosol and ozone)
               CALL DO_RRTMG_RAD_TRANSFER( am_I_Root, DAY, MONTH, 
     &                                     ICLD, 0, ISEED, Input_Opt, 
     &                                     State_Met, State_Chm, RC)

               !Second call with aerosols or ozone if requested
               DO N = 1, NSPECRADMENU
                  IF (LSPECRADMENU(N).EQ.1) THEN
                     !generate mask for species in RT
                     CALL SET_SPECMASK(N)
                     CALL DO_RRTMG_RAD_TRANSFER( am_I_Root, DAY, MONTH, 
     &                                           ICLD, N, ISEED, 
     &                                           Input_Opt, State_Met, 
     &                                           State_Chm, RC )
                  ENDIF
               ENDDO

               ! Stop the run if DO_RRTMG_RAD_TRANSFER returns an error
               IF ( RC /= GC_SUCCESS ) THEN
                  CALL ERROR_STOP('Error in DO_RRTMG_RAD_TRANSFER!',
     &                            'main.F')
               ENDIF

               ! Increment radiation timestep counter
               CALL SET_CT_RAD( INCREMENT=.TRUE. )

               !### Debug
               IF ( prtDebug ) THEN 
                  CALL DEBUG_MSG( '### MAIN: a DO_RRTMG_RAD_TRANSFER' )
               ENDIF

#if defined( USE_TIMERS )
               CALL GEOS_Timer_End( "RRTMG", RC )
#endif

            ENDIF
         ENDIF

#endif

         ! Eventually write dobson columns into diagnostics array
#if defined( DIAG_DEVEL )
         CALL CalcDobsonColumn ( am_I_Root, Input_Opt,
     &                           State_Met, State_Chm, RC )
#endif

         !==============================================================
         !   ***** I N C R E M E N T   E L A P S E D   T I M E *****
         !============================================================== 
         ! Moved before diagnostics to count the last timestep as done.
         ! Need to save timestamps for filenames.
         ! (ccc, 5/13/09)

#if defined( USE_TIMERS )
         CALL GEOS_Timer_Start( "Diagnostics", RC )
#endif

         ! Plane following diagnostic
         IF ( ND40 > 0 .and. DO_DIAG_WRITE ) THEN 
         
            ! Call SETUP_PLANEFLIGHT routine if necessary
            IF ( ITS_A_NEW_DAY() ) THEN
               
               ! If it's a full-chemistry simulation but LCHEM=F,
               ! or if it's an offline simulation, call setup routine 
               IF ( ITS_A_FULLCHEM_SIM ) THEN
                  IF ( .not. LCHEM ) THEN
                     CALL SETUP_PLANEFLIGHT( am_I_Root, Input_Opt,
     &                                       State_Chm, RC )
                  ENDIF
               ELSE
                  CALL SETUP_PLANEFLIGHT( am_I_Root, Input_Opt, 
     &                                    State_Chm, RC )
               ENDIF
            ENDIF
         ENDIF

         CALL TIMESTAMP_DIAG
         CALL SET_ELAPSED_MIN
         CALL SET_CURRENT_TIME
         IF ( prtDebug ) THEN
            CALL DEBUG_MSG( '### MAIN: after SET_ELAPSED_MIN' )
         ENDIF

         ! Plane following diagnostic, moved here to be consistent
         ! with species array and not be biased by emissions, kyu 30
         ! jan 2015
         IF ( ND40 > 0 ) THEN 
            ! Archive data along the flight track
            CALL PLANEFLIGHT
     &           ( am_I_Root, Input_Opt, State_Met, State_Chm, RC )
         ENDIF

         IF ( prtDebug ) CALL DEBUG_MSG( '### MAIN: after DIAG40' )

         !==============================================================
         !       ***** A R C H I V E   D I A G N O S T I C S *****
         !==============================================================
         IF ( ITS_TIME_FOR_DIAG() ) THEN

            !### Debug
            IF ( prtDebug ) CALL DEBUG_MSG( '### MAIN: b DIAGNOSTICS' )

            ! Accumulate several diagnostic quantities
            CALL DIAG1
     &         ( am_I_Root, Input_Opt, State_Met, State_Chm, RC )
            IF ( prtDebug ) CALL DEBUG_MSG( '### MAIN: after DIAG1' )

            ! ND41: save PBL height in 1200-1600 LT (amf)
            ! (for comparison w/ Holzworth, 1967)
            IF ( ND41 > 0 ) CALL DIAG41
            IF ( prtDebug ) CALL DEBUG_MSG( '### MAIN: after DIAG41' )

            ! ND42: SOA concentrations [ug/m3]
            IF ( ND42 > 0 ) THEN
               CALL DIAG42( Input_Opt, State_Met, State_Chm )
            ENDIF
            IF ( prtDebug ) CALL DEBUG_MSG( '### MAIN: after DIAG42' )

            ! 24-hr timeseries
            IF ( DO_SAVE_DIAG50 ) THEN
               CALL DIAG50
     &            ( am_I_Root, Input_Opt, State_Met, State_Chm, RC )
            ENDIF

            ! Increment diagnostic timestep counter. (ccc, 5/13/09)
            CALL SET_CT_DIAG( INCREMENT=.TRUE. )

            ! Planeflight diagnostic moved to be after chemistry, kyu
            !### Debug
            IF ( prtDebug ) CALL DEBUG_MSG( '### MAIN: a DIAGNOSTICS' )
         ENDIF

         !==============================================================
         !   ***** T I M E S E R I E S   D I A G N O S T I C S  *****
         !
         ! NOTE: Since we are saving soluble tracers, we must move
         !       the ND40, ND49, and ND52 timeseries diagnostics
         !       to after the call to DO_WETDEP (bmy, 4/22/04)
         !============================================================== 

         IF ( prtDebug ) CALL DEBUG_MSG( '### MAIN: before TIMESERIES' )

         ! Station timeseries
         IF ( ITS_TIME_FOR_DIAG48() ) THEN
            CALL DIAG48
     &         ( am_I_Root, Input_Opt, State_Met, State_Chm, RC )
         ENDIF
         IF ( prtDebug ) CALL DEBUG_MSG( '### MAIN: after DIAG48' )

         ! 3-D timeseries
         IF ( ITS_TIME_FOR_DIAG49() ) THEN
            CALL DIAG49
     &         ( am_I_Root, Input_Opt, State_Met, State_Chm, RC )
         ENDIF
         IF ( prtDebug ) CALL DEBUG_MSG( '### MAIN: after DIAG49' )

         ! Ship timeseries
         IF ( ITS_TIME_FOR_DIAG63() ) THEN
            CALL DIAG63( am_I_Root, Input_Opt, State_Chm, RC )
         ENDIF
         IF ( prtDebug ) CALL DEBUG_MSG( '### MAIN: after DIAG63' )

         ! Morning or afternoon timeseries
         IF ( DO_SAVE_DIAG51 ) THEN
            CALL DIAG51
     &         ( am_I_Root, Input_Opt, State_Met, State_Chm, RC )
         ENDIF
         IF ( DO_SAVE_DIAG51b ) THEN
            CALL DIAG51b
     &         ( am_I_Root, Input_Opt, State_Met, State_Chm, RC )
         ENDIF
         IF ( prtDebug ) CALL DEBUG_MSG( '### MAIN: after DIAG51' )

         ! Comment out for now 
         !! Column timeseries
         !IF ( ND52 > 0 .and. ITS_TIME_FOR_ND52() ) THEN
         !   CALL DIAG52
         !   IF ( prtDebug ) CALL DEBUG_MSG( '### MAIN: a ND52' )
         !ENDIF

!#if defined( NC_DIAG )
!         ! Write diagnostic data to the GEOS-Chem diagnostics netCDF file.
!         ! Each diagnostic has a specified averaging period specified 
!         ! in Input_Opt. (ckeller, bmy, 1/15/15)
!         CALL Diagnostics_Write( am_I_Root, Input_Opt, State_Chm,
!     &                           .FALSE.,    RC )
!#endif

         !### After diagnostics
         IF ( prtDebug ) CALL DEBUG_MSG( '### MAIN: after TIMESERIES' )

         !==============================================================
         !  ***** E N D   O F   D Y N A M I C   T I M E S T E P *****
         !==============================================================

         ! Check for NaN, Negatives, Infinities in tracers each time diag are
         ! saved. (ccc, 5/13/09)
         IF ( ITS_TIME_FOR_DIAG() ) THEN

         ! Sometimes tracers in the stratosphere can be negative at 
         ! the nested-grid domain edges. Force them to be zero before
         ! CHECK_SPC. (yxw)
#if defined( GRID05x0666 ) || defined( GRID05x0625 ) || defined( GRID025x03125)
            CALL CHECK_SPC_NESTED( State_Chm, Input_Opt,
     &                             'End of Dynamic Loop' )
#else
            CALL CHECK_SPC( State_Chm, Input_Opt, 
     &                      'End of Dynamic Loop' )
#endif

         ENDIF

#if defined( USE_TIMERS )
         CALL GEOS_Timer_End( "Diagnostics", RC )
#endif

      ENDDO

      !=================================================================
      !     ***** C O P Y   I - 3   O R   I - 6   F I E L D S *****
      !
      !     The I-3 or I-6 fields at the end of this outer timestep 
      !     become the fields at the beginning of the next timestep
      !=================================================================
      CALL COPY_I3_I6_FIELDS( State_Met )
      IF ( prtDebug ) THEN 
         CALL DEBUG_MSG( '### MAIN: after COPY_I3_I6_FIELDS' )
      ENDIF

      ENDDO

      !=================================================================
      !         ***** C L E A N U P   A N D   Q U I T *****
      !=================================================================
 9999 CONTINUE

      ! Remove all files from temporary directory 
      IF ( LUNZIP ) THEN
         
         ! Type of operation
         ZTYPE = 'remove all'

         ! Remove A3, A6, I6 fields
         CALL UNZIP_A3_FIELDS( Input_Opt, ZTYPE )
         CALL UNZIP_A6_FIELDS( Input_Opt, ZTYPE )
         CALL UNZIP_I6_FIELDS( Input_Opt, ZTYPE )

#if   defined( GCAP )
         ! Remove GCAP PHIS field (if necessary)
         CALL UNZIP_GCAP_FIELDS( Input_Opt, ZTYPE )
#endif

      ENDIF

      ! Print the mass-weighted mean OH concentration (if applicable)
      CALL PRINT_DIAG_OH( am_I_Root, Input_Opt, RC )

!#if defined( NC_DIAG )
!      ! Write diagnostic data to the GEOS-Chem restart netCDF file.
!      ! Each diagnostic has a specified averaging period specified 
!      ! in Input_Opt. (ckeller, bmy, 1/15/15)
!      CALL Diagnostics_Write( am_I_Root, Input_Opt, State_Chm,
!     &                        .FALSE.,   RC )
!#endif

      ! Save final tracer mass [kg] to disk for model benchmarking
      ! for Rn-Pb-Be and full-chem benchmark simulations
      ! NOTE: Species concentrations enter the subroutine in [kg/kg dry air]
      ! and are converted locally to [kg] (mps, 11/29/16)
      IF ( LSTDRUN ) THEN 
         CALL STDRUN( am_I_Root, Input_Opt, State_Met, State_Chm, RC,
     &                LBEGIN=.FALSE. )
      ENDIF

#if defined( USE_TIMERS )
      CALL GEOS_Timer_End ( "Timesteps", RC )
#endif

      ! Close all files
      CALL CLOSE_FILES
      IF ( prtDebug ) CALL DEBUG_MSG( '### MAIN: a CLOSE_FILES' )

      ! Deallocate fields of the Input Options object
      CALL Cleanup_Input_Opt( am_I_Root, Input_Opt, RC )
      IF ( prtDebug ) CALL DEBUG_MSG( '### MAIN: a cleanup Input_Opt' )

      ! Deallocate fields of the Chemistry State object
      CALL Cleanup_State_Chm( am_I_Root, State_Chm, RC )
      IF ( prtDebug ) CALL DEBUG_MSG( '### MAIN: a cleanup State_Chm' )

      ! Deallocate fields of the Meteorology State object
      CALL Cleanup_State_Met( am_I_Root, State_Met, RC )
      IF ( prtDebug ) CALL DEBUG_MSG( '### MAIN: a cleanup State_met' )

      ! Deallocate dynamic module arrays
      CALL CLEANUP( am_I_Root, .FALSE., RC )

      ! Deallocate the derived type object w/ mapping info
      CALL CLEANUP_MAPPING( mapping )

!-----------------------------------------------------------------------------
! TWO-WAY NESTING OPTION
! This is only invoked when compiling GEOS-Chem with COUPLE=y
!
#if defined( EXCHANGE )
      ! Finalize the two-way nesting
      CALL CLEANUP_EXCHANGE( am_I_Root )
#endif
!-----------------------------------------------------------------------------

      ! Free the shadow variables in error_mod.F
      CALL CLEANUP_ERROR()
#if defined( GTMM_Hg )
      ! Deallocate arrays from GTMM model for mercury simulation
      IF ( LGTMM ) CALL CleanupCASAarrays
#endif

      IF ( prtDebug ) CALL DEBUG_MSG( '### MAIN: a CLEANUP' )

#if defined( USE_TIMERS )
      CALL GEOS_Timer_End( "GEOS-Chem", RC )
      CALL GEOS_Timer_PrintAll( am_I_Root, RC )
#endif

      ! Print ending time of simulation
      CALL DISPLAY_END_TIME

      ! Flush the buffer to get output
      CALL FLUSH( 6 )

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
      SUBROUTINE DISPLAY_GRID_AND_MODEL
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
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! For system time stamp
      CHARACTER(LEN=16) :: STAMP

      !-------------------------------------------------
      ! Print resolution info
      !-------------------------------------------------
#if   defined( GRID4x5   )
      WRITE( 6, '(a)' )                   
     &    REPEAT( '*', 13 )                                          //
     &    '   S T A R T I N G   4 x 5   G E O S--C H E M   '         //
     &    REPEAT( '*', 13 )

#elif defined( GRID2x25  )
      WRITE( 6, '(a)' ) 
     &    REPEAT( '*', 13 )                                          // 
     &    '   S T A R T I N G   2 x 2.5   G E O S--C H E M   '       //
     &    REPEAT( '*', 13 )

#elif defined( GRID1x125 )
      WRITE( 6, '(a)' ) 
     &    REPEAT( '*', 13 )                                          // 
     &    '   S T A R T I N G   1 x 1.25   G E O S--C H E M   '      //
     &    REPEAT( '*', 13 )

#elif defined( GRID1x1 )
      WRITE( 6, '(a)' ) 
     &    REPEAT( '*', 13 )                                          // 
     &    '   S T A R T I N G   1 x 1   G E O S -- C H E M   '       //
     &    REPEAT( '*', 13 )

#elif defined( GRID05x0666 )
      WRITE( 6, '(a)' ) 
     &    REPEAT( '*', 13 )                                          // 
     &    '   S T A R T I N G   0.5 x 0.666   G E O S -- C H E M   ' //
     &    REPEAT( '*', 13 )

#elif defined( GRID05x0625 )
      WRITE( 6, '(a)' ) 
     &    REPEAT( '*', 13 )                                          // 
     &    '   S T A R T I N G   0.5 x 0.625   G E O S -- C H E M   ' //
     &    REPEAT( '*', 13 )

#elif defined( GRID025x03125 )
      WRITE( 6, '(a)' ) 
     &    REPEAT( '*', 13 )                                          // 
     &    '   S T A R T I N G  0.25 x 0.3125  G E O S -- C H E M   ' //
     &    REPEAT( '*', 13 )

#endif

      !-------------------------------------------------
      ! Print machine info
      !-------------------------------------------------

      ! Get the proper FORMAT statement for the model being used
#if   defined( COMPAQ    )
      WRITE( 6, '(a)' ) 'Created w/ HP/COMPAQ Alpha compiler'
#elif defined( IBM_AIX   )
      WRITE( 6, '(a)' ) 'Created w/ IBM-AIX compiler'
#elif defined( LINUX_PGI )
      WRITE( 6, '(a)' ) 'Created w/ LINUX/PGI compiler'
#elif defined( LINUX_IFORT )
      WRITE( 6, '(a)' ) 'Created w/ LINUX/IFORT compiler'
      WRITE( 6, '(a)' ) 'Use ifort -V to print version information'
#elif defined( SGI_MIPS  )
      WRITE( 6, '(a)' ) 'Created w/ SGI MIPSpro compiler'
#elif defined( SPARC     )
      WRITE( 6, '(a)' ) 'Created w/ Sun/SPARC compiler'
#endif

      !-------------------------------------------------
      ! Print met field info
      !-------------------------------------------------
#if   defined( GEOS_4 )
      WRITE( 6, '(a)' ) 'Using GMAO GEOS-4 met fields'
#elif defined( GEOS_5 )
      WRITE( 6, '(a)' ) 'Using GMAO GEOS-5 met fields'
#elif defined( GEOS_FP )
      WRITE( 6, '(a)' ) 'Using GMAO GEOS-FP met fields'
#elif defined( MERRA )
      WRITE( 6, '(a)' ) 'Using GMAO MERRA met fields'
#elif defined( MERRA2 )
      WRITE( 6, '(a)' ) 'Using GMAO MERRA2 met fields'
#elif defined( GCAP  )
      WRITE( 6, '(a)' ) 'Using GCAP/GISS met fields'
#endif

      !-------------------------------------------------
      ! Print msg if ISORROPIA is turned off for debug
      !-------------------------------------------------
#if   defined( NO_ISORROPIA ) 
      WRITE( 6, '(a)' ) REPEAT( '#', 39 )
      WRITE( 6, '(a)' ) '# ISORROPIA TURNED OFF FOR TESTING!!! #'
      WRITE( 6, '(a)' ) REPEAT( '#', 39 )
#endif

      !-------------------------------------------------
      ! Print msg if OpenMP is turned off for debug
      !-------------------------------------------------
#if   defined( NO_OMP ) 
      WRITE( 6, '(a)' ) REPEAT( '#', 39 )
      WRITE( 6, '(a)' ) '# OPENMP TURNED OFF FOR TESTING!!!    #'
      WRITE( 6, '(a)' ) REPEAT( '#', 39 )
#endif

      !-------------------------------------------------
      ! Print msg if DEVEL is turned on for testing
      !-------------------------------------------------
#if   defined( DEVEL ) 
      WRITE( 6, '(a)' ) REPEAT( '#', 34 )
      WRITE( 6, '(a)' ) '# USING -DDEVEL FOR TESTING GCHP #'
      WRITE( 6, '(a)' ) REPEAT( '#', 34 )
#endif

      
      !-------------------------------------------------
      ! System time stamp
      !-------------------------------------------------
      STAMP = SYSTEM_TIMESTAMP()
      WRITE( 6, 100 ) STAMP
 100  FORMAT( /, '===> SIMULATION START TIME: ', a, ' <===', / )

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
      SUBROUTINE CTM_FLUSH
! 
! !REVISION HISTORY: 
!  31 Aug 2000 - R. Yantosca - Initial version
!  13 Aug 2010 - R. Yantosca - Added ProTeX headers
!  06 Aug 2012 - R. Yantosca - IU_BPCH is only global file LUN still needed
!EOP
!------------------------------------------------------------------------------
!BOC
      CALL FLUSH( IU_BPCH )  

      END SUBROUTINE CTM_FLUSH
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
      SUBROUTINE DISPLAY_END_TIME
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

      END SUBROUTINE DISPLAY_END_TIME
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_initial_met_fields
!
! !DESCRIPTION: Internal subroutine READ\_INITIAL\_MET\_FIELDS calls the
!  various routines to read met fields at the beginning of a GEOS-Chem
!  simulation.  This code was moved out of the main routine for clarity,
!  due to the many \#if defined() blocks that are required.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE READ_INITIAL_MET_FIELDS()
! 
! !REMARKS:
!  All variables used in this routine are declared above in the main 
!  program, and as such, are visible here.
!                                                                             .
!  Also calls the following routines:
!  (1) AVGPOLE   (average pressure @ poles) when I3 or I6 fields are read
!  (2) LIGHTNING (lightning NOx emissions)  when A3 or A6 fields are read
!
! !REVISION HISTORY: 
!  07 Feb 2012 - R. Yantosca - Initial version
!  28 Feb 2012 - R. Yantosca - Removed support for GEOS-3
!  23 Oct 2013 - R. Yantosca - Now pass Input_Opt to GET_A6_FIELDS
!  23 Oct 2013 - R. Yantosca - Now pass Input_Opt to GET_MERRA_A3_FIELDS
!  24 Jun 2014 - R. Yantosca - Now pass Input_Opt to other routines
!  24 Jun 2014 - R. Yantosca - Cosmetic changes, line up arguments
!  12 Aug 2015 - R. Yantosca - Call routines for reading MERRA2 fields
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! Convenience variable for date & time
      INTEGER :: D(2)

#if   defined( GEOS_FP )

      !=================================================================
      !    *****  R E A D   G E O S -- F P   F I E L D S    *****
      !    *****  At the start of the GEOS-Chem simulation  *****
      !
      !    Handle GEOS-FP met fields separately from other met 
      !=================================================================

      ! Read time-invariant data
      CALL GEOSFP_READ_CN  (             Input_Opt, State_Met )

      ! Read 1-hr time-averaged data
      D = GET_FIRST_A1_TIME()
      CALL GEOSFP_READ_A1  ( D(1), D(2), Input_Opt, State_Met )

      ! Read 3-hr time averaged data
      D = GET_FIRST_A3_TIME()
      CALL GEOSFP_READ_A3  ( D(1), D(2), Input_Opt, State_Met )

      ! Read 3-hr time averaged data
      D = GET_FIRST_I3_TIME()
      CALL GEOSFP_READ_I3_1( D(1), D(2), Input_Opt, State_Met )

#elif defined( MERRA2 )

      !=================================================================
      !    *****  R E A D   M E R R A 2   F I E L D S    *****
      !    *****  At the start of the GEOS-Chem simulation  *****
      !
      !    Handle MERRA2 met fields separately from other met 
      !=================================================================

      ! Read time-invariant data
      CALL MERRA2_READ_CN  (             Input_Opt, State_Met )

      ! Read 1-hr time-averaged data
      D = GET_FIRST_A1_TIME()
      CALL MERRA2_READ_A1  ( D(1), D(2), Input_Opt, State_Met )

      ! Read 3-hr time averaged data
      D = GET_FIRST_A3_TIME()
      CALL MERRA2_READ_A3  ( D(1), D(2), Input_Opt, State_Met )

      ! Read 3-hr time averaged data
      D = GET_FIRST_I3_TIME()
      CALL MERRA2_READ_I3_1( D(1), D(2), Input_Opt, State_Met )

#elif defined( MERRA )

      !=================================================================
      !    *****  R E A D   M E R R A   M E T   F I E L D S  *****
      !    *****  At the start of the GEOS-Chem simulation   *****
      !
      !    Handle MERRA met fields separately from other met products
      !=================================================================

      ! Open constant fields
      D = (/ 20000101, 000000 /)
      CALL OPEN_MERRA_CN_FIELDS ( D(1), D(2), Input_Opt               )
      CALL GET_MERRA_CN_FIELDS  ( D(1), D(2), Input_Opt, State_Met    )
      IF ( LPRT ) THEN
         CALL DEBUG_MSG( '### MAIN: a 1st MERRA CN TIME' )
      ENDIF

      ! Open and read A-1 fields
      D = GET_FIRST_A1_TIME()
      CALL OPEN_MERRA_A1_FIELDS ( D(1), D(2), Input_Opt, RESET=.TRUE. )
      CALL GET_MERRA_A1_FIELDS  ( D(1), D(2), Input_Opt, State_Met    )
      IF ( LPRT ) THEN
         CALL DEBUG_MSG( '### MAIN: a 1st MERRA A1 TIME' )
      ENDIF

      ! Open and read A-3 fields
      D = GET_FIRST_A3_TIME()
      CALL OPEN_MERRA_A3_FIELDS ( D(1), D(2), Input_Opt               )
      CALL GET_MERRA_A3_FIELDS  ( D(1), D(2), Input_Opt, State_Met    )
      IF ( LPRT ) THEN
         CALL DEBUG_MSG( '### MAIN: a 1st MERRA A3 TIME' )
      ENDIF

      ! Open & read I-6 fields
      D = GET_FIRST_I6_TIME()
      CALL OPEN_MERRA_I6_FIELDS ( D(1), D(2), Input_Opt               )
      CALL GET_MERRA_I6_FIELDS_1( D(1), D(2), Input_Opt, State_Met    )
      IF ( LPRT ) THEN 
         CALL DEBUG_MSG( '### MAIN: a 1st I6 TIME' )
      ENDIF

#else

      !=================================================================
      !    *****      U N Z I P   M E T   F I E L D S        *****
      !    ***** At at the start of the GEOS-Chem simulation *****
      !
      !   Here we unzip the initial GEOS-3, GEOS-4, GEOS-5, GCAP data
      !=================================================================
      IF ( Input_Opt%LUNZIP ) THEN

         !---------------------
         ! Remove all files
         !---------------------

         ! Type of unzip operation
         ZTYPE = 'remove all'
         
         ! Remove any leftover A-3, A-6, I-6, in temp dir
         CALL UNZIP_A3_FIELDS( Input_Opt, ZTYPE )
         CALL UNZIP_A6_FIELDS( Input_Opt, ZTYPE )
         CALL UNZIP_I6_FIELDS( Input_Opt, ZTYPE )

#if   defined( GCAP )
         ! Unzip GCAP PHIS field (if necessary)
         CALL UNZIP_GCAP_FIELDS( Input_Opt, ZTYPE )
#endif

         !---------------------
         ! Unzip in foreground
         !---------------------

         ! Type of unzip operation
         ZTYPE = 'unzip foreground'

         ! Unzip A-3, A-6, I-6 files for START of run
         CALL UNZIP_A3_FIELDS( Input_Opt, ZTYPE, NYMDb )
         CALL UNZIP_A6_FIELDS( Input_Opt, ZTYPE, NYMDb )
         CALL UNZIP_I6_FIELDS( Input_Opt, ZTYPE, NYMDb )

#if   defined( GCAP )
         ! Unzip GCAP PHIS field (if necessary)
         CALL UNZIP_GCAP_FIELDS( Input_Opt, ZTYPE )
#endif

         !### Debug output
         IF ( LPRT ) CALL DEBUG_MSG( '### MAIN: a UNZIP' )
      ENDIF

      !=================================================================
      !      *****      R E A D   M E T   F I E L D S       *****
      !      ***** At the start of the GEOS-Chem simulation *****
      !
      !  Here we read in the initial GEOS-3, GEOS-4, GEOS-5, GCAP data
      !=================================================================

      ! Open and read A-3 fields
      D = GET_FIRST_A3_TIME()
      CALL OPEN_A3_FIELDS ( D(1), D(2), Input_Opt, RESET=.TRUE. )
      CALL GET_A3_FIELDS  ( D(1), D(2), Input_Opt, State_Met    )
      IF ( LPRT ) THEN
         CALL DEBUG_MSG( '### MAIN: a 1st A3 TIME' )
      ENDIF

      ! Open & read A-6 fields
      D = GET_FIRST_A6_TIME()
      CALL OPEN_A6_FIELDS ( D(1), D(2), Input_Opt               )
      CALL GET_A6_FIELDS  ( D(1), D(2), Input_Opt, State_Met    )      
      IF ( LPRT ) THEN
         CALL DEBUG_MSG( '### MAIN: a 1st A6 TIME' )
      ENDIF

      ! Open & read I-6 fields
      D = GET_FIRST_I6_TIME()
      CALL OPEN_I6_FIELDS ( D(1), D(2), Input_Opt               )
      CALL GET_I6_FIELDS_1( D(1), D(2), Input_Opt, State_Met    )
      IF ( LPRT ) THEN
         CALL DEBUG_MSG( '### MAIN: a 1st I6 TIME' )
      ENDIF

#if   defined( GCAP )

      !-----------------------------------------------------------------
      ! Read additional fields for GCAP meteorology
      !-----------------------------------------------------------------

      ! Read GCAP PHIS and LWI fields (if necessary)
      CALL OPEN_GCAP_FIELDS( Input_Opt            )
      CALL GET_GCAP_FIELDS ( Input_Opt, State_Met )

      ! Remove temporary file (if necessary)
      IF ( LUNZIP ) THEN
         CALL UNZIP_GCAP_FIELDS( Input_Opt, 'remove date' )
      ENDIF

#endif

#endif

#if defined( RRTMG )
      !  Adding the call to retrieve MODIS-derived surface
      !  albedo and emissivity here.
      !  There is a check within READ_SURFACE_RAD to make
      !  sure we dont re-read the same files
      WRITE(6,*) 'IN INITIAL_MET_FIELDS'
      !  Start day may not be one of the 8-day values
      !  Therefore we pass the init flag to ensure
      !  the closest date to the current day is read
      !  rather than nothing at all
      CALL READ_SURFACE_RAD( Input_Opt, FORCEREAD =.TRUE. )
#endif

      END SUBROUTINE READ_INITIAL_MET_FIELDS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_met_fields
!
! !DESCRIPTION: Internal subroutine READ\_MET\_FIELDS calls the
!  various routines to read met fields in the main GEOS-Chem timestepping
!  loop.  This code was moved out of the main routine for clarity, due to
!  the many \#if defined() blocks that are required.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE READ_MET_FIELDS()
!
! !REMARKS:
!  All variables used in this routine are declared above in the main 
!  program, and as such, are visible here.
!                                                                             .
!  Also calls the following routines:
!  (1) AVGPOLE   (average pressure @ poles) when I3 or I6 fields are read
!  (2) LIGHTNING (lightning NOx emissions)  when A3 or A6 fields are read
!
! !REVISION HISTORY: 
!  07 Feb 2012 - R. Yantosca - Initial version
!  28 Feb 2012 - R. Yantosca - Removed support for GEOS-3
!  04 Jan 2013 - M. Payer    - Call UPDATE_T_DAY for MERRA and GEOS-5.7.2 (tmf)
!  23 Oct 2013 - R. Yantosca - Now pass Input_Opt to GET_A6_FIELDS
!  24 Jun 2014 - R. Yantosca - Now pass Input_Opt to other routines
!  24 Jun 2014 - R. Yantosca - Cosmetic changes, line up arguments
!  12 Aug 2015 - R. Yantosca - Now call routines to read MERRA2 met data
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! Convenience variable for date & time
      INTEGER :: D(2)

      ! Set previous specific humidity to current specific humdidity 
      ! prior to MET-field update (ewl, 10/29/15)
      State_Met%SPHU_prev = State_Met%SPHU

#if   defined( GEOS_FP )

      !==============================================================
      !  ****** R E A D   G E O S -- F P   F I E L D S  *****
      !==============================================================

      !---------------------------------
      ! A-1 fields (1hr time averaged)
      !---------------------------------
      IF ( ITS_TIME_FOR_A1() ) THEN
         D = GET_A1_TIME()
         CALL GEOSFP_READ_A1  ( D(1), D(2), Input_Opt, State_Met )
      ENDIF

      !----------------------------------
      ! A-3 fields (3-hr time averaged)
      !----------------------------------
      IF ( ITS_TIME_FOR_A3() ) THEN
         D = GET_A3_TIME()
         CALL GEOSFP_READ_A3  ( D(1), D(2), Input_Opt, State_Met )
      ENDIF

      !----------------------------------
      ! I-3 fields (3-hr instantaneous
      !----------------------------------
      IF ( ITS_TIME_ FOR_I3() ) THEN
         D = GET_I3_TIME()
         CALL GEOSFP_READ_I3_2( D(1), D(2), Input_Opt, State_Met )

         ! Set dry surface pressure (PS2_DRY) from State_Met%PS2_WET
         ! and average as polar caps
         CALL SET_DRY_SURFACE_PRESSURE( State_Met, 2 )
         CALL AVGPOLE( State_Met%PS2_DRY )

         ! Compute avg moist pressure at polar caps 
         CALL AVGPOLE( State_Met%PS2_WET )

      ENDIF

#elif defined( MERRA2 )

      !==============================================================
      !  ****** R E A D   M E R R A 2   F I E L D S  *****
      !==============================================================

      !---------------------------------
      ! A-1 fields (1hr time averaged)
      !---------------------------------
      IF ( ITS_TIME_FOR_A1() ) THEN
         D = GET_A1_TIME()
         CALL MERRA2_READ_A1  ( D(1), D(2), Input_Opt, State_Met )
      ENDIF

      !----------------------------------
      ! A-3 fields (3-hr time averaged)
      !----------------------------------
      IF ( ITS_TIME_FOR_A3() ) THEN
         D = GET_A3_TIME()
         CALL MERRA2_READ_A3  ( D(1), D(2), Input_Opt, State_Met )
      ENDIF

      !----------------------------------
      ! I-3 fields (3-hr instantaneous
      !----------------------------------
      IF ( ITS_TIME_ FOR_I3() ) THEN
         D = GET_I3_TIME()
         CALL MERRA2_READ_I3_2( D(1), D(2), Input_Opt, State_Met )

         ! Set dry surface pressure (PS2_DRY) from State_Met%PS2_WET
         ! and average as polar caps
         CALL SET_DRY_SURFACE_PRESSURE( State_Met, 2 )
         CALL AVGPOLE( State_Met%PS2_DRY )

         ! Compute avg moist pressure at polar caps 
         CALL AVGPOLE( State_Met%PS2_WET )

      ENDIF

#elif defined( MERRA )

      !==============================================================
      !    ***** R E A D   M E R R A   A - 1   F I E L D S *****
      !
      !    The MERRA archive contains hourly surface data fields.
      !==============================================================
      IF ( ITS_TIME_FOR_A1() ) THEN

         ! Get the date/time for the next A-3 data block
         D = GET_A1_TIME()

         ! Open & read A-3 fields
         CALL OPEN_MERRA_A1_FIELDS ( D(1), D(2), Input_Opt            )
         CALL GET_MERRA_A1_FIELDS  ( D(1), D(2), Input_Opt, State_Met )
      ENDIf

      !==============================================================
      !    ***** R E A D   M E R R A   A - 3   F I E L D S *****
      !
      !     The MERRA archive contains 3-hourly 3-D data fields.
      !==============================================================
      IF ( ITS_TIME_FOR_A3() ) THEN
         
         ! Get the date/time for the next A-6 data block
         D = GET_A3_TIME()

         ! Open and read A-6 fields
         CALL OPEN_MERRA_A3_FIELDS ( D(1), D(2), Input_Opt            )
         CALL GET_MERRA_A3_FIELDS  ( D(1), D(2), Input_Opt, State_Met )
      ENDIF

      !==============================================================
      !    ***** R E A D   M E R R A   I - 6   F I E L D S *****
      !
      !    The MERRA archive contains 6-hourly instantaneous data.
      !==============================================================
      IF ( ITS_TIME_FOR_I6() ) THEN

         ! Get the date/time for the next I-6 data block
         D = GET_I6_TIME()

         ! Open and read files
         CALL OPEN_MERRA_I6_FIELDS ( D(1), D(2), Input_Opt            )
         CALL GET_MERRA_I6_FIELDS_2( D(1), D(2), Input_Opt, State_Met )

         ! Set dry surface pressure (PS2_DRY) from State_Met%PS2_WET
         ! and average as polar caps
         CALL SET_DRY_SURFACE_PRESSURE( State_Met, 2 )
         CALL AVGPOLE( State_Met%PS2_DRY )

         ! Compute avg moist pressure at polar caps 
         CALL AVGPOLE( State_Met%PS2_WET )

      ENDIF

#else

      !===============================================================
      !         ***** U N Z I P   M E T   F I E L D S *****
      !
      !      Some met data (except MERRA) are stored compressed.
      !===============================================================
      IF ( LUNZIP .and. ITS_TIME_FOR_UNZIP() ) THEN
         
         ! Get the date & time for 12h (720 mins) from now
         DATE = GET_TIME_AHEAD( 720 )
      
         ! If LWAIT=T then wait for the met fields to be
         ! fully unzipped before proceeding w/ the run.
         ! Otherwise, unzip fields in the background
         IF ( LWAIT ) THEN
            ZTYPE = 'unzip foreground'
         ELSE
            ZTYPE = 'unzip background'
         ENDIF
         
         ! Unzip A3, A6, I6 fields
         CALL UNZIP_A3_FIELDS( Input_Opt, ZTYPE, DATE(1) )
         CALL UNZIP_A6_FIELDS( Input_Opt, ZTYPE, DATE(1) )
         CALL UNZIP_I6_FIELDS( Input_Opt, ZTYPE, DATE(1) )
      
      ENDIF

      !===============================================================
      !        ***** R E M O V E   M E T   F I E L D S *****  
      !===============================================================
      IF ( LUNZIP .and. ITS_TIME_FOR_DEL() ) THEN

         ! Type of operation
         ZTYPE = 'remove date'

         ! Remove A-3, A-6, and I-6 files only for the current date
         CALL UNZIP_A3_FIELDS( Input_Opt, ZTYPE, NYMD )
         CALL UNZIP_A6_FIELDS( Input_Opt, ZTYPE, NYMD )
         CALL UNZIP_I6_FIELDS( Input_Opt, ZTYPE, NYMD )

      ENDIF  

      !==============================================================
      !          ***** R E A D   A - 3   F I E L D S *****
      !
      !  All met data (except MERRA) contain 3-hourly surface data.
      !==============================================================
      IF ( ITS_TIME_FOR_A3() ) THEN

         ! Get the date/time for the next A-3 data block
         D = GET_A3_TIME()

         ! Open & read A-3 fields
         CALL OPEN_A3_FIELDS( D(1), D(2), Input_Opt            ) 
         CALL GET_A3_FIELDS ( D(1), D(2), Input_Opt, State_Met )

      ENDIF

      !==============================================================
      !          ***** R E A D   A - 6   F I E L D S *****  
      !
      !      All other met fields contain 6-hourly 3-D data. 
      !==============================================================
      IF ( ITS_TIME_FOR_A6() ) THEN
         
         ! Get the date/time for the next A-6 data block
         D = GET_A6_TIME()

         ! Open and read A-6 fields
         CALL OPEN_A6_FIELDS( D(1), D(2), Input_Opt            )
         CALL GET_A6_FIELDS(  D(1), D(2), Input_Opt, State_Met )      

      ENDIF

      !==============================================================
      !          ***** R E A D   I - 6   F I E L D S *****   
      !==============================================================
      IF ( ITS_TIME_FOR_I6() ) THEN

         ! Get the date/time for the next I-6 data block
         D = GET_I6_TIME()

         ! Open and read files
         CALL OPEN_I6_FIELDS ( D(1), D(2), Input_Opt            )
         CALL GET_I6_FIELDS_2( D(1), D(2), Input_Opt, State_Met )

         ! Set dry surface pressure (PS2_DRY) from State_Met%PS2_WET
         ! and average as polar caps
         CALL SET_DRY_SURFACE_PRESSURE( State_Met, 2 )
         CALL AVGPOLE( State_Met%PS2_DRY )

         ! Compute avg pressure at polar caps 
         CALL AVGPOLE( State_Met%PS2_WET )

      ENDIF

#endif

#if defined( RRTMG )
      !  Adding the call to retrieve MODIS-derived surface
      !  albedo and emissivity here.
      !  There is a check within READ_SURFACE_RAD to make
      !  sure we dont re-read the same files
      IF ( ITS_TIME_FOR_SURFACE_RAD() ) THEN
!         WRITE(6,*) 'ITS TIME FOR SURFACE RAD'
         CALL READ_SURFACE_RAD( Input_Opt )
      ENDIF
#endif

      END SUBROUTINE READ_MET_FIELDS
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
      SUBROUTINE GET_OVERHEAD_O3_FOR_FASTJ( am_I_Root )
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
            CALL COMPUTE_OVERHEAD_O3( am_I_Root,  
     &                                DAY,
     &                                Input_Opt%USE_O3_FROM_MET, 
     &                                State_Met%TO3 )
         ENDIF
      ENDIF

      END SUBROUTINE GET_OVERHEAD_O3_FOR_FASTJ
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: initialize_regridding
!
! !DESCRIPTION: Internal subroutine Initialize\_Regridding passes several
!  variables to regrid\_a2a\_mod.F90, where they are locally shadowed.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE Initialize_Regridding( am_I_Root, Input_Opt, RC )
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
! !REMARKS:
!  This routine is a wrapper for Init_Map_A2A in regrid_a2a_mod.F90.  
!  Passing variables via Init_Map_A2A helps us to remove dependencies on
!  other GEOS-Chem routines from regrid_a2a_mod.F90.  This in turn 
!  facilitates the implementation of the HEMCO emissions package.
!
! !REVISION HISTORY: 
!  15 Jul 2014 - R. Yantosca - Initial version
!  05 Mar 2015 - R. Yantosca - Now read data w/r/t ExtData/HEMCO
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!     
      ! Scalars
      INTEGER            :: I, J
      CHARACTER(LEN=255) :: DIR

      ! Arrays
      REAL(fpp)          :: LONEDG(IIPAR+1       )  ! W longitude edges [deg]
      REAL(fpp)          :: LATSIN(       JJPAR+1)  ! SIN(Lat edges)    [1  ]
      REAL(fpp)          :: AREAS (IIPAR, JJPAR  )  ! Surface Areas     [m2 ]

      !================================================================
      ! Initialize_Regridding begins here!
      !================================================================

      ! Assume success
      RC  = GC_SUCCESS

      ! Directory where netCDF ifles are found
      DIR = TRIM( Input_Opt%DATA_DIR ) // 'HEMCO/MAP_A2A/v2014-07/'

      ! Initialize longitudes [deg]
      DO I = 1, IIPAR+1
         LONEDG(I) = GET_XEDGE( I, 1, 1 )
      ENDDO

      ! Initialize sines of latitude [1]
      DO J = 1, JJPAR+1
         LATSIN(J) = GET_YSIN( 1, J, 1 )
      ENDDO

      ! Initialize surface areas [m2]
      DO J = 1, JJPAR
      DO I = 1, IIPAR
         AREAS(I,J) = GET_AREA_M2( I, J, 1 )
      ENDDO
      ENDDO

      ! Pass to regrid_a2a_mod.F90, where these will be shadowed locally
      CALL Init_Map_A2A( IIPAR, JJPAR, LONEDG, LATSIN, AREAS, DIR )
      
      END SUBROUTINE Initialize_Regridding
!EOC
      END PROGRAM GEOS_CHEM
#endif
