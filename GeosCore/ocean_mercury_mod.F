! $Id: ocean_mercury_mod.f,v 1.11 2009/09/01 19:21:18 cdh Exp $
      MODULE OCEAN_MERCURY_MOD
!
!******************************************************************************
!  Module OCEAN_MERCURY_MOD contains variables and routines needed to compute
!  the oceanic flux of mercury.  Original code by Sarah Strode at UWA/Seattle.
!  (sas, bmy, 1/21/05, 4/17/06)
!
!  Module Variables:
!  ============================================================================
!  (1 ) Hg_RST_FILE (CHAR   )  : Name of restart file with ocean tracers
!  (2 ) USE_CHECKS  (LOGICAL)  : Flag for turning on error-checking  
!  (3 ) MAX_RELERR  (REAL*8 )  : Max error for total-tag error check [unitless]
!  (4 ) MAX_ABSERR  (REAL*8 )  : Max abs error for total-tag err chk [unitless]
!  (5 ) MAX_FLXERR  (REAL*8 )  : Max error tol for flux error check  [unitless]
!  (6 ) Hg2aq_tot   (REAL*8 )  : Total Hg2 conc. in the mixed layer  [kg      ]
!  (7 ) DD_Hg2      (REAL*8 )  : Array for Hg(II) dry dep'd to ocean [kg      ]
!  (8 ) Hgaq_tot    (REAL*8 )  : Total Hg conc. in the mixed layer   [kg      ]
!  (9 ) Hg0aq       (REAL*8 )  : Array for ocean mass of Hg(0)       [kg      ]
!  (10) Hg2aq       (REAL*8 )  : Array for ocean mass of Hg(II)      [kg      ]
!  (11) HgPaq       (REAL*8 )  : Array for ocean mass of HgP         [kg      ]
!  (12) dMLD        (REAL*8 )  : Array for Change in ocean MLD       [cm      ]
!  (13) MLD         (REAL*8 )  : Array for instantaneous ocean MLD   [cm      ]
!  (14) MLDav       (REAL*8 )  : Array for monthly mean ocean MLD    [cm      ]
!  (15) newMLD      (REAL*8 )  : Array for next month's ocean MLD    [cm      ]
!  (16) NPP         (REAL*8 )  : Array for mean net primary prod.    [unitless]
!  (17) RAD         (REAL*8 )  : Array for mean solar radiation      [W/m2    ]
!  (18) UPVEL       (REAL*8 )  : Array for ocean upwelling velocity  [m/s     ]
!  (19) WD_Hg2      (REAL*8 )  : Array for Hg(II) wet dep'd to ocean [kg      ]
!  (20) CHL         (REAL*8 )  : Chl surface concentration           [mg(m3   ]
!  (21) CDEEPATL    (REAL*8 )  : Conc. Hg0, Hg2, HgP below MLD-Atl   [pM      ]
!  (22) CDEEP       (REAL*8 )  : Conc. of Hg0, Hg2, HgP below MLD    [pM      ]
!  (23) CDEEPNAT    (REAL*8 )  : Conc. Hg0, Hg2, HgP below MLD-NAtl  [pM      ]
!  (24) CDEEPSAT    (REAL*8 )  : Conc. Hg0, Hg2, HgP below MLD-SAtl  [pM      ]
!  (25) CDEEPANT    (REAL*8 )  : Conc. Hg0, Hg2, HgP below MLD-Ant   [pM      ]
!  (26) CDEEPARC    (REAL*8 )  : Conc. Hg0, Hg2, HgP below MLD-Arc   [pM      ]
!
!
!  Module Routines:
!  ============================================================================
!  (1 ) ADD_Hg2_DD             : Archives Hg2 lost to drydep in DD_HG2
!  (2 ) ADD_Hg2_WD             : Archives Hg2 lost to wetdep in WD_HG2
!  (3 ) OCEAN_MERCURY_FLUX     : Routine to compute flux of oceanic mercury
!  (4 ) OCEAN_MERCURY_READ     : Routine to read MLD, NPP, RADSWG data fields
!  (5 ) GET_MLD_FOR_NEXT_MONTH : Routine to read MLD for the next month
!  (6 ) MLD_ADJUSTMENT         : Adjusts MLD 
!  (7 ) READ_OCEAN_Hg_RESTART  : Reads restart file with ocean Hg tracers
!  (8 ) CHECK_DIMENSIONS       : Checks dims of data blocks from restart file
!  (9 ) CHECK_DATA_BLOCKS      : Checks for missing/multiple data blocks
!  (10) MAKE_OCEAN_Hg_RESTART  : Writes new restart file with ocean Hg tracers
!  (11) CHECK_ATMOS_MERCURY    : Checks mass of total & tagged atm Hg0 & Hg2 
!  (12) CHECK_OCEAN_MERCURY    : Checks mass of total & tagged oc Hg0 & Hg2
!  (13) CHECK_OCEAN_FLUXES     : Checks mass of total & tagged DD & WD fluxes
!  (14) INIT_OCEAN_MERCURY     : Allocates and zeroes all module variables
!  (15) CLEANUP_OCEAN_MERCURY  : Deallocates all module variables
!
!  GEOS-CHEM modules referenced by ocean_mercury_mod.f
!  ============================================================================
!  (1 ) bpch2_mod.f            : Module w/ routines for binary pch file I/O
!  (2 ) dao_mod.f              : Module w/ arrays for DAO met fields
!  (3 ) diag03_mod.f           : Module w/ ND03 diagnostic arrays 
!  (2 ) file_mod.f             : Module w/ file unit numbers and error checks
!  (9 ) grid_mod.f             : Module w/ horizontal grid information
!  (10) logical_mod.f          : Module w/ GEOS-CHEM logical switches
!  (11) pressure_mod.f         : Module w/ routines to compute P(I,J,L)
!  (12) time_mod.f             : Module w/ routines to compute date & time
!  (13) tracer_mod.f           : Module w/ GEOS-CHEM tracer array STT etc.
!  (14) tracerid_mod.f         : Module w/ pointers to tracers & emissions
!  (15) transfer_mod.f         : Module w/ routines to cast & resize arrays
!
!  References:
!  ============================================================================
!  (1 ) Xu et al (1999). Formulation of bi-directional atmosphere-surface
!        exchanges of elemental mercury.  Atmospheric Environment 
!        33, 4345-4355.
!  (2 ) Nightingale et al (2000).  In situ evaluation of air-sea gas exchange
!        parameterizations using novel conservative and volatile tracers.  
!        Global Biogeochemical Cycles, 14, 373-387.
!  (3 ) Lin and Tau (2003).  A numerical modelling study on regional mercury 
!        budget for eastern North America.  Atmos. Chem. Phys. Discuss., 
!        3, 983-1015.  And other references therein.
!  (4 ) Poissant et al (2000).  Mercury water-air exchange over the upper St.
!        Lawrence River and Lake Ontario.  Environ. Sci. Technol., 34, 
!        3069-3078. And other references therein.
!  (5 ) Wangberg et al. (2001).  Estimates of air-sea exchange of mercury in 
!        the Baltic Sea.  Atmospheric Environment 35, 5477-5484.
!  (6 ) Clever, Johnson and Derrick (1985).  The Solubility of Mercury and some
!        sparingly soluble mercury salts in water and aqueous electrolyte
!        solutions.  J. Phys. Chem. Ref. Data, Vol. 14, No. 3, 1985.
!  (7 ) Sunderland, E. M. and R. Mason (2007), Human impacts on open 
!        ocean mercury concentrations, Global Biogeochemical Cycles, 21, GB4022,
!        doi:10.1029/2006GB002876, 2007.
!  (8 ) Corbitt, E.S. et al. (2011), Global source-receptor relationsihps for
!       mercury deposition under present-day and 2050 emissions scenarios,
!       Environ. Sci. Technol., 45, 10477-10484, 2011.
!  
!  Nomenclature: 
!  ============================================================================
!  (1 ) Hg(0)  a.k.a. Hg0 : Elemental   mercury
!  (2 ) Hg(II) a.k.a. Hg2 : Divalent    mercury
!  (3 ) HgP               : Particulate mercury
!
!  NOTES:
!  (1 ) Modified ocean flux w/ Sarah's new Ks value (sas, bmy, 2/24/05)
!  (2 ) Now get HALFPOLAR for GCAP or GEOS grids (bmy, 6/28/05)
!  (3 ) Now can read data for both GCAP or GEOS grids (bmy, 8/16/05)
!  (4 ) Include updates from S. Strode and C. Holmes (cdh, sas, bmy, 4/6/06)
!  (5 ) Change HgC (colloidal) to HgP (particulate) or HgPaq. (ccc, 7/20/10)
!  12 Apr 2011 - J. Fisher   - Add missing code from Holmes 2010
!  27 Sep 2011 - H. Amos     - remove LHg_WETDasHNO3 logical, it's obsolete
!  25 Oct 2011 - H. Amos     - bring in code for Hg2 gas-particle partitioning
!  07 Feb 2012 - E. Corbitt  - Imported tagged tracers
!  08 Feb 2012 - R. Yantosca - Add modifications for GEOS-5.7.2 data
!  01 Mar 2012 - R. Yantosca - Now reference new grid_mod.F90
!  19 Apr 2012 - E. Corbitt  - Added LGCBROMINE to use GEOS-Chem bromine.
!  01 Aug 2012 - R. Yantosca - Add reference to findFreeLUN from inqure_mod.F90
!  03 Aug 2012 - R. Yantosca - Move calls to findFreeLUN out of DEVEL block
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!******************************************************************************
!
      USE inquireMod, ONLY : findFreeLUN

      IMPLICIT NONE

      !=================================================================
      ! MODULE PRIVATE DECLARATIONS -- keep certain internal variables 
      ! and routines from being seen outside "ocean_mercury_mod.f"
      !=================================================================

      ! Make everything PRIVATE ...
      PRIVATE

      ! ... except these routines
      PUBLIC :: INIT_OCEAN_MERCURY
      PUBLIC :: CLEANUP_OCEAN_MERCURY
      PUBLIC :: OCEAN_MERCURY_FLUX
      PUBLIC :: READ_OCEAN_Hg_RESTART
      PUBLIC :: MAKE_OCEAN_Hg_RESTART
      PUBLIC :: LDYNSEASALT, LGCAPEMIS, LPOLARBR, LBRCHEM, LBROCHEM
      PUBLIC :: LGEIA05 !(add for 2005 GEIA, eck 10/19/11)
      PUBLIC :: LVEGEMIS ! jaf
      PUBLIC :: LRED_JNO2,   LGEOSLWC
      PUBLIC :: LHg2HalfAerosol
      PUBLIC :: STRAT_BR_FACTOR,        LAnthroHgOnly
      PUBLIC :: LOHO3CHEM,              LnoUSAemis
      PUBLIC :: LGCBROMINE !eds 4/19/12         
      PUBLIC :: READ_HG2_PARTITIONING
      PUBLIC :: Fp, Fg        ! hma 
      PUBLIC :: LNEI2005,               LInPlume     ! yzh
         
      !=================================================================
      ! MODULE VARIABLES
      !=================================================================

      ! Scalars
      LOGICAL              :: USE_CHECKS
      CHARACTER(LEN=255)   :: Hg_RST_FILE

      ! Parameters
      REAL*4,  PARAMETER   :: MAX_RELERR = 5.0d-2
      REAL*4,  PARAMETER   :: MAX_ABSERR = 5.0d-3
      REAL*4,  PARAMETER   :: MAX_FLXERR = 5.0d-1 

      REAL*8   :: CDEEP(3)  
      REAL*8   :: CDEEPATL(3)
      REAL*8   :: CDEEPNAT(3)
      REAL*8   :: CDEEPSAT(3)
      REAL*8   :: CDEEPANT(3)
      REAL*8   :: CDEEPARC(3)
      REAL*8   :: CDEEPNPA(3)

      ! Arrays
      REAL*8,  ALLOCATABLE :: dMLD(:,:)
      REAL*8,  ALLOCATABLE :: Hg0aq(:,:,:)
      REAL*8,  ALLOCATABLE :: Hg2aq(:,:,:)
!- eds 8/31/10 ------------------------------------------------------
!      REAL*8,  ALLOCATABLE :: HgPaq(:,:)
      REAL*8,  ALLOCATABLE :: HgPaq(:,:,:)
      REAL*8,  ALLOCATABLE :: HgPaq_SUNK(:,:,:) !eds 5/15/12 fix
      REAL*8,  ALLOCATABLE :: Hgaq_tot(:,:,:)
!--------------------------------------------------------------------
      REAL*8,  ALLOCATABLE :: MLD(:,:)
      REAL*8,  ALLOCATABLE :: MLDav(:,:)
      REAL*8,  ALLOCATABLE :: newMLD(:,:)
      REAL*8,  ALLOCATABLE :: prevMLD(:,:) !eds 10/19/10 debugging restart
      REAL*8,  ALLOCATABLE :: NPP(:,:)
      REAL*8,  ALLOCATABLE :: RAD(:,:)
      REAL*8,  ALLOCATABLE :: UPVEL(:,:)
      REAL*8,  ALLOCATABLE :: CHL(:,:)                              
! added by hma for Hg2 partitioning --------------------------------------
      REAL*8,  ALLOCATABLE :: BULK_CONC(:,:,:)  
      REAL*8,  ALLOCATABLE :: Fp(:,:,:)     
      REAL*8,  ALLOCATABLE :: Fg(:,:,:)   
      REAL*8,  ALLOCATABLE :: SO4_GC(:,:,:)  
      REAL*8,  ALLOCATABLE :: NH4_CONC(:,:,:)  
      REAL*8,  ALLOCATABLE :: NIT_CONC(:,:,:)
      REAL*8,  ALLOCATABLE :: OC_CONC(:,:,:) 
      REAL*8,  ALLOCATABLE :: BC_CONC(:,:,:)
      REAL*8,  ALLOCATABLE :: DST_CONC(:,:,:) 
      REAL*8,  ALLOCATABLE :: R(:,:,:) 
      REAL*8,  ALLOCATABLE :: SO4_WAC(:,:,:) 
      REAL*8,  ALLOCATABLE :: SO4_CONC(:,:,:) 
!-------------------------------------------------------------------------

      ! Logical switches for the mercury simulation, all of which are 
      ! set in INIT_MERCURY (cdh, 9/1/09)
      LOGICAL   :: LDYNSEASALT, LGCAPEMIS, LPOLARBR, LBRCHEM, LBROCHEM
      LOGICAL   :: LGEIA05 !eck, 10/19/11
      LOGICAL   :: LVEGEMIS ! jaf
      LOGICAL   :: LRED_JNO2,   LGEOSLWC
      LOGICAL   :: LHg2HalfAerosol
      LOGICAL   :: LAnthroHgOnly,          LOHO3CHEM
      LOGICAL   :: LGCBROMINE !eds 4/19/12
      LOGICAL   :: LnoUSAemis
      LOGICAL   :: LNEI2005, LInPlume    ! yzh
      REAL*8    :: STRAT_BR_FACTOR

      ! CDH Set this TRUE to use corrected area-flux relationship
      ! Set this to FALSE to use original Strode et al. (2007) model
      LOGICAL, PARAMETER :: LOCEANFIX=.TRUE. 
      ! CDH average ocean area per grid box: 1.67d11 m2/box
      ! used when eliminating AREA * FRAC_O
      REAL*8, PARAMETER :: FUDGE=1.67D11

      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement
      !=================================================================
      CONTAINS

!-------------------------------------------------------------------------

      SUBROUTINE READ_HG2_PARTITIONING( THISMONTH, State_Met )

!****************************************************************************  
!  SUBROUTINE READ_HG2_PARTITIONING calculates the fractions of Hg(II) is the particle 
!     gas phases.
!
!
!  Arguments as Input:
!  ===========================================================================
!  (1 ) THISMONTH (INTEGER) : Current month number (1-12)
!
!
!  Description of gas-particle partitioning of Hg2 
!  ===========================================================================
!
!  References:
!   (1) Yamasaki et al (1982). Effects of Ambient Temperature on Aspects of 
!       Airbone Polycyclic Aromatic Hydrocarbons, Env Sci & Tech
!   (2) Pankow (1994). An Absorption Model of Gas/Particle Partitioning of
!       Organic Comounds in the Atmosphere, Atmos Env
!   (3) Rutter and Schauer (2007). The effect of temperature on the gas-
!       particle partitioning of reactive mercury in atmospheric aerosols,
!       Atmos Env
!   (4) Vijayaraghavan et al (2008). Plume-in-grid modeling of atmopsheric 
!       mercury   
!   (5) Amos et al. (2012, ACPD). Gas-particle partitioning of Hg(II)
!       and its effect on global mercury deposition
!
!  Ratio of reactive mercury adsorbed onto particulate matter to reactive
!  mercury in the gas phase:
! 
!  (PHg,ads)/RGM = 10^(b/T - a)* PM     
! 
!      PHg,ads  =  adsorbed RGM                        (pg m^-3)
!      RGM      =  reactive gaseous mercury            (pg m^-3)
!      T        =  temperature                         (K)
!      PM       =  ambient aerosol concentration       (ug m^-3)
!      b        =  slope from simple linear regression
!      a        =  y-intercept from simple linear regression
!
!
!  Aerosol concentrations are being taken from a GEOS-Chem v8-02-03,
!  GEOS-5, 4x5, full-chem simulation run for 2007 by Lin Zhang. The units
!  reported by GEOS-Chem are ppb for aerosol mixing ratio (IJ-AVG-$).
!  Units must be converted from mol/mol to ug/m3. 
!  
!  Converting aerosol concentration ppbv --> ug/m3: 
!   ( Modeled after SUBROUTINE CONVER_UNITS in dao_mod.f )
!
!      1 ppbv = 1e-9 mol/mol
!      AIRDEN = air density,  [kg/m3]
!      CU     = aerosol molecular weight / molecular weight of air
!             = [(aero kg/mol) / (air kg/mol)]
!
!            aero mol     air kg   aero kg/mol     1e9 ug
!      PM = ---------- x ------- x ------------ x -------
!            air  mol     air m3    air kg/mol       kg
!
!      PM = (IJ-AVG-$) * AIRDEN * CU
! 
!
!  The objective of this subroutine is to determine the fraction of reactive
!  mercury in the gas-phase (Fg) and the fraction in the particle-phase (Fp).
!  That can be done easily once R is calculated, where R is the ratio of
!  particle to gas:
!
!      R  = (HgP,ads)/RGM                              (unitless)
!      Fg = 1/(R + 1)                                  (unitless)
!      Fp = 1 - Fg   *or*  = R/(R + 1)                 (unitless)
!
!  Notes:
!  ===========================================================================
!  (1 ) Scaling factor for dust emissions --> Fairlie et al., 'Impact of
!       mineral dust on nitrate, sulfate, and ozone in transpacific Asian
!       pollution plumes'. (2009) - in preparation. [hma, 22 Oct 09]
!  (2 ) SUBROUTINE READ_HG2_PARTITION moved from mercury_mod.f to ocean_mercury_mod.f!    
!       because partitioning needs to happen before vdiff_mod.f90 (in terms of
!       the order in which modules are compiled.)
!  (3 ) Use benchmark V8-02-01 geos5 runs for aerosol concentration (ppbv)
!       and aerosol surface area (cm2/cm3). hma 3 Feb 10
!  (4 ) Exclude sea-salt from BULK_CONC and BULK_SA. Chris already takes care
!       of Hg uptake on sea-salt aerosols in a different routine. 
!       (hma 6 Feb 10)
!  (5 ) R is no longer computed using aerosol surface area. (hma 20100713)
!  (6 ) Exclude coarse dust (dst2-4), sea-salt SO4 (SO4s), and sea-salt NO3 (NITS).
!         (hma 20100712)
!  23 Feb 2011 - H Amos - The aerosol mixing ratios read from file here were
!                         run with 54 tracers. The 54-tracer run is now obsolete
!                         and has been extended to 59 tracers.
!  09 Mar 2011 - H Amos - clean up comments
!  25 Oct 2011 - H Amos - brought into v9-01-02
!  09 Nov 2012 - M. Payer    - Replaced all met field arrays with State_Met
!                              derived type object
!*****************************************************************************
     
      ! References F90 modules
      USE BPCH2_MOD,          ONLY : GET_NAME_EXT ,  GET_RES_EXT
      USE BPCH2_MOD,          ONLY : GET_TAU0,       READ_BPCH2
      USE DIAG03_MOD,         ONLY : ND03                
      USE DIRECTORY_MOD,      ONLY : DATA_DIR
      USE DIRECTORY_MOD,      ONLY : DATA_DIR      
      USE ERROR_MOD,          ONLY : GEOS_CHEM_STOP, ERROR_STOP
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE TIME_MOD,           ONLY : SET_Hg2_DIAG        
      USE TRANSFER_MOD,       ONLY : TRANSFER_3D  ,  TRANSFER_3D_TROP

      USE CMN_SIZE_MOD             ! Size parameters3
      USE CMN_DIAG_MOD             ! Diagnostic arrays & parameters

      ! Arguments
      INTEGER,        INTENT(IN)  :: THISMONTH
      TYPE(MetState), INTENT(IN)  :: State_Met   ! Meteorology State object

      ! Arrays to hold aerosol concentrations  
      REAL*4               :: ARRAYso4(  IGLOB, JGLOB, LGLOB )  ! so4
      REAL*4               :: ARRAYnit(  IGLOB, JGLOB, LGLOB )  ! nit
      REAL*4               :: ARRAYnh4(  IGLOB, JGLOB, LGLOB )  ! nh4
      REAL*4               :: ARRAYbcpi( IGLOB, JGLOB, LGLOB )  ! bcpi
      REAL*4               :: ARRAYocpi( IGLOB, JGLOB, LGLOB )  ! ocpi
      REAL*4               :: ARRAYbcpo( IGLOB, JGLOB, LGLOB )  ! bcpo
      REAL*4               :: ARRAYocpo( IGLOB, JGLOB, LGLOB )  ! ocpo
      REAL*4               :: ARRAYdst1( IGLOB, JGLOB, LGLOB )  ! dst1

      ! Arrays to hold bulk concentration and surface area
      REAL*4               :: ARRAYconc(IGLOB,JGLOB,LGLOB)  ! concentration

      ! Local variables
      REAL*8               :: XTAU
      INTEGER              :: I, J, L, NN
      CHARACTER(LEN=155)   :: FILENAME
      CHARACTER(LEN=2)     :: MON                 ! for months 10-12
      CHARACTER(LEN=1)     :: SMON                ! for months 1-9
      REAL*8, PARAMETER    :: TINY  = 1d-16       ! to prevent NaN
      REAL*8, DIMENSION(5) :: AERO_MW             ! to hold aerosol molec weights 
      REAL*8               :: AIR_MW = 2.897D-02  ! avg molecular weight of air 
      REAL*8, DIMENSION(5) :: CU                  ! ratio of aerosol/air MW

      !----------------------------------------!
      ! Molecular weights, for unit conversion !
      !----------------------------------------!

      ! aerosol molecular weights, [kg/mol]
      AERO_MW(1) = 9.6D-02    ! SO4, SO4s
      AERO_MW(2) = 1.2D-02    ! OC, BC
      AERO_MW(3) = 6.2D-02    ! NIT, NITs
      AERO_MW(4) = 1.8D-02    ! NH4
      AERO_MW(5) = 2.9D-02    ! DST 

      ! ratio of aerosol molecular weight / air molecular weight
      CU = AERO_MW / AIR_MW   !  [kg aerosol / kg air]

#if defined( GRID4x5 ) || defined( GRID2x25 )
      ! tropospheric aerosols (v8-02-03, GEOS5, 4x5)
      FILENAME = TRIM( DATA_DIR )   //
     &           'mercury_201007/Hg_partitioning/' //
     &           'aero.IJ_AVG.2007.'               //
     &           GET_NAME_EXT()     // '.'         //
     &           GET_RES_EXT()
#else
      CALL ERROR_STOP( "You need to make an aerosols file. ", 
     &                 "PARTITION_Hg2")
#endif

      ! Echo some information to the standard output
      WRITE( 6, 110 )  FILENAME 
 110  FORMAT( '     -READ_HG2_PARTITIONING: Reading ', a )

      !------------------------------------------------------!
      ! Put aerosol mixing ratio (mol/mol) into arrays.     
      !------------------------------------------------------!

      !---------------------------
      ! Tropospheric aerosols
      !---------------------------
      
      ! Get the TAU0 value for the start of the given month
      ! Assume "generic" year 1985 (TAU0 = [0, 744, ... 8016])
      XTAU = GET_TAU0( THISMONTH, 1, 2007 )  


      ! so4, tracer = 27
      CALL READ_BPCH2( FILENAME, 'IJ-AVG-$', 27,      
     &                 XTAU,      IGLOB,     JGLOB,
     &                 LGLOB,     ARRAYso4,  QUIET=.TRUE. )

      ! convert REAL*4 to REAL*8 
      CALL TRANSFER_3D_TROP( ARRAYso4, SO4_GC ) 
     
      ! nh4, tracer = 31
      CALL READ_BPCH2( FILENAME, 'IJ-AVG-$', 31,      
     &                 XTAU,      IGLOB,     JGLOB,
     &                 LGLOB,     ARRAYnh4,  QUIET=.TRUE. )

      ! convert REAL*4 to REAL*8 
      CALL TRANSFER_3D_TROP( ARRAYnh4, NH4_CONC ) 

      ! nit, tracer = 32
      CALL READ_BPCH2( FILENAME, 'IJ-AVG-$', 32,      
     &                 XTAU,      IGLOB,     JGLOB,
     &                 LGLOB,     ARRAYnit,  QUIET=.TRUE. )

      ! convert REAL*4 to REAL*8 
      CALL TRANSFER_3D_TROP( ARRAYnit, NIT_CONC ) 

      ! bcpi, tracer = 34
      CALL READ_BPCH2( FILENAME, 'IJ-AVG-$', 34,      
     &                 XTAU,      IGLOB,     JGLOB,
     &                 LGLOB,     ARRAYbcpi, QUIET=.TRUE. )

      ! bcpo, tracer = 36
      CALL READ_BPCH2( FILENAME, 'IJ-AVG-$', 36,      
     &                 XTAU,      IGLOB,     JGLOB,
     &                 LGLOB,     ARRAYbcpo, QUIET=.TRUE. )

      ! convert REAL*4 to REAL*8 
      CALL TRANSFER_3D_TROP( ARRAYbcpi + ARRAYbcpo , BC_CONC ) 

      ! ocpi, tracer = 35
      CALL READ_BPCH2( FILENAME, 'IJ-AVG-$', 35,      
     &                 XTAU,      IGLOB,     JGLOB,
     &                 LGLOB,     ARRAYocpi, QUIET=.TRUE. )

      ! ocpo, tracer = 37
      CALL READ_BPCH2( FILENAME, 'IJ-AVG-$', 37,      
     &                 XTAU,      IGLOB,     JGLOB,
     &                 LGLOB,     ARRAYocpo, QUIET=.TRUE. )

      ! convert REAL*4 to REAL*8 
      CALL TRANSFER_3D_TROP( ARRAYocpi + ARRAYocpo , OC_CONC ) 

      ! dst1, tracer = 38   
      CALL READ_BPCH2( FILENAME, 'IJ-AVG-$', 49,  
     &                 XTAU,      IGLOB,     JGLOB,
     &                 LGLOB,     ARRAYdst1, QUIET=.TRUE. )
 
      ! convert REAL*4 to REAL*8 
      CALL TRANSFER_3D_TROP( ARRAYdst1 , DST_CONC ) 

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )

      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR

        ! Sum aerosol types (ppbv) to get bulk aerosol                     
        !                                                                  
        ! - Multiply species by CU factor as part of unit conversion from  
        !    ppbv --> ug/m3 
        ! Scale dust from 4x5 run (ref: Fairlie et al. (2010))            
        BULK_CONC(I,J,L) = CU(1)*SO4_CONC(I,J,L)                   +
     &                     CU(2)*(BC_CONC(I,J,L) + OC_CONC(I,J,L)) +
     &                     CU(3)*NIT_CONC(I,J,L)                   +
     &                     CU(4)*NH4_CONC(I,J,L) +
     &                     CU(5)*(DST_CONC(I,J,L)*2d0)           
      
        ! convert bulk aerosol mass concentration  to ug/m3
        ! Note: Even though both bpch files say the unit is 'ppb',
        !       4x5 gets read in as mol/mol and 2x2.5 gets read in
        !       as ppb, so only 4x5 needs to be scaled by 1d9.
        !       (hma, 20120104)
        ! 0..5x0.666 file is regridded from 4x5 file, yzh, 12/15/12
#if defined( GRID4x5 ) || defined( GRID05x0666 )
        BULK_CONC(I,J,L) = 1d9 * ( BULK_CONC(I,J,L) *
     &                     State_Met%AIRDEN(L,I,J)  )
#elif defined( GRID2x25 )
        BULK_CONC(I,J,L) = ( BULK_CONC(I,J,L) * State_Met%AIRDEN(L,I,J))
#endif 

        ! Calculate R = HgP_ads/RGM (i.e. the ratio of Hg2
        ! adsorbed onto aerosol to Hg2 in the  gas phase) 
        R(I,J,L) = BULK_CONC(I,J,L)*
     &           ( 10d0**( ( 2.5d3 / State_Met%T(I,J,L)) - 10d0 ))

        ! Fraction of Hg(II) in the gas phase (unitless)
        Fg(I,J,L) = 1d0 / (R(I,J,L) + 1d0)   
            
        ! Fraction of Hg(II) in the particle phase (unitless)
        Fp(I,J,L) = 1d0 - Fg(I,J,L)
  
      ENDDO
      ENDDO
      ENDDO 

!$OMP END PARALLEL DO

      !-------------------------------------------!
      ! Store Fg and Fp in ND03 diagnostic.       ! 
      !-------------------------------------------!

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR

            ! error check 
            IF (Fg(I,J,L)+Fp(I,J,L) > 1) THEN
               PRINT*, ' Fg + Fp > 1 @ ocean_mercury_mod.f'
               CALL GEOS_CHEM_STOP 
            ENDIF

      ENDDO
      ENDDO
      ENDDO 
!$OMP END PARALLEL DO


      ! Increment diagnostic timestep counter. 
      CALL SET_Hg2_DIAG( INCREMENT=.TRUE. )

      ! Return to calling program.
      END SUBROUTINE READ_HG2_PARTITIONING

!------------------------------------------------------------------------------

      SUBROUTINE OCEAN_MERCURY_FLUX( FLUX,      Input_Opt, 
     &                               State_Met, State_Chm )
!
!******************************************************************************
!  Subroutine OCEAN_MERCURY_FLUX calculates emissions of Hg(0) from 
!  the ocean in [kg/s].  (sas, bmy, 1/19/05, 4/17/06)
!
!  NOTE: The emitted flux may be negative when ocean conc. is very low. 
!
!  ALSO NOTE: The ocean flux was tuned with GEOS-4 4x5 met fields.  We also
!  now account for the smaller grid size if using GEOS-4 2x25 met fields.
!    
!  Arguments as Output
!  ============================================================================
!  (1 ) FLUX (REAL*8) : Flux of Hg(0) from the ocean [kg/s]
!_____________________________________________________________________________
!
!  GENERAL SOLUTION - OXIDATION, REDUCTION, SINKING, EVASION, UPWELLING
!
!  dHg0/dt  = Hg0(upw) + Hg0(ent) + Hg0(oa) -k_ox
!             + k_red * Frac_Hg2 * Reducible * HgII
!
!  dHgII/dt = HgII(dep) + HgII(up) + HgII(ent) - HgII(sink) +
!             k_ox * Hg0-k_red * Frac_Hg2 * Reducible * HgII
!____________________________________________________________________________
!
! Hg(tot)aq REDUCTION RATE CONSTANTS
!
! Hg(tot)aq reduction is split into biological and radiative reduction
!  (1.1 added to NPP for abiotic particles)
!
!   k_red     = k_red_bio + k_red_rad
!   k_red_rad = k_radbase * RADz     = ( s-1 W-1 m2 ) * ( W m-2 )
!   k_red_bio = k_biobse * NPP * 1.1 = ( s-1 mgC-1 d ) * ( mgC m-2 d-1 )
!
!
! Hg(0)aq OXIDATION RATE CONSTANTS
!
!   k_ox      = k_oxbase * RADz + k_dark
!
! k_dark is a constant dark oxidation component
!
! RADz is the integrated ligth attenuation based on Beer-Lamberts law 
! (Schwarzenbach et al. 1993)
!
!    RADz = (1/(x1-x2))(RAD/EC)(1-e**-EC * x2)
! 
! x1  = surface depth (=0) (m)
! x2  = depth of mixed layer (m)
! EC  = extinction coefficient (m-1)
! RAD = incomming radiation from GEOS5
!
! Extinction coefficient
! EC = ECwater + ECdoc * Cdoc (NPP/NPPavg) + ECchla * CHL/1000
!
! ECwater = 0.0145 m-1
! ECdoc   = 0.654 m-1
! Cdoc    = 1.5 mgL-1
! ECchla  = 31 m-1
! CHL     = amount dependent on inputfile (mg/m3) but we need 
!           mg/L so divide CHL by 1000
!____________________________________________________________________________
!
! TOTAL ORGANIC CARBON AND SUSPENDED PARTICULATE MATTER (TOTAL BIOMASS)
!
! Hg(II) - Hg(P) partitioning coefficient
!
!   Fraction of Hg2 = Frac_Hg2 =  1 / ( 1 + kd_part * SPM )
!
! Kd_part is based on Mason et al. 1998 and Mason & Fitzgerald 1993. (L/kg)
! SPM is converted to kg/L by 10E-9
!
! SPM is Suspended particulate matter (kg/L)
!
!   SPM = ( OC_tot * 10 / MLD ) * 1.1
!
! Total biomas is a proxy for SPM (mg/m3) used in Hg(II)
! partitioning. Calculated by multiplying the standing 
! stock of organic carbon (OC_tot) with 10 (exp Bundy 2004)
! 1.1 is to include abiotic particles
!
! OC_tot is the standing stock of organic carbon (mgC/m²)
!
!   OC_tot = C_tot * 80
!
! Standing stock is calculated based on C:Chl ratio of 80 (wetzel et al 2006)
!
! C_tot is the integrated pigment content in euphotic layer (mg/m2)
!
! The parameters for calculating integrated Chl is based on a
! model by Uitz et al (2006).
!
! CHL   = average Chl a conc. detected by Modis (mg/m3)
! Zm    = mixed layer depth (m)
! Ze    = euphotic depth (PAR 1% of surface value (m)
!
! C_tot differs dependent on the water being stratified 
! or well-mixed.
!___________________________________________________________________________
!
! GAS EXCHANGE
!
! Net flux from the ocean is given by the equation:
!  
!   F = Kw * ( CHg0_aq - CHg0_atm / H )    (Lis & Slanter 1974)
!
! Kw is the mass transfer coefficient (cm/h)
!   There are different possibilities for calculating Kw. The default is:
!
!   Kw = 0.25 * u^2 / SQRT ( Sc / ScCO2 )  (Nightingale et al. 2000)
!
! u^2 is the square of the wind speed (10m above ground) (m²/s²)
!
! Sc is the Schmidt # for Hg [unitless]                             
!    (ref: Poissant et al 2000; Wilke and Chang 1995)
!    to correct for seawater D0 is decreased by 6% as suggested
!    by Wanninkhof (1992)
!
!   Sc = v/D = (0.017 * exp(-0.025T))/D = kinematic viscosity/diffusivity
!
! Diffusivity is calculated by:
!   D = (7.4*10D-8 scrt(2.26 * Mw) * TK) / (vi * N**0.6)
!
!   vi = viscocity of water
!   N  = molal volumen of mercury = 14.18
!
! Viscocity is taken from Loux (2001)
!
! H is the diemensionless Henrys coefficient for elemental mercury
!
!   H = exp (-2404.3/T - 6.92) where T is sea temp in K (Andersson et al. 2008)  
!___________________________________________________________________________
!
! PARTICLE SINKING
!
! (from Sunderland & Mason 2007)
!
! JorgC_kg = 0.1 (NPP**1.77) (MLD**-0.74)
!____________________________________________________________________________
!
!  NOTES:
!  (1 ) Change Ks to make ocean flux for 2001 = 2.03e6 kg/year.
!        (sas, bmy, 2/24/05)
!  (2 ) Rewritten to include Sarah Strode's latest ocean Hg model code.
!        Also now accounts for 2x25 grid. (sas, cdh, bmy, 4/6/06)
!  (3 ) Ocean parameterizations are rewritten entirely to account for actual
!        processes in the ocean. Different subsurface conc. are included
!        (anls, 20/10/09)
!  (4 ) Now use MERRA land fraction information (jaf, 4/26/11)
!  12 Apr 2011 - J. Fisher   - Add missing code from Holmes 2010
!  27 Jul 2011 - R. Yantosca - Use the D exponent in the defnition of KD_PART, 
!                              or else the PGI compiler will crash
!  09 Nov 2011 - H. Amos     - Apply Anne's bug fix for the Hg(0) evasion 
!                              and uptake diagnostics
!  08 Feb 2012 - R. Yantosca - Now treat GEOS-5.7.2 in the same way as MERRA
!  01 Mar 2012 - R. Yantosca - Now use GET_AREA_CM2(I,J,L) from grid_mod.F90
!  01 Mar 2012 - R. Yantosca - Now use GET_XMID(I,J,L) from grid_mod.F90
!  01 Mar 2012 - R. Yantosca - Now use GET_YMID(I,J,L) from grid_mod.F90
!  09 Nov 2012 - M. Payer    - Replaced all met field arrays with State_Met
!                              derived type object
!  26 Sep 2013 - R. Yantosca - Renamed GEOS_57 Cpp switch to GEOS_FP
!******************************************************************************
!
      ! References to F90 modules
      USE DAO_MOD,            ONLY : IS_WATER 
      USE DEPO_MERCURY_MOD,   ONLY : DD_Hg2, WD_Hg2, DD_HgP, WD_HgP
      USE DIAG03_MOD,         ONLY : AD03, ND03
      USE DIRECTORY_MOD,      ONLY : DATA_DIR
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE GIGC_State_Chm_Mod, ONLY : ChmState
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE GRID_MOD,           ONLY : GET_AREA_M2, GET_XMID, GET_YMID 
      USE TIME_MOD,           ONLY : GET_TS_EMIS,     GET_MONTH 
      USE TIME_MOD,           ONLY : ITS_A_NEW_MONTH, ITS_MIDMONTH
      USE TIME_MOD,           ONLY : GET_YEAR 
!- eds 8/31/10 -------------------------------------------------------------
!      USE TRACERID_MOD,  ONLY : ID_Hg_tot,       ID_Hg_oc
      USE TRACERID_MOD,       ONLY : ID_Hg_tot
      USE TRACERID_MOD,       ONLY : ID_Hg_atl,  ID_Hg_nat,  ID_Hg_sat
      USE TRACERID_MOD,       ONLY : ID_Hg_npa,  ID_Hg_arc,  ID_Hg_ant
      USE TRACERID_MOD,       ONLY : ID_Hg_ocn
!---------------------------------------------------------------------------
      USE TRACERID_MOD,       ONLY : ID_Hg0,          N_Hg_CATS
      USE TRANSFER_MOD,       ONLY : TRANSFER_2D

      USE CMN_SIZE_MOD             ! Size parameters
!
! !INPUT PARAMETERS:
!

      TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
      TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology State object
!
! !INPUT/OUTPUT PARAMETERS:
!
      TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
!
! !OUTPUT PARAMETERS:
!
      REAL*8,         INTENT(OUT) :: FLUX(IIPAR,JJPAR,N_Hg_CATS)
!
! !LOCAL VARIABLES:
!
      ! Local variables
      LOGICAL, SAVE         :: FIRST = .TRUE.
      LOGICAL               :: IS_OCEAN_BOX ! jaf
      CHARACTER(LEN=255)    :: FILENAME
      INTEGER               :: I,         J,        NN, C
      INTEGER               :: N,         N_tot_oc
      INTEGER               :: NEXTMONTH, THISMONTH
      INTEGER               :: THISYEAR

      REAL*8                :: A_M2,     DTSRCE,   MLDCM
      REAL*8                :: CHg0aq,   CHg0,     vi,       JorgC_kg           
      REAL*8                :: TC,       TK,       Kw
      REAL*8                :: Sc,       ScCO2,    USQ,      MHg
      REAL*8                :: Hg2_RED,  Hg2_GONE, Hg2_CONV  !HgPaq_SUNK
      REAL*8                :: FRAC_L,   FRAC_O,   H,        TOTDEP
      REAL*8                :: oldMLD,   XTAU,     TOTDEPall                   
      REAL*8                :: FUP(IIPAR,JJPAR,N_Hg_CATS)
      REAL*8                :: FDOWN(IIPAR,JJPAR,N_Hg_CATS)
      REAL*8                :: X,        Y,        D                           
      REAL*8                :: NPP_tot,  A_ocean,  NPP_avg,  RADz        
      REAL*8                :: EC         
      REAL*8                :: k_red,    k_red_rad,  k_red_bio          
      REAL*8                :: k_ox    
      REAL*8                :: SPM,      Frac_Hg2, OC_tot_kg          
      !REAL*8                :: Hgaq_tot 
      REAL*8                ::Hg2aq_tot                      
      REAL*8                :: C_tot,    Ze,       OC_tot,   Hg0_OX     

      ! Parameters
      REAL*8, PARAMETER     :: EC_w      = 0.0145d0       
      REAL*8, PARAMETER     :: EC_doc    = 0.654d0
      REAL*8, PARAMETER     :: C_doc     = 1.5d0 
      !REAL*8, PARAMETER     :: k_radbase = 1.73d-6    
      !decrease photoreduction in atm & surf ocean
      ! per H. Amos 23 Sep 2011 to multiply by 0.9d0
      ! implemented by eck 10/19/11
      REAL*8, PARAMETER     :: k_radbase = 1.557d-6
      REAL*8, PARAMETER     :: k_biobase = 4.1d-10    
      REAL*8, PARAMETER     :: k_oxbase  = 6.64d-6 
#if defined ( LINUX_PGI )
      REAL*8, PARAMETER     :: Kd_part   = 316227.76601684
#else
      REAL*8, PARAMETER     :: Kd_part   = 10d0**(5.5d0)
#endif
      REAL*8, PARAMETER     :: k_ox_dark = 1d-7  
      REAL*8, PARAMETER     :: ECchla    = 31d0     

      ! Conversion factor from [cm/h * ng/L] --> [kg/m2/s]
      REAL*8,  PARAMETER    :: TO_KGM2S = 1.0D-11 / 3600D0 

      ! Small numbers to avoid dividing by zero
      REAL*8,  PARAMETER   :: SMALLNUM   = 1D-32
      REAL*8,  PARAMETER   :: NPPMINNUM   = 5D-2       
      REAL*8,  PARAMETER   :: CHLMINNUM   = 1D-1       

      ! External functions
      REAL*8,  EXTERNAL     :: SFCWINDSQR 
     
      ! For values from Input_Opt
      LOGICAL         :: LSPLIT
      INTEGER         :: N_TRACERS
      REAL*8          :: TRACER_MW_KG(Input_Opt%N_TRACERS)

      ! Pointers
      ! We need to define local arrays to hold corresponding values 
      ! from the Chemistry State (State_Chm) object. (mpayer, 12/6/12)
      REAL*8, POINTER :: STT(:,:,:,:)

      !=================================================================
      ! OCEAN_MERCURY_FLUX begins here!
      !=================================================================

      ! Copy values from Input_Opt
      LSPLIT       = Input_Opt%LSPLIT
      N_TRACERS    = Input_Opt%N_TRACERS
      TRACER_MW_KG = Input_Opt%TRACER_MW_KG(1:N_TRACERS)

      ! Loop limit for use below
      IF ( LSPLIT ) THEN
         N_tot_oc = 2
      ELSE
         N_tot_oc = 1
      ENDIF

      ! Molecular weight of Hg (applicable to all tagged tracers)
      MHg = TRACER_MW_KG(ID_Hg_tot)

      !-----------------------------------------------
      ! Check tagged & total sums (if necessary)
      !-----------------------------------------------
      IF ( USE_CHECKS .and. LSPLIT ) THEN
         CALL CHECK_ATMOS_MERCURY( State_Chm,
     &                             'start of OCEAN_MERCURY_FLUX' )
         CALL CHECK_OCEAN_MERCURY( 'start of OCEAN_MERCURY_FLUX' )
         CALL CHECK_OCEAN_FLUXES ( 'start of OCEAN_MERCURY_FLUX' )
      ENDIF

      !-----------------------------------------------
      ! Read monthly NPP, RADSW, MLD, UPVEL data
      !-----------------------------------------------
      IF ( ITS_A_NEW_MONTH() ) THEN

         ! Get current month
         THISMONTH = GET_MONTH()

         ! Get monthly MLD, NPP, CHL etc.                         
         !CALL OCEAN_MERCURY_READ( THISMONTH )
         CALL OCEAN_MERCURY_READ( THISMONTH, THISYEAR ) !jaf

      ENDIF    

      !eds 10/19/10
      IF ( FIRST ) THEN
         CALL GET_MLD_FOR_LAST_MONTH( THISMONTH )
         FIRST = .FALSE.
      ENDIF


      !-----------------------------------------------
      ! MLD and entrainment change in middle of month
      !-----------------------------------------------
      IF ( ITS_MIDMONTH() ) THEN

         ! Get current month
         THISMONTH = GET_MONTH()

!--jaf.start
         ! Get current year to check if leap year (jaf, 7/6/11)
         THISYEAR = GET_YEAR()

         ! Read next month's MLD
         !CALL GET_MLD_FOR_NEXT_MONTH( THISMONTH )
         CALL GET_MLD_FOR_NEXT_MONTH( THISMONTH, THISYEAR )
!--jaf.end
         
      ENDIF


      ! Emission timestep [s]
      DTSRCE = GET_TS_EMIS() * 60d0

      !----------------------------------------------------------------
      ! Calculate total mean NPP (mg/m2/day) for later                                                  
      !----------------------------------------------------------------                                     
      ! Initialize values
      NPP_tot = 0d0
      A_ocean = 0d0

      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Grid box surface area [m2]
         A_M2 = GET_AREA_M2( I, J, 1 )

!--jaf.start
         ! Use fractional land type information in MERRA (jaf, 4/26/11)
         ! FROCEAN is a constant, so to get correct ocean fraction we
         ! need to subtract the sea ice fraction.
#if   defined( MERRA ) || defined( GEOS_FP )
         NPP_tot = NPP_tot+ NPP(I,J) * A_M2 * 
     &             ( State_Met%FROCEAN(I,J) - State_Met%FRSEAICE(I,J) )
         A_ocean = A_ocean + A_M2 * 
     &             ( State_Met%FROCEAN(I,J) - State_Met%FRSEAICE(I,J) )
#else
         NPP_tot = NPP_tot+ NPP(I,J) * A_M2 *
     &             ( 1d0 - State_Met%FRCLND(I,J))
         A_ocean = A_ocean + A_M2 * ( 1 - State_Met%FRCLND(I,J) ) 
#endif

      ENDDO
      ENDDO

      NPP_avg = NPP_tot / A_ocean

      ! Initialize GEOS-Chem tracer array [kg] from Chemistry State object
      ! (mpayer, 12/6/12)
      STT => State_Chm%Tracers

      ! Loop over surface boxes  
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I,   vi,   A_M2,    Hg2_RED  )
!$OMP+PRIVATE( J,   NN,   k_ox,    OC_tot,  Hg2_CONV   )
!$OMP+PRIVATE( N,   TK,   CHg0,    FRAC_L,  k_red_bio  )
!$OMP+PRIVATE( C,   TC,   RADz,    Hg0_OX,  k_red_rad  )
!$OMP+PRIVATE( D,   EC,   k_red,   OLDMLD,  TOTDEPall  )
!$OMP+PRIVATE( Y,   Ze,   ScCO2,   FRAC_O,  Frac_Hg2,   Hg2aq_tot  )
!$OMP+PRIVATE( H,   Kw,   MLDCM,   TOTDEP,  OC_tot_kg  )
!$OMP+PRIVATE( X,   SPM,  CHg0aq,  Hg2_GONE   )
!$OMP+PRIVATE( Sc,  Usq,  C_tot,   JorgC_kg   )
!$OMP+PRIVATE( IS_OCEAN_BOX                   )
!$OMP+SCHEDULE( DYNAMIC )
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Grid box surface area [m2]
         A_M2       = GET_AREA_M2( I, J, 1 )

         ! Initialize values
         Kw         = 0d0
!         HgPaq_SUNK = 0d0
         Hg2_CONV   = 0d0
         TK         = 0d0
         TC         = 0d0 
         JorgC_kg   = 0d0                              
         EC         = 0d0
         RADz       = 0d0
         k_red      = 0d0
         k_red_rad  = 0d0
         k_red_bio  = 0d0
         SPM        = 0d0
         Frac_Hg2   = 0d0   
         Hg2aq_tot  = 0d0
!         Hgaq_tot   = 0d0
         Hg2_RED    = 0d0
         C_tot      = 0d0
         Ze         = 0d0
         OC_tot     = 0d0
         OC_tot_kg  = 0d0 
         Hg0_OX     = 0d0
         D          = 0d0
         TOTDEPall  = 0d0
         k_ox       = 0d0

         OLDMLD     = MLDav(I,J)
         MLDav(I,J) = MLDav(I,J) + dMLD(I,J) * DTSRCE
         MLDcm      = MLDav(I,J)
!--jaf.start
         ! Add error trap to prevent new MLD from being negative
         ! (jaf, 7/6/11)
         IF (MLDcm .LT. 0d0) MLDcm = 0d0
!--jaf.end

         ! Get fractions of land and ocean in the grid box [unitless]
         ! Use fractional land type information in MERRA. Also make sure
         ! we do not use boxes that are mostly sea ice for consistency
         ! FROCEAN is a constant, so to get correct ocean fraction we
         ! need to subtract the sea ice fraction. Don't let the fraction
         ! be less than zero (jaf, 4/26/11)
#if   defined( MERRA ) || defined( GEOS_FP )
         FRAC_L       = State_Met%FRLAND(I,J)
         FRAC_O       = MAX( State_Met%FROCEAN(I,J) -
     &                       State_Met%FRSEAICE(I,J), 0d0 )
         IS_OCEAN_BOX = ( ( FRAC_O > 0d0                     ) .and. 
     &                    ( State_Met%SEAICE00(I,J)  > 0.5d0 ) )
#else
         FRAC_L       = State_Met%FRCLND(I,J)
         FRAC_O       = 1d0 - FRAC_L
         IS_OCEAN_BOX = ( IS_WATER( I, J, State_Met ) )
#endif

!--jaf.end

         ! Change ocean mass due to mixed layer depth change
         ! Keep before next IF so that we adjust mass in ice-covered boxes
         CALL MLD_ADJUSTMENT( I, J, OLDMLD*1d-2, MLDcm*1d-2, 
     &                        Input_Opt, State_Met )

         !===========================================================
         ! Make sure we are in an ocean box
         !===========================================================
         ! Update to Holmes et al. 2010 version (jaf, 4/11/11)
!         IF ( ( ALBD(I,J) <= 0.4d0 ) .and. 
!     &        ( FRAC_L    <  0.8d0 ) .and.
!     &        ( MLDCM     > 0.99d0 )      ) THEN
!--jaf.start
          ! Use consistent criteria for Ocean/Land/Ice categories
          ! with snowpack and terrestrial emissions  !CDH 5/18/2010
          !IF ( ( IS_WATER(I,J) ) .AND. ( MLDCM > 0.99d0 ) ) THEN
         IF ( (IS_OCEAN_BOX) .and. (MLDCM > 0.99d0) ) THEN
!--jaf.end

            !===========================================================
            ! Reduction and oxidation coefficients
            !===========================================================    
            ! Avoid having NPP or CHL to be zero
            NPP(I,J) = MAX ( NPP(I,J) , NPPMINNUM )  

            CHL(I,J)    = MAX ( CHL(I,J) , CHLMINNUM )                       

            ! Light attenuation (RADz) is calculated
            EC     = (EC_w + ( EC_doc * C_doc * ( NPP(I,J) / NPP_avg ) ) 
     &               + ( ECchla * CHL(I,J) / 1000 ) )             

            RADz   = ( 1 / ( MLDcm * 1d-2 ))
     &               * (State_Met%RADSWG(I,J) / EC )    
     &               * ( 1 - EXP( -EC * ( MLDcm * 1d-2) ) )

            !--------------------------------------------------------
            ! Hg(tot)aq reduction rate constants
            !--------------------------------------------------------

            k_red_rad   = ( k_radbase * RADz )    

            k_red_bio   = ( ( k_biobase * NPP(I,J) ) * 1.1 )   !NPP is increased by 0.1

            k_red       = k_red_rad + k_red_bio


            !-------------------------------------------------------
            ! Hg(0)aq oxidation rate constants
            !------------------------------------------------------
            
            k_ox        = ( k_ox_dark + ( k_oxbase * RADz ) )  


            !=========================================================
            ! Partitioning and organic carbon
            !========================================================= 

            ! Calculation of C_tot for stratified waters
            IF (CHL(I,J) <= 1.0) THEN
               C_tot    = 36.1d0 * CHL(I,J)**0.357d0
            ELSE
               C_tot    = 37.7d0 * CHL(I,J)**0.615d0
            ENDIF

            ! Calculation of the euphotic depth
            IF (C_tot > 13.65) THEN
               Ze       = 912.0d0 * C_tot**(-0.839d0)
            ELSE
               Ze       = 426.3d0 * C_tot**(-0.547d0)
            ENDIF

            ! Recalculation of C_tot if water is shown to be well-mixed
            IF ((Ze/(MLDcm*1d-2)) < 1) THEN
               C_tot    = 42.1d0 * CHL(I,J)**0.538d0
            ENDIF

            !--------------------------------------------------------------
            ! Standing stock of organic carbon and total biomass
            !--------------------------------------------------------------
            ! Calculated based on C:Chl ratio of 80 (wetzel et al 2006)
            ! Stodk of organic carbon is in mgC/m2
            ! Then converting to OC_tot_kg in kg/grid

            OC_tot      = C_tot * 80.0d0

            OC_tot_kg   = OC_tot * 1d-6 * A_M2 * FRAC_O


            ! Total biomas is a proxy for SPM (mg/m3) used in Hg(II)
            ! partitioning. Calculated by multiplying the standing 
            ! stock of organic carbon with 10 (exp Bundy 2004)

            SPM = ( OC_tot * 10.0d0 / ( MLDcm * 1d-2 ) ) * 1.1                   

            !-------------------------------------------------------------- 
            ! Hg(II) - Hg(P) partitioning coefficient
            !--------------------------------------------------------------
            ! Kd_part is based on Mason et al. 1998 and Mason &
            ! Fitzgerald 1993. (L/kg)
            ! SPM is converted to kg/L by 10E-9

            ! SPM = Suspended particulate matter (kg/L)

            Frac_Hg2    = 1 / ( 1 + Kd_part * SPM * 1d-9)


            !--------------------------------------------------------------
            ! Sea surface temperature in both [K] and [C]
            !--------------------------------------------------------------
            ! where TSKIN is the temperature (K) at the ground/sea surface
            ! (Use as surrogate for SST, cap at freezing point)

            TK     = MAX( State_Met%TSKIN(I,J), 273.15d0 )

            TC     = TK - 273.15d0

            !==============================================================
            ! Volatilisation of Hg0
            !==============================================================
            
            ! Henry's law constant (gas->liquid) [unitless] [L water/L air]  
            ! (ref: Andersson et al. 2008)

            H      = EXP( ( -2404.3d0 / TK ) + 6.92d0 )

            ! Viscosity as a function of changing temperatures
            ! (ref: Loux 2001)
            ! The paper says the viscosity is given in cP but us really P
            ! and we therefor multiply with 100 to get cP.

            vi    = ( 10**( ( 1301.0d0 / ( 998.333d0 + 8.1855d0 
     &              * ( TC - 20.0d0 )+ 0.00585d0 * (TC - 20.0d0 )**2 ) ) 
     &              - 3.30233d0 ) ) * 100.0d0      

            ! Schmidt # for Hg [unitless]                             
            ! Sc = v/D = kinematic viscosity/diffusivity
            ! (ref: Poissant et al 2000; Wilke and Chang 1995)
            ! to correct for seawater D0 is decreased by 6% as suggested
            ! by Wanninkhof (1992)

            D = 7.4D-8 * sqrt( 2.26 * 18.0 ) * TK /
     &             ( ( 14.8**0.6 ) *vi )

            Sc   = ( 0.017d0 * EXP( -0.025d0 * TC ) ) / D                      
            
            ! Schmidt # of CO2 [unitless] for CO2 in seawater at 20 degrees C
            ! The value is set to a constant based on other ocean studies
            ! (Gardfeld et al. 2003, Rolfhus & Fitzgerald 2004, Mason et al. 2001)

            ! Correction of the Schmidt # with temperature based on Poissant
            ! et al. (2000) (for freshwatersystems).

            ScCO2  = 644.7d0 + TC * ( -6.16d0 + TC * ( 0.11d0 ) ) 

            ! Square of surface (actually 10m) wind speed [m2/s2]
            Usq    = SFCWINDSQR( I, J, State_Met%U10M, State_Met%V10M )

            !------------------------------------------------------
            ! Parameterizations for calculating water side mass trasfer coefficient 
            !------------------------------------------------------
            ! Mass transfer coefficient [cm/h], from Nightingale et al. 2000
            Kw     = ( 0.25d0 * Usq ) / SQRT( Sc / ScCO2 )             

            !-----------------------------------------------------
            ! Additional parameterizations:

            ! Nightinale et al. 2000 for instantanous winds

!            Kw     = ( 0.33d0*SQRT(usq)+0.22d0*Usq) / SQRT( Sc / ScCO2 )    

            ! Lis and Merlivat 1986
            ! Has less emphasis on windspeed as a driver for evasion
            ! Gives a less total evasion than the Nigthingale et al. 2000

!            IF (SQRT(Usq) <= 3.6d0 ) THEN                               
!               Kw = ( 0.17d0 * SQRT(Usq) * ( Sc / ScCO2 )**0.67d0 )
!            ELSE IF (SQRT(Usq) > 3.6d0 .and. SQRT(Usq) <= 13d0 ) THEN
!               Kw = ( ( 2.8d0 * SQRT(Usq))-9.6 ) * ( Sc / ScCO2 )**0.5d0
!            ELSE  
!               Kw = ( ( 5.9d0 * SQRT(Usq))-49.3 ) * ( Sc / ScCO2)**0.5d0
!            ENDIF

            ! Wanninkhof et al (1992)

!            Kw     = ( 0.31d0 * Usq ) / SQRT( Sc / ScCO2 )

            !===========================================================
            ! Particulate sinking                              
            !===========================================================
            ! HgP sinking is based on Sunderland & Mason 2007.
            ! JorgC originally in gC m-2 year-1, which is convereted 
            ! to kgC grid-1 timestep-1
            ! NPP is converted from mgC/m2/d-1 to gC/m2/year-1
            ! JorgC = 0.1 ( ( NPP * 12 )**1.77 ) *  MLD**n * M2 * Frac_O 
            !         * 10^-3 * DTSRCE / ( 365 * 24 * 60 * 60 )

            JorgC_kg  = ( ( 0.1d0 * ((( NPP(I,J) * 365) / 1000 )**1.77 )
     &                * ( ( MLDcm * 1d-2 )**(-0.74d0) ) * A_M2 * FRAC_O 
     &                * 1d-3) / ( 365.0d0 * 24.0d0 * 60.0d0 * 60.0d0 ) )
     &                * DTSRCE                     


            !-----------------------------------------------------------
            ! Physical transport for tracers, Part II:
            ! Upward current transport (Ekman pumping)
            ! Upward mass flux is:
            ! Mass = (Vol upwelling water) * (Conc. below thermocline)
            ! Mass = (VEL * AREA * TIME  ) * (C * Molar Mass )
            !-----------------------------------------------------------
 
            ! Use CDEEPATL to scale deepwater in NAtlantic            

            IF ( UPVEL(I,J) > 0d0 ) THEN
                 
            ! Loop over total Hg (and ocean Hg if necessary)
            DO C = 1, N_tot_oc

               ! Grid-box latitude [degrees]
               Y = GET_YMID( I, J, 1 )
         
               ! Grid box longitude [degrees]
               X = GET_XMID( I, J, 1 )

               ! Get Hg category #
               IF ( C == 1 ) NN = ID_Hg_tot

               !- eds 8/31/10 ----------------------------------------------
               ! On second loop, now assign regional tags individually
               !IF ( C == 2 ) NN = ID_Hg_oc
               !------------------------------------------------------------


               ! Atlantic
               IF ( ( X >= -80.0 .and. X < 25.0 )  .and.
     &              ( Y >= -25.0 .and. Y < 55.0 ) ) THEN    !(anls,100114)

                  !- eds 8/31/10 ------------------------------------------                  
                  IF ( C == 2 ) NN = ID_Hg_atl
                  !--------------------------------------------------------

                  ! Hg0 (kg)
                  Hg0aq(I,J,NN) = Hg0aq(I,J,NN) + UPVEL(I,J)
     &                 * ( MHg * A_M2 * FRAC_O * DTSRCE * CDeepatl(1) )

                  ! Hg2 
                  Hg2aq(I,J,NN) = Hg2aq(I,J,NN) + UPVEL(I,J)
     &                 * ( MHg * A_M2 * FRAC_O * DTSRCE * CDeepatl(2) )

                  ! Hg particulate
!                  IF ( C == 1 ) THEN
!                     HgC(I,J)   = HgC(I,J) + UPVEL(I,J) 
!     &                 * ( MHg * A_M2 * FRAC_O * DTSRCE * CDeepatl(3) )
                     HgPaq(I,J,NN)   = HgPaq(I,J,NN) + UPVEL(I,J) 
     &                 * ( MHg * A_M2 * FRAC_O * DTSRCE * CDeepatl(3) )
!                  ENDIF
 
               !North Pacific (west)  
               ELSE IF ( ( X >= -180.0 .and. X < -80.0 )  .and.
     &                   ( Y >=   30.0 .and. Y <  70.0 ) ) THEN

       	       	  !- eds 8/31/10 ------------------------------------------
                  IF ( C == 2 ) NN = ID_Hg_npa
       	       	  !--------------------------------------------------------

                  ! Hg0 (kg)
                  Hg0aq(I,J,NN) = Hg0aq(I,J,NN) + UPVEL(I,J)
     &                 * ( MHg * A_M2 * FRAC_O * DTSRCE * CDeepnpa(1) )

                  ! Hg2 
                  Hg2aq(I,J,NN) = Hg2aq(I,J,NN) + UPVEL(I,J)
     &                 * ( MHg * A_M2 * FRAC_O * DTSRCE * CDeepnpa(2) )

                  ! Hg particulate
!                  IF ( C == 1 ) THEN
                     HgPaq(I,J,NN)   = HgPaq(I,J,NN) + UPVEL(I,J) 
     &                 * ( MHg * A_M2 * FRAC_O * DTSRCE * CDeepnpa(3) )
!                  ENDIF

               !North Pacific (east)  
               ELSE IF ( ( X >= 25.0 .and. X < 180.0 )  .and.
     &                   ( Y >= 30.0 .and. Y <  70.0 ) ) THEN

       	       	  !- eds 8/31/10 ------------------------------------------
                  IF ( C == 2 ) NN = ID_Hg_npa
       	       	  !--------------------------------------------------------

                  ! Hg0 (kg)
                  Hg0aq(I,J,NN) = Hg0aq(I,J,NN) + UPVEL(I,J)
     &                 * ( MHg * A_M2 * FRAC_O * DTSRCE * CDeepnpa(1) )

                  ! Hg2 
                  Hg2aq(I,J,NN) = Hg2aq(I,J,NN) + UPVEL(I,J)
     &                 * ( MHg * A_M2 * FRAC_O * DTSRCE * CDeepnpa(2) )

                  ! Hg particulate
!                  IF ( C == 1 ) THEN
                     HgPaq(I,J,NN)   = HgPaq(I,J,NN) + UPVEL(I,J) 
     &                 * ( MHg * A_M2 * FRAC_O * DTSRCE * CDeepnpa(3) )
!                  ENDIF


               ! North Atlantic
               ELSE IF ( ( X >= -80.0 .and. X < 25.0 )  .and.
     &                   ( Y >=  55.0 .and. Y < 70.0 ) ) THEN

       	       	  !- eds 8/31/10 ------------------------------------------
                  IF ( C == 2 ) NN = ID_Hg_nat
       	       	  !--------------------------------------------------------

                  ! Hg0 (kg)
                  Hg0aq(I,J,NN) = Hg0aq(I,J,NN) + UPVEL(I,J)
     &                 * ( MHg * A_M2 * FRAC_O * DTSRCE * CDeepnat(1) )

                  ! Hg2 
                  Hg2aq(I,J,NN) = Hg2aq(I,J,NN) + UPVEL(I,J)
     &                 * ( MHg * A_M2 * FRAC_O * DTSRCE * CDeepnat(2) )

                  ! Hg particulate
!                  IF ( C == 1 ) THEN
                     HgPaq(I,J,NN)   = HgPaq(I,J,NN) + UPVEL(I,J) 
     &                 * ( MHg * A_M2 * FRAC_O * DTSRCE * CDeepnat(3) )
!                  ENDIF

               
               ! South Atlantic
               ELSE IF ( ( X >= -80.0 .and. X <  25.0 )  .and.
     &                   ( Y >= -65.0 .and. Y < -25.0 ) ) THEN   !(anls,100114)

       	       	  !- eds 8/31/10 ------------------------------------------
                  IF ( C == 2 ) NN = ID_Hg_sat
       	       	  !--------------------------------------------------------

                  ! Hg0 (kg)
                  Hg0aq(I,J,NN) = Hg0aq(I,J,NN) + UPVEL(I,J)
     &                 * ( MHg * A_M2 * FRAC_O * DTSRCE * CDeepsat(1) )

                  ! Hg2 
                  Hg2aq(I,J,NN) = Hg2aq(I,J,NN) + UPVEL(I,J)
     &                 * ( MHg * A_M2 * FRAC_O * DTSRCE * CDeepsat(2) )

                  ! Hg particulate
!                  IF ( C == 1 ) THEN
                     HgPaq(I,J,NN)   = HgPaq(I,J,NN) + UPVEL(I,J) 
     &                 * ( MHg * A_M2 * FRAC_O * DTSRCE * CDeepsat(3) )
!                  ENDIF


               ! Antarctic
               ELSE IF ( Y >=  -90.0 .and. Y <  -65.0 ) THEN

       	       	  !- eds 8/31/10 ------------------------------------------
                  IF ( C == 2 ) NN = ID_Hg_ant
       	       	  !--------------------------------------------------------

                  ! Hg0 (kg)
                  Hg0aq(I,J,NN) = Hg0aq(I,J,NN) + UPVEL(I,J)
     &                 * ( MHg * A_M2 * FRAC_O * DTSRCE * CDeepant(1) )

                  ! Hg2 
                  Hg2aq(I,J,NN) = Hg2aq(I,J,NN) + UPVEL(I,J)
     &                 * ( MHg * A_M2 * FRAC_O * DTSRCE * CDeepant(2) )

                  ! Hg particulate
!                  IF ( C == 1 ) THEN
                     HgPaq(I,J,NN)   = HgPaq(I,J,NN) + UPVEL(I,J) 
     &                 * ( MHg * A_M2 * FRAC_O * DTSRCE * CDeepant(3) )
!                  ENDIF


               ! Arctic
               ELSE IF ( Y >=  70.0 .and. Y <  90.0 ) THEN

       	       	  !- eds 8/31/10 ------------------------------------------
                  IF ( C == 2 ) NN = ID_Hg_arc
       	       	  !--------------------------------------------------------

                  ! Hg0 (kg)
                  Hg0aq(I,J,NN) = Hg0aq(I,J,NN) + UPVEL(I,J)
     &                 * ( MHg * A_M2 * FRAC_O * DTSRCE * CDeeparc(1) )

                  ! Hg2 
                  Hg2aq(I,J,NN) = Hg2aq(I,J,NN) + UPVEL(I,J)
     &                 * ( MHg * A_M2 * FRAC_O * DTSRCE * CDeeparc(2) )

                  ! Hg particulate
!                  IF ( C == 1 ) THEN
                     HgPaq(I,J,NN)   = HgPaq(I,J,NN) + UPVEL(I,J) 
     &                 * ( MHg * A_M2 * FRAC_O * DTSRCE * CDeeparc(3) )
!                  ENDIF
               
               ELSE

       	       	  !- eds 8/31/10 ------------------------------------------
                  IF ( C == 2 ) NN = ID_Hg_ocn
       	       	  !--------------------------------------------------------

                  ! Hg0 (kg)
                  Hg0aq(I,J,NN) = Hg0aq(I,J,NN) + UPVEL(I,J)
     &                 * ( MHg * A_M2 * FRAC_O * DTSRCE * CDeep(1) )

                  ! Hg2 
                  Hg2aq(I,J,NN) = Hg2aq(I,J,NN) + UPVEL(I,J)
     &                 * ( MHg * A_M2 * FRAC_O * DTSRCE * CDeep(2) )

                  ! Hg particulate
!                  IF ( C == 1 ) THEN
                     HgPaq(I,J,NN)   = HgPaq(I,J,NN) + UPVEL(I,J) 
     &                 * ( MHg * A_M2 * FRAC_O * DTSRCE * CDeep(3) )
!                  ENDIF
               
               ENDIF

            ENDDO
 
                  
            !----------------------------------------------------------
            ! Physical transport for TOTAL TRACERS, Part III:
            ! Downward current transport (Ekman pumping)
            ! Treated as a deposition velocity
            ! d(Mass)/dt = - VEL * Mass / BoxHeight
            !----------------------------------------------------------
            ELSE  

               ! Loop over all types of tagged tracers
               DO NN = 1, N_Hg_CATS

                  ! Hg0
                  Hg0aq(I,J,NN) = Hg0aq(I,J,NN) 
     &                * ( 1d0 + UPVEL(I,J) * DTSRCE / ( MLDcm * 1d-2 ) ) 
                  
                  ! Hg2
                  Hg2aq(I,J,NN) = Hg2aq(I,J,NN) 
     &                * ( 1d0 + UPVEL(I,J) * DTSRCE / ( MLDcm * 1d-2 ) )
               
                  ! Hg particulate
                     HgPaq(I,J,NN)  = HgPaq(I,J,NN) 
     &                * ( 1d0 + UPVEL(I,J) * DTSRCE / ( MLDcm * 1d-2 ) )

               ENDDO

            ENDIF


            !===========================================================
            ! Calculate reduction, conversion, sinking, evasion
            !
            ! (1) Hg2 <-> HgP and HgP sinks
            ! (2) Hg2 <-> Hg0 and Hg0 evades
            !
            ! NOTE: N is the GEOS-CHEM tracer # (for STT)
            !       and NN is the Hg category # (for Hg0aq, Hg2aq, HgP)
            !===========================================================

            ! Loop over all Hg categories
            DO NN = 1, N_Hg_CATS

               ! Reset flux each timestep
               FLUX(I,J,NN)  = 0d0 
               FUP(I,J,NN)   = 0d0
               FDOWN(I,J,NN) = 0d0

                  
               !--------------------------------------------------------
               ! Calculate new Hg(II) mass
               !--------------------------------------------------------

               ! Before 11/3/2009 (cdh, hamos)
               !! Total Hg(II) deposited on ocean surface [kg]
               !TOTDEP = (WD_Hg2(I,J,NN) + DD_Hg2(I,J,NN))*FRAC_O 
               !                 
               ! Total Hg(II) deposited on ocean surface [kg]
               ! Includes gaseous and particulate reactive Hg(II)
               ! plus anthropogenic primary Hg(p) (cdh, hamos 11/3/2009)
               TOTDEPall = (WD_Hg2(I,J,NN) + DD_Hg2(I,J,NN) +
     &                   WD_HgP(I,J,NN) + DD_HgP(I,J,NN) ) 

               TOTDEP        = TOTDEPall * FRAC_O

               ! Add deposited Hg(II) to the Hg(II)tot ocean mass [kg]
               Hg2aq_tot     = Hg2aq(I,J,NN) + HgPaq(I,J,NN) + TOTDEP      

               Hg2aq(I,J,NN) = Hg2aq_tot * Frac_Hg2               


               ! Mass of Hg(II)  -->  Hg(0) 
               ! Only a certain percentage of Hg(II) is considered reducible
               Hg2_RED       = Hg2aq(I,J,NN) * 0.4d0 * k_red * DTSRCE       
               ! Mass of Hg(0) --> Hg(II)
               Hg0_OX        = Hg0aq(I,J,NN) * k_ox * DTSRCE

               ! Amount of Hg(II) that is lost [kg]
               Hg2_GONE      = Hg2_RED - Hg0_OX                        

               ! Cap Hg2_GONE with available Hg2
               IF ( Hg2_GONE > Hg2aq(I,J,NN) ) THEN 
                  Hg2_GONE   = MIN( Hg2_GONE, Hg2aq(I,J,NN) )
               ENDIF

               IF ( (Hg2_GONE * (-1d0)) >  Hg0aq(I,J,NN)) THEN         
                  Hg2_GONE   = (Hg0aq(I,J,NN)*(-1))   
                  !MAX (Hg2_GONE ,(Hg0aq(I,J,NN)*(-1d0)))
               ENDIF

               ! Hg(II) ocean mass after reduction and conversion [kg]
               Hg2aq(I,J,NN) = Hg2aq(I,J,NN) - Hg2_GONE

               !--------------------------------------------------------
               ! Calculate new Hg(P) mass
               !--------------------------------------------------------

               !- eds 8/31/10 ------------------------------------------
               !IF ( NN == 1 ) THEN
               !--------------------------------------------------------

                  ! HgP ocean mass after conversion
                  HgPaq(I,J,NN)   = Hg2aq_tot * ( 1 - Frac_Hg2)

                  !----------------------------------------------------
                  ! Conversion between OC and Hg                          
                  !----------------------------------------------------
                  ! Hg/C ratio based on HgP(kg) and Stock of organic C(kg)
                  ! HgPaq_sunk funtion of C sunk and HgP/C ratio   
                  HgPaq_SUNK(I,J,NN)  = JorgC_kg * ( HgPaq(I,J,NN )  
     &                                  / OC_tot_kg)

                  ! HgP ocean mass after sinking [kg]
                  HgPaq(I,J,NN)   = HgPaq(I,J,NN) - HgPaq_SUNK(I,J,NN)

                  HgPaq(I,J,NN)   = MAX ( HgPaq(I,J,NN) , 0d0 )                 

                  ! Store carbon sinking [kgC/time]
                  IF ( ND03 > 0 ) THEN
!- eds 8/31/10 --------------------------------------------------------
!                     AD03(I,J,12) = AD03(I,J,12) + JorgC_kg
                     AD03(I,J,12,1) = AD03(I,J,12,1) + JorgC_kg
!----------------------------------------------------------------------
                  ENDIF
                
               !- eds 8/31/10 -----------------------------------------
               !ENDIF
               !-------------------------------------------------------

               !--------------------------------------------------------
               ! Calculate new Hg(0) mass
               !--------------------------------------------------------

               ! Hg0 tracer number (for STT)
               N             = ID_Hg0(NN)

               ! Add converted Hg(II) and subtract converted Hg(0) mass 
               ! to the ocean mass of Hg(0) [kg]
               Hg0aq(I,J,NN) = Hg0aq(I,J,NN) + Hg2_GONE     

               !--------------------------------------------------------
               ! Calculate oceanic and gas-phase concentration of Hg(0)
               !--------------------------------------------------------
                  
               ! Concentration of Hg(0) in the ocean [ng/L]
               CHg0aq        = ( Hg0aq(I,J,NN) * 1d11   ) /   
     &                         ( A_M2          * FRAC_O ) / MLDcm 
               
               ! Gas phase Hg(0) concentration: convert [kg] -> [ng/L]
               CHg0          = STT(I,J,1,N) * 1.0D9 /
     &                         State_Met%AIRVOL(I,J,1)
               
               !--------------------------------------------------------
               ! Compute flux of Hg(0) from the ocean to the air
               !--------------------------------------------------------

               ! Compute ocean flux of Hg0 [cm/h*ng/L]
               FLUX(I,J,NN)  = Kw * ( CHg0aq - ( CHg0 / H ) )     

               ! TURN OFF EVASION
!               FLUX(I,J,NN)= MIN(0.,FLUX(I,J,NN))

               !Prior to 09 Nov 2011, H Amos ---------------------
               !Extra diagnostic: compute flux up and flux down
               !FUP(I,J,NN)   = ( Kw * CHg0aq )
               !FDOWN(I,J,NN) = ( Kw * CHg0 / H )     
               !--------------------------------------------------


               ! Convert [cm/h*ng/L] --> [kg/m2/s] --> [kg/s]
               ! Also account for ocean fraction of grid box
               FLUX(I,J,NN)  = FLUX(I,J,NN) * TO_KGM2S * A_M2 * FRAC_O 

               !Prior to 09 Nov 2011, H Amos ---------------------             
               !FUP(I,J,NN)  = FUP(I,J,NN) * TO_KGM2S * A_M2 * FRAC_O
               !FDOWN(I,J,NN)  = FDOWN(I,J,NN) * TO_KGM2S * A_M2 * FRAC_O
               !--------------------------------------------------

               !--------------------------------------------------------
               ! Flux limited by ocean and atm Hg(0)
               !--------------------------------------------------------

               !Prior to 09 Nov 2011, H Amos --------------------- 
               ! Cap the flux w/ the available Hg(0) ocean mass
               !IF ( FLUX(I,J,NN) * DTSRCE > Hg0aq(I,J,NN) ) THEN 
               !   FLUX(I,J,NN) = Hg0aq(I,J,NN) / DTSRCE 
               !   FUP(I,J,NN)  = FLUX(I,J,NN)-FDOWN(I,J,NN)
               !ENDIF
               IF ( FLUX(I,J,NN) * DTSRCE > Hg0aq(I,J,NN) ) THEN 
                  FLUX(I,J,NN) = Hg0aq(I,J,NN) / DTSRCE 
               ENDIF
               !--------------------------------------------------

                
               ! Cap the neg flux w/ the available Hg(0) atm mass
               IF ( (-FLUX(I,J,NN) * DTSRCE ) > STT(I,J,1,N) ) THEN
                  FLUX(I,J,NN) = -STT(I,J,1,N) / DTSRCE       
               ENDIF
                
               
               ! Cap FDOWN with available Hg(0) atm mass

!               IF ((FDOWN(I,J,NN)*DTSRCE)>STT(I,J,1,N)) THEN
!                  FDOWN(I,J,NN) = STT(I,J,1,N) / DTSRCE
!               ENDIF

               ! make sure Fup and Fdown do not underflow either
               ! debug 2x2.5 diagnostic?
               FUP(I,J,NN) = MAX (FUP(I,J,NN), SMALLNUM )
               FDOWN(I,J,NN) = MAX (FDOWN(I,J,NN),SMALLNUM )
               !--------------------------------------------------------
               ! Remove amt of Hg(0) that is leaving the ocean [kg]
               !--------------------------------------------------------
               Hg0aq(I,J,NN) = Hg0aq(I,J,NN) - ( FLUX(I,J,NN) * DTSRCE ) 

               ! Make sure Hg0aq does not underflow (cdh, bmy, 3/28/06)
               Hg0aq(I,J,NN) = MAX( Hg0aq(I,J,NN), SMALLNUM )

!               Hgaq_tot = HgC(I,J) + Hg0aq(I,J,NN) + Hg2aq(I,J,NN)
               Hgaq_tot(I,J,NN) = HgPaq(I,J,NN) + Hg0aq(I,J,NN) 
     &                            + Hg2aq(I,J,NN) !eds 5/15/12

            ENDDO   

            !-----------------------------------------------------------
            ! ND03 diagnostics ("OCEAN-HG")
            !-----------------------------------------------------------

            IF ( ND03 > 0 ) THEN

              ! eds 9/9/10 added NN loop
              DO NN = 1, N_HG_CATS

               ! Aqueous Hg(0) mass [kg]
!               AD03(I,J,2)  = AD03(I,J,2)  + Hg0aq(I,J,ID_Hg_tot) 
               AD03(I,J,2,NN)  = AD03(I,J,2,NN)  + Hg0aq(I,J,NN)

               ! Aqueous Hg(II) mass [kg] 
!               AD03(I,J,7)  = AD03(I,J,7)  + Hg2aq(I,J,ID_Hg_tot) 
               AD03(I,J,7,NN)  = AD03(I,J,7,NN)  + Hg2aq(I,J,NN)

               ! Hg2 sunk deep into the ocean [kg/time]
!               AD03(I,J,8)  = AD03(I,J,8)  + HgPaq_SUNK
               AD03(I,J,8,NN)= AD03(I,J,8,NN) + HgPaq_SUNK(I,J,NN)

               ! HgTot aqua mass [kg] 
!               AD03(I,J,10) =AD03(I,J,10) + Hgaq_tot        
               AD03(I,J,10,NN) =AD03(I,J,10,NN) + Hgaq_tot(I,J,NN)
 
               ! HgP ocean mass [kg]
!               AD03(I,J,11) = AD03(I,J,11) + HgPaq(I,J) 
               AD03(I,J,11,NN) = AD03(I,J,11,NN) + HgPaq(I,J,NN)

               !Prior to 09 Nov 2011, H Amos --------------------------
               ! flux up and down (eck)
               !AD03(I,J,16) = AD03(I,J,16) + FUP(I,J,ID_Hg_tot)*DTSRCE
               !AD03(I,J,17) = AD03(I,J,17) + FDOWN(I,J,ID_Hg_tot)*DTSRCE
!              IF (FLUX(I,J,ID_Hg_tot) > 0d0) THEN
!                 AD03(I,J,16) = AD03(I,J,16) 
!    &                           + FLUX(I,J,ID_Hg_tot) * DTSRCE
!              ELSE IF (FLUX(I,J,ID_Hg_tot) < 0d0) THEN
!                 AD03(I,J,17) = AD03(I,J,17) 
!    &                           + ( abs(FLUX(I,J,ID_Hg_tot) * DTSRCE ))
!              ENDIF
               ! Loop over tagged tracers (eds)
               IF (FLUX(I,J,NN) > 0d0) THEN
                  AD03(I,J,16,NN) = AD03(I,J,16,NN) 
     &                           + FLUX(I,J,NN) * DTSRCE
               ELSE IF (FLUX(I,J,NN) < 0d0) THEN
                  AD03(I,J,17,NN) = AD03(I,J,17,NN) 
     &                           + ( abs(FLUX(I,J,NN) * DTSRCE ))
               ENDIF
               !-------------------------------------------------------

              ENDDO

            ENDIF

           
         !==============================================================
         ! If we are not in an ocean box, set Hg(0) flux to zero
         !==============================================================
         ELSE

            DO NN = 1, N_Hg_CATS 
               FLUX(I,J,NN) = 0d0
               FUP(I,J,NN)=0d0
               FDOWN(I,J,NN)=0d0
            ENDDO               

         ENDIF 
      
         !==============================================================
         ! Zero amts of deposited Hg2 for next timestep [kg]  
         !==============================================================
         
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! Free pointer
      NULLIFY( STT )

      !=================================================================
      ! Check tagged & total sums (if necessary)
      !=================================================================
      IF ( USE_CHECKS .and. LSPLIT ) THEN
         CALL CHECK_ATMOS_MERCURY(  State_Chm, 
     &                              'end of OCEAN_MERCURY_FLUX' )
         CALL CHECK_OCEAN_MERCURY(  'end of OCEAN_MERCURY_FLUX' )
         CALL CHECK_OCEAN_FLUXES (  'end of OCEAN_MERCURY_FLUX' )
         CALL CHECK_FLUX_OUT( FLUX, 'end of OCEAN_MERCURY_FLUX' )
      ENDIF

      ! Return to calling program
      END SUBROUTINE OCEAN_MERCURY_FLUX

!------------------------------------------------------------------------------

!      SUBROUTINE OCEAN_MERCURY_READ( THISMONTH )
      SUBROUTINE OCEAN_MERCURY_READ( THISMONTH, THISYEAR )
!
!******************************************************************************
!  Subroutine OCEAN_MERCURY_READ reads in the mixed layer depth, net primary 
!  productivity, upwelling and radiation climatology for each month.  
!  This is needed for the ocean flux computation. 
!  (sas, cdh, bmy, 1/20/05, 3/28/06)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) THISMONTH (INTEGER) : Month to read fields (1-12)
!
!  NOTES:
!  (1 ) Modified for S. Strode's latest ocean Hg code.  Now read files
!        from DATA_DIR_1x1/mercury_200511. (sas, cdh, bmy, 3/28/06)
!******************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD
      USE DIRECTORY_MOD, ONLY : DATA_DIR_1x1, DATA_DIR
      USE TRANSFER_MOD,  ONLY : TRANSFER_2D
      USE TIME_MOD,	 ONLY : ITS_A_LEAPYEAR ! jaf
      USE FILE_MOD,	 ONLY : IOERROR, IU_FILE ! jaf

      USE CMN_SIZE_MOD      ! Size parameters

      ! Arguments
      INTEGER, INTENT(IN)    :: THISMONTH
      INTEGER, INTENT(IN)    :: THISYEAR ! jaf

      ! Local Variables
      LOGICAL, SAVE          :: FIRST = .TRUE.
      REAL*4                 :: ARRAY(IIPAR,JJPAR,1)
      REAL*8                 :: TAU
      CHARACTER(LEN=255)     :: FILENAME
      ! Add time variables to get correct number of seconds per month
      ! and prevent negative flow later (jaf, 12/8/11)
      INTEGER	   :: DAYS_IN_MONTH, LASTMONTH
      INTEGER	   :: M(12) = (/ 31, 28, 31, 30, 31, 30,
     &                           31, 31, 30, 31, 30, 31 /)
     
      !=================================================================
      ! OCEAN_MERCURY_READ begins here!
      !=================================================================
     
      ! Get the last month (jaf, 12/8/11)
      IF ( THISMONTH == 1 ) THEN
         LASTMONTH = 12
      ELSE
         LASTMONTH = THISMONTH - 1
      ENDIF

      ! Calculate days in this month (jaf, 12/8/11)
      IF ( LASTMONTH == 2 .and. ITS_A_LEAPYEAR( THISYEAR ) ) THEN
         DAYS_IN_MONTH = M(LASTMONTH) + 1
      ELSE
         DAYS_IN_MONTH = M(LASTMONTH)
      ENDIF

      !------------------------------
      ! Mixed layer depth [cm]
      !------------------------------

      ! MLD file name
#if defined( GRID05x0666 ) && defined( NESTED_NA )
      ! Need to read emission file at nested-resolution, Y. Zhang 11/1/2011
      FILENAME = TRIM( DATA_DIR )                  // 
     &           'mercury_201203/MLD_DReqDT.'      // 
     &           GET_NAME_EXT()    // '.'          //
     &           GET_RES_EXT()
#else
      FILENAME = TRIM( DATA_DIR )                  // 
     &           'mercury_201007/MLD_DReqDT.geos.' // 
     &           GET_RES_EXT()
#endif

      ! Echo info
      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( '     - OCEAN_MERCURY_READ: Reading ', a )  

      ! TAU0 value (uses year 2003) !(uses year 1985, anls)
      TAU = GET_TAU0( THISMONTH, 1, 1985 )

      ! Read from disk; original units are [m]
      CALL READ_BPCH2( FILENAME, 'BXHGHT-$',    5,  
     &                 TAU,       IIPAR,        JJPAR,      
     &                 1,         ARRAY(:,:,1), QUIET=.TRUE. )

      ! Resize and cast to REAL*8
      CALL TRANSFER_2D( ARRAY(:,:,1), MLD )

      ! Convert [m] to [cm]
      MLD = MLD * 100d0

      ! First-time only: Set MDLav [cm] to MLD of first month
      IF ( FIRST ) THEN   
!-eds 10/19/10 causes differences dep. on restart date -----------
!         MLDav = MLD   
!         dMLD  = 0.0
!         MLDav = MLDprev
!----------------------------------------------------------------
         FIRST = .FALSE.
      ENDIF

!-------------------------------------------------      (anls, 090520) chl
! Chl from Modis [mg/m3]                                only for 4x5 for now  
!-------------------------------------------------

      ! Chl file name
#if defined( GRID05x0666 ) && defined( NESTED_NA )
      FILENAME = TRIM( DATA_DIR )                  // 
     &           'mercury_201203/Chl_2003.'        // 
     &           GET_NAME_EXT()     // '.'         //
     &           GET_RES_EXT()
#else
      ! Need to read emission file at nested-resolution, Y. Zhang 11/1/2011
      FILENAME = TRIM( DATA_DIR )                  // 
     &           'mercury_201007/Chl_2003.geos.'   //
     &           GET_RES_EXT()
#endif

      ! Echo info
      WRITE( 6, 100 ) TRIM( FILENAME )

      ! TAU0 values (uses year 2003)
      TAU = GET_TAU0( THISMONTH, 1, 2003 )      
#if defined (GRID2x25)
      ! Read data
      CALL READ_BPCH2( FILENAME, 'CHLO-A-$',    0,  
     &                 TAU,       IIPAR,        JJPAR,      
     &                 1,         ARRAY(:,:,1), QUIET=.FALSE. )
#else 

      CALL READ_BPCH2( FILENAME, 'CHLO-A-$',    1,
     &                 TAU,       IIPAR,        JJPAR,
     &                 1,         ARRAY(:,:,1), QUIET=.FALSE. )
#endif

      ! Resize and cast to REAL*8
      CALL TRANSFER_2D( ARRAY(:,:,1), CHL )


      !--------------------------------
      ! Net primary productivity [mg/m2/day]
      !--------------------------------
 
      ! NPP file name (anls, 100111)
#if defined( GRID05x0666 ) && defined( NESTED_NA )
      ! Need to read emission file at nested-resolution, Y. Zhang 11/1/2011
      FILENAME = TRIM( DATA_DIR )                  // 
     &           'mercury_201203/NPP_2003.'        // 
     &          GET_NAME_EXT()     // '.'          //
     &          GET_RES_EXT()  
#else
      FILENAME = TRIM( DATA_DIR )                  // 
     &           'mercury_201007/NPP_2003.geos.'   //
     &           GET_RES_EXT()
#endif

      ! Echo info
      WRITE( 6, 100 ) TRIM( FILENAME )

      ! TAU0 values (uses year 2003)
      TAU = GET_TAU0( THISMONTH, 1, 2003 )
 
      ! Read data
      CALL READ_BPCH2( FILENAME, 'GLOB-NPP',    1,  
     &                 TAU,       IIPAR,        JJPAR,      
     &                 1,         ARRAY(:,:,1), QUIET=.TRUE. )

      ! Resize and cast to REAL*8
      CALL TRANSFER_2D( ARRAY(:,:,1), NPP )
  

      !---------------------------------
      ! Ekman upwelling velocity [cm/s]
      !---------------------------------

      ! NPP file name
#if defined( GRID05x0666 ) && defined( NESTED_NA )
      ! Need to read emission file at nested-resolution, Y. Zhang 11/1/2011
      FILENAME = TRIM( DATA_DIR )                   // 
     &           'mercury_201203/ekman_upvel.'      // 
     &           GET_NAME_EXT()    // '.'           //
     &           GET_RES_EXT() 
#else
      FILENAME = TRIM( DATA_DIR_1x1 )               // 
     &           'mercury_200511/ekman_upvel.geos.' //
     &           GET_RES_EXT()
#endif

      ! Echo info
      WRITE( 6, 100 ) TRIM( FILENAME )

      ! TAU0 value (uses year 1985)
      TAU = GET_TAU0( THISMONTH, 1, 1985 )

      ! Read from disk; original units are [cm/s]
#if defined( GRID05x0666 ) && defined( NESTED_NA )
      CALL READ_BPCH2( FILENAME, 'BXHGHT-$',    6,  
     &                 TAU,       IIPAR,        JJPAR,      
     &                 1,         ARRAY(:,:,1), QUIET=.TRUE. )
#else
      CALL READ_BPCH2( FILENAME, 'EKMAN-V',     1,  
     &                 TAU,       IIPAR,        JJPAR,      
     &                 1,         ARRAY(:,:,1), QUIET=.TRUE. )
#endif

      ! Resize and cast to REAL*8
      CALL TRANSFER_2D( ARRAY(:,:,1), UPVEL )

      ! convert [cm/s] to [m/s]
      UPVEL = UPVEL * 1.D-2
  
      ! Return to calling program
      END SUBROUTINE OCEAN_MERCURY_READ

!------------------------------------------------------------------------------

!      SUBROUTINE GET_MLD_FOR_NEXT_MONTH( THISMONTH )
      SUBROUTINE GET_MLD_FOR_NEXT_MONTH( THISMONTH, THISYEAR )
!
!******************************************************************************
!  Subroutine GET_MLD_FOR_NEXT_MONTH reads the mixed-layer depth (MLD) 
!  values for the next month. (sas, cdh, bmy, 3/28/06)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) THISMONTH (INTEGER) : Current month number (1-12)
!  (2 ) THISYEAR (INTEGER)  : Current year 
!
!  NOTES:
!  (1 ) Now read files from DATA_DIR_1x1/mercury_200511 (bmy, 3/28/06)
!  (2 ) Now use actual number of days in given month (jaf, 7/6/11)
!******************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD,     ONLY : GET_TAU0, GET_RES_EXT, READ_BPCH2
      USE DIRECTORY_MOD, ONLY : DATA_DIR
      USE TRANSFER_MOD,  ONLY : TRANSFER_2D
      USE TIME_MOD,      ONLY : ITS_A_LEAPYEAR
      USE BPCH2_MOD,     ONLY : GET_NAME_EXT   ! Y. Zhang

      USE CMN_SIZE_MOD      ! Size parameters

      ! Arguments
      INTEGER, INTENT(IN)    :: THISMONTH
      INTEGER, INTENT(IN)    :: THISYEAR

      ! Local variables
      INTEGER                :: I, J, NEXTMONTH
      REAL*4                 :: ARRAY(IIPAR,JJPAR,1)
      REAL*8                 :: TAU
      CHARACTER(LEN=255)     :: FILENAME
      INTEGER                :: DAYS_IN_MONTH
      INTEGER                :: M(12) = (/ 31, 28, 31, 30, 31, 30,
     &                                     31, 31, 30, 31, 30, 31 /)

      !=================================================================
      ! GET_MLD_FOR_NEXT_MONTH begins here!
      !=================================================================
      
      ! Calculate days in this month (jaf, 7/6/11)
      IF ( THISMONTH == 2 .and. ITS_A_LEAPYEAR( THISYEAR ) ) THEN
         DAYS_IN_MONTH = M(THISMONTH) + 1
      ELSE
         DAYS_IN_MONTH = M(THISMONTH)
      ENDIF

      ! MLD file name
#if defined( GRID05x0666 ) && defined( NESTED_NA )
      ! Need to read emission file at nested-resolution, Y. Zhang 11/1/2011 
      FILENAME = TRIM( DATA_DIR )                  // 
     &           'mercury_201203/MLD_DReqDT.'      // 
     &           GET_NAME_EXT()    // '.'          //
     &           GET_RES_EXT()   
#else
      FILENAME = TRIM( DATA_DIR )                  // 
     &           'mercury_201007/MLD_DReqDT.geos.' //
     &            GET_RES_EXT()
#endif

      ! Echo info
      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( '     - GET_MLD_FOR_NEXT_MONTH: Reading ', a )  
      
      ! Get the next month
      NEXTMONTH = MOD( THISMONTH, 12 ) +1

      ! TAU0 value for next month (uses year 1985)
      TAU       = GET_TAU0( NEXTMONTH, 1, 1985 )

      ! Read from disk; original units are [m]
      CALL READ_BPCH2( FILENAME, 'BXHGHT-$',    5,  
     &                 TAU,       IIPAR,        JJPAR,      
     &                 1,         ARRAY(:,:,1), QUIET=.TRUE. )

      ! Resize and cast to REAL*8
      CALL TRANSFER_2D( ARRAY(:,:,1), newMLD )

      ! Convert [m] to [cm]
      newMLD = newMLD * 100d0

      ! get rate of change of MLD; convert [cm/month] -> [cm/s] 
      DO J = 1, JJPAR
      DO I = 1, IIPAR
         ! This calculation can result in negative values of MLD at the
         ! end of the the month when there are more than 30.5 days in
         ! the month. Now use actual number of days per month (jaf, 
         ! 7/6/11)
         !dMLD(I,J) = (newMLD(I,J) - MLD(I,J)) / ( 3.6d3 *24d0 * 30.5d0 )
         dMLD(I,J) = (newMLD(I,J) - MLD(I,J)) /
     &               ( 3.6d3 *24d0 * DAYS_IN_MONTH )
      ENDDO
      ENDDO

      ! Return to calling program
      END SUBROUTINE GET_MLD_FOR_NEXT_MONTH

!------------------------------------------------------------------------------

      SUBROUTINE GET_MLD_FOR_LAST_MONTH( THISMONTH )
!
!******************************************************************************
!  Subroutine GET_MLD_FOR_LAST_MONTH reads the mixed-layer depth (MLD)
!  values for the next month. Modified from GET_MLD_FOR_NEXT_MONTH (eds 10/19/10)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) THISMONTH (INTEGER) : Current month number (1-12)
!
!  NOTES:
!  (1 ) Now read files from DATA_DIR_1x1/mercury_200511 (bmy, 3/28/06)
!******************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD,     ONLY : GET_TAU0, GET_RES_EXT, READ_BPCH2
      USE BPCH2_MOD,     ONLY : GET_NAME_EXT   
      USE DIRECTORY_MOD, ONLY : DATA_DIR
      USE TRANSFER_MOD,  ONLY : TRANSFER_2D

      USE CMN_SIZE_MOD     ! Size parameters

      ! Arguments
      INTEGER, INTENT(IN)    :: THISMONTH

      ! Local variables
      INTEGER                :: I, J, LASTMONTH
      REAL*4                 :: ARRAY(IGLOB,JGLOB,1)
      REAL*8                 :: TAU
      CHARACTER(LEN=255)     :: FILENAME

      !=================================================================
      ! GET_MLD_FOR_NEXT_MONTH begins here!
      !=================================================================

      ! MLD file name
#if defined( GRID05x0666 ) && defined( NESTED_NA )
      ! Need to read emission file at nested-resolution, Y. Zhang 11/1/2011 
      FILENAME = TRIM( DATA_DIR )                  // 
     &           'mercury_201203/MLD_DReqDT.'      // 
     &           GET_NAME_EXT()    // '.'          //
     &           GET_RES_EXT()   
#else
      FILENAME = TRIM( DATA_DIR )	//
     &           'mercury_201007/MLD_DReqDT.geos.' // GET_RES_EXT()
#endif

      ! Echo info
      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( '     - GET_MLD_FOR_LAST_MONTH: Reading ', a )

      ! Get the next month
!      LASTMONTH = MOD( (THISMONTH-2), 12 ) +1
      IF (THISMONTH .EQ. 1) THEN
         LASTMONTH = 12
      ELSE 
        LASTMONTH = THISMONTH -1
      ENDIF
     
      PRINT *, 'LAST MONTH: '
      PRINT *, LASTMONTH

      ! TAU0 value for next month (uses year 1985)
      TAU	= GET_TAU0( LASTMONTH, 1, 1985 )

      ! Read from disk; original units are [m]
      CALL READ_BPCH2( FILENAME, 'BXHGHT-$',    5,
     &                 TAU,	  IGLOB,        JGLOB,
     &                 1,         ARRAY(:,:,1), QUIET=.TRUE. )

      ! Resize and cast to REAL*8
      CALL TRANSFER_2D( ARRAY(:,:,1), prevMLD )

      ! Convert [m] to [cm]
      prevMLD = prevMLD * 100d0

      ! get rate of change of MLD; convert [cm/month] -> [cm/s]
      DO J = 1, JJPAR
      DO I = 1, IIPAR
         dMLD(I,J) = (MLD(I,J) - prevMLD(I,J)) / (3.6d3 *24d0 * 30.5d0)
         ! set current MLD to intermediate b/w last month and this month
         MLDav(I,J) = prevMLD(I,J) + 0.5d0 * (MLD(I,J) - prevMLD(I,J))
      ENDDO
      ENDDO


      ! Return to calling program
      END SUBROUTINE GET_MLD_FOR_LAST_MONTH

!------------------------------------------------------------------------------

      SUBROUTINE MLD_ADJUSTMENT( I, J, MLDold, MLDnew, 
     &                           Input_Opt, State_Met )
!
!******************************************************************************
!  Subroutine MLD_ADJUSTMENT entrains new water when mixed layer depth deepens
!  and conserves concentration (leaves mass behind) when mixed layer shoals.
!  (sas, cdh, bmy, 4/18/05, 3/28/06)
!  The MLD depth is constrained so that the mean monthly concentration equals
!  the concentration in the middle of the given month. The MLD hereafter 
!  changes linearily until it reaches the middle of the next months where the 
!  process is repeted (anls, 4/30/09)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) I      (INTEGER) : GEOS-CHEM longitude index
!  (2 ) J      (INTEGER) : GEOS-CHEM latitude index
!  (3 ) MLDold (REAL*8 ) : Old ocean mixed layer depth [m]
!  (4 ) MLDnew (REAL*8 ) : New ocean mixed layer depth [m]
!
!  NOTES:
!  (1 ) Now use MERRA land fraction information (jaf, 4/26/11)
!   8 Feb 2012 - R. Yantosca - Treat GEOS-5.7.2 in the same way as MERRA
!  01 Mar 2012 - R. Yantosca - Now use GET_AREA_CM2(I,J,L) from grid_mod.F90
!  01 Mar 2012 - R. Yantosca - Now use GET_XMID(I,J,L) from grid_mod.F90
!  01 Mar 2012 - R. Yantosca - Now use GET_YMID(I,J,L) from grid_mod.F90
!  09 Nov 2012 - M. Payer    - Replaced all met field arrays with State_Met
!                              derived type object
!  25 Mar 2013 - R. Yantosca - Now use logical fields from Input_Opt
!  26 Sep 2013 - R. Yantosca - Renamed GEOS_57 Cpp switch to GEOS_FP
!******************************************************************************
!
      ! Reference to fortran90 modules
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE GRID_MOD,           ONLY : GET_AREA_M2, GET_XMID, GET_YMID
!- eds 8/31/10 -------------------------------------------------------------
!      USE TRACERID_MOD, ONLY : ID_Hg_tot, ID_Hg_oc, N_Hg_CATS
      USE TRACERID_MOD,       ONLY : ID_Hg_tot, N_Hg_CATS
      USE TRACERID_MOD,       ONLY : ID_Hg_atl,  ID_Hg_nat,  ID_Hg_sat
      USE TRACERID_MOD,       ONLY : ID_Hg_npa,  ID_Hg_arc,  ID_Hg_ant
      USE TRACERID_MOD,       ONLY : ID_Hg_ocn
!---------------------------------------------------------------------------

      USE CMN_SIZE_MOD             ! Size parameters
      
      ! Arguments
      INTEGER,        INTENT(IN)  :: I, J 
      REAL*8,         INTENT(IN)  :: MLDold, MLDnew  
      TYPE(OptInput), INTENT(IN)  :: Input_Opt   ! Input Options object
      TYPE(MetState), INTENT(IN)  :: State_Met   ! Meteorology State object

      ! Local variables
      INTEGER               :: C,    NN,     N_tot_oc    
      INTEGER               :: K, L
      REAL*8                :: A_M2, DELTAH, FRAC_O,  MHg
      REAL*8                :: X, Y                   !(added anls 01/05/09)

      ! For fields from Input_Opt
      LOGICAL               :: LSPLIT
      INTEGER               :: N_TRACERS
      REAL*8                :: TRACER_MW_KG(Input_Opt%N_TRACERS)

      !=================================================================
      ! MLD_ADJUSTMENT begins here!
      !=================================================================

      ! Copy values from Input_Opt
      LSPLIT       = Input_Opt%LSPLIT
      N_TRACERS    = Input_Opt%N_TRACERS
      TRACER_MW_KG = Input_Opt%TRACER_MW_KG(1:N_TRACERS)

      ! Loop limit for use below
      IF ( LSPLIT ) THEN
         N_tot_oc = 2
      ELSE
         N_tot_oc = 1
      ENDIF

      ! Grid box surface area [m2]
      A_M2   = GET_AREA_M2( I, J, 1 )

      ! Fraction of box that is ocean
      ! Use fractional land type information in MERRA (jaf, 4/26/11)
      ! FROCEAN is a constant, so to get correct ocean fraction we
      ! need to subtract the sea ice fraction.
#if   defined( MERRA ) || defined( GEOS_FP )
      FRAC_O     = MAX( State_Met%FROCEAN(I,J) - 
     &                  State_Met%FRSEAICE(I,J), 0d0)
#else
      FRAC_O = 1d0 - State_Met%FRCLND(I,J)
#endif

      ! Molecular weight of Hg (valid for all tagged tracers)
      MHg    = TRACER_MW_KG(ID_Hg_tot)

      ! Test if MLD increased
      IF ( MLDnew > MLDold ) THEN

         !==============================================================
         ! IF MIXED LAYER DEPTH HAS INCREASED:
         !
         ! Entrain water with a concentration specified by CDeep
         !
         ! Entrained Mass = ( Vol water entrained ) * CDeep * Molar mass
         !                = ( DELTAH * AREA * FRAC_O ) * CDeep * MHg
         !==============================================================

         ! Increase in MLD [m]
         DELTAH = MLDnew - MLDold

         ! Add Cdeepatl to North Atlantic and Cdeep to rest if the world (anls, 01/05/09)

         ! Grid-box latitude [degrees]
         Y = GET_YMID( I, J, 1 )
         
         ! Grid box longitude [degrees]
         X = GET_XMID( I, J, 1 )

         ! Loop over total Hg (and ocean Hg if necessary)
         DO C = 1, N_tot_oc

            ! Get Hg category number
            IF ( C == 1 ) NN = ID_Hg_tot

            !- eds 8/31/10-----------------------------------
            !IF ( C == 2 ) NN = ID_Hg_oc
            !------------------------------------------------

            ! Atlantic
            IF ( ( X >= -80.0 .and. X < 25.0 )  .and.
     &           ( Y >=  -25.0 .and. Y <  55.0 ) ) THEN !(anls,100114)
                        
               !- eds 8/31/10-----------------------------------
               IF ( C == 2 ) NN = ID_Hg_atl
               !------------------------------------------------

               ! Hg0
               Hg0aq(I,J,NN) = Hg0aq(I,J,NN)
     &                 + ( DELTAH * CDeepatl(1) * MHg * A_M2 * FRAC_O )

               ! Hg2
               Hg2aq(I,J,NN) = Hg2aq(I,J,NN)
     &                 + ( DELTAH * CDeepatl(2) * MHg * A_M2 * FRAC_O )

               ! HgP
!               IF ( C == 1 ) THEN
                  HgPaq(I,J,NN)   = HgPaq(I,J,NN)          
     &                 + ( DELTAH * CDeepatl(3) * MHg * A_M2 * FRAC_O )

!               ENDIF
  
            ! North Pacific (west) 
            ELSE IF ( ( X >= -180.0 .and. X < -80.0 )  .and.
     &              ( Y >=  30.0 .and. Y <  70.0 ) ) THEN
                        
               !- eds 8/31/10-----------------------------------
               IF ( C == 2 ) NN = ID_Hg_npa
               !------------------------------------------------

               ! Hg0
               Hg0aq(I,J,NN) = Hg0aq(I,J,NN)
     &                 + ( DELTAH * CDeepnpa(1) * MHg * A_M2 * FRAC_O )

               ! Hg2
               Hg2aq(I,J,NN) = Hg2aq(I,J,NN)
     &                 + ( DELTAH * CDeepnpa(2) * MHg * A_M2 * FRAC_O )

               ! HgP
!               IF ( C == 1 ) THEN
                  HgPaq(I,J,NN)   = HgPaq(I,J,NN)          
     &                 + ( DELTAH * CDeepnpa(3) * MHg * A_M2 * FRAC_O )
!               ENDIF

            ! North Pacific (east) 
            ELSE IF ( ( X >= 25.0 .and. X < 180.0 )  .and.
     &              ( Y >=  30.0 .and. Y <  70.0 ) ) THEN

               !- eds 8/31/10-----------------------------------
               IF ( C == 2 ) NN = ID_Hg_npa
               !------------------------------------------------
                        
               ! Hg0
               Hg0aq(I,J,NN) = Hg0aq(I,J,NN)
     &                 + ( DELTAH * CDeepnpa(1) * MHg * A_M2 * FRAC_O )

               ! Hg2
               Hg2aq(I,J,NN) = Hg2aq(I,J,NN)
     &                 + ( DELTAH * CDeepnpa(2) * MHg * A_M2 * FRAC_O )

               ! HgP
!               IF ( C == 1 ) THEN
                  HgPaq(I,J,NN)   = HgPaq(I,J,NN)          
     &                 + ( DELTAH * CDeepnpa(3) * MHg * A_M2 * FRAC_O )
!               ENDIF

            ! North Atlantic
            ELSE IF ( ( X >= -80.0 .and. X < 25.0 )  .and.
     &              ( Y >=  55.0 .and. Y <  70.0 ) ) THEN
                        
               !- eds 8/31/10-----------------------------------
               IF ( C == 2 ) NN = ID_Hg_nat
               !------------------------------------------------

               ! Hg0
               Hg0aq(I,J,NN) = Hg0aq(I,J,NN)
     &                 + ( DELTAH * CDeepnat(1) * MHg * A_M2 * FRAC_O )

               ! Hg2
               Hg2aq(I,J,NN) = Hg2aq(I,J,NN)
     &                 + ( DELTAH * CDeepnat(2) * MHg * A_M2 * FRAC_O )

               ! HgP
!               IF ( C == 1 ) THEN
                  HgPaq(I,J,NN)   = HgPaq(I,J,NN)          
     &                 + ( DELTAH * CDeepnat(3) * MHg * A_M2 * FRAC_O )

!               ENDIF

            ! South Atlantic
            ELSE IF ( ( X >= -80.0 .and. X < 25.0 )  .and.
     &              ( Y >=  -65.0 .and. Y <  -25.0 ) ) THEN    !(anls,100114)

               !- eds 8/31/10-----------------------------------
               IF ( C == 2 ) NN = ID_Hg_sat
               !------------------------------------------------

               ! Hg0
               Hg0aq(I,J,NN) = Hg0aq(I,J,NN)
     &                 + ( DELTAH * CDeepsat(1) * MHg * A_M2 * FRAC_O )

               ! Hg2
               Hg2aq(I,J,NN) = Hg2aq(I,J,NN)
     &                 + ( DELTAH * CDeepsat(2) * MHg * A_M2 * FRAC_O )

               ! HgP
!               IF ( C == 1 ) THEN
                  HgPaq(I,J,NN)   = HgPaq(I,J,NN)          
     &                 + ( DELTAH * CDeepsat(3) * MHg * A_M2 * FRAC_O )

!               ENDIF

            ! Antarctic
            ELSE IF ( Y >=  -90.0 .and. Y <  -65.0 ) THEN
                        
               !- eds 8/31/10-----------------------------------
               IF ( C == 2 ) NN = ID_Hg_ant
               !------------------------------------------------

               ! Hg0
               Hg0aq(I,J,NN) = Hg0aq(I,J,NN)
     &                 + ( DELTAH * CDeepant(1) * MHg * A_M2 * FRAC_O )

               ! Hg2
               Hg2aq(I,J,NN) = Hg2aq(I,J,NN)
     &                 + ( DELTAH * CDeepant(2) * MHg * A_M2 * FRAC_O )

               ! HgP
!               IF ( C == 1 ) THEN
                  HgPaq(I,J,NN)   = HgPaq(I,J,NN)          
     &                 + ( DELTAH * CDeepant(3) * MHg * A_M2 * FRAC_O )

!               ENDIF


            ! Arctic
            ELSE IF ( Y >=  70.0 .and. Y <  90.0 ) THEN
                        
               !- eds 8/31/10-----------------------------------
               IF ( C == 2 ) NN = ID_Hg_arc
               !------------------------------------------------

               ! Hg0
               Hg0aq(I,J,NN) = Hg0aq(I,J,NN)
     &                 + ( DELTAH * CDeeparc(1) * MHg * A_M2 * FRAC_O )

               ! Hg2
               Hg2aq(I,J,NN) = Hg2aq(I,J,NN)
     &                 + ( DELTAH * CDeeparc(2) * MHg * A_M2 * FRAC_O )

               ! HgP
!               IF ( C == 1 ) THEN
                  HgPaq(I,J,NN)   = HgPaq(I,J,NN)          
     &                 + ( DELTAH * CDeeparc(3) * MHg * A_M2 * FRAC_O )

!               ENDIF

            ELSE

               !- eds 8/31/10-----------------------------------
               IF ( C == 2 ) NN = ID_Hg_ocn
               !------------------------------------------------

               ! Hg0
               Hg0aq(I,J,NN) = Hg0aq(I,J,NN)
     &                    + ( DELTAH * CDeep(1) * MHg * A_M2 * FRAC_O )

               ! Hg2
               Hg2aq(I,J,NN) = Hg2aq(I,J,NN)
     &                    + ( DELTAH * CDeep(2) * MHg * A_M2 * FRAC_O )

               ! HgP
!               IF ( C == 1 ) THEN
                  HgPaq(I,J,NN)   = HgPaq(I,J,NN)          
     &                    + ( DELTAH * CDeep(3) * MHg * A_M2 * FRAC_O )
!               ENDIF

            ENDIF
         ENDDO
               
      ELSE 

         !==============================================================
         ! IF MIXED LAYER DEPTH HAS DECREASED:
         !
         ! Conserve concentration, but shed mass for ALL tracers.  
         ! Mass changes by same ratio as volume.
         !==============================================================

         ! Avoid dividing by zero
         IF ( MLDold > 0d0 ) THEN

            ! Update Hg0 and Hg2 categories
            DO NN = 1, N_Hg_CATS
               Hg0aq(I,J,NN) = Hg0aq(I,J,NN) * ( MLDnew / MLDold )
               Hg2aq(I,J,NN) = Hg2aq(I,J,NN) * ( MLDnew / MLDold )
               HgPaq(I,J,NN) = HgPaq(I,J,NN) * ( MLDnew / MLDold )
            ENDDO
            
         ENDIF

      ENDIF
      
      ! Return to calling program
      END SUBROUTINE MLD_ADJUSTMENT

!------------------------------------------------------------------------------

      SUBROUTINE READ_OCEAN_Hg_RESTART( YYYYMMDD, HHMMSS, Input_Opt ) 
!
!******************************************************************************
!  Subroutine READ_OCEAN_Hg_RESTART initializes GEOS-CHEM oceanic mercury 
!  tracer masses from a restart file. (sas, cdh, bmy, 3/28/06)
!
!  Arguments as input:
!  ============================================================================
!  (1 ) YYYYMMDD : Year-Month-Day 
!  (2 ) HHMMSS   :  and Hour-Min-Sec for which to read restart file
! 
!  NOTES:
!  01 Aug 2012 - R. Yantosca - Add reference to findFreeLUN from inqure_mod.F90
!  07 Aug 2012 - R. Yantosca - Now print LUN used to open file
!  25 Mar 2013 - R. Yantosca - Now use logical fields from Input_Opt
!******************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD,          ONLY : OPEN_BPCH2_FOR_READ
      USE CMN_SIZE_MOD
      USE DEPO_MERCURY_MOD,   ONLY : SNOW_HG, CHECK_DIMENSIONS
      USE DIRECTORY_MOD,      ONLY : RUN_DIR
      USE ERROR_MOD,          ONLY : DEBUG_MSG
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE FILE_MOD,           ONLY : IOERROR
      USE TIME_MOD,           ONLY : EXPAND_DATE
      USE TRACERID_MOD,       ONLY : GET_Hg0_CAT, GET_Hg2_CAT, N_Hg_CATS
      USE TRACERID_MOD,       ONLY : ID_Hg0,      ID_Hg2 
      !- eds 8/31/10 ---------------------------------------------------------
      USE TRACERID_MOD,       ONLY : ID_HgP
      USE TRACERID_MOD,       ONLY : GET_HgP_CAT !eds 10/18/10
      !-----------------------------------------------------------------------

      ! Arguments
      INTEGER,        INTENT(IN) :: YYYYMMDD, HHMMSS
      TYPE(OptInput), INTENT(IN) :: Input_Opt   ! Input Options object

      ! Local Variables
      INTEGER             :: IU_FILE
      INTEGER             :: I, IOS, J, L, NN, N_oc
      INTEGER             :: YEAR, MONTH, DAY
      INTEGER             :: NCOUNT(NNPAR) 
      REAL*4              :: Hg_OCEAN(IIPAR,JJPAR,1)
      CHARACTER(LEN=255)  :: FILENAME

      ! For binary punch file, version 2.0
      INTEGER             :: NI,        NJ,      NL
      INTEGER             :: IFIRST,    JFIRST,  LFIRST
      INTEGER             :: NTRACER,   NSKIP
      INTEGER             :: HALFPOLAR, CENTER180
      REAL*4              :: LONRES,    LATRES
      REAL*8              :: ZTAU0,     ZTAU1
      CHARACTER(LEN=20)   :: MODELNAME
      CHARACTER(LEN=40)   :: CATEGORY
      CHARACTER(LEN=40)   :: UNIT     
      CHARACTER(LEN=40)   :: RESERVED

      ! For fields from Input_Opt
      LOGICAL             :: LPRT
      LOGICAL             :: LSPLIT
      LOGICAL             :: LFUTURE
      INTEGER             :: N_TRACERS
      CHARACTER(LEN=14)   :: TRACER_NAME(Input_Opt%N_TRACERS)
      REAL*8              :: TRACER_MW_KG(Input_Opt%N_TRACERS)

      !=================================================================
      ! READ_OCEAN_Hg_RESTART begins here!
      !=================================================================

      ! Copy values from Input_Opt
      LPRT         = Input_Opt%LPRT
      LSPLIT       = Input_Opt%LSPLIT
      N_TRACERS    = Input_Opt%N_TRACERS
      TRACER_NAME  = Input_Opt%TRACER_NAME(1:N_TRACERS)
      TRACER_MW_KG = Input_Opt%TRACER_MW_KG(1:N_TRACERS)

      ! Initialize some variables
      NCOUNT(:)       = 0
      Hg_OCEAN(:,:,:) = 0e0

      ! Find a free file LUN
      IU_FILE         = findFreeLUN()

      ! Copy input file name to a local variable
      FILENAME        = TRIM( RUN_DIR ) // TRIM( Hg_RST_FILE )

      ! Replace YYYY, MM, DD, HH tokens in FILENAME w/ actual values
      CALL EXPAND_DATE( FILENAME, YYYYMMDD, HHMMSS )

      ! Echo some input to the screen
      WRITE( 6, '(a)' ) REPEAT( '=', 79 )
      WRITE( 6, 100   ) 
      WRITE( 6, 110   ) TRIM( FILENAME ), IU_FILE
 100  FORMAT( 'O C E A N   H g   R E S T A R T   F I L E   I N P U T' )
 110  FORMAT( /, 'READ_OCEAN_Hg_RESTART: Reading ', a, ' on unit ', i4 )

      ! Open the binary punch file for input
      CALL OPEN_BPCH2_FOR_READ( IU_FILE, FILENAME )
      
      ! Echo more output
      WRITE( 6, 120 )
 120  FORMAT( /, 'Min and Max of each tracer, as read from the file:',
     &        /, '(in volume mixing ratio units: v/v)' )
      
      !=================================================================
      ! Read concentrations -- store in the TRACER array
      !=================================================================
      DO 
         READ( IU_FILE, IOSTAT=IOS ) 
     &     MODELNAME, LONRES, LATRES, HALFPOLAR, CENTER180

         ! IOS < 0 is end-of-file, so exit
         IF ( IOS < 0 ) EXIT

         ! IOS > 0 is a real I/O error -- print error message
         IF ( IOS > 0 ) CALL IOERROR( IOS, IU_FILE, 'rd_oc_hg_rst:1' )

         READ( IU_FILE, IOSTAT=IOS ) 
     &        CATEGORY, NTRACER,  UNIT, ZTAU0,  ZTAU1,  RESERVED,
     &        NI,       NJ,       NL,   IFIRST, JFIRST, LFIRST,
     &        NSKIP

         IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_FILE, 'rd_oc_hg_rst:2' )

         READ( IU_FILE, IOSTAT=IOS ) 
     &        ( ( ( Hg_OCEAN(I,J,L), I=1,NI ), J=1,NJ ), L=1,NL )

         IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_FILE, 'rd_oc_hg_rst:3' )

         !==============================================================
         ! Assign data from the TRACER array to the STT array.
         !==============================================================
  
         ! Only process concentration data (i.e. mixing ratio)
         IF ( CATEGORY(1:8) == 'OCEAN-HG' ) THEN 

            ! Make sure array dimensions are of global size
            ! (NI=IIPAR; NJ=JJPAR, NL=LLPAR), or stop the run
            CALL CHECK_DIMENSIONS( NI, NJ, NL )

            ! Save into arrays
            IF ( ANY( ID_Hg0 == NTRACER ) ) THEN

               !----------
               ! Hg(0)
               !----------
               
               ! Get the Hg category #
               NN              = GET_Hg0_CAT( NTRACER )

               ! Store ocean Hg(0) in Hg0aq array
               Hg0aq(:,:,NN)   = Hg_OCEAN(:,:,1)
               
               ! Check for negative concentrations (jaf, 7/6/11)
               DO I = 1, IIPAR
               DO J = 1, JJPAR
                  IF ( Hg0aq(I,J,NN) .LT. 0 ) Hg0aq(I,J,NN) = 0d0
               ENDDO
               ENDDO

               ! Increment NCOUNT
               NCOUNT(NTRACER) = NCOUNT(NTRACER) + 1

            ELSE IF ( ANY( ID_Hg2 == NTRACER ) ) THEN
               
               !----------
               ! Hg(II)
               !----------

               ! Get the Hg category #
               NN              = GET_Hg2_CAT( NTRACER )

               ! Store ocean Hg(II) in Hg2_aq array
               Hg2aq(:,:,NN)   = Hg_OCEAN(:,:,1)

               ! Check for negative concentrations (jaf, 7/6/11)
               DO I = 1, IIPAR
               DO J = 1, JJPAR
                  IF ( Hg2aq(I,J,NN) .LT. 0 ) Hg2aq(I,J,NN) = 0d0
               ENDDO
               ENDDO

               ! Increment NCOUNT
               NCOUNT(NTRACER) = NCOUNT(NTRACER) + 1

            !- eds 8/31/10 ---------------------------------------
            !ELSE IF ( NTRACER == 3 ) THEN
            ELSE IF ( ANY( ID_HgP == NTRACER ) ) THEN
            !----------------------------------------------------

               !----------
               ! Hg(P)
               !----------
               ! Get the Hg category #   !eds 10/18/10
               NN              = GET_HgP_CAT( NTRACER )

               ! Particulate Hg
               HgPaq(:,:,NN)        = Hg_OCEAN(:,:,1)

               ! Increment NCOUNT
               NCOUNT(NTRACER) = NCOUNT(NTRACER) + 1

            ENDIF

            ! CDH snowpack (added following IF)
         ELSE IF ( CATEGORY(1:7) == 'SNOW-HG' ) THEN  
               !----------
               ! Hg in snow
               !----------

               ! Get the Hg category #
               NN              = GET_Hg0_CAT( NTRACER )

               ! Store ocean Hg(0) in Hg0aq array
               SNOW_HG(:,:,NN)   = Hg_OCEAN(:,:,1)
               
               ! Increment NCOUNT
!               NCOUNT(NTRACER) = NCOUNT(NTRACER) + 1
         ENDIF
      ENDDO

      ! Close file
      CLOSE( IU_FILE )      

      !=================================================================
      ! Examine data blocks, print totals, and return
      !=================================================================

      ! Tagged simulation has 17 ocean tracers; otherwise 3
      IF ( LSPLIT ) THEN
!- eds 8/31/2010 -------------------------------------------------------
!         N_oc = 17
         N_oc = 87  !Note: 87 = total Hg0, Hg2, HgP; should it be 1/3 of that? (eds)
!-----------------------------------------------------------------------
      ELSE
         N_oc = 3
      ENDIF

      ! eds temp debugging
      PRINT*, 'NNPAR', NNPAR
      PRINT*, 'NCOUNT', NCOUNT
      PRINT*, 'N_oc', N_oc

      ! Check for missing or duplicate data blocks
      CALL CHECK_DATA_BLOCKS( N_oc, NCOUNT )    

      !=================================================================
      ! Print totals
      !=================================================================

      ! Echo info
      WRITE( 6, 130 )

      ! Hg0
      DO NN = 1, N_Hg_CATS
         WRITE( 6, 140 ) ID_Hg0(NN), TRACER_NAME( Id_Hg0(NN) ), 
     &                   SUM( Hg0aq(:,:,NN) ),  'kg'
      ENDDO

      ! Hg2
      DO NN = 1, N_Hg_CATS
!     Bug fixed, Y. Zhang, 11/1/2011
!         WRITE( 6, 140 ) ID_Hg2(NN), TRACER_NAME( Id_Hg2(NN) ), 
!     &                   SUM( Hg0aq(:,:,NN) ), 'kg'
         WRITE( 6, 140 ) ID_Hg2(NN), TRACER_NAME( Id_Hg2(NN) ), 
     &                   SUM( Hg2aq(:,:,NN) ), 'kg'
      ENDDO

      ! HgP
      WRITE( 6, 140 ) 3, 'HgP       ', SUM( HgPaq ), 'kg'

      ! Format strings
 130  FORMAT( /, 'Total masses for each ocean tracer: ' ) 
 140  FORMAT( 'Tracer ', i3, ' (', a10, ') ', es12.5, 1x, a4)

      ! Fancy output
      WRITE( 6, '(a)' ) REPEAT( '=', 79 )

      ! Make sure tagged & total tracers sum up
      IF ( USE_CHECKS .and. LSPLIT ) THEN
         CALL CHECK_OCEAN_MERCURY( 'end of READ_OCEAN_Hg_RESTART' )
      ENDIF

      !### Debug
      IF ( LPRT ) CALL DEBUG_MSG( '### READ_OCEAN_Hg_RST: read file' )

      ! Return to calling program
      END SUBROUTINE READ_OCEAN_Hg_RESTART

!------------------------------------------------------------------------------

      SUBROUTINE CHECK_DATA_BLOCKS( N_TRACERS, NCOUNT )
!
!******************************************************************************
!  Subroutine CHECK_DATA_BLOCKS checks to see if we have multiple or 
!  missing data blocks for a given tracer. (sas, cdh, bmy, 3/28/06)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) N_TRACERS (INTEGER) : Number of tracers
!  (2 ) NCOUNT    (INTEGER) : Ctr array - # of data blocks found per tracer
!
!  NOTES:
!******************************************************************************
!      
      ! References to F90 modules
      USE ERROR_MOD, ONLY : GEOS_CHEM_STOP

      USE CMN_SIZE_MOD    ! Size parameters

      ! Arguments
      INTEGER, INTENT(IN) :: N_TRACERS, NCOUNT(NNPAR)
  
      ! Local variables
      INTEGER             :: N

      !=================================================================
      ! CHECK_DATA_BLOCKS begins here! 
      !=================================================================

      ! Loop over all tracers
      DO N = 1, N_TRACERS

         ! Stop if a tracer has more than one data block 
         IF ( NCOUNT(N) > 1 ) THEN 
            WRITE( 6, 100 ) N
            WRITE( 6, 120 ) 
            WRITE( 6, '(a)' ) REPEAT( '=', 79 )
            CALL GEOS_CHEM_STOP
         ENDIF
         
         ! Stop if a tracer has no data blocks 
         IF ( NCOUNT(N) == 0 ) THEN
            WRITE( 6, 110 ) N
            WRITE( 6, 120 ) 
            WRITE( 6, '(a)' ) REPEAT( '=', 79 )
            CALL GEOS_CHEM_STOP
         ENDIF
      ENDDO

      ! FORMAT statements
 100  FORMAT( 'More than one record found for tracer : ', i4 )
 110  FORMAT( 'No records found for tracer : ',           i4 ) 
 120  FORMAT( 'STOP in CHECK_DATA_BLOCKS (restart_mod.f)'    )

      ! Return to calling program
      END SUBROUTINE CHECK_DATA_BLOCKS

!------------------------------------------------------------------------------

      SUBROUTINE MAKE_OCEAN_Hg_RESTART( NYMD, NHMS, TAU, Input_Opt )
!
!******************************************************************************
!  Subroutine MAKE_OCEAN_Hg_RESTART writes an ocean mercury restart file.
!  (sas, cdh, bmy, 3/28/06)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) YYYYMMDD : Year-Month-Date 
!  (2 ) HHMMSS   :  and Hour-Min-Sec for which to create a restart file       
!  (3 ) TAU      : GEOS-CHEM TAU value corresponding to YYYYMMDD, HHMMSS
!  
!  NOTES:
!  01 Aug 2012 - R. Yantosca - Add reference to findFreeLUN from inqure_mod.F90
!  25 Mar 2013 - R. Yantosca - Now use logical fields from Input_Opt
!******************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD
      USE DEPO_MERCURY_MOD,   ONLY : SNOW_HG
      USE DIRECTORY_MOD,      ONLY : RUN_DIR
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE GRID_MOD,           ONLY : GET_XOFFSET, GET_YOFFSET
      USE TIME_MOD,           ONLY : EXPAND_DATE, GET_TAU
      USE TRACERID_MOD,       ONLY : ID_Hg_tot,   ID_Hg0
      USE TRACERID_MOD,       ONLY : ID_Hg2,      N_Hg_CATS
      USE TRACERID_MOD,       ONLY : ID_HgP !eds 9/9/10 

      USE CMN_SIZE_MOD           ! Size parameters

      ! Arguments           
      INTEGER,        INTENT(IN) :: NYMD, NHMS
      REAL*8,         INTENT(IN) :: TAU
      TYPE(OptInput), INTENT(IN) :: Input_Opt   ! Input Options object

      ! Local variables
      INTEGER             :: IU_FILE
      INTEGER             :: HALFPOLAR, CENTER180
      INTEGER             :: IFIRST,    JFIRST,   LFIRST
      INTEGER             :: N,         NN
      REAL*4              :: LONRES,    LATRES,   ARRAY(IIPAR,JJPAR,1)
      CHARACTER(LEN=20)   :: MODELNAME
      CHARACTER(LEN=40)   :: CATEGORY,  UNIT,     RESERVED
      CHARACTER(LEN=255)  :: FILENAME

      ! For values from Input_Opt
      LOGICAL             :: LSPLIT

      !=================================================================
      ! MAKE_OCEAN_Hg_RESTART begins here!
      !=================================================================

      ! Copy values from Input_Opt
      LSPLIT    = Input_Opt%LSPLIT

      ! Initialize values
      IFIRST    = GET_XOFFSET( GLOBAL=.TRUE. ) + 1
      JFIRST    = GET_YOFFSET( GLOBAL=.TRUE. ) + 1
      LFIRST    = 1
      HALFPOLAR = GET_HALFPOLAR()
      CENTER180 = 1
      LONRES    = DISIZE
      LATRES    = DJSIZE
      MODELNAME = GET_MODELNAME()
      CATEGORY  = 'OCEAN-HG'
      RESERVED  = ''
      UNIT      = 'kg'

      ! Find a free file LUN
      IU_FILE   = findFreeLUN()

      ! Expand date in filename
      FILENAME  = TRIM( RUN_DIR ) // Hg_RST_FILE
      CALL EXPAND_DATE( FILENAME, NYMD, NHMS )

      ! Echo info
      WRITE( 6, 100 ) TRIM( FILENAME ), IU_FILE
 100  FORMAT( '     - MAKE_RESTART_FILE: Writing ', a, ' on unit ', i6 )

      ! Open BPCH file for output
      CALL OPEN_BPCH2_FOR_WRITE( IU_FILE, FILENAME )

      !---------------------------
      ! Total Hg(0) in ocean
      !---------------------------
      N            = ID_Hg0(Id_Hg_tot)
      ARRAY(:,:,1) = Hg0aq(:,:,ID_Hg_tot)

      CALL BPCH2( IU_FILE,   MODELNAME, LONRES,   LATRES,
     &            HALFPOLAR, CENTER180, CATEGORY, N, 
     &            UNIT,      TAU,       TAU,      RESERVED,
     &            IIPAR,     JJPAR,     1,        IFIRST,
     &            JFIRST,    LFIRST,    ARRAY(:,:,1) )

      !---------------------------
      ! Total Hg(II) in ocean
      !---------------------------
      N            = ID_Hg2(ID_Hg_tot)
      ARRAY(:,:,1) = Hg2aq(:,:,ID_Hg_tot)

      CALL BPCH2( IU_FILE,   MODELNAME, LONRES,   LATRES,
     &            HALFPOLAR, CENTER180, CATEGORY, N, 
     &            UNIT,      TAU,       TAU,      RESERVED,
     &            IIPAR,     JJPAR,     1,        IFIRST,
     &            JFIRST,    LFIRST,    ARRAY(:,:,1) )

      !---------------------------
      ! Total HgP in ocean
      !---------------------------
      N            = ID_HgP(ID_Hg_tot)
      ARRAY(:,:,1) = HgPaq(:,:,ID_Hg_tot)   !eds 10/18/10

      CALL BPCH2( IU_FILE,   MODELNAME, LONRES,   LATRES,
     &            HALFPOLAR, CENTER180, CATEGORY, N, 
     &            UNIT,      TAU,       TAU,      RESERVED,
     &            IIPAR,     JJPAR,     1,        IFIRST,
     &            JFIRST,    LFIRST,    ARRAY(:,:,1) )

      ! Save tagged ocean tracers if present
      IF ( LSPLIT ) THEN

         !------------------------
         ! Tagged Hg(0) in ocean
         !------------------------
         DO NN = 2, N_Hg_CATS
            N            = ID_Hg0(NN)
            ARRAY(:,:,1) = Hg0aq(:,:,NN)

            CALL BPCH2( IU_FILE,   MODELNAME, LONRES,   LATRES,
     &                  HALFPOLAR, CENTER180, CATEGORY, N, 
     &                  UNIT,      TAU,       TAU,      RESERVED,
     &                  IIPAR,     JJPAR,     1,        IFIRST,
     &                  JFIRST,    LFIRST,    ARRAY(:,:,1) )
         ENDDO

         !------------------------
         ! Tagged Hg(II) in ocean
         !------------------------
         DO NN = 2, N_Hg_CATS
            N            = ID_Hg2(NN)
            ARRAY(:,:,1) = Hg2aq(:,:,NN)

            CALL BPCH2( IU_FILE,   MODELNAME, LONRES,   LATRES,
     &                  HALFPOLAR, CENTER180, CATEGORY, N, 
     &                  UNIT,      TAU,       TAU,      RESERVED,
     &                  IIPAR,     JJPAR,     1,        IFIRST,
     &                  JFIRST,    LFIRST,    ARRAY(:,:,1) )
         ENDDO

         !------------------------
         ! Tagged Hg(P) in ocean !eds 9/9/10 
         !------------------------
         DO NN = 2, N_Hg_CATS
            N            = ID_HgP(NN)  
            ARRAY(:,:,1) = HgPaq(:,:,NN)

            CALL BPCH2( IU_FILE,   MODELNAME, LONRES,   LATRES,
     &                  HALFPOLAR, CENTER180, CATEGORY, N,
     &                  UNIT,	   TAU,       TAU,	RESERVED,
     &                  IIPAR,     JJPAR,     1,        IFIRST,
     &                  JFIRST,    LFIRST,    ARRAY(:,:,1) )
         ENDDO
      ENDIF

      !---------------------------
      ! Total Hg in snowpack
      !---------------------------
         DO NN = 1, N_Hg_CATS
            CATEGORY     = 'SNOW-HG'
            N            = ID_Hg0(NN)
            ARRAY(:,:,1) = SNOW_HG(:,:,NN)

            CALL BPCH2( IU_FILE,   MODELNAME, LONRES,   LATRES,
     &                  HALFPOLAR, CENTER180, CATEGORY, N, 
     &                  UNIT,      TAU,       TAU,      RESERVED,
     &                  IIPAR,     JJPAR,     1,        IFIRST,
     &                  JFIRST,    LFIRST,    ARRAY(:,:,1) )
         ENDDO

      ! Close file
      CLOSE( IU_FILE )

      ! Make sure tagged & total tracers sum up
      IF ( USE_CHECKS .and. LSPLIT ) THEN
         CALL CHECK_OCEAN_MERCURY( 'end of MAKE_OCEAN_Hg_RESTART' )
      ENDIF

      ! Return to calling program
      END SUBROUTINE MAKE_OCEAN_Hg_RESTART

!------------------------------------------------------------------------------

      SUBROUTINE CHECK_ATMOS_MERCURY( State_Chm, LOC )
!
!******************************************************************************
!  Subroutine CHECK_ATMOS_MERCURY tests whether the total and tagged tracers 
!  the GEOS-CHEM tracer array STT sum properly within each grid box.
!  (cdh, bmy, 3/28/06)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) LOC (CHARACTER) : Name of routine where CHECK_ATMOS_MERCURY is called
!
!  NOTES:
!******************************************************************************
!
      ! References to F90 modules
      USE GIGC_State_Chm_Mod,  ONLY : ChmState
      USE ERROR_MOD,           ONLY : ERROR_STOP
      USE TRACERID_MOD,        ONLY : ID_Hg0,    ID_Hg2,   ID_HgP
      USE TRACERID_MOD,        ONLY : ID_Hg_tot, N_Hg_CATS

      USE CMN_SIZE_MOD            ! Size parameters

      ! Arguments as Input
      CHARACTER(LEN=*), INTENT(IN) :: LOC
      TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object

      ! Local variables
      LOGICAL                      :: FLAG
      INTEGER                      :: I,       J,       L
      INTEGER                      :: N,       NN
      REAL*8                       :: Hg0_tot, Hg0_tag, RELERR0, ABSERR0      
      REAL*8                       :: Hg2_tot, Hg2_tag, RELERR2, ABSERR2
      REAL*8                       :: HgP_tot, HgP_tag, RELERRP, ABSERRP

      ! Pointers
      ! We need to define local arrays to hold corresponding values 
      ! from the Chemistry State (State_Chm) object. (mpayer, 12/6/12)
      REAL*8, POINTER :: STT(:,:,:,:)

      !=================================================================
      ! CHECK_ATMOS_MERCURY begins here!
      !=================================================================

      ! Set error flags
      FLAG = .FALSE.

      ! Initialize GEOS-Chem tracer array [kg] from Chemistry State object
      ! (mpayer, 12/6/12)
      STT => State_Chm%Tracers

      ! Loop over grid boxes
! OMP PARALLEL DO
! OMP+DEFAULT( SHARED )
! OMP+PRIVATE( I,       J,       L,       N,      NN            )
! OMP+PRIVATE( Hg0_tot, RELERR0, ABSERR0                        )
! OMP+PRIVATE( Hg2_tot, RELERR2, ABSERR2                        )
! OMP+PRIVATE( HgP_tot, RELERRP, ABSERRP                        )
! OMP+REDUCTION( +:     Hg0_tag, Hg2_tag, HgP_tag               )
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Initialize
         Hg0_tot = 0d0
         Hg0_tag = 0d0
         RELERR0 = 0d0
         ABSERR0 = 0d0
         Hg2_tot = 0d0
         Hg2_tag = 0d0
         RELERR2 = 0d0
         ABSERR2 = 0d0
         HgP_tot = 0d0
         Hgp_tag = 0d0
         RELERRP = 0d0
         ABSERRP = 0d0

         !--------
         ! Hg(0)
         !--------

         ! Total Hg(0)
         N       = ID_Hg0(ID_Hg_tot)
         Hg0_tot = STT(I,J,L,N)

         ! Sum of tagged Hg(0)
         DO NN = 2, N_Hg_CATS
            N       = ID_Hg0(NN) 
            Hg0_tag = Hg0_tag + STT(I,J,L,N)
         ENDDO

         ! Absolute error for Hg0
         ABSERR0 = ABS( Hg0_tot - Hg0_tag )

         ! Relative error for Hg0 (avoid div by zero)
         IF ( Hg0_tot > 0d0 ) THEN
            RELERR0 = ABS( ( Hg0_tot - Hg0_tag ) / Hg0_tot )
         ELSE
            RELERR0 = -999d0
         ENDIF

         !--------
         ! Hg(II)
         !--------

         ! Total Hg(II)
         N       = ID_Hg2(ID_Hg_tot)
         Hg2_tot = STT(I,J,L,N)

         ! Sum of tagged Hg(II)
         DO NN = 2, N_Hg_CATS
            N       = ID_Hg2(NN) 
            Hg2_tag = Hg2_tag + STT(I,J,L,N)
         ENDDO

         ! Absolute error for Hg2
         ABSERR2 = ABS( Hg2_tot - Hg2_tag )

         ! Relative error for Hg2 (avoid div by zero)
         IF ( Hg2_tot > 0d0 ) THEN
            RELERR2 = ABS( ( Hg2_tot - Hg2_tag ) / Hg2_tot )
         ELSE
            RELERR2 = -999d0
         ENDIF

         !--------
         ! HgP
         !--------

         ! Total Hg(P)
         N       = ID_HgP(ID_Hg_tot)
         HgP_tot = STT(I,J,L,N)

         ! Sum of tagged Hg(P)
         DO NN = 2, N_Hg_CATS
            N = ID_HgP(NN)
            IF ( N > 0 ) HgP_tag = HgP_tag + STT(I,J,L,N)
         ENDDO

         ! Absolute error for HgP
         ABSERRP = ABS( HgP_tot - HgP_tag )

         ! Relative error for HgP (avoid div by zero)
         IF ( HgP_tot > 0d0 ) THEN
            RELERRP = ABS( ( HgP_tot - HgP_tag ) / HgP_tot )
         ELSE
            RELERRP = -999d0
         ENDIF

         !----------------------------
         ! Hg(0) error is too large
         !----------------------------
         IF ( RELERR0 > MAX_RELERR .and. ABSERR0 > MAX_ABSERR ) THEN
! OMP CRITICAL
            FLAG = .TRUE.
            WRITE( 6, 100 ) I, J, L, Hg0_tot, Hg0_tag, RELERR0, ABSERR0
! OMP END CRITICAL
         ENDIF

         !----------------------------
         ! Hg(0) error is too large
         !----------------------------
         IF ( RELERR2 > MAX_RELERR .and. ABSERR2 > MAX_ABSERR ) THEN
! OMP CRITICAL
            FLAG = .TRUE.
            WRITE( 6, 110 ) I, J, L, Hg2_tot, Hg2_tag, RELERR2, ABSERR2
! OMP END CRITICAL
         ENDIF

         !----------------------------
         ! HgP error is too large
         !----------------------------
         IF ( RELERRP > MAX_RELERR .and. ABSERRP > MAX_ABSERR ) THEN
! OMP CRITICAL
            FLAG = .TRUE.
            WRITE( 6, 120 ) I, J, L, HgP_tot, HgP_tag, RELERRP, ABSERRP
! OMP END CRITICAL
         ENDIF
      ENDDO
      ENDDO
      ENDDO
! OMP END PARALLEL DO

      ! Free pointer
      NULLIFY( STT )

      ! FORMAT strings
 100  FORMAT( 'Hg0 error: ', 3i5, 4es13.6 )
 110  FORMAT( 'Hg2 error: ', 3i5, 4es13.6 )
 120  FORMAT( 'HgP error: ', 3i5, 4es13.6 )
 
      ! Stop if Hg0 and Hg2 errors are too large
      IF ( FLAG ) THEN
         CALL ERROR_STOP( 'Tagged Hg0, Hg2, HgP do not add up!', LOC )
      ENDIF

      ! Return to calling program 
      END SUBROUTINE CHECK_ATMOS_MERCURY

!------------------------------------------------------------------------------

      SUBROUTINE CHECK_OCEAN_MERCURY( LOC )
!
!******************************************************************************
!  Subroutine CHECK_TAGGED_HG_OC tests whether tagged tracers in Hg0aq and
!  Hg2aq add properly within each grid box. (cdh, bmy, 3/28/06)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) LOC (CHARACTER) : Name of routine where CHECK_OCEAN_MERCURY is called
!
!  NOTES:
!******************************************************************************
!
      ! References to F90 modules
      USE ERROR_MOD,           ONLY : ERROR_STOP
      USE LOGICAL_MOD,         ONLY : LSPLIT
      USE TRACERID_MOD,        ONLY : ID_Hg_tot, N_Hg_CATS

      USE CMN_SIZE_MOD            ! Size parameters

      ! Arguments
      CHARACTER(LEN=*), INTENT(IN) :: LOC

      ! Local variables
      LOGICAL, SAVE                :: FIRST = .TRUE.
      LOGICAL                      :: FLAG
      INTEGER                      :: I,       J
      REAL*8                       :: Hg0_tot, Hg0_tag, RELERR0, ABSERR0      
      REAL*8                       :: Hg2_tot, Hg2_tag, RELERR2, ABSERR2

      !=================================================================
      ! CHECK_OCEAN_MERCURY begins here!
      !=================================================================

      ! Set error condition flag
      FLAG = .FALSE.

      ! Loop over ocean surface boxes
! OMP PARALLEL DO
! OMP+DEFAULT( SHARED )
! OMP+PRIVATE( I, J, Hg0_tot, Hg0_tag, RELERR0, ABSERR0 ) 
! OMP+PRIVATE        Hg2_tot, Hg2_tag, RELERR2, ABSERR2 )
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         !--------------------------------------
         ! Relative and absolute errors for Hg0 
         !--------------------------------------
         Hg0_tot = Hg0aq(I,J,ID_Hg_tot)
         Hg0_tag = SUM( Hg0aq(I,J,2:N_Hg_CATS) )
         ABSERR0 = ABS( Hg0_tot - Hg0_tag )

         ! Avoid div by zero
         IF ( Hg0_tot > 0d0 ) THEN
            RELERR0 = ABS( ( Hg0_tot - Hg0_tag ) / Hg0_tot )
         ELSE
            RELERR0 = -999d0
         ENDIF

         !--------------------------------------
         ! Relative and absolute errors for Hg2
         !--------------------------------------
         Hg2_tot = Hg2aq(I,J,ID_Hg_tot)
         Hg2_tag = SUM( Hg2aq(I,J,2:N_Hg_CATS) )
         ABSERR2 = ABS( Hg2_tot - Hg2_tag )

         ! Avoid div by zero
         IF ( Hg2_tot > 0d0 ) THEN
            RELERR2 = ABS( ( Hg2_tot - Hg2_tag ) / Hg2_tot )
         ELSE
            RELERR2 = -999d0
         ENDIF

         !--------------------------------------
         ! Hg(0) error is too large
         !--------------------------------------
         IF ( RELERR0 > MAX_RELERR .and. ABSERR0 > MAX_ABSERR ) THEN
! OMP CRITICAL
            FLAG = .TRUE.
            WRITE( 6, 100 ) I, J, Hg0_tot, Hg0_tag, RELERR0, ABSERR0
! OMP END CRITICAL
         ENDIF

         !--------------------------------------
         ! Hg(II) error is too large
         !--------------------------------------
         IF ( RELERR2 > MAX_RELERR .and. ABSERR2 > MAX_ABSERR ) THEN
! OMP CRITICAL
            FLAG = .TRUE.
            WRITE( 6, 110 ) I, J, Hg2_tot, Hg2_tag, RELERR2, ABSERR2
! OMP END CRITICAL
         ENDIF
      ENDDO
      ENDDO
! OMP END PARALLEL DO

      ! FORMAT strings
 100  FORMAT( 'Hg0aq error: ', 2i5, 4es13.6 )
 110  FORMAT( 'Hg2aq error: ', 2i5, 4es13.6 )

      ! Stop if Hg0 and Hg2 errors are too large
      IF ( FLAG ) THEN
         CALL ERROR_STOP( 'Tagged Hg0aq, Hg2aq do not add up!', LOC )
      ENDIF

      ! Return to calling program
      END SUBROUTINE CHECK_OCEAN_MERCURY

!------------------------------------------------------------------------------

      SUBROUTINE CHECK_OCEAN_FLUXES( LOC )
!
!******************************************************************************
!  Subroutine CHECK_OCEAN_FLUXES tests whether the drydep and wetdep fluxes in
!  DD_Hg2 and WD_Hg2 sum together in each grid box. (cdh, bmy, 3/28/06)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) LOC (CHARACTER) : Name of routine where CHECK_OCEAN_FLUXES is called
!
!  NOTES:
!******************************************************************************
!
      ! References to F90 modules
      USE ERROR_MOD,           ONLY : ERROR_STOP
      USE LOGICAL_MOD,         ONLY : LSPLIT
      USE TRACERID_MOD,        ONLY : ID_Hg_tot, N_Hg_CATS
      USE DEPO_MERCURY_MOD,    ONLY : DD_Hg2, WD_Hg2, DD_HgP, WD_HgP

      USE CMN_SIZE_MOD            ! Size parameters

      ! Arguments
      CHARACTER(LEN=*), INTENT(IN) :: LOC

      ! Local variables
      LOGICAL                      :: FLAG
      INTEGER                      :: I,         J
      REAL*8                       :: DD_tot,    DD_tag 
      REAL*8                       :: DD_RELERR, DD_ABSERR      
      REAL*8                       :: WD_tot,    WD_tag
      REAL*8                       :: WD_RELERR, WD_ABSERR

      !=================================================================
      ! CHECK_OCEAN_MERCURY begins here!
      !=================================================================

      ! Echo
      WRITE( 6, 100 )
 100  FORMAT( '     - In CHECK_OCEAN_FLUXES' )

      ! Set error condition flag
      FLAG = .FALSE.

      ! Loop over ocean surface boxes
! OMP PARALLEL DO
! OMP+DEFAULT( SHARED )
! OMP+PRIVATE( I, J, DD_tot, DD_tag, DD_RELERR, DD_ABSERR )
! OMP+PRIVATE(       WD_tot, WD_tag, WD_RELERR, WD_ABSERR )
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         !---------------------------------------
         ! Absolute & relative errors for DD_Hg2
         !---------------------------------------
         DD_tot    = DD_Hg2(I,J,1)
         DD_tag    = SUM( DD_Hg2(I,J,2:N_Hg_CATS) )
         DD_ABSERR = ABS( DD_tot - DD_tag ) 

         ! Avoid div by zero
         IF ( DD_tot > 0d0 ) THEN
            DD_RELERR = ABS( ( DD_tot - DD_tag ) / DD_tot )
         ELSE
            DD_RELERR = -999d0
         ENDIF

         !---------------------------------------
         ! Absolute & relative errors for WD_Hg2
         !---------------------------------------
         WD_tot    = WD_Hg2(I,J,1)
         WD_tag    = SUM( WD_Hg2(I,J,2:N_Hg_CATS) )
         WD_ABSERR = ABS( WD_tot - WD_tag )

         ! Avoid div by zero
         IF ( WD_tot > 0d0 ) THEN
            WD_RELERR = ABS( ( WD_tot - WD_tag ) / WD_tot )
         ELSE
            WD_RELERR = -999d0
         ENDIF

         !---------------------------------------
         ! DD flux error is too large
         !---------------------------------------
         IF ( DD_RELERR > MAX_RELERR .and. DD_ABSERR > MAX_FLXERR ) THEN
! OMP CRITICAL
            FLAG = .TRUE.
            WRITE( 6, 110 ) I, J, DD_tot, DD_tag, DD_RELERR, DD_ABSERR
! OMP END CRITICAL
         ENDIF

         !---------------------------------------
         ! WD flux error is too large
         !---------------------------------------
         IF ( WD_RELERR > MAX_RELERR .and. WD_ABSERR > MAX_FLXERR ) THEN
! OMP CRITICAL
            FLAG = .TRUE.
            WRITE( 6, 120 ) I, J, WD_tot, WD_tag, WD_RELERR, WD_ABSERR
! OMP END CRITICAL
         ENDIF
      ENDDO
      ENDDO
! OMP END PARALLEL DO

      ! FORMAT strings
 110  FORMAT( 'DD_Hg2 flux error: ', 2i5, 4es13.6 )
 120  FORMAT( 'WD_Hg2 flux error: ', 2i5, 4es13.6 )

      ! Stop if Hg0 and Hg2 errors are too large
      IF ( FLAG ) THEN
         CALL ERROR_STOP( 'Tagged DD, WD fluxes do not add up!', LOC )
      ENDIf

      ! Return to calling program
      END SUBROUTINE CHECK_OCEAN_FLUXES

!------------------------------------------------------------------------------

      SUBROUTINE CHECK_FLUX_OUT( FLUX, LOC )
!
!******************************************************************************
!  Subroutine CHECK_FLUX_OUT tests whether tagged quantities in FLUX sum 
!  together in each grid box. (cdh, bmy, 3/20/06)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) FLUX (REAL*8)   : Flux array (output of OCEAN_MERCURY_FLUX)
!  (2 ) LOC (CHARACTER) : Name of routine where CHECK_FLUX_OUT is called
!
!  NOTES:
!******************************************************************************
!
      ! References to F90 modules
      USE ERROR_MOD,           ONLY : ERROR_STOP
      USE LOGICAL_MOD,         ONLY : LSPLIT
      USE TRACERID_MOD,        ONLY : ID_Hg_tot, N_Hg_CATS

      USE CMN_SIZE_MOD            ! Size parameters

      ! Arguments
      REAL*8,           INTENT(IN) :: FLUX(IIPAR,JJPAR,N_Hg_CATS)            
      CHARACTER(LEN=*), INTENT(IN) :: LOC

      ! Local variables
      LOGICAL                      :: FLAG
      INTEGER                      :: I,          J
      REAL*8                       :: FLX_tot,    FLX_tag
      REAL*8                       :: FLX_RELERR, FLX_ABSERR

      !=================================================================
      ! CHECK_FLUX_OUT begins here!
      !=================================================================

      ! Echo
      WRITE( 6, 100 )
 100  FORMAT( '     - In CHECK_FLUX_OUT' )

      ! Set error condition flag
      FLAG = .FALSE.

      ! Loop over ocean surface boxes
! OMP PARALLEL DO
! OMP+DEFAULT( SHARED )
! OMP+PRIVATE( I, J, FLX_tot, FLX_tag, FLX_err )
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         !----------------------------------------
         ! Absolute & relative errors for FLX_Hg2
         !----------------------------------------
         FLX_tot    = FLUX(I,J,1)
         FLX_tag    = SUM( FLUX(I,J,2:N_Hg_CATS) )
         FLX_ABSERR = ABS( FLX_tot - FLX_tag )
         
         ! Avoid div by zero
         IF ( FLX_tot > 0d0 ) THEN
            FLX_RELERR = ABS( ( FLX_tot - FLX_tag ) / FLX_tot )
         ELSE
            FLX_RELERR = -999d0
         ENDIF

         !----------------------------
         ! Flux error is too large
         !----------------------------
         IF ( FLX_RELERR > MAX_RELERR  .and. 
     &        FLX_ABSERR > MAX_ABSERR ) THEN
! OMP CRITICAL
            FLAG = .TRUE.
            WRITE( 6, 110 ) I, J, FLX_tot,    FLX_tag, 
     &                            FLX_RELERR, FLX_ABSERR
! OMP END CRITICAL
         ENDIF

      ENDDO
      ENDDO
! OMP END PARALLEL DO

      ! FORMAT strings
 110  FORMAT( 'FLX_Hg2 flux error: ', 2i5, 4es13.6 )
 
      ! Stop if Hg0 and Hg2 errors are too large
      IF ( FLAG ) THEN
         CALL ERROR_STOP( 'Tagged emission fluxes do not add up!', LOC )
      ENDIf

      ! Return to calling program
      END SUBROUTINE CHECK_FLUX_OUT

!------------------------------------------------------------------------------

      SUBROUTINE INIT_OCEAN_MERCURY( THIS_Hg_RST_FILE, THIS_USE_CHECKS )
!
!******************************************************************************
!  Subroutine INIT_OCEAN_MERCURY allocates and zeroes module arrays.  
!  (sas, cdh, bmy, 1/19/05, 3/28/06)
!
!  NOTES:
!  (1 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (2 ) Now just allocates arrays.  We have moved the reading of the ocean
!        Hg restart file to READ_OCEAN_Hg_RESTART.  Now make Hg0aq and Hg2aq
!        3-D arrays. Now pass Hg_RST_FILE and USE_CHECKS from "input_mod.f"
!        via the argument list. (cdh, sas, bmy, 2/27/06)
!  29 Nov 2011 - R. Yantosca - Updated erroneous comment.  "Southern Ocean"
!                               should be "South Atlantic Ocean".
!******************************************************************************
!
      ! References to F90 modules
      USE ERROR_MOD,    ONLY : ALLOC_ERR
      USE TRACERID_MOD, ONLY : N_Hg_CATS
      USE LOGICAL_MOD,  ONLY : LPREINDHG

      USE CMN_SIZE_MOD     ! Size parameters

      ! Arguments
      CHARACTER(LEN=*), INTENT(IN) :: THIS_Hg_RST_FILE
      LOGICAL,          INTENT(IN) :: THIS_USE_CHECKS

      ! Local variables
      INTEGER                      :: AS

      !=================================================================
      ! INIT_OCEAN_MERCURY begins here!
      !=================================================================

      ! Ocean Hg restart file name
      Hg_RST_FILE = THIS_Hg_RST_FILE
      
      ! Turn on error checks for tagged & total sums?
      USE_CHECKS  = THIS_USE_CHECKS
 

      ! Set up concentrations of Hg(0), Hg(II), Hg(C) in deep ocean REDALERT
 
      IF (LPREINDHG) THEN
         CDEEP    = (/ 2d-11, 1.67d-10, 1.67d-10 /)
         CDEEPATL = (/ 2d-11, 1.67d-10, 1.67d-10 /)
         CDEEPNAT = (/ 2d-11, 1.67d-10, 1.67d-10 /)
         CDEEPSAT = (/ 2d-11, 1.67d-10, 1.67d-10 /)         
         CDEEPANT = (/ 2d-11, 1.67d-10, 1.67d-10 /)
         CDEEPARC = (/ 2d-11, 1.67d-10, 1.67d-10 /)
         CDEEPNPA = (/ 2d-11, 1.67d-10, 1.67d-10 /)
      ELSE
         CDEEP    = (/ 1.0d-10, 4.0d-10, 4.0d-10 /)
         CDEEPATL = (/ 1.4d-10, 9.3d-10, 9.3d-10 /)  
         CDEEPNAT = (/ 1.5d-10, 8.2d-10, 8.2d-10 /)
!         CDEEPSAT = (/ 1.0d-10, 5.0d-10, 5.0d-10 /)   !(anls,100114)
         CDEEPSAT = (/0.8d-10,4.1d-10,4.1d-10/)  !eck, 10/19/11
         !reduce intermediate water mercury concentration in
         !South Atlantic Ocean to 0.9pM total
         !ref: low end of uncertainty range Sunderland and Mason '07
         CDEEPANT = (/ 1.0d-10, 3.2d-10, 3.2d-10 /)
         CDEEPARC = (/ 1.2d-10, 7.5d-10, 7.5d-10 /)
         CDEEPNPA = (/ 1.0d-10, 6.0d-10, 6.0d-10 /)        
      ENDIF

      
      ! Allocate arrays     
      ALLOCATE( dMLD( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'dMLD' )
      dMLD = 0d0

      ALLOCATE( Hg0aq( IIPAR, JJPAR, N_Hg_CATS ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'Hg0aq' )
      Hg0aq = 0d0

      ALLOCATE( Hg2aq( IIPAR, JJPAR, N_Hg_CATS ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'Hg2aq' )
      Hg2aq = 0d0

!- eds 8/31/10 ---------------------------------------------
! Added NN to HgPaq throughout module
!      ALLOCATE( HgPaq( IIPAR, JJPAR ), STAT=AS )
      ALLOCATE( HgPaq (IIPAR, JJPAR, N_Hg_CATS ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'HgPaq' )
      HgPaq = 0d0

      !eds 5/15/12 fix
      ALLOCATE( Hgaq_tot (IIPAR, JJPAR, N_Hg_CATS ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'Hgaq_tot' )
      HgPaq = 0d0

      ALLOCATE( HgPaq_SUNK (IIPAR, JJPAR, N_Hg_CATS ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'HgPaq_SUNK' )
      HgPaq = 0d0

      ALLOCATE( MLD( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'MLD' )
      MLD = 0d0

      ALLOCATE( MLDav( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'MLDav' )
      MLDav = 0d0

      ALLOCATE( newMLD( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'newMLD' )
      newMLD = 0d0

      !eds 10/19/10 fixing restart bug
      ALLOCATE( prevMLD( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /=0 ) CALL ALLOC_ERR( 'prevMLD' )
      prevMLD = 0d0

      ALLOCATE( NPP( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'NPP' )
      NPP = 0d0

      ALLOCATE( CHL( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'CHL' )
      CHL = 0d0

      ALLOCATE( UPVEL( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'UPVEL' )
      UPVEL = 0d0

      ALLOCATE( Fg( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'Fg' )
      Fg = 0d0    

      ALLOCATE( Fp( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'Fp' )
      Fp = 0d0  
     
      ALLOCATE( BULK_CONC( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'BULK_CONC' )
      BULK_CONC = 0d0  

      ALLOCATE( SO4_GC( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'SO4_GC' )
      SO4_GC = 0d0 

      ALLOCATE( NIT_CONC( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'NIT_CONC' )
      NIT_CONC = 0d0 

      ALLOCATE( NH4_CONC( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'NH4_CONC' )
      NH4_CONC = 0d0 

      ALLOCATE( OC_CONC( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'OC_CONC' )
      OC_CONC = 0d0 

      ALLOCATE( BC_CONC( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'BC_CONC' )
      BC_CONC = 0d0 

      ALLOCATE( DST_CONC( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'DST_CONC' )
      DST_CONC = 0d0 

      ALLOCATE( R( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'R' )
      R = 0d0 

      ALLOCATE( SO4_CONC( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'SO4_CONC' )
      SO4_CONC = 0d0 


      ! Return to calling program
      END SUBROUTINE INIT_OCEAN_MERCURY

!------------------------------------------------------------------------------

      SUBROUTINE CLEANUP_OCEAN_MERCURY
!
!******************************************************************************
!  Subroutine CLEANUP_OCEAN_MERCURY deallocates all arrays.  
!  (sas, cdh, bmy, 1/20/05, 3/28/06)
!  
!  NOTES:
!  (1 ) Now call GET_HALFPOLAR from "bpch2_mod.f" to get the HALFPOLAR flag 
!        value for GEOS or GCAP grids. (bmy, 6/28/05)
!  (2 ) Now just deallocate arrays.  We have moved the writing of the Hg
!        restart file to MAKE_OCEAN_Hg_RESTART.  Now also deallocate HgP, dMLD
!        and MLDav arrays. (sas, cdh, bmy, 3/28/06)
!  25 Oct 2011 - H. Amos - bring Hg2 gas-particle partitioning code into
!                          v9-01-02
!******************************************************************************
!     
      !=================================================================
      ! CLEANUP_OCEAN_MERCURY begins here!
      !=================================================================
      IF ( ALLOCATED( dMLD      ) ) DEALLOCATE( dMLD      )
      IF ( ALLOCATED( Hg0aq     ) ) DEALLOCATE( Hg0aq     )  
      IF ( ALLOCATED( Hg2aq     ) ) DEALLOCATE( Hg2aq     )
      IF ( ALLOCATED( HgPaq     ) ) DEALLOCATE( HgPaq     )  
      IF ( ALLOCATED( MLD       ) ) DEALLOCATE( MLD       )
      IF ( ALLOCATED( MLDav     ) ) DEALLOCATE( MLDav     )
      IF ( ALLOCATED( newMLD    ) ) DEALLOCATE( newMLD    )
      IF ( ALLOCATED( prevMLD   ) ) DEALLOCATE( prevMLD   ) !eds 10/19/10
      IF ( ALLOCATED( NPP       ) ) DEALLOCATE( NPP       )
      IF ( ALLOCATED( CHL       ) ) DEALLOCATE( CHL       )
      IF ( ALLOCATED( UPVEL     ) ) DEALLOCATE( UPVEL     )
      IF ( ALLOCATED( BULK_CONC ) ) DEALLOCATE( BULK_CONC )
      IF ( ALLOCATED( Fg        ) ) DEALLOCATE( Fg        )
      IF ( ALLOCATED( Fp        ) ) DEALLOCATE( Fp        )
      IF ( ALLOCATED( SO4_GC    ) ) DEALLOCATE( SO4_GC    )
      IF ( ALLOCATED( NIT_CONC  ) ) DEALLOCATE( NIT_CONC  )
      IF ( ALLOCATED( NH4_CONC  ) ) DEALLOCATE( NH4_CONC  )
      IF ( ALLOCATED( OC_CONC   ) ) DEALLOCATE( OC_CONC   )
      IF ( ALLOCATED( BC_CONC   ) ) DEALLOCATE( BC_CONC   )
      IF ( ALLOCATED( DST_CONC  ) ) DEALLOCATE( DST_CONC  )
      IF ( ALLOCATED( R         ) ) DEALLOCATE( R         )
      IF ( ALLOCATED( SO4_CONC  ) ) DEALLOCATE( SO4_CONC  )

      ! Return to calling program
      END SUBROUTINE CLEANUP_OCEAN_MERCURY

!------------------------------------------------------------------------------
     
      ! End of module
      END MODULE OCEAN_MERCURY_MOD
