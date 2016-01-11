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
!  (3 ) MAX_RELERR  (REAL(fpp)): Max error for total-tag error check [unitless]
!  (4 ) MAX_ABSERR  (REAL(fpp)): Max abs error for total-tag err chk [unitless]
!  (5 ) MAX_FLXERR  (REAL(fpp)): Max error tol for flux error check  [unitless]
!  (6 ) Hg2aq_tot   (REAL(fpp)): Total Hg2 conc. in the mixed layer  [kg      ]
!  (7 ) DD_Hg2      (REAL(fpp)): Array for Hg(II) dry dep'd to ocean [kg      ]
!  (8 ) Hgaq_tot    (REAL(fpp)): Total Hg conc. in the mixed layer   [kg      ]
!  (9 ) Hg0aq       (REAL(fpp)): Array for ocean mass of Hg(0)       [kg      ]
!  (10) Hg2aq       (REAL(fpp)): Array for ocean mass of Hg(II)      [kg      ]
!  (11) HgPaq       (REAL(fpp)): Array for ocean mass of HgP         [kg      ]
!  (12) dMLD        (REAL(fpp)): Array for Change in ocean MLD       [cm      ]
!  (13) MLD         (REAL(fpp)): Array for instantaneous ocean MLD   [cm      ]
!  (14) MLDav       (REAL(fpp)): Array for monthly mean ocean MLD    [cm      ]
!  (15) newMLD      (REAL(fpp)): Array for next month's ocean MLD    [cm      ]
!  (16) NPP         (REAL(fpp)): Array for mean net primary prod.    [unitless]
!  (17) RAD         (REAL(fpp)): Array for mean solar radiation      [W/m2    ]
!  (18) UPVEL       (REAL(fpp)): Array for ocean upwelling velocity  [m/s     ]
!  (19) WD_Hg2      (REAL(fpp)): Array for Hg(II) wet dep'd to ocean [kg      ]
!  (20) CHL         (REAL(fpp)): Chl surface concentration           [mg(m3   ]
!  (21) CDEEPATL    (REAL(fpp)): Conc. Hg0, Hg2, HgP below MLD-Atl   [pM      ]
!  (22) CDEEP       (REAL(fpp)): Conc. of Hg0, Hg2, HgP below MLD    [pM      ]
!  (23) CDEEPNAT    (REAL(fpp)): Conc. Hg0, Hg2, HgP below MLD-NAtl  [pM      ]
!  (24) CDEEPSAT    (REAL(fpp)): Conc. Hg0, Hg2, HgP below MLD-SAtl  [pM      ]
!  (25) CDEEPANT    (REAL(fpp)): Conc. Hg0, Hg2, HgP below MLD-Ant   [pM      ]
!  (26) CDEEPARC    (REAL(fpp)): Conc. Hg0, Hg2, HgP below MLD-Arc   [pM      ]
!
!
!  Module Routines:
!  ============================================================================
!  (1 ) READ_HG2_PARTITIONING  : Calculates fraction of Hg2 in gas phase
!  (2 ) DELIVER_SNOW_HG        : Delivers Hg stored in snowpack to ocean
!  (3 ) OCEAN_MERCURY_FLUX     : Routine to compute flux of oceanic mercury
!  (4 ) OCEAN_MERCURY_READ     : Routine to read MLD, NPP, RADSWG data fields
!  (5 ) GET_MLD_FOR_NEXT_MONTH : Routine to read MLD for the next month
!  (6 ) GET_MLD_FOR_LAST_MONTH : Routine to read MLD for the last month
!  (7 ) MLD_ADJUSTMENT         : Adjusts MLD 
!  (8 ) READ_OCEAN_Hg_RESTART  : Reads restart file with ocean Hg tracers
!  (9 ) CHECK_DATA_BLOCKS      : Checks for missing/multiple data blocks
!  (10) MAKE_OCEAN_Hg_RESTART  : Writes new restart file with ocean Hg tracers
!  (11) CHECK_ATMOS_MERCURY    : Checks mass of total & tagged atm Hg0 & Hg2 
!  (12) CHECK_OCEAN_MERCURY    : Checks mass of total & tagged oc Hg0 & Hg2
!  (13) CHECK_OCEAN_FLUXES     : Checks mass of total & tagged DD & WD fluxes
!  (14) CHECK_FLUX_OUT         : Checks mass of total & tagged sums
!  (15) INIT_OCEAN_MERCURY     : Allocates and zeroes all module variables
!  (16) CLEANUP_OCEAN_MERCURY  : Deallocates all module variables
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
!  (9 ) Toole, J. M. et al. (2010), Influences of the ocean surface mixed layer
!        and thermohaline stratification on Arctic Sea ice in the central
!        Canada Basin, J. Geophys. Res., 115, C10018, doi:10.1029/2009JC005660.
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
!  17 Nov 2014 - M. Yannetti - Added PRECISION_MOD
!  23 Feb 2015 - J. Fisher   - Add Arctic MLD (orig. 28/11/11)
!  05 Aug 2015 - J. Fisher   - Compute ocean chemistry in all ocean boxes
!			       irrespective of ice cover and reduce radiation,
!			       Hg deposition accordingly in partial-ice boxes
!                              (orig. 28/11/11)
!  06 Aug 2015 - J. Fisher   - Added delivery of snow Hg to ocean (orig. 28/11/11)
!  11 Aug 2015 - J. Fisher   - Remove obsolete parameters LOCEANFIX, FUDGE
!  12 Aug 2015 - J. Fisher   - Added input of Hg from Arctic rivers to ocean
!                              (orig. 8/12/11)
!  17 Aug 2015 - J. Fisher   - Add Arctic NPP, Chl from Anne Soerensen
!                              (orig. 3/19/13)
!******************************************************************************
!
      USE inquireMod, ONLY : findFreeLUN

      USE PRECISION_MOD, ONLY : fpp => fp   ! For GEOS-Chem Precision (fpp)
      USE PRECISION_MOD, ONLY : f4          !  Rename to avoid conflicts
      USE PRECISION_MOD, ONLY : f8

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
      PUBLIC :: LOCEANCOEF        ! sjs 
         
      !=================================================================
      ! MODULE VARIABLES
      !=================================================================

      ! Scalars
      LOGICAL              :: USE_CHECKS

      ! Parameters
      REAL*4,  PARAMETER   :: MAX_RELERR = 5.0d-2
      REAL*4,  PARAMETER   :: MAX_ABSERR = 5.0d-3
      REAL*4,  PARAMETER   :: MAX_FLXERR = 5.0d-1

      REAL(fpp)   :: CDEEP(3)  
      REAL(fpp)   :: CDEEPATL(3)
      REAL(fpp)   :: CDEEPNAT(3)
      REAL(fpp)   :: CDEEPSAT(3)
      REAL(fpp)   :: CDEEPANT(3)
      REAL(fpp)   :: CDEEPARC(3)
      REAL(fpp)   :: CDEEPNPA(3)

      ! For Arctic rivers (jaf, 12/8/11)
      REAL(fpp)   :: RIVERFLOW(12)
      REAL(fpp)   :: dFLOW
      REAL(fpp)   :: dFLOW1
      REAL(fpp)   :: dFLOW2
      REAL(fpp)   :: FLOWNOW 

      ! Arrays
      REAL(fpp),  ALLOCATABLE :: dMLD(:,:)
      REAL(fpp),  ALLOCATABLE :: Hg0aq(:,:,:)
      REAL(fpp),  ALLOCATABLE :: Hg2aq(:,:,:)
!- eds 8/31/10 ------------------------------------------------------
!      REAL(fpp),  ALLOCATABLE :: HgPaq(:,:)
      REAL(fpp),  ALLOCATABLE :: HgPaq(:,:,:)
      REAL(fpp),  ALLOCATABLE :: HgPaq_SUNK(:,:,:) !eds 5/15/12 fix
      REAL(fpp),  ALLOCATABLE :: Hgaq_tot(:,:,:)
!--------------------------------------------------------------------
!      REAL(fpp),  ALLOCATABLE :: MLD(:,:)
      REAL(fpp),  ALLOCATABLE :: MLDav(:,:)
      REAL(fpp),  ALLOCATABLE :: newMLD(:,:)
      REAL(fpp),  ALLOCATABLE :: prevMLD(:,:) !eds 10/19/10 debugging restart
!      REAL(fpp),  ALLOCATABLE :: NPP(:,:)
      REAL(fpp),  ALLOCATABLE :: RAD(:,:)
!      REAL(fpp),  ALLOCATABLE :: UPVEL(:,:)
! added by hma for Hg2 partitioning --------------------------------------
      REAL(fpp),  ALLOCATABLE :: BULK_CONC(:,:,:)  
      REAL(fpp),  ALLOCATABLE :: Fp(:,:,:)     
      REAL(fpp),  ALLOCATABLE :: Fg(:,:,:)   
      REAL(fpp),  ALLOCATABLE :: SO4_GC(:,:,:)  
      REAL(fpp),  ALLOCATABLE :: NH4_CONC(:,:,:)  
      REAL(fpp),  ALLOCATABLE :: NIT_CONC(:,:,:)
      REAL(fpp),  ALLOCATABLE :: OC_CONC(:,:,:) 
      REAL(fpp),  ALLOCATABLE :: BC_CONC(:,:,:)
      REAL(fpp),  ALLOCATABLE :: DST_CONC(:,:,:) 
      REAL(fpp),  ALLOCATABLE :: R(:,:,:) 
      REAL(fpp),  ALLOCATABLE :: SO4_WAC(:,:,:) 
      REAL(fpp),  ALLOCATABLE :: SO4_CONC(:,:,:) 
!-------------------------------------------------------------------------

      ! For IAV in NPP (jaf, 3/19/13)
      REAL(fpp)   :: NPP_SCF

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
      LOGICAL   :: LOCEANCOEF    ! sjs
      REAL(fpp) :: STRAT_BR_FACTOR

      ! CDH Set this TRUE to use corrected area-flux relationship
      ! Set this to FALSE to use original Strode et al. (2007) model
      LOGICAL, PARAMETER :: LOCEANFIX=.TRUE. 
      ! CDH average ocean area per grid box: 1.67d11 m2/box
      ! used when eliminating AREA * FRAC_O
      REAL(fpp), PARAMETER :: FUDGE=1.67e+11_fpp

      ! Pointers to fields in the HEMCO data structure.
      ! These need to be declared REAL(f4), aka REAL*4.
      REAL(f4), POINTER :: CHL   (:,:) => NULL()
      REAL(f4), POINTER :: CHL_A (:,:) => NULL()
      REAL(f4), POINTER :: MLD   (:,:) => NULL()
      REAL(f4), POINTER :: NPP   (:,:) => NULL()
      REAL(f4), POINTER :: NPP_A (:,:) => NULL()
      REAL(f4), POINTER :: UPVEL (:,:) => NULL()
      REAL(f4), POINTER :: dMLD1 (:,:) => NULL()
      REAL(f4), POINTER :: dMLD2 (:,:) => NULL()

      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement
      !=================================================================
      CONTAINS

!-------------------------------------------------------------------------

      SUBROUTINE READ_HG2_PARTITIONING( am_I_Root, Input_Opt,
     &                                  State_Met, THISMONTH, RC )

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
!  15 Apr 2014 - R. Yantosca - Eliminate array temporaries 
!  23 Jun 2014 - R. Yantosca - Now accept am_I_Root, Input_Opt, RC
!  06 Nov 2014 - R. Yantosca - Now use State_Met%AIRDEN(I,J,L)
!  17 Dec 2014 - R. Yantosca - Leave time/date variables as 8-byte
!  05 Mar 2015 - R. Yantosca - Add Input_Opt%RES_DIR to data path
!  13 Mar 2015 - R. Yantosca - Now remove bpch input
!*****************************************************************************
!
! !USES:
!     
      USE DIAG03_MOD,         ONLY : ND03                
      USE ERROR_MOD,          ONLY : GEOS_CHEM_STOP, ERROR_STOP
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE HCO_EmisList_Mod,   ONLY : HCO_GetPtr
      USE TIME_MOD,           ONLY : SET_Hg2_DIAG        
      USE TRANSFER_MOD,       ONLY : TRANSFER_3D
      USE CMN_SIZE_MOD
      USE CMN_DIAG_MOD
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)  :: am_I_Root   ! Are we on the root CPU?
      TYPE(OptInput), INTENT(IN)  :: Input_Opt   ! Input Options object
      TYPE(MetState), INTENT(IN)  :: State_Met   ! Meteorology State object
      INTEGER,        INTENT(IN)  :: THISMONTH   ! Current month
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT) :: RC          ! Success or failure?
!
! !LOCAL VARIABLES:
!
      ! Pointers to fields in the HEMCO data structure
      ! These have to be declared REAL(f4), aka REAL*4.
      REAL(f4), POINTER       :: ARRAYso4 (:,:,:) => NULL() ! so4
      REAL(f4), POINTER       :: ARRAYnit (:,:,:) => NULL() ! nit
      REAL(f4), POINTER       :: ARRAYnh4 (:,:,:) => NULL() ! nh4
      REAL(f4), POINTER       :: ARRAYbcpi(:,:,:) => NULL() ! bcpi
      REAL(f4), POINTER       :: ARRAYocpi(:,:,:) => NULL() ! ocpi
      REAL(f4), POINTER       :: ARRAYbcpo(:,:,:) => NULL() ! bcpo
      REAL(f4), POINTER       :: ARRAYocpo(:,:,:) => NULL() ! ocpo
      REAL(f4), POINTER       :: ARRAYdst1(:,:,:) => NULL() ! dst1

      ! Arrays to hold bulk concentration and surface area
      REAL(f4)                :: ARRAYtemp(IGLOB,JGLOB,LGLOB)
      REAL(f4)                :: ARRAYconc(IGLOB,JGLOB,LGLOB)  ! concentration

      ! Local variables
      INTEGER                 :: I, J, L, NN, LL
      CHARACTER(LEN=155)      :: FILENAME
      CHARACTER(LEN=2)        :: MON                 ! for months 10-12
      CHARACTER(LEN=1)        :: SMON                ! for months 1-9
      REAL(fpp), PARAMETER    :: TINY  = 1e-16_fpp       ! to prevent NaN
      REAL(fpp), DIMENSION(5) :: AERO_MW             ! to hold aerosol molec weights 
      REAL(fpp)               :: AIR_MW = 2.897e-02_fpp  ! avg molecular weight of air 
      REAL(fpp), DIMENSION(5) :: CU                  ! ratio of aerosol/air MW

      ! Strings
      CHARACTER(LEN=255) :: LOC = 
     &     'READ_HG2_PARTITIONING (ocean_mercury_mod.F)'

      ! Assume success
      RC        =  GIGC_SUCCESS

      ! Initialize
      
      !----------------------------------------!
      ! Molecular weights, for unit conversion !
      !----------------------------------------!

      ! aerosol molecular weights, [kg/mol]
      AERO_MW(1) = 9.6e-02_fpp    ! SO4, SO4s
      AERO_MW(2) = 1.2e-02_fpp    ! OC, BC
      AERO_MW(3) = 6.2e-02_fpp    ! NIT, NITs
      AERO_MW(4) = 1.8e-02_fpp    ! NH4
      AERO_MW(5) = 2.9e-02_fpp    ! DST 

      ! ratio of aerosol molecular weight / air molecular weight
      CU = AERO_MW / AIR_MW   !  [kg aerosol / kg air]

      !------------------------------------------------------!
      ! Put aerosol mixing ratio (mol/mol) into arrays.     
      !------------------------------------------------------!

      !---------------------------
      ! Read SO4 from HEMCO
      !---------------------------
      CALL HCO_GetPtr( am_I_Root, 'AERO_SO4', ARRAYso4, RC )
      IF ( RC /= GIGC_SUCCESS ) THEN
         CALL ERROR_STOP ( 'Cannot get pointer to AERO_SO4', LOC )
      ENDIF

      ! convert REAL*4 to REAL(fpp) 
      SO4_GC(:,:,1:LLCHEM_FIX) = ARRAYso4(:,:,1:LLCHEM_FIX)

      !---------------------------
      ! Read NH4 from HEMCO
      !---------------------------
      CALL HCO_GetPtr( am_I_Root, 'AERO_NH4', ARRAYnh4, RC )
      IF ( RC /= GIGC_SUCCESS ) THEN
         CALL ERROR_STOP ( 'Cannot get pointer to AERO_SO4', LOC )
      ENDIF

      ! convert REAL*4 to REAL(fpp) 
      NH4_CONC(:,:,1:LLCHEM_FIX) = ARRAYnh4(:,:,1:LLCHEM_FIX)

      !---------------------------
      ! Read NIT from HEMCO
      !---------------------------
      CALL HCO_GetPtr( am_I_Root, 'AERO_NIT', ARRAYnit, RC )
      IF ( RC /= GIGC_SUCCESS ) THEN
         CALL ERROR_STOP ( 'Cannot get pointer to AERO_SO4', LOC )
      ENDIF
      
      ! convert REAL*4 to REAL(fpp) 
      NIT_CONC(:,:,1:LLCHEM_FIX) = ARRAYnit(:,:,1:LLCHEM_FIX)

      !---------------------------
      ! Read BCPI+BCPO from HEMCO
      !---------------------------

      ! BCPI
      CALL HCO_GetPtr( am_I_Root, 'AERO_BCPI', ARRAYbcpi, RC )
      IF ( RC /= GIGC_SUCCESS ) THEN
         CALL ERROR_STOP ( 'Cannot get pointer to AERO_BCPI', LOC )
      ENDIF

      ! BCPO
      CALL HCO_GetPtr( am_I_Root, 'AERO_BCPO', ARRAYbcpo, RC )
      IF ( RC /= GIGC_SUCCESS ) THEN
         CALL ERROR_STOP ( 'Cannot get pointer to AERO_BCPO', LOC )
      ENDIF

      ! First get the size of ARRAYbcpi.  ARRAYbcpo is the same size
      ! since they are both stored in the same file (bmy, 3/13/15)
      LL = size( ARRAYbcpi, 3 )

      ! convert REAL*4 to REAL(fpp) 
      ARRAYtemp                 = 0.e0_f4
      ARRAYtemp(:,:,1:LL)       = ARRAYbcpi + ARRAYbcpo
      BC_CONC(:,:,1:LLCHEM_FIX) = ARRAYtemp(:,:,1:LLCHEM_FIX)

      !---------------------------
      ! Read OCPI+OCPO from HEMCO
      !---------------------------

      ! OCPI
      CALL HCO_GetPtr( am_I_Root, 'AERO_OCPI', ARRAYocpi, RC )
      IF ( RC /= GIGC_SUCCESS ) THEN
         CALL ERROR_STOP ( 'Cannot get pointer to AERO_OCPI', LOC )
      ENDIF

      ! OCPO
      CALL HCO_GetPtr( am_I_Root, 'AERO_OCPO', ARRAYocpo, RC )
      IF ( RC /= GIGC_SUCCESS ) THEN
         CALL ERROR_STOP ( 'Cannot get pointer to AERO_OCPO', LOC )
      ENDIF

      ! First get the size of ARRAYbcpi.  ARRAYbcpo is the same size
      ! since they are both stored in the same file (bmy, 3/13/15)
      LL = size( ARRAYocpi, 3 )

      ! convert REAL*4 to REAL(fpp) 
      ARRAYtemp                 = 0.e0_f4
      ARRAYtemp(:,:,1:LL)       = ARRAYocpi + ARRAYocpo
      OC_CONC(:,:,1:LLCHEM_FIX) = ARRAYtemp(:,:,1:LLCHEM_FIX)

      !---------------------------
      ! Read DST1 from HEMCO
      !---------------------------
      CALL HCO_GetPtr( am_I_Root, 'AERO_DST1', ARRAYdst1, RC )
      IF ( RC /= GIGC_SUCCESS ) THEN
         CALL ERROR_STOP ( 'Cannot get pointer to AERO_DST1', LOC )
      ENDIF

      ! convert REAL*4 to REAL(fpp) 
      DST_CONC(:,:,1:LLCHEM_FIX) = ARRAYdst1(:,:,1:LLCHEM_FIX)

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
     &                     CU(5)*(DST_CONC(I,J,L)*2e+0_fpp)           
      
        ! convert bulk aerosol mass concentration  to ug/m3
        ! Note: Even though both bpch files say the unit is 'ppb',
        !       4x5 gets read in as mol/mol and 2x2.5 gets read in
        !       as ppb, so only 4x5 needs to be scaled by 1d9.
        !       (hma, 20120104)
        ! 0..5x0.666 file is regridded from 4x5 file, yzh, 12/15/12
#if defined( GRID4x5 ) || defined( GRID05x0666 )
        BULK_CONC(I,J,L) = ( BULK_CONC(I,J,L) *
     &                     State_Met%AIRDEN(I,J,L)  )
#elif defined( GRID2x25 )
        BULK_CONC(I,J,L) = ( BULK_CONC(I,J,L) * State_Met%AIRDEN(I,J,L))
#endif 

        ! Calculate R = HgP_ads/RGM (i.e. the ratio of Hg2
        ! adsorbed onto aerosol to Hg2 in the  gas phase) 
        R(I,J,L) = BULK_CONC(I,J,L)*
     &    ( 10e+0_fpp**( ( 2.5e+3_fpp / State_Met%T(I,J,L)) 
     &    - 10e+0_fpp ))

        ! Fraction of Hg(II) in the gas phase (unitless)
        Fg(I,J,L) = 1e+0_fpp / (R(I,J,L) + 1e+0_fpp)   
            
        ! Fraction of Hg(II) in the particle phase (unitless)
        Fp(I,J,L) = 1e+0_fpp - Fg(I,J,L)
  
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

      ! Free pointers
      ARRAYso4  => NULL()
      ARRAYnit  => NULL()
      ARRAYnh4  => NULL()
      ARRAYbcpi => NULL()
      ARRAYocpi => NULL()
      ARRAYbcpo => NULL()
      ARRAYocpo => NULL()
      ARRAYdst1 => NULL()

      ! Return to calling program.
      END SUBROUTINE READ_HG2_PARTITIONING

!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: deliver_snow_hg
!
! !DESCRIPTION: Subroutine DELIVER\_SNOW\_HG delivers Hg accumulated in snow
!  to the ocean.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE DELIVER_SNOW_HG( SNOW_Hg2aq, I, J, State_Met )
!
! !USES:
!
      USE DAO_MOD,            ONLY : IS_ICE,  IS_LAND
      USE DIAG03_MOD,         ONLY : AD03,    ND03 
      USE DEPO_MERCURY_MOD,   ONLY : LHGSNOW
      USE DEPO_MERCURY_MOD,   ONLY : SNOW_HG_OC,  SNOW_HG_STORED_OC
      USE DEPO_MERCURY_MOD,   ONLY : SNOW_HG_LN,  SNOW_HG_STORED_LN
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE TRACERID_MOD,       ONLY : N_Hg_CATS

      USE CMN_SIZE_MOD      ! Size parameters
!
! !INPUT PARAMETERS: 
!
      INTEGER, INTENT(IN)        :: I, J 
      TYPE(MetState), INTENT(IN) :: State_Met   ! Meteorology State object
!
! !OUTPUT PARAMETERS:
!
      REAL(fpp), INTENT(OUT)     :: SNOW_Hg2aq(N_Hg_CATS)  ! Hg2 from snow

!
! !REMARKS:
!  Snowpack Hg is delivered to the ocean instantaneously in an ionic pulse 
!  at the start of snowmelt (i.e. when the temperature exceeds 276K). At
!  this point, both the non-reducible Hg deposited to the snowpack and the
!  remaining reducible Hg are delivered to the ocean. Over land, the snow Hg
!  reservoirs are also emptied, but this Hg is not yet added to the land.
! 
! !REVISION HISTORY: 
!  17 Jun 2011 - J. Fisher - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC  
!     
! !LOCAL VARIABLES:
!
      INTEGER       :: NN
      REAL(fpp)     :: FROCSNOW, FRLNSNOW
      LOGICAL       :: IS_MELT, IS_SNOWHG_OC, IS_SNOWHG_LN
      LOGICAL       :: IS_OPEN_OCEAN, IS_OPEN_LAND

      !=================================================================
      ! DELIVER_SNOW_HG begins here!
      !=================================================================

      ! Return to calling program if snowpack model is disabled
      IF (.NOT. LHGSNOW) RETURN

      ! Intialize
      SNOW_Hg2aq = 0d0

! FRSEAICE, FRSNO only defined in MERRA, so need to treat these differently
! in older code versions (jaf, 6/19/13)
#if   defined( MERRA ) || defined( GEOS_FP ) || defined( MERRA2 )
      IS_OPEN_OCEAN = ( ( State_Met%FROCEAN(I,J) - 
     &                    State_Met%FRSEAICE(I,J) ) > 0e+0_fpp )
      IS_OPEN_LAND  = ( ( State_Met%FRLAND(I,J)  - 
     &                    State_Met%FRSNO(I,J)    ) > 0e+0_fpp )
#elif defined( GEOS_5 )
      IS_OPEN_OCEAN = ( State_Met%FROCEAN(I,J) > 0e+0_fpp .AND.
     &                  .NOT. IS_ICE( I, J, State_Met ) )
      IS_OPEN_LAND  = ( IS_LAND( I, J, State_Met ) .AND.
     &                  ( State_Met%SNOMAS(I,J) < 10e+0_fpp ) )
#else
      IS_OPEN_OCEAN = ( .NOT. IS_ICE( I, J, State_Met ) )
      IS_OPEN_LAND  = ( IS_LAND( I, J, State_Met) .AND.
     &                  ( State_Met%SNOW(I,J) < 10e+0_fpp ) )
#endif

      IS_MELT = ( State_Met%TS(I,J) >= 276e+0_fpp )

      DO NN = 1, N_Hg_CATS

         IS_SNOWHG_OC = ( (SNOW_HG_OC(I,J,NN)        > 0e+0_fpp) .OR.
     &                 (SNOW_HG_STORED_OC(I,J,NN) > 0e+0_fpp) )
         IS_SNOWHG_LN = ( (SNOW_HG_LN(I,J,NN)        > 0e+0_fpp) .OR.
     &                 (SNOW_HG_STORED_LN(I,J,NN) > 0e+0_fpp) )

         ! OCEAN
         ! Check if melt conditions reached and snow in grid box
         IF ( IS_MELT .AND. IS_SNOWHG_OC .AND. IS_OPEN_OCEAN ) THEN

            ! Add all snow Hg to aqueous Hg2 reservoir
            SNOW_Hg2aq(NN) =
     &            ( SNOW_HG_OC(I,J,NN) + SNOW_HG_STORED_OC(I,J,NN) )
                  
            ! Store diagnostic of meltwater delivery to ocean
            IF ( ND03 > 0 ) AD03(I,J,19,1) = AD03(I,J,19,1) +
     &         ( SNOW_HG_OC(I,J,NN) + SNOW_HG_STORED_OC(I,J,NN) )

            ! Zero reservoirs over ocean box, since we've added it
            ! all to the ocean
            SNOW_HG_OC(I,J,NN) = 0e+0_fpp
            SNOW_HG_STORED_OC(I,J,NN) = 0e+0_fpp
               
         ENDIF ! ocean box
               
         ! Zero reservoirs if there is exposed land
         IF ( IS_MELT .AND. IS_SNOWHG_LN .AND. IS_OPEN_LAND ) THEN
              SNOW_HG_LN(I,J,NN) = 0e+0_fpp
              SNOW_HG_STORED_LN(I,J,NN) = 0e+0_fpp
         ENDIF
      
      ENDDO ! Hg Tracers
      
      END SUBROUTINE DELIVER_SNOW_HG
!EOC
!------------------------------------------------------------------------------

      SUBROUTINE OCEAN_MERCURY_FLUX( am_I_Root, Input_Opt, State_Met, 
     &                               State_Chm, FLUX,      RC         )
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
!  (1 ) FLUX (REAL(fpp)) : Flux of Hg(0) from the ocean [kg/s]
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
!  25 Jul 2014 - R. Yantosca - Now remove reference to function SFCWINDSQR
!  06 Nov 2014 - R. Yantosca - Replace TRANSFER_2D with direct casts
!  05 Mar 2015 - R. Yantosca - Add Input_Opt%RES_DIR to data path
!  12 Aug 2015 - R. Yantosca - Add support for MERRA2 meteorology
!
!  ##########################################################################
!  #####    NOTE: BINARY PUNCH INPUT IS BEING PHASED OUT.  THIS DATA    #####
!  #####    WILL EVENTUALLY BE READ IN FROM netCDF FILES VIA HEMCO!     #####
!  #####       -- Bob Yantosca (05 Mar 2015)                            #####
!  ##########################################################################
!  05 Aug 2015 - J. Fisher   - Compute ocean chemistry in all ocean boxes
!			       irrespective of ice cover and reduce radiation,
!			       Hg deposition accordingly in partial-ice boxes
!                              (orig. 28/11/11)
!  06 Aug 2015 - J. Fisher   - Added delivery of snow Hg to ocean (orig. 28/11/11)
!  12 Aug 2015 - J. Fisher   - Added input of Hg from Arctic rivers to ocean
!                              (orig. 8/12/11)
!  14 Aug 2015 - J. Fisher   - Reduce partitioning coefficient Kd in Arctic 
!                              (orig. 3/23/12)
!  17 Aug 2015 - J. Fisher   - Add scaling of reducible HgII based on UV-B
!                              (orig. 3/24/12)
!******************************************************************************
!
      ! References to F90 modules
      USE DAO_MOD,            ONLY : IS_WATER, IS_ICE
      USE DEPO_MERCURY_MOD,   ONLY : DD_Hg2, WD_Hg2, DD_HgP, WD_HgP
      USE DIAG03_MOD,         ONLY : AD03, ND03, AD03_RIV
      USE ERROR_MOD,          ONLY : ERROR_STOP
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE GIGC_State_Chm_Mod, ONLY : ChmState
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE GRID_MOD,           ONLY : GET_AREA_M2, GET_XMID, GET_YMID 
      USE TIME_MOD
      USE TRACERID_MOD
      USE CMN_SIZE_MOD             ! Size parameters
      USE HCO_EmisList_Mod,   ONLY : HCO_GetPtr 
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
      REAL(fpp),      INTENT(OUT)   :: FLUX(IIPAR,JJPAR,N_Hg_CATS)
      INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!
! !LOCAL VARIABLES:
!
      ! Local variables
      LOGICAL, SAVE        :: FIRST = .TRUE.
      LOGICAL              :: IS_OCEAN_BOX ! jaf
      CHARACTER(LEN=255)   :: FILENAME
      INTEGER              :: I,         J,        NN, C
      INTEGER              :: N,         N_tot_oc
      INTEGER              :: NEXTMONTH, THISMONTH
      INTEGER              :: THISYEAR

      REAL(fpp)            :: A_M2,     DTSRCE,   MLDCM
      REAL(fpp)            :: CHg0aq,   CHg0,     vi,       JorgC_kg 
      REAL(fpp)            :: TC,       TK,       Kw
      REAL(fpp)            :: Sc,       ScCO2,    USQ,      MHg
      REAL(fpp)            :: Hg2_RED,  Hg2_GONE, Hg2_CONV  !HgPaq_SUNK
      REAL(fpp)            :: FRAC_L,   FRAC_O,   H,        TOTDEP
      REAL(fpp)            :: oldMLD,   XTAU,     TOTDEPall          
      REAL(fpp)            :: FUP(IIPAR,JJPAR,N_Hg_CATS)
      REAL(fpp)            :: FDOWN(IIPAR,JJPAR,N_Hg_CATS)
      REAL(fpp)            :: X,        Y,        D                  
      REAL(fpp)            :: NPP_tot,  A_ocean,  NPP_avg,  RADz     
      REAL(fpp)            :: EC         
      REAL(fpp)            :: k_red,    k_red_rad,  k_red_bio        
      REAL(fpp)            :: k_ox    
      REAL(fpp)            :: SPM,      Frac_Hg2, OC_tot_kg          
      !REAL(fpp)            :: Hgaq_tot 
      REAL(fpp)            ::Hg2aq_tot                      
      REAL(fpp)            :: C_tot,    Ze,       OC_tot,   Hg0_OX  
      REAL(fpp)            :: Kd_part, k_ox_dark ! sjs
      REAL(fpp)            :: FRAC_OPEN_OCEAN,    FRAC_OCEAN_OR_ICE
      REAL(fpp)            :: OLDFLOW,            RIVER_HG
      REAL(fpp)            :: FRAC_REDUCIBLE,     UVI_RATIO
      REAL(fpp)            :: SNOW_Hg2aq(N_Hg_CATS)

      ! Parameters
      REAL(fpp), PARAMETER :: EC_w      = 0.0145e+0_fpp       
      REAL(fpp), PARAMETER :: EC_doc    = 0.654e+0_fpp
      REAL(fpp), PARAMETER :: C_doc     = 1.5e+0_fpp 
      !REAL(fpp), PARAMETER :: k_radbase = 1.73d-6    
      !decrease photoreduction in atm & surf ocean
      ! per H. Amos 23 Sep 2011 to multiply by 0.9e+0_fpp
      ! implemented by eck 10/19/11

      REAL(fpp), PARAMETER :: k_radbase = 1.557e-6_fpp
      REAL(fpp), PARAMETER :: k_biobase = 4.1e-10_fpp    
      REAL(fpp), PARAMETER :: k_oxbase  = 6.64e-6_fpp 
      REAL(fpp), PARAMETER :: ECchla    = 31e+0_fpp     

      ! Conversion factor from [cm/h * ng/L] --> [kg/m2/s]
      REAL(fpp), PARAMETER :: TO_KGM2S = 1.0e-11_fpp / 3600e+0_fpp 

      ! Monthly Arctic river Hg concentrations (jaf, 12/8/11)
      ! These assume no river Hg flux in Nov-Apr and concentrations at 
      ! freshet (May-June) 3x higher than in summer (Leitch et al., 2009). 
      ! These values are chosen to maximize agreement with atmospheric
      ! Hg0 observations at Arctic sites (Alert, Zeppelin, Amderma).
      ! They were calculated using:
      ! Flux = CHg*Flow = CHg_may*Flow_may + CHg_jun*Flow_jun + ...
      ! If (e.g. here) CHg is 3x higher in May-June, you end up with
      ! Flux = 3CHg_jul*Flow_may + 3CHg_jul*Flow_jun + CHG_jul*Flow_jul + ...
      !      = CHg_jul * ( 3*Flow_may + 3*Flow_jun + Flow_jul + ...)
      ! You can then calculate CHg_jul and use it to calculate CHg for
      ! other months.
      ! This is in kg/km3 = ng/L
      REAL*8,  PARAMETER :: RIVER_CHg(12) = (/0d0,  0d0,  0d0,  0d0,
     &                                       80d0, 80d0, 27d0, 27d0,
     &                                       27d0, 27d0,  0d0,  0d0 /)
      REAL*8             :: A_ARCTIC_OCEAN

      ! Small numbers to avoid dividing by zero
      REAL(fpp), PARAMETER :: SMALLNUM   = 1e-32_fpp
      REAL(fpp), PARAMETER :: NPPMINNUM   = 5e-2_fpp       
      REAL(fpp), PARAMETER :: CHLMINNUM   = 1e-1_fpp       

      ! For values from Input_Opt
      LOGICAL              :: LSPLIT
      INTEGER              :: N_TRACERS
      REAL(fpp)            :: TRACER_MW_KG(Input_Opt%N_TRACERS)
      LOGICAL              :: LArcticRiv,  LKRedUV

      ! Pointers
      ! We need to define local arrays to hold corresponding values 
      ! from the Chemistry State (State_Chm) object. (mpayer, 12/6/12)
      REAL(fpp), POINTER   :: STT(:,:,:,:)

      ! Pointers to fields in the HEMCO data structure
      REAL(f4), POINTER    :: TOMS   (:,:) => NULL() ! O3   
      REAL(f4), POINTER    :: TOMS_PD(:,:) => NULL() ! present day O3
      REAL(f4), POINTER    :: TOMS_LT(:,:) => NULL() ! long-term O3

      ! Strings
      CHARACTER(LEN=255) :: 
     & LOC = 'OCEAN_MERCURY_FLUX (GeosCore/ocean_mercury_mod.f)'

      
      IF (LOCEANCOEF) THEN
         k_ox_dark = 10d0**(-5.2d0)
      ELSE
         k_ox_dark = 1d-7  
      ENDIF

      !=================================================================
      ! OCEAN_MERCURY_FLUX begins here!
      !=================================================================

      ! Copy values from Input_Opt
      LSPLIT       = Input_Opt%LSPLIT
      LArcticRiv   = Input_Opt%LArcticRiv
      LKRedUV      = Input_Opt%LKRedUV
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

      ! Get current month
      THISMONTH = GET_MONTH()

      ! Get current year to check if leap year (jaf, 8/12/11)
      THISYEAR  = GET_YEAR()

      !-----------------------------------------------
      ! Check tagged & total sums (if necessary)
      !-----------------------------------------------
      IF ( USE_CHECKS .and. LSPLIT ) THEN
         CALL CHECK_ATMOS_MERCURY( State_Chm,
     &                             'start of OCEAN_MERCURY_FLUX' )
         CALL CHECK_OCEAN_MERCURY( 'start of OCEAN_MERCURY_FLUX' )
         CALL CHECK_OCEAN_FLUXES ( 'start of OCEAN_MERCURY_FLUX' )
      ENDIF

      !-----------------------------------------------------------------
      ! Read O3 data from HEMCO
      !-----------------------------------------------------------------

      ! Now also read TOMS O3 column data
      IF ( LKRedUV ) THEN

         ! TOMS O3 columns [dobsons]
         CALL HCO_GetPtr( am_I_Root, 'TOMS_O3_COL', TOMS, RC )
         IF ( RC /= GIGC_SUCCESS ) THEN
            CALL ERROR_STOP ( 'Cannot get pointer to TOMS_O3_COL', LOC )
         ENDIF
   
         ! TOMS O3 columns, present-day [dobsons]
         CALL HCO_GetPtr( am_I_Root, 'TOMS_O3_PD', TOMS_PD, RC )
         IF ( RC /= GIGC_SUCCESS ) THEN
            CALL ERROR_STOP ( 'Cannot get pointer to TOMSPD_O3_COL', 
     &                        LOC )
         ENDIF
   
         ! TOMS O3 columns, long-term [dobsons]
         CALL HCO_GetPtr( am_I_Root, 'TOMS_O3_LT', TOMS_LT, RC )
         IF ( RC /= GIGC_SUCCESS ) THEN
            CALL ERROR_STOP ( 'Cannot get pointer to TOMSLT_O3_COL', 
     &                        LOC )
         ENDIF

      ENDIF

      !-----------------------------------------------
      ! Read monthly NPP, RADSW, MLD, UPVEL, dMLD data
      !-----------------------------------------------

      IF ( ITS_A_NEW_MONTH() ) THEN

         ! Get monthly MLD, NPP, CHL etc.                         
         CALL OCEAN_MERCURY_READ
     &      ( am_I_Root, Input_Opt, THISMONTH, THISYEAR, RC ) !jaf

      ENDIF    

      !eds 10/19/10
      IF ( FIRST ) THEN
         MLDav = MLD
         FIRST = .FALSE.
      ENDIF
      
      !-----------------------------------------------
      ! MLD and entrainment change in middle of month
      !-----------------------------------------------
      
      dMLD = dMLD1 
      dFLOW = dFLOW1 

      IF ( ITS_MIDMONTH() ) THEN
         dMLD = dMLD2
         dFLOW = dFLOW2
      ENDIF

      ! Emission timestep [s]
      DTSRCE = GET_TS_EMIS() * 60e+0_fpp

      !----------------------------------------------------------------
      ! Calculate total mean NPP (mg/m2/day) for later                                                  
      !----------------------------------------------------------------                                     
      ! Initialize values
      NPP_tot = 0e+0_fpp
      A_ocean = 0e+0_fpp
      A_ARCTIC_OCEAN = 0e+0_fpp

      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Grid box surface area [m2]
         A_M2 = GET_AREA_M2( I, J, 1 )

         ! Grid-box latitude [degrees]
         Y = GET_YMID( I, J, 1 )

         ! Use fractional land type information in MERRA (jaf, 4/26/11)
         ! FROCEAN is a constant, so to get correct ocean fraction we
         ! need to subtract the sea ice fraction.
         ! We now compute ocean chemistry for entire ocean grid box,
         ! irrespective of ice cover (jaf, 11/28/11)
#if   defined( MERRA ) || defined( GEOS_FP ) || defined( MERRA2 )
         NPP_tot = NPP_tot + NPP(I,J) * A_M2 * State_Met%FROCEAN(I,J)
         A_ocean = A_ocean + A_M2 * State_Met%FROCEAN(I,J)
         IF ( Y >= 70 ) THEN
         !IF ( (State_Met%FROCEAN(I,J) - State_Met%FRSEAICE(I,J)) > 0d0 )
         IF ( State_Met%FROCEAN(I,J) > 0d0 )
     &      A_ARCTIC_OCEAN = A_ARCTIC_OCEAN + A_M2 
     &                     * State_Met%FROCEAN(I,J)
         ENDIF
#elif defined(GEOS_5)
         NPP_tot = NPP_tot + NPP(I,J) * A_M2 * State_Met%FROCEAN(I,J)
         A_ocean = A_ocean + A_M2 * State_Met%FROCEAN(I,J)
         IF ( Y >= 70 ) THEN
            IF ( State_Met%FROCEAN(I,J) > 0d0 )
     &      A_ARCTIC_OCEAN = A_ARCTIC_OCEAN + A_M2 
     &                     * State_Met%FROCEAN(I,J)
         ENDIF
#else
         NPP_tot = NPP_tot+ NPP(I,J) * A_M2 *
     &             ( 1e+0_fpp - State_Met%FRCLND(I,J))
         A_ocean = A_ocean + A_M2 * ( 1e+0_fpp - State_Met%FRCLND(I,J) ) 
         IF ( Y >= 70 ) THEN
            IF ( State_Met%FROCEAN(I,J) > 0d0 )
     &      A_ARCTIC_OCEAN = A_ARCTIC_OCEAN + A_M2 
     &                     * ( 1e+0_fpp - State_Met%FRCLND(I,J) )
         ENDIF
#endif

      ENDDO
      ENDDO

      NPP_avg = NPP_tot / A_ocean

      !----------------------------------------------------------------
      ! Calculate Arctic river flow and Hg flux (jaf, 12/8/11)
      !---------------------------------------------------------------- 
      ! FLOW is in m3/s; RIVER_CHg in kg/km3; RIVERHG in kg/m2/s
      OLDFLOW = FLOWNOW 
      FLOWNOW = OLDFLOW + dFLOW * DTSRCE

      ! For some months / met fields, there may be no non-ice ocean
      IF ( A_ARCTIC_OCEAN == 0e+0_fpp ) THEN
         RIVER_HG = 0e+0_fpp
      ELSE
         RIVER_HG = FLOWNOW * RIVER_CHg(THISMONTH) * 1e-9_fpp /
     &              A_ARCTIC_OCEAN
      ENDIF

      ! Initialize GEOS-Chem tracer array [kg] from Chemistry State object
      ! (mpayer, 12/6/12)
      STT => State_Chm%Tracers

      ! Loop over surface boxes  
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I,   vi,   A_M2,    Hg2_RED                         )
!$OMP+PRIVATE( J,   NN,   k_ox,    OC_tot,  Hg2_CONV               )
!$OMP+PRIVATE( N,   TK,   CHg0,    k_red_bio                       )
!$OMP+PRIVATE( C,   TC,   RADz,    Hg0_OX,  k_red_rad              )
!$OMP+PRIVATE( D,   EC,   k_red,   OLDMLD,  TOTDEPall              )
!$OMP+PRIVATE( Y,   Ze,   ScCO2,   FRAC_O,  Frac_Hg2,  Hg2aq_tot   )
!$OMP+PRIVATE( H,   Kw,   MLDCM,   TOTDEP,  OC_tot_kg              )
!$OMP+PRIVATE( X,   SPM,  CHg0aq,  Hg2_GONE                        )
!$OMP+PRIVATE( Sc,  Usq,  C_tot,   JorgC_kg                        )
!$OMP+PRIVATE( IS_OCEAN_BOX                                        )
!$OMP+PRIVATE( FRAC_OPEN_OCEAN,    FRAC_OCEAN_OR_ICE,  Snow_Hg2aq  )
!$OMP+PRIVATE( FRAC_REDUCIBLE,     UVI_RATIO                       )
!$OMP+PRIVATE( Kd_part                                             )
!$OMP+SCHEDULE( DYNAMIC )
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Grid box surface area [m2]
         A_M2       = GET_AREA_M2( I, J, 1 )

         ! Grid-box latitude [degrees]
         Y = GET_YMID( I, J, 1 )

         ! Initialize values
         Kw         = 0e+0_fpp
!         HgPaq_SUNK = 0e+0_fpp
         Hg2_CONV   = 0e+0_fpp
         TK         = 0e+0_fpp
         TC         = 0e+0_fpp 
         JorgC_kg   = 0e+0_fpp                              
         EC         = 0e+0_fpp
         RADz       = 0e+0_fpp
         k_red      = 0e+0_fpp
         k_red_rad  = 0e+0_fpp
         k_red_bio  = 0e+0_fpp
         SPM        = 0e+0_fpp
         Frac_Hg2   = 0e+0_fpp   
         Hg2aq_tot  = 0e+0_fpp
!         Hgaq_tot   = 0e+0_fpp
         Hg2_RED    = 0e+0_fpp
         C_tot      = 0e+0_fpp
         Ze         = 0e+0_fpp
         OC_tot     = 0e+0_fpp
         OC_tot_kg  = 0e+0_fpp 
         Hg0_OX     = 0e+0_fpp
         D          = 0e+0_fpp
         TOTDEP     = 0e+0_fpp
         TOTDEPall  = 0e+0_fpp
         k_ox       = 0e+0_fpp
         SNOW_Hg2aq = 0e+0_fpp
         UVI_RATIO  = 0e+0_fpp 
         FRAC_REDUCIBLE = 0e+0_fpp

         OLDMLD     = MLDav(I,J)
         MLDav(I,J) = MLDav(I,J) + dMLD(I,J) * DTSRCE
         MLDcm      = MLDav(I,J)

         ! Add error trap to prevent new MLD from being negative
         ! (jaf, 7/6/11)
         IF (MLDcm .LT. 0e+0_fpp) MLDcm = 0e+0_fpp

         ! Get fractions of land and ocean in the grid box [unitless]
         ! Use fractional land type information in MERRA. Also make sure
         ! we do not use boxes that are mostly sea ice for consistency
         ! FROCEAN is a constant, so to get correct ocean fraction we
         ! need to subtract the sea ice fraction. Don't let the fraction
         ! be less than zero (jaf, 4/26/11)
         !---------------------------------------------------------------
         ! Updated to distinguish between the fraction of the box that is
         ! open ocean and that which is ocean at depth (i.e. surface is
         ! either ocean or ice). Unless MERRA is used, these are equal.
         ! For MERRA, we perform ocean chemistry whether or not there is
         ! ice cover, but ice will reduce both solar radiation input to
         ! ocean and mercury deposition to ocean (jaf, 11/28/11)
         !---------------------------------------------------------------
#if   defined( MERRA ) || defined( GEOS_FP ) || defined( MERRA2 )
!         FRAC_L       = State_Met%FRLAND(I,J)
!         FRAC_O       = MAX( State_Met%FROCEAN(I,J) -
!     &                       State_Met%FRSEAICE(I,J), 0e+0_fpp )
!         IS_OCEAN_BOX = ( ( FRAC_O > 0e+0_fpp                  ) .and.
!     &                   ( State_Met%SEAICE00(I,J)  > 0.5e+0_fpp ) )
         FRAC_OPEN_OCEAN   = MAX( State_Met%FROCEAN(I,J) - 
     &                            State_Met%FRSEAICE(I,J), 0e+0_fpp) 
         FRAC_OCEAN_OR_ICE = State_Met%FROCEAN(I,J)
         FRAC_O            = FRAC_OCEAN_OR_ICE 
         IS_OCEAN_BOX      = ( FRAC_OCEAN_OR_ICE > 0e+0_fpp ) 
#elif defined(GEOS_5)
! Update to make these more appropriate in GEOS-5 (jaf, 05/27/13)
! If it's a sea ice box, then ocean/ice fraction=1, open ocean frac=0
! If it's a non-sea ice ocean box, then ocean/ice fraction=1, open ocean
! frac = frac_o
! Note that GEOS-5 met includes FRLAND, FRLANDICE, and FROCEAN but not
! FRSEAICE
         FRAC_OCEAN_OR_ICE = State_Met%FROCEAN(I,J)
         FRAC_O            = FRAC_OCEAN_OR_ICE
         IS_OCEAN_BOX      = ( FRAC_OCEAN_OR_ICE > 0d0 ) 
         IF ( IS_ICE( I, J, State_Met ) ) THEN
           FRAC_OPEN_OCEAN = 0e+0_fpp
         ELSE
           FRAC_OPEN_OCEAN = 1e+0_fpp
         ENDIF
#else
!         FRAC_L       = State_Met%FRCLND(I,J)
!         FRAC_O       = 1e+0_fpp - FRAC_L
         FRAC_O            = 1e+0_fpp - State_Met%FRCLND(I,J)
         FRAC_OPEN_OCEAN   = FRAC_O
         FRAC_OCEAN_OR_ICE = FRAC_O
         IS_OCEAN_BOX = ( IS_WATER( I, J, State_Met ) )
#endif

         ! Change ocean mass due to mixed layer depth change
         ! Keep before next IF so that we adjust mass in ice-covered boxes
         CALL MLD_ADJUSTMENT( I, J, OLDMLD*1e-2_fpp, MLDcm*1e-2_fpp, 
     &                        Input_Opt, State_Met )

         !---------------------------------------------------------------
         ! Deliver snowpack Hg to ocean if snow has melted
         ! Call before IF statement so that snow Hg is zeroed in land
         ! boxes as well as ocean boxes (jaf, 11/29/11)
         CALL DELIVER_SNOW_HG( SNOW_Hg2aq, I, J, State_Met )

         ! Loop over total Hg (and ocean Hg if necessary)
         DO C = 1, N_tot_oc

            ! Get Hg category #
            IF ( C == 1 ) NN = ID_Hg_tot
            IF ( C == 2 ) NN = ID_Hg_ocn

            ! Add snowpack Hg to ocean reservoirs
            Hg2aq(I,J,NN) = Hg2aq(I,J,NN) + SNOW_Hg2aq(NN)

         ENDDO

         !===========================================================
         ! Make sure we are in an ocean box
         !===========================================================
         ! Use consistent criteria for Ocean/Land/Ice categories
         ! with snowpack and terrestrial emissions  !CDH 5/18/2010
         ! Use ocean criteria developed above (jaf, 11/29/11)
         IF ( (IS_OCEAN_BOX) .and. (MLDCM > 0.99e+0_fpp) ) THEN

            !===========================================================
            ! Reduction and oxidation coefficients
            !===========================================================    
            ! Avoid having NPP or CHL to be zero
            ! Moved to OCEAN_MERCURY_READ
!            !GanLuo+NPP(I,J) = MAX ( NPP(I,J) , NPPMINNUM )
!            NPP(I,J) = MAX ( NPP(I,J)*1.d0 , NPPMINNUM )
!
!            !GanLuo+CHL(I,J) = MAX ( CHL(I,J) , CHLMINNUM )
!            CHL(I,J) = MAX ( CHL(I,J)*1.d0 , CHLMINNUM )

            ! Light attenuation (RADz) is calculated
            EC     = (EC_w + ( EC_doc * C_doc * ( NPP(I,J) / NPP_avg ) )
     &               + ( ECchla * CHL(I,J) / 1000 ) )             

            ! NOTE: State_Met%RADSWG is net radiation at ground for all
            ! MET fields except GEOS-FP and MERRA2. For those data sets,
            ! net radiation was not processed and therefore we use
            ! incident radiation at ground. We split them up below for
            ! clarity, using the new SWGDN field for incident radiation.
            ! SWGDN is also available for MERRA but we do not use it
            ! here to preserve the legacy method (ewl, 9/23/15)
#if defined ( GEOS_FP ) || defined ( MERRA2 )
            RADz   = ( 1 / ( MLDcm * 1e-2_fpp ))
     &               * (State_Met%SWGDN(I,J) / EC )    
     &               * ( 1 - EXP( -EC * ( MLDcm * 1e-2_fpp) ) )
#else
            RADz   = ( 1 / ( MLDcm * 1e-2_fpp ))
     &               * (State_Met%RADSWG(I,J) / EC )    
     &               * ( 1 - EXP( -EC * ( MLDcm * 1e-2_fpp) ) )
#endif

            ! Reduce RADz in fractional sea ice boxes to only include
            ! the fraction of open ocean (i.e. assume zero light
            ! penetration through sea ice (jaf, 11/28/11)
            IF ( FRAC_OCEAN_OR_ICE > FRAC_OPEN_OCEAN )
     &         RADz = RADz * FRAC_OPEN_OCEAN

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
               C_tot    = 36.1e+0_fpp * CHL(I,J)**0.357e+0_fpp
            ELSE
               C_tot    = 37.7e+0_fpp * CHL(I,J)**0.615e+0_fpp
            ENDIF

            ! Calculation of the euphotic depth
            IF (C_tot > 13.65) THEN
               Ze       = 912.0e+0_fpp * C_tot**(-0.839e+0_fpp)
            ELSE
               Ze       = 426.3e+0_fpp * C_tot**(-0.547e+0_fpp)
            ENDIF

            ! Recalculation of C_tot if water is shown to be well-mixed
            IF ((Ze/(MLDcm*1e-2_fpp)) < 1) THEN
               C_tot    = 42.1e+0_fpp * CHL(I,J)**0.538e+0_fpp
            ENDIF

            !--------------------------------------------------------------
            ! Standing stock of organic carbon and total biomass
            !--------------------------------------------------------------
            ! Calculated based on C:Chl ratio of 80 (wetzel et al 2006)
            ! Stodk of organic carbon is in mgC/m2
            ! Then converting to OC_tot_kg in kg/grid

            OC_tot      = C_tot * 80.0e+0_fpp

            OC_tot_kg   = OC_tot * 1e-6_fpp * A_M2 * FRAC_O

            ! Total biomas is a proxy for SPM (mg/m3) used in Hg(II)
            ! partitioning. Calculated by multiplying the standing 
            ! stock of organic carbon with 10 (exp Bundy 2004)

            SPM = ( OC_tot * 10.0e+0_fpp / ( MLDcm * 1e-2_fpp ) ) * 1.1

            !-------------------------------------------------------------- 
            ! Hg(II) - Hg(P) partitioning coefficient
            !--------------------------------------------------------------
            ! Kd_part is based on Mason et al. 1998 and Mason &
            ! Fitzgerald 1993. (L/kg)
            ! SPM is converted to kg/L by 10E-9

            ! SPM = Suspended particulate matter (kg/L)

            ! For Arctic, log10(Kd)=5.0, see Fisher et al. 2012 (jaf, 3/23/12)
            ! For Arctic, now use log10(Kd)=4.5 to match Fisher et al. (2013)
            ! When not in the arctic and LOCEANCOEF is true, use Kd_part 
            ! and k_ox_dark from Song et al. ACP 2015, otherwise use ocean 
            ! coefficients from Soerensen et al. EST 2010 ! sjs,08/12/2015
#if defined ( LINUX_PGI )
            IF ( Y >= 70.0e+0_fpp ) THEN
               Kd_part = 31622.776601684
            ELSE IF ( LOCEANCOEF ) THEN
               Kd_part = 15848.93192461
            ELSE
               Kd_part = 316227.76601684
            ENDIF
#else
            IF ( Y >= 70.0e+0_fpp ) THEN
               Kd_part = 10e+0_fpp**(4.5e+0_fpp)
            ELSE IF ( LOCEANCOEF ) THEN
               Kd_part = 10e+0_fpp**(4.2e+0_fpp)
            ELSE
               Kd_part = 10e+0_fpp**(5.5e+0_fpp)
            ENDIF
#endif

            Frac_Hg2    = 1 / ( 1 + Kd_part * SPM * 1e-9_fpp)

            !--------------------------------------------------------------
            ! Sea surface temperature in both [K] and [C]
            !--------------------------------------------------------------
            ! where TSKIN is the temperature (K) at the ground/sea surface
            ! (Use as surrogate for SST, cap at freezing point)

            TK     = MAX( State_Met%TSKIN(I,J), 273.15e+0_fpp )

            TC     = TK - 273.15e+0_fpp

            !==============================================================
            ! Volatilisation of Hg0
            !==============================================================
            
            ! Henry's law constant (gas->liquid) [unitless] [L water/L air]  
            ! (ref: Andersson et al. 2008)

            H      = EXP( ( -2404.3e+0_fpp / TK ) + 6.92e+0_fpp )

            ! Viscosity as a function of changing temperatures
            ! (ref: Loux 2001)
            ! The paper says the viscosity is given in cP but us really P
            ! and we therefor multiply with 100 to get cP.

            vi    = ( 10**( ( 1301.0e+0_fpp / ( 998.333e+0_fpp 
     &              + 8.1855e+0_fpp 
     &              * ( TC - 20.0e+0_fpp )+ 0.00585e+0_fpp 
     &              * (TC - 20.0e+0_fpp )**2 ) ) 
     &              - 3.30233e+0_fpp ) ) * 100.0e+0_fpp      

            ! Schmidt # for Hg [unitless]                             
            ! Sc = v/D = kinematic viscosity/diffusivity
            ! (ref: Poissant et al 2000; Wilke and Chang 1995)
            ! to correct for seawater D0 is decreased by 6% as suggested
            ! by Wanninkhof (1992)

            D = 7.4e-8_fpp * sqrt( 2.26 * 18.0 ) * TK /
     &             ( ( 14.8**0.6 ) *vi )

            Sc   = ( 0.017e+0_fpp * EXP( -0.025e+0_fpp * TC ) ) / D  
            
            ! Schmidt # of CO2 [unitless] for CO2 in seawater at 20 degrees C
            ! The value is set to a constant based on other ocean studies
            ! (Gardfeld et al. 2003, Rolfhus & Fitzgerald 2004, Mason et al. 2001)

            ! Correction of the Schmidt # with temperature based on Poissant
            ! et al. (2000) (for freshwatersystems).

            ScCO2  = 644.7e+0_fpp + TC * ( -6.16e+0_fpp + TC * 
     &               ( 0.11e+0_fpp ) ) 

            ! Square of surface (actually 10m) wind speed [m2/s2]
            Usq    = State_Met%U10M(I,J)**2 + State_Met%V10M(I,J)**2

            !------------------------------------------------------
            ! Parameterizations for calculating water side mass trasfer coefficient 
            !------------------------------------------------------
            ! Mass transfer coefficient [cm/h], from Nightingale et al. 2000
            Kw     = ( 0.25e+0_fpp * Usq ) / SQRT( Sc / ScCO2 )         

            !-----------------------------------------------------
            ! Additional parameterizations:

            ! Nightinale et al. 2000 for instantanous winds

!            Kw     = ( 0.33e+0_fpp*SQRT(usq)+0.22e+0_fpp*Usq) / SQRT( Sc / ScCO2 )    

            ! Lis and Merlivat 1986
            ! Has less emphasis on windspeed as a driver for evasion
            ! Gives a less total evasion than the Nigthingale et al. 2000

!            IF (SQRT(Usq) <= 3.6e+0_fpp ) THEN                               
!               Kw = ( 0.17e+0_fpp * SQRT(Usq) * ( Sc / ScCO2 )**0.67e+0_fpp )
!            ELSE IF (SQRT(Usq) > 3.6e+0_fpp .and. SQRT(Usq) <= 13e+0_fpp ) THEN
!               Kw = ( ( 2.8e+0_fpp * SQRT(Usq))-9.6 ) * ( Sc / ScCO2 )**0.5e+0_fpp
!            ELSE  
!               Kw = ( ( 5.9e+0_fpp * SQRT(Usq))-49.3 ) * ( Sc / ScCO2)**0.5e+0_fpp
!            ENDIF

            ! Wanninkhof et al (1992)

!            Kw     = ( 0.31e+0_fpp * Usq ) / SQRT( Sc / ScCO2 )

            !===========================================================
            ! Particulate sinking                              
            !===========================================================
            ! HgP sinking is based on Sunderland & Mason 2007.
            ! JorgC originally in gC m-2 year-1, which is convereted 
            ! to kgC grid-1 timestep-1
            ! NPP is converted from mgC/m2/d-1 to gC/m2/year-1
            ! JorgC = 0.1 ( ( NPP * 12 )**1.77 ) *  MLD**n * M2 * Frac_O 
            !         * 10^-3 * DTSRCE / ( 365 * 24 * 60 * 60 )

            JorgC_kg  = ( ( 0.1e+0_fpp * ((( NPP(I,J) * 365) / 1000 )
     &                 **1.77) * (( MLDcm * 1e-2_fpp )**(-0.74e+0_fpp) )
     &                 * A_M2 * FRAC_O * 1e-3_fpp) / ( 365.0e+0_fpp 
     &                 * 24.0e+0_fpp * 60.0e+0_fpp * 60.0e+0_fpp ) )
     &                * DTSRCE


            !-----------------------------------------------------------
            ! Physical transport for tracers, Part II:
            ! Upward current transport (Ekman pumping)
            ! Upward mass flux is:
            ! Mass = (Vol upwelling water) * (Conc. below thermocline)
            ! Mass = (VEL * AREA * TIME  ) * (C * Molar Mass )
            !-----------------------------------------------------------
 
            ! Use CDEEPATL to scale deepwater in NAtlantic            

            IF ( UPVEL(I,J) > 0e+0_fpp ) THEN
                 
            ! Loop over total Hg (and ocean Hg if necessary)
            DO C = 1, N_tot_oc

               ! Move to start of loop for use earlier (jaf, 12/8/11)
               ! Grid-box latitude [degrees]
               !Y = GET_YMID( I, J, 1 )
         
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
     &                * ( 1e+0_fpp + UPVEL(I,J) * DTSRCE / 
     &                ( MLDcm * 1e-2_fpp ) ) 
                  
                  ! Hg2
                  Hg2aq(I,J,NN) = Hg2aq(I,J,NN) 
     &                * ( 1e+0_fpp + UPVEL(I,J) * DTSRCE / 
     &                ( MLDcm * 1e-2_fpp ) )
               
                  ! Hg particulate
                     HgPaq(I,J,NN)  = HgPaq(I,J,NN) 
     &                * ( 1e+0_fpp + UPVEL(I,J) * DTSRCE / 
     &                ( MLDcm * 1e-2_fpp ) )

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
               FLUX(I,J,NN)  = 0e+0_fpp 
               FUP(I,J,NN)   = 0e+0_fpp
               FDOWN(I,J,NN) = 0e+0_fpp

               !--------------------------------------------------------
               ! Calculate new Hg(II) mass
               !--------------------------------------------------------
               ! Add flux of Hg(II) from rivers in Arctic (jaf, 12/8/11)
               ! Convert to kg/box/timestep and add to Hg2aq
               IF ( (LArcticRiv) .AND. (Y >= 70) ) THEN
                  Hg2aq(I,J,NN) = Hg2aq(I,J,NN) + RIVER_HG *
     &                            A_M2 * DTSRCE * FRAC_O

                  IF ( ND03 .and. NN == ID_Hg_tot )
     &                AD03_riv(I,J) = AD03_riv(I,J) + RIVER_HG *
     &                                A_M2 * DTSRCE * FRAC_O
               ENDIF
                  
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

               ! When distinguishing open ocean vs ice, only deposition
               ! to ocean should be included in Hg exchange (jaf, 11/28/11)
               !TOTDEP        = TOTDEPall * FRAC_O
               TOTDEP        = TOTDEPall * FRAC_OPEN_OCEAN

               ! Add deposited Hg(II) to the Hg(II)tot ocean mass [kg]
               Hg2aq_tot     = Hg2aq(I,J,NN) + HgPaq(I,J,NN) + TOTDEP  

               Hg2aq(I,J,NN) = Hg2aq_tot * Frac_Hg2               


               ! Mass of Hg(II)  -->  Hg(0) 
               !---------------------------
               ! Only a certain percentage of Hg(II) is considered reducible

               ! Previously this value was set as 40%. Now it is determined
               ! dynamically. Based on lab studies of Hg photoreduction,
               ! UV-B radiation leads to more DGM production than UV-A. We
               ! assume that this difference reflects the ability of higher
               ! energy UV-B to break bonds that UV-A can't. In other words
               ! UV-B increases the reducible fraction of Hg(II).
               ! Changes in O3 columns are leading to changes in UV radiation
               ! especially in the UV-B portion of the spectrum. Based on
               ! realistic present-day conditions, we assume 50% of the 
               ! reducible pool (20% of total Hg(II) can be reduced by any 
               ! wavelength, but the remaining 50% by UV-B only, and this 
               ! fraction changes with time.
               !
               ! The relative change in UV is computed from the change in
               ! O3 column. O3 column data is provided from TOMS satellite
               ! data. UV is approximated using the UV Index (UVI), which
               ! is heavily weighted to the UV-B portion of the spectrum.
               ! To compute the UVI, we use the simplified analytical
               ! formula from Madronich, 2007, Photochem & Photobio.
               ! This formula does not apply to high albedo, but we use it
               ! here only for radiation entering the water column. It also
               ! assumes clear-sky, aerosol-free conditions. However, we
               ! are considering only a ratio between present day and
               ! current model time step, so we can assume on average those
               ! conditions haven't changed (or we don't know enough about
               ! how they have). Errors are induced for high solar zenith
               ! angle, but for a given grid box / date-time, we can assume
               ! no variation between years and therefore this should
               ! have minimal impact on the ratio. Finally, SZA in a
               ! given grid box at a given day/time doesn't change from
               ! year to year, so we can simplify the formula by
               ! excluding all common factors and do not need to rely on
               ! SZA (which is occasionally negative even when there is
               ! some incoming solar radiation, especially at the
               ! equinox (jaf, 2/27/12, 8/3/12)
               !
               ! Some years have missing ozone data. In this case, we
               ! replace the missing values with the longterm mean,
               ! which will be >0 unless it is a region in polar
               ! darkness. In that case, the darkness component should
               ! already be caught with the RADz term used in the
               ! reduction part of the code. (jaf, 8/2/12)
               !
               ! Use fraction=0.4 as error trapping default value
               ! to prevent catastrophically low reduction. This
               ! should basically never be used since in case of no
               ! O3 , RADSWG will also be zero and there will be no
               ! photo-reduction below. (jaf, 8/2/12)

               ! Default: 40% HgII photo-reducible
               FRAC_REDUCIBLE = 0.4e+0_fpp
               IF ( LKRedUV ) THEN

                  IF ( (TOMS(I,J) < 0) .AND. (TOMS_LT(I,J) > 0) )
     &               TOMS(I,J) = TOMS_LT(I,J)

                  IF ( (TOMS(I,J) > 0) .AND. (TOMS_PD(I,J) > 0) ) THEN 
                     UVI_RATIO = ( TOMS(I,J) / 
     &                           TOMS_PD(I,J) )**(-1.23e+0_fpp) 
                     FRAC_REDUCIBLE = 0.2e+0_fpp + 0.2e+0_fpp * 
     &                                ( UVI_RATIO )
                  ENDIF
               ENDIF

               ! Now use new FRAC_REDUCIBLE for photo-reduction only.
               ! Retain 40% for biological reduction
               !Hg2_RED       = Hg2aq(I,J,NN) * 0.4d0 * k_red * DTSRCE
               Hg2_RED = Hg2aq(I,J,NN) * DTSRCE * ( 0.4e+0_fpp * 
     &                   k_red_bio + FRAC_REDUCIBLE * k_red_rad )    

               ! Mass of Hg(0) --> Hg(II)
               Hg0_OX        = Hg0aq(I,J,NN) * k_ox * DTSRCE

               ! Amount of Hg(II) that is lost [kg]
               Hg2_GONE      = Hg2_RED - Hg0_OX                        

               ! Cap Hg2_GONE with available Hg2
               IF ( Hg2_GONE > Hg2aq(I,J,NN) ) THEN 
                  Hg2_GONE   = MIN( Hg2_GONE, Hg2aq(I,J,NN) )
               ENDIF

               IF ( (Hg2_GONE * (-1e+0_fpp)) >  Hg0aq(I,J,NN)) THEN      
                  Hg2_GONE   = (Hg0aq(I,J,NN)*(-1))   
                  !MAX (Hg2_GONE ,(Hg0aq(I,J,NN)*(-1e+0_fpp)))
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

                  HgPaq(I,J,NN)   = MAX ( HgPaq(I,J,NN) , 0e+0_fpp )    

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
               CHg0aq        = ( Hg0aq(I,J,NN) * 1e+11_fpp   ) /   
     &                         ( A_M2          * FRAC_O ) / MLDcm 
               
               ! Gas phase Hg(0) concentration: convert [kg] -> [ng/L]
               CHg0          = STT(I,J,1,N) * 1.0e+9_fpp /
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
               !FLUX(I,J,NN)  = FLUX(I,J,NN) * TO_KGM2S * A_M2 * FRAC_O 
               ! Assume fast horizontal equilibration w/in the grid box
               ! therefore no need to scale by open ocean fraction
               ! (jaf, 6/22/11)
               ! Evasive flux only if there is some open ocean in the
               ! grid box (not 100% sea ice)
               IF ( FRAC_OPEN_OCEAN > 0d0 ) THEN
                 FLUX(I,J,NN)  = FLUX(I,J,NN) * TO_KGM2S * A_M2 * FRAC_O
               ELSE
                 FLUX(I,J,NN)  = 0e+0_fp
               ENDIF


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
!              IF (FLUX(I,J,ID_Hg_tot) > 0e+0_fpp) THEN
!                 AD03(I,J,16) = AD03(I,J,16) 
!    &                           + FLUX(I,J,ID_Hg_tot) * DTSRCE
!              ELSE IF (FLUX(I,J,ID_Hg_tot) < 0e+0_fpp) THEN
!                 AD03(I,J,17) = AD03(I,J,17) 
!    &                           + ( abs(FLUX(I,J,ID_Hg_tot) * DTSRCE ))
!              ENDIF
               ! Loop over tagged tracers (eds)
               IF (FLUX(I,J,NN) > 0e+0_fpp) THEN
                  AD03(I,J,16,NN) = AD03(I,J,16,NN) 
     &                           + FLUX(I,J,NN) * DTSRCE
               ELSE IF (FLUX(I,J,NN) < 0e+0_fpp) THEN
                  AD03(I,J,17,NN) = AD03(I,J,17,NN) 
     &                           + ( abs(FLUX(I,J,NN) * DTSRCE ))
               ENDIF
               !-------------------------------------------------------

              ENDDO

               ! Total HgII/HgP dep to open ocean [kg]
               AD03(I,J,20,1) = AD03(I,J,20,1) + TOTDEP

            ENDIF

           
         !==============================================================
         ! If we are not in an ocean box, set Hg(0) flux to zero
         !==============================================================
         ELSE

            DO NN = 1, N_Hg_CATS 
               FLUX(I,J,NN) = 0e+0_fpp
               FUP(I,J,NN)=0e+0_fpp
               FDOWN(I,J,NN)=0e+0_fpp
            ENDDO               

         ENDIF 
      
         !==============================================================
         ! Zero amts of deposited Hg2 for next timestep [kg]  
         !==============================================================
         
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! Free pointers
      NULLIFY( STT     )
      NULLIFY( TOMS    )
      NULLIFY( TOMS_PD )
      NULLIFY( TOMS_LT )

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

      SUBROUTINE OCEAN_MERCURY_READ( am_I_Root, Input_Opt,
     &                               THISMONTH, THISYEAR,  RC  )
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
!  06 Nov 2014 - R. Yantosca - Replace TRANSFER_2D with direct casts
!  17 Dec 2014 - R. Yantosca - Leave time/date variables as 8-byte
!  05 Mar 2015 - R. Yantosca - Add Input_Opt%RES_DIR to data path
!  13 Mar 2015 - R. Yantosca - Now read netCDF data via HEMCO
!  23 Feb 2015 - J. Fisher   - Add Arctic MLD (orig. 28/11/11)
!  12 Aug 2015 - J. Fisher   - Added input of Hg from Arctic rivers to ocean
!                              (orig. 8/12/11)
!  17 Aug 2015 - J. Fisher   - Added Arctic NPP, Chl from Anne Soerensen
!                              (orig. 3/19/13)
!******************************************************************************
!
! !USES:
!
      USE CMN_SIZE_MOD
      USE Error_Mod,          ONLY : Error_Stop
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE HCO_EmisList_Mod,   ONLY : HCO_GetPtr
      USE TIME_MOD,	      ONLY : ITS_A_LEAPYEAR
      USE GRID_MOD,           ONLY : GET_XMID, GET_YMID 
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)  :: am_I_Root   ! Are we on the root CPU?
      TYPE(OptInput), INTENT(IN)  :: Input_Opt   ! Input Options object
      INTEGER,        INTENT(IN)  :: THISMONTH
      INTEGER,        INTENT(IN)  :: THISYEAR    
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT) :: RC          ! Success or failure?
!
! !LOCAL VARIABLES:
!
      ! Local Variables
      LOGICAL, SAVE          :: FIRST = .TRUE.
      REAL*4                 :: ARRAY(IIPAR,JJPAR,1)
      INTEGER                :: I, J, X, Y

      ! Add time variables to get correct number of seconds per month
      ! and prevent negative flow later (jaf, 12/8/11)
      INTEGER	   :: LASTMONTH, DAYS_IN_LAST_MONTH
      INTEGER	   :: NEXTMONTH, DAYS_IN_THIS_MONTH
      INTEGER	   :: M(12) = (/ 31, 28, 31, 30, 31, 30,
     &                           31, 31, 30, 31, 30, 31 /)
      REAL*8,  PARAMETER :: SEC_PER_DAY = 3.6e+3_fpp * 24e+0_fpp

      ! Moved from ocean_mercury_flux (jaf,3/20/13)
      REAL*8,  PARAMETER   :: NPPMINNUM   = 5e-2_fpp
      REAL*8,  PARAMETER   :: CHLMINNUM   = 1e-1_fpp

      ! Arctic MLD parameters
      REAL*8,  PARAMETER   :: ArcticMLD_summer = 15e+2_fpp
      REAL*8,  PARAMETER   :: ArcticMLD_other  = 20e+2_fpp

      ! Strings
      CHARACTER(LEN=255) :: 
     & LOC = 'OCEAN_MERCURY_READ (GeosCore/ocean_mercury_mod.f)'

      !=================================================================
      ! OCEAN_MERCURY_READ begins here!
      !=================================================================
     
      ! Assume success
      RC = GIGC_SUCCESS

      ! Get the last month (jaf, 12/8/11)
      IF ( THISMONTH == 1 ) THEN
         LASTMONTH = 12
      ELSE
         LASTMONTH = THISMONTH - 1
      ENDIF

      ! Get the next month
      NEXTMONTH = MOD( THISMONTH, 12 ) +1

      ! Calculate days in last month (jaf, 12/8/11)
      IF ( LASTMONTH == 2 .and. ITS_A_LEAPYEAR( THISYEAR ) ) THEN
         DAYS_IN_LAST_MONTH = M(LASTMONTH) + 1
      ELSE
         DAYS_IN_LAST_MONTH = M(LASTMONTH)
      ENDIF

      ! Calculate days in this month (jaf, 12/8/11)
      IF ( THISMONTH == 2 .and. ITS_A_LEAPYEAR( THISYEAR ) ) THEN
         DAYS_IN_THIS_MONTH = M(THISMONTH) + 1
      ELSE
         DAYS_IN_THIS_MONTH = M(THISMONTH)
      ENDIF

      !-----------------------------------------------------------------
      ! Mixed layer depth [cm]
      !-----------------------------------------------------------------
      CALL HCO_GetPtr( am_I_Root, 'OCEAN_MLD', MLD, RC )
      IF ( RC /= GIGC_SUCCESS ) THEN
         CALL ERROR_STOP ( 'Cannot get pointer to OCEAN_MLD', LOC )
      ENDIF

      ! Convert [m] to [cm]
      MLD = MLD * 100e+0_f4

      !-----------------------------------------------------------------
      ! Chl from Modis [kg/m3] --> [mg/m3]
      !-----------------------------------------------------------------
      CALL HCO_GetPtr( am_I_Root, 'OCEAN_CHLA', CHL, RC )
      IF ( RC /= GIGC_SUCCESS ) THEN
         CALL ERROR_STOP ( 'Cannot get pointer to OCEAN_CHLA', LOC )
      ENDIF

      ! Convert [kg/m3] to [mg/m3]
      CHL = CHL * 1e6_f4

      !-----------------------------------------------------------------
      ! Arctic Chl from Jin et al 2012 [mg/m3]
      !-----------------------------------------------------------------
      CALL HCO_GetPtr( am_I_Root, 'OCEAN_CHLA_A', CHL_A, RC )
      IF ( RC /= GIGC_SUCCESS ) THEN
         CALL ERROR_STOP ( 'Cannot get pointer to OCEAN_ACHLA', LOC )
      ENDIF

      !-----------------------------------------------------------------
      ! Net primary productivity [kg/m2] -> [mg/m2]
      !-----------------------------------------------------------------
      CALL HCO_GetPtr( am_I_Root, 'OCEAN_NPP', NPP, RC )
      IF ( RC /= GIGC_SUCCESS ) THEN
         CALL ERROR_STOP ( 'Cannot get pointer to OCEAN_NPP', LOC )
      ENDIF

      ! Convert [kg/m2] -> [mg/m2]
      NPP = NPP * 1e6_f4

      !-----------------------------------------------------------------
      ! Arctic Net primary productivity [mg/m2/day]
      !-----------------------------------------------------------------
      CALL HCO_GetPtr( am_I_Root, 'OCEAN_NPP_A', NPP_A, RC )
      IF ( RC /= GIGC_SUCCESS ) THEN
         CALL ERROR_STOP ( 'Cannot get pointer to OCEAN_ANPP', LOC )
      ENDIF

      !-----------------------------------------------------------------
      ! Ekman upwelling velocity [m/s]
      !-----------------------------------------------------------------
      CALL HCO_GetPtr( am_I_Root, 'OCEAN_EKMAN_V', UPVEL, RC )
      IF ( RC /= GIGC_SUCCESS ) THEN
         CALL ERROR_STOP ( 'Cannot get pointer to OCEAN_EKMAN_V', LOC )
      ENDIF

      !-----------------------------------------------------------------
      ! MLD tendency, first half of month [cm/s]
      !-----------------------------------------------------------------
      CALL HCO_GetPtr( am_I_Root, 'OCEAN_dMLD1', dMLD1, RC )
      IF ( RC /= GIGC_SUCCESS ) THEN
         CALL ERROR_STOP ( 'Cannot get pointer to dMLD1', LOC )
      ENDIF

       ! Convert [m/s] to [cm/s]
      dMLD1 = dMLD1 * 100e+0_f4

      !-----------------------------------------------------------------
      ! MLD tendency, second half of month [cm/s]
      !-----------------------------------------------------------------
      CALL HCO_GetPtr( am_I_Root, 'OCEAN_dMLD2', dMLD2, RC )
      IF ( RC /= GIGC_SUCCESS ) THEN
         CALL ERROR_STOP ( 'Cannot get pointer to dMLD2', LOC )
      ENDIF

      ! Convert [m/s] to [cm/s]
      dMLD2 = dMLD2 * 100e+0_f4

      !-----------------------------------------------------------------
      ! Overwrite fields with Arctic specific parameters
      !-----------------------------------------------------------------

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, X, Y )

      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Grid-box longitude [degrees]
         X = GET_XMID( I, J, 1 )
         ! Grid-box latitude [degrees]
         Y = GET_YMID( I, J, 1 )

          ! In central Arctic, overwrite MLD values where global MLD
          ! dataset is based on very limited information. Instead use
          ! values from climatology of Toole et al. (2010). (jaf, 11/28/11)
          IF ( (Y >= 80e+0_fpp) .OR.
     &         ( (Y >= 70e+0_fpp) .AND. ( X >= 50e+0_fpp .OR. 
     &                                    X <= 0e+0_fpp ) ) ) THEN
               SELECT CASE ( THISMONTH )
                  CASE( 7, 8, 9 ) ! summer
                      MLD(I,J) = ArcticMLD_summer
                  CASE DEFAULT    ! rest of year
                      MLD(I,J) = ArcticMLD_other
               END SELECT

               ! Also set up MLD tendency in cm/s for each half month
               ! Months with change are 6-7, 9-10
               SELECT CASE ( THISMONTH )
                   CASE( 6 ) ! june
                      dMLD1(I,J) = 0e+0_fpp
                      dMLD2(I,J) = (ArcticMLD_summer - ArcticMLD_other)/
     &                             ( SEC_PER_DAY * DAYS_IN_LAST_MONTH )
                   CASE( 7 ) ! july
                      dMLD1(I,J) = (ArcticMLD_summer - ArcticMLD_other)/
     &                             ( SEC_PER_DAY * DAYS_IN_LAST_MONTH )
                      dMLD2(I,J) = 0e+0_fpp
                  CASE( 9 ) ! september
                      dMLD1(I,J) = 0e+0_fpp 
                      dMLD2(I,J) = (ArcticMLD_other - ArcticMLD_summer)/
     &                             ( SEC_PER_DAY * DAYS_IN_LAST_MONTH )
                   CASE ( 10 ) ! october
                      dMLD1(I,J) = (ArcticMLD_other - ArcticMLD_summer)/
     &                             ( SEC_PER_DAY * DAYS_IN_LAST_MONTH )
                      dMLD2(I,J) = 0e+0_fpp
                   CASE DEFAULT ! rest of year
                      dMLD1(I,J) = 0e+0_fpp
                      dMLD2(I,J) = 0e+0_fpp
               END SELECT


          ENDIF

         ! Added NPP for the Arctic from Jin et al. 2012. (anls)
         IF ( Y >= 60 ) THEN
            IF ( NPP_A(I,J) > NPP(I,J) ) NPP(I,J) = NPP_A(I,J)
            IF ( CHL_A(I,J) > CHL(I,J) ) CHL(I,J) = CHL_A(I,J)
            ! Apply interannual scaling
            IF ( Y >= 70 ) THEN
              NPP(I,J) = NPP(I,J) * NPP_SCF
              CHL(I,J) = CHL(I,J) * NPP_SCF
            ENDIF
         ENDIF

         ! Avoid having NPP or CHL to be zero
         NPP(I,J) = MAX ( NPP(I,J) , NPPMINNUM )
         CHL(I,J) = MAX ( CHL(I,J) , CHLMINNUM )

      ENDDO
      ENDDO
! OMP END PARALLEL DO

      !-------------------------------------------------
      ! River flow
      !-------------------------------------------------
      dFLOW1 = ( RIVERFLOW(THISMONTH) - RIVERFLOW(LASTMONTH) ) /
     &         ( SEC_PER_DAY * DAYS_IN_LAST_MONTH )
      dFLOW2 = ( RIVERFLOW(NEXTMONTH) - RIVERFLOW(THISMONTH) ) /
     &         ( SEC_PER_DAY * DAYS_IN_THIS_MONTH )

      ! Set FLOWNOW first time (jaf, 12/8/11)
      IF ( FIRST ) THEN
         FLOWNOW = RIVERFLOW(LASTMONTH) + dFLOW *
     &             SEC_PER_DAY * ( DAYS_IN_LAST_MONTH - 15 )
         FIRST = .FALSE.
      ENDIF

      ! Return to calling program
      END SUBROUTINE OCEAN_MERCURY_READ

! prior to 8/26/15 (ewl) 
!------------------------------------------------------------------------------
!
!      SUBROUTINE GET_MLD_FOR_NEXT_MONTH( am_I_Root, Input_Opt,
!     &                                   THISMONTH, THISYEAR, RC )
!!
!!******************************************************************************
!!  Subroutine GET_MLD_FOR_NEXT_MONTH reads the mixed-layer depth (MLD) 
!!  values for the next month. (sas, cdh, bmy, 3/28/06)
!!
!!  Arguments as Input:
!!  ============================================================================
!!  (1 ) THISMONTH (INTEGER) : Current month number (1-12)
!!  (2 ) THISYEAR (INTEGER)  : Current year 
!!
!!  NOTES:
!!  (1 ) Now read files from DATA_DIR_1x1/mercury_200511 (bmy, 3/28/06)
!!  (2 ) Now use actual number of days in given month (jaf, 7/6/11)
!!  06 Nov 2014 - R. Yantosca - Replace TRANSFER_2D with direct casts
!!  17 Dec 2014 - R. Yantosca - Leave time/date variables as 8-byte
!!  05 Mar 2015 - R. Yantosca - Add Input_Opt%RES_DIR to data path
!!
!!  ##########################################################################
!!  #####    NOTE: BINARY PUNCH INPUT IS BEING PHASED OUT.  THIS DATA    #####
!!  #####    WILL EVENTUALLY BE READ IN FROM netCDF FILES VIA HEMCO!     #####
!!  #####       -- Bob Yantosca (05 Mar 2015)                            #####
!!  ##########################################################################
!!  23 Feb 2015 - J. Fisher   - Add Arctic MLD (orig. 28/11/11)
!!******************************************************************************
!!
!      ! References to F90 modules
!      USE BPCH2_MOD,          ONLY : GET_TAU0, GET_RES_EXT, READ_BPCH2
!      USE BPCH2_MOD,          ONLY : GET_NAME_EXT
!      USE CMN_SIZE_MOD
!      USE GIGC_ErrCode_Mod
!      USE GIGC_Input_Opt_Mod, ONLY : OptInput
!      USE TIME_MOD,           ONLY : ITS_A_LEAPYEAR
!!
!! !INPUT PARAMETERS:
!!
!      LOGICAL,        INTENT(IN)  :: am_I_Root   ! Are we on the root CPU?
!      TYPE(OptInput), INTENT(IN)  :: Input_Opt   ! Input Options object
!      INTEGER,        INTENT(IN)  :: THISMONTH
!      INTEGER,        INTENT(IN)  :: THISYEAR
!!
!! !OUTPUT PARAMETERS:
!!
!      INTEGER,        INTENT(OUT) :: RC          ! Success or failure?
!!
!! !LOCAL VARIABLES:
!!
!!
!      INTEGER                :: I, J, NEXTMONTH
!      INTEGER                :: Y, X
!      REAL*4                 :: ARRAY(IIPAR,JJPAR,1)
!      REAL(f8)               :: TAU
!      CHARACTER(LEN=255)     :: FILENAME
!      INTEGER                :: DAYS_IN_MONTH
!      INTEGER                :: M(12) = (/ 31, 28, 31, 30, 31, 30,
!     &                                     31, 31, 30, 31, 30, 31 /)
!
!      !=================================================================
!      ! GET_MLD_FOR_NEXT_MONTH begins here!
!      !=================================================================
!      
!      ! Assume success
!      RC = GIGC_SUCCESS
!
!      ! Calculate days in this month (jaf, 7/6/11)
!      IF ( THISMONTH == 2 .and. ITS_A_LEAPYEAR( THISYEAR ) ) THEN
!         DAYS_IN_MONTH = M(THISMONTH) + 1
!      ELSE
!         DAYS_IN_MONTH = M(THISMONTH)
!      ENDIF
!
!      ! MLD file name
!#if defined( GRID05x0666 ) && defined( NESTED_NA )
!      ! Need to read emission file at nested-resolution, Y. Zhang 11/1/2011 
!      FILENAME = TRIM( Input_Opt%DATA_DIR )        // 
!     &           TRIM( Input_Opt%RES_DIR  )        //
!     &           'mercury_201203/MLD_DReqDT.'      // 
!     &           GET_NAME_EXT()    // '.'          //
!     &           GET_RES_EXT()   
!#else
!      FILENAME = TRIM( Input_Opt%DATA_DIR )        // 
!     &           TRIM( Input_Opt%RES_DIR  )        //
!     &           'mercury_201007/MLD_DReqDT.geos.' //
!     &            GET_RES_EXT()
!#endif
!
!      ! Echo info
!      IF ( am_I_Root ) THEN
!         WRITE( 6, 100 ) TRIM( FILENAME )
! 100     FORMAT( '     - GET_MLD_FOR_NEXT_MONTH: Reading ', a )  
!      ENDIF
!
!
!      ! TAU0 value for next month (uses year 1985)
!      TAU       = GET_TAU0( NEXTMONTH, 1, 1985 )
!
!      ! Read from disk; original units are [m]
!      CALL READ_BPCH2( FILENAME, 'BXHGHT-$',    5,  
!     &                 TAU,       IIPAR,        JJPAR,      
!     &                 1,         ARRAY(:,:,1), QUIET=.TRUE. )
!
!      ! Resize and cast to REAL(fpp)
!      newMLD = ARRAY(:,:,1)
!
!      ! Convert [m] to [cm]
!      newMLD = newMLD * 100e+0_fpp
!
!!$OMP PARALLEL DO
!!$OMP+DEFAULT( SHARED )
!!$OMP+PRIVATE( I, J, X, Y )
!      ! get rate of change of MLD; convert [cm/month] -> [cm/s] 
!      DO J = 1, JJPAR
!      DO I = 1, IIPAR
!
!         ! Grid-box longitude [degrees] (jaf, 11/28/11)
!         X = GET_XMID( I, J, 1 )
!         ! Grid-box latitude [degrees] (jaf, 11/28/11)
!         Y = GET_YMID( I, J, 1 )
!
!         !----------------------------------------------------------------
!         ! In central Arctic, overwrite MLD values where global MLD
!         ! dataset is based on very limited information. Instead use
!         ! values from climatology of Toole et al. (2010). (jaf, 11/28/11)
!         !----------------------------------------------------------------
!          IF ( (Y >= 80d0) .OR.
!     &         ( (Y >= 70d0) .AND. (X >= 50d0 .OR. X <= 0d0) ) ) THEN
!
!               SELECT CASE ( NEXTMONTH )
!                  CASE( 7, 8, 9 ) ! summer
!                      newMLD(I,J) = 15d2
!                  CASE DEFAULT    ! rest of year
!                      newMLD(I,J) = 20d2
!               END SELECT
!          ENDIF
!         !----------------------------------------------------------------
!
!         ! This calculation can result in negative values of MLD at the
!         ! end of the the month when there are more than 30.5 days in
!         ! the month. Now use actual number of days per month (jaf, 
!         ! 7/6/11)
!         !dMLD(I,J) = (newMLD(I,J) - MLD(I,J)) / ( 3.6d3 *24e+0_fpp * 30.5e+0_fpp )
!         dMLD(I,J) = (newMLD(I,J) - MLD(I,J)) /
!     &               ( 3.6e+3_fpp *24e+0_fpp * DAYS_IN_MONTH )
!      ENDDO
!      ENDDO
!! OMP END PARALLEL DO
!
!      ! Now also get rate of change of river flow, which isn't grid box
!      ! dependent. Flow is in m3/s, so dFLOW is in m3/s2 (jaf, 12/8/11)
!      dFLOW = ( RIVERFLOW(NEXTMONTH) - RIVERFLOW(THISMONTH) ) /  
!     &        ( 3.6d3 * 24d0 * DAYS_IN_MONTH )
!
!      ! Return to calling program
!      END SUBROUTINE GET_MLD_FOR_NEXT_MONTH
!
!!------------------------------------------------------------------------------
!
!      SUBROUTINE GET_MLD_FOR_LAST_MONTH( am_I_Root, Input_Opt,
!     &                                   THISMONTH, RC         )
!!
!!******************************************************************************
!!  Subroutine GET_MLD_FOR_LAST_MONTH reads the mixed-layer depth (MLD)
!!  values for the next month. Modified from GET_MLD_FOR_NEXT_MONTH (eds 10/19/10)
!!
!!  Arguments as Input:
!!  ============================================================================
!!  (1 ) THISMONTH (INTEGER) : Current month number (1-12)
!!
!!  NOTES:
!!  (1 ) Now read files from DATA_DIR_1x1/mercury_200511 (bmy, 3/28/06)
!!  06 Nov 2014 - R. Yantosca - Replace TRANSFER_2D with direct casts
!!  17 Dec 2014 - R. Yantosca - Leave time/date variables as 8-byte
!!  05 Mar 2015 - R. Yantosca - Add Input_Opt%RES_DIR to data path
!!
!!  ##########################################################################
!!  #####    NOTE: BINARY PUNCH INPUT IS BEING PHASED OUT.  THIS DATA    #####
!!  #####    WILL EVENTUALLY BE READ IN FROM netCDF FILES VIA HEMCO!     #####
!!  #####       -- Bob Yantosca (05 Mar 2015)                            #####
!!  ##########################################################################
!!******************************************************************************
!!
!      ! References to F90 modules
!      USE BPCH2_MOD,          ONLY : GET_TAU0, GET_RES_EXT, READ_BPCH2
!      USE BPCH2_MOD,          ONLY : GET_NAME_EXT   
!      USE CMN_SIZE_MOD
!      USE GIGC_ErrCode_Mod
!      USE GIGC_Input_Opt_Mod, ONLY : OptInput
!!
!! !INPUT PARAMETERS:
!!
!      LOGICAL,        INTENT(IN)  :: am_I_Root   ! Are we on the root CPU?
!      TYPE(OptInput), INTENT(IN)  :: Input_Opt   ! Input Options object
!      INTEGER,        INTENT(IN)  :: THISMONTH
!!
!! !OUTPUT PARAMETERS:
!!
!      INTEGER,        INTENT(OUT) :: RC          ! Success or failure?
!!
!! !LOCAL VARIABLES:
!!
!      INTEGER            :: I, J, LASTMONTH
!      REAL*4             :: ARRAY(IGLOB,JGLOB,1)
!      REAL(f8)           :: TAU
!      CHARACTER(LEN=255) :: FILENAME
!
!      !=================================================================
!      ! GET_MLD_FOR_NEXT_MONTH begins here!
!      !=================================================================
!
!      ! Assume success
!      RC        =  GIGC_SUCCESS
!
!      ! MLD file name
!#if defined( GRID05x0666 ) && defined( NESTED_NA )
!      ! Need to read emission file at nested-resolution, Y. Zhang 11/1/2011 
!      FILENAME = TRIM( Input_Opt%DATA_DIR )        // 
!     &           TRIM( Input_Opt%RES_DIR  )        //
!     &           'mercury_201203/MLD_DReqDT.'      // 
!     &           GET_NAME_EXT()    // '.'          //
!     &           GET_RES_EXT()   
!#else
!      FILENAME = TRIM( Input_Opt%DATA_DIR )        //
!     &           TRIM( Input_Opt%RES_DIR  )        //
!     &           'mercury_201007/MLD_DReqDT.geos.' // GET_RES_EXT()
!#endif
!
!      ! Echo info
!      IF ( am_I_Root ) THEN
!         WRITE( 6, 100 ) TRIM( FILENAME )
! 100     FORMAT( '     - GET_MLD_FOR_LAST_MONTH: Reading ', a )
!      ENDIF
!
!      ! Get the next month
!!      LASTMONTH = MOD( (THISMONTH-2), 12 ) +1
!      IF (THISMONTH .EQ. 1) THEN
!         LASTMONTH = 12
!      ELSE 
!        LASTMONTH = THISMONTH -1
!      ENDIF
!     
!      IF ( am_I_Root ) THEN
!         PRINT *, 'LAST MONTH: '
!         PRINT *, LASTMONTH
!      ENDIF
!
!      ! TAU0 value for next month (uses year 1985)
!      TAU	= GET_TAU0( LASTMONTH, 1, 1985 )
!
!      ! Read from disk; original units are [m]
!      CALL READ_BPCH2( FILENAME, 'BXHGHT-$',    5,
!     &                 TAU,	  IGLOB,        JGLOB,
!     &                 1,         ARRAY(:,:,1), QUIET=.TRUE. )
!
!      ! Resize and cast to REAL(fpp)
!      prevMLD = ARRAY(:,:,1)
!
!      ! Convert [m] to [cm]
!      prevMLD = prevMLD * 100e+0_fpp
!
!      ! get rate of change of MLD; convert [cm/month] -> [cm/s]
!   !   DO J = 1, JJPAR
!   !   DO I = 1, IIPAR
!   !      dMLD(I,J) = (MLD(I,J) - prevMLD(I,J)) / (3.6e+3_fpp 
!   !  &             * 24e+0_fpp * 30.5e+0_fpp)
!         ! set current MLD to intermediate b/w last month and this month
!   !      MLDav(I,J) = prevMLD(I,J) + 0.5e+0_fpp * (MLD(I,J) - 
!   !  &                prevMLD(I,J))
!   !   ENDDO
!   !   ENDDO
!
!      ! add Shaojie bug-fix for hot spots
!      MLDav = MLD
!      dMlD = 0.0e+0_fpp
!
!
!      ! Return to calling program
!      END SUBROUTINE GET_MLD_FOR_LAST_MONTH
!
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
!  (1 ) I      (INTEGER)  : GEOS-CHEM longitude index
!  (2 ) J      (INTEGER)  : GEOS-CHEM latitude index
!  (3 ) MLDold (REAL(fpp)): Old ocean mixed layer depth [m]
!  (4 ) MLDnew (REAL(fpp)): New ocean mixed layer depth [m]
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
!  05 Aug 2015 - J. Fisher   - Compute MLD adjustment in all ocean boxes
!			       irrespective of ice cover and reduce radiation,
!			       Hg deposition accordingly in partial-ice boxes
!                              (orig. 28/11/11)
!  12 Aug 2015 - R. Yantosca - Add support for MERRA2 meteorology
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
      REAL(fpp),      INTENT(IN)  :: MLDold, MLDnew  
      TYPE(OptInput), INTENT(IN)  :: Input_Opt   ! Input Options object
      TYPE(MetState), INTENT(IN)  :: State_Met   ! Meteorology State object

      ! Local variables
      INTEGER               :: C,    NN,     N_tot_oc    
      INTEGER               :: K, L
      REAL(fpp)             :: A_M2, DELTAH, FRAC_O,  MHg
      REAL(fpp)             :: X, Y                   !(added anls 01/05/09)

      ! For fields from Input_Opt
      LOGICAL               :: LSPLIT
      INTEGER               :: N_TRACERS
      REAL(fpp)             :: TRACER_MW_KG(Input_Opt%N_TRACERS)

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
#if   defined( MERRA ) || defined( GEOS_FP ) || defined( MERRA2 )
      ! We now apply ocean processes, including MLD adjustment, to
      ! entire ocean grid box, irrespective of sea ice cover (jaf,
      ! 11/28/11)
!      FRAC_O     = MAX( State_Met%FROCEAN(I,J) - 
!     &                  State_Met%FRSEAICE(I,J), 0d0)
      FRAC_O     = State_Met%FROCEAN(I,J)
#else
      FRAC_O = 1e+0_fpp - State_Met%FRCLND(I,J)
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
         IF ( MLDold > 0e+0_fpp ) THEN

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
!  23 Jun 2014 - R. Yantosca - Now use Input_Opt instead of directory_mod.F
!  30 May 2013 - J. Fisher   - Add land snow Hg reservoirs
!  06 Aug 2015 - J. Fisher   - Add stored snow Hg reservoirs (orig. 11/29/11)
!******************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD,          ONLY : OPEN_BPCH2_FOR_READ
      USE CMN_SIZE_MOD
      USE DEPO_MERCURY_MOD,   ONLY : SNOW_HG_OC, SNOW_HG_LN
      USE DEPO_MERCURY_MOD,   ONLY : SNOW_HG_STORED_OC
      USE DEPO_MERCURY_MOD,   ONLY : SNOW_HG_STORED_LN
      USE DEPO_MERCURY_MOD,   ONLY : CHECK_DIMENSIONS
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
      REAL(f8)              :: ZTAU0,     ZTAU1
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
      REAL(fpp)           :: TRACER_MW_KG(Input_Opt%N_TRACERS)

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
      FILENAME        = TRIM( Input_Opt%IN_HG_RST_FILE )

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
                  IF ( Hg0aq(I,J,NN) .LT. 0 ) Hg0aq(I,J,NN) = 0e+0_fpp
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
                  IF ( Hg2aq(I,J,NN) .LT. 0 ) Hg2aq(I,J,NN) = 0e+0_fpp
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

               ! Check for negative concentrations (jaf, 7/6/11)
               DO I = 1, IIPAR
               DO J = 1, JJPAR
                  IF ( HgPaq(I,J,NN) < 0 ) HgPaq(I,J,NN) = 0e+0_fpp
               ENDDO
               ENDDO

               ! Increment NCOUNT
               NCOUNT(NTRACER) = NCOUNT(NTRACER) + 1

            ENDIF

         ! Non-reducible Hg in ocean snowpack (jaf, 6/17/11)
         ELSE IF ( CATEGORY(1:8) == 'SNOW-HGN' ) THEN
               !----------
               ! Non-reducible Hg in snow
               !----------
            
               ! Get the Hg category #
               NN              = GET_Hg0_CAT( NTRACER )
               
               ! Store ocean Hg(0) in Hg0aq array
               SNOW_HG_STORED_OC(:,:,NN)   = Hg_OCEAN(:,:,1)
            
         ELSE IF ( CATEGORY(1:7) == 'SNOW-HG' ) THEN  
               !----------
               ! Hg in snow
               !----------

               ! Get the Hg category #
               NN              = GET_Hg0_CAT( NTRACER )

               ! Store ocean Hg(0) in Hg0aq array
               SNOW_HG_OC(:,:,NN)   = Hg_OCEAN(:,:,1)
               
         ! Non-reducible Hg in land  snowpack (jaf, 6/17/11)
         ELSE IF ( CATEGORY(1:8) == 'SNHGN-LN' ) THEN
               !----------
               ! Non-reducible Hg in snow on land
               !----------

               ! Get the Hg category #
               NN              = GET_Hg0_CAT( NTRACER )

               ! Store ocean Hg(0) in Hg0aq array
               SNOW_HG_STORED_LN(:,:,NN)   = Hg_OCEAN(:,:,1)

         ELSE IF ( CATEGORY(1:7) == 'SNHG-LN' ) THEN
               !----------
               ! Hg in snow
               !----------

               ! Get the Hg category #
               NN              = GET_Hg0_CAT( NTRACER )

               ! Store ocean Hg(0) in Hg0aq array
               SNOW_HG_LN(:,:,NN)   = Hg_OCEAN(:,:,1)

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

!##############################################################################
!### DEBUG: SET OCEAN RESTART VALUES TO ZERO FOR TESTING (uncomment if needed)
!###      Hg0aq   = 0.0e+0_fpp
!###      Hg2aq   = 0.0e+0_fpp
!###      HgPaq   = 0.0e+0_fpp
!###      SNOW_Hg = 0.0e+0_fpp
!##############################################################################

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
         WRITE( 6, 140 ) ID_Hg2(NN), TRACER_NAME( Id_Hg2(NN) ), 
     &                   SUM( Hg2aq(:,:,NN) ), 'kg'
      ENDDO

      ! HgP
      WRITE( 6, 140 ) 3, 'HgP       ', SUM( HgPaq ), 'kg'

      ! Snow Hg (jaf, 6/17/11)
      WRITE( 6, 140 ) 1, 'Reducible snow Hg ocean',SUM(SNOW_HG_OC), 'kg'
      WRITE( 6, 140 ) 1, 'Non-reducible snow Hg ocean',
     &                   SUM(SNOW_HG_STORED_OC), 'kg'
      WRITE( 6, 140 ) 1, 'Reducible snow Hg land', SUM(SNOW_HG_LN), 'kg'
      WRITE( 6, 140 ) 1, 'Non-reducible snow Hg land',
     &                   SUM(SNOW_HG_STORED_LN), 'kg'

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
!  23 Jun 2014 - R. Yantosca - Now use Input_Opt instead of directory_mod.F
!  17 Dec 2014 - R. Yantosca - Leave time/date variables as 8-byte
!  30 May 2013 - J. Fisher   - Add land snow Hg reservoirs
!  06 Aug 2015 - J. Fisher   - Add stored snow Hg reservoirs (orig. 11/29/11)
!******************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD
      USE DEPO_MERCURY_MOD,   ONLY : SNOW_HG_OC, SNOW_HG_LN
      USE DEPO_MERCURY_MOD,   ONLY : SNOW_HG_STORED_OC
      USE DEPO_MERCURY_MOD,   ONLY : SNOW_HG_STORED_LN
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE GRID_MOD,           ONLY : GET_XOFFSET, GET_YOFFSET
      USE TIME_MOD,           ONLY : EXPAND_DATE, GET_TAU
      USE TRACERID_MOD,       ONLY : ID_Hg_tot,   ID_Hg0
      USE TRACERID_MOD,       ONLY : ID_Hg2,      N_Hg_CATS
      USE TRACERID_MOD,       ONLY : ID_HgP !eds 9/9/10 

      USE CMN_SIZE_MOD           ! Size parameters

      ! Arguments           
      INTEGER,        INTENT(IN) :: NYMD, NHMS
      REAL(f8),       INTENT(IN) :: TAU
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
      FILENAME  = TRIM( Input_Opt%OUT_HG_RST_FILE )
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
      ! Total reducible Hg in snowpack on ocean
      !---------------------------
         DO NN = 1, N_Hg_CATS
            CATEGORY     = 'SNOW-HG'
            N            = ID_Hg0(NN)
            ARRAY(:,:,1) = SNOW_HG_OC(:,:,NN)

            CALL BPCH2( IU_FILE,   MODELNAME, LONRES,   LATRES,
     &                  HALFPOLAR, CENTER180, CATEGORY, N, 
     &                  UNIT,      TAU,       TAU,      RESERVED,
     &                  IIPAR,     JJPAR,     1,        IFIRST,
     &                  JFIRST,    LFIRST,    ARRAY(:,:,1) )
         ENDDO

      !-----------------------------------
      ! Total non-reducible Hg in snowpack on ocean
      !-----------------------------------
         DO NN = 1, N_Hg_CATS
            CATEGORY     = 'SNOW-HGN'
            N            = ID_Hg0(NN)
            ARRAY(:,:,1) = SNOW_HG_STORED_OC(:,:,NN)

            CALL BPCH2( IU_FILE,   MODELNAME, LONRES,   LATRES,
     &                  HALFPOLAR, CENTER180, CATEGORY, N,
     &                  UNIT,      TAU,       TAU,      RESERVED,
     &                  IIPAR,     JJPAR,     1,        IFIRST,
     &                  JFIRST,    LFIRST,    ARRAY(:,:,1) )
         ENDDO

      !-----------------------------------
      ! Total reducible Hg in snowpack on land
      !-----------------------------------
         DO NN = 1, N_Hg_CATS
            CATEGORY     = 'SNHG-LN'
            N            = ID_Hg0(NN)
            ARRAY(:,:,1) = SNOW_HG_LN(:,:,NN)

            CALL BPCH2( IU_FILE,   MODELNAME, LONRES,   LATRES,
     &                  HALFPOLAR, CENTER180, CATEGORY, N,
     &                  UNIT,      TAU,       TAU,      RESERVED,
     &                  IIPAR,     JJPAR,     1,        IFIRST,
     &                  JFIRST,    LFIRST,    ARRAY(:,:,1) )
         ENDDO

      !-----------------------------------
      ! Total non-reducible Hg in snowpack on land
      !-----------------------------------
         DO NN = 1, N_Hg_CATS
            CATEGORY     = 'SNHGN-LN'
            N            = ID_Hg0(NN)
            ARRAY(:,:,1) = SNOW_HG_STORED_LN(:,:,NN)

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

      USE CMN_SIZE_MOD              ! Size parameters

      ! Arguments as Input
      CHARACTER(LEN=*), INTENT(IN)    :: LOC
      TYPE(ChmState),   INTENT(INOUT) :: State_Chm   ! Chemistry State object

      ! Local variables
      LOGICAL                      :: FLAG
      INTEGER                      :: I,       J,       L
      INTEGER                      :: N,       NN
      REAL(fpp)                    :: Hg0_tot, Hg0_tag, RELERR0, ABSERR0
      REAL(fpp)                    :: Hg2_tot, Hg2_tag, RELERR2, ABSERR2
      REAL(fpp)                    :: HgP_tot, HgP_tag, RELERRP, ABSERRP

      ! Pointers
      ! We need to define local arrays to hold corresponding values 
      ! from the Chemistry State (State_Chm) object. (mpayer, 12/6/12)
      REAL(fpp), POINTER :: STT(:,:,:,:)

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
         Hg0_tot = 0e+0_fpp
         Hg0_tag = 0e+0_fpp
         RELERR0 = 0e+0_fpp
         ABSERR0 = 0e+0_fpp
         Hg2_tot = 0e+0_fpp
         Hg2_tag = 0e+0_fpp
         RELERR2 = 0e+0_fpp
         ABSERR2 = 0e+0_fpp
         HgP_tot = 0e+0_fpp
         Hgp_tag = 0e+0_fpp
         RELERRP = 0e+0_fpp
         ABSERRP = 0e+0_fpp

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
         IF ( Hg0_tot > 0e+0_fpp ) THEN
            RELERR0 = ABS( ( Hg0_tot - Hg0_tag ) / Hg0_tot )
         ELSE
            RELERR0 = -999e+0_fpp
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
         IF ( Hg2_tot > 0e+0_fpp ) THEN
            RELERR2 = ABS( ( Hg2_tot - Hg2_tag ) / Hg2_tot )
         ELSE
            RELERR2 = -999e+0_fpp
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
         IF ( HgP_tot > 0e+0_fpp ) THEN
            RELERRP = ABS( ( HgP_tot - HgP_tag ) / HgP_tot )
         ELSE
            RELERRP = -999e+0_fpp
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
      USE TRACERID_MOD,        ONLY : ID_Hg_tot, N_Hg_CATS

      USE CMN_SIZE_MOD              ! Size parameters

      ! Arguments
      CHARACTER(LEN=*), INTENT(IN) :: LOC

      ! Local variables
      LOGICAL, SAVE                :: FIRST = .TRUE.
      LOGICAL                      :: FLAG
      INTEGER                      :: I,       J
      REAL(fpp)                    :: Hg0_tot, Hg0_tag, RELERR0, ABSERR0
      REAL(fpp)                    :: Hg2_tot, Hg2_tag, RELERR2, ABSERR2

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
         IF ( Hg0_tot > 0e+0_fpp ) THEN
            RELERR0 = ABS( ( Hg0_tot - Hg0_tag ) / Hg0_tot )
         ELSE
            RELERR0 = -999e+0_fpp
         ENDIF

         !--------------------------------------
         ! Relative and absolute errors for Hg2
         !--------------------------------------
         Hg2_tot = Hg2aq(I,J,ID_Hg_tot)
         Hg2_tag = SUM( Hg2aq(I,J,2:N_Hg_CATS) )
         ABSERR2 = ABS( Hg2_tot - Hg2_tag )

         ! Avoid div by zero
         IF ( Hg2_tot > 0e+0_fpp ) THEN
            RELERR2 = ABS( ( Hg2_tot - Hg2_tag ) / Hg2_tot )
         ELSE
            RELERR2 = -999e+0_fpp
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
      USE CMN_SIZE_MOD
      USE ERROR_MOD,           ONLY : ERROR_STOP
      USE TRACERID_MOD,        ONLY : ID_Hg_tot, N_Hg_CATS
      USE DEPO_MERCURY_MOD,    ONLY : DD_Hg2, WD_Hg2, DD_HgP, WD_HgP

      ! Arguments
      CHARACTER(LEN=*), INTENT(IN) :: LOC

      ! Local variables
      LOGICAL                      :: FLAG
      INTEGER                      :: I,         J
      REAL(fpp)                    :: DD_tot,    DD_tag 
      REAL(fpp)                    :: DD_RELERR, DD_ABSERR      
      REAL(fpp)                    :: WD_tot,    WD_tag
      REAL(fpp)                    :: WD_RELERR, WD_ABSERR

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
         IF ( DD_tot > 0e+0_fpp ) THEN
            DD_RELERR = ABS( ( DD_tot - DD_tag ) / DD_tot )
         ELSE
            DD_RELERR = -999e+0_fpp
         ENDIF

         !---------------------------------------
         ! Absolute & relative errors for WD_Hg2
         !---------------------------------------
         WD_tot    = WD_Hg2(I,J,1)
         WD_tag    = SUM( WD_Hg2(I,J,2:N_Hg_CATS) )
         WD_ABSERR = ABS( WD_tot - WD_tag )

         ! Avoid div by zero
         IF ( WD_tot > 0e+0_fpp ) THEN
            WD_RELERR = ABS( ( WD_tot - WD_tag ) / WD_tot )
         ELSE
            WD_RELERR = -999e+0_fpp
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
!  (1 ) FLUX (REAL(fpp)): Flux array (output of OCEAN_MERCURY_FLUX)
!  (2 ) LOC (CHARACTER) : Name of routine where CHECK_FLUX_OUT is called
!
!  NOTES:
!******************************************************************************
!
      ! References to F90 modules
      USE ERROR_MOD,           ONLY : ERROR_STOP
      USE TRACERID_MOD,        ONLY : ID_Hg_tot, N_Hg_CATS

      USE CMN_SIZE_MOD              ! Size parameters

      ! Arguments
      REAL(fpp),        INTENT(IN) :: FLUX(IIPAR,JJPAR,N_Hg_CATS)       
      CHARACTER(LEN=*), INTENT(IN) :: LOC

      ! Local variables
      LOGICAL                      :: FLAG
      INTEGER                      :: I,          J
      REAL(fpp)                    :: FLX_tot,    FLX_tag
      REAL(fpp)                    :: FLX_RELERR, FLX_ABSERR

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
         IF ( FLX_tot > 0e+0_fpp ) THEN
            FLX_RELERR = ABS( ( FLX_tot - FLX_tag ) / FLX_tot )
         ELSE
            FLX_RELERR = -999e+0_fpp
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

      SUBROUTINE INIT_OCEAN_MERCURY( am_I_Root, Input_Opt, RC )
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
!  23 Feb 2015 - J. Fisher   - change intermediate water Hg2/HgP speciation in
!			        Arctic so Hg2 dominates over HgP (orig. 28/11/11)
!  12 Aug 2015 - J. Fisher   - Added interannual river flow variability
!                              (orig. 23/3/12)
!  17 Aug 2015 - J. Fisher   - Added interannual Arctic NPP variability
!                              (orig. 22/3/13)
!  27 Aug 2015 - E. Lundgren - Hard-code interannual scaling of riverflow
!                              and NPP to avoid reading ASCII
!******************************************************************************
!
! !USES:
!
      USE CMN_SIZE_MOD
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE ERROR_MOD,          ONLY : ALLOC_ERR
      USE TRACERID_MOD,       ONLY : N_Hg_CATS
      USE FILE_MOD,      ONLY : IOERROR,   IU_FILE
      USE TIME_MOD,      ONLY : GET_YEAR,  EXPAND_DATE

      ! For riverflow and NPP interannual scaling (jaf, 3/26/12)
      INTEGER               :: THISYEAR, THISYEARINDEX
      REAL(fpp), DIMENSION(31) :: RSCF_05, RSCF_06, RSCF_07 
      REAL(fpp), DIMENSION(31) :: RSCF_08, RSCF_09, RSCF_10 
      REAL(fpp), DIMENSION(31) :: NSCF
!      REAL(fpp), DIMENSION(12) :: RIVERFLOW
!      REAL(fpp)                :: NPP_SCF

      !=================================================================
      ! INIT_OCEAN_MERCURY begins here!
      !=================================================================
! end

!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)    :: am_I_Root   ! Are we on the root CPU?
      TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
 
      !=================================================================
      ! INIT_OCEAN_MERCURY begins here!
      !=================================================================

      !-----------------------------------------------------------------
      ! Allocate these arrays regardless of whether you are using
      ! the dynamic ocean.  These are needed for Hg2 partitioning.
      ! (bmy, 3/30/15)
      !-----------------------------------------------------------------

      ALLOCATE( BULK_CONC( IIPAR, JJPAR, LLPAR ), STAT=RC )
      IF ( RC /= 0 ) CALL ALLOC_ERR( 'BULK_CONC' )
      BULK_CONC = 0e+0_fpp  

      ALLOCATE( SO4_GC( IIPAR, JJPAR, LLPAR ), STAT=RC )
      IF ( RC /= 0 ) CALL ALLOC_ERR( 'SO4_GC' )
      SO4_GC = 0e+0_fpp 

      ALLOCATE( NIT_CONC( IIPAR, JJPAR, LLPAR ), STAT=RC )
      IF ( RC /= 0 ) CALL ALLOC_ERR( 'NIT_CONC' )
      NIT_CONC = 0e+0_fpp 

      ALLOCATE( NH4_CONC( IIPAR, JJPAR, LLPAR ), STAT=RC )
      IF ( RC /= 0 ) CALL ALLOC_ERR( 'NH4_CONC' )
      NH4_CONC = 0e+0_fpp 

      ALLOCATE( OC_CONC( IIPAR, JJPAR, LLPAR ), STAT=RC )
      IF ( RC /= 0 ) CALL ALLOC_ERR( 'OC_CONC' )
      OC_CONC = 0e+0_fpp 

      ALLOCATE( BC_CONC( IIPAR, JJPAR, LLPAR ), STAT=RC )
      IF ( RC /= 0 ) CALL ALLOC_ERR( 'BC_CONC' )
      BC_CONC = 0e+0_fpp 

      ALLOCATE( DST_CONC( IIPAR, JJPAR, LLPAR ), STAT=RC )
      IF ( RC /= 0 ) CALL ALLOC_ERR( 'DST_CONC' )
      DST_CONC = 0e+0_fpp 

      ALLOCATE( R( IIPAR, JJPAR, LLPAR ), STAT=RC )
      IF ( RC /= 0 ) CALL ALLOC_ERR( 'R' )
      R = 0e+0_fpp 

      ALLOCATE( SO4_CONC( IIPAR, JJPAR, LLPAR ), STAT=RC )
      IF ( RC /= 0 ) CALL ALLOC_ERR( 'SO4_CONC' )
      SO4_CONC = 0e+0_fpp 

      ALLOCATE( Fg( IIPAR, JJPAR, LLPAR ), STAT=RC )
      IF ( RC /= 0 ) CALL ALLOC_ERR( 'Fg' )
      Fg = 0e+0_fpp    

      ALLOCATE( Fp( IIPAR, JJPAR, LLPAR ), STAT=RC )
      IF ( RC /= 0 ) CALL ALLOC_ERR( 'Fp' )
      Fp = 0e+0_fpp  
     
      !-----------------------------------------------------------------
      ! Exit if you are not using the dynamic ocean
      !-----------------------------------------------------------------
      IF ( .not. Input_Opt%LDYNOCEAN ) RETURN
      
      !-----------------------------------------------------------------
      ! These variables are only needed for the dynamic ocean option
      !-----------------------------------------------------------------

      ! Turn on error checks for tagged & total sums?
      USE_CHECKS  = Input_Opt%USE_CHECKS

      ! Get current year for scale factor application
      THISYEAR = GET_YEAR()

      ! Set up concentrations of Hg(0), Hg(II), Hg(C) in deep ocean REDALERT
      IF ( Input_Opt%LPREINDHG) THEN
         CDEEP    = (/ 2e-11_fpp, 1.67e-10_fpp, 1.67e-10_fpp /)
         CDEEPATL = (/ 2e-11_fpp, 1.67e-10_fpp, 1.67e-10_fpp /)
         CDEEPNAT = (/ 2e-11_fpp, 1.67e-10_fpp, 1.67e-10_fpp /)
         CDEEPSAT = (/ 2e-11_fpp, 1.67e-10_fpp, 1.67e-10_fpp /)         
         CDEEPANT = (/ 2e-11_fpp, 1.67e-10_fpp, 1.67e-10_fpp /)
         CDEEPARC = (/ 2e-11_fpp, 1.67e-10_fpp, 1.67e-10_fpp /)
         CDEEPNPA = (/ 2e-11_fpp, 1.67e-10_fpp, 1.67e-10_fpp /)
      ELSE
         CDEEP    = (/ 1.0e-10_fpp, 4.0e-10_fpp, 4.0e-10_fpp /)
         CDEEPATL = (/ 1.4e-10_fpp, 9.3e-10_fpp, 9.3e-10_fpp /)  
         CDEEPNAT = (/ 1.5e-10_fpp, 8.2e-10_fpp, 8.2e-10_fpp /)
!         CDEEPSAT = (/ 1.0e-10_fpp, 5.0e-10_fpp, 5.0e-10_fpp /)   !(anls,100114)
         CDEEPSAT = (/0.8e-10_fpp,4.1e-10_fpp,4.1e-10_fpp/)  !eck, 10/19/11
         !reduce intermediate water mercury concentration in
         !South Atlantic Ocean to 0.9pM total
         !ref: low end of uncertainty range Sunderland and Mason '07
         CDEEPANT = (/ 1.0e-10_fpp, 3.2e-10_fpp, 3.2e-10_fpp /)
         ! Redistribute Hg fractions so Hg2 /= HgP (jaf, 11/28/11)
         !CDEEPARC = (/ 1.2e-10_fpp, 7.5e-10_fpp, 7.5e-10_fpp /)
         CDEEPARC = (/ 1.2e-10_fpp, 1.2e-9_fpp,  2.8e-10_fpp /)
         CDEEPNPA = (/ 1.0e-10_fpp, 6.0e-10_fpp, 6.0e-10_fpp /)        
      ENDIF


      ! Set monthly Arctic river flow rate in m3/s (jaf, 12/8/11)
      ! These come from the 8 largest Arctic rivers, with data
      ! from UNH/GRDC
      RIVERFLOW = (/ 2.2d4,  1.9d4,  1.8d4, 1.9d4, 9.6d4, 24.6d4,
     &               14.4d4, 10.2d4, 8.4d4, 6.2d4, 3.0d4, 2.3d4   /)

      ! Initialize default interannual NPP scaling factor (jaf, 3/19/13)
      NPP_SCF = 1e+0_fpp

      !-----------------------------------------------------------------
      ! Apply riverflow and set NPP scale factors if 2009 or prior
      !-----------------------------------------------------------------

      ! Use interannual flow scaling factors for May-Oct (jaf, 3/23/11)
      ! NOTE: scale factors are hard-coded to avoid reading ascii (ewl, 8/27/15)
      IF ( THISYEAR <= 2009) THEN

         ! May 1979-2009 river scale factors
         RSCF_05 =
     &    (/ 0.63, 0.61, 1.12, 1.09, 0.50, 1.10, 0.68, 0.69, 0.59, 0.92, 
     &       1.11, 1.61, 1.02, 1.04, 0.59, 0.84, 1.13, 0.55, 2.03, 0.43,
     &       1.46, 1.54, 1.13, 1.34, 1.00, 0.58, 1.31, 0.57, 1.58, 1.20,
     &       1.00  /)  
 
         ! Jun 1979-2009 river scale factors
         RSCF_06 =
     &    (/ 1.18, 1.00, 0.89, 0.93, 1.15, 1.00, 0.94, 1.02, 0.90, 1.03, 
     &       1.17, 0.94, 1.00, 1.03, 1.00, 1.12, 0.87, 1.05, 0.79, 1.12, 
     &       0.95, 0.91, 1.14, 1.05, 0.91, 0.99, 0.85, 1.01, 1.01, 1.01, 
     &       1.06  /) 
 
         ! Jul 1979-2009 river scale factors
         RSCF_07 =
     &    (/ 0.97, 1.12, 1.03, 0.88, 1.17, 0.83, 0.92, 0.90, 1.08, 1.06, 
     &       1.11, 0.89, 0.94, 1.00, 1.00, 1.02, 0.97, 1.07, 1.03, 1.22, 
     &       1.03, 0.98, 0.88, 0.96, 0.93, 1.06, 0.89, 0.99, 1.13, 0.98, 
     &       0.97  /) 
 
         ! Aug 1979-2009 river scale factors
         RSCF_08 =
     &    (/ 1.20, 1.01, 1.00, 0.85, 1.03, 0.92, 1.05, 0.96, 0.99, 1.19, 
     &       0.95, 0.78, 0.83, 0.89, 0.96, 0.97, 0.85, 0.80, 1.09, 1.13, 
     &       1.08, 1.03, 1.13, 1.14, 0.88, 0.93, 1.02, 0.90, 1.37, 1.15, 
     &       0.91  /) 
 
         ! Sep 1979-2008 river scale factors 
         RSCF_09 =
     &    (/ 1.42, 0.91, 1.02, 1.04, 1.06, 0.78, 0.95, 0.89, 0.67, 1.09, 
     &       1.03, 0.75, 0.84, 0.71, 0.89, 0.79, 0.92, 0.85, 1.20, 0.96, 
     &       0.92, 1.07, 0.97, 1.35, 0.83, 0.95, 1.13, 1.30, 1.45, 1.26,
     &       1.00  /) ! Use 1.00 for 2009 since no scale factor available
  
         ! Oct 1979-2008 river scale factors
         RSCF_10 =
     &    (/ 0.96, 0.89, 0.89, 0.80, 1.14, 0.79, 0.88, 1.04, 0.73, 1.11, 
     &       0.92, 0.90, 0.86, 0.94, 1.04, 1.04, 1.05, 1.06, 1.12, 0.90, 
     &       0.99, 1.09, 0.93, 1.14, 0.95, 1.17, 1.09, 1.29, 1.19, 1.10,
     &       1.00  /) ! Use 1.00 for 2009 since no scale factor available

         ! 1979-2009 NPP scale factors (scaled to 2009)
         NSCF =
     &     (/ 0.857, 0.807, 0.853, 0.835, 0.831, 0.857, 0.880, 0.823, 
     &        0.827, 0.827, 0.864, 0.933, 0.913, 0.823, 0.917, 0.857, 
     &        0.939, 0.798, 0.896, 0.914, 0.868, 0.927, 0.969, 0.955, 
     &        0.927, 0.894, 0.917, 1.035, 1.153, 1.073, 1.000 /) 

         ! Calculate year index for this year given factors start in 1979
         THISYEARINDEX = THISYEAR - 1978
           
         ! Apply riverflow scale factors for May thru Oct
         RIVERFLOW(5)  = RIVERFLOW(5)  * RSCF_05( THISYEARINDEX )
         RIVERFLOW(6)  = RIVERFLOW(6)  * RSCF_06( THISYEARINDEX )
         RIVERFLOW(7)  = RIVERFLOW(7)  * RSCF_07( THISYEARINDEX )
         RIVERFLOW(8)  = RIVERFLOW(8)  * RSCF_08( THISYEARINDEX )
         RIVERFLOW(9)  = RIVERFLOW(9)  * RSCF_09( THISYEARINDEX )
         RIVERFLOW(10) = RIVERFLOW(10) * RSCF_10( THISYEARINDEX )

         ! Set NPP scaling factor
         NPP_SCF = NSCF( THISYEARINDEX )

      ENDIF

      ! Allocate arrays     
      ALLOCATE( dMLD( IIPAR, JJPAR ), STAT=RC )
      IF ( RC /= 0 ) CALL ALLOC_ERR( 'dMLD' )
      dMLD = 0e+0_fpp

      ALLOCATE( Hg0aq( IIPAR, JJPAR, N_Hg_CATS ), STAT=RC )
      IF ( RC /= 0 ) CALL ALLOC_ERR( 'Hg0aq' )
      Hg0aq = 0e+0_fpp

      ALLOCATE( Hg2aq( IIPAR, JJPAR, N_Hg_CATS ), STAT=RC )
      IF ( RC /= 0 ) CALL ALLOC_ERR( 'Hg2aq' )
      Hg2aq = 0e+0_fpp

!- eds 8/31/10 ---------------------------------------------
! Added NN to HgPaq throughout module
!      ALLOCATE( HgPaq( IIPAR, JJPAR ), STAT=RC )
      ALLOCATE( HgPaq (IIPAR, JJPAR, N_Hg_CATS ), STAT=RC )
      IF ( RC /= 0 ) CALL ALLOC_ERR( 'HgPaq' )
      HgPaq = 0e+0_fpp

      !eds 5/15/12 fix
      ALLOCATE( Hgaq_tot (IIPAR, JJPAR, N_Hg_CATS ), STAT=RC )
      IF ( RC /= 0 ) CALL ALLOC_ERR( 'Hgaq_tot' )
      HgPaq = 0e+0_fpp

      ALLOCATE( HgPaq_SUNK (IIPAR, JJPAR, N_Hg_CATS ), STAT=RC )
      IF ( RC /= 0 ) CALL ALLOC_ERR( 'HgPaq_SUNK' )
      HgPaq = 0e+0_fpp

      ALLOCATE( MLDav( IIPAR, JJPAR ), STAT=RC )
      IF ( RC /= 0 ) CALL ALLOC_ERR( 'MLDav' )
      MLDav = 0e+0_fpp

      ALLOCATE( newMLD( IIPAR, JJPAR ), STAT=RC )
      IF ( RC /= 0 ) CALL ALLOC_ERR( 'newMLD' )
      newMLD = 0e+0_fpp

      !eds 10/19/10 fixing restart bug
      ALLOCATE( prevMLD( IIPAR, JJPAR ), STAT=RC )
      IF ( RC /=0 ) CALL ALLOC_ERR( 'prevMLD' )
      prevMLD = 0e+0_fpp

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
      IF ( ALLOCATED( MLDav     ) ) DEALLOCATE( MLDav     )
      IF ( ALLOCATED( newMLD    ) ) DEALLOCATE( newMLD    )
      IF ( ALLOCATED( prevMLD   ) ) DEALLOCATE( prevMLD   ) !eds 10/19/10
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

      ! Free pointers
      CHL   => NULL()
      CHL_A => NULL()
      MLD   => NULL()
      NPP   => NULL()
      NPP_A => NULL()
      UPVEL => NULL()
      dMLD1 => NULL()
      dMLD2 => NULL()

      ! Return to calling program
      END SUBROUTINE CLEANUP_OCEAN_MERCURY

!------------------------------------------------------------------------------
     
      ! End of module
      END MODULE OCEAN_MERCURY_MOD
