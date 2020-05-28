!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: ocean_mercury_mod.F90
!
! !DESCRIPTION: Module OCEAN\_MERCURY\_MOD contains variables and routines
!  needed to compute the oceanic flux of mercury.  Original code by Sarah
!  Strode at UWA/Seattle. (sas, bmy, 1/21/05, 4/17/06)
!\\
!\\
! !INTERFACE:
!
MODULE OCEAN_MERCURY_MOD
!
! !USES:
!
  USE inquireMod,    ONLY : findFreeLUN
  USE PRECISION_MOD, ONLY : fpp => fp   ! For GEOS-Chem Precision (fpp)
  USE PRECISION_MOD, ONLY : f4          !  Rename to avoid conflicts
  USE PRECISION_MOD, ONLY : f8

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: INIT_OCEAN_MERCURY
  PUBLIC :: CHECK_OCEAN_MERCURY
  PUBLIC :: CLEANUP_OCEAN_MERCURY
  PUBLIC :: OCEAN_MERCURY_FLUX
  PUBLIC :: LDYNSEASALT, LGCAPEMIS, LPOLARBR, LBRCHEM, LBROCHEM
  PUBLIC :: L_ADD_MBL_BR
  PUBLIC :: LGEIA05
  PUBLIC :: LVEGEMIS
  PUBLIC :: LRED_JNO2,   LGEOSLWC
  PUBLIC :: LHALOGENCHEM
  PUBLIC :: LHGAQCHEM
  PUBLIC :: LRED_CLOUDONLY
  PUBLIC :: LHg2HalfAerosol
  PUBLIC :: STRAT_BR_FACTOR,        LAnthroHgOnly
  PUBLIC :: LOHO3CHEM,              LnoUSAemis
  PUBLIC :: LGCBROMINE
  PUBLIC :: READ_HG2_PARTITIONING
  PUBLIC :: Fp, Fg
  PUBLIC :: LNEI2005
  PUBLIC :: LInPlume
  PUBLIC :: LOCEANCOEF
!
! !REMARKS:
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
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !PRIVATE TYPES:
!
  !============================================================================
  ! Hg_RST_FILE : Name of restart file with ocean tracers
  ! USE_CHECKS  : Flag for turning on error-checking
  ! MAX_RELERR  : Max error for total-tag error check [unitless]
  ! MAX_ABSERR  : Max abs error for total-tag err chk [unitless]
  ! MAX_FLXERR  : Max error tol for flux error check  [unitless]
  ! Hg2aq_tot   : Total Hg2 conc. in the mixed layer  [kg      ]
  ! DD_Hg2      : Array for Hg(II) dry dep'd to ocean [kg      ]
  ! Hgaq_tot    : Total Hg conc. in the mixed layer   [kg      ]
  ! Hg0aq       : Array for ocean mass of Hg(0)       [kg      ]
  ! Hg2aq       : Array for ocean mass of Hg(II)      [kg      ]
  ! HgPaq       : Array for ocean mass of HgP         [kg      ]
  ! dMLD        : Array for Change in ocean MLD       [cm      ]
  ! MLD         : Array for instantaneous ocean MLD   [cm      ]
  ! MLDav       : Array for monthly mean ocean MLD    [cm      ]
  ! newMLD      : Array for next month's ocean MLD    [cm      ]
  ! NPP         : Array for mean net primary prod.    [unitless]
  ! RAD         : Array for mean solar radiation      [W/m2    ]
  ! UPVEL       : Array for ocean upwelling velocity  [m/s     ]
  ! WD_Hg2      : Array for Hg(II) wet dep'd to ocean [kg      ]
  ! CHL         : Chl surface concentration           [mg(m3   ]
  ! CDEEPATL    : Conc. Hg0, Hg2, HgP below MLD-Atl   [pM      ]
  ! CDEEP       : Conc. of Hg0, Hg2, HgP below MLD    [pM      ]
  ! CDEEPNAT    : Conc. Hg0, Hg2, HgP below MLD-NAtl  [pM      ]
  ! CDEEPSAT    : Conc. Hg0, Hg2, HgP below MLD-SAtl  [pM      ]
  ! CDEEPANT    : Conc. Hg0, Hg2, HgP below MLD-Ant   [pM      ]
  ! CDEEPARC    : Conc. Hg0, Hg2, HgP below MLD-Arc   [pM      ]
  !============================================================================

  ! Scalars
  LOGICAL              :: USE_CHECKS

  ! Parameters
  REAL*4,  PARAMETER   :: MAX_RELERR = 5.0d-2
  REAL*4,  PARAMETER   :: MAX_ABSERR = 5.0d-3
  REAL*4,  PARAMETER   :: MAX_FLXERR = 5.0d-1

  REAL(fpp)            :: CDEEP(3)
  REAL(fpp)            :: CDEEPATL(3)
  REAL(fpp)            :: CDEEPNAT(3)
  REAL(fpp)            :: CDEEPSAT(3)
  REAL(fpp)            :: CDEEPANT(3)
  REAL(fpp)            :: CDEEPARC(3)
  REAL(fpp)            :: CDEEPNPA(3)

  ! For Arctic rivers (jaf, 12/8/11)
  REAL(fpp)            :: RIVERFLOW(12)
  REAL(fpp)            :: dFLOW
  REAL(fpp)            :: dFLOW1
  REAL(fpp)            :: dFLOW2
  REAL(fpp)            :: FLOWNOW

  ! Private arrays
  REAL(fpp),  ALLOCATABLE :: Hgaq_tot(:,:,:)
  REAL(fpp),  ALLOCATABLE :: dMLD(:,:)
  REAL(fpp),  ALLOCATABLE :: HgPaq_SUNK(:,:,:)
  REAL(fpp),  ALLOCATABLE :: MLDav(:,:)
  REAL(fpp),  ALLOCATABLE :: newMLD(:,:)
  REAL(fpp),  ALLOCATABLE :: prevMLD(:,:)
  REAL(fpp),  ALLOCATABLE :: RAD(:,:)

  ! added by hma for Hg2 partitioning
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

  ! For IAV in NPP (jaf, 3/19/13)
  REAL(fpp) :: NPP_SCF

  ! Logical switches for the mercury simulation, all of which are
  ! set in INIT_MERCURY (cdh, 9/1/09)
  LOGICAL   :: LDYNSEASALT, LGCAPEMIS, LPOLARBR, LBRCHEM, LBROCHEM
  LOGICAL   :: L_ADD_MBL_BR, LRED_CLOUDONLY
  LOGICAL   :: LGEIA05
  LOGICAL   :: LVEGEMIS 
  LOGICAL   :: LRED_JNO2,   LGEOSLWC
  LOGICAL   :: LHALOGENCHEM
  LOGICAL   :: LHGAQCHEM
  LOGICAL   :: LHg2HalfAerosol
  LOGICAL   :: LAnthroHgOnly,          LOHO3CHEM
  LOGICAL   :: LGCBROMINE
  LOGICAL   :: LnoUSAemis
  LOGICAL   :: LNEI2005, LInPlume
  LOGICAL   :: LOCEANCOEF
  REAL(fpp) :: STRAT_BR_FACTOR

  ! CDH Set this TRUE to use corrected area-flux relationship
  ! Set this to FALSE to use original Strode et al. (2007) model
  LOGICAL,   PARAMETER :: LOCEANFIX=.TRUE.
  ! CDH average ocean area per grid box: 1.67d11 m2/box
  ! used when eliminating AREA * FRAC_O
  REAL(fpp), PARAMETER :: FUDGE=1.67e+11_fpp

  ! Pointers to fields in the HEMCO data structure.
  ! These need to be declared REAL(f4), aka REAL*4.
  ! NOTE: We can declare these NULL here because
  ! these are SAVED global pointers (bmy, 4/29/16)
  REAL(f4), POINTER :: CHL   (:,:) => NULL()
  REAL(f4), POINTER :: CHL_A (:,:) => NULL()
  REAL(f4), POINTER :: MLD   (:,:) => NULL()
  REAL(f4), POINTER :: NPP   (:,:) => NULL()
  REAL(f4), POINTER :: NPP_A (:,:) => NULL()
  REAL(f4), POINTER :: UPVEL (:,:) => NULL()
  REAL(f4), POINTER :: dMLD1 (:,:) => NULL()
  REAL(f4), POINTER :: dMLD2 (:,:) => NULL()

  !-----------------------------------------------------------------
  ! Now keep local Hg tracer indexing variables instead of getting
  ! them from tracerid_mod.f.  The tracerid_mod.f module is going
  ! to be removed for the FlexChem implementation. (bmy, 4/26/16)
  !-----------------------------------------------------------------

  ! Scalars for Hg indexing
  INTEGER           :: N_Hg_CATS
  INTEGER           :: ID_Hg_atl,  ID_Hg_nat,  ID_Hg_sat
  INTEGER           :: ID_Hg_npa,  ID_Hg_arc,  ID_Hg_ant
  INTEGER           :: ID_Hg_ocn,  ID_Hg_tot

  ! Pointers for Hg indexing
  ! NOTE: We can declare these NULL here because
  ! these are SAVED global pointers (bmy, 4/29/16)
  INTEGER, POINTER  :: Hg0_Id_List(:) => NULL()
  INTEGER, POINTER  :: Hg2_Id_List(:) => NULL()
  INTEGER, POINTER  :: HgP_Id_List(:) => NULL()

CONTAINS
!EOC
!-------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_hg2_partitioning
!
! !DESCRIPTION: Subroutine READ\_HG2\_PARTITIONING calculates the fractions of
!  Hg(II) is the particle gas phases.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE READ_HG2_PARTITIONING( Input_Opt, State_Grid, State_Met, &
                                    THISMONTH, RC )
!
! !USES:
!
#ifdef BPCH_DIAG
    USE CMN_DIAG_MOD
    USE DIAG03_MOD,         ONLY : ND03
    USE TIME_MOD,           ONLY : SET_Hg2_DIAG
#endif
    USE ErrCode_Mod
    USE ERROR_MOD,          ONLY : GEOS_CHEM_STOP, ERROR_STOP
    USE HCO_INTERFACE_MOD,  ONLY : HcoState
    USE HCO_EmisList_Mod,   ONLY : HCO_GetPtr
    USE Input_Opt_Mod,      ONLY : OptInput
    USE State_Grid_Mod,     ONLY : GrdState
    USE State_Met_Mod,      ONLY : MetState
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)  :: Input_Opt   ! Input Options object
    TYPE(MetState), INTENT(IN)  :: State_Met   ! Meteorology State object
    TYPE(GrdState), INTENT(IN)  :: State_Grid  ! Grid State object
    INTEGER,        INTENT(IN)  :: THISMONTH   ! Current month
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT) :: RC          ! Success or failure?
!
! !REMARKS:
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
!  The objective of this subroutine is to determine the fraction of reactive
!  mercury in the gas-phase (Fg) and the fraction in the particle-phase (Fp).
!  That can be done easily once R is calculated, where R is the ratio of
!  particle to gas:
!
!      R  = (HgP,ads)/RGM                              (unitless)
!      Fg = 1/(R + 1)                                  (unitless)
!      Fp = 1 - Fg   *or*  = R/(R + 1)                 (unitless)
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Pointers to fields in the HEMCO data structure
    ! These have to be declared REAL(f4), aka REAL*4.
    REAL(f4), POINTER  :: ARRAYso4 (:,:,:)
    REAL(f4), POINTER  :: ARRAYnit (:,:,:)
    REAL(f4), POINTER  :: ARRAYnh4 (:,:,:)
    REAL(f4), POINTER  :: ARRAYbcpi(:,:,:)
    REAL(f4), POINTER  :: ARRAYocpi(:,:,:)
    REAL(f4), POINTER  :: ARRAYbcpo(:,:,:)
    REAL(f4), POINTER  :: ARRAYocpo(:,:,:)
    REAL(f4), POINTER  :: ARRAYdst1(:,:,:)

    ! Arrays to hold bulk concentration and surface area
    REAL(f4)           :: ARRAYtemp(State_Grid%NX,State_Grid%NY,State_Grid%NZ)
    REAL(f4)           :: ARRAYconc(State_Grid%NX,State_Grid%NY,State_Grid%NZ)

    ! Scalars
    INTEGER                 :: I, J, L, NN, LL
    CHARACTER(LEN=155)      :: FILENAME
    CHARACTER(LEN=2)        :: MON                 ! for months 10-12
    CHARACTER(LEN=1)        :: SMON                ! for months 1-9
    REAL(fpp), PARAMETER    :: TINY  = 1e-16_fpp   ! to prevent NaN
    REAL(fpp), DIMENSION(5) :: AERO_MW             ! aerosol molec weights
    REAL(fpp)               :: AIR_MW = 2.897e-02_fpp ! avg molec weight of air
    REAL(fpp), DIMENSION(5) :: CU                  ! ratio of aerosol/air MW

    ! Strings
    CHARACTER(LEN=255) :: LOC = 'READ_HG2_PARTITIONING (ocean_mercury_mod.F90)'

    !=================================================================
    ! READ_HG2_PARTITIONING begins here!
    !=================================================================

    ! Assume success
    RC        =  GC_SUCCESS

    ! Initialize pointers
    ARRAYso4  => NULL() ! so4
    ARRAYnit  => NULL() ! nit
    ARRAYnh4  => NULL() ! nh4
    ARRAYbcpi => NULL() ! bcpi
    ARRAYocpi => NULL() ! ocpi
    ARRAYbcpo => NULL() ! bcpo
    ARRAYocpo => NULL() ! ocpo
    ARRAYdst1 => NULL() ! dst1

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
    CALL HCO_GetPtr( HcoState, 'AERO_SO4', ARRAYso4, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL ERROR_STOP ( 'Cannot get pointer to AERO_SO4', LOC )
    ENDIF

    ! convert REAL*4 to REAL(fpp)
    SO4_GC(:,:,1:State_Grid%MaxChemLev) = ARRAYso4(:,:,1:State_Grid%MaxChemLev)

    !---------------------------
    ! Read NH4 from HEMCO
    !---------------------------
    CALL HCO_GetPtr( HcoState, 'AERO_NH4', ARRAYnh4, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL ERROR_STOP ( 'Cannot get pointer to AERO_SO4', LOC )
    ENDIF

    ! convert REAL*4 to REAL(fpp)
    NH4_CONC(:,:,1:State_Grid%MaxChemLev) = ARRAYnh4(:,:,1:State_Grid%MaxChemLev)

    !---------------------------
    ! Read NIT from HEMCO
    !---------------------------
    CALL HCO_GetPtr( HcoState, 'AERO_NIT', ARRAYnit, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL ERROR_STOP ( 'Cannot get pointer to AERO_SO4', LOC )
    ENDIF

    ! convert REAL*4 to REAL(fpp)
    NIT_CONC(:,:,1:State_Grid%MaxChemLev) = ARRAYnit(:,:,1:State_Grid%MaxChemLev)

    !---------------------------
    ! Read BCPI+BCPO from HEMCO
    !---------------------------

    ! BCPI
    CALL HCO_GetPtr( HcoState, 'AERO_BCPI', ARRAYbcpi, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL ERROR_STOP ( 'Cannot get pointer to AERO_BCPI', LOC )
    ENDIF

    ! BCPO
    CALL HCO_GetPtr( HcoState, 'AERO_BCPO', ARRAYbcpo, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL ERROR_STOP ( 'Cannot get pointer to AERO_BCPO', LOC )
    ENDIF

    ! First get the size of ARRAYbcpi.  ARRAYbcpo is the same size
    ! since they are both stored in the same file (bmy, 3/13/15)
    LL = size( ARRAYbcpi, 3 )

    ! convert REAL*4 to REAL(fpp)
    ARRAYtemp           = 0.e0_f4
    ARRAYtemp(:,:,1:LL) = ARRAYbcpi + ARRAYbcpo
    BC_CONC(:,:,1:State_Grid%MaxChemLev) = ARRAYtemp(:,:,1:State_Grid%MaxChemLev)

    !---------------------------
    ! Read OCPI+OCPO from HEMCO
    !---------------------------

    ! OCPI
    CALL HCO_GetPtr( HcoState, 'AERO_OCPI', ARRAYocpi, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL ERROR_STOP ( 'Cannot get pointer to AERO_OCPI', LOC )
    ENDIF

    ! OCPO
    CALL HCO_GetPtr( HcoState, 'AERO_OCPO', ARRAYocpo, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL ERROR_STOP ( 'Cannot get pointer to AERO_OCPO', LOC )
    ENDIF

    ! First get the size of ARRAYbcpi.  ARRAYbcpo is the same size
    ! since they are both stored in the same file (bmy, 3/13/15)
    LL = size( ARRAYocpi, 3 )

    ! convert REAL*4 to REAL(fpp)
    ARRAYtemp           = 0.e0_f4
    ARRAYtemp(:,:,1:LL) = ARRAYocpi + ARRAYocpo
    OC_CONC(:,:,1:State_Grid%MaxChemLev) = ARRAYtemp(:,:,1:State_Grid%MaxChemLev)

    !---------------------------
    ! Read DST1 from HEMCO
    !---------------------------
    CALL HCO_GetPtr( HcoState, 'AERO_DST1', ARRAYdst1, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL ERROR_STOP ( 'Cannot get pointer to AERO_DST1', LOC )
    ENDIF

    ! convert REAL*4 to REAL(fpp)
    DST_CONC(:,:,1:State_Grid%MaxChemLev) = ARRAYdst1(:,:,1:State_Grid%MaxChemLev)

    !$OMP PARALLEL DO       &
    !$OMP DEFAULT( SHARED ) &
    !$OMP PRIVATE( I, J, L )
    DO L = 1, State_Grid%NZ
    DO J = 1, State_Grid%NY
    DO I = 1, State_Grid%NX

       ! Sum aerosol types (ppbv) to get bulk aerosol
       !
       ! - Multiply species by CU factor as part of unit conversion from
       !    ppbv --> ug/m3
       ! Scale dust from 4x5 run (ref: Fairlie et al. (2010))
       BULK_CONC(I,J,L) = CU(1)*SO4_CONC(I,J,L)                   + &
                          CU(2)*(BC_CONC(I,J,L) + OC_CONC(I,J,L)) + &
                          CU(3)*NIT_CONC(I,J,L)                   + &
                          CU(4)*NH4_CONC(I,J,L) +                   &
                          CU(5)*(DST_CONC(I,J,L)*2e+0_fpp)

       ! convert bulk aerosol mass concentration  to ug/m3
       BULK_CONC(I,J,L) = ( BULK_CONC(I,J,L) * State_Met%AIRDEN(I,J,L)  )

       ! Calculate R = HgP_ads/RGM (i.e. the ratio of Hg2
       ! adsorbed onto aerosol to Hg2 in the  gas phase)
       R(I,J,L) = BULK_CONC(I,J,L) * &
            ( 10e+0_fpp**( ( 2.5e+3_fpp / State_Met%T(I,J,L)) - 10e+0_fpp ))

       ! Fraction of Hg(II) in the gas phase (unitless)
       Fg(I,J,L) = 1e+0_fpp / (R(I,J,L) + 1e+0_fpp)

       ! Fraction of Hg(II) in the particle phase (unitless)
       Fp(I,J,L) = 1e+0_fpp - Fg(I,J,L)

    ENDDO
    ENDDO
    ENDDO
    !$OMP END PARALLEL DO

#ifdef BPCH_DIAG
    !-------------------------------------------!
    ! Store Fg and Fp in ND03 diagnostic.       !
    !-------------------------------------------!
    !$OMP PARALLEL DO       &
    !$OMP DEFAULT( SHARED ) &
    !$OMP PRIVATE( I, J, L )
    DO L = 1, State_Grid%NZ
    DO J = 1, State_Grid%NY
    DO I = 1, State_Grid%NX

       ! error check
       IF (Fg(I,J,L)+Fp(I,J,L) > 1) THEN
          PRINT*, ' Fg + Fp > 1 @ ocean_mercury_mod.F90'
          CALL GEOS_CHEM_STOP
       ENDIF

    ENDDO
    ENDDO
    ENDDO
    !$OMP END PARALLEL DO

    ! Increment diagnostic timestep counter.
    CALL SET_Hg2_DIAG( INCREMENT=.TRUE. )
#endif

    ! Free pointers
    ARRAYso4  => NULL()
    ARRAYnit  => NULL()
    ARRAYnh4  => NULL()
    ARRAYbcpi => NULL()
    ARRAYocpi => NULL()
    ARRAYbcpo => NULL()
    ARRAYocpo => NULL()
    ARRAYdst1 => NULL()

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
  SUBROUTINE DELIVER_SNOW_HG( SNOW_Hg2aq, I, J, State_Met, &
                              State_Chm, State_Diag )
!
! !USES:
!
#ifdef BPCH_DIAG
    USE DIAG03_MOD,         ONLY : AD03,    ND03
#endif
    USE DEPO_MERCURY_MOD,   ONLY : LHGSNOW
    USE State_Chm_Mod,      ONLY : ChmState
    USE State_Diag_Mod,     ONLY : DgnState
    USE State_Met_Mod,      ONLY : MetState
    USE Time_Mod,           ONLY : Get_Ts_Chem
!
! !INPUT PARAMETERS:
!
    INTEGER, INTENT(IN)           :: I, J
    TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology State object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
    TYPE(DgnState), INTENT(INOUT) :: State_Diag  ! Diagnostics State object
!
! !OUTPUT PARAMETERS:
!
    REAL(fpp),      INTENT(OUT)   :: SNOW_Hg2aq(N_Hg_CATS)  ! Hg2 from snow
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
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER       :: NN
    REAL(fpp)     :: FROCSNOW, FRLNSNOW, DT
    LOGICAL       :: IS_MELT, IS_SNOWHG_OC, IS_SNOWHG_LN
    LOGICAL       :: IS_OPEN_OCEAN, IS_OPEN_LAND

    ! Pointers
    REAL(fpp), POINTER :: SNOW_HG_OC(:,:,:)
    REAL(fpp), POINTER :: SNOW_HG_LN(:,:,:)
    REAL(fpp), POINTER :: SNOW_HG_STORED_OC(:,:,:)
    REAL(fpp), POINTER :: SNOW_HG_STORED_LN(:,:,:)

    !=================================================================
    ! DELIVER_SNOW_HG begins here!
    !=================================================================

    ! Return to calling program if snowpack model is disabled
    IF (.NOT. LHGSNOW) RETURN

    ! Chemistry timestep in seconds
    DT = GET_TS_CHEM()

    ! Intialize
    SNOW_Hg2aq = 0d0

    ! Point to fields in State_Chm
    SNOW_HG_OC        => State_Chm%SnowHgOcean
    SNOW_HG_LN        => State_Chm%SnowHgLand
    SNOW_HG_STORED_OC => State_Chm%SnowHgOceanStored
    SNOW_HG_STORED_LN => State_Chm%SnowHgLandStored

    IS_OPEN_OCEAN = ( ( State_Met%FROCEAN(I,J) - &
                        State_Met%FRSEAICE(I,J) ) > 0e+0_fpp )
    IS_OPEN_LAND  = ( ( State_Met%FRLAND(I,J)  - &
                        State_Met%FRSNO(I,J)    ) > 0e+0_fpp )

    IS_MELT = ( State_Met%TS(I,J) >= 276e+0_fpp )

    DO NN = 1, N_Hg_CATS

       IS_SNOWHG_OC = ( (SNOW_HG_OC(I,J,NN)        > 0e+0_fpp) .OR. &
                        (SNOW_HG_STORED_OC(I,J,NN) > 0e+0_fpp) )
       IS_SNOWHG_LN = ( (SNOW_HG_LN(I,J,NN)        > 0e+0_fpp) .OR. &
                        (SNOW_HG_STORED_LN(I,J,NN) > 0e+0_fpp) )

       ! OCEAN
       ! Check if melt conditions reached and snow in grid box
       IF ( IS_MELT .AND. IS_SNOWHG_OC .AND. IS_OPEN_OCEAN ) THEN

          ! Add all snow Hg to aqueous Hg2 reservoir
          SNOW_Hg2aq(NN) = ( SNOW_HG_OC(I,J,NN) + SNOW_HG_STORED_OC(I,J,NN) )

#ifdef BPCH_DIAG
          !===========================================================
          ! %%%%% ND03 (bpch) DIAGNOSTIC %%%%%
          !
          ! Store diagnostic of meltwater delivery to ocean [kg]
          !===========================================================
          IF ( ND03 > 0 ) THEN
             AD03(I,J,19,1) = AD03(I,J,19,1) + &
                             ( SNOW_HG_OC(I,J,NN) + SNOW_HG_STORED_OC(I,J,NN) )
          ENDIF
#endif

          !===========================================================
          ! %%%%% HISTORY (aka netCDF diagnostics) %%%%%
          !
          ! Store diagnostic of meltwater delivery to ocean [kg/s]
          !===========================================================
          IF ( State_Diag%Archive_EmisHg2snowToOcean ) THEN
             State_Diag%EmisHg2snowToOcean(I,J) = &
                  ( SNOW_HG_OC(I,J,NN) + SNOW_HG_STORED_OC(I,J,NN) ) / DT
          ENDIF

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

    ! Free pointers
    SNOW_HG_OC        => NULL()
    SNOW_HG_LN        => NULL()
    SNOW_HG_STORED_OC => NULL()
    SNOW_HG_STORED_LN => NULL()

  END SUBROUTINE DELIVER_SNOW_HG
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ocean_mercury_flux
!
! !DESCRIPTION: Subroutine OCEAND\_MERCURY\_FLUX calculates emissions of Hg(0)
!  the ocean in [kg/s].  (sas, bmy, 1/19/05, 4/17/06)
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE OCEAN_MERCURY_FLUX( Input_Opt, State_Chm, State_Diag, State_Grid, &
                                 State_Met, FLUX,      RC )
!
! !USES:
!
    USE DEPO_MERCURY_MOD,   ONLY : DD_Hg2, WD_Hg2, DD_HgP, WD_HgP
#ifdef BPCH_DIAG
    USE DIAG03_MOD,         ONLY : AD03, ND03, AD03_RIV
#endif
    USE ErrCode_Mod
    USE ERROR_MOD,          ONLY : ERROR_STOP
    USE HCO_INTERFACE_MOD,  ONLY : HcoState
    USE HCO_EmisList_Mod,   ONLY : HCO_GetPtr
    USE Input_Opt_Mod,      ONLY : OptInput
    USE State_Chm_Mod,      ONLY : ChmState
    USE State_Diag_Mod,     ONLY : DgnState
    USE State_Grid_Mod,     ONLY : GrdState
    USE State_Met_Mod,      ONLY : MetState
    USE TIME_MOD
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
    TYPE(GrdState), INTENT(IN)    :: State_Grid  ! Grid State object
    TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology State object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
    TYPE(DgnState), INTENT(INOUT) :: State_Diag  ! Diagnostics State object
!
! !OUTPUT PARAMETERS:
!
    ! Flux of Hg(0) from the ocean [kg/s]
    REAL(fpp),      INTENT(OUT)   :: FLUX(State_Grid%NX,State_Grid%NY,N_Hg_CATS)
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!
! !REMARKS:
!  NOTE: The emitted flux may be negative when ocean conc. is very low.
!
!  ALSO NOTE: The ocean flux was tuned with GEOS-4 4x5 met fields.  We also
!  now account for the smaller grid size if using GEOS-4 2x25 met fields.
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
!  Hg(tot)aq REDUCTION RATE CONSTANTS
!
!  Hg(tot)aq reduction is split into biological and radiative reduction
!   (1.1 added to NPP for abiotic particles)
!
!    k_red     = k_red_bio + k_red_rad
!    k_red_rad = k_radbase * RADz     = ( s-1 W-1 m2 ) * ( W m-2 )
!    k_red_bio = k_biobse * NPP * 1.1 = ( s-1 mgC-1 d ) * ( mgC m-2 d-1 )
!
!  Hg(0)aq OXIDATION RATE CONSTANTS
!
!    k_ox      = k_oxbase * RADz + k_dark
!
!  k_dark is a constant dark oxidation component
!
!  RADz is the integrated ligth attenuation based on Beer-Lamberts law
!  (Schwarzenbach et al. 1993)
!
!     RADz = (1/(x1-x2))(RAD/EC)(1-e**-EC * x2)
!
!  x1  = surface depth (=0) (m)
!  x2  = depth of mixed layer (m)
!  EC  = extinction coefficient (m-1)
!  RAD = incomming radiation from GEOS5
!
!  Extinction coefficient
!  EC = ECwater + ECdoc * Cdoc (NPP/NPPavg) + ECchla * CHL/1000
!
!  ECwater = 0.0145 m-1
!  ECdoc   = 0.654 m-1
!  Cdoc    = 1.5 mgL-1
!  ECchla  = 31 m-1
!  CHL     = amount dependent on inputfile (mg/m3) but we need
!            mg/L so divide CHL by 1000
!____________________________________________________________________________
!
!  TOTAL ORGANIC CARBON AND SUSPENDED PARTICULATE MATTER (TOTAL BIOMASS)
!
!  Hg(II) - Hg(P) partitioning coefficient
!
!    Fraction of Hg2 = Frac_Hg2 =  1 / ( 1 + kd_part * SPM )
!
!  Kd_part is based on Mason et al. 1998 and Mason & Fitzgerald 1993. (L/kg)
!  SPM is converted to kg/L by 10E-9
!
!  SPM is Suspended particulate matter (kg/L)
!
!    SPM = ( OC_tot * 10 / MLD ) * 1.1
!
!  Total biomas is a proxy for SPM (mg/m3) used in Hg(II)
!  partitioning. Calculated by multiplying the standing
!  stock of organic carbon (OC_tot) with 10 (exp Bundy 2004)
!  1.1 is to include abiotic particles
!
!  OC_tot is the standing stock of organic carbon (mgC/m²)
!
!    OC_tot = C_tot * 80
!
!  Standing stock is calculated based on C:Chl ratio of 80 (wetzel et al 2006)
!
!  C_tot is the integrated pigment content in euphotic layer (mg/m2)
!
!  The parameters for calculating integrated Chl is based on a
!  model by Uitz et al (2006).
!
!  CHL   = average Chl a conc. detected by Modis (mg/m3)
!  Zm    = mixed layer depth (m)
!  Ze    = euphotic depth (PAR 1% of surface value (m)
!
!  C_tot differs dependent on the water being stratified
!  or well-mixed.
!___________________________________________________________________________
!
!  GAS EXCHANGE
!
!  Net flux from the ocean is given by the equation:
!
!    F = Kw * ( CHg0_aq - CHg0_atm / H )    (Lis & Slanter 1974)
!
!  Kw is the mass transfer coefficient (cm/h)
!    There are different possibilities for calculating Kw. The default is:
!
!    Kw = 0.25 * u^2 / SQRT ( Sc / ScCO2 )  (Nightingale et al. 2000)
!
!  u^2 is the square of the wind speed (10m above ground) (m²/s²)
!
!  Sc is the Schmidt # for Hg [unitless]
!     (ref: Poissant et al 2000; Wilke and Chang 1995)
!     to correct for seawater D0 is decreased by 6% as suggested
!     by Wanninkhof (1992)
!
!    Sc = v/D = (0.017 * exp(-0.025T))/D = kinematic viscosity/diffusivity
!
!  Diffusivity is calculated by:
!    D = (7.4*10D-8 scrt(2.26 * Mw) * TK) / (vi * N**0.6)
!
!    vi = viscocity of water
!    N  = molal volumen of mercury = 14.18
!
!  Viscocity is taken from Loux (2001)
!
!  H is the diemensionless Henrys coefficient for elemental mercury
!
!    H = exp (-2404.3/T - 6.92) where T is sea temp in K (Andersson et al. 2008)
!___________________________________________________________________________
!
!  PARTICLE SINKING
!
!  (from Sunderland & Mason 2007)
!
!  JorgC_kg = 0.1 (NPP**1.77) (MLD**-0.74)
!____________________________________________________________________________
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    LOGICAL, SAVE        :: FIRST = .TRUE.
    LOGICAL              :: IS_OCEAN_BOX
    CHARACTER(LEN=255)   :: FILENAME
    INTEGER              :: nAdvect,   NA
    INTEGER              :: I,         J,        NN, C
    INTEGER              :: N,         N_tot_oc
    INTEGER              :: NEXTMONTH, THISMONTH
    INTEGER              :: THISYEAR

    REAL(fpp)            :: A_M2,     DTSRCE,   MLDCM
    REAL(fpp)            :: CHg0aq,   CHg0,     vi,       JorgC_kg
    REAL(fpp)            :: TC,       TK,       Kw
    REAL(fpp)            :: Sc,       ScCO2,    USQ,      MHg
    REAL(fpp)            :: Hg2_RED,  Hg2_GONE, Hg2_CONV
    REAL(fpp)            :: FRAC_L,   FRAC_O,   H,        TOTDEP
    REAL(fpp)            :: oldMLD,   XTAU,     TOTDEPall
    REAL(fpp)            :: FUP(State_Grid%NX,State_Grid%NY,N_Hg_CATS)
    REAL(fpp)            :: FDOWN(State_Grid%NX,State_Grid%NY,N_Hg_CATS)
    REAL(fpp)            :: X,        Y,        D
    REAL(fpp)            :: NPP_tot,  A_ocean,  NPP_avg,  RADz
    REAL(fpp)            :: EC
    REAL(fpp)            :: k_red,    k_red_rad,  k_red_bio
    REAL(fpp)            :: k_ox
    REAL(fpp)            :: SPM,      Frac_Hg2, OC_tot_kg
    REAL(fpp)            :: Hg2aq_tot
    REAL(fpp)            :: C_tot,    Ze,       OC_tot,   Hg0_OX
    REAL(fpp)            :: Kd_part, k_ox_dark
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
    REAL*8,  PARAMETER :: RIVER_CHg(12) = (/0d0,  0d0,  0d0,  0d0, &
                                           80d0, 80d0, 27d0, 27d0, &
                                           27d0, 27d0,  0d0,  0d0 /)
    REAL*8             :: A_ARCTIC_OCEAN

    ! Small numbers to avoid dividing by zero
    REAL(fpp), PARAMETER :: SMALLNUM   = 1e-32_fpp
    REAL(fpp), PARAMETER :: NPPMINNUM   = 5e-2_fpp
    REAL(fpp), PARAMETER :: CHLMINNUM   = 1e-1_fpp

    ! For values from Input_Opt
    LOGICAL              :: LSPLIT
    LOGICAL              :: LArcticRiv,  LKRedUV

    ! Pointers
    REAL(fpp), POINTER   :: Spc(:,:,:,:)
    REAL(fpp), POINTER   :: Hg0aq(:,:,:)
    REAL(fpp), POINTER   :: Hg2aq(:,:,:)
    REAL(fpp), POINTER   :: HgPaq(:,:,:)

    ! Pointers to fields in the HEMCO data structure
    REAL(f4), POINTER    :: TOMS   (:,:)   ! O3
    REAL(f4), POINTER    :: TOMS_PD(:,:)   ! present day O3
    REAL(f4), POINTER    :: TOMS_LT(:,:)   ! long-term O3

    ! Strings
    CHARACTER(LEN=255) :: LOC = 'OCEAN_MERCURY_FLUX (GeosCore/ocean_mercury_mod.F90)'

    IF (LOCEANCOEF) THEN
       k_ox_dark = 10d0**(-5.2d0)
    ELSE
       k_ox_dark = 1d-7
    ENDIF

    !=================================================================
    ! OCEAN_MERCURY_FLUX begins here!
    !=================================================================

    ! Initialize pointers
    TOMS    => NULL()
    TOMS_PD => NULL()
    TOMS_LT => NULL()

    ! Copy values from Input_Opt
    LSPLIT       = Input_Opt%LSPLIT
    LArcticRiv   = Input_Opt%LArcticRiv
    LKRedUV      = Input_Opt%LKRedUV

    ! Point to fields in State_Chm
    Spc      => State_Chm%Species
    Hg0aq    => State_Chm%OceanHg0
    Hg2aq    => State_Chm%OceanHg2
    HgPaq    => State_Chm%OceanHgP

    ! Loop limit for use below
    IF ( LSPLIT ) THEN
       N_tot_oc = 2
    ELSE
       N_tot_oc = 1
    ENDIF

    ! Molecular weight of Hg (applicable to all tagged tracers)
    MHg = State_Chm%SpcData(1)%Info%emMW_g * 1e-3_fpp

    ! Get current month
    THISMONTH = GET_MONTH()

    ! Get current year to check if leap year (jaf, 8/12/11)
    THISYEAR  = GET_YEAR()

    !-----------------------------------------------
    ! Check tagged & total sums (if necessary)
    !-----------------------------------------------
    IF ( USE_CHECKS .and. LSPLIT ) THEN
       CALL CHECK_ATMOS_MERCURY( State_Chm, State_Grid, &
                                 'start of OCEAN_MERCURY_FLUX' )
       CALL CHECK_OCEAN_MERCURY( State_Chm, State_Grid, &
                                 'start of OCEAN_MERCURY_FLUX' )
       CALL CHECK_OCEAN_FLUXES ( State_Grid, &
                                 'start of OCEAN_MERCURY_FLUX' )
    ENDIF

    !-----------------------------------------------------------------
    ! Read O3 data from HEMCO
    !-----------------------------------------------------------------

    ! Now also read TOMS O3 column data
    IF ( LKRedUV ) THEN

       ! TOMS O3 columns [dobsons]
       CALL HCO_GetPtr( HcoState, 'TOMS_O3_COL', TOMS, RC )
       IF ( RC /= GC_SUCCESS ) THEN
          CALL ERROR_STOP ( 'Cannot get pointer to TOMS_O3_COL', LOC )
       ENDIF

       ! TOMS O3 columns, present-day [dobsons]
       CALL HCO_GetPtr( HcoState, 'TOMS_O3_PD', TOMS_PD, RC )
       IF ( RC /= GC_SUCCESS ) THEN
          CALL ERROR_STOP ( 'Cannot get pointer to TOMSPD_O3_COL', LOC )
       ENDIF

       ! TOMS O3 columns, long-term [dobsons]
       CALL HCO_GetPtr( HcoState, 'TOMS_O3_LT', TOMS_LT, RC )
       IF ( RC /= GC_SUCCESS ) THEN
          CALL ERROR_STOP ( 'Cannot get pointer to TOMSLT_O3_COL', LOC )
       ENDIF

    ENDIF

    !-----------------------------------------------
    ! Read monthly NPP, RADSW, MLD, UPVEL, dMLD data
    !-----------------------------------------------

    IF ( ITS_A_NEW_MONTH() ) THEN

       ! Get monthly MLD, NPP, CHL etc.
       CALL OCEAN_MERCURY_READ( Input_Opt, State_Grid, THISMONTH, THISYEAR, RC )

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
    DTSRCE = GET_TS_EMIS()

    !----------------------------------------------------------------
    ! Calculate total mean NPP (mg/m2/day) for later
    !----------------------------------------------------------------

    ! Initialize values
    NPP_tot = 0e+0_fpp
    A_ocean = 0e+0_fpp
    A_ARCTIC_OCEAN = 0e+0_fpp

    DO J = 1, State_Grid%NY
    DO I = 1, State_Grid%NX

       ! Grid box surface area [m2]
       A_M2 = State_Grid%Area_M2(I,J)

       ! Grid-box latitude [degrees]
       Y = State_Grid%YMid(I,J)

       ! Use fractional land type information from met fields (jaf, 4/26/11)
       ! FROCEAN is a constant, so to get correct ocean fraction we
       ! need to subtract the sea ice fraction.
       ! We now compute ocean chemistry for entire ocean grid box,
       ! irrespective of ice cover (jaf, 11/28/11)
       NPP_tot = NPP_tot + NPP(I,J) * A_M2 * State_Met%FROCEAN(I,J)
       A_ocean = A_ocean + A_M2 * State_Met%FROCEAN(I,J)
       IF ( Y >= 70 ) THEN
          !IF ( (State_Met%FROCEAN(I,J) - State_Met%FRSEAICE(I,J)) > 0d0 )
          IF ( State_Met%FROCEAN(I,J) > 0d0 ) &
               A_ARCTIC_OCEAN = A_ARCTIC_OCEAN + A_M2 * State_Met%FROCEAN(I,J)
       ENDIF

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
       RIVER_HG = FLOWNOW * RIVER_CHg(THISMONTH) * 1e-9_fpp / A_ARCTIC_OCEAN
    ENDIF

    ! Loop over surface boxes
    !$OMP PARALLEL DO       &
    !$OMP DEFAULT( SHARED ) &
    !$OMP PRIVATE( I,   vi,   A_M2,    Hg2_RED                         ) &
    !$OMP PRIVATE( J,   NN,   k_ox,    OC_tot,  Hg2_CONV               ) &
    !$OMP PRIVATE( N,   TK,   CHg0,    k_red_bio                       ) &
    !$OMP PRIVATE( C,   TC,   RADz,    Hg0_OX,  k_red_rad              ) &
    !$OMP PRIVATE( D,   EC,   k_red,   OLDMLD,  TOTDEPall              ) &
    !$OMP PRIVATE( Y,   Ze,   ScCO2,   FRAC_O,  Frac_Hg2,  Hg2aq_tot   ) &
    !$OMP PRIVATE( H,   Kw,   MLDCM,   TOTDEP,  OC_tot_kg              ) &
    !$OMP PRIVATE( X,   SPM,  CHg0aq,  Hg2_GONE                        ) &
    !$OMP PRIVATE( Sc,  Usq,  C_tot,   JorgC_kg                        ) &
    !$OMP PRIVATE( IS_OCEAN_BOX                                        ) &
    !$OMP PRIVATE( FRAC_OPEN_OCEAN,    FRAC_OCEAN_OR_ICE,  Snow_Hg2aq  ) &
    !$OMP PRIVATE( FRAC_REDUCIBLE,     UVI_RATIO                       ) &
    !$OMP PRIVATE( Kd_part                                             ) &
    !$OMP SCHEDULE( DYNAMIC )
    DO J = 1, State_Grid%NY
    DO I = 1, State_Grid%NX

       ! Grid box surface area [m2]
       A_M2       = State_Grid%Area_M2(I,J)

       ! Grid-box latitude [degrees]
       Y          = State_Grid%YMid(I,J)

       ! Initialize values
       Kw         = 0e+0_fpp
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
       ! Use fractional land type information from met fields. Also make sure
       ! we do not use boxes that are mostly sea ice for consistency
       ! FROCEAN is a constant, so to get correct ocean fraction we
       ! need to subtract the sea ice fraction. Don't let the fraction
       ! be less than zero (jaf, 4/26/11)
       !---------------------------------------------------------------
       ! Updated to distinguish between the fraction of the box that is
       ! open ocean and that which is ocean at depth (i.e. surface is
       ! either ocean or ice).
       ! We perform ocean chemistry whether or not there is
       ! ice cover, but ice will reduce both solar radiation input to
       ! ocean and mercury deposition to ocean (jaf, 11/28/11)
       !---------------------------------------------------------------
       FRAC_OPEN_OCEAN   = MAX( State_Met%FROCEAN(I,J) - &
                                State_Met%FRSEAICE(I,J), 0e+0_fpp)
       FRAC_OCEAN_OR_ICE = State_Met%FROCEAN(I,J)
       FRAC_O            = FRAC_OCEAN_OR_ICE
       IS_OCEAN_BOX      = ( FRAC_OCEAN_OR_ICE > 0e+0_fpp )

       ! Change ocean mass due to mixed layer depth change
       ! Keep before next IF so that we adjust mass in ice-covered boxes
       CALL MLD_ADJUSTMENT( I, J, OLDMLD*1e-2_fpp, MLDcm*1e-2_fpp, &
                            Input_Opt, State_Chm, State_Grid, State_Met )

       !---------------------------------------------------------------
       ! Deliver snowpack Hg to ocean if snow has melted
       ! Call before IF statement so that snow Hg is zeroed in land
       ! boxes as well as ocean boxes (jaf, 11/29/11)
       CALL DELIVER_SNOW_HG( SNOW_Hg2aq, I, J, State_Met, State_Chm, &
                             State_Diag )

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
          !!GanLuo+NPP(I,J) = MAX ( NPP(I,J) , NPPMINNUM )
          !NPP(I,J) = MAX ( NPP(I,J)*1.d0 , NPPMINNUM )
          !
          !!GanLuo+CHL(I,J) = MAX ( CHL(I,J) , CHLMINNUM )
          !CHL(I,J) = MAX ( CHL(I,J)*1.d0 , CHLMINNUM )

          ! Light attenuation (RADz) is calculated
          EC     = (EC_w + ( EC_doc * C_doc * ( NPP(I,J) / NPP_avg ) ) &
                 + ( ECchla * CHL(I,J) / 1000 ) )

          RADz   = ( 1 / ( MLDcm * 1e-2_fpp )) &
                   * (State_Met%SWGDN(I,J) / EC ) &
                   * ( 1 - EXP( -EC * ( MLDcm * 1e-2_fpp) ) )

          ! Reduce RADz in fractional sea ice boxes to only include
          ! the fraction of open ocean (i.e. assume zero light
          ! penetration through sea ice (jaf, 11/28/11)
          IF ( FRAC_OCEAN_OR_ICE > FRAC_OPEN_OCEAN ) &
               RADz = RADz * FRAC_OPEN_OCEAN

          !--------------------------------------------------------
          ! Hg(tot)aq reduction rate constants
          !--------------------------------------------------------
          k_red_rad   = ( k_radbase * RADz )

          ! NPP is increased by 0.1
          k_red_bio   = ( ( k_biobase * NPP(I,J) ) * 1.1 )

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
          !
          ! SPM = Suspended particulate matter (kg/L)
          !
          ! For Arctic, log10(Kd)=5.0, see Fisher et al. 2012 (jaf, 3/23/12)
          ! For Arctic, now use log10(Kd)=4.5 to match Fisher et al. (2013)
          ! When not in the arctic and LOCEANCOEF is true, use Kd_part
          ! and k_ox_dark from Song et al. ACP 2015, otherwise use ocean
          ! coefficients from Soerensen et al. EST 2010 ! sjs,08/12/2015
          IF ( Y >= 70.0e+0_fpp ) THEN
             Kd_part = 10e+0_fpp**(4.5e+0_fpp)
          ELSE IF ( LOCEANCOEF ) THEN
             Kd_part = 10e+0_fpp**(4.2e+0_fpp)
          ELSE
             Kd_part = 10e+0_fpp**(5.5e+0_fpp)
          ENDIF

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
          vi    = ( 10**( ( 1301.0e+0_fpp / ( 998.333e+0_fpp &
                  + 8.1855e+0_fpp                            &
                  * ( TC - 20.0e+0_fpp )+ 0.00585e+0_fpp     &
                  * (TC - 20.0e+0_fpp )**2 ) )               &
                  - 3.30233e+0_fpp ) ) * 100.0e+0_fpp

          ! Schmidt # for Hg [unitless]
          ! Sc = v/D = kinematic viscosity/diffusivity
          ! (ref: Poissant et al 2000; Wilke and Chang 1995)
          ! to correct for seawater D0 is decreased by 6% as suggested
          ! by Wanninkhof (1992)
          D  = 7.4e-8_fpp * sqrt( 2.26 * 18.0 ) * TK / ( ( 14.8**0.6 ) *vi )

          Sc = ( 0.017e+0_fpp * EXP( -0.025e+0_fpp * TC ) ) / D

          ! Schmidt # of CO2 [unitless] for CO2 in seawater at 20 degrees C
          ! The value is set to a constant based on other ocean studies
          ! (Gardfeld et al. 2003, Rolfhus & Fitzgerald 2004, Mason et al. 2001)
          !
          ! Correction of the Schmidt # with temperature based on Poissant
          ! et al. (2000) (for freshwatersystems).
          ScCO2  = 644.7e+0_fpp + TC * ( -6.16e+0_fpp + TC * &
                   ( 0.11e+0_fpp ) )

          ! Square of surface (actually 10m) wind speed [m2/s2]
          Usq    = State_Met%U10M(I,J)**2 + State_Met%V10M(I,J)**2

          !------------------------------------------------------
          ! Parameterizations for calculating water side mass
          ! trasfer coefficient
          !------------------------------------------------------
          ! Mass transfer coefficient [cm/h], from Nightingale et al. 2000
          Kw     = ( 0.25e+0_fpp * Usq ) / SQRT( Sc / ScCO2 )

          !-----------------------------------------------------
          ! Additional parameterizations:

          ! Nightinale et al. 2000 for instantanous winds
          !Kw     = ( 0.33e+0_fpp*SQRT(usq)+0.22e+0_fpp*Usq) / SQRT( Sc / ScCO2 )

          ! Lis and Merlivat 1986
          ! Has less emphasis on windspeed as a driver for evasion
          ! Gives a less total evasion than the Nigthingale et al. 2000
          !IF (SQRT(Usq) <= 3.6e+0_fpp ) THEN
          !   Kw = ( 0.17e+0_fpp * SQRT(Usq) * ( Sc / ScCO2 )**0.67e+0_fpp )
          !ELSE IF (SQRT(Usq) > 3.6e+0_fpp .and. SQRT(Usq) <= 13e+0_fpp ) THEN
          !   Kw = ( ( 2.8e+0_fpp * SQRT(Usq))-9.6 ) * ( Sc / ScCO2 )**0.5e+0_fpp
          !ELSE
          !   Kw = ( ( 5.9e+0_fpp * SQRT(Usq))-49.3 ) * ( Sc / ScCO2)**0.5e+0_fpp
          !ENDIF

          ! Wanninkhof et al (1992)
          !Kw     = ( 0.31e+0_fpp * Usq ) / SQRT( Sc / ScCO2 )

          !===========================================================
          ! Particulate sinking
          !===========================================================
          ! HgP sinking is based on Sunderland & Mason 2007.
          ! JorgC originally in gC m-2 year-1, which is convereted
          ! to kgC grid-1 timestep-1
          ! NPP is converted from mgC/m2/d-1 to gC/m2/year-1
          ! JorgC = 0.1 ( ( NPP * 12 )**1.77 ) *  MLD**n * M2 * Frac_O
          !         * 10^-3 * DTSRCE / ( 365 * 24 * 60 * 60 )
          JorgC_kg  = ( ( 0.1e+0_fpp * ((( NPP(I,J) * 365) / 1000 )     &
                      **1.77) * (( MLDcm * 1e-2_fpp )**(-0.74e+0_fpp) ) &
                      * A_M2 * FRAC_O * 1e-3_fpp) / ( 365.0e+0_fpp      &
                      * 24.0e+0_fpp * 60.0e+0_fpp * 60.0e+0_fpp ) )     &
                      * DTSRCE

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
                !Y = State_Grid%YMid(I,J)

                ! Grid box longitude [degrees]
                X = State_Grid%XMid(I,J)

                ! Get Hg category #
                IF ( C == 1 ) NN = ID_Hg_tot

                !--------------------------------------------------------
                ! Atlantic
                !--------------------------------------------------------
                IF ( ( X >= -80.0 .and. X < 25.0 )  .and. &
                     ( Y >= -25.0 .and. Y < 55.0 ) ) THEN    !(anls,100114)

                   ! Assign region tag
                   IF ( C == 2 ) NN = ID_Hg_atl

                   ! Hg0 (kg)
                   Hg0aq(I,J,NN) = Hg0aq(I,J,NN) + UPVEL(I,J) &
                        * ( MHg * A_M2 * FRAC_O * DTSRCE * CDeepatl(1) )

                   ! Hg2
                   Hg2aq(I,J,NN) = Hg2aq(I,J,NN) + UPVEL(I,J) &
                        * ( MHg * A_M2 * FRAC_O * DTSRCE * CDeepatl(2) )

                   ! Hg particulate
                   !IF ( C == 1 ) THEN
                   !   HgC(I,J)   = HgC(I,J) + UPVEL(I,J) &
                   !      * ( MHg * A_M2 * FRAC_O * DTSRCE * CDeepatl(3) )
                   HgPaq(I,J,NN)   = HgPaq(I,J,NN) + UPVEL(I,J) &
                        * ( MHg * A_M2 * FRAC_O * DTSRCE * CDeepatl(3) )
                   !ENDIF

                !--------------------------------------------------------
                ! North Pacific (west)
                !--------------------------------------------------------
                ELSE IF ( ( X >= -180.0 .and. X < -80.0 )  .and. &
                          ( Y >=   30.0 .and. Y <  70.0 ) ) THEN

                   ! Assign region tag
                   IF ( C == 2 ) NN = ID_Hg_npa

                   ! Hg0 (kg)
                   Hg0aq(I,J,NN) = Hg0aq(I,J,NN) + UPVEL(I,J) &
                        * ( MHg * A_M2 * FRAC_O * DTSRCE * CDeepnpa(1) )

                   ! Hg2
                   Hg2aq(I,J,NN) = Hg2aq(I,J,NN) + UPVEL(I,J) &
                        * ( MHg * A_M2 * FRAC_O * DTSRCE * CDeepnpa(2) )

                   ! Hg particulate
                   !IF ( C == 1 ) THEN
                   HgPaq(I,J,NN)   = HgPaq(I,J,NN) + UPVEL(I,J) &
                        * ( MHg * A_M2 * FRAC_O * DTSRCE * CDeepnpa(3) )
                   !ENDIF

                ! North Pacific (east)
                ELSE IF ( ( X >= 25.0 .and. X < 180.0 )  .and. &
                          ( Y >= 30.0 .and. Y <  70.0 ) ) THEN

                   ! Assign region tag
                   IF ( C == 2 ) NN = ID_Hg_npa

                   ! Hg0 (kg)
                   Hg0aq(I,J,NN) = Hg0aq(I,J,NN) + UPVEL(I,J) &
                        * ( MHg * A_M2 * FRAC_O * DTSRCE * CDeepnpa(1) )

                   ! Hg2
                   Hg2aq(I,J,NN) = Hg2aq(I,J,NN) + UPVEL(I,J) &
                       * ( MHg * A_M2 * FRAC_O * DTSRCE * CDeepnpa(2) )

                   ! Hg particulate &
                   !IF ( C == 1 ) THEN
                   HgPaq(I,J,NN)   = HgPaq(I,J,NN) + UPVEL(I,J) &
                        * ( MHg * A_M2 * FRAC_O * DTSRCE * CDeepnpa(3) )
                   !ENDIF

                !--------------------------------------------------------
                ! North Atlantic
                !--------------------------------------------------------
                ELSE IF ( ( X >= -80.0 .and. X < 25.0 )  .and. &
                          ( Y >=  55.0 .and. Y < 70.0 ) ) THEN

                   ! Assign region tag
                   IF ( C == 2 ) NN = ID_Hg_nat

                   ! Hg0 (kg)
                   Hg0aq(I,J,NN) = Hg0aq(I,J,NN) + UPVEL(I,J) &
                        * ( MHg * A_M2 * FRAC_O * DTSRCE * CDeepnat(1) )

                   ! Hg2
                   Hg2aq(I,J,NN) = Hg2aq(I,J,NN) + UPVEL(I,J) &
                        * ( MHg * A_M2 * FRAC_O * DTSRCE * CDeepnat(2) )

                   ! Hg particulate
                   !IF ( C == 1 ) THEN
                   HgPaq(I,J,NN)   = HgPaq(I,J,NN) + UPVEL(I,J) &
                        * ( MHg * A_M2 * FRAC_O * DTSRCE * CDeepnat(3) )
                   !ENDIF

                !--------------------------------------------------------
                ! South Atlantic
                !--------------------------------------------------------
                ELSE IF ( ( X >= -80.0 .and. X <  25.0 )  .and. &
                          ( Y >= -65.0 .and. Y < -25.0 ) ) THEN   !(anls,100114)

                   ! Assign region tag
                   IF ( C == 2 ) NN = ID_Hg_sat

                   ! Hg0 (kg)
                   Hg0aq(I,J,NN) = Hg0aq(I,J,NN) + UPVEL(I,J) &
                        * ( MHg * A_M2 * FRAC_O * DTSRCE * CDeepsat(1) )

                   ! Hg2
                   Hg2aq(I,J,NN) = Hg2aq(I,J,NN) + UPVEL(I,J) &
                        * ( MHg * A_M2 * FRAC_O * DTSRCE * CDeepsat(2) )

                   ! Hg particulate
                   !IF ( C == 1 ) THEN
                   HgPaq(I,J,NN)   = HgPaq(I,J,NN) + UPVEL(I,J) &
                        * ( MHg * A_M2 * FRAC_O * DTSRCE * CDeepsat(3) )
                   !ENDIF

                !--------------------------------------------------------
                ! Antarctic
                !--------------------------------------------------------
                ELSE IF ( Y >=  -90.0 .and. Y <  -65.0 ) THEN

                   ! Assign region tag
                   IF ( C == 2 ) NN = ID_Hg_ant

                   ! Hg0 (kg)
                   Hg0aq(I,J,NN) = Hg0aq(I,J,NN) + UPVEL(I,J) &
                        * ( MHg * A_M2 * FRAC_O * DTSRCE * CDeepant(1) )

                   ! Hg2
                   Hg2aq(I,J,NN) = Hg2aq(I,J,NN) + UPVEL(I,J) &
                        * ( MHg * A_M2 * FRAC_O * DTSRCE * CDeepant(2) )

                   ! Hg particulate
                   !IF ( C == 1 ) THEN
                   HgPaq(I,J,NN)   = HgPaq(I,J,NN) + UPVEL(I,J) &
                        * ( MHg * A_M2 * FRAC_O * DTSRCE * CDeepant(3) )
                   !ENDIF

                !--------------------------------------------------------
                ! Arctic
                !--------------------------------------------------------
                ELSE IF ( Y >=  70.0 .and. Y <  90.0 ) THEN

                   ! Assign region tag
                   IF ( C == 2 ) NN = ID_Hg_arc

                   ! Hg0 (kg)
                   Hg0aq(I,J,NN) = Hg0aq(I,J,NN) + UPVEL(I,J) &
                        * ( MHg * A_M2 * FRAC_O * DTSRCE * CDeeparc(1) )

                   ! Hg2
                   Hg2aq(I,J,NN) = Hg2aq(I,J,NN) + UPVEL(I,J) &
                        * ( MHg * A_M2 * FRAC_O * DTSRCE * CDeeparc(2) )

                   ! Hg particulate
                   !IF ( C == 1 ) THEN
                   HgPaq(I,J,NN)   = HgPaq(I,J,NN) + UPVEL(I,J) &
                        * ( MHg * A_M2 * FRAC_O * DTSRCE * CDeeparc(3) )
                   !ENDIF

                ELSE

                   ! Assign region tag
                   IF ( C == 2 ) NN = ID_Hg_ocn

                   ! Hg0 (kg)
                   Hg0aq(I,J,NN) = Hg0aq(I,J,NN) + UPVEL(I,J) &
                        * ( MHg * A_M2 * FRAC_O * DTSRCE * CDeep(1) )

                   ! Hg2
                   Hg2aq(I,J,NN) = Hg2aq(I,J,NN) + UPVEL(I,J) &
                        * ( MHg * A_M2 * FRAC_O * DTSRCE * CDeep(2) )

                   ! Hg particulate
                   !IF ( C == 1 ) THEN
                   HgPaq(I,J,NN)   = HgPaq(I,J,NN) + UPVEL(I,J) &
                        * ( MHg * A_M2 * FRAC_O * DTSRCE * CDeep(3) )
                   !ENDIF

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
                Hg0aq(I,J,NN) = Hg0aq(I,J,NN) &
                     * ( 1e+0_fpp + UPVEL(I,J) * DTSRCE / &
                     ( MLDcm * 1e-2_fpp ) )

                ! Hg2
                Hg2aq(I,J,NN) = Hg2aq(I,J,NN) &
                     * ( 1e+0_fpp + UPVEL(I,J) * DTSRCE / &
                     ( MLDcm * 1e-2_fpp ) )

                ! Hg particulate
                HgPaq(I,J,NN)  = HgPaq(I,J,NN) &
                     * ( 1e+0_fpp + UPVEL(I,J) * DTSRCE / &
                     ( MLDcm * 1e-2_fpp ) )

             ENDDO

          ENDIF

          !===========================================================
          ! Calculate reduction, conversion, sinking, evasion
          !
          ! (1) Hg2 <-> HgP and HgP sinks
          ! (2) Hg2 <-> Hg0 and Hg0 evades
          !
          ! NOTE: N is the GEOS-CHEM tracer # (for Spc)
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
             !
             ! Add flux of Hg(II) from rivers in Arctic (jaf, 12/8/11)
             ! Convert to kg/box/timestep and add to Hg2aq
             !--------------------------------------------------------
             IF ( (LArcticRiv) .AND. (Y >= 70) ) THEN
                Hg2aq(I,J,NN) = Hg2aq(I,J,NN) + RIVER_HG * &
                                A_M2 * DTSRCE * FRAC_O

#ifdef BPCH_DIAG
                !-----------------------------------------------------
                ! %%%%% ND03 (bpch) DIAGNOSTIC %%%%%
                !
                ! Flux of Hg(II) from rivers to ocean in Arctic [kg/s]
                ! NOTE: RIVER_HG is kg/m2/s, convert to kg/s
                !-----------------------------------------------------
                IF ( ( ND03 > 0 ) .and. ( NN == ID_Hg_tot ) ) &
                     AD03_riv(I,J) = AD03_riv(I,J) + RIVER_HG * &
                                     A_M2 * DTSRCE * FRAC_O
#endif

                !-----------------------------------------------------
                ! %%%%% HISTORY DIAGNOSTIC %%%%%
                !
                ! Flux of Hg(II) from rivers to ocean in Arctic [kg/s]
                ! NOTE: RIVER_HG is kg/m2/s, convert to kg/s
                !-----------------------------------------------------
                IF ( State_Diag%Archive_EmisHg2rivers ) THEN
                   State_Diag%EmisHg2rivers(I,J) = ( RIVER_HG * A_M2 * FRAC_O )
                ENDIF

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
             TOTDEPall = (WD_Hg2(I,J,NN) + DD_Hg2(I,J,NN) + &
                          WD_HgP(I,J,NN) + DD_HgP(I,J,NN) )

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

                IF ( (TOMS(I,J) < 0) .AND. (TOMS_LT(I,J) > 0) ) &
                     TOMS(I,J) = TOMS_LT(I,J)

                IF ( (TOMS(I,J) > 0) .AND. (TOMS_PD(I,J) > 0) ) THEN
                   UVI_RATIO = ( TOMS(I,J) / TOMS_PD(I,J) )**(-1.23e+0_fpp)
                   FRAC_REDUCIBLE = 0.2e+0_fpp + 0.2e+0_fpp * ( UVI_RATIO )
                ENDIF
             ENDIF

             ! Now use new FRAC_REDUCIBLE for photo-reduction only.
             ! Retain 40% for biological reduction
             !Hg2_RED       = Hg2aq(I,J,NN) * 0.4d0 * k_red * DTSRCE
             Hg2_RED = Hg2aq(I,J,NN) * DTSRCE * ( 0.4e+0_fpp * &
                       k_red_bio + FRAC_REDUCIBLE * k_red_rad )

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

             ! HgP ocean mass after conversion
             HgPaq(I,J,NN)   = Hg2aq_tot * ( 1 - Frac_Hg2)

             !----------------------------------------------------
             ! Conversion between OC and Hg
             !----------------------------------------------------
             ! Hg/C ratio based on HgP(kg) and Stock of organic C(kg)
             ! HgPaq_sunk funtion of C sunk and HgP/C ratio
             HgPaq_SUNK(I,J,NN)  = JorgC_kg * ( HgPaq(I,J,NN ) / OC_tot_kg )

             ! HgP ocean mass after sinking [kg]
             HgPaq(I,J,NN)   = HgPaq(I,J,NN) - HgPaq_SUNK(I,J,NN)
             HgPaq(I,J,NN)   = MAX ( HgPaq(I,J,NN) , 0e+0_fpp )

#ifdef BPCH_DIAG
             !-----------------------------------------------------
             ! %%%%% ND03 (bpch) DIAGNOSTIC %%%%%
             !
             ! Store organic carbon sinking to ocean [kgC/time]
             !-----------------------------------------------------
             IF ( ND03 > 0 ) THEN
                AD03(I,J,12,1) = AD03(I,J,12,1) + JorgC_kg
             ENDIF
#endif

             !-----------------------------------------------------
             ! %%%%% HISTORY DIAGNOSTIC %%%%%
             !
             ! Store organic carbon sinking to ocean [kgC/s]
             !-----------------------------------------------------
             IF ( State_Diag%Archive_FluxOCtoDeepOcean ) THEN
                State_Diag%FluxOCtoDeepOcean(I,J) = JorgC_kg / DTSRCE
             ENDIF

             !--------------------------------------------------------
             ! Calculate new Hg(0) mass
             !--------------------------------------------------------

             ! Hg0 tracer number (for Spc)
             N             = Hg0_Id_List(NN)

             ! Add converted Hg(II) and subtract converted Hg(0) mass
             ! to the ocean mass of Hg(0) [kg]
             Hg0aq(I,J,NN) = Hg0aq(I,J,NN) + Hg2_GONE

             !--------------------------------------------------------
             ! Calculate oceanic and gas-phase concentration of Hg(0)
             !--------------------------------------------------------

             ! Concentration of Hg(0) in the ocean [ng/L]
             CHg0aq        = ( Hg0aq(I,J,NN) * 1e+11_fpp   ) / &
                             ( A_M2          * FRAC_O ) / MLDcm

             ! Gas phase Hg(0) concentration: convert [kg] -> [ng/L]
             CHg0          = Spc(I,J,1,N) * 1.0e+9_fpp / &
                             State_Met%AIRVOL(I,J,1)

             !--------------------------------------------------------
             ! Compute flux of Hg(0) from the ocean to the air
             !--------------------------------------------------------

             ! Compute ocean flux of Hg0 [cm/h*ng/L]
             FLUX(I,J,NN)  = Kw * ( CHg0aq - ( CHg0 / H ) )

             ! TURN OFF EVASION
             !FLUX(I,J,NN)= MIN(0.,FLUX(I,J,NN))

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
                FLUX(I,J,NN)  = 0e+0_fpp
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
             IF ( (-FLUX(I,J,NN) * DTSRCE ) > Spc(I,J,1,N) ) THEN
                FLUX(I,J,NN) = -Spc(I,J,1,N) / DTSRCE
             ENDIF

             ! Cap FDOWN with available Hg(0) atm mass
             !IF ((FDOWN(I,J,NN)*DTSRCE)>Spc(I,J,1,N)) THEN
             !   FDOWN(I,J,NN) = Spc(I,J,1,N) / DTSRCE
             !ENDIF

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

             !Hgaq_tot = HgC(I,J) + Hg0aq(I,J,NN) + Hg2aq(I,J,NN)
             Hgaq_tot(I,J,NN) = HgPaq(I,J,NN) + Hg0aq(I,J,NN) + Hg2aq(I,J,NN)

          ENDDO

#ifdef BPCH_DIAG
          !-----------------------------------------------------------
          ! %%%%% ND03 (bpch) DIAGNOSTICS
          !
          ! Oceanic Hg quantities
          !-----------------------------------------------------------
          IF ( ND03 > 0 ) THEN

             ! eds 9/9/10 added NN loop
             DO NN = 1, N_HG_CATS

                ! Aqueous Hg(0) mass [kg]
                AD03(I,J,2,NN)  = AD03(I,J,2,NN)  + Hg0aq(I,J,NN)

                ! Aqueous Hg(II) mass [kg]
                AD03(I,J,7,NN)  = AD03(I,J,7,NN)  + Hg2aq(I,J,NN)

                ! Hg2 sunk deep into the ocean [kg/time]
                AD03(I,J,8,NN)= AD03(I,J,8,NN) + HgPaq_SUNK(I,J,NN)

                ! HgTot aqua mass [kg]
                AD03(I,J,10,NN) =AD03(I,J,10,NN) + Hgaq_tot(I,J,NN)

                ! HgP ocean mass [kg]
                AD03(I,J,11,NN) = AD03(I,J,11,NN) + HgPaq(I,J,NN)

                !Prior to 09 Nov 2011, H Amos --------------------------
                ! flux up and down (eck)
                !AD03(I,J,16) = AD03(I,J,16) + FUP(I,J,ID_Hg_tot)*DTSRCE
                !AD03(I,J,17) = AD03(I,J,17) + FDOWN(I,J,ID_Hg_tot)*DTSRCE
                ! Loop over tagged tracers (eds)
                IF (FLUX(I,J,NN) > 0e+0_fpp) THEN
                   AD03(I,J,16,NN) = AD03(I,J,16,NN) + &
                                     FLUX(I,J,NN) * DTSRCE
                ELSE IF (FLUX(I,J,NN) < 0e+0_fpp) THEN
                   AD03(I,J,17,NN) = AD03(I,J,17,NN) + &
                                     ( abs(FLUX(I,J,NN) * DTSRCE ))
                ENDIF
                !-------------------------------------------------------

             ENDDO

             ! Total HgII/HgP dep to open ocean [kg]
             AD03(I,J,20,1) = AD03(I,J,20,1) + TOTDEP

          ENDIF
#endif

          !-----------------------------------------------------------
          ! %%%%% HISTORY (aka netCDF) DIAGNOSTICS %%%%%
          !
          ! Oceanic Hg quantities
          !
          ! NOTE: We have converted units from kg
          !-----------------------------------------------------------

          ! NOTE: for now just use the total category
          NN = 1

          ! Mass of oceanic Hg0 [kg]
          IF ( State_Diag%Archive_MassHg0inOcean ) THEN
             State_Diag%MassHg0inOcean(I,J) = Hg0aq(I,J,NN)
          ENDIF

          ! Mass of oceanic Hg2 [kg]
          IF ( State_Diag%Archive_MassHg2inOcean ) THEN
             State_Diag%MassHg2inOcean(I,J) = Hg2aq(I,J,NN)
          ENDIF

          ! Mass of oceanic HgP [kg]
          IF ( State_Diag%Archive_MassHgPinOcean ) THEN
             State_Diag%MassHgPinOcean(I,J) = HgPaq(I,J,NN)
          ENDIF

          ! Total oceanic mercury [kg]
          IF ( State_Diag%Archive_MassHgTotalInOcean ) THEN
             State_Diag%MassHgTotalInOcean(I,J) = Hgaq_tot(I,J,NN)
          ENDIF

          ! Flux of Hg2 sunk to deep ocean [kg/s]
          IF ( State_Diag%Archive_FluxHg2toDeepOcean ) THEN
             State_Diag%FluxHg2toDeepOcean(I,J) = HgPaq_SUNK(I,J,NN) / DTSRCE
          ENDIF

          ! Ocean-to-air and air-to-ocean fluxes
          IF ( FLUX(I,J,NN) > 0e+0_fpp ) THEN

             ! Volatization flux of Hg0 from ocean to air
             ! NOTE: Units are kg/s
             IF ( State_Diag%Archive_FluxHg0fromOceantoAir ) THEN
                State_Diag%FluxHg0fromOceanToAir(I,J) = FLUX(I,J,NN)
             ENDIF

          ELSE IF ( FLUX(I,J,NN) < 0e+0_fpp)  THEN

             ! Drydep flux of Hg0 from air to ocean
             IF ( State_Diag%Archive_FluxHg0fromAirToOcean ) THEN
                State_Diag%FluxHg0fromAirToOcean(I,J) = ABS( FLUX(I,J,NN) )
             ENDIF

          ENDIF

          ! Total Hg2/HgP deposited to ocean
          ! NOTE: Changed unit to kg/s
          IF ( State_Diag%Archive_FluxHg2HgPfromAirToOcean ) THEN
             State_Diag%FluxHg2HgPfromAirToOcean(I,J) = TOTDEP / DTSRCE
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
    Spc      => NULL()
    Hg0aq    => NULL()
    Hg2aq    => NULL()
    HgPaq    => NULL()
    TOMS     => NULL()
    TOMS_PD  => NULL()
    TOMS_LT  => NULL()

    !=================================================================
    ! Check tagged & total sums (if necessary)
    !=================================================================
    IF ( USE_CHECKS .and. LSPLIT ) THEN
       CALL CHECK_ATMOS_MERCURY( State_Chm, State_Grid, &
                                 'end of OCEAN_MERCURY_FLUX' )
       CALL CHECK_OCEAN_MERCURY( State_Chm, State_Grid, &
                                 'end of OCEAN_MERCURY_FLUX' )
       CALL CHECK_OCEAN_FLUXES ( State_Grid, &
                                 'end of OCEAN_MERCURY_FLUX' )
       CALL CHECK_FLUX_OUT( State_Grid, FLUX, &
                            'end of OCEAN_MERCURY_FLUX' )
    ENDIF

  END SUBROUTINE OCEAN_MERCURY_FLUX
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ocean_mercury_read
!
! !DESCRIPTION: Subroutine OCEAN\_MERCURY\_READ reads in the mixed layer depth,
!  net primary productivity, upwelling and radiation climatology for each month.
!  This is needed for the ocean flux computation.(sas, cdh, bmy, 1/20/05, 3/28/06)
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE OCEAN_MERCURY_READ( Input_Opt, State_Grid, &
                                 THISMONTH, THISYEAR,   RC  )
!
! !USES:
!
    USE ErrCode_Mod
    USE Error_Mod,          ONLY : Error_Stop
    USE HCO_INTERFACE_MOD,  ONLY : HcoState
    USE HCO_EmisList_Mod,   ONLY : HCO_GetPtr
    USE Input_Opt_Mod,      ONLY : OptInput
    USE State_Grid_Mod,     ONLY : GrdState
    USE TIME_MOD,	      ONLY : ITS_A_LEAPYEAR
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)  :: Input_Opt   ! Input Options object
    TYPE(GrdState), INTENT(IN)  :: State_Grid  ! Grid State object
    INTEGER,        INTENT(IN)  :: THISMONTH
    INTEGER,        INTENT(IN)  :: THISYEAR
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT) :: RC          ! Success or failure?
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    LOGICAL, SAVE          :: FIRST = .TRUE.
    REAL*4                 :: ARRAY(State_Grid%NX,State_Grid%NY,1)
    INTEGER                :: I, J, X, Y

    ! Add time variables to get correct number of seconds per month
    ! and prevent negative flow later (jaf, 12/8/11)
    INTEGER	   :: LASTMONTH, DAYS_IN_LAST_MONTH
    INTEGER	   :: NEXTMONTH, DAYS_IN_THIS_MONTH
    INTEGER	   :: M(12) = (/ 31, 28, 31, 30, 31, 30, &
                                 31, 31, 30, 31, 30, 31 /)
    REAL*8,  PARAMETER :: SEC_PER_DAY = 3.6e+3_fpp * 24e+0_fpp

    ! Moved from ocean_mercury_flux (jaf,3/20/13)
    REAL*8,  PARAMETER   :: NPPMINNUM   = 5e-2_fpp
    REAL*8,  PARAMETER   :: CHLMINNUM   = 1e-1_fpp

    ! Arctic MLD parameters
    REAL*8,  PARAMETER   :: ArcticMLD_summer = 15e+2_fpp
    REAL*8,  PARAMETER   :: ArcticMLD_other  = 20e+2_fpp

    ! Strings
    CHARACTER(LEN=255) :: LOC = 'OCEAN_MERCURY_READ (GeosCore/ocean_mercury_mod.F90)'

    !=================================================================
    ! OCEAN_MERCURY_READ begins here!
    !=================================================================

    ! Assume success
    RC = GC_SUCCESS

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
    CALL HCO_GetPtr( HcoState, 'OCEAN_MLD', MLD, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL ERROR_STOP ( 'Cannot get pointer to OCEAN_MLD', LOC )
    ENDIF

    ! Convert [m] to [cm]
    MLD = MLD * 100e+0_f4

    !-----------------------------------------------------------------
    ! Chl from Modis [kg/m3] --> [mg/m3]
    !-----------------------------------------------------------------
    CALL HCO_GetPtr( HcoState, 'OCEAN_CHLA', CHL, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL ERROR_STOP ( 'Cannot get pointer to OCEAN_CHLA', LOC )
    ENDIF

    ! Convert [kg/m3] to [mg/m3]
    CHL = CHL * 1e6_f4

    !-----------------------------------------------------------------
    ! Arctic Chl from Jin et al 2012 [mg/m3]
    !-----------------------------------------------------------------
    CALL HCO_GetPtr( HcoState, 'OCEAN_CHLA_A', CHL_A, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL ERROR_STOP ( 'Cannot get pointer to OCEAN_ACHLA', LOC )
    ENDIF

    !-----------------------------------------------------------------
    ! Net primary productivity [kg/m2] -> [mg/m2]
    !-----------------------------------------------------------------
    CALL HCO_GetPtr( HcoState, 'OCEAN_NPP', NPP, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL ERROR_STOP ( 'Cannot get pointer to OCEAN_NPP', LOC )
    ENDIF

    ! Convert [kg/m2] -> [mg/m2]
    NPP = NPP * 1e6_f4

    !-----------------------------------------------------------------
    ! Arctic Net primary productivity [mg/m2/day]
    !-----------------------------------------------------------------
    CALL HCO_GetPtr( HcoState, 'OCEAN_NPP_A', NPP_A, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL ERROR_STOP ( 'Cannot get pointer to OCEAN_ANPP', LOC )
    ENDIF

    !-----------------------------------------------------------------
    ! Ekman upwelling velocity [m/s]
    !-----------------------------------------------------------------
    CALL HCO_GetPtr( HcoState, 'OCEAN_EKMAN_V', UPVEL, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL ERROR_STOP ( 'Cannot get pointer to OCEAN_EKMAN_V', LOC )
    ENDIF

    !-----------------------------------------------------------------
    ! MLD tendency, first half of month [cm/s]
    !-----------------------------------------------------------------
    CALL HCO_GetPtr( HcoState, 'OCEAN_dMLD1', dMLD1, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL ERROR_STOP ( 'Cannot get pointer to dMLD1', LOC )
    ENDIF

    ! Convert [m/s] to [cm/s]
    dMLD1 = dMLD1 * 100e+0_f4

    !-----------------------------------------------------------------
    ! MLD tendency, second half of month [cm/s]
    !-----------------------------------------------------------------
    CALL HCO_GetPtr( HcoState, 'OCEAN_dMLD2', dMLD2, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL ERROR_STOP ( 'Cannot get pointer to dMLD2', LOC )
    ENDIF

    ! Convert [m/s] to [cm/s]
    dMLD2 = dMLD2 * 100e+0_f4

    !-----------------------------------------------------------------
    ! Overwrite fields with Arctic specific parameters
    !-----------------------------------------------------------------

    !$OMP PARALLEL DO       &
    !$OMP DEFAULT( SHARED ) &
    !$OMP PRIVATE( I, J, X, Y )
    DO J = 1, State_Grid%NY
    DO I = 1, State_Grid%NX

       ! Grid-box longitude [degrees]
       X = State_Grid%XMid(I,J)

       ! Grid-box latitude [degrees]
       Y = State_Grid%YMid(I,J)

       ! In central Arctic, overwrite MLD values where global MLD
       ! dataset is based on very limited information. Instead use
       ! values from climatology of Toole et al. (2010). (jaf, 11/28/11)
       IF ( (Y >= 80e+0_fpp) .OR. &
          ( (Y >= 70e+0_fpp) .AND. ( X >= 50e+0_fpp .OR. &
                                     X <= 0e+0_fpp ) ) ) THEN
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
             dMLD2(I,J) = (ArcticMLD_summer - ArcticMLD_other)/ &
                          ( SEC_PER_DAY * DAYS_IN_LAST_MONTH )
          CASE( 7 ) ! july
             dMLD1(I,J) = (ArcticMLD_summer - ArcticMLD_other)/ &
                          ( SEC_PER_DAY * DAYS_IN_LAST_MONTH )
             dMLD2(I,J) = 0e+0_fpp
          CASE( 9 ) ! september
             dMLD1(I,J) = 0e+0_fpp
             dMLD2(I,J) = (ArcticMLD_other - ArcticMLD_summer)/ &
                          ( SEC_PER_DAY * DAYS_IN_LAST_MONTH )
          CASE ( 10 ) ! october
             dMLD1(I,J) = (ArcticMLD_other - ArcticMLD_summer)/ &
                          ( SEC_PER_DAY * DAYS_IN_LAST_MONTH )
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
    dFLOW1 = ( RIVERFLOW(THISMONTH) - RIVERFLOW(LASTMONTH) ) / &
             ( SEC_PER_DAY * DAYS_IN_LAST_MONTH )
    dFLOW2 = ( RIVERFLOW(NEXTMONTH) - RIVERFLOW(THISMONTH) ) / &
             ( SEC_PER_DAY * DAYS_IN_THIS_MONTH )

    ! Set FLOWNOW first time (jaf, 12/8/11)
    IF ( FIRST ) THEN
       FLOWNOW = RIVERFLOW(LASTMONTH) + dFLOW1 * &
                 SEC_PER_DAY * ( DAYS_IN_LAST_MONTH - 15 )
       FIRST = .FALSE.
    ENDIF

  END SUBROUTINE OCEAN_MERCURY_READ
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: mld_adjustment
!
! !DESCRIPTION: Subroutine MLD\_ADJUSTMENT entrains new water when mixed layer
!  depth deepens and conserves concentration (leaves mass behind) when mixed
!  layer shoals. (sas, cdh, bmy, 4/18/05, 3/28/06)
!  The MLD depth is constrained so that the mean monthly concentration equals
!  the concentration in the middle of the given month. The MLD hereafter
!  changes linearily until it reaches the middle of the next months where the
!  process is repeted (anls, 4/30/09)
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE MLD_ADJUSTMENT( I, J, MLDold, MLDnew, Input_Opt, &
                             State_Chm, State_Grid, State_Met )
!
! !USES:
!
    USE Input_Opt_Mod,     ONLY : OptInput
    USE State_Chm_Mod,     ONLY : ChmState
    USE State_Grid_Mod,    ONLY : GrdState
    USE State_Met_Mod,     ONLY : MetState
!
! !INPUT PARAMETERS: 
!
    INTEGER,        INTENT(IN) :: I, J
    REAL(fpp),      INTENT(IN) :: MLDold      ! Old ocean mixed layer depth [m]
    REAL(fpp),      INTENT(IN) :: MLDnew      ! New ocean mixed layer depth [m]
    TYPE(OptInput), INTENT(IN) :: Input_Opt   ! Input Options object
    TYPE(ChmState), INTENT(IN) :: State_Chm   ! Chemistry State object
    TYPE(GrdState), INTENT(IN) :: State_Grid  ! Grid State object
    TYPE(MetState), INTENT(IN) :: State_Met   ! Meteorology State object
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER               :: C,    NN,     N_tot_oc
    INTEGER               :: K, L
    REAL(fpp)             :: A_M2, DELTAH, FRAC_O,  MHg
    REAL(fpp)             :: X, Y                   !(added anls 01/05/09)

    ! For fields from Input_Opt
    LOGICAL               :: LSPLIT

    ! Pointers
    REAL(fpp), POINTER   :: Hg0aq(:,:,:)
    REAL(fpp), POINTER   :: Hg2aq(:,:,:)
    REAL(fpp), POINTER   :: HgPaq(:,:,:)

    !=================================================================
    ! MLD_ADJUSTMENT begins here!
    !=================================================================

    ! Copy values from Input_Opt
    LSPLIT       = Input_Opt%LSPLIT

    ! Point to fields in State_Chm
    Hg0aq       => State_Chm%OceanHg0
    Hg2aq       => State_Chm%OceanHg2
    HgPaq       => State_Chm%OceanHgP

    ! Loop limit for use below
    IF ( LSPLIT ) THEN
       N_tot_oc = 2
    ELSE
       N_tot_oc = 1
    ENDIF

    ! Grid box surface area [m2]
    A_M2   = State_Grid%Area_M2(I,J)

    ! Fraction of box that is ocean
    ! Use fractional land type information from met fields (jaf, 4/26/11)
    ! FROCEAN is a constant, so to get correct ocean fraction we
    ! need to subtract the sea ice fraction.
    ! We now apply ocean processes, including MLD adjustment, to
    ! entire ocean grid box, irrespective of sea ice cover (jaf,
    ! 11/28/11)
    FRAC_O     = State_Met%FROCEAN(I,J)

    ! Molecular weight of Hg (valid for all tagged tracers)
    MHg    = State_Chm%SpcData(ID_Hg_tot)%Info%emMW_g * 1e-3_fpp

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

       ! Add Cdeepatl to North Atlantic and Cdeep to rest if the world
       ! (anls, 01/05/09)

       ! Grid-box latitude [degrees]
       Y = State_Grid%YMid(I,J)

       ! Grid box longitude [degrees]
       X = State_Grid%XMid(I,J)

       ! Loop over total Hg (and ocean Hg if necessary)
       DO C = 1, N_tot_oc

          ! Get Hg category number
          IF ( C == 1 ) NN = ID_Hg_tot

          !------------------------------------------------
          ! Atlantic
          !------------------------------------------------
          IF ( ( X >= -80.0 .and. X < 25.0 )  .and. &
               ( Y >= -25.0 .and. Y < 55.0 ) ) THEN !(anls,100114)

             ! Assign region tag
             IF ( C == 2 ) NN = ID_Hg_atl

             ! Hg0
             Hg0aq(I,J,NN) = Hg0aq(I,J,NN) + &
                             ( DELTAH * CDeepatl(1) * MHg * A_M2 * FRAC_O )

             ! Hg2
             Hg2aq(I,J,NN) = Hg2aq(I,J,NN) + &
                             ( DELTAH * CDeepatl(2) * MHg * A_M2 * FRAC_O )

             ! HgP
             !IF ( C == 1 ) THEN
             HgPaq(I,J,NN) = HgPaq(I,J,NN) + &
                             ( DELTAH * CDeepatl(3) * MHg * A_M2 * FRAC_O )
             !ENDIF

          !------------------------------------------------
          ! North Pacific (west)
          !------------------------------------------------
          ELSE IF ( ( X >= -180.0 .and. X < -80.0 )  .and. &
                    ( Y >=   30.0 .and. Y <  70.0 ) ) THEN

             ! Assign region tag
             IF ( C == 2 ) NN = ID_Hg_npa

             ! Hg0
             Hg0aq(I,J,NN) = Hg0aq(I,J,NN) + &
                             ( DELTAH * CDeepnpa(1) * MHg * A_M2 * FRAC_O )

             ! Hg2
             Hg2aq(I,J,NN) = Hg2aq(I,J,NN) + &
                             ( DELTAH * CDeepnpa(2) * MHg * A_M2 * FRAC_O )

             ! HgP
             !IF ( C == 1 ) THEN
             HgPaq(I,J,NN) = HgPaq(I,J,NN) + &
                             ( DELTAH * CDeepnpa(3) * MHg * A_M2 * FRAC_O )
             !ENDIF

          !------------------------------------------------
          ! North Pacific (east)
          !------------------------------------------------
          ELSE IF ( ( X >= 25.0 .and. X < 180.0 )  .and. &
                    ( Y >= 30.0 .and. Y <  70.0 ) ) THEN

             ! Assign region tag
             IF ( C == 2 ) NN = ID_Hg_npa

             ! Hg0
             Hg0aq(I,J,NN) = Hg0aq(I,J,NN) + &
                             ( DELTAH * CDeepnpa(1) * MHg * A_M2 * FRAC_O )

             ! Hg2
             Hg2aq(I,J,NN) = Hg2aq(I,J,NN) + &
                             ( DELTAH * CDeepnpa(2) * MHg * A_M2 * FRAC_O )

             ! HgP
             !IF ( C == 1 ) THEN
             HgPaq(I,J,NN) = HgPaq(I,J,NN) + &
                             ( DELTAH * CDeepnpa(3) * MHg * A_M2 * FRAC_O )
             !ENDIF

          !------------------------------------------------
          ! North Atlantic
          !------------------------------------------------
          ELSE IF ( ( X >= -80.0 .and. X < 25.0 )  .and. &
                    ( Y >=  55.0 .and. Y < 70.0 ) ) THEN

             ! Assign region tag
             IF ( C == 2 ) NN = ID_Hg_nat

             ! Hg0
             Hg0aq(I,J,NN) = Hg0aq(I,J,NN) + &
                             ( DELTAH * CDeepnat(1) * MHg * A_M2 * FRAC_O )

             ! Hg2
             Hg2aq(I,J,NN) = Hg2aq(I,J,NN) + &
                             ( DELTAH * CDeepnat(2) * MHg * A_M2 * FRAC_O )

             ! HgP
             !IF ( C == 1 ) THEN
             HgPaq(I,J,NN) = HgPaq(I,J,NN) + &
                             ( DELTAH * CDeepnat(3) * MHg * A_M2 * FRAC_O )
             !ENDIF

          !------------------------------------------------
          ! South Atlantic
          !------------------------------------------------
          ELSE IF ( ( X >= -80.0 .and. X <  25.0 )  .and. &
                    ( Y >= -65.0 .and. Y < -25.0 ) ) THEN    !(anls,100114)

             ! Assign region tag
             IF ( C == 2 ) NN = ID_Hg_sat

             ! Hg0
             Hg0aq(I,J,NN) = Hg0aq(I,J,NN) + &
                             ( DELTAH * CDeepsat(1) * MHg * A_M2 * FRAC_O )

             ! Hg2
             Hg2aq(I,J,NN) = Hg2aq(I,J,NN) + &
                             ( DELTAH * CDeepsat(2) * MHg * A_M2 * FRAC_O )

             ! HgP
             !IF ( C == 1 ) THEN
             HgPaq(I,J,NN) = HgPaq(I,J,NN) + &
                             ( DELTAH * CDeepsat(3) * MHg * A_M2 * FRAC_O )
             !ENDIF

          !------------------------------------------------
          ! Antarctic
          !------------------------------------------------
          ELSE IF ( Y >=  -90.0 .and. Y <  -65.0 ) THEN

             ! Assign region tag
             IF ( C == 2 ) NN = ID_Hg_ant

             ! Hg0
             Hg0aq(I,J,NN) = Hg0aq(I,J,NN) + &
                             ( DELTAH * CDeepant(1) * MHg * A_M2 * FRAC_O )

             ! Hg2
             Hg2aq(I,J,NN) = Hg2aq(I,J,NN) + &
                             ( DELTAH * CDeepant(2) * MHg * A_M2 * FRAC_O )

             ! HgP
             !IF ( C == 1 ) THEN
             HgPaq(I,J,NN) = HgPaq(I,J,NN) + &
                             ( DELTAH * CDeepant(3) * MHg * A_M2 * FRAC_O )
             !ENDIF

          !------------------------------------------------
          ! Arctic
          !------------------------------------------------
          ELSE IF ( Y >=  70.0 .and. Y <  90.0 ) THEN

             ! Assign region tag
             IF ( C == 2 ) NN = ID_Hg_arc

             ! Hg0
             Hg0aq(I,J,NN) = Hg0aq(I,J,NN) + &
                             ( DELTAH * CDeeparc(1) * MHg * A_M2 * FRAC_O )

             ! Hg2
             Hg2aq(I,J,NN) = Hg2aq(I,J,NN) + &
                             ( DELTAH * CDeeparc(2) * MHg * A_M2 * FRAC_O )

             ! HgP
             !IF ( C == 1 ) THEN
             HgPaq(I,J,NN) = HgPaq(I,J,NN) + &
                             ( DELTAH * CDeeparc(3) * MHg * A_M2 * FRAC_O )
             !ENDIF

          ELSE

             ! Assign region tag
             IF ( C == 2 ) NN = ID_Hg_ocn

             ! Hg0
             Hg0aq(I,J,NN) = Hg0aq(I,J,NN) + &
                             ( DELTAH * CDeep(1) * MHg * A_M2 * FRAC_O )

             ! Hg2
             Hg2aq(I,J,NN) = Hg2aq(I,J,NN) + &
                             ( DELTAH * CDeep(2) * MHg * A_M2 * FRAC_O )

             ! HgP
             !IF ( C == 1 ) THEN
             HgPaq(I,J,NN) = HgPaq(I,J,NN) + &
                             ( DELTAH * CDeep(3) * MHg * A_M2 * FRAC_O )
             !ENDIF

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

    ! Free pointers
    Hg0aq    => NULL()
    Hg2aq    => NULL()
    HgPaq    => NULL()

  END SUBROUTINE MLD_ADJUSTMENT
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: check_atmos_mercury
!
! !DESCRIPTION: Subroutine CHECK\_ATMOS\_MERCURY tests whether the total and
!  tagged tracers the GEOS-CHEM species array sum properly within each grid box.
!  (cdh, bmy, 3/28/06)
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE CHECK_ATMOS_MERCURY( State_Chm, State_Grid, LOC )
!
! !USES:
!
    USE ERROR_MOD,           ONLY : ERROR_STOP
    USE State_Chm_Mod,       ONLY : ChmState
    USE State_Grid_Mod,      ONLY : GrdState
!
! !INPUT PARAMETERS:
!
    CHARACTER(LEN=*), INTENT(IN)    :: LOC         ! Calling routine
    TYPE(GrdState),   INTENT(IN)    :: State_Grid  ! Grid State object
    TYPE(ChmState),   INTENT(INOUT) :: State_Chm   ! Chemistry State object
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    LOGICAL            :: FLAG
    INTEGER            :: I,       J,       L
    INTEGER            :: N,       NN
    REAL(fpp)          :: Hg0_tot, Hg0_tag, RELERR0, ABSERR0
    REAL(fpp)          :: Hg2_tot, Hg2_tag, RELERR2, ABSERR2
    REAL(fpp)          :: HgP_tot, HgP_tag, RELERRP, ABSERRP

    ! Pointers
    REAL(fpp), POINTER :: Spc(:,:,:,:)

    !=================================================================
    ! CHECK_ATMOS_MERCURY begins here!
    !=================================================================

    ! Set error flags
    FLAG = .FALSE.

    ! Point to chemical species array [kg]
    Spc => State_Chm%Species

    ! Loop over grid boxes
    !OMP PARALLEL DO       &
    !OMP DEFAULT( SHARED ) &
    !OMP PRIVATE( I,       J,       L,       N,      NN            ) &
    !OMP PRIVATE( Hg0_tot, RELERR0, ABSERR0                        ) &
    !OMP PRIVATE( Hg2_tot, RELERR2, ABSERR2                        ) &
    !OMP PRIVATE( HgP_tot, RELERRP, ABSERRP                        ) &
    !OMP REDUCTION( +:     Hg0_tag, Hg2_tag, HgP_tag               )
    DO L = 1, State_Grid%NZ
    DO J = 1, State_Grid%NY
    DO I = 1, State_Grid%NX

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
       N       = Hg0_Id_List(ID_Hg_tot)
       Hg0_tot = Spc(I,J,L,N)

       ! Sum of tagged Hg(0)
       DO NN = 2, N_Hg_CATS
          N       = Hg0_Id_List(NN)
          Hg0_tag = Hg0_tag + Spc(I,J,L,N)
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
       N       = Hg2_Id_List(ID_Hg_tot)
       Hg2_tot = Spc(I,J,L,N)

       ! Sum of tagged Hg(II)
       DO NN = 2, N_Hg_CATS
          N       = Hg2_Id_List(NN)
          Hg2_tag = Hg2_tag + Spc(I,J,L,N)
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
       N       = HgP_Id_List(ID_Hg_tot)
       HgP_tot = Spc(I,J,L,N)

       ! Sum of tagged Hg(P)
       DO NN = 2, N_Hg_CATS
          N = HgP_Id_List(NN)
          IF ( N > 0 ) HgP_tag = HgP_tag + Spc(I,J,L,N)
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
          !OMP CRITICAL
          FLAG = .TRUE.
          WRITE( 6, 100 ) I, J, L, Hg0_tot, Hg0_tag, RELERR0, ABSERR0
          !OMP END CRITICAL
       ENDIF

       !----------------------------
       ! Hg(0) error is too large
       !----------------------------
       IF ( RELERR2 > MAX_RELERR .and. ABSERR2 > MAX_ABSERR ) THEN
          !OMP CRITICAL
          FLAG = .TRUE.
          WRITE( 6, 110 ) I, J, L, Hg2_tot, Hg2_tag, RELERR2, ABSERR2
          !OMP END CRITICAL
       ENDIF

       !----------------------------
       ! HgP error is too large
       !----------------------------
       IF ( RELERRP > MAX_RELERR .and. ABSERRP > MAX_ABSERR ) THEN
          !OMP CRITICAL
          FLAG = .TRUE.
          WRITE( 6, 120 ) I, J, L, HgP_tot, HgP_tag, RELERRP, ABSERRP
          !OMP END CRITICAL
       ENDIF
    ENDDO
    ENDDO
    ENDDO
    !OMP END PARALLEL DO

    ! Free pointer
    Spc => NULL()
    
    ! FORMAT strings
100 FORMAT( 'Hg0 error: ', 3i5, 4es13.6 )
110 FORMAT( 'Hg2 error: ', 3i5, 4es13.6 )
120 FORMAT( 'HgP error: ', 3i5, 4es13.6 )

    ! Stop if Hg0 and Hg2 errors are too large
    IF ( FLAG ) THEN
       CALL ERROR_STOP( 'Tagged Hg0, Hg2, HgP do not add up!', LOC )
    ENDIF

  END SUBROUTINE CHECK_ATMOS_MERCURY
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: check_ocean_mercury
!
! !DESCRIPTION: Subroutine CHECK\_OCEAN_MERCURY tests whether tagged tracers in
!  Hg0aq and Hg2aq add properly within each grid box. (cdh, bmy, 3/28/06)
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE CHECK_OCEAN_MERCURY( State_Chm, State_Grid, LOC )
!
! !USES:
!
    USE ERROR_MOD,           ONLY : ERROR_STOP
    USE State_Chm_Mod,       ONLY : ChmState
    USE State_Grid_Mod,      ONLY : GrdState
!
! !INPUT PARAMETERS:
!
    CHARACTER(LEN=*), INTENT(IN)    :: LOC         ! Calling routine
    TYPE(GrdState),   INTENT(IN)    :: State_Grid  ! Grid State object
    TYPE(ChmState),   INTENT(INOUT) :: State_Chm   ! Chemistry State object
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
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
    !OMP PARALLEL DO       &
    !OMP DEFAULT( SHARED ) &
    !OMP PRIVATE( I, J   ) &
    !OMP PRIVATE( Hg0_tot, Hg0_tag, RELERR0, ABSERR0 )
    !OMP PRIVATE( Hg2_tot, Hg2_tag, RELERR2, ABSERR2 )
    DO J = 1, State_Grid%NY
    DO I = 1, State_Grid%NX

       !--------------------------------------
       ! Relative and absolute errors for Hg0
       !--------------------------------------
       Hg0_tot = State_Chm%OceanHg0(I,J,ID_Hg_tot)
       Hg0_tag = SUM( State_Chm%OceanHg0(I,J,2:N_Hg_CATS) )
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
       Hg2_tot = State_Chm%OceanHg2(I,J,ID_Hg_tot)
       Hg2_tag = SUM( State_Chm%OceanHg2(I,J,2:N_Hg_CATS) )
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
          !OMP CRITICAL
          FLAG = .TRUE.
          WRITE( 6, 100 ) I, J, Hg0_tot, Hg0_tag, RELERR0, ABSERR0
          !OMP END CRITICAL
       ENDIF

       !--------------------------------------
       ! Hg(II) error is too large
       !--------------------------------------
       IF ( RELERR2 > MAX_RELERR .and. ABSERR2 > MAX_ABSERR ) THEN
          !OMP CRITICAL
          FLAG = .TRUE.
          WRITE( 6, 110 ) I, J, Hg2_tot, Hg2_tag, RELERR2, ABSERR2
          !OMP END CRITICAL
       ENDIF

    ENDDO
    ENDDO
    !OMP END PARALLEL DO

    ! FORMAT strings
100 FORMAT( 'Hg0aq error: ', 2i5, 4es13.6 )
110 FORMAT( 'Hg2aq error: ', 2i5, 4es13.6 )

    ! Stop if Hg0 and Hg2 errors are too large
    IF ( FLAG ) THEN
       CALL ERROR_STOP( 'Tagged Hg0aq, Hg2aq do not add up!', LOC )
    ENDIF

  END SUBROUTINE CHECK_OCEAN_MERCURY
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: check_ocean_fluxes
!
! !DESCRIPTION: Subroutine CHECK\_OCEAN\_FLUXES tests whether the drydep and
!  wetdep fluxes in DD_Hg2 and WD_Hg2 sum together in each grid box. (cdh, bmy,
!  3/28/06)
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE CHECK_OCEAN_FLUXES( State_Grid, LOC )
!
! !USES:
!
    USE DEPO_MERCURY_MOD,    ONLY : DD_Hg2, WD_Hg2, DD_HgP, WD_HgP
    USE ERROR_MOD,           ONLY : ERROR_STOP
    USE State_Grid_Mod,      ONLY : GrdState
!
! !INPUT PARAMETERS:
!
    CHARACTER(LEN=*), INTENT(IN) :: LOC         ! Calling routine
    TYPE(GrdState),   INTENT(IN) :: State_Grid  ! Grid State object
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
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
100 FORMAT( '     - In CHECK_OCEAN_FLUXES' )

    ! Set error condition flag
    FLAG = .FALSE.

    ! Loop over ocean surface boxes
    !OMP PARALLEL DO       &
    !OMP DEFAULT( SHARED ) &
    !OMP PRIVATE( I, J   ) &
    !OMP PRIVATE( DD_tot, DD_tag, DD_RELERR, DD_ABSERR ) &
    !OMP PRIVATE( WD_tot, WD_tag, WD_RELERR, WD_ABSERR )
    DO J = 1, State_Grid%NY
    DO I = 1, State_Grid%NX

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
          !OMP CRITICAL
          FLAG = .TRUE.
          WRITE( 6, 110 ) I, J, DD_tot, DD_tag, DD_RELERR, DD_ABSERR
          !OMP END CRITICAL
       ENDIF

       !---------------------------------------
       ! WD flux error is too large
       !---------------------------------------
       IF ( WD_RELERR > MAX_RELERR .and. WD_ABSERR > MAX_FLXERR ) THEN
          !OMP CRITICAL
          FLAG = .TRUE.
          WRITE( 6, 120 ) I, J, WD_tot, WD_tag, WD_RELERR, WD_ABSERR
          !OMP END CRITICAL
       ENDIF

    ENDDO
    ENDDO
    !OMP END PARALLEL DO

    ! FORMAT strings
110 FORMAT( 'DD_Hg2 flux error: ', 2i5, 4es13.6 )
120 FORMAT( 'WD_Hg2 flux error: ', 2i5, 4es13.6 )

    ! Stop if Hg0 and Hg2 errors are too large
    IF ( FLAG ) THEN
       CALL ERROR_STOP( 'Tagged DD, WD fluxes do not add up!', LOC )
    ENDIf

  END SUBROUTINE CHECK_OCEAN_FLUXES
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: check_flux_out
!
! !DESCRIPTION: Subroutine CHECK\_FLUX\_OUT tests whether tagged quantities in
!  FLUX sum together in each grid box. (cdh, bmy, 3/20/06)
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE CHECK_FLUX_OUT( State_Grid, FLUX, LOC )
!
! !USES:
!
    USE ERROR_MOD,           ONLY : ERROR_STOP
    USE State_Grid_Mod,      ONLY : GrdState
!
! !INPUT PARAMETERS:
!
    TYPE(GrdState),   INTENT(IN) :: State_Grid
    REAL(fpp),        INTENT(IN) :: FLUX(State_Grid%NX,State_Grid%NY,N_Hg_CATS)
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
100 FORMAT( '     - In CHECK_FLUX_OUT' )

    ! Set error condition flag
    FLAG = .FALSE.

    ! Loop over ocean surface boxes
    !OMP PARALLEL DO       &
    !OMP DEFAULT( SHARED ) &
    !OMP PRIVATE( I, J, FLX_tot, FLX_tag, FLX_err )
    DO J = 1, State_Grid%NY
    DO I = 1, State_Grid%NX

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
       IF ( FLX_RELERR > MAX_RELERR  .and. FLX_ABSERR > MAX_ABSERR ) THEN
          !OMP CRITICAL
          FLAG = .TRUE.
          WRITE( 6, 110 ) I, J, FLX_tot, FLX_tag, FLX_RELERR, FLX_ABSERR
          !OMP END CRITICAL
       ENDIF

    ENDDO
    ENDDO
    !OMP END PARALLEL DO

    ! FORMAT strings
110 FORMAT( 'FLX_Hg2 flux error: ', 2i5, 4es13.6 )

    ! Stop if Hg0 and Hg2 errors are too large
    IF ( FLAG ) THEN
       CALL ERROR_STOP( 'Tagged emission fluxes do not add up!', LOC )
    ENDIf

  END SUBROUTINE CHECK_FLUX_OUT
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_land_mercury
!
! !DESCRIPTION: Subroutine INIT\_OCEAN\_MERCURY allocates and zeroes all
!  module arrays. (sas, cdh, bmy, 1/19/05, 3/28/06)
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE INIT_OCEAN_MERCURY( Input_Opt, State_Chm, State_Grid, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE ERROR_MOD,          ONLY : ALLOC_ERR
    USE FILE_MOD,           ONLY : IOERROR,   IU_FILE
    USE Input_Opt_Mod,      ONLY : OptInput
    USE Species_Mod,        ONLY : Species
    USE State_Chm_Mod,      ONLY : ChmState
    USE State_Grid_Mod,     ONLY : GrdState
    USE TIME_MOD,           ONLY : GET_YEAR,  EXPAND_DATE
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
    TYPE(GrdState), INTENT(IN)    :: State_Grid  ! Grid State object
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
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! For Hg indexing (bmy, 4/26/16)
    INTEGER                  :: N
    TYPE(Species), POINTER   :: SpcInfo

    ! For riverflow and NPP interannual scaling (jaf, 3/26/12)
    INTEGER                  :: THISYEAR, THISYEARINDEX
    REAL(fpp), DIMENSION(31) :: RSCF_05, RSCF_06, RSCF_07
    REAL(fpp), DIMENSION(31) :: RSCF_08, RSCF_09, RSCF_10
    REAL(fpp), DIMENSION(31) :: NSCF

    !=================================================================
    ! INIT_OCEAN_MERCURY begins here!
    !=================================================================

    ! Assume success
    RC = GC_SUCCESS

    ! Exit if this is a dry-run
    IF ( Input_Opt%DryRun ) RETURN

    ! Store the # of tagHg categories in a module variable
    N_Hg_CATS =  State_Chm%N_Hg_CATS

    ! Hg species index corresponding to a given Hg category number
    Hg0_Id_List    => State_Chm%Hg0_Id_List
    Hg2_Id_List    => State_Chm%Hg2_Id_List
    HgP_Id_List    => State_Chm%HgP_Id_List

    ! Loop over all Hg species
    DO N = 1, State_Chm%nSpecies

       ! Point to the Species Database entry for tracer # N
       SpcInfo => State_Chm%SpcData(N)%Info

       ! Save ocean Hg categories into module variable indices
       SELECT CASE( TRIM( SpcInfo%Name ) )
       CASE( 'Hg0'     )
          ID_Hg_tot = SpcInfo%Hg_Cat
       CASE( 'Hg0_atl' )
          ID_Hg_atl = SpcInfo%Hg_Cat
       CASE( 'Hg0_nat' )
          ID_Hg_nat = SpcInfo%Hg_Cat
       CASE( 'Hg0_sat' )
          ID_Hg_sat = SpcInfo%Hg_Cat
       CASE( 'Hg0_npa' )
          ID_Hg_npa = SpcInfo%Hg_Cat
       CASE( 'Hg0_arc' )
          ID_Hg_arc = SpcInfo%Hg_Cat
       CASE( 'Hg0_ant' )
          ID_Hg_ant = SpcInfo%Hg_Cat
       CASE( 'Hg0_ocn' )
          ID_Hg_ocn = SpcInfo%Hg_Cat
       CASE DEFAULT
          ! Do nothing
       END SELECT

       ! Free pointer
       SpcInfo => NULL()

    ENDDO

    !-----------------------------------------------------------------
    ! Allocate these arrays regardless of whether you are using
    ! the dynamic ocean.  These are needed for Hg2 partitioning.
    ! (bmy, 3/30/15)
    !-----------------------------------------------------------------

    ALLOCATE( BULK_CONC( State_Grid%NX, State_Grid%NY, State_Grid%NZ), &
              STAT=RC )
    IF ( RC /= 0 ) CALL ALLOC_ERR( 'BULK_CONC' )
    BULK_CONC = 0e+0_fpp

    ALLOCATE( SO4_GC( State_Grid%NX, State_Grid%NY, State_Grid%NZ ), &
              STAT=RC )
    IF ( RC /= 0 ) CALL ALLOC_ERR( 'SO4_GC' )
    SO4_GC = 0e+0_fpp

    ALLOCATE( NIT_CONC( State_Grid%NX, State_Grid%NY, State_Grid%NZ ), &
              STAT=RC )
    IF ( RC /= 0 ) CALL ALLOC_ERR( 'NIT_CONC' )
    NIT_CONC = 0e+0_fpp

    ALLOCATE( NH4_CONC( State_Grid%NX, State_Grid%NY, State_Grid%NZ ), &
              STAT=RC )
    IF ( RC /= 0 ) CALL ALLOC_ERR( 'NH4_CONC' )
    NH4_CONC = 0e+0_fpp

    ALLOCATE( OC_CONC( State_Grid%NX, State_Grid%NY, State_Grid%NZ ), &
              STAT=RC )
    IF ( RC /= 0 ) CALL ALLOC_ERR( 'OC_CONC' )
    OC_CONC = 0e+0_fpp

    ALLOCATE( BC_CONC( State_Grid%NX, State_Grid%NY, State_Grid%NZ ), &
              STAT=RC )
    IF ( RC /= 0 ) CALL ALLOC_ERR( 'BC_CONC' )
    BC_CONC = 0e+0_fpp

    ALLOCATE( DST_CONC( State_Grid%NX, State_Grid%NY, State_Grid%NZ ), &
              STAT=RC )
    IF ( RC /= 0 ) CALL ALLOC_ERR( 'DST_CONC' )
    DST_CONC = 0e+0_fpp

    ALLOCATE( R( State_Grid%NX, State_Grid%NY, State_Grid%NZ ), &
              STAT=RC )
    IF ( RC /= 0 ) CALL ALLOC_ERR( 'R' )
    R = 0e+0_fpp

    ALLOCATE( SO4_CONC( State_Grid%NX, State_Grid%NY, State_Grid%NZ ), &
              STAT=RC )
    IF ( RC /= 0 ) CALL ALLOC_ERR( 'SO4_CONC' )
    SO4_CONC = 0e+0_fpp

    ALLOCATE( Fg( State_Grid%NX, State_Grid%NY, State_Grid%NZ ), &
              STAT=RC )
    IF ( RC /= 0 ) CALL ALLOC_ERR( 'Fg' )
    Fg = 0e+0_fpp

    ALLOCATE( Fp( State_Grid%NX, State_Grid%NY, State_Grid%NZ ), &
              STAT=RC )
    IF ( RC /= 0 ) CALL ALLOC_ERR( 'Fp' )
    Fp = 0e+0_fpp

    !eds 5/15/12 fix
    ALLOCATE( Hgaq_tot (State_Grid%NX, State_Grid%NY, N_Hg_CATS ), &
              STAT=RC )
    IF ( RC /= 0 ) CALL ALLOC_ERR( 'Hgaq_tot' )
    Hgaq_tot = 0e+0_fpp

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
       CDEEPSAT = (/ 0.8e-10_fpp, 4.1e-10_fpp, 4.1e-10_fpp/)  !eck, 10/19/11
       !reduce intermediate water mercury concentration in
       !South Atlantic Ocean to 0.9pM total
       !ref: low end of uncertainty range Sunderland and Mason '07
       CDEEPANT = (/ 1.0e-10_fpp, 3.2e-10_fpp, 3.2e-10_fpp /)
       ! Redistribute Hg fractions so Hg2 /= HgP (jaf, 11/28/11)
       CDEEPARC = (/ 1.2e-10_fpp, 1.2e-9_fpp,  2.8e-10_fpp /)
       CDEEPNPA = (/ 1.0e-10_fpp, 6.0e-10_fpp, 6.0e-10_fpp /)
    ENDIF

    ! Set monthly Arctic river flow rate in m3/s (jaf, 12/8/11)
    ! These come from the 8 largest Arctic rivers, with data
    ! from UNH/GRDC
    RIVERFLOW = (/ 2.2d4,  1.9d4,  1.8d4, 1.9d4, 9.6d4, 24.6d4, &
                   14.4d4, 10.2d4, 8.4d4, 6.2d4, 3.0d4, 2.3d4   /)

    ! Initialize default interannual NPP scaling factor (jaf, 3/19/13)
    NPP_SCF = 1e+0_fpp

    !-----------------------------------------------------------------
    ! Apply riverflow and set NPP scale factors if 2009 or prior
    !-----------------------------------------------------------------

    ! Use interannual flow scaling factors for May-Oct (jaf, 3/23/11)
    ! NOTE: scale factors are hard-coded to avoid reading ascii (ewl, 8/27/15)
    IF ( THISYEAR <= 2009) THEN

       ! May 1979-2009 river scale factors
       RSCF_05 = &
          (/ 0.63, 0.61, 1.12, 1.09, 0.50, 1.10, 0.68, 0.69, 0.59, 0.92, &
             1.11, 1.61, 1.02, 1.04, 0.59, 0.84, 1.13, 0.55, 2.03, 0.43, &
             1.46, 1.54, 1.13, 1.34, 1.00, 0.58, 1.31, 0.57, 1.58, 1.20, &
             1.00  /)

       ! Jun 1979-2009 river scale factors
       RSCF_06 = &
          (/ 1.18, 1.00, 0.89, 0.93, 1.15, 1.00, 0.94, 1.02, 0.90, 1.03, &
             1.17, 0.94, 1.00, 1.03, 1.00, 1.12, 0.87, 1.05, 0.79, 1.12, &
             0.95, 0.91, 1.14, 1.05, 0.91, 0.99, 0.85, 1.01, 1.01, 1.01, &
             1.06  /)

       ! Jul 1979-2009 river scale factors
       RSCF_07 = &
          (/ 0.97, 1.12, 1.03, 0.88, 1.17, 0.83, 0.92, 0.90, 1.08, 1.06, &
             1.11, 0.89, 0.94, 1.00, 1.00, 1.02, 0.97, 1.07, 1.03, 1.22, &
             1.03, 0.98, 0.88, 0.96, 0.93, 1.06, 0.89, 0.99, 1.13, 0.98, &
             0.97  /)

       ! Aug 1979-2009 river scale factors
       RSCF_08 = &
          (/ 1.20, 1.01, 1.00, 0.85, 1.03, 0.92, 1.05, 0.96, 0.99, 1.19, &
             0.95, 0.78, 0.83, 0.89, 0.96, 0.97, 0.85, 0.80, 1.09, 1.13, &
             1.08, 1.03, 1.13, 1.14, 0.88, 0.93, 1.02, 0.90, 1.37, 1.15, &
             0.91  /)

       ! Sep 1979-2008 river scale factors
       RSCF_09 = &
          (/ 1.42, 0.91, 1.02, 1.04, 1.06, 0.78, 0.95, 0.89, 0.67, 1.09, &
             1.03, 0.75, 0.84, 0.71, 0.89, 0.79, 0.92, 0.85, 1.20, 0.96, &
             0.92, 1.07, 0.97, 1.35, 0.83, 0.95, 1.13, 1.30, 1.45, 1.26, &
             1.00  /) ! Use 1.00 for 2009 since no scale factor available

       ! Oct 1979-2008 river scale factors
       RSCF_10 = &
          (/ 0.96, 0.89, 0.89, 0.80, 1.14, 0.79, 0.88, 1.04, 0.73, 1.11, &
             0.92, 0.90, 0.86, 0.94, 1.04, 1.04, 1.05, 1.06, 1.12, 0.90, &
             0.99, 1.09, 0.93, 1.14, 0.95, 1.17, 1.09, 1.29, 1.19, 1.10, &
             1.00  /) ! Use 1.00 for 2009 since no scale factor available

       ! 1979-2009 NPP scale factors (scaled to 2009)
       NSCF = &
          (/ 0.857, 0.807, 0.853, 0.835, 0.831, 0.857, 0.880, 0.823, &
             0.827, 0.827, 0.864, 0.933, 0.913, 0.823, 0.917, 0.857, &
             0.939, 0.798, 0.896, 0.914, 0.868, 0.927, 0.969, 0.955, &
             0.927, 0.894, 0.917, 1.035, 1.153, 1.073, 1.000 /)

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
    ALLOCATE( dMLD( State_Grid%NX, State_Grid%NY ), STAT=RC )
    IF ( RC /= 0 ) CALL ALLOC_ERR( 'dMLD' )
    dMLD = 0e+0_fpp

    ALLOCATE( HgPaq_SUNK (State_Grid%NX, State_Grid%NY, N_Hg_CATS ), &
              STAT=RC )
    IF ( RC /= 0 ) CALL ALLOC_ERR( 'HgPaq_SUNK' )
    HgPaq_SUNK = 0e+0_fpp

    ALLOCATE( MLDav( State_Grid%NX, State_Grid%NY ), STAT=RC )
    IF ( RC /= 0 ) CALL ALLOC_ERR( 'MLDav' )
    MLDav = 0e+0_fpp

    ALLOCATE( newMLD( State_Grid%NX, State_Grid%NY ), STAT=RC )
    IF ( RC /= 0 ) CALL ALLOC_ERR( 'newMLD' )
    newMLD = 0e+0_fpp

    !eds 10/19/10 fixing restart bug
    ALLOCATE( prevMLD( State_Grid%NX, State_Grid%NY ), STAT=RC )
    IF ( RC /=0 ) CALL ALLOC_ERR( 'prevMLD' )
    prevMLD = 0e+0_fpp

  END SUBROUTINE INIT_OCEAN_MERCURY
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cleanup_ocean_mercury
!
! !DESCRIPTION: Subroutine CLEANUP\_OCEAN\_MERCURY deallocates all module
!  arrays.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE CLEANUP_OCEAN_MERCURY
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC

    IF ( ALLOCATED( Hgaq_tot  ) ) DEALLOCATE( Hgaq_tot  )
    IF ( ALLOCATED( dMLD      ) ) DEALLOCATE( dMLD      )
    IF ( ALLOCATED( MLDav     ) ) DEALLOCATE( MLDav     )
    IF ( ALLOCATED( newMLD    ) ) DEALLOCATE( newMLD    )
    IF ( ALLOCATED( prevMLD   ) ) DEALLOCATE( prevMLD   )
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
    CHL         => NULL()
    CHL_A       => NULL()
    MLD         => NULL()
    NPP         => NULL()
    NPP_A       => NULL()
    UPVEL       => NULL()
    dMLD1       => NULL()
    dMLD2       => NULL()
    Hg0_Id_List => NULL()
    Hg2_Id_List => NULL()
    HgP_Id_List => NULL()

  END SUBROUTINE CLEANUP_OCEAN_MERCURY
!EOC
END MODULE OCEAN_MERCURY_MOD
