!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: mercury_mod.F90
!
! !DESCRIPTION: Module MERCURY\_MOD contains variables and routines for the
!  GEOS-CHEM mercury simulation.  Many choices of reaction mechanism and
!  model processes can be selected with logical switches located in
!  INIT\_MERCURY.
!\\
!\\
! !INTERFACE:
!
MODULE MERCURY_MOD
!
! !USES:
!
  USE DEPO_MERCURY_MOD,  ONLY : ADD_HG2_SNOWPACK
  USE DEPO_MERCURY_MOD,  ONLY : LHGSNOW
  USE OCEAN_MERCURY_MOD, ONLY : LDYNSEASALT
  USE OCEAN_MERCURY_MOD, ONLY : LPOLARBR
  USE OCEAN_MERCURY_MOD, ONLY : L_ADD_MBL_BR
  USE OCEAN_MERCURY_MOD, ONLY : LGEIA05
  USE OCEAN_MERCURY_MOD, ONLY : LVEGEMIS
  USE OCEAN_MERCURY_MOD, ONLY : LBRCHEM
  USE OCEAN_MERCURY_MOD, ONLY : LRED_JNO2
  USE OCEAN_MERCURY_MOD, ONLY : LRED_CLOUDONLY
  USE OCEAN_MERCURY_MOD, ONLY : LHALOGENCHEM
  USE OCEAN_MERCURY_MOD, ONLY : LHGAQCHEM
  USE OCEAN_MERCURY_MOD, ONLY : LHg2HalfAerosol
  USE OCEAN_MERCURY_MOD, ONLY : STRAT_BR_FACTOR
  USE OCEAN_MERCURY_MOD, ONLY : LAnthroHgOnly
  USE OCEAN_MERCURY_MOD, ONLY : LOHO3CHEM
  USE OCEAN_MERCURY_MOD, ONLY : LGCBROMINE
  USE OCEAN_MERCURY_MOD, ONLY : LnoUSAemis
  USE OCEAN_MERCURY_MOD, ONLY : LBROCHEM
  USE OCEAN_MERCURY_MOD, ONLY : LNEI2005
  USE OCEAN_MERCURY_MOD, ONLY : LInPlume
  USE OCEAN_MERCURY_MOD, ONLY : LOCEANCOEF
  USE PhysConstants           ! Physical constants
  USE PRECISION_MOD           ! For GEOS-Chem Precision (fp)

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC  :: CHEMMERCURY
  PUBLIC  :: CLEANUP_MERCURY
  PUBLIC  :: INIT_MERCURY
  PUBLIC  :: EMISSMERCURY
  PUBLIC  :: PARTITIONHG2
  PUBLIC  :: Reset_Hg_Diags
!
! !PRIVATE MEMBER FUNCTIONS:
!
  PRIVATE :: EMITHG
  PRIVATE :: SRCHg0
  PRIVATE :: SRCHg2
  PRIVATE :: SRCHgP
  PRIVATE :: MERCURY_READYR
  PRIVATE :: OHNO3TIME
  PRIVATE :: CALC_HG2_SEASALT_LOSSRATE
  PRIVATE :: RED_INPLUME_GRID
  PRIVATE :: DO_RED_INPLUME
  PRIVATE :: GET_HGBR_RATE
  PRIVATE :: GET_O3
  PRIVATE :: GET_OH
  PRIVATE :: GET_BR
  PRIVATE :: GET_JNO2
  PRIVATE :: GET_HGCL_RATE
  PRIVATE :: GET_HGBRY_RATE
  PRIVATE :: GET_HGAQ_RATE
  PRIVATE :: GET_NO2
  PRIVATE :: GET_HO2
  PRIVATE :: GET_Cl
  PRIVATE :: GET_ClO
  PRIVATE :: GET_HOCl
!
! !REMARKS:
!  Nomenclature:
!  ============================================================================
!  (1 ) Hg(0)  a.k.a. Hg0  : Elemental   mercury
!  (2 ) Hg(II) a.k.a. Hg2  : Divalent    mercury
!  (3 ) HgP                : Particulate mercury
!                                                                             .
!  Mercury Species (1-3 are always defined; 4-87 are defined for tagged runs)
!  ============================================================================
!  (1 ) Hg(0)              : Hg(0)  - total species
!  (2 ) Hg(II)             : Hg(II) - total species
!  (3 ) HgP                : HgP    - total species
!  ------------------------+---------------------------------------------------
!  (4 ) Hg0_can            : Hg(0) - Canadian Anthropogenic
!  (5 ) Hg0_usa            : Hg(0) - USA Anthropogenic
!  (6 ) Hg0_cam            : Hg(0) - Central American  Anthropogenic
!  (7 ) Hg0_sam            : Hg(0) - South American Anthropogenic
!  (8 ) Hg0_waf            : Hg(0) - West African Anthropogenic
!  (9 ) Hg0_eaf            : Hg(0) - East African Anthropogenic
!  (10) Hg0_saf            : Hg(0) - South African Anthropogenic
!  (11) Hg0_naf            : Hg(0) - North African Anthropogenic
!  (12) Hg0_eur            : Hg(0) - OECD European Anthropogenic
!  (13) Hg0_eeu            : Hg(0) - Eastern European Anthropogenic
!  (14) Hg0_sov            : Hg(0) - Former USSR Anthropogenic
!  (15) Hg0_mde            : Hg(0) - Middle Eastern Anthropogenic
!  (16) Hg0_sas            : Hg(0) - South Asian Anthropogenic
!  (17) Hg0_eas            : Hg(0) - East Asian Anthropogenic
!  (18) Hg0_sea            : Hg(0) - Southeast Asian Anthropogenic
!  (19) Hg0_jpn            : Hg(0) - Japanese Anthropogenic
!  (20) Hg0_ocn            : Hg(0) - Oceanian Anthropogenic
!  (21) Hg0_so             : Hg(0) - Organic Soils
!  (22) Hg0_bb             : Hg(0) - Biomass Burning
!  (23) Hg0_geo            : Hg(0) - Geogenic
!  (24) Hg0_atl            : Hg(0) - Middle Atlantic Subsurface Waters
!  (25) Hg0_nat            : Hg(0) - North Atlantic Subsurface Waters
!  (26) Hg0_sat            : Hg(0) - South Atlantic Subsurface Waters
!  (27) Hg0_npa            : Hg(0) - North Pacific Subsurface Waters
!  (28) Hg0_arc            : Hg(0) - Arctic Subsurface Waters
!  (29) Hg0_ant            : Hg(0) - Antarctic Subsurface Waters
!  (30) Hg0_oce            : Hg(0) - Indo-Pacific Subsurface Waters
!  (31) Hg0_str            : Hg(0) - Stratospheric Hg from Intial Conditions
!  (32-59) Same as (4-31) but for Hg(II)
!  (60-87) Same as (4-31) but for Hg(P)
!                                                                             .
!  References:
!  ============================================================================
!  (1 ) Hall, B. (1995). "The gas phase oxidation of elemental mercury by
!        ozone.", Water, Air, and Soil Pollution 80: 301-315.
!  (2 ) Sommar, J., et al. (2001). "A kinetic study of the gas-phase
!        reaction between the hydroxyl radical and atomic mercury."
!        Atmospheric Environment 35: 3049-3054.
!  (3 ) Selin, N., et al. (2007). "Chemical cycling and deposition of
!       atmospheric mercury: Global constraints from observations."
!       J. Geophys. Res. 112.
!  (4 ) Selin, N., et al. (2008). "Global 3-D land-ocean-atmospehre model
!       for mercury: present-day versus preindustrial cycles and
!       anthropogenic enrichment factors for deposition." Global
!       Biogeochemical Cycles 22: GB2011.
!  (5 ) Allison, J.D. and T.L. Allison (2005) "Partition coefficients for
!       metals in surface water, soil and waste." Rep. EPA/600/R-05/074,
!       US EPA, Office of Research and Development, Washington, D.C.
!  (6 ) Mintz, Y and G.K. Walker (1993). "Global fields of soil moisture
!       and land surface evapotranspiration derived from observed
!       precipitation and surface air temperature." J. Appl. Meteorol. 32 (8),
!       1305-1334.
!  (7 ) Soerensen, A. et al. (2010), An improved global model for air-sea
!       exchange of mercury: High concentrations over the North Atlantic,
!       Environ. Sci. Technol., 44, 8574-8580.
!  (8 ) Corbitt, E.S. et al. (2011), Global source-receptor relationships for
!       mercury deposition under present-day and 2050 emissions scenarios,
!       Environ. Sci. Technol., 45, 10477-10484.
!  (9 ) Street, D.G. et al. (2009), Projections of global mercury emissions
!       in 2050, Environ. Sci. Technol., 43, 2983-2988.
!  (10) Holmes, C.D., et al. (2010) Global atmospheric model for mercury
!       including oxidation by bromine atoms, AC&P, 10, 12,037-12,057.
!  (11) Parrella, J. et al. (2012), Tropospheric bromine chemistry:
!       implications for present and pre-industrial ozone and mercury, ACP.
!  (12) Pohler, D. et al. (2010), Observation of halogen species in the Amundsen
!       Gulf, Arctic, by active long-path differential optical absorption
!       spectroscopy, Proc. Natl. Acad. Sci, 107(15): 6528-6587.
!  (13) Prados-Roman, C. et al. (2011), Airborne DOAS limb measurements of
!       tropospheric trace gas profiles: case studies on the profile retrieval
!       of O4 and BrO, Atmos. Meas. Tech., 4: 1241-1260.
!
! !REVISION HISTORY:
!  06 Dec 2004 - N. (Eckley) Selin - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !PRIVATE TYPES:
!
  ! Parameters
  REAL(fp),  PARAMETER  :: SMALLNUM = 1e-20_fp

  ! Arrays
  INTEGER,  ALLOCATABLE :: AN_Hg0(:,:)    ! Index array for anth Hg0 regions
  INTEGER,  ALLOCATABLE :: AN_Hg2(:,:)    ! Index array for anth Hg2 regions
  INTEGER,  ALLOCATABLE :: AN_HgP(:,:)    ! Index array for anth HgP regions
  REAL(fp), ALLOCATABLE :: COSZM(:,:)     ! Max daily solar zenith angle
  REAL(fp), ALLOCATABLE :: EHg0_an(:,:)   ! Anth Hg0 emis [kg/s] - Total
  REAL(fp), ALLOCATABLE :: EHg0_can(:,:)  ! Anth Hg0 emis [kg/s] - Canada
  REAL(fp), ALLOCATABLE :: EHg0_usa(:,:)  ! Anth Hg0 emis [kg/s] - USA
  REAL(fp), ALLOCATABLE :: EHg0_cam(:,:)  ! Anth Hg0 emis [kg/s] - C America
  REAL(fp), ALLOCATABLE :: EHg0_sam(:,:)  ! Anth Hg0 emis [kg/s] - S America
  REAL(fp), ALLOCATABLE :: EHg0_waf(:,:)  ! Anth Hg0 emis [kg/s] - W Africa
  REAL(fp), ALLOCATABLE :: EHg0_eaf(:,:)  ! Anth Hg0 emis [kg/s] - E Africa
  REAL(fp), ALLOCATABLE :: EHg0_saf(:,:)  ! Anth Hg0 emis [kg/s] - S Africa
  REAL(fp), ALLOCATABLE :: EHg0_naf(:,:)  ! Anth Hg0 emis [kg/s] - N Africa
  REAL(fp), ALLOCATABLE :: EHg0_eur(:,:)  ! Anth Hg0 emis [kg/s] - Europe
  REAL(fp), ALLOCATABLE :: EHg0_eeu(:,:)  ! Anth Hg0 emis [kg/s] - E Europe
  REAL(fp), ALLOCATABLE :: EHg0_mde(:,:)  ! Anth Hg0 emis [kg/s] - Mid. East
  REAL(fp), ALLOCATABLE :: EHg0_sov(:,:)  ! Anth Hg0 emis [kg/s] - Fmr USSR
  REAL(fp), ALLOCATABLE :: EHg0_sas(:,:)  ! Anth Hg0 emis [kg/s] - S Asia
  REAL(fp), ALLOCATABLE :: EHg0_eas(:,:)  ! Anth Hg0 emis [kg/s] - E Asia
  REAL(fp), ALLOCATABLE :: EHg0_sea(:,:)  ! Anth Hg0 emis [kg/s] - SE Asia
  REAL(fp), ALLOCATABLE :: EHg0_jpn(:,:)  ! Anth Hg0 emis [kg/s] - Japan
  REAL(fp), ALLOCATABLE :: EHg0_oce(:,:)  ! Anth Hg0 emis [kg/s] - Oceania
  REAL(fp), ALLOCATABLE :: EHg2_an(:,:)   ! Anth Hg2 emis [kg/s] - Total
  REAL(fp), ALLOCATABLE :: EHg2_can(:,:)  ! Anth Hg2 emis [kg/s] - Canada
  REAL(fp), ALLOCATABLE :: EHg2_usa(:,:)  ! Anth Hg2 emis [kg/s] - USA
  REAL(fp), ALLOCATABLE :: EHg2_cam(:,:)  ! Anth Hg2 emis [kg/s] - C America
  REAL(fp), ALLOCATABLE :: EHg2_sam(:,:)  ! Anth Hg2 emis [kg/s] - S America
  REAL(fp), ALLOCATABLE :: EHg2_waf(:,:)  ! Anth Hg2 emis [kg/s] - W Africa
  REAL(fp), ALLOCATABLE :: EHg2_eaf(:,:)  ! Anth Hg2 emis [kg/s] - E Africa
  REAL(fp), ALLOCATABLE :: EHg2_saf(:,:)  ! Anth Hg2 emis [kg/s] - S Africa
  REAL(fp), ALLOCATABLE :: EHg2_naf(:,:)  ! Anth Hg2 emis [kg/s] - N Africa
  REAL(fp), ALLOCATABLE :: EHg2_eur(:,:)  ! Anth Hg2 emis [kg/s] - Europe
  REAL(fp), ALLOCATABLE :: EHg2_eeu(:,:)  ! Anth Hg2 emis [kg/s] - E Europe
  REAL(fp), ALLOCATABLE :: EHg2_mde(:,:)  ! Anth Hg2 emis [kg/s] - Mid. East
  REAL(fp), ALLOCATABLE :: EHg2_sov(:,:)  ! Anth Hg2 emis [kg/s] - Fmr USSR
  REAL(fp), ALLOCATABLE :: EHg2_sas(:,:)  ! Anth Hg2 emis [kg/s] - S Asia
  REAL(fp), ALLOCATABLE :: EHg2_eas(:,:)  ! Anth Hg2 emis [kg/s] - E Asia
  REAL(fp), ALLOCATABLE :: EHg2_sea(:,:)  ! Anth Hg2 emis [kg/s] - SE Asia
  REAL(fp), ALLOCATABLE :: EHg2_jpn(:,:)  ! Anth Hg2 emis [kg/s] - Japan
  REAL(fp), ALLOCATABLE :: EHg2_oce(:,:)  ! Anth Hg2 emis [kg/s] - Oceania
  REAL(fp), ALLOCATABLE :: EHgP_an(:,:)   ! Anth HgP emis [kg/s] - Total
  REAL(fp), ALLOCATABLE :: EHgP_can(:,:)  ! Anth HgP emis [kg/s] - Canada
  REAL(fp), ALLOCATABLE :: EHgP_usa(:,:)  ! Anth HgP emis [kg/s] - USA
  REAL(fp), ALLOCATABLE :: EHgP_cam(:,:)  ! Anth HgP emis [kg/s] - C America
  REAL(fp), ALLOCATABLE :: EHgP_sam(:,:)  ! Anth HgP emis [kg/s] - S America
  REAL(fp), ALLOCATABLE :: EHgP_waf(:,:)  ! Anth HgP emis [kg/s] - W Africa
  REAL(fp), ALLOCATABLE :: EHgP_eaf(:,:)  ! Anth HgP emis [kg/s] - E Africa
  REAL(fp), ALLOCATABLE :: EHgP_saf(:,:)  ! Anth HgP emis [kg/s] - S Africa
  REAL(fp), ALLOCATABLE :: EHgP_naf(:,:)  ! Anth HgP emis [kg/s] - N Africa
  REAL(fp), ALLOCATABLE :: EHgP_eur(:,:)  ! Anth HgP emis [kg/s] - Europe
  REAL(fp), ALLOCATABLE :: EHgP_eeu(:,:)  ! Anth HgP emis [kg/s] - E Europe
  REAL(fp), ALLOCATABLE :: EHgP_mde(:,:)  ! Anth HgP emis [kg/s] - Mid. East
  REAL(fp), ALLOCATABLE :: EHgP_sov(:,:)  ! Anth HgP emis [kg/s] - Fmr USSR
  REAL(fp), ALLOCATABLE :: EHgP_sas(:,:)  ! Anth HgP emis [kg/s] - S Asia
  REAL(fp), ALLOCATABLE :: EHgP_eas(:,:)  ! Anth HgP emis [kg/s] - E Asia
  REAL(fp), ALLOCATABLE :: EHgP_sea(:,:)  ! Anth HgP emis [kg/s] - SE Asia
  REAL(fp), ALLOCATABLE :: EHgP_jpn(:,:)  ! Anth HgP emis [kg/s] - Japan
  REAL(fp), ALLOCATABLE :: EHgP_oce(:,:)  ! Anth HgP emis [kg/s] - Oceania
  REAL(fp), ALLOCATABLE :: EHg0_am(:,:)   ! Artisinal mining Hg0 emis [kg/s]
  REAL(fp), ALLOCATABLE :: EHg0_oc(:,:,:) ! Ocean Hg0 emis [kg/s]
  REAL(fp), ALLOCATABLE :: EHg0_ln(:,:,:) ! Hg reemission from land [kg/s]
  REAL(fp), ALLOCATABLE :: EHg0_dist(:,:) ! Spatial dist of terrestrial Hg0
                                          !  sources [unitless]
  REAL(fp), ALLOCATABLE :: EHg0_geo(:,:)  ! Geogenic Hg0 emis [kg/s]
  REAL(fp), ALLOCATABLE :: EHg0_bb(:,:)   ! Biomass burning Hg0 emis [kg/s]
  REAL(fp), ALLOCATABLE :: EHg0_vg(:,:)   ! Vegetation Hg0 emis [kg/s
  REAL(fp), ALLOCATABLE :: EHg0_so(:,:)   ! Soil Hg0 emis [kg/s]
  REAL(fp), ALLOCATABLE :: EHg0_gtm(:,:)  ! GTMM Hg0 emis [kg/s]
  REAL(fp), ALLOCATABLE :: EHg0_snow(:,:,:) !Snow Hg0 emis [kg/s]
  REAL(fp), ALLOCATABLE :: TCOSZ(:,:)     ! Sum of solar zenith angle
  REAL(fp), ALLOCATABLE :: TTDAY(:,:)     ! Total daylight time at I,J [min]
  REAL(fp), ALLOCATABLE :: ZERO_DVEL(:,:) ! Zero drydep velocity [cm/s]
  REAL(fp), ALLOCATABLE :: HG2_SEASALT_LOSSRATE(:,:)

  ! For now, we need an emission array for the HG simulation
  ! that can be passed to vdiff_mod.F90 since Trac_Tend does
  ! not exist anymore (ckeller, 10/21/2014).
  REAL(fp), ALLOCATABLE, PUBLIC :: HG_EMIS(:,:,:)

  ! Pointers to fields in the HEMCO data structure.
  ! These need to be declared REAL(f4), aka REAL*4.
  ! (NOTE: We can set them to NULL here because
  ! they are globally SAVEd variables (bmy, 4/29/16)
  REAL(f4), POINTER :: O3(:,:,:)            => NULL()
  REAL(f4), POINTER :: OH_trop(:,:,:)       => NULL()
  REAL(f4), POINTER :: OH_strat(:,:,:)      => NULL()
  REAL(f4), POINTER :: JNO2(:,:,:)          => NULL()
  REAL(f4), POINTER :: HEM_NO2(:,:,:)       => NULL()
  REAL(f4), POINTER :: HEM_NO(:,:,:)        => NULL()
  REAL(f4), POINTER :: HEM_HOCl(:,:,:)      => NULL()
  REAL(f4), POINTER :: HEM_HO2_trop(:,:,:)  => NULL()
  REAL(f4), POINTER :: HEM_HO2_strat(:,:,:) => NULL()
  REAL(f4), POINTER :: HEM_CLO(:,:,:)       => NULL()
  REAL(f4), POINTER :: HEM_CL(:,:,:)        => NULL()
  REAL(f4), POINTER :: HEM_OA(:,:,:)        => NULL()
  REAL(f4), POINTER :: OCEAN_CONC(:,:,:)    => NULL()

  ! Hg species IDs
  INTEGER           :: N_Hg_CATS
  INTEGER           :: id_Hg0,     id_Hg2,     id_HgP
  INTEGER           :: ID_Hg_tot,  ID_Hg_can,  ID_Hg_usa
  INTEGER           :: ID_Hg_cam,  ID_Hg_sam,  ID_Hg_waf
  INTEGER           :: ID_Hg_eaf,  ID_Hg_saf,  ID_Hg_naf
  INTEGER           :: ID_Hg_eur,  ID_Hg_eeu,  ID_Hg_sov
  INTEGER           :: ID_Hg_mde,  ID_Hg_sas,  ID_Hg_eas
  INTEGER           :: ID_Hg_sea,  ID_Hg_jpn,  ID_Hg_oce
  INTEGER           :: ID_Hg_so,   ID_Hg_bb,   ID_Hg_geo
  INTEGER           :: ID_Hg_atl,  ID_Hg_nat,  ID_Hg_sat
  INTEGER           :: ID_Hg_npa,  ID_Hg_arc,  ID_Hg_ant
  INTEGER           :: ID_Hg_ocn,  ID_Hg_str

  ! Pointers for Hg indexing
  ! (NOTE: We can set them to NULL here because
  ! they are globally SAVEd variables (bmy, 4/29/16)
  INTEGER, POINTER  :: Hg0_Id_List(:) => NULL()
  INTEGER, POINTER  :: Hg2_Id_List(:) => NULL()
  INTEGER, POINTER  :: HgP_Id_List(:) => NULL()

CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: chemmercury
!
! !DESCRIPTION: Subroutine CHEMMERCURY is the driver routine for mercury
!  chemistry in the GEOS-CHEM module.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE CHEMMERCURY( Input_Opt,  State_Chm, State_Diag, &
                          State_Grid, State_Met, RC )
!
! !USES:
!
    USE DEPO_MERCURY_MOD,   ONLY : ADD_HG2_DD
    USE DEPO_MERCURY_MOD,   ONLY : ADD_HgP_DD
#ifdef BPCH_DIAG
    USE CMN_DIAG_MOD
    USE DIAG03_MOD,         ONLY : AD03_Hg2_Hg0, AD03_Hg2_O3
    USE DIAG03_MOD,         ONLY : AD03_Hg2_Br,  AD03_Hg2_OH
    USE DIAG03_MOD,         ONLY : AD03_Hg2_BRY,  AD03_Hg2_CLY
    USE DIAG03_MOD,         ONLY : AD03_Hg2_Br_Y
    USE DIAG03_MOD,         ONLY : AD03_Hg2_SS,  LD03
    USE DIAG03_MOD,         ONLY : AD03_Hg2_SSR, ND03
    USE DIAG03_MOD,         ONLY : AD03_Br
#endif
    USE ErrCode_Mod
    USE ERROR_MOD,          ONLY : ERROR_STOP
    USE ERROR_MOD,          ONLY : DEBUG_MSG
    USE ERROR_MOD,          ONLY : SAFE_DIV
    USE HCO_INTERFACE_MOD,  ONLY : HcoState
    USE HCO_EmisList_Mod,   ONLY : HCO_GetPtr
    USE GLOBAL_BR_MOD,      ONLY : GET_GLOBAL_BR
    USE Input_Opt_Mod,      ONLY : OptInput
    USE Species_Mod,        ONLY : Species
    USE State_Chm_Mod,      ONLY : ChmState
    USE State_Diag_Mod,     ONLY : DgnState
    USE State_Grid_Mod,     ONLY : GrdState
    USE State_Met_Mod,      ONLY : MetState
    USE TIME_MOD,           ONLY : GET_MONTH
    USE TIME_MOD,           ONLY : GET_TS_CHEM
    USE TIME_MOD,           ONLY : ITS_A_NEW_MONTH
    USE TIME_MOD,           ONLY : ITS_TIME_FOR_A3
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
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!
! !REMARKS:
!  Description of the chemistry mechanism:
!  ============================================================================
!  The comments and code below refer to 3 mercury species:
!     1. Hg0  or Hg(0)     : gaseous elemental mercury
!     2. Hg2g or Hg(II)g   : gaseous Hg(II) (comparable to RGM observations)
!     3. Hg2p or Hg(II)p   : particle-bound Hg(II) (comparable to PBM
!        observations)
!                                                                             .
!  Following Amos et al. (2012), the Hg(II) partitioning between
!  Hg2g and Hg2p is based on temperature and aerosol surface area.
!  There is no refractory particulate Hg.
!  In some legacy code below Hg2p is called HgP and Hg2g is called Hg2.
!                                                                             .
!  The differential equations for oxidation, reduction and deposition are:
!                                                                             .
!     d[Hg0]/dt = -(k_Ox + k_Dep0) [Hg0] + k_Red_Hg2g [Hg2g]
!                                        + k_Red_Hg2p [Hg2p]
!                                                                             .
!     d[Hg2g]/dt =  k_Ox [Hg0] - (k_Red_Hg2g + k_Dep_Hg2g) [Hg2g]
!                                                                             .
!     d[Hg2p]/dt =  -(k_Red_Hg2p + k_Dep_Hg2p) [Hg2p]
!                                                                             .
!  The chemical mechanism currently assumes that the product of Hg(0)
!  oxidation is in the gas phase. This can be easily changed.
!                                                                             .
!  Equilibrium partitioning between Hg2g and Hg2p is established at
!  the beginning and end of CHEMMERCURY.
!                                                                             .
!  Notes on chemistry and deposition:
!  ===========================================================================
!                                                                             .
!  (1 ) Oxidation: Hg(0) --> Hg(II):
!     Oxidation by Br, BrO, OH, and O3 can be selectively enabled or disabled
!     with switches (LBRCHEM, LBROCHEM, LOHO3CHEM) in INIT_MERCURY.
!                                                                             .
!     Hg(0)(g) + Br(g) --> + Br/OH --> Hg(II)g, rates are selected with
!        METHOD keyword below. Recommded kinetics are 'DonohoueYBBalabanov'
!        which use rates from Donohoue et al. (2006), Goodsite et al. (2004)
!        and Balabanov et al. (2005)
!                                                                             .
!     Hg(0)(g)+ O3(g) --> Hg(II) ,  k  = 3.0e-20 cm3 molec-1 s-1
!                                       Source: Hall, 1995
!
!     Hg(0)(g)+ OH(g) --> Hg(II) ,  k  = 8.7e-14 cm3 s-1
!                                       Source: Sommar et al. 2001
!                                                                             .
!  (2 ) Reduction:
!     Reduction rates can be scaled to NO2 photolysis rates or OH
!     concentrations with the switch LRED_JNO2.  In either case,
!     aqueous-phase photochemical reduction of Hg(II) is included based
!     on estimate of rate constant and scaled to NO2 photolysis or [OH].
!     The rate is tuned to match the global Hg(0) concentration and
!     seasonal cycle.
!                                                                             .
!  (3 ) Hg(0) dry deposition:
!     The dry deposition frequency is calculated by drydep_mod. If the
!     non-local PBL mixing scheme is used, however, all dry deposition is
!     done elsewhere.  The ocean module separately cacluates Hg(0) dry
!     deposition over ocean, so CHEMMERCURY only includes Hg(0) dry
!     deposition over land.
!                                                                             .
!  (4 ) Hg(II)g and Hg(II)p dry deposition:
!     The dry deposition frequencies are calculated by drydep_mod. If the
!     non-local PBL mixing scheme is used, however, all dry deposition is
!     done elsewhere.
!                                                                             .
!  (5 ) Sea-salt uptake of Hg(II)g
!     Hg(II)g is taken up into sea-salt aerosol and deposited to the ocean
!     surface at a rate based on wind speed and relative humidity. Hg(II)p
!     uptake into sea-salt aerosol may also occur at a slower rate, but is
!     not yet treated here.
!
! !REVISION HISTORY:
!  01 Oct 1995 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! SAVEd scalars
    LOGICAL, SAVE      :: FIRST = .TRUE.

    ! Scalars
    LOGICAL            :: LDYNOCEAN
    LOGICAL            :: LGTMM
    LOGICAL            :: LNLPBL
    LOGICAL            :: prtDebug
    INTEGER            :: NA, nAdvect
    INTEGER            :: I, J, L, K, N, NN, MONTH, NSTEPS, Hg_Cat
    REAL(fp)           :: K_DRYD0, K_DRYD2G, K_DRYD2P, K_SALT
    REAL(fp)           :: C_OH, C_O3, C_BR, C_BRO, C_HOCL
    REAL(fp)           :: C_NO2, C_OA, C_HO2, C_CLO, C_CL
    REAL(fp)           :: K_OH, K_O3, K_BR, K_BRO
    REAL(fp)           :: K_CL, K_R1CL, K_R2CL
    REAL(fp)           :: K_BR2, K_BROH, K_BRHO2, K_BRNO2, K_BRCLO
    REAL(fp)           :: K_BRBRO, K_BRCL, K_R1, K_R2, K_R2A, K_R3ANO2
    REAL(fp)           :: K_AQ_OX(3), K_HOCL, K_OH_AQ, K_O3_AQ
    REAL(fp)           :: K_OX, K_RED, K_DEP0, K_DEP2G, K_DEP2P
    REAL(fp)           :: DEP_Hg0, DEP_HG2G, DEP_HG2P, K_RED2P
    REAL(fp)           :: DEP_HG2G_DRY, DEP_HG2G_SALT, DEP_HG2P_DRY
    REAL(fp)           :: DEP_HG2_DRY
    REAL(fp)           :: GROSS_OX, GROSS_OX_BR
    REAL(fp)           :: GROSS_OX_OH, GROSS_OX_O3
    REAL(fp)           :: GROSS_OX_BRY, GROSS_OX_CLY
    REAL(fp)           :: GROSS_OX_BR2, GROSS_OX_BRBRO
    REAL(fp)           :: GROSS_OX_BRNO2, GROSS_OX_BRCLO
    REAL(fp)           :: GROSS_OX_BRHO2, GROSS_OX_BROH
    REAL(fp)           :: GROSS_RED, NET_OX
    REAL(fp)           :: F_PBL, FC, Faq, LWC
    REAL(fp)           :: AREA_CM2, DEP_DRY_FLX
    REAL(fp)           :: F_HG0_DEP
    REAL(fp)           :: DTCHEM, DT
    REAL(fp)           :: KMAX
    REAL(fp)           :: K0, R_DRY, RH3, Kstar, Hg2CR, TPL
    REAL(fp)           :: A(9,9), Xold(9), Xnew(9)
    REAL(fp)           :: s1(9), s2(9), s3(9), s4(9)
    REAL(fp)           :: K_BRY(11), K_CLY(3)

    ! Pointers
    REAL(fp), POINTER  :: Spc(:,:,:,:)
    REAL(fp), POINTER  :: T(:,:,:)

    ! Strings
    CHARACTER(LEN=255) :: ThisLoc
    CHARACTER(LEN=512) :: ErrMsg
!
! !DEFINED PARAMETERS:
!
    ! K for reaction Hg0 + O3  [cm3 /molecule /s] (Source: Hall 1995)
    REAL(fp), PARAMETER     :: K_HG_O3  = 3.0e-20_fp !3.0d-20 (Gas phase)

    ! K for reaction Hg0 + OH  [cm3 /molecule /s] (Source: Sommar 2001)
    REAL(fp), PARAMETER     :: K_HG_OH  = 8.7e-14_fp !8.7d-14 (Gas phase)

    ! K for reaction Hg0 + BrO  [cm3 /molecule /s]
    ! (Source: Raofie and Ariya 2003; 2004)
    REAL(fp), PARAMETER     :: K_HG_BRO  = 1e-14_fp !1d-15 - 1d-14

    ! K for reduction [cm3 /molecule /s] (Source: Selin 2007)
    ! This variable is unused if LREDJNO2 is TRUE
    REAL(fp), PARAMETER     :: K_RED_OH = 1e-10_fp!4.2d-10 from Noelle
                               !4d-10 works well for Hg+OH/O3

    ! K for reduction, using J_NO2, [unitless scale factor]
    ! Source: Holmes et al. 2010
    ! 3.5D-3 for Hg+Br simulation; 1.3D-2 for Hg+OH/O3 simulation
    ! These K_RED_JNO2 parameters have been defined based on GEOS-5,
    ! GEOS-FP, and MERRA2 runs with the chemistry from Horowitz 2017.
    ! It is set up to fail on compilation when an untested MET-GRID
    ! pair is used. (C. Thackray, 2/2018)
    ! This should be adapated and
    ! tested by other users as appropriate (J. Fisher, 3/2016)
    ! To simplify this code, K_RED_JNO2 should be moved to HEMCO_Config.rc
    ! or input.geos since it will vary with MET/GRID/chemistry (cpt)
    REAL(fp), PARAMETER     :: K_RED_JNO2 = 16.0e-2_fp

    ! Set of Hg/Br rate constants to use
    ! Recommended: DonohoueYBBalabanov, GoodsiteY, or DonohoueYB
    ! DonohoueYBBalabanov used by Holmes et al. 2010
    CHARACTER(LEN=*), PARAMETER  :: METHOD='GoodsiteUpdate'

    ! Objects
    TYPE(Species), POINTER :: SpcInfo

    !================================================================
    ! CHEMMERCURY begins here!
    !=================================================================

    ! Assume success
    RC        = GC_SUCCESS
    ErrMsg    = ''
    ThisLoc   = ' -> at ChemMercury (in GeosCore/mercury_mod.F90)'

    ! Copy values from Input_Opt
    LDYNOCEAN = Input_Opt%LDYNOCEAN
    LGTMM     = Input_Opt%LGTMM
    LNLPBL    = Input_Opt%LNLPBL
    prtDebug  = ( Input_Opt%LPRT .and. Input_Opt%amIRoot )

    ! Initialize pointers
    Spc      => State_Chm%Species   ! Chemical species array [kg]
    T        => State_Met%T         ! Temperature [K]
    SpcInfo  => NULL()

    ! Get info about this species from the database
    SpcInfo => State_Chm%SpcData(id_Hg2)%Info

    ! Get Henry's Law constant for Hg2 [mol/L/atm] (ewl, 1/5/16)
    K0    = SpcInfo%Henry_K0
    Hg2CR = SpcInfo%Henry_CR

    ! Calculate dry air gas constant in [L atm/K/mol] (ewl, 1/5/16)
    R_DRY = Rd * AIRMW / ATM

    ! Override the default value of LGCBROMINE and LRED_JNO2 depending
    ! on the settings in the HEMCO configuration file.  NOTE: This cannot
    ! be done in INIT_MERCURY because at that time the HEMCO configuration
    ! file has not been read yet. (bmy, 3/12/15)
    IF ( FIRST ) THEN
       CALL SET_OPTIONS_FROM_HEMCO( Input_Opt, State_Grid, RC )

       FIRST = .FALSE.
    ENDIF

    !=================================================================
    ! Read monthly mean Br, OH, O3, and J_NO2 fields
    !=================================================================
    IF ( ITS_A_NEW_MONTH() ) THEN

       ! Get the current month
       MONTH = GET_MONTH()

       ! Get a pointer to the monthly mean fields from HEMCO
       CALL HCO_GetPtr( HcoState, 'GLOBAL_OH_trop', OH_trop, RC )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Cannot get pointer to GLOBAL_OH_trop!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       CALL HCO_GetPtr( HcoState, 'GLOBAL_OH_strat', OH_strat, RC )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Cannot get pointer to GLOBAL_OH_strat!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       CALL HCO_GetPtr( HcoState, 'GLOBAL_O3', O3, RC )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Cannot get pointer to GLOBAL_O3!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       ! Br and BrO not directly from HEMCO
       CALL GET_GLOBAL_BR( Input_Opt, State_Grid, State_Met, MONTH, RC )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Cannot get pointer to GLOBAL_BR!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       IF ( prtDebug ) CALL DEBUG_MSG( '### CHEMMERC: a GET_GLOBAL_BR' )

       !-------------------------------------------------------------
       ! Some oxidants only required by certain mechanisms:
       !-------------------------------------------------------------
       IF ( LHALOGENCHEM ) THEN
          CALL HCO_GetPtr( HcoState,'GLOBAL_NO2', HEM_NO2, RC )
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = 'Cannot get pointer to GLOBAL_NO2!'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF

          CALL HCO_GetPtr( HcoState,'GLOBAL_NO', HEM_NO, RC )
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = 'Cannot get pointer to GLOBAL_NO!'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF

          CALL HCO_GetPtr( HcoState,'GLOBAL_HO2_trop', HEM_HO2_trop, RC )
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = 'Cannot get pointer to GLOBAL_HO2_trop!'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF

          CALL HCO_GetPtr( HcoState,'GLOBAL_HO2_strat', HEM_HO2_strat, RC )
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = 'Cannot get pointer to GLOBAL_HO2_strat!'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF

          CALL HCO_GetPtr( HcoState,'GLOBAL_ClO', HEM_CLO, RC )
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = 'Cannot get pointer to GLOBAL_ClO!'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF

          CALL HCO_GetPtr( HcoState,'GLOBAL_Cl', HEM_CL, RC )
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = 'Cannot get pointer to GLOBAL_NO2!'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF
       ENDIF
       !-------------------------------------------------------------

       IF ( LHGAQCHEM ) THEN
          CALL HCO_GetPtr( HcoState,'GLOBAL_HOCl', HEM_HOCL, RC )
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = 'Cannot get pointer to GLOBAL_HOCl!'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF
       ENDIF

       !-------------------------------------------------------------

       IF ( LRED_JNO2 ) THEN
          CALL HCO_GetPtr( HcoState, 'JNO2', JNO2,  RC )
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = 'Cannot get pointer to JNO2!'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF

          CALL HCO_GetPtr( HcoState,'GLOBAL_OA', HEM_OA, RC )
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = 'Cannot get pointer to GLOBAL_OA!'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF
       ENDIF

       !-------------------------------------------------------------
    ENDIF

    !=================================================================
    ! Perform chemistry on Hg species
    !=================================================================

    ! Chemistry timestep [s]
    DTCHEM = GET_TS_CHEM()

    ! Compute diurnal scaling for OH
    CALL OHNO3TIME( State_Grid )
    IF ( prtDebug ) CALL DEBUG_MSG( 'CHEMMERCURY: a OHNO3TIME' )

    ! Calculate the rate of sea salt aerosol uptake of Hg2
    IF ( LDYNSEASALT .AND. ITS_TIME_FOR_A3() ) THEN
       CALL CALC_HG2_SEASALT_LOSSRATE( State_Grid, State_Met )
       IF ( prtDebug ) CALL DEBUG_MSG( 'CHEMMERCURY: a SEASALT_LOSSRATE' )
    ENDIF

    ! Partition Hg(II) between gas and aerosol
    CALL PARTITIONHG2( Input_Opt, State_Chm, State_Grid, RC )

    IF ( prtDebug ) CALL DEBUG_MSG( 'CHEMMERCURY: before do loop' )

    !$OMP PARALLEL DO                                              &
    !$OMP DEFAULT( SHARED )                                        &
    !$OMP PRIVATE( I,       J,        L,        K,       N)        &
    !$OMP PRIVATE( NN,      NSTEPS )                               &
    !$OMP PRIVATE( K_DRYD0, K_DRYD2G, K_DRYD2P, K_SALT )           &
    !$OMP PRIVATE( C_OH,    C_O3,     C_BR,     C_BRO,   C_OA )    &
    !$OMP PRIVATE( C_HO2,   C_CLO,    C_CL,     C_NO2,   C_HOCL )  &
    !$OMP PRIVATE( K_OH,    K_O3,     K_BR,     K_BRO,   K_AQ_OX ) &
    !$OMP PRIVATE( K_CLY,   K_BRY,    K_CL,     K_BROH,  K_BRHO2)  &
    !$OMP PRIVATE( K_BRNO2)                                        &
    !$OMP PRIVATE( K_OH_AQ, K_O3_AQ,  K_R1CL,   K_R2CL , K_HOCL)   &
    !$OMP PRIVATE( K_BRCLO, K_BRBRO,  K_BRCL,   K_BR2)             &
    !$OMP PRIVATE( K_R3ANO2, K_R2A,   K_R2,     K_R1)              &
    !$OMP PRIVATE( K_OX,    K_RED,    K_RED2P,  K_DEP0,  K_DEP2G)  &
    !$OMP PRIVATE( K_DEP2P )                                       &
    !$OMP PRIVATE( F_PBL,   FC,       Faq,      LWC, RH3 )         &
    !$OMP PRIVATE( DT,      KMAX )                                 &
    !$OMP PRIVATE( A,       Xold,     Xnew,     s1, s2, s3, s4 )   &
    !$OMP PRIVATE( DEP_Hg0, DEP_Hg2g, DEP_Hg2p )                   &
    !$OMP PRIVATE( DEP_HG2G_DRY,      DEP_HG2G_SALT )              &
    !$OMP PRIVATE( DEP_HG2P_DRY,      DEP_HG2_DRY )                &
    !$OMP PRIVATE( GROSS_OX,     GROSS_OX_BR,    GROSS_OX_OH )     &
    !$OMP PRIVATE( GROSS_OX_O3,  GROSS_OX_BRY,   GROSS_OX_CLY )    &
    !$OMP PRIVATE( GROSS_OX_BR2, GROSS_OX_BRBRO, GROSS_OX_BRHO2 )  &
    !$OMP PRIVATE( GROSS_OX_BRNO2, GROSS_OX_BRCLO, GROSS_OX_BROH ) &
    !$OMP PRIVATE( GROSS_RED, NET_OX )                             &
    !$OMP PRIVATE( AREA_CM2,  DEP_DRY_FLX )                        &
    !$OMP PRIVATE( F_HG0_DEP, Hg_Cat,      Kstar,       TPL )
    !%%%% NOTE: LOOP IS IN WRONG ORDER, WILL BE MORE EFFICIENT IF REVERSED %%%%
    !%%%% Putting this loop in the correct order causes differences in the %%%%
    !%%%% AD03_HG2_SS diagnostic array when using mp vs sp (mps, 9/3/15)   %%%%
    DO I=1, State_Grid%NX
    DO J=1, State_Grid%NY
    DO L=1, State_Grid%NZ

       !------------------------------------------------------
       ! Dry deposition frequencies, 1/s
       !------------------------------------------------------

       ! Initialize
       ! These deposition frequencies never become non-zero in this
       K_DRYD0  = 0e+0_fp
       K_DRYD2G = 0e+0_fp
       K_DRYD2P = 0e+0_fp

       !------------------------------------------------------
       ! Sea-salt aerosol uptake
       !------------------------------------------------------

       ! IsWater returns true for ocean boxes.
       ! Fresh water may be TRUE in
       ! principle but is not at 4x5 resolution.
       ! Sea salt will not be active in coastal boxes that have some ocean
       ! but where IsWater is FALSE.

       ! RGM uptake on sea-salt aerosol, 1/s
       ! (based on SSA production rate (wind speed)
       ! and salinity (RH)) already zero over land
       K_SALT = HG2_SEASALT_LOSSRATE(I,J)

       ! Alternatively, uncomment the following lines
       ! to use a constant uptake rate (1/s) that is tuned to
       ! fit Okinawa RGM observations (Selin et al. 2007)
       !IF (State_Met%IsWater(I,J)) THEN
       !   K_SALT = 3.8D-5
       !ELSE
       !   K_SALT = 0e+0_fp
       !ENDIF

       !------------------------------------------------------
       ! Disable Hg(0) deposition to ocean, ice, snow
       !
       ! Hg(0) exchange with the ocean is handled by ocean_mercury_mod
       ! so disable deposition over water here.
       ! Turn off Hg(0) deposition to snow and ice because we haven't yet
       ! included emission from these surfaces and most field studies
       ! suggest Hg(0) emissions exceed deposition during sunlit hours.
       ! (Holmes et al. 2010)
       !------------------------------------------------------

       ! Area fraction of box with Hg0 dep. Default is dep everywhere
       F_Hg0_Dep = 1e+0_fp

       ! Deposition occurs over areas that are not ocean, snow or ice
       ! (jaf, 4/26/11)
       F_Hg0_Dep = 1e+0_fp - MAX( 0e+0_fp, MIN( 1e+0_fp, &
                   State_Met%FROCEAN(I,J) + &
                   State_Met%FRSNO(I,J)   + &
                   State_Met%FRLANDIC(I,J) ) )

       !------------------------------------------------------
       ! Total deposition frequencies, 1/s
       !------------------------------------------------------

       ! Fraction of box (I,J,L) underneath the PBL top [dimensionless]
       F_PBL = State_Met%F_UNDER_PBLTOP(I,J,L)

       ! Do dry deposition only for the fraction of the box in the PBL
       IF (F_PBL > 0.1e+0_fp) THEN
          K_DEP0  = F_PBL * ( K_DRYD0  * F_Hg0_Dep )
          K_DEP2G = F_PBL * ( K_DRYD2G + K_SALT )
          K_DEP2P = F_PBL *   K_DRYD2P
       ELSE
          K_DEP0  = 0e+0_fp
          K_DEP2G = 0e+0_fp
          K_DEP2P = 0e+0_fp
       ENDIF

       !-------------------------------------------------
       ! Oxidation rates
       !-------------------------------------------------

       ! Concentrations of oxidants, molec/cm3
       C_O3        = GET_O3( I, J, L,        State_Met             )
       C_OH        = GET_OH( I, J, L,        State_Met             )
       C_BR        = GET_BR( I, J, L, C_BRO, State_Diag, State_Grid, State_Met )

       IF ( LHALOGENCHEM ) THEN
          C_NO2       = GET_NO2(  I, J, L, State_Grid, State_Met, &
                                  HEM_NO(I,J,L), HEM_NO2(I,J,L), C_O3 )
          C_HO2       = GET_HO2(  I, J, L, State_Met )
          C_CLO       = GET_CLO(  I, J, L, State_Grid, State_Met )
          C_CL        = GET_CL(   I, J, L, State_Grid, State_Met )
       ENDIF

       IF ( LHGAQCHEM ) THEN
          C_HOCL      = GET_HOCl( I, J, L, State_Grid, State_Met )
       ENDIF

       IF ( LRED_JNO2 ) THEN
          C_OA        = HEM_OA(   I, J, L )
          TPL = State_Met%TropLev(I, J)
          IF ( L >= TPL ) THEN
             C_OA = 0.0e0_fp
          ENDIF
       ENDIF

       ! Oxidation rate by O3 and OH, 1/s
       IF (LOHO3CHEM) THEN
          K_O3 = K_HG_O3 * C_O3
          K_OH = K_HG_OH * C_OH
       ELSE
          K_O3 = 0e+0_fp
          K_OH = 0e+0_fp
       ENDIF

       ! Oxidation rate by Br, 1/s (accounts for HgBr intermediate)
       IF (LBRCHEM) THEN
          K_BR = GET_HGBR_RATE( I, J, L, C_BR, State_Met, C_OH, METHOD )
       ELSE
          K_BR = 0e+0_fp
       ENDIF

       ! Oxidation rate by BrO
       IF (LBROCHEM) THEN
          K_BRO = K_HG_BRO * C_BRO
       ELSE
          K_BRO = 0e+0_fp
       ENDIF

       ! Oxidation rate by Cl-initiated oxidation
       IF (LHALOGENCHEM) THEN
          K_CLY = GET_HGCL_RATE( I, J, L, C_CL,  C_OH,  C_HO2, &
                                 C_NO2,   C_CLO, C_BRO, C_BR,  State_Met )
          K_R1CL = K_CLY(2)
          K_R2CL = K_CLY(3)
          K_CL = K_CLY(1)
       ELSE
          K_CL = 0d0
          K_R1CL=0d0
          K_R2CL=0d0
       ENDIF

       ! Heterogeneous oxidation by OH, O3, and HOCL
       IF (LHGAQCHEM) THEN
          K_AQ_OX = GET_HGAQ_RATE( I, J, L, C_BR, C_OH, C_O3, C_HOCL, &
                                   State_Met )
          K_HOCL  = K_AQ_OX(1)
          K_OH_AQ = K_AQ_OX(2)
          K_O3_AQ = K_AQ_OX(3)
       ELSE
          K_AQ_OX = 0e0_fp
          K_HOCL  = 0e0_fp
          K_OH_AQ = 0e0_fp
          K_O3_AQ = 0e0_fp
       ENDIF

       ! Oxidation through Hg0 + Br -> HgBr, HgBr + Y -> Hg2
       IF (LHALOGENCHEM) THEN
          K_BRY = GET_HGBRY_RATE( I, J, L, C_BR,  C_OH,  C_NO2, &
                                  C_HO2,   C_CLO, C_BRO, C_CL,  State_Met )
          K_BR2    = K_BRY(1)
          K_BROH   = K_BRY(2)
          K_BRHO2  = K_BRY(3)
          K_BRNO2  = K_BRY(4)
          K_BRCLO  = K_BRY(5)
          K_BRBRO  = K_BRY(6)
          K_BRCL   = K_BRY(7)
          ! The following rates for diagnostic purposes:
          K_R1     = K_BRY(8)
          K_R2     = K_BRY(9)
          K_R2A    = K_BRY(10)
          K_R3ANO2 = K_BRY(11)
       ELSE
          K_BR2    = 0e+0_fp
          K_BROH   = 0e+0_fp
          K_BRHO2  = 0e+0_fp
          K_BRNO2  = 0e+0_fp
          K_BRCLO  = 0e+0_fp
          K_BRBRO  = 0e+0_fp
          K_BRCL   = 0e+0_fp
          K_R1     = 0e+0_fp
          K_R2     = 0e+0_fp
          K_R2A    = 0e+0_fp
          K_R3ANO2 = 0e+0_fp
       ENDIF

       ! Total oxidation rate, 1/s
       K_OX = K_O3    + K_OH    + K_BR    + K_BRO   + K_BR2  + &
              K_BRHO2 + K_BRNO2 + K_BRCLO + K_BRCL  + K_HOCL + &
              K_OH_AQ + K_O3_AQ + K_CL    + K_BRBRO + K_BROH

       !---------------------------------------------
       ! Reduction rates (requires liquid water content)
       !---------------------------------------------

       ! Get cloud fraction from met fields [unitless]
       FC = State_Met%CLDF(I,J,L)

       ! Get grid-average liquid water content [m3 H2O/m3 air] from met
       ! Units: [kg H2O/kg air] * [kg air/m3 air] * [m3 H2O/1e3 kg H2O]
       LWC = State_Met%QL(I,J,L) * State_Met%AIRDEN(I,J,L) * 1e-3_fp

       ! LWC is a grid-box averaged quantity. To improve representation
       ! of sulfate chemistry, we divide LWC by the cloud fraction and
       ! compute sulfate chemistry based on the LWC within the cloud.  We
       ! get the appropriate grid-box averaged mass of SO2 and sulfate by
       ! multiplying these quantities by FC AFTER computing the aqueous
       ! sulfur chemistry within the cloud. (lzh, jaf, bmy, 5/27/11)
       ! Adopted for mercury chemistry (yzh, 11/9/2011)
       IF ( LRED_CLOUDONLY ) THEN
          LWC     = SAFE_DIV( LWC, FC, 0e+0_fp )
       ENDIF

       ! Define fraction of in-cloud Hg(II) which is in aqueous solution
       ! [dimensionless]

       ! Make partitioning to water temperature-dependent
       ! (cpt, 1/31/2018) Should just use COMPUTE_L2G
       Kstar = K0 * EXP( Hg2CR * ( ( 1e0_fp / T(I,J,L) ) &
               - ( 1e0_fp / 298e0_fp ) ) )
       Faq  = ( Kstar * R_DRY * T(I,J,L) * LWC )
       Faq  = Faq / ( 1e+0_fp + Faq )

       RH3 = State_Met%RH(I,J,L)

       ! Cl- in sea-salt aerosol enhances solubility 2000X in MBL
       IF( LRED_JNO2 .AND. (F_PBL >0.1) .AND. State_Met%IsWater(I,J)) THEN
          Faq = ( Kstar * 2e+3_fp * R_DRY * T(I,J,L) * LWC )
          Faq = Faq / ( 1e+0_fp + Faq )
       ENDIF

       ! Fraction of all Hg(II) in the grid box which is in cloud water.
       ! We have used the in-cloud LWC to compute the sulfate
       ! aqueous chemistry.  We get the appropriate grid-box averaged
       ! mass of SO2 and sulfate by multiplying the reaction rates
       ! by the cloud fraction after the aqueous chemistry
       ! has been done.  (lzh, jaf, bmy, 5/27/11)
       ! Adopted for mercury chemistry (yzh, 11/9/2011)
       IF ( LRED_CLOUDONLY ) THEN
          Faq =  Faq * FC
       ENDIF

       ! Define K for the reduction reaction.
       IF (LRED_JNO2) THEN
          ! Reduction happens to aqueous fraction of Hg2
          K_RED = K_RED_JNO2 * Faq &
                  * GET_JNO2( I, J, L, State_Grid, State_Met ) &
                  * C_OA

          IF (RH3 > 35.0e0_fp) THEN
             ! Reduction of HgP (35% RH is for presence of water)
             K_RED2P = K_RED_JNO2 &
                       * GET_JNO2( I, J, L, State_Grid, State_Met ) &
                       * C_OA
          ELSE
             K_RED2P = 0e0_fp
          ENDIF

       ELSE
          ! Include the fraction of
          ! Hg(II) in air within the Kreduction &
          ! scale to OH concentration [/s]
          K_RED =   K_RED_OH * Faq * C_OH
          K_RED2P = K_RED
       ENDIF

       !---------------------------------------------
       ! Round off small chemical rates
       !---------------------------------------------

       ! Round very small K_OX, K_RED to prevent numerical errors
       ! (jaf, cdh, 11/17/11)
       IF ( K_OX  < 1e-10_fp )   K_OX    = 0e+0_fp
       IF ( K_RED < 1e-10_fp )   K_RED   = 0e+0_fp
       IF ( K_RED2P < 1e-10_fp ) K_RED2P = 0e+0_fp

       !==============================================================
       ! CHEMICAL SOLVER
       !
       ! The chemical ODEs are
       ! dX/dt = A * X,
       ! where X is a column vector consisting of
       ! X = [Hg(0), RGM, PBM]
       ! plus additional terms for total deposition, oxidation,
       ! and reduction of each species.
       !
       ! We use an explicit 4th-order Runge-Kutta solver.
       ! The internal time step is chosen to be 5 times smaller than
       ! the shortest chemical or deposition time scale.
       ! Since the shortest time scales are usually a few hours for Hg,
       ! there is rarely more than 1 internal time step per hour.
       ! (C. Holmes 4/2/2013)
       !==============================================================

       ! Largest rate, 1/s
       KMAX = K_OX + K_RED + K_RED2P + K_DEP0 + K_DEP2G + K_DEP2P

       ! Number of internal time steps
       ! Choose so there are 5 steps within the shortest timescale
       ! but not fewer than 1.
       NSTEPS = MAX( CEILING( 5e+0_fp * KMAX * DTCHEM ), 1 )

       ! Internal timestep, s
       DT = DTCHEM / NSTEPS

       !------------------------------------------------
       ! Load coefficients into the A = d/dt matrix
       !------------------------------------------------

       ! Reset elements to zero
       A(:,:) = 0e+0_fp

       ! Row 1: d(Hg0)/dt = -(K_OX+K_DEP0) * Hg0 + K_RED * (Hg2g+Hg2p)
       A(1,1) = - K_DEP0 - K_OX
       A(1,2) = + K_RED
       A(1,3) = + K_RED2P

       ! Row 2: d(Hg2g)/dt = K_OX * Hg0  - (K_RED + K_DEP2G) * Hg2g
       A(2,1) = + K_OX
       A(2,2) = - K_DEP2G - K_RED

       ! Row 3: d(Hg2p)/dt = -(K_DEP2P + K_RED2P) * Hg2p
       A(3,3) = - K_DEP2P - K_RED2P

       ! Row 4-6 accumulate dry dep of Hg0, Hg2g, Hg2p, respectively
       A(4,1) = K_DEP0
       A(5,2) = K_DEP2G
       A(6,3) = K_DEP2P

       ! Row 7 accumulate gross oxidation of Hg0
       A(7,1) = K_OX

       ! Row 8-9 accumulate gross reduction of Hg2g, Hg2p, respectively
       A(8,2) = K_RED
       A(9,3) = K_RED2P

       !------------------------------------------------
       ! Run solver for each Hg regional tag
       !------------------------------------------------
       DO N = 1, N_HG_CATS

          ! Load initial concentrations into vector, kg
          Xold(:) = 0e+0_fp
          Xold(1) = MAX( Spc(I,J,L,Hg0_Id_List(N)), SMALLNUM )
          Xold(2) = MAX( Spc(I,J,L,Hg2_Id_List(N)), SMALLNUM )
          Xold(3) = MAX( Spc(I,J,L,HgP_Id_List(N)), SMALLNUM )
          ! Rows 4-9 accumulate oxidation, reduction, and deposition fluxes
          ! so these start at zero

          ! Loop over internal time steps
          DO K=1, NSTEPS

             !-----------------------------------------
             ! Runge-Kutta equations
             !-----------------------------------------
             s1 = matmul(A,Xold)
             s2 = matmul(A,Xold+DT/2e+0_fp*s1)
             s3 = matmul(A,Xold+DT/2e+0_fp*s2)
             s4 = matmul(A,Xold+DT*s3)
             Xnew = Xold + DT/6e+0_fp * (s1 + 2*s2 + 2*s3 + s4 )

             ! Make sure concentrations are not negative
             IF ( (Xnew(1) < 0e+0_fp) .or. (Xnew(2)<0) .or. (Xnew(3)<0) ) THEN
                WRITE(6,101) I,J,L,K,NSTEPS, &
                             Xnew(1), Xnew(2), Xnew(3), &
                             Xold(1), Xold(2), Xold(3), &
                             K_OX, K_RED, K_RED2P, K_DEP0, K_DEP2G, K_DEP2P
101             FORMAT( 'Negative Hg in Box (I,J,L):', 3I3, &
                        ', step',I3,' out of ',I3, &
                        ', Current conc: ',3E10.3, &
                        ', Initial conc: ',3E10.3, &
                        ', Rate coeffs: ',7E10.3 )
                CALL ERROR_STOP( 'Negative Hg concentration', &
                                 ' MERCURY_MOD: CHEMMERCURY' )
             ENDIF

             ! Set for next loop
             Xold = Xnew

          ENDDO

          !--------------------------------------------------
          ! Extract final concentrations and fluxes
          ! from solution vector
          !--------------------------------------------------

          ! Archive Concentrations, kg/box
          Spc(I,J,L,Hg0_Id_List(N)) = Xnew(1)
          Spc(I,J,L,Hg2_Id_List(N)) = Xnew(2)
          Spc(I,J,L,HgP_Id_List(N)) = Xnew(3)

          ! Accumulated fluxes, kg/box/timestep [timestep=DTCHEM]
          DEP_HG0  = Xnew(4)
          DEP_Hg2g = Xnew(5)
          DEP_Hg2p = Xnew(6)
          GROSS_OX = Xnew(7)
          GROSS_RED= Xnew(8)+Xnew(9)

          ! Apportion gross oxidation between oxidants [kg]
          IF ( (K_OX   < SMALLNUM) .OR. (GROSS_OX < SMALLNUM) ) THEN
             GROSS_OX_OH    = 0e+0_fp
             GROSS_OX_BR    = 0e+0_fp
             GROSS_OX_O3    = 0e+0_fp
             GROSS_OX_CLY   = 0e+0_fp
             GROSS_OX_BRY   = 0e+0_fp
             GROSS_OX_BR2   = 0e+0_fp
             GROSS_OX_BRBRO = 0e+0_fp
             GROSS_OX_BRCLO = 0e+0_fp
             GROSS_OX_BROH  = 0e+0_fp
             GROSS_OX_BRHO2 = 0e+0_fp
             GROSS_OX_BRNO2 = 0e+0_fp
          ELSE
             GROSS_OX_OH  = GROSS_OX * ( K_OH + K_OH_AQ ) / K_OX
             GROSS_OX_BR  = GROSS_OX * ( K_BR + K_BRO   ) / K_OX
             GROSS_OX_O3  = GROSS_OX * ( K_O3 + K_O3_AQ ) / K_OX
             GROSS_OX_CLY = GROSS_OX * ( K_CL           ) / K_OX
             GROSS_OX_BRY = GROSS_OX * ( K_BR2  + K_BRBRO + K_BRCLO + &
                                         K_BROH + K_BRHO2 + K_BRNO2 ) / K_OX
             ! specific hgbr products
             GROSS_OX_BR2   = GROSS_OX * ( K_BR2   ) / K_OX
             GROSS_OX_BRBRO = GROSS_OX * ( K_BRBRO ) / K_OX
             GROSS_OX_BRCLO = GROSS_OX * ( K_BRCLO ) / K_OX
             GROSS_OX_BROH  = GROSS_OX * ( K_BROH  ) / K_OX
             GROSS_OX_BRHO2 = GROSS_OX * ( K_BRHO2 ) / K_OX
             GROSS_OX_BRNO2 = GROSS_OX * ( K_BRNO2 ) / K_OX
          ENDIF

          ! Apportion deposition between dry deposition and sea salt [kg]
          IF ( (K_DEP2G    < SMALLNUM) .OR. (DEP_Hg2g  < SMALLNUM) ) THEN
             DEP_HG2g_SALT = 0e+0_fp
             DEP_HG2g_DRY  = 0e+0_fp
          ELSE
             DEP_HG2g_DRY  = DEP_HG2g * K_DRYD2G / K_DEP2G
             DEP_HG2g_SALT = DEP_HG2g - DEP_HG2G_DRY
          ENDIF

          ! We currently assume that Hg2p is not scavenged by sea-salt
          ! aerosol, so all Hg2p deposition is dry deposition
          DEP_Hg2p_DRY = DEP_Hg2p

          !=================================================================
          ! Add deposited Hg(II) to the ocean module. OCEAN_MERCURY_MOD
          ! determines whether the box is marine, so we don't need to here.
          ! We should add an if statement to test whether DYNAMIC LAND is
          ! active.
          !
          ! IMPORTANT NOTE: DEP_Hg2G, DEP_Hg2P, and DEP_Hg2_DRY are defined
          ! on each 3-D (I,J,K) iteration.  But in routines ADD_Hg2_DD,
          ! ADD_Hg2_WD and ADD_Hg2_SNOWPACK, these are saved to diagnostic
          ! arrays that have only 2-D (I,J) spatial size.  This can cause
          ! slight numerical differences when OpenMP parallelization is
          ! turned on, because more than one CPU is trying to write to the
          ! diagnostic arrays simultaneously.  To avoid this situation, we
          ! call  ADD_Hg2_WD and ADD_Hg2_SNOWPACK from within an
          ! !$OMP CRITICAL block, which ensures that only one CPU at a time
          ! can write to the diagnostic arrays. (bmy, 4/20/16)
          !=================================================================
          !$OMP CRITICAL
          ! Archive dry-deposited Hg2
          CALL ADD_Hg2_DD( I, J, N, DEP_Hg2G )
          ! Archive dry-deposited HgP for Hg category # N
          CALL ADD_HgP_DD( I, J, N, DEP_Hg2P )

          ! Add deposited Hg(II) to the snowpack
          IF ( LHGSNOW ) THEN
             DEP_HG2_DRY =  DEP_HG2g_DRY + DEP_HG2p_DRY
             CALL ADD_HG2_SNOWPACK( I, J, N, DEP_HG2_DRY, &
                                    State_Met, State_Chm, State_Diag )
          ENDIF
          !$OMP END CRITICAL

#ifdef BPCH_DIAG
          !==============================================================
          ! %%%%% ND03 (bpch) DIAGNOSTIC %%%%%
          !
          ! Hg(II) production
          ! Concentration of Br and BrO
          ! Loss of Hg2 by seasalt
          !==============================================================
          IF ( ND03 > 0 .AND. L <= LD03 ) THEN

             ! Store chemistry diagnostics only for total species
             IF ( N == 1 ) THEN

                ! Net oxidation [kg]
                NET_OX = GROSS_OX - GROSS_RED

                ! Production of Hg2 from Hg0 [kg]
                AD03_Hg2_Hg0(I,J,L,N)  = AD03_Hg2_Hg0(I,J,L,N) + NET_OX

                ! Production of Hg2 from Br [kg]
                AD03_Hg2_Br(I,J,L,N)   = AD03_Hg2_Br(I,J,L,N) + GROSS_OX_BR

                ! Production of Hg2 from OH [kg]
                AD03_Hg2_OH(I,J,L,N)   = AD03_Hg2_OH(I,J,L,N) + GROSS_OX_OH

                ! Production of Hg2 from O3 [kg]
                AD03_Hg2_O3(I,J,L,N)   = AD03_Hg2_O3(I,J,L,N) + GROSS_OX_O3

                ! Production of Hg2 from BrY [kg]
                AD03_Hg2_BRY(I,J,L,N)  = AD03_Hg2_BRY(I,J,L,N) + GROSS_OX_BRY

                ! Production of Hg2 from ClY [kg]
                AD03_Hg2_CLY(I,J,L,N)  = AD03_Hg2_CLY(I,J,L,N) + GROSS_OX_CLY

                ! Production of Hg2 from Br2 [kg]
                AD03_Hg2_Br_Y(I,J,L,1) = AD03_Hg2_Br_Y(I,J,L,1) + GROSS_OX_BR2

                ! Production of Hg2 from BrBrO [kg]
                AD03_Hg2_Br_Y(I,J,L,2) = AD03_Hg2_Br_Y(I,J,L,2) + GROSS_OX_BRBRO

                ! Production of Hg2 from BrHO2 [kg]
                AD03_Hg2_Br_Y(I,J,L,3) = AD03_Hg2_Br_Y(I,J,L,3) + GROSS_OX_BRHO2

                ! Production of Hg2 from BrNO2 [kg]
                AD03_Hg2_Br_Y(I,J,L,4) = AD03_Hg2_Br_Y(I,J,L,4) + GROSS_OX_BRNO2

                ! Production of Hg2 from BrClO [kg]
                AD03_Hg2_Br_Y(I,J,L,5) = AD03_Hg2_Br_Y(I,J,L,5) + GROSS_OX_BRCLO

                ! Production of Hg2 from BrOH [kg]
                AD03_Hg2_Br_Y(I,J,L,6) = AD03_Hg2_Br_Y(I,J,L,6) + GROSS_OX_BROH

                ! Concentrations of Br and BrO [molec/cm3]
                AD03_Br(I,J,L,1) = AD03_Br(I,J,L,1) + C_BR
                AD03_Br(I,J,L,2) = AD03_Br(I,J,L,2) + C_BRO

             ENDIF

             ! Sea salt diagnostic is 2-D [kg]
             AD03_Hg2_SS(I,J,N) = AD03_Hg2_SS(I,J,N) + DEP_HG2G_SALT

             ! Sea-salt loss rate diagnostic [/s]
             IF ( L == 1 ) THEN
                AD03_Hg2_SSR(I,J,N) = AD03_Hg2_SSR(I,J,N) + K_SALT
             ENDIF

          ENDIF
#endif

          !===========================================================
          ! %%%%% HISTORY (aka netCDF diagnostics) %%%%%
          !
          ! Dry deposition flux of Hg(II) in molec/cm2/s
          !
          ! NOTE: This is only updated when using the full PBL mixing
          ! option (aka TURBDAY).  When you use the non-local PBL
          ! mixing, option (aka VDIFF), the dry deposition fluxes are
          ! applied at the same time the PBL mixing is done.
          !===========================================================
          IF ( ( .not. LNLPBL ) .and. ( F_PBL > 0.1_fp ) ) THEN

             ! Grid box surface area [cm2]
             AREA_CM2 = State_Grid%Area_M2(I,J) * 1e+4_fp

             ! Archive Hg(0) drydep flux [molec/cm2/s]
             IF ( State_Diag%Archive_DryDepChm .or. &
                  State_Diag%Archive_DryDep    ) THEN

                ! Hg0 species index
                NN = Hg0_Id_List(N)

                ! Amt of Hg(0) lost to drydep [molec/cm2/s]
                DEP_DRY_FLX  = DEP_HG0 * AVO / ( 1.e-3_fp * &
                               State_Chm%SpcData(NN)%Info%emMW_g ) / &
                               ( AREA_CM2 * DTCHEM )

                ! Archive to State_Diag (sum over levels)
                State_Diag%DryDepChm(I,J,NN) = State_Diag%DryDepChm(I,J,NN) + &
                                               DEP_DRY_FLX
             ENDIF

             ! Archive Hg(II) drydep flux [molec/cm2/s]
             IF ( State_Diag%Archive_DryDepChm .or. &
                  State_Diag%Archive_DryDep    ) THEN

                ! Hg2 species index
                NN = Hg2_Id_List(N)

                ! Amt of Hg(II)g lost to drydep [molec/cm2/s]
                DEP_DRY_FLX  = DEP_HG2G_DRY * AVO / ( 1.e-3_fp * &
                               State_Chm%SpcData(NN)%Info%emMW_g ) / &
                               ( AREA_CM2 * DTCHEM )

                ! Archive to State_Diag (sum over levels)
                State_Diag%DryDepChm(I,J,NN) = State_Diag%DryDepChm(I,J,NN) + &
                                                DEP_DRY_FLX
             ENDIF

             ! Archive Hg(0) drydep flux [molec/cm2/s]
             IF ( State_Diag%Archive_DryDepChm .or. &
                  State_Diag%Archive_DryDep    ) THEN

                ! HgP species index
                NN = HgP_Id_List(N)

                ! Amt of Hg(II)p lost to drydep [molec/cm2/s]
                DEP_DRY_FLX  = DEP_HG2P_DRY * AVO / ( 1.e-3_fp * &
                               State_Chm%SpcData(NN)%Info%emMW_g ) / &
                               ( AREA_CM2 * DTCHEM )

                ! Archive to State_Diag (sum over levels)
                State_Diag%DryDepChm(I,J,NN) = State_Diag%DryDepChm(I,J,NN) + &
                                               DEP_DRY_FLX
             ENDIF
          ENDIF

          !===========================================================
          ! %%%%% HISTORY (aka netCDF diagnostics) %%%%%
          !
          ! Hg(II) production
          ! Concentration of Br and BrO
          ! Loss of Hg2 by seasalt
          !
          ! Accorrding to Colin Thackray, these are Hg2 diagnostics,
          ! that should only get updated once per timestep, for the
          ! total Hg2 species (i.e. category N==1).  Therefore, these
          ! diagnostics will not archive the tagged Hg2 species.
          !===========================================================
          IF ( N == 1 ) THEN

             ! Net oxidation [kg]
             NET_OX = GROSS_OX - GROSS_RED

             ! Concentrations of Br [molec/cm3]
             IF ( State_Diag%Archive_ConcBr ) THEN
                State_Diag%ConcBr(I,J,L) = C_Br
             ENDIF

             ! Concentrations of BrO [molec/cm3]
             IF ( State_Diag%Archive_ConcBrO) THEN
                State_Diag%ConcBrO(I,J,L) = C_BrO
             ENDIF

             ! Production of Hg2 from Br [kg/s]
             IF ( State_Diag%Archive_ProdHg2fromHg0 ) THEN
                State_Diag%ProdHg2fromBr(I,J,L) = GROSS_OX_BR / DTCHEM
             ENDIF

             ! Production of Hg2 from BrY [kg/s]
             IF ( State_Diag%Archive_ProdHg2fromBrY ) THEN
                State_Diag%ProdHg2fromBrY(I,J,L) = GROSS_OX_BRY / DTCHEM
             ENDIF

             ! Production of Hg2 from ClY [kg/s]
             IF ( State_Diag%Archive_ProdHg2fromClY ) THEN
                State_Diag%ProdHg2FromClY(I,J,L) = GROSS_OX_CLY / DTCHEM
             ENDIF

             ! Production of Hg2 from Hg0 [kg/s]
             IF ( State_Diag%Archive_ProdHg2fromHg0 ) THEN
                State_Diag%ProdHg2fromHg0(I,J,L) = NET_OX / DTCHEM
             ENDIF

             ! Production of Hg2 from O3 [kg/s]
             IF ( State_Diag%Archive_ProdHg2fromO3 ) THEN
                State_Diag%ProdHg2fromO3(I,J,L) = GROSS_OX_O3 / DTCHEM
             ENDIF

             ! Production of Hg2 from OH [kg/s]
             IF ( State_Diag%Archive_ProdHg2fromOH ) THEN
                State_Diag%ProdHg2fromOH(I,J,L) = GROSS_OX_OH / DTCHEM
             ENDIF

             ! Production of Hg2 from HgBr + Br2 [kg/s]
             IF ( State_Diag%Archive_ProdHg2fromHgBrPlusBr2 ) THEN
                State_Diag%ProdHg2fromHgBrPlusBr2(I,J,L) = GROSS_OX_BR2 / &
                                                           DTCHEM
             ENDIF

             ! Production of Hg2 from HgBr + BrBrO [kg/s]
             IF ( State_Diag%Archive_ProdHg2fromHgBrPlusBrBrO ) THEN
                State_Diag%ProdHg2fromHgBrPlusBrBrO(I,J,L) = GROSS_OX_BRBRO / &
                                                             DTCHEM
             ENDIF

             ! Production of Hg2 from HgBr + BrClO [kg/s]
             IF ( State_Diag%Archive_ProdHg2fromHgBrPlusBrClO ) THEN
                State_Diag%ProdHg2fromHgBrPlusBrClO(I,J,L) = GROSS_OX_BRCLO / &
                                                             DTCHEM
             ENDIF

             ! Production of Hg2 from HgBr + BrHO2 [kg/s]
             IF ( State_Diag%Archive_ProdHg2fromHgBrPlusBrHO2 ) THEN
                State_Diag%ProdHg2fromHgBrPlusBrHO2(I,J,L) = GROSS_OX_BRHO2 / &
                                                             DTCHEM
             ENDIF

             ! Production of Hg2 from HgBr + BrNO2 [kg/s]
             IF ( State_Diag%Archive_ProdHg2fromHgBrPlusBrNO2 ) THEN
                State_Diag%ProdHg2fromHgBrPlusBrNO2(I,J,L) = GROSS_OX_BRNO2 / &
                                                             DTCHEM
             ENDIF

             ! Production of Hg2 from HgBr + BrOH [kg/s]
             IF ( State_Diag%Archive_ProdHg2fromHgBrPlusBrOH ) THEN
                State_Diag%ProdHg2fromHgBrPlusBrOH(I,J,L) = GROSS_OX_BROH / &
                                                            DTCHEM
             ENDIF

             ! Loss Hg2 by sea salt aerosol
             ! NOTE: Sum contributions from PBL into the diagnostic
             IF ( State_Diag%Archive_LossHg2bySeaSalt ) THEN
                State_Diag%LossHg2bySeaSalt(I,J,L) = DEP_HG2G_SALT / DTCHEM
             ENDIF

             ! Sea-salt loss rate diagnostic [/s]
             IF ( L == 1 ) THEN
                IF ( State_Diag%Archive_LossRateHg2bySeaSalt ) THEN
                   State_Diag%LossRateHg2bySeaSalt(I,J) = K_SALT
                ENDIF
             ENDIF

          ENDIF
       ENDDO

    ENDDO
    ENDDO
    ENDDO
    !$OMP END PARALLEL DO

    ! Partition Hg(II) between gas and aerosol
    CALL PARTITIONHG2( Input_Opt, State_Chm, State_Grid, RC )

    ! Free pointer memory
    Spc     => NULL()
    T       => NULL()
    SpcInfo => NULL()

    IF ( prtDebug ) CALL DEBUG_MSG( 'CHEMMERCURY: a CHEM_HgP' )

  END SUBROUTINE CHEMMERCURY
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
! !IROUTINE: get_hgaq_rate
!
! !DESCRIPTION: Function GET\_HGAQ\_RATE computes the effective 1st order
!  conversion rate of Hg(0) to Hg(II) via aqueous-phase oxidation in cloud
!  droplets  for three reactions: HOCl, O3, and OH (2-step). hmh 11/19/14
!\\
!\\
! !INTERFACE:
  FUNCTION GET_HGAQ_RATE( I, J, L, BR, OH, O3, HOCL, State_Met ) &
       RESULT( K_AQ_OX )
!
! !USES:
!
    USE State_Met_Mod, ONLY : MetState
    USE ERROR_MOD,     ONLY : SAFE_DIV
!
! !INPUT PARAMETERS:
!
    INTEGER,            INTENT(IN)  :: I, J, L
    REAL(fp),           INTENT(IN)  :: BR          ! Atomic Br [molec/cm3]
    REAL(fp),           INTENT(IN)  :: OH          ! OH        [molec/cm3]
    REAL(fp),           INTENT(IN)  :: O3          ! O3        [molec/cm3]
    REAL(fp),           INTENT(IN)  :: HOCL        ! HOCl      [molec/cm3]
    TYPE(MetState),     INTENT(IN)  :: State_Met   ! Meteorology State obj
!
! !RETURN VALUE:

    ! Effective 1st order loss rate of Hg(0) [1/s] for:
    ! reaction with aqueous (1) HOCl, (2) OH, and (3) O3
    REAL(fp)                        :: K_AQ_OX(3)
!
!  References:
!  ============================================================================
!  1. Wang and Pehkonen, 2004
!  2. Pehkonen and Lin, 1997
!  3. Buxton et al., 1988
!  4. Munthe, 1992
!
!  !REVISION HISTORY:
!  19 November 2014 - H. Horowitz (HMH) - initial version finalized
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:

    REAL(fp)               :: KH_Hg0, KH_HOCl, KH_O3, delta_OH
    REAL(fp)               :: FC,     LWC
    REAL(fp)               :: F_Hg0,  F_HOCl,  F_O3
    REAL(fp)               :: R_HOCl, R_OH,    R_O3
    REAL(fp)               :: Nair, T, P
!
! !DEFINED PARAMETERS:
!
    ! Gas constant
    REAL(fp), PARAMETER    :: R_L_atm = 8.2e-2_fp ! [L atm /K /mol]
    REAL(fp), PARAMETER    :: R_m3_Pa = 8.31e0_fp ! [m3 Pa /K /mol]

    ! K for aqueous reactions  [L /mol /s]
    REAL(fp), PARAMETER    :: kHOCl   = 2.0e6_fp  ! Source: Wang and Pehkonen
    REAL(fp), PARAMETER    :: kOH     = 2.0e9_fp  ! Source: Lin and Pehkonen
    REAL(fp), PARAMETER    :: kO3     = 4.7e7_fp  ! Source: Munthe, 1992

    !=================================================================
    ! GET_HGAQ_RATE begins here!
    !=================================================================

    T = State_Met%T(I,J,L)        ! Temperature [K]
    P = State_Met%PMID_DRY(I,J,L) ! Pressure [hPa]

    !==============================================================
    ! Calculate cloud parameters for aqueous Hg chemistry
    !==============================================================

    ! Get cloud fraction from met fields
    FC      = State_Met%CLDF(I,J,L)

    ! Get liquid water content [m3 H2O/m3 air] within cloud from met flds
    ! Units: [kg H2O/kg air] * [kg air/m3 air] * [m3 H2O/1e3 kg H2O]
    LWC   = State_Met%QL(I,J,L) * State_Met%AIRDEN(I,J,L) * 1e-3_fp

    ! LWC is a grid-box averaged quantity. To improve the representation
    ! of sulfate chemistry, we divide LWC by the cloud fraction and
    ! compute sulfate chemistry based on the LWC within the cloud.  We
    ! get the appropriate grid-box averaged mass of SO2 and sulfate by
    ! multiplying these quantities by FC AFTER computing the aqueous
    ! sulfur chemistry within the cloud. (lzh, jaf, bmy, 5/27/11)
    LWC     = SAFE_DIV( LWC, FC, 0e0_fp )

    ! Henry's law constants (M per atm)
    ! Hg(0):
    KH_Hg0 = 1.28e-1_fp * exp(2482e0_fp * ((1/T) - 1/298.15e0_fp))
    ! Toyota et al. 2014 citing Sanemasa 1975
    ! or KH_Hg0 = 1.1E-1 ! Hynes et al. 2009

    ! HOCl:
    KH_HOCl = 6.6e2_fp * exp(5900e0_fp * ((1/T)-1/298.15e0_fp))
    ! NIST webbook, original citation: Huthwelker et al 1995
    ! (http://www.mpch-mainz.mpg.de/~sander/res/henry.html)
    ! fairly close to average value of other studies:
    ! others: 7.3E2 ; 4.8E2; 9.3E2; 2.6E2*exp(5100.*(1/T-1/298.15))

    ! O3:
    KH_O3 = 1.1e-2_fp * exp(2400e0_fp * ((1/T) - 1/298.15e0_fp))
    ! Jacob, 1986
    ! representative of average of values from 10 papers cited in NIST

    ! OH:
    ! assume a simple parameterization from Jacob et al. (2005) where
    ! [OH(aq)] = delta*[OH(g)] where [OH(g)] is the gas-phase concentrati
    ! calculated in GEOS-Chem without consideration of aqueous-phase cloud
    ! chemistry, and delta = 1E-19 M cm3 molecule-1 is chosen to fit the
    ! cloud chemistry model results of Jacob (1986).
    delta_OH = 1e-19_fp ! M cm3 molecule-1

    ! Aqueous oxidation reactions for Hg: HOCl and OCl
    ! Both reactions have similar rate constants within 5%.
    ! Therefore we will not explicitly model the dissociation of HOCl
    ! once dissolved in the cloud droplet, but group ([HOCl]+[OCl-])
    ! as equal to the initial dissolution of HOCl(g) based on Henry's
    ! Law.

    ! (1) Hg0(aq) + HOCl(aq) -> Hg(2+) + OH- + Cl- , k1=kHOCl

    ! rate expression: R_HOCl = kHOCl[Hg0(aq)][HOCl(aq)]
    ! we can rewrite [Hg0(aq)] and [HOCl(aq)] in terms of the partial
    ! equilibrium pressure and the Henry's law constants:
    ! kHOCl * pHg0 * KH_Hg0 * pHOCl * KH_HOCl

    ! to get the partial pressure applied only over the cloud droplets:
    ! Let us define a fraction of gas phase of A species
    ! in equilibrium with aqueous phase as
    !
    !        F_A  = 1/(1+f),
    !
    ! where  f   = hA * R * T * LWC,
    !        hA  = Henry's constant,
    !        R   = gas constant (in terms of atm not Pa!),
    !        T   = temperature in kelvin,
    !        LWC = liquid water content [m3/m3] of cloud

    ! Then, and recall that pX = c_x * P
    !   c_x is the mixing ratio in mol/mol
    !   P is the total pressure
    ! c_x = n_x / n_a
    !   n_x is the number density of x (in molec/cm3), aka [X(g)]
    !   n_a is the number density of air in the same units.
    ! n_a = P*1E2*Av/(RT)*1E-6
    !   P is the pressure in hPa
    !   Av is avogadro's number to convert from moles to molecules
    !   R is the gas constant (8.314)
    !   factors to convert hPa to Pa (1E2) and m-3 to cm-3 (1E-6).

    ! now our rate expression is as follows (converting P from hPa to atm)
    ! R_HOCl = k * pHg0              * KH_Hg0 * pHOCl       * KH_HOCl
    ! R_HOCl = k * (P/1013) * c_Hg0 * F_Hg0 * KH_Hg0 * (P/1013) * c_HOCl *
    !          F_HOCl* KH_HOCl
    ! for Hg0, similar to HOCl:
    ! c_Hg0 = nHg0/na = [Hg0(g)]/(P*1E2*Av/(RT)*1E-6)
    ! F_Hg0 = 1/(1+KH_Hg0*R*T*LWC)

    ! This rate expression R_HOCl gives us the rate of Hg0 oxidation in
    ! units of M / s, or moles per L of cloudwater / s.
    ! To convert back to units of molecules/cm3 of air/s, we need to use
    ! LWC of cloudwater and the cloud fraction:
    ! R_HOCl(mol/L cloud H2O/s)* LWC (m3 cloud H2O/m3 air in cloud)* FC
    !(vol cloud/vol air)
    !       * Av * 1000 L H2O/m3 H2O * (1E-6 m3 air/cm3 air) = rate in
    !molec/cm3/s

    ! calculate fractions
    F_Hg0  = 1e0_fp/(1e0_fp+KH_Hg0 *R_L_atm*T*LWC) ! unitless
    F_HOCl = 1e0_fp/(1e0_fp+KH_HOCl*R_L_atm*T*LWC) ! unitless
    ! or: F_Hg0 = 1+KH_Hg0*R_L_atm*T*LWC -> F_Hg0 = 1/(1+F_Hg0)

    ! calculate number density of air in molecules per cm3:
    Nair = P * 1e2_fp / (R_m3_Pa * T ) * 6.02e23_fp / 1e6_fp

    ! calculate rate of reaction (without Hg0 multiplied) in cloudwater
    R_HOCl = kHOCl * P/1013e0_fp * 1e0_fp/Nair * F_Hg0  * KH_Hg0 * &
                     P/1013e0_fp * HOCl / Nair * F_HOCl * KH_HOCl

    ! calculate first-order effective rate constant w.r.t. Hg0 for volume
    !of air:
    K_AQ_OX(1) = R_HOCl * LWC * FC * 6.02e23_fp / 1e3_fp

    !------------------------!
    ! OH aqueous             !
    !------------------------!
    !(1) Hg(0)(aq) + OH(aq) -> Hg+ + OH-     , k1=kOH
    !(2) Hg+       + OH(aq) -> Hg(2+) + OH-  , k2

    ! assume Hg+ is in steady state, then the rate of Hg(2+) production
    ! aqueous oxidation with OH is:
    ! R_OH = k1[Hg(0)(aq)][OH(aq)]

    ! we calculate the rate similarly, except that [OH(aq)] is not depende
    ! on Henry's law but on the delta factor (explained above)

    ! calculate rate of reaction (without Hg0 multiplied) in cloudwater
    R_OH = kOH * P/1013e0_fp * 1e0_fp/Nair * F_Hg0 * KH_Hg0 * delta_OH * OH

    ! calculate first-order effective rate constant w.r.t. Hg0 for volume
    ! of air
    K_AQ_OX(2) = R_OH * LWC * FC * 6.02e23_fp / 1e3_fp
    !------------------------!
    ! O3 aqueous             !
    !------------------------!
    !(1) Hg(0)(aq) + O3(aq) -> Hg(2+) , k1=kO3

    ! calculate the rate similarly to HOCl reaction but with O3:
    ! calculate fraction
    F_O3 = 1e0_fp/(1e0_fp+KH_O3*R_L_atm*T*LWC) ! unitless

    ! calculate rate of reaction (without Hg0 multiplied) in cloudwater
    R_O3 = kO3 *P/1013e0_fp *1e0_fp/Nair * F_Hg0 * KH_Hg0 *P/1013e0_fp &
                            * O3   /Nair * F_O3  * KH_O3

    ! calculate first-order effective rate constant w.r.t. Hg0 for volume
    !of air:
    K_AQ_OX(3) = R_O3 * LWC * FC * 6.02e23_fp / 1e3_fp

  END FUNCTION GET_HGAQ_RATE
!EOC

!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_hgbr_rate
!
! !DESCRIPTION: Function GET\_HGBR\_RATE computes the effective 1st order
!  conversion rate of Hg(0) to Hg(II) via two-step recombination with Br
!  and OH.
!\\
!\\
! !INTERFACE:
!
  FUNCTION GET_HGBR_RATE( I, J, L, BR, State_Met, OH, METHOD ) &
       RESULT(K_HGBR)
!
! USES:
!
    USE State_Met_Mod,      ONLY : MetState
!
! !INPUT PARAMETERS:
!
    INTEGER,            INTENT(IN)  :: I, J, L
    REAL(fp),           INTENT(IN)  :: BR           ! Atomic Br [molec/cm3]
    REAL(fp),           INTENT(IN)  :: OH           ! OH [molec/cm3]
    CHARACTER(LEN=*),   INTENT(IN)  :: METHOD       ! Set of rate constants to
                                                    !  use in calculation
    TYPE(MetState),     INTENT(IN)  :: State_Met    ! Meteorology State obj
!
! !RETURN VALUE:
!
    ! K_HgBr : Effective 1st order loss rate of Hg(0) [1/s]
    REAL(fp)                        :: K_HGBR
!
! !REMARKS:
!  ============================================================================
!  This subroutine calculates the net rate of Hg(0) oxidation to Hg(II) through
!     the following reactions. All are gas phase.
!                                                                             .
!     (1  )  Hg(0) + Br -> HgBr
!     (2  )  HgBr + M -> Hg(0) + Br
!     (2a )  HgBr     -> Hg(0) + Br
!     (3Br)  HgBr + Br -> HgBr2
!     (3OH)  HgBr + OH -> HgBrOH
!                                                                             .
!  References:
!  ============================================================================
!  1. Ariya, P. A., A. Khalizov, and A. Gidas (2002), Reaction of gaseous
!     mercury with atomic and molecular halogens: kinetics, product studies,
!     and atmospheric implications, Journal of Physical Chemistry A, 106,
!     7310-7320.
!                                                                             .
!  2. Balabanov, N. B., B. C. Shepler, and K. A. Peterson (2005), Accurate
!     global potential energy surface and reaction dynamics for the ground
!     state of HgBr2, Journal of Physical Chemistry A, 109(39), 8765-8773.
!                                                                             .
!  3. Donohoue, D. L., D. Bauer, B. Cossairt, and A. J. Hynes (2006),
!     Temperature and Pressure Dependent Rate Coefficients for the Reaction
!     of Hg with Br and the Reaction of Br with Br: A Pulsed Laser
!     Photolysis-Pulsed Laser Induced Fluorescence Study, Journal of
!     Physical Chemistry A, 110, 6623-6632.
!                                                                             .
!  4. Goodsite, M. E., J. M. C. Plane, and H. Skov (2004), A theoretical
!     study of the oxidation of Hg-0 to HgBr2 in the troposphere, Environmental
!     Science & Technology, 38(6), 1772-1776.
!                                                                             .
!  5. Holmes, C. D., et al. (2006), Global lifetime of elemental mercury
!     against oxidation by atomic bromine in the free troposphere, Geophys.
!     Res. Lett., 33(20).
!                                                                             .
!  6. Khalizov, A. F., B. Viswanathan, P. Larregaray, and P. A. Ariya (2003),
!     A theoretical study on the reactions of Hg with halogens: Atmospheric
!     implications, Journal of Physical Chemistry A, 107(33), 6360-6365.
!
! !REVISION HISTORY:
!  06 Jul 2006 - C. Holmes   - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    REAL(fp) :: k1,   k2, k2a, k3br, k3oh
    REAL(fp) :: k1G,  k2G
    REAL(fp) :: Nair, N0, T, P

    !=================================================================
    ! GET_HGBR_RATE begins here!
    !=================================================================

    T = State_Met%T(I,J,L)        ! Temperature [K]
    P = State_Met%PMID_DRY(I,J,L) ! Pressure [hPa]

    !number density of air [molec/cm3]
    Nair = P * 1e+2_fp / (AIRMW * Rd * T ) * AVO / 1e+3_fp

    !standard air density at STP: 1 atm, 273K [molec/cm3]
    N0 = ATM / (AIRMW * Rd * 273e+0_fp ) * AVO / 1e+3_fp

    SELECT CASE( METHOD )
    CASE( 'Goodsite' )
       !All rates from Goodsite et al. 2004
       !No HgBr+OH reaction
       k1   = 1.1e-12_fp * ( T / 298e+0_fp ) ** ( -2.37e+0_fp ) * ( Nair / N0 )
       k2   = 1.2e+10_fp  * exp( -8357e+0_fp / T )
       k2a  = 0e+0_fp
       k3br = 2.5e-10_fp * ( T / 298e+0_fp ) ** ( -.57e+0_fp )
       k3oh = 0e+0_fp

    CASE( 'GoodsiteY' )
       !All rates from Goodsite et al. 2004
       !Include HgBr+OH reaction
       k1   = 1.1e-12_fp * ( T / 298e+0_fp ) ** ( -2.37e+0_fp ) * ( Nair / N0 )
       k2   = 1.2e+10_fp  * exp( -8357e+0_fp / T )
       k2a  = 0e+0_fp
       k3br = 2.5e-10_fp * ( T / 298e+0_fp ) ** ( -.57e+0_fp )
       k3oh = k3br

    CASE( 'Donohoue' )
       !k1 from Donohoue et al. 2006
       !Other rates from Goodsite et al. 2004
       !No HgBr+OH reaction
       k1   = 1.46e-32_fp * ( T / 298e+0_fp ) ** (-1.86e+0_fp) * Nair
       k2   = 1.2e+10_fp  * exp( -8357e+0_fp / T )
       k2a  = 0e+0_fp
       k3br = 2.5e-10_fp * ( T / 298e+0_fp ) ** ( -.57e+0_fp )
       k3oh = 0e+0_fp

    CASE( 'DonohoueY' )
       !k1 from Donohoue et al. 2006
       !Other rates from Goodsite et al. 2004
       !Include HgBr+OH reaction
       k1   = 1.46e-32_fp * ( T / 298e+0_fp ) ** (-1.86e+0_fp) * Nair
       k2   = 1.2e+10_fp  * exp( -8357e+0_fp / T )
       k2a  = 0e+0_fp
       k3br = 2.5e-10_fp * ( T / 298e+0_fp ) ** ( -.57e+0_fp )
       k3oh = k3br

    CASE( 'DonohoueYB' )
       !k1 from Donohoue et al. 2006
       !k2 derived from Goodsite et al. 2004 rates
       !   preserving the k2/k1 equilibrium (detailed balance)
       !k3 from Goodsite et al. 2004
       !Include HgBr+OH reaction

       !Goodsite et al. 2004 rates
       k1G  = 1.1e-12_fp * ( T / 298e+0_fp ) ** ( -2.37e+0_fp ) * ( Nair / N0 )
       k2G  = 1.2e+10_fp  * exp(-8357e+0_fp/T)

       k1   = 1.46e-32_fp * ( T / 298e+0_fp ) ** (-1.86e+0_fp) * Nair
       k2   = k2G * k1 / k1G
       k2a  = 0e+0_fp
       k3br = 2.5e-10_fp * ( T / 298e+0_fp ) ** ( -.57e+0_fp )
       k3oh = k3br

    CASE( 'Balabanov' )
       !k1, k2 from Goodsite et al. 2004
       !k2a from Balabanov et al. 2005
       !k3 from the zero pressure limit in Balabanov et al. 2005
       k1   = 1.1e-12_fp * ( T / 298e+0_fp ) ** ( -2.37e+0_fp ) * ( Nair / N0 )
       k2   = 1.2e+10_fp  * exp( -8357e+0_fp / T )
       k2a  = 3.9e-11_fp
       k3br = 3e-11_fp
       k3oh = k3br

    CASE( 'KhalizovB' )
       !k1 from Khalizov et al. 2003
       !k2 derived from Goodsite et al. 2004 rates
       !   preserving the k2/k1 equilibrium (detailed balance)
       !k3 from Goodsite et al. 2004
       !Include HgBr+OH reaction

       !Goodsite et al. 2004 rates
       k1G  = 1.1e-12_fp * ( T / 298e+0_fp ) ** ( -2.37e+0_fp ) * ( Nair / N0 )
       k2G  = 1.2e+10_fp  * exp( -8357e+0_fp / T )

       k1   = 1.0e-12_fp * exp( 209e+0_fp / T ) * ( Nair / N0 )
       k2   = k2G * k1 / k1G
       k2a  = 0e+0_fp
       k3br = 2.5e-10_fp * ( T / 298e+0_fp ) ** ( -.57e+0_fp )
       k3oh = k3br

    CASE( 'AriyaB' )
       !k1 from Ariya et al. 2002
       !k2 derived from Goodsite et al. 2004 rates
       !   preserving the k2/k1 equilibrium (detailed balance)
       !k3 from Goodsite et al. 2004
       !Include HgBr+OH reaction

       !Goodsite et al. 2004 rates
       k1G  = 1.1e-12_fp * ( T / 298e+0_fp ) ** ( -2.37e+0_fp ) * ( Nair / N0 )
       k2G  = 1.2e+10_fp  * exp( -8357e+0_fp / T )

       k1   = 3.2e-12_fp * ( Nair / N0 )
       k2   = k2G * k1 / k1G
       k2a  = 0e+0_fp
       k3br = 2.5e-10_fp * ( T / 298e+0_fp ) ** ( -.57e+0_fp )
       k3oh = k3br

    CASE( 'DonohoueYBBalabanov' )
       !k1 from Donohoue et al. 2006
       !k2 derived from Goodsite et al. 2004 rates
       !   preserving the k2/k1 equilibrium (detailed balance)
       !k2a from Balabanov et al. 2005
       !k3 from Goodsite et al. 2004
       !Include HgBr+OH reaction

       !Goodsite et al. 2004 rates
       k1G  = 1.1e-12_fp * ( T / 298e+0_fp ) ** ( -2.37e+0_fp ) * ( Nair / N0 )
       k2G  = 1.2e+10_fp  * exp(-8357e+0_fp/T)

       k1   = 1.46e-32_fp * ( T / 298e+0_fp ) ** (-1.86e+0_fp) * Nair
       k2   = k2G * k1 / k1G
       k2a  = 3.9e-11_fp
       k3br = 2.5e-10_fp * ( T / 298e+0_fp ) ** ( -.57e+0_fp )
       k3oh = k3br

    CASE( 'GoodsiteUpdate' )
       !k1 from Donohoue et al. 2006
       !k2 derived from Goodsite et al. updated 2012 rates
       !   preserving the k2/k1 equilibrium (detailed balance)
       !k2a from Balabanov et al. 2005
       !k3 from Goodsite et al. 2004
       !Include HgBr+OH reaction

       !Goodsite et al. 2004 rates
       k1G  = 3.7e-13_fp * ( T / 298e+0_fp ) ** ( -2.76e+0_fp ) * ( Nair / N0 )
       k2G  = 4.0e+9_fp  * exp(-7292e+0_fp/T)

       k1   = 1.46e-32_fp * ( T / 298e+0_fp ) ** (-1.86e+0_fp) * Nair
       k2   = k2G * k1 / k1G
       k2a  = 3.9e-11_fp
       k3br = 2.5e-10_fp * ( T / 298e+0_fp ) ** ( -.57e+0_fp )
       k3oh = k3br

    CASE DEFAULT

    END SELECT

    IF ( BR > SMALLNUM ) THEN

       ! effective 1st order loss of Hg(0) [ 1/s ]
       K_HGBR = k1 * BR * (k3br * BR + k3oh * OH ) / &
              ( k2 + k2a * BR + k3br * BR + k3oh * OH )

    ELSE

       ! Avoid divide by zero in rate
       K_HGBR = 0

    ENDIF

  END FUNCTION GET_HGBR_RATE
!EOP
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
! !IROUTINE: get_hgcl_rate
!
! !DESCRIPTION: Function GET\_HGCL\_RATE computes the effective 1st order
!  conversion rate of Hg(0) to Hg(II) initiated by Cl and 2nd step with Cl,
!  OH, NO2, HO2, ClO, and BrO. hmh 11/21/14 copied from get\_hgbry\_rate.
!\\
!\\
! !INTERFACE:
  FUNCTION GET_HGCL_RATE( I, J, L, CL, OH, HO2, NO2, CLO, BRO, BR, State_Met ) &
       RESULT(K_CL)
!
! USES:
!
    USE State_Met_Mod,      ONLY : MetState
!
! !INPUT PARAMETERS:
!
    INTEGER,        INTENT(IN)  :: I, J, L
    REAL(fp),       INTENT(IN)  :: CL          ! Atomic Cl [molec/cm3]
    REAL(fp),       INTENT(IN)  :: OH          ! OH  [molec/cm3]
    REAL(fp),       INTENT(IN)  :: HO2         ! HO2 [molec/cm3]
    REAL(fp),       INTENT(IN)  :: NO2         ! NO2 [molec/cm3]
    REAL(fp),       INTENT(IN)  :: CLO         ! ClO [molec/cm3]
    REAL(fp),       INTENT(IN)  :: BRO         ! BrO [molec/cm3]
    REAL(fp),       INTENT(IN)  :: BR          ! Br  [molec/cm3]
    TYPE(MetState), INTENT(IN)  :: State_Met   ! Meteorology State obj
!
! !RETURN VALUE:
!
    ! K_Cl : Effective 1st order loss rate of Hg(0) [1/s]
    ! sum of all pathways
    ! HMH 3/8/16 ADDING PROD/LOSS HGCL
    REAL(fp)                    :: K_CL(3)
!
! !REMARKS:
!  ============================================================================
!  This subroutine calculates the net rate of Hg(0) oxidation to Hg(II) through
!     the following reactions. All are gas phase.
!                                                                             .
!     (1  )  Hg(0) + Cl -> HgCl
!     (2  )  HgCl + M -> Hg(0) + Cl
!     (2a )  HgCl + Cl-> Hg(0) + Cl2
!     (3BR)   HgCl + Cl  -> HgCl2
!     (3OH)   HgCl + OH  -> HgClOH
!     (3HO2)  HgCl + HO2 -> HgClHO2
!     (3NO2)  HgCl + NO2 -> HgClNO2
!     (3CLO)  HgCl + ClO -> HgClClO
!     (3BRO)  HgCl + BrO -> HgClBrO
!     (3BR)   HgCl + Br  -> HgClBr
!                                                                             .
!  References:
!  ============================================================================
!  1. Donohoue et al., 2005
!  2. Holmes et al., 2009
!  3. Hynes et al., 2009
!  4. Wilcox, 2009
!  5. Balabanov et al., 2005
!
!  !REVISION HISTORY:
!   21 Nov 2014 - H. Horowitz (HMH) - first version, copied from get_hgbry_rat
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    REAL(fp)               :: k1cl,   k2cl, k2acl
    REAL(fp)               :: k3cl, k3no2, k3ho2
    REAL(fp)               :: k0NO2, k0HO2, kinfNO2, kinfHO2
    REAL(fp)               :: Nair, T, P
!
! !DEFINED PARAMETERS:
!
    !molar gas constant
    REAL(fp), PARAMETER    ::  R = 287e0_fp

    !=================================================================
    ! GET_HGCL_RATE begins here!
    !=================================================================

    T = State_Met%T(I,J,L)        ! Temperature [K]
    P = State_Met%PMID_DRY(I,J,L) ! Pressure [hPa]

    !number density of air [molec/cm3]
    Nair = P * 1e2_fp / (0.029e0_fp * R * T ) * 6.02e23_fp / 1e6_fp

    !k1 from Donohoue et al. 2005
    k1cl=2.2e-32_fp*exp(6.80e2_fp*(1.0e0_fp/T-1.0e0_fp/298e0_fp))*Nair

    !k2 assumed 0 from Holmes et al., 2009; Hynes et al., 2009
    k2cl  = 0e0_fp

    !k2a from Wilcox, 2009
    k2acl  = 1.2e-11_fp * exp(-5.939e3_fp/T)
    ! OK SO:
    ! cm3/(mol*s) = 7.23E12 * exp(-11.8/(RT))
    ! 11.8 is activation energy in kcal/mol.
    !     R in kcal mol-1 K-1 = 1.987d-3
    !     exponent becomes = - 5939 (K-1) /T (K)
    ! 7.23E12 cm3/(mol*s) * 1 mol/6.022E23 molec --> 1.2d-11 cm3 / molec *s

    !k3
    !HMH changed rate to Balabanov et al. 2005 rate
    ! which is appropriate at all pressures, not just high-P
    ! and uses a more accurate level of theory than Goodsite 2004
    ! assume HgBr+Y same as HgCl+Y based on Dibble et al 2012
    k3cl = 3.0e-11_fp

    k0NO2 = 6.91e-29_fp*(300/T)**4.5e0_fp
    k0HO2 = 2.21e-29_fp*(300/T)**4.37e0_fp
    kinfNO2 = 2.21e-11_fp*EXP(4.97e2_fp/T)
    kinfHO2 = 8.42e-12_fp*EXP(6.25e2_fp/T)

    k3no2 = k0NO2*Nair/(1+k0NO2*Nair/kinfNO2) * 0.6e0_fp ** &
            ((1+ (log10(k0NO2*Nair/kinfNO2)) ** (2e0_fp) ) ** (-1e0_fp))

    ! this is for ho2, clo, cl, oh, bro
    k3ho2 = k0HO2*Nair/(1+k0HO2*Nair/kinfHO2) * 0.6e0_fp ** &
            ((1+ (log10(k0HO2*Nair/kinfHO2)) ** (2e0_fp) ) ** (-1e0_fp))

    !!! reverse reaction HgBr + NO2 - > Hg0 + BrNO2. pressure independent
    !!!k3ano2 = 3.4d-12 * exp(391d0/T)

    IF ( CL > SMALLNUM ) THEN

       ! effective 1st order loss of Hg(0) [ 1/s ]
       ! calculate steady state concentration of HgCl
       ! then calculate the total production of Hg(2) from all
       ! of the seven 'Y' pathways as k3[HgCl][Y]!

       ! HMH adding Br here! 8/10/15
       ! hmh 12/11/16 (began 12/10/16) revising based on dibble
       ! HgCl + (Cl, OH, HO2, NO2, ClO, BrO, and BR):
       K_CL(1) = k1cl * CL * (k3cl* BR + k3no2 * NO2 + k3ho2 * (CL + OH + &
                 HO2 + CLO + BRO ))/ ( k2cl + k2acl * CL + k3cl * Br + &
                 k3no2 * NO2 + k3ho2 * ( HO2 + CLO + BRO + CL + OH) )

       ! MAKE PROD LOSS HGCL REACTION RATES 3/8/16
       ! this is for Hg0 + Cl -> HgCl
       ! = [Hg0]*k1*CL
       K_CL(2) = k1cl * CL
       ! this is for HgCl + Cl -> Hg(0) + Cl2
       ! = [HgCl]*Cl*k2a
       K_CL(3) = k1cl * CL * k2acl * CL / &
               ( k2cl + k2acl * CL + k3cl * Br + &
                 k3no2 * NO2 + k3ho2 * ( HO2 + CLO + BRO + CL + OH) )

    ELSE

       ! Avoid divide by zero in rate
       K_CL(1) = 0e0_fp
       K_CL(2) = 0e0_fp
       K_CL(3) = 0e0_fp

    ENDIF

  END FUNCTION GET_HGCL_RATE
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
! !IROUTINE: get_hgbry_rate
!
! !DESCRIPTION: Function GET\_HGBRY\_RATE computes the effective 1st order
!  conversion rate of Hg(0) to Hg(II) via two-step recombination with Br
!  OH, NO2, HO2, ClO, and BrO.  hmh 6/4/14 based on eds 5/2/13, modified
! for consistency with GET\_HGBR\_RATE.
!\\
!\\
! !INTERFACE:
  FUNCTION GET_HGBRY_RATE(I, J, L, BR, OH, NO2, HO2, CLO, BRO, CL, State_Met) &
       RESULT(K_HGBRY)
!
! !USES:
!
    USE State_Met_Mod,      ONLY : MetState
!
! !INPUT PARAMETERS:
!
    INTEGER,            INTENT(IN) :: I, J, L
    REAL(fp),           INTENT(IN) :: BR          ! Atomic Br [molec/cm3]
    REAL(fp),           INTENT(IN) :: OH          ! OH  [molec/cm3]
    REAL(fp),           INTENT(IN) :: NO2         ! NO2 [molec/cm3]
    REAL(fp),           INTENT(IN) :: HO2         ! HO2 [molec/cm3]
    REAL(fp),           INTENT(IN) :: CLO         ! ClO [molec/cm3]
    REAL(fp),           INTENT(IN) :: BRO         ! BrO [molec/cm3]
    REAL(fp),           INTENT(IN) :: CL          ! Cl  [molec/cm3]
    TYPE(MetState),     INTENT(IN) :: State_Met   ! Meteorology State obj
!
!  !RETURN VALUE:
!
    ! K_HgBrY : Effective 1st order loss rate of Hg(0) [1/s]
    ! each pathway - Br, OH, HO2, NO2, ClO, BrO, Cl in that order.
    ! #s 8 - 10 are for R1, R2, and R2a!
    REAL(fp)                       :: K_HGBRY(11)
!
! !REMARKS:
!  ============================================================================
!  This subroutine calculates the net rate of Hg(0) oxidation to Hg(II) through
!     the following reactions. All are gas phase.
!                                                                             .
!     (1  )  Hg(0) + Br -> HgBr
!     (2  )  HgBr + M -> Hg(0) + Br
!     (2a )  HgBr + Br-> Hg(0) + Br2   ! HMH there was a typo in this rxn!
!     (3BR)   HgBr + Br  -> HgBr2
!     (3OH)   HgBr + OH  -> HgBrOH
!     (3HO2)  HgBr + HO2 -> HgBrHO2
!     (3NO2)  HgBr + NO2 -> HgBrNO2
!     (3CLO)  HgBr + ClO -> HgBrClO
!     (3BRO)  HgBr + BrO -> HgBrBrO
!     (3CL)   HgBr + Cl  -> HgBrCl
!                                                                             .
!  References:
!  ============================================================================
!  1. Dibble et al. 2012
!  2. Horowitz et al. in prep
!
!  !REVISION HISTORY:
!  03 May 2013 - E. Corbitt - initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    REAL(fp)             :: k1, k2, k2a, k3ano2
    REAL(fp)             :: k3br, k3no2, k3ho2
    REAL(fp)             :: Nair, N0
    REAL(fp)             :: k0NO2, k0HO2, kinfNO2, kinfHO2
    REAL(fp)             :: C_HgBr, T, P

    REAL(fp), PARAMETER  :: R = 287e+0_fp

    T = State_Met%T(I,J,L)        ! Temperature [K]
    P = State_Met%PMID_DRY(I,J,L) ! Pressure [hPa]

    !number density of air [molec/cm3]
    Nair = P * 1e+2_fp / (0.029e+0_fp * R * T ) * 6.02e+23_fp/1e+6_fp

    !standard air density at STP: 1 atm, 273K [molec/cm3]
    N0 = 1013e+2_fp/(0.029e+0_fp*R*273e+0_fp)*6.02e+23_fp/1e+6_fp

    !k1 from Donohoue et al. 2006
    k1  = 1.46e-32_fp * ( T / 298e+0_fp )**(-1.86e+0_fp) * Nair
    ! the uncertainty here is 1.46+/-0.34 and in the temperature part is
    ! -(1.86+/-1.49)(HMH, from Dibble et al., 2012)

    !k2 from Dibble et al. 2012
    k2  = 1.6e-9_fp * (T/298e+0_fp ) ** (-1.86e+0_fp) &
          * exp(-7801e+0_fp/T) * Nair
    ! uncertainty: (1.60+/-0.37)d-9 * (T/298)^(-(1.86+/-1.49)) &
    !              *exp((-7801+\/-201)/T) (HMH, Dibble et al 2012)

    !k2a from Balabanov et al. 2005
    k2a  = 3.9e-11_fp

    !k3
    !HMH changed rate to Balabanov et al. 2005 rate
    ! which is appropriate at all pressures, not just high-P
    ! and uses a more accurate level of theory than Goodsite 2004
    k3br = 3.0e-11_fp

    k0NO2 = 6.91e-29_fp*(300/T)**4.5e0_fp
    k0HO2 = 2.21e-29_fp*(300/T)**4.37e0_fp
    kinfNO2 = 2.21e-11_fp*EXP(4.97e2_fp/T)
    kinfHO2 = 8.42e-12_fp*EXP(6.25e2_fp/T)

    k3no2 = k0NO2*Nair/(1+k0NO2*Nair/kinfNO2) * 0.6e+0_fp ** &
            ((1+ (log10(k0NO2*Nair/kinfNO2))**(2e+0_fp))**(-1e+0_fp))

    ! this is for ho2, clo, cl, oh, bro
    k3ho2 = k0HO2*Nair/(1+k0HO2*Nair/kinfHO2) *0.6e+0_fp** &
            ((1+ (log10(k0HO2*Nair/kinfHO2))**(2e+0_fp))**(-1e+0_fp))

    ! reverse reaction HgBr + NO2 - > Hg0 + BrNO2. pressure independen
    k3ano2 = 3.4e-12_fp * exp(391e+0_fp/T)

    IF ( BR > SMALLNUM ) THEN

       ! effective 1st order loss of Hg(0) [ 1/s ]
       ! calculate steady state concentration of HgBr
       ! then calculate the production of Hg(2) from each
       ! of the seven 'Y' pathways as k3[HgBr][Y]!

       C_HgBr = k1 * BR   / &
                ( k2 + (k2a +k3br) * BR + (k3NO2+k3ano2) * NO2 &
                + k3HO2 * (HO2 + CLO + BRO + BR + OH + CL) )

       ! HgBr + Br
       K_HGBRY(1) = C_HgBr * k3br * BR

       ! HgBr + OH
       K_HGBRY(2) = C_HgBR *  k3HO2 * OH

       ! HgBr + HO2
       K_HGBRY(3) = C_HgBR *  k3HO2 * HO2

       ! HgBr + NO2
       K_HGBRY(4) = C_HgBR *  k3NO2 * NO2

       ! HgBr + ClO
       K_HGBRY(5) = C_HgBR *  k3HO2 * CLO

       ! HgBr + BrO
       K_HGBRY(6) = C_HgBR *  k3HO2 * BRO

       ! HgBr + Cl
       K_HGBRY(7) = C_HgBR *  k3HO2 * CL

       ! this is for Hg0 + Br -> HgBr
       ! = [Hg0]*k1*BR
       K_HGBRY(8) = k1 * BR

       ! this is for HgBr + M -> Hg(0) + Br
       ! = [HgBr]*k2
       K_HGBRY(9) = C_HgBr * k2

       ! = [HgBr]*Br*k2a
       K_HGBRY(10) = C_HgBR * k2a * BR

       ! HgBr + NO2 -> Hg(0) + BrNO2
       K_HGBRY(11) = C_HgBR * k3ano2 * NO2

    ELSE

       ! Avoid divide by zero in rate
       K_HGBRY(1) = 0
       K_HGBRY(2) = 0
       K_HGBRY(3) = 0
       K_HGBRY(4) = 0
       K_HGBRY(5) = 0
       K_HGBRY(6) = 0
       K_HGBRY(7) = 0
       K_HGBRY(8) = 0
       K_HGBRY(9) = 0
       K_HGBRY(10) = 0
       K_HGBRY(11) = 0

    ENDIF

  END FUNCTION GET_HGBRY_RATE
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: emissmercury
!
! !DESCRIPTION: Subroutine EMISSMERCURY is the driver routine for mercury
!  emissions.
!\\
!\\
! NOTE/TODO: The mercury simulation is the only GEOS-Chem emission code that
! is not yet fully compatible with HEMCO. So far, only the anthropogenic
! emissions are included in HEMCO. For all other sources, the original
! mercury code is used.
!\\
!\\
! For the non-local PBL mixing, all emissions are written into module array
! HG\_EMIS (in kg m-2 s-1). These values are then used by vdiff\_mod.F90.
! This is just a workaround to ensure backwards compatibility of the mercury
! code. Once we have added all mercury emissions to HEMCO, HG\_EMIS is not
! used any more (ckeller, 10/21/2014).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE EMISSMERCURY( Input_Opt,  State_Chm, State_Diag, &
                           State_Grid, State_Met, RC )
!
! !USES:
!
    USE DEPO_MERCURY_MOD,   ONLY : RESET_HG_DEP_ARRAYS
    USE ErrCode_Mod
    USE ERROR_MOD
    USE Input_Opt_Mod,      ONLY : OptInput
    USE LAND_MERCURY_MOD,   ONLY : LAND_MERCURY_FLUX, VEGEMIS
    USE LAND_MERCURY_MOD,   ONLY : SOILEMIS, BIOMASSHG
    USE LAND_MERCURY_MOD,   ONLY : SNOWPACK_MERCURY_FLUX
    USE OCEAN_MERCURY_MOD,  ONLY : OCEAN_MERCURY_FLUX
    USE State_Chm_Mod,      ONLY : ChmState
    USE State_Diag_Mod,     ONLY : DgnState
    USE State_Grid_Mod,     ONLY : GrdState
    USE State_Met_Mod,      ONLY : MetState
    USE TIME_MOD,           ONLY : GET_MONTH, ITS_A_NEW_MONTH
    USE UnitConv_Mod,       ONLY : Convert_Spc_Units
    
    ! Added for GTMM (ccc, 11/19/09)
    !USE LAND_MERCURY_MOD,   ONLY : GTMM_DR
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
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!
! !REMARKS:
!
!
! !REVISION HISTORY:
!  03 Jun 2013 - N. (Eckley) Selin - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    LOGICAL, SAVE     :: FIRST = .TRUE.
    INTEGER           :: THISMONTH, I, J
    CHARACTER(LEN=63) :: OrigUnit

    ! For fields from Input_Opt
    LOGICAL           :: LDYNOCEAN
    LOGICAL           :: LGTMM
    LOGICAL           :: LNLPBL
    LOGICAL           :: LPREINDHG
    LOGICAL           :: LEMIS
    LOGICAL           :: prtDebug

    ! Strings
    CHARACTER(LEN=255)   :: ErrMsg
    CHARACTER(LEN=255)   :: ThisLoc

    !=================================================================
    ! EMISSMERCURY begins here!
    !=================================================================

    ! Assume success
    RC      = GC_SUCCESS
    ErrMsg  = ''
    ThisLoc = ' -> at EMISSMERCURY (in module GeosCore/mercury_mod.F90)'

    ! Copy fields from Input_Opt
    LDYNOCEAN = Input_Opt%LDYNOCEAN
    LGTMM     = Input_Opt%LGTMM
    LNLPBL    = Input_Opt%LNLPBL
    LPREINDHG = Input_Opt%LPREINDHG
    LEMIS     = Input_Opt%LEMIS
    prtDebug  = ( Input_Opt%LPRT .and. Input_Opt%amIRoot )

    ! Convert species units to [kg] for EMISSMERCURY (ewl, 8/12/15)
    CALL Convert_Spc_Units( Input_Opt, State_Chm, State_Grid, State_Met, &
                            'kg', RC, OrigUnit=OrigUnit )

    ! Trap potential errors
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered in "Convert_Spc_Units" #1!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    ! First-time initialization
    IF ( FIRST ) THEN

       ! Check that emissions are turned on. Print error message
       ! and stop GEOS-Chem if emissions are turned off (ewl, 9/1/15)
       IF ( .not. LEMIS ) THEN
          ErrMsg = 'ERROR: Hg emissions are need for simulation ' // &
                   'but emissions are turned off in input.geos'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       ! Read anthro, ocean, land emissions of Hg from disk
       CALL MERCURY_READYR( Input_Opt, RC )

       ! Trap potential errors
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Error encountered call to "MERCURY_READYR"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       ! Reset first-time flag
       FIRST = .FALSE.
    ENDIF

    !=================================================================
    ! Call emission routines for Hg(0), Hg(II), and Hg(P)
    !=================================================================

    ! Ocean flux of Hg(0)
    IF ( LDYNOCEAN ) THEN

       ! Set to zero to clear emissions from previous time step
       ! (cdh, 4/30/09)
       EHg0_oc = 0e+0_fp

       CALL OCEAN_MERCURY_FLUX( Input_Opt,  State_Chm, State_Diag, &
                                State_Grid, State_Met, EHg0_oc,    RC )

       ! Trap potential errors
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Error encountered call to "OCEAN_MERCURY_FLUX"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       IF ( prtDebug ) CALL DEBUG_MSG( '### EMISSMERCURY: a OCEAN_FLUX' )

    ELSE
       EHg0_oc = 0e+0_fp
       CALL OFFLINEOCEAN_READMO( State_Chm, State_Diag, State_Grid, &
                                 State_Met, EHg0_oc, RC )
    ENDIF

    !==========================================================================
    ! Disable GTMM until it is brought up-to-date (mps, 3/10/19)
    !IF ( LGTMM ) THEN
    !   !--------------------------------------------------------------
    !   ! Here we are using the Global Terrstrial Mercury Model...
    !   !--------------------------------------------------------------
    !   IF ( ITS_A_NEW_MONTH() ) THEN
    !
    !      ! OLD CODE: CALL GTMM_DR( EHg0_gtm(:,:), State_Met )
    !      CALL GTMM_DR( Input_Opt, State_Chm, State_Grid, State_Met, &
    !                    EHg0_gtm,  RC )
    !
    !      IF ( prtDebug ) CALL DEBUG_MSG( '### EMISSMERCURY: a GTMM' )
    !   ENDIF
    !
    !ELSE
    !==========================================================================

    !--------------------------------------------------------------
    ! Here we are NOT using the Global Terrstrial Mercury Model...
    !--------------------------------------------------------------
    CALL LAND_MERCURY_FLUX ( EHg0_ln, LHGSNOW, State_Grid, State_Met )
    IF ( prtDebug )  CALL DEBUG_MSG( '### EMISSMERCURY: a LAND_FLUX' )

    CALL VEGEMIS( Input_Opt, State_Met, LVEGEMIS,  EHg0_dist, EHg0_vg,   RC )
    IF ( prtDebug ) CALL DEBUG_MSG( '### EMISSMERCURY: a VEGEMIS' )

    CALL SOILEMIS( EHg0_dist, EHg0_so, State_Grid, State_Met )
    IF ( prtDebug ) CALL DEBUG_MSG( '### EMISSMERCURY: a SOILEMIS' )

    !==========================================================================
    !ENDIF
    !==========================================================================

    CALL SNOWPACK_MERCURY_FLUX( EHg0_snow,  LHGSNOW,  State_Chm, &
                                State_Grid, State_Met )
    IF ( prtDebug ) CALL DEBUG_MSG( '### EMISSMERCURY: a SNOW_FLUX' )

    CALL BIOMASSHG( Input_Opt, EHg0_bb, RC )
    IF ( prtDebug ) CALL DEBUG_MSG( '### EMISSMERCURY: a BIOMASS' )

    IF ( LnoUSAemis ) THEN

       ! No Anthropogenic emissions over USA (cdh, 5/6/2009)
       DO J = 1, State_Grid%NY
       DO I = 1, State_Grid%NX

          IF ( State_Grid%XMid(I,J) > -140 .AND. &
               State_Grid%XMid(I,J) < -40  .AND. &
               State_Grid%YMid(I,J) >  10  .AND. &
               State_Grid%YMid(I,J) <  70 ) THEN

             EHg0_An(I,J) = 0e+0_fp
             EHg2_An(I,J) = 0e+0_fp
             EHgP_An(I,J) = 0e+0_fp
             EHg0_bb(I,J) = 0e+0_fp
             EHg0_am(I,J) = 0e+0_fp

          ENDIF

       ENDDO
       ENDDO

    ENDIF

    CALL RESET_HG_DEP_ARRAYS
    IF ( prtDebug ) CALL DEBUG_MSG( '### EMISSMERCURY: ' // &
         'a RESET_HG_DEP_ARRAYS' )

    ! If we are using the non-local PBL mixing,
    ! we need to initialize the EMIS_SAVE array (cdh, 08/27/09)
    ! EMIS_SAVE is now HG_EMIS (ckeller, 10/21/2014)
    IF ( LNLPBL ) HG_EMIS = 0.0e+0_fp

    ! Add Hg(0) source into State_Chm%Species [kg]
    CALL SRCHg0( Input_Opt,  State_Chm, State_Diag, State_Grid, State_Met, RC )
    IF ( prtDebug ) CALL DEBUG_MSG( '### EMISSMERCURY: a SRCHg0' )

    ! Add Hg(II) source into State_Chm%Species [kg]
    CALL SRCHg2( Input_Opt,  State_Chm, State_Diag, State_Grid, State_Met, RC )
    IF ( prtDebug ) CALL DEBUG_MSG( '### EMISSMERCURY: a SRCHg2' )

    ! Add HgP source into State_Chm%Species [kg]
    CALL SRCHgP( Input_Opt, State_Chm, State_Grid, State_Met, RC )
    IF ( prtDebug ) CALL DEBUG_MSG( '### EMISSMERCURY: a SRCHgP' )

    ! Convert species units back to original unit
    CALL Convert_Spc_Units( Input_Opt, State_Chm, State_Grid, State_Met, &
                            OrigUnit,  RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error('Unit conversion error', RC, &
                     'Routine EMISSMERCURY in mercury_mod.F90')
       RETURN
    ENDIF

  END SUBROUTINE EMISSMERCURY
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: emithg
!
! !DESCRIPTION: Subroutine EMITHG directs emission either to the chemical
!  species array (State\_Chm%Species) directly or to Hg\_EMIS for use by the
!  non-local PBL mixing. This is a programming convenience.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE EMITHG( I, J, L, ID, E_HG, Input_Opt, State_Chm, State_Grid )
!
! !USES:
!
    USE ErrCode_Mod
    USE Input_Opt_Mod,      ONLY : OptInput
    USE State_Chm_Mod,      ONLY : ChmState
    USE State_Grid_Mod,     ONLY : GrdState
    USE TIME_MOD,           ONLY : GET_TS_EMIS
!
! !INPUT PARAMETERS:
!
    INTEGER,        INTENT(IN)    :: I, J, L, ID ! Grid boxes + species #
    REAL(fp),       INTENT(IN)    :: E_Hg        ! Hg emissions
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
    TYPE(GrdState), INTENT(IN)    :: State_Grid  ! Grid State object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
!
! !REVISION HISTORY:
!  27 Aug 2009 - C. Holmes   - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    REAL(fp)          :: AM2, TS

    ! Pointers
    REAL(fp), POINTER :: Spc(:,:,:,:)

    !=================================================================
    ! EMITHG begins here!
    !=================================================================

    IF ( Input_Opt%LNLPBL ) THEN

       !--------------------------------------------------------------
       ! We are using FULL PBL MIXING (routine TURBDAY)
       !
       ! Save emissions for non-local PBL mixing or emit directly.
       ! Make sure that emitted mass is non-negative
       ! This is hear only for consistency with old code which warned
       ! of underflow error (cdh, 08/27/09)
       ! EMIS_SAVE is now HG_EMIS array. Convert kg to kg/m2/s
       ! (ckeller, 09/23/2014)
       !--------------------------------------------------------------

       ! Surface area [m2]
       AM2             = State_Grid%Area_M2(I,J)

       ! Emission timestep
       TS              = GET_TS_EMIS()

       ! Save emissions as [kg/m2/s].  These will be added
       ! to the chemical species array in routine DO_TEND
       ! (in mixing_mod.F90).
       HG_EMIS(I,J,ID) = HG_EMIS(I,J,ID) + ( MAX(E_HG,0e+0_fp) / AM2 / TS )

    ELSE

       !--------------------------------------------------------------
       ! We are using FULL PBL MIXING (routine TURBDAY)
       ! so add directly to the State_Chm%Species array
       !--------------------------------------------------------------

       ! Point to the chemical spcies array [kg]
       Spc             => State_Chm%Species

       ! Add emissions
       Spc(I,J,L,ID)   =  Spc(I,J,L,ID) + MAX( E_HG, 0e+0_fp )

       ! Free pointer
       Spc             => NULL()

    ENDIF

  END SUBROUTINE EMITHG
!EOP
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: srcHg0
!
! !DESCRIPTION: Subroutine SRCHg0 is the subroutine for Hg(0) emissions.
!  Emissions of Hg(0) will be distributed throughout the boundary layer.
!  (eck, cdh, bmy, 1/21/05, 4/6/06)
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE SRCHg0( Input_Opt,  State_Chm, State_Diag, &
                     State_Grid, State_Met, RC )
!
! !USES:
!
#ifdef BPCH_DIAG
    USE DIAG03_MOD,         ONLY : AD03, ND03, AD03_nat
#endif
    USE ErrCode_Mod
    USE Input_Opt_Mod,      ONLY : OptInput
    USE State_Chm_Mod,      ONLY : ChmState
    USE State_Diag_Mod,     ONLY : DgnState
    USE State_Grid_Mod,     ONLY : GrdState
    USE State_Met_Mod,      ONLY : MetState
    USE TIME_MOD,           ONLY : GET_TS_EMIS
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
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!
! !REVISION HISTORY:
!  21 Jan 2005 - N. (Eckley) Selin, C. Holmes - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER  :: I,      J,    L,        N,    NN,     PBL_MAX
    REAL(fp) :: DTSRCE, E_Hg, F_OF_PBL, T_Hg, T_Hg_An

    ! For fields from Input_Opt
    LOGICAL  :: LSPLIT, LPREINDHG, LGTMM

    !=================================================================
    ! SRCHg0 begins here!
    !=================================================================

    ! Assume success
    RC        = GC_SUCCESS

    ! Copy values from Input_Opt
    LSPLIT    = Input_Opt%LSPLIT
    LPREINDHG = Input_Opt%LPREINDHG
    LGTMM     = Input_Opt%LGTMM

    ! Emission timestep [s]
    DTSRCE    = GET_TS_EMIS()

    ! Maximum extent of the PBL [model levels]
    PBL_MAX   = State_Met%PBL_MAX_L

    ! Loop over grid boxes
    !$OMP PARALLEL DO       &
    !$OMP DEFAULT( SHARED ) &
    !$OMP PRIVATE( I, J, L, N, T_Hg_An, T_Hg, F_OF_PBL, E_Hg, NN)
    DO J = 1, State_Grid%NY
    DO I = 1, State_Grid%NX

       IF ( LPREINDHG ) THEN  !eds

          ! Anthropogenic emissions
          T_Hg_An = 0e+0_fp

          ! No biomass burning emissions
          EHg0_bb(I,J) = 0e+0_fp

       ELSE

          ! Compute total anthropogenic Hg(0) emissions
          T_Hg_An = EHg0_an(I,J)

          IF ( LAnthroHgOnly ) THEN
             ! No other emissions
             EHg0_bb(I,J)      = 0e+0_fp
             EHg0_oc(I,J,:)    = 0e+0_fp
             EHg0_snow(I,J,:)  = 0e+0_fp
             IF ( LGTMM ) THEN
                EHg0_gtm(I,J) = 0e+0_fp
             ELSE
                EHg0_ln(I,J,:) = 0e+0_fp
                EHg0_vg(I,J)   = 0e+0_fp
                EHg0_so(I,J)   = 0e+0_fp
             ENDIF
          ENDIF

       ENDIF

       ! Compute total Hg(0) emissions (anthro+oceans+land+natural)
       IF ( LGTMM ) THEN
          T_Hg = T_Hg_An +                &
                 EHg0_bb(I,J) +           &
                 EHg0_oc(I,J,ID_Hg_tot) + &
                 EHg0_geo(I,J) +          &
                 EHg0_gtm(I,J) +          &
                 EHg0_snow(I,J,ID_Hg_tot)
       ELSE
          T_Hg = T_Hg_An +                &
                 EHg0_bb(I,J) +           &
                 EHg0_oc(I,J,ID_Hg_tot) + &
                 EHg0_ln(I,J,ID_Hg_tot) + &
                 EHg0_geo(I,J) +          &
                 EHg0_vg(I,J) +           &
                 EHg0_so(I,J) +           &
                 EHg0_snow(I,J,ID_Hg_tot)
       ENDIF

       !==============================================================
       ! Partition Hg0 throughout PBL; store into State_Chm%Species [kg]
       ! Now make sure State_Chm%Species does not underflow (cdh, bmy, 4/6/06)
       !==============================================================

       ! Loop up to max PBL level
       DO L = 1, PBL_MAX

          ! Fraction of box (I,J,L) w/in the PBL [unitless]
          F_OF_PBL = State_Met%F_OF_PBL(I,J,L)

          !-----------------
          ! Total Hg species
          !-----------------
          N    = Hg0_Id_List(ID_Hg_tot)
          E_Hg = F_OF_PBL * T_Hg * DTSRCE
          CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

          !-----------------
          ! Tagged species
          !-----------------
          IF ( LSPLIT ) THEN

             !--------------------
             ! Primary emissions
             !--------------------

             ! Anthro Canada Hg0
             N    = Hg0_Id_List(ID_Hg_can)
             E_Hg = F_OF_PBL * EHg0_can(I,J) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Anthro United States Hg0
             N    = Hg0_Id_List(ID_Hg_usa)
             E_Hg = F_OF_PBL * EHg0_usa(I,J) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Anthro Central America Hg0
             N    = Hg0_Id_List(ID_Hg_cam)
             E_Hg = F_OF_PBL * EHg0_cam(I,J) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Anthro South America Hg0
             N    = Hg0_Id_List(ID_Hg_sam)
             E_Hg = F_OF_PBL * EHg0_sam(I,J) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Anthro West Africa Hg0
             N    = Hg0_Id_List(ID_Hg_waf)
             E_Hg = F_OF_PBL * EHg0_waf(I,J) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Anthro East Africa Hg0
             N    = Hg0_Id_List(ID_Hg_eaf)
             E_Hg = F_OF_PBL * EHg0_eaf(I,J) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Anthro South Africa Hg0
             N    = Hg0_Id_List(ID_Hg_saf)
             E_Hg = F_OF_PBL * EHg0_saf(I,J) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Anthro North Africa Hg0
             N    = Hg0_Id_List(ID_Hg_naf)
             E_Hg = F_OF_PBL * EHg0_naf(I,J) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Anthro OECD Europe Hg0
             N    = Hg0_Id_List(ID_Hg_eur)
             E_Hg = F_OF_PBL * EHg0_eur(I,J) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Anthro Eastern Europe Hg0
             N    = Hg0_Id_List(ID_Hg_eeu)
             E_Hg = F_OF_PBL * EHg0_eeu(I,J) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Anthro Middle East Hg0
             N    = Hg0_Id_List(ID_Hg_mde)
             E_Hg = F_OF_PBL * EHg0_mde(I,J) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Anthro Former USSR Hg0
             N    = Hg0_Id_List(ID_Hg_sov)
             E_Hg = F_OF_PBL * EHg0_sov(I,J) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Anthro South Asia Hg0
             N    = Hg0_Id_List(ID_Hg_sas)
             E_Hg = F_OF_PBL * EHg0_sas(I,J) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Anthro East Asia Hg0
             N    = Hg0_Id_List(ID_Hg_eas)
             E_Hg = F_OF_PBL * EHg0_eas(I,J) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Anthro Southeast Asia Hg0
             N    = Hg0_Id_List(ID_Hg_sea)
             E_Hg = F_OF_PBL * EHg0_sea(I,J) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Anthro Japan Hg0
             N    = Hg0_Id_List(ID_Hg_jpn)
             E_Hg = F_OF_PBL * EHg0_jpn(I,J) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Anthro Oceania Hg0
             N    = Hg0_Id_List(ID_Hg_oce)
             E_Hg = F_OF_PBL * EHg0_oce(I,J) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Soil emissions
             N    = Hg0_Id_List(ID_Hg_so)
             E_Hg = F_OF_PBL * EHg0_so(I,J) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Biomass burning emissions
             N    = Hg0_Id_List(ID_Hg_bb)
             E_Hg = F_OF_PBL * EHg0_bb(I,J) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Geogenic emissions
             N    = Hg0_Id_List(ID_Hg_geo)
             E_Hg = F_OF_PBL * EHg0_geo(I,J) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             !--------------------
             ! Ocean re-emissions
             !--------------------

             ! Anthro Canada re-emission from ocean
             N    = Hg0_Id_List(ID_Hg_can)
             E_Hg = F_OF_PBL * EHg0_oc(I,J,ID_Hg_can) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Anthro United States re-emission from ocean
             N    = Hg0_Id_List(ID_Hg_usa)
             E_Hg = F_OF_PBL * EHg0_oc(I,J,ID_Hg_usa) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Anthro Central America re-emission from ocean
             N    = Hg0_Id_List(ID_Hg_cam)
             E_Hg = F_OF_PBL * EHg0_oc(I,J,ID_Hg_cam) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Anthro South America re-emission from ocean
             N    = Hg0_Id_List(ID_Hg_sam)
             E_Hg = F_OF_PBL * EHg0_oc(I,J,ID_Hg_sam) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Anthro West Africa re-emission from ocean
             N    = Hg0_Id_List(ID_Hg_waf)
             E_Hg = F_OF_PBL * EHg0_oc(I,J,ID_Hg_waf) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Anthro East Africa re-emission from ocean
             N    = Hg0_Id_List(ID_Hg_eaf)
             E_Hg = F_OF_PBL * EHg0_oc(I,J,ID_Hg_eaf) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Anthro South Africa re-emission from ocean
             N    = Hg0_Id_List(ID_Hg_saf)
             E_Hg = F_OF_PBL * EHg0_oc(I,J,ID_Hg_saf) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Anthro North Africa re-emission from ocean
             N    = Hg0_Id_List(ID_Hg_naf)
             E_Hg = F_OF_PBL * EHg0_oc(I,J,ID_Hg_naf) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Anthro OECD Europe re-emission from ocean
             N    = Hg0_Id_List(ID_Hg_eur)
             E_Hg = F_OF_PBL * EHg0_oc(I,J,ID_Hg_eur) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Anthro Eastern Europe re-emission from ocean
             N    = Hg0_Id_List(ID_Hg_eeu)
             E_Hg = F_OF_PBL * EHg0_oc(I,J,ID_Hg_eeu) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Anthro Middle East re-emission from ocean
             N    = Hg0_Id_List(ID_Hg_mde)
             E_Hg = F_OF_PBL * EHg0_oc(I,J,ID_Hg_mde) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Anthro Former Soviet Union re-emission from ocean
             N    = Hg0_Id_List(ID_Hg_sov)
             E_Hg = F_OF_PBL * EHg0_oc(I,J,ID_Hg_sov) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Anthro South Asia re-emission from ocean
             N    = Hg0_Id_List(ID_Hg_sas)
             E_Hg = F_OF_PBL * EHg0_oc(I,J,ID_Hg_sas) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Anthro East Asia re-emission from ocean
             N    = Hg0_Id_List(ID_Hg_eas)
             E_Hg = F_OF_PBL * EHg0_oc(I,J,ID_Hg_eas) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Anthro Southeast Asia re-emission from ocean
             N    = Hg0_Id_List(ID_Hg_sea)
             E_Hg = F_OF_PBL * EHg0_oc(I,J,ID_Hg_sea) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Anthro Japan re-emission from ocean
             N    = Hg0_Id_List(ID_Hg_jpn)
             E_Hg = F_OF_PBL * EHg0_oc(I,J,ID_Hg_jpn) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Anthro Oceania re-emission from ocean
             N    = Hg0_Id_List(ID_Hg_oce)
             E_Hg = F_OF_PBL * EHg0_oc(I,J,ID_Hg_oce) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Soil re-emission from ocean
             N    = Hg0_Id_List(ID_Hg_so)
             E_Hg = F_OF_PBL * EHg0_oc(I,J,ID_Hg_so) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Biomass burning re-emission from ocean
             N    = Hg0_Id_List(ID_Hg_bb)
             E_Hg = F_OF_PBL * EHg0_oc(I,J,ID_Hg_bb) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Geogenic re-emission from ocean
             N    = Hg0_Id_List(ID_Hg_geo)
             E_Hg = F_OF_PBL * EHg0_oc(I,J,ID_Hg_geo) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Below Mixed Layer Mid-Atlantic re-emission from ocean
             N    = Hg0_Id_List(ID_Hg_atl)
             E_Hg = F_OF_PBL * EHg0_oc(I,J,ID_Hg_atl) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Below Mixed Layer North Atlantic re-emission from ocean
             N    = Hg0_Id_List(ID_Hg_nat)
             E_Hg = F_OF_PBL * EHg0_oc(I,J,ID_Hg_nat) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Below Mixed Layer South Atlantic re-emission from ocean
             N    = Hg0_Id_List(ID_Hg_sat)
             E_Hg = F_OF_PBL * EHg0_oc(I,J,ID_Hg_sat) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Below Mixed Layer North Pacific re-emission from ocean
             N    = Hg0_Id_List(ID_Hg_npa)
             E_Hg = F_OF_PBL * EHg0_oc(I,J,ID_Hg_npa) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Below Mixed Layer Arctic re-emission from ocean
             N    = Hg0_Id_List(ID_Hg_arc)
             E_Hg = F_OF_PBL * EHg0_oc(I,J,ID_Hg_arc) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt,  State_Chm, State_Grid )

             ! Below Mixed Layer Antarctic re-emission from ocean
             N    = Hg0_Id_List(ID_Hg_ant)
             E_Hg = F_OF_PBL * EHg0_oc(I,J,ID_Hg_ant) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Below Mixed Layer All Other Ocean re-emission from ocean
             N    = Hg0_Id_List(ID_Hg_ocn)
             E_Hg = F_OF_PBL * EHg0_oc(I,J,ID_Hg_ocn) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Residual Stratosphere from Spin-up re-emission from ocean
             N    = Hg0_Id_List(ID_Hg_str)
             E_Hg = F_OF_PBL * EHg0_oc(I,J,ID_Hg_str) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             !--------------------
             ! Land re-emissions
             !--------------------

             ! Anthro Canada re-emission from land
             N    = Hg0_Id_List(ID_Hg_can)
             E_Hg = F_OF_PBL * EHg0_ln(I,J,ID_Hg_can) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Anthro United States re-emission from land
             N    = Hg0_Id_List(ID_Hg_usa)
             E_Hg = F_OF_PBL * EHg0_ln(I,J,ID_Hg_usa) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Anthro Central America re-emission from land
             N    = Hg0_Id_List(ID_Hg_cam)
             E_Hg = F_OF_PBL * EHg0_ln(I,J,ID_Hg_cam) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Anthro South America re-emission from land
             N    = Hg0_Id_List(ID_Hg_sam)
             E_Hg = F_OF_PBL * EHg0_ln(I,J,ID_Hg_sam) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Anthro West Africa re-emission from land
             N    = Hg0_Id_List(ID_Hg_waf)
             E_Hg = F_OF_PBL * EHg0_ln(I,J,ID_Hg_waf) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Anthro East Africa re-emission from land
             N    = Hg0_Id_List(ID_Hg_eaf)
             E_Hg = F_OF_PBL * EHg0_ln(I,J,ID_Hg_eaf) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Anthro South Africa re-emission from land
             N    = Hg0_Id_List(ID_Hg_saf)
             E_Hg = F_OF_PBL * EHg0_ln(I,J,ID_Hg_saf) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Anthro North Africa re-emission from land
             N    = Hg0_Id_List(ID_Hg_naf)
             E_Hg = F_OF_PBL * EHg0_ln(I,J,ID_Hg_naf) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Anthro OECD Europe re-emission from land
             N    = Hg0_Id_List(ID_Hg_eur)
             E_Hg = F_OF_PBL * EHg0_ln(I,J,ID_Hg_eur) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Anthro Eastern Europe re-emission from land
             N    = Hg0_Id_List(ID_Hg_eeu)
             E_Hg = F_OF_PBL * EHg0_ln(I,J,ID_Hg_eeu) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Anthro Middle East re-emission from land
             N    = Hg0_Id_List(ID_Hg_mde)
             E_Hg = F_OF_PBL * EHg0_ln(I,J,ID_Hg_mde) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Anthro Former Soviet Union re-emission from land
             N    = Hg0_Id_List(ID_Hg_sov)
             E_Hg = F_OF_PBL * EHg0_ln(I,J,ID_Hg_sov) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Anthro South Asia re-emission from land
             N    = Hg0_Id_List(ID_Hg_sas)
             E_Hg = F_OF_PBL * EHg0_ln(I,J,ID_Hg_sas) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Anthro East Asia re-emission from land
             N    = Hg0_Id_List(ID_Hg_eas)
             E_Hg = F_OF_PBL * EHg0_ln(I,J,ID_Hg_eas) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Anthro Southeast Asia re-emission from land
             N    = Hg0_Id_List(ID_Hg_sea)
             E_Hg = F_OF_PBL * EHg0_ln(I,J,ID_Hg_sea) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Anthro Japan re-emission from land
             N    = Hg0_Id_List(ID_Hg_jpn)
             E_Hg = F_OF_PBL * EHg0_ln(I,J,ID_Hg_jpn) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Anthro Oceania re-emission from land
             N    = Hg0_Id_List(ID_Hg_oce)
             E_Hg = F_OF_PBL * EHg0_ln(I,J,ID_Hg_oce) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Soil re-emission from land
             N    = Hg0_Id_List(ID_Hg_so)
             E_Hg = F_OF_PBL * EHg0_ln(I,J,ID_Hg_so) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Biomass burning re-emission from land
             N    = Hg0_Id_List(ID_Hg_bb)
             E_Hg = F_OF_PBL * EHg0_ln(I,J,ID_Hg_bb) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Geogenic re-emission from land
             N    = Hg0_Id_List(ID_Hg_geo)
             E_Hg = F_OF_PBL * EHg0_ln(I,J,ID_Hg_geo) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Below Mixed Layer Mid-Atlantic re-emission from land
             N    = Hg0_Id_List(ID_Hg_atl)
             E_Hg = F_OF_PBL * EHg0_ln(I,J,ID_Hg_atl) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Below Mixed Layer North Atlantic re-emission from land
             N    = Hg0_Id_List(ID_Hg_nat)
             E_Hg = F_OF_PBL * EHg0_ln(I,J,ID_Hg_nat) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Below Mixed Layer South Atlantic re-emission from land
             N    = Hg0_Id_List(ID_Hg_sat)
             E_Hg = F_OF_PBL * EHg0_ln(I,J,ID_Hg_sat) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Below Mixed Layer North Pacific re-emission from land
             N    = Hg0_Id_List(ID_Hg_npa)
             E_Hg = F_OF_PBL * EHg0_ln(I,J,ID_Hg_npa) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Below Mixed Layer Arctic re-emission from land
             N    = Hg0_Id_List(ID_Hg_arc)
             E_Hg = F_OF_PBL * EHg0_ln(I,J,ID_Hg_arc) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Below Mixed Layer Antarctic re-emission from land
             N    = Hg0_Id_List(ID_Hg_ant)
             E_Hg = F_OF_PBL * EHg0_ln(I,J,ID_Hg_ant) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Below Mixed Layer All Other Ocean re-emission from land
             N    = Hg0_Id_List(ID_Hg_ocn)
             E_Hg = F_OF_PBL * EHg0_ln(I,J,ID_Hg_ocn) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Residual Stratosphere from Spin-up re-emission from land
             N    = Hg0_Id_List(ID_Hg_str)
             E_Hg = F_OF_PBL * EHg0_ln(I,J,ID_Hg_str) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             !--------------------
             ! Snow re-emissions
             !--------------------

             ! Anthro Canada re-emission from snow
             N    = Hg0_Id_List(ID_Hg_can)
             E_Hg = F_OF_PBL * EHg0_snow(I,J,ID_Hg_can) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Anthro United States re-emission from snow
             N    = Hg0_Id_List(ID_Hg_usa)
             E_Hg = F_OF_PBL * EHg0_snow(I,J,ID_Hg_usa) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Anthro Central America re-emission from snow
             N    = Hg0_Id_List(ID_Hg_cam)
             E_Hg = F_OF_PBL * EHg0_snow(I,J,ID_Hg_cam) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Anthro South America re-emission from snow
             N    = Hg0_Id_List(ID_Hg_sam)
             E_Hg = F_OF_PBL * EHg0_snow(I,J,ID_Hg_sam) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Anthro West Africa re-emission from snow
             N    = Hg0_Id_List(ID_Hg_waf)
             E_Hg = F_OF_PBL * EHg0_snow(I,J,ID_Hg_waf) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Anthro East Africa re-emission from snow
             N    = Hg0_Id_List(ID_Hg_eaf)
             E_Hg = F_OF_PBL * EHg0_snow(I,J,ID_Hg_eaf) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Anthro South Africa re-emission from snow
             N    = Hg0_Id_List(ID_Hg_saf)
             E_Hg = F_OF_PBL * EHg0_snow(I,J,ID_Hg_saf) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Anthro North Africa re-emission from snow
             N    = Hg0_Id_List(ID_Hg_naf)
             E_Hg = F_OF_PBL * EHg0_snow(I,J,ID_Hg_naf) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Anthro OECD Europe re-emission from snow
             N    = Hg0_Id_List(ID_Hg_eur)
             E_Hg = F_OF_PBL * EHg0_snow(I,J,ID_Hg_eur) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Anthro Eastern Europe re-emission from snow
             N    = Hg0_Id_List(ID_Hg_eeu)
             E_Hg = F_OF_PBL * EHg0_snow(I,J,ID_Hg_eeu) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Anthro Middle East re-emission from snow
             N    = Hg0_Id_List(ID_Hg_mde)
             E_Hg = F_OF_PBL * EHg0_snow(I,J,ID_Hg_mde) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Anthro Former Soviet Union re-emission from snow
             N    = Hg0_Id_List(ID_Hg_sov)
             E_Hg = F_OF_PBL * EHg0_snow(I,J,ID_Hg_sov) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Anthro South Asia re-emission from snow
             N    = Hg0_Id_List(ID_Hg_sas)
             E_Hg = F_OF_PBL * EHg0_snow(I,J,ID_Hg_sas) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Anthro East Asia re-emission from snow
             N    = Hg0_Id_List(ID_Hg_eas)
             E_Hg = F_OF_PBL * EHg0_snow(I,J,ID_Hg_eas) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Anthro Southeast Asia re-emission from snow
             N    = Hg0_Id_List(ID_Hg_sea)
             E_Hg = F_OF_PBL * EHg0_snow(I,J,ID_Hg_sea) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Anthro Japan re-emission from snow
             N    = Hg0_Id_List(ID_Hg_jpn)
             E_Hg = F_OF_PBL * EHg0_snow(I,J,ID_Hg_jpn) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Anthro Oceania re-emission from snow
             N    = Hg0_Id_List(ID_Hg_oce)
             E_Hg = F_OF_PBL * EHg0_snow(I,J,ID_Hg_oce) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Soil re-emission from snow
             N    = Hg0_Id_List(ID_Hg_so)
             E_Hg = F_OF_PBL * EHg0_snow(I,J,ID_Hg_so) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Biomass burning re-emission from snow
             N    = Hg0_Id_List(ID_Hg_bb)
             E_Hg = F_OF_PBL * EHg0_snow(I,J,ID_Hg_bb) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Geogenic re-emission from snow
             N    = Hg0_Id_List(ID_Hg_geo)
             E_Hg = F_OF_PBL * EHg0_snow(I,J,ID_Hg_geo) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Below Mixed Layer Mid-Atlantic re-emission from snow
             N    = Hg0_Id_List(ID_Hg_atl)
             E_Hg = F_OF_PBL * EHg0_snow(I,J,ID_Hg_atl) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Below Mixed Layer North Atlantic re-emission from snow
             N    = Hg0_Id_List(ID_Hg_nat)
             E_Hg = F_OF_PBL * EHg0_snow(I,J,ID_Hg_nat) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Below Mixed Layer South Atlantic re-emission from snow
             N    = Hg0_Id_List(ID_Hg_sat)
             E_Hg = F_OF_PBL * EHg0_snow(I,J,ID_Hg_sat) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Below Mixed Layer North Pacific re-emission from snow
             N    = Hg0_Id_List(ID_Hg_npa)
             E_Hg = F_OF_PBL * EHg0_snow(I,J,ID_Hg_npa) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Below Mixed Layer Arctic re-emission from snow
             N    = Hg0_Id_List(ID_Hg_arc)
             E_Hg = F_OF_PBL * EHg0_snow(I,J,ID_Hg_arc) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Below Mixed Layer Antarctic re-emission from snow
             N    = Hg0_Id_List(ID_Hg_ant)
             E_Hg = F_OF_PBL * EHg0_snow(I,J,ID_Hg_ant) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Below Mixed Layer All Other Ocean re-emission from snow
             N    = Hg0_Id_List(ID_Hg_ocn)
             E_Hg = F_OF_PBL * EHg0_snow(I,J,ID_Hg_ocn) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Residual Stratosphere from Spin-up re-emission from snow
             N    = Hg0_Id_List(ID_Hg_str)
             E_Hg = F_OF_PBL * EHg0_snow(I,J,ID_Hg_str) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

          ENDIF
       ENDDO

#ifdef BPCH_DIAG
       !==============================================================
       ! %%%%% ND03 (bpch) DIAGNOSTIC %%%%%
       !
       ! Total Hg(0) emissions [kg]
       ! 1=anthro; 3=from ocean; 4=land re-emission; 5=natural src
       !==============================================================

       IF ( ND03 > 0 ) THEN

          IF ( LGTMM ) THEN

             !eds 9/9/10
             DO NN = 1, N_HG_CATS
                AD03(I,J, 1,NN)=AD03(I,J, 1,NN)+(T_Hg_An          *DTSRCE)
                AD03(I,J, 3,NN)=AD03(I,J, 3,NN)+(EHg0_oc(I,J,NN)  *DTSRCE)
                AD03(I,J, 4,NN)=AD03(I,J, 4,NN)+(EHg0_gtm(I,J)    *DTSRCE)
                AD03(I,J, 5,NN)=AD03(I,J, 5,NN)+(EHg0_geo(I,J)    *DTSRCE)
                AD03(I,J,13,NN)=AD03(I,J,13,NN)+(EHg0_bb(I,J)     *DTSRCE)
                AD03(I,J,18,NN)=AD03(I,J,18,NN)+(EHg0_snow(I,J,NN)*DTSRCE)
             ENDDO

          ELSE

             !eds 9/9/10
             DO NN = 1, N_HG_CATS
                AD03(I,J, 1,NN)=AD03(I,J, 1,NN)+(T_Hg_An          *DTSRCE)
                AD03(I,J, 3,NN)=AD03(I,J, 3,NN)+(EHg0_oc(I,J,NN)  *DTSRCE)
                AD03(I,J, 4,NN)=AD03(I,J, 4,NN)+(EHg0_ln(I,J,NN)  *DTSRCE)
                AD03(I,J, 5,NN)=AD03(I,J, 5,NN)+(EHg0_geo(I,J)    *DTSRCE)
                AD03(I,J,13,NN)=AD03(I,J,13,NN)+(EHg0_bb(I,J)     *DTSRCE)
                AD03(I,J,14,NN)=AD03(I,J,14,NN)+(EHg0_vg(I,J)     *DTSRCE)
                AD03(I,J,15,NN)=AD03(I,J,15,NN)+(EHg0_so(I,J)     *DTSRCE)
                AD03(I,J,18,NN)=AD03(I,J,18,NN)+(EHg0_snow(I,J,NN)*DTSRCE)
             ENDDO

          ENDIF

          ! for preindustrial simulation, archive only soil,
          ! CDH- WHY ONLY SOIL??
          ! for present day archive soil, geogenic, biomass burning,
          ! vegetation, and rapid recycing
          IF ( LPREINDHG ) THEN !eds
             AD03_nat(I,J,:)=MAX(EHg0_so(I,J)*DTSRCE, SMALLNUM)
          ELSE
             IF ( LGTMM ) THEN
                DO N = 1, N_HG_CATS
                   AD03_nat(I,J,N) = DTSRCE * ( EHg0_geo(I,J)     + &
                                                EHg0_bb(I,J)      + &
                                                EHg0_gtm(I,J)     + &
                                                EHg0_snow(I,J,N)  )
                ENDDO
             ELSE
                DO N = 1, N_HG_CATS
                   AD03_nat(I,J,N) = DTSRCE * ( EHg0_ln(I,J,N)   + &
                                                EHg0_geo(I,J)    + &
                                                EHg0_bb(I,J)     + &
                                                EHg0_vg(I,J)     + &
                                                EHg0_so(I,J)     + &
                                                EHg0_snow(I,J,N) )
                ENDDO
             ENDIF
          ENDIF

       ENDIF
#endif

       !==============================================================
       ! %%%%% HISTORY (aka netCDF diagnostics) %%%%%
       !
       ! Save the various Hg0 emissions to fields of State_Diag
       ! in units of [kg/s].
       !
       ! NOTES (bmy, 10/24/18):
       ! (1) Normally these would go to the HEMCO_diagnostics*.nc
       !     file, but several of these diagnostics were defined as
       !     HEMCO manual diagnostics.  Therefore to preserve
       !     backwards compatibility with the existing Hg simulation,
       !     we save the emission fields into State_Diag.
       !
       ! (2) For now, save emissions only for total species
       !     and not the tagged species.  Colin Thackray says that
       !     the tagged Hg simulation is in some disrepair and is
       !     not widely used.
       !
       ! (3) The total Hg0 emissions diagnostic is not defined;
       !     users can sum up the individual categories in
       !     post-processing.
       !
       ! (4) In the current Hg simulation, HgP anthro emissions are
       !     added into the Hg2 species.  Therefore we only have
       !     a single emissions field of State_Diag for the
       !     combined Hg2 + HgP anthro emissions.
       !==============================================================

       ! Anthropogenic Hg0 emissions [kg/s]
       IF ( State_Diag%Archive_EmisHg0anthro ) THEN
          State_Diag%EmisHg0anthro(I,J) = EHg0_an(I,J)
       ENDIF

       ! Biomass Hg0 emissions [kg/s]
       IF ( State_Diag%Archive_EmisHg0biomass ) THEN
          State_Diag%EmisHg0biomass(I,J) = EHg0_bb(I,J)
       ENDIF

       ! Geogenic Hg0 emissions [kg/s]
       IF ( State_Diag%Archive_EmisHg0geogenic ) THEN
          State_Diag%EmisHg0geogenic(I,J) = EHg0_geo(I,J)
       ENDIF

       ! Land Hg0 emissions [kg/s]
       IF ( State_Diag%Archive_EmisHg0ocean ) THEN
          State_Diag%EmisHg0land(I,J) = EHg0_ln(I,J,1)
       ENDIF

       ! Oceanic Hg0 emissions [kg/s]
       IF ( State_Diag%Archive_EmisHg0ocean ) THEN
          State_Diag%EmisHg0ocean(I,J) = EHg0_oc(I,J,1)
       ENDIF

       ! Snow Hg0 emissions [kg/s]
       IF ( State_Diag%Archive_EmisHg0snow ) THEN
          State_Diag%EmisHg0snow(I,J) = EHg0_snow(I,J,1)
       ENDIF

       ! Soil Hg0 emissions [kg/s]
       IF ( State_Diag%Archive_EmisHg0soil ) THEN
          State_Diag%EmisHg0soil(I,J) = EHg0_so(I,J)
       ENDIF

       ! Vegetation Hg0 emissions [kg/s]
       IF ( State_Diag%Archive_EmisHg0vegetation ) THEN
          State_Diag%EmisHg0vegetation(I,J) = EHg0_vg(I,J)
       ENDIF

    ENDDO
    ENDDO
    !$OMP END PARALLEL DO

  END SUBROUTINE SRCHg0
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: srcHg2
!
! !DESCRIPTION: Subroutine SRCHg2 is the subroutine for Hg(II) emissions.
!  Emissions of Hg(II) will be distributed throughout the boundary layer.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE SRCHg2( Input_Opt,  State_Chm, State_Diag, &
                     State_Grid, State_Met, RC )
!
! !USES:
!
#ifdef BPCH_DIAG
    USE DIAG03_MOD,         ONLY : AD03, ND03
#endif
    USE ErrCode_Mod
    USE Input_Opt_Mod,      ONLY : OptInput
    USE State_Chm_Mod,      ONLY : ChmState
    USE State_Diag_Mod,     ONLY : DgnState
    USE State_Grid_Mod,     ONLY : GrdState
    USE State_Met_Mod,      ONLY : MetState
    USE TIME_MOD,           ONLY : GET_TS_EMIS
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
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!
! !REVISION HISTORY:
!  07 Dec 2004 - N. (Eckley) Selin - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER  :: I,      J,        L,    N,   PBL_MAX
    REAL(fp) :: DTSRCE, F_OF_PBL, E_Hg

    ! For values from Input_Opt
    LOGICAL  :: LSPLIT, LPREINDHG

    !=================================================================
    ! SRCHg2 begins here!
    !=================================================================

    ! Assume success
    RC        = GC_SUCCESS

    ! Copy values from Input_Opt
    LSPLIT    = Input_Opt%LSPLIT
    LPREINDHG = Input_Opt%LPREINDHG

    ! Emission timestep [s]
    DTSRCE    = GET_TS_EMIS()

    ! Maximum extent of the PBL [model levels]
    PBL_MAX   = State_Met%PBL_MAX_L

    IF (.NOT. LPREINDHG ) THEN

       ! Loop over grid boxes
       !$OMP PARALLEL DO       &
       !$OMP DEFAULT( SHARED ) &
       !$OMP PRIVATE( I, J, L, F_OF_PBL, E_Hg, N )
       DO J = 1, State_Grid%NY
       DO I = 1, State_Grid%NX

          ! Loop up to the max PBL layer
          DO L = 1, PBL_MAX

             ! Fraction of box (I,J,L) w/in the PBL [unitless]
             F_OF_PBL = State_Met%F_OF_PBL(I,J,L)

             ! Partition total Hg2 into box (I,J,L) [kg]
             E_Hg     = F_OF_PBL * EHg2_an(I,J) * DTSRCE

             !---------------------------
             ! Total anthro Hg(II) [kg]
             !---------------------------
             N    = Hg2_Id_List(ID_Hg_tot)
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             !---------------------------
             ! Tagged anthro Hg(II) [kg]
             !---------------------------
             IF ( LSPLIT ) THEN

                ! Anthro Canada Hg2
                N    = Hg2_Id_List(ID_Hg_can)
                E_Hg = F_OF_PBL * EHg2_can(I,J) * DTSRCE
                CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, &
                             State_Grid )

                ! Anthro United States Hg2
                N    = Hg2_Id_List(ID_Hg_usa)
                E_Hg = F_OF_PBL * EHg2_usa(I,J) * DTSRCE
                CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, &
                             State_Grid )

                ! Anthro Central America Hg2
                N    = Hg2_Id_List(ID_Hg_cam)
                E_Hg = F_OF_PBL * EHg2_cam(I,J) * DTSRCE
                CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, &
                             State_Grid )

                ! Anthro South America Hg2
                N    = Hg2_Id_List(ID_Hg_sam)
                E_Hg = F_OF_PBL * EHg2_sam(I,J) * DTSRCE
                CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, &
                             State_Grid )

                ! Anthro West Africa Hg2
                N    = Hg2_Id_List(ID_Hg_waf)
                E_Hg = F_OF_PBL * EHg2_waf(I,J) * DTSRCE
                CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, &
                             State_Grid )

                ! Anthro East Africa Hg2
                N    = Hg2_Id_List(ID_Hg_eaf)
                E_Hg = F_OF_PBL * EHg2_eaf(I,J) * DTSRCE
                CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, &
                             State_Grid )

                ! Anthro South Africa Hg2
                N    = Hg2_Id_List(ID_Hg_saf)
                E_Hg = F_OF_PBL * EHg2_saf(I,J) * DTSRCE
                CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, &
                             State_Grid )

                ! Anthro North Africa Hg2
                N    = Hg2_Id_List(ID_Hg_naf)
                E_Hg = F_OF_PBL * EHg2_naf(I,J) * DTSRCE
                CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, &
                             State_Grid )

                ! Anthro OECD Europe Hg2
                N    = Hg2_Id_List(ID_Hg_eur)
                E_Hg = F_OF_PBL * EHg2_eur(I,J) * DTSRCE
                CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, &
                             State_Grid )

                ! Anthro Eastern Europe Hg2
                N    = Hg2_Id_List(ID_Hg_eeu)
                E_Hg = F_OF_PBL * EHg2_eeu(I,J) * DTSRCE
                CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, &
                             State_Grid )

                ! Anthro Middle East Hg2
                N    = Hg2_Id_List(ID_Hg_mde)
                E_Hg = F_OF_PBL * EHg2_mde(I,J) * DTSRCE
                CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, &
                             State_Grid )

                ! Anthro Former USSR Hg2
                N    = Hg2_Id_List(ID_Hg_sov)
                E_Hg = F_OF_PBL * EHg2_sov(I,J) * DTSRCE
                CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, &
                             State_Grid )

                ! Anthro South Asia Hg2
                N    = Hg2_Id_List(ID_Hg_sas)
                E_Hg = F_OF_PBL * EHg2_sas(I,J) * DTSRCE
                CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, &
                             State_Grid )

                ! Anthro East Asia Hg2
                N    = Hg2_Id_List(ID_Hg_eas)
                E_Hg = F_OF_PBL * EHg2_eas(I,J) * DTSRCE
                CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, &
                             State_Grid )

                ! Anthro Southeast Asia Hg2
                N    = Hg2_Id_List(ID_Hg_sea)
                E_Hg = F_OF_PBL * EHg2_sea(I,J) * DTSRCE
                CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, &
                             State_Grid )

                ! Anthro Japan Hg2
                N    = Hg2_Id_List(ID_Hg_jpn)
                E_Hg = F_OF_PBL * EHg2_jpn(I,J) * DTSRCE
                CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, &
                             State_Grid )

                ! Anthro Oceania Hg2
                N    = Hg2_Id_List(ID_Hg_oce)
                E_Hg = F_OF_PBL * EHg2_oce(I,J) * DTSRCE
                CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, &
                             State_Grid )

             ENDIF
          ENDDO

#ifdef BPCH_DIAG
          !==============================================================
          ! %%%%% ND03 (bpch) DIAGNOSTICS %%%%%
          !
          ! Anthro Hg(II) [kg]
          ! NOTE: HgP is emitted into Hg2
          !==============================================================
          IF ( ND03 > 0 ) THEN
             AD03(I,J,6,1) = AD03(I,J,6,1) + ( EHg2_an(I,J) * DTSRCE )
          ENDIF
#endif

          !==============================================================
          ! %%%%% HISTORY (aka netCDF diagnostics) %%%%%
          !
          ! Save the combined Hg2 + HgP emissions to State_Diag
          ! in units of [kg/s].
          !
          ! NOTES (bmy, 10/24/18):
          !
          ! (1) In the current Hg simulation, HgP anthro emissions
          !     are added into the Hg2 species.  Therefore we only
          !     have defined a single emissions field of State_Diag
          !     for the combined Hg2 + HgP anthro emissions.
          !==============================================================
          IF ( State_Diag%Archive_EmisHg2HgPanthro ) THEN
             State_Diag%EmisHg2HgPanthro(I,J) = EHg2_an(I,J)
          ENDIF

       ENDDO
       ENDDO
       !$OMP END PARALLEL DO

    ENDIF

  END SUBROUTINE SRCHg2
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: srcHgp
!
! !DESCRIPTION: Subroutine SRCHgP is the subroutine for HgP emissions.
!  Emissions of HgP will be distributed throughout the boundary layer.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE SRCHgP( Input_Opt, State_Chm, State_Grid, State_Met, RC )
!
! !USES:
!
#ifdef BPCH_DIAG
    USE DIAG03_MOD,         ONLY : AD03, ND03
#endif
    USE ErrCode_Mod
    USE Input_Opt_Mod,      ONLY : OptInput
    USE State_Chm_Mod,      ONLY : ChmState
    USE State_Grid_Mod,     ONLY : GrdState
    USE State_Met_Mod,      ONLY : MetState
    USE TIME_MOD,           ONLY : GET_TS_EMIS
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
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!
! !REVISION HISTORY:
!  07 Dec 2004 - N. (Eckley) Selin - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER  :: I,      J,        L,    N,   PBL_MAX
    REAL(fp) :: DTSRCE, F_OF_PBL, E_Hg

    ! For values from Input_Opt
    LOGICAL  :: LSPLIT, LPREINDHG

    !=================================================================
    ! SRCHgP begins here!
    !=================================================================

    ! Assume success
    RC        = GC_SUCCESS

    ! Copy values from Input_Opt
    LSPLIT    = Input_Opt%LSPLIT
    LPREINDHG = Input_Opt%LPREINDHG

    ! Chemistry timestep [s]
    DTSRCE    = GET_TS_EMIS()

    ! Maximum extent of the PBL [model levels]
    PBL_MAX   = State_Met%PBL_MAX_L

    IF (.NOT. LPREINDHG) THEN

       !$OMP PARALLEL DO       &
       !$OMP DEFAULT( SHARED ) &
       !$OMP PRIVATE( I, J, L, F_OF_PBL, E_Hg, N )
       DO J = 1, State_Grid%NY
       DO I = 1, State_Grid%NX

          ! Loop up to PBL top layer
          DO L = 1, PBL_MAX

             ! Fraction of box (I,J,L) w/in the PBL [unitless]
             F_OF_PBL           = State_Met%F_OF_PBL(I,J,L)

             ! Partition HgP into box (I,J,L) [kg]
             E_Hg               = F_OF_PBL * EHgP_an(I,J) * DTSRCE

             !------------------------
             ! Total anthro HgP [kg]
             !------------------------
             N               = HgP_Id_List(ID_Hg_tot)
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, &
                          State_Grid )

             !------------------------
             ! Tagged anthro HgP [kg]
             !------------------------
             IF ( LSPLIT ) THEN

                ! Anthro Canada HgP
                N            = HgP_Id_List(ID_Hg_can)
                E_Hg         = F_OF_PBL * EHgP_can(I,J) * DTSRCE
                CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, &
                             State_Grid )

                ! Anthro United States HgP
                N            = HgP_Id_List(ID_Hg_usa)
                E_Hg         = F_OF_PBL * EHgP_usa(I,J) * DTSRCE
                CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, &
                             State_Grid )
                
                ! Anthro Central America HgP
                N            = HgP_Id_List(ID_Hg_cam)
                E_Hg         = F_OF_PBL * EHgP_cam(I,J) * DTSRCE
                CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, &
                             State_Grid )
                
                ! Anthro South America HgP
                N            = HgP_Id_List(ID_Hg_sam)
                E_Hg         = F_OF_PBL * EHgP_sam(I,J) * DTSRCE
                CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, &
                             State_Grid )
                
                ! Anthro West Africa HgP
                N            = HgP_Id_List(ID_Hg_waf)
                E_Hg         = F_OF_PBL * EHgP_waf(I,J) * DTSRCE
                CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, &
                             State_Grid )
                
                ! Anthro East Africa HgP
                N            = HgP_Id_List(ID_Hg_eaf)
                E_Hg         = F_OF_PBL * EHgP_eaf(I,J) * DTSRCE
                CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, &
                             State_Grid )
                
                ! Anthro South Africa HgP
                N            = HgP_Id_List(ID_Hg_saf)
                E_Hg         = F_OF_PBL * EHgP_saf(I,J) * DTSRCE
                CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, &
                             State_Grid )
                
                ! Anthro North Africa HgP
                N            = HgP_Id_List(ID_Hg_naf)
                E_Hg         = F_OF_PBL * EHgP_naf(I,J) * DTSRCE
                CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, &
                             State_Grid )
                
                ! Anthro OECD Europe HgP
                N            = HgP_Id_List(ID_Hg_eur)
                E_Hg         = F_OF_PBL * EHgP_eur(I,J) * DTSRCE
                CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, &
                             State_Grid )
                
                ! Anthro Eastern Europe HgP
                N            = HgP_Id_List(ID_Hg_eeu)
                E_Hg         = F_OF_PBL * EHgP_eeu(I,J) * DTSRCE
                CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, &
                             State_Grid )
                
                ! Anthro Middle East HgP
                N            = HgP_Id_List(ID_Hg_mde)
                E_Hg         = F_OF_PBL * EHgP_mde(I,J) * DTSRCE
                CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, &
                             State_Grid )
                
                ! Anthro Former USSR HgP
                N            = HgP_Id_List(ID_Hg_sov)
                E_Hg         = F_OF_PBL * EHgP_sov(I,J) * DTSRCE
                CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, &
                             State_Grid )
                
                ! Anthro South Asia HgP
                N            = HgP_Id_List(ID_Hg_sas)
                E_Hg         = F_OF_PBL * EHgP_sas(I,J) * DTSRCE
                CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, &
                             State_Grid )
                
                ! Anthro East Asia HgP
                N            = HgP_Id_List(ID_Hg_eas)
                E_Hg         = F_OF_PBL * EHgP_eas(I,J) * DTSRCE
                CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, &
                             State_Grid )
                
                ! Anthro Southeast Asia HgP
                N            = HgP_Id_List(ID_Hg_sea)
                E_Hg         = F_OF_PBL * EHgP_sea(I,J) * DTSRCE
                CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, &
                             State_Grid )
                
                ! Anthro Japan HgP
                N            = HgP_Id_List(ID_Hg_jpn)
                E_Hg         = F_OF_PBL * EHgP_jpn(I,J) * DTSRCE
                CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, &
                             State_Grid )
                
                ! Anthro Oceania HgP
                N            = HgP_Id_List(ID_Hg_oce)
                E_Hg         = F_OF_PBL * EHgP_oce(I,J) * DTSRCE
                CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, &
                             State_Grid )

             ENDIF
          ENDDO

#ifdef BPCH_DIAG
          !-----------------------------------------------------------
          ! %%%%% ND03 (bpch) diagnostics
          !
          ! Anthro HgP [kg]
          ! NOTE: This is zero, because in the current Hg simulation
          ! HgP is emitted into the Hg2 species (bmy, 10/26/18)
          !-----------------------------------------------------------
          IF ( ND03 > 0 ) THEN
             AD03(I,J,9,1) = AD03(I,J,9,1) + (EHgP_an(I,J) * DTSRCE)
          ENDIF
#endif

       ENDDO
       ENDDO
       !$OMP END PARALLEL DO

    ENDIF

  END SUBROUTINE SRCHgP
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: mercury_readyr
!
! !DESCRIPTION: Subroutine MERCURY\_READYR reads the year-invariant emissions
!  for Mercury from anthropogenic, ocean, and land sources.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE MERCURY_READYR( Input_Opt, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE HCO_ERROR_MOD
    USE HCO_INTERFACE_MOD,  ONLY : HcoState
    USE HCO_EMISLIST_MOD,   ONLY : HCO_GetPtr
    USE HCO_INTERFACE_MOD,  ONLY : GetHcoDiagn
    USE Input_Opt_Mod,      ONLY : OptInput
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!
! !REVISION HISTORY:
!  06 Dec 2004 - N. (Eckley) Selin - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS
!
    REAL(fp), PARAMETER      :: SEC_PER_YR = 365.25e+0_fp*86400e+0_fp
!
! !LOCAL VARIABLES:
!
    ! For values from Input_Opt
    LOGICAL                  :: LDYNOCEAN
    LOGICAL                  :: LPREINDHG
    LOGICAL                  :: LSPLIT

    ! HEMCO update
    CHARACTER(LEN=63)        :: DgnName
    REAL(f4),        POINTER :: Ptr2D(:,:)

    ! Strings
    CHARACTER(LEN=255)       :: ThisLoc
    CHARACTER(LEN=512)       :: ErrMsg

    !=================================================================
    ! MERCURY_READYR begins here!
    !=================================================================

    ! Assume success
    RC        = GC_SUCCESS
    ErrMsg    = ''
    ThisLoc   = '-> at MERCURY_READYR (in GeosCore/mercury_mod.F90)'

    ! Initialize pointers
    Ptr2D     => NULL()

    ! Copy values from Input_Opt
    LDYNOCEAN = Input_Opt%LDYNOCEAN
    LPREINDHG = Input_Opt%LPREINDHG
    LSPLIT    = Input_Opt%LSPLIT

    !=================================================================
    ! Anthropogenic Emissions
    !=================================================================

    !=================================================================
    ! Now get the emission fields through HEMCO. The emissions
    ! settings and source data are now set in the HEMCO configuration
    ! file (specified in input.geos) and emissions are calculated
    ! based on the content of the configuration file.
    ! Here, we just import the HEMCO diagnostics fields (defined
    ! in hcoi_gc_diagn_mod.F90) and pass the emission fields to
    ! the corresponding arrays (EHg0_an, etc.). Those internal
    ! arrays are in kg/s, so we need to convert the diagnostics
    ! from kg/m2/s to kg/s.
    !
    ! NOTE: the current implementation is a hybrid version between
    ! HEMCO and the old mercury code, primarily to make sure that
    ! all the diagnostics based on EHg* still work. Because of
    ! that, we convert all emissions from kg/m2/s to kg/s and then
    ! back to kg/m2/s when we pass them to HG_EMIS (via SRCHg0 and
    ! SRCHg2). HG_EMIS is used for the non-local PBL mixing scheme.
    ! (ckeller, 09/24/2014)
    !=================================================================

    ! Get HEMCO state object
    IF ( .NOT. ASSOCIATED(HcoState) ) THEN
       ErrMsg = 'HcoState not associated!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    ! No anthropogenic emissions in preindustrial simulation
    IF ( LPREINDHG ) THEN

       EHg0_an = 0e+0_fp  !eds
       EHg2_an = 0e+0_fp
       EHgP_an = 0e+0_fp

    ELSE

       ! ---------------
       ! Hg0 emissions
       ! ---------------

       ! Anthropogenic emissions
       DgnName = 'HG0_ANTHRO'
       CALL GetHcoDiagn( DgnName, .TRUE., RC, Ptr2D=Ptr2D )
       IF ( RC /= HCO_SUCCESS ) THEN
          ErrMsg = 'Could not get HEMCO field ' // TRIM( DgnName )
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
       EHg0_an =  Ptr2D(:,:) * HcoState%Grid%AREA_M2%Val(:,:)
       Ptr2D   => NULL()

       ! Artisanal emissions
       DgnName = 'HG0_ARTISANAL'
       CALL GetHcoDiagn( DgnName, .TRUE., RC, Ptr2D=Ptr2D )
       IF ( RC /= HCO_SUCCESS ) THEN
          ErrMsg = 'Could not get HEMCO field ' // TRIM( DgnName )
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
       EHg0_am =  Ptr2D(:,:) * HcoState%Grid%AREA_M2%Val(:,:)
       Ptr2D   => NULL()

       ! ---------------------------------------------
       ! Hg2 emissions
       ! Note: HgP emissions are added to Hg2 by HEMCO
       ! ---------------------------------------------

       ! Anthropogenic emissions
       DgnName = 'HG2_ANTHRO'
       CALL GetHcoDiagn( DgnName, .TRUE., RC, Ptr2D=Ptr2D )
       IF ( RC /= HCO_SUCCESS ) THEN
          ErrMsg = 'Could not get HEMCO field ' // TRIM( DgnName )
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
       EHg2_an =  Ptr2D(:,:) * HcoState%Grid%AREA_M2%Val(:,:)
       Ptr2D   => NULL()

    ENDIF

    ! Natural emissions
    DgnName = 'HG0_NATURAL'
    CALL GetHcoDiagn( DgnName, .TRUE., RC, Ptr2D=Ptr2D )
    IF ( RC /= HCO_SUCCESS ) THEN
       ErrMsg = 'Could not get HEMCO field ' // TRIM( DgnName )
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    EHg0_geo =  Ptr2D(:,:) * HcoState%Grid%AREA_M2%Val(:,:)
    Ptr2D    => NULL()

    ! Soil distribution
    CALL HCO_GetPtr( HcoState, 'HG0_SOILDIST', Ptr2D, RC )
    IF ( RC /= HCO_SUCCESS ) THEN
       ErrMsg = 'Could not get pointer to HEMCO field HG0_SOILDIST!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    EHg0_dist =  Ptr2D(:,:)
    Ptr2D     => NULL()

    !=================================================================
    ! Print totals to the screen in [Gg/yr]
    !=================================================================
    WRITE( 6, '(a)' ) REPEAT( '=', 79 )
    WRITE( 6, 210   )
    WRITE( 6, '(a)' )
    WRITE( 6, 211   ) SUM( EHg0_an ) * SEC_PER_YR * 1e-6_fp
    WRITE( 6, 213   ) SUM( EHg0_ln ) * SEC_PER_YR * 1e-6_fp
    WRITE( 6, 214   ) SUM( EHg0_geo ) * SEC_PER_YR * 1e-6_fp

    ! Only write ocean total if we are doing offline ocean
    IF ( .not. LDYNOCEAN ) THEN
       WRITE( 6, 217   ) SUM( EHg0_oc ) * SEC_PER_YR * 1e-6_fp
    ENDIF

    WRITE( 6, 215   ) SUM( EHg2_an ) * SEC_PER_YR * 1e-6_fp
    WRITE( 6, 216   ) SUM( EHgP_an ) * SEC_PER_YR * 1e-6_fp
    WRITE( 6, '(a)' ) REPEAT( '=', 79 )

    ! FORMAT strings
210 FORMAT( 'M E R C U R Y   E M I S S I O N S' )
211 FORMAT( 'Total Anthro     Hg(0)  : ', f7.3, ' [Gg/yr]' )
213 FORMAT( 'Total Re-Emitted Hg(0)  : ', f7.3, ' [Gg/yr]' )
214 FORMAT( 'Total Natural    Hg(0)  : ', f7.3, ' [Gg/yr]' )
215 FORMAT( 'Total Anthro     Hg(II) : ', f7.3, ' [Gg/yr]' )
216 FORMAT( 'Total Anthro     HgP    : ', f7.3, ' [Gg/yr]' )
217 FORMAT( 'Total Ocean      Hg(0)  : ', f7.3, ' [Gg/yr]' )

  END SUBROUTINE MERCURY_READYR
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_o3
!
! !DESCRIPTION: Function GET\_O3 returns monthly mean O3 for the mercury
!  simulation.
!\\
!\\
! !INTERFACE:
!
  FUNCTION GET_O3( I, J, L, State_Met ) RESULT( O3_MOLEC_CM3 )
!
! !USES:
!
    USE State_Met_Mod,      ONLY : MetState
!
! !INPUT PARAMETERS:
!
    INTEGER,        INTENT(IN)  :: I              ! Longitude index
    INTEGER,        INTENT(IN)  :: J              ! Latitude index
    INTEGER,        INTENT(IN)  :: L              ! Level index
    TYPE(MetState), INTENT(IN)  :: State_Met      ! Meteorology State object
!
! !RETURN VALUE:
!
    REAL(fp)                    :: O3_MOLEC_CM3   ! O3 conc [molec/cm3]
!
! !REVISION HISTORY:
!  06 Dec 2002 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    REAL(fp) :: BOXVL

    !=================================================================
    ! GET_O3 begins here!
    !=================================================================

    ! Grid box volume [cm3]
    BOXVL        = State_Met%AIRVOL(I,J,L) * 1e+6_fp

    ! Get ozone [v/v] for this gridbox & month
    ! and convert to [molec/cm3] (eck, 12/2/04)
    O3_MOLEC_CM3 = O3(I,J,L)       * ( AVO * 1e+3_fp / AIRMW ) * &
                   State_Met%AD(I,J,L) / BOXVL

  END FUNCTION GET_O3
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_oh
!
! !DESCRIPTION: Function GET\_OH returns monthly mean OH and imposes a
!  diurnal variation.
!\\
!\\
! !INTERFACE:
!
  FUNCTION GET_OH( I, J, L, State_Met ) RESULT( OH_MOLEC_CM3 )
!
! !USES:
!
    USE State_Met_Mod,      ONLY : MetState
    USE TIME_MOD,           ONLY : GET_TS_CHEM
!
! !INPUT PARAMETERS:
!
    INTEGER,        INTENT(IN)  :: I              ! Longitude index
    INTEGER,        INTENT(IN)  :: J              ! Latitude index
    INTEGER,        INTENT(IN)  :: L              ! Level index
    TYPE(MetState), INTENT(IN)  :: State_Met      ! Meteorology State object
!
! !RETURN VALUE:
!
    REAL(fp)                    :: OH_MOLEC_CM3   ! OH conc [molec/cm3]
!
! !REVISION HISTORY:
!  07 Dec 2004 - N. (Eckley) Selin - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    REAL(fp) :: TPL
    REAL(fp) :: HEM_OH

    !=================================================================
    ! GET_OH begins here!
    !=================================================================

    TPL = State_Met%TropLev(I, J)

    IF ( L > TPL ) THEN
       HEM_OH = OH_strat(I,J,L)
    ELSE
       HEM_OH = OH_trop(I,J,L)
    ENDIF

    ! Test for sunlight...
    IF ( State_Met%SUNCOS(I,J) > 0e+0_fp .and.  TCOSZ(I,J) > 0e+0_fp ) THEN

       ! Impose a diurnal variation on OH during the day
       ! OH from HEMCO is in kg/m3
       OH_MOLEC_CM3 = HEM_OH * 1e-6_fp * ( AVO / 0.017) * &
                      ( State_Met%SUNCOS(I,J)  / TCOSZ(I,J)    ) * &
                      ( 86400e+0_fp            / GET_TS_CHEM() )

       ! Make sure OH is not negative
       OH_MOLEC_CM3 = MAX( OH_MOLEC_CM3, 0e+0_fp )

    ELSE

       ! At night, OH goes to zero
       OH_MOLEC_CM3 = 0e+0_fp

    ENDIF

  END FUNCTION GET_OH
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_ho2
!
! !DESCRIPTION: Function GET\_HO2 returns monthly mean HO2 and imposes a
!  diurnal variation. HMH 9/24/14 adding for dibble chemistry
!\\
!\\
! !INTERFACE:
!
  FUNCTION GET_HO2( I, J, L, State_Met ) RESULT( HO2_MOLEC_CM3 )
!
! !USES:
!
    USE State_Met_Mod, ONLY : MetState
    USE TIME_MOD,      ONLY : GET_TS_CHEM
!
! !INPUT PARAMETERS:
!
    INTEGER,        INTENT(IN)  :: I              ! Longitude index
    INTEGER,        INTENT(IN)  :: J              ! Latitude index
    INTEGER,        INTENT(IN)  :: L              ! Level index
    TYPE(MetState), INTENT(IN)  :: State_Met      ! Meteorology State object
!
! !RETURN VALUE:
!
    REAL(fp)                    :: HO2_MOLEC_CM3 ! HO2 conc [molec/cm3]
!
! !REVISION HISTORY:
!  07 Dec 2004 - N. (Eckley) Selin - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    REAL(fp) :: HO2_VV
    REAL(fp) :: HEM_HO2
    REAL(fp) :: TPL

    !=================================================================
    ! GET_HO2 begins here!
    !=================================================================

    TPL = State_Met%TropLev(I, J)

    IF ( L > TPL ) THEN
       HEM_HO2 = HEM_HO2_strat(I,J,L)
    ELSE
       HEM_HO2 = HEM_HO2_trop(I,J,L)
    ENDIF

    ! Test for sunlight...
    IF ( State_Met%SUNCOS(I,J) > 0e0_fp .and. TCOSZ(I,J) > 0e0_fp ) THEN

       ! Impose a diurnal variation on HO2 during the day
       HO2_VV = HEM_HO2 * ( State_Met%SUNCOS(I,J) / TCOSZ(I,J)    ) * &
                          ( 86400e0_fp            / GET_TS_CHEM() )

       ! Make sure HO2 is not negative
       HO2_VV = MAX( HO2_VV, 0e0_fp )
       HO2_MOLEC_CM3 = HO2_VV * State_Met%AIRNUMDEN(I,J,L)

    ELSE

       ! At night, HO2 goes to zero
       HO2_MOLEC_CM3 = 0e0_fp

    ENDIF

  END FUNCTION GET_HO2
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_clo
!
! !DESCRIPTION: Function GET\_CLO imposes a diurnal scaling for ClO, hmh 10/6/
!               adding for dibble chemistry. Now including scaling of Cl and
!               HOCl for aqueous and gas-phase oxidation 11/19/14.
!\\
!\\
! !INTERFACE:
!
  FUNCTION GET_ClO( I, J, L, State_Grid, State_Met ) &
       RESULT( ClO_MOLEC_CM3 )
!
! !USES:
!
    USE State_Grid_Mod,     ONLY : GrdState
    USE State_Met_Mod,      ONLY : MetState
    USE ERROR_MOD,          ONLY : SAFE_DIV
    USE TIME_MOD,           ONLY : GET_TS_CHEM, GET_LOCALTIME
!
! !INPUT PARAMETERS:
!
    INTEGER,        INTENT(IN)  :: I              ! Longitude index
    INTEGER,        INTENT(IN)  :: J              ! Latitude index
    INTEGER,        INTENT(IN)  :: L              ! Level index
    TYPE(GrdState), INTENT(IN)  :: State_Grid     ! Grid State object
    TYPE(MetState), INTENT(IN)  :: State_Met      ! Meteorology State object
!
! !RETURN VALUE:
!
    REAL(fp)                    :: ClO_MOLEC_CM3  ! ClO  conc [molec/cm3]
!
! !REVISION HISTORY:
!  06 Oct 2014 - H. Horowitz - initial version copied from get_br.
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    REAL(fp) :: ClO_PPBV
    REAL*8   :: ClO_FAC, LOCALTIME, HOUR

    !=================================================================
    ! GET_ClO begins here! copying from get_br!
    !=================================================================

    !----------------------------------------------------------------
    ! Impose a diurnal cycle
    !
    ! This formula comes from the Holmes et al (2009) MBL box model
    ! which fit RGM observations for Okinawa and a Pacific cruise
    !----------------------------------------------------------------

    ! Test for sunlight...
    IF ( (State_Met%SUNCOS(I,J) > 0e0_fp) .and. (TTDAY(I,J) > 0e0_fp) ) THEN

       ! Use a constant function if daylight is < 2 hours or > 22hours
       IF ( ( TTDAY(I,J) < 120e0_fp  ) .OR. ( TTDAY(I,J) > 1320e0_fp ) ) THEN

          ClO_FAC = SAFE_DIV( 1440e0_fp, TTDAY(I,J), 0E0_FP )

       ELSE

          ! Local time: 0-23.999
          LOCALTIME = GET_LOCALTIME( I, J, L, State_Grid )

          ! Interpolate the real hour to lie within an equinoctal day
          ! i.e. between 6-18
          HOUR = 12E0_FP + ( LOCALTIME - 12E0_FP ) * 720E0_FP / ( TTDAY(I,J) )

          ! Multiplicative factor for the diurnal cycle
          ClO_FAC = ( 1e6_fp - 1e2_fp * ( HOUR - 6E0_FP ) ** 4E0_FP + &
                    ( 18E0_FP - HOUR ) * 1e6_fp / 12E0_FP ) / 2e6_fp

          ! Normalize the multiplicative factor to have a 24-h mean of 1
          ClO_FAC = SAFE_DIV( ClO_FAC, (4e-4_fp * TTDAY(I,J) ),0E0_FP)

       ENDIF

       ! Make sure that diurnal scaling factor is non-negative
       ClO_FAC = MAX( ClO_FAC, 0E0_FP )

       ! The instantaneous concentration is the 24-h mean times
       ! the time-varying factor from the diurnal cycle
       ClO_PPBV  = HEM_CLO(I,J,L) * ClO_FAC
       ClO_MOLEC_CM3 = ClO_PPBV * 1.0e-9_fp * State_Met%AIRNUMDEN(I,J,L)

    ELSE

       ! At night, ClO, HOCl, Cl goes to zero
       ClO_MOLEC_CM3  = 0E0_FP
    ENDIF

  END FUNCTION GET_ClO
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_cl
!
! !DESCRIPTION: Function GET\_CL imposes a diurnal scaling for Cl, hmh 11/19/14
!               adding for Cl-initiated oxidation.
!\\
!\\
! !INTERFACE:
!
  FUNCTION GET_Cl( I, J, L, State_Grid, State_Met ) &
       RESULT( Cl_MOLEC_CM3 )
!
! !USES:
!
    USE State_Grid_Mod,     ONLY : GrdState
    USE State_Met_Mod,      ONLY : MetState
    USE ERROR_MOD,          ONLY : SAFE_DIV
    USE TIME_MOD,           ONLY : GET_TS_CHEM, GET_LOCALTIME
!
! !INPUT PARAMETERS:
!
    INTEGER,        INTENT(IN)  :: I              ! Longitude index
    INTEGER,        INTENT(IN)  :: J              ! Latitude index
    INTEGER,        INTENT(IN)  :: L              ! Level index
    TYPE(GrdState), INTENT(IN)  :: State_Grid     ! Grid State object
    TYPE(MetState), INTENT(IN)  :: State_Met      ! Meteorology State object
!
! !RETURN VALUE:
!
    REAL(fp)                    :: Cl_MOLEC_CM3   ! Cl   conc [molec/cm3]
!
! !REVISION HISTORY:
!  19 Nov 2014 - H. Horowitz - initial version copied from get_clO.
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    REAL(fp) :: Cl_PPBV
    REAL*8   :: Cl_FAC, LOCALTIME, HOUR

    !=================================================================
    ! GET_Cl begins here! copying from get_ClO!
    !=================================================================

    !----------------------------------------------------------------
    ! Impose a diurnal cycle
    !
    ! This formula comes from the Holmes et al (2009) MBL box model
    ! which fit RGM observations for Okinawa and a Pacific cruise
    !----------------------------------------------------------------

    ! Test for sunlight...
    IF ( (State_Met%SUNCOS(I,J) > 0e0_fp) .and. (TTDAY(I,J) > 0e0_fp) ) THEN

       ! Use a constant function if daylight is < 2 hours or > 22hours
       IF ( ( TTDAY(I,J) < 120e0_fp  ) .OR. ( TTDAY(I,J) > 1320e0_fp ) ) THEN

          Cl_FAC = SAFE_DIV( 1440e0_fp, TTDAY(I,J), 0E0_FP )

       ELSE

          ! Local time: 0-23.999
          LOCALTIME = GET_LOCALTIME( I, J, L, State_Grid )

          ! Interpolate the real hour to lie within an equinoctal day
          ! i.e. between 6-18
          HOUR = 12E0_FP + ( LOCALTIME - 12E0_FP ) * 720E0_FP / ( TTDAY(I,J) )

          ! Multiplicative factor for the diurnal cycle
          Cl_FAC = ( 1e6_fp - 1e2_fp * ( HOUR - 6E0_FP ) ** 4E0_FP + &
                   ( 18E0_FP - HOUR ) * 1e6_fp / 12E0_FP ) / 2e6_fp

          ! Normalize the multiplicative factor to have a 24-h mean of 1
          Cl_FAC = SAFE_DIV( Cl_FAC, (4e-4_fp * TTDAY(I,J) ), 0E0_FP )

       ENDIF

       ! Make sure that diurnal scaling factor is non-negative
       Cl_FAC = MAX( Cl_FAC, 0E0_FP )

       ! The instantaneous concentration is the 24-h mean times
       ! the time-varying factor from the diurnal cycle
       Cl_PPBV   = HEM_CL(I,J,L) * Cl_FAC
       Cl_MOLEC_CM3 = Cl_PPBV * 1.0e-9_fp * State_Met%AIRNUMDEN(I,J,L)

    ELSE

       ! At night, ClO, HOCl, Cl goes to zero
       Cl_MOLEC_CM3   = 0E0_FP

    ENDIF

  END FUNCTION GET_Cl
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_HOCl
!
! !DESCRIPTION: Function GET\_HOCl returns instantaneous HOCl concentration
!  calculated from the monthly mean and an imposed diurnal variation.
!\\
!\\
! !INTERFACE:

  FUNCTION GET_HOCl( I, J, L, State_Grid, State_Met ) &
       RESULT( HOCL_INST )
!
! !USES:
!
    USE ERROR_MOD,      ONLY : SAFE_DIV
    USE TIME_MOD,       ONLY : GET_LOCALTIME
    USE State_Grid_Mod, ONLY : GrdState
    USE State_Met_Mod,  ONLY : MetState
!
! !INPUT PARAMETERS:
!
    INTEGER,        INTENT(IN) :: I, J, L     ! Grid box indices
    TYPE(GrdState), INTENT(IN) :: State_Grid  ! Grid State object
    TYPE(MetState), INTENT(IN) :: State_Met   ! Meteorology State object
!
! !RETURN VALUE:
!
    REAL(fp)                   :: HOCL_INST
!
! !REVISION HISTORY:
!  ca. 2014 - H. Horowitz - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    REAL(fp)            :: HOCl_FAC, LOCALTIME, HOUR

    !=================================================================
    ! GET_HOCl begins here! copying from get_ClO!
    !=================================================================

    !----------------------------------------------------------------
    ! Impose a diurnal cycle
    !
    ! This formula comes from the Holmes et al (2009) MBL box model
    ! which fit RGM observations for Okinawa and a Pacific cruise
    !----------------------------------------------------------------

    ! Test for sunlight...
    IF ( (State_Met%SUNCOS(I,J) > 0e0_fp) .and. (TTDAY(I,J) > 0e0_fp) ) THEN

       ! Use a constant function if daylight is < 2 hours or > 22hours
       IF ( ( TTDAY(I,J) < 120e0_fp  ) .OR. ( TTDAY(I,J) > 1320e0_fp ) ) THEN

          HOCl_FAC = SAFE_DIV( 1440e0_fp, TTDAY(I,J), 0E0_FP )

       ELSE

          ! Local time: 0-23.999
          LOCALTIME = GET_LOCALTIME( I, J, L, State_Grid )

          ! Interpolate the real hour to lie within an equinoctal day
          ! i.e. between 6-18
          HOUR = 12E0_FP + ( LOCALTIME - 12E0_FP ) * 720E0_FP / ( TTDAY(I,J) )

          ! Multiplicative factor for the diurnal cycle
          HOCl_FAC = ( 1e6_fp - 1e2_fp * ( HOUR - 6E0_FP ) ** 4E0_FP + &
                     ( 18E0_FP - HOUR ) * 1e6_fp / 12E0_FP ) / 2e6_fp

          ! Normalize the multiplicative factor to have a 24-h mean of 1
          HOCl_FAC = SAFE_DIV( HOCl_FAC, (4e-4_fp *TTDAY(I,J)),0E0_FP)

       ENDIF

       ! Make sure that diurnal scaling factor is non-negative
       HOCl_FAC = MAX( HOCl_FAC, 0E0_FP )

       HOCl_INST = HEM_HOCL(I,J,L) * HOCl_FAC ! ppbv
       HOCl_INST = HOCl_INST * 1.0e-9_fp * State_Met%AIRNUMDEN(I,J,L)

    ELSE

       ! Concentration 0 in the dark
       HOCl_INST = 0e0_fp

    ENDIF

  END FUNCTION GET_HOCl
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_NO2
!
! !DESCRIPTION: Function GET\_NO2 returns instantaneous NO2 concentration
!  calculated from the monthly mean and an imposed diurnal variation.
!\\
!\\
! !INTERFACE:
!
  FUNCTION GET_NO2( I, J, L, State_Grid, State_Met, NO_MO, NO2_MO, O3_MO ) &
       RESULT( NO2_INST )
!
! !USES:
!
    USE ERROR_MOD,      ONLY : SAFE_DIV
    USE State_Grid_Mod, ONLY : GrdState
    USE State_Met_Mod,  ONLY : MetState
!
! !INPUT PARAMETERS:
!
    INTEGER,        INTENT(IN) :: I,J,L
    TYPE(GrdState), INTENT(IN) :: State_Grid
    TYPE(MetState), INTENT(IN) :: State_Met
    REAL(4),        INTENT(IN) :: NO_MO
    REAL(4),        INTENT(IN) :: NO2_MO
    REAL(fp),       INTENT(IN) :: O3_MO
!
! !RETURN VALUE:
!
    REAL(fp) :: NO2_INST
!
! !REVISION HISTORY:
!  ca. 2014    - H. Horowitz - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    REAL(fp) :: J_NO2
    REAL(fp) :: k3
    REAL(fp) :: F_NO2
    REAL(fp) :: A
    REAL(fp) :: EdivR
    REAL(fp) :: NOx_MO
    REAL(fp) :: k3_O3_MO_prod, J_NO2_k3_O3_MO_sum

    !  =======================================================================
    !  NOTES:
    !
    !  (R1) NO2 + hv -> NO + O,  j
    !  (R2) O +O2 -> O3,         k2
    !  (R3) O3 + NO -> NO2 + O2, k3
    !
    !  [NOx] = [NO] + [NO2]
    !
    !  NOx steady state:
    !  j[NO2] = k3[O3][NO]
    !  j[NO2] = k3[O3]([NOx]-[NO2])
    !  [NO2](j+k3[O3] = k3[O3][NOx]
    !  [NO2]/[NOx] = k3[O3]/(j+k3[O3])
    !
    !   k3 = A exp(-E / RT) Arrhenius Equation
    !   A = 3.0e-12 Seinfeld & Pandis
    !   E/R = 1500
    !   k3 = 1.0e-11 @ 298k
    !
    !***************************************************************************
!
! !DEFINED PARAMETERS:
!
    A     = 3.0e-12_fp
    EdivR = 1.5e+3_fp

    !=================================================================
    ! GET_NO2 begins here!
    !=================================================================

    k3 = A*exp(-EdivR/State_Met%T(I,J,L))

    !Get total monthly NOx
    NOx_MO = NO_MO + NO2_MO

    J_NO2 = GET_JNO2( I, J, L, State_Grid, State_Met )
    k3_O3_MO_prod = k3*O3_MO
    J_NO2_k3_O3_MO_sum = J_NO2 + k3_O3_MO_prod
    F_NO2 = SAFE_DIV( k3_O3_MO_prod, J_NO2_k3_O3_MO_sum, 0E0_FP )

    IF (O3_MO .LT. SMALLNUM) THEN
       F_NO2 = SMALLNUM
    ENDIF

    IF (J_NO2 .LT. SMALLNUM) THEN
       F_NO2 = 1e+0_fp
    ENDIF

    NO2_INST = F_NO2 * NOx_MO ! ppbv
    NO2_INST = NO2_INST * 1.0e-9_fp * State_Met%AIRNUMDEN(I,J,L)

  END FUNCTION GET_NO2
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_Br
!
! !DESCRIPTION: Function GET\_BR returns instantaneous Br concentration
!  calculated from the monthly mean and an imposed diurnal variation.
!\\
!\\
! !INTERFACE:
!
  FUNCTION GET_BR( I, J, L, BRO_MOLEC_CM3, State_Diag, State_Grid, State_Met ) &
       RESULT( BR_MOLEC_CM3 )
!
! !USES:
!
    USE Calc_Met_Mod,       ONLY : GET_OBK
#ifdef BPCH_DIAG
    USE DIAG03_MOD,         ONLY : AD03_Br,   LD03,     ND03
#endif
    USE ERROR_MOD,          ONLY : SAFE_DIV
    USE GLOBAL_BR_MOD,      ONLY : BR_MERGE,  BR_TROP
    USE GLOBAL_BR_MOD,      ONLY : BR_STRAT,  J_BRO
    USE GLOBAL_BR_MOD,      ONLY : BRO_MERGE, BRO_TROP, BRO_STRAT
    USE State_Diag_Mod,     ONLY : DgnState
    USE State_Grid_Mod,     ONLY : GrdState
    USE State_Met_Mod,      ONLY : MetState
    USE TIME_MOD,           ONLY : GET_LOCALTIME
    USE TIME_MOD,           ONLY : ITS_A_NEW_DAY
    USE TIME_MOD,           ONLY : GET_MONTH
!
! !INPUT PARAMETERS:
!
    INTEGER,        INTENT(IN)    :: I            ! Longitude index
    INTEGER,        INTENT(IN)    :: J            ! Latitude index
    INTEGER,        INTENT(IN)    :: L            ! Level index
    TYPE(GrdState), INTENT(IN)    :: State_Grid   ! Grid State object
    TYPE(MetState), INTENT(IN)    :: State_Met    ! Meteorology State object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(DgnState), INTENT(INOUT) :: State_Diag
!
! OUTPUT PARAMETERS:
!
    REAL(fp),       INTENT(OUT)   :: BrO_MOLEC_CM3  ! BrO conc [molec/cm3]
!
! !RETURN VALUE:
!
    REAL(fp)                      :: BR_MOLEC_CM3   ! Br conc [molec/cm3]
!
! !REVISION HISTORY:
!  06 Jul 2006 - C. Holmes   - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    REAL(fp)              :: BR_PPTV, BR_MBL, BR_FAC
    REAL(fp)              :: LOCALTIME, HOUR, BRO_PPTV, FPBL
    REAL(fp)              :: BOXVL
!
! !DEFINED PARAMETERS:
!
    ! Constant concentration of BrO assumed for the marine boundary
    ! layer, based on finding from the Holmes et al.(2009) box model
    ! Use 0.5pptv for 24-hr mean, based on 1pptv in box model daytime
    REAL(fp), PARAMETER   :: BRO_MBL=0.5e+0_fp ! pptv

    ! Constant concentration of BrO assumed for the Polar boundary layer
    ! during springtime bromine explosion. Based on DOAS observations
    ! at Alert (Simpson et al. ACP 2007). This concentration applies
    ! only under stable atmospheric conditions.
    !REAL(fp), PARAMETER   :: BRO_POLAR=10e+0_fp
    ! Reduce concentration of BrO in Arctic during depletion events to 5 pptv
    ! Ref low end of uncertainty range in Holmes et al. 2010
    ! Ref low [BrO] observed Neuman et al 2010
    !REAL(fp), PARAMETER   :: BRO_POLAR=5e+0_fp

    ! Make polar BrO variable rather than fixed, varying as a factor of
    ! solar radiation and temperatures. (jaf, 11/29/11)
    REAL(fp)           :: BRO_POLAR, BR_POLAR, O3_POLAR, SWRAD
    LOGICAL            :: IS_MOSTLY_ICE
    ! Assume 5 ppb O3 when BrO present. Pohler et al. (2010) shows 5 ppb O3
    ! can exist for 3 ppt < [BrO] < 40 ppt. (jaf, 11/29/11)
    REAL(fp), PARAMETER  :: O3_POLAR_PPB = 5e+0_fp

    ! Parameters for calculating Br/BrO photostationary state
    ! BrO J value, /s
    ! Rate coefficient BrO + NO -> Br + NO2, cm3/molec/s
    REAL(fp), PARAMETER   :: K_BRO_NO = 2.1e-11_fp
    ! Rate coefficient Br + O3 -> BrO + O2, cm3/molec/s
    REAL(fp), PARAMETER   :: K_BR_O3  = 1.2e-12_fp
    ! Concentration of NO, based on 10pptv, molec/cm3
    REAL(fp), PARAMETER   :: C_NO     = 2.5e+8_fp

    !=================================================================
    ! GET_BR begins here!
    !=================================================================

    !----------------------------------------------------------------
    ! Get monthly-mean Br from bromocarbon and halon sources
    !----------------------------------------------------------------

    ! Zero polar BrOx, O3 concentrations (jaf, 11/29/11)
    BRO_POLAR     = 0e+0_fp
    BR_POLAR      = 0e+0_fp
    O3_POLAR      = 0e+0_fp
    SWRAD         = 0e+0_fp

    BR_PPTV  = BR_MERGE(I,J,L)
    BRO_PPTV = BRO_MERGE(I,J,L)

    ! Multiple stratospheric concentrations by Factor, if necessary
    IF ( State_Met%InStratosphere(I,J,L) .and. &
       ( ABS(STRAT_BR_FACTOR - 1e+0_fp) > 1e-2_fp) )    THEN

       BR_PPTV  = BR_PPTV  * STRAT_BR_FACTOR
       BRO_PPTV = BRO_PPTV * STRAT_BR_FACTOR

    ENDIF

    !----------------------------------------------------------------
    ! Add Br in the MBL
    !----------------------------------------------------------------
    IF ( (L_ADD_MBL_BR) .OR. (LPOLARBR) ) THEN
       ! Fraction of grid box mass which is within the mixed layer
       FPBL = State_Met%F_UNDER_PBLTOP(I,J,L)
    ENDIF

    IF ( L_ADD_MBL_BR ) THEN

       IF ( (State_Met%IsWater(I,J)) .AND. (FPBL > 0e+0_fp) ) THEN

          ! Convert BrO concentration to corresponding Br based on
          ! photochemical steady state, pptv
          ! Platt & Janssen, Faraday Discussions (1995)
          BR_MBL = BRO_MBL *( J_BRO(I,J,L) + K_BRO_NO * C_NO ) / &
                            ( K_BR_O3      * GET_O3(I,J,L,State_Met) ) * FPBL

          ! MBL concentration is the greater of the TOMCAT Br or 1pptv BrO
          BR_PPTV  = MAX( BR_PPTV,  BR_MBL )
          BrO_PPTV = MAX( BRO_PPTV, BRO_MBL*FPBL )
       ENDIF

    ENDIF

    ! Now moved to after imposed diurnal cycle, updated parameters (J. Fisher)
    !----------------------------------------------------------------
    ! OLD VERSION: Add Br in polar boundary layer
    !
    ! Bromine in the polar boundary layer requires cold temperatures
    ! (below 268K), sunlight and contact between open water and sea ice
    ! Assuming these conditions are met, we assume a fixed 5-10pptv BrO
    ! At the moment, we only have Br north of 60N, i.e. not Hudson's Bay
    !----------------------------------------------------------------
    !
    !IF ( ( LPOLARBR ) .AND. ( FPBL >  0e+0_fp   ) .AND. &
    !     ( State_Met%SUNCOS(I,J)   >  0e+0_fp   ) .AND. &
    !     ( State_Met%T(I,J,1)      <= 268e+0_fp ) .AND. &
    !     ( State_Met%LWI(I,J)      == 2         ) .AND. &
    !       GET_OBK(I,J,State_Met)  >  0e+0_fp   ) THEN
    !
    !   IF ( ( (State_Grid%YMid(I,J) > 60e+0_fp) .AND.             &
    !          (GET_MONTH() >= 3) .AND. (GET_MONTH() <= 5)  ) .OR. &
    !        ( (State_Grid%YMid(I,J) < -50e+0_fp) .AND.            &
    !          (GET_MONTH() >= 9) .AND. (GET_MONTH() <= 11)  ) ) THEN
    !
    !      ! Br concentration due to bromine explosion, pptv
    !      ! Assume [O3] is 2ppb during event
    !      BR_MBL = BRO_POLAR *( J_BRO(I,J,L) + K_BRO_NO * C_NO ) / &
    !                          ( K_BR_O3 * 5e+10_fp ) * FPBL
    !
    !      BR_PPTV  = BR_PPTV  + BR_MBL
    !      BrO_PPTV = BrO_PPTV + BRO_POLAR * FPBL
    !
    !   ENDIF
    !
    !ENDIF

    !----------------------------------------------------------------
    ! Impose a diurnal cycle
    !
    ! This formula comes from the Holmes et al (2009) MBL box model
    ! which fit RGM observations for Okinawa and a Pacific cruise
    !----------------------------------------------------------------

    ! Test for sunlight...
    IF ( (State_Met%SUNCOS(I,J) > 0e+0_fp) .and. (TTDAY(I,J) > 0e+0_fp) ) THEN

       ! Use a constant function if daylight is < 2 hours or > 22hours
       IF ( ( TTDAY(I,J) < 120e+0_fp  ) .OR. ( TTDAY(I,J) > 1320e+0_fp ) ) THEN

          BR_FAC = SAFE_DIV( 1440e+0_fp, TTDAY(I,J), 0e+0_fp )

       ELSE

          ! Local time: 0-23.999
          LOCALTIME = GET_LOCALTIME( I, J, L, State_Grid )

          ! Interpolate the real hour to lie within an equinoctal day
          ! i.e. between 6-18
          HOUR = 12e+0_fp + ( LOCALTIME - 12e+0_fp ) * 720e+0_fp / ( TTDAY(I,J))

          ! Multiplicative factor for the diurnal cycle
          BR_FAC = (1e+6_fp - 1e+2_fp * (HOUR - 6e+0_fp) ** 4e+0_fp + &
                   ( 18e+0_fp - HOUR ) * 1e+6_fp / 12e+0_fp ) / 2e+6_fp

          ! Normalize the multiplicative factor to have a 24-h mean of 1
          BR_FAC = SAFE_DIV( BR_FAC, (4e-4_fp * TTDAY(I,J)), 0e+0_fp )

       ENDIF

       ! Make sure that diurnal scaling factor is non-negative
       BR_FAC = MAX( BR_FAC, 0e+0_fp )

       ! The instantaneous concentration is the 24-h mean times
       ! the time-varying factor from the diurnal cycle
       BR_PPTV = BR_PPTV * BR_FAC

    ELSE

       ! At night, BR goes to zero
       BR_PPTV = 0e+0_fp

    ENDIF

    ! Grid box volume [cm3]
    BOXVL  = State_Met%AIRVOL(I,J,L) * 1e+6_fp

    ! Do this AFTER application of the diurnal cycle, since polar
    ! bromine has implicit diurnal cycle driven by solar radiation
    ! and represents instantaneous concentrations (jaf, 12/16/11)
    !----------------------------------------------------------------
    ! Add Br in polar boundary layer (jaf, 11/29/11)
    !
    ! Bromine in the polar boundary layer requires the following
    ! criteria be met:
    ! - In the PBL
    ! - Downward shortwave radiation > 100 W/m2 (Pohler et al. 2010)
    ! - Sea ice exists (>50% native boxes have >10% ice cover)
    ! - Breaks in sea ice exist (<100% native boxes have >90% ice cover)
    ! - Month is between Feb & June (Arctic) or Aug & Dec (Antarctic)
    !   based on http://bro.aeronomie.be/level3_monthly.php?cmd=map
    ! - Temperature is less than 0C
    !
    ! If these criteria are met, BrO is a function of ambient temp.
    ! with [BrO] based on findings from Pohler et al. (2010) and
    ! Prados-Roman et al. (2011). O3 used to convert BrO to Br is 5
    ! ppb, based on data from Pohler et al. (2010).
    !----------------------------------------------------------------
    SWRAD         = State_Met%SWGDN(I,J)
    IS_MOSTLY_ICE = ( State_Met%SEAICE00(I,J) <= 0.5e+0_fp .AND. &
                      State_Met%SEAICE90(I,J) < 1e+0_fp )

    IF ( (LPOLARBR) .AND. (FPBL > 0e+0_fp) .AND. (IS_MOSTLY_ICE) .AND. &
         (SWRAD > 1e+2_fp) .AND. (State_Met%TS(I,J) <= 273e+0_fp) ) THEN

       IF ( ((State_Grid%YMid(I,J) > 0e+0_fp) .AND.           &
            (GET_MONTH() >= 2) .AND. (GET_MONTH() <= 6)) .OR. &
            ((State_Grid%YMid(I,J) < 0e+0_fp) .AND.           &
            (GET_MONTH() >= 8) .AND. (GET_MONTH() <= 12)) ) THEN

          ! [BrO] is a linear function of temperature derived based on
          ! results from Pohler et al. (2010), Prados-Roman et al. (2011)
          ! and ability to match Hg0 seasonal cycle at Alert. (jaf,
          ! 12/24/11)
          IF ( State_Met%TS(I,J) <= 253e+0_fp ) THEN
             BRO_POLAR = 2e+1_fp
          ELSE
             BRO_POLAR = -1e+0_fp * ( State_Met%TS(I,J) - 253e+0_fp ) + 2e+1_fp
          ENDIF

          ! Convert O3 to molec/cm3
          O3_POLAR = O3_POLAR_PPB * 1e-9_fp * ( AVO * 1e+3_fp / &
                     AIRMW ) * State_Met%AD(I,J,L) / BOXVL

          ! Compute polar Br, BrO concentrations in pptv
          BrO_POLAR = BrO_POLAR * FPBL
          Br_POLAR  = BRO_POLAR * ( J_BRO(I,J,L) + K_BRO_NO * C_NO ) / &
                                  ( K_BR_O3 * O3_POLAR )

          BR_PPTV  = BR_PPTV  + Br_POLAR
          BrO_PPTV = BrO_PPTV + BrO_POLAR

#ifdef BPCH_DIAG
          !==============================================================
          ! %%%%% ND03 (bpch) DIAGNOSTIC %%%%%
          !
          ! Polar Br [pptv]; Polar Br [pptv]; Polar O3 [ppbv]
          !==============================================================
          IF ( ND03 > 0 .AND. L <= LD03 ) THEN
             AD03_Br(I,J,L,3) = AD03_Br(I,J,L,3) + BR_POLAR
             AD03_Br(I,J,L,4) = AD03_Br(I,J,L,4) + BRO_POLAR
             AD03_Br(I,J,L,5) = AD03_Br(I,J,L,5) + O3_POLAR_PPB
          ENDIF
#endif

          !==============================================================
          ! %%%%% HISTORY (aka netCDF diagnostics) %%%%%
          !
          ! Polar Br [pptv]; Polar Br [pptv]; Polar O3 [ppbv]
          !==============================================================
          IF ( State_Diag%Archive_PolarConcBr ) THEN
             State_Diag%PolarConcBr(I,J,L) = Br_Polar
          ENDIF

          IF ( State_Diag%Archive_PolarConcBrO ) THEN
             State_Diag%PolarConcBr(I,J,L) = BrO_Polar
          ENDIF

          IF ( State_Diag%Archive_PolarConcO3 ) THEN
             State_Diag%PolarConcO3(I,J,L) = O3_Polar_PPB
          ENDIF
       ENDIF ! Month

    ENDIF ! Polar BrO criteria
    !----------------------------------------------------------------

    ! Convert pptv mixing ratio -> molec/cm3
    BR_MOLEC_CM3  = BR_PPTV  * 1e-12_fp * ( AVO * 1e+3_fp / AIRMW ) * &
                    State_Met%AD(I,J,L) / BOXVL
    BRO_MOLEC_CM3 = BRO_PPTV * 1e-12_fp * ( AVO * 1e+3_fp / AIRMW ) * &
                    State_Met%AD(I,J,L) / BOXVL

  END FUNCTION GET_BR
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ohno3time
!
! !DESCRIPTION: Subroutine OHNO3TIME computes the sum of cosine of the
!  solar zenith angle over a 24 hour day, as well as the total length of
!  daylight.  This is needed to scale the offline OH and NO3 concentrations.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE OHNO3TIME( State_Grid )
!
! !USES:
!
    USE State_Grid_Mod, ONLY : GrdState
    USE TIME_MOD,       ONLY : GET_NHMSb,   GET_ELAPSED_SEC
    USE TIME_MOD,       ONLY : GET_TS_CHEM, GET_DAY_OF_YEAR, GET_GMT
!
! !INPUT PARAMETERS:
!
    TYPE(GrdState), INTENT(IN) :: State_Grid  ! Grid State object
!
! !REVISION HISTORY:
!  16 Dec 2002 - R. Park & R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    LOGICAL, SAVE :: FIRST = .TRUE.
    INTEGER       :: I, J, L, N, NT, NDYSTEP
    REAL(fp)      :: A0, A1, A2, A3, B1, B2, B3
    REAL(fp)      :: LHR0, R, AHR, DEC, TIMLOC, YMID_R
    REAL(fp)      :: SUNTMP(State_Grid%NX,State_Grid%NY)

    !=================================================================
    ! OHNO3TIME begins here!
    !=================================================================

    !  Solar declination angle (low precision formula, good enough for us):
    A0 = 0.006918
    A1 = 0.399912
    A2 = 0.006758
    A3 = 0.002697
    B1 = 0.070257
    B2 = 0.000907
    B3 = 0.000148
    R  = 2.* PI * float( GET_DAY_OF_YEAR() - 1 ) / 365.

    DEC = A0 - A1*cos(  R) + B1*sin(  R) &
             - A2*cos(2*R) + B2*sin(2*R) &
             - A3*cos(3*R) + B3*sin(3*R)

    LHR0 = int(float( GET_NHMSb() )/10000.)

    ! Only do the following at the start of a new day
    IF ( FIRST .or. GET_GMT() < 1e-5 ) THEN

       ! Zero arrays
       TTDAY(:,:) = 0e+0_fp
       TCOSZ(:,:) = 0e+0_fp
       COSZM(:,:) = 0e+0_fp

       ! NDYSTEP is # of chemistry time steps in this day
       NDYSTEP = ( 24 - INT( GET_GMT() ) ) * 3600 / GET_TS_CHEM()

       ! NT is the elapsed time [s] since the beginning of the run
       NT = GET_ELAPSED_SEC()

       ! Loop forward through NDYSTEP "fake" timesteps for this day
       DO N = 1, NDYSTEP

          ! Zero SUNTMP array
          SUNTMP = 0e+0_fp

          ! Loop over surface grid boxes
          !$OMP PARALLEL DO       &
          !$OMP DEFAULT( SHARED ) &
          !$OMP PRIVATE( I, J, YMID_R, TIMLOC, AHR )
          DO J = 1, State_Grid%NY
          DO I = 1, State_Grid%NX

             ! Grid box latitude center [radians]
             YMID_R = State_Grid%YMid_R( I, J )

             TIMLOC = real(LHR0) + real(NT)/3600.0 + State_Grid%XMid(I,J)/15.0

             DO WHILE (TIMLOC < 0)
                TIMLOC = TIMLOC + 24.0
             ENDDO

             DO WHILE (TIMLOC > 24.0)
                TIMLOC = TIMLOC - 24.0
             ENDDO

             AHR = abs(TIMLOC - 12.) * 15.0 * PI_180

             !===========================================================
             ! The cosine of the solar zenith angle (SZA) is given by:
             !
             !  cos(SZA) = sin(LAT)*sin(DEC) + cos(LAT)*cos(DEC)*cos(AHR)
             !
             ! where LAT = the latitude angle,
             !       DEC = the solar declination angle,
             !       AHR = the hour angle, all in radians.
             !
             ! If SUNCOS < 0, then the sun is below the horizon, and
             ! therefore does not contribute to any solar heating.
             !===========================================================

             ! Compute Cos(SZA)
             SUNTMP(I,J) = sin(YMID_R) * sin(DEC) + &
                           cos(YMID_R) * cos(DEC) * cos(AHR)

             ! TCOSZ is the sum of SUNTMP at location (I,J)
             ! Do not include negative values of SUNTMP
             TCOSZ(I,J) = TCOSZ(I,J) + MAX( SUNTMP(I,J), 0e+0_fp )

             ! COSZM is the peak value of SUMTMP during a day at (I,J)
             ! (rjp, bmy, 3/30/04)
             COSZM(I,J) = MAX( COSZM(I,J), SUNTMP(I,J) )

             ! TTDAY is the total daylight time at location (I,J)
             IF ( SUNTMP(I,J) > 0e+0_fp ) THEN
                TTDAY(I,J) = TTDAY(I,J) + DBLE( GET_TS_CHEM() ) / 60e+0_fp
             ENDIF
          ENDDO
          ENDDO
          !$OMP END PARALLEL DO

          ! Increment elapsed time [sec]
          NT = NT + GET_TS_CHEM()
       ENDDO

       ! Reset first-time flag
       FIRST = .FALSE.
    ENDIF

  END SUBROUTINE OHNO3TIME
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: calc_hg2_seasalt_lossrate
!
! !DESCRIPTION: Subroutine CALC\_HG2\_SEASALT\_LOSSRATE calculates the loss
!  rate of RGM (/s) by uptake of RGM into sea salt aerosol for each model
!  grid. Return value is a loss frequency (/s)
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE CALC_HG2_SEASALT_LOSSRATE( State_Grid, State_Met )
!
! !USES:
!
    USE State_Grid_Mod,     ONLY : GrdState
    USE State_Met_Mod,      ONLY : MetState
!
! !INPUT PARAMETERS:
!
    TYPE(GrdState), INTENT(IN)  :: State_Grid  ! Grid State object
    TYPE(MetState), INTENT(IN)  :: State_Met   ! Meteorology State object
!
! !REMARKS:
!  The formula used here is a least-squares fit to the full-physics model of
!  sea-salt aerosol emissions, hydroscopic growth, mass-transport limited
!  uptake of Hg(II), and aerosol deposition presented by Holmes et al. (2009)
!  See Holmes et al. 2010 for evaluation of this parameterization.
!  (cdh, 11/25/09)
!
! !REVISION HISTORY:
!  25 Nov 2009 - C. Holmes   - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    REAL(fp)      :: U10M, S
    REAL(fp)      :: LOSS_FREQ
    INTEGER       :: I, J
    REAL(fp),SAVE :: TABLE_S(21), TABLE_U10(20)
    REAL(fp),SAVE :: TABLE_UPTAKE(21,20)
    REAL(fp)      :: SFCWINDSQR

    ! Flag for first call
    LOGICAL, SAVE :: FIRST=.TRUE.

    !=================================================================
    ! HG2_SEASALT_LOSSRATE begins here!
    !=================================================================

    !$OMP PARALLEL DO       &
    !$OMP DEFAULT( SHARED ) &
    !$OMP PRIVATE( I, J, U10M, S, LOSS_FREQ, SFCWINDSQR )
    DO J = 1, State_Grid%NY
    DO I = 1, State_Grid%NX

       ! Only calculate deposition via sea salt over water
       IF ( State_Met%IsWater(I,J) ) THEN

          ! Wind speed at 10m altitude [m/s]
          SFCWINDSQR = State_Met%U10M(I,J)**2 + State_Met%V10M(I,J)**2
          U10M       = SQRT( SFCWINDSQR )

          ! Don't allow wind >20 m/s which is the limit of this parameterization
          U10M       = MAX( MIN( U10M, 20e+0_fp ), 1e+0_fp )

          ! Relative humidity as a saturation ratio
          ! Use the relative humidity of the lowest layer, although this is
          ! lower than the higher layers of the MBL
          !
          ! Don't allow supersaturation, as [Cl-] is undefined for RH>=1
          ! Cap RH at 99%, Don't allow RH < 75% as happens in coastal areas
          S = MAX( MIN( State_Met%RH(I,J,1), 99e+0_fp ), 75e+0_fp ) * 1e-2_fp

          LOSS_FREQ = 1e-10_fp * ( 1e+0_fp - EXP( -57.758e+0_fp * &
                      (1e+0_fp-S) ) ) * &
                      EXP( -1.9351e+0_fp  * U10M + &
                            9.0047e+0_fp  * SQRT( U10M ) + &
                            0.14788e+0_fp * U10M**1.5e+0_fp )

          ! Loss frequency must be positive
          LOSS_FREQ = MAX( LOSS_FREQ, 1e-10_fp )

       ELSE

          ! No loss over land
          LOSS_FREQ = 0e+0_fp

       ENDIF

       HG2_SEASALT_LOSSRATE(I,J) = LOSS_FREQ

    ENDDO
    ENDDO
    !$OMP END PARALLEL DO

  END SUBROUTINE CALC_HG2_SEASALT_LOSSRATE
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_jno2
!
! !DESCRIPTION: Function GET\_JNO2 returns monthly mean JNO2 and imposes a
!  diurnal variation.
!\\
!\\
! !INTERFACE:
!
  FUNCTION GET_JNO2( I, J, L, State_Grid, State_Met ) &
       RESULT( JNO2_NOW )
!
! !USES:
!
    USE PhysConstants
    USE State_Grid_Mod,     ONLY : GrdState
    USE State_Met_Mod,      ONLY : MetState
    USE TIME_MOD,           ONLY : GET_DAY_OF_YEAR
!
! !INPUT PARAMETERS:
!
    INTEGER,        INTENT(IN) :: I           ! Longitude index
    INTEGER,        INTENT(IN) :: J           ! Latitude index
    INTEGER,        INTENT(IN) :: L           ! Level index
    TYPE(GrdState), INTENT(IN) :: State_Grid  ! Grid State object
    TYPE(MetState), INTENT(IN) :: State_Met   ! Meteorology State object
!
! !RETURN VALUE:
!
    REAL(fp)                     :: JNO2_NOW    ! J(NO2) value [s-1]
!
! !REMARKS:
!  Impose the diurnal variation of JNO2 found by Parrish et al. (1983) under
!  clear skies. J-NO2 ~ exp( -0.360 * sec(SZA) )
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER  :: JDAY
    REAL(fp) :: CSZANOON
    REAL(fp) :: A0, A1, A2, A3, B1, B2, B3
    REAL(fp) :: R, DEC, YMID_R

    !=================================================================
    ! GET_JNO2 begins here!
    !=================================================================

    ! Test for sunlight...
    IF ( State_Met%SUNCOS(I,J) > 0e0_fp .and. TCOSZ(I,J) > 0e0_fp ) THEN

       ! Day of year
       JDAY = GET_DAY_OF_YEAR()

       ! Coefficients for solar declination angle
       A0  = 0.006918e+0_fp
       A1  = 0.399912e+0_fp
       A2  = 0.006758e+0_fp
       A3  = 0.002697e+0_fp
       B1  = 0.070257e+0_fp
       B2  = 0.000907e+0_fp
       B3  = 0.000148e+0_fp

       ! Path length of earth's orbit traversed since Jan 1 [radians]
       R   = ( 2e+0_fp * PI / 365e+0_fp ) * DBLE( JDAY - 1 )

       ! Solar declination angle (low precision formula)
       DEC = A0 - A1*COS(     R ) + B1*SIN(     R ) &
                - A2*COS( 2e+0_fp*R ) + B2*SIN( 2e+0_fp*R ) &
                - A3*COS( 3e+0_fp*R ) + B3*SIN( 3e+0_fp*R )

       ! Latitude of grid box [radians]
       YMID_R = State_Grid%YMid_R( I, J )

       ! Cosine of solar zenith angle at local noon
       CSZANOON = SIN( YMID_R ) * SIN( DEC ) + COS( YMID_R ) * COS( DEC )

       ! Parrish et al (1983) found
       ! J-NO2 ~ exp( -0.360 * sec(SZA) ), so
       ! J-NO2(now) = J-NO2(noon) *  exp( -0.360 * sec(SZANOW)  ) /
       !                             exp( -0.360 * sec(SZANOON) )
       !            = J-NO2(noon) * exp( 0.360 * [sec(SZANOON)-sec(SZANOW)] )

       ! Impose a diurnal variation on JNO2 during the day
       ! Note: We don't need to check for divide-by-zero errors
       ! because we already checked SUNCOS(I,J) and we know
       ! CSZANOON >= SUNCOS(I,J)
       JNO2_NOW = JNO2(I,J,L) * &
        EXP( 0.36e+0_fp * ( 1e+0_fp/CSZANOON - 1e+0_fp/State_Met%SUNCOS(I,J) ) )

       ! Make sure OH is not negative
       JNO2_NOW= MAX( JNO2_NOW, 0e+0_fp )

    ELSE

       ! At night, JNO2 goes to zero
       JNO2_NOW = 0e+0_fp

    ENDIF

  END FUNCTION GET_JNO2
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: partitionhg2
!
! !DESCRIPTION: Subroutine PARTITIONHG2 splits Hg(II) into gas and aerosol
!  portions  according to the thermodynamic equilibrium determined by
!  temperature and aerosol surface area.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE PARTITIONHG2( Input_Opt, State_Chm, State_Grid, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE ERROR_MOD,          ONLY : SAFE_DIV  , GEOS_CHEM_STOP
    USE Input_Opt_Mod,      ONLY : OptInput
    USE OCEAN_MERCURY_MOD,  ONLY : Fg
    USE State_Chm_Mod,      ONLY : ChmState
    USE State_Grid_Mod,     ONLY : GrdState
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
!  12-Jul-2010 - H. Amos     - Add option to partition Hg2 according to Fg/Fp
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER           :: I, J, L, N
    REAL(fp)          :: FGas, Hg2TOT

    ! Pointers
    REAL(fp), POINTER :: Spc(:,:,:,:)

    !=================================================================
    ! PARTITIONHG2 begins here!
    !=================================================================

    ! Assume success
    RC  =  GC_SUCCESS

    ! Point to the chemical species array [kg]
    Spc => State_Chm%Species

    !$OMP PARALLEL DO       &
    !$OMP DEFAULT( SHARED ) &
    !$OMP PRIVATE( I, J, L, N, FGas, Hg2TOT )
    DO L=1, State_Grid%NZ
    DO J=1, State_Grid%NY
    DO I=1, State_Grid%NX

       ! Fraction of Hg(II) as gas
       IF ( LHG2HALFAEROSOL ) THEN
          ! Assume constant partitioning
          Fgas = 0.5
       ELSE
          ! Temperature and aerosol-dependent partitioning
          Fgas = Fg(I,J,L)
       ENDIF

       ! Loop over the Hg regional tags
       DO N = 1, N_HG_CATS

          ! Total Hg(II) (gas +aerosol)
          HG2TOT = Spc(I,J,L,Hg2_Id_List(N)) + Spc(I,J,L,HgP_Id_List(N))

          ! Gas portion
          Spc(I,J,L,Hg2_Id_List(N)) = Hg2TOT * Fgas

          ! Aerosol portion
          Spc(I,J,L,HgP_Id_List(N)) = Hg2TOT * (1e+0_fp - Fgas)

       ENDDO

    ENDDO
    ENDDO
    ENDDO
    !$OMP END PARALLEL DO

    ! Free pointer memory
    Spc => NULL()

  END SUBROUTINE PARTITIONHG2
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: red_inplume_grid
!
! !DESCRIPTION: Subroutine RED\_INPLUME\_GRID conducts in plume reduction of
!  Hg2 for selected grids
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE RED_INPLUME_GRID( I, J, E_plant )
!
! !INPUT PARAMETERS:
!
    INTEGER,  INTENT(IN) :: I, J
    REAL(fp), INTENT(IN) :: E_plant
!
! !REVISION HISTORY:
!  11 Jan 2011 - Y. Zhang - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    REAL(fp) :: E_deg
!
! !DEFINED PARAMETERS:
!
    ! The source species profile (0.2 Hg0)
    ! 0.78 Hg2, 0.02 HgP, Streets, 2005), we use 75%, very high end
    REAL(fp),  PARAMETER :: K_Red_InPlume = 7.5e-1_fp

    !=================================================================
    ! RED_INPLUME_GRID begins here!
    !=================================================================

    ! Calculate the mass of Hg2 been degraded in plume
    E_deg = K_Red_InPlume * E_plant

    ! Subtract this part of Hg2 from the emission
    EHg2_an(I,J) = EHg2_an(I,J) - E_deg

    ! Degraded to Hg0
    EHg0_an(I,J) = EHg0_an(I,J) + E_deg

  END SUBROUTINE RED_INPLUME_GRID
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: do_red_inplume
!
! !DESCRIPTION: Subroutine DO\_RED\_INPLUME conducts in plume reduction of
!  Hg2 for selected grids.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE DO_RED_INPLUME( Input_Opt, State_Grid, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE HCO_INTERFACE_MOD,  ONLY : HcoState
    USE HCO_EmisList_Mod,   ONLY : HCO_GetPtr
    USE Input_Opt_Mod,      ONLY : OptInput
    USE State_Grid_Mod,     ONLY : GrdState
    USE TIME_MOD,           ONLY : EXPAND_DATE

    IMPLICIT NONE
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)  :: Input_Opt   ! Input Options object
    TYPE(GrdState), INTENT(IN)  :: State_Grid  ! Grid State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT) :: RC          ! Success or failure?
!
! !REVISION HISTORY:
!  11 Jan 2011 - Y. Zhang    - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER             :: I, J
    REAL(fp)            :: E_plant

    ! Strings
    CHARACTER(LEN=255)  :: ThisLoc
    CHARACTER(LEN=512)  :: ErrMsg

    ! Pointers
    REAL(f4), POINTER   :: E_ELEC_Hg2(:,:)
!
! !DEFINED PARAMETERS:
!
    REAL(fp), PARAMETER :: SEC_PER_YR = 365.25e+0_fp * 86400e+0_fp

    !=================================================================
    ! DO_RED_INPLUME begins here!
    !=================================================================

    ! Initialize
    RC         =  GC_SUCCESS
    ErrMsg     =  ''
    ThisLoc    =  ' -> at DO_RED_INPLUME (in GeosCore/mercury_mod.F90)'
    E_ELEC_Hg2 => NULL()

    ! Get a pointer to the monthly mean OH from HEMCO (bmy, 3/11/15)
    CALL HCO_GetPtr( HcoState, 'CFPP_NEI2005_Hg2', E_ELEC_Hg2, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Cannot get pointer to HEMCO field CFPP_NEI2005_Hg2!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    ! Do the reduction
    DO J = 1, State_Grid%NY
    DO I = 1, State_Grid%NX

       ! Hg emission in CFPP sector decreased for 45.14% during 2005-2010
       ! Convert [kg/m2/s] --> [kg/s]
       E_plant = E_ELEC_Hg2(I,J) * State_Grid%Area_M2(I,J)

       ! Reduce the Hg2 from plume
       CALL RED_INPLUME_GRID( I, J, E_plant )

    ENDDO
    ENDDO

    ! Free npointer
    E_ELEC_Hg2 => NULL()

  END SUBROUTINE DO_RED_INPLUME
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: offlineocean_readmo
!
! !DESCRIPTION: Subroutine OFFLINEOCEAN\_READMO reads the monthly varying
!     offline ocean evasion emissions if LDYNOCEAN is FALSE. Will not actually
!     need mixed layer depth when i get stuff from Yanxu
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE OFFLINEOCEAN_READMO( State_Chm, State_Diag, State_Grid, &
                                  State_Met, FLUX, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE TIME_MOD,          ONLY : EXPAND_DATE,GET_YEAR, GET_TS_EMIS
    USE TIME_MOD,          ONLY : ITS_A_NEW_MONTH, GET_MONTH
#ifdef BPCH_DIAG
    USE DIAG03_MOD,        ONLY : AD03, ND03
#endif
    USE State_Chm_Mod,     ONLY : ChmState
    USE State_Diag_Mod,    ONLY : DgnState
    USE State_Grid_Mod,    ONLY : GrdState
    USE State_Met_Mod,     ONLY : MetState
    USE HCO_INTERFACE_MOD, ONLY : HcoState
    USE HCO_EmisList_Mod,  ONLY : HCO_GetPtr
!
! !INPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(IN)    :: State_Chm    ! Chemistry State object
    TYPE(GrdState), INTENT(IN)    :: State_Grid   ! Grid State object
    TYPE(MetState), INTENT(IN)    :: State_Met    ! Meteorology State object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(DgnState), INTENT(INOUT) :: State_Diag   ! Diagnostics State Object
!
! !OUTPUT PARAMETERS:
!
    REAL*8,         INTENT(OUT)   :: FLUX(State_Grid%NX,State_Grid%NY,1)
    INTEGER,        INTENT(OUT)   :: RC           ! Success or failure?
!
! !REVISION HISTORY:
!  12 Aug 2015 - H. Horowitz - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER                 :: I, J, M, DN(12) ! M is mon., DN is days in mon
    INTEGER                 :: NNN
    INTEGER                 :: THISYEAR
    INTEGER                 :: THISMONTH

    LOGICAL, SAVE           :: FIRST = .TRUE.
    LOGICAL                 :: IS_OCEAN_BOX
    INTEGER                 :: NN, N
    REAL(fp)                :: A_M2,     DTSRCE
    REAL(fp)                :: CHg0aq,   CHg0,     vi, Hg0aqtemp
    REAL(fp)                :: TC,       TK,       Kw
    REAL(fp)                :: Sc,       ScCO2,    USQ
    REAL(fp)                :: FRAC_L,   FRAC_O,   H, D
    REAL(fp)                :: FUP(State_Grid%NX,State_Grid%NY,N_Hg_CATS)
    REAL(fp)                :: FDOWN(State_Grid%NX,State_Grid%NY,N_Hg_CATS)
    REAL(fp)                :: Hg0aq(State_Grid%NX,State_Grid%NY,N_Hg_CATS)
    REAL(fp)                :: MHg0_air

    ! Conversion factor from [cm/h * ng/L] --> [kg/m2/s]
    REAL(fp),  PARAMETER    :: TO_KGM2S = 1.0e-11_fp / 3600.0e+0_fp

    ! Small numbers to avoid dividing by zero
    REAL(fp),  PARAMETER    :: SMALLNUM = 1.0e-32_fp

    REAL(fp)                :: SFCWINDSQR

    ! Pointers
    ! We need to define local arrays to hold corresponding values
    ! from the Chemistry State (State_Chm) object. (mpayer, 12/6/12)
    REAL(fp), POINTER       :: STT(:,:,:,:)

    ! Characters
    CHARACTER(LEN=255)      :: ThisLoc
    CHARACTER(LEN=512)      :: ErrMsg

    !=================================================================
    ! OFFLINEOCEAN_READMO begins here!
    !=================================================================

    ! Initialize
    RC      = GC_SUCCESS
    ErrMsg  = ''
    ThisLoc = ' -> at OFFLINEOCEAN_READMO (in GeosCore/mercury_mod.F90)'

    ! Get month
    THISMONTH = GET_MONTH()
    M         = THISMONTH

    ! Days in each month (will use later) 9/16/15 hmh
    DN =  (/ 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /)

    !-----------------------------------------------------------------
    ! %%%%% HISTORY (aka netCDF diagnostics) %%%%%
    !
    ! Zero flux arrays of State_Diag to prevent leftover values
    ! from the last timestep from being included in the averaging
    !-----------------------------------------------------------------
    IF ( State_Diag%Archive_FluxHg0fromOceanToAir ) THEN
       State_Diag%FluxHg0fromOceanToAir = 0.0_f4
    ENDIF

    IF ( State_Diag%Archive_FluxHg0fromAirToOcean ) THEN
       State_Diag%FluxHg0fromAirToOcean = 0.0_f4
    ENDIF

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! read monthly ocean evasion  !
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    IF ( ITS_A_NEW_MONTH() ) THEN

       CALL HCO_GetPtr( HcoState, 'GLOBAL_OCEAN', OCEAN_CONC, RC )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Cannot get pointer to HEMCO field GLOBAL_OCEAN!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

    ENDIF

    ! Only doing Hg0 overall, should add trap for LSPLIT (cpt - 2017)
    DO NNN = 1, N_Hg_CATS
       IF ( NNN .eq. 1 ) THEN
          Hg0aq(:,:,NNN) = OCEAN_CONC(:,:,1)
       ELSE
          Hg0aq(:,:,NNN) = 0.0e+0_fp
       ENDIF
    ENDDO

    ! Emission timestep [s]
    DTSRCE = GET_TS_EMIS()

    STT => State_Chm%Species

    ! Loop over surface boxes
    !$OMP PARALLEL DO       &
    !$OMP DEFAULT( SHARED ) &
    !$OMP PRIVATE( I,   A_M2,    vi,     ScCO2 )               &
    !$OMP PRIVATE( J,   NN,      TC,     TK )                  &
    !$OMP PRIVATE( N,   CHg0,    FRAC_L, FRAC_O  )             &
    !$OMP PRIVATE( H,   Kw,      CHg0aq, Hg0aqtemp, MHg0_air ) &
    !$OMP PRIVATE( IS_OCEAN_BOX, Sc,     Usq, D )              &
    !$OMP SCHEDULE( DYNAMIC )
    DO J = 1, State_Grid%NY
    DO I = 1, State_Grid%NX

       ! Grid box surface area [m2]
       A_M2       = State_Grid%Area_M2( I, J )

       ! Initialize values
       Kw         = 0e0_fp
       TK         = 0e0_fp
       TC         = 0e0_fp

       ! Get fractions of land and ocean in the grid box [unitless]
       ! Use fractional land type information in MERRA. Also make sure
       ! we do not use boxes that are mostly sea ice for consistency
       ! FROCEAN is a constant, so to get correct ocean fraction we
       ! need to subtract the sea ice fraction. Don't let the fraction
       ! be less than zero (jaf, 4/26/11)
       FRAC_L       = State_Met%FRLAND(I,J)
       FRAC_O       = MAX( State_Met%FROCEAN(I,J) - &
                           State_Met%FRSEAICE(I,J), 0e0_fp )
       IS_OCEAN_BOX = ( ( FRAC_O > 0e0_fp ) .and. &
                        ( State_Met%SEAICE00(I,J)  > 0.5e0_fp ) )

       IF ( (IS_OCEAN_BOX) ) THEN

          !--------------------------------------------------------------
          ! Sea surface temperature in both [K] and [C]
          !--------------------------------------------------------------
          ! where TSKIN is the temperature (K) at the ground/sea surface
          ! (Use as surrogate for SST, cap at freezing point)
          TK     = MAX( State_Met%TSKIN(I,J), 273.15e0_fp )
          TC     = TK - 273.15e0_fp

          !==============================================================
          ! Volatilisation of Hg0
          !==============================================================

          ! Henry's law constant (gas->liquid) [unitless] [L water/L air]
          ! (ref: Andersson et al. 2008)
          H      = EXP( ( -2404.3e0_fp / TK ) + 6.92e0_fp )

          ! Viscosity as a function of changing temperatures
          ! (ref: Loux 2001)
          ! The paper says the viscosity is given in cP but us really P
          ! and we therefor multiply with 100 to get cP.
          vi    = ( 10**( ( 1301.0e0_fp / ( 998.333e0_fp + 8.1855e0_fp &
                  * ( TC - 20.0e0_fp )+ 0.00585e0_fp &
                  * ( TC - 20.0e0_fp )**2 ) ) - 3.30233e0_fp ) ) * 100.0e0_fp

          ! Schmidt # for Hg [unitless]
          ! Sc = v/D = kinematic viscosity/diffusivity
          ! (ref: Poissant et al 2000; Wilke and Chang 1995)
          ! to correct for seawater D0 is decreased by 6% as suggested
          ! by Wanninkhof (1992)
          D = 7.4e-8_fp * sqrt( 2.26e0_fp * 18.0e0_fp ) * TK / &
              ( ( 14.8e0_fp**0.6e0_fp ) *vi )
          Sc   = ( 0.017e0_fp * EXP( -0.025e0_fp * TC ) ) / D

          ! Schmidt # of CO2 [unitless] for CO2 in seawater at 20 degrees C
          ! The value is set to a constant based on other ocean studies
          ! (Gardfeld et al. 2003, Rolfhus & Fitzgerald2004,Mason et al.2001)
          !
          ! Correction of the Schmidt # with temperature based on Poissant
          ! et al. (2000) (for freshwatersystems).
          ScCO2  = 644.7e0_fp + TC * ( -6.16e0_fp + TC * ( 0.11e0_fp))

          ! Square of surface (actually 10m) wind speed [m2/s2]
          Usq    = State_Met%U10M(I,J)**2 + State_Met%V10M(I,J)**2

          !------------------------------------------------------
          ! Parameterizations for calculating water side mass trasfer
          ! coefficient
          !------------------------------------------------------
          ! Mass transfer coefficient [cm/h], from Nightingale et al. 2000
          Kw     = ( 0.25e0_fp * Usq ) / SQRT( Sc / ScCO2 )

          ! Loop over all Hg categories
          DO NN = 1, N_Hg_CATS

             ! Hg0 tracer number (for STT)
             N = Hg0_Id_List(NN)

             !--------------------------------------------------------
             ! Calculate oceanic and gas-phase concentration of Hg(0)
             !--------------------------------------------------------

             ! Concentration of Hg(0) in the ocean [ng/L]
             ! now converting from Hg0aq in mol/m3 to ng/L
             CHg0aq = Hg0aq(I,J,NN) *200.59e0_fp * 1.0e9_fp / 1.0e3_fp

             ! Gas phase Hg(0) concentration: convert [kg] -> [ng/L]
             MHg0_air = STT(I,J,1,N)
             CHg0     = MHg0_air *  1.0e9_fp /State_Met%AIRVOL(I,J,1)

             !--------------------------------------------------------
             ! Compute flux of Hg(0) from the ocean to the air
             !--------------------------------------------------------

             ! Compute ocean flux of Hg0 [cm/h*ng/L]
             FLUX(I,J,NN)  = Kw * ( CHg0aq - ( CHg0 / H ) )

             !Extra diagnostic: compute flux up and flux down
             FUP(I,J,NN)   = ( Kw * CHg0aq )
             FDOWN(I,J,NN) = ( Kw * CHg0 / H )

             !--------------------------------------------------
             ! Convert [cm/h*ng/L] --> [kg/m2/s] --> [kg/s]
             ! Also account for ocean fraction of grid box
             FLUX(I,J,NN)  = FLUX(I,J,NN) * TO_KGM2S * A_M2 * FRAC_O

             ! hmh 5/11/16 reverting to old version and uncommenting here
             FUP(I,J,NN)  = FUP(I,J,NN) * TO_KGM2S * A_M2 * FRAC_O
             FDOWN(I,J,NN)  = FDOWN(I,J,NN) * TO_KGM2S * A_M2 * FRAC_O
             !--------------------------------------------------

             !--------------------------------------------------------
             ! Flux limited by ocean and atm Hg(0)
             !--------------------------------------------------------

             ! Cap the flux w/ the available Hg(0) ocean mass
             ! THIS IS THE PROBLEM!
             Hg0aqtemp = CHg0aq * A_M2 * FRAC_O *1.0e-8_fp

             IF ( FLUX(I,J,NN) * DTSRCE > Hg0aqtemp ) THEN
                FLUX(I,J,NN) = Hg0aqtemp / DTSRCE
                FUP(I,J,NN)  = FLUX(I,J,NN)-FDOWN(I,J,NN)
             ENDIF

             ! Cap the neg flux w/ the available Hg(0) atm mass
             IF ( (-FLUX(I,J,NN) * DTSRCE ) > MHg0_air ) THEN
                FLUX(I,J,NN) = -MHg0_air / DTSRCE
             ENDIF

             ! make sure Fup and Fdown do not underflow either
             ! debug 2x2.5 diagnostic?
             FUP(I,J,NN) = MAX(FUP(I,J,NN), SMALLNUM )
             FDOWN(I,J,NN) = MAX(FDOWN(I,J,NN), SMALLNUM )

#ifdef BPCH_DIAG
             !--------------------------------------------------------
             ! %%%%% ND03 (bpch) DIAGNOSTICS %%%%%
             !
             ! Fluxes of Hg0 from air to ocean and ocean to air [kg]
             !--------------------------------------------------------
             IF ( ND03 > 0 ) THEN
                AD03(I,J,16,NN) = AD03(I,J,16,NN) + FUP(I,J,NN) * DTSRCE
                AD03(I,J,17,NN) = AD03(I,J,17,NN) + FDOWN(I,J,NN) *DTSRCE
             ENDIF
#endif

             !--------------------------------------------------------
             ! %%%%% HISTORY (aka netCDF diagnostics) %%%%%
             !
             ! Fluxes of Hg0 from air to ocean and ocean to air [kg/s]
             ! NOTE: Implement for total Hg species at ths time
             !--------------------------------------------------------
             IF ( NN == 1 ) THEN

                ! Flux of Hg0 from ocean to air [kg/s]
                IF ( State_Diag%Archive_FluxHg0fromOceanToAir ) THEN
                   State_Diag%FluxHg0fromOceanToAir(I,J) = FUP(I,J,NN)
                ENDIF

                IF ( State_Diag%Archive_FluxHg0fromAirToOcean ) THEN
                   State_Diag%FluxHg0fromAirToOcean(I,J) = FDOWN(I,J,NN)
                ENDIF

             ENDIF

          ENDDO

       ELSE

          DO NN = 1, N_Hg_CATS
             FLUX(I,J,NN)  = 0e0_fp
             FUP(I,J,NN)   = 0e0_fp
             FDOWN(I,J,NN) = 0e0_fp
          ENDDO

       ENDIF
    ENDDO
    ENDDO
    !$OMP END PARALLEL DO

    ! Free pointer
    NULLIFY( STT )

  END SUBROUTINE OFFLINEOCEAN_READMO
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Reset_Hg_Diags
!
! !DESCRIPTION: Zeroes the relevant diagnostic fields of State_Diag for
!  the Hg and tagged Hg simulations.  Some of these need to be done
!  for example, at the start of each timestep.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Reset_Hg_Diags( Input_Opt, State_Diag, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE Input_Opt_Mod,  ONLY : OptInput
    USE State_Diag_Mod, ONLY : DgnState
    USE Time_Mod,       ONLY : Its_Time_For_Chem
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(DgnState), INTENT(INOUT) :: State_Diag  ! Diagnostics State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure
!
! !REMARKS:
!  NOTE: Not all netCDF diagnostic fields need to be zeroed.  Some fields
!  of State_Diag are always
!
! !REVISION HISTORY:
!  06 Jan 2015 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    CHARACTER(LEN=255) :: ErrMsg, ThisLoc

    ! Assume success
    RC      = GC_SUCCESS
    ErrMsg  = ''
    ThisLoc = ' -> at Reset_Hg_Diags (in module GeosCore/mercury_mod.F90'

    !=================================================================
    ! Reset_Hg_Diags begins here!
    !=================================================================

    ! Exit if it's not a Hg sim
    IF ( .not. Input_Opt%ITS_A_MERCURY_SIM ) THEN
       ErrMsg = 'Routine "Reset_Hg_Diags" was called, but this ' // &
                'routine is only valid for Hg or tagHg simulations!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    !--------------------------------------------------------------
    ! Zero diagnostics in depo_mercury_mod.F90
    !--------------------------------------------------------------
    IF ( State_Diag%Archive_FluxHg2HgPfromAirToSnow ) THEN
       State_Diag%FluxHg2HgPfromAirToSnow = 0.0_fp
    ENDIF

    ! These are called once per chemistry or emissions timestep
    IF ( Its_Time_For_Chem() ) THEN

       !--------------------------------------------------------------
       ! Zero diagnostics in mercury_mod.F90
       !--------------------------------------------------------------
       IF ( State_Diag%Archive_DryDepChm .or. State_Diag%Archive_DryDep ) THEN
          State_Diag%DryDepChm = 0.0_f4
       ENDIF

       !-------------------------------------------------------------
       ! Zero diagnostics in ocean_mercury_mod.F90
       !-------------------------------------------------------------
       IF ( State_Diag%Archive_EmisHg2rivers ) THEN
          State_Diag%EmisHg2rivers = 0.0_f4
       ENDIF

       IF ( State_Diag%Archive_EmisHg2snowToOcean ) THEN
          State_Diag%EmisHg2snowToOcean = 0.0_f4
       ENDIF

       IF ( State_Diag%Archive_FluxHg0fromAirToOcean ) THEN
          State_Diag%FluxHg0fromAirToOcean = 0.0_f4
       ENDIF

       IF ( State_Diag%Archive_FluxHg0fromOceantoAir ) THEN
          State_Diag%FluxHg0fromOceanToAir = 0.0_f4
       ENDIF

       IF ( State_Diag%Archive_FluxHg2HgPfromAirToOcean ) THEN
          State_Diag%FluxHg2HgPfromAirToOcean = 0.0_f4
       ENDIF

       IF ( State_Diag%Archive_FluxOCtoDeepOcean ) THEN
          State_Diag%FluxOCtoDeepOcean = 0.0_f4
       ENDIF

    ENDIF

  END SUBROUTINE Reset_Hg_Diags
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_mercury
!
! !DESCRIPTION: Subroutine INIT\_MERCURY allocates and zeroes all
!  module arrays.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE INIT_MERCURY( Input_Opt, State_Chm, State_Grid, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE Input_Opt_Mod,      ONLY : OptInput
    USE Species_Mod,        ONLY : Species
    USE State_Chm_Mod,      ONLY : ChmState
    USE State_Grid_Mod,     ONLY : GrdState
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)  :: Input_Opt   ! Input Options object
    TYPE(ChmState), INTENT(IN)  :: State_Chm   ! Chemistry State object
    TYPE(GrdState), INTENT(IN)  :: State_Grid  ! Grid State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT) :: RC          ! Success or failure?
!
! !REVISION HISTORY:
!  02 Dec 2004 - N. (Eckley) Selin - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    LOGICAL, SAVE          :: IS_INIT = .FALSE.
    LOGICAL                :: LSPLIT
    LOGICAL                :: LDRYD
    LOGICAL                :: LNLPBL
    LOGICAL                :: LGTMM
    LOGICAL                :: LHALOGENCHEM
    INTEGER                :: nAdvect
    INTEGER                :: AS, N

    ! Strings
    CHARACTER(LEN=255)     :: ThisLoc
    CHARACTER(LEN=512)     :: ErrMsg

    ! Pointers
    TYPE(Species), POINTER :: SpcInfo

    !=================================================================
    ! INIT_MERCURY begins here!
    !=================================================================

    ! Assume success
    RC =  GC_SUCCESS

    ! Return if we have already allocated arrays
    IF ( IS_INIT ) RETURN

    ! Initialize
    SpcInfo  => NULL()
    LDRYD    = Input_Opt%LDRYD            ! Use drydep?
    LGTMM    = Input_Opt%LGTMM            ! Use GTMM model?
    LNLPBL   = Input_Opt%LNLPBL           ! Use non-local PBL?
    LSPLIT   = Input_Opt%LSPLIT           ! Tagged simulation?
    nAdvect  = State_Chm%nAdvect          ! # of Hg advected species

    ! Location string for error messages
    ErrMsg   = ''
    ThisLoc  = '-> at DEFINE_TAGGED_Hg (in GeosCore/mercury_mod.F90)'

    ! Store the # of tagged Hg categories in a module variable
    N_Hg_CATS   = State_Chm%N_Hg_CATS

    ! Hg species index corresponding to a given Hg category number
    Hg0_Id_List => State_Chm%Hg0_Id_List
    Hg2_Id_List => State_Chm%Hg2_Id_List
    HgP_Id_List => State_Chm%HgP_Id_List

    ! Locate various species flags
    DO N = 1, State_Chm%nSpecies

       ! Point to the species database entry for species # N
       SpcInfo => State_Chm%SpcData(N)%Info

       SELECT CASE( TRIM( SpcInfo%Name ) )
       CASE( 'Hg0'     )
          id_Hg0    = SpcInfo%ModelId
          ID_Hg_tot = SpcInfo%Hg_Cat
       CASE( 'Hg2'     )
          id_Hg2    = SpcInfo%ModelId
       CASE( 'HgP'     )
          id_HgP    = SpcInfo%ModelId
       CASE( 'Hg0_can' )
          ID_Hg_can = SpcInfo%Hg_Cat
       CASE( 'Hg0_usa' )
          ID_Hg_usa = SpcInfo%Hg_Cat
       CASE( 'Hg0_cam' )
          ID_Hg_cam = SpcInfo%Hg_Cat
       CASE( 'Hg0_sam' )
          ID_Hg_sam = SpcInfo%Hg_Cat
       CASE( 'Hg0_waf' )
          ID_Hg_waf = SpcInfo%Hg_Cat
       CASE( 'Hg0_eaf' )
          ID_Hg_eaf = SpcInfo%Hg_Cat
       CASE( 'Hg0_saf' )
          ID_Hg_saf = SpcInfo%Hg_Cat
       CASE( 'Hg0_naf' )
          ID_Hg_naf = SpcInfo%Hg_Cat
       CASE( 'Hg0_eur' )
          ID_Hg_eur = SpcInfo%Hg_Cat
       CASE( 'Hg0_eeu' )
          ID_Hg_eeu = SpcInfo%Hg_Cat
       CASE( 'Hg0_sov' )
          ID_Hg_sov = SpcInfo%Hg_Cat
       CASE( 'Hg0_mde' )
          ID_Hg_mde = SpcInfo%Hg_Cat
       CASE( 'Hg0_sas' )
          ID_Hg_sas = SpcInfo%Hg_Cat
       CASE( 'Hg0_eas' )
          ID_Hg_eas = SpcInfo%Hg_Cat
       CASE( 'Hg0_sea' )
          ID_Hg_sea = SpcInfo%Hg_Cat
       CASE( 'Hg0_jpn' )
          ID_Hg_jpn = SpcInfo%Hg_Cat
       CASE( 'Hg0_oce' )
          ID_Hg_oce = SpcInfo%Hg_Cat
       CASE( 'Hg0_so'  )
          ID_Hg_so  = SpcInfo%Hg_Cat
       CASE( 'Hg0_bb'  )
          ID_Hg_bb  = SpcInfo%Hg_Cat
       CASE( 'Hg0_geo' )
          ID_Hg_geo = SpcInfo%Hg_Cat
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
       CASE( 'Hg0_str' )
          ID_Hg_str = SpcInfo%Hg_cat
       CASE DEFAULT
          ! Do nothing
       END SELECT

       ! Free pointer
       SpcInfo => NULL()
    ENDDO

    !=================================================================
    ! Allocate arrays
    !=================================================================
    ALLOCATE( COSZM( State_Grid%NX, State_Grid%NY ), STAT=RC )
    CALL GC_CheckVar( 'mercury_mod.F90:COSZM', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    COSZM = 0e+0_fp

    ALLOCATE( EHg0_an( State_Grid%NX, State_Grid%NY ), STAT=RC )
    CALL GC_CheckVar( 'mercury_mod.F90:EHg0_an', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    EHg0_an = 0e+0_fp

    ALLOCATE( EHg0_am( State_Grid%NX, State_Grid%NY ), STAT=RC )
    CALL GC_CheckVar( 'mercury_mod.F90:EHg0_am', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    EHg0_am = 0e+0_fp

    ALLOCATE( EHg2_an( State_Grid%NX, State_Grid%NY ), STAT=RC )
    CALL GC_CheckVar( 'mercury_mod.F90:EHg2_an', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    EHg2_an = 0e+0_fp

    ALLOCATE( EHgP_an( State_Grid%NX, State_Grid%NY ), STAT=RC )
    CALL GC_CheckVar( 'mercury_mod.F90:EHgP_an', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    EHgP_an = 0e+0_fp

    ALLOCATE( EHg0_oc( State_Grid%NX, State_Grid%NY, N_Hg_CATS ), STAT=RC )
    CALL GC_CheckVar( 'mercury_mod.F90:EHg0_oc', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    EHg0_oc = 0e+0_fp

    ALLOCATE( EHg0_dist( State_Grid%NX, State_Grid%NY), STAT=RC )
    CALL GC_CheckVar( 'mercury_mod.F90:EHg0_dist', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    EHg0_dist = 0e+0_fp

    ALLOCATE( EHg0_geo( State_Grid%NX, State_Grid%NY ), STAT=RC )
    CALL GC_CheckVar( 'mercury_mod.F90:EHg0_Geo', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    EHg0_geo = 0e+0_fp

    ALLOCATE( EHg0_bb( State_Grid%NX, State_Grid%NY ), STAT=RC )
    CALL GC_CheckVar( 'mercury_mod.F90:EHg0_bb', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    EHg0_bb = 0e+0_fp

    ALLOCATE( EHg0_snow( State_Grid%NX, State_Grid%NY, N_Hg_CATS ), STAT=RC )
    CALL GC_CheckVar( 'mercury_mod.F90:EHg0_snow', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    EHg0_snow = 0e+0_fp

    ! Allocate the following if tagged Hg simulation
    IF ( LSPLIT ) THEN

       ! Tagged Hg0 arrays (eds 8/31/10)
       ALLOCATE( EHg0_can( State_Grid%NX, State_Grid%NY ), STAT=RC )
       CALL GC_CheckVar( 'mercury_mod.F90:EHg0_can', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       EHg0_can = 0e+0_fp

       ALLOCATE( EHg0_usa( State_Grid%NX, State_Grid%NY ), STAT=RC )
       CALL GC_CheckVar( 'mercury_mod.F90:EHg0_usa', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       EHg0_usa = 0e+0_fp

       ALLOCATE( EHg0_cam( State_Grid%NX, State_Grid%NY ), STAT=RC )
       CALL GC_CheckVar( 'mercury_mod.F90:EHg0_cam', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       EHg0_cam = 0e+0_fp

       ALLOCATE( EHg0_sam( State_Grid%NX, State_Grid%NY ), STAT=RC )
       CALL GC_CheckVar( 'mercury_mod.F90:EHg0_sam', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       EHg0_sam = 0e+0_fp

       ALLOCATE( EHg0_waf( State_Grid%NX, State_Grid%NY ), STAT=RC )
       CALL GC_CheckVar( 'mercury_mod.F90:EHg0_waf', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       EHg0_waf = 0e+0_fp

       ALLOCATE( EHg0_eaf( State_Grid%NX, State_Grid%NY ), STAT=RC )
       CALL GC_CheckVar( 'mercury_mod.F90:EHg0_eaf', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       EHg0_eaf = 0e+0_fp

       ALLOCATE( EHg0_saf( State_Grid%NX, State_Grid%NY ), STAT=RC )
       CALL GC_CheckVar( 'mercury_mod.F90:EHg0_saf', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       EHg0_saf = 0e+0_fp

       ALLOCATE( EHg0_naf( State_Grid%NX, State_Grid%NY ), STAT=RC )
       CALL GC_CheckVar( 'mercury_mod.F90:EHg0_naf', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       EHg0_naf = 0e+0_fp

       ALLOCATE( EHg0_eur( State_Grid%NX, State_Grid%NY ), STAT=RC )
       CALL GC_CheckVar( 'mercury_mod.F90:EHg0_eur', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       EHg0_eur = 0e+0_fp

       ALLOCATE( EHg0_eeu( State_Grid%NX, State_Grid%NY ), STAT=RC )
       CALL GC_CheckVar( 'mercury_mod.F90:EHg0_eeu', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       EHg0_eeu = 0e+0_fp

       ALLOCATE( EHg0_mde( State_Grid%NX, State_Grid%NY ), STAT=RC )
       CALL GC_CheckVar( 'mercury_mod.F90:EHg0_mde', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       EHg0_mde = 0e+0_fp

       ALLOCATE( EHg0_sov( State_Grid%NX, State_Grid%NY ), STAT=RC )
       CALL GC_CheckVar( 'mercury_mod.F90:EHg0_sov', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       EHg0_sov = 0e+0_fp

       ALLOCATE( EHg0_sas( State_Grid%NX, State_Grid%NY ), STAT=RC )
       CALL GC_CheckVar( 'mercury_mod.F90:EHg0_sas', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       EHg0_sas = 0e+0_fp

       ALLOCATE( EHg0_eas( State_Grid%NX, State_Grid%NY ), STAT=RC )
       CALL GC_CheckVar( 'mercury_mod.F90:EHg0_eas', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       EHg0_eas = 0e+0_fp

       ALLOCATE( EHg0_sea( State_Grid%NX, State_Grid%NY ), STAT=RC )
       CALL GC_CheckVar( 'mercury_mod.F90:EHg0_sea', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       EHg0_sea = 0e+0_fp

       ALLOCATE( EHg0_jpn( State_Grid%NX, State_Grid%NY ), STAT=RC )
       CALL GC_CheckVar( 'mercury_mod.F90:EHg0_jpn', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       EHg0_jpn = 0e+0_fp

       ALLOCATE( EHg0_oce( State_Grid%NX, State_Grid%NY ), STAT=RC )
       CALL GC_CheckVar( 'mercury_mod.F90:EHg0_oce', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       EHg0_oce = 0e+0_fp

       ! Tagged Hg2 arrays (eds 8/31/10)
       ALLOCATE( EHg2_can( State_Grid%NX, State_Grid%NY ), STAT=RC )
       CALL GC_CheckVar( 'mercury_mod.F90:EHg2_can', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       EHg2_can = 0e+0_fp

       ALLOCATE( EHg2_usa( State_Grid%NX, State_Grid%NY ), STAT=RC )
       CALL GC_CheckVar( 'mercury_mod.F90:EHg2_usa', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       EHg2_usa = 0e+0_fp

       ALLOCATE( EHg2_cam( State_Grid%NX, State_Grid%NY ), STAT=RC )
       CALL GC_CheckVar( 'mercury_mod.F90:EHg2_cam', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       EHg2_cam = 0e+0_fp

       ALLOCATE( EHg2_sam( State_Grid%NX, State_Grid%NY ), STAT=RC )
       CALL GC_CheckVar( 'mercury_mod.F90:EHg2_sam', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       EHg2_sam = 0e+0_fp

       ALLOCATE( EHg2_waf( State_Grid%NX, State_Grid%NY ), STAT=RC )
       CALL GC_CheckVar( 'mercury_mod.F90:EHg2_waf', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       EHg2_waf = 0e+0_fp

       ALLOCATE( EHg2_eaf( State_Grid%NX, State_Grid%NY ), STAT=RC )
       CALL GC_CheckVar( 'mercury_mod.F90:EHg2_eaf', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       EHg2_eaf = 0e+0_fp

       ALLOCATE( EHg2_saf( State_Grid%NX, State_Grid%NY ), STAT=RC )
       CALL GC_CheckVar( 'mercury_mod.F90:EHg2_saf', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       EHg2_saf = 0e+0_fp

       ALLOCATE( EHg2_naf( State_Grid%NX, State_Grid%NY ), STAT=RC )
       CALL GC_CheckVar( 'mercury_mod.F90:EHg2_naf', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       EHg2_naf = 0e+0_fp

       ALLOCATE( EHg2_eur( State_Grid%NX, State_Grid%NY ), STAT=RC )
       CALL GC_CheckVar( 'mercury_mod.F90:EHg2_eur', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       EHg2_eur = 0e+0_fp

       ALLOCATE( EHg2_eeu( State_Grid%NX, State_Grid%NY ), STAT=RC )
       CALL GC_CheckVar( 'mercury_mod.F90:EHg2_eeu', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       EHg2_eeu = 0e+0_fp

       ALLOCATE( EHg2_mde( State_Grid%NX, State_Grid%NY ), STAT=RC )
       CALL GC_CheckVar( 'mercury_mod.F90:EHg2_mde', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       EHg2_mde = 0e+0_fp

       ALLOCATE( EHg2_sov( State_Grid%NX, State_Grid%NY ), STAT=RC )
       CALL GC_CheckVar( 'mercury_mod.F90:EHg2_sov', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       EHg2_sov = 0e+0_fp

       ALLOCATE( EHg2_sas( State_Grid%NX, State_Grid%NY ), STAT=RC )
       CALL GC_CheckVar( 'mercury_mod.F90:EHg2_sas', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       EHg2_sas = 0e+0_fp

       ALLOCATE( EHg2_eas( State_Grid%NX, State_Grid%NY ), STAT=RC )
       CALL GC_CheckVar( 'mercury_mod.F90:Hg2_eas', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       EHg2_eas = 0e+0_fp

       ALLOCATE( EHg2_sea( State_Grid%NX, State_Grid%NY ), STAT=RC )
       CALL GC_CheckVar( 'mercury_mod.F90:EHg2_sea', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       EHg2_sea = 0e+0_fp

       ALLOCATE( EHg2_jpn( State_Grid%NX, State_Grid%NY ), STAT=RC )
       CALL GC_CheckVar( 'mercury_mod.F90:EHg2_jpn', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       EHg2_jpn = 0e+0_fp

       ALLOCATE( EHg2_oce( State_Grid%NX, State_Grid%NY ), STAT=RC )
       CALL GC_CheckVar( 'mercury_mod.F90:EHg2_oce', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       EHg2_oce = 0e+0_fp

       ! Tagged HgP arrays (eds 8/31/10)
       ALLOCATE( EHgP_can( State_Grid%NX, State_Grid%NY ), STAT=RC )
       CALL GC_CheckVar( 'mercury_mod.F90:EHgP_can', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       EHgP_can = 0e+0_fp

       ALLOCATE( EHgP_usa( State_Grid%NX, State_Grid%NY ), STAT=RC )
       CALL GC_CheckVar( 'mercury_mod.F90:EHgP_usa', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       EHgP_usa = 0e+0_fp

       ALLOCATE( EHgP_cam( State_Grid%NX, State_Grid%NY ), STAT=RC )
       CALL GC_CheckVar( 'mercury_mod.F90:EHgP_cam', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       EHgP_cam = 0e+0_fp

       ALLOCATE( EHgP_sam( State_Grid%NX, State_Grid%NY ), STAT=RC )
       CALL GC_CheckVar( 'mercury_mod.F90:EHgP_sam', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       EHgP_sam = 0e+0_fp

       ALLOCATE( EHgP_waf( State_Grid%NX, State_Grid%NY ), STAT=RC )
       CALL GC_CheckVar( 'mercury_mod.F90:EHgP_waf', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       EHgP_waf = 0e+0_fp

       ALLOCATE( EHgP_eaf( State_Grid%NX, State_Grid%NY ), STAT=RC )
       CALL GC_CheckVar( 'mercury_mod.F90:EHgP_eaf', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       EHgP_eaf = 0e+0_fp

       ALLOCATE( EHgP_saf( State_Grid%NX, State_Grid%NY ), STAT=RC )
       CALL GC_CheckVar( 'mercury_mod.F90:EHgP_saf', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       EHgP_saf = 0e+0_fp

       ALLOCATE( EHgP_naf( State_Grid%NX, State_Grid%NY ), STAT=RC )
       CALL GC_CheckVar( 'mercury_mod.F90:EHgP_naf', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       EHgP_naf = 0e+0_fp

       ALLOCATE( EHgP_eur( State_Grid%NX, State_Grid%NY ), STAT=RC )
       CALL GC_CheckVar( 'mercury_mod.F90:EHgP_eur', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       EHgP_eur = 0e+0_fp

       ALLOCATE( EHgP_eeu( State_Grid%NX, State_Grid%NY ), STAT=RC )
       CALL GC_CheckVar( 'mercury_mod.F90:EHgP_eeu', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       EHgP_eeu = 0e+0_fp

       ALLOCATE( EHgP_mde( State_Grid%NX, State_Grid%NY ), STAT=RC )
       CALL GC_CheckVar( 'mercury_mod.F90:EHgP_mde', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       EHgP_mde = 0e+0_fp

       ALLOCATE( EHgP_sov( State_Grid%NX, State_Grid%NY ), STAT=RC )
       CALL GC_CheckVar( 'mercury_mod.F90:EHgP_sov', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       EHgP_sov = 0e+0_fp

       ALLOCATE( EHgP_sas( State_Grid%NX, State_Grid%NY ), STAT=RC )
       CALL GC_CheckVar( 'mercury_mod.F90:EHgP_sas', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       EHgP_sas = 0e+0_fp

       ALLOCATE( EHgP_eas( State_Grid%NX, State_Grid%NY ), STAT=RC )
       CALL GC_CheckVar( 'mercury_mod.F90:EHgP_eas', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       EHgP_eas = 0e+0_fp

       ALLOCATE( EHgP_sea( State_Grid%NX, State_Grid%NY ), STAT=RC )
       CALL GC_CheckVar( 'mercury_mod.F90:EHgP_sea', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       EHgP_sea = 0e+0_fp

       ALLOCATE( EHgP_jpn( State_Grid%NX, State_Grid%NY ), STAT=RC )
       CALL GC_CheckVar( 'mercury_mod.F90:EHgP_jpn', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       EHgP_jpn = 0e+0_fp

       ALLOCATE( EHgP_oce( State_Grid%NX, State_Grid%NY ), STAT=RC )
       CALL GC_CheckVar( 'mercury_mod.F90:EHgP_oce', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       EHgP_oce = 0e+0_fp

    ENDIF

    IF ( LGTMM ) THEN

       ALLOCATE( EHg0_gtm( State_Grid%NX, State_Grid%NY ), STAT=RC )
       CALL GC_CheckVar( 'mercury_mod.F90:EHg0_gtm', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       EHg0_gtm = 0e+0_fp

    ELSE

       ALLOCATE( EHg0_ln( State_Grid%NX, State_Grid%NY, N_Hg_CATS ), STAT=RC )
       CALL GC_CheckVar( 'mercury_mod.F90:EHg0_ln', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       EHg0_ln = 0e+0_fp

       ALLOCATE( EHg0_vg( State_Grid%NX, State_Grid%NY ), STAT=RC )
       CALL GC_CheckVar( 'mercury_mod.F90:EHg0_vg', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       EHg0_vg = 0e+0_fp

       ALLOCATE( EHg0_so( State_Grid%NX, State_Grid%NY ), STAT=RC )
       CALL GC_CheckVar( 'mercury_mod.F90:EHg0_so', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       EHg0_so = 0e+0_fp

    ENDIF

    ALLOCATE( TCOSZ( State_Grid%NX, State_Grid%NY ), STAT=RC )
    CALL GC_CheckVar( 'mercury_mod.F90:TCOSZ', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    TCOSZ = 0e+0_fp

    ALLOCATE( TTDAY( State_Grid%NX, State_Grid%NY ), STAT=RC )
    CALL GC_CheckVar( 'mercury_mod.F90:TTDAY', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    TTDAY = 0e+0_fp

    ! Allocate ZERO_DVEL if we use non-local PBL mixing or
    ! if drydep is turned off
    IF ( LNLPBL .OR. (.not. LDRYD) ) THEN
       ALLOCATE( ZERO_DVEL( State_Grid%NX, State_Grid%NY ), STAT=RC )
       CALL GC_CheckVar( 'mercury_mod.F90:ZERO_DVEL', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       ZERO_DVEL = 0e+0_fp
    ENDIF

    ALLOCATE( HG2_SEASALT_LOSSRATE( State_Grid%NX, State_Grid%NY ), STAT=RC )
    CALL GC_CheckVar( 'mercury_mod.F90:HG2_SEASALT_LOSSRATE', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    HG2_SEASALT_LOSSRATE = 0e+0_fp

    !=================================================================
    ! Allocate & initialize arrays for tagged species indices
    !=================================================================
    IF ( LSPLIT ) THEN

       ! Species indices for tagged anthro regions
       ALLOCATE( AN_Hg0( State_Grid%NX, State_Grid%NY ), STAT=RC )
       CALL GC_CheckVar( 'mercury_mod.F90:AN_Hg0', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       AN_Hg0 = 0e+0_fp

       ALLOCATE( AN_Hg2( State_Grid%NX, State_Grid%NY ), STAT=RC )
       CALL GC_CheckVar( 'mercury_mod.F90:AN_Hg2', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       AN_Hg2 = 0e+0_fp

       ALLOCATE( AN_HgP( State_Grid%NX, State_Grid%NY ), STAT=RC )
       CALL GC_CheckVar( 'mercury_mod.F90:AN_HgP', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       AN_HgP = 0e+0_fp

    ENDIF

    ! HG_EMIS is needed for non-local PBL mixing
    ALLOCATE ( HG_EMIS( State_Grid%NX, State_Grid%NY, nAdvect ), STAT=RC )
    CALL GC_CheckVar( 'mercury_mod.F90:HG_EMIS', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    HG_EMIS = 0e+0_fp

    !=================================================================
    ! Settings
    !=================================================================

    ! Switch uses ocean rate coefficients from parameter inversion,
    ! ref. Song et al. 2015 ACP
    LOCEANCOEF=.FALSE.

    ! Switch determines whether uptake of Hg2 by sea-salt aerosol
    ! is calculated dynamically (TRUE) or uses a constant rate (FALSE)
    LDYNSEASALT = .TRUE.

    ! Use GEIA 2005 inventory
    LGEIA05=.FALSE.

    ! Switch use NEI2005 and NPRI2005 emission inventories
    LNEI2005=.TRUE.

    ! Switch modifying the speciation profile of Hg emission
    LInPlume=.FALSE.

    ! no Hg emitted through transpiration (VEGEMIS off)
    LVEGEMIS=.FALSE.

    ! Switch adds bromine in marine boundary layer
    L_ADD_MBL_BR=.FALSE.

    ! Switch adds bromine explosion in Northern springtime
    LPOLARBR=.TRUE.

    ! Switch for only doing reduction in cloud water
    LRED_CLOUDONLY = .TRUE.

    ! Switches for new reduction parameterization
    LRED_JNO2=.TRUE. ! Make propto JNO2 otherwise [OH]

    ! Switch for using GEOS-Chem tropospheric bromine fields,
    ! ref. Parrella et al. 2012, instead of older TOMCAT fields
    LGCBROMINE = .TRUE.

    ! Switch specifies that Hg2 is 50% bound to aerosol and 50% in
    ! gas phase (TRUE). If FALSE, then use temperature dependent
    ! partitioning as described in Amos et al. (2011, ACPD)
    LHg2HalfAerosol = .FALSE.

    ! Switch turns on snowpack Hg storage until snowmelt
    LHGSNOW = .TRUE.

    ! Multiplicative factor for increasing stratospheric Br and BrO
    STRAT_BR_FACTOR = 1e+0_fp

    ! Switch turns off all emissions except direct anthropogenic
    LAnthroHgOnly = .FALSE.

    ! Switch turns off all anthropogenic emissions from contiguous USA
    LnoUSAemis = .FALSE.

    !=================================================================
    ! Done
    !=================================================================

    ! Reset IS_INIT, since we have already allocated arrays
    IS_INIT = .TRUE.

  END SUBROUTINE INIT_MERCURY
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: set_options_from_hemco
!
! !DESCRIPTION: Overrides some of the Hg simulation settings depending on
!  the inputs that are specified in the HEMCO configuration file.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE SET_OPTIONS_FROM_HEMCO( Input_Opt, State_Grid, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE HCO_ERROR_MOD
    USE HCO_INTERFACE_MOD, ONLY : HcoState
    USE HCO_ExtList_Mod,   ONLY : GetExtOpt
    USE Input_Opt_Mod,     ONLY : OptInput
    USE State_Grid_Mod,    ONLY : GrdState
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)  :: Input_Opt
    TYPE(GrdState), INTENT(IN)  :: State_Grid
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT) :: RC   ! Success or failure?
!
! !REMARKS:
!
! !REVISION HISTORY:
!  06 Jan 2015 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    LOGICAL :: LGC
    LOGICAL :: LTOMCAT
    LOGICAL :: LPREINDHG
    LOGICAL :: FOUND

    ! Strings
    CHARACTER(LEN=255) :: ThisLoc
    CHARACTER(LEN=512) :: ErrMsg

    !-----------------------------------------------------------------
    ! Initialize
    !-----------------------------------------------------------------
    RC      = HCO_SUCCESS
    ErrMsg  = ''
    ThisLoc = ' -> at SET_OPTIONS_FROM_HEMCO (in GeosCore/mercury_mod.F90)'

    !-----------------------------------------------------------------
    ! Set the value of chemistry flags depending on whether or not
    ! the HEMCO collection LFLAGNAME is activated
    !-----------------------------------------------------------------
    CALL GetExtOpt( HcoState%Config, -999, 'LRED_JNO2', &
                    OptValBool=LRED_JNO2, FOUND=FOUND, RC=RC )
    IF ( RC /= HCO_SUCCESS ) THEN
       ErrMsg = 'LRED_JNO2 not found in HEMCO_Config.rc file!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    IF ( .not. FOUND ) THEN
       LRED_JNO2 = .FALSE.
    ENDIF

    CALL GetExtOpt( HcoState%Config, -999, 'LHALOGENCHEM', &
                    OptValBool=LHALOGENCHEM, FOUND=FOUND, RC=RC )
    IF ( RC /= HCO_SUCCESS ) THEN
       ErrMsg = 'LHALOGENCHEM not found in HEMCO_Config.rc file!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    IF ( .not. FOUND ) THEN
       LHALOGENCHEM = .TRUE.
    ENDIF

    CALL GetExtOpt( HcoState%Config, -999, 'LHGAQCHEM', &
                    OptValBool=LHGAQCHEM, FOUND=FOUND, RC=RC )
    IF ( RC /= HCO_SUCCESS ) THEN
       ErrMsg = 'LHGAQCHEM not found in HEMCO_Config.rc file!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    IF ( .not. FOUND ) THEN
       LHGAQCHEM = .FALSE.
    ENDIF

    CALL GetExtOpt( HcoState%Config, -999, 'LBRCHEM', &
                    OptValBool=LBRCHEM, FOUND=FOUND, RC=RC )
    IF ( RC /= HCO_SUCCESS ) THEN
       ErrMsg = 'LBRCHEM not found in HEMCO_Config.rc file!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    IF ( .not. FOUND ) THEN
       LBRCHEM = .FALSE.
    ENDIF

    CALL GetExtOpt( HcoState%Config, -999, 'LBROCHEM', &
                    OptValBool=LBROCHEM, FOUND=FOUND, RC=RC )
    IF ( RC /= HCO_SUCCESS ) THEN
       ErrMsg = 'LBROCHEM not found in HEMCO_Config.rc file!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    IF ( .not. FOUND ) THEN
       LBROCHEM = .FALSE.
    ENDIF

    CALL GetExtOpt( HcoState%Config, -999, 'LOHO3CHEM', &
                    OptValBool=LOHO3CHEM, FOUND=FOUND, RC=RC )
    IF ( RC /= HCO_SUCCESS ) THEN
       ErrMsg = 'LOHO3CHEM not found in HEMCO_Config.rc file!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    IF ( .not. FOUND ) THEN
       LOHO3CHEM = .FALSE.
    ENDIF

    !-----------------------------------------------------------------
    ! Set the value of LGCBROMINE depending on the values of the
    ! HEMCO collections BrOx_GC and BrOx_TOMCAT
    !-----------------------------------------------------------------

    ! First look for BrOx_GC
    CALL GetExtOpt( HcoState%Config, -999, 'BrOx_GC', &
                    OptValBool=LGC, FOUND=FOUND, RC=RC )
    IF ( RC /= HCO_SUCCESS ) THEN
       ErrMsg = 'BrOx_GC not found in HEMCO_Config.rc file!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    IF ( .not. FOUND ) THEN
       LGC = .FALSE.
    ENDIF

    ! Set LGCBROMINE = .TRUE. if BrOx_GC is true
    LGCBROMINE = LGC

    ! Set BrOx_TOMCAT to be the opposite of LGCBROMINE
    LTOMCAT   = ( .not. LGCBROMINE )

    ! Are we doing a preindustrial simulation?
    LPREINDHG = Input_Opt%LPREINDHG

    !-----------------------------------------------------------------
    ! Set the value of LNEI2005 depending on whether or not
    ! the HEMCO collection NEI2005 is activated
    !-----------------------------------------------------------------
    CALL GetExtOpt( HcoState%Config, -999, 'NEI2005', &
                    OptValBool=LNEI2005, FOUND=FOUND, RC=RC )
    IF ( RC /= HCO_SUCCESS ) THEN
       ErrMsg = 'NEI2005 not found in HEMCO_Config.rc file!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    IF ( .not. FOUND ) THEN
       LNEI2005 = .FALSE.
    ENDIF

    !-----------------------------------------------------------------
    ! Set the value of LInPlume depending on whether or not
    ! the HEMCO collection NEI2005 is activated
    !-----------------------------------------------------------------
    CALL GetExtOpt( HcoState%Config, -999, 'LRED_INPLUME', &
                    OptValBool=LInPlume, FOUND=FOUND, RC=RC )
    IF ( RC /= HCO_SUCCESS ) THEN
       ErrMsg = 'LRED_INPLUME not found in HEMCO_Config.rc file!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    IF ( .not. FOUND ) THEN
       LInPlume = .FALSE.
    ENDIF

    ! In plume degradation of Hg2 by SO2 in U.S. and Canada at CFPPs,
    ! (yzh 11/1/2011).  Move this here so that the HEMCO_Config file
    ! will have been already read by this point. (bmy, 10/11/16)
    IF ( LInPlume .AND. .NOT. LPREINDHG .AND. LNEI2005 ) THEN
       CALL DO_RED_INPLUME( Input_Opt, State_Grid, RC )
    ENDIF

    !-----------------------------------------------------------------
    ! Echo output
    !-----------------------------------------------------------------
    WRITE( 6, '(a)' ) REPEAT( '=', 79 )
    WRITE( 6, 100   )
    WRITE( 6, 110   ) LRED_JNO2
    WRITE( 6, 120   ) LGCBROMINE
    WRITE( 6, 130   ) LNEI2005
    WRITE( 6, 140   ) LInPlume
    WRITE( 6, 150   ) LHALOGENCHEM
    WRITE( 6, 160   ) LHGAQCHEM
    WRITE( 6, 170   ) LBRCHEM
    WRITE( 6, 180   ) LBROCHEM
    WRITE( 6, 190   ) LOHO3CHEM
    WRITE( 6, '(a)' ) REPEAT( '=', 79 )
100 FORMAT( 'Adjusting Hg simulation settings from HEMCO inputs' )
110 FORMAT( 'LRED_JNO2  is set to ', L1                          )
120 FORMAT( 'LGCBROMINE is set to ', L1                          )
130 FORMAT( 'LNEI2005   is set to ', L1                          )
140 FORMAT( 'LInPlume   is set to ', L1                          )
150 FORMAT( 'LHALOGENCHEM   is set to ', L1                          )
160 FORMAT( 'LHGAQCHEM  is set to ', L1                          )
170 FORMAT( 'LBRCHEM    is set to ', L1                          )
180 FORMAT( 'LBROCHEM   is set to ', L1                          )
190 FORMAT( 'LOHO3CHEM  is set to ', L1                          )

  END SUBROUTINE SET_OPTIONS_FROM_HEMCO
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cleanup_mercury
!
! !DESCRIPTION: Subroutine CLEANUP\_MERCURY deallocates all module arrays.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE CLEANUP_MERCURY
!
! !REVISION HISTORY:
!  06 Dec 2004 - N. (Eckley) Selin - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC

    IF ( ALLOCATED( AN_Hg0   ) ) DEALLOCATE( AN_Hg0   )
    IF ( ALLOCATED( AN_Hg2   ) ) DEALLOCATE( AN_Hg2   )
    IF ( ALLOCATED( AN_HgP   ) ) DEALLOCATE( AN_HgP   )
    IF ( ALLOCATED( COSZM    ) ) DEALLOCATE( COSZM    )
    IF ( ALLOCATED( EHg0_an  ) ) DEALLOCATE( EHg0_an  )
    IF ( ALLOCATED( EHg0_can ) ) DEALLOCATE( EHg0_can )
    IF ( ALLOCATED( EHg0_usa ) ) DEALLOCATE( EHg0_usa )
    IF ( ALLOCATED( EHg0_sam ) ) DEALLOCATE( EHg0_sam )
    IF ( ALLOCATED( EHg0_eaf ) ) DEALLOCATE( EHg0_eaf )
    IF ( ALLOCATED( EHg0_waf ) ) DEALLOCATE( EHg0_waf )
    IF ( ALLOCATED( EHg0_saf ) ) DEALLOCATE( EHg0_saf )
    IF ( ALLOCATED( EHg0_naf ) ) DEALLOCATE( EHg0_naf )
    IF ( ALLOCATED( EHg0_eur ) ) DEALLOCATE( EHg0_eur )
    IF ( ALLOCATED( EHg0_eeu ) ) DEALLOCATE( EHg0_eeu )
    IF ( ALLOCATED( EHg0_sov ) ) DEALLOCATE( EHg0_sov )
    IF ( ALLOCATED( EHg0_mde ) ) DEALLOCATE( EHg0_mde )
    IF ( ALLOCATED( EHg0_sas ) ) DEALLOCATE( EHg0_sas )
    IF ( ALLOCATED( EHg0_eas ) ) DEALLOCATE( EHg0_eas )
    IF ( ALLOCATED( EHg0_sea ) ) DEALLOCATE( EHg0_sea )
    IF ( ALLOCATED( EHg0_jpn ) ) DEALLOCATE( EHg0_jpn )
    IF ( ALLOCATED( EHg0_oce ) ) DEALLOCATE( EHg0_oce )
    IF ( ALLOCATED( EHg0_am  ) ) DEALLOCATE( EHg0_am  )
    IF ( ALLOCATED( EHg2_an  ) ) DEALLOCATE( EHg2_an  )
    IF ( ALLOCATED( EHg2_can ) ) DEALLOCATE( EHg2_can )
    IF ( ALLOCATED( EHg2_usa ) ) DEALLOCATE( EHg2_usa )
    IF ( ALLOCATED( EHg2_sam ) ) DEALLOCATE( EHg2_sam )
    IF ( ALLOCATED( EHg2_eaf ) ) DEALLOCATE( EHg2_eaf )
    IF ( ALLOCATED( EHg2_waf ) ) DEALLOCATE( EHg2_waf )
    IF ( ALLOCATED( EHg2_saf ) ) DEALLOCATE( EHg2_saf )
    IF ( ALLOCATED( EHg2_naf ) ) DEALLOCATE( EHg2_naf )
    IF ( ALLOCATED( EHg2_eur ) ) DEALLOCATE( EHg2_eur )
    IF ( ALLOCATED( EHg2_eeu ) ) DEALLOCATE( EHg2_eeu )
    IF ( ALLOCATED( EHg2_sov ) ) DEALLOCATE( EHg2_sov )
    IF ( ALLOCATED( EHg2_mde ) ) DEALLOCATE( EHg2_mde )
    IF ( ALLOCATED( EHg2_sas ) ) DEALLOCATE( EHg2_sas )
    IF ( ALLOCATED( EHg2_eas ) ) DEALLOCATE( EHg2_eas )
    IF ( ALLOCATED( EHg2_sea ) ) DEALLOCATE( EHg2_sea )
    IF ( ALLOCATED( EHg2_jpn ) ) DEALLOCATE( EHg2_jpn )
    IF ( ALLOCATED( EHg2_oce ) ) DEALLOCATE( EHg2_oce )
    IF ( ALLOCATED( EHgP_an  ) ) DEALLOCATE( EHgP_an  )
    IF ( ALLOCATED( EHgP_can ) ) DEALLOCATE( EHgP_can )
    IF ( ALLOCATED( EHgP_usa ) ) DEALLOCATE( EHgP_usa )
    IF ( ALLOCATED( EHgP_sam ) ) DEALLOCATE( EHgP_sam )
    IF ( ALLOCATED( EHgP_eaf ) ) DEALLOCATE( EHgP_eaf )
    IF ( ALLOCATED( EHgP_waf ) ) DEALLOCATE( EHgP_waf )
    IF ( ALLOCATED( EHgP_saf ) ) DEALLOCATE( EHgP_saf )
    IF ( ALLOCATED( EHgP_naf ) ) DEALLOCATE( EHgP_naf )
    IF ( ALLOCATED( EHgP_eur ) ) DEALLOCATE( EHgP_eur )
    IF ( ALLOCATED( EHgP_eeu ) ) DEALLOCATE( EHgP_eeu )
    IF ( ALLOCATED( EHgP_sov ) ) DEALLOCATE( EHgP_sov )
    IF ( ALLOCATED( EHgP_mde ) ) DEALLOCATE( EHgP_mde )
    IF ( ALLOCATED( EHgP_sas ) ) DEALLOCATE( EHgP_sas )
    IF ( ALLOCATED( EHgP_eas ) ) DEALLOCATE( EHgP_eas )
    IF ( ALLOCATED( EHgP_sea ) ) DEALLOCATE( EHgP_sea )
    IF ( ALLOCATED( EHgP_jpn ) ) DEALLOCATE( EHgP_jpn )
    IF ( ALLOCATED( EHgP_oce ) ) DEALLOCATE( EHgP_oce )
    IF ( ALLOCATED( EHg0_oc  ) ) DEALLOCATE( EHg0_oc  )
    IF ( ALLOCATED( EHg0_ln  ) ) DEALLOCATE( EHg0_ln  )
    IF ( ALLOCATED( EHg0_snow) ) DEALLOCATE( EHg0_snow)
    IF ( ALLOCATED( EHg0_geo ) ) DEALLOCATE( EHg0_geo )
    IF ( ALLOCATED( EHg0_bb  ) ) DEALLOCATE( EHg0_bb  )
    IF ( ALLOCATED( EHg0_gtm ) ) DEALLOCATE( EHg0_gtm )
    IF ( ALLOCATED( EHg0_vg  ) ) DEALLOCATE( EHg0_vg  )
    IF ( ALLOCATED( EHg0_so  ) ) DEALLOCATE( EHg0_so  )
    IF ( ALLOCATED( EHg0_dist) ) DEALLOCATE( EHg0_dist)
    IF ( ALLOCATED( TCOSZ    ) ) DEALLOCATE( TCOSZ    )
    IF ( ALLOCATED( TTDAY    ) ) DEALLOCATE( TTDAY    )
    IF ( ALLOCATED( ZERO_DVEL) ) DEALLOCATE( ZERO_DVEL)
    IF ( ALLOCATED( HG_EMIS  ) ) DEALLOCATE( HG_EMIS  )
    IF ( ALLOCATED( HG2_SEASALT_LOSSRATE ) ) DEALLOCATE( HG2_SEASALT_LOSSRATE )

    ! Free pointers to HEMCO fields
    O3            => NULL()
    OH_trop       => NULL()
    OH_strat      => NULL()
    JNO2          => NULL()
    HEM_NO        => NULL()
    HEM_NO2       => NULL()
    HEM_CLO       => NULL()
    HEM_CL        => NULL()
    HEM_OA        => NULL()
    HEM_HOCl      => NULL()
    HEM_HO2_trop  => NULL()
    HEM_HO2_strat => NULL()
    OCEAN_CONC    => NULL()

    ! Free Hg indexing pointers
    Hg0_Id_List => NULL()
    Hg2_Id_List => NULL()
    HgP_Id_List => NULL()

  END SUBROUTINE CLEANUP_MERCURY
!EOC
END MODULE MERCURY_MOD
