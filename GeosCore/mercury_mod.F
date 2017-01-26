!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: mercury_mod.F
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
      USE OCEAN_MERCURY_MOD, ONLY : LGCAPEMIS
      USE OCEAN_MERCURY_MOD, ONLY : LPOLARBR
      USE OCEAN_MERCURY_MOD, ONLY : LGEIA05
      USE OCEAN_MERCURY_MOD, ONLY : LVEGEMIS
      USE OCEAN_MERCURY_MOD, ONLY : LBRCHEM
      USE OCEAN_MERCURY_MOD, ONLY : LRED_JNO2
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
      PRIVATE :: GET_LWC
      PRIVATE :: GET_VCLDF
      PRIVATE :: GET_O3
      PRIVATE :: GET_OH
      PRIVATE :: GET_BR
      PRIVATE :: GET_JNO2
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
!  (1 ) Updated for reduction rxn and online Hg0 ocean flux.  Now use 
!        diagnostic arrays from "diag03_mod.f".  (eck, sas, bmy, 1/21/05)
!  (2 ) Now references "pbl_mix_mod.f".  Remove FEMIS array and routine
!        COMPUTE_FEMIS. (bmy, 2/15/05)
!  (3 ) Now can read data for both GEOS and GCAP grids (bmy, 8/16/05)
!  (4 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (5 ) Now references XNUMOL from "tracer_mod.f" (bmy, 10/25/05)
!  (6 ) Various updates added for tagged Hg sim. (eck, sas, cdh, bmy, 4/6/06)
!  (7 ) Now includes LPREINDHG logical switch for preindustrial simulation 
!       (eds 7/30/08)
!  13 Aug 2010 - R. Yantosca - Add modifications for MERRA (treat like GEOS-5)
!  01 Apr 2011 - H. Amos     - added LVEGEMIS, distinct from LGCAPEMIS
!  27 Sep 2011 - H. Amos     - Remove LHg_WETDasHNO3 logical, it's obsolete
!  01 Nov 2011 - Y. Zhang    - Add subroutines to do in plume reduction of Hg
!  07 Feb 2012 - E. Corbitt  - Imported tagged tracers; renamed nt to geo.
!  08 Feb 2012 - R. Yantosca - Add modifications for GEOS-5.7.x met
!  01 Mar 2012 - R. Yantosca - Now reference new grid_mod.F90
!  01 Mar 2012 - R. Yantosca - Use updated GET_LOCALTIME from time_mod.F
!  19 Apr 2012 - E. Corbitt  - Added LGCBROMINE to use GEOS-Chem bromine.
!  03 Jun 2013 - C. Holmes   - Extensively rewritten for simplicity
!  03 Jun 2013 - R. Yantosca - Added ProTeX headers
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!  17 Nov 2014 - M. Yannetti - Added PRECISION_MOD
!  12 Mar 2015 - R. Yantosca - Now get JNO2 from the HEMCO data structure
!  12 Mar 2015 - R. Yantosca - Remove GET_GLOBAL_JNO2 routine
!  10 Aug 2015 - J. Fisher   - Updated treatment of polar BrO chemistry
!                              (orig. 29/11/11)
!  05 Jan 2016 - E. Lundgren - Now use global dry air gas and Henry's Law consts
!  31 May 2016 - E. Lundgren - Replace Input_Opt%XNUMOL with AVO / (emMW_g*1e-3)
!                              where emMW_g is emitted MW from species database
!  29 Nov 2016 - R. Yantosca - grid_mod.F90 is now gc_grid_mod.F90
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
      REAL(f4), POINTER :: O3  (:,:,:) => NULL()
      REAL(f4), POINTER :: OH  (:,:,:) => NULL()
      REAL(f4), POINTER :: JNO2(:,:,:) => NULL()

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
!
! !DEFINED PARAMETERS:
!
      ! Henry's Law constant for Hg2 [mol /L /atm]
      ! Now use Henry's Law constant defined in species database (ewl, 1/5/15)
!      REAL(fp), PARAMETER    :: HL = 1.4e+6_fp

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
! !IROUTINE: chemmercury
!
! !DESCRIPTION: Subroutine CHEMMERCURY is the driver routine for mercury 
!  chemistry in the GEOS-CHEM module. 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CHEMMERCURY( am_I_Root, Input_Opt,
     &                        State_Met, State_Chm, RC )
!
! !USES:
!
      USE CMN_SIZE_MOD
      USE CMN_DIAG_MOD
      USE DAO_MOD,            ONLY : IS_ICE
      USE DAO_MOD,            ONLY : IS_LAND
      USE DAO_MOD,            ONLY : IS_WATER
      USE DEPO_MERCURY_MOD,   ONLY : ADD_HG2_DD
      USE DEPO_MERCURY_MOD,   ONLY : ADD_HgP_DD
#if defined( BPCH_DIAG )
      USE DIAG_MOD,           ONLY : AD44
#endif
      USE DIAG03_MOD,         ONLY : AD03_Hg2_Hg0, AD03_Hg2_O3
      USE DIAG03_MOD,         ONLY : AD03_Hg2_Br,  AD03_Hg2_OH
      USE DIAG03_MOD,         ONLY : AD03_Hg2_SS,  LD03
      USE DIAG03_MOD,         ONLY : AD03_Hg2_SSR, ND03 
      USE DIAG03_MOD,         ONLY : AD03_Br
      USE ErrCode_Mod
      USE ERROR_MOD,          ONLY : DEBUG_MSG
      USE ERROR_MOD,          ONLY : ERROR_STOP
      USE ERROR_MOD,          ONLY : SAFE_DIV
      USE HCO_INTERFACE_MOD,  ONLY : HcoState
      USE HCO_EmisList_Mod,   ONLY : HCO_GetPtr
      USE GLOBAL_BR_MOD,      ONLY : GET_GLOBAL_BR
      USE GC_GRID_MOD,        ONLY : GET_AREA_CM2
      USE Input_Opt_Mod,      ONLY : OptInput
      USE PBL_MIX_MOD,        ONLY : GET_FRAC_UNDER_PBLTOP
      USE Species_Mod,        ONLY : Species
      USE State_Chm_Mod,      ONLY : ChmState
      USE State_Met_Mod,      ONLY : MetState
      USE TIME_MOD,           ONLY : GET_MONTH
      USE TIME_MOD,           ONLY : GET_TS_CHEM
      USE TIME_MOD,           ONLY : ITS_A_NEW_MONTH
      USE TIME_MOD,           ONLY : ITS_TIME_FOR_A3 
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
!  08 Dec 2009 - R. Yantosca - Added ProTeX headers
!  06 Dec 2004 - N. (Eckley) Selin - Initial version
!  (1 ) Now references routine GET_PBL_MAX_L from "pbl_mix_mod.f".  Now
!        references AD44 from "diag_mod.f".  Now sum the levels from T44 into 
!        the AD44 array.  Now references N_TRACERS from "tracer_mod.f".
!        (bmy, 2/24/05)
!  (2 ) Bug fix: Set T44 to 0e0 for single precision.  Now allow for zero
!        dry deposition velocity.  Now call INIT_MERCURY from "input_mod.f"
!        (bmy, 4/6/06)
!  02 Apr 2013 - C. Holmes - CHEMMERCURY rewritten for clarity and easier
!        maintainability. The chemistry and deposition for particle-bound
!        mercury is now solved simultaneously with Hg(0) and gaseous Hg(II).
!        The subroutine now uses a 4th-order Runge-Kutta method to solve
!        the coupled chemistry. The subroutines CHEM_HG0_HG2 and CHEM_HGP 
!        are now obsolete.
!  04 Apr 2013 - C. Holmes   - Use State_Met for met fields.
!  03 Jun 2013 - R. Yantosca - Now use fields from Input_OPt instead of from
!                              logical_mod.F and tracer_mod.F
!  26 Sep 2013 - R. Yantosca - Renamed GEOS_57 Cpp switch to GEOS_FP
!  23 Oct 2013 - R. Yantosca - Now pass objects to GET_GLOBAL_OH routine
!  06 Nov 2014 - R. Yantosca - Now use State_Met%AIRDEN(I,J,L)
!  06 Nov 2014 - R. Yantosca - Now use State_Met%CLDF(I,J,L)
!  26 Feb 2015 - E. Lundgren - Replace GET_PCENTER with State_Met%PMID_DRY.
!                              Remove dependency on pressure_mod.
!  12 Aug 2015 - R. Yantosca - Added support for MERRA2 meteorology
!  05 Jan 2016 - E. Lundgren - Use Henry's constant from species database and
!                              calculate dry air gas constant from global params
!  25 Apr 2016 - R. Yantosca - Now pass Hg category # to ADD_Hg2_* routines
!  29 Apr 2016 - R. Yantosca - Don't initialize pointers in declaration stmts
!  30 Jun 2016 - R. Yantosca - Remove instances of STT.  Now get the advected
!                              species ID from State_Chm%Map_Advect.
!  01 Jul 2016 - R. Yantosca - Now rename species DB object ThisSpc to SpcInfo
!  11 Oct 2016 - R. Yantosca - Pass additional args to SET_OPTIONS_FROM_HEMCO
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:

      ! SAVEd scalars
      LOGICAL, SAVE      :: FIRST = .TRUE.

      ! Scalars
      LOGICAL            :: LDYNOCEAN
      LOGICAL            :: LGTMM 
      LOGICAL            :: LNLPBL
      LOGICAL            :: LPRT
      INTEGER            :: NA, nAdvect
      INTEGER            :: I, J, L, K, N, NN, MONTH, NSTEPS, Hg_Cat
      REAL(fp)           :: K_DRYD0, K_DRYD2G, K_DRYD2P, K_SALT
      REAL(fp)           :: C_OH, C_O3, C_BR, C_BRO
      REAL(fp)           :: K_OH, K_O3, K_BR, K_BRO
      REAL(fp)           :: K_OX, K_RED, K_DEP0, K_DEP2G, K_DEP2P
      REAL(fp)           :: DEP_Hg0, DEP_HG2G, DEP_HG2P
      REAL(fp)           :: DEP_HG2G_DRY, DEP_HG2G_SALT, DEP_HG2P_DRY
      REAL(fp)           :: DEP_HG2_DRY
      REAL(fp)           :: GROSS_OX, GROSS_OX_BR
      REAL(fp)           :: GROSS_OX_OH, GROSS_OX_O3
      REAL(fp)           :: GROSS_RED, NET_OX
      REAL(fp)           :: F_PBL, FC, Faq, LWC
      REAL(fp)           :: AREA_CM2, DEP_DRY_FLX
      REAL(fp)           :: F_HG0_DEP
      REAL(fp)           :: DTCHEM, DT
      REAL(fp)           :: KMAX
      REAL(fp)           :: A(9,9), Xold(9), Xnew(9)
      REAL(fp)           :: s1(9), s2(9), s3(9), s4(9)
      REAL(fp)           :: K0, R_DRY

      ! Pointers
      REAL(fp), POINTER  :: Spc(:,:,:,:)
      REAL(fp), POINTER  :: T(:,:,:)

      ! Strings
      CHARACTER(LEN=255) :: LOC = 'CHEMMERCURY (mercury_mod.F)'
!
! !DEFINED PARAMETERS:
!
      ! Now calculate the dry air gas constant in [L atm/K/mol] from
      ! global parameters that are in SI units (ewl, 1/5/16)
      ! Gas constant [L atm /K /mol]
      !REAL(fp), PARAMETER     :: R        = 8.2e-2_fp

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

      ! Reduce K_RED_JNO2 for GEIA 2005 inventory
      ! recommended by Helen Amos, 23 Sep 2011
      ! change implemented 10/19/11 by eck
      ! inserted ifdefs for future tuning of 2x2.5 by eck
      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      !%%% NOTE: Someone will have to add tunings for the GEOS-5.7.2
      !%%% nested grids at 0.25 x 0.3125 resolution (bmy, 6/3/13)
      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      ! These K_RED_JNO2 parameters have been defined based on GEOS-5
      ! runs and remain appropriate for GEOS-5 at 4x5 in v11-01c. Based
      ! on the GEOS-FP 4x5 run, we need to increase K_RED_JNO2 for
      ! GEOS-FP. In the absences of benchmark runs for GEOS-FP at other
      ! resolutions, we assume the same decrease in GEOS-FP relative to
      ! GEOS-5, regardless of resolution. This should be adapated and
      ! tested by other users as appropriate (J. Fisher, 3/2016)
      ! Although we have not tested other met fields, we assume same
      ! parameters apply for GEOS-FP, MERRA, and MERRA-2.
      ! To simplify this code, for now we remove the various resolution
      ! differences (which haven't been tested) and simply include
      ! different versions for GEOS-FP vs GEOS-5 and nested vs.
      ! non-nested. Users are encouraged to test and adapt. (J. Fisher,
      ! 4/2016)
#if   defined( MERRA ) || defined( GEOS_FP ) || defined( MERRA2 )
#if   defined( NESTED_NA )
      REAL(fp), PARAMETER     :: K_RED_JNO2 = 9.6e-3_fp * 1.33
#else
      REAL(fp), PARAMETER     :: K_RED_JNO2 = 9.6e-3_fp
#endif
#else
#if   defined( NESTED_NA )
      REAL(fp), PARAMETER     :: K_RED_JNO2 = 6.4e-3_fp * 1.33
#else
      REAL(fp), PARAMETER     :: K_RED_JNO2 = 6.4e-3_fp
#endif
#endif

!!!#if defined( GRID4x5 )
!!!      REAL(fp), PARAMETER     :: K_RED_JNO2 = 6.4e-3_fp * 1.5e+0_fp
!!!#elif defined( GRID2x25 )
!!!      REAL(fp), PARAMETER     :: K_RED_JNO2 = 3.5e-3_fp/2e+0_fp * 1.5e+0_fp
!!!#elif defined( GRID05x0625 ) && defined( NESTED_NA )
!!!      REAL(fp), PARAMETER     :: K_RED_JNO2 = 6.4e-3_fp * 1.33 * 1.5e+0_fp
!!!#elif defined( GRID05x0666 ) && defined( NESTED_NA )
!!!      REAL(fp), PARAMETER     :: K_RED_JNO2 = 6.4e-3_fp * 1.33 * 1.5e+0_fp
!!!#else
!!!      REAL(fp), PARAMETER     :: K_RED_JNO2 = 6.4e-3_fp  * 1.5e+0_fp
!!!#endif   
!!!#else
!!!#if defined( GRID4x5 )
!!!      !REAL(fp), PARAMETER     :: K_RED_JNO2 = 3.5d-3 / 2e+0_fp ! ref: Amos et al. (2011)
!!!      REAL(fp), PARAMETER     :: K_RED_JNO2 = 6.4e-3_fp ! yzh, consider cloud fraction 
!!!#elif defined( GRID2x25 )
!!!      REAL(fp), PARAMETER     :: K_RED_JNO2 = 3.5e-3_fp/2e+0_fp
!!!#elif defined( GRID05x0625 ) && defined( NESTED_NA )
!!!      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!      !%%% Bob Y. used same value as for GRID05x0666             %%%
!!!      !%%% Team Hg needs to supply a better value (bmy, 8/12/15) %%%
!!!      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!      REAL(fp), PARAMETER     :: K_RED_JNO2 = 6.4e-3_fp * 1.33 ! yzh
!!!#elif defined( GRID05x0666 ) && defined( NESTED_NA )
!!!      REAL(fp), PARAMETER     :: K_RED_JNO2 = 6.4e-3_fp * 1.33 ! yzh
!!!#else
!!!      REAL(fp), PARAMETER     :: K_RED_JNO2 = 6.4e-3_fp        ! yzh, need tune 
!!!#endif   
!!!#endif   

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

      ! Copy values from Input_Opt
      LDYNOCEAN = Input_Opt%LDYNOCEAN
      LGTMM     = Input_Opt%LGTMM
      LNLPBL    = Input_Opt%LNLPBL
      LPRT      = ( Input_Opt%LPRT .and. am_I_Root )

      ! Initialize pointers
      Spc      => State_Chm%Species   ! Chemical species array [kg]
      T        => State_Met%T         ! Temperature [K]
      SpcInfo  => NULL()

      ! Get info about this species from the database
      SpcInfo => State_Chm%SpcData(id_Hg2)%Info

      ! Get Henry's Law constant for Hg2 [mol/L/atm] (ewl, 1/5/16)
      K0 = SpcInfo%Henry_K0

      ! Calculate dry air gas constant in [L atm/K/mol] (ewl, 1/5/16)
      R_DRY = Rd * AIRMW / ATM

      ! Override the default value of LGCBROMINE and LRED_JNO2 depending 
      ! on the settings in the HEMCO configuration file.  NOTE: This cannot
      ! be done in INIT_MERCURY because at that time the HEMCO configuration 
      ! file has not been read yet. (bmy, 3/12/15)
      IF ( FIRST ) THEN
         CALL SET_OPTIONS_FROM_HEMCO( am_I_Root, Input_Opt, RC )

#if !defined( MERRA ) && !defined( GEOS_FP ) && !defined( MERRA2 )
         ! Warning for non-MERRA
         WRITE( 6, * ) 'WARNING: Polar BrO parameterization '//
     &                 'designed for use with incident shortwave'
         WRITE( 6, * ) 'radiation. Only net shortwave radiation '//
     &                 'available.'
#endif
         FIRST = .FALSE.
      ENDIF

      !=================================================================
      ! Read monthly mean Br, OH, O3, and J_NO2 fields
      !=================================================================
      IF ( ITS_A_NEW_MONTH() ) THEN 

         ! Get the current month
         MONTH = GET_MONTH()

         ! Get a pointer to the monthly mean OH from HEMCO (bmy, 3/11/15)
         CALL HCO_GetPtr( am_I_Root, HcoState, 'GLOBAL_OH', OH,   RC )
         IF ( RC /= GC_SUCCESS ) THEN
            CALL ERROR_STOP ( 'Cannot get pointer to GLOBAL_OH', LOC )
         ENDIF

         ! Get a pointer to the monthly mean O3 from HEMCO (bmy, 3/11/15)
         CALL HCO_GetPtr( am_I_Root, HcoState, 'GLOBAL_O3', O3,  RC )
         IF ( RC /= GC_SUCCESS ) THEN
            CALL ERROR_STOP ( 'Cannot get pointer to GLOBAL_O3', LOC )
         ENDIF

         CALL GET_GLOBAL_BR( am_I_Root, Input_Opt, State_Met, MONTH, RC)
         IF ( LPRT ) CALL DEBUG_MSG( '### CHEMMERC: a GET_GLOBAL_BR' )

         ! Get a pointer to the monthly mean OH from HEMCO (bmy, 3/11/15)
         IF ( LRED_JNO2 ) THEN
            CALL HCO_GetPtr( am_I_Root, HcoState, 'JNO2', JNO2,  RC )
            IF ( RC /= GC_SUCCESS ) THEN
               CALL ERROR_STOP ( 'Cannot get pointer to JNO2', LOC )
            ENDIF
         ENDIF

      ENDIF
      
      !=================================================================
      ! Perform chemistry on Hg species 
      !=================================================================
      
      ! Chemistry timestep [s]
      DTCHEM = GET_TS_CHEM() * 60e+0_fp

      ! Compute diurnal scaling for OH
      CALL OHNO3TIME
      IF ( LPRT ) CALL DEBUG_MSG( 'CHEMMERCURY: a OHNO3TIME' )

      ! Calculate the rate of sea salt aerosol uptake of Hg2
      IF ( LDYNSEASALT .AND. ITS_TIME_FOR_A3() ) THEN
         CALL CALC_HG2_SEASALT_LOSSRATE( State_Met )
         IF ( LPRT ) CALL DEBUG_MSG( 'CHEMMERCURY: a SEASALT_LOSSRATE' )
      ENDIF

      ! Partition Hg(II) between gas and aerosol
      CALL PARTITIONHG2( am_I_Root, Input_Opt, State_Chm, RC )

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I,       J,        L,        K,       N, NN, NSTEPS )
!$OMP+PRIVATE( K_DRYD0, K_DRYD2G, K_DRYD2P, K_SALT )
!$OMP+PRIVATE( C_OH,    C_O3,     C_BR,     C_BRO )
!$OMP+PRIVATE( K_OH,    K_O3,     K_BR,     K_BRO )
!$OMP+PRIVATE( K_OX,    K_RED,    K_DEP0,   K_DEP2G, K_DEP2P )
!$OMP+PRIVATE( F_PBL,   FC,       Faq,      LWC )
!$OMP+PRIVATE( DT,      KMAX )
!$OMP+PRIVATE( A,       Xold,     Xnew,     s1, s2, s3, s4 )
!$OMP+PRIVATE( DEP_Hg0, DEP_HG2G, DEP_HG2P )
!$OMP+PRIVATE( DEP_HG2G_DRY,      DEP_HG2G_SALT,     DEP_HG2P_DRY )
!$OMP+PRIVATE( DEP_HG2_DRY )
!$OMP+PRIVATE( GROSS_OX, GROSS_OX_BR, GROSS_OX_OH, GROSS_OX_O3 )
!$OMP+PRIVATE( GROSS_RED, NET_OX )
!$OMP+PRIVATE( AREA_CM2,  DEP_DRY_FLX )
!$OMP+PRIVATE( F_HG0_DEP, Hg_Cat )
!%%%% NOTE: LOOP IS IN WRONG ORDER, WILL BE MORE EFFICIENT IF REVERSED %%%%
!%%%% Putting this loop in the correct order causes differences in the %%%%
!%%%% AD03_HG2_SS diagnostic array when using mp vs sp (mps, 9/3/15)   %%%%
      DO I=1, IIPAR
      DO J=1, JJPAR
      DO L=1, LLPAR

         !------------------------------------------------------
         ! Dry deposition frequencies, 1/s
         !------------------------------------------------------
      
         ! Initialize
         K_DRYD0  = 0e+0_fp
         K_DRYD2G = 0e+0_fp
         K_DRYD2P = 0e+0_fp

         !------------------------------------------------------
         ! Sea-salt aerosol uptake
         !------------------------------------------------------

         ! IS_WATER returns true for ocean boxes. 
         ! Fresh water may be TRUE in 
         ! principle but is not at 4x5 resolution.
         ! Sea salt will not be active in coastal boxes that have some ocean
         ! but where IS_WATER is FALSE.

         ! RGM uptake on sea-salt aerosol, 1/s
         ! (based on SSA production rate (wind speed) 
         ! and salinity (RH)) already zero over land
         K_SALT = HG2_SEASALT_LOSSRATE(I,J)
         
         ! Alternatively, uncomment the following lines
         ! to use a constant uptake rate (1/s) that is tuned to
         ! fit Okinawa RGM observations (Selin et al. 2007)
         !IF (IS_WATER(I,J)) THEN
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
#if   defined( MERRA ) || defined( GEOS_FP ) || defined( MERRA2 )
         ! Deposition occurs over areas that are not ocean, snow or ice
         ! (jaf, 4/26/11)
         F_Hg0_Dep = 1e+0_fp - MAX( 0e+0_fp, MIN( 1e+0_fp, 
     &        State_Met%FROCEAN(I,J) + 
     &        State_Met%FRSNO(I,J) + 
     &        State_Met%FRLANDIC(I,J) ) )
#elif defined( GEOS_5 )
         ! Deposition disabled for entire box if it is ocean, ice
         ! or land with more than 10 mm water-equivalent snow depth
         IF ( (State_Met%LWI(I,J) == 0) .OR. 
     &        (IS_ICE(I,J, State_Met)) .OR. 
     &        (IS_LAND(I,J, State_Met) .AND.
     &        State_Met%SNOMAS(I,J) > 10e+0_fp) ) THEN
            F_Hg0_Dep = 0e+0_fp
         ENDIF
#else
         ! Deposition disabled for entire box if it is ocean, ice
         ! or land with more than 10 mm water-equivalent snow depth
         IF ( (State_Met%LWI(I,J) == 0) .OR. 
     &        (IS_ICE(I,J, State_Met)) .OR. 
     &        (IS_LAND(I,J, State_Met) .AND. 
     &        State_Met%SNOW(I,J) > 10e+0_fp) ) THEN
            F_Hg0_Dep = 0e+0_fp
         ENDIF                 
#endif 

         !------------------------------------------------------
         ! Total deposition frequencies, 1/s
         !------------------------------------------------------

         ! Fraction of box (I,J,L) underneath the PBL top [dimensionless]
         F_PBL = GET_FRAC_UNDER_PBLTOP( I, J, L )

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
         C_O3        = GET_O3( I, J, L,        State_Met )
         C_OH        = GET_OH( I, J, L,        State_Met )
         C_BR        = GET_BR( I, J, L, C_BRO, State_Met )

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
            K_BR = GET_HGBR_RATE( C_BR, T(I,J,L), 
     &           State_Met%PMID_DRY(I,J,L), C_OH, METHOD )
         ELSE
            K_BR = 0e+0_fp
         ENDIF

         ! Oxidation rate by BrO
         IF (LBROCHEM) THEN
            K_BRO = K_HG_BRO * C_BRO
         ELSE
            K_BRO = 0e+0_fp
         ENDIF
         
         ! Total oxidation rate, 1/s
         K_OX = K_O3 + K_OH + K_BR + K_BRO

         !---------------------------------------------
         ! Reduction rate (requires liquid water content)
         !---------------------------------------------

#if   defined( GEOS_5 ) || defined( MERRA ) || defined( GEOS_FP ) || defined( MERRA2 )

         !---------------------------------------------
         ! GEOS-5/MERRA: Get LWC, FC from met fields
         ! (skim, bmy, 1/14/10)
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
         LWC     = SAFE_DIV( LWC, FC, 0e+0_fp )

#else
         !---------------------------------------------
         ! Otherwise, compute FC, LWC as before
         !---------------------------------------------    

         ! Get volume fraction of gridbox containing 
         ! cloud [unitless], following Sundqvist et al 1989.
         FC          = GET_VCLDF( I, J, L, State_Met  )

         ! Get liquid water content in cloud (m3/m3)
         LWC         = GET_LWC( T(I,J,L) ) 

         ! There should be no liquid water when T < 258K 
         IF ( T(I,J,L) < 258e+0_fp )  LWC = 0e+0_fp

#endif

         ! Define fraction of in-cloud Hg(II) which is in aqueous solution
         ! [dimensionless]
         Faq  = ( K0 * R_DRY * T(I,J,L) * LWC )
         Faq  = Faq / ( 1e+0_fp + Faq )

         ! Cl- in sea-salt aerosol enhances solubility 2000X in MBL
         IF (LRED_JNO2 .AND. (F_PBL >0.1) .AND. 
     &        IS_WATER(I,J, State_Met )) THEN
            Faq = ( K0 * 2e+3_fp * R_DRY * T(I,J,L) * LWC )
            Faq = Faq / ( 1e+0_fp + Faq )
         ENDIF
         
         ! Fraction of all Hg(II) in the grid box which is in cloud water.
         ! We have used the in-cloud LWC to compute the sulfate
         ! aqueous chemistry.  We get the appropriate grid-box averaged 
         ! mass of SO2 and sulfate by multiplying the reaction rates
         ! by the cloud fraction after the aqueous chemistry
         ! has been done.  (lzh, jaf, bmy, 5/27/11)
         ! Adopted for mercury chemistry (yzh, 11/9/2011)
         Faq =  Faq * FC

         ! Define K for the reduction reaction. 
         IF (LRED_JNO2) THEN
            K_RED = K_RED_JNO2 * Faq * GET_JNO2(I,J,L,State_Met)
         ELSE
            ! Include the fraction of 
            ! Hg(II) in air within the Kreduction & 
            ! scale to OH concentration [/s]
            K_RED = K_RED_OH * Faq * C_OH
         ENDIF
         
         !---------------------------------------------    
         ! Round off small chemical rates
         !---------------------------------------------    
         
         ! Round very small K_OX, K_RED to prevent numerical errors
         ! (jaf, cdh, 11/17/11)
         IF ( K_OX  < 1e-10_fp ) K_OX  = 0e+0_fp
         IF ( K_RED < 1e-10_fp ) K_RED = 0e+0_fp

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
         KMAX = K_OX + K_RED + K_DEP0 + K_DEP2G + K_DEP2P

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
         A(1,1) = -K_DEP0-K_OX
         A(1,2) = +K_RED
         A(1,3) = +K_RED

         ! Row 2: d(Hg2g)/dt = K_OX * Hg0  - (K_RED + K_DEP2G) * Hg2g
         A(2,1) = +K_OX
         A(2,2) = -K_DEP2G - K_RED

         ! Row 3: d(Hg2p)/dt = -(K_DEP2P + K_RED) * Hg2p 
         A(3,3) = -K_DEP2P - K_RED

         ! Row 4-6 accumulate dry dep of Hg0, Hg2g, Hg2p, respectively
         A(4,1) = K_DEP0
         A(5,2) = K_DEP2G
         A(6,3) = K_DEP2P

         ! Row 7 accumulate gross oxidation of Hg0
         A(7,1) = K_OX
            
         ! Row 8-9 accumulate gross reduction of Hg2g, Hg2p, respectively
         A(8,2) = K_RED
         A(9,3) = K_RED

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
               Xnew = Xold + DT/6e+0_fp * (s1 +2*s2 + 2*s3 + s4 )

               ! Make sure concentrations are not negative
               IF ( (Xnew(1) < 0e+0_fp) .or. (Xnew(2)<0) .or.
     &              (Xnew(3)<0) ) THEN
                  WRITE(6,101) I,J,L,K,NSTEPS, 
     &                 Xnew(1),Xnew(2),Xnew(3),
     &                 Xold(1),Xold(2),Xold(3),
     &                 K_OX, K_RED, K_DEP0, K_DEP2G, K_DEP2P
 101              FORMAT( 'Negative Hg in Box (I,J,L):', 3I3,
     &                 ', step',I3,' out of ',I3,
     &                 ', Initial conc: ',3E10.3,
     &                 ', Current conc: ',3E10.3,
     &                 ', Rate coeffs: ',6E10.3 )
                  CALL ERROR_STOP( 'Negative Hg concentration',
     &                 ' MERCURY_MOD: CHEMMERCURY' )
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

            ! Apportion gross oxidation between O3 and OH [kg]
            IF ( (K_OX   < SMALLNUM) .OR. 
     &           (GROSS_OX < SMALLNUM) ) THEN
               GROSS_OX_OH = 0e+0_fp
               GROSS_OX_BR = 0e+0_fp
               GROSS_OX_O3 = 0e+0_fp
            ELSE
               GROSS_OX_OH = GROSS_OX * K_OH / K_OX
               GROSS_OX_BR = GROSS_OX * K_BR / K_OX
               GROSS_OX_O3 = GROSS_OX * K_O3 / K_OX
            ENDIF
               
            ! Apportion deposition between dry deposition and sea salt [kg]
            IF ( (K_DEP2G    < SMALLNUM) .OR. 
     &           (DEP_Hg2g  < SMALLNUM) ) THEN
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
            IF ( LDYNOCEAN ) THEN

               ! Archive dry deposited Hg2 for Hg category # N
               CALL ADD_Hg2_DD( I, J, N, DEP_Hg2G )

               ! Archive dry-deposited HgP for Hg category # N
               CALL ADD_HgP_DD( I, J, N, DEP_Hg2P )

            ENDIF           

            ! Add deposited Hg(II) to the snowpack
            IF ( LHGSNOW ) THEN
               DEP_HG2_DRY =  DEP_HG2g_DRY + DEP_HG2p_DRY
               CALL ADD_HG2_SNOWPACK( I, J, N, DEP_HG2_DRY, 
     &                                         State_Met    )
            ENDIF           
!$OMP END CRITICAL

#if defined( BPCH_DIAG )
            !=================================================================
            ! ND44 diagnostic: drydep flux of Hg(II) [molec/cm2/s]
            !=================================================================
            IF ( ( ND44 > 0 .OR. LGTMM ) .AND. (.NOT. LNLPBL) ) THEN

               ! Grid box surface area [cm2]
               AREA_CM2 = GET_AREA_CM2( I, J, L )

               ! Amt of Hg(0) lost to drydep [molec/cm2/s]
               DEP_DRY_FLX  = DEP_HG0 * AVO / ( 1.e-3_fp *
     &             State_Chm%SpcData(Hg0_Id_List(N))%Info%emMW_g ) / 
     &             ( AREA_CM2 * DTCHEM )

               ! Archive Hg(0) drydep flux in AD44 array [molec/cm2/s]
               AD44(I,J,Hg0_Id_List(N),1) = 
     &         AD44(I,J,Hg0_Id_List(N),1) +  DEP_DRY_FLX

               ! Amt of Hg(II)g lost to drydep [molec/cm2/s]
               DEP_DRY_FLX  = DEP_HG2G_DRY * AVO / ( 1.e-3_fp *
     &             State_Chm%SpcData(Hg2_Id_List(N))%Info%emMW_g ) / 
     &             ( AREA_CM2 * DTCHEM )

               ! Archive Hg(II) drydep flux in AD44 array [molec/cm2/s]
               AD44(I,J,Hg2_Id_List(N),1) = 
     &         AD44(I,J,Hg2_Id_List(N),1) + DEP_DRY_FLX

               ! Amt of Hg(II)p lost to drydep [molec/cm2/s]
               DEP_DRY_FLX  = DEP_HG2P_DRY * AVO / ( 1.e-3_fp *
     &             State_Chm%SpcData(HgP_Id_List(N))%Info%emMW_g ) / 
     &             ( AREA_CM2 * DTCHEM )

               ! Archive Hg(II)p drydep flux in AD44 array [molec/cm2/s]
               AD44(I,J,HgP_Id_List(N),1) = 
     &         AD44(I,J,HgP_Id_List(N),1) + DEP_DRY_FLX

            ENDIF
#endif

            !==============================================================
            ! ND03 diagnostic: Hg(II) production [kg]
            !==============================================================
            IF ( ND03 > 0 .AND. L <= LD03 ) THEN 

               ! Store chemistry diagnostics only for total species
               IF ( Hg0_Id_List(N) == ID_HG_TOT) THEN

                  ! Net oxidation [kg]
                  NET_OX = GROSS_OX - GROSS_RED

                  AD03_Hg2_Hg0(I,J,L,N)= AD03_Hg2_Hg0(I,J,L,N) +
     &                 NET_OX
                  AD03_Hg2_Br(I,J,L,N) = AD03_Hg2_Br(I,J,L,N)  + 
     &                 GROSS_OX_BR 
                  AD03_Hg2_OH(I,J,L,N) = AD03_Hg2_OH(I,J,L,N)  + 
     &                 GROSS_OX_OH
                  AD03_Hg2_O3(I,J,L,N) = AD03_Hg2_O3(I,J,L,N)  + 
     &                 GROSS_OX_O3

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

         ENDDO

      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! Partition Hg(II) between gas and aerosol
      CALL PARTITIONHG2( am_I_Root, Input_Opt, State_Chm, RC )

      ! Free pointer memory
      Spc     => NULL()
      T       => NULL()
      SpcInfo => NULL()

      IF ( LPRT ) CALL DEBUG_MSG( 'CHEMMERCURY: a CHEM_HgP' )

      END SUBROUTINE CHEMMERCURY
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
      FUNCTION GET_HGBR_RATE( BR, T, P, OH, METHOD ) RESULT( K_HGBR ) 
!
! !INPUT PARAMETERS: 
!
      ! BR     : Concentration of atomic Br [molec/cm3] 
      ! T      : Temperature [K]
      ! P      : Pressure [hPa]
      ! OH     : Concentration of OH [molec/cm3]
      ! METHOD : Set of rate constants to use in calculation
      REAL(fp),         INTENT(IN)  :: BR, OH, T, P
      CHARACTER(LEN=*), INTENT(IN)  :: METHOD 
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
!  05 Jan 2016 - E. Lundgren - Use global physical parameters
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:

      REAL(fp)            :: k1,   k2, k2a, k3br, k3oh
      REAL(fp)            :: k1G,  k2G
      REAL(fp)            :: Nair, N0
!
! !DEFINED PARAMETERS:
!

      !=================================================================
      ! GET_HGBR_RATE begins here!
      !=================================================================

      !number density of air [molec/cm3]
      Nair = P * 1e+2_fp / (AIRMW * Rd * T ) * AVO / 1e+3_fp

      !standard air density at STP: 1 atm, 273K [molec/cm3]
      N0 = ATM / (AIRMW * Rd * 273e+0_fp ) * AVO / 1e+3_fp

      SELECT CASE( METHOD )
      CASE( 'Goodsite' )
         !All rates from Goodsite et al. 2004
         !No HgBr+OH reaction
         k1   = 1.1e-12_fp * ( T / 298e+0_fp ) ** ( -2.37e+0_fp ) * 
     &          ( Nair / N0 )
         k2   = 1.2e+10_fp  * exp( -8357e+0_fp / T )
         k2a  = 0e+0_fp
         k3br = 2.5e-10_fp * ( T / 298e+0_fp ) ** ( -.57e+0_fp ) 
         k3oh = 0e+0_fp
      
      CASE( 'GoodsiteY' )
         !All rates from Goodsite et al. 2004
         !Include HgBr+OH reaction
         k1   = 1.1e-12_fp * ( T / 298e+0_fp ) ** ( -2.37e+0_fp ) * 
     &          ( Nair / N0 )
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
         k1G  = 1.1e-12_fp * ( T / 298e+0_fp ) ** ( -2.37e+0_fp ) * 
     &          ( Nair / N0 )
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
         k1   = 1.1e-12_fp * ( T / 298e+0_fp ) ** ( -2.37e+0_fp ) * 
     &          ( Nair / N0 )
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
         k1G  = 1.1e-12_fp * ( T / 298e+0_fp ) ** ( -2.37e+0_fp ) * 
     &          ( Nair / N0 )
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
         k1G  = 1.1e-12_fp * ( T / 298e+0_fp ) ** ( -2.37e+0_fp ) * 
     &          ( Nair / N0 )
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
         k1G  = 1.1e-12_fp * ( T / 298e+0_fp ) ** ( -2.37e+0_fp ) * 
     &          ( Nair / N0 )
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
         k1G  = 3.7e-13_fp * ( T / 298e+0_fp ) ** ( -2.76e+0_fp ) * 
     &          ( Nair / N0 )
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
         K_HGBR = k1 * BR * (k3br * BR + k3oh * OH ) / 
     &          ( k2 + k2a * BR + k3br * BR + k3oh * OH )  

      ELSE
         
         ! Avoid divide by zero in rate
         K_HGBR = 0
      
      ENDIF

      ! Return to calling program
      END FUNCTION GET_HGBR_RATE
!EOP
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
      SUBROUTINE EMISSMERCURY( am_I_Root, Input_Opt,
     &                         State_Met, State_Chm, RC )
!
! !USES:
!
      USE CMN_SIZE_MOD
      USE DEPO_MERCURY_MOD,   ONLY : RESET_HG_DEP_ARRAYS
      USE ErrCode_Mod
      USE ERROR_MOD
      USE GC_GRID_MOD,        ONLY : GET_XMID, GET_YMID
      USE Input_Opt_Mod,      ONLY : OptInput
      USE LAND_MERCURY_MOD,   ONLY : LAND_MERCURY_FLUX, VEGEMIS
      USE LAND_MERCURY_MOD,   ONLY : SOILEMIS, BIOMASSHG
      USE LAND_MERCURY_MOD,   ONLY : SNOWPACK_MERCURY_FLUX
      USE OCEAN_MERCURY_MOD,  ONLY : OCEAN_MERCURY_FLUX
      USE State_Chm_Mod,      ONLY : ChmState
      USE State_Met_Mod,      ONLY : MetState
      USE TIME_MOD,           ONLY : GET_MONTH, ITS_A_NEW_MONTH
      USE UnitConv_Mod

      ! Added for GTMM (ccc, 11/19/09)
      USE LAND_MERCURY_MOD,   ONLY : GTMM_DR
!
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
! !REMARKS:
! 
! 
! !REVISION HISTORY: 
!  03 Jun 2013 - N. (Eckley) Selin - Initial version
!  (1 ) Now call OCEAN_MERCURY_FLUX from "ocean_mercury_mod.f" to compute 
!        the emissions of Hg0 from the ocean instead of reading it from disk.
!        (sas, bmy, 1/20/05)
!  (2 ) Now no longer call COMPUTE_FEMIS, since we can get the same information
!        from routine GET_FRAC_OF_PBL in "pbl_mix_mod.f" (bmy, 2/22/05)
!  (3 ) Now modified for new ocean mercury module. (cdh, sas, bmy, 4/6/06)
!  01 Mar 2012 - R. Yantosca - Now use GET_XMID(I,J,L) from grid_mod.F90
!  01 Mar 2012 - R. Yantosca - Now use GET_YMID(I,J,L) from grid_mod.F90
!  25 Mar 2013 - R. Yantosca - Now accept am_I_Root, Input_Opt, State_Chm, RC
!  08 Dec 2009 - R. Yantosca - Added ProTeX headers
!  20 Jun 2014 - R. Yantosca - Now pass Input_Opt to VEGEMIS
!  23 Jun 2014 - R. Yantosca - Now pass am_I_Root to OCEAN_MERCURY_FLUX
!  05 Oct 2015 - E. Lundgren - Now convert tracer units to kg locally
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL, SAVE   :: FIRST = .TRUE. 
      INTEGER         :: THISMONTH, I, J

      ! For fields from Input_Opt
      LOGICAL         :: LDYNOCEAN
      LOGICAL         :: LGTMM
      LOGICAL         :: LNLPBL
      LOGICAL         :: LPRT
      LOGICAL         :: LPREINDHG
      LOGICAL         :: LEMIS

      !=================================================================
      ! EMISSMERCURY begins here!
      !=================================================================

      ! Assume success
      RC        = GC_SUCCESS
      
      ! Copy fields from Input_Opt
      LDYNOCEAN = Input_Opt%LDYNOCEAN
      LGTMM     = Input_Opt%LGTMM
      LNLPBL    = Input_Opt%LNLPBL
      LPREINDHG = Input_Opt%LPREINDHG
      LEMIS     = Input_Opt%LEMIS
      LPRT      = ( Input_Opt%LPRT .and. am_I_Root )

      ! Convert species units to [kg] for EMISSMERCURY (ewl, 8/12/15)
      ! NOTE: This is temporary until mercury code is in [kg/kg dry air]
      CALL ConvertSpc_KgKgDry_to_Kg( am_I_Root, State_Met,  
     &                               State_Chm, RC ) 
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error('Unit conversion error', RC, 
     &                 'Routine EMISSIONS_RUN in emissions_mod.F')
         RETURN
      ENDIF                         

      ! First-time initialization
      IF ( FIRST ) THEN

         ! Check that emissions are turned on. Print error message 
         ! and stop GEOS-Chem if emissions are turned off (ewl, 9/1/15)
         IF ( .not. LEMIS ) THEN
            PRINT *, 'ERROR: Hg emissions are need for simulation ' //
     &               'but emissions are turned off in input.geos'
            CALL GEOS_CHEM_STOP
         ENDIF

         ! Read anthro, ocean, land emissions of Hg from disk
         CALL MERCURY_READYR( am_I_Root, Input_Opt, RC )

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

         CALL OCEAN_MERCURY_FLUX( am_I_Root, Input_Opt, State_Met, 
     &                            State_Chm, EHg0_oc,   RC         )
     &                            
         IF ( LPRT ) CALL DEBUG_MSG( '### EMISSMERCURY: a OCEAN_FLUX' )

      ENDIF

      IF ( LGTMM ) THEN
         !--------------------------------------------------------------
         ! Here we are using the Global Terrstrial Mercury Model ...
         !--------------------------------------------------------------
         IF ( ITS_A_NEW_MONTH() ) THEN

            ! OLD CODE: CALL GTMM_DR( EHg0_gtm(:,:), State_Met )
            CALL GTMM_DR
     &           ( am_I_Root, Input_Opt, State_Met, EHg0_gtm, RC )

            IF ( LPRT ) CALL DEBUG_MSG( '### EMISSMERCURY: a GTMM' )
         ENDIF
  
      ELSE
         !--------------------------------------------------------------
         ! Here we are NOT using the Global Terrstrial Mercury Model ...
         !--------------------------------------------------------------
         CALL LAND_MERCURY_FLUX ( EHg0_ln, LHGSNOW, State_Met )
         IF ( LPRT ) CALL DEBUG_MSG( '### EMISSMERCURY: a LAND_FLUX' )
      
         ! Bug fix: VEGEMIS shouldn't be tied to GCAP emissions
         ! (jaf, eds, 4/1/11)
         CALL VEGEMIS( am_I_Root, Input_Opt, State_Met,
     &                 LVEGEMIS,  EHg0_dist, EHg0_vg,   RC )
         IF ( LPRT ) CALL DEBUG_MSG( '### EMISSMERCURY: a VEGEMIS' )
      
         CALL SOILEMIS( EHg0_dist, EHg0_so, State_Met )
         IF ( LPRT ) CALL DEBUG_MSG( '### EMISSMERCURY: a SOILEMIS' )

      ENDIF

      CALL SNOWPACK_MERCURY_FLUX ( EHg0_snow, LHGSNOW, State_Met )
      IF ( LPRT ) CALL DEBUG_MSG( '### EMISSMERCURY: a SNOW_FLUX' )
      
      CALL BIOMASSHG( am_I_Root, Input_Opt, EHg0_bb, RC )
      IF ( LPRT ) CALL DEBUG_MSG( '### EMISSMERCURY: a BIOMASS' )

      IF ( LnoUSAemis ) THEN

         ! No Anthropogenic emissions over USA (cdh, 5/6/2009)
         DO J = 1, JJPAR
         DO I = 1, IIPAR

!            IF ( GET_XMID( I, J, 1 ) > -125 .AND. 
!     &           GET_XMID( I, J, 1 ) < -65  .AND.
!     &           GET_YMID( I, J, 1 ) >  25  .AND. 
!     &           GET_YMID( I, J, 1 ) < 50 ) THEN
            IF ( GET_XMID( I, J, 1 ) > -140 .AND.
     &           GET_XMID( I, J, 1 ) < -40  .AND.
     &           GET_YMID( I, J, 1 ) >  10  .AND.
     &           GET_YMID( I, J, 1 ) <  70 ) THEN
               
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
      IF ( LPRT ) CALL DEBUG_MSG( '### EMISSMERCURY: ' //
     &     'a RESET_HG_DEP_ARRAYS' )
      
      ! If we are using the non-local PBL mixing,
      ! we need to initialize the EMIS_SAVE array (cdh, 08/27/09)
      ! EMIS_SAVE is now HG_EMIS (ckeller, 10/21/2014)
      IF ( LNLPBL ) HG_EMIS = 0.0e+0_fp

      ! Add Hg(0) source into State_Chm%Species [kg]
      CALL SRCHg0( am_I_Root, Input_Opt, State_Chm, RC )
      IF ( LPRT ) CALL DEBUG_MSG( '### EMISSMERCURY: a SRCHg0' )

      ! Add Hg(II) source into State_Chm%Species [kg]
      CALL SRCHg2( am_I_Root, Input_Opt, State_Chm, RC )
      IF ( LPRT ) CALL DEBUG_MSG( '### EMISSMERCURY: a SRCHg2' )

      ! Add HgP source into State_Chm%Species [kg]
      CALL SRCHgP( am_I_Root, Input_Opt, State_Chm, RC )
      IF ( LPRT ) CALL DEBUG_MSG( '### EMISSMERCURY: a SRCHgP' )
      
      ! Convert species units back to [kg/kg] after mercury emissions
      !  (ewl, 8/12/15)
      CALL ConvertSpc_Kg_to_KgKgDry( am_I_Root, State_Met,  
     &                               State_Chm, RC ) 
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error('Unit conversion error', RC, 
     &                 'Routine EMISSMERCURY in mercury_mod.F')
         RETURN
      ENDIF  

      ! Return to calling program 
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
      SUBROUTINE EMITHG( I, J, L, ID, E_HG, Input_Opt, State_Chm )
!
! !USES:
!
      USE ErrCode_Mod
      USE GC_GRID_MOD,        ONLY : GET_AREA_M2
      USE Input_Opt_Mod,      ONLY : OptInput
      USE State_Chm_Mod,      ONLY : ChmState
      USE TIME_MOD,           ONLY : GET_TS_EMIS
!
! !INPUT PARAMETERS:
!
      INTEGER,        INTENT(IN)    :: I, J, L, ID ! Grid boxes + species #
      REAL(fp),       INTENT(IN)    :: E_Hg        ! Hg emissions 
      TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
!
! !INPUT/OUTPUT PARAMETERS:
!
      TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
! 
! !REVISION HISTORY: 
!  27 Aug 2009 - C. Holmes   - Initial version
!  25 Mar 2013 - R. Yantosca - Now accept Input_Opt, State_Chm, arguments
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

         ! Save emissions for non-local PBL mixing or emit directly.
         ! Make sure that emitted mass is non-negative
         ! This is hear only for consistency with old code which warned
         ! of underflow error (cdh, 08/27/09)
         ! EMIS_SAVE is now HG_EMIS array. Convert kg to kg/m2/s
         ! (ckeller, 09/23/2014)
         !--------------------------------------------------------------

         ! Surface area [m2]
         AM2             = GET_AREA_M2(I,J,1)

         ! Emission timestep
         TS              = GET_TS_EMIS() * 60.0e+0_fp

         ! Save emissions as [kg/m2/s].  These will be added
         ! to the chemical species array in routine DO_TEND
         ! (in mixing_mod.F90).
         HG_EMIS(I,J,ID) = HG_EMIS(I,J,ID) 
     &                   + ( MAX(E_HG,0e+0_fp) / AM2 / TS )

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
      SUBROUTINE SRCHg0( am_I_Root, Input_Opt, State_Chm, RC )
!
! !USES:
!
      USE CMN_SIZE_MOD
      USE DIAG03_MOD,         ONLY : AD03, ND03, AD03_nat
      USE ErrCode_Mod
      USE ERROR_MOD,          ONLY : ERROR_STOP
      USE Input_Opt_Mod,      ONLY : OptInput
      USE PBL_MIX_MOD,        ONLY : GET_FRAC_OF_PBL, GET_PBL_MAX_L
      USE State_Chm_Mod,      ONLY : ChmState
      USE TIME_MOD,           ONLY : GET_TS_EMIS
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)    :: am_I_Root   ! Are we on the root CPU?
      TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
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
!  21 Jan 2005 - N. (Eckley) Selin, C. Holmes - Initial version
!  (1 ) Now use diagnostic arrays from "diag03_mod.f" (bmy, 1/21/05)
!  (2 ) Now references GET_FRAC_OF_PBL and GET_PBL_MAX_L from "pbl_mix_mod.f".
!        Remove reference to FEMIS. (bmy, 2/22/05)
!  (3 ) EHg0_an is now a 2-D array.  Modified for new ocean mercury module.
!        Now use ID_Hg0 index array from "tracerid_mod.f".  Now make sure
!        STT does not underflow. (cdh, bmy, 4/6/06)
!  25 Mar 2013 - R. Yantosca - Now accept am_I_Root, Input_Opt, State_Chm, RC
!  30 Jun 2016 - R. Yantosca - Remove instances of STT.  Now get the advected
!                              species ID from State_Chm%Map_Advect.
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
      DTSRCE    = GET_TS_EMIS() * 60e+0_fp

      ! Maximum extent of the PBL [model levels]
      PBL_MAX   = GET_PBL_MAX_L() 

      ! Loop over grid boxes
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, N, T_Hg_An, T_Hg, F_OF_PBL, E_Hg, NN)
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         IF ( LPREINDHG ) THEN  !eds

            ! Anthropogenic emissions
            T_Hg_An = 0e+0_fp

            ! No biomass burning emissions
            EHg0_bb(I,J) = 0e+0_fp
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
            T_Hg = T_Hg_An +
     &             EHg0_bb(I,J) +
     &             EHg0_oc(I,J,ID_Hg_tot) + 
     &             EHg0_geo(I,J) +
     &             EHg0_gtm(I,J) +
     &             EHg0_snow(I,J,ID_Hg_tot)
         ELSE
            T_Hg = T_Hg_An +
     &             EHg0_bb(I,J) +
     &             EHg0_oc(I,J,ID_Hg_tot) + 
     &             EHg0_ln(I,J,ID_Hg_tot) +
     &             EHg0_geo(I,J) +
     &             EHg0_vg(I,J) +
     &             EHg0_so(I,J) +
     &             EHg0_snow(I,J,ID_Hg_tot)
         ENDIF

         !==============================================================
         ! Partition Hg0 throughout PBL; store into STT [kg]
         ! Now make sure STT does not underflow (cdh, bmy, 4/6/06)
         !==============================================================

         ! Loop up to max PBL level
         DO L = 1, PBL_MAX

            ! Fraction of box (I,J,L) w/in the PBL [unitless]
            F_OF_PBL = GET_FRAC_OF_PBL( I, J, L )

            !-----------------
            ! Total Hg species
            !-----------------
            N    = Hg0_Id_List(ID_Hg_tot)
            E_Hg = F_OF_PBL * T_Hg * DTSRCE
            CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm )

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
               CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm )

               ! Anthro United States Hg0
               N    = Hg0_Id_List(ID_Hg_usa)
               E_Hg = F_OF_PBL * EHg0_usa(I,J) * DTSRCE
               CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm )

               ! Anthro Central America Hg0
               N    = Hg0_Id_List(ID_Hg_cam)
               E_Hg = F_OF_PBL * EHg0_cam(I,J) * DTSRCE
               CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm )

               ! Anthro South America Hg0
               N    = Hg0_Id_List(ID_Hg_sam)
               E_Hg = F_OF_PBL * EHg0_sam(I,J) * DTSRCE
               CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm )

               ! Anthro West Africa Hg0
               N    = Hg0_Id_List(ID_Hg_waf)
               E_Hg = F_OF_PBL * EHg0_waf(I,J) * DTSRCE
               CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm )

               ! Anthro East Africa Hg0
               N    = Hg0_Id_List(ID_Hg_eaf)
               E_Hg = F_OF_PBL * EHg0_eaf(I,J) * DTSRCE
               CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm )

               ! Anthro South Africa Hg0
               N    = Hg0_Id_List(ID_Hg_saf)
               E_Hg = F_OF_PBL * EHg0_saf(I,J) * DTSRCE
               CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm )

               ! Anthro North Africa Hg0
               N    = Hg0_Id_List(ID_Hg_naf)
               E_Hg = F_OF_PBL * EHg0_naf(I,J) * DTSRCE
               CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm )

               ! Anthro OECD Europe Hg0
               N    = Hg0_Id_List(ID_Hg_eur)
               E_Hg = F_OF_PBL * EHg0_eur(I,J) * DTSRCE
               CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm )

               ! Anthro Eastern Europe Hg0
               N    = Hg0_Id_List(ID_Hg_eeu)
               E_Hg = F_OF_PBL * EHg0_eeu(I,J) * DTSRCE
               CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm )

               ! Anthro Middle East Hg0
               N    = Hg0_Id_List(ID_Hg_mde)
               E_Hg = F_OF_PBL * EHg0_mde(I,J) * DTSRCE
               CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm )

               ! Anthro Former USSR Hg0
               N    = Hg0_Id_List(ID_Hg_sov)
               E_Hg = F_OF_PBL * EHg0_sov(I,J) * DTSRCE
               CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm )

               ! Anthro South Asia Hg0
               N    = Hg0_Id_List(ID_Hg_sas)
               E_Hg = F_OF_PBL * EHg0_sas(I,J) * DTSRCE
               CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm )

               ! Anthro East Asia Hg0
               N    = Hg0_Id_List(ID_Hg_eas)
               E_Hg = F_OF_PBL * EHg0_eas(I,J) * DTSRCE
               CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm )

               ! Anthro Southeast Asia Hg0
               N    = Hg0_Id_List(ID_Hg_sea)
               E_Hg = F_OF_PBL * EHg0_sea(I,J) * DTSRCE
               CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm )

               ! Anthro Japan Hg0
               N    = Hg0_Id_List(ID_Hg_jpn)
               E_Hg = F_OF_PBL * EHg0_jpn(I,J) * DTSRCE
               CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm )

               ! Anthro Oceania Hg0
               N    = Hg0_Id_List(ID_Hg_oce)
               E_Hg = F_OF_PBL * EHg0_oce(I,J) * DTSRCE
               CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm )

               ! Soil emissions
               N    = Hg0_Id_List(ID_Hg_so)
               E_Hg = F_OF_PBL * EHg0_so(I,J) * DTSRCE
               CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm )

               ! Biomass burning emissions             
               N    = Hg0_Id_List(ID_Hg_bb)
               E_Hg = F_OF_PBL * EHg0_bb(I,J) * DTSRCE
               CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm )

               ! Geogenic emissions             
               N    = Hg0_Id_List(ID_Hg_geo)
               E_Hg = F_OF_PBL * EHg0_geo(I,J) * DTSRCE
               CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm )

               !--------------------
               ! Ocean re-emissions
               !--------------------

               ! Anthro Canada re-emission from ocean
               N    = Hg0_Id_List(ID_Hg_can)
               E_Hg = F_OF_PBL * EHg0_oc(I,J,ID_Hg_can) * DTSRCE
               CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm )

               ! Anthro United States re-emission from ocean 
               N    = Hg0_Id_List(ID_Hg_usa)
               E_Hg = F_OF_PBL * EHg0_oc(I,J,ID_Hg_usa) * DTSRCE
               CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm )

               ! Anthro Central America re-emission from ocean 
               N    = Hg0_Id_List(ID_Hg_cam)
               E_Hg = F_OF_PBL * EHg0_oc(I,J,ID_Hg_cam) * DTSRCE
               CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm )

               ! Anthro South America re-emission from ocean 
               N    = Hg0_Id_List(ID_Hg_sam)
               E_Hg = F_OF_PBL * EHg0_oc(I,J,ID_Hg_sam) * DTSRCE
               CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm )

               ! Anthro West Africa re-emission from ocean 
               N    = Hg0_Id_List(ID_Hg_waf)
               E_Hg = F_OF_PBL * EHg0_oc(I,J,ID_Hg_waf) * DTSRCE
               CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm )

               ! Anthro East Africa re-emission from ocean 
               N    = Hg0_Id_List(ID_Hg_eaf)
               E_Hg = F_OF_PBL * EHg0_oc(I,J,ID_Hg_eaf) * DTSRCE
               CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm )

               ! Anthro South Africa re-emission from ocean 
               N    = Hg0_Id_List(ID_Hg_saf)
               E_Hg = F_OF_PBL * EHg0_oc(I,J,ID_Hg_saf) * DTSRCE
               CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm )

               ! Anthro North Africa re-emission from ocean 
               N    = Hg0_Id_List(ID_Hg_naf)
               E_Hg = F_OF_PBL * EHg0_oc(I,J,ID_Hg_naf) * DTSRCE
               CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm )

               ! Anthro OECD Europe re-emission from ocean 
               N    = Hg0_Id_List(ID_Hg_eur)
               E_Hg = F_OF_PBL * EHg0_oc(I,J,ID_Hg_eur) * DTSRCE
               CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm )

               ! Anthro Eastern Europe re-emission from ocean 
               N    = Hg0_Id_List(ID_Hg_eeu)
               E_Hg = F_OF_PBL * EHg0_oc(I,J,ID_Hg_eeu) * DTSRCE
               CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm )

               ! Anthro Middle East re-emission from ocean 
               N    = Hg0_Id_List(ID_Hg_mde)
               E_Hg = F_OF_PBL * EHg0_oc(I,J,ID_Hg_mde) * DTSRCE
               CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm )

               ! Anthro Former Soviet Union re-emission from ocean 
               N    = Hg0_Id_List(ID_Hg_sov)
               E_Hg = F_OF_PBL * EHg0_oc(I,J,ID_Hg_sov) * DTSRCE
               CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm )

               ! Anthro South Asia re-emission from ocean 
               N    = Hg0_Id_List(ID_Hg_sas)
               E_Hg = F_OF_PBL * EHg0_oc(I,J,ID_Hg_sas) * DTSRCE
               CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm )

               ! Anthro East Asia re-emission from ocean 
               N    = Hg0_Id_List(ID_Hg_eas)
               E_Hg = F_OF_PBL * EHg0_oc(I,J,ID_Hg_eas) * DTSRCE
               CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm )

               ! Anthro Southeast Asia re-emission from ocean 
               N    = Hg0_Id_List(ID_Hg_sea)
               E_Hg = F_OF_PBL * EHg0_oc(I,J,ID_Hg_sea) * DTSRCE
               CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm )

               ! Anthro Japan re-emission from ocean 
               N    = Hg0_Id_List(ID_Hg_jpn)
               E_Hg = F_OF_PBL * EHg0_oc(I,J,ID_Hg_jpn) * DTSRCE
               CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm )

               ! Anthro Oceania re-emission from ocean 
               N    = Hg0_Id_List(ID_Hg_oce)
               E_Hg = F_OF_PBL * EHg0_oc(I,J,ID_Hg_oce) * DTSRCE
               CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm )

               ! Soil re-emission from ocean 
               N    = Hg0_Id_List(ID_Hg_so)
               E_Hg = F_OF_PBL * EHg0_oc(I,J,ID_Hg_so) * DTSRCE
               CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm )

               ! Biomass burning re-emission from ocean 
               N    = Hg0_Id_List(ID_Hg_bb)
               E_Hg = F_OF_PBL * EHg0_oc(I,J,ID_Hg_bb) * DTSRCE
               CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm )

               ! Geogenic re-emission from ocean 
               N    = Hg0_Id_List(ID_Hg_geo)
               E_Hg = F_OF_PBL * EHg0_oc(I,J,ID_Hg_geo) * DTSRCE
               CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm )

               ! Below Mixed Layer Mid-Atlantic re-emission from ocean 
               N    = Hg0_Id_List(ID_Hg_atl)
               E_Hg = F_OF_PBL * EHg0_oc(I,J,ID_Hg_atl) * DTSRCE
               CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm )

               ! Below Mixed Layer North Atlantic re-emission from ocean 
               N    = Hg0_Id_List(ID_Hg_nat)
               E_Hg = F_OF_PBL * EHg0_oc(I,J,ID_Hg_nat) * DTSRCE
               CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm )

               ! Below Mixed Layer South Atlantic re-emission from ocean 
               N    = Hg0_Id_List(ID_Hg_sat)
               E_Hg = F_OF_PBL * EHg0_oc(I,J,ID_Hg_sat) * DTSRCE
               CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm )

               ! Below Mixed Layer North Pacific re-emission from ocean
               N    = Hg0_Id_List(ID_Hg_npa)
               E_Hg = F_OF_PBL * EHg0_oc(I,J,ID_Hg_npa) * DTSRCE
               CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm )

               ! Below Mixed Layer Arctic re-emission from ocean 
               N    = Hg0_Id_List(ID_Hg_arc)
               E_Hg = F_OF_PBL * EHg0_oc(I,J,ID_Hg_arc) * DTSRCE
               CALL EMITHG( I, J, L, N, E_Hg, Input_Opt,  State_Chm )

               ! Below Mixed Layer Antarctic re-emission from ocean 
               N    = Hg0_Id_List(ID_Hg_ant)
               E_Hg = F_OF_PBL * EHg0_oc(I,J,ID_Hg_ant) * DTSRCE
               CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm )

               ! Below Mixed Layer All Other Ocean re-emission from ocean
               N    = Hg0_Id_List(ID_Hg_ocn)
               E_Hg = F_OF_PBL * EHg0_oc(I,J,ID_Hg_ocn) * DTSRCE
               CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm )

               ! Residual Stratosphere from Spin-up re-emission from ocean 
               N    = Hg0_Id_List(ID_Hg_str)
               E_Hg = F_OF_PBL * EHg0_oc(I,J,ID_Hg_str) * DTSRCE
               CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm )

               !--------------------
               ! Land re-emissions
               !--------------------

               ! Anthro Canada re-emission from land
               N    = Hg0_Id_List(ID_Hg_can)
               E_Hg = F_OF_PBL * EHg0_ln(I,J,ID_Hg_can) * DTSRCE
               CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm )

               ! Anthro United States re-emission from land 
               N    = Hg0_Id_List(ID_Hg_usa)
               E_Hg = F_OF_PBL * EHg0_ln(I,J,ID_Hg_usa) * DTSRCE
               CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm )

               ! Anthro Central America re-emission from land 
               N    = Hg0_Id_List(ID_Hg_cam)
               E_Hg = F_OF_PBL * EHg0_ln(I,J,ID_Hg_cam) * DTSRCE
               CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm )

               ! Anthro South America re-emission from land 
               N    = Hg0_Id_List(ID_Hg_sam)
               E_Hg = F_OF_PBL * EHg0_ln(I,J,ID_Hg_sam) * DTSRCE
               CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm )

               ! Anthro West Africa re-emission from land 
               N    = Hg0_Id_List(ID_Hg_waf)
               E_Hg = F_OF_PBL * EHg0_ln(I,J,ID_Hg_waf) * DTSRCE
               CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm )

               ! Anthro East Africa re-emission from land 
               N    = Hg0_Id_List(ID_Hg_eaf)
               E_Hg = F_OF_PBL * EHg0_ln(I,J,ID_Hg_eaf) * DTSRCE
               CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm )

               ! Anthro South Africa re-emission from land 
               N    = Hg0_Id_List(ID_Hg_saf)
               E_Hg = F_OF_PBL * EHg0_ln(I,J,ID_Hg_saf) * DTSRCE
               CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm )

               ! Anthro North Africa re-emission from land 
               N    = Hg0_Id_List(ID_Hg_naf)
               E_Hg = F_OF_PBL * EHg0_ln(I,J,ID_Hg_naf) * DTSRCE
               CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm )

               ! Anthro OECD Europe re-emission from land 
               N    = Hg0_Id_List(ID_Hg_eur)
               E_Hg = F_OF_PBL * EHg0_ln(I,J,ID_Hg_eur) * DTSRCE
               CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm )

               ! Anthro Eastern Europe re-emission from land 
               N    = Hg0_Id_List(ID_Hg_eeu)
               E_Hg = F_OF_PBL * EHg0_ln(I,J,ID_Hg_eeu) * DTSRCE
               CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm )

               ! Anthro Middle East re-emission from land 
               N    = Hg0_Id_List(ID_Hg_mde)
               E_Hg = F_OF_PBL * EHg0_ln(I,J,ID_Hg_mde) * DTSRCE
               CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm )

               ! Anthro Former Soviet Union re-emission from land 
               N    = Hg0_Id_List(ID_Hg_sov)
               E_Hg = F_OF_PBL * EHg0_ln(I,J,ID_Hg_sov) * DTSRCE
               CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm )

               ! Anthro South Asia re-emission from land 
               N    = Hg0_Id_List(ID_Hg_sas)
               E_Hg = F_OF_PBL * EHg0_ln(I,J,ID_Hg_sas) * DTSRCE
               CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm )

               ! Anthro East Asia re-emission from land 
               N    = Hg0_Id_List(ID_Hg_eas)
               E_Hg = F_OF_PBL * EHg0_ln(I,J,ID_Hg_eas) * DTSRCE
               CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm )

               ! Anthro Southeast Asia re-emission from land 
               N    = Hg0_Id_List(ID_Hg_sea)
               E_Hg = F_OF_PBL * EHg0_ln(I,J,ID_Hg_sea) * DTSRCE
               CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm )

               ! Anthro Japan re-emission from land 
               N    = Hg0_Id_List(ID_Hg_jpn)
               E_Hg = F_OF_PBL * EHg0_ln(I,J,ID_Hg_jpn) * DTSRCE
               CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm )

               ! Anthro Oceania re-emission from land 
               N    = Hg0_Id_List(ID_Hg_oce)
               E_Hg = F_OF_PBL * EHg0_ln(I,J,ID_Hg_oce) * DTSRCE
               CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm )

               ! Soil re-emission from land 
               N    = Hg0_Id_List(ID_Hg_so)
               E_Hg = F_OF_PBL * EHg0_ln(I,J,ID_Hg_so) * DTSRCE
               CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm )

               ! Biomass burning re-emission from land 
               N    = Hg0_Id_List(ID_Hg_bb)
               E_Hg = F_OF_PBL * EHg0_ln(I,J,ID_Hg_bb) * DTSRCE
               CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm )

               ! Geogenic re-emission from land 
               N    = Hg0_Id_List(ID_Hg_geo)
               E_Hg = F_OF_PBL * EHg0_ln(I,J,ID_Hg_geo) * DTSRCE
               CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm )

               ! Below Mixed Layer Mid-Atlantic re-emission from land 
               N    = Hg0_Id_List(ID_Hg_atl)
               E_Hg = F_OF_PBL * EHg0_ln(I,J,ID_Hg_atl) * DTSRCE
               CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm )

               ! Below Mixed Layer North Atlantic re-emission from land 
               N    = Hg0_Id_List(ID_Hg_nat)
               E_Hg = F_OF_PBL * EHg0_ln(I,J,ID_Hg_nat) * DTSRCE
               CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm )

               ! Below Mixed Layer South Atlantic re-emission from land 
               N    = Hg0_Id_List(ID_Hg_sat)
               E_Hg = F_OF_PBL * EHg0_ln(I,J,ID_Hg_sat) * DTSRCE
               CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm )

               ! Below Mixed Layer North Pacific re-emission from land
               N    = Hg0_Id_List(ID_Hg_npa)
               E_Hg = F_OF_PBL * EHg0_ln(I,J,ID_Hg_npa) * DTSRCE
               CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm )

               ! Below Mixed Layer Arctic re-emission from land 
               N    = Hg0_Id_List(ID_Hg_arc)
               E_Hg = F_OF_PBL * EHg0_ln(I,J,ID_Hg_arc) * DTSRCE
               CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm )

               ! Below Mixed Layer Antarctic re-emission from land 
               N    = Hg0_Id_List(ID_Hg_ant)
               E_Hg = F_OF_PBL * EHg0_ln(I,J,ID_Hg_ant) * DTSRCE
               CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm )

               ! Below Mixed Layer All Other Ocean re-emission from land
               N    = Hg0_Id_List(ID_Hg_ocn)
               E_Hg = F_OF_PBL * EHg0_ln(I,J,ID_Hg_ocn) * DTSRCE
               CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm )

               ! Residual Stratosphere from Spin-up re-emission from land 
               N    = Hg0_Id_List(ID_Hg_str)
               E_Hg = F_OF_PBL * EHg0_ln(I,J,ID_Hg_str) * DTSRCE
               CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm )
               
               !--------------------
               ! Snow re-emissions
               !--------------------

               ! Anthro Canada re-emission from snow
               N    = Hg0_Id_List(ID_Hg_can)
               E_Hg = F_OF_PBL * EHg0_snow(I,J,ID_Hg_can) * DTSRCE
               CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm )

               ! Anthro United States re-emission from snow 
               N    = Hg0_Id_List(ID_Hg_usa)
               E_Hg = F_OF_PBL * EHg0_snow(I,J,ID_Hg_usa) * DTSRCE
               CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm )

               ! Anthro Central America re-emission from snow 
               N    = Hg0_Id_List(ID_Hg_cam)
               E_Hg = F_OF_PBL * EHg0_snow(I,J,ID_Hg_cam) * DTSRCE
               CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm )

               ! Anthro South America re-emission from snow 
               N    = Hg0_Id_List(ID_Hg_sam)
               E_Hg = F_OF_PBL * EHg0_snow(I,J,ID_Hg_sam) * DTSRCE
               CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm )

               ! Anthro West Africa re-emission from snow 
               N    = Hg0_Id_List(ID_Hg_waf)
               E_Hg = F_OF_PBL * EHg0_snow(I,J,ID_Hg_waf) * DTSRCE
               CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm )

               ! Anthro East Africa re-emission from snow 
               N    = Hg0_Id_List(ID_Hg_eaf)
               E_Hg = F_OF_PBL * EHg0_snow(I,J,ID_Hg_eaf) * DTSRCE
               CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm )

               ! Anthro South Africa re-emission from snow 
               N    = Hg0_Id_List(ID_Hg_saf)
               E_Hg = F_OF_PBL * EHg0_snow(I,J,ID_Hg_saf) * DTSRCE
               CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm )

               ! Anthro North Africa re-emission from snow 
               N    = Hg0_Id_List(ID_Hg_naf)
               E_Hg = F_OF_PBL * EHg0_snow(I,J,ID_Hg_naf) * DTSRCE
               CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm )

               ! Anthro OECD Europe re-emission from snow 
               N    = Hg0_Id_List(ID_Hg_eur)
               E_Hg = F_OF_PBL * EHg0_snow(I,J,ID_Hg_eur) * DTSRCE
               CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm )

               ! Anthro Eastern Europe re-emission from snow 
               N    = Hg0_Id_List(ID_Hg_eeu)
               E_Hg = F_OF_PBL * EHg0_snow(I,J,ID_Hg_eeu) * DTSRCE
               CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm )

               ! Anthro Middle East re-emission from snow 
               N    = Hg0_Id_List(ID_Hg_mde)
               E_Hg = F_OF_PBL * EHg0_snow(I,J,ID_Hg_mde) * DTSRCE
               CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm )

               ! Anthro Former Soviet Union re-emission from snow 
               N    = Hg0_Id_List(ID_Hg_sov)
               E_Hg = F_OF_PBL * EHg0_snow(I,J,ID_Hg_sov) * DTSRCE
               CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm )

               ! Anthro South Asia re-emission from snow 
               N    = Hg0_Id_List(ID_Hg_sas)
               E_Hg = F_OF_PBL * EHg0_snow(I,J,ID_Hg_sas) * DTSRCE
               CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm )

               ! Anthro East Asia re-emission from snow 
               N    = Hg0_Id_List(ID_Hg_eas)
               E_Hg = F_OF_PBL * EHg0_snow(I,J,ID_Hg_eas) * DTSRCE
               CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm )

               ! Anthro Southeast Asia re-emission from snow 
               N    = Hg0_Id_List(ID_Hg_sea)
               E_Hg = F_OF_PBL * EHg0_snow(I,J,ID_Hg_sea) * DTSRCE
               CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm )

               ! Anthro Japan re-emission from snow 
               N    = Hg0_Id_List(ID_Hg_jpn)
               E_Hg = F_OF_PBL * EHg0_snow(I,J,ID_Hg_jpn) * DTSRCE
               CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm )

               ! Anthro Oceania re-emission from snow 
               N    = Hg0_Id_List(ID_Hg_oce)
               E_Hg = F_OF_PBL * EHg0_snow(I,J,ID_Hg_oce) * DTSRCE
               CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm )

               ! Soil re-emission from snow
               N    = Hg0_Id_List(ID_Hg_so)
               E_Hg = F_OF_PBL * EHg0_snow(I,J,ID_Hg_so) * DTSRCE
               CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm )

               ! Biomass burning re-emission from snow 
               N    = Hg0_Id_List(ID_Hg_bb)
               E_Hg = F_OF_PBL * EHg0_snow(I,J,ID_Hg_bb) * DTSRCE
               CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm )

               ! Geogenic re-emission from snow 
               N    = Hg0_Id_List(ID_Hg_geo)
               E_Hg = F_OF_PBL * EHg0_snow(I,J,ID_Hg_geo) * DTSRCE
               CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm )

               ! Below Mixed Layer Mid-Atlantic re-emission from snow 
               N    = Hg0_Id_List(ID_Hg_atl)
               E_Hg = F_OF_PBL * EHg0_snow(I,J,ID_Hg_atl) * DTSRCE
               CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm )

               ! Below Mixed Layer North Atlantic re-emission from snow 
               N    = Hg0_Id_List(ID_Hg_nat)
               E_Hg = F_OF_PBL * EHg0_snow(I,J,ID_Hg_nat) * DTSRCE
               CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm )

               ! Below Mixed Layer South Atlantic re-emission from snow 
               N    = Hg0_Id_List(ID_Hg_sat)
               E_Hg = F_OF_PBL * EHg0_snow(I,J,ID_Hg_sat) * DTSRCE
               CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm )

               ! Below Mixed Layer North Pacific re-emission from snow
               N    = Hg0_Id_List(ID_Hg_npa)
               E_Hg = F_OF_PBL * EHg0_snow(I,J,ID_Hg_npa) * DTSRCE
               CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm )

               ! Below Mixed Layer Arctic re-emission from snow 
               N    = Hg0_Id_List(ID_Hg_arc)
               E_Hg = F_OF_PBL * EHg0_snow(I,J,ID_Hg_arc) * DTSRCE
               CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm )

               ! Below Mixed Layer Antarctic re-emission from snow 
               N    = Hg0_Id_List(ID_Hg_ant)
               E_Hg = F_OF_PBL * EHg0_snow(I,J,ID_Hg_ant) * DTSRCE
               CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm )

               ! Below Mixed Layer All Other Ocean re-emission from snow
               N    = Hg0_Id_List(ID_Hg_ocn)
               E_Hg = F_OF_PBL * EHg0_snow(I,J,ID_Hg_ocn) * DTSRCE
               CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm )

               ! Residual Stratosphere from Spin-up re-emission from snow 
               N    = Hg0_Id_List(ID_Hg_str)
               E_Hg = F_OF_PBL * EHg0_snow(I,J,ID_Hg_str) * DTSRCE
               CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm )

            ENDIF
         ENDDO
        
         !==============================================================
         ! ND03 diagnostic: Total Hg(0) emissions [kg]
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
                     AD03_nat(I,J,N) = DTSRCE * ( EHg0_geo(I,J)     +
     &                                            EHg0_bb(I,J)      +
     &                                            EHg0_gtm(I,J)     +
     &                                            EHg0_snow(I,J,N)  )
                  ENDDO
               ELSE
                  DO N = 1, N_HG_CATS
                     AD03_nat(I,J,N) = DTSRCE * ( EHg0_ln(I,J,N)   +
     &                                            EHg0_geo(I,J)    +
     &                                            EHg0_bb(I,J)     +
     &                                            EHg0_vg(I,J)     +
     &                                            EHg0_so(I,J)     +
     &                                            EHg0_snow(I,J,N) )
                  ENDDO
               ENDIF 
            ENDIF
            
         ENDIF
      ENDDO
      ENDDO
!$OMP END PARALLEL DO
  
      ! Return to calling program
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
      SUBROUTINE SRCHg2( am_I_Root, Input_Opt, State_Chm, RC )
!
! !USES:
!
      USE CMN_SIZE_MOD
      USE DIAG03_MOD,         ONLY : AD03, ND03
      USE ErrCode_Mod
      USE ERROR_MOD,          ONLY : ERROR_STOP
      USE Input_Opt_Mod,      ONLY : OptInput
      USE PBL_MIX_MOD,        ONLY : GET_FRAC_OF_PBL, GET_PBL_MAX_L
      USE State_Chm_Mod,      ONLY : ChmState
      USE State_Met_Mod,      ONLY : MetState
      USE TIME_MOD,           ONLY : GET_TS_EMIS 
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)    :: am_I_Root   ! Are we on the root CPU?
      TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
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
!  (1 ) Now use diagnostic arrays from "diag03_mod.f" (bmy, 1/21/05)
!  (2 ) Now references GET_FRAC_OF_PBL and GET_PBL_MAX_L from "pbl_mix_mod.f".
!        Remove reference to FEMIS. (bmy, 2/22/05)
!  (3 ) EHg2_an is now a 2-D array.  Now use ID_Hg2 index array from 
!        "tracerid_mod.f".  Now make sure STT does not underflow.
!        (eck, cdh, bmy, 4/6/06)
!  25 Mar 2013 - R. Yantosca - Now accept am_I_Root, Input_Opt, State_Chm, RC
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
      DTSRCE    = GET_TS_EMIS() * 60e+0_fp

      ! Maximum extent of the PBL [model levels]
      PBL_MAX   = GET_PBL_MAX_L()

      IF (.NOT. LPREINDHG ) THEN

      ! Loop over grid boxes
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, F_OF_PBL, E_Hg, N )
         DO J = 1, JJPAR
         DO I = 1, IIPAR

            ! Loop up to the max PBL layer
            DO L = 1, PBL_MAX

               ! Fraction of box (I,J,L) w/in the PBL [unitless]
               F_OF_PBL = GET_FRAC_OF_PBL( I, J, L )

               ! Partition total Hg2 into box (I,J,L) [kg]
               E_Hg     = F_OF_PBL * EHg2_an(I,J) * DTSRCE
                        
               !---------------------------
               ! Total anthro Hg(II) [kg]
               !---------------------------
               N    = Hg2_Id_List(ID_Hg_tot)
               CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm )

               !---------------------------
               ! Tagged anthro Hg(II) [kg]
               !---------------------------
               IF ( LSPLIT ) THEN 

                  ! Anthro Canada Hg2
                  N    = Hg2_Id_List(ID_Hg_can)
                  E_Hg = F_OF_PBL * EHg2_can(I,J) * DTSRCE
                  CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm )

                  ! Anthro United States Hg2
                  N    = Hg2_Id_List(ID_Hg_usa)
                  E_Hg = F_OF_PBL * EHg2_usa(I,J) * DTSRCE
                  CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm )

                  ! Anthro Central America Hg2
                  N    = Hg2_Id_List(ID_Hg_cam)
                  E_Hg = F_OF_PBL * EHg2_cam(I,J) * DTSRCE
                  CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm )

                  ! Anthro South America Hg2
                  N    = Hg2_Id_List(ID_Hg_sam)
                  E_Hg = F_OF_PBL * EHg2_sam(I,J) * DTSRCE
                  CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm )

                  ! Anthro West Africa Hg2
                  N    = Hg2_Id_List(ID_Hg_waf)
                  E_Hg = F_OF_PBL * EHg2_waf(I,J) * DTSRCE
                  CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm )
   
                  ! Anthro East Africa Hg2
                  N    = Hg2_Id_List(ID_Hg_eaf)
                  E_Hg = F_OF_PBL * EHg2_eaf(I,J) * DTSRCE
                  CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm )
   
                  ! Anthro South Africa Hg2
                  N    = Hg2_Id_List(ID_Hg_saf)
                  E_Hg = F_OF_PBL * EHg2_saf(I,J) * DTSRCE
                  CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm )
   
                  ! Anthro North Africa Hg2
                  N    = Hg2_Id_List(ID_Hg_naf)
                  E_Hg = F_OF_PBL * EHg2_naf(I,J) * DTSRCE
                  CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm )
   
                  ! Anthro OECD Europe Hg2
                  N    = Hg2_Id_List(ID_Hg_eur)
                  E_Hg = F_OF_PBL * EHg2_eur(I,J) * DTSRCE
                  CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm )
   
                  ! Anthro Eastern Europe Hg2
                  N    = Hg2_Id_List(ID_Hg_eeu)
                  E_Hg = F_OF_PBL * EHg2_eeu(I,J) * DTSRCE
                  CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm )
   
                  ! Anthro Middle East Hg2
                  N    = Hg2_Id_List(ID_Hg_mde)
                  E_Hg = F_OF_PBL * EHg2_mde(I,J) * DTSRCE
                  CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm )
   
                  ! Anthro Former USSR Hg2
                  N    = Hg2_Id_List(ID_Hg_sov)
                  E_Hg = F_OF_PBL * EHg2_sov(I,J) * DTSRCE
                  CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm )
   
                  ! Anthro South Asia Hg2
                  N    = Hg2_Id_List(ID_Hg_sas)
                  E_Hg = F_OF_PBL * EHg2_sas(I,J) * DTSRCE
                  CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm )
   
                  ! Anthro East Asia Hg2
                  N    = Hg2_Id_List(ID_Hg_eas)
                  E_Hg = F_OF_PBL * EHg2_eas(I,J) * DTSRCE
                  CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm )
   
                  ! Anthro Southeast Asia Hg2
                  N    = Hg2_Id_List(ID_Hg_sea)
                  E_Hg = F_OF_PBL * EHg2_sea(I,J) * DTSRCE
                  CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm )
   
                  ! Anthro Japan Hg2
                  N    = Hg2_Id_List(ID_Hg_jpn)
                  E_Hg = F_OF_PBL * EHg2_jpn(I,J) * DTSRCE
                  CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm )
   
                  ! Anthro Oceania Hg2
                  N    = Hg2_Id_List(ID_Hg_oce)
                  E_Hg = F_OF_PBL * EHg2_oce(I,J) * DTSRCE
                  CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm )

               ENDIF
            ENDDO
            
            !-------------------------------
            ! ND03 diag: Anthro Hg(II) [kg]
            !-------------------------------
            IF ( ND03 > 0 ) THEN
               AD03(I,J,6,1) = AD03(I,J,6,1) + ( EHg2_an(I,J) * DTSRCE )
            ENDIF
         ENDDO
         ENDDO
!$OMP END PARALLEL DO

      ENDIF

      ! Return to calling program
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
      SUBROUTINE SRCHgP( am_I_Root, Input_Opt, State_Chm, RC )
!
! !USES:
!
      USE CMN_SIZE_MOD
      USE DIAG03_MOD,         ONLY : AD03, ND03
      USE ErrCode_Mod
      USE ERROR_MOD,          ONLY : ERROR_STOP
      USE Input_Opt_Mod,      ONLY : OptInput
      USE PBL_MIX_MOD,        ONLY : GET_FRAC_OF_PBL, GET_PBL_MAX_L
      USE State_Chm_Mod,      ONLY : ChmState
      USE State_Met_Mod,      ONLY : MetState
      USE TIME_MOD,           ONLY : GET_TS_EMIS 
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)    :: am_I_Root   ! Are we on the root CPU?
      TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
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
!  (1 ) Now use diagnostic arrays from "diag03_mod.f" (bmy, 1/21/05)
!  (2 ) Now references GET_FRAC_OF_PBL and GET_PBL_MAX_L from "pbl_mix_mod.f".
!        Remove reference to FEMIS. (bmy, 2/22/05)
!  (3 ) EHgP_an is now a 2-D array.  Now use ID_HgP index array from 
!        "tracerid_mod.f".  Now make sure STT does not underflow. 
!        (eck, cdh, bmy, 4/6/06) 
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER :: I,      J,        L,    N,   PBL_MAX
      REAL(fp)  :: DTSRCE, F_OF_PBL, E_Hg 

      ! For values from Input_Opt
      LOGICAL :: LSPLIT, LPREINDHG

      !=================================================================
      ! SRCHgP begins here!
      !=================================================================

      ! Assume success
      RC        = GC_SUCCESS      

      ! Copy values from Input_Opt
      LSPLIT    = Input_Opt%LSPLIT
      LPREINDHG = Input_Opt%LPREINDHG

      ! Chemistry timestep [s]
      DTSRCE    = GET_TS_EMIS() * 60e+0_fp
      
      ! Maximum extent of the PBL [model levels]
      PBL_MAX   = GET_PBL_MAX_L()

      IF (.NOT. LPREINDHG) THEN

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, F_OF_PBL, E_Hg, N )
         DO J = 1, JJPAR
         DO I = 1, IIPAR

            ! Loop up to PBL top layer
            DO L = 1, PBL_MAX

               ! Fraction of box (I,J,L) w/in the PBL [unitless]
               F_OF_PBL           = GET_FRAC_OF_PBL( I, J, L )

               ! Partition HgP into box (I,J,L) [kg]
               E_Hg               = F_OF_PBL * EHgP_an(I,J) * DTSRCE
            
               !------------------------
               ! Total anthro HgP [kg]
               !------------------------
               N               = HgP_Id_List(ID_Hg_tot)
               CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm )

               !------------------------
               ! Tagged anthro HgP [kg]
               !------------------------
               IF ( LSPLIT ) THEN

               ! Anthro Canada HgP
               N            = HgP_Id_List(ID_Hg_can)
               E_Hg         = F_OF_PBL * EHgP_can(I,J) * DTSRCE
               CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm )

               ! Anthro United States HgP
               N            = HgP_Id_List(ID_Hg_usa)
               E_Hg         = F_OF_PBL * EHgP_usa(I,J) * DTSRCE
               CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm )

               ! Anthro Central America HgP
               N            = HgP_Id_List(ID_Hg_cam)
               E_Hg         = F_OF_PBL * EHgP_cam(I,J) * DTSRCE
               CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm )

               ! Anthro South America HgP
               N            = HgP_Id_List(ID_Hg_sam)
               E_Hg         = F_OF_PBL * EHgP_sam(I,J) * DTSRCE
               CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm )

               ! Anthro West Africa HgP
               N            = HgP_Id_List(ID_Hg_waf)
               E_Hg         = F_OF_PBL * EHgP_waf(I,J) * DTSRCE
               CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm )

               ! Anthro East Africa HgP
               N            = HgP_Id_List(ID_Hg_eaf)
               E_Hg         = F_OF_PBL * EHgP_eaf(I,J) * DTSRCE
               CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm )

               ! Anthro South Africa HgP
               N            = HgP_Id_List(ID_Hg_saf)
               E_Hg         = F_OF_PBL * EHgP_saf(I,J) * DTSRCE
               CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm )

               ! Anthro North Africa HgP
               N            = HgP_Id_List(ID_Hg_naf)
               E_Hg         = F_OF_PBL * EHgP_naf(I,J) * DTSRCE
               CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm )

               ! Anthro OECD Europe HgP
               N            = HgP_Id_List(ID_Hg_eur)
               E_Hg         = F_OF_PBL * EHgP_eur(I,J) * DTSRCE
               CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm )

               ! Anthro Eastern Europe HgP
               N            = HgP_Id_List(ID_Hg_eeu)
               E_Hg         = F_OF_PBL * EHgP_eeu(I,J) * DTSRCE
               CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm )

               ! Anthro Middle East HgP
               N            = HgP_Id_List(ID_Hg_mde)
               E_Hg         = F_OF_PBL * EHgP_mde(I,J) * DTSRCE
               CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm )

               ! Anthro Former USSR HgP
               N            = HgP_Id_List(ID_Hg_sov)
               E_Hg         = F_OF_PBL * EHgP_sov(I,J) * DTSRCE
               CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm )

               ! Anthro South Asia HgP
               N            = HgP_Id_List(ID_Hg_sas)
               E_Hg         = F_OF_PBL * EHgP_sas(I,J) * DTSRCE
               CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm )

               ! Anthro East Asia HgP
               N            = HgP_Id_List(ID_Hg_eas)
               E_Hg         = F_OF_PBL * EHgP_eas(I,J) * DTSRCE
               CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm )

               ! Anthro Southeast Asia HgP
               N            = HgP_Id_List(ID_Hg_sea)
               E_Hg         = F_OF_PBL * EHgP_sea(I,J) * DTSRCE
               CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm )

               ! Anthro Japan HgP
               N            = HgP_Id_List(ID_Hg_jpn)
               E_Hg         = F_OF_PBL * EHgP_jpn(I,J) * DTSRCE
               CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm )

               ! Anthro Oceania HgP
               N            = HgP_Id_List(ID_Hg_oce)
               E_Hg         = F_OF_PBL * EHgP_oce(I,J) * DTSRCE
               CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm )

               ENDIF
            ENDDO

            !----------------------------
            ! ND03 diag: Anthro HgP [kg]
            !----------------------------
            IF ( ND03 > 0 ) THEN
               AD03(I,J,9,1) = AD03(I,J,9,1) + (EHgP_an(I,J) * DTSRCE)
            ENDIF
         
         ENDDO 
         ENDDO
!$OMP END PARALLEL DO

      ENDIF
         
      ! Return to calling program
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
      SUBROUTINE MERCURY_READYR( am_I_Root, Input_Opt, RC )
!
! !USES:
!
      USE ErrCode_Mod
      USE ERROR_MOD,          ONLY : ERROR_STOP
      USE HCO_ERROR_MOD
      USE HCO_INTERFACE_MOD,  ONLY : HcoState
      USE HCO_EMISLIST_MOD,   ONLY : HCO_GetPtr
      USE HCO_INTERFACE_MOD,  ONLY : GetHcoDiagn 
      USE Input_Opt_Mod,      ONLY : OptInput
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)    :: am_I_Root   ! Are we on the root CPU?
      TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
! 
! !REVISION HISTORY: 
!  06 Dec 2004 - N. (Eckley) Selin - Initial version
!  (1 ) Now read data from mercury_200501 subdirectory.  Now compute oceanic 
!        Hg(0) emissions w/ ocean flux module instead of reading them from 
!        disk.  Now use 1985 TAU values. (sas, bmy, 1/20/05) 
!  (2 ) Now can read data for both GEOS and GCAP grids (bmy, 8/16/05)
!  (3 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (4 ) Now read anthro emissions on GEOS 1x1 grid in DATA_1x1_DIR.  Also
!        keep 2x25 and 4x5 files together in DATA_1x1_DIR.  Also now use new
!        land re-emissions files from Noelle Selin. (eck, bmy, 4/6/06)
!  (5 ) Now reads anthro emissions for 17 world regions from David Streets
!       for use in tagged simulations. (eds 2/7/12)
!  (6 ) Bug fix: emissions are read on GENERIC grid (jaf, 2/8/12)
!  01 Mar 2012 - R. Yantosca - Now use GET_XMID(I,J,L) from grid_mod.F90
!  01 Mar 2012 - R. Yantosca - Now use GET_YMID(I,J,L) from grid_mod.F90
!  13 Mar 2012 - M. Cooper   - Changed regrid algorithm to map_a2a
!  24 May 2012 - R. Yantosca - Fixed minor bugs in map_a2a algorithm
!  24 May 2012 - H. Amos     - updated file paths for soil Hg files
!  28 Jun 2012 - S. Kim      - Bug fix: you need to provide a location
!                              in the call to ERROR_STOP
!  24 Aug 2012 - R. Yantosca - DO_REGRID_A2A now reads netCDF input file
!  03 Jan 2013 - M. Payer    - Renamed PERAREA to IS_MASS in DO_REGRID_A2A
!  13 Jan 2013 - M. Payer    - Update soil Hg files for v9-02c (Y. Zhang)
!  06 Jun 2013 - R. Yantosca - Get ANTHRO_Hg_YEAR, Hg_SCENARIO from Input_Opt
!  06 Sep 2013 - R. Yantosca - Kludge: use MERRA soil file for GEOS-FP testing
!                              and debugging.  Worry about this later.
!  06 Sep 2013 - R. Yantosca - Kludge: use GEOS-5 soil file for GEOS-4 testing
!                              and debugging.  Worry about this later.
!  24 Sep 2014 - C. Keller   - HEMCO update
!  29 Apr 2016 - R. Yantosca - Don't initialize pointers in declaration stmts
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
      INTEGER                  :: ANTHRO_Hg_YEAR
      CHARACTER(LEN=10)        :: SCENARIO

      ! HEMCO update
      CHARACTER(LEN=63)        :: DgnName
      REAL(f4),        POINTER :: Ptr2D(:,:)

      ! Strings
      CHARACTER(LEN=255)       :: MSG
      CHARACTER(LEN=255)       :: LOC = 'MERCURY_READYR (mercury_mod.F)'

      !=================================================================
      ! MERCURY_READYR begins here! 
      !=================================================================

      ! Assume success
      RC             = GC_SUCCESS

      ! Initialize pointers
      Ptr2D          => NULL()  

      ! Copy values from Input_Opt
      LDYNOCEAN      = Input_Opt%LDYNOCEAN 
      LPREINDHG      = Input_Opt%LPREINDHG
      LSPLIT         = Input_Opt%LSPLIT
      ANTHRO_HG_YEAR = Input_Opt%ANTHRO_Hg_YEAR
      SCENARIO       = TRIM( Input_Opt%Hg_SCENARIO )

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
         CALL ERROR_STOP ( 'HcoState not associated!', LOC )
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
         CALL GetHcoDiagn( am_I_Root, DgnName, .TRUE., RC, Ptr2D=Ptr2D )
         IF ( RC /= HCO_SUCCESS ) CALL ERROR_STOP ( DgnName, LOC )
         EHg0_an =  Ptr2D(:,:) * HcoState%Grid%AREA_M2%Val(:,:)
         Ptr2D   => NULL()

         ! Artisanal emissions
         DgnName = 'HG0_ARTISANAL'
         CALL GetHcoDiagn( am_I_Root, DgnName, .TRUE., RC, Ptr2D=Ptr2D )
         IF ( RC /= HCO_SUCCESS ) CALL ERROR_STOP ( DgnName, LOC )
         EHg0_am =  Ptr2D(:,:) * HcoState%Grid%AREA_M2%Val(:,:)
         Ptr2D   => NULL()

         ! ---------------------------------------------
         ! Hg2 emissions
         ! Note: HgP emissions are added to Hg2 by HEMCO
         ! ---------------------------------------------

         ! Anthropogenic emissions
         DgnName = 'HG2_ANTHRO'
         CALL GetHcoDiagn( am_I_Root, DgnName, .TRUE., RC, Ptr2D=Ptr2D )
         IF ( RC /= HCO_SUCCESS ) CALL ERROR_STOP ( DgnName, LOC )
         EHg2_an =  Ptr2D(:,:) * HcoState%Grid%AREA_M2%Val(:,:)
         Ptr2D   => NULL()
      ENDIF
 
      ! Natural emissions
      DgnName = 'HG0_NATURAL'
      CALL GetHcoDiagn( am_I_Root, DgnName, .TRUE., RC, Ptr2D=Ptr2D )
      IF ( RC /= HCO_SUCCESS ) CALL ERROR_STOP ( DgnName, LOC )
      EHg0_geo =  Ptr2D(:,:) * HcoState%Grid%AREA_M2%Val(:,:)
      Ptr2D    => NULL()

      ! Soil distribution 
      CALL HCO_GetPtr( am_I_Root, HcoState, 'HG0_SOILDIST', Ptr2D, RC )
      IF ( RC /= HCO_SUCCESS ) THEN
         CALL ERROR_STOP ( 'Cannot get pointer to HG0_SOILDIST', LOC )
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
 210  FORMAT( 'M E R C U R Y   E M I S S I O N S' )
 211  FORMAT( 'Total Anthro     Hg(0)  : ', f7.3, ' [Gg/yr]' )
 213  FORMAT( 'Total Re-Emitted Hg(0)  : ', f7.3, ' [Gg/yr]' )
 214  FORMAT( 'Total Natural    Hg(0)  : ', f7.3, ' [Gg/yr]' )
 215  FORMAT( 'Total Anthro     Hg(II) : ', f7.3, ' [Gg/yr]' )
 216  FORMAT( 'Total Anthro     HgP    : ', f7.3, ' [Gg/yr]' )
 217  FORMAT( 'Total Ocean      Hg(0)  : ', f7.3, ' [Gg/yr]' )

      ! Return to calling program
      END SUBROUTINE MERCURY_READYR
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_lwc
!
! !DESCRIPTION: Function GET\_LWC returns the cloud liquid water content 
!  [m3 H2O/m3 air] at a GEOS-CHEM grid box as a function of temperature.
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_LWC( T ) RESULT( LWC )
!
! !USES:
!

!
! !INPUT PARAMETERS: 
!
      REAL(fp), INTENT(IN) :: T     ! Temperature [K]
!
! !RETURN VALUE:
!
      REAL(fp)             :: LWC   ! Liquid water content [m3 H2O/m3 air]
! 
! !REVISION HISTORY: 
!  31 Oct 2002 - R. Park     - Initial version
!  18 Jan 2011 - R. Yantosca - Updated comments
!EOP
!------------------------------------------------------------------------------
!BOC
      !=================================================================
      ! GET_LWC begins here!
      !=================================================================

      ! Compute Liquid water content in [g/m3]
      IF (T>293e+0_fp) THEN
         LWC=0.2e+0_fp
       ELSE IF ( T >= 280.e+0_fp .AND. T <= 293.e+0_fp ) THEN
         LWC = 0.32e+0_fp - 0.0060e+0_fp * ( T - 273.e+0_fp ) 
 
      ELSE IF ( T >= 248.e+0_fp .AND. T < 280.e+0_fp ) THEN
         LWC = 0.23e+0_fp + 0.0065e+0_fp * ( T - 273.e+0_fp )

      ELSE IF ( T < 248.e+0_fp ) THEN
         LWC = 0.07e+0_fp

      ENDIF

      ! Convert from [g/m3] to [m3/m3]
      ! Units: [g H2O/m3 air] * [1 kg H2O/1000g H2O] * [m3 H2O/1000kg H2O]
      LWC = LWC * 1.e-6_fp         

      ! Return to calling program
      END FUNCTION GET_LWC
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_vcldf 
!
! !DESCRIPTION: Subroutine GET\_VCLDF computes the volume cloud fraction for 
!  Hg0 and Hg2 chemistry.
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_VCLDF( I, J, L, State_Met ) RESULT( VCLDF )
!
! !USES:
!
      USE State_Met_Mod,      ONLY : MetState
!
! !INPUT PARAMETERS: 
!
      INTEGER,        INTENT(IN)  :: I           ! Longitude index
      INTEGER,        INTENT(IN)  :: J           ! Latitude index
      INTEGER,        INTENT(IN)  :: L           ! Level index
      TYPE(MetState), INTENT(IN)  :: State_Met   ! Meteorology State object
!
! !RETURN VALUE:
!
      REAL(fp)                      :: VCLDF       ! Volume cloud fraction
!
! !REMARKS:
!  References:
!  ============================================================================
!  (1  ) Sundqvist et al. [1989]
! 
! !REVISION HISTORY: 
!  06 Dec 2004 - N. (Eckley) Selin - Copied from "sulfate_mod.f" but was made 
!                                    into a function since we are already 
!                                    looping over (I,J,L) in CHEM_Hg0_Hg2
!  09 Nov 2012 - M. Payer          - Replaced all met field arrays with 
!                                    State_Met derived type object
!  26 Feb 2015 - E. Lundgren       - Replace GET_PCENTER and GET_PEDGE with 
!                                    State_Met%PMID and State_Met%PEDGE. 
!                                    Remove dependency on pressure_mod.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL(fp)             :: PRES, PSFC, RH2, R0, B0
!
! !DEFINED PARAMETERS:
!
      REAL(fp),  PARAMETER :: ZRT = 0.60e+0_fp, ZRS = 0.99e+0_fp
	
      !=================================================================
      ! GET_VCLDF begins here!
      !=================================================================
     
      ! Surface pressure
      PSFC = State_Met%PEDGE(I,J,1)         

      ! Pressure at center of grid box L
      PRES = State_Met%PMID(I,J,L)
        
      ! RH (from "dao_mod.f") is relative humidity [%]
      ! Convert to fraction and store in RH2
      RH2  = State_Met%RH(I,J,L) * 1.0e-2_fp
        
      ! Terms from Sundqvist ???
      R0   = ZRT + ( ZRS - ZRT ) * EXP( 1e+0_fp - ( PSFC / PRES )**2.5 )
      B0   = ( RH2 - R0 ) / ( 1e+0_fp - R0 )
	
      ! Force B0 into the range 0-1
      IF ( RH2 < R0  ) B0 = 0e+0_fp
      IF ( B0  > 1e+0_fp ) B0 = 1e+0_fp

      ! Volume cloud fraction
      VCLDF = 1e+0_fp - SQRT( 1e+0_fp - B0 )

      ! Return to calling program
      END FUNCTION GET_VCLDF
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
      USE CMN_SIZE_MOD
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
!  (1 ) We assume SETTRACE has been called to define IDO3. (bmy, 12/16/02)
!  (2 ) Now reference inquiry functions from "tracer_mod.f" (bmy, 7/20/04)
!  09 Nov 2012 - M. Payer    - Replaced all met field arrays with State_Met
!                              derived type object
!  24 Jul 2014 - R. Yantosca - Now compute BOXVL internally
!  11 Mar 2015 - R. Yantosca - Remove reference to GLOBAL_O3_MOD
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
      O3_MOLEC_CM3 = O3(I,J,L)       * ( AVO * 1e+3_fp / AIRMW ) * 
     &               State_Met%AD(I,J,L) / BOXVL

      ! Return to calling program
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
      USE CMN_SIZE_MOD
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
!  28 Nov 2012 - R. Yantosca - Replace SUNCOS with State_Met%SUNCOS
!  28 Nov 2012 - R. Yantosca - Now pass State_Met via the argument list
!  11 Mar 2015 - R. Yantosca - Remove reference to GLOBAL_OH_MOD

!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      !=================================================================
      ! GET_OH begins here!
      !=================================================================

      ! Test for sunlight...
      IF ( State_Met%SUNCOS(I,J) > 0e+0_fp .and. 
     &     TCOSZ(I,J) > 0e+0_fp ) THEN

         ! Impose a diurnal variation on OH during the day
         OH_MOLEC_CM3 = OH(I,J,L)                                 *    
     &                  ( State_Met%SUNCOS(I,J) / TCOSZ(I,J)    ) *
     &                  ( 1440e+0_fp            / GET_TS_CHEM() )

         ! Make sure OH is not negative
         OH_MOLEC_CM3 = MAX( OH_MOLEC_CM3, 0e+0_fp )
               
      ELSE

         ! At night, OH goes to zero
         OH_MOLEC_CM3 = 0e+0_fp

      ENDIF

      ! Return to calling program
      END FUNCTION GET_OH
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
      FUNCTION GET_BR( I, J, L, BRO_MOLEC_CM3, State_Met )
     &   RESULT( BR_MOLEC_CM3 )
!
! !USES:
!
      USE CHEMGRID_MOD,       ONLY : ITS_IN_THE_STRAT
      USE CMN_SIZE_MOD
      USE DAO_MOD,            ONLY : GET_OBK
      USE DAO_MOD,            ONLY : IS_WATER
      USE DIAG03_MOD,         ONLY : AD03_Br,   LD03,     ND03
      USE ERROR_MOD,          ONLY : SAFE_DIV
      USE GLOBAL_BR_MOD,      ONLY : BR_MERGE,  BR_TROP
      USE GLOBAL_BR_MOD,      ONLY : BR_STRAT,  J_BRO
      USE GLOBAL_BR_MOD,      ONLY : BRO_MERGE, BRO_TROP, BRO_STRAT
      USE GC_GRID_MOD,        ONLY : GET_YMID
      USE PBL_MIX_MOD,        ONLY : GET_FRAC_UNDER_PBLTOP
      USE State_Met_Mod,      ONLY : MetState
      USE TIME_MOD,           ONLY : GET_TS_CHEM, GET_LOCALTIME 
      USE TIME_MOD,           ONLY : ITS_A_NEW_DAY
      USE TIME_MOD,           ONLY : GET_MONTH
!
! !INPUT PARAMETERS: 
!
      INTEGER,        INTENT(IN)  :: I              ! Longitude index
      INTEGER,        INTENT(IN)  :: J              ! Latitude index
      INTEGER,        INTENT(IN)  :: L              ! Level index
      REAL(fp),       INTENT(OUT) :: BrO_MOLEC_CM3  ! BrO conc [molec/cm3]
      TYPE(MetState), INTENT(IN)  :: State_Met      ! Meteorology State object
!
! !RETURN VALUE:
!
      REAL(fp)                    :: BR_MOLEC_CM3   ! Br conc [molec/cm3]
! 
! !REVISION HISTORY: 
!  06 Jul 2006 - C. Holmes   - Initial version
!  01 Mar 2012 - R. Yantosca - Now use GET_YMID(I,J,L) from grid_mod.F90
!  01 Mar 2012 - R. Yantosca - Now use GET_LOCALTIME(I,J,L) from time_mod.F90
!  09 Nov 2012 - M. Payer    - Replaced all met field arrays with State_Met
!                              derived type object
!  28 Nov 2012 - R. Yantosca - Replace SUNCOS with State_Met%SUNCOS
!  24 Jul 2014 - R. Yantosca - Now compute BOXVL internally
!  26 Feb 2015 - E. Lundgren - Remove dependency on pressure_mod (not used)
!  10 Aug 2015 - J. Fisher - implement variable polar BrO concentrations based
!                            on temperature and incoming solar radiation
!  10 Aug 2015 - J. Fisher - only apply diurnal cycle to BrO from pTOMCAT
!                            daily means, NOT to polar BrO which represents
!                            instantaneous values
!  10 Aug 2015 - J. Fisher - add upper ice threshold to MERRA conditions for
!                            polar BrO release
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
      IF ( (ITS_IN_THE_STRAT( I, J, L, State_Met )) .AND.
     &     (ABS(STRAT_BR_FACTOR - 1e+0_fp) > 1e-2_fp) )    THEN

         BR_PPTV  = BR_PPTV  * STRAT_BR_FACTOR
         BRO_PPTV = BRO_PPTV * STRAT_BR_FACTOR
         
      ENDIF 

      !----------------------------------------------------------------
      ! Add Br in the MBL
      !----------------------------------------------------------------

      ! Fraction of grid box mass which is within the mixed layer
      FPBL = GET_FRAC_UNDER_PBLTOP(I,J,L)

      IF ( IS_WATER( I, J, State_Met ) .AND. (FPBL > 0e+0_fp) ) THEN

         ! Convert BrO concentration to corresponding Br based on
         ! photochemical steady state, pptv
         ! Platt & Janssen, Faraday Discussions (1995)
         BR_MBL = BRO_MBL *( J_BRO(I,J,L) + K_BRO_NO * C_NO ) / 
     &                     ( K_BR_O3      * GET_O3(I,J,L,State_Met) ) *
     &                       FPBL

         ! MBL concentration is the greater of the TOMCAT Br or 1pptv BrO
         BR_PPTV  = MAX( BR_PPTV,  BR_MBL )
         BrO_PPTV = MAX( BRO_PPTV, BRO_MBL*FPBL )


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
!      IF ( ( LPOLARBR ) .AND. ( FPBL >  0e+0_fp   ) .AND.
!     &     ( State_Met%SUNCOS(I,J)   >  0e+0_fp   ) .AND. 
!     &     ( State_Met%T(I,J,1)      <= 268e+0_fp ) .AND.
!     &     ( State_Met%LWI(I,J)      == 2     ) .AND. 
!     &       GET_OBK(I,J,State_Met)  >  0e+0_fp   ) THEN
!
!         IF ( ( (GET_YMID(I,J,L) > 60e+0_fp) .AND. 
!     &          (GET_MONTH() >= 3) .AND. (GET_MONTH() <= 5)  ) .OR.
!     &        ( (GET_YMID(I,J,L) < -50e+0_fp) .AND. 
!     &          (GET_MONTH() >= 9) .AND. (GET_MONTH() <= 11)  ) ) THEN
!
!            ! Br concentration due to bromine explosion, pptv
!            ! Assume [O3] is 2ppb during event
!            BR_MBL = BRO_POLAR *( J_BRO(I,J,L) + K_BRO_NO * C_NO ) / 
!     &                          ( K_BR_O3 * 5e+10_fp ) * FPBL
!        
!            BR_PPTV  = BR_PPTV  + BR_MBL
!            BrO_PPTV = BrO_PPTV + BRO_POLAR * FPBL
!
!         ENDIF
!
!      ENDIF

      !----------------------------------------------------------------
      ! Impose a diurnal cycle
      !
      ! This formula comes from the Holmes et al (2009) MBL box model
      ! which fit RGM observations for Okinawa and a Pacific cruise
      !----------------------------------------------------------------
            
      ! Test for sunlight...
      IF ( (State_Met%SUNCOS(I,J) > 0e+0_fp) .and. 
     &     (TTDAY(I,J) > 0e+0_fp) ) THEN
      
         ! Use a constant function if daylight is < 2 hours or > 22hours
         IF ( ( TTDAY(I,J) < 120e+0_fp  ) .OR. 
     &        ( TTDAY(I,J) > 1320e+0_fp ) ) THEN
            
            BR_FAC = SAFE_DIV( 1440e+0_fp, TTDAY(I,J), 0e+0_fp )  
            
         ELSE

            ! Local time: 0-23.999
            LOCALTIME = GET_LOCALTIME( I, J, L )
      
            ! Interpolate the real hour to lie within an equinoctal day
            ! i.e. between 6-18
            HOUR = 12e+0_fp + ( LOCALTIME - 12e+0_fp ) *
     &           720e+0_fp / ( TTDAY(I,J) ) 

            ! Multiplicative factor for the diurnal cycle 
            BR_FAC = (1e+6_fp - 1e+2_fp * (HOUR - 6e+0_fp) ** 4e+0_fp + 
     &           ( 18e+0_fp - HOUR ) * 1e+6_fp / 12e+0_fp ) / 2e+6_fp

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
      BOXVL         = State_Met%AIRVOL(I,J,L) * 1e+6_fp

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
! KLUDGE: can't get SWGDN in State_Met; add this later (jaf)
! Use incident radiation at ground (State_Met%SWGDN) where available. 
! Otherwise, use net radiation at ground (State_Met%RADSWG)
#if   defined( MERRA ) || defined( GEOS_FP ) || defined( MERRA2 )
      SWRAD         = State_Met%SWGDN(I,J)
!      SWRAD = State_Met%RADSWG(I,J) - old code (jaf)
! END KLUDGE
      IS_MOSTLY_ICE = ( State_Met%SEAICE00(I,J) <= 0.5e+0_fp .AND.
     &                  State_Met%SEAICE90(I,J) < 1e+0_fp )
#else
      SWRAD = State_Met%RADSWG(I,J)
      IS_MOSTLY_ICE = ( State_Met%LWI(I,J) == 2 )
#endif

      IF ( (LPOLARBR) .AND. (FPBL > 0e+0_fp) .AND. (IS_MOSTLY_ICE) .AND.
     &     (SWRAD > 1e+2_fp) .AND. (State_Met%TS(I,J) <= 273e+0_fp) ) 
     &     THEN 

         IF ( ((GET_YMID(I,J,L) > 0e+0_fp) .AND.
     &          (GET_MONTH() >= 2) .AND. (GET_MONTH() <= 6)) .OR.
     &        ((GET_YMID(I,J,L) < 0e+0_fp) .AND.
     &          (GET_MONTH() >= 8) .AND. (GET_MONTH() <= 12)) ) THEN

            ! [BrO] is a linear function of temperature derived based on
            ! results from Pohler et al. (2010), Prados-Roman et al. (2011)
            ! and ability to match Hg0 seasonal cycle at Alert. (jaf,
            ! 12/24/11)
            IF ( State_Met%TS(I,J) <= 253e+0_fp ) THEN 
               BRO_POLAR = 2e+1_fp 
            ELSE 
               BRO_POLAR = -1e+0_fp * ( State_Met%TS(I,J) - 
     &                                  253e+0_fp ) + 2e+1_fp
            ENDIF

            ! Convert O3 to molec/cm3
            O3_POLAR = O3_POLAR_PPB * 1e-9_fp * ( AVO * 1e+3_fp / 
     &                 AIRMW ) * State_Met%AD(I,J,L) / 
     &                 BOXVL

            ! Compute polar Br, BrO concentrations in pptv
            BrO_POLAR = BrO_POLAR * FPBL
            Br_POLAR  = BRO_POLAR * ( J_BRO(I,J,L) + K_BRO_NO * C_NO ) /
     &                              ( K_BR_O3 * O3_POLAR )
         
            BR_PPTV  = BR_PPTV  + Br_POLAR
            BrO_PPTV = BrO_PPTV + BrO_POLAR

            !==============================================================
            ! ND03 diagnostic: Polar Bromine [pptv]
            !==============================================================
            IF ( ND03 > 0 .AND. L <= LD03 ) THEN 

                  AD03_Br(I,J,L,3) = AD03_Br(I,J,L,3) + BR_POLAR
                  AD03_Br(I,J,L,4) = AD03_Br(I,J,L,4) + BRO_POLAR
                  AD03_Br(I,J,L,5) = AD03_Br(I,J,L,5) + O3_POLAR_PPB

            ENDIF ! ND03

         ENDIF ! Month

      ENDIF ! Polar BrO criteria
      !----------------------------------------------------------------

      ! Convert pptv mixing ratio -> molec/cm3
      BR_MOLEC_CM3  = BR_PPTV             * 1e-12_fp      *
     &                ( AVO * 1e+3_fp     / AIRMW ) * 
     &                State_Met%AD(I,J,L) / BOXVL
      BRO_MOLEC_CM3 = BRO_PPTV            * 1e-12_fp      *
     &                ( AVO * 1e+3_fp     / AIRMW ) * 
     &                State_Met%AD(I,J,L) / BOXVL

      ! Return to calling program
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
      SUBROUTINE OHNO3TIME
!
! !USES:
!
      USE CMN_SIZE_MOD
      USE GC_GRID_MOD, ONLY : GET_XMID,    GET_YMID_R
      USE TIME_MOD,    ONLY : GET_NHMSb,   GET_ELAPSED_SEC
      USE TIME_MOD,    ONLY : GET_TS_CHEM, GET_DAY_OF_YEAR, GET_GMT
! 
! !REVISION HISTORY: 
!  16 Dec 2002 - R. Park & R. Yantosca - Initial version
!  (1 ) Copy code from COSSZA directly for now, so that we don't get NaN
!        values.  Figure this out later (rjp, bmy, 1/10/03)
!  (2 ) Now replace XMID(I) with routine GET_XMID from "grid_mod.f".  
!        Now replace RLAT(J) with routine GET_YMID_R from "grid_mod.f". 
!        Removed NTIME, NHMSb from the arg list.  Now use GET_NHMSb,
!        GET_ELAPSED_SEC, GET_TS_CHEM, GET_DAY_OF_YEAR, GET_GMT from 
!        "time_mod.f". (bmy, 3/27/03)
!  (3 ) Now store the peak SUNCOS value for each surface grid box (I,J) in 
!        the COSZM array. (rjp, bmy, 3/30/04)
!  (4 ) Also added parallel loop over grid boxes (eck, bmy, 12/8/04)
!  01 Mar 2012 - R. Yantosca - Now use GET_XMID(I,J,L) from grid_mod.F90
!  01 Mar 2012 - R. Yantosca - Now use GET_YMID_R(I,J,L) from grid_mod.F90
!  16 May 2016 - M. Sulprizio- Remove IJLOOP and change SUNTMP array dimensions
!                              from (MAXIJ) to (IIPAR,JJPAR)
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! References to F90 modules
      ! Local variables
      LOGICAL, SAVE :: FIRST = .TRUE.
      INTEGER       :: I, J, L, N, NT, NDYSTEP
      REAL(fp)      :: A0, A1, A2, A3, B1, B2, B3
      REAL(fp)      :: LHR0, R, AHR, DEC, TIMLOC, YMID_R
      REAL(fp)      :: SUNTMP(IIPAR,JJPAR)
      
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

      DEC = A0 - A1*cos(  R) + B1*sin(  R)
     &         - A2*cos(2*R) + B2*sin(2*R)
     &         - A3*cos(3*R) + B3*sin(3*R)

      LHR0 = int(float( GET_NHMSb() )/10000.)

      ! Only do the following at the start of a new day
      IF ( FIRST .or. GET_GMT() < 1e-5 ) THEN 
      
         ! Zero arrays
         TTDAY(:,:) = 0e+0_fp
         TCOSZ(:,:) = 0e+0_fp
         COSZM(:,:) = 0e+0_fp

         ! NDYSTEP is # of chemistry time steps in this day
         NDYSTEP = ( 24 - INT( GET_GMT() ) ) * 60 / GET_TS_CHEM()     

         ! NT is the elapsed time [s] since the beginning of the run
         NT = GET_ELAPSED_SEC()

         ! Loop forward through NDYSTEP "fake" timesteps for this day 
         DO N = 1, NDYSTEP
            
            ! Zero SUNTMP array
            SUNTMP = 0e+0_fp

            ! Loop over surface grid boxes
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, YMID_R, TIMLOC, AHR )
            DO J = 1, JJPAR
            DO I = 1, IIPAR

               ! Grid box latitude center [radians]
               YMID_R = GET_YMID_R( I, J, 1 )

               TIMLOC = real(LHR0) + real(NT)/3600.0 
     &                + GET_XMID(I,J,1)/15.0
         
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
               SUNTMP(I,J) = sin(YMID_R) * sin(DEC) +
     &                       cos(YMID_R) * cos(DEC) * cos(AHR)

               ! TCOSZ is the sum of SUNTMP at location (I,J)
               ! Do not include negative values of SUNTMP
               TCOSZ(I,J) = TCOSZ(I,J) + MAX( SUNTMP(I,J), 0e+0_fp )

               ! COSZM is the peak value of SUMTMP during a day at (I,J)
               ! (rjp, bmy, 3/30/04)
               COSZM(I,J) = MAX( COSZM(I,J), SUNTMP(I,J) )

               ! TTDAY is the total daylight time at location (I,J)
               IF ( SUNTMP(I,J) > 0e+0_fp ) THEN
                  TTDAY(I,J) = TTDAY(I,J) + DBLE( GET_TS_CHEM() )
               ENDIF
            ENDDO
            ENDDO
!$OMP END PARALLEL DO

            ! Increment elapsed time [sec]
            NT = NT + ( GET_TS_CHEM() * 60 )             
         ENDDO

         ! Reset first-time flag
         FIRST = .FALSE.
      ENDIF

      ! Return to calling program
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
      SUBROUTINE CALC_HG2_SEASALT_LOSSRATE( State_Met )
!
! !USES:
!
      USE CMN_SIZE_MOD             ! Size parameters
      USE DAO_MOD,            ONLY : IS_WATER
      USE PBL_MIX_MOD,        ONLY : GET_PBL_TOP_M
      USE State_Met_Mod,      ONLY : MetState
!
! !INPUT PARAMETERS: 
!
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
!  25 Jul 2014 - R. Yantosca - Remove reference to function SFCWINDSQR
!  04 Dec 2014 - M. Yannetti - Remove reference to RNPBE_MOD
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      ! Local variables
      REAL(fp)                :: U10M, S
      REAL(fp)                :: LOSS_FREQ
      INTEGER               :: I, J
      REAL(fp),SAVE           :: TABLE_S(21), TABLE_U10(20)
      REAL(fp),SAVE           :: TABLE_UPTAKE(21,20) 
      REAL(fp)                :: SFCWINDSQR

      ! Flag for first call
      LOGICAL, SAVE         :: FIRST=.TRUE.

      !=================================================================
      ! HG2_SEASALT_LOSSRATE begins here!
      !=================================================================


!$OMP PARALLEL DO 
!$OMP+DEFAULT( SHARED ) 
!$OMP+PRIVATE( I, J, U10M, S, LOSS_FREQ, SFCWINDSQR )
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Only calculate deposition via sea salt over water
         IF ( IS_WATER( I, J, State_Met ) ) THEN

            ! Wind speed at 10m altitude [m/s]
            SFCWINDSQR = State_Met%U10M(I,J)**2 
     &                 + State_Met%V10M(I,J)**2
            U10M       = SQRT( SFCWINDSQR )

            ! Don't allow wind >20 m/s which is the limit of this 
            ! parameterization
            U10M       = MAX( MIN( U10M, 20e+0_fp ), 1e+0_fp )

            ! Relative humidity as a saturation ratio
            ! Use the relative humidity of the lowest layer, although this is 
            ! lower than the higher layers of the MBL
            !
            ! Don't allow supersaturation, as [Cl-] is undefined for RH>=1 
            ! Cap RH at 99%, Don't allow RH < 75% as happens in coastal areas
            S = MAX( MIN( State_Met%RH(I,J,1), 99e+0_fp ), 75e+0_fp ) * 
     &             1e-2_fp
            
            LOSS_FREQ = 1e-10_fp * ( 1e+0_fp - EXP( -57.758e+0_fp * 
     &           (1e+0_fp-S) ) ) *
     &           EXP( -1.9351e+0_fp  * U10M + 
     &                 9.0047e+0_fp  * SQRT( U10M ) + 
     &                 0.14788e+0_fp * U10M**1.5e+0_fp ) 
            
            ! Loss frequency must be positive
            LOSS_FREQ = MAX( LOSS_FREQ, 1e-10_fp ) 

         ELSE 

            ! No loss over land
            LOSS_FREQ = 0e+0_fp
            
         ENDIF
         
         HG2_SEASALT_LOSSRATE(I,J) = LOSS_FREQ


      ENDDO ! I 
      ENDDO ! J
!$OMP END PARALLEL DO

      ! Return to calling program
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
      FUNCTION GET_JNO2( I, J, L, State_Met ) RESULT( JNO2_NOW )
!
! !USES:
!
      USE CMN_SIZE_MOD
      USE GC_GRID_MOD,        ONLY : GET_YMID_R
      USE PhysConstants
      USE State_Met_Mod,      ONLY : MetState
      USE TIME_MOD,           ONLY : GET_TS_CHEM, GET_DAY_OF_YEAR 
!
! !INPUT PARAMETERS: 
!
      INTEGER,        INTENT(IN) :: I           ! Longitude index
      INTEGER,        INTENT(IN) :: J           ! Latitude index
      INTEGER,        INTENT(IN) :: L           ! Level index
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
!  01 Mar 2012 - R. Yantosca - Now use GET_YMID_R(I,J,L) from grid_mod.F90
!  28 Nov 2012 - R. Yantosca - Now pass State_Met via the argument list
!  28 Nov 2012 - R. Yantosca - Replace SUNCOS with State_Met%SUNCOS
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
      IF ( State_Met%SUNCOS(I,J) > TINY(1e+0_fp) ) THEN

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
         DEC = A0 - A1*COS(     R ) + B1*SIN(     R )
     &            - A2*COS( 2e+0_fp*R ) + B2*SIN( 2e+0_fp*R )
     &            - A3*COS( 3e+0_fp*R ) + B3*SIN( 3e+0_fp*R )

         ! Latitude of grid box [radians]
         YMID_R = GET_YMID_R( I, J, L )

         ! Cosine of solar zenith angle at local noon
         CSZANOON = SIN( YMID_R ) * SIN( DEC ) +
     &        COS( YMID_R ) * COS( DEC ) 

         ! Parrish et al (1983) found 
         ! J-NO2 ~ exp( -0.360 * sec(SZA) ), so
         ! J-NO2(now) = J-NO2(noon) *  exp( -0.360 * sec(SZANOW)  ) /
         !                             exp( -0.360 * sec(SZANOON) )
         !            = J-NO2(noon) * exp( 0.360 * [sec(SZANOON)-sec(SZANOW)] )

         ! Impose a diurnal variation on JNO2 during the day
         ! Note: We don't need to check for divide-by-zero errors
         ! because we already checked SUNCOS(I,J) and we know
         ! CSZANOON >= SUNCOS(I,J)
         JNO2_NOW = JNO2(I,J,L) * 
     &     EXP( 0.36e+0_fp * 
     &        ( 1e+0_fp/CSZANOON - 1e+0_fp/State_Met%SUNCOS(I,J) ) )

         ! Make sure OH is not negative
         JNO2_NOW= MAX( JNO2_NOW, 0e+0_fp )
               
      ELSE

         ! At night, JNO2 goes to zero
         JNO2_NOW = 0e+0_fp

      ENDIF

      ! Return to calling program
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
      SUBROUTINE PARTITIONHG2( am_I_Root, Input_Opt, State_Chm, RC )
!
! !USES:
!
      USE CMN_SIZE_MOD
      USE ErrCode_Mod
      USE ERROR_MOD,          ONLY : SAFE_DIV  , GEOS_CHEM_STOP
      USE Input_Opt_Mod,      ONLY : OptInput
      USE OCEAN_MERCURY_MOD,  ONLY : Fg
      USE State_Chm_Mod,      ONLY : ChmState
      USE State_Met_Mod,      ONLY : MetState
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)    :: am_I_Root   ! Are we on the root CPU?
      TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
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
!  02-Jan-2011 - H. Amos     - Add/clean up comments
!  04 Jan 2012 - H. Amos     - modify algorithms to reflect the fact that 
!                              anthropogenic Hg(p) is now emitted as Hg(II) 
!                              (i.e. it's no longer considered refractory). 
!  28-Mar-2013 - C. Holmes   - Since we now assume that all HgP and Hg2 
!                              actively partition between gas and aerosol, 
!                              we no longer need to 'reverse' partition. 
!  03 Jun 2013 - R. Yantosca - Now accept am_I_Root, Input_Opt, State_Chm, RC
!  18 Nov 2014 - M. Yannetti - Took out OCEAN_MERCURY_MOD calling Fp
!  30 Jun 2016 - R. Yantosca - Remove instances of STT.  Now get the advected
!                              species ID from State_Chm%Map_Advect.
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

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, N, FGas, Hg2TOT )
      DO L=1, LLPAR
      DO J=1, JJPAR
      DO I=1, IIPAR
         
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
            HG2TOT = Spc(I,J,L,Hg2_Id_List(N)) + 
     &               Spc(I,J,L,HgP_Id_List(N))

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

      ! Return to calling program
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
      INTEGER,INTENT(IN) :: I, J
      REAL(fp), INTENT(IN) :: E_plant
! 
! !REVISION HISTORY: 
!  11 Jan 2011 - Y. Zhang    - Initial version
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
      SUBROUTINE DO_RED_INPLUME( am_I_Root, Input_Opt, RC )
!
! !USES:
!
      USE CMN_SIZE_MOD
      USE ErrCode_Mod
      USE ERROR_MOD,          ONLY : ERROR_STOP
      USE GC_GRID_MOD,        ONLY : GET_AREA_M2
      USE HCO_INTERFACE_MOD,  ONLY : HcoState
      USE HCO_EmisList_Mod,   ONLY : HCO_GetPtr
      USE Input_Opt_Mod,      ONLY : OptInput
      USE TIME_MOD,           ONLY : EXPAND_DATE

      IMPLICIT NONE
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)  :: am_I_Root   ! Are we on the root CPU?
      TYPE(OptInput), INTENT(IN)  :: Input_Opt   ! Input Options object
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT) :: RC          ! Success or failure?
!
! !REVISION HISTORY: 
!  11 Jan 2011 - Y. Zhang    - Initial version
!  23 Jun 2014 - R. Yantosca - Now accept am_I_Root, Input_Opt, RC
!  06 Nov 2014 - R. Yantosca - Replace TRANSFER_2D with direct casts
!  17 Dec 2014 - R. Yantosca - Leave time/date variables as 8-byte
!  05 Mar 2015 - R. Yantosca - Add Input_Opt%RES_DIR to data path'
!  16 Mar 2015 - R. Yantosca - Remove bpch input
!  29 Apr 2016 - R. Yantosca - Don't initialize pointers in declaration stmts
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
      CHARACTER(LEN=255)  :: LOC = 'DO_RED_INPLUME (mercury_mod.F)'

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
      E_ELEC_Hg2 => NULL()

      
      ! Get a pointer to the monthly mean OH from HEMCO (bmy, 3/11/15)
      CALL HCO_GetPtr( am_I_Root, HcoState, 
     &                'CFPP_NEI2005_Hg2', E_ELEC_Hg2, RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL ERROR_STOP ( 
     &      'Cannot get pointer to CFPP_NEI2005_Hg2!', LOC )
      ENDIF

      ! Do the reduction
      DO J = 1, JJPAR
      DO I = 1, IIPAR
             
         ! Hg emission in CFPP sector decreased for 45.14% during 2005-2010 
         ! Convert [kg/m2/s] --> [kg/s]
         E_plant = E_ELEC_Hg2(I,J) * GET_AREA_M2(I,J,1)

         ! Reduce the Hg2 from plume
         CALL RED_INPLUME_GRID( I, J, E_plant )

      ENDDO
      ENDDO

      ! Free npointer
      E_ELEC_Hg2 => NULL()

      ! Return to calling program
      END SUBROUTINE DO_RED_INPLUME
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
      SUBROUTINE INIT_MERCURY( am_I_Root, Input_Opt, State_Chm, RC )
!
! !USES:
!
      USE CMN_DIAG_MOD
      USE CMN_SIZE_MOD
      USE ErrCode_Mod
      USE ERROR_MOD,          ONLY : ALLOC_ERR, ERROR_STOP
      USE Input_Opt_Mod,      ONLY : OptInput
      USE Species_Mod,        ONLY : Species
      USE State_Chm_Mod,      ONLY : ChmState
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)    :: am_I_Root   ! Are we on the root CPU?
      TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
      TYPE(ChmState), INTENT(IN)    :: State_Chm   ! Chemistry State object
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
! 
! !REVISION HISTORY: 
!  02 Dec 2004 - N. (Eckley) Selin - Initial version
!  (1 ) Removed reference to FEMIS array.  Now also allocates and zeroes
!        the T44 array.  Added reference to CMN_DIAG.  Now references 
!        N_TRACERS from "tracer_mod.f". (bmy, 2/24/05)
!  (2 ) EHg0_an, EHg2_an, EHgP_an are now 2-D arrays.  Now modified for 
!        updated ocean mercury module. (eck, cdh, sas, bmy, 4/6/06)
!  27 Sep 2011 - H. Amos - remove LHg_WETDasHNO3 logical, it's obsolete
!  11 Apr 2012 - R. Yantosca - Now retire lai_mod.F; remove call to INIT_LAI
!  11 Apr 2012 - R. Yantosca - Now move LAI initialization to main.F
!  03 Jun 2013 - R. Yantosca - Now use fields from Input_Opt
!  12 Mar 2015 - R. Yantosca - Now test if HEMCO_Collections are defined
!  25 Apr 2016 - R. Yantosca - Now use species database for tagHg indexing
!  29 Apr 2016 - R. Yantosca - Don't initialize pointers in declaration stmts
!  01 Jul 2016 - R. Yantosca - Now rename species DB object ThisSpc to SpcInfo
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
      INTEGER                :: nAdvect
      INTEGER                :: AS, N
      CHARACTER(LEN=255)     :: LOCATION

      ! Pointers
      TYPE(Species), POINTER :: SpcInfo

      !=================================================================
      ! INIT_MERCURY begins here!
      !=================================================================

      ! Assume success
      RC        =  GC_SUCCESS

      ! Return if we have already allocated arrays
      IF ( IS_INIT ) RETURN

      ! Initialize 
      SpcInfo        => NULL()
      LDRYD          = Input_Opt%LDRYD            ! Use drydep?
      LGTMM          = Input_Opt%LGTMM            ! Use GTMM model?
      LNLPBL         = Input_Opt%LNLPBL           ! Use non-local PBL?
      LSPLIT         = Input_Opt%LSPLIT           ! Tagged simulation?
      nAdvect        = State_Chm%nAdvect          ! # of Hg advected species

      ! Location string for ERROR_STOP
      LOCATION       = 'DEFINE_TAGGED_Hg ("mercury_mod.f")'

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
      ALLOCATE( COSZM( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'COSZM' )
      COSZM = 0e+0_fp

      ALLOCATE( EHg0_an( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EHg0_an' )
      EHg0_an = 0e+0_fp

      ALLOCATE( EHg0_am( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EHg0_am' )
      EHg0_am = 0e+0_fp

      ALLOCATE( EHg2_an( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EHg2_an' )
      EHg2_an = 0e+0_fp

      ALLOCATE( EHgP_an( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EHgP_an' )
      EHgP_an = 0e+0_fp

      ALLOCATE( EHg0_oc( IIPAR, JJPAR, N_Hg_CATS ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EHg0_oc' )
      EHg0_oc = 0e+0_fp

      ALLOCATE( EHg0_dist( IIPAR, JJPAR), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EHg0_dist' )
      EHg0_dist = 0e+0_fp

      ALLOCATE( EHg0_geo( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EHg0_geo' )
      EHg0_geo = 0e+0_fp

      ALLOCATE( EHg0_bb( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EHg0_bb' )
      EHg0_bb = 0e+0_fp

      ALLOCATE( EHg0_snow( IIPAR, JJPAR, N_Hg_CATS ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EHg0_snow' )
      EHg0_snow = 0e+0_fp

      ! Allocate the following if tagged Hg simulation
      IF ( LSPLIT ) THEN

         ! Tagged Hg0 arrays (eds 8/31/10)
         ALLOCATE( EHg0_can( IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'EHg0_can' )
         EHg0_can = 0e+0_fp

         ALLOCATE( EHg0_usa( IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'EHg0_usa' )
         EHg0_usa = 0e+0_fp

         ALLOCATE( EHg0_cam( IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'EHg0_cam' )
         EHg0_cam = 0e+0_fp

         ALLOCATE( EHg0_sam( IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'EHg0_sam' )
         EHg0_sam = 0e+0_fp

         ALLOCATE( EHg0_waf( IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'EHg0_waf' )
         EHg0_waf = 0e+0_fp

         ALLOCATE( EHg0_eaf( IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'EHg0_eaf' )
         EHg0_eaf = 0e+0_fp

         ALLOCATE( EHg0_saf( IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'EHg0_saf' )
         EHg0_saf = 0e+0_fp

         ALLOCATE( EHg0_naf( IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'EHg0_naf' )
         EHg0_naf = 0e+0_fp

         ALLOCATE( EHg0_eur( IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'EHg0_eur' )
         EHg0_eur = 0e+0_fp

         ALLOCATE( EHg0_eeu( IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'EHg0_eeu' )
         EHg0_eeu = 0e+0_fp

         ALLOCATE( EHg0_mde( IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'EHg0_mde' )
         EHg0_mde = 0e+0_fp

         ALLOCATE( EHg0_sov( IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'EHg0_sov' )
         EHg0_sov = 0e+0_fp

         ALLOCATE( EHg0_sas( IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'EHg0_sas' )
         EHg0_sas = 0e+0_fp

         ALLOCATE( EHg0_eas( IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'EHg0_eas' )
         EHg0_eas = 0e+0_fp

         ALLOCATE( EHg0_sea( IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'EHg0_sea' )
         EHg0_sea = 0e+0_fp

         ALLOCATE( EHg0_jpn( IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'EHg0_jpn' )
         EHg0_jpn = 0e+0_fp

         ALLOCATE( EHg0_oce( IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'EHg0_oce' )
         EHg0_oce = 0e+0_fp

         ! Tagged Hg2 arrays (eds 8/31/10)
         ALLOCATE( EHg2_can( IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'EHg2_can' )
         EHg2_can = 0e+0_fp

         ALLOCATE( EHg2_usa( IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'EHg2_usa' )
         EHg2_usa = 0e+0_fp

         ALLOCATE( EHg2_cam( IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'EHg2_cam' )
         EHg2_cam = 0e+0_fp

         ALLOCATE( EHg2_sam( IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'EHg2_sam' )
         EHg2_sam = 0e+0_fp

         ALLOCATE( EHg2_waf( IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'EHg2_waf' )
         EHg2_waf = 0e+0_fp

         ALLOCATE( EHg2_eaf( IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'EHg2_eaf' )
         EHg2_eaf = 0e+0_fp

         ALLOCATE( EHg2_saf( IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'EHg2_saf' )
         EHg2_saf = 0e+0_fp

         ALLOCATE( EHg2_naf( IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'EHg2_naf' )
         EHg2_naf = 0e+0_fp

         ALLOCATE( EHg2_eur( IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'EHg2_eur' )
         EHg2_eur = 0e+0_fp

         ALLOCATE( EHg2_eeu( IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'EHg2_eeu' )
         EHg2_eeu = 0e+0_fp

         ALLOCATE( EHg2_mde( IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'EHg2_mde' )
         EHg2_mde = 0e+0_fp

         ALLOCATE( EHg2_sov( IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'EHg2_sov' )
         EHg2_sov = 0e+0_fp

         ALLOCATE( EHg2_sas( IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'EHg2_sas' )
         EHg2_sas = 0e+0_fp

         ALLOCATE( EHg2_eas( IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'EHg2_eas' )
         EHg2_eas = 0e+0_fp

         ALLOCATE( EHg2_sea( IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'EHg2_sea' )
         EHg2_sea = 0e+0_fp

         ALLOCATE( EHg2_jpn( IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'EHg2_jpn' )
         EHg2_jpn = 0e+0_fp

         ALLOCATE( EHg2_oce( IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'EHg2_oce' )
         EHg2_oce = 0e+0_fp

         ! Tagged HgP arrays (eds 8/31/10)
         ALLOCATE( EHgP_can( IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'EHgP_can' )
         EHgP_can = 0e+0_fp

         ALLOCATE( EHgP_usa( IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'EHgP_usa' )
         EHgP_usa = 0e+0_fp

         ALLOCATE( EHgP_cam( IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'EHgP_cam' )
         EHgP_cam = 0e+0_fp

         ALLOCATE( EHgP_sam( IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'EHgP_sam' )
         EHgP_sam = 0e+0_fp

         ALLOCATE( EHgP_waf( IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'EHgP_waf' )
         EHgP_waf = 0e+0_fp

         ALLOCATE( EHgP_eaf( IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'EHgP_eaf' )
         EHgP_eaf = 0e+0_fp

         ALLOCATE( EHgP_saf( IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'EHgP_saf' )
         EHgP_saf = 0e+0_fp

         ALLOCATE( EHgP_naf( IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'EHgP_naf' )
         EHgP_naf = 0e+0_fp

         ALLOCATE( EHgP_eur( IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'EHgP_eur' )
         EHgP_eur = 0e+0_fp

         ALLOCATE( EHgP_eeu( IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'EHgP_eeu' )
         EHgP_eeu = 0e+0_fp

         ALLOCATE( EHgP_mde( IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'EHgP_mde' )
         EHgP_mde = 0e+0_fp

         ALLOCATE( EHgP_sov( IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'EHgP_sov' )
         EHgP_sov = 0e+0_fp

         ALLOCATE( EHgP_sas( IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'EHgP_sas' )
         EHgP_sas = 0e+0_fp

         ALLOCATE( EHgP_eas( IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'EHgP_eas' )
         EHgP_eas = 0e+0_fp

         ALLOCATE( EHgP_sea( IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'EHgP_sea' )
         EHgP_sea = 0e+0_fp

         ALLOCATE( EHgP_jpn( IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'EHgP_jpn' )
         EHgP_jpn = 0e+0_fp

         ALLOCATE( EHgP_oce( IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'EHgP_oce' )
         EHgP_oce = 0e+0_fp

      ENDIF

      IF ( LGTMM ) THEN

         ALLOCATE( EHg0_gtm( IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'EHg0_gtm' )
         EHg0_gtm = 0e+0_fp

      ELSE

         ALLOCATE( EHg0_ln( IIPAR, JJPAR, N_Hg_CATS ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'EHg0_ln' )
         EHg0_ln = 0e+0_fp

         ALLOCATE( EHg0_vg( IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'EHg0_vg' )
         EHg0_vg = 0e+0_fp

         ALLOCATE( EHg0_so( IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'EHg0_so' )
         EHg0_so = 0e+0_fp

      ENDIF

      ALLOCATE( TCOSZ( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'TCOSZ' )
      TCOSZ = 0e+0_fp

      ALLOCATE( TTDAY( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'TTDAY' )
      TTDAY = 0e+0_fp
 
      ! Allocate ZERO_DVEL if we use non-local PBL mixing or
      ! if drydep is turned off 
      IF ( LNLPBL .OR. (.not. LDRYD) ) THEN
         ALLOCATE( ZERO_DVEL( IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'ZERO_DVEL' )
         ZERO_DVEL = 0e+0_fp
      ENDIF

      ALLOCATE( HG2_SEASALT_LOSSRATE( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'HG2_SEASALT_LOSSRATE' )
      HG2_SEASALT_LOSSRATE = 0e+0_fp

      !=================================================================
      ! Allocate & initialize arrays for tagged species indices
      !=================================================================
      IF ( LSPLIT ) THEN

         ! Species indices for tagged anthro regions
         ALLOCATE( AN_Hg0( IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'AN_Hg0' )
         AN_Hg0 = 0e+0_fp

         ALLOCATE( AN_Hg2( IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'AN_Hg2' )
         AN_Hg2 = 0e+0_fp

         ALLOCATE( AN_HgP( IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'AN_HgP' )
         AN_HgP = 0e+0_fp

      ENDIF

      ! HG_EMIS is needed for non-local PBL mixing 
      ALLOCATE ( HG_EMIS( IIPAR, JJPAR, nAdvect ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'HG_EMIS' )
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

      ! Switch uses Bess' best-guess emissions if TRUE or emissions from
      ! Selin et al. 2008 (+bug fixes) if FALSE
      ! This is temporary
      ! Bess's changes: rescale GEIA to use Streets' 2006 regional totals
      ! With these emissions, Bess found the best chemistry results
      LGCAPEMIS=.TRUE.

      ! Use GEIA 2005 inventory
      LGEIA05=.FALSE.

      ! Switch use NEI2005 and NPRI2005 emission inventories
      LNEI2005=.TRUE.

      ! Switch modifying the speciation profile of Hg emission
      LInPlume=.FALSE.

      ! no Hg emitted through transpiration (VEGEMIS off)
      LVEGEMIS=.FALSE.
     
      ! Switch adds bromine explosion in Northern springtime
      LPOLARBR=.TRUE.

      ! Switch turns on Hg+Br chemistry
      LBRCHEM=.TRUE.

      ! Switch turns on Hg+BrO chemistry
      LBROCHEM=.FALSE.

      ! Switch turn on Hg+OH, Hg+O3 chemistry
      LOHO3CHEM=.FALSE.

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

      ! Return to calling program
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
      SUBROUTINE SET_OPTIONS_FROM_HEMCO( am_I_Root, Input_Opt, RC )
!
! !USES:
!
      USE ERROR_MOD,         ONLY : ERROR_STOP
      USE HCO_ERROR_MOD
      USE HCO_INTERFACE_MOD, ONLY : HcoState
      USE HCO_ExtList_Mod,   ONLY : GetExtOpt
      USE Input_Opt_Mod,     ONLY : OptInput
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)  :: am_I_Root
      TYPE(OptInput), INTENT(IN)  :: Input_Opt
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT) :: RC   ! Success or failure?
!
! !REMARKS:
!
! !REVISION HISTORY:
!  06 Jan 2015 - R. Yantosca - Initial version
!  11 Oct 2016 - R. Yantosca - Move call to DO_RED_INPLUME here
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
      CHARACTER(LEN=255) :: LOC='SET_OPTIONS_FROM_HEMCO (mercury_mod.F)'

      !-----------------------------------------------------------------
      ! Set the value of LRED_JNO2 depending on whether or not
      ! the HEMCO collection LRED_JNO2 is activated
      !-----------------------------------------------------------------
      CALL GetExtOpt( HcoState%Config, -999, 'LRED_JNO2',  
     &                OptValBool=LRED_JNO2, FOUND=FOUND, RC=RC )
      IF ( RC /= HCO_SUCCESS ) THEN
         CALL ERROR_STOP( 'GetExtOpt LRED_JNO2', LOC )
      ENDIF
      IF ( .not. FOUND ) THEN
         LRED_JNO2 = .FALSE.
      ENDIF

      !-----------------------------------------------------------------
      ! Set the value of LGCBROMINE depending on the values of the
      ! HEMCO collections BrOx_GC and BrOx_TOMCAT
      !-----------------------------------------------------------------

      ! First look for BrOx_GC
      CALL GetExtOpt( HcoState%Config, -999, 'BrOx_GC', 
     &                OptValBool=LGC, FOUND=FOUND, RC=RC )
      IF ( RC /= HCO_SUCCESS ) THEN
         CALL ERROR_STOP( 'GetExtOpt BrOx_GC', LOC )
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
      CALL GetExtOpt( HcoState%Config, -999, 'NEI2005',  
     &                OptValBool=LNEI2005, FOUND=FOUND, RC=RC )
      IF ( RC /= HCO_SUCCESS ) THEN
         CALL ERROR_STOP( 'GetExtOpt LNEI2005', LOC )
      ENDIF
      IF ( .not. FOUND ) THEN
         LNEI2005 = .FALSE.
      ENDIF

      !-----------------------------------------------------------------
      ! Set the value of LInPlume depending on whether or not
      ! the HEMCO collection NEI2005 is activated
      !-----------------------------------------------------------------
      CALL GetExtOpt( HcoState%Config, -999, 'LRED_INPLUME',  
     &                OptValBool=LInPlume, FOUND=FOUND, RC=RC )
      IF ( RC /= HCO_SUCCESS ) THEN
         CALL ERROR_STOP( 'GetExtOpt LInPlume5', LOC )
      ENDIF
      IF ( .not. FOUND ) THEN
         LInPlume = .FALSE.
      ENDIF

      ! In plume degradation of Hg2 by SO2 in U.S. and Canada at CFPPs, 
      ! (yzh 11/1/2011).  Move this here so that the HEMCO_Config file
      ! will have been already read by this point. (bmy, 10/11/16)
      IF ( LInPlume .AND. .NOT. LPREINDHG .AND. LNEI2005 ) THEN
         CALL DO_RED_INPLUME( am_I_Root, Input_Opt, RC )
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
      WRITE( 6, '(a)' ) REPEAT( '=', 79 )
 100  FORMAT( 'Adjusting Hg simulation settings from HEMCO inputs' )
 110  FORMAT( 'LRED_JNO2  is set to ', L1                          )
 120  FORMAT( 'LGCBROMINE is set to ', L1                          )
 130  FORMAT( 'LNEI2005   is set to ', L1                          )
 140  FORMAT( 'LInPlume   is set to ', L1                          )
      
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
!  (1 ) Now deallocate MLD, NPP, RAD (sas, bmy, 1/18/05)
!  (2 ) Now deallocate T44 (bmy, 2/24/05)
!EOP
!------------------------------------------------------------------------------
!BOC
!
      !=================================================================
      ! CLEANUP_MERCURY begins here!
      !=================================================================
      IF ( ALLOCATED( AN_Hg0   ) ) DEALLOCATE( AN_Hg0  )  
      IF ( ALLOCATED( AN_Hg2   ) ) DEALLOCATE( AN_Hg2  )      
      IF ( ALLOCATED( AN_HgP   ) ) DEALLOCATE( AN_HgP  )      
      IF ( ALLOCATED( COSZM    ) ) DEALLOCATE( COSZM   )      
      IF ( ALLOCATED( EHg0_an  ) ) DEALLOCATE( EHg0_an )
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
      IF ( ALLOCATED( EHg0_snow) ) DEALLOCATE( EHg0_snow )
      IF ( ALLOCATED( EHg0_geo ) ) DEALLOCATE( EHg0_geo )
      IF ( ALLOCATED( EHg0_bb  ) ) DEALLOCATE( EHg0_bb )
      IF ( ALLOCATED( EHg0_gtm  ) ) DEALLOCATE( EHg0_gtm )
      IF ( ALLOCATED( EHg0_vg  ) ) DEALLOCATE( EHg0_vg )
      IF ( ALLOCATED( EHg0_so  ) ) DEALLOCATE( EHg0_so )
      IF ( ALLOCATED( EHg0_dist) ) DEALLOCATE( EHg0_dist )
      IF ( ALLOCATED( TCOSZ    ) ) DEALLOCATE( TCOSZ   )
      IF ( ALLOCATED( TTDAY    ) ) DEALLOCATE( TTDAY   )
!      IF ( ALLOCATED( TRANSP   ) ) DEALLOCATE( TRANSP  )
      IF ( ALLOCATED( ZERO_DVEL) ) DEALLOCATE( ZERO_DVEL )
      IF ( ALLOCATED( HG2_SEASALT_LOSSRATE ) ) 
     &     DEALLOCATE( HG2_SEASALT_LOSSRATE   ) 
      IF ( ALLOCATED( HG_EMIS   ) ) DEALLOCATE( HG_EMIS   ) 
  
      ! Free pointers to HEMCO fields
      O3    => NULL()
      OH    => NULL()
      JNO2  => NULL()

      ! Free Hg indexing pointers
      Hg0_Id_List => NULL()
      Hg2_Id_List => NULL()
      HgP_Id_List => NULL()

      END SUBROUTINE CLEANUP_MERCURY
!EOC
      END MODULE MERCURY_MOD
