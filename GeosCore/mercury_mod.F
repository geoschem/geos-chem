! $Id: mercury_mod.f,v 1.24 2009/09/01 19:21:18 cdh Exp $
      MODULE MERCURY_MOD
!
!******************************************************************************
!  Module MERCURY_MOD contains variables and routines for the GEOS-CHEM 
!  mercury simulation. (eck, bmy, 12/14/04, 4/6/06)
!
!  Many choices of reaction mechanism and model processes can be selected with
!  logical switches located in INIT_MERCURY. (cdh, 11/25/09)
!
!  Module Variables:
!  ============================================================================
!  (1 ) AN_Hg0   (INTEGER) : Tracer index array for tagged anth Hg(0)  regions
!  (2 ) AN_Hg2   (INTEGER) : Tracer index array for tagged anth Hg(II) regions
!  (3 ) AN_HgP   (INTEGER) : Tracer index array for tagged anth HgP    regions
!  (4 ) COSZM    (REAL*8 ) : Max daily value of COS( S. Z. Angle ) [unitless]
!  (5 ) DRYHg0   (INTEGER) : Index for Hg0 in DEPSAV array (drydep freq)
!  (6 ) DRYHg2   (INTEGER) : Index for Hg2 in DEPSAV array (drydep freq)
!  (7 ) DRYHgP   (INTEGER) : Index for HgP in DEPSAV array (drydep freq)
!  (8 ) EHg0_an  (REAL*8 ) : Anthropogenic Hg0 emissions [kg/s]
!  (9 ) EHg2_an  (REAL*8 ) : Anthropogenic Hg2 emissions [kg/s]
!  (10) EHgP_an  (REAL*8 ) : Anthropogenic HgP emissions [kg/s]
!  (11) EHg0_am  (REAL*8 ) : Artisanal mining Hg0 emissions [kg/s]
!  (12) EHg0_oc  (REAL*8 ) : Hg0 emissions from oceans [kg/s]
!  (13) EHg0_ln  (REAL*8 ) : Re-emissions of Hg0 from land [kg/s]
!  (14) EHg0_nt  (REAL*8 ) : Hg0 emissions from natural land sources [kg/s] 
!  (15) EHg0_bb  (REAL*8 ) : Hg0 emissions fom biomass burning sources [kg/s]
!  (16) EHg0_vg  (REAL*8 ) : Hg0 emissions fom vegetation sources [kg/s] 
!  (17) EHg0_so  (REAL*8 ) : Hg0 emissions from soil sources [kg/s]
!  (18) EHg0_dist(REAL*8 ) : Spatial distribution of terrestrial Hg0 sources
!                              [dimensionless] 
!  (19) TCOSZ    (REAL*8 ) : Sum of COS( Solar Zenith Angle ) [unitless]
!  (20) TTDAY    (REAL*8 ) : Total daylight time at location (I,J) [minutes]
!  (21) T44      (REAL*4 ) : Local array for drydep diagnostic
!  (22) ZERO_DVEL(REAL*8 ) : Array with zero dry deposition velocity [cm/s]
!  (23) ANTHRO_Hg_YEAR(INT): Anthropogenic Hg emissions year (1995 or 2000)
!  (24) TRANSP   (REAL*8 ) : Plant transpiration rate [m/s]
!  (25) SMALLNUM (REAL*8 ) : A small number to prevent underflow
!
!  Module Routines:
!  ===========================================================================
!  (1 ) CHEMMERCURY        : Chemistry routine for Hg
!  (2 ) CHEM_Hg0_Hg2       : Chemistry for Hg0, Hg2 and drydep of Hg2
!  (3 ) RXN_REDOX_NODEP    : Redox chemistry of Hg(0) and Hg(II), no deposition
!  (4 ) RXN_REDOX_WITHDEP  : Redox chemistry and deposition of Hg(0) and Hg(II)
!  (5 ) CHEM_HGP           : Chemistry (via drydep loss) for HgP
!  (6 ) RXN_HgP_DRYD       : Loss of HgP via drydep
!  (7 ) EMISSMERCURY       : Emissions for Hg
!  (8 ) BIOMASSHG          : Wildfire Hg emissions
!  (9 ) VEGEMIS            : Transpiration Hg emissions 
!  (10) SOILEMIS           : Soil Hg emissions
!  (11) READ_NASA_TRANSP   : Read transpiration rates from file
!  (12) SRCHG0             : Applies emissions of Hg0
!  (13) SRCHG2             : Applies emissions of Hg2
!  (14) SRCHGP             : Applies emissions of HgP
!  (15) MERCURY_READYR     : Reads mercury emissions and converts to [kg/s]
!  (16) GET_LWC            : Computes liquid water content as a function of T
!  (17) GET_VCLDF          : Computes volume cloud fraction 
!  (18) GET_O3             : Returns monthly mean O3 field
!  (19) GET_OH             : Returns monthly mean OH field (w/ diurnal scaling)
!  (20) OHNO3TIME          : Computes diurnal scaling for monthly mean OH
!  (21) DEFINE_TAGGED_Hg   : Defines tracer number for tagged Hg tracers 
!  (22) INIT_MERCURY       : Allocates and zeroes module arrays
!  (23) CLEANUP_MERCURY    : Deallocates module arrays
!
!  GEOS-CHEM modules referenced by mercury_mod.f
!  ============================================================================
!  (1 ) bpch2_mod.f        : Module w/ routines for binary pch file I/O
!  (2 ) comode_mod.f       : Module w/ SMVGEAR allocatable arrays
!  (3 ) dao_mod.f          : Module w/ arrays for DAO met fields
!  (4 ) diag_mod.f         : Module w/ GEOS-CHEM diagnostic arrays
!  (5 ) drydep_mod.f       : Module w/ GEOS-CHEM dry deposition routines
!  (6 ) error_mod.f        : Module w/ NaN, other error check routines
!  (7 ) global_o3_mod.f    : Module w/ routines to read 3-D O3 field
!  (8 ) global_oh_mod.f    : Module w/ routines to read 3-D OH field
!  (9 ) global_br_mod.f    : Module w/ routines to read 3-D Br field
!  (10) grid_mod.f         : Module w/ horizontal grid information
!  (11) logical_mod.f      : Module w/ GEOS-CHEM logical switches
!  (12) pbl_mix_mod.f      : Module w/ routines for PBL height & mixing
!  (13) pressure_mod.f     : Module w/ routines to compute P(I,J,L)
!  (14) regrid_a2a_mod.f   : Module w/ routines to regrid data
!  (15) time_mod.f         : Module w/ routines to compute date & time
!  (16) tracer_mod.f       : Module w/ GEOS-CHEM tracer array STT etc.
!  (17) tracerid_mod.f     : Module w/ pointers to tracers & emissions
!  (18) transfer_mod.f     : Module w/ routines to cast & resize arrays
!
!  Nomenclature: 
!  ============================================================================
!  (1 ) Hg(0)  a.k.a. Hg0  : Elemental   mercury
!  (2 ) Hg(II) a.k.a. Hg2  : Divalent    mercury
!  (3 ) HgP                : Particulate mercury
!
!  Mercury Tracers (1-3 are always defined; 4-21 are defined for tagged runs)
!  ============================================================================
!  (1 ) Hg(0)              : Hg(0)  - total tracer
!  (2 ) Hg(II)             : Hg(II) - total tracer 
!  (3 ) HgP                : HgP    - total tracer
!  ------------------------+---------------------------------------------------
!  (4 ) Hg(0)_an           : Hg(0)  - North American anthropogenic
!  (5 ) Hg(0)_ae           : Hg(0)  - European Anthropogenic
!  (6 ) Hg(0)_aa           : Hg(0)  - Asian anthropogenic
!  (7 ) Hg(0)_ar           : Hg(0)  - Rest of World Anthropogenic
!  (8 ) Hg(0)_oc           : Hg(0)  - Ocean emission
!  (9 ) Hg(0)_ln           : Hg(0)  - Land reemission
!  (10) Hg(0)_nt           : Hg(0)  - Land natural emission
!  ------------------------+---------------------------------------------------
!  (11) Hg(II)_an          : Hg(II) - North American anthropogenic
!  (12) Hg(II)_ae          : Hg(II) - European Anthropogenic
!  (13) Hg(II)_aa          : Hg(II) - Asian anthropogenic
!  (14) Hg(II)_ar          : Hg(II) - Rest of World Anthropogenic
!  (15) Hg(II)_oc          : Hg(II) - Ocean emission
!  (16) Hg(II)_ln          : Hg(II) - Land reemission
!  (17) Hg(II)_nt          : Hg(II) - Land natural emission
!  ------------------------+---------------------------------------------------
!  (18) HgP_an             : HgP    - North American anthropogenic
!  (19) HgP_ae             : HgP    - European anthropogenic
!  (20) HgP_aa             : HgP    - Asian anthropogenic
!  (21) HgP_ar             : HgP    - Rest of world anthropogenic
!  ------------------------+---------------------------------------------------
!  (22) HgP_oc             : HgP    - Ocean emission        (FOR FUTURE USE)
!  (23) HgP_ln             : HgP    - Land reemission       (FOR FUTURE USE)
!  (24) HgP_nt             : HgP    - Land natural emission (FOR FUTURE USE)
!
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
!       Environ. Sci. Technol., 44, 8574-8580, 2010
!
!  NOTES:
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
!  08 Feb 2012 - R. Yantosca - Add modifications for GEOS-5.7.x met
!  01 Mar 2012 - R. Yantosca - Now reference new grid_mod.F90
!  01 Mar 2012 - R. Yantosca - Use updated GET_LOCALTIME from time_mod.F
!******************************************************************************
!
      ! References to F90 modules
      USE OCEAN_MERCURY_MOD, ONLY : LDYNSEASALT, LGCAPEMIS, LPOLARBR
      !Added 2005 GEIA inventory flag (eck, 10/19/11)
      USE OCEAN_MERCURY_MOD, ONLY : LGEIA05
      ! Added LVEGEMIS, distinct from LGCAPEMIS (jaf, eds, 4/1/11)
      USE OCEAN_MERCURY_MOD, ONLY : LVEGEMIS
      USE OCEAN_MERCURY_MOD, ONLY : LBRCHEM,     LRED_JNO2, LGEOSLWC
      USE DEPO_MERCURY_MOD , ONLY : LHGSNOW, ADD_HG2_SNOWPACK
      USE OCEAN_MERCURY_MOD, ONLY : LHg2HalfAerosol
      USE OCEAN_MERCURY_MOD, ONLY : STRAT_BR_FACTOR
      USE OCEAN_MERCURY_MOD, ONLY : LAnthroHgOnly, LOHO3CHEM
      USE OCEAN_MERCURY_MOD, ONLY : LnoUSAemis,  LBROCHEM

      IMPLICIT NONE
#     include "define.h"

      !=================================================================
      ! MODULE PRIVATE DECLARATIONS -- keep certain internal variables 
      ! and routines from being seen outside "mercury_mod.f"
      !=================================================================

      ! Make everything PRIVATE ...
      PRIVATE

      ! ... except these routines
      PUBLIC :: CHEMMERCURY
      PUBLIC :: CLEANUP_MERCURY
      PUBLIC :: INIT_MERCURY
      PUBLIC :: EMISSMERCURY
      PUBLIC :: PARTITIONHG
      PUBLIC :: HGPFRAC

      !=================================================================
      ! MODULE VARIABLES
      !=================================================================

      ! Scalars 
      INTEGER              :: ANTHRO_Hg_YEAR

      ! Parameters
      REAL*8,  PARAMETER   :: SMALLNUM = 1D-20

      ! Arrays
      INTEGER, ALLOCATABLE :: AN_Hg0(:,:)
      INTEGER, ALLOCATABLE :: AN_Hg2(:,:)
      INTEGER, ALLOCATABLE :: AN_HgP(:,:)
      REAL*8,  ALLOCATABLE :: COSZM(:,:)
      REAL*8,  ALLOCATABLE :: EHg0_an(:,:)
      REAL*8,  ALLOCATABLE :: EHg0_am(:,:)
      REAL*8,  ALLOCATABLE :: EHg2_an(:,:)
      REAL*8,  ALLOCATABLE :: EHgP_an(:,:)
      REAL*8,  ALLOCATABLE :: EHg0_oc(:,:,:)
      REAL*8,  ALLOCATABLE :: EHg0_ln(:,:,:)
      REAL*8,  ALLOCATABLE :: EHg0_dist(:,:)
      REAL*8,  ALLOCATABLE :: EHg0_nt(:,:)
      REAL*8,  ALLOCATABLE :: EHg0_bb(:,:)
      REAL*8,  ALLOCATABLE :: EHg0_vg(:,:)
      REAL*8,  ALLOCATABLE :: EHg0_so(:,:)
      REAL*8,  ALLOCATABLE :: EHg0_gtm(:,:)
      REAL*8,  ALLOCATABLE :: EHg0_snow(:,:,:)
      REAL*8,  ALLOCATABLE :: TCOSZ(:,:)
      REAL*8,  ALLOCATABLE :: TTDAY(:,:)
      REAL*8,  ALLOCATABLE :: ZERO_DVEL(:,:)
      REAL*8,  ALLOCATABLE :: HG2_SEASALT_LOSSRATE(:,:) 
      REAL*8,  ALLOCATABLE :: JNO2(:,:,:) 

      ! Henry's Law constant for Hg2 [mol /L /atm]
      REAL*8, PARAMETER     :: HL       = 1.4d6
      REAL*8, ALLOCATABLE :: HGPFRAC(:,:,:)

      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement
      !=================================================================
      CONTAINS

!------------------------------------------------------------------------------

#if defined( DEVEL )
      SUBROUTINE CHEMMERCURY( LOCAL_MET )
#else
      SUBROUTINE CHEMMERCURY
#endif
!
!******************************************************************************
!  Subroutine CHEMMERCURY is the driver routine for mercury chemistry
!  in the GEOS-CHEM module. (eck, bmy, 12/6/04, 4/6/06)
!
!  NOTES:
!  (1 ) Now references routine GET_PBL_MAX_L from "pbl_mix_mod.f".  Now
!        references AD44 from "diag_mod.f".  Now sum the levels from T44 into 
!        the AD44 array.  Now references N_TRACERS from "tracer_mod.f".
!        (bmy, 2/24/05)
!  (2 ) Bug fix: Set T44 to 0e0 for single precision.  Now allow for zero
!        dry deposition velocity.  Now call INIT_MERCURY from "input_mod.f"
!        (bmy, 4/6/06)
!******************************************************************************
!
      ! References to F90 modules
      USE DRYDEP_MOD       , ONLY : DEPSAV
      USE ERROR_MOD        , ONLY : DEBUG_MSG
      USE GLOBAL_O3_MOD    , ONLY : GET_GLOBAL_O3
      USE GLOBAL_OH_MOD    , ONLY : GET_GLOBAL_OH
      USE GLOBAL_BR_MOD    , ONLY : GET_GLOBAL_BR
      USE PBL_MIX_MOD      , ONLY : GET_PBL_MAX_L
      USE LOGICAL_MOD      , ONLY : LPRT, LGTMM, LNLPBL 
      USE TIME_MOD         , ONLY : GET_MONTH, ITS_A_NEW_MONTH
      USE TRACER_MOD       , ONLY : N_TRACERS
      USE TRACERID_MOD     , ONLY : N_HG_CATS
      USE TIME_MOD         , ONLY : ITS_TIME_FOR_A3 
      USE DRYDEP_MOD       , ONLY : DRYHg0, DRYHg2, DRYHgP
      USE OCEAN_MERCURY_MOD, ONLY : Fg, Fp  
      USE OCEAN_MERCURY_MOD, ONLY : PARTITION_Hg2

#if defined( DEVEL )
      USE GC_TYPE_MOD,       ONLY : GC_MET_LOCAL
#endif

      USE CMN_SIZE_MOD        ! Size parameters

      !Arguments
#if defined( DEVEL )
      TYPE(GC_MET_LOCAL), INTENT(IN) :: LOCAL_MET  ! Obj w/ met fields
#endif

      ! Local variables
      LOGICAL, SAVE          :: FIRST = .TRUE.
      INTEGER                :: I, J, L, MONTH, N, PBL_MAX

      REAL*8                 :: K_DRYD2(IIPAR,JJPAR)

      !=================================================================
      ! CHEMMERCURY begins here!
      !
      ! Read monthly mean OH and O3 fields
      !=================================================================
      IF ( ITS_A_NEW_MONTH() ) THEN 

         ! Get the current month
         MONTH = GET_MONTH()

         ! Read monthly mean OH and O3 from disk
         CALL GET_GLOBAL_OH( MONTH )
         IF ( LPRT ) CALL DEBUG_MSG( '### CHEMMERC: a GET_GLOBAL_OH' )

         CALL GET_GLOBAL_O3( MONTH )
         IF ( LPRT ) CALL DEBUG_MSG( '### CHEMMERC: a GET_GLOBAL_O3' )

         ! Renamed to GET_GLOBAL_BR
#if defined( DEVEL )
         CALL GET_GLOBAL_BR( MONTH, LOCAL_MET )
#else
         CALL GET_GLOBAL_BR( MONTH )
#endif
         IF ( LPRT ) CALL DEBUG_MSG( '### CHEMMERC: a GET_GLOBAL_BR' )

         IF (LRED_JNO2) THEN
            CALL GET_GLOBAL_JNO2( MONTH ) 
            IF ( LPRT ) 
     &           CALL DEBUG_MSG( '### CHEMMERC: a GET_GLOBAL_JNO2' )
         ENDIF

      ENDIF
      
      !=================================================================
      ! Perform chemistry on Hg tracers 
      !=================================================================
      
      ! Compute diurnal scaling for OH
      CALL OHNO3TIME
      IF ( LPRT ) CALL DEBUG_MSG( 'CHEMMERCURY: a OHNO3TIME' )

      ! Calculate the rate of sea salt aerosol uptake of Hg2
      IF ( LDYNSEASALT .AND. ITS_TIME_FOR_A3() ) THEN
#if defined( DEVEL )
         CALL CALC_HG2_SEASALT_LOSSRATE( LOCAL_MET )
#else
         CALL CALC_HG2_SEASALT_LOSSRATE
#endif
         IF ( LPRT ) CALL DEBUG_MSG( 'CHEMMERCURY: a SEASALT_LOSSRATE' )
      ENDIF

      ! Choose dry deposition frequency for Hg2, 1/s
      IF (LHg2HalfAerosol) THEN
         K_DRYD2 = ( DEPSAV(:,:,DRYHg2) +  DEPSAV(:,:,DRYHgP) ) / 2D0
      ELSE
         !--------------------------------------------------------
         ! Previous to 25 Oct 2011, H Amos
         !K_DRYD2 = DEPSAV(:,:,DRYHg2)
         !
         ! Weigh dry dep velocity by gas and aerosol fractions. 
         ! [ref: Amos et al., 2011, submitted to ACP]
         K_DRYD2 = ( DEPSAV(:,:,DRYHg2)*Fg(:,:,1) + 
     &               DEPSAV(:,:,DRYHgP)*Fp(:,:,1) )
         !--------------------------------------------------------
      ENDIF


      !-------------------------
      ! Hg0 and Hg2 chemistry
      !-------------------------
      IF ( LPRT ) CALL DEBUG_MSG( 'CHEMMERCURY: b CHEM_Hg0_Hg2' )
      
#if defined( DEVEL )
      ! Add option for non-local PBL (cdh, 08/27/09)
      IF ( LNLPBL ) THEN

         ! Dry deposition occurs with PBL mixing,
         ! pass zero deposition frequency
         CALL CHEM_Hg0_Hg2( ZERO_DVEL, ZERO_DVEL, LOCAL_MET )

      ELSE

         IF ( DRYHg2 > 0 .and. DRYHg0 > 0 ) THEN

            ! Dry deposition active for both Hg0 and Hg2; 
            ! pass drydep frequency to CHEM_Hg0_Hg2
            ! (NOTE: DEPSAV has units 1/s)
            CALL CHEM_Hg0_Hg2( K_DRYD2, DEPSAV(:,:,DRYHg0), LOCAL_MET )

         ELSEIF (DRYHg2 > 0 .and. DRYHg0 .le. 0) THEN

            ! Only Hg2 dry deposition is active
            CALL CHEM_Hg0_Hg2( K_DRYD2, ZERO_DVEL, LOCAL_MET )

         ELSEIF (DRYHg2 <= 0 .and. DRYHg0 > 0) THEN

            ! Only Hg0 dry deposition is active
            CALL CHEM_Hg0_Hg2( ZERO_DVEL , DEPSAV(:,:,DRYHg0),
     &                         LOCAL_MET )

         ELSE

            ! No dry deposition, pass zero deposition frequency
            CALL CHEM_Hg0_Hg2( ZERO_DVEL , ZERO_DVEL, LOCAL_MET )

#else

      ! Add option for non-local PBL (cdh, 08/27/09)
      IF ( LNLPBL ) THEN

         ! Dry deposition occurs with PBL mixing,
         ! pass zero deposition frequency
         CALL CHEM_Hg0_Hg2( ZERO_DVEL, ZERO_DVEL)

      ELSE

         IF ( DRYHg2 > 0 .and. DRYHg0 > 0 ) THEN
         
            ! Dry deposition active for both Hg0 and Hg2; 
            ! pass drydep frequency to CHEM_Hg0_Hg2 (NOTE: DEPSAV has units 1/s)
            CALL CHEM_Hg0_Hg2( K_DRYD2, DEPSAV(:,:,DRYHg0) )
            
         ELSEIF (DRYHg2 > 0 .and. DRYHg0 .le. 0) THEN

            ! Only Hg2 dry deposition is active
            CALL CHEM_Hg0_Hg2( K_DRYD2, ZERO_DVEL) 
            
         ELSEIF (DRYHg2 <= 0 .and. DRYHg0 > 0) THEN

            ! Only Hg0 dry deposition is active
            CALL CHEM_Hg0_Hg2( ZERO_DVEL , DEPSAV(:,:,DRYHg0))
            
         ELSE

            ! No dry deposition, pass zero deposition frequency
            CALL CHEM_Hg0_Hg2( ZERO_DVEL , ZERO_DVEL)
#endif

         ENDIF

      ENDIF      

      IF ( LPRT ) CALL DEBUG_MSG( 'CHEMMERCURY: a CHEM_Hg0_Hg2' )
   
      !--------------------------
      ! HgP chemistry
      !--------------------------
      IF ( LPRT ) CALL DEBUG_MSG( 'CHEMMERCURY: b CHEM_HgP' )
      
      ! Add option for non-local PBL (cdh, 08/27/09)
      ! Dry deposition done with PBL mixing if non-local selected
      IF ( DRYHgP == 0 .OR. LNLPBL ) THEN

         ! Otherwise pass zero drydep frequency
         CALL CHEM_HgP( ZERO_DVEL )

      ELSE
         
         ! If DRYHgP > 0, then drydep is active;
         ! Pass drydep frequency to CHEM_HgP
         CALL CHEM_HgP( DEPSAV(:,:,DRYHgP) )

      ENDIF

      IF ( LPRT ) CALL DEBUG_MSG( 'CHEMMERCURY: a CHEM_HgP' )

      ! Return to calling program
      END SUBROUTINE CHEMMERCURY

!------------------------------------------------------------------------------

#if defined( DEVEL )
      SUBROUTINE CHEM_Hg0_Hg2( V_DEP_Hg2, V_DEP_Hg0, LOCAL_MET )
#else
      SUBROUTINE CHEM_Hg0_Hg2( V_DEP_Hg2, V_DEP_Hg0 )
#endif
!
!******************************************************************************
!  Subroutine CHEM_Hg0_Hg2 is the chemistry subroutine for the oxidation,
!  reduction and deposition of Hg(0) and Hg(II), including tagged tracers of 
!  these species.
!  (eck, bmy, cdh, 12/6/04, 1/9/06, 7/25/08)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) V_DEP_Hg2 (REAL*8) : Dry deposition frequency for Hg(II) [/s]
!  (1 ) V_DEP_Hg0 (REAL*8) : Dry deposition frequency for Hg(0) [/s]
!
!  Description of the chemistry mechanism:
!  ============================================================================
!  The oxidation, reduction and deposition properties can be changed with 
!  logical switches found in INIT_MERCURY (below). Here we describe two 
!  alternatives.
!
!  -----------------------------------------------
!  Chemistry and deposition as described by Holmes et al. (2010):
!
!  (1 ) Oxidation: Hg(0) --> Hg(II):
!       
!     Hg(0)(g) + Br(g) --> + Br/OH --> Hg(II), rates are selected with
!        METHOD keyword below. Recommded kinetics are 'DonohoueYBBalabanov'
!        which use rates from Donohoue et al. (2006), Goodsite et al. (2004)
!        and Balabanov et al. (2005)
!           
!  (2 ) Aqueous-phase photochemical reduction of Hg(II) is included based
!        on estimate of rate constant and scaled to NO2 photolysis. The
!        rate is tuned to match the global Hg(0) concentration and
!        seasonal cycle.
! 
!  (3 ) Hg(II) is dry deposited,        kd calculated by drydep_mod [/s]     
!        The dry deposition rate is an average of the particulate and gaseous
!        deposition rates, reflecting that Hg(II) partitions between gas and
!        aerosol.
!
!  (4 ) Hg(0) is dry deposited,         kd calculated by drydep_mod [/s]
!         The ocean module separately cacluates Hg(0) dry deposition over
!         ocean, so this module only includes Hg(0) dry deposition over land.
!
!  -----------------------------------------------
!  Chemistry and deposition as described by Selin et al. (2010):
!
!  (1 ) Oxidation: Hg(0) --> Hg(II): 
!
!         Hg(0)(g)+ O3(g) --> Hg(II) ,  k  = 3.0e-20 cm3 molec-1 s-1
!                                       Source: Hall, 1995
!       
!         Hg(0)(g)+ OH(g) --> Hg(II) ,  k  = 8.7e-14 cm3 s-1 
!                                       Source: Sommar et al. 2001
!           
!  (2 ) Aqueous-phase photochemical reduction of Hg(II) is included based
!        on estimate of rate constant and scaled to OH concentrations. The
!        rate is tuned with the OH oxidation rate to match the global Hg(0) 
!        concentration and seasonal cycle.
! 
!  (3 ) Hg(II) is dry-deposited,        kd calculated by drydep_mod [/s]    
!
!  (4 ) Hg(0) is dry deposited,         kd calculated by drydep_mod [/s]
!         The ocean module separately cacluates Hg(0) dry deposition over
!         ocean, so this module only includes Hg(0) dry deposition over land.
!     
!  NOTES:
!  (1 ) Updated for reduction reaction.  Now use diagnostic arrays from
!        "diag03_mod.f" (eck, bmy, 1/21/05)
!  (2 ) Now references GET_FRAC_UNDER_PBLTOP from "pbl_mix_mod.f".  Now
!        performs drydep for all levels in the PBL.  Changed Kred to 2.1e-10
!        on advice of Noelle Eckley Selin. (bmy, 2/24/05)
!  (3 ) Now references XNUMOL from "tracer_mod.f" (bmy, 10/25/05)
!  (4 ) Now prevent divide-by-zero error.  Now use ID_Hg0 and ID_Hg2 index
!        arrays from "tracerid_mod.f".  Also modified for updated ocean
!        mercury module.  Updated some constants.  Also saves out diagnostic 
!        of Hg2 lost to rxn w/ seasalt. (eck, cdh, sas, bmy, 4/6/06)
!  (5 ) Added Hg0 dry deposition (eck)
!  (6 ) Chemistry and dry deposition now occur simultaneously. (cdh, 7/9/08)
!  13 Aug 2010 - R. Yantosca - Treat MERRA like GEOS-5 
!  14 Jan 2011 - R. Yantosca - Now get volume cloud fraction directly
!                              from GEOS-5 and MERRA met fields
!  08 Feb 2012 - R. Yantosca - Treat GEOS-5.7.x in the same way as MERRA
!  01 Mar 2012 - R. Yantosca - Now use GET_AREA_CM2(I,J,L) from grid_mod.F90
!******************************************************************************
!
      ! References to F90 modules
      USE DAO_MOD,          ONLY : T, AD, IS_WATER, IS_ICE, IS_LAND
      USE DIAG03_MOD,       ONLY : AD03_Hg2_Hg0, AD03_Hg2_O3
      USE DIAG03_MOD,       ONLY : AD03_Hg2_Br, AD03_Hg2_OH
      USE DIAG03_MOD,       ONLY : AD03_Hg2_SS,  LD03,        ND03
      USE DIAG03_MOD,       ONLY : AD03_Hg2_SSR 
      USE DIAG03_MOD,       ONLY : AD03_Br 
      USE LOGICAL_MOD,      ONLY : LSPLIT, LGTMM !ccc for GTMM
      USE PBL_MIX_MOD,      ONLY : GET_FRAC_UNDER_PBLTOP
      USE TIME_MOD,         ONLY : GET_TS_CHEM
      USE TRACER_MOD,       ONLY : STT,    XNUMOL
      USE TRACERID_MOD,     ONLY : ID_Hg0, ID_Hg2, ID_Hg_tot, N_Hg_CATS
      USE TRACERID_MOD,     ONLY : IS_Hg2
      USE GRID_MOD,         ONLY : GET_AREA_CM2
      USE LOGICAL_MOD,      ONLY : LPRT, LDYNOCEAN, LNLPBL
      USE ERROR_MOD,        ONLY : DEBUG_MSG 
      USE DEPO_MERCURY_MOD, ONLY : ADD_HG2_DD
      USE DAO_MOD,          ONLY : SNOW, SNOMAS
      USE PRESSURE_MOD,     ONLY : GET_PCENTER 
      USE DIAG_MOD,         ONLY : AD44
      USE DAO_MOD,          ONLY : AIRDEN, QL
      USE DAO_MOD,          ONLY : AD, CLDF
      USE DAO_MOD,          ONLY : LWI
      USE DAO_MOD,          ONLY : FRSNO, FRLANDIC, FROCEAN

#if defined( DEVEL )
      USE GC_TYPE_MOD,      ONLY : GC_MET_LOCAL
#endif

      USE CMN_SIZE_MOD       ! Size parameters
      USE CMN_DIAG_MOD       ! ND44

      ! Arguments
      REAL*8, INTENT(IN)    :: V_DEP_Hg2(IIPAR,JJPAR)
      REAL*8, INTENT(IN)    :: V_DEP_Hg0(IIPAR,JJPAR)
#if defined( DEVEL )
      TYPE(GC_MET_LOCAL), INTENT(IN) :: LOCAL_MET  ! Obj w/ met fields
#endif

      ! Local variables
      INTEGER               :: I, J, L, NN
      REAL*8                :: DTCHEM
      REAL*8                :: FC,          FA,          F_PBL
      REAL*8                :: LWC,         AREA_CM2
      REAL*8                :: C_O3,        C_OH,        C_BR
      REAL*8                :: C_BRO
      REAL*8                :: K_O3,        K_OH,        K_BR, K_BRO   
      REAL*8                :: K_OX,        K_RED       
      REAL*8                :: E_KOX_T,     E_KRED_T    

      REAL*8                :: K_DRYD0,     K_DRYD2,     K_SALT
      REAL*8                :: K_DEP0,      K_DEP2

      REAL*8                :: OLD_HG0,     OLD_HG2
      REAL*8                :: NEW_HG0,     NEW_HG2
      REAL*8                :: HG0_BL,      HG2_BL
      REAL*8                :: HG0_FT,      HG2_FT

      REAL*8                :: TMP_HG0,     TMP_HG2,     TMP_OX
      REAL*8                :: GROSS_OX,    NET_OX 
      REAL*8                :: GROSS_OX_OH, GROSS_OX_O3, GROSS_OX_BR
      REAL*8                :: DEP_HG0,     DEP_HG2
      REAL*8                :: DEP_HG2_DRY, DEP_HG2_SALT
      REAL*8                :: DEP_DRY_FLX
      REAL*8                :: FRAC_NO_HG0_DEP !jaf 
      LOGICAL               :: ZERO_HG0_DEP !jaf 

      ! K for reaction Hg0 + O3  [cm3 /molecule /s] (Source: Hall 1995)
      REAL*8, PARAMETER     :: K_HG_O3  = 3.0d-20 !3.0d-20 (Gas phase)

      ! K for reaction Hg2 + OH  [cm3 /molecule /s] (Source: Sommar 2001)
      REAL*8, PARAMETER     :: K_HG_OH  = 8.7d-14 !8.7d-14 (Gas phase)

      ! K for reaction Hg2 + BrO  [cm3 /molecule /s]
      ! (Source: Raofie and Ariya 2003; 2004)
      REAL*8, PARAMETER     :: K_HG_BRO  = 1d-14 !1d-15 - 1d-14

      ! Gas constant [L atm /K /mol]
      REAL*8, PARAMETER     :: R        = 8.2d-2

      ! K for reduction [cm3 /molecule /s] (Source: Selin 2007)
      ! This variable is unused if LREDJNO2 is TRUE
      REAL*8, PARAMETER     :: K_RED_OH = 1d-10!4.2d-10 from Noelle
                                !4d-10 works well for Hg+OH/O3
     
      ! K for reduction, using J_NO2, [unitless scale factor]
      ! Source: Holmes et al. 2010
      ! 3.5D-3 for Hg+Br simulation; 1.3D-2 for Hg+OH/O3 simulation

      ! Reduce K_RED_JNO2 for GEIA 2005 inventory
      ! recommended by Helen Amos, 23 Sep 2011
      ! change implemented 10/19/11 by eck
      ! inserted ifdefs for future tuning of 2x2.5 by eck
#if defined (GRID4x5)
      REAL*8, PARAMETER     :: K_RED_JNO2 = 3.5D-3/2d0 
#elif defined (GRID2x25)
      REAL*8, PARAMETER     :: K_RED_JNO2 = 3.5D-3/2d0
#else
      REAL*8, PARAMETER     :: K_RED_JNO2 = 3.5D-3/2d0
#endif      

      ! K for sea salt (eck, bmy, 4/6/06) [/s]
      REAL*8, PARAMETER     :: K_SALT_FIXED   = 3.8d-5

      ! Set of Hg/Br rate constants to use
      ! Recommended: DonohoueYBBalabanov, GoodsiteY, or DonohoueYB
      ! DonohoueYBBalabanov used by Holmes et al. 2010
      CHARACTER(LEN=*), PARAMETER  :: METHOD='DonohoueYBBalabanov' 


      !=================================================================
      ! CHEM_Hg0_Hg2 begins here!
      !=================================================================

      ! Chemistry timestep [s]
      DTCHEM = GET_TS_CHEM() * 60d0

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I,       J,        L,         NN                   )
!$OMP+PRIVATE( F_PBL,   FC,       FA,        LWC,      AREA_CM2   )
!$OMP+PRIVATE( C_O3,    C_OH,     C_BR,      C_BRO                )   
!$OMP+PRIVATE( K_O3,    K_OH,    K_BR,       K_BRO                )
!$OMP+PRIVATE( K_OX,    K_RED,   E_KOX_T,    E_KRED_T             )
!$OMP+PRIVATE( K_DEP0,  K_DEP2,  K_DRYD0,    K_DRYD2,  K_SALT     )
!$OMP+PRIVATE( OLD_HG0, OLD_HG2,  NEW_HG0,   NEW_HG2              )
!$OMP+PRIVATE( HG0_BL,  HG2_BL,   HG0_FT,    HG2_FT               )
!$OMP+PRIVATE( TMP_HG0, TMP_HG2,  TMP_OX                          )    
!$OMP+PRIVATE( DEP_HG0, DEP_HG2,  DEP_HG2_DRY, DEP_HG2_SALT       )
!!$OMP+PRIVATE( DEP_DRY_FLX,      SNOW_HT                          )
!$OMP+PRIVATE( DEP_DRY_FLX,       ZERO_HG0_DEP, FRAC_NO_HG0_DEP   )
!$OMP+PRIVATE( NET_OX,  GROSS_OX                                  )
!$OMP+PRIVATE( GROSS_OX_OH,       GROSS_OX_O3, GROSS_OX_BR        )
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR


         ! Fraction of box (I,J,L) underneath the PBL top [dimensionless]
         F_PBL = GET_FRAC_UNDER_PBLTOP( I, J, L )

         ! Monthly mean O3, OH Br concentrations [molec/cm3]
#if defined( DEVEL )
         C_BR        = GET_BR( I, J, L, LOCAL_MET, C_BRO )
         C_O3        = GET_O3( I, J, L, LOCAL_MET )
         C_OH        = GET_OH( I, J, L )
#else
         C_BR        = GET_BR( I, J, L, C_BRO )
         C_O3        = GET_O3( I, J, L )
         C_OH        = GET_OH( I, J, L )
#endif

#if   defined( GEOS_5 ) || defined( MERRA ) || defined( GEOS_57 )

         !---------------------------------------------
         ! GEOS-5/MERRA: Get LWC, FC from met fields
         ! (skim, bmy, 1/14/10)
         !---------------------------------------------

#if defined( DEVEL )
         ! Get cloud fraction from met fields [unitless]
         FC = LOCAL_MET%CLDF(L,I,J)

         IF ( LGEOSLWC ) THEN

            ! Get liquid water content [m3 H2O/m3 air] within cloud from met
            ! Units: [kg H2O/kg air] * [kg air/m3 air] * [m3 H2O/1e3 kg H2O]
            LWC = LOCAL_MET%QL(I,J,L) * LOCAL_MET%AIRDEN(L,I,J) * 1D-3

            !%%% NOTE: Someone should investigate effect of dividing LWC by the
            !%%% cloud fraction, as is done in sulfate_mod.f (bmy, 5/27/11)

         ELSE

            ! Get liquid water content for entire grid box,
            ! based on formula for LWC in cloud [m3 H2O/m3 air]
            LWC = GET_LWC( LOCAL_MET%T(I,J,L) ) * FC 

            ! There should be no liquid water when T < 258K 
            IF ( LOCAL_MET%T(I,J,L) < 258D0 )  LWC = 0D0

         ENDIF
#else
         ! Get cloud fraction from met fields [unitless]
         FC = CLDF(L,I,J)

         IF ( LGEOSLWC ) THEN

            ! Get liquid water content [m3 H2O/m3 air] within cloud from met
            ! Units: [kg H2O/kg air] * [kg air/m3 air] * [m3 H2O/1e3 kg H2O]
            LWC = QL(I,J,L) * AIRDEN(L,I,J) * 1D-3

            !%%% NOTE: Someone should investigate effect of dividing LWC by the
            !%%% cloud fraction, as is done in sulfate_mod.f (bmy, 5/27/11)

         ELSE

            ! Get liquid water content for entire grid box,
            ! based on formula for LWC in cloud [m3 H2O/m3 air]
            LWC = GET_LWC( T(I,J,L) ) * FC 

            ! There should be no liquid water when T < 258K 
            IF ( T(I,J,L) < 258D0 )  LWC = 0D0

         ENDIF
#endif

#else
         !---------------------------------------------
         ! Otherwise, compute FC, LWC as before
         !---------------------------------------------    

         ! Get volume fraction of gridbox containing 
         ! cloud [unitless], following Sundqvist et al 1989.
#if defined( DEVEL )
         FC          = GET_VCLDF( I, J, L, LOCAL_MET )

         ! Get liquid water content for entire grid box,
         ! based on formula for LWC in cloud (m3/m3)
         LWC         = GET_LWC( LOCAL_MET%T(I,J,L) ) * FC 

         ! There should be no liquid water when T < 258K 
         IF ( LOCAL_MET%T(I,J,L) < 258D0 )  LWC = 0D0
#else
         FC          = GET_VCLDF( I, J, L )

         ! Get liquid water content for entire grid box,
         ! based on formula for LWC in cloud (m3/m3)
         LWC         = GET_LWC( T(I,J,L) ) * FC 

         ! There should be no liquid water when T < 258K 
         IF ( T(I,J,L) < 258D0 )  LWC = 0D0
#endif

#endif

#if defined( DEVEL )
         ! Define fraction of Hg(II) which is in aqueous solution
         ! [dimensionless]
         FA          = ( HL * R * LOCAL_MET%T(I,J,L) * LWC )
         FA          = FA / ( 1d0 + FA )

         ! Cl- in sea-salt aerosol enhances solubility 2000X in MBL
         IF ( LRED_JNO2 .AND. ( F_PBL >0.1  )  .AND. 
     &        IS_WATER( I, J, LOCAL_MET )) THEN
            FA          = ( HL * 2D3 * R * LOCAL_MET%T(I,J,L) * LWC )
            FA          = FA / ( 1d0 + FA )
            
         ENDIF

         ! Define K's for the oxidation reactions [/s]
         K_O3        = K_HG_O3 * C_O3
         K_OH        = K_HG_OH * C_OH

         IF (LBRCHEM) THEN
            K_BR = GET_HGBR_RATE( C_BR, LOCAL_MET%T(I,J,L), 
     &             GET_PCENTER(I,J,L), C_OH, METHOD )
         ELSE
            K_BR = 0d0
         ENDIF
#else
         ! Define fraction of Hg(II) which is in aqueous solution
         ! [dimensionless]
         FA          = ( HL * R * T(I,J,L) * LWC )
         FA          = FA / ( 1d0 + FA )

         ! Cl- in sea-salt aerosol enhances solubility 2000X in MBL
         IF (LRED_JNO2 .AND. (F_PBL >0.1) .AND. IS_WATER(I,J)) THEN
            FA          = ( HL * 2D3 * R * T(I,J,L) * LWC )
            FA          = FA / ( 1d0 + FA )
            
         ENDIF

         ! Define K's for the oxidation reactions [/s]
         K_O3        = K_HG_O3 * C_O3
         K_OH        = K_HG_OH * C_OH

         IF (LBRCHEM) THEN
            K_BR = GET_HGBR_RATE( C_BR, T(I,J,L), 
     &           GET_PCENTER(I,J,L), C_OH, METHOD )
         ELSE
            K_BR = 0d0
         ENDIF
#endif

         IF (LBROCHEM) THEN
            K_BRO = K_HG_BRO * C_BRO
         ELSE
            K_BRO = 0d0
         ENDIF
         

         IF (.NOT. LOHO3CHEM) THEN
            K_O3 = 0D0
            K_OH = 0D0
         ENDIF

         K_OX        = K_O3 + K_OH + K_BR + K_BRO
         
         ! Define K for the reduction reaction. 
         IF (LRED_JNO2) THEN
            K_RED = K_RED_JNO2 * FA * GET_JNO2(I,J,L)

         ELSE
            ! Include the fraction of 
            ! Hg(II) in air within the Kreduction & 
            ! scale to OH concentration [/s]
            K_RED = K_RED_OH * FA * C_OH
         ENDIF

         ! Round very small K_OX, K_RED to prevent numerical errors
         ! (jaf, cdh, 11/17/11)
         IF ( K_OX  < 1D-10 ) K_OX  = 0d0
         IF ( K_RED < 1D-10 ) K_RED = 0d0

         ! Define K for dry deposition, [/s]
         K_DRYD0 = V_DEP_HG0(I,J)
         K_DRYD2 = V_DEP_HG2(I,J)

         ! Define K for total deposition over ocean [/s]
         ! IS_WATER returns true for ocean boxes. Fresh water may be TRUE in 
         ! principle but is not at 4x5 resolution.
         ! Sea salt will not be active in coastal boxes that have some ocean
         ! but where IS_WATER is FALSE.
         ! Note: we only need the IF structure for the fixed
         ! K_SALT, so that deposition doesn't occur over land.
         ! Remove the IF statement once HG2_SEASALT_LOSSRATE 
         ! is in the standard version. (cdh, 2/8/08)
#if defined( DEVEL )
         IF ( IS_WATER( I, J, LOCAL_MET ) ) THEN
#else
         IF ( IS_WATER(I,J) ) THEN
#endif

            IF (LDYNSEASALT) THEN
               ! Uptake based on sea-salt aerosol flux and salinity [/s]
               K_SALT = HG2_SEASALT_LOSSRATE(I,J)
            ELSE
               !Constant rate tuned for Okinawa [/s]
               K_SALT = K_SALT_FIXED
            ENDIF

            K_DEP2 = K_DRYD2 + K_SALT

         ELSE
            K_SALT = 0d0

            K_DEP2 = K_DRYD2 
         ENDIF
         
!--jaf.start
         ! Restructure to allow use of fractional land area info in MERRA
         ! jaf, 4/26/11
!         ! Disable dry deposition of Hg(0) to ice because we do not have
!         ! an ice emission model. Perennial ice should have equal emission
!         ! and deposition averaged over multiple years. (cdh, 9/11/09)
!#if   defined( GEOS_5 ) || defined( MERRA )
!         ! GEOS5 snow height (water equivalent) in mm. (Docs wrongly say m)
!         SNOW_HT = SNOMAS(I,J)
!#else
!         ! GEOS1-4 snow heigt (water equivalent) in mm
!         SNOW_HT = SNOW(I,J)
!#endif 
!         IF ( IS_ICE(I,J) .OR. (IS_LAND(I,J) .AND. SNOW_HT>10d0) ) THEN
!            K_DEP0 = 0D0
!         ENDIF

         ! Hg(0) exchange with the ocean is handled by ocean_mercury_mod
         ! so disable deposition over water here.
         ! Turn off Hg(0) deposition to snow and ice because we haven't yet
         ! included emission from these surfaces and most field studies
         ! suggest Hg(0) emissions exceed deposition during sunlit hours.
            
         ! Except in MERRA, we assume entire grid box is water or ice
         ! if conditions are met (jaf, 4/26/11)
         FRAC_NO_HG0_DEP = 1d0

#if defined( DEVEL )

#if   defined( MERRA ) || defined( GEOS_57 )
         FRAC_NO_HG0_DEP = MIN( LOCAL_MET%FROCEAN(I,J) + 
     &                          LOCAL_MET%FRSNO(I,J)   +
     &                          LOCAL_MET%FRLANDIC(I,J), 1d0 ) 
         ZERO_HG0_DEP = ( FRAC_NO_HG0_DEP > 0d0 )

#elif defined( GEOS_5 )
         ! GEOS5 snow height (water equivalent) in mm. (Docs wrongly say m)
         ZERO_HG0_DEP = (( LOCAL_MET%LWI(I,J)   == 0    )  .OR. 
     &                   ( IS_ICE ( I, J, LOCAL_MET     )) .OR.
     &                   ( IS_LAND( I, J, LOCAL_MET     )  .AND.
     &                     LOCAL_MET%SNOMAS(I,J) > 10d0 ))

#else
         ! GEOS1-4 snow heigt (water equivalent) in mm
         ZERO_HG0_DEP = (( LOCAL_MET%LWI(I,J) == 0    )  .OR.
     &                   ( IS_ICE ( I, J, LOCAL_MET   )) .OR. 
     &                   ( IS_LAND( I, J, LOCAL_MET   )  .AND.
     &                     LOCAL_MET%SNOW(I,J) > 10d0 ))
#endif 

#else

#if   defined( MERRA ) || defined( GEOS_57 )
         FRAC_NO_HG0_DEP =
     &        MIN(FROCEAN(I,J) + FRSNO(I,J) + FRLANDIC(I,J), 1d0) 
         ZERO_HG0_DEP = ( FRAC_NO_HG0_DEP > 0d0 )
#elif defined( GEOS_5 )
         ! GEOS5 snow height (water equivalent) in mm. (Docs wrongly say m)
         ZERO_HG0_DEP = ( (LWI(I,J) == 0) .OR. 
     &                   (IS_ICE(I,J)) .OR.
     &                   (IS_LAND(I,J) .AND. SNOMAS(I,J) > 10d0) )
#else
         ! GEOS1-4 snow heigt (water equivalent) in mm
         ZERO_HG0_DEP = ( (LWI(I,J) == 0) .OR.
     &                   (IS_ICE(I,J)) .OR. 
     &                   (IS_LAND(I,J) .AND. SNOW(I,J) > 10d0) )
#endif

#endif  ! DEVEL

         K_DEP0 = K_DRYD0
         IF ( ZERO_HG0_DEP ) THEN
            K_DEP0 = K_DEP0 * MAX(1d0 - FRAC_NO_HG0_DEP,0d0)
         ENDIF
!--jaf.end

         ! Precompute exponential factors [dimensionless]
         E_KOX_T   = EXP( -K_OX   * DTCHEM )
         E_KRED_T  = EXP( -K_RED  * DTCHEM )

         !==============================================================
         ! CHEMISTRY AND DEPOSITION REACTIONS
         !==============================================================
         
         ! Loop over the Hg regional tags
         DO NN = 1, N_HG_CATS

            ! Initial concentrations of Hg(0) and Hg(II) [kg]
            OLD_Hg0 = MAX( STT(I,J,L,ID_Hg0(NN)), SMALLNUM )
            OLD_Hg2 = MAX( STT(I,J,L,ID_Hg2(NN)), SMALLNUM )   

            IF ( F_PBL < 0.05D0 .OR. 
     &           (K_DEP0 < SMALLNUM .AND. K_DEP2 < SMALLNUM) ) THEN

               !==============================================================
               ! Entire box is in the free troposphere
               ! or deposition is turned off, so use RXN without deposition
               !==============================================================

               CALL RXN_REDOX_NODEP(  OLD_HG0, OLD_HG2, 
     &              K_OX,    K_RED,   DTCHEM, 
     &              E_KOX_T, E_KRED_T,
     &              NEW_HG0, NEW_HG2, GROSS_OX )

               ! No deposition occurs [kg]
               DEP_HG0 = 0D0
               DEP_HG2 = 0D0

            ELSE IF ( F_PBL > 0.95D0 ) THEN 

               !==============================================================
               ! Entire box is in the boundary layer
               ! so use RXN with deposition
               !==============================================================

               CALL RXN_REDOX_WITHDEP( OLD_HG0,   OLD_HG2,
     &              K_OX,    K_RED,    K_DEP0,    K_DEP2,   DTCHEM, 
     &              E_KOX_T, E_KRED_T,
     &              NEW_HG0, NEW_HG2,  GROSS_OX,  DEP_HG0,  DEP_HG2 )

            ELSE

               !==============================================================
               ! Box spans the top of the boundary layer
               ! Part of the mass is in the boundary layer and subject to 
               ! deposition while part is in the free troposphere and
               ! experiences no deposition.
               !
               ! We apportion the mass between the BL and FT according to the
               ! volume fraction of the box in the boundary layer.
               ! Arguably we should assume uniform mixing ratio, instead of
               ! uniform density but if the boxes are short, the air density
               ! doesn't change much.
               ! But assuming uniform mixing ratio across the inversion layer
               ! is a poor assumption anyway, so we are just using the
               ! simplest approach.
               !==============================================================

               ! Boundary layer portion of Hg [kg]
               Hg0_BL = OLD_HG0 * F_PBL
               Hg2_BL = OLD_HG2 * F_PBL

               ! Free troposphere portion of Hg [kg]
               Hg0_FT = OLD_HG0 - Hg0_BL
               Hg2_FT = OLD_HG2 - Hg2_BL
               
               ! Do chemistry with deposition on BL fraction
               CALL RXN_REDOX_WITHDEP( Hg0_BL,  Hg2_BL,
     &              K_OX,    K_RED,    K_DEP0,    K_DEP2,   DTCHEM, 
     &              E_KOX_T, E_KRED_T, 
     &              NEW_HG0, NEW_HG2,  GROSS_OX,  DEP_HG0,  DEP_HG2 )

               ! Do chemistry without deposition on the FT fraction
               CALL RXN_REDOX_NODEP(  Hg0_FT, Hg2_FT,
     &              K_OX,    K_RED,   DTCHEM, 
     &              E_KOX_T, E_KRED_T,
     &              TMP_HG0, TMP_HG2, TMP_OX )
               
               ! Recombine the boundary layer and free troposphere parts [kg]
               NEW_HG0 = NEW_HG0 + TMP_HG0
               NEW_HG2 = NEW_HG2 + TMP_HG2
               
               ! Total gross oxidation in the BL and FT [kg]
               GROSS_OX = GROSS_OX + TMP_OX

            ENDIF

            ! Ensure positive concentration [kg]
            NEW_HG0 = MAX( NEW_HG0, SMALLNUM )
            NEW_HG2 = MAX( NEW_HG2, SMALLNUM )

            ! Archive new Hg values [kg]
            STT(I,J,L,ID_Hg0(NN)) = NEW_Hg0 
            STT(I,J,L,ID_Hg2(NN)) = NEW_Hg2

            ! Net oxidation [kg]
            NET_OX = OLD_HG0 - NEW_HG0 - DEP_HG0

            ! Error check on gross oxidation [kg]
            IF ( GROSS_OX < 0D0 ) 
     &          CALL DEBUG_MSG('CHEM_HG0_HG2: negative gross oxidation')

            ! Apportion gross oxidation between O3 and OH [kg]
            IF ( (K_OX     < SMALLNUM) .OR. 
     &           (GROSS_OX < SMALLNUM) ) THEN
               GROSS_OX_OH = 0D0
               GROSS_OX_BR = 0D0
               GROSS_OX_O3 = 0D0
            ELSE
               GROSS_OX_OH = GROSS_OX * K_OH / K_OX
               GROSS_OX_BR = GROSS_OX * K_BR / K_OX
               GROSS_OX_O3 = GROSS_OX * K_O3 / K_OX
            ENDIF
               
            ! Apportion deposition between dry deposition and sea salt [kg]
            IF ( (K_DEP2  < SMALLNUM) .OR. 
     &           (DEP_HG2 < SMALLNUM) ) THEN
               DEP_HG2_SALT = 0D0
               DEP_HG2_DRY  = 0D0
            ELSE
               DEP_HG2_DRY  = DEP_HG2 * K_DRYD2 / K_DEP2
               DEP_HG2_SALT = DEP_HG2 - DEP_HG2_DRY 
            ENDIF
            
            ! Add deposited Hg(II) to the ocean module. OCEAN_MERCURY_MOD
            ! determines whether the box is marine, so we don't need to here.
            ! We should add an if statement to test whether DYNAMIC LAND is 
            ! active.
            IF ( LDYNOCEAN .AND. (DEP_Hg2 > 0d0) ) THEN
               CALL ADD_Hg2_DD( I, J, ID_Hg2(NN), DEP_HG2 )
            ENDIF           

            ! Add deposited Hg(II) to the snowpack
            IF ( LHGSNOW .AND. (DEP_Hg2>0d0) ) THEN
#if defined( DEVEL )
               CALL ADD_HG2_SNOWPACK(I,J,ID_Hg2(NN),DEP_HG2_DRY,
     &                               LOCAL_MET)
#else
               CALL ADD_HG2_SNOWPACK(I,J,ID_Hg2(NN),DEP_HG2_DRY)
#endif
            ENDIF
               

            !=================================================================
            ! ND44 diagnostic: drydep flux of Hg(II) [molec/cm2/s]
            !=================================================================
            IF ( ( ND44 > 0 .OR. LGTMM ) .AND. (.NOT. LNLPBL) ) THEN

               ! Grid box surface area [cm2]
               AREA_CM2 = GET_AREA_CM2( I, J, L )

               ! Amt of Hg(II) lost to drydep [molec/cm2/s]
               DEP_DRY_FLX  = DEP_HG2_DRY * XNUMOL(ID_Hg2(NN)) / 
     &              ( AREA_CM2 * DTCHEM )

               ! Archive Hg(II) drydep flux in AD44 array [molec/cm2/s]
               AD44(I,J,ID_HG2(NN),1) = AD44(I,J,ID_HG2(NN),1) +
     &              DEP_DRY_FLX

               ! Amt of Hg(0) lost to drydep [molec/cm2/s]
               DEP_DRY_FLX  = DEP_HG0 * XNUMOL(ID_Hg0(NN)) / 
     &              ( AREA_CM2 * DTCHEM )

               ! Archive Hg(0) drydep flux in AD44 array [molec/cm2/s]
               AD44(I,J,ID_HG0(NN),1) = AD44(I,J,ID_HG0(NN),1) +
     &              DEP_DRY_FLX

            ENDIF

            !==============================================================
            ! ND03 diagnostic: Hg(II) production [kg]
            !==============================================================
            IF ( ND03 > 0 .AND. L <= LD03 ) THEN 

               ! Store chemistry diagnostics only for total tracer
               IF ( ID_HG0(NN) == ID_HG_TOT) THEN

                  AD03_Hg2_Hg0(I,J,L)= AD03_Hg2_Hg0(I,J,L) + NET_OX
                  AD03_Hg2_Br(I,J,L) = AD03_Hg2_Br(I,J,L)  + GROSS_OX_BR 
                  AD03_Hg2_OH(I,J,L) = AD03_Hg2_OH(I,J,L)  + GROSS_OX_OH
                  AD03_Hg2_O3(I,J,L) = AD03_Hg2_O3(I,J,L)  + GROSS_OX_O3

                  AD03_Br(I,J,L,1) = AD03_Br(I,J,L,1) + C_BR
                  AD03_Br(I,J,L,2) = AD03_Br(I,J,L,2) + C_BRO


               ENDIF
               
               ! Sea salt diagnostic is 2-D [kg]
               AD03_Hg2_SS(I,J,NN) = AD03_Hg2_SS(I,J,NN) + DEP_HG2_SALT

               ! Sea-salt loss rate diagnostic [/s]
               IF ( L == 1 ) THEN
                  AD03_Hg2_SSR(I,J) = AD03_Hg2_SSR(I,J) + K_SALT
               ENDIF
               
            ENDIF               

         ENDDO               
         
      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! Return to calling program
      END SUBROUTINE CHEM_Hg0_Hg2

!------------------------------------------------------------------------------

      SUBROUTINE RXN_REDOX_NODEP( OLD_Hg0, OLD_Hg2, 
     &     K_OX,    K_RED,   DT,
     &     E_KOX_T, E_KRED_T,
     &     NEW_Hg0, NEW_Hg2, GROSS_OX  )

!
!******************************************************************************
!  Subroutine RXN_REDOX_NODEP calculates new masses of Hg0 and Hg2 for given
!  rates of oxidation and reduction, without any deposition. This is for the
!  free troposphere, or simulations with deposition turned off. 
!  The analytic forms for the solutions were derived in Mathematica. 
!  (cdh 1/17/2008)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) OLD_Hg0  (REAL*8) : Initial mass of Hg(0) [kg]
!  (2 ) OLD_Hg2  (REAL*8) : Initial mass of Hg(II) [kg]
!  (3 ) K_OX     (REAL*8) : 1st order oxidation rate [/s]
!  (4 ) K_RED    (REAL*8) : 1st order reduction rate [/s]
!  (5 ) DT       (REAL*8) : Chemistry time step [s]
!  (6 ) E_KOX_T  (REAL*8) : Precalculated exponential 
!                             exp( - K_OX * DT )  [dimensionless]
!  (7 ) E_KRED_T (REAL*8) : Precalculated exponential 
!                             exp( - K_RED * DT ) [dimensionless]
!
!  Arguments as Output:
!  ============================================================================
!  (1 ) NEW_Hg0  (REAL*8) : Final mass of Hg(0) [kg]
!  (2 ) NEW_Hg2  (REAL*8) : Final mass of Hg(II) [kg]
!  (3 ) GROSS_OX (REAL*8) : Total mass of Hg(0) oxidized [kg] 
! 
!------------------------------------------------------------------------------
!  GENERAL SOLUTION - OXIDATION AND REDUCTION
!
!  The differential equations for oxidation, reduction and deposition are:
!
!     d[Hg0]/dt = -kOx [Hg0] + kRed [Hg2]
! 
!     d[Hg2]/dt =  kOx [Hg0] - kRed [Hg2] 
!
!  The solution is:
!  
!     [Hg0](t) = [Hg0](0) ( kRed + kOx exp(-kt) ) / k
!              + [Hg2](0) ( 1 - exp(-kt) ) kRed / k
!
!     [Hg2](t) = [Hg0](0) ( 1 - exp(-kt) ) kOx / k
!              + [Hg2](0) ( kOx + kRed exp(-kt) ) / k
!
!  where
!
!     k    = kOx + kRed
!
!  In addition, we want to know the gross oxidation flux for diagnostic 
!  reasons:
!     
!     Ox(t) = kOx / k^2 *
!           (   [Hg0](0) ( kOx (1-exp(-kt)) + kRed k t ) 
!             + [Hg2](0) (-kRed(1-exp(-kt)) + kRed k t ) )
!
!  The solutions for [Hg0](t), [Hg2](t), and Ox(t) can become numerically
!  unstable if k << 1. In that case both kOx << 1 and kRed << 1, so the
!  final concentrations will be unchanged from the initial conditions.
!
!------------------------------------------------------------------------------
!  SPECIAL CASE - OXIDATION ONLY   
!
!  The solution simplifies when there is no reduction (kRed=0). Although the
!  exact solution will still apply when kRed=0, we include a separate solution
!  in this subroutine to reduce the computational burden for simulations with
!  kRed=0.
!
!     [Hg0](t) = [Hg0](0) exp(-kOx t)
!
!     [Hg2](t) = [HG2](0) + [Hg0](0) - [Hg0](t) 
!
!  In the absence of reduction, gross and net oxidation are the same
!  
!     Ox(t) = [Hg0](0) - [Hg0](t)
!
!  NOTES:
!  (1 ) Previous versions of the mercury simulation used an inexact solution to
!       the differential equations of chemistry and deposition. The assumptions
!       previously included operator splitting between chemistry and
!       deposition, as well as other numerical approximations. This version
!       simultaneously does chemistry and deposition and uses exact solutions
!       to the differential equations plus some numerical approximations where 
!       the exact solutions may be numerically unstable.
!******************************************************************************
      
      ! Arguments
      REAL*8,  INTENT(IN)  :: OLD_Hg0, OLD_Hg2
      REAL*8,  INTENT(IN)  :: K_OX,    K_RED,   DT
      REAL*8,  INTENT(IN)  :: E_KOX_T, E_KRED_T
      REAL*8,  INTENT(OUT) :: NEW_Hg0, NEW_Hg2, GROSS_OX

      ! Local variables
      REAL*8               :: K,  KT, E_KT

      !=================================================================
      ! RXN_REDOX_NODEP begins here!
      !=================================================================

      IF ( K_RED < SMALLNUM ) THEN

         !=================================================================
         ! Oxidation Only
         !=================================================================

         ! New concentration of Hg0
         NEW_HG0 = OLD_HG0 * E_KOX_T

         ! New concentration of Hg2
         NEW_HG2 = OLD_HG2 + OLD_HG0 - NEW_HG0

         ! Gross oxidation is the same as net oxidation
         GROSS_OX = OLD_HG0 - NEW_HG0

      ELSE

         !=================================================================
         ! Oxidation and Reduction
         !=================================================================

         ! Define useful terms
         K    = K_OX + K_RED
         KT   = K * DT
         E_KT = E_KOX_T * E_KRED_T

         ! Avoid a small divisor
         IF ( K < SMALLNUM ) THEN

            ! When K is very small, then there is no oxidation or reduction
            NEW_HG0 = OLD_HG0
            NEW_HG2 = OLD_HG2
            GROSS_OX = 0D0

         ELSE

            ! New concentration of Hg0
            NEW_HG0 = OLD_HG0 * ( K_RED + K_OX * E_KT ) / K +
     &                OLD_HG2 * K_RED * ( 1D0 - E_KT )  / K

            ! New concentration of Hg2
            NEW_HG2 = OLD_HG2 + OLD_HG0 - NEW_HG0

            ! Gross oxidation
            GROSS_OX = K_OX / ( K**2 ) * 
     &           ( OLD_HG0 * ( K_OX  * (  1D0 - E_KT ) + K_RED * KT ) +
     &             OLD_HG2 * ( K_RED * ( -1D0 + E_KT ) + K_RED * KT ) )
            
         ENDIF

      ENDIF
      
      ! Return to calling program
      END SUBROUTINE RXN_REDOX_NODEP

!------------------------------------------------------------------------------

      SUBROUTINE RXN_REDOX_WITHDEP( OLD_Hg0, OLD_Hg2, 
     &     K_OX,    K_RED,    K_DEP0,        K_DEP2,    DT, 
     &     E_KOX_T, E_KRED_T,
     &     NEW_Hg0, NEW_Hg2,  GROSS_OX,      DEP_HG0,   DEP_HG2 )

!
!******************************************************************************
!  Subroutine RXN_REDOX_WITHDEP calculates new masses of Hg0 and Hg2 for given
!  rates of oxidation, reduction, and deposition of Hg2. This is for the
!  boundary layer. The analytic forms for the solutions were derived in
!  Mathematica. (cdh, 1/17/2008, 7/8/2008)
! 
!  Arguments as Input:
!  ============================================================================
!  (1 ) OLD_Hg0  (REAL*8) : Initial mass of Hg(0) [kg]
!  (2 ) OLD_Hg2  (REAL*8) : Initial mass of Hg(II) [kg]
!  (3 ) K_OX     (REAL*8) : 1st order oxidation rate [/s]
!  (4 ) K_RED    (REAL*8) : 1st order reduction rate [/s]
!  (5 ) K_DEP0   (REAL*8) : Hg(0) deposition rate [/s]
!  (6 ) K_DEP2   (REAL*8) : Hg(II) deposition rate [/s]
!  (7 ) DT       (REAL*8) : Chemistry time step [s]
!  (8 ) E_KOX_T  (REAL*8) : Precalculated exponential 
!                             exp( - K_OX * DT )  [dimensionless]
!  (9 ) E_KRED_T (REAL*8) : Precalculated exponential 
!                             exp( - K_RED * DT ) [dimensionless]
!
!  Arguments as Output:
!  ============================================================================
!  (1 ) NEW_Hg0  (REAL*8) : Final mass of Hg(0) [kg]
!  (2 ) NEW_Hg2  (REAL*8) : Final mass of Hg(II) [kg]
!  (3 ) GROSS_OX (REAL*8) : Total mass of Hg(0) oxidized [kg] 
!  (4 ) DEP_HG0  (REAL*8) : Total mass of Hg(0) deposited [kg] 
!  (5 ) DEP_HG2  (REAL*8) : Total mass of Hg(II) deposited [kg] 
! 
!------------------------------------------------------------------------------
!  GENERAL SOLUTION - OXIDATION, REDUCTION, DEPOSITION
!
!  The differential equations for oxidation, reduction and deposition are:
!
!     d[Hg0]/dt = -kOx [Hg0] + kRed [Hg2] - kDep0 [Hg0]
! 
!     d[Hg2]/dt =  kOx [Hg0] - kRed [Hg2] - kDep2 [Hg2]
!
!  The solution is:
!  
!     [Hg0](t) = [Hg0](0) A(t) /(2R) *
!                         ( R(1+exp(Rt)) + (kOx+kD0-kRed-kD2)(1-exp(Rt)) )
!              + [Hg2](0) A(t) kRed (1-exp(Rt)) / R
!
!     [Hg2](t) = [Hg0](0) A(t) kOx (1-exp(Rt)) / R
!              + [Hg2](0) A(t)/(2R) *
!                         ( R(1+exp(Rt)) + (kRed+kD2-kOx-kD0)(1-exp(Rt)) )
!
!  where
!
!     k    = kOx + kRed + kD0 + kD2
!     R    = sqrt( k^2 - 4( kOx kD2 + kD0( kD2 + kRed ) )
!     A(t) = exp( -(k+R)t / 2 )
!
!  In addition, we want to know the gross oxidation flux for diagnostic
!  reasons:
!     
!     Ox(t) = [Hg0](0) A(t) kOx / (2R ( kOx kD2 + kD0( kD2 + kRed ) ) ) 
!             * ( (1-exp(Rt))( kD2^2 - kD2( kOx - 2kRed ) 
!                              -kD0( kD2 + kRed) +kR( kOx +kRed) )
!               - R ( 1 + exp(Rt) - 2/A(t) )(kD2+kRed) )
!           + [Hg2](0) A(t) kOx kRed / (2R ( kOx kD2 + kD0( kD2 + kRed ) ) )
!             * ( k(1-exp(Rt)) - R( 1 + exp(Rt) + 2/A(t) ) )
!
!     Dep0(t) = Ox(t) kD0 / kOx
!
!     Dep2(t) = [Hg0](0) + [Hg2](0) - [Hg0](t) - [Hg2](t) - Dep0(t)
!
!  The solutions for [Hg0](t), [Hg2](t), and Ox(t) can become numerically
!  unstable if R << 1. In that case we use an approximation that is accurate 
!  to second order in R:
!     (1-exp(Rt))/R ~= -t (1+Rt/2)
!
!------------------------------------------------------------------------------
!
!  The expression for Gross Oxidation has a multiplicative term
!       1/( kOx kD2 + kD0( kD2 + kRed ) ).
!  This is undefined in the following cases:
!    a) kD0 = kOx  = 0 (no Hg(0) loss)
!    b) kD2 = kRed = 0 (no Hg(II) loss)
!    c) kD0 = kD2  = 0 (no deposition; this is handled by RXN_REDOX_NODEP)
!  In addition the solutions are very simple when there is no chemistry
!    d) kOx = kRed = 0
!
!------------------------------------------------------------------------------
!  SPECIAL CASE - DEPOSITION ONLY
!  
!  The solution is very simple when (kOx=kRed=0). We include this case
!  because it reduces computation time.
!
!     [Hg0](t) = [Hg0](0) exp( -kD0 t)
!     [Hg2](t) = [Hg2](0) exp( -kD2 t)
!  
!     Ox(t)   = 0
!     Dep0(t) = [Hg0](0) - [Hg0](t)
!     Dep2(t) = [Hg2](0) - [Hg2](t)
!
!------------------------------------------------------------------------------
!  SPECIAL CASE - OXIDATION AND DEPOSITION, NO REDUCTION   
!
!  The solution simplifies when there is no reduction (kRed=0). Although the
!  exact solution will still apply when kRed=0, we include a separate solution
!  in this subroutine to reduce the computational burden for simulations with
!  kRed=0.
!
!     [Hg0](t) = [Hg0](0) exp( -(kOx+kD0) t)
!
!     [Hg2](t) = [Hg0](0) kOx (exp(-(kOx+kD0) t) - exp(-kD2 t)) / (kD2-kOx-kD0)
!              + [Hg2](0) exp(-kD2 t)
!
!     Ox(t)   = ( [Hg0](0) - [Hg0](t) ) * kOx / ( kOx + kD0 )
!     
!     Dep0(t) = ( [Hg0](0) - [Hg0](t) - Ox(t) )
!     Dep2(t) = ( [Hg0](0) + [Hg2](0) - [Hg0](t) - [Hg2](t) - Dep0(t) )
!
!  The solution for [Hg2](t) can become numerically unstable if 
!  (kD2-kOx-kD0)<<1. In that case, we use an approximation that is accurate
!  to first order in (kD2-kOx-kD0).
!  One can show
!     ( exp(-kt) - exp(-(k+e)t) ) / e ~= -t exp(-kt)
!
!  Therefore, when (kD2-kOx-kD0)<<1, we know kD2~=(kOx+kD0) and thus
!
!     [Hg2](t) ~= [Hg0](0) kOx exp(-kD2 t) (-t) 
!               + [Hg2](0) exp(-kD2 t)
!
!------------------------------------------------------------------------------
!  SPECIAL CASE - REDUCTION AND DEPOSITION, NO OXIDATION   
!
!  The solution simplifies when there is no oxidation (kOx=0). Although the
!  exact solution will still apply when kOx=0, we include a separate solution
!  in this subroutine to reduce the computational burden for simulations with
!  kOx=0.
!
!     [Hg0](t) = [Hg0](0) exp( -kD0 t)
!              + [Hg2](0) kRed (exp( -(kRed+kD2) t) - exp(-kD0 t)) / 
!                   ( kD0 - kRed - kD2 )
!
!     [Hg2](t) = [Hg2](0) exp( -(kRed+kD2) t)
!
!     Ox(t)   = 0
!     Dep2(t) = ( [Hg2](0) - [Hg2](t) ) * kD2 / ( kRed + kD2 )
!     Dep0(t) =   [Hg0](0) + [Hg2](0) - [Hg2](t) - [Hg2](t) - Dep2(t)
!    
!  The solution for [Hg0](t) can become numerically unstable if 
!  (kD0-kRed-kD2)<<1. In that case, we use an approximation that is accurate
!  to first order in (kD0-kRed-kD2).
!  One can show
!     ( exp(-kt) - exp(-(k+e)t) ) / e ~= -t exp(-kt)
!
!  Therefore, when (kD0-kRed-kD2)<<1, we know kD0~=(kRed+kD2) and thus
!
!     [Hg0](t) ~= [Hg0](0) exp(-kD0 t)
!               + [Hg2](0) kRed exp(-kD0 t) (-t) 
!
!  NOTES:
!  (1 ) Previous versions of the mercury simulation used an inexact solution to
!       the differential equations of chemistry and deposition. The assumptions
!       previously included operator splitting between chemistry and
!       deposition, as well as other numerical approximations. This version
!       simultaneously does chemistry and deposition and uses exact solutions
!       to the differential equations plus some numerical approximations where 
!       the exact solutions may be numerically unstable.
!******************************************************************************

      ! Refernces to F90 modules
      USE ERROR_MOD,   ONLY : ERROR_STOP

      ! Arguments
      REAL*8,  INTENT(IN)  :: OLD_Hg0 , OLD_Hg2  , DT
      REAL*8,  INTENT(IN)  :: K_OX    , K_RED    , K_DEP0  ,    K_DEP2
      REAL*8,  INTENT(IN)  :: E_KOX_T , E_KRED_T
      REAL*8,  INTENT(OUT) :: NEW_Hg0 , NEW_Hg2  , GROSS_OX  
      REAL*8,  INTENT(OUT) :: DEP_HG0 , DEP_HG2

      ! Local variables
      REAL*8               :: K, RAD   , E_RAD_T  , AA, QTY1, QTY2
      REAL*8               :: E_KDEP0_T, E_KDEP2_T

      !=================================================================
      ! RXN_REDOX_WITHDEP begins here!
      !=================================================================

      ! Precompute exponential factors [dimensionless]
      E_KDEP0_T = EXP( -K_DEP0 * DT )
      E_KDEP2_T = EXP( -K_DEP2 * DT )

      IF ( (K_RED < SMALLNUM) .AND. (K_OX < SMALLNUM)  ) THEN      

         !=================================================================
         ! No Chemistry, Deposition only
         !=================================================================

         ! New mass of Hg0 and Hg2, [kg]
         NEW_HG0 = OLD_HG0 * E_KDEP0_T
         NEW_HG2 = OLD_HG2 * E_KDEP2_T

         ! Oxidation of Hg0 [kg]
         GROSS_OX = 0D0

         ! Deposited Hg0 and Hg2 [kg]
         DEP_HG0 = OLD_HG0 - NEW_HG0
         DEP_HG2 = OLD_HG2 - NEW_HG2

      ELSE IF (K_RED < SMALLNUM) THEN

         !=================================================================
         ! Oxidation and Deposition only
         !=================================================================

         ! To avoid a small divisor, use an approximation, if necessary
         ! This is accurate to first order in K_OX or K_DEP
         IF ( ABS(K_DEP2 - K_OX - K_DEP0) < SMALLNUM ) THEN
            ! First order approximation
            QTY1 = -DT * E_KDEP2_T
         ELSE
            ! Exact form
            QTY1 = ( E_KOX_T * E_KDEP0_T - E_KDEP2_T ) / 
     &             ( K_DEP2 - K_OX - K_DEP0 )
         ENDIF

         ! New concentration of Hg0 [kg]
         NEW_HG0 = OLD_HG0 * E_KOX_T * E_KDEP0_T
         
         ! New concentration of Hg2 [kg]
         NEW_HG2 = OLD_HG0 * K_OX * QTY1 
     &           + OLD_HG2 * E_KDEP2_T

         ! Gross oxidation mass [kg]
         GROSS_OX = ( OLD_HG0 - NEW_HG0 ) * K_OX / ( K_OX + K_DEP0 )
         GROSS_OX = MAX( GROSS_OX, 0D0 )

         ! Hg0 deposition [kg]
         DEP_HG0 = ( OLD_HG0 - NEW_HG0 - GROSS_OX )
         DEP_HG0 = MAX( DEP_HG0, 0D0 )

         ! Hg2 deposition [kg]
         DEP_HG2 = OLD_HG0 + OLD_HG2 - NEW_HG0 - NEW_HG2 - DEP_HG0
         DEP_HG2 = MAX( DEP_HG2, 0D0 )

      ELSE IF ( K_OX < SMALLNUM ) THEN

         !=================================================================
         ! Reduction and Deposition only
         !=================================================================

         ! To avoid a small divisor, use an approximation, if necessary
         ! This is accurate to first order in K_OX or K_DEP
         IF ( ABS(K_DEP2 + K_RED - K_DEP0) < SMALLNUM ) THEN
            ! First order approximation
            QTY1 = -DT * E_KDEP0_T
         ELSE
            ! Exact form
            QTY1 = ( E_KDEP0_T - E_KDEP2_T * E_KRED_T ) / 
     &             ( K_DEP2 + K_RED - K_DEP0 )
         ENDIF

         ! New concentration of Hg0 [kg]
         NEW_HG0 = OLD_HG0 * E_KDEP0_T 
     &           + OLD_HG2 * K_RED * QTY1
         
         ! New concentration of Hg2 [kg]
         NEW_HG2 = OLD_HG2 * E_KRED_T * E_KDEP2_T 

         ! Gross oxidation mass [kg]
         GROSS_OX = 0D0

         ! Hg2 deposition [kg]
         DEP_HG2 = ( OLD_HG2 - NEW_HG2 ) * K_DEP2 / ( K_DEP2 + K_RED )
         DEP_HG2 = MAX( DEP_HG2, 0D0 )

         ! Hg0 deposition [kg]
         DEP_HG0 = OLD_HG0 + OLD_HG2 - NEW_HG0 - NEW_HG2 - DEP_HG2
         DEP_HG0 = MAX( DEP_HG0, 0D0 )

      ELSE

         !=================================================================
         ! Oxidation, Reduction, and Deposition
         !=================================================================

         ! Define common quantities
         K       = K_OX + K_RED + K_DEP2 + K_DEP0
         RAD     = SQRT( K**2 - 4D0 * 
     &                ( K_OX * K_DEP2 + K_DEP0 * ( K_DEP2 + K_RED ) ) )
         E_RAD_T = EXP( RAD * DT )
         AA      = EXP( -( K + RAD ) * DT / 2D0 )
         QTY2    = K_OX * K_DEP2 + K_DEP0 *( K_DEP2 + K_RED ) 
       
         ! To avoid a small divisor, use an approximation, if necessary
         ! This is accurate to second order in RAD
         IF ( ABS(RAD) < SMALLNUM ) THEN
            ! Second Order Approximation
            QTY1 = -DT * ( 1D0 + RAD * DT / 2D0 ) 
         ELSE
            ! Exact form
            QTY1 = ( 1D0 - E_RAD_T ) / RAD
         ENDIF

         ! New concentration of Hg0 [kg]
         NEW_HG0 = 
     &        OLD_HG0 * AA / 2D0 * 
     &          ( ( 1D0 + E_RAD_T )  
     &          + ( K_OX + K_DEP0 - K_RED - K_DEP2 ) * QTY1 )
     &      - OLD_HG2 * AA * K_RED * QTY1

         ! New concentration of Hg2 [kg]
         NEW_HG2 =
     &       -OLD_HG0 * AA * K_OX * QTY1 
     &       +OLD_HG2 * AA / 2D0 *
     &          ( ( 1D0 + E_RAD_T ) 
     &          - ( K_OX + K_DEP0 - K_RED - K_DEP2 ) * QTY1 ) 

         ! The following conditions will make the oxidation calculation
         ! unstable. The code will require revisions if these conditions
         ! occur, but I don't think they will because of the IF statements
         IF ( AA < SMALLNUM ) THEN
            CALL ERROR_STOP( 'GROSS_OX unstable when AA << 1 ', 
     &           'MERCURY_MOD: RXN_REDOX_WITHDEP')
         ENDIF
         IF ( QTY2 < SMALLNUM ) THEN
            CALL ERROR_STOP( 'GROSS_OX unstable when QTY2 << 1', 
     &           'MERCURY_MOD: RXN_REDOX_WITHDEP')
         ENDIF
         
         ! Gross oxidation mass [kg]
         GROSS_OX =
     &        OLD_HG0 * AA * K_OX / ( 2D0 * QTY2 ) *
     &          ( QTY1 * 
     &              ( K_DEP2**2 - K_DEP2 * ( K_OX - 2D0 * K_RED ) - 
     &                K_DEP0 * ( K_DEP2 + K_RED ) + 
     &                K_RED  * ( K_OX + K_RED ) ) -
     &            ( K_DEP2 + K_RED ) * ( 1D0 + E_RAD_T - 2D0 / AA ) )  
     &      + OLD_HG2 * AA * K_OX * K_RED / ( 2D0 * QTY2 ) * 
     &          ( K * QTY1 - ( 1D0 + E_RAD_T - 2D0 / AA ) )

         ! Deposition of Hg0 and Hg2 [kg]
         DEP_HG0 = GROSS_OX * K_DEP0 / K_OX
         DEP_HG2 = OLD_HG0 + OLD_HG2 - NEW_HG0 - NEW_HG2 - DEP_HG0

    
      ENDIF

      ! Return to calling program
      END SUBROUTINE RXN_REDOX_WITHDEP

!------------------------------------------------------------------------------

      FUNCTION GET_HGBR_RATE( BR, T, P, OH, METHOD ) RESULT( K_HGBR ) 

!
!******************************************************************************
!  Function GET_HGBR_RATE computes the effective 1st order conversion rate 
!  of Hg(0) to Hg(II) via two-step recombination with Br and OH. (cdh, 7/6/06)
!  
!  Arguments as Input:
!  ============================================================================
!  (1  ) BR       (REAL*8 ) : Concentration of atomic Br [molec/cm3] 
!  (2  ) T        (REAL*8 ) : Temperature [K]
!  (3  ) P        (REAL*8 ) : Pressure [hPa]
!  (4  ) OH       (REAL*8 ) : Concentration of OH [molec/cm3]
!  (5  ) METHOD   (REAL*8 ) : Set of rate constants to use in calculation
!
!  Arguments as Output:
!  ============================================================================
!  (1  ) K_HGBR   (REAL*8 ) : Effective 1st order loss rate of Hg(0) [1/s] 
!
!
!  NOTES:
!  ============================================================================
!  This subroutine calculates the net rate of Hg(0) oxidation to Hg(II) through
!     the following reactions. All are gas phase.
!     
!     (1  )  Hg(0) + Br -> HgBr 
!     (2  )  HgBr + M -> Hg(0) + Br
!     (2a )  HgBr     -> Hg(0) + Br
!     (3Br)  HgBr + Br -> HgBr2
!     (3OH)  HgBr + OH -> HgBrOH
      
!
!  References:
!  ============================================================================
!     Ariya, P. A., A. Khalizov, and A. Gidas (2002), Reaction of gaseous
!     mercury with atomic and molecular halogens: kinetics, product studies,
!     and atmospheric implications, Journal of Physical Chemistry A, 106,
!     7310-7320.
!
!     Balabanov, N. B., B. C. Shepler, and K. A. Peterson (2005), Accurate 
!     global potential energy surface and reaction dynamics for the ground
!     state of HgBr2, Journal of Physical Chemistry A, 109(39), 8765-8773.
!
!     Donohoue, D. L., D. Bauer, B. Cossairt, and A. J. Hynes (2006), 
!     Temperature and Pressure Dependent Rate Coefficients for the Reaction
!     of Hg with Br and the Reaction of Br with Br: A Pulsed Laser
!     Photolysis-Pulsed Laser Induced Fluorescence Study, Journal of
!     Physical Chemistry A, 110, 6623-6632.
! 
!     Goodsite, M. E., J. M. C. Plane, and H. Skov (2004), A theoretical
!     study of the oxidation of Hg-0 to HgBr2 in the troposphere, Environmental
!     Science & Technology, 38(6), 1772-1776.
!
!     Holmes, C. D., et al. (2006), Global lifetime of elemental mercury 
!     against oxidation by atomic bromine in the free troposphere, Geophys.
!     Res. Lett., 33(20).
!
!     Khalizov, A. F., B. Viswanathan, P. Larregaray, and P. A. Ariya (2003), 
!     A theoretical study on the reactions of Hg with halogens: Atmospheric
!     implications, Journal of Physical Chemistry A, 107(33), 6360-6365.
!
!******************************************************************************
!
      ! References to F90 modules

      ! Arguments
      REAL*8,             INTENT(IN)  :: BR, OH, T, P
      CHARACTER(LEN=*),   INTENT(IN)  :: METHOD 

      ! Local variables
      REAL*8               :: k1,   k2, k2a, k3br, k3oh
      REAL*8               :: k1G,  k2G
      REAL*8               :: Nair, N0

      !molar gas constant
      REAL*8, PARAMETER    ::  R = 287d0

      ! Function value
      REAL*8               :: K_HGBR

      !=================================================================
      ! GET_HGBR_RATE begins here!
      !=================================================================

      !number density of air [molec/cm3]
      Nair = P * 1d2 / (0.029d0 * R * T ) * 6.02d23 / 1d6

      !standard air density at STP: 1 atm, 273K [molec/cm3]
      N0 = 1013d2 / (0.029d0 * R * 273d0 ) * 6.02d23 / 1d6

      SELECT CASE( METHOD )
      CASE( 'Goodsite' )
         !All rates from Goodsite et al. 2004
         !No HgBr+OH reaction
         k1   = 1.1d-12 * ( T / 298d0 ) ** ( -2.37d0 ) * ( Nair / N0 )
         k2   = 1.2d10  * exp( -8357d0 / T )
         k2a  = 0d0
         k3br = 2.5d-10 * ( T / 298d0 ) ** ( -.57d0 ) 
         k3oh = 0d0
      
      CASE( 'GoodsiteY' )
         !All rates from Goodsite et al. 2004
         !Include HgBr+OH reaction
         k1   = 1.1d-12 * ( T / 298d0 ) ** ( -2.37d0 ) * ( Nair / N0 )
         k2   = 1.2d10  * exp( -8357d0 / T )
         k2a  = 0d0
         k3br = 2.5d-10 * ( T / 298d0 ) ** ( -.57d0 ) 
         k3oh = k3br
      
      CASE( 'Donohoue' )
         !k1 from Donohoue et al. 2006
         !Other rates from Goodsite et al. 2004
         !No HgBr+OH reaction
         k1   = 1.46d-32 * ( T / 298d0 ) ** ( -1.86d0 ) * Nair
         k2   = 1.2d10  * exp( -8357d0 / T )
         k2a  = 0d0
         k3br = 2.5d-10 * ( T / 298d0 ) ** ( -.57d0 ) 
         k3oh = 0d0
     
      CASE( 'DonohoueY' )
         !k1 from Donohoue et al. 2006
         !Other rates from Goodsite et al. 2004
         !Include HgBr+OH reaction
         k1   = 1.46d-32 * ( T / 298d0 ) ** ( -1.86d0 ) * Nair
         k2   = 1.2d10  * exp( -8357d0 / T )
         k2a  = 0d0
         k3br = 2.5d-10 * ( T / 298d0 ) ** ( -.57d0 ) 
         k3oh = k3br
     
      CASE( 'DonohoueYB' )
         !k1 from Donohoue et al. 2006
         !k2 derived from Goodsite et al. 2004 rates
         !   preserving the k2/k1 equilibrium (detailed balance)
         !k3 from Goodsite et al. 2004
         !Include HgBr+OH reaction

         !Goodsite et al. 2004 rates
         k1G  = 1.1d-12 * ( T / 298d0 ) ** ( -2.37d0 ) * ( Nair / N0 )
         k2G  = 1.2d10  * exp(-8357d0/T)

         k1   = 1.46d-32 * ( T / 298d0 ) ** ( -1.86d0 ) * Nair 
         k2   = k2G * k1 / k1G 
         k2a  = 0d0
         k3br = 2.5d-10 * ( T / 298d0 ) ** ( -.57d0 ) 
         k3oh = k3br
      
      CASE( 'Balabanov' )
         !k1, k2 from Goodsite et al. 2004
         !k2a from Balabanov et al. 2005
         !k3 from the zero pressure limit in Balabanov et al. 2005
         k1   = 1.1d-12 * ( T / 298d0 ) ** ( -2.37d0 ) * ( Nair / N0 )
         k2   = 1.2d10  * exp( -8357d0 / T )
         k2a  = 3.9d-11
         k3br = 3d-11 
         k3oh = k3br
     
      CASE( 'KhalizovB' )
         !k1 from Khalizov et al. 2003
         !k2 derived from Goodsite et al. 2004 rates
         !   preserving the k2/k1 equilibrium (detailed balance)
         !k3 from Goodsite et al. 2004
         !Include HgBr+OH reaction

         !Goodsite et al. 2004 rates
         k1G  = 1.1d-12 * ( T / 298d0 ) ** ( -2.37d0 ) * ( Nair / N0 )
         k2G  = 1.2d10  * exp( -8357d0 / T )

         k1   = 1.0d-12 * exp( 209d0 / T ) * ( Nair / N0 )
         k2   = k2G * k1 / k1G 
         k2a  = 0d0
         k3br = 2.5d-10 * ( T / 298d0 ) ** ( -.57d0 ) 
         k3oh = k3br
      
      CASE( 'AriyaB' )
         !k1 from Ariya et al. 2002
         !k2 derived from Goodsite et al. 2004 rates
         !   preserving the k2/k1 equilibrium (detailed balance)
         !k3 from Goodsite et al. 2004
         !Include HgBr+OH reaction

         !Goodsite et al. 2004 rates
         k1G  = 1.1d-12 * ( T / 298d0 ) ** ( -2.37d0 ) * ( Nair / N0 )
         k2G  = 1.2d10  * exp( -8357d0 / T )

         k1   = 3.2d-12 * ( Nair / N0 )
         k2   = k2G * k1 / k1G 
         k2a  = 0d0
         k3br = 2.5d-10 * ( T / 298d0 ) ** ( -.57d0 ) 
         k3oh = k3br

      CASE( 'DonohoueYBBalabanov' )
         !k1 from Donohoue et al. 2006
         !k2 derived from Goodsite et al. 2004 rates
         !   preserving the k2/k1 equilibrium (detailed balance)
         !k2a from Balabanov et al. 2005
         !k3 from Goodsite et al. 2004
         !Include HgBr+OH reaction

         !Goodsite et al. 2004 rates
         k1G  = 1.1d-12 * ( T / 298d0 ) ** ( -2.37d0 ) * ( Nair / N0 )
         k2G  = 1.2d10  * exp(-8357d0/T)

         k1   = 1.46d-32 * ( T / 298d0 ) ** ( -1.86d0 ) * Nair
         k2   = k2G * k1 / k1G 
         k2a  = 3.9d-11
         k3br = 2.5d-10 * ( T / 298d0 ) ** ( -.57d0 ) 
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

!------------------------------------------------------------------------------

      SUBROUTINE CHEM_HgP( V_DEP_HgP )
!
!******************************************************************************
!  Subroutine CHEM_HgP is the chemistry subroutine for HgP (particulate
!  mercury.  HgP is lost via dry deposition. (eck, bmy, 12/7/04, 4/6/06)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) V_DEP_HgP (REAL*8) : Dry deposition velocity for Hg(II) [cm/s]
!     
!  NOTES:
!  (1 ) Removed references to AD44 and "CMN_DIAG".  Now compute drydep for all
!        levels within the PBL.  Now references ND03, LD03 from "diag03_mod.f".
!        (bmy, 2/24/05)
!  (2 ) Now references XNUMOL & XNUMOLAIR from "tracer_mod.f" (bmy, 10/25/05)
!  (3 ) Now use ID_HgP index array from "tracerid_mod.f". (cdh, bmy, 1/9/06)
!******************************************************************************
!
      ! Refernces to F90 modules
      USE DIAG03_MOD,   ONLY : ND03, LD03
      USE LOGICAL_MOD,  ONLY : LSPLIT
      USE PBL_MIX_MOD,  ONLY : GET_FRAC_UNDER_PBLTOP, GET_PBL_MAX_L
      USE TIME_MOD,     ONLY : GET_TS_CHEM
      USE TRACER_MOD,   ONLY : STT,       XNUMOL
      USE TRACERID_MOD, ONLY : ID_HgP,    N_Hg_CATS

      USE CMN_SIZE_MOD  ! Size parameters

      ! Arguments
      REAL*8, INTENT(IN) :: V_DEP_HgP(IIPAR,JJPAR)

      ! Local variables
      INTEGER            :: I,      J,    L,          NN,  PBL_MAX
      REAL*8             :: DTCHEM, E_KT, F_UNDER_TOP

      !=================================================================
      ! CHEM_HgP begins here!
      !=================================================================

      ! Chemistry timestep [s]
      DTCHEM  = GET_TS_CHEM() * 60d0

      ! Maximum extent of the PBL
      PBL_MAX = GET_PBL_MAX_L()

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, F_UNDER_TOP, E_KT, NN )
!$OMP+SCHEDULE( DYNAMIC )
      DO L = 1, PBL_MAX
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Fraction of box (I,J,L) underneath the PBL top [unitless]
         F_UNDER_TOP = GET_FRAC_UNDER_PBLTOP( I, J, L )

         ! If we are in the PBL and there is nonzero drydep vel ...
         IF ( F_UNDER_TOP > 0d0 .and. V_DEP_HgP(I,J) > 0d0 ) THEN

            ! Pre-compute the exponential term for use below [unitless]
            E_KT = EXP( -V_DEP_HgP(I,J) * DTCHEM * F_UNDER_TOP  )

            ! Compute new conc of HgP tracers after drydep [kg]
            DO NN = 1, N_Hg_CATS
               
               ! Do dry deposition (skip undefined HgP tagged tracers)
               IF ( ID_HgP(NN) > 0 ) THEN
                  CALL RXN_HgP_DRYD( I, J, L, ID_HgP(NN), E_KT, DTCHEM )
               ENDIF

            ENDDO
         ENDIF
      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! Return to calling program
      END SUBROUTINE CHEM_HgP

!------------------------------------------------------------------------------

      SUBROUTINE RXN_HgP_DRYD( I, J, L, N, E_KT, DTCHEM )
!
!******************************************************************************
!  Subroutine RXN_Hg0_DRYD computes the new concentration of HgP 
!  after dry deposition.  ND44 diagnostics are also archived. 
!  (eck, bmy, 12/14/04, 4/6/06)
! 
!  Arguments as Input:
!  ============================================================================
!  (1-3) I, J, L  (INTEGER) : GEOS-CHEM lon, lat, alt grid box indices
!  (4  ) N        (INTEGER) : Index for HgP total or tagged tracers
!  (5  ) E_KT     (REAL*8 ) : Value of EXP( -KT ) for drydep rxn [unitless]
!  (6  ) DTCHEM   (INTEGER) : Chemistry timestep [s]
!
!  NOTES:
!  (1  ) Remove references to "diag_mod.f" and "CMN_DIAG".  Now save drydep
!         fluxes into T44 array. (bmy, 2/24/05)
!  (2  ) Now references XNUMOL from "tracer_mod.f" (bmy, 10/25/05)
!  (3  ) Now uses ID_HgP and ID_Hg_tot from "tracerid_mod.f" (cdh, bmy,4/6/06)
!  01 Mar 2012 - R. Yantosca - Now use GET_AREA_CM2(I,J,L) from grid_mod.F90
!******************************************************************************
!
      ! References to F90 modules
      USE DIAG_MOD        , ONLY : AD44
      USE GRID_MOD        , ONLY : GET_AREA_CM2
      USE TRACER_MOD      , ONLY : STT, XNUMOL
      USE DEPO_MERCURY_MOD, ONLY : ADD_HgP_DD
      USE TRACERID_MOD    , ONLY : IS_HgP
      USE LOGICAL_MOD     , ONLY : LNLPBL

      USE CMN_SIZE_MOD   ! Size parameters
      USE CMN_DIAG_MOD   ! ND44 

      ! Arguments
      INTEGER, INTENT(IN) :: I, J, L, N
      REAL*8,  INTENT(IN) :: DTCHEM, E_KT

      ! Local variables
      INTEGER             :: NN
      REAL*8              :: AREA_CM2, DRYDEP, OLD_HgP, NEW_HgP

      !=================================================================
      ! RXN_HgP_DRYD begins here!
      !=================================================================

      ! Error check tracer number
      IF ( N < 1 ) RETURN

      ! Initial concentration of HgP [kg]
      OLD_HgP      = MAX( STT(I,J,L,N), SMALLNUM )

      ! New concentration of HgP after drydep [kg]
      NEW_HgP      = MAX( ( OLD_HgP * E_KT ), SMALLNUM )

      ! Save new concentration of HgP in STT [kg]
      STT(I,J,L,N) = NEW_HgP

      ! Archive Hg(P) lost to drydep [kg] for the ocean mercury flux routines
      ! in "ocean_mercury_mod.f" if necessary.  Do not call ADD_HgP_DD if the 
      DRYDEP       = OLD_HgP - NEW_HgP  
      IF ( IS_HgP(N)) THEN
        CALL ADD_HgP_DD( I, J, N, DRYDEP )
      ENDIF

      !=================================================================
      ! ND44 diagnostic: drydep flux of Hg(II) [molec/cm2/s]
      !=================================================================
      IF ( ND44 > 0 .AND. (.NOT. LNLPBL) ) THEN
         
         ! Grid box surface area [cm2]
         AREA_CM2     = GET_AREA_CM2( I, J, L )

         ! Amt of Hg(II) lost to drydep [molec/cm2/s]
         DRYDEP       = OLD_HgP - NEW_HgP  
         DRYDEP       = ( DRYDEP * XNUMOL(N) ) / ( AREA_CM2 * DTCHEM )

         ! Archive Hg(II) drydep flux in AD44 array [molec/cm2/s]
         AD44(I,J,N,1) = AD44(I,J,N,1) + DRYDEP

      ENDIF
     
      ! Return to calling program
      END SUBROUTINE RXN_HgP_DRYD

!------------------------------------------------------------------------------

#if defined( DEVEL )
      SUBROUTINE EMISSMERCURY( LOCAL_MET )
#else
      SUBROUTINE EMISSMERCURY

#endif
!
!******************************************************************************
!  Subroutine EMISSMERCURY is the driver routine for mercury emissions.
!  (eck, bmy, 12/7/04, 4/6/06)
! 
!  NOTES:
!  (1 ) Now call OCEAN_MERCURY_FLUX from "ocean_mercury_mod.f" to compute 
!        the emissions of Hg0 from the ocean instead of reading it from disk.
!        (sas, bmy, 1/20/05)
!  (2 ) Now no longer call COMPUTE_FEMIS, since we can get the same information
!        from routine GET_FRAC_OF_PBL in "pbl_mix_mod.f" (bmy, 2/22/05)
!  (3 ) Now modified for new ocean mercury module. (cdh, sas, bmy, 4/6/06)
!  01 Mar 2012 - R. Yantosca - Now use GET_XMID(I,J,L) from grid_mod.F90
!  01 Mar 2012 - R. Yantosca - Now use GET_YMID(I,J,L) from grid_mod.F90
!******************************************************************************
!
      ! References to F90 modules
      USE ERROR_MOD,         ONLY : DEBUG_MSG
      USE LOGICAL_MOD,       ONLY : LPRT, LDYNOCEAN, LNLPBL 
      USE OCEAN_MERCURY_MOD, ONLY : OCEAN_MERCURY_FLUX
      USE DEPO_MERCURY_MOD,  ONLY : RESET_HG_DEP_ARRAYS
      USE TIME_MOD,          ONLY : GET_MONTH, ITS_A_NEW_MONTH
      USE TRACER_MOD,        ONLY : STT
      USE TRACERID_MOD,      ONLY : N_HG_CATS
      USE DIAG03_MOD,        ONLY : AD03_nat
      USE VDIFF_PRE_MOD,     ONLY : EMIS_SAVE 
      USE LAND_MERCURY_MOD,  ONLY : LAND_MERCURY_FLUX, VEGEMIS
      USE LAND_MERCURY_MOD,  ONLY : SOILEMIS, BIOMASSHG
      USE LAND_MERCURY_MOD,  ONLY : SNOWPACK_MERCURY_FLUX
      USE GRID_MOD,          ONLY : GET_XMID, GET_YMID

! Added for GTMM (ccc, 11/19/09)
      USE LOGICAL_MOD,       ONLY : LGTMM
      USE LAND_MERCURY_MOD,  ONLY : GTMM_DR
      USE TRACERID_MOD,      ONLY : ID_Hg_tot, ID_Hg0

#if defined( DEVEL )
      USE GC_TYPE_MOD,       ONLY : GC_MET_LOCAL
#endif

      USE CMN_SIZE_MOD          ! Size parameters

      ! Arguments
#if defined( DEVEL )
      TYPE(GC_MET_LOCAL), INTENT(IN) :: LOCAL_MET  ! Obj w/ met fields
#endif

      ! Local variables
      LOGICAL, SAVE :: FIRST = .TRUE. 
      INTEGER       :: THISMONTH, I, J, N

      !=================================================================
      ! EMISSMERCURY begins here!
      !=================================================================

      ! First-time initialization
      IF ( FIRST ) THEN

         ! Read anthro, ocean, land emissions of Hg from disk
         CALL MERCURY_READYR

         ! Reset first-time flag
         FIRST = .FALSE.
      ENDIF

      !=================================================================
      ! Call emission routines for Hg(0), Hg(II), and Hg(P)
      !=================================================================
    
#if defined( DEVEL )
      ! Ocean flux of Hg(0)
      IF ( LDYNOCEAN ) THEN
         
         ! Set to zero to clear emissions from previous time step
         ! (cdh, 4/30/09)
         EHg0_oc = 0d0

         CALL OCEAN_MERCURY_FLUX( EHg0_oc, LOCAL_MET )
         IF ( LPRT ) CALL DEBUG_MSG( '### EMISSMERCURY: a OCEAN_FLUX' )

      ENDIF

      IF ( LGTMM ) THEN
         IF ( ITS_A_NEW_MONTH() ) THEN
            CALL GTMM_DR( EHg0_gtm(:,:), LOCAL_MET )
            IF ( LPRT ) CALL DEBUG_MSG( '### EMISSMERCURY: a GTMM' )
         ENDIF
  
      ELSE
         
         CALL LAND_MERCURY_FLUX ( EHg0_ln, LHGSNOW, LOCAL_MET )
         IF ( LPRT ) CALL DEBUG_MSG( '### EMISSMERCURY: a LAND_FLUX' )
      
         ! Bug fix: VEGEMIS shouldn't be tied to GCAP emissions
         ! (jaf, eds, 4/1/11)
         CALL VEGEMIS( LVEGEMIS, EHg0_dist, EHg0_vg, LOCAL_MET )
         IF ( LPRT ) CALL DEBUG_MSG( '### EMISSMERCURY: a VEGEMIS' )
      
         CALL SOILEMIS( EHg0_dist, EHg0_so, LOCAL_MET )
         IF ( LPRT ) CALL DEBUG_MSG( '### EMISSMERCURY: a SOILEMIS' )

      ENDIF

      CALL SNOWPACK_MERCURY_FLUX ( EHg0_snow, LHGSNOW, LOCAL_MET )
      IF ( LPRT ) CALL DEBUG_MSG( '### EMISSMERCURY: a SNOW_FLUX' )

#else

      ! Ocean flux of Hg(0)
      IF ( LDYNOCEAN ) THEN
         
         ! Set to zero to clear emissions from previous time step
         ! (cdh, 4/30/09)
         EHg0_oc = 0d0

         CALL OCEAN_MERCURY_FLUX( EHg0_oc )
         IF ( LPRT ) CALL DEBUG_MSG( '### EMISSMERCURY: a OCEAN_FLUX' )

      ENDIF

      IF ( LGTMM ) THEN
         IF ( ITS_A_NEW_MONTH() ) THEN
            CALL GTMM_DR(EHg0_gtm(:,:))
            IF ( LPRT ) CALL DEBUG_MSG( '### EMISSMERCURY: a GTMM' )
         ENDIF
  
      ELSE
         
         CALL LAND_MERCURY_FLUX ( EHg0_ln, LHGSNOW )
         IF ( LPRT ) CALL DEBUG_MSG( '### EMISSMERCURY: a LAND_FLUX' )
      
         ! Bug fix: VEGEMIS shouldn't be tied to GCAP emissions
         ! (jaf, eds, 4/1/11)
         CALL VEGEMIS( LVEGEMIS, EHg0_dist, EHg0_vg )
         IF ( LPRT ) CALL DEBUG_MSG( '### EMISSMERCURY: a VEGEMIS' )
      
         CALL SOILEMIS( EHg0_dist, EHg0_so )
         IF ( LPRT ) CALL DEBUG_MSG( '### EMISSMERCURY: a SOILEMIS' )

      ENDIF

      CALL SNOWPACK_MERCURY_FLUX ( EHg0_snow, LHGSNOW )
      IF ( LPRT ) CALL DEBUG_MSG( '### EMISSMERCURY: a SNOW_FLUX' )

#endif  ! DEVEL
      
      CALL BIOMASSHG( EHg0_bb )
      IF ( LPRT ) CALL DEBUG_MSG( '### EMISSMERCURY: a BIOMASS' )

      IF ( LnoUSAemis ) THEN

         ! No Anthropogenic emissions over USA (cdh, 5/6/2009)
         DO J = 1, JJPAR
         DO I = 1, IIPAR

            IF ( GET_XMID( I, J, 1 ) > -125 .AND. 
     &           GET_XMID( I, J, 1 ) < -65  .AND.
     &           GET_YMID( I, J, 1 ) >  25  .AND. 
     &           GET_YMID( I, J, 1 ) < 50 ) THEN
               
               EHg0_An(I,J) = 0d0
               EHg2_An(I,J) = 0d0
               EHgP_An(I,J) = 0d0
               
               EHg0_bb(I,J) = 0d0
               EHg0_am(I,J) = 0d0
               
            ENDIF

         ENDDO
         ENDDO

      ENDIF

      CALL RESET_HG_DEP_ARRAYS
      IF ( LPRT ) CALL DEBUG_MSG( '### EMISSMERCURY: ' //
     &     'a RESET_HG_DEP_ARRAYS' )
      
      ! If we are using the non-local PBL mixing,
      ! we need to initialize the EMIS_SAVE array (cdh, 08/27/09)
      IF (LNLPBL) EMIS_SAVE = 0d0

      ! Add Hg(0) source into STT [kg]
      CALL SRCHg0
      IF ( LPRT ) CALL DEBUG_MSG( '### EMISSMERCURY: a SRCHg0' )

      ! Add Hg(II) source into STT [kg]
      CALL SRCHg2
      IF ( LPRT ) CALL DEBUG_MSG( '### EMISSMERCURY: a SRCHg2' )

      ! Add HgP source into STT [kg]
      CALL SRCHgP
      IF ( LPRT ) CALL DEBUG_MSG( '### EMISSMERCURY: a SRCHgP' )
      
      ! Return to calling program 
      END SUBROUTINE EMISSMERCURY

!----------------------------------------------------------------------------

      SUBROUTINE EMITHG( I, J, L, ID, E_HG )
!
!*****************************************************************************
!  Subroutine EMITHG directs emission either to STT directly or to EMIS_SAVE
!  for use by the non-local PBL mixing. This is a programming convenience.
!  (cdh, 08/27/09)
! 
!  Arguments as Input/Output:
!  ============================================================================
!  (1-3 ) I, J, L (INTEGER) : Grid box indices for lon, lat, vertical level
!  (4   ) ID (INTEGER)      : Tracer ID number in STT
!  (5   ) E_HG (REAL*8)     : Emitted mass of Hg [kg]
!
!******************************************************************************
!

      ! Reference to diagnostic arrays
      USE TRACER_MOD,    ONLY : STT
      USE LOGICAL_MOD,   ONLY : LNLPBL
      USE VDIFF_PRE_MOD, ONLY : EMIS_SAVE

      ! Local variables
      INTEGER, INTENT(IN)   :: I, J, L, ID
      REAL*8,  INTENT(IN)   :: E_Hg

      !=================================================================
      ! EMITHG begins here!
      !=================================================================

      ! Save emissions for non-local PBL mixing or emit directly.
      ! Make sure that emitted mass is non-negative
      ! This is hear only for consistency with old code which warned of
      ! underflow error (cdh, 08/27/09)
      IF (LNLPBL) THEN
         EMIS_SAVE(I,J,ID) = EMIS_SAVE(I,J,ID) + MAX( E_HG, 0D0 )
      ELSE
         STT(I,J,L,ID) = STT(I,J,L,ID) + MAX( E_HG, 0D0 )
      ENDIF

      END SUBROUTINE EMITHG

!-----------------------------------------------------------------------------

      SUBROUTINE SRCHg0
!
!******************************************************************************
!  Subroutine SRCHg0 is the subroutine for Hg(0) emissions.  
!  Emissions of Hg(0) will be distributed throughout the boundary layer. 
!  (eck, cdh, bmy, 1/21/05, 4/6/06)
! 
!  Arguments as Input/Output:
!  ============================================================================
!  (1 ) TC (REAL*8) : Tracer concentration of Hg(0) [kg]
!
!  NOTES:
!  (1 ) Now use diagnostic arrays from "diag03_mod.f" (bmy, 1/21/05)
!  (2 ) Now references GET_FRAC_OF_PBL and GET_PBL_MAX_L from "pbl_mix_mod.f".
!        Remove reference to FEMIS. (bmy, 2/22/05)
!  (3 ) EHg0_an is now a 2-D array.  Modified for new ocean mercury module.
!        Now use ID_Hg0 index array from "tracerid_mod.f".  Now make sure
!        STT does not underflow. (cdh, bmy, 4/6/06)
!******************************************************************************
!
      ! Reference to diagnostic arrays
      USE DIAG03_MOD,   ONLY : AD03, ND03, AD03_nat
      USE ERROR_MOD,    ONLY : ERROR_STOP
      USE LOGICAL_MOD,  ONLY : LSPLIT, LPREINDHG, LGTMM
      USE PBL_MIX_MOD,  ONLY : GET_FRAC_OF_PBL, GET_PBL_MAX_L
      USE TIME_MOD,     ONLY : GET_TS_EMIS
      USE TRACERID_MOD, ONLY : ID_Hg_tot, ID_Hg_na, ID_Hg_eu, ID_Hg_as
      USE TRACERID_MOD, ONLY : ID_Hg_rw,  ID_Hg_oc, ID_Hg_ln, ID_Hg_nt
      USE TRACERID_MOD, ONLY : ID_Hg0
      
      USE CMN_SIZE_MOD     ! Size parameters

      ! Local variables
      INTEGER               :: I,      J,    L,        N,    PBL_MAX
      REAL*8                :: DTSRCE, E_Hg, F_OF_PBL, T_Hg, T_Hg_An

      !=================================================================
      ! SRCHg0 begins here!
      !=================================================================

      ! Emission timestep [s]
      DTSRCE  = GET_TS_EMIS() * 60d0

      ! Maximum extent of the PBL [model levels]
      PBL_MAX = GET_PBL_MAX_L() 
           
      ! Loop over grid boxes
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, N, T_Hg_An, T_Hg, F_OF_PBL, E_Hg)
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         IF ( LPREINDHG ) THEN  !eds

            ! Anthropogenic emissions
            T_Hg_An = 0D0

            ! No biomass burning emissions
            EHg0_bb = 0D0

         ELSE

            ! Compute total anthropogenic Hg(0) emissions
            T_Hg_An = EHg0_an(I,J)

            IF ( LAnthroHgOnly ) THEN
               ! No other emissions
               EHg0_bb = 0D0
               EHg0_oc = 0D0
               EHg0_nt = 0D0
               EHg0_snow = 0D0
               IF ( LGTMM ) THEN
                  EHg0_gtm = 0D0
               ELSE
                  EHg0_ln = 0D0
                  EHg0_vg = 0D0
                  EHg0_so = 0D0
               ENDIF
            ENDIF
            
         ENDIF

         ! Compute total Hg(0) emissions (anthro+oceans+land+natural)
         IF ( LGTMM ) THEN
            T_Hg = T_Hg_An +
     &             EHg0_bb(I,J) +
     &             EHg0_oc(I,J,ID_Hg_tot) + 
     &             EHg0_nt(I,J) +
     &             EHg0_gtm(I,J) +
     &             EHg0_snow(I,J,ID_Hg_tot)
         ELSE
            T_Hg = T_Hg_An +
     &             EHg0_bb(I,J) +
     &             EHg0_oc(I,J,ID_Hg_tot) + 
     &             EHg0_ln(I,J,ID_Hg_tot) +
     &             EHg0_nt(I,J) +
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
            F_OF_PBL        = GET_FRAC_OF_PBL( I, J, L )

            !-----------------
            ! Total Hg tracer
            !-----------------
            N               = ID_Hg0(ID_Hg_tot)
            E_Hg            = F_OF_PBL * T_Hg * DTSRCE
            CALL EMITHG( I, J, L, N, E_Hg )

            !-----------------
            ! Tagged tracers
            !-----------------
            IF ( LSPLIT ) THEN 

               !--------------------
               ! Primary emissions
               !--------------------

               ! Anthro Hg0 by region
               N            = AN_Hg0(I,J)
               E_Hg         = F_OF_PBL * T_Hg_an * DTSRCE
               CALL EMITHG( I, J, L, N, E_Hg )
               

               ! Natural land sources of Hg0
               N            = ID_Hg0(ID_Hg_nt)

               E_Hg         = F_OF_PBL * EHg0_nt(I,J) * DTSRCE
               CALL EMITHG( I, J, L, N, E_Hg )

               !--------------------
               ! Ocean re-emissions
               !--------------------

               ! Anthro re-emission from ocean in N. AMERICA
               N            = ID_Hg0(ID_Hg_na)
               E_Hg         = F_OF_PBL * EHg0_oc(I,J,ID_Hg_na) * DTSRCE
               CALL EMITHG( I, J, L, N, E_Hg )

               ! Anthro re-emission from ocean in EUROPE
               N            = ID_Hg0(ID_Hg_eu)
               E_Hg         = F_OF_PBL * EHg0_oc(I,J,ID_Hg_eu) * DTSRCE
               CALL EMITHG( I, J, L, N, E_Hg )

               ! Anthro re-emission from ocean in ASIA
               N            = ID_Hg0(ID_Hg_as)
               E_Hg         = F_OF_PBL * EHg0_oc(I,J,ID_Hg_as) * DTSRCE
               CALL EMITHG( I, J, L, N, E_Hg )

               ! Anthro re-emission from ocean in REST OF WORLD
               N            = ID_Hg0(ID_Hg_rw)
               E_Hg         = F_OF_PBL * EHg0_oc(I,J,ID_Hg_rw) * DTSRCE
               CALL EMITHG( I, J, L, N, E_Hg )

               ! Re-emission from ocean in OCEAN 
               N            = ID_Hg0(ID_Hg_oc)
               E_Hg         = F_OF_PBL * EHg0_oc(I,J,ID_Hg_oc) * DTSRCE
               CALL EMITHG( I, J, L, N, E_Hg )

               ! Re-emission from ocean in LAND REEMISSION
               N            = ID_Hg0(ID_Hg_ln)
               E_Hg         = F_OF_PBL * EHg0_oc(I,J,ID_Hg_ln) * DTSRCE
               CALL EMITHG( I, J, L, N, E_Hg )

               ! Re-emission from ocean in NATURAL
               N            = ID_Hg0(ID_Hg_nt)
               E_Hg         = F_OF_PBL * EHg0_oc(I,J,ID_Hg_nt) * DTSRCE
               CALL EMITHG( I, J, L, N, E_Hg )
              

               !--------------------
               ! Land re-emissions
               !--------------------

               ! Anthro re-emission from ocean in N. AMERICA
               N            = ID_Hg0(ID_Hg_na)
               E_Hg         = F_OF_PBL * EHg0_ln(I,J,ID_Hg_na) * DTSRCE
               CALL EMITHG( I, J, L, N, E_Hg )

               ! Anthro re-emission from ocean in EUROPE
               N            = ID_Hg0(ID_Hg_eu)
               E_Hg         = F_OF_PBL * EHg0_ln(I,J,ID_Hg_eu) * DTSRCE
               CALL EMITHG( I, J, L, N, E_Hg )

               ! Anthro re-emission from ocean in ASIA
               N            = ID_Hg0(ID_Hg_as)
               E_Hg         = F_OF_PBL * EHg0_ln(I,J,ID_Hg_as) * DTSRCE
               CALL EMITHG( I, J, L, N, E_Hg )

               ! Anthro re-emission from ocean in REST OF WORLD
               N            = ID_Hg0(ID_Hg_rw)
               E_Hg         = F_OF_PBL * EHg0_ln(I,J,ID_Hg_rw) * DTSRCE
               CALL EMITHG( I, J, L, N, E_Hg )

               ! Re-emission from ocean in OCEAN 
               N            = ID_Hg0(ID_Hg_oc)
               E_Hg         = F_OF_PBL * EHg0_ln(I,J,ID_Hg_oc) * DTSRCE
               CALL EMITHG( I, J, L, N, E_Hg )

               ! Re-emission from ocean in LAND REEMISSION
               N            = ID_Hg0(ID_Hg_ln)
               E_Hg         = F_OF_PBL * EHg0_ln(I,J,ID_Hg_ln) * DTSRCE
               CALL EMITHG( I, J, L, N, E_Hg )

               ! Re-emission from ocean in NATURAL
               N            = ID_Hg0(ID_Hg_nt)
               E_Hg         = F_OF_PBL * EHg0_ln(I,J,ID_Hg_nt) * DTSRCE
               CALL EMITHG( I, J, L, N, E_Hg )
               
               !--------------------
               ! Snow re-emissions
               !--------------------

               ! Anthro re-emission from snow in N. AMERICA
               N          = ID_Hg0(ID_Hg_na)
               E_Hg       = F_OF_PBL * EHg0_snow(I,J,ID_Hg_na) * DTSRCE
               CALL EMITHG( I, J, L, N, E_Hg )

               ! Anthro re-emission from snow in EUROPE
               N          = ID_Hg0(ID_Hg_eu)
               E_Hg       = F_OF_PBL * EHg0_snow(I,J,ID_Hg_eu) * DTSRCE
               CALL EMITHG( I, J, L, N, E_Hg )

               ! Anthro re-emission from snow in ASIA
               N          = ID_Hg0(ID_Hg_as)
               E_Hg       = F_OF_PBL * EHg0_snow(I,J,ID_Hg_as) * DTSRCE
               CALL EMITHG( I, J, L, N, E_Hg )

               ! Anthro re-emission from snow in REST OF WORLD
               N          = ID_Hg0(ID_Hg_rw)
               E_Hg       = F_OF_PBL * EHg0_snow(I,J,ID_Hg_rw) * DTSRCE
               CALL EMITHG( I, J, L, N, E_Hg )

               ! Re-emission from snow in OCEAN 
               N          = ID_Hg0(ID_Hg_oc)
               E_Hg       = F_OF_PBL * EHg0_snow(I,J,ID_Hg_oc) * DTSRCE
               CALL EMITHG( I, J, L, N, E_Hg )

               ! Re-emission from snow in LAND REEMISSION
               N          = ID_Hg0(ID_Hg_ln)
               E_Hg       = F_OF_PBL * EHg0_snow(I,J,ID_Hg_ln) * DTSRCE
               CALL EMITHG( I, J, L, N, E_Hg )

               ! Re-emission from snow in NATURAL
               N          = ID_Hg0(ID_Hg_nt)
               E_Hg       = F_OF_PBL * EHg0_snow(I,J,ID_Hg_nt) * DTSRCE
               CALL EMITHG( I, J, L, N, E_Hg )
               
            ENDIF
         ENDDO
        
         !==============================================================
         ! ND03 diagnostic: Total Hg(0) emissions [kg]
         ! 1=anthro; 3=from ocean; 4=land re-emission; 5=natural src
         !==============================================================
         IF ( ND03 > 0 ) THEN
            N = ID_Hg_tot
            IF ( LGTMM ) THEN
               AD03(I,J,1) = AD03(I,J,1) + ( T_Hg_An        * DTSRCE )
               AD03(I,J,3) = AD03(I,J,3) + ( EHg0_oc(I,J,N) * DTSRCE )
               AD03(I,J,4) = AD03(I,J,4) + ( EHg0_gtm(I,J)  * DTSRCE )
               AD03(I,J,5) = AD03(I,J,5) + ( EHg0_nt(I,J) * DTSRCE )
               AD03(I,J,13)= AD03(I,J,13)+ ( EHg0_bb(I,J)*DTSRCE)
               AD03(I,J,18)= AD03(I,J,18)+ ( EHg0_snow(I,J,N)*DTSRCE)
            ELSE
               AD03(I,J,1) = AD03(I,J,1) + ( T_Hg_An        * DTSRCE )
               AD03(I,J,3) = AD03(I,J,3) + ( EHg0_oc(I,J,N) * DTSRCE )
               AD03(I,J,4) = AD03(I,J,4) + ( EHg0_ln(I,J,N) * DTSRCE )
               AD03(I,J,5) = AD03(I,J,5) + ( EHg0_nt(I,J)  * DTSRCE )
               AD03(I,J,13)= AD03(I,J,13)+ ( EHg0_bb(I,J)*DTSRCE)
               AD03(I,J,14)= AD03(I,J,14)+ ( EHg0_vg(I,J)*DTSRCE)
               AD03(I,J,15)= AD03(I,J,15)+ ( EHg0_so(I,J)*DTSRCE)
               AD03(I,J,18)= AD03(I,J,18)+ ( EHg0_snow(I,J,N)*DTSRCE)
            ENDIF

            ! for preindustrial simulation, archive only soil,
            ! CDH- WHY ONLY SOIL??
            ! for present day archive soil, geogenic, biomass burning, 
            ! vegetation, and rapid recycing 
            IF ( LPREINDHG ) THEN !eds
               AD03_nat(I,J,N)=MAX(EHg0_so(I,J)*DTSRCE, SMALLNUM)
            ELSE
               IF ( LGTMM ) THEN
                  AD03_nat(I,J,N) = DTSRCE * (
     &                 EHg0_nt(I,J)     + 
     &                 EHg0_bb(I,J)     +
     &                 EHg0_gtm(I,J)    + 
     &                 EHg0_snow(I,J,N)  ) 
               ELSE
                  AD03_nat(I,J,N) = DTSRCE * (
     &                 EHg0_ln(I,J,N) + 
     &                 EHg0_nt(I,J) + 
     &                 EHg0_bb(I,J) +
     &                 EHg0_vg(I,J) + 
     &                 EHg0_so(I,J) + 
     &                 EHg0_snow(I,J,N) ) 
               ENDIF 
            ENDIF
            
         ENDIF
      ENDDO
      ENDDO
!$OMP END PARALLEL DO
  
      ! Return to calling program
      END SUBROUTINE SRCHg0

!-----------------------------------------------------------------------------

      SUBROUTINE SRCHg2
!
!******************************************************************************
!  Subroutine SRCHg2 is the subroutine for Hg(II) emissions.  
!  Emissions of Hg(II) will be distributed throughout the boundary layer. 
!  (eck, bmy, 12/7/04, 4/6/06)
! 
!  Arguments as Input/Output:
!  ============================================================================
!  (1 ) TC (REAL*8) : Tracer concentration of Hg(II) [kg]
!
!  NOTES:
!  (1 ) Now use diagnostic arrays from "diag03_mod.f" (bmy, 1/21/05)
!  (2 ) Now references GET_FRAC_OF_PBL and GET_PBL_MAX_L from "pbl_mix_mod.f".
!        Remove reference to FEMIS. (bmy, 2/22/05)
!  (3 ) EHg2_an is now a 2-D array.  Now use ID_Hg2 index array from 
!        "tracerid_mod.f".  Now make sure STT does not underflow.
!        (eck, cdh, bmy, 4/6/06)
!******************************************************************************
!
      ! Reference to F90 modules
      USE DIAG03_MOD,   ONLY : AD03, ND03
      USE ERROR_MOD,    ONLY : ERROR_STOP
      USE LOGICAL_MOD,  ONLY : LSPLIT, LPREINDHG
      USE PBL_MIX_MOD,  ONLY : GET_FRAC_OF_PBL, GET_PBL_MAX_L
      USE TIME_MOD,     ONLY : GET_TS_EMIS
      USE TRACERID_MOD, ONLY : ID_Hg_tot, ID_Hg2

      USE CMN_SIZE_MOD     ! Size parameters

      ! Local variables
      INTEGER               :: I,      J,        L,    N,   PBL_MAX
      REAL*8                :: DTSRCE, F_OF_PBL, E_Hg 

      !=================================================================
      ! SRCHg2 begins here!
      !=================================================================
      
      ! Emission timestep [s]
      DTSRCE  = GET_TS_EMIS() * 60d0

      ! Maximum extent of the PBL [model levels]
      PBL_MAX = GET_PBL_MAX_L()

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
               F_OF_PBL           = GET_FRAC_OF_PBL( I, J, L )

               ! Partition total Hg2 into box (I,J,L) [kg]
               E_Hg               = F_OF_PBL * EHg2_an(I,J) * DTSRCE
                        
               !---------------------------
               ! Total anthro Hg(II) [kg]
               !---------------------------
               N               = ID_Hg2(ID_Hg_tot)
               CALL EMITHG( I, J, L, N, E_Hg )

               !---------------------------
               ! Tagged anthro Hg(II) [kg]
               !---------------------------
               IF ( LSPLIT ) THEN 
                  N            = AN_Hg2(I,J)
                  CALL EMITHG( I, J, L, N, E_Hg )

               ENDIF
            ENDDO
            
            !-------------------------------
            ! ND03 diag: Anthro Hg(II) [kg]
            !-------------------------------
            IF ( ND03 > 0 ) THEN
               AD03(I,J,6) = AD03(I,J,6) + ( EHg2_an(I,J) * DTSRCE )
            ENDIF
         ENDDO
         ENDDO
!$OMP END PARALLEL DO

      ENDIF

      ! Return to calling program
      END SUBROUTINE SRCHg2

!-----------------------------------------------------------------------------

      SUBROUTINE SRCHgP
!
!******************************************************************************
!  Subroutine SRCHgP is the subroutine for HgP emissions.  
!  Emissions of HgP will be distributed throughout the boundary layer. 
!  (eck, bmy, 12/7/04, 4/6/06)
! 
!  Arguments as Input/Output:
!  ============================================================================
!  (1 ) TC (REAL*8) : Tracer concentration of Hg(II) [kg]
!
!  NOTES:
!  (1 ) Now use diagnostic arrays from "diag03_mod.f" (bmy, 1/21/05)
!  (2 ) Now references GET_FRAC_OF_PBL and GET_PBL_MAX_L from "pbl_mix_mod.f".
!        Remove reference to FEMIS. (bmy, 2/22/05)
!  (3 ) EHgP_an is now a 2-D array.  Now use ID_HgP index array from 
!        "tracerid_mod.f".  Now make sure STT does not underflow. 
!        (eck, cdh, bmy, 4/6/06) 
!******************************************************************************
!
      ! Reference to diagnostic arrays
      USE DIAG03_MOD,   ONLY : AD03, ND03
      USE ERROR_MOD,    ONLY : ERROR_STOP
      USE LOGICAL_MOD,  ONLY : LSPLIT, LPREINDHG
      USE PBL_MIX_MOD,  ONLY : GET_FRAC_OF_PBL, GET_PBL_MAX_L
      USE TIME_MOD,     ONLY : GET_TS_EMIS 
      USE TRACERID_MOD, ONLY : ID_Hg_tot, ID_HgP

      USE CMN_SIZE_MOD     ! Size paramters

      ! Local variables
      INTEGER               :: I,      J,        L,    N,   PBL_MAX
      REAL*8                :: DTSRCE, F_OF_PBL, E_Hg 

      !=================================================================
      ! SRCHgP begins here!
      !=================================================================

      ! Chemistry timestep [s]
      DTSRCE  = GET_TS_EMIS() * 60d0
      
      ! Maximum extent of the PBL [model levels]
      PBL_MAX = GET_PBL_MAX_L()

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
               N               = ID_HgP(ID_Hg_tot)
!               STT(I,J,L,N)    = STT(I,J,L,N) + E_Hg
!               STT(I,J,L,N)    = MAX( STT(I,J,L,N), SMALLNUM )
               CALL EMITHG( I, J, L, N, E_Hg )

               !------------------------
               ! Tagged anthro HgP [kg]
               !------------------------
               IF ( LSPLIT ) THEN
                  N            = AN_HgP(I,J)
!                  STT(I,J,L,N) = STT(I,J,L,N) + E_Hg
!                  STT(I,J,L,N) = MAX( STT(I,J,L,N), SMALLNUM )
                  CALL EMITHG( I, J, L, N, E_Hg )

               ENDIF
            ENDDO

            !----------------------------
            ! ND03 diag: Anthro HgP [kg]
            !----------------------------
            IF ( ND03 > 0 ) THEN
               AD03(I,J,9)  = AD03(I,J,9) + ( EHgP_an(I,J) * DTSRCE )
            ENDIF
         
         ENDDO 
         ENDDO
!$OMP END PARALLEL DO

      ENDIF
         
      ! Return to calling program
      END SUBROUTINE SRCHgP

!------------------------------------------------------------------------------

      SUBROUTINE MERCURY_READYR
!
!******************************************************************************
!  Subroutine MERCURY_READYR reads the year-invariant emissions for Mercury
!  from anthropogenic, ocean, and land sources. (eck, bmy, 12/6/04, 4/6/06)
!  
!  NOTES:
!  (1 ) Now read data from mercury_200501 subdirectory.  Now compute oceanic 
!        Hg(0) emissions w/ ocean flux module instead of reading them from 
!        disk.  Now use 1985 TAU values. (sas, bmy, 1/20/05) 
!  (2 ) Now can read data for both GEOS and GCAP grids (bmy, 8/16/05)
!  (3 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (4 ) Now read anthro emissions on GEOS 1x1 grid in DATA_1x1_DIR.  Also
!        keep 2x25 and 4x5 files together in DATA_1x1_DIR.  Also now use new
!        land re-emissions files from Noelle Selin. (eck, bmy, 4/6/06)
!  (5 ) Bug fix: emissions are read on GENERIC grid (jaf, 2/8/12)
!  01 Mar 2012 - R. Yantosca - Now use GET_XMID(I,J,L) from grid_mod.F90
!  01 Mar 2012 - R. Yantosca - Now use GET_YMID(I,J,L) from grid_mod.F90
!  13 Mar 2012 - M. Cooper   - Changed regrid algorithm to map_a2a
!  24 May 2012 - R. Yantosca - Fixed minor bugs in map_a2a algorithm
!  24 May 2012 - H. Amos     - updated file paths for soil Hg files
!******************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD,      ONLY : GET_NAME_EXT_2D, GET_RES_EXT
      USE BPCH2_MOD,      ONLY : GET_TAU0,        READ_BPCH2
      USE DIRECTORY_MOD,  ONLY : DATA_DIR_1x1, DATA_DIR
      USE LOGICAL_MOD,    ONLY : LDYNOCEAN, LPREINDHG
      USE REGRID_A2A_MOD, ONLY : DO_REGRID_A2A
      USE TIME_MOD,       ONLY : EXPAND_DATE
      USE TRANSFER_MOD,   ONLY : TRANSFER_2D
      USE GRID_MOD,       ONLY : GET_XMID,    GET_YMID
      USE ERROR_MOD,      ONLY : ERROR_STOP

      USE CMN_SIZE_MOD         ! Size parameters

      ! Local variables
      INTEGER                 :: NYMD, I, J
      REAL*4                  :: ARRAY(IIPAR,JJPAR,1)
      REAL*4                  :: ARRAY1(I1x1,J1x1-1,1)
      REAL*8, TARGET          :: GEN_1x1(I1x1,J1x1-1)
      REAL*8                  :: XTAU
      REAL*8                  :: X, Y
      REAL*8, PARAMETER       :: SEC_PER_YR = 365.25d0 * 86400d0
      CHARACTER(LEN=255)      :: FILENAME
      CHARACTER(LEN=255)      :: LLFILENAME
      REAL*8, POINTER         :: INGRID(:,:) => NULL()

      !=================================================================
      ! MERCURY_READYR begins here! 
      !=================================================================


      !=================================================================
      ! Anthropogenic Emissions
      !=================================================================

      ! No anthropogenic emissions in preindustrial simulation
      IF ( LPREINDHG ) THEN

         EHg0_an = 0D0  !eds
         EHg2_an = 0D0
         EHgP_an = 0D0

      ELSE

         !---------------------------
         ! Hg(0) emissions, except mining [kg/s]
         !---------------------------

         ! File with lat/lon edges for regridding
         LLFILENAME = TRIM( DATA_DIR_1x1) //
     &                'MAP_A2A_Regrid_201203/latlon_generic_new.txt'

         IF (LGCAPEMIS) THEN
            ! Use Bess' GCAP emissions
   
            ! Mercury data is either for 2006
            NYMD     = ( 2006 * 10000 ) + 0101
            XTAU     = GET_TAU0( 1, 1, 2006 )

            ! Filename for anthropogenic mercury source
            FILENAME = TRIM( DATA_DIR_1x1 )       // 
     &           'mercury_201002/GEIA_Streets_Hg0.geos.1x1.YYYY'

            ! Add year to the filename
            CALL EXPAND_DATE( FILENAME, NYMD, 000000 )

            ! Echo info
            WRITE( 6, 100 ) TRIM( FILENAME )
 100        FORMAT( '     - MERCURY_READYR: Reading ', a )

            ! Read data in [kg/yr]
            CALL READ_BPCH2( FILENAME, 'HG-SRCE', 1,
     &                       XTAU,      I1x1,     J1x1-1,
     &                       1,         ARRAY1,   QUIET=.TRUE. )

            ! Cast to REAL*8 before regridding
            GEN_1x1(:,:) = ARRAY1(:,:,1)

            ! Regrid from GENERIC 1x1 to the current grid
            INGRID  => GEN_1x1(:,:)
            CALL DO_REGRID_A2A( LLFILENAME, I1x1,    J1x1-1,
     &                          INGRID,     EHg0_an, PERAREA=1 )

            ! Convert from [kg/yr] to [kg/s]
            EHg0_an = EHg0_an / SEC_PER_YR            
        
            ! Free pointer
            NULLIFY( INGRID )

         ELSEIF (LGEIA05) THEN
            ! Use GEIA 2005 emissions inventory (eck, 10/19/11)
   
            ! Mercury data is for 2005
            NYMD     = ( 2005 * 10000 ) + 0101
            XTAU     = GET_TAU0( 1, 1, 2005 )

            ! Filename for anthropogenic mercury source
            FILENAME = TRIM( DATA_DIR_1x1 )       // 
     &           'mercury_201002/GEIA_Hg0.geos.1x1.YYYY'

            ! Add year to the filename
            CALL EXPAND_DATE( FILENAME, NYMD, 000000 )

            ! Echo info
            WRITE( 6, 100 ) TRIM( FILENAME )

            ! Read data in [kg/yr]
            CALL READ_BPCH2( FILENAME, 'HG-SRCE', 1,
     &                       XTAU,      I1x1,     J1x1-1,
     &                       1,         ARRAY1,   QUIET=.TRUE. )

            ! Cast to REAL*8 before regridding
            GEN_1x1(:,:) = ARRAY1(:,:,1)

            ! Regrid from GENERIC 1x1 to the current grid
            INGRID  => GEN_1x1(:,:)
            CALL DO_REGRID_A2A( LLFILENAME, I1x1,    J1x1-1,
     &                          INGRID,     EHg0_an, PERAREA=1 )

            ! Convert from [kg/yr] to [kg/s]
            EHg0_an = EHg0_an / SEC_PER_YR    

            ! For GEIA 2005, change Hg0:Hg2:HgP ratio
            ! from 50:40:10 to 86.5:9.9:3.6 for fossil fuel emissions.
            ! Fossil fuels make up ~47% of total global anthro emis
            ! Refs: Amos et al. (2011, ACPD), Y. Zhang et al. (2011, in prep)
            ! Implemented 10/19/11 by eck
            EHg0_an=1.24d0*EHg0_an
            
            ! Free pointer
            NULLIFY( INGRID )

         ELSE
            ! Otherwise use Selin et al. 2008 emissions

            ! Mercury data is either for 1995 or 2000
            NYMD     = ( ANTHRO_Hg_YEAR * 10000 ) + 0101
            XTAU     = GET_TAU0( 1, 1, ANTHRO_Hg_YEAR )

            ! Filename for anthropogenic mercury source
            FILENAME = TRIM( DATA_DIR_1x1 ) // 
     &           'mercury_200511/GEIA_Hg0.geos.1x1.YYYY'

            ! Add year to the filename
            CALL EXPAND_DATE( FILENAME, NYMD, 000000 )

            ! Echo info
            WRITE( 6, 100 ) TRIM( FILENAME )

            ! Read data in [kg/yr]
            CALL READ_BPCH2( FILENAME, 'HG-SRCE', 1,
     &                       XTAU,      I1x1,     J1x1-1,
     &                       1,         ARRAY1,   QUIET=.TRUE. )
            
            ! Cast to REAL*8 before regridding
            GEN_1x1(:,:) = ARRAY1(:,:,1)

            ! Regrid from GENERIC 1x1 to the current grid
            INGRID  => GEN_1x1(:,:)
            CALL DO_REGRID_A2A( LLFILENAME, I1x1,    J1x1-1,
     &                          INGRID,     EHg0_an, PERAREA=1 )

            ! Convert from [kg/yr] to [kg/s]
            EHg0_an = EHg0_an / SEC_PER_YR

            ! Scale up the inventory, as explained in Selin et al., GBC 2008
            DO J = 1, JJPAR
            DO I = 1, IIPAR

               ! Grid box longitude [degrees]
               X = GET_XMID( I, J, 1 )

               ! Grid-box latitude [degrees]
               Y = GET_YMID( I, J, 1 )

               ! Asian Anthro Hg increase 50%
               IF ( ( X >= 70.0 .and. X < 152.5 )  .and.
     &              ( Y >=  8.0 .and. Y <  45.0 ) ) THEN

                  EHg0_an(I,J) = EHg0_an(I,J) * 1.5d0

               ELSE

                  ! otherwise, 30% increase
                  EHg0_an(I,J) = EHg0_an(I,J) * 1.3d0

               ENDIF 

            ENDDO
            ENDDO 

            ! Free pointer
            NULLIFY( INGRID )

         ENDIF

         !---------------------------
         ! Artisanal mining Hg(0) emissions [kg/s]
         !---------------------------

         ! Read artisanal mining and add to anthropogenic
         ! Use "generic" year 1985 for TAU values
         XTAU     = GET_TAU0( 1, 1, 1985 )

         ! Filename for artisanal mining emissions (4x5 only for now)
         FILENAME=TRIM( DATA_DIR ) // 'mercury_201007/artisanal.bpch'

         ! Echo info
         WRITE( 6, 100 ) TRIM( FILENAME )

         ! Read data in [kg/yr]
         CALL READ_BPCH2( FILENAME, 'HG-SRCE',     1,  
     &       XTAU,      IIPAR,        JJPAR,    
     &        1,         ARRAY(:,:,1), QUIET=.TRUE. )

         ! Cast to REAL*8 and resize
         CALL TRANSFER_2D( ARRAY(:,:,1), EHg0_am(:,:) )
     
         ! Convert from Mg/yr to kg/s; multiply by 50%
         ! need more source: try 75%
         EHg0_am = Ehg0_am * 1D3 / SEC_PER_YR * 0.75d0

         !Artisanal mining is included in the 2006 Streets and 2005
         !GEIA emissions inventories. Here, we don't add it if those
         !flags are set. (eck, 10/19/11)
         IF (LGCAPEMIS) THEN
            EHg0_an = EHg0_an
         ELSEIF (LGEIA05) THEN
            EHg0_an = EHg0_an
         ELSE
         ! Add artisanal mining to anthropogenic
         EHg0_an = EHg0_an + EHg0_am
         ENDIF

         !---------------------------
         ! Anthropogenic Hg(II) emissions [kg/s]
         !---------------------------

         IF (LGCAPEMIS) THEN
            ! Use Bess' GCAP emissions

            ! Filename for anthropogenic mercury source
            FILENAME = TRIM( DATA_DIR_1x1 )               // 
     &           'mercury_201002/GEIA_Streets_Hg2.geos.1x1.YYYY'

            ! Add year to the filename
            CALL EXPAND_DATE( FILENAME, NYMD, 000000 )
            XTAU     = GET_TAU0( 1, 1, 2006 )

            ! Echo info
            WRITE( 6, 100 ) TRIM( FILENAME )

            ! Read data in [kg/yr]
            CALL READ_BPCH2( FILENAME, 'HG-SRCE', 6,
     &                       XTAU,      I1x1,     J1x1-1,
     &                       1,         ARRAY1,   QUIET=.TRUE. )

            ! Cast to REAL*8 before regridding
            GEN_1x1(:,:) = ARRAY1(:,:,1)

            ! Regrid from GENERIC 1x1 to the current grid
            INGRID  => GEN_1x1(:,:)
            CALL DO_REGRID_A2A( LLFILENAME, I1x1,    J1x1-1,
     &                          INGRID,     EHg2_an, PERAREA=1 )

            ! Convert from [kg/yr] to [kg/s]
            EHg2_an = EHg2_an / SEC_PER_YR
        
            ! Free pointer
            NULLIFY( INGRID )

         ELSEIF (LGEIA05) THEN
            ! Use GEIA 2005 emissions (eck, 10/19/11)

            ! Filename for anthropogenic mercury source
            FILENAME = TRIM( DATA_DIR_1x1 )               // 
     &           'mercury_201002/GEIA_Hg2.geos.1x1.YYYY'

            ! Add year to the filename
            CALL EXPAND_DATE( FILENAME, NYMD, 000000 )
            XTAU     = GET_TAU0( 1, 1, 2005 )

            ! Echo info
            WRITE( 6, 100 ) TRIM( FILENAME )

            ! Read data in [kg/yr]
            CALL READ_BPCH2( FILENAME, 'HG-SRCE', 6,
     &                       XTAU,      I1x1,     J1x1-1,
     &                       1,         ARRAY1,   QUIET=.TRUE. )

            ! Cast to REAL*8 before regridding
            GEN_1x1(:,:) = ARRAY1(:,:,1)

            ! Regrid from GENERIC 1x1 to the current grid
            INGRID  => GEN_1x1(:,:)
            CALL DO_REGRID_A2A( LLFILENAME, I1x1,    J1x1-1,
     &                          INGRID,     EHg2_an, PERAREA=1 )

            ! Convert from [kg/yr] to [kg/s]
            EHg2_an = EHg2_an / SEC_PER_YR

            ! For GEIA 2005, change Hg0:Hg2:HgP ratio
            ! from 50:40:10 to 86.5:9.9:3.6 for fossil fuel emissions.
            ! Fossil fuels make up ~47% of total global anthro emis
            ! Refs: Amos et al. (2011, ACPD), Y. Zhang et al. (2011, in prep)
            ! Implemented 10/19/11 by eck
            EHg2_an=0.45d0*EHg2_an

            ! Free pointer
            NULLIFY( INGRID )

         ELSE
            ! Otherwise use Selin et al. 2008 emissions

            ! Filename for anthropogenic mercury source
            FILENAME = TRIM( DATA_DIR_1x1 ) // 
     &           'mercury_200511/GEIA_Hg2.geos.1x1.YYYY'

            ! Add year to the filename
            CALL EXPAND_DATE( FILENAME, NYMD, 000000 )
            XTAU     = GET_TAU0( 1, 1, ANTHRO_Hg_YEAR )

            ! Echo info
            WRITE( 6, 100 ) TRIM( FILENAME )

            ! Read data in [kg/yr]
            CALL READ_BPCH2( FILENAME, 'HG-SRCE', 6,
     &                       XTAU,      I1x1,     J1x1-1,
     &                       1,         ARRAY1,   QUIET=.TRUE. )

            ! Cast to REAL*8 before regridding
            GEN_1x1(:,:) = ARRAY1(:,:,1)

            ! Regrid from GENERIC 1x1 to the current grid
            INGRID  => GEN_1x1(:,:)
            CALL DO_REGRID_A2A( LLFILENAME, I1x1,    J1x1-1,
     &                          INGRID,     EHg2_an, PERAREA=1 )

            ! Convert from [kg/yr] to [kg/s]; increase Pacyna by 30%
            EHg2_an = EHg2_an * 1.3d0 / SEC_PER_YR

            ! Free pointer
            NULLIFY( INGRID )

         ENDIF

         !---------------------------
         ! HgP emissions [kg/s]
         !---------------------------

         IF (LGCAPEMIS) THEN
            ! Use Bess' GCAP emissions

            ! Filename for anthropogenic mercury source
            FILENAME = TRIM( DATA_DIR_1x1)   //
     &           'mercury_201002/GEIA_Streets_HgP.geos.1x1.YYYY'

            ! Add year to the filename
            CALL EXPAND_DATE( FILENAME, NYMD, 000000 )

            ! Echo info
            WRITE( 6, 100 ) TRIM( FILENAME )

            ! Read data in [kg/yr]
            CALL READ_BPCH2( FILENAME, 'HG-SRCE', 9,
     &                       XTAU,      I1x1,     J1x1-1,
     &                       1,         ARRAY1,   QUIET=.TRUE. )

            ! Cast to REAL*8 before regridding
            GEN_1x1(:,:) = ARRAY1(:,:,1)

            ! Regrid from GENERIC 1x1 to the current grid
            INGRID  => GEN_1x1(:,:)
            CALL DO_REGRID_A2A( LLFILENAME, I1x1,    J1x1-1,
     &                          INGRID,     EHgP_an, PERAREA=1 )

            ! Convert from [kg/yr] to [kg/s]
            EHgP_an = EHgP_an / SEC_PER_YR

            ! Free pointer
            NULLIFY( INGRID )

         ELSEIF (LGEIA05) THEN
            ! Use GEIA 2005 emissions (eck, 10/19/11)

            ! Filename for anthropogenic mercury source
            FILENAME = TRIM( DATA_DIR_1x1)   //
     &           'mercury_201002/GEIA_HgP.geos.1x1.YYYY'

            ! Add year to the filename
            CALL EXPAND_DATE( FILENAME, NYMD, 000000 )

            ! Echo info
            WRITE( 6, 100 ) TRIM( FILENAME )

            ! Read data in [kg/yr]
            CALL READ_BPCH2( FILENAME, 'HG-SRCE', 9,
     &                       XTAU,      I1x1,     J1x1-1,
     &                       1,         ARRAY1,   QUIET=.TRUE. )

            ! Cast to REAL*8 before regridding
            GEN_1x1(:,:) = ARRAY1(:,:,1)

            ! Regrid from GENERIC 1x1 to the current grid
            INGRID  => GEN_1x1(:,:)
            CALL DO_REGRID_A2A( LLFILENAME, I1x1,    J1x1-1,
     &                          INGRID,     EHgP_an, PERAREA=1 )

            ! Convert from [kg/yr] to [kg/s]
            EHgP_an = EHgP_an / SEC_PER_YR

            ! For GEIA 2005, change Hg0:Hg2:HgP ratio
            ! from 50:40:10 to 86.5:9.9:3.6 for fossil fuel emissions.
            ! Fossil fuels make up ~47% of total global anthro emis
            ! Refs: Amos et al. (2011, ACPD), Y. Zhang et al. (2011, in prep)
            ! Implemented 10/19/11 by eck
            EHgP_an=0.55d0*EHgP_an

            ! Free pointer
            NULLIFY( INGRID )

         ELSE
            ! Otherwise use Selin et al. 2008 emissions

            ! Filename for anthropogenic mercury source
            FILENAME = TRIM( DATA_DIR_1x1) //
     &           'mercury_200511/GEIA_HgP.geos.1x1.YYYY'

            ! Add year to the filename
            CALL EXPAND_DATE( FILENAME, NYMD, 000000 )

            ! Echo info
            WRITE( 6, 100 ) TRIM( FILENAME )

            ! Read data in [kg/yr]
            CALL READ_BPCH2( FILENAME, 'HG-SRCE', 9,
     &                       XTAU,      I1x1,     J1x1-1,
     &                       1,         ARRAY1,   QUIET=.TRUE. )

            ! Cast to REAL*8 before regridding
            GEN_1x1(:,:) = ARRAY1(:,:,1)

            ! Regrid from GENERIC 1x1 to the current grid
            INGRID  => GEN_1x1(:,:)
            CALL DO_REGRID_A2A( LLFILENAME, I1x1,    J1x1-1,
     &                          INGRID,     EHgP_an, PERAREA=1 )

            ! Convert from [kg/yr] to [kg/s];increase Pacyna by 30%
            EHgP_an = EHgP_an * 1.3D0 / SEC_PER_YR

            ! Free pointer
            NULLIFY( INGRID )

         ENDIF

      ENDIF
 
      ! Assume HgP is part of the pool of Hg2 
      ! available for gas-particle partitioning
      EHg2_an = EHg2_an + EHgP_an
      EHgP_an = 0d0 * EHgP_an

      !=================================================================
      ! Read distribution of Hg on land.
      ! See explanation in Selin et al., GBC 2008 or the GEOS-Chem Hg
      !  wiki.
      !
      ! EHg0_dist is a dimensionless spatial scaling factor that distributes
      ! the soil Hg content according to deposition. It also accounts for
      ! the global anthropogenic enrichment in the present day.
      !=================================================================

      ! Use "generic" year 1985 for TAU values
      XTAU     = GET_TAU0( 1, 1, 1985 )

      IF ( LPREINDHG ) THEN
         !-----------------------
         ! Pre-industrial 
         !-----------------------
#if   defined( GEOS_5 ) && defined( GRID4x5 )
        FILENAME= TRIM( DATA_DIR )            //
     &            'mercury_201110/soilHg.102311.preind.bpch'
#elif defined( GEOS_5 ) && defined( GRID2x25 )
        FILENAME= TRIM( DATA_DIR )            //
     &            'mercury_201110/soilHg.102311.preind.bpch'
#elif defined( MERRA ) && defined( GRID4x5 )
        FILENAME= TRIM( DATA_DIR )            //
     &            'mercury_201110/soilHg.merra.preind.bpch'
#else
      CALL ERROR_STOP('You need to make a Hg soil file. ')
#endif

      ELSE
         !-----------------------
         ! Present day
         !-----------------------
#if   defined( GEOS_5 ) && defined( GRID4x5 )
        FILENAME= TRIM( DATA_DIR )            //
     &            'mercury_201110/soilHg.102311.presentday.bpch'
#elif defined( GEOS_5 ) && defined( GRID2x25 )
        FILENAME= TRIM( DATA_DIR )            //
     &            'mercury_201110/soilHg.102311.presentday.bpch'
#elif defined( MERRA ) && defined( GRID4x5 )
        FILENAME= TRIM( DATA_DIR )            //
     &            'mercury_201110/soilHg.merra.presentday.bpch'
#else
      CALL ERROR_STOP('You need to make a Hg soil file. ')
#endif

      ENDIF

      ! Echo info
      WRITE( 6, 100 ) TRIM( FILENAME )

      ! Read data in [kg/yr]
      CALL READ_BPCH2( FILENAME, 'HG-SRCE',     4,  
     &                 XTAU,      IIPAR,        JJPAR,    
     &                 1,         ARRAY(:,:,1), QUIET=.TRUE. )

      ! Cast to REAL*8 and resize
      CALL TRANSFER_2D( ARRAY(:,:,1), EHg0_dist(:,:) )

      !=================================================================
      ! Read annual emissions of Hg(0) from natural land sources [kg/s]
      !=================================================================

      ! Use "generic" year 1985 for TAU values
      XTAU     = GET_TAU0( 1, 1, 1985 )

      ! Filename for natural land-source mercury
      FILENAME= TRIM( DATA_DIR )            //
     &     'mercury_201007/newnatural.bpch'

      ! Echo info
      WRITE( 6, 100 ) TRIM( FILENAME )

      ! TAU value corresponding to the data 
      XTAU   = GET_TAU0( 1, 1, 1995 )
      
      ! Read data in [kg/yr]
      CALL READ_BPCH2( FILENAME, 'HG-SRCE',     5, 
     &                 XTAU,      IIPAR,        JJPAR,    
     &                 1,         ARRAY(:,:,1), QUIET=.TRUE. )

      ! Cast to REAL*8 and resize
      CALL TRANSFER_2D( ARRAY(:,:,1), EHg0_nt )

      ! Convert from [kg/yr] to [kg/s]
      EHg0_nt = EHg0_nt / SEC_PER_YR 

      !Recommended for v8-03-02 and later: 
      !Reduce Geogenic emissions by 50%
      !Ref: Soerensen et al., 2010
      EHg0_nt = EHg0_nt *0.5d0

      !=================================================================
      ! Read offline ocean Hg0 emissions (if LDYNOCEAN = .FALSE.)
      ! 
      ! Note from cdh (8/13/2008):
      ! These static ocean emissions were developed for GEOS-Chem before
      ! v7-04-06 and may not be useful any longer. They were
      ! tuned with older, lower emissions, and with different chemical rate
      ! constants.
      !=================================================================
      IF ( .not. LDYNOCEAN ) THEN

         ! File name
         FILENAME = TRIM( DATA_DIR_1x1 )            // 
     &              'mercury_200511/Hg_ocean.geos.' // GET_RES_EXT()    

         ! Echo info
         WRITE( 6, 100 ) TRIM( FILENAME )

         ! TAU value corresponding to the data 
         XTAU  = GET_TAU0( 1, 1, 1985 )

         ! Read data in [kg/yr]
         CALL READ_BPCH2( FILENAME, 'HG-SRCE',     3,  
     &                    XTAU,      IIPAR,        JJPAR,    
     &                    1,         ARRAY(:,:,1), QUIET=.TRUE. )

         ! Cast to REAL*8 and resize
         CALL TRANSFER_2D( ARRAY(:,:,1), EHg0_oc(:,:,1) )

         ! Convert from [kg/yr] to [kg/s]
         EHg0_oc = EHg0_oc / SEC_PER_YR

      ENDIF

      !=================================================================
      ! Print totals to the screen in [Gg/yr]
      !=================================================================
      WRITE( 6, '(a)' ) REPEAT( '=', 79 )
      WRITE( 6, 110   )
      WRITE( 6, '(a)' )
      WRITE( 6, 111   ) SUM( EHg0_an ) * SEC_PER_YR * 1d-6
      WRITE( 6, 113   ) SUM( EHg0_ln ) * SEC_PER_YR * 1d-6
      WRITE( 6, 114   ) SUM( EHg0_nt ) * SEC_PER_YR * 1d-6

      ! Only write ocean total if we are doing offline ocean
      IF ( .not. LDYNOCEAN ) THEN
         WRITE( 6, 117   ) SUM( EHg0_oc ) * SEC_PER_YR * 1d-6 
      ENDIF

      WRITE( 6, 115   ) SUM( EHg2_an ) * SEC_PER_YR * 1d-6
      WRITE( 6, 116   ) SUM( EHgP_an ) * SEC_PER_YR * 1d-6
      WRITE( 6, '(a)' ) REPEAT( '=', 79 )

      ! FORMAT strings
 110  FORMAT( 'M E R C U R Y   E M I S S I O N S' )
 111  FORMAT( 'Total Anthro     Hg(0)  : ', f7.3, ' [Gg/yr]' )
 113  FORMAT( 'Total Re-Emitted Hg(0)  : ', f7.3, ' [Gg/yr]' )
 114  FORMAT( 'Total Natural    Hg(0)  : ', f7.3, ' [Gg/yr]' )
 115  FORMAT( 'Total Anthro     Hg(II) : ', f7.3, ' [Gg/yr]' )
 116  FORMAT( 'Total Anthro     HgP    : ', f7.3, ' [Gg/yr]' )
 117  FORMAT( 'Total Ocean      Hg(0)  : ', f7.3, ' [Gg/yr]' )
 
      ! Return to calling program
      END SUBROUTINE MERCURY_READYR

!------------------------------------------------------------------------------

      FUNCTION GET_LWC( T ) RESULT( LWC )
!
!******************************************************************************
!  Function GET_LWC returns the cloud liquid water content [m3 H2O/m3 air] at 
!  a GEOS-CHEM grid box as a function of temperature. (rjp, bmy, 10/31/02, 
!  12/7/04)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) T (REAL*8) : Temperature value at a GEOS-CHEM grid box [K]
!
!  NOTES:
!  18 Jan 2011 - R. Yantosca - Updated comments
!******************************************************************************
!
      ! Arguments
      REAL*8, INTENT(IN) :: T

      ! Function value
      REAL*8             :: LWC

      !=================================================================
      ! GET_LWC begins here!
      !=================================================================

      ! Compute Liquid water content in [g/m3]
      IF (T>293d0) THEN
         LWC=0.2d0
       ELSE IF ( T >= 280.d0 .AND. T <= 293.d0 ) THEN
         LWC = 0.32d0 - 0.0060d0 * ( T - 273.D0 ) 
 
      ELSE IF ( T >= 248.d0 .AND. T < 280.d0 ) THEN
         LWC = 0.23d0 + 0.0065d0 * ( T - 273.D0 )

      ELSE IF ( T < 248.d0 ) THEN
         LWC = 0.07d0

      ENDIF

      ! Convert from [g/m3] to [m3/m3]
      ! Units: [g H2O/m3 air] * [1 kg H2O/1000g H2O] * [m3 H2O/1000kg H2O]
      LWC = LWC * 1.D-6         

      ! Return to calling program
      END FUNCTION GET_LWC

!------------------------------------------------------------------------------

#if defined( DEVEL )
      FUNCTION GET_VCLDF( I, J, L, LOCAL_MET ) RESULT( VCLDF )
#else
      FUNCTION GET_VCLDF( I, J, L ) RESULT( VCLDF )
#endif
! 
!******************************************************************************
!  Subroutine GET_VCLDF computes the volume cloud fraction for Hg0 and Hg2
!  chemistry. (eck, bmy, 12/6/04)
!
!  Arguments as Input:
!  ============================================================================
!  (1-3) I, J, L (INTEGER) : GEOS-CHEM lon, lat, alt indices 
!
!  References:
!  ============================================================================
!  (1  ) Sundqvist et al. [1989]
!
!  NOTES:
!  (1  ) Copied from "sulfate_mod.f" but was made into a function since we are
!         already looping over (I,J,L) in CHEM_Hg0_Hg2 (eck, bmy, 12/6/04)
!******************************************************************************
!
      ! References to F90 modules 
      USE DAO_MOD,      ONLY : RH
      USE PRESSURE_MOD, ONLY : GET_PCENTER, GET_PEDGE
#if defined( DEVEL )
      USE GC_TYPE_MOD,  ONLY : GC_MET_LOCAL
#endif
      
      ! Arguments
      INTEGER, INTENT(IN)  :: I, J, L
#if defined( DEVEL )
      TYPE(GC_MET_LOCAL), INTENT(IN) :: LOCAL_MET  ! Obj w/ met fields
#endif

      ! Local variables
      REAL*8               :: PRES, PSFC, RH2, R0, B0
      REAL*8,  PARAMETER   :: ZRT = 0.60d0, ZRS = 0.99d0
	
      ! Function value
      REAL*8               :: VCLDF
	
      !=================================================================
      ! GET_VCLDF begins here!
      !=================================================================
     
      ! Surface pressure
      PSFC = GET_PEDGE(I,J,1)         

      ! Pressure at center of grid box L
      PRES = GET_PCENTER(I,J,L)
        
      ! RH (from "dao_mod.f") is relative humidity [%]
      ! Convert to fraction and store in RH2
#if defined( DEVEL )
      RH2  = LOCAL_MET%RH(I,J,L) * 1.0d-2
#else
      RH2  = RH(I,J,L) * 1.0d-2
#endif
        
      ! Terms from Sundqvist ???
      R0   = ZRT + ( ZRS - ZRT ) * EXP( 1d0 - ( PSFC / PRES )**2.5 )
      B0   = ( RH2 - R0 ) / ( 1d0 - R0 )
	
      ! Force B0 into the range 0-1
      IF ( RH2 < R0  ) B0 = 0d0
      IF ( B0  > 1d0 ) B0 = 1d0

      ! Volume cloud fraction
      VCLDF = 1d0 - SQRT( 1d0 - B0 )

      ! Return to calling program
      END FUNCTION GET_VCLDF

!------------------------------------------------------------------------------

#if defined( DEVEL )
      FUNCTION GET_O3( I, J, L, LOCAL_MET ) RESULT( O3_MOLEC_CM3 )
#else
      FUNCTION GET_O3( I, J, L ) RESULT( O3_MOLEC_CM3 )
#endif
!
!******************************************************************************
!  Function GET_O3 returns monthly mean O3 for offline sulfate aerosol
!  simulations. (bmy, 12/16/02)
!
!  Arguments as Input:
!  ============================================================================
!  (1-3) I, J, L   (INTEGER) : Grid box indices for lon, lat, vertical level
!
!  NOTES:
!  (1 ) We assume SETTRACE has been called to define IDO3. (bmy, 12/16/02)
!  (2 ) Now reference inquiry functions from "tracer_mod.f" (bmy, 7/20/04)
!******************************************************************************
!
      ! References to F90 modules
      USE DAO_MOD      , ONLY : AD
      USE GLOBAL_O3_MOD, ONLY : O3
#if defined( DEVEL )
      USE GC_TYPE_MOD,   ONLY : GC_MET_LOCAL
#endif

      USE CMN_SIZE_MOD  ! Size parameters

      ! Arguments
      INTEGER, INTENT(IN) :: I, J, L
#if defined( DEVEL )
      TYPE(GC_MET_LOCAL), INTENT(IN) :: LOCAL_MET  ! Obj w/ met fields
#endif

      ! Local variables
      REAL*8              :: O3_MOLEC_CM3

      ! External functions
      REAL*8, EXTERNAL    :: BOXVL
      
      !=================================================================
      ! GET_O3 begins here!
      !=================================================================

      ! Get ozone [v/v] for this gridbox & month
      ! and convert to [molec/cm3] (eck, 12/2/04)
#if defined( DEVEL )
      O3_MOLEC_CM3 = O3(I,J,L)           * ( 6.022d23 / 28.97d-3 ) * 
     &               LOCAL_MET%AD(I,J,L) / BOXVL(I,J,L,LOCAL_MET)
#else
      O3_MOLEC_CM3 = O3(I,J,L) * ( 6.022d23 / 28.97d-3 ) * 
     &               AD(I,J,L)  /  BOXVL(I,J,L)
#endif

      ! Return to calling program
      END FUNCTION GET_O3

!------------------------------------------------------------------------------

      FUNCTION GET_OH( I, J, L ) RESULT( OH_MOLEC_CM3 )
!
!******************************************************************************
!  Function GET_OH returns monthly mean OH and imposes a diurnal variation.
!  (eck, bmy, 12/7/04)
!
!  Arguments as Input:
!  ============================================================================
!  (1-3) I, J, L (INTEGER) : Grid box indices for lon, lat, vertical level
!
!  NOTES:
!******************************************************************************
!
      ! References to F90 modules
      USE DAO_MOD      , ONLY : SUNCOS
      USE GLOBAL_OH_MOD, ONLY : OH
      USE TIME_MOD     , ONLY : GET_TS_CHEM

      USE CMN_SIZE_MOD  ! Size parameters

      ! Arguments
      INTEGER, INTENT(IN) :: I, J, L

      ! Local variables
      INTEGER             :: JLOOP
      REAL*8              :: OH_MOLEC_CM3
 
      !=================================================================
      ! GET_OH begins here!
      !=================================================================

      ! 1-D grid box index for SUNCOS
      JLOOP = ( (J-1) * IIPAR ) + I

      ! Test for sunlight...
      IF ( SUNCOS(JLOOP) > 0d0 .and. TCOSZ(I,J) > 0d0 ) THEN

         ! Impose a diurnal variation on OH during the day
         OH_MOLEC_CM3 = OH(I,J,L)                      *           
     &                  ( SUNCOS(JLOOP) / TCOSZ(I,J) ) *
     &                  ( 1440d0        / GET_TS_CHEM() )

         ! Make sure OH is not negative
         OH_MOLEC_CM3 = MAX( OH_MOLEC_CM3, 0d0 )
               
      ELSE

         ! At night, OH goes to zero
         OH_MOLEC_CM3 = 0d0

      ENDIF

      ! Return to calling program
      END FUNCTION GET_OH

!------------------------------------------------------------------------------

#if defined( DEVEL )
      FUNCTION GET_BR( I, J, L, LOCAL_MET, BRO_MOLEC_CM3 )
     &   RESULT( BR_MOLEC_CM3 )
#else
      FUNCTION GET_BR( I, J, L, BRO_MOLEC_CM3 ) RESULT( BR_MOLEC_CM3 )
#endif
!
!******************************************************************************
!  Function GET_BR returns instantaneous Br concentration calculated from
!  the monthly mean and an imposed diurnal variation.
!  (cdh, 07/06/06)
!
!  Arguments as Input:
!  ============================================================================
!  (1-3) I, J, L (INTEGER) : Grid box indices for lon, lat, vertical level
!
!  Arguments as Output:
!  ============================================================================
!  (1-3) BrO (REAL*8)      : Concentration of BrO
!
!  NOTES:
!  01 Mar 2012 - R. Yantosca - Now use GET_YMID(I,J,L) from grid_mod.F90
!  01 Mar 2012 - R. Yantosca - Now use GET_LOCALTIME(I,J,L) from time_mod.F90
!******************************************************************************
!
      ! References to F90 modules
      USE DAO_MOD       , ONLY : SUNCOS,    AD,       T
      USE DAO_MOD       , ONLY : LWI,       GET_OBK
      USE DAO_MOD       , ONLY : IS_WATER
      USE GLOBAL_BR_MOD , ONLY : BR_MERGE,  BR_TROP,  BR_STRAT, J_BRO
      USE GLOBAL_BR_MOD , ONLY : BRO_MERGE, BRO_TROP, BRO_STRAT
      USE TIME_MOD      , ONLY : GET_TS_CHEM, GET_LOCALTIME 
      USE TIME_MOD      , ONLY : ITS_A_NEW_DAY
      USE TROPOPAUSE_MOD, ONLY : ITS_IN_THE_STRAT
      USE PBL_MIX_MOD   , ONLY : GET_FRAC_UNDER_PBLTOP
      USE GRID_MOD      , ONLY : GET_YMID
      USE TIME_MOD      , ONLY : GET_MONTH
      USE PRESSURE_MOD  , ONLY : GET_PCENTER
      USE ERROR_MOD     , ONLY : SAFE_DIV

#if defined( DEVEL )
      USE GC_TYPE_MOD,    ONLY : GC_MET_LOCAL
#endif

      USE CMN_SIZE_MOD  ! Size parameters

      ! Arguments
      INTEGER, INTENT(IN) :: I, J, L
      REAL*8, INTENT(OUT) :: BrO_MOLEC_CM3 
#if defined( DEVEL )
      TYPE(GC_MET_LOCAL), INTENT(IN) :: LOCAL_MET  ! Obj w/ met fields
#endif

      ! Local variables
      INTEGER             :: JLOOP
      REAL*8              :: BR_MOLEC_CM3, BR_PPTV, BR_MBL, BR_FAC
      REAL*8              :: LOCALTIME, HOUR, BRO_PPTV, FPBL

      ! External functions
      REAL*8, EXTERNAL    :: BOXVL           

      ! Constant concentration of BrO assumed for the marine boundary
      ! layer, based on finding from the Holmes et al.(2009) box model
      ! Use 0.5pptv for 24-hr mean, based on 1pptv in box model daytime
      REAL*8, PARAMETER   :: BRO_MBL=0.5D0 ! pptv

      ! Constant concentration of BrO assumed for the Polar boundary layer
      ! during springtime bromine explosion. Based on DOAS observations
      ! at Alert (Simpson et al. ACP 2007). This concentration applies
      ! only under stable atmospheric conditions.
      !REAL*8, PARAMETER   :: BRO_POLAR=10D0
      ! Reduce concentration of BrO in Arctic during depletion events to 5 pptv
      ! Ref low end of uncertainty range in Holmes et al. 2010
      ! Ref low [BrO] observed Neuman et al 2010
      REAL*8, PARAMETER   :: BRO_POLAR=5D0

      ! Parameters for calculating Br/BrO photostationary state
      ! BrO J value, /s 
      ! Rate coefficient BrO + NO -> Br + NO2, cm3/molec/s
      REAL*8, PARAMETER   :: K_BRO_NO = 2.1D-11
      ! Rate coefficient Br + O3 -> BrO + O2, cm3/molec/s
      REAL*8, PARAMETER   :: K_BR_O3  = 1.2D-12
      ! Concentration of NO, based on 10pptv, molec/cm3
      REAL*8, PARAMETER   :: C_NO     = 2.5D8

      !=================================================================
      ! GET_BR begins here!
      !=================================================================

      !----------------------------------------------------------------
      ! Get monthly-mean Br from bromocarbon and halon sources 
      !----------------------------------------------------------------
    
      BR_PPTV  = BR_MERGE(I,J,L)
      BRO_PPTV = BRO_MERGE(I,J,L)

      ! Multiple stratospheric concentrations by Factor, if necessary
#if defined( DEVEL )
      IF ( (ITS_IN_THE_STRAT( I, J, L, LOCAL_MET )) .AND.
     &     (ABS(STRAT_BR_FACTOR - 1d0) > 1d-2) )    THEN
#else
      IF ( (ITS_IN_THE_STRAT(I,J,L)) .AND.
     &     (ABS(STRAT_BR_FACTOR - 1d0) > 1d-2) ) THEN
#endif

         BR_PPTV  = BR_PPTV  * STRAT_BR_FACTOR
         BRO_PPTV = BRO_PPTV * STRAT_BR_FACTOR
         
      ENDIF 

      !----------------------------------------------------------------
      ! Add Br in the MBL
      !----------------------------------------------------------------

      ! Fraction of grid box mass which is within the mixed layer
      FPBL = GET_FRAC_UNDER_PBLTOP(I,J,L)

#if defined( DEVEL )
      IF ( IS_WATER( I, J, LOCAL_MET ) .AND. (FPBL > 0d0) ) THEN
#else
      IF ( IS_WATER(I,J) .AND. (FPBL > 0d0) ) THEN
#endif

         ! Convert BrO concentration to corresponding Br based on
         ! photochemical steady state, pptv
         ! Platt & Janssen, Faraday Discussions (1995)
#if defined( DEVEL )
         BR_MBL = BRO_MBL *( J_BRO(I,J,L) + K_BRO_NO * C_NO ) / 
     &                     ( K_BR_O3      * GET_O3(I,J,L,LOCAL_MET) ) *
     &                       FPBL
#else
         BR_MBL = BRO_MBL *( J_BRO(I,J,L) + K_BRO_NO * C_NO ) / 
     &                     ( K_BR_O3 * GET_O3(I,J,L) ) *
     &                     FPBL
#endif

         ! MBL concentration is the greater of the TOMCAT Br or 1pptv BrO
         BR_PPTV  = MAX( BR_PPTV,  BR_MBL )
         BrO_PPTV = MAX( BRO_PPTV, BRO_MBL*FPBL )


      ENDIF            

      !----------------------------------------------------------------
      ! Add Br in polar boundary layer
      ! 
      ! Bromine in the polar boundary layer requires cold temperatures
      ! (below 268K), sunlight and contact between open water and sea ice
      ! Assuming these conditions are met, we assume a fixed 5-10pptv BrO 
      ! At the moment, we only have Br north of 60N, i.e. not Hudson's Bay
      !----------------------------------------------------------------

      ! 1-D grid box index for SUNCOS
      JLOOP = ( (J-1) * IIPAR ) + I

#if defined( DEVEL )
      IF ( ( LPOLARBR ) .AND. ( FPBL >  0d0   ) .AND.
     &     ( SUNCOS(JLOOP)           >  0D0   ) .AND. 
     &     ( LOCAL_MET%T(I,J,1)      <= 268D0 ) .AND.
     &     ( LOCAL_MET%LWI(I,J)      == 2     ) .AND. 
     &       GET_OBK(I,J,LOCAL_MET)  >  0D0   ) THEN
#else
      IF ( (LPOLARBR)         .AND. (FPBL > 0d0) .AND.
     &     (SUNCOS(JLOOP) > 0D0) .AND. (T(I,J,1) <= 268D0) .AND.
     &     ( LWI(I,J)==2 )  .AND. GET_OBK(I,J) > 0D0 ) THEN
#endif

         IF ( ( (GET_YMID(I,J,L) > 60D0) .AND. 
     &          (GET_MONTH() >= 3) .AND. (GET_MONTH() <= 5)  ) .OR.
     &        ( (GET_YMID(I,J,L) < -50D0) .AND. 
     &          (GET_MONTH() >= 9) .AND. (GET_MONTH() <= 11)  ) ) THEN

            ! Br concentration due to bromine explosion, pptv
            ! Assume [O3] is 2ppb during event
            BR_MBL = BRO_POLAR *( J_BRO(I,J,L) + K_BRO_NO * C_NO ) / 
     &                          ( K_BR_O3 * 5D10 ) * FPBL
        
            BR_PPTV  = BR_PPTV  + BR_MBL
            BrO_PPTV = BrO_PPTV + BRO_POLAR * FPBL

         ENDIF

      ENDIF

      !----------------------------------------------------------------
      ! Impose a diurnal cycle
      !
      ! This formula comes from the Holmes et al (2009) MBL box model
      ! which fit RGM observations for Okinawa and a Pacific cruise
      !----------------------------------------------------------------
            
      ! 1-D grid box index for SUNCOS
      JLOOP = ( (J-1) * IIPAR ) + I
      
      ! Test for sunlight...
      IF ( (SUNCOS(JLOOP) > 0d0) .and. (TTDAY(I,J) > 0d0) ) THEN
      
         ! Use a constant function if daylight is < 2 hours or > 22hours
         IF ( ( TTDAY(I,J) < 120d0  ) .OR. 
     &        ( TTDAY(I,J) > 1320d0 ) ) THEN
            
            BR_FAC = SAFE_DIV( 1440d0, TTDAY(I,J), 0D0 )  
            
         ELSE

            ! Local time: 0-23.999
            LOCALTIME = GET_LOCALTIME( I, J, L )
      
            ! Interpolate the real hour to lie within an equinoctal day
            ! i.e. between 6-18
            HOUR = 12D0 + ( LOCALTIME - 12D0 ) *
     &           720D0 / ( TTDAY(I,J) ) 

            ! Multiplicative factor for the diurnal cycle 
            BR_FAC = ( 1D6 - 1D2 * ( HOUR - 6D0 ) ** 4D0 + 
     &           ( 18D0 - HOUR ) * 1D6 / 12D0 ) / 2D6

            ! Normalize the multiplicative factor to have a 24-h mean of 1
            BR_FAC = SAFE_DIV( BR_FAC, ( 4D-4 * TTDAY(I,J) ), 0D0 )

         ENDIF


         ! Make sure that diurnal scaling factor is non-negative
         BR_FAC = MAX( BR_FAC, 0D0 )

         ! The instantaneous concentration is the 24-h mean times
         ! the time-varying factor from the diurnal cycle
         BR_PPTV = BR_PPTV * BR_FAC


      ELSE

         ! At night, BR goes to zero
         BR_PPTV = 0D0

      ENDIF


      ! Convert pptv mixing ratio -> molec/cm3
#if defined( DEVEL )
      BR_MOLEC_CM3  = BR_PPTV             * 1D-12      *
     &                ( 6.022D23          / 28.97D-3 ) * 
     &                LOCAL_MET%AD(I,J,L) / BOXVL(I,J,L,LOCAL_MET)
      BRO_MOLEC_CM3 = BRO_PPTV            * 1D-12      *
     &                ( 6.022D23          / 28.97D-3 ) * 
     &                LOCAL_MET%AD(I,J,L) / BOXVL(I,J,L,LOCAL_MET)
#else
      BR_MOLEC_CM3  = BR_PPTV * 1D-12 * ( 6.022D23 / 28.97D-3 ) * 
     &                AD(I,J,L)  /  BOXVL(I,J,L)
      BRO_MOLEC_CM3 = BRO_PPTV * 1D-12 * ( 6.022D23 / 28.97D-3 ) * 
     &                AD(I,J,L)  /  BOXVL(I,J,L)
#endif

      ! Return to calling program
      END FUNCTION GET_BR

!------------------------------------------------------------------------------

      SUBROUTINE OHNO3TIME
!
!******************************************************************************
!  Subroutine OHNO3TIME computes the sum of cosine of the solar zenith
!  angle over a 24 hour day, as well as the total length of daylight. 
!  This is needed to scale the offline OH and NO3 concentrations.
!  (rjp, bmy, 12/16/02, 12/8/04)
!
!  NOTES:
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
!******************************************************************************
!
      ! References to F90 modules
      USE GRID_MOD, ONLY : GET_XMID,    GET_YMID_R
      USE TIME_MOD, ONLY : GET_NHMSb,   GET_ELAPSED_SEC
      USE TIME_MOD, ONLY : GET_TS_CHEM, GET_DAY_OF_YEAR, GET_GMT

      USE CMN_SIZE_MOD  ! Size parameters
      USE CMN_GCTM_MOD  ! Physical constants

      ! Local variables
      LOGICAL, SAVE       :: FIRST = .TRUE.
      INTEGER             :: I, IJLOOP, J, L, N, NT, NDYSTEP
      REAL*8              :: A0, A1, A2, A3, B1, B2, B3
      REAL*8              :: LHR0, R, AHR, DEC, TIMLOC, YMID_R
      REAL*8              :: SUNTMP(MAXIJ)
      
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
         TTDAY(:,:) = 0d0
         TCOSZ(:,:) = 0d0
         COSZM(:,:) = 0d0

         ! NDYSTEP is # of chemistry time steps in this day
         NDYSTEP = ( 24 - INT( GET_GMT() ) ) * 60 / GET_TS_CHEM()         

         ! NT is the elapsed time [s] since the beginning of the run
         NT = GET_ELAPSED_SEC()

         ! Loop forward through NDYSTEP "fake" timesteps for this day 
         DO N = 1, NDYSTEP
            
            ! Zero SUNTMP array
            SUNTMP(:) = 0d0

            ! Loop over surface grid boxes
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, YMID_R, IJLOOP, TIMLOC, AHR )
            DO J = 1, JJPAR
            DO I = 1, IIPAR

               ! Grid box latitude center [radians]
               YMID_R = GET_YMID_R( I, J, 1 )

               ! Increment IJLOOP
               IJLOOP = ( (J-1) * IIPAR ) + I
               TIMLOC = real(LHR0) + real(NT)/3600.0 
     &                + GET_XMID(I,J,1)/15.0
         
               DO WHILE (TIMLOC .lt. 0)
                  TIMLOC = TIMLOC + 24.0
               ENDDO

               DO WHILE (TIMLOC .gt. 24.0)
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
               SUNTMP(IJLOOP) = sin(YMID_R) * sin(DEC) +
     &                          cos(YMID_R) * cos(DEC) * cos(AHR)

               ! TCOSZ is the sum of SUNTMP at location (I,J)
               ! Do not include negative values of SUNTMP
               TCOSZ(I,J) = TCOSZ(I,J) + MAX( SUNTMP(IJLOOP), 0d0 )

               ! COSZM is the peak value of SUMTMP during a day at (I,J)
               ! (rjp, bmy, 3/30/04)
               COSZM(I,J) = MAX( COSZM(I,J), SUNTMP(IJLOOP) )

               ! TTDAY is the total daylight time at location (I,J)
               IF ( SUNTMP(IJLOOP) > 0d0 ) THEN
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

!------------------------------------------------------------------------------

      SUBROUTINE DEFINE_TAGGED_Hg
!
!******************************************************************************
!  Subroutine DEFINE_TAGGED_Hg defines the tagged tracer numbers for
!  anthropogenic (by geographic region) Hg0, Hg2, and HgP.  The position of 
!  Hg2 and HgP in the DEPSAV array is also computed for future use.  This 
!  routine only has to be called once at the start of the simulation. 
!  (eck, cdh, bmy, 12/15/04, 4/6/06)
!
!  NOTES:
!  (1 ) Now only define AN_Hg0, AN_Hg2, AN_HgP.  Now use ID_Hg0, ID_Hg2, and
!        ID_HgP index arrays from "tracerid_mod.f". (eck, cdh, bmy, 4/6/06)
!  01 Mar 2012 - R. Yantosca - Now use GET_XMID(I,J,L) from grid_mod.F90
!  01 Mar 2012 - R. Yantosca - Now use GET_YMID(I,J,L) from grid_mod.F90
!******************************************************************************
!
      ! References to F90 modules
      USE GRID_MOD    , ONLY : GET_XMID  , GET_YMID
      USE ERROR_MOD   , ONLY : ERROR_STOP
      USE TRACER_MOD  , ONLY : N_TRACERS
      USE TRACERID_MOD, ONLY : ID_Hg0    , ID_Hg2  , ID_HgP
      USE TRACERID_MOD, ONLY : ID_Hg_tot , ID_Hg_na, ID_Hg_eu, ID_Hg_as
      USE TRACERID_MOD, ONLY : ID_Hg_rw  , ID_Hg_oc, ID_Hg_ln, ID_Hg_nt

      USE CMN_SIZE_MOD     ! Size parameters

      ! Local variables
      INTEGER               :: I, J
      REAL*8                :: X, Y
      CHARACTER(LEN=255)    :: LOCATION

      !=================================================================
      ! DEFINE_TAGGED_Hg begins here!
      !=================================================================

      ! Location string for ERROR_STOP
      LOCATION = 'DEFINE_TAGGED_Hg ("mercury_mod.f")'

      !---------------------------------
      ! Anthropogenic tracer indices
      !---------------------------------

      ! Loop over latitudes
!$OMP PARALLEL DO 
!$OMP+DEFAULT( SHARED ) 
!$OMP+PRIVATE( I, J, X, Y )
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Grid box longitude and latitude [degrees]
         X = GET_XMID( I, J, 1 )
         Y = GET_YMID( I, J, 1 )
         
         ! North American Anthro Hg 
         IF ( ( X >= -172.5 .and. X < -17.5 )   .and. 
     &        ( Y >=   24.0 .and. Y <  88.0 ) ) THEN
            AN_Hg0(I,J) = ID_Hg0(ID_Hg_na)
            AN_Hg2(I,J) = ID_Hg2(ID_Hg_na)
            AN_HgP(I,J) = ID_HgP(ID_Hg_na)
         
         ! European Anthro Hg (1st sub-box)
         ELSE IF ( ( X >= -17.5 .and. X < 72.5 )  .and. 
     &             ( Y >=  36.0 .and. Y < 45.0 ) ) THEN
            AN_Hg0(I,J) = ID_Hg0(ID_Hg_eu)
            AN_Hg2(I,J) = ID_Hg2(ID_Hg_eu)
            AN_HgP(I,J) = ID_HgP(ID_Hg_eu)

         ! European Anthro Hg (2nd sub-box)
         ELSE IF ( ( X >= -17.5 .and. X < 172.5 )  .and. 
     &             ( Y >=  45.0 .and. Y <  88.0 ) ) THEN
            AN_Hg0(I,J) = ID_Hg0(ID_Hg_eu)
            AN_Hg2(I,J) = ID_Hg2(ID_Hg_eu)
            AN_HgP(I,J) = ID_HgP(ID_Hg_eu)

         ! Asian Anthro Hg 
         ELSE IF ( ( X >= 70.0 .and. X < 152.5 )  .and.
     &             ( Y >=  8.0 .and. Y <  45.0 ) ) THEN
            AN_Hg0(I,J) = ID_Hg0(ID_Hg_as)
            AN_Hg2(I,J) = ID_Hg2(ID_Hg_as)
            AN_HgP(I,J) = ID_HgP(ID_Hg_as)

         ! Rest-of-world Anthro Hg 
         ELSE
            AN_Hg0(I,J) = ID_Hg0(ID_Hg_rw)
            AN_Hg2(I,J) = ID_Hg2(ID_Hg_rw)
            AN_HgP(I,J) = ID_HgP(ID_Hg_rw)

         ENDIF
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      !=================================================================
      ! Error check tracers: make sure they are not zero
      ! since this can cause array-out-of-bounds errors
      !=================================================================

      ! Tagged Hg0
      IF ( ANY( AN_Hg0 == 0 ) ) THEN
         CALL ERROR_STOP( 'AN_Hg0 tracers are undefined!', LOCATION )   
      ENDIF
      
      ! Tagged Hg2
      IF ( ANY( AN_Hg2 == 0 ) ) THEN
         CALL ERROR_STOP( 'AN_Hg2 tracers are undefined!', LOCATION )   
      ENDIF

      ! Tagged Hg2
      IF ( ANY( AN_HgP == 0 ) ) THEN
         CALL ERROR_STOP( 'AN_HgP tracer are undefined!', LOCATION )   
      ENDIF

      !---------------------------------
      ! Error check # of tracers
      !---------------------------------
      IF ( N_TRACERS < 21 ) THEN 
         CALL ERROR_STOP( 'Too few Hg tagged tracers!',
     &                    'DEFINE_TAGGED_Hg ("mercury_mod.f")' )
      ENDIF

      ! Return to calling program
      END SUBROUTINE DEFINE_TAGGED_Hg

!------------------------------------------------------------------------------

#if defined( DEVEL )
      SUBROUTINE CALC_HG2_SEASALT_LOSSRATE( LOCAL_MET )
#else
      SUBROUTINE CALC_HG2_SEASALT_LOSSRATE
#endif
!
!******************************************************************************
!  Subroutine HG2_SEASALT_LOSSRATE calculates the loss rate of RGM (/s) by 
!  uptake of RGM into sea salt aerosol for each model grid
!  
!  The formula used here is a least-squares fit to the full-physics model of
!  sea-salt aerosol emissions, hydroscopic growth, mass-transport limited
!  uptake of Hg(II), and aerosol deposition presented by Holmes et al. (2009)
!  See Holmes et al. 2010 for evaluation of this parameterization. 
!  (cdh, 11/25/09)
!  
!  Return value is a loss frequency (/s)
!
!  NOTES:
!******************************************************************************
!

      ! References to F90 modules
      USE DAO_MOD    , ONLY : IS_WATER, RH
      USE DAO_MOD    , ONLY : U => U10M, V => V10M
      USE PBL_MIX_MOD, ONLY : GET_PBL_TOP_M
      USE FILE_MOD   , ONLY : IU_FILE
      USE RNPBBE_MOD , ONLY : SLQ

#if defined( DEVEL )
      USE GC_TYPE_MOD, ONLY : GC_MET_LOCAL
#endif

      USE CMN_SIZE_MOD     ! Size parameters

      !Arguments
#if defined( DEVEL )
      TYPE(GC_MET_LOCAL), INTENT(IN) :: LOCAL_MET  ! Obj w/ met fields
#endif

      ! Local variables
      REAL*8                :: U10M, S
      REAL*8                :: LOSS_FREQ
      INTEGER               :: I, J
      REAL*8,SAVE           :: TABLE_S(21), TABLE_U10(20)
      REAL*8,SAVE           :: TABLE_UPTAKE(21,20) 
      

      ! External functions
      REAL*8,  EXTERNAL     :: SFCWINDSQR

      ! Flag for first call
      LOGICAL,SAVE          :: FIRST=.TRUE.

      !=================================================================
      ! HG2_SEASALT_LOSSRATE begins here!
      !=================================================================


!$OMP PARALLEL DO 
!$OMP+DEFAULT( SHARED ) 
!$OMP+PRIVATE( I, J, U10M, S, LOSS_FREQ )
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Only calculate deposition via sea salt over water
#if defined( DEVEL )
         IF ( IS_WATER( I, J, LOCAL_MET ) ) THEN
#else
         IF ( IS_WATER(I,J) ) THEN
#endif

            ! Wind speed at 10m altitude [m/s]
#if defined( DEVEL )
            U10M = SQRT( SFCWINDSQR( I, J,
     &                               LOCAL_MET%U10M, LOCAL_MET%V10M ) )
#else
            U10M = SQRT( SFCWINDSQR( I, J, U, V ) )
#endif
            
            ! Don't allow wind >20 m/s which is the limit of this 
            ! parameterization
            U10M = MAX( MIN( U10M, 20d0 ), 1D0 )

            ! Relative humidity as a saturation ratio
            ! Use the relative humidity of the lowest layer, although this is 
            ! lower than the higher layers of the MBL
            !
            ! Don't allow supersaturation, as [Cl-] is undefined for RH>=1 
            ! Cap RH at 99%, Don't allow RH < 75% as happens in coastal areas
#if defined( DEVEL )
            S = MAX( MIN( LOCAL_MET%RH(I,J,1), 99D0 ), 75D0 ) * 1D-2
#else
            S = MAX( MIN( RH(I,J,1), 99D0 ), 75D0 ) * 1D-2
#endif
            
            LOSS_FREQ = 1D-10 * ( 1D0 - EXP( -57.758D0 * (1D0-S) ) ) *
     &           EXP( -1.9351D0  * U10M + 
     &                 9.0047D0  * SQRT( U10M ) + 
     &                 0.14788D0 * U10M**1.5D0 ) 
            
            ! Loss frequency must be positive
            LOSS_FREQ = MAX( LOSS_FREQ, 1D-10 ) 

         ELSE 

            ! No loss over land
            LOSS_FREQ = 0D0
            
         ENDIF
         
         HG2_SEASALT_LOSSRATE(I,J) = LOSS_FREQ


      ENDDO ! I 
      ENDDO ! J
!$OMP END PARALLEL DO

      ! Return to calling program
      END SUBROUTINE CALC_HG2_SEASALT_LOSSRATE

!------------------------------------------------------------------------------

      SUBROUTINE GET_GLOBAL_JNO2( THISMONTH )
!
!******************************************************************************
!  Subroutine GET_GLOBAL_JNO2 reads monthly mean JNO2 data fields.  
!  These are needed for Hg sim. (cdh, 11/25/09)
!
!  Arguments as Input:
!  ===========================================================================
!  (1 ) THISMONTH (INTEGER) : Current month number (1-12)
!
!  NOTES:
!******************************************************************************

      ! References to F90 modules
      USE BPCH2_MOD,     ONLY : GET_NAME_EXT, GET_RES_EXT
      USE BPCH2_MOD,     ONLY : GET_TAU0,     READ_BPCH2
      USE DIRECTORY_MOD, ONLY : DATA_DIR
      USE TRANSFER_MOD,  ONLY : TRANSFER_3D

      USE CMN_SIZE_MOD   ! Size parameters

      IMPLICIT NONE
#     include "define.h"

      ! Arguments
      INTEGER, INTENT(IN)  :: THISMONTH

      ! Local variables
      REAL*4               :: ARRAY(IIPAR,JJPAR,LGLOB)
      REAL*4               :: ARRAY2(IIPAR,JJPAR,LLPAR)
      REAL*8               :: XTAU
      CHARACTER(LEN=255)   :: FILENAME

      !=================================================================
      ! GET_GLOBAL_JNO2 begins here!
      !=================================================================

      !CDH -WE EVENTUALLY NEED TO SUPPORT OTHER MODEL RESOLUTIONS
      ! Filename for full-level GEOS5 model
      FILENAME = 
     &           TRIM( DATA_DIR )  // 'mercury_201007/jvalues.noon.' //
     &           GET_NAME_EXT() // '.' // GET_RES_EXT()


      ! Echo some information to the standard output
      WRITE( 6, 110 ) TRIM( FILENAME )
 110  FORMAT( '     - GET_GLOBAL_JNO2: Reading ', a )

      ! Get the TAU0 value for the start of the given month
      ! Assume "generic" year 1985 (TAU0 = [0, 744, ... 8016])
      XTAU = GET_TAU0( THISMONTH, 1, 2005 )

#if   defined( GRIDREDUCED )
 
      ! Read O3 data (v/v) from the binary punch file (tracer #51)
      CALL READ_BPCH2( FILENAME, 'JV-MAP-$', 1,     
     &                 XTAU,      IIPAR,     JJPAR,      
     &                 LGLOB,     ARRAY,     QUIET=.TRUE. )

      ! Assign data from ARRAY to the module variable O3
      CALL TRANSFER_3D( ARRAY, JNO2 )

#else

      ! Read O3 data (v/v) from the binary punch file (tracer #51)
      CALL READ_BPCH2( FILENAME, 'JV-MAP-$', 1,      
     &                 XTAU,      IIPAR,     JJPAR,
     &                 LLPAR,     ARRAY2,    QUIET=.TRUE. )

      ! Assign data from ARRAY to the module variable O3
      ! (don't have to fold layers in the stratosphere)
      JNO2 = ARRAY2

#endif


      ! Return to calling program
      END SUBROUTINE GET_GLOBAL_JNO2

!------------------------------------------------------------------------------

      FUNCTION GET_JNO2( I, J, L ) RESULT( JNO2_NOW )
!
!******************************************************************************
!  Function GET_JNO2 returns monthly mean JNO2 and imposes a diurnal 
!
!  Impose the diurnal variation of JNO2 found by Parrish et al. (1983) under 
!  clear skies. J-NO2 ~ exp( -0.360 * sec(SZA) )
!
!  Arguments as Input:
!  ============================================================================
!  (1-3) I, J, L (INTEGER) : Grid box indices for lon, lat, vertical level
!
!  NOTES:
!  01 Mar 2012 - R. Yantosca - Now use GET_YMID_R(I,J,L) from grid_mod.F90
!******************************************************************************
!
      ! References to F90 modules
      USE DAO_MOD,       ONLY : SUNCOS
      USE TIME_MOD,      ONLY : GET_TS_CHEM, GET_DAY_OF_YEAR 
      USE GRID_MOD,      ONLY : GET_YMID_R

      USE CMN_SIZE_MOD  ! Size parameters
      USE CMN_GCTM_MOD    ! Physical constants

      ! Arguments
      INTEGER, INTENT(IN) :: I, J, L

      ! Local variables
      INTEGER             :: JLOOP, JDAY
      REAL*8              :: JNO2_NOW, CSZANOON
      REAL*8              :: A0, A1, A2, A3, B1, B2, B3
      REAL*8              :: R, DEC, YMID_R  
 
      !=================================================================
      ! GET_JNO2 begins here!
      !=================================================================
      
      ! 1-D grid box index for SUNCOS
      JLOOP = ( (J-1) * IIPAR ) + I

      ! Test for sunlight...
      IF ( SUNCOS(JLOOP) > TINY(1D0) ) THEN

         ! Day of year
         JDAY = GET_DAY_OF_YEAR()

         ! Coefficients for solar declination angle
         A0  = 0.006918d0
         A1  = 0.399912d0
         A2  = 0.006758d0
         A3  = 0.002697d0
         B1  = 0.070257d0
         B2  = 0.000907d0
         B3  = 0.000148d0
 
         ! Path length of earth's orbit traversed since Jan 1 [radians]
         R   = ( 2d0 * PI / 365d0 ) * DBLE( JDAY - 1 ) 

         ! Solar declination angle (low precision formula)
         DEC = A0 - A1*COS(     R ) + B1*SIN(     R )
     &            - A2*COS( 2d0*R ) + B2*SIN( 2d0*R )
     &            - A3*COS( 3d0*R ) + B3*SIN( 3d0*R )

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
         ! because we already checked SUNCOS(JLOOP) and we know
         ! CSZANOON >= SUNCOS(JLOOP)
         JNO2_NOW = JNO2(I,J,L) * 
     &        EXP( 0.36D0 * ( 1D0/CSZANOON - 1D0/SUNCOS(JLOOP) ) )

         ! Make sure OH is not negative
         JNO2_NOW= MAX( JNO2_NOW, 0d0 )
               
      ELSE

         ! At night, JNO2 goes to zero
         JNO2_NOW = 0d0

      ENDIF

      ! Return to calling program
      END FUNCTION GET_JNO2

!------------------------------------------------------------------------------

      SUBROUTINE PARTITIONHG( DIRECTION, STT_COPY )

!******************************************************************************
!  SUBROUTINE PARTITIONHG partitions the Hg2 tracer between gas and particle
!  phases for the purpose of wet/dry deposition.
!
!  Arguments as Input/Output:
!  ============================================================================
!  DIRECTION  : 1 for forward , split Hg2 into gas and particle phases
!             : 2 for backward, return gas and particle phases to Hg2 tracer
!  STT_COPY   : STT array of tracer concentrations
!
!  NOTES:
!   12-Jul-2010 H. Amos - Add option to partition Hg2 according to Fg/Fp
!   02-Jan-2011 H. Amos - Add/clean up comments
!   04 Jan 2012 H. Amos - modify algorithms to reflect the fact that 
!                         anthropogenic Hg(p) is now emitted as Hg(II) (i.e. it's
!                         no longer considered refractory). 
!******************************************************************************

      ! References to F90 modules
      USE ERROR_MOD        , ONLY : SAFE_DIV  , GEOS_CHEM_STOP
      USE TRACERID_MOD     , ONLY : ID_HG2    , ID_HGP
      USE TRACER_MOD       , ONLY : N_TRACERS
      USE OCEAN_MERCURY_MOD, ONLY : Fp        , Fg

      USE CMN_SIZE_MOD  ! Size parameters

      INTEGER, INTENT(IN)     :: DIRECTION
      REAL*8, INTENT(INOUT)   :: STT_COPY(IIPAR,JJPAR,LLPAR,N_TRACERS)

      INTEGER                 :: I,J,L


      !=================================================================
      ! PARTITIONHG begins here!
      !=================================================================
      
      SELECT CASE (DIRECTION) 
         CASE (1) ! Split Hg2 into gas and particle phases

            ! Note: Now that anthropogenic Hg(p) is emitted as Hg(II)
            !       and not treated as refractory particulate Hg, the
            !       Hg(p) tracer is empty. (hma, 20120104)


!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
            DO L=1, LLPAR
            DO J=1, JJPAR
            DO I=1, IIPAR

              IF ( LHG2HALFAEROSOL ) THEN 
                !----------------------------
                ! Partition Hg(II) 50/50
                !----------------------------

                ! Put particulate Hg(II) into the Hg(p) tracer
                STT_COPY(I,J,L,ID_HGP(1)) = 
     &               0.5D0 * STT_COPY(I,J,L,ID_HG2(1))
                STT_COPY(I,J,L,ID_HG2(1)) = 0.5D0 *
     &               STT_COPY(I,J,L,ID_HG2(1))
              ELSE 
               !--------------------------------------------
               ! temperature dependent Hg(II) partitioning
               !--------------------------------------------

                ! Put particulate Hg(II) into the Hg(p) tracer
                STT_COPY(I,J,L,ID_HGP(1)) = 
     &               Fp(I,J,L) * STT_COPY(I,J,L,ID_HG2(1))
                STT_COPY(I,J,L,ID_HG2(1)) = Fg(I,J,L) *
     &               STT_COPY(I,J,L,ID_HG2(1))
              ENDIF
                  
            ENDDO
            ENDDO
            ENDDO
!$OMP END PARALLEL DO
      
         CASE(2) ! Recombine particulate Hg2 with  gaseous Hg2
            
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
            DO L=1, LLPAR
            DO J=1, JJPAR
            DO I=1, IIPAR
               STT_COPY(I,J,L,ID_HG2(1)) = STT_COPY(I,J,L,ID_HG2(1)) +
     &                                     STT_COPY(I,J,L,ID_HGP(1))
               STT_COPY(I,J,L,ID_HGP(1)) = 0d0
            ENDDO
            ENDDO
            ENDDO               

!$OMP END PARALLEL DO

         CASE DEFAULT
            WRITE( 6,'(A)') 'Unrecognized conversion direction'
            CALL GEOS_CHEM_STOP

         END SELECT

      ! Return to calling program
      END SUBROUTINE PARTITIONHG

!-----------------------------------------------------------------------------


      SUBROUTINE INIT_MERCURY( THIS_ANTHRO_Hg_YEAR )
!
!******************************************************************************
!  Subroutine INIT_MERCURY allocates and zeroes all module arrays.
!  (eck, cdh, sas, bmy, 12/2/04, 4/6/06)
!  
!  NOTES:
!  (1 ) Removed reference to FEMIS array.  Now also allocates and zeroes
!        the T44 array.  Added reference to CMN_DIAG.  Now references 
!        N_TRACERS from "tracer_mod.f". (bmy, 2/24/05)
!  (2 ) EHg0_an, EHg2_an, EHgP_an are now 2-D arrays.  Now modified for 
!        updated ocean mercury module. (eck, cdh, sas, bmy, 4/6/06)
!  27 Sep 2011 - H. Amos - remove LHg_WETDasHNO3 logical, it's obsolete
!  11 Apr 2012 - R. Yantosca - Now retire lai_mod.F; remove call to INIT_LAI
!  11 Apr 2012 - R. Yantosca - Now move LAI initialization to main.F
!******************************************************************************
!
      ! References to F90 modules
      USE DRYDEP_MOD  , ONLY : DEPNAME  , NUMDEP
      USE ERROR_MOD   , ONLY : ALLOC_ERR, ERROR_STOP
      USE LOGICAL_MOD , ONLY : LSPLIT   , LDRYD     , LNLPBL
      USE LOGICAL_MOD , ONLY : LGTMM
      USE TRACER_MOD  , ONLY : N_TRACERS
      USE TRACERID_MOD, ONLY : N_Hg_CATS
!------------------------------------------------------------------------------
! Prior to 4/11/12:
! Replace lai_mod.F with modis_lai_mod.F.  Also now initialize the LAI
! arrays from main.F instead of from within this module. (bmy, 4/11/12)
!      USE LAI_MOD     , ONLY : INIT_LAI
!------------------------------------------------------------------------------

      USE CMN_SIZE_MOD     ! Size parameters
      USE CMN_DIAG_MOD     ! ND44

      ! Arguments
      INTEGER, INTENT(IN)   :: THIS_ANTHRO_Hg_YEAR

      ! Local variables
      LOGICAL, SAVE         :: IS_INIT = .FALSE. 
      INTEGER               :: AS, N
      CHARACTER(LEN=255)    :: LOCATION

      !=================================================================
      ! INIT_MERCURY begins here!
      !=================================================================

      ! Return if we have already allocated arrays
      IF ( IS_INIT ) RETURN

      ! Anthropogenic Hg emissions year
      ANTHRO_Hg_YEAR = THIS_ANTHRO_Hg_YEAR

      ! Location string for ERROR_STOP
      LOCATION       = 'DEFINE_TAGGED_Hg ("mercury_mod.f")'

      !=================================================================
      ! Allocate arrays
      !=================================================================
      ALLOCATE( COSZM( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'COSZM' )
      COSZM = 0d0

      ALLOCATE( EHg0_an( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EHg0_an' )
      EHg0_an = 0d0

      ALLOCATE( EHg0_am( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EHg0_am' )
      EHg0_am = 0d0

      ALLOCATE( EHg2_an( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EHg2_an' )
      EHg2_an = 0d0

      ALLOCATE( EHgP_an( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EHgP_an' )
      EHgP_an = 0d0

      ALLOCATE( EHg0_oc( IIPAR, JJPAR, N_Hg_CATS ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EHg0_oc' )
      EHg0_oc = 0d0

      ALLOCATE( EHg0_dist( IIPAR, JJPAR), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EHg0_dist' )
      EHg0_dist = 0d0

      ALLOCATE( EHg0_nt( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EHg0_nt' )
      EHg0_nt = 0d0

      ALLOCATE( EHg0_bb( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EHg0_bb' )
      EHg0_bb = 0d0

      ALLOCATE( EHg0_snow( IIPAR, JJPAR, N_Hg_CATS ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EHg0_snow' )
      EHg0_snow = 0d0

      IF ( LGTMM ) THEN

         ALLOCATE( EHg0_gtm( IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'EHg0_gtm' )
         EHg0_gtm = 0d0

      ELSE

         ALLOCATE( EHg0_ln( IIPAR, JJPAR, N_Hg_CATS ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'EHg0_ln' )
         EHg0_ln = 0d0

         ALLOCATE( EHg0_vg( IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'EHg0_vg' )
         EHg0_vg = 0d0

         ALLOCATE( EHg0_so( IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'EHg0_so' )
         EHg0_so = 0d0

      ENDIF

      ALLOCATE( TCOSZ( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'TCOSZ' )
      TCOSZ = 0d0

      ALLOCATE( TTDAY( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'TTDAY' )
      TTDAY = 0d0
 
      ! Allocate ZERO_DVEL if we use non-local PBL mixing or
      ! if drydep is turned off 
      IF ( LNLPBL .OR. (.not. LDRYD) ) THEN
         ALLOCATE( ZERO_DVEL( IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'ZERO_DVEL' )
         ZERO_DVEL = 0d0
      ENDIF

      ALLOCATE( HG2_SEASALT_LOSSRATE( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'HG2_SEASALT_LOSSRATE' )
      HG2_SEASALT_LOSSRATE = 0d0

      ALLOCATE( JNO2( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'JNO2' )
      JNO2 = 0d0

      ALLOCATE( HGPFRAC( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'HGPFRAC' )
      HGPFRAC = 0d0

!------------------------------------------------------------------------------
! Prior to 4/11/12:
! Replace lai_mod.F with modis_lai_mod.F.  Also now initialize the LAI
! arrays from main.F instead of from within this module. (bmy, 4/11/12)
!      CALL INIT_LAI
!------------------------------------------------------------------------------

      !=================================================================
      ! Allocate & initialize arrays for tagged tracer indices
      !=================================================================
      IF ( LSPLIT ) THEN

         ! Tracer indices for tagged anthro regions
         ALLOCATE( AN_Hg0( IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'AN_Hg0' )
         AN_Hg0 = 0d0

         ALLOCATE( AN_Hg2( IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'AN_Hg2' )
         AN_Hg2 = 0d0

         ALLOCATE( AN_HgP( IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'AN_HgP' )
         AN_HgP = 0d0

         ! Define the tagged tracer indices 
         CALL DEFINE_TAGGED_Hg

      ENDIF

      !=================================================================
      ! Settings
      !=================================================================

      ! Switch determines whether uptake of Hg2 by sea-salt aerosol
      ! is calculated dynamically (TRUE) or uses a constant rate (FALSE)
      LDYNSEASALT = .TRUE.

      ! Switch uses Bess' best-guess emissions if TRUE or emissions from
      ! Selin et al. 2008 (+bug fixes) if FALSE
      ! This is temporary
      ! Bess's changes: rescale GEIA to use Streets' 2006 regional totals
      ! With these emissions, Bess found the best chemistry results
      LGCAPEMIS=.FALSE.
      ! Use GEIA 2005 inventory
      LGEIA05=.TRUE.

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
      LGEOSLWC=.TRUE. ! Use GEOS LWC

      ! Switch specifies that Hg2 is 50% bound to aerosol and 50% in
      ! gas phase (TRUE). If FALSE, then use temperature dependent 
      ! partitioning as described in Amos et al. (2011, ACPD)
      LHg2HalfAerosol = .FALSE.

      ! Switch turns on snowpack Hg storage until snowmelt
      LHGSNOW = .TRUE.

      ! Multiplicative factor for increasing stratospheric Br and BrO
      STRAT_BR_FACTOR = 1d0

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

!------------------------------------------------------------------------------

      SUBROUTINE CLEANUP_MERCURY
!
!******************************************************************************
!  Subroutine CLEANUP_MERCURY deallocates all module arrays.
!  (eck, bmy, 12/6/04, 2/24/05)
!
!  NOTES:
!  (1 ) Now deallocate MLD, NPP, RAD (sas, bmy, 1/18/05)
!  (2 ) Now deallocate T44 (bmy, 2/24/05)
!******************************************************************************
!
      !=================================================================
      ! CLEANUP_MERCURY begins here!
      !=================================================================
      IF ( ALLOCATED( AN_Hg0    ) ) DEALLOCATE( AN_Hg0    )  
      IF ( ALLOCATED( AN_Hg2    ) ) DEALLOCATE( AN_Hg2    )      
      IF ( ALLOCATED( AN_HgP    ) ) DEALLOCATE( AN_HgP    )      
      IF ( ALLOCATED( COSZM     ) ) DEALLOCATE( COSZM     )       
      IF ( ALLOCATED( EHg0_an   ) ) DEALLOCATE( EHg0_an   )
      IF ( ALLOCATED( EHg0_am   ) ) DEALLOCATE( EHg0_am   )
      IF ( ALLOCATED( EHg2_an   ) ) DEALLOCATE( EHg2_an   )
      IF ( ALLOCATED( EHgP_an   ) ) DEALLOCATE( EHgP_an   )
      IF ( ALLOCATED( EHg0_oc   ) ) DEALLOCATE( EHg0_oc   )
      IF ( ALLOCATED( EHg0_ln   ) ) DEALLOCATE( EHg0_ln   )
      IF ( ALLOCATED( EHg0_snow ) ) DEALLOCATE( EHg0_snow )
      IF ( ALLOCATED( EHg0_nt   ) ) DEALLOCATE( EHg0_nt   )
      IF ( ALLOCATED( EHg0_bb   ) ) DEALLOCATE( EHg0_bb   )
      IF ( ALLOCATED( EHg0_gtm  ) ) DEALLOCATE( EHg0_gtm  )
      IF ( ALLOCATED( EHg0_vg   ) ) DEALLOCATE( EHg0_vg   )
      IF ( ALLOCATED( EHg0_so   ) ) DEALLOCATE( EHg0_so   )
      IF ( ALLOCATED( EHg0_dist ) ) DEALLOCATE( EHg0_dist )
      IF ( ALLOCATED( TCOSZ     ) ) DEALLOCATE( TCOSZ     )
      IF ( ALLOCATED( TTDAY     ) ) DEALLOCATE( TTDAY     )
      IF ( ALLOCATED( ZERO_DVEL ) ) DEALLOCATE( ZERO_DVEL )
      IF ( ALLOCATED( HG2_SEASALT_LOSSRATE ) ) 
     &     DEALLOCATE( HG2_SEASALT_LOSSRATE   ) 
      IF ( ALLOCATED( JNO2      ) ) DEALLOCATE( JNO2      ) 
  
      ! Return to calling program
      END SUBROUTINE CLEANUP_MERCURY

!------------------------------------------------------------------------------

      ! End of module
      END MODULE MERCURY_MOD
