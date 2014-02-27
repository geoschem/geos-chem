!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: sulfate_mod
!
! !DESCRIPTION: Module SULFATE\_MOD contains arrays and routines for performing
!  either a coupled chemistry/aerosol run or an offline sulfate aerosol
!  simulation. Original code taken from Mian Chin's GOCART model and modified
!  accordingly. (rjp, bdf, bmy, 6/22/00, 8/26/10)
!\\
!\\
! !INTERFACE: 
!
      MODULE SULFATE_MOD
!
! !USES:
!
      USE VDIFF_PRE_MOD, ONLY : emis_save ! (Lin, 03/31/09)
      USE inquireMod,    ONLY : findFreeLUN

      IMPLICIT NONE
      PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS: 
!
      PUBLIC :: CHEMSULFATE       
      PUBLIC :: EMISSSULFATE      
      PUBLIC :: CLEANUP_SULFATE   
      PUBLIC :: INIT_SULFATE
!
! 
! !REMARKS:
!  References:
!  ============================================================================
!  (1 ) Andreae, M.O. & P. Merlet, "Emission of trace gases and aerosols from
!        biomass burning", Global Biogeochem. Cycles, 15, 955-966, 2001.
!  (2 ) Nightingale et al [2000a], J. Geophys. Res, 14, 373-387
!  (3 ) Nightingale et al [2000b], Geophys. Res. Lett, 27, 2117-2120
!  (4 ) Wanninkhof, R., "Relation between wind speed and gas exchange over
!        the ocean", J. Geophys. Res, 97, 7373-7382, 1992.
! 
! !REVISION HISTORY: 
!  (1 ) All module variables are declared PRIVATE (i.e., they can only
!        be seen from within this module (bmy, 6/2/00)
!  (2 ) The routines in "sulfate_mod.f" assume that we are doing chemistry
!        over the global region (e.g. IIPAR=IIPAR, JJPAR=JJPAR). (bmy, 6/8/00)
!  (3 ) Removed obsolete code from DRYDEP_SULFATE (bmy, 12/21/00)
!  (4 ) Removed obsolete commented-out code from module routines (bmy, 4/23/01)
!  (5 ) Now read data files from DATA_DIR/sulfate_sim_200106/ (bmy, 6/19/01)
!  (6 ) Updated comments (bmy, 9/4/01)
!  (7 ) XTRA2(IREF,JREF,5) is now XTRA2(I,J).  Now reference COSSZA from
!        "dao_mod.f". (bmy, 9/27/01)
!  (8 ) Removed obsolete commented out code from 9/01 (bmy, 10/24/01)
!  (9 ) Minor fixes to facilitate compilation on ALPHA (bmy, 11/15/01)
!  (11) Now divide module header into MODULE PRIVATE, MODULE VARIABLES, and
!        MODULE ROUTINES sections.  Updated comments (bmy, 5/28/02)
!  (12) Replaced all instances of IM with IIPAR and JM with JJPAR, in order
!        to prevent namespace confusion for the new TPCORE (bmy, 6/25/02)
!  (13) Now reference "file_mod.f" (bmy, 6/27/02)
!  (14) Now references GET_PEDGE from "pressure_mod.f", which computes P at
!        the bottom edge of grid box (I,J,L).  Also deleted obsolete,
!        commented-out code. (dsa, bdf, bmy, 8/21/02)
!  (15) Added updated code from Rokjin Park and Brendan Field, in order to
!        perform coupled chemistry-aerosol simulations.  Also added parallel
!        DO-loops in several subroutines.  Updated comments, cosmetic
!        changes.  Now reference "error_mod.f" and "wetscav_mod.f".  
!        Now only do chemistry below the tropopause. (rjp, bdf, bmy, 12/6/02)
!  (16) Added ENH3_na array to hold natural source NH3 emissions.  Also now
!        facilitate passing DMS, SO2, SO4, NH3 to SMVGEAR for fullchem
!        simulations.  Added subroutine READ_NATURAL_NH3. (rjp, bmy, 3/23/03)
!  (17) Now references "grid_mod.f" and "time_mod.f".  Also made other minor
!        cosmetic changes. (bmy, 3/27/03)
!  (18) Updated chemistry routines to apply drydep losses throughout the
!        entire PBL. (rjp, bmy, 8/1/03)
!  (19) Now accounts for GEOS-4 PBL being in meters (bmy, 1/15/04)
!  (20) Fix ND44 diag so that we get same results for sp or mp (bmy, 3/24/04)
!  (21) Added COSZM array.  Now use diurnal varying JH2O2 in CHEM_H2O2. 
!        (rjp, bmy, 3/39/04)
!  (22) Added more parallel DO-loops (bmy, 4/14/04)
!  (23) Now add SO2 from ships (bec, bmy, 5/20/04)
!  (24) Now references "directory_mod.f", "logical_mod.f" and "tracer_mod.f".
!        Now removed IJSURF. (bmy, 7/20/04)
!  (25) Can overwrite USA with EPA/NEI99 emissions (rjp, rch, bmy, 11/16/04)
!  (26) Modified for AS, AHS, LET, SO4aq, NH4aq (cas, bmy, 1/11/05)
!  (27) Now also references "pbl_mix_mod.f".  NOTE: Comment out phase 
!        transition  code for now since it is still under development and
!        will take a while to be rewritten. (bmy, 3/15/05)
!  (28) Modified for SO4s, NITs chemistry (bec, 4/13/05)
!  (29) Now reads updated files for SST and offline chemistry.  Now read data
!        for both GCAP and GEOS grids.  Now references "tropopause_mod.f".
!        (bmy, 8/22/05)
!  (30) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (31) Now references XNUMOL & XNUMOLAIR from "tracer_mod.f" (bmy, 10/25/05)
!  (32) Now read int'annual SST data on GEOS 1x1 grid (bmy, 11/17/05)
!  (33) Bug fix for offline aerosol sim in SEASALT_CHEM (bec, bmy, 3/29/06)
!  (34) Bug fix in INIT_DRYDEP (bmy, 5/23/06)
!  (35) Now references "bravo_mod.f" (rjp, kfb, bmy, 6/26/06)
!  (36) Now references "streets_anthro_mod.f" (yxw, bmy, 8/17/06)
!  (37) Now references "biomass_mod.f" (bmy, 9/27/06)
!  (38) Now prevent seg fault error in READ_BIOFUEL_SO2 (bmy, 11/3/06)
!  (39) Bug fix in SEASALT_CHEM (havala, bec, bmy, 12/8/06)
!  (40) Extra error check for low RH in GRAV_SETTLING (phs, 6/11/08)
!  (41) Now references "cac_anthro_mod.f".  And apply SO2 yearly scale factor
!        to SO2 from GEIA (amv, phs, 3/11/08)  
!  (41) Bug fixes in reading EDGAR data w/ the right tracer number, 
!        when we are doing offline or nonstd simulations (dkh, 10/31/08)
!  (42) Bug fix for AD13_SO2_sh in SRCSO2 (phs, 2/27/09)
!  (43) Bug fix: need to add CAC_AN to PRIVATE statements (bmy, 5/27/09)
!  (44) Constrain surface emissions to the first level and save them into
!        emis_save (lin, 5/29/09)
!  (45) Last year of SST data is now 2008 (see READ_SST) (bmy, 7/13/09)
!  (46) Updated rxns in CHEM_DMS and CHEM_SO2 to JPL 2006 (jaf, bmy, 10/15/09)
!  (47) Added new volcanic emissions of SO2 (jaf, bmy, 10/15/09)
!  (48) Now accounts for NEI 2005 emissions, and multilevels SOxan emissions
!        (amv, phs, 10/15/2009) 
!  (49) Fixes in SRCSO2 for SunStudio compiler (bmy, 12/3/09)
!  (50) Add new subroutine SRCSF30 for emission to 30bin sulfate (win, 1/25/10)
!  (51) Add new array PSO4_SO2AQ for SO4 produced via aqueous chemistry of SO2 
!        excluding that from heterogeous reaction on sea-salt. (win, 1/25/10)
!  (52) Standardized patch in READ_ANTHRO_NH3 (dkh, bmy, 3/5/10)
!  (53) Use LWC from GEOS-5 met fields (jaf, bmy, 6/30/10)
!  (54) Add module parameters MNYEAR_VOLC and MXYEAR_VOLC to define the 1st 
!       and last year with data for volcanic emissions. (ccc, 9/30/10)
!  (55) Use updated volcanic emissions from 1979 to 2009
!  26 Aug 2010 - R. Yantosca   - Add modifications for MERRA
!  12 Nov 2010 - R. Yantosca   - Avoid div-by-zero when computing L2S, L3S
!  07 Sep 2011 - P. Kasibathla - Modified to include GFED3
!  22 Dec 2011 - M. Payer      - Added ProTeX headers 
!  08 Feb 2012 - R. Yantosca   - Add modifications for GEOS-5.7.2 met
!  01 Mar 2012 - R. Yantosca   - Now reference new grid_mod.F90
!  13 Mar 2012 - M. Cooper     - Changed regrid algorithm to map_a2a
!  28 Nov 2012 - R. Yantosca   - Use SUNCOS fields from the State_Met object
!  04 Mar 2013 - R. Yantosca   - Now call INIT_SULFATE from the init stage
!                                which facilitates connection to GEOS-5 GCM
!  05 Mar 2013 - R. Yantosca   - Now use Input_Opt%LNLPBL instead of LNLPBL
!                                from logical_mod.F
!  13 Mar 2013 - R. Yantosca   - Bug fix: make sure we pass values to the
!                                SOIL_DRYDEP routine even when ND44 is off
!  30 May 2013 - S. Farina     - Merged TOMAS code into sulfate_mod.F
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
      !========================================================================
      ! MODULE PARAMETERS:
      !
      ! XNUMOL_OH  : Molecules OH  per kg OH          [molec/kg]
      ! XNUMOL_O3  : Molecules O3  per kg O3          [molec/kg]
      ! XNUMOL_NO3 : Molecules NO3 per kg NO3         [molec/kg]
      ! TCVV_S     : Ratio: Molwt air / Molwt S       [unitless]
      !=======================================================================
      REAL*8,  PARAMETER   :: XNUMOL_OH  = 6.022d23 / 17d-3
      REAL*8,  PARAMETER   :: XNUMOL_O3  = 6.022d23 / 48d-3
      REAL*8,  PARAMETER   :: XNUMOL_NO3 = 6.022d23 / 62d-3
      REAL*8,  PARAMETER   :: TCVV_S     = 28.97d0  / 32d0 
      REAL*8,  PARAMETER   :: SMALLNUM   = 1d-20

#if defined( TOMAS )
      !---------------------------------------------------------------
      ! For TOMAS microphysics: Add parameter for scaling anthro SO2
      !---------------------------------------------------------------
      REAL*8,  PARAMETER   :: scaleanthso2 = 1.0d0
#endif
!
! !PRIVATE TYPES:
!
      !========================================================================
      ! MODULE VARIABLES:
      !
      ! DMSo       : DMS oceanic emissions            [v/v/timestep]
      ! DRYH2O2    : Pointer to H2O2  in DEPVEL array [unitless] 
      ! DRYSO2     : Pointer to SO2   in DEPVEL array [unitless]
      ! DRYSO4     : Pointer to SO4   in DEPVEL array [unitless]
      ! DRYSO4s    : Pointer to SO4s  in DEPVEL array [unitless]
      ! DRYMSA     : Pointer to MSA   in DEPVEL array [unitless]
      ! DRYNH3     : Pointer to NH3   in DEPVEL array [unitless]
      ! DRYNH4     : Pointer to NH4   in DEPVEL array [unitless]
      ! DRYNIT     : Pointer to NIT   in DEPVEL array [unitless]
      ! DRYNITs    : Pointer to NITs  in DEPVEL array [unitless]
      ! DRYSO4aq   : Pointer to SO4aq in DEPVEL array [unitless]
      ! DRYAS      : Pointer to AS    in DEPVEL array [unitless]  
      ! DRYAHS     : Pointer to AHS   in DEPVEL array [unitless]
      ! DRYLET     : Pointer to LET   in DEPVEL array [unitless]
      ! DRYNH4aq   : Pointer to NH4aq in DEPVEL array [unitless]
      ! ENH3_an    : NH3 anthropogenic emissions      [kg NH3/box/s]
      ! ENH3_bb    : NH3 biomass emissions            [kg NH3/box/s]
      ! ENH3_bf    : NH3 biofuel emissions            [kg NH3/box/s]
      ! ENH3_na    : NH73 natural source emissions    [kg NH3/box/s]
      ! ESO2_ac    : SO2 aircraft emissions           [kg SO2/box/s]
      ! ESO2_an    : SO2 anthropogenic emissions      [kg SO2/box/s]
      ! ESO2_ev    : SO2 eruptive volcanic em.        [kg SO2/box/s]
      ! ESO2_nv    : SO2 non-eruptive volcanic em.    [kg SO2/box/s]
      ! ESO2_bb    : SO2 biomass burning emissions    [kg SO2/box/s]
      ! ESO2_bf    : SO2 biofuel burning emissions    [kg SO2/box/s]
      ! ESO2_sh    : SO2 ship emissions               [kg SO2/box/s]
      ! ESO4_an    : SO4 anthropogenic emissions      [kg SO4/box/s]
      ! JH2O2      : Monthly mean J(H2O2) values      [s-1]
      ! O3m        : Monthly mean O3 concentration    [v/v]
      ! PH2O2m     : Monthly mean P(H2O2)             [molec/cm3/s]
      ! PMSA_DMS   : P(MSA) from DMS                  [v/v/timestep]
      ! PSO2_DMS   : P(SO2) from DMS                  [v/v/timestep]
      ! PSO4_SO2   : P(SO4) from SO2                  [v/v/timestep]
      ! SSTEMP     : Sea surface temperatures         [K]
      ! VCLDF      : Volume cloud frac. for SO2 aq.   [unitless]
      ! Eev        : SO2 em. from eruptive volcanoes  [kg SO2/box/s]
      ! Env        : SO2 em. from non-erup volcanoes  [kg SO2/box/s]
      ! TCOSZ      : Sum of cos(SZA) for offline run  [unitless] 
      ! TTDAY      : Total daylight length at (I,J)   [minutes]
      ! SMALLNUM   : Small number - prevent underflow [unitless]
      ! COSZM      : Array for MAX(cos(SZA)) at (I,J) [unitless]
      ! LVOLC      : Number of volcanic levels (20)   [unitless]
      !========================================================================

      ! Time variable
      INTEGER              :: ELAPSED_SEC

      ! Allocatable arrays
      REAL*8,  ALLOCATABLE :: DMSo(:,:) 
      REAL*8,  ALLOCATABLE :: ENH3_an(:,:)
      REAL*8,  ALLOCATABLE :: ENH3_bb(:,:)
      REAL*8,  ALLOCATABLE :: ENH3_bf(:,:)
      REAL*8,  ALLOCATABLE :: ENH3_na(:,:)
      REAL*8,  ALLOCATABLE :: ESO2_ac(:,:,:) 
      REAL*8,  ALLOCATABLE :: ESO2_an(:,:,:)
      REAL*8,  ALLOCATABLE :: ESO2_bb(:,:)     
      REAL*8,  ALLOCATABLE :: ESO2_bf(:,:)
      REAL*8,  ALLOCATABLE :: ESO2_ev(:,:,:)
      REAL*8,  ALLOCATABLE :: ESO2_nv(:,:,:)
      REAL*8,  ALLOCATABLE :: ESO2_sh(:,:) 
      REAL*8,  ALLOCATABLE :: ESO4_an(:,:,:) 
      REAL*8,  ALLOCATABLE :: JH2O2(:,:,:)
! not used for now - comment it (phs, 6/29/09)      
!      REAL*8,  ALLOCATABLE :: LSO2_AQ(:,:,:)
      REAL*8,  ALLOCATABLE :: O3m(:,:,:)
      REAL*8,  ALLOCATABLE :: PH2O2m(:,:,:)
      REAL*8,  ALLOCATABLE :: PMSA_DMS(:,:,:)
      REAL*8,  ALLOCATABLE :: PSO2_DMS(:,:,:)
      REAL*8,  ALLOCATABLE :: PSO4_SO2(:,:,:)
      REAL*8,  ALLOCATABLE :: PSO4_SS(:,:,:)
      REAL*8,  ALLOCATABLE :: PNITs(:,:,:)
      REAL*4,  ALLOCATABLE :: SOx_SCALE(:,:)
      REAL*8,  ALLOCATABLE :: SSTEMP(:,:)
      REAL*8,  ALLOCATABLE :: TCOSZ(:,:)
      REAL*8,  ALLOCATABLE :: TTDAY(:,:)
      REAL*8,  ALLOCATABLE :: VCLDF(:,:,:)
      REAL*8,  ALLOCATABLE :: COSZM(:,:)

#if   defined( TOMAS )
      !---------------------------------------------------------------
      ! For TOMAS microphysics: Define PSO4_SO2aq array
      !---------------------------------------------------------------
      REAL*8,  ALLOCATABLE :: PSO4_SO2AQ(:,:,:)
#endif

      ! Volcanic emissions parameters (ccc, 9/30/10)
      INTEGER, PARAMETER   :: LVOLC=20            ! Volcano levels
      INTEGER, PARAMETER   :: MNYEAR_VOLC = 1979  ! Min year of data
      INTEGER, PARAMETER   :: MXYEAR_VOLC = 2009  ! Max year of data

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%% NOTE: This is for the SEAC4RS volcano emissions.  Leave code commented
!%%% out here so that we can restore it later. (bmy, 9/25/13)
!
!      ! Additional variables for volcano emissions (jaf, skim, 1/11/13)
!      !    Includes extra info for eruptive & noneruptive volcanos
!      INTEGER              :: VOLCYEAR            ! Year of data to use
!      INTEGER              :: N_ERUP, N_NONERUP   ! # of each volc events
!      LOGICAL              :: DO_ERUPTIVE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      ! Eruptive volcanoes
      REAL*8,  ALLOCATABLE :: EEV(:,:,:)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%% NOTE: This is for the SEAC4RS volcano emissions.  Leave code commented
!%%% out here so that we can restore it later. (bmy, 9/25/13)
!      REAL*8,  ALLOCATABLE :: VEV_ELEV(:), VEV_CLOUDTOP(:), VEV_SO2(:)
!      INTEGER, ALLOCATABLE :: VEV_DATE(:), VEV_II(:), VEV_JJ(:)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      ! Non-eruptive volcanoes 
      REAL*8,  ALLOCATABLE :: ENV(:,:,:)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%% NOTE: This is for the SEAC4RS volcano emissions.  Leave code commented
!%%% out here so that we can restore it later. (bmy, 9/25/13)
!      REAL*8,  ALLOCATABLE :: VNV_ELEV(:), VNV_SO2(:)
!      INTEGER, ALLOCATABLE :: VNV_DATE(:), VNV_II(:), VNV_JJ(:)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      ! Pointers to drydep species w/in DEPSAV
      INTEGER              :: DRYSO2,  DRYSO4,   DRYMSA,  DRYNH3  
      INTEGER              :: DRYNH4,  DRYNIT,   DRYSO4s, DRYNITs
      INTEGER              :: DRYH2O2, DRYSO4aq, DRYAS,   DRYAHS
      INTEGER              :: DRYLET,  DRYNH4aq

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
! !IROUTINE: get_vcldf
!
! !DESCRIPTION: Subroutine GET\_VCLDF computes the volume cloud fraction for
!  SO2 chemistry. (rjp, bdf, bmy, 9/23/02)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE GET_VCLDF( am_I_Root, State_Met, RC )
!
! !USES:
!
      USE CMN_SIZE_MOD
      USE GIGC_ErrCode_Mod
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE PRESSURE_MOD,       ONLY : GET_PCENTER
      USE PRESSURE_MOD,       ONLY : GET_PEDGE
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)  :: am_I_Root   ! Are we on the root CPU?
      TYPE(MetState), INTENT(IN)  :: State_Met   ! Meteorology State object
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT) :: RC          ! Success or failure?
!
! !REMARKS:
!  References:
!  ============================================================================
!  (1) Sundqvist et al. [1989]
!
! !REVISION HISTORY: 
!  14 Jan 2011 - R. Yantosca - Return if VCLDF is not allocated
!  22 Dec 2011 - M. Payer    - Added ProTeX headers
!  14 Nov 2012 - R. Yantosca - Added am_I_Root, RC arguments
!  15 Nov 2012 - M. Payer    - Replaced all met field arrays with State_Met
!                              derived type object
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
      REAL*8,  PARAMETER   :: ZRT = 0.60d0, ZRS = 0.99d0
!
! !LOCAL VARIABLES:
!
      INTEGER              :: I,    J,    L
      REAL*8               :: PRES, PSFC, RH2, R0, B0
	
      !=================================================================
      ! GET_VCLDF begins here!
      !=================================================================

      ! Assume success
      RC  = GIGC_SUCCESS

      ! Exit if VCLDF is not allocated.  We will now get the cloud
      ! fraction from the GEOS-5 or MERRA met fields. (skim, bmy, 1/14/10)
      IF ( .not. ALLOCATED( VCLDF ) ) RETURN

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, PSFC, PRES, RH2, R0, B0 )
      DO L = 1, LLTROP
      DO J = 1, JJPAR 
      DO I = 1, IIPAR
	
         ! Surface pressure
         PSFC = GET_PEDGE(I,J,1)

         ! Pressure at the center of the grid box
         PRES = GET_PCENTER(I,J,L)

         ! RH (from "dao_mod.f") is relative humidity [%]
         ! Convert to fraction and store in RH2
         RH2  = State_Met%RH(I,J,L) * 1.0d-2

         ! Terms from Sundqvist ???
         R0   = ZRT + ( ZRS - ZRT ) * EXP( 1d0 - ( PSFC / PRES )**2.5 )
         B0   = ( RH2 - R0 ) / ( 1d0 - R0 )
	   
         ! Force B0 into the range 0-1
         IF ( RH2 < R0  ) B0 = 0d0
         IF ( B0  > 1d0 ) B0 = 1d0

         ! Volume cloud fraction
         VCLDF(I,J,L) = 1d0 - SQRT( 1d0 - B0 )

      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      END SUBROUTINE GET_VCLDF
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_lwc
!
! !DESCRIPTION: Function GET\_LWC returns the cloud liquid water content
!  [m3 H2O/m3 air] at a  GEOS-CHEM grid box as a function of temperature.
!  (rjp, bmy, 10/31/02, 1/14/03)
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_LWC( T ) RESULT( LWC )
!
! !INPUT PARAMETERS: 
!
      REAL*8, INTENT(IN) :: T ! Temperature value at a GEOS-CHEM grid box [K]
!
! !RETURN VALUE:
!
      REAL*8             :: LWC
! 
! !REVISION HISTORY: 
!  18 Jan 2011 - R. Yantosca - Updated comments 
!  22 Dec 2011 - M. Payer    - Added ProTeX header
!EOP
!------------------------------------------------------------------------------
!BOC
!

      !=================================================================
      ! GET_LWC begins here!
      !=================================================================

      ! Compute Liquid water content in [g/m3]
      IF ( T > 293d0 ) THEN
         LWC = 0.2d0

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

      END FUNCTION GET_LWC
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: chemsulfate
!
! !DESCRIPTION: Subroutine CHEMSULFATE is the interface between the GEOS-CHEM
!  main program and the sulfate chemistry routines.  The user has the option of
!  running a coupled chemistry-aerosols simulation or an offline aerosol
!  simulation. (rjp, bdf, bmy, 5/31/00, 3/16/06)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CHEMSULFATE( am_I_Root, Input_Opt, 
     &                        State_Met, State_Chm, RC )
!
! !USES:
!
      USE CMN_SIZE_MOD                   
      USE DAO_MOD,            ONLY : CONVERT_UNITS
      USE DRYDEP_MOD,         ONLY : DEPSAV
      USE ERROR_MOD,          ONLY : DEBUG_MSG
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE GIGC_State_Chm_Mod, ONLY : ChmState
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE GLOBAL_OH_MOD,      ONLY : GET_GLOBAL_OH
      USE GLOBAL_NO3_MOD,     ONLY : GET_GLOBAL_NO3
      USE TIME_MOD,           ONLY : GET_MONTH
      USE TIME_MOD,           ONLY : GET_TS_CHEM
      USE TIME_MOD,           ONLY : GET_ELAPSED_SEC
      USE TIME_MOD,           ONLY : ITS_A_NEW_MONTH
      USE TRACERID_MOD,       ONLY : IDTNITs
      USE TRACERID_MOD,       ONLY : IDTSO4s
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)    :: am_I_Root   ! Is this the root CPU?
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
!  (1 ) Now reference all arguments except FIRSTCHEM and RH from either F90 
!        modules or from common block header files.  Updated comments, 
!        cosmetic changes.  Added NH3, NH4, NITRATE chemistry routines.   
!        Also call MAKE_RH and CONVERT_UNITS from "dao_mod.f".  Now references
!        IDTDMS, IDTSO2 etc. from "tracerid_mod.f".  Now make FIRSTCHEM a 
!        local SAVEd variable.  Now reference DEPSAV from "drydep_mod.f".
!        Also get rid of extraneous dimensions of DEPSAV.  Added NTIME,
!        NHMSb arrays for OHNO3TIME.  (rjp, bdf, bmy, 12/16/02)
!  (2 ) CHEM_DMS is now only called for offline sulfate simulations.  
!        (rjp, bmy, 3/23/03)
!  (3 ) Now remove NTIME, NHMSb from the arg list and call to OHNO3TIME.
!        Now references functions GET_MONTH, GET_TS_CHEM, and GET_ELAPSED_SEC
!        from the new "time_mod.f". (bmy, 3/27/03)
!  (4 ) Now reference STT, TCVV, N_TRACERS, ITS_AN_AEROSOL_SIM from
!        "tracer_mod.f".  Now reference ITS_A_NEW_MONTH from "time_mod.f".
!        Now references LPRT from "logical_mod.f". (bmy, 7/20/04)
!  (5 ) Updated for AS, AHS, LET, SO4aq, NH4aq.  Now references LCRYST from
!        logical_mod.f.  Now locate species in the DEPSAV array w/in 
!        INIT_SULFATE. (bmy, 12/21/04)
!  (6 ) Now handle gravitational settling of SO4s, NITs (bec, bmy, 4/13/05)
!  (7 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (8 ) Remove reference to MAKE_RH, it's not needed here (bmy, 3/16/06)
!  (9 ) Reference to LTOMAS and add call CHEM_SO4_AQ using aqueous oxidation
!        which is one of the TOMAS microphysics subroutine  (win, 1/25/10)
!  05 Oct 2011 - R. Yantosca - SUNCOS is no longer needed here
!  22 Dec 2011 - M. Payer    - Added ProTeX headers
!  30 Jul 2012 - R. Yantosca - Now accept am_I_Root as an argument when
!                              running with the traditional driver main.F
!  14 Nov 2012 - R. Yantosca - Add Input_Opt, RC as arguments
!  15 Nov 2012 - M. Payer    - Replaced all met field arrays with State_Met
!                              derived type object
!  04 Mar 2013 - R. Yantosca - Remove call to INIT_SULFATE
!  19 Mar 2013 - R. Yantosca - Now copy Input_Opt%TCVV(1:N_TRACERS)
!  25 Mar 2013 - M. Payer    - Now pass State_Chm object via the arg list
!  23 Apr 2013 - R. Yantosca - Remove LTOMAS logical, since we now invoke TOMAS
!                              with either TOMAS=yes or TOMAS40=yes
!  31 May 2013 - R. Yantosca - Now pass am_I_root, Input_Opt, State_Chm
!                              and RC to TOMAS routine CHEM_SO4_AQ
!  23 Oct 2013 - R. Yantosca - Now pass objects to GET_GLOBAL_OH routine
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! SAVEd scalars
      LOGICAL, SAVE     :: FIRSTCHEM = .TRUE.
      INTEGER, SAVE     :: LASTMONTH = -99

      ! Non-SAVEd scalars
      LOGICAL           :: IT_IS_AN_AEROSOL_SIM
      LOGICAL           :: prtDebug
      INTEGER           :: I, J, L, N, MONTH
      REAL*8            :: DTCHEM

      ! For fields from Input_Opt
      LOGICAL           :: LCRYST
      LOGICAL           :: LPRT
      INTEGER           :: N_TRACERS
      REAL*8            :: TCVV(Input_Opt%N_TRACERS)

      ! Pointers
      ! We need to define local arrays to hold corresponding values 
      ! from the Chemistry State (State_Chm) object. (mpayer, 12/6/12)
      REAL*8, POINTER   :: STT(:,:,:,:)
!
! !EXTERNAL FUNCTIONS:   
!
      REAL*8,  EXTERNAL :: BOXVL

      !=================================================================
      ! CHEMSULFATE begins here!
      !=================================================================

      ! Assume success
      RC                   = GIGC_SUCCESS

      ! Copy fields from INPUT_OPT to local variables for use below
      LCRYST               = Input_Opt%LCRYST
      LPRT                 = Input_Opt%LPRT
      N_TRACERS            = Input_Opt%N_TRACERS
      TCVV                 = Input_Opt%TCVV(1:N_TRACERS)
      IT_IS_AN_AEROSOL_SIM = Input_Opt%ITS_AN_AEROSOL_SIM

      ! Initialize GEOS-Chem tracer array [kg] from Chemistry State object
      ! (mpayer, 12/6/12)
      STT => State_Chm%Tracers

      ! Should we print debug output?
      prtDebug             = ( LPRT .and. am_I_Root )

      ! Get current month
      MONTH                = GET_MONTH()

      ! If it's an offline simulation ...
      IF ( IT_IS_AN_AEROSOL_SIM ) THEN

         ! Then read monthly data files ...
         IF ( ITS_A_NEW_MONTH() ) THEN 
!-----------------------------------------------------------------------------
! Prior to 10/23/13:
!            CALL GET_GLOBAL_OH( MONTH, State_Met )
!-----------------------------------------------------------------------------
            CALL GET_GLOBAL_OH
     &           ( am_I_Root, Input_Opt, State_Met, MONTH, RC )
            CALL GET_GLOBAL_NO3( MONTH )
         ENDIF

         ! And compute time scaling arrays for offline OH, NO3
         CALL OHNO3TIME
         
      ENDIF

      ! Store NTIME in a shadow variable
      ELAPSED_SEC = GET_ELAPSED_SEC()

      ! DTCHEM is the chemistry timestep in seconds
      DTCHEM = GET_TS_CHEM() * 60d0

      ! Initialize module arrays
      PSO2_DMS = 0d0
      PMSA_DMS = 0d0
      PSO4_SO2 = 0d0
      PSO4_SS  = 0d0
      PNITs    = 0d0
#if   defined( TOMAS )
      PSO4_SO2AQ = 0d0     ! For TOMAS microphysics
#endif
                  
      !================================================================= 
      ! Call individual chemistry routines for sulfate/aerosol tracers
      !=================================================================

      ! SO4s [kg] gravitational settling
      CALL GRAV_SETTLING( am_I_Root,          Input_Opt, State_Met,
     &                    STT(:,:,:,IDTSO4s), 1,         RC         )
      IF ( prtDebug ) THEN
         CALL DEBUG_MSG( '### CHEMSULFATE: GRAV_SET, SO4S' )
      ENDIF

      ! NITs [kg] gravitational settling
      CALL GRAV_SETTLING( am_I_Root,          Input_Opt, State_Met, 
     &                    STT(:,:,:,IDTNITs), 2,         RC )
      IF ( prtDebug ) THEN
         CALL DEBUG_MSG( '### CHEMSULFATE: GRAV_SET, NITS' )
      ENDIF

      ! Convert all tracers in STT from [kg] -> [v/v]
      CALL CONVERT_UNITS( 1, N_TRACERS, TCVV, State_Met%AD, STT )
      IF ( prtDebug ) THEN
         CALL DEBUG_MSG( '### CHEMSULFATE: a CONVERT UNITS' )
      ENDIF

      ! For offline runs only ...
      IF ( IT_IS_AN_AEROSOL_SIM ) THEN

         ! DMS (offline only)
         CALL CHEM_DMS( am_I_Root, Input_Opt, State_Met, State_Chm, RC )
         IF ( prtDebug ) THEN
            CALL DEBUG_MSG( '### CHEMSULFATE: a CHEM_DMS' )
         ENDIF

         ! H2O2 (offline only)
         CALL CHEM_H2O2( am_I_Root, Input_Opt, State_Met, State_Chm,
     &                   RC )
         IF ( prtDebug ) THEN
            CALL DEBUG_MSG( '### CHEMSULFATE: a CHEM_H2O2' )
         ENDIF

      ENDIF

      ! SO2
      CALL GET_VCLDF( am_I_Root, State_Met, RC )
      IF ( prtDebug ) THEN
         CALL DEBUG_MSG( '### CHEMSULFATE: a get VCLDF' )
      ENDIF

      CALL CHEM_SO2( am_I_Root, Input_Opt, State_Met, State_Chm, RC )
      IF ( prtDebug ) THEN
         CALL DEBUG_MSG( '### CHEMSULFATE: a CHEM_SO2' )
      ENDIF

      ! SO4
      CALL CHEM_SO4( am_I_Root, Input_Opt, State_Met, State_Chm, RC )
      IF ( prtDebug ) THEN
         CALL DEBUG_MSG( '### CHEMSULFATE: a CHEM_SO4' )
      ENDIF

#if   defined( TOMAS )
      !-----------------------------------------------------------------
      ! For TOMAS microphysics:
      !
      ! SO4 from aqueous chemistry of SO2 (in-cloud oxidation)
      !
      ! SO4 produced via aqueous chemistry is distributed onto 30-bin
      ! aerosol by TOMAS subroutine AQOXID.   NOTE: This may be moved
      ! to tomas_mod.f in the future, but for now it still needs to get
      ! the PSO4_SO2AQ value while CHEMSULFATE is called
      !-----------------------------------------------------------------
      CALL CHEM_SO4_AQ( am_I_Root, Input_Opt, State_Met, State_Chm, RC )
      IF ( LPRT ) CALL DEBUG_MSG( '### CHEMSULFATE: a CHEM_SO4_AQ' )
#endif

      ! MSA
      CALL CHEM_MSA( am_I_Root, Input_Opt, State_Met, State_Chm, RC )
      IF ( prtDebug ) THEN
         CALL DEBUG_MSG( '### CHEMSULFATE: a CHEM_MSA' )
      ENDIF

      ! NH3
      CALL CHEM_NH3( am_I_Root, Input_Opt, State_Met, State_Chm, RC )
      IF ( prtDebug ) THEN
         CALL DEBUG_MSG( '### CHEMSULFATE: a CHEM_NH3' )
      ENDIF

      ! NH4 (gas-phase)
      CALL CHEM_NH4( am_I_Root, Input_Opt, State_Met, State_Chm, RC )
      IF ( prtDebug ) THEN
         CALL DEBUG_MSG( '### CHEMSULFATE: a CHEM_NH4' )
      ENDIF

      ! Sulfur Nitrate
      CALL CHEM_NIT( am_I_Root, Input_Opt, State_Met, State_Chm, RC )
      IF ( prtDebug ) THEN
         CALL DEBUG_MSG( '### CHEMSULFATE: a CHEM_NIT' )
      ENDIF

      ! Convert STT from [v/v] -> [kg]
      CALL CONVERT_UNITS( 2, N_TRACERS, TCVV, State_Met%AD, STT )

      ! Free pointer
      NULLIFY( STT )

      ! We have already gone thru one chemistry iteration
      FIRSTCHEM = .FALSE. 
         
      END SUBROUTINE CHEMSULFATE
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: grav_settling
!
! !DESCRIPTION: Subroutine GRAV\_SETTLING performs gravitational settling of
!  sulfate and nitrate in coarse sea salt (SO4S and NITS).
!  (bec, rjp, bmy, 4/20/04, 7/20/04, 10/25/05)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE GRAV_SETTLING( am_I_Root, Input_Opt, State_Met, 
     &                          TC,        N,         RC         )
!
! !USES:
!
      USE CMN_GCTM_MOD        
      USE CMN_DIAG_MOD        
      USE CMN_SIZE_MOD        
      USE DIAG_MOD,           ONLY : AD44
      USE DRYDEP_MOD,         ONLY : DEPSAV
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE GRID_MOD,           ONLY : GET_AREA_CM2
      USE PRESSURE_MOD,       ONLY : GET_PCENTER
      USE TRACERID_MOD,       ONLY : IDTSO4s
      USE TRACERID_MOD,       ONLY : IDTNITs
      USE TIME_MOD,           ONLY : GET_ELAPSED_SEC
      USE TIME_MOD,           ONLY : GET_TS_CHEM
!
! !INPUT PARAMETERS: 
!
      LOGICAL,        INTENT(IN)    :: am_I_Root   ! Are we on the root CPU?
      TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
      TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology State object
      INTEGER,        INTENT(IN)    :: N           ! N=1 is SO4S; N=2 is NITS
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!
! !INPUT/OUTPUT PARAMETERS: 
!
      REAL*8,         INTENT(INOUT) :: TC(IIPAR,JJPAR,LLPAR) ! Tracer [kg]
! 
! !REVISION HISTORY:
!  (1 ) Now references SALA_REDGE_um and SALC_REDGE_um from "tracer_mod.f"
!        (bmy, 7/20/04)
!  (2 ) Now references XNUMOL from "tracer_mod.f" (bmy, 10/25/05)
!  (3 ) Now limit relative humidity to [tiny(real*8),0.99] range for DLOG
!         argument (phs, 5/1/08)
!  (4 ) Bug fixes to the Gerber hygroscopic growth for sea salt aerosols
!       (jaegle, 5/5/11)
!  (5 ) Update hygroscopic growth to Lewis and Schwartz formulation (2006) and
!       density calculation based on Tang et al. (1997) (bec, jaegle 5/5/11)
!  22 Dec 2011 - M. Payer    - Added ProTeX headers
!  01 Mar 2012 - R. Yantosca - Now use GET_AREA_CM2(I,J,L) from grid_mod.F90
!  14 Nov 2012 - R. Yantosca - Now pass am_I_Root, Input_Opt, RC as arguments
!  15 Nov 2012 - M. Payer    - Replaced all met field arrays with State_Met
!                              derived type object
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
      REAL*8,  PARAMETER     :: C1   =  0.7674d0 
      REAL*8,  PARAMETER     :: C2   =  3.079d0 
      REAL*8,  PARAMETER     :: C3   =  2.573d-11
      REAL*8,  PARAMETER     :: C4   = -1.424d0
      REAL*8,  PARAMETER     :: DEN  = 2200.0d0 ! [kg/m3] sea-salt density

      ! Parameters for polynomial coefficients to derive seawater
      ! density. From Tang et al. (1997) (bec, jaegle, 5/11/11)
      REAL*8,  PARAMETER     :: A1   =  7.93d-3
      REAL*8,  PARAMETER     :: A2   = -4.28d-5
      REAL*8,  PARAMETER     :: A3   =  2.52d-6
      REAL*8,  PARAMETER     :: A4   = -2.35d-8
      REAL*8,  PARAMETER     :: EPSI = 1.0D-4 
!
! !LOCAL VARIABLES:
!
      INTEGER                :: I,      J,     L,        DTCHEM
      REAL*8                 :: DELZ,   DELZ1, REFF
      REAL*8                 :: P,      DP,    PDP,      TEMP        
      REAL*8                 :: CONST,  SLIP,  VISC,     FAC1
      REAL*8                 :: FAC2,   FLUX,  AREA_CM2, RHB
      ! replace RCM (radius in CM with RUM radius in microns) jaegle 5/11/11
      REAL*8                 :: RUM,    RWET,  RATIO_R,  RHO 
      REAL*8                 :: TOT1,   TOT2
      REAL*8                 :: VTS(LLPAR)  
      REAL*8                 :: TC0(LLPAR)
      ! added variables for density calculation (jaegle, bec 5/11/11)
      REAL*8                 :: RHO1, WTP
      
      ! Arrays
      INTEGER                :: IDDEP(2)
      INTEGER                :: IDTRC(2)	
      REAL*8                 :: SALA_REDGE_um(2)
      REAL*8                 :: SALC_REDGE_um(2)
      REAL*8                 :: XNUMOL(Input_Opt%N_TRACERS)

      !=================================================================
      ! GRAV_SETTLING begins here!
      !=================================================================

      ! Return if tracers are undefined
      IF ( IDTSO4s == 0 .and. IDTNITs == 0 ) RETURN

      ! Return if it's the start of the run
      IF ( GET_ELAPSED_SEC() == 0 ) RETURN

      ! Assume success
      RC            = GIGC_SUCCESS

      ! Copy fields from INPUT_OPT to local variables for use below
      SALA_REDGE_um = Input_Opt%SALA_REDGE_um
      SALC_REDGE_um = Input_Opt%SALC_REDGE_um
      XNUMOL        = Input_Opt%XNUMOL

      ! Chemistry timestep [s]
      DTCHEM = GET_TS_CHEM() * 60d0

      ! Store in IDDEP array
      IDDEP(1) = DRYSO4s
      IDDEP(2) = DRYNITs

      ! Tracer array
      IDTRC(1) = IDTSO4s
      IDTRC(2) = IDTNITs

      ! Coarse mode
      REFF = 0.5d-6 * ( SALC_REDGE_um(1) + SALC_REDGE_um(2) )
            
      ! Sea salt radius [cm]
      ! The Gerber formula for hygroscopic growth uses the radius in
      ! micrometers instead of centimeters. This fix is implemented by using
      ! RUM instead of RCM (jaegle 5/5/11)
      RUM  = REFF * 1d6


      ! Exponential factors
      ! replace with radius in microns (jaegle 5/5/11)
      FAC1 = C1 * ( RUM**C2 )
      FAC2 = C3 * ( RUM**C4 )

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I,       J,     L,    VTS,  P,        TEMP, RHB,  RWET ) 
!$OMP+PRIVATE( RATIO_R, RHO,   DP,   PDP,  CONST,    SLIP, VISC, TC0  )
!$OMP+PRIVATE( DELZ,    DELZ1, TOT1, TOT2, AREA_CM2, FLUX             )
!$OMP+PRIVATE( RHO1,    WTP                                           ) 
!$OMP+SCHEDULE( DYNAMIC )
      DO J = 1, JJPAR
      DO I = 1, IIPAR       

         ! Initialize 
         DO L = 1, LLPAR
            VTS(L) = 0d0
         ENDDO

         ! Loop over levels
         DO L = 1, LLPAR

            ! Pressure at center of the level [kPa]
            P       = GET_PCENTER(I,J,L) * 0.1d0

            ! Temperature [K]
            TEMP    = State_Met%T(I,J,L)

            ! Cap RH at 0.99 
            RHB     = MIN( 0.99d0, State_Met%RH(I,J,L) * 1d-2 )

            ! Safety check (phs, 5/1/08)
            RHB     = MAX( TINY(RHB), RHB           )

            ! Aerosol growth with relative humidity in radius [m] 
            ! (Gerber, 1985)
            ! Several bug fixes to the Gerber formulation: a log10 (instead of
            ! ln) should be used and the dry radius should be expressed in
            ! micrometers (instead of cm) also add more significant digits to
            ! the exponent (jaegle 5/5/11)
            !RWET    = 1d-6*(FAC1/(FAC2-LOG10(RHB))+RUM**3.d0)**0.33333d0

            ! Use equation 5 in Lewis and Schwartz (2006) [m] for sea salt
            ! growth (jaegle 5/11/11)
            RWET = REFF * (4.d0 / 3.7d0) *
     &                  ( (2.d0 - RHB)/(1.d0 - RHB) )**(1.d0/3.d0)


            ! Ratio dry over wet radii at the cubic power
            RATIO_R = ( REFF / RWET )**3.d0

            ! Density of the wet aerosol (kg/m3)
            !RHO     = RATIO_R * DEN + ( 1.d0 - RATIO_R ) * 1000.d0

            ! Above density calculation is chemically unsound because it
            ! ignores chemical solvation.  
            ! Iteratively solve Tang et al., 1997 equation 5 to calculate
            ! density of wet aerosol (kg/m3) 
            ! (bec, jaegle 5/11/11)
            RATIO_R = ( REFF / RWET )
            ! Assume an initial density of 1000 kg/m3
            RHO  = 1000.D0
            RHO1 = 0.d0 !initialize (bec, 6/21/10)
            DO WHILE ( ABS( RHO1-RHO ) .gt. EPSI )
                ! First calculate weight percent of aerosol (kg_RH=0.8/kg_wet) 
                WTP    = 100.d0 * DEN/RHO * RATIO_R**3.d0
                ! Then calculate density of wet aerosol using equation 5 
                ! in Tang et al., 1997 [kg/m3]
                RHO1   = ( 0.9971d0 + (A1 * WTP) + (A2 * WTP**2.d0) + 
     $               (A3 * WTP**3.d0) + (A4 * WTP**4.d0) ) * 1000.d0
                ! Now calculate new weight percent using above density
                ! calculation
                WTP    = 100.d0 * DEN/RHO1 * RATIO_R**3.d0
                ! Now recalculate new wet density [kg/m3]
                RHO   = ( 0.9971d0 + (A1 * WTP) + (A2 * WTP**2.d0) + 
     $              (A3 * WTP**3.d0) + (A4 * WTP**4.d0) ) * 1000.d0
            ENDDO

            ! Dp = particle diameter [um]
            DP      = 2.d0 * RWET * 1.d6        

            ! PdP = P * dP [hPa * um]
            PDp     = P * Dp

            ! Constant
            CONST   = 2.d0 * RHO * RWET**2 * g0 / 9.d0

            !===========================================================
            ! NOTE: Slip correction factor calculations following 
            ! Seinfeld, pp464 which is thought to be more accurate 
            ! but more computation required. (rjp, 1/24/02)
            !
            ! # air molecule number density
            ! num = P * 1d3 * 6.023d23 / (8.314 * Temp) 
            !
            ! # gas mean free path
            ! lamda = 1.d6/( 1.41421 * num * 3.141592 * (3.7d-10)**2 ) 
            !
            ! # Slip correction
            ! Slip = 1. + 2. * lamda * (1.257 + 0.4 * exp( -1.1 * Dp     
            !     &     / (2. * lamda))) / Dp
            !
            ! NOTE: Eq) 3.22 pp 50 in Hinds (Aerosol Technology)
            ! which produces slip correction factore with small error
            ! compared to the above with less computation.
            !===========================================================  
          
            ! Slip correction factor (as function of P*dp)
            Slip = 1.d0+(15.60d0 + 7.0d0 * EXP(-0.059d0 * PDp)) / PDp

            ! Viscosity [Pa*s] of air as a function of temperature 
            VISC = 1.458d-6 * (Temp)**(1.5d0) / ( Temp + 110.4d0 )

            ! Settling velocity [m/s]
            VTS(L) = CONST * Slip / VISC
         ENDDO

         ! Method is to solve bidiagonal matrix which is
         ! implicit and first order accurate in z (rjp, 1/24/02)

         ! Save initial tracer concentration in column
         DO L = 1, LLPAR
            TC0(L) = TC(I,J,L)
         ENDDO

         ! We know the boundary condition at the model top
         L    = LLTROP
         DELZ = State_Met%BXHEIGHT(I,J,L)

         TC(I,J,L) = TC(I,J,L) / ( 1.d0 + DTCHEM * VTS(L) / DELZ )

         DO L = LLTROP-1, 1, -1
            DELZ  = State_Met%BXHEIGHT(I,J,L)
            DELZ1 = State_Met%BXHEIGHT(I,J,L+1)
            TC(I,J,L) = 1.d0 / ( 1.d0 + DTCHEM * VTS(L) / DELZ )
     &                * ( TC(I,J,L) + DTCHEM * VTS(L+1) / DELZ1
     &                *  TC(I,J,L+1) )
         ENDDO
         
         !==============================================================
         ! ND44 diagnostic: sea salt loss [molec/cm2/s]
         !==============================================================
         IF ( ND44 > 0 ) THEN

            ! Initialize
            TOT1 = 0d0
            TOT2 = 0d0
            
            ! Compute column totals of TCO(:) and TC(I,J,:,N)
            DO L = 1, LLPAR
               TOT1 = TOT1 + TC0(L)
               TOT2 = TOT2 + TC(I,J,L)
            ENDDO

            ! Surface area [cm2]
            AREA_CM2 = GET_AREA_CM2( I, J, 1 )

            ! Convert sea salt flux from [kg/s] to [molec/cm2/s]
            FLUX     = ( TOT1 - TOT2 ) / DTCHEM
            FLUX     = FLUX * XNUMOL(IDTRC(N)) / AREA_CM2 
   
            ! Store in AD44 array
            AD44(I,J,IDDEP(N),1) = AD44(I,J,IDDEP(N),1) + FLUX
         ENDIF
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      END SUBROUTINE GRAV_SETTLING
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: chem_dms
!
! !DESCRIPTION: Subroutine CHEM\_DMS is the DMS chemistry subroutine from Mian
!  Chin's GOCART model, modified for use with the GEOS-CHEM model.
!  (rjp, bdf, bmy, 5/31/00, 10/15/09)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CHEM_DMS( am_I_Root, Input_Opt, 
     &                     State_Met, State_Chm, RC )
!
! !USES:
!
      USE CMN_GCTM_MOD         
      USE CMN_DIAG_MOD         
      USE CMN_SIZE_MOD         
      USE DIAG_MOD,           ONLY : AD05
      USE DRYDEP_MOD,         ONLY : DEPSAV
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE GIGC_State_Chm_Mod, ONLY : ChmState
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE TIME_MOD,           ONLY : GET_TS_CHEM
      USE TRACERID_MOD,       ONLY : IDTDMS
      USE TROPOPAUSE_MOD,     ONLY : ITS_IN_THE_STRAT
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
!  Reaction List (by Mian Chin, chin@rondo.gsfc.nasa.gov)                  
!  ============================================================================
!                                                                             .
!  R1:    DMS + OH  -> a*SO2 + b*MSA                OH addition channel    
!         k1 = { 1.7e-42*exp(7810/T)*[O2] / (1+5.5e-31*exp(7460/T)*[O2] }  
!         a = 0.75, b = 0.25                                               
!                                                                             .
!  R2:    DMS + OH  ->   SO2 + ...                  OH abstraction channel 
!         k2 = 1.2e-11*exp(-260/T)                                         
!                                                                             .
!         DMS_OH = DMS0 * exp(-(r1+r2)* NDT1)                                  
!         where DMS0 is the DMS concentration at the beginning,            
!         r1 = k1*[OH], r2 = k2*[OH].                                      
!                                                                             .
!  R3:    DMS + NO3 ->   SO2 + ...                                         
!         k3 = 1.9e-13*exp(500/T)                                          
!                                                                             .
!         DMS = DMS_OH * exp(-r3*NDT1)                                         
!         where r3 = k3*[NO3].                                             
!                                                                             .
!  R4:    DMS + X   ->   SO2 + ...                                         
!         assume to be at the rate of DMS+OH and DMS+NO3 combined.         
!                                                                             .
!  The production of SO2 and MSA here, PSO2_DMS and PMSA_DMS, are saved    
!  for use in CHEM_SO2 and CHEM_MSA subroutines as a source term.  They    
!  are in unit of [v/v/timestep]. 
!
! !REVISION HISTORY:
!  (1 ) Now reference AD, AIRDEN, and SUNCOS from "dao_mod.f".  Added 
!        parallel DO-loops.  Also now extract OH and NO3 from SMVGEAR
!        for coupled chemistry-aerosol runs. (rjp, bdf, bmy, 9/16/02)
!  (2 ) Bug fix: remove duplicate definition of RK3 (bmy, 3/23/03)
!  (3 ) Now use function GET_TS_CHEM from "time_mod.f".  (bmy, 3/27/03)
!  (4 ) Now reference STT and ITS_A_FULLCHEM_SIM from "tracer_mod.f"
!        Now replace IJSURF w/ an analytic function. (bmy, 7/20/04)
!  (5 ) Shift rows 8,9 in AD05 to 9,10 in to make room for P(SO4) from O3 
!        oxidation in sea-salt aerosols (bec, bmy, 4/13/05)
!  (6 ) Now remove reference to CMN, it's obsolete.  Now reference 
!        ITS_IN_THE_STRAT from "tropopause_mod.f". (bmy, 8/22/05)
!  (7 ) Now references XNUMOL from "tracer_mod.f" (bmy, 10/25/05)
!  (8 ) Now correctly records P(SO2) from OH in AD05 (pjh)
!  (9 ) Update reaction rate to match JPL06 and full chem (jaf, bmy, 10/15/09)
!  22 Dec 2011 - M. Payer    - Added ProTeX headers
!  31 Jul 2012 - R. Yantosca - Now loop from 1..LLPAR for GIGC compatibility
!  14 Nov 2012 - R. Yantosca - Add am_I_Root, Input_Opt, RC as arguments
!  15 Nov 2012 - M. Payer    - Replaced all met field arrays with State_Met
!                              derived type object
!  28 Nov 2012 - R. Yantosca - Replace SUNCOS with State_Met%SUNCOS
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
      REAL*8, PARAMETER :: FX = 1.0d0
      REAL*8, PARAMETER :: A  = 0.75d0
      REAL*8, PARAMETER :: B  = 0.25d0

      ! From D4: only 0.8 efficiency, also some goes to DMSO and lost.  
      ! So we assume 0.75 efficiency for DMS addtion channel to form     
      ! products.                                                        
      REAL*8, PARAMETER :: EFF = 1d0
!
! !LOCAL VARIABLES:
!
      ! Scalars
      INTEGER           :: I,   J,    L
      REAL*8            :: TK,  O2,   RK1,    RK2,    RK3,   F  
      REAL*8            :: DMS, DMS0, DMS_OH, DTCHEM, XOH,   XN3 
      REAL*8            :: XX,  OH,   OH0,    XNO3,   XNO30, LOH
      REAL*8            :: LNO3
                        
      ! For fields from Input_Opt
      LOGICAL           :: IS_FULLCHEM
      REAL*8            :: XNUMOL(Input_Opt%N_TRACERS)

      ! Pointers
      ! We need to define local arrays to hold corresponding values 
      ! from the Chemistry State (State_Chm) object. (mpayer, 12/6/12)
      REAL*8, POINTER :: STT(:,:,:,:)
!
! !EXTERNAL FUNCTIONS:   
!
      REAL*8,  EXTERNAL      :: BOXVL
      
      !=================================================================
      ! CHEM_DMS begins here!
      !=================================================================
      IF ( IDTDMS == 0 ) RETURN

      ! Assume success
      RC          = GIGC_SUCCESS

      ! Copy fields from INPUT_OPT to local variables for use below
      IS_FULLCHEM = Input_Opt%ITS_A_FULLCHEM_SIM
      XNUMOL      = Input_Opt%XNUMOL

      ! Initialize GEOS-Chem tracer array [kg] from Chemistry State object
      ! (mpayer, 12/6/12)
      STT => State_Chm%Tracers

      ! DTCHEM is the chemistry timestep in seconds
      DTCHEM      = GET_TS_CHEM() * 60d0

      ! Factor to convert AIRDEN from kgair/m3 to molecules/cm3:
      f           = 1000.d0 / AIRMW * 6.022d23 * 1.d-6
      
      !=================================================================
      ! Do the chemistry over all tropospheric grid boxes!
      !=================================================================
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, TK, O2, DMS0,OH, XNO3, RK1, RK2 )
!$OMP+PRIVATE( RK3, DMS_OH, DMS, OH0, XNO30, XOH, XN3, XX, LOH, LNO3  )
!$OMP+SCHEDULE( DYNAMIC )
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Skip stratospheric boxes
         IF ( ITS_IN_THE_STRAT( I, J, L, State_Met ) ) CYCLE

         ! Temperature [K]
         TK     = State_Met%T(I,J,L)

         ! Get O2 [molec/cm3], DMS [v/v], OH [molec/cm3], NO3 [molec/cm3]
         O2     = State_Met%AIRDEN(L,I,J) * f * 0.21d0
         DMS0   = STT(I,J,L,IDTDMS)
         OH     = GET_OH(  I, J, L, Input_Opt, State_Met )
         XNO3   = GET_NO3( I, J, L, Input_Opt, State_Met )

         !==============================================================
         ! (1) DMS + OH:  RK1 - addition channel  
         !                RK2 - abstraction channel   
         !==============================================================
         RK1 = 0.d0
         RK2 = 0.d0
         RK3 = 0.d0

         IF ( OH > 0.d0 ) THEN
            RK1 = ( 1.7d-42 * EXP( 7810.d0 / TK ) * O2 ) /
     &            ( 1.d0 + 5.5d-31 * EXP( 7460.d0 / TK ) * O2 ) * OH

            ! Update reaction rate to match JPL06 and full chem
            ! (jaf, bmy, 10/15/09)
            RK2 = 1.1d-11 * EXP( -240.d0 / TK ) * OH 
         ENDIF
            
         !==============================================================
         ! (2) DMS + NO3 (only happens at night):  
         !==============================================================
         IF ( State_Met%SUNCOS(I,J) <= 0d0 ) THEN
            RK3 = 1.9d-13 * EXP( 500.d0 / TK ) * XNO3
         ENDIF

         !==============================================================
         ! Update DMS concentrations after reaction with OH and NO3, 
         ! and also account for DMS + X assuming at a rate as 
         ! (DMS+OH)*Fx in the day and (DMS+NO3)*Fx at night:   
         ! 
         ! DMS_OH :  DMS concentration after reaction with OH  
         ! DMS    :  DMS concentration after reaction with NO3       
         !           (min(DMS) = 1.0E-32)       
         !
         ! NOTE: If we are doing a coupled fullchem/aerosol run, then
         ! also modify OH and NO3 concentrations after rxn w/ DMS.
         !==============================================================
         DMS_OH = DMS0   * EXP( -( RK1 + RK2 ) * Fx * DTCHEM )
         DMS    = DMS_OH * EXP( -( RK3       ) * Fx * DTCHEM ) 
         IF ( DMS < SMALLNUM ) DMS = 0d0

         ! Archive initial OH and NO3 for diagnostics
         OH0    = OH
         XNO30  = XNO3

         IF ( IS_FULLCHEM ) THEN
         
            ! Update OH after rxn w/ DMS (coupled runs only)
            OH    = OH0 - ( ( DMS0 - DMS_OH ) *
     &              State_Met%AIRDEN(L,I,J) * f )
            IF ( OH < SMALLNUM ) OH = 0d0

            ! Update NO3 after rxn w/ DMS (coupled runs only)
            XNO3  = XNO30 - ( ( DMS_OH - DMS ) *
     &              State_Met%AIRDEN(L,I,J) * f )
            IF ( XNO3 < SMALLNUM ) XNO3 = 0d0

         ENDIF 

         ! Save DMS back to the tracer array
         STT(I,J,L,IDTDMS) = DMS

         !==============================================================
         ! Save SO2 and MSA production from DMS oxidation 
         ! in [mixing ratio/timestep]:    
         !
         ! SO2 is formed in DMS+OH addition (0.85) and abstraction 
         ! (1.0) channels as well as DMS + NO3 reaction.  We also 
         ! assume that SO2 yield from DMS + X is 1.0.  
         !
         ! MSA is formed in DMS + OH addition (0.15) channel. 
         !==============================================================
         IF ( ( RK1 + RK2 ) == 0.d0 ) THEN
            PMSA_DMS(I,J,L) = 0.d0
         ELSE
            PMSA_DMS(I,J,L) = ( DMS0 - DMS_OH ) * 
     &                          B*RK1 / ( ( RK1 + RK2 ) * Fx ) * EFF
         ENDIF

         PSO2_DMS(I,J,L) =  DMS0 - DMS - PMSA_DMS(I,J,L)

         !==============================================================
         ! ND05 diagnostic: production and loss  
         !
         ! For the offline run, we are reading in monthly mean OH, NO3 
         ! from disk.  We don't modify these, so LOH = 0 and LNO3 = 0.
         !==============================================================
         IF ( ND05 > 0 .and. L <= LD05 ) THEN

            ! P(SO2) from DMS+OH, DMS+NO3, and DMS+X
            XOH  = ( DMS0   - DMS_OH - PMSA_DMS(I,J,L) ) /
     &                 Fx   * State_Met%AD(I,J,L) / TCVV_S
            XN3  = ( DMS_OH - DMS ) / Fx * State_Met%AD(I,J,L) / TCVV_S
            XX   = ( ( DMS0 - DMS ) * State_Met%AD(I,J,L) / TCVV_S )
     &               - XOH  - XN3
        
            ! Convert L(OH) and L(NO3) from [molec/cm3] to [kg/timestep]
            LOH  = ( OH0   - OH  ) * BOXVL(I,J,L,State_Met) / XNUMOL_OH
            LNO3 = ( XNO30 - XNO3) * BOXVL(I,J,L,State_Met) / XNUMOL_NO3 

            ! Store P(SO2) from DMS + OH [kg S/timestep]
            AD05(I,J,L,1) = AD05(I,J,L,1) + XOH

            ! Store P(SO2) from DMS + NO3 [kg S/timestep]
            AD05(I,J,L,2) = AD05(I,J,L,2) + XN3

            ! Store total P(SO2) from DMS [kg S/timestep]
            AD05(I,J,L,3) = AD05(I,J,L,3)   + 
     &                    ( PSO2_DMS(I,J,L) * State_Met%AD(I,J,L) /
     &                      TCVV_S )

            ! Store P(MSA) from DMS [kg S/timestep]
            AD05(I,J,L,4) = AD05(I,J,L,4)   + 
     &                    ( PMSA_DMS(I,J,L) * State_Met%AD(I,J,L) /
     &                      TCVV_S )

            ! Store L(OH) by DMS [kg OH/timestep]
            AD05(I,J,L,9) = AD05(I,J,L,9) + LOH
            
            ! Store L(NO3) by DMS [kg NO3/timestep]
            AD05(I,J,L,10) = AD05(I,J,L,10) + LNO3

         ENDIF

         !==============================================================
         ! For a coupled fullchem/aerosol run, save OH [molec/cm3] 
         ! and NO3 [molec/cm3] back into the CSPEC array of SMVGEAR
         !==============================================================
         IF ( IS_FULLCHEM ) THEN
            CALL SET_OH( I, J, L, OH )
            CALL SET_NO3( I, J, L, XNO3 )
         ENDIF
      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! Free pointer
      NULLIFY( STT )

      END SUBROUTINE CHEM_DMS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: chem_h2o2
!
! !DESCRIPTION: Subroutine CHEM\_H2O2 is the H2O2 chemistry subroutine for
!  offline sulfate simulations.  For coupled runs, H2O2 chemistry is already
!  computed by the SMVGEAR module. (rjp, bmy, 11/26/02, 10/25/05)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CHEM_H2O2( am_I_Root, Input_Opt, 
     &                      State_Met, State_Chm, RC )
!
! !USES:
!
      USE BPCH2_MOD,          ONLY : GET_NAME_EXT
      USE BPCH2_MOD,          ONLY : GET_RES_EXT
      USE BPCH2_MOD,          ONLY : GET_TAU0
      USE BPCH2_MOD,          ONLY : READ_BPCH2
      USE CMN_SIZE_MOD
      USE CMN_FJ_MOD,         ONLY : JPMAX, JPPJ
      USE CMN_DIAG_MOD
      USE CMN_GCTM_MOD
      USE DIAG_MOD,           ONLY : AD44 
      USE DIRECTORY_MOD,      ONLY : DATA_DIR
      USE DRYDEP_MOD,         ONLY : DEPSAV
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE GIGC_State_Chm_Mod, ONLY : ChmState
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE GRID_MOD,           ONLY : GET_AREA_CM2
      USE PBL_MIX_MOD,        ONLY : GET_FRAC_UNDER_PBLTOP
      USE TIME_MOD,           ONLY : GET_MONTH
      USE TIME_MOD,           ONLY : GET_TS_CHEM
      USE TIME_MOD,           ONLY : ITS_A_NEW_MONTH
      USE TRACERID_MOD,       ONLY : IDTH2O2
      USE TRANSFER_MOD,       ONLY : TRANSFER_3D_TROP
      USE TROPOPAUSE_MOD,     ONLY : ITS_IN_THE_STRAT
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)    :: am_I_Root   ! Is this the root CPU?
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
!  (1 ) Bug fix: need to multiply DXYP by 1d4 for cm2 (bmy, 3/23/03)
!  (2 ) Now replace DXYP(JREF)*1d4 with routine GET_AREA_CM2 of "grid_mod.f"
!        Now use functions GET_MONTH and GET_TS_CHEM from "time_mod.f".
!        (bmy, 3/27/03)
!  (3 ) Now references PBLFRAC from "drydep_mod.f".  Now apply dry deposition 
!        throughout the entire PBL.  Added FREQ variable. (bmy, 8/1/03)
!  (4 ) Now use ND44_TMP array to store vertical levels of drydep flux, then 
!        sum into AD44 array.  This preents numerical differences when using
!        multiple processors. (bmy, 3/24/04)    
!  (5 ) Now use diurnally-varying JO1D.  Now use new unit conversion for
!        the ND44 diagnostic. (rjp, bmy, 3/30/04)
!  (6 ) Now use parallel DO-loop to zero ND44_TMP.  Now uses ITS_A_NEW_MONTH
!        from time_mod.f. (bmy, 4/14/04)
!  (7 ) Now reference STT & TCVV from "tracer_mod.f".  Also replace IJSURF
!        with an analytic function.  Now references DATA_DIR from 
!        "directory_mod.f". (bmy, 7/20/04)
!  (8 ) Now suppress output from READ_BPCH with QUIET keyword (bmy, 1/25/05)
!  (9 ) Replace PBLFRAC from "drydep_mod.f" with GET_FRAC_UNDER_PBLTOP
!        from "pbl_mix_mod.f" (bmy, 2/22/05)
!  (10) Now read offline files from "sulfate_sim_200508/offline".  Now remove 
!        reference to CMN, it's obsolete.  Now reference ITS_IN_THE_STRAT from 
!        "tropopause_mod.f". (bmy, 8/22/05)
!  (11) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (12) Now references XNUMOL from "tracer_mod.f" (bmy, 10/25/05)
!  22 Dec 2011 - M. Payer    - Added ProTeX headers
!  01 Mar 2012 - R. Yantosca - Now use GET_AREA_CM2(I,J,L) from grid_mod.F90
!  30 Jul 2012 - R. Yantosca - Now accept am_I_Root as an argument when
!                              running with the traditional driver main.F
!  31 Jul 2012 - R. Yantosca - Now loop from 1..LLPAR for GIGC compatibility
!  31 Jul 2012 - R. Yantosca - Declare temp drydep arrays w/ LLPAR (not LLTROP)
!  14 Nov 2012 - R. Yantosca - Add am_I_Root, Input_Opt, RC as arguments
!  15 Nov 2012 - M. Payer    - Replaced all met field arrays with State_Met
!                              derived type object
!  26 Nov 2012 - R. Yantosca - Dimension ND44_TMP array with LLPAR, not LLTROP
!  28 Nov 2012 - R. Yantosca - Replace SUNCOS with State_Met%SUNCOS
!  05 Mar 2013 - R. Yantosca - Now use Input_Opt%LNLPBL
!  19 Mar 2013 - R. Yantosca - Now copy Input_Opt%TCVV(1:N_TRACERS) and
!                              Input_Opt%XNUMOL(1:N_TRACERS) -- avoid OOB errs
!  25 Mar 2013 - M. Payer    - Now pass State_Chm object via the arg list
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
! 
      REAL*8,  PARAMETER :: A = 2.9d-12
!
! !LOCAL VARIABLES:
!
      ! Scalars
      LOGICAL            :: LNLPBL
      LOGICAL            :: FIRST     = .TRUE.
      INTEGER, SAVE      :: LASTMONTH = -99
      INTEGER            :: I, J, L, JLOOP, N_TRACERS
      REAL*4             :: ARRAY(IIPAR,JJPAR,LLTROP)
      REAL*8             :: DT,    Koh,  DH2O2, M,    F ,   XTAU   
      REAL*8             :: H2O20, H2O2, ALPHA, FLUX, FREQ, PHOTJ
      CHARACTER(LEN=255) :: FILENAME

      ! Arrays
      REAL*8             :: TCVV  (Input_Opt%N_TRACERS)
      REAL*8             :: XNUMOL(Input_Opt%N_TRACERS)
      REAL*8             :: ND44_TMP(IIPAR,JJPAR,LLPAR)

      ! Pointers
      ! We need to define local arrays to hold corresponding values 
      ! from the Chemistry State (State_Chm) object. (mpayer, 12/6/12)
      REAL*8, POINTER :: STT(:,:,:,:)

      !=================================================================
      ! CHEM_H2O2 begins here!
      !=================================================================
      IF ( IDTH2O2 == 0 .or. DRYH2O2 == 0 ) RETURN 

      ! Assume success
      RC        = GIGC_SUCCESS

      ! Copy fields from INPUT_OPT to local variables for use below
      N_TRACERS = Input_Opt%N_TRACERS
      TCVV      = Input_Opt%TCVV  (1:N_TRACERS)
      XNUMOL    = Input_Opt%XNUMOL(1:N_TRACERS)
      LNLPBL    = Input_Opt%LNLPBL

      ! Initialize GEOS-Chem tracer array [kg] from Chemistry State object
      ! (mpayer, 12/6/12)
      STT    => State_Chm%Tracers

      ! Chemistry timestep [s]
      DT        = GET_TS_CHEM() * 60d0

      ! Factor to convert AIRDEN from kgair/m3 to molecules/cm3:
      F         = 1000.d0 / AIRMW * 6.022d23 * 1.d-6
      
      ! Zero ND44_TMP array
      IF ( ND44 > 0 ) THEN
         ND44_TMP = 0d0
      ENDIF

      !=================================================================
      ! For offline run: read J(H2O2) from disk below
      !=================================================================
      IF ( ITS_A_NEW_MONTH() ) THEN

         ! File name to read data 
         FILENAME = TRIM( DATA_DIR )                       // 
     &              'sulfate_sim_200508/offline/JH2O2.'    // 
     &              GET_NAME_EXT() // '.' // GET_RES_EXT()
           
         ! Print filename
         IF ( am_I_Root ) THEN
            WRITE( 6, 100 ) TRIM( FILENAME )
         ENDIF
 100     FORMAT( '     - CHEM_H2O2: Reading ', a )

         ! Get TAU0 value for this month in "generic" year 1985
         XTAU = GET_TAU0( GET_MONTH(), 1, 1985 )

         ! Read J(H2O2) [s-1]  from disk (only up to tropopause)
         ! limit array 3d dimension to LLTROP_FIX, i.e, case of annual mean
         ! tropopause. This is backward compatibility with 
         ! offline data set.
         CALL READ_BPCH2( FILENAME, 'JV-MAP-$', 3,      
     &     XTAU,        IIPAR,                    JJPAR,      
     &     LLTROP_FIX,  ARRAY(:,:,1:LLTROP_FIX),  QUIET=.TRUE. )

         ! Cast to REAL*8 and resize if necessary
         CALL TRANSFER_3D_TROP( ARRAY, JH2O2 )
            
         ! Reset LASTMONTH
         !LASTMONTH = GET_MONTH()
      ENDIF

      !=================================================================
      ! Loop over tropopsheric grid boxes and do chemistry
      !=================================================================
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, M, H2O20, KOH, FREQ, ALPHA, DH2O2, H2O2, FLUX )
!$OMP+PRIVATE( JLOOP, PHOTJ )
!$OMP+SCHEDULE( DYNAMIC )
      DO L  = 1, LLPAR
      DO J  = 1, JJPAR
      DO I  = 1, IIPAR

         ! Initialize for safety's sake 
         FLUX = 0d0
         FREQ = 0d0

         ! Skip stratospheric boxes
         IF ( ITS_IN_THE_STRAT( I, J, L, State_Met ) ) CYCLE

         ! Density of air [molec/cm3]
         M     = State_Met%AIRDEN(L,I,J) * f  

         ! Initial H2O2 [v/v]
         H2O20 = STT(I,J,L,IDTH2O2)

         ! Loss frequenty due to OH oxidation [s-1]
         KOH   = A * EXP( -160.d0 / State_Met%T(I,J,L) ) * 
     &           GET_OH( I, J, L, Input_Opt, State_Met )

         ! H2O2 drydep frequency [1/s].  Account for the fraction
         ! of grid box (I,J,L) that is located beneath the PBL top.
         FREQ  = DEPSAV(I,J,DRYH2O2) * GET_FRAC_UNDER_PBLTOP( I, J, L ) 

         ! Add option for non-local PBL (Lin, 03/31/09)
         IF ( LNLPBL ) FREQ = 0.D0

         ! 1-D grid box index for SUNCOS
         JLOOP = ( (J-1) * IIPAR ) + I

         ! Impose a diurnal variation of jH2O2 by multiplying COS of 
         ! solar zenith angle normalized by maximum solar zenith angle 
         ! because the archived JH2O2 is for local noon time
         IF ( COSZM(I,J) > 0.d0 ) THEN
            PHOTJ = JH2O2(I,J,L) * State_Met%SUNCOS(I,J) / COSZM(I,J)
            PHOTJ = MAX( PHOTJ, 0d0 )
         ELSE
            PHOTJ = 0d0
         ENDIF

         ! Compute loss fraction from OH, photolysis, drydep [unitless].  
         ALPHA = 1.D0 + ( KOH + PHOTJ + FREQ ) * DT 

         ! Delta H2O2 [v/v]
         DH2O2 = PH2O2m(I,J,L) * DT / ( ALPHA * M )
         
         ! Final H2O2 [v/v]
         H2O2  = ( H2O20 / ALPHA + DH2O2 )
         IF ( H2O2 < SMALLNUM ) H2O2 = 0d0

         ! Store final H2O2 in STT
         STT(I,J,L,IDTH2O2) = H2O2

         !==============================================================
         ! ND44 diagnostics: H2O2 drydep loss [molec/cm2/s]
         !==============================================================
         IF ( ND44 > 0 .AND. FREQ > 0d0 ) THEN

            ! Convert H2O2 from [v/v] to H2O2 [molec/cm2/s]
            FLUX = H2O20 * FREQ * DT / ( 1.D0 + FREQ * DT )
            FLUX = FLUX * State_Met%AD(I,J,L) / TCVV(IDTH2O2)
            FLUX = FLUX * XNUMOL(IDTH2O2) 
     &            / ( GET_AREA_CM2( I, J, L ) * DT )

            ! Save dryd flx in ND44_TMP as a placeholder
            ND44_TMP(I,J,L) = ND44_TMP(I,J,L) + FLUX
         ENDIF
      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO
      
      !===============================================================
      ! ND44: Sum drydep fluxes by level into the AD44 array in
      ! order to ensure that  we get the same results w/ sp or mp 
      !===============================================================
      IF ( ND44 > 0 ) THEN 
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
         DO J = 1, JJPAR
         DO I = 1, IIPAR
         DO L = 1, LLPAR
            AD44(I,J,DRYH2O2,1) = AD44(I,J,DRYH2O2,1) + ND44_TMP(I,J,L)
         ENDDO
         ENDDO
         ENDDO
!$OMP END PARALLEL DO
      ENDIF

      ! Free pointer
      NULLIFY( STT )

      END SUBROUTINE CHEM_H2O2
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: chem_so2
!
! !DESCRIPTION: Subroutine CHEM\_SO2 is the SO2 chemistry subroutine. 
!  (rjp, bmy, 11/26/02, 8/26/10) 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CHEM_SO2( am_I_Root, Input_Opt, 
     &                     State_Met, State_Chm, RC )
!
! !USES:
!
      USE CMN_GCTM_MOD
      USE CMN_DIAG_MOD
      USE CMN_SIZE_MOD
      USE DAO_MOD,            ONLY : IS_WATER
      USE DIAG_MOD,           ONLY : AD05
      USE DIAG_MOD,           ONLY : AD44
      USE DRYDEP_MOD,         ONLY : DEPSAV
      USE DIRECTORY_MOD,      ONLY : DATA_DIR
      USE ERROR_MOD,          ONLY : IS_SAFE_EXP
      USE ERROR_MOD,          ONLY : SAFE_DIV
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE GIGC_State_Chm_Mod, ONLY : ChmState
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE GLOBAL_HNO3_MOD,    ONLY : GET_GLOBAL_HNO3
      USE GLOBAL_HNO3_MOD,    ONLY : GET_HNO3_VV
      USE GRID_MOD,           ONLY : GET_AREA_CM2
      USE PBL_MIX_MOD,        ONLY : GET_FRAC_UNDER_PBLTOP
      USE PRESSURE_MOD,       ONLY : GET_PCENTER
      USE TIME_MOD,           ONLY : GET_TS_CHEM, GET_MONTH
      USE TIME_MOD,           ONLY : ITS_A_NEW_MONTH
      USE TRACERID_MOD,       ONLY : IDTH2O2
      USE TRACERID_MOD,       ONLY : IDTSO2
      USE TRACERID_MOD,       ONLY : IDTSO4,  IDTNH3,  IDTNH4,  IDTHNO3
      USE TRACERID_MOD,       ONLY : IDTNIT,  IDTDST1, IDTDST2, IDTDST3
      USE TRACERID_MOD,       ONLY : IDTDST4, IDTSALA, IDTSALC
      USE SEASALT_MOD,        ONLY : GET_ALK
      USE WETSCAV_MOD,        ONLY : H2O2s
      USE WETSCAV_MOD,        ONLY : SO2s
      USE TROPOPAUSE_MOD,     ONLY : ITS_IN_THE_STRAT
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)    :: am_I_Root   ! Is this the root CPU?
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
!  Reaction List (by Rokjin Park, rjp@io.harvard.edu)                      
!  ============================================================================
!  (1 ) SO2 production:                                                      
!       DMS + OH, DMS + NO3 (saved in CHEM_DMS)                               
!                                                                             .
!  (2 ) SO2 loss:                                                         
!       (a) SO2 + OH  -> SO4                                               
!       (b) SO2       -> drydep                                             
!       (c) SO2 + H2O2 or O3 (aq) -> SO4                         
!                                                                             .
!  (3 ) SO2 = SO2_0 * exp(-bt) +  PSO2_DMS/bt * [1-exp(-bt)]   
!                                                                             .
!       where b is the sum of the reaction rate of SO2 + OH and the dry       
!       deposition rate of SO2, PSO2_DMS is SO2 production from DMS in        
!       MixingRatio/timestep.                                                 
!                                                                             .
!  If there is cloud in the gridbox (fraction = fc), then the aqueous      
!  phase chemistry also takes place in cloud. The amount of SO2 oxidized   
!  by H2O2 in cloud is limited by the available H2O2; the rest may be      
!  oxidized due to additional chemistry, e.g, reaction with O3 or O2       
!  (catalyzed by trace metal).                                             
!                                                                          
! !REVISION HISTORY: 
!  (1 ) Removed duplicate definition of Ki (bmy, 11/15/01)     
!  (2 ) Eliminate duplicate HPLUS definition.  Make adjustments to facilitate 
!        SMVGEAR chemistry for fullchem runs (rjp, bmy, 3/23/03)
!  (3 ) Now replace DXYP(J+J0)*1d4 with routine GET_AREA_CM2 of "grid_mod.f"
!        Now use function GET_TS_CHEM from "time_mod.f".
!  (4 ) Now apply dry deposition to entire PBL.  Now references PBLFRAC array
!        from "drydep_mod.f". (bmy, 8/1/03)  
!  (5 ) Now use ND44_TMP array to store vertical levels of drydep flux, then 
!        sum into AD44 array.  This preents numerical differences when using
!        multiple processors. (bmy, 3/24/04)
!  (6 ) Now use parallel DO-loop to zero ND44_TMP (bmy, 4/14/04)
!  (7 ) Now reference STT, TCVV, & ITS_AN_AEROSOL_SIM from "tracer_mod.f".
!        Now reference DATA_DIR from "directory_mod.f" (bmy, 7/20/04)
!  (8 ) Replace PBLFRAC from "drydep_mod.f" with GET_FRAC_UNDER_PBLTOP from 
!        "pbl_mix_mod.f" (bmy, 2/22/05)
!  (9 ) Modified for SO4s, NITs.  Also modified for alkalinity w/in the
!        seasalt chemistry. (bec, bmy, 4/13/05)
!  (10) Now remove reference to CMN, it's obsolete.  Now reference 
!        ITS_IN_THE_STRAT from "tropopause_mod.f" (bmy, 8/22/05)
!  (11) Now references XNUMOL from "tracer_mod.f" (bmy, 10/25/05)
!  (12) Updated to match JPL 2006 + full chem (jaf, bmy, 10/15/09)
!  (13) Now prevent floating-point exceptions when taking the exponential
!        terms. (win, bmy, 1/4/10)
!  (14) Save aqueous production rate to PSO4_SO2AQ for TOMAS microphyics
!        (win, 1/25/10)
!  (15) Added extra error checks to prevent negative L2S, L3S (bmy, 4/28/10)
!  (16) Use liq. water content from met fields in GEOS-5 (jaf, bmy, 6/30/10)
!  26 Aug 2010 - R. Yantosca - Use liquid water content from MERRA
!  12 Nov 2010 - R. Yantosca - Prevent div-by-zero when computing L2S and L3S
!  27 May 2011 - L. Zhang    - Divide LWC by cloud fraction for GEOS/MERRA
!                              and adjust the L2S and L3S rates accordingly
!  22 Dec 2011 - M. Payer    - Added ProTeX headers
!  08 Feb 2012 - R. Yantosca - Treat GEOS-5.7.2 in the same way as MERRA
!  01 Mar 2012 - R. Yantosca - Now use GET_AREA_CM2(I,J,L) from grid_mod.F90
!  31 Jul 2012 - R. Yantosca - Now loop over 1..LLPAR for GIGC compatibility
!  31 Jul 2012 - R. Yantosca - Declare temp drydep arrays w/ LLPAR (not LLTROP)
!  14 Nov 2012 - R. Yantosca - Add am_I_Root, Input_Opt, RC as arguments
!  15 Nov 2012 - M. Payer    - Replaced all met field arrays with State_Met
!                              derived type object
!  05 Mar 2013 - R. Yantosca - Now use Input_Opt%LNLPBL
!  19 Mar 2013 - R. Yantosca - Now copy Input_Opt%TCVV(1:N_TRACERS) and
!                              Input_Opt%XNUMOL(1:N_TRACERS) -- avoid OOB errs
!  25 Mar 2013 - M. Payer    - Now pass State_Chm object via the arg list
!  05 Sep 2013 - M. Sulprizio- Add modifications for cloud pH (B. Alexander)
!   6 Sep 2013 - M. Sulprizio- Bug fix: Prevent divide-by-zero if LWC=0. Only
!                              do aqueous SO2 chemistry when LWC>0.
!  26 Sep 2013 - R. Yantosca - Renamed GEOS_57 Cpp switch to GEOS_FP
!  28 Jan 2014 - R. Yantosca - Bug fix for TOMAS. Set ALKdst=0 since TOMAS
!                              carries its own dust tracers instead of DST1-4.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
!      REAL*8,  PARAMETER    :: HPLUS  = 3.16227766016837953d-5  !pH = 4.5
      REAL*8,  PARAMETER    :: HPLUS_45  = 3.16227766016837953d-5  !pH = 4.5
      REAL*8,  PARAMETER    :: HPLUS_50  = 1.0d-5  !pH = 5.0
      REAL*8,  PARAMETER    :: MINDAT = 1.d-20
!
! !LOCAL VARIABLES:
!
      ! Scalars
      LOGICAL               :: LNLPBL
      LOGICAL               :: IS_OFFLINE
      LOGICAL               :: IS_FULLCHEM
      INTEGER               :: I,      J,       L       !,      I1,   I2
      INTEGER               :: II,     NSTEP,   N_TRACERS
      INTEGER               :: BULK,   SIZE_RES 
      REAL*8                :: K0,     Ki,      KK,     M,    L1
      REAL*8                :: L2,     L3,      Ld,     F,    Fc
      REAL*8                :: RK,     RKT,     DTCHEM, DT_T, TK
      REAL*8                :: F1,     RK1,     RK2,    RK3,  SO20
      REAL*8                :: SO2_cd, H2O20,   O3,     L2S,  L3S
      REAL*8                :: LWC,    KaqH2O2, KaqO3,  PATM, FLUX
      REAL*8                :: ALK,    ALK1,    ALK2,    SO2_ss
      REAL*8                :: Kt1,    Kt2,     AREASS1, AREASS2
      REAL*8                :: PSO4E,  PSO4F,   Kt1N,    Kt2N
      REAL*8                :: XX,     AREA_CM2
      REAL*8                :: HPLUS,  SO4nss, TNH3,   TNO3,  GNO3, ANIT
      REAL*8                :: LSTOT,  ALKdst, ALKss,  ALKds, NH3
      REAL*8                :: SSCvv,  aSO4,   SO2_sr, SR,    TANIT

      ! Arrays
      REAL*8                :: ND44_TMP(IIPAR,JJPAR,LLPAR)
      REAL*8                :: TCVV  (Input_Opt%N_TRACERS)
      REAL*8                :: XNUMOL(Input_Opt%N_TRACERS)

      ! Pointers
      ! We need to define local arrays to hold corresponding values 
      ! from the Chemistry State (State_Chm) object. (mpayer, 12/6/12)
      REAL*8, POINTER :: STT(:,:,:,:)

      !=================================================================
      ! CHEM_SO2 begins here!
      !=================================================================

      IF ( IDTH2O2 == 0 .or. IDTSO2 == 0 .or. DRYSO2 == 0 ) RETURN

      ! Assume success
      RC          = GIGC_SUCCESS

      ! Copy fields from INPUT_OPT to local variables for use below
      IS_FULLCHEM = Input_Opt%ITS_A_FULLCHEM_SIM
      IS_OFFLINE  = Input_Opt%ITS_AN_AEROSOL_SIM
      N_TRACERS   = Input_Opt%N_TRACERS
      TCVV        = Input_Opt%TCVV  (1:N_TRACERS)
      XNUMOL      = Input_Opt%XNUMOL(1:N_TRACERS)
      LNLPBL      = Input_Opt%LNLPBL

      ! Initialize GEOS-Chem tracer array [kg] from Chemistry State object
      ! (mpayer, 12/6/12)
      STT => State_Chm%Tracers

      ! Read HNO3 for offline simulation
      IF ( IS_OFFLINE ) THEN
         IF ( ITS_A_NEW_MONTH() ) THEN
            CALL GET_GLOBAL_HNO3( GET_MONTH() )
         ENDIF
      ENDIF

      ! DTCHEM is the chemistry timestep in seconds
      DTCHEM   = GET_TS_CHEM() * 60d0

      ! Factor to convert AIRDEN from [kg air/m3] to [molec air/cm3]
      F        = 1000.d0 / AIRMW * 6.022d23 * 1.d-6

      ! Zero ND44_TMP array
      IF ( ND44 > 0 ) ND44_TMP = 0d0
      
      ! Loop over tropospheric grid boxes
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, SO20, H2O20, O3, PATM, TK, K0, M, KK, F1, RK1  )
!$OMP+PRIVATE( RK2, RK, RKT, SO2_cd, L1, Ld, L2, L2S, L3, L3S, FC, LWC )
!$OMP+PRIVATE( KaqH2O2, KaqO3, AREA_CM2, FLUX, ALK, ALK1, ALK2         )
!$OMP+PRIVATE( Kt1, Kt2, AREASS1, AREASS2, SO2_ss, Kt1N, Kt2N          )
!$OMP+PRIVATE( PSO4E, PSO4F, XX                                        )
!$OMP+PRIVATE( HPLUS, SO4nss, TNH3, TNO3,  GNO3, ANIT,   LSTOT, ALKdst )
!$OMP+PRIVATE( ALKds, ALKss,  NH3,  SSCvv, aSO4, SO2_sr, SR,    TANIT  )
!$OMP+PRIVATE( BULK,  SIZE_RES                                         )
!$OMP+SCHEDULE( DYNAMIC )
      DO L = 1, LLPAR 
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Initialize for safety's sake 
         AREA_CM2 = 0d0
         FLUX     = 0d0
         Ld       = 0d0

         ! Skip stratospheric boxes
         IF ( ITS_IN_THE_STRAT( I, J, L, State_Met ) ) CYCLE

         ! Initial SO2, H2O2 and O3 [v/v]
         SO20   = STT(I,J,L,IDTSO2)         
         H2O20  = STT(I,J,L,IDTH2O2)
         O3     = GET_O3( I, J, L, Input_Opt, State_Met )

         ! PATM  : Atmospheric pressure in atm
         PATM   = GET_PCENTER( I, J, L ) / 1013.25d0

         ! TK : Temperature [K]
         TK     = State_Met%T(I,J,L)

         IF ( IS_OFFLINE ) THEN

            ! Gas phase SO4 production is done here in offline run only 
            ! Updated to match JPL 2006 + full chem (jaf, 10/14/09)
            K0  = 3.3d-31 * ( 300.d0 / TK )**4.3d0
            Ki  = 1.6d-12
            M   = State_Met%AIRDEN(L,I,J) * F
            KK  = K0 * M / Ki
            F1  = ( 1.d0 + ( LOG10( KK ) )**2 )**( -1 )
            RK1 = ( K0 * M / ( 1.d0 + KK ) ) * 0.6d0**F1 * 
     &            GET_OH( I, J, L, Input_Opt, State_Met )

         ELSE 

            ! For online runs, SMVGEAR deals w/ this computation,
            ! so we can simply set RK1 = 0 (rjp, bmy, 3/23/03)
            K0  = 0.d0
            M   = 0.d0
            KK  = 0.d0
            F1  = 0.d0
            RK1 = 0.d0

         ENDIF

         ! SO2 drydep frequency [1/s].  Also accounts for the fraction
         ! of grid box (I,J,L) that is located beneath the PBL top.
         RK2    = DEPSAV(I,J,DRYSO2) * GET_FRAC_UNDER_PBLTOP( I, J, L )

         ! Add option for non-local PBL (Lin, 03/31/09)
         IF ( LNLPBL ) RK2 = 0.D0

         ! RK: total reaction rate [1/s]
         RK     = ( RK1 + RK2 )
       
         ! RKT: RK * DTCHEM [unitless] (bmy, 6/1/00)
         RKT    =  RK * DTCHEM

         !==============================================================
         ! Update SO2 conc. after gas phase chemistry and deposition
         !==============================================================
         IF ( RK > 0.d0 ) THEN
            SO2_cd = ( SO20  * EXP( -RKT ) ) +
     &               ( PSO2_DMS(I,J,L) * ( 1.d0 - EXP( -RKT ) ) / RKT )

            L1     = ( SO20 - SO2_cd + PSO2_DMS(I,J,L) ) * RK1/RK
             
            Ld     = ( SO20 - SO2_cd + PSO2_DMS(I,J,L) ) * RK2/RK
            
         ELSE
            SO2_cd = SO20
            L1     = 0.d0
         ENDIF

         !==============================================================
         ! Update SO2 conc. after seasalt chemistry (bec, 12/7/04)
         !==============================================================

         ! Get alkalinity of accum (ALK1) and coarse (ALK2) [kg]
         CALL GET_ALK( I, J, L, ALK1, ALK2, Kt1, Kt2, Kt1N, Kt2N,
     &                 State_Met )

         ! Total alkalinity [kg]
         ALK = ALK1 + ALK2

         ! If (1) there is alkalinity, (2) there is SO2 present, and 
         ! (3) O3 is in excess, then compute seasalt SO2 chemistry
         IF  ( ( ALK    > MINDAT )  .AND.
     &         ( SO2_cd > MINDAT )  .AND. 
     &         ( SO2_cd < O3     ) ) THEN

            ! Compute oxidation of SO2 -> SO4 and condensation of
            ! HNO3 -> nitrate within the seasalt aerosol
            CALL SEASALT_CHEM( I,         J,         L,        ALK1, 
     &                         ALK2,      SO2_cd,    Kt1,      Kt2,   
     &                         Kt1N,      Kt2N,      SO2_ss,   PSO4E,     
     &                         PSO4F,     am_I_Root, Input_Opt, 
     &                         State_Met, State_Chm, RC               )

         ELSE

            ! Otherwise set equal to zero
            SO2_ss       = SO2_cd
            PSO4E        = 0.d0
            PSO4F        = 0.d0
            PNITS(I,J,L) = 0.d0

         ENDIF

         !==============================================================
         ! Update SO2 concentration after cloud chemistry          
         ! SO2 chemical loss rate = SO4 production rate [v/v/timestep]
         !==============================================================
#if   defined ( GEOS_5 ) || defined( MERRA ) || defined( GEOS_FP )

         !---------------------------------------------
         ! GEOS-5/MERRA: Get LWC, FC from met fields
         ! (jaf, bmy, 6/30/10)
         !---------------------------------------------

         ! Get cloud fraction from met fields
         FC      = State_Met%CLDF(L,I,J)

         ! Get liquid water content [m3 H2O/m3 air] within cloud from met flds
         ! Units: [kg H2O/kg air] * [kg air/m3 air] * [m3 H2O/1e3 kg H2O]
         LWC     = State_Met%QL(I,J,L) * State_Met%AIRDEN(L,I,J) * 1D-3

         ! LWC is a grid-box averaged quantity. To improve the representation 
         ! of sulfate chemistry, we divide LWC by the cloud fraction and 
         ! compute sulfate chemistry based on the LWC within the cloud.  We 
         ! get the appropriate grid-box averaged mass of SO2 and sulfate by 
         ! multiplying these quantities by FC AFTER computing the aqueous 
         ! sulfur chemistry within the cloud. (lzh, jaf, bmy, 5/27/11)
         LWC     = SAFE_DIV( LWC, FC, 0d0 )

#else
         !---------------------------------------------
         ! Otherwise, compute FC, LWC as before
         !---------------------------------------------

         ! Volume cloud fraction (Sundqvist et al 1989) [unitless]
         FC      = VCLDF(I,J,L)

         ! Liquid water content in cloudy area of grid box [m3/m3]
         ! LWC as returned from the GET_LWC function is the in-cloud liquid
         ! water content.  To improve the representation of sulfate chemistry, 
         ! we use this LWC to compute the aqueous sulfate chemistry.  We then
         ! get the appropriate grid-box averaged mass of SO2 and sulfate by 
         ! multiplying by FC AFTER computing the aqueous sulfur chemistry
         ! within the cloud. (lzh, jaf, bmy, 5/27/11)
         LWC = GET_LWC( TK )

#endif

         ! Zero variables
         KaqH2O2 = 0.d0
         KaqO3   = 0.d0
         L2      = 0.d0
         L3      = 0.d0
         L2S     = 0.d0
         L3S     = 0.d0
         
         ! If (1) there is cloud, (2) there is SO2 present, and 
         ! (3) the T > -15 C, then compute aqueous SO2 chemistry
         ! Prevent divide-by-zero if LWC=0 (mpayer, 9/6/13)
         IF ( ( FC     > 0.d0   )  .AND. 
     &        ( SO2_ss > MINDAT )  .AND. 
     &        ( TK     > 258.0  )  .AND.
     &        ( LWC    > 0.d0   ) ) THEN

            !===========================================================
            ! NOTE...Sulfate production from aquatic reactions of SO2 
            ! with H2O2 & O3 is computed here and followings are 
            ! approximations or method used for analytical (integral) 
            ! solution of these computations. Please email us 
            ! (rjp@io.harvard.edu or bmy@io.harvard.edu) if you find
            ! anything wrong or questionable. 
            ! 
            ! 1) with H2O2(aq)
            !      [HSO3-] + [H+] + [H2O2(aq)] => [SO4=]     (rxn)
            !      d[SO4=]/dt = k[H+][HSO3-][H2O2(aq)] (M/s) (rate)
            !
            ! we can rewrite k[H+][HSO3-] as K1 pSO2 hSO2, 
            ! where pSO2 is equilibrium vapor pressure of SO2(g) 
            ! in atm, and hSO2 is henry's law constant for SO2
            !
            ! Therefore, rate can be written as 
            !
            !       k * K1 * pSO2 * hSO2 * pH2O2 * hH2O2,
            !
            ! where pH2O2 is equilibrium vapor pressure of H2O2(g), 
            ! and hH2O2 is henry's law constant for H2O2. Detailed 
            ! values are given in AQCHEM_SO2 routine.
            ! 
            ! Let us define a fraction of gas phase of A species 
            ! in equilibrium with aqueous phase as 
            !
            !        xA  = 1/(1+f), 
            !
            ! where  f   = hA * R * T * LWC, 
            !        hA  = Henry's constant,
            !        R   = gas constant, 
            !        T   = temperature in kelvin, 
            !        LWC = liquid water content [m3/m3]
            !
            ! As a result, the rate would become:
            !
            !    d[SO4=]   
            !    ------- = k K1 hSO2 hH2O2 xSO2 xH2O2 P P [SO2][H2O2]
            !      dt      
            !      ^       ^                            ^   ^    ^
            !      |       |____________________________|   |    |
            !
            !   mole/l/s               mole/l/s            v/v  v/v
            !
            !
            ! And we multiply rate by (LWC * R * T / P) in order to 
            ! convert unit from mole/l/s to v/v/s
            !
            ! Finally we come to 
            !
            !    d[SO4=]  
            !    ------- = KaqH2O2 [SO2][H2O2], 
            !      dt 
            !
            ! where
            !
            !   KaqH2O2 = k K1 hSO2 hH2O2 xSO2 xH2O2 P LWC R T, 
            !
            ! this new rate corresponds to a typical second order 
            ! reaction of which analytical (integral) solution is 
            !
            !   X  = A0 B0 ( exp[(A0-B0) Ka t] - 1 ) 
            !      / ( A0 exp[(A0-B0) Ka t] - B0 ) 
            !
            ! inserting variables into solution then we get
            ! [SO4=] =  [SO2][H2O2](exp[([SO2]-[H2O2]) KaqH2O2 t] - 1 )
            !        / ( [SO2] exp[([SO2]-[H2O2]) KaqH2O2 t] - [H2O2] )
            !
            ! Note...Exactly same method can be applied to O3 reaction 
            ! in aqueous phase with different rate constants. 
            !===========================================================

	    ! Get concentrations for cloud pH calculation (bec, 12/23/11)

	    ! Get sulfate concentration and convert from [v/v] to
            ! [moles/liter]
	    ! Use a cloud scavenging ratio of 0.7 
	    SO4nss  =  STT(I,J,L,IDTSO4) * State_Met%AIRDEN(L,I,J) *
     &                 0.7d0 / ( 28.97d0 * LWC )

	    ! Get total ammonia (NH3 + NH4+) concentration [v/v]
	    ! Use a cloud scavenging ratio of 0.7 for NH4+
            TNH3     = ( STT(I,J,L,IDTNH4) * 0.7d0 ) + 
     &                   STT(I,J,L,IDTNH3)

	    ! Get total nitrate (HNO3 + NIT) concentrations [v/v] 
     	    ! Use a cloud scavenging ratio of 0.7 for NIT
	    IF ( IS_FULLCHEM ) THEN
               TNO3 = STT(I,J,L,IDTHNO3) + 
     &              ( STT(I,J,L,IDTNIT) * 0.7d0 )
	       GNO3 = STT(I,J,L,IDTHNO3) !For Fahey & Pandis decision algorithm
            ELSE IF ( IS_OFFLINE ) THEN
               TANIT = STT(I,J,L,IDTNIT) !aerosol nitrate [v/v]
	       GNO3  = GET_HNO3_VV(I,J,L) - TANIT ! gas-phase nitric acid [v/v]
   	       ANIT  = TANIT * 0.7d0 ! aerosol nitrate in the cloud drops [v/v]
	       TNO3  = GNO3 + ANIT   ! total nitrate for cloud pH calculations
	    ENDIF
		
            ! Calculate cloud pH
	    CALL GET_HPLUS( SO4nss, TNH3, TNO3,     SO2_ss, TK, 
     &                        PATM,   LWC,  HPLUS_45, HPLUS )

            ! Compute aqueous rxn rates for SO2
            CALL AQCHEM_SO2( LWC, TK,    PATM,    SO2_ss, H2O20, 
     &                       O3,  HPLUS, KaqH2O2, KaqO3 ) 

            !----------------------------------------------------------
            ! Compute loss by H2O2.  Prevent floating-point exception
            ! by not allowing the exponential to go to infinity if 
            ! the argument is too large.  (win, bmy, 1/4/09)
            !----------------------------------------------------------

            ! Argument of the exponential
            XX  = ( SO2_ss - H2O20 ) * KaqH2O2 * DTCHEM

            ! Test if EXP(XX) can be computed w/o numerical exception
            IF ( IS_SAFE_EXP( XX ) .and. ABS( XX ) > 0d0 ) THEN

               ! Aqueous phase SO2 loss rate w/ H2O2 [v/v/timestep]
               L2  = EXP( XX )

               ! Loss by H2O2
               L2S = SO2_ss * H2O20 * ( L2 - 1.D0 ) / 
     &               ( (SO2_ss * L2) - H2O20 )  
            ELSE

               ! NOTE from Jintai Lin (4/28/10):
               ! However, in the case of a negative XX, L2S should be 
               ! approximated as SO2_ss, instead of H2O20. In other words, 
               ! L2S = SO2_ss * H2O20 * ( L2 - 1.D0 ) / ( (SO2_ss*L2) - H2O20 )
               ! reaches different limits when XX reaches positive infinity 
               ! and negative infinity.
               IF ( XX > 0.d0 ) THEN 
                  L2S = H2O20 
               ELSE 
                  L2S = SO2_ss
               ENDIF

            ENDIF

            !----------------------------------------------------------
            ! Compute loss by O3.  Prevent floating-point exception
            ! by not allowing the exponential to go to infinity if 
            ! the argument is too large. (win, bmy, 1/4/09)
            !----------------------------------------------------------

            ! Argument of the exponential
            XX = ( SO2_ss - O3 ) * KaqO3 * DTCHEM 

            ! Test if EXP(XX) can be computed w/o numerical exception
            IF ( IS_SAFE_EXP( XX ) .and. ABS( XX ) > 0d0 ) THEN

               ! Aqueous phase SO2 loss rate w/ O3 [v/v/timestep]
               L3  = EXP( XX )

               ! Loss by O3
               L3S = SO2_ss * O3 * (L3 - 1.D0) / ((SO2_ss * L3) - O3)  

            ELSE
 
               ! Follow the same logic for L3S as described in
               ! Jintai Lin's note above (bmy, 4/28/10)
               IF ( XX > 0.d0 ) THEN 
                  L3S = O3 
               ELSE 
                  L3S = SO2_ss 
               ENDIF
            ENDIF

            ! Decide whether or not it is necessary to use heterogeneous cloud
            ! pH calculations based on the Fahey and Pandis, 2001 decision
            ! algorithm (bec, 12/23/11)

            ! Add up total seasalt and dust and convert to ug/m3
            ! Note that it is better to use dust and sea-salt alkalinity
            ! tracers if these are being transported (bec, 12/23/11)

#if defined( TOMAS )
            !%%%%%%%%%%%%%%%%% BUG FIX FOR TOMAS %%%%%%%%%%%%%%%%%%%%%%%
            ! NOTE: TOMAS uses its own dust tracers and does not
            ! carry ALKdst.  Set ALKdst to zero here. (bmy, 1/28/14)
            ALKdst = 0d0
#else
            ! For other simulations, Sum up the contributions from
            ! DST1 thru DST4 tracers into ALKdst. (bmy, 1/28/14)
            ALKdst = ( STT(I,J,L,IDTDST1) + STT(I,J,L,IDTDST2) +
     &                 STT(I,J,L,IDTDST3) + STT(I,J,L,IDTDST4) ) *
     &                 1.D9 * State_Met%AD(I,J,L) / TCVV(IDTDST1) /
     &                 State_Met%AIRVOL(I,J,L)
#endif

            ALKss  = ( STT(I,J,L,IDTSALA  ) + STT(I,J,L,IDTSALC) ) *
     &                 1.D9 * State_Met%AD(I,J,L) / TCVV(IDTSALA) /
     &                 State_Met%AIRVOL(I,J,L)

            ALKds = ALKdst + ALKss

            ! Get NH3 concentrations (v/v)
            NH3 = STT(I,J,L,IDTNH3) 

            ! Initialize
            BULK = 0
	    SIZE_RES = 0

            ! Fahey and Seinfeld decision algorithm
	    IF ( H2O20 > SO2_ss + 1D-9 ) THEN
               BULK = 1
            ELSEIF( LWC < 0.1D6 ) THEN !10^6 coversion from m3/m3 --> g/m3
               SIZE_RES = 1
            ELSEIF( gno3 > NH3 ) THEN
               IF ( SO2_ss >= 5.D-9          .and. 
     &              H2O20  >= SO20   ) 
     &              BULK    = 1
               IF ( LWC    >= 0.3D6          .and. 
     &              SO2_ss >= 3.D-9          .and. 
     &              H2O20  >= SO2_ss ) 
     &              BULK    = 1
               IF ( ALKds  >= 5.D0           .and. 
     &              LWC    >= 0.5D6          .and. 
     &              H2O20  >= SO2_ss ) 
     &              BULK    = 1
               IF ( LWC    >= 0.1D6          .and. 
     &              gno3   <= (NH3 + 2.D-9) ) 
     &              BULK    = 1
            ELSEIF( LWC    >= 0.1d6 ) THEN
               IF ( NH3    <= 1.d-9          .and. 
     &              ALKds  >= 5.d0   ) 
     &              BULK    = 1
               IF ( NH3    <= 5.d-9          .and. 
     &              ALKds  >= 10.d0  ) 
     &              BULK    = 1
               IF ( gno3   <= 1.d-9          .and. 
     &              NH3    >= (gno3 + 2.d-9) .and. 
     &              SO2_ss <= 7.d-9 )
     &              BULK    = 1
               IF ( gno3   <= 1.d-9          .and. 
     &              NH3    >= (gno3 + 2.d-9) .and. 
     &              ALKds  >= 2.d0 )  BULK = 1
               IF ( gno3   <= 3.d-9          .and. 
     &              NH3    >= (gno3 + 4.d-9) )
     &              BULK    = 1
               IF ( gno3   <= 7.d-9          .and. 
     &              NH3    >= (gno3 + 3.d-9) .and. 
     &              SO2_ss <= 5.d-9 )
     &              BULK    = 1
               IF ( gno3   <= 7.d-9          .and. 
     &              NH3    >= (gno3 + 3.d-9) .and. 
     &		    ALKds  >= 4.d0           .and. 
     &              SO2_ss <= 9.d-9  ) 
     &              BULK    = 1
               IF ( ALKds  >= 3.d0           .and. 
     &              NH3    <= 3.d-9          .and. 
     &              SO2_ss <= 4.d-9 ) 
     &              BULK    = 1
               IF ( ALKds  >= 5.d0           .and. 
     &              SO2_ss <= 5.d-9          .and. 
     &              NH3    <= 7.d-9 )  
     &              BULK    = 1
               IF ( NH3    >= (gno3 + 2.d-9) .and. 
     &              SO2_ss <= 5.d-9 )  
     &              BULK    = 1
               IF ( NH3    >= (gno3 + 4.d-9) .and. 
     &              SO2_ss <= 10.d-9 )
     &              BULK    = 1
               IF ( ALKds  >= 2.d0           .and. 
     &              NH3    <= 10.d-9         .and. 
     &              H2O20  >= SO2_ss ) 
     &              BULK    = 1
               IF ( NH3    <= 1.d-9          .and. 
     &              SO2_ss >= 3.d-9          .and. 
     &              H2O20  >= SO2_ss ) 
     &              BULK    = 1
            ELSEIF ( LWC   >= 0.3d6 ) THEN
               IF ( NH3    >= (gno3 + 5.d-9) .and. 
     &              SO2_ss <= 10.d-9 )
     &              BULK    = 1
               IF ( gno3   <= 1.d-9          .and. 
     &              NH3    >= (gno3 + 2.d-9) )
     &              BULK    = 1
               IF ( gno3   <= 7.d-9          .and. 
     &              NH3    >= (gno3 + 3.d-9) )
     &              BULK    = 1
               IF ( ALKds  >= 3.d0           .and. 
     &              NH3 <= 10d-9             .and. 
     &              SO2_ss <= 5d-9 ) 
     &              BULK    = 1
               IF ( ALKds  >= 5.d0           .and. 
     &              NH3    <= 10.d-9         .and. 
     &              SO2_ss <= 5.d-9 ) 
     &              BULK    = 1
               IF ( SO2_ss >= 1.5d-9         .and. 
     &              H2O20  >= SO2_ss ) 
     &              BULK    = 1
               IF ( NH3    <= 12.d-9         .and. 
     &              ALKds  >=10.d0 ) 
     &              BULK    = 1
               IF ( NH3    <= 1.d-9          .and. 
     &              ALKds  >= 4.d0           .and. 
     &              SO2_ss <= 10.d-9 ) 
     &              BULK    = 1
               IF ( NH3    <= 5.d-9          .and. 
     &              ALKds  >= 6.d0           .and. 
     &              SO2_ss <= 10.d-9 ) 
     &              BULK    = 1
               IF ( NH3    <= 7.d-9          .and. 
     &              ALKds   >-8.d0           .and. 
     &              SO2_ss <= 10.d-9 ) 
     &              BULK    = 1
            ELSEIF ( LWC   >= 0.5d6 ) THEN
               IF ( H2O20  >= (0.9d0 * SO2_ss) ) 
     &              BULK    = 1
               IF ( NH3    <= 1.d-9          .and. 
     &              ALKds  >= 5.d0           .and. 
     &              SO2_ss <= 10.d-9 ) 
     &              BULK    = 1
            ELSE
               SIZE_RES = 1
            ENDIF

            ! Decide whether or not to perform sulfate production rate
            ! enhancement due to cloud drop heterogenity in pH over the oceans
            ! (bec, 12/23/11)
            IF ( SIZE_RES == 1 .AND. IS_WATER( I, J, State_MEt) .AND. 
     &           TK > 268.15 ) THEN

	       ! Get total in-cloud sulfate production based on bulk cloud pH
               ! calculations for use in HET_DROP_CHEM
               LSTOT = L2S + L3S

	       ! Get coarse-mode sea-salt concentration for use in
               ! HET_DROP_CHEM [v/v]
               ! Note that it is better to use coarse sea salt alkalinity
               ! tracer if it is being transported (bec, 12/23/11)
               SSCvv  = STT(I,J,L,IDTSALC)

	       ! Get sulfate concentrations for use in HET_DROP_CHEM [v/v]
               aSO4  =  STT(I,J,L,IDTSO4) 

               ! This is to make sure HET_DROP_CHEM does not compute more
               ! sulfate then there is SO2
               SO2_sr = MAX( SO2_ss - ( ( L2S + L3S ) * FC ), MINDAT )

               CALL HET_DROP_CHEM( I,    J,   L,      LSTOT, SSCvv, 
     &                             aSO4, NH3, SO2_sr, H2O20, GNO3,  SR,
     &                             State_Met )

	    ELSE
		SR = 0.d0
	    ENDIF

            ! We have used the in-cloud LWC to compute the sulfate
            ! aqueous chemistry.  We get the appropriate grid-box averaged 
            ! mass of SO2 and sulfate by multiplying the reaction rates
            ! L2S and L3s by the cloud fraction after the aqueous chemistry
            ! has been done.  (lzh, jaf, bmy, 5/27/11)
            L2S =  L2S * FC
            L3S =  L3S * FC    
	    SR  =  SR  * FC

            ! Make sure SO2_ss and H2O20 are in the proper range
!            SO2_ss = MAX( SO2_ss - ( L2S + L3S ), MINDAT )
            SO2_ss = MAX( SO2_ss - ( L2S + L3S + SR ), MINDAT )
            H2O20  = MAX( H2O20  - L2S,                MINDAT )

            ! Update SO2 level, save SO2[ppv], H2O2[ppv] for WETDEP
            SO2s( I,J,L) = SO2_ss
            H2O2s(I,J,L) = H2O20

         ELSE

            ! Otherwise, don't do aqueous chemistry, and
            ! save the original concentrations into SO2 and H2O2
            H2O2s(I,J,L) = MAX( H2O20,  1.0d-32 )
            SO2s(I,J,L ) = MAX( SO2_ss, 1.0d-32 )
            L2S          = 0.d0
            L3S          = 0.d0
	    SR		 = 0.d0

         ENDIF

         ! Store updated SO2, H2O2 back to the tracer arrays 
         STT(I,J,L,IDTSO2)  = SO2s( I,J,L)
         STT(I,J,L,IDTH2O2) = H2O2s(I,J,L)

         ! SO2 chemical loss rate  = SO4 production rate [v/v/timestep]
!         PSO4_SO2(I,J,L) = L1 + L2S + L3S + PSO4E
         PSO4_SO2(I,J,L) = L1 + L2S + L3S + PSO4E + SR
         PSO4_ss (I,J,L) = PSO4F

#if   defined( TOMAS )
         PSO4_SO2AQ(I,J,L) = L2S + L3S  ! For TOMAS microphysics
#endif

         !=================================================================
         ! ND05 Diagnostics [kg S/timestep]
         !=================================================================
         IF ( ND05 > 0 .and. L <= LD05 ) THEN
           
            ! P(SO4) from gas-phase oxidation [kg S/timestep]
            AD05(I,J,L,5) = AD05(I,J,L,5) +
     &                      ( L1  * State_Met%AD(I,J,L) / TCVV_S )

            ! P(SO4) from aqueous-phase oxidation with H2O2 [kg S/timestep]
            AD05(I,J,L,6) = AD05(I,J,L,6) +
     &                      ( L2S * State_Met%AD(I,J,L) / TCVV_S )

            ! P(SO4) from aqueous-phase oxidation with O3 [kg S/timestep]
!            AD05(I,J,L,7) = AD05(I,J,L,7) +
!     &                      ( L3S * State_Met%AD(I,J,L) / TCVV_S )
            AD05(I,J,L,7) = AD05(I,J,L,7) +
     &                  ( ( L3S + SR ) * State_Met%AD(I,J,L) / TCVV_S )

            ! P(SO4) from O3 oxidation in sea-salt aerosols [kg S/timestep]
            AD05(I,J,L,8) = AD05(I,J,L,8) +
     &                      ( (PSO4E + PSO4F) * State_Met%AD(I,J,L) /
     &                      TCVV_S )

         ENDIF

         !=================================================================
         ! ND44 Diagnostic: Drydep flux of SO2 [molec/cm2/s]
         !=================================================================
         IF ( ND44 > 0 .AND. Ld > 0d0 ) THEN

            ! Surface area [cm2]
            AREA_CM2 = GET_AREA_CM2( I, J, L )

            ! Convert [v/v/timestep] to [molec/cm2/s]
            FLUX = Ld   * State_Met%AD(I,J,L)      / TCVV(IDTSO2)
            FLUX = FLUX * XNUMOL(IDTSO2) / AREA_CM2 / DTCHEM
            
            ! Store dryd flx in ND44_TMP as a placeholder
            ND44_TMP(I,J,L) = ND44_TMP(I,J,L) + FLUX
         ENDIF
      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      !===============================================================
      ! ND44: Sum drydep fluxes by level into the AD44 array in
      ! order to ensure that  we get the same results w/ sp or mp 
      !===============================================================
      IF ( ND44 > 0 ) THEN 
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
         DO J = 1, JJPAR
         DO I = 1, IIPAR
         DO L = 1, LLPAR
            AD44(I,J,DRYSO2,1) = AD44(I,J,DRYSO2,1) + ND44_TMP(I,J,L)
         ENDDO
         ENDDO
         ENDDO
!$OMP END PARALLEL DO
      ENDIF

      ! Free pointer
      NULLIFY( STT )

      END SUBROUTINE CHEM_SO2
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: seasalt_chem
!
! !DESCRIPTION: Subroutine SEASALT\_CHEM computes SO4 formed from S(IV) + O3 on
!  seasalt aerosols as a function of seasalt alkalinity. (bec, bmy, 4/13/05,
!  10/7/08)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE SEASALT_CHEM ( I,         J,         L,
     &                          ALK1,      ALK2,      SO2_cd,
     &                          Kt1,       Kt2,       Kt1N,
     &                          Kt2N,      SO2_ss,    PSO4E,
     &                          PSO4F,     am_I_Root, Input_Opt,
     &                          State_Met, State_Chm, RC )
!
! !USES:
!
      !---------------------------------------------------------------
      ! DIAGNOSTICS -- leave commented out for now (bec, bmy, 4/13/05)
      !USE CMN_DIAG_MOD          ! ND19
      !USE DIAG_MOD,        ONLY : AD09
      !---------------------------------------------------------------
      USE CMN_GCTM_MOD          
      USE CMN_SIZE_MOD          
      USE COMODE_MOD,         ONLY : CSPEC
      USE COMODE_MOD,         ONLY : JLOP
      USE COMODE_MOD,         ONLY : VOLUME
      USE ERROR_MOD,          ONLY : GEOS_CHEM_STOP
      USE ERROR_MOD,          ONLY : IT_IS_NAN
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE GIGC_State_Chm_Mod, ONLY : ChmState
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE GLOBAL_HNO3_MOD,    ONLY : GET_HNO3_UGM3
      USE TIME_MOD,           ONLY : GET_TS_CHEM
      USE TIME_MOD,           ONLY : GET_ELAPSED_SEC
      USE TIME_MOD,           ONLY : GET_MONTH 
      USE TIME_MOD,           ONLY : ITS_A_NEW_MONTH
      USE TRACER_MOD,         ONLY : XNUMOLAIR
      USE TRACERID_MOD,       ONLY : IDTHNO3
      USE TRACERID_MOD,       ONLY : IDTNITs
      USE TRACERID_MOD,       ONLY : IDTSO2
      USE TRACERID_MOD,       ONLY : IDTSO4
      USE TRACERID_MOD,       ONLY : IDTSO4s
!
! !INPUT PARAMETERS: 
!

      INTEGER,        INTENT(IN)   :: I, J, L    ! Grid box indices
      REAL*8,         INTENT(IN)   :: SO2_cd     ! SO2 mixing ratio [v/v] after 
                                                 !  gas phase chemistry and 
                                                 !  dry deposition
      REAL*8,         INTENT(IN)   :: Kt1, Kt2   ! Rate constant [s-1] for 
                                                 !  sulfate formation on sea
                                                 !  salt aerosols from GET_ALK 
                                                 !  (1=fine; 2=coarse)
      REAL*8,         INTENT(IN)   :: Kt1N, Kt2N 
      REAL*8,         INTENT(IN)   :: ALK1, ALK2 ! Alkalinity [kg] from 
                                                 !  seasalt_mod
      TYPE(MetState), INTENT(IN)   :: State_Met  ! Meteorology State object
      TYPE(OptInput), INTENT(IN)   :: Input_Opt  ! Input Options object
      LOGICAL,        INTENT(IN)   :: am_I_Root  ! Are we on the root CPU?
!
! !INPUT/OUTPUT PARAMETERS: 
!
      TYPE(ChmState), INTENT(INOUT) :: State_Chm ! Chemistry State object
!
! !OUTPUT PARAMETERS: 
!
      REAL*8,         INTENT(OUT)   :: SO2_ss    ! SO2 mixing ratio [v/v] 
                                                 !  after sea salt chemistry
      REAL*8,         INTENT(OUT)   :: PSO4E     ! SO4E (sulfate produced by 
                                                 !  S(IV)+O3 on fine seasalt) 
                                                 !  mixing ratio [v/v]
      REAL*8,         INTENT(OUT)   :: PSO4F     ! SO4F (sulfate produced by 
                                                 !  S(IV)+O3 on coarse seasalt)
      INTEGER,        INTENT(OUT)   :: RC        ! Success or failure?
!
! !REMARKS:
!  Chemical reactions:
!  ============================================================================
!  (R1) SO2 + O3 + ALK => SO4 + O2
!       Modeled after Chamedies and Stelson, 1992?
!
! !REVISION HISTORY:
!  (1 ) Now references XNUMOLAIR from "tracer_mod.f" (bmy, 10/25/05)
!  (2 ) Bug fix: now avoid seg fault error if IDTHNO3 is zero, as it would
!        be for an offline aerosol simulation. (bmy, 3/29/06)
!  (3 ) Fixed typo in FALK_A_SO2 equation: C_FLUX_C should be C_FLUX_A.
!        (havala, bec, bmy, 12/8/06)
!  (4 ) Bug fix for mass balance, replace TITR_HNO3 w/ HNO3_SSC in the
!        expression for HNO3_ss.  Bug fix: now do equivalent computation 
!        for GET_GNO3, which is now no longer called because it's in 
!        "isoropia_mod.f". (bec, bmy, 7/30/08)
!  22 Dec 2011 - M. Payer    - Added ProTeX headers
!  09 Nov 2012 - M. Payer    - Replaced all met field arrays with State_Met
!                              derived type object
!  25 Mar 2013 - M. Payer    - Now pass State_Chm object via the arg list
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
      REAL*8,  PARAMETER :: MINDAT    = 1.0d-20 
      REAL*8,  PARAMETER :: TCVV_HNO3 = 28.97d0 / 63.0d0 
!
! !LOCAL VARIABLES:
!
      ! Scalars
      LOGICAL            :: IT_IS_A_FULLCHEM_SIM
      INTEGER            :: JLOOP,       N_TRACERS
      REAL*8             :: SO2_chem,    DTCHEM
      REAL*8             :: O3_cspec,    O3_lost
      REAL*8             :: EQ_1_C,      EQ_2_C
      REAL*8             :: SO4E,        SO2_new,    SO4F
      REAL*8             :: SO2_eq,      N_FLUX_A,   N_FLUX_C
      REAL*8             :: END_ALK,     L5A,        L5C
      REAL*8             :: EQ1,         EQ2,        TITR_SO2
      REAL*8             :: TITR_HNO3,   NIT_vv,     NITs_vv
      REAL*8             :: NIT0,        NITS0
      REAL*8             :: F_SO2,       FALK_A_SO2, FALK_C_SO2
      REAL*8             :: EQ_BEG,      F_SO2_A,    F_SO2_C
      REAL*8             :: ALKA,        ALKC,       TOTAL_ACID_FLUX
      REAL*8             :: HNO3_EQ,     TOT_FLUX_A, TOT_FLUX_C
      REAL*8             :: FALK_A_HNO3, HNO3_vv
      REAL*8             :: FALK_C_HNO3, F_HNO3_A,   F_HNO3_C
      REAL*8             :: EQ_1_N,      EQ_2_N,     F_HNO3
      REAL*8             :: HNO3_SSA,    HNO3_SSC,   N_FLUX
      REAL*8             :: HNO3_EQ_C,   L6A,        L6C   
      REAL*8             :: C_FLUX_A,    C_FLUX_C,   C_FLUX
      REAL*8             :: HNO3_ss,     HNO3_kg

      ! Arrays
      REAL*8             :: TCVV(Input_Opt%N_TRACERS)

      ! Pointers
      REAL*8, POINTER    :: STT(:,:,:,:)
      REAL*8, POINTER    :: AD(:,:,:)
      REAL*8, POINTER    :: AIRVOL(:,:,:)

      !=================================================================
      ! SEASALT_CHEM begins here!
      !=================================================================

      ! Assume success
      RC                   = GIGC_SUCCESS

      ! Copy fields from INPUT_OPT to local variables for use below
      N_TRACERS            = Input_Opt%N_TRACERS
      TCVV                 = Input_Opt%TCVV(1:N_TRACERS)
      IT_IS_A_FULLCHEM_SIM = Input_Opt%ITS_A_FULLCHEM_SIM

      ! Initialize pointers
      STT                 => State_Chm%Tracers
      AD                  => State_Met%AD
      AIRVOL              => State_Met%AIRVOL

      ! DTCHEM is the chemistry timestep in seconds
      DTCHEM = GET_TS_CHEM() * 60d0

      ! Convert SO2 [v/v] to  [eq/gridbox]
      SO2_eq = ( 2.d0 * SO2_cd  * AD(I,J,L) ) / 
     &         ( TCVV(IDTSO2) * 0.064d0 )
      SO2_eq = MAX( SO2_eq, MINDAT )

      IF ( IT_IS_A_FULLCHEM_SIM ) THEN

 	 ! Convert HNO3 [v/v] to [equivalents]
         HNO3_vv = STT(I,J,L,IDTHNO3)
         HNO3_eq = HNO3_vv * AD(I,J,L) / ( 28.97d0 / 63d0 ) / 63.d-3

      ELSE
         
         ! Get HNO3 in ug/m3, then multiply by volume in m3
         ! and 1e-6 kg/ug to get HNO3 in kg
         HNO3_kg = GET_HNO3_UGM3( I, J, L, State_Met ) *
     &             AIRVOL(I,J,L) * 1e-6

	 ! Convert HNO3 [kg] first to [v/v] 
         HNO3_vv = HNO3_kg * ( 28.97d0 / 63d0 ) / AD(I,J,L)
 
         ! Then convert HNO3 [kg] to [equivalents]
         HNO3_eq = HNO3_kg / 63d-3

      ENDIF

      !-----------
      ! SO2
      !-----------

      ! Available flux of SO2 to accum sea salt aerosols [v/v/timestep]
      L5A      = EXP( -Kt1 * DTCHEM )
      F_SO2_A  = SO2_cd * ( 1.d0 - L5A )
      F_SO2_A  = MAX( F_SO2_A, 1.d-32 )

      ! Convert to [eq/timestep] 
      C_FLUX_A = 2.d0 * F_SO2_A * AD(I,J,L) / TCVV(IDTSO2) / 0.064d0

      ! Available flux of SO2 to coarse sea salt aerosols [v/v/timestep]
      L5C      = EXP( - Kt2 * DTCHEM )
      F_SO2_C  = SO2_cd * ( 1.d0 - L5C )
      F_SO2_C  = MAX( F_SO2_C, 1.0d-32 )

      ! Convert to [eq/timestep] 
      C_FLUX_C = 2.d0 * F_SO2_C * AD(I,J,L) / TCVV(IDTSO2) / 0.064d0

      ! Total flux of SO2 [v/v/timestep]
      F_SO2    = F_SO2_A + F_SO2_C 

      ! Total flux of SO2 [eq/timestep]
      C_FLUX   = C_FLUX_A + C_FLUX_C 

      !-----------
      ! HNO3
      !-----------

      ! Available flux of HNO3 to accum sea salt aerosols [v/v/timestep]
      L6A = EXP( - Kt1N * DTCHEM )
      F_HNO3_A = HNO3_vv * ( 1.D0 - L6A )
      F_HNO3_A = MAX( F_HNO3_A, 1.0D-32 )

      ! Convert to [eq/timestep] 
      N_FLUX_A = F_HNO3_A * AD(I,J,L) / ( 28.97d0 / 63d0 ) / 0.063d0

      ! Available flux of HNO3 to coarse sea salt aerosols [v/v/timestep]
      L6C = EXP( - Kt2N * DTCHEM )
      F_HNO3_C = HNO3_vv * ( 1.D0 - L6C )
      F_HNO3_C = MAX( F_HNO3_C, 1.0D-32 )

      ! convert to [eq/timestep] 
      N_FLUX_C = F_HNO3_C * AD(I,J,L) / ( 28.97d0 / 63d0 ) / 0.063d0

      ! Total flux of HNO3
      F_HNO3 = F_HNO3_A + F_HNO3_C ![v/v/timestep]
      N_FLUX = N_FLUX_A + N_FLUX_C ![eq/timestep]

      !-----------
      ! Acid
      !-----------

      ! Total acid flux to accum sea-salt aerosols [eq/box/timestep]
      TOT_FLUX_A = C_FLUX_A + N_FLUX_A 
      TOT_FLUX_A = MAX( TOT_FLUX_A, MINDAT )

      ! Total acid flux to coarse sea-salt aerosols [eq/box/timestep]
      TOT_FLUX_C = C_FLUX_C + N_FLUX_C 
      TOT_FLUX_C = MAX( TOT_FLUX_C, MINDAT )

      ! Total  acid flux to sea salt aerosols
      TOTAL_ACID_FLUX = TOT_FLUX_A + TOT_FLUX_C

      ! Total available alkalinity [eq]
      EQ1 = ALK1 * 0.07d0
      EQ2 = ALK2 * 0.07d0

      !----------------------------------------------------------------------
      ! NOTE: This was a sensitivity simulation, keep for future reference
      !       cf Alexander et al 2005 (bec, bmy, 4/13/05)
      !! Total available alkalinity [eq] doubled for Sievering run
      !EQ1 = ALK1 * 0.14d0
      !EQ2 = ALK2 * 0.14d0
      !----------------------------------------------------------------------

      !----------------------------------------------------------------------
      ! DIAGNOSTIC -- leave uncommented for now (bec, bmy, 4/13/05)
      !! Write out beginning alkalinity [eq SO4]
      !EQ_BEG = EQ1 + EQ2
      !IF ( ND09 > 0 ) AD09(I,J,L,1) = AD09(I,J,L,1) + EQ_BEG
      !----------------------------------------------------------------------

      IF ( TOT_FLUX_A > EQ1 ) THEN

	 ! Fraction of alkalinity available for each acid
         FALK_A_SO2  = C_FLUX_A / TOT_FLUX_A
	 FALK_A_HNO3 = N_FLUX_A / TOT_FLUX_A
         FALK_A_SO2  = MAX( FALK_A_SO2, MINDAT )
	 FALK_A_HNO3 = MAX( FALK_A_HNO3, MINDAT )

      ELSE

	 FALK_A_SO2  = 1.0d0
	 FALK_A_HNO3 = 1.0d0

      ENDIF
      
      IF ( TOT_FLUX_C > EQ2 ) THEN

         ! Fraction of flkalinity available for each acid
	 FALK_C_SO2  = C_FLUX_C/TOT_FLUX_C
	 FALK_C_HNO3 = N_FLUX_C/TOT_FLUX_C
         FALK_C_SO2  = MAX( FALK_C_SO2, MINDAT )
	 FALK_C_HNO3 = MAX( FALK_C_HNO3, MINDAT )

      ELSE

	 FALK_C_SO2  = 1.0d0
	 FALK_C_HNO3 = 1.0d0

      ENDIF

      ! Alkalinity available for S(IV) --> S(VI)
      EQ_1_C       = EQ1 * FALK_A_SO2
      EQ_1_C       = MAX( EQ_1_C, MINDAT )
      EQ_1_N       = EQ1 * FALK_A_HNO3
      EQ_1_N       = MAX( EQ_1_N, MINDAT )
                  
      EQ_2_C       = EQ2 * FALK_C_SO2
      EQ_2_C       = MAX( EQ_2_C, MINDAT )
      EQ_2_N       = EQ2 * FALK_C_HNO3
      EQ_2_N       = MAX( EQ_2_N, MINDAT )

      !-----------------
      ! Fine Seasalt
      !-----------------

      ! don't produce more SO4 than available ALK or SO2
      SO4E         = MIN( C_FLUX_A, EQ_1_C, SO2_eq ) 
      SO4E         = MAX( SO4E, MINDAT )

      ! Update SO2 concentration [eq/box] 
      SO2_new      = SO2_eq - SO4E
      SO2_new      = MAX( SO2_new, MINDAT )

      !-----------------
      ! Coarse Seasalt
      !-----------------     
      IF ( SO2_new > MINDAT ) THEN

 	 ! don't produce more SO4 than available ALK or SO2
	 SO4F      = MIN( C_FLUX_C, SO2_new, EQ_2_C ) 
	 SO4F      = MAX( SO4F, MINDAT )

	 !Update SO2 concentration [eq] 
	 SO2_chem  = SO2_new - SO4F
	 SO2_chem  = MAX( SO2_chem, MINDAT )
      ELSE
	 SO4F      = MINDAT
	 SO2_chem  = MINDAT
      ENDIF

      ! Alkalinity titrated by S(IV) --> S(VI) [eq]
      TITR_SO2     = SO4E + SO4F

      !-------------------------------------------------------------------
      ! DIAGNOSTIC -- leave uncommented for now
      !! write out in diagnostic
      !IF ( ND09 > 0 ) AD09(I,J,L,2) = AD09(I,J,L,2) + TITR_SO2
      !-------------------------------------------------------------------

      !Modified SO2 [eq] converted back to [v/v]
      SO2_ss       = SO2_chem * 0.064 * TCVV(IDTSO2) /
     &               AD(I,J,L) / 2.0d0
      SO2_ss       = MAX( SO2_ss, MINDAT )

      !SO4E produced converted from [eq/timestep] to [v/v/timestep] 
      PSO4E        = SO4E * 0.096 * TCVV(IDTSO4) /
     &               AD(I,J,L) / 2.0d0

      !SO4F produced converted from [eq/timestep] to [v/v/timestep] 
      PSO4F        = SO4F * 0.096 * TCVV(IDTSO4S) /
     &               AD(I,J,L) / 2.0d0

      ! Alkalinity titrated by HNO3
      HNO3_SSA     = MIN(N_FLUX_A, HNO3_EQ, EQ_1_N)
      HNO3_SSA     = MAX(HNO3_SSA, MINDAT)
      HNO3_EQ_C    = HNO3_EQ - HNO3_SSA
      HNO3_EQ_C    = MAX(HNO3_EQ_C, MINDAT)
      HNO3_SSC     = MIN(N_FLUX_C, HNO3_EQ_C, EQ_2_N)
      HNO3_SSC     = MAX(HNO3_SSC, MINDAT)
      TITR_HNO3    = HNO3_SSA + HNO3_SSC

      !----------------------------------------------------------------------
      ! DIAGNOSTIC -- leave commented out for now
      ! !write out alkalinity titrated by HNO3(g)
      !IF ( ND09 > 0 ) AD09(I,J,L,3) = AD09(I,J,L,3) + TITR_HNO3
      !----------------------------------------------------------------------

      ! HNO3 lost [eq/timestep] converted back to [v/v/timestep]
      IF ( IDTHNO3 > 0 ) THEN

         ! Coupled sim: IDTHNO3 is defined, so use it
         HNO3_ss = HNO3_SSC * 0.063 * TCVV(IDTHNO3) / AD(I,J,L)
         STT(I,J,L,IDTHNO3) = MAX( HNO3_vv - HNO3_ss, MINDAT )

      ELSE

         ! Offline aerosol sim: IDTHNO3 isn't defined, use TCVV_HNO3
         HNO3_ss = TITR_HNO3 * 0.063 * TCVV_HNO3 /AD(I,J,L)

      ENDIF

      ! NITS produced converted from [eq/timestep] to [v/v/timestep] 
      PNITs(I,J,L) = HNO3_SSC * 0.063 * TCVV(IDTNITS) / AD(I,J,L)
	 
      ! Modified accum alkalinity 
      ALKA         = EQ1 - (SO4E + HNO3_SSA)
      ALKA         = MAX( ALKA, MINDAT )

      !------------------------------------------------------------------------
      ! Uncomment this if you want to transport alkalinity (bec, bmy, 4/13/05)
      ![eq] --> [kg] --> [v/v] use this only if transporting alkalinity
      !ALKAvv = (ALKA * TCVV(IDTSAL1))/(7.0d-2 * AD(I,J,L))
      !ALKAvv = MAX( ALKAvv, MINDAT )
      !------------------------------------------------------------------------

      ! Modified accum alkalinity 
      ALKC         = EQ2 - (SO4F + HNO3_SSC)
      ALKC         = MAX( ALKC, MINDAT )
      
      !------------------------------------------------------------------------
      ! Uncomment this if you want to transport alkalinity (bec, bmy, 4/13/05)
      !! [eq] --> [kg] --> [v/v] use this only if transporting alkalinity
      !ALKCvv = (ALKC * TCVV(IDTSAL2))/(7.0d-2 * AD(I,J,L))
      !ALKCvv = MAX(ALKCvv, MINDAT)
      !------------------------------------------------------------------------

      !------------------------------------------------------------------------
      ! DIAGNOSTIC -- leave commented out for now (bec, bmy, 4/13/05)
      !! write out ending alkalinity
      !END_ALK = ALKA + ALKC
      !IF ( ND09 > 0 ) AD09(I,J,L,4) = AD09(I,J,L,4) + END_ALK
      !------------------------------------------------------------------------

      ! Free pointer
      NULLIFY( STT    )
      NULLIFY( AD     )
      NULLIFY( AIRVOL )

      END SUBROUTINE SEASALT_CHEM
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_hplus
!
! !DESCRIPTION: Subroutine GET\_HPLUS computes H+ concentrations in cloud
!  liquid water for pH dependent cloud chemistry. (bec, 4/11/11)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE GET_HPLUS( SO4nss, TNH3, TNO3, SO2,
     &                      T,      PRES, LWC,  iHPLUS, HPLUS )
!
! !USES:
!
      USE ERROR_MOD,       ONLY : IT_IS_NAN, GEOS_CHEM_STOP

!
! !INPUT PARAMETERS: 
!
      REAL*8,  INTENT(IN)    :: SO4nss ! Total nss sulfate mixing ratio [M]
      REAL*8,  INTENT(IN)    :: TNO3   ! Total nitrate (gas+particulate) mixing
                                       ! ratio [v/v]
      REAL*8,  INTENT(IN)    :: TNH3   ! NH3 mixing ratio [v/v]
      REAL*8,  INTENT(IN)    :: SO2    ! SO2 mixing ratio [v/v] 
      REAL*8,  INTENT(IN)    :: T      ! Temperature [K]
      REAL*8,  INTENT(IN)    :: PRES   ! Pressure [atm]
      REAL*8,  INTENT(IN)    :: LWC    ! Cloud liquid water content [m3/m3]
      REAL*8,  INTENT(IN)    :: iHPLUS ! Initial [H+] [M]
!
! !OUTPUT PARAMETERS: 
!
      REAL*8,  INTENT(OUT)   :: HPLUS  ! Calculated [H+] [M]
! !REMARKS:
!  Calculation:
!  ============================================================================
!  Solve the following electroneutrality equation:
!  [H+] = 2[SO4]nss + [Cl] + [OH] + [HCO3] + 2[CO3] + [HSO3] + 2[SO3] + [NO3]
!        - [Na] - 2[Ca] - [K] - 2[Mg] - [NH4]
!
!  Aqueous concentrations of [Cl], [Na], [Ca], [K], and [Mg] come from 
!  ISORROPIA II
!
!  Let concentrations of [HCO3], [CO3], [HSO3], [SO3], [NO3] and [NH4] evolve
!  according to Henry's law equilibrium.
!
!  Assume [S(VI)] = [SO4]nss (this applies for pH > 3)
!            
! !REVISION HISTORY:
!  25 Jan 2012 - M. Payer    - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!

      ! Water dissociation constants
      REAL*8,  PARAMETER     :: Kw   = 1.0d-14
      REAL*8,  PARAMETER     :: DhrKw = -6710.d0
      REAL*8,  PARAMETER     :: MINVAL = 0.01
!
! !LOCAL VARIABLES:
!
      REAL*8                 :: D, Kw_T, ipH, newpH, nHPLUS
      REAL*8                 :: kCO2p, kCO2p2
      REAL*8                 :: kSO2p, kSO2p2
      REAL*8                 :: kHNO3p, kNH3p, kHClp
      INTEGER                :: count
      REAL*8                 :: E, F, G, H, P, Q, R, A, B
      REAL*8                 :: A2, A1, A0, X, Y
      REAL*8                 :: CRUTES( 3 ) ! Coeff and roots of cubic equation
      INTEGER                :: NR          ! Number of roots to cubic equation

      !=================================================================
      ! GET_HPLUS begins here!
      !=================================================================

      ! Initial pH guess 
      ipH = -log10(iHPLUS)

      ! Non-volatile aerosol concentration [M] 
      ! For now sulfate is the only non-volatile species
      D = (2.d0*SO4nss) 

      ! Temperature dependent water equilibrium constant
      Kw_T = Kw*exp(DhrKw*((1.d0/T)-(1.d0/298.d0)))

      ! Initialize
      newpH   = 0.0 
      COUNT = 0

      DO WHILE ( ABS(ipH-newpH) .gt. MINVAL )

      COUNT = COUNT+1
 
         IF ( COUNT .EQ. 1 ) THEN
            ipH = ipH
	 ELSE
            ipH = newpH
	 ENDIF

         nHPLUS = 10.d0**(-ipH)

         kCO2p  = kCO21 ( PRES, T, LWC, nHPLUS ) 

	 kCO2p2 = kCO22 ( PRES, T, LWC, nHPLUS ) 

	 kSO2p  = kSO21 ( PRES, T, LWC, nHPLUS, SO2 )

	 kSO2p2 = kSO22 ( PRES, T, LWC, nHPLUS, SO2 )

	 kHNO3p = kHNO3 ( PRES, T, LWC, nHPLUS, TNO3 )

	 kNH3p  = kNH3  ( PRES, T, LWC, nHPLUS, TNH3, Kw_T )

         ! Save this for later in case we want to include HCl in cloud pH
         ! calculations
         !kHClp  = kHCl  ( PRES, T, LWC, nHPLUS, CL  ) 

         !E = KHClp + Kw_T + KCO2p + KSO2p + KHNO3p !This would be used if HCl is included
	 E = Kw_T + KCO2p + KSO2p + KHNO3p 

	 F = KCO2p2 + KSO2p2

	 G = KNH3p

	 H = 1 + G

	 P = -(D/H)

	 Q = -E/H

	 R = F/H

	 A = (1.d0/3.d0)*((3.d0*Q)-(P*P))

  	 B = (1.d0/27.d0)*((2.0d0*P*P*P)-(9.0d0*P*Q)+(27.0d0*R))

 	 A2 = 0.
 	 A1 = A
 	 A0 = B

	 !write(6,*) 'calling CUBIC'
         CALL CUBIC ( A2, A1, A0, NR, CRUTES )
	 !write(6,*) 'after CUBIC'

         ! Code assumes the smallest positive root is in CRUTES(1)
         X = CRUTES( 1 )

	 ! Y = [H+]
 	 Y = X - (P/3.0d0)

  	 ! Set minimum [H+] = 1.d-14 (pH = 14)
    	 Y = MAX(Y,1.0d-14)

  	 ! Set maximum [H+] = 1.d-1 (pH = 1)
    	 Y = MIN(Y,1.0d-1)

  	 ! If solution does not converge after 5 iterations
  	 ! average last 2 pH calculations
  	 IF (count > 5) THEN
	    newpH = ((-log10(Y)) + (-log10(nHPLUS))) / 2.0d0

            IF (IT_IS_NAN( newpH )) THEN
               write(6,*) 'newpH = ', newpH
               write(6,*) 'Y = ', Y
               write(6,*) 'nHPLUS = ', nHPLUS
               CALL GEOS_CHEM_STOP
            ENDIF

            EXIT
         ELSE 
            newpH = -log10(Y)

            IF (IT_IS_NAN( newpH )) THEN
               write(6,*) 'newpH = ', newpH
               write(6,*) 'Y = ', Y
               CALL GEOS_CHEM_STOP
            ENDIF

         ENDIF

      ENDDO

      HPLUS = 10.0d0**(-newpH)

      END SUBROUTINE GET_HPLUS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: kCO21
!
! !DESCRIPTION: Function kCO21
!\\
!\\
! !INTERFACE:
!
      FUNCTION kCO21 ( P, T, LWC, HPLUS ) RESULT ( KCO2p )
!
! !INPUT PARAMETERS: 
!
      REAL*8,  INTENT(IN) :: T, P, LWC, HPLUS
!
! !OUTPUT PARAMETERS: 
!
      REAL*8              :: KCO2p, KCO2p2
!
! !REVISION HISTORY:
!  25 Jan 2012 - M. Payer    - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
      ! CO2 dissociation constants
      REAL*8,  PARAMETER  :: Kc1 = 4.3e-7
      REAL*8,  PARAMETER  :: Kc2 =4.68e-11
      REAL*8,  PARAMETER  :: DhrKc1 = -1000.
      REAL*8,  PARAMETER  :: DhrKc2 = -1760.
      REAL*8,  PARAMETER  :: Hco2 = 3.4e-2
      REAL*8,  PARAMETER  :: Dhco2 = -4.85/0.04
      ! CO2 concentration [v/v]
      REAL*8,  PARAMETER  :: CO2 = 380.0d-6 !380 ppmv
!
! !LOCAL VARIABLES:
!
      REAL*8              :: Hco2_T, Kc1_T, Kc2_T
      REAL*8              :: Hco2eff, xCO2, pCO2

      !=================================================================
      ! kCO2 begins here!
      !=================================================================

      !CO2 dissolution constants
      Hco2_T = Hco2*exp(Dhco2*((1.d0/T)-(1.d0/298.d0)))
      Kc1_T = Kc1*exp(DhrKc1*((1.d0/T)-(1.d0/298.d0)))
      Kc2_T = Kc2*exp(DhrKc1*((1.d0/T)-(1.d0/298.d0)))

      !CO2 dissolution  
      Hco2eff = Hco2_T*(1.d0+(Kc1_T/HPLUS)+((Kc1_T*Kc2_T)/
     &         (HPLUS*HPLUS)))
      xCO2 = 1.d0 / ( 1.d0 + ( Hco2eff * 0.08205d0 * T * LWC ) )
      pCO2 = CO2 * P * xCO2

      KCO2p  = Hco2_T * Kc1_T * pCO2

      END FUNCTION kCO21
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: kCO22
!
! !DESCRIPTION: Function kCO22
!\\
!\\
! !INTERFACE:
!
      FUNCTION kCO22 ( P, T, LWC, HPLUS ) RESULT ( KCO2p2 )
!
! !INPUT PARAMETERS: 
!
      REAL*8,  INTENT(IN) :: T, P, LWC, HPLUS
!
! !OUTPUT PARAMETERS: 
!
      REAL*8              :: KCO2p, KCO2p2
!
! !REVISION HISTORY:
!  25 Jan 2012 - M. Payer    - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
      ! CO2 dissociation constants
      REAL*8,  PARAMETER  :: Kc1 = 4.3e-7
      REAL*8,  PARAMETER  :: Kc2 =4.68e-11
      REAL*8,  PARAMETER  :: DhrKc1 = -1000.
      REAL*8,  PARAMETER  :: DhrKc2 = -1760.
      REAL*8,  PARAMETER  :: Hco2 = 3.4e-2
      REAL*8,  PARAMETER  :: Dhco2 = -4.85/0.04
      ! CO2 concentration [v/v]
      REAL*8,  PARAMETER  :: CO2 = 380.0d-6 
!
! !LOCAL VARIABLES:
!
      REAL*8              :: Hco2_T, Kc1_T, Kc2_T
      REAL*8              :: Hco2eff, xCO2, pCO2

      !=================================================================
      ! kCO2 begins here!
      !=================================================================

      !CO2 dissolution constants
      Hco2_T = Hco2*exp(Dhco2*((1.d0/T)-(1.d0/298.d0)))
      Kc1_T = Kc1*exp(DhrKc1*((1.d0/T)-(1.d0/298.d0)))
      Kc2_T = Kc2*exp(DhrKc1*((1.d0/T)-(1.d0/298.d0)))

      !CO2 dissolution  
      Hco2eff = Hco2_T*(1.d0+(Kc1_T/HPLUS)+((Kc1_T*Kc2_T)/
     &         (HPLUS*HPLUS)))
      xCO2 = 1.d0 / ( 1.d0 + ( Hco2eff * 0.08205d0 * T * LWC ) )
      pCO2 = CO2 * P * xCO2

      KCO2p  = Hco2_T * Kc1_T * pCO2
      KCO2p2 = 2.d0 * KCO2p *  Kc2_T

      END FUNCTION kCO22
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: kSO21
!
! !DESCRIPTION: Function kSO21
!\\
!\\
! !INTERFACE:
!
      FUNCTION kSO21 ( P, T, LWC, HPLUS, SO2 ) RESULT ( KSO2p )
!
! !INPUT PARAMETERS: 
!
      REAL*8,  INTENT(IN) :: T, P, LWC, HPLUS, SO2
!
! !OUTPUT PARAMETERS: 
!
      REAL*8              :: KSO2p, KSO2p2
!
! !REVISION HISTORY:
!  25 Jan 2012 - M. Payer    - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
      ! SO2 dissociation constants
      REAL*8,  PARAMETER  :: Ks1 = 1.3e-2
      REAL*8,  PARAMETER  :: Ks2 = 6.6e-8
      REAL*8,  PARAMETER  :: Hso2 = 1.23
      REAL*8,  PARAMETER  :: Dhso2 = -6.25/0.04
      REAL*8,  PARAMETER  :: DhrKso21 = 1960.
      REAL*8,  PARAMETER  :: DhrKso22 = 1500.
!
! !LOCAL VARIABLES:
!
      REAL*8              :: Hso2_T, Ks1_T, Ks2_T
      REAL*8              :: Hso2eff, xSO2, pSO2

      !=================================================================
      ! kSO2 begins here!
      !=================================================================


      ! SO2 dissolution constants
      Hso2_T = Hso2*exp(Dhso2*((1.d0/T)-(1.d0/298.d0)))
      Ks1_T = Ks1*exp(DhrKso21*((1.d0/T)-(1.d0/298.d0)))
      Ks2_T = Ks2*exp(DhrKso22*((1.d0/T)-(1.d0/298.d0)))

      ! SO2 dissolution
      Hso2eff = Hso2_T*(1.d0+(Ks1_T/HPLUS)+((Ks1_T*Ks2_T)/
     &         (HPLUS*HPLUS)))
      xSO2 = 1.d0 / ( 1.d0 + ( Hso2eff * 0.08205d0 * T * LWC ) )
      pSO2 = SO2 * P * xSO2

      KSO2p  = Hso2_T * Ks1_T * pSO2

      END FUNCTION kSO21
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: kSO22
!
! !DESCRIPTION: Function kSO22
!\\
!\\
! !INTERFACE:
!
      FUNCTION kSO22 ( P, T, LWC, HPLUS, SO2 ) RESULT ( KSO2p2 )
!
! !INPUT PARAMETERS: 
!
      REAL*8,  INTENT(IN) :: T, P, LWC, HPLUS, SO2
!
! !OUTPUT PARAMETERS: 
!
      REAL*8              :: KSO2p, KSO2p2
!
! !REVISION HISTORY:
!  25 Jan 2012 - M. Payer    - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
      ! SO2 dissociation constants
      REAL*8,  PARAMETER  :: Ks1 = 1.3e-2
      REAL*8,  PARAMETER  :: Ks2 = 6.6e-8
      REAL*8,  PARAMETER  :: Hso2 = 1.23
      REAL*8,  PARAMETER  :: Dhso2 = -6.25/0.04
      REAL*8,  PARAMETER  :: DhrKso21 = 1960.
      REAL*8,  PARAMETER  :: DhrKso22 = 1500.
!
! !LOCAL VARIABLES:
!
      REAL*8              :: Hso2_T, Ks1_T, Ks2_T
      REAL*8              :: Hso2eff, xSO2, pSO2

      !=================================================================
      ! kSO2 begins here!
      !=================================================================


      ! SO2 dissolution constants
      Hso2_T = Hso2*exp(Dhso2*((1.d0/T)-(1.d0/298.d0)))
      Ks1_T  = Ks1 *exp(DhrKso21*((1.d0/T)-(1.d0/298.d0)))
      Ks2_T  = Ks2 *exp(DhrKso22*((1.d0/T)-(1.d0/298.d0)))

      !SO2 dissolution
      Hso2eff = Hso2_T*(1.d0+(Ks1_T/HPLUS)+((Ks1_T*Ks2_T)/
     &         (HPLUS*HPLUS)))
      xSO2 = 1.d0 / ( 1.d0 + ( Hso2eff * 0.08205d0 * T * LWC ) )
      pSO2 = SO2 * P * xSO2

      KSO2p  = Hso2_T * Ks1_T * pSO2
      KSO2p2 = 2.d0 * KSO2p * Ks2_T

      END FUNCTION kSO22
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: kHNO3
!
! !DESCRIPTION: Function kNO3
!\\
!\\
! !INTERFACE:
!
      FUNCTION kHNO3 ( P, T, LWC, HPLUS, HNO3 ) RESULT ( KHNO3p )
!
! !INPUT PARAMETERS: 
!
      REAL*8,  INTENT(IN) :: T, P, LWC, HPLUS, HNO3
!
! !OUTPUT PARAMETERS: 
!
      REAL*8              :: KHNO3p
!
! !REVISION HISTORY:
!  25 Jan 2012 - M. Payer    - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
      ! HNO3 dissociation constants
      REAL*8,  PARAMETER  :: Kn1 = 15.4
      REAL*8,  PARAMETER  :: Hhno3 = 2.1e5
      REAL*8,  PARAMETER  :: Dhhno3 = 0.
      REAL*8,  PARAMETER  :: DhrKn1 = 8700.
!
! !LOCAL VARIABLES:
!
      REAL*8              :: Hhno3_T, Kn1_T
      REAL*8              :: Hhno3eff, xHNO3, pHNO3

      !=================================================================
      ! kHNO3 begins here!
      !=================================================================

      ! HNO3 dissolution constants
      Hhno3_T = Hhno3*exp(Dhhno3*((1.d0/T)-(1.d0/298.d0)))
      Kn1_T = Kn1*exp(DhrKn1*((1.d0/T)-(1.d0/298.d0)))

      !HNO3 dissolution
      Hhno3eff = 3.2e6/HPLUS
      xHNO3 = 1.d0 / ( 1.d0 + ( Hhno3eff * 0.08205d0 * T * LWC ) )
      pHNO3 = HNO3 * P * xHNO3

      kHNO3p = Hhno3_T * Kn1_T * pHNO3

      END FUNCTION kHNO3
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: kHCl
!
! !DESCRIPTION: Function kHCl
!\\
!\\
! !INTERFACE:
!
      FUNCTION kHCl ( P, T, LWC, HPLUS, Cl ) RESULT ( KHClp )
!
! !INPUT PARAMETERS: 
!
      REAL*8,  INTENT(IN) :: T, P, LWC, HPLUS, Cl
!
! !OUTPUT PARAMETERS: 
!
      REAL*8              :: KHClp
!
! !REVISION HISTORY:
!  25 Jan 2012 - M. Payer    - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
      ! HNO3 dissociation constants
      REAL*8,  PARAMETER  :: Kcl = 1.74d6
      REAL*8,  PARAMETER  :: Hcl = 1.1d0
      REAL*8,  PARAMETER  :: Dhcl = -2.3158d0
      REAL*8,  PARAMETER  :: DhrKcl = 6900.d0
!
! !LOCAL VARIABLES:
!
      REAL*8              :: Hcl_T, Kcl_T
      REAL*8              :: Hcleff, xCl, pHCl

      !=================================================================
      ! kHCl begins here!
      !=================================================================

      ! HNO3 dissolution constants
      HCl_T = Hcl*exp(Dhcl*((1.0d0/T)-(1.0d0/298.0d0)))
      Kcl_T = Kcl*exp(DhrKcl*((1.0d0/T)-(1.0d0/298.0d0)))

      !HCl dissolution
      Hcleff = Hcl_T*(1.0d0+(Kcl_T/HPLUS))
      xCl = 1.0d0 / ( 1.0d0 + ( Hcleff * 0.08205d0 * T * LWC ) )
      pHCl = Cl * P * xCl

      kHClp = Hcl_T * Kcl_T * pHCl

      END FUNCTION kHCl
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: kNH3
!
! !DESCRIPTION: Function kNH3
!\\
!\\
! !INTERFACE:
!
      FUNCTION kNH3 ( P, T, LWC, HPLUS, NH3, Kw ) RESULT ( KNH3p )
!
! !INPUT PARAMETERS: 
!
      REAL*8,  INTENT(IN) :: T, P, LWC, HPLUS, NH3, Kw
!
! !OUTPUT PARAMETERS: 
!
      REAL*8              :: KNH3p
!
! !REVISION HISTORY:
!  25 Jan 2012 - M. Payer    - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
      ! NH3 dissociation contants
      REAL*8,  PARAMETER  :: Ka1 = 1.7e-5
      REAL*8,  PARAMETER  :: Hnh3 = 62.
      REAL*8,  PARAMETER  :: Dhnh3 = -8.17/0.04
      REAL*8,  PARAMETER  :: DhrKa1 = -450.

      ! Variables
      REAL*8              :: Hnh3_T, Ka1_T
      REAL*8              :: Hnh3eff, xNH3, pNH3

      !=================================================================
      ! kNH3 begins here!
      !=================================================================

      !NH3 dissolution constants
      Hnh3_T = Hnh3*exp(Dhnh3*((1.d0/T)-(1.d0/298.d0)))
      Ka1_T = Ka1*exp(DhrKa1*((1.d0/T)-(1.d0/298.d0)))

      !NH3 dissolution
      Hnh3eff = Hnh3_T*(1.d0+((Ka1_T* HPLUS) / Kw))
      xNH3 = 1.d0 / ( 1.d0 + ( Hnh3eff * 0.08205d0 * T * LWC ) )
      pNH3 = NH3 * P * xNH3

      KNH3p = Hnh3_T * Ka1_T * pNH3 / Kw

      END FUNCTION kNH3
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cubic
!
! !DESCRIPTION: Subroutine CUBIC finds the roots of a cubic equation / 3rd 
! order polynomial
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CUBIC( A2, A1, A0, NR, CRUTES )
!
! !USES:
!
!      USE ERROR_MOD, ONLY : GEOS_CHEM_STOP !ERROR_STOP
!
! !INPUT PARAMETERS: 
!
      INTEGER           :: NR
      REAL*8            :: A2, A1, A0
      REAL*8            :: CRUTES(3)
!
! !REMARKS:
! Formulae can be found in numer. recip.  on page 145
!   kiran  developed  this version on 25/4/1990
!   Dr. Francis S. Binkowski modified the routine on 6/24/91, 8/7/97
! ***
! *** modified 2/23/98 by fsb to incorporate Dr. Ingmar Ackermann's
!     recommendations for setting a0, a1,a2 as real*8 variables.
!
! Modified by Bob Yantosca (10/15/02) 
! - Now use upper case / white space
! - force double precision with "D" exponents
! - updated comments / cosmetic changes
! - now call ERROR_STOP from "error_mod.f" to stop the run safely
!
! !REVISION HISTORY:
!  25 Jan 2012 - M. Payer    - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
      REAL*8, PARAMETER :: ONE    = 1.0d0
      REAL*8, PARAMETER :: SQRT3  = 1.732050808d0
      REAL*8, PARAMETER :: ONE3RD = 0.333333333d0
!
! !LOCAL VARIABLES:
!
      REAL*8            :: QQ,    RR,    A2SQ,  THETA, DUM1, DUM2
      REAL*8            :: PART1, PART2, PART3, RRSQ,  PHI,  YY1
      REAL*8            :: YY2,   YY3,   COSTH, SINTH

      !=================================================================
      ! CUBIC begins here!
      !=================================================================
      A2SQ = A2 * A2
      QQ   = ( A2SQ - 3.d0*A1 ) / 9.d0
      RR   = ( A2*( 2.d0*A2SQ - 9.d0*A1 ) + 27.d0*A0 ) / 54.d0

      ! CASE 1 THREE REAL ROOTS or  CASE 2 ONLY ONE REAL ROOT
      DUM1 = QQ * QQ * QQ
      RRSQ = RR * RR
      DUM2 = DUM1 - RRSQ

      IF ( DUM2 .GE. 0.d0 ) THEN

         ! Now we have three real roots
         PHI = SQRT( DUM1 )

         IF ( ABS( PHI ) .LT. 1.d-20 ) THEN
            CRUTES(1) = 0.0d0
            CRUTES(2) = 0.0d0
            CRUTES(3) = 0.0d0
            NR        = 0
         ENDIF
         
         THETA = ACOS( RR / PHI ) / 3.0d0
         COSTH = COS( THETA )
         SINTH = SIN( THETA )

         ! Use trig identities to simplify the expressions
         ! Binkowski's modification
         PART1     = SQRT( QQ )
         YY1       = PART1 * COSTH
         YY2       = YY1 - A2/3.0d0
         YY3       = SQRT3 * PART1 * SINTH
         CRUTES(3) = -2.0d0*YY1 - A2/3.0d0
         CRUTES(2) = YY2 + YY3
         CRUTES(1) = YY2 - YY3

         ! Set negative roots to a large positive value
         IF ( CRUTES(1) .LT. 0.0d0 ) CRUTES(1) = 1.0d9
         IF ( CRUTES(2) .LT. 0.0d0 ) CRUTES(2) = 1.0d9
         IF ( CRUTES(3) .LT. 0.0d0 ) CRUTES(3) = 1.0d9

         ! Put smallest positive root in crutes(1)
         CRUTES(1) = MIN( CRUTES(1), CRUTES(2), CRUTES(3) )
         NR        = 3

      ELSE  

         ! Now here we have only one real root
         PART1     = SQRT( RRSQ - DUM1 )
         PART2     = ABS( RR )
         PART3     = ( PART1 + PART2 )**ONE3RD
         CRUTES(1) = -SIGN(ONE,RR) * ( PART3 + (QQ/PART3) ) - A2/3.D0
         CRUTES(2) = 0.D0
         CRUTES(3) = 0.D0
         NR        = 1

      ENDIF
      
      ! Return to calling program
      END SUBROUTINE CUBIC
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: aqchem_so2
!
! !DESCRIPTION: Subroutine AQCHEM\_SO2 computes the reaction rates for aqueous
! SO2 chemistry. (rjp, bmy, 10/31/02, 12/12/02)  
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE AQCHEM_SO2( LWC, T,     P,       SO2, H2O2, 
     &                       O3,  Hplus, KaqH2O2, KaqO3 ) 
!
! !INPUT PARAMETERS: 
!
      REAL*8, INTENT(IN)  :: LWC     ! Liq water content [m3/m3]=1.E-6*L [g/m3]
      REAL*8, INTENT(IN)  :: T       ! Temperature [K]
      REAL*8, INTENT(IN)  :: P       ! Pressure [atm]
      REAL*8, INTENT(IN)  :: SO2     ! SO2  mixing ratio [v/v]
      REAL*8, INTENT(IN)  :: H2O2    ! H2O2 mixing ratio [v/v]
      REAL*8, INTENT(IN)  :: O3      ! O3   mixing ratio [v/v]
      REAL*8, INTENT(IN)  :: HPLUS   ! Concentration of H+ ion (i.e. pH) [v/v]
!
! !OUTPUT PARAMETERS: 
!
      REAL*8, INTENT(OUT) :: KaqH2O2 ! Reaction rate for H2O2
      REAL*8, INTENT(OUT) :: KaqO3   ! Reaction rate for O3
!
! !REMARKS:
!  Chemical Reactions:
!  ============================================================================
!  (R1) HSO3- + H2O2(aq) + H+ => SO4-- + 2H+ + H2O [Jacob, 1986]   
!                                                                             .
!      d[S(VI)]/dt = k[H+][H2O2(aq)][HSO3-]/(1 + K[H+]) 
!      [Seinfeld and Pandis, 1998, page 366]
!                                                                             .
!  (R2) SO2(aq) + O3(aq) =>                                        
!       HSO3-   + O3(aq) =>  
!       SO3--   + O3(aq) =>
!       [Jacob, 1986; Jacobson, 1999]
!                                                                             .
!       d[S(VI)]/dt = (k0[SO2(aq)] + k1[HSO3-] + K2[SO3--])[O3(aq)]
!       [Seinfeld and Pandis, 1998, page 363]
!                                                                             .
!  Reaction rates can be given as
!       Ra     = k [H2O2(ag)] [S(IV)]  [mole/liter*s]  OR
!       Krate  = Ra LWC R T / P        [1/s]
!                                                                             .
!  Where:
!       LWC = Liquid water content(g/m3)*10-6 [m3(water)/m3(gas)]
!       R   = 0.08205  (atm L / mol-K), Universal gas const.
!       T   = Temperature (K)
!       P   = Pressure (atm)
!                                                                             .
!  Procedure:
!  ============================================================================
!  (a ) Given [SO2] which is assumed to be total SO2 (gas+liquid) in 
!        equilibrium between gas and liquid phase. 
!                                                                             .
!  (b ) We can compute SO2(g) using Henry's law
!          P(so2(g)) = Xg * [SO2]
!          Xg = 1/(1 + Faq), Fraction of SO2 in gas
!       where: 
!          Faq   = Kheff * R * T * LWC, 
!          KHeff = Effective Henry's constant
!                                                                             .
!  (c ) Then Calculate Aquous phase, S[IV] concentrations
!        S[IV] = Kheff * P(so2(g) in atm) [M]
!                                                                             .
!  (d ) The exact same procedure is applied to calculate H2O2(aq)
!
! !REVISION HISTORY:
!  (1 ) Updated by Rokjin Park (rjp, bmy, 12/12/02)
!  22 Dec 2011 - M. Payer    - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
      REAL*8, PARAMETER   :: R = 0.08205d0 
!
! !LOCAL VARIABLES:
!
      REAL*8              :: KH2O2,  RA,     KS1, KS2,    HCSO2
      REAL*8              :: FHCSO2, XSO2G,  SIV, HSO3,   XSO2AQ
      REAL*8              :: XHSO3,  XSO3,   KH1, HCH2O2, FHCH2O2
      REAL*8              :: XH2O2G, H2O2aq, KO0, KO1,    KO2
      REAL*8              :: HCO3,   XO3g,   O3aq

      !=================================================================
      ! AQCHEM_SO2 begins here!
      !
      ! Aqueous reaction rate
      ! HSO3- + H2O2 + H+ => SO4-- + 2H+ + H2O [Jacob, 1986]
      !=================================================================

      ! [Jacob, 1986]
      KH2O2 = 6.31d14 * EXP( -4.76d3 / T )  

!      ! [Jacobson, 1999]
!      KH2O2 = 7.45d07 * EXP( -15.96d0 * ( (298.15/T) - 1.) ) / 
!     &                  ( 1.d0 + 13.d0 * Hplus)

      !=================================================================
      ! Equilibrium reaction of SO2-H2O
      !    SO2 + H2O = SO2(aq)        (s0)
      !    SO2(ag)   = HSO3- + H+     (s1)
      !    HSO3-     = SO3-- + H+     (s2)
      !
      ! Reaction constant for Aqueous chemistry -- No big difference 
      ! between Jacob and Jacobson, choose one of them.
      !
      ! Reaction rate dependent on Temperature is given
      !   H = A exp ( B (T./T - 1) ) 
      !
      ! For equilibrium reactions of SO2:
      !            As1      Bs1   As2      Bs2  
      !  Seinfeld  1.30d-2  7.02  6.60d-8  3.76   [1998]
      !  Jacob     1.30d-2  6.75  6.31d-8  5.05   [1986]
      !  Jacobson  1.71d-2  7.04  5.99d-8  3.74   [1996]
      !=================================================================
      Ks1    = 1.30d-2 * EXP( 6.75d0 * ( 298.15d0 / T - 1.d0 ) )
      Ks2    = 6.31d-8 * EXP( 5.05d0 * ( 298.15d0 / T - 1.d0 ) )

      ! SIV Fraction
      XSO2aq = 1.d0/(1.d0 + Ks1/Hplus + Ks1*Ks2/(Hplus*Hplus))
      XHSO3  = 1.d0/(1.d0 + Hplus/Ks1 + Ks2/Hplus)
      XSO3   = 1.d0/(1.d0 + Hplus/Ks2 + Hplus*Hplus/(Ks1*Ks2))

      ! Henry's constant [mol/l-atm] and Effective Henry's constant for SO2
      HCSO2  = 1.22d0 * EXP( 10.55d0 * ( 298.15d0 / T - 1.d0) )         
      FHCSO2 = HCSO2 * (1.d0 + (Ks1/Hplus) + (Ks1*Ks2 / (Hplus*Hplus)))
      
      XSO2g  = 1.d0 / ( 1.d0 + ( FHCSO2 * R * T * LWC ) )
      SIV    = FHCSO2 * XSO2g * SO2 * P
!      HSO3   = Ks1 * HCSO2 * XSO2g * SO2 * P

      !=================================================================
      ! H2O2 equilibrium reaction
      ! H2O2 + H2O = H2O2.H2O
      ! H2O2.H2O   = HO2- + H+   1)
      !
      ! Reaction rate dependent on Temperature is given
      !   H = A exp ( B (T./T - 1) ) 
      !
      ! For equilibrium reactions of SO2
      !            Ah1       Bh1
      !  Jacob     1.58E-12  -12.49  [1986]
      !  Jacobson  2.20E-12  -12.52  [1996]
      !=================================================================
      Kh1 = 2.20d-12 * EXP( -12.52d0 * ( 298.15d0 / T - 1.d0 ) )

      ! Henry's constant [mol/l-atm] and Effective Henry's constant for H2O2
      ! [Seinfeld and Pandis, 1998]
      ! HCH2O2  = 7.45D4 * EXP( 24.48d0 * ( 298.15d0 / T - 1.d0) ) 

      ! [Jacobson,1999]
      HCH2O2  = 7.45D4 * EXP( 22.21d0 * (298.15d0 / T - 1.d0) )
      FHCH2O2 = HCH2O2 * (1.d0 + (Kh1 / Hplus))

      XH2O2g  = 1.d0 / ( 1.d0 + ( FHCH2O2 * R * T * LWC ) )
!      H2O2aq  = FHCH2O2 * XH2O2g * H2O2 * P

      ! Conversion rate from SO2 to SO4 via reaction with H2O2
      KaqH2O2  = kh2o2 * Ks1 * FHCH2O2 * HCSO2 * XH2O2g * XSO2g
     &         * P * LWC * R * T            ! [v/v/s]

      !=================================================================
      !  Aqueous reactions of SO2 with O3
      !  SO2(aq) + O3 =>                       (0)
      !  HSO3-   + O3 => SO4-- + H+ + O2       (1)
      !  SO3--   + O3 => SO4-- + O2            (2)
      !
      ! NOTE
      ! [Jacob, 1986]
      !    KO1  = 3.49E12 * EXP( -4.83E3 / T )  
      !    KO2  = 7.32E14 * EXP( -4.03E3 / T )
      !
      ! [Jacobson, 1999]
      !    KO0  = 2.40E+4
      !    KO1  = 3.70E+5 * EXP( -18.56 * ((298.15/T) - 1.))
      !    KO2  = 1.50E+9 * EXP( -17.72 * ((298.15/T) - 1.))
      !
      ! Rate constants from Jacobson is larger than those of Jacob
      ! and results in faster conversion from S(IV) to S(VI)
      ! We choose Jacob 1) 2) and Jacobson 0) here
      !=================================================================
      KO0 = 2.40d+4
      KO1 = 3.49d12 * EXP( -4.83d3 / T )  
      KO2 = 7.32d14 * EXP( -4.03d3 / T )

      !=================================================================
      ! H2O2 equilibrium reaction
      ! O3 + H2O = O3.H2O
      !  HCO3  = 1.13E-2 * EXP( 8.51 * (298.15/T -1.) ), S & P
      !  HCO3  = 1.13E-2 * EXP( 7.72 * (298.15/T -1.) ), Jacobson
      !=================================================================

      ! Calculate Henry's Law constant for atmospheric temperature
      HCO3  = 1.13d-2 * EXP( 8.51d0 * ( 298.15d0 / T - 1.d0 ) )

      XO3g  = 1.d0 / ( 1.d0 + ( HCO3 * R * T * LWC ) )
!      O3aq  = HCO3 * XO3g * O3 * P
      
      ! Conversion rate from SO2 to SO4 via reaction with O3
      KaqO3 = (KO0*XSO2AQ + KO1*XHSO3 + KO2*XSO3) * FHCSO2 * XSO2g
     &      * P * HCO3 * XO3g * LWC * R * T   ! [v/v/s]

      END SUBROUTINE AQCHEM_SO2
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: het_drop_chem
!
! !DESCRIPTION: Subroutine HET\_DROP\_CHEM estimates the in-cloud sulfate
!  production rate in heterogeneous cloud droplets based on the Yuen et al.,
!  1996 parameterization. (bec, 6/16/11)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE HET_DROP_CHEM( I, J, L, LSTOT, SSCvv, aSO4, GNH3,
     &                          SO2_sr, H2O20, GNO3, SR, State_Met )
!
! !USES:
!
      USE ERROR_MOD,          ONLY : IT_IS_FINITE, GEOS_CHEM_STOP
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE TIME_MOD,           ONLY : GET_TS_CHEM
      USE TRACER_MOD,         ONLY : TCVV
      USE TRACERID_MOD,       ONLY : IDTSO4, IDTSALC
      USE TRACERID_MOD,       ONLY : IDTDST2, IDTDST3, IDTDST4
!
! !INPUT PARAMETERS: 
!
      INTEGER,        INTENT(IN)  :: I, J, L
      REAL*8,         INTENT(IN)  :: LSTOT
      REAL*8,         INTENT(IN)  :: SSCvv
      REAL*8,         INTENT(IN)  :: aSO4
      REAL*8,         INTENT(IN)  :: GNH3
      REAL*8,         INTENT(IN)  :: SO2_sr
      REAL*8,         INTENT(IN)  :: H2O20
      REAL*8,         INTENT(IN)  :: GNO3
      TYPE(MetState), INTENT(IN)  :: State_Met   ! Meteorology State object
!
! !OUTPUT PARAMETERS: 
!
      REAL*8,         INTENT(OUT) :: SR ! Sulfate production rate
!
! !REVISION HISTORY:
!  25 Jan 2012 - M. Payer    - Added ProTeX headers
!  05 Sep 2013 - M. Sulprizio- Now pass met fields using the State_Met object
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
      REAL*8, PARAMETER   :: PI    = 3.14159D0
      REAL*8, PARAMETER   :: SS_DEN = 2200.d0 !dry sea-salt density [kg/m3]
      ! sigma of the size distribution for sea-salt (Jaegle et al., 2011)
      REAL*8, PARAMETER   :: SIG_S = 1.8d0
      ! geometric dry mean diameters [m] for computing lognormal size distribution
      REAL*8, PARAMETER   :: RG_S = 0.4d-6 !(Jaegle et a., 2011)
      REAL*8, PARAMETER   :: RG_D2 = 1.5d-6!(Ginoux et al., 2001)
      REAL*8, PARAMETER   :: RG_D3 = 2.5d-6
      REAL*8, PARAMETER   :: RG_D4 = 4.d-6
!
! !LOCAL VARIABLES:
!
      REAL*8              :: alpha_NH3, alpha_SO2, alpha_H2O2
      REAL*8              :: alpha_HNO3, alpha_B, alpha_CN
      REAL*8              :: alpha_W, alpha_SO4, sum_gas, H
      REAL*8              :: NDss, NDd2, NDd3, NDd4, CN, DEN, REFF, W
      REAL*8              :: DTCHEM, APV, DSVI
      REAL*8              :: B, NH3, SO2, H2O2, HNO3, SO4
      REAL*8              :: CNss, CNd2, CNd3, CNd4

      ! Pointers
      REAL*8, POINTER     :: AD(:,:,:)
      REAL*8, POINTER     :: AIRDEN(:,:,:)
      REAL*8, POINTER     :: AIRVOL(:,:,:)
      REAL*8, POINTER     :: U(:,:,:)
      REAL*8, POINTER     :: V(:,:,:)

      !=================================================================
      ! HET_DROP_CHEM begins here!
      !=================================================================

      ! Initialize pointers
      AD     => State_Met%AD
      AIRDEN => State_Met%AIRDEN
      AIRVOL => State_Met%AIRVOL
      U      => State_Met%U
      V      => State_Met%V

      ! Convert gas phase concentrations from [v/v] to [pptv]
      NH3  = GNH3  * 1.0d12
      SO2  = SO2_sr  * 1.0d12
      H2O2 = H2O20 * 1.0d12
      HNO3 = GNO3 * 1.0d12

      ! Convert sulfate aerosol concentrations from [v/v] to [ug/m3]
      SO4 = ( aSO4         * AD(I,J,L)     * 1.0d9 ) / 
     &      ( TCVV(IDTSO4) * AIRVOL(I,J,L) )

      ! Convert in cloud sulfate production rate from [v/v/timestep] to
      ! [ug/m3/timestep]
      B  = ( LSTOT        * AD(I,J,L)     * 1.0d9 ) / 
     &     ( TCVV(IDTSO4) * AIRVOL(I,J,L) )

      ! Convert coarse-mode aerosol concentrations from [v/v] to [#/cm3]
      ! based on equation in Hofmann, Science, 1990.
      ! First convert from [v/v] to [kg/m3 air]
      CNss = SSCvv         * AD(I,J,L)     /
     &     ( TCVV(IDTSALC) * AIRVOL(I,J,L) )

      ! Now convert from [kg/m3 air] to [#/cm3 air]
      ! Sea-salt
      NDss = ( (3.d0/4.d0) * CNss ) / 
     &       (PI * SS_DEN * RG_S**3.d0 * exp( (9.d0/2.d0) * 
     &       (LOG(SIG_S)) ** 2.d0 ) )

      ! Total coarse mode number concentration [#/cm3]
      CN = NDss ! sea-salt 

      ! Determine regression coefficients based on the local SO2 concentration
      IF ( SO2 <= 200.0d0 ) THEN
         alpha_B    = 0.5318d0
         alpha_NH3  = -1.67d-7
         alpha_SO2  = 2.59d-6
         alpha_H2O2 = -1.77d-7
         alpha_HNO3 = -1.72d-7
         alpha_W    = 1.22d-6
         alpha_CN   = 4.58d-6
         alpha_SO4  = -1.00d-5
      ELSE IF ( SO2 > 200.0d0 .AND. SO2 <= 500.0d0 ) THEN
         alpha_B    = 0.5591d0
         alpha_NH3  = 3.62d-6
         alpha_SO2  = 1.66d-6
         alpha_H2O2 = 1.06d-7
         alpha_HNO3 = -5.45d-7
         alpha_W    = -5.79d-7
         alpha_CN   = 1.63d-5
         alpha_SO4  = -7.40d-6
      ELSE IF ( SO2 > 500.0d0 .AND. SO2 < 1000.0d0 ) THEN
         alpha_B    = 1.1547d0
         alpha_NH3  = -4.28d-8
         alpha_SO2  = -1.23d-7
         alpha_H2O2 = -9.05d-7
         alpha_HNO3 = 1.73d-7
         alpha_W    = 7.22d-6
         alpha_CN   = 2.44d-5
         alpha_SO4  = 3.25d-5
      ELSE IF ( SO2 >= 1000.0d0 ) THEN
         alpha_B    = 1.1795d0
         alpha_NH3  = 2.57d-7
         alpha_SO2  = -5.54d-7
         alpha_H2O2 = -1.08d-6
         alpha_HNO3 = 1.95d-6
         alpha_W    = 6.14d-6
         alpha_CN   = 1.64d-5
         alpha_SO4  = 2.48d-6
      ENDIF

      ! Updraft velocity over the oceans [cm/s]
      W = 500d0

      ! Compute H (integration time interval * air parcel velocity) [m]
      ! DTCHEM is the chemistry timestep in seconds
      DTCHEM = GET_TS_CHEM() * 60d0

      ! Compute air parcel velocity [m/s]
      APV = SQRT( (U(I,J,L) * U(I,J,L)) + (V(I,J,L) * V(I,J,L)) )

      H   = DTCHEM * APV          ![m]

      sum_gas = (alpha_NH3 * NH3) + (alpha_SO2 * SO2) +
     &          (alpha_H2O2 * H2O2) + (alpha_HNO3 * HNO3)

      DSVI = ( alpha_B * B ) +
     &       ( ( (alpha_CN * CN) + (alpha_W * W) + (alpha_SO4 * SO4) +
     &           sum_gas ) * H )

      ! additional sulfate production that can be attributed to
      ! ozone [ug/m3/timestep]
      SR = DSVI - B 

      ! Convert SR from [ug/m3/timestep] to [v/v/timestep]
      SR = SR *TCVV(IDTSO4) * 1.D-9 / AIRDEN(L,I,J)
	
      ! Don't allow SR to be negative
      SR = MAX( SR, 0.d0 )

      ! Don't produce more SO4 than SO2 available after AQCHEM_SO2
      SR = MIN( SR, SO2_sr )

      ! Free pointers
      NULLIFY( AD     )
      NULLIFY( AIRDEN )
      NULLIFY( AIRVOL )
      NULLIFY( U      )
      NULLIFY( V      )

      END SUBROUTINE HET_DROP_CHEM
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: chem_so4
!
! !DESCRIPTION: Subroutine CHEM\_SO4 is the SO4 chemistry subroutine from Mian
!  Chin's GOCART model, modified for the GEOS-CHEM model.  Now also modified to
!  account for production of crystalline and aqueous sulfur tracers. 
!  (rjp, bdf, cas, bmy, 5/31/00, 5/23/06) 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CHEM_SO4( am_I_Root, Input_Opt, 
     &                     State_Met, State_Chm, RC )
!
! !USES:
!
      USE CMN_SIZE_MOD
      USE CMN_DIAG_MOD
      USE DIAG_MOD,           ONLY : AD44
      USE DRYDEP_MOD,         ONLY : DEPSAV
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE GIGC_State_Chm_Mod, ONLY : ChmState
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE GRID_MOD,           ONLY : GET_AREA_CM2
      USE PBL_MIX_MOD,        ONLY : GET_FRAC_UNDER_PBLTOP
      USE TIME_MOD,           ONLY : GET_TS_CHEM
      USE TRACERID_MOD,       ONLY : IDTSO4
      USE TRACERID_MOD,       ONLY : IDTSO4s
      USE TRACERID_MOD,       ONLY : IDTAS
      USE TRACERID_MOD,       ONLY : IDTAHS 
      USE TRACERID_MOD,       ONLY : IDTLET
      USE TRACERID_MOD,       ONLY : IDTSO4aq
      USE TRACERID_MOD,       ONLY : IDTNH4aq
      USE TROPOPAUSE_MOD,     ONLY : ITS_IN_THE_STRAT
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)    :: am_I_Root   ! Is this the root CPU?
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
!  Reaction List (by Mian Chin, chin@rondo.gsfc.nasa.gov)                  
!  ============================================================================
!  The Only production is from SO2 oxidation (save in CHEM_SO2), and the only  
!  loss is dry depsition here.  Wet deposition will be treated in "wetdep.f".
!                                                                             .
!  SO4 = SO4_0 * exp(-kt) + PSO4_SO2/kt * (1.-exp(-kt))                    
!    where k = dry deposition.                                             
!                                                                             .
! !REVISION HISTORY:
!  (1 ) Now reference AD from "dao_mod.f".  Added parallel DO-loops.  
!        Updated comments, cosmetic changes. (rjp, bdf, bmy, 9/16/02)
!  (2 ) Now replace DXYP(JREF)*1d4 with routine GET_AREA_CM2 of "grid_mod.f"
!        Now use function GET_TS_CHEM from "time_mod.f" (bmy, 3/27/03)
!  (3 ) Now reference PBLFRAC from "drydep_mod.f".  Now apply dry deposition
!        to the entire PBL. (rjp, bmy, 8/1/03)
!  (4 ) Now use ND44_TMP array to store vertical levels of drydep flux, then 
!        sum into AD44 array.  This preents numerical differences when using
!        multiple processors. (bmy, 3/24/04)
!  (5 ) Now use parallel DO-loop to zero ND44_TMP (bmy, 4/14/04)
!  (6 ) Now reference STT & TCVV from "tracer_mod.f" (bmy, 7/20/04)
!  (7 ) Now references LCRYST from "logical_mod.f".  Modified for crystalline
!        and aqueous sulfate2 tracers: AS, AHS, LET, SO4aq.  Also changed name
!        of ND44_TMP to T44 to save space. (cas, bmy, 12/21/04)
!  (8 ) Replace PBLFRAC from "drydep_mod.f" with GET_FRAC_UNDER_PBLTOP from 
!        "pbl_mix_mod.f" (bmy, 2/22/05)
!  (9 )  Now remove reference to CMN, it's obsolete.  Now reference 
!         ITS_IN_THE_STRAT from "tropopause_mod.f" (bmy, 8/22/05)
!  (10) Now references XNUMOL from "tracer_mod.f" (bmy, 10/25/05)
!  (11) Rearrange error check to avoid SEG FAULTS (bmy, 5/23/06)
!  22 Dec 2011 - M. Payer    - Added ProTeX headers
!  01 Mar 2012 - R. Yantosca - Now use GET_AREA_CM2(I,J,L) from grid_mod.F90
!  30 Jul 2012 - R. Yantosca - Now accept am_I_Root as an argument when
!                              running with the traditional driver main.F
!  31 Jul 2012 - R. Yantosca - Now loop from 1..LLPAR for GIGC compatibility
!  31 Jul 2012 - R. Yantosca - Declare temp drydep arrays w/ LLPAR (not LLTROP)
!  14 Nov 2012 - R. Yantosca - Add am_I_Root, Input_Opt, RC as arguments
!  15 Nov 2012 - M. Payer    - Replaced all met field arrays with State_Met
!                              derived type object
!  05 Mar 2013 - R. Yantosca - Now use Input_Opt%LNLPBL
!  19 Mar 2013 - R. Yantosca - Now copy Input_Opt%TCVV(1:N_TRACERS) and
!                              Input_Opt%XNUMOL(1:N_TRACERS) -- avoid OOB errs
!  25 Mar 2013 - M. Payer    - Now pass State_Chm object via the arg list
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL               :: LNLPBL, LCRYST, LSSALT
      INTEGER               :: N_TRACERS
      INTEGER               :: I,      J,      L,        N,     N_ND44
      REAL*8                :: AS,     AS0,    AHS,      AHS0,  LET   
      REAL*8                :: LET0,   SO4,    SO40,     SO4aq, SO4aq0
      REAL*8                :: SO4s,   SO40s,  RKT,      RKTs,  E_RKT
      REAL*8                :: E_RKTs, DTCHEM, AREA_CM2, FLUX

      ! Arrays
      REAL*8                :: TCVV  (Input_Opt%N_TRACERS)
      REAL*8                :: XNUMOL(Input_Opt%N_TRACERS)
      REAL*8                :: T44(IIPAR,JJPAR,LLPAR,6) 

      ! Pointers
      REAL*8, POINTER       :: STT(:,:,:,:)
      REAL*8, POINTER       :: AD(:,:,:)

      !=================================================================
      ! CHEM_SO4 begins here!
      !=================================================================

      ! Return if tracers are not defined
      IF ( IDTSO4 == 0 .or. IDTSO4s == 0 ) RETURN
      IF ( DRYSO4 == 0 .or. DRYSO4s == 0 ) RETURN

      ! Assume success
      RC        = GIGC_SUCCESS

      ! Copy fields from INPUT_OPT to local variables for use below
      LCRYST    = Input_Opt%LCRYST
      LNLPBL    = Input_Opt%LNLPBL
      LSSALT    = Input_Opt%LSSALT
      N_TRACERS = Input_Opt%N_TRACERS
      TCVV      = Input_Opt%TCVV  (1:N_TRACERS) 
      XNUMOL    = Input_Opt%XNUMOL(1:N_TRACERS)

      ! Initialize pointers
      STT      => State_Chm%Tracers
      AD       => State_Met%AD

      ! DTCHEM is the chemistry timestep in seconds
      DTCHEM    = GET_TS_CHEM() * 60d0

      ! Number of drydep tracers to save
      N_ND44    = 2

!------------------------------------------------------------------------------
!%%% Currently under development (rjp, bmy, 4/13/05)
!%%%       ! Number of drydep tracers to save
!%%%       IF ( LCRYST ) THEN 
!%%%          N_ND44 = 6
!%%%       ELSE
!%%%          N_ND44 = 2
!%%%       ENDIF
!------------------------------------------------------------------------------

      ! Zero T44 array
      IF ( ND44 > 0 ) THEN
         T44 = 0d0
      ENDIF

      ! Loop over tropospheric grid boxes
!------------------------------------------------------------------------------
!%%% Currently under development (rjp, bmy, 3/15/05)
!%%% $OMP PARALLEL DO
!%%% $OMP+DEFAULT( SHARED )
!%%% $OMP+PRIVATE( I,     J,      L,    AREA_CM2, RKT,   RKTs   )
!%%% $OMP+PRIVATE( E_RKT, E_RKTs, FLUX, SO4,      SO4s,  SO4aq  )
!%%% $OMP+PRIVATE( AS,    AHS,    LET,  SO40,     SO40s, SO4aq0 )
!%%% $OMP+PRIVATE( AS0,   AHS0,   LET0                          )
!%%% $OMP+SCHEDULE( DYNAMIC ) 
!------------------------------------------------------------------------------
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I,      J,    L,   AREA_CM2, RKT,   RKTs, E_RKT )
!$OMP+PRIVATE( E_RKTs, FLUX, SO4, SO4s,     SO40,  SO40s       )
!$OMP+SCHEDULE( DYNAMIC ) 
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Initialize for safety's sake
         AREA_CM2 = 0d0
         RKT      = 0d0
         RKTs     = 0d0
         E_RKT    = 0d0
         E_RKTs   = 0d0
         FLUX     = 0d0
         SO4      = 0d0
         SO4s     = 0d0
!------------------------------------------------------------------------------
!%%% Currently under development (rjp, bmy, 3/15/05)
!%%%         SO4aq    = 0d0
!%%%         AS       = 0d0
!%%%         AHS      = 0d0
!%%%         LET      = 0d0
!------------------------------------------------------------------------------

         ! Skip stratospheric boxes
         IF ( ITS_IN_THE_STRAT( I, J, L, State_Met ) ) CYCLE

         !==============================================================
         ! Initial concentrations before chemistry
         !==============================================================

         ! SO4 [v/v]
         SO40  = STT(I,J,L,IDTSO4)
         
         ! SO4 within coarse seasalt aerosol [v/v]
         SO40s = STT(I,J,L,IDTSO4s)

!------------------------------------------------------------------------------
!%%% Currently under development (rjp, bmy, 3/15/05)
!%%%          ! SO4aq, AS, LET, AHS (if necessary) [v/v]
!%%%          IF ( LCRYST ) THEN
!%%%             SO4aq0 = STT(I,J,L,IDTSO4aq)
!%%%             AS0    = STT(I,J,L,IDTAS)
!%%%             AHS0   = STT(I,J,L,IDTAHS)
!%%%             LET0   = STT(I,J,L,IDTLET)          
!%%%          ENDIF
!------------------------------------------------------------------------------

         !==============================================================
         ! SO4 chemistry: 
         !
         ! (CASE 1) SO4 production from SO2 and loss by drydep
         !          --> see equation in header notes above
         !
         ! (CASE 2) SO4 production from SO2 with no SO4 loss by drydep
         !==============================================================

         ! SO4 drydep frequency [1/s].  Also accounts for the fraction
         ! of each vertical level that is located below the PBL top
         RKT  = DEPSAV(I,J,DRYSO4)  * GET_FRAC_UNDER_PBLTOP( I, J, L )

         ! Add option for non-local PBL (Lin, 03/31/09) 
         IF ( LNLPBL ) RKT = 0.D0

         ! RKT > 0 denotes that we have SO4 drydep occurring
         IF ( RKT > 0d0 ) THEN
            
            !-----------------------------------------------------------
            ! CASE 1: SO4 production from SO2 and SO4 loss by drydep
            !-----------------------------------------------------------

            ! Fraction of SO4 lost to drydep [unitless]
            RKT   = RKT * DTCHEM
            
            ! Pre-compute exponential term for use below
            E_RKT = EXP( -RKT ) 

            ! Updated SO4 (gas phase) [v/v]
            SO4   = ( SO40                *          E_RKT ) + 
     &              ( PSO4_SO2(I,J,L)/RKT * ( 1.d0 - E_RKT ) )

!------------------------------------------------------------------------------
!%%% Currently under development (rjp, bmy, 3/15/05)
!%%%             IF ( LCRYST ) THEN
!%%% 
!%%%                ! Updated SO4 (aqueous phase) [v/v]
!%%%                SO4aq = ( SO4aq0              *          E_RKT   ) + 
!%%%      &                 ( PSO4_SO2(I,J,L)/RKT * ( 1.d0 - E_RKT ) ) 
!%%% 
!%%%                ! Updated AS, AHS, LET [v/v] 
!%%%                AS    = AS0  * E_RKT
!%%%                AHS   = AHS0 * E_RKT
!%%%                LET   = LET0 * E_RKT
!%%% 
!%%%             ENDIF
!------------------------------------------------------------------------------

         ELSE

            !-----------------------------------------------------------
            ! CASE 2: Production of SO4 from SO2; no SO4 drydep loss
            !-----------------------------------------------------------

!-----------------------------------------------------------------------------
!%%% Currently under development (rjp, bmy, 3/15/05)
!%%%            IF ( LCRYST ) THEN
!%%%
!%%%               ! SO4 production from SO2 (both gas-phase & aqueous)
!%%%               SO4   = SO40   + ( 2d0 * PSO4_SO2(I,J,L) )
!%%%               SO4aq = SO4aq0 + ( 2d0 * PSO4_SO2(I,J,L) )
!%%%
!%%%               ! No production from AS, AHS, LET
!%%%               AS    = AS0
!%%%               AHS   = AHS0
!%%%               LET   = LET0
!%%%
!%%%            ELSE
!%%%
!%%%               ! SO4 production from SO2 [v/v/timestep]
!%%%               SO4 = SO40 + PSO4_SO2(I,J,L)
!%%%
!%%%            ENDIF
!-----------------------------------------------------------------------------

               ! SO4 production from SO2 [v/v/timestep]
               SO4 = SO40 + PSO4_SO2(I,J,L)

         ENDIF

         !==============================================================
         ! SO4s (SO4 w/in seasalt aerosol) chemistry: 
         !
         ! (CASE 3) SO4s production from seasalt and loss by drydep
         !          --> see equation in header notes above
         !
         ! (CASE 4) SO4s prod from seasalt w/ no SO4s loss by drydep
         !==============================================================

         ! SO4s drydep frequency [1/s].   Also accounts for the fraction
         ! of each vertical level that is located below the PBL top
         RKTs = DEPSAV(I,J,DRYSO4s) * GET_FRAC_UNDER_PBLTOP( I, J, L )

         ! Add option for non-local PBL (Lin, 03/31/09)
         IF ( LNLPBL ) RKTs = 0.D0
        
         ! RKTs > 0 indicates that SO4s drydep occurs
         IF ( RKTs > 0d0 ) THEN
            
            !-----------------------------------------------------------
            ! CASE 3: SO4s prod from seasalt SO4s loss by drydep
            !-----------------------------------------------------------

            ! Fraction of SO4s lost to drydep [unitless]
            RKTs   = RKTs * DTCHEM
            
            ! Pre-compute exponential term for use below
            E_RKTs = EXP( -RKTs ) 
               
            ! Updated SO4 (gas phase) [v/v]
            SO4s   = ( SO40s               *          E_RKTs ) + 
     &               ( PSO4_ss(I,J,L)/RKTs * ( 1.d0 - E_RKTs ) )

         ELSE

            !--------------------------------------------------------
            ! CASE 4: Prod of SO4s from seasalt; no SO4s drydep loss
            !--------------------------------------------------------

            ! SO4 production from SO2 [v/v/timestep]
            SO4s = SO40s + PSO4_ss(I,J,L)

         ENDIF

         !==============================================================
         ! Final concentrations after chemistry
         !==============================================================

         ! Error check
         IF ( SO4  < SMALLNUM ) SO4  = 0d0
         IF ( SO4s < SMALLNUM ) SO4s = 0d0

         ! Final concentrations [v/v]
         STT(I,J,L,IDTSO4)  = SO4
         STT(I,J,L,IDTSO4s) = SO4s

!-----------------------------------------------------------------------------
!%%% Currently under development (bmy, 3/15/05)
!%%%          ! SO4aq, AS, AHS, LET (if necessary)
!%%%          IF ( LCRYST ) THEN
!%%% 
!%%%             ! Error check
!%%%             IF ( SO4aq < SMALLNUM ) SO4aq = 0d0
!%%%             IF ( AS    < SMALLNUM ) AS    = 0d0
!%%%             IF ( AHS   < SMALLNUM ) AHS   = 0d0
!%%%             IF ( LET   < SMALLNUM ) LET   = 0d0
!%%% 
!%%%             ! Final SO4aq, AS, AHS, LET [v/v]
!%%%             STT(I,J,L,IDTSO4aq) = SO4aq
!%%%             STT(I,J,L,IDTAS)    = AS
!%%%             STT(I,J,L,IDTAHS)   = AHS
!%%%             STT(I,J,L,IDTLET)   = LET
!%%% 
!%%%          ENDIF
!-----------------------------------------------------------------------------

         !==============================================================
         ! ND44 Diagnostic: Drydep flux of SO4 and the crystalline & 
         ! aqueous tracers (AS, AHS, LET, SO4aq) in [molec/cm2/s]
         !==============================================================
         IF ( ND44 > 0 ) THEN

            ! Surface area [cm2]
            AREA_CM2 = GET_AREA_CM2( I, J, L )

            ! SO4 drydep flux [molec/cm2/s]
            FLUX = SO40  - SO4 + PSO4_SO2(I,J,L) 
            FLUX = FLUX  * AD(I,J,L) / TCVV(IDTSO4)
            FLUX = FLUX  * XNUMOL(IDTSO4)      / AREA_CM2 / DTCHEM
            T44(I,J,L,1) = T44(I,J,L,1) + FLUX

            ! SO4s drydep flux [molec/cm2/s]
            FLUX = SO40s - SO4s + PSO4_ss(I,J,L) 
            FLUX = FLUX  * AD(I,J,L) / TCVV(IDTSO4s)
            FLUX = FLUX  * XNUMOL(IDTSO4s)     / AREA_CM2 / DTCHEM
            T44(I,J,L,2) = T44(I,J,L,2) + FLUX

!------------------------------------------------------------------------------
!%%% Currently under development (rjp, bmy, 3/15/05)
!%%%             ! SO4aq, AS, AHS, LET drydep fluxes (if necessary)
!%%%             IF ( LCRYST ) THEN
!%%% 
!%%%                ! SO4aq drydep flux [molec/cm2/s]
!%%%                FLUX = SO4aq0 - SO4aq + PSO4_SO2(I,J,L)     
!%%%                FLUX = FLUX  * AD(I,J,L) / TCVV(IDTSO4aq)
!%%%                FLUX = FLUX  * XNUMOL(IDTSO4aq)    / AREA_CM2 / DTCHEM
!%%%                T44(I,J,L,3) = T44(I,J,L,3) + FLUX
!%%% 
!%%%                ! AS drydep flux [molec/cm2/s]
!%%%                FLUX = AS0   - AS
!%%%                FLUX = FLUX  * AD(I,J,L) / TCVV(IDTAS)
!%%%                FLUX = FLUX  * XNUMOL(IDTAS)       / AREA_CM2 / DTCHEM
!%%%                T44(I,J,L,4) = T44(I,J,L,4) + FLUX
!%%% 
!%%%                ! AHS drydep flux [molec/cm2/s]
!%%%                FLUX = AHS0  - AHS
!%%%                FLUX = FLUX  * AD(I,J,L) / TCVV(IDTAHS)
!%%%                FLUX = FLUX  * XNUMOL(IDTAHS)      / AREA_CM2 / DTCHEM
!%%%                T44(I,J,L,5) = T44(I,J,L,5) + FLUX
!%%% 
!%%%                ! LET drydep flux [molec/cm2/s]
!%%%                FLUX = LET0  - LET
!%%%                FLUX = FLUX  * AD(I,J,L) / TCVV(IDTLET)
!%%%                FLUX = FLUX  * XNUMOL(IDTLET)      / AREA_CM2 / DTCHEM
!%%%                T44(I,J,L,6) = T44(I,J,L,6) + FLUX
!%%% 
!%%%             ENDIF
!------------------------------------------------------------------------------
         ENDIF
      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      !===============================================================
      ! ND44: Sum drydep fluxes by level into the AD44 array in
      ! order to ensure that we get the same results w/ sp or mp 
      !===============================================================
      IF ( ND44 > 0 ) THEN 
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
         DO J = 1, JJPAR
         DO I = 1, IIPAR
         DO L = 1, LLPAR

            ! Sum SO4, SO4s drydep fluxes in the vertical [molec/cm2/s]
            AD44(I,J,DRYSO4, 1) = AD44(I,J,DRYSO4, 1) + T44(I,J,L,1)
            AD44(I,J,DRYSO4s,1) = AD44(I,J,DRYSO4s,1) + T44(I,J,L,2)

!------------------------------------------------------------------------------
!%%% Currently under development (rjp, bmy, 3/15/05)
!%%%             ! Sum SO4aq, AS, AHS, LET drydep fluxes (if necessary)
!%%%             IF ( LCRYST ) THEN
!%%%                AD44(I,J,DRYSO4aq,1) = AD44(I,J,DRYSO4aq,1)+T44(I,J,L,3)
!%%%                AD44(I,J,DRYAS,   1) = AD44(I,J,DRYAS,   1)+T44(I,J,L,4)
!%%%                AD44(I,J,DRYAHS,  1) = AD44(I,J,DRYAHS,  1)+T44(I,J,L,5)
!%%%                AD44(I,J,DRYLET,  1) = AD44(I,J,DRYLET,  1)+T44(I,J,L,6)
!%%%             ENDIF
!------------------------------------------------------------------------------

         ENDDO
         ENDDO
         ENDDO
!$OMP END PARALLEL DO
      ENDIF

      ! Free pointers
      NULLIFY( STT )
      NULLIFY( AD  )

      END SUBROUTINE CHEM_SO4
!EOC
#if   defined( TOMAS )
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: chem_so4_aq
!
! !DESCRIPTION: Subroutine CHEM\_SO4\_AQ takes the SO4 produced via aqueous 
!  chemistry of SO2 and distribute onto the size-resolved aerosol number and 
!  sulfate mass as a part of the TOMAS aerosol microphysics module 
!  (win, 1/25/10)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CHEM_SO4_AQ( am_I_Root, Input_Opt, 
     &                        State_Met, State_Chm, RC )
!
! !USES:
!
      USE CMN_SIZE_MOD
      USE DAO_MOD,            ONLY : CONVERT_UNITS 
      USE GIGC_Input_Opt_Mod, ONLY : OptINput 
      USE GIGC_State_Chm_Mod, ONLY : ChmState
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE PBL_MIX_MOD,        ONLY : GET_FRAC_UNDER_PBLTOP
      USE TOMAS_MOD,          ONLY : AQOXID, GETACTBIN
      USE TRACERID_MOD,       ONLY : IDTSO4, IDTNK10
      USE TROPOPAUSE_MOD,     ONLY : ITS_IN_THE_STRAT
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
!  NOTE: This subroutine is ignored unless we compile for TOMAS microphysics.
!
! !REVISION HISTORY:
!  (1 ) As of now the SO4 produced via heterogeneous reaction on the 2-mode
!        seasalt is not include in this treatment (win, 7/23/07)
!  (2 ) Change a fixed kmin = 8 (corresponding to the assumed activation dia. 
!        of 55nm to be varying with current chemical composition.  Take average
!        of the activating bin for LS and CONV rains. (win, 9/25/07)
!  16 Feb 2012 - R. Yantosca - Added ProTeX headers
!  09 Nov 2012 - M. Payer    - Replaced all met field arrays with State_Met
!                              derived type object
!  31 May 2013 - R. Yantosca - Now accept am_I_Root, Input_Opt, State_Chm
!                              and RC arguments
!  31 May 2013 - R. Yantosca - Now pass State_Chm to TOMAS routines
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER :: I, J, L
      integer :: k, binact1, binact2
      INTEGER :: KMIN
      REAL*8  :: SO4OXID

      ! For values from Input_Opt
      INTEGER :: N_TRACERS
      REAL*8  :: TCVV(Input_Opt%N_TRACERS)
      REAL*8  :: XNUMOL(Input_Opt%N_TRACERS)

      ! Pointers
      REAL*8, POINTER :: AD(:,:,:)

      !=================================================================
      ! CHEM_SO4_AQ begins here!
      !=================================================================
 
      ! Copy values from Input_Opt
      N_TRACERS = Input_Opt%N_TRACERS
      TCVV      = Input_Opt%TCVV  (1:N_TRACERS)
      XNUMOL    = Input_Opt%XNUMOL(1:N_TRACERS)

      ! Initialize pointers
      AD     => State_Met%AD

      ! Convert STT from [v/v] -> [kg] 
      CALL CONVERT_UNITS( 2, N_TRACERS, TCVV, AD, State_Chm%Tracers )

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
!$OMP+PRIVATE( KMIN, SO4OXID, BINACT1, BINACT2 )
!$OMP+SCHEDULE( DYNAMIC )
      DO L = 1, LLTROP
      DO J = 1, JJPAR
      DO I = 1, IIPAR
                  
         ! Skip stratospheric boxes
         IF (ITS_IN_THE_STRAT( I, J, L, State_Met ) ) CYCLE

         SO4OXID = PSO4_SO2AQ(I,J,L) * AD(I,J,L)/ TCVV(IDTSO4)

         IF ( SO4OXID > 0d0 ) THEN

            ! The activating bin of the composition of BIN10
            ! is, in most any case, the same as activating bin of other bins.
            CALL GETACTBIN
     &         ( I, J, L, IDTNK10, .TRUE. , BINACT1, State_Chm )

            CALL GETACTBIN
     &         ( I, J, L, IDTNK10, .FALSE., BINACT2, State_Chm )

            KMIN = ( BINACT1 + BINACT2 )/ 2.

            CALL AQOXID( SO4OXID, KMIN, I, J, L, State_Met, State_Chm )
         ENDIF
      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

       ! Convert STT from [kg] -> [v/v] 
      CALL CONVERT_UNITS( 1, N_TRACERS, TCVV, AD, State_Chm%Tracers )
      
      ! Free pointers
      NULLIFY( AD )

      END SUBROUTINE CHEM_SO4_AQ
!EOC
#endif
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: chem_msa
!
! !DESCRIPTION: Subroutine CHEM\_MSA is the SO4 chemistry subroutine from Mian
!  Chin's GOCART model, modified for the GEOS-CHEM model. (rjp, bdf, bmy,
!  5/31/00, 10/25/05)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CHEM_MSA( am_I_Root, Input_Opt, 
     &                     State_Met, State_Chm, RC )
!
! !USES:
!
      USE CMN_GCTM_MOD
      USE CMN_DIAG_MOD
      USE CMN_SIZE_MOD
      USE DIAG_MOD,           ONLY : AD44
      USE DRYDEP_MOD,         ONLY : DEPSAV
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE GIGC_State_Chm_Mod, ONLY : ChmState
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE GRID_MOD,           ONLY : GET_AREA_CM2
      USE PBL_MIX_MOD,        ONLY : GET_FRAC_UNDER_PBLTOP
      USE PBL_MIX_MOD,        ONLY : GET_PBL_MAX_L
      USE TIME_MOD,           ONLY : GET_TS_CHEM
      USE TRACERID_MOD,       ONLY : IDTMSA
      USE TROPOPAUSE_MOD,     ONLY : ITS_IN_THE_STRAT
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
!  Reaction List (by Mian Chin, chin@rondo.gsfc.nasa.gov)                  
!  ============================================================================
!  The Only production is from DMS oxidation (saved in CHEM_DMS), and the only
!  loss is dry depsition here.  Wet deposition will be treated in "wetdep.f".
!                                                                             .
!  MSA = MSA_0 * exp(-dt) + PMSA_DMS/kt * (1.-exp(-kt))                    
!    where k = dry deposition.                                             
!        
! !REVISION HISTORY:
!  (1 ) Now reference AD from "dao_mod.f".  Added parallel DO-loops.  
!        Updated comments, cosmetic changes. (rjp, bmy, bdf, 9/16/02)
!  (2 ) Now replace DXYP(JREF)*1d4 with routine GET_AREA_CM2 of "grid_mod.f"
!        Now use function GET_TS_CHEM from "time_mod.f" (bmy, 3/27/03)
!  (3 ) Now reference PBLFRAC from "drydep_mod.f".  Now apply dry deposition
!        to the entire PBL. (rjp, bmy, 8/1/03) 
!  (4 ) Now use ND44_TMP array to store vertical levels of drydep flux, then 
!        sum into AD44 array.  This preents numerical differences when using
!        multiple processors. (bmy, 3/24/04) 
!  (5 ) Now use parallel DO-loop to zero ND44_TMP (bmy, 4/14/04)
!  (6 ) Now references STT & TCVV from "tracer_mod.f" (bmy, 7/20/04)
!  (7 ) Replace PBLFRAC from "drydep_mod.f" with GET_FRAC_UNDER_PBLTOP from 
!        "pbl_mix_mod.f".  Also reference GET_PBL_MAX_L from "pbl_mix_mod.f"
!        Vertical DO-loops can run up to PBL_MAX and not LLTROP.   Also
!        remove reference to header file CMN. (bmy, 2/22/05)
!  (8 ) Now references XNUMOL from "tracer_mod.f" (bmy, 10/25/05)
!  (9 ) Change loop back to over entire troposphere to correctly add production
!        of MSA (PMSA_dms) to the MSA tracer array.
!       Added reference USE_TROPOPAUSE_MOD, ONLY : ITS_IN_THE_STRAT
!        as a precaution.  (pjh, 8/19/2009)
!  22 Dec 2011 - M. Payer    - Added ProTeX headers
!  01 Mar 2012 - R. Yantosca - Now use GET_AREA_CM2(I,J,L) from grid_mod.F90
!  31 Jul 2012 - R. Yantosca - Now loop from 1..LLPAR for GIGC compatibility
!  31 Jul 2012 - R. Yantosca - Declare temp drydep arrays w/ LLPAR (not LLTROP)
!  14 Nov 2012 - R. Yantosca - Add am_I_Root, Input_Opt, RC as arguments
!  15 Nov 2012 - M. Payer    - Replaced all met field arrays with State_Met
!                              derived type object
!  05 Mar 2013 - R. Yantosca - Now use Input_Opt%LNLPBL
!  19 Mar 2013 - R. Yantosca - Now copy Input_Opt%TCVV(1:N_TRACERS) and
!                              Input_Opt%XNUMOL(1:N_TRACERS) -- avoid OOB errs
!  25 Mar 2013 - M. Payer    - Now pass State_Chm object via the arg list
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! Scalars
      LOGICAL               :: LNLPBL
      INTEGER               :: N_TRACERS
      INTEGER               :: I,      J,    L,        PBL_MAX
      REAL*8                :: DTCHEM, MSA0, MSA,      RK       
      REAL*8                :: RKT,    FLUX, AREA_CM2, F_UNDER_TOP

      ! Arrays
      REAL*8                :: ND44_TMP(IIPAR,JJPAR,LLPAR)
      REAL*8                :: TCVV  (Input_Opt%N_TRACERS)
      REAL*8                :: XNUMOL(Input_Opt%N_TRACERS)

      ! Pointers
      ! We need to define local arrays to hold corresponding values 
      ! from the Chemistry State (State_Chm) object. (mpayer, 12/6/12)
      REAL*8, POINTER :: STT(:,:,:,:)

      !=================================================================
      ! CHEM_MSA begins here!
      !=================================================================
      IF ( IDTMSA == 0 .or. DRYMSA == 0 ) RETURN

      ! Assume success
      RC        = GIGC_SUCCESS

      ! Copy fields from INPUT_OPT to local variables for use below
      N_TRACERS = Input_Opt%N_TRACERS
      TCVV      = Input_Opt%TCVV  (1:N_TRACERS)
      XNUMOL    = Input_Opt%XNUMOL(1:N_TRACERS)
      LNLPBL    = Input_Opt%LNLPBL

      ! Initialize GEOS-Chem tracer array [kg] from Chemistry State object
      ! (mpayer, 12/6/12)
      STT => State_Chm%Tracers

      ! DTCHEM is the chemistry interval in seconds
      DTCHEM  = GET_TS_CHEM() * 60d0 

      ! Zero ND44_TMP array
      IF ( ND44 > 0 ) THEN
         ND44_TMP = 0d0
      ENDIF

      ! Loop over tropospheric grid boxes
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, F_UNDER_TOP, MSA0, RKT, MSA, AREA_CM2, FLUX )
!$OMP+SCHEDULE( DYNAMIC )
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Skip stratospheric boxes
         IF ( ITS_IN_THE_STRAT( I, J, L, State_Met ) ) CYCLE      

         ! Fraction of box (I,J,L) underneath the PBL top [unitless]
         F_UNDER_TOP = GET_FRAC_UNDER_PBLTOP( I, J, L )

         ! Initial MSA [v/v]
         MSA0 = STT(I,J,L,IDTMSA) 

         ! MSA drydep frequency [1/s].  Also accounts for the fraction
         ! of each grid box (I,J,L) that is located beneath the PBL top
         RKT = DEPSAV(I,J,DRYMSA) * F_UNDER_TOP

         ! Add option for non-local PBL (Lin, 03/31/09)
         IF ( LNLPBL ) RKT = 0.D0
         
         ! RKT > 0 denotes that we have drydep occurring
         IF ( RKT > 0.d0 ) THEN

            ! Fraction of MSA lost to drydep [unitless]
            RKT = RKT * DTCHEM
            
            ! Modified MSA concentration 
            MSA = ( MSA0 * EXP( -RKT )                        ) +
     &            ( PMSA_DMS(I,J,L)/RKT * ( 1d0 - EXP( -RKT ) ) )

         ELSE
            
            ! MSA production from DMS [v/v/timestep]
            MSA = MSA0 + PMSA_DMS(I,J,L)

         ENDIF

         ! Final MSA [v/v]
         IF ( MSA < SMALLNUM ) MSA = 0d0
         STT(I,J,L,IDTMSA) = MSA

         !===========================================================
         ! ND44 Diagnostic: Drydep flux of MSA [molec/cm2/s]
         !===========================================================
         IF ( ND44 > 0 .and. RKT > 0d0 ) THEN
                     
            ! Surface area [cm2]
            AREA_CM2 = GET_AREA_CM2( I, J, L )

            ! Convert [v/v/timestep] to [molec/cm2/s]
            FLUX = MSA0 - MSA + PMSA_DMS(I,J,L)                    
            FLUX = FLUX * State_Met%AD(I,J,L) / TCVV(IDTMSA)            
            FLUX = FLUX * XNUMOL(IDTMSA) / AREA_CM2 / DTCHEM    

            ! Store dryd flux in ND44_TMP as a placeholder
            ND44_TMP(I,J,L) = ND44_TMP(I,J,L) + FLUX
         ENDIF
      
      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO
      
      !===============================================================
      ! ND44: Sum drydep fluxes by level into the AD44 array in
      ! order to ensure that  we get the same results w/ sp or mp 
      !===============================================================
      IF ( ND44 > 0 ) THEN 
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
         DO J = 1, JJPAR
         DO I = 1, IIPAR
         DO L = 1, LLPAR
            AD44(I,J,DRYMSA,1) = AD44(I,J,DRYMSA,1) + ND44_TMP(I,J,L)
         ENDDO
         ENDDO
         ENDDO
!$OMP END PARALLEL DO
      ENDIF

      ! Free pointer
      NULLIFY( STT )

      END SUBROUTINE CHEM_MSA
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: chem_nh3
!
! !DESCRIPTION:  Subroutine CHEM\_NH3 removes NH3 from the surface via dry
!  deposition. (rjp, bdf, bmy, 1/2/02, 10/25/05)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CHEM_NH3( am_I_Root, Input_Opt, 
     &                     State_Met, State_Chm, RC )
!
! !USES:
!
      USE CMN_DIAG_MOD
      USE CMN_SIZE_MOD
      USE DIAG_MOD,           ONLY : AD44
      USE DRYDEP_MOD,         ONLY : DEPSAV
      USE GET_NDEP_MOD,       ONLY : SOIL_DRYDEP
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE GIGC_State_Chm_Mod, ONLY : ChmState
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE GRID_MOD,           ONLY : GET_AREA_CM2
      USE PBL_MIX_MOD,        ONLY : GET_FRAC_UNDER_PBLTOP
      USE PBL_MIX_MOD,        ONLY : GET_PBL_MAX_L
      USE TIME_MOD,           ONLY : GET_TS_CHEM
      USE TRACERID_MOD,       ONLY : IDTNH3
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
!  If you are using the non-local PBL mixing scheme (VDIFF), then routine
!  SOIL_DRYDEP and the ND44 diagnostic updating will be done there.  
!                                                                             .
!  If you are using the full TURBDAY PBL mixing scheme, then we have to 
!  call the SOIL_DRYDEP and archive the ND44 diagnostics here.
!                                                                             .
!  Reaction List:
!  ============================================================================
!  (1 ) NH3 = NH3_0 * EXP( -dt )  where d = dry deposition rate [s-1]
!        
! !REVISION HISTORY:
!  (1 ) Now reference AD from "dao_mod.f".  Added parallel DO-loops.  
!        Updated comments, cosmetic changes. (rjp, bmy, bdf, 9/16/02)
!  (2 ) Now replace DXYP(J+J0)*1d4 with routine GET_AREA_CM2 from "grid_mod.f"
!        Now use function GET_TS_CHEM from "time_mod.f" (bmy, 3/27/03)
!  (3 ) Now reference PBLFRAC from "drydep_mod.f".  Now apply dry deposition
!        to the entire PBL.  Added L and FREQ variables.  Recode to avoid 
!        underflow from the EXP() function. (rjp, bmy, 8/1/03) 
!  (4 ) Now use ND44_TMP array to store vertical levels of drydep flux, then 
!        sum into AD44 array.  This preents numerical differences when using
!        multiple processors. (bmy, 3/24/04)    
!  (5 ) Now use parallel DO-loop to zero ND44_TMP (bmy, 4/14/04)
!  (6 ) Now references STT & TCVV from "tracer_mod.f" Also remove reference to
!        CMN, it's not needed(bmy, 7/20/04)
!  (7 ) Replace PBLFRAC from "drydep_mod.f" with GET_FRAC_UNDER_PBLTOP from 
!        "pbl_mix_mod.f".  Also reference GET_PBL_MAX_L from "pbl_mix_mod.f"
!        Vertical DO-loops can run up to PBL_MAX and not LLTROP. (bmy, 2/22/05)
!  (8 ) Now references XNUMOL from "tracer_mod.f" (bmy, 10/25/05)
!  22 Dec 2011 - M. Payer    - Added ProTeX headers
!  01 Mar 2012 - R. Yantosca - Now use GET_AREA_CM2(I,J,L) from grid_mod.F90
!  31 Jul 2012 - R. Yantosca - Now loop from 1..LLPAR for GIGC compatibility
!  31 Jul 2012 - R. Yantosca - Declare temp drydep arrays w/ LLPAR (not LLTROP)
!  14 Nov 2012 - R. Yantosca - Add am_I_Root, Input_Opt, RC as arguments
!  15 Nov 2012 - M. Payer    - Replaced all met field arrays with State_Met
!                              derived type object
!  05 Mar 2013 - R. Yantosca - Now use Input_Opt%LNLPBL
!  13 Mar 2013 - R. Yantosca - Bug fix: make sure we pass values to the
!                              SOIL_DRYDEP routine even when ND44 is off
!                              (this happens when LNLPBL = F)
!  19 Mar 2013 - R. Yantosca - Now copy Input_Opt%TCVV(1:N_TRACERS) and
!                              Input_Opt%XNUMOL(1:N_TRACERS) -- avoid OOB errs
!  25 Mar 2013 - M. Payer    - Now pass State_Chm object via the arg list
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! Scalars
      LOGICAL :: LNLPBL, LSOILNOX
      INTEGER :: I,      J,        L,    N_TRACERS
      REAL*8  :: DTCHEM, NH30,     NH3
      REAL*8  :: FREQ,   AREA_CM2, FLUX, F_UNDER_TOP

      ! Arrays
      REAL*8  :: DRYD_FLX(IIPAR,JJPAR,LLPAR)
      REAL*8  :: TCVV  (Input_Opt%N_TRACERS)
      REAL*8  :: XNUMOL(Input_Opt%N_TRACERS)

      ! Pointers
      ! We need to define local arrays to hold corresponding values 
      ! from the Chemistry State (State_Chm) object. (mpayer, 12/6/12)
      REAL*8, POINTER :: STT(:,:,:,:)

      !=================================================================
      ! CHEM_NH3 begins here!
      !=================================================================
      IF ( IDTNH3 == 0 .or. DRYNH3 == 0 ) RETURN

      ! Assume success
      RC        = GIGC_SUCCESS

      ! Copy fields from INPUT_OPT to local variables for use below
      LNLPBL    = Input_Opt%LNLPBL
      LSOILNOX  = Input_Opt%LSOILNOX
      N_TRACERS = Input_Opt%N_TRACERS
      TCVV      = Input_Opt%TCVV  (1:N_TRACERS)
      XNUMOL    = Input_Opt%XNUMOL(1:N_TRACERS)

      ! Initialize GEOS-Chem tracer array [kg] from Chemistry State object
      ! (mpayer, 12/6/12)
      STT => State_Chm%Tracers

      ! DTCHEM is the chemistry interval in seconds
      DTCHEM    = GET_TS_CHEM() * 60d0

      ! Zero DRYD_FLX array.  This will store drydep values even if 
      ! the ND44 diagnostic has been turned off (bmy, 3/13/13)
      DRYD_FLX  = 0d0

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, F_UNDER_TOP, FREQ, NH30, NH3, AREA_CM2, FLUX )
!$OMP+SCHEDULE( DYNAMIC )
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Fraction of box (I,J,L) underneath the PBL top [unitless]
         F_UNDER_TOP = GET_FRAC_UNDER_PBLTOP( I, J, L )
         
         ! Only apply drydep to boxes w/in the PBL
         IF ( F_UNDER_TOP > 0d0 ) THEN

            ! NH3 drydep frequency [1/s].  Also accounts for the fraction
            ! of each grid box (I,J,L) that is located beneath the PBL top
            FREQ = DEPSAV(I,J,DRYNH3) * F_UNDER_TOP

            ! Add option for non-local PBL (Lin, 03/31/09)
            IF ( LNLPBL ) FREQ = 0.D0

            ! Only compute drydep loss if FREQ is nonzero
            IF ( FREQ > 0d0 ) THEN

               ! Initial NH3 [v/v]
               NH30 = STT(I,J,L,IDTNH3)
            
               ! Amount of NH3 lost to drydep [v/v]
               NH3 = NH30 * ( 1d0 - EXP( -FREQ * DTCHEM ) )

               ! Prevent underflow condition
               IF ( NH3 < SMALLNUM ) NH3 = 0d0

               ! Subtract NH3 lost to drydep from initial NH3 [v/v]
               STT(I,J,L,IDTNH3) = NH30 - NH3

               !========================================================
               ! Archive drydep flux of NH3 [molec/cm2/s]
               ! This is needed for SOIL_DRYDEP and the ND44 diagnostic
               !========================================================
               IF ( NH3 > 0d0 ) THEN

                  ! Surface area [cm2]
                  AREA_CM2 = GET_AREA_CM2( I, J, L )
                  
                  ! Convert drydep loss from [v/v/timestep] to [molec/cm2/s]
                  FLUX = NH3  * State_Met%AD(I,J,L) / TCVV(IDTNH3)
                  FLUX = FLUX * XNUMOL(IDTNH3) / AREA_CM2 / DTCHEM

                  ! Store dryd flx in DRYD_TMP for use below                  
                  DRYD_FLX(I,J,L) = DRYD_FLX(I,J,L) + FLUX
               ENDIF
            ENDIF
         ENDIF
      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      !===============================================================
      ! If we are not using the non-local PBL scheme, then:
      !===============================================================
      IF ( .not. LNLPBL ) THEN 

         !------------------------------------------------------------
         ! (1) Call SOIL_DRYDEP to archive the NH3 drydep fluxes.
         !------------------------------------------------------------
         IF ( LSOILNOX ) THEN
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, IDTNH3 )
            DO L = 1, LLPAR
            DO J = 1, JJPAR
            DO I = 1, IIPAR
               CALL SOIL_DRYDEP( I, J, L, IDTNH3, DRYD_FLX(I,J,L) )
            ENDDO
            ENDDO
            ENDDO
!$OMP END PARALLEL DO
         ENDIF

         !------------------------------------------------------------
         ! (2) Save archived NH3 drydep fluxes into the AD44 array
         !------------------------------------------------------------
         IF ( ND44 > 0 ) THEN 
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
            DO J = 1, JJPAR
            DO I = 1, IIPAR
            DO L = 1, LLPAR
               AD44(I,J,DRYNH3,1) = AD44(I,J,DRYNH3,1) + DRYD_FLX(I,J,L)
            ENDDO
            ENDDO
            ENDDO
!$OMP END PARALLEL DO
         ENDIF
      ENDIF

      ! Free pointer
      NULLIFY( STT )

      END SUBROUTINE CHEM_NH3
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: chem_nh4
!
! !DESCRIPTION: Subroutine CHEM\_NH4 removes NH4 from the surface via dry
!  deposition. (rjp, bdf, bmy, 1/2/02, 10/25/05) 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CHEM_NH4( am_I_Root, Input_Opt, 
     &                     State_Met, State_Chm, RC )
!
! !USES:
!
      USE CMN_DIAG_MOD
      USE CMN_SIZE_MOD
      USE DIAG_MOD,           ONLY : AD44
      USE DRYDEP_MOD,         ONLY : DEPSAV
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE GIGC_State_Chm_Mod, ONLY : ChmState
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE GRID_MOD,           ONLY : GET_AREA_CM2
      USE PBL_MIX_MOD,        ONLY : GET_FRAC_UNDER_PBLTOP
      USE PBL_MIX_MOD,        ONLY : GET_PBL_MAX_L
      USE TIME_MOD,           ONLY : GET_TS_CHEM
      USE TRACERID_MOD,       ONLY : IDTNH4
      USE GET_NDEP_MOD,       ONLY : SOIL_DRYDEP
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
!  If you are using the non-local PBL mixing scheme (VDIFF), then routine
!  SOIL_DRYDEP and the ND44 diagnostic updating will be done there.  
!                                                                             .
!  If you are using the full TURBDAY PBL mixing scheme, then we have to 
!  call the SOIL_DRYDEP and archive the ND44 diagnostics here.
!                                                                             .
!  Reaction List:
!  ============================================================================
!  (1 ) NH4 = NH4_0 * EXP( -dt )  where d = dry deposition rate [s-1]
!        
! !REVISION HISTORY:
!  (1 ) Now reference AD from "dao_mod.f".  Added parallel DO-loops.  
!        Updated comments, cosmetic changes. (rjp, bmy, bdf, 9/16/02)
!  (2 ) Now replace DXYP(JREF)*1d4 with routine GET_AREA_CM2 of "grid_mod.f".
!        Now use function GET_TS_CHEM from "time_mod.f" (bmy, 3/27/03)
!  (3 ) Now reference PBLFRAC from "drydep_mod.f".  Now apply dry deposition
!        to the entire PBL.  Added L and FREQ variables.  Recode to avoid 
!        underflow from EXP(). (rjp, bmy, 8/1/03) 
!  (4 ) Now use ND44_TMP array to store vertical levels of drydep flux, then 
!        sum into AD44 array.  This preents numerical differences when using
!        multiple processors. (bmy, 3/24/04)    
!  (5 ) Now use parallel DO-loop to zero ND44_TMP (bmy, 4/14/04)
!  (6 ) Now reference STT & TCVV from "tracer_mod.f".   Also remove reference 
!        to CMN, it's not needed (bmy, 7/20/04)
!  (7 ) Replace PBLFRAC from "drydep_mod.f" with GET_FRAC_UNDER_PBLTOP from 
!        "pbl_mix_mod.f".  Also reference GET_PBL_MAX_L from "pbl_mix_mod.f"
!        Vertical DO-loops can run up to PBL_MAX and not LLTROP. (bmy, 2/22/05)
!  (8 ) Now references XNUMOL from "tracer_mod.f" (bmy, 10/25/05)
!  22 Dec 2011 - M. Payer    - Added ProTeX headers
!  01 Mar 2012 - R. Yantosca - Now use GET_AREA_CM2(I,J,L) from grid_mod.F90
!  31 Jul 2012 - R. Yantosca - Now loop from 1..LLPAR for GIGC compatibility
!  31 Jul 2012 - R. Yantosca - Declare temp drydep arrays w/ LLPAR (not LLTROP)
!  14 Nov 2012 - R. Yantosca - Add am_I_Root, Input_Opt, RC as arguments
!  15 Nov 2012 - M. Payer    - Replaced all met field arrays with State_Met
!                              derived type object
!  05 Mar 2013 - R. Yantosca - Now use Input_Opt%LNLPBL
!  13 Mar 2013 - R. Yantosca - Bug fix: make sure we pass values to the
!                              SOIL_DRYDEP routine even when ND44 is off
!                              (this happens when LNLPBL = F)
!  19 Mar 2013 - R. Yantosca - Now copy Input_Opt%TCVV(1:N_TRACERS) and
!                              Input_Opt%XNUMOL(1:N_TRACERS) -- avoid OOB errs
!  25 Mar 2013 - M. Payer    - Now pass State_Chm object via the arg list
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! Scalars
      LOGICAL :: LNLPBL, LSOILNOX
      INTEGER :: I,      J,    L,        N_TRACERS
      REAL*8  :: DTCHEM, NH4,  NH40
      REAL*8  :: FREQ,   FLUX, AREA_CM2, F_UNDER_TOP

      ! Arrays
      REAL*8  :: DRYD_FLX(IIPAR,JJPAR,LLPAR)
      REAL*8  :: TCVV  (Input_Opt%N_TRACERS)
      REAL*8  :: XNUMOL(Input_Opt%N_TRACERS)

      ! Pointers
      ! We need to define local arrays to hold corresponding values 
      ! from the Chemistry State (State_Chm) object. (mpayer, 12/6/12)
      REAL*8, POINTER :: STT(:,:,:,:)

      !=================================================================
      ! CHEM_NH4 begins here!
      !=================================================================
      IF ( IDTNH4 == 0 .or. DRYNH4 == 0 ) RETURN

      ! Assume success
      RC        = GIGC_SUCCESS

      ! Copy fields from INPUT_OPT to local variables for use below
      LNLPBL    = Input_Opt%LNLPBL
      LSOILNOX  = Input_Opt%LSOILNOX
      N_TRACERS = Input_Opt%N_TRACERS
      TCVV      = Input_Opt%TCVV  (1:N_TRACERS)
      XNUMOL    = Input_Opt%XNUMOL(1:N_TRACERS)

      ! Initialize GEOS-Chem tracer array [kg] from Chemistry State object
      ! (mpayer, 12/6/12)
      STT => State_Chm%Tracers

      ! DTCHEM is the chemistry interval in seconds
      DTCHEM    = GET_TS_CHEM() * 60d0 

      ! Zero DRYD_FLX array.  This will store drydep values even if 
      ! the ND44 diagnostic has been turned off (bmy, 3/13/13)
      DRYD_FLX = 0d0

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, F_UNDER_TOP, FREQ, NH40, NH4, AREA_CM2, FLUX )
!$OMP+SCHEDULE( DYNAMIC )
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Fraction of box (I,J,L) underneath the PBL top [unitless]
         F_UNDER_TOP = GET_FRAC_UNDER_PBLTOP( I, J, L )       

         ! Only apply drydep to boxes w/in the PBL
         IF ( F_UNDER_TOP > 0d0 ) THEN

            ! NH4 drydep frequency [1/s].  Also accounts for the fraction
            ! of each grid box (I,J,L) that is located beneath the PBL top
            FREQ = DEPSAV(I,J,DRYNH4) * F_UNDER_TOP

            ! Add option for non-local PBL (Lin, 03/31/09)
            IF ( LNLPBL ) FREQ = 0.D0

            ! Only apply drydep loss if FREQ is nonzero
            IF ( FREQ > 0d0 ) THEN

               ! Initial NH4 [v/v]
               NH40 = STT(I,J,L,IDTNH4)
         
               ! Amount of NH4 lost to drydep [v/v]
               NH4 = NH40 * ( 1d0 - EXP( -FREQ * DTCHEM ) )

               ! Prevent underflow condition
               IF ( NH4 < SMALLNUM ) NH4 = 0d0

               ! Subtract NH4 lost to drydep from initial NH4 [v/v]
               STT(I,J,L,IDTNH4) = NH40 - NH4

               !========================================================
               ! Archive drydep flux of NH4 [molec/cm2/s]
               !========================================================
               IF ( NH4 > 0d0 ) THEN
         
                  ! Surface area [cm2]
                  AREA_CM2 = GET_AREA_CM2( I, J, L )
                  
                  ! Convert drydep loss from [v/v/timestep] to [molec/cm2/s]
                  FLUX = NH4  * State_Met%AD(I,J,L) / TCVV(IDTNH4)
                  FLUX = FLUX * XNUMOL(IDTNH4) / AREA_CM2 / DTCHEM

                  ! Store dryd flx in DRYD_FLX as a placeholder
                  DRYD_FLX(I,J,L) = DRYD_FLX(I,J,L) + FLUX
               ENDIF
            ENDIF
         ENDIF
      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      !===============================================================
      ! If we are not using the non-local PBL scheme, then:
      !===============================================================
      IF ( .not. LNLPBL ) THEN 

         !------------------------------------------------------------
         ! (1) Call SOIL_DRYDEP to archive the NH4 drydep fluxes.
         !------------------------------------------------------------
         IF ( LSOILNOX ) THEN
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, IDTNH4 )
            DO L = 1, LLPAR
            DO J = 1, JJPAR
            DO I = 1, IIPAR
               CALL SOIL_DRYDEP( I, J, L, IDTNH4, DRYD_FLX(I,J,L) )
            ENDDO
            ENDDO
            ENDDO
!$OMP END PARALLEL DO
         ENDIF

         !------------------------------------------------------------
         ! (2) Save archived NH4 drydep fluxes into the AD44 array
         !------------------------------------------------------------
         IF ( ND44 > 0 ) THEN 
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
            DO J = 1, JJPAR
            DO I = 1, IIPAR
            DO L = 1, LLPAR
               AD44(I,J,DRYNH4,1) = 
     &         AD44(I,J,DRYNH4,1) + DRYD_FLX(I,J,L)
            ENDDO
            ENDDO
            ENDDO
!$OMP END PARALLEL DO
         ENDIF
      ENDIF

      ! Free pointer
      NULLIFY( STT )

      END SUBROUTINE CHEM_NH4
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: chem_nit
!
! !DESCRIPTION: Subroutine CHEM\_NIT removes SULFUR NITRATES (NIT) from the
!  surface via dry deposition. (rjp, bdf, bmy, 1/2/02, 5/23/06)  
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CHEM_NIT( am_I_Root, Input_Opt, 
     &                     State_Met, State_Chm, RC )
!
! !USES:
!
      USE CMN_DIAG_MOD
      USE CMN_SIZE_MOD
      USE DIAG_MOD,           ONLY : AD44
      USE DRYDEP_MOD,         ONLY : DEPSAV
      USE GET_NDEP_MOD,       ONLY : SOIL_DRYDEP
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE GIGC_State_Chm_Mod, ONLY : ChmState
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE GRID_MOD,           ONLY : GET_AREA_CM2
      USE PBL_MIX_MOD,        ONLY : GET_FRAC_UNDER_PBLTOP
      USE PBL_MIX_MOD,        ONLY : GET_PBL_MAX_L
      USE TIME_MOD,           ONLY : GET_TS_CHEM
      USE TRACERID_MOD,       ONLY : IDTNIT
      USE TRACERID_MOD,       ONLY : IDTNITs
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
!#endif
!  
!  If you are using the non-local PBL mixing scheme (VDIFF), then routine
!  SOIL_DRYDEP and the ND44 diagnostic updating will be done there.  
!                                                                             .
!  If you are using the full TURBDAY PBL mixing scheme, then we have to 
!  call the SOIL_DRYDEP and archive the ND44 diagnostics here.
!                                                                             .
! !REMARKS:
!  Reaction List:
!  ============================================================================
!  (1 ) NIT = NIT_0 * EXP( -dt )  where d = dry deposition rate [s-1]
!        
! !REVISION HISTORY:
!  (1 ) Now reference AD from "dao_mod.f".  Added parallel DO-loops.  
!        Updated comments, cosmetic changes. (rjp, bmy, bdf, 9/20/02)
!  (2 ) Now replace DXYP(J+J0)*1d4 with routine GET_AREA_CM2 from "grid_mod.f".
!        Now use function GET_TS_CHEM from "time_mod.f" (bmy, 3/27/03)
!  (3 ) Now reference PBLFRAC from "drydep_mod.f".  Now apply dry deposition
!        to the entire PBL.  Added L and FREQ variables.  Recode to avoid
!        underflow from EXP(). (rjp, bmy, 8/1/03) 
!  (4 ) Now use ND44_TMP array to store vertical levels of drydep flux, then 
!        sum into AD44 array.  This preents numerical differences when using
!        multiple processors. (bmy, 3/24/04)    
!  (5 ) Now use parallel DO-loop to zero ND44_TMP (bmy, 4/14/04)
!  (6 ) Now reference STT & TCVV from "tracer_mod.f".  Also remove reference
!        to CMN, it's not needed anymore. (bmy, 7/20/04)
!  (7 ) Replace PBLFRAC from "drydep_mod.f" with GET_FRAC_UNDER_PBLTOP from 
!        "pbl_mix_mod.f".  Also reference GET_PBL_MAX_L from "pbl_mix_mod.f"
!        Vertical DO-loops can run up to PBL_MAX and not LLTROP. (bmy, 2/22/05)
!  (8 ) Now references XNUMOL from "tracer_mod.f" (bmy, 10/25/05)
!  (9 ) Rearrange error check to avoid SEG FAULTS (bmy, 5/23/06)
!  22 Dec 2011 - M. Payer    - Added ProTeX headers
!  01 Mar 2012 - R. Yantosca - Now use GET_AREA_CM2(I,J,L) from grid_mod.F90
!  31 Jul 2012 - R. Yantosca - Now loop from 1..LLPAR for GIGC compatibility
!  31 Jul 2012 - R. Yantosca - Declare temp drydep arrays w/ LLPAR (not LLTROP)
!  09 Nov 2012 - M. Payer    - Replaced all met field arrays with State_Met
!                              derived type object
!  05 Mar 2013 - R. Yantosca - Now use Input_Opt%LNLPBL
!  13 Mar 2013 - R. Yantosca - Bug fix: make sure we pass values to the
!                              SOIL_DRYDEP routine even when ND44 is off
!                              (this happens when LNLPBL = F)
!  19 Mar 2013 - R. Yantosca - Now copy Input_Opt%TCVV(1:N_TRACERS) and
!                              Input_Opt%XNUMOL(1:N_TRACERS) -- avoid OOB errs
!  25 Mar 2013 - M. Payer    - Now pass State_Chm object via the arg list
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! Scalars
      LOGICAL :: LNLPBL,   LSOILNOX,     LSSALT
      INTEGER :: I,        J,      L,    N,     N_TRACERS
      REAL*8  :: DTCHEM,   NIT,    NITs, NIT0,  NIT0s,  E_RKT
      REAL*8  :: E_RKTs,   FLUX,   FREQ, FREQs, RKT,    RKTs   
      REAL*8  :: AREA_CM2, F_UNDER_TOP

      ! Arrays
      REAL*8  :: DRYD_FLX(IIPAR,JJPAR,LLPAR,2)
      REAL*8  :: TCVV  (Input_Opt%N_TRACERS)
      REAL*8  :: XNUMOL(Input_Opt%N_TRACERS)

      ! Pointers
      ! We need to define local arrays to hold corresponding values 
      ! from the Chemistry State (State_Chm) object. (mpayer, 12/6/12)
      REAL*8, POINTER :: STT(:,:,:,:)

      !=================================================================
      ! CHEM_NIT begins here!
      !=================================================================

      ! Return if tracers are not defined
      IF ( IDTNIT == 0 .or. IDTNITs == 0 ) RETURN
      IF ( DRYNIT == 0 .or. DRYNITs == 0 ) RETURN

      ! Assume success
      RC        = GIGC_SUCCESS

      ! Copy fields from INPUT_OPT to local variables for use below
      LNLPBL    = Input_Opt%LNLPBL
      LSOILNOX  = Input_Opt%LSOILNOX
      LSSALT    = Input_Opt%LSSALT
      N_TRACERS = Input_Opt%N_TRACERS
      TCVV      = Input_Opt%TCVV  (1:N_TRACERS)
      XNUMOL    = Input_Opt%XNUMOL(1:N_TRACERS)

      ! Initialize GEOS-Chem tracer array [kg] from Chemistry State object
      ! (mpayer, 12/6/12)
      STT => State_Chm%Tracers

      ! DTCHEM is the chemistry interval in seconds
      DTCHEM    = GET_TS_CHEM() * 60d0 

      ! Zero DRYD_FLX array.  This will store drydep values even if 
      ! the ND44 diagnostic has been turned off (bmy, 3/13/13)
      DRYD_FLX  = 0d0

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I,    J,      L,        NIT0,        NIT0s, NIT   )
!$OMP+PRIVATE( NITs, FREQ,   FREQs,    F_UNDER_TOP, RKT,   E_RKT ) 
!$OMP+PRIVATE( RKTs, E_RKTs, AREA_CM2, FLUX                      )
!$OMP+SCHEDULE( DYNAMIC )
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Initial NITRATE [v/v]
         NIT0  = STT(I,J,L,IDTNIT)

         ! Initial NITRATE w/in seasalt [v/v]
         NIT0s = STT(I,J,L,IDTNITs)

         ! Initialize variables
         NIT   = 0d0
         NITs  = 0d0
         FREQ  = 0d0
         FREQs = 0d0

         ! Fraction of box (I,J,L) underneath the PBL top [unitless]
         F_UNDER_TOP = GET_FRAC_UNDER_PBLTOP( I, J, L )     

         ! Only apply drydep to boxes w/in the PBL
         IF ( F_UNDER_TOP > 0d0 ) THEN 

            !===========================================================
            ! NIT chemistry
            !===========================================================

            ! NIT drydep frequency [1/s].  Also accounts for the fraction
            ! of each vertical level that is located below the PBL top
            FREQ  = DEPSAV(I,J,DRYNIT)  * F_UNDER_TOP

            ! Add option for non-local PBL (Lin, 03/31/09)
            IF ( LNLPBL ) FREQ = 0.D0

            ! If there is drydep ...
            IF ( FREQ > 0d0 ) THEN

               ! Fraction of NIT lost to drydep [unitless] (bec, 12/15/04)
               RKT  = FREQ  * DTCHEM

               ! Pre-compute the exponential term
               E_RKT = EXP( -RKT )

               ! Amount of NITRATE lost to drydep [v/v]
               NIT = NIT0 * ( 1d0 - E_RKT )

               ! Prevent underflow condition
               IF ( NIT < SMALLNUM ) NIT = 0d0

               ! Subtract NITRATE lost to drydep from initial NITRATE [v/v]
               STT(I,J,L,IDTNIT) = NIT0 - NIT

            ELSE
	
               ! No deposition occurs
               NIT = 0d0

            ENDIF

            !===========================================================
            ! NITs chemistry
            !===========================================================

            ! NITs drydep frequency [1/s].  Also accounts for the fraction
            ! of each vertical level that is located below the PBL top
            FREQs = DEPSAV(I,J,DRYNITs) * F_UNDER_TOP
               
            ! Add option for non-local PBL (Lin, 03/31/09)
            IF ( LNLPBL ) FREQs = 0.D0

            ! If there is drydep ...
            IF ( FREQs > 0d0 ) THEN

               ! Fraction of NIT lost to drydep [unitless] (bec, 12/15/04)
               RKTs   = FREQs * DTCHEM

               ! Pre-compute the exponential term
               E_RKTs = EXP( -RKTs )

               ! Compute new NIT concentration [v/v],
               ! updated for seasalt chemistry
               NITs   = ( NIT0s             *          E_RKTs ) + 
     &                  ( PNITs(I,J,L)/RKTs * ( 1.d0 - E_RKTs ) )

            ELSE

               ! NIT prod from HNO3 uptake on fine sea-salt [v/v/timestep]
               NITs = NIT0s + PNITs(I,J,L)

            ENDIF
            
            ! Store final concentration in STT [v/v]
            STT(I,J,L,IDTNITs) = NITs
            
            !========================================================
            ! Archive drydep flux of NIT and NITs [molec/cm2/s]
            !========================================================
         
            ! Surface area [cm2]
            AREA_CM2 = GET_AREA_CM2( I, J, L )
            
            !-------------------
            ! NIT drydep flux
            !-------------------
             
            ! If NIT drydep occurs ...
            IF ( FREQ > 0d0 ) THEN

               ! Convert from [v/v/timestep] to [molec/cm2/s]
               FLUX = NIT * State_Met%AD(I,J,L) / TCVV(IDTNIT)
               FLUX = FLUX * XNUMOL(IDTNIT) / AREA_CM2 / DTCHEM
 
               ! Store NITs drydep flx [molec/cm2/s] for use below
               DRYD_FLX(I,J,L,1) = DRYD_FLX(I,J,L,1) + FLUX

            ENDIF
                  
            !-------------------
            ! NITs drydep flux
            !-------------------
            
            ! NOTE: if drydep doesn't occur then we still have
            ! production from seasalt (bec, bmy, 4/13/05)

            ! Convert from [v/v/timestep] to [molec/cm2/s]
            FLUX = NIT0s - NITs + PNITs(I,J,L) 
            FLUX = FLUX * State_Met%AD(I,J,L) / TCVV(IDTNITs)
            FLUX = FLUX * XNUMOL(IDTNITs) / AREA_CM2 / DTCHEM

            ! Store NITs drydep flx [molec/cm2/s] for use below
            DRYD_FLX(I,J,L,2) = DRYD_FLX(I,J,L,2) + FLUX

         ENDIF
      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      !===============================================================
      ! If we are not using the non-local PBL scheme, then:
      !===============================================================
      IF ( .not. LNLPBL ) THEN 

         !------------------------------------------------------------
         ! (1) Call SOIL_DRYDEP to archive NIT & NITs drydep fluxes
         !------------------------------------------------------------
         IF ( LSOILNOX ) THEN
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, IDTNIT, IDTNITs )
            DO L = 1, LLPAR
            DO J = 1, JJPAR
            DO I = 1, IIPAR
               CALL SOIL_DRYDEP( I, J, L, IDTNIT,  DRYD_FLX(I,J,L,1) )
               CALL SOIL_DRYDEP( I, J, L, IDTNITs, DRYD_FLX(I,J,L,2) )
            ENDDO
            ENDDO
            ENDDO
!$OMP END PARALLEL DO
         ENDIF

         !------------------------------------------------------------
         ! (2) Save archived NIT & NITs drydep fluxes into AD44 array
         !------------------------------------------------------------
         IF ( ND44 > 0 ) THEN 
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
            DO J = 1, JJPAR
            DO I = 1, IIPAR
            DO L = 1, LLPAR
               AD44(I,J,DRYNIT, 1) = 
     &         AD44(I,J,DRYNIT, 1) + DRYD_FLX(I,J,L,1)

               AD44(I,J,DRYNITs,1) = 
     &         AD44(I,J,DRYNITs,1) + DRYD_FLX(I,J,L,2)
            ENDDO
            ENDDO
            ENDDO
!$OMP END PARALLEL DO
         ENDIF         
      ENDIF

      ! Free pointer
      NULLIFY( STT )

      END SUBROUTINE CHEM_NIT
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: emisssulfate
!
! !DESCRIPTION: Subroutine EMISSSULFATE is the interface between the GEOS-CHEM
!  model and the sulfate emissions routines in "sulfate\_mod.f" (bmy, 6/7/00,
!  10/15/09)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE EMISSSULFATE( am_I_Root, Input_Opt, 
     &                         State_Met, State_Chm, RC ) 
!
! !USES:
!
      USE CMN_SIZE_MOD
      USE ERROR_MOD,          ONLY : DEBUG_MSG
      USE GFED2_BIOMASS_MOD,  ONLY : GFED2_IS_NEW
      USE GFED3_BIOMASS_MOD,  ONLY : GFED3_IS_NEW
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE GIGC_State_Chm_Mod, ONLY : ChmState
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE TIME_MOD,           ONLY : GET_SEASON, GET_MONTH
      USE TIME_MOD,           ONLY : GET_YEAR,   ITS_A_NEW_MONTH
      USE TIME_MOD,           ONLY : GET_DAY,    ITS_A_NEW_DAY
      USE TRACERID_MOD,       ONLY : IDTNITs,    IDTSO4s
      USE TRACERID_MOD,       ONLY : IDTDMS,     IDTSO2 
      USE TRACERID_MOD,       ONLY : IDTSO4,     IDTNH3
#if   defined( TOMAS )
      USE TRACERID_MOD,       ONLY : IDTSF1,     IDTAW1,  IDTNH4
      USE TRACERID_MOD,       ONLY : IDTNK1
      USE TOMAS_MOD,          ONLY : IBINS,      ICOMP,   IDIAG
      USE TOMAS_MOD,          ONLY : NH4BULKTOBIN
      USE TOMAS_MOD,          ONLY : SRTNH4
#endif
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)    :: am_I_Root   ! Are we on the root CPU?
      TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
!
! !INPUT/OUTPUT PARAMETERS:
!
      TYPE(MetState), INTENT(INOUT) :: State_Met   ! Meteorology State object
      TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
! 
! !REVISION HISTORY:
!  (1 ) BXHEIGHT is now dimensioned IIPAR,JJPAR,LLPAR (bmy, 9/26/01)
!  (2 ) Removed obsolete commented out code from 9/01 (bmy, 10/24/01)
!  (3 ) Now reference all arguments except FIRSTEMISS, LENV, LEEV from 
!        header files or F90 modules.  Removed NSRCE,  MONTH, JDAY, 
!        LWI, BXHEIGHT, DXYP, AD, PTOP, SIGE, PS, PBL, XTRA2, STT, DATA_DIR, 
!        JYEAR from the arg list.  Now reference GET_PEDGE from F90 module
!        "pressure_mod.f" to compute grid box edge pressures.  Now uses
!        GET_SEASON from "time_mod.f" to get the season.  Now references
!        IDTDMS, IDTSO2, etc from "tracerid_mod.f".  Now make FIRSTEMISS
!        a local SAVEd variable.  Now call READ_BIOMASS_NH3 to read NH3
!        biomass and biofuel emissions. (bmy, 12/13/02)
!  (4 ) Now call READ_NATURAL_NH3 to read the NH3 source from natural
!        emissions. (rjp, bmy, 3/23/03)
!  (5 ) Now use functions GET_SEASON and GET_MONTH from the new "time_mod.f"
!        (bmy, 3/27/03)
!  (6 ) Added first-time printout message (bmy, 4/6/04)
!  (7 ) Now references CMN_SETUP.  Now read ship SO2 if LSHIPSO2=T.  Also
!        references ITS_A_NEW_MONTH from "time_mod.f". (bec, bmy, 5/20/04)
!  (8 ) Now references STT and ITS_AN_AEROSOL_SIM from "tracer_mod.f".  
!        Now references LSHIPSO2 from "logical_mod.f" (bmy, 7/20/04)
!  (9 ) Now references GET_YEAR from "time_mod.f". (bmy, 8/1/05)
!  (10) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (11) Now check if GFED2 has been updated (yc, phs, 12/23/08)
!  (12) Add LANTHRO switch to properly turn off the anthropogenic emissions,
!        READ_AIRCRAFT_SO2, READ_ANTHRO_SOx, READ_ANTHRO_NH3 (ccc, 4/15/09)
!  (13) Now read new volcanic SO2 emissions daily (jaf, bmy, 10/15/09)
!  (14) Add LBIOFUEL switch to properly turn off the biofuel emissions,
!       READ_BIOFUEL_SO2, READ_BIOFUEL_NH3. (ccc, 7/16/10)
!  (14a) Now call SRCSF30 to emit sulfate mass and aerosol number into the 
!        size-resolved TOMAS aerosol tracers. Reference to tracer IDs of 
!        the TOMAS aerosol from tracerid_mod  (win, 1/25/10)
!  22 Dec 2011 - M. Payer    - Added ProTeX headers
!  04 Mar 2013 - R. Yantosca - Remove call to INIT_SULFATE
!  31 Jul 2013 - M. Sulprizio- Now set SO2 aircraft emissions to zero if using
!                              AEIC aircraft emissions (S. Eastham)
!  29 Jan 2014 - R. Yantosca - Avoid array temporaries in call to NH4BULKTOBIN
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL, SAVE :: FIRSTEMISS = .TRUE. 
      INTEGER       :: NSEASON, DAY, MONTH, YEAR   

#if   defined( TOMAS )
      ! Fields for TOMAS simulation
      REAL*8        :: BINMASS(IIPAR,JJPAR,LLPAR,IBINS*ICOMP)
      REAL*8        :: tempnh4(ibins)
      INTEGER       :: TID, I, J, L
      REAL*8        :: MK_TEMP(IBINS)
#endif

      ! For fields from Input_Opt
      LOGICAL       :: LANTHRO
      LOGICAL       :: LBIOFUEL
      LOGICAL       :: LBIOMASS
      LOGICAL       :: LPRT
      LOGICAL       :: LSHIPSO2
      LOGICAL       :: LAEIC
      LOGICAL       :: ITS_AN_AEROSOL_SIM
      INTEGER       :: N_TRACERS

      ! Pointers
      ! We need to define local arrays to hold corresponding values 
      ! from the Chemistry State (State_Chm) object. (mpayer, 12/6/12)
      REAL*8, POINTER :: STT(:,:,:,:)

      !=================================================================
      ! EMISSSULFATE begins here!
      !=================================================================

      ! Assume success
      RC                 = GIGC_SUCCESS

      ! Copy values from Input_Opt
      LANTHRO            = Input_Opt%LANTHRO
      LBIOFUEL           = Input_Opt%LBIOFUEL
      LBIOMASS           = Input_Opt%LBIOMASS
      LPRT               = Input_Opt%LPRT
      LSHIPSO2           = Input_Opt%LSHIPSO2
      LAEIC              = Input_Opt%LAEIC
      ITS_AN_AEROSOL_SIM = Input_Opt%ITS_AN_AEROSOL_SIM
      N_TRACERS          = Input_Opt%N_TRACERS

      ! Initialize GEOS-Chem tracer array [kg] from Chemistry State object
      ! (mpayer, 12/6/12)
      STT => State_Chm%Tracers

      ! Do only on the first timestep
      IF ( FIRSTEMISS ) THEN

         ! Echo info
         WRITE( 6, '(a)' ) REPEAT( '=', 79 )
         WRITE( 6, 100   )
         WRITE( 6, 110   )
         WRITE( 6, 120   )
         WRITE( 6, 130   ) 
         WRITE( 6, '(a)' ) REPEAT( '=', 79 )
       
         ! FORMAT strings
 100     FORMAT( 'S U L F A T E   A E R O S O L   E M I S S I O N S'   )
 110     FORMAT( 'Routines originally by Mian Chin''s GOCART model'    )
 120     FORMAT( 'Modified for GEOS-CHEM by R. Park and R. Yantosca'   ) 
 130     FORMAT( 'Last Modification Date: 4/6/04'                      )

         ! We have now gone thru the first timestep
         FIRSTEMISS  = .FALSE.
      ENDIF

      ! Get the season and month
      NSEASON = GET_SEASON()
      DAY     = GET_DAY()
      MONTH   = GET_MONTH()
      YEAR    = GET_YEAR()

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%% NOTE: This is for the SEAC4RS volcano emissions.  Leave code commented
!%%% out here so that we can restore it later. (bmy, 9/25/13)
!      ! Set year to use for noneruptive volcanic emissions
!      ! (jaf, skim, 1/11/13)
!      VOLCYEAR = YEAR
!      VOLCYEAR = MAX( YEAR, MNYEAR_VOLC )
!      VOLCYEAR = MIN( YEAR, MXYEAR_VOLC )
!
!      ! Only do eruptive emissions if current year data is available
!      ! (jaf, skim, 1/11/13)
!      DO_ERUPTIVE = .FALSE.
!      IF ( ( YEAR .GE. MNYEAR_VOLC ) .AND.
!     &     ( YEAR .LE. MXYEAR_VOLC ) ) DO_ERUPTIVE = .TRUE.
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      !=================================================================
      ! If this is a new day, read in the volcanic emissions 
      ! (jaf, bmy, 10/15/09)
      !=================================================================
      IF ( ITS_A_NEW_DAY() ) THEN 

         ! Read volcanic data
         CALL READ_NONERUP_VOLC( DAY, MONTH, YEAR )
         CALL READ_ERUP_VOLC(    DAY, MONTH, YEAR )

      ENDIF

      !=================================================================
      ! If this is a new month, read in the monthly mean quantities
      !=================================================================
      IF ( ITS_A_NEW_MONTH() ) THEN 

         ! Read monthly mean data
         CALL READ_SST( MONTH, YEAR, State_Met )
         CALL READ_OCEAN_DMS( MONTH )
         CALL READ_NATURAL_NH3( MONTH )

         ! Add LBIOFUEL szitch to turn off biofuel emissions.
         ! Warning: Streets 2006 inventory combines both anthropogenic
         ! and biofuel emissions, so with this inventory it is not
         ! possible to turn off biofuel emissions.
         !(ccc, 7/16/10)
         IF ( LBIOFUEL ) THEN
            CALL READ_BIOFUEL_SO2( MONTH )
#if defined( TOMAS )
            ! Scale anthro SO2 for TOMAS
            ESO2_bf = ESO2_bf * scaleanthso2
#endif
            CALL READ_BIOFUEL_NH3( MONTH )
         ENDIF

         ! Add LANTHRO switch to turn off anthropogenic emissions.
         ! (ccc, 4/15/09)
         IF ( LANTHRO ) THEN

            IF ( LAEIC) THEN
               ESO2_ac = 0d0
            ELSE
               CALL READ_AIRCRAFT_SO2( MONTH, State_Met )
            ENDIF
#if defined( TOMAS )
            ! Scale anthro SO2 for TOMAS
            ESO2_ac = ESO2_ac * scaleanthso2
#endif

            CALL READ_ANTHRO_SOx( Input_Opt, MONTH, NSEASON )
#if defined( TOMAS )
            ! Scale anthro SO2 for TOMAS
            ESO2_an  = ESO2_an * scaleanthso2
            ESO4_an  = ESO4_an * scaleanthso2
#endif

            CALL READ_ANTHRO_NH3( Input_Opt, MONTH )

         ENDIF

         ! Also read ship exhaust SO2 if necessary 
         !CALL READ_SHIP_SO2( MONTH )
         CALL READ_SHIP_SO2( Input_Opt, MONTH )
#if defined( TOMAS )
         ! Scale anthro SO2 for TOMAS
         ESO2_sh  = ESO2_sh * scaleanthso2
#endif

         ! Read oxidants for the offline simulation only
         IF ( ITS_AN_AEROSOL_SIM ) CALL READ_OXIDANT( MONTH )

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%% NOTE: This is for the SEAC4RS volcano emissions.  Leave code commented
!%%% out here so that we can restore it later. (bmy, 9/25/13)
!         ! Read volcanic emission data (jaf, skim, 1/11/13)
!         CALL READ_VOLC( MONTH, VOLCYEAR )
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      ENDIF

      IF(  ( GFED2_IS_NEW() .or. GFED3_IS_NEW() .or.
     &       ITS_A_NEW_MONTH() ) .AND. 
     &     ( LBIOMASS )  ) THEN 
         CALL GET_BIOMASS_SO2
         IF ( LPRT ) CALL DEBUG_MSG( '### EMISSSULFATE: GET_BM_SO2')
         
         CALL GET_BIOMASS_NH3
         IF ( LPRT ) CALL DEBUG_MSG( '### EMISSSULFATE: GET_BM_NH3')         
      ENDIF

      !=================================================================
      ! Add emissions into the STT tracer array
      !=================================================================
      IF ( IDTDMS /= 0 ) THEN
         CALL SRCDMS( am_I_Root, Input_Opt, State_Met, State_Chm, RC )
      ENDIF

      IF ( IDTSO2 /= 0 ) THEN
         CALL SRCSO2( am_I_Root, NSEASON, Input_Opt, State_Chm, RC )
      ENDIF

      IF ( IDTSO4 /= 0 ) THEN 
         CALL SRCSO4( am_I_Root, Input_Opt, State_Chm, RC ) 
      ENDIF

      IF ( IDTNH3 /= 0 ) THEN
         CALL SRCNH3( am_I_Root, Input_Opt, State_Chm, RC ) 
      ENDIF

#if   defined( TOMAS )
      !=================================================================
      ! For TOMAS microphysics
      !
      ! To emit size-resolved sulfate, need to bring the NH4 amount
      ! as one of the aerosol mass component (win, 1/25/10)
      !=================================================================
      IF (IDTSF1 > 0 .and. IDTNK1 > 0 ) THEN

      !!!  Get NH4 and aerosol water into the same array 
         BINMASS(:,:,:,1:IBINS*(ICOMP-IDIAG)) = 
     &        STT(:,:,:,IDTSF1:IDTSF1+IBINS*(ICOMP-IDIAG) - 1)

         IF ( SRTNH4 > 0 ) THEN
            TID = IBINS*(ICOMP-IDIAG) + 1

!%%%OMP PARALLEL DO
!%%%OMP+DEFAULT( SHARED )
!%%%OMP+PRIVATE( I, J, L, TMPNH4, MKTEMP )
!%%%OMP+SCHEDULE( DYNAMIC )
            DO L=1,LLPAR
            DO J=1,JJPAR
            DO I=1,IIPAR

!-----------------------------------------------------------------------------
! Prior to 1/29/14:
! Rewrite subroutine call to avoid array temporaries (bmy, 1/29/14)
!               CALL NH4BULKTOBIN( STT(I,J,L,IDTSF1:IDTSF1-1+IBINS), 
!     &              STT(I,J,L,IDTNH4), TEMPNH4(1:IBINS) )
!-----------------------------------------------------------------------------

               ! Take a slice of size IBINS from STT
               MK_TEMP = STT(I,J,L,IDTSF1:IDTSF1-1+IBINS)

               ! Avoid array tmeporaries in subroutine call (bmy, 1/29/14)
               CALL NH4BULKTOBIN( MK_TEMP, STT(I,J,L,IDTNH4), TEMPNH4 )

               BINMASS(I,J,L,TID:TID+IBINS-1) = TEMPNH4(1:IBINS)
            ENDDO
            ENDDO
            ENDDO
!%%%OMP END PARALLEL DO

         ENDIF

         TID = IBINS*(ICOMP-1) +1
         BINMASS(:,:,:,TID:TID+IBINS-1) = 
     &           STT(:,:,:,IDTAW1:IDTAW1+IBINS-1)


      !IF ( IDTSF1 > 0 ) THEN
         CALL SRCSF30(  STT(:,:,:,IDTNK1:IDTNK1+IBINS-1), 
     &                  BINMASS(:,:,:,:), Input_Opt, State_Met )

         ! Return the aerosol mass after emission subroutine to STT
         ! excluding the NH4 aerosol and aerosol water (win, 9/27/08) 
         STT(:,:,:,IDTSF1:IDTSF1+IBINS*(ICOMP-IDIAG)-1) = 
     &        BINMASS(:,:,:,1:IBINS*(ICOMP-IDIAG))
      ENDIF
#endif

      ! Free pointer
      NULLIFY( STT )

      END SUBROUTINE EMISSSULFATE
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: SULFATE_PBL_MIX
!
! !DESCRIPTION: Subroutine SULFATE\_PBL\_MIX partitions the total
! anthro sulfate emissions thru the entire boundary layer. Emissions
! above the PBL are not used, and left in their level, regardless of
! the mixing scheme. For non-local mixing scheme, all emissions
! within the PBL are put in the first level.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE SULFATE_PBL_MIX ( EMISS, SULFATE, FRAC_OF_PBL,
     $                             PBL_TOP, IS_LOCAL )
!
! !USES:
!                 
      USE ERROR_MOD,    ONLY : ERROR_STOP
      IMPLICIT NONE
!
! !INPUT PARAMETERS:
!             
      INTEGER, INTENT(IN)    :: PBL_TOP ! Top level of boundary layer
      LOGICAL, INTENT(IN)    :: IS_LOCAL ! mixing scheme
      REAL*8,  INTENT(IN)    :: FRAC_OF_PBL(:) ! 
      REAL*8,  INTENT(IN)    :: EMISS(:)  
!     
! !OUTPUT PARAMETERS:
!             
      REAL*8,  INTENT(INOUT) :: SULFATE(:) ! partitioned emissions
!
! !REVISION HISTORY: 
!  27 Oct 2009 - P. Le Sager - initial
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER                :: TOPMIX, TOPEMISS
      REAL*8                 :: TSULFATE

      ! Zero  at all levels 
      SULFATE  = 0.0

      ! Number of emission levels
      TOPEMISS = SIZE( EMISS )
      
      ! Higher level of emiss to be partitionned
      TOPMIX   = MIN( PBL_TOP, TOPEMISS )        
        
      ! Get total emiss SULFATE in PBL
      TSULFATE = SUM( EMISS(1:TOPMIX) ) 

      
      ! Partition if local scheme
      IF ( IS_LOCAL ) THEN

         ! Fraction of total SULFATE in each layer
         SULFATE( 1:PBL_TOP ) = FRAC_OF_PBL( 1:PBL_TOP ) * TSULFATE
         
      ELSE
         
         SULFATE(1) = TSULFATE

      ENDIF

      ! Do not touch emissions above PBL, regardless of mixing scheme
      IF ( TOPEMISS > TOPMIX )
     $     SULFATE( TOPMIX+1 : TOPEMISS ) = EMISS( TOPMIX+1 : TOPEMISS )
        
!#### DEBUG
!      IF ( ABS( SUM( SULFATE(1:TOPMIX) ) - TSULFATE ) > 1.D-5 ) THEN
!         PRINT*, '### ERROR in SULFATE_PBL_MIX!'
!         PRINT*, '### SUM(SULFATE) : ', SUM( SULFATE(1:TOPMIX) )
!         PRINT*, '### TSULFATE     : ', TSULFATE
!         CALL ERROR_STOP( 'Check SULFATE PBL EMISSIONS MIXING',
!     &                    'SULFATE_PBL_MIX (sulfate_mod.f)' )
!      ENDIF
!#### DEBUG

      END SUBROUTINE SULFATE_PBL_MIX
!EOC    
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: srcdms
!
! !DESCRIPTION: Subroutine SRCDMS, from Mian Chin's GOCART model, add DMS
!  emissions to the tracer array.  Modified for use with the GEOS-CHEM model.
!  (bmy, 6/2/00, 8/16/05)
!\\
!\\
! !INTERFACE:
!

      SUBROUTINE SRCDMS( am_I_Root, Input_Opt,
     &                   State_Met, State_Chm, RC )
!
! !USES:
!
      USE DIAG_MOD,           ONLY : AD13_DMS
      USE DAO_MOD,            ONLY : IS_WATER
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE GIGC_State_Chm_Mod, ONLY : ChmState
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE GRID_MOD,           ONLY : GET_AREA_M2
      USE PBL_MIX_MOD,        ONLY : GET_FRAC_OF_PBL, GET_PBL_TOP_L
      USE TIME_MOD,           ONLY : GET_TS_EMIS
      USE TRACERID_MOD,       ONLY : IDTDMS ! (Lin, 03/31/09)

      USE CMN_SIZE_MOD             ! Size parameters
      USE CMN_DIAG_MOD             ! ND13 (for now)
      USE CMN_GCTM_MOD             ! SCALE_HEIGHT
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
!  (1 ) Now reference NSRCE, LWI, DXYP, XTRA2 from either header files or
!        F90 modules.  Now use routines from "pressure_mod.f" to compute
!        grid box surface pressures. (bmy, 9/18/02)
!  (2 ) Now replace DXYP(J) with routine GET_AREA_M2 of "grid_mod.f"
!        Now use routine GET_TS_EMIS from the new "time_mod.f". (bmy, 3/27/03)
!  (3 ) For GEOS-4, convert PBL from [m] to [hPa] w/ the hydrostatic law.
!        Now references SCALE_HEIGHT from "CMN_GCTM".  Added BLTHIK variable
!        for PBL thickness in [hPa]. (bmy, 1/15/04)
!  (4 ) Remove reference to "pressure_mod.f".  Now reference GET_FRAC_OF_PBL
!        and GET_PBL_TOP_L from "pbl_mix_mod.f". (bmy, 2/22/05)
!  (5 ) Switch from Liss & Merlivat to Nightingale formulation for DMS
!        emissions. (swu, bmy, 8/16/05)
!  22 Dec 2011 - M. Payer    - Added ProTeX headers
!  09 Nov 2012 - M. Payer    - Replaced all met field arrays with State_Met
!                              derived type object
!  25 Mar 2013 - R. Yantosca - Now accept am_I_Root, Input_Opt, State_Chm, RC
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
! 
      ! Molecular weight of DMS, kg/mole
      REAL*8,  PARAMETER     :: DMS_MW = 62d0

      ! Ratio of molecular weights: S/DMS
      REAL*8,  PARAMETER     :: S_DMS = 32d0 / 62d0
!
! !LOCAL VARIABLES:
!
      INTEGER                :: I,      J,    L,     NTOP
      REAL*8                 :: DTSRCE, SST,  Sc,    CONC,   W10 
      REAL*8                 :: ScCO2,  AKw,  ERATE, DMSSRC, FEMIS

      ! For fields from Input_Opt
      LOGICAL                :: LNLPBL

      ! Pointers
      ! We need to define local arrays to hold corresponding values 
      ! from the Chemistry State (State_Chm) object. (mpayer, 12/6/12)
      REAL*8,  POINTER       :: TC(:,:,:)
!
! !EXTERNAL FUNCTIONS:   
!
      REAL*8,  EXTERNAL      :: SFCWINDSQR

      !=================================================================
      ! SRCDMS begins here!
      !=================================================================

      ! Assume success
      RC     =  GIGC_SUCCESS

      ! Copy values from Input_Opt
      LNLPBL = Input_Opt%LNLPBL

      ! Chemistry timestep in seconds
      DTSRCE = GET_TS_EMIS() * 60d0

      ! Point to the DMS tracer slice
      TC     => State_Chm%TRACERS(:,:,:,IDTDMS)

      !=================================================================      
      ! Compute DMS emissions = seawater DMS * transfer velocity
      !=================================================================
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I,   J,     SST,    Sc,   CONC, W10,  ScCO2 )
!$OMP+PRIVATE( AKw, ERATE, DMSSRC, NTOP, L,    FEMIS       )
!$OMP+SCHEDULE( DYNAMIC )
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Sea surface temperature in Celsius
         SST = State_Met%SST(I,J) - 273.15d0

         ! Only do the following for water boxes
         IF ( IS_WATER( I, J, State_Met ) ) THEN

            ! Schmidt number for DMS (Saltzman et al., 1993) 
            Sc = 2674.0d0         - 147.12d0*SST + 
     &           3.726d0*(SST**2) - 0.038d0*(SST**3)
 
            !===========================================================
            ! Calculate transfer velocity in cm/hr  (AKw)  
            !                                      
            ! Tans et al. transfer velocity (1990) for CO2 at 
            ! 25oC (Erickson, 1993)
            !                                                                 
            ! Tans et al. assumed AKW=0 when W10<=3. I modified it 
            ! to let DMS emit at low windseeds too.  Chose 3.6m/s as 
            ! the threshold.        
            !  
            ! Schmidt number for CO2: Sc = 600  (20oC, fresh water)          
            !                         Sc = 660  (20oC, seawater)             
            !                         Sc = 428  (25oC, Erickson 93)   
            !===========================================================
            CONC = DMSo(I,J)    
            W10  = SQRT( SFCWINDSQR( I, J,
     &                               State_Met%U10M,
     &                               State_Met%V10M  ) )

            !-----------------------------------------------------------
            ! Tans et al. (1990) 
            !ScCO2 = 428.d0
            !IF (W10 .LE. 3.6) THEN
            !   AKw = 1.0667d0 * W10
            !ELSE
            !   AKw = 6.4d0 * (W10 - 3.d0)
            !ENDIF
            !-----------------------------------------------------------
            ! Wanninkhof (1992) 
            !ScCO2 = 660.d0
            !AKw = 0.31d0 * W10**2
            !-----------------------------------------------------------
            !! Liss and Merlivat (1986) 
            !ScCO2 = 600.d0
            !IF ( W10 <= 3.6d0 ) then
            !   AKw = 0.17d0 * W10
            !   
            !ELSE IF ( W10 <= 13.d0 ) THEN
            !   AKw = 2.85d0 * W10 - 9.65d0
            !   
            !ELSE
            !   AKw = 5.90d0 * W10 - 49.3d0
            !   
            !ENDIF
            !-----------------------------------------------------------
            ! NOTE: Also need to uncomment this section if using
            !       Tans, Wanninkhof, or Liss & Merlivat
            !IF ( W10 <= 3.6d0 ) THEN
            !   AKw = AKw * ( (ScCO2/Sc)**0.667 )
            !ELSE
            !   AKw = AKw * SQRT(ScCO2/Sc)
            !ENDIF
            !-----------------------------------------------------------
            ! Nightingale [2000] (swu, bmy, 8/16/05)
            !
            ! Note that from Nightingale et al [2000a], 
            ! the best fit formulation should be:
            !
            !   AKw = ( 0.222*W10*W10 + 0.333*W10 ) * sqrt( ScCO2/Sc ) 
            !
            ! But from Nightingale et al [2000b], which reported that 
            ! more measurements were incorported, they claimed that 
            ! the following is the best fit:
            !
            ScCO2 = 600.d0            
            AKw   = ( 0.24d0*W10*W10 + 0.061d0*W10 ) * SQRT( ScCO2/Sc )  
            !-----------------------------------------------------------

            !===========================================================
            ! Calculate emission flux in kg/box/timestep   
            !
            ! AKw    is in cm/hr         : AKw/100/3600    is m/sec.    
            ! CONC   is in nM/L (nM/dm3) : CONC*1E-12*1000 is kmole/m3. 
            ! DMS_MW is in g DMS/mol = kg/kmole                          
            ! ERATE  is in kg DMS/m2/timestep   
            ! DMSSRC is in kg DMS/box/timestep  
            !===========================================================
            ERATE = ( AKw  / 100.d0 / 3600.d0 ) * 
     &              ( CONC * 1.d-12 * 1000.d0 ) * DMS_MW * DTSRCE  

            DMSSRC = ERATE * GET_AREA_M2( I, J, 1 )

#if defined( DEVEL )
            State_Chm%TRAC_TEND(I,J,1,IDTDMS) = 
     &      State_Chm%TRAC_TEND(I,J,1,IDTDMS) + DMSSRC/DTSRCE ! kg DMS/grid/s
#endif
            !===========================================================
            ! Add DMS emissions [kg DMS/box] into the tracer array
            !===========================================================

            ! Top layer of the PBL
            NTOP = CEILING( GET_PBL_TOP_L( I, J ) )
            
            ! Loop thru the boundary layer
            DO L = 1, NTOP

               ! Fraction of PBL spanned by grid box (I,J,L) [unitless]
               FEMIS     = GET_FRAC_OF_PBL( I, J, L )

               ! Add option for non-local PBL (Lin, 03/31/09)
               IF ( LNLPBL ) FEMIS = 0.D0

               ! DMS in box (I,J,L) plus emissions [kg]
               TC(I,J,L) = TC(I,J,L) + ( FEMIS * DMSSRC )

            ENDDO

         ELSE                   

            ! If we are not over water, then there is no DMS source
            DMSSRC = 0.d0

         ENDIF                  
 
         ! Save surface emis for non-local PBL mixing (vdiff_mod.f) (units: kg)
         ! (Lin, 06/09/08)
         if ( LNLPBL ) emis_save(I,J,IDTDMS) = DMSSRC

         !==============================================================
         ! ND13 diagnostic:  DMS emissions [kg S/box/timestep]
         !==============================================================
         IF ( ND13 > 0 ) THEN
            AD13_DMS(I,J) = AD13_DMS(I,J) + ( DMSSRC * S_DMS ) ! / DTSRCE
         ENDIF
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! Free pointer
      NULLIFY( TC )

      END SUBROUTINE SRCDMS 
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: srcso2
!
! !DESCRIPTION: Subroutine SRCSO2 (originally from Mian Chin) computes SO2
!  emissons from aircraft, biomass, and anthro sources. (rjp, bdf, bmy, 6/2/00,
!  12/3/09)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE SRCSO2( am_I_Root, NSEASON, Input_Opt, State_Chm, RC )
!
! !USES:
!
      USE BRAVO_MOD,      ONLY     : GET_BRAVO_ANTHRO, GET_BRAVO_MASK
      USE CAC_ANTHRO_MOD, ONLY     : GET_CANADA_MASK,  GET_CAC_ANTHRO
      USE DIAG_MOD,       ONLY     : AD13_SO2_an,      AD13_SO2_ac
      USE DIAG_MOD,       ONLY     : AD13_SO2_bb,      AD13_SO2_nv
      USE DIAG_MOD,       ONLY     : AD13_SO2_ev,      AD13_SO2_bf
      USE DIAG_MOD,       ONLY     : AD13_SO2_sh
      USE EPA_NEI_MOD,    ONLY     : GET_EPA_ANTHRO,   GET_EPA_BIOFUEL
      USE EPA_NEI_MOD,    ONLY     : GET_USA_MASK
      USE ERROR_MOD,      ONLY     : ERROR_STOP,       GEOS_CHEM_STOP
      USE ERROR_MOD,      ONLY     : IS_SAFE_DIV
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE GIGC_State_Chm_Mod, ONLY : ChmState
      USE GRID_MOD,       ONLY     : GET_AREA_CM2
      USE GRID_MOD,       ONLY     : GET_XOFFSET, GET_YOFFSET
      USE NEI2005_ANTHRO_MOD, ONLY : GET_NEI2005_ANTHRO
      USE NEI2005_ANTHRO_MOD, ONLY : NEI05_MASK => USA_MASK
      USE PBL_MIX_MOD,    ONLY     : GET_FRAC_OF_PBL,  GET_PBL_TOP_L
      USE PRESSURE_MOD,   ONLY     : GET_PEDGE
      USE TIME_MOD,       ONLY     : GET_TS_EMIS,      GET_DAY_OF_YEAR 
      USE TIME_MOD,       ONLY     : GET_DAY_OF_WEEK_LT
      USE TRACERID_MOD,   ONLY     : IDTSO2

      USE CMN_SIZE_MOD             ! Size parameters
      USE CMN_DIAG_MOD             ! ND13, LD13 (for now)
      USE CMN_GCTM_MOD             ! SCALE_HEIGHT
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)    :: am_I_Root   ! Are we on the root CPU?
      INTEGER,        INTENT(IN)    :: NSEASON     ! Season #: 1=DJF; 2=MAM; 
                                                   !           3=JJA; 4=SON
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
!  (1 ) Now reference NSRCE, JDAY, PBL, XTRA2, BXHEIGHT from either header
!        files or F90 modules.  Also use routines from "pressure_mod.f" to
!        compute grid box pressures. (bmy, 9/18/02)
!  (2 ) Now use routines GET_TS_EMIS and GET_DAY_OF_YEAR from the new 
!        "time_mod.f" (bmy, 3/27/03)
!  (3 ) For GEOS-4, convert PBL from [m] to [hPa] w/ the hydrostatic law.
!        Now references SCALE_HEIGHT from "CMN_GCTM".  Added BLTHIK variable
!        to hold PBL thickness in [hPa]. (bmy, 1/15/04)
!  (4 ) Now references AD13_SO2_sh array from "diag_mod.f".  Also references
!        LSHIPSO2 from "CMN_SETUP" (bec, bmy, 5/20/04) 
!  (5 ) Now references LSHIPSO2 from "logical_mod.f" (bmy, 7/20/04)
!  (6 ) Now references routines GET_EPA_ANTHRO and GET_USA_MASK from 
!        "epa_nei_mod.f".  Now references GET_AREA_CM2 from "grid_mod.f".  
!        Now references GET_DAY_OF_WEEK from "time_mod.f"  Now references 
!        LNEI99 from "logical_mod.f".  Now can overwrite the anthro SOx 
!        emissions over the continental US if LNEI99=T.  Now references IDTSO2
!        from "tracerid_mod.f. (rch, rjp, bmy, 11/16/04)
!  (7 ) Remove reference to "pressure_mod.f".  Now reference GET_FRAC_OF_PBL 
!        and GET_PBL_TOP_L from "pbl_mix_mod.f".  Removed reference to header
!        file CMN. (bmy, 2/22/05)
!  (8 ) Now references XNUMOL from "tracer_mod.f" (bmy, 10/25/05)
!  (9 ) Now references GET_BRAVO_ANTHRO and GET_BRAVO_MASK from "bravo_mod.f" 
!        for BRAVO Mexican emissions. (rjp, kfb, bmy, 6/26/06)
!  (10) Bug fix: EPA emissions were overwritten by regular ones when both BRAVO
!         and EPA were used. (phs, 10/4/07)
!  (11) Now use CAC Canadian emissions, if necessary (amv, 1/10/08)
!  (12) Bug fix: Always fill the diagnostic array AD13_SO2_sh because it 
!        is allocated anyway (phs, 2/27/09) 
!  (13) Changed processing of volcanic SO2 emissions (jaf, bmy, 10/15/09)
!  (14) Read NEI now (amv, 10/07/2009)
!  (15) Now calls SULFATE_PBL_MIX to do the PBL mixing of
!        emissions (phs, 10/27/09)
!  (16) Rewrite Aerocom SO2 emissions section to avoid errors on SunStudio
!        compiler. Also avoid division by zero. (bmy, 12/3/09)
!  22 Dec 2011 - M. Payer    - Added ProTeX headers
!  01 Mar 2012 - R. Yantosca - Now use GET_AREA_CM2(I,J,L) from grid_mod.F90
!  14 Jun 2013 - R. Yantosca - Now determine weekday/weekend with respect to
!                              the local time at each grid box.  (Formerly,
!                              this had been done w/r/t the GMT time.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
      ! Ratio of molecular weights: S/SO2
      REAL*8,  PARAMETER     :: S_SO2 = 32d0 / 64d0
!
! !LOCAL VARIABLES:
!
      LOGICAL                :: WEEKDAY, IS_LOCAL
      INTEGER                :: I, J, K, L, NTOP
      INTEGER                :: DOW_LT
      REAL*8                 :: SO2(LLPAR)
      REAL*8                 :: DTSRCE,      SO2SRC
      REAL*8                 :: TSO2,        FEMIS(LLPAR)
      REAL*8                 :: AREA_CM2,    AN,        BF
      REAL*8                 :: SO2an(IIPAR,JJPAR,NOXLEVELS)
      REAL*8                 :: SO2bf(IIPAR,JJPAR)
      REAL*8                 :: TEMPEMISS(NOXLEVELS)
      REAL*8                 :: T_DIFF, T_SUM

      ! New variables for volcanic emissions (jaf, 10/15/09)
      INTEGER                :: IREF, JREF, I0, J0
      REAL*8                 :: PLOW, PHIGH, PVOLCLOW, PVOLCHIGH, FRAC

      ! For fields from Input_Opt
      LOGICAL                :: LHIST,    LNLPBL, LBRAVO, LNEI99  
      LOGICAL                :: LSHIPSO2, LCAC,   LNEI05
      INTEGER                :: N_TRACERS
      REAL*8                 :: XNUMOL(Input_Opt%N_TRACERS)

      ! Pressure edges of volcanic 0.5-km grid. Edges come from
      ! converting 500 m increments to pressure using the U.S. standard
      ! atmosphere since volcanic emissions are provided as meters above
      ! sea level (jaf, 10/15/09)
      REAL*8                 :: VOLCPRESS(LVOLC+1) = (/
     &  1013.25, 954.61, 898.74, 845.56, 794.95, 746.83, 701.08,
     &   657.64, 616.40, 577.28, 540.20, 505.07, 471.81, 440.35,
     &   410.61, 382.52, 356.00, 330.99, 307.42, 285.24, 264.36 /)

      ! Pointers
      ! We need to define local arrays to hold corresponding values 
      ! from the Chemistry State (State_Chm) object. (mpayer, 12/6/12)
      REAL*8, POINTER        :: TC(:,:,:)

      !=================================================================
      ! SRCSO2 begins here!
      !================================================================

      ! Assume success
      RC        =  GIGC_SUCCESS

      ! Copy values from Input_Opt
      LHIST     = Input_Opt%LHIST
      LNLPBL    = Input_Opt%LNLPBL
      LBRAVO    = Input_Opt%LBRAVO
      LNEI99    = Input_Opt%LNEI99
      LSHIPSO2  = Input_Opt%LSHIPSO2
      LCAC      = Input_Opt%LCAC
      LNEI05    = Input_Opt%LNEI05
      N_TRACERS = Input_Opt%N_TRACERS
      XNUMOL    = Input_Opt%XNUMOL(1:N_TRACERS)

      ! Point to the SO2 tracer slice
      TC        => State_Chm%TRACERS(:,:,:,IDTSO2)

      ! DTSRCE is the emission timestep in seconds
      DTSRCE    = GET_TS_EMIS() * 60d0

      ! is it local PBL mixing scheme?
      IS_LOCAL  = ( .NOT. LNLPBL )
      
      !=================================================================
      ! SO2 emissions from volcanoes [kg SO2/box/s].
      ! 
      ! SO2 emissions are saved separately for eruptive and non-eruptive
      ! volcanoes, but in both cases they have been gridded onto a
      ! 1 degree x 1 degree by 0.5 km grid. Horizontal regridding is
      ! done in READ_NONERUP_VOLC and READ_ERUP_VOLC. Here we do the
      ! vertical regridding, with the scheme based on that in AIREMISS
      ! in aircraft_nox_mod.f. Since both types of volcanoes are on the
      ! same grid and have been similarly pre-processed, we can treat
      ! them simultaneously here.
      !
      ! As before, non-eruptive volcanic emissions are injected into the
      ! model layer corresponding to the altitude of the volcanic
      ! crater. Eruptive emissions are injected into the top third of
      ! the volcanic plume. These steps were performed in the pre-
      ! processing and hence are not done here. 
      !
      ! (jaf, bmy, 10/15/09)
      !=================================================================
      ESO2_ev = 0d0
      ESO2_nv = 0d0

      ! Get nested-grid offsets
      I0 = GET_XOFFSET()
      J0 = GET_YOFFSET()
    
      ! Loop over surface grid boxes
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, K )
!$OMP+PRIVATE( PLOW, PHIGH, PVOLCLOW, PVOLCHIGH, FRAC, T_SUM, T_DIFF )
      DO J = 1, JJPAR
      DO I = 1, IIPAR
    
         !==============================================================
         ! Loop over tropospheric GEOS-CHEM levels 
         !==============================================================
         DO L = 1, LLTROP

            ! PLOW  is the pressure at the bottom edge of sigma level L
            ! PHIGH is the pressure at the top    edge of sigma level L 
            PLOW  = GET_PEDGE(I,J,L)
            PHIGH = GET_PEDGE(I,J,L+1)

            ! Make sure PLOW is not smaller than VOLCPRESS(1)
            IF ( L == 1 .AND. PLOW < VOLCPRESS(1) ) PLOW = VOLCPRESS(1)

            ! Exit loop if GEOS-Chem level is higher than top of grid
            IF ( PLOW .LE. VOLCPRESS(LVOLC+1) ) EXIT

            !===========================================================
            ! Loop over the native 1-km volcanic SO2 grid layers
            !===========================================================
            DO K = 1, LVOLC

               ! PVOLCLOW  is pressure at the bottom of 1-km grid layer K
               ! PVOLCHIGH is pressure at the top    of 1-km grid layer K
               PVOLCLOW = VOLCPRESS(K)
               PVOLCHIGH = VOLCPRESS(K+1)

               ! Compute the fraction of each 1-km layer K that
               ! lies within the given GEOS-CHEM layer L
               FRAC = 0.0d0

               ! GC grid box is completely below volcanic grid box
               IF ( PHIGH >= PVOLCLOW ) THEN
                  GOTO 10

               ! Bottom of GC grid box is in volcanic grid box
               ELSE IF ( PLOW < PVOLCLOW .AND. PLOW > PVOLCHIGH ) THEN
                  ! Only part of GC box is in volcanic box
                  IF ( PHIGH < PVOLCHIGH ) THEN
                     FRAC = (PLOW-PVOLCHIGH) / (PVOLCLOW-PVOLCHIGH)
                  ! All of GC box is in volcanic box
                  ELSE
                     FRAC =  (PLOW-PHIGH) / (PVOLCLOW-PVOLCHIGH)
                  ENDIF

               ! Top of GC grid box is in volcanic box (bottom isn't)
               ELSE IF ( PHIGH <  PVOLCLOW .AND. 
     &                   PHIGH >  PVOLCHIGH .AND.
     &                   PLOW  >= PVOLCLOW) THEN
                  FRAC = (PVOLCLOW-PHIGH) / (PVOLCLOW-PVOLCHIGH)

               ! All of volcanic box in GC box
               ELSE IF ( PHIGH <= PVOLCHIGH .AND. 
     &                   PLOW >= PVOLCLOW ) THEN
                  FRAC = 1.0d0

               ENDIF

               ! Add contribution from this layer of volcanic grid (K)
               ! to non-eruptive SO2 in model grid layer L.
               ESO2_nv(I,J,L) = ESO2_nv(I,J,L) +
     &                          FRAC * ENV(I,J,K)
  
               ! Add contribution from this layer of volcanic grid (K)
               ! to eruptive SO2 in model grid layer L.
               ESO2_ev(I,J,L) = ESO2_ev(I,J,L) +
     &                          FRAC * EEV(I,J,K)

            ENDDO ! K

10          CONTINUE
         ENDDO    ! L

         !===========================================================
         ! Consistency check: columns before and after regridding 
         ! should be equal.  Be careful not to divide by zero!
         !===========================================================

         ! Non-eruptive volcanoes
         T_SUM  = SUM( ENV(I,J,:) )
         T_DIFF = ABS( SUM( ESO2_nv(I,J,:) ) - T_SUM )         

         IF ( IS_SAFE_DIV( T_DIFF, T_SUM ) ) THEN
            IF ( T_DIFF / T_SUM > 1d-5 ) THEN
               PRINT*, 'Non-eruptive volcanic SO2 emissions before ' //
     &                 'and after regridding are not equivalent!'
               CALL FLUSH(6)
               CALL GEOS_CHEM_STOP
            ENDIF
         ENDIF

         ! Eruptive volcanoes
         T_SUM  = SUM( EEV(I,J,:) )
         T_DIFF = ABS( SUM( ESO2_ev(I,J,:) ) - T_SUM )         

         IF ( IS_SAFE_DIV( T_DIFF, T_SUM ) ) THEN
            IF ( T_DIFF / T_SUM > 1d-5 ) THEN
               PRINT*, 'Eruptive volcanic SO2 emissions before ' //
     &                 'and after regridding are not equivalent!'
               CALL FLUSH(6)
               CALL GEOS_CHEM_STOP
            ENDIF
         ENDIF

      ENDDO       ! I
      ENDDO       ! J
!$OMP END PARALLEL DO

      !=================================================================
      ! Overwrite USA    w/ EPA/NEI99 (anthro+biofuel) SO2 emissions 
      ! Overwrite MEXICO w/ BRAVO     (anthro only   ) SO2 emissions
      ! Overwrite CANADA w/ CAC       (anthro only   ) SO2 emissions
      !-----------------------------------------------------------------
      ! Note that we:
      ! Overwrite ASIA w/ STREETS and EUROPE w/ EMEP 
      !  in READ_ANTHRO_SOx.
      !
      ! In both cases, SO4 is a fraction of provided SO2 (except for
      ! EPA). It is done in READ_ANTHRO_SOX in the 1st case, and in
      ! SRCSO4 for inventories dealt with here. EPA is the only one to
      ! provide direct SO4, which is why we deal with it here, even
      ! though it does not have to be like that. Historical.
      !
      ! So, since we have EPA here, we have to deal with BRAVO and
      ! CAC here to deal with their overlaping mask.
      !
      ! (amv, phs, 3/12/08, 8/24/09)
      !=================================================================
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, AREA_CM2, AN, BF, DOW_LT, WEEKDAY )
!$OMP+SCHEDULE( DYNAMIC )
      DO J = 1, JJPAR

         ! Initialize
         AN       = 0d0
         BF       = 0d0

         DO I = 1, IIPAR

            ! Grid box surface area [cm2]
            AREA_CM2 = GET_AREA_CM2( I, J, 1 )

!------- Prior to 3/8/07 (phs)         
            !!-----------------------------------------------------------
            !! If we are using EPA/NEI99 (anthro + biofuel) ...
            !!-----------------------------------------------------------
            !IF ( LNEI99 ) THEN
            !
            !   ! If we are over the USA ...
            !   IF ( GET_USA_MASK( I, J ) > 0d0 ) THEN
            !
            !      ! Determine if we should use weekday or weekend NEI
            !      ! emissions at grid box (I,J,L).  Since NEI is over
            !      ! the US, then weekend is Sat/Sun.
            !      DOW_LT       = GET_DAY_OF_WEEK_LT( I, J, 1 )
            !      WEEKDAY      = ( DOW_LT > 0 .and. DOW_LT < 6 )
            !
            !      ! Read EPA/NEI99 SO2 emissions in [molec/cm2/s]
            !      AN           = GET_EPA_ANTHRO( I, J, IDTSO2, WEEKDAY )
            !      BF           = GET_EPA_BIOFUEL(I, J, IDTSO2, WEEKDAY )
            !   
            !      ! Convert anthro SO2 from [molec/cm2/s] to [kg/box/s] 
            !      ! Place all anthro SO2 into surface layer
            !      SO2an(I,J,1) = AN * AREA_CM2 / XNUMOL(IDTSO2)
            !      SO2an(I,J,2) = 0d0
            !   
            !      ! Convert anthro SO2 from [molec/cm2/s] to [kg/box/s] 
            !      SO2bf(I,J)   = BF * AREA_CM2 / XNUMOL(IDTSO2)
            !
            !   ELSE
            !
            !      ! If we are not over the USA, then just use the regular 
            !      ! emissions from ESO2_an and ESO2_bf (bmy, 11/16/04)
            !      SO2an(I,J,1) = ESO2_an(I,J,1)
            !      SO2an(I,J,2) = ESO2_an(I,J,2)
            !      SO2bf(I,J)   = ESO2_bf(I,J)
            !
            !   ENDIF
            !
            !ELSE
            ! 
            !   ! If we are not using EPA/NEI99 emissions, then just copy 
            !   ! ESO2_an and ESO2_bf into local arrays (bmy, 11/16/04)
            !   SO2an(I,J,1) = ESO2_an(I,J,1)
            !   SO2an(I,J,2) = ESO2_an(I,J,2)
            !   SO2bf(I,J)   = ESO2_bf(I,J)
            !
            !ENDIF
            !
            !!-----------------------------------------------------------
            !! If we are using BRAVO emissions over Mexico ...
            !!-----------------------------------------------------------
            !IF ( LBRAVO ) THEN
            !
            !   ! If we are over Mexico ...
            !   IF ( GET_BRAVO_MASK( I, J ) > 0d0 ) THEN
            !
            !      ! Read BRAVO SO2 emissions in [molec/cm2/s]
            !      AN           = GET_BRAVO_ANTHRO( I, J, IDTSO2 )
            !   
            !      ! Convert anthro SO2 from [molec/cm2/s] to [kg/box/s] 
            !      ! Place all anthro SO2 into surface layer
            !      SO2an(I,J,1) = AN * AREA_CM2 / XNUMOL(IDTSO2)
            !      SO2an(I,J,2) = 0d0
            !
            !   ELSE
            !
            !      ! If we are not over MEXICO, then just use 
            !      ! the regular emissions from ESO2_an
            !      SO2an(I,J,1) = ESO2_an(I,J,1)
            !      SO2an(I,J,2) = ESO2_an(I,J,2)
            !
            !   ENDIF
            !
            !ELSE
            ! 
            !   ! If we are not using BRAVO emissions, then just copy 
            !   ! ESO2_an and ESO2_bf into local arrays
            !   SO2an(I,J,1) = ESO2_an(I,J,1)
            !   SO2an(I,J,2) = ESO2_an(I,J,2)
            !   SO2bf(I,J)   = ESO2_bf(I,J)
            !
            !ENDIF

            !-----------------------------------------------------------
            ! Default SO2 from GEIA or EDGAR (w/ optional STREETS for 
            ! ASIA, and EMEP for Europe)
            !-----------------------------------------------------------
            DO L = 1, NOXLEVELS
               SO2an(I,J,L) = ESO2_an(I,J,L)
            ENDDO
            SO2bf(I,J)   = ESO2_bf(I,J)

            !-----------------------------------------------------------
            ! If we are using EPA/NEI99 over the USA (anthro + biofuel)
            !-----------------------------------------------------------
![eml Don't overwrite US if historical
            IF ( LNEI99 .and. .not. LHIST ) THEN
!eml]
            IF ( GET_USA_MASK( I, J ) > 0d0 ) THEN

               ! Determine if we should use weekday or weekend NEI
               ! emissions at grid box (I,J,L).  Since NEI is over
               ! the US, then weekend is Sat/Sun.
               DOW_LT       = GET_DAY_OF_WEEK_LT( I, J, 1 )
               WEEKDAY      = ( DOW_LT > 0 .and. DOW_LT < 6 )
            
               ! Read EPA/NEI99 SO2 emissions in [molec/cm2/s]
               AN           = GET_EPA_ANTHRO( I, J, IDTSO2, WEEKDAY )
               BF           = GET_EPA_BIOFUEL(I, J, IDTSO2, WEEKDAY )
               
               ! Convert anthro SO2 from [molec/cm2/s] to [kg/box/s] 
               ! Place all anthro SO2 into surface layer
               SO2an(I,J,1) = AN * AREA_CM2 / XNUMOL(IDTSO2)
#if defined( TOMAS )
     &                      * scaleanthso2
#endif

               DO L = 2, NOXLEVELS
                  SO2an(I,J,L) = 0d0
               ENDDO
               
               ! Convert biofuel SO2 from [molec/cm2/s] to [kg/box/s] 
               SO2bf(I,J)   = BF * AREA_CM2 / XNUMOL(IDTSO2)
#if defined( TOMAS )
     &                      * scaleanthso2
#endif

            ENDIF
            ENDIF

            !-----------------------------------------------------------
            ! If we are using BRAVO emissions over Mexico (anthro)
            !-----------------------------------------------------------
![eml
            IF ( LBRAVO .and. .not. LHIST )  THEN 
!eml]
            IF ( GET_BRAVO_MASK( I, J ) > 0d0 )  THEN
            
               ! Read BRAVO SO2 emissions in [molec/cm2/s]
               AN           = GET_BRAVO_ANTHRO( I, J, IDTSO2 )
               
               ! Convert anthro SO2 from [molec/cm2/s] to [kg/box/s] 
               ! Place all anthro SO2 into surface layer.
               ! Add to USA emissions if on the border. 
               IF ( LNEI99 ) then
                  IF ( GET_USA_MASK( I, J) > 0d0 ) THEN

                  SO2an(I,J,1) = SO2an(I,J,1) + 
     &                           AN * AREA_CM2 / XNUMOL(IDTSO2)
#if defined( TOMAS )
     &                           *scaleanthso2
#endif
                 ENDIF 
               ELSE
                  SO2an(I,J,1) = AN * AREA_CM2 / XNUMOL(IDTSO2)
#if defined( TOMAS )
     &                           *scaleanthso2
#endif
               ENDIF

               DO L = 2, NOXLEVELS
                  SO2an(I,J,L) = 0d0
               ENDDO

            ENDIF
            ENDIF


            !-----------------------------------------------------------
            ! If we are using CAC emissions over Canada ...
            !-----------------------------------------------------------
![eml
            IF ( LCAC .and. .not. LHIST ) THEN
!eml]
            IF ( GET_CANADA_MASK( I, J ) > 0d0 ) THEN

               ! Read CAC SO2 emissions in [molec/cm2/s]
               AN           = GET_CAC_ANTHRO( I, J, IDTSO2,
     &                            MOLEC_CM2_s=.TRUE. )

               ! Convert anthro SO2 from [molec/cm2/s] to [kg/box/s] 
               ! Place all anthro SO2 into surface layer.
               ! Add to USA emissions if on the border. 
               IF ( LNEI99 ) then
                  if ( GET_USA_MASK( I, J) > 0d0 ) THEN

                     SO2an(I,J,1) = SO2an(I,J,1) + 
     &                    AN * AREA_CM2 / XNUMOL(IDTSO2)
#if defined( TOMAS )
     &                           *scaleanthso2
#endif
                  endif
               ELSE
                  SO2an(I,J,1) = AN * AREA_CM2 / XNUMOL(IDTSO2)
#if defined( TOMAS )
     &                           *scaleanthso2
#endif
               ENDIF

               DO L = 2, NOXLEVELS
                  SO2an(I,J,L) = 0d0
               ENDDO

            ENDIF
            ENDIF


            !-----------------------------------------------------------
            ! If we are using EPA/NEI 2005 over USA.
            ! Must be called after CAC and BRAVO to simply overwrite
            ! where they overlap
            !-----------------------------------------------------------  
![eml
            IF ( LNEI05 .and. .not. LHIST ) THEN
!eml]
            IF ( NEI05_MASK( I, J ) > 0d0 ) THEN

               ! Determine if we should use weekday or weekend NEI
               ! emissions at grid box (I,J,L).  Since NEI is over
               ! the US, then weekend is Sat/Sun.
               DOW_LT  = GET_DAY_OF_WEEK_LT( I, J, 1 )
               WEEKDAY = ( DOW_LT > 0 .and. DOW_LT < 6 )

               ! Read USA SO2 emissions in [molec/cm2/s]
               DO L = 1, NOXLEVELS

                  AN = GET_NEI2005_ANTHRO( I, J, L, IDTSO2, WEEKDAY,
     &                            MOLEC_CM2_s=.TRUE. )

                  ! Convert anthro SO2 from [molec/cm2/s] to [kg/box/s]
                  SO2an(I,J,L) =  AN * AREA_CM2 / XNUMOL(IDTSO2)
               ENDDO

            ENDIF
            ENDIF

         ENDDO

#if defined( DEVEL )
         ! kg/m2/s
         State_Chm%TRAC_TEND(:,J,1,IDTSO2) = 
     &   State_Chm%TRAC_TEND(:,J,1,IDTSO2) + 
     &       ( SO2an(:,J,1) + Eso2_bb(:,J) +
     &         SO2bf(:,J)   + ESO2_sh(:,J) )

         State_Chm%TRAC_TEND(:,J,2:NOXLEVELS,IDTSO2) = 
     &   State_Chm%TRAC_TEND(:,J,2:NOXLEVELS,IDTSO2) + SO2an(:,J,2:)

#endif
      ENDDO
!$OMP END PARALLEL DO

      !=================================================================
      ! Add SO2 emissions into model levels
      !=================================================================
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, NTOP, L, SO2, FEMIS, SO2SRC, TEMPEMISS )
!$OMP+SCHEDULE( DYNAMIC )
      DO J = 1, JJPAR
      DO I = 1, IIPAR

          ! Top of the boundary layer
          NTOP = CEILING( GET_PBL_TOP_L( I, J ) ) 

          
          ! Fraction of PBL spanned by grid box (I,J,L) [unitless]
          DO L = 1, LLPAR
             FEMIS(L)  = GET_FRAC_OF_PBL( I, J, L )
          ENDDO

         !=============================================================
         ! Do PBL mixing of anthro (multilevel) + biomass + biofuel +
         ! ship exhaust SO2 at (I,J)
         !=============================================================
         TEMPEMISS(1) = SO2an(I,J,1) + Eso2_bb(I,J) +
     $                  SO2bf(I,J)   + ESO2_sh(I,J)

         TEMPEMISS(2:)= SO2an(I,J,2:)
         
         call SULFATE_PBL_MIX( TEMPEMISS, SO2, 
     $                         FEMIS,     NTOP, IS_LOCAL )

         
         !==============================================================
         ! Add anthro SO2, aircraft SO2, volcano SO2, and biomass SO2
         ! Convert from [kg SO2/box/s] -> [kg SO2/box/timestep]
         !==============================================================
         DO L = 1, LLPAR
            
            ! SO2 emissions [kg/box/s]
            SO2SRC = SO2(L)         + ESO2_ac(I,J,L) + 
     &               ESO2_nv(I,J,L) + ESO2_ev(I,J,L) 

            ! Add SO2 to TC array [kg/box/timestep]
            ! Add option for non-local PBL (Lin, 03/31/09)
            IF (.NOT. LNLPBL) THEN
               TC(I,J,L) = TC(I,J,L) + ( SO2SRC * DTSRCE )
            ELSE
               IF (L /= 1) then
                  TC(I,J,L) = TC(I,J,L) + ( SO2SRC * DTSRCE )
               ELSE
                  ! Save surface emis for non-local PBL mixing (vdiff_mod.f) 
                  ! (units: kg) (Lin, 06/09/08)
                  EMIS_SAVE(I,J,IDTSO2) = SO2SRC * DTSRCE
               ENDIF
            ENDIF
         ENDDO

         !==============================================================
         ! ND13 Diagnostic: SO2 emissions in [kg S/box/timestep]
         !==============================================================
         IF ( ND13 > 0 ) THEN 

            ! Anthropogenic SO2 -- Levels 1-NOXLEVELS
            DO L = 1, NOXLEVELS
               AD13_SO2_an(I,J,L) = AD13_SO2_an(I,J,L) +
     &                              ( SO2an(I,J,L) * S_SO2 * DTSRCE )
            ENDDO

            ! SO2 from biomass burning
            AD13_SO2_bb(I,J)      = AD13_SO2_bb(I,J) +
     &                              ( ESO2_bb(I,J) * S_SO2 * DTSRCE )
 
            ! SO2 from biofuel burning
            AD13_SO2_bf(I,J)      = AD13_SO2_bf(I,J) +
     &                              ( SO2bf(I,J)   * S_SO2 * DTSRCE )

            ! SO2 from ship emissions (bec, bmy, 5/20/04)
            ! Always fill the diagnostic array since
            ! it is allocated anyway (phs, 2/27/09) 
            AD13_SO2_sh(I,J)      = AD13_SO2_sh(I,J) +
     &                              ( ESO2_sh(I,J) * S_SO2 * DTSRCE )

            ! Loop thru LD13 levels
            DO L = 1, LD13 

               ! SO2 from aircraft emissions
               AD13_SO2_ac(I,J,L) = AD13_SO2_ac(I,J,L) +
     &                              ( ESO2_ac(I,J,L) * S_SO2 * DTSRCE )

               ! SO2 from non-eruptive volcanoes
               AD13_SO2_nv(I,J,L) = AD13_SO2_nv(I,J,L) +
     &                              ( ESO2_nv(I,J,L) * S_SO2 * DTSRCE )

               ! SO2 from eruptive volcanoes
               AD13_SO2_ev(I,J,L) = AD13_SO2_ev(I,J,L) +
     &                              ( ESO2_ev(I,J,L) * S_SO2 * DTSRCE )
            ENDDO
         ENDIF
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! Free pointer
      NULLIFY( TC )

      END SUBROUTINE SRCSO2
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: srcso4
!
! !DESCRIPTION: Subroutine SRCSO4 (originally from Mian Chin) computes SO4 
!  emissions from anthropogenic sources (rjp, bdf, bmy, 6/2/00, 5/27/09)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE SRCSO4( am_I_Root, Input_Opt, State_Chm, RC )
!
! !USES:
!
      USE BRAVO_MOD,          ONLY : GET_BRAVO_ANTHRO, GET_BRAVO_MASK
      USE CAC_ANTHRO_MOD,     ONLY : GET_CANADA_MASK,  GET_CAC_ANTHRO
      USE NEI2005_ANTHRO_MOD, ONLY : GET_NEI2005_ANTHRO
      USE NEI2005_ANTHRO_MOD, ONLY : NEI05_MASK => USA_MASK
      USE DIAG_MOD,           ONLY : AD13_SO4_an,     AD13_SO4_bf
      USE EPA_NEI_MOD,        ONLY : GET_EPA_ANTHRO,  GET_EPA_BIOFUEL
      USE EPA_NEI_MOD,        ONLY : GET_USA_MASK
      USE ERROR_MOD,          ONLY : ERROR_STOP
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE GIGC_State_Chm_Mod, ONLY : ChmState
      USE GRID_MOD,           ONLY : GET_AREA_CM2
      USE PBL_MIX_MOD,        ONLY : GET_FRAC_OF_PBL, GET_PBL_TOP_L
      USE TIME_MOD,           ONLY : GET_TS_EMIS
      USE TIME_MOD,           ONLY : GET_DAY_OF_WEEK_LT
      USE TRACERID_MOD,       ONLY : IDTSO4, IDTSO2

      USE CMN_SIZE_MOD             ! Size parameters
      USE CMN_DIAG_MOD             ! ND13 (for now)
      USE CMN_GCTM_MOD             ! SCALE_HEIGHT
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
! 
! !REVISION HISTORY:
!  (1 ) Emission of SO4 is read in SULFATE_READYR, in [kg/box/s]. 
!        It is converted to [kg/box/timestep] here.  
!  (2 ) Now use routine GET_TS_EMIS from the new "time_mod.f" (bmy, 3/27/03)
!  (3 ) For GEOS-4, convert PBL from [m] to [hPa] w/ the barometric law.
!        Now references SCALE_HEIGHT from "CMN_GCTM".  Added BLTHIK variable
!        to hold PBL thickness in [hPa]. (bmy, 1/15/04)
!  (4 ) Now references GET_EPA_ANTHRO, GET_EPA_BIOFUEL, and GET_USA_MASK from 
!        "epa_nei_mod.f".  Now references AD13_SO4_bf from "diag_mod.f".  Now 
!        references GET_AREA_CM2 from "grid_mod.f".  Now references 
!        GET_DAY_OF_WEEK from "time_mod.f".  Now references LNEI99 from 
!        "logical_mod.f".  Now can overwrite the anthro SOx emissions over 
!        the continental US if LNEI99=T.  Now references IDTSO4 from 
!        "tracerid_mod.f". (rch, rjp, bmy, 11/16/04)
!  (5 ) Remove reference to "pressure_mod.f".  Now reference GET_FRAC_OF_PBL 
!        and GET_PBL_TOP_L from "pbl_mix_mod.f".  Removed reference to header 
!        file CMN. (bmy, 2/22/05)
!  (6 ) Now references XNUMOL & XNUMOLAIR from "tracer_mod.f" (bmy, 10/25/05)
!  (7 ) Now overwrite CAC emissions over Canada, if necessary (amv,  1/10/08)
!  (8 ) Need to add CAC_AN to the PRIVATE statement (bmy, 5/27/09)
!  (9 ) Now account for BRAVO SO4. Fix typo for CAC (phs, 8/24/09)
!  (10) Now account for NEI 2005 inventory (amv, 10/07/2009)
!  (11) Now calls SULFATE_PBL_MIX to do the PBL mixing of
!        emissions (phs, 10/27/09)
!  22 Dec 2011 - M. Payer    - Added ProTeX headers
!  01 Mar 2012 - R. Yantosca - Now use GET_AREA_CM2(I,J,L) from grid_mod.F90
!  25 Mar 2013 - R. Yantosca - Now accept am_I_Root, Input_Opt, State_Chm, RC
!  14 Jun 2013 - R. Yantosca - Now determine weekday/weekend with respect to
!                              the local time at each grid box.  (Formerly,
!                              this had been done w/r/t the GMT time.
!  13 Dec 2013 - M. Sulprizio- Bug fix: Make sure we apply NEI05 emissions when
!                              NEI05_MASK > 0, not NEI05_MASK < 0 (J. Fisher)
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
      ! Ratio of molecular weights: S/SO4
      REAL*8,  PARAMETER     :: S_SO4 = 32d0 / 96d0
!
! !LOCAL VARIABLES:
!
      LOGICAL                :: WEEKDAY, IS_LOCAL
      INTEGER                :: I, J, K, L, DOW_LT, NTOP
      REAL*8                 :: SO4(LLPAR), DTSRCE  
      REAL*8                 :: TSO4,       FEMIS(LLPAR)
      REAL*8                 :: AREA_CM2,   EPA_AN,  EPA_BF
      REAL*8                 :: AN
      REAL*8                 :: SO4an(IIPAR,JJPAR,NOXLEVELS)
      REAL*8                 :: SO4bf(IIPAR,JJPAR)
      REAL*8                 :: TEMPEMISS(NOXLEVELS)

      ! For fields from Input_Opt
      LOGICAL                :: LHIST,  LNLPBL, LBRAVO  
      LOGICAL                :: LNEI99, LCAC,   LNEI05
      INTEGER                :: N_TRACERS
      REAL*8                 :: XNUMOL(Input_Opt%N_TRACERS)

      ! Pointers
      ! We need to define local arrays to hold corresponding values 
      ! from the Chemistry State (State_Chm) object. (mpayer, 12/6/12)
      REAL*8, POINTER        :: TC(:,:,:)

      !=================================================================
      ! SRCSO4 begins here!
      !=================================================================

      ! Assume success
      RC        =  GIGC_SUCCESS

      ! Copy values from Input_Opt
      LHIST     = Input_Opt%LHIST
      LNLPBL    = Input_Opt%LNLPBL
      LBRAVO    = Input_Opt%LBRAVO
      LNEI99    = Input_Opt%LNEI99
      LCAC      = Input_Opt%LCAC
      LNEI05    = Input_Opt%LNEI05
      N_TRACERS = Input_Opt%N_TRACERS
      XNUMOL    = Input_Opt%XNUMOL(1:N_TRACERS)

      ! Point to the SO4 tracer slice
      TC        => State_Chm%TRACERS(:,:,:,IDTSO4)

      ! DTSRCE is the emission timestep in seconds
      DTSRCE    = GET_TS_EMIS() * 60d0

      ! Get current day of the week
      ! is it local PBL mixing scheme?
      IS_LOCAL  = .NOT. LNLPBL

      !=================================================================
      ! Overwrite USA    w/ EPA/NEI99 SO4 emissions       (if necessary)
      ! Overwrite CANADA w/ CAC       SO2-fraction emiss. (if necessary)
      ! Overwrite MEXICO w/ BRAVO     SO2-fraction emiss. (if necessary)
      ! Store emissions into local arrays SO4an, SO4bf
      !=================================================================
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, AREA_CM2, EPA_AN, EPA_BF, AN )
!$OMP+PRIVATE( NTOP, FEMIS, SO4, TEMPEMISS, DOW_LT, WEEKDAY )
      DO J = 1, JJPAR
         
         DO I = 1, IIPAR

            ! Grid box surface area [cm2]
            AREA_CM2 = GET_AREA_CM2( I, J, 1 )

            ! Top level of boundary layer at (I,J)
            NTOP = CEILING( GET_PBL_TOP_L( I, J ) )

            ! Fraction of PBL spanned by grid box (I,J,L) [unitless]
            DO L = 1, LLPAR
               FEMIS(L)  = GET_FRAC_OF_PBL( I, J, L )
            ENDDO

            !-----------------------------------------------------------
            ! Default SO4 from GEIA, or (as a fraction of SO2) from 
            ! EDGAR w/ optional STREETS for S.E.-ASIA, and optional
            ! EMEP for Europe
            !-----------------------------------------------------------
            DO L = 1, NOXLEVELS
               SO4an(I,J,L) = ESO4_an(I,J,L)
            ENDDO
            SO4bf(I,J)   = 0d0

            
            !-----------------------------------------------------------
            ! If we are using EPA/NEI99 emissions and over the USA
![eml don't overwrite if historical
            IF ( LNEI99 .and. .not. LHIST ) THEN
!eml]
            IF ( GET_USA_MASK( I, J ) > 0d0 )  THEN

                  ! Determine if we should use weekday or weekend NEI
                  ! emissions at grid box (I,J,L).  Since NEI is over
                  ! the US, then weekend is Sat/Sun.
                  DOW_LT       = GET_DAY_OF_WEEK_LT( I, J, 1 )
                  WEEKDAY      = ( DOW_LT > 0 .and. DOW_LT < 6 )

                  ! Read SO4 emissions in [molec/cm2/s]
                  EPA_AN       = GET_EPA_ANTHRO( I, J, IDTSO4, WEEKDAY )
                  EPA_BF       = GET_EPA_BIOFUEL(I, J, IDTSO4, WEEKDAY ) 

                  ! Convert anthro SO4 from [molec/cm2/s] to [kg/box/s] 
                  ! Place all EPA/NEI99 anthro SO4 into surface layer
                  SO4an(I,J,1) = EPA_AN * AREA_CM2 / XNUMOL(IDTSO4)
#if defined( TOMAS )
     &                           *scaleanthso2
#endif
                  DO L = 2,NOXLEVELS
                     SO4an(I,J,L) = 0d0
                  ENDDO
               
                  ! Convert biofuel SO4 from [molec/cm2/s] to [kg/box/s]
                  SO4bf(I,J)   = EPA_BF * AREA_CM2 / XNUMOL(IDTSO4)
#if defined( TOMAS )
     &                           *scaleanthso2
#endif
            ENDIF
            ENDIF

            !-----------------------------------------------------------
            ! If we are using CAC emissions and over CANADA ...
![eml
            IF ( LCAC .and. .not. LHIST ) THEN
!eml]
            IF ( GET_CANADA_MASK( I, J) > 0d0 ) THEN

               ! Read SO4 emissions in [molec/cm2/s]
               AN = GET_CAC_ANTHRO( I, J, IDTSO2, MOLEC_CM2_S=.TRUE. )

               ! Convert anthro SO2 to SO4 and from [molec/cm2/s] to 
               ! [kg/box/s]
               ! Place all CAC anthro SO4 into surface layer
               IF ( LNEI99 ) then
                  if (GET_USA_MASK( I, J) > 0d0 ) THEN
                     SO4an(I,J,1) = SO4an(I,J,1) + AN * 0.014/ 0.986 
     &                    * AREA_CM2 / XNUMOL(IDTSO4)
#if defined( TOMAS )
     &                           *scaleanthso2
#endif
                  endif
               ELSE
                  SO4an(I,J,1) = AN * 0.014 / 0.986 * AREA_CM2 
     &                           / XNUMOL(IDTSO4)
#if defined( TOMAS )
     &                           *scaleanthso2
#endif
                  SO4bf(I,J)   = 0d0
               ENDIF

               DO L = 2, NOXLEVELS
                  SO4an(I,J,L) = 0d0
               ENDDO

            ENDIF
            ENDIF

            !-----------------------------------------------------------
            ! If we are using BRAVO emissions and over MEXICO ...
![eml
            IF ( LBRAVO .and. .not. LHIST ) THEN
!eml]
            IF ( GET_BRAVO_MASK( I, J) > 0d0 ) THEN

               ! Read SO4 emissions in [molec/cm2/s]
               AN = GET_BRAVO_ANTHRO( I, J, IDTSO2 )
               
               ! Convert anthro SO2 to SO4 and from [molec/cm2/s] to 
               ! [kg/box/s]
               ! Place all BRAVO anthro SO4 into surface layer
               IF ( LNEI99 ) then
                  if (GET_USA_MASK( I, J) > 0d0 ) THEN
                     SO4an(I,J,1) = SO4an(I,J,1) + AN * 0.014/ 0.986 
     &                    * AREA_CM2 / XNUMOL(IDTSO4)
                  endif
               ELSE
                  SO4an(I,J,1) = AN * 0.014 / 0.986 * AREA_CM2 
     &                           / XNUMOL(IDTSO4)
                  SO4bf(I,J)   = 0d0
               ENDIF

               DO L = 2, NOXLEVELS
                  SO4an(I,J,L) = 0d0
               ENDDO

            ENDIF
            ENDIF

            !-----------------------------------------------------------
            ! If we are using NEI 2005 over the USA ...
            ! Must be called after CAC and BRAVO to simply overwrite
            ! where they overlap
![eml
            IF ( LNEI05 .and. .not. LHIST ) THEN
!eml]
               IF ( NEI05_MASK( I, J) > 0d0 ) THEN

                  ! Determine if we should use weekday or weekend NEI
                  ! emissions at grid box (I,J,L).  Since NEI is over
                  ! the US, then weekend is Sat/Sun.
                  DOW_LT  = GET_DAY_OF_WEEK_LT( I, J, 1 )
                  WEEKDAY = ( DOW_LT > 0 .and. DOW_LT < 6 )

                  ! Read SO4 emissions in [molec/cm2/s]
                  DO L = 1, NOXLEVELS
                     EPA_AN    = GET_NEI2005_ANTHRO( I, J, L, IDTSO4,
     &                                  WEEKDAY,  MOLEC_CM2_S=.TRUE. )
                     SO4an(I,J,L) = EPA_AN * AREA_CM2 / XNUMOL(IDTSO4)

                  ENDDO
                  SO4bf(I,J)   = 0d0
                  
               ENDIF
            ENDIF


            !=============================================================
            ! Do PBL mixing of emissions
            !=============================================================
            TEMPEMISS(1) = SO4an(I,J,1) + SO4bf(I,J)

            TEMPEMISS(2:)= SO4an(I,J,2:)

            call SULFATE_PBL_MIX( TEMPEMISS,   SO4, 
     $                            FEMIS, NTOP, IS_LOCAL )


            
            !=============================================================
            ! Add SO4 emissions to tracer array 
            ! Convert from [kg SO4/box/s] -> [kg SO4/box/timestep]
            !=============================================================
            IF ( IS_LOCAL ) THEN
               TC(I,J,:)  = TC(I,J,:)  + SO4(:)  * DTSRCE
            ELSE
               TC(I,J,2:) = TC(I,J,2:) + SO4(2:) * DTSRCE
               EMIS_SAVE(I,J,IDTSO4) = SO4(1) * DTSRCE
            ENDIF
               

            !==============================================================
            ! ND13 Diagnostic: SO4 emission in [kg S/box/timestep]       
            !==============================================================
            IF ( ND13 > 0 ) THEN 

               ! Anthro SO4
               AD13_SO4_an(I,J,:) = AD13_SO4_an(I,J,:) +
     $                              ( SO4an(I,J,:) * S_SO4 * DTSRCE )

               ! Biofuel SO4
               AD13_SO4_bf(I,J) = AD13_SO4_bf(I,J) + 
     &                            ( SO4bf(I,J) * S_SO4 * DTSRCE ) 
            ENDIF
            
            
         ENDDO
#if defined( DEVEL )
         State_Chm%TRAC_TEND(:,J,1,IDTSO4)           = 
     &   State_Chm%TRAC_TEND(:,J,1,IDTSO4)           + 
     &        ( SO4an(:,J,1) + SO4bf(:,J) )

         State_Chm%TRAC_TEND(:,J,2:NOXLEVELS,IDTSO4) = 
     &   State_Chm%TRAC_TEND(:,J,2:NOXLEVELS,IDTSO4) + SO4an(:,J,2:)
#endif
      ENDDO
!$OMP END PARALLEL DO

      ! Free pointer
      NULLIFY( TC )

      END SUBROUTINE SRCSO4
!EOC

#if   defined( TOMAS )
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: srcsf30
!
! !DESCRIPTION: !  Subroutine SRCSF30 handles 30bin sulfate SF1-SF30 and
!  number NK1-NK30 emissions into the GEOS-CHEM tracer array.(win, 1/25/10)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE SRCSF30( TC1, TC2, Input_Opt, State_Met )
!
! !USES:
!
      USE CAC_ANTHRO_MOD,     ONLY : GET_CANADA_MASK
      USE CAC_ANTHRO_MOD,     ONLY : GET_CAC_ANTHRO
      USE DIAG_MOD,           ONLY : AD59_SULF,     AD59_NUMB
      USE EPA_NEI_MOD,        ONLY : GET_EPA_ANTHRO,  GET_EPA_BIOFUEL
      USE EPA_NEI_MOD,        ONLY : GET_USA_MASK
      USE ERROR_MOD,          ONLY : ERROR_STOP,  IT_IS_NAN
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE GRID_MOD,           ONLY : GET_AREA_CM2
      USE PBL_MIX_MOD,        ONLY : GET_FRAC_OF_PBL, GET_PBL_TOP_L
      USE PRESSURE_MOD,       ONLY : GET_PCENTER
      USE TIME_MOD,           ONLY : GET_DAY_OF_WEEK
      USE TIME_MOD,           ONLY : GET_DAY_OF_WEEK_LT
      USE TIME_MOD,           ONLY : GET_TS_EMIS
      USE TOMAS_MOD,          ONLY : IBINS, AVGMASS, ICOMP
      USE TOMAS_MOD,          ONLY : Xk,  SUBGRIDCOAG, MNFIX
      USE TOMAS_MOD,          ONLY : SRTSO4, SRTNH4,  DEBUGPRINT
      USE TRACER_MOD,         ONLY : XNUMOL
      USE TRACERID_MOD        

      USE CMN_SIZE_MOD             ! Size parameters
      USE CMN_DIAG_MOD             ! ND13 (for now)
      USE CMN_GCTM_MOD             ! SCALE_HEIGHT
!
! !INPUT PARAMETERS: 
!
      TYPE(OptInput), INTENT(IN)  :: Input_Opt   ! Input Options object
      TYPE(MetState), INTENT(IN)  :: State_Met   ! Meteorology State object
!
! !INPUT/OUTPUT PARAMETERS:
!
      ! TC1 (REAL*8 ) : Array for 30-bin NK1-NK30 tracer number [no.]
      ! TC2 (REAL*8)  : Array for 30-bin SF1-SF30 tracer mass [kg] 
      REAL*8,  INTENT(INOUT) :: TC1(IIPAR,JJPAR,LLPAR,IBINS)
      REAL*8,  INTENT(INOUT) :: TC2(IIPAR,JJPAR,LLPAR,IBINS*ICOMP)
!
! !REMARKS:
!  NOTE: This subroutine is ignored unless you compile for TOMAS microphysics.
!                                                                             .
!  Variable comments
!  ============================================================================
!  NdistInit(IBINS)    : Initial size dist of emission
!  NdistFinal(IBINS)   : Final size dist of emissions after subgrid coag 
!  Ndist(IBINS)        : Number size dist of grid cell  
!  Mdist(IBINS,ICOMP)  : Mass size dist of grid cell
!  MaddFinal(NBINS)    : madd added to bins due to coagulation
!  tscale              : mixing time (s) currently assumed 10 hr. (win, 10/4/07)
!
! !REVISION HISTORY:
!  (1 ) For now use the call subroutine SRCSO4 as the source of bulk mass
!        emitted (win, 7/17/07) ---> did not work.  I need the SO4an and SO4bf
!  (2 ) Copy the subroutine SRCSO4 to use here. It may seem unwise to have the
!        same code in two subroutine (if have to modify in the future, have to
!        do at both places), but let's have it this way for now in case I can
!        run the microphysics run without needing to have bulk SO4 tracer
!        (win, 7/18/07)
!  (3 ) Overwrite SO4 emission of EPA/NEI99 inventory with a calculation using
!        an assumption that EPA/NEI99 SO2 emission is 99% of total S emission
!        (by mole) then calculate SO4 to be 1% of the total (win, 10/4/07)
!  (4 ) Add subgrid coagulation calculation (win, 10/4/07)
!  (5 ) Import SRCSF30 from geoschem v.7-04-11 to v.8-02-02 and update the code
!        for CAC emissions and Lin's PBL scheme (win, 1/25/10)
!  16 Feb 2012 - R. Yantosca - Added ProTeX headers 
!  01 Mar 2012 - R. Yantosca - Now use GET_AREA_CM2(I,J,L) from grid_mod.F90
!  09 Nov 2012 - M. Payer    - Replaced all met field arrays with State_Met
!                              derived type object
!  14 Jun 2013 - R. Yantosca - Now determine weekday/weekend with respect to
!                              the local time at each grid box.  (Formerly,
!                              this had been done w/r/t the GMT time.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL                :: WEEKDAY
      INTEGER                :: I, J, K, L, DOW_LT, NTOP, C
      REAL*8                 :: SO4(LLPAR), DTSRCE  , EFRAC(LLPAR)
      REAL*8                 :: TSO4,       FEMIS
      REAL*8                 :: AREA_CM2,   EPA_AN,  EPA_BF
      REAL*8                 :: CAC_AN
      REAL*8                 :: SO4an(IIPAR,JJPAR,2)
      REAL*8                 :: SO4bf(IIPAR,JJPAR)
      REAL*8                 :: SO4anbf(IIPAR,JJPAR,2)
      REAL*8 BFRAC(IBINS)     ! Mass fraction emitted to each bin

#if  defined( TOMAS12) || defined( TOMAS15)
      DATA BFRAC/
# if  defined( TOMAS15)
     &           0.0d0     , 0.0d0     , 0.0d0, 
# endif
     &           4.3760E-02, 6.2140E-02, 3.6990E-02, 1.8270E-02,
     &           4.2720E-02, 1.1251E-01, 1.9552E-01, 2.2060E-01,
     &           1.6158E-01, 7.6810E-02, 2.8884E-02, 2.0027E-04/

#else
!else we are using 30 or 40 bin model

      DATA BFRAC/
# if  defined( TOMAS40)
     &           0.000d00 , 0.000d00 , 0.000d00 , 0.000d00 , 0.000d00 ,
     &           0.000d00 , 0.000d00 , 0.000d00 , 0.000d00 , 0.000d00 ,
# endif
     &           1.728d-02, 2.648d-02, 3.190d-02, 3.024d-02, 2.277d-02,
     &           1.422d-02, 9.029d-03, 9.241d-03, 1.531d-02, 2.741d-02,
     &           4.529d-02, 6.722d-02, 8.932d-02, 1.062d-01, 1.130d-01,
     &           1.076d-01, 9.168d-02, 6.990d-02, 4.769d-02, 2.912d-02,
     &           1.591d-02, 7.776d-03, 3.401d-03, 1.331d-03, 4.664d-04,
     &           1.462d-04, 4.100d-05, 1.029d-05, 2.311d-06, 4.645d-07/
#endif

      REAL*8                 :: NDISTINIT(IBINS) 
      REAL*8                 :: NDISTFINAL(IBINS) 
      REAL*8                 :: MADDFINAL(IBINS) 
      REAL*8                 :: NDIST(IBINS),  MDIST(IBINS,ICOMP)
      REAL*8                 :: NDIST2(IBINS), MDIST2(IBINS,ICOMP)
      REAL*4                 :: TSCALE, BOXVOL, TEMP, PRES

      REAL*8                 :: N0(LLPAR,IBINS), M0(LLPAR,IBINS)

      LOGICAL                :: ERRORSWITCH, SGCOAG = .TRUE.
      logical                :: pdbug !(temporary) win, 10/24/07
      integer                :: ii, jj, ll
      data ii, jj, ll / 61, 1, 7 /

      ! Ratio of molecular weights: S/SO4
      REAL*8,  PARAMETER     :: S_SO4 = 32d0 / 96d0

      ! debugging 
      real*8   dummy

      ! For fields from Input_Opt
      LOGICAL :: LNEI99,  LPRT, LCAC, LNLPBL

      !=================================================================
      ! SRCSF30 begins here!
      !=================================================================

      ! COpy values from Input_Opt
      LNEI99  = Input_Opt%LNEI99
      LPRT    = Input_Opt%LPRT
      LCAC    = Input_Opt%LCAC
      LNLPBL  = Input_Opt%LNLPBL

      ! DTSRCE is the emission timestep in seconds
      DTSRCE  = GET_TS_EMIS() * 60d0

      !=================================================================
      ! Overwrite USA w/ EPA/NEI99 SO4 emissions (if necessary)
      ! Overwrite CANADA w/ CAC       SO2-fraction emiss. (if necessary)
      ! Store emissions into local arrays SO4an, SO4bf
      !=================================================================
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, AREA_CM2, EPA_AN, EPA_BF, CAC_AN, DOW_LT, WEEKDAY )
      DO J = 1, JJPAR
         
         DO I = 1, IIPAR

            ! Grid box surface area [cm2]
            AREA_CM2 = GET_AREA_CM2( I, J, 1 )

            !-----------------------------------------------------------
            ! Default SO4 from GEIA, or (as a fraction of SO2) from 
            ! EDGAR w/ optional STREETS for S.E.-ASIA, and optional
            ! EMEP for Europe
            !-----------------------------------------------------------
            SO4an(I,J,1) = ESO4_an(I,J,1)
            SO4an(I,J,2) = ESO4_an(I,J,2)
            SO4bf(I,J)   = 0d0

            ! If we are using EPA/NEI99 emissions
            IF ( LNEI99 ) THEN
            IF ( GET_USA_MASK( I, J ) > 0d0 ) THEN

               ! Determine if we should use weekday or weekend NEI
               ! emissions at grid box (I,J,L).  Since NEI is over
               ! the US, then weekend is Sat/Sun.
               DOW_LT  = GET_DAY_OF_WEEK_LT( I, J, 1 )
               WEEKDAY = ( DOW_LT > 0 .and. DOW_LT < 6 )

               ! Calculate SO4 emissions from assuming the SO2 emitted
               ! is 99% of total mass (win, 10/4/07)
               EPA_AN = GET_EPA_ANTHRO(I,J,IDTSO2, WEEKDAY)/  
     &                  0.99d0*0.01d0
               EPA_BF = GET_EPA_BIOFUEL(I,J,IDTSO2, WEEKDAY)/
     &                  0.990*0.01d0

               ! Convert anthro SO4 from [molec/cm2/s] to [kg/box/s] 
               ! Place all EPA/NEI99 anthro SO4 into surface layer
               SO4an(I,J,1) = EPA_AN * AREA_CM2 / XNUMOL(IDTSO4)
               SO4an(I,J,2) = 0d0
               
               ! Convert biofuel SO4 from [molec/cm2/s] to [kg/box/s]
               SO4bf(I,J)   = EPA_BF * AREA_CM2 / XNUMOL(IDTSO4)

            ENDIF
            ENDIF

            ! If we are using CAC emissions and over CANADA ...
            IF ( LCAC ) THEN
            IF ( GET_CANADA_MASK( I, J) < 0d0 ) THEN

               ! Read SO4 emissions in [molec/cm2/s]
               CAC_AN       = GET_CAC_ANTHRO( I, J, IDTSO2, 
     &                                        MOLEC_CM2_S=.TRUE. )

               ! Calculate SO4 emissions from assuming the SO2 emitted
               ! is 99% of total mass
               CAC_AN  = CAC_AN / 0.99d0*0.01d0

               ! Convert anthro SO2 to SO4 and from [molec/cm2/s] to 
               ! [kg/box/s]
               ! Place all CAC anthro SO4 into surface layer
               IF ( LNEI99 .and. GET_USA_MASK( I, J) > 0d0 ) THEN
                  SO4an(I,J,1) = SO4an(I,J,1) + CAC_AN * 0.014/ 0.986 
     &                           * AREA_CM2 / XNUMOL(IDTSO4)
               ELSE
                  SO4an(I,J,1) = CAC_AN * 0.014 / 0.986 * AREA_CM2 
     &                           / XNUMOL(IDTSO4)
                  SO4bf(I,J)   = 0d0
               ENDIF

               SO4an(I,J,2) = 0d0

            ENDIF
            ENDIF 

         ENDDO
      ENDDO
!$OMP END PARALLEL DO

      !=================================================================
      ! Compute SO4 emissions 
      !=================================================================
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, NTOP, SO4, TSO4, L, FEMIS, EFRAC, N0, M0, K )
!$OMP+PRIVATE( NDISTINIT, NDIST, MDIST, NDISTFINAL, MADDFINAL )
!$OMP+PRIVATE( NDIST2, MDIST2, C )
!$OMP+PRIVATE( BOXVOL, TEMP, PRES, TSCALE, pdbug )
!$OMP+SCHEDULE( DYNAMIC )
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Top level of boundary layer at (I,J)
         NTOP = CEILING( GET_PBL_TOP_L( I, J ) )

         ! Zero SO4 array at all levels 
         DO L = 1, LLPAR
            SO4(L) = 0.0
         ENDDO

         ! Compute total anthro SO4 (surface + 100m) plus biofuel SO4
         TSO4 = 0.d0
         TSO4 = SUM( SO4an(I,J,:) ) + SO4bf(I,J)
         IF ( TSO4 == 0d0 ) CYCLE


!debug         if(i==60.and.j==35) print *,'TSO4',TSO4
         !=============================================================
         ! First calculate emission distribution vertically within PBL
         !=============================================================
         ! EFRAC(30) = fraction of total emission splitted for each
         !             level until reaching PBL top.
         EFRAC = 0d0

         !==============================================================
         ! Partition the total anthro SO4 emissions thru the entire 
         ! boundary layer (if PBL top is higher than level 2)
         !==============================================================
         ! Add option for non-local PBL (Lin, 03/31/09)
         IF (.NOT. LNLPBL) THEN
         IF ( NTOP > 2 ) THEN

            ! Loop thru boundary layer
            DO L = 1, NTOP
                 
               ! Fraction of PBL spanned by grid box (I,J,L) [unitless]
               EFRAC(L)  = GET_FRAC_OF_PBL( I, J, L )
!               print *, 'L, EFRAC(L)',L, EFRAC(L)
            
            ENDDO
         ELSE
            EFRAC(1) = ( SO4an(I,J,1) + SO4bf(I,J) ) / TSO4
            EFRAC(2) = SO4an(I,J,2) / TSO4
         ENDIF

         IF ( ABS( SUM( EFRAC(:)) - 1.d0 ) > 1.D-5 ) THEN
            PRINT*, '### ERROR in SRCSF30!'
            PRINT*, '### I, J : ', I, J
            print*, 'EFRAC',EFRAC(:)
            PRINT*, '### SUM(EFRAC) : ', SUM( EFRAC(:) )
            PRINT*, '### This should exactly 1.00'
            CALL ERROR_STOP( 'Check SO4 redistribution',
     &                       'SRCSF30 (sulfate_mod.f)' )
         ENDIF

         ELSE 
            ! stop the program for now since I don't totally implement
            ! the subgrid coagulation option w/ Lin's new PBL scheme (win, 1/25/10)
            print *,'If the program stops here, that means you are ',
     &           'running TOMAS simulation with the new PBL scheme ',
     &           'implemented since GEOS-Chem v.8-02-01.', 
     &           '-----> Try not using the non-local PBL option'
            CALL ERROR_STOP( 'Code does not support new PBL scheme',
     &                       'SRCSF30 (sulfate_mod.f)')
         
         ENDIF  ! .not. LNLPBL

         !=============================================================
         ! Add the size-resolved SO4 emission to tracer array
         ! Having the options to do sub-grid coagulation or simply 
         ! emit.  
         ! Sub-grid coagulation reduces the number being emitted
         ! and modifies the mass size distribution of existing particle
         ! as well as the size distribution being emitted.
         ! (win, 10/4/07)
         !=============================================================
         IF ( SGCOAG ) THEN

            N0(:,:) = TC1(I,J,:,:)
            M0(:,:) = TC2(I,J,:,1:IBINS)

          DO L = 1, LLPAR
            SO4(L) = TSO4 * EFRAC(L) * DTSRCE
            IF ( SO4(L) == 0.d0 ) CYCLE
            DO K = 1, IBINS
               NDISTINIT(K) = SO4(L) * BFRAC(K) /
     &                                ( SQRT( XK(K)*XK(K+1) ) )
               NDIST(K) = TC1(I,J,L,K)
               DO C = 1, ICOMP
                  MDIST(K,C) = TC2(I,J,L,K+(C-1)*IBINS)
                  IF( IT_IS_NAN( MDIST(K,C) ) ) THEN
                     PRINT *,'+++++++ Found NaN in SRCSF30  +++++++'
                     PRINT *,'Location (I,J,L):',I,J,L,'Bin',K,'comp',C
                     CALL  ERROR_STOP('SRCSF30 SGCCOAG','sulfate_mod.f')
                  ENDIF
               ENDDO
               NDISTFINAL(K) = 0
               MADDFINAL(K)  = 0
            ENDDO  

!      print *, 'mnfix in sulfate_mod:6346'
            CALL MNFIX( NDIST, MDIST, ERRORSWITCH )
            IF( ERRORSWITCH ) PRINT *,'SRCSF30: MNFIX found error ',
     &              'before SUBGRIDCOAG at ',I,J,L
            ERRORSWITCH = .FALSE.

            ! Overwrite number and mass before emission for diagnostic 
            ! just in case there was any change by MNFIX (win, 10/27/08)
            N0(L,:) = NDIST(:)
            M0(L,:) = MDIST(:,SRTSO4)

!debug            if (I== ii .and. J==jj .and. L==ll)  
!     &           CALL DEBUGPRINT(NDIST,MDIST,I,J,L,
!     &           'SRCSF30 Before SUBGRIDCOAG ')


      ! Define subgrid coagulation timescale (win, 10/28/08)
#if   defined( GRID4x5   ) 
            TSCALE = 10.*3600.  ! 10 hours
#elif defined( GRID2x25  )
            TSCALE = 5.*3600.    
#elif defined( GRID1x125 )
            TSCALE = 2.*3600.    
#elif defined( GRID1x1   ) 
            TSCALE = 2.*3600.    
#elif defined( GRID05x0666 ) 
            TSCALE = 1.*3600.    
#elif defined( GRID025x03125 )
            !%%% KLUDGE, just copied the 0.5 x 0.666 value
            !%%% Someone needs to add the right value (bmy, 2/16/12)
            TSCALE = 1.*3600.    
#endif
            !Prior to 10/28/08 (win)     
            !TSCALE = 10.*3600.
            BOXVOL  = State_Met%AIRVOL(I,J,L) * 1.e6 !convert from m3 -> cm3 
            TEMP    = State_Met%T(I,J,L)
            PRES    = GET_PCENTER(i,j,l)*100.0 ! in Pa

!            print *, 'Now doing subgrid coag with timescale ',
!     &           TSCALE/3600.,'hr'

            pdbug=.false.
!debug            if (I.eq.II.and.J.eq.JJ .and. L==LL) pdbug=.true.
                             
            CALL SUBGRIDCOAG( NDISTINIT, NDIST, MDIST, BOXVOL,TEMP,
     &                        PRES, TSCALE, NDISTFINAL, MADDFINAL,pdbug) 

!debug            if (I== ii .and. J==jj .and. L==ll)then
!               print*,'PS'
!               print*,'SO4 emis',SO4(L)
!               print*,'Mdistnit',SO4(L)*BFRAC(:)
!               print*,'tscale=',tscale
!               print*,'ndistinit=',ndistinit
!               print*,'ndistfinal=',ndistfinal
!               print*,'maddfinal=',maddfinal
!            endif
!debug---
            DO K = 1, IBINS
               NDIST(K) = NDIST(K) + NDISTFINAL(K)
               MDIST(K,SRTSO4) = MDIST(K,SRTSO4) +
     &                           NDISTFINAL(K) * (SQRT(XK(K)*XK(K+1)))+
     &                           MADDFINAL(K)   

!debug               if (I.eq.II.and.J.eq.JJ .and. L==ll) then
!                  print *,(NDISTFINAL(K) * (SQRT(XK(K)*XK(K+1)))+
!     &                 MADDFINAL(K) )
!               endif
!debug---
            ENDDO
            DO K= 1, IBINS
               NDIST2(K) = NDIST(K)
               DO C = 1, ICOMP
                  MDIST2(K,C) = MDIST(K,C)
               ENDDO
            ENDDO
!debug            if(i==ii.and.j==jj.and.l==ll  ) 
!     &          CALL DEBUGPRINT(NDIST2,MDIST2,I,J,L, 
!     &                         'SRCSF30 : After SUBGRIDCOAG')

 
            ERRORSWITCH = .FALSE.
!debug            if(i==ii.and.j==ii.and.l==ll) errorswitch=.true.
!      print *, 'mnfix in sulfate_mod:6430'
            CALL MNFIX( NDIST2, MDIST2, ERRORSWITCH )

            IF( ERRORSWITCH ) PRINT *,'SRCSF30: MNFIX found error ',
     &           'after SUBGRIDCOAG at ',I,J,L

!debug            if(i==ii .and.j==jj .and.l==ll) 
!     &          CALL DEBUGPRINT(NDIST2,MDIST2,I,J,L, 
!     &                         'SRCSF30 : After MNFIX')

            DO K = 1, IBINS
               TC1(I,J,L,K) = NDIST2(K)
               DO C=1,ICOMP
                  TC2(I,J,L,K+(C-1)*IBINS) = MDIST2(K,C)
               ENDDO
            ENDDO

         ENDDO ! L loop


           !==============================================================
           ! ND59 Diagnostic: Size-resolved primary sulfate emission in 
           !                 [kg S/box/timestep] and the corresponding
           !                  number emission [no./box/timestep]
           !==============================================================
            IF ( ND59 > 0 ) THEN 
               DO L = 1, LLPAR
               DO K = 1, IBINS
                  if(TC2(I,J,L,K)-M0(L,K) < 0d0) 
     &              print *,'Negative SF emis ',TC2(I,J,L,K)-M0(L,K),
     &                 'at',I,J,L,K
                  if(TC1(I,J,L,K)-N0(L,K) < 0d0) 
     &              print *,'Negative NK emis ',TC1(I,J,L,K)-N0(L,K),
     &                 'at',I,J,L,K
                  AD59_SULF(I,J,1,K) = AD59_SULF(I,J,1,K) + 
     &                                  (TC2(I,J,L,K)-M0(L,K))*S_SO4
                  AD59_NUMB(I,J,1,K) = AD59_NUMB(I,J,1,K) + 
     &                                  TC1(I,J,L,K)-N0(L,K)
               ENDDO
               ENDDO
            ENDIF
           
         ELSE

         ! Distributing primary emission without sub-grid coagulation
         !=============================================================
         ! Add SO4 emissions to tracer array 
         ! For SF: Convert from [kg SO4/box/s] -> [kg SO4/box/timestep]
         ! For NK: Convert from [kg SO4/box/s] -> [No.   /box/timestep]
         !=============================================================
            DO L = 1, LLPAR
               SO4(L) = TSO4 * EFRAC(L)
               DO K = 1, IBINS
                  TC1(I,J,L,K) = TC1(I,J,L,K) + 
     &                    ( SO4(L) * DTSRCE * BFRAC(K) / AVGMASS(K) )
                  TC2(I,J,L,K) = TC2(I,J,L,K) + 
     &                    ( SO4(L) * DTSRCE * BFRAC(K)               )
               ENDDO
            ENDDO

           !==============================================================
           ! ND59 Diagnostic: Size-resolved primary sulfate emission in 
           !                 [kg S/box/timestep] and the corresponding
           !                  number emission [no./box/timestep]
           !==============================================================
            IF ( ND59 > 0 ) THEN 
               SO4anbf(:,:,1) = SO4an(:,:,1) + SO4bf(:,:)
               SO4anbf(:,:,2) = SO4an(:,:,2)
            
               DO L = 1, 2
               DO K = 1, IBINS
                  AD59_SULF(I,J,L,K) = AD59_SULF(I,J,L,K) + 
     &                                ( SO4anbf(I,J,L) * BFRAC(K) 
     &                                    * S_SO4 * DTSRCE         )
                  AD59_NUMB(I,J,L,K) = AD59_NUMB(I,J,L,K) + 
     &                                ( SO4anbf(I,J,L) * BFRAC(K) 
     &                                    / AVGMASS(K) * DTSRCE   )
               ENDDO
               ENDDO

            ENDIF

         ENDIF !SGCOAG

      ENDDO
      ENDDO
!$OMP END PARALLEL DO
      IF ( LPRT ) print *,'   ### Finish SRCSF30'

      END SUBROUTINE SRCSF30
!EOC
#endif
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: srcnh3
!
! !DESCRIPTION: Subroutine SRCNH3 handles NH3 emissions into the GEOS-CHEM
!  tracer array. (rjp, bmy, 12/17/01, 5/27/09)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE SRCNH3( am_I_Root, Input_Opt, State_Chm, RC )
!
! !USES:
!
      USE CAC_ANTHRO_MOD,     ONLY : GET_CANADA_MASK
      USE CAC_ANTHRO_MOD,     ONLY : GET_CAC_ANTHRO
      USE NEI2005_ANTHRO_MOD, ONLY : GET_NEI2005_ANTHRO
      USE NEI2005_ANTHRO_MOD, ONLY : NEI05_MASK => USA_MASK
      USE DIAG_MOD,           ONLY : AD13_NH3_an, AD13_NH3_bb
      USE DIAG_MOD,           ONLY : AD13_NH3_bf, AD13_NH3_na
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE GIGC_State_Chm_Mod, ONLY : ChmState
      USE GRID_MOD,           ONLY : GET_AREA_CM2
      USE EPA_NEI_MOD,        ONLY : GET_EPA_ANTHRO, GET_EPA_BIOFUEL
      USE EPA_NEI_MOD,        ONLY : GET_USA_MASK
      USE ERROR_MOD,          ONLY : ERROR_STOP
      USE PBL_MIX_MOD,        ONLY : GET_FRAC_OF_PBL, GET_PBL_TOP_L
      USE TIME_MOD,           ONLY : GET_DAY_OF_WEEK_LT
      USE TIME_MOD,           ONLY : GET_TS_EMIS
      USE TRACERID_MOD,       ONLY : IDTNH3

      USE CMN_SIZE_MOD             ! Size parameters
      USE CMN_DIAG_MOD             ! ND13
      USE CMN_GCTM_MOD             ! SCALE_HEIGHT
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
!  (1 ) Now save NH3 emissions to ND13 diagnostic (bmy, 12/13/02)
!  (2 ) Now reference AD13_NH3_na from "diag_mod.f", and archive natural 
!       source NH3 diagnostics for ND13.  Also consider natural source NH3
!       when partitioning by level into the STT array. (rjp, bmy, 3/23/03)
!  (3 ) Now use routine GET_TS_EMIS from the new "time_mod.f" (bmy, 3/27/03)
!  (4 ) For GEOS-4, convert PBL from [m] to [hPa] w/ the barometric law.
!       Now references SCALE_HEIGHT from "CMN_GCTM".  Added BLTHIK variable
!       to hold PBL thickness in [hPa]. (bmy, 1/15/04)
!  (5 ) Now references GET_EPA_ANTHRO, GET_EPA_BIOFUEL, and GET_USA_MASK from 
!       "epa_nei_mod.f".  Now references GET_DAY_OF_WEEK from "time_mod.f".  
!       Now references LNEI99 from "logical_mod.f".  Now references 
!       GET_AREA_CM2 from "grid_mod.f".  Now references IDTNH3 from 
!       "tracerid_mod.f".  Now references XNUMOL from CMN_O3.  Now can 
!       overwrite the anthro & biofuel NH3 emissions over the continental US 
!       if LNEI99=T.  Now references IDTNH3 from "tracerid_mod.f". 
!       (rjp, rch, bmy, 11/16/04)
!  (6 ) Remove reference to "pressure_mod.f".  Now reference GET_FRAC_OF_PBL 
!       and GET_PBL_TOP_L from "pbl_mix_mod.f".  Removed reference to header 
!       file CMN. (bmy, 2/22/05)
!  (7 ) Now references XNUMOL from "tracer_mod.f" (bmy, 10/25/05)
!  (8 ) Need to add CAC_AN to the PRIVATE loop (bmy, 5/27/09)
!  (9 ) Added NIE 2005 (amv, 10/07/2009) 
!  (10) Made NH3an 3D; Calls SULFATE_PBL_MIX to do the PBL mixing of
!       emissions, and allows for emissions above the PBL (phs, 10/27/09)
!  22 Dec 2011 - M. Payer    - Added ProTeX headers
!  01 Mar 2012 - R. Yantosca - Now use GET_AREA_CM2(I,J,L) from grid_mod.F90
!  25 Mar 2013 - R. Yantosca - Now accept am_I_Root, Input_Opt, State_Chm, RC
!  14 Jun 2013 - R. Yantosca - Now determine weekday/weekend with respect to
!                              the local time at each grid box.  (Formerly,
!                              this had been done w/r/t the GMT time.
!  22 Jul 2013 - M. Sulprizio- Now copy LRCP from Input_Opt
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL                :: WEEKDAY,  IS_LOCAL
      INTEGER                :: I, J, L,  K, NTOP, DOW_LT
      REAL*8                 :: DTSRCE,   AREA_CM2    
      REAL*8                 :: CAC_AN,   EPA_AN,  EPA_BF
      REAL*8                 :: FEMIS(LLPAR), NH3(LLPAR)
      REAL*8                 :: NH3an(IIPAR,JJPAR,NOXLEVELS)
      REAL*8                 :: NEI2005(IIPAR,JJPAR,NOXLEVELS)
      REAL*8                 :: NH3bf(IIPAR,JJPAR)
      REAL*8                 :: TEMPEMISS(NOXLEVELS)
      
      ! For fields from Input_Opt
      LOGICAL                :: LNEI99, LCAC, LNEI05, LNLPBL, LRCP
      INTEGER                :: N_TRACERS
      REAL*8                 :: XNUMOL(Input_Opt%N_TRACERS)

      ! Pointers
      ! We need to define local arrays to hold corresponding values 
      ! from the Chemistry State (State_Chm) object. (mpayer, 12/6/12)
      REAL*8, POINTER        :: TC(:,:,:)

      !=================================================================
      ! SRCNH3 begins here!
      !=================================================================

      ! Assume success
      RC        =  GIGC_SUCCESS

      ! Copy values from Input_Opt
      LNLPBL    = Input_Opt%LNLPBL
      LNEI99    = Input_Opt%LNEI99
      LCAC      = Input_Opt%LCAC
      LNEI05    = Input_Opt%LNEI05
      LRCP      = Input_Opt%LRCP
      N_TRACERS = Input_Opt%N_TRACERS
      XNUMOL    = Input_Opt%XNUMOL(1:N_TRACERS)

      ! Point to the NH3 tracer slice
      TC        => State_Chm%TRACERS(:,:,:,IDTNH3)

      ! Emission timestep [s]
      DTSRCE    = GET_TS_EMIS() * 60d0
 
      ! init
      NH3an     = 0D0
      
      ! is it local PBL mixing scheme?
      IS_LOCAL   = .NOT. LNLPBL

      !=================================================================
      ! Overwrite USA with EPA/NEI NH3 emissions (if necessary)
      ! Store emissions into local arrays NH3an, NH3bf
      !=================================================================
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, AREA_CM2, EPA_AN, EPA_BF, CAC_AN, L )
!$OMP+PRIVATE( NTOP, FEMIS, NH3, TEMPEMISS, DOW_LT, WEEKDAY )
      DO J = 1, JJPAR

         DO I = 1, IIPAR

            !----------------------------------------------------------------
            ! NOTE: There seems to be some problems with the EPA/NEI NH3 
            ! emissions.  Therefore we will use the existing emissions for NH3 
            ! until further notice.  Comment out the lines below until 
            ! further notice.  (bmy, 11/17/04) 
            !! Grid box surface area [cm2]
            AREA_CM2 = GET_AREA_CM2( I, J, 1 )
            !----------------------------------------------------------------

            ! Top of the boundary layer
            NTOP = CEILING( GET_PBL_TOP_L( I, J ) ) 

            ! Fraction of PBL spanned by grid box (I,J,L) [unitless]
            DO L = 1, LLPAR
               FEMIS(L)  = GET_FRAC_OF_PBL( I, J, L )
            ENDDO

            !-----------------------------------------------------------------
            ! NOTE: There seems to be some problems with the EPA/NEI NH3 
            ! emissions.  Therefore we will use the existing emissions for 
            ! NH3 until further notice.  Comment out the lines below until 
            ! further notice.  (bmy, 11/17/04) 
            !
            !! If we are using EPA/NEI99 emissions ...
            !IF ( LNEI99 ) THEN
            !
            !   ! If we are over the USA ...
            !   IF ( GET_USA_MASK( I, J ) > 0d0 ) THEN
            !
            !      ! Determine if we should use weekday or weekend NEI
            !      ! emissions at grid box (I,J,L).  Since NEI is over
            !      ! the US, then weekend is Sat/Sun.
            !      DOW_LT     = GET_DAY_OF_WEEK_LT( I, J, 1 )
            !      WEEKDAY    = ( DOW_LT > 0 .and. DOW_LT < 6 )
            !
            !      ! Read NH3 anthro emissions in [molec NH3/cm2/s]
            !      EPA_AN     = GET_EPA_ANTHRO(  I, J, IDTNH3, WEEKDAY )
            !      EPA_BF     = GET_EPA_BIOFUEL( I, J, IDTNH3, WEEKDAY )
            !
            !      ! Convert from [molec NH3/cm2/s] to [kg NH3/box/sec]
            !      NH3an(I,J) = EPA_AN * AREA_CM2 / XNUMOL(IDTNH3)
            !      NH3bf(I,J) = EPA_BF * AREA_CM2 / XNUMOL(IDTNH3)
            !
            !   ELSE
            !
            !      ! If we are not over the USA, just use the regular 
            !      ! emissions in NH3_an and NH3bf (bmy, 11/16/04)
            !      NH3an(I,J) = ENH3_an(I,J)
            !      NH3bf(I,J) = ENH3_bf(I,J)
            !
            !   ENDIF
            !
            !ELSE
            !-----------------------------------------------------------------

            ! If we are not using the EPA/NEI emissions, just copy the 
            ! regular ENH3_an and ENH3_bf to local arrays. (bmy, 11/16/04)
            NH3an(I,J,1) = ENH3_an(I,J)
            NH3bf(I,J)   = ENH3_bf(I,J)
            
            ! The RCP anthropogenic emission inventory already includes
            ! biofuel emissions, so we turn off the additional biofuel 
            ! inventory (cdh, 8/31/12)
            IF ( LRCP ) NH3bf(I,J) = 0d0 
            
            !-----------------------------------------------------------------
            ! NOTE: There seems to be some problems with the EPA/NEI NH3 
            ! emissions.  Therefore we will use the existing emissions for 
            ! NH3 until further notice.  Comment out the lines below until 
            ! further notice.  (bmy, 11/17/04) 
            !ENDIF
            !-----------------------------------------------------------------

            ! If we are using CAC emissions and over CANADA
            IF ( LCAC ) THEN
            IF ( GET_CANADA_MASK( I, J ) > 0d0 ) THEN

               ! Read NH3 anthro emissions in [molec NH3/cm2/s]
               CAC_AN = GET_CAC_ANTHRO( I, J, IDTNH3, 
     &                                  MOLEC_CM2_S=.TRUE.)

               ! Convert from [molec NH3/cm2/c] to [kg NH3/box/sec]
               NH3an(I,J,1) = CAC_AN * AREA_CM2 / XNUMOL(IDTNH3)

            ENDIF
            ENDIF

            

            ! If we are using NEI 2005 over North America
            IF ( LNEI05 ) THEN
            IF ( NEI05_MASK( I, J ) > 0d0 ) THEN
             
               ! Determine if we should use weekday or weekend NEI
               ! emissions at grid box (I,J,L).  Since NEI is over
               ! the US, then weekend is Sat/Sun.
               DOW_LT  = GET_DAY_OF_WEEK_LT( I, J, 1 )
               WEEKDAY = ( DOW_LT > 0 .and. DOW_LT < 6 )

               DO L = 1, NOXLEVELS 

                  ! Read NH3 anthro emissions in [molec NH3/cm2/s]                  
                  EPA_AN = GET_NEI2005_ANTHRO( I, J, L, IDTNH3, 
     &                                   WEEKDAY, MOLEC_CM2_S=.TRUE.)
     
                  ! Convert from [molec NH3/cm2/c] to [kg NH3/box/sec]
                  NH3an(I,J,L) = EPA_AN * AREA_CM2 / XNUMOL(IDTNH3)

               ENDDO

            ENDIF
            ENDIF



            !=============================================================
            ! Do PBL mixing of all types of NH3 emission [kg/box/s]
            !=============================================================
            TEMPEMISS(1) = NH3an(I,J,1) + ENH3_bb(I,J) + 
     &                     NH3bf(I,J)   + ENH3_na(I,J)

            TEMPEMISS(2:)= NH3an(I,J,2:)
            
            call SULFATE_PBL_MIX( TEMPEMISS, NH3, 
     $                            FEMIS,     NTOP, IS_LOCAL )


            !=============================================================
            ! Add NH3 emissions to tracer/emission arrays, and
            ! convert from [kg NH3/box/s] -> [kg NH3/box/timestep]
            !=============================================================
            IF ( IS_LOCAL ) THEN
               TC(I,J,:)  = TC(I,J,:)  + NH3(:)  * DTSRCE
            ELSE
               TC(I,J,2:) = TC(I,J,2:) + NH3(2:) * DTSRCE
               EMIS_SAVE(I,J,IDTNH3) = NH3(1) * DTSRCE
            ENDIF
            
            !============================================================
            ! ND13 diagnostics: NH3 emissions [kg NH3/box/timestep]
            !============================================================
            IF ( ND13 > 0 ) THEN

               ! Anthro NH3
               DO L = 1, NOXLEVELS      
                  AD13_NH3_an(I,J,L) = AD13_NH3_an(I,J,L) + 
     &                              ( NH3an(I,J,L)   * DTSRCE )            
               ENDDO

               ! Biomass NH3
               AD13_NH3_bb(I,J) = AD13_NH3_bb(I,J) + 
     &                         ( ENH3_bb(I,J) * DTSRCE )

               ! Biofuel NH3
               AD13_NH3_bf(I,J) = AD13_NH3_bf(I,J) +
     &                         ( NH3bf(I,J)   * DTSRCE )

               ! Natural source NH3
               AD13_NH3_na(I,J) = AD13_NH3_na(I,J) + 
     &                         ( ENH3_na(I,J) * DTSRCE )


            ENDIF
            
         ENDDO
#if defined( DEVEL )
         State_Chm%TRAC_TEND(:,J,1,IDTNH3) = 
     &   State_Chm%TRAC_TEND(:,J,1,IDTNH3) + 
     &        ( NH3an(:,J,1) + ENH3_bb(:,J) + 
     &          NH3bf(:,J)   + ENH3_na(:,J) )

         State_Chm%TRAC_TEND(:,J,2:NOXLEVELS,IDTNH3) = 
     &   State_Chm%TRAC_TEND(:,J,2:NOXLEVELS,IDTNH3) + NH3an(:,J,2:)
#endif
      ENDDO
!$OMP END PARALLEL DO

      ! Free pointer
      NULLIFY( TC )

      END SUBROUTINE SRCNH3
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_oh
!
! !DESCRIPTION: Function GET\_OH returns OH from SMVGEAR's CSPEC array (for
!  coupled runs) or monthly mean OH (for offline runs).  Imposes a diurnal 
!  variation on OH for offline simulations. (bmy, 12/16/02, 7/20/04)
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_OH( I, J, L, Input_Opt, State_Met ) 
     &         RESULT( OH_MOLEC_CM3 )
!
! !USES:
!
      USE CMN_SIZE_MOD
      USE COMODE_MOD,         ONLY : CSPEC, JLOP
      USE ERROR_MOD,          ONLY : ERROR_STOP
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE GLOBAL_OH_MOD,      ONLY : OH
      USE TIME_MOD,           ONLY : GET_TS_CHEM
      USE TRACERID_MOD,       ONLY : IDOH
!
! !INPUT PARAMETERS: 
!
      INTEGER,        INTENT(IN) :: I, J, L     ! Lon, lat, level indices
      TYPE(OptInput), INTENT(IN) :: Input_Opt   ! Input Options object
      TYPE(MetState), INTENT(IN) :: State_Met   ! Meteorology State object
!
! !REVISION HISTORY: 
!  (1 ) We assume SETTRACE has been called to define IDOH (bmy, 11/1/02)
!  (2 ) Now use function GET_TS_CHEM from "time_mod.f" (bmy, 3/27/03)
!  (3 ) Now reference ITS_A_FULLCHEM_SIM, ITS_AN_AEROSOL_SIM from 
!        "tracer_mod.f".  Also replace IJSURF w/ an analytic function. 
!        (bmy, 7/20/04)
!  22 Dec 2011 - M. Payer    - Added ProTeX headers
!  28 Nov 2012 - R. Yantosca - Replace SUNCOS with State_Met%SUNCOS
!  28 Nov 2012 - R. Yantosca - Add State_Met to argument list
!  04 Mar 2013 - R. Yantosca - Now pass Input_Opt%ITS_A_FULLCHEM_SIM and
!                              Input_Opt%ITS_AN_AEROSOL_SIM
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER             :: JLOOP
      REAL*8              :: OH_MOLEC_CM3
 
      !=================================================================
      ! GET_OH begins here!
      !=================================================================

      IF ( Input_Opt%ITS_A_FULLCHEM_SIM ) THEN

         !---------------------
         ! Coupled simulation
         !---------------------

         ! JLOOP = SMVGEAR 1-D grid box index
         JLOOP = JLOP(I,J,L)

         ! Take OH from the SMVGEAR array CSPEC
         ! OH is defined only in the troposphere
         IF ( JLOOP > 0 ) THEN
            OH_MOLEC_CM3 = CSPEC(JLOOP,IDOH)
         ELSE
            OH_MOLEC_CM3 = 0d0
         ENDIF

      ELSE IF ( Input_Opt%ITS_AN_AEROSOL_SIM ) THEN

         !---------------------
         ! Offline simulation
         !---------------------

         ! 1-D grid box index for SUNCOS
         JLOOP = ( (J-1) * IIPAR ) + I

         ! Test for sunlight...
         IF ( State_Met%SUNCOS(I,J) > 0d0 .and. TCOSZ(I,J) > 0d0 ) THEN

            ! Impose a diurnal variation on OH during the day
            OH_MOLEC_CM3 = OH(I,J,L)                              *           
     &                     ( State_Met%SUNCOS(I,J) / TCOSZ(I,J) ) *
     &                     ( 1440d0                / GET_TS_CHEM() )

            ! Make sure OH is not negative
            OH_MOLEC_CM3 = MAX( OH_MOLEC_CM3, 0d0 )
               
         ELSE

            ! At night, OH goes to zero
            OH_MOLEC_CM3 = 0d0

         ENDIF

      ELSE

         !---------------------
         ! Invalid simulation
         !---------------------         
         CALL ERROR_STOP( 'Invalid NSRCX!', 'GET_OH (sulfate_mod.f)')

      ENDIF

      END FUNCTION GET_OH
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: set_oh
!
! !DESCRIPTION: Subroutine SET\_OH saves the modified OH value back to
!  SMVGEAR's CSPEC array for coupled sulfate/aerosol simulations. (bmy, 
!  12/16/02)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE SET_OH( I, J, L, OH ) 
!
! !USES:
!
      USE CMN_SIZE_MOD
      USE COMODE_MOD,   ONLY : CSPEC, JLOP
      USE TRACERID_MOD, ONLY : IDOH
!
! !INPUT PARAMETERS: 
!
      INTEGER, INTENT(IN) :: I, J, L ! Indices for lon, lat, vertical level
      REAL*8,  INTENT(IN) :: OH      ! OH at (I,J,L) to be saved into CSPEC
!
! !REVISION HISTORY:
!  (1 ) We assume SETTRACE has been called to define IDOH (bmy, 12/16/02)
!  22 Dec 2011 - M. Payer    - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER             :: JLOOP

      !=================================================================
      ! SET_OH begins here!
      !=================================================================

      ! JLOOP = SMVGEAR 1-D grid box index
      JLOOP = JLOP(I,J,L) 

      ! Replace OH into CSPEC(troposphere only)
      IF ( JLOOP > 0 ) THEN
         CSPEC(JLOOP,IDOH) = OH
      ENDIF

      END SUBROUTINE SET_OH
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_no3
!
! !DESCRIPTION: Function GET\_NO3 returns NO3 from SMVGEAR's CSPEC array (for
!  coupled runs) or monthly mean OH (for offline runs).  For offline runs, the
!  concentration of NO3 is set to zero during the day. (rjp, bmy, 12/16/02)
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_NO3( I, J, L, Input_Opt, State_Met ) 
     &         RESULT( NO3_MOLEC_CM3 ) 
!
! !USES:
!

      USE CMN_SIZE_MOD
      USE COMODE_MOD,         ONLY : CSPEC, JLOP
      USE ERROR_MOD,          ONLY : ERROR_STOP
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE GLOBAL_NO3_MOD,     ONLY : NO3
      USE TRACERID_MOD,       ONLY : IDNO3
!
! !INPUT PARAMETERS: 
!
      INTEGER,        INTENT(IN)  :: I, J, L     ! Lon, lat, vertical level
      TYPE(OptInput), INTENT(IN)  :: Input_Opt   ! Input Options object
      TYPE(MetState), INTENT(IN)  :: State_Met   ! Meteorology State object
!
! !REVISION HISTORY:
!  (1 ) Now references ERROR_STOP from "error_mod.f".  We also assume that
!        SETTRACE has been called to define IDNO3.  Now also set NO3 to 
!        zero during the day. (rjp, bmy, 12/16/02)
!  (2 ) Now reference ITS_A_FULLCHEM_SIM and ITS_AN_AEROSOL_SIM from 
!        "tracer_mod.f".  Also remove reference to CMN.   Also replace
!        IJSURF with an analytic function. (bmy, 7/20/04)
!  22 Dec 2011 - M. Payer    - Added ProTeX headers
!  09 Nov 2012 - M. Payer    - Replaced all met field arrays with State_Met
!                              derived type object
!  28 Nov 2012 - R. Yantosca - Replace SUNCOS with State_Met%SUNCOS
!  04 Mar 2013 - R. Yantosca - Now pass Input_Opt%ITS_A_FULLCHEM_SIM and
!                              Input_Opt%ITS_AN_AEROSOL_SIM
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER             :: JLOOP
      REAL*8              :: NO3_MOLEC_CM3
!
! !EXTERNAL FUNCTIONS:   
!
      REAL*8,  EXTERNAL   :: BOXVL

      !=================================================================
      ! GET_NO3 begins here!
      !=================================================================

      IF ( Input_Opt%ITS_A_FULLCHEM_SIM ) THEN

         !--------------------
         ! Coupled simulation
         !--------------------
            
         ! 1-D SMVGEAR grid box index
         JLOOP = JLOP(I,J,L)

         ! Take NO3 from the SMVGEAR array CSPEC
         ! NO3 is defined only in the troposphere
         IF ( JLOOP > 0 ) THEN
            NO3_MOLEC_CM3 = CSPEC(JLOOP,IDNO3)
         ELSE
            NO3_MOLEC_CM3 = 0d0
         ENDIF

      ELSE IF ( Input_Opt%ITS_AN_AEROSOL_SIM ) THEN

         !==============================================================  
         ! Offline simulation: Read monthly mean GEOS-CHEM NO3 fields
         ! in [v/v].  Convert these to [molec/cm3] as follows:
         !
         !  vol NO3   moles NO3    kg air     kg NO3/mole NO3
         !  ------- = --------- * -------- * ---------------- =  kg NO3 
         !  vol air   moles air      1        kg air/mole air
         !
         ! And then we convert [kg NO3] to [molec NO3/cm3] by:
         !  
         !  kg NO3   molec NO3   mole NO3     1     molec NO3
         !  ------ * --------- * -------- * ----- = --------- 
         !     1     mole NO3     kg NO3     cm3       cm3
         !          ^                    ^
         !          |____________________|  
         !            this is XNUMOL_NO3
         !
         ! If at nighttime, use the monthly mean NO3 concentration from
         ! the NO3 array of "global_no3_mod.f".  If during the daytime,
         ! set the NO3 concentration to zero.  We don't have to relax to 
         ! the monthly mean concentration every 3 hours (as for HNO3) 
         ! since NO3 has a very short lifetime. (rjp, bmy, 12/16/02) 
         !==============================================================

         ! 1-D grid box index for SUNCOS
         JLOOP = ( (J-1) * IIPAR ) + I

         ! Test if daylight
         IF ( State_Met%SUNCOS(I,J) > 0d0 ) THEN

            ! NO3 goes to zero during the day
            NO3_MOLEC_CM3 = 0d0
              
         ELSE

            ! At night: Get NO3 [v/v] and convert it to [kg]
            NO3_MOLEC_CM3 = NO3(I,J,L)     * State_Met%AD(I,J,L) *
     &                      ( 62d0/28.97d0 ) 

            ! Convert NO3 from [kg] to [molec/cm3]
            NO3_MOLEC_CM3 = NO3_MOLEC_CM3  * XNUMOL_NO3 /
     &                      BOXVL(I,J,L,State_Met)
               
         ENDIF
            
         ! Make sure NO3 is not negative
         NO3_MOLEC_CM3  = MAX( NO3_MOLEC_CM3, 0d0 )

      ELSE

         !--------------------
         ! Invalid simulation
         !--------------------
         CALL ERROR_STOP( 'Invalid NSRCX!','GET_NO3 (sulfate_mod.f)')

      ENDIF

      END FUNCTION GET_NO3
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: set_no3
!
! !DESCRIPTION: Subroutine SET\_NO3 saves the modified NO3 value back to 
!  SMVGEAR's CSPEC array for coupled lfate/aerosol simulations. (rjp, bmy, 
!  12/16/02, 7/20/04)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE SET_NO3( I, J, L, NO3 ) 
!
! !USES:
!
      USE CMN_SIZE_MOD
      USE COMODE_MOD,   ONLY : CSPEC, JLOP
      USE TRACERID_MOD, ONLY : IDNO3
!
! !INPUT PARAMETERS: 
!
      INTEGER, INTENT(IN) :: I, J, L ! Indices for lon, lat, vertical level
      REAL*8,  INTENT(IN) :: NO3     ! NO3 at (I,J,L) to be saved into CSPEC
!
! !REVISION HISTORY:
!  (1 ) We assume SETTRACE has been called to define IDNO3. (bmy, 12/16/02)
!  (2 ) Remove references to "error_mod.f" and CMN (bmy, 7/20/04)
!  22 Dec 2011 - M. Payer    - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER             :: JLOOP

      !=================================================================
      ! SET_NO3 begins here!
      !=================================================================

      ! 1-D grid box index for CSPEC
      JLOOP = JLOP(I,J,L) 

      ! Replace OH into CSPEC (troposphere only)
      IF ( JLOOP > 0 ) THEN
         CSPEC(JLOOP,IDNO3) = NO3
      ENDIF

      END SUBROUTINE SET_NO3
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_o3
!
! !DESCRIPTION: Function GET\_O3 returns monthly mean O3 for offline sulfate
!  aerosol simulations. (bmy, 12/16/02, 10/25/05)
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_O3( I, J, L, Input_Opt, State_Met ) RESULT( O3_VV )
!
! !USES:
!
      USE CMN_SIZE_MOD
      USE COMODE_MOD,         ONLY : CSPEC, JLOP, VOLUME
      USE ERROR_MOD,          ONLY : ERROR_STOP
      USE GIGC_Input_Opt_MOd, ONLY : OptInput
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE TRACER_MOD,         ONLY : XNUMOLAIR
      USE TRACERID_MOD,       ONLY : IDO3
!
! !INPUT PARAMETERS: 
!
      INTEGER,        INTENT(IN)  :: I, J, L     ! Lon, lat, vertical level
      TYPE(OptInput), INTENT(IN)  :: Input_Opt  
      TYPE(MetState), INTENT(IN)  :: State_Met   ! Meteorology State object
!
! !REVISION HISTORY: 
!  (1 ) We assume SETTRACE has been called to define IDO3. (bmy, 12/16/02)
!  (2 ) Now reference inquiry functions from "tracer_mod.f" (bmy, 7/20/04)
!  (3 ) Now remove reference to CMN, it's obsolete.  (bmy, 8/22/05)
!  (4 ) Now references XNUMOLAIR from "tracer_mod.f" (bmy, 10/25/05)
!  22 Dec 2011 - M. Payer    - Added ProTeX headers
!  09 Nov 2012 - M. Payer    - Replaced all met field arrays with State_Met
!                              derived type object
!  04 Mar 2013 - R. Yantosca - Now pass Input_Opt%ITS_A_FULLCHEM_SIM and
!                              Input_Opt%ITS_AN_AEROSOL_SIM
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER               :: JLOOP
      REAL*8                :: O3_VV
 
      !=================================================================
      ! GET_O3 begins here!
      !=================================================================

      IF ( Input_Opt%ITS_A_FULLCHEM_SIM ) THEN

         !--------------------
         ! Coupled simulation
         !--------------------

         ! JLOOP = SMVGEAR 1-D grid box index
         JLOOP = JLOP(I,J,L)

         ! Get O3 from CSPEC [molec/cm3] and convert it to [v/v]
         ! O3 data will only be defined below the tropopause
         IF ( JLOOP  > 0 ) THEN
            O3_VV = ( CSPEC(JLOOP,IDO3)       * 1d6       ) / 
     &              ( State_Met%AIRDEN(L,I,J) * XNUMOLAIR )
         ELSE
            O3_VV = 0d0
         ENDIF

      ELSE IF ( Input_Opt%ITS_AN_AEROSOL_SIM ) THEN
         
         !--------------------
         ! Offline simulation
         !--------------------

         ! Get O3 [v/v] for this gridbox & month
         ! O3 data will only be defined below the tropopause
         IF ( L <= LLTROP ) THEN
            O3_VV = O3m(I,J,L)
         ELSE
            O3_VV = 0d0
         ENDIF

      ELSE

         !--------------------
         ! Invalid simulation
         !--------------------         
         CALL ERROR_STOP( 'Invalid NSRCX!', 'GET_O3 (sulfate_mod.f)')

      ENDIF

      END FUNCTION GET_O3
!EOC
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!%%% NOTE: This is for the SEAC4RS volcano emissions.  Leave code commented
!!%%% out here so that we can restore it later. (bmy, 9/25/13)
!!------------------------------------------------------------------------------
!!                  GEOS-Chem Global Chemical Transport Model                  !
!!------------------------------------------------------------------------------
!!BOP
!!
!! !IROUTINE: read_volc
!!
!! !DESCRIPTION: Subroutine READ\_VOLC reads SO2 emissions from
!!  volcanoes. (bmy, jaf, skim, 1/11/13)
!!\\
!!\\
!! !INTERFACE:
!!
!      SUBROUTINE READ_VOLC( INMONTH, INYEAR )
!!
!! !USES:
!!
!      ! Modules for netCDF read
!      USE m_netcdf_io_open
!      USE m_netcdf_io_get_dimlen
!      USE m_netcdf_io_read
!      USE m_netcdf_io_readattr
!      USE m_netcdf_io_close
!
!      ! GEOS-Chem modules
!      USE ERROR_MOD,     ONLY : ALLOC_ERR
!      USE DIRECTORY_MOD, ONLY : DATA_DIR_1x1
!      USE GRID_MOD,      ONLY : GET_XOFFSET, GET_YOFFSET
!      USE TIME_MOD,      ONLY : EXPAND_DATE
!
!      USE CMN_SIZE_MOD         ! Size parameters
!
!      IMPLICIT NONE
!
!#     include "netcdf.inc"
!!
!! !INPUT PARAMETERS: 
!!
!      INTEGER, INTENT(IN)     :: INMONTH ! Current month number (1-12)
!      INTEGER, INTENT(IN)     :: INYEAR  ! Current 4-digit year
!!
!! !REMARKS:
!!  1. Assumes that you have a netCDF library (either v3 or v4) installed on 
!!  your system.
!!  2. Uses most recent year of non-eruptive emissions. Eruptive emissions 
!!  are only used if simulation year is within year of available data 
!!
!! !REVISION HISTORY: 
!!  16 Jun 2012 - J. Fisher   - Initial version (replaces READ_ERUP_VOLC
!!                              and READ_NONERUP_VOLC)
!!  11 Jan 2013 - S. Kim      - Remove unneeded for loops
!!                              Minor bug fixes for reading netCDF files
!!EOP
!!------------------------------------------------------------------------------
!!BOC
!!
!! !LOCAL VARIABLES:
!!
!      !=================================================================
!      ! Variable declarations
!      !=================================================================
!
!      ! Data arrays
!      INTEGER              :: NEVENTS
!      INTEGER, ALLOCATABLE :: v_date(:)
!      INTEGER, ALLOCATABLE :: v_elevation(:)
!      INTEGER, ALLOCATABLE :: v_cloud_column_height(:)
!      REAL*4,  ALLOCATABLE :: v_lon(:),v_lat(:)
!      REAL*4,  ALLOCATABLE :: v_so2(:)
!
!      ! Character strings
!      CHARACTER(LEN=255) :: nc_dir             ! netCDF directory name
!      CHARACTER(LEN=255) :: nc_file            ! netCDF file name
!      CHARACTER(LEN=255) :: v_name             ! netCDF variable name 
!      CHARACTER(LEN=255) :: a_name             ! netCDF attribute name
!      CHARACTER(LEN=255) :: a_val              ! netCDF attribute value
!
!      ! Arrays for netCDF start and count values
!      INTEGER            :: st1d(1), ct1d(1)   ! For 1D arrays    
!      INTEGER            :: st2d(2), ct2d(2)   ! For 2D arrays 
!      INTEGER            :: st3d(3), ct3d(3)   ! For 3D arrays 
!      INTEGER            :: st4d(4), ct4d(4)   ! For 4D arrays 
!      INTEGER            :: st5d(5), ct5d(5)   ! For 5D arrays 
!      INTEGER            :: st6d(6), ct6d(6)   ! For 6D arrays 
!
!      ! Variables for reading, etc.
!      INTEGER            :: YYYYMMDD
!      INTEGER            :: I0, J0
!      INTEGER            :: AS, fId, N, n_nv, n_ev, I, J
!      REAL*8             :: XEDGE(IIPAR+1), YEDGE(JJPAR+1)
!
!      !=================================================================
!      ! Open and read data from the netCDF file
!      !=================================================================
!
!      ! Set reading date
!      YYYYMMDD = 10000 * InYear + 100 * InMonth + 1
!
!      ! Replace time & date tokens in the directory name
!      nc_dir  = TRIM( DATA_DIR_1x1 ) // 'volcano_SO2_201206/YYYY'
!      CALL Expand_Date( nc_dir, YYYYMMDD, 000000 )
!
!      ! Replace time & date tokens in the file name
!      nc_file = '/Volcano_SO2_YYYYMM.nc'
!      CALL Expand_Date( nc_file, YYYYMMDD, 000000 )
!
!      ! Construct complete file path
!      nc_file = TRIM( nc_dir ) // TRIM( nc_file )
!
!      ! Echo output
!      WRITE( 6, 150 ) TRIM( nc_file )
! 150  FORMAT( '     - READ_VOLC: Reading ', a )
!
!      ! Open netCDF file
!      CALL NcOp_Rd( fId, TRIM( nc_file ) )
!
!      ! Read the dimensions from the netCDF file
!      CALL NcGet_DimLen( fId, 'nevents', NEVENTS )
!
!      ! Allocate arrays for reading
!      ALLOCATE( v_date ( NEVENTS ), STAT=AS )
!      IF ( AS /= 0 ) CALL ALLOC_ERR( 'v_date' )
!      v_date = 0d0
!      
!      ALLOCATE( v_elevation ( NEVENTS ), STAT=AS )
!      IF ( AS /= 0 ) CALL ALLOC_ERR( 'v_elevation' )
!      v_elevation = 0d0
!      
!      ALLOCATE( v_cloud_column_height ( NEVENTS ), STAT=AS )
!      IF ( AS /= 0 ) CALL ALLOC_ERR( 'v_cloud_column_height' )
!      v_cloud_column_height = 0d0
!      
!      ALLOCATE( v_lon ( NEVENTS ), STAT=AS )
!      IF ( AS /= 0 ) CALL ALLOC_ERR( 'v_lon' )
!      v_lon = 0d0
!      
!      ALLOCATE( v_lat ( NEVENTS ), STAT=AS )
!      IF ( AS /= 0 ) CALL ALLOC_ERR( 'v_lat' )
!      v_lat = 0d0
!      
!      ALLOCATE( v_so2 ( NEVENTS ), STAT=AS )
!      IF ( AS /= 0 ) CALL ALLOC_ERR( 'v_so2' )
!      v_so2 = 0d0
!      
!      !----------------------------------------
!      ! VARIABLE: date
!      !----------------------------------------
!
!      ! Variable name
!      v_name = "date"
!
!      ! Read date from file
!      st1d   = (/ 1 /)
!      ct1d   = (/ NEVENTS /)
!      CALL NcRd( v_date, fId, TRIM(v_name), st1d, ct1d )
!
!      ! Read the date:long_name attribute
!      a_name = "long_name"
!      CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )
!
!      ! Echo info to stdout
!      WRITE( 6, 130 ) TRIM(v_name), TRIM(a_val)
!
!      !----------------------------------------
!      ! VARIABLE: elevation
!      !----------------------------------------
!
!      ! Variable name
!      v_name = "elevation"
!
!      ! Read elevation from file
!      st1d   = (/ 1 /)
!      ct1d   = (/ NEVENTS /)
!      CALL NcRd( v_elevation, fId, TRIM(v_name), st1d, ct1d )
!
!      ! Read the elevation:long_name attribute
!      a_name = "long_name"
!      CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )
!
!      ! Read the elevation:units attribute
!      a_name = "units"
!      CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )
!
!      ! Echo info to stdout
!      WRITE( 6, 130 ) TRIM(v_name), TRIM(a_val)
!
!      !----------------------------------------
!      ! VARIABLE: cloud_column_height
!      !----------------------------------------
!
!      ! Variable name
!      v_name = "cloud_column_height"
!
!      ! Read cloud_column_height from file
!      st1d   = (/ 1 /)
!      ct1d   = (/ NEVENTS /)
!      CALL NcRd( v_cloud_column_height, fId, TRIM(v_name), st1d, ct1d )
!
!      ! Read the cloud_column_height:long_name attribute
!      a_name = "long_name"
!      CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )
!
!      ! Read the cloud_column_height:units attribute
!      a_name = "units"
!      CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )
!
!      ! Echo info to stdout
!      WRITE( 6, 130 ) TRIM(v_name), TRIM(a_val)
!
!      !----------------------------------------
!      ! VARIABLE: lon
!      !----------------------------------------
!
!      ! Variable name
!      v_name = "lon"
!
!      ! Read lon from file
!      st1d   = (/ 1 /)
!      ct1d   = (/ NEVENTS /)
!      CALL NcRd( v_lon, fId, TRIM(v_name), st1d, ct1d )
!
!      ! Read the lon:long_name attribute
!      a_name = "long_name"
!      CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )
!
!      ! Read the lon:units attribute
!      a_name = "units"
!      CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )
!
!      ! Echo info to stdout
!      WRITE( 6, 130 ) TRIM(v_name), TRIM(a_val)
!
!      !----------------------------------------
!      ! VARIABLE: lat
!      !----------------------------------------
!
!      ! Variable name
!      v_name = "lat"
!
!      ! Read lat from file
!      st1d   = (/ 1 /)
!      ct1d   = (/ NEVENTS /)
!      CALL NcRd( v_lat, fId, TRIM(v_name), st1d, ct1d )
!
!      ! Read the lat:long_name attribute
!      a_name = "long_name"
!      CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )
!
!      ! Read the lat:units attribute
!      a_name = "units"
!      CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )
!
!      ! Echo info to stdout
!      WRITE( 6, 130 ) TRIM(v_name), TRIM(a_val)
!
!      !----------------------------------------
!      ! VARIABLE: so2
!      !----------------------------------------
!
!      ! Variable name
!      v_name = "so2"
!
!      ! Read so2 from file
!      st1d   = (/ 1 /)
!      ct1d   = (/ NEVENTS /)
!      CALL NcRd( v_so2, fId, TRIM(v_name), st1d, ct1d )
!
!      ! Read the so2:long_name attribute
!      a_name = "long_name"
!      CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )
!
!      ! Read the so2:units attribute
!      a_name = "units"
!      CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )
!
!      ! Echo info to stdout
!      WRITE( 6, 130 ) TRIM(v_name), TRIM(a_val)
!
!      ! Close netCDF file
!      CALL NcCl( fId )
!
!      ! Echo info to stdout
!      WRITE( 6, 140 )
!      WRITE( 6, 100 ) REPEAT( '%', 79 )
!
!      ! FORMAT statements
! 100  FORMAT( a                                              )
! 130  FORMAT( '%% Successfully read ',       a, ' [', a, ']' )
! 140  FORMAT( '%% Successfully closed file!'                 )
!
!      !=================================================================
!      ! Separate eruptive and non-eruptive events
!      !=================================================================
!
!      ! Count event types
!      N_NONERUP = 0
!      N_ERUP = 0
!      DO N = 1, NEVENTS
!         IF ( v_elevation(N) == v_cloud_column_height(N) ) THEN
!            N_NONERUP = N_NONERUP + 1
!         ELSE
!            N_ERUP = N_ERUP + 1
!         ENDIF
!      ENDDO
!
!      ! If already allocated from a previous month
!      !    deallocate to prevent allocation error
!      !    (skim, 1/14/13)
!      IF ( ALLOCATED( VEV_ELEV     ) ) DEALLOCATE( VEV_ELEV     )
!      IF ( ALLOCATED( VEV_SO2      ) ) DEALLOCATE( VEV_SO2      )
!      IF ( ALLOCATED( VEV_DATE     ) ) DEALLOCATE( VEV_DATE     )
!      IF ( ALLOCATED( VEV_II       ) ) DEALLOCATE( VEV_II       )
!      IF ( ALLOCATED( VEV_JJ       ) ) DEALLOCATE( VEV_JJ       )
!      IF ( ALLOCATED( VEV_CLOUDTOP ) ) DEALLOCATE( VEV_CLOUDTOP )
!      IF ( ALLOCATED( VNV_ELEV     ) ) DEALLOCATE( VNV_ELEV     )
!      IF ( ALLOCATED( VNV_SO2      ) ) DEALLOCATE( VNV_SO2      )
!      IF ( ALLOCATED( VNV_DATE     ) ) DEALLOCATE( VNV_DATE     )
!      IF ( ALLOCATED( VNV_II       ) ) DEALLOCATE( VNV_II       )
!      IF ( ALLOCATED( VNV_JJ       ) ) DEALLOCATE( VNV_JJ       )
!
!      ! Allocate arrays for eruptive events
!      ALLOCATE( vev_date ( N_ERUP ), STAT=AS )
!      IF ( AS /= 0 ) CALL ALLOC_ERR( 'vev_date' )
!      vev_date = 0d0
!      ALLOCATE( vev_ii ( N_ERUP ), STAT=AS )
!      IF ( AS /= 0 ) CALL ALLOC_ERR( 'vev_ii' )
!      vev_ii = 0
!      ALLOCATE( vev_jj ( N_ERUP ), STAT=AS )
!      IF ( AS /= 0 ) CALL ALLOC_ERR( 'vev_jj' )
!      vev_jj = 0
!      ALLOCATE( vev_so2 ( N_ERUP ), STAT=AS )
!      IF ( AS /= 0 ) CALL ALLOC_ERR( 'vev_so2' )
!      vev_so2 = 0d0
!      ALLOCATE( vev_elev ( N_ERUP ), STAT=AS )
!      IF ( AS /= 0 ) CALL ALLOC_ERR( 'vev_elev' )
!      vev_elev = 0d0
!      ALLOCATE( vev_cloudtop ( N_ERUP ), STAT=AS )
!      IF ( AS /= 0 ) CALL ALLOC_ERR( 'vev_cloudtop' )
!      vev_cloudtop = 0d0
!      
!      ! Allocate arrays for non-eruptive events
!      ALLOCATE( vnv_date ( N_NONERUP ), STAT=AS )
!      IF ( AS /= 0 ) CALL ALLOC_ERR( 'vnv_date' )
!      vnv_date = 0d0
!      ALLOCATE( vnv_ii ( N_NONERUP ), STAT=AS )
!      IF ( AS /= 0 ) CALL ALLOC_ERR( 'vnv_ii' )
!      vnv_ii = 0
!      ALLOCATE( vnv_jj ( N_NONERUP ), STAT=AS )
!      IF ( AS /= 0 ) CALL ALLOC_ERR( 'vnv_jj' )
!      vnv_jj = 0
!      ALLOCATE( vnv_so2 ( N_NONERUP ), STAT=AS )
!      IF ( AS /= 0 ) CALL ALLOC_ERR( 'vnv_so2' )
!      vnv_so2 = 0d0
!      ALLOCATE( vnv_elev ( N_NONERUP ), STAT=AS )
!      IF ( AS /= 0 ) CALL ALLOC_ERR( 'vnv_elev' )
!      vnv_elev = 0d0
!
!      I0 = GET_XOFFSET( GLOBAL=.TRUE. )
!      J0 = GET_YOFFSET( GLOBAL=.TRUE. )
!
!      ! Re-assign values to eruptive or non-eruptive arrays
!      ! Change the counters to start at zero (skim, 8/17/13)
!      n_nv = 1
!      n_ev = 1
!      DO N = 1, NEVENTS
!
!         ! Determine which grid box to place the emissions
!         I = INT( ( v_lon(N) + 180d0 - (I0 * DISIZE) )/DISIZE + 1.5d0 )
!         J = INT( ( v_lat(N) +  90d0 - (J0 * DJSIZE) )/DJSIZE + 1.5d0 )
!
!         ! Only check for the date line for global simulations 
!         ! (skim, 8/17/13)
!#if defined( GRID4x5 ) || defined( GRID2x25 ) 
!         ! Handle the date line correctly
!         IF ( I > IIPAR ) I = I - IIPAR
!#endif
!
!         ! Non-eruptive
!         IF ( v_elevation(N) == v_cloud_column_height(N) ) THEN
!            vnv_date(n_nv) = v_date(n)
!            vnv_elev(n_nv) = v_elevation(n)
!            vnv_so2(n_nv)  = v_so2(n)
!            vnv_ii(n_nv)   = I
!            vnv_jj(n_nv)   = J
!            ! Augment non-eruptive counter
!            n_nv = n_nv + 1
!         ELSE
!         ! Eruptive
!         IF ( DO_ERUPTIVE ) THEN
!            vev_date(n_ev)     = v_date(n)
!            vev_elev(n_ev)     = v_elevation(n)
!            vev_so2(n_ev)      = v_so2(n)
!            vev_cloudtop(n_ev) = v_cloud_column_height(n)
!            vev_ii(n_ev)       = I
!            vev_jj(n_ev)       = J
!            ! Augment eruptive counter
!            n_ev = n_ev + 1
!         ENDIF
!         ENDIF
!      ENDDO
!
!      !=================================================================
!      ! Convert units (kt/event/day --> kg/event/s)
!      !=================================================================
!      IF ( DO_ERUPTIVE ) vev_so2 = vev_so2 * 1d6 / ( 24d0 * 3600d0 )
!      vnv_so2 = vnv_so2 * 1d6 / ( 24d0 * 3600d0 )
!
!      END SUBROUTINE READ_VOLC
!!EOC
!!------------------------------------------------------------------------------
!!                  GEOS-Chem Global Chemical Transport Model                  !
!!------------------------------------------------------------------------------
!!BOP
!!
!! !IROUTINE: grid_volc
!!
!! !DESCRIPTION: Subroutine GRID\_VOLC grids SO2 emissions from
!!  volcanoes. (jaf, 6/17/12)
!!\\
!!\\
!! !INTERFACE:
!!
!      SUBROUTINE GRID_VOLC( INDAY, INMONTH, INYEAR, State_Met )
!!
!! !USES:
!!
!      USE ERROR_MOD,          ONLY : GEOS_CHEM_STOP, IS_SAFE_DIV
!      USE GIGC_State_Met_Mod, ONLY : MetState
!
!      USE CMN_SIZE_MOD         ! Size parameters
!
!      IMPLICIT NONE
!
!!
!! !INPUT PARAMETERS: 
!!
!      INTEGER, INTENT(IN)        :: INDAY     ! Current day number (1-31)
!      INTEGER, INTENT(IN)        :: INMONTH   ! Current month number (1-12)
!      INTEGER, INTENT(IN)        :: INYEAR    ! Current 4-digit year
!      TYPE(MetState), INTENT(IN) :: State_Met ! Meteorology state object
!!
!! !REVISION HISTORY: 
!!  17 Jun 2012 - J. Fisher   - Initial version (replaces READ_ERUP_VOLC
!!                              and READ_NONERUP_VOLC)
!!  11 Jan 2013 - S. Kim      - Update to use GIGC conventions
!!  26 Jul 2013 - R. Yantosca - Add trap for OOB errors
!!EOP
!!------------------------------------------------------------------------------
!!BOC
!!
!! !LOCAL VARIABLES:
!!
!      INTEGER            :: YYYYMMDD, N, L, L_BOT, L_TOP
!      REAL*8             :: BOX_BOT, BOX_TOP
!      REAL*8             :: PLUME_BOT, PLUME_TOP, DTHICK(LLPAR)
!      REAL*8             :: CHECKSUM_NV, CHECKSUM_EV, T_DIFF
!
!      !=================================================================
!      ! Set-up
!      !=================================================================
!
!      ! Reset emissions to zero each day
!      ESO2_nv = 0d0
!      ESO2_ev = 0d0
!      CHECKSUM_NV = 0d0
!      CHECKSUM_EV = 0d0
!
!      ! Set date
!      YYYYMMDD = 10000 * InYear + 100 * InMonth + InDay
!
!      !=================================================================
!      ! Loop over non-eruptive events
!      !=================================================================
!
!      DO N = 1, N_NONERUP
!         IF ( vnv_date(N) < YYYYMMDD ) CYCLE
!         IF ( vnv_date(N) > YYYYMMDD ) EXIT
!
!         ! Extra check to avoid OOB errors (bmy, 7/26/13)
!         IF ( vnv_ii(N) < 1 .or. vnv_ii(N) > IIPAR ) CYCLE
!         IF ( vnv_jj(N) < 1 .or. vnv_jj(N) > JJPAR ) CYCLE
!
!         ! For consistency check
!         CHECKSUM_NV = CHECKSUM_NV + vnv_so2(N)
!
!         ! Loops over levels (non-eruptive emissions at crater elevation)
!         BOX_BOT = 0d0
!         DO L = 1, LLPAR
!            BOX_TOP = State_Met%BXHEIGHT(vnv_ii(N),vnv_jj(N),L) 
!     &                  + BOX_BOT
!            IF ( ( vnv_elev(N) >= BOX_BOT ) .and. 
!     &           ( vnv_elev(N)  < BOX_TOP ) ) THEN
!               ESO2_nv(vnv_ii(N),vnv_jj(N),L) = 
!     &            ESO2_nv(vnv_ii(N),vnv_jj(N),L) + vnv_so2(N)
!               EXIT
!            ENDIF
!            BOX_BOT = BOX_TOP
!         ENDDO 
!         ! Above model top layer
!         IF ( vnv_elev(N) >= BOX_TOP ) 
!     &      ESO2_nv(vnv_ii(N),vnv_jj(N),LLPAR) = 
!     &         ESO2_nv(vnv_ii(N),vnv_jj(N),LLPAR) + vnv_so2(N)
!
!      ENDDO
!
!      !=================================================================
!      ! Loop over eruptive events
!      !=================================================================
!
!      IF ( DO_ERUPTIVE ) THEN 
!      DO N = 1, N_ERUP
!         IF ( vev_date(N) < YYYYMMDD ) CYCLE
!         IF ( vev_date(N) > YYYYMMDD ) EXIT
!
!         ! Extra check to avoid OOB erros (bmy, 7/26/13)
!         IF ( vev_ii(N) < 1 .or. vev_ii(N) > IIPAR ) CYCLE
!         IF ( vev_jj(N) < 1 .or. vev_jj(N) > JJPAR ) CYCLE
!
!         ! For consistency check
!         CHECKSUM_EV = CHECKSUM_EV + vev_so2(N)
!
!         ! Define plume edges (eruptive emissions in top 1/3 of cloud)
!         PLUME_TOP = vev_cloudtop(N)
!         PLUME_BOT = PLUME_TOP - (PLUME_TOP - vev_elev(N))/3d0
!
!         ! Loop over levels
!         L_BOT = 0
!         L_TOP = 0
!         BOX_BOT = 0d0
!         DO L = 1, LLPAR
!
!            ! New box top
!            BOX_TOP = State_Met%BXHEIGHT(vev_ii(N),vev_jj(N),L) 
!     &                   + BOX_BOT
!
!            ! Default thickness (used for plumes that extend over
!            ! multiple levels)
!            DTHICK(L) = State_Met%BXHEIGHT(vev_ii(N),vev_jj(N),L)
!
!            ! Case 1: If top edge of level < bottom of plume, go to next
!            ! level
!            IF ( BOX_TOP <= PLUME_BOT ) THEN
!               BOX_BOT = BOX_TOP
!               CYCLE
!            ENDIF
!
!            ! Case 2: If the plume is only in the current level
!            IF ( BOX_BOT <= PLUME_BOT .and. BOX_TOP >= PLUME_TOP ) THEN
!               L_BOT = L
!               L_TOP = L
!               DTHICK(L) = PLUME_TOP - PLUME_BOT ! Entire plume
!               EXIT
!            ENDIF
!
!            ! Case 3: The plume extends more than one level. Find the
!            ! lowest (L_BOT) highest (L_TOP), and the thickness (DTHICK)
!            ! in each level.
!            IF (BOX_BOT <= PLUME_BOT .and. BOX_TOP > PLUME_BOT) THEN
!               L_BOT = L
!               DTHICK(L) = BOX_TOP - PLUME_BOT
!            ELSE IF (BOX_BOT<PLUME_TOP .and. BOX_TOP>=PLUME_TOP) THEN
!               L_TOP = L
!               DTHICK(L) = PLUME_TOP - BOX_BOT
!               EXIT
!            ENDIF
!
!            ! New box bottom
!            BOX_BOT = BOX_TOP
!
!         ENDDO
!
!         ! Above top model layer
!         IF ( L_TOP == 0 ) L_TOP = LLPAR
!
!         ! Inject SO2 emissions into correct levels, distributed evenly
!         ! with altitude across the plume (dthick/plume_thickness)
!         DO L = L_BOT, L_TOP
!            ESO2_ev(vev_ii(N),vev_jj(N),L) = 
!     &         ESO2_ev(vev_ii(N),vev_jj(N),L) + 
!     &         vev_so2(N) * DTHICK(L) / (PLUME_TOP-PLUME_BOT)
!         ENDDO
!
!      ENDDO
!      ENDIF
!
!      !=================================================================
!      ! Consistency check
!      !=================================================================
!
!      ! Non-eruptive volcanoes
!      T_DIFF = ABS( SUM( ESO2_nv ) - CHECKSUM_NV )         
!      IF ( IS_SAFE_DIV( T_DIFF, CHECKSUM_NV ) ) THEN
!         IF ( T_DIFF / CHECKSUM_NV > 1d-5 ) THEN
!            PRINT*, 'Non-eruptive volcanic SO2 emissions before ' //
!     &              'and after regridding are not equivalent!'
!            CALL FLUSH(6)
!            CALL GEOS_CHEM_STOP
!         ENDIF
!      ENDIF
!      IF ( DO_ERUPTIVE ) THEN 
!      ! Non-eruptive volcanoes
!      T_DIFF = ABS( SUM( ESO2_ev ) - CHECKSUM_EV )         
!      IF ( IS_SAFE_DIV( T_DIFF, CHECKSUM_EV ) ) THEN
!         IF ( T_DIFF / CHECKSUM_EV > 1d-5 ) THEN
!            PRINT*, 'Eruptive volcanic SO2 emissions before ' //
!     &              'and after regridding are not equivalent!'
!            CALL FLUSH(6)
!            CALL GEOS_CHEM_STOP
!         ENDIF
!      ENDIF
!      ENDIF
!
!      END SUBROUTINE GRID_VOLC
!!EOC
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_nonerup_volc
!
! !DESCRIPTION: Subroutine READ\_NONERUP\_VOLC reads SO2 emissions from
!  non-eruptive volcanoes. (rjp, bdf, bmy, jaf, 9/19/02, 10/3/05, 10/15/09)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE READ_NONERUP_VOLC( INDAY, INMONTH, INYEAR )
!
! !USES:
!
      USE BPCH2_MOD,      ONLY : GET_TAU0, READ_BPCH2
      USE CMN_SIZE_MOD
      USE DIRECTORY_MOD,  ONLY : DATA_DIR_1x1
      USE REGRID_A2A_MOD, ONLY : DO_REGRID_A2A
      USE TIME_MOD,       ONLY : EXPAND_DATE
!
! !INPUT PARAMETERS: 
!
      INTEGER, INTENT(IN)     :: INDAY   ! Current 2-digit day 
      INTEGER, INTENT(IN)     :: INMONTH ! Current month number (1-12)
      INTEGER, INTENT(IN)     :: INYEAR  ! Current 4-digit year
!
! !REVISION HISTORY: 
!  (1 ) Split off from old module routine "sulfate_readyr" (bmy, 9/19/02)
!  (2 ) Now references DATA_DIR from "directory_mod.f" (bmy, 7/20/04)
!  (3 ) Now read files from "sulfate_sim_200508/" (bmy, 7/28/05)
!  (4 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (5 ) Complete re-write as volcanic emissions are now monthly and
!	stored as BPCH files (jaf, bmy, 10/15/09)
!  (6 ) Now use MNYEAR_VOLC and MXYEAR_VOLC as 1st and last year of emissions.
!       (ccc, 9/30/10)
!  (7 ) Volcanic data have been updated. Use a new directory. (ccc, 9/30/10)
!  22 Dec 2011 - M. Payer    - Added ProTeX headers
!  13 Mar 2012 - M. Cooper   - Changed regrid algorithm to map_a2a
!  24 May 2012 - R. Yantosca - Fixed minor bugs in map_a2a implementation
!  24 Aug 2012 - R. Yantosca - DO_REGRID_A2A now reads netCDF input file
!  03 Jan 2013 - M. Payer    - Renamed PERAREA to IS_MASS in DO_REGRID_A2A
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL*4               :: ARRAY_GEN(I1x1,J1x1-1,LVOLC)
      REAL*8, TARGET       :: ARRAY_GEN8(I1x1,J1x1-1)
      REAL*8, TARGET       :: ARRAY(IIPAR,JJPAR,LVOLC)
      REAL*8               :: XTAU
      INTEGER              :: THISYEAR, THISDAY, YYYYMMDD, L
      CHARACTER(LEN=255)   :: FILENAME
      CHARACTER(LEN=255)   :: LLFILENAME
      REAL*8, POINTER      :: OUTGRID(:,:) => NULL()
      REAL*8, POINTER      :: INGRID(:,:)  => NULL()

      !=================================================================
      ! READ_NONERUP_VOLC begins here!
      !=================================================================

      ! Set year based on availability of volcanic emission files
      THISYEAR = INYEAR
      THISYEAR = MAX( INYEAR, MNYEAR_VOLC )
      THISYEAR = MIN( INYEAR, MXYEAR_VOLC )

      ! Need to deal with leap years for which there is no data (i.e.
      ! 2008). Assume emissions on Feb. 29th are identical to emissions
      ! on Feb. 28th (jaf, 10/15/09)
      THISDAY = INDAY
      IF ( (INMONTH == 2) .and. (INDAY == 29) ) THISDAY = 28

      ! Set date
      YYYYMMDD = 10000 * THISYEAR + 100 * INMONTH + THISDAY

      ! File name
      FILENAME = TRIM( DATA_DIR_1x1 )  //
     &           'volcano_SO2_201010/' //
     &           'YYYY/SO2_volc.nonerup.YYYYMM.generic.1x1'

      ! Replace YYYY/MM in the file name
      CALL EXPAND_DATE( FILENAME, YYYYMMDD, 000000 )

      ! Echo output
      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( '     - READ_NONERUP_VOLC: Reading ', a )

      ! TAU value (use this year)
      XTAU = GET_TAU0( INMONTH, THISDAY, THISYEAR )

      ! Read non-eruptive emissions (kg/box/day)
      CALL READ_BPCH2( FILENAME, 'SO2-NV-$',        26,
     &                 XTAU,      I1x1,             J1x1-1,
     &                 LVOLC,     ARRAY_GEN(:,:,:), QUIET=.TRUE. )

      ! File with lat/lon edges for regridding
      LLFILENAME = TRIM( DATA_DIR_1x1 ) //
     &             'MAP_A2A_Regrid_201203/MAP_A2A_latlon_generic1x1.nc'

      ! Array regridding has to be done level by level since the array
      ! has a different vertical dimension than expected
      DO L = 1, LVOLC

         ! Cast from REAL*4 to REAL*8
         ARRAY_GEN8 = ARRAY_GEN(:,:,L)

         ! Point to input & output arrays
         INGRID  => ARRAY_GEN8                     !Units kg 
         OUTGRID => ARRAY(:,:,L)

         ! Regrid
         CALL DO_REGRID_A2A( LLFILENAME, I1x1,    J1x1-1,
     &                       INGRID,     OUTGRID, IS_MASS=1,
     &                       netCDF=.TRUE.                   )

         ! Free pointers
         NULLIFY( INGRID, OUTGRID )

      ENDDO

      ! Unit conversion: [kg SO2/box/day] -> [kg SO2/box/s]
      ENV = ARRAY / ( 24.d0 * 3600.d0 )

      END SUBROUTINE READ_NONERUP_VOLC
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_erup_volc
!
! !DESCRIPTION: Subroutine READ\_ERUP\_VOLC reads SO2 emissions from eruptive
!  volcanoes. (rjp, bdf, bmy, jaf, 9/19/02, 10/3/05, 10/15/09)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE READ_ERUP_VOLC( INDAY, INMONTH, INYEAR )
!
! !USES:
!
      USE BPCH2_MOD,      ONLY : GET_TAU0, READ_BPCH2
      USE CMN_SIZE_MOD
      USE DIRECTORY_MOD,  ONLY : DATA_DIR_1x1
      USE REGRID_A2A_MOD, ONLY : DO_REGRID_A2A
      USE TIME_MOD,       ONLY : EXPAND_DATE
!
! !INPUT PARAMETERS: 
!
      INTEGER, INTENT(IN)     :: INDAY   ! Current 2-digit day
      INTEGER, INTENT(IN)     :: INMONTH ! Current month number (1-12)
      INTEGER, INTENT(IN)     :: INYEAR  ! Current 4-digit year
!
! !REVISION HISTORY: 
!  (1 ) Split off from old module routine "sulfate_readyr" (bmy, 9/19/02)
!  (2 ) Now references DATA_DIR from "directory_mod.f" (bmy, 7/20/04)
!  (3 ) Now read files from "sulfate_sim_200508/" (bmy, 7/28/05)
!  (4 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (5 ) Complete re-write as volcanic emissions are now monthly and
!	stored as BPCH files (jaf, bmy, 10/15/09)
!  (6 ) Now use MNYEAR_VOLC and MXYEAR_VOLC as 1st and last year of emissions.
!       (ccc, 9/30/10)
!  (7 ) Volcanic data have been updated. Use a new directory. (ccc, 9/30/10)
!  22 Dec 2011 - M. Payer    - Added ProTeX headers
!  13 Mar 2012 - M. Cooper   - Changed regrid algorithm to map_a2a
!  24 May 2012 - R. Yantosca - Fixed minor bugs in map_a2a implementation
!  24 Aug 2012 - R. Yantosca - DO_REGRID_A2A now reads netCDF input file
!  03 Jan 2013 - M. Payer    - Renamed PERAREA to IS_MASS in DO_REGRID_A2A
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL*4               :: ARRAY_GEN(I1x1,J1x1-1,LVOLC)
      REAL*8, TARGET       :: ARRAY_GEN8(I1x1,J1x1-1)
      REAL*8, TARGET       :: ARRAY(IIPAR,JJPAR,LVOLC)
      REAL*8               :: XTAU
      INTEGER              :: YYYYMMDD, L
      CHARACTER(LEN=255)   :: FILENAME
      CHARACTER(LEN=255)   :: LLFILENAME
      REAL*8, POINTER      :: OUTGRID(:,:) => NULL()
      REAL*8, POINTER      :: INGRID (:,:) => NULL()
     
      !==================================================================
      ! READ_ERUP_VOLC begins here
      !==================================================================

      ! If the current year falls within the range of available data,
      ! get eruptive volcanic emissions (jaf, 10/15/09)
      IF ( ( INYEAR .GE. MNYEAR_VOLC )  .AND. 
     &     ( INYEAR .LE. MXYEAR_VOLC ) ) THEN

         ! Set date
         YYYYMMDD = 10000 * INYEAR + 100 * INMONTH + INDAY

         ! File name
         FILENAME = TRIM( DATA_DIR_1x1 )  //
     &              'volcano_SO2_201010/' //
     &              'YYYY/SO2_volc.erup.YYYYMM.generic.1x1'

         ! Replace YYYY/MM in the file name
         CALL EXPAND_DATE( FILENAME, YYYYMMDD, 000000 )

         ! TAU value (use this year)
         XTAU = GET_TAU0( INMONTH, INDAY, INYEAR )

         ! Echo output
         WRITE( 6, 100 ) TRIM( FILENAME )
 100     FORMAT( '     - READ_ERUP_VOLC: Reading ', a )

         ! Read eruptive emissions (kg/box/day)
         CALL READ_BPCH2( FILENAME, 'SO2-EV-$',        26,
     &                    XTAU,      I1x1,             J1x1-1,
     &                    LVOLC,     ARRAY_GEN(:,:,:), QUIET=.TRUE. )

         ! File with lat/lon edges for regridding
         LLFILENAME = TRIM( DATA_DIR_1x1 )           //
     &                'MAP_A2A_Regrid_201203/'       // 
     &                'MAP_A2A_latlon_generic1x1.nc'

         ! Array regridding has to be done level by level since the array
         ! has a different vertical dimension than expected
         DO L = 1, LVOLC

            ! Cast from REAL*4 to REAL*8
            ARRAY_GEN8 = ARRAY_GEN(:,:,L)

            ! Point to input & output arrays
            INGRID  => ARRAY_GEN8                      !Units kg 
            OUTGRID => ARRAY(:,:,L)
            
            ! Regrid
            CALL DO_REGRID_A2A( LLFILENAME, I1x1,    J1x1-1,
     &                          INGRID,     OUTGRID, IS_MASS=1,
     &                          netCDF=.TRUE.                   )

            ! Free pointer
            NULLIFY( INGRID, OUTGRID )

         ENDDO

         ! Unit conversion: [kg SO2/box/day] -> [kg SO2/box/s]
         EEV = ARRAY / ( 24.d0 * 3600.d0 )

      ! If no data are available for current year, set eruptive
      ! emissions to zero and print a warning message.
      ! Current data range is MNYEAR_VOLC - MXYEAR_VOLC (jaf, 10/15/09)
      ! Year parameters are defined at the beginning of the module.
      ! (ccc, 9/30/10)
      ELSE

         WRITE(6, 110) MNYEAR_VOLC, MXYEAR_VOLC
         WRITE(6, *) 'Eruptive SO2 emissions are being set to zero!'

 110     FORMAT( 'WARNING: Eruptive SO2 emissions only available for ', 
     &        i4,'-',i4,'. You are outside the window.')


         EEV = 0d0

      ENDIF

      END SUBROUTINE READ_ERUP_VOLC
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_anthro_sox
!
! !DESCRIPTION: Suborutine READ\_ANTHRO\_SOx reads the anthropogenic SOx from
!  disk,  and partitions it into anthropogenic SO2 and SO4. (rjp, bdf, bmy,
!  9/20/02, 10/31/08)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE READ_ANTHRO_SOx( Input_Opt, THISMONTH, NSEASON )
!
! !USES:
!
      USE BPCH2_MOD,            ONLY : GET_NAME_EXT_2D, GET_RES_EXT
      USE BPCH2_MOD,            ONLY : GET_TAU0,        READ_BPCH2
      USE DIRECTORY_MOD,        ONLY : DATA_DIR
      USE EDGAR_MOD,            ONLY : GET_EDGAR_ANTH_SO2
      USE EMEP_MOD,             ONLY : GET_EMEP_ANTHRO
      USE EMEP_MOD,             ONLY : GET_EUROPE_MASK
      USE FUTURE_EMISSIONS_MOD, ONLY : GET_FUTURE_SCALE_SO2ff
      USE GIGC_Input_Opt_Mod,   ONLY : OptInput
      USE GRID_MOD,             ONLY : GET_XMID,        GET_YMID
      USE GRID_MOD,             ONLY : GET_AREA_CM2
      USE RCP_MOD,              ONLY : GET_RCP_EMISSION
      USE STREETS_ANTHRO_MOD,   ONLY : GET_SE_ASIA_MASK
      USE STREETS_ANTHRO_MOD,   ONLY : GET_STREETS_ANTHRO
      USE TIME_MOD,             ONLY : GET_YEAR
      USE TRACER_MOD,           ONLY : XNUMOL
      USE TRACERID_MOD,         ONLY : IDTSO2, IDTSO4
      USE TRANSFER_MOD,         ONLY : TRANSFER_2D
      USE SCALE_ANTHRO_MOD,     ONLY : GET_ANNUAL_SCALAR
![eml
      USE LOGICAL_MOD,          ONLY : LHIST
!eml]
#if   defined( TOMAS )
      USE TRACERID_MOD,         ONLY : IDTNK1    ! For TOMAS microphysics
#endif

      USE CMN_SIZE_MOD               ! Size parameters
      USE CMN_O3_MOD                 ! FSCALYR
!
! !INPUT PARAMETERS: 
!
      TYPE(OptInput), INTENT(IN)    :: Input_Opt ! Input Options object
      INTEGER,        INTENT(IN)    :: THISMONTH ! Current month number (1-12)
      INTEGER,        INTENT(IN)    :: NSEASON   ! Season #: 1=DJF; 2=MAM; 
                                                 !           3=JJA; 4=SON

!  
! !REVISION HISTORY:
!  (1 ) Now use functions GET_XMID and GET_YMID to compute lon and lat
!        centers of grid box (I,J).  Now replace DXYP(JREF)*1d4 with routine
!        GET_AREA_CM2 of "grid_mod.f".  Now use functions GET_MONTH and
!        GET_YEAR of time_mod.f".  Now call READ_BPCH2 with QUIET=.TRUE. 
!        (bmy, 3/27/03)
!  (2 ) Now references DATA_DIR from "directory_mod.f".  Also removed
!        reference to CMN, it's not needed. (bmy, 7/20/04)
!  (3 ) Now read files from "sulfate_sim_200508/".  Now read data for both
!        GCAP and GEOS grids (bmy, 8/16/05)
!  (4 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (5 ) Now references XNUMOL from "tracer_mod.f" (bmy, 10/25/05)
!  (6 ) Now computes future SOx emissions (swu, bmy, 5/30/06)
!  (7 ) Now can read either EDGAR or GEIA emissions (avd, bmy, 7/14/06)
!  (8 ) Now overwrite David Streets' SO2, if necessary (yxw, bmy, 8/14/06)
!  (9 ) Now accounts for FSCLYR (phs, 3/17/08)
!  (9 ) Bug fix: Using tracer #30 in the call to GET_STREETS_ANTHRO can cause
!        problems when adding or removing species.  Replace w/ IDTNH3.
!        (dkh, 10/31/08)
!  (10) Account for multilevels emissions (amv, 10/07/2009)
!  (11 ) Use 1% SO4 fraction (out of total SOx) everywhere when running
!         with TOMAS aerosols. (win, 1/25/10)
!  22 Dec 2011 - M. Payer    - Added ProTeX headers
!  01 Mar 2012 - R. Yantosca - Now use GET_AREA_CM2(I,J,L) from grid_mod.F90
!  01 Mar 2012 - R. Yantosca - Now use GET_XMID(I,J,L) from grid_mod.F90
!  01 Mar 2012 - R. Yantosca - Now use GET_YMID(I,J,L) from grid_mod.F90
!  30 May 2013 - S. Farina   - Merge in TOMAS updates
!  22 Jul 2013 - M. Sulprizio- Now copy LRCP from Input_Opt
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER                       :: I, J, L, LL, IX, JX, IOS
      INTEGER, SAVE                 :: LASTYEAR = -99
      INTEGER                       :: SCALEYEAR
      REAL*4                        :: E_SOx(IIPAR,JJPAR,2)
      REAL*4                        :: ARRAY(IIPAR,JJPAR,1)
      REAL*8                        :: XTAU, Fe, X, Y, AREA_CM2
      REAL*8                        :: EDG_SO2
      CHARACTER(LEN=4)              :: SYEAR
      CHARACTER(LEN=255)            :: FILENAME

      ! For fields from Input_Opt
      LOGICAL                       :: LFUTURE
      LOGICAL                       :: LEDGARSOx
      LOGICAL                       :: LSTREETS
      LOGICAL                       :: LEMEP
      LOGICAL                       :: LRCP

      !=================================================================
      ! READ_ANTHRO_SOx begins here!
      !=================================================================

      ! Copy values from Input_Opt
      LFUTURE        = Input_Opt%LFUTURE
      LEDGARSOx      = Input_Opt%LEDGARSOx
      LSTREETS       = Input_Opt%LSTREETS
      LEMEP          = Input_Opt%LEMEP
      LRCP           = Input_Opt%LRCP
         
      IF ( FSCALYR < 0 ) THEN
         SCALEYEAR = GET_YEAR()
      ELSE
         SCALEYEAR = FSCALYR
      ENDIF


      IF ( LEDGARSOx .OR. LRCP ) THEN

         !==============================================================
         ! Use EDGAR SOx emissions

         ! NOTE: When TOMAS is being used, then do this: 
         !   Accoroding to a recent finding from GISS GCM-II', emitting 
         !   3% of SOx as primary SO4 results in overprediction of 
         !   aerosol number even in sulfate-seasalt simulation, which 
         !   should be rather lower estimate in absence of carbonaceous
         !   aerosols.  So try this 1% for now. (win, 1/25/10)
         !
         ! Otherwise, if TOMAS is not used, then partition SOx into 
         ! SO2 and SO4, according to the following fractions
         ! (cf Chin et al, 2000):
         ! 
         ! Europe     [ 36N-78N,  12.5W-60.0E ]:  5.0% of SOx is SO4
         !                                       95.0% of SOx is SO2   
         !
         ! N. America [ 26N-74N, 167.5W-52.5W ]:  1.4% of SOx is SO4
         !                                       98.6% of SOx is SO2
         !
         ! Everywhere else:                       3.0% of SOx is SO4
         !                                       97.0% of SOx is SO2
         !==============================================================

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, X, Y, AREA_CM2, EDG_SO2, Fe, LL ) 

         ! Loop over latitudes
         DO J = 1, JJPAR

            ! Loop over longitudes
            DO I = 1, IIPAR

               ! Longitude and latitude [degrees]
               X       = GET_XMID( I, J, 1 )
               Y       = GET_YMID( I, J, 1 )

               ! Get EDGAR SO2 emissions [kg/s]
               ! NOTE: Future emissions are already applied!
               EDG_SO2 = GET_EDGAR_ANTH_SO2( I, J, KG_S=.TRUE. ) 

               ! Use RCP SO2 emissions
               IF ( LRCP ) THEN

                  ! Get RCP SO2 emissions, molec/cm2/s
                  ! NOTE: Future emissions assumed to already be included
                  EDG_SO2 = GET_RCP_EMISSION( I, J, IDTSO2, 
     &                 LAND=.TRUE., SHIP=.FALSE. )

                  ! Convert molec/cm2/s -> kg/s
                  EDG_SO2 = EDG_SO2 * AREA_CM2 / XNUMOL(IDTSO2) 

                  ! RCP includes biofuels
                  ! However, if regional inventories without biofuels 
                  ! are used, then there will be no biofuel emissions 
                  ! from that region
                  ESO2_bf(I,J) = 0d0

               ENDIF

               ! If we are using David Streets' emissions ...
![eml Not if historical
               IF ( LSTREETS .and. .not. LHIST ) THEN
!eml]
                  ! If we are over the SE Asia region ...
                  IF ( GET_SE_ASIA_MASK( I, J ) > 0d0 ) THEN

                     ! Overwrite EDGAR SO2 w/ David Streets' [kg SO2/s]
                     EDG_SO2 = GET_STREETS_ANTHRO( I,      J, 
     &                                             IDTSO2, KG_S=.TRUE. )
                     
                     ! Streets 2006 includes biofuels.
                     IF ( SCALEYEAR >= 2006 ) ESO2_bf(I,J) = 0d0

                  ENDIF
               ENDIF

               ! If we are using EMEP over Europe...
![eml
               IF ( LEMEP .and. .not. LHIST ) THEN
!eml]

                  IF (GET_EUROPE_MASK(I,J) > 0d0) THEN

                     EDG_SO2 = GET_EMEP_ANTHRO( I,      J,
     $                                          IDTSO2, KG_S=.TRUE. )

                  ENDIF

               ENDIF

           
#if   defined( TOMAS )
               !--------------------------------------------------------
               ! For TOMAS microphysics: use 1% SO4 fraction
               !--------------------------------------------------------
               Fe = 0.01d0 
#else
               !--------------------------------------------------------
               ! For std GEOS-Chem: preserve existing code
               !--------------------------------------------------------

               ! Compute SO4/SOx fraction for EUROPE
               IF ( ( X >= -12.5 .and. X <= 60.0 )  .and.
     &              ( Y >=  36.0 .and. Y <= 78.0 ) ) THEN
                  Fe = 0.05d0

               ! Compute SO4/SOx fraction for NORTH AMERICA
               ELSE IF ( ( X >= -167.5 .and. X <= -52.5 )  .and.
     &                   ( Y >=   26.0 .and. Y <=  74.0 ) ) THEN
                  Fe = 0.014d0

               ! Compute SO4/SOx fraction for EVERYWHERE ELSE
               ELSE
                  Fe = 0.03d0

               ENDIF
#endif

               ! Store SO2 emission [kg SO2/s]
               ESO2_an(I,J,1) = EDG_SO2
               DO LL = 2,NOXLEVELS
                  ESO2_an(I,J,LL) = 0d0
               ENDDO

               ! Compute SO4 from SO2 [kg SO4/s]
               ESO4_an(I,J,1) = EDG_SO2 * Fe / ( 1.d0 - Fe )
               DO LL = 2, NOXLEVELS
                  ESO4_an(I,J,LL) = 0d0
               ENDDO
            ENDDO
         ENDDO
!$OMP END PARALLEL DO

      ELSE

         !==============================================================
         ! Use GEIA SOx emissions
         !==============================================================

         ! Define filename
         FILENAME = TRIM( DATA_DIR )                          //
     &              'fossil_200104/merge_nobiofuels.'         //
     &              GET_NAME_EXT_2D() // '.' // GET_RES_EXT() 
     
         ! Echo output
         WRITE( 6, 100 ) TRIM( FILENAME )
 100     FORMAT( '     - READ_ANTHRO_SOx: Reading ', a )

         ! Pick the right TAU value for the given season
         ! Seasons: 1=DJF, 2=MAM, 3=JJA, 4=SON
         SELECT CASE ( NSEASON )
            CASE ( 1 )
               XTAU = -744d0
            CASE ( 2 )
               XTAU = 1416d0
            CASE ( 3 )
               XTAU = 3624d0
            CASE ( 4 )
               XTAU = 5832d0
         END SELECT

         ! Read anthropogenic SOx [molec/cm2/s] 
         CALL READ_BPCH2( FILENAME, 'ANTHSRCE', 27, 
     &                    XTAU,      IIPAR,     JJPAR,     
     &                    2,         E_SOx,     QUIET=.TRUE. )

         !=================================================================
         ! Read in yearly SO2 scale factors here
         ! (For now we only have 1998, deal w/ other years later)
         !=================================================================

         ! Get annual scalar factor (amv, 08/24/07)
         CALL GET_ANNUAL_SCALAR( 73, 1985, SCALEYEAR, SOx_SCALE )

         !==============================================================
         ! Partition SOx into SO2 and SO4, according to the following
         ! fractions (cf Chin et al, 2000):
         ! 
         ! Europe     [ 36N-78N,  12.5W-60.0E ]:  5.0% of SOx is SO4
         !                                       95.0% of SOx is SO2   
         !
         ! N. America [ 26N-74N, 167.5W-52.5W ]:  1.4% of SOx is SO4
         !                                       98.6% of SOx is SO2
         !
         ! Everywhere else:                       3.0% of SOx is SO4
         !                                       97.0% of SOx is SO2
         !==============================================================
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, AREA_CM2, Y, X, Fe, LL )
         DO L = 1, 2
         DO J = 1, JJPAR
         DO I = 1, IIPAR

            ! Grid box surface area [cm2]
            AREA_CM2 = GET_AREA_CM2( I, J, L )
            
            ! Longitude and latitude [degrees]
            X = GET_XMID( I, J, L )
            Y = GET_YMID( I, J, L )

            ! First scale SOx to the given fossil fuel year
            E_SOx(I,J,L) = E_SOx(I,J,L) * SOx_SCALE(I,J)
            
            ! Compute future SOx emissions (if necessary)
            IF ( LFUTURE ) THEN
               E_SOx(I,J,L) = E_SOx(I,J,L)                  * 
     &                        GET_FUTURE_SCALE_SO2ff( I, J )
            ENDIF

            ! Compute SO4/SOx fraction for EUROPE
            IF ( ( X >= -12.5 .and. X <= 60.0 )  .and. 
     &           ( Y >=  36.0 .and. Y <= 78.0 ) ) THEN
               Fe = 0.05d0

            ! Compute SO4/SOx fraction for NORTH AMERICA
            ELSE IF ( ( X >= -167.5 .and. X <= -52.5 )  .and.   
     &                ( Y >=   26.0 .and. Y <=  74.0 ) ) THEN
               Fe = 0.014d0
 
            ! Compute SO4/SOx fraction for EVERYWHERE ELSE
            ELSE
               Fe = 0.03d0
             
            ENDIF
         
            ! Compute SO2 (tracer #2) from SOx
            ! Convert from [molec SOx/cm2/s] to [kg SO2/box/s]
            ESO2_an(I,J,L) = E_SOx(I,J,L) * ( 1.d0 - Fe ) * 
     &                       AREA_CM2 / XNUMOL(IDTSO2)            

            ! If we are using David Streets' emissions
            ! Remember: those include BF if Year is GE 2006
            IF ( LSTREETS ) THEN

               ! If we are over the SE Asia region
               IF ( GET_SE_ASIA_MASK( I, J ) > 0d0 ) THEN

                  ! Overwrite GEIA SO2 w/ David Streets' SO2 [kg SO2/s]
                  ESO2_an(I,J,1) = GET_STREETS_ANTHRO( I, J, IDTSO2, 
     &                                                 KG_S=.TRUE. )

                  ! Zero 2nd level of emissions
                  ESO2_an(I,J,2) = 0d0

               ENDIF
            ENDIF

            IF ( LEMEP ) THEN

               IF (GET_EUROPE_MASK(I,J) > 0d0 ) THEN

                  ESO2_an(I,J,1) = GET_EMEP_ANTHRO( I, J, IDTSO2,
     &                                              KG_S=.TRUE. ) 

                  ESO2_an(I,J,2) = 0d0
               ENDIF

            ENDIF

            ESO4_an(I,J,L) = ESO2_an(I,J,L) * Fe / (1.d0-Fe)

            IF ( NOXLEVELS > 2 ) THEN
               DO LL = 3,NOXLEVELS
                  ESO2_an(I,J,LL) = 0d0
                  ESO4_an(I,J,LL) = 0d0
               ENDDO
            ENDIF

         ENDDO
         ENDDO
         ENDDO
!$OMP END PARALLEL DO
      ENDIF

      END SUBROUTINE READ_ANTHRO_SOx
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_ocean_dms
!
! !DESCRIPTION: Subroutine READ\_OCEAN\_DMS reads seawater concentrations of
!  DMS (nmol/L). (rjp, bdf, bmy, 9/20/02, 10/3/05)
!\\
!\\
! !INTERFACE:
! 
      SUBROUTINE READ_OCEAN_DMS( THISMONTH )
!
! !USES:
!
      USE CMN_SIZE_MOD
      USE BPCH2_MOD,     ONLY : GET_NAME_EXT_2D, GET_RES_EXT
      USE BPCH2_MOD,     ONLY : GET_TAU0,        READ_BPCH2
      USE DIRECTORY_MOD, ONLY : DATA_DIR
      USE TRANSFER_MOD,  ONLY : TRANSFER_2D
!
! !INPUT PARAMETERS: 
!
      INTEGER, INTENT(IN) :: THISMONTH ! Current month number (1-12)
!
! !REVISION HISTORY:
!  (1 ) Extracted from old module routine SULFATE_READMON (bmy, 9/18/02)
!  (2 ) Now call READ_BPCH2 with QUIET=.TRUE. (bmy, 3/27/03)
!  (3 ) Now references DATA_DIR from "directory_mod.f" (bmy, 7/20/04)
!  (4 ) Now read files from "sulfate_sim_200508/".  Now read data for both
!        GCAP and GEOS grids (bmy, 8/16/05) 
!  (5 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  22 Dec 2011 - M. Payer    - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL*4              :: ARRAY(IIPAR,JJPAR,1)
      REAL*8              :: XTAU
      CHARACTER(LEN=255)  :: FILENAME

      !==================================================================
      ! READ_OCEAN_DMS begins here!
      !==================================================================

      ! File name
      FILENAME = TRIM( DATA_DIR )                         // 
     &           'sulfate_sim_200508/DMS_seawater.'       //
     &           GET_NAME_EXT_2D() // '.' // GET_RES_EXT()

      ! Echo output
      WRITE( 6, 100 ) TRIM( FILENAME )  
 100  FORMAT( '     - READ_OCEAN_DMS: Reading ', a ) 

      ! TAU value (use generic year 1985)
      XTAU = GET_TAU0( THISMONTH, 1, 1985 )

      ! Read ocean DMS [nmol/L]
      CALL READ_BPCH2( FILENAME, 'BIOGSRCE',    25, 
     &                 XTAU,      IIPAR,        JJPAR,      
     &                 1,         ARRAY(:,:,1), QUIET=.TRUE. ) 

      ! Cast from REAL*4 to REAL*8 
      CALL TRANSFER_2D( ARRAY(:,:,1), DMSo )
      
      END SUBROUTINE READ_OCEAN_DMS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_sst
!
! !DESCRIPTION: Subroutine READ\_SST reads monthly mean sea surface
!  temperatures. (rjp, bdf, bmy, 9/18/02, 7/13/09)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE READ_SST( THISMONTH, THISYEAR, State_Met )
!
! !USES:
!
      USE BPCH2_MOD,          ONLY : GET_NAME_EXT_2D, GET_RES_EXT
      USE BPCH2_MOD,          ONLY : GET_TAU0,        READ_BPCH2
      USE CMN_SIZE_MOD
      USE DIRECTORY_MOD,      ONLY : DATA_DIR,        DATA_DIR_1x1
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE REGRID_A2A_MOD,     ONLY : DO_REGRID_A2A
      USE TRANSFER_MOD,       ONLY : TRANSFER_2D
!
! !INPUT PARAMETERS: 
!
      INTEGER,        INTENT(IN)    :: THISMONTH ! Current month number (1-12)
      INTEGER,        INTENT(IN)    :: THISYEAR  ! Current 4-digit year
!
! !INPUT/OUTPUT PARAMETERS:
!
      TYPE(MetState), INTENT(INOUT) :: State_Met ! Meteorology State object
!
! !REVISION HISTORY:
!  (1 ) Extracted from old module routine SULFATE_READMON (bmy, 9/18/02)
!  (2 ) Now call READ_BPCH2 with QUIET=.TRUE. (bmy, 3/27/03)
!  (3 ) Now references DATA_DIR from "directory_mod.f" (bmy, 7/20/04)
!  (4 ) Now use interannual SST data from NOAA if present; otherwise use
!        climatological SST data.  Now read data for both GCAP and GEOS 
!        grids (bmy, 8/16/05) 
!  (5 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (6 ) Now read int'annual SST data on the GEOS 1x1 grid (bmy, 11/17/05)
!  (7 ) Last year of data is now 2008 (bmy, 7/13/09)
!  22 Dec 2011 - M. Payer    - Added ProTeX headers
!  13 Mar 2012 - M. Cooper   - Changed regrid algorithm to map_a2a
!  24 May 2012 - R. Yantosca - Fixed minor bugs in map_a2a implementation
!  24 Aug 2012 - R. Yantosca - DO_REGRID_A2A now reads netCDF input file
!  03 Jan 2013 - M. Payer    - Renamed PERAREA to IS_MASS in DO_REGRID_A2A 
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL*4                  :: ARRAY(IIPAR,JJPAR,1)
      REAL*4                  :: ARRAY1(I1x1,J1x1,1)
      REAL*8                  :: XTAU
      CHARACTER(LEN=4)        :: SYEAR
      CHARACTER(LEN=255)      :: FILENAME
      CHARACTER(LEN=255)      :: LLFILENAME
      REAL*8                  :: INGRID(I1x1,J1x1)

      !==================================================================
      ! READ_SST begins here!
      !==================================================================

      IF ( THISYEAR >= 1985 .and. THISYEAR <= 2008 ) THEN 
         
         !------------------------------------
         ! Use interannual SST data from NOAA
         ! Data exists for 1985 - 2004,
         ! Add other years as necessary
         !------------------------------------

         ! Make a string for THISYEAR
         WRITE( SYEAR, '(i4)' ) THISYEAR

         ! File name
         FILENAME = TRIM( DATA_DIR_1x1 )       // 
     &              'SST_200508/SST.geos.1x1.' // SYEAR

         ! Echo output
         WRITE( 6, 100 ) TRIM( FILENAME )  
 100     FORMAT( '     - READ_SST: Reading ', a ) 

         ! TAU value (use this year)
         XTAU = GET_TAU0( THISMONTH, 1, THISYEAR )

         ! Read sea surface temperature [K]
         CALL READ_BPCH2( FILENAME, 'GMAO-2D',      69, 
     &                    XTAU,      I1x1,          J1x1,     
     &                    1,         ARRAY1(:,:,1), QUIET=.TRUE. ) 

         ! File with lat/lon edges for regridding
         LLFILENAME = TRIM( DATA_DIR_1x1) //
     &                'MAP_A2A_Regrid_201203/MAP_A2A_latlon_geos1x1.nc'

         ! Regrid from 1x1 and cast to REAL*8
         INGRID = ARRAY1(:,:,1)
         CALL DO_REGRID_A2A( LLFILENAME, I1x1,    J1x1,
     &                       INGRID,     SSTEMP,  IS_MASS=0,
     &                       netCDF=.TRUE.                   )

      ELSE

         !-------------------------------
         ! Use climatological SST data 
         !-------------------------------

         ! File name
         FILENAME = TRIM( DATA_DIR )          // 
     &              'sulfate_sim_200508/SST.' // GET_NAME_EXT_2D() //
     &              '.'                       // GET_RES_EXT()

         ! Echo output
         WRITE( 6, 100 ) TRIM( FILENAME )  

         ! TAU value (use generic year 1985)
         XTAU = GET_TAU0( THISMONTH, 1, 1985 )

         ! Read sea surface temperature [K]
         CALL READ_BPCH2( FILENAME, 'DAO-FLDS',    5, 
     &                    XTAU,      IIPAR,        JJPAR,     
     &                    1,         ARRAY(:,:,1), QUIET=.TRUE. ) 

         ! Cast from REAL*4 to REAL*8 
         CALL TRANSFER_2D( ARRAY(:,:,1), SSTEMP )

      ENDIF

      ! Copy met field to State_Met
      State_Met%SST = SSTEMP

      END SUBROUTINE READ_SST
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_biofuel_so2
!
! !DESCRIPTION: Subroutine READ\_BIOFUEL\_SO2 reads monthly mean biomass
!  burning  emissions for SO2.  SOx is read in, and converted to SO2. 
!  (rjp, bdf, bmy, phs, 1/16/03, 12/23/08)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE READ_BIOFUEL_SO2( THISMONTH )
!
! !USES:
!
      USE BIOMASS_MOD,          ONLY : BIOMASS
      USE BPCH2_MOD,            ONLY : GET_NAME_EXT_2D, GET_RES_EXT
      USE BPCH2_MOD,            ONLY : GET_TAU0,        READ_BPCH2
      USE DIRECTORY_MOD,        ONLY : DATA_DIR,        DATA_DIR_1x1
      USE FUTURE_EMISSIONS_MOD, ONLY : GET_FUTURE_SCALE_SO2bf
      USE LOGICAL_MOD,          ONLY : LBIOMASS,        LFUTURE
      USE TIME_MOD,             ONLY : ITS_A_LEAPYEAR
      USE TRACER_MOD,           ONLY : XNUMOL
      USE TRACERID_MOD,         ONLY : IDTSO2
      USE TRANSFER_MOD,         ONLY : TRANSFER_2D
      USE REGRID_A2A_MOD,       ONLY : DO_REGRID_A2A
![eml
      USE LOGICAL_MOD,          ONLY : LHIST
      USE TIME_MOD,             ONLY : GET_HISTYR
!eml]
      
      USE CMN_SIZE_MOD               ! Size parameters
!
! !INPUT PARAMETERS: 
!
      INTEGER, INTENT(IN)           :: THISMONTH ! Current month number (1-12)

!
! !REVISION HISTORY:
!  (1 ) Extracted from old module routine SULFATE_READMON (bmy, 9/18/02)
!  (2 ) Modified molar ratio of biomass burning SO2 per CO.  Added SO2
!        emission from biofuel burning. (rjp, bmy, 1/16/03)
!  (3 ) Now replace DXYP(J+J0)*1d4 with routine GET_AREA_CM2 of "grid_mod.f"
!        Now replace MONTH with the argument THISMONTH.  Now call READ_BPCH2
!        with QUIET=.TRUE. (bmy, 3/27/03)
!  (4 ) Now references DATA_DIR from "directory_mod.f".  Also removed
!        references to CMN and CMN_SETUP. (bmy, 7/20/04)
!  (5 ) Now can read either seasonal or interannual biomass burning emissions.
!        Now references routines from both "logical_mod.f" and "time_mod.f".
!        Now reads SO2 biomass emissions directly rather than computing
!        it by mole fraction from CO. (rjp, bmy, 1/11/05)
!  (6 ) Now read data for both GCAP and GEOS grids (bmy, 8/16/05) 
!  (7 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (8 ) Now computes future biomass emissions, if necessary (swu, bmy, 5/30/06)
!  (9 ) Now only read the biofuel, we have moved the biomass-reading code
!        to "gc_biomass_mod.f" for compatibility with GFED2 biomass emissions
!        (bmy, 9/27/06)
!  (10) Now prevent seg fault if BIOMASS emissions are turned off.
!        (bmy, 10/3/06)
!  (11) Renamed READ_BIOFUEL_SO2, and move all biomass code to GET_BIOMASS_SO2
!        to account for several GFED2 products (yc, phs, 12/23/08)
!  (12) IDBSO2 is not used anymore (ccc, 01/29/10)
!  22 Dec 2011 - M. Payer    - Added ProTeX headers
!  01 Mar 2012 - R. Yantosca - Removed reference to GET_AREA_CM2
!  08 Mar 2012 - M. Payer    - Added modifications for historical emissions of
!                              SO2 (E. Leibensperger)
!  06 Apr 2012 - M. Payer    - Changed regrid algorithm to map_a2a (M. Cooper)
!  24 May 2012 - R. Yantosca - Fixed minor bugs in map_a2a implementation
!  24 Aug 2012 - R. Yantosca - DO_REGRID_A2A now reads netCDF input file
!  03 Jan 2013 - M. Payer    - Bug fix for regridding of historical emissions.
!                              Changed to PERAREA=1 since units are kg/yr.
!  03 Jan 2013 - M. Payer    - Renamed PERAREA to IS_MASS in DO_REGRID_A2A
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER                       :: I, J, THISYEAR
      REAL*4                        :: ARRAY(IIPAR,JJPAR,1)
      REAL*8                        :: BIOCO(IIPAR,JJPAR)
      REAL*8                        :: XTAU
      CHARACTER(LEN=4  )            :: CYEAR
      CHARACTER(LEN=255)            :: FILENAME
      CHARACTER(LEN=255)            :: LLFILENAME
      REAL*8, POINTER               :: OUTGRID(:,:) => NULL()
      REAL*8, POINTER               :: INGRID (:,:) => NULL()
![eml
      CHARACTER(LEN=4  )            :: CHISTYR
      REAL*4                        :: ESO2_GEN(I1x1,J1x1-1)
      REAL*8, TARGET                :: ESO2_GEN8(I1x1,J1x1-1)
      REAL*8, TARGET                :: TEMP_ED(IIPAR,JJPAR,1)
!eml]
      ! Days per month              
      REAL*8 :: NMDAY(12) = (/ 31d0, 28d0, 31d0, 30d0, 31d0, 30d0, 
     &                         31d0, 31d0, 30d0, 31d0, 30d0, 31d0/)

      !=================================================================
      ! READ_BIOFUEL_SO2 begins here!
      !=================================================================

      !=================================================================
      ! Compute biofuel SO2 from biofuel CO.  Use a molar 
      ! ratio of 0.0015 moles SO2/mole CO from biofuel burning. 
      ! (Table 2, [Andreae and Merlet, 2001])
      !=================================================================
      
      ! File name for biofuel burning file
![eml
      IF(.not. LHIST) THEN
!eml]
         FILENAME = TRIM( DATA_DIR )          // 
     &              'biofuel_200202/biofuel.' // GET_NAME_EXT_2D() //
     &              '.'                       // GET_RES_EXT()

         ! Echo info
         WRITE( 6, 100 ) TRIM( FILENAME )
 100     FORMAT( '     - READ_BIOFUEL_SO2: Reading ', a )

         ! Get TAU0 value (use generic year 1985)
         XTAU = GET_TAU0( 1, 1, 1985 )

         ! Read Biofuel burning of CO [kg/yr]
         CALL READ_BPCH2( FILENAME, 'BIOFSRCE',    4, 
     &                    XTAU,      IIPAR,        JJPAR,     
     &                    1,         ARRAY(:,:,1), QUIET=.TRUE.  ) 

         ! Cast from REAL*4 to REAL*8
         CALL TRANSFER_2D( ARRAY(:,:,1), BIOCO )

         !=================================================================
         ! Unit conversion to [kg SO2/s]
         !=================================================================
  
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J )
         ! Loop over longitudes
         DO J = 1, JJPAR

            ! Loop over latitudes
            DO I = 1, IIPAR

               ! Convert biofuel SO2 from [kg CO/yr] to [kg SO2/s]
               ESO2_bf(I,J) = ( BIOCO(I,J) * 64d-3 * 0.0015d0 /
     &                        ( 28d-3 * 86400.d0 * 365.25d0 ) )

               ! Apply future emissions to biofuel SO2, if necessary
               IF ( LFUTURE ) THEN
                  ESO2_bf(I,J) = ESO2_bf(I,J) * 
     &                           GET_FUTURE_SCALE_SO2bf( I, J )
               ENDIF
            ENDDO
         ENDDO
!$OMP END PARALLEL DO

![eml
      ELSE

         !-----------------------------------------------------------------
         ! Historical emissions of SO2
         !-----------------------------------------------------------------

         WRITE(CHISTYR,'(I4)') GET_HISTYR()

         ! Get TAU0 for this year
         XTAU = GET_TAU0( 1, 1, GET_HISTYR() )

         ! Read biofuel SO2
         FILENAME = TRIM( DATA_DIR_1x1 )                       // 
     &              'historical_emissions_201203/SO2/'         //
     &              'EDGAR.BCSO2.' // CHISTYR // 'generic.1x1'
        
         WRITE( 6, 101 ) TRIM( FILENAME )
 101     FORMAT( '     - READ_EDGAR_NOx: Reading ', a )
         CALL READ_BPCH2( TRIM(FILENAME), 'EDGAR-2D', 26,
     &                    XTAU,           I1x1,       J1x1-1,
     &                    1,              ESO2_GEN,   QUIET=.TRUE. )

         ! Cast from REAL*4 to REAL*8
         ESO2_GEN8 = ESO2_GEN

         ! File with lat/lon edges for regridding
         LLFILENAME = TRIM( DATA_DIR_1x1 )           //
     &                'MAP_A2A_Regrid_201203/'       //
     &                'MAP_A2A_latlon_generic1x1.nc'

         ! Point to input & output arrays
         INGRID  => ESO2_GEN8
         OUTGRID => TEMP_ED(:,:,1)

         ! Regrid to current model resolution [kg/yr]
         CALL DO_REGRID_A2A( LLFILENAME, I1x1,    J1x1-1, 
     &                       INGRID,     OUTGRID, IS_MASS=1,
     &                       netCDF=.TRUE.                   )

         ! Free pointers
         NULLIFY( INGRID, OUTGRID )

         ! Convert from kg SO2/yr to kg SO2/s
         ESO2_bf = TEMP_ED(:,:,1) / (365.25d0 * 86400.d0)

      ENDIF
!eml]

      END SUBROUTINE READ_BIOFUEL_SO2
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_biomass_so2
!
! !DESCRIPTION: Subroutine GET\_BIOMASS\_SO2 retrieve monthly/8-day/3hr biomass
!  burning emissions for SO2.  (yc, phs, 12/23/08)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE GET_BIOMASS_SO2
!
! !USES:
!
      USE BIOMASS_MOD,          ONLY : BIOMASS
      USE CMN_SIZE_MOD
      USE TRACERID_MOD,         ONLY : IDBSO2
      USE GRID_MOD,             ONLY : GET_AREA_CM2
      USE TRACER_MOD,           ONLY : XNUMOL
      USE TRACERID_MOD,         ONLY : IDTSO2
      USE TRANSFER_MOD,         ONLY : TRANSFER_2D
! 
! !REVISION HISTORY:
!  (1 ) Extracted from old module subroutine READ_BIOMASS_SO2
!        (yc, phs, 12/23/08)
!  (2 ) IDBSO2 is now in tracerid_mod.f (fp, 6/2009)
!  22 Dec 2011 - M. Payer    - Added ProTeX headers
!  01 Mar 2012 - R. Yantosca - Now use GET_AREA_CM2(I,J,L) from grid_mod.F90
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER                       :: I, J
      REAL*8                        :: CONV
                                    
      !=================================================================
      ! GET_BIOMASS_SO2 begins here!
      !=================================================================

      ! Unit conversion to [kg SO2/s]
  
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, CONV )

      ! Loop over surface boxes
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Conversion factor for [cm2 * kg/molec]
         CONV = GET_AREA_CM2( I, J, 1 ) / XNUMOL(IDTSO2)

         ! Convert biomass SO2 from [molec SO2/cm2/s] -> [kg SO2/s] 
         ! NOTE: Future scale has already been applied by this point
         ESO2_bb(I,J) = BIOMASS(I,J,IDBSO2) * CONV
            
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      END SUBROUTINE GET_BIOMASS_SO2
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_aircraft_so2
!
! !DESCRIPTION: Subroutine READ\_AIRCRAFT\_SO2 reads monthly mean aircraft fuel
!  emissions and converts them to SO2 emissions. (rjp, bdf, bmy, 9/18/02,
!  10/3/05)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE READ_AIRCRAFT_SO2( THISMONTH, State_Met )
!
! !USES:
!
      USE BPCH2_MOD,          ONLY : GET_RES_EXT, GET_TAU0, READ_BPCH2
      USE CMN_SIZE_MOD
      USE DIRECTORY_MOD,      ONLY : DATA_DIR 
      USE FILE_MOD,           ONLY : IOERROR
      USE GIGC_State_Met_Mod, ONLY : MetState
!
! !INPUT PARAMETERS: 
!
      INTEGER,        INTENT(IN) :: THISMONTH   ! Current month number (1-12)
      TYPE(MetState), INTENT(IN) :: State_Met   ! Meteorology State object
!
! !REVISION HISTORY:
!  (1 ) Extracted from old module routine SULFATE_READMON (bmy, 9/18/02)
!  (2 ) Now references DATA_DIR from "directory_mod.f" (bmy, 7/20/04)
!  (3 ) Now read files from "sulfate_sim_200508/".  Now read data for both 
!        GCAP and GEOS grids (bmy, 8/16/05)
!  (4 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  22 Dec 2011 - M. Payer    - Added ProTeX headers
!  02 Jul 2012 - R. Yantosca - Add simple kludge that prevents OOB error when
!                              reading SO2 aircraft data for NA nested grid
!  01 Aug 2012 - R. Yantosca - Add reference to findFreeLUN from inqure_mod.F90
!  03 Aug 2012 - R. Yantosca - Move calls to findFreeLUN out of DEVEL block
!  09 Nov 2012 - M. Payer    - Replaced all met field arrays with State_Met
!                              derived type object
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER            :: I, IOS, J, K, L
      INTEGER            :: IU_FILE
      REAL*8             :: ACSO2(IIPAR,JJPAR,20)
      REAL*8             :: FAC, FUEL, DZ(LLPAR), ZH(0:LLPAR)
      CHARACTER(LEN=255) :: FILENAME

      ! Month names
      CHARACTER(LEN=3)   :: CMONTH(12) = (/'jan', 'feb', 'mar', 'apr', 
     &                                     'may', 'jun', 'jul', 'aug',
     &                                     'sep', 'oct', 'nov', 'dec'/)

      !=================================================================
      ! READ_AIRCRAFT_SO2 begins here!
      !=================================================================
      
      ! Zero arrays
      ESO2_ac = 0d0
      ACSO2   = 0d0
      
      ! File name
      FILENAME = TRIM( DATA_DIR )               // 
     &           'sulfate_sim_200508/aircraft.' // GET_RES_EXT() //
     &           '.1992.'                       // CMONTH(THISMONTH)

      ! Echo output
      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( '     - READ_AIRCRAFT_SO2: Reading ', a )     

      !=================================================================
      ! Read aircraft emissions.  These are fuel burned in [kg/box/day],
      ! from AEAP for 1992.  SO2 emission is calculated by assuming    
      ! an emission index EI of 1.0, i.e., 1g of SO2 emitted per kg    
      ! of fuel burned.  It is also assumed that there is no diurnal   
      ! variation of emission rate. Convert to [kg SO2/box/s]. 
      !=================================================================

      ! Find a free file LUN
      IU_FILE = findFreeLUN()

      ! Open file 
      OPEN( IU_FILE, FILE=FILENAME, STATUS='OLD', IOSTAT=IOS )
      IF ( IOS > 0 ) CALL IOERROR( IOS, IU_FILE, 'read_aircraft_so2:1' )

      ! Read header line
      READ( IU_FILE, '(/)', IOSTAT=IOS )
      IF ( IOS > 0 ) CALL IOERROR( IOS, IU_FILE, 'read_aircraft_so2:2' )
      
      ! Read data values until an EOF is found
      DO 
         READ( IU_FILE, '(3i4,e11.3)', IOSTAT=IOS ) I, J, L, FUEL

         ! EOF encountered
         IF ( IOS < 0 ) EXIT

         ! I/O error condition
         IF ( IOS > 0 ) THEN
            CALL IOERROR( IOS, IU_FILE, 'read_aircraft_so2:3' )
         ENDIF

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%% KLUDGE FOR NORTH AMERICAN NESTED GRID (bmy, 7/2/12)
!%%% For some reason there seems to be bad data within the N. American
!%%% nested grid data files.  Some of the longitude indices are out of bounds
!%%% for the grid.  Put in a simple kludge to just skip over these boxes.
!%%% We are going to be reinventing how emissions get done in GEOS-Chem,
!%%% so it's probably not worth recreating the the file at this time.
!%%%
#if   defined( NESTED_NA )        
         IF ( I > IIPAR ) CYCLE
         IF ( J > JJPAR ) CYCLE
#endif         
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

         ! Unit conversion: [kg Fuel/box/day] -> [kg SO2/box/s]
         ! Assuming an emission index of 1.0, 
         ! 1 g SO2 / kg fuel burned [Weisenstein et al., 1996]
         ACSO2(I,J,L+1) = 1.d-3 * FUEL / ( 24.d0 * 3600d0 )
      ENDDO

      ! Close file
      CLOSE( IU_FILE )

      !=================================================================
      ! Interpolate from the 1-km grid to the given GEOS-CHEM grid
      ! NOTE: we need to account for window grids (bmy, 9/20/02)
      !=================================================================
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! ACSO2 is the aircraft SO2 on the 1-km vertical grid
         FUEL = SUM( ACSO2(I,J,:) )
         IF ( FUEL < 1d-20 ) CYCLE

         ! There are 20 1-km levels
         DO K = 1, 20

            ! Initialize
            ZH(0) = 0.d0

            ! Loop over levels
            DO L = 1, LLPAR

               ! Altitude of top edge of level L, from ground [km]
               ZH(L) = ZH(L-1) + ( State_Met%BXHEIGHT(I,J,L) * 1d-3 )
               
               IF ( ZH(L-1) > DBLE(K)   ) EXIT
               IF ( ZH(L  ) < DBLE(K-1) ) CYCLE
               
               IF ( ZH(L) < DBLE(K) ) THEN
                  FAC            = ZH(L) - MAX( ZH(L-1), DBLE(K-1) )
                  ESO2_ac(I,J,L) = ESO2_ac(I,J,L) + ACSO2(I,J,K) * FAC
               ELSE
                  FAC            = DBLE(K) - MAX( ZH(L-1), DBLE(K-1) )
                  ESO2_ac(I,J,L) = ESO2_ac(I,J,L) + ACSO2(I,J,K) * FAC
                  EXIT
               ENDIF		     
            ENDDO
         ENDDO     
      ENDDO
      ENDDO

      END SUBROUTINE READ_AIRCRAFT_SO2
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_ship_so2
!
! !DESCRIPTION: Subroutine READ\_SHIP\_SO2 reads in ship SO2 emissions, from
!  either Corbett et al or EDGAR inventories. (bec, qli, 10/01/03, 7/14/06)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE READ_SHIP_SO2( Input_Opt, THISMONTH )
!
! !USES:
!
      USE ARCTAS_SHIP_EMISS_MOD,ONLY : GET_ARCTAS_SHIP
      USE ICOADS_SHIP_MOD,      ONLY : GET_ICOADS_SHIP   !(cklee, 7/09/09)
      USE BPCH2_MOD,            ONLY : GET_NAME_EXT_2D, GET_RES_EXT
      USE BPCH2_MOD,            ONLY : GET_TAU0,        READ_BPCH2
      USE CMN_SIZE_MOD
      USE EDGAR_MOD,            ONLY : GET_EDGAR_SHIP_SO2
      USE EMEP_MOD,             ONLY : GET_EMEP_ANTHRO, GET_EUROPE_MASK
      USE FUTURE_EMISSIONS_MOD, ONLY : GET_FUTURE_SCALE_SO2ff
      USE GIGC_Input_Opt_Mod,   ONLY : OptInput
      USE GRID_MOD,             ONLY : GET_AREA_CM2
      USE RCP_MOD,              ONLY : GET_RCP_EMISSION
      USE TRACER_MOD,           ONLY : XNUMOL
      USE TRACERID_MOD,         ONLY : IDTSO2
      USE TRANSFER_MOD,         ONLY : TRANSFER_2D
!
! !INPUT PARAMETERS: 
!
      TYPE(OptInput), INTENT(IN) :: Input_Opt
      INTEGER,        INTENT(IN) :: THISMONTH ! Current month (1-12)
!
! !REVISION HISTORY:
!  (1 ) Now references DATA_DIR from "directory_mod.f" (bmy, 7/20/04)
!  (2 ) Now read files from "sulfate_sim_200508/".  Now read data for both 
!       GCAP and GEOS grids. (bmy, 8/16/05)
!  (3 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (4 ) Now references XNUMOL from "tracer_mod.f" (bmy, 10/25/05)
!  (5 ) Now get EDGAR ship SO2 emissions if necessary.  Also apply future
!       emissions scale factors to the default Corbett et al ship emissions.
!       (avd, bmy, 7/14/06)
!  (6 ) Now references GET_ARCTAS_HIP from 'arctas_ship_emiss_mod.f" and
!       GET_EMEP_ANTHRO to get ARCTAS and EMEP SO2 ship emissions (phs,12/5/08)
!  (7 ) Now get ICOADS ship SO2 if necessary (phs, cklee, 6/30/09) 
!  22 Dec 2011 - M. Payer    - Added ProTeX headers
!  01 Mar 2012 - R. Yantosca - Now use GET_AREA_CM2(I,J,L) from grid_mod.F90
!  17 Jun 2013 - R. Yantosca - Now take values from Input_Opt
!  22 Jul 2013 - M. Sulprizio- Now copy LRCPSHIP from Input_Opt
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER                       :: I, J
      REAL*4                        :: ARRAY(IIPAR,JJPAR,1)
      REAL*4                        :: SHIPSO2(IIPAR,JJPAR)
      REAL*8                        :: XTAU, AREA_CM2
      CHARACTER (LEN=255)           :: FILENAME

      ! For fields from Input_Opt
      CHARACTER(LEN=255)            :: DATA_DIR
      LOGICAL                       :: LEDGARSHIP
      LOGICAL                       :: LEMEPSHIP
      LOGICAL                       :: LFUTURE
      LOGICAL                       :: LARCSHIP
      LOGICAL                       :: LSHIPSO2
      LOGICAL                       :: LICOADSSHIP
      LOGICAL                       :: LRCPSHIP

      !=================================================================
      ! READ_SHIP_SO2 begins here!
      !=================================================================

      ! Copy values from Input_Opt
      DATA_DIR    = Input_Opt%DATA_DIR
      LEDGARSHIP  = Input_Opt%LEDGARSHIP
      LEMEPSHIP   = Input_Opt%LEMEPSHIP
      LFUTURE     = Input_Opt%LFUTURE
      LARCSHIP    = Input_Opt%LARCSHIP
      LSHIPSO2    = Input_Opt%LSHIPSO2
      LICOADSSHIP = Input_Opt%LICOADSSHIP
      LRCPSHIP    = Input_Opt%LRCPSHIP

      ! Reset
      ESO2_sh = 0D0

      ! Test for EDAGR last, since this is default inventory by design.
      ! So we can still use EDGAR SHIP to get ship-NOX and CO, and
      ! overwrite ship-SO2 with ARCTAS or Colbert (phs, 12/5/08)
      ! Does the same with ICOADS (phs, 7/21/09)

      !-----------------------------------------------------------
      ! Use ARCTAS SHIP emissions (EDGAR 2006 update) 
      !-----------------------------------------------------------
      IF ( LARCSHIP ) THEN
            
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J )
         DO J = 1, JJPAR
         DO I = 1, IIPAR

         ! Read ARCTAS SO2 emissions in [kg SO2/BOX/s]
         ESO2_sh(I,J) = GET_ARCTAS_SHIP( I, J, IDTSO2,  KG_S=.TRUE. )
               
            IF ( LEMEPSHIP ) THEN
               IF ( GET_EUROPE_MASK(I,J) > 0d0 )
     $           ESO2_sh(I,J) = GET_EMEP_ANTHRO(I, J, IDTSO2, 
     &                                          KG_S=.TRUE.,
     $                                          SHIP=.TRUE.)
            ENDIF

         ENDDO 
         ENDDO
!$OMP END PARALLEL DO

      !-----------------------------------------------------------
      ! Or RCP SHIP emissions 
      ! (cdh, 10/14/11)
      !-----------------------------------------------------------
      ELSE IF ( LRCPSHIP ) THEN
            
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J )
         DO J = 1, JJPAR
         DO I = 1, IIPAR

               ! Grid box surface area [cm2]
               AREA_CM2 = GET_AREA_CM2( I, J, 1 )

               ! Read RCP SO2 emissions in [molec/cm2/s]
               ESO2_sh(I,J) = GET_RCP_EMISSION( I, J, IDTSO2,  
     &                                          SHIP=.TRUE. )
               
               ! Convert molec/cm2/s -> kg SO2/BOX/s
               ESO2_sh(I,J) = ESO2_sh(I,J) * AREA_CM2 / XNUMOL(IDTSO2)

               IF ( LEMEPSHIP ) THEN
                  IF ( GET_EUROPE_MASK(I,J) > 0d0 )
     $                 ESO2_sh(I,J) = GET_EMEP_ANTHRO(I, J, IDTSO2, 
     &                 KG_S=.TRUE.,
     $                 SHIP=.TRUE.)
               ENDIF

            ENDDO 
         ENDDO
!$OMP END PARALLEL DO

      !----------------------------------------
      ! Or Corbett et al ship SO2 emissions
      !----------------------------------------
      ELSE IF ( LSHIPSO2 ) THEN

         ! Filename
         FILENAME = TRIM( DATA_DIR )           // 
     &           'sulfate_sim_200508/shipSOx.' // GET_NAME_EXT_2D() //
     &           '.'                           // GET_RES_EXT()

         ! Echo some information to the standard output
         WRITE( 6, 110 ) TRIM( FILENAME )
 110     FORMAT( '     - READ_SHIP_SO2 ', a )
      
         ! TAU value at the beginning of this month
         XTAU = GET_TAU0( THISMONTH, 1, 1985 )
      
         ! Read in this month's ship SO2 emissions [molec SO2/cm2/s]
         CALL READ_BPCH2( FILENAME, 'SO2-SHIP',     26,  
     &                    XTAU,      IIPAR,         JJPAR, 
     &                    1,         ARRAY(:,:,1),  QUIET=.TRUE. )

         ! Cast from REAL*4 to REAL*8
         CALL TRANSFER_2D( ARRAY(:,:,1), SHIPSO2 )

         ! Loop over surface boxes
         DO J = 1, JJPAR
         DO I = 1, IIPAR

            ! Grid box surface area [cm2]
            AREA_CM2 = GET_AREA_CM2( I, J, 1 )

            ! Convert [molec SO2/cm2/s] to [kg SO2/box/s]
            ESO2_sh(I,J) = SHIPSO2(I,J) * AREA_CM2 / XNUMOL(IDTSO2)

            ! Apply future emissions (if necessary)
            IF ( LFUTURE ) THEN
               ESO2_sh(I,J) = ESO2_sh(I,J) *
     &                        GET_FUTURE_SCALE_SO2ff( I, J ) 
            ENDIF

            IF ( LEMEPSHIP ) THEN
               IF ( GET_EUROPE_MASK(I,J) > 0d0 )
     &              ESO2_sh(I,J) = GET_EMEP_ANTHRO(I, J, IDTSO2,
     &                                     KG_S=.TRUE., SHIP=.TRUE.)
            ENDIF               
             
         ENDDO 
         ENDDO

      !-----------------------------------------------------------
      ! Test for EDGAR ship emissions
      !-----------------------------------------------------------
      ELSE IF ( LEDGARSHIP ) THEN 

         !----------------------------------------
         ! Use EDGAR ship SO2 emissions
         !----------------------------------------

         ! Get EDGAR ship SO2 [kg SO2/box/s]
         ! NOTE: Future emissions have already been applied! 
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J )
         DO J = 1, JJPAR
         DO I = 1, IIPAR

            ESO2_sh(I,J) = GET_EDGAR_SHIP_SO2( I, J, KG_S=.TRUE. )

            IF ( LEMEPSHIP ) THEN
               IF ( GET_EUROPE_MASK(I,J) > 0d0 )
     $           ESO2_sh(I,J) = GET_EMEP_ANTHRO(I, J, IDTSO2, 
     &                                          KG_S=.TRUE.,
     $                                          SHIP=.TRUE.)
            ENDIF
         ENDDO 
         ENDDO
!$OMP END PARALLEL DO
         

      !----------------------------------------
      ! Use ICOADS ship SOx emissions   !(cklee, 6/24/09)
      ! Replace the above ship emissions 
      !----------------------------------------
      ELSE IF ( LICOADSSHIP ) THEN

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J )
         DO J = 1, JJPAR
         DO I = 1, IIPAR

         ! Read ICOADS emissions in [kg SO2/box/s]
         ESO2_sh(I,J) = GET_ICOADS_SHIP( I, J, IDTSO2,  KG_S=.TRUE. )
               
            IF ( LEMEPSHIP ) THEN
               IF ( GET_EUROPE_MASK(I,J) > 0d0 )
     $           ESO2_sh(I,J) = GET_EMEP_ANTHRO(I, J, 26, KG_S=.TRUE.,
     $                                          SHIP=.TRUE.)
            ENDIF

         ENDDO 
         ENDDO
!$OMP END PARALLEL DO
         
      ENDIF

      END SUBROUTINE READ_SHIP_SO2
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_anthro_nh3
!
! !DESCRIPTION: Subroutine READ\_ANTHRO\_NH3 reads the monthly mean 
!  anthropogenic NH3 emissions from disk and converts to [kg NH3/box/s]. 
!  (rjp, bdf, bmy, 9/20/02, 3/5/10)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE READ_ANTHRO_NH3( Input_Opt, THISMONTH )
!
! !USES:
!
      USE BPCH2_MOD,            ONLY : GET_NAME_EXT_2D, GET_RES_EXT
      USE BPCH2_MOD,            ONLY : GET_TAU0,        READ_BPCH2
      USE CMN_SIZE_MOD
      USE DIRECTORY_MOD,        ONLY : DATA_DIR
      USE EMEP_MOD,             ONLY : GET_EMEP_ANTHRO
      USE EMEP_MOD,             ONLY : GET_EUROPE_MASK
      USE FUTURE_EMISSIONS_MOD, ONLY : GET_FUTURE_SCALE_NH3an
      USE GIGC_Input_Opt_Mod,   ONLY : OptInput
      USE GRID_MOD,             ONLY : GET_AREA_CM2
      USE RCP_MOD,              ONLY : GET_RCP_EMISSION
      USE STREETS_ANTHRO_MOD,   ONLY : GET_SE_ASIA_MASK
      USE STREETS_ANTHRO_MOD,   ONLY : GET_STREETS_ANTHRO
      USE TRACERID_MOD,         ONLY : IDTNH3
      USE TRANSFER_MOD,         ONLY : TRANSFER_2D
!
! !INPUT PARAMETERS: 
!
      TYPE(OptInput), INTENT(IN)    :: Input_Opt ! Input Options object
      INTEGER, INTENT(IN)           :: THISMONTH ! Current month number (1-12)
!
! !REVISION HISTORY:
!  (1 ) Renamed from NH3_READ to READ_ANTHRO_NH3.  Also updated comments,
!        made cosmetic changes. (bmy, 9/20/02)
!  (2 ) Changed filename to NH3_anthsrce.geos.*.  Also now reads data under
!        category name "NH3-ANTH". (rjp, bmy, 3/23/03)
!  (3 ) Now reads from NH3emis.monthly.geos.* files.  Now call READ_BPCH2
!        with QUIET=.TRUE. (bmy, 3/27/03)
!  (4 ) Now references DATA_DIR from "directory_mod.f" (bmy, 7/20/04)
!  (5 ) Now read files from "sulfate_sim_200508/". Now read data for both 
!        GCAP and GEOS grids. (bmy, 8/16/05)
!  (5 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (6 ) Now compute future emissions, if necessary (swu, bmy, 5/30/06)
!  (7 ) Now overwrite w/ David Streets' NH3, if necessary (yxw, bmy, 8/17/06)
!  (8 ) Bug fix: Using tracer #30 in the call to GET_STREETS_ANTHRO can cause
!        problems when adding or removing species.  Replace w/ IDTNH3.
!        (dkh, 10/31/08)
!  (9 ) Now check if NH3 Streets is available (phs, 12/10/08)  
!  (10) Bug fix -- STREETS needs to be PRIVATE (dkh, bmy, 3/5/10)
!  22 Dec 2011 - M. Payer    - Added ProTeX headers
!  22 Jul 2013 - M. Sulprizio- Now copy LRCP from Input_Opt
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER                       :: I, J, DAY_NUM
      REAL*4                        :: ARRAY(IIPAR,JJPAR,1)
      REAL*8                        :: AREA_CM2, EPA_NEI, XTAU
      REAL*8                        :: XNUMOL(Input_Opt%N_TRACERS)

      ! Days per month
      REAL*8              :: NMDAY(12) = (/ 31d0, 28d0, 31d0, 30d0,
     &                                      31d0, 30d0, 31d0, 31d0, 
     &                                      30d0, 31d0, 30d0, 31d0 /)

      CHARACTER(LEN=255)            :: FILENAME
      REAL*8                        :: STREETS

      ! For fields from Input_Opt
      LOGICAL                :: LFUTURE
      LOGICAL                :: LSTREETS
      LOGICAL                :: LEMEP
      LOGICAL                :: LRCP

      !=================================================================
      ! READ_ANTHRO_NH3 begins here!
      !=================================================================

      ! Copy fields from INPUT_OPT to local variables for use below
      LFUTURE       = Input_Opt%LFUTURE
      LSTREETS      = Input_Opt%LSTREETS
      LEMEP         = Input_Opt%LEMEP
      LRCP          = Input_Opt%LRCP
      XNUMOL        = Input_Opt%XNUMOL

      !-----------------------------------------------------------
      ! Use RCP NH3 emissions 
      ! (mpayer, cdh, 6/28/12)
      !-----------------------------------------------------------   
      IF ( LRCP ) THEN

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, AREA_CM2 )
         DO J = 1, JJPAR
         DO I = 1, IIPAR

               ! Grid box surface area [cm2]
               AREA_CM2 = GET_AREA_CM2( I, J, 1 )

               ! Read RCP NH3 emissions in [molec/cm2/s]
               ENH3_an(I,J) = GET_RCP_EMISSION( I, J, IDTNH3,  
     &                                          LAND=.TRUE. )
               
               ! Convert molec/cm2/s -> kg NH3/BOX/s
               ENH3_an(I,J) = ENH3_an(I,J) * AREA_CM2 / XNUMOL(IDTNH3)
               
            ENDDO 
         ENDDO
!$OMP END PARALLEL DO

      ELSE

         ! File name
         FILENAME = TRIM( DATA_DIR )                         //
     &              'sulfate_sim_200508/NH3_anthsrce.'       //
     &              GET_NAME_EXT_2D() // '.' // GET_RES_EXT()

         ! Echo output
         WRITE( 6, 100 ) TRIM( FILENAME )
 100     FORMAT( '     - READ_ANTHRO_NH3: Reading ', a )
      
         ! Get TAU value (use year 1990, the year of the data!)
         XTAU = GET_TAU0( THISMONTH, 1, 1990 )
	
         ! Read 1990 NH3 emissions [kg N/box/mon]
         CALL READ_BPCH2( FILENAME, 'NH3-ANTH',    29,  
     &                    XTAU,      IIPAR,        JJPAR,       
     &                    1,         ARRAY(:,:,1), QUIET=.TRUE. )

         ! Cast from REAL*4 to REAL*8
         CALL TRANSFER_2D( ARRAY(:,:,1), ENH3_an )

         ! Convert from [kg N/box/mon] to [kg NH3/box/s]
         ENH3_an = ENH3_an * ( 17.d0 / 14.d0 ) 
     &             / ( NMDAY(THISMONTH) * 86400.d0 )


      ENDIF

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, STREETS )
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Compute future NH3an emissions, if necessary
         ! Moved here since Streets and EMEP should have already
         ! applied FUTURE scale factors if needed 
         IF ( LFUTURE ) THEN
            ENH3_an(I,J) = ENH3_an(I,J) *
     &                     GET_FUTURE_SCALE_NH3an( I, J )
         ENDIF

         ! If we are using David Streets' emissions ...
         IF ( LSTREETS ) THEN

            ! If we are over the SE Asia region ...
            IF ( GET_SE_ASIA_MASK( I, J ) > 0d0  ) THEN

               ! Overwrite with David Streets emissions if available
               ! [kg NH3/s]
               STREETS = GET_STREETS_ANTHRO( I,      J, 
     &                                       IDTNH3, KG_S=.TRUE.)

               IF ( .not. ( STREETS < 0d0 ) )
     $              ENH3_an(I,J) = STREETS
               
            ENDIF
         ENDIF

         IF ( LEMEP ) THEN
            IF ( GET_EUROPE_MASK(I,J) > 0d0) THEN
               ENH3_an(I,J) = GET_EMEP_ANTHRO(I,J,IDTNH3,KG_S=.TRUE.)
            ENDIF
         ENDIF

      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      END SUBROUTINE READ_ANTHRO_NH3
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_natural_nh3
!
! !DESCRIPTION: Subroutine READ\_NATURAL\_NH3 reads the monthly mean natural 
!  NH3 emissions from disk and converts to [kg NH3/box/s]. 
!  (rjp, bdf, bmy, 9/20/02, 10/3/05)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE READ_NATURAL_NH3( THISMONTH )
!
! !USES:
!
      USE BPCH2_MOD,     ONLY : GET_NAME_EXT_2D, GET_RES_EXT
      USE BPCH2_MOD,     ONLY : GET_TAU0,        READ_BPCH2
      USE CMN_SIZE_MOD
      USE DIRECTORY_MOD, ONLY : DATA_DIR
      USE TRANSFER_MOD,  ONLY : TRANSFER_2D
!
! !INPUT PARAMETERS: 
!
      INTEGER, INTENT(IN) :: THISMONTH ! Current month number (1-12)
!
! !REVISION HISTORY: 
!  (1 ) Updated FORMAT string.  Now also call READ_BPCH2 with QUIET=.TRUE.
!        (bmy, 4/8/03)
!  (2 ) Now references DATA_DIR from "directory_mod.f" (bmy, 7/20/04)
!  (3 ) Now read files from "sulfate_sim_200508/".  Now read data for both 
!        GCAP and GEOS grids. (bmy, 8/16/05)
!  (4 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  22 Dec 2011 - M. Payer    - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL*4              :: ARRAY(IIPAR,JJPAR,1)
      REAL*8              :: XTAU
      REAL*8              :: NMDAY(12) = (/ 31d0, 28d0, 31d0, 30d0,
     &                                      31d0, 30d0, 31d0, 31d0, 
     &                                      30d0, 31d0, 30d0, 31d0 /)
      CHARACTER(LEN=255)  :: FILENAME

      !=================================================================
      ! READ_NATURAL_NH3 begins here!
      !=================================================================

      ! File name
      FILENAME = TRIM( DATA_DIR )                         //
     &           'sulfate_sim_200508/NH3_natusrce.'       //
     &           GET_NAME_EXT_2D() // '.' // GET_RES_EXT()

      ! Echo output
      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( '     - READ_NATURAL_NH3: Reading ', a )
      
      ! Get TAU value (use year 1990, the year of the data!)
      XTAU = GET_TAU0( THISMONTH, 1, 1990 )
	
      ! Read 1990 NH3 emissions [kg N/box/mon]
      CALL READ_BPCH2( FILENAME, 'NH3-NATU',    29,  
     &                 XTAU,      IIPAR,        JJPAR,       
     &                 1,         ARRAY(:,:,1), QUIET=.TRUE. )

      ! Cast from REAL*4 to REAL*8
      CALL TRANSFER_2D( ARRAY(:,:,1), ENH3_na )

      ! Convert from [kg N/box/mon] to [kg NH3/box/s]
      ENH3_na = ENH3_na * ( 17.d0 / 14.d0 ) /
     &          ( NMDAY(THISMONTH) * 86400.d0 ) 
 
      END SUBROUTINE READ_NATURAL_NH3
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_biofuel_nh3
!
! !DESCRIPTION: Subroutine READ\_BIOFUEL\_NH3 reads the monthly mean biomass
!  NH3 and biofuel emissions from disk and converts to [kg NH3/box/s]. 
!  (rjp, bdf, bmy, phs, 9/20/02, 12/23/08)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE READ_BIOFUEL_NH3( THISMONTH ) 
!
! !USES:
!
      USE BIOMASS_MOD,          ONLY : BIOMASS
      USE BPCH2_MOD,            ONLY : GET_NAME_EXT_2D, GET_RES_EXT
      USE BPCH2_MOD,            ONLY : GET_TAU0,        READ_BPCH2
      USE CMN_SIZE_MOD
      USE DIRECTORY_MOD,        ONLY : DATA_DIR
      USE FUTURE_EMISSIONS_MOD, ONLY : GET_FUTURE_SCALE_NH3bf
      USE LOGICAL_MOD,          ONLY : LBIOMASS,        LFUTURE
      USE TIME_MOD,             ONLY : ITS_A_LEAPYEAR
      USE TRACER_MOD,           ONLY : XNUMOL
      USE TRACERID_MOD,         ONLY : IDTNH3
      USE TRANSFER_MOD,         ONLY : TRANSFER_2D

!
! !INPUT PARAMETERS: 
!
      INTEGER, INTENT(IN)           :: THISMONTH ! Current month number (1-12)

!
! !REVISION HISTORY: 
!  (1 ) Renamed from NH3_READ to READ_BIOMASS_NH3.  Also updated comments,
!        made cosmetic changes.  Now reads in both biomass and biofuel
!        emissions. (rjp, bmy, 12/13/02)
!  (2 ) Now replace DXYP(J+J0) with routine GET_AREA_M2 of "grid_mod.f"
!        Now use function GET_YEAR from "time_mod.f".  Replace MONTH with 
!        THISMONTH when referencing the NMDAY variable.  Now call READ_BPCH2
!        with QUIET=.TRUE. (bmy, 3/27/03)
!  (3 ) If using interannual biomass emissions, substitute seasonal emissions
!        for years where internannual emissions do not exist.  Now also
!        reference GET_TAU from "time_mod.f" (bmy, 5/15/03)
!  (4 ) Now use ENCODE statement for PGI/F90 on Linux (bmy, 9/29/03)
!  (5 ) Changed cpp switch name from LINUX to LINUX_PGI (bmy, 12/2/03)
!  (6 ) Now references DATA_DIR from "directory_mod.f".  Now references LBBSEA
!        from "logical_mod.f".  Removed references to CMN and CMN_SETUP.
!        (bmy, 7/20/04)
!  (7 ) Now can read either seasonal or interannual biomass burning emissions.
!        Now references routines from both and "time_mod.f".  Now reads SO2 
!        biomass emissions directly rather than computing it by mole fraction 
!        from CO. (rjp, bmy, 1/11/05)
!  (8 ) Now read files from "sulfate_sim_200508/".  Now read data for both 
!        GCAP and GEOS grids. (bmy, 8/16/05)
!  (9 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (10) Now only read the biofuel, we have moved the biomass-reading code to 
!        "gc_biomass_mod.f" for compatibility with GFED2 biomass emissions
!        (bmy, 9/27/06)
!  (11) Prevent seg fault error when LBIOMASS=F (bmy, 11/3/06)
!  (12) Renamed READ_BIOFUEL_NH3, and move all biomass code to GET_BIOMASS_NH3
!        to account for several GFED2 products (yc, phs, 12/23/08)
!  (13) IDBNH3 is not used anymore (ccc, 01/29/10)   
!  22 Dec 2011 - M. Payer    - Added ProTeX headers
!  01 Mar 2012 - R. Yantosca - Removed reference to GET_AREA_CM2
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER                       :: I, J, THISYEAR
      REAL*4                        :: ARRAY(IIPAR,JJPAR,1)
      REAL*8                        :: XTAU, DMASS!, CONV
      REAL*8                        :: NMDAY(12) = (/31d0, 28d0, 31d0, 
     &                                               30d0, 31d0, 30d0, 
     &                                               31d0, 31d0, 30d0, 
     &                                               31d0, 30d0, 31d0/)
      CHARACTER(LEN=4  )            :: CYEAR
      CHARACTER(LEN=255)            :: FILENAME

      !=================================================================
      ! READ_BIOFUEL_NH3 begins here!
      !=================================================================

      !=================================================================
      ! Read NH3 biofuel emissions
      !=================================================================

      ! File name
      FILENAME = TRIM( DATA_DIR )                         // 
     &           'sulfate_sim_200508/NH3_biofuel.'        // 
     &           GET_NAME_EXT_2D() // '.' // GET_RES_EXT()
   
      ! Echo output
      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( '     - READ_BIOFUEL_NH3: Reading ', a )

      ! Get TAU0 value for 1998
      XTAU  = GET_TAU0( THISMONTH, 1, 1998 )

      ! Read NH3 biofuel data [kg NH3/box/month]
      CALL READ_BPCH2( FILENAME, 'BIOFSRCE',    29, 
     &                 XTAU,      IIPAR,        JJPAR,       
     &                 1,         ARRAY(:,:,1), QUIET=.TRUE. )
 
      ! Cast from REAL*4 to REAL*8 and resize if necesary
      CALL TRANSFER_2D( ARRAY(:,:,1), ENH3_bf )

      ! Store NH3 in ENH3_bf array [kg NH3/box/s]
      ENH3_bf = ENH3_bf / ( NMDAY(THISMONTH) * 86400.d0 )

      !=================================================================
      ! Convert units and apply IPCC future emissions (if necessary)
      !=================================================================

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J )

      ! Loop over latitudes
      DO J = 1, JJPAR
 
         ! Loop over longitudes
         DO I = 1, IIPAR

            ! Scale biofuel NH3 to IPCC future scenario (if necessary)
            IF ( LFUTURE ) THEN
               ENH3_bf(I,J) = ENH3_bf(I,J) * 
     &                        GET_FUTURE_SCALE_NH3bf( I, J )
            ENDIF
         ENDDO
      ENDDO
!$OMP END PARALLEL DO
      
      END SUBROUTINE READ_BIOFUEL_NH3
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_biomass_nh3
!
! !DESCRIPTION: Subroutine GET\_BIOMASS\_NH3 retrieve the monthly/8days/3hr
!  mean biomass NH3 (yc, phs, 12/23/08)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE GET_BIOMASS_NH3
!
! !USES:
!
      USE BIOMASS_MOD,          ONLY : BIOMASS
      USE CMN_SIZE_MOD
      USE GRID_MOD,             ONLY : GET_AREA_CM2
      USE TRACER_MOD,           ONLY : XNUMOL
      USE TRACERID_MOD,         ONLY : IDTNH3, IDBNH3
!
! !REVISION HISTORY: 
!  (1 ) Extracted from old module subroutine READ_BIOMASS_NH3
!        (yc, phs, 12/23/08)
!  (2 ) IDBNH3 is in tracerid_mod.f now (fp, 6/2009)
!  22 Dec 2011 - M. Payer    - Added ProTeX headers
!  01 Mar 2012 - R. Yantosca - Now use GET_AREA_CM2(I,J,L) from grid_mod.F90
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER                       :: I, J
      REAL*8                        :: CONV
      
      !=================================================================
      ! READ_BIOMASSBURN_NH3 begins here!
      !=================================================================
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, CONV )

      ! Loop over surface grid boxes
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Conversion factor for [cm2 * kg/molec]
         CONV = GET_AREA_CM2( I, J, 1 ) / XNUMOL(IDTNH3)

         ! Convert biomass NH3 from [molec NH3/cm2/s] -> [kg NH3/s]
         ! NOTE: Future scale is applied by this point (if necessary)
         ENH3_bb(I,J) = BIOMASS(I,J,IDBNH3) * CONV

      ENDDO
      ENDDO
!$OMP END PARALLEL DO
      
      END SUBROUTINE GET_BIOMASS_NH3
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_oxidant
!
! !DESCRIPTION: Subroutine READ\_OXIDANT reads in monthly mean H2O2 and O3
!  fields for the offline sulfate + aerosol simulation. (rjp, bdf, bmy, 
!  11/1/02, 10/3/05)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE READ_OXIDANT( MONTH )
!
! !USES:
!
      USE BPCH2_MOD,     ONLY : GET_NAME_EXT, GET_RES_EXT
      USE BPCH2_MOD,     ONLY : GET_TAU0,     READ_BPCH2
      USE CMN_SIZE_MOD
      USE DIRECTORY_MOD, ONLY : DATA_DIR
      USE TRANSFER_MOD,  ONLY : TRANSFER_3D_TROP
!
! !INPUT PARAMETERS: 
!
      INTEGER, INTENT(IN)    :: MONTH ! Emission timestep in minutes
!
! !REVISION HISTORY: 
!  (1 ) Now call READ_BPCH2 with QUIET=.TRUE. (bmy, 3/27/03)
!  (2 ) Now references DATA_DIR from "directory_mod.f" (bmy, 7/20/04)
!  (3 ) Now read files from "sulfate_sim_200508/offline/".  Now read data
!        for both GEOS and GCAP grids (bmy, 8/16/05)
!  (4 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  22 Dec 2011 - M. Payer    - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER             :: I, J, L, K      
      REAL*4              :: ARRAY(IIPAR,JJPAR,LLTROP)
      REAL*8              :: XTAU
      CHARACTER(LEN=255)  :: FILENAME

      !=================================================================
      ! READ_OXIDANT begins here !
      !
      ! Oxidant fields were computed for 1998 using coupled aerosol
      ! and gas chemistry GEOS-CHEM by Brendan Field (bdf, 5/23/02).  
      ! Bob Yantosca has regridded these fields to all GEOS-CHEM grids.  
      ! Data is saved from the surface to the tropopause. 
      !=================================================================

      ! Use generic year 1985
      XTAU = GET_TAU0( MONTH, 1, 1985 )

      !=================================================================
      ! Read monthly mean PH2O2 (from HO2 + HO2 = H2O2) [molec/cm3/s]
      !=================================================================
      FILENAME = TRIM( DATA_DIR )                      // 
     &           'sulfate_sim_200508/offline/PH2O2.'   // 
     &           GET_NAME_EXT() //  '.' // GET_RES_EXT()

      ! Echo filename
      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( '     - READ_OXIDANT: Reading ', a ) 

      ! Read data
      CALL READ_BPCH2( FILENAME, 'PORL-L=$', 5,     
      ! limit array 3d dimension to LLTROP_FIX, i.e, case of annual mean
      ! tropopause. This is backward compatibility with 
      ! offline data set.
     &     XTAU,        IIPAR,                    JJPAR,      
     &     LLTROP_FIX,  ARRAY(:,:,1:LLTROP_FIX),  QUIET=.TRUE. )
!     &                 XTAU,      IIPAR,     JJPAR,     
!     &                 LLTROP,    ARRAY,     QUIET=.TRUE. )

      ! Cast to REAL*8 and resize if necessary
      CALL TRANSFER_3D_TROP( ARRAY, PH2O2m )
            
      !=================================================================
      ! Read monthly mean O3 [v/v]
      !=================================================================
      FILENAME = TRIM( DATA_DIR )                      // 
     &           'sulfate_sim_200508/offline/O3.'      //
     &           GET_NAME_EXT() // '.' // GET_RES_EXT()

      ! Echo filename
      WRITE( 6, 100 ) TRIM( FILENAME )

      ! Read data
      ! limit array 3d dimension to LLTROP_FIX, i.e, case of annual mean
      ! tropopause. This is backward compatibility with 
      ! offline data set.
      CALL READ_BPCH2( FILENAME, 'IJ-AVG-$', 51,     
     &     XTAU,        IIPAR,                    JJPAR,      
     &     LLTROP_FIX,  ARRAY(:,:,1:LLTROP_FIX),  QUIET=.TRUE. )
!     &                 XTAU,      IIPAR,     JJPAR,     
!     &                 LLTROP,    ARRAY,     QUIET=.TRUE. )

      ! Cast to REAL*8 and resize if necessary
      CALL TRANSFER_3D_TROP( ARRAY, O3m ) 
      
      END SUBROUTINE READ_OXIDANT
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ohno3time
!
! !DESCRIPTION: Subroutine OHNO3TIME computes the sum of cosine of the solar
!  zenith angle over a 24 hour day, as well as the total length of daylight. 
!  This is needed to scale the offline OH and NO3 concentrations.
!  (rjp, bmy, 12/16/02, 3/30/04)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE OHNO3TIME
!
! !USES:
!
      USE CMN_SIZE_MOD
      USE CMN_GCTM_MOD
      USE GRID_MOD, ONLY : GET_XMID,    GET_YMID_R
      USE TIME_MOD, ONLY : GET_NHMSb,   GET_ELAPSED_SEC
      USE TIME_MOD, ONLY : GET_TS_CHEM, GET_DAY_OF_YEAR, GET_GMT
!
! !REVISION HISTORY: 
!  (1 ) Copy code from COSSZA directly for now, so that we don't get NaN
!        values.  Figure this out later (rjp, bmy, 1/10/03)
!  (2 ) Now replace XMID(I) with routine GET_XMID from "grid_mod.f".  
!        Now replace RLAT(J) with routine GET_YMID_R from "grid_mod.f". 
!        Removed NTIME, NHMSb from the arg list.  Now use GET_NHMSb,
!        GET_ELAPSED_SEC, GET_TS_CHEM, GET_DAY_OF_YEAR, GET_GMT from 
!        "time_mod.f". (bmy, 3/27/03)
!  (3 ) Now store the peak SUNCOS value for each surface grid box (I,J) in 
!        the COSZM array. (rjp, bmy, 3/30/04)
!  22 Dec 2011 - M. Payer    - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
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

            ! IJLOOP is the 1-D loop index for SUNCOS
            IJLOOP = 0

            ! Loop over surface grid boxes
            DO J = 1, JJPAR
            DO I = 1, IIPAR

               ! Grid box latitude center [radians]
               YMID_R = GET_YMID_R( I, J, 1 )

               ! Increment IJLOOP
               IJLOOP = IJLOOP + 1
               TIMLOC = real(LHR0) + real(NT)/3600.0 + 
     &                  GET_XMID( I, J, 1 )/15.0
         
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

            !### Debug
            !PRINT*, '### IN OHNO3TIME'
            !PRINT*, '### N       : ', N
            !PRINT*, '### NDYSTEP : ', NDYSTEP
            !PRINT*, '### NT      : ', NT
            !PRINT*, '### JDAY    : ', JDAY
            !PRINT*, '### RLAT    : ', RLAT
            !PRINT*, '### XMID    : ', XMID
            !PRINT*, '### SUNTMP  : ', SUNTMP
            !PRINT*, '### TCOSZ   : ', MINVAL( TCOSZ ), MAXVAL( TCOSZ )
            !PRINT*, '### TTDAY   : ', MINVAL( TCOSZ ), MAXVAL( TCOSZ )

            ! Increment elapsed time [sec]
            NT = NT + ( GET_TS_CHEM() * 60 )             
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
! !IROUTINE: init_sulfate
!
! !DESCRIPTION: Subroutine INIT\_SULFATE initializes and zeros all allocatable
!  arrays declared in "sulfate\_mod.f" (bmy, 6/2/00, 10/15/09)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE INIT_SULFATE( am_I_Root, Input_Opt, RC )
!
! !USES:
!
      USE CMN_SIZE_MOD
      USE DRYDEP_MOD,         ONLY : DEPNAME, NUMDEP
      USE ERROR_MOD,          ONLY : ALLOC_ERR
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)  :: am_I_Root   ! Are we on the root CPU?
      TYPE(OptInput), INTENT(IN)  :: Input_Opt   ! Input Options object
!
! !INPUT/OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT) :: RC          ! Success or failure?
!
! !REVISION HISTORY: 
!  (1 ) Only allocate some arrays for the standalone simulation (NSRCX==10).
!        Also reference NSRCX from "CMN".  Now eferences routine ALLOC_ERR
!        from "error_mod.f" ((rjp, bdf, bmy, 10/15/02)
!  (2 ) Now also allocate the IJSURF array to keep the 1-D grid box indices
!        for SUNCOS (for both coupled & offline runs).  Now allocate PH2O2m 
!        and O3m for offline runs.  Also allocate ESO2_bf (bmy, 1/16/03)
!  (3 ) Now allocate ENH3_na array (rjp, bmy, 3/23/03)
!  (4 ) Now allocate COSZM for offline runs (bmy, 3/30/04)
!  (5 ) Now allocate ESO2_sh array (bec, bmy, 5/20/04)
!  (6 ) Now allocates ITS_AN_AEROSOL_SIM from "tracer_mod.f".  Now remove 
!        IJSURF (bmy, 7/20/04)
!  (7 ) Now locate species in the DEPSAV array here instead of in CHEMSULFATE.
!        Now reference LDRYD from "logical_mod.f".  Updated for AS, AHS, LET, 
!        SO4aq, NH4aq. (bmy, 1/6/06)
!  (8 ) Now allocates PSO4_ss, PNITs (bec, bmy, 4/13/05)
!  (9 ) Initialize drydep flags outside of IF block (bmy, 5/23/06)
!  (10) Now redimension EEV & NEV arrays for new SO2 volcanic emissions
!        inventory.  Deleted obsolete arrays from older SO2 volcanic
!        emissions inventory. (jaf, bmy, 10/15/09)
!  (11) Now alllocate PSO4_SO2AQ (win, 1/25/10)
!  22 Dec 2011 - M. Payer    - Added ProTeX headers
!  04 Mar 2013 - R. Yantosca - Now accept am_I_Root, Input_Opt, RC arguments
!  05 Mar 2013 - R. Yantosca - Now use Input_Opt%ITS_AN_AEROSOL_SIM
!  30 May 2013 - S. Farina   - Allocate PSO4_SO2AQ for TOMAS
!  26 Sep 2013 - R. Yantosca - Renamed GEOS_57 Cpp switch to GEOS_FP
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER            :: AS, I, J, N

      !=================================================================
      ! INIT_SULFATE begins here!
      !=================================================================

      ! Assume success
      RC = GIGC_SUCCESS

      ! Allocate arrays
      ALLOCATE( SSTEMP( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'SSTEMP' )
      SSTEMP = 0d0

      ALLOCATE( DMSo( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'DMSo' )
      DMSo = 0d0

      ALLOCATE( EEV( IIPAR, JJPAR, LVOLC ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EEV' )
      EEV = 0d0

      ALLOCATE( ENV( IIPAR, JJPAR, LVOLC ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'ENV' )
      ENV = 0d0

      ALLOCATE( ENH3_an( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'ENH3_an' )
      ENH3_an = 0d0

      ALLOCATE( ENH3_bb( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'ENH3_bb' )
      ENH3_bb = 0d0

      ALLOCATE( ENH3_bf( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'ENH3_bf' )
      ENH3_bf = 0d0

      ALLOCATE( ENH3_na( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'ENH3_na' )
      ENH3_na = 0d0

      ALLOCATE( ESO2_ac( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'ESO2_ac' )
      ESO2_ac = 0d0

      ALLOCATE( ESO2_an( IIPAR, JJPAR, NOXLEVELS ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'ESO2_an' )
      ESO2_an = 0d0

      ALLOCATE( ESO2_bb( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'ESO2_bb' )
      ESO2_bb = 0d0

      ALLOCATE( ESO2_bf( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'ESO2_bf' )
      ESO2_bf = 0d0

      ALLOCATE( ESO2_ev( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'ESO2_ev' )
      ESO2_ev = 0d0

      ALLOCATE( ESO2_nv( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'ESO2_nv' )
      ESO2_nv = 0d0

      ALLOCATE( ESO2_sh( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'ESO2_sh' )
      ESO2_sh = 0d0

      ALLOCATE( ESO4_an( IIPAR, JJPAR, NOXLEVELS ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'ESO4_an' )
      ESO4_an = 0d0

      ALLOCATE( PMSA_DMS( IIPAR, JJPAR, LLTROP ), STAT=AS ) 
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'PMSA_DMS' )
      PMSA_DMS = 0d0

      ALLOCATE( PSO2_DMS( IIPAR, JJPAR, LLTROP ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'PSO2_DMS' )
      PSO2_DMS = 0d0

      ALLOCATE( PSO4_SO2( IIPAR, JJPAR, LLTROP ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'PSO4_SO2' )
      PSO4_SO2 = 0d0

      ALLOCATE( PSO4_ss( IIPAR, JJPAR, LLTROP ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'PSO4_ss' )
      PSO4_ss = 0d0

      ALLOCATE( PNITs( IIPAR, JJPAR, LLTROP ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'PNITs' )
      PNITs = 0d0

      ALLOCATE( SOx_SCALE( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'SOx_SCALE' )
      SOx_SCALE = 0d0

#if   !defined( GEOS_5 ) && !defined( MERRA ) && !defined( GEOS_FP )
      ! If we are using GEOS-5 or MERRA met, then get the cloud fraction 
      ! directly from the met fields.  (skim, bmy, 1/14/10)
      ALLOCATE( VCLDF( IIPAR, JJPAR, LLTROP ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'VCLDF' )
      VCLDF = 0d0
#endif

#if   defined( TOMAS )
      ! Allocate for TOMAS microphysics
      ALLOCATE( PSO4_SO2AQ( IIPAR, JJPAR, LLTROP ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'PSO4_SO2AQ' )
      PSO4_SO2AQ = 0d0
#endif

      !=================================================================
      ! Only initialize the following for offline aerosol simulations
      !=================================================================
      IF ( Input_Opt%ITS_AN_AEROSOL_SIM ) THEN

         ALLOCATE( PH2O2m( IIPAR, JJPAR, LLTROP ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'PH2O2m' )
         PH2O2m = 0d0

         ALLOCATE( O3m( IIPAR, JJPAR, LLTROP ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'O3m' )
         O3m = 0d0

         ALLOCATE( JH2O2( IIPAR, JJPAR, LLTROP ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'JH2O' )
         JH2O2 = 0d0

         ALLOCATE( TCOSZ( IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'TCOSZ' )
         TCOSZ = 0d0

         ALLOCATE( TTDAY( IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'TTDAY' )
         TTDAY = 0d0
 
         ALLOCATE( COSZM( IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'COSZM' )
         COSZM = 0d0
      ENDIF

      !================================================================
      ! Find drydep species in the DEPSAV array
      !=================================================================

      ! Initialize flags
      DRYH2O2  = 0
      DRYSO2   = 0
      DRYSO4   = 0
      DRYSO4s  = 0
      DRYMSA   = 0
      DRYNH3   = 0
      DRYNH4   = 0
      DRYNIT   = 0
      DRYSO4s  = 0
      DRYAS    = 0
      DRYAHS   = 0
      DRYLET   = 0
      DRYSO4aq = 0
      DRYNH4aq = 0  

      ! Locate position of each tracer in DEPSAV
      IF ( Input_Opt%LDRYD ) THEN
         DO N = 1, Input_Opt%NUMDEP
            SELECT CASE ( TRIM( Input_Opt%DEPNAME(N) ) )
               CASE ( 'H2O2'   )
                  DRYH2O2  = N
               CASE ( 'SO2'   )
                  DRYSO2   = N
               CASE ( 'SO4'   )
                  DRYSO4   = N
               CASE ( 'SO4S'   )
                  DRYSO4s  = N
               CASE ( 'MSA'   )
                  DRYMSA   = N
               CASE ( 'NH3'   )
                  DRYNH3   = N
               CASE ( 'NH4'   )
                  DRYNH4   = N
               CASE ( 'NIT'   )
                  DRYNIT   = N
               CASE ( 'NITS'   )
                  DRYNITs  = N
               CASE ( 'AS'    )
                  DRYAS    = N
               CASE ( 'AHS'   )
                  DRYAHS   = N
               CASE ( 'LET'   )
                  DRYLET   = N
               CASE ( 'SO4aq' )
                  DRYSO4aq = N
               CASE ( 'NH4aq' )
                  DRYNH4aq = N
               CASE DEFAULT
                  ! Nothing
            END SELECT        
         ENDDO
      ENDIF

      END SUBROUTINE INIT_SULFATE
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cleanup_sulfate
!
! !DESCRIPTION: Subroutine CLEANUP\_SULFATE deallocates all previously 
!  allocated arrays for sulfate emissions -- call at the end of the run
!  (bmy, 6/1/00, 10/15/09)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CLEANUP_SULFATE  
! 
! !REVISION HISTORY: 
!  (1 ) Now also deallocates IJSURF. (bmy, 11/12/02)
!  (2 ) Now also deallocates ENH3_na (rjp, bmy, 3/23/03)
!  (3 ) Now also deallocates COSZM (rjp, bmy, 3/30/04)
!  (4 ) Now also deallocates ESO4_sh (bec, bmy, 5/20/04)
!  (5 ) Now remove IJSURF (bmy, 7/20/04)
!  (6 ) Bug fix: now deallocate PSO4_ss, PNITs (bmy, 5/3/06)
!  (7 ) Deleted obsolete arrays from older SO2 volcanic emissions 
!        inventory (jaf, bmy, 10/15/09)
!  (8 ) Deallocate PSO4_SO2AQ (win, 1/25/10)
!  22 Dec 2011 - M. Payer    - Added ProTeX headers
!  30 May 2013 - S. Farina   - Deallocate PSO4_SO2AQ for TOMAS
!EOP
!------------------------------------------------------------------------------
!BOC
!

      !=================================================================
      ! CLEANUP_SULFATE begins here!
      !=================================================================
      IF ( ALLOCATED( DMSo       ) ) DEALLOCATE( DMSo       )
      IF ( ALLOCATED( EEV        ) ) DEALLOCATE( EEV        )
      IF ( ALLOCATED( ENV        ) ) DEALLOCATE( ENV        )
      IF ( ALLOCATED( ENH3_an    ) ) DEALLOCATE( ENH3_an    )
      IF ( ALLOCATED( ENH3_bb    ) ) DEALLOCATE( ENH3_bb    )
      IF ( ALLOCATED( ENH3_bf    ) ) DEALLOCATE( ENH3_bf    )
      IF ( ALLOCATED( ENH3_na    ) ) DEALLOCATE( ENH3_na    )
      IF ( ALLOCATED( ESO2_ac    ) ) DEALLOCATE( ESO2_ac    )
      IF ( ALLOCATED( ESO2_an    ) ) DEALLOCATE( ESO2_an    )
      IF ( ALLOCATED( ESO2_nv    ) ) DEALLOCATE( ESO2_nv    )
      IF ( ALLOCATED( ESO2_ev    ) ) DEALLOCATE( ESO2_ev    )
      IF ( ALLOCATED( ESO2_bb    ) ) DEALLOCATE( ESO2_bb    )
      IF ( ALLOCATED( ESO2_bf    ) ) DEALLOCATE( ESO2_bf    )
      IF ( ALLOCATED( ESO2_sh    ) ) DEALLOCATE( ESO2_sh    )
      IF ( ALLOCATED( ESO4_an    ) ) DEALLOCATE( ESO4_an    )
      IF ( ALLOCATED( JH2O2      ) ) DEALLOCATE( JH2O2      )
      IF ( ALLOCATED( O3m        ) ) DEALLOCATE( O3m        )
      IF ( ALLOCATED( PH2O2m     ) ) DEALLOCATE( PH2O2m     )
      IF ( ALLOCATED( PMSA_DMS   ) ) DEALLOCATE( PMSA_DMS   )
      IF ( ALLOCATED( PNITs      ) ) DEALLOCATE( PNITs      )
      IF ( ALLOCATED( PSO2_DMS   ) ) DEALLOCATE( PSO2_DMS   )
      IF ( ALLOCATED( PSO4_SO2   ) ) DEALLOCATE( PSO4_SO2   )
#if   defined( TOMAS )
      IF ( ALLOCATED( PSO4_SO2AQ ) ) DEALLOCATE( PSO4_SO2AQ ) 
#endif
      IF ( ALLOCATED( PSO4_ss    ) ) DEALLOCATE( PSO4_ss    )
      IF ( ALLOCATED( SOx_SCALE  ) ) DEALLOCATE( SOx_SCALE  )
      IF ( ALLOCATED( SSTEMP     ) ) DEALLOCATE( SSTEMP     )
      IF ( ALLOCATED( TCOSZ      ) ) DEALLOCATE( TCOSZ      )
      IF ( ALLOCATED( TTDAY      ) ) DEALLOCATE( TTDAY      )          
      IF ( ALLOCATED( VCLDF      ) ) DEALLOCATE( VCLDF      )
      IF ( ALLOCATED( COSZM      ) ) DEALLOCATE( COSZM      )
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%% NOTE: This is for the SEAC4RS volcano emissions.  Leave code commented
!%%% out here so that we can restore it later. (bmy, 9/25/13)
!      IF ( ALLOCATED( VEV_ELEV  ) ) DEALLOCATE( VEV_ELEV  )
!      IF ( ALLOCATED( VEV_SO2   ) ) DEALLOCATE( VEV_SO2   )
!      IF ( ALLOCATED( VEV_DATE  ) ) DEALLOCATE( VEV_DATE  )
!      IF ( ALLOCATED( VEV_II    ) ) DEALLOCATE( VEV_II    )
!      IF ( ALLOCATED( VEV_JJ    ) ) DEALLOCATE( VEV_JJ    )
!      IF ( ALLOCATED( VNV_ELEV  ) ) DEALLOCATE( VNV_ELEV  )
!      IF ( ALLOCATED( VNV_SO2   ) ) DEALLOCATE( VNV_SO2   )
!      IF ( ALLOCATED( VNV_DATE  ) ) DEALLOCATE( VNV_DATE  )
!      IF ( ALLOCATED( VNV_II    ) ) DEALLOCATE( VNV_II    )
!      IF ( ALLOCATED( VNV_JJ    ) ) DEALLOCATE( VNV_JJ    )
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      END SUBROUTINE CLEANUP_SULFATE
!EOC
      END MODULE SULFATE_MOD
