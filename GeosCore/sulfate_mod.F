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
      USE HCO_ERROR_MOD    ! For HEMCO error handling
      USE PRECISION_MOD    ! For GEOS-Chem Precision (fp, f4, f8)

      IMPLICIT NONE
      PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS: 
!
      PUBLIC :: CHEMSULFATE       
      PUBLIC :: CLEANUP_SULFATE   
      PUBLIC :: INIT_SULFATE
#if defined( TOMAS )
      PUBLIC :: EMISSSULFATETOMAS ! JKodros (6/2/15 - this is to connect TOMAS to HEMCO)
#endif
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
!  20 Aug 2013 - R. Yantosca   - Removed "define.h", this is now obsolete
!  18 Sep 2014 - M. Sulprizio  - Get oxidant fields for offline aerosol
!                                simulation from HEMCO 
!  03 Nov 2014 - C. Keller     - Incorporated GET_ALK from seasalt_mod.F
!  20 Nov 2014 - M. Yannetti   - Added PRECISION_MOD
!  04 Mar 2015 - R. Yantosca   - Remove obsolete, commented-out code
!  04 Mar 2015 - R. Yantosca   - Use REAL(f4) for pointer args to HCO_GetPtr
!  22 May 2015 - R. Yantosca   - Remove variables made obsolete by HEMCO
!  12 Jun 2015 - R. Yantosca   - Now remove orphaned ND44 variables
!  12 Jun 2015 - R. Yantosca   - Remove CHEM_NH3 and CHEM_NH4 routines
!                                because drydep is now done in mixing_mod.F90
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
      REAL(fp),  PARAMETER   :: XNUMOL_OH   = 6.022e+23_fp / 17e-3_fp
      REAL(fp),  PARAMETER   :: XNUMOL_O3   = 6.022e+23_fp / 48e-3_fp
      REAL(fp),  PARAMETER   :: XNUMOL_NO3  = 6.022e+23_fp / 62e-3_fp
      REAL(fp),  PARAMETER   :: XNUMOL_H2O2 = 6.022e+23_fp / 34e-3_fp
      REAL(fp),  PARAMETER   :: TCVV_S      = 28.97e+0_fp  / 32e+0_fp 
      REAL(fp),  PARAMETER   :: SMALLNUM    = 1e-20_fp
      REAL(fp),  PARAMETER   :: CM3PERM3    = 1.e6_fp

#if defined( TOMAS )
      !---------------------------------------------------------------
      ! For TOMAS microphysics: Add parameter for scaling anthro SO2
      !---------------------------------------------------------------
      REAL(fp),  PARAMETER   :: scaleanthso2 = 1.0e+0_fp
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
      !
      !%%% NOTE: THESE ARE NOW OBTAINED VIA HEMCO (bmy, 5/22/15) %%%%%%%%%%%
      !% ENH3_an    : NH3 anthropogenic emissions      [kg NH3/box/s]
      !% ENH3_bb    : NH3 biomass emissions            [kg NH3/box/s]
      !% ENH3_bf    : NH3 biofuel emissions            [kg NH3/box/s]
      !% ENH3_na    : NH73 natural source emissions    [kg NH3/box/s]
      !% ESO2_ac    : SO2 aircraft emissions           [kg SO2/box/s]
      !% ESO2_an    : SO2 anthropogenic emissions      [kg SO2/box/s]
      !% ESO2_ev    : SO2 eruptive volcanic em.        [kg SO2/box/s]
      !% ESO2_nv    : SO2 non-eruptive volcanic em.    [kg SO2/box/s]
      !% ESO2_bb    : SO2 biomass burning emissions    [kg SO2/box/s]
      !% ESO2_bf    : SO2 biofuel burning emissions    [kg SO2/box/s]
      !% ESO2_sh    : SO2 ship emissions               [kg SO2/box/s]
      !% ESO4_an    : SO4 anthropogenic emissions      [kg SO4/box/s]
      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      !
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
      REAL(fp),  ALLOCATABLE :: DMSo(:,:) 
      REAL(fp),  ALLOCATABLE :: PMSA_DMS(:,:,:)
      REAL(fp),  ALLOCATABLE :: PSO2_DMS(:,:,:)
      REAL(fp),  ALLOCATABLE :: PSO4_SO2(:,:,:)
      REAL(fp),  ALLOCATABLE :: PSO4_SS(:,:,:)
      REAL(fp),  ALLOCATABLE :: PNITs(:,:,:)
      REAL(f4),  ALLOCATABLE :: SOx_SCALE(:,:)
      REAL(fp),  ALLOCATABLE :: SSTEMP(:,:)
      REAL(fp),  ALLOCATABLE :: TCOSZ(:,:)
      REAL(fp),  ALLOCATABLE :: TTDAY(:,:)
      REAL(fp),  ALLOCATABLE :: VCLDF(:,:,:)
      REAL(fp),  ALLOCATABLE :: COSZM(:,:)

#if   defined( TOMAS )
      !---------------------------------------------------------------
      ! For TOMAS microphysics: Define PSO4_SO2aq array
      !---------------------------------------------------------------
      REAL(fp),  ALLOCATABLE :: PSO4_SO2AQ(:,:,:)
      REAL(fp),  ALLOCATABLE :: SO4_ANTH(:,:,:)
      REAL(fp),  ALLOCATABLE :: SO4_BIOF(:,:)
#endif

      ! Pointers to drydep species w/in DEPSAV
      INTEGER                :: DRYSO2,  DRYSO4,   DRYMSA,  DRYNH3  
      INTEGER                :: DRYNH4,  DRYNIT,   DRYSO4s, DRYNITs
      INTEGER                :: DRYH2O2, DRYSO4aq, DRYAS,   DRYAHS
      INTEGER                :: DRYLET,  DRYNH4aq

      ! These are pointers to fields in the HEMCO data structure.
      ! Declare these with REAL(fp), aka REAL*4. (bmy, 3/4/15)
      REAL(f4), POINTER      :: O3m(:,:,:)      => NULL()
      REAL(f4), POINTER      :: PH2O2m(:,:,:)   => NULL()
      REAL(f4), POINTER      :: JH2O2(:,:,:)    => NULL()
      REAL(f4), POINTER      :: OH(:,:,:)       => NULL()
      REAL(f4), POINTER      :: NO3(:,:,:)      => NULL()
      REAL(f4), POINTER      :: HNO3(:,:,:)     => NULL()
      REAL(f4), POINTER      :: NDENS_SALA(:,:) => NULL()
      REAL(f4), POINTER      :: NDENS_SALC(:,:) => NULL()

      ! Emission timestep (imported from HEMCO)
      REAL(fp)               :: TS_EMIS

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
      REAL(fp),  PARAMETER   :: ZRT = 0.60e+0_fp, ZRS = 0.99e+0_fp
!
! !LOCAL VARIABLES:
!
      INTEGER              :: I,    J,    L
      REAL(fp)               :: PRES, PSFC, RH2, R0, B0
	
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
      DO L = 1, LLCHEM
      DO J = 1, JJPAR 
      DO I = 1, IIPAR
	
         ! Surface pressure
         PSFC = GET_PEDGE(I,J,1)

         ! Pressure at the center of the grid box
         PRES = GET_PCENTER(I,J,L)

         ! RH (from "dao_mod.f") is relative humidity [%]
         ! Convert to fraction and store in RH2
         RH2  = State_Met%RH(I,J,L) * 1.0e-2_fp

         ! Terms from Sundqvist ???
         R0   = ZRT + ( ZRS - ZRT ) * EXP( 1e+0_fp 
     &      - ( PSFC / PRES )**2.5 )
         B0   = ( RH2 - R0 ) / ( 1e+0_fp - R0 )
	   
         ! Force B0 into the range 0-1
         IF ( RH2 < R0  ) B0 = 0e+0_fp
         IF ( B0  > 1e+0_fp ) B0 = 1e+0_fp

         ! Volume cloud fraction
         VCLDF(I,J,L) = 1e+0_fp - SQRT( 1e+0_fp - B0 )

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
      REAL(fp), INTENT(IN) :: T ! Temperature value at a GEOS-CHEM grid box [K]
!
! !RETURN VALUE:
!
      REAL(fp)             :: LWC
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
      IF ( T > 293e+0_fp ) THEN
         LWC = 0.2e+0_fp

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
      USE ERROR_MOD,          ONLY : DEBUG_MSG
      USE ERROR_MOD,          ONLY : ERROR_STOP
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE GIGC_State_Chm_Mod, ONLY : ChmState
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE HCO_EMISLIST_MOD,   ONLY : HCO_GetPtr
      USE HCOI_GC_MAIN_MOD,   ONLY : GetHcoState
      USE HCO_STATE_MOD,      ONLY : HCO_STATE
      USE TIME_MOD,           ONLY : GET_MONTH
      USE TIME_MOD,           ONLY : GET_TS_CHEM
      USE TIME_MOD,           ONLY : GET_ELAPSED_SEC
      USE TIME_MOD,           ONLY : ITS_A_NEW_MONTH
      USE TRACERID_MOD,       ONLY : IDTNITs
      USE TRACERID_MOD,       ONLY : IDTSO4s
      USE UCX_MOD,            ONLY : SETTLE_STRAT_AER
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
!  18 Sep 2014 - M. Sulprizio- Get oxidant fields for offline aerosol
!                              simulation from HEMCO
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
      REAL(fp)            :: DTCHEM

      ! For fields from Input_Opt
      LOGICAL           :: LCRYST
      LOGICAL           :: LPRT
      LOGICAL           :: LGRAVSTRAT
      INTEGER           :: N_TRACERS
      REAL(fp)          :: TCVV(Input_Opt%N_TRACERS)

      ! Pointers
      ! We need to define local arrays to hold corresponding values 
      ! from the Chemistry State (State_Chm) object. (mpayer, 12/6/12)
      REAL(fp), POINTER  :: STT(:,:,:,:)

      ! For getting fields from HEMCO 
      LOGICAL            :: aIR
      CHARACTER(LEN=255) :: LOC = 'CHEMSULFATE (sulfate_mod.F)'

      ! To get emission timestep
      TYPE(HCO_STATE), POINTER :: HcoState => NULL()

      !=================================================================
      ! CHEMSULFATE begins here!
      !=================================================================

      ! Assume success
      RC                   = GIGC_SUCCESS

      ! am I root? 
      aIR                  = am_I_Root

      ! Copy fields from INPUT_OPT to local variables for use below
      LCRYST               = Input_Opt%LCRYST
      LPRT                 = Input_Opt%LPRT
      LGRAVSTRAT           = Input_Opt%LGRAVSTRAT
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

         ! Get offline oxidant fields from HEMCO (mps, 9/18/14)
         IF ( FIRSTCHEM ) THEN

            CALL HCO_GetPtr( aIR, 'O3',    O3m,    RC )
            IF ( RC /= GIGC_SUCCESS )
     &      CALL ERROR_STOP( 'Cannot get pointer to O3',    LOC )

            CALL HCO_GetPtr( aIR, 'PH2O2', PH2O2m, RC )
            IF ( RC /= GIGC_SUCCESS )
     &      CALL ERROR_STOP( 'Cannot get pointer to PH2O2', LOC )

            CALL HCO_GetPtr( aIR, 'JH2O2', JH2O2,  RC )
            IF ( RC /= GIGC_SUCCESS )
     &      CALL ERROR_STOP( 'Cannot get pointer to JH2O2', LOC )

            CALL HCO_GetPtr( aIR, 'GLOBAL_OH',  OH,    RC )
            IF ( RC /= GIGC_SUCCESS )
     &      CALL ERROR_STOP( 'Cannot get pointer to GLOBAL_OH',   LOC )

            CALL HCO_GetPtr( aIR, 'GLOBAL_NO3', NO3,   RC )
            IF ( RC /= GIGC_SUCCESS )
     &      CALL ERROR_STOP( 'Cannot get pointer to GLOBAL_NO3',  LOC )

            CALL HCO_GetPtr( aIR, 'GLOBAL_HNO3', HNO3, RC )
            IF ( RC /= GIGC_SUCCESS )
     &      CALL ERROR_STOP( 'Cannot get pointer to GLOBAL_HNO3', LOC )

         ENDIF

         ! And compute time scaling arrays for offline OH, NO3
         CALL OHNO3TIME
         
      ENDIF

      ! Store NTIME in a shadow variable
      ELAPSED_SEC = GET_ELAPSED_SEC()

      ! DTCHEM is the chemistry timestep in seconds
      DTCHEM = GET_TS_CHEM() * 60e+0_fp

      ! TS_EMIS is the emission timestep (in seconds). This is a module
      ! variable, hence define only on first call.
      IF ( FIRSTCHEM ) THEN
         CALL GetHcoState ( HcoState )
         IF ( .NOT. ASSOCIATED(HcoState) ) 
     &        CALL ERROR_STOP ( 'Cannot get HcoState', LOC )
         TS_EMIS = HcoState%TS_EMIS
         HcoState => NULL()
      ENDIF

      ! Initialize module arrays
      PSO2_DMS = 0e+0_fp
      PMSA_DMS = 0e+0_fp
      PSO4_SO2 = 0e+0_fp
      PSO4_SS  = 0e+0_fp
      PNITs    = 0e+0_fp
#if   defined( TOMAS )
      PSO4_SO2AQ = 0e+0_fp     ! For TOMAS microphysics
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

#if defined( UCX )
      ! Stratospheric aerosol graviational settling
      IF ( LGRAVSTRAT ) THEN
         CALL SETTLE_STRAT_AER( am_I_Root, Input_Opt, 
     &                          State_Met, State_Chm, RC )
         IF ( prtDebug ) THEN
            CALL DEBUG_MSG( '### CHEMSULFATE: GRAV_SET, STRAT' )
         ENDIF
      ENDIF
#endif

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

      ! Sulfur Nitrate.
      ! CHEM_NIT includes a source term from sea salt aerosols, so keep
      ! here.
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
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: emisssulfatetomas
!
! !DESCRIPTION: Subroutine EMISSSULFATETOMAS connects HEMCO bulk emissions to 
! the TOMAS tracers. Only use this for TOMAS sims. This should be quite similar 
! to the TOMAS relevant parts of 'emisssulfate' in v9 (Jkodros 6/2/15)
!
!\\
!\\
#if defined( TOMAS )
      SUBROUTINE EMISSSULFATETOMAS( am_I_Root, Input_Opt,
     &                         State_Met, State_Chm, RC )

! !USES:
      USE CMN_DIAG_MOD
      USE CMN_SIZE_MOD
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE GIGC_State_Chm_Mod, ONLY : ChmState
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE TRACERID_MOD,       ONLY : IDTSF1,     IDTAW1,  IDTNH4
      USE TRACERID_MOD,       ONLY : IDTNK1
      USE TOMAS_MOD,          ONLY : IBINS,      ICOMP,   IDIAG
      USE TOMAS_MOD,          ONLY : NH4BULKTOBIN
      USE TOMAS_MOD,          ONLY : SRTNH4
!#endif
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)    :: am_I_Root   ! Are we on the root CPU?
      TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
      TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology State object

! !INPUT/OUTPUT PARAMETERS:
!
      TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State objectt
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?

! Local variables
      ! Fields for TOMAS simulation
      REAL*8        :: BINMASS(IIPAR,JJPAR,LLPAR,IBINS*ICOMP)
      REAL*8        :: tempnh4(ibins)
      INTEGER       :: TID, I, J, L
      REAL*8        :: MK_TEMP(IBINS)
      INTEGER       :: ii=53, jj=29, ll=1

      ! Pointers
      ! We need to define local arrays to hold corresponding values
      ! from the Chemistry State (State_Chm) object. (mpayer, 12/6/12)
      REAL*8, POINTER :: STT(:,:,:,:)

      !=================================================================
      ! EMISSSULFATETOMAS begins here!
      !=================================================================
      ! Initialize GEOS-Chem tracer array [kg] from Chemistry State object
      STT => State_Chm%Tracers

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

               ! Take a slice of size IBINS from STT
               MK_TEMP = STT(I,J,L,IDTSF1:IDTSF1-1+IBINS)

               ! Avoid array tmeporaries in subroutine call (bmy, 1/29/14)
               CALL NH4BULKTOBIN( MK_TEMP, STT(I,J,L,IDTNH4), TEMPNH4 )

               BINMASS(I,J,L,TID:TID+IBINS-1) = TEMPNH4(1:IBINS)
            Enddo
            ENDDO
            ENDDO
!%%%OMP END PARALLEL DO

         ENDIF

         TID = IBINS*(ICOMP-1) +1
         BINMASS(:,:,:,TID:TID+IBINS-1) =
     &           STT(:,:,:,IDTAW1:IDTAW1+IBINS-1)

      !IF ( IDTSF1 > 0 ) THEN
         CALL SRCSF30(  STT(:,:,:,IDTNK1:IDTNK1+IBINS-1),
     &         BINMASS(:,:,:,:), am_I_Root, Input_Opt, State_Met, RC )

         ! Return the aerosol mass after emission subroutine to STT
         ! excluding the NH4 aerosol and aerosol water (win, 9/27/08)
         STT(:,:,:,IDTSF1:IDTSF1+IBINS*(ICOMP-IDIAG)-1) =
     &        BINMASS(:,:,:,1:IBINS*(ICOMP-IDIAG))
      ENDIF


      ! Free pointer
      NULLIFY( STT )

      END SUBROUTINE EMISSSULFATETOMAS
#endif

!
!-----------------------------------------------------------------------------
!                  Jack Kodros re-writing this
!-----------------------------------------------------------------------------
#if defined( TOMAS )
      SUBROUTINE SRCSF30( TC1, TC2,
     &     am_I_Root, Input_Opt, State_Met, RC )

      USE DIAG_MOD,           ONLY : AD59_SULF,     AD59_NUMB
      USE ERROR_MOD,          ONLY : ERROR_STOP,  IT_IS_NAN
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE PBL_MIX_MOD,        ONLY : GET_FRAC_OF_PBL, GET_PBL_TOP_L
      USE PRESSURE_MOD,       ONLY : GET_PCENTER

      USE TOMAS_MOD,          ONLY : IBINS, AVGMASS, ICOMP
      USE TOMAS_MOD,          ONLY : Xk,  SUBGRIDCOAG, MNFIX
      USE TOMAS_MOD,          ONLY : SRTSO4, SRTNH4,  DEBUGPRINT
      USE TRACERID_MOD

      USE CMN_SIZE_MOD             ! Size parameters
      USE CMN_DIAG_MOD             ! ND13 (for now)
      USE CMN_GCTM_MOD             ! SCALE_HEIGHT

      ! HEMCO update
      USE HCOI_GC_MAIN_MOD,   ONLY : GetHcoState, GetHcoDiagn
      USE HCO_STATE_MOD,      ONLY : HCO_STATE
      USE HCO_STATE_MOD,      ONLY : HCO_STATE

!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)    :: am_I_Root   ! Are we on the root CPU?
      TYPE(OptInput), INTENT(IN)  :: Input_Opt   ! Input Options object
      TYPE(MetState), INTENT(IN)  :: State_Met   ! Meteorology State object
! !INPUT/OUTPUT PARAMETERS:
!
      REAL*8,  INTENT(INOUT) :: TC1(IIPAR,JJPAR,LLPAR,IBINS)
      REAL*8,  INTENT(INOUT) :: TC2(IIPAR,JJPAR,LLPAR,IBINS*ICOMP)
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?

! Local variables
      INTEGER                :: I, J, K, L, DOW_LT, NTOP, C
      REAL*8                 :: SO4(LLPAR), DTSRCE  , EFRAC(LLPAR)
      REAL*8                 :: TSO4,       FEMIS
      REAL*8                 :: AREA_CM2
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

      REAL(fp)                 :: AREA(IIPAR, JJPAR)
      REAL(fp)                 :: AREA3D(IIPAR, JJPAR,2)

      ! Pointers
      TYPE(HCO_State), POINTER :: HcoState   => NULL()
      REAL(f4),        POINTER :: Ptr2D(:,:) => NULL()
      REAL(f4),        POINTER :: Ptr3D(:,:,:) => NULL()

      INTEGER                  :: N_TRACERS
      REAL(fp)                 :: XNUMOL(Input_Opt%N_TRACERS)

      LOGICAL                :: ERRORSWITCH, SGCOAG = .TRUE.
      INTEGER                  :: FLAG, ERR
      logical                :: pdbug !(temporary) win, 10/24/07
      !integer                :: ii, jj, ll
      !data ii, jj, ll / 61, 1, 7 /
      INTEGER                :: ii=53, jj=29, ll=1

      ! Ratio of molecular weights: S/SO4
      REAL*8,  PARAMETER     :: S_SO4 = 32d0 / 96d0

      ! debugging
      real*8   dummy

      ! For fields from Input_Opt
      LOGICAL :: LPRT, LNLPBL
      LOGICAL                 :: jkdbg=.true.

      ! Strings
      CHARACTER(LEN= 63)       :: DgnName
      CHARACTER(LEN=255)       :: MSG
      CHARACTER(LEN=255)       :: LOC='srcsf30 (sulfate_mod.F)'

      !=================================================================
      ! SRCSF30 begins here!
      !=================================================================

      ! COpy values from Input_Opt
      LPRT    = Input_Opt%LPRT
      LNLPBL  = Input_Opt%LNLPBL

      ! Import emissions from HEMCO (through HEMCO state)
      CALL GetHcoState( HcoState )
      IF ( .NOT. ASSOCIATED(HcoState) ) THEN
         CALL ERROR_STOP ( 'HcoState not defined!', LOC )
      ENDIF

      ! Emission timestep [seconds]
      DTSRCE = HcoState%TS_EMIS

      ! Grid box aarea
      AREA = HcoState%Grid%AREA_M2%Val(:,:)
      AREA3D(:,:,1) = AREA(:,:)
      AREA3D(:,:,2) = AREA(:,:)

      !================================================================
      ! READ IN HEMCO EMISSIONS
      !================================================================
      DgnName = 'SO4_ANTH'
      CALL GetHcoDiagn( am_I_Root, DgnName, .FALSE., ERR,
     &  Ptr3D=Ptr3D )
      IF ( .NOT. ASSOCIATED(Ptr3D) ) THEN
         CALL HCO_WARNING('Not found: '//TRIM(DgnName),ERR,THISLOC=LOC)
      ELSE
         SO4_ANTH = Ptr3D(:,:,:)
      ENDIF
      Ptr3D => NULL()

      DgnName = 'SO4_BIOF'
      CALL GetHcoDiagn( am_I_Root, DgnName, .FALSE., ERR,
     &  Ptr2D=Ptr2D )
      IF ( .NOT. ASSOCIATED(Ptr2D) ) THEN
         CALL HCO_WARNING('Not found: '//TRIM(DgnName),ERR,THISLOC=LOC)
      ELSE
         SO4_BIOF = Ptr2D(:,:)
      ENDIF
      Ptr2D => NULL()

      ! convert to kg/box/sec
      SO4an = SO4_ANTH(:,:,:) * AREA3d(:,:,:)
      SO4bf = SO4_BIOF(:,:) * AREA(:,:)

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
               !print*, 'JACK IN ND59 SULFATE'
               DO L = 1, LLPAR
               DO K = 1, IBINS
                  if(TC2(I,J,L,K)-M0(L,K) < 0d0)
     &              print *,'Negative SF emis ',TC2(I,J,L,K)-M0(L,K),
     &                 'at',I,J,L,K
                  if(TC1(I,J,L,K)-N0(L,K) < 0d0)
     &              print *,'Negative NK emis ',TC1(I,J,L,K)-N0(L,K),
     &                 'at',I,J,L,K
               !   if(I==ii .and. J==jj) then
               !      print*, 'TC2: ', TC2(ii,jj,0,:)
               !      print*, 'M0: ', M0(0,:)
               !      print*, 'S_SO4: ', S_SO4
               !   ENDIF
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
#endif

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
#if !defined( NO_BPCH )
      USE DIAG_MOD,           ONLY : AD44
#endif
!      USE DRYDEP_MOD,         ONLY : DEPSAV
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
      REAL(fp),         INTENT(INOUT) :: TC(IIPAR,JJPAR,LLPAR) ! Tracer [kg]
! 
! !REVISION HISTORY:
!  (1 ) Now references SALA_REDGE_um and SALC_REDGE_um from "tracer_mod.f"
!        (bmy, 7/20/04)
!  (2 ) Now references XNUMOL from "tracer_mod.f" (bmy, 10/25/05)
!  (3 ) Now limit relative humidity to [tiny(real(fp)),0.99] range for DLOG
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
!  06 Jan 2015 - M. Yannetti - Changed some variables to f8 as needed
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
      REAL(fp),  PARAMETER     :: C1   =  0.7674e+0_fp 
      REAL(fp),  PARAMETER     :: C2   =  3.079e+0_fp 
      REAL(fp),  PARAMETER     :: C3   =  2.573e-11_fp
      REAL(fp),  PARAMETER     :: C4   = -1.424e+0_fp
      REAL(fp),  PARAMETER     :: DEN  = 2200.0e+0_fp ! [kg/m3] sea-salt density

      ! Parameters for polynomial coefficients to derive seawater
      ! density. From Tang et al. (1997) (bec, jaegle, 5/11/11)
      REAL(fp),  PARAMETER     :: A1   =  7.93e-3_fp
      REAL(fp),  PARAMETER     :: A2   = -4.28e-5_fp
      REAL(fp),  PARAMETER     :: A3   =  2.52e-6_fp
      REAL(fp),  PARAMETER     :: A4   = -2.35e-8_fp
      REAL(f8),  PARAMETER     :: EPSI = 1.0e-4_f8 
!
! !LOCAL VARIABLES:
!
      INTEGER                :: I,      J,     L,        DTCHEM
      REAL(fp)                 :: DELZ,   DELZ1, REFF
      REAL(fp)                 :: P,      DP,    PDP,      TEMP        
      REAL(fp)                 :: CONST,  SLIP,  VISC,     FAC1
      REAL(fp)                 :: FAC2,   FLUX,  AREA_CM2, RHB
      ! replace RCM (radius in CM with RUM radius in microns) jaegle 5/11/11
      REAL(fp)                 :: RUM,    RWET,  RATIO_R
      REAL(fp)                 :: TOT1,   TOT2
      REAL(fp)                 :: VTS(LLPAR)  
      REAL(fp)                 :: TC0(LLPAR)
      ! added variables for density calculation (jaegle, bec 5/11/11)
      REAL(f8)                 :: RHO1, WTP, RHO
 
      ! Arrays
      INTEGER                :: IDDEP(2)
      INTEGER                :: IDTRC(2)	
      REAL(fp)                 :: SALA_REDGE_um(2)
      REAL(fp)                 :: SALC_REDGE_um(2)
      REAL(fp)                 :: XNUMOL(Input_Opt%N_TRACERS)

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
      DTCHEM = GET_TS_CHEM() * 60e+0_fp

      ! Store in IDDEP array
      IDDEP(1) = DRYSO4s
      IDDEP(2) = DRYNITs

      ! Tracer array
      IDTRC(1) = IDTSO4s
      IDTRC(2) = IDTNITs

      ! Coarse mode
      REFF = 0.5e-6_fp * ( SALC_REDGE_um(1) + SALC_REDGE_um(2) )
            
      ! Sea salt radius [cm]
      ! The Gerber formula for hygroscopic growth uses the radius in
      ! micrometers instead of centimeters. This fix is implemented by using
      ! RUM instead of RCM (jaegle 5/5/11)
      RUM  = REFF * 1e+6_fp


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
            VTS(L) = 0e+0_fp
         ENDDO

         ! Loop over levels
         DO L = 1, LLPAR

            ! Pressure at center of the level [kPa]
            P       = GET_PCENTER(I,J,L) * 0.1e+0_fp

            ! Temperature [K]
            TEMP    = State_Met%T(I,J,L)

            ! Cap RH at 0.99 
            RHB    = MIN( 0.99e+0_fp, State_Met%RH(I,J,L) * 1e-2_fp )

            ! Safety check (phs, 5/1/08)
            RHB     = MAX( TINY(RHB), RHB           )

            ! Aerosol growth with relative humidity in radius [m] 
            ! (Gerber, 1985)
            ! Several bug fixes to the Gerber formulation: a log10 (instead of
            ! ln) should be used and the dry radius should be expressed in
            ! micrometers (instead of cm) also add more significant digits to
            ! the exponent (jaegle 5/5/11)
            !RWET    = 1d-6*(FAC1/(FAC2-LOG10(RHB))+RUM**3.e+0_fp)**0.33333e+0_fp

            ! Use equation 5 in Lewis and Schwartz (2006) [m] for sea salt
            ! growth (jaegle 5/11/11)
            RWET = REFF * (4.e+0_fp / 3.7e+0_fp) *
     &        ( (2.e+0_fp - RHB)/(1.e+0_fp - RHB) )**(1.e+0_fp/3.e+0_fp)


            ! Ratio dry over wet radii at the cubic power
            RATIO_R = ( REFF / RWET )**3.e+0_fp

            ! Density of the wet aerosol (kg/m3)
            !RHO     = RATIO_R * DEN + ( 1.e+0_fp - RATIO_R ) * 1000.e+0_fp

            ! Above density calculation is chemically unsound because it
            ! ignores chemical solvation.  
            ! Iteratively solve Tang et al., 1997 equation 5 to calculate
            ! density of wet aerosol (kg/m3) 
            ! (bec, jaegle 5/11/11)
            RATIO_R = ( REFF / RWET )
            ! Assume an initial density of 1000 kg/m3
            RHO  = 1000.e+0_f8
            RHO1 = 0.e+0_f8 !initialize (bec, 6/21/10)
            DO WHILE ( ABS( RHO1-RHO ) .gt. EPSI )
                ! First calculate weight percent of aerosol (kg_RH=0.8/kg_wet) 
                WTP    = 100.e+0_f8 * DEN/RHO * RATIO_R**3.e+0_f8
                ! Then calculate density of wet aerosol using equation 5 
                ! in Tang et al., 1997 [kg/m3]
                RHO1   = ( 0.9971e+0_f8 + (A1 * WTP) 
     $          + (A2 * WTP**2.e+0_f8)
     $          + (A3 * WTP**3.e+0_f8) 
     $          + (A4 * WTP**4.e+0_f8) ) * 1000.e+0_f8
                ! Now calculate new weight percent using above density
                ! calculation
                WTP    = 100.e+0_f8 * DEN/RHO1 * RATIO_R**3.e+0_f8
                ! Now recalculate new wet density [kg/m3]
                RHO   = ( 0.9971e+0_f8 + (A1 * WTP) 
     $          + (A2 * WTP**2.e+0_f8)
     $          + (A3 * WTP**3.e+0_f8) 
     $          + (A4 * WTP**4.e+0_f8) ) * 1000.e+0_f8
            ENDDO

            ! Dp = particle diameter [um]
            DP      = 2.e+0_fp * RWET * 1.e+6_fp        

            ! PdP = P * dP [hPa * um]
            PDp     = P * Dp

            ! Constant
            CONST   = 2.e+0_fp * RHO * RWET**2 * g0 / 9.e+0_fp

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
            Slip = 1.e+0_fp+(15.60e+0_fp + 7.0e+0_fp * 
     &             EXP(-0.059e+0_fp * PDp)) / PDp

            ! Viscosity [Pa*s] of air as a function of temperature 
            VISC = 1.458e-6_fp * (Temp)**(1.5e+0_fp) / 
     &             ( Temp + 110.4e+0_fp )

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
         L    = LLCHEM
         DELZ = State_Met%BXHEIGHT(I,J,L)

         TC(I,J,L) = TC(I,J,L) / ( 1.e+0_fp + DTCHEM * VTS(L) / DELZ )

         DO L = LLCHEM-1, 1, -1
            DELZ  = State_Met%BXHEIGHT(I,J,L)
            DELZ1 = State_Met%BXHEIGHT(I,J,L+1)
            TC(I,J,L) = 1.e+0_fp / ( 1.e+0_fp + DTCHEM * VTS(L) / DELZ )
     &                * ( TC(I,J,L) + DTCHEM * VTS(L+1) / DELZ1
     &                *  TC(I,J,L+1) )
         ENDDO
         
#if !defined( NO_BPCH )
         !==============================================================
         ! ND44 diagnostic: sea salt loss [molec/cm2/s]
         !==============================================================
         IF ( ND44 > 0 ) THEN

            ! Initialize
            TOT1 = 0e+0_fp
            TOT2 = 0e+0_fp
            
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
#endif
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
      USE CHEMGRID_MOD,       ONLY : ITS_IN_THE_NOCHEMGRID
      USE CMN_GCTM_MOD         
      USE CMN_DIAG_MOD         
      USE CMN_SIZE_MOD         
      USE DIAG_MOD,           ONLY : AD05
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE GIGC_State_Chm_Mod, ONLY : ChmState
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE TIME_MOD,           ONLY : GET_TS_CHEM
      USE TRACERID_MOD,       ONLY : IDTDMS
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
!  24 Jul 2014 - R. Yantosca - Now compute BOXVL internally
!  06 Nov 2014 - R. Yantosca - Now use State_Met%AIRDEN(I,J,L)
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
      REAL(fp), PARAMETER :: FX = 1.0e+0_fp
      REAL(fp), PARAMETER :: A  = 0.75e+0_fp
      REAL(fp), PARAMETER :: B  = 0.25e+0_fp

      ! From D4: only 0.8 efficiency, also some goes to DMSO and lost.  
      ! So we assume 0.75 efficiency for DMS addtion channel to form     
      ! products.                                                        
      REAL(fp), PARAMETER :: EFF = 1e+0_fp
!
! !LOCAL VARIABLES:
!
      ! Scalars
      INTEGER           :: I,    J,    L
      REAL(fp)            :: TK,   O2,   RK1,    RK2,    RK3,   F  
      REAL(fp)            :: DMS,  DMS0, DMS_OH, DTCHEM, XOH,   XN3 
      REAL(fp)            :: XX,   OH,   OH0,    XNO3,   XNO30, LOH
      REAL(fp)            :: LNO3, BOXVL
                        
      ! For fields from Input_Opt
      LOGICAL           :: IS_FULLCHEM
      REAL(fp)            :: XNUMOL(Input_Opt%N_TRACERS)

      ! Pointers
      ! We need to define local arrays to hold corresponding values 
      ! from the Chemistry State (State_Chm) object. (mpayer, 12/6/12)
      REAL(fp), POINTER :: STT(:,:,:,:)
      
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
      DTCHEM      = GET_TS_CHEM() * 60e+0_fp

      ! Factor to convert AIRDEN from kgair/m3 to molecules/cm3:
      f           = 1000.e+0_fp / AIRMW * 6.022e+23_fp * 1.e-6_fp
      
      !=================================================================
      ! Do the chemistry over all chemically-active grid boxes!
      !=================================================================
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, TK, O2, DMS0,OH, XNO3, RK1, RK2, BOXVL )
!$OMP+PRIVATE( RK3, DMS_OH, DMS, OH0, XNO30, XOH, XN3, XX, LOH, LNO3  )
!$OMP+SCHEDULE( DYNAMIC )
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Skip non-chemistry boxes
         IF ( ITS_IN_THE_NOCHEMGRID( I, J, L, State_Met ) ) CYCLE

         ! Temperature [K]
         TK     = State_Met%T(I,J,L)

         ! Get O2 [molec/cm3], DMS [v/v], OH [molec/cm3], NO3 [molec/cm3]
         O2     = State_Met%AIRDEN(I,J,L) * f * 0.21e+0_fp
         DMS0   = STT(I,J,L,IDTDMS)
         OH     = GET_OH(  I, J, L, Input_Opt, State_Met )
         XNO3   = GET_NO3( I, J, L, Input_Opt, State_Met )

         !==============================================================
         ! (1) DMS + OH:  RK1 - addition channel  
         !                RK2 - abstraction channel   
         !==============================================================
         RK1 = 0.e+0_fp
         RK2 = 0.e+0_fp
         RK3 = 0.e+0_fp

         IF ( OH > 0.e+0_fp ) THEN
            RK1 = ( 1.7e-42_fp * EXP( 7810.e+0_fp / TK ) * O2 ) /
     &            ( 1.e+0_fp + 5.5e-31_fp * EXP( 7460.e+0_fp / TK ) 
     &            * O2 ) * OH

            ! Update reaction rate to match JPL06 and full chem
            ! (jaf, bmy, 10/15/09)
            RK2 = 1.1e-11_fp * EXP( -240.e+0_fp / TK ) * OH 
         ENDIF
            
         !==============================================================
         ! (2) DMS + NO3 (only happens at night):  
         !==============================================================
         IF ( State_Met%SUNCOS(I,J) <= 0e+0_fp ) THEN
            RK3 = 1.9e-13_fp * EXP( 500.e+0_fp / TK ) * XNO3
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
         IF ( DMS < SMALLNUM ) DMS = 0e+0_fp

         ! Archive initial OH and NO3 for diagnostics
         OH0    = OH
         XNO30  = XNO3

         IF ( IS_FULLCHEM ) THEN
         
            ! Update OH after rxn w/ DMS (coupled runs only)
            OH    = OH0 - ( ( DMS0 - DMS_OH ) *
     &              State_Met%AIRDEN(I,J,L) * f )
            IF ( OH < SMALLNUM ) OH = 0e+0_fp

            ! Update NO3 after rxn w/ DMS (coupled runs only)
            XNO3  = XNO30 - ( ( DMS_OH - DMS ) *
     &              State_Met%AIRDEN(I,J,L) * f )
            IF ( XNO3 < SMALLNUM ) XNO3 = 0e+0_fp

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
         IF ( ( RK1 + RK2 ) == 0.e+0_fp ) THEN
            PMSA_DMS(I,J,L) = 0.e+0_fp
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
        
            ! Grid box volume [cm3]
            BOXVL = State_Met%AIRVOL(I,J,L) * 1e+6_fp

            ! Convert L(OH) and L(NO3) from [molec/cm3] to [kg/timestep]
            LOH  = ( OH0   - OH  ) * BOXVL / XNUMOL_OH
            LNO3 = ( XNO30 - XNO3) * BOXVL / XNUMOL_NO3 

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
      USE CHEMGRID_MOD,       ONLY : ITS_IN_THE_NOCHEMGRID
      USE CMN_SIZE_MOD
      USE CMN_DIAG_MOD
      USE CMN_GCTM_MOD
#if !defined( NO_BPCH )
      USE DIAG_MOD,           ONLY : AD44 
#endif
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE GIGC_State_Chm_Mod, ONLY : ChmState
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE PBL_MIX_MOD,        ONLY : GET_FRAC_UNDER_PBLTOP
      USE TIME_MOD,           ONLY : GET_MONTH
      USE TIME_MOD,           ONLY : GET_TS_CHEM
      USE TIME_MOD,           ONLY : ITS_A_NEW_MONTH
      USE TRACERID_MOD,       ONLY : IDTH2O2
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
!  18 Sep 2014 - M. Sulprizio- Now get J(H2O2) and PH2O2m from HEMCO
!  06 Nov 2014 - R. Yantosca - Now use State_Met%AIRDEN(I,J,L)
!  12 Jun 2015 - R. Yantosca - Now remove orphaned ND44 variables
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
! 
      REAL(fp),  PARAMETER :: A = 2.9e-12_fp
!
! !LOCAL VARIABLES:
!
      ! Scalars
      LOGICAL            :: FIRST     = .TRUE.
      INTEGER, SAVE      :: LASTMONTH = -99
      INTEGER            :: I, J, L, JLOOP, N_TRACERS
      REAL(fp)           :: DT,    Koh,  DH2O2, M,    F ,   XTAU   
      REAL(fp)           :: H2O20, H2O2, ALPHA, FREQ, PHOTJ
      CHARACTER(LEN=255) :: FILENAME

      ! Arrays
      REAL*4             :: ARRAY(IIPAR,JJPAR,LLCHEM)
      REAL(fp)           :: TCVV  (Input_Opt%N_TRACERS)
      REAL(fp)           :: XNUMOL(Input_Opt%N_TRACERS)

      ! Pointers
      ! We need to define local arrays to hold corresponding values 
      ! from the Chemistry State (State_Chm) object. (mpayer, 12/6/12)
      REAL(fp), POINTER :: STT(:,:,:,:)

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

      ! Initialize GEOS-Chem tracer array [kg] from Chemistry State object
      ! (mpayer, 12/6/12)
      STT    => State_Chm%Tracers

      ! Chemistry timestep [s]
      DT        = GET_TS_CHEM() * 60e+0_fp

      ! Factor to convert AIRDEN from kgair/m3 to molecules/cm3:
      F         = 1000.e+0_fp / AIRMW * 6.022d23 * 1.e-6_fp
      
      !=================================================================
      ! Loop over tropopsheric grid boxes and do chemistry
      !=================================================================
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, M, H2O20, KOH, FREQ, ALPHA, DH2O2, H2O2 )
!$OMP+PRIVATE( JLOOP, PHOTJ )
!$OMP+SCHEDULE( DYNAMIC )
      DO L  = 1, LLPAR
      DO J  = 1, JJPAR
      DO I  = 1, IIPAR

         ! Initialize for safety's sake 
         FREQ = 0e+0_fp

         ! Skip non-chemistry boxes
         IF ( ITS_IN_THE_NOCHEMGRID( I, J, L, State_Met ) ) CYCLE

         ! Density of air [molec/cm3]
         M     = State_Met%AIRDEN(I,J,L) * f  

         ! Initial H2O2 [v/v]
         H2O20 = STT(I,J,L,IDTH2O2)

         ! Loss frequenty due to OH oxidation [s-1]
         KOH   = A * EXP( -160.e+0_fp / State_Met%T(I,J,L) ) * 
     &           GET_OH( I, J, L, Input_Opt, State_Met )

         ! Now do all dry deposition in mixing_mod.F90 (ckeller, 3/5/15)
         FREQ = 0.e+0_fp

         ! 1-D grid box index for SUNCOS
         JLOOP = ( (J-1) * IIPAR ) + I

         ! Impose a diurnal variation of jH2O2 by multiplying COS of 
         ! solar zenith angle normalized by maximum solar zenith angle 
         ! because the archived JH2O2 is for local noon time
         IF ( COSZM(I,J) > 0.e+0_fp ) THEN
            PHOTJ = JH2O2(I,J,L) * State_Met%SUNCOS(I,J) / COSZM(I,J)
            PHOTJ = MAX( PHOTJ, 0e+0_fp )
         ELSE
            PHOTJ = 0e+0_fp
         ENDIF

         ! Compute loss fraction from OH, photolysis, drydep [unitless].  
         ALPHA = 1.e+0_fp + ( KOH + PHOTJ + FREQ ) * DT 

         ! Delta H2O2 [v/v]
         ! PH2O2m is in kg/m3 (from HEMCO), convert to molec/cm3/s (mps,9/18/14)
         DH2O2 = ( PH2O2m(I,J,L) / TS_EMIS * XNUMOL_H2O2 / CM3PERM3 )
     &           * DT / ( ALPHA * M )
         
         ! Final H2O2 [v/v]
         H2O2  = ( H2O20 / ALPHA + DH2O2 )
         IF ( H2O2 < SMALLNUM ) H2O2 = 0e+0_fp

         ! Store final H2O2 in STT
         STT(I,J,L,IDTH2O2) = H2O2

      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

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
      USE CHEMGRID_MOD,       ONLY : ITS_IN_THE_NOCHEMGRID
      USE CMN_GCTM_MOD
      USE CMN_DIAG_MOD
      USE CMN_SIZE_MOD
      USE DAO_MOD,            ONLY : IS_WATER
      USE DIAG_MOD,           ONLY : AD05
      USE ERROR_MOD,          ONLY : IS_SAFE_EXP
      USE ERROR_MOD,          ONLY : SAFE_DIV
      USE GIGC_ErrCode_Mod
      USE ERROR_MOD,          ONLY : ERROR_STOP
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE GIGC_State_Chm_Mod, ONLY : ChmState
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE PRESSURE_MOD,       ONLY : GET_PCENTER
      USE TIME_MOD,           ONLY : GET_TS_CHEM, GET_MONTH
      USE TIME_MOD,           ONLY : ITS_A_NEW_MONTH
      USE TRACERID_MOD,       ONLY : IDTH2O2
      USE TRACERID_MOD,       ONLY : IDTSO2
      USE TRACERID_MOD,       ONLY : IDTSO4,  IDTNH3,  IDTNH4,  IDTHNO3
      USE TRACERID_MOD,       ONLY : IDTNIT,  IDTDST1, IDTDST2, IDTDST3
      USE TRACERID_MOD,       ONLY : IDTDST4, IDTSALA, IDTSALC
      USE WETSCAV_MOD,        ONLY : H2O2s
      USE WETSCAV_MOD,        ONLY : SO2s
      USE HCOI_GC_MAIN_MOD,   ONLY : GetHcoDiagn
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
      INTEGER,        INTENT(INOUT) :: RC          ! Success or failure?
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
!  06 Sep 2013 - M. Sulprizio- Bug fix: Prevent divide-by-zero if LWC=0. Only
!                              do aqueous SO2 chemistry when LWC>0.
!  26 Sep 2013 - R. Yantosca - Renamed GEOS_57 Cpp switch to GEOS_FP
!  28 Jan 2014 - R. Yantosca - Bug fix for TOMAS. Set ALKdst=0 since TOMAS
!                              carries its own dust tracers instead of DST1-4.
!  25 Jun 2014 - R. Yantosca - Now pass Input_Opt to GET_ALK
!  18 Sep 2014 - M. Sulprizio- Now get HNO3 for offline aerosol sim from HEMCO
!  06 Nov 2014 - R. Yantosca - Now use State_Met%AIRDEN(I,J,L)
!  06 Nov 2014 - R. Yantosca - Now use State_Met%CLDF(I,J,L)
!  12 Jan 2015 - C. Keller   - Now allow NDENS_SALA and NDENS_SALC to be empty.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
!      REAL(fp),  PARAMETER    :: HPLUS  = 3.16227766016837953d-5  !pH = 4.5
      REAL(fp),  PARAMETER  :: HPLUS_45  = 3.16227766016837953e-5_fp  !pH = 4.5
      REAL(fp),  PARAMETER  :: HPLUS_50  = 1.0e-5_fp  !pH = 5.0
      REAL(fp),  PARAMETER  :: MINDAT = 1.e-20_fp
!
! !LOCAL VARIABLES:
!
      ! Scalars
      LOGICAL               :: IS_OFFLINE
      LOGICAL               :: IS_FULLCHEM
      INTEGER               :: I,      J,       L       !,      I1,   I2
      INTEGER               :: II,     NSTEP,   N_TRACERS
      INTEGER               :: BULK,   SIZE_RES 
      REAL(fp)              :: K0,     Ki,      KK,     M,    L1
      REAL(fp)              :: L2,     L3,      Ld,     F,    Fc
      REAL(fp)              :: RK,     RKT,     DTCHEM, DT_T, TK
      REAL(fp)              :: F1,     RK1,     RK2,    RK3,  SO20
      REAL(fp)              :: SO2_cd, H2O20,   O3,     L2S,  L3S
      REAL(fp)              :: LWC,    KaqH2O2, KaqO3,  PATM
      REAL(fp)              :: ALK,    ALK1,    ALK2,    SO2_ss
      REAL(fp)              :: Kt1,    Kt2,     AREASS1, AREASS2
      REAL(fp)              :: PSO4E,  PSO4F,   Kt1N,    Kt2N
      REAL(fp)              :: XX
      REAL(fp)              :: HPLUS,  SO4nss, TNH3,   TNO3,  GNO3, ANIT
      REAL(fp)              :: LSTOT,  ALKdst, ALKss,  ALKds, NH3
      REAL(fp)              :: SSCvv,  aSO4,   SO2_sr, SR,    TANIT

      ! Arrays
      REAL(fp)              :: TCVV  (Input_Opt%N_TRACERS)
      REAL(fp)              :: XNUMOL(Input_Opt%N_TRACERS)

      ! Pointers
      ! We need to define local arrays to hold corresponding values 
      ! from the Chemistry State (State_Chm) object. (mpayer, 12/6/12)
      REAL(fp), POINTER :: STT(:,:,:,:)

      ! For HEMCO update
      LOGICAL, SAVE   :: FIRST = .TRUE.

      CHARACTER(LEN=255), PARAMETER :: LOC = 'CHEM_SO2 (sulfate_mod.F)'

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

      ! Initialize GEOS-Chem tracer array [kg] from Chemistry State object
      ! (mpayer, 12/6/12)
      STT => State_Chm%Tracers

      ! DTCHEM is the chemistry timestep in seconds
      DTCHEM   = GET_TS_CHEM() * 60e+0_fp

      ! Factor to convert AIRDEN from [kg air/m3] to [molec air/cm3]
      F        = 1000.e+0_fp / AIRMW * 6.022e+23_fp * 1.e-6_fp

      ! On first call, get pointers to HEMCO diagnostics arrays.
      ! These are the sea salt aerosol number densities for the fine
      ! and coarse mode, respectively. Values are in # / surface grid
      ! box. These values are needed in the GET_ALK call below.
      ! If the diagnostics are not being found, e.g. because the
      ! sea salt emissions extension is turned off (or LEMIS is 
      ! disabled), the passed pointers NDENS_SALA and NDENS_SALC
      ! will stay nullified. Values of zero will be used in this
      ! case! (ckeller, 01/12/2015)
      IF ( FIRST ) THEN

         ! Sea salt density, fine mode
         CALL GetHcoDiagn( am_I_Root, 'SEASALT_DENS_FINE', 
     &          StopIfNotFound=.FALSE., RC=RC, Ptr2D=NDENS_SALA )
         IF ( RC /= HCO_SUCCESS )
     &      CALL ERROR_STOP( 'Cannot get SEASALT_DENS_FINE', LOC )

         ! Sea salt density, coarse mode 
         CALL GetHcoDiagn( am_I_Root, 'SEASALT_DENS_COARSE', 
     &          StopIfNotFound=.FALSE., RC=RC, Ptr2D=NDENS_SALC )
         IF ( RC /= HCO_SUCCESS )
     &      CALL ERROR_STOP( 'Cannot get SEASALT_DENS_COARSE', LOC )
         
         ! Adjust first flag
         FIRST = .FALSE.
      ENDIF 
      
      ! Loop over chemistry grid boxes
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, SO20, H2O20, O3, PATM, TK, K0, M, KK, F1, RK1  )
!$OMP+PRIVATE( RK2, RK, RKT, SO2_cd, L1, Ld, L2, L2S, L3, L3S, FC, LWC )
!$OMP+PRIVATE( KaqH2O2, KaqO3, ALK, ALK1, ALK2                         )
!$OMP+PRIVATE( Kt1, Kt2, AREASS1, AREASS2, SO2_ss, Kt1N, Kt2N          )
!$OMP+PRIVATE( PSO4E, PSO4F, XX                                        )
!$OMP+PRIVATE( HPLUS, SO4nss, TNH3, TNO3,  GNO3, ANIT,   LSTOT, ALKdst )
!$OMP+PRIVATE( ALKds, ALKss,  NH3,  SSCvv, aSO4, SO2_sr, SR,    TANIT  )
!$OMP+PRIVATE( BULK,  SIZE_RES,     RC                                 )
!$OMP+SCHEDULE( DYNAMIC )
      DO L = 1, LLPAR 
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Initialize for safety's sake 
         Ld       = 0e+0_fp

         ! Skip non-chemistry boxes
         IF ( ITS_IN_THE_NOCHEMGRID( I, J, L, State_Met ) ) CYCLE

         ! Initial SO2, H2O2 and O3 [v/v]
         SO20   = STT(I,J,L,IDTSO2)         
         H2O20  = STT(I,J,L,IDTH2O2)
         O3     = GET_O3( I, J, L, Input_Opt, State_Met )

         ! PATM  : Atmospheric pressure in atm
         PATM   = GET_PCENTER( I, J, L ) / 1013.25e+0_fp

         ! TK : Temperature [K]
         TK     = State_Met%T(I,J,L)

         IF ( IS_OFFLINE ) THEN

            ! Gas phase SO4 production is done here in offline run only 
            ! Updated to match JPL 2006 + full chem (jaf, 10/14/09)
            K0  = 3.3e-31_fp * ( 300.e+0_fp / TK )**4.3e+0_fp
            Ki  = 1.6e-12_fp
            M   = State_Met%AIRDEN(I,J,L) * F
            KK  = K0 * M / Ki
            F1  = ( 1.e+0_fp + ( LOG10( KK ) )**2 )**( -1 )
            RK1 = ( K0 * M / ( 1.e+0_fp + KK ) ) * 0.6e+0_fp**F1 * 
     &            GET_OH( I, J, L, Input_Opt, State_Met )

         ELSE 

            ! For online runs, SMVGEAR deals w/ this computation,
            ! so we can simply set RK1 = 0 (rjp, bmy, 3/23/03)
            K0  = 0.e+0_fp
            M   = 0.e+0_fp
            KK  = 0.e+0_fp
            F1  = 0.e+0_fp
            RK1 = 0.e+0_fp

         ENDIF

         ! Now do all dry deposition in mixing_mod.F90 (ckeller, 3/5/15)
         RK2  = 0.e+0_fp

         ! RK: total reaction rate [1/s]
         RK     = ( RK1 + RK2 )
       
         ! RKT: RK * DTCHEM [unitless] (bmy, 6/1/00)
         RKT    =  RK * DTCHEM

         !==============================================================
         ! Update SO2 conc. after gas phase chemistry and deposition
         !==============================================================
         IF ( RK > 0.e+0_fp ) THEN
            SO2_cd = ( SO20  * EXP( -RKT ) ) +
     &          ( PSO2_DMS(I,J,L) * ( 1.e+0_fp - EXP( -RKT ) ) / RKT )

            L1     = ( SO20 - SO2_cd + PSO2_DMS(I,J,L) ) * RK1/RK
             
            Ld     = ( SO20 - SO2_cd + PSO2_DMS(I,J,L) ) * RK2/RK
            
         ELSE
            SO2_cd = SO20
            L1     = 0.e+0_fp
         ENDIF

         !==============================================================
         ! Update SO2 conc. after seasalt chemistry (bec, 12/7/04)
         !==============================================================

         ! Get alkalinity of accum (ALK1) and coarse (ALK2) [kg]
         CALL GET_ALK( am_I_Root, I, J, L, ALK1, ALK2, Kt1, Kt2,  
     &                 Kt1N, Kt2N, Input_Opt, State_Met, RC )

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
            PSO4E        = 0.e+0_fp
            PSO4F        = 0.e+0_fp
            PNITS(I,J,L) = 0.e+0_fp

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
         FC      = State_Met%CLDF(I,J,L)

         ! Get liquid water content [m3 H2O/m3 air] within cloud from met flds
         ! Units: [kg H2O/kg air] * [kg air/m3 air] * [m3 H2O/1e3 kg H2O]
         LWC     = State_Met%QL(I,J,L) * State_Met%AIRDEN(I,J,L) *
     &             1e-3_fp

         ! LWC is a grid-box averaged quantity. To improve the representation 
         ! of sulfate chemistry, we divide LWC by the cloud fraction and 
         ! compute sulfate chemistry based on the LWC within the cloud.  We 
         ! get the appropriate grid-box averaged mass of SO2 and sulfate by 
         ! multiplying these quantities by FC AFTER computing the aqueous 
         ! sulfur chemistry within the cloud. (lzh, jaf, bmy, 5/27/11)
         LWC     = SAFE_DIV( LWC, FC, 0e+0_fp )

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
         KaqH2O2 = 0.e+0_fp
         KaqO3   = 0.e+0_fp
         L2      = 0.e+0_fp
         L3      = 0.e+0_fp
         L2S     = 0.e+0_fp
         L3S     = 0.e+0_fp
         
         ! If (1) there is cloud, (2) there is SO2 present, and 
         ! (3) the T > -15 C, then compute aqueous SO2 chemistry
         ! Prevent divide-by-zero if LWC=0 (mpayer, 9/6/13)
         IF ( ( FC     > 0.e+0_fp   )  .AND. 
     &        ( SO2_ss > MINDAT )  .AND. 
     &        ( TK     > 258.0  )  .AND.
     &        ( LWC    > 0.e+0_fp   ) ) THEN

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
	    SO4nss  =  STT(I,J,L,IDTSO4) * State_Met%AIRDEN(I,J,L) *
     &                 0.7e+0_fp / ( 28.97e+0_fp * LWC )

	    ! Get total ammonia (NH3 + NH4+) concentration [v/v]
	    ! Use a cloud scavenging ratio of 0.7 for NH4+
            TNH3     = ( STT(I,J,L,IDTNH4) * 0.7e+0_fp ) + 
     &                   STT(I,J,L,IDTNH3)

	    ! Get total nitrate (HNO3 + NIT) concentrations [v/v] 
     	    ! Use a cloud scavenging ratio of 0.7 for NIT
	    IF ( IS_FULLCHEM ) THEN
               TNO3 = STT(I,J,L,IDTHNO3) + 
     &              ( STT(I,J,L,IDTNIT) * 0.7e+0_fp )
	       GNO3 = STT(I,J,L,IDTHNO3) !For Fahey & Pandis decision algorithm
            ELSE IF ( IS_OFFLINE ) THEN
               TANIT = STT(I,J,L,IDTNIT) !aerosol nitrate [v/v]
	       GNO3  = HNO3(I,J,L) - TANIT ! gas-phase nitric acid [v/v]
   	       ANIT  = TANIT * 0.7e+0_fp ! aerosol nitrate in the cloud drops [v/v]
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
            IF ( IS_SAFE_EXP( XX ) .and. ABS( XX ) > 0e+0_fp ) THEN

               ! Aqueous phase SO2 loss rate w/ H2O2 [v/v/timestep]
               L2  = EXP( XX )

               ! Loss by H2O2
               L2S = SO2_ss * H2O20 * ( L2 - 1.e+0_fp ) / 
     &               ( (SO2_ss * L2) - H2O20 )  
            ELSE

               ! NOTE from Jintai Lin (4/28/10):
               ! However, in the case of a negative XX, L2S should be 
               ! approximated as SO2_ss, instead of H2O20. In other words, 
               ! L2S = SO2_ss * H2O20 * ( L2 - 1.D0 ) / ( (SO2_ss*L2) - H2O20 )
               ! reaches different limits when XX reaches positive infinity 
               ! and negative infinity.
               IF ( XX > 0.e+0_fp ) THEN 
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
            IF ( IS_SAFE_EXP( XX ) .and. ABS( XX ) > 0e+0_fp ) THEN

               ! Aqueous phase SO2 loss rate w/ O3 [v/v/timestep]
               L3  = EXP( XX )

               ! Loss by O3
               L3S = SO2_ss * O3 * (L3 - 1.e+0_fp)/((SO2_ss * L3) - O3)

            ELSE
 
               ! Follow the same logic for L3S as described in
               ! Jintai Lin's note above (bmy, 4/28/10)
               IF ( XX > 0.e+0_fp ) THEN 
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
            ALKdst = 0e+0_fp
#else
            ! For other simulations, Sum up the contributions from
            ! DST1 thru DST4 tracers into ALKdst. (bmy, 1/28/14)
            ALKdst = ( STT(I,J,L,IDTDST1) + STT(I,J,L,IDTDST2) +
     &                 STT(I,J,L,IDTDST3) + STT(I,J,L,IDTDST4) ) *
     &                 1.e+9_fp * State_Met%AD(I,J,L) / TCVV(IDTDST1) /
     &                 State_Met%AIRVOL(I,J,L)
#endif

            ALKss  = ( STT(I,J,L,IDTSALA  ) + STT(I,J,L,IDTSALC) ) *
     &                 1.e+9_fp * State_Met%AD(I,J,L) / TCVV(IDTSALA) /
     &                 State_Met%AIRVOL(I,J,L)

            ALKds = ALKdst + ALKss

            ! Get NH3 concentrations (v/v)
            NH3 = STT(I,J,L,IDTNH3) 

            ! Initialize
            BULK = 0
	    SIZE_RES = 0

            ! Fahey and Seinfeld decision algorithm
	    IF ( H2O20 > SO2_ss + 1e-9_fp ) THEN
               BULK = 1
            ELSEIF( LWC < 0.1e+6_fp ) THEN !10^6 coversion from m3/m3 --> g/m3
               SIZE_RES = 1
            ELSEIF( gno3 > NH3 ) THEN
               IF ( SO2_ss >= 5.e-9_fp          .and. 
     &              H2O20  >= SO20   ) 
     &              BULK    = 1
               IF ( LWC    >= 0.3e+6_fp          .and. 
     &              SO2_ss >= 3.e-9_fp          .and. 
     &              H2O20  >= SO2_ss ) 
     &              BULK    = 1
               IF ( ALKds  >= 5.e+0_fp           .and. 
     &              LWC    >= 0.5e+6_fp          .and. 
     &              H2O20  >= SO2_ss ) 
     &              BULK    = 1
               IF ( LWC    >= 0.1e+6_fp          .and. 
     &              gno3   <= (NH3 + 2.e-9_fp) ) 
     &              BULK    = 1
            ELSEIF( LWC    >= 0.1e+6_fp ) THEN
               IF ( NH3    <= 1.e-9_fp          .and. 
     &              ALKds  >= 5.e+0_fp   ) 
     &              BULK    = 1
               IF ( NH3    <= 5.e-9_fp          .and. 
     &              ALKds  >= 10.e+0_fp  ) 
     &              BULK    = 1
               IF ( gno3   <= 1.e-9_fp          .and. 
     &              NH3    >= (gno3 + 2.e-9_fp) .and. 
     &              SO2_ss <= 7.e-9_fp )
     &              BULK    = 1
               IF ( gno3   <= 1.e-9_fp          .and. 
     &              NH3    >= (gno3 + 2.e-9_fp) .and. 
     &              ALKds  >= 2.e+0_fp )  BULK = 1
               IF ( gno3   <= 3.e-9_fp          .and. 
     &              NH3    >= (gno3 + 4.e-9_fp) )
     &              BULK    = 1
               IF ( gno3   <= 7.e-9_fp          .and. 
     &              NH3    >= (gno3 + 3.e-9_fp) .and. 
     &              SO2_ss <= 5.e-9_fp )
     &              BULK    = 1
               IF ( gno3   <= 7.e-9_fp          .and. 
     &              NH3    >= (gno3 + 3.e-9_fp) .and. 
     &		    ALKds  >= 4.e+0_fp           .and. 
     &              SO2_ss <= 9.e-9_fp  ) 
     &              BULK    = 1
               IF ( ALKds  >= 3.e+0_fp           .and. 
     &              NH3    <= 3.e-9_fp          .and. 
     &              SO2_ss <= 4.e-9_fp ) 
     &              BULK    = 1
               IF ( ALKds  >= 5.e+0_fp           .and. 
     &              SO2_ss <= 5.e-9_fp          .and. 
     &              NH3    <= 7.e-9_fp )  
     &              BULK    = 1
               IF ( NH3    >= (gno3 + 2.e-9_fp) .and. 
     &              SO2_ss <= 5.e-9_fp )  
     &              BULK    = 1
               IF ( NH3    >= (gno3 + 4.e-9_fp) .and. 
     &              SO2_ss <= 10.e-9_fp )
     &              BULK    = 1
               IF ( ALKds  >= 2.e+0_fp           .and. 
     &              NH3    <= 10.e-9_fp         .and. 
     &              H2O20  >= SO2_ss ) 
     &              BULK    = 1
               IF ( NH3    <= 1.e-9_fp          .and. 
     &              SO2_ss >= 3.e-9_fp          .and. 
     &              H2O20  >= SO2_ss ) 
     &              BULK    = 1
            ELSEIF ( LWC   >= 0.3d6 ) THEN
               IF ( NH3    >= (gno3 + 5.e-9_fp) .and. 
     &              SO2_ss <= 10.e-9_fp )
     &              BULK    = 1
               IF ( gno3   <= 1.e-9_fp          .and. 
     &              NH3    >= (gno3 + 2.e-9_fp) )
     &              BULK    = 1
               IF ( gno3   <= 7.e-9_fp          .and. 
     &              NH3    >= (gno3 + 3.e-9_fp) )
     &              BULK    = 1
               IF ( ALKds  >= 3.e+0_fp           .and. 
     &              NH3 <= 10e-9_fp             .and. 
     &              SO2_ss <= 5e-9_fp ) 
     &              BULK    = 1
               IF ( ALKds  >= 5.e+0_fp           .and. 
     &              NH3    <= 10.e-9_fp         .and. 
     &              SO2_ss <= 5.e-9_fp ) 
     &              BULK    = 1
               IF ( SO2_ss >= 1.5e-9_fp         .and. 
     &              H2O20  >= SO2_ss ) 
     &              BULK    = 1
               IF ( NH3    <= 12.e-9_fp         .and. 
     &              ALKds  >=10.e+0_fp ) 
     &              BULK    = 1
               IF ( NH3    <= 1.e-9_fp          .and. 
     &              ALKds  >= 4.e+0_fp           .and. 
     &              SO2_ss <= 10.e-9_fp ) 
     &              BULK    = 1
               IF ( NH3    <= 5.e-9_fp          .and. 
     &              ALKds  >= 6.e+0_fp           .and. 
     &              SO2_ss <= 10.e-9_fp ) 
     &              BULK    = 1
               IF ( NH3    <= 7.e-9_fp          .and. 
     &              ALKds   >-8.e+0_fp           .and. 
     &              SO2_ss <= 10.e-9_fp ) 
     &              BULK    = 1
            ELSEIF ( LWC   >= 0.5e+6_fp ) THEN
               IF ( H2O20  >= (0.9e+0_fp * SO2_ss) ) 
     &              BULK    = 1
               IF ( NH3    <= 1.e-9_fp          .and. 
     &              ALKds  >= 5.e+0_fp           .and. 
     &              SO2_ss <= 10.e-9_fp ) 
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
     &                             Input_Opt, State_Met )

	    ELSE
		SR = 0.e+0_fp
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
            H2O2s(I,J,L) = MAX( H2O20,  1.0e-32_fp )
            SO2s(I,J,L ) = MAX( SO2_ss, 1.0e-32_fp )
            L2S          = 0.e+0_fp
            L3S          = 0.e+0_fp
	    SR		 = 0.e+0_fp

         ENDIF

         ! Store updated SO2, H2O2 back to the tracer arrays 
         STT(I,J,L,IDTSO2)  = SO2s( I,J,L)
         STT(I,J,L,IDTH2O2) = H2O2s(I,J,L)

         ! SO2 chemical loss rate  = SO4 production rate [v/v/timestep]
!         PSO4_SO2(I,J,L) = L1 + L2S + L3S + PSO4E
         PSO4_SO2(I,J,L) = L1 + L2S + L3S + PSO4E + SR
         PSO4_ss (I,J,L) = PSO4F

#if   defined( TOMAS )
         PSO4_SO2AQ(I,J,L) = L2S + L3S + SR ! For TOMAS microphysics
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
      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

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
      REAL(fp),         INTENT(IN)   :: SO2_cd     ! SO2 mixing ratio [v/v] after 
                                                 !  gas phase chemistry and 
                                                 !  dry deposition
      REAL(fp),         INTENT(IN)   :: Kt1, Kt2   ! Rate constant [s-1] for 
                                                 !  sulfate formation on sea
                                                 !  salt aerosols from GET_ALK 
                                                 !  (1=fine; 2=coarse)
      REAL(fp),         INTENT(IN)   :: Kt1N, Kt2N 
      REAL(fp),         INTENT(IN)   :: ALK1, ALK2 ! Alkalinity [kg] from 
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
      REAL(fp),         INTENT(OUT)   :: SO2_ss    ! SO2 mixing ratio [v/v] 
                                                 !  after sea salt chemistry
      REAL(fp),         INTENT(OUT)   :: PSO4E     ! SO4E (sulfate produced by 
                                                 !  S(IV)+O3 on fine seasalt) 
                                                 !  mixing ratio [v/v]
      REAL(fp),         INTENT(OUT)   :: PSO4F     ! SO4F (sulfate produced by 
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
!  18 Sep 2014 - M. Sulprizio- Now get HNO3 for offline aerosol sim from HEMCO
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
      REAL(fp),  PARAMETER :: MINDAT    = 1.0e-20_fp 
      REAL(fp),  PARAMETER :: TCVV_HNO3 = 28.97e+0_fp / 63.0e+0_fp 
!
! !LOCAL VARIABLES:
!
      ! Scalars
      LOGICAL            :: IT_IS_A_FULLCHEM_SIM
      INTEGER            :: JLOOP,       N_TRACERS
      REAL(fp)             :: SO2_chem,    DTCHEM
      REAL(fp)             :: O3_cspec,    O3_lost
      REAL(fp)             :: EQ_1_C,      EQ_2_C
      REAL(fp)             :: SO4E,        SO2_new,    SO4F
      REAL(fp)             :: SO2_eq,      N_FLUX_A,   N_FLUX_C
      REAL(fp)             :: END_ALK,     L5A,        L5C
      REAL(fp)             :: EQ1,         EQ2,        TITR_SO2
      REAL(fp)             :: TITR_HNO3,   NIT_vv,     NITs_vv
      REAL(fp)             :: NIT0,        NITS0
      REAL(fp)             :: F_SO2,       FALK_A_SO2, FALK_C_SO2
      REAL(fp)             :: EQ_BEG,      F_SO2_A,    F_SO2_C
      REAL(fp)             :: ALKA,        ALKC,       TOTAL_ACID_FLUX
      REAL(fp)             :: HNO3_EQ,     TOT_FLUX_A, TOT_FLUX_C
      REAL(fp)             :: FALK_A_HNO3, HNO3_vv
      REAL(fp)             :: FALK_C_HNO3, F_HNO3_A,   F_HNO3_C
      REAL(fp)             :: EQ_1_N,      EQ_2_N,     F_HNO3
      REAL(fp)             :: HNO3_SSA,    HNO3_SSC,   N_FLUX
      REAL(fp)             :: HNO3_EQ_C,   L6A,        L6C   
      REAL(fp)             :: C_FLUX_A,    C_FLUX_C,   C_FLUX
      REAL(fp)             :: HNO3_ss,     HNO3_kg

      ! Arrays
      REAL(fp)             :: TCVV(Input_Opt%N_TRACERS)

      ! Pointers
      REAL(fp), POINTER    :: STT(:,:,:,:)
      REAL(fp), POINTER    :: AD(:,:,:)
      REAL(fp), POINTER    :: AIRVOL(:,:,:)

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
      DTCHEM = GET_TS_CHEM() * 60e+0_fp

      ! Convert SO2 [v/v] to  [eq/gridbox]
      SO2_eq = ( 2.e+0_fp * SO2_cd  * AD(I,J,L) ) / 
     &         ( TCVV(IDTSO2) * 0.064e+0_fp )
      SO2_eq = MAX( SO2_eq, MINDAT )

      IF ( IT_IS_A_FULLCHEM_SIM ) THEN

 	 ! Convert HNO3 [v/v] to [equivalents]
         HNO3_vv = STT(I,J,L,IDTHNO3)
         HNO3_eq = HNO3_vv * AD(I,J,L) / ( 28.97e+0_fp / 63e+0_fp ) 
     &           / 63.e-3_fp

      ELSE
         
         ! HNO3 is in v/v (from HEMCO)
         HNO3_vv = HNO3(I,J,L)

         ! Convert HNO3 [v/v] to [equivalents]
         HNO3_eq = HNO3_vv * AD(I,J,L) / ( 28.97e+0_fp / 63e+0_fp ) 
     &           / 63.e-3_fp

      ENDIF

      !-----------
      ! SO2
      !-----------

      ! Available flux of SO2 to accum sea salt aerosols [v/v/timestep]
      L5A      = EXP( -Kt1 * DTCHEM )
      F_SO2_A  = SO2_cd * ( 1.e+0_fp - L5A )
      F_SO2_A  = MAX( F_SO2_A, 1.e-32_fp )

      ! Convert to [eq/timestep] 
      C_FLUX_A = 2.e+0_fp * F_SO2_A * AD(I,J,L) / TCVV(IDTSO2) 
     &         / 0.064e+0_fp

      ! Available flux of SO2 to coarse sea salt aerosols [v/v/timestep]
      L5C      = EXP( - Kt2 * DTCHEM )
      F_SO2_C  = SO2_cd * ( 1.e+0_fp - L5C )
      F_SO2_C  = MAX( F_SO2_C, 1.0e-32_fp )

      ! Convert to [eq/timestep] 
      C_FLUX_C = 2.e+0_fp * F_SO2_C * AD(I,J,L) / TCVV(IDTSO2) 
     &         / 0.064e+0_fp

      ! Total flux of SO2 [v/v/timestep]
      F_SO2    = F_SO2_A + F_SO2_C 

      ! Total flux of SO2 [eq/timestep]
      C_FLUX   = C_FLUX_A + C_FLUX_C 

      !-----------
      ! HNO3
      !-----------

      ! Available flux of HNO3 to accum sea salt aerosols [v/v/timestep]
      L6A = EXP( - Kt1N * DTCHEM )
      F_HNO3_A = HNO3_vv * ( 1.e+0_fp - L6A )
      F_HNO3_A = MAX( F_HNO3_A, 1.0e-32_fp )

      ! Convert to [eq/timestep] 
      N_FLUX_A = F_HNO3_A * AD(I,J,L) / ( 28.97e+0_fp / 63e+0_fp ) 
     &         / 0.063e+0_fp

      ! Available flux of HNO3 to coarse sea salt aerosols [v/v/timestep]
      L6C = EXP( - Kt2N * DTCHEM )
      F_HNO3_C = HNO3_vv * ( 1.e+0_fp - L6C )
      F_HNO3_C = MAX( F_HNO3_C, 1.0e-32_fp )

      ! convert to [eq/timestep] 
      N_FLUX_C = F_HNO3_C * AD(I,J,L) / ( 28.97e+0_fp / 63e+0_fp ) 
     &         / 0.063e+0_fp

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
      EQ1 = ALK1 * 0.07e+0_fp
      EQ2 = ALK2 * 0.07e+0_fp

      !----------------------------------------------------------------------
      ! NOTE: This was a sensitivity simulation, keep for future reference
      !       cf Alexander et al 2005 (bec, bmy, 4/13/05)
      !! Total available alkalinity [eq] doubled for Sievering run
      !EQ1 = ALK1 * 0.14e+0_fp
      !EQ2 = ALK2 * 0.14e+0_fp
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

	 FALK_A_SO2  = 1.0e+0_fp
	 FALK_A_HNO3 = 1.0e+0_fp

      ENDIF
      
      IF ( TOT_FLUX_C > EQ2 ) THEN

         ! Fraction of flkalinity available for each acid
	 FALK_C_SO2  = C_FLUX_C/TOT_FLUX_C
	 FALK_C_HNO3 = N_FLUX_C/TOT_FLUX_C
         FALK_C_SO2  = MAX( FALK_C_SO2, MINDAT )
	 FALK_C_HNO3 = MAX( FALK_C_HNO3, MINDAT )

      ELSE

	 FALK_C_SO2  = 1.0e+0_fp
	 FALK_C_HNO3 = 1.0e+0_fp

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
     &               AD(I,J,L) / 2.0e+0_fp
      SO2_ss       = MAX( SO2_ss, MINDAT )

      !SO4E produced converted from [eq/timestep] to [v/v/timestep] 
      PSO4E        = SO4E * 0.096 * TCVV(IDTSO4) /
     &               AD(I,J,L) / 2.0e+0_fp

      !SO4F produced converted from [eq/timestep] to [v/v/timestep] 
      PSO4F        = SO4F * 0.096 * TCVV(IDTSO4S) /
     &               AD(I,J,L) / 2.0e+0_fp

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
      REAL(fp),  INTENT(IN)    :: SO4nss ! Total nss sulfate mixing ratio [M]
      REAL(fp),  INTENT(IN)    :: TNO3   ! Total nitrate (gas+particulate) mixing
                                       ! ratio [v/v]
      REAL(fp),  INTENT(IN)    :: TNH3   ! NH3 mixing ratio [v/v]
      REAL(fp),  INTENT(IN)    :: SO2    ! SO2 mixing ratio [v/v] 
      REAL(fp),  INTENT(IN)    :: T      ! Temperature [K]
      REAL(fp),  INTENT(IN)    :: PRES   ! Pressure [atm]
      REAL(fp),  INTENT(IN)    :: LWC    ! Cloud liquid water content [m3/m3]
      REAL(fp),  INTENT(IN)    :: iHPLUS ! Initial [H+] [M]
!
! !OUTPUT PARAMETERS: 
!
      REAL(fp),  INTENT(OUT)   :: HPLUS  ! Calculated [H+] [M]
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
!  06 Jan 2015 - M. Yannetti - Set some variables to f8 that required it
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!

      ! Water dissociation constants
      REAL(fp),  PARAMETER     :: Kw   = 1.0e-14_fp
      REAL(fp),  PARAMETER     :: DhrKw = -6710.e+0_fp
      REAL(fp),  PARAMETER     :: MINVAL = 0.01
!
! !LOCAL VARIABLES:
!
      REAL(fp)                 :: D, Kw_T, ipH, newpH, nHPLUS
      REAL(fp)                 :: kCO2p, kCO2p2
      REAL(fp)                 :: kSO2p, kSO2p2
      REAL(fp)                 :: kHNO3p, kNH3p, kHClp
      INTEGER                :: count
      REAL(fp)                 :: E, F, G, H, Q, R
      REAL(fp)                 :: X, Y
      REAL(fp)                 :: CRUTES( 3 ) ! Coeff and roots of cubic equation
      REAL(f8)                 :: A, B, A1, A2, A0, P
      INTEGER                :: NR          ! Number of roots to cubic equation

      !=================================================================
      ! GET_HPLUS begins here!
      !=================================================================

      ! Initial pH guess 
      ipH = -log10(iHPLUS)

      ! Non-volatile aerosol concentration [M] 
      ! For now sulfate is the only non-volatile species
      D = (2.e+0_fp*SO4nss) 

      ! Temperature dependent water equilibrium constant
      Kw_T = Kw*exp(DhrKw*((1.e+0_fp/T)-(1.e+0_fp/298.e+0_fp)))

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

         nHPLUS = 10.e+0_fp**(-ipH)

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

	 A = (1.e+0_f8/3.e+0_f8)*((3.e+0_f8*Q)-(P*P))

  	 B = (1.e+0_f8/27.e+0_f8)*((2.0e+0_f8*P*P*P)
     &        -(9.0e+0_f8*P*Q)+(27.0e+0_f8*R))

 	 A2 = 0.
 	 A1 = A
 	 A0 = B

	 !write(6,*) 'calling CUBIC'
         CALL CUBIC ( A2, A1, A0, NR, CRUTES )
	 !write(6,*) 'after CUBIC'

         ! Code assumes the smallest positive root is in CRUTES(1)
         X = CRUTES( 1 )

	 ! Y = [H+]
 	 Y = X - (P/3.0e+0_fp)

  	 ! Set minimum [H+] = 1.d-14 (pH = 14)
    	 Y = MAX(Y,1.0e-14_fp)

  	 ! Set maximum [H+] = 1.d-1 (pH = 1)
    	 Y = MIN(Y,1.0e-1_fp)

  	 ! If solution does not converge after 5 iterations
  	 ! average last 2 pH calculations
  	 IF (count > 5) THEN
	    newpH = ((-log10(Y)) + (-log10(nHPLUS))) / 2.0e+0_fp

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

      HPLUS = 10.0e+0_fp**(-newpH)

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
      REAL(fp),  INTENT(IN) :: T, P, LWC, HPLUS
!
! !OUTPUT PARAMETERS: 
!
      REAL(fp)              :: KCO2p, KCO2p2
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
      REAL(fp),  PARAMETER  :: Kc1 = 4.3e-7
      REAL(fp),  PARAMETER  :: Kc2 =4.68e-11
      REAL(fp),  PARAMETER  :: DhrKc1 = -1000.
      REAL(fp),  PARAMETER  :: DhrKc2 = -1760.
      REAL(fp),  PARAMETER  :: Hco2 = 3.4e-2
      REAL(fp),  PARAMETER  :: Dhco2 = -4.85/0.04
      ! CO2 concentration [v/v]
      REAL(fp),  PARAMETER  :: CO2 = 380.0e-6_fp !380 ppmv
!
! !LOCAL VARIABLES:
!
      REAL(fp)              :: Hco2_T, Kc1_T, Kc2_T
      REAL(fp)              :: Hco2eff, xCO2, pCO2

      !=================================================================
      ! kCO2 begins here!
      !=================================================================

      !CO2 dissolution constants
      Hco2_T = Hco2*exp(Dhco2*((1.e+0_fp/T)-(1.e+0_fp/298.e+0_fp)))
      Kc1_T = Kc1*exp(DhrKc1*((1.e+0_fp/T)-(1.e+0_fp/298.e+0_fp)))
      Kc2_T = Kc2*exp(DhrKc1*((1.e+0_fp/T)-(1.e+0_fp/298.e+0_fp)))

      !CO2 dissolution  
      Hco2eff = Hco2_T*(1.e+0_fp+(Kc1_T/HPLUS)+((Kc1_T*Kc2_T)/
     &         (HPLUS*HPLUS)))
      xCO2 = 1.e+0_fp / ( 1.e+0_fp 
     &     + ( Hco2eff * 0.08205e+0_fp * T * LWC ) )
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
      REAL(fp),  INTENT(IN) :: T, P, LWC, HPLUS
!
! !OUTPUT PARAMETERS: 
!
      REAL(fp)              :: KCO2p, KCO2p2
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
      REAL(fp),  PARAMETER  :: Kc1 = 4.3e-7
      REAL(fp),  PARAMETER  :: Kc2 =4.68e-11
      REAL(fp),  PARAMETER  :: DhrKc1 = -1000.
      REAL(fp),  PARAMETER  :: DhrKc2 = -1760.
      REAL(fp),  PARAMETER  :: Hco2 = 3.4e-2
      REAL(fp),  PARAMETER  :: Dhco2 = -4.85/0.04
      ! CO2 concentration [v/v]
      REAL(fp),  PARAMETER  :: CO2 = 380.0e-6_fp 
!
! !LOCAL VARIABLES:
!
      REAL(fp)              :: Hco2_T, Kc1_T, Kc2_T
      REAL(fp)              :: Hco2eff, xCO2, pCO2

      !=================================================================
      ! kCO2 begins here!
      !=================================================================

      !CO2 dissolution constants
      Hco2_T = Hco2*exp(Dhco2*((1.e+0_fp/T)-(1.e+0_fp/298.e+0_fp)))
      Kc1_T = Kc1*exp(DhrKc1*((1.e+0_fp/T)-(1.e+0_fp/298.e+0_fp)))
      Kc2_T = Kc2*exp(DhrKc1*((1.e+0_fp/T)-(1.e+0_fp/298.e+0_fp)))

      !CO2 dissolution  
      Hco2eff = Hco2_T*(1.e+0_fp+(Kc1_T/HPLUS)+((Kc1_T*Kc2_T)/
     &         (HPLUS*HPLUS)))
      xCO2 = 1.e+0_fp / ( 1.e+0_fp 
     &     + ( Hco2eff * 0.08205e+0_fp * T * LWC ) )
      pCO2 = CO2 * P * xCO2

      KCO2p  = Hco2_T * Kc1_T * pCO2
      KCO2p2 = 2.e+0_fp * KCO2p *  Kc2_T

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
      REAL(fp),  INTENT(IN) :: T, P, LWC, HPLUS, SO2
!
! !OUTPUT PARAMETERS: 
!
      REAL(fp)              :: KSO2p, KSO2p2
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
      REAL(fp),  PARAMETER  :: Ks1 = 1.3e-2
      REAL(fp),  PARAMETER  :: Ks2 = 6.6e-8
      REAL(fp),  PARAMETER  :: Hso2 = 1.23
      REAL(fp),  PARAMETER  :: Dhso2 = -6.25/0.04
      REAL(fp),  PARAMETER  :: DhrKso21 = 1960.
      REAL(fp),  PARAMETER  :: DhrKso22 = 1500.
!
! !LOCAL VARIABLES:
!
      REAL(fp)              :: Hso2_T, Ks1_T, Ks2_T
      REAL(fp)              :: Hso2eff, xSO2, pSO2

      !=================================================================
      ! kSO2 begins here!
      !=================================================================


      ! SO2 dissolution constants
      Hso2_T = Hso2*exp(Dhso2*((1.e+0_fp/T)-(1.e+0_fp/298.e+0_fp)))
      Ks1_T = Ks1*exp(DhrKso21*((1.e+0_fp/T)-(1.e+0_fp/298.e+0_fp)))
      Ks2_T = Ks2*exp(DhrKso22*((1.e+0_fp/T)-(1.e+0_fp/298.e+0_fp)))

      ! SO2 dissolution
      Hso2eff = Hso2_T*(1.e+0_fp+(Ks1_T/HPLUS)+((Ks1_T*Ks2_T)/
     &         (HPLUS*HPLUS)))
      xSO2 = 1.e+0_fp / ( 1.e+0_fp 
     &     + ( Hso2eff * 0.08205e+0_fp * T * LWC ) )
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
      REAL(fp),  INTENT(IN) :: T, P, LWC, HPLUS, SO2
!
! !OUTPUT PARAMETERS: 
!
      REAL(fp)              :: KSO2p, KSO2p2
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
      REAL(fp),  PARAMETER  :: Ks1 = 1.3e-2
      REAL(fp),  PARAMETER  :: Ks2 = 6.6e-8
      REAL(fp),  PARAMETER  :: Hso2 = 1.23
      REAL(fp),  PARAMETER  :: Dhso2 = -6.25/0.04
      REAL(fp),  PARAMETER  :: DhrKso21 = 1960.
      REAL(fp),  PARAMETER  :: DhrKso22 = 1500.
!
! !LOCAL VARIABLES:
!
      REAL(fp)              :: Hso2_T, Ks1_T, Ks2_T
      REAL(fp)              :: Hso2eff, xSO2, pSO2

      !=================================================================
      ! kSO2 begins here!
      !=================================================================


      ! SO2 dissolution constants
      Hso2_T = Hso2*exp(Dhso2*((1.e+0_fp/T)-(1.e+0_fp/298.e+0_fp)))
      Ks1_T  = Ks1 *exp(DhrKso21*((1.e+0_fp/T)-(1.e+0_fp/298.e+0_fp)))
      Ks2_T  = Ks2 *exp(DhrKso22*((1.e+0_fp/T)-(1.e+0_fp/298.e+0_fp)))

      !SO2 dissolution
      Hso2eff = Hso2_T*(1.e+0_fp+(Ks1_T/HPLUS)+((Ks1_T*Ks2_T)/
     &         (HPLUS*HPLUS)))
      xSO2 = 1.e+0_fp / ( 1.e+0_fp 
     &     + ( Hso2eff * 0.08205e+0_fp * T * LWC ) )
      pSO2 = SO2 * P * xSO2

      KSO2p  = Hso2_T * Ks1_T * pSO2
      KSO2p2 = 2.e+0_fp * KSO2p * Ks2_T

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
      REAL(fp),  INTENT(IN) :: T, P, LWC, HPLUS, HNO3
!
! !OUTPUT PARAMETERS: 
!
      REAL(fp)              :: KHNO3p
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
      REAL(fp),  PARAMETER  :: Kn1 = 15.4
      REAL(fp),  PARAMETER  :: Hhno3 = 2.1e5
      REAL(fp),  PARAMETER  :: Dhhno3 = 0.
      REAL(fp),  PARAMETER  :: DhrKn1 = 8700.
!
! !LOCAL VARIABLES:
!
      REAL(fp)              :: Hhno3_T, Kn1_T
      REAL(fp)              :: Hhno3eff, xHNO3, pHNO3

      !=================================================================
      ! kHNO3 begins here!
      !=================================================================

      ! HNO3 dissolution constants
      Hhno3_T = Hhno3*exp(Dhhno3*((1.e+0_fp/T)-(1.e+0_fp/298.e+0_fp)))
      Kn1_T = Kn1*exp(DhrKn1*((1.e+0_fp/T)-(1.e+0_fp/298.e+0_fp)))

      !HNO3 dissolution
      Hhno3eff = 3.2e6/HPLUS
      xHNO3 = 1.e+0_fp / ( 1.e+0_fp 
     &      + ( Hhno3eff * 0.08205e+0_fp * T * LWC ) )
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
      REAL(fp),  INTENT(IN) :: T, P, LWC, HPLUS, Cl
!
! !OUTPUT PARAMETERS: 
!
      REAL(fp)              :: KHClp
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
      REAL(fp),  PARAMETER  :: Kcl = 1.74e+6_fp
      REAL(fp),  PARAMETER  :: Hcl = 1.1e+0_fp
      REAL(fp),  PARAMETER  :: Dhcl = -2.3158e+0_fp
      REAL(fp),  PARAMETER  :: DhrKcl = 6900.e+0_fp
!
! !LOCAL VARIABLES:
!
      REAL(fp)              :: Hcl_T, Kcl_T
      REAL(fp)              :: Hcleff, xCl, pHCl

      !=================================================================
      ! kHCl begins here!
      !=================================================================

      ! HNO3 dissolution constants
      HCl_T = Hcl*exp(Dhcl*((1.0e+0_fp/T)-(1.0e+0_fp/298.0e+0_fp)))
      Kcl_T = Kcl*exp(DhrKcl*((1.0e+0_fp/T)-(1.0e+0_fp/298.0e+0_fp)))

      !HCl dissolution
      Hcleff = Hcl_T*(1.0e+0_fp+(Kcl_T/HPLUS))
      xCl = 1.0e+0_fp / ( 1.0e+0_fp 
     &    + ( Hcleff * 0.08205e+0_fp * T * LWC ) )
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
      REAL(fp),  INTENT(IN) :: T, P, LWC, HPLUS, NH3, Kw
!
! !OUTPUT PARAMETERS: 
!
      REAL(fp)              :: KNH3p
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
      REAL(fp),  PARAMETER  :: Ka1 = 1.7e-5
      REAL(fp),  PARAMETER  :: Hnh3 = 62.
      REAL(fp),  PARAMETER  :: Dhnh3 = -8.17/0.04
      REAL(fp),  PARAMETER  :: DhrKa1 = -450.

      ! Variables
      REAL(fp)              :: Hnh3_T, Ka1_T
      REAL(fp)              :: Hnh3eff, xNH3, pNH3

      !=================================================================
      ! kNH3 begins here!
      !=================================================================

      !NH3 dissolution constants
      Hnh3_T = Hnh3*exp(Dhnh3*((1.e+0_fp/T)-(1.e+0_fp/298.e+0_fp)))
      Ka1_T = Ka1*exp(DhrKa1*((1.e+0_fp/T)-(1.e+0_fp/298.e+0_fp)))

      !NH3 dissolution
      Hnh3eff = Hnh3_T*(1.e+0_fp+((Ka1_T* HPLUS) / Kw))
      xNH3 = 1.e+0_fp / ( 1.e+0_fp 
     &     + ( Hnh3eff * 0.08205e+0_fp * T * LWC ) )
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
      REAL(f8)            :: A2, A1, A0
      REAL(fp)            :: CRUTES(3)
!
! !REMARKS:
! Formulae can be found in numer. recip.  on page 145
!   kiran  developed  this version on 25/4/1990
!   Dr. Francis S. Binkowski modified the routine on 6/24/91, 8/7/97
! ***
! *** modified 2/23/98 by fsb to incorporate Dr. Ingmar Ackermann's
!     recommendations for setting a0, a1,a2 as real(fp) variables.
!
! Modified by Bob Yantosca (10/15/02) 
! - Now use upper case / white space
! - force double precision with "D" exponents
! - updated comments / cosmetic changes
! - now call ERROR_STOP from "error_mod.f" to stop the run safely
!
! !REVISION HISTORY:
!  25 Jan 2012 - M. Payer    - Added ProTeX headers
!  06 Jan 2015 - M. Yannetti - Manual changes of some variables to f8
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
      REAL(fp), PARAMETER :: ONE    = 1.0e+0_fp
      REAL(fp), PARAMETER :: SQRT3  = 1.732050808e+0_fp
      REAL(fp), PARAMETER :: ONE3RD = 0.333333333e+0_fp
!
! !LOCAL VARIABLES:
!
      REAL(fp)            :: A2SQ,  THETA 
      REAL(f8)            :: PART1, PART2, PART3, YY1
      REAL(fp)            :: YY2,   YY3,   COSTH, SINTH
      REAL(f8)            :: PHI,   DUM1,  DUM2,  RRSQ
      REAL(f8)            :: QQ,    RR

      !=================================================================
      ! CUBIC begins here!
      !=================================================================
      A2SQ = A2 * A2
      QQ   = ( A2SQ - 3.e+0_f8 * A1 ) / 9.e+0_f8
      RR   = ( A2*( 2.e+0_f8*A2SQ - 9.e+0_f8*A1 ) + 27.e+0_f8*A0 ) 
     &       / 54.e+0_f8

      ! CASE 1 THREE REAL ROOTS or  CASE 2 ONLY ONE REAL ROOT
      DUM1 = QQ * QQ * QQ
      RRSQ = RR * RR
      DUM2 = DUM1 - RRSQ

      IF ( DUM2 .GE. 0.e+0_f8 ) THEN

         ! Now we have three real roots
         PHI = SQRT( DUM1 )

         IF ( ABS( PHI ) .LT. 1.e-20_f8 ) THEN
            CRUTES(1) = 0.0e+0_fp
            CRUTES(2) = 0.0e+0_fp
            CRUTES(3) = 0.0e+0_fp
            NR        = 0
         ENDIF
         
         THETA = ACOS( RR / PHI ) / 3.0e+0_fp
         COSTH = COS( THETA )
         SINTH = SIN( THETA )

         ! Use trig identities to simplify the expressions
         ! Binkowski's modification
         PART1     = SQRT( QQ )
         YY1       = PART1 * COSTH
         YY2       = YY1 - A2/3.0e+0_fp
         YY3       = SQRT3 * PART1 * SINTH
         CRUTES(3) = -2.0e+0_fp*YY1 - A2/3.0e+0_fp
         CRUTES(2) = YY2 + YY3
         CRUTES(1) = YY2 - YY3

         ! Set negative roots to a large positive value
         IF ( CRUTES(1) .LT. 0.0e+0_fp ) CRUTES(1) = 1.0e+9_fp
         IF ( CRUTES(2) .LT. 0.0e+0_fp ) CRUTES(2) = 1.0e+9_fp
         IF ( CRUTES(3) .LT. 0.0e+0_fp ) CRUTES(3) = 1.0e+9_fp

         ! Put smallest positive root in crutes(1)
         CRUTES(1) = MIN( CRUTES(1), CRUTES(2), CRUTES(3) )
         NR        = 3

      ELSE  

         ! Now here we have only one real root
         PART1     = SQRT( RRSQ - DUM1 )
         PART2     = ABS( RR )
         PART3     = ( PART1 + PART2 )**ONE3RD
         CRUTES(1) = -SIGN(ONE,RR) * ( PART3 + (QQ/PART3) ) - A2/3.e0_fp
         CRUTES(2) = 0.e+0_fp
         CRUTES(3) = 0.e+0_fp
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
      REAL(fp), INTENT(IN)  :: LWC     ! Liq water content [m3/m3]=1.E-6*L [g/m3]
      REAL(fp), INTENT(IN)  :: T       ! Temperature [K]
      REAL(fp), INTENT(IN)  :: P       ! Pressure [atm]
      REAL(fp), INTENT(IN)  :: SO2     ! SO2  mixing ratio [v/v]
      REAL(fp), INTENT(IN)  :: H2O2    ! H2O2 mixing ratio [v/v]
      REAL(fp), INTENT(IN)  :: O3      ! O3   mixing ratio [v/v]
      REAL(fp), INTENT(IN)  :: HPLUS   ! Concentration of H+ ion (i.e. pH) [v/v]
!
! !OUTPUT PARAMETERS: 
!
      REAL(fp), INTENT(OUT) :: KaqH2O2 ! Reaction rate for H2O2
      REAL(fp), INTENT(OUT) :: KaqO3   ! Reaction rate for O3
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
      REAL(fp), PARAMETER   :: R = 0.08205e+0_fp 
!
! !LOCAL VARIABLES:
!
      REAL(fp)              :: KH2O2,  RA,     KS1, KS2,    HCSO2
      REAL(fp)              :: FHCSO2, XSO2G,  SIV, HSO3,   XSO2AQ
      REAL(fp)              :: XHSO3,  XSO3,   KH1, HCH2O2, FHCH2O2
      REAL(fp)              :: XH2O2G, H2O2aq, KO0, KO1,    KO2
      REAL(fp)              :: HCO3,   XO3g,   O3aq

      !=================================================================
      ! AQCHEM_SO2 begins here!
      !
      ! Aqueous reaction rate
      ! HSO3- + H2O2 + H+ => SO4-- + 2H+ + H2O [Jacob, 1986]
      !=================================================================

      ! [Jacob, 1986]
      KH2O2 = 6.31e+14_fp * EXP( -4.76e+3_fp / T )  

!      ! [Jacobson, 1999]
!      KH2O2 = 7.45e+0_fp7 * EXP( -15.96e+0_fp * ( (298.15/T) - 1.) ) / 
!     &                  ( 1.e+0_fp + 13.e+0_fp * Hplus)

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
      Ks1    = 1.30e-2_fp * EXP( 6.75e+0_fp 
     &                      * ( 298.15e+0_fp / T - 1.e+0_fp ) )
      Ks2    = 6.31e-8_fp * EXP( 5.05e+0_fp 
     &                      * ( 298.15e+0_fp / T - 1.e+0_fp ) )

      ! SIV Fraction
      XSO2aq = 1.e+0_fp/(1.e+0_fp + Ks1/Hplus + Ks1*Ks2/(Hplus*Hplus))
      XHSO3  = 1.e+0_fp/(1.e+0_fp + Hplus/Ks1 + Ks2/Hplus)
      XSO3   = 1.e+0_fp/(1.e+0_fp + Hplus/Ks2 + Hplus*Hplus/(Ks1*Ks2))

      ! Henry's constant [mol/l-atm] and Effective Henry's constant for SO2
      HCSO2  = 1.22e+0_fp * EXP( 10.55e+0_fp 
     &                      * ( 298.15e+0_fp / T - 1.e+0_fp) )         
      FHCSO2 = HCSO2 * (1.e+0_fp + (Ks1/Hplus) + 
     &        (Ks1*Ks2 / (Hplus*Hplus)))
      
      XSO2g  = 1.e+0_fp / ( 1.e+0_fp + ( FHCSO2 * R * T * LWC ) )
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
      Kh1 = 2.20e-12_fp * EXP( -12.52e+0_fp * 
     &     ( 298.15e+0_fp / T - 1.e+0_fp ) )

      ! Henry's constant [mol/l-atm] and Effective Henry's constant for H2O2
      ! [Seinfeld and Pandis, 1998]
      ! HCH2O2  = 7.45D4 * EXP( 24.48e+0_fp * ( 298.15e+0_fp / T - 1.e+0_fp) ) 

      ! [Jacobson,1999]
      HCH2O2  = 7.45e+4_fp * EXP( 22.21e+0_fp * 
     &         (298.15e+0_fp / T - 1.e+0_fp) )
      FHCH2O2 = HCH2O2 * (1.e+0_fp + (Kh1 / Hplus))

      XH2O2g  = 1.e+0_fp / ( 1.e+0_fp + ( FHCH2O2 * R * T * LWC ) )
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
      KO0 = 2.40e+4_fp
      KO1 = 3.49e+12_fp * EXP( -4.83e+3_fp / T )  
      KO2 = 7.32e+14_fp * EXP( -4.03e+3_fp / T )

      !=================================================================
      ! H2O2 equilibrium reaction
      ! O3 + H2O = O3.H2O
      !  HCO3  = 1.13E-2 * EXP( 8.51 * (298.15/T -1.) ), S & P
      !  HCO3  = 1.13E-2 * EXP( 7.72 * (298.15/T -1.) ), Jacobson
      !=================================================================

      ! Calculate Henry's Law constant for atmospheric temperature
      HCO3  = 1.13e-2_fp * EXP( 8.51e+0_fp 
     &      * ( 298.15e+0_fp / T - 1.e+0_fp ) )

      XO3g  = 1.e+0_fp / ( 1.e+0_fp + ( HCO3 * R * T * LWC ) )
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
      SUBROUTINE HET_DROP_CHEM( I,    J,    L,      LSTOT, SSCvv, 
     &                          aSO4, GNH3, SO2_sr, H2O20, GNO3, 
     &                          SR,   Input_Opt,    State_Met    )
!
! !USES:
!
      USE ERROR_MOD,          ONLY : IT_IS_FINITE, GEOS_CHEM_STOP
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE TIME_MOD,           ONLY : GET_TS_CHEM
      USE TRACERID_MOD,       ONLY : IDTSO4, IDTSALC
      USE TRACERID_MOD,       ONLY : IDTDST2, IDTDST3, IDTDST4
!
! !INPUT PARAMETERS: 
!
      INTEGER,        INTENT(IN)  :: I, J, L
      REAL(fp),         INTENT(IN)  :: LSTOT
      REAL(fp),         INTENT(IN)  :: SSCvv
      REAL(fp),         INTENT(IN)  :: aSO4
      REAL(fp),         INTENT(IN)  :: GNH3
      REAL(fp),         INTENT(IN)  :: SO2_sr
      REAL(fp),         INTENT(IN)  :: H2O20
      REAL(fp),         INTENT(IN)  :: GNO3
      TYPE(OptInput), INTENT(IN)  :: Input_Opt   ! Input Options object
      TYPE(MetState), INTENT(IN)  :: State_Met   ! Meteorology State object
!
! !OUTPUT PARAMETERS: 
!
      REAL(fp),         INTENT(OUT) :: SR ! Sulfate production rate
!
! !REVISION HISTORY:
!  25 Jan 2012 - M. Payer    - Added ProTeX headers
!  05 Sep 2013 - M. Sulprizio- Now pass met fields using the State_Met object
!  06 Nov 2014 - R. Yantosca - Now use State_Met%AIRDEN(I,J,L)
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
      REAL(fp), PARAMETER   :: PI    = 3.14159e+0_fp
      REAL(fp), PARAMETER   :: SS_DEN = 2200.e+0_fp !dry sea-salt density [kg/m3]
      ! sigma of the size distribution for sea-salt (Jaegle et al., 2011)
      REAL(fp), PARAMETER   :: SIG_S = 1.8e+0_fp
      ! geometric dry mean diameters [m] for computing lognormal size distribution
      REAL(fp), PARAMETER   :: RG_S = 0.4e-6_fp !(Jaegle et a., 2011)
      REAL(fp), PARAMETER   :: RG_D2 = 1.5e-6_fp!(Ginoux et al., 2001)
      REAL(fp), PARAMETER   :: RG_D3 = 2.5e-6_fp
      REAL(fp), PARAMETER   :: RG_D4 = 4.e-6_fp
!
! !LOCAL VARIABLES:
!
      REAL(fp)              :: alpha_NH3, alpha_SO2, alpha_H2O2
      REAL(fp)              :: alpha_HNO3, alpha_B, alpha_CN
      REAL(fp)              :: alpha_W, alpha_SO4, sum_gas, H
      REAL(fp)              :: NDss, NDd2, NDd3, NDd4, CN, DEN, REFF, W
      REAL(fp)              :: DTCHEM, APV, DSVI
      REAL(fp)              :: B, NH3, SO2, H2O2, HNO3, SO4
      REAL(fp)              :: CNss, CNd2, CNd3, CNd4

      ! Pointers
      REAL(fp), POINTER     :: AD(:,:,:)
      REAL(fp), POINTER     :: AIRDEN(:,:,:)
      REAL(fp), POINTER     :: AIRVOL(:,:,:)
      REAL(fp), POINTER     :: U(:,:,:)
      REAL(fp), POINTER     :: V(:,:,:)

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
      NH3  = GNH3  * 1.0e+12_fp
      SO2  = SO2_sr  * 1.0e+12_fp
      H2O2 = H2O20 * 1.0e12_fp
      HNO3 = GNO3 * 1.0e12_fp

      ! Convert sulfate aerosol concentrations from [v/v] to [ug/m3]
      SO4 = ( aSO4                   * AD(I,J,L)     * 1.0e+9_fp ) / 
     &      ( Input_Opt%TCVV(IDTSO4) * AIRVOL(I,J,L) )

      ! Convert in cloud sulfate production rate from [v/v/timestep] to
      ! [ug/m3/timestep]
      B  = ( LSTOT                  * AD(I,J,L)     * 1.0e+9_fp ) / 
     &     ( Input_Opt%TCVV(IDTSO4) * AIRVOL(I,J,L) )

      ! Convert coarse-mode aerosol concentrations from [v/v] to [#/cm3]
      ! based on equation in Hofmann, Science, 1990.
      ! First convert from [v/v] to [kg/m3 air]
      CNss = SSCvv                   * AD(I,J,L)     /
     &     ( Input_Opt%TCVV(IDTSALC) * AIRVOL(I,J,L) )

      ! Now convert from [kg/m3 air] to [#/cm3 air]
      ! Sea-salt
      NDss = ( (3.e+0_fp/4.e+0_fp) * CNss ) / 
     &       (PI * SS_DEN * RG_S**3.e+0_fp * 
     &       exp( (9.e+0_fp/2.e+0_fp) * 
     &       (LOG(SIG_S)) ** 2.e+0_fp ) )

      ! Total coarse mode number concentration [#/cm3]
      CN = NDss ! sea-salt 

      ! Determine regression coefficients based on the local SO2 concentration
      IF ( SO2 <= 200.0e+0_fp ) THEN
         alpha_B    = 0.5318e+0_fp
         alpha_NH3  = -1.67e-7_fp
         alpha_SO2  = 2.59e-6_fp
         alpha_H2O2 = -1.77e-7_fp
         alpha_HNO3 = -1.72e-7_fp
         alpha_W    = 1.22e-6_fp
         alpha_CN   = 4.58e-6_fp
         alpha_SO4  = -1.00e-5_fp
      ELSE IF ( SO2 > 200.0e+0_fp .AND. SO2 <= 500.0e+0_fp ) THEN
         alpha_B    = 0.5591e+0_fp
         alpha_NH3  = 3.62e-6_fp
         alpha_SO2  = 1.66e-6_fp
         alpha_H2O2 = 1.06e-7_fp
         alpha_HNO3 = -5.45e-7_fp
         alpha_W    = -5.79e-7_fp
         alpha_CN   = 1.63e-5_fp
         alpha_SO4  = -7.40e-6_fp
      ELSE IF ( SO2 > 500.0e+0_fp .AND. SO2 < 1000.0e+0_fp ) THEN
         alpha_B    = 1.1547e+0_fp
         alpha_NH3  = -4.28e-8_fp
         alpha_SO2  = -1.23e-7_fp
         alpha_H2O2 = -9.05e-7_fp
         alpha_HNO3 = 1.73e-7_fp
         alpha_W    = 7.22e-6_fp
         alpha_CN   = 2.44e-5_fp
         alpha_SO4  = 3.25e-5_fp
      ELSE IF ( SO2 >= 1000.0e+0_fp ) THEN
         alpha_B    = 1.1795e+0_fp
         alpha_NH3  = 2.57e-7_fp
         alpha_SO2  = -5.54e-7_fp
         alpha_H2O2 = -1.08e-6_fp
         alpha_HNO3 = 1.95e-6_fp
         alpha_W    = 6.14e-6_fp
         alpha_CN   = 1.64e-5_fp
         alpha_SO4  = 2.48e-6_fp
      ENDIF

      ! Updraft velocity over the oceans [cm/s]
      W = 500e+0_fp

      ! Compute H (integration time interval * air parcel velocity) [m]
      ! DTCHEM is the chemistry timestep in seconds
      DTCHEM = GET_TS_CHEM() * 60e+0_fp

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
      SR = SR * Input_Opt%TCVV(IDTSO4) * 1.e-9_fp / AIRDEN(I,J,L)
	
      ! Don't allow SR to be negative
      SR = MAX( SR, 0.e+0_fp )

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
      USE CHEMGRID_MOD,       ONLY : ITS_IN_THE_NOCHEMGRID
      USE CMN_SIZE_MOD
      USE CMN_DIAG_MOD
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE GIGC_State_Chm_Mod, ONLY : ChmState
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE TIME_MOD,           ONLY : GET_TS_CHEM
      USE TRACERID_MOD,       ONLY : IDTSO4
      USE TRACERID_MOD,       ONLY : IDTSO4s
      USE TRACERID_MOD,       ONLY : IDTAS
      USE TRACERID_MOD,       ONLY : IDTAHS 
      USE TRACERID_MOD,       ONLY : IDTLET
      USE TRACERID_MOD,       ONLY : IDTSO4aq
      USE TRACERID_MOD,       ONLY : IDTNH4aq
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
!  The only production is from SO2 oxidation (save in CHEM_SO2).  Dry
!  deposition is now handled in mixing_mod.F90, so we can must add the
!  production from SO2 into the SO4 tracers.
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
!  12 Jun 2015 - R. Yantosca - Now remove orphaned ND44 variables
!  12 Jun 2015 - R. Yantosca - Drydep is now handled in mixing_mod.F90,
!                              so we can greatly collapse this code
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER           :: I,   J,    L,    N
      REAL(fp)          :: SO4, SO40, SO4s, SO40s

      ! Pointers
      REAL(fp), POINTER :: STT(:,:,:,:)

      !=================================================================
      ! CHEM_SO4 begins here!
      !=================================================================

      ! Return if tracers are not defined
      IF ( IDTSO4 == 0 .or. IDTSO4s == 0 ) RETURN
      IF ( DRYSO4 == 0 .or. DRYSO4s == 0 ) RETURN

      ! Assume success
      RC        = GIGC_SUCCESS

      ! Initialize pointers
      STT      => State_Chm%Tracers

      ! Loop over chemistry grid boxes
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, SO4, SO4s, SO40, SO40s )
!$OMP+SCHEDULE( DYNAMIC ) 
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Initialize for safety's sake
         SO4  = 0e+0_fp
         SO4s = 0e+0_fp

         ! Skip non-chemistry boxes
         IF ( ITS_IN_THE_NOCHEMGRID( I, J, L, State_Met ) ) CYCLE

         !==============================================================
         ! Initial concentrations before chemistry
         !==============================================================

         ! SO4 [v/v]
         SO40  = STT(I,J,L,IDTSO4)
         
         ! SO4 within coarse seasalt aerosol [v/v]
         SO40s = STT(I,J,L,IDTSO4s)

         !==============================================================
         ! SO4 chemistry
         !==============================================================

         ! SO4 production from SO2 [v/v/timestep]
         SO4 = SO40 + PSO4_SO2(I,J,L)

         !==============================================================
         ! SO4s (SO4 w/in seasalt aerosol) chemistry: 
         !==============================================================

         ! SO4 production from SO2 [v/v/timestep]
         SO4s = SO40s + PSO4_ss(I,J,L)

         !==============================================================
         ! Final concentrations after chemistry
         !==============================================================

         ! Error check
         IF ( SO4  < SMALLNUM ) SO4  = 0e+0_fp
         IF ( SO4s < SMALLNUM ) SO4s = 0e+0_fp

         ! Final concentrations [v/v]
         STT(I,J,L,IDTSO4)  = SO4
         STT(I,J,L,IDTSO4s) = SO4s

      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! Free pointers
      NULLIFY( STT )

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
      USE CHEMGRID_MOD,       ONLY : ITS_IN_THE_NOCHEMGRID
      USE CMN_SIZE_MOD
      USE DAO_MOD,            ONLY : CONVERT_UNITS 
      USE GIGC_Input_Opt_Mod, ONLY : OptINput 
      USE GIGC_State_Chm_Mod, ONLY : ChmState
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE TOMAS_MOD,          ONLY : AQOXID, GETACTBIN
! JKodros (6/2/15 - Allow for different TOMAS bin lengths)
# if   defined( TOMAS12 )
      USE TRACERID_MOD,       ONLY : IDTSO4, IDTNK5
# elif defined( TOMAS15 )
      USE TRACERID_MOD,       ONLY : IDTSO4, IDTNK8
# elif defined( TOMAS30 )
      USE TRACERID_MOD,       ONLY : IDTSO4, IDTNK10
# else
      USE TRACERID_MOD,       ONLY : IDTSO4, IDTNK20
# endif
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
      INTEGER           :: I, J, L
      integer           :: k, binact1, binact2
      INTEGER           :: KMIN
      REAL(fp)          :: SO4OXID

      ! For values from Input_Opt
      INTEGER           :: N_TRACERS
      REAL(fp)          :: TCVV(Input_Opt%N_TRACERS)
      REAL(fp)          :: XNUMOL(Input_Opt%N_TRACERS)

      ! Pointers
      REAL(fp), POINTER :: AD(:,:,:)

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
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR
                  
         ! Skip non-chemistry boxes
         IF ( ITS_IN_THE_NOCHEMGRID( I, J, L, State_Met ) ) CYCLE

         SO4OXID = PSO4_SO2AQ(I,J,L) * AD(I,J,L)/ TCVV(IDTSO4)

         IF ( SO4OXID > 0e+0_fp ) THEN
! JKodros (6/2/15 - Set activating bin based on which TOMAS bin length being used)
# if   defined( TOMAS12 )
                CALL GETACTBIN
     &         ( I, J, L, IDTNK5, .TRUE. , BINACT1, State_Chm )

            CALL GETACTBIN
     &         ( I, J, L, IDTNK5, .FALSE., BINACT2, State_Chm )
# elif defined ( TOMAS15 )
                CALL GETACTBIN
     &         ( I, J, L, IDTNK8, .TRUE. , BINACT1, State_Chm )

            CALL GETACTBIN
     &         ( I, J, L, IDTNK8, .FALSE., BINACT2, State_Chm )
# elif defined ( TOMAS30 )
            CALL GETACTBIN
     &         ( I, J, L, IDTNK10, .TRUE. , BINACT1, State_Chm )

            CALL GETACTBIN
     &         ( I, J, L, IDTNK10, .FALSE., BINACT2, State_Chm )
# else
            CALL GETACTBIN
     &         ( I, J, L, IDTNK20, .TRUE. , BINACT1, State_Chm )

            CALL GETACTBIN
     &         ( I, J, L, IDTNK20, .FALSE., BINACT2, State_Chm )
# endif

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
      USE CHEMGRID_MOD,       ONLY : ITS_IN_THE_NOCHEMGRID
      USE CMN_GCTM_MOD
      USE CMN_DIAG_MOD
      USE CMN_SIZE_MOD
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE GIGC_State_Chm_Mod, ONLY : ChmState
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE TIME_MOD,           ONLY : GET_TS_CHEM
      USE TRACERID_MOD,       ONLY : IDTMSA
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
!  The only production is from DMS oxidation (saved in CHEM_DMS).
!  Dry deposition is now treaded in mixing_mod.F90.
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
!  12 Jun 2015 - R. Yantosca - Now remove orphaned ND44 variables
!  12 Jun 2015 - R. Yantosca - Drydep is now handled in mixing_mod.F90,
!                              so we can greatly collapse this code
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! Scalars
      INTEGER           :: I,    J,    L
      REAL(fp)          :: MSA0, MSA

      ! Pointers
      REAL(fp), POINTER :: STT(:,:,:,:)

      !=================================================================
      ! CHEM_MSA begins here!
      !=================================================================
      IF ( IDTMSA == 0 .or. DRYMSA == 0 ) RETURN

      ! Assume success
      RC  = GIGC_SUCCESS

      ! Initialize GEOS-Chem tracer array [kg] from Chemistry State object
      ! (mpayer, 12/6/12)
      STT => State_Chm%Tracers

      ! Loop over chemistry grid boxes
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, MSA0, MSA )
!$OMP+SCHEDULE( DYNAMIC )
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Skip non-chemistry boxes
         IF ( ITS_IN_THE_NOCHEMGRID( I, J, L, State_Met ) ) CYCLE      

         ! Initial MSA [v/v]
         MSA0 = STT(I,J,L,IDTMSA) 

         ! MSA production from DMS [v/v/timestep]
         MSA = MSA0 + PMSA_DMS(I,J,L)

         ! Final MSA [v/v]
         IF ( MSA < SMALLNUM ) MSA = 0e+0_fp
         STT(I,J,L,IDTMSA) = MSA

      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! Free pointer
      NULLIFY( STT )

      END SUBROUTINE CHEM_MSA
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
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE GIGC_State_Chm_Mod, ONLY : ChmState
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE PBL_MIX_MOD,        ONLY : GET_FRAC_UNDER_PBLTOP
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

!
! !REMARKS:
!  Dry deposition is now applied in mixing_mod.F90.  Therefore we can 
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
!  12 Jun 2015 - R. Yantosca - Now remove orphaned ND44 variables
!  12 Jun 2015 - R. Yantosca - Drydep is now handled in mixing_mod.F90,
!                              so we can greatly collapse this code
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! Scalars
      INTEGER           :: I,   J,    L,    N
      REAL(fp)          :: NIT, NITs, NIT0, NIT0s, F_UNDER_TOP

      ! Pointers
      REAL(fp), POINTER :: STT(:,:,:,:)

      !=================================================================
      ! CHEM_NIT begins here!
      !=================================================================

      ! Return if tracers are not defined
      IF ( IDTNIT == 0 .or. IDTNITs == 0 ) RETURN
      IF ( DRYNIT == 0 .or. DRYNITs == 0 ) RETURN

      ! Assume success
      RC  = GIGC_SUCCESS

      ! Initialize GEOS-Chem tracer array [kg] from Chemistry State object
      ! (mpayer, 12/6/12)
      STT => State_Chm%Tracers

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, NIT0, NIT0s, NIT, NITs, F_UNDER_TOP )
!$OMP+SCHEDULE( DYNAMIC )
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Initial NITRATE [v/v]
         NIT0  = STT(I,J,L,IDTNIT)

         ! Initial NITRATE w/in seasalt [v/v]
         NIT0s = STT(I,J,L,IDTNITs)

         ! Initialize variables
         NIT   = 0e+0_fp
         NITs  = 0e+0_fp

         ! Fraction of box (I,J,L) underneath the PBL top [unitless]
         F_UNDER_TOP = GET_FRAC_UNDER_PBLTOP( I, J, L )     

         ! Only apply drydep to boxes w/in the PBL
         IF ( F_UNDER_TOP > 0e+0_fp ) THEN 

            !===========================================================
            ! NITs chemistry
            !===========================================================

            ! NIT prod from HNO3 uptake on fine sea-salt [v/v/timestep]
            NITs = NIT0s + PNITs(I,J,L)

            ! Store final concentration in STT [v/v]
            STT(I,J,L,IDTNITs) = NITs
            
         ENDIF
      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! Free pointer
      NULLIFY( STT )

      END SUBROUTINE CHEM_NIT
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
      REAL(fp),  INTENT(IN)    :: FRAC_OF_PBL(:) ! 
      REAL(fp),  INTENT(IN)    :: EMISS(:)  
!     
! !OUTPUT PARAMETERS:
!             
      REAL(fp),  INTENT(INOUT) :: SULFATE(:) ! partitioned emissions
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
      REAL(fp)                 :: TSULFATE

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
!  18 Sep 2014 - M. Sulprizio- Now get OH for offline aerosol sim from HEMCO
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER             :: JLOOP
      REAL(fp)              :: OH_MOLEC_CM3
 
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
            OH_MOLEC_CM3 = 0e+0_fp
         ENDIF

      ELSE IF ( Input_Opt%ITS_AN_AEROSOL_SIM ) THEN

         !---------------------
         ! Offline simulation
         !---------------------

         ! 1-D grid box index for SUNCOS
         JLOOP = ( (J-1) * IIPAR ) + I

         ! Test for sunlight...
         IF ( State_Met%SUNCOS(I,J) > 0e+0_fp 
     &    .and. TCOSZ(I,J) > 0e+0_fp ) THEN

            ! Impose a diurnal variation on OH during the day
            OH_MOLEC_CM3 = OH(I,J,L)                              * 
     &                     ( State_Met%SUNCOS(I,J) / TCOSZ(I,J) ) *
     &                     ( 1440e+0_fp                / GET_TS_CHEM() )

            ! OH is in kg/m3 (from HEMCO), convert to molec/cm3 (mps, 9/18/14)
            OH_MOLEC_CM3 = OH_MOLEC_CM3 * XNUMOL_OH / CM3PERM3

            ! Make sure OH is not negative
            OH_MOLEC_CM3 = MAX( OH_MOLEC_CM3, 0e+0_fp )
               
         ELSE

            ! At night, OH goes to zero
            OH_MOLEC_CM3 = 0e+0_fp

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
      REAL(fp),  INTENT(IN) :: OH      ! OH at (I,J,L) to be saved into CSPEC
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
!  24 Jul 2014 - R. Yantosca - Now compute BOXVL internally
!  18 Sep 2014 - M. Sulprizio- Now get NO3 for offline aerosol sims from HEMCO
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER :: JLOOP
      REAL(fp)  :: NO3_MOLEC_CM3
      REAL(fp)  :: BOXVL

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
            NO3_MOLEC_CM3 = 0e+0_fp
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
         IF ( State_Met%SUNCOS(I,J) > 0e+0_fp ) THEN

            ! NO3 goes to zero during the day
            NO3_MOLEC_CM3 = 0e+0_fp
              
         ELSE

            ! Grid box volume [cm3]
            BOXVL         = State_Met%AIRVOL(I,J,L) * 1e+6_fp

            ! At night: Get NO3 [v/v] and convert it to [kg]
            NO3_MOLEC_CM3 = NO3(I,J,L)     * State_Met%AD(I,J,L) *
     &                      ( 62e+0_fp/28.97e+0_fp ) 

            ! Convert NO3 from [kg] to [molec/cm3]
            NO3_MOLEC_CM3 = NO3_MOLEC_CM3  * XNUMOL_NO3 / BOXVL
               
         ENDIF
            
         ! Make sure NO3 is not negative
         NO3_MOLEC_CM3  = MAX( NO3_MOLEC_CM3, 0e+0_fp )

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
      REAL(fp),  INTENT(IN) :: NO3     ! NO3 at (I,J,L) to be saved into CSPEC
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
!  06 Nov 2014 - R. Yantosca - Now use State_Met%AIRDEN(I,J,L)
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER               :: JLOOP
      REAL(fp)                :: O3_VV
 
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
            O3_VV = ( CSPEC(JLOOP,IDO3)       * 1e+6_fp   ) / 
     &              ( State_Met%AIRDEN(I,J,L) * XNUMOLAIR )
         ELSE
            O3_VV = 0e+0_fp
         ENDIF

      ELSE IF ( Input_Opt%ITS_AN_AEROSOL_SIM ) THEN
         
         !--------------------
         ! Offline simulation
         !--------------------

         ! Get O3 [v/v] for this gridbox & month
         ! O3 data will only be defined within the chemistry grid
         IF ( L <= LLCHEM ) THEN
            O3_VV = O3m(I,J,L)
         ELSE
            O3_VV = 0e+0_fp
         ENDIF

      ELSE

         !--------------------
         ! Invalid simulation
         !--------------------         
         CALL ERROR_STOP( 'Invalid NSRCX!', 'GET_O3 (sulfate_mod.f)')

      ENDIF

      END FUNCTION GET_O3
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
      REAL(fp)              :: A0, A1, A2, A3, B1, B2, B3
      REAL(fp)              :: LHR0, R, AHR, DEC, TIMLOC, YMID_R
      REAL(fp)              :: SUNTMP(MAXIJ)
      
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
            SUNTMP(:) = 0e+0_fp

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
               TCOSZ(I,J) = TCOSZ(I,J) + MAX( SUNTMP(IJLOOP), 0e+0_fp )

               ! COSZM is the peak value of SUMTMP during a day at (I,J)
               ! (rjp, bmy, 3/30/04)
               COSZM(I,J) = MAX( COSZM(I,J), SUNTMP(IJLOOP) )

               ! TTDAY is the total daylight time at location (I,J)
               IF ( SUNTMP(IJLOOP) > 0e+0_fp ) THEN
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
! !IROUTINE: get_alk
!
! !DESCRIPTION: Subroutine GET\_ALK returns the seasalt alkalinity emitted at
!  each timestep to sulfate\_mod.f for chemistry on seasalt aerosols.
!  (bec, 12/7/04, 11/23/09)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE GET_ALK( am_I_Root, I,J,L, ALK1, ALK2, Kt1, Kt2, 
     &                    Kt1N, Kt2N, Input_Opt, State_Met, RC)
!
! !USES:
!
      USE GIGC_ErrCode_Mod
      USE ERROR_MOD,          ONLY : IT_IS_NAN, ERROR_STOP
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE PBL_MIX_MOD,        ONLY : GET_FRAC_OF_PBL, GET_PBL_TOP_L
      USE HCOI_GC_MAIN_MOD,   ONLY : GetHcoVal
      USE TRACERID_MOD,       ONLY : IDTSALA, IDTSALC
!
! !INPUT PARAMETERS: 
!
      LOGICAL,        INTENT(IN)    :: am_I_Root   ! Root CPU? 
      INTEGER,        INTENT(IN)    :: I, J, L     ! Lon-lat-alt indices
      TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
      TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology State object
!
! !OUTPUT PARAMETERS: 
!
      REAL(fp),         INTENT(OUT)   :: ALK1, ALK2  ! [kg]
      REAL(fp),         INTENT(OUT)   :: Kt1, Kt2, Kt1N, Kt2N ! [s-1]
!
! !INPUT/OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(INOUT) :: RC
!
! !REVISION HISTORY: 
!  (1 ) Becky Alexander says we can remove AREA1, AREA2 (bec, bmy, 9/5/06)
!  (2 ) Bug fix to remove a double-substitution.  Replace code lines for 
!        TERM{123}A, TERM{123}B, TERM{123}AN, TERM{123}BN. (bec, bmy, 7/18/08)
!  (3 ) Updated hygroscopic growth parameters (bec, bmy, 11/23/09)
!  22 Dec 2011 - M. Payer    - Added ProTeX headers
!  09 Nov 2012 - M. Payer    - Replaced all met field arrays with State_Met
!                              derived type object
!  25 Jun 2014 - R. Yantosca - Now accept Input_Opt via the arg list
!  25 Jun 2014 - R. Yantosca - Removed references to tracer_mod.F     
!  02 Nov 2014 - C. Keller   - Now get alkalinity and number density from HEMCO
!  04 Nov 2014 - C. Keller   - Moved from seasalt_mod.F
!  12 Nov 2014 - C. Keller   - Now weight # dens. and alk. by fraction of PBL
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL(fp),  PARAMETER :: PI = 3.14159265
      REAL(fp)             :: N1, N2, Kt
      REAL(fp)             :: HGF, ALK
      REAL(fp)             :: RAD1, RAD2, RAD3
      REAL(fp)             :: term1a, term2a, term3a
      REAL(fp)             :: term1b, term2b, term3b
      REAL(fp)             :: term1aN, term2aN, term3aN
      REAL(fp)             :: term1bN, term2bN, term3bN
      REAL(fp)             :: const1, const2, const1N, const2N
      REAL(fp)             :: a1, a2, b1, b2, a1N, a2N, b1N, b2N
      REAL(fp),  PARAMETER :: MINDAT = 1.e-20_fp
      INTEGER            :: IRH
      REAL(fp),  PARAMETER :: gamma_SO2 = 0.11e+0_fp !from Worsnop et al. (1989)
      REAL(fp),  PARAMETER :: gamma_HNO3 = 0.2e+0_fp !from JPL [2001] 
      REAL(fp),  PARAMETER :: Dg = 0.2e+0_fp !gas phase diffusion coeff. [cm2/s]
      REAL(fp),  PARAMETER :: v = 3.0e+4_fp  !cm/s

      ! HEMCO update
      REAL(fp)                       :: FEMIS
      INTEGER                      :: NTOP
      LOGICAL                      :: FOUND
 
      !=================================================================
      ! GET_ALK begins here!
      !=================================================================

      ! Zero variables
      KT1   = 0.e+0_fp
      KT2   = 0.e+0_fp
      KT1N  = 0.e+0_fp
      KT2N  = 0.e+0_fp
      N1    = 0.e+0_fp
      N2    = 0.e+0_fp
      ALK1  = 0.e+0_fp
      ALK2  = 0.e+0_fp

      !-----------------------------------------------------------------------
      ! Get alkalinity from HEMCO. This is just the current emissions,
      ! converted from kg/m2/s to kg. In the original seasalt code, the 
      ! alkalinity was set to the total surface flux for every layer within
      ! the PBL, and to zero above. This seems unrealistic, and in the code
      ! below the alkalinity (and number density) are scaled by the fraction
      ! of PBL. This approach ensures that the total number density and 
      ! alkalinity is preserved.
      !                                                 (ckeller, 10/31/2014)
      !-----------------------------------------------------------------------
      ! [kg] use this when not transporting alk
!      ALK1  = ALK_EMIS(I,J,L,1)
!      ALK2  = ALK_EMIS(I,J,L,2)

      ! Layer in which the PBL top occurs
      NTOP = CEILING( GET_PBL_TOP_L( I, J ) )
      
      ! Do the following only if we are within the PBL
      IF ( L <= NTOP ) THEN      

         ! Fraction of the PBL spanned by box (I,J,L) [unitless]
         FEMIS = GET_FRAC_OF_PBL( I, J, L )
         
         ! Uncomment the following line to reproduce pre-HEMCO code.
         ! FEMIS = 1.0d0

         ! Get ALK1 and ALK2 surface emissions from HEMCO. These are in 
         ! kg/m2/s. 
         CALL GetHcoVal( IDTSALA, I, J, 1, FOUND, Emis=ALK1 )
         CALL GetHcoVal( IDTSALC, I, J, 1, FOUND, Emis=ALK2 )

         ! kg/m2/s --> kg. Weight by fraction of PBL
         ALK1 = MAX(ALK1,0.0e+0_fp) 
     &        * State_Met%AREA_M2(I,J,1) * TS_EMIS * FEMIS 
         ALK2 = MAX(ALK2,0.0e+0_fp) 
     &        * State_Met%AREA_M2(I,J,1) * TS_EMIS * FEMIS 

         ! Get number densities in # / cm3. Weight by fraction of PBL
         IF ( ASSOCIATED(NDENS_SALA) ) THEN
            N1 = NDENS_SALA(I,J) / State_Met%AIRVOL(I,J,L) 
     &         * 1.e-6_fp        * FEMIS
         ENDIF
         IF ( ASSOCIATED(NDENS_SALC) ) THEN
            N2 = NDENS_SALC(I,J) / State_Met%AIRVOL(I,J,L) 
     &         * 1.e-6_fp        * FEMIS
         ENDIF

      ENDIF

      !-----------------------------------------------------------------------
      ! NOTE: If you want to transport alkalinity then uncomment this section
      ! (bec, bmy, 4/13/05)
      ! 
      !! alkalinity [v/v] to [kg] use this when transporting alk
      !! or using Liao et al [2004] assumption of a continuous supply of
      ! alkalinity based on Laskin et al. [2003]
      !ALK1 = STT(I,J,L,IDTSALA) * State_Met%AD(I,J,L)/TCVV(IDTSALA)
      !ALK2 = STT(I,J,L,IDTSALC) * State_Met%AD(I,J,L)/TCVV(IDTSALC)
      !-----------------------------------------------------------------------

      ! Conversion from [m-3] --> [cm-3]
!      N1 = N_DENS(I,J,L,1) * 1.d-6
!      N2 = N_DENS(I,J,L,2) * 1.d-6

      ALK = ALK1 + ALK2

      ! If there is any alkalinity ...
      IF ( ALK > MINDAT ) THEN

         ! set humidity index IRH as a percent
         IRH = State_Met%RH(I,J,L)
         IRH = MAX(  1, IRH )
         IRH = MIN( 99, IRH )

         ! Hygroscopic growth factor for sea-salt from Chin et al. (2002)
         ! Updated (bec, bmy, 11/23/09)
         IF ( IRH < 100 ) HGF = 4.8e+0_fp
         IF ( IRH < 99  ) HGF = 2.9e+0_fp
         IF ( IRH < 95  ) HGF = 2.4e+0_fp
         IF ( IRH < 90  ) HGF = 2.0e+0_fp
         IF ( IRH < 80  ) HGF = 1.8e+0_fp
         IF ( IRH < 70  ) HGF = 1.6e+0_fp
         IF ( IRH < 50  ) HGF = 1.0e+0_fp

         ! radius of sea-salt aerosol size bins [cm] accounting for 
         ! hygroscopic growth
         RAD1 = Input_Opt%SALA_REDGE_um(1) * HGF * 1.e-4_fp 
         RAD2 = Input_Opt%SALA_REDGE_um(2) * HGF * 1.e-4_fp 
         RAD3 = Input_Opt%SALC_REDGE_um(2) * HGF * 1.e-4_fp 

         !----------------------------------
         ! SO2 uptake onto fine particles 
         !----------------------------------

	 ! calculate gas-to-particle rate constant for uptake of 
	 ! SO2 onto fine sea-salt aerosols [Jacob, 2000] analytical solution
         CONST1 = 4.e+0_fp/(V*GAMMA_SO2)
         A1     = (RAD1/DG)+CONST1
         B1     = (RAD2/DG)+CONST1
         TERM1A = ((B1**2)/2.0e+0_fp) - ((A1**2)/2.0e+0_fp)
         TERM2A = 2.e+0_fp*CONST1*(B1-A1)
         TERM3A = (CONST1**2)*LOG(B1/A1)
         KT1    = 4.e+0_fp*PI*N1*(DG**3)*(TERM1A - TERM2A + TERM3A)

         !----------------------------------
         ! SO2 uptake onto coarse particles 
         !----------------------------------
         
	 ! calculate gas-to-particle rate constant for uptake of 
	 ! SO2 onto coarse sea-salt aerosols [Jacob, 2000] analytical solution
         CONST2 = 4.e+0_fp/(V*GAMMA_SO2)
         A2     = (RAD2/DG)+CONST2
         B2     = (RAD3/DG)+CONST2
         TERM1B = ((B2**2)/2.0e+0_fp) - ((A2**2)/2.0e+0_fp)
         TERM2B = 2.e+0_fp*CONST2*(B2-A2)
         TERM3B = (CONST2**2)*LOG(B2/A2)
         KT2    = 4.e+0_fp*PI*N2*(DG**3)*(TERM1B - TERM2B + TERM3B)
         KT     = KT1 + KT2

         !----------------------------------
         ! HNO3 uptake onto fine particles 
         !----------------------------------

         ! calculate gas-to-particle rate constant for uptake of 
         ! HNO3 onto fine sea-salt aerosols [Jacob, 2000] analytical solution
         CONST1N = 4.e+0_fp/(V*GAMMA_HNO3)
         A1N     = (RAD1/DG)+CONST1N
         B1N     = (RAD2/DG)+CONST1N
         TERM1AN = ((B1N**2)/2.0e+0_fp) - ((A1N**2)/2.0e+0_fp)
         TERM2AN = 2.e+0_fp*CONST1N*(B1N-A1N)
         TERM3AN = (CONST1N**2)*LOG(B1N/A1N)
         KT1N    = 4.e+0_fp*PI*N1*(DG**3)*(TERM1AN - TERM2AN + TERM3AN)

         !----------------------------------
         ! HNO3 uptake onto coarse particles 
         !----------------------------------

	 ! calculate gas-to-particle rate constant for uptake of 
	 ! HNO3 onto coarse sea-salt aerosols [Jacob, 2000] analytical solution
         CONST2N = 4.e+0_fp/(V*GAMMA_HNO3)
         A2N     = (RAD2/DG)+CONST2N
         B2N     = (RAD3/DG)+CONST2N
         TERM1BN = ((B2N**2)/2.0e+0_fp) - ((A2N**2)/2.0e+0_fp)
         TERM2BN = 2.e+0_fp*CONST2N*(B2N-A2N)
         TERM3BN = (CONST2N**2)*LOG(B2N/A2N)
         KT2N    = 4.e+0_fp*PI*N2*(DG**3)*(TERM1BN - TERM2BN + TERM3BN)


      ELSE

         ! If no alkalinity, set everything to zero
         ALK1 = 0.e+0_fp
         ALK2 = 0.e+0_fp
         KT1  = 0.e+0_fp
         KT1N = 0.e+0_fp
         KT2  = 0.e+0_fp
         KT2N = 0.e+0_fp

      ENDIF

      ! Return w/ success
      RC = GIGC_SUCCESS

      END SUBROUTINE GET_ALK
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
!  22 May 2015 - R. Yantosca - Remove variables made obsolete by HEMCO
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
      SSTEMP = 0e+0_fp

      ALLOCATE( DMSo( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'DMSo' )
      DMSo = 0e+0_fp

      ALLOCATE( PMSA_DMS( IIPAR, JJPAR, LLCHEM ), STAT=AS ) 
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'PMSA_DMS' )
      PMSA_DMS = 0e+0_fp

      ALLOCATE( PSO2_DMS( IIPAR, JJPAR, LLCHEM ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'PSO2_DMS' )
      PSO2_DMS = 0e+0_fp

      ALLOCATE( PSO4_SO2( IIPAR, JJPAR, LLCHEM ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'PSO4_SO2' )
      PSO4_SO2 = 0e+0_fp

      ALLOCATE( PSO4_ss( IIPAR, JJPAR, LLCHEM ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'PSO4_ss' )
      PSO4_ss = 0e+0_fp

      ALLOCATE( PNITs( IIPAR, JJPAR, LLCHEM ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'PNITs' )
      PNITs = 0e+0_fp

      ALLOCATE( SOx_SCALE( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'SOx_SCALE' )
      SOx_SCALE = 0e+0_fp

#if   !defined( GEOS_5 ) && !defined( MERRA ) && !defined( GEOS_FP )
      ! If we are using GEOS-5 or MERRA met, then get the cloud fraction 
      ! directly from the met fields.  (skim, bmy, 1/14/10)
      ALLOCATE( VCLDF( IIPAR, JJPAR, LLCHEM ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'VCLDF' )
      VCLDF = 0e+0_fp
#endif

#if   defined( TOMAS )
      ! Allocate for TOMAS microphysics
      ALLOCATE( PSO4_SO2AQ( IIPAR, JJPAR, LLCHEM ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'PSO4_SO2AQ' )
      PSO4_SO2AQ = 0e+0_fp

      ALLOCATE( SO4_ANTH( IIPAR, JJPAR, LLCHEM ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'SO4_ANTH' )
      SO4_ANTH = 0e+0_fp

      ALLOCATE( SO4_BIOF( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'SO4_BIOF' )
      SO4_BIOF = 0e+0_fp
#endif

      !=================================================================
      ! Only initialize the following for offline aerosol simulations
      !=================================================================
      IF ( Input_Opt%ITS_AN_AEROSOL_SIM ) THEN

         ALLOCATE( TCOSZ( IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'TCOSZ' )
         TCOSZ = 0e+0_fp

         ALLOCATE( TTDAY( IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'TTDAY' )
         TTDAY = 0e+0_fp
 
         ALLOCATE( COSZM( IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'COSZM' )
         COSZM = 0e+0_fp
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
      SUBROUTINE CLEANUP_SULFATE()  
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
!  04 Mar 2015 - R. Yantosca - Remove EEV, NEV.  Volcano eruptions are now
!                              handled via HEMCO.
!  22 May 2015 - R. Yantosca - Remove variables made obsolete by HEMCO
!EOP
!------------------------------------------------------------------------------
!BOC
!

      !=================================================================
      ! CLEANUP_SULFATE begins here!
      !=================================================================
      IF ( ALLOCATED( DMSo       ) ) DEALLOCATE( DMSo       )
      IF ( ALLOCATED( PMSA_DMS   ) ) DEALLOCATE( PMSA_DMS   )
      IF ( ALLOCATED( PNITs      ) ) DEALLOCATE( PNITs      )
      IF ( ALLOCATED( PSO2_DMS   ) ) DEALLOCATE( PSO2_DMS   )
      IF ( ALLOCATED( PSO4_SO2   ) ) DEALLOCATE( PSO4_SO2   )
#if   defined( TOMAS )
      IF ( ALLOCATED( PSO4_SO2AQ ) ) DEALLOCATE( PSO4_SO2AQ ) 
      IF ( ALLOCATED( SO4_ANTH   ))  DEALLOCATE( SO4_ANTH   )
      IF ( ALLOCATED( SO4_BIOF   ))  DEALLOCATE( SO4_BIOF   )
#endif
      IF ( ALLOCATED( PSO4_ss    ) ) DEALLOCATE( PSO4_ss    )
      IF ( ALLOCATED( SOx_SCALE  ) ) DEALLOCATE( SOx_SCALE  )
      IF ( ALLOCATED( SSTEMP     ) ) DEALLOCATE( SSTEMP     )
      IF ( ALLOCATED( TCOSZ      ) ) DEALLOCATE( TCOSZ      )
      IF ( ALLOCATED( TTDAY      ) ) DEALLOCATE( TTDAY      )          
      IF ( ALLOCATED( VCLDF      ) ) DEALLOCATE( VCLDF      )
      IF ( ALLOCATED( COSZM      ) ) DEALLOCATE( COSZM      )

      ! Free pointers
      IF ( ASSOCIATED( O3m        ) ) O3m        => NULL()
      IF ( ASSOCIATED( PH2O2m     ) ) PH2O2m     => NULL()
      IF ( ASSOCIATED( JH2O2      ) ) JH2O2      => NULL()
      IF ( ASSOCIATED( OH         ) ) OH         => NULL()
      IF ( ASSOCIATED( NO3        ) ) NO3        => NULL()
      IF ( ASSOCIATED( HNO3       ) ) HNO3       => NULL()
      IF ( ASSOCIATED( NDENS_SALA ) ) NDENS_SALA => NULL()
      IF ( ASSOCIATED( NDENS_SALC ) ) NDENS_SALC => NULL()

      END SUBROUTINE CLEANUP_SULFATE
!EOC
      END MODULE SULFATE_MOD
