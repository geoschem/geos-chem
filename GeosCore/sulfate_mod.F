! $Id: sulfate_mod.f,v 1.11 2010/03/15 19:33:20 ccarouge Exp $
      MODULE SULFATE_MOD
!
!******************************************************************************
!  Module SULFATE_MOD contains arrays and routines for performing either a
!  coupled chemistry/aerosol run or an offline sulfate aerosol simulation.
!  Original code taken from Mian Chin's GOCART model and modified accordingly.
!  (rjp, bdf, bmy, 6/22/00, 8/26/10)
!
!  Module Variables:
!  ============================================================================
!  (1 ) XNUMOL_OH  (REAL*8 ) : Molecules OH  per kg OH          [molec/kg]
!  (2 ) XNUMOL_O3  (REAL*8 ) : Molecules O3  per kg O3          [molec/kg]
!  (3 ) XNUMOL_NO3 (REAL*8 ) : Molecules NO3 per kg NO3         [molec/kg]
!  (4 ) TCVV_S     (REAL*8 ) : Ratio: Molwt air / Molwt S       [unitless]
!  (5 ) DMSo       (REAL*8 ) : DMS oceanic emissions            [v/v/timestep]
!  (6 ) DRYH2O2    (INTEGER) : Pointer to H2O2  in DEPVEL array [unitless] 
!  (7 ) DRYSO2     (INTEGER) : Pointer to SO2   in DEPVEL array [unitless]
!  (8 ) DRYSO4     (INTEGER) : Pointer to SO4   in DEPVEL array [unitless]
!  (9 ) DRYSO4s    (INTEGER) : Pointer to SO4s  in DEPVEL array [unitless]
!  (10) DRYMSA     (INTEGER) : Pointer to MSA   in DEPVEL array [unitless]
!  (11) DRYNH3     (INTEGER) : Pointer to NH3   in DEPVEL array [unitless]
!  (12) DRYNH4     (INTEGER) : Pointer to NH4   in DEPVEL array [unitless]
!  (13) DRYNIT     (INTEGER) : Pointer to NIT   in DEPVEL array [unitless]
!  (14) DRYNITs    (INTEGER) : Pointer to NITs  in DEPVEL array [unitless]
!  (15) DRYSO4aq   (INTEGER) : Pointer to SO4aq in DEPVEL array [unitless]
!  (16) DRYAS      (INTEGER) : Pointer to AS    in DEPVEL array [unitless]  
!  (17) DRYAHS     (INTEGER) : Pointer to AHS   in DEPVEL array [unitless]
!  (18) DRYLET     (INTEGER) : Pointer to LET   in DEPVEL array [unitless]
!  (19) DRYNH4aq   (INTEGER) : Pointer to NH4aq in DEPVEL array [unitless]
!  (20) ENH3_an    (REAL*8 ) : NH3 anthropogenic emissions      [kg NH3/box/s]
!  (21) ENH3_bb    (REAL*8 ) : NH3 biomass emissions            [kg NH3/box/s]
!  (22) ENH3_bf    (REAL*8 ) : NH3 biofuel emissions            [kg NH3/box/s]
!  (23) ENH3_na    (REAL*8 ) : NH73 natural source emissions    [kg NH3/box/s]
!  (24) ESO2_ac    (REAL*8 ) : SO2 aircraft emissions           [kg SO2/box/s]
!  (25) ESO2_an    (REAL*8 ) : SO2 anthropogenic emissions      [kg SO2/box/s]
!  (26) ESO2_ev    (REAL*8 ) : SO2 eruptive volcanic em.        [kg SO2/box/s]
!  (27) ESO2_nv    (REAL*8 ) : SO2 non-eruptive volcanic em.    [kg SO2/box/s]
!  (28) ESO2_bb    (REAL*8 ) : SO2 biomass burning emissions    [kg SO2/box/s]
!  (29) ESO2_bf    (REAL*8 ) : SO2 biofuel burning emissions    [kg SO2/box/s]
!  (30) ESO2_sh    (REAL*8 ) : SO2 ship emissions               [kg SO2/box/s]
!  (31) ESO4_an    (REAL*8 ) : SO4 anthropogenic emissions      [kg SO4/box/s]
!  (32) JH2O2      (REAL*8 ) : Monthly mean J(H2O2) values      [s-1]
!  (33) O3m        (REAL*8 ) : Monthly mean O3 concentration    [v/v]
!  (34) PH2O2m     (REAL*8 ) : Monthly mean P(H2O2)             [molec/cm3/s]
!  (35) PMSA_DMS   (REAL*8 ) : P(MSA) from DMS                  [v/v/timestep]
!  (36) PSO2_DMS   (REAL*8 ) : P(SO2) from DMS                  [v/v/timestep]
!  (37) PSO4_SO2   (REAL*8 ) : P(SO4) from SO2                  [v/v/timestep]
!  (38) SSTEMP     (REAL*8 ) : Sea surface temperatures         [K]
!  (39) VCLDF      (REAL*8 ) : Volume cloud frac. for SO2 aq.   [unitless]
!  (40) Eev        (REAL*8 ) : SO2 em. from eruptive volcanoes  [kg SO2/box/s]
!  (41) Env        (INTEGER) : SO2 em. from non-erup volcanoes  [kg SO2/box/s]
!  (42) TCOSZ      (REAL*8 ) : Sum of cos(SZA) for offline run  [unitless] 
!  (43) TTDAY      (REAL*8 ) : Total daylight length at (I,J)   [minutes]
!  (44) SMALLNUM   (REAL*8 ) : Small number - prevent underflow [unitless]
!  (45) COSZM      (REAL*8 ) : Array for MAX(cos(SZA)) at (I,J) [unitless]
!  (46) LVOLC      (INTEGER) : Number of volcanic levels (20)   [unitless]
!  
!  Module Routines:
!  ===========================================================================
!  (1 ) GET_VCLDF         : Computes volume cloud fraction for SO2 chemistry 
!  (2 ) GET_LWC           : Computes liquid water content as a function of T
!  (3 ) CHEMSULFATE       : Driver routine for sulfate/aerosol chemistry
!  (4 ) GRAV_SETTLING     : Routine to compute settling of SO4s and NITs
!  (5 ) CHEM_DMS          : Chemistry routine for DMS tracer
!  (6 ) CHEM_H2O2         : Chemistry routine for H2O2 tracer
!  (7 ) CHEM_SO2          : Chemistry routine for SO2 tracer
!  (8 ) SEASALT_CHEM      : Computes SO2->SO4 and HNO3->nitrate w/in seasalt
!  (9 ) AQCHEM_SO2        : Computes reaction rates for aqueous SO2 chemistry
!  (10) CHEM_SO4          : Chemistry routine for SO4 tracer
!  (11) PHASE_SO4         : Computes phase transition for crystalline tracers 
!  (12) PHASE_RADIATIVE   : Computes radiative forcing for crystalline tracers
!  (13) CHEM_MSA          : Chemistry routine for MSA tracer
!  (14) CHEM_NH3          : Chemistry routine for ammonia tracer
!  (15) CHEM_NH4          : Chemistry routine for ammonium tracer
!  (16) CHEM_NIT          : Chemistry routine for nitrates tracer
!  (17) EMISSSULFATE      : Driver routine for sulfate/aerosol emissions
!  (18) SRCDMS            : Emission routine for DMS tracer
!  (19) SRCSO2            : Emission routine for SO2 tracer
!  (20) SRCSO4            : Emission routine for SO4 tracer
!  (21) SRCNH3            : Emission routine for NH3 tracer
!  (22) GET_OH            : Returns OH for coupled or offline simulations
!  (23) SET_OH            : Resets modified OH in SMVGEAR's CSPEC array
!  (24) GET_NO3           : Returns NO3 for coupled or offline simulations
!  (25) SET_NO3           : Resets modified OH in SMVGEAR's CSPEC array
!  (26) GET_O3            : Returns O3 for coupled or offline simulations
!  (27) READ_NONERUP_VOLC : Reads SO2 emissions from non-eruptive volcanoes
!  (28) READ_ERUP_VOLC    : Reads SO2 emissions from eruptive volcanoes 
!  (29) READ_ANTHRO_SOx   : Reads anthropogenic SO2 and SO4 emissions
!  (30) READ_OCEAN_DMS    : Reads biogenic DMS emissions from oceans
!  (31) READ_SST          : Reads monthly mean sea-surface temperatures
!  (32) READ_BIOFUEL_SO2  : Reads SO2 emissions from biomass burning
!  (33) READ_AIRCRAFT_SO2 : Reads SO2 emissions from aircraft exhaust
!  (34) READ_SHIP_SO2     : Reads SO2 emissions from ship exhaust
!  (35) READ_ANTHRO_NH3   : Reads NH3 emissions from anthropogenic sources
!  (36) READ_NATURAL_NH3  : Reads NH3 emissions from natural sources
!  (37) READ_BIOMASS_NH3  : Reads NH3 biomass burning emissions
!  (38) READ_OXIDANT      : Reads monthly mean O3 and H2O2 for offline run
!  (39) OHNO3TIME         : Computes time arrays for scaling offline OH, NO3
!  (40) INIT_SULFATE      : Allocates & zeroes module arrays
!  (41) CLEANUP_SULFATE   : Deallocates module arrays
!
!  GEOS-CHEM modules referenced by sulfate_mod.f
!  ============================================================================
!  (1 ) biomass_mod.f          : Module w/ routines for biomass burning
!  (2 ) bpch2_mod.f            : Module w/ routines for binary pch file I/O
!  (3 ) bravo_mod.f            : Module w/ routines to read BRAVO emissions
!  (4 ) comode_mod.f           : Module w/ SMVGEAR allocatable arrays
!  (5 ) dao_mod.f              : Module w/ DAO met field arrays
!  (6 ) diag_mod.f             : Module w/ GEOS-Chem diagnostic arrays
!  (7 ) directory_mod.f        : Module w/ GEOS-Chem data & met field dirs
!  (8 ) drydep_mod.f           : Module w/ GEOS-Chem dry deposition routines
!  (9 ) epa_nei_mod.f          : Module w/ routines to read EPA/NEI99 data
!  (10) error_mod.f            : Module w/ NaN, other error check routines
!  (11) file_mod.f             : Module w/ file unit numbers & error checks
!  (12) future_emissions_mod.f : Module w/ routines for IPCC future emissions
!  (13) grid_mod.f             : Module w/ horizontal grid information
!  (14) global_no3_mod.f       : Module w/ routines to read 3-D NO3 field
!  (15) global_oh_mod.f        : Module w/ routines to read 3-D OH field
!  (16) isoropia_mod.f         : Module w/ ISORROPIA routines for aer thermo eq
!  (17) logical_mod.f          : Module w/ GEOS-Chem logical switches
!  (18) pbl_mix_mod.f          : Module w/ routines for PBL height & mixing
!  (19) pressure_mod.f         : Module w/ routines to compute P(I,J,L)
!  (20) seasalt_mod.f          : Module w/ routines for seasalt chemistry
!  (21) streets_anthro_mod.f   : Module w/ routines for David Streets' emiss
!  (22) time_mod.f             : Module w/ routines to compute time & date
!  (23) tracer_mod.f           : Module w/ GEOS-Chem tracer array STT etc.
!  (24) tracerid_mod.f         : Module w/ pointers to tracers & emissions
!  (25) transfer_mod.f         : Module w/ routines to cast & resize arrays
!  (26) tropopause_mod.f       : Module w/ routines to read ann mean tropopause
!  (27) uvalbedo_mod.f         : Module w/ UV albedo array and reader
!  (28) wetscav_mod.f          : Module w/ routines for wetdep & scavenging
!
!  References:
!  ============================================================================
!  (1 ) Andreae, M.O. & P. Merlet, "Emission of trace gases and aerosols from
!        biomass burning", Global Biogeochem. Cycles, 15, 955-966, 2001.
!  (2 ) Nightingale et al [2000a], J. Geophys. Res, 14, 373-387
!  (3 ) Nightingale et al [2000b], Geophys. Res. Lett, 27, 2117-2120
!  (4 ) Wanninkhof, R., "Relation between wind speed and gas exchange over
!        the ocean", J. Geophys. Res, 97, 7373-7382, 1992.
!
!  NOTES:
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
!  (50) Standardized patch in READ_ANTHRO_NH3 (dkh, bmy, 3/5/10)
!  (51) Use LWC from GEOS-5 met fields (jaf, bmy, 6/30/10)
!  (52) Add module parameters MNYEAR_VOLC and MXYEAR_VOLC to define the 1st 
!       and last year with data for volcanic emissions. (ccc, 9/30/10)
!  (53) Use updated volcanic emissions from 1979 to 2009
!  26 Aug 2010 - R. Yantosca - Add modifications for MERRA
!  12 Nov 2010 - R. Yantosca - Avoid div-by-zero when computing L2S, L3S
!******************************************************************************
!
      USE LOGICAL_MOD,   ONLY : LNLPBL ! (Lin, 03/31/09)
      USE VDIFF_PRE_MOD, ONLY : emis_save ! (Lin, 03/31/09)

      IMPLICIT NONE
      
      !=================================================================
      ! MODULE PRIVATE DECLARATIONS -- keep certain internal variables 
      ! and routines from being seen outside "sulfate_mod.f"
      !=================================================================

      ! Make everything PRIVATE ...
      PRIVATE

      ! ... except these routines
      PUBLIC :: CHEMSULFATE       
      PUBLIC :: EMISSSULFATE      
      PUBLIC :: CLEANUP_SULFATE   
     
      !=================================================================
      ! MODULE VARIABLES (see descriptions listed above)
      !=================================================================

      ! Time variable
      INTEGER              :: ELAPSED_SEC

      ! Parameters
      REAL*8,  PARAMETER   :: XNUMOL_OH  = 6.022d23 / 17d-3
      REAL*8,  PARAMETER   :: XNUMOL_O3  = 6.022d23 / 48d-3
      REAL*8,  PARAMETER   :: XNUMOL_NO3 = 6.022d23 / 62d-3
      REAL*8,  PARAMETER   :: TCVV_S     = 28.97d0  / 32d0
      REAL*8,  PARAMETER   :: SMALLNUM   = 1d-20

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

      ! Volcanic emissions parameters (ccc, 9/30/10)
      INTEGER, PARAMETER   :: LVOLC=20            ! Volcano levels
      INTEGER, PARAMETER   :: MNYEAR_VOLC = 1979  ! Min year of data
      INTEGER, PARAMETER   :: MXYEAR_VOLC = 2009  ! Max year of data


      ! Eruptive volcanoes
      REAL*8,  ALLOCATABLE :: EEV(:,:,:)

      ! Non-eruptive volcanoes 
      REAL*8,  ALLOCATABLE :: ENV(:,:,:)
      
      ! Pointers to drydep species w/in DEPSAV
      INTEGER              :: DRYSO2,  DRYSO4,   DRYMSA,  DRYNH3  
      INTEGER              :: DRYNH4,  DRYNIT,   DRYSO4s, DRYNITs
      INTEGER              :: DRYH2O2, DRYSO4aq, DRYAS,   DRYAHS
      INTEGER              :: DRYLET,  DRYNH4aq

      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement
      !=================================================================
      CONTAINS

!------------------------------------------------------------------------------

      SUBROUTINE GET_VCLDF
!
!******************************************************************************
!  Subroutine GET_VCLDF computes the volume cloud fraction for SO2 chemistry.
!  (rjp, bdf, bmy, 9/23/02)
!
!  References:
!  ============================================================================
!  (1) Sundqvist et al. [1989]
!
!  NOTES:
!  14 Jan 2011 - R. Yantosca - Return if VCLDF is not allocated
!******************************************************************************
!
      ! References to F90 modules 
      USE DAO_MOD,      ONLY : RH
      USE PRESSURE_MOD, ONLY : GET_PCENTER, GET_PEDGE

#     include "F77_CMN_SIZE"   ! Size parameters

      ! Local variables
      INTEGER              :: I,    J,    L
      REAL*8               :: PRES, PSFC, RH2, R0, B0

      ! Parameters
      REAL*8,  PARAMETER   :: ZRT = 0.60d0, ZRS = 0.99d0
		
      !=================================================================
      ! GET_VCLDF begins here!
      !=================================================================

      ! Exit if we VCLDF is not allocated.  We will now get the cloud
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
         RH2  = RH(I,J,L) * 1.0d-2

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

      ! Return to calling program
      END SUBROUTINE GET_VCLDF

!------------------------------------------------------------------------------

      FUNCTION GET_LWC( T ) RESULT( LWC )
!
!******************************************************************************
!  Function GET_LWC returns the cloud liquid water content [m3 H2O/m3 air] at 
!  a  GEOS-CHEM grid box as a function of temperature. (rjp, bmy, 10/31/02, 
!  1/14/03)
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

      ! Return to calling program
      END FUNCTION GET_LWC

!------------------------------------------------------------------------------
      
      SUBROUTINE CHEMSULFATE
!
!******************************************************************************
!  Subroutine CHEMSULFATE is the interface between the GEOS-CHEM main program
!  and the sulfate chemistry routines.  The user has the option of running
!  a coupled chemistry-aerosols simulation or an offline aerosol simulation.
!  (rjp, bdf, bmy, 5/31/00, 3/16/06)
!
!  NOTES:
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
!******************************************************************************
!
      ! References to F90 modules
      USE DAO_MOD,        ONLY : AD,     AIRDEN,  CLDF
      USE DAO_MOD,        ONLY : SUNCOS, CONVERT_UNITS
      USE DRYDEP_MOD,     ONLY : DEPSAV
      USE ERROR_MOD,      ONLY : DEBUG_MSG
      USE GLOBAL_OH_MOD,  ONLY : GET_GLOBAL_OH
      USE GLOBAL_NO3_MOD, ONLY : GET_GLOBAL_NO3
      USE LOGICAL_MOD,    ONLY : LCRYST,          LPRT
      USE TIME_MOD,       ONLY : GET_MONTH,       GET_TS_CHEM
      USE TIME_MOD,       ONLY : GET_ELAPSED_SEC, ITS_A_NEW_MONTH
      USE TRACER_MOD,     ONLY : STT,             TCVV 
      USE TRACER_MOD,     ONLY : N_TRACERS,       ITS_AN_AEROSOL_SIM
      USE TRACERID_MOD,   ONLY : IDTNITs,         IDTSO4s
 
#     include "F77_CMN_SIZE"     ! Size parameters 

      ! Local variables
      LOGICAL, SAVE       :: FIRSTCHEM = .TRUE.
      INTEGER, SAVE       :: LASTMONTH = -99
      INTEGER             :: I, J, L, N, MONTH
      REAL*8              :: DTCHEM

      ! External functions   
      REAL*8,  EXTERNAL   :: BOXVL

      !=================================================================
      ! CHEMSULFATE begins here!
      !=================================================================

      ! Get current month
      MONTH = GET_MONTH()

      ! Establish indices w/in DEPSAV array
      IF ( FIRSTCHEM ) THEN

         ! Initialize arrays (if not already done before)
         CALL INIT_SULFATE

         ! Reset first-time flag
         FIRSTCHEM = .FALSE.
      ENDIF

      ! If it's an offline simulation ...
      IF ( ITS_AN_AEROSOL_SIM() ) THEN

         ! Then read monthly data files ...
         IF ( ITS_A_NEW_MONTH() ) THEN 
            CALL GET_GLOBAL_OH( MONTH )
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
                  
      !================================================================= 
      ! Call individual chemistry routines for sulfate/aerosol tracers
      !=================================================================

      ! SO4s [kg] gravitational settling 
      CALL GRAV_SETTLING( STT(:,:,:,IDTSO4s), 1 )
      IF ( LPRT ) CALL DEBUG_MSG( '### CHEMSULFATE: GRAV_SET, SO4S' )

      ! NITs [kg] gravitational settling 
      CALL GRAV_SETTLING( STT(:,:,:,IDTNITs), 2 )
      IF ( LPRT ) CALL DEBUG_MSG( '### CHEMSULFATE: GRAV_SET, NITS' )

      ! Convert all tracers in STT from [kg] -> [v/v] 
      CALL CONVERT_UNITS( 1, N_TRACERS, TCVV, AD, STT )
      IF ( LPRT ) CALL DEBUG_MSG( '### CHEMSULFATE: a CONVERT UNITS' )

      ! For offline runs only ...
      IF ( ITS_AN_AEROSOL_SIM() ) THEN
         
         ! DMS (offline only)
         CALL CHEM_DMS
         IF ( LPRT ) CALL DEBUG_MSG( '### CHEMSULFATE: a CHEM_DMS' ) 
      
         ! H2O2 (offline only)
         CALL CHEM_H2O2
         IF ( LPRT ) CALL DEBUG_MSG( '### CHEMSULFATE: a CHEM_H2O2' )

      ENDIF

      ! SO2 
      CALL GET_VCLDF
      IF ( LPRT ) CALL DEBUG_MSG( '### CHEMSULFATE: a get VCLDF' )
      CALL CHEM_SO2
      IF ( LPRT ) CALL DEBUG_MSG( '### CHEMSULFATE: a CHEM_SO2' )

      ! SO4 
      CALL CHEM_SO4
      IF ( LPRT ) CALL DEBUG_MSG( '### CHEMSULFATE: a CHEM_SO4' )

!-----------------------------------------------------------------------------
!%%% Currently under development (rjp, bmy, 3/15/05)
!%%%      ! Only do the following if the crystalline sulfate & aqueous 
!%%%      ! tracers (AS, AHS, LET, SO4aq, NH4aq) are defined
!%%%      IF ( LCRYST ) THEN
!%%%         
!%%%         ! Phase change
!%%%         CALL PHASE_SO4
!%%%         IF ( LPRT ) CALL DEBUG_MSG( '### CHEMSULFATE: a PHASE_SO4' )
!%%%      
!%%%         ! Radiative forcing
!%%%         CALL PHASE_RADIATIVE
!%%%         IF ( LPRT ) CALL DEBUG_MSG( '### CHEMSULFATE: a PHASE_RAD' )
!%%%      
!%%%      ENDIF
!-----------------------------------------------------------------------------

      ! MSA 
      CALL CHEM_MSA
      IF ( LPRT ) CALL DEBUG_MSG( '### CHEMSULFATE: a CHEM_MSA' )

      ! NH3 
      CALL CHEM_NH3
      IF ( LPRT ) CALL DEBUG_MSG( '### CHEMSULFATE: a CHEM_NH3' )

      ! NH4 (gas-phase)
      CALL CHEM_NH4
      IF ( LPRT ) CALL DEBUG_MSG( '### CHEMSULFATE: a CHEM_NH4' )

!-----------------------------------------------------------------------------
!%%% Currently under development (rjp, bmy, 3/15/05)
!%%%      ! NH4 (aqueous phase)
!%%%      IF ( LCRYST ) THEN
!%%%         CALL CHEM_NH4aq
!%%%         IF ( LPRT ) CALL DEBUG_MSG( '### CHEMSULFATE: a CHEM_NH4aq' )
!%%%      ENDIF
!-----------------------------------------------------------------------------

      ! Sulfur Nitrate 
      CALL CHEM_NIT
      IF ( LPRT ) CALL DEBUG_MSG( '### CHEMSULFATE: a CHEM_NIT' )
 
      ! Convert STT from [v/v] -> [kg]
      CALL CONVERT_UNITS( 2, N_TRACERS, TCVV, AD, STT )

      ! We have already gone thru one chemistry iteration
      FIRSTCHEM = .FALSE. 
         
      ! Return to calling program
      END SUBROUTINE CHEMSULFATE

!------------------------------------------------------------------------------

      SUBROUTINE GRAV_SETTLING( TC, N )
!
!******************************************************************************
!  Subroutine GRAV_SETTLING performs gravitational settling of sulfate
!  and nitrate in coarse sea salt (SO4S and NITS).
!  (bec, rjp, bmy, 4/20/04, 7/20/04, 10/25/05)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) TC (REAL*8 ) : Tracer [kg]
!  (2 ) N  (INTEGER) : N=1 is SO4S; N=2 is NITS
!
!  Arguments as Output:
!  ============================================================================
!  (1 ) TC (REAL*8 ) : Contains modified tracer
!
!  NOTES:
!  (1 ) Now references SALA_REDGE_um and SALC_REDGE_um from "tracer_mod.f"
!        (bmy, 7/20/04)
!  (2 ) Now references XNUMOL from "tracer_mod.f" (bmy, 10/25/05)
!  (3 ) Now limit relative humidity to [tiny(real*8),0.99] range for DLOG
!         argument (phs, 5/1/08)
!******************************************************************************
!
      ! References to F90 modules
      USE DAO_MOD,       ONLY : T, BXHEIGHT, RH
      USE DIAG_MOD,      ONLY : AD44
      USE DRYDEP_MOD,    ONLY : DEPSAV
      USE PRESSURE_MOD,  ONLY : GET_PCENTER
      USE TRACER_MOD,    ONLY : SALA_REDGE_um,   SALC_REDGE_um,  XNUMOL
      USE TRACERID_MOD,  ONLY : IDTSO4s,         IDTNITs
      USE TIME_MOD,      ONLY : GET_ELAPSED_SEC, GET_TS_CHEM
      USE GRID_MOD,      ONLY : GET_AREA_CM2

#     include "F77_CMN_SIZE"      ! Size parameters
#     include "F77_CMN_GCTM"      ! g0
#     include "F77_CMN_DIAG"      ! ND44

      ! Arguments
      INTEGER, INTENT(IN)    :: N
      REAL*8,  INTENT(INOUT) :: TC(IIPAR,JJPAR,LLPAR)

      ! Local variables
      INTEGER                :: I,      J,     L,        DTCHEM
      REAL*8                 :: DELZ,   DELZ1, REFF
      REAL*8                 :: P,      DP,    PDP,      TEMP        
      REAL*8                 :: CONST,  SLIP,  VISC,     FAC1
      REAL*8                 :: FAC2,   FLUX,  AREA_CM2, RHB
      REAL*8                 :: RCM,    RWET,  RATIO_R,  RHO
      REAL*8                 :: TOT1,   TOT2
      REAL*8                 :: VTS(LLPAR)  
      REAL*8                 :: TC0(LLPAR)
      
      ! Parameters
      REAL*8,  PARAMETER     :: C1 =  0.7674d0 
      REAL*8,  PARAMETER     :: C2 =  3.079d0 
      REAL*8,  PARAMETER     :: C3 =  2.573d-11
      REAL*8,  PARAMETER     :: C4 = -1.424d0
      REAL*8,  PARAMETER     :: DEN = 2200.0d0 ! [kg/m3] sea-salt density

      ! Arrays
      INTEGER              :: IDDEP(2)
      INTEGER              :: IDTRC(2)	

      !=================================================================
      ! GRAV_SETTLING begins here!
      !=================================================================

      ! Return if tracers are undefined
      IF ( IDTSO4s == 0 .and. IDTNITs == 0 ) RETURN

      ! Return if it's the start of the run
      IF ( GET_ELAPSED_SEC() == 0 ) RETURN

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
      RCM  = REFF * 100d0  

      ! Exponential factors
      FAC1 = C1 * ( RCM**C2 )
      FAC2 = C3 * ( RCM**C4 )

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I,       J,     L,    VTS,  P,        TEMP, RHB,  RWET ) 
!$OMP+PRIVATE( RATIO_R, RHO,   DP,   PDP,  CONST,    SLIP, VISC, TC0  )
!$OMP+PRIVATE( DELZ,    DELZ1, TOT1, TOT2, AREA_CM2, FLUX             )
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
            TEMP    = T(I,J,L)

            ! Cap RH at 0.99 
            RHB     = MIN( 0.99d0, RH(I,J,L) * 1d-2 )

            ! Safety check (phs, 5/1/08)
            RHB     = MAX( TINY(RHB), RHB           )

            ! Aerosol growth with relative humidity in radius [m] 
            ! (Gerber, 1985)
            RWET    = 0.01d0*(FAC1/(FAC2-DLOG(RHB))+RCM**3.d0)**0.33d0

            ! Ratio dry over wet radii at the cubic power
            RATIO_R = ( REFF / RWET )**3.d0

            ! Density of the wet aerosol (kg/m3)
            RHO     = RATIO_R * DEN + ( 1.d0 - RATIO_R ) * 1000.d0

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
         DELZ = BXHEIGHT(I,J,L)

         TC(I,J,L) = TC(I,J,L) / ( 1.d0 + DTCHEM * VTS(L) / DELZ )

         DO L = LLTROP-1, 1, -1
            DELZ  = BXHEIGHT(I,J,L)
            DELZ1 = BXHEIGHT(I,J,L+1)
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
            AREA_CM2 = GET_AREA_CM2( J )

            ! Convert sea salt flux from [kg/s] to [molec/cm2/s]
            FLUX     = ( TOT1 - TOT2 ) / DTCHEM
            FLUX     = FLUX * XNUMOL(IDTRC(N)) / AREA_CM2 
   
            ! Store in AD44 array
            AD44(I,J,IDDEP(N),1) = AD44(I,J,IDDEP(N),1) + FLUX
         ENDIF
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! Return to calling program
      END SUBROUTINE GRAV_SETTLING

!------------------------------------------------------------------------------

      SUBROUTINE CHEM_DMS
!
!******************************************************************************
!  Subroutine CHEM_DMS is the DMS chemistry subroutine from Mian Chin's    
!  GOCART model, modified for use with the GEOS-CHEM model.
!  (rjp, bdf, bmy, 5/31/00, 10/15/09)  
!                                                                           
!  Module Variables used:                                                     
!  ============================================================================
!  (1 ) PSO2_DMS (REAL*8 ) : Array for P(SO2) from DMS [v/v]                
!  (2 ) PMSA_DMS (REAL*8 ) : Array for P(MSA) from DMS [v/v]                
!                                                                            
!  Reaction List (by Mian Chin, chin@rondo.gsfc.nasa.gov)                  
!  ============================================================================
!                                                                           
!  R1:    DMS + OH  -> a*SO2 + b*MSA                OH addition channel    
!         k1 = { 1.7e-42*exp(7810/T)*[O2] / (1+5.5e-31*exp(7460/T)*[O2] }  
!         a = 0.75, b = 0.25                                               
!                                                                           
!  R2:    DMS + OH  ->   SO2 + ...                  OH abstraction channel 
!         k2 = 1.2e-11*exp(-260/T)                                         
!                                                                           
!         DMS_OH = DMS0 * exp(-(r1+r2)* NDT1)                                  
!         where DMS0 is the DMS concentration at the beginning,            
!         r1 = k1*[OH], r2 = k2*[OH].                                      
!                                                                           
!  R3:    DMS + NO3 ->   SO2 + ...                                         
!         k3 = 1.9e-13*exp(500/T)                                          
!                                                                           
!         DMS = DMS_OH * exp(-r3*NDT1)                                         
!         where r3 = k3*[NO3].                                             
!                                                                           
!  R4:    DMS + X   ->   SO2 + ...                                         
!         assume to be at the rate of DMS+OH and DMS+NO3 combined.         
!                                                                           
!  The production of SO2 and MSA here, PSO2_DMS and PMSA_DMS, are saved    
!  for use in CHEM_SO2 and CHEM_MSA subroutines as a source term.  They    
!  are in unit of [v/v/timestep]. 
!
!  NOTES: 
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
!******************************************************************************
!
      ! Reference to F90 modules
      USE DAO_MOD,        ONLY : AD, AIRDEN, SUNCOS, T
      USE DIAG_MOD,       ONLY : AD05
      USE DRYDEP_MOD,     ONLY : DEPSAV
      USE TIME_MOD,       ONLY : GET_TS_CHEM
      USE TRACER_MOD,     ONLY : STT, ITS_A_FULLCHEM_SIM, XNUMOL
      USE TRACERID_MOD,   ONLY : IDTDMS
      USE TROPOPAUSE_MOD, ONLY : ITS_IN_THE_STRAT

#     include "F77_CMN_SIZE"     ! Size parameters
#     include "F77_CMN_GCTM"     ! AIRMW
#     include "F77_CMN_DIAG"     ! ND05, LD05

      ! Local variables
      LOGICAL                :: IS_FULLCHEM
      INTEGER                :: I,   J,    L,      IJLOOP
      REAL*8                 :: TK,  O2,   RK1,    RK2,    RK3,   F  
      REAL*8                 :: DMS, DMS0, DMS_OH, DTCHEM, XOH,   XN3 
      REAL*8                 :: XX,  OH,   OH0,    XNO3,   XNO30, LOH
      REAL*8                 :: LNO3

      ! Parameters
      REAL*8, PARAMETER      :: FX = 1.0d0
      REAL*8, PARAMETER      :: A  = 0.75d0
      REAL*8, PARAMETER      :: B  = 0.25d0

      ! From D4: only 0.8 efficiency, also some goes to DMSO and lost.  
      ! So we assume 0.75 efficiency for DMS addtion channel to form     
      ! products.                                                        
      REAL*8, PARAMETER      :: EFF = 1d0
      
      ! External functions   
      REAL*8,  EXTERNAL      :: BOXVL
      
      !=================================================================
      ! CHEM_DMS begins here!
      !=================================================================
      IF ( IDTDMS == 0 ) RETURN

      ! Flag for fullchem simulation
      IS_FULLCHEM = ITS_A_FULLCHEM_SIM()

      ! DTCHEM is the chemistry timestep in seconds
      DTCHEM = GET_TS_CHEM() * 60d0

      ! Factor to convert AIRDEN from kgair/m3 to molecules/cm3:
      f  = 1000.d0 / AIRMW * 6.022d23 * 1.d-6
      
      !=================================================================
      ! Do the chemistry over all tropospheric grid boxes!
      !=================================================================
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, IJLOOP, J, L, TK, O2, DMS0, OH, XNO3, RK1, RK2 )
!$OMP+PRIVATE( RK3, DMS_OH, DMS, OH0, XNO30, XOH, XN3, XX, LOH, LNO3  )
!$OMP+SCHEDULE( DYNAMIC )
      DO L = 1, LLTROP  
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Skip stratospheric boxes
         IF ( ITS_IN_THE_STRAT( I, J, L ) ) CYCLE

         ! IJLOOP is the 1-D grid box index for SUNCOS
         IJLOOP = ( (J-1) * IIPAR ) + I

         ! Temperature [K]
         TK     = T(I,J,L)

         ! Get O2 [molec/cm3], DMS [v/v], OH [molec/cm3], NO3 [molec/cm3]
         O2     = AIRDEN(L,I,J) * f * 0.21d0
         DMS0   = STT(I,J,L,IDTDMS)
         OH     = GET_OH( I, J, L )
         XNO3   = GET_NO3( I, J, L )

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
         IF ( SUNCOS(IJLOOP) <= 0d0 ) THEN
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
            OH    = OH0 - ( ( DMS0 - DMS_OH ) * AIRDEN(L,I,J) * f )
            IF ( OH < SMALLNUM ) OH = 0d0

            ! Update NO3 after rxn w/ DMS (coupled runs only)
            XNO3  = XNO30 - ( ( DMS_OH - DMS ) * AIRDEN(L,I,J) * f )
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
     $           Fx * AD(I,J,L) / TCVV_S
            XN3  = ( DMS_OH - DMS    ) / Fx * AD(I,J,L) / TCVV_S
            XX   = ( ( DMS0 - DMS ) * AD(I,J,L) / TCVV_S ) - XOH - XN3
        
            ! Convert L(OH) and L(NO3) from [molec/cm3] to [kg/timestep]
            LOH  = ( OH0   - OH   ) * BOXVL(I,J,L) / XNUMOL_OH
            LNO3 = ( XNO30 - XNO3 ) * BOXVL(I,J,L) / XNUMOL_NO3 

            ! Store P(SO2) from DMS + OH [kg S/timestep]
            AD05(I,J,L,1) = AD05(I,J,L,1) + XOH

            ! Store P(SO2) from DMS + NO3 [kg S/timestep]
            AD05(I,J,L,2) = AD05(I,J,L,2) + XN3

            ! Store total P(SO2) from DMS [kg S/timestep]
            AD05(I,J,L,3) = AD05(I,J,L,3) + 
     &                      ( PSO2_DMS(I,J,L) * AD(I,J,L) / TCVV_S )

            ! Store P(MSA) from DMS [kg S/timestep]
            AD05(I,J,L,4) = AD05(I,J,L,4) + 
     &                      ( PMSA_DMS(I,J,L) * AD(I,J,L) / TCVV_S )

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

      ! Return to calling program
      END SUBROUTINE CHEM_DMS

!------------------------------------------------------------------------------

      SUBROUTINE CHEM_H2O2
!
!******************************************************************************
!  Subroutine CHEM_H2O2 is the H2O2 chemistry subroutine for offline sulfate
!  simulations.  For coupled runs, H2O2 chemistry is already computed by
!  the SMVGEAR module. (rjp, bmy, 11/26/02, 10/25/05)
!                                                                           
!  NOTES:
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
!******************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD,      ONLY : GET_NAME_EXT, GET_RES_EXT
      USE BPCH2_MOD,      ONLY : GET_TAU0,     READ_BPCH2
      USE DAO_MOD,        ONLY : AD, AIRDEN, OPTD, SUNCOS, T
      USE DIAG_MOD,       ONLY : AD44 
      USE DIRECTORY_MOD,  ONLY : DATA_DIR
      USE DRYDEP_MOD,     ONLY : DEPSAV
      USE GRID_MOD,       ONLY : GET_AREA_CM2
      USE PBL_MIX_MOD,    ONLY : GET_FRAC_UNDER_PBLTOP
      USE TIME_MOD,       ONLY : GET_MONTH, GET_TS_CHEM, ITS_A_NEW_MONTH
      USE TRACER_MOD,     ONLY : STT,       TCVV,        XNUMOL
      USE TRACERID_MOD,   ONLY : IDTH2O2
      USE TRANSFER_MOD,   ONLY : TRANSFER_3D_TROP
      USE UVALBEDO_MOD,   ONLY : UVALBEDO
      USE TROPOPAUSE_MOD, ONLY : ITS_IN_THE_STRAT
 
#     include "cmn_fj.h"       ! IPAR, JPAR, LPAR, F77_CMN_SIZE
#     include "F77_CMN_GCTM"       ! AIRMW
#     include "F77_CMN_DIAG"       ! ND44
      
      ! Local variables
      LOGICAL                 :: FIRST     = .TRUE.
      INTEGER, SAVE           :: LASTMONTH = -99
      INTEGER                 :: I, J, L, JLOOP
      REAL*4                  :: ARRAY(IIPAR,JJPAR,LLTROP)
      REAL*8                  :: ND44_TMP(IIPAR,JJPAR,LLTROP)
      REAL*8                  :: DT,    Koh,  DH2O2, M,    F ,   XTAU   
      REAL*8                  :: H2O20, H2O2, ALPHA, FLUX, FREQ, PHOTJ
      REAL*8,  PARAMETER      :: A = 2.9d-12
      CHARACTER(LEN=255)      :: FILENAME

      !=================================================================
      ! CHEM_H2O2 begins here!
      !=================================================================
      IF ( IDTH2O2 == 0 .or. DRYH2O2 == 0 ) RETURN 

      ! Chemistry timestep [s]
      DT = GET_TS_CHEM() * 60d0

      ! Factor to convert AIRDEN from kgair/m3 to molecules/cm3:
      F  = 1000.d0 / AIRMW * 6.022d23 * 1.d-6
      
      ! Zero ND44_TMP array
      IF ( ND44 > 0 ) THEN
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
         DO L = 1, LLTROP
         DO J = 1, JJPAR
         DO I = 1, IIPAR
            ND44_TMP(I,J,L) = 0d0
         ENDDO
         ENDDO
         ENDDO
!$OMP END PARALLEL DO
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
         WRITE( 6, 100 ) TRIM( FILENAME )
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
!     &                    XTAU,      IIPAR,     JJPAR,      
!     &                    LLTROP,    ARRAY,     QUIET=.TRUE. )



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
      DO L  = 1, LLTROP
      DO J  = 1, JJPAR
      DO I  = 1, IIPAR

         ! Initialize for safety's sake 
         FLUX = 0d0
         FREQ = 0d0

         ! Skip stratospheric boxes
         IF ( ITS_IN_THE_STRAT( I, J, L ) ) CYCLE

         ! Density of air [molec/cm3]
         M     = AIRDEN(L,I,J) * f  

         ! Initial H2O2 [v/v]
         H2O20 = STT(I,J,L,IDTH2O2)

         ! Loss frequenty due to OH oxidation [s-1]
         KOH   = A * EXP( -160.d0 / T(I,J,L) ) * GET_OH(I,J,L)  

         ! H2O2 drydep frequency [1/s].  Account for the fraction
         ! of grid box (I,J,L) that is located beneath the PBL top.
         FREQ  = DEPSAV(I,J,DRYH2O2) * GET_FRAC_UNDER_PBLTOP( I, J, L ) 

         ! Add option for non-local PBL (Lin, 03/31/09)
         IF (LNLPBL) FREQ = 0.D0

         ! 1-D grid box index for SUNCOS
         JLOOP = ( (J-1) * IIPAR ) + I

         ! Impose a diurnal variation of jH2O2 by multiplying COS of 
         ! solar zenith angle normalized by maximum solar zenith angle 
         ! because the archived JH2O2 is for local noon time
         IF ( COSZM(I,J) > 0.d0 ) THEN
            PHOTJ = JH2O2(I,J,L) * SUNCOS(JLOOP) / COSZM(I,J)
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
            FLUX = FLUX * AD(I,J,L) / TCVV(IDTH2O2)
            FLUX = FLUX * XNUMOL(IDTH2O2) / ( GET_AREA_CM2( J ) * DT )

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
         DO L = 1, LLTROP
            AD44(I,J,DRYH2O2,1) = AD44(I,J,DRYH2O2,1) + ND44_TMP(I,J,L)
         ENDDO
         ENDDO
         ENDDO
!$OMP END PARALLEL DO
      ENDIF

      ! Return to calling program
      END SUBROUTINE CHEM_H2O2

!------------------------------------------------------------------------------

      SUBROUTINE CHEM_SO2
!
!******************************************************************************
!  Subroutine CHEM_SO2 is the SO2 chemistry subroutine 
!  (rjp, bmy, 11/26/02, 8/26/10) 
!                                                                          
!  Module variables used:
!  ============================================================================
!  (1 ) PSO2_DMS (REAL*8 ) : Array for P(SO2) from DMS          [v/v/timestep]
!  (2 ) PSO4_SO2 (REAL*8 ) : Array for P(SO4) from SO2          [v/v/timestep]
!  (3 ) LSO2_AQ  (REAL*8 ) : Array for L(SO2) from Aqueuos chem [v/v/timestep]
!                                                                           
!  Reaction List (by Rokjin Park, rjp@io.harvard.edu)                      
!  ============================================================================
!  (1 ) SO2 production:                                                      
!       DMS + OH, DMS + NO3 (saved in CHEM_DMS)                               
!                                                                          
!  (2 ) SO2 loss:                                                         
!       (a) SO2 + OH  -> SO4                                               
!       (b) SO2       -> drydep                                             
!       (c) SO2 + H2O2 or O3 (aq) -> SO4                         
!                                                                          
!  (3 ) SO2 = SO2_0 * exp(-bt) +  PSO2_DMS/bt * [1-exp(-bt)]   
! 
!       where b is the sum of the reaction rate of SO2 + OH and the dry       
!       deposition rate of SO2, PSO2_DMS is SO2 production from DMS in        
!       MixingRatio/timestep.                                                 
!                                                                          
!  If there is cloud in the gridbox (fraction = fc), then the aqueous      
!  phase chemistry also takes place in cloud. The amount of SO2 oxidized   
!  by H2O2 in cloud is limited by the available H2O2; the rest may be      
!  oxidized due to additional chemistry, e.g, reaction with O3 or O2       
!  (catalyzed by trace metal).                                             
!                                                                          
!  NOTES:                                                                   
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
!  (14) Added extra error checks to prevent negative L2S, L3S (bmy, 4/28/10)
!  (15) Use liq. water content from met fields in GEOS-5 (jaf, bmy, 6/30/10)
!  26 Aug 2010 - R. Yantosca - Use liquid water content from MERRA
!  12 Nov 2010 - R. Yantosca - Prevent div-by-zero when computing L2S and L3S
!  27 May 2011 - L. Zhang    - Divide LWC by cloud fraction for GEOS/MERRA
!                              and adjust the L2S and L3S rates accordingly
!******************************************************************************
!
      ! Reference to diagnostic arrays
      USE DAO_MOD,         ONLY : AD,      AIRDEN, T
      USE DIAG_MOD,        ONLY : AD05,    AD44
      USE DRYDEP_MOD,      ONLY : DEPSAV
      USE DIRECTORY_MOD,   ONLY : DATA_DIR
      USE ERROR_MOD,       ONLY : IS_SAFE_EXP
      USE ERROR_MOD,       ONLY : SAFE_DIV
      USE GLOBAL_HNO3_MOD, ONLY : GET_GLOBAL_HNO3
      USE GRID_MOD,        ONLY : GET_AREA_CM2
      USE PBL_MIX_MOD,     ONLY : GET_FRAC_UNDER_PBLTOP
      USE PRESSURE_MOD,    ONLY : GET_PCENTER
      USE TIME_MOD,        ONLY : GET_TS_CHEM, GET_MONTH
      USE TIME_MOD,        ONLY : ITS_A_NEW_MONTH
      USE TRACER_MOD,      ONLY : STT,     TCVV,  ITS_AN_AEROSOL_SIM
      USE TRACER_MOD,      ONLY : XNUMOL
      USE TRACERID_MOD,    ONLY : IDTH2O2, IDTSO2
      USE SEASALT_MOD,     ONLY : GET_ALK
      USE WETSCAV_MOD,     ONLY : H2O2s,   SO2s
      USE TROPOPAUSE_MOD,  ONLY : ITS_IN_THE_STRAT

      ! For LWC from met fields in GEOS-5 (jaf, 6/30/10)
      USE DAO_MOD,         ONLY : AIRDEN, QL, CLDF 

#     include "F77_CMN_SIZE"    ! Size parameters
#     include "F77_CMN_GCTM"    ! AIRMW
#     include "F77_CMN_DIAG"    ! LD05, ND05, ND44

      ! Local variables
      LOGICAL               :: IS_OFFLINE
      INTEGER               :: I,      J,       L,      I1,   I2
      INTEGER               :: II,     NSTEP
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
      REAL*8                :: ND44_TMP(IIPAR,JJPAR,LLTROP)

      ! Parameters
      REAL*8,  PARAMETER    :: HPLUS  = 3.16227766016837953d-5  !pH = 4.5
      REAL*8,  PARAMETER    :: MINDAT = 1.d-20

      !=================================================================
      ! CHEM_SO2 begins here!
      !=================================================================
      IF ( IDTH2O2 == 0 .or. IDTSO2 == 0 .or. DRYSO2 == 0 ) RETURN

      ! Is it an offline simulation?
      IS_OFFLINE = ITS_AN_AEROSOL_SIM()

      ! Read HNO3 for offline simulation
      IF ( IS_OFFLINE ) THEN
         IF ( ITS_A_NEW_MONTH() ) THEN
            CALL GET_GLOBAL_HNO3( GET_MONTH() )
         ENDIF
      ENDIF

      ! DTCHEM is the chemistry timestep in seconds
      DTCHEM = GET_TS_CHEM() * 60d0

      ! Factor to convert AIRDEN from [kg air/m3] to [molec air/cm3]
      F      = 1000.d0 / AIRMW * 6.022d23 * 1.d-6

      ! Zero ND44_TMP array
      ND44_TMP = 0d0
      
      ! Loop over tropospheric grid boxes
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, SO20, H2O20, O3, PATM, TK, K0, M, KK, F1, RK1  )
!$OMP+PRIVATE( RK2, RK, RKT, SO2_cd, L1, Ld, L2, L2S, L3, L3S, FC, LWC )
!$OMP+PRIVATE( KaqH2O2, KaqO3, AREA_CM2, FLUX, ALK, ALK1, ALK2         )
!$OMP+PRIVATE( Kt1, Kt2, AREASS1, AREASS2, SO2_ss, Kt1N, Kt2N          )
!$OMP+PRIVATE( PSO4E, PSO4F, XX                                        )
!$OMP+SCHEDULE( DYNAMIC )
      DO L = 1, LLTROP  
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Initialize for safety's sake 
         AREA_CM2 = 0d0
         FLUX     = 0d0
         Ld       = 0d0

         ! Skip stratospheric boxes
         IF ( ITS_IN_THE_STRAT( I, J, L ) ) CYCLE

         ! Initial SO2, H2O2 and O3 [v/v]
         SO20   = STT(I,J,L,IDTSO2)         
         H2O20  = STT(I,J,L,IDTH2O2)
         O3     = GET_O3(I,J,L)

         ! PATM  : Atmospheric pressure in atm
         PATM   = GET_PCENTER( I, J, L ) / 1013.25d0

         ! TK : Temperature [K]
         TK     = T(I,J,L)

         IF ( IS_OFFLINE ) THEN

            ! Gas phase SO4 production is done here in offline run only 
            ! Updated to match JPL 2006 + full chem (jaf, 10/14/09)
            K0  = 3.3d-31 * ( 300.d0 / TK )**4.3d0
            Ki  = 1.6d-12
            M   = AIRDEN(L,I,J) * F
            KK  = K0 * M / Ki
            F1  = ( 1.d0 + ( LOG10( KK ) )**2 )**( -1 )
            RK1 = ( K0 * M / ( 1.d0 + KK ) ) * 0.6d0**F1 * GET_OH(I,J,L)

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
         IF (LNLPBL) RK2 = 0.D0

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
         CALL GET_ALK( I, J, L, ALK1, ALK2, Kt1, Kt2, Kt1N, Kt2N )

         ! Total alkalinity [kg]
         ALK = ALK1 + ALK2

         ! If (1) there is alkalinity, (2) there is SO2 present, and 
         ! (3) O3 is in excess, then compute seasalt SO2 chemistry
         IF  ( ( ALK    > MINDAT )  .AND.
     &         ( SO2_cd > MINDAT )  .AND. 
     &         ( SO2_cd < O3     ) ) THEN

            ! Compute oxidation of SO2 -> SO4 and condensation of
            ! HNO3 -> nitrate within the seasalt aerosol 
            CALL SEASALT_CHEM( I,      J,     L,   ALK1, ALK2,
     &                         SO2_cd, Kt1,   Kt2, Kt1N, Kt2N,
     &                         SO2_ss, PSO4E, PSO4F ) 

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
#if   defined ( GEOS_5 ) || defined( MERRA )

         !---------------------------------------------
         ! GEOS-5/MERRA: Get LWC, FC from met fields
         ! (jaf, bmy, 6/30/10)
         !---------------------------------------------

         ! Get cloud fraction from met fields
         FC      = CLDF(L,I,J)

         ! Get liquid water content [m3 H2O/m3 air] within cloud from met flds
         ! Units: [kg H2O/kg air] * [kg air/m3 air] * [m3 H2O/1e3 kg H2O]
         LWC     = QL(I,J,L) * AIRDEN(L,I,J) * 1D-3

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
         !---------------------------------------------------------------------
         ! Prior to 5/21/11:
         !LWC     = GET_LWC( TK ) * FC
         !---------------------------------------------------------------------
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
         IF ( ( FC     > 0.d0   )  .AND. 
     &        ( SO2_ss > MINDAT )  .AND. 
     &        ( TK     > 258.0  ) ) THEN

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

            ! We have used the in-cloud LWC to compute the sulfate
            ! aqueous chemistry.  We get the appropriate grid-box averaged 
            ! mass of SO2 and sulfate by multiplying the reaction rates
            ! L2S and L3s by the cloud fraction after the aqueous chemistry
            ! has been done.  (lzh, jaf, bmy, 5/27/11)
            L2S =  L2S * FC
            L3S =  L3S * FC    

            ! Make sure SO2_ss and H2O20 are in the proper range
            SO2_ss = MAX( SO2_ss - ( L2S + L3S ), MINDAT )
            H2O20  = MAX( H2O20  - L2S,           MINDAT )

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

         ENDIF

         ! Store updated SO2, H2O2 back to the tracer arrays 
         STT(I,J,L,IDTSO2)  = SO2s( I,J,L)
         STT(I,J,L,IDTH2O2) = H2O2s(I,J,L)

         ! SO2 chemical loss rate  = SO4 production rate [v/v/timestep]
         PSO4_SO2(I,J,L) = L1 + L2S + L3S + PSO4E
         PSO4_ss (I,J,L) = PSO4F

         !=================================================================
         ! ND05 Diagnostics [kg S/timestep]
         !=================================================================
         IF ( ND05 > 0 .and. L <= LD05 ) THEN
           
            ! P(SO4) from gas-phase oxidation [kg S/timestep]
            AD05(I,J,L,5) = AD05(I,J,L,5) +
     &                      ( L1  * AD(I,J,L) / TCVV_S )

            ! P(SO4) from aqueous-phase oxidation with H2O2 [kg S/timestep]
            AD05(I,J,L,6) = AD05(I,J,L,6) +
     &                      ( L2S * AD(I,J,L) / TCVV_S )

            ! P(SO4) from aqueous-phase oxidation with O3 [kg S/timestep]
            AD05(I,J,L,7) = AD05(I,J,L,7) +
     &                      ( L3S * AD(I,J,L) / TCVV_S )

            ! P(SO4) from O3 oxidation in sea-salt aerosols [kg S/timestep]
            AD05(I,J,L,8) = AD05(I,J,L,8) +
     &                      ( (PSO4E + PSO4F) * AD(I,J,L) / TCVV_S )

         ENDIF

         !=================================================================
         ! ND44 Diagnostic: Drydep flux of SO2 [molec/cm2/s]
         !=================================================================
         IF ( ND44 > 0 .AND. Ld > 0d0 ) THEN

            ! Surface area [cm2]
            AREA_CM2 = GET_AREA_CM2( J )

            ! Convert [v/v/timestep] to [molec/cm2/s]
            FLUX = Ld   * AD(I,J,L)      / TCVV(IDTSO2)
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
         DO L = 1, LLTROP
            AD44(I,J,DRYSO2,1) = AD44(I,J,DRYSO2,1) + ND44_TMP(I,J,L)
         ENDDO
         ENDDO
         ENDDO
!$OMP END PARALLEL DO
      ENDIF

      ! Return to calling program
      END SUBROUTINE CHEM_SO2

!------------------------------------------------------------------------------

      SUBROUTINE SEASALT_CHEM ( I,      J,     L,   ALK1, ALK2,
     &                          SO2_cd, Kt1,   Kt2, Kt1N, Kt2N,
     &                          SO2_ss, PSO4E, PSO4F )
!
!******************************************************************************
!  Function SEASALT_CHEM computes SO4 formed from S(IV) + O3 on seasalt 
!  aerosols as a function of seasalt alkalinity. (bec, bmy, 4/13/05, 10/7/08)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) I      (INTEGER) :
!  (2 ) J      (INTEGER) :
!  (3 ) L      (INTEGER) :
!  (4 ) O30    (REAL*8)  : Initial O3 mixing ratio (v/v]
!  (5 ) ALK    (REAL*8)  : Alkalinity [kg] from seasalt from seasalt_mod
!  (6 ) SO2_cd (REAL*8)  : SO2 mixing ratio [v/v] after gas phase chemistry
!                           and dry deposition
!  (7 ) Kt1    (REAL*8)  : Rate constant [s-1] for sulfate formation on 
!                           fine sea-salt aerosols from GET_ALK
!  (8 ) Kt2    (REAL*8)  : Rate constant [s-1] for sulfate formation on 
!                           coarse sea-salt aerosols from GET_ALK
!
!  Arguments as Output:
!  ============================================================================
!  (9 ) SO2_ss 	(REAL*8) : SO2 mixing ratio [v/v] updated after SS chemistry
!  (10) SO4E    (REAL*8) : SO4E (sulfate produced by S(IV)+O3 on fine seasalt)
!                          mixing ratio [v/v]
!  (11) SO4F    (REAL*8) : SO4F(sulfate produced by S(IV)+O3 on coarse seasalt)
!                          mixing ratio [v/v]
!  (12) O3      (REAL*8) : Updated O3 mixing ratio [v/v] for fullchem runs 
!                           only. Otherwise O30=O3.
!
!  Chemical reactions:
!  ============================================================================
!  (R1) SO2 + O3 + ALK => SO4 + O2
!       Modeled after Chamedies and Stelson, 1992?
!
!  NOTES:
!  (1 ) Now references XNUMOLAIR from "tracer_mod.f" (bmy, 10/25/05)
!  (2 ) Bug fix: now avoid seg fault error if IDTHNO3 is zero, as it would
!        be for an offline aerosol simulation. (bmy, 3/29/06)
!  (3 ) Fixed typo in FALK_A_SO2 equation: C_FLUX_C should be C_FLUX_A.
!        (havala, bec, bmy, 12/8/06)
!  (4 ) Bug fix for mass balance, replace TITR_HNO3 w/ HNO3_SSC in the
!        expression for HNO3_ss.  Bug fix: now do equivalent computation 
!        for GET_GNO3, which is now no longer called because it's in 
!        "isoropia_mod.f". (bec, bmy, 7/30/08)
!******************************************************************************
!
      ! References to F90 modules
      USE COMODE_MOD,      ONLY : CSPEC, JLOP, VOLUME
      USE DAO_MOD,         ONLY : AD, AIRDEN, AIRVOL
      USE TRACERID_MOD
      !----------------------------------------------------------------
      ! DIAGNOSTIC -- leave commented out for now (bec, bmy, 4/13/05)
      !USE DIAG_MOD,        ONLY : AD09
      !----------------------------------------------------------------
      USE ERROR_MOD,       ONLY : GEOS_CHEM_STOP
      USE TIME_MOD,        ONLY : GET_TS_CHEM,        GET_ELAPSED_SEC
      USE ERROR_MOD,       ONLY : IT_IS_NAN
      USE TRACER_MOD,      ONLY : ITS_A_FULLCHEM_SIM, STT
      USE TRACER_MOD,      ONLY : TCVV,               XNUMOLAIR
      USE GLOBAL_HNO3_MOD, ONLY : GET_HNO3_UGM3
      USE TIME_MOD,        ONLY : GET_ELAPSED_SEC,    GET_MONTH 
      USE TIME_MOD,        ONLY : ITS_A_NEW_MONTH

      ! Add these for GET_GNO3 fix (lyj, bmy, 10/7/08)
      USE GLOBAL_HNO3_MOD, ONLY : GET_HNO3_UGM3
      USE DAO_MOD,         ONLY : AIRVOL

#     include "F77_CMN_SIZE"  ! Size parameters
!---------------------------------------------------------------
! DIAGNOSTIC -- leave commented out for now (bec, bmy, 4/13/05)
!#     include "F77_CMN_DIAG"  ! ND19
!---------------------------------------------------------------
#     include "F77_CMN_GCTM"  ! AIRMW

      ! Arguments
      INTEGER, INTENT(IN)  :: I, J, L 
      REAL*8,  INTENT(IN)  :: SO2_cd, Kt1, Kt2, Kt1N, Kt2N
      REAL*8,  INTENT(IN)  :: ALK1, ALK2
      REAL*8,  INTENT(OUT) :: SO2_ss, PSO4E, PSO4F

      ! Local variables
      INTEGER              :: JLOOP
      REAL*8               :: SO2_chem,    DTCHEM
      REAL*8               :: O3_cspec,    O3_lost
      REAL*8               :: EQ_1_C,      EQ_2_C
      REAL*8               :: SO4E,        SO2_new,    SO4F
      REAL*8               :: SO2_eq,      N_FLUX_A,   N_FLUX_C
      REAL*8               :: END_ALK,     L5A,        L5C
      REAL*8               :: EQ1,         EQ2,        TITR_SO2
      REAL*8               :: TITR_HNO3,   NIT_vv,     NITs_vv
      REAL*8               :: NIT0,        NITS0
      REAL*8               :: F_SO2,       FALK_A_SO2, FALK_C_SO2
      REAL*8               :: EQ_BEG,      F_SO2_A,    F_SO2_C
      REAL*8               :: ALKA,        ALKC,       TOTAL_ACID_FLUX
      REAL*8               :: HNO3_EQ,     TOT_FLUX_A, TOT_FLUX_C
      REAL*8               :: FALK_A_HNO3, HNO3_vv
      REAL*8               :: FALK_C_HNO3, F_HNO3_A,   F_HNO3_C
      REAL*8               :: EQ_1_N,      EQ_2_N,     F_HNO3
      REAL*8               :: HNO3_SSA,    HNO3_SSC,   N_FLUX
      REAL*8               :: HNO3_EQ_C,   L6A,        L6C   
      REAL*8               :: C_FLUX_A,    C_FLUX_C,   C_FLUX
      REAL*8               :: HNO3_ss,     HNO3_kg
      REAL*8,  PARAMETER   :: MINDAT    = 1.0d-20 
      REAL*8,  PARAMETER   :: TCVV_HNO3 = 28.97d0 / 63.0d0 

      !=================================================================
      ! SEASALT_CHEM begins here!
      !=================================================================

      ! DTCHEM is the chemistry timestep in seconds
      DTCHEM = GET_TS_CHEM() * 60d0

      ! Convert SO2 [v/v] to  [eq/gridbox]
      SO2_eq = ( 2.d0 * SO2_cd * AD(I,J,L) ) / ( TCVV(IDTSO2) * 0.064d0)
      SO2_eq = MAX( SO2_eq, MINDAT )

      IF ( ITS_A_FULLCHEM_SIM() ) THEN

 	 ! Convert HNO3 [v/v] to [equivalents]
         HNO3_vv = STT(I,J,L,IDTHNO3)
         HNO3_eq = HNO3_vv * AD(I,J,L) / ( 28.97d0 / 63d0 ) / 63.d-3

      ELSE
         
         ! Get HNO3 in ug/m3, then multiply by volume in m3
         ! and 1e-6 kg/ug to get HNO3 in kg
         HNO3_kg = GET_HNO3_UGM3( I, J, L ) * AIRVOL(I,J,L) * 1e-6

	 ! Convert HNO3 [kg] first to [v/v] 
         HNO3_vv = HNO3_kg * ( 28.97d0 / 63d0 )   / AD(I,J,L)
 
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
      N_FLUX_A = F_HNO3_A * AD(I,J,L)/( 28.97d0 / 63d0 )/0.063d0

      ! Available flux of HNO3 to coarse sea salt aerosols [v/v/timestep]
      L6C = EXP( - Kt2N * DTCHEM )
      F_HNO3_C = HNO3_vv * ( 1.D0 - L6C )
      F_HNO3_C = MAX( F_HNO3_C, 1.0D-32 )

      ! convert to [eq/timestep] 
      N_FLUX_C = F_HNO3_C * AD(I,J,L)/( 28.97d0 / 63d0 )/0.063d0

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
      SO2_ss       = SO2_chem * 0.064 * TCVV(IDTSO2)/AD(I,J,L)/2.0d0
      SO2_ss       = MAX( SO2_ss, MINDAT )

      !SO4E produced converted from [eq/timestep] to [v/v/timestep] 
      PSO4E        = SO4E * 0.096 * TCVV(IDTSO4)/AD(I,J,L)/2.0d0

      !SO4F produced converted from [eq/timestep] to [v/v/timestep] 
      PSO4F        = SO4F * 0.096 * TCVV(IDTSO4S)/AD(I,J,L)/2.0d0

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
         HNO3_ss = HNO3_SSC * 0.063 * TCVV(IDTHNO3)/AD(I,J,L)
         STT(I,J,L,IDTHNO3) = MAX( HNO3_vv - HNO3_ss, MINDAT )

      ELSE

         ! Offline aerosol sim: IDTHNO3 isn't defined, use TCVV_HNO3
         HNO3_ss = TITR_HNO3 * 0.063 * TCVV_HNO3 / AD(I,J,L)

      ENDIF

      ! NITS produced converted from [eq/timestep] to [v/v/timestep] 
      PNITs(I,J,L) = HNO3_SSC * 0.063 * TCVV(IDTNITS)/AD(I,J,L)
	 
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

      ! Return to calling program
      END SUBROUTINE SEASALT_CHEM

!------------------------------------------------------------------------------

      SUBROUTINE AQCHEM_SO2( LWC, T,     P,       SO2, H2O2, 
     &                       O3,  Hplus, KaqH2O2, KaqO3 ) 
!
!******************************************************************************
!  Function AQCHEM_SO2 computes the reaction rates for aqueous SO2 chemistry.
!  (rjp, bmy, 10/31/02, 12/12/02)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) LWC     (REAL*8) : Liquid water content [m3/m3] = 1.E-6*L [g/m3]
!  (2 ) T       (REAL*8) : Temperature [K]
!  (3 ) P       (REAL*8) : Pressure [atm]
!  (4 ) SO2     (REAL*8) : SO2  mixing ratio [v/v]
!  (5 ) H2O2    (REAL*8) : H2O2 mixing ratio [v/v]
!  (6 ) O3      (REAL*8) : O3   mixing ratio [v/v]
!  (7 ) HPLUS   (REAL*8) : Concentration of H+ ion (i.e. the pH) [v/v]
!
!  Arguments as Output:
!  ============================================================================
!  (8 ) KaqH2O2 (REAL*8) : Reaction rate for H2O2
!  (9 ) KaqO3   (REAL*8) : Reaction rate for O3
!
!  Chemical Reactions:
!  ============================================================================
!  (R1) HSO3- + H2O2(aq) + H+ => SO4-- + 2H+ + H2O [Jacob, 1986]   
!
!      d[S(VI)]/dt = k[H+][H2O2(aq)][HSO3-]/(1 + K[H+]) 
!      [Seinfeld and Pandis, 1998, page 366]
!
!  (R2) SO2(aq) + O3(aq) =>                                        
!       HSO3-   + O3(aq) =>  
!       SO3--   + O3(aq) =>
!       [Jacob, 1986; Jacobson, 1999]
!
!       d[S(VI)]/dt = (k0[SO2(aq)] + k1[HSO3-] + K2[SO3--])[O3(aq)]
!       [Seinfeld and Pandis, 1998, page 363]
!
!  Reaction rates can be given as
!       Ra     = k [H2O2(ag)] [S(IV)]  [mole/liter*s]  OR
!       Krate  = Ra LWC R T / P        [1/s]
!
!  Where:
!       LWC = Liquid water content(g/m3)*10-6 [m3(water)/m3(gas)]
!       R   = 0.08205  (atm L / mol-K), Universal gas const.
!       T   = Temperature (K)
!       P   = Pressure (atm)
!
!  Procedure:
!  ============================================================================
!  (a ) Given [SO2] which is assumed to be total SO2 (gas+liquid) in 
!        equilibrium between gas and liquid phase. 
!
!  (b ) We can compute SO2(g) using Henry's law 
!          P(so2(g)) = Xg * [SO2]
!          Xg = 1/(1 + Faq), Fraction of SO2 in gas
!       where: 
!          Faq   = Kheff * R * T * LWC, 
!          KHeff = Effective Henry's constant
!
!  (c ) Then Calculate Aquous phase, S[IV] concentrations
!        S[IV] = Kheff * P(so2(g) in atm) [M]
!
!  (d ) The exact same procedure is applied to calculate H2O2(aq)
!
!  NOTES:
!  (1 ) Updated by Rokjin Park (rjp, bmy, 12/12/02)
!******************************************************************************
!
      ! Arguments
      REAL*8, INTENT(IN)  :: LWC, T, P, SO2, H2O2, O3, HPLUS
      REAL*8, INTENT(OUT) :: KaqH2O2, KaqO3

      ! Local variables
      REAL*8, PARAMETER   :: R = 0.08205d0 
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

      ! Return to calling program
      END SUBROUTINE AQCHEM_SO2

!------------------------------------------------------------------------------

      SUBROUTINE CHEM_SO4
!
!******************************************************************************
!  Subroutine CHEM_SO4 is the SO4 chemistry subroutine from Mian Chin's GOCART
!  model, modified for the GEOS-CHEM model.  Now also modified to account
!  for production of crystalline & aqueous sulfur tracers. 
!  (rjp, bdf, cas, bmy, 5/31/00, 5/23/06) 
!                                                                          
!  Module Variables Used:
!  ============================================================================
!  (1 ) PSO4_SO2 (REAL*8 ) : Array for P(SO4) from SO2 [v/v/timestep]
!  (2 ) PSO4_ss  (REAL*8 ) : Array for P(SO4) from SO2 
!                            (coarse sea-salt aerosols) [v/v/timestep]
!                                                                           
!  Reaction List (by Mian Chin, chin@rondo.gsfc.nasa.gov)                  
!  ============================================================================
!  The Only production is from SO2 oxidation (save in CHEM_SO2), and the only  
!  loss is dry depsition here.  Wet deposition will be treated in "wetdep.f".
!                                                                          
!  SO4 = SO4_0 * exp(-kt) + PSO4_SO2/kt * (1.-exp(-kt))                    
!    where k = dry deposition.                                             
!                      
!  NOTES:              
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
!******************************************************************************
!
      ! References to F90 modules
      USE DAO_MOD,        ONLY : AD
      USE DIAG_MOD,       ONLY : AD44
      USE DRYDEP_MOD,     ONLY : DEPSAV
      USE GRID_MOD,       ONLY : GET_AREA_CM2
      USE LOGICAL_MOD,    ONLY : LCRYST, LSSALT
      USE PBL_MIX_MOD,    ONLY : GET_FRAC_UNDER_PBLTOP
      USE TIME_MOD,       ONLY : GET_TS_CHEM
      USE TRACER_MOD,     ONLY : STT,    TCVV,     XNUMOL
      USE TRACERID_MOD,   ONLY : IDTSO4, IDTSO4s,  IDTAS,   IDTAHS 
      USE TRACERID_MOD,   ONLY : IDTLET, IDTSO4aq, IDTNH4aq
      USE TROPOPAUSE_MOD, ONLY : ITS_IN_THE_STRAT

#     include "F77_CMN_SIZE"     ! Size parameters
#     include "F77_CMN_DIAG"     ! ND44

      ! Local variables
      INTEGER               :: I,      J,      L,        N,     N_ND44
      REAL*8                :: AS,     AS0,    AHS,      AHS0,  LET   
      REAL*8                :: LET0,   SO4,    SO40,     SO4aq, SO4aq0
      REAL*8                :: SO4s,   SO40s,  RKT,      RKTs,  E_RKT
      REAL*8                :: E_RKTs, DTCHEM, AREA_CM2, FLUX
      REAL*8                :: T44(IIPAR,JJPAR,LLTROP,6) 

      !=================================================================
      ! CHEM_SO4 begins here!
      !=================================================================

      ! Return if tracers are not defined
      IF ( IDTSO4 == 0 .or. IDTSO4s == 0 ) RETURN
      IF ( DRYSO4 == 0 .or. DRYSO4s == 0 ) RETURN

      ! DTCHEM is the chemistry timestep in seconds
      DTCHEM = GET_TS_CHEM() * 60d0

      ! Number of drydep tracers to save
      N_ND44 = 2

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
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, N )
         DO N = 1, N_ND44
         DO L = 1, LLTROP 
         DO J = 1, JJPAR
         DO I = 1, IIPAR
            T44(I,J,L,N) = 0d0
         ENDDO
         ENDDO
         ENDDO
         ENDDO
!$OMP END PARALLEL DO
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
      DO L = 1, LLTROP 
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
         IF ( ITS_IN_THE_STRAT( I, J, L ) ) CYCLE

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
         IF (LNLPBL) RKT = 0.D0

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
         IF (LNLPBL) RKTs = 0.D0
        
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
            AREA_CM2 = GET_AREA_CM2( J )

            ! SO4 drydep flux [molec/cm2/s]
            FLUX = SO40  - SO4 + PSO4_SO2(I,J,L) 
            FLUX = FLUX  * AD(I,J,L)       / TCVV(IDTSO4)
            FLUX = FLUX  * XNUMOL(IDTSO4)  / AREA_CM2 / DTCHEM
            T44(I,J,L,1) = T44(I,J,L,1) + FLUX

            ! SO4s drydep flux [molec/cm2/s]
            FLUX = SO40s - SO4s + PSO4_ss(I,J,L) 
            FLUX = FLUX  * AD(I,J,L)       / TCVV(IDTSO4s)
            FLUX = FLUX  * XNUMOL(IDTSO4s) / AREA_CM2 / DTCHEM
            T44(I,J,L,2) = T44(I,J,L,2) + FLUX

!------------------------------------------------------------------------------
!%%% Currently under development (rjp, bmy, 3/15/05)
!%%%             ! SO4aq, AS, AHS, LET drydep fluxes (if necessary)
!%%%             IF ( LCRYST ) THEN
!%%% 
!%%%                ! SO4aq drydep flux [molec/cm2/s]
!%%%                FLUX = SO4aq0 - SO4aq + PSO4_SO2(I,J,L)     
!%%%                FLUX = FLUX  * AD(I,J,L)        / TCVV(IDTSO4aq)
!%%%                FLUX = FLUX  * XNUMOL(IDTSO4aq) / AREA_CM2 / DTCHEM
!%%%                T44(I,J,L,3) = T44(I,J,L,3) + FLUX
!%%% 
!%%%                ! AS drydep flux [molec/cm2/s]
!%%%                FLUX = AS0   - AS
!%%%                FLUX = FLUX  * AD(I,J,L)     / TCVV(IDTAS)
!%%%                FLUX = FLUX  * XNUMOL(IDTAS) / AREA_CM2 / DTCHEM
!%%%                T44(I,J,L,4) = T44(I,J,L,4) + FLUX
!%%% 
!%%%                ! AHS drydep flux [molec/cm2/s]
!%%%                FLUX = AHS0  - AHS
!%%%                FLUX = FLUX  * AD(I,J,L)      / TCVV(IDTAHS)
!%%%                FLUX = FLUX  * XNUMOL(IDTAHS) / AREA_CM2 / DTCHEM
!%%%                T44(I,J,L,5) = T44(I,J,L,5) + FLUX
!%%% 
!%%%                ! LET drydep flux [molec/cm2/s]
!%%%                FLUX = LET0  - LET
!%%%                FLUX = FLUX  * AD(I,J,L)      / TCVV(IDTLET)
!%%%                FLUX = FLUX  * XNUMOL(IDTLET) / AREA_CM2 / DTCHEM
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
         DO L = 1, LLTROP

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

      ! Return to calling program
      END SUBROUTINE CHEM_SO4

!------------------------------------------------------------------------------

      !SUBROUTINE PHASE_SO4
      !
      ! *** Currently under development ***
      !
      !END SUBROUTINE PHASE_SO4

!------------------------------------------------------------------------------

      !SUBROUTINE PHASE_RADIATIVE
      !
      ! *** Currently under development ***
      !
      !END SUBROUTINE PHASE_RADIATIVE

!------------------------------------------------------------------------------

      SUBROUTINE CHEM_MSA
!
!******************************************************************************
!  Subroutine CHEM_MSA is the SO4 chemistry subroutine from Mian Chin's GOCART
!  model, modified for the GEOS-CHEM model. (rjp, bdf, bmy, 5/31/00, 10/25/05)
!                                                                          
!  Module Variables Used:
!  ============================================================================
!  (1 ) PMSA_DMS (REAL*8 ) : Array for P(MSA) from DMS [v/v/timestep]
!                                                                          
!  Reaction List (by Mian Chin, chin@rondo.gsfc.nasa.gov)                  
!  ============================================================================
!  The Only production is from DMS oxidation (saved in CHEM_DMS), and the only
!  loss is dry depsition here.  Wet deposition will be treated in "wetdep.f".
!                                                                          
!  MSA = MSA_0 * exp(-dt) + PMSA_DMS/kt * (1.-exp(-kt))                    
!    where k = dry deposition.                                             
!        
!  NOTES:
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
!******************************************************************************
!
      ! References to F90 modules
      USE DAO_MOD,        ONLY : AD
      USE DIAG_MOD,       ONLY : AD44
      USE DRYDEP_MOD,     ONLY : DEPSAV
      USE GRID_MOD,       ONLY : GET_AREA_CM2
      USE PBL_MIX_MOD,    ONLY : GET_FRAC_UNDER_PBLTOP, GET_PBL_MAX_L
      USE TIME_MOD,       ONLY : GET_TS_CHEM
      USE TRACER_MOD,     ONLY : STT, TCVV, XNUMOL
      USE TRACERID_MOD,   ONLY : IDTMSA
      USE TROPOPAUSE_MOD, ONLY : ITS_IN_THE_STRAT

#     include "F77_CMN_SIZE"     ! Size parameters
#     include "F77_CMN_GCTM"     ! AIRMW
#     include "F77_CMN_DIAG"     ! ND44

      ! Local variables
      INTEGER               :: I,      J,    L,        PBL_MAX
      REAL*8                :: DTCHEM, MSA0, MSA,      RK       
      REAL*8                :: RKT,    FLUX, AREA_CM2, F_UNDER_TOP
      REAL*8                :: ND44_TMP(IIPAR,JJPAR,LLTROP)

      !=================================================================
      ! CHEM_MSA begins here!
      !=================================================================
      IF ( IDTMSA == 0 .or. DRYMSA == 0 ) RETURN

      ! DTCHEM is the chemistry interval in seconds
      DTCHEM  = GET_TS_CHEM() * 60d0 

      ! Model level where maximum PBL height occurs 
      PBL_MAX = GET_PBL_MAX_L()

      ! Zero ND44_TMP array
      IF ( ND44 > 0 ) THEN
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
!$OMP+SCHEDULE( DYNAMIC )
         DO L = 1, LLTROP 
         DO J = 1, JJPAR
         DO I = 1, IIPAR
            ND44_TMP(I,J,L) = 0d0
         ENDDO
         ENDDO
         ENDDO
!$OMP END PARALLEL DO
      ENDIF

      ! Loop over tropospheric grid boxes
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, F_UNDER_TOP, MSA0, RKT, MSA, AREA_CM2, FLUX )
!$OMP+SCHEDULE( DYNAMIC )
      DO L = 1, LLTROP
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Skip stratospheric boxes
         IF ( ITS_IN_THE_STRAT( I, J, L ) ) CYCLE      

         ! Fraction of box (I,J,L) underneath the PBL top [unitless]
         F_UNDER_TOP = GET_FRAC_UNDER_PBLTOP( I, J, L )

         ! Initial MSA [v/v]
         MSA0 = STT(I,J,L,IDTMSA) 

         ! MSA drydep frequency [1/s].  Also accounts for the fraction
         ! of each grid box (I,J,L) that is located beneath the PBL top
         RKT = DEPSAV(I,J,DRYMSA) * F_UNDER_TOP

         ! Add option for non-local PBL (Lin, 03/31/09)
         IF (LNLPBL) RKT = 0.D0
         
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
            AREA_CM2 = GET_AREA_CM2( J )

            ! Convert [v/v/timestep] to [molec/cm2/s]
            FLUX = MSA0 - MSA + PMSA_DMS(I,J,L)                    
            FLUX = FLUX * AD(I,J,L)      / TCVV(IDTMSA)            
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
         DO L = 1, PBL_MAX
            AD44(I,J,DRYMSA,1) = AD44(I,J,DRYMSA,1) + ND44_TMP(I,J,L)
         ENDDO
         ENDDO
         ENDDO
!$OMP END PARALLEL DO
      ENDIF

      ! Return to calling program
      END SUBROUTINE CHEM_MSA

!------------------------------------------------------------------------------

      SUBROUTINE CHEM_NH3
!
!******************************************************************************
!  Subroutine CHEM_NH3 removes NH3 from the surface via dry deposition.
!  (rjp, bdf, bmy, 1/2/02, 10/25/05)  
!                                                                          
!  Reaction List:
!  ============================================================================
!  (1 ) NH3 = NH3_0 * EXP( -dt )  where d = dry deposition rate [s-1]
!        
!  NOTES:
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
!******************************************************************************
!
      ! References to F90 modules
      USE DAO_MOD,      ONLY : AD
      USE DIAG_MOD,     ONLY : AD44
      USE DRYDEP_MOD,   ONLY : DEPSAV
      USE GRID_MOD,     ONLY : GET_AREA_CM2
      USE PBL_MIX_MOD,  ONLY : GET_FRAC_UNDER_PBLTOP, GET_PBL_MAX_L
      USE TIME_MOD,     ONLY : GET_TS_CHEM
      USE TRACER_MOD,   ONLY : STT, TCVV, XNUMOL
      USE TRACERID_MOD, ONLY : IDTNH3

#     include "F77_CMN_SIZE"     ! Size parameters
#     include "F77_CMN_DIAG"     ! ND44

      ! Local variables
      INTEGER :: I,      J,        L,    PBL_MAX
      REAL*8  :: DTCHEM, NH30,     NH3
      REAL*8  :: FREQ,   AREA_CM2, FLUX, F_UNDER_TOP
      REAL*8  :: ND44_TMP(IIPAR,JJPAR,LLTROP)

      !=================================================================
      ! CHEM_NH3 begins here!
      !=================================================================
      IF ( IDTNH3 == 0 .or. DRYNH3 == 0 ) RETURN

      ! DTCHEM is the chemistry interval in seconds
      DTCHEM  = GET_TS_CHEM() * 60d0

      ! Model level where maximum PBL height occurs 
      PBL_MAX = GET_PBL_MAX_L()

      ! Zero ND44_TMP array
      IF ( ND44 > 0 ) THEN
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
!$OMP+SCHEDULE( DYNAMIC )
         DO L = 1, LLTROP 
         DO J = 1, JJPAR
         DO I = 1, IIPAR
            ND44_TMP(I,J,L) = 0d0
         ENDDO
         ENDDO
         ENDDO
!$OMP END PARALLEL DO
      ENDIF

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, F_UNDER_TOP, FREQ, NH30, NH3, AREA_CM2, FLUX )
!$OMP+SCHEDULE( DYNAMIC )
      DO L = 1, PBL_MAX
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
         IF (LNLPBL) FREQ = 0.D0

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
               ! ND44 diagnostic: Drydep flux of NH3 [molec/cm2/s]
               !========================================================
               IF ( ND44 > 0 .and. NH3 > 0d0 ) THEN

                  ! Surface area [cm2]
                  AREA_CM2 = GET_AREA_CM2( J )
                  
                  ! Convert drydep loss from [v/v/timestep] to [molec/cm2/s]
                  FLUX = NH3  * AD(I,J,L)      / TCVV(IDTNH3)
                  FLUX = FLUX * XNUMOL(IDTNH3) / AREA_CM2 / DTCHEM

                  ! Store dryd flx in ND44_TMP as a placeholder
                  ND44_TMP(I,J,L) = ND44_TMP(I,J,L) + FLUX
               ENDIF
            ENDIF
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
         DO L = 1, PBL_MAX
            AD44(I,J,DRYNH3,1) = AD44(I,J,DRYNH3,1) + ND44_TMP(I,J,L)
         ENDDO
         ENDDO
         ENDDO
!$OMP END PARALLEL DO
      ENDIF

      ! Return to calling program
      END SUBROUTINE CHEM_NH3

!------------------------------------------------------------------------------

      SUBROUTINE CHEM_NH4
!
!******************************************************************************
!  Subroutine CHEM_NH4 removes NH4 from the surface via dry deposition.
!  (rjp, bdf, bmy, 1/2/02, 10/25/05)  
!                                                                          
!  Reaction List:
!  ============================================================================
!  (1 ) NH4 = NH4_0 * EXP( -dt )  where d = dry deposition rate [s-1]
!        
!  NOTES:
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
!******************************************************************************
!
      ! References to F90 modules
      USE DAO_MOD,      ONLY : AD
      USE DIAG_MOD,     ONLY : AD44
      USE DRYDEP_MOD,   ONLY : DEPSAV
      USE GRID_MOD,     ONLY : GET_AREA_CM2
      USE PBL_MIX_MOD,  ONLY : GET_FRAC_UNDER_PBLTOP, GET_PBL_MAX_L
      USE TIME_MOD,     ONLY : GET_TS_CHEM
      USE TRACER_MOD,   ONLY : STT, TCVV, XNUMOL
      USE TRACERID_MOD, ONLY : IDTNH4

#     include "F77_CMN_SIZE"  ! Size parameters
#     include "F77_CMN_DIAG"  ! ND44

      ! Local variables
      INTEGER :: I,      J,    L,        PBL_MAX
      REAL*8  :: DTCHEM, NH4,  NH40
      REAL*8  :: FREQ,   FLUX, AREA_CM2, F_UNDER_TOP
      REAL*8  :: ND44_TMP(IIPAR,JJPAR,LLTROP)

      !=================================================================
      ! CHEM_NH4 begins here!
      !=================================================================
      IF ( IDTNH4 == 0 .or. DRYNH4 == 0 ) RETURN

      ! DTCHEM is the chemistry interval in seconds
      DTCHEM = GET_TS_CHEM() * 60d0 

      ! Model level where maximum PBL height occurs 
      PBL_MAX = GET_PBL_MAX_L()

      ! Zero ND44_TMP array
      IF ( ND44 > 0 ) THEN
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
         DO L = 1, LLTROP 
         DO J = 1, JJPAR
         DO I = 1, IIPAR
            ND44_TMP(I,J,L) = 0d0
         ENDDO
         ENDDO
         ENDDO
!$OMP END PARALLEL DO
      ENDIF

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, F_UNDER_TOP, FREQ, NH40, NH4, AREA_CM2, FLUX )
!$OMP+SCHEDULE( DYNAMIC )
      DO L = 1, PBL_MAX
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
         IF (LNLPBL) FREQ = 0.D0

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
               ! ND44 diagnostic: Drydep flux of NH4 [molec/cm2/s]
               !========================================================
               IF ( ND44 > 0 .and. NH4 > 0d0 ) THEN
         
                  ! Surface area [cm2]
                  AREA_CM2 = GET_AREA_CM2( J )
                  
                  ! Convert drydep loss from [v/v/timestep] to [molec/cm2/s]
                  FLUX = NH4  * AD(I,J,L)      / TCVV(IDTNH4)
                  FLUX = FLUX * XNUMOL(IDTNH4) / AREA_CM2 / DTCHEM

                  ! Store dryd flx in ND44_TMP as a placeholder
                  ND44_TMP(I,J,L) = ND44_TMP(I,J,L) + FLUX
               ENDIF
            ENDIF
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
         DO L = 1, PBL_MAX
            AD44(I,J,DRYNH4,1) = AD44(I,J,DRYNH4,1) + ND44_TMP(I,J,L)
         ENDDO
         ENDDO
         ENDDO
!$OMP END PARALLEL DO
      ENDIF

      ! Return to calling program
      END SUBROUTINE CHEM_NH4

!------------------------------------------------------------------------------

      SUBROUTINE CHEM_NH4aq
!
!******************************************************************************
!  Subroutine CHEM_NH4aq removes NH4aq from the surface via dry deposition.
!  (cas, bmy, 1/6/05, 10/25/05)  
!                                                                          
!  Reaction List:
!  ============================================================================
!  (1 ) NH4aq = NH4_0aq * EXP( -dt )  where d = dry deposition rate [s-1]
!        
!  NOTES:
!  (1 ) Replace PBLFRAC from "drydep_mod.f" with GET_FRAC_UNDER_PBLTOP from 
!        "pbl_mix_mod.f".  Also reference GET_PBL_MAX_L from "pbl_mix_mod.f"
!        Vertical DO-loops can run up to PBL_MAX and not LLTROP. (bmy, 2/22/05)
!  (31) Now references XNUMOL from "tracer_mod.f" (bmy, 10/25/05)
!******************************************************************************
!
      ! References to F90 modules
      USE DAO_MOD,      ONLY : AD
      USE DIAG_MOD,     ONLY : AD44
      USE DRYDEP_MOD,   ONLY : DEPSAV
      USE GRID_MOD,     ONLY : GET_AREA_CM2
      USE PBL_MIX_MOD,  ONLY : GET_FRAC_UNDER_PBLTOP, GET_PBL_MAX_L
      USE TIME_MOD,     ONLY : GET_TS_CHEM
      USE TRACER_MOD,   ONLY : STT, TCVV, XNUMOL
      USE TRACERID_MOD, ONLY : IDTNH4aq

#     include "F77_CMN_SIZE"  ! Size parameters
#     include "F77_CMN_DIAG"  ! ND44

      ! Local variables
      INTEGER :: I,      J,     L,        PBL_MAX
      REAL*8  :: DTCHEM, NH4aq, NH4aq0
      REAL*8  :: FREQ,   FLUX,  AREA_CM2, F_UNDER_TOP
      REAL*8  :: T44(IIPAR,JJPAR,LLTROP)

      !=================================================================
      ! CHEM_NH4 begins here!
      !=================================================================
      IF ( IDTNH4aq == 0 .or. DRYNH4aq == 0 ) RETURN

      ! DTCHEM is the chemistry interval in seconds
      DTCHEM  = GET_TS_CHEM() * 60d0 

      ! Model level where maximum PBL height occurs 
      PBL_MAX = GET_PBL_MAX_L()      

      ! Zero T44 array
      IF ( ND44 > 0 ) THEN
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
         DO L = 1, LLTROP 
         DO J = 1, JJPAR
         DO I = 1, IIPAR
            T44(I,J,L) = 0d0
         ENDDO
         ENDDO
         ENDDO
!$OMP END PARALLEL DO
      ENDIF

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I,      J,     L,        F_UNDER_TOP, FREQ  )
!$OMP+PRIVATE( NH4aq0, NH4aq, AREA_CM2, FLUX               )
!$OMP+SCHEDULE( DYNAMIC )
      DO L = 1, PBL_MAX
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Fraction of box (I,J,L) underneath the PBL top [unitless]
         F_UNDER_TOP = GET_FRAC_UNDER_PBLTOP( I, J, L )     

         ! Only apply drydep to boxes w/in the PBL
         IF ( F_UNDER_TOP > 0d0 ) THEN

            ! NH4 drydep frequency [1/s] -- PBLFRAC accounts for the fraction
            ! of each grid box (I,J,L) that is located beneath the PBL top
            FREQ = DEPSAV(I,J,DRYNH4aq) * F_UNDER_TOP

         ! Add option for non-local PBL (Lin, 03/31/09)
         IF (LNLPBL) FREQ = 0.D0

            ! Only apply drydep loss if FREQ is nonzero
            IF ( FREQ > 0d0 ) THEN

               ! Initial NH4 [v/v]
               NH4aq0 = STT(I,J,L,IDTNH4aq)
         
               ! Amount of NH4 lost to drydep [v/v]
               NH4aq = NH4aq0 * ( 1d0 - EXP( -FREQ * DTCHEM ) )

               ! Prevent underflow condition
               IF ( NH4aq < SMALLNUM ) NH4aq = 0d0

               ! Subtract NH4 lost to drydep from initial NH4 [v/v]
               STT(I,J,L,IDTNH4aq) = NH4aq0 - NH4aq

               !========================================================
               ! ND44 diagnostic: Drydep flux of NH4 [molec/cm2/s]
               !========================================================
               IF ( ND44 > 0 .and. NH4aq > 0d0 ) THEN
         
                  ! Surface area [cm2]
                  AREA_CM2 = GET_AREA_CM2( J )
                  
                  ! Convert drydep loss from [v/v/timestep] to [molec/cm2/s]
                  FLUX = NH4aq * AD(I,J,L)        / TCVV(IDTNH4aq)
                  FLUX = FLUX  * XNUMOL(IDTNH4aq) / AREA_CM2 / DTCHEM

                  ! Store dryd flx in ND44_TMP as a placeholder
                  T44(I,J,L) = T44(I,J,L) + FLUX
               ENDIF
            ENDIF
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
         DO L = 1, PBL_MAX
            AD44(I,J,DRYNH4aq,1) = AD44(I,J,DRYNH4aq,1) + T44(I,J,L)
         ENDDO
         ENDDO
         ENDDO
!$OMP END PARALLEL DO
      ENDIF

      ! Return to calling program
      END SUBROUTINE CHEM_NH4aq

!------------------------------------------------------------------------------

      SUBROUTINE CHEM_NIT
!
!******************************************************************************
!  Subroutine CHEM_NIT removes SULFUR NITRATES (NIT) from the surface 
!  via dry deposition. (rjp, bdf, bmy, 1/2/02, 5/23/06)  
!                                                                          
!  Reaction List:
!  ============================================================================
!  (1 ) NIT = NIT_0 * EXP( -dt )  where d = dry deposition rate [s-1]
!        
!  NOTES:
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
!******************************************************************************
!
      ! References to F90 modules
      USE DAO_MOD,      ONLY : AD
      USE DIAG_MOD,     ONLY : AD44
      USE DRYDEP_MOD,   ONLY : DEPSAV
      USE GRID_MOD,     ONLY : GET_AREA_CM2
      USE LOGICAL_MOD,  ONLY : LSSALT
      USE PBL_MIX_MOD,  ONLY : GET_FRAC_UNDER_PBLTOP, GET_PBL_MAX_L
      USE TIME_MOD,     ONLY : GET_TS_CHEM
      USE TRACER_MOD,   ONLY : STT, TCVV, XNUMOL
      USE TRACERID_MOD, ONLY : IDTNIT, IDTNITs

#     include "F77_CMN_SIZE"     ! Size parameters
#     include "F77_CMN_DIAG"     ! ND44

      ! Local variables
      INTEGER :: I,        J,     L,    N,     N_ND44, PBL_MAX
      REAL*8  :: DTCHEM,   NIT,   NITs, NIT0,  NIT0s,  E_RKT
      REAL*8  :: E_RKTs,   FLUX,  FREQ, FREQs, RKT,    RKTs   
      REAL*8  :: AREA_CM2, F_UNDER_TOP
      REAL*8  :: T44(IIPAR,JJPAR,LLTROP,2)

      !=================================================================
      ! CHEM_NIT begins here!
      !=================================================================

      ! Return if tracers are not defined
      IF ( IDTNIT == 0 .or. IDTNITs == 0 ) RETURN
      IF ( DRYNIT == 0 .or. DRYNITs == 0 ) RETURN

      ! DTCHEM is the chemistry interval in seconds
      DTCHEM = GET_TS_CHEM() * 60d0 

      ! Model level where maximum PBL height occurs 
      PBL_MAX = GET_PBL_MAX_L()      
      
      ! Number of tracers for ND44
      N_ND44 = 2 

      ! Zero ND44 array
      IF ( ND44 > 0 ) THEN
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, N )
         DO N = 1, N_ND44
         DO L = 1, LLTROP 
         DO J = 1, JJPAR
         DO I = 1, IIPAR
            T44(I,J,L,N) = 0d0
         ENDDO
         ENDDO
         ENDDO
         ENDDO
!$OMP END PARALLEL DO
      ENDIF

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I,    J,      L,        NIT0,        NIT0s, NIT   )
!$OMP+PRIVATE( NITs, FREQ,   FREQs,    F_UNDER_TOP, RKT,   E_RKT ) 
!$OMP+PRIVATE( RKTs, E_RKTs, AREA_CM2, FLUX                      )
!$OMP+SCHEDULE( DYNAMIC )
      DO L = 1, PBL_MAX
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
         IF (LNLPBL) FREQ = 0.D0

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
         IF (LNLPBL) FREQs = 0.D0

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
            ! ND44 Diagnostic: Drydep flux of NIT [molec/cm2/s]
            !========================================================
            IF ( ND44 > 0 ) THEN
         
               ! Surface area [cm2]
               AREA_CM2 = GET_AREA_CM2( J )
                  
               !-------------
               ! NIT drydep
               !-------------
               
               ! If NIT drydep occurs ...
               IF ( FREQ > 0d0 ) THEN

                  ! Convert from [v/v/timestep] to [molec/cm2/s]
                  FLUX = NIT * AD(I,J,L)      / TCVV(IDTNIT)
                  FLUX = FLUX * XNUMOL(IDTNIT) / AREA_CM2 / DTCHEM
 
                  ! Store dryd flx in ND44_TMP as a placeholder
                  T44(I,J,L,1) = T44(I,J,L,1) + FLUX

               ENDIF
                  
               !-------------
               ! NITs drydep
               !-------------

               ! NOTE: if drydep doesn't occur then we still have
               ! production from seasalt (bec, bmy, 4/13/05)

               ! Convert from [v/v/timestep] to [molec/cm2/s]
               FLUX = NIT0s - NITs + PNITs(I,J,L) 
               FLUX = FLUX * AD(I,J,L)       / TCVV(IDTNITs)
               FLUX = FLUX * XNUMOL(IDTNITs) / AREA_CM2 / DTCHEM

               ! Store dryd flx in ND44_TMP as a placeholder
               T44(I,J,L,2) = T44(I,J,L,2) + FLUX

            ENDIF
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
         DO L = 1, PBL_MAX
            AD44(I,J,DRYNIT, 1) = AD44(I,J,DRYNIT, 1) + T44(I,J,L,1)
            AD44(I,J,DRYNITs,1) = AD44(I,J,DRYNITs,1) + T44(I,J,L,2)
         ENDDO
         ENDDO
         ENDDO
!$OMP END PARALLEL DO
      ENDIF

      ! Return to calling program
      END SUBROUTINE CHEM_NIT

!-----------------------------------------------------------------------------

      SUBROUTINE EMISSSULFATE
!
!******************************************************************************
!  Subroutine EMISSSULFATE is the interface between the GEOS-CHEM model and
!  the sulfate emissions routines in "sulfate_mod.f" (bmy, 6/7/00, 10/15/09)
! 
!  NOTES:
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
!  (7 ) Now references F77_CMN_SETUP.  Now read ship SO2 if LSHIPSO2=T.  Also
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
!******************************************************************************
!
      ! References to F90 modules
      USE ERROR_MOD,         ONLY : DEBUG_MSG
      USE LOGICAL_MOD,       ONLY : LSHIPSO2,   LPRT,   LBIOMASS !(win,5/1/09)
      USE TIME_MOD,          ONLY : GET_SEASON, GET_MONTH
      USE TIME_MOD,          ONLY : GET_YEAR,   ITS_A_NEW_MONTH
      USE TIME_MOD,          ONLY : GET_DAY,    ITS_A_NEW_DAY
      USE TRACER_MOD,        ONLY : STT,        ITS_AN_AEROSOL_SIM
      USE TRACERID_MOD,      ONLY : IDTNITs,    IDTSO4s
      USE TRACERID_MOD,      ONLY : IDTDMS,     IDTSO2 
      USE TRACERID_MOD,      ONLY : IDTSO4,     IDTNH3
      USE GFED2_BIOMASS_MOD, ONLY : GFED2_IS_NEW
      USE LOGICAL_MOD,       ONLY : LANTHRO, LBIOFUEL

#     include "F77_CMN_SIZE"  ! Size parameters

      ! Local variables
      LOGICAL, SAVE      :: FIRSTEMISS = .TRUE. 
      INTEGER            :: NSEASON, DAY, MONTH, YEAR   

      !=================================================================
      ! EMISSSULFATE begins here!
      !=================================================================

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

         ! Initialize arrays
         CALL INIT_SULFATE

         ! We have now gone thru the first timestep
         FIRSTEMISS  = .FALSE.
      ENDIF

      ! Get the season and month
      NSEASON = GET_SEASON()
      DAY     = GET_DAY()
      MONTH   = GET_MONTH()
      YEAR    = GET_YEAR()

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
         CALL READ_SST( MONTH, YEAR )
         CALL READ_OCEAN_DMS( MONTH )
         CALL READ_NATURAL_NH3( MONTH )

         ! Add LBIOFUEL szitch to turn off biofuel emissions.
         ! Warning: Streets 2006 inventory combines both anthropogenic
         ! and biofuel emissions, so with this inventory it is not
         ! possible to turn off biofuel emissions.
         !(ccc, 7/16/10)
         IF ( LBIOFUEL ) THEN
            CALL READ_BIOFUEL_SO2( MONTH )
            CALL READ_BIOFUEL_NH3( MONTH )
         ENDIF

         ! Add LANTHRO switch to turn off anthropogenic emissions.
         ! (ccc, 4/15/09)
         IF ( LANTHRO ) THEN
            CALL READ_AIRCRAFT_SO2( MONTH )
            CALL READ_ANTHRO_SOx( MONTH, NSEASON )
            CALL READ_ANTHRO_NH3( MONTH )
         ENDIF

         ! Also read ship exhaust SO2 if necessary 
         CALL READ_SHIP_SO2( MONTH )

         ! Read oxidants for the offline simulation only
         IF ( ITS_AN_AEROSOL_SIM() ) CALL READ_OXIDANT( MONTH )

      ENDIF

      IF(  ( GFED2_IS_NEW() .or. ITS_A_NEW_MONTH() ) .AND. 
     &     ( LBIOMASS )  ) THEN 
         CALL GET_BIOMASS_SO2
         IF ( LPRT ) CALL DEBUG_MSG( '### EMISSSULFATE: GET_BM_SO2')
         
         CALL GET_BIOMASS_NH3
         IF ( LPRT ) CALL DEBUG_MSG( '### EMISSSULFATE: GET_BM_NH3')         
      ENDIF

      !=================================================================
      ! Add emissions into the STT tracer array
      !=================================================================
      IF ( IDTDMS /= 0 ) CALL SRCDMS( STT(:,:,:,IDTDMS)          )
      IF ( IDTSO2 /= 0 ) CALL SRCSO2( STT(:,:,:,IDTSO2), NSEASON )
      IF ( IDTSO4 /= 0 ) CALL SRCSO4( STT(:,:,:,IDTSO4)          )
      IF ( IDTNH3 /= 0 ) CALL SRCNH3( STT(:,:,:,IDTNH3)          )

      ! Return to calling program
      END SUBROUTINE EMISSSULFATE
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
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
!   27 Oct 2009 - P. Le Sager - initial
!
! !REMARKS:
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
!---------------------------------------------------------------------------
      
      SUBROUTINE SRCDMS( TC )
!
!***************************************************************************** 
!  Subroutine SRCDMS, from Mian Chin's GOCART model, add DMS emissions 
!  to the tracer array.  Modified for use with the GEOS-CHEM model.
!  (bmy, 6/2/00, 8/16/05)
!
!  Arguments as Input/Output:
!  ===========================================================================
!  (1 ) TC (REAL*8 ) : Initial tracer mass [kg], plus DMS emissions
!
!  NOTES:
!  (1 ) Now reference NSRCE, LWI, DXYP, XTRA2 from either header files or
!        F90 modules.  Now use routines from "pressure_mod.f" to compute
!        grid box surface pressures. (bmy, 9/18/02)
!  (2 ) Now replace DXYP(J) with routine GET_AREA_M2 of "grid_mod.f"
!        Now use routine GET_TS_EMIS from the new "time_mod.f". (bmy, 3/27/03)
!  (3 ) For GEOS-4, convert PBL from [m] to [hPa] w/ the hydrostatic law.
!        Now references SCALE_HEIGHT from "F77_CMN_GCTM".  Added BLTHIK variable
!        for PBL thickness in [hPa]. (bmy, 1/15/04)
!  (4 ) Remove reference to "pressure_mod.f".  Now reference GET_FRAC_OF_PBL
!        and GET_PBL_TOP_L from "pbl_mix_mod.f". (bmy, 2/22/05)
!  (5 ) Switch from Liss & Merlivat to Nightingale formulation for DMS
!        emissions. (swu, bmy, 8/16/05)
!******************************************************************************
!
      ! Reference to diagnostic arrays
      USE DIAG_MOD,     ONLY : AD13_DMS
      USE DAO_MOD,      ONLY : IS_WATER, LWI, PBL
      USE GRID_MOD,     ONLY : GET_AREA_M2
      USE PBL_MIX_MOD,  ONLY : GET_FRAC_OF_PBL, GET_PBL_TOP_L
      USE TIME_MOD,     ONLY : GET_TS_EMIS
      USE TRACERID_MOD,   ONLY : IDTDMS ! (Lin, 03/31/09)

#     include "F77_CMN_SIZE"     ! Size parameters
#     include "F77_CMN_DIAG"     ! ND13 (for now)
#     include "F77_CMN_GCTM"     ! SCALE_HEIGHT

      ! Arguments
      REAL*8,  INTENT(INOUT) :: TC(IIPAR,JJPAR,LLPAR)

      ! Local variables
      INTEGER                :: I,      J,    L,     NTOP
      REAL*8                 :: DTSRCE, SST,  Sc,    CONC,   W10 
      REAL*8                 :: ScCO2,  AKw,  ERATE, DMSSRC, FEMIS

      ! Molecular weight of DMS, kg/mole
      REAL*8,  PARAMETER     :: DMS_MW = 62d0

      ! Ratio of molecular weights: S/DMS
      REAL*8,  PARAMETER     :: S_DMS = 32d0 / 62d0

      ! External functions
      REAL*8,  EXTERNAL      :: SFCWINDSQR

      !=================================================================
      ! SRCDMS begins here!
      !=================================================================

      ! Chemistry timestep in seconds
      DTSRCE = GET_TS_EMIS() * 60d0

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
         SST = SSTEMP(I,J) - 273.15d0

         ! Only do the following for water boxes
         IF ( IS_WATER(I,J) ) THEN

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
            W10  = SQRT( SFCWINDSQR(I,J) )

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

            DMSSRC = ERATE * GET_AREA_M2( J )

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
               IF (LNLPBL) FEMIS = 0.D0

               ! DMS in box (I,J,L) plus emissions [kg]
               TC(I,J,L) = TC(I,J,L) + ( FEMIS * DMSSRC )

            ENDDO

         ELSE                   

            ! If we are not over water, then there is no DMS source
            DMSSRC = 0.d0

         ENDIF                  
 
         ! Save surface emis for non-local PBL mixing (vdiff_mod.f) (units: kg)
         ! (Lin, 06/09/08)
         if (LNLPBL) emis_save(I,J,IDTDMS) = DMSSRC

         !==============================================================
         ! ND13 diagnostic:  DMS emissions [kg S/box/timestep]
         !==============================================================
         IF ( ND13 > 0 ) THEN
            AD13_DMS(I,J) = AD13_DMS(I,J) + ( DMSSRC * S_DMS ) ! / DTSRCE
         ENDIF
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! Return to calling program
      END SUBROUTINE SRCDMS 

!------------------------------------------------------------------------------

      SUBROUTINE SRCSO2( TC, NSEASON )
!
!******************************************************************************
!  Subroutine SRCSO2 (originally from Mian Chin) computes SO2 emissons from 
!  aircraft, biomass, and anthro sources. (rjp, bdf, bmy, 6/2/00, 12/3/09)
!
!  Arguments as Input/Output:
!  ===========================================================================
!  (1 ) NSEASON (INTEGER) : Season number: 1=DJF; 2=MAM; 3=JJA; 4=SON
!  (2 ) TC      (REAL*8 ) : SO2 tracer mass [kg]
! 
!  NOTES:
!  (1 ) Now reference NSRCE, JDAY, PBL, XTRA2, BXHEIGHT from either header
!        files or F90 modules.  Also use routines from "pressure_mod.f" to
!        compute grid box pressures. (bmy, 9/18/02)
!  (2 ) Now use routines GET_TS_EMIS and GET_DAY_OF_YEAR from the new 
!        "time_mod.f" (bmy, 3/27/03)
!  (3 ) For GEOS-4, convert PBL from [m] to [hPa] w/ the hydrostatic law.
!        Now references SCALE_HEIGHT from "F77_CMN_GCTM".  Added BLTHIK variable
!        to hold PBL thickness in [hPa]. (bmy, 1/15/04)
!  (4 ) Now references AD13_SO2_sh array from "diag_mod.f".  Also references
!        LSHIPSO2 from "F77_CMN_SETUP" (bec, bmy, 5/20/04) 
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
!******************************************************************************
!
      ! Reference to diagnostic arrays
      USE BRAVO_MOD,      ONLY : GET_BRAVO_ANTHRO, GET_BRAVO_MASK
      USE CAC_ANTHRO_MOD, ONLY : GET_CANADA_MASK,  GET_CAC_ANTHRO
      USE DIAG_MOD,       ONLY : AD13_SO2_an,      AD13_SO2_ac
      USE DIAG_MOD,       ONLY : AD13_SO2_bb,      AD13_SO2_nv
      USE DIAG_MOD,       ONLY : AD13_SO2_ev,      AD13_SO2_bf
      USE DIAG_MOD,       ONLY : AD13_SO2_sh
      USE DAO_MOD,        ONLY : BXHEIGHT, PBL
      USE EPA_NEI_MOD,    ONLY : GET_EPA_ANTHRO,   GET_EPA_BIOFUEL
      USE EPA_NEI_MOD,    ONLY : GET_USA_MASK
      USE ERROR_MOD,      ONLY : ERROR_STOP,       GEOS_CHEM_STOP
      USE ERROR_MOD,      ONLY : IS_SAFE_DIV
      USE GRID_MOD,       ONLY : GET_AREA_CM2
      USE GRID_MOD,       ONLY : GET_XOFFSET, GET_YOFFSET
      USE LOGICAL_MOD,    ONLY : LBRAVO, LNEI99,   LSHIPSO2
      USE LOGICAL_MOD,    ONLY : LCAC, LNEI05
      USE NEI2005_ANTHRO_MOD, ONLY : GET_NEI2005_ANTHRO
      USE NEI2005_ANTHRO_MOD, ONLY : NEI05_MASK => USA_MASK
      USE PBL_MIX_MOD,    ONLY : GET_FRAC_OF_PBL,  GET_PBL_TOP_L
      USE PRESSURE_MOD,   ONLY : GET_PEDGE
      USE TIME_MOD,       ONLY : GET_TS_EMIS,      GET_DAY_OF_YEAR 
      USE TIME_MOD,       ONLY : GET_DAY_OF_WEEK
      USE TRACER_MOD,     ONLY : XNUMOL
      USE TRACERID_MOD,   ONLY : IDTSO2
![eml
      USE LOGICAL_MOD,    ONLY : LHIST
!eml]

#     include "F77_CMN_SIZE"     ! Size parameters
#     include "F77_CMN_DIAG"     ! ND13, LD13 (for now)
#     include "F77_CMN_GCTM"     ! SCALE_HEIGHT

      ! Arguments
      INTEGER, INTENT(IN)    :: NSEASON
      REAL*8,  INTENT(INOUT) :: TC(IIPAR,JJPAR,LLPAR)

      ! Local variables
      LOGICAL                :: WEEKDAY, IS_LOCAL
      INTEGER                :: I, J, K, L, NTOP
      INTEGER                :: DAY_NUM
      REAL*8                 :: SO2(LLPAR)
      REAL*8                 :: DTSRCE,      SO2SRC
      REAL*8                 :: TSO2,        FEMIS(LLPAR)
      REAL*8                 :: AREA_CM2,    AN,        BF
      REAL*8                 :: SO2an(IIPAR,JJPAR,NOXLEVELS)
      REAL*8                 :: SO2bf(IIPAR,JJPAR)
      REAL*8                 :: TEMPEMISS(NOXLEVELS)
      REAL*8                 :: T_DIFF, T_SUM

      ! Ratio of molecular weights: S/SO2
      REAL*8,  PARAMETER     :: S_SO2 = 32d0 / 64d0

      ! New variables for volcanic emissions (jaf, 10/15/09)
      INTEGER                :: IREF, JREF, I0, J0
      REAL*8                 :: PLOW, PHIGH, PVOLCLOW, PVOLCHIGH, FRAC

      ! Pressure edges of volcanic 0.5-km grid. Edges come from
      ! converting 500 m increments to pressure using the U.S. standard
      ! atmosphere since volcanic emissions are provided as meters above
      ! sea level (jaf, 10/15/09)
      REAL*8                 :: VOLCPRESS(LVOLC+1) = (/
     &  1013.25, 954.61, 898.74, 845.56, 747.95, 746.83, 701.08,
     &   657.64, 616.40, 577.28, 540.20, 505.07, 471.81, 440.35,
     &   410.61, 382.52, 356.00, 330.99, 307.42, 285.24, 264.36 /)

      !=================================================================
      ! SRCSO2 begins here!
      !================================================================

      ! DTSRCE is the emission timestep in seconds
      DTSRCE  = GET_TS_EMIS() * 60d0

      ! Get current day of the week
      DAY_NUM = GET_DAY_OF_WEEK()

      ! Is it a weekday?
      WEEKDAY = ( DAY_NUM > 0 .and. DAY_NUM < 6 )

      ! is it local PBL mixing scheme?
      IS_LOCAL = .NOT. LNLPBL      
      
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
!$OMP+PRIVATE( I, J, L, AREA_CM2, AN, BF )
!$OMP+SCHEDULE( DYNAMIC )
      DO J = 1, JJPAR

         ! Initialize
         AN       = 0d0
         BF       = 0d0

         ! Grid box surface area [cm2]
         AREA_CM2 = GET_AREA_CM2(J)

         DO I = 1, IIPAR
!------- Prior to 3/8/07 (phs)         
            !!-----------------------------------------------------------
            !! If we are using EPA/NEI99 (anthro + biofuel) ...
            !!-----------------------------------------------------------
            !IF ( LNEI99 ) THEN
            !
            !   ! If we are over the USA ...
            !   IF ( GET_USA_MASK( I, J ) > 0d0 ) THEN
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
            
               ! Read EPA/NEI99 SO2 emissions in [molec/cm2/s]
               AN           = GET_EPA_ANTHRO( I, J, IDTSO2, WEEKDAY )
               BF           = GET_EPA_BIOFUEL(I, J, IDTSO2, WEEKDAY )
               
               ! Convert anthro SO2 from [molec/cm2/s] to [kg/box/s] 
               ! Place all anthro SO2 into surface layer
               SO2an(I,J,1) = AN * AREA_CM2 / XNUMOL(IDTSO2)
               DO L = 2, NOXLEVELS
                  SO2an(I,J,L) = 0d0
               ENDDO
               
               ! Convert biofuel SO2 from [molec/cm2/s] to [kg/box/s] 
               SO2bf(I,J)   = BF * AREA_CM2 / XNUMOL(IDTSO2)

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
                 ENDIF 
               ELSE
                  SO2an(I,J,1) = AN * AREA_CM2 / XNUMOL(IDTSO2)
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
                  endif
               ELSE
                  SO2an(I,J,1) = AN * AREA_CM2 / XNUMOL(IDTSO2)
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

      ! Return to calling program
      END SUBROUTINE SRCSO2

!-----------------------------------------------------------------------------

      SUBROUTINE SRCSO4( TC )
!
!******************************************************************************
!  Subroutine SRCSO4 (originally from Mian Chin) computes SO4 emissions from 
!  anthropogenic sources (rjp, bdf, bmy, 6/2/00, 5/27/09)
!
!  Arguments as Input/Output:
!  ===========================================================================
!  (2) TC     (REAL*8 ) : Array for SO4 mass [kg]
! 
!  NOTES:
!  (1 ) Emission of SO4 is read in SULFATE_READYR, in [kg/box/s]. 
!        It is converted to [kg/box/timestep] here.  
!  (2 ) Now use routine GET_TS_EMIS from the new "time_mod.f" (bmy, 3/27/03)
!  (3 ) For GEOS-4, convert PBL from [m] to [hPa] w/ the barometric law.
!        Now references SCALE_HEIGHT from "F77_CMN_GCTM".  Added BLTHIK variable
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
!******************************************************************************
!
      ! Reference to diagnostic arrays
      USE BRAVO_MOD,          ONLY : GET_BRAVO_ANTHRO, GET_BRAVO_MASK
      USE CAC_ANTHRO_MOD,     ONLY : GET_CANADA_MASK,  GET_CAC_ANTHRO
      USE NEI2005_ANTHRO_MOD, ONLY : GET_NEI2005_ANTHRO
      USE NEI2005_ANTHRO_MOD, ONLY : NEI05_MASK => USA_MASK
      USE DAO_MOD,            ONLY : PBL
      USE DIAG_MOD,           ONLY : AD13_SO4_an,     AD13_SO4_bf
      USE EPA_NEI_MOD,        ONLY : GET_EPA_ANTHRO,  GET_EPA_BIOFUEL
      USE EPA_NEI_MOD,        ONLY : GET_USA_MASK
      USE ERROR_MOD,          ONLY : ERROR_STOP
      USE GRID_MOD,           ONLY : GET_AREA_CM2
      USE LOGICAL_MOD,        ONLY : LNEI99, LCAC, LBRAVO, LNEI05
      USE PBL_MIX_MOD,        ONLY : GET_FRAC_OF_PBL, GET_PBL_TOP_L
      USE TIME_MOD,           ONLY : GET_DAY_OF_WEEK, GET_TS_EMIS
      USE TRACER_MOD,         ONLY : XNUMOL
      USE TRACERID_MOD,       ONLY : IDTSO4, IDTSO2
![eml
      USE LOGICAL_MOD,        ONLY : LHIST
!eml]

#     include "F77_CMN_SIZE"     ! Size parameters
#     include "F77_CMN_DIAG"     ! ND13 (for now)
#     include "F77_CMN_GCTM"     ! SCALE_HEIGHT

      ! Arguments      
      REAL*8,  INTENT(INOUT) :: TC(IIPAR,JJPAR,LLPAR)

      ! Local variables
      LOGICAL                :: WEEKDAY, IS_LOCAL
      INTEGER                :: I, J, K, L, DAY_NUM, NTOP
      REAL*8                 :: SO4(LLPAR), DTSRCE  
      REAL*8                 :: TSO4,       FEMIS(LLPAR)
      REAL*8                 :: AREA_CM2,   EPA_AN,  EPA_BF
      REAL*8                 :: AN
      REAL*8                 :: SO4an(IIPAR,JJPAR,NOXLEVELS)
      REAL*8                 :: SO4bf(IIPAR,JJPAR)
      REAL*8                 :: TEMPEMISS(NOXLEVELS)

      ! Ratio of molecular weights: S/SO4
      REAL*8,  PARAMETER     :: S_SO4 = 32d0 / 96d0

      !=================================================================
      ! SRCSO4 begins here!
      !=================================================================

      ! DTSRCE is the emission timestep in seconds
      DTSRCE  = GET_TS_EMIS() * 60d0

      ! Get current day of the week
      DAY_NUM = GET_DAY_OF_WEEK()

      ! Is it a weekday?
      WEEKDAY = ( DAY_NUM > 0 .and. DAY_NUM < 6 )

      ! is it local PBL mixing scheme?
      IS_LOCAL = .NOT. LNLPBL

      !=================================================================
      ! Overwrite USA    w/ EPA/NEI99 SO4 emissions       (if necessary)
      ! Overwrite CANADA w/ CAC       SO2-fraction emiss. (if necessary)
      ! Overwrite MEXICO w/ BRAVO     SO2-fraction emiss. (if necessary)
      ! Store emissions into local arrays SO4an, SO4bf
      !=================================================================
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, AREA_CM2, EPA_AN, EPA_BF, AN )
!$OMP+PRIVATE( NTOP, FEMIS, SO4, TEMPEMISS )
      DO J = 1, JJPAR
         
         ! Grid box surface area [cm2]
         AREA_CM2 = GET_AREA_CM2(J)

         DO I = 1, IIPAR

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
            
                  ! Read SO4 emissions in [molec/cm2/s]
                  EPA_AN       = GET_EPA_ANTHRO( I, J, IDTSO4, WEEKDAY )
                  EPA_BF       = GET_EPA_BIOFUEL(I, J, IDTSO4, WEEKDAY ) 

                  ! Convert anthro SO4 from [molec/cm2/s] to [kg/box/s] 
                  ! Place all EPA/NEI99 anthro SO4 into surface layer
                  SO4an(I,J,1) = EPA_AN * AREA_CM2 / XNUMOL(IDTSO4)
                  DO L = 2,NOXLEVELS
                     SO4an(I,J,L) = 0d0
                  ENDDO
               
                  ! Convert biofuel SO4 from [molec/cm2/s] to [kg/box/s]
                  SO4bf(I,J)   = EPA_BF * AREA_CM2 / XNUMOL(IDTSO4)

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
               IF ( NEI05_MASK( I, J) < 0d0 ) THEN

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
      ENDDO
!$OMP END PARALLEL DO

      ! Return to calling program
      END SUBROUTINE SRCSO4

!------------------------------------------------------------------------------

      SUBROUTINE SRCNH3( TC )
!
!******************************************************************************
!  Subroutine SRCNH3 handles NH3 emissions into the GEOS-CHEM tracer array.
!  (rjp, bmy, 12/17/01, 5/27/09)
!
!  Arguments as Input/Output
!  ============================================================================
!  (1 ) TC (REAL*8 ) : Array for NH3 tracer mass in kg
!
!  NOTES:
!  (1 ) Now save NH3 emissions to ND13 diagnostic (bmy, 12/13/02)
!  (2 ) Now reference AD13_NH3_na from "diag_mod.f", and archive natural 
!        source NH3 diagnostics for ND13.  Also consider natural source NH3
!        when partitioning by level into the STT array. (rjp, bmy, 3/23/03)
!  (3 ) Now use routine GET_TS_EMIS from the new "time_mod.f" (bmy, 3/27/03)
!  (4 ) For GEOS-4, convert PBL from [m] to [hPa] w/ the barometric law.
!        Now references SCALE_HEIGHT from "F77_CMN_GCTM".  Added BLTHIK variable
!        to hold PBL thickness in [hPa]. (bmy, 1/15/04)
!  (5 ) Now references GET_EPA_ANTHRO, GET_EPA_BIOFUEL, and GET_USA_MASK from 
!        "epa_nei_mod.f".  Now references GET_DAY_OF_WEEK from "time_mod.f".  
!        Now references LNEI99 from "logical_mod.f".  Now references 
!        GET_AREA_CM2 from "grid_mod.f".  Now references IDTNH3 from 
!        "tracerid_mod.f".  Now references XNUMOL from F77_CMN_O3.  Now can 
!        overwrite the anthro & biofuel NH3 emissions over the continental US 
!        if LNEI99=T.  Now references IDTNH3 from "tracerid_mod.f". 
!        (rjp, rch, bmy, 11/16/04)
!  (6 ) Remove reference to "pressure_mod.f".  Now reference GET_FRAC_OF_PBL 
!        and GET_PBL_TOP_L from "pbl_mix_mod.f".  Removed reference to header 
!        file CMN. (bmy, 2/22/05)
!  (7 ) Now references XNUMOL from "tracer_mod.f" (bmy, 10/25/05)
!  (8 ) Need to add CAC_AN to the PRIVATE loop (bmy, 5/27/09)
!  (9 ) Added NIE 2005 (amv, 10/07/2009) 
!  (10) Made NH3an 3D; Calls SULFATE_PBL_MIX to do the PBL mixing of
!        emissions, and allows for emissions above the PBL (phs, 10/27/09)        
!******************************************************************************
!
      ! References to F90 modules
      USE CAC_ANTHRO_MOD,     ONLY : GET_CANADA_MASK
      USE CAC_ANTHRO_MOD,     ONLY : GET_CAC_ANTHRO
      USE NEI2005_ANTHRO_MOD, ONLY : GET_NEI2005_ANTHRO
      USE NEI2005_ANTHRO_MOD, ONLY : NEI05_MASK => USA_MASK
      USE DIAG_MOD,           ONLY : AD13_NH3_an, AD13_NH3_bb
      USE DIAG_MOD,           ONLY : AD13_NH3_bf, AD13_NH3_na
      USE DAO_MOD,            ONLY : PBL
      USE GRID_MOD,           ONLY : GET_AREA_CM2
      USE EPA_NEI_MOD,        ONLY : GET_EPA_ANTHRO, GET_EPA_BIOFUEL
      USE EPA_NEI_MOD,        ONLY : GET_USA_MASK
      USE ERROR_MOD,          ONLY : ERROR_STOP
      USE LOGICAL_MOD,        ONLY : LNEI99, LCAC, LNEI05
      USE PBL_MIX_MOD,        ONLY : GET_FRAC_OF_PBL, GET_PBL_TOP_L
      USE TIME_MOD,           ONLY : GET_DAY_OF_WEEK, GET_TS_EMIS
      USE TRACER_MOD,         ONLY : XNUMOL
      USE TRACERID_MOD,       ONLY : IDTNH3

#     include "F77_CMN_SIZE"     ! Size parameters
#     include "F77_CMN_DIAG"     ! ND13
#     include "F77_CMN_GCTM"     ! SCALE_HEIGHT
      
      ! Argumetns
      REAL*8,  INTENT(INOUT) :: TC(IIPAR,JJPAR,LLPAR)

      ! Local variables
      LOGICAL                :: WEEKDAY,  IS_LOCAL
      INTEGER                :: I, J, L,  K, NTOP, DAY_NUM
      REAL*8                 :: DTSRCE,   AREA_CM2    
      REAL*8                 :: CAC_AN,   EPA_AN,  EPA_BF
      REAL*8                 :: FEMIS(LLPAR), NH3(LLPAR)
      REAL*8                 :: NH3an(IIPAR,JJPAR,NOXLEVELS)
      REAL*8                 :: NEI2005(IIPAR,JJPAR,NOXLEVELS)
      REAL*8                 :: NH3bf(IIPAR,JJPAR)
      REAL*8                 :: TEMPEMISS(NOXLEVELS)
      
      !=================================================================
      ! SRCNH3 begins here!
      !=================================================================

      ! Emission timestep [s]
      DTSRCE  = GET_TS_EMIS() * 60d0
 
      ! Get current day of the week
      DAY_NUM = GET_DAY_OF_WEEK()

      ! Is it a weekday?
      WEEKDAY = ( DAY_NUM > 0 .and. DAY_NUM < 6 )

      ! init
      NH3an = 0D0
      
      ! is it local PBL mixing scheme?
      IS_LOCAL = .NOT. LNLPBL

      !=================================================================
      ! Overwrite USA with EPA/NEI NH3 emissions (if necessary)
      ! Store emissions into local arrays NH3an, NH3bf
      !=================================================================
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, AREA_CM2, EPA_AN, EPA_BF, CAC_AN, L )
!$OMP+PRIVATE( NTOP, FEMIS, NH3, TEMPEMISS )
      DO J = 1, JJPAR
         
         !-------------------------------------------------------------------
         ! NOTE: There seems to be some problems with the EPA/NEI NH3 
         ! emissions.  Therefore we will use the existing emissions for NH3 
         ! until further notice.  Comment out the lines below until 
         ! further notice.  (bmy, 11/17/04) 
         !! Grid box surface area [cm2]
         AREA_CM2 = GET_AREA_CM2( J )
         !-------------------------------------------------------------------

         DO I = 1, IIPAR

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
      ENDDO
!$OMP END PARALLEL DO

      ! Return to calling program
      END SUBROUTINE SRCNH3

!------------------------------------------------------------------------------

      FUNCTION GET_OH( I, J, L ) RESULT( OH_MOLEC_CM3 )
!
!******************************************************************************
!  Function GET_OH returns OH from SMVGEAR's CSPEC array (for coupled runs)
!  or monthly mean OH (for offline runs).  Imposes a diurnal variation on
!  OH for offline simulations. (bmy, 12/16/02, 7/20/04)
!
!  Arguments as Input:
!  ============================================================================
!  (1-3) I, J, L (INTEGER) : Grid box indices for lon, lat, vertical level
!
!  NOTES:
!  (1 ) We assume SETTRACE has been called to define IDOH (bmy, 11/1/02)
!  (2 ) Now use function GET_TS_CHEM from "time_mod.f" (bmy, 3/27/03)
!  (3 ) Now reference ITS_A_FULLCHEM_SIM, ITS_AN_AEROSOL_SIM from 
!        "tracer_mod.f".  Also replace IJSURF w/ an analytic function. 
!        (bmy, 7/20/04)
!******************************************************************************
!
      ! References to F90 modules
      USE COMODE_MOD,    ONLY : CSPEC, JLOP
      USE DAO_MOD,       ONLY : SUNCOS
      USE ERROR_MOD,     ONLY : ERROR_STOP
      USE GLOBAL_OH_MOD, ONLY : OH
      USE TIME_MOD,      ONLY : GET_TS_CHEM
      USE TRACER_MOD,    ONLY : ITS_A_FULLCHEM_SIM, ITS_AN_AEROSOL_SIM
      USE TRACERID_MOD,  ONLY : IDOH

#     include "F77_CMN_SIZE"  ! Size parameters

      ! Arguments
      INTEGER, INTENT(IN) :: I, J, L

      ! Local variables
      INTEGER             :: JLOOP
      REAL*8              :: OH_MOLEC_CM3
 
      !=================================================================
      ! GET_OH begins here!
      !=================================================================
      IF ( ITS_A_FULLCHEM_SIM() ) THEN

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

      ELSE IF ( ITS_AN_AEROSOL_SIM() ) THEN

         !---------------------
         ! Offline simulation
         !---------------------

         ! 1-D grid box index for SUNCOS
         JLOOP = ( (J-1) * IIPAR ) + I

         ! Test for sunlight...
         IF ( SUNCOS(JLOOP) > 0d0 .and. TCOSZ(I,J) > 0d0 ) THEN

            ! Impose a diurnal variation on OH during the day
            OH_MOLEC_CM3 = OH(I,J,L)                      *           
     &                     ( SUNCOS(JLOOP) / TCOSZ(I,J) ) *
     &                     ( 1440d0        / GET_TS_CHEM() )

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

      ! Return to calling program
      END FUNCTION GET_OH

!------------------------------------------------------------------------------

      SUBROUTINE SET_OH( I, J, L, OH ) 
!
!******************************************************************************
!  Function SET_OH saves the modified OH value back to SMVGEAR's CSPEC array
!  for coupled sulfate/aerosol simulations. (bmy, 12/16/02)
!
!  Arguments as Input:
!  ============================================================================
!  (1-3) I, J, L   (INTEGER) : Grid box indices for lon, lat, vertical level
!  (4  ) OH        (REAL*8 ) : OH at grid box (I,J,L) to be saved into CSPEC
!
!  NOTES:
!  (1 ) We assume SETTRACE has been called to define IDOH (bmy, 12/16/02)
!******************************************************************************
!
      ! References to F90 modules
      USE COMODE_MOD,   ONLY : CSPEC, JLOP
      USE TRACERID_MOD, ONLY : IDOH
 
#     include "F77_CMN_SIZE"  ! Size parameters

      ! Arguments
      INTEGER, INTENT(IN) :: I, J, L
      REAL*8,  INTENT(IN) :: OH

      ! Local variables
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

      ! Return to calling program
      END SUBROUTINE SET_OH

!------------------------------------------------------------------------------

      FUNCTION GET_NO3( I, J, L ) RESULT( NO3_MOLEC_CM3 ) 
!
!******************************************************************************
!  Function GET_NO3 returns NO3 from SMVGEAR's CSPEC array (for coupled runs)
!  or monthly mean OH (for offline runs).  For offline runs, the concentration
!  of NO3 is set to zero during the day. (rjp, bmy, 12/16/02)
!
!  Arguments as Input:
!  ============================================================================
!  (1-3) I, J, L (INTEGER) : Grid box indices for lon, lat, vertical level
!
!  NOTES:
!  (1 ) Now references ERROR_STOP from "error_mod.f".  We also assume that
!        SETTRACE has been called to define IDNO3.  Now also set NO3 to 
!        zero during the day. (rjp, bmy, 12/16/02)
!  (2 ) Now reference ITS_A_FULLCHEM_SIM and ITS_AN_AEROSOL_SIM from 
!        "tracer_mod.f".  Also remove reference to CMN.   Also replace
!        IJSURF with an analytic function. (bmy, 7/20/04)
!******************************************************************************
!
      ! References to F90 modules
      USE COMODE_MOD,     ONLY : CSPEC, JLOP
      USE DAO_MOD,        ONLY : AD,    SUNCOS
      USE ERROR_MOD,      ONLY : ERROR_STOP
      USE GLOBAL_NO3_MOD, ONLY : NO3
      USE TRACER_MOD,     ONLY : ITS_A_FULLCHEM_SIM, ITS_AN_AEROSOL_SIM
      USE TRACERID_MOD,   ONLY : IDNO3

#     include "F77_CMN_SIZE"  ! Size parameters

      ! Arguments
      INTEGER, INTENT(IN) :: I, J, L

      ! Local variables
      INTEGER             :: JLOOP
      REAL*8              :: NO3_MOLEC_CM3
 
      ! External functions
      REAL*8,  EXTERNAL   :: BOXVL

      !=================================================================
      ! GET_NO3 begins here!
      !=================================================================
      IF ( ITS_A_FULLCHEM_SIM() ) THEN

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

      ELSE IF ( ITS_AN_AEROSOL_SIM() ) THEN

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
         IF ( SUNCOS(JLOOP) > 0d0 ) THEN

            ! NO3 goes to zero during the day
            NO3_MOLEC_CM3 = 0d0
              
         ELSE

            ! At night: Get NO3 [v/v] and convert it to [kg]
            NO3_MOLEC_CM3 = NO3(I,J,L) * AD(I,J,L) * ( 62d0/28.97d0 ) 
               
            ! Convert NO3 from [kg] to [molec/cm3]
            NO3_MOLEC_CM3 = NO3_MOLEC_CM3 * XNUMOL_NO3 / BOXVL(I,J,L)
                  
         ENDIF
            
         ! Make sure NO3 is not negative
         NO3_MOLEC_CM3  = MAX( NO3_MOLEC_CM3, 0d0 )

      ELSE

         !--------------------
         ! Invalid simulation
         !--------------------
         CALL ERROR_STOP( 'Invalid NSRCX!','GET_NO3 (sulfate_mod.f)')

      ENDIF

      ! Return to calling program
      END FUNCTION GET_NO3

!------------------------------------------------------------------------------

      SUBROUTINE SET_NO3( I, J, L, NO3 ) 
!
!******************************************************************************
!  Function SET_NO3 saves the modified NO3 value back to SMVGEAR's CSPEC array
!  for coupled lfate/aerosol simulations. (rjp, bmy, 12/16/02, 7/20/04)
!
!  Arguments as Input:
!  ============================================================================
!  (1-3) I, J, L (INTEGER) : Grid box indices for lon, lat, vertical level
!  (4  ) NO3     (REAL*8 ) : OH at grid box (I,J,L) to be saved into CSPEC
!
!  NOTES:
!  (1 ) We assume SETTRACE has been called to define IDNO3. (bmy, 12/16/02)
!  (2 ) Remove references to "error_mod.f" and CMN (bmy, 7/20/04)
!******************************************************************************
!
      ! References to F90 modules
      USE COMODE_MOD,   ONLY : CSPEC, JLOP
      USE TRACERID_MOD, ONLY : IDNO3
      
#     include "F77_CMN_SIZE"  ! Size parameters 

      ! Arguments
      INTEGER, INTENT(IN) :: I, J, L
      REAL*8,  INTENT(IN) :: NO3

      ! Local variables
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

      ! Return to calling program
      END SUBROUTINE SET_NO3
      
!------------------------------------------------------------------------------

      FUNCTION GET_O3( I, J, L ) RESULT( O3_VV )
!
!******************************************************************************
!  Function GET_O3 returns monthly mean O3 for offline sulfate aerosol
!  simulations. (bmy, 12/16/02, 10/25/05)
!
!  Arguments as Input:
!  ============================================================================
!  (1-3) I, J, L   (INTEGER) : Grid box indices for lon, lat, vertical level
!
!  NOTES:
!  (1 ) We assume SETTRACE has been called to define IDO3. (bmy, 12/16/02)
!  (2 ) Now reference inquiry functions from "tracer_mod.f" (bmy, 7/20/04)
!  (3 ) Now remove reference to CMN, it's obsolete.  (bmy, 8/22/05)
!  (4 ) Now references XNUMOLAIR from "tracer_mod.f" (bmy, 10/25/05)
!******************************************************************************
!
      ! References to F90 modules
      USE COMODE_MOD,   ONLY : CSPEC, JLOP, VOLUME
      USE DAO_MOD,      ONLY : AIRDEN
      USE ERROR_MOD,    ONLY : ERROR_STOP
      USE TRACER_MOD,   ONLY : ITS_A_FULLCHEM_SIM, ITS_AN_AEROSOL_SIM  
      USE TRACER_MOD,   ONLY : XNUMOLAIR
      USE TRACERID_MOD, ONLY : IDO3

#     include "F77_CMN_SIZE"     ! Size parameters

      ! Arguments
      INTEGER, INTENT(IN)   :: I, J, L

      ! Local variables
      INTEGER               :: JLOOP
      REAL*8                :: O3_VV
 
      !=================================================================
      ! GET_O3 begins here!
      !=================================================================
      IF ( ITS_A_FULLCHEM_SIM() ) THEN

         !--------------------
         ! Coupled simulation
         !--------------------

         ! JLOOP = SMVGEAR 1-D grid box index
         JLOOP = JLOP(I,J,L)

         ! Get O3 from CSPEC [molec/cm3] and convert it to [v/v]
         ! O3 data will only be defined below the tropopause
         IF ( JLOOP  > 0 ) THEN
            O3_VV = ( CSPEC(JLOOP,IDO3) * 1d6       ) / 
     &              ( AIRDEN(L,I,J)     * XNUMOLAIR )
         ELSE
            O3_VV = 0d0
         ENDIF

      ELSE IF ( ITS_AN_AEROSOL_SIM() ) THEN
         
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
         CALL ERROR_STOP( 'Invalid NSRCX!', 'GET_OH (sulfate_mod.f)')

      ENDIF

      ! Return to calling program
      END FUNCTION GET_O3

!------------------------------------------------------------------------------
      
      SUBROUTINE READ_NONERUP_VOLC( INDAY, INMONTH, INYEAR )
!
!******************************************************************************
!  Subroutine READ_NONERUP_VOLC reads SO2 emissions from non-eruptive
!  volcanoes. (rjp, bdf, bmy, jaf, 9/19/02, 10/3/05, 10/15/09)
!
!  Arguments as input:
!  ============================================================================
!  (1 ) INMONTH (INTEGER) : Current month number (1-12)
!  (2 ) INYEAR  (INTEGER) : Current 4-digit year 
!  (3 ) INDAY   (INTEGER) : Current 2-digit day 
!
!  NOTES:
!  (1 ) Split off from old module routine "sulfate_readyr" (bmy, 9/19/02)
!  (2 ) Now references DATA_DIR from "directory_mod.f" (bmy, 7/20/04)
!  (3 ) Now read files from "sulfate_sim_200508/" (bmy, 7/28/05)
!  (4 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (5 ) Complete re-write as volcanic emissions are now monthly and
!	stored as BPCH files (jaf, bmy, 10/15/09)
!  (6 ) Now use MNYEAR_VOLC and MXYEAR_VOLC as 1st and last year of emissions.
!       (ccc, 9/30/10)
!  (7 ) Volcanic data have been updated. Use a new directory. (ccc, 9/30/10)
!******************************************************************************
! 
      ! References to F90 modules
      USE BPCH2_MOD,      ONLY : GET_TAU0, READ_BPCH2
      USE DIRECTORY_MOD,  ONLY : DATA_DIR_1x1
      USE REGRID_1x1_MOD, ONLY : DO_REGRID_1x1, DO_REGRID_G2G_1x1
      USE TIME_MOD,       ONLY : EXPAND_DATE

#     include "F77_CMN_SIZE"  ! Size parameters

      ! Arguments
      INTEGER, INTENT(IN)     :: INDAY, INMONTH, INYEAR
     
      ! Local variables
      REAL*4               :: ARRAY_GEN(I1x1,J1x1-1,LVOLC)
      REAL*8               :: ARRAY_GEN8(I1x1,J1x1-1)
      REAL*8               :: ARRAY_1x1(I1x1,J1x1,1)
      REAL*8               :: ARRAY_1LEV(IIPAR,JJPAR)
      REAL*8               :: ARRAY(IIPAR,JJPAR,LVOLC)
      REAL*8               :: XTAU
      INTEGER              :: THISYEAR, THISDAY, YYYYMMDD, L
      CHARACTER(LEN=255)   :: FILENAME

      !=================================================================
      ! READ_NONERUP_VOLC begins here!
      !=================================================================

      ! Set year based on availability of volcanic emission files
      THISYEAR = INYEAR
!-- Previous to (ccc, 9/30/10). Add parameters for edge years.
!      IF ( INYEAR < 1985 ) THISYEAR = 1985
!      IF ( INYEAR > 2007 ) THISYEAR = 2007
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
!-- Previous to (ccc, 9/30/10). New directory.
!      FILENAME = TRIM( DATA_DIR_1x1 )  //
!     &           'volcano_SO2_200909/' //
!     &           'YYYY/SO2_volc.nonerup.YYYYMM.generic.1x1'

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

      ! Array regridding has to be done level by level since the array
      ! has a different vertical dimension than expected
      DO L = 1, LVOLC
         ! Cast from REAL*4 to REAL*8
         ARRAY_GEN8 = ARRAY_GEN(:,:,L)

         ! Regrid from GENERIC 1x1 GRID to GEOS-Chem 1x1 GRID.
         CALL DO_REGRID_G2G_1x1( 'kg', ARRAY_GEN8, ARRAY_1x1(:,:,1) )

         ! Regrid from 1x1 to correct grid
         CALL DO_REGRID_1x1( 'kg', ARRAY_1x1, ARRAY_1LEV )

         ! Save into correct level of ARRAY
         ARRAY(:,:,L) = ARRAY_1LEV
      ENDDO

      ! Unit conversion: [kg SO2/box/day] -> [kg SO2/box/s]
      ENV = ARRAY / ( 24.d0 * 3600.d0 )

      ! Return to calling program
      END SUBROUTINE READ_NONERUP_VOLC

!------------------------------------------------------------------------------
      
      SUBROUTINE READ_ERUP_VOLC( INDAY, INMONTH, INYEAR )
!
!***************************************************************************** 
!  Subroutine READ_ERUP_VOLC reads SO2 emissions from eruptive
!  volcanoes. (rjp, bdf, bmy, jaf, 9/19/02, 10/3/05, 10/15/09)
!
!  Arguments as input:
!  ============================================================================
!  (1 ) INMONTH (INTEGER) : Current month number (1-12)
!  (2 ) INYEAR  (INTEGER) : Current 4-digit year 
!  (3 ) INDAY   (INTEGER) : Current 2-digit day 
!
!  NOTES:
!  (1 ) Split off from old module routine "sulfate_readyr" (bmy, 9/19/02)
!  (2 ) Now references DATA_DIR from "directory_mod.f" (bmy, 7/20/04)
!  (3 ) Now read files from "sulfate_sim_200508/" (bmy, 7/28/05)
!  (4 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (5 ) Complete re-write as volcanic emissions are now monthly and
!	stored as BPCH files (jaf, bmy, 10/15/09)
!  (6 ) Now use MNYEAR_VOLC and MXYEAR_VOLC as 1st and last year of emissions.
!       (ccc, 9/30/10)
!  (7 ) Volcanic data have been updated. Use a new directory. (ccc, 9/30/10)
!*****************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD,      ONLY : GET_TAU0, READ_BPCH2
      USE DIRECTORY_MOD,  ONLY : DATA_DIR_1x1
      USE REGRID_1x1_MOD, ONLY : DO_REGRID_1x1, DO_REGRID_G2G_1x1
      USE TIME_MOD,       ONLY : EXPAND_DATE
      
#     include "F77_CMN_SIZE"   ! Size parameters

      ! Arguments
      INTEGER, INTENT(IN)     :: INDAY, INMONTH, INYEAR
 
      ! Local variables
      REAL*4               :: ARRAY_GEN(I1x1,J1x1-1,LVOLC)
      REAL*8               :: ARRAY_GEN8(I1x1,J1x1-1)
      REAL*8               :: ARRAY_1x1(I1x1,J1x1,1)
      REAL*8               :: ARRAY_1LEV(IIPAR,JJPAR)
      REAL*8               :: ARRAY(IIPAR,JJPAR,LVOLC)
      REAL*8               :: XTAU
      INTEGER              :: YYYYMMDD, L
      CHARACTER(LEN=255)   :: FILENAME
     
      !==================================================================
      ! READ_ERUP_VOLC begins here
      !==================================================================

      ! If the current year falls within the range of available data,
      ! get eruptive volcanic emissions (jaf, 10/15/09)
!-- Previous to (ccc, 9/30/10). Add parameters for edge years
!      IF ( ( INYEAR .GE. 1985 ) .AND. ( INYEAR .LE. 2007 ) ) THEN
      IF ( ( INYEAR .GE. MNYEAR_VOLC )  .AND. 
     &     ( INYEAR .LE. MXYEAR_VOLC ) ) THEN

         ! Set date
         YYYYMMDD = 10000 * INYEAR + 100 * INMONTH + INDAY

         ! File name
!-- Previous to (ccc, 9/30/10). New directory.
!         FILENAME = TRIM( DATA_DIR_1x1 )  //
!     &              'volcano_SO2_200909/' //
!     &              'YYYY/SO2_volc.erup.YYYYMM.generic.1x1'

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

         ! Array regridding has to be done level by level since the array
         ! has a different vertical dimension than expected
         DO L = 1, LVOLC
            ! Cast from REAL*4 to REAL*8
            ARRAY_GEN8 = ARRAY_GEN(:,:,L)

            ! Regrid from GENERIC 1x1 GRID to GEOS-Chem 1x1 GRID.
            CALL DO_REGRID_G2G_1x1( 'kg', ARRAY_GEN8, ARRAY_1x1(:,:,1) )

            ! Regrid from 1x1 to correct grid
            CALL DO_REGRID_1x1( 'kg', ARRAY_1x1, ARRAY_1LEV )

            ! Save into correct level of ARRAY
            ARRAY(:,:,L) = ARRAY_1LEV
         ENDDO

         ! Unit conversion: [kg SO2/box/day] -> [kg SO2/box/s]
         EEV = ARRAY / ( 24.d0 * 3600.d0 )

      ! If no data are available for current year, set eruptive
      ! emissions to zero and print a warning message.
      ! Current data range is MNYEAR_VOLC - MXYEAR_VOLC (jaf, 10/15/09)
      ! Year parameters are defined at the beginning of the module.
      ! (ccc, 9/30/10)
      ELSE

!-- Previous to (ccc, 9/30/10). Add parameters for edge years
!         WRITE(6, *) 'WARNING: Eruptive SO2 emissions only available'//
!     &               ' for 1985-2007. You are outside the window.'
!         WRITE(6, *) 'Eruptive SO2 emissions are being set to zero!'


         WRITE(6, 110) MNYEAR_VOLC, MXYEAR_VOLC
         WRITE(6, *) 'Eruptive SO2 emissions are being set to zero!'

 110     FORMAT( 'WARNING: Eruptive SO2 emissions only available for ', 
     &        i4,'-',i4,'. You are outside the window.')


         EEV = 0d0

      ENDIF

      ! Return to calling program
      END SUBROUTINE READ_ERUP_VOLC
      
!------------------------------------------------------------------------------

      SUBROUTINE READ_ANTHRO_SOx( THISMONTH, NSEASON )
!
!******************************************************************************
!  Suborutine READ_ANTHRO_SOx reads the anthropogenic SOx from disk, 
!  and partitions it into anthropogenic SO2 and SO4. 
!  (rjp, bdf, bmy, 9/20/02, 10/31/08)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) THISMONTH (INTEGER) : Current month number (1-12)
!  
!  NOTES:
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
!******************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD,            ONLY : GET_NAME_EXT_2D, GET_RES_EXT
      USE BPCH2_MOD,            ONLY : GET_TAU0,        READ_BPCH2
      USE DIRECTORY_MOD,        ONLY : DATA_DIR
      USE EDGAR_MOD,            ONLY : GET_EDGAR_ANTH_SO2
      USE EMEP_MOD,             ONLY : GET_EMEP_ANTHRO
      USE EMEP_MOD,             ONLY : GET_EUROPE_MASK
      USE FUTURE_EMISSIONS_MOD, ONLY : GET_FUTURE_SCALE_SO2ff
      USE GRID_MOD,             ONLY : GET_XMID,        GET_YMID
      USE GRID_MOD,             ONLY : GET_AREA_CM2
      USE LOGICAL_MOD,          ONLY : LFUTURE,         LEDGARSOx
      USE LOGICAL_MOD,          ONLY : LSTREETS,        LEMEP
      USE STREETS_ANTHRO_MOD,   ONLY : GET_SE_ASIA_MASK
      USE STREETS_ANTHRO_MOD,   ONLY : GET_STREETS_ANTHRO
      USE TIME_MOD,             ONLY : GET_YEAR
      USE TRACER_MOD,           ONLY : XNUMOL
      USE TRACERID_MOD,         ONLY : IDTSO2, IDTSO4
      USE TRANSFER_MOD,         ONLY : TRANSFER_2D
      USE SCALE_ANTHRO_MOD,     ONLY : GET_ANNUAL_SCALAR
![eml
      USE LOGICAL_MOD,          ONLY : LHIST
      USE TIME_MOD,             ONLY : GET_HISTYR
!eml]

#     include "F77_CMN_SIZE"             ! Size parameters
#     include "F77_CMN_O3"               ! FSCALYR
                                    
      ! Arguments                   
      INTEGER, INTENT(IN)           :: THISMONTH, NSEASON
                                    
      ! Local variables             
      INTEGER                       :: I, J, L, LL, IX, JX, IOS
      INTEGER, SAVE                 :: LASTYEAR = -99
      INTEGER                       :: SCALEYEAR
      REAL*4                        :: E_SOx(IIPAR,JJPAR,2)
      REAL*4                        :: ARRAY(IIPAR,JJPAR,1)
      REAL*8                        :: XTAU, Fe, X, Y, AREA_CM2
      REAL*8                        :: EDG_SO2
      CHARACTER(LEN=4)              :: SYEAR
      CHARACTER(LEN=255)            :: FILENAME

      !=================================================================
      ! READ_ANTHRO_SOx begins here!
      !=================================================================
         
      IF ( FSCALYR < 0 ) THEN
         SCALEYEAR = GET_YEAR()
      ELSE
         SCALEYEAR = FSCALYR
      ENDIF


      IF ( LEDGARSOx ) THEN

         !==============================================================
         ! Use EDGAR SOx emissions
         !
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
!$OMP+PRIVATE( I, J, X, Y, EDG_SO2, Fe, LL ) 

         ! Loop over latitudes
         DO J = 1, JJPAR

            ! Latitude [degrees]
            Y = GET_YMID( J )

            ! Loop over longitudes
            DO I = 1, IIPAR

               ! Longitude [degrees]
               X  = GET_XMID( I )

               ! Get EDGAR SO2 emissions [kg/s]
               ! NOTE: Future emissions are already applied!
               EDG_SO2 = GET_EDGAR_ANTH_SO2( I, J, KG_S=.TRUE. ) 

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

            ! Grid box surface area [cm2]
            AREA_CM2 = GET_AREA_CM2( J )

            ! Latitude [degrees]
            Y = GET_YMID( J )

            DO I = 1, IIPAR

               ! Longitude [degrees]  
               X = GET_XMID( I )

               ! First scale SOx to the given fossil fuel year
               E_SOx(I,J,L) = E_SOx(I,J,L) * SOx_SCALE(I,J)
            
               ! Compute future SOx emissions (if necessary)
               IF ( LFUTURE ) THEN
                  E_SOx(I,J,L) = E_SOx(I,J,L)                  * 
     &                           GET_FUTURE_SCALE_SO2ff( I, J )
               ENDIF

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
         
               ! Compute SO2 (tracer #2) from SOx
               ! Convert from [molec SOx/cm2/s] to [kg SO2/box/s]
               ESO2_an(I,J,L) = E_SOx(I,J,L) * ( 1.d0 - Fe ) * 
     &                          AREA_CM2 / XNUMOL(IDTSO2)            

               ! If we are using David Streets' emissions
               ! Remember: those include BF if Year is GE 2006
               IF ( LSTREETS ) THEN

                  ! If we are over the SE Asia region
                  IF ( GET_SE_ASIA_MASK( I, J ) > 0d0 ) THEN

                     ! Overwrite GEIA SO2 w/ David Streets' SO2 [kg SO2/s]
                     ESO2_an(I,J,1) = GET_STREETS_ANTHRO( I, J, IDTSO2, 
     &                                                    KG_S=.TRUE. )

                     ! Zero 2nd level of emissions
                     ESO2_an(I,J,2) = 0d0

                  ENDIF
               ENDIF

               IF ( LEMEP ) THEN

                  IF (GET_EUROPE_MASK(I,J) > 0d0 ) THEN

                     ESO2_an(I,J,1) = GET_EMEP_ANTHRO( I, J, IDTSO2,
     &                                                 KG_S=.TRUE. ) 

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

      ! Return to calling program
      END SUBROUTINE READ_ANTHRO_SOx
 
!------------------------------------------------------------------------------

      SUBROUTINE READ_OCEAN_DMS( THISMONTH )
!
!***************************************************************************** 
!  Subroutine READ_OCEAN_DMS reads seawater concentrations of DMS (nmol/L).
!  (rjp, bdf, bmy, 9/20/02, 10/3/05)
!
!  Arguments as input:
!  ===========================================================================
!  (1 ) THISMONTH (INTEGER) : Current month number (1-12)
!
!  NOTES:
!  (1 ) Extracted from old module routine SULFATE_READMON (bmy, 9/18/02)
!  (2 ) Now call READ_BPCH2 with QUIET=.TRUE. (bmy, 3/27/03)
!  (3 ) Now references DATA_DIR from "directory_mod.f" (bmy, 7/20/04)
!  (4 ) Now read files from "sulfate_sim_200508/".  Now read data for both
!        GCAP and GEOS grids (bmy, 8/16/05) 
!  (5 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!*****************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD,     ONLY : GET_NAME_EXT_2D, GET_RES_EXT
      USE BPCH2_MOD,     ONLY : GET_TAU0,        READ_BPCH2
      USE DIRECTORY_MOD, ONLY : DATA_DIR
      USE TRANSFER_MOD,  ONLY : TRANSFER_2D

#     include "F77_CMN_SIZE"  ! Size parameters 

      ! Arguments
      INTEGER, INTENT(IN) :: THISMONTH
      
      ! Local variables
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
      
      ! Return to calling program
      END SUBROUTINE READ_OCEAN_DMS

!------------------------------------------------------------------------------

      SUBROUTINE READ_SST( THISMONTH, THISYEAR )
!
!***************************************************************************** 
!  Subroutine READ_SST reads monthly mean sea surface temperatures.
!  (rjp, bdf, bmy, 9/18/02, 7/13/09)
!
!  Arguments as input:
!  ===========================================================================
!  (1 ) THISMONTH (INTEGER) : Current month number (1-12)
!  (2 ) THISYEAR  (INTEGER) : Current 4-digit year 
!
!  NOTES:
!  (1 ) Extracted from old module routine SULFATE_READMON (bmy, 9/18/02)
!  (2 ) Now call READ_BPCH2 with QUIET=.TRUE. (bmy, 3/27/03)
!  (3 ) Now references DATA_DIR from "directory_mod.f" (bmy, 7/20/04)
!  (4 ) Now use interannual SST data from NOAA if present; otherwise use
!        climatological SST data.  Now read data for both GCAP and GEOS 
!        grids (bmy, 8/16/05) 
!  (5 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (6 ) Now read int'annual SST data on the GEOS 1x1 grid (bmy, 11/17/05)
!  (7 ) Last year of data is now 2008 (bmy, 7/13/09)
!*****************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD,      ONLY : GET_NAME_EXT_2D, GET_RES_EXT
      USE BPCH2_MOD,      ONLY : GET_TAU0,        READ_BPCH2
      USE DIRECTORY_MOD,  ONLY : DATA_DIR,        DATA_DIR_1x1
      USE REGRID_1x1_MOD, ONLY : DO_REGRID_1x1
      USE TRANSFER_MOD,   ONLY : TRANSFER_2D

#     include "F77_CMN_SIZE"       ! Size parameters

      ! Arguments
      INTEGER, INTENT(IN)     :: THISMONTH, THISYEAR
      
      ! Local variables
      REAL*4                  :: ARRAY(IIPAR,JJPAR,1)
      REAL*4                  :: ARRAY1(I1x1,J1x1,1)
      REAL*8                  :: XTAU
      CHARACTER(LEN=4)        :: SYEAR
      CHARACTER(LEN=255)      :: FILENAME

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

         ! Regrid from 1x1 and cast to REAL*8
         CALL DO_REGRID_1x1( 'K', ARRAY1, SSTEMP )

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

      ! Return to calling program
      END SUBROUTINE READ_SST

!------------------------------------------------------------------------------

      SUBROUTINE READ_BIOFUEL_SO2( THISMONTH )
!
!******************************************************************************
!  Subroutine READ_BIOFUEL_SO2 reads monthly mean biomass burning 
!  emissions for SO2.  SOx is read in, and converted to SO2. 
!  (rjp, bdf, bmy, phs, 1/16/03, 12/23/08)
!
!  Arguments as input:
!  ===========================================================================
!  (1 ) THISMONTH (INTEGER) : Current month number (1-12)
!
!  NOTES:
!  (1 ) Extracted from old module routine SULFATE_READMON (bmy, 9/18/02)
!  (2 ) Modified molar ratio of biomass burning SO2 per CO.  Added SO2
!        emission from biofuel burning. (rjp, bmy, 1/16/03)
!  (3 ) Now replace DXYP(J+J0)*1d4 with routine GET_AREA_CM2 of "grid_mod.f"
!        Now replace MONTH with the argument THISMONTH.  Now call READ_BPCH2
!        with QUIET=.TRUE. (bmy, 3/27/03)
!  (4 ) Now references DATA_DIR from "directory_mod.f".  Also removed
!        references to CMN and F77_CMN_SETUP. (bmy, 7/20/04)
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
!******************************************************************************
!
      ! References to F90 modules
      USE BIOMASS_MOD,          ONLY : BIOMASS
      ! IDBSO2 is not used anymore (ccc, 01/29/10)
      !USE BIOMASS_MOD,          ONLY : BIOMASS,         IDBSO2
      USE BPCH2_MOD,            ONLY : GET_NAME_EXT_2D, GET_RES_EXT
      USE BPCH2_MOD,            ONLY : GET_TAU0,        READ_BPCH2
      USE DIRECTORY_MOD,        ONLY : DATA_DIR
      USE FUTURE_EMISSIONS_MOD, ONLY : GET_FUTURE_SCALE_SO2bf
      USE GRID_MOD,             ONLY : GET_AREA_CM2
      USE LOGICAL_MOD,          ONLY : LBIOMASS,        LFUTURE
      USE TIME_MOD,             ONLY : ITS_A_LEAPYEAR
      USE TRACER_MOD,           ONLY : XNUMOL
      USE TRACERID_MOD,         ONLY : IDTSO2
      USE TRANSFER_MOD,         ONLY : TRANSFER_2D
![eml
      USE LOGICAL_MOD,          ONLY : LHIST
      USE TIME_MOD,             ONLY : GET_HISTYR
      USE REGRID_1x1_MOD,       ONLY : DO_REGRID_G2G_1x1, DO_REGRID_1x1
!eml]

#     include "F77_CMN_SIZE"             ! Size parameters
                                    
      ! Arguments                   
      INTEGER, INTENT(IN)           :: THISMONTH
                                    
      ! Local variables             
      INTEGER                       :: I, J, THISYEAR
      REAL*4                        :: ARRAY(IIPAR,JJPAR,1)
      REAL*8                        :: BIOCO(IIPAR,JJPAR)
      REAL*8                        :: XTAU
      CHARACTER(LEN=4  )            :: CYEAR
      CHARACTER(LEN=255)            :: FILENAME
![eml
      CHARACTER(LEN=4  )            :: CHISTYR
      CHARACTER(LEN=255)            :: DATADIR
      REAL*4                        :: ESO2_GEN(I1x1,J1x1-1)
      REAL*8                        :: ESO2_GEN8(I1x1,J1x1-1)
      REAL*8                        :: ESO2_GEO(I1x1,J1x1,1)
      REAL*8                        :: TEMP_ED(IIPAR,JJPAR,1)
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
         WRITE(CHISTYR,'(I4)') GET_HISTYR()
         DATADIR ='/home/eml/data/EMISSIONS/SO2/proper/'
         FILENAME=TRIM(DATADIR)//'EDGAR.BCSO2.'//CHISTYR//'generic.1x1'
         
         XTAU = GET_TAU0( 1, 1, GET_HISTYR() )
         WRITE( 6, 101 ) TRIM( FILENAME )
 101     FORMAT( '     - READ_EDGAR_NOx: Reading ', a )
         CALL READ_BPCH2( TRIM(FILENAME), 'EDGAR-2D', 26,
     &                    XTAU,           I1x1,       J1x1-1,
     &                    1,              ESO2_GEN,   QUIET=.TRUE. )

         ESO2_GEN8 = ESO2_GEN

         CALL DO_REGRID_G2G_1x1( 'kg/yr', ESO2_GEN8, ESO2_GEO(:,:,1) )
         CALL DO_REGRID_1x1(1, 'kg/yr', ESO2_GEO, TEMP_ED )
         ! Convert from kg SO2/yr to kg SO2/s
         ESO2_bf = TEMP_ED(:,:,1) / (365.25d0 * 86400.d0)
      ENDIF
!eml]

      ! Return to calling program
      END SUBROUTINE READ_BIOFUEL_SO2

!------------------------------------------------------------------------------

      SUBROUTINE GET_BIOMASS_SO2
!
!******************************************************************************
!  Subroutine GET_BIOMASS_SO2 retrieve monthly/8-day/3hr biomass burning 
!  emissions for SO2.  (yc, phs, 12/23/08)
!
!  Arguments as input:
!  ===========================================================================
!  (1 ) 
!
!  NOTES:
!  (1 ) Extracted from old module subroutine READ_BIOMASS_SO2
!        (yc, phs, 12/23/08)
!  (2 ) IDBSO2 is now in tracerid_mod.f (fp, 6/2009)
!******************************************************************************
!
      ! References to F90 modules
      !IDBSO2 now in tracerid_mod (fp, 6/2009)
      !USE BIOMASS_MOD,          ONLY : BIOMASS,         IDBSO2
      USE BIOMASS_MOD,          ONLY : BIOMASS
      USE TRACERID_MOD,         ONLY : IDBSO2
      USE GRID_MOD,             ONLY : GET_AREA_CM2
      USE TRACER_MOD,           ONLY : XNUMOL
      USE TRACERID_MOD,         ONLY : IDTSO2
      USE TRANSFER_MOD,         ONLY : TRANSFER_2D
      
#     include "F77_CMN_SIZE"             ! Size parameters
                                    
      ! Local variables             
      INTEGER                       :: I, J
      REAL*8                        :: CONV
                                    
      !=================================================================
      ! GET_BIOMASS_SO2 begins here!
      !=================================================================
      ! Unit conversion to [kg SO2/s]
  
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, CONV )

      ! Loop over longitudes
      DO J = 1, JJPAR

         ! Conversion factor for [cm2 * kg/molec]
         CONV = GET_AREA_CM2( J ) / XNUMOL(IDTSO2)

         ! Loop over latitudes
         DO I = 1, IIPAR

            ! Convert biomass SO2 from [molec SO2/cm2/s] -> [kg SO2/s] 
            ! NOTE: Future scale has already been applied by this point
            ESO2_bb(I,J) = BIOMASS(I,J,IDBSO2) * CONV

         ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! Return to calling program
      END SUBROUTINE GET_BIOMASS_SO2

!------------------------------------------------------------------------------

      SUBROUTINE READ_AIRCRAFT_SO2( THISMONTH )
!
!******************************************************************************
!  Subroutine READ_AIRCRAFT_SO2 reads monthly mean aircraft fuel emissions 
!  and converts them to SO2 emissions. (rjp, bdf, bmy, 9/18/02, 10/3/05)
!
!  Arguments as input:
!  ===========================================================================
!  (1 ) THISMONTH (INTEGER) : Current month number (1-12)
!
!  NOTES:
!  (1 ) Extracted from old module routine SULFATE_READMON (bmy, 9/18/02)
!  (2 ) Now references DATA_DIR from "directory_mod.f" (bmy, 7/20/04)
!  (3 ) Now read files from "sulfate_sim_200508/".  Now read data for both 
!        GCAP and GEOS grids (bmy, 8/16/05)
!  (4 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!******************************************************************************
!
      ! Reference to F90 modules
      USE BPCH2_MOD,     ONLY : GET_RES_EXT, GET_TAU0, READ_BPCH2
      USE DAO_MOD,       ONLY : BXHEIGHT
      USE DIRECTORY_MOD, ONLY : DATA_DIR 
      USE FILE_MOD,      ONLY : IU_FILE, IOERROR

#     include "F77_CMN_SIZE"  ! Size parameters

      ! Arguments
      INTEGER, INTENT(IN) :: THISMONTH

      ! Local variables
      INTEGER             :: I, IOS, J, K, L
      REAL*8              :: ACSO2(IIPAR,JJPAR,20)
      REAL*8              :: FAC, FUEL, DZ(LLPAR), ZH(0:LLPAR)
      CHARACTER(LEN=255)  :: FILENAME

      ! Month names
      CHARACTER(LEN=3)    :: CMONTH(12) = (/'jan', 'feb', 'mar', 'apr', 
     &                                      'may', 'jun', 'jul', 'aug',
     &                                      'sep', 'oct', 'nov', 'dec'/)

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
               ZH(L) = ZH(L-1) + ( BXHEIGHT(I,J,L) * 1d-3 )
               
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

      ! Return to calling program
      END SUBROUTINE READ_AIRCRAFT_SO2

!------------------------------------------------------------------------------

      SUBROUTINE READ_SHIP_SO2( THISMONTH )
!
!******************************************************************************
!  Subroutine READ_SHIP_SO2 reads in ship SO2 emissions, from either Corbett
!  et al or EDGAR inventories. (bec, qli, 10/01/03, 7/14/06)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) THISMONTH (INTEGER) : Current month (1-12)
!
!  NOTES:
!  (1 ) Now references DATA_DIR from "directory_mod.f" (bmy, 7/20/04)
!  (2 ) Now read files from "sulfate_sim_200508/".  Now read data for both 
!        GCAP and GEOS grids. (bmy, 8/16/05)
!  (3 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (4 ) Now references XNUMOL from "tracer_mod.f" (bmy, 10/25/05)
!  (5 ) Now get EDGAR ship SO2 emissions if necessary.  Also apply future
!        emissions scale factors to the default Corbett et al ship emissions.
!        (avd, bmy, 7/14/06)
!  (6 ) Now references GET_ARCTAS_HIP from 'arctas_ship_emiss_mod.f" and
!        GET_EMEP_ANTHRO to get ARCTAS and EMEP SO2 ship emissions (phs, 12/5/08)
!  (7 ) Now get ICOADS ship SO2 if necessary (phs, cklee, 6/30/09) 
!******************************************************************************
!
      ! References to F90 modules
      USE ARCTAS_SHIP_EMISS_MOD,ONLY : GET_ARCTAS_SHIP
      USE ICOADS_SHIP_MOD,      ONLY : GET_ICOADS_SHIP   !(cklee, 7/09/09)
      USE BPCH2_MOD,            ONLY : GET_NAME_EXT_2D, GET_RES_EXT
      USE BPCH2_MOD,            ONLY : GET_TAU0,        READ_BPCH2
      USE DIRECTORY_MOD,        ONLY : DATA_DIR 
      USE EDGAR_MOD,            ONLY : GET_EDGAR_SHIP_SO2
      USE EMEP_MOD,             ONLY : GET_EMEP_ANTHRO, GET_EUROPE_MASK
      USE FUTURE_EMISSIONS_MOD, ONLY : GET_FUTURE_SCALE_SO2ff
      USE GRID_MOD,             ONLY : GET_AREA_CM2
      USE LOGICAL_MOD,          ONLY : LEDGARSHIP,      LFUTURE, 
     &                                 LARCSHIP,        LSHIPSO2,
     $                                 LEMEPSHIP
      USE LOGICAL_MOD,          ONLY : LICOADSSHIP  !(cklee, 6/30/09)
      USE TRACER_MOD,           ONLY : XNUMOL
      USE TRACERID_MOD,         ONLY : IDTSO2
      USE TRANSFER_MOD,         ONLY : TRANSFER_2D

#     include "F77_CMN_SIZE"             ! Size parameters 

      ! Arguments
      INTEGER, INTENT(IN)           :: THISMONTH

      ! Local variables
      INTEGER                       :: I, J
      REAL*4                        :: ARRAY(IIPAR,JJPAR,1)
      REAL*4                        :: SHIPSO2(IIPAR,JJPAR)
      REAL*8                        :: XTAU, AREA_CM2
      CHARACTER (LEN=255)           :: FILENAME

      !=================================================================
      ! READ_SHIP_SO2 begins here!
      !=================================================================

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

         ! Loop over latitudes
         DO J = 1, JJPAR

            ! Grid box surface area [cm2]
            AREA_CM2 = GET_AREA_CM2( J )

            ! Loop over longitudes
            DO I = 1, IIPAR
            
               ! Convert [molec SO2/cm2/s] to [kg SO2/box/s]
               ESO2_sh(I,J) = SHIPSO2(I,J) * AREA_CM2 / XNUMOL(IDTSO2)

               ! Apply future emissions (if necessary)
               IF ( LFUTURE ) THEN
                   ESO2_sh(I,J) = ESO2_sh(I,J) *
     &                            GET_FUTURE_SCALE_SO2ff( I, J ) 
               ENDIF

               IF ( LEMEPSHIP ) THEN
                  IF ( GET_EUROPE_MASK(I,J) > 0d0 )
     $                 ESO2_sh(I,J) = GET_EMEP_ANTHRO(I, J, IDTSO2,
     $                                         KG_S=.TRUE., SHIP=.TRUE.)
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


      ! Return to calling program
      END SUBROUTINE READ_SHIP_SO2

!------------------------------------------------------------------------------

      SUBROUTINE READ_ANTHRO_NH3( THISMONTH )
!
!******************************************************************************
!  Subroutine READ_ANTHRO_NH3 reads the monthly mean anthropogenic 
!  NH3 emissions from disk and converts to [kg NH3/box/s]. 
!  (rjp, bdf, bmy, 9/20/02, 3/5/10)
!
!  Arguments as input:
!  ===========================================================================
!  (1 ) THISMONTH (INTEGER) : Current month number (1-12)
!
!  NOTES:
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
!******************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD,            ONLY : GET_NAME_EXT_2D, GET_RES_EXT
      USE BPCH2_MOD,            ONLY : GET_TAU0,        READ_BPCH2
      USE DIRECTORY_MOD,        ONLY : DATA_DIR
      USE EMEP_MOD,             ONLY : GET_EMEP_ANTHRO
      USe EMEP_MOD,             ONLY : GET_EUROPE_MASK
      USE FUTURE_EMISSIONS_MOD, ONLY : GET_FUTURE_SCALE_NH3an
      USE LOGICAL_MOD,          ONLY : LFUTURE,         LSTREETS
      USE LOGICAL_MOD,          ONLY : LEMEP
      USE STREETS_ANTHRO_MOD,   ONLY : GET_SE_ASIA_MASK
      USE STREETS_ANTHRO_MOD,   ONLY : GET_STREETS_ANTHRO
      USE TRACERID_MOD,         ONLY : IDTNH3
      USE TRANSFER_MOD,         ONLY : TRANSFER_2D

#     include "F77_CMN_SIZE"             ! Size parameters

      ! Arguments
      INTEGER, INTENT(IN)           :: THISMONTH

      ! Local variables
      LOGICAL                       :: WEEKDAY
      INTEGER                       :: I, J, DAY_NUM
      REAL*4                        :: ARRAY(IIPAR,JJPAR,1)
      REAL*8                        :: AREA_CM2, EPA_NEI, XTAU
      REAL*8              :: NMDAY(12) = (/ 31d0, 28d0, 31d0, 30d0,
     &                                      31d0, 30d0, 31d0, 31d0, 
     &                                      30d0, 31d0, 30d0, 31d0 /)
      CHARACTER(LEN=255)            :: FILENAME
      REAL*8                        :: STREETS

      !=================================================================
      ! READ_ANTHRO_NH3 begins here!
      !=================================================================

      ! File name
      FILENAME = TRIM( DATA_DIR )                         //
     &           'sulfate_sim_200508/NH3_anthsrce.'       //
     &           GET_NAME_EXT_2D() // '.' // GET_RES_EXT()

      ! Echo output
      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( '     - READ_ANTHRO_NH3: Reading ', a )
      
      ! Get TAU value (use year 1990, the year of the data!)
      XTAU = GET_TAU0( THISMONTH, 1, 1990 )
	
      ! Read 1990 NH3 emissions [kg N/box/mon]
      CALL READ_BPCH2( FILENAME, 'NH3-ANTH',    29,  
     &                 XTAU,      IIPAR,        JJPAR,       
     &                 1,         ARRAY(:,:,1), QUIET=.TRUE. )

      ! Cast from REAL*4 to REAL*8
      CALL TRANSFER_2D( ARRAY(:,:,1), ENH3_an )

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, STREETS )
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Convert from [kg N/box/mon] to [kg NH3/box/s]
         ENH3_an(I,J) = ENH3_an(I,J) * ( 17.d0 / 14.d0 ) 
     &                / ( NMDAY(THISMONTH) * 86400.d0 )

         ! Compute future NH3an emissions, if necessary
         ! Moved here since Streets and EMEP should have already
         ! applied FUTURE scale factors if needed 
         IF ( LFUTURE ) THEN
            ENH3_an(I,J) = ENH3_an(I,J) * GET_FUTURE_SCALE_NH3an( I, J )
         ENDIF

         ! If we are using David Streets' emissions ...
         IF ( LSTREETS ) THEN

            ! If we are over the SE Asia region ...
            IF ( GET_SE_ASIA_MASK( I, J ) > 0d0  ) THEN

               ! Overwrite with David Streets emissions if available [kg NH3/s]
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

      ! Return to calling program
      END SUBROUTINE READ_ANTHRO_NH3

!------------------------------------------------------------------------------

      SUBROUTINE READ_NATURAL_NH3( THISMONTH )
!
!******************************************************************************
!  Subroutine READ_NATURAL_NH3 reads the monthly mean natural 
!  NH3 emissions from disk and converts to [kg NH3/box/s]. 
!  (rjp, bdf, bmy, 9/20/02, 10/3/05)
!
!  Arguments as input:
!  ===========================================================================
!  (1 ) THISMONTH (INTEGER) : Current month number (1-12)
!
!  NOTES:
!  (1 ) Updated FORMAT string.  Now also call READ_BPCH2 with QUIET=.TRUE.
!        (bmy, 4/8/03)
!  (2 ) Now references DATA_DIR from "directory_mod.f" (bmy, 7/20/04)
!  (3 ) Now read files from "sulfate_sim_200508/".  Now read data for both 
!        GCAP and GEOS grids. (bmy, 8/16/05)
!  (4 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!******************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD,     ONLY : GET_NAME_EXT_2D, GET_RES_EXT
      USE BPCH2_MOD,     ONLY : GET_TAU0,        READ_BPCH2
      USE DIRECTORY_MOD, ONLY : DATA_DIR
      USE TRANSFER_MOD,  ONLY : TRANSFER_2D

#     include "F77_CMN_SIZE"  ! Size parameters

      ! Arguments
      INTEGER, INTENT(IN) :: THISMONTH

      ! Local variables
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
 
      ! Return to calling program
      END SUBROUTINE READ_NATURAL_NH3

!------------------------------------------------------------------------------

      SUBROUTINE READ_BIOFUEL_NH3( THISMONTH ) 
!
!******************************************************************************
!  Subroutine READ_BIOFUEL_NH3 reads the monthly mean biomass NH3 
!  and biofuel emissions from disk and converts to [kg NH3/box/s]. 
!  (rjp, bdf, bmy, phs, 9/20/02, 12/23/08)
!
!  Arguments as input:
!  ===========================================================================
!  (1 ) THISMONTH (INTEGER) : Current month number (1-12)
!
!  NOTES:
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
!        from "logical_mod.f".  Removed references to CMN and F77_CMN_SETUP.
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
!******************************************************************************
!
      ! References to F90 modules
      !USE BIOMASS_MOD,          ONLY : BIOMASS,         IDBNH3
      USE BIOMASS_MOD,          ONLY : BIOMASS
      USE BPCH2_MOD,            ONLY : GET_NAME_EXT_2D, GET_RES_EXT
      USE BPCH2_MOD,            ONLY : GET_TAU0,        READ_BPCH2
      USE DIRECTORY_MOD,        ONLY : DATA_DIR
      USE FUTURE_EMISSIONS_MOD, ONLY : GET_FUTURE_SCALE_NH3bf
      USE LOGICAL_MOD,          ONLY : LBIOMASS,        LFUTURE
      USE GRID_MOD,             ONLY : GET_AREA_CM2
      USE TIME_MOD,             ONLY : ITS_A_LEAPYEAR
      USE TRACER_MOD,           ONLY : XNUMOL
      USE TRACERID_MOD,         ONLY : IDTNH3
      USE TRANSFER_MOD,         ONLY : TRANSFER_2D

#     include "F77_CMN_SIZE"             ! Size parameters

      ! Arguments
      INTEGER, INTENT(IN)           :: THISMONTH

      ! Local variables
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
      
      ! Return to calling program
      END SUBROUTINE READ_BIOFUEL_NH3

!------------------------------------------------------------------------------

      SUBROUTINE GET_BIOMASS_NH3
!
!******************************************************************************
!  Subroutine GET_BIOMASS_NH3 retrieve the monthly/8days/3hr mean biomass NH3 
!  (yc, phs, 12/23/08)
!
!  Arguments as input:
!  ===========================================================================
!  (1 ) 
!
!  NOTES:
!  (1 ) Extracted from old module subroutine READ_BIOMASS_NH3
!        (yc, phs, 12/23/08)
!  (2 ) IDBNH3 is in tracerid_mod.f now (fp, 6/2009)
!******************************************************************************
!
      ! References to F90 modules
      ! IDBNH3 now in tracerid_mod (fp, 6/2009)
      !USE BIOMASS_MOD,          ONLY : BIOMASS,         IDBNH3
      USE BIOMASS_MOD,          ONLY : BIOMASS
      USE GRID_MOD,             ONLY : GET_AREA_CM2
      USE TRACER_MOD,           ONLY : XNUMOL
      USE TRACERID_MOD,         ONLY : IDTNH3, IDBNH3

#     include "F77_CMN_SIZE"             ! Size parameters

      ! Local variables
      INTEGER                       :: I, J
      REAL*8                        :: CONV
      
      !=================================================================
      ! READ_BIOMASSBURN_NH3 begins here!
      !=================================================================
      ! Convert units

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, CONV )

      ! Loop over latitudes
      DO J = 1, JJPAR
 
         ! Conversion factor for [cm2 * kg/molec]
         CONV = GET_AREA_CM2( J ) / XNUMOL(IDTNH3)

         ! Loop over longitudes
         DO I = 1, IIPAR

            ! Convert biomass NH3 from [molec NH3/cm2/s] -> [kg NH3/s]
            ! NOTE: Future scale is applied by this point (if necessary)
            ENH3_bb(I,J) = BIOMASS(I,J,IDBNH3) * CONV

         ENDDO
      ENDDO
!$OMP END PARALLEL DO
      
      ! Return to calling program
      END SUBROUTINE GET_BIOMASS_NH3
      
!------------------------------------------------------------------------------

      SUBROUTINE READ_OXIDANT( MONTH )
!
!******************************************************************************
!  Subroutine READ_OXIDANT reads in monthly mean H2O2 and O3 fields for the
!  offline sulfate + aerosol simulation. (rjp, bdf, bmy, 11/1/02, 10/3/05)
!
!  Arguments as input:
!  ============================================================================
!  (1 ) MONTH    (INTEGER  ) : Emission timestep in minutes
!
!  NOTES:
!  (1 ) Now call READ_BPCH2 with QUIET=.TRUE. (bmy, 3/27/03)
!  (2 ) Now references DATA_DIR from "directory_mod.f" (bmy, 7/20/04)
!  (3 ) Now read files from "sulfate_sim_200508/offline/".  Now read data
!        for both GEOS and GCAP grids (bmy, 8/16/05)
!  (4 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!******************************************************************************
!  
      ! References to F90 modules
      USE BPCH2_MOD,     ONLY : GET_NAME_EXT, GET_RES_EXT
      USE BPCH2_MOD,     ONLY : GET_TAU0,     READ_BPCH2
      USE DIRECTORY_MOD, ONLY : DATA_DIR
      USE TRANSFER_MOD,  ONLY : TRANSFER_3D_TROP

#     include "F77_CMN_SIZE"  ! Size parameters

      ! Arguments
      INTEGER, INTENT(IN) :: MONTH

      ! Local variables 
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
      
      ! Return to calling program
      END SUBROUTINE READ_OXIDANT

!------------------------------------------------------------------------------

      SUBROUTINE OHNO3TIME
!
!******************************************************************************
!  Subroutine OHNO3TIME computes the sum of cosine of the solar zenith
!  angle over a 24 hour day, as well as the total length of daylight. 
!  This is needed to scale the offline OH and NO3 concentrations.
!  (rjp, bmy, 12/16/02, 3/30/04)
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
!******************************************************************************
!
      ! References to F90 modules
      USE GRID_MOD, ONLY : GET_XMID,    GET_YMID_R
      USE TIME_MOD, ONLY : GET_NHMSb,   GET_ELAPSED_SEC
      USE TIME_MOD, ONLY : GET_TS_CHEM, GET_DAY_OF_YEAR, GET_GMT

#     include "F77_CMN_SIZE"  ! Size parameters
#     include "F77_CMN_GCTM"

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

            ! IJLOOP is the 1-D loop index for SUNCOS
            IJLOOP = 0

            ! Loop over surface grid boxes
            DO J = 1, JJPAR

               ! Grid box latitude center [radians]
               YMID_R = GET_YMID_R( J )

            DO I = 1, IIPAR

               ! Increment IJLOOP
               IJLOOP = IJLOOP + 1
               TIMLOC = real(LHR0) + real(NT)/3600.0 + GET_XMID(I)/15.0
         
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

      ! Return to calling program
      END SUBROUTINE OHNO3TIME

!------------------------------------------------------------------------------
          
      SUBROUTINE INIT_SULFATE
!
!******************************************************************************
!  Subroutine INIT_SULFATE initializes and zeros all allocatable arrays
!  declared in "sulfate_mod.f" (bmy, 6/2/00, 10/15/09)
!
!  NOTES:
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
!******************************************************************************
!
      ! References to F90 modules
      USE DRYDEP_MOD,    ONLY : DEPNAME, NUMDEP
      USE ERROR_MOD,     ONLY : ALLOC_ERR
      USE LOGICAL_MOD,   ONLY : LDRYD
      USE TRACER_MOD,    ONLY : ITS_AN_AEROSOL_SIM

#     include "F77_CMN_SIZE" ! Size parameters

      ! Local variables
      LOGICAL, SAVE      :: IS_INIT = .FALSE.
      INTEGER            :: AS, I, J, N, IJLOOP

      !=================================================================
      ! INIT_SULFATE begins here!
      !=================================================================

      ! Return if we have already initialized
      IF ( IS_INIT ) RETURN
      
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

#if   !defined( GEOS_5 ) && !defined( MERRA )
      ! If we are using GEOS-5 or MERRA met, then get the cloud fraction 
      ! directly from the met fields.  (skim, bmy, 1/14/10)
      ALLOCATE( VCLDF( IIPAR, JJPAR, LLTROP ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'VCLDF' )
      VCLDF = 0d0
#endif

      !=================================================================
      ! Only initialize the following for offline aerosol simulations
      !=================================================================
      IF ( ITS_AN_AEROSOL_SIM() ) THEN

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

      IF ( LDRYD ) THEN
         
         ! Locate position of each tracer in DEPSAV
         DO N = 1, NUMDEP
            SELECT CASE ( TRIM( DEPNAME(N) ) )
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

      ! Reset IS_INIT so we do not allocate arrays again
      IS_INIT = .TRUE.

      ! Return to calling program
      END SUBROUTINE INIT_SULFATE

!-----------------------------------------------------------------------------

      SUBROUTINE CLEANUP_SULFATE
!
!******************************************************************************
!  Subroutine CLEANUP_SULFATE deallocates all previously allocated arrays 
!  for sulfate emissions -- call at the end of the run (bmy, 6/1/00, 10/15/09)
! 
!  NOTES:
!  (1 ) Now also deallocates IJSURF. (bmy, 11/12/02)
!  (2 ) Now also deallocates ENH3_na (rjp, bmy, 3/23/03)
!  (3 ) Now also deallocates COSZM (rjp, bmy, 3/30/04)
!  (4 ) Now also deallocates ESO4_sh (bec, bmy, 5/20/04)
!  (5 ) Now remove IJSURF (bmy, 7/20/04)
!  (6 ) Bug fix: now deallocate PSO4_ss, PNITs (bmy, 5/3/06)
!  (7 ) Deleted obsolete arrays from older SO2 volcanic emissions 
!        inventory (jaf, bmy, 10/15/09)
!******************************************************************************
! 
      !=================================================================
      ! CLEANUP_SULFATE begins here!
      !=================================================================
      IF ( ALLOCATED( DMSo      ) ) DEALLOCATE( DMSo      )
      IF ( ALLOCATED( EEV       ) ) DEALLOCATE( EEV       )
      IF ( ALLOCATED( ENV       ) ) DEALLOCATE( ENV       )
      IF ( ALLOCATED( ENH3_an   ) ) DEALLOCATE( ENH3_an   )
      IF ( ALLOCATED( ENH3_bb   ) ) DEALLOCATE( ENH3_bb   )
      IF ( ALLOCATED( ENH3_bf   ) ) DEALLOCATE( ENH3_bf   )
      IF ( ALLOCATED( ENH3_na   ) ) DEALLOCATE( ENH3_na   )
      IF ( ALLOCATED( ESO2_ac   ) ) DEALLOCATE( ESO2_ac   )
      IF ( ALLOCATED( ESO2_an   ) ) DEALLOCATE( ESO2_an   )
      IF ( ALLOCATED( ESO2_nv   ) ) DEALLOCATE( ESO2_nv   )
      IF ( ALLOCATED( ESO2_ev   ) ) DEALLOCATE( ESO2_ev   )
      IF ( ALLOCATED( ESO2_bb   ) ) DEALLOCATE( ESO2_bb   )
      IF ( ALLOCATED( ESO2_bf   ) ) DEALLOCATE( ESO2_bf   )
      IF ( ALLOCATED( ESO2_sh   ) ) DEALLOCATE( ESO2_sh   )
      IF ( ALLOCATED( ESO4_an   ) ) DEALLOCATE( ESO4_an   )
      IF ( ALLOCATED( JH2O2     ) ) DEALLOCATE( JH2O2     )
      IF ( ALLOCATED( O3m       ) ) DEALLOCATE( O3m       )
      IF ( ALLOCATED( PH2O2m    ) ) DEALLOCATE( PH2O2m    )
      IF ( ALLOCATED( PMSA_DMS  ) ) DEALLOCATE( PMSA_DMS  )
      IF ( ALLOCATED( PNITs     ) ) DEALLOCATE( PNITs     )
      IF ( ALLOCATED( PSO2_DMS  ) ) DEALLOCATE( PSO2_DMS  )
      IF ( ALLOCATED( PSO4_SO2  ) ) DEALLOCATE( PSO4_SO2  )
      IF ( ALLOCATED( PSO4_ss   ) ) DEALLOCATE( PSO4_ss   )
      IF ( ALLOCATED( SOx_SCALE ) ) DEALLOCATE( SOx_SCALE )
      IF ( ALLOCATED( SSTEMP    ) ) DEALLOCATE( SSTEMP    )
      IF ( ALLOCATED( TCOSZ     ) ) DEALLOCATE( TCOSZ     )
      IF ( ALLOCATED( TTDAY     ) ) DEALLOCATE( TTDAY     )          
      IF ( ALLOCATED( VCLDF     ) ) DEALLOCATE( VCLDF     )
      IF ( ALLOCATED( COSZM     ) ) DEALLOCATE( COSZM     )

      ! Return to calling program
      END SUBROUTINE CLEANUP_SULFATE

!------------------------------------------------------------------------------

      END MODULE SULFATE_MOD
